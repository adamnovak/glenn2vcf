// main.cpp: Main file for core graph merger

#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <set>
#include <utility>
#include <algorithm>
#include <getopt.h>

#include "ekg/vg/vg.hpp"
#include "ekg/vg/index.hpp"
#include "ekg/vg/vcflib/src/Variant.h"

// Note: THIS CODE IS TERRIBLE
// TODO:
//  - Decide if we need to have sibling alts detect (somehow) and coordinate with each other
//  - Parallelize variant generation
//  - Make variant stamping out some kind of function, don't duplicate the same variant construction code 6 times

// How many bases may we put in an allele in VCF if we expect GATK to be able to
// parse it?
const static int MAX_ALLELE_LENGTH = 4096;

/**
 * Holds indexes of the reference: position to node, node to position and
 * orientation, and the full reference string.
 */
struct ReferenceIndex {
    // Index from node ID to first position on the reference string and
    // orientation it occurs there.
    std::map<int64_t, std::pair<size_t, bool>> byId;
    
    // Index from start position on the reference to the oriented node that
    // begins there.  Some nodes may be backward (orientation true) at their
    // canonical reference positions. In this case, the last base of the node
    // occurs at the given position.
    std::map<size_t, vg::NodeTraversal> byStart;
    
    // The actual sequence of the reference.
    std::string sequence;
};

/**
 * Make a letter into a full string because apparently that's too fancy for the
 * standard library.
 */
std::string char_to_string(const char& letter) {
    std::string toReturn;
    toReturn.push_back(letter);
    return toReturn;
}

/**
 * Write a minimal VCF header for a single-sample file.
 */
void write_vcf_header(std::ostream& stream, std::string& sample_name, std::string& contig_name, size_t contig_size) {
    stream << "##fileformat=VCFv4.2" << std::endl;
    stream << "##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location\">" << std::endl;
    stream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
    if(!contig_name.empty()) {
        // Announce the contig as well.
        stream << "##contig=<ID=" << contig_name << ",length=" << contig_size << ">" << std::endl;
    }
    stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sample_name << std::endl;
}

/**
 * Create the reference allele for an empty vcflib Variant, since apaprently
 * there's no method for that already. Must be called before any alt alleles are
 * added.
 */
void create_ref_allele(vcflib::Variant& variant, const std::string& allele) {
    // Set the ref allele
    variant.ref = allele;
    
    for(size_t i = 0; i < variant.ref.size(); i++) {
        // Look at all the bases
        if(variant.ref[i] != 'A' && variant.ref[i] != 'C' && variant.ref[i] != 'G' && variant.ref[i] != 'T') {
            // Correct anything bogus (like "X") to N
            variant.ref[i] = 'N';
        }
    }
    
    // Make it 0 in the alleles-by-index list
    variant.alleles.push_back(allele);
    // Build the reciprocal index-by-allele mapping
    variant.updateAlleleIndexes();
}

/**
 * Return true if a variant may be output, or false if this variant is valid but
 * the GATK might choke on it.
 *
 * Mostly used to throw out variants with very logn alleles, because GATK has an
 * allele length limit. How alleles that really *are* 1 megabase deletions are
 * to be specified to GATK is left as an exercise to the reader.
 */
bool can_write_alleles(vcflib::Variant& variant) {
    for(auto& allele : variant.alleles) {
        if(allele.size() > MAX_ALLELE_LENGTH) {
            return false;
        }
    }
    return true;
}

/**
 * Add a new alt allele to a vcflib Variant, since apaprently there's no method
 * for that already.
 *
 * If that allele already exists in the variant, does not add it again.
 *
 * Retuerns the allele number (0, 1, 2, etc.) corresponding to the given allele
 * string in the given variant. 
 */
int add_alt_allele(vcflib::Variant& variant, const std::string& allele) {
    // Copy the allele so we can throw out bad characters
    std::string fixed(allele);
    
    for(size_t i = 0; i < fixed.size(); i++) {
        // Look at all the bases
        if(fixed[i] != 'A' && fixed[i] != 'C' && fixed[i] != 'G' && fixed[i] != 'T') {
            // Correct anything bogus (like "X") to N
            fixed[i] = 'N';
        }
    }
    
    for(int i = 0; i < variant.alleles.size(); i++) {
        if(variant.alleles[i] == fixed) {
            // Already exists
            return i;
        }
    }

    // Add it as an alt
    variant.alt.push_back(fixed);
    // Make it next in the alleles-by-index list
    variant.alleles.push_back(fixed);
    // Build the reciprocal index-by-allele mapping
    variant.updateAlleleIndexes();

    // We added it in at the end
    return variant.alleles.size() - 1;
}

/**
 * Return true if a mapping is a perfect match, and false if it isn't.
 */
bool mapping_is_perfect_match(const vg::Mapping& mapping) {
    for (auto edit : mapping.edit()) {
        if (edit.from_length() != edit.to_length() || !edit.sequence().empty()) {
            // This edit isn't a perfect match
            return false;
        }
    }
    
    // If we get here, all the edits are perfect matches.
    // Note that Mappings with no edits at all are full-length perfect matches.
    return true;
}

/**
 * Do a breadth-first search left from the given node traversal, and return
 * maximal-copy-number, minimal-length-in-nodes node list paths starting at the
 * given node and ending on the indexed reference path. The paths are broken up
 * into lists in order of descending copy number.
 */
std::vector<std::list<std::list<vg::NodeTraversal>>> bfs_left(vg::VG& graph,
    vg::NodeTraversal node, const ReferenceIndex& index,
    const std::map<vg::Node*, size_t>& nodeCopyNumbers, int64_t maxDepth = 10) {

    // We'll fill this in based on copy number we can push. Right now we fill it
    // in with two empty lists, because we don't support copy numbers greater
    // than 2. Slot 0 corresponds to copy number 1, slot 1 corresponds to copy
    // number 2.
    std::vector<std::list<std::list<vg::NodeTraversal>>> toReturn { {}, {} };
    
    // Do a BFS
    
    // This holds the paths to get to NodeTraversals to visit (all of which will
    // end with the node we're starting with) and the max copy number that can
    // be pushed through each.
    std::list<std::pair<std::list<vg::NodeTraversal>, size_t>> toExtend;
    
    // This keeps a set of all the oriented nodes we already got to and don't
    // need to queue again, stratified by copy number. Slot 0 corresponds to
    // copy number 1, slot 1 corresponds to copy number 2.
    std::vector<std::set<vg::NodeTraversal>> alreadyQueued { {}, {} };
    
    // Start at this node at depth 0, with the copy number it has
    toExtend.emplace_back(std::make_pair(std::list<vg::NodeTraversal> {node}, nodeCopyNumbers.at(node.node)));
    if(nodeCopyNumbers.at(node.node) == 0) {
        // We can't find any paths for this node since it is absent.
        return toReturn;
    }
    assert(nodeCopyNumbers.at(node.node) <= 2);
    // Find the set for this copy number and put in the node
    alreadyQueued.at(nodeCopyNumbers.at(node.node) - 1).insert(node);
    
    while(!toExtend.empty()) {
        // Keep going until we've visited every node up to our max search depth.
        
        // Dequeue a path to extend and the copy number we're looking for from it.
        std::pair<std::list<vg::NodeTraversal>, size_t> nextToExtend = toExtend.front();
        toExtend.pop_front();
        auto& path = nextToExtend.first;
        auto& flow = nextToExtend.second;
        
        // We can't just throw out longer paths, because shorter paths may need
        // to visit a node twice (in opposite orientations) and thus might get
        // rejected later. Or they might overlap with paths on the other side.
        
        // Look up and see if the front node on the path is on our reference
        // path
        if(index.byId.count(path.front().node->id())) {
            // This node is on the reference path. TODO: we don't care if it
            // lands in a place that is itself deleted.
            
            // Say we got to the right place
            toReturn.at(flow - 1).push_back(path);
            
            // Don't bother looking for extensions, we already got there. If
            // this path gets disqualified, extensions will also get
            // disqualified.
        } else if(path.size() <= maxDepth) {
            // We haven't hit the reference path yet, but we also haven't hit
            // the max depth. Extend with all the possible extensions.
            
            // Look left
            vector<vg::NodeTraversal> prevNodes;
            graph.nodes_prev(path.front(), prevNodes);
            
            for(auto prevNode : prevNodes) {
                // For each node we can get to
                
                // Is this an anchoring node?
                bool isPrimaryPath = index.byId.count(prevNode.node->id());
            
                // What's its copy number effect? We just say we can push all
                // the flow if we hit the reference.
                auto nodeCopyNumber = isPrimaryPath ? flow : nodeCopyNumbers.at(prevNode.node);
                
                // How much flow can we push if we add it to the path
                auto newMaxPushable = std::min(nodeCopyNumber, flow);
                
                if(newMaxPushable == 0) {
                    // Skip this node; we can't use any copies of it.
                    continue;
                }
                // We aren't built to handle copy numbers > 2.
                assert(newMaxPushable <= 2);
            
                // Can we already get here with this copy number or more?
                bool canAlreadyReachWithThisCopyNumberOrGreater = false;
                for(size_t i = alreadyQueued.size(); i >= newMaxPushable; i--) {
                    // Run the loop from high copy numbers to low copy numbers
                    if(alreadyQueued.at(i - 1).count(prevNode)) {
                        canAlreadyReachWithThisCopyNumberOrGreater = true;
                    }
                }
                if(canAlreadyReachWithThisCopyNumberOrGreater) {
                    // If we can, we don't need to add this node. TODO:
                    // shouldn't we care about finding all the paths? Since some
                    // path combinations can be invalid?
                    continue;
                }
            
            
                // Make a new path extended left with the node
                std::list<vg::NodeTraversal> extended(path);
                extended.push_front(prevNode);
                toExtend.push_back(std::make_pair(extended, newMaxPushable));
                
                // Remember we found a way to this node, so we don't try and
                // visit it other ways.
                alreadyQueued.at(newMaxPushable - 1).insert(prevNode);
            }
        }
        
    }
    
    // When we get here, we've found at least one of the paths to reference
    // nodes at each possible length, with max pushable copy number.
    return toReturn;
}

/**
 * Flip a NodeTraversal around and return the flipped copy.
 */
vg::NodeTraversal flip(vg::NodeTraversal toFlip) {
    return vg::NodeTraversal(toFlip.node, !toFlip.backward);
}

/**
 * Do a breadth-first search right from the given node traversal, and return
 * maximal-copy-number, minimal-length-in-nodes node list paths starting at the
 * given node and ending on the indexed reference path. The paths are broken up
 * into lists in order of descending copy number.
 */
std::vector<std::list<std::list<vg::NodeTraversal>>> bfs_right(vg::VG& graph,
    vg::NodeTraversal node, const ReferenceIndex& index,
    const std::map<vg::Node*, size_t>& nodeCopyNumbers, int64_t maxDepth = 10) {

    // Look left from the backward version of the node.
    std::vector<std::list<std::list<vg::NodeTraversal>>> toReturn = bfs_left(graph, flip(node), index, nodeCopyNumbers, maxDepth);
    
    for(auto& pathGroup : toReturn) {
        for(auto& path : pathGroup) {
            // Invert the order of every path in palce
            path.reverse();
            
            for(auto& traversal : path) {
                // And invert the orientation of every node in the path in place.
                traversal = flip(traversal);
            }
        }
    }
    
    return toReturn;
}

#define debug
/**
 * Given a vg graph, a node in the graph, and an index for the reference path,
 * look out from the node in both directions to find a bubble relative to the
 * path, with a consistent orientation, such that the bubble has the maximum
 * copy number pushable along it, and such that it is the shortest bubble on
 * each side for that copy number. The bubble may not visit the same node twice.
 * 
 * Return the ordered and oriented nodes in the bubble, with the outer nodes
 * being oriented forward along the named path, and with the first node coming
 * before the last node in the reference. Needs a set of nodes identified as
 * entirely present, so it can search paths using only those nodes.
 */
std::pair<std::vector<vg::NodeTraversal>, size_t>
find_bubble(vg::VG& graph, vg::Node* node, const ReferenceIndex& index,
const std::map<vg::Node*, size_t>& nodeCopyNumbers) {

    // Find paths on both sides, with nodes on the primary path at the outsides
    // and this node in the middle. Returns paths stratified by copy number, so
    // you want to comapre the first pair of lists first, then things in each
    // first list with the other second list, then things in each first and
    // second list with the other third list, and so on (if we ever went that
    // far).
    auto leftPaths = bfs_left(graph, vg::NodeTraversal(node), index, nodeCopyNumbers);
    auto rightPaths = bfs_right(graph, vg::NodeTraversal(node), index, nodeCopyNumbers);
    
    // Find a combination of two paths which gets us to the reference in a
    // consistent orientation (meaning that when you look at the ending nodes'
    // Mappings in the reference path, the ones with minimal ranks have the same
    // orientations) and which doesn't use the same nodes on both sides.
    
    // We need to look in different combinations of lists.
    auto testCombination = [&](const std::list<std::list<vg::NodeTraversal>>& leftList,
        const std::list<std::list<vg::NodeTraversal>>& rightList) {

        for(auto leftPath : leftList) {
            // Figure out the relative orientation for the leftmost node.
#ifdef debug        
            std::cerr << "Left path: " << std::endl;
            for(auto traversal : leftPath ) {
                std::cerr << "\t" << traversal << std::endl;
            }
#endif    
            // Split out its node pointer and orientation
            auto leftNode = leftPath.front().node;
            auto leftOrientation = leftPath.front().backward;
            
            // Get where it falls in the reference as a position, orientation pair.
            auto leftRefPos = index.byId.at(leftNode->id());
            
            // We have a backward orientation relative to the reference path if we
            // were traversing the anchoring node backwards, xor if it is backwards
            // in the reference path.
            bool leftRelativeOrientation = leftOrientation != leftRefPos.second;
            
            // Make a set of all the nodes in the left path
            std::set<int64_t> leftPathNodes;
            for(auto visit : leftPath) {
                leftPathNodes.insert(visit.node->id());
            }
            
            for(auto rightPath : rightList) {
                // Figure out the relative orientation for the rightmost node.
#ifdef debug            
                std::cerr << "Right path: " << std::endl;
                for(auto traversal : rightPath ) {
                    std::cerr << "\t" << traversal << std::endl;
                }
#endif            
                // Split out its node pointer and orientation
                // Remember it's at the end of this path.
                auto rightNode = rightPath.back().node;
                auto rightOrientation = rightPath.back().backward;
                
                // Get where it falls in the reference as a position, orientation pair.
                auto rightRefPos = index.byId.at(rightNode->id());
                
                // We have a backward orientation relative to the reference path if we
                // were traversing the anchoring node backwards, xor if it is backwards
                // in the reference path.
                bool rightRelativeOrientation = rightOrientation != rightRefPos.second;
                
                if(leftRelativeOrientation == rightRelativeOrientation &&
                    ((!leftRelativeOrientation && leftRefPos.first < rightRefPos.first) ||
                    (leftRelativeOrientation && leftRefPos.first > rightRefPos.first))) {
                    // We found a pair of paths that get us to and from the
                    // reference without turning around, and that don't go back to
                    // the reference before they leave.
                    
                    // Start with the left path
                    std::vector<vg::NodeTraversal> fullPath{leftPath.begin(), leftPath.end()};
                    
                    // We need to detect overlap with the left path
                    bool overlap = false;
                    
                    for(auto it = ++(rightPath.begin()); it != rightPath.end(); ++it) {
                        // For all but the first node on the right path, add that in
                        fullPath.push_back(*it);
                        
                        if(leftPathNodes.count((*it).node->id())) {
                            // We already visited this node on the left side. Try
                            // the next right path instead.
                            overlap = true;
                        }
                    }
                    
                    if(overlap) {
                        // Can't combine this right with this left, as they share
                        // nodes and we can't handle the copy number implications.
                        // Try the next right.
                        // TODO: handle the copy number implications.
                        continue;
                    }
                    
                    if(leftRelativeOrientation) {
                        // Turns out our anchored path is backwards.
                        
                        // Reorder everything the other way
                        std::reverse(fullPath.begin(), fullPath.end());
                        
                        for(auto& traversal : fullPath) {
                            // Flip each traversal
                            traversal = flip(traversal);
                        }
                    }
                    
                    // Count up max copy number pushable.
                    // TODO: we check the main node twice.
                    size_t maxPushable = nodeCopyNumbers.at(node);
                    for(auto traversal : fullPath) {
                        // Limit to the narrowest node on the path, not counting the reference nodes.
                        if(!index.byId.count(traversal.node->id())) {
                            maxPushable = std::min(maxPushable, nodeCopyNumbers.at(traversal.node));
                        }
                    }
                    
                    // Just give the first valid path we find.
#ifdef debug        
                    std::cerr << "Merged path:" << std::endl;
                    for(auto traversal : fullPath) {
                        std::cerr << "\t" << traversal << std::endl;
                    }
#endif
                    return std::make_pair(fullPath, maxPushable);
                }
                
            }
        }
        
        // Return no copy number through no path if we can't find anything.
        return std::make_pair(std::vector<vg::NodeTraversal>(), (size_t) 0);
        
    };
    
    // Try all the groups of lists stratified by copy number and return the best
    // thing we find. TODO: since copy numbers can only go up to 2 we're just
    // going to hardcode this. We really should generate it.
    assert(leftPaths.size() <= 2);
    assert(rightPaths.size() <= 2);
    
    std::vector<std::pair<int, int>> ordering = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
    for(auto indexes : ordering) {
        if(indexes.first >= leftPaths.size() || indexes.second >= rightPaths.size()) {
            // We don't have that many path groups
            continue;
        }
        
        // Look for a valid combination in this tranche
        auto result = testCombination(leftPaths[indexes.first], rightPaths[indexes.second]);
        
        if(result.second > 0) {
            // We found one!
            return result;
        }
    }
    
    // No combinations found in any tranche.
    return std::make_pair(std::vector<vg::NodeTraversal>(), 0);
}
#undef debug

/**
 * Trace out the reference path in the given graph named by the given name.
 * Returns a structure with useful indexes of the reference.
 */
ReferenceIndex trace_reference_path(vg::VG& vg, std::string refPathName) {
    // Make sure the reference path is present
    assert(vg.paths.has_path(refPathName));
    
    // We'll fill this in and then return it.
    ReferenceIndex index;
    
    // We're also going to build the reference sequence string
    std::stringstream refSeqStream;
    
    // What base are we at in the reference
    size_t referenceBase = 0;
    
    // What was the last rank? Ranks must always go up.
    int64_t lastRank = -1;
    
    for(auto mapping : vg.paths.get_path(refPathName)) {
        // All the mappings need to be perfect matches.
        assert(mapping_is_perfect_match(mapping));
    
        if(!index.byId.count(mapping.position().node_id())) {
            // This is the first time we have visited this node in the reference
            // path.
            
            // Add in a mapping.
            index.byId[mapping.position().node_id()] = 
                std::make_pair(referenceBase, mapping.is_reverse());
#ifdef debug
            std::cerr << "Node " << mapping.position().node_id() << " rank " << mapping.rank()
                << " starts at base " << referenceBase << " with "
                << vg.get_node(mapping.position().node_id())->sequence() << std::endl;
#endif
            
            // Make sure ranks are monotonically increasing along the path.
            assert(mapping.rank() > lastRank);
            lastRank = mapping.rank();
        }
        
        // Find the node's sequence
        std::string sequence = vg.get_node(mapping.position().node_id())->sequence();
        
        while(referenceBase == 0 && sequence.size() > 0 &&
            (sequence[0] != 'A' && sequence[0] != 'T' && sequence[0] != 'C' &&
            sequence[0] != 'G' && sequence[0] != 'N')) {
            
            // If the path leads with invalid characters (like "X"), throw them
            // out when computing reference path positions.
            
            // TODO: this is a hack to deal with the debruijn-brca1-k63 graph,
            // which leads with an X.
            
            std::cerr << "Warning: dropping invalid leading character "
                << sequence[0] << " from node " << mapping.position().node_id()
                << std::endl;
                
            sequence.erase(sequence.begin());
        }
        
        if(mapping.is_reverse()) {
            // Put the reverse sequence in the reference path
            refSeqStream << vg::reverse_complement(sequence);
        } else {
            // Put the forward sequence in the reference path
            refSeqStream << sequence;
        }
            
        // Say that this node appears here along the reference in this
        // orientation.
        index.byStart[referenceBase] = vg::NodeTraversal(
            vg.get_node(mapping.position().node_id()), mapping.is_reverse()); 
            
        // Whether we found the right place for this node in the reference or
        // not, we still need to advance along the reference path. We assume the
        // whole node (except any leading bogus characters) is included in the
        // path (since it sort of has to be, syntactically, unless it's the
        // first or last node).
        referenceBase += sequence.size();
        
        // TODO: handle leading bogus characters in calls on the first node.
    }
    
    // Create the actual reference sequence we will use
    index.sequence = refSeqStream.str();
    
    // Announce progress.
    std::cerr << "Traced " << referenceBase << " bp reference path " << refPathName << "." << std::endl;
    
    if(index.sequence.size() < 100) {
        std::cerr << "Reference sequence: " << index.sequence << std::endl;
    }
    
    // Give back the indexes we have been making
    return index;
}

void help_main(char** argv) {
    std::cerr << "usage: " << argv[0] << " [options] VGFILE GLENNFILE" << std::endl
        << "Convert a Glenn-format vg graph and variant file pair to a VCF." << std::endl
        << std::endl
        << "There are three objects in play: the reference (a single path), "
        << "the graph (containing the reference as a path) and the sample "
        << "(which is a set of calls on the graph, with some substitutions, "
        << "defined by the Glenn file)."
        << std::endl
        << "options:" << std::endl
        << "    -r, --ref PATH      use the given path name as the reference path" << std::endl
        << "    -c, --contig NAME   use the given name as the VCF contig name" << std::endl
        << "    -g, --gvcf          include lines for non-variant positions" << std::endl
        << "    -s, --sample NAME   name the sample in the VCF with the given name" << std::endl
        << "    -o, --offset INT    offset variant positions by this amount" << std::endl
        << "    -h, --help          print this help message" << std::endl;
}

#define debug
int main(int argc, char** argv) {
    
    if(argc == 1) {
        // Print the help
        help_main(argv);
        return 1;
    }
    
    // Option variables
    // What's the name of the reference path in the graph?
    std::string refPathName = "";
    // What name should we give the contig in the VCF file?
    std::string contigName = "";
    // What name should we use for the sample in the VCF file?
    std::string sampleName = "SAMPLE";
    // Should we output lines for all the reference positions that do exist?
    bool announceNonVariant = false;
    // How far should we offset positions of variants?
    int64_t variantOffset = 0;
    
    optind = 1; // Start at first real argument
    bool optionsRemaining = true;
    while(optionsRemaining) {
        static struct option longOptions[] = {
            {"ref", required_argument, 0, 'r'},
            {"contig", required_argument, 0, 'c'},
            {"gvcf", no_argument, 0, 'g'},
            {"sample", required_argument, 0, 's'},
            {"offset", required_argument, 0, 'o'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int optionIndex = 0;

        char option = getopt_long(argc, argv, "r:c:gs:o:h", longOptions, &optionIndex);
        switch(option) {
        // Option value is in global optarg
        case 'r':
            // Set the reference path name
            refPathName = optarg;
            break;
        case 'c':
            // Set the contig name
            contigName = optarg;
            break;
        case 'g':
            // Say we need to announce non-variant ref positions
            announceNonVariant = true;
            break;
        case 's':
            // Set the sample name
            sampleName = optarg;
            break;
        case 'o':
            // Offset variants
            variantOffset = std::stoll(optarg);
            break;
        case -1:
            optionsRemaining = false;
            break;
        case 'h': // When the user asks for help
        case '?': // When we get options we can't parse
            help_main(argv);
            exit(1);
            break;
        default:
            std::cerr << "Illegal option: " << option << std::endl;
            exit(1);
        }
    }
    
    if(argc - optind < 2) {
        // We don't have two positional arguments
        // Print the help
        help_main(argv);
        return 1;
    }
    
    // Pull out the file names
    std::string vgFile = argv[optind++];
    std::string glennFile = argv[optind++];
    
    // Open the vg file
    std::ifstream vgStream(vgFile);
    if(!vgStream.good()) {
        std::cerr << "Could not read " << vgFile << std::endl;
        exit(1);
    }
    
    // Load up the VG file
    vg::VG vg(vgStream);
    
    vg.paths.sort_by_mapping_rank();
    vg.paths.rebuild_mapping_aux();
    
    if(refPathName.empty()) {
        std:cerr << "Graph has " << vg.paths.size() << " paths to choose from."
            << std::endl;
        if(vg.paths.size() == 1) {
            // Autodetect the reference path name as the name of the only path
            refPathName = (*vg.paths._paths.begin()).first;
        } else {
            refPathName = "ref";
        }
        
        std::cerr << "Guessed reference path name of " << refPathName
            << std::endl;
    }
    
    // Follow the reference path and extract indexes we need: index by node ID,
    // index by node start, and the reconstructed path sequence.
    ReferenceIndex index = trace_reference_path(vg, refPathName);
    
    // Open up the Glenn-file
    std::ifstream glennStream(glennFile);
    
    // Parse it into an internal format, where we track status and copy number
    // for nodes and edges.
    
    // This holds copy numbers for all the nodes we have copy numbers called
    // for, by the node pointer in the vg graph.
    std::map<vg::Node*, size_t> nodeCopyNumbers;
    // This holds all the edges that are deletions, by the pointer to the stored
    // Edge object in the VG graph
    std::set<vg::Edge*> deletionEdges;
    
    // Loop through all the lines
    std::string line;
    while(std::getline(glennStream, line)) {
        // For each line
        
        if(line == "") {
            // Skip blank lines
            continue;
        }
        
        // Make a stringstream to read out tokens
        std::stringstream tokens(line);
        
        // Read the kind of line this is ("N"ode or "E"dge)
        std::string lineType;
        tokens >> lineType; 
        
        if(lineType == "N") {
            // This is a line about a node
            
            // Read the node ID
            int64_t nodeId;
            tokens >> nodeId;
            
            if(!vg.has_node(nodeId)) {
                throw std::runtime_error("Invalid node: " + std::to_string(nodeId));
            }
            
            // What kind of call is it? We only care that it isn't "U"ncalled
            std::string callType;
            tokens >> callType;
            
            if(callType == "U") {
                // This node has no called copy number
#ifdef debug
                std::cerr << "Uncalled node: " << nodeId << endl;
#endif
                // Handle the next line
                continue;
            }
            
            // Read the copy number
            size_t copyNumber;
            tokens >> copyNumber;
            
            assert(0 <= copyNumber);
            assert(copyNumber <= 2);
            
#ifdef debug
            std::cerr << "Node " << nodeId << " has copy number " << copyNumber << endl;
#endif
            
            // Save it
            nodeCopyNumbers[vg.get_node(nodeId)] = copyNumber;
            
        } else if(lineType == "E") {
        
            // Read the edge data
            std::string edgeDescription;
            tokens >> edgeDescription;
            
            // Split on commas. We'd just iterate the regex iterator ourselves,
            // but it seems to only split on the first comma if we do that.
            std::vector<string> parts;
            std::copy(std::sregex_token_iterator(edgeDescription.begin(), edgeDescription.end(), std::regex(","), -1), std::sregex_token_iterator(), std::back_inserter(parts));
            
            // We need the four fields to describe an edge.
            assert(parts.size() == 4);
            
            // Parse the from node
            int64_t from = std::stoll(parts[0]);
            // And the from_start flag
            bool fromStart = std::stoi(parts[1]);
            // Make a NodeSide for the from side
            vg::NodeSide fromSide(from, !fromStart);
            // Parse the to node
            int64_t to = std::stoll(parts[2]);
            // And the to_end flag
            bool toEnd = std::stoi(parts[3]);
            // Make a NodeSide for the to side
            vg::NodeSide toSide(to, toEnd);
            
            if(!vg.has_edge(std::make_pair(fromSide, toSide))) {
                // Ensure we really have that edge
                throw std::runtime_error("Edge " + edgeDescription + " not in graph.");
            }
            
            // Parse the mode
            std::string mode;
            tokens >> mode;
            
            if(mode == "L") {
                // This is a deletion edge
#ifdef debug
                std::cerr << "Edge " << edgeDescription << " is a deletion." << endl;
#endif
                
                // Say it's a deletion
                deletionEdges.insert(vg.get_edge(std::make_pair(fromSide, toSide)));

            }
        
        } else {
            // This is not a real kind of line
            throw std::runtime_error("Unknown line type: " + lineType);
        }
        
        
        
        
    }
    
    
    // Generate a vcf header. We can't make Variant records without a
    // VariantCallFile, because the variants need to know which of their
    // available info fields or whatever are defined in the file's header, so
    // they know what to output.
    std::stringstream headerStream;
    write_vcf_header(headerStream, sampleName, contigName, index.sequence.size() + variantOffset);
    
    // Load the headers into a new VCF file object
    vcflib::VariantCallFile vcf;
    std::string headerString = headerStream.str();
    assert(vcf.openForOutput(headerString));
    
    // Spit out the header
    std::cout << headerStream.str();
    
    // Then go through it from the graph's point of view: first over alt nodes
    // backending into the reference (creating things occupying ranges to which
    // we can attribute copy number) and then over reference nodes.
    
    // We need to track the bases lost.
    size_t basesLost = 0;
    
    // TODO: look at every deletion edge and spit out variants for them.
    // Complain and maybe remember if they don't connect two primary path nodes.
    
    vg.for_each_node([&](vg::Node* node) {
        // Look at every node in the graph and spit out variants for the ones
        // that are non-reference, but for which we can push copy number to the
        // reference path, greedily.
    
        // Ensure this node is nonreference
        if(index.byId.count(node->id())) {
            // Skip reference nodes
            return;
        }
        
        while(nodeCopyNumbers.at(node) > 0) {
            // We still have copy number on this node.
            
            // Find a path to the primary reference from here that can use maximal copy number
            auto found = find_bubble(vg, node, index, nodeCopyNumbers);
            
            // Break out the actual path of NodeTraversals and the max copy number we can push.
            auto& path = found.first;
            auto& canPush = found.second;
            
            if(canPush == 0) {
                // If we can't find one, complain we discarded this node's worth
                // of variation (* copy number?)
                std::cerr << "Can't account for " << nodeCopyNumbers.at(node) <<
                    " copies of node " << node->id() << " length " <<
                    node->sequence().size() << std::endl;
                
                basesLost += node->sequence().size();
                
                // Skip this node and move on to the next
                break;
            }
            
            // Turn it into a substitution/insertion
            
            // The position we have stored for this start node is the first
            // position along the reference at which it occurs. Our bubble
            // goes forward in the reference, so we must come out of the
            // opposite end of the node from the one we have stored.
            auto referenceIntervalStart = index.byId.at(path.front().node->id()).first +
                path.front().node->sequence().size();
            
            // The position we have stored for the end node is the first
            // position it occurs in the reference, and we know we go into
            // it in a reference-concordant direction, so we must have our
            // past-the-end position right there.
            auto referenceIntervalPastEnd = index.byId.at(path.back().node->id()).first;
            
            
            // We'll fill in this stream with all the node sequences we visit on
            // the path, except for the first and last.
            std::stringstream altStream;
            
            
            if(referenceIntervalPastEnd - referenceIntervalStart == 0) {
                // If this is an insert, make sure we have the 1 base before it.
                
                // TODO: we should handle an insert at the very beginning
                assert(referenceIntervalStart > 0);
                
                // Budge jeft and add that character to the alt as well
                referenceIntervalStart--;
                altStream << index.sequence[referenceIntervalStart];
            }
            
            // We also need a list of all the alt node IDs for anming the
            // variant.
            std::stringstream idStream;
            
            for(int64_t i = 1; i < path.size() - 1; i++) {
                // For all but the first and last nodes, grab their sequences in
                // the correct orientation.
                
                std::string addedSequence = path[i].node->sequence();
            
                if(path[i].backward) {
                    // If the node is traversed backward, we need to flip its sequence.
                    addedSequence = vg::reverse_complement(addedSequence);
                }
                
                // Stick the sequence
                altStream << addedSequence;
                
                // Debit copy number.
                nodeCopyNumbers[path[i].node] -= canPush;
                
                // Record ID
                idStream << std::to_string(path[i].node->id());
                if(i != path.size() - 2) {
                    // Add a separator (-2 since the last thing is path is an
                    // anchoring reference node)
                    idStream << "_";
                }
            }

            

            // Make the variant and emit it.
            std::string refAllele = index.sequence.substr(
                referenceIntervalStart, referenceIntervalPastEnd - referenceIntervalStart);
            std::string altAllele = altStream.str();
            
            // Make a Variant
            vcflib::Variant variant;
            variant.sequenceName = contigName;
            variant.setVariantCallFile(vcf);
            variant.quality = 0;
            variant.position = referenceIntervalStart + 1 + variantOffset;
            variant.id = idStream.str();
            
            // Initialize the ref allele
            create_ref_allele(variant, refAllele);
            
            // Add the alt allele
            int altNumber = add_alt_allele(variant, altAllele);
            
            // Say we're going to spit out the genotype for this sample.        
            variant.format.push_back("GT");
            auto& genotype = variant.samples[sampleName]["GT"];
            
            if(canPush == 1) {
                // We're allele alt and ref heterozygous.
                genotype.push_back(std::to_string(altNumber) + "/0");
            } else if(canPush == 2) {
                // We're alt homozygous, other overlapping variants notwithstanding.
                genotype.push_back(std::to_string(altNumber) + "/" + std::to_string(altNumber));
            } else {
                // We're something weird
                throw std::runtime_error("Invalid copy number for alt: " + std::to_string(canPush));
            }
            
#ifdef debug
            std::cerr << "Found variant " << refAllele << " -> " << altAllele
                << " caused by nodes " <<  variant.id
                << " at 1-based reference position " << variant.position
                << std::endl;
#endif

            if(can_write_alleles(variant)) {
                // Output the created VCF variant.
                std::cout << variant << std::endl;
            } else {
                std::cerr << "Variant is too large" << std::endl;
                // TODO: account for the 1 base we added extra if it was a pure
                // insert.
                basesLost += altAllele.size();
            }
            
            
        }
        
        
    });
    
    // Announce how much we can't show.
    std::cerr << "Had to drop " << basesLost << " bp of unrepresentable variation." << std::endl;
    
    return 0;
}



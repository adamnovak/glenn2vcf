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
 * minimal-length-in-nodes node list paths starting on the given named path and
 * ending at the given node, using only the given present nodes as non-primary-
 * path parts of the bubble.
 */
std::vector<std::list<vg::NodeTraversal>> bfs_left(vg::VG& graph,
    vg::NodeTraversal node, const std::string& pathName,
    const std::set<int64_t>& nodesPresent, int64_t maxDepth = 10) {

    // Do a BFS => Use a queue of paths to extend
    
    
    // This holds the paths to get to NodeTraversals to visit (all of which will
    // end with the node we're starting with).
    std::list<std::list<vg::NodeTraversal>> toExtend;
    
    // This keeps a set of all the oriented nodes we already got to and don't
    // need to queue again. TODO: what if one way to get to them fails due to
    // not being called as existing or something?
    std::set<vg::NodeTraversal> alreadyQueued;
    
    // Start at this node at depth 0.
    toExtend.emplace_back(std::list<vg::NodeTraversal> {node});
    alreadyQueued.insert(node);
    
    // Fill in this result list
    std::vector<std::list<vg::NodeTraversal>> toReturn;
    
    while(!toExtend.empty()) {
        // Keep going until we've visited every node up to our max search depth.
        
        // Dequeue a path to extend
        std::list<vg::NodeTraversal> path = toExtend.front();
        toExtend.pop_front();
        
        if(!toReturn.empty() && path.size() > toReturn.front().size()) {
            // We already know of a path that gets to the reference path in
            // fewer steps. Drop this one and any possible descendants.
            continue;
        }
        
        // Look up and see if the front node on the path is on our reference
        // path
        if(graph.paths.has_node_mapping(path.front().node) &&
            graph.paths.get_node_mapping(path.front().node).count(pathName) &&
            graph.paths.get_node_mapping(path.front().node).at(pathName).size()) {
            // There are paths that visit this node, and one of them is the
            // reference path.
            
            // Say we got to the right place
            toReturn.push_back(path);
            
            // Don't bother looking for extensions, we already got there.
        } else if(path.size() <= maxDepth) {
            // We haven't hit the reference path yet, but we also haven't hit
            // the max depth. Extend with all the possible extensions.
            
            // Look left
            vector<vg::NodeTraversal> prevNodes;
            graph.nodes_prev(path.front(), prevNodes);
            
            for(auto prevNode : prevNodes) {
                // For each node we can get to
                
                // Is this an anchoring node?
                bool isPrimaryPath = graph.paths.has_node_mapping(prevNode.node) &&
                    graph.paths.get_node_mapping(prevNode.node).count(pathName) &&
                    graph.paths.get_node_mapping(prevNode.node).at(pathName).size();
            
                if((nodesPresent.count(prevNode.node->id()) || isPrimaryPath) && !alreadyQueued.count(prevNode)) {
                    // This node is one we're allowed to visit (called as
                    // present or on the primary path), and we haven't already
                    // found a way to get to it.
            
                    // Make a new path extended left with each of these
                    std::list<vg::NodeTraversal> extended(path);
                    extended.push_front(prevNode);
                    toExtend.push_back(extended);
                    
                    // Remember we found a way to this node, so we don't try and
                    // visit it other ways.
                    alreadyQueued.insert(prevNode);
                }
            }
        }
        
    }
    
    // When we get here, we've filled in the vector with all the paths from the
    // primary path to our designated node of minimal node count. TODO: sort
    // them somehow.
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
 * minimal-length-in-nodes node list paths starting at the given node and ending
 * on the given named path, using only the given present nodes as non-primary-
 * path parts of the bubble.
 */
std::vector<std::list<vg::NodeTraversal>> bfs_right(vg::VG& graph,
    vg::NodeTraversal node, const std::string& pathName,
    const std::set<int64_t>& nodesPresent, int64_t maxDepth = 10) {

    // Look left from the backward version of the node.
    std::vector<std::list<vg::NodeTraversal>> toReturn = bfs_left(graph, flip(node), pathName, nodesPresent, maxDepth);
    
    for(auto& path : toReturn) {
        // Invert the order of every path in palce
        path.reverse();
        
        for(auto& traversal : path) {
            // And invert the orientation of every node in the path in place.
            traversal = flip(traversal);
        }
    }
    
    return toReturn;
}

/**
 * Given a vg graph, a node in the graph, and a path in the graph, look out from
 * the node in both directions to find a bubble relative to the path, with a
 * consistent orientation. Return the ordered and oriented nodes in the bubble,
 * with the outer nodes being oriented forward along the named path, and with
 * the first node coming before the last node in the reference. Needs a set of
 * nodes identified as entirely present, so it can search paths using only those
 * nodes.
 */
std::vector<vg::NodeTraversal>
find_bubble(vg::VG& graph, vg::Node* node, const std::string& pathName,
const std::set<int64_t>& nodesPresent) {

    // Find paths on both sides, with nodes on the primary path at the outsides
    // and this node in the middle.
    auto leftPaths = bfs_left(graph, vg::NodeTraversal(node), pathName, nodesPresent);
    auto rightPaths = bfs_right(graph, vg::NodeTraversal(node), pathName, nodesPresent);
    
    // Find a combination of two paths which gets us to the reference in a
    // consistent orientation (meaning that when you look at the ending nodes'
    // Mappings in the reference path, the ones with minimal ranks have the same
    // orientations).

    for(auto leftPath : leftPaths) {
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
        
        // Get all the Mappings in a nonempty set
        std::set<vg::Mapping*> leftMappings = graph.paths.get_node_mapping(leftNode).at(pathName);
        
        vg::Mapping* firstMapping = nullptr;
        
        for(vg::Mapping* mapping : leftMappings) {
            // Look at all the mappings of this node
            if(firstMapping == nullptr || firstMapping->rank() > mapping->rank()) {
                // This one is the leftmost one we have seen so far
                firstMapping = mapping;
            }
        }
        
        assert(firstMapping != nullptr);
        
        // We have a backward orientation relative to the reference path if we
        // were traversing the anchoring node backwards, xor if it is backwards
        // in the reference path.
        bool leftRelativeOrientation = leftOrientation != firstMapping->is_reverse();
        
        for(auto rightPath : rightPaths) {
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
            
            // Get all the Mappings in a nonempty set
            std::set<vg::Mapping*> rightMappings = graph.paths.get_node_mapping(rightNode).at(pathName);
            
            vg::Mapping* lastMapping = nullptr;
            
            for(vg::Mapping* mapping : rightMappings) {
                // Look at all the mappings of this node
                if(lastMapping == nullptr || lastMapping->rank() > mapping->rank()) {
                    // This one is the leftmost one we have seen so far.
                    // We do infact want the leftmost mapping we can find on the right.
                    lastMapping = mapping;
                }
            }
            
            assert(lastMapping != nullptr);
            
            // We have a backward orientation relative to the reference path if we
            // were traversing the anchoring node backwards, xor if it is backwards
            // in the reference path.
            bool rightRelativeOrientation = rightOrientation != lastMapping->is_reverse();
            
            if(leftRelativeOrientation == rightRelativeOrientation &&
                ((!leftRelativeOrientation && firstMapping->rank() < lastMapping->rank()) ||
                (leftRelativeOrientation && firstMapping->rank() > lastMapping->rank()))) {
                // We found a pair of paths that get us to and from the
                // reference without turning around, and that don't go back to
                // the reference before they leave.
                
                // Start with the left path
                std::vector<vg::NodeTraversal> toReturn{leftPath.begin(), leftPath.end()};
                
                for(auto it = ++(rightPath.begin()); it != rightPath.end(); ++it) {
                    // For all but the first node on the right path, add that in
                    toReturn.push_back(*it);
                }
                
                if(leftRelativeOrientation) {
                    // Turns out our anchored path is backwards.
                    
                    // Reorder everything the other way
                    std::reverse(toReturn.begin(), toReturn.end());
                    
                    for(auto& traversal : toReturn) {
                        // Flip each traversal
                        traversal = flip(traversal);
                    }
                }
                
                // Just give the first valid path we find. TODO: search for
                // paths where the whole thing is annotated present or whatever.
#ifdef debug        
                std::cerr << "Merged path:" << std::endl;
                for(auto traversal : toReturn) {
                    std::cerr << "\t" << traversal << std::endl;
                }
#endif
                return toReturn;
            }
            
        }
    }
    
    // Return the empty vector if we can't find compatible paths.
    return std::vector<vg::NodeTraversal>();
}

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
    std::map<vg::Node*, int> nodeCopyNumbers;
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
            int copyNumber;
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
        
        
    });
    
    // Announce how much we can't show.
    std::cerr << "Had to drop " << basesLost << " bp of unrepresentable variation." << std::endl;
    
    return 0;
}



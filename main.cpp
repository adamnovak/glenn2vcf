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
 * Represents our opinion of a particular base in a node in the graph.
 */
struct BaseCall {
    // How many alts are allowed?
    static const int MAX_ALTS = 2;

    // Is the default base here peresnt?
    bool graphBasePresent = false;
    // How many alts are here?
    char numberOfAlts = 0;
    // What are the actual alt bases?
    char alts[MAX_ALTS];
    
    /**
     * Create a new BaseCall representing what's going on at this position in
     * the graph. Uses the calls from the Glenn file (up to two one-character
     * strings). This constructor is responsible for interpreting the "-" and
     * "." special call characters.
     */
    BaseCall(const std::set<std::string>& altSet) {
        
        // We start with no alts.
        numberOfAlts = 0;
        // And with no indication that this base is present in the graph.
        graphBasePresent = false;
        for(auto alt : altSet) {
            if(alt == "-") {
                // This isn't a real alt base. It just means "same as the other
                // character". Skip it.
                continue;
            } else if(alt == ".") {
                // The occurrence of this character means that the graph's
                // normal base is actually present.
                graphBasePresent = true;
                continue;
            }
            // Otherwise we got a real letter.
            // Make sure we're not going out of bounds
            assert(alt.size() == 1);
            assert(numberOfAlts < MAX_ALTS);
            // Store each alt at the next position in  the array.
            alts[numberOfAlts] = alt[0];
            numberOfAlts++;
        }
    }
    
    /**
     * Default constructor for being put in a vector.
     */
    BaseCall() {
        // Nothing to do!
    }
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
    for(int i = 0; i < variant.alleles.size(); i++) {
        if(variant.alleles[i] == allele) {
            // Already exists
            return i;
        }
    }
    
    // Copy the allele so we can throw out bad characters
    std::string fixed(allele);
    
    for(size_t i = 0; i < fixed.size(); i++) {
        // Look at all the bases
        if(fixed[i] != 'A' && fixed[i] != 'C' && fixed[i] != 'G' && fixed[i] != 'T') {
            // Correct anything bogus (like "X") to N
            fixed[i] = 'N';
        }
    }

    // Add it as an alt
    variant.alt.push_back(allele);
    // Make it next in the alleles-by-index list
    variant.alleles.push_back(allele);
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
 * ending at the given node, using only the given present nodes.
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
                if(nodesPresent.count(prevNode.node->id()) && !alreadyQueued.count(prevNode)) {
                    // This node is one we're allowed to visit, and we haven't
                    // already found a way to get to it.
            
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
 * on the given named path, using only the given present nodes.
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
        << "    -d, --deletions     include reference deletions or replacements with <NON_REF>" << std::endl
        << "    -b, --bubbles       include nodes on larger inferred bubbles" << std::endl
        << "    -s, --sampe NAME    name the sample in the VCF with the given name" << std::endl
        << "    -o, --offset INT    offset variant positions by this amount" << std::endl
        << "    -h, --help          print this help message" << std::endl;
}

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
    // Should we include deletions from the reference when the reference base
    // isn't noted as present?
    bool announceDeletions = false;
    // Should arbitrarily find long expanded bubbles for non-anchored nodes?
    bool useBubbles = false;
    // How far should we offset positions of variants?
    int64_t variantOffset = 0;
    
    optind = 1; // Start at first real argument
    bool optionsRemaining = true;
    while(optionsRemaining) {
        static struct option longOptions[] = {
            {"ref", required_argument, 0, 'r'},
            {"contig", required_argument, 0, 'c'},
            {"gvcf", no_argument, 0, 'g'},
            {"deletions", no_argument, 0, 'd'},
            {"bubbles", no_argument, 0, 'b'},
            {"sample", required_argument, 0, 's'},
            {"offset", required_argument, 0, 'o'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int optionIndex = 0;

        char option = getopt_long(argc, argv, "r:c:gdbs:o:h", longOptions, &optionIndex);
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
        case 'd':
            // Say we need to announce deleted ref positions
            announceDeletions = true;
            break;
         case 'b':
            // Say we should admit inferred large bubbles
            useBubbles = true;
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
    
    // Make sure the reference path is present
    assert(vg.paths.has_path(refPathName));
    
    // Trace the reference path, and assign each node a canonical reference
    // range. The first base of the node occurs at the given position in the
    // reference. Some nodes may be backward (orientation true) at their
    // canonical reference positions. In this case, the last base of the node
    // occurs at the given position.
    std::map<int64_t, std::pair<size_t, bool>> referencePositionAndOrientation;
    
    // We're also going to build the reference sequence string
    std::stringstream refSeqStream;
    
    // We also need to be able to map any node start position in the reference
    // to the NodeTraversal that lives there.
    std::map<size_t, vg::NodeTraversal> nodesByReference;
    
    // We're going to need to note which reference nodes have had copy number
    // taken from them by being bypassed by a present alt. We're not going to do
    // it in a way that's clever enough to resolve multiple reference
    // occurrences smartly, but we're going to do it. This maps from reference
    // node ID to total copy number of that node used by bypassing alts.
    // We also use this to note alts that have had copy number used by other
    // bubbles.
    std::map<int64_t, int64_t> copynumberUsedByAlts;
    
    // What base are we at in the reference
    size_t referenceBase = 0;
    
    // What was the last rank? Ranks must always go up.
    int64_t lastRank = -1;
    
    for(auto mapping : vg.paths.get_path(refPathName)) {
        // All the mappings need to be perfect matches.
        assert(mapping_is_perfect_match(mapping));
    
        if(!referencePositionAndOrientation.count(mapping.position().node_id())) {
            // This is the first time we have visited this node in the reference
            // path.
            
            // Add in a mapping.
            referencePositionAndOrientation[mapping.position().node_id()] = 
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
        nodesByReference[referenceBase] = vg::NodeTraversal(
            vg.get_node(mapping.position().node_id()), mapping.is_reverse()); 
            
        // All reference nodes start out with no copy number used by alts.
        copynumberUsedByAlts[mapping.position().node_id()] = 0;
        
        // Whether we found the right place for this node in the reference or
        // not, we still need to advance along the reference path. We assume the
        // whole node (except any leading bogus characters) is included in the
        // path (since it sort of has to be, syntactically, unless it's the
        // first or last node).
        referenceBase += sequence.size();
        
        // TODO: handle leading bogus characters in calls on the first node.
    }
    
    // Create the actual reference sequence we will use
    std::string refSeq(refSeqStream.str());
    
    // Announce progress.
    std::cerr << "Traced " << referenceBase << " bp reference path " << refPathName << "." << std::endl;
    
    if(refSeq.size() < 100) {
        std::cerr << "Reference sequence: " << refSeq << std::endl;
    }
    
    // Open up the Glenn-file
    std::ifstream glennStream(glennFile);
    
    // Parse it into an internal format, where we keep track of whether bases
    // exist or not.
    // Stores call info for a position by graph node and index in the node.    
    std::map<int64_t, std::vector<BaseCall>> callsByNodeOffset;
    
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
        
        // Read the node id
        int64_t nodeId;
        tokens >> nodeId;
        
        // Read the offset
        size_t offset;
        tokens >> offset;
        offset--; // Make it 0-based
        
        // Read the base that the graph has at this position
        char graphBase;
        tokens >> graphBase;
        
        if(vg.get_node(nodeId)->sequence().size() <= offset) {
            throw std::runtime_error("Node " + std::to_string(nodeId) + " offset " +
                std::to_string(offset) + " is out of range on " +
                std::to_string(vg.get_node(nodeId)->sequence().size()) + " bp node");
        }
        
        // TODO: make sure the graph we're using agrees with the reference base.
        if(vg.get_node(nodeId)->sequence()[offset] != graphBase) {
            throw std::runtime_error("Node " + std::to_string(nodeId) + " offset " +
                std::to_string(offset) + " is " + char_to_string(vg.get_node(nodeId)->sequence()[offset]) +
                " and not " + graphBase);
        }
        
        // Read the call string
        std::string call;
        tokens >> call;
        
        // Split that out into a set of call character strings. This is how we
        // split on commas in C++. Ask for the non-matched parts of a regex
        // iterator that matches commas.
        std::set<std::string> callCharacters(std::sregex_token_iterator(
            call.begin(), call.end(), std::regex(","), -1),
            std::sregex_token_iterator());
        
        // Fill in the callsByNodeOffset map for this base of this node.
        
        if(callsByNodeOffset[nodeId].size() <= offset) {
            // Make sure there's room in the vector
            callsByNodeOffset[nodeId].resize(offset + 1);
        }
        
        // Interpret the meaning of the -,. or A,C type character pairs that
        // Glenn is using.
        callsByNodeOffset[nodeId][offset] = BaseCall(callCharacters);

#ifdef debug
        std::cerr << "Node " << nodeId << " base " << offset << " status: "
            << (callsByNodeOffset[nodeId][offset].graphBasePresent ? "Present" : "Absent")
            << std::endl;
#endif
    }
    
    
    vg.for_each_node([&](vg::Node* node) {
        // Make sure we have all the calls we need (one per base)
        if(!callsByNodeOffset.count(node->id())) {
            // This node isn't called. Complain and compensate.
            std::cerr << "WARNING: Uncalled node: " << std::to_string(node->id()) << std::endl;
            callsByNodeOffset[node->id()].resize(node->sequence().size());
        }
        
        if(callsByNodeOffset[node->id()].size() != node->sequence().size()) {
            // All nodes are called to sufficient length
            std::cerr << "WARNING: Calls for node are too short: " << std::to_string(node->id()) << std::endl;
            callsByNodeOffset[node->id()].resize(node->sequence().size());
        }
    });
    
    // Scan through all the nodes and mark the ones that are entirely present.
    std::set<int64_t> nodesPresent;
    vg.for_each_node([&](vg::Node* node) {
        // Determine how many bases of this node are present.
        size_t basesFound = 0;
        for(auto& call : callsByNodeOffset[node->id()]) {
            // For every base in the node, note if it's graph base is present.
            basesFound += call.graphBasePresent;
        }
        
        if(basesFound == node->sequence().size()) {
            // All the bases were called present. The whole node is present.
            nodesPresent.insert(node->id());
        }
    });
    
    
    // Generate a vcf header. We can't make Variant records without a
    // VariantCallFile, because the variants need to know which of their
    // available info fields or whatever are defined in the file's header, so
    // they know what to output.
    std::stringstream headerStream;
    write_vcf_header(headerStream, sampleName, contigName, refSeq.size() + variantOffset);
    
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
    
    vg.for_each_node([&](vg::Node* node) {
        // Look at every node in the graph and spit out variants for the ones
        // that are non-reference, but attached to two reference nodes and are
        // called as present.
    
        // Ensure this node is nonreference
        if(referencePositionAndOrientation.count(node->id())) {
            // Skip reference nodes
            return;
        }
        
        // Determine how many bases of this node are present.
        size_t graphBasesPresent = 0;
        // Determine how many alt calls on the node are also present.
        int totalAltsPresent = 0;
        for(auto& call : callsByNodeOffset[node->id()]) {
            // For every base in the node, note if it's graph base is present.
            // TODO: We can't just use the node presence set because we
            // distinguish later between nodes with no bases present and nodes
            // with some bases present.
            graphBasesPresent += call.graphBasePresent;
            
            // And note how many alts are also present.
            totalAltsPresent += call.numberOfAlts;
        }
        
        if(totalAltsPresent > 0) {
            // This node is present, but it has alts we should call on it and
            // won't.
            std::cerr << "Node " << node->id() <<" has " << totalAltsPresent << " bp additional novel alts!" << std::endl;
                
            // TODO: we leave the node in, because at least one copy of it
            // exists, but we might end up calling it homozygous when really we
            // have one of it and one of a modified version of it.
            // We lose these extra alt bases no matter what happens now.
            basesLost += totalAltsPresent;
        }
        
        if(graphBasesPresent == 0) {
            // This node isn't used at all in this sample, so ignore it.
            // We lose alts if it somehow had any.
            return;
        }
        
        if(graphBasesPresent < node->sequence().size()) {
            // We shouldn't call this as a variant; they're not even
            // heterozygous this alt.
            std::cerr << "Node " << node->id() << " has only "
                << graphBasesPresent << "/" << node->sequence().size()
                << " present. Skipping!" << std::endl;
            // We lose the bases that were there.
            basesLost += graphBasesPresent;
            return;
        }
        
        // Ensure this node attaches to two reference nodes, with correct
        // orientations.
        
        // Get all the oriented nodes to the left of this one's local forward.
        vector<vg::NodeTraversal> prevNodes;
        vg.nodes_prev(vg::NodeTraversal(node), prevNodes);
        
        // Find the leftmost reference node we're attached to at our start.
        size_t leftmostInPosition = (size_t) -1;
        vg::NodeTraversal leftmostInNode;
        
        for(auto& candidate : prevNodes) {
            // Look it up in the reference
            auto found = referencePositionAndOrientation.find(candidate.node->id());
            if(found == referencePositionAndOrientation.end()) {
                // We only care about nodes that actually are in the reference
                continue;
            }
            
            if((*found).second.first < leftmostInPosition) {
                // Take this as our new leftmost way into this node from the
                // reference. We know the reference positions are strictly
                // ordered, so orientation in the reference doesn't matter here.
                leftmostInNode = candidate;
                leftmostInPosition = (*found).second.first;
            }
        }
        
        // Get all the oriented nodes to the right of this one's local forward.
        vector<vg::NodeTraversal> nextNodes;
        vg.nodes_next(vg::NodeTraversal(node), nextNodes);
        
        // Find the leftmost reference node we're attached to at our end.
        size_t leftmostOutPosition = (size_t) -1;
        vg::NodeTraversal leftmostOutNode;
        
        for(auto& candidate : nextNodes) {
            // Look it up in the reference
            auto found = referencePositionAndOrientation.find(candidate.node->id());
            if(found == referencePositionAndOrientation.end()) {
                // We only care about nodes that actually are in the reference
                continue;
            }
            
            if((*found).second.first < leftmostOutPosition) {
                // Take this as our new leftmost way into this node from the
                // reference. We know the reference positions are strictly
                // ordered, so orientation in the reference doesn't matter here.
                leftmostOutNode = candidate;
                leftmostOutPosition = (*found).second.first;
            }
        }
        
        // We're going to fill in these variables to describe a single- or
        // multi-node bubble.
        // Start of the interval we replace on the reference
        size_t referenceIntervalStart;
        // Past-the-end position of the interval we replace on the reference
        size_t referenceIntervalPastEnd;
        
        // Oriented nodes we're replacing that interval with (from which alleles
        // will be derived)
        std::vector<vg::NodeTraversal> newNodes;
        
        // Now check the above to make sure we're actually placed in a
        // consistent place in the reference. We need to be able to read along
        // the reference forward, into this node, and out the other end into the
        // reference later in the same orientation.
        if(leftmostInNode.node == nullptr || leftmostOutNode.node == nullptr) {
            if(useBubbles) {
                // Try and find a bubble that lets us place the node, and then
                // say the whole bubble exists due to the node. TODO: this can
                // have us asserting nodes we shouldn't assert when the nodes we
                // do assert would allow a different bubble.
                
                if(copynumberUsedByAlts.count(node->id()) && copynumberUsedByAlts.at(node->id()) > 0) {
                    // This node has already had some copy number accounted for
                    // in another bubble. TODO: Just assume it's genotype was
                    // correct there. There's no way to go back and assert it as
                    // homozygous when it was already called het in a bubble,
                    // for example.
                    return;
                }
                
                // Look for the bubble. We already know that the node we're starting from is present.
                std::vector<vg::NodeTraversal> bubble = find_bubble(vg, node, refPathName, nodesPresent);
                
                if(bubble.size() == 0) {
                    // We couldn't find a way to make a bubble.
                    std::cerr << "Node " << node->id() << " not part of bubble." << std::endl;
                    // We lose the bases we wanted to represent.
                    basesLost += graphBasesPresent;
                    return;
                }
                
                // Otherwise we have a bubble we can use.
                
                // TODO: ensure all the nodes we're tagging here are supposed to
                // be present, and try other bubbles if they aren't. TODO: once
                // we call the whole bubble for one of these nodes, don't call
                // it again for other nodes.
                
                // Where do we start along the reference?
                
                // The position we have stored for this start node is the first
                // position along the reference at which it occurs. Our bubble
                // goes forward in the reference, so we must come out of the
                // opposite end of the node from the one we have stored.
                referenceIntervalStart = referencePositionAndOrientation.at(bubble.front().node->id()).first +
                    bubble.front().node->sequence().size();
                
                // The position we have stored for the end node is the first
                // position it occurs in the reference, and we know we go into
                // it in a reference-concordant direction, so we must have our
                // past-the-end position right there.
                referenceIntervalPastEnd = referencePositionAndOrientation.at(bubble.back().node->id()).first;
                
                for(size_t i = 1; i < bubble.size() - 1; i++) {
                    // Stick in all the non-anchoring nodes, dropping the first and last anchoring nodes.
                    newNodes.push_back(bubble[i]);
                }
                
            } else {
                // We're missing a reference node on one side, and we can't do longer bubbles.
                std::cerr << "Node " << node->id() << " not anchored to reference." << std::endl;
                // We lose the bases we wanted to represent.
                basesLost += graphBasesPresent;
                return;
            }
        } else {
            // We can proceed with a single-node bubble.
        
            // Determine if we read into this node forward along the reference
            // (true) or backward along the reference (false). If we found the node
            // to our left in the same orientation as it occurs in the reference,
            // then we do read in forward.
            bool readInForward = leftmostInNode.backward == referencePositionAndOrientation.at(leftmostInNode.node->id()).second;
            
            // If we found the node to our right in the same orientation as it
            // occurs in the reference, then we do read out forward as well.
            bool readOutForward = leftmostOutNode.backward == referencePositionAndOrientation.at(leftmostOutNode.node->id()).second;
            
            if(readInForward != readOutForward) {
                // Going through this node would cause us to invert the direction
                // we're traversing the reference in.
                std::cerr << "Node " << node->id() << " inverts reference path." << std::endl;
                // We lose the bases we wanted to represent.
                basesLost += graphBasesPresent;
                return;
            }
            
            // We need to work out what orientation we have relative to the
            // reference.
            vg::NodeTraversal altNode(node);
            
            if(!readInForward) {
                // We have a consistent orientation, but it's backward!
                // Swap the in and out nodes, and traverse our node in reverse.
                altNode.backward = true;
                std::swap(leftmostInNode, leftmostOutNode);
            }
            
            // Now we know that the in node really is where we come into the alt,
            // and the out node really is where we leave the alt, when reading along
            // the reference path. Either may still be backward in the reference
            // path, though.
            
            // Work out where and how they are positioned in the reference
            auto& inNodePlacement = referencePositionAndOrientation.at(leftmostInNode.node->id());
            auto& outNodePlacement = referencePositionAndOrientation.at(leftmostOutNode.node->id());
            
            if(outNodePlacement.first <= inNodePlacement.first) {
                // We're perfectly fine, orientation-wise, except we let you time
                // travel and leave before you arrived.
                std::cerr << "Node " << node->id() << " allows duplication." << std::endl;
                // We lose the bases we wanted to represent.
                basesLost += graphBasesPresent;
                return;
            }
            
            // So what are the actual bounds of the reference interval covered by
            // the node? Since the node placement positions are just the first bases
            // along the reference at which the nodes occur, we don't care about
            // orientation of the anchoring node sequences.
            referenceIntervalStart = inNodePlacement.first + leftmostInNode.node->sequence().size();
            referenceIntervalPastEnd = outNodePlacement.first;
            
            // Just one node
            newNodes.push_back(altNode);
        }
        
        
        // Make sure we got reasonable bounds whatever kind of bubble we have.
        assert(referenceIntervalPastEnd >= referenceIntervalStart);
        
        // How long is this interval in the reference?
        size_t referenceIntervalSize = referenceIntervalPastEnd - referenceIntervalStart;
        
        // Trace the reference between our in node and our out node.
        size_t refPosition = referenceIntervalStart;
        
        // We want to know if the reference path opposite us is ever called as
        // present or has a novel SNP. If so, since we're present, we know we
        // must be heterozygous here. If not, we'll call ourselves homozygous
        // here. TODO: catch conflicts between homozygous non-reference mutually
        // exclusive variants.
        // This is false by default; we assume it's missing and can be proven wrong.
        // TODO: this makes us call insertions as homozygous.
        bool refPathExists = false;
        
        while(refPosition < referenceIntervalPastEnd) {
            // While we aren't at the start of the reference node that comes
            // after this alt...
            
            // Get the node starting here in the reference.
            if(!nodesByReference.count(refPosition)) {
                std::cerr << "Missing node at reference position " << refPosition << std::endl;
                throw std::runtime_error("Reference position incorrect");
            }
            vg::Node* refNode = nodesByReference.at(refPosition).node;
            
            // We know we can iterate over the whole reference node, because it
            // must start immediatelty after the previous node ends.
            for(int i = 0; i < refNode->sequence().size(); i++) {
                // See if the reference node is ever called as absent with no
                // novel SNP alt.
                
                auto& call = callsByNodeOffset[refNode->id()][i];
                if(call.graphBasePresent || call.numberOfAlts > 0) {
                    // We found evidence the reference exists in alternation
                    // with this allele.
                    refPathExists = true;
                }
            }
            
            // Advance to the start of the next reference node
            refPosition += refNode->sequence().size();
        }
        
        // Tell the reference nodes that we're using their copy number, so if
        // they turn up missing but were bypassed by an alt they don't act
        // deleted.
        refPosition = referenceIntervalStart;
        while(refPosition < referenceIntervalPastEnd) {
            // While we aren't at the start of the reference node that comes
            // after this alt...
            
            // Get the node starting here in the reference.
            if(!nodesByReference.count(refPosition)) {
                std::cerr << "Missing node at reference position " << refPosition << std::endl;
                throw std::runtime_error("Reference position incorrect");
            }
            vg::Node* refNode = nodesByReference.at(refPosition).node;
            
            // If we saw any of the reference path existing, we use 1 copy
            // number. Otherwise we use 2. Note that this can end up using >2
            // copy number from a node if multiple alts bypass it, or if one alt
            // bypasses it twice.
            copynumberUsedByAlts[refNode->id()] += (refPathExists ? 1 : 2);
                        
            // Advance to the start of the next reference node
            refPosition += refNode->sequence().size();
        }
        
        // Make a Variant
        vcflib::Variant variant;
        variant.sequenceName = contigName;
        variant.setVariantCallFile(vcf);
        variant.quality = 0;
        
        // Pull out the string for the reference allele
        std::string refAllele = refSeq.substr(referenceIntervalStart, referenceIntervalSize);
        // And for the alt allele
        std::stringstream altAlleleStream;
        
        for(auto& addedNode : newNodes) {
            // For each node we're adding in, stick in the sequence in the
            // correct orientation
            std::string addedSequence = addedNode.node->sequence();
            
            if(addedNode.backward) {
                // If the node is traversed backward, we need to flip its sequence.
                addedSequence = vg::reverse_complement(addedSequence);
            }
            
            // Stick the sequence
            altAlleleStream << addedSequence;
            
            // Also account for use of this node's copy number, so we won't go
            // and use it to create its own bubble.
            if(!copynumberUsedByAlts.count(addedNode.node->id())) {
                // If the node has no copy number tally, make some space.
                copynumberUsedByAlts[addedNode.node->id()] = 0;
            }
            
            // Whether the reference path exists or not, we use 2 copy number
            // here between the ref and the bubble
            copynumberUsedByAlts[addedNode.node->id()] += 2;
            
            if(copynumberUsedByAlts[addedNode.node->id()] > 2) {
                // If we end up doing this multiple times for a node, complain
                // to the user, because we may not quite be correct in our
                // conversion.
                std::cerr << "Warning: Node " << addedNode.node->id() << " used by multiple bubbles!" << std::endl;
            }
        }
        
        // Convert to a proper string
        std::string altAllele(altAlleleStream.str());
        
        if(refAllele.size() == 0) {
            // Shift everybody left by 1 base for the anchoring base that VCF
            // requires for insertions.
            assert(referenceIntervalStart != 0);
            referenceIntervalStart--;
            referenceIntervalSize++;
            // Add that base to the start of both alleles.
            refAllele = refSeq[referenceIntervalStart] + refAllele;
            altAllele = refSeq[referenceIntervalStart] + altAllele;
        }
        
        // Alt allele size can't be 0, no need to do the same shift for
        // deletions
        
        // Set the variant position. Convert to 1-based.
        variant.position = referenceIntervalStart + 1 + variantOffset;
        
        // Name it with the node numbers we're saying exist.
        std::stringstream idStream;
        for(size_t i = 0; i < newNodes.size(); i++) {
            idStream << std::to_string(newNodes[i].node->id());
            if(i != newNodes.size() - 1) {
                // Add a separator
                idStream << "_";
            }
        }
        variant.id = idStream.str();
        
        // Initialize the ref allele
        create_ref_allele(variant, refAllele);
        
        // Add the graph version
        int altNumber = add_alt_allele(variant, altAllele);
        
        // Say we're going to spit out the genotype for this sample.        
        variant.format.push_back("GT");
        auto& genotype = variant.samples[sampleName]["GT"];
        
        // Make it hom/het as appropriate
        if(refPathExists) {
            // We're allele alt and ref heterozygous.
            genotype.push_back(std::to_string(altNumber) + "/0");
        } else {
            // We're alt homozygous, other overlapping variants notwithstanding.
            genotype.push_back(std::to_string(altNumber) + "/" + std::to_string(altNumber));
        }

#ifdef debug
        std::cerr << "Found variant " << refAllele << " -> " << altAllele
            << " caused by node " << altNode.node->id()
            << " at 1-based reference position " << variant.position
            << std::endl;
#endif

        if(can_write_alleles(variant)) {
            // Output the created VCF variant.
            std::cout << variant << std::endl;
        } else {
            std::cerr << "Variant is too large" << std::endl;
            basesLost += altAllele.size();
        }
    });
    
    vg.for_each_node([&](vg::Node* node) {
        // Now we go through all the nodes on the reference path, and add in
        // SNPs on them.
        
        if(!referencePositionAndOrientation.count(node->id())) {
            // Skip non-reference nodes
            return;
        }
        
        // How many bases on this node are present?
        size_t basesPresent = 0;
        for(int i = 0; i < node->sequence().size(); i++) {
            auto& call = callsByNodeOffset[node->id()][i];
            if(call.graphBasePresent) {
                basesPresent++;
            }
        }
        
        // This gets set if there's a running deletion of several bases in a row
        // from a ref node that is otherwise present. It makes it less obvious
        // that you can't just declare a deletion without considering the
        // surrounding bases (we still might say G -> C and GATTTT -> G about
        // the same G).
        int64_t runningDelStart = -1;
        
        for(int i = 0; i < node->sequence().size(); i++) {
            // Are we continuing a running deletion?
            bool runningDelContinues = false;
        
            // Where are we in node local coordinates?
            int positionOnNode;
        
            // Work out where it is in the reference
            size_t referencePosition = referencePositionAndOrientation.at(node->id()).first + i;
            if(referencePositionAndOrientation.at(node->id()).second) {
                // We're backward in the reference, so incrementing i goes
                // towards the node's index 0
                positionOnNode = (node->sequence().size() - i - 1);
            } else {
                // We're forward in the reference, so incrementing i goes
                // towards the node's max index.
                positionOnNode = i;
            }
        
            // We need to keep track of the alts present.
            std::set<int> altNumbersPresent;
        
            // For each position along the node, grab the call there.
            auto& call = callsByNodeOffset[node->id()][positionOnNode];
            if(call.numberOfAlts > 0) {
                // At least one alt is present here.
                // Make the variant.
                
                // Make a Variant
                vcflib::Variant variant;
                variant.sequenceName = contigName;
                variant.setVariantCallFile(vcf);
                variant.quality = 0;
                
                // Grab its reference base
                std::string refAllele = char_to_string(refSeq[referencePosition]);
                // Initialize the ref allele
                create_ref_allele(variant, refAllele);
                
                // Add in alt bases, with the right orientation
                for(int j = 0; j < call.numberOfAlts; j++) {
                    std::string altAllele = char_to_string(call.alts[j]);
                    if(referencePositionAndOrientation.at(node->id()).second) {
                        // We need to flip the orientation to reference
                        // orientation
                        altAllele = vg::reverse_complement(altAllele);
                    }
                    // Add the novel SNP allele
                    altNumbersPresent.insert(add_alt_allele(variant, altAllele));
                }
                
                if(altNumbersPresent.count(0)) {
                    // We found an alt that actually matches ref. Make sure it's represented as ref instead.
                    altNumbersPresent.erase(altNumbersPresent.find(0));
                    call.graphBasePresent = true;
                    std::cerr << "ERROR: Found alt equal to reference. Correcting." << std::endl;
                    // TODO: see if we have something silly like an off by 1 error? Throw an error here?
                }
                
                // Set the variant position. Convert to 1-based.
                variant.position = referencePosition + 1 + variantOffset;
                
                // Cram the node and offset into the variant ID column, for these point variants
                variant.id = std::to_string(node->id()) + "." + std::to_string(positionOnNode);
                
                // Say we're going to spit out the genotype for this sample.        
                variant.format.push_back("GT");
                auto& genotype = variant.samples[sampleName]["GT"];
                
                // Make it hom/het as appropriate
                if(call.graphBasePresent) {
                    if(altNumbersPresent.size() > 0) {
                        // We have the ref and we also have an alt.
                        genotype.push_back("1/0");
                    } else {
                        // We have the ref and no alts distinct from it showed
                        // up. Should not happen because we already made sure
                        // call.numberOfAlts > 0
                        std::cerr << "ERROR: Found no alts when alts should exist" << std::endl;
                        genotype.push_back("0/0");
                    }
                } else if(altNumbersPresent.size() == 1) {
                    // We have only one alt allele, and no reference.
                    if(copynumberUsedByAlts[node->id()] == 0) {
                        // We still have all our copy number, but we aren't
                        // present. Give it all to the alt.
                        genotype.push_back("1/1");
                    } else {
                        // We're in alternation with one or more alts that took
                        // our copy number. Call this base as het, where the
                        // "ref" side might already have been replaced.
                        // TODO: should we make haploid calls here?
                        genotype.push_back("1/0");
                    }
                } else if(altNumbersPresent.size() == 2) {
                    // We have two alt alleles and no reference. We must have
                    // both present.
                    genotype.push_back("1/2");
                } else {
                    // This should never happen
                    throw std::runtime_error("Semantically invalid BaseCall");
                }

#ifdef debug
                std::cerr << "Found  alt-bearing variant " << refAllele << " -> ";
                for(int j = 0; j < call.numberOfAlts; j++) {
                    std::cerr << call.alts[j] << ",";
                }
                std::cerr << " on node " << node->id()
                    << " at 1-based reference position " << variant.position
                    << std::endl;
#endif
                    
                if(can_write_alleles(variant)) {
                    // Output the created VCF variant.
                    std::cout << variant << std::endl;
                } else {
                    std::cerr << "Variant is too large" << std::endl;
                    basesLost += altNumbersPresent.size();
                }
            } else if(!call.graphBasePresent) {
                // This reference base isn't present at all!
                if(copynumberUsedByAlts[node->id()] < 2 && announceDeletions) {
                    // All the copy number we would expect to see hasn't been
                    // used up by alts bypassing this base, and we want to hear
                    // about deletions. We need to say it's missing.
                    
                    if(basesPresent == 0) {
                        // The whole node is missing! And we don't know where it
                        // went! Give a generic NON_REF GVCF answer.
                        
                        
                        // We need to work out what the ref (present) and
                        // alt (missing) alleles will be.
                        std::string refAllele = refSeq.substr(referencePosition, 1);
                        std::string altAllele = "<NON_REF>";
                        
                        // Make a Variant
                        vcflib::Variant variant;
                        variant.sequenceName = contigName;
                        variant.setVariantCallFile(vcf);
                        variant.quality = 0;
                        
                        // Cram the node and offset into the variant ID column, for these point variants
                        variant.id = std::to_string(node->id()) + "." + std::to_string(positionOnNode);
                        
                        // Initialize the ref allele
                        create_ref_allele(variant, refAllele);
                        
                        // Add the novel deletion allele
                        int altNumber = add_alt_allele(variant, altAllele);
                        
                        // Say it's homozygous alt
                        variant.format.push_back("GT");
                        auto& genotype = variant.samples[sampleName]["GT"];
                        genotype.push_back(std::to_string(altNumber) + "/" + std::to_string(altNumber));
                        
                        // Set the variant position. Convert to 1-based.
                        variant.position = referencePosition + 1 + variantOffset;
                        
#ifdef debug
                        std::cerr << "Found NR variant " << refAllele << " -> "
                            << altAllele << " on node " << node->id()
                            << " at 1-based reference position " << variant.position
                            << std::endl;
#endif
                            
                        if(can_write_alleles(variant)) {
                            // Output the created VCF variant.
                            std::cout << variant << std::endl;
                        } else {
                            std::cerr << "Variant is too large" << std::endl;
                            basesLost += refAllele.size();
                        }
                        
                    } else {
                        // One or a few bases of the node are missing. Call them
                        // as deletions.
                        
                        if(runningDelStart == -1) {
                            // Start a new deletion here. Note that here may be
                            // base 0.
                            runningDelStart = referencePosition;
                        }
                        // Otherwise we get used by the running del
                        runningDelContinues = true;
                    }
                            
                }
            } else {
                // This base is present.
                // Is that described in an alt? Or do we want to account for it now?
                if(copynumberUsedByAlts[node->id()] == 0 && announceNonVariant) {
                    // Nope. Nobody has ever heard of this reference base. We
                    // need to talk about it.
                    
                    std::string refAllele = refSeq.substr(referencePosition, 1);
                    std::string altAllele = "<NON_REF>";
                    
                    // Make a Variant
                    vcflib::Variant variant;
                    variant.sequenceName = contigName;
                    variant.setVariantCallFile(vcf);
                    variant.quality = 0;
                    
                    // Cram the node and offset into the variant ID column, for these point variants
                    variant.id = std::to_string(node->id()) + "." + std::to_string(positionOnNode);
                    
                    // Initialize the ref allele
                    create_ref_allele(variant, refAllele);
                    
                    // Add the novel deletion allele
                    int altNumber = add_alt_allele(variant, altAllele);
                    
                    // Say it's homozygous ref
                    variant.format.push_back("GT");
                    auto& genotype = variant.samples[sampleName]["GT"];
                    genotype.push_back("0/0");
                    
                    // Set the variant position. Convert to 1-based.
                    variant.position = referencePosition + 1 + variantOffset;
                    
                    if(can_write_alleles(variant)) {
                        // Output the created VCF variant.
                        std::cout << variant << std::endl;
                    } else {
                        std::cerr << "Variant is too large" << std::endl;
                        basesLost += refAllele.size();
                    }
                }
            }
            
            if(runningDelStart != -1 && !runningDelContinues) {
                // We're right after the end of a running deletion.
                
                // How many bases were deleted?
                size_t runningDelLength = referencePosition - runningDelStart;
                
#ifdef debug
                std::cerr << "Running del length " << runningDelLength << " internal to node." << std::endl;
#endif
                
                std::string refAllele;
                std::string altAllele;
                
               
                
                if(runningDelStart == 0) {
                    // Special handling of deletion of the first base
                    
                    // Grab its n reference bases
                    refAllele = refSeq.substr(runningDelStart, runningDelLength + 1);
                    
                    // Grab its one alt base
                    altAllele = refSeq.substr(referencePosition, 1);
                    
                } else {
                
                    // We want the ref allele to start a base early and include
                    // the base before the deleted bases.
                    runningDelStart -= 1;
                    
                    // Grab its two reference bases
                    refAllele = refSeq.substr(runningDelStart, runningDelLength + 1);
                    
                    // Grab its one alt base
                    altAllele = refSeq.substr(runningDelStart, 1);
                }
                
                
                // Make a Variant
                vcflib::Variant variant;
                variant.sequenceName = contigName;
                variant.setVariantCallFile(vcf);
                variant.quality = 0;
                
                // Initialize the ref allele
                create_ref_allele(variant, refAllele);
                
                // Add the novel deletion allele
                int altNumber = add_alt_allele(variant, altAllele);
                
                // Say it's homozygous deleted for the whole run.
                variant.format.push_back("GT");
                auto& genotype = variant.samples[sampleName]["GT"];
                genotype.push_back(std::to_string(altNumber) + "/" + std::to_string(altNumber));
                
                // Set the variant position. Convert to 1-based.
                variant.position = runningDelStart + 1 + variantOffset;
                
#ifdef debug
                std::cerr << "Found variant " << refAllele << " -> "
                    << altAllele << " on node " << node->id()
                    << " at 1-based reference position " << variant.position
                    << std::endl;
#endif
                    
                if(can_write_alleles(variant)) {
                    // Output the created VCF variant.
                    std::cout << variant << std::endl;
                } else {
                    std::cerr << "Variant is too large" << std::endl;
                    basesLost += refAllele.size();
                }
                    
                // Clear the running deletion
                runningDelStart = -1;
                
            }
        }
        
        // We're now at the 1-past-the-end position on the node.
        if(runningDelStart != -1) {
            // We were deleting up to the end of the node.
            
            // We know we didn't delete the whole node, no need to check for pos
            // 0 in ref.
            
            // How many bases were deleted? It's the end ref position of the
            // node minus the ref position where the del started.
            size_t runningDelLength = (referencePositionAndOrientation.at(node->id()).first +
                node->sequence().size() - runningDelStart);
                
#ifdef debug
            std::cerr << "Running del length " << runningDelLength << " at end of node." << std::endl;
#endif
            
            // We want the ref allele to start a base early and include
            // the base before the deleted bases.
            runningDelStart -= 1;
            
            // Grab its two reference bases
            std::string refAllele = refSeq.substr(runningDelStart, runningDelLength + 1);
            
            // Grab its one alt base
            std::string altAllele = refSeq.substr(runningDelStart, 1);
            
            // Make a Variant
            vcflib::Variant variant;
            variant.sequenceName = contigName;
            variant.setVariantCallFile(vcf);
            variant.quality = 0;
            
            // Initialize the ref allele
            create_ref_allele(variant, refAllele);
            
            // Add the novel deletion allele
            int altNumber = add_alt_allele(variant, altAllele);
            
            // Say it's homozygous deleted for the whole run.
            variant.format.push_back("GT");
            auto& genotype = variant.samples[sampleName]["GT"];
            genotype.push_back(std::to_string(altNumber) + "/" + std::to_string(altNumber));
            
            // Set the variant position. Convert to 1-based.
            variant.position = runningDelStart + 1 + variantOffset;
            
#ifdef debug
            std::cerr << "Found variant " << refAllele << " -> "
                << altAllele << " on node " << node->id()
                << " at 1-based reference position " << variant.position
                << std::endl;
#endif
                
            if(can_write_alleles(variant)) {
                // Output the created VCF variant.
                std::cout << variant << std::endl;
            } else {
                std::cerr << "Variant is too large" << std::endl;
                basesLost += refAllele.size();
            }
                
            // Clear the running deletion
            runningDelStart = -1;
        }
        
    });
    
    // Announce how much we can't show.
    std::cerr << "Had to drop " << basesLost << " bp of unrepresentable variation." << std::endl;
    
    return 0;
}



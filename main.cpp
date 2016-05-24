// main.cpp: Main file for core graph merger

#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <set>
#include <utility>
#include <algorithm>
#include <getopt.h>

#include "ekg/vg/src/vg.hpp"
#include "ekg/vg/src/index.hpp"
#include "ekg/vg/deps/vcflib/src/Variant.h"

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
    stream << "##INFO=<ID=XREF,Number=0,Type=Flag,Description=\"Present in original graph\">" << std::endl;
    stream << "##INFO=<ID=XSEE,Number=.,Type=String,Description=\"Original graph node:offset cross-references\">" << std::endl;
    stream << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">" << std::endl;
    stream << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << std::endl;
    stream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
    stream << "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">" << std::endl;
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
 * Return true if a variant may be output, or false if this variant is valid but
 * the GATK might choke on it.
 *
 * Mostly used to throw out variants with very long alleles, because GATK has an
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
 * Get the length of a path through nodes, in base pairs.
 */
size_t bp_length(const std::list<vg::NodeTraversal>& path) {
    size_t length = 0;
    for(auto& traversal : path) {
        // Sum up length of each node's sequence
        length += traversal.node->sequence().size();
    }
    return length;
}

/**
 * Do a breadth-first search left from the given node traversal, and return
 * lengths and paths starting at the given node and ending on the indexed
 * reference path.
 */
std::set<std::pair<size_t, std::list<vg::NodeTraversal>>> bfs_left(vg::VG& graph,
    vg::NodeTraversal node, const ReferenceIndex& index, int64_t maxDepth = 10,
    bool stopIfVisited = false) {

    // Holds partial paths we want to return, with their lengths in bp.
    std::set<std::pair<size_t, std::list<vg::NodeTraversal>>> toReturn;
    
    // Do a BFS
    
    // This holds the paths to get to NodeTraversals to visit (all of which will
    // end with the node we're starting with).
    std::list<std::list<vg::NodeTraversal>> toExtend;
    
    // This keeps a set of all the oriented nodes we already got to and don't
    // need to queue again.
    std::set<vg::NodeTraversal> alreadyQueued;
    
    // Start at this node at depth 0
    toExtend.emplace_back(std::list<vg::NodeTraversal> {node});
    // Mark this traversal as already queued
    alreadyQueued.insert(node);
    
    // How many ticks have we spent searching?
    size_t searchTicks = 0;
    // Track how many options we have because size may be O(n).
    size_t stillToExtend = toExtend.size();
    
    while(!toExtend.empty()) {
        // Keep going until we've visited every node up to our max search depth.
        
#ifdef debug
        searchTicks++;
        if(searchTicks % 100 == 0) {
            // Report on how much searching we are doing.
            std::cerr << "Search tick " << searchTicks << ", " << stillToExtend << " options." << std::endl;
        }
#endif
        
        // Dequeue a path to extend.
        // Make sure to move out of the list to avoid a useless copy.
        std::list<vg::NodeTraversal> path(std::move(toExtend.front()));
        toExtend.pop_front();
        stillToExtend--;
        
        // We can't just throw out longer paths, because shorter paths may need
        // to visit a node twice (in opposite orientations) and thus might get
        // rejected later. Or they might overlap with paths on the other side.
        
        // Look up and see if the front node on the path is on our reference
        // path
        if(index.byId.count(path.front().node->id())) {
            // This node is on the reference path. TODO: we don't care if it
            // lands in a place that is itself deleted.
            
            // Say we got to the right place
            toReturn.emplace(bp_length(path), std::move(path));
            
            // Don't bother looking for extensions, we already got there.
        } else if(path.size() <= maxDepth) {
            // We haven't hit the reference path yet, but we also haven't hit
            // the max depth. Extend with all the possible extensions.
            
            // Look left
            vector<vg::NodeTraversal> prevNodes;
            graph.nodes_prev(path.front(), prevNodes);
            
            for(auto prevNode : prevNodes) {
                // For each node we can get to
                
                if(stopIfVisited && alreadyQueued.count(prevNode)) {
                    // We already have a way to get here.
                    break;
                }
            
                // Make a new path extended left with the node
                std::list<vg::NodeTraversal> extended(path);
                extended.push_front(prevNode);
                toExtend.emplace_back(std::move(extended));
                stillToExtend++;
                
                // Remember we found a way to this node, so we don't try and
                // visit it other ways.
                alreadyQueued.insert(prevNode);
            }
        }
        
    }
    
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
 * lengths and paths starting at the given node and ending on the indexed
 * reference path.
 */
std::set<std::pair<size_t, std::list<vg::NodeTraversal>>> bfs_right(vg::VG& graph,
    vg::NodeTraversal node, const ReferenceIndex& index, int64_t maxDepth = 10,
    bool stopIfVisited = false) {

    // Look left from the backward version of the node.
    auto toConvert = bfs_left(graph, flip(node), index, maxDepth, stopIfVisited);
    
    // Since we can't modify set records in place, we need to do a copy
    std::set<std::pair<size_t, std::list<vg::NodeTraversal>>> toReturn;
    
    for(auto lengthAndPath : toConvert) {
        // Flip every path to run the other way
        lengthAndPath.second.reverse();
        for(auto& traversal : lengthAndPath.second) {
            // And invert the orientation of every node in the path in place.
            traversal = flip(traversal);
        }
        // Stick it in the new set
        toReturn.emplace(std::move(lengthAndPath));
    }
    
    return toReturn;
}

/**
 * Given a vg graph, a node in the graph, and an index for the reference path,
 * look out from the node in both directions to find a shortest bubble relative
 * to the path, with a consistent orientation. The bubble may not visit the same
 * node twice.
 *
 * Takes a max depth for the searches producing the paths on each side.
 * 
 * Return the ordered and oriented nodes in the bubble, with the outer nodes
 * being oriented forward along the named path, and with the first node coming
 * before the last node in the reference.
 */
std::vector<vg::NodeTraversal>
find_bubble(vg::VG& graph, vg::Node* node, const ReferenceIndex& index,
    int64_t maxDepth = 10) {

    // Find paths on both sides, with nodes on the primary path at the outsides
    // and this node in the middle. Returns path lengths and paths in pairs in a
    // set.
    auto leftPaths = bfs_left(graph, vg::NodeTraversal(node), index, maxDepth);
    auto rightPaths = bfs_right(graph, vg::NodeTraversal(node), index, maxDepth);
    
    // Find a combination of two paths which gets us to the reference in a
    // consistent orientation (meaning that when you look at the ending nodes'
    // Mappings in the reference path, the ones with minimal ranks have the same
    // orientations) and which doesn't use the same nodes on both sides.
    
    // We need to look in different combinations of lists.
    auto testCombinations = [&](const std::list<std::list<vg::NodeTraversal>>& leftList,
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
                    
                    // Just give the first valid path we find.
#ifdef debug        
                    std::cerr << "Merged path:" << std::endl;
                    for(auto traversal : fullPath) {
                        std::cerr << "\t" << traversal << std::endl;
                    }
#endif
                    return fullPath;
                }
                
            }
        }
        
        // Return the empty path if we can't find anything.
        return std::vector<vg::NodeTraversal>();
        
    };
    
    // Convert sets to lists, which requires a copy again...
    std::list<std::list<vg::NodeTraversal>> leftConverted;
    for(auto lengthAndPath : leftPaths) {
        leftConverted.emplace_back(std::move(lengthAndPath.second));
    }
    std::list<std::list<vg::NodeTraversal>> rightConverted;
    for(auto lengthAndPath : rightPaths) {
        rightConverted.emplace_back(std::move(lengthAndPath.second));
    }
    
    // Look for a valid combination, or return an empty path if one iesn't
    // found.
    return testCombinations(leftConverted, rightConverted);
    
}


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
                std::make_pair(referenceBase, mapping.position().is_reverse());
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
        
        if(mapping.position().is_reverse()) {
            // Put the reverse sequence in the reference path
            refSeqStream << vg::reverse_complement(sequence);
        } else {
            // Put the forward sequence in the reference path
            refSeqStream << sequence;
        }
            
        // Say that this node appears here along the reference in this
        // orientation.
        index.byStart[referenceBase] = vg::NodeTraversal(
            vg.get_node(mapping.position().node_id()), mapping.position().is_reverse()); 
            
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

/**
 * Given a collection of pileups by original node ID, and a set of original node
 * id:offset cross-references in both ref and alt categories, produce a VCF
 * comment line giving the pileup for each of those positions on those nodes.
 * Includes a trailing newline if nonempty.
 *
 * TODO: VCF comments aren't really a thing.
 */
std::string get_pileup_line(const std::map<int64_t, vg::NodePileup>& nodePileups,
    const std::set<std::pair<int64_t, size_t>>& refCrossreferences,
    const std::set<std::pair<int64_t, size_t>>& altCrossreferences) {
    // We'll make a stringstream to write to.
    std::stringstream out;
    
    out << "#";
    
    for(const auto& xref : refCrossreferences) {
        // For every cross-reference
        if(nodePileups.count(xref.first) && nodePileups.at(xref.first).base_pileup_size() > xref.second) {
            // If we have that base pileup, grab it
            auto basePileup = nodePileups.at(xref.first).base_pileup(xref.second);
            
            out << xref.first << ":" << xref.second << " (ref) " << basePileup.bases() << "\t";
        }
        // Nodes with no pileups (either no pileups were provided or they didn't
        // appear/weren't visited by reads) will not be mentioned on this line
    }
    
    for(const auto& xref : altCrossreferences) {
        // For every cross-reference
        if(nodePileups.count(xref.first) && nodePileups.at(xref.first).base_pileup_size() > xref.second) {
            // If we have that base pileup, grab it
            auto basePileup = nodePileups.at(xref.first).base_pileup(xref.second);
            
            out << xref.first << ":" << xref.second << " (alt) " << basePileup.bases() << "\t";
        }
        // Nodes with no pileups (either no pileups were provided or they didn't
        // appear/weren't visited by reads) will not be mentioned on this line
    }
    // TODO: make these nearly-identical loops a loop or a lambda or something.
    
    if(out.str().size() > 1) {
        // We actually found something. Send it out with a trailing newline
        out << std::endl;
        return out.str();
    } else {
        // Give an empty string.
        return "";
    }
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
        << "    -s, --sample NAME   name the sample in the VCF with the given name" << std::endl
        << "    -o, --offset INT    offset variant positions by this amount" << std::endl
        << "    -l, --length INT    override total sequence length" << std::endl
        << "    -d, --depth INT     maximum depth for path search (default 10 nodes)" << std::endl
        << "    -p, --pileup FILE   filename for a pileup to use to annotate variants" << std::endl
        << "    -f, --min_fraction  min fraction of average coverage at which to call" << std::endl
        << "    -b, --max_het_bias  max imbalance factor between alts to call heterozygous" << std::endl
        << "    -n, --min_count     min total supporting read count to call a variant" << std::endl
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
    // How far should we offset positions of variants?
    int64_t variantOffset = 0;
    // How many nodes should we be willing to look at on our path back to the
    // primary path? Keep in mind we need to look at all valid paths (and all
    // combinations thereof) until we find a valid pair.
    int64_t maxDepth = 10;
    // What should the total sequence length reported in the VCF header be?
    int64_t lengthOverride = -1;
    // Should we load a pileup and print out pileup info as comments after
    // variants?
    std::string pileupFilename;
    // What fraction of average coverage should be the minimum to call a variant (or a single copy)?
    // Default to 0 because vg call is still applying depth thresholding
    double minFractionForCall = 0;
    // What fraction of the reads supporting an alt are we willing to discount?
    // At 2, if twice the reads support one allele as the other, we'll call
    // homozygous instead of heterozygous. At infinity, every call will be
    // heterozygous if even one read supports each allele.
    double maxHetBias = 20;
    // What's the minimum integer number of reads that must support a call? We
    // don't necessarily want to call a SNP as het because we have a single
    // supporting read, even if there are only 10 reads on the site.
    size_t minTotalSupportForCall = 2;
    
    optind = 1; // Start at first real argument
    bool optionsRemaining = true;
    while(optionsRemaining) {
        static struct option longOptions[] = {
            {"ref", required_argument, 0, 'r'},
            {"contig", required_argument, 0, 'c'},
            {"sample", required_argument, 0, 's'},
            {"offset", required_argument, 0, 'o'},
            {"depth", required_argument, 0, 'd'},
            {"length", required_argument, 0, 'l'},
            {"pileup", required_argument, 0, 'p'},
            {"min_fraction", required_argument, 0, 'f'},
            {"max_het_bias", required_argument, 0, 'b'},
            {"min_count", required_argument, 0, 'n'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int optionIndex = 0;

        char option = getopt_long(argc, argv, "r:c:s:o:d:l:p:f:b:n:h", longOptions, &optionIndex);
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
        case 's':
            // Set the sample name
            sampleName = optarg;
            break;
        case 'o':
            // Offset variants
            variantOffset = std::stoll(optarg);
            break;
        case 'd':
            // Limit max depth for pathing to primary path
            maxDepth = std::stoll(optarg);
            break;
        case 'l':
            // Set a length override
            lengthOverride = std::stoll(optarg);
            break;
        case 'p':
            // Set a pileup filename
            pileupFilename = optarg;
            break;
        case 'f':
            // Set min fraction of average coverage for a call
            minFractionForCall = std::stod(optarg);
            break;
        case 'b':
            // Set max factor between reads on one alt and reads on the other
            // alt for calling a het.
            maxHetBias = std::stod(optarg);
            break;
        case 'n':
            // How many reads need to touch an allele before we are willing to
            // call it?
            minTotalSupportForCall = std::stoll(optarg);
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
    
    // This holds read support for all the nodes we have read support provided
    // for, by the node pointer in the vg graph.
    std::map<vg::Node*, size_t> nodeReadSupport;
    // And read support for the edges
    std::map<vg::Edge*, size_t> edgeReadSupport;
    // This holds all the edges that are deletions, by the pointer to the stored
    // Edge object in the VG graph
    std::set<vg::Edge*> deletionEdges;
    
    // This holds where nodes came from (node and offset) in the original, un-
    // augmented graph. For pieces of original nodes, this is where the piece
    // started. For novel nodes, this is where the piece that thois is an
    // alternative to started.
    std::map<vg::Node*, std::pair<int64_t, size_t>> nodeSources;
    
    // We also need to track what edges and nodes are reference (i.e. already
    // known)
    std::set<vg::Node*> knownNodes;
    std::set<vg::Edge*> knownEdges;
    
    // Loop through all the lines
    std::string line;
    size_t lineNumber = 0;
    while(std::getline(glennStream, line)) {
        // For each line
        
        lineNumber++;
        
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
                throw std::runtime_error("Line " + std::to_string(lineNumber) + ": Invalid node: " + std::to_string(nodeId));
            }
            
            // Retrieve the node we're talking about 
            vg::Node* nodePointer = vg.get_node(nodeId);
            
            // What kind of call is it? Could be "U"ncalled, or "R"eference
            // (i.e. known in the original graph), which we have special
            // handling for.
            std::string callType;
            tokens >> callType;
            
            if(callType == "U") {
                // This node has no called copy number
#ifdef debug
                std::cerr << "Line " << std::to_string(lineNumber) << ": Uncalled node: " << nodeId << endl;
#endif
                // Put it down as copy number 0 so we don't try and path through
                // it. TODO: just make the pathing code treat no-CN-stored nodes
                // as unusable.
                nodeReadSupport[vg.get_node(nodeId)] = 0;

            }
            
            // Read the read support
            size_t readSupport;
            if(tokens >> readSupport) {
                // For nodes with the number there, actually process the read support
            
#ifdef debug
                std::cerr << "Line " << std::to_string(lineNumber) << ": Node " << nodeId
                    << " has read support " << readSupport << endl;
#endif
                
                // Save it
                nodeReadSupport[nodePointer] = readSupport;
                
                if(callType == "R") {
                    // Note that this is a reference node
                    knownNodes.insert(nodePointer);
                }
            }
            
            // Load the original node ID and offset for this node, if present.
            int64_t originalId;
            size_t originalOffset;
            
            if(tokens >> originalId && tokens >> originalOffset && originalId != 0) {
                nodeSources[nodePointer] = std::make_pair(originalId, originalOffset);
            }
            
        } else if(lineType == "E") {
        
            // Read the edge data
            std::string edgeDescription;
            tokens >> edgeDescription;
            
            // Split on commas. We'd just iterate the regex iterator ourselves,
            // but it seems to only split on the first comma if we do that.
            std::vector<string> parts;
            std::regex comma_re(",");
            std::copy(std::sregex_token_iterator(edgeDescription.begin(), edgeDescription.end(), comma_re, -1), std::sregex_token_iterator(), std::back_inserter(parts));
            
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
                throw std::runtime_error("Line " + std::to_string(lineNumber) + ": Edge " + edgeDescription + " not in graph.");
            }
            
            // Get the edge
            vg::Edge* edgePointer = vg.get_edge(std::make_pair(fromSide, toSide));
            
            // Parse the mode
            std::string mode;
            tokens >> mode;
            
            if(mode == "L" || mode == "R") {
                // This is a deletion edge, or an edge in the primary path that
                // may describe a nonzero-length deletion.
#ifdef debug
                std::cerr << "Line " << std::to_string(lineNumber) << ": Edge "
                    << edgeDescription << " may describe a deletion." << endl;
#endif

                // Say it's a deletion
                deletionEdges.insert(edgePointer);
                
                if(mode == "R") {
                    // The reference edges also get marked as such
                    knownEdges.insert(edgePointer);
#ifdef debug
                    std::cerr << "Line " << std::to_string(lineNumber) << ": Edge "
                        << edgeDescription << " is reference." << endl;
#endif
                }

            }
            
            // Read the read support
            size_t readSupport;
            if(tokens >> readSupport) {
                // For nodes with the number there, actually process the read support
            
#ifdef debug
                std::cerr << "Line " << std::to_string(lineNumber) << ": Edge " << edgeDescription
                    << " has read support " << readSupport << endl;
#endif
                
                // Save it
                edgeReadSupport[edgePointer] = readSupport;
            }
        
        } else {
            // This is not a real kind of line
            throw std::runtime_error("Line " + std::to_string(lineNumber) + ": Unknown line type: " + lineType);
        }
        
        
        
        
    }
    
    std::cerr << "Loaded " << lineNumber << " lines from " << glennFile << endl;
    
    // Crunch the numbers on the reference and its read support. How much read
    // support in total (node length * aligned reads) does the primary path get?
    size_t primaryPathTotalSupport = 0;
    for(auto& pointerAndSupport : nodeReadSupport) {
        if(index.byId.count(pointerAndSupport.first->id())) {
            // This is a primary path node. Add in the total read bases supporting it
            primaryPathTotalSupport += pointerAndSupport.first->sequence().size() * pointerAndSupport.second;
        }
    }
    // Calculate average support in reads per base
    double primaryPathAverageSupport = (double)primaryPathTotalSupport / index.sequence.size();
    
    std::cerr << "Primary path average coverage: " << primaryPathAverageSupport << endl;
    
    // If applicable, load the pileup.
    // This will hold pileup records by node ID.
    std::map<int64_t, vg::NodePileup> nodePileups;
    
    std::function<void(vg::Pileup&)> handlePileup = [&](vg::Pileup& p) { 
        // Handle each pileup chunk
        for(size_t i = 0; i < p.node_pileups_size(); i++) {
            // Pull out every node pileup
            auto& pileup = p.node_pileups(i);
            // Save the pileup under its node's pointer.
            nodePileups[pileup.node_id()] = pileup;
        }
    };
    if(!pileupFilename.empty()) {
        // We have to load some pileups
        std::ifstream in;
        in.open(pileupFilename.c_str());
        stream::for_each(in, handlePileup);
    }
    
    // Generate a vcf header. We can't make Variant records without a
    // VariantCallFile, because the variants need to know which of their
    // available info fields or whatever are defined in the file's header, so
    // they know what to output.
    // Handle length override if specified.
    std::stringstream headerStream;
    write_vcf_header(headerStream, sampleName, contigName,
        lengthOverride != -1 ? lengthOverride : (index.sequence.size() + variantOffset));
    
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
        
        if(nodeReadSupport.at(node) > 0) {
            // We have copy number on this node.
            
            // Find a path to the primary reference from here
            auto path = find_bubble(vg, node, index, maxDepth);
            
            if(path.empty()) {
                // We couldn't find a path back to the primary path. Discard
                // this material.
                basesLost += node->sequence().size();
                return;
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
                
                // Budge left and add that character to the alt as well
                referenceIntervalStart--;
                altStream << index.sequence[referenceIntervalStart];
            }
            
            // Variants should be reference if most of their bases are
            // reference, and novel otherwise. A single SNP base on a known
            // insert should not make it a novel insert.
            size_t knownAltBases = 0;
            
            // How many alt bases are there overall, both known and novel?
            size_t altBases = 0;
            
            // We also need a list of all the alt node IDs for anming the
            // variant.
            std::stringstream idStream;
            
            // And collections of all the involved node pointers, on both the
            // ref and alt sides
            std::set<vg::Node*> refInvolvedNodes;
            std::set<vg::Node*> altInvolvedNodes;
            
            // And we want to know how much read support in total ther alt has
            // (support * node length)
            size_t altReadSupportTotal = 0;
            
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
                
                // Record ID
                idStream << std::to_string(path[i].node->id());
                if(i != path.size() - 2) {
                    // Add a separator (-2 since the last thing is path is an
                    // anchoring reference node)
                    idStream << "_";
                }
                
                // Record involvement
                altInvolvedNodes.insert(path[i].node);
                
                if(nodeReadSupport.count(path[i].node)) {
                    // We have read support for this node. Add it in to the total support for the alt.
                    altReadSupportTotal += path[i].node->sequence().size() * nodeReadSupport.at(path[i].node);
                }
                
                if(knownNodes.count(path[i].node)) {
                    // This is a reference node.
                    knownAltBases += path[i].node->sequence().size();
                }
                // We always need to add in the length of the node to the total
                // length
                altBases += path[i].node->sequence().size();
            }


            // Find the primary path nodes that are being skipped over/bypassed
            
            // First collect all the IDs of nodes we aren't skipping because
            // they're in the alt. Don't count those.
            std::set<int64_t> altIds;
            for(auto& visit : path) {
                altIds.insert(visit.node->id());
            }

            // Holds total primary path base readings observed (read support * node length).
            double refReadSupportTotal = 0;
            // And total bases of material we looked at on the primary path and
            // not this alt
            size_t refBases = 0;
            int64_t refNodeStart = referenceIntervalStart;
            
            while(refNodeStart < referenceIntervalPastEnd) {
            
                // Find the reference node starting here or later. Remember that
                // a variant anchored at its left base to a reference position
                // may have no node starting right where it starts.
                auto found = index.byStart.lower_bound(refNodeStart);
                if(found == index.byStart.end()) {
                    // No reference nodes here! That's a bit weird. But stop the
                    // loop.
                    break;
                }
                if((*found).first >= referenceIntervalPastEnd) {
                    // The next reference node we can find is out of the space
                    // being replaced. We're done.
                    break;
                }
                
                // Pull out the reference node we located
                auto* refNode = (*found).second.node;
                
                // Record involvement
                refInvolvedNodes.insert(refNode);
            
                // Next iteration look where this node ends.
                refNodeStart = (*found).first + refNode->sequence().size();
            
                if(altIds.count(refNode->id())) {
                    // This node is also involved in the alt we did take, so
                    // skip it. TODO: work out how to deal with shared nodes.
#ifdef debug
                    std::cerr << "Node " << refNode->id() << " also used in alt" << std::endl;
#endif
                    continue;
                }
                
                // Say we saw these bases, which may or may not have been called present
                refBases += refNode->sequence().size();
#ifdef debug
                std::cerr << "Node " << refNode->id() << " has " << nodeReadSupport.at(refNode) << " copies" << std::endl;
#endif
                
                // Count the bases we see not deleted
                refReadSupportTotal += refNode->sequence().size() * nodeReadSupport.at(refNode);
            }
            
#ifdef debug
            std::cerr << idStream.str() << " ref alternative: " << refReadSupportTotal << "/" << refBases
                << " from " << referenceIntervalStart << " to " << referenceIntervalPastEnd << std::endl;
#endif

            // We divide the read support of stuff passed over by the total
            // bases of stuff passed over to get the average read support for
            // the primary path allele.
            double refReadSupportAverage = refBases == 0 ? 0 : (double)refReadSupportTotal / refBases;
            
            // And similarly for the alt
            double altReadSupportAverage = altBases == 0 ? 0 : (double)altReadSupportTotal / altBases;

            if(refBases == 0) {
                // There's no reference node; we're a pure insert. Like with
                // deletions, we should look at the edge that bypasses us to see
                // if there's any support for it.
                
                // We eant an edge from the end of the node before us to the
                // start of the node after us.
                std::pair<vg::NodeSide, vg::NodeSide> edgeWanted = std::make_pair(
                    vg::NodeSide(path.front().node->id(), true),
                    vg::NodeSide(path.back().node->id()));
                
                if(vg.has_edge(edgeWanted)) {
                    // We found it!
                    vg::Edge* bypass = vg.get_edge(edgeWanted);
                    
                    // Any reads supporting the edge bypassing the insert are
                    // really ref support reads, and should count as supporting
                    // the whole ref allele.
                    refReadSupportTotal = edgeReadSupport.count(bypass) ? edgeReadSupport.at(bypass) : 0; 
                    refReadSupportAverage = refReadSupportTotal;
                }
            }
            // Otherwise if there's no edge or no support for that edge, the ref support should stay 0.
            
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
            
            if(knownAltBases >= altBases / 2) {
                // Flag the variant as reference. Don't put in a false entry if
                // it isn't, because vcflib will spit out the flag anyway...
                variant.infoFlags["XREF"] = true;
            }
            
            // We need to deduplicate the corss-references, because multiple
            // involved nodes may cross-reference the same original node and
            // offset, and because the same original node and offset can be
            // referenced by both ref and alt paths.
            std::set<std::pair<int64_t, size_t>> refCrossreferences;
            std::set<std::pair<int64_t, size_t>> altCrossreferences;
            std::set<std::pair<int64_t, size_t>> crossreferences;
            
            for(auto* node : refInvolvedNodes) {
                // Every involved node gets its original node:offset recorded as
                // an XSEE cross-reference.
                
                if(nodeSources.count(node)) {
                    // It has a source. Find it
                    auto& source = nodeSources.at(node);
                    // Then add it to be referenced.
                    refCrossreferences.insert(source);
                    crossreferences.insert(source);
                }
            }
            for(auto* node : altInvolvedNodes) {
                // Every involved node gets its original node:offset recorded as
                // an XSEE cross-reference.
                
                if(nodeSources.count(node)) {
                    // It has a source. Find it
                    auto& source = nodeSources.at(node);
                    // Then add it to be referenced.
                    altCrossreferences.insert(source);
                    crossreferences.insert(source);
                }
            }
            for(auto& crossreference : crossreferences) {
                variant.info["XSEE"].push_back(std::to_string(crossreference.first) + ":" +
                    std::to_string(crossreference.second));
            }
            
            // Initialize the ref allele
            create_ref_allele(variant, refAllele);
            
            // Add the alt allele
            int altNumber = add_alt_allele(variant, altAllele);
            
            // Say we're going to spit out the genotype for this sample.        
            variant.format.push_back("GT");
            auto& genotype = variant.samples[sampleName]["GT"];
            
            
            // We're going to make some really bad calls at low depth. We can
            // pull them out with a depth filter, but for now just elide them.
            if(refReadSupportAverage + altReadSupportAverage >= primaryPathAverageSupport * minFractionForCall) {
                if(refReadSupportAverage > maxHetBias * altReadSupportAverage &&
                    refReadSupportTotal >= minTotalSupportForCall) {
                    // Biased enough towards ref, and ref has enough total reads.
                    // Say it's hom ref
                    genotype.push_back("0/0");
                } else if(altReadSupportAverage > maxHetBias * refReadSupportAverage
                    && altReadSupportTotal >= minTotalSupportForCall) {
                    // Say it's hom alt
                    genotype.push_back(std::to_string(altNumber) + "/" + std::to_string(altNumber));
                } else if(refReadSupportTotal >= minTotalSupportForCall &&
                    altReadSupportTotal >= minTotalSupportForCall) {
                    // Say it's het
                    genotype.push_back("0/" + std::to_string(altNumber));
                } else {
                    // We can't really call this as anything.
                    genotype.push_back("./.");
                }
            } else {
                // Depth too low. Say we have no idea.
                // TODO: elide variant?
                genotype.push_back("./.");
            }
            // TODO: use legit thresholds here.
            
            // Add depth for the variant and the samples
            std::string depthString = std::to_string((int64_t)round(refReadSupportAverage + altReadSupportAverage));
            variant.format.push_back("DP");
            variant.samples[sampleName]["DP"].push_back(depthString);
            variant.info["DP"].push_back(depthString); // We only have one sample, so variant depth = sample depth
            
            // Also allelic depths
            variant.format.push_back("AD");
            variant.samples[sampleName]["AD"].push_back(std::to_string((int64_t)round(refReadSupportAverage)));
            variant.samples[sampleName]["AD"].push_back(std::to_string((int64_t)round(altReadSupportAverage)));
            
            
#ifdef debug
            std::cerr << "Found variant " << refAllele << " -> " << altAllele
                << " caused by nodes " <<  variant.id
                << " at 1-based reference position " << variant.position
                << std::endl;
#endif

            if(can_write_alleles(variant)) {
                // Output the created VCF variant.
                std::cout << variant << std::endl;
                
                // Output the pileup line, which will be nonempty if we have pileups
                std::cout << get_pileup_line(nodePileups, refCrossreferences, altCrossreferences);
            
            } else {
                std::cerr << "Variant is too large" << std::endl;
                // TODO: account for the 1 base we added extra if it was a pure
                // insert.
                basesLost += altAllele.size();
            }
            
            
        }
        
        
    });
    
    for(vg::Edge* deletion : deletionEdges) {
        // Make deletion variants for each deletion edge
        
        // Make a string naming the edge
        std::string edgeName = std::to_string(deletion->from()) +
            (deletion->from_start() ? "L" : "R") + "->" +
            std::to_string(deletion->to()) + (deletion->to_end() ? "R" : "L");
        
        if(!index.byId.count(deletion->from()) || !index.byId.count(deletion->to())) {
            // This deletion edge does not cover a reference interval.
            // TODO: take into account its presence when pushing copy number.
#ifdef debug
            std::cerr << "Deletion edge " << edgeName << " does not cover a reference interval. Skipping!" << endl;
#endif
            continue;
        }
        
        // Where are we from and to in the reference (leftmost position and
        // relative orientation)
        auto& fromPlacement = index.byId.at(deletion->from());
        auto& toPlacement = index.byId.at(deletion->to());
        
#ifdef debug
        std::cerr << "Node " << deletion->from() << " is at ref position " 
            << fromPlacement.first << " orientation " << fromPlacement.second << std::endl;
        std::cerr << "Node " << deletion->to() << " is at ref position "
            << toPlacement.first << " orientation " << toPlacement.second << std::endl;
#endif
        
        // Are we attached to the reference-relative left or right of our from
        // base?
        bool fromFirst = fromPlacement.second != deletion->from_start();
        
        // And our to base?
        bool toLast = toPlacement.second != deletion->to_end();
        
        // What base should the from end really be on? This is the non-deleted
        // base outside the deletion on the from end.
        int64_t fromBase = fromPlacement.first + (fromFirst ? 0 : vg.get_node(deletion->from())->sequence().size() - 1);
        
        // And the to end?
        int64_t toBase = toPlacement.first + (toLast ? vg.get_node(deletion->to())->sequence().size() - 1 : 0);

        if(toBase <= fromBase) {
            // Our edge ought to be running backward.
            if(!(fromFirst && toLast)) {
                // We're not a proper deletion edge in the backwards spelling
                // Discard the edge
                std::cerr << "Improper deletion edge " << edgeName << std::endl;
                basesLost += toBase - fromBase;
                continue;
            } else {
                // Just invert the from and to bases.
                std::swap(fromBase, toBase);
#ifdef debug
                std::cerr << "Inverted deletion edge " << edgeName << std::endl;
#endif
            }
        } else if(fromFirst || toLast) {
            // We aren't a proper deletion edge in the forward spelling either.
            std::cerr << "Improper deletion edge " << edgeName << std::endl;
            basesLost += fromBase - toBase;
            continue;
        }
        
        
        if(toBase <= fromBase + 1) {
            // No bases were actually deleted. Maybe this is just a normal reference edge.
            continue;
        }
        
#ifdef debug
        std::cerr << "Deletion " << edgeName << " bookended by " << fromBase << " and " << toBase << std::endl;
#endif
        
        // What original node:offset places do we care about?
        std::set<std::pair<int64_t, size_t>> crossreferences;
        
        // Guess the copy number of the deletion.
        // Holds total base copies (node length * read support) observed as not deleted.
        double refReadSupportTotal = 0;
        int64_t deletedNodeStart = fromBase + 1;
        while(deletedNodeStart != toBase) {
#ifdef debug
            std::cerr << "Next deleted node starts at " << deletedNodeStart << std::endl;
#endif
        
            // Find the deleted node starting here in the reference
            auto* deletedNode = index.byStart.at(deletedNodeStart).node;
            // We know the next reference node should start just after this one.
            // Even if it previously existed in the reference.
            deletedNodeStart += deletedNode->sequence().size();
            
            // Count the read observations we see not deleted
            refReadSupportTotal += deletedNode->sequence().size() * nodeReadSupport.at(deletedNode);
            
            // Add the beginning of this node as a see also. The deletion is
            // stored at the beginning of the first node, and the other
            // locations will help us get an idea of how much support the
            // deleted stuff has.
            crossreferences.insert(nodeSources.at(deletedNode));
        }
        
        // We divide the total read support by the total bases deleted to get
        // the average read support for the reference.
        double refReadSupportAverage = (double)refReadSupportTotal / (toBase - fromBase - 1);
        
        // Get the support for the edge itself?
        size_t altReadSupportTotal = edgeReadSupport.count(deletion) ? edgeReadSupport[deletion] : 0;
        
        // No sense averaging the deletion edge read support because there are no bases.
        
        
        // What copy number do we call for the deletion?
        int64_t copyNumberCall = 2;
        
        // We're going to make some really bad calls at low depth. We can
        // pull them out with a depth filter, but for now just elide them.
        if(refReadSupportAverage + altReadSupportTotal >= primaryPathAverageSupport * minFractionForCall) {
            if(refReadSupportAverage > maxHetBias * altReadSupportTotal &&
                refReadSupportTotal >= minTotalSupportForCall) {
                // Say it's hom ref
                copyNumberCall = 0;
            } else if(altReadSupportTotal > maxHetBias * refReadSupportAverage &&
                altReadSupportTotal >= minTotalSupportForCall) {
                // Say it's hom alt
                copyNumberCall = 2;
            } else if(refReadSupportTotal >= minTotalSupportForCall &&
                altReadSupportTotal >= minTotalSupportForCall) {
                // Say it's het
                copyNumberCall = 1;
            } else {
                // We're not biased enough towards either homozygote, but we
                // don't have enough support for each allele to call het.
                // TODO: we just don't call
                copyNumberCall = 0;
            }
        } else {
            // Depth too low. Don't call.
            copyNumberCall = 0;
        }
        // TODO: use legit thresholds here.
        
        
        if(copyNumberCall == 0) {
            // Actually don't call a deletion
            continue;
        }
        
        // Now we know fromBase is the last non-deleted base and toBase is the
        // first non-deleted base. We'll make an alt replacing the first non-
        // deleted base plus the deletion with just the first non-deleted base.
        // Rename everything to the same names we were using before.
        size_t referenceIntervalStart = fromBase;
        size_t referenceIntervalPastEnd = toBase;
        
        // Make the variant and emit it.
        std::string refAllele = index.sequence.substr(
            referenceIntervalStart, referenceIntervalPastEnd - referenceIntervalStart);
        std::string altAllele = index.sequence.substr(referenceIntervalStart, 1);
        
        // Make a Variant
        vcflib::Variant variant;
        variant.sequenceName = contigName;
        variant.setVariantCallFile(vcf);
        variant.quality = 0;
        variant.position = referenceIntervalStart + 1 + variantOffset;
        variant.id = edgeName;
        
        if(knownEdges.count(deletion)) {
            // Mark it as reference if it is a reference edge. Apparently vcflib
            // does not check the flag value when serializing, so don't put in a
            // flase entry if it's not a reference edge.œ
            variant.infoFlags["XREF"] = true;
#ifdef debug
            std::cerr << edgeName << " is a reference deletion" << std::endl;
#endif
        }
        
        for(auto& crossreference : crossreferences) {
            // Add in all the deduplicated cross-references
            variant.info["XSEE"].push_back(std::to_string(crossreference.first) + ":" +
                std::to_string(crossreference.second));
        }
        
        // Initialize the ref allele
        create_ref_allele(variant, refAllele);
        
        // Add the alt allele
        int altNumber = add_alt_allele(variant, altAllele);
        
        // Say we're going to spit out the genotype for this sample.        
        variant.format.push_back("GT");
        auto& genotype = variant.samples[sampleName]["GT"];
        
        if(copyNumberCall == 1) {
            // We're allele alt and ref heterozygous.
            genotype.push_back("0/" + std::to_string(altNumber));
        } else if(copyNumberCall == 2) {
            // We're alt homozygous, other overlapping variants notwithstanding.
            genotype.push_back(std::to_string(altNumber) + "/" + std::to_string(altNumber));
        } else {
            // We're something weird
            throw std::runtime_error("Invalid copy number for deletion: " + std::to_string(copyNumberCall));
        }
        
        // Add depth for the variant and the samples
        std::string depthString = std::to_string((int64_t)round(refReadSupportAverage + altReadSupportTotal));
        variant.format.push_back("DP");
        variant.samples[sampleName]["DP"].push_back(depthString);
        variant.info["DP"].push_back(depthString); // We only have one sample, so variant depth = sample depth
        
        // Also allelic depths
        variant.format.push_back("AD");
        variant.samples[sampleName]["AD"].push_back(std::to_string((int64_t)round(refReadSupportAverage)));
        variant.samples[sampleName]["AD"].push_back(std::to_string(altReadSupportTotal));
        
#ifdef debug
        std::cerr << "Found variant " << refAllele << " -> " << altAllele
            << " caused by edge " <<  variant.id
            << " at 1-based reference position " << variant.position
            << std::endl;
#endif

        if(can_write_alleles(variant)) {
            // Output the created VCF variant.
            std::cout << variant << std::endl;
            
            // Output the pileup line, which will be nonempty if we have pileups
            // We only have ref crossreferences here. TODO: make the ref and alt
            // labels make sense for deletions/re-design the way labeling works.
            std::cout << get_pileup_line(nodePileups, crossreferences, std::set<std::pair<int64_t, size_t>>());
            
        } else {
            std::cerr << "Variant is too large" << std::endl;
            // TODO: Drop the anchoring base that doesn't really belong to the
            // deletion, when we can be consistent with inserts.
            basesLost += altAllele.size();
        }
        
    }
    
    // Announce how much we can't show.
    std::cerr << "Had to drop " << basesLost << " bp of unrepresentable variation." << std::endl;
    
    return 0;
}



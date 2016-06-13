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

// We represent support as a pair, but we define math for it.
// We use doubles because we may need fractional math.
typedef std::pair<double, double> Support;

/**
 * Add two Support values together, accounting for strand.
 */
Support operator+(const Support& one, const Support& other) {
    return std::make_pair(one.first + other.first, one.second + other.second);
}

/**
 * Add in a Support to another.
 */
Support& operator+=(Support& one, const Support& other) {
    one.first += other.first;
    one.second += other.second;
    return one;
}


/**
 * Scale a support by a factor.
 */
template<typename Scalar>
Support operator*(const Support& support, const Scalar& scale) {
    return std::make_pair(support.first * scale, support.second * scale);
}

/**
 * Scale a support by a factor, the other way
 */
template<typename Scalar>
Support operator*(const Scalar& scale, const Support& support) {
    return std::make_pair(support.first * scale, support.second * scale);
}

/**
 * Divide a support by a factor.
 */
template<typename Scalar>
Support operator/(const Support& support, const Scalar& scale) {
    return std::make_pair(support.first / scale, support.second / scale);
}

/**
 * Allow printing a support.
 */
std::ostream& operator<<(std::ostream& stream, const Support& support) {
    return stream << support.first << "," << support.second;
}

/**
 * Get the total read support in a support.
 */
double total(const Support& support) {
    return support.first + support.second;
}

/**
 * Get the strand bias of a support.
 */
double strand_bias(const Support& support) {
    return std::max(support.first, support.second) / (support.first + support.second);
}

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
    stream << "##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Forward and reverse support for ref and alt alleles.\">" << std::endl;
    // We need this field to stratify on for VCF comparison. The info is in SB but vcfeval can't pull it out
    stream << "##FORMAT=<ID=XAAD,Number=1,Type=Integer,Description=\"Alt allele read count.\">" << std::endl;
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
 * Given a pair of nodes in order (constituting a superbubble), return unique
 * traversals through the graph from the first node to the last node.
 *
 * The start node must be upstream of the end node in the given orientations,
 * because we will try to leave the start traversal on its right, to enter the
 * end traversal on its left.
 */
std::vector<std::list<vg::NodeTraversal>> find_alleles(vg::VG& graph,
    const vg::NodeTraversal& start, const vg::NodeTraversal& end,
    const std::map<vg::Node*, Support>& nodeReadSupport, int64_t maxDepth = 10) {
    
    // Holds paths we want to return.
    // TODO: does this need to be a set, really?
    std::vector<std::list<vg::NodeTraversal>> toReturn;
    
    // Do a BFS
    
    // This holds the paths to extend from
    std::list<std::list<vg::NodeTraversal>> toExtend;
    
    // Start at the start node
    toExtend.emplace_back(std::list<vg::NodeTraversal> {start});
    
    while(!toExtend.empty()) {
    
        // Dequeue a path to extend.
        // Make sure to move out of the list to avoid a useless copy.
        std::list<vg::NodeTraversal> path(std::move(toExtend.front()));
        toExtend.pop_front();
        
        if(path.back() == end) {
            // We have finally reached the end!
            
            // Say we got to the right place
            toReturn.emplace_back(std::move(path));
            
            // Don't bother looking for extensions, we already got there.
        } else if(path.size() < maxDepth) {
            // We haven't hit the end yet, but we also haven't hit the max
            // depth. Extend with all the possible extensions.
            
            // Look right
            vector<vg::NodeTraversal> nextNodes = graph.nodes_next(path.back());
            
            for(auto nextNode : nextNodes) {
                // For each node we can get to
                
                if(!nodeReadSupport.empty() && (!nodeReadSupport.count(nextNode.node) ||
                    total(nodeReadSupport.at(nextNode.node)) == 0)) {
                    
                    // We have no support at all for visiting this node (but we
                    // do have some node read support data). Don't visit it.
                    continue;
                }
                
                // Make a new path extended right with the node
                std::list<vg::NodeTraversal> extended(path);
                extended.push_back(nextNode);
                toExtend.emplace_back(std::move(extended));
            }
        }
    }
    
    // Return the set of alleles for the bubble, with leading and trailing start
    // and end traversals.
    return toReturn;
    
    
}

/**
 * Do a breadth-first search left from the given node traversal, and return
 * lengths and paths starting at the given node and ending on the indexed
 * reference path. Refuses to visit nodes with no support.
 */
std::set<std::pair<size_t, std::list<vg::NodeTraversal>>> bfs_left(vg::VG& graph,
    vg::NodeTraversal node, const ReferenceIndex& index,
    const std::map<vg::Node*, Support>& nodeReadSupport, int64_t maxDepth = 10,
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
                
                if(!nodeReadSupport.empty() && (!nodeReadSupport.count(prevNode.node) ||
                    total(nodeReadSupport.at(prevNode.node)) == 0)) {
                    
                    // We have no support at all for visiting this node (but we
                    // do have some node read support data)
                    continue;
                }
                
                if(stopIfVisited && alreadyQueued.count(prevNode)) {
                    // We already have a way to get here.
                    continue;
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
    vg::NodeTraversal node, const ReferenceIndex& index,
    const std::map<vg::Node*, Support>& nodeReadSupport, int64_t maxDepth = 10,
    bool stopIfVisited = false) {

    // Look left from the backward version of the node.
    auto toConvert = bfs_left(graph, flip(node), index, nodeReadSupport, maxDepth, stopIfVisited);
    
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
    const std::map<vg::Node*, Support>& nodeReadSupport, int64_t maxDepth = 10) {

    // Find paths on both sides, with nodes on the primary path at the outsides
    // and this node in the middle. Returns path lengths and paths in pairs in a
    // set.
    auto leftPaths = bfs_left(graph, vg::NodeTraversal(node), index, nodeReadSupport, maxDepth);
    auto rightPaths = bfs_right(graph, vg::NodeTraversal(node), index, nodeReadSupport, maxDepth);
    
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
    
    // Fix up doubly reversing edges, which upset the sduperbubble-finding
    // algorithm
    vg.flip_doubly_reversed_edges();
    
    if(refPathName.empty()) {
        std::cerr << "Graph has " << vg.paths.size() << " paths to choose from."
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
    
    // Copy the graph and dagify the copy. We need to hold this translation
    // table from new node ID to old node and relative orientation.
    map<vg::id_t, pair<vg::id_t, bool>> dagTranslation;
    vg::VG* dag = new vg::VG(std::move(vg.dagify(1, dagTranslation)));
    
    std::cerr << "Computing superbubbles..." << std::endl;
    
    // Find the superbubbles in the DAG
    std::vector<std::pair<vg::id_t, vg::id_t>> superbubbles = dag->get_superbubbles();
    
    for(auto& superbubble : superbubbles) {
        // Translate them back to the original node ID space.
        // Relative orientations don't really matter here.
        superbubble.first = dagTranslation[superbubble.first].first;
        superbubble.second = dagTranslation[superbubble.second].first;
    }
    
    // Discard the copy.
    delete dag;
    dagTranslation.clear();
    
    std::cerr << "Graph has " << superbubbles.size() << " superbubbles." << std::endl;
    
    // Follow the reference path and extract indexes we need: index by node ID,
    // index by node start, and the reconstructed path sequence.
    ReferenceIndex index = trace_reference_path(vg, refPathName);
    
    // Open up the Glenn-file
    std::ifstream glennStream(glennFile);
    
    // Parse it into an internal format, where we track status and copy number
    // for nodes and edges.
    
    // This holds read support, on each strand, for all the nodes we have read
    // support provided for, by the node pointer in the vg graph.
    std::map<vg::Node*, Support> nodeReadSupport;
    // And read support for the edges
    std::map<vg::Edge*, Support> edgeReadSupport;
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
                nodeReadSupport[vg.get_node(nodeId)] = std::make_pair(0.0, 0.0);

            }
            
            // Read the read support
            Support readSupport;
            if((tokens >> readSupport.first) && (tokens >> readSupport.second)) {
                // For nodes with the number there, actually process the read support
            
#ifdef debug
                std::cerr << "Line " << std::to_string(lineNumber) << ": Node " << nodeId
                    << " has read support " << readSupport.first << "," << readSupport.second << endl;
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
            Support readSupport;
            if((tokens >> readSupport.first) && (tokens >> readSupport.second)) {
                // For nodes with the number there, actually process the read support
            
#ifdef debug
                std::cerr << "Line " << std::to_string(lineNumber) << ": Edge " << edgeDescription
                    << " has read support " << readSupport.first << "," << readSupport.second << endl;
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
    Support primaryPathTotalSupport = std::make_pair(0.0, 0.0);
    for(auto& pointerAndSupport : nodeReadSupport) {
        if(index.byId.count(pointerAndSupport.first->id())) {
            // This is a primary path node. Add in the total read bases supporting it
            primaryPathTotalSupport += pointerAndSupport.first->sequence().size() * pointerAndSupport.second;
        }
    }
    // Calculate average support in reads per base
    auto primaryPathAverageSupport = primaryPathTotalSupport / index.sequence.size();
    
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
    
    for(auto startAndEnd : superbubbles) {
        // Get where the two nodes at the ends of the bubble are in the
        // reference, and in what orientations they occur
        
#ifdef debug
        std::cerr << "Evaluating superbubble " << startAndEnd.first << " to " << 
                startAndEnd.second << std::endl;
#endif
        
        if(!index.byId.count(startAndEnd.first) || !index.byId.count(startAndEnd.second)) {
            // We aren't anchored to the primary path, so there's not much we
            // can do.
            std::cerr << "Superbubble " << startAndEnd.first << " to " << 
                startAndEnd.second << " is not anchored to the reference on both ends!" << std::endl;
            // Try the next superbubble
            continue;
        }
        
        auto startPositionAndOrientation = index.byId.at(startAndEnd.first);     
        auto endPositionAndOrientation = index.byId.at(startAndEnd.second);
        
        if(startPositionAndOrientation.first >= endPositionAndOrientation.first) {
            // In the wrong order! Swap them!
            std::cerr << "Warning: swapping nodes " << startAndEnd.first << 
                " and " << startAndEnd.second << " in superbubble!" << std::endl;
            std::swap(startAndEnd.first, startAndEnd.second);
            std::swap(startPositionAndOrientation, endPositionAndOrientation);
        }
        
        // Make traversals in the correct orientations for the start and end nodes.
        vg::NodeTraversal start(vg.get_node(startAndEnd.first), startPositionAndOrientation.second);
        vg::NodeTraversal end(vg.get_node(startAndEnd.second), endPositionAndOrientation.second);
        
        // Look for all the alleles
        std::vector<std::list<vg::NodeTraversal>> alleles = find_alleles(vg, start, end, nodeReadSupport, maxDepth);
        
        if(alleles.size() < 2) {
            // We don't want to try and emit variants with less than 2 alleles.
            std::cerr << "Warning! Fewer than 2 alleles at superbubble " << startAndEnd.first << 
                " to " << startAndEnd.second << std::endl;
            // TODO: we get lots of non-variable superbubbles by default it seems.
            // Account for dropped bases
            basesLost += endPositionAndOrientation.first - startPositionAndOrientation.first;
            continue;
        }
        
        // Determine the sequence of the ref allele, by cutting out the right bit of the reference
        // The position we have stored for this start node is the first
        // position along the reference at which it occurs. Our bubble
        // goes forward in the reference, so we must come out of the
        // opposite end of the node from the one we have stored.
        auto referenceIntervalStart = index.byId.at(start.node->id()).first +
            start.node->sequence().size();
        
        // The position we have stored for the end node is the first
        // position it occurs in the reference, and we know we go into
        // it in a reference-concordant direction, so we must have our
        // past-the-end position right there.
        auto referenceIntervalPastEnd = index.byId.at(end.node->id()).first;
        
        // Find the reference allele, if it's not already in the alleles we
        // found. We begin with the start node.
        std::list<vg::NodeTraversal> refAllele{start};
        
        // This tracks where the next node in the ref allele will start
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
            
            // Next iteration look where this node ends.
            refNodeStart = (*found).first + refNode->sequence().size();
            
            // Stick the traversal in the allele.
            refAllele.push_back((*found).second);
        }
        // The ref allele will always end with the end node, which we won't find
        // in the loop.
        refAllele.push_back(end);
        
        for(size_t i = 0; i < alleles.size(); i++) {
            if(alleles[i] == refAllele) {
                // Promote the ref allele to position 0
                std::swap(alleles[0], alleles[i]);
                break;
            }
        }
        
        if(alleles[0] != refAllele) {
            // We didn't manage to find the ref allele in our BFS
            std::cerr << "Warning! No ref allele found!" << std::endl;
            basesLost += referenceIntervalPastEnd - referenceIntervalStart;
            
            // Skip the rest of the superbubble
            continue;
        }
        
        // Otherwise we now have the ref allele, which we know falls at index 0,
        // and some alts.
        
        // For each allele, compute sequence and support and length.
        std::vector<string> alleleSequence(alleles.size());
        std::vector<Support> alleleSupport(alleles.size());
        std::vector<size_t> alleleLength(alleles.size());
        // Also average support
        std::vector<Support> alleleAverageSupport(alleles.size());
        // We also want to know if any alleles are empty
        bool haveEmptyAllele = false;
        
        // We also record all the IDs involved
        std::stringstream idStream;
        
        // We also need to find the most and second most supported alleles.
        int64_t mostSupportedAllele = -1;
        Support mostSupportedAlleleSupport = std::make_pair(0.0, 0.0);
        
        int64_t secondMostSupportedAllele = -1;
        Support secondMostSupportedAlleleSupport = std::make_pair(0.0, 0.0);
        
        for(size_t i = 0; i < alleles.size(); i++) {
            // Copy to a vector for random access, so skipping the first and last nodes is easier.
            std::vector<vg::NodeTraversal> allele{alleles[i].begin(), alleles[i].end()};
        
            // We stream the allele's sequence together
            std::stringstream seqStream;
        
            for(int64_t j = 1; j + 1 < allele.size(); j++) {
                // For all but the first and last nodes, grab their sequences in
                // the correct orientation.
                
                std::string addedSequence = allele[j].node->sequence();
            
                if(allele[j].backward) {
                    // If the node is traversed backward, we need to flip its sequence.
                    addedSequence = vg::reverse_complement(addedSequence);
                }
                
                // Stick the sequence
                seqStream << addedSequence;
                
                // Record ID
                idStream << std::to_string(allele[j].node->id());
                // And separator
                idStream << "_";
                
                if(nodeReadSupport.count(allele[j].node)) {
                    // We have read support for this node. Add it in to the total support for the alt.
                    alleleSupport[i] += allele[j].node->sequence().size() * nodeReadSupport.at(allele[j].node);
                }
                
                // We always need to add in the length of the node to the total
                // length
                alleleLength[i] += allele[j].node->sequence().size();
            }
            // Finish the allele sequence
            alleleSequence[i] = seqStream.str();
            
            assert(alleleSequence[i].size() == alleleLength[i]);
            if(alleleLength[i] == 0) {
                // Remember that at least one allele is empty (and so we'll have to add a reference-anchoring base)
                haveEmptyAllele = true;
            }
            
#ifdef debug
            std::cerr << "Allele " << i << ": " << alleleSequence[i] << std::endl;
#endif
            
            if(allele.size() == 2) {
                // This is a pure deletion/non-insertion allele. Its support is just the support for its edge.
                
                // We want an edge from the end of the left anchoring node to
                // the start of the right anchoring node.
                std::pair<vg::NodeSide, vg::NodeSide> edgeWanted = std::make_pair(
                    vg::NodeSide(allele.front().node->id(), !allele.front().backward),
                    vg::NodeSide(allele.back().node->id(), allele.back().backward));
                
                if(vg.has_edge(edgeWanted)) {
                    // We found it!
                    vg::Edge* bypass = vg.get_edge(edgeWanted);
                    
                    // Any reads supporting the edge bypassing the insert are
                    // really ref support reads, and should count as supporting
                    // the whole ref allele.
                    alleleSupport[i] = edgeReadSupport.count(bypass) ? edgeReadSupport.at(bypass) : std::make_pair(0.0, 0.0); 
                    // Hack the allele length to 1 so we average to this value
                    alleleLength[i] = 1;
                }
            }
            
            // Compute the average support (total bases of support over length)
            alleleAverageSupport[i] = alleleSupport[i] / alleleLength[i];
            
            if(mostSupportedAllele == -1 || total(alleleAverageSupport[i]) > total(mostSupportedAlleleSupport)) {
                // New most supported allele!
                // Demote the old one
                secondMostSupportedAllele = mostSupportedAllele;
                secondMostSupportedAlleleSupport = mostSupportedAlleleSupport;
                
                // And replace it
                mostSupportedAllele = i;
                mostSupportedAlleleSupport = alleleAverageSupport[i];
            } else if(secondMostSupportedAllele == -1 || 
                total(alleleAverageSupport[i]) > total(secondMostSupportedAlleleSupport)) {
                // Replace the old secomd most supported allele instead
                secondMostSupportedAllele = i;
                secondMostSupportedAlleleSupport = alleleAverageSupport[i];
            }
        }
        
        // We need to have both of these alleles found.
        assert(mostSupportedAllele != -1);
        assert(secondMostSupportedAllele != -1);
        assert(mostSupportedAllele != secondMostSupportedAllele);
        
        if(haveEmptyAllele) {
            // Apply anchoring bases if ref is empty. We know we won't be
            // replacing the very first base of the reference, so this is easy.
            assert(referenceIntervalStart > 0);
            
            // Adjust the reference interval start, from which the variant
            // position is derived
            referenceIntervalStart--;
            // Find the base before this variation
            char anchoringBase = index.sequence[referenceIntervalStart];
            for(size_t i = 0; i < alleles.size(); i++) {
                // Extend each allele on the left with the anchoring base.
                alleleSequence[i].insert(alleleSequence[i].begin(), anchoringBase);
                alleleLength[i]++;
            }
        }
        
        // Now we only consider the most and second most supported alleles when
        // deciding the genotype.
        
        // Make a Variant
        vcflib::Variant variant;
        variant.sequenceName = contigName;
        variant.setVariantCallFile(vcf);
        variant.quality = 0;
        variant.position = referenceIntervalStart + 1 + variantOffset;
        variant.id = idStream.str();
        // Lop off the trailing underscore that we placed after the last node
        // ID.
        assert(variant.id.size() > 0);
        variant.id.pop_back();
        
        // Initialize the ref allele sequence
        create_ref_allele(variant, alleleSequence[0]);
        
        // Add the alt alleles.
        // Record alt number for the most supported allele.
        int mostSupportedAlt = 0;
        if(mostSupportedAllele > 0) {
            // We need an alt for the most supported allele
            mostSupportedAlt = add_alt_allele(variant, alleleSequence[mostSupportedAllele]);
        }
        // Otherwise it's ref
        
        // And for the second most supported allele.
        int secondMostSupportedAlt = 0;
        if(secondMostSupportedAllele > 0) {
            // We need an alt for the second most supported allele
            secondMostSupportedAlt = add_alt_allele(variant, alleleSequence[secondMostSupportedAllele]);
        }
        // Otherwise it's ref
        
        // Say we're going to spit out the genotype for this sample.        
        variant.format.push_back("GT");
        auto& genotype = variant.samples[sampleName]["GT"];
        
        // Make some references for convenience
        auto& majorSupport = alleleAverageSupport[mostSupportedAllele];
        auto& minorSupport = alleleAverageSupport[secondMostSupportedAllele];
        
        // We're going to make some really bad calls at low depth. We can
        // pull them out with a depth filter, but for now just elide them.
        if(total(majorSupport + minorSupport) >= total(primaryPathAverageSupport) * minFractionForCall) {
            if(total(majorSupport) > maxHetBias * total(minorSupport) &&
                total(majorSupport) >= minTotalSupportForCall) {
                // Biased enough towards major allele, and it has enough total
                // reads. Say it's hom major
                genotype.push_back(std::to_string(mostSupportedAlt) + "/" + std::to_string(mostSupportedAlt));
            } else if(total(majorSupport) >= minTotalSupportForCall &&
                total(minorSupport) >= minTotalSupportForCall) {
                // Say it's het
                genotype.push_back(std::to_string(mostSupportedAlt) + "/" + std::to_string(secondMostSupportedAlt));
            } else if(total(majorSupport) >= minTotalSupportForCall) {
                // We have enough for one of the major allele, and not any of
                // the minor allele, and we can't call it het. Say we have one major allele copy and one... something.
                genotype.push_back(std::to_string(mostSupportedAlt) + "/.");
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
        
        // Add depth for the variant and the samples. It ends up being the total of all the supports.
        double totalDepth = 0.0;
        for(auto& support : alleleAverageSupport) {
            totalDepth += total(support);
        }
        std::string depthString = std::to_string((int64_t)round(totalDepth));
        variant.format.push_back("DP");
        variant.samples[sampleName]["DP"].push_back(depthString);
        variant.info["DP"].push_back(depthString); // We only have one sample, so variant depth = sample depth
        
        // Also allelic depths for all alleles
        variant.format.push_back("AD");
        for(auto& support : alleleAverageSupport) {
            variant.samples[sampleName]["AD"].push_back(std::to_string((int64_t)round(total(support))));
        }
                
        // Also strand biases
        // TODO: this is sort of silly for multiallelic sites.
        // Pick an alt allele, which is the most supported allele unles that's the reference.
        int altAllele = mostSupportedAllele == 0 ? secondMostSupportedAllele : mostSupportedAllele;
        variant.format.push_back("SB");
        variant.samples[sampleName]["SB"].push_back(std::to_string((int64_t)round(alleleAverageSupport[0].first)));
        variant.samples[sampleName]["SB"].push_back(std::to_string((int64_t)round(alleleAverageSupport[0].second)));
        variant.samples[sampleName]["SB"].push_back(std::to_string((int64_t)round(alleleAverageSupport[altAllele].first)));
        variant.samples[sampleName]["SB"].push_back(std::to_string((int64_t)round(alleleAverageSupport[altAllele].second)));
        
        // And total alt allele depth
        // TODO: total minor allele depth?
        Support altAlleleDepth = std::make_pair(0.0, 0.0);
        if(mostSupportedAllele != 0) {
            // Most supported allele involved isn't ref, so add it in.
            altAlleleDepth += majorSupport;
        }
        if(secondMostSupportedAllele != 0) {
            // Second most supported allele involved isn't ref, so add it in.
            altAlleleDepth += minorSupport;
        }
        // Then report total support for used alt alleles.
        variant.format.push_back("XAAD");
        variant.samples[sampleName]["XAAD"].push_back(std::to_string((int64_t)round(total(altAlleleDepth))));
        
        
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
            basesLost += alleleLength[mostSupportedAllele] + alleleLength[secondMostSupportedAllele];
        }
    }
    
    // This handles all substitutions, inserts, and deletes. No need for special
    // handling later.
    
    // Announce how much we can't show.
    std::cerr << "Had to drop " << basesLost << " bp of unrepresentable variation." << std::endl;
    
    return 0;
}



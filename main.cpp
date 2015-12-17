// main.cpp: Main file for core graph merger

#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <set>
#include <getopt.h>

#include "ekg/vg/vg.hpp"
#include "ekg/vg/index.hpp"
#include "ekg/vg/vcflib/src/Variant.h"

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
void write_vcf_header(std::ostream& stream, std::string sample_name) {
   stream << "##fileformat=VCFv4.2" << std::endl;
   stream << "##FORMAT=<ID=GT,Number=1,Type=Integer,Description=\"Genotype\">" << std::endl;
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
    // Make it 0 in the alleles-by-index list
    variant.alleles.push_back(allele);
    // Build the reciprocal index-by-allele mapping
    variant.updateAlleleIndexes();
}

/**
 * Add a new alt allele to a vcflib Variant, since apaprently there's no method
 * for that already.
 */
void add_alt_allele(vcflib::Variant& variant, const std::string& allele) {
    // Add it as an alt
    variant.alt.push_back(allele);
    // Make it next in the alleles-by-index list
    variant.alleles.push_back(allele);
    // Build the reciprocal index-by-allele mapping
    variant.updateAlleleIndexes();
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
        << "    -h, --help          print this help message" << std::endl;
}

int main(int argc, char** argv) {
    
    if(argc == 1) {
        // Print the help
        help_main(argv);
        return 1;
    }
    
    optind = 1; // Start at first real argument
    bool optionsRemaining = true;
    while(optionsRemaining) {
        static struct option longOptions[] = {
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int optionIndex = 0;

        switch(getopt_long(argc, argv, "h", longOptions, &optionIndex)) {
        // Option value is in global optarg
        case 'h': // When the user asks for help
        case '?': // When we get options we can't parse
            help_main(argv);
            exit(1);
            break;
        default:
            // TODO: keep track of the option
            std::cerr << "Illegal option" << std::endl;
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
    
    // Open up the Glenn-file
    std::ifstream glennStream(glennFile);
    
    // We can't make Variant records without a VariantCallFile, because the
    // variants need to know which of their available info fields or whatever
    // are defined in the file's header, so they know what to output.
    
    // Generate a header
    std::stringstream headerStream;
    // TODO: get sample name from file or a command line option.
    write_vcf_header(headerStream, "SAMPLE");
    
    // Load the headers into a new VCF
    vcflib::VariantCallFile vcf;
    std::string headerString = headerStream.str();
    assert(vcf.openForOutput(headerString));
    
    // Spit out the header
    std::cout << headerStream.str();
    
    // Loop through all the lines
    std::string line;
    while(std::getline(glennStream, line)) {
        // For each line
        
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
        
        // Read the call string
        std::string call;
        tokens >> call;
        
        // Split that out into a set of call character strings. This is how we
        // split on commas in C++. Ask for the non-matched parts of a regex
        // iterator that matches commas.
        std::set<std::string> callCharacters(std::sregex_token_iterator(
            call.begin(), call.end(), std::regex(","), -1),
            std::sregex_token_iterator());
        
        // TODO: real conversion
        
        // Spit out a nonsense VCF line
        vcflib::Variant variant;
        variant.setVariantCallFile(vcf);
        
        // Initialize the ref allele
        create_ref_allele(variant, char_to_string(graphBase));
        
        for(auto character : callCharacters) {
            // For each character (in its own string) that we called
            if(character != "-" && character != ".") {
                // It's an alt character
                // Add it to the variant
                add_alt_allele(variant, character);
            }
        }
        
        // Output the created VCF variant.
        std::cout << variant << std::endl;
    }

    return 0;
}



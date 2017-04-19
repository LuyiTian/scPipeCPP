// sc_detect_bc.cpp

#include "detect_barcode.h"
#include "timing.h"

int main(int argc, char* argv[]) {
    if (argc < 4) 
    { 
        std::cout << "Usage: sc_detect_bc.  \n" <<\
            "to find barcode sequence enrichment when barcode is unknown a priori, like in Drop-seq \n" <<\
            "\t-I <infq> the input fastq file, this should be the output of ./sc_trim_barcode (required)\n"<<\
            "\t-O <outcsv> the output csv file, which can be used as the cell barcode annotation (required)\n"<<\
            "\t-L <bc_len> the total length of barcode (required)\n"<<\
            "\t-S <surfix> the surfix added in cell name (default: `CELL_`)\n"<<\
            "\t-R <max_reads> the maximum number of reads processed (default: 1000000)\n"<<\
            "\t-N <min_count> barcode that have smaller count will be discarded (default: 10)\n"<<\
            "\t-M <max_mismatch> maximum mismatch allowed when merging barcode (default: 1)\n";

        exit(0);
    } 
    else 
    {
        std::string infq, outcsv, surfix;
        surfix = "CELL_";
        int max_mismatch = 1;
        int bc_len = 0;
        int max_reads = 1000000;
        int min_count = 10;
        for (int i = 1; i < argc; i++) 
        {
            std::string arg = argv[i];

            if (i + 1 != argc) 
            {
                if (arg == "-I") 
                {
                    infq = argv[i + 1];
                } 
                else if (arg == "-O") 
                {
                    outcsv = argv[i + 1];
                }
                else if (arg == "-S") 
                {
                    surfix = argv[i + 1];
                }
                else if (arg == "-L")
                {
                    bc_len = std::atoi(argv[i + 1]);
                }
                else if (arg == "-R")
                {
                    max_reads = std::atoi(argv[i + 1]);
                }
                else if (arg == "-N")
                {
                    min_count = std::atoi(argv[i + 1]);
                }
                else if (arg == "-M")
                {
                    max_mismatch = std::atoi(argv[i + 1]);
                }
            }
        }
        std::cout << "######### detect barcode:" << std::endl;
        std::cout << "parameters:" << std::endl;
        std::cout << "\t" << "out csv file: " << outcsv << std::endl;
        std::cout << "\t" << "input fastq file: " << infq << std::endl;
        std::cout << "\t" << "barcode length: " << bc_len << std::endl;
        std::cout << "\t" << "max reads: " << max_reads << std::endl;
        std::cout << "\t" << "min count: " << min_count << std::endl;
        std::cout << "\t" << "max mismatch " << max_mismatch << std::endl;

        std::unordered_map<std::string, int> counter = summarize_barcode(infq, bc_len, max_reads, max_mismatch, min_count);
        write_barcode_summary(outcsv, surfix, counter);
        return 0;
    }
}
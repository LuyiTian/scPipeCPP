//detect_barcode.cpp

#include "utils.h"

void summarize_barcode(std::string fn, std::string outfn, int bc_len, int max_reads, int max_mismatch)
{
    BGZF *fp = bgzf_open(fn.c_str(), "r");  // input file
    std::ofstream o_file(outfn);  // output file
    if (max_mismatch<1000)
    {
        std::cerr << "max_mismatch should be larger than 1000." << std::endl;
        exit(1);
    }
    if (bc_len < 4)
    {
        std::cerr << "bc_len should be larger than 3." << std::endl;
        exit(1);
    }

    kseq_t *seq1;
    seq =  kseq_init(fq1);
    int cnt = 0;
    while((cnt < max_mismatch) && (kseq_read(seq)>= 0))
    {
        
    }
}
//detect_barcode.cpp

#include "utils.h"

void merge_barcode(std::unordered_map<std::string, int> &counter, int max_mismatch, int min_count)
{
    for (auto bc1 = counter.begin(); bc1 != counter.end();) // use normal iterator in order to use `erase`
    {
        found = false;
        for (auto const& bc2: counter) // use range based
        {
            if (hamming_distance(bc1->first, bc2.first) == 1)
            {
                if (bc1->second == max_mismatch || bc1->second < bc2.second*0.5)
                {
                    found = true;
                    // merge two barcodes
                    counter[bc2.first] += counter[bc1->first];
                    if (__DEBUG){std::cout << "merge: " <<  bc1->first << "::" << bc2.first << "\t" << bc1->second << "::" << bc2.second << std::endl;}
                    break;
                }
            }
        }
        if (found)
        {
            // delete bc1
            bc1 = counter.erase(bc1);
        }
        else
        {
            bc1++;
        }
    }

    for (auto bc1 = counter.begin(); bc1 != counter.end();) // use normal iterator in order to use `erase`
    {
        if (bc1->second < min_count)
        {
            bc1 = counter.erase(bc1);
        }
        else
        {
            bc1++;
        }
    }
}



std::unordered_map<std::string, int> summarize_barcode(std::string fn, int bc_len, int max_reads, int max_mismatch, int min_count)
{
    BGZF *fp = bgzf_open(fn.c_str(), "r");  // input file
    std::unordered_map<std::string, int> counter;
    std::string tmp_bc;
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
        tmp_bc = string(seq->name.s).substr(0, bc_len);
        if (counter.find(tmp_bc) != counter.end())
        {
            counter[tmp_bc] ++;
        }
        else
        {
            counter[tmp_bc] = 1;
        }
    }

    merge_barcode(counter, max_mismatch, min_count);



    return counter;
}

void write_barcode_summary(std::string outfn, std::string surfix, std::unordered_map<std::string, int> counter)
{
    std::ofstream o_file(outfn);  // output file
    int cnt = 0;
    int dig = std::to_string(counter.size()).length()+1;
    for (auto const& bc: counter)
    {
        o_file << surfix << padding(cnt, dig) << "," << bc.first << "," << bc.second << std::endl;
        cnt++;
    }
}
    
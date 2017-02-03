#include <iostream>
#include <fstream>
#include <cstring>
#include <unordered_map>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include "utils.h"
#include "cellbarcode.h"

#ifndef PARSEBAM_H
#define PARSEBAM_H

// overall_count_stat: count statistics for all cells
// chr_aligned_stat: per chromosome count
// per cell statistics:
// cell_mapped: mapped to exon
// cell_mapped_intron: mapped to intron
// cell_mapped_ambiguous: ambiguous map to multiple exon
// cell_align_unmapped: aligned to genome but not mapped to any feature
// cell_unaligned: not aligned
//
class Bamdemultiplex
{
public:
    Barcode bar;
    std::string c_tag;
    std::string m_tag;
    std::string g_tag;
    std::string a_tag;
    std::string out_dir;
    std::string mt_tag;

    std::unordered_map<std::string, int> overall_count_stat;
    std::unordered_map<std::string, int> chr_aligned_stat;
    std::unordered_map<std::string, int> cell_mapped_exon;
    std::unordered_map<std::string, int> cell_mapped_intron;
    std::unordered_map<std::string, int> cell_mapped_ambiguous;
    std::unordered_map<std::string, int> cell_align_unmapped;
    std::unordered_map<std::string, int> cell_unaligned;
    std::unordered_map<std::string, int> cell_ERCC;
    std::unordered_map<std::string, int> cell_MT;

    Bamdemultiplex(std::string odir, Barcode b, std::string cellular_tag, std::string molecular_tag, std::string gene_tag, std::string map_tag, std::string MT_tag);
    ~Bamdemultiplex();
    int barcode_demultiplex(std::string bam_path, int max_mismatch);

    void write_statistics(std::string overall_stat_f, std::string chr_stat_f, std::string cell_stat_f);

    /* data */
};

#endif
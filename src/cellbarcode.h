#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include "utils.h"

#ifndef CELLBARCODE_H
#define CELLBARCODE_H

// a class that stores cellular barcode annotation and
// find close barcode for a given sequence
class Barcode
{
public:
    std::unordered_map<std::string, std::string> barcode_dict;
    std::vector<std::string> cellid_list;
    std::vector<std::string> barcode_list;

    // if annotation is given
    void read_anno(std::string fn);

    std::unordered_map<std::string, std::ofstream> get_count_file_w(std::string out_dir);

    std::unordered_map<std::string, std::ifstream> get_count_file_r(std::string in_dir);

    std::string get_closest_match(std::string bc_seq, int max_mismatch);

    friend std::ostream& operator<< (std::ostream& out, const Barcode& obj); 
};

#endif
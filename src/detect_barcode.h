#include <zlib.h> // for reading compressed .fq file
#include <string>
#include <stdio.h>
#include <iostream>
#include "config_hts.h"
#include "utils.h"

#ifndef DETECTBARCODE_H
#define DETECTBARCODE_H

std::unordered_map<std::string, int> summarize_barcode(std::string fn, int bc_len, int max_reads, int max_mismatch);
#endif
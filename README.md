<<<<<<< HEAD
# scPipeCPP
a C++ pipeline for scRNA-seq processing
=======
[![Standard](https://img.shields.io/badge/c%2B%2B-11-blue.svg)](https://en.wikipedia.org/wiki/C%2B%2B#Standardization)

# scPipe
This is a fork of the C++ portion of the scPipe software from https://github.com/LuyiTian/scPipe

This pipeline provides tools for summarising gene counts from single or paired end reads containing unique molecular identifiers (UMI). Gene counting takes into account the UMI such that genes are only counted once per UMI for each exon.

**Currently this software only supports macOS and Linux.**

## Installation
```
make
make install
```
## Getting started
This pipeline contains 4 programs:

* sc_demultiplex
* sc_gene_counting
* sc_trim_barcode
* sc_exon_mapping

Run any of the programs with no arguments to see usage instructions.

## Dependencies
* [htslib](https://github.com/samtools/htslib) (Included with scPipe)
* [zlib](https://github.com/madler/zlib) (Should be present on all modern macOS or Linux versions)

## Acknowledgements
This package is inspired by `scater` and `scran` packages. Output from this pipeline can be converted to `SCESet` class in `scater` for downstream analysis.
>>>>>>> scPipe/master

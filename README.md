

# scPipeCPP
[![Standard](https://img.shields.io/badge/c%2B%2B-11-blue.svg)](https://en.wikipedia.org/wiki/C%2B%2B#Standardization)

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

after `make install`. you can use `make test` to test the package. some test data and script will be in the `test` folder.

Run any of the programs with no arguments to see usage instructions.

## Dependencies
* [zlib](https://github.com/madler/zlib) (Should be present on all modern macOS or Linux versions)

## TODO

* update UMI correction methods
* update `basic_test.cpp` and remove local path dependencies.
* simulate Drop-seq and MAR-seq data. and use them to test our pipeline.

## Acknowledgements
This package is inspired by `scater` and `scran` packages. Output from this pipeline can be converted to `SCESet` class in `scater` for downstream analysis.

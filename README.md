# scPipe
This is a fork of the C++ portion of the scPipe software from https://github.com/LuyiTian/scPipe

This pipeline provides tools for summarising gene counts from single or paired end reads containing unique molecular identifiers (UMI). Gene counting takes into account the UMI such that genes are only counted once per UMI for each exon.

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

## Acknowledgements
This package is inspired by `scater` and `scran` packages. Output from this pipeline can be converted to `SCESet` class in `scater` for downstream analysis.

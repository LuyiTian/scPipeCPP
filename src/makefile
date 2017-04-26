################################################################
## scPipe Makefile
## Version 0.1
## Author: Luyi Tian (tian.l@wehi.edu.au)
################################################################

DIR = $(shell pwd)
prefix ?= /usr/local

INSTALL = install -p
INSTALL_PROGRAM = $(INSTALL)

CPPFLAGS = -std=c++11 -O2 -pthread -static-libstdc++
DEVEL_FLAGS = -std=c++11 -O2 -pthread -Wall -Wno-unused-function
CFLAGS = -O2
LIBS = -lz

HTS_OBJS = sam.hts.o bgzf.hts.o hfile.hts.o hts.hts.o kstring.hts.o
OBJS = utils.o cellbarcode.o parsebam.o parsecount.o trimbarcode.o transcriptmapping.o detect_barcode.o fq_simulator.o timing.o
RM = rm -rf

.PHONY: all clean buildtest runtest test buildPipe install devel
.SUFFIXES: .o .hts.o

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(LIBS) -o $@ $<

%.hts.o: %.c
	$(CC) -c $(CFLAGS) -o $@ $<

all: buildPipe buildtest

buildtest:
	$(CXX) -o test/test_simulator test/test_simulator.cpp $(LIBS) $(OBJS) $(HTS_OBJS) $(CPPFLAGS)
	#$(CXX) -o test/basic_test test/basic_test.cpp $(LIBS) $(OBJS) $(HTS_OBJS) $(CPPFLAGS)

testpipe:
	cd test && chmod u+x ./test_celseq2_pipe.sh && ./test_celseq2_pipe.sh 
	cd test && chmod u+x ./test_dropseq_pipe.sh && ./test_dropseq_pipe.sh

runtest:
	cd test && ./test_simulator
	#./test/basic_test

test: runtest

buildPipe: $(HTS_OBJS) $(OBJS)
	mkdir -p bin
	$(CXX) -o bin/sc_demultiplex sc_demultiplex.cpp $(LIBS) $(OBJS) $(HTS_OBJS) $(CPPFLAGS)
	$(CXX) -o bin/sc_gene_counting sc_gene_counting.cpp $(LIBS) $(OBJS) $(HTS_OBJS) $(CPPFLAGS)
	$(CXX) -o bin/sc_trim_barcode sc_trim_barcode.cpp $(LIBS) $(OBJS) $(HTS_OBJS) $(CPPFLAGS)
	$(CXX) -o bin/sc_exon_mapping sc_exon_mapping.cpp $(LIBS) $(OBJS) $(HTS_OBJS) $(CPPFLAGS)
	$(CXX) -o bin/sc_detect_bc sc_detect_bc.cpp $(LIBS) $(OBJS) $(HTS_OBJS) $(CPPFLAGS)

devel: $(HTS_OBJS) $(OBJS)
	mkdir -p bin
	$(CXX) -o bin/sc_demultiplex sc_demultiplex.cpp $(LIBS) $(OBJS) $(HTS_OBJS) $(DEVEL_FLAGS)
	$(CXX) -o bin/sc_gene_counting sc_gene_counting.cpp $(LIBS) $(OBJS) $(HTS_OBJS) $(DEVEL_FLAGS)
	$(CXX) -o bin/sc_trim_barcode sc_trim_barcode.cpp $(LIBS) $(OBJS) $(HTS_OBJS) $(DEVEL_FLAGS)
	$(CXX) -o bin/sc_exon_mapping sc_exon_mapping.cpp $(LIBS) $(OBJS) $(HTS_OBJS) $(DEVEL_FLAGS)
	$(CXX) -o bin/sc_detect_bc sc_detect_bc.cpp $(LIBS) $(OBJS) $(HTS_OBJS) $(CPPFLAGS)

clean:
	$(RM) $(OBJS)
	$(RM) $(HTS_OBJS)
	$(RM) bin/*
	$(RM) test/basic_test
	$(RM) test/test_simulator
	$(RM) test/celseq2/*
	$(RM) test/dropseq/*

install:
	$(INSTALL_PROGRAM) bin/sc_demultiplex $(prefix)/bin/sc_demultiplex
	$(INSTALL_PROGRAM) bin/sc_gene_counting $(prefix)/bin/sc_gene_counting
	$(INSTALL_PROGRAM) bin/sc_trim_barcode $(prefix)/bin/sc_trim_barcode
	$(INSTALL_PROGRAM) bin/sc_exon_mapping $(prefix)/bin/sc_exon_mapping
	$(INSTALL_PROGRAM) bin/sc_detect_bc $(prefix)/bin/sc_detect_bc

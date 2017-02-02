################################################################
## scPipe Makefile
## Version 0.1
## Author: Luyi Tian (tian.l@wehi.edu.au)
################################################################

DIR = $(shell pwd)
HTS_DIR = $(DIR)/hts
INSTALL_PATH ?= /usr/local

CPPFLAGS = -std=c++11 -static-libstdc++ -O2

LIBS = -I$(HTS_DIR)/include -L$(HTS_DIR)/lib -Wl,-rpath,$(HTS_DIR)/lib -lhts -lz

OBJS = utils.o cellbarcode.o parsebam.o parsecount.o trimbarcode.o transcriptmapping.o
RM = rm -rf

.PHONY: all buildhts buildtest runtest test buildPipe install

%.o: %.cpp $(HTS_DIR)
	$(CXX) -c $(CPPFLAGS) $(LIBS) -o $@ $<

all: buildPipe

$(HTS_DIR):
	mkdir -p $(HTS_DIR) &&\
	cd htslib-1.3 && \
	make clean && \
	./configure --disable-libcurl --prefix=$(HTS_DIR) CC=$(CC) CXX=$(CXX) && \
	make && \
	make install && \
	cd ..

buildtest: $(OBJS)
	$(CXX) -o test/basic_test test/basic_test.cpp $(LIBS) $(OBJS) $(CPPFLAGS)

runtest:
	./test/basic_test

test: buildtest runtest

buildhts: $(HTS_DIR)

buildPipe: $(OBJS)
	mkdir -p bin
	$(CXX) -o bin/sc_demultiplex sc_demultiplex.cpp $(LIBS) $(OBJS) $(CPPFLAGS)
	$(CXX) -o bin/sc_gene_counting sc_gene_counting.cpp $(LIBS) $(OBJS) $(CPPFLAGS)
	$(CXX) -o bin/sc_trim_barcode sc_trim_barcode.cpp $(LIBS) $(OBJS) $(CPPFLAGS)
	$(CXX) -o bin/sc_exon_mapping sc_exon_mapping.cpp $(LIBS) $(OBJS) $(CPPFLAGS)

clean:
	$(RM) $(OBJS)
	$(RM) test/basic_test
	$(RM) ./bin
	$(RM) hts

install:
	cp bin/sc_demultiplex $INSTALL_PATH/bin/
	cp bin/sc_gene_counting $INSTALL_PATH/bin/
	cp bin/sc_trim_barcode $INSTALL_PATH/bin/
	cp bin/sc_exon_mapping $INSTALL_PATH/bin/

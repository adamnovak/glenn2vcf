.PHONY: all clean

CXX=g++
INCLUDES=-Iekg/vg/src -Iekg/vg/include
CXXFLAGS=-O3 -std=c++11 -fopenmp -g $(INCLUDES)
VGLIBDIR=ekg/vg/lib
LDSEARCH=-Lekg/vg/src -L$(VGLIBDIR)
LDFLAGS=-lm -lpthread -lz -lbz2 -lsnappy -ldivsufsort -ldivsufsort64 -ljansson $(LDSEARCH)
LIBVG=ekg/vg/lib/libvg.a
LIBXG=ekg/vg/lib/libxg.a
LIBPROTOBUF=ekg/vg/lib/libprotobuf.a
LIBSDSL=ekg/vg/lib/libsdsl.a
LIBGSSW=ekg/vg/lib/libgssw.a
LIBSNAPPY=ekg/vg/lib/libsnappy.a
LIBROCKSDB=ekg/vg/lib/librocksdb.a
LIBHTS=ekg/vg/lib/libhts.a
LIBGCSA2=ekg/vg/lib/libgcsa2.a
LIBVCFLIB=ekg/vg/lib/libvcflib.a
LIBRAPTOR=ekg/vg/lib/libraptor2.a
LIBGFAKLUGE=ekg/vg/lib/libgfakluge.a
LIBSUPBUB=ekg/vg/lib/libsupbub.a
VGLIBS=$(LIBVG) $(LIBXG) $(LIBVCFLIB) $(LIBGSSW) $(LIBSNAPPY) $(LIBROCKSDB) $(LIBHTS) $(LIBGCSA2) $(LIBPROTOBUF) $(LIBRAPTOR) $(LIBGFAKLUGE) $(LIBSUPBUB) $(LIBSDSL)

#Some little adjustments to build on OSX
#(tested with gcc4.9 and jansson installed from MacPorts)
SYS=$(shell uname -s)
ifeq (${SYS},Darwin)
	LDFLAGS:=$(LDFLAGS) -L/opt/local/lib/ # needed for macports jansson
else
	LDFLAGS:=$(LDFLAGS) -lrt
endif

all: glenn2vcf

$(LIBSDSL): $(LIBVG)

$(LIBPROTOBUF): $(LIBVG)

$(LIBXG): $(LIBVG)

$(LIBVG): ekg/vg/src/*
	cd ekg/vg && . ./source_me.sh && $(MAKE)

main.o: $(LIBVG)

glenn2vcf: main.o $(LIBSONLIB) $(VGLIBS) 
	$(CXX) $^ -o $@ $(CXXFLAGS) $(LDFLAGS)

clean:
	rm -f glenn2vcf
	rm -f *.o
	cd ekg/vg && $(MAKE) clean

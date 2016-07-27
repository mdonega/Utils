# Basic Makefile for ROOT
# Code in waves.C -> int main()

CXX=`root-config --cxx`
CXXFLAGS=`root-config --cflags`
LDFLAGS=`root-config --ldflags`
LDLIBS=`root-config --glibs`
ROOTLIBS='-lRooFit -lHtml -lMinuit -lRooFitCore -lRooStats -lHistFactory -lCore -lCint -lGraf -lGraf3d -lHist -lMatrix -lPostscript -lProof -lTree -lGpad -lGui'


all: waves

% : %.C
g++ -g -Wall `root-config --cflags --libs`-L$(ROOTSYS)/lib $(ROOTLIBS) $? -o $@

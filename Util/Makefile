CXX = g++

ROOTFLAGS = `root-config --cflags`
ROOTLIBS  = `root-config --libs`

all: $(shell ls *.cpp | sed "s/\.cpp/\.o/g")

.cpp.o:
	$(CXX) $(ROOTFLAGS) -c $< -o $@
Style.o       : const.h Style.h Style.cpp

.PHONY: clean
clean:
	rm $(shell ls *.cpp | sed "s/\.cpp/\.o/g")

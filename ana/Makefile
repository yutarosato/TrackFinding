CXX = g++

ROOTFLAGS = `root-config --cflags`
ROOTLIBS  = `root-config --libs`

all: util $(shell ls *.cpp | sed "s|\.cpp||g")
other = $(shell ls ../Util/*.cpp | sed "s|\.cpp|.o|g")

util:
	cd ../Util; make

tracking:tracking.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

tracking2:tracking2.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

.cpp.o:
	$(CXX) $(ROOTFLAGS) -c $< -o $@
tracking.o        : setting.h tracking.cpp     $(other)
tracking2.o        : setting.h tracking2.cpp     $(other)

.PHONY: clean
clean:
	rm $(shell ls *.cpp | sed "s|\.cpp|\.o|g") $(shell ls *.cpp | sed "s|\.cpp||g")

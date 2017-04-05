# target
#TARGET = main # main_v0.1 main_v0.1.2 main_v0.1.2.2 main_v0.2
#TARGET = main_v0.1
#TARGET = main_v0.1.2
#TARGET = main_v0.1.2.2
#TARGET = main_v0.2
TARGET = main_my

# suffix
ObjSuf          = .o
SrcSuf          = .cxx
ExeSuf          =

# set suffixs
TARGETOBJ = $(TARGET)$(ObjSuf)
TARGETSRC = $(TARGET)$(SrcSuf)
TARGETEXE = $(TARGET)$(ExeSuf)

# for ROOT
ROOTCFLAGS      = $(shell root-config --cflags)
ROOTLIBS        = $(shell root-config --libs)
ROOTGLIBS       = $(shell root-config --glibs)

# for linux
CXX     = g++ -g
CXXFLAGS= -O2 -Wall
LD      = g++ -g
LDFLAGS = -O2 -Wall

# create parameters
CXXFLAGS += $(ROOTCFLAGS)
LIBS      = $(ROOTLIBS) $(SYSLIBS)
GLIBS     = $(ROOTGLIBS) $(SYSLIBS)

all: $(TARGETEXE)
OBJ = root.o io.o algo.o Display.o

$(TARGETEXE): $(TARGETOBJ) $(OBJ)
	$(LD) $(LDFLAGS) $(LIBS) -o $(TARGETEXE) $(TARGETOBJ) $(OBJ)

$(TARGETOBJ): $(TARGETSRC)
	$(CXX) $(CXXFLAGS) -c -o $(TARGETOBJ) $(TARGETSRC)

root.o : root.cxx
	$(CXX) $(CXXFLAGS) -c -o root.o root.cxx

io.o : io.cxx
	$(CXX) $(CXXFLAGS) -c -o io.o io.cxx

algo.o : algo.cxx
	$(CXX) $(CXXFLAGS) -c -o algo.o algo.cxx

Display.o : Display.cxx
	$(CXX) $(CXXFLAGS) -c -o Display.o Display.cxx

clean:
	rm -f *.o
	rm -f *~


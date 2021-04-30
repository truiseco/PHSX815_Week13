# C++ compiler
CXX = g++

# necessary compiler flags for using ROOT (root.cern.ch) - remove these if you're not using root
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTGLIBS     = $(shell root-config --glibs)

# ROOT shared library flags
#GLIBS = $(filter-out -stdlib=libc++ -pthread , $(ROOTGLIBS))
GLIBS = $(ROOTGLIBS)

# some compiler flags
CXXFLAGS = -std=c++17 -g -Wall -O2
# ROOT flags
CXXFLAGS += -fPIC $(filter-out -stdlib=libc++ -pthread , $(ROOTCFLAGS))
#CXXFLAGS = $(ROOTCFLAGS)

# location of source code
SRCDIR = ./src/

#location of header files
INCLUDEDIR = ./include/

CXXFLAGS += -I$(INCLUDEDIR)

# location of object files (from compiled library files)
OUTOBJ = ./obj/

CPP_FILES := $(wildcard src/*.cpp)
H_FILES := $(wildcard include/*.h)
OBJ_FILES := $(addprefix $(OUTOBJ),$(notdir $(CPP_FILES:.cpp=.o)))

# targets to make
all: KMeans.x

# recipe for building KMeans.x
KMeans.x:  $(SRCDIR)KMeans.C $(OBJ_FILES) $(H_FILES)
	$(CXX) $(CXXFLAGS) $ $< $(GLIBS) -o KMeans.x
	touch KMeans.x

$(OUTOBJ)%.o: src/%.cpp include/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# clean-up target (make clean)
clean:
	rm -f *.x
	rm -rf *.dSYM
	rm -f $(OUTOBJ)*.o

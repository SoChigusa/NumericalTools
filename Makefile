
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs) -lMinuit2

# for iMac
CXX = icpc
CXXFLAGS = -g -Wall -std=c++14 -Ofast
NUMERICALFLAGS = -I/Users/SoChigusa/works/NumericalTools/

# for laptop
# CXX := g++
# CXXFLAGS := -g -Wall -std=c++11 -O2

all:		test.out

clean:
		$(RM) *.out

%.out: 	%.cpp
		$(CXX) $(CXXFLAGS) $(NUMERICALFLAGS) -o $@ $^ $(ROOTLIBS) $(ROOTFLAGS)

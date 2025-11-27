# build the example code 
# to build your own program, add the EXE name in the EXES list below
# then include a target section following the example of one of the
# other programs below

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs)
ROOTFLAGS   = $(ROOTCFLAGS) $(ROOTLIBS) $(ROOTGLIBS) 
CXXFLAGS  += $(ROOTCFLAGS) -I$(ODELIB) -Wall -O3
LDFLAGS    = $(ROOTLIBS) $(ROOTGLIBS) # ,-rpath #,$(ODELIB) -L$(ODELIB) -lode
GSLFLAGS   = -lgsl -lgslcblas
GXX	   = g++ $(CXXFLAGS)


EXES=baseball1 vterm baseball2


all: $(EXES)

vterm: vterm.cpp RungeKutta.hpp
	$(GXX) $(CXXFLAGS) RungeKutta.hpp $@.cpp -o$@ $(LDFLAGS) 

baseball1: baseball1.cpp RungeKutta.hpp
	$(GXX) $(CXXFLAGS) RungeKutta.hpp $@.cpp -o$@ $(LDFLAGS) 

baseball2: baseball2.cpp RungeKutta.hpp
	$(GXX) $(CXXFLAGS) RungeKutta.hpp $@.cpp -o$@ $(LDFLAGS) 


clean:
	rm -f *.o *.so *.d *~ $(EXES)



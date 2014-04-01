CXX = $(shell root-config --cxx)
CPPFLAGS = -isystem$(shell root-config --incdir) -I inc
CXXFLAGS = -Wall -Wextra -pedantic -O2 -Wshadow -fPIC $(shell root-config --cflags)
LD = $(shell root-config --ld)
LDFLAGS = $(shell root-config --ldflags)
LDLIBS =  $(shell root-config --glibs) -lMinuit2 -lMathMore

VPATH = inc:src:obj

OBJECTS = obj/TopMass.o

COMPILE = $(CXX) $(CXXFLAGS) $(CPPFLAGS) -c
LINK = $(LD) $(LDFLAGS)
LINKEND = $(OBJECTS) $(LDLIBS)

# add your executable names here
all : DoFit

DoFit : DoFit.o $(OBJECTS)
	$(LINK) -o DoFit DoFit.o $(LINKEND)
DoFit.o: DoFit.C TopMass.h
	$(COMPILE) DoFit.C

clean:
	-rm -f DoFit obj/*.o *.o

obj/TopMass.o : TopMass.C TopMass.h
	$(COMPILE) src/TopMass.C -o obj/TopMass.o

.PHONY : clean


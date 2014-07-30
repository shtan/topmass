CXX = $(shell root-config --cxx)
CPPFLAGS = -isystem$(shell root-config --incdir) -I inc 
CXXFLAGS = -Wall -Wextra -pedantic -O2 -Wshadow -fPIC $(shell root-config --cflags)
LD = $(shell root-config --ld)
LDFLAGS = $(shell root-config --ldflags)
LDLIBS =  $(shell root-config --glibs) -lMinuit2 -lMathMore

VPATH = inc:src:obj

OBJECTS = obj/TopMass.o obj/Mt2Calculator.o obj/Diagnostics.o obj/Shapes.o obj/Plotting.o

COMPILE = $(CXX) $(CXXFLAGS) $(CPPFLAGS) -c
LINK = $(LD) $(LDFLAGS)
LINKEND = $(OBJECTS) $(LDLIBS)

# add your executable names here
all : DoFit

DoFit : DoFit.o $(OBJECTS)
	$(LINK) -o DoFit DoFit.o $(LINKEND)
DoFit.o: DoFit.C TopMass.h Shapes.h Plotting.h
	$(COMPILE) DoFit.C

clean:
	-rm -f DoFit obj/*.o *.o

obj/TopMass.o : TopMass.C TopMass.h Mt2Calculator.h Shapes.h
	$(COMPILE) src/TopMass.C -o obj/TopMass.o

obj/Mt2Calculator.o : Mt2Calculator.C Mt2Calculator.h
	$(COMPILE) src/Mt2Calculator.C -o obj/Mt2Calculator.o

obj/Diagnostics.o : Diagnostics.C TopMass.h Shapes.h
	$(COMPILE) src/Diagnostics.C -o obj/Diagnostics.o

obj/Shapes.o : Shapes.C Shapes.h
	$(COMPILE) src/Shapes.C -o obj/Shapes.o

obj/Plotting.o : Plotting.C Plotting.h Mt2Calculator.h TopMass.h
	$(COMPILE) src/Plotting.C -o obj/Plotting.o

.PHONY : clean


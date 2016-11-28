ROOTCFLAGS    = $(shell /usr/local/bin/root-config --cflags)
ROOTLIBS      = $(shell /usr/local/bin/root-config --libs)
ROOTGLIBS     = $(shell /usr/local/bin/root-config --glibs)

CXX           = g++
CXXFLAGS += -O -g -Wall -fPIC $(ROOTCFLAGS) -I/usr/local/include/
LIBS          = $(ROOTLIBS)
GLIBS         = $(ROOTGLIBS)

EXEC = Generator
OBJS = JGenBeamEnergy.o Generator.o
OBJSL = JGenBeamEnergy.o
LDLIBS =$(ROOTLIBS)


all:    $(OBJS) 
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(OBJS) $(LDLIBS)
$(OBJS): %.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@
clean:
	-rm $(OBJS) $(EXEC)

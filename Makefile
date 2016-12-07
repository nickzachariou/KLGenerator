ROOTCFLAGS    = $(shell /scratch/mbashkan/GlueX/gluex_install/gluex_top/root/root-6.06.08/bin/root-config --cflags)
ROOTLIBS      = $(shell /scratch/mbashkan/GlueX/gluex_install/gluex_top/root/root-6.06.08/bin/root-config --libs)
ROOTGLIBS     = $(shell /scratch/mbashkan/GlueX/gluex_install/gluex_top/root/root-6.06.08/bin/root-config --glibs)

CXX           = g++
CXXFLAGS += -O -g -Wall -fPIC $(ROOTCFLAGS) -I/usr/local/include/ -I$(HALLD_HOME)/src/.Linux_RH-x86_64-gcc4.8.5/libraries/
LIBS          = $(ROOTLIBS) $(HALLD_HOME)/src/.Linux_RH-x86_64-gcc4.8.5/libraries/
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

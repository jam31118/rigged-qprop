
include $(QPROP_HOME)/src/main/variables.mk

# The option -march=native may allow the compiler to use fancy instructions (SSE, AVX and the like).
# It also means you have to recompile on each CPU architecture you are using.
cpu_type = generic
arch_flags = -mtune=$(cpu_type) #-march=$(arch)
gxxoptions = -g -std=c++0x -O3 -Wno-deprecated $(arch_flags) -funroll-all-loops -finline-functions -fexpensive-optimizations -ffast-math -fPIC -Xlinker -defsym -Xlinker MAIN__=main -I. $(extgxxopts)
#gxxoptions = -g -std=c++0x -O8 -Wno-deprecated $(arch_flags) -funroll-all-loops -finline-functions -fexpensive-optimizations -ffast-math -fPIC -Xlinker -defsym -Xlinker MAIN__=main -I. $(extgxxopts)
#gxxoptions = -g -std=c++0x -O8 -Wno-deprecated -march=native -funroll-all-loops -finline-functions -fexpensive-optimizations -ffast-math -fPIC -Xlinker -defsym -Xlinker MAIN__=main -I.
libpath = ../../lib/x86_64

# $@ - means the target
# $^ - means all prerequisites
# $< - means just the first prerequisite

all: install

bar.o: bar.cc
	$(CXX) $(gxxoptions) -c -o $@ $<

ylm.o: ylm.cc
	$(CXX) $(gxxoptions) -c -o $@ $<

wavefunction.o: wavefunction.cc
	$(CXX) $(gxxoptions) -c -o $@ $<

hamop.o: hamop.cc
	$(CXX) $(gxxoptions) -c -o $@ $<

grid.o: grid.cc
	$(CXX) $(gxxoptions) -c -o $@ $<

fluid.o: fluid.cc
	$(CXX) $(gxxoptions) -c -o $@ $<

factorial.o: factorial.cc
	$(CXX) $(gxxoptions) -c -o $@ $<

winop.o: winop.cc
	$(CXX) $(gxxoptions) -c -o $@ $<

smallHelpers.o: smallHelpers.cc
	$(CXX) $(gxxoptions) -c -o $@ $<

parameter.o: parameter.cc
	$(CXX) $(gxxoptions) -c -o $@ $<

vecpot.o: vecpot.cc
	$(CXX) $(gxxoptions) -c -o $@ $<

potentials.o: potentials.cc
	$(CXX) $(gxxoptions) -c -o $@ $<

# ar rvs ---> replace files, be verbose, Write an object-file index into the archive
libqprop.a: bar.o ylm.o wavefunction.o hamop.o grid.o fluid.o factorial.o winop.o parameter.o smallHelpers.o vecpot.o potentials.o
	ar rvs $@ $^

clean:
	rm -f *.o *.a

install: $(libpath)/libqprop.a

# solution with mkdir -p is ugly but it works
$(libpath)/libqprop.a: libqprop.a
	mkdir -p $(libpath); cp libqprop.a $(libpath)


# The option -march=native may allow the compiler to use fancy instructions (SSE, AVX and the like).
# It also means you have to recompile on each CPU architecture you are using.
gxxoptions = -g -std=c++0x -O8 -Wno-deprecated -march=native -funroll-all-loops -finline-functions -fexpensive-optimizations -ffast-math  -Xlinker -defsym -Xlinker MAIN__=main -I.
libpath = ../../lib/x86_64

# $@ - means the target
# $^ - means all prerequisites
# $< - means just the first prerequisite

all: install

bar.o: bar.cc
	g++ $(gxxoptions) -c -o $@ $<

ylm.o: ylm.cc
	g++ $(gxxoptions) -c -o $@ $<

wavefunction.o: wavefunction.cc
	g++ $(gxxoptions) -c -o $@ $<

hamop.o: hamop.cc
	g++ $(gxxoptions) -c -o $@ $<

grid.o: grid.cc
	g++ $(gxxoptions) -c -o $@ $<

fluid.o: fluid.cc
	g++ $(gxxoptions) -c -o $@ $<

factorial.o: factorial.cc
	g++ $(gxxoptions) -c -o $@ $<

winop.o: winop.cc
	g++ $(gxxoptions) -c -o $@ $<

parameter.o: parameter.cc
	g++ $(gxxoptions) -c -o $@ $<

# ar rvs ---> replace files, be verbose, Write an object-file index into the archive
libqprop.a: bar.o ylm.o wavefunction.o hamop.o grid.o fluid.o factorial.o winop.o parameter.o
	ar rvs $@ $^

clean:
	rm -f *.o *.a

install: $(libpath)/libqprop.a

# solution with mkdir -p is ugly but it works
$(libpath)/libqprop.a: libqprop.a
	mkdir -p $(libpath); cp libqprop.a $(libpath)

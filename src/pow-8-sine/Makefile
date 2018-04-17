include ../common-variables-for-make.mk

all: hydrogen_im hydrogen_re eval-tsurff

hydrogen_im : libqprop.a hydrogen_im.cc
	$(CXX) $(optimargs) -std=c++0x hydrogen_im.cc -o hydrogen_im -I. -I../base -lqprop -lm -L../../lib/x86_64

hydrogen_re : libqprop.a hydrogen_re.cc
	$(CXX) $(optimargs) -std=c++0x hydrogen_re.cc -o hydrogen_re -I. -I../base -lqprop -lm $(gslargs) -L../../lib/x86_64  

eval-tsurff : libqprop.a eval-tsurff.cc
	$(CXX) $(optimargs) -std=c++0x eval-tsurff.cc -o eval-tsurff -I. -I../base -lqprop -lm $(gslargs) -L../../lib/x86_64

eval-tsurff-mpi : libqprop.a eval-tsurff-mpi.cc
	$(CXX) $(optimargs) -std=c++0x eval-tsurff-mpi.cc -o eval-tsurff-mpi -I. -I../base -lqprop -lm  $(gslargs) $(mpiargs) -L../../lib/x86_64

libqprop.a:
	cd ../base/ && make 

clean:
	cd ../base/ && make clean
	rm -f ../../lib/x86_64/libqprop.a

include $(QPROP_SRC)/common-variables-for-make.mk

all: hydrogen_im hydrogen_re eval-tsurff

hydrogen_im : libqprop.a hydrogen_im.cc
	$(CXX) $(optimargs) -std=c++0x hydrogen_im.cc -o hydrogen_im -I. -I$(QPROP_SRC)/base -lqprop -lm -L$(QPROP_SRC)/../lib/x86_64

hydrogen_re : libqprop.a hydrogen_re.cc
	$(CXX) $(optimargs) -std=c++0x hydrogen_re.cc -o hydrogen_re -I. -I$(QPROP_SRC)/base -lqprop -lm $(gslargs) -L$(QPROP_SRC)/../lib/x86_64 $(boostargs) -DHAVE_BOOST 

eval-tsurff : libqprop.a eval-tsurff.cc
	$(CXX) $(optimargs) -std=c++0x eval-tsurff.cc -o eval-tsurff -I. -I$(QPROP_SRC)/base -lqprop -lm $(gslargs) -L$(QPROP_SRC)/../lib/x86_64

eval-tsurff-mpi : libqprop.a eval-tsurff-mpi.cc
	$(CXX) $(optimargs) -std=c++0x eval-tsurff-mpi.cc -o eval-tsurff-mpi -I. -I$(QPROP_SRC)/base -lqprop -lm  $(gslargs) $(mpiargs) -L$(QPROP_SRC)/../lib/x86_64

libqprop.a:
	cd $(QPROP_SRC)/base/ && make 

clean:
	cd $(QPROP_SRC)/base/ && make clean
	rm -f $(QPROP_SRC)/lib/x86_64/libqprop.a


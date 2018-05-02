include $(QPROP_SRC)/common-variables-for-make.mk

imag_bin = imag-prop
real_bin = real-prop
eval_bin = eval-tsurff
eval_mpi_bin = eval-tsurff-mpi

all: $(imag_bin) $(real_bin) $(eval_bin)

$(imag_bin) : libqprop.a imag-prop.cc
	$(CXX) $(optimargs) -std=c++0x imag-prop.cc -o $@ -I. -I$(QPROP_SRC)/base -lqprop -lm -L$(QPROP_SRC)/../lib/x86_64

$(real_bin) : libqprop.a real-prop.cc
	$(CXX) $(optimargs) -std=c++0x real-prop.cc -o $@ -I. -I$(QPROP_SRC)/base -lqprop -lm $(gslargs) -L$(QPROP_SRC)/../lib/x86_64 $(boostargs) -DHAVE_BOOST 

$(eval_bin) : libqprop.a eval-tsurff.cc
	$(CXX) $(optimargs) -std=c++0x eval-tsurff.cc -o $@ -I. -I$(QPROP_SRC)/base -lqprop -lm $(gslargs) -L$(QPROP_SRC)/../lib/x86_64

$(eval_mpi_bin) : libqprop.a eval-tsurff-mpi.cc
	$(CXX) $(optimargs) -std=c++0x eval-tsurff-mpi.cc -o $@ -I. -I$(QPROP_SRC)/base -lqprop -lm  $(gslargs) $(mpiargs) -L$(QPROP_SRC)/../lib/x86_64

libqprop.a:
	cd $(QPROP_SRC)/base/ && make 

clean:
	cd $(QPROP_SRC)/base/ && make clean
	rm -f $(QPROP_SRC)/lib/x86_64/libqprop.a

clean-bin: clean
	rm -f $(imag_bin) $(real_bin) $(eval_bin) $(eval_mpi_bin)

clean-data:
	rm -f *.dat *.raw *.log *.bin

clean-all: clean clean-bin clean-data


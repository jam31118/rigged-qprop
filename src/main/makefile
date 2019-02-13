QPROP_SRC = $(QPROP_HOME)/src
QPROP_BIN = $(QPROP_HOME)/bin

include $(QPROP_SRC)/main/variables.mk
include $(CU_PROP_HOME)/src/var.mk

imag_bin = imag-prop
real_bin = real-prop
eval_bin = eval-tsurff
eval_mpi_bin = eval-tsurff-mpi

ppp_bin = ppp
ppp_mpi_bin = ppp-mpi
ppp_gpu_bin = ppp-gpu

bin_list = $(imag_bin) $(real_bin) $(eval_bin) $(eval_mpi_bin) $(ppp_bin) $(ppp_mpi_bin) 
ext_bin_list = $(ppp_gpu_bin)
total_bin_list = $(bin_list) $(ext_bin_list)

MPICXX ?= $(CXX)

all: install

$(imag_bin) : libqprop.a imag-prop.cc imag-prop-main.cc
	$(CXX) $(optimargs) -std=c++0x imag-prop-main.cc imag-prop.cc -o $@ -I. -I$(QPROP_SRC)/base -L$(QPROP_SRC)/../lib/x86_64 -lqprop -lm 

$(real_bin) : libqprop.a real-prop.cc real-prop-main.cc
	$(CXX) $(optimargs) -std=c++0x real-prop-main.cc real-prop.cc -o $@ -I. -I$(QPROP_SRC)/base -lqprop -lm -I$(GSLHOME)/include -L$(QPROP_SRC)/../lib/x86_64 -I$(BOOSTHOME)/include -DHAVE_BOOST 

$(eval_bin) : libqprop.a eval-tsurff.cc eval-tsurff-main.cc
	$(CXX) $(optimargs) -std=c++0x eval-tsurff-main.cc eval-tsurff.cc -o $@ -I. -I$(QPROP_SRC)/base -lqprop -lm $(gslargs) -L$(QPROP_SRC)/../lib/x86_64

$(eval_mpi_bin) : libqprop.a eval-tsurff-mpi.cc eval-tsurff-mpi-main.cc
	$(MPICXX) $(optimargs) -std=c++0x eval-tsurff-mpi-main.cc eval-tsurff-mpi.cc -o $@ -I. -I$(QPROP_SRC)/base -lqprop -lm  $(gslargs) $(mpiargs) -DHAVE_MPI -L$(QPROP_SRC)/../lib/x86_64

$(ppp_bin): ppp.cc libqprop.a libmatrix.a 
	$(CXX) -std=c++0x $< -o $@ $(optimargs) $(QPROP_FLAGS) $(MATRIX_FLAGS) -DHAVE_BOOST

$(ppp_mpi_bin): ppp.cc libqprop.a libmatrix.a
	$(CXX) -std=c++0x $< -o $@ $(optimargs) $(QPROP_FLAGS) $(MATRIX_FLAGS) $(mpiargs) -DHAVE_MPI -DHAVE_BOOST

$(ppp_gpu_bin): ppp.cc libqprop.a libmatrix.a libcu_propagator.so
	$(CXX) -std=c++0x $< -o $@ $(optimargs) $(QPROP_FLAGS) $(MATRIX_FLAGS) $(CU_PROP_FLAGS) $(CUDA_LIB_FLAGS) -DHAVE_CUDA -DHAVE_BOOST
	$(MAKE) move-bin

install: $(bin_list)
	$(MAKE) move-bin

move-bin:
	if [ ! -d "$(QPROP_BIN)" ]; then mkdir $(QPROP_BIN); fi
	for bin in $(total_bin_list); do if [ -f "$$bin" ]; then mv $$bin $(QPROP_BIN); fi; done

libqprop.a:
	+$(MAKE) -C $(QPROP_SRC)/base/ 

libmatrix.a:
	+$(MAKE) install -C $(MATRIX_HOME)

libcu_propagator.so:
	+$(MAKE) install -C $(CU_PROP_HOME)

clean: clean-bin
	cd $(QPROP_SRC)/base/ && $(MAKE) clean
	rm -f $(QPROP_SRC)/lib/x86_64/libqprop.a

clean-bin:
	for bin_name in $(total_bin_list); do if [ -f "$$bin_name" ]; then rm $$bin_name; fi; done
	if [ -d "$(QPROP_BIN)" ]; then cd $(QPROP_BIN) && rm -f $(total_bin_list); fi

clean-data:
	rm -f *.dat *.raw *.log *.bin

clean-all: clean clean-bin clean-data


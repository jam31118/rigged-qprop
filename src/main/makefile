QPROP_SRC = $(QPROP_HOME)/src
QPROP_BIN = $(QPROP_HOME)/bin

include $(QPROP_SRC)/main/variables.mk

imag_bin = imag-prop
real_bin = real-prop
eval_bin = eval-tsurff
eval_mpi_bin = eval-tsurff-mpi

bin_list = $(imag_bin) $(real_bin) $(eval_bin) $(eval_mpi_bin)

MPICXX ?= g++

all: install

$(imag_bin) : libqprop.a imag-prop.cc imag-prop-main.cc
	$(CXX) $(optimargs) -std=c++0x imag-prop-main.cc imag-prop.cc -o $@ -I. -I$(QPROP_SRC)/base -lqprop -lm -L$(QPROP_SRC)/../lib/x86_64

$(real_bin) : libqprop.a real-prop.cc real-prop-main.cc
	$(CXX) $(optimargs) -std=c++0x real-prop-main.cc real-prop.cc -o $@ -I. -I$(QPROP_SRC)/base -lqprop -lm $(gslargs) -L$(QPROP_SRC)/../lib/x86_64 $(boostargs) -DHAVE_BOOST 

$(eval_bin) : libqprop.a eval-tsurff.cc eval-tsurff-main.cc
	$(CXX) $(optimargs) -std=c++0x eval-tsurff-main.cc eval-tsurff.cc -o $@ -I. -I$(QPROP_SRC)/base -lqprop -lm $(gslargs) -L$(QPROP_SRC)/../lib/x86_64

$(eval_mpi_bin) : libqprop.a eval-tsurff-mpi.cc eval-tsurff-mpi-main.cc
	$(MPICXX) $(optimargs) -std=c++0x eval-tsurff-mpi-main.cc eval-tsurff-mpi.cc -o $@ -I. -I$(QPROP_SRC)/base -lqprop -lm  $(gslargs) $(mpiargs) -L$(QPROP_SRC)/../lib/x86_64


install: $(bin_list)
	if [ ! -d "$(QPROP_BIN)" ]; then mkdir $(QPROP_BIN); fi
	mv $^ $(QPROP_BIN)

libqprop.a:
	+make -C $(QPROP_SRC)/base/ 

clean:
	cd $(QPROP_SRC)/base/ && make clean
	rm -f $(QPROP_SRC)/lib/x86_64/libqprop.a

clean-bin:
	if [ -d "$(QPROP_BIN)" ]; then cd $(QPROP_BIN) && rm -f $(bin_list); fi

clean-data:
	rm -f *.dat *.raw *.log *.bin

clean-all: clean clean-bin clean-data


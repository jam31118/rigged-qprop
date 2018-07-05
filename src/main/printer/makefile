QPROP_SRC = $(QPROP_HOME)/src
QPROP_BIN = $(QPROP_HOME)/bin

common_flags = -std=c++0x 
common_inc = -I$(QPROP_SRC)/base 
common_lib = -L$(QPROP_SRC)/../lib/x86_64 -lqprop

targets = print-imagpot

all: install

print-imagpot: libqprop.a print-imagpot.cc print-imagpot-main.cc print-imagpot.hh
	$(CXX) $(common_flags) -o $@ print-imagpot.cc print-imagpot-main.cc $(common_inc) $(common_lib)

install: $(targets)
	if [ ! -d "$(QPROP_BIN)" ]; then mkdir $(QPROP_BIN); fi
	mv $^ $(QPROP_BIN)

libqprop.a:
	+$(MAKE) -C $(QPROP_SRC)/base/ 


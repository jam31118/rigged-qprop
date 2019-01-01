QPROP_FLAGS = -L$(QPROP_HOME)/lib/x86_64 -I$(QPROP_HOME)/src/base -lqprop
PPP_FLAGS = -I..

MATRIX_HOME = ../../dep/matrix
MATRIX_FLAGS = -L$(MATRIX_HOME)/lib -I$(MATRIX_HOME)/include -lmatrix

OPTIMIZE_FLAGS = -O3
DEBUG_FLAGS = 

all: main

main: main.cc ../propagator.cc
	$(CXX) -o $@ $^ $(OPTIMIZE_FLAGS) $(DEBUG_FLAGS) $(QPROP_FLAGS) $(PPP_FLAGS) $(MATRIX_FLAGS)

libqprop.a: 
	$(MAKE) -C $(QPROP_HOME)/src/base

clean:
	$(RM) -f main



QPROP_FLAGS = -L$(QPROP_HOME)/lib/x86_64 -I$(QPROP_HOME)/src/base -lqprop
PPP_FLAGS = -I..

MATRIX_HOME = ../../dep/matrix
MATRIX_FLAGS = -L$(MATRIX_HOME)/lib -I$(MATRIX_HOME)/include -lmatrix

all: main

main: main.cc ../propagator.cc
	$(CXX) -o $@ $^ $(QPROP_FLAGS) $(PPP_FLAGS) $(MATRIX_FLAGS)

clean:
	$(RM) -f main



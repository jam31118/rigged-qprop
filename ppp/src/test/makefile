QPROP_FLAGS = -L$(QPROP_HOME)/lib/x86_64 -I$(QPROP_HOME)/src/base -lqprop
PPP_FLAGS = -I..

MATRIX_HOME = ../../dep/matrix
MATRIX_FLAGS = -L$(MATRIX_HOME)/lib -I$(MATRIX_HOME)/include -lmatrix

MPI_HOME = $(QPROP_DEP_DIR)/openmpi
MPI_FLAGS = -I$(MPI_HOME)/include -L$(MPI_HOME)/lib -lmpi -lmpi_cxx

OPTIMIZE_FLAGS = -O3
DEBUG_FLAGS = 

all: main

main: main.cc ../propagator.cc
	$(CXX) -o $@ $^ $(OPTIMIZE_FLAGS) $(DEBUG_FLAGS) $(QPROP_FLAGS) $(PPP_FLAGS) $(MATRIX_FLAGS) $(MPI_FLAGS)

libqprop.a: 
	$(MAKE) -C $(QPROP_HOME)/src/base

clean:
	$(RM) -f main



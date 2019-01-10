cpu_type = generic
arch_flags = -mtune=$(cpu_type) 
optimargs=-O3 $(arch_flags) -ffast-math

LIBHOME=$(QPROP_DEP_DIR)
GSLHOME=$(LIBHOME)/gsl
MPIHOME=$(LIBHOME)/openmpi
BOOSTHOME=$(LIBHOME)/boost

gslargs=-I$(GSLHOME)/include -L$(GSLHOME)/lib -lgsl -lgslcblas
mpiargs=-I$(MPIHOME)/include -L$(MPIHOME)/lib -lmpi -lmpi_cxx
boostargs=-I$(BOOSTHOME)/include -L$(BOOSTHOME)/lib

QPROP_FLAGS = -L$(QPROP_HOME)/lib/x86_64 -I$(QPROP_HOME)/src/base -lqprop

MATRIX_HOME = $(QPROP_HOME)/dep/matrix
MATRIX_FLAGS = -L$(MATRIX_HOME)/lib -I$(MATRIX_HOME)/include -lmatrix

CU_PROP_HOME = $(QPROP_HOME)/dep/cu-tridiag
CU_PROP_FLAGS = -L$(CU_PROP_HOME)/lib -I$(CU_PROP_HOME)/include -lcu_propagator



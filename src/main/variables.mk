cpu_type = generic
arch_flags = -mtune=$(cpu_type) 
optimargs=-O3 $(arch_flags) -ffast-math

LIBHOME=$(QPROP_DEP_DIR)
GSLHOME=$(LIBHOME)/gsl/
MPIHOME=$(LIBHOME)/openmpi/
BOOSTHOME=$(LIBHOME)/boost/

gslargs=-I$(GSLHOME)/include -L$(GSLHOME)/lib -lgsl -lgslcblas
mpiargs=-I$(MPIHOME)/include -L$(MPIHOME)/lib -lmpi -lmpi_cxx -DHAVE_MPI
boostargs=-I$(BOOSTHOME)/include -L$(BOOSTHOME)/lib


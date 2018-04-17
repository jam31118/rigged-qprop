optimargs=-O3 -ffast-math -march=native

LIBHOME=../prereq/tool/install
GSLHOME=$(LIBHOME)/gsl-2.4
MPIHOME=$(LIBHOME)/openmpi-1.10.7
BOOSTHOME=$(LIBHOME)/boost_1_65_1

gslargs=-I$(GSLHOME)/include -L$(GSLHOME)/lib -lgsl -lgslcblas
mpiargs=-I$(MPIHOME)/include -L$(MPIHOME)/lib -lmpi -lmpi_cxx -DHAVE_MPI
boostargs=-I$(BOOSTHOME)/include -L$(BOOSTHOME)/lib


#gslargs=-I/usr/include -L/usr/lib -lgsl -lgslcblas
#mpiargs=-I/usr/include -L/usr/lib -lmpi -lmpi_cxx -DHAVE_MPI
# for Ubuntu this seems to work for Open MPI
#mpiargs=-I/usr/include/openmpi -L/usr/lib -lmpi -lmpi_cxx -DHAVE_MPI

cpu_type = generic
arch_flags = -mtune=$(cpu_type) #-march=$(arch)
optimargs=-O3 $(arch_flags) #-ffast-math
#optimargs=-O3 -ffast-math -march=native

LIBHOME=$(QPROP_DEP_DIR)
GSLHOME=$(LIBHOME)/gsl/
MPIHOME=$(LIBHOME)/openmpi/
BOOSTHOME=$(LIBHOME)/boost/

gslargs=-I$(GSLHOME)/include -L$(GSLHOME)/lib -lgsl -lgslcblas
mpiargs=-I$(MPIHOME)/include -L$(MPIHOME)/lib -lmpi -lmpi_cxx -DHAVE_MPI
boostargs=-I$(BOOSTHOME)/include -L$(BOOSTHOME)/lib


#TOOL_HOME = /home/sjahn/tool
#GMP_HOME = $(TOOL_HOME)/gmp/gmp
#GMP_FLAGS = -I$(GMP_HOME)/include -L$(GMP_HOME)/lib
#MPC_HOME = $(TOOL_HOME)/mpc/mpc
#MPC_FLAGS = -I$(MPC_HOME)/include -L$(MPC_HOME)/lib
#GLIBC_HOME = $(TOOL_HOME)/glibc/glibc
#GLIBC_FLAGS = -I$(GLIBC_HOME)/include -L$(GLIBC_HOME)/lib
#MPFR_HOME = $(TOOL_HOME)/mpfr/mpfr
#MPFR_FLAGS = -I$(MPFR_HOME)/include -L$(MPFR_HOME)/lib
#ISL_HOME = $(TOOL_HOME)/isl/isl/
#ISL_FLAGS = -I$(ISL_HOME)/include -L$(ISL_HOME)/lib
#GCC_HOME = $(TOOL_HOME)/gcc/gcc
#GCC_FLAGS = -I$(GCC_HOME)/include -L$(GCC_HOME)/lib -L$(GCC_HOME)/lib64 -I$(GCC_HOME)/lib/gcc/x86_64-unknown-linux-gnu/5.4.0/include/ -I$(GCC_HOME)/include/c++/5.4.0 -I/home/sjahn/tool/gcc/gcc-5.4.0/include/c++/5.4.0/x86_64-unknown-linux-gnu -I/home/sjahn/tool/gcc/gcc-5.4.0/include/c++/5.4.0/x86_64-unknown-linux-gnu/32 -I/home/sjahn/tool/gcc/gcc-5.4.0/gcc-5.4.0/libstdc++-v3/include/c_compatibility


# external gcc paths
#suppress_default_flags = -nostdlib -lgcc -nostdinc -nostdinc++ -nodefaultlibs #-print-search-dirs
#extgxxopts = $(GMP_FLAGS) $(MPC_FLAGS) $(GLIBC_FLAGS) $(MPFR_FLAGS) $(ISL_FLAGS) $(GCC_FLAGS) $(suppress_default_flags) 
#extgxxopts = -L/home/sjahn/tool/gcc/gcc-5.4.0/lib/

#extgxxopts =  -L/home/sjahn/tool/mpfr/lib -L/home/sjahn/tool/mpc/mpc/lib -L/home/sjahn/tool/isl/isl/lib -L/home/sjahn/tool/glibc/glibc -L/home/sjahn/tool/gcc/gcc/lib -I/home/sjahn/tool/gcc/gcc/include 

#CXX = g++
#CXX = /home/sjahn/tool/gcc/gcc-7.3.0/bin/g++
#CXX = /home/sjahn/tool/gcc/gcc-5.4.0/bin/g++
#CXX = icpc


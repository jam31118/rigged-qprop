# Construct flags for using OpenMPI C++ binding
OMPI_INCLUDE_PATH = $(QPROP_DEP_DIR)/openmpi/include
OMPI_LIB_PATH = $(QPROP_DEP_DIR)/openmpi/lib
OMPI_FLAGS = -I$(OMPI_INCLUDE_PATH) -L$(OMPI_LIB_PATH) -lmpi -lmpi_cxx


include ../variable.mk

all: io-read

io-read: io-read.cc
	$(CXX) -o $@ $^ $(OMPI_FLAGS)

clean:
	$(RM) -f io-read


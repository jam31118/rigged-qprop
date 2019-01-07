include ../variable.mk

all: hello

hello: hello.cc
	$(CXX) -o $@ $^ $(OMPI_FLAGS)

clean:
	$(RM) -f hello


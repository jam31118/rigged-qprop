QPROP_FLAGS = -L$(QPROP_HOME)/lib/x86_64 -I$(QPROP_HOME)/src/base -lqprop
PPP_FLAGS = -I..

all: main

main: main.cc ../propagator.cc
	$(CXX) -o $@ $^ $(QPROP_FLAGS) $(PPP_FLAGS)

clean:
	$(RM) -f main



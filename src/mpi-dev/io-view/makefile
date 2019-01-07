include ../variable.mk

target = io-view

all: $(target)

$(target): $(target).cc
	$(CXX) -o $@ $^ $(OPTIMIZE_FLAGS) $(DEBUG_FLAGS) $(OMPI_FLAGS)

clean:
	$(RM) -f $(target)


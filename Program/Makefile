UNAME_S := $(shell uname -s)
CCC = g++-14
TARGETDIR=.
CCFLAGS = -fopenmp -O3 -Wall -std=c++20 -Wno-unused-but-set-variable

debug: CCFLAGS += -g 
release: CCFLAGS += -march=native
 
OBJS = \
		$(TARGETDIR)/Main/main.o 

$(TARGETDIR)/runTest: $(OBJS)
	$(CCC) $(CCFLAGS) -o $(TARGETDIR)/runTest $(OBJS)

$(TARGETDIR)/Main/main.o: $(TARGETDIR)/Main/main.cpp
	$(CCC) $(CCFLAGS) -c $< -o $@

test: runTest
	./runTest

clean:
	$(RM) \
	$(TARGETDIR)/Main/main.o \
	$(TARGETDIR)/runTest
    
debug: clean
	$(MAKE) runTest CCC=g++-14

rebuild: 
	$(MAKE) debug
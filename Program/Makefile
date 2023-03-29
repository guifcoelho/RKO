UNAME_S := $(shell uname -s)
CCC = g++
TARGETDIR=.
CCFLAGS = -O3 -Wall -std=c++11
ifeq ($(UNAME_S),Linux)
    CCFLAGS += -fopenmp
endif
ifeq ($(UNAME_S),Darwin)
    CCFLAGS += -openmp
endif

debug: CCFLAGS += -g
release: CCFLAGS += -O3 -march=native

debug: runTest
release: runTest
 
OBJS2 = \
		$(TARGETDIR)/main.o \
        $(TARGETDIR)/MultiStart.o \
		$(TARGETDIR)/SA.o \
		$(TARGETDIR)/ILS.o \
		$(TARGETDIR)/BRKGA.o \
		$(TARGETDIR)/BRKGA_QL.o \
        $(TARGETDIR)/Problem.o 

$(TARGETDIR)/runTest: $(OBJS2)
	$(CCC) $(CCFLAGS) -o $(TARGETDIR)/runTest $(OBJS2)
    
$(TARGETDIR)/Problem.o: Problem.h Problem.cpp
	$(CCC) $(CCFLAGS) -c Problem.cpp -o $(TARGETDIR)/Problem.o

$(TARGETDIR)/MultiStart.o: MultiStart.h MultiStart.cpp
	$(CCC) $(CCFLAGS) -c MultiStart.cpp -o $(TARGETDIR)/MultiStart.o

$(TARGETDIR)/SA.o: SA.h SA.cpp
	$(CCC) $(CCFLAGS) -c SA.cpp -o $(TARGETDIR)/SA.o

$(TARGETDIR)/ILS.o: ILS.h ILS.cpp
	$(CCC) $(CCFLAGS) -c ILS.cpp -o $(TARGETDIR)/ILS.o

$(TARGETDIR)/BRKGA.o: BRKGA.h BRKGA.cpp
	$(CCC) $(CCFLAGS) -c BRKGA.cpp -o $(TARGETDIR)/BRKGA.o

$(TARGETDIR)/BRKGA_QL.o: BRKGA_QL.h BRKGA_QL.cpp
	$(CCC) $(CCFLAGS) -c BRKGA_QL.cpp -o $(TARGETDIR)/BRKGA_QL.o

$(TARGETDIR)/main.o: main.cpp
	$(CCC) $(CCFLAGS) -c main.cpp -o $(TARGETDIR)/main.o


test: runTest
	./runTest

clean:
	$(RM) \
	$(TARGETDIR)/main.o \
	$(TARGETDIR)/MultiStart.o \
	$(TARGETDIR)/SA.o \
	$(TARGETDIR)/ILS.o \
	$(TARGETDIR)/BRKGA.o \
	$(TARGETDIR)/BRKGA_QL.o \
	$(TARGETDIR)/Problem.o
    
rebuild:
	$(MAKE) clean
	$(MAKE) debug
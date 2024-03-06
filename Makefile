##### Makefile for 2d ins code #####

# $(APlusPlus) is the location of the A++ library, left general here

all = ins2d

# compilers
CC = gcc
CXX = g++

opt = -O3
omp = -fopenmp

CCFLAGS= -fPIC $(opt) $(omp) -I$(APlusPlus)/include 

# List of libraries for A++
AppLibraries = -Wl,-rpath,$(APlusPlus)/lib -L$(APlusPlus)/lib -lApp -lApp_static
ompLibraries = -lgomp

# List of all libraries
LIBS = $(AppLibraries) $(ompLibraries)

# implicit (generic) rule to compile .C files
# $@ = file name of the target
# $< = name of the first prerequistite
%.o : %.C
	$(CXX) $(CCFLAGS) -o $@ -c $<

# 2D ins
ins2dFiles = ins2d.o
ins2d: $(ins2dFiles)
	$(CXX) $(CCFLAGS) -o $@ $(ins2dFiles) $(LIBS)


clean:; rm *.o




CC = g++
CPP_FLAGS = -std=c++0x -I.
L_FLAGS = -lgsl -lgslcblas

all: main debug

main: bernstein.h  bezier.h  converged.h  cubature.h  min.h  lj.h  point.h  vwrapper.h
	$(CC) $(CPP_FLAGS) -o main main.cxx hcubature.c pcubature.c $(L_FLAGS)

debug: bernstein.h  bezier.h  converged.h  cubature.h  min.h  lj.h  point.h  vwrapper.h
	$(CC) $(CPP_FLAGS) -o debug debug.cxx hcubature.c pcubature.c $(L_FLAGS)

clean:
	rm -rf main debug

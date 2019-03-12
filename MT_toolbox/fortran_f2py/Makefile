SHELL = /bin/bash

FC  = gfortran
FPY = f2py3.5

ROOT            = /Users/larnier/Documents/dias/mt/
FPYFLAGS        = -c --opt='-fcheck=all -O3 -fexternal-blas '


matrix:
	$(FPY) $(FPYFLAGS) matrixOperations.f90 -m libMatrix

all: matrix

clean: rm *.so

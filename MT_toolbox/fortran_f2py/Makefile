SHELL = /bin/bash

FC  = gfortran
FPY = f2py3.5

ROOT            = C:/Users/johnk/devProjects/Python\PyDev/MT_toolbox/
FPYFLAGS        = -c --opt='-fcheck=all -O3 -fexternal-blas '


matrix:
	$(FPY) $(FPYFLAGS) matrixOperations.f90 -m libMatrix

all: matrix

clean: rm *.so

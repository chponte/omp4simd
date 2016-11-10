# Compiler full path
CC=/home/christian.ponte/compilers/intel/bin/icc

# Intel MIC library path
CL=-L /home/christian.ponte/compilers/intel/lib/mic

# Compilation flags
CFLAGS=-mcmodel=medium -qopenmp -O2 -no-vec -mmic
#-qopt-report=4 -qopt-report-phase=all -opt-report-file=stdout

# Output directory
ODIR=$(shell pwd)/bin

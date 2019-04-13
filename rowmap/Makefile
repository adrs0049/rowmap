#
# Makefile to produce python bindings for rowmap
#

V = 0

MKDIR_P = mkdir -p
OUT_DIR = build

MODULE=_rowmap

SOURCES=$(MODULE).pyf rowmap.f90

ifndef MKLROOT
LIBS=-lblas
COMPILER=--fcompiler=gnu95
else
export CC=icc
export CFLAGS=-O3 -fPIC
export FLAGS=-O3 -fPIC
INCLUDES=-I${MKLROOT}/include
LIBS=-L${MKLROOT}/lib/intel64
LIBS+=-lmkl_rt -lmkl_intel_thread
COMPILER=--fcompiler=intelem --opt="-O3 -fopenmp"
endif

python_exec=python2

CWD := $(shell pwd)

.DELETE_ON_ERROR:

all: check-env $(MODULE)

.refresh:
	touch .refresh

check-env:
ifndef MKLROOT
	@echo "Intel-mkl not found. Compiling against system BLAS."
endif

$(MODULE): $(SOURCES)
	$(python_exec) -m numpy.f2py $(COPTS) $(COMPILER) $(INCLUDES) $(LIBS) -c $(MODULE).pyf -m $(MODULE) rowmap.f90
	#f2py $(INCLUDES) $(LIBS) -c $(MODULE).pyf -m $(MODULE) rowmap.f90

.PHONY: setup clean force all

clean:
	rm -f *.so

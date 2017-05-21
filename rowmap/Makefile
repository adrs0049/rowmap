#
# Makefile to produce python bindings for rowmap
#

V = 0

MKDIR_P = mkdir -p
OUT_DIR = build

MODULE=_rowmap

SOURCES=$(MODULE).pyf rowmap.f90
INCLUDES=
LIBS=-L/opt/intel/compilers_and_libraries_2016.3.210/linux/mkl/lib/intel64 -lmkl_rt


CWD := $(shell pwd)

.DELETE_ON_ERROR:

all: $(MODULE)

.refresh:
	touch .refresh

$(MODULE): $(SOURCES)
	f2py $(INCLUDES) $(LIBS) -c $(MODULE).pyf -m $(MODULE) rowmap.f90

.PHONY: setup clean force all

clean:
	rm $(MODULE).so

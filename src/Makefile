# ARCH should be defined by caller
ARCH ?=
arch ?= x86_64
COPT ?= -O2
GCC ?= gcc
EXEEXT ?=
CFLAGS := -D_FILE_OFFSET_BITS=64 -Wall $(COPT)

all: ../bin$(ARCH)/polyt$(EXEEXT)

clean:
	rm -f ../bin$(ARCH)/*

../bin$(ARCH)/polyt$(EXEEXT): polyt.c $(TOOLS)
	$(GCC) polyt.c $(CFLAGS) -o ../bin$(ARCH)/polyt$(EXEEXT)

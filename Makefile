#
# Makefile to build stencil MPI
#

CC=mpicc

COMP=GNU
ifeq ($(COMP), GNU)
  CFLAGS=-std=c99 -qopenmp -pthread -Wall -Ofast
endif

EXE1=stencil.exe
EXES=$(EXE1) $(EXE2)

all: $(EXES)

$(EXES): %.exe : %.c
	$(CC) $(CFLAGS) -o $@ $^

.PHONY: clean all

clean:
	\rm -f $(EXES)
	\rm -f *.o


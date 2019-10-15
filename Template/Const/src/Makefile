#
# Makefile for test
#

.SUFFIXES:
.SUFFIXES: .f90 .o

.PHONY: clean help

# FC = f90
# FC = g95
FC = gfortran
FC = g95 -mmacosx-version-min=10.13

# FFLAGS = -Wall
FFLAGS =


.f90.o:
	$(FC) $(FFLAGS) -c $<


OBJS = 	const.o \
	main.o

all: runme
	./runme

runme: main.o
	$(FC) $(FFLAGS) -o runme $(OBJS)

main.o: const.o

help:
	echo 'Usage: make'

clean:
	rm -f runme core $(OBJS) *.mod

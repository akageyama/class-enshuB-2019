#
# Makefile for test
#

.SUFFIXES:
.SUFFIXES: .f90 .o

.PHONY: clean help

# FC = f90
# FC = g95
FC = gfortran

# FFLAGS = -Wall
FFLAGS =


.f90.o:
	$(FC) $(FFLAGS) -c $<


OBJS = 	const.o \
	main.o

runme: main.o
	$(FC) $(FFLAGS) -o runme $(OBJS)

main.o: const.o

help:
	echo 'Usage: make && ./runme'

clean:
	rm -f runme core $(OBJS) *.mod

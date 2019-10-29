#
# Makefile for test
#

.SUFFIXES:
.SUFFIXES: .f90 .o

.PHONY: clean help

# FC = f90
FC = g95
# FC = gfortran

# FFLAGS = -Wall
FFLAGS =


.f90.o:
	$(FC) $(FFLAGS) -c $<


OBJS = const.o \
       pbm.o \
       ppm.o \
       main.o

temp_grayscale.ppm temp_kobe_rev.pbm: runme
	./runme

runme: main.o
	$(FC) $(FFLAGS) -o runme $(OBJS)

main.o: const.o ppm.o pbm.o

help:
	echo 'Usage: make'

clean:
	rm -f runme core $(OBJS) *.mod
	rm -rf temp_grayscale.ppm temp_kobe_rev.pbm

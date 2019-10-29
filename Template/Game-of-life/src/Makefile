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
       gol.o \
       main.o

movie.gif: runme kobe_bitmap.pbm pentomino_bitmap.pbm start.pbm
	./runme
	convert -delay 5 ../data/*.pbm movie.gif

runme: main.o
	$(FC) $(FFLAGS) -o runme $(OBJS)

main.o: const.o pbm.o gol.o

gol.o: const.o pbm.o

help:
	echo 'Usage: make && ./runme'

clean:
	rm -f runme core $(OBJS) *.mod
	rm -f ../data/*.pbm movie.gif

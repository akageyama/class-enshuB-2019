#
# Makefile for SmoothLife
#

.SUFFIXES:
.SUFFIXES: .f90 .o

.PHONY: clean help

# FC = f90
# FC = gfortran
FC = g95 -Wall -mmacosx-version-min=10.13

# FFLAGS = -Wall
FFLAGS = -O5 -fopenmp 
# FFLAGS = -Wall -fbounds-check  -O -Wuninitialized


.f90.o:
	$(FC) $(FFLAGS) -c $<


# OBJS = const.o \
#        pbm.o \
#        gol.o \
#        main.o

OBJS = const.o \
       pgm.o \
       sml.o \
       main.o

# movie.gif: runme kobe_gray.pgm kobe_bitmap.pbm pentomino_bitmap.pbm start.pbm
# 	./runme
# 	convert -delay 5 ../data/*.pbm movie.gif

movie.gif: runme sample_face.pgm
	./runme
	convert -delay 5 ../data/*.pgm movie.gif

runme: main.o
	$(FC) $(FFLAGS) -o runme $(OBJS)

# main.o: const.o pbm.o gol.o
main.o: const.o pgm.o sml.o

# gol.o: const.o pbm.o

sml.o: const.o pgm.o

help:
	echo 'Usage: make && ./runme'

clean:
	rm -f runme core $(OBJS) *.mod
	rm -f ../data/*.pgm ../data/*.pbm movie.gif

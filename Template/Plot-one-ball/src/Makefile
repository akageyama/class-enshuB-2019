#
# Makefile for test
#

.SUFFIXES:
.SUFFIXES: .f90 .o

.PHONY: clean help

# FC = f90
# FC = gfortran
FC = g95 -mmacosx-version-min=10.13

# FFLAGS = -Wall
FFLAGS =


.f90.o:
	$(FC) $(FFLAGS) -c $<


OBJS = 	const.o  \
        particle.o  \
        vect2d.o  \
	main.o

fig: particle.pos.data makefig.plt
	gnuplot makefig.plt

runme: main.o
	$(FC) $(FFLAGS) -o runme $(OBJS)
        
particle.pos.data: runme
	./runme  # will save to 'particle.pos.data'

main.o: const.o  \
        particle.o  \
        vect2d.o

particle.o: const.o  \
            vect2d.o

help:
	@echo 'Usage: make fig'

clean:
	rm -f runme core $(OBJS) *.mod
	rm -f particle.pos.data

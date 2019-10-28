#
# Makefile for test
#
#  make or make FC=gfortran
#

.SUFFIXES:
.SUFFIXES: .f90 .o

.PHONY: clean help

# THISMACHINE = $(shell hostname | cut -d. -f3)
# ENSHUSHITSU = center

# CENT = a

# FC = f90
FC = g95 -mmacosx-version-min=10.13


.f90.o:
	$(FC) $(FFLAGS) -c $<


OBJS = 	const.o  \
        particles.o  \
        vect2d.o  \
	main.o

fig: particles.pos.data.00000 makefig.plt
	gnuplot makefig.plt

runme: main.o
	$(FC) $(FFLAGS) -o runme $(OBJS)
        
particles.pos.data.00000: runme
	./runme  # will save to 'particles.pos.data.?????'

main.o: const.o  \
        particles.o  \
        vect2d.o

particles.o: const.o  \
             vect2d.o

help:
	@echo 'Usage: Just type "make"'
	@echo ' Note: Default compiler is "g95 -mmacosx-version-min=10.13"'
	@echo '       If you have gfortran, type "make FC=gfortran"'

clean:
	rm -f runme core $(OBJS) *.mod
	rm -f particles.pos.data

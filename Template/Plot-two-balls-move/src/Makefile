#
# Makefile for test
#

.SUFFIXES:
.SUFFIXES: .f90 .o

.PHONY: clean help

# THISMACHINE = $(shell hostname | cut -d. -f3)
# ENSHUSHITSU = center

# CENT = a

# FC = f90
FC = gfortran
# ifeq ($(THISMACHINE),$(ENSHUSHITSU))
#   FC = g95 -mmacosx-version-min=10.13
# endif
ifdef C
  FC = g95 -mmacosx-version-min=10.13
endif

test:
	echo $(FC)

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
	@echo 'Usage: make fig'

clean:
	rm -f runme core $(OBJS) *.mod
	rm -f particles.pos.data.?????

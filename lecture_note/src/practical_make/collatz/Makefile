.PHONY: all clean

all: collatz
	./collatz

collatz: collatz.o graph.o
	gcc -o collatz collatz.o graph.o

collatz.o: collatz.c collatz.h
	gcc -c collatz.c

graph.o: graph.c
	gcc -c graph.c

clean:
	rm -rf collatz.o graph.o collatz


.PHONY: clean all

all: collatz01 collatz02 collatz03 collatz04 collatz05


%: %.c collatz.h
	gcc $< -o $@


clean:
	rm -f collatz0? collatz0?.o *~


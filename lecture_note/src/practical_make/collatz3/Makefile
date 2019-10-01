.PHONY: max count clean


max: out.txt
	./extract_top_max.sh out.txt

count: out.txt
	./extract_top_count.sh out.txt

out.txt: collatz3
	./collatz3 > out.txt

collatz3: collatz3.c
	gcc -o collatz3 collatz3.c

clean:
	rm -rf out.txt collatz3.o collatz3


#include <stdio.h>

typedef struct {
  long initial_n;
  long repeat_count;
  long max_n;
} Record;

void collatz(long n) {
  Record rec;

  rec.initial_n = n;
  rec.repeat_count = 0;
  rec.max_n = n;

  while ( n > 1 ) {
    rec.repeat_count++;
    if ( n % 2 == 0 ) {
      n = n / 2;
    }
    else {
      n = 3 * n + 1;
      if ( n > rec.max_n ) rec.max_n = n;
    }
    // printf("n = %d\n", n);
  }
  printf("%ld %ld %ld\n", rec.max_n, rec.repeat_count, rec.initial_n);
}

int main(void) {
  for (long i=1; i<=10000; i++) {
    collatz(i);
  }
}

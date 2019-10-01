#include <stdio.h>
#include <unistd.h>
#include "collatz.h"

int main(void) {
  int n = N;
  printf(" n = %d\n", n);
  while ( n > 1 ) {
    if ( n % 2 == 0 ) {
      n = n / 2;
    }
    else {
      n = 3 * n + 1;
    }
    printf(" n = %d\n", n);
    sleep(3);
  }
}

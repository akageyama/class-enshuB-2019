#include <stdio.h>
#include <unistd.h>

#include "collatz.h"

extern void graph();

int main(void) {
  int n = N_INITIAL;
  //printf(" n = %d\n", n);
  graph(n);

  while ( n > 1 ) {
    if ( n % 2 == 0 ) {
      n = n / 2;
    }
    else {
      n = 3 * n + 1;
    }
    //printf(" n = %d\n", n);
    graph(n);
    // sleep(1);
  }
}

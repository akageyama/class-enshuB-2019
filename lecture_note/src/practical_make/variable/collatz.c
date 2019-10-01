#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

int main(int argc,char *argv[]) {
  int n = 100; // default
  if (argc>1) n = atoi(argv[1]);
  printf(" n = %d\n", n);

  while ( n > 1 ) {
    if ( n % 2 == 0 ) {
      n = n / 2;
    }
    else {
      n = 3 * n + 1;
    }
    printf(" n = %d\n", n);
    sleep(1);
  }
}

#include <stddef.h>
#include <stdio.h>
#include <math.h>

#include "crandom.h"

#define T (10000000)
#define N (1000)
#define A (-5)
#define B (5)


size_t histogram[N + 1];


int main() {
  size_t i = 0;
  struct cRandom * crandom = dSFMTRandomNew();

  for(i = 0; i < T; ++i) {
    double x = normal(crandom, 0, 1);

    if( x < A )
      x = A;
    if( B < x )
      x = B;

    histogram[ (size_t) floor(N * (x - A) / (B - A)) ]++;
  }

  crandom->release(crandom);

  for(i = 0; i <= N; ++i) {
    if( histogram[i] )
      printf("%.3f %.3f\n", A + ((double)i) * (B - A) / N, ((double)histogram[i]) * N / (T * (B - A)));
  }

  return 0;
}

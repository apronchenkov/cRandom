/** 
 * Copyright (C) 2009 Alexander G. Pronchenkov. All rights reserved.
 *
 * The new BSD License is applied to this software.
 */

/* This code is based on the library rvgs.c (Steve Park & Dave Geyer) */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "crandom.h"
#include "dSFMT/dSFMT.h"


/**
 * Interface for random variates generator form uniformly distribute
 */
struct dSFMTRandom {
  struct cRandom crandom;

  dsfmt_t dsfmt;
};


/**
 * Returns the next pseudorandom, uniformly distributed double value
 * between 0.0 and 1.0 from this random number generator's sequence.
 *
 * Range: 0 <= x < 1
 */
double dSFMTRandomNext(void * that) {
  struct dSFMTRandom * random = (struct dSFMTRandom *) that;

  return dsfmt_genrand_close_open(&random->dsfmt);
}


/**
 * Create a new cRandom object (dSFMT based)
 *
 * Initialize it by 32-bit integer.
 */
struct cRandom * dSFMTRandomNewBySeed(int seed) {
  struct dSFMTRandom * random = (struct dSFMTRandom *) malloc(sizeof(*random));

  if( random == NULL )
    return NULL;

  dsfmt_init_gen_rand(&random->dsfmt, seed);

  random->crandom.next = &dSFMTRandomNext;
  random->crandom.release = &free;

  return (struct cRandom *)random;
}


/**
 * Create a new cRandom object (dSFMT based)
 *
 * Initialize it by array.
 */
struct cRandom * dSFMTRandomNewByArray(int * array, int arrayLength) {
  struct dSFMTRandom * random = (struct dSFMTRandom *) malloc(sizeof(*random));

  if( random == NULL )
    return NULL;

  dsfmt_init_by_array(&random->dsfmt, (uint32_t*)array, (arrayLength * sizeof(int)) / sizeof(uint32_t));

  random->crandom.next = &dSFMTRandomNext;
  random->crandom.release = &free;

  return (struct cRandom *)random;
}


/**
 * Create a new cRandom object (dSFMT based)
 */
struct cRandom * dSFMTRandomNew(void) {
  return dSFMTRandomNewBySeed( (int) time(NULL));
}



/***********************
 * With finite support *
 ***********************/


/**
 * Returns 1 with probability p or 0 with probability 1 - p. 
 * NOTE: use 0.0 < p < 1.0
 *
 * Range:    0, 1
 * Mean:     p
 * Variance: p * (1 - p)
 */
int bernoulli(struct cRandom * crandom, double p) {
  assert( 0.0 < p && p < 1.0 );

  return (crandom->next(crandom) > p);
}


/**
 * Returns a binomial distributed integer between 0 and n inclusive. 
 * NOTE: use n > 0 and 0.0 < p < 1.0
 *
 * Range:    0, ..., n
 * Mean:     n * p
 * Variance: n * p * (1 - p)
 */
int binomial(struct cRandom * crandom, int n, double p) {
  int i, x = 0;

  assert( 0 < n );
  assert( 0.0 < p && p < 1.0 );

  for(i = 0; i < n; ++i)
    x += (crandom->next(crandom) > p); /* x += Bernoulli(crandom, p); */

  return x;
}


/**
 * Returns an discrete uniform distributed integer between a and b inclusive. 
 * NOTE: use a < b
 *         
 * Range:    a, ..., b
 * Mean:     (a + b) / 2
 * Variance: (sqr(b - a + 1) - 1) / 12
 */
int equilikely(struct cRandom * crandom, int a, int b) {
  assert( a < b );

  return (a + (int) ((b - a + 1) * crandom->next(crandom)));
}


/**
 * Returns a geometric distributed non-negative integer.
 * NOTE: use 0.0 < p < 1.0
 *
 * Range:    0, ...
 * Mean:     p / (1 - p)
 * Variance: p / sqr(1 - p)
 */
int geometric(struct cRandom * crandom, double p) {
  assert( 0.0 < p && p < 1.0 );

  return ((int) (log(1.0 - crandom->next(crandom)) / log(p)));
}


/**
 * Returns a Pascal distributed non-negative integer. 
 * NOTE: use n > 0 and 0.0 < p < 1.0
 *
 * Range:    0, ...
 * Mean:     n * p / (1 - p)
 * Variance: n * p / sqr(1 - p)
 */
int pascal(struct cRandom * crandom, int n, double p) {
  const double log_p = log(p);
  int i, x = 0;

  assert( 0 < n );
  assert( 0.0 < p && p < 1.0 );

  for (i = 0; i < n; i++)
    x += ((int) (log(1.0 - crandom->next(crandom)) / log_p)); /* x += geometric(crandom, p); */

  return (x);
}


/**
 * Returns a Poisson distributed non-negative integer. 
 * NOTE: use m > 0.0
 *
 * Range:    0, ...
 * Mean:     m
 * Variance: m
 */
int Poisson(struct cRandom * crandom, double m) {
  double t = 0.0;
  int x = -1;

  assert( 0.0 < m );

  while (t < m) {
    t -= m * log(1.0 - crandom->next(crandom)); /* t += exponential(crandom, 1.0); */
    x++;
  }

  return x;
}


/*************************
 * With infinite support *
 *************************/


/**
 * Returns a uniformly distributed real number between a and b. 
 * NOTE: use a < b
 *
 * Range:    a < x < b
 * Mean:     (a + b) / 2
 * Variance: sqr(b - a) / 12 
 */
double uniform(struct cRandom * crandom, double a, double b) {
  assert( a < b );

  return a + (b - a) * crandom->next(crandom);
}


/**
 * Returns an exponentially distributed positive real number. 
 * NOTE: use m > 0.0
 *
 * Range:    0 < x
 * Mean:     m
 * Variance: sqr(m)
 */
double exponential(struct cRandom * crandom, double m) {
  assert( 0.0 < m );

  return - m * log(1.0 - crandom->next(crandom));
}


/**
 * Returns an Erlang distributed positive real number.
 * NOTE: use n > 0 and b > 0.0
 *
 * Range:    0 < x
 * Mean:     n * b
 * Variance: n * sqr(b)
 */
double erlang(struct cRandom * crandom, int n, double b) {
  int i;
  double x = 0.0;

  assert( 0 < n );
  assert( 0.0 < b );

  for (i = 0; i < n; i++) 
    x -= b * log(1.0 - crandom->next(crandom)); /* x += exponential(crandom, b); */

  return (x);

}


/**
 * Returns a normal (Gaussian) distributed real number.
 * NOTE: use s > 0.0
 *
 * Range:    all x
 * Mean:     m
 * Variance: sqr(s)
 */
double normal(struct cRandom * crandom, double m, double s) {
  /*
   * Uses a very accurate approximation of the normal idf due to Odeh & Evans, 
   * J. Applied Statistics, 1974, vol 23, pp 96-97.
   */
  const double p0 = 0.322232431088;     const double q0 = 0.099348462606;
  const double p1 = 1.0;                const double q1 = 0.588581570495;
  const double p2 = 0.342242088547;     const double q2 = 0.531103462366;
  const double p3 = 0.204231210245e-1;  const double q3 = 0.103537752850;
  const double p4 = 0.453642210148e-4;  const double q4 = 0.385607006340e-2;
  double u, t, p, q, z;

  assert( 0.0 < s );

  u = crandom->next(crandom);
  if( u < 0.5 )
    t = sqrt(-2.0 * log(u));
  else
    t = sqrt(-2.0 * log(1.0 - u));
  p   = p0 + t * (p1 + t * (p2 + t * (p3 + t * p4)));
  q   = q0 + t * (q1 + t * (q2 + t * (q3 + t * q4)));
  if( u < 0.5 )
    z = (p / q) - t;
  else
    z = t - (p / q);

  return (m + s * z);
}


/**
 * Returns a lognormal distributed positive real number. 
 * NOTE: use b > 0.0
 *
 * Range:    0 < x
 * Mean:     exp(a + sqr(b) / 2)
 * Variance: (exp(sqr(b) - 1) * exp(2 * a + sqr(b))
 */
double lognormal(struct cRandom * crandom, double a, double b) {
  assert( 0.0 < b );

  return (exp(a + b * normal(crandom, 0.0, 1.0)));
}


/**
 * Returns a chi-square distributed positive real number. 
 * NOTE: use n > 0
 *
 * Range:    0 < x
 * Mean:     n
 * Variance: 2 * n
 */
double chisquare(struct cRandom * crandom, int n) {
  long   i;
  double z, x = 0.0;

  assert( 0 < n );

  for (i = 0; i < n; ++i) {
    z = normal(crandom, 0.0, 1.0);
    x += z * z;
  }

  return x;

}


/**
 * Returns a student-t distributed real number.
 * NOTE: use n > 0
 *
 * Range:    all x
 * Mean:     0           (when n > 1)
 * Variance: n / (n - 2) (when n > 2)
 */
double student(struct cRandom * crandom, int n) {
  assert( 0 < n );

  return (normal(crandom, 0.0, 1.0) / sqrt(chisquare(crandom, n) / n));
}

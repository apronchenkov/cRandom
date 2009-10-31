/** 
 * Copyright (C) 2009 Alexander G. Pronchenkov. All rights reserved.
 *
 * The new BSD License is applied to this software.
 */

/* This code is based on the library rvgs.c (Steve Park & Dave Geyer) */

#ifndef __crandom_h__
#define __crandom_h__

#ifdef __cplusplus
extern "C" {
#endif


/**
 * Interface for random variates generator form uniformly distribute
 */
struct cRandom {
  /**
   * Returns the next pseudorandom, uniformly distributed double value
   * between 0.0 and 1.0 from this random number generator's sequence.
   *
   * Range: 0 <= x < 1
   */
  double (* next)(void * that);


  /**
   * Releases resources of a cRandom object
   */
  void (* release)(void * that);
};


/**
 * Create a new cRandom object (dSFMT based)
 */
struct cRandom * dSFMTRandomNew(void);


/**
 * Create a new cRandom object (dSFMT based)
 *
 * Initialize it by 32-bit integer.
 */
struct cRandom * dSFMTRandomNewBySeed(int seed);


/**
 * Create a new cRandom object (dSFMT based)
 *
 * Initialize it by array.
 */
struct cRandom * dSFMTRandomNewByArray(int * array, int arrayLength);




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
int bernoulli(struct cRandom * crandom, double p);


/**
 * Returns a binomial distributed integer between 0 and n inclusive. 
 * NOTE: use n > 0 and 0.0 < p < 1.0
 *
 * Range:    0, ..., n
 * Mean:     n * p
 * Variance: n * p * (1 - p)
 */
int binomial(struct cRandom * crandom, int n, double p);


/**
 * Returns an discrete uniform distributed integer between a and b inclusive. 
 * NOTE: use a < b
 *         
 * Range:    a, ..., b
 * Mean:     (a + b) / 2
 * Variance: (sqr(b - a + 1) - 1) / 12
 */
int equilikely(struct cRandom * crandom, int a, int b);


/**
 * Returns a geometric distributed non-negative integer.
 * NOTE: use 0.0 < p < 1.0
 *
 * Range:    0, ...
 * Mean:     p / (1 - p)
 * Variance: p / sqr(1 - p)
 */
int geometric(struct cRandom * crandom, double p);


/**
 * Returns a Pascal distributed non-negative integer. 
 * NOTE: use n > 0 and 0.0 < p < 1.0
 *
 * Range:    0, ...
 * Mean:     n * p / (1 - p)
 * Variance: n * p / sqr(1 - p)
 */
int pascal(struct cRandom * crandom, int n, double p);


/**
 * Returns a Poisson distributed non-negative integer. 
 * NOTE: use m > 0
 *
 * Range:    0, ...
 * Mean:     m
 * Variance: m
 */
int Poisson(struct cRandom * crandom, double m);


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
double uniform(struct cRandom * crandom, double a, double b);


/**
 * Returns an exponentially distributed positive real number. 
 * NOTE: use m > 0.0
 *
 * Range:    0 < x
 * Mean:     m
 * Variance: sqr(m)
 */
double exponential(struct cRandom * crandom, double m);


/**
 * Returns an Erlang distributed positive real number.
 * NOTE: use n > 0 and b > 0.0
 *
 * Range:    0 < x
 * Mean:     n * b
 * Variance: n * sqr(b)
 */
double erlang(struct cRandom * crandom, int n, double b);


/**
 * Returns a normal (Gaussian) distributed real number.
 * NOTE: use s > 0.0
 *
 * Range:    all x
 * Mean:     m
 * Variance: sqr(s)
 */
double normal(struct cRandom * crandom, double m, double s);


/**
 * Returns a lognormal distributed positive real number. 
 * NOTE: use b > 0.0
 *
 * Range:    0 < x
 * Mean:     exp(a + sqr(b) / 2)
 * Variance: (exp(sqr(b) - 1) * exp(2 * a + sqr(b))
 */
double lognormal(struct cRandom * crandom, double a, double b);


/**
 * Returns a chi-square distributed positive real number. 
 * NOTE: use n > 0
 *
 * Range:    0 < x
 * Mean:     n
 * Variance: 2 * n
 */
double chisquare(struct cRandom * crandom, int n);


/**
 * Returns a student-t distributed real number.
 * NOTE: use n > 0
 *
 * Range:    all x
 * Mean:     0           (when n > 1)
 * Variance: n / (n - 2) (when n > 2)
 */
double student(struct cRandom * crandom, int n);


#ifdef __cplusplus
}
#endif


#endif /*_crandom_h__*/

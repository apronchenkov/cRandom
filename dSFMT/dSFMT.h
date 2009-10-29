/** 
 * @file dSFMT.h 
 *
 * @brief double precision SIMD oriented Fast Mersenne Twister(dSFMT)
 * pseudorandom number generator based on IEEE 754 format.
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (C) 2007, 2008 Mutsuo Saito, Makoto Matsumoto and
 * Hiroshima University. All rights reserved.
 *
 * The new BSD License is applied to this software.
 * see LICENSE.txt
 *
 * @note We assume that your system has inttypes.h.
 */

#ifndef DSFMT_H
#define DSFMT_H

#include <stdio.h>
#include <assert.h>

#if !defined(DSFMT_MEXP)
#ifdef __GNUC__
  #warning "DSFMT_MEXP is not defined. I assume DSFMT_MEXP is 19937."
#endif
  #define DSFMT_MEXP 19937
#endif
/*-----------------
  BASIC DEFINITIONS
  -----------------*/
/* Mersenne Exponent. The period of the sequence 
 *  is a multiple of 2^DSFMT_MEXP-1. */
/** DSFMT generator has an internal state array of 128-bit integers,
 * and N is its size. */
#define DSFMT_N ((DSFMT_MEXP - 128) / 104 + 1)
/** N64 is the size of internal state array when regarded as an array
 * of 64-bit integers.*/
#define DSFMT_N64 (DSFMT_N * 2)

#if !defined(DSFMT_BIG_ENDIAN)
#  if defined(__BYTE_ORDER) && defined(__BIG_ENDIAN)
#    if __BYTE_ORDER == __BIG_ENDIAN
#      define DSFMT_BIG_ENDIAN 1
#    endif
#  elif defined(_BYTE_ORDER) && defined(_BIG_ENDIAN)
#    if _BYTE_ORDER == _BIG_ENDIAN
#      define DSFMT_BIG_ENDIAN 1
#    endif
#  elif defined(__BYTE_ORDER__) && defined(__BIG_ENDIAN__)
#    if __BYTE_ORDER__ == __BIG_ENDIAN__
#      define DSFMT_BIG_ENDIAN 1
#    endif
#  elif defined(BYTE_ORDER) && defined(BIG_ENDIAN)
#    if BYTE_ORDER == BIG_ENDIAN
#      define DSFMT_BIG_ENDIAN 1
#    endif
#  elif defined(__BIG_ENDIAN) || defined(_BIG_ENDIAN) \
    || defined(__BIG_ENDIAN__) || defined(BIG_ENDIAN)
#      define DSFMT_BIG_ENDIAN 1
#  endif
#endif

#if defined(DSFMT_BIG_ENDIAN) && defined(__amd64)
#  undef DSFMT_BIG_ENDIAN
#endif

#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)
#  include <inttypes.h>
#elif defined(_MSC_VER) || defined(__BORLANDC__)
#  if !defined(DSFMT_UINT32_DEFINED) && !defined(SFMT_UINT32_DEFINED)
typedef unsigned int uint32_t;
typedef unsigned __int64 uint64_t;
#    define UINT64_C(v) (v ## ui64)
#    define DSFMT_UINT32_DEFINED
#    if !defined(inline)
#      define inline __inline
#    endif
#  endif
#else
#  include <inttypes.h>
#  if !defined(inline)
#    if defined(__GNUC__)
#      define inline __inline__
#    else
#      define inline
#    endif
#  endif
#endif

#ifndef UINT64_C
#  define UINT64_C(v) (v ## ULL) 
#endif

/*------------------------------------------
  128-bit SIMD like data type for standard C
  ------------------------------------------*/
#if defined(HAVE_ALTIVEC)
#  if !defined(__APPLE__)
#    include <altivec.h>
#  endif
/** 128-bit data structure */
union W128_T {
    vector unsigned int s;
    uint64_t u[2];
    uint32_t u32[4];
    double d[2];
};

#elif defined(HAVE_SSE2)
#  include <emmintrin.h>

/** 128-bit data structure */
union W128_T {
    __m128i si;
    __m128d sd;
    uint64_t u[2];
    uint32_t u32[4];
    double d[2];
};
#else  /* standard C */
/** 128-bit data structure */
union W128_T {
    uint64_t u[2];
    uint32_t u32[4];
    double d[2];
};
#endif

/** 128-bit data type */
typedef union W128_T w128_t;

/** the 128-bit internal state array */
struct DSFMT_T {
    w128_t status[DSFMT_N + 1];
    int idx;
};
typedef struct DSFMT_T dsfmt_t;


void dsfmt_gen_rand_all(dsfmt_t *dsfmt);
void dsfmt_fill_array_open_close(dsfmt_t *dsfmt, double array[], int size);
void dsfmt_fill_array_close_open(dsfmt_t *dsfmt, double array[], int size);
void dsfmt_fill_array_open_open(dsfmt_t *dsfmt, double array[], int size);
void dsfmt_fill_array_close1_open2(dsfmt_t *dsfmt, double array[], int size);
void dsfmt_init_gen_rand(dsfmt_t *dsfmt, uint32_t seed);
void dsfmt_init_by_array(dsfmt_t *dsfmt, uint32_t init_key[],
			     int key_length);
const char *dsfmt_get_idstring(void);
int dsfmt_get_min_array_size(void);

#if defined(__GNUC__)
#  define DSFMT_PRE_INLINE inline static
#  define DSFMT_PST_INLINE __attribute__((always_inline))
#elif defined(_MSC_VER) && _MSC_VER >= 1200
#  define DSFMT_PRE_INLINE __forceinline static
#  define DSFMT_PST_INLINE
#else
#  define DSFMT_PRE_INLINE inline static
#  define DSFMT_PST_INLINE
#endif
DSFMT_PRE_INLINE double dsfmt_genrand_close1_open2(dsfmt_t *dsfmt)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double dsfmt_genrand_close_open(dsfmt_t *dsfmt)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double dsfmt_genrand_open_close(dsfmt_t *dsfmt)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double dsfmt_genrand_open_open(dsfmt_t *dsfmt)
    DSFMT_PST_INLINE;


/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range [1, 2).  This is
 * the primitive and faster than generating numbers in other ranges.
 * dsfmt_init_gen_rand() or dsfmt_init_by_array() must be called
 * before this function.
 * @param dsfmt dsfmt internal state date
 * @return double precision floating point pseudorandom number
 */
inline static double dsfmt_genrand_close1_open2(dsfmt_t *dsfmt) {
    double r;
    double *psfmt64 = &dsfmt->status[0].d[0];

    if (dsfmt->idx >= DSFMT_N64) {
	dsfmt_gen_rand_all(dsfmt);
	dsfmt->idx = 0;
    }
    r = psfmt64[dsfmt->idx++];
    return r;
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range [0, 1).
 * dsfmt_init_gen_rand() or dsfmt_init_by_array() must be called
 * before this function.
 * @param dsfmt dsfmt internal state date
 * @return double precision floating point pseudorandom number
 */
inline static double dsfmt_genrand_close_open(dsfmt_t *dsfmt) {
    return dsfmt_genrand_close1_open2(dsfmt) - 1.0;
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range (0, 1].
 * dsfmt_init_gen_rand() or dsfmt_init_by_array() must be called
 * before this function. 
 * @param dsfmt dsfmt internal state date
 * @return double precision floating point pseudorandom number
 */
inline static double dsfmt_genrand_open_close(dsfmt_t *dsfmt) {
    return 2.0 - dsfmt_genrand_close1_open2(dsfmt);
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range (0, 1).
 * dsfmt_init_gen_rand() or dsfmt_init_by_array() must be called
 * before this function.
 * @param dsfmt dsfmt internal state date
 * @return double precision floating point pseudorandom number
 */
inline static double dsfmt_genrand_open_open(dsfmt_t *dsfmt) {
    double *dsfmt64 = &dsfmt->status[0].d[0];
    union {
	double d;
	uint64_t u;
    } r;

    if (dsfmt->idx >= DSFMT_N64) {
	dsfmt_gen_rand_all(dsfmt);
	dsfmt->idx = 0;
    }
    r.d = dsfmt64[dsfmt->idx++];
    r.u |= 1;
    return r.d - 1.0;
}

#endif /* DSFMT_H */

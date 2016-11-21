/* 128 bit UUID computation / Rabin fingerprinting algorithm                     */
/*                                                                               */
/* F. Prill, DWD                                                                 */
/*                                                                               */
/* @par Copyright and License                                                    */
/*                                                                               */
/* This code is subject to the DWD and MPI-M-Software-License-Agreement in       */
/* its most recent form.                                                         */
/* Please see the file LICENSE in the root of the source tree for this code.     */
/* Where software is supplied by third parties, it is indicated in the           */
/* headers of the routines.                                                      */

#ifndef UTIL_UUID_H
#define UTIL_UUID_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

typedef struct
{
  unsigned char data[16];
} uuid_t;

typedef enum {
  UUID_EQUAL,
  UUID_EQUAL_LIMITED_ACCURACY,
  UUID_UNEQUAL
} cmp_UUID_t;


typedef struct 
{
  unsigned int initialized;

  unsigned int f48[2]; // 48 bit fingerprint
  unsigned int f64[2]; // 64 bit fingerprint

  // pre-computed tables for fingerprint algorithm
  unsigned int 
    TA64[256][2], TB64[256][2], TC64[256][2], TD64[256][2],
    TC48[256][2], TD48[256][2], TE48[256][2], TF48[256][2];

} context_t;


void            uuid_generate(const double* val, const int nval, uuid_t* uuid);
cmp_UUID_t      compare_UUID(const uuid_t uuid_A, const uuid_t uuid_B, 
                             double* min_difference);
void            uuid_unparse(char *buffer, const uuid_t *uuid);

/* private routines */

void            longRightShift(unsigned int* sWord, unsigned int s);
void            longLeftShift(unsigned int* sWord, unsigned int s);
void            generate_polynomial(unsigned int k, unsigned int* W);
void            multiT8_64(unsigned int* s, unsigned int* in_s, const unsigned int* p);
void            multiT8_48(unsigned int* s, unsigned int* in_s, const unsigned int* p);
void            compute_table_64(const unsigned int*  p, unsigned int TA[256][2], unsigned int TB[256][2], 
                                 unsigned int TC[256][2], unsigned int TD[256][2]);
void            compute_table_48(const unsigned int*  p, 
                                 unsigned int TC[256][2], unsigned int TD[256][2], unsigned int TE[256][2], 
                                 unsigned int TF[256][2]);
void            fp4_64(unsigned int* s, const unsigned char* w, unsigned int TA[256][2], unsigned int TB[256][2], 
                       unsigned int TC[256][2], unsigned int TD[256][2]);
void            fp4_48(unsigned int* s, const unsigned char* w, 
                       unsigned int TC[256][2], unsigned int TD[256][2], unsigned int TE[256][2], 
                       unsigned int TF[256][2]);
void            fp1_64(unsigned int* s, unsigned char w, const unsigned int* p);
void            fp1_48(unsigned int* s, unsigned char w, const unsigned int* p);
int             extract_exponent64(const double val);
int             extract_exponent32(const float val);
void            convert_to_fixed64(const double val, unsigned long long* fp, int* izero);
void            decode_UUID(const uuid_t uuid_in, long long unsigned int* f64, long long unsigned int* f48, 
                            unsigned char* izero);

#endif /* UTIL_UUID_H */

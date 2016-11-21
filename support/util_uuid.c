/* 128 bit UUID computation / Rabin fingerprinting algorithm                                        */
/*                                                                                                  */
/* ("A fingerprinting algorithm is a procedure that maps an arbitrarily large data item             */
/* to a much shorter bit string" [Wikipedia-en]).                                                   */
/*                                                                                                  */
/* The fingerprint of string A is computed as                                                       */
/*   f(A) = A(t) mod P(t)   where P(t) is an irreducible polynomial                                 */
/*                          in the ring of integers modulus 2 ("Z2").                               */
/*                                                                                                  */
/* @author F. Prill, DWD (2015-01-20)                                                               */
/*                                                                                                  */
/* Literature:                                                                                      */
/* [1] Rabin, M. O.:      "Fingerprinting by Random Polynomials."                                   */
/*                          Center for Research in Computing Technology, Harvard University., 1981  */
/* [2] Broder, A. Z.:     "Some applications of Rabin's fingerprinting method".                     */
/*                         Sequences II: Methods in Communications, Security, and Computer Science, */  
/*                         Springer-Verlag, 1993, 143-152                                           */
/* [3] Chan, C. & Lu, H.: "Fingerprinting using Polynomial (Rabin's method)"                        */
/*                         CMPUT690 Term Project, 2001                                              */
/* [4] IEEE754 floating point numbers: see, e.g., lecture notes by M. Hill,                         */
/*     http://pages.cs.wisc.edu/~markhill/cs354/Fall2008/notes/flpt.apprec.html                     */
/*                                                                                                  */
/* Remarks:                                                                                         */
/*                                                                                                  */
/* 1. The fingerprint ("UUID") can be used to distinguish between sequences of double precision     */
/*    64 bit values with a probability of failure of                                                */
/*      Pr(f(A) = f(B) | A != B) <= |length of sequence|/2^57                                       */
/*                                                                                                  */
/* 2. The fingerprint can also be used to determine if number sequences are equal up to a           */
/*    precision of 2^-24 ~ 5.96e-8, provided the numbers are known to be smaller than 2^40.         */
/*    The error probability for this is bound by                                                    */
/*      Pr(f(A) = f(B) | A != B) <= |length of sequence|/2^41                                       */
/*                                                                                                  */
/* 3. For identical sequences the algorithm deterministically generates the same UUID.              */
/*                                                                                                  */
/* 4. The algorithm generates 128 bit fingerprints that do not contain NaN or Inf values in their   */
/*    representation as four 32 bit floating point numbers.                                         */
/*                                                                                                  */
/* 5. As a slight drawback, the fingerprint is not cryptographically secure, i.e. an adversary      */
/*    is able to alter the data but to maintain the signature.                                      */
/*                                                                                                  */
/* 6. This implementation assumes a little-endian hardware architecture.                            */
/*                                                                                                  */
/* 7. The fingerprinting technique could easily allow for a cascade of accuracy tests, if more      */
/*    than 128 bits were available for the UUID.                                                    */
/*                                                                                                  */
/*                                                                                                  */
/* @par Copyright and License                                                                       */
/*                                                                                                  */
/* This code is subject to the DWD and MPI-M-Software-License-Agreement in                          */
/* its most recent form.                                                                            */
/* Please see the file LICENSE in the root of the source tree for this code.                        */
/* Where software is supplied by third parties, it is indicated in the                              */
/* headers of the routines.                                                                         */

#include "util_uuid.h"

/* -------------------- General Utility Routines -------------------- */

#define BYTETOBINARYPATTERN "%d%d%d%d%d%d%d%d"
#define BYTETOBINARY(byte)  \
  (byte & 0x80 ? 1 : 0), \
  (byte & 0x40 ? 1 : 0), \
  (byte & 0x20 ? 1 : 0), \
  (byte & 0x10 ? 1 : 0), \
  (byte & 0x08 ? 1 : 0), \
  (byte & 0x04 ? 1 : 0), \
  (byte & 0x02 ? 1 : 0), \
  (byte & 0x01 ? 1 : 0) 


// Display integer number as a binary string.
//
void printBinary(const unsigned int* sWord, unsigned int len) 
{
  unsigned char* s = ((unsigned char*) sWord)+4*len-1;
  for (unsigned int i=0; i<4*len; i++,s--)
    printf(BYTETOBINARYPATTERN"  ", BYTETOBINARY(*s));
  printf("\n");
}

// Display integer number as a binary string.
//
void printBinary64(const unsigned long long val) 
{
  unsigned char* s = ((unsigned char*) &val)+7;
  for (int i=0; i<8; i++,s--)
    printf(BYTETOBINARYPATTERN"  ", BYTETOBINARY(*s));
  printf("\n");
}


/* ---------- Auxiliary Routines for Z2 Arithmetics ----------------- */

// Perform 64 bit right-shift operation.
//
void longRightShift(unsigned int* sWord, unsigned int s) 
{
  if (s > 8) {
    longRightShift(sWord, 8);
    longRightShift(sWord, s-8);
  } else {
    long long unsigned int*  w = (long long unsigned int*) sWord;
    (*w) = (*w) >> s;
  }
}


// Perform 64 bit left-shift operation.
//
void longLeftShift(unsigned int* sWord, unsigned int s) 
{
  if (s > 8) {
    longLeftShift(sWord, 8);
    longLeftShift(sWord, s-8);
  } else {
    long long unsigned int*  w = (long long unsigned int*) sWord;
    (*w) = (*w) << s;
  }
}


// Generate irreducible polynomial of degree k
//
// see [Chan&Lu], Section 2.2.2
//
void generate_polynomial(unsigned int k, unsigned int* W)
{
  srand(time(NULL));

  unsigned int r = rand() % ((k-1)/2);
  r = 2*r + 1;
  unsigned int A[k];
  A[0] = 1;
  for (unsigned int i=1; i<k; i++) A[i] = 0;

  for (unsigned int i=1; i<r; i++) {
    unsigned int j = 0;
    while (A[j] != 0)  j = rand() % (k-1);
    A[j+1] = 1;
  }
  const unsigned int len = (k+31)/32;
  for (unsigned int i=0; i<len; i++) W[i] = 0;
  for (int i=(k-1); i>=0; i--) {
    longLeftShift(W, 1);
    W[0] += A[i];
  }
}


// Compute i x t^8 mod P(t) in the ring of integers modulus 2.
//
// This variant operates on 64 bit numbers, represented by two 32 bit
// integers.
//
void multiT8_64(unsigned int* s, 
                  unsigned int* in_s, 
                  const unsigned int* p) 
{
  s[0] = in_s[0];
  s[1] = in_s[1];
  for (int i=0; i<8; i++) {
    int needXOR = (s[1] & 0x80000000) == 0x80000000;  // conditional XOR 

    // left shift 1 bit: extract highest bit from lower word, insert
    // as lowest bit in higher word
    unsigned int temp = (s[0] >> 31) & 0x00000001;
    s[1] = (s[1] << 1) | temp;
    s[0] =  s[0] << 1;

    if (needXOR) {
      s[0] ^= p[0];
      s[1] ^= p[1];
    }
  }
}


// Compute i x t^8 mod P(t) in the ring of integers modulus 2.
//
// This variant operates on 48 bit numbers, represented by "one and a
// half" 32 bit integer. Therefore, the left shift operation explicitly
// sets the output bit to zero.
//
void multiT8_48(unsigned int* s, 
                  unsigned int* in_s, 
                  const unsigned int* p) 
{
  s[0] = in_s[0];
  s[1] = in_s[1];
  for (int i=0; i<8; i++) {
    int needXOR = (s[1] & 0x00008000) == 0x00008000;  // conditional XOR 
    // left shift 1 bit: extract highest bit from lower word, insert
    // as lowest bit in higher word
    unsigned int temp = (s[0] >> 31) & 0x00000001;
    s[1]  = (s[1] << 1) | temp;
    s[1] &= 0xFFFF; // clear old leading bit
    s[0]  =  s[0] << 1;

    if (needXOR) {
      s[0] ^= p[0];
      s[1] ^= p[1];
    }
  }
}


/* ---------- Fingerprint Routines, 48 bit and 64 bit --------------- */

// Pre-compute tables to process 4 bytes at once for 64 bit
// fingerprinting.
//
void compute_table_64(const unsigned int*  p, 
		      unsigned int TA[256][2], 
		      unsigned int TB[256][2], 
		      unsigned int TC[256][2], 
		      unsigned int TD[256][2]) 
{
  for (int i=0; i<256; i++) {
    // compute TD[i] = i x t^(k-8)
    //         TD[i] = i x t^8 mod P(t)
    TD[i][0] = 0x00000000;
    TD[i][1] = i << 24;
    multiT8_64(TD[i], TD[i], p);  
    multiT8_64(TC[i], TD[i], p);
    multiT8_64(TB[i], TC[i], p);
    multiT8_64(TA[i], TB[i], p);
  }
}


// Pre-compute tables to process 4 bytes at once for 32 bit
// fingerprinting.
//
// Compared to the 64 bit fingerprint we need two additional tables: 
//    TE(t) = F_c5(t)*t^56 mod P(t)  and 
//    TF(t) = F_c6(t)*t^48 mod P(t)
//
void compute_table_48(const unsigned int*  p, 
		      unsigned int TC[256][2], 
		      unsigned int TD[256][2],
		      unsigned int TE[256][2], 
		      unsigned int TF[256][2]) 
{
  for (int i=0; i<256; i++) {
    // compute TF[i] = i x t^40
    //         TF[i] = i x t^8 mod P(t)
    TF[i][0] = 0x00000000;
    TF[i][1] = i << 8;
    multiT8_48(TF[i], TF[i], p);  
    multiT8_48(TE[i], TF[i], p);  
    multiT8_48(TD[i], TE[i], p);  
    multiT8_48(TC[i], TD[i], p);
  }
}


// Fingerprint 4 characters (bytes) at once, 64 bit fingerprint.
//
void fp4_64(unsigned int* s, const unsigned char* w, 
	    unsigned int TA[256][2], 
	    unsigned int TB[256][2], 
	    unsigned int TC[256][2], 
	    unsigned int TD[256][2]) 
{
  unsigned int f[4];
  unsigned int j = 24, ww = 0;
  for (int i=0; i<4; i++,j-=8) {
    f[i]  =  (s[1] >> j) & 0xff;
    ww    =  (ww << 8) | w[i];
  }
  s[1] = s[0] ^ TA[f[0]][1] ^ TB[f[1]][1] ^ TC[f[2]][1] ^ TD[f[3]][1];
  s[0] = ww   ^ TA[f[0]][0] ^ TB[f[1]][0] ^ TC[f[2]][0] ^ TD[f[3]][0];
}


// Fingerprint 4 characters at once, 48 bit fingerprint.
//
void fp4_48(unsigned int* s, const unsigned char* w, 
	    unsigned int TC[256][2], 
	    unsigned int TD[256][2],
	    unsigned int TE[256][2], 
	    unsigned int TF[256][2]) 
{
  unsigned int f[6];
  unsigned int j = 24, ww = 0;
  for (int i=0; i<4; i++,j-=8) {
    f[i]  =  (s[1] >> j) & 0xff;
    ww    =  (ww << 8) | w[i];
  }
  f[4] = (s[0] >> 24) & 0xff;
  f[5] = (s[0] >> 16) & 0xff;
  unsigned int s0 = s[0] & 0xFFFF;
  s[1] = s0  ^ TC[f[2]][1] ^ TD[f[3]][1] ^ TE[f[4]][1] ^ TF[f[5]][1];
  s[0] = ww  ^ TC[f[2]][0] ^ TD[f[3]][0] ^ TE[f[4]][0] ^ TF[f[5]][0];
}


// Fingerprint 1 character, 64 bit fingerprint version.
//
// Computing the fingerprint extended by one bit consists of one shift
// left operation with the new bit as input bit and r1 as output bit and
// then, conditioned upon r1=1, a bit-wise exclusive operation, the
// second operand being P with the leading coefficient removed.
//
void fp1_64(unsigned int* s, unsigned char w, 
             const unsigned int* p) 
{
  for (int i=0; i<8; i++) {
    int needXOR = (s[1] & 0x80000000) == 0x80000000;

    // left shift 1 bit
    unsigned int temp = (s[0] >> 31) & 0x00000001;
    s[1] = (s[1] << 1) | temp;
    s[0] =  s[0] << 1;

    s[0] |= ((w & 0x80) >> 7);
    w = w << 1;
    if (needXOR) {
      s[0] ^= p[0];
      s[1] ^= p[1];
    }
  }
}


// Fingerprint 1 character, 48 bit fingerprint version.
//
void fp1_48(unsigned int* s, unsigned char w, 
             const unsigned int* p) 
{
  for (int i=0; i<8; i++) {
    int needXOR = s[1] & 0x00008000;

    // left shift 1 bit
    unsigned int temp = (s[0] >> 31) & 0x00000001;
    s[1]  = (s[1] << 1) | temp;
    s[1] &= 0xFFFF; // clear old leading bit
    s[0]  =  s[0] << 1;

    s[0] |= ((w & 0x80) >> 7);
    w = w << 1;
    if (needXOR) {
      s[0] ^= p[0];
      s[1] ^= p[1];
    }
  }
}


/* ----------------IEEE float to fixed point conversion ------------- */

// Extract unbiased exponent from 64 bit IEEE floating point number
//
int extract_exponent64(const double val)
{
  const unsigned int EXP_BIAS = 1023;
  unsigned char*     word     = ((unsigned char*) &val);
  return (((unsigned int) (word[7] & 0x7F)) << 4) + (word[6] >> 4) - EXP_BIAS;
}


// Extract unbiased exponent from 32 bit IEEE floating point number
//
int extract_exponent32(const float val)
{
  const unsigned int EXP_BIAS = 127;
  unsigned char*     word     = ((unsigned char*) &val);
  return (((int) word[3] & 0x7F) << 1) + (word[2] >> 7) - EXP_BIAS;
}


// Convert IEEE float to a 64 bit fixed point number.
//
// *IMPORTANT Assumption: Little-Endianness*
//
// IEEE:         [   1 Sign bit                                   ]
//               [ + 11-bit, biased-1023 integer for the exponent ]
//               [ + 52-bit mantissa                              ]
//
// Fixed-Point:  [ Sign bit + 2^39, ..., 2^-23 ]
//
// Parameter "izero": Consider the part of the mantissa that is
// truncated by converting a number~$x$ to its fixed-point
// representation and let $q(x)$ denote the largest bit position with a
// zero coefficient.
//
void convert_to_fixed64(const double val, unsigned long long* fp, int* izero)
{
  // split floating point number into its sign, exponent and mantissa:
  unsigned char*     word     = ((unsigned char*) &val);
  unsigned int       sign     = ((unsigned int) (word[7] & 0x80));  // note: sign is at MSB position
  int                exponent = extract_exponent64(val);

  // re-order normalized little-endian mantissa
  unsigned long long mantissa = 
    0x0010000000000000ull                         + // hidden bit
    (((unsigned long long) (word[6] & 0x0F)) << 48) + 
    (((unsigned long long) word[5])          << 40) + 
    (((unsigned long long) word[4])          << 32) + 
    (((unsigned long long) word[3])          << 24) + 
    (((unsigned long long) word[2])          << 16) + 
    (((unsigned long long) word[1])          <<  8) +
    ((unsigned long long) word[0]);

  // shift mantissa to the right position in the fixed point representation
  const int smallest_exp = -23;
  //                sign bit    shift to 2^0 location 
  const int shift = 1         + 63-12+smallest_exp        - exponent;

  (*izero) = 0;
  if (shift > 0) 
    {
      // determine highest-value zero in the truncated part of the
      // mantissa; performance could certainly be improved by using a
      // lookup table.
      unsigned long long bitmask = 0, testbit=1;
      for (int j=0; j<shift; ++j, testbit<<=1)
	bitmask = (bitmask << 1) + 1;
      for (*izero=0; *izero<shift; ++(*izero)) {
	testbit >>= 1;	
	if ((mantissa & testbit) == 0) break;
      }
      longRightShift((unsigned int*) &mantissa, shift);
    }
  else if (shift < 0)
    longLeftShift( (unsigned int*) &mantissa, -1*shift);
  (*fp) = mantissa;

  // set sign bit
  *((unsigned char*) fp+7) |= sign;
  (*izero) -= smallest_exp - 1;
}


/* ----------------------- UUID Computation ------------------------- */

// irreducible polynomial x^64 + x^4 + x^3 +   x + 1
// with the leading coefficient removed.
const unsigned int P64[2] = { 0x0000001B, 0x00000000 };
const unsigned int P48[2] = { 0xaf5433a5, 0xb09a };


// Print UUID in standard hexadecimal format.
//
void uuid_unparse(char *buffer, const uuid_t *uuid)
{
  const unsigned char *d = uuid->data;

  sprintf(buffer,
      "%02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x", 
           d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7], d[8], d[9], d[10], d[11], d[12], d[13], d[14], d[15]);
}


unsigned char parse_hexpair(const char *s)
{
  int result;
  if (isdigit(*s))
    result = (*s - '0') << 4;
  else if (isupper(*s))
    result = (*s - 'A' + 10) << 4;
  else
    result = (*s - 'a' + 10) << 4;
  ++s;
  if (isdigit(*s))
    result |= (*s - '0');
  else if (isupper(*s))
    result |= (*s - 'A' + 10);
  else
    result |= (*s - 'a' + 10);
  return (unsigned char) result;
}


// Convert a hex-string representation of a UUID into a uuid byte
// array.
//
int uuid_parse(uuid_t *uuid, const char *uuid_str)
{
  int i;
  unsigned char *d = uuid->data;

  for (i = 0; i < 36; ++i) {
      char c = uuid_str[i];
      if (!isxdigit(c)
          && !(c == '-' && (i == 8 || i == 13 || i == 18 || i == 23)))
        return -1;
    }
  if (uuid_str[36] != '\0')  return -1;

  d[0] = parse_hexpair(&uuid_str[0]);
  d[1] = parse_hexpair(&uuid_str[2]);
  d[2] = parse_hexpair(&uuid_str[4]);
  d[3] = parse_hexpair(&uuid_str[6]);

  d[4] = parse_hexpair(&uuid_str[9]);
  d[5] = parse_hexpair(&uuid_str[11]);

  d[6] = parse_hexpair(&uuid_str[14]);
  d[7] = parse_hexpair(&uuid_str[16]);

  d[8] = parse_hexpair(&uuid_str[19]);
  d[9] = parse_hexpair(&uuid_str[21]);

  for (i = 6; i--;)
    d[10 + i] = parse_hexpair(&uuid_str[i * 2 + 24]);
  return 0;
}


// Auxiliary routine for "uuid_generate": add data to 48 bit and 64
// bit fingerprint.
//
void uuid_scan_data(const double* val, const int nval, context_t* fprnt, int *max_zero)
{
  // pre-compute tables for fingerprint algorithm
  if (fprnt->initialized == 0) {
    compute_table_64(P64, fprnt->TA64, fprnt->TB64, fprnt->TC64, fprnt->TD64);
    compute_table_48(P48, fprnt->TC48, fprnt->TD48, fprnt->TE48, fprnt->TF48);
    fprnt->f48[0] = fprnt->f48[1] = 0; // 48 bit fingerprint
    fprnt->f64[0] = fprnt->f64[1] = 0; // 64 bit fingerprint

    fprnt->initialized = 1;
  }

  // --- loop over the sequence of double precision floating point
  //     numbers, let "t" denote the current float
  unsigned long long val_fp;
  int izero;
  for (int i=0; i<nval; i++)
    {
      // --- --- update the 64 bit fingerprint F1 by the current float t
      fp4_64(fprnt->f64, (unsigned char*) &val[i],     fprnt->TA64, fprnt->TB64, fprnt->TC64, fprnt->TD64);
      fp4_64(fprnt->f64, (unsigned char*) &val[i] + 4, fprnt->TA64, fprnt->TB64, fprnt->TC64, fprnt->TD64);
      // --- --- translate float t into 64 bit *fixed* point number t'
      convert_to_fixed64(val[i], &val_fp, &izero);
      if ((*max_zero) <= izero)  (*max_zero)=izero;
      // --- --- update 48 bit fingerprint F2 with t'
      fp4_48(fprnt->f48, (unsigned char*) &val_fp,     fprnt->TC48, fprnt->TD48, fprnt->TE48, fprnt->TF48);
      fp4_48(fprnt->f48, (unsigned char*) &val_fp + 4, fprnt->TC48, fprnt->TD48, fprnt->TE48, fprnt->TF48);
    }
}


// Given a sequence of 64 bit double precision floating point numbers,
// the UUID (Universally Unique Identifier) is computed here as follows.
//
// The 128 bit (16 byte UUID or four 32 bit REALs) consists of three parts:
//
// bits   0- 63 : 64 bit fingerprint
//                of the floating point number sequence.
// bits  64-111 : 48 bit fingerprint
//                of a truncated fixed-point number sequence.
// bits 112-119 : largest bit position with a zero coefficient (max over all numbers)
// bits 120-123 : 4 bits to avoid NaN/Inf in the representation
//                of the UUID as four single precision floating point
//                numbers (these bits are moved out of the REALs' exponents 
//                where they are set accordingly if a NaN/Inf needs to be 
//                avoided).
//
void uuid_generate(const double* val, const int nval, uuid_t* uuid) 
{
  context_t context;
  context.initialized = 0;

  // --- loop over the sequence of floating point numbers
  int max_zero = 0;
  uuid_scan_data(val, nval, &context, &max_zero);

  // --- concatenate the two fingerprints
  for (int i=0; i<16; ++i) uuid->data[i] = 0x00;
  unsigned char* words64 = (unsigned char*) context.f64;
  unsigned char* words48 = (unsigned char*) context.f48;
  for (int i=0; i<8; ++i) uuid->data[i]   = words64[i];
  for (int i=0; i<6; ++i) uuid->data[i+8] = words48[i];

  // --- cast (F1,F2) into four 32 bit floating point numbers,
  //     insert four bits, each into the exponent of the
  //     4-floats representation of the 128 bit string.
  unsigned int old_bits = 0;
  for (int i=0; i<4; i++) {
    unsigned char new_bit = (uuid->data[4*i+3] & 0x40) >> 6;
    old_bits = (old_bits << 1) + new_bit;
    float* data = (float*) &uuid->data[i];
    int exponent = extract_exponent32(*data);
    if (exponent == 0)
      uuid->data[4*i+3] |=  0x40;
    else
      uuid->data[4*i+3] &= ~0x40;
  }

  // maximum of highest-value zero in the truncated part of the
  // mantissa
  uuid->data[14] = max_zero;

  //     The original bits are appended at bit position 120-123
  //     fingerprints.
  uuid->data[15] = old_bits << 4;
}


/* ------------------------ UUID Decoding --------------------------- */

// Splits the UUID into its fingerprint parts.
//
void decode_UUID(const uuid_t uuid_in, 
		 long long unsigned int* f64,
		 long long unsigned int* f48,
		 unsigned char* izero)
{
  // --- insert bits 112 - 115 into locations 0x40 of bytes 3,7,11,15
  uuid_t uuid = uuid_in;
  unsigned int old_bits = uuid.data[15];
  for (int i=0; i<4; i++) {
    if ((old_bits & 0x80) == 0x80)
      uuid.data[4*i+3] |=  0x40;  // set bit
    else
      uuid.data[4*i+3] &= ~0x40;  // unset bit
    old_bits = old_bits << 1;
  }
  // --- split UUID into 64 bit and 48 bit fingerprint
  (*f64) = 0;
  (*f48) = 0;
  unsigned char* words64 = (unsigned char*) f64;
  unsigned char* words48 = (unsigned char*) f48;
  for (int i=0; i<8; ++i) words64[i] = uuid.data[i];
  for (int i=0; i<6; ++i) words48[i] = uuid.data[i+8];
  (*izero) = uuid.data[14];
}


// Compare two UUIDs.
//
cmp_UUID_t compare_UUID(const uuid_t uuid_A, const uuid_t uuid_B, 
			double* min_difference)
{
  // decode the full-precision and the low-accuracy
  // fingerprint from the UUID:
  long long unsigned int f64_A, f48_A, f64_B, f48_B;
  unsigned char          q_A, q_B;
  decode_UUID(uuid_A, &f64_A, &f48_A, &q_A);
  decode_UUID(uuid_B, &f64_B, &f48_B, &q_B);

  (*min_difference) = -1.;

  if (f64_A == f64_B) 
    return UUID_EQUAL;
  else if (f48_A == f48_B) 
    return UUID_EQUAL_LIMITED_ACCURACY;
  else {
    unsigned char q = (q_A > q_B) ? q_A : q_B;
    (*min_difference) = exp(-1.*q*log(2.));
    return UUID_UNEQUAL;
  }
}


void
uuid_format(char *buffer, const uuid_t *uuid)
{
  const unsigned char *d = uuid->data;

  sprintf(buffer,
      "%02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x", d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7], d[8], d[9], d[10], d[11], d[12], d[13], d[14], d[15]);
}

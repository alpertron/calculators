//
// This file is part of Alpertron Calculators.
//
// Copyright 2026 Dario Alejandro Alpern
//
// Alpertron Calculators is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Alpertron Calculators is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Alpertron Calculators.  If not, see <http://www.gnu.org/licenses/>.
//

#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include "bignbr.h"

#define MAX_FFT_BITS 18
#define MAX_FFT_LEN (1 << MAX_FFT_BITS)

/* ---------- 3 NTT primes ---------- */

#define NTTPrime1 998244353U
#define NTTPrime2 1004535809U
#define NTTPrime3 469762049U

#define pInv1 998244351U
#define pInv2 1004535807U
#define pInv3 469762047U

#define R2Prime1 932051910U
#define R2Prime2 542374313U
#define R2Prime3 460175152U

#define oneMont1 301989884U
#define oneMont2 276824060U
#define oneMont3 67108855U

#define primitiveRoot1 905969652U
#define primitiveRoot2 830472180U
#define primitiveRoot3 201326565U

/* ---------- Triple struct: one element per prime, padded to 16 bytes ---------- */

typedef struct
{
  uint32_t v1;
  uint32_t v2;
  uint32_t v3;
} triple_t;

static triple_t firstFactor[MAX_FFT_LEN];
static triple_t secondFactor[MAX_FFT_LEN];
static triple_t twiddleFactor[(MAX_FFT_LEN / 2) + 1];
static triple_t inv_n[MAX_FFT_BITS];
static triple_t MontgomeryMultNTransf[MAX_FFT_LEN];
static triple_t CustomNbrTransf[MAX_FFT_LEN];
static triple_t TestNbrTransf[MAX_FFT_LEN];
static int constantsInitializedFFTSize = 0;

/* ---------- Montgomery multiplication ---------- */

static inline uint32_t mont_mul_p1(uint32_t a, uint32_t b)
{
  uint64_t t = (uint64_t)a * b;
  uint32_t m0 = (uint32_t)t * pInv1;
  return (uint32_t)((t + (uint64_t)m0 * NTTPrime1) >> 32);  // Range: [0, 2p).
}

static inline uint32_t mont_mul_p2(uint32_t a, uint32_t b)
{
  uint64_t t = (uint64_t)a * b;
  uint32_t m0 = (uint32_t)t * pInv2;
  return (uint32_t)((t + (uint64_t)m0 * NTTPrime2) >> 32);  // Range: [0, 2p).
}

static inline uint32_t mont_mul_p3(uint32_t a, uint32_t b)
{
  uint64_t t = (uint64_t)a * b;
  uint32_t m0 = (uint32_t)t * pInv3;
  return (uint32_t)((t + (uint64_t)m0 * NTTPrime3) >> 32);  // Range: [0, 2p).
}

/* ---------- Power mod ---------- */

static uint32_t mod_pow_p1(uint32_t base, uint32_t expon)
{
  uint32_t power = oneMont1;
  while (expon != 0)
  {
    if (expon & 1)
    {
      power = mont_mul_p1(power, base);
    }
    base = mont_mul_p1(base, base);
    expon >>= 1;
  }
  return power;
}

static uint32_t mod_pow_p2(uint32_t base, uint32_t expon)
{
  uint32_t power = oneMont2;
  while (expon != 0)
  {
    if (expon & 1)
    {
      power = mont_mul_p2(power, base);
    }
    base = mont_mul_p2(base, base);
    expon >>= 1;
  }
  return power;
}

static uint32_t mod_pow_p3(uint32_t base, uint32_t expon)
{
  uint32_t power = oneMont3;
  while (expon != 0)
  {
    if (expon & 1)
    {
      power = mont_mul_p3(power, base);
    }
    base = mont_mul_p3(base, base);
    expon >>= 1;
  }
  return power;
}

/* ---------- Init all three primes together ---------- */

static void init_ntt(void)
{
  /* inv_n[k-1] = inverse of 2^k, outside Montgomery domain */
  for (int k = 1; k <= MAX_FFT_BITS; k++)
  {
    int len = 1 << k;
    uint32_t ml1 = mont_mul_p1(len, R2Prime1);
    uint32_t ml2 = mont_mul_p2(len, R2Prime2);
    uint32_t ml3 = mont_mul_p3(len, R2Prime3);
    uint32_t invMont1 = mod_pow_p1(ml1, NTTPrime1 - 2);
    inv_n[k - 1].v1 = mont_mul_p1(invMont1, 1U);
    uint32_t invMont2 = mod_pow_p2(ml2, NTTPrime2 - 2);
    inv_n[k - 1].v2 = mont_mul_p2(invMont2, 1U);
    uint32_t invMont3 = mod_pow_p3(ml3, NTTPrime3 - 2);
    inv_n[k - 1].v3 = mont_mul_p3(invMont3, 1U);
  }

  /* twiddle table: root of order MAX_FFT_LEN/2 */
  uint32_t root1 = mod_pow_p1(primitiveRoot1, (NTTPrime1 - 1) / (MAX_FFT_LEN / 2));
  uint32_t root2 = mod_pow_p2(primitiveRoot2, (NTTPrime2 - 1) / (MAX_FFT_LEN / 2));
  uint32_t root3 = mod_pow_p3(primitiveRoot3, (NTTPrime3 - 1) / (MAX_FFT_LEN / 2));
  uint32_t w1 = oneMont1;
  uint32_t w2 = oneMont2;
  uint32_t w3 = oneMont3;

  triple_t* ptrTwiddle = twiddleFactor;
  for (int j = 0; j <= MAX_FFT_LEN / 2; j++)
  {
    ptrTwiddle->v1 = w1;
    ptrTwiddle->v2 = w2;
    ptrTwiddle->v3 = w3;
    // Perform Montgomery multiplications on all three primes.
    uint64_t prod1 = (uint64_t)w1 * root1;
    uint64_t prod2 = (uint64_t)w2 * root2;
    uint64_t prod3 = (uint64_t)w3 * root3;
    uint32_t m1 = (uint32_t)prod1 * pInv1;
    uint32_t m2 = (uint32_t)prod2 * pInv2;
    uint32_t m3 = (uint32_t)prod3 * pInv3;
    w1 = (uint32_t)((prod1 + (uint64_t)m1 * NTTPrime1) >> 32);
    w2 = (uint32_t)((prod2 + (uint64_t)m2 * NTTPrime2) >> 32);
    w3 = (uint32_t)((prod3 + (uint64_t)m3 * NTTPrime3) >> 32);
    w1 -= NTTPrime1;
    w2 -= NTTPrime2;
    w3 -= NTTPrime3;
    w1 += ((int32_t)w1 >> 31) & NTTPrime1;
    w2 += ((int32_t)w2 >> 31) & NTTPrime2;
    w3 += ((int32_t)w3 >> 31) & NTTPrime3;
    ptrTwiddle++;
  }
  assert(twiddleFactor[MAX_FFT_LEN / 2].v1 == oneMont1);
  assert(twiddleFactor[MAX_FFT_LEN / 2].v2 == oneMont2);
  assert(twiddleFactor[MAX_FFT_LEN / 2].v3 == oneMont3);
}

/* ---------- DIF forward NTT, all three primes together ---------- */

static void ntt_dif(triple_t* signal, int convolutionLen)
{
  int step = (MAX_FFT_LEN / 2) / convolutionLen;
  for (int len = convolutionLen; len >= 2; len >>= 1)
  {
    int halfLen = len >> 1;

    for (int i = 0; i < convolutionLen; i += len)
    {
      triple_t* p1 = signal + i;
      triple_t* p2 = p1 + halfLen;
      const triple_t* ptrTwiddle = twiddleFactor;

      for (int j = 0; j < halfLen; j++)
      {
        uint32_t u1 = p1->v1;   // Range: [0, 2p) for all three primes.
        uint32_t u2 = p1->v2;
        uint32_t u3 = p1->v3;
        uint32_t t1 = p2->v1;   // Range: [0, 2p) for all three primes.
        uint32_t t2 = p2->v2;
        uint32_t t3 = p2->v3;

        /* sum */
        uint32_t s1 = u1 + t1 - (2 * NTTPrime1); // Range [-2p, 2p)
        p1->v1 = s1 + ((uint32_t)((int32_t)s1 >> 31) & (2 * NTTPrime1));  // Range [0, 2p)
        uint32_t s2 = u2 + t2 - (2 * NTTPrime2); // Range [-2p, 2p)
        p1->v2 = s2 + ((uint32_t)((int32_t)s2 >> 31) & (2 * NTTPrime2));  // Range [0, 2p)
        uint32_t s3 = u3 + t3 - (2 * NTTPrime3); // Range [-2p, 2p)
        p1->v3 = s3 + ((uint32_t)((int32_t)s3 >> 31) & (2 * NTTPrime3));  // Range [0, 2p)

        /* diff then twiddle */
        uint32_t d1 = u1 + (2 * NTTPrime1) - t1; // Range [0, 4p)
        uint32_t d2 = u2 + (2 * NTTPrime2) - t2; // Range [0, 4p)
        uint32_t d3 = u3 + (2 * NTTPrime3) - t3; // Range [0, 4p)
    
        uint64_t mul1 = (uint64_t)d1 * ptrTwiddle->v1;
        uint64_t mul2 = (uint64_t)d2 * ptrTwiddle->v2;
        uint64_t mul3 = (uint64_t)d3 * ptrTwiddle->v3;
        uint32_t m1 = (uint32_t)mul1 * pInv1;
        uint32_t m2 = (uint32_t)mul2 * pInv2;
        uint32_t m3 = (uint32_t)mul3 * pInv3;
        p2->v1 = (uint32_t)((mul1 + (uint64_t)m1 * NTTPrime1) >> 32); // Range [0, 2p).
        p2->v2 = (uint32_t)((mul2 + (uint64_t)m2 * NTTPrime2) >> 32); // Range [0, 2p).
        p2->v3 = (uint32_t)((mul3 + (uint64_t)m3 * NTTPrime3) >> 32); // Range [0, 2p).

        ptrTwiddle += step;
        p1++;
        p2++;
      }
    }
    step <<= 1;
  }
}

/* ---------- DIT inverse NTT, all three primes together ---------- */

static void getTransform(const limb* origNumber, triple_t* transf,
  int length, int convolutionLen)
{
  const limb* ptrOrigNumber = origNumber;
  triple_t* ptrTransf = transf;
  for (int i = 0; i < length; i++)
  {
    uint32_t x = ptrOrigNumber->x;
    // Perform Montgomery multiplications for all three primes. Range is [0, 2p).
    uint64_t prod1 = (uint64_t)(x % NTTPrime1) * R2Prime1;
    uint64_t prod2 = (uint64_t)(x % NTTPrime2) * R2Prime2;
    uint64_t prod3 = (uint64_t)(x % NTTPrime3) * R2Prime3;
    uint32_t m1 = (uint32_t)prod1 * pInv1;
    uint32_t m2 = (uint32_t)prod2 * pInv2;
    uint32_t m3 = (uint32_t)prod3 * pInv3;
    ptrTransf->v1 = (uint32_t)((prod1 + (uint64_t)m1 * NTTPrime1) >> 32);
    ptrTransf->v2 = (uint32_t)((prod2 + (uint64_t)m2 * NTTPrime2) >> 32);
    ptrTransf->v3 = (uint32_t)((prod3 + (uint64_t)m3 * NTTPrime3) >> 32);
    ptrOrigNumber++;
    ptrTransf++;
  }
  for (int i = length; i < convolutionLen; i++)
  {
    ptrTransf->v1 = 0;
    ptrTransf->v2 = 0;
    ptrTransf->v3 = 0;
    ptrTransf++;
  }
  ntt_dif(transf, convolutionLen);
}

static int ntt_dit_inverse(triple_t* signal, int convolutionLen)
{
  int level = 0;
  int step = MAX_FFT_LEN / 2;
  assert(convolutionLen >= 2);
  for (int len = 2; len <= convolutionLen; len <<= 1)
  {
    int halfLen = len >> 1;
    step >>= 1;

    for (int i = 0; i < convolutionLen; i += len)
    {
      triple_t* p1 = signal + i;
      triple_t* p2 = p1 + halfLen;
      const triple_t* ptrTwiddle = &twiddleFactor[MAX_FFT_LEN / 2];

      for (int j = 0; j < halfLen; j++)
      {
        uint32_t u1 = p1->v1;   // Range [0, 2p).
        uint32_t u2 = p1->v2;   // Range [0, 2p).
        uint32_t u3 = p1->v3;   // Range [0, 2p).

        uint64_t mul1 = (uint64_t)p2->v1 * ptrTwiddle->v1;
        uint64_t mul2 = (uint64_t)p2->v2 * ptrTwiddle->v2;
        uint64_t mul3 = (uint64_t)p2->v3 * ptrTwiddle->v3;
        uint32_t m1 = (uint32_t)mul1 * pInv1;
        uint32_t m2 = (uint32_t)mul2 * pInv2;
        uint32_t m3 = (uint32_t)mul3 * pInv3;
        uint32_t v1 = (uint32_t)((mul1 + (uint64_t)m1 * NTTPrime1) >> 32);  // Range [0, 2p).
        uint32_t v2 = (uint32_t)((mul2 + (uint64_t)m2 * NTTPrime2) >> 32);  // Range [0, 2p).
        uint32_t v3 = (uint32_t)((mul3 + (uint64_t)m3 * NTTPrime3) >> 32);  // Range [0, 2p).

        uint32_t s1 = u1 + v1 - (2 * NTTPrime1);  // Range [-2p, 2p).
        uint32_t s2 = u2 + v2 - (2 * NTTPrime2);  // Range [-2p, 2p).
        uint32_t s3 = u3 + v3 - (2 * NTTPrime3);  // Range [-2p, 2p).
  
        p1->v1 = s1 + ((uint32_t)((int32_t)s1 >> 31) & (2 * NTTPrime1));    // Range [0, 2p).
        p1->v2 = s2 + ((uint32_t)((int32_t)s2 >> 31) & (2 * NTTPrime2));    // Range [0, 2p).
        p1->v3 = s3 + ((uint32_t)((int32_t)s3 >> 31) & (2 * NTTPrime3));    // Range [0, 2p).

        uint32_t d1 = u1 - v1;  // Range [-2p, 2p).
        uint32_t d2 = u2 - v2;  // Range [-2p, 2p).
        uint32_t d3 = u3 - v3;  // Range [-2p, 2p).
        p2->v1 = d1 + ((uint32_t)((int32_t)d1 >> 31) & (2 * NTTPrime1));    // Range [0, 2p).
        p2->v2 = d2 + ((uint32_t)((int32_t)d2 >> 31) & (2 * NTTPrime2));    // Range [0, 2p).
        p2->v3 = d3 + ((uint32_t)((int32_t)d3 >> 31) & (2 * NTTPrime3));    // Range [0, 2p).
     
        ptrTwiddle -= step;
        p1++;
        p2++;
      }
    }
    level++;
  }
  return level;
}

/* ---------- Merged convolution for all three primes ---------- */

static void convolution_merged(const limb* first, const limb* second,
  int len1, int len2, int convolutionLen)
{
  const triple_t* ptrSecondFactor;
  getTransform(first, firstFactor, len1, convolutionLen);
  bool accessingCustomNbr = (CustomNbrCached == NBR_CACHED) && (second == CustomNbrAddr);
  bool accessingTestNbr = (TestNbrCached == NBR_CACHED) && (second == TestNbr);
  bool accessingMontgomeryMultN = (MontgomeryMultNCached == NBR_CACHED) && (second == MontgomeryMultN);
  int nbrBytes = convolutionLen * (int)sizeof(triple_t);

  if (first != second)
  {
    if (accessingTestNbr)
    {
      ptrSecondFactor = TestNbrTransf;
    }
    else if (accessingMontgomeryMultN)
    {
      ptrSecondFactor = MontgomeryMultNTransf;
    }
    else if (accessingCustomNbr)
    {
      ptrSecondFactor = CustomNbrTransf;
    }
    else
    {
      getTransform(second, secondFactor, len2, convolutionLen);
      if ((TestNbrCached == NBR_READY_TO_BE_CACHED) && (second == TestNbr))
      {
        (void)memcpy(TestNbrTransf, secondFactor, nbrBytes);
        TestNbrCached = NBR_CACHED;
      }
      else if ((MontgomeryMultNCached == NBR_READY_TO_BE_CACHED) && (second == MontgomeryMultN))
      {
        (void)memcpy(MontgomeryMultNTransf, secondFactor, nbrBytes);
        MontgomeryMultNCached = NBR_CACHED;
      }
      else if ((CustomNbrCached == NBR_READY_TO_BE_CACHED) && (second == CustomNbrAddr))
      {
        (void)memcpy(CustomNbrTransf, secondFactor, nbrBytes);
        CustomNbrCached = NBR_CACHED;
      }
      else
      {  // No more conditions.
      }
      ptrSecondFactor = secondFactor;
    }
  }
  else
  {
    ptrSecondFactor = firstFactor;
  }

  triple_t *ptrFirstFactor = firstFactor;
  for (int i = 0; i < convolutionLen; i++)
  {
    // Perform Montgomery multiplications for all three primes. Output range is [0, 2p).
    uint64_t prod1 = (uint64_t)ptrFirstFactor->v1 * ptrSecondFactor->v1;
    uint64_t prod2 = (uint64_t)ptrFirstFactor->v2 * ptrSecondFactor->v2;
    uint64_t prod3 = (uint64_t)ptrFirstFactor->v3 * ptrSecondFactor->v3;
    uint32_t m1 = (uint32_t)prod1 * pInv1;
    uint32_t m2 = (uint32_t)prod2 * pInv2;
    uint32_t m3 = (uint32_t)prod3 * pInv3;
    ptrFirstFactor->v1 = (uint32_t)((prod1 + (uint64_t)m1 * NTTPrime1) >> 32);
    ptrFirstFactor->v2 = (uint32_t)((prod2 + (uint64_t)m2 * NTTPrime2) >> 32);
    ptrFirstFactor->v3 = (uint32_t)((prod3 + (uint64_t)m3 * NTTPrime3) >> 32);
    ptrFirstFactor++;
    ptrSecondFactor++;
  }

  int level = ntt_dit_inverse(firstFactor, convolutionLen);
  triple_t inv = inv_n[level - 1];
  ptrFirstFactor = firstFactor;
  for (int i = 0; i < convolutionLen; i++)
  {
    // Perform Montgomery multiplications for all three primes. Output range is [0, p).
    uint64_t prod1 = (uint64_t)ptrFirstFactor->v1 * inv.v1;
    uint64_t prod2 = (uint64_t)ptrFirstFactor->v2 * inv.v2;
    uint64_t prod3 = (uint64_t)ptrFirstFactor->v3 * inv.v3;
    uint32_t m1 = (uint32_t)prod1 * pInv1;
    uint32_t m2 = (uint32_t)prod2 * pInv2;
    uint32_t m3 = (uint32_t)prod3 * pInv3;
    uint32_t v1 = (uint32_t)((prod1 + (uint64_t)m1 * NTTPrime1) >> 32);
    uint32_t v2 = (uint32_t)((prod2 + (uint64_t)m2 * NTTPrime2) >> 32);
    uint32_t v3 = (uint32_t)((prod3 + (uint64_t)m3 * NTTPrime3) >> 32);
    v1 -= NTTPrime1;
    v2 -= NTTPrime2;
    v3 -= NTTPrime3;
    ptrFirstFactor->v1 = v1 + (((int32_t)v1 >> 31) & NTTPrime1);
    ptrFirstFactor->v2 = v2 + (((int32_t)v2 >> 31) & NTTPrime2);
    ptrFirstFactor->v3 = v3 + (((int32_t)v3 >> 31) & NTTPrime3);
    assert(ptrFirstFactor->v1 < NTTPrime1);
    assert(ptrFirstFactor->v2 < NTTPrime2);
    assert(ptrFirstFactor->v3 < NTTPrime3);
    ptrFirstFactor++;
  }
}

/* ---------- Chinese Remainder (Garner simplified) ---------- */

static void crt3_96(uint32_t r1, uint32_t r2, uint32_t r3,
  uint64_t* lo, uint64_t* hi)
{
  const uint64_t p1 = NTTPrime1;
  const uint64_t p2 = NTTPrime2;
  const uint64_t p3 = NTTPrime3;

  const uint32_t inv_p1_mod_p2 = 669690699U;
  const uint32_t inv_p12_mod_p3 = 354521948U;

  uint64_t t1 = r1;

  uint64_t tmp = r2 + p2 - (t1 % p2);
  if (tmp >= p2) tmp -= p2;
  uint64_t t2 = (tmp * inv_p1_mod_p2) % p2;

  uint64_t middle = t1 + p1 * t2;
  uint64_t middle_mod_p3 = middle % p3;

  tmp = r3 + p3 - middle_mod_p3;
  if (tmp >= p3) tmp -= p3;
  uint64_t t3 = (tmp * inv_p12_mod_p3) % p3;

  uint64_t lo2 = p1 * t2;
  uint64_t hi2 = (p1 * t2) >> 63 >> 1;

  uint64_t lo3 = p1 * p2;
  uint64_t a_lo = (uint32_t)lo3;
  uint64_t a_hi = lo3 >> 32;
  uint64_t b_lo = (uint32_t)t3;
  uint64_t b_hi = t3 >> 32;
  uint64_t p0 = a_lo * b_lo;
  uint64_t p1_ = a_lo * b_hi;
  uint64_t p2_ = a_hi * b_lo;
  uint64_t p3_ = a_hi * b_hi;
  uint64_t mid = p1_ + p2_;
  uint64_t carry = (mid < p1_);
  uint64_t lo_prod = p0 + (mid << 32);
  carry += (lo_prod < p0);
  uint64_t hi_prod = p3_ + (mid >> 32) + carry;

  uint64_t lo_sum = t1 + lo2;
  uint64_t carry2 = (lo_sum < t1);
  lo_sum += lo_prod;
  carry2 += (lo_sum < lo_prod);

  *lo = lo_sum;
  *hi = hi2 + hi_prod + carry2;
}

/* ============================================================ */
/*                   MAIN ENTRY FUNCTION                        */
/* ============================================================ */

void nttMultiplication(const limb* factor1,
  const limb* factor2,
  limb* result,
  int len1,
  int len2,
  int* pResultLen)
{
  assert(len1 >= 1);
  assert(len2 >= 1);
  int convolutionLen = 1;
  while (convolutionLen < len1 + len2)
  {
    convolutionLen <<= 1;
  }
  assert(convolutionLen <= MAX_FFT_LEN);

  if (convolutionLen != constantsInitializedFFTSize)
  {
    init_ntt();
    constantsInitializedFFTSize = convolutionLen;
  }

  convolution_merged(factor1, factor2, len1, len2, convolutionLen);

  uint64_t carry = 0;
  int total = len1 + len2;
  const triple_t* ptrProductTransf = firstFactor;
  limb* ptrResult = result;
  for (int i = 0; i < total; i++)
  {
    uint64_t lo;
    uint64_t hi;
    crt3_96(ptrProductTransf->v1, ptrProductTransf->v2, ptrProductTransf->v3, &lo, &hi);
    uint64_t value = lo + carry;
    carry = (value >> BITS_PER_GROUP) + (hi << (64 - BITS_PER_GROUP));
    ptrResult->x = value & MAX_INT_NBR_U;
    ptrProductTransf++;
    ptrResult++;
  }
  if (pResultLen != NULL)
  {
    *pResultLen = total;
  }
}
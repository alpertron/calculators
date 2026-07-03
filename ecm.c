//
// This file is part of Alpertron Calculators.
//
// Copyright 2015-2021 Dario Alejandro Alpern
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
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include "string/strings.h"
#include "bignbr.h"
#include "expression.h"
#include "factor.h"
#include "commonstruc.h"
#include "skiptest.h"
#include "copyStr.h"
#include "ecmInternal.h"

#define SIEVE_COMPOSITE        0
#define SIEVE_PROBABLE_PRIME   1

#define NOT_COPRIME_2310   SIEVE_COMPOSITE
#define COPRIME_2310       SIEVE_PROBABLE_PRIME

#ifdef __EMSCRIPTEN__
extern int intPrime;
extern int intStep1Bound;
extern int64_t longPrime;
extern int64_t longStep2Bound;
extern char* ptrLowerText;
extern char lowerText[MAX_LEN * 16];
#endif
extern int StepECM;
extern int maxIndexM;
extern int indexM;
static int SmallPrime[1335]; /* Primes < 11000 */
static int boundStep1;
static int64_t boundStep2;
static int sqrtBoundStep1;
static int NumberSizeBytes;
static int groupSize;
static int sieveSize;
static int halfSieveSize;

struct sBounds
{
  int digitLevel;
  int nbrCurves;
  int boundStep1;
  int sqrtBoundStep1;
};

const struct sBounds stEcmBounds[] =
{
  { 15, 25, 2000, 45 },                // ECM bounds for 15 digits
  { 20, 90, 11000, 105 },              // ECM bounds for 20 digits
  { 25, 300, 50000, 224 },             // ECM bounds for 25 digits
  { 30, 700, 250000, 501 },            // ECM bounds for 30 digits
  { 35, 1800, 1000000, 1001 },         // ECM bounds for 35 digits
  { 40, 5100, 3000000, 1733 },         // ECM bounds for 40 digits
  { 45, 10600, 11000000, 3317 },       // ECM bounds for 45 digits
  { 50, 19300, 43000000, 6558 },       // ECM bounds for 50 digits
  { 55, 49000, 110000000, 10489 },     // ECM bounds for 55 digits
};

/* ECM limits for 30, 35, ..., 95 digits */
static int limits[] = { 10, 10, 10, 10, 10, 15, 22, 26, 60, 130, 200, 270, 350 };


#define ADD 6  /* number of multiplications in an addition */
#define DUP 5  /* number of multiplications in a duplicate */

/* returns the number of modular multiplications */
static int lucas_cost(int multiplier, double v)
{
  int nbrMultiplications;
  int e;
  int d = multiplier;
  double dr = ((double)d / v) + 0.5;
  int r = (int)dr;
  if (r >= multiplier)
  {
    return (ADD * multiplier);
  }
  d = multiplier - r;
  e = (2 * r) - multiplier;
  nbrMultiplications = DUP + ADD; /* initial duplicate and final addition */
  while (d != e)
  {
    if (d < e)
    {
      r = d;
      d = e;
      e = r;
    }
    if (((4 * d) <= (5 * e)) && (((d + e) % 3) == 0))
    { /* condition 1 */
      r = ((2 * d) - e) / 3;
      e = ((2 * e) - d) / 3;
      d = r;
      nbrMultiplications += 3 * ADD; /* 3 additions */
    }
    else if (((4 * d) <= (5 * e)) && (((d - e) % 6) == 0))
    { /* condition 2 */
      d = (d - e) / 2;
      nbrMultiplications += ADD + DUP; /* one addition, one duplicate */
    }
    else if (d <= (4 * e))
    { /* condition 3 */
      d -= e;
      nbrMultiplications += ADD; /* one addition */
    }
    else if (((d + e) % 2) == 0)
    { /* condition 4 */
      d = (d - e) / 2;
      nbrMultiplications += ADD + DUP; /* one addition, one duplicate */
    }
    else if ((d % 2) == 0)
    { /* condition 5 */
      d /= 2;
      nbrMultiplications += ADD + DUP; /* one addition, one duplicate */
    }
    else if ((d % 3) == 0)
    { /* condition 6 */
      d = (d / 3) - e;
      nbrMultiplications += (3 * ADD) + DUP; /* three additions, one duplicate */
    }
    else if (((d + e) % 3) == 0)
    { /* condition 7 */
      d = (d - (2 * e)) / 3;
      nbrMultiplications += (3 * ADD) + DUP; /* three additions, one duplicate */
    }
    else if (((d - e) % 3) == 0)
    { /* condition 8 */
      d = (d - e) / 3;
      nbrMultiplications += (3 * ADD) + DUP; /* three additions, one duplicate */
    }
    else if ((e % 2) == 0)
    { /* condition 9 */
      e /= 2;
      nbrMultiplications += ADD + DUP; /* one addition, one duplicate */
    }
    else
    { /* no more conditions */
    }
  }
  return nbrMultiplications;
}

/* computes nP from P=(x:z) and puts the result in (x:z). Assumes n>2. */
static void prac(int multiplier, limb* x, limb *z)
{
  int d;
  int e;
  int r;
  int i;
  double dr;
  limb* t;
  limb* xA = x;
  limb* zA = z;
  limb* xB = common.ecm.Aux1;
  limb* zB = common.ecm.Aux2;
  limb* xC = common.ecm.Aux3;
  limb* zC = common.ecm.Aux4;
  limb* xT = common.ecm.W1;
  limb* zT = common.ecm.W2;
  limb* xT2 = common.ecm.W3;
  limb* zT2 = common.ecm.W4;
  const double v[] =
  {
    1.61803398875,
    1.72360679775,
    1.618347119656,
    1.617914406529,
    1.612429949509,
    1.632839806089,
    1.620181980807,
    1.580178728295,
    1.617214616534,
    1.38196601125 };

  /* chooses the best value of v */
  r = lucas_cost(multiplier, v[0]);
  i = 0;
  for (d = 1; d < 10; d++)
  {
    e = lucas_cost(multiplier, v[d]);
    if (e < r)
    {
      r = e;
      i = d;
    }
  }
  d = multiplier;
  dr = ((double)d / v[i]) + 0.5;
  r = (int)dr;
  /* first iteration always begins by Condition 3, then a swap */
  d = multiplier - r;
  e = (2 * r) - multiplier;
  (void)memcpy(xB, xA, NumberSizeBytes);   // B <- A
  (void)memcpy(zB, zA, NumberSizeBytes);
  (void)memcpy(xC, xA, NumberSizeBytes);   // C <- A
  (void)memcpy(zC, zA, NumberSizeBytes);
  duplicate(xA, zA, xA, zA); /* A=2*A */
  while (d != e)
  {
    if (d < e)
    {
      r = d; d = e; e = r;
      t = xA; xA = xB; xB = t;
      t = zA; zA = zB; zB = t;
    }
    /* do the first line of Table 4 whose condition qualifies */
    if (((4 * d) <= (5 * e)) && (((d + e) % 3) == 0))
    { /* condition 1 */
      r = ((2 * d) - e) / 3;
      e = ((2 * e) - d) / 3;
      d = r;
      add3(xT, zT, xA, zA, xB, zB, xC, zC); /* T = f(A,B,C) */
      add3(xT2, zT2, xT, zT, xA, zA, xB, zB); /* T2 = f(T,A,B) */
      add3(xB, zB, xB, zB, xT, zT, xA, zA); /* B = f(B,T,A) */
      t = xA; xA = xT2; xT2 = t;
      t = zA; zA = zT2; zT2 = t; /* swap A and T2 */
    }
    else if (((4 * d) <= (5 * e)) && (((d - e) % 6) == 0))
    { /* condition 2 */
      d = (d - e) / 2;
      add3(xB, zB, xA, zA, xB, zB, xC, zC); /* B = f(A,B,C) */
      duplicate(xA, zA, xA, zA); /* A = 2*A */
    }
    else if (d <= (4 * e))
    { /* condition 3 */
      d -= e;
      add3(xT, zT, xB, zB, xA, zA, xC, zC); /* T = f(B,A,C) */
      t = xB; xB = xT; xT = xC; xC = t;
      t = zB; zB = zT; zT = zC; zC = t; /* circular permutation (B,T,C) */
    }
    else if (((d + e) % 2) == 0)
    { /* condition 4 */
      d = (d - e) / 2;
      add3(xB, zB, xB, zB, xA, zA, xC, zC); /* B = f(B,A,C) */
      duplicate(xA, zA, xA, zA); /* A = 2*A */
    }
    else if ((d % 2) == 0)
    { /* condition 5 */
      d /= 2;
      add3(xC, zC, xC, zC, xA, zA, xB, zB); /* C = f(C,A,B) */
      duplicate(xA, zA, xA, zA); /* A = 2*A */
    }
    else if ((d % 3) == 0)
    { /* condition 6 */
      d = (d / 3) - e;
      duplicate(xT, zT, xA, zA); /* T1 = 2*A */
      add3(xT2, zT2, xA, zA, xB, zB, xC, zC); /* T2 = f(A,B,C) */
      add3(xA, zA, xT, zT, xA, zA, xA, zA); /* A = f(T1,A,A) */
      add3(xT, zT, xT, zT, xT2, zT2, xC, zC); /* T1 = f(T1,T2,C) */
      t = xC; xC = xB; xB = xT; xT = t;
      t = zC; zC = zB; zB = zT; zT = t; /* circular permutation (C,B,T) */
    }
    else if (((d + e) % 3) == 0)
    { /* condition 7 */
      d = (d - (2 * e)) / 3;
      add3(xT, zT, xA, zA, xB, zB, xC, zC); /* T1 = f(A,B,C) */
      add3(xB, zB, xT, zT, xA, zA, xB, zB); /* B = f(T1,A,B) */
      duplicate(xT, zT, xA, zA);
      add3(xA, zA, xA, zA, xT, zT, xA, zA); /* A = 3*A */
    }
    else if (((d - e) % 3) == 0)
    { /* condition 8 */
      d = (d - e) / 3;
      add3(xT, zT, xA, zA, xB, zB, xC, zC); /* T1 = f(A,B,C) */
      add3(xC, zC, xC, zC, xA, zA, xB, zB); /* C = f(A,C,B) */
      t = xB; xB = xT; xT = t;
      t = zB; zB = zT; zT = t; /* swap B and T */
      duplicate(xT, zT, xA, zA);
      add3(xA, zA, xA, zA, xT, zT, xA, zA); /* A = 3*A */
    }
    else if ((e % 2) == 0)
    { /* condition 9 */
      e /= 2;
      add3(xC, zC, xC, zC, xB, zB, xA, zA); /* C = f(C,B,A) */
      duplicate(xB, zB, xB, zB); /* B = 2*B */
    }
    else
    {
      /* no more conditions */
    }
  }
  add3(x, z, xA, zA, xB, zB, xC, zC);
}

// Compute the point (sum_x:sum_z) <- (Q_x:Q_z) + (R_x:R_z)
// where (diff_x:diff_z) equals (Q_x:Q_z) - (R_x:R_z)
void add3(limb* sum_x, limb* sum_z, const limb* Q_x, const limb* Q_z,
  const limb* R_x, const limb* R_z, const limb* diff_x, const limb* diff_z)
{
  limb* t = common.ecm.Aux5;
  limb* u = common.ecm.Aux6;
  limb* v = common.ecm.Aux7;
  limb* w = common.ecm.Aux8;
  SubtBigNbrModN(Q_x, Q_z, v, TestNbr, NumberLength);  // v = x2-z2
  AddBigNbrModN(R_x, R_z, w, TestNbr, NumberLength);   // w = x1+z1
  modmult(v, w, u);                                    // u = (x2-z2)*(x1+z1)
  AddBigNbrModN(Q_x, Q_z, w, TestNbr, NumberLength);   // w = x2+z2
  SubtBigNbrModN(R_x, R_z, t, TestNbr, NumberLength);  // t = x1-z1
  modmult(t, w, v);                                    // v = (x2+z2)*(x1-z1)
  AddBigNbrModN(u, v, t, TestNbr, NumberLength);       // t = 2*(x1*x2-z1*z2)
  modmult(t, t, w);                                    // w = 4*(x1*x2-z1*z2)^2
  SubtBigNbrModN(u, v, t, TestNbr, NumberLength);      // t = 2*(x2*z1-x1*z2)
  modmult(t, t, v);                                    // v = 4*(x2*z1-x1*z2)^2
  if (!memcmp(diff_x, sum_x, NumberSizeBytes))
  {
    (void)memcpy(u, diff_x, NumberSizeBytes);
    (void)memcpy(t, w, NumberSizeBytes);
    modmult(diff_z, t, w);
    modmult(v, u, sum_z);
    (void)memcpy(sum_x, w, NumberSizeBytes);
  }
  else
  {
    modmult(w, diff_z, sum_x);                    // sum_x = 4*z*(x1*x2-z1*z2)^2
    modmult(diff_x, v, sum_z);                    // sum_z = 4*x*(x2*z1-x1*z2)^2
  }
}

/* computes 2P=(x2:z2) from P=(x1:z1), with 5 mul, 4 add/sub, 5 mod.
Uses the following global variables:
- n : number to factor
- b : (a+2)/4 mod n
- u, v, w : auxiliary variables
Modifies: x2, z2, u, v, w
*/
void duplicate(limb* dupl_x, limb* dupl_z, const limb* P_x, const limb* P_z)
{
  limb* u = common.ecm.Aux5;
  limb* v = common.ecm.Aux6;
  limb* w = common.ecm.Aux7;
  AddBigNbrModN(P_x, P_z, w, TestNbr, NumberLength);   // w = x1+z1
  modmult(w, w, u);                                    // u = (x1+z1)^2
  SubtBigNbrModN(P_x, P_z, w, TestNbr, NumberLength);  // w = x1-z1
  modmult(w, w, v);                                    // v = (x1-z1)^2
  modmult(u, v, dupl_x);                               // x2 = u*v = (x1^2 - z1^2)^2
  SubtBigNbrModN(u, v, w, TestNbr, NumberLength);      // w = u-v = 4*x1*z1
  modmult(common.ecm.AA, w, u);
  AddBigNbrModN(u, v, u, TestNbr, NumberLength);       // u = (v+b*w)
  modmult(w, u, dupl_z);                               // z2 = (w*u)
}

int gcdIsOne(const limb* value)
{
  UncompressLimbsBigInteger(value, &common.ecm.Temp1);
  UncompressLimbsBigInteger(TestNbr, &common.ecm.Temp2);
  // Return zero if value is zero or both numbers are equal.
  if (BigIntIsZero(&common.ecm.Temp1))
  {
    return 0;
  }
  BigIntSubt(&common.ecm.Temp1, &common.ecm.Temp2, &common.ecm.Temp3);
  if (BigIntIsZero(&common.ecm.Temp3))
  {
    return 0;
  }
  BigIntGcd(&common.ecm.Temp1, &common.ecm.Temp2, &common.ecm.Temp3);
  CompressLimbsBigInteger(common.ecm.GD, &common.ecm.Temp3);
  if ((common.ecm.Temp3.nbrLimbs == 1) && (common.ecm.Temp3.limbs[0].x < 2))
  {
    return common.ecm.Temp3.limbs[0].x;    // GCD is less than 2.
  }
  return 2;      // GCD is greater than one.
}

// This routine requires the array isCoprime210_2310 to be initialized with 1 for numbers
// that are coprime with 2310 and 0 otherwise.
void GenerateSieve(int initial, int indexFirstPrime)
{
  int indexPrime;
  int i;
  int prime;
  for (int startBlockOffset = 0; startBlockOffset < MAX_SIEVE_PRIME; startBlockOffset += SIEVE_SIZE)
  {
    (void)memcpy(&common.ecm.sieve[startBlockOffset], common.ecm.isCoprime210_2310, SIEVE_SIZE);
  }
  indexPrime = indexFirstPrime;
  prime = SmallPrime[indexPrime];
  do
  {
    if (initial > (prime * prime))
    {
      int initModPrime = initial % prime;
      if ((initModPrime & 1) != 0)
      {    // initModPrime is odd
        i = (prime - initModPrime) / 2;
      }
      else if (initModPrime == 0)
      {
        i = 0;
      }
      else
      {    // initModPrime is even
        i = prime - (initModPrime / 2);
      }
      for (; i < MAX_SIEVE_PRIME; i += prime)
      {
        common.ecm.sieve[i] = SIEVE_COMPOSITE; /* Composite */
      }
    }
    else
    {
      i = (prime * prime) - initial;
      if (i < (2 * MAX_SIEVE_PRIME))
      {
        for (i = i / 2; i < MAX_SIEVE_PRIME; i += prime)
        {
          common.ecm.sieve[i] = SIEVE_COMPOSITE; /* Composite */
        }
      }
      else
      {
        break;
      }
    }
    indexPrime++;
    prime = SmallPrime[indexPrime];
  } while (prime < MAX_SIEVE_PRIME);
}

/*******************************/
/* First step of ECM algorithm */
/*******************************/
static enum eEcmResult ecmStep1(void)
{
  int prime;
  int sieveIndex;
  int sievedNumber;
  int retcode;
  int bufSize = (NumberLength + 1) * (int)sizeof(limb);
  (void)memcpy(common.ecm.Xbak, common.ecm.X, bufSize);
  (void)memcpy(common.ecm.Zbak, common.ecm.Z, bufSize);
  (void)memcpy(common.ecm.GcdAccumulated, MontgomeryMultR1, bufSize);
  StepECM = 1;
  for (int pass = 0; pass < 2; pass++)
  {
    /* For powers of 2 and 3 */
    for (prime = 2; prime <= 3; prime++)
    {
      for (int powerPrime = prime; powerPrime <= boundStep1; powerPrime *= prime)
      {
        if (prime == 2)
        {    // Multiply point by 2.
          duplicate(common.ecm.X, common.ecm.Z, common.ecm.X, common.ecm.Z);
        }
        else if (prime == 3)
        {    // Multiply point by 3.
          duplicate(common.ecm.W1, common.ecm.W2, common.ecm.X, common.ecm.Z);
          add3(common.ecm.X, common.ecm.Z, common.ecm.X, common.ecm.Z,
            common.ecm.W1, common.ecm.W2, common.ecm.X, common.ecm.Z);
        }
        if (pass == 1)
        {
          retcode = gcdIsOne(common.ecm.Z);
          if (retcode == 0)
          {
            return FACTOR_NOT_FOUND_GCD;
          }
          if (retcode > 1)
          {
            return FACTOR_FOUND;
          }
        }
      }
    }
    /* for powers of odd primes */

    indexM = 1;
    do
    {
      prime = SmallPrime[indexM];
#ifdef __EMSCRIPTEN__
      intPrime = prime;
#endif
      for (int64_t largePrimePower = prime; largePrimePower <= boundStep1;
        largePrimePower *= prime)
      {
        prac(prime, common.ecm.X, common.ecm.Z);  // Multiply (X:Z) by prime.
      }
      indexM++;
      if (pass == 1)
      {
        retcode = gcdIsOne(common.ecm.Z);
        if (retcode == 0)
        {
          return FACTOR_NOT_FOUND_GCD;
        }
        if (retcode > 1)
        {
          return FACTOR_FOUND;
        }
      }
    } while (SmallPrime[indexM - 1] <= sqrtBoundStep1);
    int startSieve = prime + 2;

    /* Initialize isCoprime210_2310[n]: 1 if gcd(P+2n,2310) == 1, 0 otherwise */
    sievedNumber = startSieve;
    for (sieveIndex = 0; sieveIndex < SIEVE_SIZE; sieveIndex++)
    {
      common.ecm.isCoprime210_2310[sieveIndex] =
        ((((sievedNumber % 3) == 0)
          || ((sievedNumber % 5) == 0)
          || ((sievedNumber % 7) == 0)
          || ((sievedNumber % 11) == 0)
          ) ? (unsigned char)NOT_COPRIME_2310 : (unsigned char)COPRIME_2310);
      sievedNumber += 2;
    }    
    do
    {
      int currentPrime = startSieve;
      /* Generate sieve */
      GenerateSieve(startSieve, 5);

      /* Walk through sieve */
      for (sieveIndex = 0; sieveIndex < MAX_SIEVE_PRIME; sieveIndex++)
      {
        if (common.ecm.sieve[sieveIndex] == SIEVE_COMPOSITE)
        {
          currentPrime += 2;
          continue; /* Do not process composites */
        }
        if (currentPrime > boundStep1)
        {
          break;
        }
#ifdef __EMSCRIPTEN__
        intPrime = currentPrime;
#endif
        prac(currentPrime, common.ecm.X, common.ecm.Z);  // Multiply (X:Z) by prime.
        if (pass == 1)
        {
          retcode = gcdIsOne(common.ecm.Z);
          if (retcode == 0)
          {
            return FACTOR_NOT_FOUND_GCD;
          }
          if (retcode > 1)
          {
            return FACTOR_FOUND;
          }
        }
        currentPrime += 2;
      }
      startSieve = currentPrime;
    } while (startSieve < boundStep1);
    if (pass == 0)
    {
      int result = gcdIsOne(common.ecm.Z);
      if (result == 1)
      {
        break;         // GCD is 1 so factor is not found.
      }
      if (result == 2)
      {                // GCD greater than zero. Factor found.
        return FACTOR_FOUND;
      }
      // Factor could not be found because GCD is zero.
      // Repeat this curve but performing GCD on every step.
      (void)memcpy(common.ecm.X, common.ecm.Xbak, bufSize);
      (void)memcpy(common.ecm.Z, common.ecm.Zbak, bufSize);
    }
  } /* end for Pass */
  return FACTOR_NOT_FOUND;
}

/******************************************************/
/* Second step (using improved standard continuation) */
/******************************************************/
static enum eEcmResult ecmStep2(void)
{
  int coprimeIndex = 0;
  int index = 0;
  int indexFirstPrime;
  StepECM = 2;
  if (GROUP_SIZE * NumberLength > sizeof(common.ecm.root) / sizeof(common.ecm.root[0]))
  {  // Numbers are big. Use short group.
    groupSize = SHORT_GROUP_SIZE;
    sieveSize = SHORT_SIEVE_SIZE;
    halfSieveSize = SHORT_SIEVE_SIZE / 2;
    indexFirstPrime = 4;
  }
  else
  {  // Numbers are small. Use long group.
    groupSize = GROUP_SIZE;
    sieveSize = SIEVE_SIZE;
    halfSieveSize = SIEVE_SIZE / 2;
    indexFirstPrime = 5;
  }
  for (int sievedNumber = 1; sievedNumber < sieveSize; sievedNumber += 2)
  {
    if (((sievedNumber % 3) == 0) || ((sievedNumber % 5) == 0) || ((sievedNumber % 7) == 0)
       || ((indexFirstPrime == 5) && ((sievedNumber % 11) == 0)))
    {
      common.ecm.isCoprime210_2310[index] = NOT_COPRIME_2310;
    }
    else
    {
      common.ecm.sieveidx[coprimeIndex] = index;
      common.ecm.isCoprime210_2310[index] = COPRIME_2310;
      coprimeIndex++;
    }
    index++;
  }
  // At this moment half of the isCoprime210_2310 array is filled.
  // Fill the other half with the same values because gcd(a, 2310) = gcd(a + 2310, 2310).
  memcpy(&common.ecm.isCoprime210_2310[halfSieveSize], common.ecm.isCoprime210_2310, halfSieveSize);
  for (int pass = 0; pass < 2; pass++)
  {
    limb* ptrRoot = common.ecm.root;
    int firstIndexM;
    (void)memcpy(common.ecm.Xbak, common.ecm.X, NumberSizeBytes);  // (X:Z) -> Q (output
    (void)memcpy(common.ecm.Zbak, common.ecm.Z, NumberSizeBytes);  //         from step 1)
    (void)memcpy(common.ecm.GcdAccumulated, MontgomeryMultR1, NumberSizeBytes);
    (void)memcpy(common.ecm.UX, common.ecm.X, NumberSizeBytes);
    (void)memcpy(common.ecm.UZ, common.ecm.Z, NumberSizeBytes);    // (UX:UZ) <- Q 
    (void)ModInvBigNbr(common.ecm.Z, common.ecm.Aux1, TestNbr, NumberLength);
    modmult(common.ecm.Aux1, common.ecm.X, ptrRoot);               // root[0] <- X/Z (Q)
    ptrRoot += NumberLength;
    duplicate(common.ecm.TX, common.ecm.TZ, common.ecm.X, common.ecm.Z);  // (TX:TZ) <- 2Q
    // Compute 3Q, 5Q, 7Q, 11Q, ... up to (sieveSize-1)*Q.
    // Store kQ where k is coprime with 2310 in root[].
    for (int sieveIndex = 3; sieveIndex < sieveSize; sieveIndex += 2)
    {
      // At this moment (X:Z) = (k-2)Q, (UX:UZ) = (k-4)Q and (TX:TZ) = 2Q.
      // The invariant is (X:Z) = (UX:UZ) + (TX:TZ) = (k-2)Q.
      add3(common.ecm.WX, common.ecm.WZ, common.ecm.X, common.ecm.Z,
        common.ecm.TX, common.ecm.TZ, common.ecm.UX, common.ecm.UZ); // (WX:WZ) <- kQ
      (void)memcpy(common.ecm.UX, common.ecm.X, NumberSizeBytes);    // (UX:UZ) <- (k-2)Q
      (void)memcpy(common.ecm.UZ, common.ecm.Z, NumberSizeBytes);
      (void)memcpy(common.ecm.X, common.ecm.WX, NumberSizeBytes);    // (X:Z) <- kQ
      (void)memcpy(common.ecm.Z, common.ecm.WZ, NumberSizeBytes);
      if (sieveIndex == halfSieveSize)
      {      // halfSieveSize is odd.
        (void)memcpy(common.ecm.W3, common.ecm.X, NumberSizeBytes);
        (void)memcpy(common.ecm.W4, common.ecm.Z, NumberSizeBytes);  // (W3:W4) <- (sieveSize/2)*Q
      }
      if (common.ecm.isCoprime210_2310[(sieveIndex-1)/2] == COPRIME_2310)
      {    // sieveIndex is coprime with 2310. Store X/Z in root[].
        (void)ModInvBigNbr(common.ecm.Z, common.ecm.Aux1, TestNbr, NumberLength);
        modmult(common.ecm.Aux1, common.ecm.X, ptrRoot);             // root[J] <- X/Z
        ptrRoot += NumberLength;
      }
    } /* end for sieveIndex */
    assert(ptrRoot - &common.ecm.root[0] == groupSize * NumberLength);
    duplicate(common.ecm.UX, common.ecm.UZ, common.ecm.W3, common.ecm.W4);  // (UX:UZ) <- sieveSize*Q
    duplicate(common.ecm.TX, common.ecm.TZ, common.ecm.UX, common.ecm.UZ);  // (TX:TZ) <- 2*sieveSize*Q
    (void)memcpy(common.ecm.X, common.ecm.UX, NumberSizeBytes);             // (X:Z) <- sieveSize*Q
    (void)memcpy(common.ecm.Z, common.ecm.UZ, NumberSizeBytes);
    firstIndexM = boundStep1 / (2 * sieveSize);
    maxIndexM = (int)(boundStep2 / (2 * sieveSize));
#ifdef __EMSCRIPTEN__
    longPrime = 1;
#endif
    for (indexM = 0; indexM <= maxIndexM; indexM++)
    {
      if (indexM >= firstIndexM)
      { // If inside step 2 range... 
        bool rc = ModInvBigNbr(common.ecm.Z, common.ecm.Aux3, TestNbr, NumberLength);
        if (rc == false)
        {
          (void)memcpy(common.ecm.GD, common.ecm.Z, NumberSizeBytes);
          return FACTOR_FOUND;
        }
        // Compute Aux as X/Z for m*sieveSize*Q.
        modmult(common.ecm.X, common.ecm.Aux3, common.ecm.Aux1);

          /* Generate sieve */
        if (((indexM % 10) == 0) || (indexM == firstIndexM))
        {  // Generate sieve for next 10 blocks of sieveSize numbers or the first time the step 2 executes.
          GenerateSieve((indexM / 10) * (20 * sieveSize) + 1, indexFirstPrime);
        }
        /* Walk through sieve */
        int startSieveBlock = halfSieveSize + (indexM % 10) * sieveSize;
        ptrRoot = common.ecm.root;
        for (int groupIndex = 0; groupIndex < groupSize; groupIndex++)
        {
          int delta = common.ecm.sieveidx[groupIndex]; // 0 < delta < sieveSize
          if ((common.ecm.sieve[startSieveBlock + delta] == SIEVE_PROBABLE_PRIME) ||
            (common.ecm.sieve[startSieveBlock - 1 - delta] == SIEVE_PROBABLE_PRIME))
          {    // At least one of the two numbers is probable prime. Compute the GCD.
            SubtBigNbrModN(common.ecm.Aux1, ptrRoot, common.ecm.Aux2, TestNbr, NumberLength);
            modmult(common.ecm.GcdAccumulated, common.ecm.Aux2, common.ecm.GcdAccumulated);
          }
          ptrRoot += NumberLength;
        }
        if (pass != 0)
        {
          if (BigNbrIsZero(common.ecm.GcdAccumulated))
          {           // This curve cannot factor the number.
            return FACTOR_NOT_FOUND_GCD;
          }
          if (gcdIsOne(common.ecm.GcdAccumulated) > 1)
          {
            return FACTOR_FOUND;
          }
        }
      }   // End if inside step 2 range.
      // At this moment (X:Z) = (2m+1)*sieveSize*Q, (UX:UZ) = (2m-1)*sieveSize*Q and (TX:TZ) = 2*sieveSize*Q.
      // The invariant is (X:Z) = (UX:UZ) + (TX:TZ) = (2m+1)*sieveSize*Q.
      add3(common.ecm.WX, common.ecm.WZ, common.ecm.X, common.ecm.Z,
           common.ecm.TX, common.ecm.TZ, common.ecm.UX, common.ecm.UZ); // (WX:WZ) <- (2m+3)*sieveSize*Q
      (void)memcpy(common.ecm.UX, common.ecm.X, NumberSizeBytes);       // (UX:UZ) <- (2m+1)*sieveSize*Q
      (void)memcpy(common.ecm.UZ, common.ecm.Z, NumberSizeBytes);
      (void)memcpy(common.ecm.X, common.ecm.WX, NumberSizeBytes);       // (X:Z) <- (2m+3)*sieveSize*Q
      (void)memcpy(common.ecm.Z, common.ecm.WZ, NumberSizeBytes);
#ifdef __EMSCRIPTEN__
      longPrime += 2 * sieveSize;
#endif
    } // end for indexM
    if (pass == 0)
    {
      int rc;
      if (BigNbrIsZero(common.ecm.GcdAccumulated))
      { // If GcdAccumulated is zero
        (void)memcpy(common.ecm.X, common.ecm.Xbak, NumberSizeBytes);
        (void)memcpy(common.ecm.Z, common.ecm.Zbak, NumberSizeBytes);
        continue; // multiple of TestNbr, continue.
      }
      rc = gcdIsOne(common.ecm.GcdAccumulated);
      if (rc == 1)
      {
        break;    // GCD is one, so this curve does not find a factor.
      }
      if (rc == 0)
      {
        continue;
      }
      // Compute GD as GCD(GcdAccumulated, TestNbr)
      if (memcmp(common.ecm.GD, TestNbr, NumberSizeBytes) != 0)
      {           // GCD is not 1 or TestNbr
        return FACTOR_FOUND;
      }
    }
  } /* end for Pass */
  return FACTOR_NOT_FOUND;
}

static void initSmallPrimeArray(void)
{
  int potentialPrime = 3;
  int numPrimes = (int)sizeof(SmallPrime) / (int)sizeof(SmallPrime[0]);
  SmallPrime[0] = 2;
  for (indexM = 1; indexM < numPrimes; indexM++)
  {     // Loop that fills the SmallPrime array.
    int divisor;
    SmallPrime[indexM] = potentialPrime; /* Store prime */
    do
    {
      potentialPrime += 2;
      for (divisor = 3; (divisor * divisor) <= potentialPrime; divisor += 2)
      { /* Check if potentialPrime is prime */
        if ((potentialPrime % divisor) == 0)
        {
          break;  /* Composite */
        }
      }
    } while ((divisor * divisor) <= potentialPrime);
  }
}

enum eEcmResult ecmCurve(int *pEC, int *pNextEC)
{
  int EC = *pEC;
  int NextEC = *pNextEC;
  enum eEcmResult result;
  NumberSizeBytes = NumberLength * (int)sizeof(limb);
#ifdef __EMSCRIPTEN__
  char text[20];
#endif
  if (SmallPrime[0] != 2)
  {    // Not initialized yet.
    initSmallPrimeArray();
  }
  EC %= 50000000;   // Convert to curve number.
  do
  {
#ifdef __EMSCRIPTEN__
    char* ptrText;
#endif
    int nbrDigits;
    if (NextEC > 0)
    {
      EC = NextEC;
      NextEC = -1;
      if (EC >= TYP_SIQS)
      {
        int bufSize = (NumberLength - 1) * (int)sizeof(limb);
        common.ecm.GD[0].x = 1;   // Set GD to 1.
        (void)memset(&common.ecm.GD[1], 0, bufSize);
        *pEC = EC;
        *pNextEC = NextEC;
        return FACTOR_FOUND;
      }
    }
    else
    {
      EC++;
#ifdef __EMSCRIPTEN__
      ptrText = text;
      formatString(&ptrText, "7$1d", EC);
      databack(text);
#endif
      nbrDigits = NumberLength * 9;          // Get number of digits.
      if ((NextEC == 0) && (nbrDigits >= 30) && (nbrDigits <= 110))
      {          // Force switch to SIQS and number not too large.
        EC += TYP_SIQS;
        *pEC = EC;
        *pNextEC = NextEC;
        (void)memcpy(common.ecm.GD, TestNbr, NumberSizeBytes);
        return CHANGE_TO_SIQS;
      }
      if ((nbrDigits > 30) && (nbrDigits <= 90))  // If between 30 and 90 digits...         
      {                             // Switch to SIQS.
        int limit = limits[(nbrDigits - 31) / 5];
        if ((EC % 50000000) >= limit)
        {                           // Switch to SIQS.
          EC += TYP_SIQS;
          *pEC = EC;
          *pNextEC = NextEC;
          (void)memcpy(common.ecm.GD, TestNbr, NumberSizeBytes);
          return CHANGE_TO_SIQS;
        }
      }
    }
    // Compute bounds according to the curve number.
    int curveNbr = EC;
    const struct sBounds* pstBounds = &stEcmBounds[0];
    do
    {
      if (curveNbr < pstBounds->nbrCurves)
      {
        break;
      }
      curveNbr -= pstBounds->nbrCurves;
      pstBounds++;
    } while (pstBounds->boundStep1 != 110000000);
    boundStep1 = pstBounds->boundStep1;
    boundStep2 = (int64_t)boundStep1 * 100;
    sqrtBoundStep1 = pstBounds->sqrtBoundStep1;
#ifdef __EMSCRIPTEN__
    intStep1Bound = boundStep1;
    longStep2Bound = boundStep2;
    ptrText = ptrLowerText;  // Point after number that is being factored.
    copyStr(&ptrText, "<p>");
    formatString(&ptrText, LITERAL_ECM1, pstBounds->digitLevel);
    formatString(&ptrText, " <meter min=\"0\" max=\"$1d\" value=\"$2d\">$3d%</meter></p><p>",
          pstBounds->nbrCurves, curveNbr, curveNbr * 100 / pstBounds->nbrCurves);
    formatString(&ptrText, LITERAL_ECM2, EC, boundStep1, boundStep2);
    copyStr(&ptrText, "</p>");
    databack(lowerText);
#endif

    //  Compute W1 <- 2 * (EC+1)*modinv(3 * (EC+1) ^ 2 - 1, N) mod N
                                               // Aux2 <- 1 in Montgomery notation.
    (void)memcpy(common.ecm.Aux2, MontgomeryMultR1, NumberSizeBytes);
    // Compute Aux2 as EC + 1.
    modmultInt(common.ecm.Aux2, EC + 1, common.ecm.Aux2);
    // Compute Aux1 as 2*(EC+1).
    modmultInt(common.ecm.Aux2, 2, common.ecm.Aux1);
    // Compute Aux3 as (EC + 1)^2.
    modmultInt(common.ecm.Aux2, EC + 1, common.ecm.Aux3);
    // Compute Aux3 as 3*(EC + 1)^2.
    modmultInt(common.ecm.Aux3, 3, common.ecm.Aux3);
    // Compute Aux2 as 3*(EC + 1)^2 - 1.
    SubtBigNbrModN(common.ecm.Aux3, MontgomeryMultR1, common.ecm.Aux2, TestNbr, NumberLength);
    (void)ModInvBigNbr(common.ecm.Aux2, common.ecm.Aux2, TestNbr, NumberLength);
    // Compute W1 as 2*(EC+1)/(3*(EC+1)^2 - 1)
    modmult(common.ecm.Aux1, common.ecm.Aux2, common.ecm.W1);

    //  if W1*(W1 ^ 2 - 1)*(9 * W1 ^ 2 - 1) mod N=0 then select another curve.
    modmult(common.ecm.W1, common.ecm.W1, common.ecm.W2);          // W2 <- W1^2
    modmult(common.ecm.W2, common.ecm.W1, common.ecm.W3);         // W3 <- W1^3
    SubtBigNbrModN(common.ecm.W3, common.ecm.W1, common.ecm.Aux1, TestNbr, NumberLength);  // Aux1 <- W1^3 - W1
    modmultInt(common.ecm.W2, 9, common.ecm.Aux2);      // Aux2 <- 9*W1^2
    SubtBigNbrModN(common.ecm.Aux2, MontgomeryMultR1, common.ecm.Aux2, TestNbr, NumberLength); // Aux2 <- 9*W1^2-1
    modmult(common.ecm.Aux1, common.ecm.Aux2, common.ecm.Aux3);
  } while (BigNbrIsZero(common.ecm.Aux3));
  //   Compute Z as 4 * W1 mod N
  modmultInt(common.ecm.W1, 4, common.ecm.Z);
  //   Compute A as(-3 * W1 ^ 4 - 6 * W1 ^ 2 + 1)*modinv(4 * W1 ^ 3, N) mod N
  modmultInt(common.ecm.W2, 6, common.ecm.Aux1);      // Aux1 <- 6*W1^2
  SubtBigNbrModN(MontgomeryMultR1, common.ecm.Aux1, common.ecm.Aux1, TestNbr, NumberLength);
  modmult(common.ecm.W2, common.ecm.W2, common.ecm.Aux2);       // Aux2 <- W1^4
  modmultInt(common.ecm.Aux2, 3, common.ecm.Aux2);     // Aux2 <- 3*W1^4
  SubtBigNbrModN(common.ecm.Aux1, common.ecm.Aux2, common.ecm.Aux1, TestNbr, NumberLength);
  modmultInt(common.ecm.W3, 4, common.ecm.Aux2);      // Aux2 <- 4*W1^3
  (void)ModInvBigNbr(common.ecm.Aux2, common.ecm.Aux3, TestNbr, NumberLength);
  modmult(common.ecm.Aux1, common.ecm.Aux3, common.ecm.W1);
  //   Compute AA as (A + 2)*modinv(4, N) mod N
  modmultInt(MontgomeryMultR1, 2, common.ecm.Aux2);  // Aux2 <- 2
  AddBigNbrModN(common.ecm.W1, common.ecm.Aux2, common.ecm.Aux1, TestNbr, NumberLength); // Aux1 <- W1+2
  modmultInt(MontgomeryMultR1, 4, common.ecm.Aux2);  // Aux2 <- 4
  (void)ModInvBigNbr(common.ecm.Aux2, common.ecm.Aux2, TestNbr, NumberLength);
  modmult(common.ecm.Aux1, common.ecm.Aux2, common.ecm.AA);
  //   Compute X as (3 * W1 ^ 2 + 1) mod N
  modmultInt(common.ecm.W2, 3, common.ecm.Aux1);    // Aux1 <- 3*W1^2
  AddBigNbrModN(common.ecm.Aux1, MontgomeryMultR1, common.ecm.X, TestNbr, NumberLength);
  result = ecmStep1();
  if (result == FACTOR_FOUND)
  {
    *pEC = EC;
    *pNextEC = NextEC;
    return result;
  }
  if (result != FACTOR_NOT_FOUND_GCD)
  {
    result = ecmStep2();
  }
  if (result == FACTOR_FOUND)
  {
    *pEC = EC;
    *pNextEC = NextEC;
    return result;
  }
  if (result == FACTOR_NOT_FOUND_GCD)
  {
    *pEC = EC + 1;
    *pNextEC = NextEC;
    return FACTOR_NOT_FOUND_GCD;
  }
  *pEC = EC;
  *pNextEC = NextEC;
  return FACTOR_NOT_FOUND;
}

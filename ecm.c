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
#include <math.h>
#include "bignbr.h"
#include "expression.h"
#include "factor.h"
#include "commonstruc.h"
#include "skiptest.h"

#ifdef __EMSCRIPTEN__
extern int nbrPrimes;
extern int indexPrimes;
extern char* ptrLowerText;
extern char lowerText[];
#endif
extern int StepECM;
extern int maxIndexM;
extern int indexM;
static int SmallPrime[670]; /* Primes < 5000 */
static int boundStep1;
static int boundStep2;
static int sqrtBoundStep1;
static int NumberSizeBytes;

/* ECM limits for 30, 35, ..., 95 digits */
static int limits[] = { 10, 10, 10, 10, 10, 15, 22, 26, 35, 50, 100, 150, 250 };

static void add3(limb* x3, limb* z3, const limb* x2, const limb* z2,
  const limb* x1, const limb* z1, const limb* x, const limb* z);
static void duplicate(limb* x2, limb* z2, const limb* x1, const limb* z1);

/******************************************************/
/* Start of code adapted from Paul Zimmermann's ECM4C */
/******************************************************/
#define ADD 6  /* number of multiplications in an addition */
#define DUP 5  /* number of multiplications in a duplicate */

/* returns the number of modular multiplications */
static int lucas_cost(int n, double v)
{
  int c;
  int e;
  int d = n;
  double dr = ((double)d / v) + 0.5;
  int r = (int)dr;
  if (r >= n)
  {
    return (ADD * n);
  }
  d = n - r;
  e = (2 * r) - n;
  c = DUP + ADD; /* initial duplicate and final addition */
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
      c += 3 * ADD; /* 3 additions */
    }
    else if (((4 * d) <= (5 * e)) && (((d - e) % 6) == 0))
    { /* condition 2 */
      d = (d - e) / 2;
      c += ADD + DUP; /* one addition, one duplicate */
    }
    else if (d <= (4 * e))
    { /* condition 3 */
      d -= e;
      c += ADD; /* one addition */
    }
    else if (((d + e) % 2) == 0)
    { /* condition 4 */
      d = (d - e) / 2;
      c += ADD + DUP; /* one addition, one duplicate */
    }
    else if ((d % 2) == 0)
    { /* condition 5 */
      d /= 2;
      c += ADD + DUP; /* one addition, one duplicate */
    }
    else if ((d % 3) == 0)
    { /* condition 6 */
      d = (d / 3) - e;
      c += (3 * ADD) + DUP; /* three additions, one duplicate */
    }
    else if (((d + e) % 3) == 0)
    { /* condition 7 */
      d = ((d - 2) * e) / 3;
      c += (3 * ADD) + DUP; /* three additions, one duplicate */
    }
    else if (((d - e) % 3) == 0)
    { /* condition 8 */
      d = (d - e) / 3;
      c += (3 * ADD) + DUP; /* three additions, one duplicate */
    }
    else if ((e % 2) == 0)
    { /* condition 9 */
      e /= 2;
      c += ADD + DUP; /* one addition, one duplicate */
    }
  }
  return c;
}

/* computes nP from P=(x:z) and puts the result in (x:z). Assumes n>2. */
void prac(int n, limb* x, limb* z, limb* xT, limb* zT, limb* xT2, limb* zT2)
{
  int d;
  int e;
  int r;
  int i;
  limb* t;
  limb* xA = x;
  limb* zA = z;
  limb* xB = common.ecm.Aux1;
  limb* zB = common.ecm.Aux2;
  limb* xC = common.ecm.Aux3;
  limb* zC = common.ecm.Aux4;
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
  r = lucas_cost(n, v[0]);
  i = 0;
  for (d = 1; d < 10; d++)
  {
    e = lucas_cost(n, v[d]);
    if (e < r)
    {
      r = e;
      i = d;
    }
  }
  d = n;
  r = (int)(((double)d / v[i]) + 0.5);
  /* first iteration always begins by Condition 3, then a swap */
  d = n - r;
  e = (2 * r) - n;
  (void)memcpy(xB, xA, NumberSizeBytes);   // B <- A
  (void)memcpy(zB, zA, NumberSizeBytes);
  (void)memcpy(xC, xA, NumberSizeBytes);   // C <- A
  (void)memcpy(zC, zA, NumberSizeBytes);
  duplicate(xA, zA, xA, zA); /* A=2*A */
  while (d != e)
  {
    if (d < e)
    {
      r = d;
      d = e;
      e = r;
      t = xA;
      xA = xB;
      xB = t;
      t = zA;
      zA = zB;
      zB = t;
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
      t = xA;
      xA = xT2;
      xT2 = t;
      t = zA;
      zA = zT2;
      zT2 = t; /* swap A and T2 */
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
      t = xB;
      xB = xT;
      xT = xC;
      xC = t;
      t = zB;
      zB = zT;
      zT = zC;
      zC = t; /* circular permutation (B,T,C) */
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
      t = xC;
      xC = xB;
      xB = xT;
      xT = t;
      t = zC;
      zC = zB;
      zB = zT;
      zT = t; /* circular permutation (C,B,T) */
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
      t = xB;
      xB = xT;
      xT = t;
      t = zB;
      zB = zT;
      zT = t; /* swap B and T */
      duplicate(xT, zT, xA, zA);
      add3(xA, zA, xA, zA, xT, zT, xA, zA); /* A = 3*A */
    }
    else if ((e % 2) == 0)
    { /* condition 9 */
      e /= 2;
      add3(xC, zC, xC, zC, xB, zB, xA, zA); /* C = f(C,B,A) */
      duplicate(xB, zB, xB, zB); /* B = 2*B */
    }
  }
  add3(x, z, xA, zA, xB, zB, xC, zC);
}

/* adds Q=(x2:z2) and R=(x1:z1) and puts the result in (x3:z3),
using 5/6 mul, 6 add/sub and 6 mod. One assumes that Q-R=P or R-Q=P where P=(x:z).
Uses the following global variables:
- n : number to factor
- x, z : coordinates of P
- u, v, w : auxiliary variables
Modifies: x3, z3, u, v, w.
(x3,z3) may be identical to (x2,z2) and to (x,z)
*/
static void add3(limb* x3, limb* z3, const limb* x2, const limb* z2,
  const limb* x1, const limb* z1, const limb* x, const limb* z)
{
  limb* t = common.ecm.fieldTX;
  limb* u = common.ecm.fieldTZ;
  limb* v = common.ecm.fieldUX;
  limb* w = common.ecm.fieldUZ;
  SubtBigNbrModN(x2, z2, v, TestNbr, NumberLength); // v = x2-z2
  AddBigNbrModN(x1, z1, w, TestNbr, NumberLength);      // w = x1+z1
  modmult(v, w, u);       // u = (x2-z2)*(x1+z1)
  AddBigNbrModN(x2, z2, w, TestNbr, NumberLength);      // w = x2+z2
  SubtBigNbrModN(x1, z1, t, TestNbr, NumberLength); // t = x1-z1
  modmult(t, w, v);       // v = (x2+z2)*(x1-z1)
  AddBigNbrModN(u, v, t, TestNbr, NumberLength);        // t = 2*(x1*x2-z1*z2)
  modmult(t, t, w);       // w = 4*(x1*x2-z1*z2)^2
  SubtBigNbrModN(u, v, t, TestNbr, NumberLength);   // t = 2*(x2*z1-x1*z2)
  modmult(t, t, v);       // v = 4*(x2*z1-x1*z2)^2
  if (!memcmp(x, x3, NumberSizeBytes))
  {
    (void)memcpy(u, x, NumberSizeBytes);
    (void)memcpy(t, w, NumberSizeBytes);
    modmult(z, t, w);
    modmult(v, u, z3);
    (void)memcpy(x3, w, NumberSizeBytes);
  }
  else
  {
    modmult(w, z, x3); // x3 = 4*z*(x1*x2-z1*z2)^2
    modmult(x, v, z3); // z3 = 4*x*(x2*z1-x1*z2)^2
  }
}

/* computes 2P=(x2:z2) from P=(x1:z1), with 5 mul, 4 add/sub, 5 mod.
Uses the following global variables:
- n : number to factor
- b : (a+2)/4 mod n
- u, v, w : auxiliary variables
Modifies: x2, z2, u, v, w
*/
static void duplicate(limb* x2, limb* z2, const limb* x1, const limb* z1)
{
  limb* u = common.ecm.fieldUZ;
  limb* v = common.ecm.fieldTX;
  limb* w = common.ecm.fieldTZ;
  AddBigNbrModN(x1, z1, w, TestNbr, NumberLength);      // w = x1+z1
  modmult(w, w, u);       // u = (x1+z1)^2
  SubtBigNbrModN(x1, z1, w, TestNbr, NumberLength); // w = x1-z1
  modmult(w, w, v);       // v = (x1-z1)^2
  modmult(u, v, x2);      // x2 = u*v = (x1^2 - z1^2)^2
  SubtBigNbrModN(u, v, w, TestNbr, NumberLength);   // w = u-v = 4*x1*z1
  modmult(common.ecm.fieldAA, w, u);
  AddBigNbrModN(u, v, u, TestNbr, NumberLength);        // u = (v+b*w)
  modmult(w, u, z2);      // z2 = (w*u)
}
/* End of code adapted from Paul Zimmermann's ECM4C */

static int gcdIsOne(const limb* value)
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

static void GenerateSieve(int initial)
{
  int i;
  int j;
  int Q;
  int initModQ;
  for (i = 0; i < (10 * SIEVE_SIZE); i += SIEVE_SIZE)
  {
    (void)memcpy(&common.ecm.sieve[i], common.ecm.sieve2310, SIEVE_SIZE);
  }
#if MAX_PRIME_SIEVE == 11
  j = 5;
  Q = 13; /* Point to prime 13 */
#else
  j = 4;
  Q = 11; /* Point to prime 11 */
#endif
  do
  {
    if (initial > (Q * Q))
    {
      initModQ = initial % Q;
      if (initModQ & 1)
      {    // initModQ is odd
        i = (Q - initModQ) >> 1;
      }
      else if (initModQ == 0)
      {
        i = 0;
      }
      else
      {    // initModQ is even
        i = Q - (initModQ >> 1);
      }
      for (; i < (10 * SIEVE_SIZE); i += Q)
      {
        common.ecm.sieve[i] = 1; /* Composite */
      }
    }
    else
    {
      i = (Q * Q) - initial;
      if (i < (20 * SIEVE_SIZE))
      {
        for (i = i / 2; i < (10 * SIEVE_SIZE); i += Q)
        {
          common.ecm.sieve[i] = 1; /* Composite */
        }
      }
      else
      {
        break;
      }
    }
    Q = SmallPrime[++j];
#if MAX_PRIME_SIEVE == 11
  } while (Q < 5000);
#else
} while (Q < (10 * SIEVE_SIZE));
#endif
}

/*******************************/
/* First step of ECM algorithm */
/*******************************/
enum eEcmResult ecmStep1(void)
{
  int I;
  int P;
  int i;
  int u;
  (void)memcpy(common.ecm.Xaux, common.ecm.X, NumberSizeBytes);
  (void)memcpy(common.ecm.Zaux, common.ecm.Z, NumberSizeBytes);
  (void)memcpy(common.ecm.GcdAccumulated, MontgomeryMultR1, (NumberLength + 1) * sizeof(limb));
  for (int pass = 0; pass < 2; pass++)
  {
    /* For powers of 2 */
#ifdef __EMSCRIPTEN__
    indexPrimes = 0;
#endif
    StepECM = 1;
    for (I = 1; I <= boundStep1; I <<= 1)
    {
      duplicate(common.ecm.X, common.ecm.Z, common.ecm.X, common.ecm.Z);
    }
    for (I = 3; I <= boundStep1; I *= 3)
    {
      duplicate(common.ecm.W1, common.ecm.W2, common.ecm.X, common.ecm.Z);
      add3(common.ecm.X, common.ecm.Z, common.ecm.X, common.ecm.Z, common.ecm.W1, common.ecm.W2, common.ecm.X, common.ecm.Z);
    }

    if (pass == 0)
    {
      modmult(common.ecm.GcdAccumulated, common.ecm.Z, common.ecm.Aux1);
      (void)memcpy(common.ecm.GcdAccumulated, common.ecm.Aux1, NumberSizeBytes);
    }
    else
    {
      if (gcdIsOne(common.ecm.Z) > 1)
      {
        return FACTOR_FOUND;
      }
    }

    /* for powers of odd primes */

    indexM = 1;
    do
    {
#ifdef __EMSCRIPTEN__
      indexPrimes++;
#endif
      P = SmallPrime[indexM];
      for (long long IP = P; IP <= boundStep1; IP *= P)
      {
        prac(P, common.ecm.X, common.ecm.Z, common.ecm.W1, common.ecm.W2, common.ecm.W3, common.ecm.W4);
      }
      indexM++;
      if (pass == 0)
      {
        modmult(common.ecm.GcdAccumulated, common.ecm.Z, common.ecm.Aux1);
        (void)memcpy(common.ecm.GcdAccumulated, common.ecm.Aux1, NumberSizeBytes);
      }
      else
      {
        if (gcdIsOne(common.ecm.Z) > 1)
        {
          return FACTOR_FOUND;
        }
      }
    } while (SmallPrime[indexM - 1] <= sqrtBoundStep1);
    P += 2;

    /* Initialize sieve2310[n]: 1 if gcd(P+2n,2310) > 1, 0 otherwise */
    u = P;
    for (i = 0; i < SIEVE_SIZE; i++)
    {
      common.ecm.sieve2310[i] =
        ((((u % 3) == 0)
          || ((u % 5) == 0)
          || ((u % 7) == 0)
#if MAX_PRIME_SIEVE == 11
          || ((u % 11) == 0)
#endif
          ) ? (unsigned char)1 : (unsigned char)0);
      u += 2;
    }
    do
    {
      /* Generate sieve */
      GenerateSieve(P);

      /* Walk through sieve */

      for (i = 0; i < (10 * SIEVE_SIZE); i++)
      {
        if (common.ecm.sieve[i] != 0)
        {
          continue; /* Do not process composites */
        }
        if ((P + (2 * i)) > boundStep1)
        {
          break;
        }
#ifdef __EMSCRIPTEN__
        indexPrimes++;
#endif
        prac(P + (2 * i), common.ecm.X, common.ecm.Z, common.ecm.W1,
          common.ecm.W2, common.ecm.W3, common.ecm.W4);
        if (pass == 0)
        {
          modmult(common.ecm.GcdAccumulated, common.ecm.Z, common.ecm.Aux1);
          (void)memcpy(common.ecm.GcdAccumulated, common.ecm.Aux1, NumberSizeBytes);
        }
        else
        {
          if (gcdIsOne(common.ecm.Z) > 1)
          {
            return FACTOR_FOUND;
          }
        }
      }
      P += 20 * SIEVE_SIZE;
    } while (P < boundStep1);
    if (pass == 0)
    {
      if (BigNbrIsZero(common.ecm.GcdAccumulated))
      { // If GcdAccumulated is
        (void)memcpy(common.ecm.X, common.ecm.Xaux, NumberSizeBytes);
        (void)memcpy(common.ecm.Z, common.ecm.Zaux, NumberSizeBytes);
        continue; // multiple of TestNbr, continue.
      }
      if (gcdIsOne(common.ecm.GcdAccumulated) > 1)
      {
        return FACTOR_FOUND;
      }
      break;
    }
  } /* end for Pass */
  return FACTOR_NOT_FOUND;
}

/******************************************************/
/* Second step (using improved standard continuation) */
/******************************************************/
enum eEcmResult ecmStep2(void)
{
  int I;
  int j;
  StepECM = 2;
  j = 0;
  for (int u = 1; u < SIEVE_SIZE; u += 2)
  {
    if (((u % 3) == 0) || ((u % 5) == 0) || ((u % 7) == 0)
#if MAX_PRIME_SIEVE == 11
      || ((u % 11) == 0)
#endif
      )
    {
      common.ecm.sieve2310[u / 2] = (unsigned char)1;
    }
    else
    {
      common.ecm.sieve2310[common.ecm.sieveidx[j++] = u / 2] = (unsigned char)0;
    }
  }
  (void)memcpy(&common.ecm.sieve2310[HALF_SIEVE_SIZE], &common.ecm.sieve2310[0], HALF_SIEVE_SIZE);
  (void)memcpy(common.ecm.Xaux, common.ecm.X, NumberSizeBytes);  // (X:Z) -> Q (output
  (void)memcpy(common.ecm.Zaux, common.ecm.Z, NumberSizeBytes);  //         from step 1)
  for (int pass = 0; pass < 2; pass++)
  {
    int Qaux;
    int J;
    (void)memcpy(common.ecm.GcdAccumulated, MontgomeryMultR1, NumberSizeBytes);
    (void)memcpy(common.ecm.UX, common.ecm.X, NumberSizeBytes);
    (void)memcpy(common.ecm.UZ, common.ecm.Z, NumberSizeBytes);  // (UX:UZ) -> Q 
    ModInvBigNbr(common.ecm.Z, common.ecm.Aux1, TestNbr, NumberLength);
    modmult(common.ecm.Aux1, common.ecm.X, common.ecm.root[0]); // root[0] <- X/Z (Q)
    J = 0;
    AddBigNbrModN(common.ecm.X, common.ecm.Z, common.ecm.Aux1, TestNbr, NumberLength);
    modmult(common.ecm.Aux1, common.ecm.Aux1, common.ecm.W1);
    SubtBigNbrModN(common.ecm.X, common.ecm.Z, common.ecm.Aux1, TestNbr, NumberLength);
    modmult(common.ecm.Aux1, common.ecm.Aux1, common.ecm.W2);
    modmult(common.ecm.W1, common.ecm.W2, common.ecm.TX);
    SubtBigNbrModN(common.ecm.W1, common.ecm.W2, common.ecm.Aux1, TestNbr, NumberLength);
    modmult(common.ecm.Aux1, common.ecm.AA, common.ecm.Aux2);
    AddBigNbrModN(common.ecm.Aux2, common.ecm.W2, common.ecm.Aux3, TestNbr, NumberLength);
    modmult(common.ecm.Aux1, common.ecm.Aux3, common.ecm.TZ); // (TX:TZ) -> 2Q
    SubtBigNbrModN(common.ecm.X, common.ecm.Z, common.ecm.Aux1, TestNbr, NumberLength);
    AddBigNbrModN(common.ecm.TX, common.ecm.TZ, common.ecm.Aux2, TestNbr, NumberLength);
    modmult(common.ecm.Aux1, common.ecm.Aux2, common.ecm.W1);
    AddBigNbrModN(common.ecm.X, common.ecm.Z, common.ecm.Aux1, TestNbr, NumberLength);
    SubtBigNbrModN(common.ecm.TX, common.ecm.TZ, common.ecm.Aux2, TestNbr, NumberLength);
    modmult(common.ecm.Aux1, common.ecm.Aux2, common.ecm.W2);
    AddBigNbrModN(common.ecm.W1, common.ecm.W2, common.ecm.Aux1, TestNbr, NumberLength);
    modmult(common.ecm.Aux1, common.ecm.Aux1, common.ecm.Aux2);
    modmult(common.ecm.Aux2, common.ecm.UZ, common.ecm.X);
    SubtBigNbrModN(common.ecm.W1, common.ecm.W2, common.ecm.Aux1, TestNbr, NumberLength);
    modmult(common.ecm.Aux1, common.ecm.Aux1, common.ecm.Aux2);
    modmult(common.ecm.Aux2, common.ecm.UX, common.ecm.Z); // (X:Z) -> 3Q
    for (I = 5; I < SIEVE_SIZE; I += 2)
    {
      (void)memcpy(common.ecm.WX, common.ecm.X, NumberSizeBytes);
      (void)memcpy(common.ecm.WZ, common.ecm.Z, NumberSizeBytes);
      SubtBigNbrModN(common.ecm.X, common.ecm.Z, common.ecm.Aux1, TestNbr, NumberLength);
      AddBigNbrModN(common.ecm.TX, common.ecm.TZ, common.ecm.Aux2, TestNbr, NumberLength);
      modmult(common.ecm.Aux1, common.ecm.Aux2, common.ecm.W1);
      AddBigNbrModN(common.ecm.X, common.ecm.Z, common.ecm.Aux1, TestNbr, NumberLength);
      SubtBigNbrModN(common.ecm.TX, common.ecm.TZ, common.ecm.Aux2, TestNbr, NumberLength);
      modmult(common.ecm.Aux1, common.ecm.Aux2, common.ecm.W2);
      AddBigNbrModN(common.ecm.W1, common.ecm.W2, common.ecm.Aux1, TestNbr, NumberLength);
      modmult(common.ecm.Aux1, common.ecm.Aux1, common.ecm.Aux2);
      modmult(common.ecm.Aux2, common.ecm.UZ, common.ecm.X);
      SubtBigNbrModN(common.ecm.W1, common.ecm.W2, common.ecm.Aux1, TestNbr, NumberLength);
      modmult(common.ecm.Aux1, common.ecm.Aux1, common.ecm.Aux2);
      modmult(common.ecm.Aux2, common.ecm.UX, common.ecm.Z); // (X:Z) -> 5Q, 7Q, ...
      if (pass == 0)
      {
        modmult(common.ecm.GcdAccumulated, common.ecm.Aux1, common.ecm.Aux2);
        (void)memcpy(common.ecm.GcdAccumulated, common.ecm.Aux2, NumberSizeBytes);
      }
      else
      {
        if (gcdIsOne(common.ecm.Aux1) > 1)
        {
          return FACTOR_FOUND;
        }
      }
      if (I == HALF_SIEVE_SIZE)
      {
        (void)memcpy(common.ecm.DX, common.ecm.X, NumberSizeBytes);
        (void)memcpy(common.ecm.DZ, common.ecm.Z, NumberSizeBytes);  // (DX:DZ) -> HALF_SIEVE_SIZE*Q
      }
      if (((I % 3) != 0) && ((I % 5) != 0) && ((I % 7) != 0)
#if MAX_PRIME_SIEVE == 11
        && ((I % 11) != 0)
#endif
        )
      {
        J++;
        ModInvBigNbr(common.ecm.Z, common.ecm.Aux1, TestNbr, NumberLength);
        modmult(common.ecm.Aux1, common.ecm.X, common.ecm.root[J]); // root[J] <- X/Z
      }
      (void)memcpy(common.ecm.UX, common.ecm.WX, NumberSizeBytes);  // (UX:UZ) <-
      (void)memcpy(common.ecm.UZ, common.ecm.WZ, NumberSizeBytes);  // Previous (X:Z)
    } /* end for I */
    AddBigNbrModN(common.ecm.DX, common.ecm.DZ, common.ecm.Aux1, TestNbr, NumberLength);
    modmult(common.ecm.Aux1, common.ecm.Aux1, common.ecm.W1);
    SubtBigNbrModN(common.ecm.DX, common.ecm.DZ, common.ecm.Aux1, TestNbr, NumberLength);
    modmult(common.ecm.Aux1, common.ecm.Aux1, common.ecm.W2);
    modmult(common.ecm.W1, common.ecm.W2, common.ecm.X);
    SubtBigNbrModN(common.ecm.W1, common.ecm.W2, common.ecm.Aux1, TestNbr, NumberLength);
    modmult(common.ecm.Aux1, common.ecm.AA, common.ecm.Aux2);
    AddBigNbrModN(common.ecm.Aux2, common.ecm.W2, common.ecm.Aux3, TestNbr, NumberLength);
    modmult(common.ecm.Aux1, common.ecm.Aux3, common.ecm.Z);
    (void)memcpy(common.ecm.UX, common.ecm.X, NumberSizeBytes);
    (void)memcpy(common.ecm.UZ, common.ecm.Z, NumberSizeBytes);    // (UX:UZ) -> SIEVE_SIZE*Q
    AddBigNbrModN(common.ecm.X, common.ecm.Z, common.ecm.Aux1, TestNbr, NumberLength);
    modmult(common.ecm.Aux1, common.ecm.Aux1, common.ecm.W1);
    SubtBigNbrModN(common.ecm.X, common.ecm.Z, common.ecm.Aux1, TestNbr, NumberLength);
    modmult(common.ecm.Aux1, common.ecm.Aux1, common.ecm.W2);
    modmult(common.ecm.W1, common.ecm.W2, common.ecm.TX);
    SubtBigNbrModN(common.ecm.W1, common.ecm.W2, common.ecm.Aux1, TestNbr, NumberLength);
    modmult(common.ecm.Aux1, common.ecm.AA, common.ecm.Aux2);
    AddBigNbrModN(common.ecm.Aux2, common.ecm.W2, common.ecm.Aux3, TestNbr, NumberLength);
    modmult(common.ecm.Aux1, common.ecm.Aux3, common.ecm.TZ); // (TX:TZ) -> 2*SIEVE_SIZE*Q
    SubtBigNbrModN(common.ecm.X, common.ecm.Z, common.ecm.Aux1, TestNbr, NumberLength);
    AddBigNbrModN(common.ecm.TX, common.ecm.TZ, common.ecm.Aux2, TestNbr, NumberLength);
    modmult(common.ecm.Aux1, common.ecm.Aux2, common.ecm.W1);
    AddBigNbrModN(common.ecm.X, common.ecm.Z, common.ecm.Aux1, TestNbr, NumberLength);
    SubtBigNbrModN(common.ecm.TX, common.ecm.TZ, common.ecm.Aux2, TestNbr, NumberLength);
    modmult(common.ecm.Aux1, common.ecm.Aux2, common.ecm.W2);
    AddBigNbrModN(common.ecm.W1, common.ecm.W2, common.ecm.Aux1, TestNbr, NumberLength);
    modmult(common.ecm.Aux1, common.ecm.Aux1, common.ecm.Aux2);
    modmult(common.ecm.Aux2, common.ecm.UZ, common.ecm.X);
    SubtBigNbrModN(common.ecm.W1, common.ecm.W2, common.ecm.Aux1, TestNbr, NumberLength);
    modmult(common.ecm.Aux1, common.ecm.Aux1, common.ecm.Aux2);
    modmult(common.ecm.Aux2, common.ecm.UX, common.ecm.Z); // (X:Z) -> 3*SIEVE_SIZE*Q
    Qaux = (int)(boundStep1 / (2 * SIEVE_SIZE));
    maxIndexM = (int)(boundStep2 / (2 * SIEVE_SIZE));
    for (indexM = 0; indexM <= maxIndexM; indexM++)
    {
      if (indexM >= Qaux)
      { // If inside step 2 range... 
        if (indexM == 0)
        {
          ModInvBigNbr(common.ecm.UZ, common.ecm.Aux3, TestNbr, NumberLength);
          modmult(common.ecm.UX, common.ecm.Aux3, common.ecm.Aux1); // Aux1 <- X/Z (SIEVE_SIZE*Q)
        }
        else
        {
          ModInvBigNbr(common.ecm.Z, common.ecm.Aux3, TestNbr, NumberLength);
          modmult(common.ecm.X, common.ecm.Aux3, common.ecm.Aux1); // Aux1 <- X/Z (3,5,*
        }                         //         SIEVE_SIZE*Q)

          /* Generate sieve */
        if (((indexM % 10) == 0) || (indexM == Qaux))
        {
          GenerateSieve(indexM / 10 * (20 * SIEVE_SIZE) + 1);
        }
        /* Walk through sieve */
        J = HALF_SIEVE_SIZE + (indexM % 10) * SIEVE_SIZE;
        for (int i = 0; i < GROUP_SIZE; i++)
        {
          j = common.ecm.sieveidx[i]; // 0 < J < HALF_SIEVE_SIZE
          if ((common.ecm.sieve[J + j] != 0) && (common.ecm.sieve[J - 1 - j] != 0))
          {
            continue; // Do not process if both are composite numbers.
          }
          SubtBigNbrModN(common.ecm.Aux1, common.ecm.root[i], common.ecm.M, TestNbr, NumberLength);
          modmult(common.ecm.GcdAccumulated, common.ecm.M, common.ecm.Aux2);
          (void)memcpy(common.ecm.GcdAccumulated, common.ecm.Aux2, NumberSizeBytes);
        }
        if (pass != 0)
        {
          if (BigNbrIsZero(common.ecm.GcdAccumulated))
          {           // This curve cannot factor the number.
            break;
          }
          if (gcdIsOne(common.ecm.GcdAccumulated) > 1)
          {
            return FACTOR_FOUND;
          }
        }
      }   // End for.
      if (indexM != 0)
      { // Update (X:Z)
        (void)memcpy(common.ecm.WX, common.ecm.X, NumberSizeBytes);
        (void)memcpy(common.ecm.WZ, common.ecm.Z, NumberSizeBytes);
        SubtBigNbrModN(common.ecm.X, common.ecm.Z, common.ecm.Aux1, TestNbr, NumberLength);
        AddBigNbrModN(common.ecm.TX, common.ecm.TZ, common.ecm.Aux2, TestNbr, NumberLength);
        modmult(common.ecm.Aux1, common.ecm.Aux2, common.ecm.W1);
        AddBigNbrModN(common.ecm.X, common.ecm.Z, common.ecm.Aux1, TestNbr, NumberLength);
        SubtBigNbrModN(common.ecm.TX, common.ecm.TZ, common.ecm.Aux2, TestNbr, NumberLength);
        modmult(common.ecm.Aux1, common.ecm.Aux2, common.ecm.W2);
        AddBigNbrModN(common.ecm.W1, common.ecm.W2, common.ecm.Aux1, TestNbr, NumberLength);
        modmult(common.ecm.Aux1, common.ecm.Aux1, common.ecm.Aux2);
        modmult(common.ecm.Aux2, common.ecm.UZ, common.ecm.X);
        SubtBigNbrModN(common.ecm.W1, common.ecm.W2, common.ecm.Aux1, TestNbr, NumberLength);
        modmult(common.ecm.Aux1, common.ecm.Aux1, common.ecm.Aux2);
        modmult(common.ecm.Aux2, common.ecm.UX, common.ecm.Z);
        (void)memcpy(common.ecm.UX, common.ecm.WX, NumberSizeBytes);
        (void)memcpy(common.ecm.UZ, common.ecm.WZ, NumberSizeBytes);
      }
    } // end for Q
    if (pass == 0)
    {
      int rc;
      if (BigNbrIsZero(common.ecm.GcdAccumulated))
      { // If GcdAccumulated is zero
        (void)memcpy(common.ecm.X, common.ecm.Xaux, NumberSizeBytes);
        (void)memcpy(common.ecm.Z, common.ecm.Zaux, NumberSizeBytes);
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
      // GD <- GCD(GcdAccumulated, TestNbr)
      if (memcmp(common.ecm.GD, TestNbr, NumberSizeBytes))
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
  SmallPrime[0] = 2;
  indexM = 1;
  for (indexM = 1; indexM < sizeof(SmallPrime) / sizeof(SmallPrime[0]); indexM++)
  {     // Loop that fills the SmallPrime array.
    int divisor;
    SmallPrime[indexM] = potentialPrime; /* Store prime */
    do
    {
      potentialPrime += 2;
      for (divisor = 3; divisor * divisor <= potentialPrime; divisor += 2)
      { /* Check if P is prime */
        if (potentialPrime % divisor == 0)
        {
          break;  /* Composite */
        }
      }
    } while (divisor * divisor <= potentialPrime);
  }
}

enum eEcmResult ecmCurve(int *pEC, int *pNextEC)
{
  int EC = *pEC;
  int NextEC = *pNextEC;
  enum eEcmResult result;
  NumberSizeBytes = NumberLength * sizeof(limb);
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
        common.ecm.GD[0].x = 1;   // Set GD to 1.
        (void)memset(&common.ecm.GD[1], 0, (NumberLength - 1) * sizeof(limb));
        *pEC = EC;
        *pNextEC = NextEC;
        return FACTOR_FOUND;
      }
    }
    else
    {
      EC++;
#ifdef __EMSCRIPTEN__
      text[0] = '7';
      ptrText = &text[1];
      int2dec(&ptrText, EC);
      *ptrText = 0;                 // Add string terminator.
      databack(text);
#endif
      nbrDigits = NumberLength * 9;          // Get number of digits.
      if ((NextEC == 0) && (nbrDigits <= 110)) // Force switch to SIQS and number not too large.
      {
        EC += TYP_SIQS;
        *pEC = EC;
        *pNextEC = NextEC;
        return CHANGE_TO_SIQS;
      }
      if ((nbrDigits > 30) && (nbrDigits <= 90))  // If between 30 and 90 digits...         
      {                             // Switch to SIQS.
        int limit = limits[(nbrDigits - 31) / 5];
        if (EC % 50000000 >= limit)
        {                           // Switch to SIQS.
          EC += TYP_SIQS;
          *pEC = EC;
          *pNextEC = NextEC;
          return CHANGE_TO_SIQS;
        }
      }
    }
    boundStep1 = 2000;
    boundStep2 = 200000;
    sqrtBoundStep1 = 45;
#ifdef __EMSCRIPTEN__
    nbrPrimes = 303; /* Number of primes less than 2000 */
#endif
    if (EC > 25)
    {
      if (EC < 326)
      {
        boundStep1 = 50000;
        boundStep2 = 5000000;
        sqrtBoundStep1 = 224;
#ifdef __EMSCRIPTEN__
        nbrPrimes = 5133; /* Number of primes less than 50000 */
#endif
      }
      else
      {
        if (EC < 2000)
        {
          boundStep1 = 1000000;
          boundStep2 = 100000000;
          sqrtBoundStep1 = 1001;
#ifdef __EMSCRIPTEN__
          nbrPrimes = 78498; /* Number of primes less than 1000000 */
#endif
        }
        else
        {
          boundStep1 = 11000000;
          boundStep2 = 1100000000;
          sqrtBoundStep1 = 3316;
#ifdef __EMSCRIPTEN__
          nbrPrimes = 726517; /* Number of primes less than 11000000 */
#endif
        }
      }
    }
#ifdef __EMSCRIPTEN__
    ptrText = ptrLowerText;  // Point after number that is being factored.
    (void)strcpy(ptrText, lang ? "<p>Curva " : "<p>Curve ");
    ptrText += strlen(ptrText);
    int2dec(&ptrText, EC);   // Show curve number.
    (void)strcpy(ptrText, lang ? " usando límites B1=" : " using bounds B1=");
    ptrText += strlen(ptrText);
    int2dec(&ptrText, boundStep1);   // Show first bound.
    (void)strcpy(ptrText, lang ? " y B2=" : " and B2=");
    ptrText += strlen(ptrText);
    int2dec(&ptrText, boundStep2);   // Show second bound.
    (void)strcpy(ptrText, "</p>");
    databack(lowerText);
#endif

    //  Compute A0 <- 2 * (EC+1)*modinv(3 * (EC+1) ^ 2 - 1, N) mod N
                                               // Aux2 <- 1 in Montgomery notation.
    (void)memcpy(common.ecm.Aux2, MontgomeryMultR1, NumberSizeBytes);
    modmultInt(common.ecm.Aux2, EC + 1, common.ecm.Aux2);            // Aux2 <- EC + 1.
    modmultInt(common.ecm.Aux2, 2, common.ecm.Aux1);                 // Aux1 <- 2*(EC+1)
    modmultInt(common.ecm.Aux2, EC + 1, common.ecm.Aux3);            // Aux3 <- (EC + 1)^2
    modmultInt(common.ecm.Aux3, 3, common.ecm.Aux3);                 // Aux3 <- 3*(EC + 1)^2
                                               // Aux2 <- 3*(EC + 1)^2 - 1 
    SubtBigNbrModN(common.ecm.Aux3, MontgomeryMultR1, common.ecm.Aux2, TestNbr, NumberLength);
    ModInvBigNbr(common.ecm.Aux2, common.ecm.Aux2, TestNbr, NumberLength);
    modmult(common.ecm.Aux1, common.ecm.Aux2, common.ecm.A0);                   // A0 <- 2*(EC+1)/(3*(EC+1)^2 - 1)

    //  if A0*(A0 ^ 2 - 1)*(9 * A0 ^ 2 - 1) mod N=0 then select another curve.
    modmult(common.ecm.A0, common.ecm.A0, common.ecm.A02);          // A02 <- A0^2
    modmult(common.ecm.A02, common.ecm.A0, common.ecm.A03);         // A03 <- A0^3
    SubtBigNbrModN(common.ecm.A03, common.ecm.A0, common.ecm.Aux1, TestNbr, NumberLength);  // Aux1 <- A0^3 - A0
    modmultInt(common.ecm.A02, 9, common.ecm.Aux2);      // Aux2 <- 9*A0^2
    SubtBigNbrModN(common.ecm.Aux2, MontgomeryMultR1, common.ecm.Aux2, TestNbr, NumberLength); // Aux2 <- 9*A0^2-1
    modmult(common.ecm.Aux1, common.ecm.Aux2, common.ecm.Aux3);
  } while (BigNbrIsZero(common.ecm.Aux3));
  //   Z <- 4 * A0 mod N
  modmultInt(common.ecm.A0, 4, common.ecm.Z);
  //   A = (-3 * A0 ^ 4 - 6 * A0 ^ 2 + 1)*modinv(4 * A0 ^ 3, N) mod N
  modmultInt(common.ecm.A02, 6, common.ecm.Aux1);      // Aux1 <- 6*A0^2
  SubtBigNbrModN(MontgomeryMultR1, common.ecm.Aux1, common.ecm.Aux1, TestNbr, NumberLength);
  modmult(common.ecm.A02, common.ecm.A02, common.ecm.Aux2);       // Aux2 <- A0^4
  modmultInt(common.ecm.Aux2, 3, common.ecm.Aux2);     // Aux2 <- 3*A0^4
  SubtBigNbrModN(common.ecm.Aux1, common.ecm.Aux2, common.ecm.Aux1, TestNbr, NumberLength);
  modmultInt(common.ecm.A03, 4, common.ecm.Aux2);      // Aux2 <- 4*A0^3
  ModInvBigNbr(common.ecm.Aux2, common.ecm.Aux3, TestNbr, NumberLength);
  modmult(common.ecm.Aux1, common.ecm.Aux3, common.ecm.A0);
  //   AA <- (A + 2)*modinv(4, N) mod N
  modmultInt(MontgomeryMultR1, 2, common.ecm.Aux2);  // Aux2 <- 2
  AddBigNbrModN(common.ecm.A0, common.ecm.Aux2, common.ecm.Aux1, TestNbr, NumberLength); // Aux1 <- A0+2
  modmultInt(MontgomeryMultR1, 4, common.ecm.Aux2);  // Aux2 <- 4
  ModInvBigNbr(common.ecm.Aux2, common.ecm.Aux2, TestNbr, NumberLength);
  modmult(common.ecm.Aux1, common.ecm.Aux2, common.ecm.AA);
  //   X <- (3 * A0 ^ 2 + 1) mod N
  modmultInt(common.ecm.A02, 3, common.ecm.Aux1);    // Aux1 <- 3*A0^2
  AddBigNbrModN(common.ecm.Aux1, MontgomeryMultR1, common.ecm.X, TestNbr, NumberLength);
  result = ecmStep1();
  if (result == FACTOR_FOUND)
  {
    *pEC = EC;
    *pNextEC = NextEC;
    return result;
  }
  result = ecmStep2();
  if (result == FACTOR_FOUND)
  {
    *pEC = EC;
    *pNextEC = NextEC;
    return result;
  }
  *pEC = EC;
  *pNextEC = NextEC;
  return FACTOR_NOT_FOUND;
}


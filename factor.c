/*
This file is part of Alpertron Calculators.

Copyright 2015 Dario Alejandro Alpern

Alpertron Calculators is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Alpertron Calculators is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Alpertron Calculators.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "bignbr.h"
#include "expression.h"
#include "factor.h"
#include "commonstruc.h"
#include "skiptest.h"

int yieldFreq;
static int oldNbrFactors;
#ifdef __EMSCRIPTEN__
char upperText[MAX_LEN*16];
char lowerText[MAX_LEN*16];
char *ptrLowerText;
extern mmCback modmultCallback;
extern int64_t lModularMult;
extern char *ptrInputText;
#endif

#define TYP_AURIF    100000000
#define TYP_TABLE    150000000
#define TYP_SIQS     200000000
#define TYP_LEHMAN   250000000
#define TYP_RABIN    300000000
#define TYP_DIVISION 350000000
#define TYP_EC       400000000

union uCommon common;
int64_t primeModMult;
int StepECM;
int skipPrimality;
int nbrECM;
int nbrPrimalityTests;
int nbrSIQS;
int timeECM;
int timePrimalityTests;
int timeSIQS;
static int nbrPrimes, indexPrimes, DegreeAurif, NextEC;
static BigInteger power, prime;
int *factorArr[FACTOR_ARRSIZE];
static int indexM, maxIndexM;
static int foundByLehman, performLehman;
static int EC;
static int SmallPrime[670]; /* Primes < 5000 */
static void add3(limb *x3, limb *z3, limb *x2, limb *z2, limb *x1, limb *z1, limb *x, limb *z);
static void duplicate(limb *x2, limb *z2, limb *x1, limb *z1);
static BigInteger Temp1, Temp2, Temp3, Temp4;
BigInteger factorValue, tofactor;
char prettyprint, cunningham, hexadecimal;
long long Gamma[386];
long long Delta[386];
long long AurifQ[386];
char tofactorDec[MAX_LEN*12];
extern int q[MAX_LEN];
int nbrToFactor[MAX_LEN];
struct sFactors astFactorsMod[5000];
int factorsMod[20000];
static void insertBigFactor(struct sFactors *pstFactors, BigInteger *divisor,
  int type);
static char *findChar(char *str, char c);
enum eEcmResult
{
  FACTOR_NOT_FOUND = 0,
  FACTOR_FOUND,
  CHANGE_TO_SIQS,
};

/* ECM limits for 30, 35, ..., 95 digits */
static int limits[] = { 10, 10, 10, 10, 10, 15, 22, 26, 35, 50, 100, 150, 250 };
/******************************************************/
/* Start of code adapted from Paul Zimmermann's ECM4C */
/******************************************************/
#define ADD 6  /* number of multiplications in an addition */
#define DUP 5  /* number of multiplications in a duplicate */

static void GetYieldFrequency(void)
{
  yieldFreq = 1000000 / (NumberLength * NumberLength) + 1;
  if (yieldFreq > 100000)
  {
    yieldFreq = yieldFreq / 100000 * 100000;
  }
  else if (yieldFreq > 10000)
  {
    yieldFreq = yieldFreq / 10000 * 10000;
  }
  else if (yieldFreq > 1000)
  {
    yieldFreq = yieldFreq / 1000 * 1000;
  }
  else if (yieldFreq > 100)
  {
    yieldFreq = yieldFreq / 100 * 100;
  }
}
                          /* returns the number of modular multiplications */
static int lucas_cost(int n, double v)
{
  int c, d, e, r;

  d = n;
  r = (int)((double)d / v + 0.5);
  if (r >= n)
  {
    return (ADD * n);
  }
  d = n - r;
  e = 2 * r - n;
  c = DUP + ADD; /* initial duplicate and final addition */
  while (d != e)
  {
    if (d < e)
    {
      r = d;
      d = e;
      e = r;
    }
    if (4 * d <= 5 * e && ((d + e) % 3) == 0)
    { /* condition 1 */
      r = (2 * d - e) / 3;
      e = (2 * e - d) / 3;
      d = r;
      c += 3 * ADD; /* 3 additions */
    }
    else if (4 * d <= 5 * e && (d - e) % 6 == 0)
    { /* condition 2 */
      d = (d - e) / 2;
      c += ADD + DUP; /* one addition, one duplicate */
    }
    else if (d <= (4 * e))
    { /* condition 3 */
      d -= e;
      c += ADD; /* one addition */
    }
    else if ((d + e) % 2 == 0)
    { /* condition 4 */
      d = (d - e) / 2;
      c += ADD + DUP; /* one addition, one duplicate */
    }
    else if (d % 2 == 0)
    { /* condition 5 */
      d /= 2;
      c += ADD + DUP; /* one addition, one duplicate */
    }
    else if (d % 3 == 0)
    { /* condition 6 */
      d = d / 3 - e;
      c += 3 * ADD + DUP; /* three additions, one duplicate */
    }
    else if ((d + e) % 3 == 0)
    { /* condition 7 */
      d = (d - 2 * e) / 3;
      c += 3 * ADD + DUP; /* three additions, one duplicate */
    }
    else if ((d - e) % 3 == 0)
    { /* condition 8 */
      d = (d - e) / 3;
      c += 3 * ADD + DUP; /* three additions, one duplicate */
    }
    else if (e % 2 == 0)
    { /* condition 9 */
      e /= 2;
      c += ADD + DUP; /* one addition, one duplicate */
    }
  }
  return c;
}

/* computes nP from P=(x:z) and puts the result in (x:z). Assumes n>2. */
void prac(int n, limb *x, limb *z, limb *xT, limb *zT, limb *xT2, limb *zT2)
{
  int d, e, r, i;
  limb *t;
  limb *xA = x, *zA = z;
  limb *xB = common.ecm.Aux1, *zB = common.ecm.Aux2;
  limb *xC = common.ecm.Aux3, *zC = common.ecm.Aux4;
  double v[] =
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
  r = (int)((double)d / v[i] + 0.5);
  /* first iteration always begins by Condition 3, then a swap */
  d = n - r;
  e = 2 * r - n;
  memcpy(xB, xA, NumberLength * sizeof(limb));   // B <- A
  memcpy(zB, zA, NumberLength * sizeof(limb));
  memcpy(xC, xA, NumberLength * sizeof(limb));   // C <- A
  memcpy(zC, zA, NumberLength * sizeof(limb));
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
    if (4 * d <= 5 * e && ((d + e) % 3) == 0)
    { /* condition 1 */
      r = (2 * d - e) / 3;
      e = (2 * e - d) / 3;
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
    else if (4 * d <= 5 * e && (d - e) % 6 == 0)
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
    else if ((d + e) % 2 == 0)
    { /* condition 4 */
      d = (d - e) / 2;
      add3(xB, zB, xB, zB, xA, zA, xC, zC); /* B = f(B,A,C) */
      duplicate(xA, zA, xA, zA); /* A = 2*A */
    }
    else if (d % 2 == 0)
    { /* condition 5 */
      d /= 2;
      add3(xC, zC, xC, zC, xA, zA, xB, zB); /* C = f(C,A,B) */
      duplicate(xA, zA, xA, zA); /* A = 2*A */
    }
    else if (d % 3 == 0)
    { /* condition 6 */
      d = d / 3 - e;
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
    else if ((d + e) % 3 == 0)
    { /* condition 7 */
      d = (d - 2 * e) / 3;
      add3(xT, zT, xA, zA, xB, zB, xC, zC); /* T1 = f(A,B,C) */
      add3(xB, zB, xT, zT, xA, zA, xB, zB); /* B = f(T1,A,B) */
      duplicate(xT, zT, xA, zA);
      add3(xA, zA, xA, zA, xT, zT, xA, zA); /* A = 3*A */
    }
    else if ((d - e) % 3 == 0)
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
    else if (e % 2 == 0)
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
static void add3(limb *x3, limb *z3, limb *x2, limb *z2,
                 limb *x1, limb *z1, limb *x, limb *z)
{
  limb *t = common.ecm.fieldTX;
  limb *u = common.ecm.fieldTZ;
  limb *v = common.ecm.fieldUX;
  limb *w = common.ecm.fieldUZ;
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
  if (!memcmp(x, x3, NumberLength*sizeof(limb)))
  {
    memcpy(u, x, NumberLength * sizeof(int));
    memcpy(t, w, NumberLength * sizeof(int));
    modmult(z, t, w);
    modmult(v, u, z3);
    memcpy(x3, w, NumberLength * sizeof(int));
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
static void duplicate(limb *x2, limb *z2, limb *x1, limb *z1)
{
  limb *u = common.ecm.fieldUZ;
  limb *v = common.ecm.fieldTX;
  limb *w = common.ecm.fieldTZ;
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

static int gcdIsOne(limb *value)
{
  UncompressLimbsBigInteger(value, &Temp1);
  UncompressLimbsBigInteger(TestNbr, &Temp2);
  // Return zero if value is zero or both numbers are equal.
  if (Temp1.nbrLimbs == 1 && Temp1.limbs[0].x == 0)
  {
    return 0;
  }
  BigIntSubt(&Temp1, &Temp2, &Temp3);
  if (Temp3.nbrLimbs == 1 && Temp3.limbs[0].x == 0)
  {
    return 0;
  }
  BigIntGcd(&Temp1, &Temp2, &Temp3);
  CompressLimbsBigInteger(common.ecm.GD, &Temp3);
  if (Temp3.nbrLimbs == 1 && Temp3.limbs[0].x < 2)
  {
    return Temp3.limbs[0].x;    // GCD is less than 2.
  }
  return 2;      // GCD is greater than one.
}

static void GenerateSieve(int initial)
{
  int i, j, Q, initModQ;
  for (i = 0; i < 10*SIEVE_SIZE; i += SIEVE_SIZE)
  {
    memcpy(&common.ecm.sieve[i], common.ecm.sieve2310, SIEVE_SIZE);
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
    if (initial > Q * Q)
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
      for (; i < 10*SIEVE_SIZE; i += Q)
      {
        common.ecm.sieve[i] = 1; /* Composite */
      }
    }
    else
    {
      i = Q * Q - initial;
      if (i < 20*SIEVE_SIZE)
      {
        for (i = i / 2; i < 10*SIEVE_SIZE; i += Q)
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
  } while (Q < 10*SIEVE_SIZE);
#endif
}

static int Cos(int N)
{
  switch (N % 8)
  {
  case 0:
    return 1;
  case 4:
    return -1;
  }
  return 0;
}

static int intTotient(int N)
{
  int totient, argumentDivisor, trialDivisor;

  totient = argumentDivisor = N;
  if (argumentDivisor % 2 == 0)
  {
    totient /= 2;
    do
    {
      argumentDivisor /= 2;
    } while (argumentDivisor % 2 == 0);
  }
  if (argumentDivisor % 3 == 0)
  {
    totient = totient * 2 / 3;
    do
    {
      argumentDivisor /= 3;
    } while (argumentDivisor % 3 == 0);
  }
  trialDivisor = 5;
  while (trialDivisor * trialDivisor <= argumentDivisor)
  {
    if (trialDivisor % 3 != 0 && argumentDivisor % trialDivisor == 0)
    {
      totient = totient * (trialDivisor - 1) / trialDivisor;
      do
      {
        argumentDivisor /= trialDivisor;
      } while (argumentDivisor % trialDivisor == 0);
    }
    trialDivisor += 2;
  }
  if (argumentDivisor > 1)
  {
    totient = totient * (argumentDivisor - 1) / argumentDivisor;
  }
  return totient;
}

int Moebius(int N)
{
  int moebius, argumentDivisor, trialDivisor;

  moebius = 1;
  argumentDivisor = N;
  if (argumentDivisor % 2 == 0)
  {
    moebius = -moebius;
    argumentDivisor /= 2;
    if (argumentDivisor % 2 == 0)
    {
      return 0;
    }
  }
  if (argumentDivisor % 3 == 0)
  {
    moebius = -moebius;
    argumentDivisor /= 3;
    if (argumentDivisor % 3 == 0)
    {
      return 0;
    }
  }
  trialDivisor = 5;
  while (trialDivisor * trialDivisor <= argumentDivisor)
  {
    if (trialDivisor % 3 != 0)
    {
      while (argumentDivisor % trialDivisor == 0)
      {
        moebius = -moebius;
        argumentDivisor /= trialDivisor;
        if (argumentDivisor % trialDivisor == 0)
        {
          return 0;
        }
      }
    }
    trialDivisor += 2;
  }
  if (argumentDivisor > 1)
  {
    moebius = -moebius;
  }
  return moebius;
}

void GetAurifeuilleFactor(struct sFactors *pstFactors, int L, BigInteger *BigBase)
{
  static BigInteger x, Csal, Dsal, Nbr1;
  int k;

  BigIntPowerIntExp(BigBase, L, &x);   // x <- BigBase^L.
  intToBigInteger(&Csal, 1);
  intToBigInteger(&Dsal, 1);
  for (k = 1; k < DegreeAurif; k++)
  {
    longToBigInteger(&Nbr1, Gamma[k]);
    BigIntMultiply(&Csal, &x, &Csal);
    BigIntAdd(&Csal, &Nbr1, &Csal);      // Csal <- Csal * x + Gamma[k]
    longToBigInteger(&Nbr1, Delta[k]);
    BigIntMultiply(&Dsal, &x, &Dsal);
    BigIntAdd(&Dsal, &Nbr1, &Dsal);      // Dsal <- Dsal * x + Gamma[k]
  }
  longToBigInteger(&Nbr1, Gamma[k]);
  BigIntMultiply(&Csal, &x, &Csal);
  BigIntAdd(&Csal, &Nbr1, &Csal);        // Csal <- Csal * x + Gamma[k]
  BigIntPowerIntExp(BigBase, (L + 1) / 2, &Nbr1);   // Nbr1 <- Dsal * base^((L+1)/2)
  BigIntMultiply(&Dsal, &Nbr1, &Nbr1);
  BigIntAdd(&Csal, &Nbr1, &Dsal);
  insertBigFactor(pstFactors, &Dsal, TYP_AURIF);
  BigIntSubt(&Csal, &Nbr1, &Dsal);
  insertBigFactor(pstFactors, &Dsal, TYP_AURIF);
}

// Get Aurifeuille factors.
void InsertAurifFactors(struct sFactors *pstFactors, BigInteger *BigBase, int Expon, int Incre)
{
  int Base = BigBase->limbs[0].x;
  if (BigBase->nbrLimbs != 1 || Base >= 386)
  {
    return;    // Base is very big, so go out.
  }
  if (Expon % 2 == 0 && Incre == -1)
  {
    do
    {
      Expon /= 2;
    } while (Expon % 2 == 0);
    Incre = Base % 4 - 2;
  }
  if (Expon % Base == 0
    && Expon / Base % 2 != 0
    && ((Base % 4 != 1 && Incre == 1) || (Base % 4 == 1 && Incre == -1)))
  {
    int N1, q, L, k;
    int N = Base;
    if (N % 4 == 1)
    {
      N1 = N;
    }
    else
    {
      N1 = 2 * N;
    }
    DegreeAurif = intTotient(N1) / 2;
    for (k = 1; k <= DegreeAurif; k += 2)
    {
      AurifQ[k] = JacobiSymbol(N, k);
    }
    for (k = 2; k <= DegreeAurif; k += 2)
    {
      int t1 = k; // Calculate t2 = gcd(k, N1)
      int t2 = N1;
      while (t1 != 0)
      {
        int t3 = t2 % t1;
        t2 = t1;
        t1 = t3;
      }
      t1 = Moebius(N1 / t2) * intTotient(t2) * Cos((N - 1) * k);
      AurifQ[k] = t1;
    }
    Gamma[0] = Delta[0] = 1;
    for (k = 1; k <= DegreeAurif / 2; k++)
    {
      int j;
      Gamma[k] = Delta[k] = 0;
      for (j = 0; j < k; j++)
      {
        Gamma[k] =
          Gamma[k]
          + N * AurifQ[2 * k
          - 2 * j
          - 1] * Delta[j]
          - AurifQ[2 * k
          - 2 * j] * Gamma[j];
        Delta[k] =
          Delta[k]
          + AurifQ[2 * k
          + 1
          - 2 * j] * Gamma[j]
          - AurifQ[2 * k
          - 2 * j] * Delta[j];
      }
      Gamma[k] /= 2 * k;
      Delta[k] = (Delta[k] + Gamma[k]) / (2 * k + 1);
    }
    for (k = DegreeAurif / 2 + 1; k <= DegreeAurif; k++)
    {
      Gamma[k] = Gamma[DegreeAurif - k];
    }
    for (k = (DegreeAurif + 1) / 2; k < DegreeAurif; k++)
    {
      Delta[k] = Delta[DegreeAurif - k - 1];
    }
    q = Expon / Base;
    L = 1;
    while (L * L <= q)
    {
      if (q % L == 0)
      {
        GetAurifeuilleFactor(pstFactors,L, BigBase);
        if (q != L * L)
        {
          GetAurifeuilleFactor(pstFactors, q / L, BigBase);
        }
      }
      L += 2;
    }
  }
}

void copyString(char *textFromServer)
{
  strcpy(common.saveFactors.text, textFromServer);
}

void Cunningham(struct sFactors *pstFactors, BigInteger *BigBase, int Expon,
                int increment, BigInteger *BigOriginal)
{
#ifdef __EMSCRIPTEN__
  char url[200];
  char *ptrUrl;
#endif
  int Expon2, k;
  char *ptrFactorsAscii;
  BigInteger Nbr1, Nbr2;

  common.saveFactors.text[0] = 0;    // Indicate no new factor found in advance.
  Expon2 = Expon;
  if (cunningham && BigOriginal->nbrLimbs > 4)
  {   // Enter here on numbers of more than 40 digits if the user selected
      // get Cunningham factors from server.
#ifdef __EMSCRIPTEN__
    databack(lang ? "4<p>Obteniendo los factores primitivos conocidos del servidor Web.</p>":
                    "4<p>Requesting known primitive factors from Web server.</p>");
    // Format URL.
    ptrUrl = url;
    strcpy(ptrUrl, "factors.pl?base=");
    ptrUrl += strlen(ptrUrl);
    int2dec(&ptrUrl, BigBase->limbs[0].x);
    strcpy(ptrUrl, "&expon=");
    ptrUrl += strlen(ptrUrl);
    int2dec(&ptrUrl, Expon);
    strcpy(ptrUrl, "&type=");
    ptrUrl += strlen(ptrUrl);
    *ptrUrl++ = (increment > 0? 'p': 'm');
    *ptrUrl = 0;
    getCunn(url, common.saveFactors.text);
#endif
  }
  ptrFactorsAscii = common.saveFactors.text;
  while (*ptrFactorsAscii > ' ')
  { // Loop through factors found in server.
    int nbrDigits;
    char *ptrEndFactor = findChar(ptrFactorsAscii, '*');
    if (ptrEndFactor == NULL)
    {
      nbrDigits = (int)strlen(ptrFactorsAscii);
      ptrEndFactor = ptrFactorsAscii + nbrDigits;
    }
    else
    {
      nbrDigits = (int)(ptrEndFactor - ptrFactorsAscii);
      ptrEndFactor++;
    }
    Dec2Bin(ptrFactorsAscii, Nbr1.limbs, nbrDigits, &Nbr1.nbrLimbs);
    Nbr1.sign = SIGN_POSITIVE;
    ptrFactorsAscii = ptrEndFactor;
    insertBigFactor(pstFactors, &Nbr1, TYP_TABLE);
  }
  while (Expon2 % 2 == 0 && increment == -1)
  {
    Expon2 /= 2;
    BigIntPowerIntExp(BigBase, Expon2, &Nbr1);
    addbigint(&Nbr1, increment);
    insertBigFactor(pstFactors, &Nbr1, TYP_TABLE);
    InsertAurifFactors(pstFactors,BigBase, Expon2, 1);
  }
  k = 1;
  while (k * k <= Expon)
  {
    if (Expon % k == 0)
    {
      if (k % 2 != 0)
      { /* Only for odd exponent */
        BigIntPowerIntExp(BigBase, Expon / k, &Nbr1);
        addbigint(&Nbr1, increment);
        BigIntGcd(&Nbr1, BigOriginal, &Nbr2);   // Nbr2 <- gcd(Base^(Expon/k)+incre, original)
        insertBigFactor(pstFactors, &Nbr2, TYP_TABLE);
        CopyBigInt(&Temp1, BigOriginal);
        BigIntDivide(&Temp1, &Nbr2, &Nbr1);
        insertBigFactor(pstFactors, &Nbr1, TYP_TABLE);
        InsertAurifFactors(pstFactors, BigBase, Expon / k, increment);
      }
      if ((Expon / k) % 2 != 0)
      { /* Only for odd exponent */
        BigIntPowerIntExp(BigBase, k, &Nbr1);
        addbigint(&Nbr1, increment); 
        BigIntGcd(&Nbr1, BigOriginal, &Nbr2);   // Nbr2 <- gcd(Base^k+incre, original)
        insertBigFactor(pstFactors, &Nbr2, TYP_TABLE);
        CopyBigInt(&Temp1, BigOriginal);
        BigIntDivide(&Temp1, &Nbr2, &Nbr1);
        insertBigFactor(pstFactors, &Nbr1, TYP_TABLE);
        InsertAurifFactors(pstFactors, BigBase, k, increment);
      }
    }
    k++;
  }
}

static boolean ProcessExponent(struct sFactors *pstFactors, BigInteger *nbrToFactor, int Exponent)
{
#ifdef __EMSCRIPTEN__
  char status[200];
  char *ptrStatus;
#endif
  BigInteger NFp1, NFm1, nthRoot, rootN1, rootN, rootbak;
  BigInteger nextroot, dif;
  double log2N;
#ifdef __EMSCRIPTEN__
  int elapsedTime = (int)(tenths() - originalTenthSecond);
  if (elapsedTime / 10 != oldTimeElapsed / 10)
  {
    oldTimeElapsed = elapsedTime;
    ptrStatus = status;
    strcpy(ptrStatus, lang ? "4<p>Transcurrió " : "4<p>Time elapsed: ");
    ptrStatus += strlen(ptrStatus);
    GetDHMS(&ptrStatus, elapsedTime / 10);
    strcpy(ptrStatus, lang ? "&nbsp;&nbsp;&nbsp;Exponente potencia +/- 1: " :
      "&nbsp;&nbsp;&nbsp;Power +/- 1 exponent: ");
    ptrStatus += strlen(ptrStatus);
    int2dec(&ptrStatus, Exponent);
    *ptrStatus = 0;                       // Add string terminator.
    databack(status);
}
#endif
  CopyBigInt(&NFp1, nbrToFactor);
  addbigint(&NFp1, 1);                    // NFp1 <- NumberToFactor + 1
  CopyBigInt(&NFm1, nbrToFactor);
  addbigint(&NFm1, -1);                   // NFm1 <- NumberToFactor - 1
  log2N = logBigNbr(&NFp1) / Exponent;    // Find nth root of number to factor.
  expBigNbr(&nthRoot, log2N);
  rootbak = nthRoot;
  for (;;)
  {
    BigIntPowerIntExp(&nthRoot, Exponent - 1, &rootN1); // rootN1 <- nthRoot ^ (Exponent-1)
    BigIntMultiply(&nthRoot, &rootN1, &rootN);     // rootN <- nthRoot ^ Exponent
    BigIntSubt(&NFp1, &rootN, &dif);            // dif <- NFp1 - rootN
    if (dif.nbrLimbs == 1 && dif.limbs[0].x == 0)
    { // Perfect power
      Cunningham(pstFactors, &nthRoot, Exponent, -1, nbrToFactor);
      return TRUE;
    }
    addbigint(&dif, 1);                         // dif <- dif + 1
    BigIntDivide(&dif, &rootN1, &Temp1);        // Temp1 <- dif / rootN1
    subtractdivide(&Temp1, 0, Exponent);        // Temp1 <- Temp1 / Exponent
    BigIntAdd(&Temp1, &nthRoot, &nextroot);     // nextroot <- Temp1 + nthRoot
    addbigint(&nextroot, -1);                   // nextroot <- nextroot - 1
    BigIntSubt(&nextroot, &nthRoot, &nthRoot);  // nthRoot <- nextroot - nthRoot
    if (nthRoot.sign == SIGN_POSITIVE)
    {
      break; // Not a perfect power
    }
    CopyBigInt(&nthRoot, &nextroot);
  }
  nthRoot = rootbak;
  for (;;)
  {
    BigIntPowerIntExp(&nthRoot, Exponent - 1, &rootN1); // rootN1 <- nthRoot ^ (Exponent-1)
    BigIntMultiply(&nthRoot, &rootN1, &rootN);     // rootN <- nthRoot ^ Exponent
    BigIntSubt(&NFm1, &rootN, &dif);            // dif <- NFm1 - rootN
    if (dif.nbrLimbs == 1 && dif.limbs[0].x == 0)
    { // Perfect power
      Cunningham(pstFactors, &nthRoot, Exponent, 1, nbrToFactor);
      return TRUE;
    }
    addbigint(&dif, 1);                         // dif <- dif + 1
    BigIntDivide(&dif, &rootN1, &Temp1);        // Temp1 <- dif / rootN1
    subtractdivide(&Temp1, 0, Exponent);        // Temp1 <- Temp1 / Exponent
    BigIntAdd(&Temp1, &nthRoot, &nextroot);     // nextroot <- Temp1 + nthRoot
    addbigint(&nextroot, -1);                   // nextroot <- nextroot - 1
    BigIntSubt(&nextroot, &nthRoot, &nthRoot);  // nthRoot <- nextroot - nthRoot
    if (nthRoot.sign == SIGN_POSITIVE)
    {
      break;                                    // Not a perfect power
    }
    CopyBigInt(&nthRoot, &nextroot);
  }
  return FALSE;
}

static void PowerPM1Check(struct sFactors *pstFactors, BigInteger *nbrToFactor)
{
  char plus1 = FALSE;
  char minus1 = FALSE;
  int Exponent = 0;
  int i, j;
  int modulus;
  int mod9 = getRemainder(nbrToFactor, 9);
  int maxExpon = nbrToFactor->nbrLimbs * BITS_PER_GROUP;
  int numPrimes = 2 * maxExpon + 3;
  double logar = logBigNbr(nbrToFactor);
  // 33219 = logarithm base 2 of max number supported = 10^10000.
  // Let n = a^b +/- 1 (n = number to factor).
  // If n!=1 or n!=0 or n!=7 (mod 8), then b cannot be even.
  modulus = nbrToFactor->limbs[0].x & 7;
  if (modulus == 0 || modulus == 1 || modulus == 7)
  {       // b can be even
    memset(common.ecm.ProcessExpon, 0xFF, sizeof(common.ecm.ProcessExpon));
  }
  else
  {       // b cannot be even
    memset(common.ecm.ProcessExpon, 0xAA, sizeof(common.ecm.ProcessExpon));
  }
  memset(common.ecm.primes, 0xFF, sizeof(common.ecm.primes));
  for (i = 2; i * i < numPrimes; i++)
  {       // Generation of primes using sieve of Eratosthenes.
    if (common.ecm.primes[i >> 3] & (1 << (i & 7)))
    {     // Number i is prime.
      for (j = i * i; j < numPrimes; j += i)
      {   // Mark multiple of i as composite.
        common.ecm.primes[j >> 3] &= ~(1 << (j & 7));
      }
    }
  }
  // Let n = a^b +/- 1 (n = number to factor).
  // If -1<=n<=2 (mod p) does not hold, b cannot be multiple of p-1.
  // If -2<=n<=2 (mod p) does not hold, b cannot be multiple of (p-1)/2.
  for (i = 2; i < numPrimes; i++)
  {
    if (common.ecm.primes[i>>3] & (1 << (i & 7)))
    {      // i is prime according to sieve.
           // If n+/-1 is multiple of p, then it must be multiple
           // of p^2, otherwise it cannot be a perfect power.
      uint64_t remainder;
      int index;
      int rem = getRemainder(nbrToFactor, i);
      longToBigInteger(&Temp1, (uint64_t)i*(uint64_t)i);
      BigIntRemainder(nbrToFactor, &Temp1, &Temp2);     // Temp2 <- nbrToFactor % (i*i)
      remainder = (uint64_t)Temp2.limbs[0].x;
      if (rem == 1 || rem == i - 1)
      {
        if (Temp2.nbrLimbs > 1)
        {
          remainder += (uint64_t)Temp2.limbs[1].x << BITS_PER_GROUP;
        }
        // NumberFactor cannot be a power + 1 if condition holds.
        plus1 = (rem == 1 && remainder != 1);
        // NumberFactor cannot be a power - 1 if condition holds.
        minus1 = (rem == i - 1 && remainder != (uint64_t)i*(uint64_t)i - 1);
      }
      index = i / 2;
      if (!(common.ecm.ProcessExpon[index >> 3] & (1<<(index&7))))
      {
        continue;
      }
      modulus = remainder % i;
      if (modulus > (plus1 ? 1 : 2) && modulus < (minus1 ? i - 1 : i - 2))
      {
        for (j = index; j <= maxExpon; j += index)
        {
          common.ecm.ProcessExpon[j >> 3] &= ~(1 << (j&7));
        }
      }
      else
      {
        if (modulus == i - 2)
        {
          for (j = i - 1; j <= maxExpon; j += i - 1)
          {
            common.ecm.ProcessExpon[j >> 3] &= ~(1 << (j & 7));
          }
        }
      }
    }
  }
  for (j = 2; j < 100; j++)
  {
    double u = logar / log(j) + .000005;
    Exponent = (int)floor(u);
    if (u - Exponent > .00001)
      continue;
    if (Exponent % 3 == 0 && mod9 > 2 && mod9 < 7)
      continue;
    if (!(common.ecm.ProcessExpon[Exponent >> 3] & (1 << (Exponent & 7))))
      continue;
    if (ProcessExponent(pstFactors, nbrToFactor, Exponent))
      return;
  }
  for (; Exponent >= 2; Exponent--)
  {
    if (Exponent % 3 == 0 && mod9 > 2 && mod9 < 7)
    {
      continue;
    }
    if (!(common.ecm.ProcessExpon[Exponent >> 3] & (1 << (Exponent & 7))))
    {
      continue;
    }
    if (ProcessExponent(pstFactors, nbrToFactor, Exponent))
    {
      return;
    }
  }
}

// Perform Lehman algorithm
static void Lehman(BigInteger *nbr, int k, BigInteger *factor)
{
  int bitsSqrLow[] =
  {
    0x00000003, // 3
    0x00000013, // 5
    0x00000017, // 7
    0x0000023B, // 11
    0x0000161B, // 13
    0x0001A317, // 17
    0x00030AF3, // 19
    0x0005335F, // 23
    0x13D122F3, // 29
    0x121D47B7, // 31
    0x5E211E9B, // 37
    0x82B50737, // 41
    0x83A3EE53, // 43
    0x1B2753DF, // 47
    0x3303AED3, // 53
    0x3E7B92BB, // 59
    0x0A59F23B, // 61
  };
  int bitsSqrHigh[] =
  {
    0x00000000, // 3
    0x00000000, // 5
    0x00000000, // 7
    0x00000000, // 11
    0x00000000, // 13
    0x00000000, // 17
    0x00000000, // 19
    0x00000000, // 23
    0x00000000, // 29
    0x00000000, // 31
    0x00000016, // 37
    0x000001B3, // 41
    0x00000358, // 43
    0x00000435, // 47
    0x0012DD70, // 53
    0x022B6218, // 59
    0x1713E694, // 61
  };
  int primes[] = { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61 };
  int nbrs[17];
  int diffs[17];
  int i, j, m, r;
  int nbrLimbs, nbrIterations;
  BigInteger sqrRoot, nextroot;
  BigInteger a, c, sqr, val;
  if ((nbr->limbs[0].x & 1) == 0)
  { // nbr Even
    r = 0;
    m = 1;
  }
  else
  {
    if (k % 2 == 0)
    { // k Even
      r = 1;
      m = 2;
    }
    else
    { // k Odd
      r = (k + nbr->limbs[0].x) & 3;
      m = 4;
    }
  }
  intToBigInteger(&sqr, k<<2);
  BigIntMultiply(&sqr, nbr, &sqr);
  squareRoot(sqr.limbs, sqrRoot.limbs, sqr.nbrLimbs, &sqrRoot.nbrLimbs);
  sqrRoot.sign = SIGN_POSITIVE;
  CopyBigInt(&a, &sqrRoot);
  for (;;)
  {
    if ((a.limbs[0].x & (m-1)) == r)
    {
      BigIntMultiply(&a, &a, &nextroot);
      BigIntSubt(&nextroot, &sqr, &nextroot);
      if (nextroot.sign == SIGN_POSITIVE)
      {
        break;
      }
    }
    addbigint(&a, 1);                         // a <- a + 1
  }
  BigIntMultiply(&a, &a, &nextroot);
  BigIntSubt(&nextroot, &sqr, &c);
  for (i = 0; i < 17; i++)
  {
    int pr = primes[i];
    nbrs[i] = getRemainder(&c, pr);    // nbrs[i] <- c % primes[i]
    diffs[i] = m * (getRemainder(&a, pr) * 2 + m) % pr;
  }
  nbrLimbs = factor->nbrLimbs;
  if (nbrLimbs > 10)
  {
    nbrIterations = 10000;
  }
  else
  {
    nbrIterations = 1000 * factor->nbrLimbs;
  }
  for (j = 0; j < nbrIterations; j++)
  {
    for (i = 0; i < 17; i++)
    {
      int shiftBits = nbrs[i];
      if (shiftBits < 32)
      {
        if ((bitsSqrLow[i] & (1 << shiftBits)) == 0)
        { // Not a perfect square
          break;
        }
      }
      else if ((bitsSqrHigh[i] & (1 << (shiftBits-32))) == 0)
      { // Not a perfect square
        break;
      }
    }
    if (i == 17)
    { // Test for perfect square
      intToBigInteger(&c, m * j);           // c <- m * j
      BigIntAdd(&a, &c, &val);
      BigIntMultiply(&val, &val, &c);       // c <- val * val
      BigIntSubt(&c, &sqr, &c);             // c <- val * val - sqr
      squareRoot(c.limbs, sqrRoot.limbs, c.nbrLimbs, &sqrRoot.nbrLimbs);
      sqrRoot.sign = SIGN_POSITIVE;         // sqrRoot <- sqrt(c)
      BigIntAdd(&sqrRoot, &val, &sqrRoot);
      BigIntGcd(&sqrRoot, nbr, &c);         // Get GCD(sqrRoot + val, nbr)
      if (c.nbrLimbs > 1)
      {    // Non-trivial factor has been found.
        CopyBigInt(factor, &c);
        return;
      }
    }
    for (i = 0; i < 17; i++)
    {
      nbrs[i] = (nbrs[i] + diffs[i]) % primes[i];
      diffs[i] = (diffs[i] + 2 * m * m) % primes[i];
    }
  }
  intToBigInteger(factor, 1);   // Factor not found.
}

static enum eEcmResult ecmCurve(BigInteger *N)
{
  BigInteger potentialFactor;
#ifdef __EMSCRIPTEN__
  char text[20];
#endif
  EC %= 50000000;   // Convert to curve number.
  for (;;)
  {
#ifdef __EMSCRIPTEN__
    char *ptrText;
#endif
    int I, Pass;
    int i, j, u;
    long long L1, L2, LS, P, IP, Paux = 1;
    if (NextEC > 0)
    {
      EC = NextEC;
      NextEC = -1;
      if (EC >= TYP_SIQS)
      {
        common.ecm.GD[0].x = 1;   // Set GD to 1.
        memset(&common.ecm.GD[1], 0, (NumberLength - 1) * sizeof(limb));
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
      L1 = NumberLength*9;          // Get number of digits.
      if (NextEC == 0 && L1 <= 110) // Force switch to SIQS and number not too large.
      {
        EC += TYP_SIQS;
        return CHANGE_TO_SIQS;
      }
      if (L1 > 30 && L1 <= 90)      // If between 30 and 90 digits...         
      {                             // Switch to SIQS.
        int limit = limits[((int)L1 - 31) / 5];
        if (EC % 50000000 >= limit)
        {                           // Switch to SIQS.
          EC += TYP_SIQS;
          return CHANGE_TO_SIQS;
        }
      }
    }
    // Try to factor BigInteger N using Lehman algorithm. Result in potentialFactor.
    Lehman(N, EC % 50000000, &potentialFactor);
    if (potentialFactor.nbrLimbs > 1)
    {                // Factor found.
      memcpy(common.ecm.GD, potentialFactor.limbs, NumberLength * sizeof(limb));
      memset(&common.ecm.GD[potentialFactor.nbrLimbs], 0, 
             (NumberLength - potentialFactor.nbrLimbs) * sizeof(limb));
      foundByLehman = TRUE;
      return FACTOR_FOUND;
    }
    L1 = 2000;
    L2 = 200000;
    LS = 45;
    Paux = EC;
    nbrPrimes = 303; /* Number of primes less than 2000 */
    if (EC > 25)
    {
      if (EC < 326)
      {
        L1 = 50000;
        L2 = 5000000;
        LS = 224;
        Paux = EC - 24;
        nbrPrimes = 5133; /* Number of primes less than 50000 */
      }
      else
      {
        if (EC < 2000)
        {
          L1 = 1000000;
          L2 = 100000000;
          LS = 1001;
          Paux = EC - 299;
          nbrPrimes = 78498; /* Number of primes less than 1000000 */
        }
        else
        {
          L1 = 11000000;
          L2 = 1100000000;
          LS = 3316;
          Paux = EC - 1900;
          nbrPrimes = 726517; /* Number of primes less than 11000000 */
        }
      }
    }
#ifdef __EMSCRIPTEN__
    ptrText = ptrLowerText;  // Point after number that is being factored.
    strcpy(ptrText, lang ? "<p>Curva ": "<p>Curve ");
    ptrText += strlen(ptrText);
    int2dec(&ptrText, EC);   // Show curve number.
    strcpy(ptrText, lang?" usando límites B1=": " using bounds B1=");
    ptrText += strlen(ptrText);
    int2dec(&ptrText, L1);   // Show first bound.
    strcpy(ptrText, lang? " y B2=" : " and B2=");
    ptrText += strlen(ptrText);
    int2dec(&ptrText, L2);   // Show second bound.
    strcpy(ptrText, "</p>");
    ptrText += strlen(ptrText);
    databack(lowerText);
#if 0
    primalityString =
      textAreaContents
      + StringToLabel
      + "\nLimit (B1="
      + L1
      + "; B2="
      + L2
      + ")    Curve ";
    UpperLine = "Digits in factor:   ";
    LowerLine = "Probability:        ";
    for (I = 0; I < 6; I++)
    {
      UpperLine += "    >= " + (I * 5 + 15);
      Prob =
        (int)Math.round(
          100
          * (1 - exp(-((double)L1 * (double)Paux) / ProbArray[I])));
      if (Prob == 100)
      {
        LowerLine += "    100% ";
      }
      else
      {
        if (Prob >= 10)
        {
          LowerLine += "     " + Prob + "% ";
        }
        else
        {
          LowerLine += "      " + Prob + "% ";
        }
      }
    } /* end for */
    lowerTextArea.setText(
      primalityString + EC + "\n" + UpperLine + "\n" + LowerLine);
#endif
#endif

    //  Compute A0 <- 2 * (EC+1)*modinv(3 * (EC+1) ^ 2 - 1, N) mod N
                                               // Aux2 <- 1 in Montgomery notation.
    memcpy(common.ecm.Aux2, MontgomeryMultR1, NumberLength * sizeof(limb));                                               
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
    if (BigNbrIsZero(common.ecm.Aux3))
    {
      continue;
    }
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
    /**************/
    /* First step */
    /**************/
    memcpy(common.ecm.Xaux, common.ecm.X, NumberLength * sizeof(limb));
    memcpy(common.ecm.Zaux, common.ecm.Z, NumberLength * sizeof(limb));
    memcpy(common.ecm.GcdAccumulated, MontgomeryMultR1, (NumberLength+1) * sizeof(limb));
    for (Pass = 0; Pass < 2; Pass++)
    {
      /* For powers of 2 */
      indexPrimes = 0;
      StepECM = 1;
      for (I = 1; I <= L1; I <<= 1)
      {
        duplicate(common.ecm.X, common.ecm.Z, common.ecm.X, common.ecm.Z);
      }
      for (I = 3; I <= L1; I *= 3)
      {
        duplicate(common.ecm.W1, common.ecm.W2, common.ecm.X, common.ecm.Z);
        add3(common.ecm.X, common.ecm.Z, common.ecm.X, common.ecm.Z, common.ecm.W1, common.ecm.W2, common.ecm.X, common.ecm.Z);
      }

      if (Pass == 0)
      {
        modmult(common.ecm.GcdAccumulated, common.ecm.Z, common.ecm.Aux1);
        memcpy(common.ecm.GcdAccumulated, common.ecm.Aux1, NumberLength * sizeof(limb));
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
        indexPrimes++;
        P = SmallPrime[indexM];
        for (IP = P; IP <= L1; IP *= P)
        {
          prac((int)P, common.ecm.X, common.ecm.Z, common.ecm.W1, common.ecm.W2, common.ecm.W3, common.ecm.W4);
        }
        indexM++;
        if (Pass == 0)
        {
          modmult(common.ecm.GcdAccumulated, common.ecm.Z, common.ecm.Aux1);
          memcpy(common.ecm.GcdAccumulated, common.ecm.Aux1, NumberLength * sizeof(limb));
        }
        else
        {
          if (gcdIsOne(common.ecm.Z) > 1)
          {
            return FACTOR_FOUND;
          }
        }
      } while (SmallPrime[indexM - 1] <= LS);
      P += 2;

      /* Initialize sieve2310[n]: 1 if gcd(P+2n,2310) > 1, 0 otherwise */
      u = (int)P;
      for (i = 0; i < SIEVE_SIZE; i++)
      {
        common.ecm.sieve2310[i] =
          (u % 3 == 0
            || u % 5 == 0
            || u % 7 == 0
#if MAX_PRIME_SIEVE == 11
            || u % 11 == 0
#endif
            ? (unsigned char)1 : (unsigned char)0);
        u += 2;
      }
      do
      {
        /* Generate sieve */
        GenerateSieve((int)P);

        /* Walk through sieve */

        for (i = 0; i < 10*SIEVE_SIZE; i++)
        {
          if (common.ecm.sieve[i] != 0)
          {
            continue; /* Do not process composites */
          }
          if (P + 2 * i > L1)
          {
            break;
          }
          indexPrimes++;
          prac((int)(P + 2 * i), common.ecm.X, common.ecm.Z, common.ecm.W1, common.ecm.W2, common.ecm.W3, common.ecm.W4);
          if (Pass == 0)
          {
            modmult(common.ecm.GcdAccumulated, common.ecm.Z, common.ecm.Aux1);
            memcpy(common.ecm.GcdAccumulated, common.ecm.Aux1, NumberLength * sizeof(limb));
          }
          else
          {
            if (gcdIsOne(common.ecm.Z) > 1)
            {
              return FACTOR_FOUND;
            }
          }
        }
        P += 20*SIEVE_SIZE;
      } while (P < L1);
      if (Pass == 0)
      {
        if (BigNbrIsZero(common.ecm.GcdAccumulated))
        { // If GcdAccumulated is
          memcpy(common.ecm.X, common.ecm.Xaux, NumberLength * sizeof(limb));
          memcpy(common.ecm.Z, common.ecm.Zaux, NumberLength * sizeof(limb));
          continue; // multiple of TestNbr, continue.
        }
        if (gcdIsOne(common.ecm.GcdAccumulated) > 1)
        {
          return FACTOR_FOUND;
        }
        break;
      }
    } /* end for Pass */

      /******************************************************/
      /* Second step (using improved standard continuation) */
      /******************************************************/
    StepECM = 2;
    j = 0;
    for (u = 1; u < SIEVE_SIZE; u += 2)
    {
      if (u % 3 == 0 || u % 5 == 0 || u % 7 == 0
#if MAX_PRIME_SIEVE == 11
        || u % 11 == 0
#endif
        )
      {
        common.ecm.sieve2310[u / 2] = (unsigned char)1;
      }
      else
      {
        common.ecm.sieve2310[(common.ecm.sieveidx[j++] = u / 2)] = (unsigned char)0;
      }
    }
    memcpy(&common.ecm.sieve2310[HALF_SIEVE_SIZE], &common.ecm.sieve2310[0], HALF_SIEVE_SIZE);
    memcpy(common.ecm.Xaux, common.ecm.X, NumberLength * sizeof(limb));  // (X:Z) -> Q (output
    memcpy(common.ecm.Zaux, common.ecm.Z, NumberLength * sizeof(limb));  //         from step 1)
    for (Pass = 0; Pass < 2; Pass++)
    {
      int Qaux, J;
      memcpy(common.ecm.GcdAccumulated, MontgomeryMultR1, NumberLength * sizeof(limb));
      memcpy(common.ecm.UX, common.ecm.X, NumberLength * sizeof(limb));
      memcpy(common.ecm.UZ, common.ecm.Z, NumberLength * sizeof(limb));  // (UX:UZ) -> Q 
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
        memcpy(common.ecm.WX, common.ecm.X, NumberLength * sizeof(limb));
        memcpy(common.ecm.WZ, common.ecm.Z, NumberLength * sizeof(limb));
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
        if (Pass == 0)
        {
          modmult(common.ecm.GcdAccumulated, common.ecm.Aux1, common.ecm.Aux2);
          memcpy(common.ecm.GcdAccumulated, common.ecm.Aux2, NumberLength * sizeof(limb));
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
          memcpy(common.ecm.DX, common.ecm.X, NumberLength * sizeof(limb));
          memcpy(common.ecm.DZ, common.ecm.Z, NumberLength * sizeof(limb));  // (DX:DZ) -> HALF_SIEVE_SIZE*Q
        }
        if (I % 3 != 0 && I % 5 != 0 && I % 7 != 0
#if MAX_PRIME_SIEVE == 11
          && I % 11 != 0
#endif
          )
        {
          J++;
          ModInvBigNbr(common.ecm.Z, common.ecm.Aux1, TestNbr, NumberLength);
          modmult(common.ecm.Aux1, common.ecm.X, common.ecm.root[J]); // root[J] <- X/Z
        }
        memcpy(common.ecm.UX, common.ecm.WX, NumberLength * sizeof(limb));  // (UX:UZ) <-
        memcpy(common.ecm.UZ, common.ecm.WZ, NumberLength * sizeof(limb));  // Previous (X:Z)
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
      memcpy(common.ecm.UX, common.ecm.X, NumberLength * sizeof(limb));
      memcpy(common.ecm.UZ, common.ecm.Z, NumberLength * sizeof(limb));    // (UX:UZ) -> SIEVE_SIZE*Q
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
      Qaux = (int)(L1 / (2*SIEVE_SIZE));
      maxIndexM = (int)(L2 / (2 * SIEVE_SIZE));
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
          if (indexM % 10 == 0 || indexM == Qaux)
          {
            GenerateSieve(indexM / 10 * (20*SIEVE_SIZE) + 1);
          }
          /* Walk through sieve */
          J = HALF_SIEVE_SIZE + (indexM % 10) * SIEVE_SIZE;
          for (i = 0; i < GROUP_SIZE; i++)
          {
            j = common.ecm.sieveidx[i]; // 0 < J < HALF_SIEVE_SIZE
            if (common.ecm.sieve[J + j] != 0 && common.ecm.sieve[J - 1 - j] != 0)
            {
              continue; // Do not process if both are composite numbers.
            }
            SubtBigNbrModN(common.ecm.Aux1, common.ecm.root[i], common.ecm.M, TestNbr, NumberLength);
            modmult(common.ecm.GcdAccumulated, common.ecm.M, common.ecm.Aux2);
            memcpy(common.ecm.GcdAccumulated, common.ecm.Aux2, NumberLength * sizeof(limb));
          }
          if (Pass != 0)
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
          memcpy(common.ecm.WX, common.ecm.X, NumberLength * sizeof(limb));
          memcpy(common.ecm.WZ, common.ecm.Z, NumberLength * sizeof(limb));
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
          memcpy(common.ecm.UX, common.ecm.WX, NumberLength * sizeof(limb));
          memcpy(common.ecm.UZ, common.ecm.WZ, NumberLength * sizeof(limb));
        }
      } // end for Q
      if (Pass == 0)
      {
        int rc;
        if (BigNbrIsZero(common.ecm.GcdAccumulated))
        { // If GcdAccumulated is zero
          memcpy(common.ecm.X, common.ecm.Xaux, NumberLength * sizeof(limb));
          memcpy(common.ecm.Z, common.ecm.Zaux, NumberLength * sizeof(limb));
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
        if (memcmp(common.ecm.GD, TestNbr, NumberLength*sizeof(limb)))
        {           // GCD is not 1 or TestNbr
          return FACTOR_FOUND;
        }
      }
    } /* end for Pass */
    performLehman = TRUE;
  }       /* End curve calculation */
}
#ifdef __EMSCRIPTEN__
char *ShowFactoredPart(BigInteger *pNbr, void *vFactors)
{
  struct sFactors *pstFactors = (struct sFactors *)vFactors;
  ptrLowerText = lowerText;
  *ptrLowerText++ = '3';
  if (vFactors != NULL && pstFactors->multiplicity > 1)
  {    // Some factorization known.
    int NumberLengthBak = NumberLength;
    strcpy(ptrLowerText, "<p class=\"blue\">");
    ptrLowerText += strlen(ptrLowerText);
    SendFactorizationToOutput(pstFactors, &ptrLowerText, 1);
    strcpy(ptrLowerText, "</p>");
    ptrLowerText += strlen(ptrLowerText);
    NumberLength = NumberLengthBak;
  }
  if (StepECM == 3)
  {
    strcpy(ptrLowerText, lang ? "<p>Comprobando si es primo " : "<p>Testing primality of ");
  }
  else
  {
    strcpy(ptrLowerText, lang ? "<p>Factorizando " : "<p>Factoring ");
  }
  ptrLowerText += strlen(ptrLowerText);
  if (hexadecimal)
  {
    Bin2Hex(pNbr->limbs, ptrLowerText, pNbr->nbrLimbs, groupLen);
  }
  else
  {
    Bin2Dec(pNbr->limbs, ptrLowerText, pNbr->nbrLimbs, groupLen);
  }
  ptrLowerText += strlen(ptrLowerText);
  strcpy(ptrLowerText, "</p>");
  ptrLowerText += strlen(ptrLowerText);
  return ptrLowerText;
}
void ShowLowerText(void)
{
  databack(lowerText);
}
#endif

static void ecm(BigInteger *N, struct sFactors *pstFactors)
{
  int P, Q;
#ifndef __EMSCRIPTEN__
  (void)pstFactors;     // Ignore parameter.
#endif
  common.ecm.fieldTX = common.ecm.TX;
  common.ecm.fieldTZ = common.ecm.TZ;
  common.ecm.fieldUX = common.ecm.UX;
  common.ecm.fieldUZ = common.ecm.UZ;
  //  int Prob;
  //  BigInteger NN;

  common.ecm.fieldAA = common.ecm.AA;
  NumberLength = N->nbrLimbs;
  memcpy(TestNbr, N->limbs, NumberLength * sizeof(limb));
  GetYieldFrequency();
  GetMontgomeryParms(NumberLength);
  memset(common.ecm.M, 0, NumberLength * sizeof(limb));
  memset(common.ecm.DX, 0, NumberLength * sizeof(limb));
  memset(common.ecm.DZ, 0, NumberLength * sizeof(limb));
  memset(common.ecm.W3, 0, NumberLength * sizeof(limb));
  memset(common.ecm.W4, 0, NumberLength * sizeof(limb));
  memset(common.ecm.GD, 0, NumberLength * sizeof(limb));
#ifdef __EMSCRIPTEN__
  ptrLowerText = ShowFactoredPart(N, pstFactors);
#endif
  EC--;
  SmallPrime[0] = 2;
  P = 3;
  indexM = 1;
  for (indexM = 1; indexM < sizeof(SmallPrime)/sizeof(SmallPrime[0]); indexM++)
  {     // Loop that fills the SmallPrime array.
    SmallPrime[indexM] = (int)P; /* Store prime */
    do
    {
      P += 2;
      for (Q = 3; Q * Q <= P; Q += 2)
      { /* Check if P is prime */
        if (P % Q == 0)
        {
          break;  /* Composite */
        }
      }
    } while (Q * Q <= P);
  }
  foundByLehman = FALSE;
  do
  {
    enum eEcmResult ecmResp = ecmCurve(N);
    if (ecmResp == CHANGE_TO_SIQS)
    {    // Perform SIQS
#ifdef __EMSCRIPTEN__
      double originalTenths = tenths();
      uint64_t oldModularMult = lModularMult;
#endif
      FactoringSIQS(TestNbr, common.ecm.GD);
#ifdef __EMSCRIPTEN__
      SIQSModMult += lModularMult - oldModularMult;
      timeSIQS += (int)(tenths() - originalTenths);
      nbrSIQS++;
#endif

      break;
    }
    else if (ecmResp == FACTOR_FOUND)
    {
      break;
    }
  } while (!memcmp(common.ecm.GD, TestNbr, NumberLength*sizeof(limb)));
  StepECM = 0; /* do not show pass number on screen */
}

void SendFactorizationToOutput(struct sFactors *pstFactors, char **pptrOutput, int doFactorization)
{
  char *ptrOutput = *pptrOutput;
  strcpy(ptrOutput, tofactorDec);
  ptrOutput += strlen(ptrOutput);
  if (doFactorization)
  {
    struct sFactors *pstFactor;
    pstFactor = pstFactors+1;
    if (tofactor.sign == SIGN_POSITIVE && pstFactors->multiplicity == 1 && pstFactor->multiplicity == 1 &&
      (*pstFactor->ptrFactor > 1 || *(pstFactor->ptrFactor + 1) > 1))
    {    // Do not show zero or one as prime.
      strcpy(ptrOutput, lang ? " es primo" : " is prime");
      ptrOutput += strlen(ptrOutput);
    }
    else
    {
      int i = 0;
      strcpy(ptrOutput, " = ");
      ptrOutput += strlen(ptrOutput);
      if (tofactor.sign == SIGN_NEGATIVE)
      {
        *ptrOutput++ = '-';
        if (tofactor.nbrLimbs > 1 || tofactor.limbs[0].x > 1)
        {
          if (prettyprint)
          {
            strcpy(ptrOutput, "1 &times; ");
          }
          else
          {
            strcpy(ptrOutput, "1 * ");
          }
          ptrOutput += strlen(ptrOutput);
        }
      }
      for (;;)
      {
        NumberLength = *pstFactor->ptrFactor;
        IntArray2BigInteger(pstFactor->ptrFactor, &factorValue);
        if (hexadecimal)
        {
          Bin2Hex(factorValue.limbs, ptrOutput, factorValue.nbrLimbs, groupLen);
        }
        else
        {
          Bin2Dec(factorValue.limbs, ptrOutput, factorValue.nbrLimbs, groupLen);
        }
        ptrOutput += strlen(ptrOutput);
        if (pstFactor->multiplicity > 1)
        {
          if (prettyprint)
          {
            strcpy(ptrOutput, "<sup>");
            ptrOutput += strlen(ptrOutput);
            int2dec(&ptrOutput, pstFactor->multiplicity);
            strcpy(ptrOutput, "</sup>");
            ptrOutput += strlen(ptrOutput);
          }
          else
          {
            *ptrOutput++ = '^';
            int2dec(&ptrOutput, pstFactor->multiplicity);
          }
        }
#ifdef ENABLE_VERBOSE
        int type = pstFactor->type;
        int isPrime = (pstFactor->upperBound == 0);
        if (type > 0)
        {
          int compositeType = type / 50000000 * 50000000;
          strcpy(ptrOutput, " <span class=\"verbose\">(");
          ptrOutput += strlen(ptrOutput);
          if (compositeType == TYP_AURIF)
          {
            strcpy(ptrOutput, "Aurifeuille");
            ptrOutput += strlen(ptrOutput);
            if (!isPrime)
            {
              strcpy(ptrOutput, lang ? " - Compuesto" : " - Composite");
            }
          }
          else if (compositeType == TYP_TABLE)
          {
            strcpy(ptrOutput, lang ? "Tabla" : "Table");
            ptrOutput += strlen(ptrOutput);
            if (!isPrime)
            {
              strcpy(ptrOutput, lang ? " - Compuesto" : " - Composite");
            }
          }
          else if (compositeType == TYP_SIQS)
          {
            strcpy(ptrOutput, lang? "<abbr title=\"Criba cuadrática autoinicializada\">SIQS</abbr>":
                                    "<abbr title=\"Self-Initializing Quadratic Sieve\">SIQS</abbr>");
            ptrOutput += strlen(ptrOutput);
            if (!isPrime)
            {
              strcpy(ptrOutput, lang ? " - Compuesto" : " - Composite");
            }
          }
          else if (compositeType == TYP_LEHMAN)
          {
            strcpy(ptrOutput, "Lehman");
            ptrOutput += strlen(ptrOutput);
            if (!isPrime)
            {
              strcpy(ptrOutput, lang ? " - Compuesto" : " - Composite");
            }
          }
          else if (compositeType == TYP_RABIN)
          {
            strcpy(ptrOutput, lang ? "Miller y Rabin" : "Miller &amp; Rabin");
            ptrOutput += strlen(ptrOutput);
            if (!isPrime)
            {
              strcpy(ptrOutput, lang ? " - Compuesto" : " - Composite");
            }
          }
          else if (compositeType == TYP_DIVISION)
          {
            strcpy(ptrOutput, lang ? "División" : "Division");
            ptrOutput += strlen(ptrOutput);
            if (!isPrime)
            {
              strcpy(ptrOutput, lang ? " - Compuesto" : " - Composite");
            }
          }
          else if (type > TYP_EC)
          {
            if (isPrime)
            {
              strcpy(ptrOutput, lang ? "<abbr title=\"Método de curvas elípticas\">ECM</abbr>, curva " :
                "<abbr title=\"Elliptic curve method\">ECM</abbr>, curve ");
              ptrOutput += strlen(ptrOutput);
              int2dec(&ptrOutput, type - TYP_EC);
              *ptrOutput = 0;     // Add string terminator.
            }
            else
            {
              strcpy(ptrOutput, lang ? "Compuesto" : "Composite");
            }
          }
          else if (!isPrime)
          {
            strcpy(ptrOutput, lang ? "Compuesto": "Composite");
          }
          ptrOutput += strlen(ptrOutput);
          strcpy(ptrOutput, ")</span>");
          ptrOutput += strlen(ptrOutput);
        }
        else if (!isPrime)
        {
          strcpy(ptrOutput, "<span class=\"terse\"> (");
          ptrOutput += strlen(ptrOutput);
          strcpy(ptrOutput, lang ? "Compuesto" : "Composite");
          ptrOutput += strlen(ptrOutput);
          strcpy(ptrOutput, ")</span>");
          ptrOutput += strlen(ptrOutput);
        }
        if (type < 0)
        {
          strcpy(ptrOutput, " (Unknown)");
          ptrOutput += strlen(ptrOutput);
        }
#endif
        if (++i == pstFactors->multiplicity)
        {
          break;
        }
        if (prettyprint)
        {
          strcpy(ptrOutput, " &times; ");
        }
        else
        {
          strcpy(ptrOutput, " * ");
        }
        ptrOutput += strlen(ptrOutput);
        pstFactor++;
      }
    }
  }
  *pptrOutput = ptrOutput;
}

static void SortFactors(struct sFactors *pstFactors)
{
  int factorNumber, factorNumber2, ctr;
  struct sFactors *pstCurFactor = pstFactors + 1;
  struct sFactors stTempFactor;
  int *ptrNewFactor;
  struct sFactors *pstNewFactor;
  for (factorNumber = 1; factorNumber <= pstFactors->multiplicity; factorNumber++, pstCurFactor++)
  {
    pstNewFactor = pstCurFactor + 1;
    for (factorNumber2 = factorNumber + 1; factorNumber2 <= pstFactors->multiplicity; factorNumber2++, pstNewFactor++)
    {
      int *ptrFactor = pstCurFactor->ptrFactor;
      int *ptrFactor2 = pstNewFactor->ptrFactor;
      if (*ptrFactor < *ptrFactor2)
      {     // Factors already in correct order.
        continue;
      }
      if (*ptrFactor == *ptrFactor2)
      {
        for (ctr = *ptrFactor; ctr > 1; ctr--)
        {
          if (*(ptrFactor + ctr) != *(ptrFactor2 + ctr))
          {
            break;
          }
        }
        if (*(ptrFactor + ctr) < *(ptrFactor2 + ctr))
        {     // Factors already in correct order.
          continue;
        }
        if (*(ptrFactor + ctr) == *(ptrFactor2 + ctr))
        {     // Factors are the same.
          pstCurFactor->multiplicity += pstNewFactor->multiplicity;
          ctr = pstFactors->multiplicity - factorNumber2;
          if (ctr > 0)
          {
            memmove(pstNewFactor, pstNewFactor + 1, ctr * sizeof(struct sFactors));
          }
          pstFactors->multiplicity--;   // Indicate one less known factor.
          continue;
        }
      }
      // Exchange both factors.
      memcpy(&stTempFactor, pstCurFactor, sizeof(struct sFactors));
      memcpy(pstCurFactor, pstNewFactor, sizeof(struct sFactors));
      memcpy(pstNewFactor, &stTempFactor, sizeof(struct sFactors));
    }
  }
  // Find location for new factors.
  ptrNewFactor = 0;
  pstCurFactor = pstFactors + 1;
  for (factorNumber = 1; factorNumber <= pstFactors->multiplicity; factorNumber++, pstCurFactor++)
  {
    int *ptrPotentialNewFactor = pstCurFactor->ptrFactor + *(pstCurFactor->ptrFactor) + 1;
    if (ptrPotentialNewFactor > ptrNewFactor)
    {
      ptrNewFactor = ptrPotentialNewFactor;
    }
  }
  pstFactors->ptrFactor = ptrNewFactor;
}

// Insert new factor found into factor array. This factor array must be sorted.
static void insertIntFactor(struct sFactors *pstFactors, struct sFactors *pstFactorDividend, int divisor, 
  int expon, BigInteger *cofactor)
{
  struct sFactors *pstCurFactor;
  int multiplicity;
  int factorNumber;
  int *ptrFactor = pstFactorDividend->ptrFactor;
  int nbrLimbs = *ptrFactor;
  int *ptrValue;
  pstFactorDividend->upperBound = divisor;
  // Divide number by factor just found.
  if (cofactor == NULL)
  {        // Find cofactor.
    DivBigNbrByInt(ptrFactor + 1, divisor, ptrFactor + 1, nbrLimbs);
    if (*(ptrFactor + nbrLimbs) == 0)
    {
      (*ptrFactor)--;
    }
  }
  else
  {        // Cofactor given as a parameter.
    NumberLength = cofactor->nbrLimbs;
    BigInteger2IntArray(ptrFactor, cofactor);
  }
  // Check whether prime is already in factor list.
  pstCurFactor = pstFactors+1;
  for (factorNumber = 1; factorNumber <= pstFactors->multiplicity; factorNumber++, pstCurFactor++)
  {
    ptrValue = pstCurFactor->ptrFactor;  // Point to factor in factor array.
    if (*ptrValue == 1 && *(ptrValue+1) == divisor)
    {  // Prime already found: increment multiplicity and go out.
      pstCurFactor->multiplicity += pstFactorDividend->multiplicity * expon;
      ptrValue = pstFactorDividend->ptrFactor;
      if (*ptrValue == 1 && *(ptrValue + 1) == 1)
      {    // Dividend is 1 now so discard it.
        *pstFactorDividend = *(pstFactors + pstFactors->multiplicity--);
      }
      SortFactors(pstFactors);
      return;
    }
    if (*ptrValue > 1 || *(ptrValue + 1) > divisor)
    {   // Factor in factor list is greater than factor to insert. Exit loop.
      break;
    }
  }
  pstFactors->multiplicity++; // Indicate new known factor.
  multiplicity = pstFactorDividend->multiplicity;
  // Move all elements.
  ptrValue = pstFactorDividend->ptrFactor;
  if (pstFactors->multiplicity > factorNumber)
  {
    memmove(pstCurFactor + 1, pstCurFactor,
      (pstFactors->multiplicity - factorNumber) * sizeof(struct sFactors));
  }
  if (*ptrValue == 1 && *(ptrValue + 1) == 1)
  {
    pstCurFactor = pstFactorDividend;
  }
  else
  {
    ptrValue = pstFactors->ptrFactor;
    pstCurFactor->ptrFactor = ptrValue;
    pstCurFactor->multiplicity = multiplicity * expon;
    pstFactors->ptrFactor += 2;  // Next free memory.
  }
  pstCurFactor->upperBound = 0;
  *ptrValue = 1;  // Number of limbs.
  *(ptrValue + 1) = divisor;
  SortFactors(pstFactors);
}

// Insert new factor found into factor array. This factor array must be sorted.
// The divisor must be also sorted.
static void insertBigFactor(struct sFactors *pstFactors, BigInteger *divisor, int type)
{
  struct sFactors *pstCurFactor;
  int factorNumber;
  int lastFactorNumber = pstFactors->multiplicity;
  struct sFactors *pstNewFactor = pstFactors + lastFactorNumber + 1;
  int *ptrNewFactorLimbs = pstFactors->ptrFactor;
  pstCurFactor = pstFactors + 1;
  for (factorNumber = 1; factorNumber <= lastFactorNumber; factorNumber++, pstCurFactor++)
  {     // For each known factor...
    int *ptrFactor = pstCurFactor->ptrFactor;
    NumberLength = *ptrFactor;
    IntArray2BigInteger(ptrFactor, &Temp2);    // Convert known factor to Big Integer.
    BigIntGcd(divisor, &Temp2, &Temp3);         // Temp3 is the GCD between known factor and divisor.
    if (Temp3.nbrLimbs == 1 && Temp3.limbs[0].x < 2)
    {                                           // divisor is not a new factor (GCD = 0 or 1).
      continue;
    }
    if (TestBigNbrEqual(&Temp2, &Temp3))
    {                                           // GCD is equal to known factor.
      continue;
    }
    // At this moment both GCD and known factor / GCD are new known factors. Replace the known factor by
    // known factor / GCD and generate a new known factor entry.
    NumberLength = Temp3.nbrLimbs;
    BigInteger2IntArray(ptrNewFactorLimbs, &Temp3);      // Append new known factor.
    BigIntDivide(&Temp2, &Temp3, &Temp4);                // Divide by this factor.
    NumberLength = Temp4.nbrLimbs;
    BigInteger2IntArray(ptrFactor, &Temp4);              // Overwrite old known factor.
    pstNewFactor->multiplicity = pstCurFactor->multiplicity;
    pstNewFactor->ptrFactor = ptrNewFactorLimbs;
    pstNewFactor->upperBound = pstCurFactor->upperBound;
    if (type < 50000000)
    {          // Factor found using ECM.
      pstNewFactor->type = TYP_EC + EC;
      type = pstCurFactor->type / 50000000 * 50000000;
      if (type == 0)
      {
        pstCurFactor->type = TYP_DIVISION + EC;
      }
      else
      {
        pstCurFactor->type = type + EC;
      }
    }
    else
    {          // Found otherwise.
      pstNewFactor->type = type;
    }
    pstNewFactor++;
    pstFactors->multiplicity++;
    ptrNewFactorLimbs += 1 + Temp3.nbrLimbs;
  }
  // Sort factors in ascending order. If two factors are equal, coalesce them.
  // Divide number by factor just found.
  SortFactors(pstFactors);
}
#ifdef __EMSCRIPTEN__
static void showECMStatus(void)
{
  char status[200];
  int elapsedTime;
  char *ptrStatus;
  if ((lModularMult % yieldFreq) != 0)
  {
    return;
  }
  elapsedTime = (int)(tenths() - originalTenthSecond);
  if (elapsedTime / 10 == oldTimeElapsed / 10)
  {
    return;
  }
  oldTimeElapsed = elapsedTime;
  ptrStatus = status;
  strcpy(ptrStatus, lang ? "4<p>Transcurrió " : "4<p>Time elapsed: ");
  ptrStatus += strlen(ptrStatus);
  GetDHMS(&ptrStatus, elapsedTime / 10);
  strcpy(ptrStatus, "&nbsp;&nbsp;&nbsp;");  // Separate with three spaces.
  ptrStatus += strlen(ptrStatus);
  switch (StepECM)
  {
  case 1:
    strcpy(ptrStatus, lang ? "Paso 1: " : "Step 1: ");
    ptrStatus += strlen(ptrStatus);
    int2dec(&ptrStatus, indexPrimes / (nbrPrimes / 100));
    *ptrStatus++ = '%';
    break;
  case 2:
    strcpy(ptrStatus, lang ? "Paso 2: " : "Step 2: ");
    ptrStatus += strlen(ptrStatus);
    int2dec(&ptrStatus, maxIndexM == 0 ? 0 : indexM / (maxIndexM / 100));
    *ptrStatus++ = '%';
    break;
  case 3:
    strcpy(ptrStatus, lang ? "Progreso: " : "Progress: ");
    ptrStatus += strlen(ptrStatus);
    int2dec(&ptrStatus, percentageBPSW);
    *ptrStatus++ = '%';
    break;
  }
  strcpy(ptrStatus, "</p>");
  databack(status);
}

// Save factors in Web Storage so factorization can continue the next time the application runs.
static void SaveFactors(struct sFactors *pstFactors)
{
#ifdef FACTORIZATION_APP
  struct sFactors *pstCurFactor = pstFactors + 1;
  int factorNbr;
  BigInteger bigint;
  char *ptrText = common.saveFactors.text;
  if (oldNbrFactors == pstFactors->multiplicity)
  {
    return;
  }
  oldNbrFactors = pstFactors->multiplicity;
  *ptrText++ = '8';
  strcpy(ptrText, ptrInputText);
  ptrText += strlen(ptrText);
  *ptrText++ = '=';
  for (factorNbr = 1; factorNbr <= pstFactors->multiplicity; factorNbr++, pstCurFactor++)
  {
    if (factorNbr > 1)
    {
      *ptrText++ = '*';
    }
    NumberLength = *pstCurFactor->ptrFactor;
    IntArray2BigInteger(pstCurFactor->ptrFactor, &bigint);
    BigInteger2Dec(&bigint, ptrText, -100000);   // Factors are saved in decimal.
    ptrText += strlen(ptrText);
    *ptrText++ = '^';
    int2dec(&ptrText, pstCurFactor->multiplicity);
    *ptrText++ = '(';
    int2dec(&ptrText, pstCurFactor->upperBound);
    *ptrText++ = ',';
    int2dec(&ptrText, pstCurFactor->type);
    *ptrText++ = ')';
  }
  *ptrText++ = 0;
  databack(common.saveFactors.text);
#endif
}

#endif

static char *findChar(char *str, char c)
{
  while (*str != 0)
  {
    if (*str == c)
    {
      return str;
    }
    str++;
  }
  return NULL;
}

static int getNextInteger(char **ppcFactors, int *result, char delimiter)
{
  char *pcFactors = *ppcFactors;
  char *ptrCharFound = findChar(pcFactors, delimiter);
  if (ptrCharFound == NULL)
  {
    return 1;
  }
  *ptrCharFound = 0;
  Dec2Bin(pcFactors, Temp1.limbs, (int)(ptrCharFound - pcFactors), &Temp1.nbrLimbs);
  if (Temp1.nbrLimbs != 1)
  {   // Exponent is too large.
    return 1;
  }
  *result = (int)Temp1.limbs[0].x;
  *ppcFactors = ptrCharFound+1;
  return 0;
}

// Return: 0 = No factors found.
//         1 = Factors found.
// Use: Xaux for square root of -1.
//      Zaux for square root of 1.
static int factorCarmichael(BigInteger *pValue, struct sFactors *pstFactors)
{
  int randomBase = 0;
#ifdef __EMSCRIPTEN__
  timePrimalityTests = 0;  // Reset to zero all timings.
  timeSIQS = 0;
  timeECM = 0;
  nbrPrimalityTests = 0;   // Reset to zero all counters.
  nbrSIQS = 0;
  nbrECM = 0;
  SIQSModMult = 0;
#endif
  int factorsFound = FALSE;
  int nbrLimbsQ, countdown, ctr;
  int nbrLimbs = pValue->nbrLimbs;
  int sqrtOneFound = FALSE;
  int sqrtMinusOneFound = FALSE;
  int Aux1Len;
  limb *pValueLimbs = pValue->limbs;
  (pValueLimbs + nbrLimbs)->x = 0;
  memcpy(q, pValueLimbs, (nbrLimbs + 1) * sizeof(limb));
  nbrLimbsQ = nbrLimbs;
  q[0]--;                     // q = p - 1 (p is odd, so there is no carry).
  memcpy(common.ecm.Aux1, q, (nbrLimbsQ + 1) * sizeof(q[0]));
  Aux1Len = nbrLimbs;
  DivideBigNbrByMaxPowerOf2(&ctr, common.ecm.Aux1, &Aux1Len);
  memcpy(TestNbr, pValueLimbs, nbrLimbs * sizeof(limb));
  TestNbr[nbrLimbs].x = 0;
  GetMontgomeryParms(nbrLimbs);
  for (countdown = 20; countdown > 0; countdown--)
  {
    int i;
    NumberLength = nbrLimbs;
    randomBase = (int)((uint64_t)randomBase * 89547121 + 1762281733) & MAX_INT_NBR;
    modPowBaseInt(randomBase, common.ecm.Aux1, Aux1Len, common.ecm.Aux2); // Aux2 = base^Aux1.
                                                 // If Mult1 = 1 or Mult1 = TestNbr-1, then try next base.
    if (checkOne(common.ecm.Aux2, nbrLimbs) != 0 || checkMinusOne(common.ecm.Aux2, nbrLimbs) != 0)
    {
      continue;    // This base cannot find a factor. Try another one.
    }
    for (i = 0; i < ctr; i++)
    {              // Loop that squares number.
      modmult(common.ecm.Aux2, common.ecm.Aux2, common.ecm.Aux3);
      if (checkOne(common.ecm.Aux3, nbrLimbs) != 0)
      {            // Non-trivial square root of 1 found.
        if (!sqrtOneFound)
        {          // Save it to perform GCD later.
          memcpy(common.ecm.Zaux, common.ecm.Aux2, nbrLimbs*sizeof(limb));
          sqrtOneFound = TRUE;
        }
        else
        {          // Try to find non-trivial factor by doing GCD.
          SubtBigNbrMod(common.ecm.Aux2, common.ecm.Zaux, common.ecm.Aux4);
          UncompressLimbsBigInteger(common.ecm.Aux4, &Temp2);
          BigIntGcd(pValue, &Temp2, &Temp4);
          if ((Temp4.nbrLimbs != 1 || Temp4.limbs[0].x > 1) &&
            (Temp4.nbrLimbs != NumberLength ||
              memcmp(pValue->limbs, Temp4.limbs, NumberLength * sizeof(limb))))
          {          // Non-trivial factor found.
            insertBigFactor(pstFactors, &Temp4, TYP_RABIN);
            factorsFound = TRUE;
          }
        }
                   // Try to find non-trivial factor by doing GCD.
        NumberLength = nbrLimbs;
        AddBigNbrMod(common.ecm.Aux2, MontgomeryMultR1, common.ecm.Aux4);
        UncompressLimbsBigInteger(common.ecm.Aux4, &Temp2);
        BigIntGcd(pValue, &Temp2, &Temp4);
        if ((Temp4.nbrLimbs != 1 || Temp4.limbs[0].x > 1) &&
          (Temp4.nbrLimbs != NumberLength ||
            memcmp(pValue->limbs, Temp4.limbs, NumberLength * sizeof(limb))))
        {          // Non-trivial factor found.
          insertBigFactor(pstFactors, &Temp4, TYP_RABIN);
          factorsFound = TRUE;
        }
        i = ctr;
        continue;  // Find more factors.
      }
      if (checkMinusOne(common.ecm.Aux3, nbrLimbs) != 0)
      {            // Square root of 1 found.
        if (!sqrtMinusOneFound)
        {          // Save it to perform GCD later.
          memcpy(common.ecm.Xaux, common.ecm.Aux2, nbrLimbs * sizeof(limb));
          sqrtOneFound = TRUE;
        }
        else
        {          // Try to find non-trivial factor by doing GCD.
          SubtBigNbrMod(common.ecm.Aux3, common.ecm.Xaux, common.ecm.Aux4);
          UncompressLimbsBigInteger(common.ecm.Aux4, &Temp2);
          BigIntGcd(pValue, &Temp2, &Temp4);
          if ((Temp4.nbrLimbs != 1 || Temp4.limbs[0].x > 1) &&
            (Temp4.nbrLimbs != NumberLength ||
              memcmp(pValue->limbs, Temp4.limbs, NumberLength * sizeof(limb))))
          {          // Non-trivial factor found.
            insertBigFactor(pstFactors, &Temp4, TYP_RABIN);
            factorsFound = TRUE;
          }
        }
        i = ctr;
        continue;  // Find more factors.
      }
      memcpy(common.ecm.Aux2, common.ecm.Aux3, nbrLimbs * sizeof(limb));
    }
  }
  return factorsFound;
}

static void factorSmallInt(int toFactor, int* factors, struct sFactors* pstFactors)
{
  int factorsFound = 0;
  int primeFactor;
  int multiplicity;
  struct sFactors* ptrFactor = pstFactors + 1;
  int* ptrFactorLimbs = factors;
  if (toFactor <= 3)
  {     // Only one factor.
    ptrFactor->ptrFactor = ptrFactorLimbs;
    ptrFactor->multiplicity = 1;
    ptrFactor->type = 0;
    ptrFactor->upperBound = 0;
    *ptrFactorLimbs++ = 1;
    *ptrFactorLimbs++ = toFactor;
    pstFactors->multiplicity = 1;
    return;
  }
  // Divide by 2 and 3.
  for (primeFactor = 2; primeFactor <= 3; primeFactor++)
  {
    multiplicity = 0;
    while (toFactor % primeFactor == 0)
    {
      toFactor /= primeFactor;
      multiplicity++;
    }
    if (multiplicity > 0)
    {
      factorsFound++;
      ptrFactor->ptrFactor = ptrFactorLimbs;
      ptrFactor->multiplicity = multiplicity;
      ptrFactor->type = 0;
      ptrFactor->upperBound = 0;
      *ptrFactorLimbs++ = 1;
      *ptrFactorLimbs++ = primeFactor;
      ptrFactor++;
    }
  }
  for (primeFactor = 5;
    (unsigned int)primeFactor * (unsigned int)primeFactor <= (unsigned int)toFactor;
    primeFactor += 2)
  {
    if (primeFactor % 3 == 0)
    {
      continue;
    }
    multiplicity = 0;
    while (toFactor % primeFactor == 0)
    {
      toFactor /= primeFactor;
      multiplicity++;
    }
    if (multiplicity > 0)
    {
      factorsFound++;
      ptrFactor->ptrFactor = ptrFactorLimbs;
      ptrFactor->multiplicity = multiplicity;
      ptrFactor->type = 0;
      ptrFactor->upperBound = 0;
      *ptrFactorLimbs++ = 1;
      *ptrFactorLimbs++ = primeFactor;
      ptrFactor++;
    }
  }
  if (toFactor > 1)
  {
    factorsFound++;
    ptrFactor->ptrFactor = ptrFactorLimbs;
    ptrFactor->multiplicity = 1;
    ptrFactor->type = 0;
    ptrFactor->upperBound = 0;
    *ptrFactorLimbs++ = 1;
    *ptrFactorLimbs++ = toFactor;
    ptrFactor++;
  }
  pstFactors->multiplicity = factorsFound;
}

void factor(BigInteger* toFactor, int* number, int* factors, struct sFactors* pstFactors)
{
  factorExt(toFactor, number, factors, pstFactors, NULL);
}

// pstFactors -> ptrFactor points to end of factors.
// pstFactors -> multiplicity indicates the number of different factors.
void factorExt(BigInteger *toFactor, int *number, int *factors, struct sFactors *pstFactors, char *pcKnownFactors)
{
  struct sFactors *pstCurFactor;
  int factorNbr, expon;
  int remainder, nbrLimbs, ctr;
  int *ptrFactor;
  int dividend;
  char *ptrCharFound;
  int result;
  initializeSmallPrimes(smallPrimes);
  if (toFactor->nbrLimbs == 1)
  {
    factorSmallInt(toFactor->limbs[0].x, factors, pstFactors);
    return;
  }
  oldNbrFactors = 0;
  NextEC = -1;
  EC = 1;
  NumberLength = toFactor->nbrLimbs;
  GetYieldFrequency();
#ifdef __EMSCRIPTEN__
  oldTimeElapsed = 0;
  lModularMult = 0;
  SIQSModMult = 0;
  primeModMult = 0;
  nbrECM = 0;
  nbrSIQS = 0;
  nbrPrimalityTests = 0;
  timeECM = 0;
  timePrimalityTests = 0;
  timeSIQS = 0;
  originalTenthSecond = tenths();
  modmultCallback = showECMStatus;   // Set callback.
#endif
  pstCurFactor = pstFactors + 1;
  if (pcKnownFactors == NULL)
  {   // No factors known.
    memcpy(factors, number, (1 + *number) * sizeof(int));
    pstFactors->multiplicity = 1;
    pstFactors->ptrFactor = factors + 1 + *factors;
    pstFactors->upperBound = 0;
    pstCurFactor->multiplicity = 1;
    pstCurFactor->ptrFactor = factors;
    pstCurFactor->upperBound = 2;
#ifdef __EMSCRIPTEN__
    SaveFactors(pstFactors);
#endif
  }
  else
  {   // Insert factors saved on Web Storage.
    pstFactors->multiplicity = 0;
    pstFactors->ptrFactor = factors;
    while (*pcKnownFactors != 0)
    {
      ptrCharFound = findChar(pcKnownFactors, '^');
      if (ptrCharFound == NULL)
      {
        break;
      }
      *ptrCharFound = 0;
      Dec2Bin(pcKnownFactors, prime.limbs, (int)(ptrCharFound - pcKnownFactors), &prime.nbrLimbs);
      BigInteger2IntArray(pstFactors->ptrFactor, &prime);
      pcKnownFactors = ptrCharFound + 1;
      if (getNextInteger(&pcKnownFactors, &pstCurFactor->multiplicity, '('))
      {     // Error on processing exponent.
        break;
      }
      char *ptrCharFound = findChar(pcKnownFactors, ',');
      if (ptrCharFound != NULL)
      {
        if (getNextInteger(&pcKnownFactors, &pstCurFactor->upperBound, ','))
        {     // Error on processing upper bound.
          break;
        }
        if (getNextInteger(&pcKnownFactors, &pstCurFactor->type, ')'))
        {     // Error on processing upper bound.
          break;
        }
      }
      else
      {
        if (getNextInteger(&pcKnownFactors, &pstCurFactor->upperBound, ')'))
        {     // Error on processing upper bound.
          break;
        }
        pstCurFactor->type = 0;
      }
      pstFactors->multiplicity++;
      pstCurFactor->ptrFactor = pstFactors->ptrFactor;
      pstFactors->ptrFactor += 1 + *pstFactors->ptrFactor;
      pstCurFactor++;
      if (*pcKnownFactors == '*')
      {
        pcKnownFactors++;  // Skip multiplication sign.
      }
      if (*pcKnownFactors == ';')
      {
        pcKnownFactors++;  // Skip separation between known factors and factor entered by user.
        Dec2Bin(pcKnownFactors, prime.limbs, (int)strlen(pcKnownFactors), &prime.nbrLimbs);
        if (foundByLehman)
        {
          insertBigFactor(pstFactors, &prime, TYP_LEHMAN);
        }
        else if (EC >= TYP_SIQS)
        {
          insertBigFactor(pstFactors, &prime, TYP_SIQS);
        }
        else
        {
          insertBigFactor(pstFactors, &prime, 0);
        }
#ifdef __EMSCRIPTEN__
        SaveFactors(pstFactors);
#endif
        break;
      }
      if (*pcKnownFactors == ',')
      {
        NextEC = 0;        // Curve number
        while (*++pcKnownFactors != 0)
        {
          NextEC = NextEC * 10 + (*pcKnownFactors & 0x0F);
        }
        oldNbrFactors = pstFactors->multiplicity;
        break;
      }
    }
  }
  if (tofactor.nbrLimbs > 1)
  {
    PowerPM1Check(pstFactors, toFactor);
  }
  pstCurFactor = pstFactors + 1;
  for (factorNbr = 1; factorNbr <= pstFactors->multiplicity; factorNbr++, pstCurFactor++)
  {
    int upperBoundIndex;
    int upperBound = pstCurFactor->upperBound;
    int delta;
    int restartFactoring = FALSE;
    // If number is prime, do not process it.
    if (upperBound == 0)
    {     // Factor is prime.
      continue;
    }
    // Get upperBoundIndex from upperBound.
    upperBoundIndex = 0;
    if (upperBound < 100000)
    {
      for (delta = 8192; delta > 0; delta >>= 1)
      {
        int tempUpperBoundIndex = upperBoundIndex + delta;
        if (tempUpperBoundIndex >= SMALL_PRIMES_ARRLEN)
        {           // Index too large.
          continue;
        }
        if (upperBound == smallPrimes[tempUpperBoundIndex])
        {
          upperBoundIndex = tempUpperBoundIndex;
          break;
        }
        if (upperBound > smallPrimes[tempUpperBoundIndex])
        {
          upperBoundIndex = tempUpperBoundIndex;
        }
      }
    }
    ptrFactor = pstCurFactor->ptrFactor;
    nbrLimbs = *ptrFactor;
    NumberLength = *pstCurFactor->ptrFactor;
    IntArray2BigInteger(pstCurFactor->ptrFactor, &power);
    NumberLength = power.nbrLimbs;
#ifdef __EMSCRIPTEN__
    char *ptrText = ShowFactoredPart(&prime, pstFactors);
    if (!skipPrimality)
    {
#endif
      expon = PowerCheck(&power, &prime);
      if (expon > 1)
      {
        NumberLength = prime.nbrLimbs;
        BigInteger2IntArray(pstCurFactor->ptrFactor, &prime);
        pstCurFactor->multiplicity *= expon;
        nbrLimbs = *pstCurFactor->ptrFactor;
      }
#ifdef __EMSCRIPTEN__
      strcpy(ptrText, lang ? "<p>División por primos menores que 100000.</p>" :
        "<p>Trial division by primes less than 100000.</p>");
      ShowLowerText();
#endif
      while (upperBound < 100000 && nbrLimbs > 1)
      {        // Number has at least 2 limbs: Trial division by small numbers.
        while (pstCurFactor->upperBound != 0)
        {            // Factor found.
          int expon, index, delta;
          ptrFactor = pstCurFactor->ptrFactor;
          remainder = RemDivBigNbrByInt(ptrFactor + 1, upperBound, nbrLimbs);
          if (remainder != 0)
          {    // Factor not found. Use new divisor.
            break;
          }
          // Small factor found. Find the exponent.
          expon = 1;
          index = 0;
          delta = 1;
          CopyBigInt(&common.trialDiv.cofactor, &prime);
          subtractdivide(&common.trialDiv.cofactor, 0, upperBound);
          intToBigInteger(&common.trialDiv.power[0], upperBound);
          for (;;)
          {      // Test whether the cofactor is multiple of power.
            BigIntDivide(&common.trialDiv.cofactor, &common.trialDiv.power[index], &common.trialDiv.quotient);
            BigIntMultiply(&common.trialDiv.quotient, &common.trialDiv.power[index], &common.trialDiv.temp);
            if (!BigIntEqual(&common.trialDiv.temp, &common.trialDiv.cofactor))
            {    // Not a multiple, so exit loop.
              break;
            }
            CopyBigInt(&common.trialDiv.cofactor, &common.trialDiv.quotient);
            BigIntMultiply(&common.trialDiv.power[index], &common.trialDiv.power[index], &common.trialDiv.power[index+1]);
            expon += delta;
            delta <<= 1;
            index++;
          }
          for (index--; index >= 0; index--)
          {
            delta >>= 1;
            BigIntDivide(&common.trialDiv.cofactor, &common.trialDiv.power[index], &common.trialDiv.quotient);
            BigIntMultiply(&common.trialDiv.quotient, &common.trialDiv.power[index], &common.trialDiv.temp);
            if (BigIntEqual(&common.trialDiv.temp, &common.trialDiv.cofactor))
            {    // It is a multiple.
              CopyBigInt(&common.trialDiv.cofactor, &common.trialDiv.quotient);
              expon += delta;
            }
          }
          insertIntFactor(pstFactors, pstCurFactor, upperBound, expon, &common.trialDiv.cofactor);
          restartFactoring = TRUE;
          break;
        }
        if (restartFactoring)
        {
          factorNbr = 0;
          pstCurFactor = pstFactors;
          break;
        }
        if (nbrLimbs == 1)
        {     // Number completely factored.
          break;
        }
        upperBound = smallPrimes[++upperBoundIndex];
      }
      if (restartFactoring)
      {
        continue;
      }
      if (nbrLimbs == 1)
      {
        dividend = *(ptrFactor + 1);
        while (upperBound < 65535 &&
                (unsigned int)upperBound * (unsigned int)upperBound <= (unsigned int)dividend)
        {              // Trial division by small numbers.
          if (dividend % upperBound == 0)
          {            // Factor found.
            insertIntFactor(pstFactors, pstCurFactor, upperBound, 1, NULL);
            restartFactoring = TRUE;
            break;
          }
          upperBound = smallPrimes[++upperBoundIndex];
        }
        if (restartFactoring)
        {
          factorNbr = 0;
          pstCurFactor = pstFactors;
          continue;
        }
        pstCurFactor->upperBound = 0;   // Number is prime.
        continue;
      }
      // No small factor. Check whether the number is prime or prime power.
#ifdef __EMSCRIPTEN__
#ifdef FACTORIZATION_APP
      StepECM = 0;
      strcpy(ptrText, lang ? "<p>Verificando si el número es potencia perfecta.<p>" :
        "<p>Testing whether the number is perfect power or not.</p>");
      ShowLowerText();
#else
      databack(lang ? "3<p>Verificando si el número es potencia perfecta.</p>" :
        "3<p>Testing whether the number is perfect power or not.</p>");
#endif
#endif
#ifdef __EMSCRIPTEN__
    }
#endif
#ifdef __EMSCRIPTEN__
    SaveFactors(pstFactors);
#endif
#ifdef FACTORIZATION_APP
    if (skipPrimality)
    {
      skipPrimality = FALSE;
      result = 1;
    }
    else
    {
#ifdef __EMSCRIPTEN__
      if (prime.nbrLimbs > 3322 / BITS_PER_GROUP)  // 3322 = 1000*log_2(10) -> 1000 digits
      {
        startSkipTest();
      }
      uint64_t oldModularMult = lModularMult;
#endif
      result = BpswPrimalityTest(&prime, pstFactors);
#ifdef __EMSCRIPTEN__
      nbrPrimalityTests++;
      primeModMult += lModularMult - oldModularMult;
      if (prime.nbrLimbs > 3322 / BITS_PER_GROUP)  // 3322 = 1000*log_2(10) -> 1000 digits
      {
        endSkipTest();
      }
#endif
    }
#else
    result = BpswPrimalityTest(&prime);
#endif
    if (result == 0)
    {   // Number is prime power.
      pstCurFactor->upperBound = 0;   // Indicate that number is prime.
      continue;                       // Check next factor.
    }
    if (result > 1)
    {    // Number is 2-Fermat probable prime. Try to factor it.
      if (factorCarmichael(&prime, pstFactors))
      {                               // Factors found.
        factorNbr--;                  // Test whether factor found is prime.
        pstCurFactor--;
        continue;
      }
    }
#ifdef __EMSCRIPTEN__
    double originalTenths = tenths();
#endif
    ecm(&prime, pstFactors);          // Factor number.
#ifdef __EMSCRIPTEN__
    nbrECM++;
    timeECM += (int)(tenths() - originalTenths);
#endif
    // Check whether GD is not one. In this case we found a proper factor.
    for (ctr = 1; ctr < NumberLength; ctr++)
    {
      if (common.ecm.GD[ctr].x != 0)
      {
        break;
      }
    }
    if (ctr != NumberLength || common.ecm.GD[0].x != 1)
    {
      int numLimbs;
      Temp1.sign = SIGN_POSITIVE;
      numLimbs = NumberLength;
      while (numLimbs > 1)
      {
        if (common.ecm.GD[numLimbs-1].x != 0)
        {
          break;
        }
        numLimbs--;
      }
      memcpy(Temp1.limbs, common.ecm.GD, numLimbs * sizeof(limb));
      Temp1.nbrLimbs = numLimbs;
      if (foundByLehman)
      {
        insertBigFactor(pstFactors, &Temp1, TYP_LEHMAN + EC);
      }
      else
      {
        insertBigFactor(pstFactors, &Temp1, EC);
      }
#ifdef __EMSCRIPTEN__
      SaveFactors(pstFactors);
#endif
      factorNbr = 0;
      pstCurFactor = pstFactors;
    }    // End if
  }      // End for
#ifdef __EMSCRIPTEN__
  SaveFactors(pstFactors);
#endif
}

EXTERNALIZE char *getFactorsAsciiPtr(void)
{
  return common.saveFactors.text;
}

static void intArrayToBigInteger(int *ptrValues, BigInteger *bigint)
{
  bigint->sign = SIGN_POSITIVE;
  bigint->nbrLimbs = *ptrValues;
  memcpy(bigint->limbs, ptrValues + 1, *ptrValues * sizeof(int));
}

// Find Euler's Totient as the product of p^(e-1)*(p-1) where p=prime and e=exponent.
void Totient(BigInteger *result)
{
  static BigInteger TempVar;
  struct sFactors *pstFactor;
  int factorNumber;
  intToBigInteger(result, 1);  // Set result to 1.
  pstFactor = &astFactorsMod[1];
  for (factorNumber = 1; factorNumber <= astFactorsMod[0].multiplicity; factorNumber++)
  {
    intArrayToBigInteger(pstFactor->ptrFactor, &factorValue);
    if (factorValue.nbrLimbs == 1 && factorValue.limbs[0].x == 1)
    {   // If factor is 1 do not do anything.
      continue;
    }
    BigIntPowerIntExp(&factorValue, pstFactor->multiplicity - 1, &TempVar);   // p^(e-1)
    BigIntMultiply(result, &TempVar, result);
    intArrayToBigInteger(pstFactor->ptrFactor, &TempVar);
    addbigint(&TempVar, -1);   // p-1
    BigIntMultiply(result, &TempVar, result);
    pstFactor++;
  }
}

void NumFactors(BigInteger* result)
{
  intToBigInteger(result, astFactorsMod[0].multiplicity);
}

void MinFactor(BigInteger* result)
{
  int factorNumber;
  struct sFactors* pstFactor = &astFactorsMod[1];
  intArrayToBigInteger(pstFactor->ptrFactor, result);
  for (factorNumber = 2; factorNumber <= astFactorsMod[0].multiplicity; factorNumber++)
  {
    intArrayToBigInteger((++pstFactor)->ptrFactor, &factorValue);
    BigIntSubt(&factorValue, result, &factorValue);
    if (factorValue.sign == SIGN_NEGATIVE)
    {
      intArrayToBigInteger(pstFactor->ptrFactor, result);
    }
  }
}

void MaxFactor(BigInteger* result)
{
  int factorNumber;
  struct sFactors* pstFactor = &astFactorsMod[1];
  intArrayToBigInteger(pstFactor->ptrFactor, result);
  for (factorNumber = 2; factorNumber <= astFactorsMod[0].multiplicity; factorNumber++)
  {
    intArrayToBigInteger((++pstFactor)->ptrFactor, &factorValue);
    BigIntSubt(&factorValue, result, &factorValue);
    if (factorValue.sign == SIGN_POSITIVE)
    {
      intArrayToBigInteger(pstFactor->ptrFactor, result);
    }
  }
}

// Find sum of divisors as the product of (p^(e+1)-1)/(p-1) where p=prime and e=exponent.
void SumOfDivisors(BigInteger *result)
{
  struct sFactors *pstFactor;
  int factorNumber;
  intToBigInteger(result, 1);  // Set result to 1.
  pstFactor = &astFactorsMod[1];
  for (factorNumber = 1; factorNumber <= astFactorsMod[0].multiplicity; factorNumber++)
  {
    intArrayToBigInteger(pstFactor->ptrFactor, &factorValue);
    if (factorValue.nbrLimbs == 1 && factorValue.limbs[0].x == 1)
    {   // If factor is 1 do not do anything.
      continue;
    }
    BigIntPowerIntExp(&factorValue, pstFactor->multiplicity + 1, &Temp1);   // p^(e+1)
    addbigint(&Temp1, -1);   // p^(e+1)-1
    BigIntMultiply(result, &Temp1, &Temp2);
    intArrayToBigInteger(pstFactor->ptrFactor, &Temp1);
    addbigint(&Temp1, -1);   // p-1
    BigIntDivide(&Temp2, &Temp1, result);
    pstFactor++;
  }
}

// Find number of divisors as the product of e+1 where p=prime and e=exponent.
void NumberOfDivisors(BigInteger *result)
{
  struct sFactors *pstFactor;
  int factorNumber;
  intToBigInteger(result, 1);  // Set result to 1.
  pstFactor = &astFactorsMod[1];
  for (factorNumber = 1; factorNumber <= astFactorsMod[0].multiplicity; factorNumber++)
  {
    intArrayToBigInteger(pstFactor->ptrFactor, &factorValue);
    if (factorValue.nbrLimbs == 1 && factorValue.limbs[0].x == 1)
    {   // If factor is 1 do not do anything.
      continue;
    }
    multint(result, result, pstFactor->multiplicity + 1);
    pstFactor++;
  }
}

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

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bignbr.h"
#include "expression.h"
#include "highlevel.h"
#include "showtime.h"
#include "batch.h"
#include "tsquares.h"

#define NBR_PRIMES_QUADR_SIEVE  17

const char* square = "<span class=\"bigger\">²</span>";
bool Computing3Squares;
bool hexadecimal;
int iMult3;
int iMult4;
int sieve[MAX_SIEVE];
int primediv[256];
int primeexp[256];
limb Mult2[MAX_LEN];
limb valueP[MAX_LEN];
limb number[MAX_LEN];
BigInteger SquareMult1;
BigInteger SquareMult2;
BigInteger SquareMult3;
BigInteger SquareMult4;
int Mult1Len;
int Mult2Len;
int Mult3Len;
int Mult4Len;
int nbrLimbsP;
int sum;
int nbrLimbs;
int nbrModExp;
BigInteger biMult1;
BigInteger biMult2;
BigInteger biMult3;
BigInteger biMult4;

bool FindTwoSquaresNoNumTheory(void)
{
  // In the following arrays, the bit n is set if n is a square mod p.
  // Bits 31-0
  const unsigned int bitsSqrLow[NBR_PRIMES_QUADR_SIEVE] =
  {
    0x00000003U, // squares mod 3
    0x00000013U, // squares mod 5
    0x00000017U, // squares mod 7
    0x0000023BU, // squares mod 11
    0x0000161BU, // squares mod 13
    0x0001A317U, // squares mod 17
    0x00030AF3U, // squares mod 19
    0x0005335FU, // squares mod 23
    0x13D122F3U, // squares mod 29
    0x121D47B7U, // squares mod 31
    0x5E211E9BU, // squares mod 37
    0x82B50737U, // squares mod 41
    0x83A3EE53U, // squares mod 43
    0x1B2753DFU, // squares mod 47
    0x3303AED3U, // squares mod 53
    0x3E7B92BBU, // squares mod 59
    0x0A59F23BU, // squares mod 61
  };
  // Bits 63-32
  const unsigned int bitsSqrHigh[NBR_PRIMES_QUADR_SIEVE] =
  {
    0x00000000U, // squares mod 3
    0x00000000U, // squares mod 5
    0x00000000U, // squares mod 7
    0x00000000U, // squares mod 11
    0x00000000U, // squares mod 13
    0x00000000U, // squares mod 17
    0x00000000U, // squares mod 19
    0x00000000U, // squares mod 23
    0x00000000U, // squares mod 29
    0x00000000U, // squares mod 31
    0x00000016U, // squares mod 37
    0x000001B3U, // squares mod 41
    0x00000358U, // squares mod 43
    0x00000435U, // squares mod 47
    0x0012DD70U, // squares mod 53
    0x022B6218U, // squares mod 59
    0x1713E694U, // squares mod 61
  };
  const int primes[NBR_PRIMES_QUADR_SIEVE] =
  { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61 };
  int nbrs[NBR_PRIMES_QUADR_SIEVE];
  int diffs[NBR_PRIMES_QUADR_SIEVE];
  int i;
  int nbrIterations;
  static BigInteger nextroot;
  static BigInteger c;
  static BigInteger val;
  if (sum > 0)
  {      // Do not attempt to find sum of two squares using this method
         // if number requires 3 or 4 squares.
    return false;
  }
  int lenBytes = nbrLimbsP * (int)sizeof(limb);
  (void)memcpy(biMult3.limbs, valueP, lenBytes);  // Mult3 <- number to find
  biMult3.sign = SIGN_POSITIVE;         // the decomposition into
  biMult3.nbrLimbs = nbrLimbsP;         // two squares.
  squareRoot(biMult3.limbs, biMult1.limbs, biMult3.nbrLimbs, &biMult1.nbrLimbs);
  biMult1.sign = SIGN_POSITIVE;
  (void)BigIntMultiply(&biMult1, &biMult1, &nextroot);
  BigIntSubt(&biMult3, &nextroot, &c);
  for (i = 0; i < NBR_PRIMES_QUADR_SIEVE; i++)
  {
    int pr = primes[i];
    nbrs[i] = getRemainder(&c, pr);    // nbrs[i] <- c % primes[i]
    // Get difference between the square of biMult1 and the next one.
    diffs[i] = ((getRemainder(&biMult1, pr) * 2) + 1) % pr;
  }
  if (nbrLimbsP > 10)
  {
    nbrIterations = 10000;
  }
  else
  {
    nbrIterations = 1000 * nbrLimbsP;
  }
  if ((biMult1.nbrLimbs == 1) && (nbrIterations > biMult1.limbs[0].x))
  {     // The number to be squared must be always positive.
    nbrIterations = biMult1.limbs[0].x - 1;
  }
  for (int iterNbr = 0; iterNbr < nbrIterations; iterNbr++)
  {
    for (i = 0; i < NBR_PRIMES_QUADR_SIEVE; i++)
    {
      unsigned int shiftBits = (unsigned int)nbrs[i];
      unsigned int bitsSqr;
      unsigned int bitsToShift;
      if (shiftBits < 32U)
      {
        bitsSqr = bitsSqrLow[i];
        bitsToShift = shiftBits;
      }
      else
      {
        bitsSqr = bitsSqrHigh[i];
        bitsToShift = shiftBits - 32U;
      }
      if (((bitsSqr >> bitsToShift) & 0x01U) == 0U)
      { // Not a perfect square
        break;
      }
    }
    if (i == NBR_PRIMES_QUADR_SIEVE)
    { // Test for perfect square
      intToBigInteger(&c, iterNbr);         // c <- iterNbr
      BigIntSubt(&biMult1, &c, &val);
      (void)BigIntMultiply(&val, &val, &c); // c <- val * val
      BigIntSubt(&biMult3, &c, &c);         // c <- nbr - val * val - sqr
      squareRoot(c.limbs, biMult2.limbs, c.nbrLimbs, &biMult2.nbrLimbs);
      biMult2.sign = SIGN_POSITIVE;   // sqrRoot <- sqrt(c)
      (void)BigIntMultiply(&biMult2, &biMult2, &val); // val <- sqroot^2
      if (BigIntEqual(&c, &val))
      {
        intToBigInteger(&c, iterNbr);       // c <- m * iterNbr
        BigIntSubt(&biMult1, &c, &biMult1);
        int nbrBytes = (nbrLimbs + 1) * (int)sizeof(limb);
        (void)memset(Mult1, 0, nbrBytes);
        (void)memset(Mult2, 0, nbrBytes);
        nbrBytes = biMult1.nbrLimbs * (int)sizeof(limb);
        (void)memcpy(Mult1, biMult1.limbs, nbrBytes);
        nbrBytes = biMult2.nbrLimbs * (int)sizeof(limb);
        (void)memcpy(Mult2, biMult2.limbs, nbrBytes);
        return true;          // Decomposition in two squares was found.
      }
    }
    for (i = 0; i < NBR_PRIMES_QUADR_SIEVE; i++)
    {
      int pr = primes[i];
      diffs[i] = (diffs[i] + pr - 2) % pr;
      nbrs[i] = (nbrs[i] + diffs[i]) % pr;
    }
  }
  return false;               // Decomposition in two squares not found.
}

void FillSieveArray(const limb *value)
{
  for (int i = 0; i < (MAX_SIEVE / 2); i++)
  {
    if (sieve[i] >= 0)                        // If prime...
    {
      int Rem;
      int LimbModQ;
      int Q;
      const limb* ptrValue = value + nbrLimbs;
      Rem = 0;
      Q = (2 * i) + 3;                        // Prime
      LimbModQ = UintToInt(LIMB_RANGE % (unsigned int)Q);
      for (int j = nbrLimbs - 1; j >= 0; j--)
      {
        ptrValue--;
        int Divid = ptrValue->x + (Rem * LimbModQ);
        Rem = UintToInt((unsigned int)Divid % (unsigned int)Q);
      }
      sieve[i] = Rem;
    }
  }
}

bool FindSumTwoSquaresNoSmallDivisors(void)
{
  bool sqrtFound;
  int nbrLimbsQ;
  int shRightMult3;
  int lenBytes;
  int base;
  // At this moment p should be prime, otherwise we must try another number.
  valueP[nbrLimbsP].x = 0;
  lenBytes = (nbrLimbsP + 1) * (int)sizeof(valueP[0]);
  (void)memcpy(valueQ, valueP, lenBytes);
  nbrLimbsQ = nbrLimbsP;
  valueQ[0].x--;                     // q = p - 1 (p is odd, so there is no carry).
  lenBytes = (nbrLimbsQ + 1) * (int)sizeof(valueQ[0]);
  (void)memcpy(Mult3, valueQ, lenBytes);
  Mult3Len = nbrLimbsP;
  DivideBigNbrByMaxPowerOf2(&shRightMult3, Mult3, &Mult3Len);
  base = 1;
  if (FindTwoSquaresNoNumTheory())
  {                  // Two squares have been found.
    return true;
  }
  lenBytes = (nbrLimbsP + 1) * (int)sizeof(valueP[0]);
  (void)memcpy(TestNbr, valueP, lenBytes);
  GetMontgomeryParms(nbrLimbsP);
  sqrtFound = false;
  do
  {                // Compute Mult1 = sqrt(-1) (mod p).
    ShowStatus();
    base++;
    modPowBaseInt(base, Mult3, Mult3Len, Mult1); // Mult1 = base^Mult3.
    nbrModExp++;   // Increment number of modular exponentiations.    
    for (int i = 0; i < shRightMult3; i++)
    {              // Loop that squares number.
      modmult(Mult1, Mult1, Mult4);
      if (checkMinusOne(Mult4, nbrLimbsP) != 0)
      {
        sqrtFound = true;
        break;      // Mult1^2 = -1 (mod p), so exit loop.
      }
      lenBytes = nbrLimbsP * (int)sizeof(limb);
      (void)memcpy(Mult1, Mult4, lenBytes);
    }
    // If power (Mult4) is 1, that means that number is at least PRP,
    // so continue loop trying to find square root of -1.
    lenBytes = nbrLimbsP * (int)sizeof(limb);
  } while (memcmp(Mult4, MontgomeryMultR1, lenBytes) == 0);
  if (!sqrtFound)
  {            // Cannot find sqrt(-1) (mod p), try next candidate.
    return false;
  }
  // Convert Mult1 from Montgomery notation to standard number
  // by multiplying by 1 in Montgomery notation.
  lenBytes = nbrLimbsP * (int)sizeof(limb);
  (void)memset(Mult2, 0, lenBytes);
  Mult2[0].x = 1;
  modmult(Mult2, Mult1, Mult1);  // Mult1 = sqrt(-1) mod p.
  // Find the sum of two squares that equals this PRP x^2 + y^2 = p using
  // Gaussian GCD as: x + iy = gcd(s + i, p)
  // where s = sqrt(-1) mod p

  // Result is stored in biMult1 and biMult2.
  // Initialize real part to square root of (-1).
  (void)memcpy(biMult1.limbs, Mult1, lenBytes);
  biMult1.nbrLimbs = nbrLimbsP;
  biMult1.sign = SIGN_POSITIVE;
  intToBigInteger(&biMult2, 1);   // Initialize imaginary part to 1.
  while ((biMult1.nbrLimbs > 1) && (biMult1.limbs[biMult1.nbrLimbs - 1].x == 0))
  {
    biMult1.nbrLimbs--;
  }
  // Initialize real part to prime.
  (void)memcpy(biMult3.limbs, TestNbr, lenBytes);
  biMult3.nbrLimbs = nbrLimbsP;
  biMult3.sign = SIGN_POSITIVE;
  while ((biMult3.nbrLimbs > 1) && (biMult3.limbs[biMult3.nbrLimbs - 1].x == 0))
  {
    biMult3.nbrLimbs--;
  }
  intToBigInteger(&biMult4, 0);   // Initialize imaginary part to 0.
  // Find gcd of (biMult1 + biMult2 * i) and 
  // (biMult3 + biMult4 * i)
  GaussianGCD(&biMult1, &biMult2,          // First number.
    &biMult3, &biMult4,          // Second number.
    &SquareMult1, &SquareMult2,  // GCD.
    &SquareMult3, &SquareMult4); // Temporary variables.
  nbrLimbs = biMult1.nbrLimbs;
  if (nbrLimbs < biMult2.nbrLimbs)
  {
    nbrLimbs = biMult2.nbrLimbs;
  }
  lenBytes = (nbrLimbs + 1) * (int)sizeof(limb);
  (void)memset(Mult1, 0, lenBytes);
  (void)memset(Mult2, 0, lenBytes);
  lenBytes = SquareMult1.nbrLimbs * (int)sizeof(limb);
  (void)memcpy(Mult1, SquareMult1.limbs, lenBytes);
  lenBytes = SquareMult2.nbrLimbs * (int)sizeof(limb);
  (void)memcpy(Mult2, SquareMult2.limbs, lenBytes);
  return true;
}

void FindSumOfTwoSquaresUsingCrt(int nbrDivisors)
{
  int nbrBytes;
  for (int i = 0; i < nbrDivisors; i++)
  {
    int index;
    int divisor = primediv[i];
    for (index = primeexp[i] - 1; index > 0; index -= 2)
    {
      MultBigNbrByInt(Mult1, divisor, Mult1, nbrLimbs + 1);
      MultBigNbrByInt(Mult2, divisor, Mult2, nbrLimbs + 1);
      if ((Mult1[nbrLimbs].x != 0) || (Mult2[nbrLimbs].x != 0))
      {
        nbrLimbs++;
        Mult1[nbrLimbs].x = 0;
        Mult2[nbrLimbs].x = 0;
      }
    }
    if (index == 0)
    {
      // Since the value primediv[i] is very low, it is faster to use
      // trial and error than the general method in order to find
      // the sum of two squares.
      int r;
      int j = 1;
      for (;;)
      {
        double dR;
        r = divisor - (j * j);
        dR = sqrt((double)r);
        r = (int)dR;
        if ((r * r) + (j * j) == divisor)
        {
          break;
        }
        j++;
      }
      // Compute Mult1 <- Previous Mult1 * j + Previous Mult2 * r
      // Compute Mult2 <- Previous Mult1 * r - Previous Mult2 * j
      // Use SquareMult1, SquareMult2, SquareMult3, SquareMult4 as temporary storage.
      MultBigNbrByInt(Mult1, j, SquareMult1.limbs, nbrLimbs + 1);
      MultBigNbrByInt(Mult2, r, SquareMult2.limbs, nbrLimbs + 1);
      MultBigNbrByInt(Mult1, r, SquareMult3.limbs, nbrLimbs + 1);
      MultBigNbrByInt(Mult2, j, SquareMult4.limbs, nbrLimbs + 1);
      AddBigNbr(SquareMult1.limbs, SquareMult2.limbs,
        Mult1, nbrLimbs + 1);
      SubtractBigNbrB(SquareMult3.limbs, SquareMult4.limbs,
        Mult2, nbrLimbs + 1);
      if ((unsigned int)Mult2[nbrLimbs].x >= LIMB_RANGE)
      {   // Since Mult2 is a difference of products, it can be
          // negative. In this case replace it by its absolute value.
        ChSignBigNbr(Mult2, nbrLimbs + 1);
      }
      if ((Mult1[nbrLimbs].x != 0) || (Mult2[nbrLimbs].x != 0))
      {
        nbrLimbs++;
        Mult1[nbrLimbs].x = 0;
        Mult2[nbrLimbs].x = 0;
      }
    }
  }
  nbrBytes = nbrLimbs * (int)sizeof(int);
  (void)memcpy(biMult1.limbs, Mult1, nbrBytes);
  (void)memcpy(biMult2.limbs, Mult2, nbrBytes);
  biMult1.sign = SIGN_POSITIVE;
  biMult2.sign = SIGN_POSITIVE;
  biMult1.nbrLimbs = nbrLimbs;
  biMult2.nbrLimbs = nbrLimbs;
  while ((biMult1.nbrLimbs > 1) && (biMult1.limbs[biMult1.nbrLimbs - 1].x == 0))
  {
    biMult1.nbrLimbs--;
  }
  while ((biMult2.nbrLimbs > 1) && (biMult2.limbs[biMult2.nbrLimbs - 1].x == 0))
  {
    biMult2.nbrLimbs--;
  }
}

bool isSumOfTwoSquares(void)
{
  int nbrDivisors;
  int divisor;
  int shRight;
  bool multipleOf4kPlus3 = false;
  limb carry;

  if (FindTwoSquaresNoNumTheory())
  {   // Sum of two squares was found.
    return true;
  }
  // Now compute modulus over the product of all odd divisors
  // less than maxsieve.
  nbrDivisors = 0;
  divisor = 3;
  for (int i = 0; i < (MAX_SIEVE / 2); i++)
  {
    if ((sieve[i] >= 0) && (((sieve[i] - sum) % divisor) == 0))
    {                          // Divisor found.
      primediv[nbrDivisors] = divisor;   // Store divisor.
      primeexp[nbrDivisors] = 0;         // Store exponent.
      for (;;)
      {
        // Compute the remainder of p and divisor.
        int LimbModDivisor = UintToInt(LIMB_RANGE % (unsigned int)divisor);
        carry.x = 0;
        for (int index = nbrLimbsP - 1; index >= 0; index--)
        {
          carry.x = UintToInt((((unsigned int)carry.x * (unsigned int)LimbModDivisor) +
            (unsigned int)valueP[index].x) % (unsigned int)divisor);
        }
        if (carry.x != 0)
        {     // Number is not a multiple of "divisor".
          break;
        }
        // Divide by divisor.
        DivBigNbrByInt(valueP, divisor, valueP, nbrLimbsP);
        if ((nbrLimbsP > 1) && (valueP[nbrLimbsP - 1].x == 0))
        {
          nbrLimbsP--;
        }
        primeexp[nbrDivisors]++;       // Increment exponent.
      }    // End for
      if (((divisor & 0x03) == 3) && ((primeexp[nbrDivisors] & 1) != 0))
      {                   // Divisor of the form 4k+3 not square.
        multipleOf4kPlus3 = true;
        break;            // Discard this value.
      }
      nbrDivisors++;
      if (FindTwoSquaresNoNumTheory())
      {   // Sum of two squares was found.
        FindSumOfTwoSquaresUsingCrt(nbrDivisors);
        return true;
      }
    }
    divisor += 2;
  }         // End for
  if (multipleOf4kPlus3)
  {         // Number cannot be expressed as sum of two squares.
    return false;
  }
  DivideBigNbrByMaxPowerOf2(&shRight, valueP, &nbrLimbsP);
  if ((valueP[0].x & 0x03) != 0x01)
  {                   // p is not congruent to 1 mod 4
    return false;     // so it cannot be a sum of two squares.
  }
  if (shRight != 0)
  {
    primediv[nbrDivisors] = 2;          // Store divisor.
    primeexp[nbrDivisors] = shRight;    // Store exponent.
    nbrDivisors++;
  }
  if ((valueP[0].x == 1) && (nbrLimbsP == 1))
  {         // Number is the product of only small primes 4k+1 and
            // squares of primes 4k+3.
    Mult1[0].x = 1;
    Mult2[0].x = 0;
    Mult1[1].x = 0;
    Mult2[1].x = 0;
    nbrLimbs = 1;
  }
  else
  {
    if (!FindSumTwoSquaresNoSmallDivisors())
    {       // Number cannot be expressed as sum of two squares.
      return false;
    }
  }
  // At this moment the big prime is equal to Mult1^2 + Mult2^2.
  // Use the other divisors of modulus in order to get Mult1 and Mult2.
  FindSumOfTwoSquaresUsingCrt(nbrDivisors);
  return true;
}

void InitSieveArray(void)
{
  int i;

  for (i = 0; i < (MAX_SIEVE / 2); i++)
  {
    sieve[i] = 0;
  }
  for (i = 3; (i * i) < (MAX_SIEVE / 2); i += 2)
  {
    int j = (i * i) - 3;
    if ((j % 2) == 0)
    {
      j = j / 2;
    }
    else
    {
      j = (j + i) / 2;
    }
    for (; j < (MAX_SIEVE / 2); j += i)
    {
      sieve[j] = -1;             // Indicate number is composite.
    }
  }
}

void ShowStatus(void)
{
#ifdef __EMSCRIPTEN__
  char status[200];
  char* ptrStatus;
  int elapsedTime = (int)(tenths() - originalTenthSecond);
  if ((elapsedTime / 10) == (oldTimeElapsed / 10))
  {
    return;
  }
  oldTimeElapsed = elapsedTime;
  ptrStatus = status;
  copyStr(&ptrStatus, lang ? "4<p>Transcurrió " : "4<p>Time elapsed: ");
  GetDHMS(&ptrStatus, elapsedTime / 10);
  copyStr(&ptrStatus, lang ? "&nbsp;&nbsp;&nbsp;Intentando suma de dos cuadrados de n &minus; " : "&nbsp;&nbsp;&nbsp;Attempting sum of two squares of <var>n</var> &minus; ");
  if (hexadecimal)
  {
    int2hex(&ptrStatus, iMult3);
  }
  else
  {
    int2dec(&ptrStatus, iMult3);
  }
  copyStr(&ptrStatus, square);
  if (!Computing3Squares)
  {
    copyStr(&ptrStatus, " &minus; ");
    if (hexadecimal)
    {
      int2hex(&ptrStatus, iMult4);
    }
    else
    {
      int2dec(&ptrStatus, iMult4);
    }
    copyStr(&ptrStatus, square);
  }
  databack(status);
#endif
}


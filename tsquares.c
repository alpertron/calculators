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
#include "string/strings.h"
#include "bignbr.h"
#include "expression.h"
#include "highlevel.h"
#include "showtime.h"
#include "batch.h"
#include "tsquares.h"
#include "copyStr.h"
#include "commonstruc.h"

#ifndef __ANDROID__
union uCommon common;
#endif

#define NBR_PRIMES_QUADR_SIEVE  17

const char* square = "<span class=\"bigger\">²</span>";

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
  int lenBytes = common.misc.nbrLimbsP * (int)sizeof(limb);
  (void)memcpy(common.misc.biMult3.limbs, common.misc.valueP, lenBytes);  // Mult3 <- number to find
  common.misc.biMult3.sign = SIGN_POSITIVE;         // the decomposition into
  common.misc.biMult3.nbrLimbs = common.misc.nbrLimbsP;         // two squares.
  squareRoot(common.misc.biMult3.limbs, common.misc.biMult1.limbs, common.misc.biMult3.nbrLimbs, &common.misc.biMult1.nbrLimbs);
  common.misc.biMult1.sign = SIGN_POSITIVE;
  (void)BigIntMultiply(&common.misc.biMult1, &common.misc.biMult1, &nextroot);
  BigIntSubt(&common.misc.biMult3, &nextroot, &c);
  for (i = 0; i < NBR_PRIMES_QUADR_SIEVE; i++)
  {
    int pr = primes[i];
    nbrs[i] = getRemainder(&c, pr);    // nbrs[i] <- c % primes[i]
    // Get difference between the square of biMult1 and the next one.
    diffs[i] = ((getRemainder(&common.misc.biMult1, pr) * 2) + 1) % pr;
  }
  if (common.misc.nbrLimbsP > 10)
  {
    nbrIterations = 10000;
  }
  else
  {
    nbrIterations = 1000 * common.misc.nbrLimbsP;
  }
  if ((common.misc.biMult1.nbrLimbs == 1) && (nbrIterations > common.misc.biMult1.limbs[0].x))
  {     // The number to be squared must be always positive.
    nbrIterations = common.misc.biMult1.limbs[0].x - 1;
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
      }
      else
      {
        bitsSqr = bitsSqrHigh[i];
      }
      bitsToShift = shiftBits & 0x1FU;
      if (((bitsSqr >> bitsToShift) & 0x01U) == 0U)
      { // Not a perfect square
        break;
      }
    }
    if (i == NBR_PRIMES_QUADR_SIEVE)
    { // Test for perfect square
      intToBigInteger(&c, iterNbr);         // c <- iterNbr
      BigIntSubt(&common.misc.biMult1, &c, &val);
      (void)BigIntMultiply(&val, &val, &c); // c <- val * val
      BigIntSubt(&common.misc.biMult3, &c, &c);         // c <- nbr - val * val - sqr
      squareRoot(c.limbs, common.misc.biMult2.limbs, c.nbrLimbs, &common.misc.biMult2.nbrLimbs);
      common.misc.biMult2.sign = SIGN_POSITIVE;   // sqrRoot <- sqrt(c)
      (void)BigIntMultiply(&common.misc.biMult2, &common.misc.biMult2, &val); // val <- sqroot^2
      if (BigIntEqual(&c, &val))
      {
        intToBigInteger(&c, iterNbr);       // c <- m * iterNbr
        BigIntSubt(&common.misc.biMult1, &c, &common.misc.biMult1);
        int nbrBytes = (common.misc.nbrLimbs + 1) * (int)sizeof(limb);
        (void)memset(common.misc.Mult1, 0, nbrBytes);
        (void)memset(common.misc.Mult2, 0, nbrBytes);
        nbrBytes = common.misc.biMult1.nbrLimbs * (int)sizeof(limb);
        (void)memcpy(common.misc.Mult1, common.misc.biMult1.limbs, nbrBytes);
        nbrBytes = common.misc.biMult2.nbrLimbs * (int)sizeof(limb);
        (void)memcpy(common.misc.Mult2, common.misc.biMult2.limbs, nbrBytes);
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

void FillSieveArray(const limb *value, int *sieveArr)
{
  for (int i = 0; i < (MAX_SIEVE / 2); i++)
  {
    if (*(sieveArr + i) >= 0)                        // If prime...
    {
      int Rem;
      int LimbModQ;
      int Q;
      const limb* ptrValue = value + common.misc.nbrLimbs;
      Rem = 0;
      Q = (2 * i) + 3;                        // Prime
      LimbModQ = UintToInt(LIMB_RANGE % (unsigned int)Q);
      for (int j = common.misc.nbrLimbs - 1; j >= 0; j--)
      {
        ptrValue--;
        int Divid = ptrValue->x + (Rem * LimbModQ);
        Rem = UintToInt((unsigned int)Divid % (unsigned int)Q);
      }
      *(sieveArr + i) = Rem;
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
  common.misc.valueP[common.misc.nbrLimbsP].x = 0;
  lenBytes = (common.misc.nbrLimbsP + 1) * (int)sizeof(common.misc.valueP[0]);
  (void)memcpy(common.misc.valueQ, common.misc.valueP, lenBytes);
  nbrLimbsQ = common.misc.nbrLimbsP;
  common.misc.valueQ[0].x--;                     // q = p - 1 (p is odd, so there is no carry).
  lenBytes = (nbrLimbsQ + 1) * (int)sizeof(common.misc.valueQ[0]);
  (void)memcpy(common.misc.Mult3, common.misc.valueQ, lenBytes);
  common.misc.Mult3Len = common.misc.nbrLimbsP;
  DivideBigNbrByMaxPowerOf2(&shRightMult3, common.misc.Mult3, &common.misc.Mult3Len);
  base = 1;
  if (FindTwoSquaresNoNumTheory())
  {                  // Two squares have been found.
    return true;
  }
  lenBytes = (common.misc.nbrLimbsP + 1) * (int)sizeof(common.misc.valueP[0]);
  (void)memcpy(TestNbr, common.misc.valueP, lenBytes);
  GetMontgomeryParms(common.misc.nbrLimbsP);
  sqrtFound = false;
  do
  {                // Compute Mult1 = sqrt(-1) (mod p).
    ShowStatus();
    base++;
    modPowBaseInt(base, common.misc.Mult3, common.misc.Mult3Len, common.misc.Mult1); // Mult1 = base^Mult3.
    common.misc.nbrModExp++;   // Increment number of modular exponentiations.    
    for (int i = 0; i < shRightMult3; i++)
    {              // Loop that squares number.
      modmult(common.misc.Mult1, common.misc.Mult1, common.misc.Mult4);
      if (checkMinusOne(common.misc.Mult4, common.misc.nbrLimbsP) != 0)
      {
        sqrtFound = true;
        break;      // Mult1^2 = -1 (mod p), so exit loop.
      }
      lenBytes = common.misc.nbrLimbsP * (int)sizeof(limb);
      (void)memcpy(common.misc.Mult1, common.misc.Mult4, lenBytes);
    }
    // If power (Mult4) is 1, that means that number is at least PRP,
    // so continue loop trying to find square root of -1.
    lenBytes = common.misc.nbrLimbsP * (int)sizeof(limb);
  } while (memcmp(common.misc.Mult4, MontgomeryMultR1, lenBytes) == 0);
  if (!sqrtFound)
  {            // Cannot find sqrt(-1) (mod p), try next candidate.
    return false;
  }
  // Convert Mult1 from Montgomery notation to standard number
  // by multiplying by 1 in Montgomery notation.
  lenBytes = common.misc.nbrLimbsP * (int)sizeof(limb);
  (void)memset(common.misc.Mult2, 0, lenBytes);
  common.misc.Mult2[0].x = 1;
  modmult(common.misc.Mult2, common.misc.Mult1, common.misc.Mult1);  // Mult1 = sqrt(-1) mod p.
  // Find the sum of two squares that equals this PRP x^2 + y^2 = p using
  // Gaussian GCD as: x + iy = gcd(s + i, p)
  // where s = sqrt(-1) mod p

  // Result is stored in biMult1 and biMult2.
  // Initialize real part to square root of (-1).
  (void)memcpy(common.misc.biMult1.limbs, common.misc.Mult1, lenBytes);
  common.misc.biMult1.nbrLimbs = common.misc.nbrLimbsP;
  common.misc.biMult1.sign = SIGN_POSITIVE;
  intToBigInteger(&common.misc.biMult2, 1);   // Initialize imaginary part to 1.
  while ((common.misc.biMult1.nbrLimbs > 1) &&
    (common.misc.biMult1.limbs[common.misc.biMult1.nbrLimbs - 1].x == 0))
  {
    common.misc.biMult1.nbrLimbs--;
  }
  // Initialize real part to prime.
  (void)memcpy(common.misc.biMult3.limbs, TestNbr, lenBytes);
  common.misc.biMult3.nbrLimbs = common.misc.nbrLimbsP;
  common.misc.biMult3.sign = SIGN_POSITIVE;
  while ((common.misc.biMult3.nbrLimbs > 1) &&
    (common.misc.biMult3.limbs[common.misc.biMult3.nbrLimbs - 1].x == 0))
  {
    common.misc.biMult3.nbrLimbs--;
  }
  intToBigInteger(&common.misc.biMult4, 0);   // Initialize imaginary part to 0.
  // Find gcd of (biMult1 + biMult2 * i) and 
  // (biMult3 + biMult4 * i)
  GaussianGCD(&common.misc.biMult1, &common.misc.biMult2,          // First number.
    &common.misc.biMult3, &common.misc.biMult4,          // Second number.
    &common.misc.SquareMult1, &common.misc.SquareMult2,  // GCD.
    &common.misc.SquareMult3, &common.misc.SquareMult4); // Temporary variables.
  common.misc.nbrLimbs = common.misc.biMult1.nbrLimbs;
  if (common.misc.nbrLimbs < common.misc.biMult2.nbrLimbs)
  {
    common.misc.nbrLimbs = common.misc.biMult2.nbrLimbs;
  }
  lenBytes = (common.misc.nbrLimbs + 1) * (int)sizeof(limb);
  (void)memset(common.misc.Mult1, 0, lenBytes);
  (void)memset(common.misc.Mult2, 0, lenBytes);
  lenBytes = common.misc.SquareMult1.nbrLimbs * (int)sizeof(limb);
  (void)memcpy(common.misc.Mult1, common.misc.SquareMult1.limbs, lenBytes);
  lenBytes = common.misc.SquareMult2.nbrLimbs * (int)sizeof(limb);
  (void)memcpy(common.misc.Mult2, common.misc.SquareMult2.limbs, lenBytes);
  return true;
}

void FindSumOfTwoSquaresUsingCrt(int nbrDivisors)
{
  int nbrBytes;
  for (int i = 0; i < nbrDivisors; i++)
  {
    int index;
    int divisor = common.misc.primediv[i];
    for (index = common.misc.primeexp[i] - 1; index > 0; index -= 2)
    {
      MultBigNbrByInt(common.misc.Mult1, divisor, common.misc.Mult1, common.misc.nbrLimbs + 1);
      MultBigNbrByInt(common.misc.Mult2, divisor, common.misc.Mult2, common.misc.nbrLimbs + 1);
      if ((common.misc.Mult1[common.misc.nbrLimbs].x != 0) || (common.misc.Mult2[common.misc.nbrLimbs].x != 0))
      {
        common.misc.nbrLimbs++;
        common.misc.Mult1[common.misc.nbrLimbs].x = 0;
        common.misc.Mult2[common.misc.nbrLimbs].x = 0;
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
      MultBigNbrByInt(common.misc.Mult1, j, common.misc.SquareMult1.limbs, common.misc.nbrLimbs + 1);
      MultBigNbrByInt(common.misc.Mult2, r, common.misc.SquareMult2.limbs, common.misc.nbrLimbs + 1);
      MultBigNbrByInt(common.misc.Mult1, r, common.misc.SquareMult3.limbs, common.misc.nbrLimbs + 1);
      MultBigNbrByInt(common.misc.Mult2, j, common.misc.SquareMult4.limbs, common.misc.nbrLimbs + 1);
      AddBigNbr(common.misc.SquareMult1.limbs, common.misc.SquareMult2.limbs,
        common.misc.Mult1, common.misc.nbrLimbs + 1);
      SubtractBigNbrB(common.misc.SquareMult3.limbs, common.misc.SquareMult4.limbs,
        common.misc.Mult2, common.misc.nbrLimbs + 1);
      if ((unsigned int)common.misc.Mult2[common.misc.nbrLimbs].x >= LIMB_RANGE)
      {   // Since Mult2 is a difference of products, it can be
          // negative. In this case replace it by its absolute value.
        ChSignBigNbr(common.misc.Mult2, common.misc.nbrLimbs + 1);
      }
      if ((common.misc.Mult1[common.misc.nbrLimbs].x != 0) || (common.misc.Mult2[common.misc.nbrLimbs].x != 0))
      {
        common.misc.nbrLimbs++;
        common.misc.Mult1[common.misc.nbrLimbs].x = 0;
        common.misc.Mult2[common.misc.nbrLimbs].x = 0;
      }
    }
  }
  nbrBytes = common.misc.nbrLimbs * (int)sizeof(int);
  (void)memcpy(common.misc.biMult1.limbs, common.misc.Mult1, nbrBytes);
  (void)memcpy(common.misc.biMult2.limbs, common.misc.Mult2, nbrBytes);
  common.misc.biMult1.sign = SIGN_POSITIVE;
  common.misc.biMult2.sign = SIGN_POSITIVE;
  common.misc.biMult1.nbrLimbs = common.misc.nbrLimbs;
  common.misc.biMult2.nbrLimbs = common.misc.nbrLimbs;
  while ((common.misc.biMult1.nbrLimbs > 1) && 
    (common.misc.biMult1.limbs[common.misc.biMult1.nbrLimbs - 1].x == 0))
  {
    common.misc.biMult1.nbrLimbs--;
  }
  while ((common.misc.biMult2.nbrLimbs > 1) &&
    (common.misc.biMult2.limbs[common.misc.biMult2.nbrLimbs - 1].x == 0))
  {
    common.misc.biMult2.nbrLimbs--;
  }
}

bool isSumOfTwoSquares(void)
{
  int nbrDivisors;
  int divisor;
  int shRight;
  bool multipleOf4kPlus3 = false;
  limb carry;
  int mod8 = common.misc.valueP[0].x & 7;
  int mod3squared = common.misc.sieve[((3 * 3) - 3) / 2];
  int mod7squared = common.misc.sieve[((7 * 7) - 3) / 2];
  int mod11squared = common.misc.sieve[((11 * 11) - 3) / 2];
  int mod19squared = common.misc.sieve[((19 * 19) - 3) / 2];

  if ((mod8 == 3) || (mod8 == 6) || (mod8 == 7) ||
    (mod3squared == 3) || (mod3squared == 6) ||
    ((mod7squared != 0) && ((mod7squared % 7) == 0)) ||
    ((mod11squared != 0) && ((mod11squared % 11) == 0)) ||
    ((mod19squared != 0) && ((mod19squared % 19) == 0)))

  {    // Number cannot be a sum of two squares in this case.
    return false;
  }
  common.misc.attempts++;
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
    if ((common.misc.sieve[i] >= 0) && ((common.misc.sieve[i] % divisor) == 0))
    {                          // Divisor found.
      common.misc.primediv[nbrDivisors] = divisor;   // Store divisor.
      common.misc.primeexp[nbrDivisors] = 0;         // Store exponent.
      for (;;)
      {
        // Compute the remainder of p and divisor.
        int LimbModDivisor = UintToInt(LIMB_RANGE % (unsigned int)divisor);
        carry.x = 0;
        for (int index = common.misc.nbrLimbsP - 1; index >= 0; index--)
        {
          carry.x = UintToInt((((unsigned int)carry.x * (unsigned int)LimbModDivisor) +
            (unsigned int)common.misc.valueP[index].x) % (unsigned int)divisor);
        }
        if (carry.x != 0)
        {     // Number is not a multiple of "divisor".
          break;
        }
        // Divide by divisor.
        DivBigNbrByInt(common.misc.valueP, divisor, common.misc.valueP, common.misc.nbrLimbsP);
        if ((common.misc.nbrLimbsP > 1) && (common.misc.valueP[common.misc.nbrLimbsP - 1].x == 0))
        {
          common.misc.nbrLimbsP--;
        }
        common.misc.primeexp[nbrDivisors]++;       // Increment exponent.
      }    // End for
      if (((divisor & 0x03) == 3) && ((common.misc.primeexp[nbrDivisors] & 1) != 0))
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
  DivideBigNbrByMaxPowerOf2(&shRight, common.misc.valueP, &common.misc.nbrLimbsP);
  if ((common.misc.valueP[0].x & 0x03) != 0x01)
  {                   // p is not congruent to 1 mod 4
    return false;     // so it cannot be a sum of two squares.
  }
  if (shRight != 0)
  {
    common.misc.primediv[nbrDivisors] = 2;          // Store divisor.
    common.misc.primeexp[nbrDivisors] = shRight;    // Store exponent.
    nbrDivisors++;
  }
  if ((common.misc.valueP[0].x == 1) && (common.misc.nbrLimbsP == 1))
  {         // Number is the product of only small primes 4k+1 and
            // squares of primes 4k+3.
    common.misc.Mult1[0].x = 1;
    common.misc.Mult2[0].x = 0;
    common.misc.Mult1[1].x = 0;
    common.misc.Mult2[1].x = 0;
    common.misc.nbrLimbs = 1;
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
    common.misc.sieve[i] = 0;
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
      common.misc.sieve[j] = -1;             // Indicate number is composite.
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
  copyStr(&ptrStatus, "4<p>");
  // Time elapsed:
  formatString(&ptrStatus, LITERAL_SHOW_STATUS1, elapsedTime / 10);
  // &nbsp;&nbsp;&nbsp;Attempting sum of two squares. Attempt #$1d
  formatString(&ptrStatus, LITERAL_SHOW_STATUS2, common.misc.attempts);
  *ptrStatus = 0;
  databack(status);
#endif
}


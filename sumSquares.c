//
// This file is part of Alpertron Calculators.
//
// Copyright 2023 Dario Alejandro Alpern
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
#include "factor.h"
#include "showtime.h"
#include "output.h"
#include "commonstruc.h"

extern BigInteger tofactor;
static BigInteger Quad1;
static BigInteger Quad2;
static BigInteger Quad3;
static BigInteger Quad4;
static BigInteger valueP;
static BigInteger bigExpon;
static BigInteger bigBase;
static BigInteger Mult1;
static BigInteger Mult2;
static BigInteger Mult3;
static BigInteger Mult4;
static BigInteger Tmp;
static BigInteger Tmp1;
static BigInteger Tmp2;
static BigInteger Tmp3;
static BigInteger Tmp4;
static BigInteger M1;
static BigInteger M2;
static BigInteger M3;
static BigInteger M4;
static BigInteger M5;
static BigInteger M6;
static BigInteger M7;
static BigInteger M8;
void heapsort(int count, int** pptrBigInts);
#ifdef __EMSCRIPTEN__
extern int64_t lModularMult;
#endif

static void modPowShowStatus(const limb* base, const limb* exp, int nbrGroupsExp, limb* power)
{
  int lenBytes = (NumberLength + 1) * (int)sizeof(*power);
  (void)memcpy(power, MontgomeryMultR1, lenBytes);  // power <- 1
  for (int index = nbrGroupsExp - 1; index >= 0; index--)
  {
    int groupExp = (exp + index)->x;
#ifdef __EMSCRIPTEN__
    percentageBPSW = (nbrGroupsExp - index) * 100 / nbrGroupsExp;
#endif
    for (unsigned int mask = HALF_INT_RANGE_U; mask > 0U; mask >>= 1)
    {
      modmult(power, power, power);
      if (((unsigned int)groupExp & mask) != 0U)
      {
        modmult(power, base, power);
      }
    }
  }
}

// A number is a sum of 3 squares if it has not the form 4^n*(8k+7) and
// there is a prime factor of the form 4k+3 with odd multiplicity.

static bool isSumOfThreeSquares(const struct sFactors* pstFactors, BigInteger* pCandidate)
{
  const struct sFactors* pstFactor = pstFactors + 1; // Point to first factor in array of factors.
  int shRight;
  bool evenExponentOf2 = true;
  bool allExponentsEven = true;  // Assume all exponents of primes are even.
  for (int indexPrimes = pstFactors->multiplicity - 1; indexPrimes >= 0; indexPrimes--)
  {
    if ((pstFactor->multiplicity % 2) == 0)
    {                            // Do not consider the case of even multiplicity.
      pstFactor++;
      continue;
    }
    // Prime factor multiplicity is odd.
    if ((*pstFactor->ptrFactor == 1) && (*(pstFactor->ptrFactor + 1) == 2))
    {       // Number is power of 2 but not of 4.
      evenExponentOf2 = false;   // We will determine later if we will use 3 or 4 squares.
    }
    else
    {
      allExponentsEven = false;  // Do not consider the factor 2 on exponent odd.
      if ((*(pstFactor->ptrFactor + 1) % 4) == 1)
      {     // Prime has the form 4k+1, so the candidate can be expressed as sum
            // of two squares.
        return false;
      }
    }
    pstFactor++;
  }
  if (allExponentsEven)
  {         // All exponents are even, so the candidate can by expressed as a
            // sum of two squares.
    return false;
  }
  if (!evenExponentOf2)
  {                        // Number is multiple of power of 2 but not power of 4.
    return true;           // Sum of three squares.
  }
  CopyBigInt(pCandidate, &tofactor);                 // Divide by power of 4.
  if (BigIntIsZero(pCandidate))
  {
    return false;    // Do not express zero as sum of three squares.
  }
  DivideBigNbrByMaxPowerOf4(&shRight, pCandidate->limbs, &pCandidate->nbrLimbs);
  // If candidate divided by power of 4 does not equal 7 (mod 8), then
  // it is a sum of three squares, all of them different from zero.
  return (pCandidate->limbs[0].x % 8) != 7;
}

// Check whether the number is prime.
// pTmp2->limbs: number to check (p).
// pTmp3->limbs: exponent (q).
// pTmp4->limbs: power.
// pTmp1->limbs: temporary storage.
// If we can find the square root of -1, the routine returns true.
static bool TryToFindSqRootMinus1(BigInteger* pTmp1, BigInteger* pTmp2,
  BigInteger* pTmp3, BigInteger* pTmp4)
{
  int shRightPower;
  int powerLen;
  int base;
  int lenBytes;
  int nbrLimbs = pTmp2->nbrLimbs;
  pTmp2->limbs[nbrLimbs].x = 0;
  lenBytes = (nbrLimbs + 1) * (int)sizeof(limb);
  (void)memcpy(pTmp3->limbs, pTmp2->limbs, lenBytes);
  pTmp3->limbs[0].x--;      // q = p - 1 (p is odd, so there is no carry).
  powerLen = nbrLimbs;
  DivideBigNbrByMaxPowerOf2(&shRightPower, pTmp3->limbs, &powerLen);
  base = 1;
  (void)memcpy(TestNbr, pTmp2->limbs, lenBytes);
  GetMontgomeryParms(nbrLimbs);
  do
  {                 // Compute Mult1 = sqrt(-1) (mod p).
    base++;
    // Mult1 = base^pTmp3.
    modPowBaseInt(base, pTmp3->limbs, powerLen, pTmp1->limbs);
#ifdef __EMSCRIPTEN__
    lModularMult++;   // Increment number of modular exponentiations.    
#endif
    for (int i = 0; i < shRightPower; i++)
    {              // Loop that squares number.
      modmult(pTmp1->limbs, pTmp1->limbs, pTmp4->limbs);
      if (checkMinusOne(pTmp4->limbs, nbrLimbs) != 0)
      {
        return true;  // Mult1^2 = -1 (mod p), so exit loop.
      }
      lenBytes = nbrLimbs * (int)sizeof(limb);
      (void)memcpy(pTmp1->limbs, pTmp4->limbs, lenBytes);
    }
    // If power (Mult4) is 1, that means that number is at least PRP,
    // so continue loop trying to find square root of -1.
    lenBytes = nbrLimbs * (int)sizeof(limb);
  } while (memcmp(pTmp4->limbs, MontgomeryMultR1, lenBytes) == 0);
  return false;
}

static void GenerateSumOfTwoSquaresOfPQ(const int* ptrArrFactors,
  BigInteger* pTmp1, BigInteger* pTmp2)
{
  int i;
  int prime = *ptrArrFactors;
  int expon = *(ptrArrFactors + 1);
  for (i = expon - 1; i > 0; i -= 2)
  {
    multint(&Quad1, &Quad1, prime);
    multint(&Quad2, &Quad2, prime);
  }
  if (i == 0)
  {
    // Since the prime is very low, it is faster to use
    // trial and error than the general method in order to find
    // the sum of two squares.
    int j = 1;
    int r;
    for (;;)
    {
      double dR;
      r = prime - (j * j);
      dR = sqrt((double)r);
      r = (int)dR;
      if ((r * r) + (j * j) == prime)
      {
        break;
      }
      j++;
    }
    // Compute Mult1 <- Previous Mult1 * j + Previous Mult2 * r
    // Compute Mult2 <- Previous Mult1 * r - Previous Mult2 * j
    // Use SquareMult1, SquareMult2, SquareMult3, SquareMult4 as temporary storage.
    addmult(pTmp2, &Quad1, j, &Quad2, r);
    addmult(pTmp1, &Quad1, r, &Quad2, -j);
    CopyBigInt(&Quad1, pTmp2);
    CopyBigInt(&Quad2, pTmp1);
  }
}

static bool DivideCandidateBySmallPrimes(BigInteger* pCandidate, int** pptrArrFactors)
{
  int prime = 2;
  int primeIndex = 0;
  int* ptrArrFactors = *pptrArrFactors;
  do
  {
    int expon = 0;
    // Divide Tmp2 by small primes.
    while (getRemainder(pCandidate, prime) == 0)
    {
      subtractdivide(pCandidate, 0, prime);
      expon++;
    }
    if ((((expon & 0x01) == 0x01) && ((prime % 4) == 3)) ||
      ((pCandidate->limbs[0].x % 4) == 3))
    {  // Number cannot be expressed as a sum of three squares.
      return false;
    }
    if (expon > 0)
    {
      *ptrArrFactors = prime;
      ptrArrFactors++;
      *ptrArrFactors = expon;
      ptrArrFactors++;
    }
    primeIndex++;
    prime = smallPrimes[primeIndex];
  } while (prime < 32768);
  *pptrArrFactors = ptrArrFactors;
  return true;
}

// Divide the number by the maximum power of 2, so it is odd.
// Subtract 1, 4, 9, etc. Let the result be n.
// Divide n by small primes. Let n=a1^e1*a2^e2*a3^e3*...*x where a,b,c,d are
// small primes.
// If e_i is odd and a_i=4k+3, then use next value of n.
// If x is not prime, use next value of n.
static void ComputeThreeSquares(BigInteger* pTmp,
  BigInteger* pTmp1, BigInteger* pTmp2, BigInteger* pTmp3, BigInteger* pTmp4,
  BigInteger* pM1, BigInteger* pM2)
{
  int nbrLimbs;
  int arrFactors[400];
  int diff = 1;
  int shRight;
  int* ptrArrFactors;
  CopyBigInt(pTmp, &tofactor);
  DivideBigNbrByMaxPowerOf4(&shRight, pTmp->limbs, &pTmp->nbrLimbs);
  diff = 0;
  for (;;)
  {
    const int* ptrArrFactorsBak;
    diff++;
    ptrArrFactors = arrFactors;
    intToBigInteger(pTmp1, diff * diff);
    BigIntSubt(pTmp, pTmp1, pTmp2);
    if (!DivideCandidateBySmallPrimes(pTmp2, &ptrArrFactors))
    {      // Number cannot be expressed as a sum of three squares.
      continue;
    }
    ptrArrFactorsBak = ptrArrFactors;
    if ((pTmp2->nbrLimbs == 1) && (pTmp2->limbs[0].x == 1))
    {      // Cofactor equals 1.
      intToBigInteger(&Quad1, 1);
      intToBigInteger(&Quad2, 0);
    }
    else if ((pTmp2->nbrLimbs == 1) && (pTmp2->limbs[0].x == 2))
    {      // Cofactor equals 2.
      intToBigInteger(&Quad1, 1);
      intToBigInteger(&Quad2, 1);
    }
    else
    {
      int lenBytes;
      nbrLimbs = pTmp2->nbrLimbs;
      if (!TryToFindSqRootMinus1(pTmp1, pTmp2, pTmp3, pTmp4))
      {            // Cannot find sqrt(-1) (mod p), go to next candidate.
        continue;
      }
      // Convert pTmp1->limbs from Montgomery notation to standard number
      // by multiplying by 1 in Montgomery notation.
      lenBytes = nbrLimbs * (int)sizeof(limb);
      (void)memset(pTmp4->limbs, 0, lenBytes);
      pTmp4->limbs[0].x = 1;
      // pTmp1->limbs = sqrt(-1) mod p.
      modmult(pTmp4->limbs, pTmp1->limbs, pTmp1->limbs);
      pTmp1->nbrLimbs = nbrLimbs;
      while ((pTmp1->nbrLimbs > 0) && (pTmp1->limbs[pTmp1->nbrLimbs - 1].x == 0))
      {
        pTmp1->nbrLimbs--;
      }
      pTmp1->sign = SIGN_POSITIVE;
      // Find the sum of two squares that equals this PRP x^2 + y^2 = p using
      // Gaussian GCD as: x + iy = gcd(s + i, p)
      // where s = sqrt(-1) mod p

      // Result is stored in biMult1 and biMult2.
      // Initialize real part to square root of (-1).
      intToBigInteger(pTmp2, 1);   // Initialize imaginary part to 1.
      // Initialize real part to prime.
      lenBytes = NumberLength * (int)sizeof(limb);
      (void)memcpy(pTmp3->limbs, TestNbr, lenBytes);
      pTmp3->sign = SIGN_POSITIVE;
      pTmp3->nbrLimbs = NumberLength;
      intToBigInteger(pTmp4, 0);   // Initialize imaginary part to 0.
      // Find gcd of (biMult1 + biMult2 * i) and
      // (biMult3 + biMult4 * i)
      GaussianGCD(pTmp1, pTmp2, pTmp3, pTmp4, &Quad1, &Quad2, pM1, pM2);
    }
    for (ptrArrFactors = arrFactors; ptrArrFactors < ptrArrFactorsBak; ptrArrFactors += 2)
    {
      GenerateSumOfTwoSquaresOfPQ(ptrArrFactors, pTmp1, pTmp2);
    }            /* end for */
    Quad1.sign = SIGN_POSITIVE;
    Quad2.sign = SIGN_POSITIVE;
    intToBigInteger(&Quad3, diff);
    intToBigInteger(&Quad4, 0);
    // Perform shift left.
    for (int count = 0; count < shRight; count++)
    {
      BigIntAdd(&Quad1, &Quad1, &Quad1);
      BigIntAdd(&Quad2, &Quad2, &Quad2);
      BigIntAdd(&Quad3, &Quad3, &Quad3);
    }
    BigIntSubt(&Quad1, &Quad2, pTmp);
    if (pTmp->sign == SIGN_NEGATIVE)
    {   // Quad1 < Quad2, so exchange them.
      CopyBigInt(pTmp, &Quad1);
      CopyBigInt(&Quad1, &Quad2);
      CopyBigInt(&Quad2, pTmp);
    }
    return;
  }
}

// Compute prime valueP as Mult1^2 + Mult2^2.
static void ComputeSumOfTwoSquaresForPrime(void)
{
  static limb minusOneMont[MAX_LEN];
  int lenBytes = NumberLength * (int)sizeof(limb);
  (void)memset(minusOneMont, 0, lenBytes);
  SubtBigNbrModN(minusOneMont, MontgomeryMultR1, minusOneMont, TestNbr, NumberLength);
  CopyBigInt(&bigExpon, &valueP);
  subtractdivide(&bigExpon, 1, 4);     // q = (prime-1)/4
  bigBase.limbs[0].x = 1;
  do
  {    // Loop that finds mult1 = sqrt(-1) mod prime in Montgomery notation.
    bigBase.limbs[0].x++;
    modPowShowStatus(bigBase.limbs, bigExpon.limbs, bigExpon.nbrLimbs, Mult1.limbs);
  } while (!memcmp(Mult1.limbs, MontgomeryMultR1, lenBytes) ||
    !memcmp(Mult1.limbs, minusOneMont, lenBytes));
  Mult1.sign = SIGN_POSITIVE;
  lenBytes = valueP.nbrLimbs * (int)sizeof(limb);
  (void)memset(Mult2.limbs, 0, lenBytes);
  Mult2.limbs[0].x = 1;
  Mult2.nbrLimbs = 1;
  Mult2.sign = SIGN_POSITIVE;
  // Convert Mult1 to standard notation by multiplying by 1 in
  // Montgomery notation.
  modmult(Mult1.limbs, Mult2.limbs, Mult3.limbs);
  (void)memcpy(Mult1.limbs, Mult3.limbs, lenBytes);
  for (Mult1.nbrLimbs = valueP.nbrLimbs; Mult1.nbrLimbs > 1; Mult1.nbrLimbs--)
  {  // Adjust number of limbs so the most significant limb is not zero.
    if (Mult1.limbs[Mult1.nbrLimbs - 1].x != 0)
    {
      break;
    }
  }
  // Initialize real part to square root of (-1).
  intToBigInteger(&Mult2, 1);   // Initialize imaginary part to 1.
  // Initialize real part to prime.
  (void)memcpy(Mult3.limbs, TestNbr, lenBytes);
  Mult3.nbrLimbs = valueP.nbrLimbs;
  Mult3.sign = SIGN_POSITIVE;
  while ((Mult3.nbrLimbs > 1) && (Mult3.limbs[Mult3.nbrLimbs - 1].x == 0))
  {
    Mult3.nbrLimbs--;
  }
  intToBigInteger(&Mult4, 0);   // Initialize imaginary part to 0.
  // Find gcd of (Mult1 + Mult2 * i) and (Mult3 + Mult4 * i)
  GaussianGCD(&Mult1, &Mult2, &Mult3, &Mult4, &Tmp1, &Tmp2, &Tmp3, &Tmp4);
  CopyBigInt(&Mult1, &Tmp1);
  CopyBigInt(&Mult2, &Tmp2);
}

// Compute prime p as Mult1^2 + Mult2^2 + Mult3^2 + Mult4^2.
static void ComputeSumOfFourSquaresForPrime(void)
{
  int lenBytes;
  int mult1 = 0;
  // Compute Mult1 and Mult2 so Mult1^2 + Mult2^2 = -1 (mod p)
  intToBigInteger(&Tmp, -1);
  intToBigInteger(&Mult2, -1);
  while (BigIntJacobiSymbol(&Tmp, &valueP) <= 0)
  {     // Not a quadratic residue. Compute next value of -1 - Mult1^2 in variable Tmp.
    BigIntAdd(&Tmp, &Mult2, &Tmp);
    Mult2.limbs[0].x += 2;
    mult1++;
  }
  // After the loop finishes, Tmp = (-1 - Mult1^2) is a quadratic residue mod p.
  // Convert base to Montgomery notation.
  BigIntAdd(&Tmp, &valueP, &Tmp);
  Tmp.limbs[NumberLength].x = 0;
  modmult(Tmp.limbs, MontgomeryMultR2, Tmp.limbs);

  intToBigInteger(&Mult1, mult1);
  CopyBigInt(&bigExpon, &valueP);
  subtractdivide(&bigExpon, -1, 4);  // q <- (p+1)/4.
  // Find Mult2 <- square root of Tmp = Tmp^q (mod p) in Montgomery notation.
  modPowShowStatus(Tmp.limbs, bigExpon.limbs, valueP.nbrLimbs, Mult2.limbs);
  // Convert Mult2 from Montgomery notation to standard notation.
  lenBytes = valueP.nbrLimbs * (int)sizeof(limb);
  (void)memset(Tmp.limbs, 0, lenBytes);
  Tmp.limbs[0].x = 1;
  intToBigInteger(&Mult3, 1);
  intToBigInteger(&Mult4, 0);
  // Convert Mult2 to standard notation by multiplying by 1 in
  // Montgomery notation.
  modmult(Mult2.limbs, Tmp.limbs, Mult2.limbs);
  for (Mult2.nbrLimbs = valueP.nbrLimbs; Mult2.nbrLimbs > 1; Mult2.nbrLimbs--)
  {  // Adjust number of limbs so the most significant limb is not zero.
    if (Mult2.limbs[Mult2.nbrLimbs - 1].x != 0)
    {
      break;
    }
  }
  Mult2.sign = SIGN_POSITIVE;
  MultiplyQuaternionBy2(&Mult1, &Mult2, &Mult3, &Mult4);
  CopyBigInt(&M1, &valueP);
  BigIntMultiplyBy2(&M1);
  intToBigInteger(&M2, 0);
  intToBigInteger(&M3, 0);
  intToBigInteger(&M4, 0);
  QuaternionGCD(&Mult1, &Mult2, &Mult3, &Mult4, &M1, &M2, &M3, &M4, &M5, &M6, &M7, &M8,
    &Tmp, &Tmp1, &Tmp2, &Tmp3);
  DivideQuaternionBy2(&M5, &M6, &M7, &M8);
  CopyBigInt(&Mult1, &M5);
  CopyBigInt(&Mult2, &M6);
  CopyBigInt(&Mult3, &M7);
  CopyBigInt(&Mult4, &M8);
}

// If p = Mult1^2 + Mult2^2 + Mult3^2 + Mult4^2 and 
// q = Quad1^2 + Quad2^2 + Quad3^2 + Quad4^2,
// compute new values of Quad1, Quad2, Quad3 and Quad4 such that:
// pq = Quad1^2 + Quad2^2 + Quad3^2 + Quad4^2.
static void GenerateSumOfFourSquaresOfPQ(void)
{
  // Compute Tmp1 as Mult1*Quad1 + Mult2*Quad2 + Mult3*Quad3 + Mult4*Quad4
  (void)BigIntMultiply(&Mult1, &Quad1, &Tmp);
  (void)BigIntMultiply(&Mult2, &Quad2, &Tmp4);
  BigIntAdd(&Tmp, &Tmp4, &Tmp);
  (void)BigIntMultiply(&Mult3, &Quad3, &Tmp4);
  BigIntAdd(&Tmp, &Tmp4, &Tmp);
  (void)BigIntMultiply(&Mult4, &Quad4, &Tmp4);
  BigIntAdd(&Tmp, &Tmp4, &Tmp1);

  // Compute Tmp2 as Mult1*Quad2 - Mult2*Quad1 + Mult3*Quad4 - Mult4*Quad3
  (void)BigIntMultiply(&Mult1, &Quad2, &Tmp);
  (void)BigIntMultiply(&Mult2, &Quad1, &Tmp4);
  BigIntSubt(&Tmp, &Tmp4, &Tmp);
  (void)BigIntMultiply(&Mult3, &Quad4, &Tmp4);
  BigIntAdd(&Tmp, &Tmp4, &Tmp);
  (void)BigIntMultiply(&Mult4, &Quad3, &Tmp4);
  BigIntSubt(&Tmp, &Tmp4, &Tmp2);

  // Compute Tmp3 as Mult1*Quad3 - Mult3*Quad1 - Mult2*Quad4 + Mult4*Quad2
  (void)BigIntMultiply(&Mult1, &Quad3, &Tmp);
  (void)BigIntMultiply(&Mult3, &Quad1, &Tmp4);
  BigIntSubt(&Tmp, &Tmp4, &Tmp);
  (void)BigIntMultiply(&Mult2, &Quad4, &Tmp4);
  BigIntSubt(&Tmp, &Tmp4, &Tmp);
  (void)BigIntMultiply(&Mult4, &Quad2, &Tmp4);
  BigIntAdd(&Tmp, &Tmp4, &Tmp3);

  // Compute Quad4 as Mult1*Quad4 - Mult4*Quad1 + Mult2*Quad3 - Mult3*Quad2
  (void)BigIntMultiply(&Mult1, &Quad4, &Tmp);
  (void)BigIntMultiply(&Mult4, &Quad1, &Tmp4);
  BigIntSubt(&Tmp, &Tmp4, &Tmp);
  (void)BigIntMultiply(&Mult2, &Quad3, &Tmp4);
  BigIntAdd(&Tmp, &Tmp4, &Tmp);
  (void)BigIntMultiply(&Mult3, &Quad2, &Tmp4);
  BigIntSubt(&Tmp, &Tmp4, &Quad4);

  CopyBigInt(&Quad3, &Tmp3);
  CopyBigInt(&Quad2, &Tmp2);
  CopyBigInt(&Quad1, &Tmp1);
}

// Sort squares so Quad1 is the greatest value
// and Quad4 the lowest one.
static void sortSquares(void)
{
  BigIntSubt(&Quad1, &Quad2, &Tmp);
  if (Tmp.sign == SIGN_NEGATIVE)
  {   // Quad1 < Quad2, so exchange them.
    CopyBigInt(&Tmp, &Quad1);
    CopyBigInt(&Quad1, &Quad2);
    CopyBigInt(&Quad2, &Tmp);
  }
  BigIntSubt(&Quad1, &Quad3, &Tmp);
  if (Tmp.sign == SIGN_NEGATIVE)
  {   // Quad1 < Quad3, so exchange them.
    CopyBigInt(&Tmp, &Quad1);
    CopyBigInt(&Quad1, &Quad3);
    CopyBigInt(&Quad3, &Tmp);
  }
  BigIntSubt(&Quad1, &Quad4, &Tmp);
  if (Tmp.sign == SIGN_NEGATIVE)
  {   // Quad1 < Quad4, so exchange them.
    CopyBigInt(&Tmp, &Quad1);
    CopyBigInt(&Quad1, &Quad4);
    CopyBigInt(&Quad4, &Tmp);
  }
  BigIntSubt(&Quad2, &Quad3, &Tmp);
  if (Tmp.sign == SIGN_NEGATIVE)
  {   // Quad2 < Quad3, so exchange them.
    CopyBigInt(&Tmp, &Quad2);
    CopyBigInt(&Quad2, &Quad3);
    CopyBigInt(&Quad3, &Tmp);
  }
  BigIntSubt(&Quad2, &Quad4, &Tmp);
  if (Tmp.sign == SIGN_NEGATIVE)
  {   // Quad2 < Quad4, so exchange them.
    CopyBigInt(&Tmp, &Quad2);
    CopyBigInt(&Quad2, &Quad4);
    CopyBigInt(&Quad4, &Tmp);
  }
  BigIntSubt(&Quad3, &Quad4, &Tmp);
  if (Tmp.sign == SIGN_NEGATIVE)
  {   // Quad3 < Quad4, so exchange them.
    CopyBigInt(&Tmp, &Quad3);
    CopyBigInt(&Quad3, &Quad4);
    CopyBigInt(&Quad4, &Tmp);
  }
}

void ComputeFourSquares(const struct sFactors* pstFactors)
{
  int indexPrimes;
  int index4k1 = 0;
  const struct sFactors* pstFactor;

  common.k.sumSquares.manyDecompositions = false;
  common.k.sumSquares.twoOddMultiplicity = false;
  intToBigInteger(&Quad1, 1);      // 1 = 1^2 + 0^2 + 0^2 + 0^2
  intToBigInteger(&Quad2, 0);
  intToBigInteger(&Quad3, 0);
  intToBigInteger(&Quad4, 0);
  if ((tofactor.nbrLimbs < 25) && isSumOfThreeSquares(pstFactors, &Tmp))
  {                                // Decompose in sum of 3 squares if less than 200 digits.
    ComputeThreeSquares(&Tmp, &Tmp1, &Tmp2, &Tmp3, &Tmp4, &M1, &M2);
    return;
  }
  pstFactor = pstFactors + 1;      // Point to first factor in array of factors.
  if ((pstFactors->multiplicity == 1) && (*pstFactor->ptrFactor == 1))
  {
    if (*(pstFactor->ptrFactor + 1) == 1)
    {                              // Number to factor is 1.
      return;
    }
    if (*(pstFactor->ptrFactor + 1) == 0)
    {                              // Number to factor is 0.
      intToBigInteger(&Quad1, 0);  // 0 = 0^2 + 0^2 + 0^2 + 0^2
      return;
    }
  }
  intToBigInteger(&common.k.sumSquares.productOtherDivisors, 1);
  for (indexPrimes = 0; indexPrimes < pstFactors->multiplicity; indexPrimes++)
  {
    NumberLength = *pstFactor->ptrFactor;
    IntArray2BigInteger(pstFactor->ptrFactor, &valueP);
    if ((valueP.nbrLimbs == 1) && (valueP.limbs[0].x == 2))
    {
      intToBigInteger(&Mult1, 1);  // 2 = 1^2 + 1^2 + 0^2 + 0^2
      intToBigInteger(&Mult2, 1);
      intToBigInteger(&Mult3, 0);
      intToBigInteger(&Mult4, 0);
      (void)BigIntMultiplyPower2(&common.k.sumSquares.productOtherDivisors,
        pstFactor->multiplicity / 2);
      if ((pstFactor->multiplicity & 1) != 0)
      {
        common.k.sumSquares.twoOddMultiplicity = true;
      }
    }
    else
    { /* Prime not 2 */
      int lenBytes = NumberLength * (int)sizeof(limb);
      NumberLength = valueP.nbrLimbs;
      (void)memcpy(&TestNbr, valueP.limbs, lenBytes);
      TestNbr[NumberLength].x = 0;
      GetMontgomeryParms(NumberLength);
      (void)memset(bigBase.limbs, 0, lenBytes);
      if ((valueP.limbs[0].x & 3) == 1)
      { /* if p = 1 (mod 4) */
        ComputeSumOfTwoSquaresForPrime();
        intToBigInteger(&Mult3, 0);
        intToBigInteger(&Mult4, 0);
        if (index4k1 < MAX_NBR_PRIMES_4K1)
        {
          CopyBigInt(&common.k.sumSquares.primeDivisors[index4k1].real, &Mult1);
          CopyBigInt(&common.k.sumSquares.primeDivisors[index4k1].imag, &Mult2);
          common.k.sumSquares.indexes[index4k1] = indexPrimes + 1;
          index4k1++;
          if ((pstFactor->multiplicity > 1) || (index4k1 > 1))
          {  // Not only one prime of the form 4k+1 with multiplicity 1.
            common.k.sumSquares.manyDecompositions = true;
          }
        }
        else
        {
          (void)BigIntPowerIntExp(&valueP, pstFactor->multiplicity, &Mult1);
          (void)BigIntMultiply(&common.k.sumSquares.productOtherDivisors,
            &Mult1, &common.k.sumSquares.productOtherDivisors);
        }
      } /* end p = 1 (mod 4) */
      else
      { /* if p = 3 (mod 4) */
        (void)BigIntPowerIntExp(&valueP, pstFactor->multiplicity / 2, &Mult1);
        (void)BigIntMultiply(&common.k.sumSquares.productOtherDivisors,
          &Mult1, &common.k.sumSquares.productOtherDivisors);
        ComputeSumOfFourSquaresForPrime();
      } /* end if p = 3 (mod 4) */
    } /* end prime not 2 */
    if ((pstFactor->multiplicity % 2) != 0)
    {
      GenerateSumOfFourSquaresOfPQ();
    }
    pstFactor++;
  } /* end for indexPrimes */
  common.k.sumSquares.nbrIndexes = index4k1;
  pstFactor = pstFactors + 1;      // Point to first factor in array of factors.
  for (indexPrimes = pstFactors->multiplicity - 1; indexPrimes >= 0; indexPrimes--)
  {
    NumberLength = *pstFactor->ptrFactor;
    IntArray2BigInteger(pstFactor->ptrFactor, &valueP);
    (void)BigIntPowerIntExp(&valueP, pstFactor->multiplicity / 2, &Tmp1);
    (void)BigIntMultiply(&Quad1, &Tmp1, &Quad1);
    (void)BigIntMultiply(&Quad2, &Tmp1, &Quad2);
    (void)BigIntMultiply(&Quad3, &Tmp1, &Quad3);
    (void)BigIntMultiply(&Quad4, &Tmp1, &Quad4);
    pstFactor++;
  }
  Quad1.sign = SIGN_POSITIVE;
  Quad2.sign = SIGN_POSITIVE;
  Quad3.sign = SIGN_POSITIVE;
  Quad4.sign = SIGN_POSITIVE;
  sortSquares();
}

static void varSquared(char** pptrOutput, char letter, char sign)
{
  char* ptrOutput = *pptrOutput;
  *ptrOutput = ' ';
  ptrOutput++;
  *ptrOutput = letter;
  ptrOutput++;
  copyStr(&ptrOutput, (prettyprint ? "&sup2;" : "^2"));
  *ptrOutput = ' ';
  ptrOutput++;
  *ptrOutput = sign;
  ptrOutput++;
  *pptrOutput = ptrOutput;
}

static void valueVar(char** pptrOutput, char letter, const BigInteger* value)
{
  char* ptrOutput = *pptrOutput;
  beginLine(&ptrOutput);
  *ptrOutput = letter;
  ptrOutput++;
  *ptrOutput = ' ';
  ptrOutput++;
  *ptrOutput = '=';
  ptrOutput++;
  *ptrOutput = ' ';
  ptrOutput++;
  if (hexadecimal)
  {
    BigInteger2Hex(&ptrOutput, value, groupLen);
  }
  else
  {
    BigInteger2Dec(&ptrOutput, value, groupLen);
  }
  finishLine(&ptrOutput);
  *pptrOutput = ptrOutput;
}

// Find product of (c1r + c1i * i) times (c2r + c2i * i):
// The real part is: c1r * c2r - c1i * c2i.
// The imaginary part is: c1r * c2i + c1i * c2r.
// When signIsPositive is false, change sign of c2i.
static void complexMultiply(struct stBigComplex* pstComp1,
  const struct stBigComplex* pstComp2, bool signIsPositive)
{
  (void)BigIntMultiply(&pstComp1->real, &pstComp2->real,
    &common.k.sumSquares.temp1);
  (void)BigIntMultiply(&pstComp1->imag, &pstComp2->imag,
    &common.k.sumSquares.temp2);
  if (signIsPositive)
  {
    BigIntSubt(&common.k.sumSquares.temp1, &common.k.sumSquares.temp2,
      &common.k.sumSquares.temp1);
  }
  else
  {
    BigIntAdd(&common.k.sumSquares.temp1, &common.k.sumSquares.temp2,
      &common.k.sumSquares.temp1);
  }
  (void)BigIntMultiply(&pstComp1->real, &pstComp2->imag,
    &common.k.sumSquares.temp2);
  (void)BigIntMultiply(&pstComp1->imag, &pstComp2->real,
    &common.k.sumSquares.temp3);
  CopyBigInt(&pstComp1->real, &common.k.sumSquares.temp1);
  if (signIsPositive)
  {
    BigIntAdd(&common.k.sumSquares.temp3, &common.k.sumSquares.temp2,
      &pstComp1->imag);
  }
  else
  {
    BigIntSubt(&common.k.sumSquares.temp3, &common.k.sumSquares.temp2,
      &pstComp1->imag);
  }
}

// Find quotient of (c1r + c1i * i) divided by (c2r + c2i * i):
// The norm of c2 is: c2r^2 + c2i^2
// The real part is: (c1r * c2r + c1i * c2i) / norm.
// The imaginary part is: (c1r * c2i - c1i * c2r) / norm.
static void complexDivide(struct stBigComplex* pstComp1,
  const struct stBigComplex* pstComp2, bool signIsPositive)
{
  // Compute dividends.
  complexMultiply(pstComp1, pstComp2, !signIsPositive);
  // Compute norm (divisor).
  (void)BigIntMultiply(&pstComp2->real, &pstComp2->real,
    &common.k.sumSquares.temp1);
  (void)BigIntMultiply(&pstComp2->imag, &pstComp2->imag,
    &common.k.sumSquares.temp2);
  BigIntAdd(&common.k.sumSquares.temp1, &common.k.sumSquares.temp2,
    &common.k.sumSquares.temp1);
  // Perform divisions.
  (void)BigIntDivide(&pstComp1->real, &common.k.sumSquares.temp1,
    &pstComp1->real);
  (void)BigIntDivide(&pstComp1->imag, &common.k.sumSquares.temp1,
    &pstComp1->imag);
}

static void initializeSumOfTwoSquares(void)
{
  int nbrIndexes2 = 0;
  int nbrExponents = common.k.sumSquares.nbrIndexes;
  // Compute sum of squares using divisor = 1.
  CopyBigInt(&common.k.sumSquares.currentValue.real,
    &common.k.sumSquares.productOtherDivisors);
  if (common.k.sumSquares.twoOddMultiplicity)
  {    // Number 2 with odd multiplicity.
    CopyBigInt(&common.k.sumSquares.currentValue.imag,
      &common.k.sumSquares.currentValue.real);
  }
  else
  {    // Number 2 with even multiplicity.
    intToBigInteger(&common.k.sumSquares.currentValue.imag, 0);
  }
  for (int exponentNbr = 0; exponentNbr < nbrExponents; exponentNbr++)
  {     // Multiply by power of M + Ni where M^2 + N^2 = 4k+1.
    bool significantBits;
    int expon;
    int index = common.k.sumSquares.indexes[exponentNbr];
    const struct sFactors* pstFactors = &astFactorsMod[index];
    struct stBigComplex* ptrDivisor = &common.k.sumSquares.divisors[exponentNbr];
    const struct stBigComplex* ptrPrimeDivisor = &common.k.sumSquares.primeDivisors[exponentNbr];
    CopyBigInt(&ptrDivisor->real, &ptrPrimeDivisor->real);
    CopyBigInt(&ptrDivisor->imag, &ptrPrimeDivisor->imag);
    expon = pstFactors->multiplicity - common.k.sumSquares.currentExp[exponentNbr];
    significantBits = false;
    if (expon == 0)
    {
      continue;
    }
    common.k.sumSquares.indexes2[nbrIndexes2] = index;
    common.k.sumSquares.indexes2toIndexes[nbrIndexes2] = exponentNbr;
    nbrIndexes2++;
    for (int mask = 0x1000000; mask > 0; mask /= 2)
    {
      if (significantBits)
      {
        complexMultiply(ptrDivisor, ptrDivisor, true);
      }
      if (expon & mask)
      {
        if (significantBits)
        {
          complexMultiply(ptrDivisor, ptrPrimeDivisor, true);
        }
        significantBits = true;
      }
    }
    complexMultiply(&common.k.sumSquares.currentValue,
      ptrDivisor, true);
  }
  common.k.sumSquares.GrayCode = 0;
  common.k.sumSquares.nbrIndexes2 = nbrIndexes2;
}

static void showSumSqButton(char** pptrOutput)
{
#ifndef __EMSCRIPTEN__
  (void)pptrOutput;
#endif
  if (common.k.sumSquares.manyDecompositions)
  {    // Sum of two squares. Show button and init sum of two squares data.
    common.k.sumSquares.initPending = true;
#ifdef __EMSCRIPTEN__
    // Show button.
    copyStr(pptrOutput,
      "<div id=\"sumSquares\"><p><button type=\"button\" id=\"showSumSq\">");
    if (lang)
    {
      copyStr(pptrOutput, "Todas las sumas de dos cuadrados");
    }
    else
    {
      copyStr(pptrOutput, "All sums of two squares");
    }
    copyStr(pptrOutput, "</button> <span class=\"new\">");
    copyStr(pptrOutput, lang ? "¡Nuevo!" : "New!");
    copyStr(pptrOutput, "</span></p></div>");
#endif
  }
}

void ShowFourSquares(char** pptrOutput)
{
  beginLine(pptrOutput);
  copyStr(pptrOutput, "n =");
  if (BigIntIsZero(&Quad4))
  {          // Quad4 equals zero.
    if (BigIntIsZero(&Quad3))
    {        // Quad3 and Quad4 equal zero.
      if (BigIntIsZero(&Quad2))
      {      // Quad2, Quad3 and Quad4 equal zero.
        varSquared(pptrOutput, 'a', ' ');
        finishLine(pptrOutput);
        valueVar(pptrOutput, 'a', &Quad1);
        showSumSqButton(pptrOutput);
        return;
      }
      varSquared(pptrOutput, 'a', '+');
      varSquared(pptrOutput, 'b', ' ');
      finishLine(pptrOutput);
      valueVar(pptrOutput, 'a', &Quad1);
      valueVar(pptrOutput, 'b', &Quad2);
      showSumSqButton(pptrOutput);
      return;
    }
    varSquared(pptrOutput, 'a', '+');
    varSquared(pptrOutput, 'b', '+');
    varSquared(pptrOutput, 'c', ' ');
    finishLine(pptrOutput);
    valueVar(pptrOutput, 'a', &Quad1);
    valueVar(pptrOutput, 'b', &Quad2);
    valueVar(pptrOutput, 'c', &Quad3);
    return;
  }
  varSquared(pptrOutput, 'a', '+');
  varSquared(pptrOutput, 'b', '+');
  varSquared(pptrOutput, 'c', '+');
  varSquared(pptrOutput, 'd', ' ');
  copyStr(pptrOutput, "</p>");
  valueVar(pptrOutput, 'a', &Quad1);
  valueVar(pptrOutput, 'b', &Quad2);
  valueVar(pptrOutput, 'c', &Quad3);
  valueVar(pptrOutput, 'd', &Quad4);
}

// In order to minimize multiplications of big numbers, Gray code is used.
// Let S = 2^r / product of all divisors of the form 4k+3.
// Divide the number by different even powers of divisors of the form 4k+1.
// Calling this product T, we get: R = N / S / T.
// Then R = (4k_1+1)^e_1 * (4k_2+1)^e_2 * ... * (4k_n+1)^e_n.
// Let 4k_s+1 = u_s^2 + v_s^2.
// Compute R = M^2 + N^2, where M + Ni = prod (u_s +/- i*v_s)^e_s.
// Use Gray code to select between plus and minus sign. The sign of the
// first element of the product will be always plus.
// Select M and N so they are both positive and M >= N.

void showSumTwoSquares(void)
{
  bool showMoreSumSquares = true;
  int nbrDivisors;
  int sumSquaresNbr;
  int nbrExponents = common.k.sumSquares.nbrIndexes;
  int* ptrFoundSumSquares = common.k.sumSquares.foundSumSquares;
  char* ptrOutput = output;
  if (common.k.sumSquares.initPending)
  {
    (void)memset(&common.k.sumSquares.currentExp, 0,
      sizeof(common.k.sumSquares.currentExp));
    (void)memset(&common.k.sumSquares.currentExpGray, 0,
      sizeof(common.k.sumSquares.currentExpGray));
    common.k.sumSquares.initPending = false;
    initializeSumOfTwoSquares();
  }
  *ptrOutput = 'T';      // Indicate this output is the sum of squares.
  ptrOutput++;
  copyStr(&ptrOutput, lang ? "<p>Suma de dos cuadrados:</p><ul>" :
    "<p>Sum of two squares:</p><ul>");
  for (sumSquaresNbr = 0; sumSquaresNbr < 1000; sumSquaresNbr++)
  {
    int exponentNbr;
    int mask;
    if (ptrFoundSumSquares > &common.k.sumSquares.foundSumSquares[900000])
    {
      break;             // No more space for holding sum of squares.
    }
    // Save current sum of squares into output array.
    // Save first the largest coefficient (real or imaginary).
    // If coefficient is negative, convert it to positive.
    CopyBigInt(&Tmp, &common.k.sumSquares.currentValue.real);
    Tmp.sign = SIGN_POSITIVE;
    CopyBigInt(&Tmp1, &common.k.sumSquares.currentValue.imag);
    Tmp1.sign = SIGN_POSITIVE;
    BigIntSubt(&Tmp, &Tmp1, &Tmp2);
    common.k.sumSquares.ptrFoundSumSquares[sumSquaresNbr] = ptrFoundSumSquares;
    for (int coeffNbr = 0; coeffNbr < 2; coeffNbr++)
    {  // 0: select maximum coefficient, 1: select minimum coefficient
      int arrLen;
      const BigInteger* ptrCurrentCoeff;
      if (coeffNbr == ((Tmp2.sign == SIGN_POSITIVE) ? 0 : 1))
      {   // Select real coefficient
        ptrCurrentCoeff = &Tmp;
      }
      else
      {   // Select imaginary coefficient.
        ptrCurrentCoeff = &Tmp1;
      }
      NumberLength = ptrCurrentCoeff -> nbrLimbs;
      BigInteger2IntArray(ptrFoundSumSquares, ptrCurrentCoeff);
      arrLen = 1 + ptrCurrentCoeff->nbrLimbs;
      ptrFoundSumSquares += arrLen;
    }
    common.k.sumSquares.GrayCode++;
    // Get the index of lowest bit of GrayCode.
    mask = 1;
    for (exponentNbr = 0; exponentNbr < common.k.sumSquares.nbrIndexes2; exponentNbr++)
    {
      if (common.k.sumSquares.GrayCode & mask)
      {
        break;
      }
      mask *= 2;
    }
    if (exponentNbr < (common.k.sumSquares.nbrIndexes2 - 1))
    {   // Test whether bit changes from 0 to 1 or from 1 to 0.
      int index;
      int leastSignificantBit = common.k.sumSquares.GrayCode &
        -common.k.sumSquares.GrayCode;
      int grayCode = common.k.sumSquares.GrayCode ^
        (common.k.sumSquares.GrayCode >> 1);
      bool signIsPositive;
      if ((grayCode & leastSignificantBit) != 0)
      {
        signIsPositive = true;
      }
      else
      {
        signIsPositive = false;
      }
      index = common.k.sumSquares.indexes2toIndexes[exponentNbr];
      complexDivide(&common.k.sumSquares.currentValue,
        &common.k.sumSquares.divisors[index], signIsPositive);
      complexMultiply(&common.k.sumSquares.currentValue,
        &common.k.sumSquares.divisors[index], !signIsPositive);
      continue;
    }
    // Found all sums of squares of previous divisor.
    // Find next divisor of number to be factored.
    for (exponentNbr = 0; exponentNbr < nbrExponents; exponentNbr++)
    {
      int index = common.k.sumSquares.indexes[exponentNbr];
      const struct sFactors* pstFactors = &astFactorsMod[index];
      if (pstFactors->multiplicity < 2)
      {
        continue;
      }
      if (common.k.sumSquares.currentExpGray[exponentNbr] == 0)
      {  // Ascending exponents.
        if ((common.k.sumSquares.currentExp[exponentNbr] + 1) <
          pstFactors->multiplicity)
        {  // Multiply by square of divisor of the form 4k+1.
          common.k.sumSquares.currentExp[exponentNbr] += 2;
          NumberLength = *pstFactors->ptrFactor;
          IntArray2BigInteger(pstFactors->ptrFactor, &Tmp);
          (void)BigIntMultiply(&common.k.sumSquares.productOtherDivisors, &Tmp,
            &common.k.sumSquares.productOtherDivisors);
          initializeSumOfTwoSquares();
          break;
        }
      }
      else
      {  // Descending exponents.
        if (common.k.sumSquares.currentExp[exponentNbr] > 1)
        {    // Divide by square of divisor of the form 4k+1.
          common.k.sumSquares.currentExp[exponentNbr] -= 2;
          NumberLength = *pstFactors->ptrFactor;
          IntArray2BigInteger(pstFactors->ptrFactor, &Tmp);
          (void)BigIntDivide(&common.k.sumSquares.productOtherDivisors, &Tmp,
            &common.k.sumSquares.productOtherDivisors);
          initializeSumOfTwoSquares();
          break;
        }
      }
      common.k.sumSquares.currentExpGray[exponentNbr] ^= 1;
    }
    if (exponentNbr == nbrExponents)
    {           // All exponents processed.
      showMoreSumSquares = false;
      sumSquaresNbr++;
      break;
    }
  }
  // Sort factors.
  nbrDivisors = sumSquaresNbr;
  heapsort(nbrDivisors, common.k.sumSquares.ptrFoundSumSquares);
  // Generate output from sorted sum of squares.
  // Both components are in consecutive addresses.
  for (sumSquaresNbr = 0; sumSquaresNbr < nbrDivisors; sumSquaresNbr++)
  {
    const int* ptrIntArray = common.k.sumSquares.ptrFoundSumSquares[sumSquaresNbr];
    int arrLen;
    copyStr(&ptrOutput, "<li>");
    for (int component = 0; component < 2; component++)
    {
      NumberLength = *ptrIntArray;
      IntArray2BigInteger(ptrIntArray, &Tmp);
      arrLen = *ptrIntArray;
      ptrIntArray++;
      if (hexadecimal)
      {
        Bin2Hex(&ptrOutput, (const limb*)ptrIntArray, arrLen, groupLen);
      }
      else
      {
        Bin2Dec(&ptrOutput, (const limb*)ptrIntArray, arrLen, groupLen);
      }
      copyStr(&ptrOutput, "&sup2;");
      if (component == 0)
      {
        copyStr(&ptrOutput, " + ");
        ptrIntArray += arrLen;
      }
    }
    copyStr(&ptrOutput, "</li>");
  }
  copyStr(&ptrOutput, "</ul>");
  if (showMoreSumSquares)
  {
#ifdef __EMSCRIPTEN__
    copyStr(&ptrOutput, "<p><button type=\"button\" id=\"showSumSq\">");
    copyStr(&ptrOutput, lang ? "Más sumas de dos cuadrados" : "More sums of two squares");
    copyStr(&ptrOutput, "</button></p>");
#endif
    output[0] = 'S';    // Indicate button present.
  }
}


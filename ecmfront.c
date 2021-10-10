//
// This file is part of Alpertron Calculators.
//
// Copyright 2016-2021 Dario Alejandro Alpern
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
#include "batch.h"
#include "output.h"
#include "commonstruc.h"

#ifdef FACTORIZATION_APP
#ifdef __EMSCRIPTEN__
extern bool skipPrimality;
extern int64_t lModularMult;
#endif
extern bool fromFile;
extern bool lineEndingCRLF;
extern BigInteger tofactor;
static BigInteger Quad1;
static BigInteger Quad2;
static BigInteger Quad3;
static BigInteger Quad4;
extern BigInteger factorValue;
static BigInteger result;
static BigInteger p;
static BigInteger q;
static BigInteger K;
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
static void ComputeFourSquares(const struct sFactors *pstFactors);
static void GetEulerTotient(char **pptrOutput);
static void GetMobius(char **pptrOutput);
static void GetNumberOfDivisors(char **pptrOutput);
static void GetSumOfDivisors(char **pptrOutput);
static void ShowFourSquares(char **pptrOutput);
static bool doFactorization;
static char *knownFactors;

void showDivisors(void);    // DEBUG BORRAR
#ifdef FACTORIZATION_APP
void batchCallback(char **pptrOutput)
{
  char *ptrFactorDec = tofactorDec;
  NumberLength = tofactor.nbrLimbs;
  if (tofactor.sign == SIGN_NEGATIVE)
  {
    *ptrFactorDec = '-';
    ptrFactorDec++;
  }
  BigInteger2IntArray(nbrToFactor, &tofactor);
  if (*nbrToFactor < 0)
  {    // If number is negative, make it positive.
    *nbrToFactor = -*nbrToFactor;
  }
  if (hexadecimal)
  {
    Bin2Hex(&ptrFactorDec, tofactor.limbs, tofactor.nbrLimbs, groupLen);
  }
  else
  {
    Bin2Dec(&ptrFactorDec, tofactor.limbs, tofactor.nbrLimbs, groupLen);
  }
  if (doFactorization)
  {
    factorExt(&tofactor, nbrToFactor, factorsMod, astFactorsMod, knownFactors);
    knownFactors = NULL;
  }
  SendFactorizationToOutput(astFactorsMod, pptrOutput, doFactorization);
}
#endif

#ifdef __EMSCRIPTEN__
static void startList(char** pptrOutput)
{
  if (!fromFile)
  {
    copyStr(pptrOutput, "<ul>");
  }
}

static void startListLine(char** pptrOutput)
{
  if (!fromFile)
  {
    copyStr(pptrOutput, "<li>");
  }
}

static void endListLine(char** pptrOutput)
{
  if (fromFile)
  {
    copyStr(pptrOutput, (lineEndingCRLF ? "\r\n" : "\n"));
  }
  else
  {
    copyStr(pptrOutput, "</li>");
  }
}

static void endList(char** pptrOutput)
{
  if (!fromFile)
  {
    copyStr(pptrOutput, "</ul>");
  }
}
#endif

static void ExponentToBigInteger(int exponent, BigInteger *bigint)
{
  if (exponent > (int)MAX_VALUE_LIMB)
  {
    bigint->limbs[0].x = exponent - (int)MAX_VALUE_LIMB;
    bigint->limbs[1].x = 1;
    bigint->nbrLimbs = 2;
  }
  else
  {
    bigint->limbs[0].x = exponent + 1;
    bigint->nbrLimbs = 1;
  }
  bigint->sign = SIGN_POSITIVE;
}

// Find number of divisors as the product of all exponents plus 1.
static void GetNumberOfDivisors(char **pptrOutput)
{
  char *ptrOutput = *pptrOutput;
  const struct sFactors *pstFactor;
  result.limbs[0].x = 1;    // Set result to 1.
  result.nbrLimbs = 1;
  result.sign = SIGN_POSITIVE;
  pstFactor = &astFactorsMod[1];
  for (int factorNumber = 1; factorNumber <= astFactorsMod[0].multiplicity; factorNumber++)
  {
    const int *ptrFactor = pstFactor->ptrFactor;
    if ((*ptrFactor == 1) && (*(ptrFactor + 1) < 2))
    {                        // Factor is 1.
      break;
    }
    ExponentToBigInteger(pstFactor->multiplicity, &factorValue);
    (void)BigIntMultiply(&factorValue, &result, &result);
    pstFactor++;
  }
  beginLine(&ptrOutput);
  copyStr(&ptrOutput, lang ? "Cantidad de divisores: " : "Number of divisors: ");
  if (hexadecimal)
  {
    BigInteger2Hex(&ptrOutput, &result, groupLen);
  }
  else
  {
    BigInteger2Dec(&ptrOutput, &result, groupLen);
  }
  finishLine(&ptrOutput);
  *pptrOutput = ptrOutput;
}

static void GetSumOfDivisors(char **pptrOutput)
{
  char *ptrOutput = *pptrOutput;
  SumOfDivisors(&result);
  beginLine(&ptrOutput);
  copyStr(&ptrOutput, lang ? "Suma de divisores: " : "Sum of divisors: ");
  if (hexadecimal)
  {
    BigInteger2Hex(&ptrOutput, &result, groupLen);
  }
  else
  {
    BigInteger2Dec(&ptrOutput, &result, groupLen);
  }
  finishLine(&ptrOutput);
  *pptrOutput = ptrOutput;
}

static void GetEulerTotient(char **pptrOutput)
{
  char *ptrOutput = *pptrOutput;
  Totient(&result);
  beginLine(&ptrOutput);
  copyStr(&ptrOutput, lang ? "Phi de Euler: " : "Euler's totient: ");
  if (hexadecimal)
  {
    BigInteger2Hex(&ptrOutput, &result, groupLen);
  }
  else
  {
    BigInteger2Dec(&ptrOutput, &result, groupLen);
  }
  finishLine(&ptrOutput);
  *pptrOutput = ptrOutput;
}

// Find Mobius as zero if some exponent is > 1, 1 if the number of factors is even, -1 if it is odd.
static void GetMobius(char **pptrOutput)
{
  char *ptrOutput = *pptrOutput;
  const struct sFactors *pstFactor;
  int mobius = 1;
  pstFactor = &astFactorsMod[1];
  if ((astFactorsMod[0].multiplicity > 1) || (*pstFactor->ptrFactor != 1) ||
    (*(pstFactor->ptrFactor + 1) != 1))
  {                                // Number to factor is not 1.
    for (int factorNumber = 1; factorNumber <= astFactorsMod[0].multiplicity; factorNumber++)
    {
      if (pstFactor->multiplicity == 1)
      {
        mobius = -mobius;
      }
      else
      {
        mobius = 0;
      }
      pstFactor++;
    }
  }
  beginLine(&ptrOutput);
  copyStr(&ptrOutput, "Möbius: ");
  if (mobius < 0)
  {
    mobius = -mobius;
    *ptrOutput = '-';
    ptrOutput++;
  }
  int2dec(&ptrOutput, mobius);
  finishLine(&ptrOutput);
  *pptrOutput = ptrOutput;
}

static void modPowShowStatus(const limb *base, const limb *exp, int nbrGroupsExp, limb *power)
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
  const struct sFactors *pstFactor = pstFactors + 1; // Point to first factor in array of factors.
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
static bool TryToFindSqRootMinus1(BigInteger *pTmp1, BigInteger *pTmp2,
  BigInteger *pTmp3, BigInteger *pTmp4)
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

static bool DivideCandidateBySmallPrimes(BigInteger* pCandidate, int **pptrArrFactors)
{
  int prime = 2;
  int primeIndex = 0;
  int *ptrArrFactors = *pptrArrFactors;
  do
  {
    int expon = 0;
    // Divide Tmp2 by small primes.
    while (getRemainder(pCandidate, prime) == 0)
    {
      subtractdivide(pCandidate, 0, prime);
      expon++;
    }
    if ((((expon % 2) == 1) && ((prime % 4) == 3)) || ((pCandidate->limbs[0].x % 4) == 3))
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
static void ComputeThreeSquares(BigInteger *pTmp,
  BigInteger *pTmp1, BigInteger *pTmp2, BigInteger *pTmp3, BigInteger *pTmp4,
  BigInteger *pM1, BigInteger *pM2)
{
  int nbrLimbs;
  int arrFactors[400];
  int diff = 1;
  int shRight;
  int *ptrArrFactors;
  CopyBigInt(pTmp, &tofactor);
  DivideBigNbrByMaxPowerOf4(&shRight, pTmp->limbs, &pTmp->nbrLimbs);
  diff = 0;
  for (;;)
  {
    const int* ptrArrFactorsBak;
    diff++;
    ptrArrFactors = arrFactors;
    intToBigInteger(pTmp1, diff*diff);
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
    for (ptrArrFactors = arrFactors; ptrArrFactors < ptrArrFactorsBak; ptrArrFactors+=2)
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

// Compute prime p as Mult1^2 + Mult2^2.
static void ComputeSumOfTwoSquaresForPrime(void)
{
  static limb minusOneMont[MAX_LEN];
  int lenBytes = NumberLength * (int)sizeof(limb);
  (void)memset(minusOneMont, 0, lenBytes);
  SubtBigNbrModN(minusOneMont, MontgomeryMultR1, minusOneMont, TestNbr, NumberLength);
  CopyBigInt(&q, &p);
  subtractdivide(&q, 1, 4);     // q = (prime-1)/4
  K.limbs[0].x = 1;
  do
  {    // Loop that finds mult1 = sqrt(-1) mod prime in Montgomery notation.
    K.limbs[0].x++;
    modPowShowStatus(K.limbs, q.limbs, q.nbrLimbs, Mult1.limbs);
  } while (!memcmp(Mult1.limbs, MontgomeryMultR1, lenBytes) ||
    !memcmp(Mult1.limbs, minusOneMont, lenBytes));
  Mult1.sign = SIGN_POSITIVE;
  lenBytes = p.nbrLimbs * (int)sizeof(limb);
  (void)memset(Mult2.limbs, 0, lenBytes);
  Mult2.limbs[0].x = 1;
  Mult2.nbrLimbs = 1;
  Mult2.sign = SIGN_POSITIVE;
  // Convert Mult1 to standard notation by multiplying by 1 in
  // Montgomery notation.
  modmult(Mult1.limbs, Mult2.limbs, Mult3.limbs);
  (void)memcpy(Mult1.limbs, Mult3.limbs, lenBytes);
  for (Mult1.nbrLimbs = p.nbrLimbs; Mult1.nbrLimbs > 1; Mult1.nbrLimbs--)
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
  Mult3.nbrLimbs = p.nbrLimbs;
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
  while (BigIntJacobiSymbol(&Tmp, &p) <= 0)
  {     // Not a quadratic residue. Compute next value of -1 - Mult1^2 in variable Tmp.
    BigIntAdd(&Tmp, &Mult2, &Tmp);
    Mult2.limbs[0].x += 2;
    mult1++;
  }
  // After the loop finishes, Tmp = (-1 - Mult1^2) is a quadratic residue mod p.
  // Convert base to Montgomery notation.
  BigIntAdd(&Tmp, &p, &Tmp);
  Tmp.limbs[NumberLength].x = 0;
  modmult(Tmp.limbs, MontgomeryMultR2, Tmp.limbs);

  intToBigInteger(&Mult1, mult1);
  CopyBigInt(&q, &p);
  subtractdivide(&q, -1, 4);  // q <- (p+1)/4.
  // Find Mult2 <- square root of Tmp = Tmp^q (mod p) in Montgomery notation.
  modPowShowStatus(Tmp.limbs, q.limbs, p.nbrLimbs, Mult2.limbs);
  // Convert Mult2 from Montgomery notation to standard notation.
  lenBytes = p.nbrLimbs * (int)sizeof(limb);
  (void)memset(Tmp.limbs, 0, lenBytes);
  Tmp.limbs[0].x = 1;
  intToBigInteger(&Mult3, 1);
  intToBigInteger(&Mult4, 0);
  // Convert Mult2 to standard notation by multiplying by 1 in
  // Montgomery notation.
  modmult(Mult2.limbs, Tmp.limbs, Mult2.limbs);
  for (Mult2.nbrLimbs = p.nbrLimbs; Mult2.nbrLimbs > 1; Mult2.nbrLimbs--)
  {  // Adjust number of limbs so the most significant limb is not zero.
    if (Mult2.limbs[Mult2.nbrLimbs - 1].x != 0)
    {
      break;
    }
  }
  Mult2.sign = SIGN_POSITIVE;
  MultiplyQuaternionBy2(&Mult1, &Mult2, &Mult3, &Mult4);
  CopyBigInt(&M1, &p);
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
  // Compute Tmp1 <- Mult1*Quad1 + Mult2*Quad2 + Mult3*Quad3 + Mult4*Quad4
  (void)BigIntMultiply(&Mult1, &Quad1, &Tmp);
  (void)BigIntMultiply(&Mult2, &Quad2, &Tmp4);
  BigIntAdd(&Tmp, &Tmp4, &Tmp);
  (void)BigIntMultiply(&Mult3, &Quad3, &Tmp4);
  BigIntAdd(&Tmp, &Tmp4, &Tmp);
  (void)BigIntMultiply(&Mult4, &Quad4, &Tmp4);
  BigIntAdd(&Tmp, &Tmp4, &Tmp1);

  // Compute Tmp2 <- Mult1*Quad2 - Mult2*Quad1 + Mult3*Quad4 - Mult4*Quad3
  (void)BigIntMultiply(&Mult1, &Quad2, &Tmp);
  (void)BigIntMultiply(&Mult2, &Quad1, &Tmp4);
  BigIntSubt(&Tmp, &Tmp4, &Tmp);
  (void)BigIntMultiply(&Mult3, &Quad4, &Tmp4);
  BigIntAdd(&Tmp, &Tmp4, &Tmp);
  (void)BigIntMultiply(&Mult4, &Quad3, &Tmp4);
  BigIntSubt(&Tmp, &Tmp4, &Tmp2);

  // Compute Tmp3 <- Mult1*Quad3 - Mult3*Quad1 - Mult2*Quad4 + Mult4*Quad2
  (void)BigIntMultiply(&Mult1, &Quad3, &Tmp);
  (void)BigIntMultiply(&Mult3, &Quad1, &Tmp4);
  BigIntSubt(&Tmp, &Tmp4, &Tmp);
  (void)BigIntMultiply(&Mult2, &Quad4, &Tmp4);
  BigIntSubt(&Tmp, &Tmp4, &Tmp);
  (void)BigIntMultiply(&Mult4, &Quad2, &Tmp4);
  BigIntAdd(&Tmp, &Tmp4, &Tmp3);

  // Compute Quad4 <- Mult1*Quad4 - Mult4*Quad1 + Mult2*Quad3 - Mult3*Quad2
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

static void ComputeFourSquares(const struct sFactors *pstFactors)
{
  int indexPrimes;
  const struct sFactors *pstFactor;

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
    {                            // Number to factor is 1.
      return;
    }
    if (*(pstFactor->ptrFactor + 1) == 0)
    {                              // Number to factor is 0.
      intToBigInteger(&Quad1, 0);  // 0 = 0^2 + 0^2 + 0^2 + 0^2
      return;
    }
  }
  for (indexPrimes = pstFactors -> multiplicity - 1; indexPrimes >= 0; indexPrimes--)
  {
    if ((pstFactor -> multiplicity % 2) == 0)
    {                              // Prime factor appears twice.
      pstFactor++;
      continue;
    }
    NumberLength = *pstFactor->ptrFactor;
    IntArray2BigInteger(pstFactor->ptrFactor, &p);
    CopyBigInt(&q, &p);
    addbigint(&q, -1);             // q <- p-1
    if ((p.nbrLimbs == 1) && (p.limbs[0].x == 2))
    {
      intToBigInteger(&Mult1, 1);  // 2 = 1^2 + 1^2 + 0^2 + 0^2
      intToBigInteger(&Mult2, 1);
      intToBigInteger(&Mult3, 0);
      intToBigInteger(&Mult4, 0);
    }
    else
    { /* Prime not 2 */
      int lenBytes = NumberLength * (int)sizeof(limb);
      NumberLength = p.nbrLimbs;
      (void)memcpy(&TestNbr, p.limbs, lenBytes);
      TestNbr[NumberLength].x = 0;
      GetMontgomeryParms(NumberLength);
      (void)memset(K.limbs, 0, lenBytes);
      if ((p.limbs[0].x & 3) == 1)
      { /* if p = 1 (mod 4) */
        ComputeSumOfTwoSquaresForPrime();
        intToBigInteger(&Mult3, 0);
        intToBigInteger(&Mult4, 0);
      } /* end p = 1 (mod 4) */
      else
      { /* if p = 3 (mod 4) */
        ComputeSumOfFourSquaresForPrime();
      } /* end if p = 3 (mod 4) */
    } /* end prime not 2 */
    GenerateSumOfFourSquaresOfPQ();
    pstFactor++;
  } /* end for indexPrimes */
  pstFactor = pstFactors + 1;      // Point to first factor in array of factors.
  for (indexPrimes = pstFactors->multiplicity - 1; indexPrimes >= 0; indexPrimes--)
  {
    NumberLength = *pstFactor->ptrFactor;
    IntArray2BigInteger(pstFactor->ptrFactor, &p);
    (void)BigIntPowerIntExp(&p, pstFactor->multiplicity / 2, &K);
    (void)BigIntMultiply(&Quad1, &K, &Quad1);
    (void)BigIntMultiply(&Quad2, &K, &Quad2);
    (void)BigIntMultiply(&Quad3, &K, &Quad3);
    (void)BigIntMultiply(&Quad4, &K, &Quad4);
    pstFactor++;
  }
  Quad1.sign = SIGN_POSITIVE;
  Quad2.sign = SIGN_POSITIVE;
  Quad3.sign = SIGN_POSITIVE;
  Quad4.sign = SIGN_POSITIVE;
  sortSquares();
}

static void varSquared(char **pptrOutput, char letter, char sign)
{
  char *ptrOutput = *pptrOutput;
  *ptrOutput = ' ';
  ptrOutput++;
  *ptrOutput = letter;
  ptrOutput++;
  copyStr(&ptrOutput, (prettyprint? "&sup2;": "^2"));
  *ptrOutput = ' ';
  ptrOutput++;
  *ptrOutput = sign;
  ptrOutput++;
  *pptrOutput = ptrOutput;
}

static void valueVar(char **pptrOutput, char letter, const BigInteger *value)
{
  char *ptrOutput = *pptrOutput;
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

static void ShowFourSquares(char **pptrOutput)
{
  char *ptrOutput = *pptrOutput;
  beginLine(&ptrOutput);
  copyStr(&ptrOutput, "n =");
  if (BigIntIsZero(&Quad4))
  {          // Quad4 equals zero.
    if (BigIntIsZero(&Quad3))
    {        // Quad3 and Quad4 equal zero.
      if (BigIntIsZero(&Quad2))
      {      // Quad2, Quad3 and Quad4 equal zero.
        varSquared(&ptrOutput, 'a', ' ');
        finishLine(&ptrOutput);
        valueVar(&ptrOutput, 'a', &Quad1);
        *pptrOutput = ptrOutput;
        return;
      }
      varSquared(&ptrOutput, 'a', '+');
      varSquared(&ptrOutput, 'b', ' ');
      finishLine(&ptrOutput);
      valueVar(&ptrOutput, 'a', &Quad1);
      valueVar(&ptrOutput, 'b', &Quad2);
      *pptrOutput = ptrOutput;
      return;
    }
    varSquared(&ptrOutput, 'a', '+');
    varSquared(&ptrOutput, 'b', '+');
    varSquared(&ptrOutput, 'c', ' ');
    finishLine(&ptrOutput);
    valueVar(&ptrOutput, 'a', &Quad1);
    valueVar(&ptrOutput, 'b', &Quad2);
    valueVar(&ptrOutput, 'c', &Quad3);
    *pptrOutput = ptrOutput;
    return;
  }
  varSquared(&ptrOutput, 'a', '+');
  varSquared(&ptrOutput, 'b', '+');
  varSquared(&ptrOutput, 'c', '+');
  varSquared(&ptrOutput, 'd', ' ');
  copyStr(&ptrOutput, "</p>");
  valueVar(&ptrOutput, 'a', &Quad1);
  valueVar(&ptrOutput, 'b', &Quad2);
  valueVar(&ptrOutput, 'c', &Quad3);
  valueVar(&ptrOutput, 'd', &Quad4);
  *pptrOutput = ptrOutput;
}

void ecmFrontText(char *tofactorText, bool performFactorization, char *factors)
{
  char *ptrOutput;
  bool isBatch;
  knownFactors = factors;
  if (valuesProcessed == 0)
  {
    doFactorization = performFactorization;
  }
  enum eExprErr rc = BatchProcessing(tofactorText, &tofactor, &ptrOutput, &isBatch);
  if (!isBatch && (rc == EXPR_OK) && (counterC == -1) && doFactorization)
  {
#ifdef __EMSCRIPTEN__
    int64_t sumSquaresModMult = 0;
#endif
    if (tofactor.sign == SIGN_POSITIVE)
    {        // Number to factor is non-negative.
#ifdef __EMSCRIPTEN__
      char* ptrText;
#endif
      if (!BigIntIsZero(&tofactor))
      {      // Number to factor is not zero.
        GetNumberOfDivisors(&ptrOutput);
        GetSumOfDivisors(&ptrOutput);
        GetEulerTotient(&ptrOutput);
        GetMobius(&ptrOutput);
      }
#ifdef __EMSCRIPTEN__
      StepECM = 3;   // Show progress (in percentage) of sum of squares.
      ptrText = ShowFactoredPart(&tofactor, astFactorsMod);
      beginLine(&ptrText);
      copyStr(&ptrText, lang ? "<p>Hallando suma de cuadrados.</p>" :
        "<p>Searching for sum of squares.</p>");
      ShowLowerText();
      sumSquaresModMult = lModularMult;
#endif
      ComputeFourSquares(astFactorsMod);
      ShowFourSquares(&ptrOutput);
#ifdef __EMSCRIPTEN__
      sumSquaresModMult = lModularMult - sumSquaresModMult;
      StepECM = 0;   // Do not show progress.
#endif
    }
#ifdef __EMSCRIPTEN__
    copyStr(&ptrOutput, "<div id=\"divisors\"><p><button type=\"button\" id=\"showdiv\">");
    copyStr(&ptrOutput, lang ? "Mostrar divisores": "Show divisors");
    copyStr(&ptrOutput, "</button></p>");
#endif
    beginLine(&ptrOutput);
    showElapsedTime(&ptrOutput);
    finishLine(&ptrOutput);
#ifdef __EMSCRIPTEN__
    if (lModularMult > 0)
    {
      beginLine(&ptrOutput);
      copyStr(&ptrOutput, lang ? "Multiplicaciones modulares:" :
        "Modular multiplications:");
      finishLine(&ptrOutput);
      startList(&ptrOutput);
      if ((lModularMult - primeModMult - SIQSModMult - sumSquaresModMult) > 0)
      {
        startListLine(&ptrOutput);
        copyStr(&ptrOutput, "ECM: ");
        long2dec(&ptrOutput, lModularMult - primeModMult - SIQSModMult - sumSquaresModMult);
        endListLine(&ptrOutput);
      }
      if (primeModMult > 0)
      {
        startListLine(&ptrOutput);
        copyStr(&ptrOutput, lang ? "Verificación de números primos probables: " :
          "Probable prime checking: ");
        long2dec(&ptrOutput, primeModMult);
        endListLine(&ptrOutput);
      }
      if (SIQSModMult > 0)
      {
        startListLine(&ptrOutput);
        copyStr(&ptrOutput, "SIQS: ");
        long2dec(&ptrOutput, SIQSModMult);
        endListLine(&ptrOutput);
      }
      if (sumSquaresModMult > 0)
      {
        startListLine(&ptrOutput);
        copyStr(&ptrOutput, lang ? "Suma de cuadrados: " : "Sum of squares: ");
        long2dec(&ptrOutput, sumSquaresModMult);
        endListLine(&ptrOutput);
      }
      endList(&ptrOutput);
    }
    if (nbrSIQS > 0)
    {
      uint64_t quotient;
      beginLine(&ptrOutput);
      copyStr(&ptrOutput, "SIQS:");
      startList(&ptrOutput);
      startListLine(&ptrOutput);
      int2dec(&ptrOutput, polynomialsSieved);
      copyStr(&ptrOutput, lang ? " polinomios utilizados" : " polynomials sieved");
      endListLine(&ptrOutput);
      startListLine(&ptrOutput);
      int2dec(&ptrOutput, trialDivisions);
      copyStr(&ptrOutput, lang ? " conjuntos de divisiones de prueba" : " sets of trial divisions");
      endListLine(&ptrOutput);
      startListLine(&ptrOutput);
      int2dec(&ptrOutput, smoothsFound);
      copyStr(&ptrOutput, lang ? " congruencias completas (1 de cada " :
        " smooth congruences found (1 out of every ");
      quotient = ValuesSieved / (uint64_t)smoothsFound;
      int2dec(&ptrOutput, (int)quotient);
      copyStr(&ptrOutput, lang ? " valores)</li><li>" : " values)</li><li>");
      int2dec(&ptrOutput, totalPartials);
      copyStr(&ptrOutput, lang ? " congruencias parciales (1 de cada " :
        " partial congruences found (1 out of every ");
      quotient = ValuesSieved / (uint64_t)totalPartials;
      int2dec(&ptrOutput, (int)quotient);
      copyStr(&ptrOutput, lang ? " valores)" : " values)");
      endListLine(&ptrOutput);
      startListLine(&ptrOutput);
      int2dec(&ptrOutput, partialsFound);
      copyStr(&ptrOutput, lang ? " congruencias parciales útiles</li><li>Tamaño de la matriz binaria: " :
        " useful partial congruences</li><li>Size of binary matrix: ");
      int2dec(&ptrOutput, matrixRows);
      copyStr(&ptrOutput, " &times; ");
      int2dec(&ptrOutput, matrixCols);
      endListLine(&ptrOutput);
      endList(&ptrOutput);
    }
    if ((nbrSIQS > 0) || (nbrECM > 0) || (nbrPrimalityTests > 0))
    {
      beginLine(&ptrOutput);
      copyStr(&ptrOutput, lang ? "Tiempos:" : "Timings:");
      startList(&ptrOutput);
      if (nbrPrimalityTests > 0)
      {
        startListLine(&ptrOutput);
        copyStr(&ptrOutput, lang ? "Test de primo probable de " : "Probable prime test of ");
        int2dec(&ptrOutput, nbrPrimalityTests);
        copyStr(&ptrOutput, lang ? " número" : " number");
        if (nbrPrimalityTests != 1)
        {
          *ptrOutput = 's';
          ptrOutput++;
        }
        *ptrOutput = ':';
        ptrOutput++;
        *ptrOutput = ' ';
        ptrOutput++;
        GetDHMSt(&ptrOutput, timePrimalityTests);
        endListLine(&ptrOutput);
      }
      if (nbrECM > 0)
      {
        startListLine(&ptrOutput);
        copyStr(&ptrOutput, lang ? "Factorización " : "Factoring ");
        int2dec(&ptrOutput, nbrECM);
        copyStr(&ptrOutput, lang ? " número" : " number");
        if (nbrECM != 1)
        {
          *ptrOutput = 's';
          ptrOutput++;
        }
        copyStr(&ptrOutput, lang ? " mediante ECM" : " using ECM:");
        *ptrOutput = ' ';
        ptrOutput++;
        GetDHMSt(&ptrOutput, timeECM - timeSIQS);
        endListLine(&ptrOutput);
      }
      if (nbrSIQS > 0)
      {
        startListLine(&ptrOutput);
        copyStr(&ptrOutput, lang ? "Factorización " : "Factoring ");
        int2dec(&ptrOutput, nbrSIQS);
        copyStr(&ptrOutput, lang ? " número" : " number");
        if (nbrSIQS != 1)
        {
          *ptrOutput = 's';
          ptrOutput++;
        }
        copyStr(&ptrOutput, lang ? " mediante SIQS" : " using SIQS:");
        *ptrOutput = ' ';
        ptrOutput++;
        GetDHMSt(&ptrOutput, timeSIQS);
        endListLine(&ptrOutput);
      }
      endList(&ptrOutput);
    }
    copyStr(&ptrOutput, "</div>");
#endif
  }
  beginLine(&ptrOutput);
  copyStr(&ptrOutput, lang ? COPYRIGHT_SPANISH : COPYRIGHT_ENGLISH );
  finishLine(&ptrOutput);
  // Initialize all exponents to zero.
  memset(common.divisors.currentExp, 0, sizeof(common.divisors.currentExp));
  memset(common.divisors.currentExpGray, 0, sizeof(common.divisors.currentExpGray));
  intToBigInteger(&common.divisors.divisor, 1);   // First divisor will be 1.
}

// Implementation of Heapsort as in Wikipedia.
static void siftDown(int start, int end)
{
  int root = start;
  while ((2 * root + 1) <= end)
  {
    int child = 2 * root + 1;
    int swap = root;
    // Compare a[swap] vs. a[child]
    if (IntArrayCompare(common.divisors.ptrFoundDivisors[swap],
      common.divisors.ptrFoundDivisors[child]) < 0)
    {        // a[swap] < a[child]
      swap = child;
    }
    if (((child + 1) <= end) &&
      (IntArrayCompare(common.divisors.ptrFoundDivisors[swap],
      common.divisors.ptrFoundDivisors[child + 1]) < 0))
    {             // a[swap] < a[child + 1]
      swap = child + 1;
    }
    if (swap == root)
    {
      return;
    }
    int* ptrTemp = common.divisors.ptrFoundDivisors[root];
    common.divisors.ptrFoundDivisors[root] = common.divisors.ptrFoundDivisors[swap];
    common.divisors.ptrFoundDivisors[swap] = ptrTemp;
    root = swap;
  }
}

static void heapify(int count)
{
  int start = (count - 2) / 2;
  while (start >= 0)
  {
    siftDown(start, count - 1);
    start--;
  }
}

static void heapsort(int count)
{
  int end = count - 1;
  heapify(count);
  while (end > 0)
  {
    int* ptrTemp = common.divisors.ptrFoundDivisors[end];
    common.divisors.ptrFoundDivisors[end] = common.divisors.ptrFoundDivisors[0];
    common.divisors.ptrFoundDivisors[0] = ptrTemp;
    end--;
    siftDown(0, end);
  }
}

// In order to minimize multiplications of big numbers, Gray code is used.
void showDivisors(void)
{
  bool showMoreDivisors = true;
  int nbrDivisors;
  int divisorNbr;
  int nbrExponents = astFactorsMod[0].multiplicity;
  int* ptrFoundDivisors = common.divisors.foundDivisors;
  char *ptrOutput = output;
  *ptrOutput = 'D';      // Indicate this output is the list of divisors.
  ptrOutput++;
  copyStr(&ptrOutput, lang ? "<p>Lista de divisores:</p><ul>" :
    "<p>List of divisors:</p><ul>");
  if ((tofactor.nbrLimbs == 1) && (tofactor.limbs[0].x <= 1))
  {
    if (tofactor.limbs[0].x == 0)
    {
      copyStr(&ptrOutput, lang? "<li>Cualquier número natural</li></ul>":
        "<li>Any natural number</li></ul>");
    }
    else
    {
      copyStr(&ptrOutput, "<li>1</li></ul>");
    }
    return;
  }
  if (nbrExponents > 50)
  {                      // Process only the first 50 exponents.
    nbrExponents = 50;
  }
  for (divisorNbr = 0; divisorNbr < 1000; divisorNbr++)
  {
    int exponentNbr;
    int arrLen;
    const struct sFactors* pstFactors;
    if (ptrFoundDivisors - &common.divisors.foundDivisors[0] > 900000)
    {
      break;             // Divisors are very large.
    }
    common.divisors.ptrFoundDivisors[divisorNbr] = ptrFoundDivisors;
    NumberLength = common.divisors.divisor.nbrLimbs;
    BigInteger2IntArray(ptrFoundDivisors, &common.divisors.divisor);
    arrLen = 1 + common.divisors.divisor.nbrLimbs;
    ptrFoundDivisors += arrLen;

    pstFactors = &astFactorsMod[1];
    // Find next divisor.
    for (exponentNbr = 0; exponentNbr < nbrExponents; exponentNbr++)
    {
      if (common.divisors.currentExpGray[exponentNbr] == 0)
      {  // Ascending exponents.
        if (common.divisors.currentExp[exponentNbr] <
          pstFactors->multiplicity)
        {
          common.divisors.currentExp[exponentNbr]++;
          NumberLength = *pstFactors->ptrFactor;
          IntArray2BigInteger(pstFactors->ptrFactor, &Tmp);
          BigIntMultiply(&common.divisors.divisor, &Tmp,
            &common.divisors.divisor);
          NumberLength = common.divisors.divisor.nbrLimbs;
          BigInteger2IntArray(ptrFoundDivisors, &common.divisors.divisor);
          break;
        }
      }
      else
      {  // Descending exponents.
        if (common.divisors.currentExp[exponentNbr] > 0)
        {
          common.divisors.currentExp[exponentNbr]--;
          NumberLength = *pstFactors->ptrFactor;
          IntArray2BigInteger(pstFactors->ptrFactor, &Tmp);
          BigIntDivide(&common.divisors.divisor, &Tmp,
            &common.divisors.divisor);
          NumberLength = common.divisors.divisor.nbrLimbs;
          BigInteger2IntArray(ptrFoundDivisors, &common.divisors.divisor);
          break;
        }
      }
      common.divisors.currentExpGray[exponentNbr] ^= 1;
      pstFactors++;
    }
    if (exponentNbr == nbrExponents)
    {           // All exponents processed.
      showMoreDivisors = false;
      divisorNbr++;
      break;
    }
  }
  // Sort factors.
  nbrDivisors = divisorNbr;
  heapsort(nbrDivisors);
  // Generate output from sorted divisors.
  for (divisorNbr = 0; divisorNbr < nbrDivisors; divisorNbr++)
  {
    const int* ptrIntArray = common.divisors.ptrFoundDivisors[divisorNbr];
    int arrLen;
    NumberLength = *common.divisors.ptrFoundDivisors[divisorNbr];
    IntArray2BigInteger(common.divisors.ptrFoundDivisors[divisorNbr], &Tmp);
    copyStr(&ptrOutput, "<li>");
    arrLen = *ptrIntArray;
    if (hexadecimal)
    {
      Bin2Hex(&ptrOutput, (const limb *)ptrIntArray + 1, arrLen, groupLen);
    }
    else
    {
      Bin2Dec(&ptrOutput, (const limb*)ptrIntArray + 1, arrLen, groupLen);
    }
    copyStr(&ptrOutput, "</li>");
  }
  copyStr(&ptrOutput, "</ul>");
  if (showMoreDivisors)
  {
    copyStr(&ptrOutput, "<p><button type=\"button\" id=\"showdiv\">");
    copyStr(&ptrOutput, lang ? "Mostrar más divisores" : "Show more divisors");
    copyStr(&ptrOutput, "</button></p>");
    output[0] = 'E';    // Indicate button present.
  }
}

#ifndef _MSC_VER
EXTERNALIZE void doWork(void)
{
  int flags;
  char *ptrData = inputString;
  char *ptrWebStorage;
  char *ptrKnownFactors;
#ifdef __EMSCRIPTEN__
  originalTenthSecond = tenths();
#endif
  if (*ptrData == 'C')
  {    // User pressed Continue button.
    ecmFrontText(NULL, false, NULL); // The 3rd parameter includes known factors.
#ifdef __EMSCRIPTEN__
    databack(output);
#endif
    return;
  }
  if (*ptrData == 'D')
  {    // User pressed SHow Divisors button.
    showDivisors();
#ifdef __EMSCRIPTEN__
    databack(output);
#endif
    return;
  }
  valuesProcessed = 0;
  groupLen = 0;
  while (*ptrData != ',')
  {
    groupLen = (groupLen * 10) + (*ptrData - '0');
    ptrData++;
  }
  ptrData++;             // Skip comma.
  flags = *ptrData;
  if (flags == '-')
  {
    ptrData++;
    flags = -*ptrData;
  }
#ifndef lang  
  lang = ((flags & 1)? true: false);
#endif
#ifdef __EMSCRIPTEN__
  if ((flags & (-2)) == '4')
  {
    skipPrimality = true;
    flags = 2;           // Do factorization.
  }
#endif
  ptrData += 2;          // Skip app number and second comma.
  fromFile = (*ptrData == '1');
  ptrData++;
  prettyprint = (*(ptrData + 1) == '1');
  if (fromFile)
  {
    prettyprint = false;
  }
  cunningham = (*(ptrData + 2) == '1');
  hexadecimal = (*(ptrData + 3) == '1');
  ptrData += 4;
  ptrWebStorage = ptrData + (int)strlen(ptrData) + 1;
  ptrKnownFactors = findChar(ptrWebStorage, '=');
  if (!prettyprint)
  {
    groupLen = -groupLen;  // Do not show number of digts.
  }
  if (ptrKnownFactors != NULL)
  {
    ptrKnownFactors++;
  }
  if ((flags & 0x80) && (ptrKnownFactors != NULL))
  {
    flags = 2;  // Do factorization.
  }
  if ((flags & 2) != 0)
  {               // Do factorization.
    ecmFrontText(ptrData, true, ptrKnownFactors); // The 3rd parameter includes known factors.
  }
  else
  {               // Do not perform factorization.
    ecmFrontText(ptrData, false, ptrKnownFactors); // The 3rd parameter includes known factors.
  }
#ifdef __EMSCRIPTEN__
  databack(output);
#endif
}
#endif
#endif
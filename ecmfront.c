/*
This file is part of Alpertron Calculators.

Copyright 2016 Dario Alejandro Alpern

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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bignbr.h"
#include "expression.h"
#include "factor.h"
#include "showtime.h"
#include "batch.h"

#ifdef __EMSCRIPTEN__
extern int skipPrimality;
extern long long lModularMult;
#endif
extern BigInteger tofactor;
static BigInteger Quad1, Quad2, Quad3, Quad4;
extern BigInteger factorValue;
static BigInteger result;
static void ComputeFourSquares(struct sFactors *pstFactors);
static void GetEulerTotient(char **pptrOutput);
static void GetMobius(char **pptrOutput);
static void GetNumberOfDivisors(char **pptrOutput);
static void GetSumOfDivisors(char **pptrOutput);
static void ShowFourSquares(char **pptrOutput);
static int doFactorization;
static char *knownFactors;

#ifdef FACTORIZATION_APP
void batchCallback(char **pptrOutput)
{
  char *ptrFactorDec = tofactorDec;
  NumberLength = tofactor.nbrLimbs;
  if (tofactor.sign == SIGN_NEGATIVE)
  {
    *ptrFactorDec++ = '-';
  }
  CompressBigInteger(nbrToFactor, &tofactor);
  if (hexadecimal)
  {
    Bin2Hex(tofactor.limbs, ptrFactorDec, tofactor.nbrLimbs, groupLen);
  }
  else
  {
    Bin2Dec(tofactor.limbs, ptrFactorDec, tofactor.nbrLimbs, groupLen);
  }
  if (doFactorization)
  {
    factorExt(&tofactor, nbrToFactor, factorsMod, astFactorsMod, knownFactors);
    knownFactors = NULL;
  }
  SendFactorizationToOutput(astFactorsMod, pptrOutput, doFactorization);
}
#endif

static void ExponentToBigInteger(int exponent, BigInteger *bigint)
{
  if (exponent > MAX_VALUE_LIMB)
  {
    bigint->limbs[0].x = exponent - MAX_VALUE_LIMB;
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
  struct sFactors *pstFactor;
  int factorNumber;
  result.limbs[0].x = 1;    // Set result to 1.
  result.nbrLimbs = 1;
  result.sign = SIGN_POSITIVE;
  pstFactor = &astFactorsMod[1];
  for (factorNumber = 1; factorNumber <= astFactorsMod[0].multiplicity; factorNumber++)
  {
    int *ptrFactor = pstFactor->ptrFactor;
    if (*ptrFactor == 1 && *(ptrFactor + 1) < 2)
    {                        // Factor is 1.
      break;
    }
    ExponentToBigInteger(pstFactor->multiplicity, &factorValue);
    BigIntMultiply(&factorValue, &result, &result);
    pstFactor++;
  }
  strcpy(ptrOutput, lang ? "<p>Cantidad de divisores: " : "<p>Number of divisors: ");
  ptrOutput += strlen(ptrOutput);
  if (hexadecimal)
  {
    BigInteger2Hex(&result, ptrOutput, groupLen);
  }
  else
  {
    BigInteger2Dec(&result, ptrOutput, groupLen);
  }
  ptrOutput += strlen(ptrOutput);
  strcpy(ptrOutput, "</p>");
  ptrOutput += strlen(ptrOutput);
  *pptrOutput = ptrOutput;
}

static void GetSumOfDivisors(char **pptrOutput)
{
  char *ptrOutput = *pptrOutput;
  SumOfDivisors(&result);
  strcpy(ptrOutput, lang ? "<p>Suma de divisores: " : "<p>Sum of divisors: ");
  ptrOutput += strlen(ptrOutput);
  if (hexadecimal)
  {
    BigInteger2Hex(&result, ptrOutput, groupLen);
  }
  else
  {
    BigInteger2Dec(&result, ptrOutput, groupLen);
  }
  ptrOutput += strlen(ptrOutput);
  strcpy(ptrOutput, "</p>");
  ptrOutput += strlen(ptrOutput);
  *pptrOutput = ptrOutput;
}

static void GetEulerTotient(char **pptrOutput)
{
  char *ptrOutput = *pptrOutput;
  Totient(&result);
  strcpy(ptrOutput, lang ? "<p>Phi de Euler: " : "<p>Euler's totient: ");
  ptrOutput += strlen(ptrOutput);
  if (hexadecimal)
  {
    BigInteger2Hex(&result, ptrOutput, groupLen);
  }
  else
  {
    BigInteger2Dec(&result, ptrOutput, groupLen);
  }
  ptrOutput += strlen(ptrOutput);
  strcpy(ptrOutput, "</p>");
  ptrOutput += strlen(ptrOutput);
  *pptrOutput = ptrOutput;
}

// Find Mobius as zero if some exponent is > 1, 1 if the number of factors is even, -1 if it is odd.
static void GetMobius(char **pptrOutput)
{
  char *ptrOutput = *pptrOutput;
  struct sFactors *pstFactor;
  int mobius = 1;
  pstFactor = &astFactorsMod[1];
  if (astFactorsMod[0].multiplicity > 1 || *pstFactor->ptrFactor != 1 ||
    *(pstFactor->ptrFactor + 1) != 1)
  {                                // Number to factor is not 1.
    int factorNumber;
    for (factorNumber = 1; factorNumber <= astFactorsMod[0].multiplicity; factorNumber++)
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
  strcpy(ptrOutput, "<p>Möbius: ");
  ptrOutput += strlen(ptrOutput);
  if (mobius < 0)
  {
    mobius = -mobius;
    *ptrOutput++ = '-';
  }
  int2dec(&ptrOutput, mobius);
  strcpy(ptrOutput, "</p>");
  ptrOutput += strlen(ptrOutput);
  *pptrOutput = ptrOutput;
}

static void modPowShowStatus(limb *base, limb *exp, int nbrGroupsExp, limb *power)
{
  int mask, index;
  memcpy(power, MontgomeryMultR1, (NumberLength + 1) * sizeof(*power));  // power <- 1
  for (index = nbrGroupsExp - 1; index >= 0; index--)
  {
    int groupExp = (int)(exp + index)->x;
#ifdef __EMSCRIPTEN__
    percentageBPSW = (nbrGroupsExp - index) * 100 / nbrGroupsExp;
#endif
    for (mask = 1 << (BITS_PER_GROUP - 1); mask > 0; mask >>= 1)
    {
      modmult(power, power, power);
      if ((groupExp & mask) != 0)
      {
        modmult(power, base, power);
      }
    }
  }
}

static void ComputeFourSquares(struct sFactors *pstFactors)
{
  int indexPrimes;
  static BigInteger p, q, K, Mult1, Mult2, Mult3, Mult4;
  static BigInteger Tmp, Tmp1, Tmp2, Tmp3, Tmp4, M1, M2, M3, M4, M5, M6, M7, M8;
  struct sFactors *pstFactor;
  static limb minusOneMont[MAX_LEN];

  intToBigInteger(&Quad1, 1);     // 1 = 1^2 + 0^2 + 0^2 + 0^2
  intToBigInteger(&Quad2, 0);
  intToBigInteger(&Quad3, 0);
  intToBigInteger(&Quad4, 0);
  pstFactor = pstFactors + 1;      // Point to first factor in array of factors.
  if (pstFactors->multiplicity == 1 && *pstFactor->ptrFactor == 1)
  {
    if (*(pstFactor->ptrFactor + 1) == 1)
    {                              // Number to factor is 1.
      return;
    }
    if (*(pstFactor->ptrFactor + 1) == 0)
    {                             // Number to factor is 0.
      intToBigInteger(&Quad1, 0);     // 0 = 0^2 + 0^2 + 0^2 + 0^2
      return;
    }
  }
  for (indexPrimes = pstFactors -> multiplicity - 1; indexPrimes >= 0; 
       indexPrimes--, pstFactor++)
  {
    if (pstFactor -> multiplicity % 2 == 0)
    {                              // Prime factor appears twice.
      continue;
    }
    NumberLength = *pstFactor->ptrFactor;
    UncompressBigInteger(pstFactor->ptrFactor, &p);
    CopyBigInt(&q, &p);
    addbigint(&q, -1);             // q <- p-1
    if (p.nbrLimbs == 1 && p.limbs[0].x == 2)
    {
      intToBigInteger(&Mult1, 1); // 2 = 1^2 + 1^2 + 0^2 + 0^2
      intToBigInteger(&Mult2, 1);
      intToBigInteger(&Mult3, 0);
      intToBigInteger(&Mult4, 0);
    }
    else
    { /* Prime not 2 */
      int nbrLimbs;
      NumberLength = p.nbrLimbs;
      memcpy(&TestNbr, p.limbs, NumberLength * sizeof(limb));
      TestNbr[NumberLength].x = 0;
      GetMontgomeryParms(NumberLength);
      memset(minusOneMont, 0, NumberLength * sizeof(limb));
      SubtBigNbrModN(minusOneMont, MontgomeryMultR1, minusOneMont, TestNbr, NumberLength);
      memset(K.limbs, 0, NumberLength * sizeof(limb));
      if ((p.limbs[0].x & 3) == 1)
      { /* if p = 1 (mod 4) */
        CopyBigInt(&q, &p);
        subtractdivide(&q, 1, 4);     // q = (prime-1)/4
        K.limbs[0].x = 1;
        do
        {    // Loop that finds mult1 = sqrt(-1) mod prime in Montgomery notation.
          K.limbs[0].x++;
          modPowShowStatus(K.limbs, q.limbs, q.nbrLimbs, Mult1.limbs);
        } while (!memcmp(Mult1.limbs, MontgomeryMultR1, NumberLength * sizeof(limb)) ||
          !memcmp(Mult1.limbs, minusOneMont, NumberLength * sizeof(limb)));
        Mult1.sign = SIGN_POSITIVE;
        memset(Mult2.limbs, 0, p.nbrLimbs * sizeof(limb));
        Mult2.limbs[0].x = 1;
        Mult2.nbrLimbs = 1;
        Mult2.sign = SIGN_POSITIVE;
        // Convert Mult1 to standard notation by multiplying by 1 in
        // Montgomery notation.
        modmult(Mult1.limbs, Mult2.limbs, Mult3.limbs);
        memcpy(Mult1.limbs, Mult3.limbs, p.nbrLimbs * sizeof(limb));
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
        memcpy(Mult3.limbs, TestNbr, p.nbrLimbs * sizeof(limb));
        Mult3.nbrLimbs = p.nbrLimbs;
        Mult3.sign = SIGN_POSITIVE;
        while (Mult3.nbrLimbs > 1 && Mult3.limbs[Mult3.nbrLimbs - 1].x == 0)
        {
          Mult3.nbrLimbs--;
        }
        intToBigInteger(&Mult4, 0);   // Initialize imaginary part to 0.
                                      // Find gcd of (Mult1 + Mult2 * i) and (Mult3 + Mult4 * i)
        GaussianGCD(&Mult1, &Mult2, &Mult3, &Mult4, &Tmp1, &Tmp2, &Tmp3, &Tmp4);
        nbrLimbs = Mult1.nbrLimbs;
        if (nbrLimbs < Mult2.nbrLimbs)
        {
          nbrLimbs = Mult2.nbrLimbs;
        }
        CopyBigInt(&Mult1, &Tmp1);
        CopyBigInt(&Mult2, &Tmp2);
        intToBigInteger(&Mult3, 0);
        intToBigInteger(&Mult4, 0);
      } /* end p = 1 (mod 4) */
      else
      { /* if p = 3 (mod 4) */
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
        memset(Tmp.limbs, 0, p.nbrLimbs * sizeof(limb));
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
      } /* end if p = 3 (mod 4) */
    } /* end prime not 2 */

    // Compute Tmp1 <- Mult1*Quad1 + Mult2*Quad2 + Mult3*Quad3 + Mult4*Quad4
    BigIntMultiply(&Mult1, &Quad1, &Tmp);
    BigIntMultiply(&Mult2, &Quad2, &Tmp4);
    BigIntAdd(&Tmp, &Tmp4, &Tmp);
    BigIntMultiply(&Mult3, &Quad3, &Tmp4);
    BigIntAdd(&Tmp, &Tmp4, &Tmp);
    BigIntMultiply(&Mult4, &Quad4, &Tmp4);
    BigIntAdd(&Tmp, &Tmp4, &Tmp1);

    // Compute Tmp2 <- Mult1*Quad2 - Mult2*Quad1 + Mult3*Quad4 - Mult4*Quad3
    BigIntMultiply(&Mult1, &Quad2, &Tmp);
    BigIntMultiply(&Mult2, &Quad1, &Tmp4);
    BigIntSubt(&Tmp, &Tmp4, &Tmp);
    BigIntMultiply(&Mult3, &Quad4, &Tmp4);
    BigIntAdd(&Tmp, &Tmp4, &Tmp);
    BigIntMultiply(&Mult4, &Quad3, &Tmp4);
    BigIntSubt(&Tmp, &Tmp4, &Tmp2);

    // Compute Tmp3 <- Mult1*Quad3 - Mult3*Quad1 - Mult2*Quad4 + Mult4*Quad2
    BigIntMultiply(&Mult1, &Quad3, &Tmp);
    BigIntMultiply(&Mult3, &Quad1, &Tmp4);
    BigIntSubt(&Tmp, &Tmp4, &Tmp);
    BigIntMultiply(&Mult2, &Quad4, &Tmp4);
    BigIntSubt(&Tmp, &Tmp4, &Tmp);
    BigIntMultiply(&Mult4, &Quad2, &Tmp4);
    BigIntAdd(&Tmp, &Tmp4, &Tmp3);

    // Compute Quad4 <- Mult1*Quad4 - Mult4*Quad1 + Mult2*Quad3 - Mult3*Quad2
    BigIntMultiply(&Mult1, &Quad4, &Tmp);
    BigIntMultiply(&Mult4, &Quad1, &Tmp4);
    BigIntSubt(&Tmp, &Tmp4, &Tmp);
    BigIntMultiply(&Mult2, &Quad3, &Tmp4);
    BigIntAdd(&Tmp, &Tmp4, &Tmp);
    BigIntMultiply(&Mult3, &Quad2, &Tmp4);
    BigIntSubt(&Tmp, &Tmp4, &Quad4);

    CopyBigInt(&Quad3, &Tmp3);
    CopyBigInt(&Quad2, &Tmp2);
    CopyBigInt(&Quad1, &Tmp1);
  } /* end for indexPrimes */
  pstFactor = pstFactors + 1;      // Point to first factor in array of factors.
  for (indexPrimes = pstFactors->multiplicity - 1; indexPrimes >= 0; indexPrimes--, pstFactor++)
  {
    NumberLength = *pstFactor->ptrFactor;
    UncompressBigInteger(pstFactor->ptrFactor, &p);
    BigIntPowerIntExp(&p, pstFactor->multiplicity / 2, &K);
    BigIntMultiply(&Quad1, &K, &Quad1);
    BigIntMultiply(&Quad2, &K, &Quad2);
    BigIntMultiply(&Quad3, &K, &Quad3);
    BigIntMultiply(&Quad4, &K, &Quad4);
  }
  Quad1.sign = SIGN_POSITIVE;
  Quad2.sign = SIGN_POSITIVE;
  Quad3.sign = SIGN_POSITIVE;
  Quad4.sign = SIGN_POSITIVE;
  // Sort squares
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

static void varSquared(char **pptrOutput, char letter, char sign)
{
  char *ptrOutput = *pptrOutput;
  *ptrOutput++ = ' ';
  *ptrOutput++ = letter;
  strcpy(ptrOutput, (prettyprint? "&sup2;": "^2"));
  ptrOutput += strlen(ptrOutput);
  *ptrOutput++ = ' ';
  *ptrOutput++ = sign;
  *pptrOutput = ptrOutput;
}

static void valueVar(char **pptrOutput, char letter, BigInteger *value)
{
  char *ptrOutput = *pptrOutput;
  strcpy(ptrOutput, "<p>");
  ptrOutput += strlen(ptrOutput);
  *ptrOutput++ = letter;
  *ptrOutput++ = ' ';
  *ptrOutput++ = '=';
  *ptrOutput++ = ' ';
  if (hexadecimal)
  {
    BigInteger2Hex(value, ptrOutput, groupLen);
  }
  else
  {
    BigInteger2Dec(value, ptrOutput, groupLen);
  }
  ptrOutput += strlen(ptrOutput);
  strcpy(ptrOutput, "</p>");
  ptrOutput += strlen(ptrOutput);
  *pptrOutput = ptrOutput;
}

static void ShowFourSquares(char **pptrOutput)
{
  char *ptrOutput = *pptrOutput;
  strcpy(ptrOutput, "<p>n =");
  ptrOutput += strlen(ptrOutput);
  if (Quad4.nbrLimbs == 1 && Quad4.limbs[0].x == 0)
  {          // Quad4 equals zero.
    if (Quad3.nbrLimbs == 1 && Quad3.limbs[0].x == 0)
    {        // Quad3 and Quad4 equal zero.
      if (Quad2.nbrLimbs == 1 && Quad2.limbs[0].x == 0)
      {      // Quad2, Quad3 and Quad4 equal zero.
        varSquared(&ptrOutput, 'a', ' ');
        strcpy(ptrOutput, "</p>");
        ptrOutput += strlen(ptrOutput);
        valueVar(&ptrOutput, 'a', &Quad1);
        *pptrOutput = ptrOutput;
        return;
      }
      varSquared(&ptrOutput, 'a', '+');
      varSquared(&ptrOutput, 'b', ' ');
      strcpy(ptrOutput, "</p>");
      ptrOutput += strlen(ptrOutput);
      valueVar(&ptrOutput, 'a', &Quad1);
      valueVar(&ptrOutput, 'b', &Quad2);
      *pptrOutput = ptrOutput;
      return;
    }
    varSquared(&ptrOutput, 'a', '+');
    varSquared(&ptrOutput, 'b', '+');
    varSquared(&ptrOutput, 'c', ' ');
    strcpy(ptrOutput, "</p>");
    ptrOutput += strlen(ptrOutput);
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
  strcpy(ptrOutput, "</p>");
  ptrOutput += strlen(ptrOutput);
  valueVar(&ptrOutput, 'a', &Quad1);
  valueVar(&ptrOutput, 'b', &Quad2);
  valueVar(&ptrOutput, 'c', &Quad3);
  valueVar(&ptrOutput, 'd', &Quad4);
  *pptrOutput = ptrOutput;
}

void ecmFrontText(char *tofactorText, int performFactorization, char *factors)
{
  char *ptrOutput, *ptrText;
  knownFactors = factors;
  if (valuesProcessed == 0)
  {
    doFactorization = performFactorization;
  }
  enum eExprErr rc = BatchProcessing(tofactorText, &tofactor, &ptrOutput);
  if (valuesProcessed == 1)
  {
    if (rc == EXPR_OK && doFactorization)
    {
      if (tofactor.sign == SIGN_POSITIVE)
      {        // Number to factor is non-negative.
        if (tofactor.nbrLimbs > 1 || tofactor.limbs[0].x > 0)
        {      // Number to factor is not zero.
          GetNumberOfDivisors(&ptrOutput);
          GetSumOfDivisors(&ptrOutput);
          GetEulerTotient(&ptrOutput);
          GetMobius(&ptrOutput);
        }
#ifdef __EMSCRIPTEN__
        StepECM = 3;   // Show progress (in percentage) of sum of squares.
        ptrText = ShowFactoredPart(&tofactor, astFactorsMod);
        strcpy(ptrText, lang ? "<p>Hallando suma de cuadrados.</p>" :
          "<p>Searching for sum of squares.</p>");
        ShowLowerText();
#endif
        ComputeFourSquares(astFactorsMod);
        ShowFourSquares(&ptrOutput);
#ifdef __EMSCRIPTEN__
        StepECM = 0;   // Do not show progress.
#endif
      }
      showElapsedTime(&ptrOutput);
    }
  }
  strcpy(ptrOutput, lang ? "<p>" COPYRIGHT_SPANISH "</p>" :
    "<p>" COPYRIGHT_ENGLISH "</p>");
}

void doWork(void)
{
  int flags;
  char *ptrData = inputString;
  char *ptrWebStorage, *ptrKnownFactors;
#ifdef __EMSCRIPTEN__
  originalTenthSecond = tenths();
#endif
  if (*ptrData == 'C')
  {    // User pressed Continue button.
    ecmFrontText(NULL, 0, NULL); // The 3rd parameter includes known factors.
#ifdef __EMSCRIPTEN__
    databack(output);
#endif
    return;
  }
  valuesProcessed = 0;
  groupLen = 0;
  while (*ptrData != ',')
  {
    groupLen = groupLen * 10 + (*ptrData++ - '0');
  }
  ptrData++;             // Skip comma.
  flags = *ptrData;
  if (flags == '-')
  {
    flags = -*(++ptrData);
  }
  lang = flags & 1;
#ifdef __EMSCRIPTEN__
  if ((flags & (-2)) == '4')
  {
    skipPrimality = TRUE;
    flags = 2;           // Do factorization.
  }
#endif
  ptrData += 2;          // Skip app number and second comma.
  verbose = (*(ptrData + 1) == '1');
  prettyprint = (*(ptrData + 2) == '1');
  cunningham = (*(ptrData + 3) == '1');
  hexadecimal = (*(ptrData + 4) == '1');
  ptrData += 5;
  ptrWebStorage = ptrData + strlen(ptrData) + 1;
  ptrKnownFactors = findChar(ptrWebStorage, '=');
  if (prettyprint == 0)
  {
    groupLen = -groupLen;  // Do not show number of digts.
  }
  if (ptrKnownFactors)
  {
    ptrKnownFactors++;
  }
  if (flags & 0x80)
  {
    if (ptrKnownFactors)
    {
      flags = 2;  // do factorization.
    }
  }
  ecmFrontText(ptrData, flags & 2, ptrKnownFactors); // The 3rd parameter includes known factors.
#ifdef __EMSCRIPTEN__
  databack(output);
#endif
}

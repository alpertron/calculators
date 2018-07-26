/*
This file is part of Alpertron Calculators.

Copyright 2017 Dario Alejandro Alpern

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

static BigInteger Solution1[400];
static BigInteger Solution2[400];
static BigInteger Increment[400];
static int Exponents[400];
static BigInteger Aux[400];
BigInteger ValA, ValB, ValC, ValN;
static int SolNbr;
static char textExp[1000];
static BigInteger LastModulus;
static int nbrFactors;
char *ptrOutput;

static int Show(BigInteger *num, char *str, int t)
{
  if (num->nbrLimbs > 1 || num->limbs[0].x != 0)
  {     // num is not zero.
    if ((t & 1) != 0 && num->sign == SIGN_POSITIVE)
    {
      *ptrOutput++ = ' ';
      *ptrOutput++ = '+';
    }
    if (num->sign == SIGN_NEGATIVE)
    {
      *ptrOutput++ = ' ';
      *ptrOutput++ = '-';
    }
    if (num->nbrLimbs != 1 || num->limbs[0].x != 1)
    {    // num is not 1 or -1.
      int signTemp = num->sign;
      num->sign = SIGN_POSITIVE;
      *ptrOutput++ = ' ';
      BigInteger2Dec(num, ptrOutput, groupLen);
      num->sign = signTemp;
      ptrOutput += strlen(ptrOutput);
    }
    strcpy(ptrOutput, str);
    ptrOutput += strlen(ptrOutput);
    return t | 1;
  }
  return t;
}

void Show1(BigInteger *num, int t)
{
  int u = Show(num, "", t);
  if ((u & 1) == 0 || (num->nbrLimbs == 1 && num->limbs[0].x == 1))
  {
    BigInteger2Dec(num, ptrOutput, groupLen);
    ptrOutput += strlen(ptrOutput);
  }
}

void Solution(BigInteger *value)
{
  SolNbr++;
  strcpy(ptrOutput, "<li>x = ");
  ptrOutput += strlen(ptrOutput);
  BigInteger2Dec(value, ptrOutput, groupLen);
  ptrOutput += strlen(ptrOutput);
  strcpy(ptrOutput, "</li>");
  ptrOutput += strlen(ptrOutput);
}

void SolveEquation(void)
{
  int factorIndex, expon, T1, E;
  BigInteger GcdAll;
  BigInteger ValNn;
  BigInteger z, Mult, currentSolution;
  BigInteger U, SqrtDisc, prime, K1;
  BigInteger discriminant, V, K, Q, L;
  struct sFactors *pstFactor;

  if (ValN.nbrLimbs == 1 && ValN.limbs[0].x == 0)
  {        // Mod 0 => Equation in integer numbers
    if (ValA.nbrLimbs == 1 && ValA.limbs[0].x == 0)
    {      // Not Quadratic Equation
      if (ValB.nbrLimbs == 1 && ValB.limbs[0].x == 0)
      {    // Constant Equation
        if (ValC.nbrLimbs == 1 && ValC.limbs[0].x == 0)
        {  // 0 = 0
          strcpy(ptrOutput, "<p>The equation is satisfied by any integer <var>x</var>.</p>");
          ptrOutput += strlen(ptrOutput);
        }
        else
        {
          return;
        }
      }
      else
      {        // Linear Equation: ValB * x + ValC = 0.
        BigIntRemainder(&ValC, &ValB, &Aux[0]);
        if (Aux[0].nbrLimbs == 1 && Aux[0].limbs[0].x == 0)
        {      // ValC is multiple of ValB: solution is -ValC / ValB.
          BigIntDivide(&ValC, &ValB, &Aux[0]);
          BigIntNegate(&Aux[0], &Aux[0]);
          Solution(&Aux[0]);
        }
        else
        {      // No solutions on integers.
          return;
        }
      }
    }
    else
    {                       // Quadratic Equation
      // The discriminant is ValB * ValB - 4 * ValA * ValC.
      BigIntMultiply(&ValB, &ValB, &Aux[0]);
      BigIntMultiply(&ValA, &ValC, &Aux[1]);
      multint(&Aux[1], &Aux[1], 4);
      BigIntSubt(&Aux[0], &Aux[1], &discriminant);
      if (discriminant.sign == SIGN_NEGATIVE)
      {        // No integer solutions.
        return;
      }
      // Set V to square root of discriminant.
      squareRoot(discriminant.limbs, V.limbs, discriminant.nbrLimbs, &V.nbrLimbs);
      V.sign = SIGN_POSITIVE;
      BigIntMultiply(&V, &V, &Aux[0]);
      BigIntSubt(&Aux[0], &discriminant, &Aux[0]);
      if (Aux[0].nbrLimbs != 1 || Aux[0].limbs[0].x != 0)
      {  // discriminant has no integer square root.
        return;
      }
      multint(&Aux[0], &ValA, 2);    // Denominator
      BigIntSubt(&V, &ValB, &Aux[1]);
      BigIntRemainder(&Aux[1], &Aux[0], &Aux[2]);
      if (Aux[2].nbrLimbs == 1 && Aux[2].limbs[0].x == 0)
      {      // (V-ValB)/(2*ValA) is integer: it is a solution.
        BigIntDivide(&Aux[1], &Aux[0], &Aux[2]);
        Solution(&Aux[2]);
      }
      BigIntNegate(&V, &V);
      BigIntSubt(&V, &ValB, &Aux[1]);
      BigIntRemainder(&Aux[1], &Aux[0], &Aux[2]);
      if (Aux[2].nbrLimbs == 1 && Aux[2].limbs[0].x == 0)
      {      // (-V-ValB)/(2*ValA) is integer: it is a solution.
        BigIntDivide(&Aux[1], &Aux[0], &Aux[2]);
        Solution(&Aux[2]);
      }
    }
    SolNbr = 1;
    return;
  }
  // Modulus is not zero.
  ValN.sign = SIGN_POSITIVE;
  BigIntRemainder(&ValA, &ValN, &ValA);
  BigIntRemainder(&ValB, &ValN, &ValB);
  BigIntRemainder(&ValC, &ValN, &ValC);
  BigIntGcd(&ValA, &ValB, &Aux[0]);
  BigIntGcd(&ValC, &Aux[0], &GcdAll);
  BigIntRemainder(&ValC, &GcdAll, &Aux[0]);
  if (Aux[0].nbrLimbs != 1 || Aux[0].limbs[0].x != 0)
  {  // ValC must be multiple of gcd(ValA, ValB).
     // Otherwise go out because there are no solutions.
    return;
  }
  BigIntGcd(&ValN, &GcdAll, &Aux[0]);
  CopyBigInt(&GcdAll, &Aux[0]);
  // Divide all coefficients by gcd(ValA, ValB).
  BigIntDivide(&ValA, &GcdAll, &Aux[0]);
  CopyBigInt(&ValA, &Aux[0]);
  BigIntDivide(&ValB, &GcdAll, &Aux[0]);
  CopyBigInt(&ValB, &Aux[0]);
  BigIntDivide(&ValC, &GcdAll, &Aux[0]);
  CopyBigInt(&ValC, &Aux[0]);
  BigIntDivide(&ValN, &GcdAll, &Aux[0]);
  CopyBigInt(&ValN, &Aux[0]);
  CopyBigInt(&ValNn, &ValN);
  if (ValNn.nbrLimbs == 1 && ValNn.limbs[0].x == 1)
  {     // All values from 0 to GcdAll - 1 are solutions.
    if (GcdAll.nbrLimbs > 1 || GcdAll.limbs[0].x > 5)
    {
      strcpy(ptrOutput, "<p>All values of <var>x</var> between 0 and ");
      ptrOutput += strlen(ptrOutput);
      addbigint(&GcdAll, -1);
      BigInteger2Dec(&GcdAll, ptrOutput, groupLen);
      ptrOutput += strlen(ptrOutput);
      strcpy(ptrOutput, " are solutions.</p>");
      ptrOutput += strlen(ptrOutput);
    }
    else
    {
      int ctr;
      for (ctr = 0; ctr < GcdAll.limbs[0].x; ctr++)
      {
        intToBigInteger(&Aux[0], ctr);
        Solution(&Aux[0]);
      }
    }
    return;
  }
  BigIntRemainder(&ValA, &ValN, &Aux[0]);
  if (Aux[0].nbrLimbs == 1 && Aux[0].limbs[0].x == 0)
  {           // Linear equation.
    BigIntGcd(&ValB, &ValN, &Aux[0]);
    if (Aux[0].nbrLimbs != 1 || Aux[0].limbs[0].x != 1)
    {         // ValB and ValN are not coprime. Go out.
      return;
    }
    // Calculate z <- -ValC / ValB (mod ValN)
    NumberLength = ValN.nbrLimbs;
    memcpy(TestNbr, ValN.limbs, NumberLength * sizeof(limb));
    TestNbr[NumberLength].x = 0;
    GetMontgomeryParms(NumberLength);
    BigIntModularDivision(&ValC, &ValB, &ValN, &z);
    BigIntNegate(&z, &z);
    if (z.sign == SIGN_NEGATIVE)
    {
      BigIntAdd(&z, &ValN, &z);
    }
    BigIntMultiply(&ValNn, &GcdAll, &Aux[0]);
    do
    {
      Solution(&z);
      BigIntAdd(&z, &ValN, &z);
      BigIntSubt(&z, &Aux[0], &Aux[1]);
    } while (Aux[1].sign == SIGN_NEGATIVE);
    return;
  }
  BigIntSubt(&LastModulus, &ValN, &Aux[0]);
  if (Aux[0].nbrLimbs > 1 || Aux[0].limbs[0].x != 0)
  {     // Last modulus is different from ValN.
    CopyBigInt(&LastModulus, &ValN);
    NumberLength = ValN.nbrLimbs;
    CompressBigInteger(nbrToFactor, &ValN);
    Bin2Dec(ValN.limbs, tofactorDec, ValN.nbrLimbs, groupLen);
    factor(&ValN, nbrToFactor, factorsMod, astFactorsMod);
  }
  intToBigInteger(&Q, 0);
  nbrFactors = astFactorsMod[0].multiplicity;
  pstFactor = &astFactorsMod[1];
  for (factorIndex = 0; factorIndex<nbrFactors; factorIndex++)
  {
    int sol1Invalid = 0;
    int sol2Invalid = 0;
    NumberLength = *pstFactor->ptrFactor;
    UncompressBigInteger(pstFactor->ptrFactor, &prime);
    expon = pstFactor->multiplicity;
    BigIntPowerIntExp(&prime, expon, &V);
    BigIntRemainder(&ValA, &V, &L);
    if (L.nbrLimbs == 1 && L.limbs[0].x == 0)
    {     // ValA multiple of prime, means linear equation mod prime.
      if (ValB.nbrLimbs == 1 && ValB.limbs[0].x == 0 &&
        !(ValC.nbrLimbs == 1 && ValC.limbs[0].x == 0))
      {
        return;  // There are no solutions: ValB=0 and ValC!=0
      }
      BigIntGcd(&ValB, &V, &U);
      BigIntDivide(&V, &U, &Q);
      BigIntDivide(&ValB, &U, &V);
      BigIntDivide(&ValC, &U, &L);
      BigIntNegate(&V, &V);
      if (prime.nbrLimbs == 1 && prime.limbs[0].x == 2)
      {      // If prime is 2.
        BigIntModularDivisionPower2(&L, &V, &Q, &Solution1[factorIndex]);
      }
      else
      {
        NumberLength = Q.nbrLimbs;
        memcpy(TestNbr, Q.limbs, NumberLength * sizeof(limb));
        TestNbr[NumberLength].x = 0;
        GetMontgomeryParms(NumberLength);
        BigIntModularDivision(&L, &V, &Q, &Solution1[factorIndex]);
      }
      CopyBigInt(&Solution2[factorIndex], &Solution1[factorIndex]);
    }
    else
    {                   /* If quadratic equation mod p */
      BigInteger squareRoot;
      BigInteger tmp1, tmp2, ValAOdd;
      int nbrBitsSquareRoot, correctBits, nbrLimbs, bitsAZero, mask, ctr;
      int deltaZeros, deltaIsZero = 0;
      // Compute discriminant = ValB^2 - 4*ValA*ValC.
      BigIntMultiply(&ValB, &ValB, &Aux[0]);
      BigIntMultiply(&ValA, &ValC, &discriminant);
      multint(&discriminant, &discriminant, 4);
      BigIntSubt(&Aux[0], &discriminant, &discriminant);
      CopyBigInt(&ValAOdd, &ValA);
      if (prime.nbrLimbs == 1 && prime.limbs[0].x == 2)
      {         /* Prime p is 2 */
        // Number of bits of square root of discriminant to compute: expon + bits_a + 1,
        // where bits_a is the number of least significant bits of a set to zero.
        // To compute the square root, compute the inverse of sqrt, so only multiplications are used.
        // f(x) = invsqrt(x), f_{n+1}(x) = f_n * (3 - x*f_n^2)/2
        // Get odd part of A and number of bits to zero.
        DivideBigNbrByMaxPowerOf2(&bitsAZero, ValAOdd.limbs, &ValAOdd.nbrLimbs);
        bitsAZero++;   // Include the number 2 in the denominator.
        // Compute inverse of -A (mod 2^expon).
        nbrLimbs = (expon + BITS_PER_GROUP - 1) / BITS_PER_GROUP;
        if (nbrLimbs > ValAOdd.nbrLimbs)
        {
          memset(&ValAOdd.limbs[ValAOdd.nbrLimbs], 0, (nbrLimbs - ValAOdd.nbrLimbs) * sizeof(limb));
        }
        if (ValAOdd.sign == ValB.sign)
        {
          ChSignBigNbr((int *)ValAOdd.limbs, nbrLimbs);
        }
        NumberLength = nbrLimbs;
        ComputeInversePower2(ValAOdd.limbs, tmp2.limbs, tmp1.limbs);
        memcpy(ValAOdd.limbs, tmp2.limbs, nbrLimbs*sizeof(limb));
        if (BigIntIsZero(&discriminant))
        {     // discriminant is zero.
          memset(squareRoot.limbs, 0, nbrLimbs*sizeof(limb));
          deltaIsZero = 1;
          deltaZeros = 0;
          nbrBitsSquareRoot = expon + bitsAZero;
          correctBits = expon / 2;
          if (correctBits == 0)
          {
            correctBits = 1;
          }
        }
        else
        {
          int LSLimb;
          CopyBigInt(&tmp1, &discriminant);
          DivideBigNbrByMaxPowerOf4(&deltaZeros, discriminant.limbs, &discriminant.nbrLimbs);
          // Find number of bits of square root to compute.
          nbrBitsSquareRoot = expon + bitsAZero - deltaZeros;
          nbrLimbs = (nbrBitsSquareRoot + BITS_PER_GROUP - 1) / BITS_PER_GROUP;
          if (nbrLimbs > discriminant.nbrLimbs)
          {
            memset(&discriminant.limbs[nbrLimbs], 0, (nbrLimbs - discriminant.nbrLimbs) * sizeof(limb));
          }
          if (discriminant.sign == SIGN_NEGATIVE)
          {
            ChSignBigNbr((int *)discriminant.limbs, nbrLimbs);
          }
          LSLimb = discriminant.limbs[0].x;
          if (expon >= deltaZeros && expon - deltaZeros < 8)
          {
            LSLimb &= (1 << (expon - deltaZeros)) - 1;
          }
          if (nbrBitsSquareRoot == 2 && (LSLimb & 3) != 1)
          {
            return;    // Square root does not exist. Go out.
          }
          if (nbrBitsSquareRoot > 2 && (LSLimb & 7) != 1)
          {
            return;    // Square root does not exist. Go out.
          }
          // First approximation to inverse of square root.
          squareRoot.limbs[0].x = ((LSLimb & 15) == 1 ? 1 : 3);
          correctBits = 2;
          while (correctBits < nbrBitsSquareRoot)
          {   // Compute f(x) = invsqrt(x), f_{n+1}(x) = f_n * (3 - x*f_n^2)/2
            correctBits *= 2;
            nbrLimbs = correctBits / BITS_PER_GROUP + 1;
            MultBigNbr((int *)squareRoot.limbs, (int *)squareRoot.limbs, (int *)tmp2.limbs, nbrLimbs);
            MultBigNbr((int *)tmp2.limbs, (int *)discriminant.limbs, (int *)tmp2.limbs, nbrLimbs);
            ChSignBigNbr((int *)tmp2.limbs, nbrLimbs);
            memset(tmp1.limbs, 0, nbrLimbs * sizeof(limb));
            tmp1.limbs[0].x = 3;
            AddBigNbr((int *)tmp1.limbs, (int *)tmp2.limbs, (int *)tmp2.limbs, nbrLimbs);
            MultBigNbr((int *)tmp2.limbs, (int *)squareRoot.limbs, (int *)tmp1.limbs, nbrLimbs);
            memcpy(squareRoot.limbs, tmp1.limbs, nbrLimbs * sizeof(limb));
            DivBigNbrByInt((int *)tmp1.limbs, 2, (int *)squareRoot.limbs, nbrLimbs);
          }
          // Get square root of discriminant from its inverse by multiplying by discriminant.
          MultBigNbr((int *)discriminant.limbs, (int *)squareRoot.limbs, (int *)tmp1.limbs, nbrLimbs);
          memcpy(squareRoot.limbs, tmp1.limbs, nbrLimbs*sizeof(limb));
          // Multiply by square root of discriminant by 2^deltaZeros.
          squareRoot.limbs[nbrLimbs].x = 0;
          for (ctr = 0; ctr < deltaZeros; ctr++)
          {
            MultBigNbrByInt((int*)squareRoot.limbs, 2, (int*)squareRoot.limbs, nbrLimbs + 1);
            if (squareRoot.limbs[nbrLimbs].x != 0)
            {
              nbrLimbs++;
              squareRoot.limbs[nbrLimbs].x = 0;
            }
          }
          squareRoot.sign = SIGN_POSITIVE;
          squareRoot.nbrLimbs = nbrLimbs;
          correctBits = expon - deltaZeros;
          if (nbrBitsSquareRoot < 2)
          {
            correctBits = nbrBitsSquareRoot;
          }
        }
        BigIntPowerOf2(&Q, correctBits);      // Store increment.
        mask = 1 << (correctBits % BITS_PER_GROUP);
        // Compute x = (b + sqrt(discriminant)) / (-2a)
        nbrLimbs = (nbrBitsSquareRoot + BITS_PER_GROUP - 1) / BITS_PER_GROUP;
        if (ValB.sign == SIGN_POSITIVE)
        {
          AddBigNbr((int *)ValB.limbs, (int *)squareRoot.limbs, (int *)tmp1.limbs, nbrLimbs);
        }
        else
        {
          SubtractBigNbr((int *)ValB.limbs, (int *)squareRoot.limbs, (int *)tmp1.limbs, nbrLimbs);
        }
        for (ctr = 0; ctr < bitsAZero; ctr++)
        {
          if (tmp1.limbs[0].x & 1)
          {   // Cannot divide by 2 because number is odd.
            sol1Invalid = 1;
            break;
          }
          DivBigNbrByInt((int *)tmp1.limbs, 2, (int *)tmp1.limbs, nbrLimbs);
        }
        nbrLimbs = (correctBits + BITS_PER_GROUP - 1) / BITS_PER_GROUP;
        MultBigNbr((int *)tmp1.limbs, (int *)ValAOdd.limbs, (int *)tmp2.limbs, nbrLimbs);
        if (mask != 1)
        {
          tmp2.limbs[nbrLimbs - 1].x &= mask - 1;
        }
        tmp2.nbrLimbs = nbrLimbs;
        tmp2.sign = SIGN_POSITIVE;
        CopyBigInt(&Solution1[factorIndex], &tmp2);
        // Compute x = (b - sqrt(discriminant)) / (-2a)
        nbrLimbs = (nbrBitsSquareRoot + BITS_PER_GROUP - 1) / BITS_PER_GROUP;
        if (ValB.sign == SIGN_NEGATIVE)
        {
          AddBigNbr((int *)ValB.limbs, (int *)squareRoot.limbs, (int *)tmp1.limbs, nbrLimbs);
        }
        else
        {
          SubtractBigNbr((int *)ValB.limbs, (int *)squareRoot.limbs, (int *)tmp1.limbs, nbrLimbs);
        }
        for (ctr = 0; ctr < bitsAZero; ctr++)
        {
          if (tmp1.limbs[0].x & 1)
          {   // Cannot divide by 2 because number is odd.
            sol2Invalid = 1;
            break;
          }
          DivBigNbrByInt((int *)tmp1.limbs, 2, (int *)tmp1.limbs, nbrLimbs);
        }
        nbrLimbs = (correctBits + BITS_PER_GROUP - 1) / BITS_PER_GROUP;
        MultBigNbr((int *)tmp1.limbs, (int *)ValAOdd.limbs, (int *)tmp2.limbs, nbrLimbs);
        if (mask != 1)
        {
          tmp2.limbs[nbrLimbs - 1].x &= (mask - 1);
        }
        tmp2.nbrLimbs = nbrLimbs;
        tmp2.sign = SIGN_POSITIVE;
        CopyBigInt(&Solution2[factorIndex], &tmp2);
      }
      else
      {                        // Prime is not 2
        // Number of bits of square root of discriminant to compute: expon + bits_a + 1,
        // where bits_a is the number of least significant bits of a set to zero.
        // To compute the square root, compute the inverse of sqrt, so only multiplications are used.
        // f(x) = invsqrt(x), f_{n+1}(x) = f_n * (3 - x*f_n^2)/2
        // Get maximum power of prime which divide ValA.
        bitsAZero = 0;
        for (;;)
        {
          BigIntRemainder(&ValAOdd, &prime, &tmp1);
          if (tmp1.nbrLimbs > 1 || tmp1.limbs[0].x != 0)
          {
            break;
          }
          BigIntDivide(&ValAOdd, &prime, &ValAOdd);
          bitsAZero++;
        }
        BigIntRemainder(&discriminant, &V, &discriminant);
        // Get maximum power of prime which divide discriminant.
        if (BigIntIsZero(&discriminant))
        {      // Discriminant is zero.
          deltaZeros = expon;
        }
        else
        {      // Discriminant is not zero.
          deltaZeros = 0;
          for (;;)
          {
            BigIntRemainder(&discriminant, &prime, &tmp1);
            if (tmp1.nbrLimbs > 1 || tmp1.limbs[0].x != 0)
            {
              break;
            }
            BigIntDivide(&discriminant, &prime, &discriminant);
            deltaZeros++;
          }
        }
        if ((deltaZeros & 1) && deltaZeros < expon)
        {          // If delta is of type m*prime^n where m is not multiple of prime
                   // and n is odd, there is no solution, so go out.
          return;
        }
        deltaZeros >>= 1;
        // Compute inverse of -2*A (mod prime^(expon - deltaZeros)).
        BigIntAdd(&ValAOdd, &ValAOdd, &ValAOdd);
        BigIntPowerIntExp(&prime, expon - deltaZeros, &tmp1);
        BigIntRemainder(&ValAOdd, &tmp1, &ValAOdd);
        nbrLimbs = tmp1.nbrLimbs;
        if (ValAOdd.sign == SIGN_NEGATIVE)
        {
          BigIntAdd(&tmp1, &ValAOdd, &ValAOdd);
        }
        else if (ValAOdd.nbrLimbs > 1 || ValAOdd.limbs[0].x != 0)
        {
          BigIntSubt(&tmp1, &ValAOdd, &ValAOdd);
        }
        intToBigInteger(&tmp2, 1);
        NumberLength = tmp1.nbrLimbs;
        memcpy(TestNbr, tmp1.limbs, NumberLength * sizeof(limb));
        TestNbr[NumberLength].x = 0;
        GetMontgomeryParms(NumberLength);
        BigIntModularDivision(&tmp2, &ValAOdd, &tmp1, &Aux[0]);
        CopyBigInt(&ValAOdd, &Aux[0]);
        if (discriminant.nbrLimbs == 1 && discriminant.limbs[0].x == 0)
        {     // Discriminant is zero.
          memset(squareRoot.limbs, 0, nbrLimbs * sizeof(limb));
          deltaIsZero = 1;
          nbrBitsSquareRoot = expon + bitsAZero;
        }
        else
        {      // Discriminant is not zero.
          // Find number of digits of square root to compute.
          nbrBitsSquareRoot = expon + bitsAZero - deltaZeros;
          BigIntPowerIntExp(&prime, nbrBitsSquareRoot, &tmp1);
          nbrLimbs = tmp1.nbrLimbs;
          BigIntRemainder(&discriminant, &tmp1, &discriminant);
          if (discriminant.sign == SIGN_NEGATIVE)
          {
            BigIntAdd(&discriminant, &tmp1, &discriminant);
          }
          if (nbrLimbs > discriminant.nbrLimbs)
          {
            memset(&discriminant.limbs[nbrLimbs], 0, (nbrLimbs - discriminant.nbrLimbs) * sizeof(limb));
          }
          BigIntRemainder(&discriminant, &prime, &Aux[3]);
          if (Aux[3].sign == SIGN_NEGATIVE)
          {
            BigIntAdd(&Aux[3], &prime, &Aux[3]);
          }
          if (BigIntJacobiSymbol(&Aux[3], &prime) != 1)
          {
            return;         // Not a quadratic residue, so go out.
          }
          // Compute square root of discriminant.
          NumberLength = prime.nbrLimbs;
          memcpy(TestNbr, prime.limbs, NumberLength * sizeof(limb));
          TestNbr[NumberLength].x = 0;
          GetMontgomeryParms(NumberLength);
          CopyBigInt(&Q, &prime);
          if ((prime.limbs[0].x & 3) == 3)
          {                 // prime mod 4 = 3
            subtractdivide(&Q, -1, 4);   // Q <- (prime+1)/4.
            BigIntModularPower(&Aux[3], &Q, &SqrtDisc);
          }
          else
          {
            limb *toConvert;
            // Convert discriminant to Montgomery notation.
            CompressLimbsBigInteger(Aux[5].limbs, &Aux[3]);
            modmult(Aux[5].limbs, MontgomeryMultR2, Aux[6].limbs);  // u
            if ((prime.limbs[0].x & 7) == 5)
            {              // prime mod 8 = 5: use Atkin's method for modular square roots.
               // Step 1. v <- (2u)^((p-5)/8) mod p
               // Step 2. i <- (2uv^2) mod p
               // Step 3. square root of u <- uv (i-1)
              // Step 1.
              subtractdivide(&Q, 5, 8);   // Q <- (prime-5)/8.
              AddBigNbrMod(Aux[6].limbs, Aux[6].limbs, Aux[7].limbs); // 2u
              modPow(Aux[7].limbs, Q.limbs, Q.nbrLimbs, Aux[8].limbs);
              // At this moment Aux[7].limbs is v in Montgomery notation.

              // Step 2.
              modmult(Aux[8].limbs, Aux[8].limbs, Aux[9].limbs);  // v^2
              modmult(Aux[7].limbs, Aux[9].limbs, Aux[9].limbs);  // i

              // Step 3.
              SubtBigNbrMod(Aux[9].limbs, MontgomeryMultR1, Aux[9].limbs); // i-1
              modmult(Aux[8].limbs, Aux[9].limbs, Aux[9].limbs);           // v*(i-1)
              modmult(Aux[6].limbs, Aux[9].limbs, Aux[9].limbs);           // u*v*(i-1)
              toConvert = Aux[9].limbs;
            }
            else
            {          // prime = 1 (mod 8). Use Shanks' method for modular square roots.
              // Step 1. Select e >= 3, q odd such that p = 2^e * q + 1.
              // Step 2. Choose x at random in the range 1 < x < p such that jacobi (x, p) = -1.
              // Step 3. Set z <- x^q (mod p).
              // Step 4. Set y <- z, r <- e, x <- a^((q-1)/2) mod p, v <- ax mod p, w <- vx mod p.
              // Step 5. If w = 1, the computation ends with +/-v as the square root.
              // Step 6. Find the smallest value of k such that w^(2^k) = 1 (mod p)
              // Step 7. Set d <- y^(2^(r-k-1)) mod p, y <- d^2 mod p, r <- k, v <- dv mod p, w <- wy mod p.
              // Step 8. Go to step 5.
              int x, e, r, k;
              // Step 1.
              subtractdivide(&Q, 1, 1);   // Q <- (prime-1).
              DivideBigNbrByMaxPowerOf2(&e, Q.limbs, &Q.nbrLimbs);
              // Step 2.
              x = 1;
              do
              {
                intToBigInteger(&Aux[3], ++x);
              } while (BigIntJacobiSymbol(&Aux[3], &prime) >= 0);
              // Step 3.
              // Get z <- x^q (mod p) in Montgomery notation.
              modPowBaseInt(x, Q.limbs, Q.nbrLimbs, Aux[4].limbs);  // z
              // Step 4.
              memcpy(Aux[5].limbs, Aux[4].limbs, NumberLength * sizeof(limb)); // y
              r = e;
              CopyBigInt(&K1, &Q);
              subtractdivide(&K1, 1, 2);
              modPow(Aux[6].limbs, K1.limbs, K1.nbrLimbs, Aux[7].limbs); // x
              modmult(Aux[6].limbs, Aux[7].limbs, Aux[8].limbs);         // v
              modmult(Aux[8].limbs, Aux[7].limbs, Aux[9].limbs);         // w
              // Step 5
              while (memcmp(Aux[9].limbs, MontgomeryMultR1, NumberLength * sizeof(limb)))
              {
                // Step 6
                k = 0;
                memcpy(Aux[10].limbs, Aux[9].limbs, NumberLength * sizeof(limb));
                do
                {
                  k++;
                  modmult(Aux[10].limbs, Aux[10].limbs, Aux[10].limbs);
                } while (memcmp(Aux[10].limbs, MontgomeryMultR1, NumberLength * sizeof(limb)));
                // Step 7
                memcpy(Aux[11].limbs, Aux[5].limbs, NumberLength * sizeof(limb)); // d
                for (ctr = 0; ctr < r - k - 1; ctr++)
                {
                  modmult(Aux[11].limbs, Aux[11].limbs, Aux[11].limbs);
                }
                modmult(Aux[11].limbs, Aux[11].limbs, Aux[5].limbs);   // y
                r = k;
                modmult(Aux[8].limbs, Aux[11].limbs, Aux[8].limbs);    // v
                modmult(Aux[9].limbs, Aux[5].limbs, Aux[9].limbs);     // w
              }
              toConvert = Aux[8].limbs;
            }
            // Convert from Montgomery to standard notation.
            memset(Aux[4].limbs, 0, NumberLength * sizeof(limb)); // Convert power to standard notation.
            Aux[4].limbs[0].x = 1;
            modmult(Aux[4].limbs, toConvert, toConvert);
            UncompressLimbsBigInteger(toConvert, &SqrtDisc);
          }
          // Obtain inverse of square root stored in SqrtDisc (mod prime).
          intToBigInteger(&tmp2, 1);
          BigIntModularDivision(&tmp2, &SqrtDisc, &prime, &squareRoot);
          correctBits = 1;
          CopyBigInt(&Q, &prime);
          // Obtain nbrBitsSquareRoot correct digits of inverse square root.
          while (correctBits < nbrBitsSquareRoot)
          {   // Compute f(x) = invsqrt(x), f_{n+1}(x) = f_n * (3 - x*f_n^2)/2
            correctBits *= 2;
            BigIntMultiply(&Q, &Q, &Q);           // Square Q.
            BigIntMultiply(&squareRoot, &squareRoot, &tmp1);
            BigIntRemainder(&tmp1, &Q, &tmp2);
            BigIntMultiply(&tmp2, &discriminant, &tmp1);
            BigIntRemainder(&tmp1, &Q, &tmp2);
            intToBigInteger(&tmp1, 3);
            BigIntSubt(&tmp1, &tmp2, &tmp2);
            BigIntMultiply(&tmp2, &squareRoot, &tmp1);
            if (tmp1.limbs[0].x & 1)
            {
              BigIntAdd(&tmp1, &Q, &tmp1);
            }
            BigIntDivide2(&tmp1);
            BigIntRemainder(&tmp1, &Q, &squareRoot);
          }
          // Get square root of discriminant from its inverse by multiplying by discriminant.
          if (squareRoot.sign == SIGN_NEGATIVE)
          {
            BigIntAdd(&squareRoot, &Q, &squareRoot);
          }
          BigIntMultiply(&squareRoot, &discriminant, &squareRoot);
          BigIntRemainder(&squareRoot, &Q, &squareRoot);
          // Multiply by square root of discriminant by prime^deltaZeros.
          for (ctr = 0; ctr < deltaZeros; ctr++)
          {
            BigIntMultiply(&squareRoot, &prime, &squareRoot);
          }
        }
        correctBits = expon - deltaZeros;
        BigIntPowerIntExp(&prime, correctBits, &Q);      // Store increment.
        // Compute x = (b + sqrt(discriminant)) / (-2a) and x = (b - sqrt(discriminant)) / (-2a)
        BigIntAdd(&ValB, &squareRoot, &tmp1);
        for (ctr = 0; ctr < bitsAZero; ctr++)
        {
          BigIntRemainder(&tmp1, &prime, &tmp2);
          if (tmp2.nbrLimbs > 1 || tmp2.limbs[0].x != 0)
          {   // Cannot divide by prime, so go out.
            sol1Invalid = 1;
            break;
          }
          BigIntDivide(&tmp1, &prime, &tmp1);
        }
        BigIntMultiply(&tmp1, &ValAOdd, &tmp1);
        BigIntRemainder(&tmp1, &Q, &Solution1[factorIndex]);
        if (Solution1[factorIndex].sign == SIGN_NEGATIVE)
        {
          BigIntAdd(&Solution1[factorIndex], &Q, &Solution1[factorIndex]);
        }
        BigIntSubt(&ValB, &squareRoot, &tmp1);
        for (ctr = 0; ctr < bitsAZero; ctr++)
        {
          BigIntRemainder(&tmp1, &prime, &tmp2);
          if (tmp2.nbrLimbs > 1 || tmp2.limbs[0].x != 0)
          {   // Cannot divide by prime, so go out.
            sol2Invalid = 1;
            break;
          }
          BigIntDivide(&tmp1, &prime, &tmp1);
        }
        BigIntMultiply(&tmp1, &ValAOdd, &tmp1);
        BigIntRemainder(&tmp1, &Q, &Solution2[factorIndex]);
        if (Solution2[factorIndex].sign == SIGN_NEGATIVE)
        {
          BigIntAdd(&Solution2[factorIndex], &Q, &Solution2[factorIndex]);
        }
      }
      if (sol1Invalid && sol2Invalid)
      {     // Both solutions are invalid. Go out.
        return;
      }
      if (sol1Invalid)
      {     // Solution1 is invalid. Overwrite it with Solution2.
        CopyBigInt(&Solution1[factorIndex], &Solution2[factorIndex]);
      }
      else if (sol2Invalid)
      {     // Solution2 is invalid. Overwrite it with Solution1.
        CopyBigInt(&Solution2[factorIndex], &Solution1[factorIndex]);
      }
      BigIntSubt(&Solution2[factorIndex], &Solution1[factorIndex], &Aux[0]);
      if (Aux[0].sign == SIGN_NEGATIVE)
      {     // Solution2 is less than Solution1, so exchange them.
        CopyBigInt(&Aux[0], &Solution1[factorIndex]);
        CopyBigInt(&Solution1[factorIndex], &Solution2[factorIndex]);
        CopyBigInt(&Solution2[factorIndex], &Aux[0]);
      }
    }
    CopyBigInt(&Increment[factorIndex], &Q);
    Exponents[factorIndex] = 0;
    pstFactor++;
  }   // end for
  // Use Chinese remainder theorem to obtain the solutions.
  do
  {
    multint(&Aux[0], &Increment[0], Exponents[0] >> 1);
    if (Exponents[0] & 1)
    {
      BigIntAdd(&Aux[0], &Solution2[0], &Aux[0]);
    }
    else
    {
      BigIntAdd(&Aux[0], &Solution1[0], &Aux[0]);
    }
    CopyBigInt(&currentSolution, &Aux[0]);
    nbrFactors = astFactorsMod[0].multiplicity;
    pstFactor = &astFactorsMod[1];
    UncompressBigInteger(pstFactor->ptrFactor, &prime);
    BigIntPowerIntExp(&prime, pstFactor->multiplicity, &Mult);
    for (T1 = 1; T1<nbrFactors; T1++)
    {
      pstFactor++;
      expon = Exponents[T1];
      multint(&Aux[T1], &Increment[T1], expon >> 1);
      if (expon & 1)
      {
        BigIntAdd(&Aux[T1], &Solution2[T1], &Aux[T1]);
      }
      else
      {
        BigIntAdd(&Aux[T1], &Solution1[T1], &Aux[T1]);
      }
      NumberLength = *pstFactor->ptrFactor;
      UncompressBigInteger(pstFactor->ptrFactor, &prime);
      BigIntPowerIntExp(&prime, pstFactor->multiplicity, &K1);
      CopyBigInt(&prime, &K1);
      for (E = 0; E<T1; E++)
      {
        BigIntSubt(&Aux[T1], &Aux[E], &Q);
        UncompressBigInteger(astFactorsMod[E+1].ptrFactor, &K);
        BigIntPowerIntExp(&K, astFactorsMod[E+1].multiplicity, &L);
        NumberLength = prime.nbrLimbs;
        memcpy(TestNbr, prime.limbs, NumberLength * sizeof(limb));
        TestNbr[NumberLength].x = 0;
        GetMontgomeryParms(NumberLength);
        BigIntModularDivision(&Q, &L, &prime, &Aux[T1]);
      }
      BigIntRemainder(&Aux[T1], &prime, &L);
      CopyBigInt(&Aux[T1], &L);
      // currentSolution <- Aux[T1] * Mult + currentSolution
      BigIntMultiply(&Aux[T1], &Mult, &L);
      BigIntAdd(&currentSolution, &L, &currentSolution);
      BigIntMultiply(&K1, &Mult, &Mult);
    }   /* end for */
    intToBigInteger(&V, 0);
    for (;;)
    {
      // if V >= GcdAll, exit loop.
      BigIntSubt(&V, &GcdAll, &K1);
      if (K1.sign == SIGN_POSITIVE)
      {
        break;
      }
      // The solution is V*ValNn + currentSolution
      BigIntMultiply(&V, &ValNn, &K1);
      BigIntAdd(&K1, &currentSolution, &K1);
      Solution(&K1);
      addbigint(&V, 1);  // V <- V + 1
    }
    for (T1 = nbrFactors - 1; T1 >= 0; T1--)
    {
      UncompressBigInteger(astFactorsMod[T1+1].ptrFactor, &K);
      BigIntPowerIntExp(&K, astFactorsMod[T1+1].multiplicity, &prime);
      BigIntSubt(&Solution1[T1], &Solution2[T1], &K1);
      if (K1.nbrLimbs == 1 && K1.limbs[0].x == 0)
      {     // Solution1[T1] == Solution2[T1]
        Exponents[T1] += 2;
      }
      else
      {     // Solution1[T1] != Solution2[T1]
        Exponents[T1]++;
      }
      multadd(&L, Exponents[T1], &Increment[T1], 0);   // L <- Exponents[T1] * Increment[T1]
      multadd(&K1, 2, &prime, 0);                      // K1 <- 2 * prime
      BigIntSubt(&L, &K1, &L);
      if (L.sign == SIGN_NEGATIVE)
      {
        break;
      }
      Exponents[T1] = 0;
    }   /* end for */
  } while (T1 >= 0);
}

void textErrorQuadMod(char *pOutput, enum eExprErr rc)
{
  switch (rc)
  {
  case EXPR_MODULUS_MUST_BE_NONNEGATIVE:
    strcpy(pOutput, lang ? "No debe ser negativo" :
      "Must not be negative");
    break;
  default:
    textError(pOutput, rc);
  }
}

void quadmodText(char *quadrText, char *linearText, char *constText, char *modText, int groupLength)
{
  int u;
  char *ptrBeginSol;
  enum eExprErr rc;
  ptrOutput = output;
  strcpy(ptrOutput, "2<p>");
  ptrOutput = output + strlen(output);
  rc = ComputeExpression(quadrText, 1, &ValA);
  if (rc != EXPR_OK)
  {
    strcpy(ptrOutput, lang ? "Coeficiente cuadrático: ": "Quadratic coefficient: ");
    ptrOutput = output + strlen(output);
    textErrorQuadMod(ptrOutput, rc);
    ptrOutput = output + strlen(output);
    strcpy(ptrOutput, "</p>");
    ptrOutput = output + strlen(output);
  }
  else
  {
    rc = ComputeExpression(linearText, 1, &ValB);
    if (rc != EXPR_OK)
    {
      strcpy(ptrOutput, lang ? "Coeficiente lineal: " : "Linear coefficient: ");
      ptrOutput = output + strlen(output);
      textErrorQuadMod(ptrOutput, rc);
      ptrOutput = output + strlen(output);
      strcpy(ptrOutput, "</p>");
      ptrOutput = output + strlen(output);
    }
    else
    {
      rc = ComputeExpression(constText, 1, &ValC);
      if (rc != EXPR_OK)
      {
        strcpy(ptrOutput, lang ? "Término independiente: " : "Constant coefficient: ");
        ptrOutput = output + strlen(output);
        textErrorQuadMod(ptrOutput, rc);
        ptrOutput = output + strlen(output);
        strcpy(ptrOutput, "</p>");
        ptrOutput = output + strlen(output);
      }
      else
      {
        rc = ComputeExpression(modText, 1, &ValN);
        if (rc == EXPR_OK && ValN.sign == SIGN_NEGATIVE)
        {
          rc = EXPR_MODULUS_MUST_BE_NONNEGATIVE;
        }
        if (rc != EXPR_OK)
        {
          strcpy(ptrOutput, lang ? "Módulo: " : "Modulus: ");
          ptrOutput = output + strlen(output);
          textErrorQuadMod(ptrOutput, rc);
          ptrOutput = output + strlen(output);
          strcpy(ptrOutput, "</p>");
          ptrOutput = output + strlen(output);
        }
        else
        {
          u = Show(&ValA, " x&sup2;", 2);
          u = Show(&ValB, " x", u);
          Show1(&ValC, u);
          strcpy(ptrOutput, " &equiv; 0 (mod ");
          ptrOutput += strlen(ptrOutput);
          BigInteger2Dec(&ValN, ptrOutput, groupLen);
          ptrOutput += strlen(ptrOutput);
          strcpy(ptrOutput, ")</p>");
          ptrOutput += strlen(ptrOutput);
          SolNbr = 0;
          ptrBeginSol = ptrOutput;
          strcpy(ptrOutput, "<ol>");
          ptrOutput += strlen(ptrOutput);
          SolveEquation();
          if (SolNbr == 0)
          {
            ptrOutput = ptrBeginSol;
            strcpy(ptrOutput, lang? "<p>No hay soluciones.</p>": "<p>There are no solutions.</p>");
          }
          else
          {
            strcpy(ptrOutput, "</ol>");
          }
          ptrOutput += strlen(ptrOutput);
        }
      }
    }
  }
  strcpy(ptrOutput, lang ? "<p>" COPYRIGHT_SPANISH "</p>" :
    "<p>" COPYRIGHT_ENGLISH "</p>");
  ptrOutput += strlen(ptrOutput);
  *ptrOutput = 0;   // Add string terminator.
}

#ifdef __EMSCRIPTEN__
void doWork(void)
{
  int flags;
  char *ptrData = inputString;
  char *ptrQuadrCoeff, *ptrLinearCoeff, *ptrConstCoeff, *ptrMod;
  groupLen = 0;
  while (*ptrData != ',')
  {
    groupLen = groupLen * 10 + (*ptrData++ - '0');
  }
  ptrData++;                    // Skip comma.
  flags = *ptrData;
  lang = flags & 1;
  ptrQuadrCoeff = ptrData + 2;  // Skip flags and comma.
  ptrLinearCoeff = ptrQuadrCoeff + strlen(ptrQuadrCoeff) + 1;
  ptrConstCoeff = ptrLinearCoeff + strlen(ptrLinearCoeff) + 1;
  ptrMod = ptrConstCoeff + strlen(ptrConstCoeff) + 1;
  quadmodText(ptrQuadrCoeff, ptrLinearCoeff, ptrConstCoeff, ptrMod, groupLen);
  databack(output);
}
#endif

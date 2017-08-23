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
static int nbrToFactor[MAX_LEN];
static int SolNbr;
static char textExp[1000];
static BigInteger LastModulus;
static int nbrFactors;
extern struct sFactors astFactorsMod[1000];
extern int factorsMod[10000];
char *ptrOutput;
struct sFactors astFactorsMod[1000];
int factorsMod[10000];

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
      *ptrOutput++ = ' ';
      BigInteger2Dec(num, ptrOutput, groupLen);
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
  int factorIndex, expon, T1, E, Expo;
  BigInteger GcdAll;
  BigInteger ValNn;
  BigInteger z, Mult, currentSolution;
  BigInteger U, SqrtDisc, prime, K1, Den;
  BigInteger discriminant, V, K, Q, R, Sol1, L;
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
  BigIntGcd(&ValA, &ValB, &Aux[0]);
  BigIntGcd(&ValC, &Aux[0], &GcdAll);
  BigIntRemainder(&ValC, &GcdAll, &Aux[0]);
  if (Aux[0].nbrLimbs != 1 || Aux[0].limbs[0].x != 0)
  {  // ValC must be multiple of gcd(ValA, ValB).
     // Otherwise go out because there are no solutions.
    return;
  }
  // Divide all coefficients and modulus by gcd(ValA, ValB).
  BigIntDivide(&ValA, &GcdAll, &Aux[0]);
  CopyBigInt(&ValA, &Aux[0]);
  BigIntDivide(&ValB, &GcdAll, &Aux[0]);
  CopyBigInt(&ValB, &Aux[0]);
  BigIntDivide(&ValC, &GcdAll, &Aux[0]);
  CopyBigInt(&ValC, &Aux[0]);
  BigIntDivide(&ValN, &GcdAll, &Aux[0]);
  CopyBigInt(&ValN, &Aux[0]);
  CopyBigInt(&ValNn, &Aux[0]);
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
    factor(&ValN, nbrToFactor, factorsMod, astFactorsMod, NULL);
  }
  intToBigInteger(&Q, 0);
  nbrFactors = astFactorsMod[0].multiplicity;
  pstFactor = &astFactorsMod[1];
  for (factorIndex = 0; factorIndex<nbrFactors; factorIndex++)
  {
    NumberLength = *pstFactor->ptrFactor;
    UncompressBigInteger(pstFactor->ptrFactor, &prime);
    expon = pstFactor->multiplicity;
    BigIntRemainder(&ValA, &prime, &V);
    if (V.nbrLimbs == 1 && V.limbs[0].x == 0)
    {     // ValA multiple of prime, means linear equation mod prime.
      if (ValB.nbrLimbs == 1 && ValB.limbs[0].x == 0 &&
        !(ValC.nbrLimbs == 1 && ValC.limbs[0].x == 0))
      {
        return;  // There are no solutions: ValB=0 and ValC!=0
      }
      BigIntPowerIntExp(&prime, expon, &V);
      BigIntRemainder(&ValA, &V, &L);
      if (L.nbrLimbs == 1 && L.limbs[0].x == 0)
      {    // Linear equation mod prime^exp.
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
          BigIntModularDivision(&L, &V, &Q, &Solution1[factorIndex]);
        }
      }
      else
      {
        BigIntGcd(&ValA, &V, &Q);
        if (prime.nbrLimbs == 1 && prime.limbs[0].x == 2)
        {      // If prime is 2.
          BigIntModularDivisionPower2(&ValC, &ValB, &Q, &Sol1);
        }
        else
        {
          BigIntModularDivision(&ValC, &ValB, &Q, &Sol1);
        }
        T1 = 0;
        CopyBigInt(&V, &prime);
        for (;;)
        {
          BigIntRemainder(&ValA, &V, &L);
          if (L.nbrLimbs > 1 || L.limbs[0].x != 0)
          {
            break;
          }
          BigIntMultiply(&V, &prime, &V);
          T1++;
        }
        do
        {
          // L <- (ValA * (Sol1)^2 + ValB * Sol1 + ValC) / Q
          BigIntMultiply(&ValA, &Sol1, &K1);
          BigIntAdd(&K1, &ValB, &L);
          BigIntMultiply(&L, &Sol1, &L);
          BigIntAdd(&L, &ValC, &L);
          BigIntDivide(&L, &Q, &K);
          CopyBigInt(&L, &K);
          // K1 <- 2 * ValA * Sol1 + ValB
          BigIntAdd(&K1, &K1, &K1);
          BigIntAdd(&K1, &ValB, &K1);
          // K <- gcd(L, K1)
          BigIntGcd(&L, &K1, &K);
          // L <- L/K
          BigIntDivide(&L, &K, &R);
          CopyBigInt(&L, &R);
          // K1 <- K1/K
          BigIntDivide(&K1, &K, &R);
          CopyBigInt(&K1, &R);
          // R <- L/K1 (mod Q)
          if (prime.nbrLimbs == 1 && prime.limbs[0].x == 2)
          {      // If prime is 2.
            BigIntModularDivisionPower2(&L, &K1, &Q, &R);
          }
          else
          {
            BigIntModularDivision(&L, &K1, &Q, &R);
          }
          // Sol1 <- Sol1 + R*Q
          BigIntMultiply(&R, &Q, &K);
          BigIntAdd(&Sol1, &K, &Sol1);
          T1 = T1 + T1;
          BigIntMultiply(&Q, &Q, &Q);
        } while (T1 < expon);
        BigIntPowerIntExp(&prime, expon, &Q);
        BigIntRemainder(&Sol1, &Q, &SqrtDisc);
        intToBigInteger(&Mult, 1);
        Solution1[factorIndex] = SqrtDisc;
      }
      Solution2[factorIndex] = Solution1[factorIndex];
    }
    else
    {                   /* If quadratic equation mod p */
      Expo = 0;
      if (prime.nbrLimbs == 1 && prime.limbs[0].x == 2)
      {         /* Prime p is 2 */
        BigInteger coefX2, coefX1, coefX0;
        BigInteger increm;
        BigInteger deltaincrem;
        CopyBigInt(&coefX2, &ValA);
        CopyBigInt(&coefX1, &ValB);
        CopyBigInt(&coefX0, &ValC);
        intToBigInteger(&increm, 0);
        intToBigInteger(&deltaincrem, 1);
        Expo = 0;
        while (Expo < expon - 1 && (coefX1.limbs[0].x & 1) == 0)
        {
          if ((coefX0.limbs[0].x & 1) == 0)
          {                   // constant term is even
            if ((coefX0.limbs[0].x & 2) == 2)
            {                 // constant term = 2 (mod 4)
              return;         // No solutions
            }
            multadd(&coefX1, 2, &coefX1, 0);   // coefX1 <- 2*coefX1
            multadd(&coefX0, 4, &coefX0, 0);   // coefX0 <- 4*coefX0
            Expo += 2;
            multadd(&deltaincrem, 2, &deltaincrem, 0);
          }
          else {              // constant term is odd
            BigIntAdd(&coefX0, &coefX1, &coefX0);
            BigIntAdd(&coefX0, &coefX2, &coefX0);
            BigIntAdd(&coefX1, &coefX2, &coefX1);
            BigIntAdd(&coefX1, &coefX2, &coefX1);
            BigIntAdd(&increm, &deltaincrem, &increm);
          }
        }
        if (Expo == expon - 1)
        {                     // Equation modulo 2 
          if ((coefX1.limbs[0].x & 1) == 0)
          {                   // Linear coefficient is even.
            BigIntPowerOf2(&Q, (expon + 1) / 2);
            if (coefX0.limbs[0].x & 1)
            {     // Constant coefficient is odd.
              BigIntPowerOf2(&Solution1[factorIndex], Expo / 2);
              BigIntAdd(&Solution1[factorIndex], &increm, &Solution1[factorIndex]);
            }
            else
            {     // Constant coefficient is even.
              CopyBigInt(&Solution1[factorIndex], &increm);
            }
          }
          else
          {
            if (coefX0.limbs[0].x & 1)
            {                 // All three coefficients odd
              return;         // No solution
            }
            BigIntPowerOf2(&Q, Expo / 2);
            CopyBigInt(&Solution1[factorIndex], &increm);
          }
          CopyBigInt(&Solution2[factorIndex], &Solution1[factorIndex]);
        }
        else
        {
          if (Expo == expon)
          {            // equation modulo 1.
            BigIntPowerOf2(&Q, Expo / 2);
            CopyBigInt(&Solution1[factorIndex], &increm);
            CopyBigInt(&Solution2[factorIndex], &increm);
          }
          else
          {            // coefX1 odd -> 1 solution odd, 1 solution even
            // Compute discriminant <- (coefX1)^2 - 4 * coefX0 * coefX2.
            BigIntMultiply(&coefX1, &coefX1, &Aux[0]);
            BigIntMultiply(&coefX0, &coefX2, &Aux[1]);
            multadd(&Aux[1], 4, &Aux[1], 0);
            BigIntSubt(&Aux[0], &Aux[1], &discriminant);
            if (coefX0.limbs[0].x & 1)
            {
              return;          // The three numbers cannot be all odd.
            }
            intToBigInteger(&Sol1, 1);
            T1 = 1;
            intToBigInteger(&Q, 2);
            if (expon != 1)
            {
              do
              {
                // L <- (ValA * (Sol1)^2 + coefX1 * Sol1 + CoefX0) / Q
                BigIntMultiply(&ValA, &Sol1, &L);
                BigIntAdd(&L, &coefX1, &L);
                BigIntMultiply(&L, &Sol1, &L);
                BigIntAdd(&L, &coefX0, &L);
                BigIntDivide(&L, &Q, &Aux[0]);
                BigIntNegate(&Aux[0], &L);
                // R <- -L/(2*Sol1*coefX2 + coefX1) (mod Q)
                BigIntAdd(&Sol1, &Sol1, &Aux[0]);
                BigIntMultiply(&Aux[0], &coefX2, &Aux[0]);
                BigIntAdd(&Aux[0], &coefX1, &Aux[0]);
                BigIntModularDivisionPower2(&L, &Aux[0], &Q, &R);
                BigIntMultiply(&Q, &R, &Aux[0]);
                BigIntAdd(&Sol1, &Aux[0], &Sol1);
                T1 = T1 + T1;
                BigIntMultiply(&Q, &Q, &Q);
              } while (T1 < expon - Expo / 2);
            }
            BigIntPowerOf2(&Q, expon - Expo);
            // Solution1[factorIndex] <- Sol1 % Q * 2**(Expo / 2) + increm
            BigIntPowerOf2(&Aux[0], Expo/2);
            BigIntRemainder(&Sol1, &Q, &Aux[1]);
            BigIntMultiply(&Aux[0], &Aux[1], &Aux[1]);
            BigIntAdd(&Aux[1], &increm, &Solution1[factorIndex]);
            // Solution2[factorIndex] <- -(coefX1 / coefX2 + Sol1) % Q * 2**(Expo / 2) + increm
            BigIntModularDivisionPower2(&coefX1, &coefX2, &Q, &Aux[1]);
            BigIntAdd(&Aux[1], &Sol1, &Aux[1]);
            BigIntRemainder(&Aux[1], &Q, &Aux[2]);
            BigIntMultiply(&Aux[0], &Aux[2], &Aux[1]);
            BigIntSubt(&increm, &Aux[2], &Solution2[factorIndex]);
            if (Solution2[factorIndex].sign == SIGN_NEGATIVE)
            {
              BigIntAdd(&Q, &Solution2[factorIndex], &Solution2[factorIndex]);
            }
            BigIntSubt(&Solution1[factorIndex], &Solution2[factorIndex], &Aux[1]);
            if (Aux[1].sign == SIGN_POSITIVE)
            {   // Make Solution2 the greatest one.
              CopyBigInt(&Sol1, &Solution1[factorIndex]);
              CopyBigInt(&Solution1[factorIndex], &Solution2[factorIndex]);
              CopyBigInt(&Solution2[factorIndex], &Sol1);
            }
            BigIntMultiply(&Q, &Aux[0], &Q);
          }
        }
      }
      else
      {                        // Prime is not 2
        // discriminant <- (ValB)^2 - 4 * ValA * ValC
        BigIntMultiply(&ValB, &ValB, &Aux[0]);
        BigIntMultiply(&ValA, &ValC, &Aux[1]);
        multadd(&Aux[1], 4, &Aux[1], 0);
        BigIntSubt(&Aux[0], &Aux[1], &discriminant);
        for (; Expo < expon; Expo++)
        {   // Find exponent
          BigIntRemainder(&discriminant, &prime, &Aux[0]);
          if (Aux[0].nbrLimbs == 1 && Aux[0].limbs[0].x == 0)
          {    // discriminant is multiple of prime.
            BigIntDivide(&discriminant, &prime, &Aux[0]);
            CopyBigInt(&discriminant, &Aux[0]);
          }
          else
          {
            break;
          }
        }
        if (Expo == expon)
        {
          BigIntPowerIntExp(&prime, (Expo + 1) / 2, &Q);
          BigIntNegate(&ValB, &Aux[2]);
          BigIntAdd(&ValA, &ValA, &Aux[3]);
          NumberLength = Q.nbrLimbs;
          memcpy(TestNbr, Q.limbs, NumberLength * sizeof(limb));
          TestNbr[NumberLength].x = 0;
          GetMontgomeryParms(NumberLength);
          BigIntModularDivision(&Aux[2], &Aux[3], &Q, &Solution1[factorIndex]);
          CopyBigInt(&Solution2[factorIndex], &Solution1[factorIndex]);
        }
        else
        {
          if ((Expo & 1) == 1)
          {
            return;         // No solutions if odd exponent.
          }
          BigIntPowerIntExp(&prime, Expo / 2, &Mult);
          if (BigIntJacobiSymbol(&discriminant, &prime) != 1)
          {
            return;         // Not a quadratic residue, so go out.
          }
          NumberLength = prime.nbrLimbs;
          memcpy(TestNbr, prime.limbs, NumberLength * sizeof(limb));
          TestNbr[NumberLength].x = 0;
          GetMontgomeryParms(NumberLength);
          CopyBigInt(&Q, &prime);
                            // Obtain discriminant mod prime.
          BigIntRemainder(&discriminant, &prime, &Aux[3]);
          if (Aux[3].sign == SIGN_NEGATIVE)
          {
            BigIntAdd(&Aux[3], &prime, &Aux[3]);
          }
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
              int x, e, r, k, ctr;
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
          CopyBigInt(&Sol1, &SqrtDisc);
          T1 = 1;
          CopyBigInt(&Q, &prime);
          if (expon != 1)
          {            // Use Hensel lifting for solutions mod prime^n.
            do
            {
              NumberLength = Q.nbrLimbs;
              memcpy(TestNbr, Q.limbs, NumberLength * sizeof(limb));
              TestNbr[NumberLength].x = 0;
              GetMontgomeryParms(NumberLength);
              // L <- ((Sol1)^2 - discriminant) / Q
              BigIntMultiply(&Sol1, &Sol1, &K1);
              BigIntSubt(&K1, &discriminant, &K1);
              BigIntDivide(&K1, &Q, &L);
              // K <- gcd(L, 2*Sol1)
              BigIntAdd(&Sol1, &Sol1, &K1);
              BigIntGcd(&K1, &L, &K);
              // R <- (-L/K) / (2*Sol1/K) (mod Q)
              BigIntNegate(&L, &Aux[4]);
              BigIntDivide(&Aux[4], &K, &Aux[5]);
              BigIntAdd(&Sol1, &Sol1, &Aux[6]);
              BigIntDivide(&Aux[6], &K, &Aux[7]);
              BigIntModularDivision(&Aux[5], &Aux[7], &Q, &R);
              // Sol1 <- R*Q + Sol1
              BigIntMultiply(&R, &Q, &Aux[4]);
              BigIntAdd(&Sol1, &Aux[4], &Sol1);
              T1 = T1 + T1;
              // Q <- Q*Q
              BigIntMultiply(&Q, &Q, &Q);
            } while (T1 < expon - Expo / 2);
            BigIntPowerIntExp(&prime, expon - Expo / 2, &Q);
          }
          NumberLength = Q.nbrLimbs;
          memcpy(TestNbr, Q.limbs, NumberLength * sizeof(limb));
          TestNbr[NumberLength].x = 0;
          GetMontgomeryParms(NumberLength);
          // SqrtDisc <- Sol1 % Q * Mult
          BigIntRemainder(&Sol1, &Q, &K1);
          BigIntMultiply(&K1, &Mult, &SqrtDisc);
          // Den <- ValA * 2
          BigIntAdd(&ValA, &ValA, &Den);
          BigIntNegate(&ValB, &K1);
          BigIntSubt(&K1, &SqrtDisc, &Aux[4]);
          BigIntModularDivision(&Aux[4], &Den, &Q, &Solution1[factorIndex]);
          BigIntAdd(&K1, &SqrtDisc, &Aux[4]);
          BigIntModularDivision(&Aux[4], &Den, &Q, &Solution2[factorIndex]);
          BigIntSubt(&Solution1[factorIndex], &Solution2[factorIndex], &Aux[4]);
          if (Aux[4].sign == SIGN_POSITIVE)
          {   // Exchange solutions so Solution1 is always less than Solution2.
            CopyBigInt(&Sol1, &Solution1[factorIndex]);
            CopyBigInt(&Solution1[factorIndex], &Solution2[factorIndex]);
            CopyBigInt(&Solution2[factorIndex], &Sol1);
          }
        }
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

void textErrorQuadMod(char *ptrOutput, enum eExprErr rc)
{
  char text[150];

  switch (rc)
  {
  case EXPR_MODULUS_MUST_BE_NONNEGATIVE:
    strcpy(text, lang ? "No debe ser negativo" :
      "Must not be negative");
    break;
  default:
    textError(text, rc);
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

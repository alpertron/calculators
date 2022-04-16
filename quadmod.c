//
// This file is part of Alpertron Calculators.
//
// Copyright 2017-2021 Dario Alejandro Alpern
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

static BigInteger sqrRoot;
static BigInteger Quadr;
static BigInteger Linear;
static BigInteger Const;
static BigInteger tmp1;
static BigInteger tmp2;
static BigInteger ValAOdd;
static BigInteger ValBOdd;
static BigInteger ValCOdd;
static BigInteger GcdAll;
static BigInteger ValNn;
static BigInteger z;
static BigInteger Mult;
static BigInteger currentSolution;
static BigInteger SqrtDisc;
static BigInteger prime;
static BigInteger K1;
static BigInteger discriminant;
static BigInteger V;
static BigInteger K;
static BigInteger Q;
static BigInteger L;
static BigInteger Solution1[400];
static BigInteger Solution2[400];
static BigInteger Increment[400];
static int Exponents[400];
static BigInteger Aux[400];
BigInteger ValA;
BigInteger ValB;
BigInteger ValC;
BigInteger ValN;
static int SolNbr;
static BigInteger LastModulus;
static int nbrFactors;
char *ptrOutput;
static bool sol1Invalid;
static bool sol2Invalid;
extern int factorsMod[20000];

static int Show(const BigInteger *num, const char *str, int t)
{
  if (!BigIntIsZero(num))
  {     // num is not zero.
    if (((t & 1) != 0) && (num->sign == SIGN_POSITIVE))
    {
      copyStr(&ptrOutput, " +");
    }
    if (num->sign == SIGN_NEGATIVE)
    {
      copyStr(&ptrOutput, " &minus;");
    }
    if ((num->nbrLimbs != 1) || (num->limbs[0].x != 1))
    {    // num is not 1 or -1.
      *ptrOutput = ' ';
      ptrOutput++;
      Bin2Dec(&ptrOutput, num->limbs, num->nbrLimbs, groupLen);
    }
    copyStr(&ptrOutput, str);
    return t | 1;
  }
  return t;
}

void Show1(const BigInteger *num, int t)
{
  int u = Show(num, "", t);
  if (((u & 1) == 0) || ((num->nbrLimbs == 1) && (num->limbs[0].x == 1)))
  {
    *ptrOutput = ' ';
    ptrOutput++;
    CopyBigInt(&z, num);
    z.sign = SIGN_POSITIVE;
    BigInteger2Dec(&ptrOutput, &z, groupLen);
  }
}

void Solution(const BigInteger *value)
{
  SolNbr++;
  copyStr(&ptrOutput, "<li>x = ");
  BigInteger2Dec(&ptrOutput, value, groupLen);
  copyStr(&ptrOutput, "</li>");
}

static void setNbrLimbs(BigInteger* pBigNbr)
{
  pBigNbr->nbrLimbs = NumberLength;
  pBigNbr->sign = SIGN_POSITIVE;
  while (pBigNbr->nbrLimbs > 1)
  {
    if (pBigNbr->limbs[pBigNbr->nbrLimbs - 1].x != 0)
    {
      break;
    }
    pBigNbr->nbrLimbs--;
  }
}

// Find quadratic solution of Quadr*x^2 + Linear*x + Const = 0 (mod 2^expon)
// when Quadr is even and Linear is odd. In this case there is unique solution.
static void findQuadraticSolution(BigInteger* pSolution, int exponent)
{
  int bytesLen;
  int expon = exponent;
  int bitMask = 1;
  limb *ptrSolution = pSolution->limbs;
  BigIntPowerOf2(&Q, expon);
  bytesLen = Q.nbrLimbs * (int)sizeof(limb);
  (void)memset(pSolution->limbs, 0, bytesLen);
  while (expon > 0)
  {
    expon--;
    BigIntPowerOf2(&K1, expon);
    addbigint(&K1, -1);
    if ((Const.limbs[0].x & 1) != 0)
    {        // Const is odd.
      ptrSolution->x |= bitMask;
  // Const <- Quadr/2 + floor(Linear/2) + floor(Const/2) + 1
      if (Const.sign == SIGN_NEGATIVE)
      {
        addbigint(&Const, -1);
      }
      BigIntDivideBy2(&Const);          // floor(Const/2)
      addbigint(&Const, 1);             // floor(Const/2) + 1
      CopyBigInt(&tmp1, &Linear);
      if (tmp1.sign == SIGN_NEGATIVE)
      {
        addbigint(&tmp1, -1);
      }
      BigIntDivideBy2(&tmp1);           // floor(Linear/2)
      BigIntAdd(&Const, &tmp1, &Const);
      CopyBigInt(&tmp1, &Quadr);
      BigIntDivideBy2(&tmp1);           // Quadr/2
      BigIntAdd(&Const, &tmp1, &Const);

      // Linear <- 2*Quadr + Linear and Quadr <- 2*Quadr.
      BigIntMultiplyBy2(&Quadr);          // Quadr*2
      BigIntAdd(&Linear, &Quadr, &Linear);
      BigIntAnd(&Linear, &K1, &Linear);   // Reduce mod 2^expon
    }
    else
    {        // Const is even.
      BigIntDivideBy2(&Const);          // Const/2
      BigIntMultiplyBy2(&Quadr);        // Quadr*2
    }
    BigIntAnd(&Const, &K1, &Const);     // Reduce mod 2^expon
    BigIntAnd(&Quadr, &K1, &Quadr);     // Reduce mod 2^expon
    bitMask *= 2;
    if (bitMask < 0)
    {
      bitMask = 1;
      ptrSolution++;
    }
  }
  NumberLength = Q.nbrLimbs;
  setNbrLimbs(pSolution);
}

// Solve Ax^2 + Bx + C = 0 in integers.
static void SolveIntegerEquation(void)
{
  if (BigIntIsZero(&ValA))
  {      // Not Quadratic Equation
    if (BigIntIsZero(&ValB))
    {    // Constant Equation
      if (BigIntIsZero(&ValC))
      {  // 0 = 0
        copyStr(&ptrOutput, "<p>The equation is satisfied by any integer <var>x</var>.</p>");
      }
      else
      {
        return;
      }
    }
    else
    {        // Linear Equation: ValB * x + ValC = 0.
      (void)BigIntRemainder(&ValC, &ValB, &Aux[0]);
      if (BigIntIsZero(&Aux[0]))
      {      // ValC is multiple of ValB: solution is -ValC / ValB.
        (void)BigIntDivide(&ValC, &ValB, &Aux[0]);
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
    (void)BigIntMultiply(&ValB, &ValB, &Aux[0]);
    (void)BigIntMultiply(&ValA, &ValC, &Aux[1]);
    multint(&Aux[1], &Aux[1], 4);
    BigIntSubt(&Aux[0], &Aux[1], &discriminant);
    if (discriminant.sign == SIGN_NEGATIVE)
    {        // No integer solutions.
      return;
    }
    // Set V to square root of discriminant.
    squareRoot(discriminant.limbs, V.limbs, discriminant.nbrLimbs, &V.nbrLimbs);
    V.sign = SIGN_POSITIVE;
    (void)BigIntMultiply(&V, &V, &Aux[0]);
    BigIntSubt(&Aux[0], &discriminant, &Aux[0]);
    if (!BigIntIsZero(&Aux[0]))
    {  // discriminant has no integer square root.
      return;
    }
    multint(&Aux[0], &ValA, 2);    // Denominator
    BigIntSubt(&V, &ValB, &Aux[1]);
    (void)BigIntRemainder(&Aux[1], &Aux[0], &Aux[2]);
    if (BigIntIsZero(&Aux[2]))
    {      // (V-ValB)/(2*ValA) is integer: it is a solution.
      (void)BigIntDivide(&Aux[1], &Aux[0], &Aux[2]);
      Solution(&Aux[2]);
    }
    BigIntNegate(&V, &V);
    BigIntSubt(&V, &ValB, &Aux[1]);
    (void)BigIntRemainder(&Aux[1], &Aux[0], &Aux[2]);
    if (BigIntIsZero(&Aux[2]))
    {      // (-V-ValB)/(2*ValA) is integer: it is a solution.
      (void)BigIntDivide(&Aux[1], &Aux[0], &Aux[2]);
      Solution(&Aux[2]);
    }
  }
  SolNbr = 1;
}

// Use Chinese remainder theorem to obtain the solutions.
static void PerformChineseRemainderTheorem(void)
{
  int T1;
  int expon;
  do
  {
    const struct sFactors* pstFactor;
    multint(&Aux[0], &Increment[0], Exponents[0] / 2);
    if ((Exponents[0] & 1) != 0)
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
    IntArray2BigInteger(pstFactor->ptrFactor, &prime);
    (void)BigIntPowerIntExp(&prime, pstFactor->multiplicity, &Mult);
    for (T1 = 1; T1 < nbrFactors; T1++)
    {
      pstFactor++;
      expon = Exponents[T1];
      multint(&Aux[T1], &Increment[T1], expon / 2);
      if ((expon & 1) != 0)
      {
        BigIntAdd(&Aux[T1], &Solution2[T1], &Aux[T1]);
      }
      else
      {
        BigIntAdd(&Aux[T1], &Solution1[T1], &Aux[T1]);
      }
      NumberLength = *pstFactor->ptrFactor;
      IntArray2BigInteger(pstFactor->ptrFactor, &prime);
      (void)BigIntPowerIntExp(&prime, pstFactor->multiplicity, &K1);
      CopyBigInt(&prime, &K1);
      for (int E = 0; E < T1; E++)
      {
        int NumberLengthBytes;
        BigIntSubt(&Aux[T1], &Aux[E], &Q);
        IntArray2BigInteger(astFactorsMod[E + 1].ptrFactor, &K);
        (void)BigIntPowerIntExp(&K, astFactorsMod[E + 1].multiplicity, &L);
        NumberLength = prime.nbrLimbs;
        NumberLengthBytes = NumberLength * (int)sizeof(limb);
        (void)memcpy(TestNbr, prime.limbs, NumberLengthBytes);
        TestNbr[NumberLength].x = 0;
        GetMontgomeryParms(NumberLength);
        BigIntModularDivision(&Q, &L, &prime, &Aux[T1]);
      }
      (void)BigIntRemainder(&Aux[T1], &prime, &L);
      CopyBigInt(&Aux[T1], &L);
      // currentSolution <- Aux[T1] * Mult + currentSolution
      (void)BigIntMultiply(&Aux[T1], &Mult, &L);
      BigIntAdd(&currentSolution, &L, &currentSolution);
      (void)BigIntMultiply(&K1, &Mult, &Mult);
    }   /* end for */
    intToBigInteger(&V, 0);
    BigIntSubt(&V, &GcdAll, &K1);
    // Perform loop while V < GcdAll.
    while (K1.sign == SIGN_NEGATIVE)
    {      
      // The solution is V*ValNn + currentSolution
      (void)BigIntMultiply(&V, &ValNn, &K1);
      BigIntAdd(&K1, &currentSolution, &K1);
      Solution(&K1);
      addbigint(&V, 1);  // V <- V + 1
      BigIntSubt(&V, &GcdAll, &K1);
    }
    for (T1 = nbrFactors - 1; T1 >= 0; T1--)
    {
      IntArray2BigInteger(astFactorsMod[T1 + 1].ptrFactor, &K);
      (void)BigIntPowerIntExp(&K, astFactorsMod[T1 + 1].multiplicity, &prime);
      BigIntSubt(&Solution1[T1], &Solution2[T1], &K1);
      if ((K1.nbrLimbs == 1) && (K1.limbs[0].x == 0))
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

// Solve Bx + C = 0 (mod N).
static void SolveModularLinearEquation(void)
{
  int NumberLengthBytes;
  int powerOf2;
  int solutionNbr = 0;
  int* ptrFactorsMod = factorsMod;
  struct sFactors* pstFactor = &astFactorsMod[1];
  BigInteger* ptrSolution1 = Solution1;
  BigInteger* ptrSolution2 = Solution2;
  BigIntGcd(&ValB, &ValN, &Aux[0]);
  if ((Aux[0].nbrLimbs != 1) || (Aux[0].limbs[0].x != 1))
  {         // ValB and ValN are not coprime. Go out.
    return;
  }
  // Calculate z <- -ValC / ValB (mod ValN)
  // Modular division routines used work for power of 2 or odd numbers.
  // This requires to compute the quotient in two steps.
  // N = r*2^k (r = odd)
  DivideBigNbrByMaxPowerOf2(&powerOf2, ValN.limbs, &ValN.nbrLimbs);
  NumberLength = ValN.nbrLimbs;
  NumberLengthBytes = NumberLength * (int)sizeof(limb);
  if ((ValN.nbrLimbs != 1) || (ValN.limbs[0].x != 1))
  {       // ValN is not 1.
    CopyBigInt(&Increment[solutionNbr], &ValN);
    Exponents[solutionNbr] = 1;
    (void)memcpy(TestNbr, ValN.limbs, NumberLengthBytes);
    TestNbr[NumberLength].x = 0;
    // Perform division using odd modulus r.
    GetMontgomeryParms(NumberLength);
    // ptrSolution1 <- ValC / |ValB|
    BigIntModularDivision(&ValC, &ValB, &ValN, ptrSolution1);
    // ptrSolution1 <- -ValC / ValB
    if (!BigIntIsZero(ptrSolution1))
    {
      BigIntSubt(&ValN, ptrSolution1, ptrSolution1);
    }
    CopyBigInt(ptrSolution2, ptrSolution1);
    BigInteger2IntArray(ptrFactorsMod, &ValN);
    pstFactor->ptrFactor = ptrFactorsMod;
    pstFactor->multiplicity = 1;
    pstFactor++;
    ptrFactorsMod += *ptrFactorsMod;
    ptrFactorsMod++;
    ptrSolution1++;
    ptrSolution2++;
    solutionNbr++;
  }
  // Perform division using power of 2.
  if (powerOf2 > 0)
  {
    BigIntPowerOf2(ptrSolution1, powerOf2);
    CopyBigInt(&Increment[solutionNbr], ptrSolution1);
    Exponents[solutionNbr] = 1;
    BigInteger2IntArray(ptrFactorsMod, ptrSolution1);
    pstFactor->ptrFactor = ptrFactorsMod;
    pstFactor->multiplicity = 1;
    GetMontgomeryParmsPowerOf2(powerOf2);
    // Use ValA (which is zero for linear equations) as a temporary area.
    // ptrSolution1 <- 1 / |ValB|
    ComputeInversePower2(ValB.limbs, ptrSolution1->limbs, ValA.limbs);
    // ptrSolution1 <- |ValC| / |ValB|
    modmult(ptrSolution1->limbs, ValC.limbs, ptrSolution1->limbs);
    NumberLengthBytes = NumberLength * (int)sizeof(int);
    // ptrSolution1 <- -ValC / ValB
    if (ValB.sign == ValC.sign)
    {
      (void)memset(&ValA.limbs, 0, NumberLengthBytes);
      SubtractBigNbr(ValA.limbs, ptrSolution1->limbs, ptrSolution1->limbs, NumberLength);
    }
    // Discard bits outside number in most significant limb.
    ptrSolution1->limbs[NumberLength - 1].x &= (1 << (powerOf2 % BITS_PER_GROUP)) - 1;
    ptrSolution1->nbrLimbs = NumberLength;
    ptrSolution1->sign = SIGN_POSITIVE;
    CopyBigInt(ptrSolution2, ptrSolution1);
    solutionNbr++;
  }
  astFactorsMod[0].multiplicity = solutionNbr;
  PerformChineseRemainderTheorem();
}

// To compute the square root, compute the inverse of sqrt,
// so only multiplications are used.
// f(x) = invsqrt(x), f_{n+1}(x) = f_n * (3 - x*f_n^2)/2
static void ComputeSquareRootModPowerOf2(int expon, int bitsCZero)
{
  int lenBytes;
  int correctBits;
  int nbrLimbs;
  // First approximation to inverse of square root.
  sqrRoot.limbs[0].x = (((ValCOdd.limbs[0].x & 15) == 1) ? 1 : 3);
  correctBits = 2;
  nbrLimbs = 1;
  while (correctBits < expon)
  {   // Compute f(x) = invsqrt(x), f_{n+1}(x) = f_n * (3 - x*f_n^2)/2
    correctBits *= 2;
    nbrLimbs = (correctBits / BITS_PER_GROUP) + 1;
    MultBigNbr(sqrRoot.limbs, sqrRoot.limbs, tmp2.limbs, nbrLimbs);
    MultBigNbr(tmp2.limbs, ValCOdd.limbs, tmp2.limbs, nbrLimbs);
    ChSignBigNbr(tmp2.limbs, nbrLimbs);
    lenBytes = nbrLimbs * (int)sizeof(limb);
    (void)memset(tmp1.limbs, 0, lenBytes);
    tmp1.limbs[0].x = 3;
    AddBigNbr(tmp1.limbs, tmp2.limbs, tmp2.limbs, nbrLimbs);
    MultBigNbr(tmp2.limbs, sqrRoot.limbs, tmp1.limbs, nbrLimbs);
    (void)memcpy(sqrRoot.limbs, tmp1.limbs, lenBytes);
    DivBigNbrByInt(tmp1.limbs, 2, sqrRoot.limbs, nbrLimbs);
  }
  // Get square root of ValCOdd from its inverse by multiplying by ValCOdd.
  MultBigNbr(ValCOdd.limbs, sqrRoot.limbs, tmp1.limbs, nbrLimbs);
  lenBytes = nbrLimbs * (int)sizeof(limb);
  (void)memcpy(sqrRoot.limbs, tmp1.limbs, lenBytes);
  setNbrLimbs(&sqrRoot);
  for (int ctr = 0; ctr < (bitsCZero / 2); ctr++)
  {
    BigIntMultiplyBy2(&sqrRoot);
  }
}

// Solve Ax^2 + Bx + C = 0 (mod 2^expon).
static bool SolveQuadraticEqModPowerOf2(int exponent, int factorIndex)
{
  int expon = exponent;
  int bitsAZero;
  int bitsBZero;
  int bitsCZero;
  // ax^2 + bx + c = 0 (mod 2^expon)
  // This follows the paper Complete solving the quadratic equation mod 2^n
  // of Dehnavi, Shamsabad and Rishakani.
  // Get odd part of A, B and C and number of bits to zero.
  CopyBigInt(&ValAOdd, &ValA);
  DivideBigNbrByMaxPowerOf2(&bitsAZero, ValAOdd.limbs, &ValAOdd.nbrLimbs);
  CopyBigInt(&ValBOdd, &ValB);
  DivideBigNbrByMaxPowerOf2(&bitsBZero, ValBOdd.limbs, &ValBOdd.nbrLimbs);
  CopyBigInt(&ValCOdd, &ValC);
  DivideBigNbrByMaxPowerOf2(&bitsCZero, ValCOdd.limbs, &ValCOdd.nbrLimbs);
  if ((bitsAZero > 0) && (bitsBZero > 0) && (bitsCZero > 0))
  {
    int minExpon = bitsAZero;
    if (minExpon < bitsBZero)
    {
      minExpon = bitsBZero;
    }
    if (minExpon < bitsCZero)
    {
      minExpon = bitsCZero;
    }
    bitsAZero -= minExpon;
    bitsBZero -= minExpon;
    bitsCZero -= minExpon;
    expon -= minExpon;
  }
  if (((bitsAZero == 0) && (bitsBZero == 0) && (bitsCZero == 0)) ||
    ((bitsAZero > 0) && (bitsBZero > 0) && (bitsCZero == 0)))
  {
    return false;   // No solutions, so go out.
  }
  if ((bitsAZero == 0) && (bitsBZero > 0))
  {           // The solution in this case requires square root.
    // compute s = ((b/2)^2 - a*c)/a^2, q = odd part of s,
    // r = maximum exponent of power of 2 that divides s.
    CopyBigInt(&tmp1, &ValB);
    BigIntDivideBy2(&tmp1);
    (void)BigIntMultiply(&tmp1, &tmp1, &tmp1);  // (b/2)^2
    (void)BigIntMultiply(&ValA, &ValC, &tmp2);  // a*c
    BigIntSubt(&tmp1, &tmp2, &tmp1);      // (b/2)^2 - a*c
    BigIntPowerOf2(&K1, expon);
    addbigint(&K1, -1);
    BigIntAnd(&tmp1, &K1, &ValCOdd);      // (b/2) - a*c mod 2^n
    NumberLength = K1.nbrLimbs;
    if (NumberLength > ValAOdd.nbrLimbs)
    {
      int lenBytes = (NumberLength - ValAOdd.nbrLimbs) * (int)sizeof(int);
      (void)memset(&ValAOdd.limbs[ValAOdd.nbrLimbs], 0, lenBytes);
    }
    if (NumberLength > tmp2.nbrLimbs)
    {
      int lenBytes = (NumberLength - tmp2.nbrLimbs) * (int)sizeof(int);
      (void)memset(&tmp2.limbs[tmp2.nbrLimbs], 0, lenBytes);
    }
    ComputeInversePower2(ValAOdd.limbs, tmp2.limbs, tmp1.limbs);
    (void)BigIntMultiply(&ValCOdd, &tmp2, &ValCOdd);
    BigIntAnd(&ValCOdd, &K1, &ValCOdd);      // ((b/2) - a*c)/a mod 2^n
    (void)BigIntMultiply(&ValCOdd, &tmp2, &ValCOdd);
    BigIntAnd(&ValCOdd, &K1, &ValCOdd);      // s = ((b/2) - a*c)/a^2 mod 2^n
    if (BigIntIsZero(&ValCOdd))
    {         // s = 0, so its square root is also zero.
      intToBigInteger(&sqrRoot, 0);
      expon -= expon / 2;
    }
    else
    {
      DivideBigNbrByMaxPowerOf2(&bitsCZero, ValCOdd.limbs, &ValCOdd.nbrLimbs);
      // At this moment, bitsCZero = r and ValCOdd = q.
      if (((ValCOdd.limbs[0].x & 7) != 1) || (bitsCZero & 1))
      {
        return false;          // q != 1 or p2(r) == 0, so go out.
      }
      if (expon < 2)
      {                        // Modulus is 2.
        intToBigInteger(&sqrRoot, (bitsCZero > 0) ? 0 : 1);
      }
      else
      {
        // Find square root of ValCOdd.
        ComputeSquareRootModPowerOf2(expon - (bitsCZero / 2) - 1, bitsCZero);
        expon--;
      }
    }
    // x = sqrRoot - b/2a.
    BigIntPowerOf2(&K1, expon);
    addbigint(&K1, -1);
    NumberLength = K1.nbrLimbs;
    ComputeInversePower2(ValAOdd.limbs, tmp2.limbs, tmp1.limbs);
    setNbrLimbs(&tmp2);
    CopyBigInt(&tmp1, &ValB);
    BigIntDivideBy2(&tmp1);               // b/2
    (void)BigIntMultiply(&tmp1, &tmp2, &tmp1);  // b/2a
    BigIntChSign(&tmp1);                  // -b/2a
    BigIntAnd(&tmp1, &K1, &tmp1);         // -b/2a mod 2^expon
    BigIntAdd(&tmp1, &sqrRoot, &tmp2);
    BigIntAnd(&tmp2, &K1, &Solution1[factorIndex]);
    BigIntSubt(&tmp1, &sqrRoot, &tmp2);
    BigIntAnd(&tmp2, &K1, &Solution2[factorIndex]);
  }
  else if ((bitsAZero == 0) && (bitsBZero == 0))
  {
    CopyBigInt(&Quadr, &ValA);
    BigIntMultiplyBy2(&Quadr);         // 2a
    CopyBigInt(&Linear, &ValB);        // b
    CopyBigInt(&Const, &ValC);
    BigIntDivideBy2(&Const);           // c/2
    findQuadraticSolution(&Solution1[factorIndex], expon - 1);
    BigIntMultiplyBy2(&Solution1[factorIndex]);

    CopyBigInt(&Quadr, &ValA);
    BigIntMultiplyBy2(&Quadr);         // 2a
    BigIntAdd(&Quadr, &ValB, &Linear); // 2a+b
    CopyBigInt(&Const, &ValA);
    BigIntAdd(&Const, &ValB, &Const);
    BigIntAdd(&Const, &ValC, &Const);
    BigIntDivideBy2(&Const);           // (a+b+c)/2
    findQuadraticSolution(&Solution2[factorIndex], expon - 1);
    BigIntMultiplyBy2(&Solution2[factorIndex]);
    addbigint(&Solution2[factorIndex], 1);
  }
  else
  {
    CopyBigInt(&Quadr, &ValA);
    CopyBigInt(&Linear, &ValB);
    CopyBigInt(&Const, &ValC);
    findQuadraticSolution(&Solution1[factorIndex], expon);
    sol2Invalid = true;
  }
  BigIntPowerOf2(&Q, expon);         // Store increment.
  return true;
}

static void ComputeSquareRootModPowerOfP(int nbrBitsSquareRoot)
{
  int correctBits;
  NumberLength = prime.nbrLimbs;
  int NumberLengthBytes = NumberLength * (int)sizeof(limb);
  (void)memcpy(TestNbr, prime.limbs, NumberLengthBytes);
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
    limb* toConvert;
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
      int x;
      int e;
      int r;
      // Step 1.
      subtractdivide(&Q, 1, 1);   // Q <- (prime-1).
      DivideBigNbrByMaxPowerOf2(&e, Q.limbs, &Q.nbrLimbs);
      // Step 2.
      x = 1;
      do
      {
        x++;
        intToBigInteger(&Aux[3], x);
      } while (BigIntJacobiSymbol(&Aux[3], &prime) >= 0);
      // Step 3.
      // Get z <- x^q (mod p) in Montgomery notation.
      modPowBaseInt(x, Q.limbs, Q.nbrLimbs, Aux[4].limbs);  // z
      // Step 4.
      NumberLengthBytes = NumberLength * (int)sizeof(limb);
      (void)memcpy(Aux[5].limbs, Aux[4].limbs, NumberLengthBytes); // y
      r = e;
      CopyBigInt(&K1, &Q);
      subtractdivide(&K1, 1, 2);
      modPow(Aux[6].limbs, K1.limbs, K1.nbrLimbs, Aux[7].limbs); // x
      modmult(Aux[6].limbs, Aux[7].limbs, Aux[8].limbs);         // v
      modmult(Aux[8].limbs, Aux[7].limbs, Aux[9].limbs);         // w
      // Step 5
      while (memcmp(Aux[9].limbs, MontgomeryMultR1, NumberLengthBytes) != 0)
      {
        // Step 6
        int k = 0;
        (void)memcpy(Aux[10].limbs, Aux[9].limbs, NumberLengthBytes);
        do
        {
          k++;
          modmult(Aux[10].limbs, Aux[10].limbs, Aux[10].limbs);
        } while (memcmp(Aux[10].limbs, MontgomeryMultR1, NumberLengthBytes) != 0);
        // Step 7
        (void)memcpy(Aux[11].limbs, Aux[5].limbs, NumberLengthBytes); // d
        for (int ctr = 0; ctr < (r - k - 1); ctr++)
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
    NumberLengthBytes = NumberLength * (int)sizeof(limb);
    (void)memset(Aux[4].limbs, 0, NumberLengthBytes); // Convert power to standard notation.
    Aux[4].limbs[0].x = 1;
    modmult(Aux[4].limbs, toConvert, toConvert);
    UncompressLimbsBigInteger(toConvert, &SqrtDisc);
  }
  // Obtain inverse of square root stored in SqrtDisc (mod prime).
  intToBigInteger(&tmp2, 1);
  BigIntModularDivision(&tmp2, &SqrtDisc, &prime, &sqrRoot);
  correctBits = 1;
  CopyBigInt(&Q, &prime);
  // Obtain nbrBitsSquareRoot correct digits of inverse square root.
  while (correctBits < nbrBitsSquareRoot)
  {   // Compute f(x) = invsqrt(x), f_{n+1}(x) = f_n * (3 - x*f_n^2)/2
    correctBits *= 2;
    (void)BigIntMultiply(&Q, &Q, &Q);           // Square Q.
    (void)BigIntMultiply(&sqrRoot, &sqrRoot, &tmp1);
    (void)BigIntRemainder(&tmp1, &Q, &tmp2);
    (void)BigIntMultiply(&tmp2, &discriminant, &tmp1);
    (void)BigIntRemainder(&tmp1, &Q, &tmp2);
    intToBigInteger(&tmp1, 3);
    BigIntSubt(&tmp1, &tmp2, &tmp2);
    (void)BigIntMultiply(&tmp2, &sqrRoot, &tmp1);
    if ((tmp1.limbs[0].x & 1) != 0)
    {
      BigIntAdd(&tmp1, &Q, &tmp1);
    }
    BigIntDivide2(&tmp1);
    (void)BigIntRemainder(&tmp1, &Q, &sqrRoot);
  }
  // Get square root of discriminant from its inverse by multiplying by discriminant.
  if (sqrRoot.sign == SIGN_NEGATIVE)
  {
    BigIntAdd(&sqrRoot, &Q, &sqrRoot);
  }
  (void)BigIntMultiply(&sqrRoot, &discriminant, &sqrRoot);
  (void)BigIntRemainder(&sqrRoot, &Q, &sqrRoot);
}

// Solve Ax^2 + Bx + C = 0 (mod p^expon).
static bool SolveQuadraticEqModPowerOfP(int expon, int factorIndex)
{
  int correctBits;
  int nbrLimbs;
  int bitsAZero;
  int ctr;
  int deltaZeros;
  int NumberLengthBytes;
  // Number of bits of square root of discriminant to compute: expon + bits_a + 1,
  // where bits_a is the number of least significant bits of a set to zero.
  // To compute the square root, compute the inverse of sqrt, so only multiplications are used.
  // f(x) = invsqrt(x), f_{n+1}(x) = f_n * (3 - x*f_n^2)/2
  // Get maximum power of prime which divide ValA.
  CopyBigInt(&ValAOdd, &ValA);
  bitsAZero = 0;
  for (;;)
  {
    (void)BigIntRemainder(&ValAOdd, &prime, &tmp1);
    if (ValAOdd.sign == SIGN_NEGATIVE)
    {
      BigIntAdd(&ValAOdd, &prime, &ValAOdd);
    }
    if (!BigIntIsZero(&tmp1))
    {
      break;
    }
    (void)BigIntDivide(&ValAOdd, &prime, &ValAOdd);
    bitsAZero++;
  }
  (void)BigIntRemainder(&discriminant, &V, &discriminant);
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
      (void)BigIntRemainder(&discriminant, &prime, &tmp1);
      if (!BigIntIsZero(&tmp1))
      {
        break;
      }
      (void)BigIntDivide(&discriminant, &prime, &discriminant);
      deltaZeros++;
    }
  }
  if (((deltaZeros & 1) != 0) && (deltaZeros < expon))
  {          // If delta is of type m*prime^n where m is not multiple of prime
             // and n is odd, there is no solution, so go out.
    return false;
  }
  deltaZeros >>= 1;
  // Compute inverse of -2*A (mod prime^(expon - deltaZeros)).
  BigIntAdd(&ValA, &ValA, &ValAOdd);
  (void)BigIntPowerIntExp(&prime, expon - deltaZeros, &tmp1);
  (void)BigIntRemainder(&ValAOdd, &tmp1, &ValAOdd);
  nbrLimbs = tmp1.nbrLimbs;
  if (ValAOdd.sign == SIGN_NEGATIVE)
  {
    ValAOdd.sign = SIGN_POSITIVE;           // Negate 2*A
  }
  else if (!BigIntIsZero(&ValAOdd))
  {
    BigIntSubt(&tmp1, &ValAOdd, &ValAOdd);  // Negate 2*A
  }
  else
  {            // Nothing to do.
  }
  intToBigInteger(&tmp2, 1);
  NumberLength = tmp1.nbrLimbs;
  NumberLengthBytes = NumberLength * (int)sizeof(limb);
  (void)memcpy(TestNbr, tmp1.limbs, NumberLengthBytes);
  TestNbr[NumberLength].x = 0;
  GetMontgomeryParms(NumberLength);
  BigIntModularDivision(&tmp2, &ValAOdd, &tmp1, &Aux[0]);
  CopyBigInt(&ValAOdd, &Aux[0]);
  if (BigIntIsZero(&discriminant))
  {     // Discriminant is zero.
    int lenBytes = nbrLimbs * (int)sizeof(limb);
    (void)memset(sqrRoot.limbs, 0, lenBytes);
  }
  else
  {      // Discriminant is not zero.
    // Find number of digits of square root to compute.
    int nbrBitsSquareRoot = expon + bitsAZero - deltaZeros;
    (void)BigIntPowerIntExp(&prime, nbrBitsSquareRoot, &tmp1);
    nbrLimbs = tmp1.nbrLimbs;
    (void)BigIntRemainder(&discriminant, &tmp1, &discriminant);
    if (discriminant.sign == SIGN_NEGATIVE)
    {
      BigIntAdd(&discriminant, &tmp1, &discriminant);
    }
    if (nbrLimbs > discriminant.nbrLimbs)
    {
      int lenBytes = (nbrLimbs - discriminant.nbrLimbs) * (int)sizeof(limb);
      (void)memset(&discriminant.limbs[nbrLimbs], 0, lenBytes);
    }
    (void)BigIntRemainder(&discriminant, &prime, &Aux[3]);
    if (Aux[3].sign == SIGN_NEGATIVE)
    {
      BigIntAdd(&Aux[3], &prime, &Aux[3]);
    }
    if (BigIntJacobiSymbol(&Aux[3], &prime) != 1)
    {
      return false;         // Not a quadratic residue, so go out.
    }
    // Compute square root of discriminant.
    ComputeSquareRootModPowerOfP(nbrBitsSquareRoot);
    // Multiply by square root of discriminant by prime^deltaZeros.
    for (ctr = 0; ctr < deltaZeros; ctr++)
    {
      (void)BigIntMultiply(&sqrRoot, &prime, &sqrRoot);
    }
  }
  correctBits = expon - deltaZeros;
  (void)BigIntPowerIntExp(&prime, correctBits, &Q);      // Store increment.
  // Compute x = (b + sqrt(discriminant)) / (-2a) and x = (b - sqrt(discriminant)) / (-2a)
  BigIntAdd(&ValB, &sqrRoot, &tmp1);
  for (ctr = 0; ctr < bitsAZero; ctr++)
  {
    (void)BigIntRemainder(&tmp1, &prime, &tmp2);
    if (!BigIntIsZero(&tmp2))
    {   // Cannot divide by prime, so go out.
      sol1Invalid = true;
      break;
    }
    (void)BigIntDivide(&tmp1, &prime, &tmp1);
  }
  (void)BigIntMultiply(&tmp1, &ValAOdd, &tmp1);
  (void)BigIntRemainder(&tmp1, &Q, &Solution1[factorIndex]);
  if (Solution1[factorIndex].sign == SIGN_NEGATIVE)
  {
    BigIntAdd(&Solution1[factorIndex], &Q, &Solution1[factorIndex]);
  }
  BigIntSubt(&ValB, &sqrRoot, &tmp1);
  for (ctr = 0; ctr < bitsAZero; ctr++)
  {
    (void)BigIntRemainder(&tmp1, &prime, &tmp2);
    if (!BigIntIsZero(&tmp2))
    {   // Cannot divide by prime, so go out.
      sol2Invalid = true;
      break;
    }
    (void)BigIntDivide(&tmp1, &prime, &tmp1);
  }
  (void)BigIntMultiply(&tmp1, &ValAOdd, &tmp1);
  (void)BigIntRemainder(&tmp1, &Q, &Solution2[factorIndex]);
  if (Solution2[factorIndex].sign == SIGN_NEGATIVE)
  {
    BigIntAdd(&Solution2[factorIndex], &Q, &Solution2[factorIndex]);
  }
  return true;
}

static void QuadraticTermMultipleOfP(int expon, int factorIndex)
{
  // Perform Newton approximation:
  // x_{n+1} = x_n - (a*x_n^2 + b*x_n + c) / (2*a_x + b)      
  BigInteger* ptrSolution = &Solution1[factorIndex];
  int NumberLengthBytes;
  NumberLength = prime.nbrLimbs;
  NumberLengthBytes = NumberLength * (int)sizeof(limb);
  (void)memcpy(TestNbr, prime.limbs, NumberLengthBytes);
  TestNbr[NumberLength].x = 0;
  GetMontgomeryParms(NumberLength);
  BigIntModularDivision(&ValC, &ValB, &prime, ptrSolution);
  BigIntNegate(ptrSolution, ptrSolution);
  if (ptrSolution->sign == SIGN_NEGATIVE)
  {
    BigIntAdd(ptrSolution, &prime, ptrSolution);
  }
  for (int currentExpon = 2; currentExpon < (2 * expon); currentExpon *= 2)
  {
    (void)BigIntPowerIntExp(&prime, currentExpon, &V);
    (void)BigIntMultiply(&ValA, ptrSolution, &Q);// a*x_n
    CopyBigInt(&L, &Q);
    BigIntAdd(&Q, &ValB, &Q);                    // a*x_n + b
    (void)BigIntRemainder(&Q, &V, &Q);
    (void)BigIntMultiply(&Q, ptrSolution, &Q);   // a*x_n^2 + b*x_n
    BigIntAdd(&Q, &ValC, &Q);                    // a*x_n^2 + b*x_n + c
    (void)BigIntRemainder(&Q, &V, &Q);           // Numerator. 
    multint(&L, &L, 2);                          // 2*a*x_n
    BigIntAdd(&L, &ValB, &L);                    // 2*a*x_n + b
    (void)BigIntRemainder(&L, &V, &L);           // Denominator
    NumberLength = V.nbrLimbs;
    NumberLengthBytes = NumberLength * (int)sizeof(limb);
    (void)memcpy(TestNbr, V.limbs, NumberLengthBytes);
    TestNbr[NumberLength].x = 0;
    GetMontgomeryParms(NumberLength);
    BigIntModularDivision(&Q, &L, &V, &z);
    BigIntSubt(ptrSolution, &z, ptrSolution);
    (void)BigIntRemainder(ptrSolution, &V, ptrSolution);
    if (ptrSolution->sign == SIGN_NEGATIVE)
    {
      BigIntAdd(ptrSolution, &V, ptrSolution);
    }
  }
  (void)BigIntPowerIntExp(&prime, expon, &Q);
  (void)BigIntRemainder(ptrSolution, &Q, ptrSolution);
  CopyBigInt(&Solution2[factorIndex], &Solution1[factorIndex]);
}

static bool QuadraticTermNotMultipleOfP(int expon, int factorIndex)
{
  bool solutions;
  sol1Invalid = false;
  sol2Invalid = false;
  // Compute discriminant = ValB^2 - 4*ValA*ValC.
  (void)BigIntMultiply(&ValB, &ValB, &Aux[0]);
  (void)BigIntMultiply(&ValA, &ValC, &discriminant);
  multint(&discriminant, &discriminant, 4);
  BigIntSubt(&Aux[0], &discriminant, &discriminant);
  if ((prime.nbrLimbs == 1) && (prime.limbs[0].x == 2))
  {         /* Prime p is 2 */
    solutions = SolveQuadraticEqModPowerOf2(expon, factorIndex);
  }
  else
  {                        // Prime is not 2
    solutions = SolveQuadraticEqModPowerOfP(expon, factorIndex);
  }
  if (!solutions || (sol1Invalid && sol2Invalid))
  {     // Both solutions are invalid. Go out.
    return false;
  }
  if (sol1Invalid)
  {     // Solution1 is invalid. Overwrite it with Solution2.
    CopyBigInt(&Solution1[factorIndex], &Solution2[factorIndex]);
  }
  else if (sol2Invalid)
  {     // Solution2 is invalid. Overwrite it with Solution1.
    CopyBigInt(&Solution2[factorIndex], &Solution1[factorIndex]);
  }
  else
  {              // Nothing to do.
  }
  BigIntSubt(&Solution2[factorIndex], &Solution1[factorIndex], &Aux[0]);
  if (Aux[0].sign == SIGN_NEGATIVE)
  {     // Solution2 is less than Solution1, so exchange them.
    CopyBigInt(&Aux[0], &Solution1[factorIndex]);
    CopyBigInt(&Solution1[factorIndex], &Solution2[factorIndex]);
    CopyBigInt(&Solution2[factorIndex], &Aux[0]);
  }
  return true;
}

// Solve Ax^2 + Bx + C = 0 (mod N).
static void SolveEquation(void)
{
  const struct sFactors *pstFactor;

  if (BigIntIsZero(&ValN))
  {        // Mod zero => Equation in integer numbers
    SolveIntegerEquation();
    return;
  }
  // Modulus is not zero.
  ValN.sign = SIGN_POSITIVE;
  (void)BigIntRemainder(&ValA, &ValN, &ValA);
  (void)BigIntRemainder(&ValB, &ValN, &ValB);
  (void)BigIntRemainder(&ValC, &ValN, &ValC);
  BigIntGcd(&ValA, &ValB, &Aux[0]);
  Aux[0].sign = SIGN_POSITIVE;
  BigIntGcd(&ValC, &Aux[0], &GcdAll);
  GcdAll.sign = SIGN_POSITIVE;
  (void)BigIntRemainder(&ValC, &GcdAll, &Aux[0]);
  if (!BigIntIsZero(&Aux[0]))
  {  // ValC must be multiple of gcd(ValA, ValB).
     // Otherwise go out because there are no solutions.
    return;
  }
  BigIntGcd(&ValN, &GcdAll, &Aux[0]);
  CopyBigInt(&GcdAll, &Aux[0]);
  // Divide all coefficients by gcd(ValA, ValB).
  (void)BigIntDivide(&ValA, &GcdAll, &ValA);
  (void)BigIntDivide(&ValB, &GcdAll, &ValB);
  (void)BigIntDivide(&ValC, &GcdAll, &ValC);
  (void)BigIntDivide(&ValN, &GcdAll, &ValN);
  CopyBigInt(&ValNn, &ValN);
  if ((ValNn.nbrLimbs == 1) && (ValNn.limbs[0].x == 1))
  {     // All values from 0 to GcdAll - 1 are solutions.
    if ((GcdAll.nbrLimbs > 1) || (GcdAll.limbs[0].x > 5))
    {
      copyStr(&ptrOutput, "<p>All values of <var>x</var> between 0 and ");
      addbigint(&GcdAll, -1);
      BigInteger2Dec(&ptrOutput, &GcdAll, groupLen);
      copyStr(&ptrOutput, " are solutions.</p>");
    }
    else
    {
      for (int ctr = 0; ctr < GcdAll.limbs[0].x; ctr++)
      {
        intToBigInteger(&Aux[0], ctr);
        Solution(&Aux[0]);
      }
    }
    return;
  }
  (void)BigIntRemainder(&ValA, &ValN, &Aux[0]);
  if (BigIntIsZero(&Aux[0]))
  {           // Linear equation.
    SolveModularLinearEquation();
    return;
  }
  if (!BigIntEqual(&LastModulus, &ValN))
  {     // Last modulus is different from ValN.
    char* ptrFactorDec;
    CopyBigInt(&LastModulus, &ValN);
    NumberLength = ValN.nbrLimbs;
    BigInteger2IntArray(nbrToFactor, &ValN);
    ptrFactorDec = tofactorDec;
    Bin2Dec(&ptrFactorDec, ValN.limbs, ValN.nbrLimbs, groupLen);
    factor(&ValN, nbrToFactor, factorsMod, astFactorsMod);
  }
  intToBigInteger(&Q, 0);
  nbrFactors = astFactorsMod[0].multiplicity;
  pstFactor = &astFactorsMod[1];
  for (int factorIndex = 0; factorIndex<nbrFactors; factorIndex++)
  {
    int expon = pstFactor->multiplicity;
    NumberLength = *pstFactor->ptrFactor;
    IntArray2BigInteger(pstFactor->ptrFactor, &prime);
    (void)BigIntPowerIntExp(&prime, expon, &V);
    (void)BigIntRemainder(&ValA, &prime, &L);
    if (BigIntIsZero(&L) && !((prime.nbrLimbs == 1) && (prime.limbs[0].x == 2)))
    {     // ValA multiple of prime, means linear equation mod prime. Also prime is not 2.
      if ((BigIntIsZero(&ValB)) && !(BigIntIsZero(&ValC)))
      {
        return;  // There are no solutions: ValB=0 and ValC!=0
      }
      QuadraticTermMultipleOfP(expon, factorIndex);
    }
    else
    {                   /* If quadratic equation mod p */
      if (!QuadraticTermNotMultipleOfP(expon, factorIndex))
      {
        return;
      }
    }
    CopyBigInt(&Increment[factorIndex], &Q);
    Exponents[factorIndex] = 0;
    pstFactor++;
  }   // end for
  PerformChineseRemainderTheorem();
}

void textErrorQuadMod(char **pptrOutput, enum eExprErr rc)
{
  if (rc == EXPR_MODULUS_MUST_BE_NONNEGATIVE)
  {
    copyStr(pptrOutput, lang ? "No debe ser negativo" :
      "Must not be negative");
  }
  else
  {
    textError(pptrOutput, rc);
  }
}

void quadmodText(const char *quadrText, const char *linearText, const char *constText,
  const char *modText, int groupLength)
{
  groupLen = groupLength;
  char *ptrBeginSol;
  enum eExprErr rc;
  ptrOutput = output;
  *ptrOutput = '2';  // Prefix for sending output to screen.
  ptrOutput++;
  rc = ComputeExpression(quadrText, &ValA, false);
  if (rc != EXPR_OK)
  {
    copyStr(&ptrOutput, lang ? "<p>Coeficiente cuadrático: ": "<p>Quadratic coefficient: ");
    textErrorQuadMod(&ptrOutput, rc);
    copyStr(&ptrOutput, "</p>");
  }
  rc = ComputeExpression(linearText, &ValB, false);
  if (rc != EXPR_OK)
  {
    copyStr(&ptrOutput, lang ? "<p>Coeficiente lineal: " : "<p>Linear coefficient: ");
    textErrorQuadMod(&ptrOutput, rc);
    copyStr(&ptrOutput, "</p>");
  }
  rc = ComputeExpression(constText, &ValC, false);
  if (rc != EXPR_OK)
  {
    copyStr(&ptrOutput, lang ? "<p>Término independiente: " : "<p>Constant coefficient: ");
    textErrorQuadMod(&ptrOutput, rc);
    copyStr(&ptrOutput, "</p>");
  }
  rc = ComputeExpression(modText, &ValN, false);
  if ((rc == EXPR_OK) && (ValN.sign == SIGN_NEGATIVE))
  {
    rc = EXPR_MODULUS_MUST_BE_NONNEGATIVE;
  }
  if (rc != EXPR_OK)
  {
    copyStr(&ptrOutput, lang ? "Módulo: " : "Modulus: ");
    textErrorQuadMod(&ptrOutput, rc);
    copyStr(&ptrOutput, "</p>");
  }
  if (ptrOutput == &output[1])
  {    // No errors found.
    copyStr(&ptrOutput, "<p>");
    int u = Show(&ValA, " x&sup2;", 2);
    u = Show(&ValB, " x", u);
    Show1(&ValC, u);
    copyStr(&ptrOutput, " &equiv; 0 (mod ");
    BigInteger2Dec(&ptrOutput, &ValN, groupLen);
    copyStr(&ptrOutput, ")</p>");
    SolNbr = 0;
    ptrBeginSol = ptrOutput;
    copyStr(&ptrOutput, "<ol>");
    SolveEquation();
    if (SolNbr == 0)
    {
      ptrOutput = ptrBeginSol;
      copyStr(&ptrOutput, lang? "<p>No hay soluciones.</p>": "<p>There are no solutions.</p>");
    }
    else
    {
      copyStr(&ptrOutput, "</ol>");
    }
  }
  copyStr(&ptrOutput, "<p>");
  copyStr(&ptrOutput, lang ? COPYRIGHT_SPANISH: COPYRIGHT_ENGLISH);
  copyStr(&ptrOutput, "</p>");
}

#if defined __EMSCRIPTEN__ && !defined _MSC_VER
EXTERNALIZE void doWork(void)
{
  int flags;
  char *ptrData = inputString;
  char* ptrQuadrCoeff;
  char* ptrLinearCoeff;
  char* ptrConstCoeff;
  char* ptrMod;
  int groupLength = 0;
  while (*ptrData != ',')
  {
    groupLength = (groupLength * 10) + (*ptrData - '0');
    ptrData++;
  }
  ptrData++;                    // Skip comma.
  flags = *ptrData;
#ifndef lang  
  lang = ((flags & 1)? true: false);
#endif
  ptrQuadrCoeff = ptrData + 2;  // Skip flags and comma.
  ptrLinearCoeff = ptrQuadrCoeff + (int)strlen(ptrQuadrCoeff) + 1;
  ptrConstCoeff = ptrLinearCoeff + (int)strlen(ptrLinearCoeff) + 1;
  ptrMod = ptrConstCoeff + (int)strlen(ptrConstCoeff) + 1;
  quadmodText(ptrQuadrCoeff, ptrLinearCoeff, ptrConstCoeff, ptrMod, groupLength);
  databack(output);
}
#endif

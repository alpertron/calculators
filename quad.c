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

enum eLinearSolution
{
  SOLUTION_FOUND = 0,
  NO_SOLUTIONS,
  INFINITE_SOLUTIONS,
};

static enum eShowSolution
{
  ONE_SOLUTION = 0,
  TWO_SOLUTIONS,
  FIRST_SOLUTION,
  SECOND_SOLUTION,
} showSolution;

static enum eCallbackQuadModType
{
  CBACK_QMOD_PARABOLIC = 0,
  CBACK_QMOD_ELLIPTIC,
  CBACK_QMOD_HYPERBOLIC,
} callbackQuadModType;

static int counters[400];
static char isDescending[400];
static int originalMultiplicities[400];
static BigInteger Solution1[400];
static BigInteger Solution2[400];
static BigInteger Increment[400];
static int Exponents[400];
static BigInteger Aux[400];
static BigInteger ValAbak, ValBbak, ValCbak;
static BigInteger ValA, ValB, ValC, ValD, ValE, ValF;
static BigInteger ValH, ValI, ValL, ValM, ValN, ValO, ValP, ValQ;
static BigInteger ValU, ValV, ValG, ValR, ValS, ValJ, ValK, ValZ;
static BigInteger ValAlpha, ValBeta, ValDen, ValDiv;
static BigInteger ValUBak, ValVBak;
static BigInteger ValGcdHomog;
static BigInteger Tmp1, Tmp2, Tmp3;
static int SolNbr;
static int showRecursiveSolution;
static char textExp[1000];
static BigInteger Xind, Yind, Xlin, Ylin;
static int nbrFactors, solFound;
static char teach = 0, also;
static char ExchXY;
static char *ptrOutput;
static char *divgcd = "<p>Dividing the equation by the greatest common divisor we obtain:</p>";
static char *varT = "t";
static char *varTprime = "t'";
static char *squareText = "&sup2;";
static BigInteger discr;
static BigInteger U1, U2, U3, V1, V2, V3, bigTmp;
static BigInteger startPeriodU, startPeriodV;
static BigInteger coeffQuadr, coeffLinear, coeffIndep, modulus;
static BigInteger currentFactor, prime;
static BigInteger Xplus, Xminus, Yplus, Yminus;
static BigInteger *Xbak, *Ybak;
static int Show(BigInteger *num, char *str, enum eLinearSolution t);
static void Show1(BigInteger *num, enum eLinearSolution t);
static void callbackQuadModParabolic(BigInteger *value);
static void callbackQuadModElliptic(BigInteger *value);
static void callbackQuadModHyperbolic(BigInteger *value);
static void ContFracPell(void);
static BigInteger Tmp[12];
static int indexEvenMultiplicity[400];
static int nbrPrimesEvenMultiplicity;
static BigInteger K, L;
#define NBR_COEFF 6
struct stValidateCoeff
{
  char *expression;
  BigInteger *bigint;
  char *textSpanish;
  char *textEnglish;
};
static struct stValidateCoeff astValidateCoeff[NBR_COEFF] =
{
  { NULL, &ValA, "Coeficiente <var>a</var>: ", "Coefficient <var>a</var>: " },
  { NULL, &ValB, "Coeficiente <var>b</var>: ", "Coefficient <var>b</var>: " },
  { NULL, &ValC, "Coeficiente <var>c</var>: ", "Coefficient <var>c</var>: " },
  { NULL, &ValD, "Coeficiente <var>d</var>: ", "Coefficient <var>d</var>: " },
  { NULL, &ValE, "Coeficiente <var>e</var>: ", "Coefficient <var>e</var>: " },
  { NULL, &ValF, "Coeficiente <var>f</var>: ", "Coefficient <var>f</var>: " },
};

static void showMinus(void)
{
  strcpy(ptrOutput, "&#8209;");
  ptrOutput += strlen(ptrOutput);
}

static void shownbr(BigInteger *value)
{
  BigInteger2Dec(value, ptrOutput, groupLen);
  ptrOutput += strlen(ptrOutput);
}

static void showInt(int value)
{
  int2dec(&ptrOutput, value);
  ptrOutput += strlen(ptrOutput);
}

static void showSquare(void)
{
  strcpy(ptrOutput, "&sup2;");
  ptrOutput += strlen(ptrOutput);
}

static void w(char *text)
{
  strcpy(ptrOutput, text);
  ptrOutput += strlen(ptrOutput);
}

static void showAlso(void)
{
  if (also)
  {
    w("and also:<BR>");
  }
  else
  {
    also = 1;
  }
}

static void ShowLin(BigInteger *coeffX, BigInteger *coeffY, BigInteger *coeffInd, char *x, char *y)
{
  enum eLinearSolution t;
  t = Show(coeffX, x, SOLUTION_FOUND);
  t = Show(coeffY, y, t);
  Show1(coeffInd, t);
}

static void ShowLinInd(BigInteger *lin, BigInteger *ind, char *var)
{
  if (BigIntIsZero(ind) && BigIntIsZero(lin))
  {
    w("0");
  }
  if (!BigIntIsZero(ind))
  {
    shownbr(ind);
  }
  *ptrOutput++ = ' ';
  if (lin->sign == SIGN_NEGATIVE)
  {
    showMinus();
  }
  else if (!BigIntIsZero(lin) && !BigIntIsZero(ind))
  {
    *ptrOutput++ = '+';
  }
  *ptrOutput++ = ' ';
  if (!BigIntIsZero(lin))
  {
    if (lin->nbrLimbs != 1 || lin->limbs[0].x != 1)
    {     // abs(lin) is not 1
      CopyBigInt(&Aux[0], lin);
      shownbr(&Aux[0]);
    }
    *ptrOutput++ = ' ';
    strcpy(ptrOutput, var);
    ptrOutput += strlen(ptrOutput);
  }
}

int PrintLinear(enum eLinearSolution Ret, char *var)
{
  if (Ret == NO_SOLUTIONS)
  {
    return 0;
  }
  if (var == varT)
  {
    showAlso();
  }
  if (Ret == INFINITE_SOLUTIONS)
  {
    w("<p>x, y: any integer</p>");
    return 1;
  }
  if (ExchXY)
  {
    // Exchange Xind and Yind
    CopyBigInt(&bigTmp, &Xind);
    CopyBigInt(&Xind, &Yind);
    CopyBigInt(&Yind, &bigTmp);
    // Exchange Xlin and Ylin
    CopyBigInt(&bigTmp, &Xlin);
    CopyBigInt(&Xlin, &Ylin);
    CopyBigInt(&Ylin, &bigTmp);
  }
  w("<p>x = ");
  ShowLinInd(&Xlin, &Xind, var);
  w("<br>y = ");
  ShowLinInd(&Ylin, &Yind, var);
  w("</p>");
  return 0;
}

static void PrintQuad(BigInteger *coeffT2, BigInteger *coeffT, BigInteger *coeffInd)
{
  if (coeffT2->nbrLimbs == 1 && coeffT2->limbs[0].x == 1)
  {             // abs(coeffT2) = 1
    if (coeffT2->sign == SIGN_POSITIVE)
    {           // coeffT2 = 1
      *ptrOutput++ = ' ';
    }
    else
    {           // coeffT2 = -1
      showMinus();
    }
    *ptrOutput++ = 't';
    showSquare();
  }
  else if (!BigIntIsZero(coeffT2))
  {             // coeffT2 is not zero.
    shownbr(coeffT2);
    *ptrOutput++ = ' ';
    *ptrOutput++ = 't';
    showSquare();
  }
  if (coeffT->sign == SIGN_NEGATIVE)
  {
    w(" - ");
  }
  else if (!BigIntIsZero(coeffT) && !BigIntIsZero(coeffT2))
  {
    w(" + ");
  }
  if (coeffT->nbrLimbs == 1 && coeffT->limbs[0].x == 1)
  {     // abs(coeffT) = 1
    w("t ");
  }
  else if (!BigIntIsZero(coeffT))
  {
    if (coeffT->sign == SIGN_NEGATIVE)
    {
      coeffT->sign = SIGN_POSITIVE;
      shownbr(coeffT);
      coeffT->sign = SIGN_NEGATIVE;
    }
    else
    {
      shownbr(coeffT);
    }
    w(" t");
  }
  if (coeffInd->sign == SIGN_NEGATIVE)
  {
    w(" - ");
    coeffInd->sign = SIGN_POSITIVE;
    shownbr(coeffInd);
    coeffInd->sign = SIGN_NEGATIVE;
  }
  else if (!BigIntIsZero(coeffInd))
  {
    if (!BigIntIsZero(coeffT) || !BigIntIsZero(coeffT2))
    {
      w(" + ");
    }  
    shownbr(coeffInd);
  }
  w("<br>");
}

void ShowXY(BigInteger *X, BigInteger *Y)
{
  enum eSign signX, signY;
  if (showSolution == TWO_SOLUTIONS)
  {
    solFound = 1;
    if (Xbak->nbrLimbs == 0)
    {
      CopyBigInt(Xbak, X);
      CopyBigInt(Ybak, Y);
      return;
    }
    // Use the lowest of |X| + |Y| and |Xbak| + |Ybak|
    signX = Xbak->sign;
    signY = Ybak->sign;
    Xbak->sign = SIGN_POSITIVE;
    Ybak->sign = SIGN_POSITIVE;
    BigIntAdd(Xbak, Ybak, &bigTmp);
    Xbak->sign = signX;
    Ybak->sign = signY;
    signX = X->sign;
    signY = Y->sign;
    X->sign = SIGN_POSITIVE;
    Y->sign = SIGN_POSITIVE;
    BigIntSubt(&bigTmp, X, &bigTmp);
    BigIntSubt(&bigTmp, Y, &bigTmp);
    X->sign = signX;
    Y->sign = signY;
    if (bigTmp.sign == SIGN_POSITIVE)
    {       // |x| + |y| <= |xbak| + |ybak|
      CopyBigInt(Xbak, X);
      CopyBigInt(Ybak, Y);
    }
  }
  else
  {
    // ONE_SOLUTION: Show it.
    showAlso();
    w("<p>x = ");
    shownbr(ExchXY ? Y : X);
    w("<BR>y = ");
    shownbr(ExchXY ? X : Y);
    w("</p>");
  }
}

static int Show(BigInteger *num, char *str, enum eLinearSolution t)
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
      strcpy(ptrOutput, "&nbsp;&#8290;");
    }
    else
    {
      strcpy(ptrOutput, "&nbsp;");
    }
    if (*str)
    {
      ptrOutput += strlen(ptrOutput);
    }
    strcpy(ptrOutput, str);
    ptrOutput += strlen(ptrOutput);
    return t | 1;
  }
  return t;
}

static void Show1(BigInteger *num, enum eLinearSolution t)
{
  int u = Show(num, "", t);
  *ptrOutput++ = ' ';
  if ((u & 1) == 0 || (num->nbrLimbs == 1 && num->limbs[0].x == 1))
  {
    enum eSign sign = num->sign;
    num->sign = SIGN_POSITIVE;
    BigInteger2Dec(num, ptrOutput, groupLen);
    num->sign = sign;
    ptrOutput += strlen(ptrOutput);
  }
}

static void ShowEq(BigInteger *coeffA, BigInteger *coeffB, BigInteger *coeffC, 
  BigInteger *coeffD, BigInteger *coeffE, BigInteger *coeffF, char *x, char *y)
{
  char var[30];
  char *ptrVar;
  enum eLinearSolution t;
  ptrVar = var;
  strcpy(ptrVar, x);
  ptrVar += strlen(ptrVar);
  strcpy(ptrVar, squareText);
  ptrVar += strlen(ptrVar);
  t = Show(coeffA, var, SOLUTION_FOUND);

  ptrVar = var;
  strcpy(ptrVar, x);
  ptrVar += strlen(ptrVar);
  strcpy(ptrVar, "&#8290;");
  ptrVar += strlen(ptrVar);
  strcpy(ptrVar, y);
  ptrVar += strlen(ptrVar);
  t = Show(coeffB, var, t);

  ptrVar = var;
  strcpy(ptrVar, y);
  ptrVar += strlen(ptrVar);
  strcpy(ptrVar, squareText);
  ptrVar += strlen(ptrVar);
  t = Show(coeffC, var, t);

  ptrVar = var;
  strcpy(ptrVar, x);
  ptrVar += strlen(ptrVar);
  t = Show(coeffD, var, t);

  ptrVar = var;
  strcpy(ptrVar, y);
  ptrVar += strlen(ptrVar);
  t = Show(coeffE, var, t);
  Show1(coeffF, t);
}

static void SolutionX(BigInteger *value)
{
  SolNbr++;
  if (teach)
  {
    strcpy(ptrOutput, "<li>t = ");
    ptrOutput += strlen(ptrOutput);
    BigInteger2Dec(value, ptrOutput, groupLen);
    ptrOutput += strlen(ptrOutput);
    strcpy(ptrOutput, "</li>");
    ptrOutput += strlen(ptrOutput);
  }
  if (callbackQuadModType == CBACK_QMOD_PARABOLIC)
  {
    callbackQuadModParabolic(value);
  }
  else if (callbackQuadModType == CBACK_QMOD_ELLIPTIC)
  {
    callbackQuadModElliptic(value);
  }
  else if (callbackQuadModType == CBACK_QMOD_HYPERBOLIC)
  {
    callbackQuadModHyperbolic(value);
  }
}

// Solve congruence an^2 + bn + c = 0 (mod n) where n is different from zero.
void SolveQuadModEquation(void)
{
  int factorIndex, expon, T1, E;
  BigInteger GcdAll;
  BigInteger ValNn;
  BigInteger z, Mult, currentSolution;
  BigInteger U, SqrtDisc, K1;
  BigInteger discriminant, V, Q;
  struct sFactors *pstFactor;

  modulus.sign = SIGN_POSITIVE;
  BigIntRemainder(&coeffQuadr, &modulus, &coeffQuadr);
  if (coeffQuadr.sign == SIGN_NEGATIVE)
  {
    BigIntAdd(&coeffQuadr, &modulus, &coeffQuadr);
  }
  BigIntRemainder(&coeffLinear, &modulus, &coeffLinear);
  if (coeffLinear.sign == SIGN_NEGATIVE)
  {
    BigIntAdd(&coeffLinear, &modulus, &coeffLinear);
  }
  BigIntRemainder(&coeffIndep, &modulus, &coeffIndep);
  if (coeffIndep.sign == SIGN_NEGATIVE)
  {
    BigIntAdd(&coeffIndep, &modulus, &coeffIndep);
  }
  BigIntGcd(&coeffQuadr, &coeffLinear, &Tmp[0]);
  BigIntGcd(&coeffIndep, &Tmp[0], &GcdAll);
  BigIntRemainder(&coeffIndep, &GcdAll, &Tmp[0]);
  if (!BigIntIsZero(&Tmp[0]))
  {  // ValC must be multiple of gcd(ValA, ValB).
     // Otherwise go out because there are no solutions.
    return;
  }
  BigIntGcd(&modulus, &GcdAll, &Tmp[0]);
  CopyBigInt(&GcdAll, &Tmp[0]);
  // Divide all coefficients by gcd(ValA, ValB).
  BigIntDivide(&coeffQuadr, &GcdAll, &Tmp[0]);
  CopyBigInt(&coeffQuadr, &Tmp[0]);
  BigIntDivide(&coeffLinear, &GcdAll, &Tmp[0]);
  CopyBigInt(&coeffLinear, &Tmp[0]);
  BigIntDivide(&coeffIndep, &GcdAll, &Tmp[0]);
  CopyBigInt(&coeffIndep, &Tmp[0]);
  BigIntDivide(&modulus, &GcdAll, &Tmp[0]);
  CopyBigInt(&modulus, &Tmp[0]);
  CopyBigInt(&ValNn, &modulus);
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
        intToBigInteger(&Tmp[0], ctr);
        SolutionX(&Tmp[0]);
      }
    }
    return;
  }
  BigIntRemainder(&coeffQuadr, &modulus, &Tmp[0]);
  if (BigIntIsZero(&Tmp[0]))
  {           // Linear equation.
    BigIntGcd(&coeffLinear, &modulus, &Tmp[0]);
    if (Tmp[0].nbrLimbs != 1 || Tmp[0].limbs[0].x != 1)
    {         // ValB and ValN are not coprime. Go out.
      return;
    }
    // Calculate z <- -ValC / ValB (mod ValN)
    NumberLength = modulus.nbrLimbs;
    memcpy(TestNbr, modulus.limbs, NumberLength * sizeof(limb));
    TestNbr[NumberLength].x = 0;
    GetMontgomeryParms(NumberLength);
    BigIntModularDivision(&coeffIndep, &coeffLinear, &modulus, &z);
    BigIntNegate(&z, &z);
    if (z.sign == SIGN_NEGATIVE)
    {
      BigIntAdd(&z, &modulus, &z);
    }
    BigIntMultiply(&ValNn, &GcdAll, &Tmp[0]);
    do
    {
      SolutionX(&z);
      BigIntAdd(&z, &modulus, &z);
      BigIntSubt(&z, &Tmp[0], &Tmp[1]);
    } while (Tmp[1].sign == SIGN_NEGATIVE);
    return;
  }
  if (callbackQuadModType == CBACK_QMOD_PARABOLIC)
  {    // For elliptic case, the factorization is already done.
    NumberLength = modulus.nbrLimbs;
    CompressBigInteger(nbrToFactor, &modulus);
    Bin2Dec(modulus.limbs, tofactorDec, modulus.nbrLimbs, groupLen);
    factor(&modulus, nbrToFactor, factorsMod, astFactorsMod, NULL);
  }
  intToBigInteger(&Q, 0);
  nbrFactors = astFactorsMod[0].multiplicity;   // Get number of different prime factors.
  pstFactor = &astFactorsMod[1];                // Point to first prime factor information.
  for (factorIndex = 0; factorIndex<nbrFactors; factorIndex++)
  {
    int sol1Invalid = 0;
    int sol2Invalid = 0;
    expon = pstFactor->multiplicity;            // Get exponent multiplicity.
    if (expon == 0)
    {
      pstFactor++;                              // Do not process prime if the exponent is zero.
      continue;
    }
    NumberLength = *pstFactor->ptrFactor;
    UncompressBigInteger(pstFactor->ptrFactor, &prime);   // Get value of prime factor.
    BigIntPowerIntExp(&prime, expon, &V);       // Get value of prime power.
    BigIntRemainder(&coeffQuadr, &V, &L);       // Check whether quadratic coefficient is multiple of prime power.
    if (BigIntIsZero(&L))
    {     // ValA multiple of prime, means linear equation mod prime.
      if (BigIntIsZero(&coeffLinear) && !BigIntIsZero(&coeffIndep))
      {
        return;  // There are no solutions: ValB is equal to zero and ValC is different from zero.
      }
      BigIntGcd(&coeffLinear, &V, &U);
      BigIntDivide(&V, &U, &Q);
      BigIntDivide(&coeffLinear, &U, &V);
      BigIntDivide(&coeffIndep, &U, &L);
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
    {                   // Quadratic equation mod prime
      BigInteger squareRoot;
      BigInteger tmp1, tmp2, ValAOdd;
      int nbrBitsSquareRoot, correctBits, nbrLimbs, bitsAZero, mask, ctr;
      int deltaZeros, deltaIsZero = 0;
      // Compute discriminant = ValB^2 - 4*ValA*ValC.
      BigIntMultiply(&coeffLinear, &coeffLinear, &Tmp[0]);
      BigIntMultiply(&coeffQuadr, &coeffIndep, &discriminant);
      multint(&discriminant, &discriminant, 4);
      BigIntSubt(&Tmp[0], &discriminant, &discriminant);
      CopyBigInt(&ValAOdd, &coeffQuadr);
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
        if (ValAOdd.sign == SIGN_POSITIVE)
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
        if (coeffLinear.sign == SIGN_POSITIVE)
        {
          AddBigNbr((int *)coeffLinear.limbs, (int *)squareRoot.limbs, (int *)tmp1.limbs, nbrLimbs);
        }
        else
        {
          SubtractBigNbr((int *)coeffLinear.limbs, (int *)squareRoot.limbs, (int *)tmp1.limbs, nbrLimbs);
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
        if (coeffLinear.sign == SIGN_NEGATIVE)
        {
          AddBigNbr((int *)coeffLinear.limbs, (int *)squareRoot.limbs, (int *)tmp1.limbs, nbrLimbs);
        }
        else
        {
          SubtractBigNbr((int *)coeffLinear.limbs, (int *)squareRoot.limbs, (int *)tmp1.limbs, nbrLimbs);
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
          for (;;)
          {
            deltaZeros = 0;
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
        BigIntModularDivision(&tmp2, &ValAOdd, &tmp1, &Tmp[0]);
        CopyBigInt(&ValAOdd, &Tmp[0]);
        if (BigIntIsZero(&discriminant))
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
          BigIntRemainder(&discriminant, &prime, &Tmp[3]);
          if (Tmp[3].sign == SIGN_NEGATIVE)
          {
            BigIntAdd(&Tmp[3], &prime, &Tmp[3]);
          }
          if (BigIntJacobiSymbol(&Tmp[3], &prime) != 1)
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
            BigIntModularPower(&Tmp[3], &Q, &SqrtDisc);
          }
          else
          {
            limb *toConvert;
            // Convert discriminant to Montgomery notation.
            CompressLimbsBigInteger(Tmp[5].limbs, &Tmp[3]);
            modmult(Tmp[5].limbs, MontgomeryMultR2, Tmp[6].limbs);  // u
            if ((prime.limbs[0].x & 7) == 5)
            {              // prime mod 8 = 5: use Atkin's method for modular square roots.
               // Step 1. v <- (2u)^((p-5)/8) mod p
               // Step 2. i <- (2uv^2) mod p
               // Step 3. square root of u <- uv (i-1)
              // Step 1.
              subtractdivide(&Q, 5, 8);   // Q <- (prime-5)/8.
              AddBigNbrMod(Tmp[6].limbs, Tmp[6].limbs, Tmp[7].limbs); // 2u
              modPow(Tmp[7].limbs, Q.limbs, Q.nbrLimbs, Tmp[8].limbs);
              // At this moment Aux[7].limbs is v in Montgomery notation.

              // Step 2.
              modmult(Tmp[8].limbs, Tmp[8].limbs, Tmp[9].limbs);  // v^2
              modmult(Tmp[7].limbs, Tmp[9].limbs, Tmp[9].limbs);  // i

              // Step 3.
              SubtBigNbrMod(Tmp[9].limbs, MontgomeryMultR1, Tmp[9].limbs); // i-1
              modmult(Tmp[8].limbs, Tmp[9].limbs, Tmp[9].limbs);           // v*(i-1)
              modmult(Tmp[6].limbs, Tmp[9].limbs, Tmp[9].limbs);           // u*v*(i-1)
              toConvert = Tmp[9].limbs;
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
                intToBigInteger(&Tmp[3], ++x);
              } while (BigIntJacobiSymbol(&Tmp[3], &prime) >= 0);
              // Step 3.
              // Get z <- x^q (mod p) in Montgomery notation.
              modPowBaseInt(x, Q.limbs, Q.nbrLimbs, Tmp[4].limbs);  // z
              // Step 4.
              memcpy(Tmp[5].limbs, Tmp[4].limbs, NumberLength * sizeof(limb)); // y
              r = e;
              CopyBigInt(&K1, &Q);
              subtractdivide(&K1, 1, 2);
              modPow(Tmp[6].limbs, K1.limbs, K1.nbrLimbs, Tmp[7].limbs); // x
              modmult(Tmp[6].limbs, Tmp[7].limbs, Tmp[8].limbs);         // v
              modmult(Tmp[8].limbs, Tmp[7].limbs, Tmp[9].limbs);         // w
              // Step 5
              while (memcmp(Tmp[9].limbs, MontgomeryMultR1, NumberLength * sizeof(limb)))
              {
                // Step 6
                k = 0;
                memcpy(Tmp[10].limbs, Tmp[9].limbs, NumberLength * sizeof(limb));
                do
                {
                  k++;
                  modmult(Tmp[10].limbs, Tmp[10].limbs, Tmp[10].limbs);
                } while (memcmp(Tmp[10].limbs, MontgomeryMultR1, NumberLength * sizeof(limb)));
                // Step 7
                memcpy(Tmp[11].limbs, Tmp[5].limbs, NumberLength * sizeof(limb)); // d
                for (ctr = 0; ctr < r - k - 1; ctr++)
                {
                  modmult(Tmp[11].limbs, Tmp[11].limbs, Tmp[11].limbs);
                }
                modmult(Tmp[11].limbs, Tmp[11].limbs, Tmp[5].limbs);   // y
                r = k;
                modmult(Tmp[8].limbs, Tmp[11].limbs, Tmp[8].limbs);    // v
                modmult(Tmp[9].limbs, Tmp[5].limbs, Tmp[9].limbs);     // w
              }
              toConvert = Tmp[8].limbs;
            }
            // Convert from Montgomery to standard notation.
            memset(Tmp[4].limbs, 0, NumberLength * sizeof(limb)); // Convert power to standard notation.
            Tmp[4].limbs[0].x = 1;
            modmult(Tmp[4].limbs, toConvert, toConvert);
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
        BigIntAdd(&coeffLinear, &squareRoot, &tmp1);
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
        BigIntSubt(&coeffLinear, &squareRoot, &tmp1);
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
      BigIntSubt(&Solution2[factorIndex], &Solution1[factorIndex], &Tmp[0]);
      if (Tmp[0].sign == SIGN_NEGATIVE)
      {     // Solution2 is less than Solution1, so exchange them.
        CopyBigInt(&Tmp[0], &Solution1[factorIndex]);
        CopyBigInt(&Solution1[factorIndex], &Solution2[factorIndex]);
        CopyBigInt(&Solution2[factorIndex], &Tmp[0]);
      }
    }
    CopyBigInt(&Increment[factorIndex], &Q);
    Exponents[factorIndex] = 0;
    pstFactor++;
  }   // end for
  // Use Chinese remainder theorem to obtain the solutions.
  do
  {
    nbrFactors = astFactorsMod[0].multiplicity;
    pstFactor = &astFactorsMod[1];
    for (T1 = 0; T1 < nbrFactors; T1++)
    {
      if (pstFactor->multiplicity)
      {
        break;   // Do not process prime factors for which the multiplicity is zero.
      }
      pstFactor++;
    }
    multint(&Aux[T1], &Increment[T1], Exponents[T1] >> 1);
    if (Exponents[T1] & 1)
    {
      BigIntAdd(&Aux[T1], &Solution2[T1], &Aux[T1]);
    }
    else
    {
      BigIntAdd(&Aux[T1], &Solution1[T1], &Aux[T1]);
    }
    CopyBigInt(&currentSolution, &Aux[T1]);
    UncompressBigInteger(pstFactor->ptrFactor, &prime);
    BigIntPowerIntExp(&prime, pstFactor->multiplicity, &Mult);
    for (T1++; T1<nbrFactors; T1++)
    {
      pstFactor++;
      if (pstFactor->multiplicity == 0)
      {
        continue;   // Do not process prime factors for which the multiplicity is zero.
      }
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
        if (astFactorsMod[E + 1].multiplicity == 0)
        {
          continue;
        }
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
      SolutionX(&K1);
      addbigint(&V, 1);  // V <- V + 1
    }
    for (T1 = nbrFactors - 1; T1 >= 0; T1--)
    {
      if (astFactorsMod[T1 + 1].multiplicity == 0)
      {
        continue;     // Do not process prime factors for which the multiplicity is zero.
      }
      UncompressBigInteger(astFactorsMod[T1+1].ptrFactor, &K);         // Get prime factor.
      BigIntPowerIntExp(&K, astFactorsMod[T1+1].multiplicity, &prime); // Get prime power.
      if (BigIntEqual(&Solution1[T1], &Solution2[T1]))
      {     // Solution1[T1] == Solution2[T1]
        Exponents[T1] += 2;
      }
      else
      {     // Solution1[T1] != Solution2[T1]
        Exponents[T1]++;
      }
      multadd(&L, Exponents[T1], &Increment[T1], 0);   // L <- Exponents[T1] * Increment[T1]
      multadd(&K1, 2, &prime, 0);                      // K1 <- 2 * prime power
      BigIntSubt(&L, &K1, &L);                         // L <- Exponents[T1] * Increment[T1] - 2 * prime power.
      if (L.sign == SIGN_NEGATIVE)
      {
        break;
      }
      Exponents[T1] = 0;
    }   /* end for */
  } while (T1 >= 0);
}

static void paren(BigInteger *num)
{
  if (num->sign == SIGN_NEGATIVE)
  {
    *ptrOutput++ = '(';
    shownbr(num);
    *ptrOutput++ = ')';
    *ptrOutput = 0;
  }
  else
  {
    BigInteger2Dec(num, ptrOutput, groupLen);
    ptrOutput += strlen(ptrOutput);
  }
}

enum eLinearSolution Linear(BigInteger *coeffX, BigInteger *coeffY, BigInteger *coeffInd)
{
  BigInteger q;
  int showSteps, stepNbr;
  if (teach)
  {
    w("<p>This is a linear equation ");
    ShowLin(coeffX, coeffY, coeffInd, "x", "y");
    w(" = 0</p>");
  }
  if (BigIntIsZero(coeffX))
  {
    if (BigIntIsZero(coeffY))
    {
      if (!BigIntIsZero(coeffInd))
      {
        return NO_SOLUTIONS;           // No solutions
      }
      else
      {
        return INFINITE_SOLUTIONS;     // Infinite number of solutions
      }
    }
    BigIntRemainder(coeffInd, coeffY, &Aux[0]);
    if (!BigIntIsZero(&Aux[0]))
    {
      return NO_SOLUTIONS;             // No solutions
    }
    else
    {
      intToBigInteger(&Xind, 0);
      intToBigInteger(&Xlin, 1);
      BigIntDivide(coeffInd, coeffY, &Yind);
      BigIntNegate(&Yind, &Yind);
      intToBigInteger(&Ylin, 0);
      return SOLUTION_FOUND;           // Solution found
    }
  }
  if (BigIntIsZero(coeffY))
  {
    BigIntRemainder(coeffInd, coeffX, &Aux[0]);
    if (!BigIntIsZero(&Aux[0]))
    {
      return NO_SOLUTIONS;             // No solutions
    }
    else
    {
      intToBigInteger(&Yind, 0);
      intToBigInteger(&Ylin, 1);
      BigIntDivide(coeffInd, coeffX, &Xind);
      BigIntNegate(&Xind, &Xind);
      intToBigInteger(&Xlin, 0);
      return SOLUTION_FOUND;           // Solution found
    }
  }
  BigIntGcd(coeffX, coeffY, &U1);
  if (U1.nbrLimbs > 1 || U1.limbs[0].x != 1)
  {                  // GCD is not 1.
    if (teach)
    {
      w("<p>To solve it, we first find the gcd of the linear coefficients, that is: gcd(");
      shownbr(coeffX);
      w(", ");
      shownbr(coeffY);
      w(") = ");
      shownbr(&U1);
      w(".</p>");
    }
    BigIntRemainder(coeffInd, &U1, &U2);
    if (!BigIntIsZero(&U2))
    {
      if (teach)
      {
        w("<p>The independent coefficient is not multiple of gcd, so there are no solutions.</p>");
      }
      return NO_SOLUTIONS;            // No solutions
    }
    BigIntDivide(coeffX, &U1, coeffX);   // Divide all coefficients by the gcd.
    BigIntDivide(coeffY, &U1, coeffY);
    BigIntDivide(coeffInd, &U1, coeffInd);
  }
  if (teach)
  {
    w(divgcd);
    ShowLin(coeffX, coeffY, coeffInd, "x", "y");
    w(" = 0</p>");
    w("<p>Now we must apply the Generalized Euclidean algorithm:</p>");
  }
  intToBigInteger(&U1, 1);    // U1 <- 1
  intToBigInteger(&U2, 0);    // U2 <- 0
  CopyBigInt(&U3, coeffX);    // U3 <- coeffX
  intToBigInteger(&V1, 0);    // V1 <- 0
  intToBigInteger(&V2, 1);    // V2 <- 1
  CopyBigInt(&V3, coeffY);    // V3 <- coeffY
  showSteps = (coeffX->nbrLimbs == 1 && coeffY->nbrLimbs == 1);
  stepNbr = 1;
  if (teach && showSteps)
  {
    w("<p>");
  }
  for (;;)
  {
    if (BigIntIsZero(&V3))
    {                    // V3 equals zero, so go out of loop.
      break;
    }
    if (teach && showSteps)
    {
      w("Step ");
      showInt(stepNbr);
      w(": ");
      paren(&U1);    // U1
      w(" &times; ");
      paren(coeffX);
      w(" + ");
      paren(&U2);    // U2
      w(" &times; ");
      paren(coeffY);
      w(" = ");
      paren(&U3);    // U3
      w("<br>");
    }
    floordiv(&U3, &V3, &q);                 // q <- floor(U3 / V3).
    BigIntMultiply(&q, &V1, &bigTmp);
    BigIntSubt(&U1, &bigTmp, &bigTmp);      // T <- U1 - q * V1
    CopyBigInt(&U1, &V1);                   // U1 <- V1
    CopyBigInt(&V1, &bigTmp);               // V1 <- T
    BigIntMultiply(&q, &V2, &bigTmp);
    BigIntSubt(&U2, &bigTmp, &bigTmp);      // T <- U2 - q * V2
    CopyBigInt(&U2, &V2);                   // U2 <- V2
    CopyBigInt(&V2, &bigTmp);               // V2 <- T
    BigIntMultiply(&q, &V3, &bigTmp);
    BigIntSubt(&U3, &bigTmp, &bigTmp);      // T <- U3 - q * V3
    CopyBigInt(&U3, &V3);                   // U3 <- V3
    CopyBigInt(&V3, &bigTmp);               // V3 <- T
    stepNbr++;
  }
  BigIntDivide(coeffInd, &U3, &q);
  BigIntChSign(&q);                         // q <- -coeffInd / U3
  // Xind <- -U1 * coeffInd / U3
  BigIntMultiply(&U1, &q, &Xind);
  // Xlin <- coeffY
  CopyBigInt(&Xlin, coeffY);
  // Yind <- -U2 * coeffInd / U3
  BigIntMultiply(&U2, &q, &Yind);
  // Ylin <- -coeffX
  CopyBigInt(&Ylin, coeffX);
  BigIntChSign(&Ylin);
  if (teach)
  {
    w("Step ");
    showInt(stepNbr);
    w(": ");
    paren(&U1);    // U1
    w(" &times; ");
    paren(coeffX);
    w(" + ");
    paren(&U2);    // U2
    w(" &times; ");
    paren(coeffY);
    w(" = ");
    paren(&U3);    // U3
    CopyBigInt(&bigTmp, coeffInd);
    BigIntChSign(&bigTmp);
    if (q.sign != SIGN_POSITIVE || q.nbrLimbs != 1 || q.limbs[0].x != 1)
    {    // Multiplier is not 1.
      w("</p><p>Multiplying the last equation by ");
      paren(&q);
      w(" we obtain:<br>");
      paren(&Xind);
      w(" &times; ");
      paren(coeffX);
      w(" + ");
      paren(&Yind);
      w(" &times; ");
      paren(coeffY);
      w(" = ");
      shownbr(&bigTmp);
    }
    w("</p><p>Adding and subtracting ");
    paren(coeffX);
    w(" &times; ");
    paren(coeffY);
    w(" t' we obtain:</p><p>(");
    shownbr(&Xind);
    w(" + ");
    paren(coeffY);
    w(" t') &times; ");
    paren(coeffX);
    w(" + (");
    shownbr(&Yind);
    w(" - ");
    paren(coeffX);
    w(" t') &times; ");
    paren(coeffY);
    w(" = ");
    shownbr(&bigTmp);
    w("</p><p>So, the solution is given by the set:</p>");
    PrintLinear(0, "t'");
    w("</p>");
  }

  // Substitute variables so the independent coefficients can be minimized.
  // Reuse variables U1, U2, U3, V1, V2, V3.
  BigIntMultiply(coeffX, coeffX, &U1);
  BigIntMultiply(coeffY, coeffY, &q);
  BigIntAdd(&U1, &q, &U1);                  // U1 <- coeffX^2 + coeffY^2
  CopyBigInt(&U2, &U1);
  subtractdivide(&U2, 0, 2);                // U2 <- (coeffX^2 + coeffY^2)/2
  BigIntMultiply(coeffX, &Yind, &q);
  BigIntAdd(&U2, &q, &U2);
  BigIntMultiply(coeffY, &Xind, &q);
  BigIntSubt(&U2, &q, &U2);
  floordiv(&U2, &U1, &U1);                  // U1 <- delta to add to t'
  if (teach)
  {
    w("<p>By making the substitution t = ");
    shownbr(&U1);
    w(" + t' we finally obtain:</p>");
  }
  BigIntMultiply(&U1, coeffY, &q);
  BigIntAdd(&Xind, &q, &Xind);        // Xind <- Xind + coeffY * delta
  BigIntMultiply(&U1, coeffX, &q);
  BigIntSubt(&Yind, &q, &Yind);       // Yind <- Yind - coeffX * delta
  if (Xlin.sign == SIGN_NEGATIVE && Ylin.sign == SIGN_NEGATIVE)
  {    // If both coefficients are negative, make them positive.
    BigIntChSign(&Xlin);
    BigIntChSign(&Ylin);
  }
  return SOLUTION_FOUND;
}

static void EnsureCoeffAIsNotZero(void)
{
  if (BigIntIsZero(&ValA))
  {      // Next algorithm does not work if A = 0. In this case, exchange x and y.
    ExchXY = 1;
    CopyBigInt(&bigTmp, &ValA);   // Exchange coefficients of x^2 and y^2.
    CopyBigInt(&ValA, &ValC);
    CopyBigInt(&ValC, &bigTmp);
    CopyBigInt(&bigTmp, &ValD);   // Exchange coefficients of x and y.
    CopyBigInt(&ValD, &ValE);
    CopyBigInt(&ValE, &bigTmp);
  }
}

static void DiscriminantIsZero(void)
{
  EnsureCoeffAIsNotZero();
  // ax^2 + bxy + cx^2 + dx + ey + f = 0 (1)
  // Multiplying by 4a:
  // (2ax + by)^2 + 4adx + 4aey + 4af = 0
  // Let t = 2ax + by. So (1) becomes: (t + d)^2 = uy + v.
  // Compute u <- 2(bd - 2ae)
  BigIntMultiply(&ValB, &ValD, &Aux[0]);
  BigIntMultiply(&ValA, &ValE, &Aux[1]);
  multint(&Aux[2], &Aux[1], 2);
  BigIntSubt(&Aux[0], &Aux[2], &Aux[1]);
  multint(&ValU, &Aux[1], 2);

  // Compute v <- d^2 - 4af
  BigIntMultiply(&ValD, &ValD, &Aux[2]);
  BigIntMultiply(&ValA, &ValF, &Aux[1]);
  multint(&Aux[3], &Aux[1], 4);
  BigIntSubt(&Aux[2], &Aux[3], &ValV);

  if (BigIntIsZero(&ValU))
  { // u equals zero, so (t+d)^2 = v.
    enum eLinearSolution ret;
    if (ValV.sign == SIGN_NEGATIVE)
    { // There are no solutions when v is negative.
      return;
    }
    if (BigIntIsZero(&ValV))
    { // v equals zero, so (1) becomes 2ax + by + d = 0
      multint(&Aux[3], &ValA, 2);
      ret = Linear(&Aux[3], &ValB, &ValD);
      if (teach)
      {
        w("<div class=\"result\">");
      }
      PrintLinear(ret, "t");
      if (teach)
      {
        w("</span>");
      }
      return;
    }
    // u equals zero but v does not.
    // v must be a perfect square, otherwise there are no solutions.
    squareRoot(ValV.limbs, ValG.limbs, ValV.nbrLimbs, &ValG.nbrLimbs);
    ValG.sign = SIGN_POSITIVE;          // g <- sqrt(v).
    BigIntMultiply(&ValG, &ValG, &Aux[3]);
    if (!BigIntEqual(&ValV, &Aux[3]))
    {   // v is not perfect square.
      return;
    }
    // The original equation is now: 2ax + by + (d +/- g) = 0
    multint(&Aux[3], &ValA, 2);
    BigIntAdd(&ValD, &ValG, &Aux[4]);
    ret = Linear(&Aux[3], &ValB, &Aux[4]);
    if (teach)
    {
      w("<div class=\"result\">");
    }
    PrintLinear(ret, "t");
    if (teach)
    {
      w("</span>");
    }
    multint(&Aux[3], &ValA, 2);
    BigIntSubt(&ValD, &ValG, &Aux[4]);
    ret = Linear(&Aux[3], &ValB, &Aux[4]);
    if (teach)
    {
      w("<div class=\"result\">");
    }
    PrintLinear(ret, "t");
    if (teach)
    {
      w("</span>");
    }
    return;
  }
  // At this moment u does not equal zero.
  // We have to solve the congruence T^2 = v (mod u) where T = t+d and t = 2ax+by.
  callbackQuadModType = CBACK_QMOD_PARABOLIC;
  intToBigInteger(&coeffQuadr, 1);
  intToBigInteger(&coeffLinear, 0);
  CopyBigInt(&coeffIndep, &ValV);
  BigIntChSign(&coeffIndep);
  CopyBigInt(&modulus, &ValU);
  SolveQuadModEquation();
}

static void callbackQuadModParabolic(BigInteger *value)
{  // The argument of this function is T. t = T - d + uk (k arbitrary).
   // Compute ValR <- (T^2 - v)/u
  BigIntMultiply(value, value, &bigTmp);
  BigIntSubt(&bigTmp, &ValV, &bigTmp);
  BigIntDivide(&bigTmp, &ValU, &ValR);
   // Compute ValS <- 2*T
  BigIntAdd(value, value, &ValS);
   // Find k from the congruence jk = K (mod z) where j = u-bs, K = d+br-T, z = 2a
   // Compute j <- u-bs
  BigIntMultiply(&ValB, &ValS, &bigTmp);
  BigIntSubt(&ValU, &bigTmp, &ValJ);
   // Compute K <- d+br-T
  BigIntMultiply(&ValB, &ValR, &bigTmp);
  BigIntAdd(&ValD, &bigTmp, &bigTmp);
  BigIntSubt(&bigTmp, value, &ValK);
   // Compute z <- 2a
  BigIntAdd(&ValA, &ValA, &ValZ);
  // If K is not multiple of gcd(g, z) there is no solution.
  BigIntGcd(&ValJ, &ValZ, &bigTmp);
  BigIntRemainder(&ValK, &bigTmp, &U1);
  if (!BigIntIsZero(&U1))
  {
    return;
  }
  // Compute g = gcd(j, K, z), then recalculate j <- j/g, K <- K/g, z <- z/g
  BigIntGcd(&bigTmp, &ValK, &U1);
  BigIntDivide(&ValJ, &U1, &U2);    // U2 <- j
  BigIntDivide(&ValK, &U1, &U3);    // U3 <- K
  BigIntDivide(&ValZ, &U1, &ValZ);
  ValZ.sign = SIGN_POSITIVE;        // Use positive sign for modulus.
  BigIntRemainder(&U2, &ValZ, &U2);    
  if (U2.sign == SIGN_NEGATIVE)
  {
    BigIntAdd(&U2, &ValZ, &U2);
  }
  BigIntRemainder(&U3, &ValZ, &U3);
  if (U3.sign == SIGN_NEGATIVE)
  {
    BigIntAdd(&U3, &ValZ, &U3);
  }
  if (BigIntIsZero(&U2))
  {    // M and N equal zero.
    intToBigInteger(&ValZ, 1);     // In this case 0*k = 0 (mod z) means any k is valid.
  }
  else
  {
    BigIntGeneralModularDivision(&U2, &U3, &ValZ, &U2);   // U2 <- x'
  }
   // Compute coefficients of x: V1 + V2 * w + V3 * w^2
   // Let m = 2be - 4cd
  BigIntMultiply(&ValB, &ValE, &U3);
  BigIntMultiply(&ValC, &ValD, &bigTmp);
  BigIntAdd(&bigTmp, &bigTmp, &bigTmp);
  BigIntSubt(&U3, &bigTmp, &U3);
  BigIntAdd(&U3, &U3, &U3);        // U3 <- m
   // Compute V1 <- (x'j - k)/2a + mx'^2
  BigIntMultiply(&U2, &ValJ, &bigTmp);
  BigIntSubt(&bigTmp, &ValK, &bigTmp);
  BigIntDivide(&bigTmp, &ValA, &V1);
  subtractdivide(&V1, 0, 2);
  BigIntMultiply(&U3, &U2, &bigTmp);
  BigIntMultiply(&bigTmp, &U2, &bigTmp);
  BigIntAdd(&V1, &bigTmp, &V1);
   // Compute V2 <- (j/2a + 2mx')z
  BigIntDivide(&ValJ, &ValA, &V2);
  subtractdivide(&V2, 0, 2);
  BigIntMultiply(&U3, &U2, &bigTmp);
  BigIntAdd(&bigTmp, &bigTmp, &bigTmp);
  BigIntAdd(&V2, &bigTmp, &V2);
  BigIntMultiply(&V2, &ValZ, &V2);
   // Compute V3 <- m*z^2
  BigIntMultiply(&U3, &ValZ, &V3);
  BigIntMultiply(&V3, &ValZ, &V3);
  showAlso();
  w("x = ");
  PrintQuad(&V3, &V2, &V1);
   // Compute coefficients of y: V1 + V2 * w + V3 * w^2
   // Compute V1 <- r + sx' + ux'^2
  BigIntMultiply(&ValU, &U2, &V1);
  BigIntAdd(&V1, &ValS, &V1);
  BigIntMultiply(&V1, &U2, &V1);
  BigIntAdd(&V1, &ValR, &V1);
   // Compute V2 <- (s + 2ux')z
  BigIntMultiply(&ValU, &U2, &V2);
  BigIntAdd(&V2, &V2, &V2);
  BigIntAdd(&V2, &ValS, &V2);
  BigIntMultiply(&V2, &ValZ, &V2);
   // Compute V3 <- uz^2
  BigIntMultiply(&ValU, &ValZ, &V3);
  BigIntMultiply(&V3, &ValZ, &V3);
  w("y = ");
  PrintQuad(&V3, &V2, &V1);
}

static void ShowPoint(BigInteger *X, BigInteger *Y)
{
  // Check first that (X+alpha) and (Y+beta) are multiple of D.
  BigIntAdd(X, &ValAlpha, &Tmp1);
  BigIntRemainder(&Tmp1, &ValDiv, &Tmp2);
  if (!BigIntIsZero(&Tmp2))
  {
    return;
  }
  BigIntAdd(Y, &ValBeta, &Tmp2);
  BigIntRemainder(&Tmp2, &ValDiv, &Tmp3);
  if (!BigIntIsZero(&Tmp3))
  {
    return;
  }
  BigIntDivide(&Tmp1, &ValDiv, &Tmp1);
  BigIntDivide(&Tmp2, &ValDiv, &Tmp2);
  ShowXY(&Tmp1, &Tmp2);
  showRecursiveSolution = 1; // Show recursive solution if it exists.
}

// Solve ax^2+bxy+cy^2 = K
// The quadratic modular equation algorithm requires that gcd(a, n) = 1.
// At this point gcd(a, b, c) = 1
// The possibilities are:
// - gcd(a, K) = 1. There is nothing to do.
// Otherwise perform the transformation x = PX + QY, y = RX + SY with PS - QR = 1.
// In particular perform: x = mX + (m-1)Y, y = X + Y
// We get: (am^2+bm+c)*X^2 + (2(am^2+bm+c) - (2am+b))*XY + ((am^2+bm+c) - (2am+b) + a)*Y^2
// Also perform: x = X + Y, y = (m-1)X + mY
// We get: (a+(m-1)*b+(m-1)^2*c)*X^2 + (2a + b(2m-1) + 2cm(m-1))*X*Y + (a+bm+cm^2)*Y^2
// The discriminant of the new formula does not change.
// Compute m=1, 2, 3,... until gcd(am^2+bm+c, K) = 1.
// When using the second formula, change sign of m so we know the formula used when
// undoing the unimodular transformation later.

// Since the algorithm discovers only primitive solutions, i.e. solutions (x,y) where
// gcd(x,y) = 1, we need to solve ax'^2+bx'y'+cy'^2 = K/R^2 where R^2 is a divisor of K.
// Then we get x = Rx', y = Ry'.
static void NonSquareDiscriminant(void)
{
  int factorNbr, numFactors;
  struct sFactors *pstFactor;
             // Find GCD(a,b,c)
  BigIntGcd(&ValA, &ValB, &bigTmp);
  BigIntGcd(&ValC, &bigTmp, &ValGcdHomog);
  CopyBigInt(&ValDiv, &discr);
    // Divide A, B, C and K by this GCD.
  BigIntDivide(&ValA, &ValGcdHomog, &ValA);
  BigIntDivide(&ValB, &ValGcdHomog, &ValB);
  BigIntDivide(&ValC, &ValGcdHomog, &ValC);
  BigIntDivide(&ValK, &ValGcdHomog, &ValK);
  // Divide discriminant by the square of GCD.
  BigIntDivide(&discr, &ValGcdHomog, &discr);
  BigIntDivide(&discr, &ValGcdHomog, &discr);
  if (BigIntIsZero(&ValK))
  {          // If k=0, the only solution is (X, Y) = (0, 0)
    ShowPoint(&ValK, &ValK);
    return;
  }
  // Factor independent term.
  NumberLength = ValK.nbrLimbs;
  CompressBigInteger(nbrToFactor, &ValK);
  Bin2Dec(ValK.limbs, tofactorDec, ValK.nbrLimbs, groupLen);
  factor(&ValK, nbrToFactor, factorsMod, astFactorsMod, NULL);
  // Find all indexes of prime factors with even multiplicity.
  nbrPrimesEvenMultiplicity = 0;
  numFactors = astFactorsMod[0].multiplicity;
  pstFactor = &astFactorsMod[1];
  for (factorNbr = 1; factorNbr <= numFactors; factorNbr++)
  {
    if (pstFactor->multiplicity > 1)
    {      // At least prime is squared.
      indexEvenMultiplicity[nbrPrimesEvenMultiplicity] = factorNbr;
      originalMultiplicities[nbrPrimesEvenMultiplicity++] = pstFactor->multiplicity & (~1); // Convert to even.
    }
    pstFactor++;
  }
  memset(counters, 0, sizeof(counters));
  memset(isDescending, 0, sizeof(isDescending));
  intToBigInteger(&ValE, 1);  // Initialize multiplier to 1.
  // Loop that cycles through all square divisors of the independent term.
  CopyBigInt(&ValAbak, &ValA);
  CopyBigInt(&ValBbak, &ValB);
  CopyBigInt(&ValCbak, &ValC);
  for (;;)
  {
    int index;
    CopyBigInt(&ValA, &ValAbak);
    CopyBigInt(&ValB, &ValBbak);
    CopyBigInt(&ValC, &ValCbak);
    BigIntGcd(&ValA, &ValK, &bigTmp);
    intToBigInteger(&ValM, 0);
    if (bigTmp.nbrLimbs != 1 || bigTmp.limbs[0].x != 1)
    {                 // gcd(a, K) is not equal to 1.
      intToBigInteger(&ValM, 0);
      do
      {                          // Compute cm^2 + bm + a and exit loop if this value is not coprime with K.
        BigIntMultiply(&ValC, &ValM, &U2);
        BigIntAdd(&U2, &ValB, &U1);
        BigIntMultiply(&U1, &ValM, &U1);
        BigIntAdd(&U1, &ValA, &U1);
        BigIntGcd(&U1, &ValK, &bigTmp);
        if (bigTmp.nbrLimbs == 1 && bigTmp.limbs[0].x == 1)
        {
          addbigint(&ValM, 1);  // Increment M.
          BigIntChSign(&ValM);  // Change sign to indicate type.
          break;
        }
        addbigint(&ValM, 1);    // Increment M.
                                // Compute am^2 + bm + c and loop while this value is not coprime with K.
        BigIntMultiply(&ValA, &ValM, &U2);
        BigIntAdd(&U2, &ValB, &U1);
        BigIntMultiply(&U1, &ValM, &U1);
        BigIntAdd(&U1, &ValC, &U1);
        BigIntGcd(&U1, &ValK, &bigTmp);
      } while (bigTmp.nbrLimbs != 1 || bigTmp.limbs[0].x != 1);
      // Compute 2am + b or 2cm + b as required.
      BigIntAdd(&U2, &U2, &U2);
      BigIntAdd(&U2, &ValB, &U2);
      if (ValM.sign == SIGN_POSITIVE)
      {
        // Compute c.
        BigIntSubt(&U1, &U2, &ValB);
        BigIntAdd(&ValB, &ValA, &ValC);
        // Compute b.
        BigIntAdd(&ValB, &U1, &ValB);
        // Compute a.
        CopyBigInt(&ValA, &U1);
      }
      else
      {
        // Compute a.
        BigIntSubt(&U1, &U2, &ValB);
        BigIntAdd(&ValB, &ValC, &ValA);
        // Compute b.
        BigIntAdd(&ValB, &U1, &ValB);
        // Compute c.
        CopyBigInt(&ValC, &U1);
      }
    }
    CopyBigInt(&coeffQuadr, &ValA);
    CopyBigInt(&coeffLinear, &ValB);
    CopyBigInt(&coeffIndep, &ValC);
    CopyBigInt(&modulus, &ValK);
    SolveQuadModEquation();
    // Adjust counters.
    for (index = 0; index < nbrPrimesEvenMultiplicity; index++)
    {            // Loop that increments counters.
      if (isDescending[index] == 0)
      {          // Ascending.
        pstFactor = &astFactorsMod[indexEvenMultiplicity[index]];
        if (counters[index] == originalMultiplicities[index])
        {
          isDescending[index] = 1;    // Next time it will be descending.
          continue;
        }
        else
        {
          NumberLength = *pstFactor->ptrFactor;
          UncompressBigInteger(pstFactor->ptrFactor, &bigTmp);
          BigIntMultiply(&bigTmp, &bigTmp, &U3);
          pstFactor->multiplicity -= 2;
          BigIntDivide(&ValK, &U3, &ValK);       // Divide by square of prime.
          BigIntMultiply(&ValE, &bigTmp, &ValE); // Multiply multiplier by prime.counters[index]++;
          counters[index] += 2;
          break;
        }
      }
      else
      {         // Descending.
        if (counters[index] <= 1)
        {
          isDescending[index] = 0;    // Next time it will be ascending.
          continue;
        }
        else
        {
          pstFactor = &astFactorsMod[indexEvenMultiplicity[index]];
          NumberLength = *pstFactor->ptrFactor;
          UncompressBigInteger(pstFactor->ptrFactor, &bigTmp);
          BigIntMultiply(&bigTmp, &bigTmp, &U3);
          pstFactor->multiplicity += 2;
          BigIntMultiply(&ValK, &U3, &ValK);     // Multiply by square of prime.
          BigIntDivide(&ValE, &bigTmp, &ValE);   // Divide multiplier by prime.counters[index]++;
          counters[index] -= 2;
          break;
        }
      }
    }
    if (index == nbrPrimesEvenMultiplicity)
    {               // All factors have been found. Exit loop.
      break;
    }
  }
  if (showRecursiveSolution && callbackQuadModType == CBACK_QMOD_HYPERBOLIC)
  {   // Show recursive solution.
    CopyBigInt(&ValA, &ValAbak);
    CopyBigInt(&ValB, &ValBbak);
    CopyBigInt(&ValC, &ValCbak);
    ContFracPell();
  }
}

static void NegativeDiscriminant(void)
{
  callbackQuadModType = CBACK_QMOD_ELLIPTIC;
  NonSquareDiscriminant();
}

static void UnimodularSubstitution(void)
{
  if (ValM.sign == SIGN_NEGATIVE)
  {     // Perform the substitution: x = X + Y, y = (|m|-1)X + |m|Y
    ValM.sign = SIGN_POSITIVE;
    BigIntAdd(&ValZ, &ValO, &Tmp[1]);     // x
    BigIntMultiply(&Tmp[1], &ValM, &Tmp[0]);
    BigIntSubt(&Tmp[0], &ValZ, &Tmp[0]);  // y
    ShowPoint(&Tmp[1], &Tmp[0]);
  }
  else if (BigIntIsZero(&ValM))
  {
    ShowPoint(&ValZ, &ValO);
  }
  else
  {     // Perform the substitution: x = mX + (m-1)Y, y = X + Y
    BigIntAdd(&ValZ, &ValO, &Tmp[1]);     // y
    BigIntMultiply(&Tmp[1], &ValM, &Tmp[0]);
    BigIntSubt(&Tmp[0], &ValO, &Tmp[0]);  // x
    ShowPoint(&Tmp[0], &Tmp[1]);
  }
}

// On input: ValH: value of u, ValI: value of v.
// Output: ((tu - nv)*E, u*E) and ((-tu + nv)*E, -u*E)
// If m is greater than zero, perform the substitution: x = mX + (m-1)Y, y = X + Y
// If m is less than zero, perform the substitution: x = X + Y, y = (|m|-1)X + |m|Y
// Do not substitute if m equals zero.

static void NonSquareDiscrSolution(BigInteger *value)
{
  // Get value of tu - Kv
  BigIntMultiply(value, &ValH, &ValZ);    // tu
  CopyBigInt(&bigTmp, &ValK);
  if (bigTmp.sign == SIGN_NEGATIVE)
  {
    BigIntChSign(&bigTmp);                // Get |K|
  }
  BigIntMultiply(&bigTmp, &ValI, &bigTmp);// |K|v
  BigIntSubt(&ValZ, &bigTmp, &ValZ);      // tu + |K|v
  BigIntMultiply(&ValZ, &ValE, &ValZ);    // X = (tu - |K|v)*E
  BigIntMultiply(&ValH, &ValE, &ValO);    // Y = u*E
  Xbak = &Xplus;
  Ybak = &Yplus;
  UnimodularSubstitution();               // Undo unimodular substitution
  BigIntChSign(&ValZ);                    // (-tu - |K|v)*E
  BigIntChSign(&ValO);                    // -u*E
  Xbak = &Xminus;
  Ybak = &Yminus;
  UnimodularSubstitution();               // Undo unimodular substitution
}

// Obtain next convergent of continued fraction of ValU/ValV
// Previous convergents U1/V1, U2/V2, U3/V3.
static void getNextConvergent(void)
{
  floordiv(&ValU, &ValV, &bigTmp);
  // Values of U3 and V3 are not used, so they can be overwritten now.
  // Compute new value of ValU and ValV.
  BigIntMultiply(&bigTmp, &ValV, &U3);
  BigIntSubt(&ValU, &U3, &U3);
  CopyBigInt(&ValU, &ValV);
  CopyBigInt(&ValV, &U3);
  // Compute new convergents: h_n = a_n*h_{n-1} + h_{n-2}, k_n = k_n*k_{n-1} + k_{n-2}
  BigIntMultiply(&bigTmp, &U1, &V3);
  BigIntAdd(&V3, &U2, &V3);
  CopyBigInt(&U3, &U2);
  CopyBigInt(&U2, &U1);
  CopyBigInt(&U1, &V3);
  BigIntMultiply(&bigTmp, &V1, &bigTmp);
  BigIntAdd(&bigTmp, &V2, &bigTmp);
  CopyBigInt(&V3, &V2);
  CopyBigInt(&V2, &V1);
  CopyBigInt(&V1, &bigTmp);
}

static void callbackQuadModElliptic(BigInteger *value)
{
  // Compute P <- (at^2+bt+c)/K
  BigIntMultiply(&ValA, value, &ValQ);
  BigIntAdd(&ValQ, &ValB, &ValP);
  BigIntMultiply(&ValP, value, &ValP);
  BigIntAdd(&ValP, &ValC, &ValP);
  BigIntDivide(&ValP, &ValK, &ValP);
  // Compute Q <- -(2at + b).
  BigIntAdd(&ValQ, &ValQ, &ValQ);
  BigIntAdd(&ValQ, &ValB, &ValQ);
  BigIntChSign(&ValQ);
  // Compute R <- aK
  BigIntMultiply(&ValA, &ValK, &ValR);
  if (ValP.sign == SIGN_POSITIVE && ValP.nbrLimbs == 1)
  {
    int Plow = ValP.limbs[0].x;
    if ((discr.nbrLimbs > 1 || discr.limbs[0].x > 4) && Plow == 1)
    {      // Discriminant is less than -4 and P equals 1.
      intToBigInteger(&ValH, 1);
      intToBigInteger(&ValI, 0);
      NonSquareDiscrSolution(value);   // (1, 0)
      return;
    }
    if (discr.nbrLimbs == 1 && discr.limbs[0].x == 4)
    {      // Discriminant is equal to -4.
      CopyBigInt(&ValG, &ValQ);
      subtractdivide(&ValG, 0, 2);
      if (Plow == 1)
      {
        intToBigInteger(&ValH, 1);
        intToBigInteger(&ValI, 0);
        NonSquareDiscrSolution(value);   // (1, 0)
        CopyBigInt(&ValH, &ValG);
        intToBigInteger(&ValI, -1);
        NonSquareDiscrSolution(value);   // (G, -1)
        return;
      }
      if (Plow == 2)
      {
        intToBigInteger(&ValI, -1);
        CopyBigInt(&ValH, &ValG);
        subtractdivide(&ValH, 1, 2);
        NonSquareDiscrSolution(value);   // ((G-1)/2, -1)
        intToBigInteger(&ValI, -1);
        CopyBigInt(&ValH, &ValG);
        subtractdivide(&ValH, -1, 2);
        NonSquareDiscrSolution(value);   // ((G+1)/2, -1)
        return;       
      }
    }
    if (discr.nbrLimbs == 1 && discr.limbs[0].x == 3)
    {      // Discriminant is equal to -3.
      if (Plow == 1)
      {
        intToBigInteger(&ValH, 1);
        intToBigInteger(&ValI, 0);
        NonSquareDiscrSolution(value);   // (1, 0)
        CopyBigInt(&ValH, &ValG);
        intToBigInteger(&ValI, -1);
        NonSquareDiscrSolution(value);   // (G, -1)
        intToBigInteger(&ValI, -1);
        BigIntAdd(&ValG, &ValI, &ValH);
        NonSquareDiscrSolution(value);   // (G-1, -1)
        return;
      }
      if (Plow == 3)
      {
        intToBigInteger(&ValI, -1);
        CopyBigInt(&ValH, &ValG);
        subtractdivide(&ValH, 1, 3);
        NonSquareDiscrSolution(value);   // ((G-1)/3, -1)
        intToBigInteger(&ValI, -1);
        CopyBigInt(&ValH, &ValG);
        subtractdivide(&ValH, -2, 3);
        NonSquareDiscrSolution(value);   // ((G+2)/3, -1)
        intToBigInteger(&ValI, -2);
        BigIntAdd(&ValG, &ValG, &ValH);
        subtractdivide(&ValH, -1, 3);
        NonSquareDiscrSolution(value);   // ((2G+1)/3, -2)
        return;
      }
    }
  }
  // Compute bound L = sqrt(4P/(-D))
  multint(&U1, &ValP, 4);
  BigIntDivide(&U1, &discr, &U1);
  BigIntChSign(&U1);               // 4P/(-D)
  squareRoot(U1.limbs, ValL.limbs, U1.nbrLimbs, &ValL.nbrLimbs);  // sqrt(4P/(-D))

  intToBigInteger(&U1, 1);         // Initial value of last convergent: 1/0.
  intToBigInteger(&V1, 0);
  intToBigInteger(&U2, 0);         // Initial value of next to last convergent: 0/1.
  intToBigInteger(&V2, 1);
  // Compute continued fraction expansion of U/V = -Q/2P.
  CopyBigInt(&ValU, &ValQ);
  BigIntChSign(&ValU);
  BigIntAdd(&ValP, &ValP, &ValV);
  for (;;)
  {
    getNextConvergent();
    BigIntSubt(&ValL, &V1, &bigTmp);    // Check whether the denominator of convergent exceeds bound.
    if (bigTmp.sign == SIGN_NEGATIVE)
    {
      break;                            // Bound exceeded, so go out.
    }
    // Test whether P*U1^2 + Q*U1*V1 + R*V1^2 = 1.
    BigIntMultiply(&ValP, &U1, &ValO);      // P*U1
    BigIntMultiply(&ValQ, &V1, &bigTmp);    // Q*V1
    BigIntAdd(&ValO, &bigTmp, &ValO);       // P*U1 + Q*V1
    BigIntMultiply(&ValO, &U1, &ValO);      // P*U1^2 + Q*U1*V1
    BigIntMultiply(&ValR, &V1, &bigTmp);    // R*V1
    BigIntMultiply(&bigTmp, &V1, &bigTmp);  // R*V1^2
    BigIntAdd(&ValO, &bigTmp, &ValO);       // P*U1^2 + Q*U1*V1 + R*V1^2
    if (ValO.sign == SIGN_POSITIVE && ValO.nbrLimbs == 1 && ValO.limbs[0].x == 1)
    {                                       // a*U1^2 + b*U1*V1 + c*V1^2 = 1.
      int D;
      CopyBigInt(&ValH, &U1);
      CopyBigInt(&ValI, &V1);
      NonSquareDiscrSolution(value);        // (U1, V1)
      D = discr.limbs[0].x;
      if (discr.nbrLimbs > 1 || D > 4)
      {                                     // Discriminant is less than -4, go out.
        return;
      }
      if (D == 3 || D == 4)
      {                                     // Discriminant is equal to -3 or -4.
        getNextConvergent();
        CopyBigInt(&ValH, &U1);
        CopyBigInt(&ValI, &V1);
        NonSquareDiscrSolution(value);      // (U1, V1)
        if (D == 3)
        {
          getNextConvergent();
          CopyBigInt(&ValH, &U1);
          CopyBigInt(&ValI, &V1);
          NonSquareDiscrSolution(value);    // (U1, V1)
        }
        return;
      }
    }
  }
}

static void CheckSolutionSquareDiscr(void)
{
  BigIntAdd(&ValK, &currentFactor, &ValJ);     // Get value of J.
  BigIntDivide(&ValZ, &currentFactor, &ValN);
  BigIntAdd(&ValO, &ValN, &ValN);              // Get value of N.
                                               // Compute x = (NI - JM) / den
  BigIntMultiply(&ValN, &ValI, &U1);           // NI
  BigIntMultiply(&ValJ, &ValM, &bigTmp);       // JM
  BigIntSubt(&U1, &bigTmp, &U1);               // NI - JM
  BigIntRemainder(&U1, &ValDen, &U2);
  if (BigIntIsZero(&U2))
  {
    BigIntDivide(&U1, &ValDen, &U1);           // x
                                               // Compute y = (JL - NH) / den
    BigIntMultiply(&ValJ, &ValL, &U2);         // JL
    BigIntMultiply(&ValN, &ValH, &bigTmp);     // NH
    BigIntSubt(&U2, &bigTmp, &U2);             // JL - NH
    BigIntRemainder(&U2, &ValDen, &U3);
    if (BigIntIsZero(&U3))
    {
      BigIntDivide(&U2, &ValDen, &U2);         // y
      ShowXY(&U1, &U2);                        // Show results.
    }
  }
}

static void PerfectSquareDiscriminant(void)
{
  enum eLinearSolution ret;
  struct sFactors *pstFactor;
  int index;

  if (BigIntIsZero(&ValA))
  { // Let R = gcd(b, c)
    BigIntGcd(&ValB, &ValC, &ValR);
  }
  else
  {
    // Multiplying by 4a we get (2aX + (b+g)Y)(2aX + (b-g)Y) = 4ak
    // Let R = gcd(2a, b+g)
    BigIntAdd(&ValA, &ValA, &U1);
    BigIntAdd(&ValB, &ValG, &U2);
    BigIntGcd(&U1, &U2, &ValR);
    // Let S = gcd(2a, b-g)
    BigIntSubt(&ValB, &ValG, &U2);
    BigIntGcd(&U1, &U2, &ValS);
  }
  if (BigIntIsZero(&ValK))
  {     // k equals zero.
    if (BigIntIsZero(&ValA))
    {    // Coefficient a does equals zero.
      // Solve Dy + beta = 0
      intToBigInteger(&Aux[0], 0);
      ret = Linear(&Aux[0], &discr, &ValBeta);
      if (teach)
      {
        w("<div class=\"result\">");
      }
      PrintLinear(ret, "t");
      if (teach)
      {
        w("</span>");
      }
      // Solve bDx + cDy + b*alpha + c*beta = 0
      BigIntMultiply(&ValB, &discr, &Aux[0]);
      BigIntMultiply(&ValC, &discr, &Aux[1]);
      BigIntMultiply(&ValB, &ValAlpha, &Aux[2]);
      BigIntMultiply(&ValC, &ValBeta, &bigTmp);
      BigIntAdd(&Aux[2], &bigTmp, &Aux[2]);
    }
    else
    {    // Coefficient a does not equal zero.
      // Solve 2aD x + (b+g)D y = 2a*alpha + (b+g)*beta
      BigIntMultiply(&ValA, &discr, &Aux[0]);
      BigIntAdd(&Aux[0], &Aux[0], &Aux[0]);
      BigIntAdd(&ValB, &ValG, &Aux[1]);
      BigIntMultiply(&Aux[1], &discr, &Aux[1]);
      BigIntMultiply(&ValA, &ValAlpha, &Aux[2]);
      BigIntAdd(&Aux[2], &Aux[2], &Aux[2]);
      BigIntAdd(&ValB, &ValG, &bigTmp);
      BigIntMultiply(&bigTmp, &ValBeta, &bigTmp);
      BigIntAdd(&Aux[2], &bigTmp, &Aux[2]);
      BigIntChSign(&Aux[2]);
      ret = Linear(&Aux[0], &Aux[1], &Aux[2]);
      if (teach)
      {
        w("<div class=\"result\">");
      }
      PrintLinear(ret, "t");
      if (teach)
      {
        w("</span>");
      }
      // Solve 2aD x + (b-g)D y = 2a*alpha + (b-g)*beta
      BigIntMultiply(&ValA, &discr, &Aux[0]);
      BigIntAdd(&Aux[0], &Aux[0], &Aux[0]);
      BigIntSubt(&ValB, &ValG, &Aux[1]);
      BigIntMultiply(&Aux[1], &discr, &Aux[1]);
      BigIntMultiply(&ValA, &ValAlpha, &Aux[2]);
      BigIntAdd(&Aux[2], &Aux[2], &Aux[2]);
      BigIntSubt(&ValB, &ValG, &bigTmp);
      BigIntMultiply(&bigTmp, &ValBeta, &bigTmp);
      BigIntAdd(&Aux[2], &bigTmp, &Aux[2]);
      BigIntChSign(&Aux[2]);
    }
    ret = Linear(&Aux[0], &Aux[1], &Aux[2]);
    if (teach)
    {
      w("<div class=\"result\">");
    }
    PrintLinear(ret, "t");
    if (teach)
    {
      w("</span>");
    }
    return;
  }
  // k does not equal zero.
  if (BigIntIsZero(&ValA))
  { // If R does not divide k, there is no solution.
    CopyBigInt(&U3, &ValK);
    CopyBigInt(&U1, &ValR);
  }
  else
  {
    // If R*S does not divide 4ak, there is no solution.
    BigIntMultiply(&ValR, &ValS, &U1);
    BigIntMultiply(&ValA, &ValK, &U2);
    multadd(&U3, 4, &U2, 0);
  }
  BigIntRemainder(&U3, &U1, &U2);
  if (!BigIntIsZero(&U2))
  {
    return;
  }
  BigIntDivide(&U3, &U1, &ValZ);
  // Compute all factors of Z = 4ak/RS
  NumberLength = ValZ.nbrLimbs;
  CompressBigInteger(nbrToFactor, &ValZ);
  Bin2Dec(ValZ.limbs, tofactorDec, ValZ.nbrLimbs, groupLen);
  factor(&ValZ, nbrToFactor, factorsMod, astFactorsMod, NULL);
  // x = (NI - JM) / D(IL - MH) and y = (JL - NH) / D(IL - MH)
  // The denominator cannot be zero here.
  // H = 2a/R, I = (b+g)/R, J = F + H * alpha + I * beta
  // L = 2a/S, M = (b-g)/S, N = Z/F + L * alpha + M * beta
  // F is any factor of Z (positive or negative).
  if (BigIntIsZero(&ValA))
  {
    intToBigInteger(&ValH, 0);              // H <- 0
    intToBigInteger(&ValI, 1);              // I <- 1
    BigIntDivide(&ValB, &ValR, &ValL);      // L <- b/R
    BigIntDivide(&ValC, &ValR, &ValM);      // M <- c/R
  }
  else
  {
    BigIntAdd(&ValA, &ValA, &U3);           // 2a
    BigIntDivide(&U3, &ValR, &ValH);        // H <- 2a/R
    BigIntDivide(&U3, &ValS, &ValL);        // L <- 2a/S
    BigIntAdd(&ValB, &ValG, &ValI);
    BigIntDivide(&ValI, &ValR, &ValI);      // I <- (b+g)/R
    BigIntSubt(&ValB, &ValG, &ValM);
    BigIntDivide(&ValM, &ValS, &ValM);      // M <- (b-g)/S
  }
  BigIntMultiply(&ValH, &ValAlpha, &ValK);  // H * alpha
  BigIntMultiply(&ValI, &ValBeta, &bigTmp); // I * beta
  BigIntAdd(&ValK, &bigTmp, &ValK);         // K <- H * alpha + I * beta
  BigIntMultiply(&ValL, &ValAlpha, &ValO);  // L * alpha
  BigIntMultiply(&ValM, &ValBeta, &bigTmp); // M * beta
  BigIntAdd(&ValO, &bigTmp, &ValO);         // O <- L * alpha + M * beta
  // Compute denominator: D(IL - MH)
  BigIntMultiply(&ValI, &ValL, &ValDen);    // IL
  BigIntMultiply(&ValM, &ValH, &bigTmp);    // MH
  BigIntSubt(&ValDen, &bigTmp, &ValDen);    // IL - MH
  BigIntMultiply(&ValDen, &discr, &ValDen); // D(IL - MH)

  // Loop that finds all factors of Z.
  // Use Gray code to use only one big number.
  // Gray code: 0->000, 1->001, 2->011, 3->010, 4->110, 5->111, 6->101, 7->100.
  // Change from zero to one means multiply, otherwise divide.
  memset(counters, 0, sizeof(counters));
  memset(isDescending, 0, sizeof(isDescending));
  intToBigInteger(&currentFactor, 1);
  nbrFactors = astFactorsMod[0].multiplicity;
  for (;;)
  {
    CheckSolutionSquareDiscr();       // Process positive divisor.
    BigIntChSign(&currentFactor);
    CheckSolutionSquareDiscr();       // Process negative divisor.
    BigIntChSign(&currentFactor);
    pstFactor = &astFactorsMod[1];
    for (index = 0; index < nbrFactors; index++, pstFactor++)
    {            // Loop that increments counters.
      if (isDescending[index] == 0)
      {          // Ascending.
        if (counters[index] == pstFactor->multiplicity)
        {
          isDescending[index] = 1;    // Next time it will be descending.
          continue;
        }
        else
        {
          UncompressBigInteger(pstFactor->ptrFactor, &prime);
          BigIntMultiply(&currentFactor, &prime, &currentFactor);
          counters[index]++;
          break;
        }
      }
      else
      {         // Descending.
        if (counters[index] == 0)
        {
          isDescending[index] = 0;    // Next time it will be ascending.
          continue;
        }
        else
        {
          UncompressBigInteger(pstFactor->ptrFactor, &prime);
          BigIntDivide(&currentFactor, &prime, &currentFactor);
          counters[index]--;
          break;
        }
      }
    }
    if (index == nbrFactors)
    {               // All factors have been found. Exit loop.
      break;
    }
  }
}

static void PositiveDiscriminant(void)
{
  callbackQuadModType = CBACK_QMOD_HYPERBOLIC;
  NonSquareDiscriminant();
}

// PQa algorithm for (P+G)/Q where G = sqrt(discriminant):
// U1 <- 1, U2 <- 0
// V1 <- 0, V2 <- 1
// Perform loop:
// a = floor((U + G)/V)
// U3 <- U2, U2 <- U1, U1 <- a*U2 + U3
// V3 <- V2, V2 <- V1, V1 <- a*V2 + V3
// U <- a*V - U
// V <- (D - U^2)/V
// Inside period when: 0 <= G - U < V
static void ContFrac(BigInteger *value, enum eShowSolution solutionNbr)
{
  int index = 0;
  int periodIndex = 0;
  char Beven = ((ValB.limbs[0].x & 1) == 0);
  // Initialize variables.
  intToBigInteger(&U1, 1);
  intToBigInteger(&U2, 0);
  intToBigInteger(&V1, 0);
  intToBigInteger(&V2, 1);
  intToBigInteger(&startPeriodU, -1);      // Less than zero means outside period.
  intToBigInteger(&startPeriodV, -1);
  if (solutionNbr == SECOND_SOLUTION)
  {
    index++;
  }
  for (;;)
  {
    if (ValV.nbrLimbs == 1 && ValV.limbs[0].x == (Beven ? 1 : 2) &&
      (index & 1) == (ValK.sign == ValV.sign ? 0 : 1))
    {         // Found solution.
      if (discr.nbrLimbs == 1 && discr.limbs[0].x == 5 && ValA.sign != ValK.sign && 
        solutionNbr == FIRST_SOLUTION)
      {       // Determinant is 5 and aK < 0. Use exceptional solution (U1-U2)/(V1-V2).
        BigIntSubt(&V1, &V2, &ValH);
        BigIntSubt(&U1, &U2, &ValI);
      }
      else
      {       // Determinant is not 5 or aK > 0. Use convergent U1/V1 as solution.
        CopyBigInt(&ValH, &V1);
        CopyBigInt(&ValI, &U1);
      }
      BigIntChSign(&ValI);
      showSolution = TWO_SOLUTIONS;
      solFound = 0;
      NonSquareDiscrSolution(value);
      if (solFound)
      {
        break;                             // Solution found. Exit loop.
      }
    }
    if (startPeriodU.sign == SIGN_POSITIVE)
    {               // Already inside period.
      periodIndex++;
      if (BigIntEqual(&ValU, &startPeriodU) && BigIntEqual(&ValV, &startPeriodV))
      {             // New period started.
        if ((periodIndex & 1) == 0)
        {           // Two periods of period length is odd, one period if even.
          break;  // Go out in this case.
        }
      }
    }
    else
    {             // Check if periodic part of continued fraction has started.
      BigIntSubt(&ValG, &ValU, &bigTmp);          // G - U must be >= 0 and < V.
      if (bigTmp.sign == SIGN_POSITIVE)
      {
        BigIntSubt(&bigTmp, &ValV, &bigTmp);
        if (bigTmp.sign == SIGN_NEGATIVE)
        {                                         // Periodic part of continued fraction started.
          CopyBigInt(&startPeriodU, &ValU);       // Save U and V to check period end.
          CopyBigInt(&startPeriodV, &ValV);
        }
      }
    }
    BigIntAdd(&ValU, &ValG, &bigTmp);
    if (ValV.sign == SIGN_NEGATIVE)
    {   // If denominator is negative, round square root upwards.
      addbigint(&bigTmp, 1);
    }
    floordiv(&bigTmp, &ValV, &Tmp1);       // Tmp1 = Term of continued fraction.
    CopyBigInt(&U3, &U2);                  // U3 <- U2, U2 <- U1, U1 <- a*U2 + U3
    CopyBigInt(&U2, &U1);
    BigIntMultiply(&Tmp1, &U2, &U1);
    BigIntAdd(&U1, &U3, &U1);
    CopyBigInt(&V3, &V2);                  // V3 <- V2, V2 <- V1, V1 <- a*V2 + V3
    CopyBigInt(&V2, &V1);
    BigIntMultiply(&Tmp1, &V2, &V1);
    BigIntAdd(&V1, &V3, &V1);
    BigIntMultiply(&Tmp1, &ValV, &bigTmp); // U <- a*V - U
    BigIntSubt(&bigTmp, &ValU, &ValU);
    BigIntMultiply(&ValU, &ValU, &bigTmp); // V <- (D - U^2)/V
    BigIntSubt(&ValL, &bigTmp, &bigTmp);
    BigIntDivide(&bigTmp, &ValV, &Tmp1);
    CopyBigInt(&ValV, &Tmp1);
    index++;
  }
}

static void ShowRecSol(char variable, BigInteger *coefX,
                       BigInteger *coefY, BigInteger *coefInd)
{
  enum eLinearSolution t;
  *ptrOutput++ = variable;
  w("<sub>n+1</sub> = ");
  t = Show(coefX, "x<sub>n</sub>", SOLUTION_FOUND);
  t = Show(coefY, "y<sub>n</sub>", t);
  Show1(coefInd, t);
}

static void ShowResult(char *text, BigInteger *value)
{
  w(text);
  w(" = ");
  shownbr(value);
  w("<br>");
}

static void ShowAllRecSols(void)
{
  if (ValP.nbrLimbs > 2 || ValQ.nbrLimbs > 2)
  {
    if (BigIntIsZero(&ValAlpha) && BigIntIsZero(&ValBeta))
    {
      w("x<sub>n+1</sub> = P&nbsp;&#8290;x<sub>n</sub> + Q&nbsp;&#8290;y<sub>n</sub><br>"
        "y<sub>n+1</sub> = R&nbsp;&#8290;x<sub>n</sub> + S&nbsp;&#8290;y<sub>n</sub></p><p>");
    }
    else
    {
      w("x<sub>n+1</sub> = P&nbsp;&#8290;x<sub>n</sub> + Q&nbsp;&#8290;y<sub>n</sub> + K<br>"
        "y<sub>n+1</sub> = R&nbsp;&#8290;x<sub>n</sub> + S&nbsp;&#8290;y<sub>n</sub> + L</p><p>");
    }
    w(lang? "donde:</p><p>": "where:</p><p>");
    ShowResult("P", &ValP);
    ShowResult("Q", &ValQ);
    if (!BigIntIsZero(&ValAlpha) || !BigIntIsZero(&ValBeta))
    {
      ShowResult("K", &ValK);
    }
    ShowResult("R", &ValR);
    ShowResult("S", &ValS);
    if (!BigIntIsZero(&ValAlpha) || !BigIntIsZero(&ValBeta))
    {
      ShowResult("L", &ValL);
    }
  }
  else
  {
    ShowRecSol('x', &ValP, &ValQ, &ValK);
    *ptrOutput++ = '<';
    *ptrOutput++ = 'b';
    *ptrOutput++ = 'r';
    *ptrOutput++ = '>';
    ShowRecSol('y', &ValR, &ValS, &ValL);
  }
  // Compute x_{n-1} from x_n and y_n
  // Compute new value of K and L as: Knew <- L*Q - K*S and Lnew <- K*R - L*P
  BigIntMultiply(&ValL, &ValQ, &Tmp1);
  BigIntMultiply(&ValK, &ValS, &bigTmp);
  BigIntSubt(&Tmp1, &bigTmp, &Tmp1);
  BigIntMultiply(&ValK, &ValR, &Tmp2);
  BigIntMultiply(&ValL, &ValP, &bigTmp);
  BigIntSubt(&Tmp2, &bigTmp, &ValL);
  CopyBigInt(&ValK, &Tmp1);
  // Compute new values of P, Q, R and S as: Pnew <- S, Qnew <- -Q, Rnew <- -R, Snew <- P
  BigIntChSign(&ValQ);
  BigIntChSign(&ValR);
  CopyBigInt(&bigTmp, &ValP);
  CopyBigInt(&ValP, &ValS);
  CopyBigInt(&ValS, &bigTmp);
  w(lang? "<p>y tambi&eacute;n:</p>": "<p>and also:</p>");
  if (ValP.nbrLimbs > 2 || ValQ.nbrLimbs > 2)
  {
    w("<p>");
    ShowResult("P", &ValP);
    ShowResult("Q", &ValQ);
    if (!BigIntIsZero(&ValAlpha) || !BigIntIsZero(&ValBeta))
    {
      ShowResult("K", &ValK);
    }
    ShowResult("R", &ValR);
    ShowResult("S", &ValS);
    if (!BigIntIsZero(&ValAlpha) || !BigIntIsZero(&ValBeta))
    {
      ShowResult("L", &ValL);
    }
  }
  else
  {
    ShowRecSol('x', &ValP, &ValQ, &ValK);
    *ptrOutput++ = '<';
    *ptrOutput++ = 'b';
    *ptrOutput++ = 'r';
    *ptrOutput++ = '>';
    ShowRecSol('y', &ValR, &ValS, &ValL);
  }
  *ptrOutput++ = '<';
  *ptrOutput++ = '/';
  *ptrOutput++ = 'p';
  *ptrOutput++ = '>';
}
// Use continued fraction of sqrt(B^2-4AC)
// If the convergent is r/s we get:
// x(n+1) = Px(n) + Qy(n) + K
// y(n+1) = Rx(n) + Sy(n) + L
// where if b is odd: 
//        P = (r - bs)/2, Q = -cs, R = as, S = (r + bs)/2,
// if b is even:
//        P = r - (b/2)s, Q = -cs, R = as, S = r + (b/2)s,
// in any case:
//        K = (alpha*(1-P) - beta*Q) / D, L = (-alpha*R + beta*(1-S)) / D.
static void ContFracPell(void)
{
  enum eSign sign = SIGN_POSITIVE;
  char Beven = ((ValB.limbs[0].x & 1) == 0);
  int periodNbr = 0;
  int periodLength;
  // Initialize variables.
  intToBigInteger(&ValU, 0);
  intToBigInteger(&ValV, 1);
  CopyBigInt(&ValH, &discr);
  if (Beven)
  {
    subtractdivide(&ValH, 0, 4);
  }
  squareRoot(ValH.limbs, ValG.limbs, ValH.nbrLimbs, &ValG.nbrLimbs);
  ValG.sign = SIGN_POSITIVE;          // g <- sqrt(discr).
  BigIntMultiply(&discr, &ValGcdHomog, &discr);
  BigIntMultiply(&discr, &ValGcdHomog, &discr);  // Obtain original discriminant.
  intToBigInteger(&U1, 1);
  intToBigInteger(&U2, 0);
  intToBigInteger(&V1, 0);
  intToBigInteger(&V2, 1);
  periodLength = 1;
  if (ValGcdHomog.limbs[0].x != 1)
  {
    periodLength = -1;
    do
    {
      BigIntAdd(&ValU, &ValG, &bigTmp);
      if (ValV.sign == SIGN_NEGATIVE)
      {   // If denominator is negative, round square root upwards.
        addbigint(&bigTmp, 1);
      }
      floordiv(&bigTmp, &ValV, &Tmp1);       // Tmp1 = Term of continued fraction.
      BigIntMultiply(&Tmp1, &ValV, &bigTmp); // U <- a*V - U
      BigIntSubt(&bigTmp, &ValU, &ValU);
      BigIntMultiply(&ValU, &ValU, &bigTmp); // V <- (D - U^2)/V
      BigIntSubt(&ValL, &bigTmp, &bigTmp);
      BigIntDivide(&bigTmp, &ValV, &Tmp1);
      CopyBigInt(&ValV, &Tmp1);
      if (periodLength < 0)
      {
        CopyBigInt(&ValUBak, &ValU);
        CopyBigInt(&ValVBak, &ValV);
      }
      periodLength++;
    } while (periodLength == 1 || !BigIntEqual(&ValU, &ValUBak) || !BigIntEqual(&ValV, &ValVBak));
    intToBigInteger(&ValU, 0);    // Reset values of U and V.
    intToBigInteger(&ValV, 1);
  }
  w(lang ? "<p>Soluciones recursivas:</p><p>" : "<p>Recursive solutions:</p><p>");
  for (;;)
  {
    int limbValue;
    BigIntAdd(&ValU, &ValG, &bigTmp);
    if (ValV.sign == SIGN_NEGATIVE)
    {   // If denominator is negative, round square root upwards.
      addbigint(&bigTmp, 1);
    }
    floordiv(&bigTmp, &ValV, &Tmp1);       // Tmp1 = Term of continued fraction.
    CopyBigInt(&U3, &U2);                  // U3 <- U2, U2 <- U1, U1 <- a*U2 + U3
    CopyBigInt(&U2, &U1);
    BigIntMultiply(&Tmp1, &U2, &U1);
    BigIntAdd(&U1, &U3, &U1);
    CopyBigInt(&V3, &V2);                  // V3 <- V2, V2 <- V1, V1 <- a*V2 + V3
    CopyBigInt(&V2, &V1);
    BigIntMultiply(&Tmp1, &V2, &V1);
    BigIntAdd(&V1, &V3, &V1);
    BigIntMultiply(&Tmp1, &ValV, &bigTmp); // U <- a*V - U
    BigIntSubt(&bigTmp, &ValU, &ValU);
    BigIntMultiply(&ValU, &ValU, &bigTmp); // V <- (D - U^2)/V
    BigIntSubt(&ValH, &bigTmp, &bigTmp);
    BigIntDivide(&bigTmp, &ValV, &Tmp1);
    CopyBigInt(&ValV, &Tmp1);
    if (sign == SIGN_POSITIVE)
    {
      sign = SIGN_NEGATIVE;
    }
    else
    {
      sign = SIGN_POSITIVE;
    }
    // Expecting denominator to be 1 (B even or odd) or 4 (B odd) with correct sign.
    if (ValV.nbrLimbs != 1 || ValV.sign != sign)
    {   
      continue;
    }
    limbValue = ValV.limbs[0].x;
    if (limbValue != 1 && (Beven || limbValue != 4))
    {
      continue;
    }
    if (++periodNbr*periodLength % ValGcdHomog.limbs[0].x != 0)
    {
      continue;
    }
    // Found solution.
    if (Beven)
    {
      CopyBigInt(&ValQ, &ValB);
      subtractdivide(&ValQ, 0, 2);
      BigIntMultiply(&ValQ, &V1, &bigTmp);
      BigIntSubt(&U1, &bigTmp, &ValP);        // P <- r - (b/2)s
      BigIntAdd(&U1, &bigTmp, &ValS);         // S <- r + (b/2)s
    }
    else
    {
      BigIntMultiply(&ValB, &V1, &bigTmp);
      BigIntSubt(&U1, &bigTmp, &ValP);        // P <- r - bs
      BigIntAdd(&U1, &bigTmp, &ValS);         // S <- r + bs
      if (limbValue == 4)
      {
        subtractdivide(&ValP, 0, 2);          // P <- (r - bs)/2
        subtractdivide(&ValS, 0, 2);          // S <- (r + bs)/2
      }
    }
    BigIntMultiply(&ValC, &V1, &ValQ);
    BigIntChSign(&ValQ);                      // Q <- -cs
    BigIntMultiply(&ValA, &V1, &ValR);        // R <- as
    if (!Beven && limbValue == 1)
    {
      BigIntAdd(&ValQ, &ValQ, &ValQ);         // Q <- -2cs
      BigIntAdd(&ValR, &ValR, &ValR);         // R <- 2as
    }
    BigIntMultiply(&ValAlpha, &ValP, &ValK);
    BigIntMultiply(&ValBeta, &ValQ, &bigTmp);
    BigIntAdd(&ValK, &bigTmp, &ValK);
    BigIntMultiply(&ValAlpha, &ValR, &ValL);
    BigIntMultiply(&ValBeta, &ValS, &bigTmp);
    BigIntAdd(&ValL, &bigTmp, &ValL);
    // Check whether alpha - K and beta - L are multiple of discriminant.
    BigIntSubt(&ValAlpha, &ValK, &bigTmp);
    BigIntRemainder(&bigTmp, &discr, &bigTmp);
    if (BigIntIsZero(&bigTmp))
    {
      BigIntSubt(&ValBeta, &ValL, &bigTmp);
      BigIntRemainder(&bigTmp, &discr, &bigTmp);
      if (BigIntIsZero(&bigTmp))
      {    // Solution found.
        BigIntSubt(&ValAlpha, &ValK, &ValK);
        BigIntDivide(&ValK, &discr, &ValK);
        BigIntSubt(&ValBeta, &ValL, &ValL);
        BigIntDivide(&ValL, &discr, &ValL);
        ShowAllRecSols();
        return;
      }
    }
    // Check whether alpha + K and beta + L are multiple of discriminant.
    BigIntAdd(&ValAlpha, &ValK, &bigTmp);
    BigIntRemainder(&bigTmp, &discr, &bigTmp);
    if (BigIntIsZero(&bigTmp))
    {
      BigIntAdd(&ValBeta, &ValL, &bigTmp);
      BigIntRemainder(&bigTmp, &discr, &bigTmp);
      if (BigIntIsZero(&bigTmp))
      {    // Solution found.
        BigIntAdd(&ValAlpha, &ValK, &ValK);
        BigIntDivide(&ValK, &discr, &ValK);
        BigIntAdd(&ValBeta, &ValL, &ValL);
        BigIntDivide(&ValL, &discr, &ValL);
        BigIntChSign(&ValP);
        BigIntChSign(&ValQ);
        BigIntChSign(&ValR);
        BigIntChSign(&ValS);
        ShowAllRecSols();
        return;
      }
    }
  }
}

static void callbackQuadModHyperbolic(BigInteger *value)
{
  // Compute theta = (t-b)/(2a) mod K. (t = value).
  char Beven = ((ValB.limbs[0].x & 1) == 0);
  // Compute P = floor((2*a*theta + b)/2)
  BigIntAdd(&ValA, &ValA, &ValP);
  BigIntMultiply(&ValP, value, &ValP);
  BigIntAdd(&ValP, &ValB, &ValP);
  subtractdivide(&ValP, ValP.limbs[0].x & 1, 2);
  // Compute Q = a*abs(K)
  CopyBigInt(&ValQ, &ValK);
  ValQ.sign = SIGN_POSITIVE;
  BigIntMultiply(&ValQ, &ValA, &ValQ);
  // Find U, V, L so we can compute the continued fraction expansion of (U+sqrt(L))/V.
  CopyBigInt(&ValL, &discr);
  if (Beven)
  {
    CopyBigInt(&ValU, &ValP);       // U <- P
    subtractdivide(&ValL, 0, 4);    // Argument of square root is discriminant/4.
    CopyBigInt(&ValV, &ValQ);
  }
  else
  {
    BigIntAdd(&ValP, &ValP, &ValU);
    addbigint(&ValU, 1);            // U <- 2P+1
    BigIntAdd(&ValQ, &ValQ, &ValV); // V <- 2Q
  }
  BigIntChSign(&ValU);
  // If L-U^2 is not multiple of V, there is no solution, so go out.
  BigIntMultiply(&ValU, &ValU, &bigTmp);
  BigIntSubt(&ValL, &bigTmp, &bigTmp);
  BigIntRemainder(&bigTmp, &ValV, &bigTmp);
  if (!BigIntIsZero(&bigTmp))
  {
    return;
  }
  // Set G <- floor(sqrt(L))
  squareRoot(ValL.limbs, ValG.limbs, ValL.nbrLimbs, &ValG.nbrLimbs);
  ValG.sign = SIGN_POSITIVE;          // g <- sqrt(discr).
  Xplus.nbrLimbs = 0;                 // Invalidate solutions.
  Xminus.nbrLimbs = 0;
  Yplus.nbrLimbs = 0;
  Yminus.nbrLimbs = 0;
  CopyBigInt(&ValUBak, &ValU);
  CopyBigInt(&ValVBak, &ValV);
  ContFrac(value, FIRST_SOLUTION);    // Continued fraction of (U+G)/V
  CopyBigInt(&ValU, &ValUBak);
  CopyBigInt(&ValV, &ValVBak);
  BigIntChSign(&ValU);
  BigIntChSign(&ValV);
  ContFrac(value, SECOND_SOLUTION);   // Continued fraction of (-U+G)/(-V)
  showSolution = ONE_SOLUTION;
  if (Xplus.nbrLimbs)
  {
    ShowXY(&Xplus, &Yplus);
  }
  if (Xminus.nbrLimbs)
  {
    ShowXY(&Xminus, &Yminus);
  }
}

void SolveQuadEquation(void)
{
  also = 0;
  showSolution = ONE_SOLUTION;
  showRecursiveSolution = 0;    // Do not show recursive solution by default.
  // Get GCD of A, B, C, D and E.
  BigIntGcd(&ValA, &ValB, &Aux[0]);
  BigIntGcd(&Aux[0], &ValC, &Aux[1]);
  BigIntGcd(&Aux[1], &ValD, &Aux[0]);
  BigIntGcd(&Aux[0], &ValE, &Aux[1]);
  BigIntRemainder(&ValF, &Aux[1], &Aux[0]);
  if (!BigIntIsZero(&Aux[0]))
  {   // F is not multiple of GCD(A, B, C, D, E) so there are no solutions.
    w(lang ? "<p>No hay soluciones.</p>" : "<p>There are no solutions.</p>");
    ptrOutput += strlen(ptrOutput);
  }
  // Divide all coefficients by GCD(A, B, C, D, E)
  BigIntDivide(&ValA, &Aux[1], &ValA);
  BigIntDivide(&ValB, &Aux[1], &ValB);
  BigIntDivide(&ValC, &Aux[1], &ValC);
  BigIntDivide(&ValD, &Aux[1], &ValD);
  BigIntDivide(&ValE, &Aux[1], &ValE);
  BigIntDivide(&ValF, &Aux[1], &ValF);
  // Test whether the equation is linear. A = B = C = 0.
  if (BigIntIsZero(&ValA) && BigIntIsZero(&ValB) && BigIntIsZero(&ValC))
  {
    enum eLinearSolution ret = Linear(&ValD, &ValE, &ValF);
    if (teach)
    {
      w("<div class=\"result\">");
    }
    PrintLinear(ret, "t");
    if (teach)
    {
      w("</span>");
    }
    return;
  }
  // Compute discriminant: b^2 - 4ac.
  BigIntMultiply(&ValB, &ValB, &Aux[0]);
  BigIntMultiply(&ValA, &ValC, &Aux[1]);
  multint(&Aux[2], &Aux[1], 4);
  BigIntSubt(&Aux[0], &Aux[2], &discr);
  if (BigIntIsZero(&discr))
  {        // Discriminant is zero.
    DiscriminantIsZero();
    return;
  }
  // Discriminant is not zero. Translate the origin (x, y) by (alpha, beta).
  // Compute alpha = 2cd - be
  BigIntMultiply(&ValC, &ValD, &ValAlpha);
  BigIntAdd(&ValAlpha, &ValAlpha, &ValAlpha);
  BigIntMultiply(&ValB, &ValE, &bigTmp);
  BigIntSubt(&ValAlpha, &bigTmp, &ValAlpha);
  // Compute beta = 2ae - bd
  BigIntMultiply(&ValA, &ValE, &ValBeta);
  BigIntAdd(&ValBeta, &ValBeta, &ValBeta);
  BigIntMultiply(&ValB, &ValD, &bigTmp);
  BigIntSubt(&ValBeta, &bigTmp, &ValBeta);
  // We get the equation ax^2 + bxy + cy^2 = k
  // where k = -D (ae^2 - bed + cd^2 + fD)
  BigIntMultiply(&ValA, &ValE, &ValK);     // ae
  BigIntMultiply(&ValK, &ValE, &ValK);     // ae^2
  BigIntMultiply(&ValB, &ValE, &bigTmp);   // be
  BigIntMultiply(&bigTmp, &ValD, &bigTmp); // bed
  BigIntSubt(&ValK, &bigTmp, &ValK);       // ae^2 - bed
  BigIntMultiply(&ValC, &ValD, &bigTmp);   // cd
  BigIntMultiply(&bigTmp, &ValD, &bigTmp); // cd^2
  BigIntAdd(&ValK, &bigTmp, &ValK);        // ae^2 - bed + cd^2
  BigIntMultiply(&ValF, &discr, &bigTmp);  // fD
  BigIntAdd(&ValK, &bigTmp, &ValK);        // ae^2 - bed + cd^2 + fD
  BigIntMultiply(&ValK, &discr, &ValK);    // D (ae^2 - bed + cd^2 + fD)
  BigIntChSign(&ValK);                     // k
                                           // Let t=gcd(a,b,c). If k is not multiple of t, there are no solutions.
  BigIntGcd(&ValA, &ValB, &bigTmp);
  BigIntGcd(&bigTmp, &ValC, &U1);
  BigIntRemainder(&ValK, &U1, &bigTmp);
  if (!BigIntIsZero(&bigTmp))
  {
    return;     // There are no solutions.
  }
  if (discr.sign == SIGN_NEGATIVE)
  {
    NegativeDiscriminant();
    return;
  }
  squareRoot(discr.limbs, ValG.limbs, discr.nbrLimbs, &ValG.nbrLimbs);
  ValG.sign = SIGN_POSITIVE;
  BigIntMultiply(&ValG, &ValG, &bigTmp);
  if (BigIntEqual(&bigTmp, &discr))
  {   // Discriminant is a perfect square.
    PerfectSquareDiscriminant();
    return;
  }
  PositiveDiscriminant();
}

void quadText(char *coefAText, char *coefBText, char *coefCText,
              char *coefDText, char *coefEText, char *coefFText)
{
  int coeffNbr;
  char *ptrBeginSol;
  enum eExprErr rc;
  struct stValidateCoeff *pstValidateCoeff = astValidateCoeff;
  astValidateCoeff[0].expression = coefAText;
  astValidateCoeff[1].expression = coefBText;
  astValidateCoeff[2].expression = coefCText;
  astValidateCoeff[3].expression = coefDText;
  astValidateCoeff[4].expression = coefEText;
  astValidateCoeff[5].expression = coefFText;
  ptrOutput = output;
  strcpy(ptrOutput, "2<p>");
  ptrOutput = output + strlen(output);
  for (coeffNbr = 0; coeffNbr < NBR_COEFF; coeffNbr++)
  {
    rc = ComputeExpression(pstValidateCoeff->expression, 1,
                           pstValidateCoeff->bigint);
    if (rc != EXPR_OK)
    {
      strcpy(ptrOutput, lang ? pstValidateCoeff->textSpanish : pstValidateCoeff->textEnglish);
      ptrOutput = output + strlen(output);
      textError(ptrOutput, rc);
      ptrOutput = output + strlen(output);
      strcpy(ptrOutput, "</p>");
      ptrOutput = output + strlen(output);
      break;
    }
    pstValidateCoeff++;
  }
  if (coeffNbr == NBR_COEFF)
  {
    ShowEq(&ValA, &ValB, &ValC, &ValD, &ValE, &ValF, "x", "y");
    strcpy(ptrOutput, " = 0</p>");
    ptrOutput += strlen(ptrOutput);
    SolNbr = 0;
    ptrBeginSol = ptrOutput;
#if 0
    strcpy(ptrOutput, "<ol>");
    ptrOutput += strlen(ptrOutput);
#endif
    SolveQuadEquation();
#if 0
    if (SolNbr == 0)
    {
      ptrOutput = ptrBeginSol;
      strcpy(ptrOutput, lang ? "<p>No hay soluciones.</p>" : "<p>There are no solutions.</p>");
    }
    else
    {
      strcpy(ptrOutput, "</ol>");
    }
    ptrOutput += strlen(ptrOutput);
#endif
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
  char *ptrCoeffA, *ptrCoeffB, *ptrCoeffC;
  char *ptrCoeffD, *ptrCoeffE, *ptrCoeffF;
  groupLen = 0;
  while (*ptrData != ',')
  {
    groupLen = groupLen * 10 + (*ptrData++ - '0');
  }
  ptrData++;                    // Skip comma.
  flags = *ptrData;
  lang = flags & 1;
  ptrCoeffA = ptrData + 2;  // Skip flags and comma.
  ptrCoeffB = ptrCoeffA + strlen(ptrCoeffA) + 1;
  ptrCoeffC = ptrCoeffB + strlen(ptrCoeffB) + 1;
  ptrCoeffD = ptrCoeffC + strlen(ptrCoeffC) + 1;
  ptrCoeffE = ptrCoeffD + strlen(ptrCoeffD) + 1;
  ptrCoeffF = ptrCoeffE + strlen(ptrCoeffE) + 1;
  quadText(ptrCoeffA, ptrCoeffB, ptrCoeffC, ptrCoeffD, ptrCoeffE, ptrCoeffF);
  databack(output);
}
#endif

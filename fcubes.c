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
#include <string.h>
#include "bignbr.h"
#include "expression.h"
#include "highlevel.h"
#include "showtime.h"
#include "batch.h"
static BigInteger value;
static BigInteger Base1;
static BigInteger Base2;
static BigInteger Base3;
static BigInteger Base4;
static BigInteger P;
static BigInteger Q;
static BigInteger R;
static BigInteger S;
static BigInteger a;
static BigInteger b;
static BigInteger P1;
static BigInteger Q1;
static BigInteger R1;
static BigInteger S1;
static BigInteger tmpP1;
static BigInteger tmpQ1;
static BigInteger tmpR1;
static BigInteger tmpS1;
static BigInteger toProcess;
static int groupLength;
static const char *cube = "<span class=\"bigger\">³</span>";
static bool hexadecimal;
static void batchCubesCallback(char** pptrOutput);
static int sums[] =
{
  6, 0, 1, -1, -1, 0, -1, 0, 1, 1,
  6, 3, 1, 0, -1, 4, 2, -5, -2, 4,
  18, 1, 2, 14, -2, -23, -3, -26, 3, 30,
  18, 7, 1, 2, 6, -1, 8, -2, -9, 2,
  18, 8, 1, -5, -1, 14, -3, 29, 3, -30,
  18, 10, 1, 6, -1, -15, -3, -32, 3, 33,
  18, 11, 1, -1, 6, 7, 8, 10, -9, -11,
  18, 17, 2, -12, -2, 21, -3, 23, 3, -27,
  54, 20, 3, -11, -3, 10, 1, 2, -1, 7,
  72, 56, -9, 4, 1, 4, 6, -2, 8, -4,
  108, 2, -1, -22, 1, 4, -3, -41, 3, 43,
  216, 92, 3, -164, -3, 160, 1, -35, -1, 71,
  270, 146, -60, 91, -3, 13, 22, -37, 59, -89,
  270, 200, 3, 259, -3, -254, 1, 62, -1, -107,
  270, 218, -3, -56, 3, 31, -5, -69, 5, 78,
  432, 380, -3, 64, 3, -80, 2, -29, -2, 65,
  540, 38, 5, -285, -5, 267, 3, -140, -3, 190,
  810, 56, 5, -755, -5, 836, 9, -1445, -9, 1420,
  1080, 380, -1, -1438, 1, 1258, -3, -4037, 3, 4057,
  1620, 1334, -5, -3269, 5, 3107, -9, -5714, 9, 5764,
  1620, 1352, -5, 434, 5, -353, 9, -722, -9, 697,
  2160, 362, -5, -180, 5, 108, -6, -149, 6, 199,
  6480, 794, -5, -83, 5, 11, -6, -35, 6, 85,
};

   // Compute tmpP1 <- Base1^3 + Base2^3 + Base3^3 + Base4^3
   // Use tmpQ1 as temporary variable.
static void getSumOfCubes(void)
{
  (void)BigIntMultiply(&Base1, &Base1, &tmpQ1);
  (void)BigIntMultiply(&tmpQ1, &Base1, &tmpP1);

  (void)BigIntMultiply(&Base2, &Base2, &tmpQ1);
  (void)BigIntMultiply(&tmpQ1, &Base2, &tmpQ1);
  BigIntAdd(&tmpP1, &tmpQ1, &tmpP1);

  (void)BigIntMultiply(&Base3, &Base3, &tmpQ1);
  (void)BigIntMultiply(&tmpQ1, &Base3, &tmpQ1);
  BigIntAdd(&tmpP1, &tmpQ1, &tmpP1);

  (void)BigIntMultiply(&Base4, &Base4, &tmpQ1);
  (void)BigIntMultiply(&tmpQ1, &Base4, &tmpQ1);
  BigIntAdd(&tmpP1, &tmpQ1, &tmpP1);
}

 // If |value1| < |value2|, exchange both numbers.
static void SortBigIntegers(BigInteger *pValue1, BigInteger *pValue2)
{
  enum eSign tmpSign;
  int index;
  limb *ptr1;
  limb *ptr2;
  int nbrLimbs1 = pValue1->nbrLimbs;
  int nbrLimbs2 = pValue2->nbrLimbs;
  if (nbrLimbs1 > nbrLimbs2)
  {
    return;    // Base1 > Base2, so nothing to do.
  }
  if (nbrLimbs1 == nbrLimbs2)
  {
    ptr1 = &pValue1->limbs[nbrLimbs1];
    ptr2 = &pValue2->limbs[nbrLimbs1];
    for (index = nbrLimbs1 - 1; index >= 0; index--)
    {
      ptr1--;
      ptr2--;
      if (ptr1->x > ptr2->x)
      {
        return;    // Base1 > Base2, so nothing to do.
      }
      if (ptr1->x < ptr2->x)
      {
        break;     // Base1 < Base2, so exchange them.
      }
    }
  }
    // Exchange lengths.
  pValue1->nbrLimbs = nbrLimbs2;
  pValue2->nbrLimbs = nbrLimbs1;
    // Exchange signs.
  tmpSign = pValue1->sign;
  pValue1->sign = pValue2->sign;
  pValue2->sign = tmpSign;
    // Exchange bytes that compose the numbers.
  ptr1 = pValue1->limbs;
  ptr2 = pValue2->limbs;
  for (index = 0; index < nbrLimbs2; index++)
  {
    limb tmp = *ptr1;
    *ptr1 = *ptr2;
    *ptr2 = tmp;
    ptr1++;
    ptr2++;
  }
}

static void EvaluateQuadraticPoly(BigInteger *pResult, const BigInteger *pValue,
  int quad, int linear, int constant)
{
  multadd(pResult, quad, pValue, linear);
  // Multiply result by value.
  (void)BigIntMultiply(pResult, pValue, pResult);
  addbigint(pResult, constant);
}

// Perform Pell solution of Demjanenko's theorem
// Using these values of P, Q, R and S, a and b will be
// always one and zero (mod 6) respectively.
// P <- -112488782561 = -(52*2^31+819632865)
// Q <- -6578430178320 = -(3063*2^31+687764496)
// R <- -1923517596 = -(0*2^31+1923517596)
// S <- P
// P1 <- 1
// Q1 <- 0
// R1 <- 0
// S1 <- 1
static void Demjanenko(void)
{
  int mod83;
  int pow;
  int exp;
  int mask;
  P.limbs[1].x = 52;
  P.limbs[0].x = 819632865;
  Q.limbs[1].x = 3063;
  Q.limbs[0].x = 687764496;
  Q.nbrLimbs = 2;
  P.nbrLimbs = 2;
  Q.sign = SIGN_NEGATIVE;
  P.sign = SIGN_NEGATIVE;
  CopyBigInt(&S, &P);
  intToBigInteger(&R, -1923517596);
  intToBigInteger(&P1, 1);
  intToBigInteger(&Q1, 0);
  intToBigInteger(&R1, 0);
  intToBigInteger(&S1, 1);
  mod83 = getRemainder(&value, 83);
  pow = 71;
  exp = 0;
  while (pow != mod83)
  {
    exp++;
    pow = (pow * 50) % 83;
  }
  if (exp > (82 / 2))
  {
    exp = 82 - exp;
    Q.sign = SIGN_POSITIVE;
    R.sign = SIGN_POSITIVE;
  }  // Now exp is in range 0-41.
  mask = 32;
  while (mask > 0)
  {
    // Compute tmpP1 as P1*P1 + Q1*R1
    // Compute tmpQ1 as (P1+S1) * Q1
    // Compute tmpR1 as (P1+S1) * R1
    // Compute tmpS1 as S1*S1 + Q1*R1
    // Compute P1 as tmpP1
    // Compute Q1 as tmpQ1
    // Compute R1 as tmpR1
    // Compute S1 as tmpS1
    (void)BigIntMultiply(&P1, &P1, &tmpP1);
    (void)BigIntMultiply(&Q1, &R1, &tmpQ1);
    (void)BigIntMultiply(&S1, &S1, &tmpS1);
    BigIntAdd(&P1, &S1, &tmpR1);
    BigIntAdd(&tmpP1, &tmpQ1, &P1);
    BigIntAdd(&tmpS1, &tmpQ1, &S1);
    (void)BigIntMultiply(&tmpR1, &Q1, &Q1);
    (void)BigIntMultiply(&tmpR1, &R1, &R1);
    if ((exp & mask) != 0)
    {
      // Compute tmpP1 as P*P1 + Q*R1
      // Compute tmpQ1 as P*Q1 + Q*S1
      // Compute tmpR1 as R*P1 + S*R1
      // Compute tmpS1 as R*Q1 + S*S1
      // Compute P1 as tmpP1
      // Compute Q1 as tmpQ1
      // Compute R1 as tmpR1
      // Compute S1 as tmpS1
      (void)BigIntMultiply(&P, &P1, &tmpP1);
      (void)BigIntMultiply(&Q, &R1, &tmpQ1);
      (void)BigIntMultiply(&R, &P1, &tmpR1);
      (void)BigIntMultiply(&S, &R1, &tmpS1);
      BigIntAdd(&tmpP1, &tmpQ1, &P1);
      BigIntAdd(&tmpR1, &tmpS1, &R1);
      (void)BigIntMultiply(&P, &Q1, &tmpP1);
      (void)BigIntMultiply(&Q, &S1, &tmpQ1);
      (void)BigIntMultiply(&R, &Q1, &tmpR1);
      (void)BigIntMultiply(&S, &S1, &tmpS1);
      BigIntAdd(&tmpP1, &tmpQ1, &Q1);
      BigIntAdd(&tmpR1, &tmpS1, &S1);
    }
    mask >>= 1;
  }
  addmult(&a, &P1, -3041, &Q1, -52);   // a <- -3041*P1 - 52*Q1
  addmult(&b, &R1, -3041, &S1, -52);   // b <- -3041*R1 - 52*S1
  addmult(&Base1, &a, 27, &b, -928);   // Base1 <- 27*a - 928*b
  addmult(&Base2, &a, -9, &b, -602);   // Base2 <- -9*a - 602*b
  addmult(&Base3, &a, 25, &b, -2937);  // Base3 <- 25*a - 2937*b
  addmult(&Base4, &a, -19, &b, 2746);  // Base4 <- -19*a - 2746*b
  // a <- (value - Base1^3 - Base2^3 - Base3^3 - Base4^3)/(18*83)
  getSumOfCubes();  // tmpP1 = Base1^3 + Base2^3 + Base3^3 + Base4^3
  BigIntSubt(&value, &tmpP1, &a);
  subtractdivide(&a, 0, 18 * 83);      // Divide a by 18*83.
  multint(&tmpP1, &a, 10);             // Base1 <- Base1 + 10*a
  BigIntAdd(&tmpP1, &Base1, &Base1);
  multint(&tmpP1, &a, -19);            // Base2 <- Base2 - 19*a
  BigIntAdd(&tmpP1, &Base2, &Base2);
  multint(&tmpP1, &a, -24);            // Base3 <- Base3 - 24*a
  BigIntAdd(&tmpP1, &Base3, &Base3);
  multint(&tmpP1, &a, 27);             // Base4 <- Base4 + 27*a
  BigIntAdd(&tmpP1, &Base4, &Base4);
}

static int fcubes(const BigInteger *pArgument)
{
  int mod18;
  int i;
  bool converted = false;
  CopyBigInt(&value, pArgument);
  // Compute argument mod 18.
  mod18 = getRemainder(pArgument, 18);
  if ((mod18 == 4) || (mod18 == 5) || (mod18 == 13) || (mod18 == 14))
  {
    return -1;  // Applet does not work if the number is congruent to 4 or 5 (mod 9)
  }
  if (mod18 == 16)
  {             // Change sign.
    converted = true;
    BigIntNegate(&value, &value);
  }
  for (i = ((int)sizeof(sums) / (int)sizeof(sums[0]))-10; i>=0; i -= 10)
  {
    int modulus = sums[i];
    if (((getRemainder(&value, modulus) + modulus)% modulus) == sums[i + 1])
    {
      subtractdivide(&value, sums[i + 1], modulus);      // value <- (value-sums[i+1])/modulus
      multadd(&Base1, sums[i + 2], &value, sums[i + 3]); // Base1 <- sums[i+2]*value+sums[i+3]
      multadd(&Base2, sums[i + 4], &value, sums[i + 5]); // Base2 <- sums[i+4]*value+sums[i+5]
      multadd(&Base3, sums[i + 6], &value, sums[i + 7]); // Base3 <- sums[i+6]*value+sums[i+7]
      multadd(&Base4, sums[i + 8], &value, sums[i + 9]); // Base4 <- sums[i+8]*value+sums[i+9]
      break;
    }
  }
  if (i < 0)
  {
    if (getRemainder(&value, 54) == 2)
    {           // If value == 2 (mod 54)...
      subtractdivide(&value, 2, 54);   // value <- (value-2)/54
      EvaluateQuadraticPoly(&Base1, &value, 29484, 2211, 43);
      EvaluateQuadraticPoly(&Base2, &value, -29484, -2157, -41);
      EvaluateQuadraticPoly(&Base3, &value, 9828, 485, 4);
      EvaluateQuadraticPoly(&Base4, &value, -9828, -971, -22);
    }
    else if (getRemainder(&value, 83 * 108) == (83 * 46))
    {           // If value == 83*46 (mod 83*108)...
      subtractdivide(&value, 83 * 46, (83 * 108)); // value <-(value - (83*46)) / (83*108)
      EvaluateQuadraticPoly(&Base1, &value, 29484, 25143, 5371);
      EvaluateQuadraticPoly(&Base2, &value, -29484, -25089, -5348);
      EvaluateQuadraticPoly(&Base3, &value, 9828, 8129, 1682);
      EvaluateQuadraticPoly(&Base4, &value, -9828, -8615, -1889);
    }
    else
    {
      Demjanenko();
    }
  }
  if (converted)
  {
    BigIntNegate(&Base1, &Base1);
    BigIntNegate(&Base2, &Base2);
    BigIntNegate(&Base3, &Base3);
    BigIntNegate(&Base4, &Base4);
  }

      // Sort cubes
  SortBigIntegers(&Base1, &Base2);
  SortBigIntegers(&Base1, &Base3);
  SortBigIntegers(&Base1, &Base4);
  SortBigIntegers(&Base2, &Base3);
  SortBigIntegers(&Base2, &Base4);
  SortBigIntegers(&Base3, &Base4);

  // Validate

  getSumOfCubes();  // tmpP1 = Base1^3 - Base2^3 - Base3^3 - Base4^3
  BigIntSubt(&tmpP1, pArgument, &tmpQ1);
  if ((tmpQ1.nbrLimbs != 1) || (tmpQ1.limbs[0].x != 0))
  {
    return 1;       // Result does not validate.
  }
  return 0;
}

void fcubesText(char *input, int grpLen)
{
  char *ptrOutput;
  if (valuesProcessed == 0)
  {
    groupLength = grpLen;
  }
  (void)BatchProcessing(input, &toProcess, &ptrOutput, NULL, batchCubesCallback);
#ifdef __EMSCRIPTEN__
  copyStr(&ptrOutput, lang ? "<p>Transcurrió " : "<p>Time elapsed: ");
  int elapsedTime = (int)(tenths() - originalTenthSecond);
  GetDHMSt(&ptrOutput, elapsedTime);
#endif
  copyStr(&ptrOutput, "</p><p>");
  copyStr(&ptrOutput, (lang ? COPYRIGHT_SPANISH: COPYRIGHT_ENGLISH));
  copyStr(&ptrOutput, "</p>");
}

// Show cube number. Use parentheses for negative numbers.
static void showCube(char** pptrOutput, const BigInteger* pBase)
{
  char* ptrOutput = *pptrOutput;
  if (pBase->sign == SIGN_NEGATIVE)
  {
    *ptrOutput = '(';
    ptrOutput++;
  }
  if (hexadecimal)
  {
    BigInteger2Hex(&ptrOutput, pBase, groupLength);
  }
  else
  {
    BigInteger2Dec(&ptrOutput, pBase, groupLength);
  }
  if (pBase->sign == SIGN_NEGATIVE)
  {
    *ptrOutput = ')';
    ptrOutput++;
  }
  copyStr(&ptrOutput, cube);
  *pptrOutput = ptrOutput;
}

static void batchCubesCallback(char **pptrOutput)
{
  int result;
  char *ptrOutput = *pptrOutput;
  NumberLength = toProcess.nbrLimbs;
  result = fcubes(&toProcess);
  // Show the number to be decomposed into sum of cubes.
  copyStr(&ptrOutput, "<p>");
  if (hexadecimal)
  {
    BigInteger2Hex(&ptrOutput, &toProcess, groupLength);
  }
  else
  {
    BigInteger2Dec(&ptrOutput, &toProcess, groupLength);
  }
  switch (result)
  {
  case -1:
    copyStr(&ptrOutput, (lang?": El applet no funciona si el número es congruente a 4 o 5 (mod 9)</p>":
      ": This applet does not work if the number is congruent to 4 or 5 (mod 9)</p>"));
    *pptrOutput = ptrOutput;
    return;
  case 1:
    copyStr(&ptrOutput, (lang?": ¡Error interno! Por favor envíe este número al autor del applet.</p>":
      ": Internal error! Please send the number to the author of the applet.</p>"));
    *pptrOutput = ptrOutput;
    return;
  case 2:
    copyStr(&ptrOutput, (lang?": El usuario detuvo el cálculo</p>": ": User stopped the calculation</p>"));
    *pptrOutput = ptrOutput;
    return;
  default:
    break;
  }
  // Show decomposition in sum of 1, 2, 3 or 4 cubes.
  copyStr(&ptrOutput, " = ");
  showCube(&ptrOutput, &Base1);
  if (!BigIntIsZero(&Base2))
  {
    copyStr(&ptrOutput, " + ");
    showCube(&ptrOutput, &Base2);
  }
  if (!BigIntIsZero(&Base3))
  {
    copyStr(&ptrOutput, " + ");
    showCube(&ptrOutput, &Base3);
  }
  if (!BigIntIsZero(&Base4))
  {
    copyStr(&ptrOutput, " + ");
    showCube(&ptrOutput, &Base4);
  }
  *pptrOutput = ptrOutput;
}

#if defined __EMSCRIPTEN__ && !defined _MSC_VER
EXTERNALIZE void doWork(void)
{
  int app;
  int grpLen = 0;
  char* ptrData = inputString;
  originalTenthSecond = tenths();
  if (*ptrData == 'C')
  {    // User pressed Continue button.
    fcubesText(NULL, 0); // Routine does not use parameters in this case.
    databack(output);
    return;
  }
  valuesProcessed = 0;
  while (*ptrData != ',')
  {
    grpLen = (grpLen * 10) + (*ptrData - '0');
    ptrData++;
  }
  ptrData++;             // Skip comma.
  app = *ptrData - '0';
  if (*(ptrData + 1) != ',')
  {
    ptrData++;
    app = (app * 10) + *ptrData - '0';
  }
#ifndef lang  
  lang = ((app & 1) ? true : false);
#endif
  if ((app & 0x40) != 0)
  {
    hexadecimal = true;
  }
  else
  {
    hexadecimal = false;
  }
  fcubesText(ptrData + 2, grpLen);
  databack(output);
}
#endif

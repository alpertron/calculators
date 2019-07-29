/*
This file is part of Alpertron Calculators.

Copyright 2015 Dario Alejandro Alpern

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
#include <string.h>
#include <stdio.h>
#include "bignbr.h"
#include "highlevel.h"
#include "batch.h"
static BigInteger value;
static BigInteger Base1, Base2, Base3, Base4;
static BigInteger P, Q, R, S, a, b;
static BigInteger P1, Q1, R1, S1;
static BigInteger tmpP1, tmpQ1, tmpR1, tmpS1;
static BigInteger toProcess;
static int groupLength;
static char *cube = "<span class=\"bigger\">³</span>";
extern int lang;
extern char hexadecimal;
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
  limb tmp;
  int index;
  limb *ptr1, *ptr2;
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
    tmp = *ptr1;
    *ptr1 = *ptr2;
    *ptr2 = tmp;
    ptr1++;
    ptr2++;
  }
}

static void EvaluateQuadraticPoly(BigInteger *pResult, BigInteger *pValue, int quad, int linear, int constant)
{
  multadd(pResult, quad, pValue, linear);
  // Multiply result by value.
  (void)BigIntMultiply(pResult, pValue, pResult);
  addbigint(pResult, constant);
}

static int fcubes(BigInteger *pArgument)
{
  int mod18, modulus, i, mod83, mask;
  int pow, exp;
  boolean converted = FALSE;
  CopyBigInt(&value, pArgument);
  // Compute argument mod 18.
  mod18 = getRemainder(pArgument, 18);
  if (mod18 == 4 || mod18 == 5 || mod18 == 13 || mod18 == 14)
  {
    return -1;  // Applet does not work if the number is congruent to 4 or 5 (mod 9)
  }
  if (mod18 == 16)
  {             // Change sign.
    converted = TRUE;
    BigIntNegate(&value, &value);
  }
  for (i = 0; i<(int)(sizeof(sums)/sizeof(sums[0])); i += 10)
  {
    modulus = sums[i];
    if ((getRemainder(&value, modulus) + modulus)% modulus == sums[i + 1])
    {
      break;
    }
  }
  if (i<(int)(sizeof(sums) / sizeof(sums[0])))
  {
    subtractdivide(&value, sums[i + 1], modulus);      // value <- (value-sums[i+1])/modulus
    multadd(&Base1, sums[i + 2], &value, sums[i + 3]); // Base1 <- sums[i+2]*value+sums[i+3]
    multadd(&Base2, sums[i + 4], &value, sums[i + 5]); // Base2 <- sums[i+4]*value+sums[i+5]
    multadd(&Base3, sums[i + 6], &value, sums[i + 7]); // Base3 <- sums[i+6]*value+sums[i+7]
    multadd(&Base4, sums[i + 8], &value, sums[i + 9]); // Base4 <- sums[i+8]*value+sums[i+9]
  }
  else if (getRemainder(&value, 54) == 2)
  {           // If value == 2 (mod 54)...
    subtractdivide(&value, 2, 54);   // value <- (value-2)/54
    EvaluateQuadraticPoly(&Base1, &value, 29484, 2211, 43);
    EvaluateQuadraticPoly(&Base2, &value, -29484, -2157, -41);
    EvaluateQuadraticPoly(&Base3, &value, 9828, 485, 4);
    EvaluateQuadraticPoly(&Base4, &value, -9828, -971, -22);
  }
  else if (getRemainder(&value, 83 * 108) == 83*46)
  {           // If value == 83*46 (mod 83*108)...
    subtractdivide(&value, 83*46, 83*108); // value <-(value - (83*46)) / (83*108)
    EvaluateQuadraticPoly(&Base1, &value, 29484, 25143, 5371);
    EvaluateQuadraticPoly(&Base2, &value, -29484, -25089, -5348);
    EvaluateQuadraticPoly(&Base3, &value, 9828, 8129, 1682);
    EvaluateQuadraticPoly(&Base4, &value, -9828, -8615, -1889);
  }
  else
  {
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
    P.limbs[1].x = S.limbs[1].x = 52;
    P.limbs[0].x = S.limbs[0].x = 819632865;
    Q.limbs[1].x = 3063;
    Q.limbs[0].x = 687764496;
    R.limbs[1].x = 0;
    R.limbs[0].x = 1923517596;
    P.nbrLimbs = Q.nbrLimbs = R.nbrLimbs = S.nbrLimbs = 3;
    P.sign = Q.sign = R.sign = S.sign = SIGN_NEGATIVE;
    P1.limbs[0].x = S1.limbs[0].x = 1;
    Q1.limbs[0].x = R1.limbs[0].x = 0;
    P1.nbrLimbs = Q1.nbrLimbs = R1.nbrLimbs = S1.nbrLimbs = 1;
    P1.sign = Q1.sign = R1.sign = S1.sign = SIGN_POSITIVE;
    mod83 = getRemainder(&value, 83);
    pow = 71;
    exp = 0;
    while (pow != mod83)
    {
      exp++;
      pow = (pow * 50) % 83;
    }
    if (exp > 82 / 2)
    {
      exp = 82 - exp;
      Q.sign = R.sign = SIGN_POSITIVE;
    }  // Now exp is in range 0-41.
    mask = 32;
    while (mask > 0)
    {
      // tmpP1 <- P1*P1 + Q1*R1
      // tmpQ1 <- (P1+S1) * Q1
      // tmpR1 <- (P1+S1) * R1
      // tmpS1 <- S1*S1 + Q1*R1
      // P1 <- tmpP1
      // Q1 <- tmpQ1
      // R1 <- tmpR1
      // S1 <- tmpS1
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
        // tmpP1 <- P*P1 + Q*R1
        // tmpQ1 <- P*Q1 + Q*S1
        // tmpR1 <- R*P1 + S*R1
        // tmpS1 <- R*Q1 + S*S1
        // P1 <- tmpP1
        // Q1 <- tmpQ1
        // R1 <- tmpR1
        // S1 <- tmpS1
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

  if (converted != 0)
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
  if (tmpQ1.nbrLimbs != 1 || tmpQ1.limbs[0].x != 0)
  {
    return 1;       // Result does not validate.
  }
  return 0;
}

void fcubesText(char *input, int groupLen)
{
  char *ptrOutput;
  if (valuesProcessed == 0)
  {
    groupLength = groupLen;
  }
  BatchProcessing(input, &toProcess, &ptrOutput);
  strcpy(ptrOutput, (lang ? "</p><p>" COPYRIGHT_SPANISH "</p>" :
    "</p><p>" COPYRIGHT_ENGLISH "</p>"));
}

void batchCubesCallback(char **pptrOutput)
{
  int result;
  char *ptrOutput = *pptrOutput;
  NumberLength = toProcess.nbrLimbs;
  result = fcubes(&toProcess);
  // Show the number to be decomposed into sum of cubes.
  strcpy(ptrOutput, "<p>");
  ptrOutput += strlen(ptrOutput);
  if (hexadecimal)
  {
    BigInteger2Hex(&toProcess, ptrOutput, groupLength);
  }
  else
  {
    BigInteger2Dec(&toProcess, ptrOutput, groupLength);
  }
  ptrOutput += strlen(ptrOutput);
  switch (result)
  {
  case -1:
    strcpy(ptrOutput, (lang==0?": This applet does not work if the number is congruent to 4 or 5 (mod 9)</p>":
      ": El applet no funciona si el número es congruente a 4 o 5 (mod 9)</p>"));
    *pptrOutput = ptrOutput + strlen(ptrOutput);
    return;
  case 1:
    strcpy(ptrOutput, (lang==0?": Internal error! Please send the number to the author of the applet.</p>":
      ": ¡Error interno!Por favor envíe este número al autor del applet.</p>"));
    *pptrOutput = ptrOutput + strlen(ptrOutput);
    return;
  case 2:
    strcpy(ptrOutput, (lang==0?": User stopped the calculation</p>":": El usuario detuvo el cálculo</p>"));
    *pptrOutput = ptrOutput + strlen(ptrOutput);
    return;
  }
  // Show decomposition in sum of 1, 2, 3 or 4 cubes.
  strcpy(ptrOutput, " = ");
  ptrOutput += strlen(ptrOutput);
  if (Base1.sign == SIGN_NEGATIVE)
  {
    *ptrOutput++ = '(';
  }
  if (hexadecimal)
  {
    BigInteger2Hex(&Base1, ptrOutput, groupLength);
  }
  else
  {
    BigInteger2Dec(&Base1, ptrOutput, groupLength);
  }
  ptrOutput += strlen(ptrOutput);
  if (Base1.sign == SIGN_NEGATIVE)
  {
    *ptrOutput++ = ')';
  }
  strcpy(ptrOutput, cube);
  ptrOutput += strlen(ptrOutput);
  if (Base2.nbrLimbs != 1 || Base2.limbs[0].x != 0)
  {
    strcpy(ptrOutput, " + ");
    ptrOutput += strlen(ptrOutput);
    if (Base2.sign == SIGN_NEGATIVE)
    {
      *ptrOutput++ = '(';
    }
    if (hexadecimal)
    {
      BigInteger2Hex(&Base2, ptrOutput, groupLength);
    }
    else
    {
      BigInteger2Dec(&Base2, ptrOutput, groupLength);
    }
    ptrOutput += strlen(ptrOutput);
    if (Base2.sign == SIGN_NEGATIVE)
    {
      *ptrOutput++ = ')';
    }
    strcpy(ptrOutput, cube);
    ptrOutput += strlen(ptrOutput);
  }
  if (Base3.nbrLimbs != 1 || Base3.limbs[0].x != 0)
  {
    strcpy(ptrOutput, " + ");
    ptrOutput += strlen(ptrOutput);
    if (Base3.sign == SIGN_NEGATIVE)
    {
      *ptrOutput++ = '(';
    }
    if (hexadecimal)
    {
      BigInteger2Hex(&Base3, ptrOutput, groupLength);
    }
    else
    {
      BigInteger2Dec(&Base3, ptrOutput, groupLength);
    }
    ptrOutput += strlen(ptrOutput);
    if (Base3.sign == SIGN_NEGATIVE)
    {
      *ptrOutput++ = ')';
    }
    strcpy(ptrOutput, cube);
    ptrOutput += strlen(ptrOutput);
  }
  if (Base4.nbrLimbs != 1 || Base4.limbs[0].x != 0)
  {
    strcpy(ptrOutput, " + ");
    ptrOutput += strlen(ptrOutput);
    if (Base4.sign == SIGN_NEGATIVE)
    {
      *ptrOutput++ = '(';
    }
    if (hexadecimal)
    {
      BigInteger2Hex(&Base4, ptrOutput, groupLength);
    }
    else
    {
      BigInteger2Dec(&Base4, ptrOutput, groupLength);
    }
    ptrOutput += strlen(ptrOutput);
    if (Base4.sign == SIGN_NEGATIVE)
    {
      *ptrOutput++ = ')';
    }
    strcpy(ptrOutput, cube);
    ptrOutput += strlen(ptrOutput);
  }
  *pptrOutput = ptrOutput;
}

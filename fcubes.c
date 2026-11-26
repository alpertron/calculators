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
#include <assert.h>
#include "string/strings.h"
#include "bignbr.h"
#include "expression.h"
#include "highlevel.h"
#include "showtime.h"
#include "batch.h"
#include "copyStr.h"
#include "output.h"
#include "commonstruc.h"

#ifndef __ANDROID__
union uCommon common;
#endif

extern bool hexadecimal;
static const char *cube = "<span class=\"bigger\">³</span>";
static void batchCubesCallback(char** pptrOutput, int type);
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
  810, 164, 3, -26, -3, 1, 5, -21, -5, 30,
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
  (void)BigIntMultiply(&common.misc.Base1, &common.misc.Base1, &common.misc.tmpQ1);
  (void)BigIntMultiply(&common.misc.tmpQ1, &common.misc.Base1, &common.misc.tmpP1);

  (void)BigIntMultiply(&common.misc.Base2, &common.misc.Base2, &common.misc.tmpQ1);
  (void)BigIntMultiply(&common.misc.tmpQ1, &common.misc.Base2, &common.misc.tmpQ1);
  BigIntAdd(&common.misc.tmpP1, &common.misc.tmpQ1, &common.misc.tmpP1);

  (void)BigIntMultiply(&common.misc.Base3, &common.misc.Base3, &common.misc.tmpQ1);
  (void)BigIntMultiply(&common.misc.tmpQ1, &common.misc.Base3, &common.misc.tmpQ1);
  BigIntAdd(&common.misc.tmpP1, &common.misc.tmpQ1, &common.misc.tmpP1);

  (void)BigIntMultiply(&common.misc.Base4, &common.misc.Base4, &common.misc.tmpQ1);
  (void)BigIntMultiply(&common.misc.tmpQ1, &common.misc.Base4, &common.misc.tmpQ1);
  BigIntAdd(&common.misc.tmpP1, &common.misc.tmpQ1, &common.misc.tmpP1);
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
  assert((nbrLimbs1 > 0) && (nbrLimbs2 > 0));
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
  common.misc.P.limbs[1].x = 52;
  common.misc.P.limbs[0].x = 819632865;
  common.misc.Q.limbs[1].x = 3063;
  common.misc.Q.limbs[0].x = 687764496;
  common.misc.Q.nbrLimbs = 2;
  common.misc.P.nbrLimbs = 2;
  common.misc.Q.sign = SIGN_NEGATIVE;
  common.misc.P.sign = SIGN_NEGATIVE;
  CopyBigInt(&common.misc.S, &common.misc.P);
  intToBigInteger(&common.misc.R, -1923517596);
  intToBigInteger(&common.misc.P1, 1);
  intToBigInteger(&common.misc.Q1, 0);
  intToBigInteger(&common.misc.R1, 0);
  intToBigInteger(&common.misc.S1, 1);
  mod83 = getRemainder(&common.misc.value, 83);
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
    common.misc.Q.sign = SIGN_POSITIVE;
    common.misc.R.sign = SIGN_POSITIVE;
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
    (void)BigIntMultiply(&common.misc.P1, &common.misc.P1, &common.misc.tmpP1);
    (void)BigIntMultiply(&common.misc.Q1, &common.misc.R1, &common.misc.tmpQ1);
    (void)BigIntMultiply(&common.misc.S1, &common.misc.S1, &common.misc.tmpS1);
    BigIntAdd(&common.misc.P1, &common.misc.S1, &common.misc.tmpR1);
    BigIntAdd(&common.misc.tmpP1, &common.misc.tmpQ1, &common.misc.P1);
    BigIntAdd(&common.misc.tmpS1, &common.misc.tmpQ1, &common.misc.S1);
    (void)BigIntMultiply(&common.misc.tmpR1, &common.misc.Q1, &common.misc.Q1);
    (void)BigIntMultiply(&common.misc.tmpR1, &common.misc.R1, &common.misc.R1);
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
      (void)BigIntMultiply(&common.misc.P, &common.misc.P1, &common.misc.tmpP1);
      (void)BigIntMultiply(&common.misc.Q, &common.misc.R1, &common.misc.tmpQ1);
      (void)BigIntMultiply(&common.misc.R, &common.misc.P1, &common.misc.tmpR1);
      (void)BigIntMultiply(&common.misc.S, &common.misc.R1, &common.misc.tmpS1);
      BigIntAdd(&common.misc.tmpP1, &common.misc.tmpQ1, &common.misc.P1);
      BigIntAdd(&common.misc.tmpR1, &common.misc.tmpS1, &common.misc.R1);
      (void)BigIntMultiply(&common.misc.P, &common.misc.Q1, &common.misc.tmpP1);
      (void)BigIntMultiply(&common.misc.Q, &common.misc.S1, &common.misc.tmpQ1);
      (void)BigIntMultiply(&common.misc.R, &common.misc.Q1, &common.misc.tmpR1);
      (void)BigIntMultiply(&common.misc.S, &common.misc.S1, &common.misc.tmpS1);
      BigIntAdd(&common.misc.tmpP1, &common.misc.tmpQ1, &common.misc.Q1);
      BigIntAdd(&common.misc.tmpR1, &common.misc.tmpS1, &common.misc.S1);
    }
    mask >>= 1;
  }
  addmult(&common.misc.a, &common.misc.P1, -3041, &common.misc.Q1, -52);   // a <- -3041*P1 - 52*Q1
  addmult(&common.misc.b, &common.misc.R1, -3041, &common.misc.S1, -52);   // b <- -3041*R1 - 52*S1
  addmult(&common.misc.Base1, &common.misc.a, 27, &common.misc.b, -928);   // Base1 <- 27*a - 928*b
  addmult(&common.misc.Base2, &common.misc.a, -9, &common.misc.b, -602);   // Base2 <- -9*a - 602*b
  addmult(&common.misc.Base3, &common.misc.a, 25, &common.misc.b, -2937);  // Base3 <- 25*a - 2937*b
  addmult(&common.misc.Base4, &common.misc.a, -19, &common.misc.b, 2746);  // Base4 <- -19*a - 2746*b
  // a <- (value - Base1^3 - Base2^3 - Base3^3 - Base4^3)/(18*83)
  getSumOfCubes();  // tmpP1 = Base1^3 + Base2^3 + Base3^3 + Base4^3
  BigIntSubt(&common.misc.value, &common.misc.tmpP1, &common.misc.a);
  subtractdivide(&common.misc.a, 0, 18 * 83);      // Divide a by 18*83.
  multint(&common.misc.tmpP1, &common.misc.a, 10);             // Base1 <- Base1 + 10*a
  BigIntAdd(&common.misc.tmpP1, &common.misc.Base1, &common.misc.Base1);
  multint(&common.misc.tmpP1, &common.misc.a, -19);            // Base2 <- Base2 - 19*a
  BigIntAdd(&common.misc.tmpP1, &common.misc.Base2, &common.misc.Base2);
  multint(&common.misc.tmpP1, &common.misc.a, -24);            // Base3 <- Base3 - 24*a
  BigIntAdd(&common.misc.tmpP1, &common.misc.Base3, &common.misc.Base3);
  multint(&common.misc.tmpP1, &common.misc.a, 27);             // Base4 <- Base4 + 27*a
  BigIntAdd(&common.misc.tmpP1, &common.misc.Base4, &common.misc.Base4);
}

static void showLinear(char **pptrOutput, int linear, int indep)
{
  if (linear == 1)
  {
    copyStr(pptrOutput, "x");
  }
  else if (linear == -1)
  {
    copyStr(pptrOutput, "&minus;x");
  }
  else
  {
    if (linear >= 0)
    {
      int2dec(pptrOutput, linear);
    }
    else
    {
      copyStr(pptrOutput, "&minus;");
      int2dec(pptrOutput, -linear);
    }
    copyStr(pptrOutput, "&#8290;x");
  }
  if (indep < 0)
  {
    copyStr(pptrOutput, " &minus; ");
    int2dec(pptrOutput, -indep);
  }
  else if (indep > 0)
  {
    copyStr(pptrOutput, " + ");
    int2dec(pptrOutput, indep);
  }
}

static void showCubeOfLinear(char** pptrOutput, int linear, int indep)
{
  copyStr(pptrOutput, "(");
  showLinear(pptrOutput, linear, indep);
  copyStr(pptrOutput, ")³");
}

static void showCubeOfQuadratic(char** pptrOutput, int quad, int linear, int indep)
{
  copyStr(pptrOutput, "(");
  if (quad >= 0)
  {
    int2dec(pptrOutput, quad);
  }
  else
  {
    copyStr(pptrOutput, "&minus;");
    int2dec(pptrOutput, -quad);
  }
  if (linear >= 0)
  {
    copyStr(pptrOutput, " + ");
    int2dec(pptrOutput, linear);
  }
  else
  {
    copyStr(pptrOutput, " &minus; ");
    int2dec(pptrOutput, -linear);
  }
  if (indep >= 0)
  {
    copyStr(pptrOutput, " + ");
    int2dec(pptrOutput, indep);
  }
  else
  {
    copyStr(pptrOutput, " &minus; ");
    int2dec(pptrOutput, -indep);
  }
  copyStr(pptrOutput, ")³");
}

static int fcubes(char **pptrOutput, const BigInteger *pArgument)
{
  int mod18;
  int i;
  bool converted = false;
  CopyBigInt(&common.misc.value, pArgument);
  // Compute argument mod 18.
  mod18 = getRemainder(pArgument, 18);
  if ((mod18 == 4) || (mod18 == 5) || (mod18 == 13) || (mod18 == 14))
  {
    return -1;  // Applet does not work if the number is congruent to 4 or 5 (mod 9)
  }
  if (mod18 == 16)
  {             // Change sign.
    converted = true;
    BigIntNegate(&common.misc.value, &common.misc.value);
  }
  for (i = ((int)sizeof(sums) / (int)sizeof(sums[0]))-10; i>=0; i -= 10)
  {
    int modulus = sums[i];
    if (((getRemainder(&common.misc.value, modulus) + modulus)% modulus) == sums[i + 1])
    {
      subtractdivide(&common.misc.value, sums[i + 1], modulus);      // value <- (value-sums[i+1])/modulus
      multadd(&common.misc.Base1, sums[i + 2], &common.misc.value, sums[i + 3]); // Base1 <- sums[i+2]*value+sums[i+3]
      multadd(&common.misc.Base2, sums[i + 4], &common.misc.value, sums[i + 5]); // Base2 <- sums[i+4]*value+sums[i+5]
      multadd(&common.misc.Base3, sums[i + 6], &common.misc.value, sums[i + 7]); // Base3 <- sums[i+6]*value+sums[i+7]
      multadd(&common.misc.Base4, sums[i + 8], &common.misc.value, sums[i + 9]); // Base4 <- sums[i+8]*value+sums[i+9]
      formatString(pptrOutput, "<p>$1s ", LITERAL_FCUBES1);  // Using the formula
      showLinear(pptrOutput, modulus, sums[i + 1]);
      copyStr(pptrOutput, " = ");
      showCubeOfLinear(pptrOutput, sums[i + 2], sums[i + 3]);
      copyStr(pptrOutput, " + ");
      showCubeOfLinear(pptrOutput, sums[i + 4], sums[i + 5]);
      copyStr(pptrOutput, " + ");
      showCubeOfLinear(pptrOutput, sums[i + 6], sums[i + 7]);
      copyStr(pptrOutput, " + ");
      showCubeOfLinear(pptrOutput, sums[i + 8], sums[i + 9]);
      copyStr(pptrOutput, ":</p>");
      break;
    }
  }
  if (i < 0)
  {
    if (getRemainder(&common.misc.value, 54) == 2)
    {           // If value == 2 (mod 54)...
      subtractdivide(&common.misc.value, 2, 54);   // value <- (value-2)/54
      EvaluateQuadraticPoly(&common.misc.Base1, &common.misc.value, 29484, 2211, 43);
      EvaluateQuadraticPoly(&common.misc.Base2, &common.misc.value, -29484, -2157, -41);
      EvaluateQuadraticPoly(&common.misc.Base3, &common.misc.value, 9828, 485, 4);
      EvaluateQuadraticPoly(&common.misc.Base4, &common.misc.value, -9828, -971, -22);
      copyStr(pptrOutput, LITERAL_FCUBES1);
      copyStr(pptrOutput, " ");
      showLinear(pptrOutput, 54, 2);
      copyStr(pptrOutput, " = ");
      showCubeOfQuadratic(pptrOutput, 29484, 2211, 43);
      copyStr(pptrOutput, " + ");
      showCubeOfQuadratic(pptrOutput, -29484, -2157, -41);
      copyStr(pptrOutput, " + ");
      showCubeOfQuadratic(pptrOutput, 9828, 485, 4);
      copyStr(pptrOutput, " + ");
      showCubeOfQuadratic(pptrOutput, -9828, -971, -22);
      copyStr(pptrOutput, ":</p>");
    }
    else if (getRemainder(&common.misc.value, 83 * 108) == (83 * 46))
    {           // If value == 83*46 (mod 83*108)...
      subtractdivide(&common.misc.value, 83 * 46, (83 * 108)); // value <-(value - (83*46)) / (83*108)
      EvaluateQuadraticPoly(&common.misc.Base1, &common.misc.value, 29484, 25143, 5371);
      EvaluateQuadraticPoly(&common.misc.Base2, &common.misc.value, -29484, -25089, -5348);
      EvaluateQuadraticPoly(&common.misc.Base3, &common.misc.value, 9828, 8129, 1682);
      EvaluateQuadraticPoly(&common.misc.Base4, &common.misc.value, -9828, -8615, -1889);
      copyStr(pptrOutput, LITERAL_FCUBES1);
      copyStr(pptrOutput, " ");
      showLinear(pptrOutput, 83*108, 83*46);
      copyStr(pptrOutput, " = ");
      showCubeOfQuadratic(pptrOutput, 29484, 25143, 5371);
      copyStr(pptrOutput, " + ");
      showCubeOfQuadratic(pptrOutput, -29484, -25089, -5348);
      copyStr(pptrOutput, " + ");
      showCubeOfQuadratic(pptrOutput, 9828, 8129, 1682);
      copyStr(pptrOutput, " + ");
      showCubeOfQuadratic(pptrOutput, -9828, -8615, -1889);
      copyStr(pptrOutput, ":</p>");
    }
    else
    {
      copyStr(pptrOutput, LITERAL_FCUBES2);
      Demjanenko();
    }
  }
  if (converted)
  {
    BigIntNegate(&common.misc.Base1, &common.misc.Base1);
    BigIntNegate(&common.misc.Base2, &common.misc.Base2);
    BigIntNegate(&common.misc.Base3, &common.misc.Base3);
    BigIntNegate(&common.misc.Base4, &common.misc.Base4);
  }

      // Sort cubes
  SortBigIntegers(&common.misc.Base1, &common.misc.Base2);
  SortBigIntegers(&common.misc.Base1, &common.misc.Base3);
  SortBigIntegers(&common.misc.Base1, &common.misc.Base4);
  SortBigIntegers(&common.misc.Base2, &common.misc.Base3);
  SortBigIntegers(&common.misc.Base2, &common.misc.Base4);
  SortBigIntegers(&common.misc.Base3, &common.misc.Base4);
  return 0;
}

void fcubesText(char *input, int grpLen)
{
  char *ptrOutput;
  if (valuesProcessed == 0)
  {
    common.misc.groupLength = grpLen;
  }
  (void)BatchProcessing(input, &common.misc.toProcess, &ptrOutput, NULL, batchCubesCallback);
#ifdef __EMSCRIPTEN__
  int elapsedTime = (int)(tenths() - originalTenthSecond);
  copyStr(&ptrOutput, "<p>");
  formatString(&ptrOutput, LITERAL_FCUBES3, elapsedTime);
  copyStr(&ptrOutput, "</p>");
#endif
  copyStr(&ptrOutput, "<p>");
  showCopyright(&ptrOutput);
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
    BigInteger2Hex(&ptrOutput, pBase, common.misc.groupLength);
  }
  else
  {
    BigInteger2Dec(&ptrOutput, pBase, common.misc.groupLength);
  }
  if (pBase->sign == SIGN_NEGATIVE)
  {
    *ptrOutput = ')';
    ptrOutput++;
  }
  copyStr(&ptrOutput, cube);
  *pptrOutput = ptrOutput;
}

static void batchCubesCallback(char **pptrOutput, int type)
{
  int result;
  char *ptrOutput = *pptrOutput;
  NumberLength = common.misc.toProcess.nbrLimbs;
  if ((type == BATCH_NO_PROCESS_DEC) || (type == BATCH_NO_PROCESS_HEX))
  {         // Do not compute sum of cubes.
    if (hexadecimal)
    {
      BigInteger2Hex(&ptrOutput, &common.misc.toProcess, common.misc.groupLength);
    }
    else
    {
      BigInteger2Dec(&ptrOutput, &common.misc.toProcess, common.misc.groupLength);
    }
    *pptrOutput = ptrOutput;
    return;
  }
  result = fcubes(&ptrOutput, &common.misc.toProcess);
  // Show the number to be decomposed into sum of cubes.
  if (type == BATCH_NO_QUOTE)
  {
    copyStr(&ptrOutput, "<p>");
    if (hexadecimal)
    {
      BigInteger2Hex(&ptrOutput, &common.misc.toProcess, common.misc.groupLength);
    }
    else
    {
      BigInteger2Dec(&ptrOutput, &common.misc.toProcess, common.misc.groupLength);
    }
  }
  if ((type == BATCH_NO_QUOTE) && (result != 0))
  {
    *ptrOutput = ':';
    ptrOutput++;
    *ptrOutput = ' ';
    ptrOutput++;
  }
  switch (result)
  {
  case -1:
    copyStr(&ptrOutput, LITERAL_FCUBES_BATCH_CBACK1);
    break;
  case 2:
    copyStr(&ptrOutput, LITERAL_FCUBES_BATCH_CBACK2);
    break;
  default:
    break;
  }
  if (result != 0)
  {
    if (type == BATCH_NO_QUOTE)
    {
      copyStr(&ptrOutput, "</p>");
    }
    *pptrOutput = ptrOutput;
    return;
  }
  // Show decomposition in sum of 1, 2, 3 or 4 cubes.
  if (type == BATCH_NO_QUOTE)
  {
    copyStr(&ptrOutput, " = ");
  }
  showCube(&ptrOutput, &common.misc.Base1);
  if (!BigIntIsZero(&common.misc.Base2))
  {
    copyStr(&ptrOutput, " + ");
    showCube(&ptrOutput, &common.misc.Base2);
  }
  if (!BigIntIsZero(&common.misc.Base3))
  {
    copyStr(&ptrOutput, " + ");
    showCube(&ptrOutput, &common.misc.Base3);
  }
  if (!BigIntIsZero(&common.misc.Base4))
  {
    copyStr(&ptrOutput, " + ");
    showCube(&ptrOutput, &common.misc.Base4);
  }
  if (type == BATCH_NO_QUOTE)
  {
    copyStr(&ptrOutput, "</p>");
  }
  *pptrOutput = ptrOutput;
}

#if defined(__EMSCRIPTEN__) && !defined(_MSC_VER)
#ifdef __ANDROID__
EXTERNALIZE void doWorkFCubes(void)
#else
EXTERNALIZE void doWork(void)
#endif
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
  while (*(ptrData + 1) != ',')
  {
    ptrData++;
    app = (app * 10) + *ptrData - '0';
  }
  if ((app & 0x40) != 0)
  {
    hexadecimal = true;
  }
  else
  {
    hexadecimal = false;
  }
  groupSeparator = app >> 8;
  fcubesText(ptrData + 2, grpLen);
  databack(output);
}
#endif

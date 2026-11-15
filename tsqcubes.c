//
// This file is part of Alpertron Calculators.
//
// Copyright 2022 Dario Alejandro Alpern
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
#include "string/strings.h"
#include "bignbr.h"
#include "expression.h"
#include "highlevel.h"
#include "showtime.h"
#include "batch.h"
#include "tsquares.h"
#include "output.h"
#include "copyStr.h"
#include "commonstruc.h"

static void batchSqCubesCallback(char** pptrOutput, int type);
static int Exponent = 3;
static char bigExponText[20];

static bool checkSumOfTwoSquares(const BigInteger *pArgument, int expon)
{
  int lenBytes;
  int numLimbs;
  int mod8;
  common.misc.power4 = 0;
  (void)BigIntPowerIntExp(&common.misc.Base3, expon, &common.misc.powerN);
  BigIntSubt(pArgument, &common.misc.powerN, &common.misc.powerN);
  // Test whether powerN is a perfect square.
  numLimbs = common.misc.powerN.nbrLimbs;
  lenBytes = numLimbs * (int)sizeof(limb);
  (void)memcpy(common.misc.Base2.limbs, common.misc.powerN.limbs, lenBytes);
  squareRoot(common.misc.Base2.limbs, common.misc.Base1.limbs,
    numLimbs, &common.misc.Base1.nbrLimbs);
  multiply(common.misc.Base1.limbs, common.misc.Base1.limbs,
    common.misc.Base2.limbs, common.misc.Base1.nbrLimbs, NULL);
  if (memcmp(common.misc.Base2.limbs, common.misc.powerN.limbs, lenBytes) == 0)
  {     // powerN can be expressed as a square of Base1.
    intToBigInteger(&common.misc.Base2, 0);
    return true;
  }
  DivideBigNbrByMaxPowerOf4(&common.misc.power4, common.misc.powerN.limbs, &common.misc.powerN.nbrLimbs);
  mod8 = common.misc.powerN.limbs[0].x & 7;
  if ((mod8 > 2) && (mod8 != 5))
  {     // If number mod 8 is not 1, 2 or 5,
        // it cannot be expressed as sum of two squares.
    return false;
  }
  if ((common.misc.powerN.nbrLimbs == 1) && (common.misc.powerN.limbs[0].x == 2))
  {
    intToBigInteger(&common.misc.Base1, 1);
    intToBigInteger(&common.misc.Base2, 1);
    return true;
  }
  common.misc.nbrLimbsP = common.misc.powerN.nbrLimbs;
  common.misc.nbrLimbs = common.misc.nbrLimbsP;
  (void)memcpy(common.misc.valueP, common.misc.powerN.limbs, common.misc.nbrLimbsP * sizeof(limb));
  FillSieveArray(common.misc.valueP, common.misc.sieve);
  if (!isSumOfTwoSquares())
  {
    return false;
  }
  CopyBigInt(&common.misc.Base1, &common.misc.biMult1);
  CopyBigInt(&common.misc.Base2, &common.misc.biMult2);
  return true;
}

// Compute Base1, Base2 and Base3 such that 
// Base1^2 + Base2^2 + Base3^exp = *pArgument.
// Compute k = floor(root(pArgument)).
// Then for i = k, k - 1, k - 2, ... try to represent
// pArgument - i^exp as a sum of two squares.
static int tsqcubes(const BigInteger *pArgument, int expon)
{
  // Compute Base3 <- floor(root(pArgument)).
  (void)BigIntRoot(pArgument, &common.misc.Base3, expon);
  InitSieveArray();
  while (!checkSumOfTwoSquares(pArgument, expon))
  { // If we cannot express powerN as a sum of two squares,
    // try previous value of Base3.
    addbigint(&common.misc.Base3, -1); // Use previous value of Base3.
  }
  common.misc.Base1.sign = SIGN_POSITIVE;
  common.misc.Base2.sign = SIGN_POSITIVE;
  (void)BigIntMultiplyPower2(&common.misc.Base1, common.misc.power4);
  (void)BigIntMultiplyPower2(&common.misc.Base2, common.misc.power4);

  // Validate

  (void)BigIntMultiply(&common.misc.Base1, &common.misc.Base1, &common.misc.valueN);
  (void)BigIntMultiply(&common.misc.Base2, &common.misc.Base2, &common.misc.powerN);
  BigIntAdd(&common.misc.powerN, &common.misc.valueN, &common.misc.valueN);
  (void)BigIntPowerIntExp(&common.misc.Base3, expon, &common.misc.powerN);
  BigIntAdd(&common.misc.powerN, &common.misc.valueN, &common.misc.valueN);
  if (!TestBigNbrEqual(&common.misc.valueN, pArgument))
  {
    return 1;       // Result does not validate.
  }
  return 0;
}

void tsqcubesText(char *input, int grpLen, int expon)
{
  char *ptrOutput;
  if (valuesProcessed == 0)
  {
    common.misc.groupLength = grpLen;
  }
  if (expon > 0)
  {
    Exponent = expon;
  }
  if (Exponent == 3)
  {
    common.misc.powerStr = "<span class=\"bigger\">³</span>";
  }
  else
  {
    char* ptrBigExponText = bigExponText;
    copyStr(&ptrBigExponText, "<sup>");
    int2dec(&ptrBigExponText, Exponent);
    copyStr(&ptrBigExponText, "</sup>");
    *ptrBigExponText = 0;
    common.misc.powerStr = bigExponText;
  }
  (void)BatchProcessing(input, &common.misc.toProcess, &ptrOutput, NULL, batchSqCubesCallback);
  showElapsedTime(&ptrOutput);
  copyStr(&ptrOutput, "<p>");
  showCopyright(&ptrOutput);
  copyStr(&ptrOutput, "</p>");
}

// Show cube number. Use parentheses for negative numbers.
static void showTerm(char** pptrOutput, const BigInteger* pBase, const char *expon)
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
  copyStr(&ptrOutput, expon);
  *pptrOutput = ptrOutput;
}

static void batchSqCubesCallback(char **pptrOutput, int type)
{
  int result;
  char *ptrOutput = *pptrOutput;
  if ((type == BATCH_NO_PROCESS_DEC) || (type == BATCH_NO_PROCESS_HEX))
  {         // Do not compute sum of squares.
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
  NumberLength = common.misc.toProcess.nbrLimbs;
  result = tsqcubes(&common.misc.toProcess, Exponent);
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
  if (result == 1)
  {
    // Internal error! Please send the number to the author of the calculator.
    copyStr(&ptrOutput, LITERAL_BATCH_SQCUBE_CBACK);
    copyStr(&ptrOutput, "</p>");
    *pptrOutput = ptrOutput;
    return;
  }
  // Use powerN as a temporary variable.
  BigIntSubt(&common.misc.Base1, &common.misc.Base2, &common.misc.powerN);
  if (common.misc.powerN.sign == SIGN_NEGATIVE)
  {    // Swap Base1 and Base2 so the first is always greater than the second.
    CopyBigInt(&common.misc.powerN, &common.misc.Base1);
    CopyBigInt(&common.misc.Base1, &common.misc.Base2);
    CopyBigInt(&common.misc.Base2, &common.misc.powerN);
  }
  // Show decomposition in sum of up to two squares and a cube.
  if (type == BATCH_NO_QUOTE)
  {
    copyStr(&ptrOutput, " = ");
  }
  if (BigIntIsZero(&common.misc.Base2) && BigIntIsZero(&common.misc.Base1))
  {
    showTerm(&ptrOutput, &common.misc.Base3, common.misc.powerStr);
  }
  else
  {
    showTerm(&ptrOutput, &common.misc.Base1, square);
    if (!BigIntIsZero(&common.misc.Base2))
    {
      copyStr(&ptrOutput, " + ");
      showTerm(&ptrOutput, &common.misc.Base2, square);
    }
    if (!BigIntIsZero(&common.misc.Base3))
    {
      copyStr(&ptrOutput, " + ");
      showTerm(&ptrOutput, &common.misc.Base3, common.misc.powerStr);
    }
  }
  if (type == BATCH_NO_QUOTE)
  {
    copyStr(&ptrOutput, "</p>");
  }
  *pptrOutput = ptrOutput;
}

#if defined(__EMSCRIPTEN__) && !defined(_MSC_VER)
#ifdef __ANDROID__
EXTERNALIZE void doWorkTSqCubes(void)
#else
EXTERNALIZE void doWork(void)
#endif
{
  int expon;
  int app;
  int grpLen = 0;
  char* ptrData = inputString;
  originalTenthSecond = tenths();
  if (*ptrData == 'C')
  {    // User pressed Continue button.
    tsqcubesText(NULL, 0, -1); // Routine does not use parameters in this case.
    databack(output);
    return;
  }
  expon = 0;
  while (*ptrData != ',')
  {
    expon = (expon * 10) + (*ptrData - '0');
    ptrData++;
  }
  ptrData++;        // Skip comma.
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
  if ((app & 0x40) != 0)
  {
    hexadecimal = true;
  }
  else
  {
    hexadecimal = false;
  }
  tsqcubesText(ptrData + 2, grpLen, expon);
  databack(output);
}
#endif

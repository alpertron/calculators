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

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "string/strings.h"
#include "bignbr.h"
#include "output.h"
#include "expression.h"
#include "highlevel.h"
#include "showtime.h"
#include "batch.h"
#include "tsquares.h"
#include "copyStr.h"
#include "commonstruc.h"

static void batchSquaresCallback(char** pptrOutput, int type);

 // If Mult1 < Mult2, exchange both numbers.
static void SortBigNbrs(BigInteger *pbiMult1, BigInteger *pbiMult2)
{   // biFirstTerm not used at this time. Use it as temp variable.
  BigIntSubt(pbiMult1, pbiMult2, &common.misc.biFirstTerm);
  if (common.misc.biFirstTerm.sign == SIGN_NEGATIVE)
  {                // mult1 < mult2 -> exchange them.
    CopyBigInt(&common.misc.biFirstTerm, pbiMult1);
    CopyBigInt(pbiMult1, pbiMult2);
    CopyBigInt(pbiMult2, &common.misc.biFirstTerm);
  }
}

static void UpdateSieveArray(void)
{
  for (int index = 0; index < MAX_SIEVE; index++)
  {
    if (common.misc.sieve[index] >= 0)
    {
      int modulus = (2 * index) + 3;
      common.misc.sieve[index] = (common.misc.sieve[index] + common.misc.delta[index]) % modulus;
      common.misc.delta[index] -= 2;
      if (common.misc.delta[index] < 0)
      {
        common.misc.delta[index] += modulus;
      }
    }
  }
}

// Try to decompose a number as a sum of three or four squares.
// Mult1 = sqrt(toProcess).
// common.misc.biMult4 = toProcess - biFirstTerm^2.
// Use common.misc.biMult1, common.misc.biMult2 and common.misc.biMult3 as temporary storage.

static void SumOfSquaresNumber(void)
{
  common.misc.attempts = 0;
  FillSieveArray(common.misc.toProcess.limbs, common.misc.sieve);
  CopyBigInt(&common.misc.biMult3, &common.misc.toProcess);
  int nbrBytes = common.misc.toProcess.nbrLimbs * (int)sizeof(limb);
  if ((common.misc.toProcess.limbs[0].x & 7) != 7)
  {              // n!=7 (mod 8) => Sum of three squares
                 // Compute common.misc.biMult4 as a sum of two squares.
    CopyBigInt(&common.misc.biSecondTerm, &common.misc.biFirstTerm);
    intToBigInteger(&common.misc.biFirstTerm, 0);
    (void)memcpy(common.misc.valueP, common.misc.toProcess.limbs, nbrBytes);
    common.misc.nbrLimbsP = common.misc.toProcess.nbrLimbs;
    if (isSumOfTwoSquares())
    {            // Number is a sum of two squares: common.misc.biMult1^2 + common.misc.biMult2^2.
      intToBigInteger(&common.misc.biFirstTerm, 0);
      intToBigInteger(&common.misc.biSecondTerm, 0);
      CopyBigInt(&common.misc.biMult3, &common.misc.biMult1);
      CopyBigInt(&common.misc.biMult4, &common.misc.biMult2);
      return;
    }
  }
  else
  {              // n==7 (mod 8) => Sum of four squares
                 // Compute common.misc.biMult4 as a sum of three squares.
    for (;;)
    {
      int pwr4;
      CopyBigInt(&common.misc.biMult3, &common.misc.biMult4);
      DivideBigNbrByMaxPowerOf4(&pwr4, common.misc.biMult3.limbs, &common.misc.biMult3.nbrLimbs);
      if ((common.misc.biMult3.limbs[0].x & 7) != 7)
      {          // common.misc.biMult4 is a sum of three squares. Compute them.
        break;
      }
      addbigint(&common.misc.biFirstTerm, -1);  // Next candidate for largest term.
      BigIntAdd(&common.misc.biMult4, &common.misc.biFirstTerm, &common.misc.biMult4);
      BigIntAdd(&common.misc.biMult4, &common.misc.biFirstTerm, &common.misc.biMult4);
      addbigint(&common.misc.biMult4, 1);
    }
    squareRoot(common.misc.biMult4.limbs, common.misc.biSecondTerm.limbs, common.misc.biMult4.nbrLimbs,
      &common.misc.biSecondTerm.nbrLimbs);
    (void)BigIntMultiply(&common.misc.biSecondTerm, &common.misc.biSecondTerm, &common.misc.biMult3);
    BigIntSubt(&common.misc.biMult4, &common.misc.biMult3, &common.misc.biMult4);
  }
  common.misc.nbrLimbs = common.misc.biMult4.nbrLimbs;
  FillSieveArray(common.misc.biMult4.limbs, common.misc.sieve);
  multint(&common.misc.biMult3, &common.misc.biSecondTerm, 2);
  addbigint(&common.misc.biMult3, -1);
  common.misc.nbrLimbs = common.misc.biMult3.nbrLimbs;
  FillSieveArray(common.misc.biMult3.limbs, common.misc.delta);
  // Compute common.misc.biMult4 as a sum of two squares:
  // If the result is not the product of a power of 2, small primes
  // of the form 4k+1 and squares of primes of the form (4k+3)^2
  // and a (big) prime, it is not a sum of two squares (or it is
  // very difficult to compute it).
  // In this case, decrement Mult2 and try again.
  for (;;)
  {
    // p should be the product of power of 2,
    // powers of small primes of form 4k+1 and
    // powers of squares of small primes of form 4k+3.
    nbrBytes = common.misc.biMult4.nbrLimbs * (int)sizeof(limb);
    (void)memcpy(common.misc.valueP, common.misc.biMult4.limbs, nbrBytes);
    common.misc.nbrLimbsP = common.misc.biMult4.nbrLimbs;
    if (isSumOfTwoSquares())
    {   // Found common.misc.biMult1^2 + common.misc.biMult2^2.
      CopyBigInt(&common.misc.biMult3, &common.misc.biMult1);
      CopyBigInt(&common.misc.biMult4, &common.misc.biMult2);
      break;
    }
    addbigint(&common.misc.biSecondTerm, -1);
    BigIntAdd(&common.misc.biMult4, &common.misc.biSecondTerm, &common.misc.biMult4);
    BigIntAdd(&common.misc.biMult4, &common.misc.biSecondTerm, &common.misc.biMult4);
    addbigint(&common.misc.biMult4, 1);
    UpdateSieveArray();
  }
}

// Variable to split in up to four squares: toProcess.
int fsquares(void)
{
#ifdef __EMSCRIPTEN__
  char *ptrOutput;
#endif
  // Use biSecondTerm as temporary variable to display the original number.
  CopyBigInt(&common.misc.biSecondTerm, &common.misc.toProcess);
  DivideBigNbrByMaxPowerOf4(&common.misc.power4, common.misc.toProcess.limbs, &common.misc.toProcess.nbrLimbs);
  // Get biFirstTerm <- square root of origNbr.
  common.misc.nbrLimbs = common.misc.toProcess.nbrLimbs;
  squareRoot(common.misc.toProcess.limbs, common.misc.biFirstTerm.limbs, common.misc.nbrLimbs,
    &common.misc.biFirstTerm.nbrLimbs);
  (void)BigIntMultiply(&common.misc.biFirstTerm, &common.misc.biFirstTerm, &common.misc.biMult4);
  BigIntSubt(&common.misc.toProcess, &common.misc.biMult4, &common.misc.biMult4);
  if (BigIntIsZero(&common.misc.biMult4))
  {          // number is a perfect square.
    CopyBigInt(&common.misc.biMult1, &common.misc.biFirstTerm);
    intToBigInteger(&common.misc.biMult2, 0);
    intToBigInteger(&common.misc.biMult3, 0);
  }
  else
  {          // number is not a perfect square.
    common.misc.nbrModExp = 0;
    InitSieveArray();
#ifdef __EMSCRIPTEN__
    // Show number to be processed on screen.
    ptrOutput = common.misc.tmpOutput;
    copyStr(&ptrOutput, "1<p><var>n</var> = ");
    if (hexadecimal)
    {
      Bin2Hex(&ptrOutput, common.misc.biSecondTerm.limbs, common.misc.biSecondTerm.nbrLimbs, common.misc.groupLength);
    }
    else
    {
      Bin2Dec(&ptrOutput, common.misc.biSecondTerm.limbs, common.misc.biSecondTerm.nbrLimbs, common.misc.groupLength);
    }
    copyStr(&ptrOutput, "</p>");
    databack(common.misc.tmpOutput);
#endif
    SumOfSquaresNumber();
    CopyBigInt(&common.misc.biMult1, &common.misc.biFirstTerm);
    CopyBigInt(&common.misc.biMult2, &common.misc.biSecondTerm);
    // Sort squares
    SortBigNbrs(&common.misc.biMult1, &common.misc.biMult2);
    SortBigNbrs(&common.misc.biMult1, &common.misc.biMult3);
    SortBigNbrs(&common.misc.biMult1, &common.misc.biMult4);
    SortBigNbrs(&common.misc.biMult2, &common.misc.biMult3);
    SortBigNbrs(&common.misc.biMult2, &common.misc.biMult4);
    SortBigNbrs(&common.misc.biMult3, &common.misc.biMult4);
  }
  (void)BigIntMultiplyPower2(&common.misc.toProcess, 2* common.misc.power4);
  (void)BigIntMultiplyPower2(&common.misc.biMult1, common.misc.power4);
  (void)BigIntMultiplyPower2(&common.misc.biMult2, common.misc.power4);
  (void)BigIntMultiplyPower2(&common.misc.biMult3, common.misc.power4);
  (void)BigIntMultiplyPower2(&common.misc.biMult4, common.misc.power4);
  // Validate result.
  BigIntMultiply(&common.misc.biMult1, &common.misc.biMult1, &common.misc.biFirstTerm);
  BigIntMultiply(&common.misc.biMult2, &common.misc.biMult2, &common.misc.biSecondTerm);
  BigIntAdd(&common.misc.biFirstTerm, &common.misc.biSecondTerm, &common.misc.biFirstTerm);
  BigIntMultiply(&common.misc.biMult3, &common.misc.biMult3, &common.misc.biSecondTerm);
  BigIntAdd(&common.misc.biFirstTerm, &common.misc.biSecondTerm, &common.misc.biFirstTerm);
  BigIntMultiply(&common.misc.biMult4, &common.misc.biMult4, &common.misc.biSecondTerm);
  BigIntAdd(&common.misc.biFirstTerm, &common.misc.biSecondTerm, &common.misc.biFirstTerm);
  if (BigIntEqual(&common.misc.biFirstTerm, &common.misc.toProcess))
  {
    return 0;    // Validation is OK.
  }
  return 1;      // Validation failed.
}

void fsquaresText(char *input, int grpLen)
{
  char *ptrOutput;
  if (valuesProcessed == 0)
  {
    common.misc.groupLength = grpLen;
  }
  (void)BatchProcessing(input, &common.misc.toProcess, &ptrOutput, NULL, batchSquaresCallback);
  copyStr(&ptrOutput, "<p>");
  showElapsedTime(&ptrOutput);
  copyStr(&ptrOutput, "</p><p>");
  showCopyright(&ptrOutput);
  copyStr(&ptrOutput, "</p>");
}

static void batchSquaresCallback(char **pptrOutput, int type)
{
  int rc;
  char *ptrOutput;
  ptrOutput = *pptrOutput;
  NumberLength = common.misc.toProcess.nbrLimbs;
  BigInteger2IntArray((int *)common.misc.number, &common.misc.toProcess);
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
  rc = fsquares();
  // Show the number to be decomposed into sum of squares.
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
  if (common.misc.toProcess.sign == SIGN_NEGATIVE)
  {
    if (type == BATCH_NO_QUOTE)
    {
      *ptrOutput = ':';
      ptrOutput++;
      *ptrOutput = ' ';
      ptrOutput++;
    }
    textError(&ptrOutput, EXPR_NUMBER_TOO_LOW);
    if (type == BATCH_NO_QUOTE)
    {
      copyStr(&ptrOutput, "</p>");
    }
    *pptrOutput = ptrOutput;
    return;
  }
  switch (rc)
  {
  case 1:
    formatString(&ptrOutput, "<p>$1s</p><p>$2s</p>",
      LITERAL_BATCH_FSQUARES_CBACK1, LITERAL_BATCH_FSQUARES_CBACK2);
    *pptrOutput = ptrOutput;
    return;
  case 2:
    copyStr(&ptrOutput, LITERAL_BATCH_FSQUARES_CBACK3);
    *pptrOutput = ptrOutput;
    return;
  default:
    break;
  }
  // Show the decomposition.
  if (type == BATCH_NO_QUOTE)
  {
    copyStr(&ptrOutput, " = ");
  }
  if (hexadecimal)
  {
    Bin2Hex(&ptrOutput, common.misc.biMult1.limbs, common.misc.biMult1.nbrLimbs, common.misc.groupLength);
  }
  else
  {
    Bin2Dec(&ptrOutput, common.misc.biMult1.limbs, common.misc.biMult1.nbrLimbs, common.misc.groupLength);
  }
  copyStr(&ptrOutput, square);
  if (!BigIntIsZero(&common.misc.biMult2))
  {
    copyStr(&ptrOutput, " + ");
    if (hexadecimal)
    {
      Bin2Hex(&ptrOutput, common.misc.biMult2.limbs, common.misc.biMult2.nbrLimbs, common.misc.groupLength);
    }
    else
    {
      Bin2Dec(&ptrOutput, common.misc.biMult2.limbs, common.misc.biMult2.nbrLimbs, common.misc.groupLength);
    }
    copyStr(&ptrOutput, square);
  }
  if (!BigIntIsZero(&common.misc.biMult3))
  {
    copyStr(&ptrOutput, " + ");
    if (hexadecimal)
    {
      Bin2Hex(&ptrOutput, common.misc.biMult3.limbs, common.misc.biMult3.nbrLimbs, common.misc.groupLength);
    }
    else
    {
      Bin2Dec(&ptrOutput, common.misc.biMult3.limbs, common.misc.biMult3.nbrLimbs, common.misc.groupLength);
    }
    copyStr(&ptrOutput, square);
  }
  if (!BigIntIsZero(&common.misc.biMult4))
  {
    copyStr(&ptrOutput, " + ");
    if (hexadecimal)
    {
      Bin2Hex(&ptrOutput, common.misc.biMult4.limbs, common.misc.biMult4.nbrLimbs, common.misc.groupLength);
    }
    else
    {
      Bin2Dec(&ptrOutput, common.misc.biMult4.limbs, common.misc.biMult4.nbrLimbs, common.misc.groupLength);
    }
    copyStr(&ptrOutput, square);
  }
  if (type == BATCH_NO_QUOTE)
  {
    copyStr(&ptrOutput, "</p>");
  }
  *pptrOutput = ptrOutput;
}

#if defined(__EMSCRIPTEN__) && !defined(_MSC_VER)
#ifdef __ANDROID__
EXTERNALIZE void doWorkFSquares(void)
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
    fsquaresText(NULL, 0); // Routine does not use parameters in this case.
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
  fsquaresText(ptrData + 2, grpLen);
  databack(output);
}
#endif

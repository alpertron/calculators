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
#include "bignbr.h"
#include "expression.h"
#include "highlevel.h"
#include "showtime.h"
#include "batch.h"
#include "tsquares.h"

static int delta[MAX_SIEVE];
static int power4;
static int groupLength;
static void batchSquaresCallback(char** pptrOutput, int type);
#ifdef __EMSCRIPTEN__
  static char tmpOutput[MAX_LEN*12];
#endif
static BigInteger toProcess;
static BigInteger biFirstTerm;
static BigInteger biSecondTerm;
extern limb TestNbr[MAX_LEN];
extern limb MontgomeryMultR1[MAX_LEN];

 // If Mult1 < Mult2, exchange both numbers.
static void SortBigNbrs(BigInteger *pbiMult1, BigInteger *pbiMult2)
{   // biFirstTerm not used at this time. Use it as temp variable.
  BigIntSubt(pbiMult1, pbiMult2, &biFirstTerm);
  if (biFirstTerm.sign == SIGN_NEGATIVE)
  {                // mult1 < mult2 -> exchange them.
    CopyBigInt(&biFirstTerm, pbiMult1);
    CopyBigInt(pbiMult1, pbiMult2);
    CopyBigInt(pbiMult2, &biFirstTerm);
  }
}

static void UpdateSieveArray(void)
{
  for (int index = 0; index < MAX_SIEVE; index++)
  {
    if (sieve[index] >= 0)
    {
      int modulus = (2 * index) + 3;
      sieve[index] = (sieve[index] + delta[index]) % modulus;
      delta[index] -= 2;
      if (delta[index] < 0)
      {
        delta[index] += modulus;
      }
    }
  }
}

// Try to decompose a number as a sum of three or four squares.
// Mult1 = sqrt(toProcess).
// biMult4 = toProcess - biFirstTerm^2.
// Use biMult1, biMult2 and biMult3 as temporary storage.

static void SumOfSquaresNumber(void)
{
  attempts = 0;
  FillSieveArray(toProcess.limbs, sieve);
  CopyBigInt(&biMult3, &toProcess);
  int nbrBytes = toProcess.nbrLimbs * (int)sizeof(int);
  if ((toProcess.limbs[0].x & 7) != 7)
  {              // n!=7 (mod 8) => Sum of three squares
                 // Compute biMult4 as a sum of two squares.
    CopyBigInt(&biSecondTerm, &biFirstTerm);
    intToBigInteger(&biFirstTerm, 0);
    (void)memcpy(valueP, toProcess.limbs, nbrBytes);
    nbrLimbsP = toProcess.nbrLimbs;
    if (isSumOfTwoSquares())
    {            // Number is a sum of two squares: biMult1^2 + biMult2^2.
      intToBigInteger(&biFirstTerm, 0);
      intToBigInteger(&biSecondTerm, 0);
      CopyBigInt(&biMult3, &biMult1);
      CopyBigInt(&biMult4, &biMult2);
      return;
    }
  }
  else
  {              // n==7 (mod 8) => Sum of four squares
                 // Compute biMult4 as a sum of three squares.
    for (;;)
    {
      int pwr4;
      CopyBigInt(&biMult3, &biMult4);
      DivideBigNbrByMaxPowerOf4(&pwr4, biMult3.limbs, &biMult3.nbrLimbs);
      if ((biMult3.limbs[0].x & 7) != 7)
      {          // biMult4 is a sum of three squares. Compute them.
        break;
      }
      addbigint(&biFirstTerm, -1);  // Next candidate for largest term.
      BigIntAdd(&biMult4, &biFirstTerm, &biMult4);
      BigIntAdd(&biMult4, &biFirstTerm, &biMult4);
      addbigint(&biMult4, 1);
    }
    squareRoot(biMult4.limbs, biSecondTerm.limbs, biMult4.nbrLimbs,
      &biSecondTerm.nbrLimbs);
    (void)BigIntMultiply(&biSecondTerm, &biSecondTerm, &biMult3);
    BigIntSubt(&biMult4, &biMult3, &biMult4);
  }
  nbrLimbs = biMult4.nbrLimbs;
  FillSieveArray(biMult4.limbs, sieve);
  multint(&biMult3, &biSecondTerm, 2);
  addbigint(&biMult3, -1);
  nbrLimbs = biMult3.nbrLimbs;
  FillSieveArray(biMult3.limbs, delta);
  // Compute biMult4 as a sum of two squares:
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
    nbrBytes = biMult4.nbrLimbs * (int)sizeof(int);
    (void)memcpy(valueP, biMult4.limbs, nbrBytes);
    nbrLimbsP = biMult4.nbrLimbs;
    if (isSumOfTwoSquares())
    {   // Found biMult1^2 + biMult2^2.
      CopyBigInt(&biMult3, &biMult1);
      CopyBigInt(&biMult4, &biMult2);
      break;
    }
    addbigint(&biSecondTerm, -1);
    BigIntAdd(&biMult4, &biSecondTerm, &biMult4);
    BigIntAdd(&biMult4, &biSecondTerm, &biMult4);
    addbigint(&biMult4, 1);
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
  CopyBigInt(&biSecondTerm, &toProcess);
  DivideBigNbrByMaxPowerOf4(&power4, toProcess.limbs, &toProcess.nbrLimbs);
  // Get biFirstTerm <- square root of origNbr.
  nbrLimbs = toProcess.nbrLimbs;
  squareRoot(toProcess.limbs, biFirstTerm.limbs, nbrLimbs,
    &biFirstTerm.nbrLimbs);
  (void)BigIntMultiply(&biFirstTerm, &biFirstTerm, &biMult4);
  BigIntSubt(&toProcess, &biMult4, &biMult4);
  if (BigIntIsZero(&biMult4))
  {          // number is a perfect square.
    CopyBigInt(&biMult1, &biFirstTerm);
    intToBigInteger(&biMult2, 0);
    intToBigInteger(&biMult3, 0);
  }
  else
  {          // number is not a perfect square.
    nbrModExp = 0;
    InitSieveArray();
#ifdef __EMSCRIPTEN__
    // Show number to be processed on screen.
    ptrOutput = tmpOutput;
    copyStr(&ptrOutput, "1<p><var>n</var> = ");
    if (hexadecimal)
    {
      Bin2Hex(&ptrOutput, biSecondTerm.limbs, biSecondTerm.nbrLimbs, groupLength);
    }
    else
    {
      Bin2Dec(&ptrOutput, biSecondTerm.limbs, biSecondTerm.nbrLimbs, groupLength);
    }
    copyStr(&ptrOutput, "</p>");
    databack(tmpOutput);
#endif
    SumOfSquaresNumber();
    CopyBigInt(&biMult1, &biFirstTerm);
    CopyBigInt(&biMult2, &biSecondTerm);
    // Sort squares
    SortBigNbrs(&biMult1, &biMult2);
    SortBigNbrs(&biMult1, &biMult3);
    SortBigNbrs(&biMult1, &biMult4);
    SortBigNbrs(&biMult2, &biMult3);
    SortBigNbrs(&biMult2, &biMult4);
    SortBigNbrs(&biMult3, &biMult4);
  }
  (void)BigIntMultiplyPower2(&toProcess, 2*power4);
  (void)BigIntMultiplyPower2(&biMult1, power4);
  (void)BigIntMultiplyPower2(&biMult2, power4);
  (void)BigIntMultiplyPower2(&biMult3, power4);
  (void)BigIntMultiplyPower2(&biMult4, power4);
  // Validate result.
  BigIntMultiply(&biMult1, &biMult1, &biFirstTerm);
  BigIntMultiply(&biMult2, &biMult2, &biSecondTerm);
  BigIntAdd(&biFirstTerm, &biSecondTerm, &biFirstTerm);
  BigIntMultiply(&biMult3, &biMult3, &biSecondTerm);
  BigIntAdd(&biFirstTerm, &biSecondTerm, &biFirstTerm);
  BigIntMultiply(&biMult4, &biMult4, &biSecondTerm);
  BigIntAdd(&biFirstTerm, &biSecondTerm, &biFirstTerm);
  if (BigIntEqual(&biFirstTerm, &toProcess))
  {
    return 0;    // Validation is OK.
  }
  return 1;      // Validation failed.
}

void fsquaresText(char *input, int grpLen)
{
  char *ptrOutput;
#ifdef __EMSCRIPTEN__
  int elapsedTime;
#endif
  if (valuesProcessed == 0)
  {
    groupLength = grpLen;
  }
  (void)BatchProcessing(input, &toProcess, &ptrOutput, NULL, batchSquaresCallback);
#ifdef __EMSCRIPTEN__
  copyStr(&ptrOutput, lang ? "<p>Transcurrió " : "<p>Time elapsed: ");
  elapsedTime = (int)(tenths() - originalTenthSecond);
  GetDHMSt(&ptrOutput, elapsedTime);
#endif
  copyStr(&ptrOutput, "<p>");
  copyStr(&ptrOutput, (lang ? COPYRIGHT_SPANISH: COPYRIGHT_ENGLISH));
  copyStr(&ptrOutput, "</p>");
}

static void batchSquaresCallback(char **pptrOutput, int type)
{
  int rc;
  char *ptrOutput;
  ptrOutput = *pptrOutput;
  NumberLength = toProcess.nbrLimbs;
  BigInteger2IntArray((int *)number, &toProcess);
  if ((type == BATCH_NO_PROCESS_DEC) || (type == BATCH_NO_PROCESS_HEX))
  {         // Do not compute sum of squares.
    if (hexadecimal)
    {
      BigInteger2Hex(&ptrOutput, &toProcess, groupLength);
    }
    else
    {
      BigInteger2Dec(&ptrOutput, &toProcess, groupLength);
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
      BigInteger2Hex(&ptrOutput, &toProcess, groupLength);
    }
    else
    {
      BigInteger2Dec(&ptrOutput, &toProcess, groupLength);
    }
  }
  if (toProcess.sign == SIGN_NEGATIVE)
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
    copyStr(&ptrOutput, (lang ? ": ¡Error interno!\n\nPor favor envíe este número al autor del applet.</p>":
      ": Internal error!\n\nPlease send the number to the author of the applet.</p>"));
    *pptrOutput = ptrOutput;
    return;
  case 2:
    copyStr(&ptrOutput, (lang?": El usuario detuvo el cálculo": ": User stopped the calculation"));
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
    Bin2Hex(&ptrOutput, biMult1.limbs, biMult1.nbrLimbs, groupLength);
  }
  else
  {
    Bin2Dec(&ptrOutput, biMult1.limbs, biMult1.nbrLimbs, groupLength);
  }
  copyStr(&ptrOutput, square);
  if (!BigIntIsZero(&biMult2))
  {
    copyStr(&ptrOutput, " + ");
    if (hexadecimal)
    {
      Bin2Hex(&ptrOutput, biMult2.limbs, biMult2.nbrLimbs, groupLength);
    }
    else
    {
      Bin2Dec(&ptrOutput, biMult2.limbs, biMult2.nbrLimbs, groupLength);
    }
    copyStr(&ptrOutput, square);
  }
  if (!BigIntIsZero(&biMult3))
  {
    copyStr(&ptrOutput, " + ");
    if (hexadecimal)
    {
      Bin2Hex(&ptrOutput, biMult3.limbs, biMult3.nbrLimbs, groupLength);
    }
    else
    {
      Bin2Dec(&ptrOutput, biMult3.limbs, biMult3.nbrLimbs, groupLength);
    }
    copyStr(&ptrOutput, square);
  }
  if (!BigIntIsZero(&biMult4))
  {
    copyStr(&ptrOutput, " + ");
    if (hexadecimal)
    {
      Bin2Hex(&ptrOutput, biMult4.limbs, biMult4.nbrLimbs, groupLength);
    }
    else
    {
      Bin2Dec(&ptrOutput, biMult4.limbs, biMult4.nbrLimbs, groupLength);
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
EXTERNALIZE void doWork(void)
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
  fsquaresText(ptrData + 2, grpLen);
  databack(output);
}
#endif

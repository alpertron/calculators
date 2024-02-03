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
void DivideBigNbrByMaxPowerOf2(int *pShRight, limb *number, int *pNbrLimbs);

 // If Mult1 < Mult2, exchange both numbers.
static void SortBigNbrs(limb *mult1, int *mult1Len, limb *mult2, int *mult2Len)
{
  limb tmp;
  int index;
  limb *ptr1;
  limb *ptr2;
  int len = *mult1Len;
  if (len > *mult2Len)
  {
    return;    // mult1 > mult2, so nothing to do.
  }
  if (len == *mult2Len)
  {
    ptr1 = &mult1[len-1];
    ptr2 = &mult2[len-1];
    for (index = *mult1Len - 1; index >= 0; index--)
    {
      if (ptr1->x > ptr2->x)
      {
        return;    // mult1 > mult2, so nothing to do.
      }
      if (ptr1->x < ptr2->x)
      {
        break;     // mult1 < mult2, so exchange them.
      }
      ptr1--;
      ptr2--;
    }
  }
    // Exchange lengths.
  len = *mult2Len;
  *mult2Len = *mult1Len;
  *mult1Len = len;
    // Exchange bytes that compose the numbers.
  ptr1 = mult1;
  ptr2 = mult2;
  for (index = 0; index < len; index++)
  {
    tmp.x = ptr1->x;
    ptr1->x = ptr2->x;
    ptr2->x = tmp.x;
    ptr1++;
    ptr2++;
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
  if (nbrLimbs == 1)
  {
    if (number[0].x == 3)
    {      // Number is 3 = 1^2 + 1^2 + 1^2.
           // At this moment Mult1 = 1.
      Mult2[0].x = 1;
      Mult2Len = 1;
      Mult3[0].x = 1;
      Mult3Len = 1;
      Mult4[0].x = 0;
      Mult4Len = 1;
      return;
    }
    if (number[0].x == 2)
    {      // Number is 2 = 1^2 + 1^2.
           // At this moment Mult1 = 1.
      Mult2[0].x = 1;
      Mult2Len = 1;
      Mult3[0].x = 0;
      Mult3Len = 1;
      Mult4[0].x = 0;
      Mult4Len = 1;
      return;
    }
  }
  attempts = 0;
  FillSieveArray(toProcess.limbs, sieve);
  CopyBigInt(&biMult3, &toProcess);
  DivideBigNbrByMaxPowerOf4(&power4, biMult3.limbs, &biMult3.nbrLimbs);
  int nbrBytes = toProcess.nbrLimbs * (int)sizeof(int);
  if ((biMult3.limbs[0].x & 7) != 7)
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
      Mult3Len = biMult1.nbrLimbs;
      Mult4Len = biMult2.nbrLimbs;
      int mult1LenBytes = Mult3Len * (int)sizeof(int);
      int mult2LenBytes = Mult4Len * (int)sizeof(int);
      (void)memcpy(Mult3, biMult1.limbs, mult1LenBytes);
      (void)memcpy(Mult4, biMult2.limbs, mult2LenBytes);
      return;
    }
  }
  else
  {              // n==7 (mod 8) => Sum of four squares
                 // Compute biMult4 as a sum of three squares.
    for (;;)
    {
      CopyBigInt(&biMult3, &biMult4);
      DivideBigNbrByMaxPowerOf4(&power4, biMult3.limbs, &biMult3.nbrLimbs);
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
      Mult3Len = biMult1.nbrLimbs;
      Mult4Len = biMult2.nbrLimbs;
      int mult1LenBytes = Mult3Len * (int)sizeof(int);
      int mult2LenBytes = Mult4Len * (int)sizeof(int);
      (void)memcpy(Mult3, biMult1.limbs, mult1LenBytes);
      (void)memcpy(Mult4, biMult2.limbs, mult2LenBytes);
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
  int tmp;
  int idx;
  int lenBytes;
  nbrLimbs = toProcess.nbrLimbs;
  // Get Mult1 <- square root of origNbr.
  squareRoot(toProcess.limbs, biFirstTerm.limbs, nbrLimbs,
    &biFirstTerm.nbrLimbs);
  (void)BigIntMultiply(&biFirstTerm, &biFirstTerm, &biMult4);
  BigIntSubt(&toProcess, &biMult4, &biMult4);
  if (BigIntIsZero(&biMult4))
  {          // number is a perfect square.
    Mult1Len = biFirstTerm.nbrLimbs;
    int Mult1LenBytes = Mult1Len * (int)sizeof(int);
    (void)memcpy(Mult1, biFirstTerm.limbs, Mult1LenBytes);
    Mult2[0].x = 0;
    Mult2Len = 1;
    Mult3[0].x = 0;
    Mult3Len = 1;
    Mult4[0].x = 0;
    Mult4Len = 1;
    power4 = 0;
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
      Bin2Hex(&ptrOutput, toProcess.limbs, toProcess.nbrLimbs, groupLength);
    }
    else
    {
      Bin2Dec(&ptrOutput, toProcess.limbs, toProcess.nbrLimbs, groupLength);
    }
    copyStr(&ptrOutput, "</p>");
    databack(tmpOutput);
#endif
    SumOfSquaresNumber();
    Mult1Len = biFirstTerm.nbrLimbs;
    int Mult1LenBytes = Mult1Len * (int)sizeof(int);
    (void)memcpy(Mult1, biFirstTerm.limbs, Mult1LenBytes);
    Mult2Len = biSecondTerm.nbrLimbs;
    int Mult2LenBytes = Mult2Len * (int)sizeof(int);
    (void)memcpy(Mult2, biSecondTerm.limbs, Mult2LenBytes);
    while ((Mult3[Mult3Len-1].x == 0) && (Mult3Len > 1))
    {
      Mult3Len--;
    }
    while ((Mult4[Mult4Len-1].x == 0) && (Mult4Len > 1))
    {
      Mult4Len--;
    }
    // Sort squares
    SortBigNbrs(Mult1, &Mult1Len, Mult2, &Mult2Len);
    SortBigNbrs(Mult1, &Mult1Len, Mult3, &Mult3Len);
    SortBigNbrs(Mult1, &Mult1Len, Mult4, &Mult4Len);
    SortBigNbrs(Mult2, &Mult2Len, Mult3, &Mult3Len);
    SortBigNbrs(Mult2, &Mult2Len, Mult4, &Mult4Len);
    SortBigNbrs(Mult3, &Mult3Len, Mult4, &Mult4Len);
  }
    // Validate result.
  idx = Mult1Len * 2;
  multiply(Mult1, Mult1, SquareMult1.limbs, Mult1Len, &tmp);
  SquareMult1.limbs[idx].x = 0;
  multiply(Mult2, Mult2, SquareMult2.limbs, Mult2Len, &tmp);
  lenBytes = ((Mult1Len - Mult2Len) * 2) * (int)sizeof(limb);
  (void)memset(&SquareMult2.limbs[Mult2Len * 2], 0, lenBytes);
  SquareMult2.limbs[idx].x = 0;
  multiply(Mult3, Mult3, SquareMult3.limbs, Mult3Len, &tmp);
  lenBytes = ((Mult1Len - Mult3Len) * 2) * (int)sizeof(limb);
  (void)memset(&SquareMult3.limbs[Mult3Len * 2], 0, lenBytes);
  SquareMult3.limbs[idx].x = 0;
  multiply(Mult4, Mult4, SquareMult4.limbs, Mult4Len, &tmp);
  lenBytes = ((Mult1Len - Mult4Len) * 2) * (int)sizeof(limb);
  (void)memset(&SquareMult4.limbs[Mult4Len * 2], 0, lenBytes);
  SquareMult4.limbs[idx].x = 0;
  idx++;
  AddBigInt(SquareMult1.limbs, SquareMult2.limbs, SquareMult1.limbs, idx);
  AddBigInt(SquareMult1.limbs, SquareMult3.limbs, SquareMult1.limbs, idx);
  AddBigInt(SquareMult1.limbs, SquareMult4.limbs, SquareMult1.limbs, idx);
  while ((idx > 1) && (SquareMult1.limbs[idx - 1].x == 0))
  {
    idx--;
  }
  if (idx != toProcess.nbrLimbs)
  {           // Invalid length.
    return 1;
  }
  for (int index = 0; index < idx; index++)
  {
    if (SquareMult1.limbs[index].x != toProcess.limbs[index].x)
    {
      return 1;
    }
  }
  return 0;
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
    Bin2Hex(&ptrOutput, Mult1, Mult1Len, groupLength);
  }
  else
  {
    Bin2Dec(&ptrOutput, Mult1, Mult1Len, groupLength);
  }
  copyStr(&ptrOutput, square);
  if ((Mult2Len != 1) || (Mult2[0].x != 0))
  {
    copyStr(&ptrOutput, " + ");
    if (hexadecimal)
    {
      Bin2Hex(&ptrOutput, Mult2, Mult2Len, groupLength);
    }
    else
    {
      Bin2Dec(&ptrOutput, Mult2, Mult2Len, groupLength);
    }
    copyStr(&ptrOutput, square);
  }
  if ((Mult3Len != 1) || (Mult3[0].x != 0))
  {
    copyStr(&ptrOutput, " + ");
    if (hexadecimal)
    {
      Bin2Hex(&ptrOutput, Mult3, Mult3Len, groupLength);
    }
    else
    {
      Bin2Dec(&ptrOutput, Mult3, Mult3Len, groupLength);
    }
    copyStr(&ptrOutput, square);
  }
  if ((Mult4Len != 1) || (Mult4[0].x != 0))
  {
    copyStr(&ptrOutput, " + ");
    if (hexadecimal)
    {
      Bin2Hex(&ptrOutput, Mult4, Mult4Len, groupLength);
    }
    else
    {
      Bin2Dec(&ptrOutput, Mult4, Mult4Len, groupLength);
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

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

static limb origNbr[MAX_LEN];
static int power4;
static limb result[MAX_LEN];
static int origNbrLimbs;
static int groupLength;
static void batchSquaresCallback(char** pptrOutput, int type);
#ifdef __EMSCRIPTEN__
  static char tmpOutput[MAX_LEN*12];
#endif
static BigInteger toProcess;
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

// Try to decompose a number as a sum of two squares.
// p = Mult1^2 + Mult2^2
// Use biMult1, biMult2 and biMult3 as temporary storage.

static void SumOfSquaresNumberGreaterThan3(void)
{
  limb carry;
  FillSieveArray(number);
  if ((number[0].x & 7) != 7)
  {              // n!=7 (mod 8) => Sum of three squares
    Computing3Squares = true;
    iMult4 = 0;
    iMult3 = -1;
  }
  else
  {              // n==7 (mod 8) => Sum of four squares
    Computing3Squares = false;
    iMult3 = -1;
    iMult4 = 0;
  }
  // If number is a sum of three squares, subtract a small square.
  // If number is a sum of four squares, subtract two small squares.
  // If the result is not the product of a power of 2, small primes
  // of the form 4k+1 and squares of primes of the form (4k+3)^2
  // and a (big) prime, try with other squares.
  do
  {
    if (!Computing3Squares && (iMult3 >= iMult4))
    {
      iMult3 = 1;
      iMult4++;
    }
    else
    {
      iMult3++;
    }
    sum = (iMult3 * iMult3) + (iMult4 * iMult4);
    carry.x = number[0].x - sum;
    valueP[0].x = UintToInt((unsigned int)carry.x & MAX_VALUE_LIMB);
    carry.x >>= BITS_PER_GROUP;
    if (nbrLimbs > 1)
    {
      carry.x += number[1].x;
      valueP[1].x = UintToInt((unsigned int)carry.x & MAX_VALUE_LIMB);
      carry.x >>= BITS_PER_GROUP;
      for (int index = 2; index < nbrLimbs; index++)
      {
        carry.x += number[index].x;
        valueP[index].x = UintToInt((unsigned int)carry.x & MAX_VALUE_LIMB);
        carry.x >>= BITS_PER_GROUP;
      }
    }
    // p should be the product of power of 2,
    // powers of small primes of form 4k+1 and
    // powers of squares of small primes of form 4k+3.
    nbrLimbsP = nbrLimbs;
  } while (!isSumOfTwoSquares());        /* end while */
}

// Variable to split in up to four squares: number
int fsquares(void)
{
#ifdef __EMSCRIPTEN__
  char *ptrOutput;
#endif
  int tmp;
  int count;
  int idx;
  int lenBytes;
  nbrLimbs = origNbrLimbs;
  lenBytes = nbrLimbs * (int)sizeof(limb);
  (void)memcpy(number, origNbr, lenBytes);
  squareRoot(number, Mult1, nbrLimbs, &Mult1Len);
  multiply(Mult1, Mult1, result, Mult1Len, NULL);
  if (memcmp(result, number, lenBytes) == 0)
  {          // number is a perfect square.
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
    ptrOutput = tmpOutput;
    copyStr(&ptrOutput, "1<p><var>n</var> = ");
    if (hexadecimal)
    {
      Bin2Hex(&ptrOutput, origNbr, nbrLimbs, groupLength);
    }
    else
    {
      Bin2Dec(&ptrOutput, origNbr, nbrLimbs, groupLength);
    }
    copyStr(&ptrOutput, "</p>");
    databack(tmpOutput);
#endif
    DivideBigNbrByMaxPowerOf4(&power4, number, &nbrLimbs);
    Mult1Len = 1;
    Mult2Len = 1;
    if ((nbrLimbs == 1) && (number[0].x < 4))
    {
      iMult3 = 0;
      iMult4 = 0;
      switch (number[0].x)
      {
      case 3:
        iMult3 = 1;
        Mult2[0].x = 1;
        Mult1[0].x = 1;
        break;
      case 2:
        Mult2[0].x = 1;
        Mult1[0].x = 1;
        break;
      case 1:
        Mult1[0].x = 1;
        break;
      default:
        break;
      }
    }
    else
    {     // Number greater than 3.
      SumOfSquaresNumberGreaterThan3();
    }
    // Shift left the number of bits that the original number was divided by 4.
    for (count = 0; count < power4; count++)
    {
      MultBigNbrByInt(Mult1, 2, Mult1, nbrLimbs + 1);
      MultBigNbrByInt(Mult2, 2, Mult2, nbrLimbs + 1);
      if ((Mult1[nbrLimbs].x != 0) || (Mult2[nbrLimbs].x != 0))
      {
        nbrLimbs++;
        Mult1[nbrLimbs].x = 0;
        Mult2[nbrLimbs].x = 0;
      }
    }
    nbrLimbsP = 1;
    Mult3[0].x = iMult3;
    Mult4[0].x = iMult4;
    Mult3[1].x = 0;
    Mult4[1].x = 0;
    for (count = 0; count < power4; count++)
    {
      MultBigNbrByInt(Mult3, 2, Mult3, nbrLimbsP + 1);
      MultBigNbrByInt(Mult4, 2, Mult4, nbrLimbsP + 1);
      if ((Mult3[nbrLimbsP].x != 0) || (Mult4[nbrLimbsP].x != 0))
      {
        nbrLimbsP++;
      }
      Mult3[nbrLimbsP].x = 0;
      Mult4[nbrLimbsP].x = 0;
    }
    Mult1Len = nbrLimbs;
    Mult2Len = nbrLimbs;
    while ((Mult1[Mult1Len-1].x == 0) && (Mult1Len > 1))
    {
      Mult1Len--;
    }
    while ((Mult2[Mult2Len-1].x == 0) && (Mult2Len > 1))
    {
      Mult2Len--;
    }
    Mult3Len = nbrLimbsP;
    Mult4Len = nbrLimbsP;
    while ((Mult3[Mult3Len-1].x == 0) && (Mult3Len > 1))
    {
      Mult3Len--;
    }
    while ((Mult1[Mult4Len-1].x == 0) && (Mult4Len > 1))
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
  if (idx != origNbrLimbs)
  {           // Invalid length.
    return 1;
  }
  for (int index = 0; index < idx; index++)
  {
    if (SquareMult1.limbs[index].x != origNbr[index].x)
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
  copyStr(&ptrOutput, "</p><p>");
  copyStr(&ptrOutput, (lang ? COPYRIGHT_SPANISH: COPYRIGHT_ENGLISH));
  copyStr(&ptrOutput, "</p>");
}

static void batchSquaresCallback(char **pptrOutput, int type)
{
  int rc;
  char *ptrOutput;
  int lenBytes;
  ptrOutput = *pptrOutput;
  NumberLength = toProcess.nbrLimbs;
  BigInteger2IntArray((int *)number, &toProcess);
  origNbrLimbs = toProcess.nbrLimbs;
  lenBytes = origNbrLimbs * (int)sizeof(limb);
  (void)memcpy(origNbr, toProcess.limbs, lenBytes);
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

#if defined __EMSCRIPTEN__ && !defined _MSC_VER
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

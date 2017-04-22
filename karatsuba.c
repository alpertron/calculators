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

/* Perform Karatsuba multiplication in little endian order */

#include <string.h>
#include "bignbr.h"
#include <stdio.h>
#include <math.h>
#include <stdint.h>

#define KARATSUBA_CUTOFF 8
static limb arr[MAX_LEN];
static limb arrayAux[MAX_LEN];
static int karatLength;
static void Karatsuba(int idxFactor1, int length, int diffIndex);

void multiply(limb *factor1, limb *factor2, limb *result, int len, int *pResultLen)
{
  int length = len;
    // Compute length of numbers for each recursion.
  if (length > KARATSUBA_CUTOFF)
  {
    int div = 1;
    while (length > KARATSUBA_CUTOFF)
    {
      div *= 2;
      length = (length + 1) / 2;
    }
    length *= div;
  }
  karatLength = length;
  memset(arr, 0, 2 * length*sizeof(limb));
  memcpy(&arr[0], factor1, len*sizeof(limb));
  memcpy(&arr[length], factor2, len*sizeof(limb));
  Karatsuba(0, length, 2 * length);
  memcpy(result, &arr[2 * (karatLength - length)], 2 * length * sizeof(limb));
  if (pResultLen != NULL)
  {
    memcpy(result, &arr[2 * (karatLength - length)], 2 * length * sizeof(limb));
    if (karatLength > length && arr[2 * (karatLength - length)-1].x == 0)
    {
      *pResultLen = length * 2 - 1;
    }
    else
    {
      *pResultLen = length * 2;
    }
  }
}

// The return value is the sign: true: negative.
// In result the absolute value of the difference is computed.
static int absSubtract(int idxMinuend, int idxSubtrahend,
                       int idxResult, int nbrLen)
{
  int sign = 0;
  limb carry;
  int i;
  for (i = nbrLen-1; i>=0; i--)
  {
    if (arr[idxMinuend + i].x != arr[idxSubtrahend + i].x)
    {
      break;
    }
  }
  if (i>=0 && arr[idxMinuend + i].x < arr[idxSubtrahend + i].x)
  {
    sign = 1;
    i = idxMinuend;    // Exchange minuend and subtrahend.
    idxMinuend = idxSubtrahend;
    idxSubtrahend = i;
  }
  carry.x = 0;
  for (i = 0; i < nbrLen; i++)
  {
    carry.x += arr[idxMinuend + i].x - arr[idxSubtrahend + i].x;
    arr[idxResult + i].x = carry.x & MAX_VALUE_LIMB;
    carry.x >>= BITS_PER_GROUP;
  }
  return sign;
}

// Multiply two groups of nbrLen limbs. The first one starts at idxFactor1
// and the second one at idxFactor2. The 2*nbrLen limb result is stored
// starting at idxFactor1. Use arrayAux as temporary storage.
// Accumulate products by result limb.
#ifdef _USING64BITS_
static void ClassicalMult2Limbs(int idxFactor1, int idxFactor2)
{
  int i;
  int32_t factor2_i;
  int32_t factor1_0 = arr[idxFactor1].x;
  int32_t factor1_1 = arr[idxFactor1+1].x;
  int32_t prod_iPlus0 = 0;
  int32_t prod_iPlus1 = 0;
  int64_t Pr;
  for (i = 0; i < 2; i++)
  {
    factor2_i = arr[idxFactor2 + i].x;
    Pr = prod_iPlus0 + (int64_t)factor2_i * factor1_0;
    arrayAux[i].x = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus1 + (int64_t)factor2_i * factor1_1 + (Pr >> BITS_PER_GROUP);
    prod_iPlus0 = (int32_t)Pr & MAX_INT_NBR;
    prod_iPlus1 = (int32_t)(Pr >> BITS_PER_GROUP);
  }
  arrayAux[2].x = prod_iPlus0;
  arrayAux[3].x = prod_iPlus1;
}

static void ClassicalMult3Limbs(int idxFactor1, int idxFactor2)
{
  int i;
  int32_t factor2_i;
  int32_t factor1_0 = arr[idxFactor1].x;
  int32_t factor1_1 = arr[idxFactor1 + 1].x;
  int32_t factor1_2 = arr[idxFactor1 + 2].x;
  int32_t prod_iPlus0 = 0;
  int32_t prod_iPlus1 = 0;
  int32_t prod_iPlus2 = 0;
  int64_t Pr;
  for (i = 0; i < 3; i++)
  {
    factor2_i = arr[idxFactor2 + i].x;
    Pr = prod_iPlus0 + (int64_t)factor2_i * factor1_0;
    arrayAux[i].x = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus1 + (int64_t)factor2_i * factor1_1 + (Pr >> BITS_PER_GROUP);
    prod_iPlus0 = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus2 + (int64_t)factor2_i * factor1_2 + (Pr >> BITS_PER_GROUP);
    prod_iPlus1 = (int32_t)Pr & MAX_INT_NBR;
    prod_iPlus2 = (int32_t)(Pr >> BITS_PER_GROUP);
  }
  arrayAux[3].x = prod_iPlus0;
  arrayAux[4].x = prod_iPlus1;
  arrayAux[5].x = prod_iPlus2;
}

static void ClassicalMult4Limbs(int idxFactor1, int idxFactor2)
{
  int i;
  int32_t factor2_i;
  int32_t factor1_0 = arr[idxFactor1].x;
  int32_t factor1_1 = arr[idxFactor1 + 1].x;
  int32_t factor1_2 = arr[idxFactor1 + 2].x;
  int32_t factor1_3 = arr[idxFactor1 + 3].x;
  int32_t prod_iPlus0 = 0;
  int32_t prod_iPlus1 = 0;
  int32_t prod_iPlus2 = 0;
  int32_t prod_iPlus3 = 0;
  int64_t Pr;
  for (i = 0; i < 4; i++)
  {
    factor2_i = arr[idxFactor2 + i].x;
    Pr = prod_iPlus0 + (int64_t)factor2_i * factor1_0;
    arrayAux[i].x = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus1 + (int64_t)factor2_i * factor1_1 + (Pr >> BITS_PER_GROUP);
    prod_iPlus0 = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus2 + (int64_t)factor2_i * factor1_2 + (Pr >> BITS_PER_GROUP);
    prod_iPlus1 = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus3 + (int64_t)factor2_i * factor1_3 + (Pr >> BITS_PER_GROUP);
    prod_iPlus2 = (int32_t)Pr & MAX_INT_NBR;
    prod_iPlus3 = (int32_t)(Pr >> BITS_PER_GROUP);
  }
  arrayAux[4].x = prod_iPlus0;
  arrayAux[5].x = prod_iPlus1;
  arrayAux[6].x = prod_iPlus2;
  arrayAux[7].x = prod_iPlus3;
}

static void ClassicalMult5Limbs(int idxFactor1, int idxFactor2)
{
  int i;
  int32_t factor2_i;
  int32_t factor1_0 = arr[idxFactor1].x;
  int32_t factor1_1 = arr[idxFactor1 + 1].x;
  int32_t factor1_2 = arr[idxFactor1 + 2].x;
  int32_t factor1_3 = arr[idxFactor1 + 3].x;
  int32_t factor1_4 = arr[idxFactor1 + 4].x;
  int32_t prod_iPlus0 = 0;
  int32_t prod_iPlus1 = 0;
  int32_t prod_iPlus2 = 0;
  int32_t prod_iPlus3 = 0;
  int32_t prod_iPlus4 = 0;
  int64_t Pr;
  for (i = 0; i < 5; i++)
  {
    factor2_i = arr[idxFactor2 + i].x;
    Pr = prod_iPlus0 + (int64_t)factor2_i * factor1_0;
    arrayAux[i].x = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus1 + (int64_t)factor2_i * factor1_1 + (Pr >> BITS_PER_GROUP);
    prod_iPlus0 = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus2 + (int64_t)factor2_i * factor1_2 + (Pr >> BITS_PER_GROUP);
    prod_iPlus1 = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus3 + (int64_t)factor2_i * factor1_3 + (Pr >> BITS_PER_GROUP);
    prod_iPlus2 = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus4 + (int64_t)factor2_i * factor1_4 + (Pr >> BITS_PER_GROUP);
    prod_iPlus3 = (int32_t)Pr & MAX_INT_NBR;
    prod_iPlus4 = (int32_t)(Pr >> BITS_PER_GROUP);
  }
  arrayAux[5].x = prod_iPlus0;
  arrayAux[6].x = prod_iPlus1;
  arrayAux[7].x = prod_iPlus2;
  arrayAux[8].x = prod_iPlus3;
  arrayAux[9].x = prod_iPlus4;
}

static void ClassicalMult6Limbs(int idxFactor1, int idxFactor2)
{
  int i;
  int32_t factor2_i;
  int32_t factor1_0 = arr[idxFactor1].x;
  int32_t factor1_1 = arr[idxFactor1 + 1].x;
  int32_t factor1_2 = arr[idxFactor1 + 2].x;
  int32_t factor1_3 = arr[idxFactor1 + 3].x;
  int32_t factor1_4 = arr[idxFactor1 + 4].x;
  int32_t factor1_5 = arr[idxFactor1 + 5].x;
  int32_t prod_iPlus0 = 0;
  int32_t prod_iPlus1 = 0;
  int32_t prod_iPlus2 = 0;
  int32_t prod_iPlus3 = 0;
  int32_t prod_iPlus4 = 0;
  int32_t prod_iPlus5 = 0;
  int64_t Pr;
  for (i = 0; i < 6; i++)
  {
    factor2_i = arr[idxFactor2 + i].x;
    Pr = prod_iPlus0 + (int64_t)factor2_i * factor1_0;
    arrayAux[i].x = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus1 + (int64_t)factor2_i * factor1_1 + (Pr >> BITS_PER_GROUP);
    prod_iPlus0 = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus2 + (int64_t)factor2_i * factor1_2 + (Pr >> BITS_PER_GROUP);
    prod_iPlus1 = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus3 + (int64_t)factor2_i * factor1_3 + (Pr >> BITS_PER_GROUP);
    prod_iPlus2 = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus4 + (int64_t)factor2_i * factor1_4 + (Pr >> BITS_PER_GROUP);
    prod_iPlus3 = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus5 + (int64_t)factor2_i * factor1_5 + (Pr >> BITS_PER_GROUP);
    prod_iPlus4 = (int32_t)Pr & MAX_INT_NBR;
    prod_iPlus5 = (int32_t)(Pr >> BITS_PER_GROUP);
  }
  arrayAux[6].x = prod_iPlus0;
  arrayAux[7].x = prod_iPlus1;
  arrayAux[8].x = prod_iPlus2;
  arrayAux[9].x = prod_iPlus3;
  arrayAux[10].x = prod_iPlus4;
  arrayAux[11].x = prod_iPlus5;
}

static void ClassicalMult7Limbs(int idxFactor1, int idxFactor2)
{
  int i;
  int32_t factor2_i;
  int32_t factor1_0 = arr[idxFactor1].x;
  int32_t factor1_1 = arr[idxFactor1 + 1].x;
  int32_t factor1_2 = arr[idxFactor1 + 2].x;
  int32_t factor1_3 = arr[idxFactor1 + 3].x;
  int32_t factor1_4 = arr[idxFactor1 + 4].x;
  int32_t factor1_5 = arr[idxFactor1 + 5].x;
  int32_t factor1_6 = arr[idxFactor1 + 6].x;
  int32_t prod_iPlus0 = 0;
  int32_t prod_iPlus1 = 0;
  int32_t prod_iPlus2 = 0;
  int32_t prod_iPlus3 = 0;
  int32_t prod_iPlus4 = 0;
  int32_t prod_iPlus5 = 0;
  int32_t prod_iPlus6 = 0;
  int64_t Pr;
  for (i = 0; i < 7; i++)
  {
    factor2_i = arr[idxFactor2 + i].x;
    Pr = prod_iPlus0 + (int64_t)factor2_i * factor1_0;
    arrayAux[i].x = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus1 + (int64_t)factor2_i * factor1_1 + (Pr >> BITS_PER_GROUP);
    prod_iPlus0 = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus2 + (int64_t)factor2_i * factor1_2 + (Pr >> BITS_PER_GROUP);
    prod_iPlus1 = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus3 + (int64_t)factor2_i * factor1_3 + (Pr >> BITS_PER_GROUP);
    prod_iPlus2 = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus4 + (int64_t)factor2_i * factor1_4 + (Pr >> BITS_PER_GROUP);
    prod_iPlus3 = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus5 + (int64_t)factor2_i * factor1_5 + (Pr >> BITS_PER_GROUP);
    prod_iPlus4 = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus6 + (int64_t)factor2_i * factor1_6 + (Pr >> BITS_PER_GROUP);
    prod_iPlus5 = (int32_t)Pr & MAX_INT_NBR;
    prod_iPlus6 = (int32_t)(Pr >> BITS_PER_GROUP);
  }
  arrayAux[7].x = prod_iPlus0;
  arrayAux[8].x = prod_iPlus1;
  arrayAux[9].x = prod_iPlus2;
  arrayAux[10].x = prod_iPlus3;
  arrayAux[11].x = prod_iPlus4;
  arrayAux[12].x = prod_iPlus5;
  arrayAux[13].x = prod_iPlus6;
}

static void ClassicalMult8Limbs(int idxFactor1, int idxFactor2)
{
  int i;
  int32_t factor2_i;
  int32_t factor1_0 = arr[idxFactor1].x;
  int32_t factor1_1 = arr[idxFactor1 + 1].x;
  int32_t factor1_2 = arr[idxFactor1 + 2].x;
  int32_t factor1_3 = arr[idxFactor1 + 3].x;
  int32_t factor1_4 = arr[idxFactor1 + 4].x;
  int32_t factor1_5 = arr[idxFactor1 + 5].x;
  int32_t factor1_6 = arr[idxFactor1 + 6].x;
  int32_t factor1_7 = arr[idxFactor1 + 7].x;
  int32_t prod_iPlus0 = 0;
  int32_t prod_iPlus1 = 0;
  int32_t prod_iPlus2 = 0;
  int32_t prod_iPlus3 = 0;
  int32_t prod_iPlus4 = 0;
  int32_t prod_iPlus5 = 0;
  int32_t prod_iPlus6 = 0;
  int32_t prod_iPlus7 = 0;
  int64_t Pr;
  for (i = 0; i < 8; i++)
  {
    factor2_i = arr[idxFactor2 + i].x;
    Pr = prod_iPlus0 + (int64_t)factor2_i * factor1_0;
    arrayAux[i].x = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus1 + (int64_t)factor2_i * factor1_1 + (Pr >> BITS_PER_GROUP);
    prod_iPlus0 = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus2 + (int64_t)factor2_i * factor1_2 + (Pr >> BITS_PER_GROUP);
    prod_iPlus1 = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus3 + (int64_t)factor2_i * factor1_3 + (Pr >> BITS_PER_GROUP);
    prod_iPlus2 = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus4 + (int64_t)factor2_i * factor1_4 + (Pr >> BITS_PER_GROUP);
    prod_iPlus3 = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus5 + (int64_t)factor2_i * factor1_5 + (Pr >> BITS_PER_GROUP);
    prod_iPlus4 = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus6 + (int64_t)factor2_i * factor1_6 + (Pr >> BITS_PER_GROUP);
    prod_iPlus5 = (int32_t)Pr & MAX_INT_NBR;
    Pr = prod_iPlus7 + (int64_t)factor2_i * factor1_7 + (Pr >> BITS_PER_GROUP);
    prod_iPlus6 = (int32_t)Pr & MAX_INT_NBR;
    prod_iPlus7 = (int32_t)(Pr >> BITS_PER_GROUP);
  }
  arrayAux[8].x = prod_iPlus0;
  arrayAux[9].x = prod_iPlus1;
  arrayAux[10].x = prod_iPlus2;
  arrayAux[11].x = prod_iPlus3;
  arrayAux[12].x = prod_iPlus4;
  arrayAux[13].x = prod_iPlus5;
  arrayAux[14].x = prod_iPlus6;
  arrayAux[15].x = prod_iPlus7;
}

#endif
static void ClassicalMult(int idxFactor1, int idxFactor2, int nbrLen)
{
#ifdef _USING64BITS_
  uint64_t product;
  switch (nbrLen)
  {
  case 1:
    product = (int64_t)arr[idxFactor1].x * arr[idxFactor2].x;
    // Store least significant limb.
    arr[idxFactor1].x = (int32_t)product & MAX_INT_NBR;
    // Store most significant limb.
    arr[idxFactor2].x = (int32_t)(product >> BITS_PER_GROUP) & MAX_INT_NBR;
    return;
  case 2:
    ClassicalMult2Limbs(idxFactor1, idxFactor2);
    break;
  case 3:
    ClassicalMult3Limbs(idxFactor1, idxFactor2);
    break;
  case 4:
    ClassicalMult4Limbs(idxFactor1, idxFactor2);
    break;
  case 5:
    ClassicalMult5Limbs(idxFactor1, idxFactor2);
    break;
  case 6:
    ClassicalMult6Limbs(idxFactor1, idxFactor2);
    break;
  case 7:
    ClassicalMult7Limbs(idxFactor1, idxFactor2);
    break;
  case 8:
    ClassicalMult8Limbs(idxFactor1, idxFactor2);
    break;
  }
#else
  limb *ptrFactor1, *ptrFactor2;
  int prodCol, fact1Col;
  double dRangeLimb = (double)LIMB_RANGE;
  double dInvRangeLimb = 1 / dRangeLimb;
  int low = 0;              // Low limb of sums of multiplications.
  double dAccumulator = 0;  // Approximation to the sum of multiplications.
  int factor1, factor2;
  for (prodCol = 0; prodCol < 2 * nbrLen - 1; prodCol++)
  {    // Process each limb of product (least to most significant limb).
    if (prodCol < nbrLen)
    {   // Processing first half (least significant) of product.
      ptrFactor2 = &arr[idxFactor2 + prodCol];
      ptrFactor1 = &arr[idxFactor1];
      fact1Col = prodCol;
    }
    else
    {  // Processing second half (most significant) of product.
      ptrFactor2 = &arr[idxFactor2 + nbrLen - 1];
      ptrFactor1 = &arr[idxFactor1 + prodCol - nbrLen + 1];
      fact1Col = 2 * (nbrLen - 1) - prodCol;
    }
    for (; fact1Col>=0; fact1Col--)
    {
      factor1 = (ptrFactor1++)->x;
      factor2 = (ptrFactor2--)->x;
      low += factor1 * factor2;
      dAccumulator += (double)factor1 * (double)factor2;
    }
    low &= MAX_VALUE_LIMB;    // Trim extra bits.
    arrayAux[prodCol].x = low;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    if (low < HALF_INT_RANGE)
    {
      dAccumulator = floor(dAccumulator * dInvRangeLimb + 0.25);
    }
    else
    {
      dAccumulator = floor(dAccumulator * dInvRangeLimb - 0.25);
    }
    low = (unsigned int)(dAccumulator - floor(dAccumulator * dInvRangeLimb) * dRangeLimb);
  }
  arrayAux[prodCol].x = low;
#endif
  memcpy(&arr[idxFactor1], &arrayAux[0], 2 * nbrLen * sizeof(limb));
  return;
}

// Recursive Karatsuba function.
static void Karatsuba(int idxFactor1, int nbrLen, int diffIndex)
{
  int idxFactor2 = idxFactor1 + nbrLen;
  int i;
  unsigned int carry1First, carry1Second;
  unsigned int carry2Second;
  limb *ptrResult, *ptrHigh, tmp;
  int middle;
  int sign;
  int halfLength;
  if (nbrLen <= KARATSUBA_CUTOFF)
  {
    // Check if one of the factors is equal to zero.
    ptrResult = &arr[idxFactor1];
    for (i = nbrLen; i > 0; i--)
    {
      if ((ptrResult++)->x != 0)
      {
        break;
      }
    }
    if (i > 0)
    {     // First factor is not zero. Check second.
      ptrResult = &arr[idxFactor2];
      for (i = nbrLen; i > 0; i--)
      {
        if ((ptrResult++)->x != 0)
        {
          break;
        }
      }
    }
    if (i==0)
    {    // One of the factors is equal to zero.
      for (i = nbrLen - 1; i >= 0; i--)
      {
        arr[idxFactor1 + i].x = arr[idxFactor2 + i].x = 0;
      }
      return;
    }
         // Below cutoff: perform standard classical multiplcation.
    ClassicalMult(idxFactor1, idxFactor2, nbrLen);
    return;
  }
  // Length > KARATSUBA_CUTOFF: Use Karatsuba multiplication.
  // It uses three half-length multiplications instead of four.
  //  x*y = (xH*b + xL)*(yH*b + yL)
  //  x*y = (b + 1)*(xH*yH*b + xL*yL) + (xH - xL)*(yL - yH)*b
  // The length of b is stored in variable halfLength.
  // Since the absolute values of (xH - xL) and (yL - yH) fit in
  // a single limb, there will be no overflow.

  // At this moment the order is: xL, xH, yL, yH.
  // Exchange high part of first factor with low part of 2nd factor.
  halfLength = nbrLen >> 1;
  for (i = idxFactor1 + halfLength; i<idxFactor2; i++)
  {
    tmp.x = arr[i].x;
    arr[i].x = arr[i + halfLength].x;
    arr[i + halfLength].x = tmp.x;
  }
  // At this moment the order is: xL, yL, xH, yH.
  // Get absolute values of (xH-xL) and (yL-yH) and the signs.
  sign = absSubtract(idxFactor1, idxFactor2, diffIndex, halfLength);
  sign ^= absSubtract(idxFactor2 + halfLength, idxFactor1 + halfLength,
    diffIndex + halfLength, halfLength);
  middle = diffIndex;
  diffIndex += nbrLen;
  Karatsuba(idxFactor1, halfLength, diffIndex); // Multiply both low parts.
  Karatsuba(idxFactor2, halfLength, diffIndex); // Multiply both high parts.
  Karatsuba(middle, halfLength, diffIndex);     // Multiply the differences.
     // Process all carries at the end.
     // Obtain (b+1)(xH*yH*b + xL*yL) = xH*yH*b^2 + (xL*yL+xH*yH)*b + xL*yL
     // The first and last terms are already in correct locations.
  ptrResult = &arr[idxFactor1+halfLength];
  carry1First = carry1Second = carry2Second = 0;
  for (i = halfLength; i > 0; i--)
  {
    // The sum of three ints overflows an unsigned int variable,
    // so two adds are required. Also carries must be separated in
    // order to avoid overflow:
    // 00000001 + 7FFFFFFF + 7FFFFFFF = FFFFFFFF
    unsigned int accum1Lo = carry1First + ptrResult->x + (ptrResult + halfLength)->x;
    unsigned int accum2Lo;
    carry1First = accum1Lo >> BITS_PER_GROUP;
    accum2Lo = carry2Second + (accum1Lo & MAX_VALUE_LIMB) +
               (ptrResult - halfLength)->x;
    carry2Second = accum2Lo >> BITS_PER_GROUP;
    accum1Lo = carry1Second + (accum1Lo & MAX_VALUE_LIMB) +
               (ptrResult + nbrLen)->x;
    carry1Second = accum1Lo >> BITS_PER_GROUP;
    (ptrResult + halfLength)->x = accum1Lo & MAX_VALUE_LIMB;
    ptrResult->x = accum2Lo & MAX_VALUE_LIMB;
    ptrResult++;
  }
  (ptrResult + halfLength)->x += carry1First + carry1Second;
  ptrResult->x += carry1First + carry2Second;
  // Process carries.
  ptrResult = &arr[idxFactor1];
  carry1First = 0;
  for (i = 2*nbrLen; i > 0; i--)
  {
    carry1First += ptrResult->x;
    (ptrResult++)->x = carry1First & MAX_VALUE_LIMB;
    carry1First >>= BITS_PER_GROUP;
  }
  // Compute final product.
  ptrHigh = &arr[middle];
  ptrResult = &arr[idxFactor1 + halfLength];
  if (sign != 0)
  {            // (xH-xL) * (yL-yH) is negative.
    int borrow = 0;
    for (i = nbrLen; i > 0; i--)
    {
      borrow += ptrResult->x - (ptrHigh++)->x;
      (ptrResult++)->x = borrow & MAX_VALUE_LIMB;
      borrow >>= BITS_PER_GROUP;
    }
    for (i = halfLength; i > 0; i--)
    {
      borrow += ptrResult->x;
      (ptrResult++)->x = borrow & MAX_VALUE_LIMB;
      borrow >>= BITS_PER_GROUP;
    }
  }
  else
  {            // (xH-xL) * (yL-yH) is positive or zero.
    unsigned int carry = 0;
    for (i = nbrLen; i > 0; i--)
    {
      carry += (unsigned int)ptrResult->x + (unsigned int)(ptrHigh++)->x;
      (ptrResult++)->x = (int)(carry & MAX_VALUE_LIMB);
      carry >>= BITS_PER_GROUP;
    }
    for (i = halfLength; i > 0; i--)
    {
      carry += (unsigned int)ptrResult->x;
      (ptrResult++)->x = (int)(carry & MAX_VALUE_LIMB);
      carry >>= BITS_PER_GROUP;
    }
  }
}

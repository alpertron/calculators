//
// This file is part of Alpertron Calculators.
//
// Copyright 2015-2025 Dario Alejandro Alpern
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

/* Perform Karatsuba multiplication in little endian order */

#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include "bignbr.h"

#define KARATSUBA_CUTOFF 16
#define FFT_THRESHOLD 100
limb arr[MAX_LEN/*3* FFT_THRESHOLD*/];
limb arrayAux[MAX_LEN/*3* FFT_THRESHOLD*/];
static limb accumulatedProd[MAX_LEN_MULT];
static limb partialProd[MAX_LEN_MULT];
static int karatLength;
static void Karatsuba(int numLen);

void multiply(const limb* factor1, const limb* factor2, limb* result,
  int len, int* pResultLen)
{
  multiplyWithBothLen(factor1, factor2, result, len, len, pResultLen);
}

static inline void multiplyWithBothLenLL(const limb *factor1, const limb *factor2,
  limb *result, int len1, int len2, int *pResultLen)
{
  int length = len1;
  int lenBytes;
  // Compute the maximum length.
  if (length < len2)
  {
    length = len2;
  }
  if (length > FFT_THRESHOLD)
  {
    fftMultiplication(factor1, factor2, result, len1, len2, pResultLen);
    return;
  }
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
  lenBytes = 2 * length * (int)sizeof(limb);
  (void)memset(arr, 0, lenBytes);
  lenBytes = len1 * (int)sizeof(limb);
  (void)memcpy(&arr[0], factor1, lenBytes);
  lenBytes = len2 * (int)sizeof(limb);
  (void)memcpy(&arr[length], factor2, lenBytes);
  Karatsuba(length);
  lenBytes = 2 * length * (int)sizeof(limb);
  (void)memcpy(result, &arr[2 * (karatLength - length)], lenBytes);
  if (pResultLen != NULL)
  {
  	(void)memcpy(result, &arr[2 * (karatLength - length)], lenBytes);
    length = len1 + len2;
    while (length > 1 && (arr[length - 1].x == 0))
    {
      length--;
    }
    *pResultLen = length;
  }
}

void multiplyWithBothLen(const limb* factor1, const limb* factor2, limb* result,
  int len1, int len2, int* pResultLen)
{
  const limb* minFact;
  const limb* maxFact;
  int minLen;
  int maxLen;
  int lenProd = len1 + len2;
  int lenBytes;
  int offset;
  assert(len1 >= 1);
  assert(len2 >= 1);
  if (len1 < len2)
  {
    minFact = factor1;
    maxFact = factor2;
    minLen = len1;
    maxLen = len2;
  }
  else
  {
    minFact = factor2;
    maxFact = factor1;
    minLen = len2;
    maxLen = len1;
  }
  if ((minLen == 1) && (minFact->x == 0))
  {     // Multiply by zero. Set product to zero.
    lenBytes = lenProd * (int)sizeof(int);
    (void)memset(result, 0, lenBytes);
    if (pResultLen != NULL)
    {
      *pResultLen = 1;
    }
    return;
  }
  if (minLen == maxLen)
  {
    multiplyWithBothLenLL(minFact, maxFact, result, minLen, maxLen, pResultLen);
    return;
  }
  // Perform several multiplications and add all products.
  lenBytes = lenProd * (int)sizeof(int);
  (void)memset(accumulatedProd, 0, lenBytes);
  for (offset = 0; offset < (maxLen - minLen); offset += minLen)
  {
    multiplyWithBothLenLL(minFact, maxFact + offset, partialProd,
      minLen, minLen, NULL);
    // Add partial product to accumulated product.
    AddBigNbr(accumulatedProd + offset, partialProd,
      accumulatedProd + offset, 2 * minLen);
  }
  multiplyWithBothLenLL(minFact, maxFact + offset, partialProd,
    minLen, maxLen - offset, NULL);
  // Add partial product to accumulated product.
  AddBigNbr(accumulatedProd + offset, partialProd,
    accumulatedProd + offset, lenProd - offset);
  lenBytes = lenProd * (int)sizeof(int);
  (void)memcpy(result, accumulatedProd, lenBytes);
  // Copy accumulatedProd to result.
  while ((lenProd > 1) && (accumulatedProd[lenProd - 1].x == 0))
  {
    lenProd--;
  }
  if (pResultLen != NULL)
  {
    *pResultLen = lenProd;
  }
}

// The return value is the sign: true: negative.
// In result the absolute value of the difference is computed.
static int absSubtract(int idxMinuend, int idxSubtrahend,
                       int idxResult, int nbrLen)
{
  int indexMinuend = idxMinuend;
  int indexSubtrahend = idxSubtrahend;
  int sign = 0;
  unsigned int borrow;
  int i;
  limb *ptrArray;
  const limb* ptrMinuend = &arr[indexMinuend];
  const limb* ptrSubtrahend = &arr[indexSubtrahend];
  for (i = nbrLen-1; i>=0; i--)
  {
    if ((ptrMinuend + i)->x != (ptrSubtrahend + i)->x)
    {
      break;
    }
  }
  if ((i>=0) && ((ptrMinuend + i)->x < (ptrSubtrahend + i)->x))
  {
    sign = 1;
    i = indexMinuend;    // Exchange minuend and subtrahend.
    indexMinuend = indexSubtrahend;
    indexSubtrahend = i;
  }
  ptrArray = &arr[idxResult];
  indexMinuend -= idxResult;
  indexSubtrahend -= idxResult;
  borrow = 0U;
  for (i = nbrLen; i > 0; i -= 2)
  {
    borrow = (unsigned int)(ptrArray + indexMinuend)->x - (unsigned int)(ptrArray + indexSubtrahend)->x -
      (borrow >> BITS_PER_GROUP);
    ptrArray->x = UintToInt(borrow & MAX_VALUE_LIMB);
    ptrArray++;
    borrow = (unsigned int)(ptrArray + indexMinuend)->x - (unsigned int)(ptrArray + indexSubtrahend)->x -
      (borrow >> BITS_PER_GROUP);
    ptrArray->x = UintToInt(borrow & MAX_VALUE_LIMB);
    ptrArray++;
  }
  return sign;
}

struct stKaratsubaStack
{
  int idxFactor1;
  int sign;
  int stage;
};

static struct stKaratsubaStack astKaratsubaStack[16];

inline static void computeFinalProduct(int nbrLen, int idxFactor1, int sign, int diffIndex)
{
  int i;
  int halfLength = nbrLen / 2;
  // Process all carries at the end.
  // Obtain (b+1)(xH*yH*b + xL*yL) = xH*yH*b^2 + (xL*yL+xH*yH)*b + xL*yL
  // The first and last terms are already in correct locations.
  limb* ptrResult = &arr[idxFactor1 + halfLength];
  unsigned int carry1First = 0U;
  unsigned int carry1Second = 0U;
  unsigned int carry2Second = 0U;
  int doubleLength = 2 * nbrLen;
  const limb* ptrHigh;
  limb* ptrLimb;

  for (i = halfLength; i > 0; i--)
  {
    // The sum of three ints overflows an unsigned int variable,
    // so two adds are required. Also carries must be separated in
    // order to avoid overflow:
    // 00000001 + 7FFFFFFF + 7FFFFFFF = FFFFFFFF
    unsigned int accum1Lo = carry1First + (unsigned int)ptrResult->x +
      (unsigned int)(ptrResult + halfLength)->x;
    unsigned int accum2Lo;
    carry1First = accum1Lo >> BITS_PER_GROUP;
    accum2Lo = carry2Second + (accum1Lo & MAX_VALUE_LIMB) +
      (unsigned int)(ptrResult - halfLength)->x;
    carry2Second = accum2Lo >> BITS_PER_GROUP;
    accum1Lo = carry1Second + (accum1Lo & MAX_VALUE_LIMB) +
      (unsigned int)(ptrResult + nbrLen)->x;
    carry1Second = accum1Lo >> BITS_PER_GROUP;
    (ptrResult + halfLength)->x = UintToInt(accum1Lo & MAX_VALUE_LIMB);
    ptrResult->x = UintToInt(accum2Lo & MAX_VALUE_LIMB);
    ptrResult++;
  }
  // Process carries.
  ptrLimb = ptrResult + halfLength;
  carry1Second += carry1First + (unsigned int)ptrLimb->x;
  ptrLimb->x = (int)carry1Second & MAX_INT_NBR;
  for (i = nbrLen + halfLength + 1; i < doubleLength; i++)
  {
    if (carry1Second < MAX_VALUE_LIMB)
    {
      break;
    }
    carry1Second >>= BITS_PER_GROUP;
    ptrLimb++;
    carry1Second += (unsigned int)ptrLimb->x;
    ptrLimb->x = (int)carry1Second & MAX_INT_NBR;
  }
  carry2Second += carry1First + (unsigned int)ptrResult->x;
  ptrResult->x = (int)carry2Second & MAX_INT_NBR;
  for (i = nbrLen + 1; i < doubleLength; i++)
  {
    if (carry2Second < MAX_VALUE_LIMB)
    {
      break;
    }
    carry2Second >>= BITS_PER_GROUP;
    ptrResult++;
    carry2Second += (unsigned int)ptrResult->x;
    ptrResult->x = (int)carry2Second & MAX_INT_NBR;
  }
  // Compute final product.
  ptrHigh = &arr[diffIndex];
  ptrResult = &arr[idxFactor1 + halfLength];
  if (sign != 0)
  {            // (xH-xL) * (yL-yH) is negative.
    unsigned int borrow = 0U;
    for (i = halfLength; i > 0; i--)
    {
      borrow = (unsigned int)ptrResult->x - (unsigned int)ptrHigh->x - borrow;
      ptrHigh++;
      ptrResult->x = UintToInt(borrow & MAX_VALUE_LIMB);
      ptrResult++;
      borrow = (unsigned int)ptrResult->x - (unsigned int)ptrHigh->x -
        (borrow >> BITS_PER_GROUP);
      ptrHigh++;
      ptrResult->x = UintToInt(borrow & MAX_VALUE_LIMB);
      ptrResult++;
      borrow >>= BITS_PER_GROUP;
    }
    for (i = halfLength; i > 0; i--)
    {
      if (borrow == 0U)
      {
        break;
      }
      borrow = (unsigned int)ptrResult->x - borrow;
      ptrResult->x = UintToInt(borrow & MAX_VALUE_LIMB);
      ptrResult++;
      borrow >>= BITS_PER_GROUP;
    }
  }
  else
  {            // (xH-xL) * (yL-yH) is positive or zero.
    unsigned int carry = 0;
    for (i = halfLength; i > 0; i--)
    {
      carry += (unsigned int)ptrResult->x + (unsigned int)ptrHigh->x;
      ptrHigh++;
      ptrResult->x = UintToInt(carry & MAX_VALUE_LIMB);
      ptrResult++;
      carry = (carry >> BITS_PER_GROUP) +
        (unsigned int)ptrResult->x + (unsigned int)ptrHigh->x;
      ptrHigh++;
      ptrResult->x = UintToInt(carry & MAX_VALUE_LIMB);
      ptrResult++;
      carry >>= BITS_PER_GROUP;
    }
    for (i = halfLength; i > 0; i--)
    {
      if (carry == 0U)
      {
        break;
      }
      carry += (unsigned int)ptrResult->x;
      ptrResult->x = UintToInt(carry & MAX_VALUE_LIMB);
      ptrResult++;
      carry >>= BITS_PER_GROUP;
    }
  }
}

inline static void belowKaratsubaCutoff(int nbrLen, int idxFactor1, int idxFactor2)
{
  const limb* ptrResult;
  int i;
  // Check if one of the factors is equal to zero.
  ptrResult = &arr[idxFactor1];
  for (i = nbrLen; i > 0; i--)
  {
    if (ptrResult->x != 0)
    {
      break;
    }
    ptrResult++;
  }
  if (i > 0)
  {     // First factor is not zero. Check second.
    ptrResult = &arr[idxFactor2];
    for (i = nbrLen; i > 0; i--)
    {
      if (ptrResult->x != 0)
      {
        break;
      }
      ptrResult++;
    }
  }
  if (i == 0)
  {    // One of the factors is equal to zero.
    for (i = nbrLen - 1; i >= 0; i--)
    {
      arr[idxFactor1 + i].x = 0;
      arr[idxFactor2 + i].x = 0;
    }
  }
  else
  {   // Below cutoff: perform standard classical multiplication.
    ClassicalMult(idxFactor1, idxFactor2, nbrLen);
  }
}

static void Karatsuba(int numLen)
{
  int nbrLen = numLen;
  int idxFactor1 = 0;
  int idxFactor2;
  limb* ptrLimb;
  const limb* endPtrLimb;
  int sign = 0;
  int halfLength;
  int diffIndex = 2 * nbrLen;
  static struct stKaratsubaStack *pstKaratsubaStack = astKaratsubaStack;
  int stage = 0;
  // Save current parameters in stack.
  pstKaratsubaStack->idxFactor1 = idxFactor1;
  pstKaratsubaStack->sign = 0;
  pstKaratsubaStack->stage = -1;
  pstKaratsubaStack++;
  do
  {
    switch (stage)
    {
    case 0:
      idxFactor2 = idxFactor1 + nbrLen;
      if (nbrLen <= KARATSUBA_CUTOFF)
      {
        belowKaratsubaCutoff(nbrLen, idxFactor1, idxFactor2);
        pstKaratsubaStack--;
        assert(pstKaratsubaStack >= &astKaratsubaStack[0]);
        idxFactor1 = pstKaratsubaStack->idxFactor1;
        nbrLen *= 2;
        diffIndex -= nbrLen;
        stage = pstKaratsubaStack->stage;
        sign = pstKaratsubaStack->sign;
        break;
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
      halfLength = nbrLen / 2;
      ptrLimb = &arr[idxFactor1 + halfLength];
      endPtrLimb = &arr[idxFactor2];
      for (; ptrLimb < endPtrLimb; ptrLimb++)
      {
        limb tmp;
        tmp.x = ptrLimb->x;
        ptrLimb->x = (ptrLimb + halfLength)->x;
        (ptrLimb + halfLength)->x = tmp.x;
      }
      
      // At this moment the order is: xL, yL, xH, yH.
      // Get absolute values of (xH-xL) and (yL-yH) and the signs.
      sign = absSubtract(idxFactor1, idxFactor2, diffIndex, halfLength);
      sign ^= absSubtract(idxFactor2 + halfLength, idxFactor1 + halfLength,
        diffIndex + halfLength, halfLength);
      // Save current parameters in stack.
      pstKaratsubaStack->idxFactor1 = idxFactor1;
      pstKaratsubaStack->sign = sign;
      pstKaratsubaStack->stage = 1;
      pstKaratsubaStack++;
      // Multiply both low parts.
      diffIndex += nbrLen;
      nbrLen = halfLength;
      break;
    case 1:
      // Multiply both high parts.
      idxFactor1 += nbrLen;
      diffIndex += nbrLen;
      nbrLen /= 2;
      pstKaratsubaStack->stage = 2;
      pstKaratsubaStack++;
      stage = 0;         // Start new Karatsuba multiplication.
      break;
    case 2:
      // Multiply the differences.
      idxFactor1 = diffIndex;
      diffIndex += nbrLen;
      nbrLen >>= 1;
      pstKaratsubaStack->stage = 3;
      pstKaratsubaStack++;
      stage = 0;         // Start new Karatsuba multiplication.
      break;
    default:
      computeFinalProduct(nbrLen, idxFactor1, sign, diffIndex);
      nbrLen *= 2;
      diffIndex -= nbrLen;
      pstKaratsubaStack--;
      assert(pstKaratsubaStack >= &astKaratsubaStack[0]);
      idxFactor1 = pstKaratsubaStack->idxFactor1;
      stage = pstKaratsubaStack->stage;
      sign = pstKaratsubaStack->sign;
      break;
    }     // End switch
  } while (stage >= 0);
}

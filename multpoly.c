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

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "bignbr.h"
#include "highlevel.h"
#include "polynomial.h"
#include "showtime.h"

#define KARATSUBA_POLY_CUTOFF      16
#define SQRT_KARATSUBA_POLY_CUTOFF 4

struct stKaratsubaStack
{
  int idxFactor1;
  int stage;
};

//#define ADJUST_MODULUS(sum, modulus) sum += modulus & (sum >> 31);
#define ADJUST_MODULUS(sum, modulus) if (sum < 0) {sum += modulus;}
static struct stKaratsubaStack astKaratsubaStack[10];
static BigInteger coeff[2 * KARATSUBA_POLY_CUTOFF];
int polyInv[COMPRESSED_POLY_MAX_LENGTH];
int polyInvCached;
static int polyMultM[COMPRESSED_POLY_MAX_LENGTH];
static int polyMultT[COMPRESSED_POLY_MAX_LENGTH];
extern int poly4[COMPRESSED_POLY_MAX_LENGTH];
extern int poly5[COMPRESSED_POLY_MAX_LENGTH];

// Multiply two groups of nbrLen coefficients. The first one starts at
// idxFactor1 and the second one at idxFactor2. The 2*nbrLen coefficient
// result is stored starting at idxFactor1. Use arrAux as temporary storage.
// Accumulate products by result coefficient.
static void ClassicalPolyMult(int idxFactor1, int idxFactor2, int coeffLen, int nbrLimbs)
{
  int i;
  int j;
  int* ptrFactor1;
  const int* ptrFactor2;
#ifndef _USING64BITS_
  limb result;
#endif
  for (i = 0; i < (2 * coeffLen) - 1; i++)
  {    // Process each limb of product (least to most significant limb).
    if (i < coeffLen)
    {   // Processing first half (least significant) of product.
      ptrFactor2 = &polyMultTemp[(idxFactor2 + i) * nbrLimbs];
      ptrFactor1 = &polyMultTemp[idxFactor1 * nbrLimbs];
      j = i;
    }
    else
    {  // Processing second half (most significant) of product.
      ptrFactor2 = &polyMultTemp[(idxFactor2 + coeffLen - 1) * nbrLimbs];
      ptrFactor1 = &polyMultTemp[(idxFactor1 + i - coeffLen + 1) * nbrLimbs];
      j = 2 * (coeffLen - 1) - i;
    }
    (void)memset(coeff[i].limbs, 0, nbrLimbs * sizeof(limb));
    if (nbrLimbs == 2)
    {    // Optimization for the case when there is only one limb.
      int modulus = TestNbr[0].x;
      int sum = 0;
#ifdef _USING64BITS_
      uint64_t dSum;
#endif
      ptrFactor1++;
      ptrFactor2++;
      if (modulus < (32768 / SQRT_KARATSUBA_POLY_CUTOFF))
      {          // Sum of products fits in one limb.
        for (; j >= 3; j -= 4)
        {
          sum += (*ptrFactor1 * *ptrFactor2) +
            (*(ptrFactor1 + 2) * *(ptrFactor2 - 2)) +
            (*(ptrFactor1 + 4) * *(ptrFactor2 - 4)) +
            (*(ptrFactor1 + 6) * *(ptrFactor2 - 6));
          ptrFactor1 += 8;
          ptrFactor2 -= 8;
        }
        while (j >= 0)
        {
          sum += (*ptrFactor1 * *ptrFactor2);
          ptrFactor1 += 2;
          ptrFactor2 -= 2;
          j--;
        }
        sum %= modulus;
      }
      else if (modulus <= 32768)
      {         // Product fits in one limb.
#ifdef _USING64BITS_
        dSum = 0;
#else
        double dSum = 0;
#endif
        for (; j >= 3; j -= 4)
        {
#ifdef _USING64BITS_
          dSum += (uint64_t)(*ptrFactor1 * *ptrFactor2) +
            (uint64_t)(*(ptrFactor1 + 2) * *(ptrFactor2 - 2)) +
            (uint64_t)(*(ptrFactor1 + 4) * *(ptrFactor2 - 4)) +
            (uint64_t)(*(ptrFactor1 + 6) * *(ptrFactor2 - 6));
#else
          dSum += (double)(*ptrFactor1 * *ptrFactor2) +
            (double)(*(ptrFactor1 + 2) * *(ptrFactor2 - 2)) +
            (double)(*(ptrFactor1 + 4) * *(ptrFactor2 - 4)) +
            (double)(*(ptrFactor1 + 6) * *(ptrFactor2 - 6));
#endif
          ptrFactor1 += 8;
          ptrFactor2 -= 8;
        }
        while (j >= 0)
        {
#ifdef _USING64BITS_
          dSum += (uint64_t)(*ptrFactor1 * *ptrFactor2);
#else
          dSum += (double)(*ptrFactor1 * *ptrFactor2);
#endif
          ptrFactor1 += 2;
          ptrFactor2 -= 2;
          j--;
        }
#ifdef _USING64BITS_
        sum = (int)(dSum % modulus);
#else
        sum = (int)(dSum - floor(dSum / modulus) * modulus);
#endif
      }
      else
      {
#ifdef _USING64BITS_
        dSum = 0;
#endif
        for (; j >= 0; j--)
        {
#ifdef _USING64BITS_
          dSum += (int64_t)*ptrFactor1 * *ptrFactor2;
          if ((int64_t)dSum < 0)
          {
            dSum %= modulus;
          }
#else
          smallmodmult(*ptrFactor1, *ptrFactor2, &result, modulus);
          sum += result.x - modulus;
          // If sum < 0 do sum <- sum + modulus else do nothing.
          ADJUST_MODULUS(sum, modulus);
#endif
          ptrFactor1 += 2;
          ptrFactor2 -= 2;
        }
#ifdef _USING64BITS_
        sum = (int)(dSum % modulus);
#endif
      }
      coeff[i].limbs[0].x = sum;
    }
    else
    {          // General case.
      for (; j >= 0; j--)
      {
        LenAndLimbs2ArrLimbs(ptrFactor1, operand3.limbs, nbrLimbs);
        LenAndLimbs2ArrLimbs(ptrFactor2, operand2.limbs, nbrLimbs);
        modmult(operand3.limbs, operand2.limbs, operand3.limbs);
        AddBigNbrMod(coeff[i].limbs, operand3.limbs, coeff[i].limbs);
        ptrFactor1 += nbrLimbs;
        ptrFactor2 -= nbrLimbs;
      }
    }
  }
  ptrFactor1 = &polyMultTemp[idxFactor1 * nbrLimbs];
  if (nbrLimbs == 2)
  {    // Optimization for the case when there is only one limb.
    for (i = 0; i < (2 * coeffLen) - 1; i++)
    {
      *ptrFactor1 = 1;
      ptrFactor1++;
      *ptrFactor1 = coeff[i].limbs[0].x;
      ptrFactor1++;
    }
  }
  else
  {
    for (i = 0; i < (2 * coeffLen) - 1; i++)
    {
      ArrLimbs2LenAndLimbs(ptrFactor1, coeff[i].limbs, nbrLimbs);
      ptrFactor1 += nbrLimbs;
    }
  }
  *ptrFactor1 = 1;
  *(ptrFactor1 + 1) = 0;
  return;
}

// Recursive Karatsuba function.
static void KaratsubaPoly(int idxFactor1, int nbrLen, int nbrLimbs)
{
  int i;
  int nbrLength = nbrLen;
  int idxFactor2;
  int* ptrResult;
  const int* ptrHigh;
  const int* ptr1;
  const int* ptr2;
  int sum;
  int modulus;
  int halfLength;
  int diffIndex = 2 * nbrLength;
  static struct stKaratsubaStack* pstKaratsubaStack = astKaratsubaStack;
  static int coeffic[MAX_LEN];
  int stage = 0;
  // Save current parameters in stack.
  pstKaratsubaStack->idxFactor1 = idxFactor1;
  pstKaratsubaStack->stage = -1;
  pstKaratsubaStack++;
  do
  {
    switch (stage)
    {
    case 0:
      idxFactor2 = idxFactor1 + nbrLength;
      if (nbrLength <= KARATSUBA_POLY_CUTOFF)
      {
        // Check if one of the factors is equal to zero.
        ptrResult = &polyMultTemp[idxFactor1 * nbrLimbs];
        i = nbrLength;
        if (nbrLimbs == 2)
        {
          ptrResult++;
          i--;
          while (((i & 0x80000000) | *ptrResult) == 0)
          {   // Loop not finished and coefficient is not zero.
            ptrResult += nbrLimbs;
            i--;
          }
        }
        else
        {
          i--;
          while (((i & 0x80000000) | (*ptrResult - 1) | *(ptrResult + 1)) == 0)
          {   // Loop not finished and coefficient is not zero.
            ptrResult += nbrLimbs;
            i--;
          }
        }
        if (i < 0)
        {      // First factor is zero. Initialize second to zero.
          ptrResult = &polyMultTemp[idxFactor2 * nbrLimbs];
          if (nbrLimbs == 2)
          {
            ptrResult++;
            for (i = nbrLength; i > 0; i--)
            {
              *ptrResult = 0;
              ptrResult += nbrLimbs;
            }
          }
          else
          {
            for (i = nbrLength; i > 0; i--)
            {
              *ptrResult = 1;
              *(ptrResult + 1) = 0;
              ptrResult += nbrLimbs;
            }
          }
        }
        else
        {     // First factor is not zero. Check second.
          ptrResult = &polyMultTemp[idxFactor2 * nbrLimbs];
          i = nbrLength;
          if (nbrLimbs == 2)
          {
            ptrResult++;
            i--;
            while (((i & 0x80000000) | *ptrResult) == 0)
            {   // Loop not finished and coefficient is not zero.
              ptrResult += nbrLimbs;
              i--;
            }
          }
          else
          {
            i--;
            while (((--i & 0x80000000) | (*ptrResult - 1) | *(ptrResult + 1)) == 0)
            {   // Loop not finished and coefficient is not zero.
              ptrResult += nbrLimbs;
              i--;
            }
          }
          if (i < 0)
          {    // Second factor is zero. Initialize first to zero.
            ptrResult = &polyMultTemp[idxFactor1 * nbrLimbs];
            if (nbrLimbs == 2)
            {
              ptrResult++;
              for (i = nbrLength; i > 0; i--)
              {
                *ptrResult = 0;
                ptrResult += nbrLimbs;
              }
            }
            else
            {
              for (i = nbrLength; i > 0; i--)
              {
                *ptrResult = 1;
                *(ptrResult + 1) = 0;
                ptrResult += nbrLimbs;
              }
            }
          }
          else
          {   // Both factors not zero: perform standard classical polynomial multiplcation.
            ClassicalPolyMult(idxFactor1, idxFactor2, nbrLength, nbrLimbs);
          }
        }
        pstKaratsubaStack--;
        idxFactor1 = pstKaratsubaStack->idxFactor1;
        nbrLength *= 2;
        diffIndex -= nbrLength;
        stage = pstKaratsubaStack->stage;
        break;
      }
      // Length > KARATSUBA_CUTOFF: Use Karatsuba multiplication.
      // It uses three half-length multiplications instead of four.
      //  x*y = (xH*b + xL)*(yH*b + yL)
      //  x*y = (b + 1)*(xH*yH*b + xL*yL) + (xH - xL)*(yL - yH)*b
      // The length of b is stored in variable halfLength.

      // At this moment the order is: xL, xH, yL, yH.
      // Exchange high part of first factor with low part of 2nd factor.
      halfLength = nbrLength / 2;
      int* ptrHighFirstFactor = &polyMultTemp[(idxFactor1 + halfLength) * nbrLimbs];
      int* ptrLowSecondFactor = &polyMultTemp[idxFactor2 * nbrLimbs];
      if (nbrLimbs == 2)
      {                        // Coefficients have only one limb.
        ptrHighFirstFactor++;  // Point to limb to be exchanged.
        ptrLowSecondFactor++;  // Point to limb to be exchanged.
        for (i = 0; i < halfLength; i++)
        {
          int coefficient = *ptrHighFirstFactor;
          *ptrHighFirstFactor = *ptrLowSecondFactor;
          *ptrLowSecondFactor = coefficient;
          ptrHighFirstFactor += 2;
          ptrLowSecondFactor += 2;
        }
      }
      else
      {
        int sizeCoeffInBytes = nbrLimbs * sizeof(int);
        for (i = 0; i < halfLength; i++)
        {
          (void)memcpy(coeffic, ptrHighFirstFactor, sizeCoeffInBytes);
          (void)memcpy(ptrHighFirstFactor, ptrLowSecondFactor, sizeCoeffInBytes);
          (void)memcpy(ptrLowSecondFactor, coeffic, sizeCoeffInBytes);
          ptrHighFirstFactor += nbrLimbs;
          ptrLowSecondFactor += nbrLimbs;
        }
      }
      // At this moment the order is: xL, yL, xH, yH.
      // Compute (xH-xL) and (yL-yH) and store them starting from index diffIndex.
      ptr1 = &polyMultTemp[idxFactor1 * nbrLimbs];
      ptr2 = &polyMultTemp[idxFactor2 * nbrLimbs];
      ptrResult = &polyMultTemp[diffIndex * nbrLimbs];
      if (nbrLimbs == 2)
      {    // Small modulus.
        modulus = TestNbr[0].x;
        ptr1++;
        ptr2++;
        ptrResult++;
        for (i = 0; i < halfLength; i++)
        {
          sum = *ptr2 - *ptr1;
          // If sum < 0 do sum <- sum + modulus else do nothing.
          ADJUST_MODULUS(sum, modulus);
          *ptrResult = sum;
          ptr1 += 2;
          ptr2 += 2;
          ptrResult += 2;
        }
        for (i = 0; i < halfLength; i++)
        {
          sum = *ptr1 - *ptr2;
          // If sum < 0 do sum <- sum + modulus else do nothing.
          ADJUST_MODULUS(sum, modulus);
          *ptrResult = sum;
          ptr1 += 2;
          ptr2 += 2;
          ptrResult += 2;
        }
      }
      else
      {    // General case.
        for (i = 0; i < halfLength; i++)
        {
          LenAndLimbs2ArrLimbs(ptr1, operand3.limbs, nbrLimbs);
          LenAndLimbs2ArrLimbs(ptr2, operand2.limbs, nbrLimbs);
          SubtBigNbrMod(operand2.limbs, operand3.limbs, operand3.limbs);
          ArrLimbs2LenAndLimbs(ptrResult, operand3.limbs, nbrLimbs);
          ptr1 += nbrLimbs;
          ptr2 += nbrLimbs;
          ptrResult += nbrLimbs;
        }
        for (i = 0; i < halfLength; i++)
        {
          LenAndLimbs2ArrLimbs(ptr1, operand3.limbs, nbrLimbs);
          LenAndLimbs2ArrLimbs(ptr2, operand2.limbs, nbrLimbs);
          SubtBigNbrMod(operand3.limbs, operand2.limbs, operand3.limbs);
          ArrLimbs2LenAndLimbs(ptrResult, operand3.limbs, nbrLimbs);
          ptr1 += nbrLimbs;
          ptr2 += nbrLimbs;
          ptrResult += nbrLimbs;
        }
      }
      // Save current parameters in stack.
      pstKaratsubaStack->idxFactor1 = idxFactor1;
      pstKaratsubaStack->stage = 1;
      pstKaratsubaStack++;
      // Multiply both low parts.
      diffIndex += nbrLength;
      nbrLength = halfLength;
      break;
    case 1:
      // Multiply both high parts.
      idxFactor1 += nbrLength;
      diffIndex += nbrLength;
      nbrLength >>= 1;
      pstKaratsubaStack->stage = 2;
      pstKaratsubaStack++;
      stage = 0;         // Start new Karatsuba multiplication.
      break;
    case 2:
      // Multiply the differences.
      idxFactor1 = diffIndex;
      diffIndex += nbrLength;
      nbrLength >>= 1;
      pstKaratsubaStack->stage = 3;
      pstKaratsubaStack++;
      stage = 0;         // Start new Karatsuba multiplication.
      break;
    default:
      halfLength = nbrLength >> 1;
      // Obtain (b+1)(xH*yH*b + xL*yL) = xH*yH*b^2 + (xL*yL+xH*yH)*b + xL*yL
      // The first and last terms are already in correct locations.
      // Add (xL*yL+xH*yH)*b.
      ptrResult = &polyMultTemp[(idxFactor1 + halfLength) * nbrLimbs];
      if (nbrLimbs == 2)
      {        // Optimization for small numbers.
        int nbrLen2 = nbrLength * 2;
        modulus = TestNbr[0].x;
        ptrResult++;
        for (i = halfLength; i > 0; i--)
        {
          // First addend is the coefficient from xH*yH*b^2 + xL*yL
          // Second addend is the coefficient from xH*yH
          int coeffTmp = *ptrResult + *(ptrResult + nbrLength) - modulus;
          // If sum < 0 do sum <- sum + modulus else do nothing.
          ADJUST_MODULUS(coeffTmp, modulus);
          // Third addend is xL*yL. Add all three coefficients.
          sum = coeffTmp + *(ptrResult - nbrLength) - modulus;
          // If sum < 0 do sum <- sum + modulus else do nothing.
          ADJUST_MODULUS(sum, modulus);
          // Store coefficient of xH*yH*b^2 + (xL*yL+xH*yH)*b + xL*yL
          *ptrResult = sum;

          // First addend is the coefficient from xL*yL
          // Second addend is the coefficient from xH*yH
          sum = coeffTmp + *(ptrResult + nbrLen2) - modulus;
          // If sum < 0 do sum <- sum + modulus else do nothing.
          ADJUST_MODULUS(sum, modulus);
          // Store coefficient of xH*yH*b^2 + (xL*yL+xH*yH)*b + xL*yL
          *(ptrResult + nbrLength) = sum;
          // Point to next address.
          ptrResult += 2;
        }
      }
      else
      {        // General case.
        for (i = halfLength; i > 0; i--)
        {
          // Obtain coefficient from xH*yH*b^2 + xL*yL
          LenAndLimbs2ArrLimbs(ptrResult, operand3.limbs, nbrLimbs);
          // Obtain coefficient from xH*yH
          LenAndLimbs2ArrLimbs(ptrResult + (halfLength * nbrLimbs), operand1.limbs, nbrLimbs);
          // Add these coefficients.
          AddBigNbrMod(operand3.limbs, operand1.limbs, operand3.limbs);
          // Obtain coefficient from xL*yL
          LenAndLimbs2ArrLimbs(ptrResult - (halfLength * nbrLimbs), operand2.limbs, nbrLimbs);
          AddBigNbrMod(operand2.limbs, operand3.limbs, operand2.limbs);
          // Store coefficient of xH*yH*b^2 + (xL*yL+xH*yH)*b + xL*yL
          ArrLimbs2LenAndLimbs(ptrResult, operand2.limbs, nbrLimbs);
          // Obtain coefficient from xH*yH
          LenAndLimbs2ArrLimbs(ptrResult + (nbrLength * nbrLimbs), operand2.limbs, nbrLimbs);
          // Add coefficient from xL*yL
          AddBigNbrMod(operand3.limbs, operand2.limbs, operand3.limbs);
          // Store coefficient of xH*yH*b^2 + (xL*yL+xH*yH)*b + xL*yL
          ArrLimbs2LenAndLimbs(ptrResult + (halfLength * nbrLimbs), operand3.limbs, nbrLimbs);
          // Point to next address.
          ptrResult += nbrLimbs;
        }
      }
      // Compute final product by adding (xH - xL)*(yL - yH)*b.
      ptrHigh = &polyMultTemp[diffIndex * nbrLimbs];
      ptrResult = &polyMultTemp[(idxFactor1 + halfLength) * nbrLimbs];
      if (nbrLimbs == 2)
      {        // Optimization for small numbers.
        modulus = TestNbr[0].x;
        ptrResult++;           // Point to limb to process.
        ptrHigh++;             // Point to limb to process.
        for (i = nbrLength; i >= 2; i -= 2)
        {
          sum = *ptrResult + *ptrHigh - modulus;
          // If sum < 0 do sum <- sum + modulus else do nothing.
          ADJUST_MODULUS(sum, modulus);
          *ptrResult = sum;
          sum = *(ptrResult + 2) + *(ptrHigh + 2) - modulus;
          // If sum < 0 do sum <- sum + modulus else do nothing.
          ADJUST_MODULUS(sum, modulus);
          *(ptrResult + 2) = sum;
          ptrHigh += 4;
          ptrResult += 4;
        }
        if (i > 0)
        {
          sum = *ptrResult + *ptrHigh - modulus;
          // If sum < 0 do sum <- sum + modulus else do nothing.
          ADJUST_MODULUS(sum, modulus);
          *ptrResult = sum;
        }
      }
      else
      {        // General case.
        for (i = nbrLength; i > 0; i--)
        {
          LenAndLimbs2ArrLimbs(ptrResult, operand3.limbs, nbrLimbs);
          LenAndLimbs2ArrLimbs(ptrHigh, operand2.limbs, nbrLimbs);
          AddBigNbrMod(operand3.limbs, operand2.limbs, operand3.limbs);
          ArrLimbs2LenAndLimbs(ptrResult, operand3.limbs, nbrLimbs);
          ptrHigh += nbrLimbs;
          ptrResult += nbrLimbs;
        }
      }
      nbrLength *= 2;
      diffIndex -= nbrLength;
      pstKaratsubaStack--;
      idxFactor1 = pstKaratsubaStack->idxFactor1;
      stage = pstKaratsubaStack->stage;
    }     // End switch
  } while (stage >= 0);
}

// Multiply factor1 by factor2.The result will be stored in polyMultTemp.
static void MultIntegerPolynomial(int degree1, int degree2,
  /*@in@*/int* factor1, /*@in@*/int* factor2)
{
  int indexes[2][MAX_DEGREE + 1];
  int* ptrIndex;
  int* piTemp;
  int* piDest;
  int currentDegree;
  int index;
  int degreeF2;

  if (degree1 < degree2)
  {       // Force degree1 >= degree2.
    int tmp = degree1;
    degree1 = degree2;
    degree2 = tmp;
    piTemp = factor1;
    factor1 = factor2;
    factor2 = piTemp;
  }
  // Fill indexes to start of each coefficient.
  ptrIndex = &indexes[0][0];
  index = 0;
  for (currentDegree = 0; currentDegree <= degree1; currentDegree++)
  {
    *ptrIndex = index;
    ptrIndex++;
    index += numLimbs(factor1 + index) + 1;
  }
  ptrIndex = &indexes[1][0];
  index = 0;
  for (currentDegree = 0; currentDegree <= degree2; currentDegree++)
  {
    *ptrIndex = index;
    ptrIndex++;
    index += numLimbs(factor2 + index) + 1;
  }
  piDest = polyMultTemp;
  for (currentDegree = 0; currentDegree < degree2; currentDegree++)
  {
    intToBigInteger(&operand4, 0);
    for (degreeF2 = 0; degreeF2 <= currentDegree; degreeF2++)
    {
      UncompressBigIntegerB(factor2 + indexes[1][degreeF2], &operand1);
      UncompressBigIntegerB(factor1 + indexes[0][currentDegree - degreeF2], &operand2);
      (void)BigIntMultiply(&operand1, &operand2, &operand3);
      BigIntAdd(&operand4, &operand3, &operand4);
      NumberLength = operand4.nbrLimbs;
    }
    BigInteger2IntArray(piDest, &operand4);
    piDest += 1 + numLimbs(piDest);
  }
  for (; currentDegree <= degree1; currentDegree++)
  {
    intToBigInteger(&operand4, 0);
    for (degreeF2 = 0; degreeF2 <= degree2; degreeF2++)
    {
      UncompressBigIntegerB(factor2 + indexes[1][degreeF2], &operand1);
      UncompressBigIntegerB(factor1 + indexes[0][currentDegree - degreeF2], &operand2);
      (void)BigIntMultiply(&operand1, &operand2, &operand3);
      BigIntAdd(&operand4, &operand3, &operand4);
      NumberLength = operand4.nbrLimbs;
    }
    BigInteger2IntArray(piDest, &operand4);
    piDest += 1 + numLimbs(piDest);
  }
  for (; currentDegree <= (degree1 + degree2); currentDegree++)
  {
    intToBigInteger(&operand4, 0);
    for (degreeF2 = currentDegree - degree1; degreeF2 <= degree2; degreeF2++)
    {
      UncompressBigIntegerB(factor2 + indexes[1][degreeF2], &operand1);
      UncompressBigIntegerB(factor1 + indexes[0][currentDegree - degreeF2], &operand2);
      (void)BigIntMultiply(&operand1, &operand2, &operand3);
      BigIntAdd(&operand4, &operand3, &operand4);
      NumberLength = operand4.nbrLimbs;
    }
    BigInteger2IntArray(piDest, &operand4);
    piDest += 1 + numLimbs(piDest);
  }
}

// Multiply factor1 by factor2. The result will be stored in polyMultTemp.
void MultPolynomial(int degree1, int degree2, /*@in@*/int* factor1, /*@in@*/int* factor2)
{
  int* ptrValue1;
  int nbrLimbs;
  int karatDegree;
  if (modulusIsZero)
  {
    MultIntegerPolynomial(degree1, degree2, factor1, factor2);
    return;
  }
  nbrLimbs = NumberLength + 1;
  if (((degree1 * degree1) < degree2) || ((degree2 * degree2) < degree1))
  {    // One of the factors is a lot smaller than the other.
       // Use classical multiplication of polynomials.
    const int *ptrSrc1;
    const int *ptrSrc2;
    int degreeProd = degree1 + degree2;
    int currentDegree1;
    int currentDegree2;
    // Initialize product to zero.
    int* ptrDest = polyMultTemp;
    for (currentDegree1 = 0; currentDegree1 <= degreeProd; currentDegree1++)
    {
      *ptrDest = 1;
      *(ptrDest + 1) = 0;
      ptrDest += nbrLimbs;
    }
    ptrSrc1 = factor1;
    for (currentDegree1 = 0; currentDegree1 <= degree1; currentDegree1++)
    {
      if ((*ptrSrc1 != 1) || *(ptrSrc1 + 1) != 0)
      {       // Only process factor if it is not zero.
        IntArray2BigInteger(ptrSrc1, &operand3);
        ptrSrc2 = factor2;
        ptrDest = &polyMultTemp[currentDegree1 * nbrLimbs];
        if ((NumberLength == 1) && (TestNbr[0].x <= 32768))
        {
          int mod = TestNbr[0].x;
          ptrSrc2++;
          ptrDest++;
          for (currentDegree2 = 0; currentDegree2 <= degree2; currentDegree2++)
          {
            *ptrDest = (*ptrDest + (operand3.limbs[0].x * *ptrSrc2)) % mod;
            ptrSrc2 += 2;
            ptrDest += 2;
          }
        }
        else
        {
          for (currentDegree2 = 0; currentDegree2 <= degree2; currentDegree2++)
          {
            if ((*ptrSrc2 != 1) || (*(ptrSrc2 + 1) != 0))
            {       // Only process factor if it is not zero.
              IntArray2BigInteger(ptrSrc2, &operand2);
              modmult(operand2.limbs, operand3.limbs, operand2.limbs);
              IntArray2BigInteger(ptrDest, &operand1);
              AddBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
              BigInteger2IntArray(ptrDest, &operand1);
            }
            ptrSrc2 += nbrLimbs;
            ptrDest += nbrLimbs;
          }
        }
      }
      ptrSrc1 += nbrLimbs;
    }
    return;
  }
  // Find the least power of 2 greater or equal than the maximum of factor1 and factor2.
  karatDegree = ((degree1 > degree2)? degree1 : degree2) + 1;
  if ((NumberLength == 1) && (karatDegree > 50) &&
    (karatDegree < (1000000 / TestNbr[0].x / TestNbr[0].x)))
  {
    fftPolyMult(factor1, factor2, polyMultTemp, degree1+1, degree2+1);
    return;
  }
  // Compute length of numbers for each recursion.
  if (karatDegree > KARATSUBA_POLY_CUTOFF)
  {
    int div = 1;
    while (karatDegree > KARATSUBA_POLY_CUTOFF)
    {
      div *= 2;
      karatDegree = (karatDegree + 1) / 2;
    }
    karatDegree *= div;
  }
  // Initialize Karatsuba polynomial.
  ptrValue1 = polyMultTemp;
  for (int currentDegree = 3 * karatDegree; currentDegree > 0; currentDegree--)
  {
    *ptrValue1 = 1;        // Initialize coefficient to zero.
    *(ptrValue1 + 1) = 0;
    ptrValue1 += nbrLimbs;
  }
  (void)memcpy(polyMultTemp, factor1, (degree1 + 1) * nbrLimbs * sizeof(limb));
  (void)memcpy(&polyMultTemp[karatDegree * nbrLimbs], factor2, (degree2 + 1) * nbrLimbs * sizeof(limb));
  KaratsubaPoly(0, karatDegree, nbrLimbs);
}

// Compute the polynomial polyInv as:
// polyInv = x^(2*polyDegree) / polymod
void GetPolyInvParm(int polyDegree, /*@in@*/int* polyMod)
{
  int degrees[15];
  int nbrDegrees = 0;
  int newtonDegree;
  int nbrLimbs = NumberLength + 1;

  newtonDegree = polyDegree;
  // Compute degrees to use in Newton loop.
  while (newtonDegree > 1)
  {
    degrees[nbrDegrees] = newtonDegree;
    nbrDegrees++;
    newtonDegree = (newtonDegree + 1) / 2;
  }

  // Point to leading coefficient of polyMod.
  // Initialize polyInv as x - K
  // where K is the coefficient of x^(polyDegree-1) of polyMod.
  SetNumberToOne(&polyInv[nbrLimbs]);
  LenAndLimbs2ArrLimbs(polyMod + (polyDegree - 1) * nbrLimbs, operand1.limbs, nbrLimbs);
  (void)memset(operand2.limbs, 0, nbrLimbs * sizeof(limb));
  SubtBigNbrMod(operand2.limbs, operand1.limbs, operand1.limbs);
  ArrLimbs2LenAndLimbs(polyInv, operand1.limbs, nbrLimbs);
  newtonDegree = 1;

  // Use Newton iteration: F_{2n}(x) = F_n(x)*(2-D(x)*F_n(x))
  // where F = polyInv (degree newtonDegree)
  // and D = polyMod (degree nextDegree).
  // Use poly5 as temporary polynomial 2-D(x)*F_n(x) (degree nextDegree).
  nbrDegrees--;
  while (nbrDegrees >= 0)
  {  
    int* ptrCoeff; 
    const int* ptrCoeff2;
    const int *ptrPolyMod;
    int nextDegree = degrees[nbrDegrees];
    // Initialize poly4 with the nextDegree most significant coefficients.
    ptrPolyMod = polyMod + (nbrLimbs * (polyDegree - nextDegree));
    (void)memcpy(poly4, ptrPolyMod, nextDegree * nbrLimbs * sizeof(limb));
    SetNumberToOne(&poly4[nextDegree * nbrLimbs]);
    polyInvCached = NBR_READY_TO_BE_CACHED;
    MultPolynomial(nextDegree, newtonDegree, poly4, polyInv);
    ptrCoeff = poly5;   // Destination of 2-D(x)*F_n(x).
    ptrCoeff2 = &polyMultTemp[newtonDegree * nbrLimbs];  // Source of D(x)*F_n(x).
    (void)memset(operand2.limbs, 0, nbrLimbs * sizeof(limb));
    for (int deg = 0; deg < nextDegree; deg++)
    {
      LenAndLimbs2ArrLimbs(ptrCoeff2, operand3.limbs, nbrLimbs);
      SubtBigNbrMod(operand2.limbs, operand3.limbs, operand3.limbs);
      ArrLimbs2LenAndLimbs(ptrCoeff, operand3.limbs, nbrLimbs);
      ptrCoeff += nbrLimbs;
      ptrCoeff2 += nbrLimbs;
    }
    SetNumberToOne(ptrCoeff);
    MultPolynomial(nextDegree, newtonDegree, poly5, polyInv);
    (void)memcpy(polyInv, &polyMultTemp[newtonDegree * nbrLimbs],
        (nextDegree+1) * nbrLimbs * sizeof(limb));
    newtonDegree = nextDegree;
    nbrDegrees--;
  }
  polyInvCached = NBR_READY_TO_BE_CACHED;
}

// Multiply two polynomials mod polyMod using inverse polynomial polyInv.
// The algorithm is:
// m <- (T*polyInv)/x^polyDegree
// return T - m*polyMod
void multUsingInvPolynomial(/*@in@*/int* polyFact1, /*@in@*/int* polyFact2,
  /*@out@*/int* polyProduct,
  int polyDegree, /*@in@*/int* polyMod)
{
  int currentDegree;
  int index;
  int nbrLimbs = NumberLength + 1;
  // Compute T
  MultPolynomial(polyDegree, polyDegree, polyFact1, polyFact2);
  (void)memcpy(polyMultT, polyMultTemp, ((2 * polyDegree) + 1) * nbrLimbs * sizeof(limb));
  // Compute m
  MultPolynomial(polyDegree, polyDegree, &polyMultT[polyDegree * nbrLimbs], polyInv);
  (void)memcpy(polyMultM, &polyMultTemp[polyDegree * nbrLimbs],
    (polyDegree + 1) * nbrLimbs * sizeof(limb));
  // Compute m*polyMod
  MultPolynomial(polyDegree, polyDegree, polyMultM, polyMod);
  // Compute T - mN.
  if (NumberLength == 1)
  {
    int mod = TestNbr[0].x;
    index = 1;
    for (currentDegree = 0; currentDegree < polyDegree; currentDegree++)
    {
      int temp = polyMultT[index] - polyMultTemp[index];
      temp += mod & (temp >> 31);
      *(polyProduct + 1) = temp;
      index += nbrLimbs;
      polyProduct += nbrLimbs;
    }
  }
  else
  {
    index = 0;
    for (currentDegree = 0; currentDegree < polyDegree; currentDegree++)
    {
      LenAndLimbs2ArrLimbs(&polyMultTemp[index], operand1.limbs, nbrLimbs);
      LenAndLimbs2ArrLimbs(&polyMultT[index], operand2.limbs, nbrLimbs);
      SubtBigNbrMod(operand2.limbs, operand1.limbs, operand1.limbs);
      ArrLimbs2LenAndLimbs(polyProduct, operand1.limbs, nbrLimbs);
      index += nbrLimbs;
      polyProduct += nbrLimbs;
    }
  }
}

// Multiply two polynomials mod polyMod.
void multPolynomialModPoly(const int* polyFact1, const int* polyFact2,
  /*@out@*/int* polyProduct, int polyDegree, const int* polyMod)
{
  int index1;
  int index2;
  int nbrLimbs = NumberLength + 1;
  const int* ptrPoly1;
  const int* ptrPoly2;
  int* ptrPolyTemp;
  // Initialize polyMultTemp with the most significant half of product.
  ptrPolyTemp = polyMultTemp + (polyDegree - 1) * nbrLimbs;
  for (index1 = polyDegree - 1; index1 >= 0; index1--)
  {
    ptrPoly1 = polyFact1 + (index1 * nbrLimbs);
    ptrPoly2 = polyFact2 + (polyDegree - 1) * nbrLimbs;
    IntArray2BigInteger(ptrPoly1, &operand1);
    IntArray2BigInteger(ptrPoly2, &operand2);
    modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    for (index2 = polyDegree - 2; index2 >= index1; index2--)
    {
      ptrPoly2 -= nbrLimbs;
      ptrPoly1 += nbrLimbs;
      IntArray2BigInteger(ptrPoly1, &operand2);
      IntArray2BigInteger(ptrPoly2, &operand3);
      modmult(operand2.limbs, operand3.limbs, operand2.limbs);
      AddBigNbrMod(operand2.limbs, operand1.limbs, operand1.limbs);
    }
    BigInteger2IntArray(ptrPolyTemp, &operand1);
    ptrPolyTemp -= nbrLimbs;
  }
  // Get remainder of long division by polyMod and append next limbs of the product.
  for (index1 = polyDegree - 2; index1 >= 0; index1--)
  {
    int* ptrPoly3 = &polyMultTemp[(polyDegree - 1) * nbrLimbs];
    // Back up leading coefficient.
    (void)memcpy(ptrPoly3 + nbrLimbs, ptrPoly3, nbrLimbs * sizeof(int));
    IntArray2BigInteger(ptrPoly3, &operand3);
    ptrPoly1 = polyMod + (polyDegree * nbrLimbs);
    for (index2 = polyDegree - 2; index2 >= 0; index2--)
    {
      ptrPoly1 -= nbrLimbs;
      ptrPoly3 -= nbrLimbs;
      IntArray2BigInteger(ptrPoly1, &operand1);
      modmult(operand3.limbs, operand1.limbs, operand1.limbs);
      IntArray2BigInteger(ptrPoly3, &operand2);
      SubtBigNbrMod(operand2.limbs, operand1.limbs, operand1.limbs);
      BigInteger2IntArray(ptrPoly3 + nbrLimbs, &operand1);
    }
    ptrPoly1 = polyFact1 + (index1 * nbrLimbs);
    ptrPoly2 = polyFact2;
    IntArray2BigInteger(ptrPoly1, &operand1);
    IntArray2BigInteger(ptrPoly2, &operand2);
    modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    for (index2 = 1; index2 <= index1; index2++)
    {
      ptrPoly1 -= nbrLimbs;
      ptrPoly2 += nbrLimbs;
      IntArray2BigInteger(ptrPoly1, &operand2);
      IntArray2BigInteger(ptrPoly2, &operand3);
      modmult(operand2.limbs, operand3.limbs, operand2.limbs);
      AddBigNbrMod(operand2.limbs, operand1.limbs, operand1.limbs);
    }
    IntArray2BigInteger(polyMod, &operand2);
    IntArray2BigInteger(&polyMultTemp[polyDegree * nbrLimbs], &operand3);
    modmult(operand3.limbs, operand2.limbs, operand2.limbs);
    SubtBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
    BigInteger2IntArray(polyMultTemp, &operand1);
  }
  (void)memcpy(polyProduct, polyMultTemp, polyDegree * nbrLimbs * sizeof(int));
}

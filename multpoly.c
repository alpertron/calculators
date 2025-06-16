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
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include "bignbr.h"
#include "expression.h"
#include "highlevel.h"
#include "polynomial.h"
#include "showtime.h"
#include "fft.h"

#define KARATSUBA_POLY_CUTOFF      16
#define SQRT_KARATSUBA_POLY_CUTOFF 4

struct stKaratsubaStack
{
  int idxFactor1;
  int stage;
};

#define ADJUST_MODULUS(sum, modulus) (sum += modulus & (sum >> 31))
static struct stKaratsubaStack astKaratsubaStack[10];
static BigInteger coeff[2 * KARATSUBA_POLY_CUTOFF];
int polyInvCached;
extern int powerOf2Exponent;
limb subst1Limbs[MAX_LEN_MULT];
limb subst2Limbs[MAX_LEN_MULT];
static limb tmp1[MAX_LEN];

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
  assert(idxFactor1 >= 0);
  assert(idxFactor2 >= 0);
  assert(nbrLimbs >= 1);
#ifndef _USING64BITS_
  limb result;
#endif
  for (i = 0; i < (2 * coeffLen) - 1; i++)
  {    // Process each limb of product (least to most significant limb).
    int lenBytes;
    if (i < coeffLen)
    {   // Processing first half (least significant) of product.
      ptrFactor2 = &common.poly.polyMultTemp[(idxFactor2 + i) * nbrLimbs];
      ptrFactor1 = &common.poly.polyMultTemp[idxFactor1 * nbrLimbs];
      j = i;
    }
    else
    {  // Processing second half (most significant) of product.
      ptrFactor2 = &common.poly.polyMultTemp[(idxFactor2 + coeffLen - 1) * nbrLimbs];
      ptrFactor1 = &common.poly.polyMultTemp[(idxFactor1 + i - coeffLen + 1) * nbrLimbs];
      j = 2 * (coeffLen - 1) - i;
    }
    lenBytes = nbrLimbs * (int)sizeof(limb);
    (void)memset(coeff[i].limbs, 0, lenBytes);
    if (nbrLimbs == 2)
    {    // Optimization for the case when there is only one limb.
      int modulus = TestNbr[0].x;
      int sum = 0;
#ifdef _USING64BITS_
      uint64_t ui64Sum;
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
        ui64Sum = 0;
#else
        double dSum = 0;
#endif
        for (; j >= 3; j -= 4)
        {
          int prod1 = *ptrFactor1 * *ptrFactor2;
          int prod2 = *(ptrFactor1 + 2) * *(ptrFactor2 - 2);
          int prod3 = *(ptrFactor1 + 4) * *(ptrFactor2 - 4);
          int prod4 = *(ptrFactor1 + 6) * *(ptrFactor2 - 6);
#ifdef _USING64BITS_
          ui64Sum += (uint64_t)prod1 + (uint64_t)prod2 +
            (uint64_t)prod3 + (uint64_t)prod4;
#else
          dSum += (double)prod1 + (double)prod2 +
            (double)prod3 + (double)prod4;
#endif
          ptrFactor1 += 8;
          ptrFactor2 -= 8;
        }
        while (j >= 0)
        {
          int prod1 = *ptrFactor1 * *ptrFactor2;
#ifdef _USING64BITS_
          ui64Sum += (uint64_t)prod1;
#else
          dSum += (double)prod1;
#endif
          ptrFactor1 += 2;
          ptrFactor2 -= 2;
          j--;
        }
#ifdef _USING64BITS_
        ui64Sum = ui64Sum % (uint64_t)modulus;
        sum = (int)ui64Sum;
#else
        dSum = dSum - floor(dSum / (double)modulus) * (double)modulus;
        sum = (int)dSum;
#endif
      }
      else
      {
#ifdef _USING64BITS_
        ui64Sum = 0;
#endif
        for (; j >= 0; j--)
        {
#ifdef _USING64BITS_
          ui64Sum += (uint64_t)*ptrFactor1 * (uint64_t)*ptrFactor2;
          if ((int64_t)ui64Sum < 0)
          {
            ui64Sum %= modulus;
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
        ui64Sum %= (unsigned int)modulus;
        sum = (int)ui64Sum;
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
  ptrFactor1 = &common.poly.polyMultTemp[idxFactor1 * nbrLimbs];
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

#define checkLimbForNbrLimbs2(nbr)  \
   case nbr:             \
   if (*ptrResult != 0)  \
   {                     \
     return false;       \
   }                     \
   ptrResult += nbrLimbs;

#define checkLimbForNbrLimbsNot2(nbr)  \
   case nbr:             \
   if (((*(ptrResult - 1) - 1) | *ptrResult) != 0)  \
   {                     \
     return false;       \
   }                     \
   ptrResult += nbrLimbs;

static inline bool isFactorEqualToZero(const int* pResult, int nbrLimbs, int nbrLength)
{
  const int* ptrResult = pResult;
  if (nbrLimbs == 2)
  {
    switch (nbrLength)
    {
      checkLimbForNbrLimbs2(16)
      checkLimbForNbrLimbs2(15)
      checkLimbForNbrLimbs2(14)
      checkLimbForNbrLimbs2(13)
      checkLimbForNbrLimbs2(12)
      checkLimbForNbrLimbs2(11)
      checkLimbForNbrLimbs2(10)
      checkLimbForNbrLimbs2(9)
      checkLimbForNbrLimbs2(8)
      checkLimbForNbrLimbs2(7)
      checkLimbForNbrLimbs2(6)
      checkLimbForNbrLimbs2(5)
      checkLimbForNbrLimbs2(4)
      checkLimbForNbrLimbs2(3)
      checkLimbForNbrLimbs2(2)
        /* No break */
    default:     // length equals to 1.
      if (*ptrResult != 0)
      {
        return false;
      }
      break;
    }
  }
  else
  {
    switch (nbrLength)
    {
      checkLimbForNbrLimbsNot2(16)
      checkLimbForNbrLimbsNot2(15)
      checkLimbForNbrLimbsNot2(14)
      checkLimbForNbrLimbsNot2(13)
      checkLimbForNbrLimbsNot2(12)
      checkLimbForNbrLimbsNot2(11)
      checkLimbForNbrLimbsNot2(10)
      checkLimbForNbrLimbsNot2(9)
      checkLimbForNbrLimbsNot2(8)
      checkLimbForNbrLimbsNot2(7)
      checkLimbForNbrLimbsNot2(6)
      checkLimbForNbrLimbsNot2(5)
      checkLimbForNbrLimbsNot2(4)
      checkLimbForNbrLimbsNot2(3)
      checkLimbForNbrLimbsNot2(2)
        /* No break */
    default:     // length equals to 1.
      if (((*(ptrResult - 1) - 1) | *ptrResult) != 0)
      {
        return false;
      }
      break;
    }
  }
  return true;
}

// Recursive Karatsuba function.
static void KaratsubaPoly(int idxFact1, int nbrLen, int nbrLimbs)
{
  int i;
  int idxFactor1 = idxFact1;
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
  assert(nbrLimbs >= 1);
  assert(idxFact1 >= 0);
  assert(nbrLen >= 1);
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
      assert(idxFactor2 >= 0);
      if (nbrLength <= KARATSUBA_POLY_CUTOFF)
      {
        // Check if one of the factors is equal to zero.
        ptrResult = &common.poly.polyMultTemp[(idxFactor1 * nbrLimbs) + 1];
        if (isFactorEqualToZero(ptrResult, nbrLimbs, nbrLength))
        {      // First factor is zero. Initialize second to zero.
          ptrResult = &common.poly.polyMultTemp[(idxFactor2 * nbrLimbs) + 1];
          if (nbrLimbs == 2)
          {
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
              *(ptrResult - 1) = 1;
              *ptrResult = 0;
              ptrResult += nbrLimbs;
            }
          }
        }
        else
        {     // First factor is not zero. Check second.
          ptrResult = &common.poly.polyMultTemp[(idxFactor2 * nbrLimbs) + 1];
          assert(idxFactor2 >= 0);
          if (isFactorEqualToZero(ptrResult, nbrLimbs, nbrLength))
          {    // Second factor is zero. Initialize first to zero.
            ptrResult = &common.poly.polyMultTemp[(idxFactor1 * nbrLimbs) + 1];
            if (nbrLimbs == 2)
            {
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
                *(ptrResult - 1) = 1;
                *ptrResult = 0;
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
        assert(pstKaratsubaStack >= &astKaratsubaStack[0]);
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
      int* ptrHighFirstFactor = &common.poly.polyMultTemp[(idxFactor1 + halfLength) * nbrLimbs];
      int* ptrLowSecondFactor = &common.poly.polyMultTemp[idxFactor2 * nbrLimbs];
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
        int sizeCoeffInBytes = nbrLimbs * (int)sizeof(int);
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
      ptr1 = &common.poly.polyMultTemp[idxFactor1 * nbrLimbs];
      ptr2 = &common.poly.polyMultTemp[idxFactor2 * nbrLimbs];
      ptrResult = &common.poly.polyMultTemp[diffIndex * nbrLimbs];
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
      halfLength = nbrLength / 2;
      // Obtain (b+1)(xH*yH*b + xL*yL) = xH*yH*b^2 + (xL*yL+xH*yH)*b + xL*yL
      // The first and last terms are already in correct locations.
      // Add (xL*yL+xH*yH)*b.
      ptrResult = &common.poly.polyMultTemp[(idxFactor1 + halfLength) * nbrLimbs];
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
      ptrHigh = &common.poly.polyMultTemp[diffIndex * nbrLimbs];
      assert(diffIndex >= 0);
      ptrResult = &common.poly.polyMultTemp[(idxFactor1 + halfLength) * nbrLimbs];
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
      assert(pstKaratsubaStack >= &astKaratsubaStack[0]);
      idxFactor1 = pstKaratsubaStack->idxFactor1;
      stage = pstKaratsubaStack->stage;
      break;
    }     // End switch
  } while (stage >= 0);
}

// Multiply factor1 by factor2.The result will be stored in common.poly.polyMultTemp.
static void MultIntegerPolynomial(int deg1, int deg2,
  const int* fact1, const int* fact2)
{
  int degree1;
  int degree2;
  const int* factor1;
  const int* factor2;
  int indexes[2][MAX_DEGREE + 1];
  int* ptrIndex;
  int* piDest;
  int currentDegree;
  int index;
  int degreeF2;
  assert(deg1 >= 0);
  assert(deg2 >= 0);
  // Force degree1 >= degree2.
  if (deg1 < deg2)
  {
    degree1 = deg2;
    degree2 = deg1;
    factor1 = fact2;
    factor2 = fact1;
  }
  else
  {
    degree1 = deg1;
    degree2 = deg2;
    factor1 = fact1;
    factor2 = fact2;
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
  piDest = common.poly.polyMultTemp;
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
    piDest += numLimbs(piDest);
    piDest++;
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
    piDest += numLimbs(piDest);
    piDest++;
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
    piDest += numLimbs(piDest);
    piDest++;
  }
}

// Compact polynomial poly so each coefficient has bit size lenBitsDest
// in BigInteger substed.
static void KroneckerSubstitution(const int* poly, int polyDegree,
  int bitsDegree, int bitsModulus, limb* substedLimbs, int *substedNbrLimbs)
{
  const int* ptrSrc = poly;
  int nbrLimbs = NumberLength + 1;
  int lenBitsDest = (2 * bitsModulus) + bitsDegree;
  int currentBitOffset = 0;
  for (int currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {
    unsigned int shLeft = (unsigned int)currentBitOffset % 
      (unsigned int)BITS_PER_GROUP;
    unsigned int shRight = (unsigned int)BITS_PER_GROUP - shLeft;
    limb *ptrDest = substedLimbs + (currentBitOffset / BITS_PER_GROUP);
    currentBitOffset += lenBitsDest;
    const limb* ptrDestNew = substedLimbs + (currentBitOffset / BITS_PER_GROUP);
    unsigned int carry;
    const int* ptrSrcBak;
    if (currentBitOffset == lenBitsDest)
    {
      carry = 0U;
    }
    else
    {
      carry = (unsigned int)ptrDest->x;
    }    
    ptrSrcBak = ptrSrc;
    ptrSrc++;
    for (int lenLimbsSrc = *(ptrSrc-1); lenLimbsSrc > 0; lenLimbsSrc--)
    {     // Shift left currentBitOffset bits.
      ptrDest->x = UintToInt((((unsigned int)*ptrSrc << shLeft) + carry) & MAX_VALUE_LIMB);
      carry = UintToInt((unsigned int)*ptrSrc >> shRight);
      ptrDest++;
      ptrSrc++;      
    }
    ptrDest->x = carry;
    while (ptrDest < ptrDestNew)
    {
      ptrDest++;
      ptrDest->x = 0;
    }
    ptrSrc = ptrSrcBak + nbrLimbs;
  }
  nbrLimbs = (currentBitOffset+BITS_PER_GROUP_MINUS_1) / BITS_PER_GROUP;
  while (nbrLimbs > 1)
  {
    if ((substedLimbs + nbrLimbs - 1)->x != 0)
    {
      break;
    }
    nbrLimbs--;
  }
  *substedNbrLimbs = nbrLimbs;
}

// Kronecker substitution replaces polynomial multiplication by
// big integer multiplication by inserting zeros between coefficients
// so there is no overflow when performing the multiplication.
// Then the product is converted back to a polynomial.
static bool MultiplyUsingKroneckerSubst(int degree1, int degree2,
  const int* factor1, const int* factor2)
{
  int degreeProd = degree1 + degree2;
  int* ptrDest;
  int lenBitsCoeff;
  int lenBitsCoeffToProcess;
  unsigned int mostSignificantLimbMask;
  int currentBitOffset = 0;
  int nbrLimbs;
  int nbrLimbsProduct;
  int nbrBitsProduct;
  unsigned int bitsDegree;
  unsigned int bitsModulus;
  int mostSignificantLimb;
  int len1Limbs;
  int len2Limbs;
  int lenProdLimbs;
  limb* prod;

  // Get number of bits of degree.
  for (bitsDegree = 0U; bitsDegree < 16U; bitsDegree++)
  {
    if ((1U << bitsDegree) > (unsigned int)degreeProd)
    {
      break;
    }
  }
  // Get number of bits of modulus.
  bitsModulus = NumberLength * BITS_PER_GROUP;
  mostSignificantLimb = TestNbr[NumberLength-1].x;
  for (int mask = HALF_INT_RANGE; mask > 0; mask /= 2)
  {
    if ((mask & mostSignificantLimb) != 0)
    {
      break;
    }
    bitsModulus--;
  }
  // With this choice of bit size, there will be no overflow
  // of coefficients when performing the multiplication.
  lenBitsCoeff = (int)((2U * bitsModulus) + bitsDegree);
  nbrBitsProduct = lenBitsCoeff * (degreeProd + 1);
  nbrLimbsProduct = (nbrBitsProduct + BITS_PER_GROUP_MINUS_1) / BITS_PER_GROUP;
  if (nbrLimbsProduct >= MAX_LEN_MULT)
  {       // Number would be very large. Cannot perform substitution.
    return false;
  }
  // Compact each polynomial factor so each coefficient has bit size lenBitsDest.
  KroneckerSubstitution(factor1, degree1, bitsDegree, bitsModulus,
    subst1Limbs, &len1Limbs);
  KroneckerSubstitution(factor2, degree2, bitsDegree, bitsModulus,
    subst2Limbs, &len2Limbs);
  multiplyWithBothLen(subst1Limbs, subst2Limbs, subst2Limbs,
    len1Limbs, len2Limbs, &lenProdLimbs);
  while (lenProdLimbs < nbrLimbsProduct)
  {  // Loop that generates most significant zeros.
    (subst2Limbs + lenProdLimbs)->x = 0;
    lenProdLimbs++;
  }
  // At this moment subst2Limbs has the product of the two polynomials
  // where each coefficient has bit size lenBitsDest.
  // Get each coefficient and get their remainders modulo TestNbr.
  // These numbers will be the coefficients of the polynomial multiplication.
  operand2.nbrLimbs = NumberLength;
  ptrDest = common.poly.polyMultTemp;
  if (powerOf2Exponent != 0)
  {          // TestNbr is even.
    lenBitsCoeffToProcess = powerOf2Exponent;
  }
  else
  {
    lenBitsCoeffToProcess = lenBitsCoeff;
  }
  mostSignificantLimbMask = lenBitsCoeffToProcess % BITS_PER_GROUP;
  if (mostSignificantLimbMask == 0U)
  {
    mostSignificantLimbMask = (unsigned int)BITS_PER_GROUP;
  }
  mostSignificantLimbMask = (1U << mostSignificantLimbMask) - 1U;
  for (int currentDegree = 0; currentDegree <= degreeProd; currentDegree++)
  {
    prod = (limb *)ptrDest + 1;
    limb *ptrDividend = prod;
    // Shift right coefficient to generate coefficient in prod.
    const limb* ptrSrc = subst2Limbs + (currentBitOffset / BITS_PER_GROUP);
    unsigned int shRight = (unsigned int)currentBitOffset %
      (unsigned int)BITS_PER_GROUP;
    unsigned int shLeft = (unsigned int)BITS_PER_GROUP - shRight;
    currentBitOffset += lenBitsCoeff;
    for (int currentBitNbr = 0; currentBitNbr < lenBitsCoeffToProcess;
      currentBitNbr += BITS_PER_GROUP)
    {
      ptrDividend->x = UintToInt((((unsigned int)ptrSrc->x >> shRight) |
        ((unsigned int)(ptrSrc + 1)->x << shLeft)) & MAX_VALUE_LIMB);
      ptrSrc++;
      ptrDividend++;
    }
    (ptrDividend - 1)->x &= mostSignificantLimbMask;
    if (powerOf2Exponent == 0)
    {
      ptrDividend->x = 0;
      (ptrDividend + 1)->x = 0;
      // Montgomery reduction works if the number has up to
      // bitsModulus + NumberLength * BITS_PER_GROUP bits,
      // but the length is 2*bitsModulus + bitsDegree.
      // So we have to execute AdjustModN before the Montgomery reduction.
      if (lenBitsCoeff >= ((int)bitsModulus + (NumberLength * BITS_PER_GROUP)))
      {
        AdjustModN(&prod[NumberLength], TestNbr, NumberLength);
      }
      if (NumberLength == 1)
      {
        if (TestNbr[0].x <= 32768)
        {
          prod->x = prod->x % TestNbr[0].x;
        }
        else
        {
#ifdef _USING64BITS_
          prod->x = (prod->x +
            ((int64_t)(prod+1)->x << BITS_PER_GROUP)) % TestNbr[0].x;
#else
          // Round up quotient.
          int quotient = (int)floor(((double)prod->x +
            ((double)(prod+1)->x * (double)MAX_VALUE_LIMB)) / (double)TestNbr[0].x + 0.5);
          int remainder = (prod->x +
            ((prod+1)->x << BITS_PER_GROUP)) - (quotient * TestNbr[0].x);
          if (remainder < 0)
          {    // Quotient was 1 more than expected. Adjust remainder.
            remainder += TestNbr[0].x;
          }
          prod->x = remainder;
#endif
        }
      }
      else
      {
        // Perform Montgomery reduction.
        // Compute m
        multiply(prod, MontgomeryMultN, tmp1, NumberLength, NULL);
        // Compute mN
        multiply(tmp1, TestNbr, tmp1, NumberLength, NULL);
        endBigModmult(tmp1, prod, prod);
      }
    }
    // Find number of limbs of coefficients.
    for (nbrLimbs = NumberLength; nbrLimbs > 1; nbrLimbs--)
    {
      if ((prod + nbrLimbs - 1)->x != 0)
      {
        break;
      }
    }
    *ptrDest = nbrLimbs;
    nbrLimbs = NumberLength + 1;
    ptrDest += nbrLimbs;
  }
  return true;
}

// Multiply factor1 by factor2. The result will be stored in common.poly.polyMultTemp.
void MultPolynomial(int degree1, int degree2,
  const int* factor1, const int* factor2)
{
  int* ptrValue1;
  int nbrLimbs;
  int karatDegree;
  int lenBytes;
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
    int* ptrDest = common.poly.polyMultTemp;
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
        ptrDest = &common.poly.polyMultTemp[currentDegree1 * nbrLimbs];
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
    fftPolyMult(factor1, factor2, common.poly.polyMultTemp, degree1+1, degree2+1);
    return;
  }
  if ((NumberLength > 1) && (degree1 > 5) && (degree2 > 5) &&
    MultiplyUsingKroneckerSubst(degree1, degree2, factor1, factor2))
  {             // Multiply polynomials using Kronecker substitution.
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
  ptrValue1 = common.poly.polyMultTemp;
  for (int currentDegree = 3 * karatDegree; currentDegree > 0; currentDegree--)
  {
    *ptrValue1 = 1;        // Initialize coefficient to zero.
    *(ptrValue1 + 1) = 0;
    ptrValue1 += nbrLimbs;
  }
  lenBytes = (degree1 + 1) * nbrLimbs * (int)sizeof(limb);
  (void)memcpy(common.poly.polyMultTemp, factor1, lenBytes);
  lenBytes = (degree2 + 1) * nbrLimbs * (int)sizeof(limb);
  (void)memcpy(&common.poly.polyMultTemp[karatDegree * nbrLimbs], factor2, lenBytes);
  KaratsubaPoly(0, karatDegree, nbrLimbs);
}

// Compute the polynomial common.poly.polyInv as:
// common.poly.polyInv = x^(2*polyDegree) / polymod
void GetPolyInvParm(int polyDegree, const int* polyMod)
{
  int degrees[15];
  int nbrDegrees = 0;
  int newtonDegree;
  int nbrLimbs = NumberLength + 1;
  int lenBytes;

  newtonDegree = polyDegree;
  // Compute degrees to use in Newton loop.
  while (newtonDegree > 1)
  {
    degrees[nbrDegrees] = newtonDegree;
    nbrDegrees++;
    newtonDegree = (newtonDegree + 1) / 2;
  }

  // Point to leading coefficient of polyMod.
  // Initialize common.poly.polyInv as x - K
  // where K is the coefficient of x^(polyDegree-1) of polyMod.
  SetNumberToOne(&common.poly.polyInv[nbrLimbs]);
  LenAndLimbs2ArrLimbs(polyMod + (polyDegree - 1) * nbrLimbs, operand1.limbs, nbrLimbs);
  lenBytes = nbrLimbs * (int)sizeof(limb);
  (void)memset(operand2.limbs, 0, lenBytes);
  SubtBigNbrMod(operand2.limbs, operand1.limbs, operand1.limbs);
  ArrLimbs2LenAndLimbs(common.poly.polyInv, operand1.limbs, nbrLimbs);
  newtonDegree = 1;

  // Use Newton iteration: F_{2n}(x) = F_n(x)*(2-D(x)*F_n(x))
  // where F = common.poly.polyInv (degree newtonDegree)
  // and D = polyMod (degree nextDegree).
  // Use common.poly.poly5 as temporary polynomial 2-D(x)*F_n(x) (degree nextDegree).
  nbrDegrees--;
  while (nbrDegrees >= 0)
  {  
    int* ptrCoeff; 
    const int* ptrCoeff2;
    const int *ptrPolyMod;
    int nextDegree = degrees[nbrDegrees];
    // Initialize common.poly.poly4 with the nextDegree most significant coefficients.
    ptrPolyMod = polyMod + (nbrLimbs * (polyDegree - nextDegree));
    lenBytes = nextDegree * nbrLimbs * (int)sizeof(limb);
    (void)memcpy(common.poly.poly4, ptrPolyMod, lenBytes);
    SetNumberToOne(&common.poly.poly4[nextDegree * nbrLimbs]);
    polyInvCached = NBR_READY_TO_BE_CACHED;
    MultPolynomial(nextDegree, newtonDegree, common.poly.poly4, common.poly.polyInv);
    ptrCoeff = common.poly.poly5;   // Destination of 2-D(x)*F_n(x).
    ptrCoeff2 = &common.poly.polyMultTemp[newtonDegree * nbrLimbs];  // Source of D(x)*F_n(x).
    int nbrLimbsBytes = nbrLimbs * sizeof(limb);
    (void)memset(operand2.limbs, 0, nbrLimbsBytes);
    for (int deg = 0; deg < nextDegree; deg++)
    {
      LenAndLimbs2ArrLimbs(ptrCoeff2, operand3.limbs, nbrLimbs);
      SubtBigNbrMod(operand2.limbs, operand3.limbs, operand3.limbs);
      ArrLimbs2LenAndLimbs(ptrCoeff, operand3.limbs, nbrLimbs);
      ptrCoeff += nbrLimbs;
      ptrCoeff2 += nbrLimbs;
    }
    SetNumberToOne(ptrCoeff);
    MultPolynomial(nextDegree, newtonDegree, common.poly.poly5, common.poly.polyInv);
    lenBytes = (nextDegree + 1) * nbrLimbs * (int)sizeof(limb);
    (void)memcpy(common.poly.polyInv, &common.poly.polyMultTemp[newtonDegree * nbrLimbs], lenBytes);
    newtonDegree = nextDegree;
    nbrDegrees--;
  }
  polyInvCached = NBR_READY_TO_BE_CACHED;
}

// Multiply two polynomials mod polyMod using inverse polynomial common.poly.polyInv.
// The algorithm is:
// T <- polyFact1 * polyFact2
// m <- (T*common.poly.polyInv)/x^polyDegree
// return T - m*polyMod
void multUsingInvPolynomial(const int* polyFact1, const int* polyFact2,
  /*@out@*/int* polyProd,
  int polyDegree, const int* polyMod)
{
  int* polyProduct = polyProd;
  int currentDegree;
  int index;
  int nbrLimbs = NumberLength + 1;
  int lenBytes;
  // Compute T
  MultPolynomial(polyDegree, polyDegree, polyFact1, polyFact2);
  lenBytes = ((2 * polyDegree) + 1) * nbrLimbs * (int)sizeof(limb);
  (void)memcpy(common.poly.polyMultT, common.poly.polyMultTemp, lenBytes);
  // Compute m
  MultPolynomial(polyDegree, polyDegree, &common.poly.polyMultT[polyDegree * nbrLimbs], common.poly.polyInv);
  lenBytes = (polyDegree + 1) * nbrLimbs * (int)sizeof(limb);
  (void)memcpy(common.poly.polyMultM, &common.poly.polyMultTemp[polyDegree * nbrLimbs], lenBytes);
  // Compute m*polyMod
  MultPolynomial(polyDegree, polyDegree, common.poly.polyMultM, polyMod);
  // Compute T - mN.
  if (NumberLength == 1)
  {
    int mod = TestNbr[0].x;
    index = 1;
    for (currentDegree = 0; currentDegree < polyDegree; currentDegree++)
    {
      int temp = common.poly.polyMultT[index] - common.poly.polyMultTemp[index];
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
      LenAndLimbs2ArrLimbs(&common.poly.polyMultTemp[index], operand1.limbs, nbrLimbs);
      LenAndLimbs2ArrLimbs(&common.poly.polyMultT[index], operand2.limbs, nbrLimbs);
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
  int lenBytes;
  // Initialize common.poly.polyMultTemp with the most significant half of product.
  ptrPolyTemp = common.poly.polyMultTemp + (polyDegree - 1) * nbrLimbs;
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
    int* ptrPoly3 = &common.poly.polyMultTemp[(polyDegree - 1) * nbrLimbs];
    // Back up leading coefficient.
    lenBytes = nbrLimbs * (int)sizeof(int);
    (void)memcpy(ptrPoly3 + nbrLimbs, ptrPoly3, lenBytes);
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
    IntArray2BigInteger(&common.poly.polyMultTemp[polyDegree * nbrLimbs], &operand3);
    modmult(operand3.limbs, operand2.limbs, operand2.limbs);
    SubtBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
    BigInteger2IntArray(common.poly.polyMultTemp, &operand1);
  }
  lenBytes = polyDegree * nbrLimbs * (int)sizeof(int);
  (void)memcpy(polyProduct, common.poly.polyMultTemp, lenBytes);
}

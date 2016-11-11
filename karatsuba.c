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

#define KARATSUBA_CUTOFF 8
static limb arr[MAX_LEN];
static limb arrAux[MAX_LEN];
static int length;
static int karatLength;
static void Karatsuba(int idxFactor1, int length, int endIndex);
int multCtr, karatCtr;

void multiply(limb *factor1, limb *factor2, limb *result, int len, int *pResultLen)
{
	length = len;
    // Compute length of numbers for each recursion.
	if (length <= KARATSUBA_CUTOFF)
	{
		karatLength = length;
	}
	else
	{
    int div;

		div = 1;
		while (length > KARATSUBA_CUTOFF)
		{
			div *= 2;
			length = (length + 1) / 2;
		}
		length *= div;
		karatLength = length;
	}
  memset(arr, 0, 2 * length*sizeof(limb));
	memcpy(&arr[0], factor1, len*sizeof(limb));
	memcpy(&arr[length], factor2, len*sizeof(limb));
	Karatsuba(0, length, 2 * length);
  memcpy(result, &arr[2 * (karatLength - length)], 2 * length * sizeof(limb));
  if (pResultLen != NULL)
  {
    memcpy(result, &arr[2 * (karatLength - length)], 2 * length * sizeof(limb));
    if (arr[2 * (karatLength - length)-1].x == 0)
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
	                     int idxResult, int length)
{
	int sign = 0;
	limb carry;
	int i;
	for (i = length-1; i>=0; i--)
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
  for (i = 0; i < length; i++)
  {
    carry.x += arr[idxMinuend + i].x - arr[idxSubtrahend + i].x;
    arr[idxResult + i].x = carry.x & MAX_VALUE_LIMB;
    carry.x >>= BITS_PER_GROUP;
  }
	return sign;
}

static void ClassicalMult(int idxFactor1, int idxFactor2, int length)
{
  int i, j;
  limb carry, sum, hiProd;
  limb *ptrFactor1, *ptrFactor2;
  limb prod;
  multCtr++;
  sum.x = 0;  // Initialize sum of product of limbs in the same column.
  for (i = 0; i < 2 * length - 1; i++)
  {    // Process each limb of product (least to most significant limb).
    carry.x = 0;
    if (i < length)
    {   // Processing first half (least significant) of product.
      ptrFactor2 = &arr[idxFactor2 + i];
      ptrFactor1 = &arr[idxFactor1];
      j = i;
    }
    else
    {  // Processing second half (most significant) of product.
      ptrFactor2 = &arr[idxFactor2 + length - 1];
      ptrFactor1 = &arr[idxFactor1 + i - length + 1];
      j = 2 * (length - 1) - i;
    }
    for (; j>0; j -= 2)
    {
      prod.x = ptrFactor1->x * ptrFactor2->x +
               (ptrFactor1 + 1)->x * (ptrFactor2 - 1)->x;
      ptrFactor1 += 2;
      ptrFactor2 -= 2;
      hiProd.x = prod.x >> BITS_PER_GROUP;
      sum.x += prod.x & MAX_VALUE_LIMB;
      carry.x += hiProd.x;
    }
    if (j == 0)
    {   // Odd number of limbs: process last one.
      prod.x = ptrFactor1->x * ptrFactor2->x;
      hiProd.x = prod.x >> BITS_PER_GROUP;
      sum.x += prod.x & MAX_VALUE_LIMB;
      carry.x += hiProd.x;
    }
    arrAux[i].x = sum.x & MAX_VALUE_LIMB;
    sum.x = carry.x + (sum.x >> BITS_PER_GROUP);
  }
  arrAux[i].x = sum.x;
  memcpy(&arr[idxFactor1], &arrAux[0], 2 * length * sizeof(limb));
  return;
}

// Recursive Karatsuba function.
static void Karatsuba(int idxFactor1, int length, int endIndex)
{
	int idxFactor2 = idxFactor1 + length;
  int i;
  limb tmp, cy;
  limb *ptrResult, *ptrHigh;
	int middle;
	int sign;
	int halfLength;
	if (length <= KARATSUBA_CUTOFF)
	{
		// Check if one of the factors is equal to zero.
    ptrResult = &arr[idxFactor1];
		for (i = length; i > 0; i--)
		{
			if ((ptrResult++)->x != 0)
			{
				break;
			}
		}
		if (i > 0)
		{     // First factor is not zero. Check second.
      ptrResult = &arr[idxFactor2];
      for (i = length; i > 0; i--)
			{
        if ((ptrResult++)->x != 0)
				{
					break;
				}
			}
		}
		if (i==0)
		{    // One of the factors is equal to zero.
			for (i = length - 1; i >= 0; i--)
			{
				arr[idxFactor1 + i].x = arr[idxFactor2 + i].x = 0;
			}
			return;
		}
         // Below cutoff: perform standard classical multiplcation.
    ClassicalMult(idxFactor1, idxFactor2, length);
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
  karatCtr++;
  halfLength = length >> 1;
	for (i = idxFactor1 + halfLength; i<idxFactor2; i++)
	{
		tmp.x = arr[i].x;
		arr[i].x = arr[i + halfLength].x;
		arr[i + halfLength].x = tmp.x;
	}
	// At this moment the order is: xL, yL, xH, yH.
	// Get absolute values of (xH-xL) and (yL-yH) and the signs.
	sign = absSubtract(idxFactor1, idxFactor2, endIndex, halfLength);
	sign ^= absSubtract(idxFactor2 + halfLength, idxFactor1 + halfLength,
		endIndex + halfLength, halfLength);
	middle = endIndex;
	endIndex += length;
	Karatsuba(idxFactor1, halfLength, endIndex); // Multiply both low parts.
	Karatsuba(idxFactor2, halfLength, endIndex); // Multiply both high parts.
	Karatsuba(middle, halfLength, endIndex);     // Multiply the differences.
     // Process all carries at the end.
     // Obtain (b+1)(xH*yH*b + xL*yL) = xH*yH*b^2 + (xL*yL+xH*yH)*b + xL*yL
     // The first and last terms are already in correct locations.
  ptrResult = &arr[idxFactor1+halfLength];
  for (i = halfLength; i > 0; i--)
  {
    tmp.x = ptrResult->x;
    ptrResult->x += (ptrResult - halfLength)->x + (ptrResult + halfLength)->x;
    (ptrResult + halfLength)->x += tmp.x + (ptrResult + length)->x;
    ptrResult++;
  }
    // Compute final product.
  ptrHigh = &arr[middle];
  ptrResult = &arr[idxFactor1 + halfLength];
  cy.x = 0;
  if (sign != 0)
  {            // (xH-xL) * (yL-yH) is negative.
    for (i = length; i > 0; i--)
    {
      cy.x += ptrResult->x - (ptrHigh++)->x;
      (ptrResult++)->x = cy.x & MAX_VALUE_LIMB;
      cy.x >>= BITS_PER_GROUP;
    }
  }
  else
  {            // (xH-xL) * (yL-yH) is positive or zero.
    for (i = length; i > 0; i--)
    {
      cy.x += ptrResult->x + (ptrHigh++)->x;
      (ptrResult++)->x = cy.x & MAX_VALUE_LIMB;
      cy.x >>= BITS_PER_GROUP;
    }
  }
  for (i = halfLength; i > 0; i--)
  {
    cy.x += ptrResult->x;
    (ptrResult++)->x = cy.x & MAX_VALUE_LIMB;
    cy.x >>= BITS_PER_GROUP;
  }
}

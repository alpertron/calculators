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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bignbr.h"
#include <math.h>

static limb approxInvSqrt[MAX_LEN];
limb approxInv[MAX_LEN];
limb adjustedArgument[MAX_LEN];
limb arrAux[MAX_LEN];
int bitLengthCycle[20];

// This routine uses Newton iteration: if x is an approximate inverse square root of N,
// a better approximation is: x(3-Nxx)/2. After the inverse square root is computed,
// the square root is found just by multiplying by N.
// The argument is multiplied by a power of 4 so the most significant limb is
// between LIMB_RANGE/4 and LIMB_RANGE-1 and there is an even number of limbs.
// At the end of the calculation, the result is divided by the power of 2.

// All computations are done in little-endian notation.
// Find power of 4 that divides the number.
// output: pNbrLimbs = pointer to number of limbs
//         pPower4 = pointer to power of 4.
static void MultiplyBigNbrByMinPowerOf4(/*@in@*/ int *pPower4, /*@in@*/limb *number, int len, /*@out@*/limb *dest)
{
  int power4;
  int power2 = 0;
  limb mostSignficLimb;
  int index2, mask, shLeft;
  limb prevLimb, currLimb;
  limb *ptrDest, *ptrSrc;

  ptrSrc = number;
  if ((len & 1) != 0)
  {
    power2 = BITS_PER_GROUP;
  }
  else
  {
    power2 = 0;
  }
  shLeft = 0;
  mostSignficLimb.x = (number + len - 1)->x;
  for (mask = LIMB_RANGE/2; mask > 0; mask >>= 1)
  {
    if ((mostSignficLimb.x & mask) != 0)
    {
      break;
    }
    shLeft++;
  }
  power2 = (power2 + shLeft) & (-2);
  power4 = power2 >> 1;
  ptrDest = dest;
  if (power2 > BITS_PER_GROUP)
  {
    shLeft = power2 - BITS_PER_GROUP;
    (ptrDest++)->x = 0;
  }
  else
  {
    shLeft = power2;
  }
  // Multiply number by this power.
  prevLimb.x = 0;
  for (index2 = len; index2 >= 0; index2--)
  {
    currLimb.x = (ptrSrc++)->x;
    (ptrDest++)->x = ((currLimb.x << shLeft) | (prevLimb.x >> (BITS_PER_GROUP - shLeft))) & MAX_VALUE_LIMB;
    prevLimb.x = currLimb.x;
  }
  *pPower4 = power4;
}

void squareRoot(/*@in@*/limb *argument, /*@out@*/limb *sqRoot, int len, /*@out@*/int *pLenSqRoot)
{
  int index, lenInvSqrt, lenInvSqrt2;
  int bitLength;
  int bitLengthNbrCycles;
  int idx;
  limb prev;
  limb *ptrApproxInvSqrt;
  limb *ptrArrAux;
  int shRight;
  limb *ptrDest, *ptrSrc;
  double invSqrt;
  int currLimb;

  // If the number of limbs is even and the upper half of the number has all
  // limbs set to MAX_VALUE_LIMB, there could be overflow. Set the square root directly.
  if ((len & 1) == 0)
  {      // Even number of limbs.
    for (index = len / 2; index < len; index++)
    {
      if ((argument + index)->x != MAX_VALUE_LIMB)
      {
        break;
      }
    }
    if (index == len)
    {   // Set square root and go out.
      for (index = len / 2-1; index >= 0; index--)
      {
        (sqRoot + index)->x = MAX_VALUE_LIMB;
      }
      *pLenSqRoot = len / 2;
      return;
    }
  }
  // Obtain logarithm of 2 of argument.
  for (index=len - 1; index>2; index--)
  {
    if ((argument+index)->x != 0)
    {
      break;
    }
  }
  memset(sqRoot, 0, (len+1)/2*sizeof(limb));
  if (index <= 1)
  {                // Argument is small, so compute directly its square root.
    if (index == 0)
    {
      sqRoot->x = (int)floor(sqrt(argument->x)+0.000001);
    }
    else
    {
      int square[3];   // MultBigNbr routine uses an extra limb for result.
      double dArg = argument->x + (double)(argument + 1)->x*(double)LIMB_RANGE;
      dArg = floor(sqrt(dArg + 0.5));
      if (dArg == (double)LIMB_RANGE)
      {
        dArg = (double)MAX_VALUE_LIMB;
      }
      sqRoot->x = (int)dArg;
      MultBigNbr(&sqRoot->x, &sqRoot->x, square, 1);
      if (square[1] > (argument + 1)->x ||
        (square[1] == (argument + 1)->x && square[0] > argument->x))
      {
        sqRoot->x--;
      }
    }
    *pLenSqRoot = 1;    // Only one limb for result.
    return;
  }
  // Adjust argument.
  MultiplyBigNbrByMinPowerOf4(&shRight, argument, len, adjustedArgument);
  if ((len & 1) != 0)
  {
    len++;   // Make number of limbs even.
  }
  lenInvSqrt = (len + 5)>> 1;
  memset(approxInvSqrt, 0, lenInvSqrt*sizeof(limb));
  // Initialize approximate inverse square root.
  invSqrt = LIMB_RANGE / sqrt(getMantissa(adjustedArgument+len, len)*LIMB_RANGE+1);
  approxInvSqrt[lenInvSqrt - 1].x = 1;
  invSqrt = (invSqrt - 1)*LIMB_RANGE;
  if (invSqrt > MAX_INT_NBR)
  {
    approxInvSqrt[lenInvSqrt - 2].x = MAX_INT_NBR;
  }
  else
  {
    approxInvSqrt[lenInvSqrt - 2].x = (int)invSqrt;
  }
               // Perform Newton approximation loop.
               // Get bit length of each cycle.
  bitLengthNbrCycles = 0;
  bitLength = lenInvSqrt*BITS_PER_GROUP;
  while (bitLength > 7)
  {
    bitLengthCycle[bitLengthNbrCycles++] = bitLength;
    bitLength = (bitLength+1) >> 1;
  }
               // Each loop increments precision.
  while (--bitLengthNbrCycles >= 0)
  {
    int limbLength, prevLimb;
    bitLength = bitLengthCycle[bitLengthNbrCycles];
    limbLength = (bitLength + 3*(BITS_PER_GROUP)-1) / BITS_PER_GROUP;
    if (limbLength > lenInvSqrt)
    {
      limbLength = lenInvSqrt;
    }     
    // Point to least significant limb for this cycle.
    ptrApproxInvSqrt = &approxInvSqrt[lenInvSqrt - limbLength];
    // Compute x(3-Nxx)/2. Start by squaring.
    multiply(ptrApproxInvSqrt, ptrApproxInvSqrt, approxInv, limbLength, NULL);
    // Multiply by argument.
    multiply(&approxInv[limbLength-1], &adjustedArgument[len - limbLength], arrAux, limbLength, NULL);
    // Subtract arrAux from 3.
    ptrArrAux = &arrAux[limbLength];
    for (idx = limbLength-1; idx > 0; idx--)
    {
      ptrArrAux->x ^= MAX_VALUE_LIMB;
      ptrArrAux++;
    }
    ptrArrAux->x = 2 - ptrArrAux->x;
    // Divide arrAux by 2.
    prevLimb = 0;
    for (idx = limbLength; idx > 0; idx--)
    {
      currLimb = ptrArrAux->x;
      (ptrArrAux--)->x = ((currLimb >> 1) | (prevLimb << (BITS_PER_GROUP - 1))) & MAX_VALUE_LIMB;
      prevLimb = currLimb;
    }
    // Multiply arrAux by approxInvSqrt.
    multiply(ptrArrAux+1, &approxInvSqrt[lenInvSqrt-limbLength], approxInv, limbLength, NULL);
    memcpy(&approxInvSqrt[lenInvSqrt - limbLength], &approxInv[limbLength - 1], limbLength*sizeof(limb));
  }
  // Multiply approxInvSqrt by argument to obtain the square root.
  if (len > lenInvSqrt)
  {
    multiply(approxInvSqrt, &adjustedArgument[len - lenInvSqrt], approxInv, lenInvSqrt, NULL);
  }
  else
  {
    multiply(approxInvSqrt, adjustedArgument, approxInv, len, NULL);
  }             // approxInv holds the square root.
  lenInvSqrt2 = ((len+1) >> 1);
  if (approxInv[2*lenInvSqrt - lenInvSqrt2-2].x > (7 << (BITS_PER_GROUP-3)))
  {                   // Increment square root.
    for (idx = 2 * lenInvSqrt - lenInvSqrt2-1; idx < 2*lenInvSqrt-1; idx++)
    {
      approxInv[idx].x = (approxInv[idx].x + 1) & MAX_INT_NBR;
      if (approxInv[idx].x != 0)
      {
        break;
      }
    }
    if (idx == 2*lenInvSqrt-1)
    {                // Roll back on overflow.
      for (idx = 2 * lenInvSqrt - lenInvSqrt2-1; idx < 2 * lenInvSqrt-1; idx++)
      {
        if (--approxInv[idx].x >= 0)
        {
          break;
        }
        approxInv[idx].x = MAX_VALUE_LIMB;
      }
    }
                    // Test whether square root is correct.
                    // It is correct only if when squared, it is <= than the argument.
    ptrApproxInvSqrt = &approxInv[2 * lenInvSqrt - lenInvSqrt2-1];
    multiply(ptrApproxInvSqrt, ptrApproxInvSqrt, approxInvSqrt, lenInvSqrt2, NULL);
    for (idx = lenInvSqrt-1; idx > 0; idx--)
    {
      if (adjustedArgument[idx].x != approxInvSqrt[idx].x)
      {
        break;
      }
    }
    if (adjustedArgument[idx].x < approxInvSqrt[idx].x)
    {                // Incorrect square root: roll back.
      for (idx = 2 * lenInvSqrt - lenInvSqrt2-1; idx < 2 * lenInvSqrt-1; idx++)
      {
        if (--approxInv[idx].x >= 0)
        {
          break;
        }
        approxInv[idx].x = MAX_VALUE_LIMB;
      }
    }
  }
  len = (len + 1) >> 1;
  *pLenSqRoot = len;
  // Shift right shRight bits into result.
  ptrSrc = &approxInv[2 * lenInvSqrt - *pLenSqRoot - 1 + len - 1];
  ptrDest = sqRoot + len - 1;
  prev.x = 0;
  for (index = len; index > 0; index--)
  {
    (ptrDest--)->x = ((prev.x << (BITS_PER_GROUP-shRight)) |
                      (ptrSrc->x >> shRight)) & MAX_VALUE_LIMB;
    prev.x = (ptrSrc--)->x;
  }
}

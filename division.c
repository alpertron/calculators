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
#include "expression.h"
#include <math.h>

extern limb approxInv[MAX_LEN];
extern limb adjustedArgument[MAX_LEN];
extern limb arrAux[MAX_LEN];
extern int bitLengthCycle[20];

// This routine uses Newton iteration: if x is an approximate inverse square root of N,
// a better approximation is: x(3-Nxx)/2. After the inverse square root is computed,
// the square root is found just by multiplying by N.
// The argument is multiplied by a power of 4 so the most significant limb is
// between LIMB_RANGE/4 and LIMB_RANGE - 1 and there is an even number of limbs.
// At the end of the calculation, the result is divided by the power of 2.

// All computations are done in little-endian notation.
// Find power of 2 that divides the number.
// output: pNbrLimbs = pointer to number of limbs
//         pPower2 = pointer to power of 2.
static void MultiplyBigNbrByMinPowerOf2(int *pPower2, limb *number, int len, limb *dest)
{
  limb mostSignficLimb, oldLimb, newLimb;
  int index2, mask, shLeft;
  limb *ptrDest, *ptrSrc;

  ptrSrc = number;
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
  ptrDest = dest;
  // Multiply number by this power.
  oldLimb.x = 0;
  for (index2 = len; index2 > 0; index2--)
  {
    newLimb.x = ptrDest->x;
    (ptrDest++)->x = ((newLimb.x << shLeft) |
      (oldLimb.x >> (BITS_PER_GROUP - shLeft))) & MAX_VALUE_LIMB;
    oldLimb.x = newLimb.x;
  }
  ptrDest->x = oldLimb.x >> (BITS_PER_GROUP - shLeft);
  *pPower2 = shLeft;
}

// After computing the number of limbs of the results, this routine finds the inverse
// of the divisor and then multiplies it by the dividend using nbrLimbs+1 limbs.
// After that, the quotient is adjusted.
enum eExprErr BigIntDivide(BigInteger *pDividend, BigInteger *pDivisor, BigInteger *pQuotient)
{
  limb carry;
  double inverse;
  limb oldLimb, newLimb;
  int nbrLimbs, nbrLimbsDividend, nbrLimbsDivisor;

  // Check whether the divisor is zero.
  if (pDivisor->limbs[0].x == 0 && pDivisor->nbrLimbs == 1)
  {  // Indicate overflow if divisor is zero.
    return EXPR_NUMBER_TOO_HIGH;
  }
  // Get number of limbs of quotient.
  nbrLimbsDividend = pDividend->nbrLimbs;
  nbrLimbsDivisor = pDivisor->nbrLimbs;
  nbrLimbs = nbrLimbsDividend - nbrLimbsDivisor;
  if (nbrLimbs < 0)
  {   // Absolute value of dividend is less than absolute value of divisor.
    pQuotient->limbs[0].x = 0;
    pQuotient->nbrLimbs = 1;
    pQuotient->sign = SIGN_POSITIVE;
    return EXPR_OK;
  }
  if (nbrLimbs == 0)
  {   // Both divisor and dividend have the same number of limbs.
    for (nbrLimbs = nbrLimbsDividend - 1; nbrLimbs > 0; nbrLimbs--)
    {
      if (pDividend->limbs[nbrLimbs].x != pDivisor->limbs[nbrLimbs].x)
      {
        break;
      }
    }
    if (pDividend->limbs[nbrLimbs].x < pDivisor->limbs[nbrLimbs].x)
    {   // Dividend is less than divisor, so quotient is zero.
      pQuotient->limbs[0].x = 0;
      pQuotient->nbrLimbs = 1;
      pQuotient->sign = SIGN_POSITIVE;
      return EXPR_OK;
    }
  }
  if (nbrLimbsDividend == 1)
  {   // If dividend is small, perform the division directly.
    pQuotient->limbs[0].x = pDividend->limbs[0].x / pDivisor->limbs[0].x;
    pQuotient->nbrLimbs = 1;
  }
  else
  {
    int index;
    int bitLength;
    int bitLengthNbrCycles;
    int idx;
    int nbrLimbsQuotient;
    int power2;
    limb *ptrDest;
    limb *ptrDivisor, *ptrDividend, *ptrQuotient;

    nbrLimbs += 2;    // Use this number of limbs for intermediate calculations.
    if (nbrLimbs > nbrLimbsDivisor)
    {
      memset(&adjustedArgument[0], 0, (nbrLimbs - nbrLimbsDivisor)*sizeof(limb));
      memcpy(&adjustedArgument[nbrLimbs - nbrLimbsDivisor], &pDivisor->limbs[0], nbrLimbsDivisor*sizeof(limb));
    }
    else
    {
      memcpy(&adjustedArgument[0], &pDivisor->limbs[nbrLimbsDivisor - nbrLimbs], nbrLimbs*sizeof(limb));
    }
    MultiplyBigNbrByMinPowerOf2(&power2, adjustedArgument, nbrLimbs, adjustedArgument);
    // Initialize approximate inverse.
    inverse = MAX_VALUE_LIMB / ((double)adjustedArgument[nbrLimbs - 1].x + 1);
    approxInv[nbrLimbs-1].x = 1;
    approxInv[nbrLimbs-2].x = (int)floor((inverse-1)*MAX_VALUE_LIMB);
    // Perform Newton approximation loop.
    // Get bit length of each cycle.
    bitLengthNbrCycles = 0;
    bitLength = nbrLimbs*BITS_PER_GROUP;
    while (bitLength > BITS_PER_GROUP)
    {
      bitLengthCycle[bitLengthNbrCycles++] = bitLength;
      bitLength = (bitLength + 1) >> 1;
    }
    // Each loop increments precision.
    // Use Newton iteration: x_{n+1} = x_n(2 - x_n)
    while (--bitLengthNbrCycles >= 0)
    {
      limb *ptrArrAux;
      int limbLength;

      bitLength = bitLengthCycle[bitLengthNbrCycles];
      limbLength = (bitLength + 3 * (BITS_PER_GROUP)-1) / BITS_PER_GROUP;
      if (limbLength > nbrLimbs)
      {
        limbLength = nbrLimbs;
      }
      // Compute x(2-Nx).
      // Multiply by argument.
      multiply(&approxInv[nbrLimbs-limbLength], &adjustedArgument[nbrLimbs - limbLength], arrAux, limbLength, NULL);
      // Subtract arrAux from 2.
      carry.x = 0;
      ptrArrAux = &arrAux[limbLength];
      for (idx = limbLength - 1; idx > 0; idx--)
      {
        carry.x = -ptrArrAux->x + carry.x;
        (ptrArrAux++)->x = carry.x & MAX_VALUE_LIMB;
        carry.x >>= BITS_PER_GROUP;
      }
      ptrArrAux->x = 2 - ptrArrAux->x + carry.x;
      // Multiply arrAux by approxInv.
      multiply(&arrAux[limbLength], &approxInv[nbrLimbs - limbLength], approxInv, limbLength, NULL);
      memmove(&approxInv[nbrLimbs - limbLength], &approxInv[limbLength - 1], limbLength*sizeof(limb));
    }
    // Multiply approxInv by argument to obtain the quotient.
    if (nbrLimbsDividend >= nbrLimbs)
    {
      multiply(&pDividend->limbs[nbrLimbsDividend - nbrLimbs], approxInv, approxInv, nbrLimbs, NULL);
    }
    else
    {
      memset(arrAux, 0, (nbrLimbs - nbrLimbsDividend)*sizeof(limb));
      memcpy(&arrAux[nbrLimbs - nbrLimbsDividend], pDividend->limbs, nbrLimbsDividend*sizeof(limb));
      multiply(arrAux, approxInv, approxInv, nbrLimbs, NULL);
    }             // approxInv holds the quotient.
    // Shift left quotient power2 bits into result.
    ptrDest = &approxInv[nbrLimbs-1];
    oldLimb.x = 0;
    for (index = nbrLimbs; index >= 0; index--)
    {
      newLimb.x = ptrDest->x;
      (ptrDest++)->x = ((newLimb.x << power2) |
                        (oldLimb.x >> (BITS_PER_GROUP-power2))) & MAX_VALUE_LIMB;
      oldLimb.x = newLimb.x;
    }
    // Determine number of limbs of quotient.
    nbrLimbsQuotient = nbrLimbsDividend - nbrLimbsDivisor;
    ptrDivisor = &pDivisor->limbs[nbrLimbsDivisor];
    ptrDividend = &pDividend->limbs[nbrLimbsDividend];
    for (idx = nbrLimbsDivisor; idx > 0; idx--)
    {
      if ((--ptrDividend)->x != (--ptrDivisor)->x)
      {
        break;
      }
    }
    ptrQuotient = &approxInv[2 * nbrLimbs - nbrLimbsQuotient - 1];
    if ((idx > 0 && ptrDividend->x > ptrDivisor->x) || idx <= 0)
    {
      nbrLimbsQuotient++;
    }
    if ((ptrQuotient - 1)->x > (7 << (BITS_PER_GROUP - 3)))
    {                   // Increment quotient.
      for (idx = 0; idx < nbrLimbsQuotient; idx++)
      {
        if ((unsigned int)++((ptrQuotient+idx)->x) < LIMB_RANGE)
        {
          break;
        }
        (ptrQuotient + idx)->x = 0;
      }
      if (idx == nbrLimbsQuotient)
      {                // Roll back on overflow.
        for (idx = 0; idx < nbrLimbsQuotient; idx++)
        {
          if (--((ptrQuotient + idx)->x) >= 0)
          {
            break;
          }
          (ptrQuotient + idx)->x = MAX_VALUE_LIMB;
        }
      }
      // Test whether the quotient is correct.
      // It is correct only if multiplied by the divisor, it is <= than the dividend.
      if (nbrLimbsQuotient > nbrLimbsDivisor)
      {
        memcpy(&approxInv[0], pDivisor->limbs, nbrLimbsDivisor*sizeof(limb));
        memset(&approxInv[nbrLimbsDivisor], 0, (nbrLimbsQuotient - nbrLimbsDivisor)*sizeof(limb));
        multiply(&approxInv[0], &approxInv[2 * nbrLimbs - nbrLimbsQuotient], arrAux, nbrLimbsQuotient, NULL);
      }
      else
      {
        memset(&approxInv[2 * nbrLimbs], 0, (nbrLimbsDivisor - nbrLimbsQuotient)*sizeof(limb));
        multiply(pDivisor->limbs, &approxInv[2 * nbrLimbs - nbrLimbsQuotient], arrAux, nbrLimbsDivisor, NULL);
      }
      ptrDividend = &pDividend->limbs[0];
      ptrDest = &arrAux[0];
      carry.x = 0;
      for (idx = pDividend->nbrLimbs; idx > 0; idx--)
      {
        carry.x += (ptrDividend++)->x - (ptrDest++)->x;
        carry.x >>= BITS_PER_GROUP;
      }
      ptrDest = &pQuotient->limbs[0];
      if (carry.x != 0)
      {  // Decrement quotient.
        carry.x = -1;
        for (idx = 0; idx < nbrLimbsQuotient; idx++)
        {
          carry.x += (ptrQuotient++)->x;
          (ptrDest++)->x = carry.x & MAX_VALUE_LIMB;
          carry.x >>= BITS_PER_GROUP;
        }
        if (carry.x != 0)
        {
          pQuotient->nbrLimbs = nbrLimbsQuotient - 1;
        }
        else
        {
          pQuotient->nbrLimbs = nbrLimbsQuotient;
        }
      }
      else
      {
        memcpy(ptrDest, ptrQuotient, nbrLimbsQuotient*sizeof(limb));
        pQuotient->nbrLimbs = nbrLimbsQuotient;
      }
    }
    else
    {
      memcpy(&pQuotient->limbs[0], ptrQuotient, nbrLimbsQuotient*sizeof(limb));
      pQuotient->nbrLimbs = nbrLimbsQuotient;
    }
  }
  if (pDividend->sign == pDivisor->sign || (pQuotient->limbs[0].x == 0 && pQuotient->nbrLimbs == 1))
  {
    pQuotient->sign = SIGN_POSITIVE;
  }
  else
  {
    pQuotient->sign = SIGN_NEGATIVE;
  }
  return EXPR_OK;
}

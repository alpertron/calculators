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
// Multiply by a power of 2 such that the most significant bit of most significant limb is one.
// output: pNbrLimbs = pointer to number of limbs
//         pPower2 = pointer to power of 2.
static void MultiplyBigNbrByMinPowerOf2(int *pPower2, const limb *number, int len, limb *dest)
{
  limb mostSignficLimb;
  limb oldLimb;
  limb newLimb;
  unsigned int shLeft;
  unsigned int shRight;
  limb *ptrDest;
  unsigned int uiOld;
  unsigned int uiDest;

  shLeft = 0U;
  mostSignficLimb.x = (number + len - 1)->x;
  for (unsigned int mask = LIMB_RANGE/2U; mask > 0U; mask >>= 1)
  {
    if (((unsigned int)mostSignficLimb.x & mask) != 0U)
    {
      break;
    }
    shLeft++;
  }
  shRight = (unsigned int)BITS_PER_GROUP - shLeft;
  ptrDest = dest;
  // Multiply number by this power.
  oldLimb.x = 0;
  for (int index2 = len; index2 > 0; index2--)
  {
    unsigned int uiNew;
    newLimb.x = ptrDest->x;
    uiNew = (unsigned int)newLimb.x;
    uiOld = (unsigned int)oldLimb.x;
    uiDest = ((uiNew << shLeft) | (uiOld >> shRight)) & MAX_VALUE_LIMB;
    ptrDest->x = (int)uiDest;
    ptrDest++;
    oldLimb.x = newLimb.x;
  }
  uiOld = (unsigned int)oldLimb.x;
  uiDest = uiOld >> shRight;
  ptrDest->x = (int)uiDest;
  *pPower2 = shLeft;
}

// After computing the number of limbs of the results, this routine finds the inverse
// of the divisor and then multiplies it by the dividend using nbrLimbs+1 limbs.
// After that, the quotient is adjusted.
enum eExprErr BigIntDivide(const BigInteger *pDividend, const BigInteger *pDivisor, BigInteger *pQuotient)
{
  double inverse;
  limb oldLimb;
  limb newLimb;
  int nbrLimbs;
  int nbrLimbsDividend;
  int nbrLimbsDivisor;

  // Check whether the divisor is zero.
  if (BigIntIsZero(pDivisor))
  {  // Indicate overflow if divisor is zero.
    return EXPR_DIVIDE_BY_ZERO;
  }
  // Get number of limbs of quotient.
  nbrLimbsDividend = pDividend->nbrLimbs;
  nbrLimbsDivisor = pDivisor->nbrLimbs;
  nbrLimbs = nbrLimbsDividend - nbrLimbsDivisor;
  if (nbrLimbs < 0)
  {   // Absolute value of dividend is less than absolute value of divisor.
    intToBigInteger(pQuotient, 0);
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
      intToBigInteger(pQuotient, 0);
      return EXPR_OK;
    }
  }
  if (nbrLimbsDividend == 1)
  {   // If dividend is small, perform the division directly.
    pQuotient->limbs[0].x = pDividend->limbs[0].x / pDivisor->limbs[0].x;
    pQuotient->nbrLimbs = 1;
  }
  else if (nbrLimbsDivisor == 1)
  {   // Divisor is small: use divide by int.
      // Sign of quotient is determined later.
      // pDivisor may be overwritten by dividend in next copy.
    int divisor = pDivisor->limbs[0].x;
    if (pQuotient != pDividend)
    {
      CopyBigInt(pQuotient, pDividend);
    }
    subtractdivide(pQuotient, 0, divisor);
  }
  else if (nbrLimbsDivisor < 64)
  {   // It is faster to perform classical division than
      // using Newton algorithm.
      // Use adjustedArgument to hold the remainder.
#ifdef _USING64BITS_
    int64_t carry;
#else
    int carry;
    double dVal = 1.0 / (double)LIMB_RANGE;
    double dSquareLimb = (double)LIMB_RANGE * (double)LIMB_RANGE;
#endif
    int i;
    int TrialQuotient;
    double dNbr;
    double dInvDivisor;
#ifndef _USING64BITS_
    double dTrialQuotient;
    double dDelta;
#endif
    limb* ptrDividend;
    const limb* ptrDivisor;
    limb* ptrQuotient;
    int lenBytes = nbrLimbsDividend * (int)sizeof(limb);

    (void)memcpy(adjustedArgument, pDividend->limbs, lenBytes);
    adjustedArgument[nbrLimbsDividend].x = 0;
    pQuotient->nbrLimbs = nbrLimbsDividend - nbrLimbsDivisor + 1;
    ptrQuotient = &pQuotient->limbs[nbrLimbsDividend - nbrLimbsDivisor];
    dInvDivisor = 1/getMantissa(&pDivisor->limbs[nbrLimbsDivisor], nbrLimbsDivisor);
    for (; nbrLimbsDividend >= nbrLimbsDivisor; nbrLimbsDividend--)
    {
#ifndef _USING64BITS_
      int low;
      double dAccumulator;
#endif
      dNbr = getMantissa(&adjustedArgument[nbrLimbsDividend+1], nbrLimbsDividend+1)*LIMB_RANGE;
      TrialQuotient = (int)(unsigned int)floor((dNbr * dInvDivisor) + 0.5);
      if ((unsigned int)TrialQuotient >= LIMB_RANGE)
      {   // Maximum value for limb.
        TrialQuotient = (int)MAX_VALUE_LIMB;
      }
      // Compute Nbr <- Nbr - TrialQuotient * Modulus
#ifndef _USING64BITS_
      dTrialQuotient = (double)TrialQuotient;
      dDelta = 0.0;
#endif
      carry = 0;
      ptrDividend = &adjustedArgument[nbrLimbsDividend - nbrLimbsDivisor];
      ptrDivisor = pDivisor->limbs;
      for (i = 0; i < nbrLimbsDivisor; i++)
      {
#ifdef _USING64BITS_
        carry += (int64_t)ptrDividend->x - ((ptrDivisor->x * (int64_t)TrialQuotient));
        ptrDividend->x = UintToInt((unsigned int)carry & MAX_VALUE_LIMB);
        carry >>= BITS_PER_GROUP;
#else
        low = (ptrDividend->x - (ptrDivisor->x * TrialQuotient) + carry) & MAX_INT_NBR;
        // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
        // In that case, there would be an error of +/- 1.
        dAccumulator = (double)ptrDividend->x - ((double)ptrDivisor->x * dTrialQuotient) +
          (double)carry + dDelta;
        dDelta = 0.0;
        if (dAccumulator < 0.0)
        {
          dAccumulator += dSquareLimb;
          dDelta = -(double)LIMB_RANGE;
        }
        if (low < HALF_INT_RANGE)
        {
          carry = (int)floor((dAccumulator + (double)FOURTH_INT_RANGE) * dVal);
        }
        else
        {
          carry = (int)floor((dAccumulator - (double)FOURTH_INT_RANGE) * dVal);
        }
        ptrDividend->x = low;
#endif
        ptrDividend++;
        ptrDivisor++;
      }
#ifdef _USING64BITS_
      carry += (int64_t)ptrDividend->x;
      ptrDividend->x = UintToInt((unsigned int)carry & MAX_VALUE_LIMB);
      carry >>= BITS_PER_GROUP;
#else
      low = (ptrDividend->x + carry) & MAX_INT_NBR;
      // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
      // In that case, there would be an error of +/- 1.
      dAccumulator = (double)ptrDividend->x + (double)carry + dDelta;
      if (dAccumulator < 0.0)
      {
        dAccumulator += dSquareLimb;
      }
      if (low < HALF_INT_RANGE)
      {
        carry = (int)floor((dAccumulator + (double)FOURTH_INT_RANGE) * dVal);
      }
      else
      {
        carry = (int)floor((dAccumulator - (double)FOURTH_INT_RANGE) * dVal);
      }
      ptrDividend->x = low;
#endif
      ptrDividend++;
      ptrDividend->x = carry & MAX_INT_NBR;
      if (((unsigned int)adjustedArgument[nbrLimbsDividend].x & MAX_VALUE_LIMB) != 0U)
      {
        unsigned int cy = 0;
        ptrDividend = &adjustedArgument[nbrLimbsDividend - nbrLimbsDivisor];
        ptrDivisor = pDivisor->limbs;
        for (i = 0; i < nbrLimbsDivisor; i++)
        {
          cy += (unsigned int)(ptrDividend->x) + (unsigned int)(ptrDivisor->x);
          ptrDividend->x = UintToInt(cy & MAX_VALUE_LIMB);
          cy >>= BITS_PER_GROUP;
          ptrDivisor++;
          ptrDividend++;
        }
        TrialQuotient--;
        adjustedArgument[nbrLimbsDividend].x = 0;
        adjustedArgument[nbrLimbsDividend + 1].x = 0;
      }
      ptrQuotient->x = TrialQuotient;
      ptrQuotient--;
    }
    while (pQuotient->nbrLimbs > 1)
    {
      if (pQuotient->limbs[pQuotient->nbrLimbs - 1].x != 0)
      {
        break;
      }
      pQuotient->nbrLimbs--;
    }
  }
  else
  {        // Divisor has more than 64 limbs. Use Newton algorithm 
           // to find the inverse and then multiply by dividend.
    int bitLength;
    int bitLengthNbrCycles;
    int idx;
    int nbrLimbsQuotient;
    int power2;
    limb *ptrDest;
    const limb *ptrDivisor;
    const limb *ptrDividend;
    limb *ptrQuotient;
    limb *ptrQuot;
    int lenBytes;
    unsigned int power2U;
    unsigned int power2Complement;

    nbrLimbs += 3;    // Use this number of limbs for intermediate calculations.
    if (nbrLimbs > nbrLimbsDivisor)
    {
      lenBytes = (nbrLimbs - nbrLimbsDivisor) * (int)sizeof(limb);
      (void)memset(&adjustedArgument[0], 0, lenBytes);
      lenBytes = nbrLimbsDivisor * (int)sizeof(limb);
      (void)memcpy(&adjustedArgument[nbrLimbs - nbrLimbsDivisor], &pDivisor->limbs[0], lenBytes);
    }
    else
    {
      lenBytes = nbrLimbs * (int)sizeof(limb);
      (void)memcpy(&adjustedArgument[0], &pDivisor->limbs[nbrLimbsDivisor - nbrLimbs], lenBytes);
    }
    MultiplyBigNbrByMinPowerOf2(&power2, adjustedArgument, nbrLimbs, adjustedArgument);
    // Initialize approximate inverse.
    inverse = LIMB_RANGE / ((double)adjustedArgument[nbrLimbs - 1].x + 
          ((double)adjustedArgument[nbrLimbs - 2].x) / LIMB_RANGE);
    approxInv[nbrLimbs-1].x = 1;
    if (inverse <= 1.0)
    {
      approxInv[nbrLimbs - 2].x = 0;
    }
    else if (inverse == 2.0)
    {
      approxInv[nbrLimbs - 2].x = MAX_VALUE_LIMB;
      approxInv[nbrLimbs - 3].x = MAX_VALUE_LIMB;
    }
    else
    {
      double t = (inverse - 1.0) * (double)LIMB_RANGE;
      double floor_t = floor(t);
      approxInv[nbrLimbs - 2].x = (int)floor_t;
      approxInv[nbrLimbs - 3].x = (int)floor((t - floor_t) * (double)LIMB_RANGE);
    }
    // Perform Newton approximation loop.
    // Get bit length of each cycle.
    bitLengthNbrCycles = 0;
    bitLength = nbrLimbs*BITS_PER_GROUP;
    while (bitLength >= BITS_PER_GROUP)
    {
      bitLengthCycle[bitLengthNbrCycles] = bitLength;
      bitLengthNbrCycles++;
      bitLength = (bitLength + 1) / 2;
    }
    // Each loop increments precision.
    // Use Newton iteration: x_{n+1} = x_n * (2 - x_n)
    bitLengthNbrCycles--;
    while (bitLengthNbrCycles >= 0)
    {
      limb *ptrArrAux;
      int limbLength;

      bitLength = bitLengthCycle[bitLengthNbrCycles];
      limbLength = (bitLength + (3 * BITS_PER_GROUP)-1) / BITS_PER_GROUP;
      if (limbLength > nbrLimbs)
      {
        limbLength = nbrLimbs;
      }
      // Compute x(2-Nx).
      // Multiply by divisor.
      multiply(&approxInv[nbrLimbs-limbLength], &adjustedArgument[nbrLimbs - limbLength], arrAux, limbLength, NULL);
      // Subtract arrAux from 2.
      ptrArrAux = &arrAux[limbLength];
      for (idx = limbLength - 1; idx > 0; idx--)
      {
        ptrArrAux->x = MAX_INT_NBR - ptrArrAux->x;
        ptrArrAux++;
      }
      ptrArrAux->x = 1 - ptrArrAux->x;
      // Multiply arrAux by approxInv.
      multiply(&arrAux[limbLength], &approxInv[nbrLimbs - limbLength], approxInv, limbLength, NULL);
      lenBytes = limbLength * (int)sizeof(limb);
      (void)memmove(&approxInv[nbrLimbs - limbLength], &approxInv[limbLength - 1], lenBytes);
      bitLengthNbrCycles--;
    }
    // Multiply approxInv by argument to obtain the quotient.
    if (nbrLimbsDividend >= nbrLimbs)
    {
      multiply(&pDividend->limbs[nbrLimbsDividend - nbrLimbs], approxInv, approxInv, nbrLimbs, NULL);
    }
    else
    {
      lenBytes = (nbrLimbs - nbrLimbsDividend) * (int)sizeof(limb);
      (void)memset(arrAux, 0, lenBytes);
      lenBytes = nbrLimbsDividend * (int)sizeof(limb);
      (void)memcpy(&arrAux[nbrLimbs - nbrLimbsDividend], pDividend->limbs, lenBytes);
      multiply(arrAux, approxInv, approxInv, nbrLimbs, NULL);
    }             // approxInv holds the quotient.
    // Shift left quotient power2 bits into result.
    ptrDest = &approxInv[nbrLimbs - 1];
    oldLimb.x = 0;
    power2U = (unsigned int)power2;
    power2Complement = (unsigned int)BITS_PER_GROUP - power2U;
    for (int index = nbrLimbs; index >= 0; index--)
    {
      newLimb.x = ptrDest->x;
      ptrDest->x = UintToInt((((unsigned int)newLimb.x << power2U) |
        ((unsigned int)oldLimb.x >> power2Complement)) & MAX_VALUE_LIMB);
      ptrDest++;
      oldLimb.x = newLimb.x;
    }

    // Determine number of limbs of quotient.
    nbrLimbsQuotient = nbrLimbsDividend - nbrLimbsDivisor;
    ptrDivisor = &pDivisor->limbs[nbrLimbsDivisor - 1];
    ptrDividend = &pDividend->limbs[nbrLimbsDividend - 1];
    for (idx = nbrLimbsDivisor - 1; idx > 0; idx--)
    {
      if (ptrDividend->x != ptrDivisor->x)
      {
        break;
      }
      ptrDividend--;
      ptrDivisor--;
    }
    if (ptrDividend->x >= ptrDivisor->x)
    {
      nbrLimbsQuotient++;
    }
    ptrQuotient = &approxInv[(2 * nbrLimbs) - nbrLimbsQuotient];
    if (approxInv[(2 * nbrLimbs) - 1].x == 0)
    {  // Most significant byte is zero, so it is not part of the quotient. 
      ptrQuotient--;
    }
    ptrQuot = ptrQuotient;
    if ((unsigned int)(ptrQuotient - 1)->x > (LIMB_RANGE / 8U * 7U))
    {                   // Increment quotient.
      for (idx = 0; idx <= nbrLimbsQuotient; idx++)
      {
        (ptrQuotient + idx)->x ++;
        if (((unsigned int)(ptrQuotient + idx)->x & MAX_INT_NBR_U) != 0U)
        {
          break;
        }
        (ptrQuotient + idx)->x = 0;
      }
      if ((idx >= nbrLimbsQuotient) && (ptrDividend->x < ptrDivisor->x))
      {                // Roll back on overflow.
        for (idx = 0; idx <= nbrLimbsQuotient; idx++)
        {
          ((ptrQuotient + idx)->x)--;
          if (((ptrQuotient + idx)->x) >= 0)
          {
            break;
          }
          (ptrQuotient + idx)->x = MAX_VALUE_LIMB;
        }
      }
      if (approxInv[(2 * nbrLimbs) - 1].x != 0)
      {    // Most significant byte is not zero, so it is part of the quotient.
        ptrQuot = &approxInv[(2 * nbrLimbs) - nbrLimbsQuotient];
      }
      // Test whether the quotient is correct.
      // It is correct only if multiplied by the divisor, it is <= than the dividend.
      if (nbrLimbsQuotient > nbrLimbsDivisor)
      {
        lenBytes = nbrLimbsDivisor * (int)sizeof(limb);
        (void)memcpy(&approxInv[0], pDivisor->limbs, lenBytes);
        lenBytes = (nbrLimbsQuotient - nbrLimbsDivisor) * (int)sizeof(limb);
        (void)memset(&approxInv[nbrLimbsDivisor], 0, lenBytes);
        multiply(&approxInv[0], ptrQuot, arrAux, nbrLimbsQuotient, NULL);
      }
      else
      {
        lenBytes = (nbrLimbsDivisor - nbrLimbsQuotient) * (int)sizeof(limb);
        (void)memset(&approxInv[2 * nbrLimbs], 0, lenBytes);
        multiply(pDivisor->limbs, ptrQuot, arrAux, nbrLimbsDivisor, NULL);
      }
      ptrDividend = &pDividend->limbs[pDividend->nbrLimbs - 1];
      ptrDest = &arrAux[pDividend->nbrLimbs - 1];
      for (idx = pDividend->nbrLimbs - 1; idx > 0; idx--)
      {
        if (ptrDividend->x != ptrDest->x)
        {
          break;
        }
        ptrDividend--;
        ptrDest--;
      }
      if (ptrDividend->x < ptrDest->x)
      {  // Decrement quotient.
        ptrQuotient = ptrQuot;
        for (idx = 0; idx < nbrLimbsQuotient; idx++)
        {
          (ptrQuotient->x)--;
          if (ptrQuotient->x >= 0)
          {
            break;
          }
          ptrQuotient->x = MAX_VALUE_LIMB;
          ptrQuotient++;
        }
        if (idx == nbrLimbsQuotient)
        {
          nbrLimbsQuotient--;
        }
      }
    }
    lenBytes = nbrLimbsQuotient * (int)sizeof(limb);
    (void)memcpy(&pQuotient->limbs[0], ptrQuot, lenBytes);
    pQuotient->nbrLimbs = nbrLimbsQuotient;
  }
  if ((pDividend->sign == pDivisor->sign) ||
    ((pQuotient->limbs[0].x == 0) && (pQuotient->nbrLimbs == 1)))
  {
    pQuotient->sign = SIGN_POSITIVE;
  }
  else
  {
    pQuotient->sign = SIGN_NEGATIVE;
  }
  return EXPR_OK;
}

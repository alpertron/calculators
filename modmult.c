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
#include <math.h>
#include <stdint.h>
#include "bignbr.h"
limb MontgomeryR1[MAX_LEN];
limb TestNbr[MAX_LEN];
limb MontgomeryMultN[MAX_LEN];
limb MontgomeryMultR1[MAX_LEN];
limb MontgomeryMultR2[MAX_LEN];
static limb aux[MAX_LEN], aux2[MAX_LEN];
static limb aux3[MAX_LEN], aux4[MAX_LEN];
static int NumberLength2;
int NumberLength, NumberLengthR1;
long long lModularMult;
mmCback modmultCallback;
static limb U[MAX_LEN], V[MAX_LEN], R[MAX_LEN], S[MAX_LEN];
static BigInteger tmpDen, tmpNum;

int getNbrLimbs(limb *bigNbr)
{
  limb *ptrLimb = bigNbr + NumberLength;
  while (ptrLimb > bigNbr)
  {
    if ((--ptrLimb)->x != 0)
    {
      return (int)(ptrLimb - &bigNbr[0] + 1);
    }
  }
  return 1;
}

// Find the inverse of value mod 2^(NumberLength*BITS_PER_GROUP)
void ComputeInversePower2(limb *value, limb *result, limb *tmp)
{
  int N, x, j;
  limb Cy;
  int currLen;
  x = N = (int)value->x;       // 2 least significant bits of inverse correct.
  x = x * (2 - N * x);         // 4 least significant bits of inverse correct.
  x = x * (2 - N * x);         // 8 least significant bits of inverse correct.
  x = x * (2 - N * x);         // 16 least significant bits of inverse correct.
  x = x * (2 - N * x);         // 32 least significant bits of inverse correct.
  result->x = x & MAX_VALUE_LIMB;
  for (currLen = 2; currLen <= NumberLength * 2; currLen <<= 1)
  {
    multiply(value, result, tmp, currLen, NULL);    // tmp <- N * x
    Cy.x = 2 - tmp[0].x;
    tmp[0].x = Cy.x & MAX_VALUE_LIMB;
    for (j = 1; j < currLen; j++)
    {
      Cy.x = (Cy.x >> BITS_PER_GROUP) - tmp[j].x;
      tmp[j].x = Cy.x & MAX_VALUE_LIMB;
    }                                               // tmp <- 2 - N * x
    multiply(result, tmp, result, currLen, NULL);   // tmp <- x * (2 - N * x)
  }
}

// Let R be a power of 2 of at least len limbs.
// Compute R1 = MontgomeryR1 and N = MontgomeryN using the formulas:
// R1 = R mod M
// N = -M^(-1) mod R

void GetMontgomeryParms(int len)
{
  int j;
  limb Cy;
  TestNbr[len].x = 0;
  NumberLength = len;
  NumberLength2 = len + len;
  if (NumberLength == 1)
  {
    MontgomeryMultR1[0].x = 1;
    MontgomeryMultR2[0].x = 1;
    NumberLengthR1 = 1;
    return;
  }
  // Compute MontgomeryMultN as -1/TestNbr (mod 2^k) using Newton method,
  // which doubles the precision for each iteration.
  // In the formula above: k = BITS_PER_GROUP * NumberLength.
  if (NumberLength >= 8)
  {
    limb *ptrResult;
    ComputeInversePower2(TestNbr, MontgomeryMultN, aux);
    ptrResult = &MontgomeryMultN[0];
    Cy.x = 0;          // Change sign.
    for (j = 0; j < NumberLength; j++)
    {
      Cy.x = (Cy.x >> BITS_PER_GROUP) - ptrResult->x;
      ptrResult->x = Cy.x & MAX_VALUE_LIMB;
      ptrResult++;
    }
    ptrResult->x = 0;
  }
  else
  {
    int x, N;
    x = N = (int)TestNbr[0].x;   // 2 least significant bits of inverse correct.
    x = x * (2 - N * x);         // 4 least significant bits of inverse correct.
    x = x * (2 - N * x);         // 8 least significant bits of inverse correct.
    x = x * (2 - N * x);         // 16 least significant bits of inverse correct.
    x = x * (2 - N * x);         // 32 least significant bits of inverse correct.
    MontgomeryMultN[0].x = (-x) & MAX_VALUE_LIMB;    // Change sign
  }
  // Compute MontgomeryMultR1 as 1 in Montgomery notation,
  // this is 2^(NumberLength*BITS_PER_GROUP) % TestNbr.
  j = NumberLength;
  MontgomeryMultR1[j].x = 1;
  do
  {
    MontgomeryMultR1[--j].x = 0;
  } while (j > 0);
  AdjustModN(MontgomeryMultR1, TestNbr, len);
  MontgomeryMultR1[NumberLength].x = 0;
  memcpy(MontgomeryMultR2, MontgomeryMultR1, (NumberLength+1)*sizeof(limb));
  for (NumberLengthR1 = NumberLength; NumberLengthR1 > 0; NumberLengthR1--)
  {
    if (MontgomeryMultR1[NumberLengthR1 - 1].x != 0)
    {
      break;
    }
  }
  // Compute MontgomeryMultR2 as 2^(2*NumberLength*BITS_PER_GROUP) % TestNbr.
  for (j = NumberLength; j > 0; j--)
  {
    memmove(&MontgomeryMultR2[1], &MontgomeryMultR2[0], NumberLength*sizeof(limb));
    MontgomeryMultR2[0].x = 0;
    AdjustModN(MontgomeryMultR2, TestNbr, len);
  }
}

// Compute Nbr <- Nbr mod TestNbr.
// TestNbr has NumberLength limbs.
void AdjustModN(limb *Nbr, limb *Modulus, int nbrLen)
{
  int i, carry;
  int TrialQuotient;
  double dNbr, dModulus, dTrialQuotient;
  double dAccumulator, dDelta;
  double dVal = 1 / (double)LIMB_RANGE;
  double dSquareLimb = (double)LIMB_RANGE * (double)LIMB_RANGE;

  dModulus = getMantissa(Modulus+nbrLen, nbrLen);
  dNbr = getMantissa(Nbr + nbrLen + 1, nbrLen + 1) * LIMB_RANGE;
  TrialQuotient = (int)(unsigned int)floor(dNbr / dModulus + 0.5);
  if ((unsigned int)TrialQuotient >= LIMB_RANGE)
  {   // Maximum value for limb.
    TrialQuotient = MAX_VALUE_LIMB;
  }
  // Compute Nbr <- Nbr - TrialQuotient * Modulus
  dTrialQuotient = (double)TrialQuotient;
  carry = 0;
  dAccumulator = 0;
  dDelta = 0;
  for (i = 0; i <= nbrLen; i++)
  {
    int low = (Nbr[i].x - Modulus[i].x * TrialQuotient + carry) & MAX_INT_NBR;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    dAccumulator = Nbr[i].x - Modulus[i].x * dTrialQuotient + carry + dDelta;
    dDelta = 0;
    if (dAccumulator < 0)
    {
      dAccumulator += dSquareLimb;
      dDelta = -(double)LIMB_RANGE;
    }
    if (low < HALF_INT_RANGE)
    {
      carry = (int)floor((dAccumulator + HALF_INT_RANGE / 2)*dVal);
    }
    else
    {
      carry = (int)floor((dAccumulator - HALF_INT_RANGE / 2)*dVal);
    }
    Nbr[i].x = low;
  }
  Nbr[i].x = carry & MAX_INT_NBR;
  if ((Nbr[nbrLen].x & MAX_VALUE_LIMB) != 0)
  {
    unsigned int cy = 0;
    for (i = 0; i < nbrLen; i++)
    {
      cy += (unsigned int)Nbr[i].x + (unsigned int)Modulus[i].x;
      Nbr[i].x = (int)(cy & MAX_VALUE_LIMB);
      cy >>= BITS_PER_GROUP;
    }
    Nbr[nbrLen].x = 0;
  }
}

void AddBigNbrModN(limb *Nbr1, limb *Nbr2, limb *Sum, limb *mod, int nbrLen)
{
  unsigned int carry;
  int borrow;
  int i;

  carry = 0;
  for (i = 0; i < nbrLen; i++)
  {
    carry = (carry >> BITS_PER_GROUP) +
      (unsigned int)(Nbr1 + i)->x + (unsigned int)(Nbr2 + i)->x;
    Sum[i].x = (int)(carry & MAX_VALUE_LIMB);
  }
  borrow = 0;
  for (i = 0; i < nbrLen; i++)
  {
    borrow = (borrow >> BITS_PER_GROUP) +
       (unsigned int)Sum[i].x - (mod + i)->x;
    Sum[i].x = (int)(borrow & MAX_VALUE_LIMB);
  }

  if (carry < LIMB_RANGE && borrow < 0)
  {
    carry = 0;
    for (i = 0; i < nbrLen; i++)
    {
      carry = (carry >> BITS_PER_GROUP) +
          (unsigned int)(Sum+i)->x + (unsigned int)(mod+i)->x;
      Sum[i].x = (int)(carry & MAX_VALUE_LIMB);
    }
  }
}

void AddBigNbrMod(limb *Nbr1, limb *Nbr2, limb *Sum)
{
  AddBigNbrModN(Nbr1, Nbr2, Sum, TestNbr, NumberLength);
}

void SubtBigNbrModN(limb *Nbr1, limb *Nbr2, limb *Diff, limb *mod, int nbrLen)
{
  int i;
  int borrow = 0;
  for (i = 0; i < nbrLen; i++)
  {
    borrow = (borrow >> BITS_PER_GROUP) + (Nbr1 + i)->x - (Nbr2 + i)->x;
    Diff[i].x = (int)(borrow & MAX_VALUE_LIMB);
  }
  if (borrow < 0)
  {
    unsigned int carry = 0;
    for (i = 0; i < nbrLen; i++)
    {
      carry = (carry >> BITS_PER_GROUP) + 
          (unsigned int)(Diff + i)->x + (unsigned int)(mod + i)->x;
      Diff[i].x = (int)(carry & MAX_VALUE_LIMB);
    }
  }
}

void SubtBigNbrMod(limb *Nbr1, limb *Nbr2, limb *Sum)
{
  SubtBigNbrModN(Nbr1, Nbr2, Sum, TestNbr, NumberLength);
}

static void smallmodmult(int factor1, int factor2, limb *product, int mod)
{
  if (mod < SMALL_NUMBER_BOUND)
  {
    product->x = factor1 * factor2 % mod;
  }
  else
  {   // TestNbr has one limb but it is not small.
#ifdef _USING64BITS_
    product->x = (int64_t)factor1 * factor2 % mod;
#else
      // Round up quotient.
    int quotient = (int)floor((double)factor1 * (double)factor2 / (double)mod + 0.5);
    int remainder = factor1 * factor2 - quotient * mod;
    if (remainder < 0)
    {    // Quotient was 1 more than expected. Adjust remainder.
      remainder += mod;
    }
    product->x = remainder;
#endif
  }
}

// Multiply two numbers in Montgomery notation.
//
// For large numbers the REDC algorithm is:
// m <- ((T mod R)N') mod R
// t <- (T + mN) / R
// if t >= N then
//   return t - N
// else
//   return t
// end if

void modmult(limb *factor1, limb *factor2, limb *product)
{
  limb carry;
  int count;
  limb Prod[13];
  unsigned int cy;
  int index;
#ifdef __EMSCRIPTEN__
  if (modmultCallback)
  {
    modmultCallback();
    lModularMult++;
  }
#endif
  if (NumberLength == 1)
  {
    smallmodmult(factor1->x, factor2->x, product, TestNbr[0].x);
    return;
  }
  if (NumberLength <= 12)
  {     // Small numbers.
    int i, j;
#ifdef _USING64BITS_
    int32_t MontDig, Nbr;
    int64_t Pr;
    memset(Prod, 0, NumberLength * sizeof(limb));
    for (i = 0; i < NumberLength; i++)
    {
      Pr = (Nbr = (factor1 + i)->x) * (int64_t)factor2->x + (uint32_t)Prod[0].x;
      MontDig = ((int32_t)Pr * MontgomeryMultN[0].x) & MAX_VALUE_LIMB;
      Prod[0].x = (Pr = (((int64_t)MontDig * TestNbr[0].x + Pr) >> BITS_PER_GROUP) +
        (int64_t)MontDig * TestNbr[1].x + (int64_t)Nbr * (factor2 + 1)->x + (uint32_t)Prod[1].x) & MAX_VALUE_LIMB;
      for (j = 2; j < NumberLength; j++)
      {
        Prod[j - 1].x = ((Pr = (Pr >> BITS_PER_GROUP) +
          (int64_t)MontDig * TestNbr[j].x + (int64_t)Nbr * (factor2 + j)->x + (uint32_t)Prod[j].x) & MAX_VALUE_LIMB);
      }
      Prod[j - 1].x = (int32_t)(Pr >> BITS_PER_GROUP);
    }
#else
    double dLimbRange = (double)LIMB_RANGE;
    double dInvLimbRange = (double)1 / dLimbRange;
    memset(Prod, 0, NumberLength*sizeof(limb));
    for (i = 0; i < NumberLength; i++)
    {
      int Nbr = (factor1 + i)->x;
      double dNbr = (double)Nbr;
      int low = Nbr * factor2->x + Prod[0].x;
      double dAccum = dNbr * (double)factor2->x + (double)Prod[0].x;
      int MontDig = (low * MontgomeryMultN[0].x) & MAX_VALUE_LIMB;
      double dMontDig = (double)MontDig;
      dAccum += dMontDig * (double)TestNbr[0].x;
      // At this moment dAccum is multiple of LIMB_RANGE.
      dAccum = floor(dAccum*dInvLimbRange + 0.5);
      low = ((unsigned int)dAccum + MontDig * TestNbr[1].x +
                   Nbr * (factor2 + 1)->x + Prod[1].x) & MAX_VALUE_LIMB;
      dAccum += dMontDig * TestNbr[1].x + dNbr * (factor2 + 1)->x + (unsigned int)Prod[1].x;
      Prod[0].x = low;
      for (j = 2; j < NumberLength; j++)
      {
        // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
        // In that case, there would be an error of +/- 1.
        if (low < HALF_INT_RANGE)
        {
          dAccum = ((dAccum + HALF_INT_RANGE / 2)*dInvLimbRange);
        }
        else
        {
          dAccum = ((dAccum - HALF_INT_RANGE / 2)*dInvLimbRange);
        }
        low = (int)(dAccum - floor(dAccum * dInvLimbRange) * dLimbRange);
        dAccum += dMontDig * TestNbr[j].x + dNbr * (factor2 + j)->x + (unsigned int)Prod[j].x;
        low = (low + MontDig * TestNbr[j].x +
               Nbr * (factor2 + j)->x + Prod[j].x) & MAX_VALUE_LIMB;
        Prod[j - 1].x = low;
      }
      if (low < HALF_INT_RANGE)
      {
        dAccum = ((dAccum + HALF_INT_RANGE / 2)*dInvLimbRange);
      }
      else
      {
        dAccum = ((dAccum - HALF_INT_RANGE / 2)*dInvLimbRange);
      }
      Prod[j - 1].x = (unsigned int)dAccum;  // Most significant limb can be greater than LIMB_RANGE
    }
#endif
    for (j = NumberLength - 1; j >= 0; j--)
    {
      if (Prod[j].x != TestNbr[j].x)
      {
        break;
      }
    }
    if (j<0 || (unsigned int)Prod[j].x >= (unsigned int)TestNbr[j].x)
    {        // Prod >= TestNbr, so perform Prod <- Prod - TestNbr
      carry.x = 0;
      for (count = 0; count < NumberLength; count++)
      {
        carry.x += Prod[count].x - TestNbr[count].x;
        Prod[count].x = carry.x & MAX_VALUE_LIMB;
        carry.x >>= BITS_PER_GROUP;
      }
    }
    memcpy(product, Prod, NumberLength*sizeof(limb));
    return;
  }
  // Compute T
  multiply(factor1, factor2, product, NumberLength, NULL);
  // Compute m
  multiply(product, MontgomeryMultN, aux, NumberLength, NULL);
  // Compute mN
  multiply(TestNbr, aux, aux2, NumberLength, NULL);
  // Check if lowest half of mN is not zero
  for (count = NumberLength-1; count >= 0; count--)
  {
    if (aux2[count].x != 0)
    {    // Lowest half of mN is not zero
      break;
    }
  }
  // If lowest half of mN is zero, compute hi(T) + hi(mN)
  // else compute hi(T) + hi(mN) + 1
  // Where hi(number) is the high half of number.
  cy = (count >= 0 ? LIMB_RANGE : 0);
  index = NumberLength;
  for (count = 0; count < NumberLength; count++)
  {
    cy = (cy >> BITS_PER_GROUP) +
      (unsigned int)(product + index)->x + (unsigned int)aux2[index].x;
      (product + count)->x = (int)(cy & MAX_VALUE_LIMB);
    index++;
  }
  // Check whether this number is greater than TestNbr.
  if (cy < LIMB_RANGE)
  {
    for (count = NumberLength - 1; count > 0; count--)
    {
      if ((product + count)->x != TestNbr[count].x)
      {
        break;
      }
    }
  }
  if (cy >= LIMB_RANGE || (product + count)->x >= TestNbr[count].x)
  {  // The number is greater or equal than Testnbr. Subtract it.
    int borrow = 0;
    for (count = 0; count < NumberLength; count++)
    {
      borrow = (borrow >> BITS_PER_GROUP) +
        (product + count)->x - TestNbr[count].x;
      (product + count)->x = (int)(borrow & MAX_VALUE_LIMB);
    }
  }
}

// Multiply big number in Montgomery notation by integer.
void modmultIntExtended(limb *factorBig, int factorInt, limb *result, limb *pTestNbr, int nbrLen)
{
  int i, low;
  int TrialQuotient;
  double dTrialQuotient, dAccumulator, dFactorInt;
  limb *ptrFactorBig, *ptrTestNbr;
  double dTestNbr, dFactorBig;
  double dInvLimbRange = 1 / (double)LIMB_RANGE;
  if (nbrLen == 1)
  {
    smallmodmult(factorBig->x, factorInt, result, pTestNbr->x);
    return;
  }
  (factorBig + nbrLen)->x = 0;
  dFactorInt = (double)factorInt;
  dTestNbr = getMantissa(pTestNbr + nbrLen, nbrLen);
  dFactorBig = getMantissa(factorBig + nbrLen, nbrLen);
  TrialQuotient = (int)(unsigned int)floor(dFactorBig * factorInt / dTestNbr + 0.5);
  if ((unsigned int)TrialQuotient >= LIMB_RANGE)
  {   // Maximum value for limb.
    TrialQuotient = MAX_VALUE_LIMB;
  }
  // Compute result <- factorBig * factorInt - TrialQuotient * TestNbr
  dTrialQuotient = (double)TrialQuotient;
  low = 0;
  dAccumulator = 0;
  ptrFactorBig = factorBig;
  ptrTestNbr = pTestNbr;
  for (i = 0; i <= nbrLen; i++)
  {
    dAccumulator += ptrFactorBig->x * dFactorInt - dTrialQuotient * ptrTestNbr->x;
    low += ptrFactorBig->x * factorInt - TrialQuotient * ptrTestNbr->x;
    low &= MAX_VALUE_LIMB;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    (result + i)->x = low;
    if (low < HALF_INT_RANGE)
    {
      dAccumulator = floor(dAccumulator*dInvLimbRange + 0.25);
    }
    else
    {
      dAccumulator = floor(dAccumulator*dInvLimbRange - 0.25);
    }
    low = (int)dAccumulator & MAX_VALUE_LIMB;
    ptrFactorBig++;
    ptrTestNbr++;
  }
  while (((result+nbrLen)->x & MAX_VALUE_LIMB) != 0)
  {
    ptrFactorBig = result;
    ptrTestNbr = pTestNbr;
    unsigned int cy = 0;
    for (i = 0; i <= nbrLen; i++)
    {
      cy += (unsigned int)ptrTestNbr->x + (unsigned int)ptrFactorBig->x;
      ptrFactorBig->x = (int)(cy & MAX_VALUE_LIMB);
      cy >>= BITS_PER_GROUP;
      ptrFactorBig++;
      ptrTestNbr++;
    }
  }
}

void modmultInt(limb *factorBig, int factorInt, limb *result)
{
  modmultIntExtended(factorBig, factorInt, result, TestNbr, NumberLength);
}

// Input: base = base in Montgomery notation.
//        exp  = exponent.
//        nbrGroupsExp = number of limbs of exponent.
// Output: power = power in Montgomery notation.
void modPow(limb *base, limb *exp, int nbrGroupsExp, limb *power)
{
  int mask, index;
  memcpy(power, MontgomeryMultR1, (NumberLength + 1)*sizeof(*power));  // power <- 1
  for (index = nbrGroupsExp - 1; index >= 0; index--)
  {
    int groupExp = (int)(exp + index)->x;
    for (mask = 1 << (BITS_PER_GROUP - 1); mask > 0; mask >>= 1)
    {
      modmult(power, power, power);
      if ((groupExp & mask) != 0)
      {
        modmult(power, base, power);
      }
    }
  }
}

// Input: base = base in Montgomery notation.
//        exp  = exponent.
// Output: power = power in Montgomery notation.
void modPowLimb(limb *base, limb *exp, limb *power)
{
  int mask, groupExp;
  memcpy(power, MontgomeryMultR1, (NumberLength + 1)*sizeof(*power));  // power <- 1
  groupExp = (int)exp->x;
  for (mask = 1 << (BITS_PER_GROUP - 1); mask > 0; mask >>= 1)
  {
    modmult(power, power, power);
    if ((groupExp & mask) != 0)
    {
      modmult(power, base, power);
    }
  }
}

void modPowBaseInt(int base, limb *exp, int nbrGroupsExp, limb *power)
{
  int mask, index;
  memcpy(power, MontgomeryMultR1, (NumberLength+1)*sizeof(limb));  // power <- 1
  for (index = nbrGroupsExp-1; index>=0; index--)
  {
    int groupExp = (int)(exp+index)->x;
    for (mask = 1<<(BITS_PER_GROUP-1); mask > 0; mask >>= 1)
    {
      modmult(power, power, power);
      if ((groupExp & mask) != 0)
      {
        modmultInt(power, base, power);
      }
    }
  }
}

/* U' <- eU + fV, V' <- gU + hV                                        */
/* U <- U', V <- V'                                                    */
static void AddMult(limb *firstBig, int e, int f, limb *secondBig, int g, int h, int nbrLen)
{
  double dVal = 1 / (double)LIMB_RANGE;
  int ctr, carryU, carryV;
  double dFactorE = (double)e;
  double dFactorF = (double)f;
  double dFactorG = (double)g;
  double dFactorH = (double)h;
  carryU = carryV = 0;
  for (ctr = 0; ctr <= nbrLen; ctr++)
  {
    int u = firstBig->x;
    int v = secondBig->x;
    int lowU = (carryU + u * e + v * f) & MAX_INT_NBR;
    int lowV = (carryV + u * g + v * h) & MAX_INT_NBR;
    // Subtract or add 0.25 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    double dCarry = ((double)carryU + (double)u * dFactorE +
                     (double)v * dFactorF)*dVal;
    if (lowU < HALF_INT_RANGE)
    {
      carryU = (int)floor(dCarry + 0.25);
    }
    else
    {
      carryU = (int)floor(dCarry - 0.25);
    }
    dCarry = ((double)carryV + (double)u * dFactorG +
              (double)v * dFactorH)*dVal;
    if (lowV < HALF_INT_RANGE)
    {
      carryV = (int)floor(dCarry + 0.25);
    }
    else
    {
      carryV = (int)floor(dCarry - 0.25);
    }
    (firstBig++)->x = lowU;
    (secondBig++)->x = lowV;
  }
}

// Perform first <- (first - second) / 2
// first must be greater than second.
static int HalveDifference(limb *first, limb *second, int len)
{
  int i;
  int borrow, prevLimb;
    // Perform first <- (first - second)/2.
  borrow = first->x - second->x;
  prevLimb = borrow & MAX_VALUE_LIMB;
  borrow >>= BITS_PER_GROUP;
  for (i = 1; i < len; i++)
  {
    int currLimb;
    borrow += (first + i)->x - (second + i)->x;
    currLimb = borrow & MAX_VALUE_LIMB;
    borrow >>= BITS_PER_GROUP;
    (first + i - 1)->x = ((prevLimb >> 1) |
      (currLimb << (BITS_PER_GROUP - 1))) & MAX_VALUE_LIMB;
    prevLimb = currLimb;
  }
  (first + i - 1)->x = prevLimb >> 1;
    // Get length of result.
  for (len--; len > 0; len--)
  {
    if ((first + len)->x != 0)
    {
      break;
    }
  }
  return len+1;
}

int modInv(int NbrMod, int currentPrime)
{
  int QQ, T1, T3;
  int V1 = 1;
  int V3 = NbrMod;
  int U1 = 0;
  int U3 = currentPrime;
  while (V3 != 0)
  {
    if (U3 < V3 + V3)
    {               // QQ = 1
      T1 = U1 - V1;
      T3 = U3 - V3;
    }
    else
    {
      QQ = U3 / V3;
      T1 = U1 - V1 * QQ;
      T3 = U3 - V3 * QQ;
    }
    U1 = V1;
    U3 = V3;
    V1 = T1;
    V3 = T3;
  }
  return U1 + (currentPrime & (U1 >> 31));
}

/***********************************************************************/
/* NAME: ModInvBigNbr                                                  */
/*                                                                     */
/* PURPOSE: Find the inverse multiplicative modulo M.                  */
/* The algorithm terminates with inv = X^(-1) mod M.                   */
/*                                                                     */
/* This routine uses Kaliski Montgomery inverse algorithm              */
/* with changes by E. Savas and C. K. Koc.                             */
/*  1. U <- M, V <- X, R <- 0, S <- 1, k <- 0                          */
/*  2. while V > 0 do                                                  */
/*  3.   if U even then U <- U / 2, S <- 2S                            */
/*  4.   elsif V even then V <- V / 2, R <- 2R                         */
/*  5.   elsif U > V  then U <- (U - V) / 2, R <- R + S, S <- 2S       */
/*  6.   else V <- (V - U) / 2, S <- S + R, R <- 2R                    */
/*  7.   k <- k + 1                                                    */
/*  8. if R >= M then R <- R - M                                       */
/*  9. R <- M - R                                                      */
/* 10. R <- MonPro(R, R2)                                              */
/* 11. return MonPro(R, 2^(m-k))                                       */
/*                                                                     */
/*  In order to reduce the calculations, several single precision      */
/*  variables are added:                                               */
/*                                                                     */
/* R' <- aR + bS, S' <- cR + dS                                        */
/***********************************************************************/
void ModInvBigNbr(limb *num, limb *inv, limb *mod, int nbrLen)
{
  int len;
  int k;
  int a, b, c, d;
  int size, i;
  int bitCount;
  int lenU, lenV;
  int borrow;
  if (nbrLen == 1)
  {
    inv->x = modInv(num->x, mod->x);
    return;
  }
  //  1. U <- M, V <- X, R <- 0, S <- 1, k <- 0
  size = (nbrLen+1)*sizeof(limb);
  (mod + nbrLen)->x = 0;
  (num + nbrLen)->x = 0;
  memcpy(U, mod, size);
  memcpy(V, num, size);
    // Maximum value of R and S can be up to 2*M, so one more limb is needed.
  memset(R, 0, size);   // R <- 0
  memset(S, 0, size);   // S <- 1
  S[0].x = 1;
  k = 0;
  // R' <- aR + bS, S' <- cR + dS
  a = d = 1;  // R' = R, S' = S.
  b = c = 0;
  len = nbrLen;
  // Find length of U.
  for (lenU = nbrLen - 1; lenU > 0; lenU--)
  {
    if (U[lenU].x != 0)
    {
      break;
    }
  }
  lenU++;
  // Find length of V.
  for (lenV = nbrLen - 1; lenV > 0; lenV--)
  {
    if (V[lenV].x != 0)
    {
      break;
    }
  }
  lenV++;
  //  2. while V > 0 do
  while (lenV > 1 || V[0].x > 0)
  {
    //  3.   if U even then U <- U / 2, S <- 2S
    if ((U[0].x & 1) == 0)
    {     // U is even.
      for (i = 0; i < lenU; i++)
      {  // Loop that divides U by 2.
        U[i].x = ((U[i].x >> 1) | (U[i + 1].x << (BITS_PER_GROUP - 1))) & MAX_VALUE_LIMB;
      }
      if (U[lenU - 1].x == 0)
      {
        lenU--;
      }
      // R' <- aR + bS, S' <- cR + dS
      c *= 2; d *= 2;  // Multiply S by 2.
    }
    //  4.   elsif V even then V <- V / 2, R <- 2R
    else if ((V[0].x & 1) == 0)
    {    // V is even.
      for (i = 0; i < lenV; i++)
      {  // Loop that divides V by 2.
        V[i].x = ((V[i].x >> 1) | (V[i+1].x << (BITS_PER_GROUP-1))) & MAX_VALUE_LIMB;
      }
      if (V[lenV - 1].x == 0)
      {
        lenV--;
      }
      // R' <- aR + bS, S' <- cR + dS
      a *= 2; b *= 2;  // Multiply R by 2.
    }
    //  5.   elsif U >= V  then U <- (U - V) / 2, R <- R + S, S <- 2S
    else
    {
      len = (lenU > lenV ? lenU : lenV);
      for (i = len - 1; i > 0; i--)
      {
        if (U[i].x != V[i].x)
        {
          break;
        }
      }
      if (U[i].x > V[i].x)
      {     // U > V
        lenU = HalveDifference(U, V, len); // U <- (U - V) / 2
        // R' <- aR + bS, S' <- cR + dS
        a += c; b += d;  // R <- R + S
        c *= 2; d *= 2;  // S <- 2S
      }
      //  6.   elsif V >= U then V <- (V - U) / 2, S <- S + R, R <- 2R
      else
      {    // V >= U
        lenV = HalveDifference(V, U, len); // V <- (V - U) / 2
        // R' <- aR + bS, S' <- cR + dS
        c += a; d += b;  // S <- S + R
        a *= 2; b *= 2;  // R <- 2R
      }
    }
    //  7.   k <- k + 1
    k++;
    // Adjust variables.
    if (a >= HALF_INT_RANGE/2 || b >= HALF_INT_RANGE/2 ||
        c >= HALF_INT_RANGE/2 || d >= HALF_INT_RANGE/2)
      {  // Possible overflow of a...d in next cycle, so
       // compute now R and S and reset a, b, c and d.
      AddMult(R, a, b, S, c, d, nbrLen+1);
      // R' <- aR + bS, S' <- cR + dS
      a = d = 1;  // R' = R, S' = S.
      b = c = 0;
    }
  }
  AddMult(R, a, b, S, c, d, nbrLen+1);
  //  8. if R >= M then R <- R - M
  for (i = nbrLen; i > 0; i--)
  {
    if (R[i].x != (mod + i)->x)
    {
      break;
    }
  }
  if ((unsigned int)R[i].x >= (unsigned int)(mod + i)->x)
  {      // R >= M.
    borrow = 0;
    for (i = 0; i <= nbrLen; i++)
    {
      borrow += R[i].x - (mod + i)->x;
      R[i].x = borrow & MAX_VALUE_LIMB;
      borrow >>= BITS_PER_GROUP;
    }
  }
  //  9. R <- M - R
  borrow = 0;
  for (i = 0; i <= nbrLen; i++)
  {
    borrow += (mod + i)->x - R[i].x;
    R[i].x = borrow & MAX_VALUE_LIMB;
    borrow >>= BITS_PER_GROUP;
  }
  R[nbrLen].x = 0;
  // At this moment R = x^(-1)*2^k
  // 10. R <- MonPro(R, R2)
  modmult(R, MontgomeryMultR2, R);
  R[nbrLen].x = 0;
  // At this moment R = x^(-1)*2^(k+m)
  // 11. return MonPro(R, 2^(m-k))
  memset(S, 0, size);
  bitCount = nbrLen*BITS_PER_GROUP - k;
  if (bitCount < 0)
  {
    bitCount += nbrLen*BITS_PER_GROUP;
    S[bitCount / BITS_PER_GROUP].x = 1 << (bitCount % BITS_PER_GROUP);
    modmult(R, S, inv);
  }
  else
  {
    S[bitCount / BITS_PER_GROUP].x = 1 << (bitCount % BITS_PER_GROUP);
    modmult(R, S, inv);
    modmult(inv, MontgomeryMultR2, inv);
  }
}
void BigIntModularDivision(BigInteger *Num, BigInteger *Den, BigInteger *mod, BigInteger *quotient)
{
  NumberLength = mod -> nbrLimbs;
  memcpy(TestNbr, mod -> limbs, NumberLength * sizeof(limb));
  TestNbr[NumberLength].x = 0;
  GetMontgomeryParms(NumberLength);
  BigIntRemainder(Num, mod, &tmpNum);
  BigIntRemainder(Den, mod, &tmpDen);
  if (tmpNum.sign == SIGN_NEGATIVE)
  {
    BigIntAdd(&tmpNum, mod, &tmpNum);
  }
  if (tmpDen.sign == SIGN_NEGATIVE)
  {
    BigIntAdd(&tmpDen, mod, &tmpDen);
  }
  CompressLimbsBigInteger(aux3, &tmpDen);
  modmult(aux3, MontgomeryMultR2, aux3);  // aux3 <- Den in Montgomery notation
  ModInvBigNbr(aux3, aux3, TestNbr, NumberLength); // aux3 <- 1 / Den in Montg notation.
  CompressLimbsBigInteger(aux4, &tmpNum);
  modmult(aux3, aux4, aux3);              // aux3 <- Num / Dev in standard notation.
  UncompressLimbsBigInteger(aux3, quotient);  // Get Num/Den
}

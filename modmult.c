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
#include <stdint.h>
#include "bignbr.h"
char MontgomeryMultNCached;
char TestNbrCached;
limb MontgomeryR1[MAX_LEN];
limb TestNbr[MAX_LEN];
limb MontgomeryMultN[MAX_LEN];
limb MontgomeryMultR1[MAX_LEN];
limb MontgomeryMultR2[MAX_LEN];
int powerOf2Exponent;
static limb aux[MAX_LEN];
static limb aux2[MAX_LEN];
static limb aux3[MAX_LEN];
static limb aux4[MAX_LEN];
static limb aux5[MAX_LEN];
static limb aux6[MAX_LEN];
static limb resultModOdd[MAX_LEN];
static limb resultModPower2[MAX_LEN];
static int NumberLength2;
int NumberLength;
int NumberLengthR1;
#ifdef __EMSCRIPTEN__
int64_t lModularMult;
#endif
mmCback modmultCallback;
static limb U[MAX_LEN];
static limb V[MAX_LEN];
static limb R[MAX_LEN];
static limb S[MAX_LEN];
static limb Ubak[MAX_LEN];
static limb Vbak[MAX_LEN];
static BigInteger tmpDen;
static BigInteger tmpNum;
static BigInteger oddValue;
static BigInteger tmpFact1;
static BigInteger tmpFact2;

int getNbrLimbs(const limb *bigNbr)
{
  const limb *ptrLimb = bigNbr + NumberLength;
  while (ptrLimb > bigNbr)
  {
    ptrLimb--;
    if (ptrLimb->x != 0)
    {
      return (int)(ptrLimb - bigNbr + 1);
    }
  }
  return 1;
}

// Find the inverse of value mod 2^(NumberLength*BITS_PER_GROUP)
void ComputeInversePower2(const limb *value, limb *result, limb *tmp)
{
  int N;
  int x;
  int j;
  limb Cy;
  N = value->x;                // 2 least significant bits of inverse correct.
  x = N;
  x = x * (2 - (N * x));       // 4 least significant bits of inverse correct.
  x = x * (2 - (N * x));       // 8 least significant bits of inverse correct.
  x = x * (2 - (N * x));       // 16 least significant bits of inverse correct.
  x = x * (2 - (N * x));       // 32 least significant bits of inverse correct.
  result->x = x & MAX_VALUE_LIMB;
  for (int currLen = 2; currLen < NumberLength; currLen <<= 1)
  {
    multiply(value, result, tmp, currLen, NULL);    // tmp <- N * x
    Cy.x = 2 - tmp[0].x;
    tmp[0].x = Cy.x & MAX_VALUE_LIMB;
    for (j = 1; j < currLen; j++)
    {
      Cy.x = (Cy.x >> BITS_PER_GROUP) - tmp[j].x;
      tmp[j].x = Cy.x & MAX_VALUE_LIMB;
    }                                                  // tmp <- 2 - N * x
    multiply(result, tmp, result, currLen, NULL);      // tmp <- x * (2 - N * x)
  }
  // Perform last approximation to inverse.
  multiply(value, result, tmp, NumberLength, NULL);    // tmp <- N * x
  Cy.x = 2 - tmp[0].x;
  tmp[0].x = Cy.x & MAX_VALUE_LIMB;
  for (j = 1; j < NumberLength; j++)
  {
    Cy.x = (Cy.x >> BITS_PER_GROUP) - tmp[j].x;
    tmp[j].x = Cy.x & MAX_VALUE_LIMB;
  }                                                    // tmp <- 2 - N * x
  multiply(result, tmp, result, NumberLength, NULL);   // tmp <- x * (2 - N * x)
}

// Let R be a power of 2 of at least len limbs.
// Compute R1 = MontgomeryR1 and N = MontgomeryN using the formulas:
// R1 = R mod M
// N = -M^(-1) mod R

void GetMontgomeryParms(int len)
{
  int j;
  limb Cy;
  MontgomeryMultNCached = NBR_NOT_CACHED;
  TestNbrCached = NBR_NOT_CACHED;
  TestNbr[len].x = 0;
  NumberLength = len;
  NumberLength2 = len + len;
  powerOf2Exponent = 0;    // Indicate not power of 2 in advance.
  if ((NumberLength == 1) && ((TestNbr[0].x & 1) != 0))
  {
    MontgomeryMultR1[0].x = 1;
    MontgomeryMultR2[0].x = 1;
    NumberLengthR1 = 1;
    return;
  }
  // Check whether TestNbr is a power of 2.
  for (j = 0; j < (NumberLength-1); j++)
  {
    if (TestNbr[j].x != 0)
    {
      break;
    }
  }
  if (j == (NumberLength - 1))
  {
    int value = TestNbr[NumberLength - 1].x;
    for (j = 0; j < BITS_PER_GROUP; j++)
    {
      if (value == 1)
      {
        powerOf2Exponent = ((NumberLength - 1)*BITS_PER_GROUP) + j;
        (void)memset(MontgomeryMultR1, 0, NumberLength*sizeof(limb));
        (void)memset(MontgomeryMultR2, 0, NumberLength * sizeof(limb));
        MontgomeryMultR1[0].x = 1;
        MontgomeryMultR2[0].x = 1;
        return;
      }
      value >>= 1;
    }
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
    int x;
    int N;
    N = TestNbr[0].x;            // 2 least significant bits of inverse correct.
    x = N;
    x = x * (2 - (N * x));       // 4 least significant bits of inverse correct.
    x = x * (2 - (N * x));       // 8 least significant bits of inverse correct.
    x = x * (2 - (N * x));       // 16 least significant bits of inverse correct.
    x = x * (2 - (N * x));       // 32 least significant bits of inverse correct.
    MontgomeryMultN[0].x = (-x) & MAX_VALUE_LIMB;    // Change sign
  }
  // Compute MontgomeryMultR1 as 1 in Montgomery notation,
  // this is 2^(NumberLength*BITS_PER_GROUP) % TestNbr.
  j = NumberLength;
  MontgomeryMultR1[j].x = 1;
  do
  {
    j--;
    MontgomeryMultR1[j].x = 0;
  } while (j > 0);
  AdjustModN(MontgomeryMultR1, TestNbr, len);
  MontgomeryMultR1[NumberLength].x = 0;
  (void)memcpy(MontgomeryMultR2, MontgomeryMultR1, (NumberLength+1)*sizeof(limb));
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
    (void)memmove(&MontgomeryMultR2[1], &MontgomeryMultR2[0], NumberLength*sizeof(limb));
    MontgomeryMultR2[0].x = 0;
    AdjustModN(MontgomeryMultR2, TestNbr, len);
  }
  MontgomeryMultNCached = NBR_READY_TO_BE_CACHED;
  TestNbrCached = NBR_READY_TO_BE_CACHED;
}

// Compute Nbr <- Nbr mod Modulus.
// Modulus has NumberLength limbs.
void AdjustModN(limb *Nbr, const limb *Modulus, int nbrLen)
{
#ifdef _USING64BITS_
  int64_t carry;
#else
  int carry;
  double dVal = 1.0 / (double)LIMB_RANGE;
  double dSquareLimb = (double)LIMB_RANGE * (double)LIMB_RANGE;
  double dDelta;
  double dTrialQuotient;
#endif
  int i;
  int TrialQuotient;
  double dNbr;
  double dInvModulus;

  dInvModulus = 1/getMantissa(Modulus+nbrLen, nbrLen);
  dNbr = getMantissa(Nbr + nbrLen + 1, nbrLen + 1) * LIMB_RANGE;
  TrialQuotient = (int)(unsigned int)floor((dNbr * dInvModulus) + 0.5);
  if ((unsigned int)TrialQuotient >= LIMB_RANGE)
  {   // Maximum value for limb.
    TrialQuotient = MAX_VALUE_LIMB;
  }
#ifndef _USING64BITS_
  dTrialQuotient = (double)TrialQuotient;
  dDelta = 0;
#endif
  // Compute Nbr <- Nbr - TrialQuotient * Modulus
  carry = 0;
  for (i = 0; i <= nbrLen; i++)
  {
#ifdef _USING64BITS_
    carry += (int64_t)Nbr[i].x - (Modulus[i].x * (int64_t)TrialQuotient);
    Nbr[i].x = (int)carry & MAX_VALUE_LIMB;
    carry >>= BITS_PER_GROUP;
#else
    int low = (Nbr[i].x - (Modulus[i].x * TrialQuotient) + carry) & MAX_INT_NBR;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    double dAccumulator = (double)Nbr[i].x - ((double)Modulus[i].x * dTrialQuotient) +
      (double)carry + dDelta;
    dDelta = 0;
    if (dAccumulator < 0)
    {
      dAccumulator += dSquareLimb;
      dDelta = -(double)LIMB_RANGE;
    }
    if (low < HALF_INT_RANGE)
    {
      carry = (int)floor((dAccumulator + (double)FOURTH_INT_RANGE)*dVal);
    }
    else
    {
      carry = (int)floor((dAccumulator - (double)FOURTH_INT_RANGE)*dVal);
    }
    Nbr[i].x = low;
#endif
  }
  Nbr[i].x = (int)carry & MAX_INT_NBR;
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

void AddBigNbrModN(const limb *Nbr1, const limb *Nbr2, limb *Sum, const limb *mod, int nbrLen)
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

  if ((carry < LIMB_RANGE) && (borrow < 0))
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

void SubtBigNbrModN(const limb *Nbr1, const limb *Nbr2, limb *Diff, const limb *mod, int nbrLen)
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

void smallmodmult(int factor1, int factor2, limb *product, int mod)
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
    int remainder = (factor1 * factor2) - (quotient * mod);
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

#ifdef _USING64BITS_
static void MontgomeryMult2(const limb *pNbr1, const limb *pNbr2, limb *pProd)
{
  uint64_t Pr;
  uint32_t Prod0 = 0;
  uint32_t Prod1 = 0;
  uint32_t TestNbr0 = TestNbr[0].x;
  uint32_t TestNbr1 = TestNbr[1].x;
  uint32_t Nbr2_0 = pNbr2->x;
  uint32_t Nbr2_1 = (pNbr2 + 1)->x;
  for (int i = 0; i<2; i++)
  {
    uint32_t Nbr;
    uint32_t MontDig;
    Nbr = (pNbr1 + i)->x;
    Pr = (Nbr * (uint64_t)Nbr2_0) + Prod0;
    MontDig = ((uint32_t)Pr * MontgomeryMultN[0].x) & MAX_INT_NBR;
    Pr = ((((uint64_t)MontDig * TestNbr0) + Pr) >> BITS_PER_GROUP) +
      ((uint64_t)MontDig * TestNbr1) + ((uint64_t)Nbr * Nbr2_1) + Prod1;
    Prod0 = Pr & MAX_INT_NBR;
    Prod1 = (uint32_t)(Pr >> BITS_PER_GROUP);
  }
  if ((Pr >= ((uint64_t)(TestNbr1 + 1) << BITS_PER_GROUP)) ||
     ((Prod1 == TestNbr1) && (Prod0 >= TestNbr0)))
  {
    int32_t borrow = (int32_t)Prod0 - (int32_t)TestNbr0;
    Prod0 = borrow & MAX_INT_NBR;
    Prod1 = ((borrow >> BITS_PER_GROUP) + (int32_t)Prod1 - (int32_t)TestNbr1) & MAX_INT_NBR;
  }
  pProd->x = Prod0;
  (pProd + 1)->x = Prod1;
}

static void MontgomeryMult3(const limb *pNbr1, const limb *pNbr2, limb *pProd)
{
  uint64_t Pr;
  uint32_t Prod0 = 0;
  uint32_t Prod1 = 0;
  uint32_t Prod2 = 0;
  uint32_t TestNbr0 = TestNbr[0].x;
  uint32_t TestNbr1 = TestNbr[1].x;
  uint32_t TestNbr2 = TestNbr[2].x;
  uint32_t Nbr2_0 = pNbr2->x;
  uint32_t Nbr2_1 = (pNbr2 + 1)->x;
  uint32_t Nbr2_2 = (pNbr2 + 2)->x;
  for (int i = 0; i<3; i++)
  {
    uint32_t MontDig;
    uint32_t Nbr;
    Pr = (Nbr = (pNbr1 + i)->x) * (uint64_t)Nbr2_0 + Prod0;
    MontDig = ((uint32_t)Pr * MontgomeryMultN[0].x) & MAX_INT_NBR;
    Prod0 = (Pr = ((((uint64_t)MontDig * TestNbr0) + Pr) >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr1 + (uint64_t)Nbr * Nbr2_1 + Prod1) & MAX_INT_NBR;
    Prod1 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr2 + (uint64_t)Nbr * Nbr2_2 + Prod2) & MAX_INT_NBR;
    Prod2 = (uint32_t)(Pr >> BITS_PER_GROUP);
  }
  if ((Pr >= ((uint64_t)(TestNbr2 + 1) << BITS_PER_GROUP))
    || ((Prod2 == TestNbr2) && ((Prod1 > TestNbr1) || ((Prod1 == TestNbr1) && (Prod0 >= TestNbr0)))))
  {
    int32_t borrow;
    Prod0 = (borrow = (int32_t)Prod0 - (int32_t)TestNbr0) & MAX_INT_NBR;
    Prod1 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod1 - (int32_t)TestNbr1) & MAX_INT_NBR;
    Prod2 = ((borrow >> BITS_PER_GROUP) + (int32_t)Prod2 - (int32_t)TestNbr2) & MAX_INT_NBR;
  }
  pProd->x = Prod0;
  (pProd + 1)->x = Prod1;
  (pProd + 2)->x = Prod2;
}

static void MontgomeryMult4(const limb *pNbr1, const limb *pNbr2, limb *pProd)
{
  uint64_t Pr;
  uint32_t Prod0 = 0;
  uint32_t Prod1 = 0;
  uint32_t Prod2 = 0;
  uint32_t Prod3 = 0;
  uint32_t TestNbr0 = TestNbr[0].x;
  uint32_t TestNbr1 = TestNbr[1].x;
  uint32_t TestNbr2 = TestNbr[2].x;
  uint32_t TestNbr3 = TestNbr[3].x;
  uint32_t Nbr2_0 = pNbr2->x;
  uint32_t Nbr2_1 = (pNbr2 + 1)->x;
  uint32_t Nbr2_2 = (pNbr2 + 2)->x;
  uint32_t Nbr2_3 = (pNbr2 + 3)->x;
  for (int i = 0; i<4; i++)
  {
    uint32_t Nbr;
    uint32_t MontDig;
    Pr = (Nbr = (pNbr1 + i)->x) * (uint64_t)Nbr2_0 + Prod0;
    MontDig = ((uint32_t)Pr * MontgomeryMultN[0].x) & MAX_INT_NBR;
    Prod0 = (Pr = ((((uint64_t)MontDig * TestNbr0) + Pr) >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr1 + (uint64_t)Nbr * Nbr2_1 + Prod1) & MAX_INT_NBR;
    Prod1 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr2 + (uint64_t)Nbr * Nbr2_2 + Prod2) & MAX_INT_NBR;
    Prod2 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr3 + (uint64_t)Nbr * Nbr2_3 + Prod3) & MAX_INT_NBR;
    Prod3 = (uint32_t)(Pr >> BITS_PER_GROUP);
  }
  if (Pr >= ((uint64_t)(TestNbr3 + 1) << BITS_PER_GROUP)
    || (Prod3 == TestNbr3
      && (Prod2 > TestNbr2
        || (Prod2 == TestNbr2
          && (Prod1 > TestNbr1 || (Prod1 == TestNbr1 && (Prod0 >= TestNbr0)))))))
  {
    int32_t borrow;
    Prod0 = (borrow = (int32_t)Prod0 - (int32_t)TestNbr0) & MAX_INT_NBR;
    Prod1 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod1 - (int32_t)TestNbr1) & MAX_INT_NBR;
    Prod2 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod2 - (int32_t)TestNbr2) & MAX_INT_NBR;
    Prod3 = ((borrow >> BITS_PER_GROUP) + (int32_t)Prod3 - (int32_t)TestNbr3) & MAX_INT_NBR;
  }
  pProd->x = Prod0;
  (pProd + 1)->x = Prod1;
  (pProd + 2)->x = Prod2;
  (pProd + 3)->x = Prod3;
}

static void MontgomeryMult5(const limb *pNbr1, const limb *pNbr2, limb *pProd)
{
  uint64_t Pr;
  uint32_t Prod0 = 0;
  uint32_t Prod1 = 0;
  uint32_t Prod2 = 0;
  uint32_t Prod3 = 0;
  uint32_t Prod4 = 0;
  uint32_t TestNbr0 = TestNbr[0].x;
  uint32_t TestNbr1 = TestNbr[1].x;
  uint32_t TestNbr2 = TestNbr[2].x;
  uint32_t TestNbr3 = TestNbr[3].x;
  uint32_t TestNbr4 = TestNbr[4].x;
  uint32_t Nbr2_0 = pNbr2->x;
  uint32_t Nbr2_1 = (pNbr2 + 1)->x;
  uint32_t Nbr2_2 = (pNbr2 + 2)->x;
  uint32_t Nbr2_3 = (pNbr2 + 3)->x;
  uint32_t Nbr2_4 = (pNbr2 + 4)->x;
  for (int i = 0; i<5; i++)
  {
    uint32_t Nbr;
    uint32_t MontDig;
    Pr = (Nbr = (pNbr1 + i)->x) * (uint64_t)Nbr2_0 + Prod0;
    MontDig = ((uint32_t)Pr * MontgomeryMultN[0].x) & MAX_INT_NBR;
    Prod0 = (Pr = (((uint64_t)MontDig * TestNbr0 + Pr) >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr1 + (uint64_t)Nbr * Nbr2_1 + Prod1) & MAX_INT_NBR;
    Prod1 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr2 + (uint64_t)Nbr * Nbr2_2 + Prod2) & MAX_INT_NBR;
    Prod2 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr3 + (uint64_t)Nbr * Nbr2_3 + Prod3) & MAX_INT_NBR;
    Prod3 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr4 + (uint64_t)Nbr * Nbr2_4 + Prod4) & MAX_INT_NBR;
    Prod4 = (uint32_t)(Pr >> BITS_PER_GROUP);
  }
  if (Pr >= ((uint64_t)(TestNbr4 + 1) << BITS_PER_GROUP)
    || (Prod4 == TestNbr4
      && (Prod3 > TestNbr3
        || (Prod3 == TestNbr3
          && (Prod2 > TestNbr2
            || (Prod2 == TestNbr2
              && (Prod1 > TestNbr1 || (Prod1 == TestNbr1 && (Prod0 >= TestNbr0)))))))))
  {
    int32_t borrow;
    Prod0 = (borrow = (int32_t)Prod0 - (int32_t)TestNbr0) & MAX_INT_NBR;
    Prod1 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod1 - (int32_t)TestNbr1) & MAX_INT_NBR;
    Prod2 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod2 - (int32_t)TestNbr2) & MAX_INT_NBR;
    Prod3 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod3 - (int32_t)TestNbr3) & MAX_INT_NBR;
    Prod4 = ((borrow >> BITS_PER_GROUP) + (int32_t)Prod4 - (int32_t)TestNbr4) & MAX_INT_NBR;
  }
  pProd->x = Prod0;
  (pProd + 1)->x = Prod1;
  (pProd + 2)->x = Prod2;
  (pProd + 3)->x = Prod3;
  (pProd + 4)->x = Prod4;
}

static void MontgomeryMult6(const limb *pNbr1, const limb *pNbr2, limb *pProd)
{
  uint64_t Pr;
  uint32_t Prod0 = 0;
  uint32_t Prod1 = 0;
  uint32_t Prod2 = 0;
  uint32_t Prod3 = 0;
  uint32_t Prod4 = 0;
  uint32_t Prod5 = 0;
  uint32_t TestNbr0 = TestNbr[0].x;
  uint32_t TestNbr1 = TestNbr[1].x;
  uint32_t TestNbr2 = TestNbr[2].x;
  uint32_t TestNbr3 = TestNbr[3].x;
  uint32_t TestNbr4 = TestNbr[4].x;
  uint32_t TestNbr5 = TestNbr[5].x;
  uint32_t Nbr2_0 = pNbr2->x;
  uint32_t Nbr2_1 = (pNbr2 + 1)->x;
  uint32_t Nbr2_2 = (pNbr2 + 2)->x;
  uint32_t Nbr2_3 = (pNbr2 + 3)->x;
  uint32_t Nbr2_4 = (pNbr2 + 4)->x;
  uint32_t Nbr2_5 = (pNbr2 + 5)->x;
  for (int i = 0; i<6; i++)
  {
    uint32_t Nbr;
    uint32_t MontDig;
    Pr = (Nbr = (pNbr1 + i)->x) * (uint64_t)Nbr2_0 + Prod0;
    MontDig = ((uint32_t)Pr * MontgomeryMultN[0].x) & MAX_INT_NBR;
    Prod0 = (Pr = (((uint64_t)MontDig * TestNbr0 + Pr) >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr1 + (uint64_t)Nbr * Nbr2_1 + Prod1) & MAX_INT_NBR;
    Prod1 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr2 + (uint64_t)Nbr * Nbr2_2 + Prod2) & MAX_INT_NBR;
    Prod2 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr3 + (uint64_t)Nbr * Nbr2_3 + Prod3) & MAX_INT_NBR;
    Prod3 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr4 + (uint64_t)Nbr * Nbr2_4 + Prod4) & MAX_INT_NBR;
    Prod4 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr5 + (uint64_t)Nbr * Nbr2_5 + Prod5) & MAX_INT_NBR;
    Prod5 = (uint32_t)(Pr >> BITS_PER_GROUP);
  }
  if (Pr >= ((uint64_t)(TestNbr5 + 1) << BITS_PER_GROUP)
    || (Prod5 == TestNbr5
      && (Prod4 > TestNbr4
        || (Prod4 == TestNbr4
          && (Prod3 > TestNbr3
            || (Prod3 == TestNbr3
              && (Prod2 > TestNbr2
                || (Prod2 == TestNbr2
                  && (Prod1 > TestNbr1
                    || (Prod1 == TestNbr1
                      && (Prod0 >= TestNbr0)))))))))))
  {
    int32_t borrow;
    Prod0 = (borrow = (int32_t)Prod0 - (int32_t)TestNbr0) & MAX_INT_NBR;
    Prod1 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod1 - (int32_t)TestNbr1) & MAX_INT_NBR;
    Prod2 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod2 - (int32_t)TestNbr2) & MAX_INT_NBR;
    Prod3 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod3 - (int32_t)TestNbr3) & MAX_INT_NBR;
    Prod4 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod4 - (int32_t)TestNbr4) & MAX_INT_NBR;
    Prod5 = ((borrow >> BITS_PER_GROUP) + (int32_t)Prod5 - (int32_t)TestNbr5) & MAX_INT_NBR;
  }
  pProd->x = Prod0;
  (pProd + 1)->x = Prod1;
  (pProd + 2)->x = Prod2;
  (pProd + 3)->x = Prod3;
  (pProd + 4)->x = Prod4;
  (pProd + 5)->x = Prod5;
}

static void MontgomeryMult7(const limb *pNbr1, const limb *pNbr2, limb *pProd)
{
  uint64_t Pr;
  uint32_t Prod0 = 0;
  uint32_t Prod1 = 0;
  uint32_t Prod2 = 0;
  uint32_t Prod3 = 0;
  uint32_t Prod4 = 0;
  uint32_t Prod5 = 0;
  uint32_t Prod6 = 0;
  uint32_t TestNbr0 = TestNbr[0].x;
  uint32_t TestNbr1 = TestNbr[1].x;
  uint32_t TestNbr2 = TestNbr[2].x;
  uint32_t TestNbr3 = TestNbr[3].x;
  uint32_t TestNbr4 = TestNbr[4].x;
  uint32_t TestNbr5 = TestNbr[5].x;
  uint32_t TestNbr6 = TestNbr[6].x;
  int Nbr2_0 = pNbr2->x;
  int Nbr2_1 = (pNbr2 + 1)->x;
  int Nbr2_2 = (pNbr2 + 2)->x;
  int Nbr2_3 = (pNbr2 + 3)->x;
  int Nbr2_4 = (pNbr2 + 4)->x;
  int Nbr2_5 = (pNbr2 + 5)->x;
  int Nbr2_6 = (pNbr2 + 6)->x;
  for (int i = 0; i<7; i++)
  {
    uint32_t Nbr;
    uint32_t MontDig;
    Pr = (Nbr = (pNbr1 + i)->x) * (uint64_t)Nbr2_0 + Prod0;
    MontDig = ((uint32_t)Pr * MontgomeryMultN[0].x) & MAX_INT_NBR;
    Prod0 = (Pr = (((uint64_t)MontDig * TestNbr0 + Pr) >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr1 + (uint64_t)Nbr * Nbr2_1 + Prod1) & MAX_INT_NBR;
    Prod1 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr2 + (uint64_t)Nbr * Nbr2_2 + Prod2) & MAX_INT_NBR;
    Prod2 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr3 + (uint64_t)Nbr * Nbr2_3 + Prod3) & MAX_INT_NBR;
    Prod3 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr4 + (uint64_t)Nbr * Nbr2_4 + Prod4) & MAX_INT_NBR;
    Prod4 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr5 + (uint64_t)Nbr * Nbr2_5 + Prod5) & MAX_INT_NBR;
    Prod5 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr6 + (uint64_t)Nbr * Nbr2_6 + Prod6) & MAX_INT_NBR;
    Prod6 = (uint32_t)(Pr >> BITS_PER_GROUP);
  }
  if (Pr >= ((uint64_t)(TestNbr6 + 1) << BITS_PER_GROUP)
    || (Prod6 == TestNbr6
      && (Prod5 > TestNbr5
        || (Prod5 == TestNbr5
          && (Prod4 > TestNbr4
            || (Prod4 == TestNbr4
              && (Prod3 > TestNbr3
                || (Prod3 == TestNbr3
                  && (Prod2 > TestNbr2
                    || (Prod2 == TestNbr2
                      && (Prod1 > TestNbr1
                        || (Prod1 == TestNbr1
                          && (Prod0 >= TestNbr0)))))))))))))
  {
    int32_t borrow;
    Prod0 = (borrow = (int32_t)Prod0 - TestNbr0) & MAX_INT_NBR;
    Prod1 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod1 - (int32_t)TestNbr1) & MAX_INT_NBR;
    Prod2 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod2 - (int32_t)TestNbr2) & MAX_INT_NBR;
    Prod3 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod3 - (int32_t)TestNbr3) & MAX_INT_NBR;
    Prod4 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod4 - (int32_t)TestNbr4) & MAX_INT_NBR;
    Prod5 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod5 - (int32_t)TestNbr5) & MAX_INT_NBR;
    Prod6 = ((borrow >> BITS_PER_GROUP) + (int32_t)Prod6 - TestNbr6) & MAX_INT_NBR;
  }
  pProd->x = Prod0;
  (pProd + 1)->x = Prod1;
  (pProd + 2)->x = Prod2;
  (pProd + 3)->x = Prod3;
  (pProd + 4)->x = Prod4;
  (pProd + 5)->x = Prod5;
  (pProd + 6)->x = Prod6;
}

static void MontgomeryMult8(const limb *pNbr1, const limb *pNbr2, limb *pProd)
{
  uint64_t Pr;
  uint32_t Prod0 = 0;
  uint32_t Prod1 = 0;
  uint32_t Prod2 = 0;
  uint32_t Prod3 = 0;
  uint32_t Prod4 = 0;
  uint32_t Prod5 = 0;
  uint32_t Prod6 = 0;
  uint32_t Prod7 = 0;
  uint32_t TestNbr0 = TestNbr[0].x;
  uint32_t TestNbr1 = TestNbr[1].x;
  uint32_t TestNbr2 = TestNbr[2].x;
  uint32_t TestNbr3 = TestNbr[3].x;
  uint32_t TestNbr4 = TestNbr[4].x;
  uint32_t TestNbr5 = TestNbr[5].x;
  uint32_t TestNbr6 = TestNbr[6].x;
  uint32_t TestNbr7 = TestNbr[7].x;
  uint32_t Nbr2_0 = pNbr2->x;
  uint32_t Nbr2_1 = (pNbr2 + 1)->x;
  uint32_t Nbr2_2 = (pNbr2 + 2)->x;
  uint32_t Nbr2_3 = (pNbr2 + 3)->x;
  uint32_t Nbr2_4 = (pNbr2 + 4)->x;
  uint32_t Nbr2_5 = (pNbr2 + 5)->x;
  uint32_t Nbr2_6 = (pNbr2 + 6)->x;
  uint32_t Nbr2_7 = (pNbr2 + 7)->x;
  for (int i = 0; i<8; i++)
  {
    uint32_t Nbr;
    uint32_t MontDig;
    Pr = (Nbr = (pNbr1 + i)->x) * (uint64_t)Nbr2_0 + Prod0;
    MontDig = ((uint32_t)Pr * MontgomeryMultN[0].x) & MAX_INT_NBR;
    Prod0 = (Pr = (((uint64_t)MontDig * TestNbr0 + Pr) >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr1 + (uint64_t)Nbr * Nbr2_1 + Prod1) & MAX_INT_NBR;
    Prod1 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr2 + (uint64_t)Nbr * Nbr2_2 + Prod2) & MAX_INT_NBR;
    Prod2 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr3 + (uint64_t)Nbr * Nbr2_3 + Prod3) & MAX_INT_NBR;
    Prod3 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr4 + (uint64_t)Nbr * Nbr2_4 + Prod4) & MAX_INT_NBR;
    Prod4 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr5 + (uint64_t)Nbr * Nbr2_5 + Prod5) & MAX_INT_NBR;
    Prod5 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr6 + (uint64_t)Nbr * Nbr2_6 + Prod6) & MAX_INT_NBR;
    Prod6 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr7 + (uint64_t)Nbr * Nbr2_7 + Prod7) & MAX_INT_NBR;
    Prod7 = (uint32_t)(Pr >> BITS_PER_GROUP);
  }
  if (Pr >= ((uint64_t)(TestNbr7 + 1) << BITS_PER_GROUP)
    || (Prod7 == TestNbr7
      && (Prod6 > TestNbr6
        || (Prod6 == TestNbr6
          && (Prod5 > TestNbr5
            || (Prod5 == TestNbr5
              && (Prod4 > TestNbr4
                || (Prod4 == TestNbr4
                  && (Prod3 > TestNbr3
                    || (Prod3 == TestNbr3
                      && (Prod2 > TestNbr2
                        || (Prod2 == TestNbr2
                          && (Prod1 > TestNbr1
                            || (Prod1 == TestNbr1
                              && (Prod0 >= TestNbr0)))))))))))))))
  {
    int32_t borrow;
    Prod0 = (borrow = (int32_t)Prod0 - (int32_t)TestNbr0) & MAX_INT_NBR;
    Prod1 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod1 - (int32_t)TestNbr1) & MAX_INT_NBR;
    Prod2 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod2 - (int32_t)TestNbr2) & MAX_INT_NBR;
    Prod3 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod3 - (int32_t)TestNbr3) & MAX_INT_NBR;
    Prod4 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod4 - (int32_t)TestNbr4) & MAX_INT_NBR;
    Prod5 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod5 - (int32_t)TestNbr5) & MAX_INT_NBR;
    Prod6 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod6 - (int32_t)TestNbr6) & MAX_INT_NBR;
    Prod7 = ((borrow >> BITS_PER_GROUP) + (int32_t)Prod7 - (int32_t)TestNbr7) & MAX_INT_NBR;
  }
  pProd->x = (int)Prod0;
  (pProd + 1)->x = (int)Prod1;
  (pProd + 2)->x = (int)Prod2;
  (pProd + 3)->x = (int)Prod3;
  (pProd + 4)->x = (int)Prod4;
  (pProd + 5)->x = (int)Prod5;
  (pProd + 6)->x = (int)Prod6;
  (pProd + 7)->x = (int)Prod7;
}
static void MontgomeryMult9(const limb *pNbr1, const limb *pNbr2, limb *pProd)
{
  uint64_t Pr;
  uint32_t Prod0 = 0;
  uint32_t Prod1 = 0;
  uint32_t Prod2 = 0;
  uint32_t Prod3 = 0;
  uint32_t Prod4 = 0;
  uint32_t Prod5 = 0;
  uint32_t Prod6 = 0;
  uint32_t Prod7 = 0;
  uint32_t Prod8 = 0;
  uint32_t TestNbr0 = TestNbr[0].x;
  uint32_t TestNbr1 = TestNbr[1].x;
  uint32_t TestNbr2 = TestNbr[2].x;
  uint32_t TestNbr3 = TestNbr[3].x;
  uint32_t TestNbr4 = TestNbr[4].x;
  uint32_t TestNbr5 = TestNbr[5].x;
  uint32_t TestNbr6 = TestNbr[6].x;
  uint32_t TestNbr7 = TestNbr[7].x;
  uint32_t TestNbr8 = TestNbr[8].x;
  uint32_t Nbr2_0 = pNbr2->x;
  uint32_t Nbr2_1 = (pNbr2 + 1)->x;
  uint32_t Nbr2_2 = (pNbr2 + 2)->x;
  uint32_t Nbr2_3 = (pNbr2 + 3)->x;
  uint32_t Nbr2_4 = (pNbr2 + 4)->x;
  uint32_t Nbr2_5 = (pNbr2 + 5)->x;
  uint32_t Nbr2_6 = (pNbr2 + 6)->x;
  uint32_t Nbr2_7 = (pNbr2 + 7)->x;
  uint32_t Nbr2_8 = (pNbr2 + 8)->x;
  for (int i = 0; i<9; i++)
  {
    uint32_t Nbr;
    uint32_t MontDig;
    Pr = (Nbr = (pNbr1 + i)->x) * (uint64_t)Nbr2_0 + Prod0;
    MontDig = ((uint32_t)Pr * MontgomeryMultN[0].x) & MAX_INT_NBR;
    Prod0 = (Pr = (((uint64_t)MontDig * TestNbr0 + Pr) >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr1 + (uint64_t)Nbr * Nbr2_1 + Prod1) & MAX_INT_NBR;
    Prod1 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr2 + (uint64_t)Nbr * Nbr2_2 + Prod2) & MAX_INT_NBR;
    Prod2 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr3 + (uint64_t)Nbr * Nbr2_3 + Prod3) & MAX_INT_NBR;
    Prod3 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr4 + (uint64_t)Nbr * Nbr2_4 + Prod4) & MAX_INT_NBR;
    Prod4 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr5 + (uint64_t)Nbr * Nbr2_5 + Prod5) & MAX_INT_NBR;
    Prod5 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr6 + (uint64_t)Nbr * Nbr2_6 + Prod6) & MAX_INT_NBR;
    Prod6 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr7 + (uint64_t)Nbr * Nbr2_7 + Prod7) & MAX_INT_NBR;
    Prod7 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr8 + (uint64_t)Nbr * Nbr2_8 + Prod8) & MAX_INT_NBR;
    Prod8 = (uint32_t)(Pr >> BITS_PER_GROUP);
  }
  if (Pr >= ((uint64_t)(TestNbr8 + 1) << BITS_PER_GROUP)
    || (Prod8 == TestNbr8
      && (Prod7 > TestNbr7
        || (Prod7 == TestNbr7
          && (Prod6 > TestNbr6
            || (Prod6 == TestNbr6
              && (Prod5 > TestNbr5
                || (Prod5 == TestNbr5
                  && (Prod4 > TestNbr4
                    || (Prod4 == TestNbr4
                      && (Prod3 > TestNbr3
                        || (Prod3 == TestNbr3
                          && (Prod2 > TestNbr2
                            || (Prod2 == TestNbr2
                              && (Prod1 > TestNbr1
                                || (Prod1 == TestNbr1
                                  && (Prod0 >= TestNbr0)))))))))))))))))
  {
    int32_t borrow;
    Prod0 = (borrow = (int32_t)Prod0 - (int32_t)TestNbr0) & MAX_INT_NBR;
    Prod1 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod1 - (int32_t)TestNbr1) & MAX_INT_NBR;
    Prod2 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod2 - (int32_t)TestNbr2) & MAX_INT_NBR;
    Prod3 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod3 - (int32_t)TestNbr3) & MAX_INT_NBR;
    Prod4 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod4 - (int32_t)TestNbr4) & MAX_INT_NBR;
    Prod5 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod5 - (int32_t)TestNbr5) & MAX_INT_NBR;
    Prod6 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod6 - (int32_t)TestNbr6) & MAX_INT_NBR;
    Prod7 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod7 - (int32_t)TestNbr7) & MAX_INT_NBR;
    Prod8 = ((borrow >> BITS_PER_GROUP) + (int32_t)Prod8 - (int32_t)TestNbr8) & MAX_INT_NBR;
  }
  pProd->x = Prod0;
  (pProd + 1)->x = Prod1;
  (pProd + 2)->x = Prod2;
  (pProd + 3)->x = Prod3;
  (pProd + 4)->x = Prod4;
  (pProd + 5)->x = Prod5;
  (pProd + 6)->x = Prod6;
  (pProd + 7)->x = Prod7;
  (pProd + 8)->x = Prod8;
}

static void MontgomeryMult10(const limb *pNbr1, const limb *pNbr2, limb *pProd)
{
  uint64_t Pr;
  uint32_t Prod0 = 0;
  uint32_t Prod1 = 0;
  uint32_t Prod2 = 0;
  uint32_t Prod3 = 0;
  uint32_t Prod4 = 0;
  uint32_t Prod5 = 0;
  uint32_t Prod6 = 0;
  uint32_t Prod7 = 0;
  uint32_t Prod8 = 0;
  uint32_t Prod9 = 0;
  uint32_t TestNbr0 = TestNbr[0].x;
  uint32_t TestNbr1 = TestNbr[1].x;
  uint32_t TestNbr2 = TestNbr[2].x;
  uint32_t TestNbr3 = TestNbr[3].x;
  uint32_t TestNbr4 = TestNbr[4].x;
  uint32_t TestNbr5 = TestNbr[5].x;
  uint32_t TestNbr6 = TestNbr[6].x;
  uint32_t TestNbr7 = TestNbr[7].x;
  uint32_t TestNbr8 = TestNbr[8].x;
  uint32_t TestNbr9 = TestNbr[9].x;
  uint32_t Nbr2_0 = pNbr2->x;
  uint32_t Nbr2_1 = (pNbr2 + 1)->x;
  uint32_t Nbr2_2 = (pNbr2 + 2)->x;
  uint32_t Nbr2_3 = (pNbr2 + 3)->x;
  uint32_t Nbr2_4 = (pNbr2 + 4)->x;
  uint32_t Nbr2_5 = (pNbr2 + 5)->x;
  uint32_t Nbr2_6 = (pNbr2 + 6)->x;
  uint32_t Nbr2_7 = (pNbr2 + 7)->x;
  uint32_t Nbr2_8 = (pNbr2 + 8)->x;
  uint32_t Nbr2_9 = (pNbr2 + 9)->x;
  for (int i = 0; i<10; i++)
  {
    uint32_t Nbr;
    uint32_t MontDig;
    Pr = (Nbr = (pNbr1 + i)->x) * (uint64_t)Nbr2_0 + Prod0;
    MontDig = ((uint32_t)Pr * MontgomeryMultN[0].x) & MAX_INT_NBR;
    Prod0 = (Pr = (((uint64_t)MontDig * TestNbr0 + Pr) >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr1 + (uint64_t)Nbr * Nbr2_1 + Prod1) & MAX_INT_NBR;
    Prod1 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr2 + (uint64_t)Nbr * Nbr2_2 + Prod2) & MAX_INT_NBR;
    Prod2 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr3 + (uint64_t)Nbr * Nbr2_3 + Prod3) & MAX_INT_NBR;
    Prod3 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr4 + (uint64_t)Nbr * Nbr2_4 + Prod4) & MAX_INT_NBR;
    Prod4 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr5 + (uint64_t)Nbr * Nbr2_5 + Prod5) & MAX_INT_NBR;
    Prod5 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr6 + (uint64_t)Nbr * Nbr2_6 + Prod6) & MAX_INT_NBR;
    Prod6 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr7 + (uint64_t)Nbr * Nbr2_7 + Prod7) & MAX_INT_NBR;
    Prod7 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr8 + (uint64_t)Nbr * Nbr2_8 + Prod8) & MAX_INT_NBR;
    Prod8 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr9 + (uint64_t)Nbr * Nbr2_9 + Prod9) & MAX_INT_NBR;
    Prod9 = (uint32_t)(Pr >> BITS_PER_GROUP);
  }
  if (Pr >= ((uint64_t)(TestNbr9 + 1) << BITS_PER_GROUP)
    || (Prod9 == TestNbr9
      && (Prod8 > TestNbr8
        || (Prod8 == TestNbr8
          && (Prod7 > TestNbr7
            || (Prod7 == TestNbr7
              && (Prod6 > TestNbr6
                || (Prod6 == TestNbr6
                  && (Prod5 > TestNbr5
                    || (Prod5 == TestNbr5
                      && (Prod4 > TestNbr4
                        || (Prod4 == TestNbr4
                          && (Prod3 > TestNbr3
                            || (Prod3 == TestNbr3
                              && (Prod2 > TestNbr2
                                || (Prod2 == TestNbr2
                                  && (Prod1 > TestNbr1
                                    || (Prod1 == TestNbr1
                                      && (Prod0 >= TestNbr0)))))))))))))))))))
  {
    int32_t borrow;
    Prod0 = (borrow = (int32_t)Prod0 - (int32_t)TestNbr0) & MAX_INT_NBR;
    Prod1 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod1 - (int32_t)TestNbr1) & MAX_INT_NBR;
    Prod2 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod2 - (int32_t)TestNbr2) & MAX_INT_NBR;
    Prod3 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod3 - (int32_t)TestNbr3) & MAX_INT_NBR;
    Prod4 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod4 - (int32_t)TestNbr4) & MAX_INT_NBR;
    Prod5 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod5 - (int32_t)TestNbr5) & MAX_INT_NBR;
    Prod6 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod6 - (int32_t)TestNbr6) & MAX_INT_NBR;
    Prod7 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod7 - (int32_t)TestNbr7) & MAX_INT_NBR;
    Prod8 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod8 - (int32_t)TestNbr8) & MAX_INT_NBR;
    Prod9 = ((borrow >> BITS_PER_GROUP) + (int32_t)Prod9 - (int32_t)TestNbr9) & MAX_INT_NBR;
  }
  pProd->x = Prod0;
  (pProd + 1)->x = Prod1;
  (pProd + 2)->x = Prod2;
  (pProd + 3)->x = Prod3;
  (pProd + 4)->x = Prod4;
  (pProd + 5)->x = Prod5;
  (pProd + 6)->x = Prod6;
  (pProd + 7)->x = Prod7;
  (pProd + 8)->x = Prod8;
  (pProd + 9)->x = Prod9;
}

static void MontgomeryMult11(const limb *pNbr1, const limb *pNbr2, limb *pProd)
{
  uint64_t Pr;
  uint32_t Prod0 = 0;
  uint32_t Prod1 = 0;
  uint32_t Prod2 = 0;
  uint32_t Prod3 = 0;
  uint32_t Prod4 = 0;
  uint32_t Prod5 = 0;
  uint32_t Prod6 = 0;
  uint32_t Prod7 = 0;
  uint32_t Prod8 = 0;
  uint32_t Prod9 = 0;
  uint32_t Prod10 = 0;
  uint32_t TestNbr0 = TestNbr[0].x;
  uint32_t TestNbr1 = TestNbr[1].x;
  uint32_t TestNbr2 = TestNbr[2].x;
  uint32_t TestNbr3 = TestNbr[3].x;
  uint32_t TestNbr4 = TestNbr[4].x;
  uint32_t TestNbr5 = TestNbr[5].x;
  uint32_t TestNbr6 = TestNbr[6].x;
  uint32_t TestNbr7 = TestNbr[7].x;
  uint32_t TestNbr8 = TestNbr[8].x;
  uint32_t TestNbr9 = TestNbr[9].x;
  uint32_t TestNbr10 = TestNbr[10].x;
  uint32_t Nbr2_0 = pNbr2->x;
  uint32_t Nbr2_1 = (pNbr2 + 1)->x;
  uint32_t Nbr2_2 = (pNbr2 + 2)->x;
  uint32_t Nbr2_3 = (pNbr2 + 3)->x;
  uint32_t Nbr2_4 = (pNbr2 + 4)->x;
  uint32_t Nbr2_5 = (pNbr2 + 5)->x;
  uint32_t Nbr2_6 = (pNbr2 + 6)->x;
  uint32_t Nbr2_7 = (pNbr2 + 7)->x;
  uint32_t Nbr2_8 = (pNbr2 + 8)->x;
  uint32_t Nbr2_9 = (pNbr2 + 9)->x;
  uint32_t Nbr2_10 = (pNbr2 + 10)->x;
  for (int i = 0; i<11; i++)
  {
    uint32_t Nbr;
    uint32_t MontDig;
    Pr = (Nbr = (pNbr1 + i)->x) * (uint64_t)Nbr2_0 + Prod0;
    MontDig = ((uint32_t)Pr * MontgomeryMultN[0].x) & MAX_INT_NBR;
    Prod0 = (Pr = (((uint64_t)MontDig * TestNbr0 + Pr) >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr1 + (uint64_t)Nbr * Nbr2_1 + Prod1) & MAX_INT_NBR;
    Prod1 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr2 + (uint64_t)Nbr * Nbr2_2 + Prod2) & MAX_INT_NBR;
    Prod2 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr3 + (uint64_t)Nbr * Nbr2_3 + Prod3) & MAX_INT_NBR;
    Prod3 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr4 + (uint64_t)Nbr * Nbr2_4 + Prod4) & MAX_INT_NBR;
    Prod4 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr5 + (uint64_t)Nbr * Nbr2_5 + Prod5) & MAX_INT_NBR;
    Prod5 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr6 + (uint64_t)Nbr * Nbr2_6 + Prod6) & MAX_INT_NBR;
    Prod6 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr7 + (uint64_t)Nbr * Nbr2_7 + Prod7) & MAX_INT_NBR;
    Prod7 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr8 + (uint64_t)Nbr * Nbr2_8 + Prod8) & MAX_INT_NBR;
    Prod8 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr9 + (uint64_t)Nbr * Nbr2_9 + Prod9) & MAX_INT_NBR;
    Prod9 = (Pr = (Pr >> BITS_PER_GROUP) +
      (uint64_t)MontDig * TestNbr10 + (uint64_t)Nbr * Nbr2_10 + Prod10) & MAX_INT_NBR;
    Prod10 = (uint32_t)(Pr >> BITS_PER_GROUP);
  }
  if (Pr >= ((uint64_t)(TestNbr10 + 1) << BITS_PER_GROUP)
    || (Prod10 == TestNbr10
      && (Prod9 > TestNbr9
        || (Prod9 == TestNbr9
          && (Prod8 > TestNbr8
            || (Prod8 == TestNbr8
              && (Prod7 > TestNbr7
                || (Prod7 == TestNbr7
                  && (Prod6 > TestNbr6
                    || (Prod6 == TestNbr6
                      && (Prod5 > TestNbr5
                        || (Prod5 == TestNbr5
                          && (Prod4 > TestNbr4
                            || (Prod4 == TestNbr4
                              && (Prod3 > TestNbr3
                                || (Prod3 == TestNbr3
                                  && (Prod2 > TestNbr2
                                    || (Prod2 == TestNbr2
                                      && (Prod1 > TestNbr1
                                        || (Prod1 == TestNbr1
                                          && (Prod0 >= TestNbr0)))))))))))))))))))))
  {
    int32_t borrow;
    Prod0 = (borrow = (int32_t)Prod0 - (int32_t)TestNbr0) & MAX_INT_NBR;
    Prod1 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod1 - (int32_t)TestNbr1) & MAX_INT_NBR;
    Prod2 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod2 - (int32_t)TestNbr2) & MAX_INT_NBR;
    Prod3 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod3 - (int32_t)TestNbr3) & MAX_INT_NBR;
    Prod4 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod4 - (int32_t)TestNbr4) & MAX_INT_NBR;
    Prod5 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod5 - (int32_t)TestNbr5) & MAX_INT_NBR;
    Prod6 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod6 - (int32_t)TestNbr6) & MAX_INT_NBR;
    Prod7 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod7 - (int32_t)TestNbr7) & MAX_INT_NBR;
    Prod8 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod8 - (int32_t)TestNbr8) & MAX_INT_NBR;
    Prod9 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod9 - (int32_t)TestNbr9) & MAX_INT_NBR;
    Prod10 = ((borrow >> BITS_PER_GROUP) + (int32_t)Prod10 - (int32_t)TestNbr10) & MAX_INT_NBR;
  }
  pProd->x = Prod0;
  (pProd + 1)->x = Prod1;
  (pProd + 2)->x = Prod2;
  (pProd + 3)->x = Prod3;
  (pProd + 4)->x = Prod4;
  (pProd + 5)->x = Prod5;
  (pProd + 6)->x = Prod6;
  (pProd + 7)->x = Prod7;
  (pProd + 8)->x = Prod8;
  (pProd + 9)->x = Prod9;
  (pProd + 10)->x = Prod10;
}
#endif

void modmult(const limb *factor1, const limb *factor2, limb *product)
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
  if (powerOf2Exponent != 0)
  {    // TestNbr is a power of 2.
    UncompressLimbsBigInteger(factor1, &tmpFact1);
    UncompressLimbsBigInteger(factor2, &tmpFact2);
    BigIntMultiply(&tmpFact1, &tmpFact2, &tmpFact1);
    CompressLimbsBigInteger(product, &tmpFact1);
    (product + powerOf2Exponent / BITS_PER_GROUP)->x &= (1 << (powerOf2Exponent % BITS_PER_GROUP)) - 1;
    return;
  }
  if (NumberLength == 1)
  {
    if (TestNbr[0].x <= 32768)
    {
      product->x = factor1->x * factor2->x % TestNbr[0].x;
      return;
    }
    smallmodmult(factor1->x, factor2->x, product, TestNbr[0].x);
    return;
  }
  if (NumberLength <= 12)
  {     // Small numbers.
    int j;
#ifdef _USING64BITS_
    switch (NumberLength)
    {
    case 2:
      MontgomeryMult2(factor1, factor2, product);
      return;
    case 3:
      MontgomeryMult3(factor1, factor2, product);
      return;
    case 4:
      MontgomeryMult4(factor1, factor2, product);
      return;
    case 5:
      MontgomeryMult5(factor1, factor2, product);
      return;
    case 6:
      MontgomeryMult6(factor1, factor2, product);
      return;
    case 7:
      MontgomeryMult7(factor1, factor2, product);
      return;
    case 8:
      MontgomeryMult8(factor1, factor2, product);
      return;
    case 9:
      MontgomeryMult9(factor1, factor2, product);
      return;
    case 10:
      MontgomeryMult10(factor1, factor2, product);
      return;
    case 11:
      MontgomeryMult11(factor1, factor2, product);
      return;
    default:
      break;
    }
    (void)memset(Prod, 0, NumberLength * sizeof(limb));
    for (int i = 0; i < NumberLength; i++)
    {
      int32_t MontDig;
      int32_t Nbr;
      int64_t Pr;
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
    double dInvLimbRange = 1.0 / dLimbRange;
    (void)memset(Prod, 0, NumberLength*sizeof(limb));
    for (int i = 0; i < NumberLength; i++)
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
          dAccum = ((dAccum + (double)FOURTH_INT_RANGE)*dInvLimbRange);
        }
        else
        {
          dAccum = ((dAccum - (double)FOURTH_INT_RANGE)*dInvLimbRange);
        }
        low = (int)(dAccum - floor(dAccum * dInvLimbRange) * dLimbRange);
        dAccum += dMontDig * TestNbr[j].x + dNbr * (factor2 + j)->x + (unsigned int)Prod[j].x;
        low = (low + MontDig * TestNbr[j].x +
               Nbr * (factor2 + j)->x + Prod[j].x) & MAX_VALUE_LIMB;
        Prod[j - 1].x = low;
      }
      if (low < HALF_INT_RANGE)
      {
        dAccum = ((dAccum + (double)FOURTH_INT_RANGE)*dInvLimbRange);
      }
      else
      {
        dAccum = ((dAccum - (double)FOURTH_INT_RANGE)*dInvLimbRange);
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
    if ((j<0) || ((unsigned int)Prod[j].x >= (unsigned int)TestNbr[j].x))
    {        // Prod >= TestNbr, so perform Prod <- Prod - TestNbr
      carry.x = 0;
      for (count = 0; count < NumberLength; count++)
      {
        carry.x += Prod[count].x - TestNbr[count].x;
        Prod[count].x = carry.x & MAX_VALUE_LIMB;
        carry.x >>= BITS_PER_GROUP;
      }
    }
    (void)memcpy(product, Prod, NumberLength*sizeof(limb));
    return;
  }
  // Compute T
  multiply(factor1, factor2, product, NumberLength, NULL);
  // Compute m
  multiply(product, MontgomeryMultN, aux, NumberLength, NULL);
  // Compute mN
  multiply(aux, TestNbr, aux2, NumberLength, NULL);
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
  if ((cy >= LIMB_RANGE) || ((product + count)->x >= TestNbr[count].x))
  {  // The number is greater or equal than TestNbr. Subtract it.
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
void modmultIntExtended(limb *factorBig, int factorInt, limb *result, const limb *pTestNbr, int nbrLen)
{
#ifdef _USING64BITS_
  int64_t carry;
#else
  double dTrialQuotient;
  double dAccumulator;
  double dFactorInt;
  double dInvLimbRange = 1.0 / (double)LIMB_RANGE;
  int low;
#endif
  int i;
  int TrialQuotient;
  limb *ptrFactorBig;
  const limb *ptrTestNbr;
  double dTestNbr;
  double dFactorBig;
  if (nbrLen == 1)
  {
    smallmodmult(factorBig->x, factorInt, result, pTestNbr->x);
    return;
  }
  (factorBig + nbrLen)->x = 0;
  dTestNbr = getMantissa(pTestNbr + nbrLen, nbrLen);
  dFactorBig = getMantissa(factorBig + nbrLen, nbrLen);
  TrialQuotient = (int)(unsigned int)floor(dFactorBig * (double)factorInt / dTestNbr + 0.5);
  if ((unsigned int)TrialQuotient >= LIMB_RANGE)
  {   // Maximum value for limb.
    TrialQuotient = MAX_VALUE_LIMB;
  }
  // Compute result <- factorBig * factorInt - TrialQuotient * TestNbr
  ptrFactorBig = factorBig;
  ptrTestNbr = pTestNbr;
#ifdef _USING64BITS_
  carry = 0;
  for (i = 0; i <= nbrLen; i++)
  {
    carry += (int64_t)ptrFactorBig->x * factorInt -
             (int64_t)TrialQuotient * ptrTestNbr->x;
    (result + i)->x = (int)carry & MAX_INT_NBR;
    carry >>= BITS_PER_GROUP;
    ptrFactorBig++;
    ptrTestNbr++;
  }
#else
  dFactorInt = (double)factorInt;
  dTrialQuotient = (double)TrialQuotient;
  low = 0;
  dAccumulator = 0;
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
#endif
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

// Compute power = base^exponent (mod modulus)
// Assumes GetMontgomeryParms routine for modulus already called.
// This works only for odd moduli.
void BigIntModularPower(const BigInteger *base, const BigInteger *exponent, BigInteger *power)
{
  CompressLimbsBigInteger(aux5, base);
  modmult(aux5, MontgomeryMultR2, aux6);   // Convert base to Montgomery notation.
  modPow(aux6, exponent->limbs, exponent->nbrLimbs, aux5);
  (void)memset(aux4, 0, NumberLength * sizeof(limb)); // Convert power to standard notation.
  aux4[0].x = 1;
  modmult(aux4, aux5, aux6);
  UncompressLimbsBigInteger(aux6, power);
}

// Input: base = base in Montgomery notation.
//        exp  = exponent.
//        nbrGroupsExp = number of limbs of exponent.
// Output: power = power in Montgomery notation.
void modPow(const limb *base, const limb *exp, int nbrGroupsExp, limb *power)
{
  (void)memcpy(power, MontgomeryMultR1, (NumberLength + 1)*sizeof(*power));  // power <- 1
  for (int index = nbrGroupsExp - 1; index >= 0; index--)
  {
    int groupExp = (exp + index)->x;
    for (int mask = 1 << (BITS_PER_GROUP - 1); mask > 0; mask >>= 1)
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
void modPowLimb(const limb *base, const limb *exp, limb *power)
{
  int groupExp;
  (void)memcpy(power, MontgomeryMultR1, (NumberLength + 1)*sizeof(*power));  // power <- 1
  groupExp = exp->x;
  for (int mask = 1 << (BITS_PER_GROUP - 1); mask > 0; mask >>= 1)
  {
    modmult(power, power, power);
    if ((groupExp & mask) != 0)
    {
      modmult(power, base, power);
    }
  }
}

void modPowBaseInt(int base, const limb *exp, int nbrGroupsExp, limb *power)
{
  (void)memcpy(power, MontgomeryMultR1, (NumberLength+1)*sizeof(limb));  // power <- 1
  for (int index = nbrGroupsExp-1; index>=0; index--)
  {
    int groupExp = (exp+index)->x;
    for (int mask = 1<<(BITS_PER_GROUP-1); mask > 0; mask >>= 1)
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
#ifdef _USING64BITS_
  int64_t carryU = 0;
  int64_t carryV = 0;
  for (int ctr = 0; ctr <= nbrLen; ctr++)
  {
    int u = firstBig->x;
    int v = secondBig->x;
    carryU += u*(int64_t)e + v*(int64_t)f;
    carryV += u*(int64_t)g + v*(int64_t)h;
    firstBig->x = (int)(carryU & MAX_INT_NBR);
    secondBig->x = (int)(carryV & MAX_INT_NBR);
    firstBig++;
    secondBig++;
    carryU >>= BITS_PER_GROUP;
    carryV >>= BITS_PER_GROUP;
  }
#else
  double dVal = 1.0 / (double)LIMB_RANGE;
  int carryU;
  int carryV;
  double dFactorE = (double)e;
  double dFactorF = (double)f;
  double dFactorG = (double)g;
  double dFactorH = (double)h;
  carryU = carryV = 0;
  for (int ctr = 0; ctr <= nbrLen; ctr++)
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
    firstBig->x = lowU;
    secondBig->x = lowV;
    firstBig++;
    secondBig++;
  }
#endif
}

// Perform first <- (first - second) / 2
// first must be greater than second.
static int HalveDifference(limb *first, const limb *second, int len)
{
  int i;
  int borrow;
  int prevLimb;
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
  return len + 1;
}

int modInv(int NbrMod, int currentPrime)
{
  int QQ;
  int T1;
  int T3;
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

static void InitHighUandV(int lenU, int lenV, double *pHighU, double *pHighV)
{
  double highU;
  double highV;
  double dLimbRange = (double)LIMB_RANGE;
  if (lenV >= lenU)
  {
    highV = ((double)V[lenV - 1].x * dLimbRange) + (double)V[lenV - 2].x;
    if (lenV >= 3)
    {
      highV += (double)V[lenV - 3].x / dLimbRange;
    }
    if (lenV == lenU)
    {
      highU = ((double)U[lenV - 1].x * dLimbRange) + (double)U[lenV - 2].x;
    }
    else if (lenV == (lenU + 1))
    {
      highU = (double)U[lenV - 2].x;
    }
    else
    {
      highU = 0;
    }
    if ((lenV <= (lenU + 2)) && (lenV >= 3))
    {
      highU += (double)U[lenV - 3].x / dLimbRange;
    }
  }
  else
  {
    highU = ((double)U[lenU - 1].x * (double)LIMB_RANGE) + (double)U[lenU - 2].x;
    if (lenU >= 3)
    {
      highU += (double)U[lenU - 3].x / dLimbRange;
    }
    if (lenU == (lenV + 1))
    {
      highV = (double)V[lenU - 2].x;
    }
    else
    {
      highV = 0;
    }
    if ((lenU <= (lenV + 2)) && (lenU >= 3))
    {
      highV += (double)V[lenU - 3].x / dLimbRange;
    }
  }
  *pHighU = highU;
  *pHighV = highV;
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
/* R' <- aR + bS, S' <-  cR + dS                                       */
/* U' <- aU - bV, V' <- -cU + dV                                       */
/***********************************************************************/
void ModInvBigNbr(limb *num, limb *inv, limb *mod, int nbrLen)
{
  int k;
  int steps;
  int a;
  int b;
  int c;
  int d;  // Coefficients used to update variables R, S, U, V.
  int size;
  int i;
  int bitCount;
  int lenRS;
  int lenU;
  int lenV;
  int lowU;
  int lowV;
  int borrow;
  if (nbrLen == 1)
  {
    inv->x = modInv(num->x, mod->x);
    return;
  }
  if (powerOf2Exponent != 0)
  {    // TestNbr is a power of 2.
    ComputeInversePower2(num, inv, aux);
    (inv + powerOf2Exponent / BITS_PER_GROUP)->x &= (1 << (powerOf2Exponent % BITS_PER_GROUP)) - 1;
    return;
  }
  //  1. U <- M, V <- X, R <- 0, S <- 1, k <- 0
  size = (nbrLen+1)*sizeof(limb);
  (mod + nbrLen)->x = 0;
  (num + nbrLen)->x = 0;
  (void)memcpy(U, mod, size);
  (void)memcpy(V, num, size);
    // Maximum value of R and S can be up to 2*M, so one more limb is needed.
  (void)memset(R, 0, size);   // R <- 0
  (void)memset(S, 0, size);   // S <- 1
  S[0].x = 1;
  lenRS = 1;
  k = 0;
  steps = 0;
  // R' <- aR + bS, S' <- cR + dS
  a = 1;  // R' = R, S' = S.
  d = 1;
  b = 0;
  c = 0;
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
  lowU = U[0].x;
  lowV = V[0].x;
  // Initialize highU and highV.
  if ((lenU > 1) || (lenV > 1))
  {
    double highU;
    double highV;
    InitHighUandV(lenU, lenV, &highU, &highV);
    //  2. while V > 0 do
    for (;;)
    {
      //  3.   if U even then U <- U / 2, S <- 2S
      if ((lowU & 1) == 0)
      {     // U is even.
        lowU >>= 1;
        highV += highV;
        // R' <- aR + bS, S' <- cR + dS
        c *= 2;
        d *= 2;  // Multiply S by 2.
      }
      //  4.   elsif V even then V <- V / 2, R <- 2R
      else if ((lowV & 1) == 0)
      {    // V is even.
        lowV >>= 1;
        highU += highU;
        // R' <- aR + bS, S' <- cR + dS
        a *= 2;
        b *= 2;  // Multiply R by 2.
      }
      else
      {
        //  5.   elsif U >= V  then U <- (U - V) / 2, R <- R + S, S <- 2S
        if (highU > highV)
        {     // U > V. Perform U <- (U - V) / 2
          lowU = (lowU - lowV) >> 1;
          highU -= highV;
          highV += highV;
          // R' <- aR + bS, S' <- cR + dS
          a += c;
          b += d;  // R <- R + S
          c *= 2;
          d *= 2;  // S <- 2S
        }
        //  6.   elsif V >= U then V <- (V - U) / 2, S <- S + R, R <- 2R
        else
        {    // V >= U. Perform V <- (V - U) / 2
          lowV = (lowV - lowU) >> 1;
          highV -= highU;
          highU += highU;
          // R' <- aR + bS, S' <- cR + dS
          c += a;
          d += b;  // S <- S + R
          a *= 2;
          b *= 2;  // R <- 2R
        }
      }
      //  7.   k <- k + 1
      // Adjust variables.
      steps++;
      if (steps == (BITS_PER_GROUP - 1))
      {  // compute now U and V and reset e, f, g and h.
         // U' <- eU + fV, V' <- gU + hV
        int len = (lenU > lenV ? lenU : lenV);
        (void)memset(&U[lenU].x, 0, (len - lenU + 1) * sizeof(limb));
        (void)memset(&V[lenV].x, 0, (len - lenV + 1) * sizeof(limb));
        (void)memcpy(Ubak, U, (len + 1) * sizeof(limb));
        (void)memcpy(Vbak, V, (len + 1) * sizeof(limb));
        AddMult(U, a, -b, V, -c, d, len);
        if ((U[lenU].x | V[lenV].x) & (1<<(BITS_PER_GROUP - 2)))
        {    // Complete expansion of U and V required for all steps.
            //  2. while V > 0 do
          (void)memcpy(U, Ubak, (len+1) * sizeof(limb));
          (void)memcpy(V, Vbak, (len+1) * sizeof(limb));
          b = 0;
          c = 0;  // U' = U, V' = V.
          a = 1;
          d = 1;
          while ((lenV > 1) || (V[0].x > 0))
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
              c *= 2;
              d *= 2;  // Multiply S by 2.
            }
            //  4.   elsif V even then V <- V / 2, R <- 2R
            else if ((V[0].x & 1) == 0)
            {    // V is even.
              for (i = 0; i < lenV; i++)
              {  // Loop that divides V by 2.
                V[i].x = ((V[i].x >> 1) | (V[i + 1].x << (BITS_PER_GROUP - 1))) & MAX_VALUE_LIMB;
              }
              if (V[lenV - 1].x == 0)
              {
                lenV--;
              }
              // R' <- aR + bS, S' <- cR + dS
              a *= 2;
              b *= 2;  // Multiply R by 2.
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
                a += c;
                b += d;  // R <- R + S
                c *= 2;
                d *= 2;  // S <- 2S
              }
              //  6.   elsif V >= U then V <- (V - U) / 2, S <- S + R, R <- 2R
              else
              {    // V >= U
                lenV = HalveDifference(V, U, len); // V <- (V - U) / 2
                                                   // R' <- aR + bS, S' <- cR + dS
                c += a;
                d += b;  // S <- S + R
                a *= 2;
                b *= 2;  // R <- 2R
              }
            }
            //  7.   k <- k + 1
            k++;
            if (k % (BITS_PER_GROUP - 1) == 0)
            {
              break;
            }
          }
          if ((lenV == 1) && (V[0].x == 0))
          {
            break;
          }
        }
        else
        {
          k += steps;
          for (i = 0; i < lenU; i++)
          {  // Loop that divides U by 2^(BITS_PER_GROUP - 1).
            U[i].x = ((U[i].x >> (BITS_PER_GROUP - 1)) | (U[i + 1].x << 1)) & MAX_VALUE_LIMB;
          }
          U[lenU].x = 0;
          while ((lenU > 0) && (U[lenU - 1].x == 0))
          {
            lenU--;
          }
          for (i = 0; i < lenV; i++)
          {  // Loop that divides V by 2^(BITS_PER_GROUP - 1).
            V[i].x = ((V[i].x >> (BITS_PER_GROUP - 1)) | (V[i + 1].x << 1)) & MAX_VALUE_LIMB;
          }
          V[lenV].x = 0;
          while ((lenV > 0) && (V[lenV - 1].x == 0))
          {
            lenV--;
          }
        }
        steps = 0;
        AddMult(R, a, b, S, c, d, lenRS);
        if ((R[lenRS].x != 0) || (S[lenRS].x != 0))
        {
          lenRS++;
        }
        lowU = U[0].x;
        lowV = V[0].x;
        b = 0;
        c = 0;  // U' = U, V' = V.
        a = 1;
        d = 1;
        if ((lenU == 0) || (lenV == 0) || ((lenV == 1) && (lenU == 1)))
        {
          break;
        }
        InitHighUandV(lenU, lenV, &highU, &highV);
      }
    }
  }
  if (lenU > 0)
  {
  //  2. while V > 0 do
    while (lowV > 0)
    {
      //  3.   if U even then U <- U / 2, S <- 2S
      if ((lowU & 1) == 0)
      {     // U is even.
        lowU >>= 1;
        // R' <- aR + bS, S' <- cR + dS
        c *= 2;
        d *= 2;  // Multiply S by 2.
      }
      //  4.   elsif V even then V <- V / 2, R <- 2R
      else if ((lowV & 1) == 0)
      {    // V is even.
        lowV >>= 1;
        // R' <- aR + bS, S' <- cR + dS
        a *= 2;
        b *= 2;  // Multiply R by 2.
      }
      //  5.   elsif U >= V  then U <- (U - V) / 2, R <- R + S, S <- 2S
      else if (lowU > lowV)
      {     // U > V. Perform U <- (U - V) / 2
        lowU = (lowU - lowV) >> 1;
        // R' <- aR + bS, S' <- cR + dS
        a += c;
        b += d;  // R <- R + S
        c *= 2;
        d *= 2;  // S <- 2S
      }
      //  6.   elsif V >= U then V <- (V - U) / 2, S <- S + R, R <- 2R
      else
      {    // V >= U. Perform V <- (V - U) / 2
        lowV = (lowV - lowU) >> 1;
        // R' <- aR + bS, S' <- cR + dS
        c += a;
        d += b;  // S <- S + R
        a *= 2;
        b *= 2;  // R <- 2R
      }
      //  7.   k <- k + 1
      steps++;
      if (steps == (BITS_PER_GROUP - 1))
      {  // compute now R and S and reset a, b, c and d.
         // R' <- aR + bS, S' <- cR + dS
        AddMult(R, a, b, S, c, d, nbrLen + 1);
        b = 0;     // R' = R, S' = S.
        c = 0;
        a = 1;
        d = 1;
        k += steps;
        steps = 0;
      }
    }
  }
  AddMult(R, a, b, S, c, d, nbrLen+1);
  k += steps;
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
  (void)memset(S, 0, size);
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

// Compute modular division for odd moduli.
void BigIntModularDivision(const BigInteger *Num, const BigInteger *Den,
  const BigInteger *mod, BigInteger *quotient)
{
  NumberLength = mod->nbrLimbs;
  // Reduce Num modulo mod.
  BigIntRemainder(Num, mod, &tmpNum);
  if (tmpNum.sign == SIGN_NEGATIVE)
  {
    BigIntAdd(&tmpNum, mod, &tmpNum);
  }
  // Reduce Den modulo mod.
  BigIntRemainder(Den, mod, &tmpDen);
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

// Modular division when modulus is a power of 2.
void BigIntModularDivisionPower2(const BigInteger *Num, const BigInteger *Den,
   const BigInteger *mod, BigInteger *quotient)
{
  int NumberLengthBak = NumberLength;
  NumberLength = mod->nbrLimbs;
  // Compute aux3 as inverse of Den mod mod.
  ComputeInversePower2(Den->limbs, aux3, aux4);
  if (Num->sign != Den->sign)
  {   // Sign of numerator is different from divisor, so negate aux3.
    int Cy = 0;
    for (int idx = 0; idx < NumberLength; idx++)
    {
      Cy -= aux3[idx].x;
      aux3[idx].x = Cy & MAX_VALUE_LIMB;
      Cy >>= BITS_PER_GROUP;
    }
  }
  multiply(aux3, Num->limbs, quotient->limbs, NumberLength, NULL);    // quotient <- Den * aux3
  quotient->limbs[NumberLength - 1].x &= mod->limbs[NumberLength - 1].x - 1;
  // Adjust number of length of quotient so the most significant limb is not zero.
  while (NumberLength > 1)
  {
    if (quotient->limbs[NumberLength - 1].x)
    {
      break;
    }
    NumberLength--;
  }
  quotient->nbrLimbs = NumberLength;
  NumberLength = NumberLengthBak;
}

void BigIntModularDivisionSaveTestNbr(const BigInteger *Num, const BigInteger *Den,
  const BigInteger *mod, BigInteger *quotient)
{
  int NumberLengthBak = NumberLength;
  (void)memcpy(U, TestNbr, (NumberLength + 1) * sizeof(limb));
  NumberLength = mod->nbrLimbs;
  (void)memcpy(TestNbr, mod->limbs, NumberLength * sizeof(limb));
  TestNbr[NumberLength].x = 0;
  GetMontgomeryParms(NumberLength);
  BigIntModularDivision(Num, Den, mod, quotient);
  NumberLength = NumberLengthBak;
  (void)memcpy(TestNbr, U, (NumberLength+1) * sizeof(limb));
}

// On input: 
// oddValue = odd modulus.
// resultModOdd = result mod odd value
// resultModPower2 = result mod 2^shRight
// result = pointer to result.
// From Knuth's TAOCP Vol 2, section 4.3.2:
// If c = result mod odd, d = result mod 2^k:
// compute result = c + (d-c)*modinv(odd,2^k)*odd
static void ChineseRemainderTheorem(int shRight, BigInteger *result)
{
  if (shRight == 0)
  {
    NumberLength = oddValue.nbrLimbs;
    UncompressLimbsBigInteger(resultModOdd, result);
    return;
  }
  if (NumberLength > oddValue.nbrLimbs)
  {
    (void)memset(&oddValue.limbs[oddValue.nbrLimbs], 0, (NumberLength - oddValue.nbrLimbs) * sizeof(limb));
  }
  SubtractBigNbr((int *)resultModPower2, (int *)resultModOdd, (int *)aux3, NumberLength);
  ComputeInversePower2(oddValue.limbs, aux4, aux);
  modmult(aux4, aux3, aux5);
  (aux5 + shRight / BITS_PER_GROUP)->x &= (1 << (shRight % BITS_PER_GROUP)) - 1;
  UncompressLimbsBigInteger(aux5, result);
  BigIntMultiply(result, &oddValue, result);
  NumberLength = oddValue.nbrLimbs;
  UncompressLimbsBigInteger(resultModOdd, &tmpDen);
  BigIntAdd(result, &tmpDen, result);
}

// Compute modular division. ModInvBigNbr does not support even moduli,
// so the division is done separately by calculating the division modulo
// n/2^k (n odd) and 2^k and then merge the results using Chinese Remainder
// Theorem.
void BigIntGeneralModularDivision(const BigInteger *Num, const BigInteger *Den, 
   const BigInteger *mod, BigInteger *quotient)
{
  int shRight;
  CopyBigInt(&oddValue, mod);
  DivideBigNbrByMaxPowerOf2(&shRight, oddValue.limbs, &oddValue.nbrLimbs);
  // Reduce Num modulo oddValue.
  BigIntRemainder(Num, &oddValue, &tmpNum);
  if (tmpNum.sign == SIGN_NEGATIVE)
  {
    BigIntAdd(&tmpNum, &oddValue, &tmpNum);
  }
  // Reduce Den modulo oddValue.
  BigIntRemainder(Den, &oddValue, &tmpDen);
  if (tmpDen.sign == SIGN_NEGATIVE)
  {
    BigIntAdd(&tmpDen, &oddValue, &tmpDen);
  }
  NumberLength = oddValue.nbrLimbs;
  (void)memcpy(TestNbr, oddValue.limbs, NumberLength * sizeof(limb));
  TestNbr[NumberLength].x = 0;
  GetMontgomeryParms(NumberLength);
  CompressLimbsBigInteger(aux3, &tmpDen);
  modmult(aux3, MontgomeryMultR2, aux3);      // aux3 <- Den in Montgomery notation
  ModInvBigNbr(aux3, aux3, TestNbr, NumberLength); // aux3 <- 1 / Den in Montg notation.
  CompressLimbsBigInteger(aux4, &tmpNum);
  modmult(aux3, aux4, resultModOdd);          // resultModOdd <- Num / Dev in standard notation.

  // Compute inverse mod power of 2.
  NumberLength = (shRight + BITS_PER_GROUP-1) / BITS_PER_GROUP;
  CompressLimbsBigInteger(aux3, Den);
  ComputeInversePower2(aux3, aux4, aux);
  powerOf2Exponent = shRight;
  modmult(Num->limbs, aux4, resultModPower2); // resultModPower2 <- Num / Dev modulus 2^k.
  ChineseRemainderTheorem(shRight, quotient);
  powerOf2Exponent = 0;
}

// Compute modular division. ModInvBigNbr does not support even moduli,
// so the division is done separately by calculating the division modulo
// n/2^k (n odd) and 2^k and then merge the results using Chinese Remainder
// Theorem.
enum eExprErr BigIntGeneralModularPower(const BigInteger *base, const BigInteger *exponent, 
  const BigInteger *mod, BigInteger *power)
{
  int shRight;
  if ((mod->nbrLimbs == 1) && (mod->limbs[0].x == 0))
  {            // Modulus is zero.
    return BigIntPower(base, exponent, power);
  }
  CopyBigInt(&oddValue, mod);
  oddValue.sign = SIGN_POSITIVE;
  DivideBigNbrByMaxPowerOf2(&shRight, oddValue.limbs, &oddValue.nbrLimbs);
  // Reduce base modulo oddValue.
  BigIntRemainder(base, &oddValue, &tmpNum);
  if (tmpNum.sign == SIGN_NEGATIVE)
  {
    BigIntAdd(&tmpNum, &oddValue, &tmpNum);
  }
  NumberLength = oddValue.nbrLimbs;
  (void)memcpy(TestNbr, oddValue.limbs, NumberLength * sizeof(limb));
  TestNbr[NumberLength].x = 0;
  GetMontgomeryParms(NumberLength);
  BigIntModularPower(&tmpNum, exponent, &tmpDen);
  (void)memcpy(resultModOdd, tmpDen.limbs, tmpDen.nbrLimbs * sizeof(limb));
  if (shRight > 0)
  {
    // Compute power mod power of 2.
    NumberLength = (shRight + BITS_PER_GROUP - 1) / BITS_PER_GROUP;
    CompressLimbsBigInteger(aux3, base);
    powerOf2Exponent = shRight;
    modPowLimb(aux3, exponent->limbs, resultModPower2);
    ChineseRemainderTheorem(shRight, power);
    powerOf2Exponent = 0;
  }
  else
  {
    CopyBigInt(power, &tmpDen);
  }
  return EXPR_OK;
}

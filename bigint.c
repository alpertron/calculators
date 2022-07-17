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
#include <string.h>
#include <math.h>
#include "bignbr.h"

#define MAX_LIMBS_SIQS 15
static BigInteger Numerator;
static BigInteger Denominator;
static BigInteger Modulus;
static BigInteger Quotient;
static BigInteger BigInt1;
static BigInteger BigInt2;
static BigInteger BigGcd;

limb NbrBak[MAX_LIMBS_SIQS];
void ChSignBigNbr(limb *nbr, int length)
{
  int carry = 0;
  const limb *ptrEndNbr = nbr + length;
  for (limb *ptrNbr = nbr; ptrNbr < ptrEndNbr; ptrNbr++)
  {
    carry -= ptrNbr->x;
    ptrNbr->x = carry & MAX_INT_NBR;
    carry >>= BITS_PER_GROUP;
  }
}

void ChSignBigNbrB(limb *nbr, int length)
{
  int carry = 0;
  limb *ptrNbr;
  const limb *ptrEndNbr = nbr + length - 1;
  for (ptrNbr = nbr; ptrNbr < ptrEndNbr; ptrNbr++)
  {
    carry -= ptrNbr->x;
    ptrNbr->x = carry & MAX_INT_NBR;
    carry >>= BITS_PER_GROUP;
  }
  ptrNbr->x = carry - ptrNbr->x;
}

void AddBigNbr(const limb *pNbr1, const limb*pNbr2, limb*pSum, int nbrLen)
{
  unsigned int carry = 0U;
  const limb*ptrNbr1 = pNbr1;
  const limb*ptrNbr2 = pNbr2;
  const limb*ptrEndSum = pSum + nbrLen;
  for (limb*ptrSum = pSum; ptrSum < ptrEndSum; ptrSum++)
  {
    unsigned int tmp;
    carry = (carry >> BITS_PER_GROUP) + (unsigned int)ptrNbr1->x +
      (unsigned int)ptrNbr2->x;
    tmp = carry & MAX_INT_NBR_U;
    ptrSum->x = (int)tmp;
    ptrNbr1++;
    ptrNbr2++;
  }
}

void SubtractBigNbr(const limb *pNbr1, const limb*pNbr2, limb*pDiff, int nbrLen)
{
  unsigned int borrow = 0U;
  const limb*ptrNbr1 = pNbr1;
  const limb*ptrNbr2 = pNbr2;
  const limb*ptrEndDiff = pDiff + nbrLen;
  for (limb*ptrDiff = pDiff; ptrDiff < ptrEndDiff; ptrDiff++)
  {
    unsigned int tmp;
    borrow = (unsigned int)ptrNbr1->x - (unsigned int)ptrNbr2->x -
      (borrow >> BITS_PER_GROUP);
    tmp = borrow & MAX_INT_NBR_U;
    ptrDiff->x = (int)tmp;
    ptrNbr1++;
    ptrNbr2++;
  }
}

void AddBigNbrB(const limb*pNbr1, const limb*pNbr2, limb*pSum, int nbrLen)
{
  unsigned int carry = 0U;
  const limb*ptrNbr1 = pNbr1;
  const limb*ptrNbr2 = pNbr2;
  limb *ptrSum;
  const limb* ptrEndSum = pSum + nbrLen - 1;
  for (ptrSum = pSum; ptrSum < ptrEndSum; ptrSum++)
  {
    unsigned int tmp;
    carry = (carry >> BITS_PER_GROUP) + (unsigned int)ptrNbr1->x +
      (unsigned int)ptrNbr2->x;
    tmp = carry & MAX_INT_NBR_U;
    ptrSum->x = (int)tmp;
    ptrNbr1++;
    ptrNbr2++;
  }
  carry = (carry >> BITS_PER_GROUP) + (unsigned int)ptrNbr1->x +
    (unsigned int)ptrNbr2->x;
  ptrSum->x = (int)carry;
}

void SubtractBigNbrB(const limb *pNbr1, const limb *pNbr2, limb *pDiff, int nbrLen)
{
  unsigned int borrow = 0U;
  const limb *ptrNbr1 = pNbr1;
  const limb *ptrNbr2 = pNbr2;
  limb *ptrDiff;
  const limb* ptrEndDiff = pDiff + nbrLen - 1;
  for (ptrDiff = pDiff; ptrDiff < ptrEndDiff; ptrDiff++)
  {
    unsigned int tmp;
    borrow = (unsigned int)ptrNbr1->x - (unsigned int)ptrNbr2->x -
      (borrow >> BITS_PER_GROUP);
    tmp = borrow & MAX_INT_NBR_U;
    ptrDiff->x = (int)tmp;
    ptrNbr1++;
    ptrNbr2++;
  }
  borrow = (unsigned int)ptrNbr1->x - (unsigned int)ptrNbr2->x -
    (borrow >> BITS_PER_GROUP);
  ptrDiff->x = borrow;
}

void AddBigIntModN(const limb *pNbr1, const limb *pNbr2, limb *pSum, const limb *pMod,
  int nbrLen)
{
  const limb* ptrNbr1 = pNbr1;
  const limb* ptrNbr2 = pNbr2;
  limb* ptrSum = pSum;
  const limb* ptrMod = pMod;
  unsigned int borrow = 0U;
  unsigned int carry = 0U;
  unsigned int tmp;
  int i;
  for (i = 0; i < nbrLen; i++)
  {
    carry = (carry >> BITS_PER_GROUP) + (unsigned int)ptrNbr1->x +
      (unsigned int)ptrNbr2->x;
    tmp = carry & MAX_INT_NBR_U;
    ptrSum->x = (int)tmp;
    ptrNbr1++;
    ptrNbr2++;
    ptrSum++;
  }
  carry >>= BITS_PER_GROUP;
  ptrSum -= nbrLen;
  for (i = 0; i < nbrLen; i++)
  {
    borrow = (unsigned int)ptrSum->x - (unsigned int)ptrMod->x -
      (borrow >> BITS_PER_GROUP);
    tmp = borrow & MAX_INT_NBR_U;
    ptrSum->x = (int)tmp;
    ptrMod++;
    ptrSum++;
  }
  borrow >>= BITS_PER_GROUP;
  if (borrow != carry)
  {    // Sum is less than zero. Add Mod again.
    ptrSum -= nbrLen;
    ptrMod -= nbrLen;
    carry = 0;
    for (i = 0; i < nbrLen; i++)
    {
      carry = (carry >> BITS_PER_GROUP) + (unsigned int)ptrSum->x +
        (unsigned int)ptrMod->x;
      tmp = carry & MAX_INT_NBR_U;
      ptrSum->x = (int)tmp;
      ptrMod++;
      ptrSum++;
    }
  }
}

void SubtractBigNbrModN(const limb *pNbr1, const limb *pNbr2, limb *pDiff, const limb *pMod,
  int nbrLen)
{
  const limb* ptrNbr1 = pNbr1;
  const limb* ptrNbr2 = pNbr2;
  limb* ptrDiff = pDiff;
  const limb* ptrMod = pMod;
  unsigned int borrow = 0U;
  unsigned int tmp;
  int i;
  for (i = 0; i < nbrLen; i++)
  {
    borrow = (unsigned int)ptrNbr1->x - (unsigned int)ptrNbr2->x -
      (borrow >> BITS_PER_GROUP);
    tmp = borrow & MAX_INT_NBR_U;
    ptrDiff->x = (int)tmp;
    ptrNbr1++;
    ptrNbr2++;
    ptrDiff++;
  }
  if (borrow != 0U)
  {
    unsigned int carry = 0U;
    ptrDiff -= nbrLen;
    for (i = 0; i < nbrLen; i++)
    {
      carry = (carry >> BITS_PER_GROUP) + (unsigned int)ptrDiff->x +
        (unsigned int)ptrMod->x;
      tmp = carry & MAX_INT_NBR_U;
      ptrDiff->x = (int)tmp;
      ptrMod++;
      ptrDiff++;
    }
  }
}

void MultBigNbrByInt(const limb *pBigFactor, int factor, limb *bigProd, int nbrLen)
{
  int secondFactor = factor;
  const limb* ptrBigFactor = pBigFactor;
  limb *bigProduct = bigProd;
  double dFactor;
  double dVal = 1.0 / (double)(1U << BITS_PER_GROUP);
  int factorPositive = 1;
  int carry;
  if (secondFactor < 0)
  {     // If factor is negative, indicate it and compute its absolute value.
    factorPositive = 0;
    secondFactor = -secondFactor;
  }
  dFactor = (double)secondFactor;
  carry = 0;
  for (int ctr = 0; ctr < nbrLen; ctr++)
  {
    double dCarry;
    int low = ((ptrBigFactor->x * secondFactor) + carry) & MAX_INT_NBR;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    dCarry = ((double)ptrBigFactor->x * dFactor) + (double)carry;
    if (low < HALF_INT_RANGE)
    {
      dCarry = floor((dCarry + (double)FOURTH_INT_RANGE) * dVal);
    }
    else
    {
      dCarry = floor((dCarry - (double)FOURTH_INT_RANGE) * dVal);
    }
    carry = (int)dCarry;
    bigProduct->x = low;
    bigProduct++;
    ptrBigFactor++;
  }
  if (factorPositive == 0)
  {         // If factor is negative, change sign of product.
    ChSignBigNbr(bigProd, nbrLen);
  }
}

void MultBigNbrByIntB(const limb *bigFactor, int factor, limb *bigProd, int nbrLen)
{
  const limb* bigFact = bigFactor;
  int fact = factor;
  limb *bigProduct = bigProd;
  double dFactor;
  double dVal = 1.0 / (double)(1U << BITS_PER_GROUP);
  int factorPositive = 1;
  int carry;
  if (fact < 0)
  {     // If factor is negative, indicate it and compute its absolute value.
    factorPositive = 0;
    fact = -fact;
  }
  dFactor = (double)fact;
  carry = 0;
  for (int ctr = 0; ctr < (nbrLen-1); ctr++)
  {
    double dCarry;
    int low = ((bigFact->x * fact) + carry) & MAX_INT_NBR;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    dCarry = ((double)bigFact->x * dFactor) + (double)carry;
    if (low < HALF_INT_RANGE)
    {
      dCarry = floor((dCarry + (double)FOURTH_INT_RANGE) * dVal);
    }
    else
    {
      dCarry = floor((dCarry - (double)FOURTH_INT_RANGE) * dVal);
    }
    carry = (int)dCarry;
    bigProduct->x = low;
    bigProduct++;
    bigFact++;
  }
  bigProduct->x = (bigFact->x * fact) + carry;
  if (factorPositive == 0)
  {         // If factor is negative, change sign of product.
    ChSignBigNbrB(bigProd, nbrLen);
  }
}

void DivBigNbrByInt(const limb *pDividend, int divisor, limb *pQuotient, int nbrLen)
{
  const limb* ptrDividend = pDividend;
  limb* ptrQuotient = pQuotient;
  unsigned int remainder = 0U;
  double dDivisor = (double)divisor;
  double dLimb = 2147483648.0;
  int nbrLenMinus1 = nbrLen - 1;
  ptrDividend += nbrLenMinus1;
  ptrQuotient += nbrLenMinus1;
  for (int ctr = nbrLenMinus1; ctr >= 0; ctr--)
  {
    unsigned int dividend = (remainder << BITS_PER_GROUP) +
      (unsigned int)ptrDividend->x;
    double dDividend = ((double)remainder * dLimb) + (double)ptrDividend->x;
    // quotient has correct value or 1 more.
    unsigned int quotient = (unsigned int)((dDividend / dDivisor) + 0.5);
    remainder = dividend - (quotient * (unsigned int)divisor);
    if (remainder >= (unsigned int)divisor)
    {     // remainder not in range 0 <= remainder < divisor. Adjust.
      quotient--;
      remainder += (unsigned int)divisor;
    }
    ptrQuotient->x = (int)quotient;
    ptrQuotient--;
    ptrDividend--;
  }
}

int RemDivBigNbrByInt(const limb *pDividend, int divisor, int nbrLen)
{
  const limb* ptrDividend = pDividend;
  unsigned int remainder = 0U;
  double dDivisor = (double)divisor;
  double dLimb = 0x80000000;
  int nbrLenMinus1 = nbrLen - 1;
  ptrDividend += nbrLenMinus1;
  for (int ctr = nbrLenMinus1; ctr >= 0; ctr--)
  {
    unsigned int dividend = (remainder << BITS_PER_GROUP) +
      (unsigned int)ptrDividend->x;
    double dDividend = ((double)remainder * dLimb) + (double)ptrDividend->x;
         // quotient has correct value or 1 more.
    unsigned int quotient = (unsigned int)((dDividend / dDivisor) + 0.5);
    remainder = dividend - (quotient * (unsigned int)divisor);
    if (remainder >= (unsigned int)divisor)
    {     // remainder not in range 0 <= remainder < divisor. Adjust.
      quotient--;
      remainder += (unsigned int)divisor;
    }
    ptrDividend--;
  }
  return remainder;
}

void MultBigNbr(const limb *pFactor1, const limb *pFactor2, limb *pProd, int nbrLen)
{
  limb* ptrProd = pProd;
  double dRangeLimb = (double)(1U << BITS_PER_GROUP);
  double dInvRangeLimb = 1.0 / dRangeLimb;
  int low = 0;
  int factor1;
  int factor2;
  double dAccumulator = 0.0;
  for (int i = 0; i < nbrLen; i++)
  {
    for (int j = 0; j <= i; j++)
    {
      factor1 = (pFactor1 + j)->x;
      factor2 = (pFactor2 + i - j)->x;
      low += factor1*factor2;
      dAccumulator += (double)factor1 * (double)factor2;
    }
    low &= MAX_INT_NBR;    // Trim extra bits.
    ptrProd->x = low;
    ptrProd++;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    if (low < HALF_INT_RANGE)
    {
      dAccumulator = floor((dAccumulator + (double)FOURTH_INT_RANGE)*dInvRangeLimb);
    }
    else
    {
      dAccumulator = floor((dAccumulator - (double)FOURTH_INT_RANGE)*dInvRangeLimb);
    }
    low = (int)(dAccumulator - floor(dAccumulator * dInvRangeLimb) * dRangeLimb);
  }
  ptrProd->x = low;
  (ptrProd+1)->x = (int)floor(dAccumulator/dRangeLimb);
}

// On input: pFactor1 and pFactor2: pointers to factors.
//           pProd: pointer to product (length = 2*nbrLen)
//           nbrLen: number of limbs of factors.
void MultBigNbrComplete(const limb *pFactor1, const limb *pFactor2, limb *pProd, int nbrLen)
{
  limb* ptrProd = pProd;
  double dRangeLimb = (double)(1U << BITS_PER_GROUP);
  double dInvRangeLimb = 1.0 / dRangeLimb;
  int low = 0;
  int i;
  int j;
  int factor1;
  int factor2;
  double dAccumulator = 0.0;
  for (i = 0; i < nbrLen; i++)
  {
    for (j = 0; j <= i; j++)
    {
      factor1 = (pFactor1 + j)->x;
      factor2 = (pFactor2 + i - j)->x;
      low += factor1 * factor2;
      dAccumulator += (double)factor1 * (double)factor2;
    }
    low &= MAX_INT_NBR;    // Trim extra bits.
    ptrProd->x = low;
    ptrProd++;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    if (low < HALF_INT_RANGE)
    {
      dAccumulator = floor((dAccumulator + (double)FOURTH_INT_RANGE)*dInvRangeLimb);
    }
    else
    {
      dAccumulator = floor((dAccumulator - (double)FOURTH_INT_RANGE)*dInvRangeLimb);
    }
    low = (int)(dAccumulator - floor(dAccumulator * dInvRangeLimb) * dRangeLimb);
  }
  for (; i < (2*nbrLen); i++)
  {
    for (j = i-nbrLen+1; j < nbrLen; j++)
    {
      factor1 = (pFactor1 + j)->x;
      factor2 = (pFactor2 + i - j)->x;
      low += factor1 * factor2;
      dAccumulator += (double)factor1 * (double)factor2;
    }
    low &= MAX_INT_NBR;    // Trim extra bits.
    ptrProd->x = low;
    ptrProd++;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    if (low < HALF_INT_RANGE)
    {
      dAccumulator = floor((dAccumulator + (double)FOURTH_INT_RANGE)*dInvRangeLimb);
    }
    else
    {
      dAccumulator = floor((dAccumulator - (double)FOURTH_INT_RANGE)*dInvRangeLimb);
    }
    low = (int)(dAccumulator - floor(dAccumulator * dInvRangeLimb) * dRangeLimb);
  }
  ptrProd->x = low;
  (ptrProd + 1)->x = (int)floor(dAccumulator / dRangeLimb);
}

void IntToBigNbr(int value, limb *bigNbr, int nbrLength)
{
  limb* ptrBigNbr = bigNbr;
  int signExtended;
  if (value >= 0)
  {     // value is positive.
    ptrBigNbr->x = value;
    signExtended = 0;
  }
  else
  {     // value is negative.
    ptrBigNbr->x = value & MAX_INT_NBR;
    signExtended = MAX_INT_NBR;
  }
  for (int index = 1; index < nbrLength; index++)
  {
    ptrBigNbr++;
    ptrBigNbr->x = signExtended;
  }
}

int BigNbrToBigInt(const BigInteger *pBigNbr, limb *pBigInt)
{
  int nbrLenBigNbr = pBigNbr->nbrLimbs;
  int lenBigNbrBytes = nbrLenBigNbr * (int)sizeof(limb);
  (void)memcpy(pBigInt, pBigNbr->limbs, lenBigNbrBytes);
  return nbrLenBigNbr;
}

void BigIntToBigNbr(BigInteger *pBigNbr, const limb *pBigInt, int nbrLenBigInt)
{
  int nbrLimbs;
  int lenBigIntBytes = nbrLenBigInt * (int)sizeof(limb);
  (void)memcpy(pBigNbr->limbs, pBigInt, lenBigIntBytes);
  const limb *ptrLimb = pBigNbr->limbs + nbrLenBigInt;
  nbrLimbs = nbrLenBigInt;
  do
  {
    ptrLimb--;
    if (ptrLimb->x != 0)
    {
      break;
    }
    nbrLimbs--;
  } while (nbrLimbs > 1);
  pBigNbr->nbrLimbs = nbrLimbs;
}

void GcdBigNbr(const limb *pNbr1, const limb *pNbr2, limb *pGcd, int nbrLen)
{
  int lenBytes;
  BigIntToBigNbr(&BigInt1, pNbr1, nbrLen);
  BigIntToBigNbr(&BigInt2, pNbr2, nbrLen);
  BigIntGcd(&BigInt1, &BigInt2, &BigGcd);
  lenBytes = NumberLength * (int)sizeof(int);
  (void)memset(pGcd, 0, lenBytes);
  (void)BigNbrToBigInt(&BigGcd, pGcd);
}

void AdjustBigIntModN(limb *Nbr, const limb *Mod, int nbrLen)
{
  AdjustModN(Nbr, Mod, nbrLen);
}

void MultBigNbrModN(const limb *Nbr1, limb *Nbr2, limb *Prod, const limb *Mod, int nbrLen)
{
  int i = nbrLen;
  limb arr[MAX_LIMBS_SIQS];
  int lenBytes;

  if ((i >= 2) && ((Mod + i - 1)->x == 0))
  {
    i--;
  }
  (Nbr2+i)->x = 0;
  lenBytes = nbrLen * (int)sizeof(*Prod);
  (void)memset(Prod, 0, lenBytes);
  do
  {
    i--;
    int Nbr = (Nbr1 + i)->x;
    int j = nbrLen;
    do
    {
      *(Prod+j) = *(Prod + j - 1);
      j--;
    } while (j > 0);
    Prod[0].x = 0;
    AdjustBigIntModN(Prod, Mod, nbrLen);
    MultBigNbrByIntModN(Nbr2, Nbr, arr, Mod, nbrLen);
    AddBigIntModN(arr, Prod, Prod, Mod, nbrLen);
  } while (i > 0);
}

void MultBigNbrByIntModN(limb *Nbr1, int Nbr2, limb *Prod, const limb *Mod, int nbrLen)
{
  int nbrLength = nbrLen;
  if ((nbrLength >= 2) && ((Mod + nbrLength - 1)->x == 0))
  {
    nbrLength--;
  }
  (Nbr1+nbrLength)->x = 0;
  modmultIntExtended(Nbr1, Nbr2, Prod, Mod, nbrLength);
}

int intDoubleModPow(int NbrMod, int Expon, int currentPrime)
{
  int exponent = Expon;
  double dPower = 1;
  double dSquare = NbrMod;
  double dModulus = currentPrime;
  while (exponent != 0)
  {
    if ((exponent & 1) == 1)
    {
      dPower *= dSquare;
      dPower -= floor(dPower / dModulus)*dModulus;
    }
    dSquare *= dSquare;
    dSquare -= floor(dSquare / dModulus)*dModulus;
    exponent >>= 1;
  }
  return (int)dPower;
}

void ModInvBigInt(const limb *num, limb *inv, const limb *mod, int nbrLenBigInt)
{
  int NumberLengthBigInt;
  int NumberLengthBak = NumberLength;
  int lenBytes = nbrLenBigInt * (int)sizeof(limb);
  (void)memset(inv, 0, lenBytes);
  NumberLength = nbrLenBigInt;
  while (NumberLength > 1)
  {
    if ((mod + NumberLength - 1)->x != 0)
    {
      break;
    }
    NumberLength--;
  }
  lenBytes = NumberLength * (int)sizeof(limb);
  (void)memcpy(TestNbr, mod, lenBytes);
  TestNbr[NumberLength].x = 0;
  GetMontgomeryParms(NumberLength);
  BigIntToBigNbr(&Denominator, num, NumberLength);
  BigIntToBigNbr(&Modulus, mod, NumberLength);
  Numerator.sign = SIGN_POSITIVE;
  Numerator.nbrLimbs = 1;
  Numerator.limbs[0].x = 1;    // Numerator <- 1.
  BigIntModularDivision(&Numerator, &Denominator, &Modulus, &Quotient);
  NumberLengthBigInt = BigNbrToBigInt(&Quotient, inv);
  NumberLength = NumberLengthBak;
  if (NumberLengthBigInt < NumberLength)
  {
    lenBytes = (NumberLength - NumberLengthBigInt) * (int)sizeof(int);
    (void)memset(inv + NumberLengthBigInt, 0, lenBytes);
  }
}

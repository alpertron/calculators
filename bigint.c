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

int NbrBak[MAX_LIMBS_SIQS];
void ChSignBigNbr(int *nbr, int length)
{
  int carry = 0;
  const int *ptrEndNbr = nbr + length;
  for (int *ptrNbr = nbr; ptrNbr < ptrEndNbr; ptrNbr++)
  {
    carry -= *ptrNbr;
    *ptrNbr = carry & MAX_INT_NBR;
    carry >>= BITS_PER_INT_GROUP;
  }
}

void ChSignBigNbrB(int *nbr, int length)
{
  int carry = 0;
  int *ptrNbr = nbr;
  const int *ptrEndNbr = nbr + length - 1;
  for (ptrNbr = nbr; ptrNbr < ptrEndNbr; ptrNbr++)
  {
    carry -= *ptrNbr;
    *ptrNbr = carry & MAX_INT_NBR;
    carry >>= BITS_PER_INT_GROUP;
  }
  *ptrNbr = carry - *ptrNbr;
}

void AddBigNbr(const int *pNbr1, const int *pNbr2, int *pSum, int nbrLen)
{
  unsigned int carry = 0;
  const int *ptrNbr1 = pNbr1;
  const int *ptrNbr2 = pNbr2;
  const int *ptrEndSum = pSum + nbrLen;
  for (int *ptrSum = pSum; ptrSum < ptrEndSum; ptrSum++)
  {
    carry = (carry >> BITS_PER_INT_GROUP) + (unsigned int)*ptrNbr1 + (unsigned int)*ptrNbr2;
    *ptrSum = (int)(carry & MAX_INT_NBR);
    ptrNbr1++;
    ptrNbr2++;
  }
}

void SubtractBigNbr(const int *pNbr1, const int *pNbr2, int *pDiff, int nbrLen)
{
  int borrow = 0;
  const int *ptrNbr1 = pNbr1;
  const int *ptrNbr2 = pNbr2;
  const int *ptrEndDiff = pDiff + nbrLen;
  for (int *ptrDiff = pDiff; ptrDiff < ptrEndDiff; ptrDiff++)
  {
    borrow = (borrow >> BITS_PER_INT_GROUP) + *ptrNbr1 - *ptrNbr2;
    *ptrDiff = borrow & MAX_INT_NBR;
    ptrNbr1++;
    ptrNbr2++;
  }
}

void AddBigNbrB(const int *pNbr1, const int *pNbr2, int *pSum, int nbrLen)
{
  unsigned int carry = 0;
  const int *ptrNbr1 = pNbr1;
  const int *ptrNbr2 = pNbr2;
  int *ptrSum;
  const int* ptrEndSum = pSum + nbrLen - 1;
  for (ptrSum = pSum; ptrSum < ptrEndSum; ptrSum++)
  {
    carry = (carry >> BITS_PER_INT_GROUP) + (unsigned int)*ptrNbr1 + (unsigned int)*ptrNbr2;
    *ptrSum = (int)(carry & MAX_INT_NBR);
    ptrNbr1++;
    ptrNbr2++;
  }
  carry = (carry >> BITS_PER_INT_GROUP) + (unsigned int)*ptrNbr1 + (unsigned int)*ptrNbr2;
  *ptrSum = (int)carry;
}

void SubtractBigNbrB(const int *pNbr1, const int *pNbr2, int *pDiff, int nbrLen)
{
  int borrow = 0;
  const int *ptrNbr1 = pNbr1;
  const int *ptrNbr2 = pNbr2;
  int *ptrDiff;
  const int* ptrEndDiff = pDiff + nbrLen - 1;
  for (ptrDiff = pDiff; ptrDiff < ptrEndDiff; ptrDiff++)
  {
    borrow = (borrow >> BITS_PER_INT_GROUP) + *ptrNbr1 - *ptrNbr2;
    *ptrDiff = borrow & MAX_INT_NBR;
    ptrNbr1++;
    ptrNbr2++;
  }
  borrow = (borrow >> BITS_PER_INT_GROUP) + *ptrNbr1 - *ptrNbr2;
  *ptrDiff = borrow;
}

void AddBigIntModN(const int *pNbr1, const int *pNbr2, int *pSum, const int *pMod, int nbrLen)
{
  int borrow = 0;
  unsigned int carry = 0;
  int i;
  for (i = 0; i < nbrLen; i++)
  {
    carry = (carry >> BITS_PER_INT_GROUP) + (unsigned int)*pNbr1 + (unsigned int)*pNbr2;
    *pSum = (int)(carry & MAX_INT_NBR);
    pNbr1++;
    pNbr2++;
    pSum++;
  }
  carry >>= BITS_PER_INT_GROUP;
  pSum -= nbrLen;
  for (i = 0; i < nbrLen; i++)
  {
    borrow = (borrow >> BITS_PER_INT_GROUP) + *pSum - *pMod;
    *pSum = borrow & MAX_INT_NBR;
    pMod++;
    pSum++;
  }
  borrow >>= BITS_PER_INT_GROUP;
  if (borrow + (int)carry != 0)
  {    // Sum is less than zero. Add Mod again.
    pSum -= nbrLen;
    pMod -= nbrLen;
    carry = 0;
    for (i = 0; i < nbrLen; i++)
    {
      carry = (carry >> BITS_PER_INT_GROUP) + (unsigned int)*pSum + (unsigned int)*pMod;
      *pSum = (int)(carry & MAX_INT_NBR);
      pMod++;
      pSum++;
    }
  }
}

void SubtractBigNbrModN(const int *pNbr1, const int *pNbr2, int *pDiff, const int *pMod, int nbrLen)
{
  int borrow = 0;
  int i;
  for (i = 0; i < nbrLen; i++)
  {
    borrow = (borrow >> BITS_PER_INT_GROUP) + *pNbr1 - *pNbr2;
    *pDiff = borrow & MAX_INT_NBR;
    pNbr1++;
    pNbr2++;
    pDiff++;
  }
  if (borrow != 0)
  {
    unsigned int carry = 0;
    pDiff -= nbrLen;
    for (i = 0; i < nbrLen; i++)
    {
      carry = (carry >> BITS_PER_INT_GROUP) + (unsigned int)*pDiff + (unsigned int)*pMod;
      *pDiff = (int)(carry & MAX_INT_NBR);
      pMod++;
      pDiff++;
    }
  }
}

void MultBigNbrByInt(const int *bigFactor, int factor, int *bigProd, int nbrLen)
{
  int *bigProduct = bigProd;
  double dFactor;
  double dVal = 1 / (double)(1U<<BITS_PER_INT_GROUP);
  int factorPositive = 1;
  int carry;
  if (factor < 0)
  {     // If factor is negative, indicate it and compute its absolute value.
    factorPositive = 0;
    factor = -factor;
  }
  dFactor = (double)factor;
  carry = 0;
  for (int ctr = 0; ctr < nbrLen; ctr++)
  {
    int low = (*bigFactor * factor + carry) & MAX_INT_NBR;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    if (low < HALF_INT_RANGE)
    {
      carry = (int)floor(((double)*bigFactor * dFactor + (double)carry + HALF_INT_RANGE/2)*dVal);
    }
    else
    {
      carry = (int)floor(((double)*bigFactor * dFactor + (double)carry - HALF_INT_RANGE/2)*dVal);
    }
    *bigProduct = low;
    bigProduct++;
    bigFactor++;
  }
  if (factorPositive == 0)
  {         // If factor is negative, change sign of product.
    ChSignBigNbr(bigProd, nbrLen);
  }
}

void MultBigNbrByIntB(const int *bigFactor, int factor, int *bigProd, int nbrLen)
{
  int *bigProduct = bigProd;
  double dFactor;
  double dVal = 1 / (double)(1U << BITS_PER_INT_GROUP);
  int factorPositive = 1;
  int carry;
  if (factor < 0)
  {     // If factor is negative, indicate it and compute its absolute value.
    factorPositive = 0;
    factor = -factor;
  }
  dFactor = (double)factor;
  carry = 0;
  for (int ctr = 0; ctr < nbrLen-1; ctr++)
  {
    int low = (*bigFactor * factor + carry) & MAX_INT_NBR;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    if (low < HALF_INT_RANGE)
    {
      carry = (int)floor(((double)*bigFactor * dFactor + (double)carry + HALF_INT_RANGE / 2)*dVal);
    }
    else
    {
      carry = (int)floor(((double)*bigFactor * dFactor + (double)carry - HALF_INT_RANGE / 2)*dVal);
    }
    *bigProduct = low;
    bigProduct++;
    bigFactor++;
  }
  *bigProduct = *bigFactor * factor + carry;
  if (factorPositive == 0)
  {         // If factor is negative, change sign of product.
    ChSignBigNbrB(bigProd, nbrLen);
  }
}

void DivBigNbrByInt(const int *pDividend, int divisor, int *pQuotient, int nbrLen)
{
  int remainder = 0;
  double dDivisor = (double)divisor;
  double dLimb = 0x80000000;
  pDividend += nbrLen - 1;
  pQuotient += nbrLen - 1;
  for (int ctr = nbrLen - 1; ctr >= 0; ctr--)
  {
    int dividend = (remainder << BITS_PER_INT_GROUP) + *pDividend;
    double dDividend = (double)remainder * dLimb + *pDividend;
    // quotient has correct value or 1 more.
    int quotient = (unsigned int)(dDividend / dDivisor + 0.5);
    remainder = dividend - quotient * divisor;
    if ((unsigned int)remainder >= (unsigned int)divisor)
    {     // remainder not in range 0 <= remainder < divisor. Adjust.
      quotient--;
      remainder += divisor;
    }
    *pQuotient = quotient;
    pQuotient--;
    pDividend--;
  }
}

int RemDivBigNbrByInt(const int *pDividend, int divisor, int nbrLen)
{
  int remainder = 0;
  double dDivisor = (double)divisor;
  double dLimb = 0x80000000;
  pDividend += nbrLen - 1;
  for (int ctr = nbrLen - 1; ctr >= 0; ctr--)
  {
    unsigned int dividend = (remainder << BITS_PER_INT_GROUP) + *pDividend;
    double dDividend = (double)remainder * dLimb + *pDividend;
         // quotient has correct value or 1 more.
    unsigned int quotient = (unsigned int)(dDividend / dDivisor + 0.5);
    remainder = dividend - quotient * divisor;
    if ((unsigned int)remainder >= (unsigned int)divisor)
    {     // remainder not in range 0 <= remainder < divisor. Adjust.
      quotient--;
      remainder += divisor;
    }
    pDividend--;
  }
  return remainder;
}

void MultBigNbr(const int *pFactor1, const int *pFactor2, int *pProd, int nbrLen)
{
  double dRangeLimb = (double)(1U << BITS_PER_INT_GROUP);
  double dInvRangeLimb = 1 / dRangeLimb;
  int low = 0;
  int factor1;
  int factor2;
  double dAccumulator = 0;
  for (int i = 0; i < nbrLen; i++)
  {
    for (int j = 0; j <= i; j++)
    {
      factor1 = *(pFactor1 + j);
      factor2 = *(pFactor2 + i - j);
      low += factor1*factor2;
      dAccumulator += (double)factor1 * (double)factor2;
    }
    low &= MAX_INT_NBR;    // Trim extra bits.
    *pProd = low;
    pProd++;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    if (low < HALF_INT_RANGE)
    {
      dAccumulator = floor((dAccumulator + HALF_INT_RANGE/2)*dInvRangeLimb);
    }
    else
    {
      dAccumulator = floor((dAccumulator - HALF_INT_RANGE/2)*dInvRangeLimb);
    }
    low = (int)(dAccumulator - floor(dAccumulator * dInvRangeLimb) * dRangeLimb);
  }
  *pProd = low;
  *(pProd+1) = (int)floor(dAccumulator/dRangeLimb);
}

// On input: pFactor1 and pFactor2: pointers to factors.
//           pProd: pointer to product (length = 2*nbrLen)
//           nbrLen: number of limbs of factors.
void MultBigNbrComplete(const int *pFactor1, const int *pFactor2, int *pProd, int nbrLen)
{
  double dRangeLimb = (double)(1U << BITS_PER_INT_GROUP);
  double dInvRangeLimb = 1 / dRangeLimb;
  int low = 0;
  int i;
  int j;
  int factor1;
  int factor2;
  double dAccumulator = 0;
  for (i = 0; i < nbrLen; i++)
  {
    for (j = 0; j <= i; j++)
    {
      factor1 = *(pFactor1 + j);
      factor2 = *(pFactor2 + i - j);
      low += factor1 * factor2;
      dAccumulator += (double)factor1 * (double)factor2;
    }
    low &= MAX_INT_NBR;    // Trim extra bits.
    *pProd = low;
    pProd++;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    if (low < HALF_INT_RANGE)
    {
      dAccumulator = floor((dAccumulator + HALF_INT_RANGE / 2)*dInvRangeLimb);
    }
    else
    {
      dAccumulator = floor((dAccumulator - HALF_INT_RANGE / 2)*dInvRangeLimb);
    }
    low = (int)(dAccumulator - floor(dAccumulator * dInvRangeLimb) * dRangeLimb);
  }
  for (; i < 2*nbrLen; i++)
  {
    for (j = i-nbrLen+1; j < nbrLen; j++)
    {
      factor1 = *(pFactor1 + j);
      factor2 = *(pFactor2 + i - j);
      low += factor1 * factor2;
      dAccumulator += (double)factor1 * (double)factor2;
    }
    low &= MAX_INT_NBR;    // Trim extra bits.
    *pProd = low;
    pProd++;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    if (low < HALF_INT_RANGE)
    {
      dAccumulator = floor((dAccumulator + HALF_INT_RANGE / 2)*dInvRangeLimb);
    }
    else
    {
      dAccumulator = floor((dAccumulator - HALF_INT_RANGE / 2)*dInvRangeLimb);
    }
    low = (int)(dAccumulator - floor(dAccumulator * dInvRangeLimb) * dRangeLimb);
  }
  *pProd = low;
  *(pProd + 1) = (int)floor(dAccumulator / dRangeLimb);
}

void IntToBigNbr(int value, int *bigNbr, int nbrLength)
{
  if (value >= 0)
  {     // value is positive.
    *bigNbr = value;
    value = 0;
  }
  else
  {     // value is negative.
    *bigNbr = value & MAX_INT_NBR;
    value = MAX_INT_NBR;
  }
  for (; nbrLength > 1; nbrLength--)
  {
    bigNbr++;
    *bigNbr = value;
  }
}

int BigNbrToBigInt(const BigInteger *pBigNbr, int *pBigInt)
{
  int nbrLenBigNbr = pBigNbr->nbrLimbs;
  int *ptrBigInt = pBigInt;
  const limb *ptrLimb = pBigNbr->limbs;
  for (int ctr = 0; ctr < nbrLenBigNbr; ctr++)
  {
    *ptrBigInt = ptrLimb->x;
    ptrBigInt++;
    ptrLimb++;
  }
  return nbrLenBigNbr;
}

void BigIntToBigNbr(BigInteger *pBigNbr, const int *pBigInt, int nbrLenBigInt)
{
  int nbrLimbs;
  const int *ptrBigInt = pBigInt;
  limb *ptrLimb = pBigNbr->limbs;
  pBigNbr->sign = SIGN_POSITIVE;
  for (int ctr = 0; ctr < nbrLenBigInt; ctr++)
  {
    ptrLimb->x = *ptrBigInt;
    ptrLimb++;
    ptrBigInt++;
  }
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

void GcdBigNbr(const int *pNbr1, const int *pNbr2, int *pGcd, int nbrLen)
{
  BigIntToBigNbr(&BigInt1, pNbr1, nbrLen);
  BigIntToBigNbr(&BigInt2, pNbr2, nbrLen);
  BigIntGcd(&BigInt1, &BigInt2, &BigGcd);
  (void)memset(pGcd, 0, NumberLength * sizeof(int));
  BigNbrToBigInt(&BigGcd, pGcd);
}

void AdjustBigIntModN(int *Nbr, int *Mod, int nbrLen)
{
  AdjustModN((limb *)Nbr, (limb *)Mod, nbrLen);
}

void MultBigNbrModN(const int *Nbr1, int *Nbr2, int *Prod, int *Mod, int nbrLen)
{
  int i;
  int arr[MAX_LIMBS_SIQS];

  if ((nbrLen >= 2) && (*(Mod + nbrLen - 1) == 0))
  {
    nbrLen--;
  }
  *(Nbr2+nbrLen) = 0;
  (void)memset(Prod, 0, nbrLen * sizeof(*Prod));
  i = nbrLen;
  do
  {
    i--;
    int Nbr = *(Nbr1 + i);
    int j = nbrLen;
    do
    {
      *(Prod+j) = *(Prod + j - 1);
      j--;
    } while (j > 0);
    Prod[0] = 0;
    AdjustBigIntModN(Prod, Mod, nbrLen);
    MultBigNbrByIntModN(Nbr2, Nbr, arr, Mod, nbrLen);
    AddBigIntModN(arr, Prod, Prod, Mod, nbrLen);
  } while (i > 0);
}

void MultBigNbrByIntModN(int *Nbr1, int Nbr2, int *Prod, int *Mod, int nbrLen)
{
  if ((nbrLen >= 2) && (*(Mod + nbrLen - 1) == 0))
  {
    nbrLen--;
  }
  *(Nbr1+nbrLen) = 0;
  modmultIntExtended((limb *)Nbr1, Nbr2, (limb *)Prod, (limb *)Mod, nbrLen);
}

int intDoubleModPow(int NbrMod, int Expon, int currentPrime)
{
  double dPower = 1;
  double dSquare = NbrMod;
  double dModulus = currentPrime;
  while (Expon != 0)
  {
    if ((Expon & 1) == 1)
    {
      dPower *= dSquare;
      dPower -= floor(dPower / dModulus)*dModulus;
    }
    dSquare *= dSquare;
    dSquare -= floor(dSquare / dModulus)*dModulus;
    Expon >>= 1;
  }
  return (int)dPower;
}

void ModInvBigInt(const int *num, int *inv, const int *mod, int nbrLenBigInt)
{
  int NumberLengthBigInt;
  int NumberLengthBak = NumberLength;
  (void)memset(inv, 0, nbrLenBigInt*sizeof(int));
  while (nbrLenBigInt > 1)
  {
    if (*(mod + nbrLenBigInt - 1) != 0)
    {
      break;
    }
    nbrLenBigInt--;
  }
  NumberLength = nbrLenBigInt;
  (void)memcpy(TestNbr, mod, NumberLength * sizeof(limb));
  TestNbr[NumberLength].x = 0;
  GetMontgomeryParms(NumberLength);
  BigIntToBigNbr(&Denominator, num, nbrLenBigInt);
  BigIntToBigNbr(&Modulus, mod, nbrLenBigInt);
  Numerator.sign = SIGN_POSITIVE;
  Numerator.nbrLimbs = 1;
  Numerator.limbs[0].x = 1;    // Numerator <- 1.
  BigIntModularDivision(&Numerator, &Denominator, &Modulus, &Quotient);
  NumberLengthBigInt = BigNbrToBigInt(&Quotient, inv);
  NumberLength = NumberLengthBak;
  if (NumberLengthBigInt < NumberLength)
  {
    (void)memset(inv + NumberLengthBigInt, 0, (NumberLength - NumberLengthBigInt) * sizeof(int));
  }
}

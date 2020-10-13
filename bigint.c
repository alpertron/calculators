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
#include <string.h>
#include <math.h>
#include "bignbr.h"

#define MAX_LIMBS_SIQS 15
static BigInteger Numerator, Denominator, Modulus, Quotient;
static BigInteger BigInt1, BigInt2, BigGcd;

int NbrBak[MAX_LIMBS_SIQS];
void ChSignBigNbr(int *nbr, int length)
{
  int carry = 0;
  int ctr;
  for (ctr = 0; ctr < length; ctr++)
  {
    carry -= *nbr;
    *nbr++ = carry & MAX_INT_NBR;
    carry >>= BITS_PER_INT_GROUP;
  }
}

void ChSignBigNbrB(int *nbr, int length)
{
  int carry = 0;
  int ctr;
  for (ctr = 0; ctr < length-1; ctr++)
  {
    carry -= *nbr;
    *nbr++ = carry & MAX_INT_NBR;
    carry >>= BITS_PER_INT_GROUP;
  }
  *nbr = carry - *nbr;
}

void AddBigNbr(int *pNbr1, int *pNbr2, int *pSum, int nbrLen)
{
  unsigned int carry = 0;
  int i;
  for (i = 0; i < nbrLen; i++)
  {
    carry = (carry >> BITS_PER_INT_GROUP) + (unsigned int)*pNbr1++ + (unsigned int)*pNbr2++;
    *pSum++ = (int)(carry & MAX_INT_NBR);
  }
}

void SubtractBigNbr(int *pNbr1, int *pNbr2, int *pDiff, int nbrLen)
{
  int borrow = 0;
  int i;
  for (i = 0; i < nbrLen; i++)
  {
    borrow = (borrow >> BITS_PER_INT_GROUP) + *pNbr1++ - *pNbr2++;
    *pDiff++ = borrow & MAX_INT_NBR;
  }
}

void AddBigNbrB(int *pNbr1, int *pNbr2, int *pSum, int nbrLen)
{
  unsigned int carry = 0;
  int i;
  for (i = 0; i < nbrLen-1; i++)
  {
    carry = (carry >> BITS_PER_INT_GROUP) + (unsigned int)*pNbr1++ + (unsigned int)*pNbr2++;
    *pSum++ = (int)(carry & MAX_INT_NBR);
  }
  carry = (carry >> BITS_PER_INT_GROUP) + (unsigned int)*pNbr1++ + (unsigned int)*pNbr2++;
  *pSum = (int)carry;
}

void SubtractBigNbrB(int *pNbr1, int *pNbr2, int *pDiff, int nbrLen)
{
  int borrow = 0;
  int i;
  for (i = 0; i < nbrLen-1; i++)
  {
    borrow = (borrow >> BITS_PER_INT_GROUP) + *pNbr1++ - *pNbr2++;
    *pDiff++ = borrow & MAX_INT_NBR;
  }
  borrow = (borrow >> BITS_PER_INT_GROUP) + *pNbr1++ - *pNbr2++;
  *pDiff++ = borrow;
}

void AddBigIntModN(int *pNbr1, int *pNbr2, int *pSum, int *pMod, int nbrLen)
{
  int borrow = 0;
  unsigned int carry = 0;
  int i;
  for (i = 0; i < nbrLen; i++)
  {
    carry = (carry >> BITS_PER_INT_GROUP) + (unsigned int)*pNbr1++ + (unsigned int)*pNbr2++;
    *pSum++ = (int)(carry & MAX_INT_NBR);
  }
  carry >>= BITS_PER_INT_GROUP;
  pSum -= nbrLen;
  for (i = 0; i < nbrLen; i++)
  {
    borrow = (borrow >> BITS_PER_INT_GROUP) + *pSum - *pMod++;
    *pSum++ = borrow & MAX_INT_NBR;
  }
  borrow >>= BITS_PER_INT_GROUP;
  if (borrow + (int)carry != 0)
  {    // Sum is less than zero. Add Mod again.
    pSum -= nbrLen;
    pMod -= nbrLen;
    carry = 0;
    for (i = 0; i < nbrLen; i++)
    {
      carry = (carry >> BITS_PER_INT_GROUP) + (unsigned int)*pSum + (unsigned int)*pMod++;
      *pSum++ = (int)(carry & MAX_INT_NBR);
    }
  }
}

void SubtractBigNbrModN(int *pNbr1, int *pNbr2, int *pDiff, int *pMod, int nbrLen)
{
  int borrow = 0;
  int i;
  for (i = 0; i < nbrLen; i++)
  {
    borrow = (borrow >> BITS_PER_INT_GROUP) + *pNbr1++ - *pNbr2++;
    *pDiff++ = borrow & MAX_INT_NBR;
  }
  if (borrow != 0)
  {
    unsigned int carry = 0;
    pDiff -= nbrLen;
    for (i = 0; i < nbrLen; i++)
    {
      carry = (carry >> BITS_PER_INT_GROUP) + (unsigned int)*pDiff + (unsigned int)*pMod++;
      *pDiff++ = (int)(carry & MAX_INT_NBR);
    }
  }
}

void MultBigNbrByInt(int *bigFactor, int factor, int *bigProd, int nbrLen)
{
  int *bigProduct = bigProd;
  double dFactor;
  double dVal = 1 / (double)(1U<<BITS_PER_INT_GROUP);
  int factorPositive = 1;
  int ctr, carry;
  if (factor < 0)
  {     // If factor is negative, indicate it and compute its absolute value.
    factorPositive = 0;
    factor = -factor;
  }
  dFactor = (double)factor;
  carry = 0;
  for (ctr = 0; ctr < nbrLen; ctr++)
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
    *bigProduct++ = low;
    bigFactor++;
  }
  if (factorPositive == 0)
  {         // If factor is negative, change sign of product.
    ChSignBigNbr(bigProd, nbrLen);
  }
}

void MultBigNbrByIntB(int *bigFactor, int factor, int *bigProd, int nbrLen)
{
  int *bigProduct = bigProd;
  double dFactor;
  double dVal = 1 / (double)(1U << BITS_PER_INT_GROUP);
  int factorPositive = 1;
  int ctr, carry;
  if (factor < 0)
  {     // If factor is negative, indicate it and compute its absolute value.
    factorPositive = 0;
    factor = -factor;
  }
  dFactor = (double)factor;
  carry = 0;
  for (ctr = 0; ctr < nbrLen-1; ctr++)
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
    *bigProduct++ = low;
    bigFactor++;
  }
  *bigProduct = *bigFactor * factor + carry;
  if (factorPositive == 0)
  {         // If factor is negative, change sign of product.
    ChSignBigNbrB(bigProd, nbrLen);
  }
}

void DivBigNbrByInt(int *pDividend, int divisor, int *pQuotient, int nbrLen)
{
  int ctr;
  int remainder = 0;
  double dDivisor = (double)divisor;
  double dLimb = 0x80000000;
  pDividend += nbrLen - 1;
  pQuotient += nbrLen - 1;
  for (ctr = nbrLen - 1; ctr >= 0; ctr--)
  {
    double dDividend;
    int quotient, dividend;
    dividend = (remainder << BITS_PER_INT_GROUP) + *pDividend;
    dDividend = (double)remainder * dLimb + *pDividend;
    // quotient has correct value or 1 more.
    quotient = (unsigned int)(dDividend / dDivisor + 0.5);
    remainder = dividend - quotient * divisor;
    if ((unsigned int)remainder >= (unsigned int)divisor)
    {     // remainder not in range 0 <= remainder < divisor. Adjust.
      quotient--;
      remainder += divisor;
    }
    *pQuotient-- = quotient;
    pDividend--;
  }
}

int RemDivBigNbrByInt(int *pDividend, int divisor, int nbrLen)
{
  int ctr;
  int remainder = 0;
  double dDivisor = (double)divisor;
  double dLimb = 0x80000000;
  pDividend += nbrLen - 1;
  for (ctr = nbrLen - 1; ctr >= 0; ctr--)
  {
    unsigned int quotient, dividend;
    double dDividend;
    dividend = (remainder << BITS_PER_INT_GROUP) + *pDividend;
    dDividend = (double)remainder * dLimb + *pDividend;
         // quotient has correct value or 1 more.
    quotient = (unsigned int)(dDividend / dDivisor + 0.5);
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

void MultBigNbr(int *pFactor1, int *pFactor2, int *pProd, int nbrLen)
{
  double dRangeLimb = (double)(1U << BITS_PER_INT_GROUP);
  double dInvRangeLimb = 1 / dRangeLimb;
  int low = 0;
  int i, j;
  int factor1, factor2;
  double dAccumulator = 0;
  for (i = 0; i < nbrLen; i++)
  {
    for (j = 0; j <= i; j++)
    {
      factor1 = *(pFactor1 + j);
      factor2 = *(pFactor2 + i - j);
      low += factor1*factor2;
      dAccumulator += (double)factor1 * (double)factor2;
    }
    low &= MAX_INT_NBR;    // Trim extra bits.
    *pProd++ = low;
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
  int i, j;
  int factor1, factor2;
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
    *pProd++ = low;
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
    *pProd++ = low;
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
    *++bigNbr = value;
  }
}

int BigNbrToBigInt(BigInteger *pBigNbr, int *pBigInt)
{
  int ctr;
  int nbrLenBigNbr = pBigNbr->nbrLimbs;
  int *ptrBigInt = pBigInt;
  limb *ptrLimb = pBigNbr->limbs;
  for (ctr = 0; ctr < nbrLenBigNbr; ctr++)
  {
    *ptrBigInt++ = ptrLimb++->x;
  }
  return nbrLenBigNbr;
}

void BigIntToBigNbr(BigInteger *pBigNbr, int *pBigInt, int nbrLenBigInt)
{
  int ctr, nbrLimbs;
  int *ptrBigInt = pBigInt;
  limb *ptrLimb = pBigNbr->limbs;
  pBigNbr->sign = SIGN_POSITIVE;
  for (ctr = 0; ctr < nbrLenBigInt; ctr++)
  {
    ptrLimb++->x = *ptrBigInt++;
  }
  nbrLimbs = nbrLenBigInt;
  do
  {
    if ((--ptrLimb)->x != 0)
    {
      break;
    }
  } while (--nbrLimbs > 1);
  pBigNbr->nbrLimbs = nbrLimbs;
}

void GcdBigNbr(int *pNbr1, int *pNbr2, int *pGcd, int nbrLen)
{
  BigIntToBigNbr(&BigInt1, pNbr1, nbrLen);
  BigIntToBigNbr(&BigInt2, pNbr2, nbrLen);
  BigIntGcd(&BigInt1, &BigInt2, &BigGcd);
  memset(pGcd, 0, NumberLength * sizeof(int));
  BigNbrToBigInt(&BigGcd, pGcd);
}

void AdjustBigIntModN(int *Nbr, int *Mod, int nbrLen)
{
  AdjustModN((limb *)Nbr, (limb *)Mod, nbrLen);
}

void MultBigNbrModN(int *Nbr1, int *Nbr2, int *Prod, int *Mod, int nbrLen)
{
  int i;
  int arr[MAX_LIMBS_SIQS];

  if (nbrLen >= 2 && *(Mod + nbrLen - 1) == 0)
  {
    nbrLen--;
  }
  *(Nbr2+nbrLen) = 0;
  memset(Prod, 0, nbrLen * sizeof(*Prod));
  i = nbrLen;
  do
  {
    int Nbr = *(Nbr1 + --i);
    int j = nbrLen;
    do
    {
      *(Prod+j) = *(Prod + j - 1);
    } while (--j > 0);
    Prod[0] = 0;
    AdjustBigIntModN(Prod, Mod, nbrLen);
    MultBigNbrByIntModN(Nbr2, Nbr, arr, Mod, nbrLen);
    AddBigIntModN(arr, Prod, Prod, Mod, nbrLen);
  } while (i > 0);
}

void MultBigNbrByIntModN(int *Nbr1, int Nbr2, int *Prod, int *Mod, int nbrLen)
{
  if (nbrLen >= 2 && *(Mod + nbrLen - 1) == 0)
  {
    nbrLen--;
  }
  *(Nbr1+nbrLen) = 0;
  modmultIntExtended((limb *)Nbr1, Nbr2, (limb *)Prod, (limb *)Mod, nbrLen);
}

int intDoubleModPow(int NbrMod, int Expon, int currentPrime)
{
  double Power = 1;
  double Square = NbrMod;
  double Modulus = currentPrime;
  while (Expon != 0)
  {
    if ((Expon & 1) == 1)
    {
      Power *= Square;
      Power -= floor(Power / Modulus)*Modulus;
    }
    Square *= Square;
    Square -= floor(Square / Modulus)*Modulus;
    Expon >>= 1;
  }
  return (int)Power;
}

void ModInvBigInt(int *num, int *inv, int *mod, int nbrLenBigInt)
{
  int NumberLengthBigInt;
  int NumberLengthBak = NumberLength;
  memset(inv, 0, nbrLenBigInt*sizeof(int));
  while (nbrLenBigInt > 1)
  {
    if (*(mod + nbrLenBigInt - 1) != 0)
    {
      break;
    }
    nbrLenBigInt--;
  }
  NumberLength = nbrLenBigInt;
  memcpy(TestNbr, mod, NumberLength * sizeof(limb));
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
    memset(inv + NumberLengthBigInt, 0, (NumberLength - NumberLengthBigInt) * sizeof(int));
  }
}

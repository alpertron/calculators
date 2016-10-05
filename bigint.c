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

#define MAX_LIMBS_SIQS 25

int NbrBak[MAX_LIMBS_SIQS];
void ChSignBigNbr(int *nbr, int length)
{
  int carry = 0;
  int ctr;
  for (ctr = 0; ctr < length; ctr++)
  {
    carry -= *nbr;
    *nbr &= -LIMB_RANGE;
    carry >>= BITS_PER_GROUP;
  }
}

void AddBigNbr(int *pNbr1, int *pNbr2, int *pSum, int nbrLen)
{
  unsigned int carry = 0;
  int i;
  for (i = 0; i < nbrLen; i++)
  {
    carry = (carry >> 31) + (unsigned int)*pNbr1++ + (unsigned int)*pNbr2++;
    *pSum++ = (int)(carry & 0x7FFFFFFF);
  }
}

void SubtractBigNbr(int *pNbr1, int *pNbr2, int *pDiff, int nbrLen)
{
  int carry = 0;
  int i;
  for (i = 0; i < nbrLen; i++)
  {
    carry = (carry >> 31) + *pNbr1++ - *pNbr2++;
    *pDiff++ = carry & 0x7FFFFFFF;
  }
}

void AddBigIntModN(int *pNbr1, int *pNbr2, int *pSum, int *pMod, int nbrLen)
{
  int borrow = 0;
  unsigned int carry = 0;
  int i;
  for (i = 0; i < nbrLen; i++)
  {
    carry = (carry >> 31) + (unsigned int)*pNbr1++ + (unsigned int)*pNbr2++;
    *pSum++ = (int)(carry & 0x7FFFFFFF);
  }
  pSum -= nbrLen;
  for (i = 0; i < nbrLen; i++)
  {
    borrow = (borrow >> 31) + *pSum - *pMod++;
    *pSum++ = borrow & 0x7FFFFFFF;
  }
  if (borrow + (int)carry != 0)
  {    // Sum is less than zero. Add Mod again.
    pSum -= nbrLen;
    pMod -= nbrLen;
    carry = 0;
    for (i = 0; i < nbrLen; i++)
    {
      carry = (carry >> 31) + (unsigned int)*pSum + (unsigned int)*pMod++;
      *pSum++ = (int)(carry & 0x7FFFFFFF);
    }
  }
}

void SubtractBigNbrModN(int *pNbr1, int *pNbr2, int *pDiff, int *pMod, int nbrLen)
{
  int borrow = 0;
  unsigned int carry;
  int i;
  for (i = 0; i < nbrLen; i++)
  {
    borrow = (borrow >> 31) + *pNbr1++ - *pNbr2++;
    *pDiff++ = borrow & 0x7FFFFFFF;
  }
  if (borrow != 0)
  {
    pDiff -= nbrLen;
    carry = 0;
    for (i = 0; i < nbrLen; i++)
    {
      carry = (carry >> 31) + (unsigned int)*pDiff + (unsigned int)*pMod++;
      *pDiff++ = (int)(carry & 0x7FFFFFFF);
    }
  }
}

void MultBigNbrByInt(int *bigFactor, int factor, int *bigProduct, int nbrLen)
{
  double dVal = 1 / 0x80000000e0;
  double dFactor = (double)factor;
  int ctr, low;
  int carry = 0;
  for (ctr = 0; ctr < nbrLen; ctr++)
  {
    low = *bigFactor * factor + carry;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    if (low < 0x40000000)
    {
      carry = (int)(((double)*bigFactor * dFactor + (double)carry + 0x20000000)*dVal);
    }
    else
    {
      carry = (int)(((double)*bigFactor * dFactor + (double)carry - 0x20000000)*dVal);
    }
    *bigProduct++ = low;
    bigFactor++;
  }
}

void DivBigNbrByInt(int *pDividend, int divisor, int *pQuotient, int nbrLen)
{
  int ctr, quotient, dividend;
  int remainder = 0;
  double dDivisor = (double)divisor;
  double dLimb = 0x80000000;
  double dDividend, dQuotient;
  pDividend += nbrLen - 1;
  pQuotient += nbrLen - 1;
  for (ctr = nbrLen - 1; ctr >= 0; ctr--)
  {
    dividend = (remainder << 31) + *pDividend;
    dDividend = (double)remainder * dLimb + *pDividend;
    dQuotient = dDividend / dDivisor + 0.5;
    quotient = (int)dQuotient;   // quotient has correct value or 1 more.
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
  int ctr, quotient, dividend;
  int remainder = 0;
  double dDivisor = (double)divisor;
  double dLimb = 0x80000000;
  double dDividend, dQuotient;
  pDividend += nbrLen - 1;
  for (ctr = nbrLen - 1; ctr >= 0; ctr--)
  {
    dividend = (remainder << 31) + *pDividend;
    dDividend = (double)remainder * dLimb + *pDividend;
    dQuotient = dDividend / dDivisor + 0.5;
    quotient = (int)dQuotient;   // quotient has correct value or 1 more.
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
  double dVal = 1 / 0x80000000e0;
  double dAccumulator;
  int low;
  int i, j;
  int factor1, factor2;
  low = 0;
  for (i = 0; i < nbrLen; i++)
  {
    dAccumulator = 0;
    for (j = 0; j <= i; j++)
    {
      factor1 = *(pFactor1 + j);
      factor2 = *(pFactor2 + i - j);
      low += factor1*factor2;
      dAccumulator += (double)factor1 * (double)factor2;
    }
    *pProd++ = low;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    if (low < 0x40000000)
    {
      dAccumulator = floor((dAccumulator + 0x20000000)*dVal);
    }
    else
    {
      dAccumulator = floor((dAccumulator - 0x20000000)*dVal);
    }
  }
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
    *bigNbr = value & 0x7FFFFFFF;
    value = 0x7FFFFFFF;
  }
  for (; nbrLength > 1; nbrLength--)
  {
    *++bigNbr = value;
  }
}

int BigNbrToBigInt(BigInteger *pBigNbr, int *pBigInt)
{
  int ctr;
  int currentValue = 0;
  int nbrLenBigNbr = pBigNbr->nbrLimbs;
  int curShift = 0;
  int *ptrBigInt = pBigInt;
  limb *ptrLimb = pBigNbr->limbs;
  for (ctr = 0; ctr < nbrLenBigNbr; ctr++)
  {
    currentValue += ptrLimb->x << curShift;
    curShift += BITS_PER_GROUP;
    if (curShift >= 31)
    {
      *ptrBigInt++ = currentValue & 0x7FFFFFFF;
      curShift -= 31;
      currentValue = ptrLimb->x >> (BITS_PER_GROUP - curShift);
    }
    ptrLimb++;
  }
  if (curShift == 0)
  {
    return (int)(ptrBigInt - pBigInt);
  }
  *ptrBigInt = currentValue;
  return (int)(ptrBigInt - pBigInt) + 1;
}

void BigIntToBigNbr(BigInteger *pBigNbr, int *pBigInt, int nbrLenBigInt)
{
  int ctr, nbrLimbs;
  int currentValue = 0;
  int curShift = 0;
  int *ptrBigInt = pBigInt;
  limb *ptrLimb = pBigNbr->limbs;
  pBigNbr->sign = SIGN_POSITIVE;
  for (ctr = 0; ctr < nbrLenBigInt; ctr++)
  {
    curShift += 31;
    while (curShift >= BITS_PER_GROUP)
    {
      currentValue |= *ptrBigInt << curShift;
      ptrLimb++->x = currentValue & MAX_VALUE_LIMB;
      curShift -= BITS_PER_GROUP;
      currentValue >>= BITS_PER_GROUP;
    }
    ptrBigInt++;
  }
  if (curShift == 0)
  {
    nbrLimbs = (int)(ptrBigInt - pBigInt);
  }
  else
  {
    nbrLimbs = (int)(ptrBigInt - pBigInt) + 1;
    ptrLimb++->x = currentValue;
  }
  for (nbrLimbs = NumberLength; nbrLimbs > 1; nbrLimbs--)
  {
    if (--ptrLimb->x != 0)
    {
      break;
    }
  }
  pBigNbr->nbrLimbs = nbrLimbs;
}

void GcdBigNbr(int *pNbr1, int *pNbr2, int *pGcd, int nbrLen)
{
  BigInteger BigInt1, BigInt2, BigGcd;
  BigIntToBigNbr(&BigInt1, pNbr1, nbrLen);
  BigIntToBigNbr(&BigInt2, pNbr2, nbrLen);
  BigIntGcd(&BigInt1, &BigInt2, &BigGcd);
  memset(pGcd, 0, NumberLength * sizeof(int));
  BigNbrToBigInt(&BigGcd, pGcd);
}

void AdjustBigIntModN(int Nbr[], int Mod[], int nbrLen)
{
  int MaxUInt = 0x7FFFFFFFL;
  int TrialQuotient;
  double dAux, dN;
  double dMaxInt = (double)(1U << 31);

  dN = (double)Mod[nbrLen - 1];
  if (nbrLen > 1)
  {
    dN += (double)Mod[nbrLen - 2] / dMaxInt;
  }
  if (nbrLen > 2)
  {
    dN += (double)Mod[nbrLen - 3] / (dMaxInt*dMaxInt);
  }
  dAux =
    (double)Nbr[nbrLen] * dMaxInt + (double)Nbr[nbrLen - 1];
  if (nbrLen > 1)
  {
    dAux += (double)Nbr[nbrLen - 2] / dMaxInt;
  }
  TrialQuotient = (long)(dAux / dN) + 3;
  if (TrialQuotient < 0)
  {     // Exceeded maximum value for TrialQuotient.
    TrialQuotient = MaxUInt;
  }
  MultBigNbrByInt(Mod, TrialQuotient, NbrBak, nbrLen + 1);
  SubtractBigNbr(Nbr, NbrBak, Nbr, nbrLen + 1);
  while (Nbr[nbrLen] != 0)
  {
    AddBigNbr(Nbr, Mod, Nbr, nbrLen + 1);
  }
}

void MultBigNbrModN(int Nbr1[], int Nbr2[], int Prod[], int Mod[], int nbrLen)
{
  int i, j;
  int Nbr;

  if (nbrLen >= 2 && Mod[nbrLen - 1] == 0 && Mod[nbrLen - 2]<0x40000000)
  {
    nbrLen--;
  }
  i = nbrLen;
  do
  {
    Prod[--i] = 0;
  } while (i > 0);
  i = nbrLen;
  do
  {
    Nbr = Nbr1[--i];
    j = nbrLen;
    do
    {
      Prod[j] = Prod[j - 1];
      j--;
    } while (j > 0);
    MultBigNbrByInt(Nbr2, Nbr, Prod, nbrLen + 1);
    AdjustBigIntModN(Prod, Mod, nbrLen);
  } while (i > 0);
}

void MultBigNbrByIntModN(int Nbr1[], int Nbr2, int Prod[], int Mod[], int nbrLen)
{
  if (nbrLen >= 2 && Mod[nbrLen - 1] == 0 && Mod[nbrLen - 2]<0x40000000)
  {
    nbrLen--;
  }
  MultBigNbrByInt(Nbr1, Nbr2, Prod, nbrLen + 1);
  AdjustBigIntModN(Prod, Mod, nbrLen);
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
  BigInteger Numerator, Denominator, Modulus, Quotient;
  BigIntToBigNbr(&Denominator, num, nbrLenBigInt);
  BigIntToBigNbr(&Modulus, mod, nbrLenBigInt);
  Numerator.sign = SIGN_POSITIVE;
  Numerator.nbrLimbs = 1;
  Numerator.limbs[0].x = 1;    // Numerator <- 1.
  BigIntModularDivision(&Numerator, &Denominator, &Modulus, &Quotient);
  BigNbrToBigInt(&Quotient, inv);
}

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
#include "expression.h"

static BigInteger Temp, Temp2, Temp3, Base, Power;
static char ProcessExpon[5000];
static char primes[5000];
extern limb Mult1[MAX_LEN];
extern limb Mult3[MAX_LEN];
extern limb Mult4[MAX_LEN];
extern int q[MAX_LEN];
extern limb TestNbr[MAX_LEN];
extern limb MontgomeryMultR1[MAX_LEN];
int groupLen = 6;

void CopyBigInt(BigInteger *pDest, BigInteger *pSrc)
{
  pDest->sign = pSrc->sign;
  pDest->nbrLimbs = pSrc->nbrLimbs;
  memcpy(pDest->limbs, pSrc->limbs, (pSrc->nbrLimbs)*sizeof(limb));
}

void AddBigInt(limb *pAddend1, limb *pAddend2, limb *pSum, int nbrLimbs)
{
  unsigned int carry = 0;
  int i;
  for (i = 0; i < nbrLimbs; i++)
  {
    carry = (carry >> BITS_PER_GROUP) + (unsigned int)(pAddend1++)->x +
                                        (unsigned int)(pAddend2++)->x;
    (pSum++)->x = (int)(carry & MAX_INT_NBR);
  }
}

void SubtractBigInt(limb *pMinuend, limb *pSubtrahend, limb *pDiff, int nbrLimbs)
{
  int borrow = 0;
  int i;
  for (i = 0; i < nbrLimbs; i++)
  {
    borrow = (borrow >> BITS_PER_INT_GROUP) + (pMinuend++)->x - (pSubtrahend++)->x;
    (pDiff++)->x = borrow & MAX_INT_NBR;
  }
}

void BigIntAdd(BigInteger *pAddend1, BigInteger *pAddend2, BigInteger *pSum)
{
  int ctr, nbrLimbs;
  limb *ptrAddend1, *ptrAddend2, *ptrSum;
  BigInteger *pTemp;
  if (pAddend1->nbrLimbs < pAddend2->nbrLimbs)
  {
    pTemp = pAddend1;
    pAddend1 = pAddend2;
    pAddend2 = pTemp;
  }           // At this moment, the absolute value of addend1 is greater than
              // or equal than the absolute value of addend2.
  else if (pAddend1->nbrLimbs == pAddend2->nbrLimbs)
  {
    for (ctr = pAddend1->nbrLimbs - 1; ctr >= 0; ctr--)
    {
      if (pAddend1->limbs[ctr].x != pAddend2->limbs[ctr].x)
      {
        break;
      }
    }
    if (ctr >= 0 && pAddend1->limbs[ctr].x < pAddend2->limbs[ctr].x)
    {
      pTemp = pAddend1;
      pAddend1 = pAddend2;
      pAddend2 = pTemp;
    }           // At this moment, the absolute value of addend1 is greater than
                // or equal than the absolute value of addend2.
  }
  nbrLimbs = pAddend2->nbrLimbs;
  ptrAddend1 = pAddend1->limbs;
  ptrAddend2 = pAddend2->limbs;
  ptrSum = pSum->limbs;
  if (pAddend1->sign == pAddend2->sign)
  {             // Both addends have the same sign. Sum their absolute values.
    unsigned int carry = 0;
    for (ctr = 0; ctr < nbrLimbs; ctr++)
    {
      carry = (carry >> BITS_PER_GROUP) + (unsigned int)(ptrAddend1++)->x +
        (unsigned int)(ptrAddend2++)->x;
      (ptrSum++)->x = (int)(carry & MAX_INT_NBR);
    }
    nbrLimbs = pAddend1->nbrLimbs;
    for (; ctr < nbrLimbs; ctr++)
    {
      carry = (carry >> BITS_PER_GROUP) + (unsigned int)(ptrAddend1++)->x;
      (ptrSum++)->x = (int)(carry & MAX_INT_NBR);
    }
    if (carry >= LIMB_RANGE)
    {
      ptrSum->x = 1;
      nbrLimbs++;
    }
  }
  else
  {           // Both addends have different sign. Subtract their absolute values.
    int borrow = 0;
    for (ctr = 0; ctr < nbrLimbs; ctr++)
    {
      borrow = (borrow >> BITS_PER_INT_GROUP) + (ptrAddend1++)->x - (ptrAddend2++)->x;
      (ptrSum++)->x = borrow & MAX_INT_NBR;
    }
    nbrLimbs = pAddend1->nbrLimbs;
    for (; ctr < nbrLimbs; ctr++)
    {
      borrow = (borrow >> BITS_PER_INT_GROUP) + (ptrAddend1++)->x;
      (ptrSum++)->x = borrow & MAX_INT_NBR;
    }
    while (nbrLimbs > 1 && pSum->limbs[nbrLimbs - 1].x == 0)
    {     // Loop that deletes non-significant zeros.
      nbrLimbs--;
    }
  }
  pSum->nbrLimbs = nbrLimbs;
  pSum->sign = pAddend1->sign;
}

void BigIntNegate(BigInteger *pSrc, BigInteger *pDest)
{
  if (pSrc != pDest)
  {
    CopyBigInt(pDest, pSrc);
  }
  if (pSrc->sign == SIGN_POSITIVE && (pSrc->nbrLimbs!=1 || pSrc->limbs[0].x != 0))
  {
    pDest->sign = SIGN_NEGATIVE;
  }
  else
  {
    pDest->sign = SIGN_POSITIVE;
  }
}

void BigIntSubt(BigInteger *pMinuend, BigInteger *pSubtrahend, BigInteger *pDifference)
{
  BigIntNegate(pSubtrahend, pSubtrahend);
  BigIntAdd(pMinuend, pSubtrahend, pDifference);
  if (pSubtrahend != pDifference)
  {
    BigIntNegate(pSubtrahend, pSubtrahend);    // Put original sign back.
  }
}

enum eExprErr BigIntMultiply(BigInteger *pFactor1, BigInteger *pFactor2, BigInteger *pProduct)
{
  int nbrLimbsFactor1 = pFactor1->nbrLimbs;
  int nbrLimbsFactor2 = pFactor2->nbrLimbs;
  int nbrLimbs;
  if ((pFactor1->nbrLimbs == 1 && pFactor1->limbs[0].x == 0) ||
    (pFactor2->nbrLimbs == 1 && pFactor2->limbs[0].x == 0))
  {    // Any factor is zero.
    pProduct->nbrLimbs = 1;      // Product is zero.
    pProduct->limbs[0].x = 0;
    pProduct->sign = SIGN_POSITIVE;
    return EXPR_OK;
  }
  if (pFactor1->nbrLimbs + pFactor2->nbrLimbs > 66438 / BITS_PER_GROUP)  // 2^66438 ~ 10^20000
  {
    return EXPR_INTERM_TOO_HIGH;
  }
  if (nbrLimbsFactor1<nbrLimbsFactor2)
  {
    memset(&pFactor1->limbs[nbrLimbsFactor1], 0, (nbrLimbsFactor2 - nbrLimbsFactor1)*sizeof(limb));
    nbrLimbs = nbrLimbsFactor2;
  }
  else
  {
    memset(&pFactor2->limbs[nbrLimbsFactor2], 0, (nbrLimbsFactor1 - nbrLimbsFactor2)*sizeof(limb));
    nbrLimbs = nbrLimbsFactor1;
  }
  multiply(&pFactor1->limbs[0], &pFactor2->limbs[0], &pProduct->limbs[0], nbrLimbs, &nbrLimbs);
  nbrLimbs = nbrLimbsFactor1 + nbrLimbsFactor2;
  if (pProduct->limbs[nbrLimbs - 1].x == 0)
  {
    nbrLimbs--;
  }
  pProduct->nbrLimbs = nbrLimbs;
  if (nbrLimbs == 1 && pProduct->limbs[0].x == 0)
  {
    pProduct->sign = SIGN_POSITIVE;
  }
  else
  {
    if (pFactor1->sign == pFactor2->sign)
    {
      pProduct->sign = SIGN_POSITIVE;
    }
    else
    {
      pProduct->sign = SIGN_NEGATIVE;
    }
  }
  return EXPR_OK;
}

enum eExprErr BigIntRemainder(BigInteger *pDividend, BigInteger *pDivisor, BigInteger *pRemainder)
{
  enum eExprErr rc;
  if (pDivisor->limbs[0].x == 0 && pDivisor->nbrLimbs == 1)
  {   // If divisor = 0, then remainder is the dividend.
    return EXPR_OK;
  }
  CopyBigInt(&Temp2, pDividend);
  rc = BigIntDivide(pDividend, pDivisor, &Base);   // Get quotient of division.
  if (rc != EXPR_OK)
  {
    return rc;
  }
  rc = BigIntMultiply(&Base, pDivisor, &Base);
  if (rc != EXPR_OK)
  {
    return rc;
  }
  BigIntSubt(&Temp2, &Base, pRemainder);
  return EXPR_OK;
}

void intToBigInteger(BigInteger *bigint, int value)
{
  if (value >= 0)
  {
    bigint->limbs[0].x = value;
    bigint->sign = SIGN_POSITIVE;
  }
  else
  {
    bigint->limbs[0].x = -value;
    bigint->sign = SIGN_NEGATIVE;
  }
  bigint->nbrLimbs = 1;
}

void longToBigInteger(BigInteger *bigint, long long value)
{
  int nbrLimbs = 0;
  bigint->sign = SIGN_POSITIVE;
  if (value < 0)
  {
    bigint->sign = SIGN_NEGATIVE;
    value = -value;
  }
  do
  {
    bigint->limbs[nbrLimbs++].x = (int)value & MAX_VALUE_LIMB;
    value >>= BITS_PER_GROUP;
  } while (value != 0);
  bigint->nbrLimbs = nbrLimbs;
}

void expBigNbr(BigInteger *bignbr, double logar)
{
  int mostSignificantLimb;
  logar /= log(2);
  bignbr->sign = SIGN_POSITIVE;
  bignbr->nbrLimbs = (int)floor(logar / BITS_PER_GROUP);
  mostSignificantLimb = (int)floor(exp((logar - BITS_PER_GROUP*bignbr->nbrLimbs) * log(2)) + 0.5);
  if (mostSignificantLimb == LIMB_RANGE)
  {
    mostSignificantLimb = 1;
    bignbr->nbrLimbs++;
  }
  bignbr->nbrLimbs++;
  memset(bignbr->limbs, 0, bignbr->nbrLimbs * sizeof(limb));
  bignbr->limbs[bignbr->nbrLimbs-1].x = mostSignificantLimb;
}

double logBigNbr(BigInteger *pBigNbr)
{
  int nbrLimbs;
  double logar;
  nbrLimbs = pBigNbr->nbrLimbs;
  if (nbrLimbs == 1)
  {
    logar = log((double)(pBigNbr->limbs[0].x));
  }
  else
  {
    double value = pBigNbr->limbs[nbrLimbs - 2].x +
                  (double)pBigNbr->limbs[nbrLimbs - 1].x * LIMB_RANGE;
    if (nbrLimbs == 2)
    {
      logar = log(value);
    }
    else
    {
      logar = log(value + (double)pBigNbr->limbs[nbrLimbs - 3].x / LIMB_RANGE);
    }
    logar += (double)((nbrLimbs - 2)*BITS_PER_GROUP)*log(2);
  }
  return logar;
}

double logLimbs(limb *pBigNbr, int nbrLimbs)
{
  double logar;
  if (nbrLimbs > 1)
  {
    logar = log((double)((pBigNbr + nbrLimbs - 2)->x +
                ((double)(pBigNbr + nbrLimbs - 1)->x * LIMB_RANGE))) +
      (double)(nbrLimbs - 2)*log((double)LIMB_RANGE);
  }
  else
  {
    logar = log((double)((pBigNbr+nbrLimbs - 1)->x)) +
      (double)(nbrLimbs - 1)*log((double)LIMB_RANGE);
  }
  return logar;
}

enum eExprErr BigIntPowerIntExp(BigInteger *pBase, int expon, BigInteger *pPower)
{
  int mask;
  double base;
  if (pBase->nbrLimbs == 1 && pBase->limbs[0].x == 0)
  {     // Base = 0 -> power = 0
    pPower->limbs[0].x = 0;
    pPower->nbrLimbs = 1;
    pPower->sign = SIGN_POSITIVE;
    return EXPR_OK;
  }
  base = logBigNbr(pBase);
  if (base*(double)expon > 46051)
  {   // More than 20000 digits. 46051 = log(10^20000)
    return EXPR_INTERM_TOO_HIGH;
  }
  CopyBigInt(&Base, pBase);
  pPower->sign = SIGN_POSITIVE;
  pPower->nbrLimbs = 1;
  pPower->limbs[0].x = 1;
  for (mask = 1 << 30; mask != 0; mask >>= 1)
  {
    if ((expon & mask) != 0)
    {
      for (; mask != 0; mask >>= 1)
      {
        (void)BigIntMultiply(pPower, pPower, pPower);
        if ((expon & mask) != 0)
        {
          (void)BigIntMultiply(pPower, &Base, pPower);
        }
      }
      break;
    }
  }
  return EXPR_OK;
}

enum eExprErr BigIntPower(BigInteger *pBase, BigInteger *pExponent, BigInteger *pPower)
{
  int expon;
  if (pExponent->sign == SIGN_NEGATIVE)
  {    // Negative exponent not accepted.
    return EXPR_INVALID_PARAM;
  }
  if (pExponent->nbrLimbs > 1)
  {     // Exponent too high.
    if (pBase->nbrLimbs == 1 && pBase->limbs[0].x < 2)
    {     // Base = 0 -> power = 0
      pPower->limbs[0].x = pBase->limbs[0].x;
      pPower->nbrLimbs = 1;
      pPower->sign = SIGN_POSITIVE;
      return EXPR_OK;
    }
    return EXPR_INTERM_TOO_HIGH;
  }
  if (pExponent->nbrLimbs == 2)
  {
    expon = pExponent->limbs[0].x + (pExponent->limbs[1].x << BITS_PER_GROUP);
  }
  else
  {
    expon = pExponent->limbs[0].x;
  }
  return BigIntPowerIntExp(pBase, expon, pPower);
}



// GCD of two numbers:
// Input: a, b positive integers
// Output : g and d such that g is odd and gcd(a, b) = g×2d
//   d : = 0
//   while a and b are both even do
//     a : = a / 2
//     b : = b / 2
//     d : = d + 1
//   while a != b do
//     if a is even then a : = a / 2
//     else if b is even then b : = b / 2
//     else if a > b then a : = (a – b) / 2
//     else b : = (b – a) / 2
//   g : = a
//     output g, d

void BigIntDivide2(BigInteger *pArg)
{
  int nbrLimbs = pArg->nbrLimbs;
  int ctr = nbrLimbs - 1;
  unsigned int carry;
  limb *ptrLimb = &pArg->limbs[ctr];
  carry = 0;
  for (; ctr >= 0; ctr--)
  {
    carry = (carry << BITS_PER_GROUP) + (unsigned int)ptrLimb->x;
    (ptrLimb--)->x = (int)(carry >> 1);
    carry &= 1;
  }
  if (nbrLimbs > 1 && pArg->limbs[nbrLimbs - 1].x == 0)
  {     // Most significant limb is zero, so reduce size by one limb.
    pArg->nbrLimbs--;
  }
}

static void BigIntMutiplyPower2(BigInteger *pArg, int power2)
{
  int ctr;
  int nbrLimbs = pArg->nbrLimbs;
  limb *ptrLimbs = pArg->limbs;
  for (; power2 > 0; power2--)
  {
    unsigned int carry = 0;
    for (ctr = 0; ctr < nbrLimbs; ctr++)
    {
      carry += (unsigned int)(ptrLimbs + ctr)->x << 1;
      (ptrLimbs + ctr)->x = (int)(carry & MAX_VALUE_LIMB);
      carry >>= BITS_PER_GROUP;
    }
    if (carry != 0)
    {
      (ptrLimbs + ctr)->x = (int)carry;
      nbrLimbs++;
    }
  }
  pArg->nbrLimbs = nbrLimbs;
}

boolean TestBigNbrEqual(BigInteger *pNbr1, BigInteger *pNbr2)
{
  int ctr;
  limb *ptrLimbs1 = pNbr1->limbs;
  limb *ptrLimbs2 = pNbr2->limbs;
  if (pNbr1->nbrLimbs != pNbr2->nbrLimbs)
  {        // Sizes of numbers are different.
    return FALSE;
  }
  if (pNbr1->sign != pNbr2->sign)
  {        // Sign of numbers are different.
    if (pNbr1->nbrLimbs == 1 && pNbr1->limbs[0].x == 0 && pNbr2->limbs[0].x == 0)
    {              // Both numbers are zero.
      return TRUE;
    }
    return FALSE;
  }

           // Check whether both numbers are equal.
  for (ctr = pNbr1->nbrLimbs - 1; ctr >= 0; ctr--)
  {
    if ((ptrLimbs1 + ctr)->x != (ptrLimbs2 + ctr)->x)
    {      // Numbers are different.
      return FALSE;
    }
  }        // Numbers are equal.
  return TRUE;
}

void BigIntGcd(BigInteger *pArg1, BigInteger *pArg2, BigInteger *pResult)
{
  int nbrLimbs1 = pArg1->nbrLimbs;
  int nbrLimbs2 = pArg2->nbrLimbs;
  int power2;
  if ((nbrLimbs1 == 1 && pArg1->limbs[0].x == 0) || (nbrLimbs2 == 1 && pArg2->limbs[0].x == 0))
  {               // Some argument is zero, so the GCD is zero.
    pResult->limbs[0].x = 0;
    pResult->nbrLimbs = 1;
    pResult->sign = SIGN_POSITIVE;
    return;
  }
  // Reuse Base and Power temporary variables.
  CopyBigInt(&Base, pArg1);
  CopyBigInt(&Power, pArg2);
  Base.sign = SIGN_POSITIVE;
  Power.sign = SIGN_POSITIVE;
  power2 = 0;
  while (((Base.limbs[0].x | Power.limbs[0].x) & 1) == 0)
  {  // Both values are even
    BigIntDivide2(&Base);
    BigIntDivide2(&Power);
    power2++;
  }
  while (TestBigNbrEqual(&Base, &Power) == 0)
  {    // Main GCD loop.
    if ((Base.limbs[0].x & 1) == 0)
    {          // Number is even. Divide it by 2.
      BigIntDivide2(&Base);
      continue;
    }
    if ((Power.limbs[0].x & 1) == 0)
    {          // Number is even. Divide it by 2.
      BigIntDivide2(&Power);
      continue;
    }
    BigIntSubt(&Base, &Power, pResult);
    if (pResult->sign == SIGN_POSITIVE)
    {
      CopyBigInt(&Base, pResult);
      BigIntDivide2(&Base);
    }
    else
    {
      CopyBigInt(&Power, pResult);
      Power.sign = SIGN_POSITIVE;
      BigIntDivide2(&Power);
    }
  }
  CopyBigInt(pResult, &Base);
  BigIntMutiplyPower2(pResult, power2);
}

static void addToAbsValue(limb *pLimbs, int *pNbrLimbs, int addend)
{
  int ctr;
  int nbrLimbs = *pNbrLimbs;
  pLimbs->x += addend;
  if ((unsigned int)pLimbs->x < LIMB_RANGE)
  {     // No overflow. Go out of routine.
    return;
  }
    pLimbs->x -= LIMB_RANGE;
    for (ctr = 1; ctr < nbrLimbs; ctr++)
    {
    pLimbs++;        // Point to next most significant limb.
    if (pLimbs->x != MAX_INT_NBR)
    {   // No overflow. Go out of routine.
      (pLimbs->x)++;   // Add carry.
      return;
    }
    pLimbs->x = 0;
  }
  (*pNbrLimbs)++;        // Result has an extra limb.
  (pLimbs + 1)->x = 1;   // Most significant limb must be 1.
}

static void subtFromAbsValue(limb *pLimbs, int *pNbrLimbs, int subt)
{
  int nbrLimbs = *pNbrLimbs;
  pLimbs->x -= subt;
  if (pLimbs->x < 0)
  {
    int ctr;
    for (ctr = 1; ctr < nbrLimbs; ctr++)
    {
      (pLimbs + ctr - 1)->x += LIMB_RANGE;
      if (--((pLimbs + ctr)->x) >= 0)
      {
        break;
      }
    }
    if (nbrLimbs > 1 && (pLimbs + nbrLimbs - 1)->x == 0)
    {
      nbrLimbs--;
    }
  }
  *pNbrLimbs = nbrLimbs;
}

void subtractdivide(BigInteger *pBigInt, int subt, int divisor)
{
  int nbrLimbs = pBigInt->nbrLimbs;
  // Point to most significant limb.
  limb *pLimbs;
  int ctr;
  int remainder = 0;
  double dDivisor = (double)divisor;
  double dLimb = (double)LIMB_RANGE;
  if (subt >= 0)
  {
    if (pBigInt->sign == SIGN_POSITIVE)
    {               // Subtract subt to absolute value.
      subtFromAbsValue(pBigInt->limbs, &nbrLimbs, subt);
    }
    else
    {               // Add subt to absolute value.
      addToAbsValue(pBigInt->limbs, &nbrLimbs, subt);
    }
  }
  else
  {
    if (pBigInt->sign == SIGN_POSITIVE)
    {               // Subtract subt to absolute value.
      addToAbsValue(pBigInt->limbs, &nbrLimbs, -subt);
    }
    else
    {               // Add subt to absolute value.
      subtFromAbsValue(pBigInt->limbs, &nbrLimbs, -subt);
    }
  }
  pLimbs = pBigInt->limbs + nbrLimbs - 1;
  // Divide number by divisor.
  for (ctr = nbrLimbs - 1; ctr >= 0; ctr--)
  {
    double dDividend, dQuotient;
    unsigned int quotient, dividend;
    dividend = (remainder << BITS_PER_INT_GROUP) + pLimbs->x;
    dDividend = (double)remainder * dLimb + pLimbs->x;
    dQuotient = floor(dDividend / dDivisor + 0.5);
    quotient = (unsigned int)dQuotient;   // quotient has correct value or 1 more.
    remainder = dividend - quotient * divisor;
    if ((unsigned int)remainder >= (unsigned int)divisor)
    {     // remainder not in range 0 <= remainder < divisor. Adjust.
      quotient--;
      remainder += divisor;
    }
    (pLimbs--)->x = (int)quotient;
  }
  if (nbrLimbs > 1 && pBigInt->limbs[nbrLimbs - 1].x == 0)
  {   // Most significant limb is now zero, so discard it.
    nbrLimbs--;
  }
  pBigInt->nbrLimbs = nbrLimbs;
}

int getRemainder(BigInteger *pBigInt, int divisor)
{
  int ctr;
  int remainder = 0;
  int nbrLimbs = pBigInt->nbrLimbs;
  double dDivisor = (double)divisor;
  double dLimb = 0x80000000;
  limb *pLimb = &pBigInt->limbs[nbrLimbs - 1];
  for (ctr = nbrLimbs - 1; ctr >= 0; ctr--)
  {
    int quotient, dividend;
    double dQuotient, dDividend;
    dividend = (remainder << BITS_PER_INT_GROUP) + pLimb->x;
    dDividend = (double)remainder * dLimb + pLimb->x;
    dQuotient = floor(dDividend / dDivisor + 0.5);
    quotient = (int)(unsigned int)dQuotient;   // quotient has correct value or 1 more.
    remainder = dividend - quotient * divisor;
    if ((unsigned int)remainder >= (unsigned int)divisor)
    {     // remainder not in range 0 <= remainder < divisor. Adjust.
      quotient--;
      remainder += divisor;
    }
    pLimb--;
  }
  if (pBigInt->sign == SIGN_NEGATIVE && remainder != 0)
  {
    remainder = divisor - remainder;
  }
  return remainder;
}

void addbigint(BigInteger *pResult, int addend)
{
  int sign;
  int nbrLimbs = pResult->nbrLimbs;
  limb *pResultLimbs = pResult->limbs;
  sign = pResult->sign;
  if (addend < 0)
  {
    addend = -addend;
    if (sign == SIGN_POSITIVE)
    {
      sign = SIGN_NEGATIVE;
    }
    else
    {
      sign = SIGN_POSITIVE;
    }
  }
  if (sign == SIGN_POSITIVE)
  {   // Add addend to absolute value of pResult.
    addToAbsValue(pResultLimbs, &nbrLimbs, addend);
  }
  else
  {  // Subtract addend from absolute value of pResult.
    if (nbrLimbs == 1)
    {
      pResultLimbs->x -= addend;
      if (pResultLimbs->x < 0)
      {
        pResultLimbs->x = -pResultLimbs->x;
        BigIntNegate(pResult, pResult);
      }
    }
    else
    {     // More than one limb.
      subtFromAbsValue(pResultLimbs, &nbrLimbs, addend);
    }
  }
  pResult->nbrLimbs = nbrLimbs;
}

void multint(BigInteger *pResult, BigInteger *pMult, int factor)
{
  double dFactor;
  double dVal = 1 / (double)LIMB_RANGE;
  int factorPositive = 1;
  int ctr, carry;
  int nbrLimbs = pMult->nbrLimbs;
  limb *pLimb = pMult->limbs;
  limb *pResultLimb = pResult->limbs;
  if (factor < 0)
  {     // If factor is negative, indicate it and compute its absolute value.
    factorPositive = 0;
    factor = -factor;
  }
  dFactor = (double)factor;
  carry = 0;
  for (ctr = 0; ctr < nbrLimbs; ctr++)
  {
    int low = (pLimb->x * factor + carry) & MAX_INT_NBR;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    if (low < HALF_INT_RANGE)
    {
      carry = (int)(((double)(pLimb->x) * dFactor + (double)carry + HALF_INT_RANGE / 2)*dVal);
    }
    else
    {
      carry = (int)(((double)(pLimb->x) * dFactor + (double)carry - HALF_INT_RANGE / 2)*dVal);
    }
    (pResultLimb++)->x = low;
    pLimb++;
  }
  if (carry != 0)
  {
    pResultLimb->x = carry;
    nbrLimbs++;
  }
  pResult->nbrLimbs = nbrLimbs;
  pResult->sign = pMult->sign;
  if (factorPositive == 0)
  {
    BigIntNegate(pResult, pResult);
  }
}

void multadd(BigInteger *pResult, int iMult, BigInteger *pMult, int addend)
{
  multint(pResult, pMult, iMult);
  addbigint(pResult, addend);
}

// Compute *pResult -> *pMult1 * iMult1 + *pMult2 * iMult2
// Use Temp as temporary variable.
void addmult(BigInteger *pResult, BigInteger *pMult1, int iMult1, BigInteger *pMult2, int iMult2)
{
  multint(pResult, pMult1, iMult1);
  multint(&Temp, pMult2, iMult2);
  BigIntAdd(pResult, &Temp, pResult);
}

// Get number of bits of given big integer.
int bitLength(BigInteger *pBigNbr)
{
  unsigned int mask;
  int bitCount;
  int lastLimb = pBigNbr->nbrLimbs-1;
  int bitLen = lastLimb*BITS_PER_GROUP;
  unsigned int limb = (unsigned int)(pBigNbr->limbs[lastLimb].x);
  mask = 1;
  for (bitCount = 0; bitCount < BITS_PER_GROUP; bitCount++)
  {
    if (limb < mask)
    {
      break;
    }
    mask *= 2;
  }
  return bitLen + bitCount;
}

int intModPow(int NbrMod, int Expon, int currentPrime)
{
  unsigned int power = 1;
  unsigned int square = (unsigned int)NbrMod;
  while (Expon != 0)
  {
    if ((Expon & 1) == 1)
    {
      power = (power * square) % (unsigned int)currentPrime;
    }
    square = (square * square) % (unsigned int)currentPrime;
    Expon >>= 1;
  }
  return (int)power;
}

static void InitTempFromInt(int value)
{
  Temp.nbrLimbs = 1;
  Temp.limbs[0].x = value;
  Temp.sign = SIGN_POSITIVE;
}

void UncompressBigInteger(/*@in@*/int *ptrValues, /*@out@*/BigInteger *bigint)
{
  limb *destLimb = bigint->limbs;
  if (NumberLength == 1)
  {
    destLimb->x = *(ptrValues + 1);
    bigint->nbrLimbs = 1;
  }
  else
  {
    int ctr, nbrLimbs;
    nbrLimbs = *ptrValues++;
    bigint->nbrLimbs = nbrLimbs;
    for (ctr = 0; ctr < nbrLimbs; ctr++)
    {
      (destLimb++)->x = *ptrValues++;
    }
    for (; ctr < NumberLength; ctr++)
    {
      (destLimb++)->x = 0;
    }
  }
}

void CompressBigInteger(/*@out@*/int *ptrValues, /*@in@*/BigInteger *bigint)
{
  limb *destLimb = bigint->limbs;
  if (NumberLength == 1)
  {
    *ptrValues = 1;
    *(ptrValues + 1) = (int)(destLimb->x);
  }
  else
  {
    int ctr, nbrLimbs;
    nbrLimbs = getNbrLimbs(bigint->limbs);
    *ptrValues++ = nbrLimbs;
    for (ctr = 0; ctr < nbrLimbs; ctr++)
    {
      *ptrValues++ = (int)((destLimb++)->x);
    }
  }
}

void UncompressLimbsBigInteger(/*@in@*/limb *ptrValues, /*@out@*/BigInteger *bigint)
{
  if (NumberLength == 1)
  {
    bigint->limbs[0].x = ptrValues->x;
    bigint->nbrLimbs = 1;
  }
  else
  {
    int nbrLimbs;
    limb *ptrValue1;
    memcpy(bigint->limbs, ptrValues, NumberLength*sizeof(limb));
    ptrValue1 = ptrValues + NumberLength;
    for (nbrLimbs = NumberLength; nbrLimbs > 1; nbrLimbs--)
    {
      --ptrValue1;
      if (ptrValue1->x != 0)
      {
        break;
      }
    }
    bigint->nbrLimbs = nbrLimbs;
  }
}

void CompressLimbsBigInteger(/*@out@*/limb *ptrValues, /*@in@*/BigInteger *bigint)
{
  if (NumberLength == 1)
  {
    ptrValues->x = bigint->limbs[0].x;
  }
  else
  {
    int nbrLimbs = bigint->nbrLimbs;
    if (nbrLimbs > NumberLength)
    {
      memcpy(ptrValues, bigint->limbs, NumberLength*sizeof(limb));
    }
    else
    {
      memcpy(ptrValues, bigint->limbs, nbrLimbs * sizeof(limb));
      memset(ptrValues + nbrLimbs, 0, (NumberLength - nbrLimbs) * sizeof(limb));
    }
  }
}

void UncompressIntLimbs(/*@in@*/int *ptrValues, /*@out@*/limb *bigint, int nbrLen)
{
  int nbrLimbs = *ptrValues;
  memcpy(bigint, ptrValues+1, nbrLimbs*sizeof(limb));
  memset(bigint + nbrLimbs + 1, 0, (nbrLen - nbrLimbs) * sizeof(limb));
}

void CompressIntLimbs(/*@out@*/int *ptrValues, /*@in@*/limb *bigint, int nbrLen)
{
  int nbrLimbs;
  memcpy(ptrValues+1, bigint, (nbrLen-1) * sizeof(limb));
  for (nbrLimbs = nbrLen-1; nbrLimbs > 1; nbrLimbs--)
  {
    if (*(ptrValues + nbrLimbs) != 0)
    {
      break;
    }
  }
  *ptrValues = nbrLimbs;
}

// This routine checks whether the number pointed by pNbr is
// a perfect power. If it is not, it returns one.
// If it is a perfect power, it returns the exponent and 
// it fills the buffer pointed by pBase with the base.
int PowerCheck(BigInteger *pBigNbr, BigInteger *pBase)
{
  limb *ptrLimb;
  double dN;
  int nbrLimbs = pBigNbr->nbrLimbs;
  int maxExpon = bitLength(pBigNbr);
  int h, j;
  int modulus;
  int intLog2root;
  int primesLength, Exponent;
  double log2N, log2root;
  int prime2310x1[] =
  { 2311, 4621, 9241, 11551, 18481, 25411, 32341, 34651, 43891, 50821 };
  // Primes of the form 2310x+1.
  boolean expon2 = TRUE, expon3 = TRUE, expon5 = TRUE;
  boolean expon7 = TRUE, expon11 = TRUE;
  for (h = 0; h < sizeof(prime2310x1) / sizeof(prime2310x1[0]); h++)
  {
    int testprime = prime2310x1[h];
    int mod = getRemainder(pBigNbr, testprime);
    if (expon2 && intModPow(mod, testprime / 2, testprime) > 1)
    {
      expon2 = FALSE;
    }
    if (expon3 && intModPow(mod, testprime / 3, testprime) > 1)
    {
      expon3 = FALSE;
    }
    if (expon5 && intModPow(mod, testprime / 5, testprime) > 1)
    {
      expon5 = FALSE;
    }
    if (expon7 && intModPow(mod, testprime / 7, testprime) > 1)
    {
      expon7 = FALSE;
    }
    if (expon11 && intModPow(mod, testprime / 11, testprime) > 1)
    {
      expon11 = FALSE;
    }
  }
  primesLength = 2 * maxExpon + 3;
  for (h = 2; h <= maxExpon; h++)
  {
    ProcessExpon[h] = TRUE;
  }
  for (h = 2; h < primesLength; h++)
  {
    primes[h] = TRUE;
  }
  for (h = 2; h * h < primesLength; h++)
  { // Generation of primes
    for (j = h * h; j < primesLength; j += h)
    { // using Eratosthenes sieve
      primes[j] = FALSE;
    }
  }
  for (h = 13; h < primesLength; h++)
  {
    if (primes[h])
    {
      int processed = 0;
      for (j = 2 * h + 1; j < primesLength; j += 2 * h)
      {
        if (primes[j])
        {
          modulus = getRemainder(pBigNbr, j);
          if (intModPow(modulus, j / h, j) > 1)
          {
            for (j = h; j <= maxExpon; j += h)
            {
              ProcessExpon[j] = FALSE;
            }
            break;
          }
        }
        if (++processed > 10)
        {
          break;
        }
      }
    }
  }
  log2N = logBigNbr(pBigNbr) / log(2);
  for (Exponent = maxExpon; Exponent >= 2; Exponent--)
  {
    if (Exponent % 2 == 0 && !expon2)
    {
      continue; // Not a square
    }
    if (Exponent % 3 == 0 && !expon3)
    {
      continue; // Not a cube
    }
    if (Exponent % 5 == 0 && !expon5)
    {
      continue; // Not a fifth power
    }
    if (Exponent % 7 == 0 && !expon7)
    {
      continue; // Not a 7th power
    }
    if (Exponent % 11 == 0 && !expon11)
    {
      continue; // Not an 11th power
    }
    if (!ProcessExpon[Exponent])
    {
      continue;
    }
    // Initialize approximation to n-th root (n = Exponent).
    log2root = log2N / Exponent;
    intLog2root = (int)floor(log2root/ BITS_PER_GROUP);
    nbrLimbs = intLog2root + 1;
    ptrLimb = &pBase->limbs[nbrLimbs - 1];
    dN = exp((log2root - intLog2root*BITS_PER_GROUP) * log(2));
    // All approximations must be >= than true answer.
    if (nbrLimbs == 1)
    {
      ptrLimb->x = (int)(unsigned int)ceil(dN);
      if ((unsigned int)ptrLimb->x == LIMB_RANGE)
      {
        nbrLimbs = 2;
        ptrLimb->x = 0;
        (ptrLimb+1)->x = 1;
      }
    }
    else
    {
      dN += 1 / (double)LIMB_RANGE;
      ptrLimb->x = (int)trunc(dN);
      dN -= trunc(dN);
      (ptrLimb - 1)->x = (int)trunc(dN*LIMB_RANGE);
    }
    pBase->nbrLimbs = nbrLimbs;
    // Perform Newton iteration for n-th root.
    for (;;)
    {   // Check whether the approximate root is actually exact.
      BigIntPowerIntExp(pBase, Exponent-1, &Temp3); // Temp3 <- x^(e-1)
      BigIntMultiply(&Temp3, pBase, &Temp2);        // Temp2 <- x^e 
      BigIntSubt(pBigNbr, &Temp2, &Temp2);            // Compare to radicand.
      if (Temp2.nbrLimbs == 1 && Temp2.limbs[0].x == 0)
      {                     // Perfect power, so go out.
        return Exponent;
      }
      if (Temp2.sign == SIGN_POSITIVE)
      {                     // x^e > radicand -> not perfect power, so go out.
        break;
      }
      BigIntDivide(pBigNbr, &Temp3, &Temp);         // Temp -> N/x^(e-1)
      BigIntSubt(&Temp, pBase, &Temp2);             // Temp2 -> N/x^(e-1) - x
      if (Temp2.nbrLimbs == 1 && Temp2.limbs[0].x == 0)
      {     // New approximation will be the same as previous. Go out.
        break;
      }
      InitTempFromInt(Exponent - 1);
      BigIntSubt(&Temp2, &Temp, &Temp2);
      InitTempFromInt(Exponent);
      BigIntDivide(&Temp2, &Temp, &Temp2);
      BigIntAdd(&Temp2, pBase, pBase);
    }
  }
  pBase->nbrLimbs = pBigNbr->nbrLimbs;
  memcpy(pBase->limbs, pBigNbr->limbs, pBase->nbrLimbs*sizeof(limb));
  return 1;
}

static int checkOne(limb *value, int nbrLimbs)
{
  int idx;
  for (idx = 0; idx < nbrLimbs; idx++)
  {
    if ((value++)->x != MontgomeryMultR1[idx].x)
    {
      return 0;    // Go out if value is not 1 (mod p)
    }
  }
  return 1;
}

int checkMinusOne(limb *value, int nbrLimbs)
{
  int idx;
  unsigned int carry;
  carry = 0;
  for (idx = 0; idx < nbrLimbs; idx++)
  {
    carry += (unsigned int)(value++)->x + (unsigned int)MontgomeryMultR1[idx].x;
    if ((carry & MAX_VALUE_LIMB) != (unsigned int)TestNbr[idx].x)
    {
      return 0;    // Go out if value is not -1 (mod p)
    }
    carry >>= BITS_PER_GROUP;
  }
  return 1;
}

// Find power of 2 that divides the number.
// output: pNbrLimbs = pointer to number of limbs
//         pShRight = pointer to power of 2.
void DivideBigNbrByMaxPowerOf2(int *pShRight, limb *number, int *pNbrLimbs)
{
  int power2 = 0;
  int index, index2, mask, shRg;
  int nbrLimbs = *pNbrLimbs;
  // Start from least significant limb (number zero).
  for (index = 0; index < nbrLimbs; index++)
  {
    if (number[index].x != 0)
    {
      break;
    }
    power2 += BITS_PER_GROUP;
  }
  for (mask = 0x1; mask <= MAX_VALUE_LIMB; mask *= 2)
  {
    if ((number[index].x & mask) != 0)
    {
      break;
    }
    power2++;
  }
  // Divide number by this power.
  shRg = power2 % BITS_PER_GROUP; // Shift right bit counter
  if ((number[nbrLimbs - 1].x & (-(1 << shRg))) != 0)
  {   // Most significant bits set.
    *pNbrLimbs = nbrLimbs - index;
  }
  else
  {   // Most significant bits not set.
    *pNbrLimbs = nbrLimbs - index - 1;
  }
      // Move number shRg bits to the right.
  mask = (1 << shRg) - 1;
  for (index2 = index; index2 < nbrLimbs; index2++)
  {
    number[index2].x = ((number[index2].x >> shRg) |
                        (number[index2+1].x << (BITS_PER_GROUP - shRg))) &
                        MAX_VALUE_LIMB;
  }
  if (index > 0)
  {   // Move limbs to final position.
    memmove(number, &number[index], (nbrLimbs - index) * sizeof(limb));
  }
  *pShRight = power2;
}

// Check if number is pseudoprime base 3
int isPseudoprime(BigInteger *pResult)
{
  int base, delta, baseNbr, ctr, i, nbrLimbsQ, Mult3Len;
  limb largeVal;
  int nbrLimbs = pResult->nbrLimbs;
  limb *pResultLimbs = pResult->limbs;
  if (nbrLimbs == 1)
  {
    largeVal.x = pResultLimbs->x;
    if (largeVal.x >= 9)
    {
      int Q;
      for (Q = 3; Q*Q <= largeVal.x; Q += 2)
      {     // Check if base is prime
        if (largeVal.x % Q == 0)
        {
          break;     // Composite
        }
      }
      if (Q*Q <= largeVal.x)
      {
        return FALSE;       // Composite. Try next candidate. 
      }
    }
    return TRUE;
  }
  base = 3;
  delta = 2;
  for (baseNbr = 100; baseNbr > 0; baseNbr--)
  {    // Compute value mod base. If it is zero, the number is composite.
    if (getRemainder(pResult, base) == 0)
    {                      // Number is composite.
      return FALSE;
    }
    base += delta;         // Skip multiples of 3.
    if (base > 5)
    {
      delta = 6 - delta;   // Exchange delta between 2 and 4.
    }
  }
  base = 3;
  (pResultLimbs + nbrLimbs)->x = 0;
  memcpy(q, pResultLimbs, (nbrLimbs + 1)*sizeof(limb));
  nbrLimbsQ = nbrLimbs;
  q[0]--;                     // q = p - 1 (p is odd, so there is no carry).
  memcpy(Mult3, q, (nbrLimbsQ + 1)*sizeof(q[0]));
  Mult3Len = nbrLimbs;
  DivideBigNbrByMaxPowerOf2(&ctr, Mult3, &Mult3Len);
  memcpy(TestNbr, pResultLimbs, nbrLimbs*sizeof(limb));
  TestNbr[nbrLimbs].x = 0;
  GetMontgomeryParms(nbrLimbs);
  for (baseNbr = 20; baseNbr > 0; baseNbr--)
  {    // Try up to 20 bases.
    modPowBaseInt(base, Mult3, Mult3Len, Mult1); // Mult1 = base^Mult3.
                                                 // If Mult1 = 1 or Mult1 = TestNbr-1, then try next base.
    if (checkOne(Mult1, nbrLimbs) != 0 || checkMinusOne(Mult1, nbrLimbs) != 0)
    {
      base += 2;
      continue;
    }
    for (i = 0; i < ctr; i++)
    {              // Loop that squares number.
      modmult(Mult1, Mult1, Mult4);
      if (checkOne(Mult4, nbrLimbs) != 0)
      {
        return FALSE;      // Number is composite
      }
      if (checkMinusOne(Mult4, nbrLimbs) != 0)
      {
        i = -1;     // Number is strong pseudoprime.
        break;
      }
      memcpy(Mult1, Mult4, nbrLimbs*sizeof(limb));
    }
    if (i != -1)
    {
      return FALSE;        // Composite. Go out SPRP loop.
    }
    // If power (Mult4) is 1, that means that number is at least PRP,
    // so continue loop trying to find square root of -1.
    base += 2;
  }
  return TRUE;
}

void NbrToLimbs(int nbr, /*@out@*/limb *limbs, int len)
{
  if (nbr >= MAX_VALUE_LIMB)
  {
    limbs->x = nbr % MAX_VALUE_LIMB;
    (limbs+1)->x = nbr / MAX_VALUE_LIMB;
    memset(limbs + 2, 0, (len-2) * sizeof(limb));
  }
  else
  {
    limbs->x = nbr;
    memset(limbs + 1, 0, (len-1) * sizeof(limb));
  }
}

int BigNbrIsZero(limb *value)
{
  int ctr;
  for (ctr = 0; ctr < NumberLength; ctr++)
  {
    if (value->x != 0)
    {
      return 0;  // Number is not zero.
    }
    value++;
  }
  return 1;      // Number is zero
}

double getMantissa(limb *ptrLimb, int nbrLimbs)
{
  double dN = (double)(ptrLimb - 1)->x;
  double dInvLimb = 1 / (double)LIMB_RANGE;
  if (nbrLimbs > 1)
  {
    dN += (double)(ptrLimb - 2)->x * dInvLimb;
  }
  if (nbrLimbs > 2)
  {
    dN += (double)(ptrLimb - 3)->x * dInvLimb * dInvLimb;
  }
  return dN;
}

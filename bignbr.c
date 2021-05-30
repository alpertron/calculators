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
#include "factor.h"
#include "expression.h"
#include "skiptest.h"

enum eOper
{
  OPERATION_AND = 0,
  OPERATION_OR,
  OPERATION_XOR,
};

static BigInteger Temp;
static BigInteger Temp2;
static BigInteger Temp3;
static BigInteger Base;
static BigInteger Power;
static BigInteger expon;
static bool ProcessExpon[(MAX_LEN*BITS_PER_GROUP) + 1000];
static bool primes[(MAX_LEN*BITS_PER_GROUP) + 1000];
extern limb Mult1[MAX_LEN];
extern limb Mult2[MAX_LEN];
extern limb Mult3[MAX_LEN];
extern limb Mult4[MAX_LEN];
extern int q[MAX_LEN];
extern limb TestNbr[MAX_LEN];
extern limb MontgomeryMultR1[MAX_LEN];
int groupLen = 6;
int smallPrimes[SMALL_PRIMES_ARRLEN+1];
#ifdef __EMSCRIPTEN__
int percentageBPSW;
#endif

void CopyBigInt(BigInteger *pDest, const BigInteger *pSrc)
{
  if (pDest != pSrc)
  {
    pDest->sign = pSrc->sign;
    pDest->nbrLimbs = pSrc->nbrLimbs;
    (void)memcpy(pDest->limbs, pSrc->limbs, (pSrc->nbrLimbs) * sizeof(limb));
  }
}

void AddBigInt(const limb *pAddend1, const limb *pAddend2, limb *pSum, int nbrLimbs)
{
  const limb *ptrAddend1 = pAddend1;
  const limb *ptrAddend2 = pAddend2;
  limb *ptrSum = pSum;
  unsigned int carry = 0;
  for (int i = 0; i < nbrLimbs; i++)
  {
    carry = (carry >> BITS_PER_GROUP) + (unsigned int)ptrAddend1->x +
                                        (unsigned int)ptrAddend2->x;
    ptrAddend1++;
    ptrAddend2++;
    ptrSum->x = (int)carry & MAX_INT_NBR;
    ptrSum++;
  }
}

void SubtractBigInt(const limb *pMinuend, const limb *pSubtrahend, limb *pDiff, int nbrLimbs)
{
  const limb *ptrMinuend = pMinuend;
  const limb *ptrSubtrahend = pSubtrahend;
  limb *ptrDiff = pDiff;
  int borrow = 0;
  for (int i = 0; i < nbrLimbs; i++)
  {
    borrow = (borrow >> BITS_PER_INT_GROUP) + ptrMinuend->x - ptrSubtrahend->x;
    ptrMinuend++;
    ptrSubtrahend++;
    ptrDiff->x = borrow & MAX_INT_NBR;
    ptrDiff++;
  }
}

// If address of num and result match, BigIntDivide will overwrite num, so it must be executed after processing num.
void floordiv(const BigInteger *num, const BigInteger *den, BigInteger *result)
{
  static BigInteger rem;
  (void)BigIntRemainder(num, den, &rem);
  if ((((num->sign == SIGN_NEGATIVE) && (den->sign == SIGN_POSITIVE)) ||
    ((num->sign == SIGN_POSITIVE) && !BigIntIsZero(num) && (den->sign == SIGN_NEGATIVE))) && !BigIntIsZero(&rem))
  {
    (void)BigIntDivide(num, den, result);
    addbigint(result, -1);
  }
  else
  {
    (void)BigIntDivide(num, den, result);
  }
}

void ceildiv(const BigInteger*num, const BigInteger *den, BigInteger *result)
{
  static BigInteger rem;
  (void)BigIntDivide(num, den, result);
  (void)BigIntRemainder(num, den, &rem);
  if ((((num->sign == SIGN_POSITIVE) && !BigIntIsZero(num) && (den->sign == SIGN_POSITIVE)) ||
    ((num->sign == SIGN_NEGATIVE) && (den->sign == SIGN_NEGATIVE))) && !BigIntIsZero(&rem))
  {
    addbigint(result, 1);
  }
}

void BigIntChSign(BigInteger *value)
{
  if ((value->nbrLimbs == 1) && (value->limbs[0].x == 0))
  {    // Value is zero. Do not change sign.
    return;
  }
  if (value->sign == SIGN_POSITIVE)
  {
    value->sign = SIGN_NEGATIVE;
  }
  else
  {
    value->sign = SIGN_POSITIVE;
  }
}

static void InternalBigIntAdd(const BigInteger *pAdd1, const BigInteger *pAdd2, 
  BigInteger *pSum, enum eSign addend2sign)
{
  const BigInteger* pAddend1 = pAdd1;
  const BigInteger* pAddend2 = pAdd2;
  int ctr;
  int nbrLimbs;
  const limb *ptrAddend1;
  const limb *ptrAddend2;
  limb *ptrSum;
  const BigInteger *pTemp;
  enum eSign addend1Sign = pAddend1->sign;
  enum eSign addend2Sign = addend2sign;
  enum eSign tmpSign;
  if (pAddend1->nbrLimbs < pAddend2->nbrLimbs)
  {
    tmpSign = addend1Sign;
    addend1Sign = addend2Sign;
    addend2Sign = tmpSign;
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
    if ((ctr >= 0) && (pAddend1->limbs[ctr].x < pAddend2->limbs[ctr].x))
    {
      tmpSign = addend1Sign;
      addend1Sign = addend2Sign;
      addend2Sign = tmpSign;
      pTemp = pAddend1;
      pAddend1 = pAddend2;
      pAddend2 = pTemp;
    }           // At this moment, the absolute value of addend1 is greater than
                // or equal than the absolute value of addend2.
  }
  else
  {             // Nothing to do.
  }
  nbrLimbs = pAddend2->nbrLimbs;
  ptrAddend1 = pAddend1->limbs;
  ptrAddend2 = pAddend2->limbs;
  ptrSum = pSum->limbs;
  if (addend1Sign == addend2Sign)
  {             // Both addends have the same sign. Sum their absolute values.
    unsigned int carry = 0;
    unsigned int limbValue;
    for (ctr = 0; ctr < nbrLimbs; ctr++)
    {
      carry = (carry >> BITS_PER_GROUP) + (unsigned int)ptrAddend1->x +
        (unsigned int)ptrAddend2->x;
      ptrAddend1++;
      ptrAddend2++;
      limbValue = carry & MAX_INT_NBR_U;
      ptrSum->x = (int)limbValue;
      ptrSum++;
    }
    nbrLimbs = pAddend1->nbrLimbs;
    for (; ctr < nbrLimbs; ctr++)
    {
      carry = (carry >> BITS_PER_GROUP) + (unsigned int)ptrAddend1->x;
      ptrAddend1++;
      ptrSum->x = (int)(carry & MAX_INT_NBR_U);
      ptrSum++;
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
      borrow = (borrow >> BITS_PER_INT_GROUP) + ptrAddend1->x - ptrAddend2->x;
      ptrSum->x = borrow & MAX_INT_NBR;
      ptrAddend1++;
      ptrAddend2++;
      ptrSum++;
    }
    nbrLimbs = pAddend1->nbrLimbs;
    for (; ctr < nbrLimbs; ctr++)
    {
      borrow = (borrow >> BITS_PER_INT_GROUP) + ptrAddend1->x;
      ptrSum->x = borrow & MAX_INT_NBR;
      ptrAddend1++;
      ptrSum++;
    }
    while ((nbrLimbs > 1) && (pSum->limbs[nbrLimbs - 1].x == 0))
    {     // Loop that deletes non-significant zeros.
      nbrLimbs--;
    }
  }
  pSum->nbrLimbs = nbrLimbs;
  pSum->sign = addend1Sign;
  if ((pSum->nbrLimbs == 1) && (pSum->limbs[0].x == 0))
  {          // Result is zero.
    pSum->sign = SIGN_POSITIVE;
  }
}

void BigIntAdd(const BigInteger* pAddend1, const BigInteger* pAddend2, BigInteger* pSum)
{
  InternalBigIntAdd(pAddend1, pAddend2, pSum, pAddend2->sign);
}

void BigIntNegate(const BigInteger *pSrc, BigInteger *pDest)
{
  if (pSrc != pDest)
  {
    CopyBigInt(pDest, pSrc);
  }
  BigIntChSign(pDest);
}

void BigIntSubt(const BigInteger *pMinuend, const BigInteger *pSubtrahend, BigInteger *pDifference)
{
  if (pSubtrahend->sign == SIGN_POSITIVE)
  {
    InternalBigIntAdd(pMinuend, pSubtrahend, pDifference, SIGN_NEGATIVE);
  }
  else
  {
    InternalBigIntAdd(pMinuend, pSubtrahend, pDifference, SIGN_POSITIVE);
  }
}

enum eExprErr BigIntMultiply(const BigInteger *pFact1, const BigInteger *pFact2, BigInteger *pProduct)
{
  const BigInteger* pFactor1 = pFact1;
  const BigInteger* pFactor2 = pFact2;
  int nbrLimbsFactor1 = pFactor1->nbrLimbs;
  int nbrLimbsFactor2 = pFactor2->nbrLimbs;
  int nbrLimbs;
  const BigInteger *temp;
  if ((pFactor1->nbrLimbs == 1) || (pFactor2->nbrLimbs == 1))
  {       // At least one the factors has only one limb.
    int factor2;
    if (pFactor1->nbrLimbs == 1)
    {     // Force the second factor to have only one limb.
      temp = pFactor1;
      pFactor1 = pFactor2;
      pFactor2 = temp;
    }
      // Multiply BigInteger by integer.
    factor2 = ((pFactor2->sign == SIGN_POSITIVE)? pFactor2->limbs[0].x : -pFactor2->limbs[0].x);
    multint(pProduct, pFactor1, factor2);
    return EXPR_OK;
  }
#ifdef FACTORIZATION_APP
  // 2^664380 ~ 10^200000
  if ((pFactor1->nbrLimbs + pFactor2->nbrLimbs) > ((664380 / BITS_PER_GROUP) + 1))
#else
  // 2^66438 ~ 10^20000
  if ((pFactor1->nbrLimbs + pFactor2->nbrLimbs) > ((66438 / BITS_PER_GROUP) + 1))
#endif
  {
    return EXPR_INTERM_TOO_HIGH;
  }
  if (((nbrLimbsFactor1 * nbrLimbsFactor1) < nbrLimbsFactor2) ||
    ((nbrLimbsFactor2 * nbrLimbsFactor2) < nbrLimbsFactor1))
  {    // One of the factors is a lot smaller than the other.
       // Use classical multiplication.
   // return EXPR_OK
  }
  multiplyWithBothLen(&pFactor1->limbs[0], &pFactor2->limbs[0], &pProduct->limbs[0],
    nbrLimbsFactor1, nbrLimbsFactor2, &nbrLimbs);
  nbrLimbs = nbrLimbsFactor1 + nbrLimbsFactor2;
  if (pProduct->limbs[nbrLimbs - 1].x == 0)
  {
    nbrLimbs--;
  }
  pProduct->nbrLimbs = nbrLimbs;
  if ((nbrLimbs == 1) && (pProduct->limbs[0].x == 0))
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

enum eExprErr BigIntRemainder(const BigInteger *pDividend,
  const BigInteger *pDivisor, BigInteger *pRemainder)
{
  enum eExprErr rc;
  if (BigIntIsZero(pDivisor))
  {   // If divisor = 0, then remainder is the dividend.
    CopyBigInt(pRemainder, pDividend);
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
  unsigned long long ullValue;
  int nbrLimbs = 0;
  unsigned int limbValue;
  if (value >= 0)
  {
    bigint->sign = SIGN_POSITIVE;
    ullValue = (unsigned long long)value;
  }
  else
  {
    bigint->sign = SIGN_NEGATIVE;
    ullValue = (unsigned long long)(-value);
  }
  do
  {
    limbValue = (unsigned int)ullValue & MAX_VALUE_LIMB;
    bigint->limbs[nbrLimbs].x = (int)limbValue;
    nbrLimbs++;
    ullValue >>= BITS_PER_GROUP;
  } while (ullValue != 0U);
  bigint->nbrLimbs = nbrLimbs;
}

void expBigNbr(BigInteger *bignbr, double logar)
{
  unsigned int mostSignificantLimb;
  int nbrBits;
  int sizeBytes;
  bignbr->sign = SIGN_POSITIVE;
  bignbr->nbrLimbs = (int)floor(logar / (double)BITS_PER_GROUP);
  nbrBits = BITS_PER_GROUP * bignbr->nbrLimbs;
  mostSignificantLimb = (unsigned int)floor(exp(logar - ((double)nbrBits * LOG_2)) + 0.5);
  if (mostSignificantLimb == LIMB_RANGE)
  {
    mostSignificantLimb = 1;
    bignbr->nbrLimbs++;
  }
  bignbr->nbrLimbs++;
  sizeBytes = bignbr->nbrLimbs * (int)sizeof(limb);
  (void)memset(bignbr->limbs, 0, sizeBytes);
  bignbr->limbs[bignbr->nbrLimbs-1].x = mostSignificantLimb;
}

double logBigNbr(const BigInteger *pBigNbr)
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
    int nbrBits;
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
    nbrBits = (nbrLimbs - 2) * BITS_PER_GROUP;
    logar += (double)nbrBits * LOG_2;
  }
  return logar;
}

double logLimbs(const limb *pBigNbr, int nbrLimbs)
{
  double logar;
  int nbrBits;
  if (nbrLimbs > 1)
  {
    nbrBits = (nbrLimbs - 2) * BITS_PER_GROUP;
    logar = log((double)(pBigNbr + nbrLimbs - 2)->x +
      ((double)(pBigNbr + nbrLimbs - 1)->x * (double)LIMB_RANGE)) +
      (double)nbrBits * LOG_2;
  }
  else
  {
    nbrBits = (nbrLimbs - 1) * BITS_PER_GROUP;
    logar = log((double)((pBigNbr + nbrLimbs - 1)->x)) +
      (double)nbrBits * LOG_2;
  }
  return logar;
}

enum eExprErr BigIntPowerIntExp(const BigInteger *pBase, int exponent, BigInteger *pPower)
{
  double base;
  enum eExprErr rc;
  if (BigIntIsZero(pBase))
  {     // Base = 0 -> power = 0
    pPower->limbs[0].x = 0;
    pPower->nbrLimbs = 1;
    pPower->sign = SIGN_POSITIVE;
    return EXPR_OK;
  }
  base = logBigNbr(pBase);
#ifdef FACTORIZATION_APP
  if (base*(double)exponent > 460510)
#else
  if (base*(double)exponent > 46051)
#endif
  {   // More than 20000 digits. 46051 = log(10^20000)
    return EXPR_INTERM_TOO_HIGH;
  }
  CopyBigInt(&Base, pBase);
  pPower->sign = SIGN_POSITIVE;
  pPower->nbrLimbs = 1;
  pPower->limbs[0].x = 1;
  for (unsigned int mask = HALF_INT_RANGE_U; mask != 0U; mask >>= 1)
  {
    if (((unsigned int)exponent & mask) != 0U)
    {
      for (unsigned int mask2 = mask; mask2 != 0U; mask2 >>= 1)
      {
        rc = BigIntMultiply(pPower, pPower, pPower);
        if (rc != EXPR_OK)
        {
          return rc;
        }
        if (((unsigned int)exponent & mask2) != 0U)
        {
          rc = BigIntMultiply(pPower, &Base, pPower);
          if (rc != EXPR_OK)
          {
            return rc;
          }
        }
      }
      break;
    }
  }
  return EXPR_OK;
}

enum eExprErr BigIntPower(const BigInteger *pBase, const BigInteger *pExponent, BigInteger *pPower)
{
  int exponent;
  if (pExponent->sign == SIGN_NEGATIVE)
  {    // Negative exponent not accepted.
    return EXPR_INVALID_PARAM;
  }
  if (pExponent->nbrLimbs > 2)
  {     // Exponent too high.
    if ((pBase->nbrLimbs == 1) && (pBase->limbs[0].x < 2))
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
    exponent = pExponent->limbs[0].x + (pExponent->limbs[1].x << BITS_PER_GROUP);
  }
  else
  {
    exponent = pExponent->limbs[0].x;
  }
  return BigIntPowerIntExp(pBase, exponent, pPower);
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
    ptrLimb->x = (int)(carry >> 1);
    ptrLimb--;
    carry &= 1;
  }
  if ((nbrLimbs > 1) && (pArg->limbs[nbrLimbs - 1].x == 0))
  {     // Most significant limb is zero, so reduce size by one limb.
    pArg->nbrLimbs--;
  }
}

static void BigIntMutiplyPower2(BigInteger *pArg, int powerOf2)
{
  int power2 = powerOf2;
  int ctr;
  int nbrLimbs = pArg->nbrLimbs;
  limb *ptrLimbs = pArg->limbs;
  for (; power2 > 0; power2--)
  {
    unsigned int carry = 0U;
    for (ctr = 0; ctr < nbrLimbs; ctr++)
    {
      unsigned int unsignedLimb;
      carry += (unsigned int)(ptrLimbs + ctr)->x << 1;
      unsignedLimb = carry & MAX_VALUE_LIMB;
      (ptrLimbs + ctr)->x = (int)unsignedLimb;
      carry >>= BITS_PER_GROUP;
    }
    if (carry != 0U)
    {
      (ptrLimbs + ctr)->x = (int)carry;
      nbrLimbs++;
    }
  }
  pArg->nbrLimbs = nbrLimbs;
}

bool TestBigNbrEqual(const BigInteger *pNbr1, const BigInteger *pNbr2)
{
  const limb *ptrLimbs1 = pNbr1->limbs;
  const limb *ptrLimbs2 = pNbr2->limbs;
  if (pNbr1->nbrLimbs != pNbr2->nbrLimbs)
  {        // Sizes of numbers are different.
    return false;
  }
  if (pNbr1->sign != pNbr2->sign)
  {        // Sign of numbers are different.
    if ((pNbr1->nbrLimbs == 1) && (pNbr1->limbs[0].x == 0) && (pNbr2->limbs[0].x == 0))
    {              // Both numbers are zero.
      return true;
    }
    return false;
  }

           // Check whether both numbers are equal.
  for (int ctr = pNbr1->nbrLimbs - 1; ctr >= 0; ctr--)
  {
    if ((ptrLimbs1 + ctr)->x != (ptrLimbs2 + ctr)->x)
    {      // Numbers are different.
      return false;
    }
  }        // Numbers are equal.
  return true;
}

void BigIntGcd(const BigInteger *pArg1, const BigInteger *pArg2, BigInteger *pResult)
{
  int power2;
  if (BigIntIsZero(pArg1))
  {               // First argument is zero, so the GCD is second argument.
    CopyBigInt(pResult, pArg2);
    return;
  }
  if (BigIntIsZero(pArg2))
  {               // Second argument is zero, so the GCD is first argument.
    CopyBigInt(pResult, pArg1);
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
  limb* ptrLimbs = pLimbs;
  int nbrLimbs = *pNbrLimbs;
  ptrLimbs->x += addend;
  if ((unsigned int)ptrLimbs->x < LIMB_RANGE)
  {     // No overflow. Go out of routine.
    return;
  }
  ptrLimbs->x -= (int)LIMB_RANGE;
  for (int ctr = 1; ctr < nbrLimbs; ctr++)
  {
    ptrLimbs++;        // Point to next most significant limb.
    if (ptrLimbs->x != MAX_INT_NBR)
    {   // No overflow. Go out of routine.
      (ptrLimbs->x)++;   // Add carry.
      return;
    }
    ptrLimbs->x = 0;
  }
  (*pNbrLimbs)++;        // Result has an extra limb.
  (ptrLimbs + 1)->x = 1;   // Most significant limb must be 1.
}

static void subtFromAbsValue(limb *pLimbs, int *pNbrLimbs, int subt)
{
  int nbrLimbs = *pNbrLimbs;
  pLimbs->x -= subt;
  if (pLimbs->x < 0)
  {
    int ctr = 0;
    do
    {      // Loop that adjust number if there is borrow.
      (pLimbs + ctr)->x += (int)LIMB_RANGE;
      ctr++;
      if (ctr == nbrLimbs)
      {    // All limbs processed. Exit loop.
        break;
      }
    } while (--((pLimbs + ctr)->x) < 0);   // Continue loop if there is borrow.
    if ((nbrLimbs > 1) && ((pLimbs + nbrLimbs - 1)->x == 0))
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
  double dDivisor = (double)divisor;
  double dInvDivisor = 1.0 / dDivisor;
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
  if (divisor == 2)
  {      // Use shifts for divisions by 2.
    limb* ptrDest = pBigInt->limbs;
    int curLimb = ptrDest->x;
    for (int ctr = 1; ctr < nbrLimbs; ctr++)
    {  // Process starting from least significant limb.
      int nextLimb = (ptrDest + 1)->x;
      ptrDest->x = ((curLimb >> 1) | (nextLimb << (BITS_PER_GROUP - 1))) & MAX_INT_NBR;
      ptrDest++;
      curLimb = nextLimb;
    }
    ptrDest->x = (curLimb >> 1) & MAX_INT_NBR;
  }
  else
  {
    int remainder = 0;
    limb* pLimbs = pBigInt->limbs + nbrLimbs - 1;
    // Divide number by divisor.
    for (int ctr = nbrLimbs - 1; ctr >= 0; ctr--)
    {
      unsigned int dividend = ((unsigned int)remainder << BITS_PER_INT_GROUP) + 
        (unsigned int)pLimbs->x;
      double dDividend = ((double)remainder * dLimb) + (double)pLimbs->x;
      double dQuotient = (dDividend * dInvDivisor) + 0.5;
      unsigned int quotient = (unsigned int)dQuotient;   // quotient has correct value or 1 more.
      remainder = dividend - (quotient * divisor);
      if (remainder < 0)
      {     // remainder not in range 0 <= remainder < divisor. Adjust.
        quotient--;
        remainder += divisor;
      }
      pLimbs->x = (int)quotient;
      pLimbs--;
    }
  }
  if ((nbrLimbs > 1) && (pBigInt->limbs[nbrLimbs - 1].x == 0))
  {   // Most significant limb is now zero, so discard it.
    nbrLimbs--;
  }
  pBigInt->nbrLimbs = nbrLimbs;
}

int getRemainder(const BigInteger *pBigInt, int divisor)
{
  int remainder = 0;
  int nbrLimbs = pBigInt->nbrLimbs;
  double dDivisor = (double)divisor;
  double dLimb = 0x80000000;
  const limb *pLimb = &pBigInt->limbs[nbrLimbs - 1];
  for (int ctr = nbrLimbs - 1; ctr >= 0; ctr--)
  {
    int dividend = (remainder << BITS_PER_INT_GROUP) + pLimb->x;
    double dDividend = ((double)remainder * dLimb) + (double)pLimb->x;
    double dQuotient = floor((dDividend / dDivisor) + 0.5);
    int quotient = (int)(unsigned int)dQuotient;   // quotient has correct value or 1 more.
    remainder = dividend - (quotient * divisor);
    if ((unsigned int)remainder >= (unsigned int)divisor)
    {     // remainder not in range 0 <= remainder < divisor. Adjust.
      quotient--;
      remainder += divisor;
    }
    pLimb--;
  }
  if ((pBigInt->sign == SIGN_NEGATIVE) && (remainder != 0))
  {
    remainder = divisor - remainder;
  }
  return remainder;
}

void addbigint(BigInteger *pResult, int addend)
{
  int intAddend = addend;
  enum eSign sign;
  int nbrLimbs = pResult->nbrLimbs;
  limb *pResultLimbs = pResult->limbs;
  sign = pResult->sign;
  if (intAddend < 0)
  {
    intAddend = -intAddend;
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
    addToAbsValue(pResultLimbs, &nbrLimbs, intAddend);
  }
  else
  {  // Subtract addend from absolute value of pResult.
    if (nbrLimbs == 1)
    {
      pResultLimbs->x -= intAddend;
      if (pResultLimbs->x < 0)
      {
        pResultLimbs->x = -pResultLimbs->x;
        BigIntNegate(pResult, pResult);
      }
    }
    else
    {     // More than one limb.
      subtFromAbsValue(pResultLimbs, &nbrLimbs, intAddend);
    }
  }
  pResult->nbrLimbs = nbrLimbs;
}

void multint(BigInteger *pResult, const BigInteger *pMult, int factor)
{
#ifdef _USING64BITS_
  int64_t carry;
#else
  int carry;
  double dFactor;
  double dVal = 1.0 / (double)LIMB_RANGE;
#endif
  int intMult = factor;
  bool factorPositive = true;
  int nbrLimbs = pMult->nbrLimbs;
  const limb *pLimb = pMult->limbs;
  limb *pResultLimb = pResult->limbs;
  if (intMult == 0)
  {   // Any number multiplied by zero is zero.
    intToBigInteger(pResult, 0);
    return;
  }
  if (intMult < 0)
  {     // If factor is negative, indicate it and compute its absolute value.
    factorPositive = false;
    intMult = -intMult;
  }
#ifndef _USING64BITS_
  dFactor = (double)intMult;
#endif
  carry = 0;
  for (int ctr = 0; ctr < nbrLimbs; ctr++)
  {
    unsigned int unsignedLimb;
#ifdef _USING64BITS_
    carry += (int64_t)pLimb->x * (int64_t)intMult;
    unsignedLimb = (unsigned int)carry & MAX_VALUE_LIMB;
    pResultLimb->x = (int)unsignedLimb;
    pResultLimb++;
    carry >>= BITS_PER_GROUP;
#else
    int low = ((pLimb->x * factor) + carry) & MAX_INT_NBR;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    if (low < HALF_INT_RANGE)
    {
      carry = (int)((((double)(pLimb->x) * dFactor) + (double)carry +
        (double)FOURTH_INT_RANGE)*dVal);
    }
    else
    {
      carry = (int)((((double)(pLimb->x) * dFactor) + (double)carry -
        (double)FOURTH_INT_RANGE)*dVal);
    }
    pResultLimb->x = low;
    pResultLimb++;
#endif
    pLimb++;
  }
  if (carry != 0)
  {
    pResultLimb->x = (int)carry;
    nbrLimbs++;
  }
  pResult->nbrLimbs = nbrLimbs;
  pResult->sign = pMult->sign;
  if (!factorPositive)
  {
    BigIntNegate(pResult, pResult);
  }
}

void multadd(BigInteger *pResult, int iMult, const BigInteger *pMult, int addend)
{
  multint(pResult, pMult, iMult);
  addbigint(pResult, addend);
}

// Compute *pResult -> *pMult1 * iMult1 + *pMult2 * iMult2
// Use Temp as temporary variable.
void addmult(BigInteger *pResult, const BigInteger *pMult1, int iMult1, 
  const BigInteger *pMult2, int iMult2)
{
  multint(pResult, pMult1, iMult1);
  multint(&Temp, pMult2, iMult2);
  BigIntAdd(pResult, &Temp, pResult);
}

// Get number of bits of given big integer.
int bitLength(const BigInteger *pBigNbr)
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
  int exponent = Expon;
  unsigned int power = 1;
  unsigned int square = (unsigned int)NbrMod;
  while (exponent != 0)
  {
    if ((exponent & 1) == 1)
    {
      power = (power * square) % (unsigned int)currentPrime;
    }
    square = (square * square) % (unsigned int)currentPrime;
    exponent >>= 1;
  }
  return (int)power;
}

void IntArray2BigInteger(const int *ptrValues, BigInteger *bigint)
{
  const int* piValues = ptrValues;
  limb *destLimb = bigint->limbs;
  int nbrLimbs = *piValues;
  piValues++;
  if (nbrLimbs > 0)
  {
    bigint->sign = SIGN_POSITIVE;
  }
  else
  {
    bigint->sign = SIGN_NEGATIVE;
    nbrLimbs = -nbrLimbs;
  }
  if (NumberLength == 1)
  {
    destLimb->x = *piValues;
    bigint->nbrLimbs = 1;
  }
  else
  {
    int ctr;
    bigint->nbrLimbs = nbrLimbs;
    for (ctr = 0; ctr < nbrLimbs; ctr++)
    {
      destLimb->x = *piValues;
      destLimb++;
      piValues++;
    }
    for (; ctr < NumberLength; ctr++)
    {
      destLimb->x = 0;
      destLimb++;
    }
  }
}

void BigInteger2IntArray(/*@out@*/int *ptrValues, const BigInteger *bigint)
{
  int* pValues = ptrValues;
  const limb *srcLimb = bigint->limbs;
  if (NumberLength == 1)
  {
    *pValues = ((bigint->sign == SIGN_POSITIVE)? 1: -1);
    *(pValues + 1) = srcLimb->x;
  }
  else
  {
    int nbrLimbs;
    nbrLimbs = getNbrLimbs(srcLimb);
    *pValues = ((bigint->sign == SIGN_POSITIVE)? nbrLimbs : -nbrLimbs);
    pValues++;
    for (int ctr = 0; ctr < nbrLimbs; ctr++)
    {
      *pValues = srcLimb->x;
      pValues++;
      srcLimb++;
    }
  }
}

void UncompressLimbsBigInteger(const limb *ptrValues, /*@out@*/BigInteger *bigint)
{
  if (NumberLength == 1)
  {
    bigint->limbs[0].x = ptrValues->x;
    bigint->nbrLimbs = 1;
  }
  else
  {
    int nbrLimbs;
    const limb *ptrValue1;
    int numberLengthBytes = NumberLength * (int)sizeof(limb);
    (void)memcpy(bigint->limbs, ptrValues, numberLengthBytes);
    ptrValue1 = ptrValues + NumberLength;
    for (nbrLimbs = NumberLength; nbrLimbs > 1; nbrLimbs--)
    {
      ptrValue1--;
      if (ptrValue1->x != 0)
      {
        break;
      }
    }
    bigint->nbrLimbs = nbrLimbs;
  }
}

void CompressLimbsBigInteger(/*@out@*/limb *ptrValues, const BigInteger *bigint)
{
  if (NumberLength == 1)
  {
    ptrValues->x = bigint->limbs[0].x;
  }
  else
  {
    int numberLengthBytes = NumberLength * (int)sizeof(limb);
    int nbrLimbs = bigint->nbrLimbs;
    if (nbrLimbs > NumberLength)
    {
      (void)memcpy(ptrValues, bigint->limbs, numberLengthBytes);
    }
    else
    {
      int nbrLimbsBytes = nbrLimbs * (int)sizeof(limb);
      (void)memcpy(ptrValues, bigint->limbs, nbrLimbsBytes);
      nbrLimbsBytes = numberLengthBytes - nbrLimbsBytes;
      (void)memset(ptrValues + nbrLimbs, 0, nbrLimbsBytes);
    }
  }
}

void LenAndLimbs2ArrLimbs(const int *ptrValues, /*@out@*/limb *bigint, int nbrLen)
{
  int nbrLimbs = *ptrValues;
  int nbrLimbsBytes;
  if (nbrLimbs < 0)
  {
    nbrLimbs = -nbrLimbs;
  }
  nbrLimbsBytes = nbrLimbs * (int)sizeof(limb);
  (void)memcpy(bigint, ptrValues+1, nbrLimbsBytes);
  if (nbrLen > nbrLimbs)
  {
    int nbrLenBytes = nbrLen * (int)sizeof(limb);
    nbrLimbsBytes = nbrLenBytes - nbrLimbsBytes;
    (void)memset(bigint + nbrLimbs, 0, nbrLimbsBytes);
  }
}

void ArrLimbs2LenAndLimbs(/*@out@*/int *ptrValues, const limb *bigint, int nbrLen)
{
  int nbrLimbs;
  int nbrLimbsBytes = (nbrLen - 1) * (int)sizeof(limb);
  (void)memcpy(ptrValues+1, bigint, nbrLimbsBytes);
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
int PowerCheck(const BigInteger *pBigNbr, BigInteger *pBase)
{
  limb *ptrLimb;
  double dN;
  int nbrLimbs;
  int maxExpon;
  int h;
  int j;
  int modulus;
  int intLog2root;
  int primesLength;
  int Exponent;
  int base = 0;
  double log2N;
  double log2root;
  double dMaxExpon;
  int prime2310x1[] =
  { 2311, 4621, 9241, 11551, 18481, 25411, 32341, 34651, 43891, 50821 };
  // Primes of the form 2310x+1.
  bool expon2 = true;
  bool expon3 = true;
  bool expon5 = true;
  bool expon7 = true;
  bool expon11 = true;
  double dLogBigNbr = logBigNbr(pBigNbr);
  if (pBigNbr->nbrLimbs > 10)
  {
    for (base = 2; base <= 100; base++)
    {     // Check whether pBigNbr is perfect power of these bases.
      int nbrLimbsBytes;
      double dProd;
      double dLogBase = log(base);
      double dExponent = dLogBigNbr / dLogBase;
      Exponent = (int)dExponent - 1;
      dProd = dLogBigNbr - ((double)Exponent * dLogBase);
      if ((dProd > 0.00000000001) || (dProd < -0.00000000001))
      {
        Exponent++;
        dProd -= dLogBase;
        if ((dProd > 0.00000000001) || (dProd < -0.00000000001))
        {
          Exponent++;
          dProd -= dLogBase;
          if ((dProd > 0.00000000001) || (dProd < -0.00000000001))
          {
            continue;           // Test next base.
          }
        }
      }
      intToBigInteger(pBase, base);
      (void)BigIntPowerIntExp(pBase, Exponent, &Temp3);
      if (BigIntEqual(&Temp3, pBigNbr))
      {
        return Exponent;
      }
      pBase->nbrLimbs = pBigNbr->nbrLimbs;
      nbrLimbsBytes = pBase->nbrLimbs * (int)sizeof(limb);
      (void)memcpy(pBase->limbs, pBigNbr->limbs, nbrLimbsBytes);
      return 1;
    }
    dMaxExpon = (dLogBigNbr / log(101)) + 0.5;
  }
  else
  {
    dMaxExpon = (dLogBigNbr / LOG_2) + 0.5;
  }
  maxExpon = (int)dMaxExpon;
  for (h = 0; h < sizeof(prime2310x1) / sizeof(prime2310x1[0]); h++)
  {
    int testprime = prime2310x1[h];
    int mod = getRemainder(pBigNbr, testprime);
    if (expon2 && (intModPow(mod, testprime / 2, testprime) > 1))
    {
      expon2 = false;
    }
    if (expon3 && (intModPow(mod, testprime / 3, testprime) > 1))
    {
      expon3 = false;
    }
    if (expon5 && (intModPow(mod, testprime / 5, testprime) > 1))
    {
      expon5 = false;
    }
    if (expon7 && (intModPow(mod, testprime / 7, testprime) > 1))
    {
      expon7 = false;
    }
    if (expon11 && (intModPow(mod, testprime / 11, testprime) > 1))
    {
      expon11 = false;
    }
  }
  primesLength = (2 * maxExpon) + 3;
  for (h = 2; h <= maxExpon; h++)
  {
    ProcessExpon[h] = true;
  }
  for (h = 2; h < primesLength; h++)
  {
    primes[h] = true;
  }
  for (h = 2; (h * h) < primesLength; h++)
  { // Generation of primes
    for (j = h * h; j < primesLength; j += h)
    { // using Eratosthenes sieve
      primes[j] = false;
    }
  }
  for (h = 13; h < primesLength; h++)
  {
    if (primes[h])
    {
      int processed = 0;
      for (j = (2 * h) + 1; j < primesLength; j += 2 * h)
      {
        if (primes[j])
        {
          modulus = getRemainder(pBigNbr, j);
          if (intModPow(modulus, j / h, j) > 1)
          {
            for (j = h; j <= maxExpon; j += h)
            {
              ProcessExpon[j] = false;
            }
            break;
          }
        }
        processed++;
        if (processed > 10)
        {
          break;
        }
      }
    }
  }
  log2N = dLogBigNbr / LOG_2;
  for (Exponent = maxExpon; Exponent >= 2; Exponent--)
  {
    int k;
    int prime;
    int nbrBits;
    if (((Exponent % 2) == 0) && !expon2)
    {
      continue; // Not a square
    }
    if (((Exponent % 3) == 0) && !expon3)
    {
      continue; // Not a cube
    }
    if (((Exponent % 5) == 0) && !expon5)
    {
      continue; // Not a fifth power
    }
    if (((Exponent % 7) == 0) && !expon7)
    {
      continue; // Not a 7th power
    }
    if (((Exponent % 11) == 0) && !expon11)
    {
      continue; // Not an 11th power
    }
    if (!ProcessExpon[Exponent])
    {
      continue;
    }
    // Initialize approximation to n-th root (n = Exponent).
    log2root = log2N / (double)Exponent;
    intLog2root = (int)floor(log2root/ (double)BITS_PER_GROUP);
    nbrLimbs = intLog2root + 1;
    ptrLimb = &pBase->limbs[nbrLimbs - 1];
    nbrBits = intLog2root * BITS_PER_GROUP;
    dN = exp((log2root - (double)nbrBits) * LOG_2);
    if (nbrLimbs == 1)
    {
      double dQuot;
      base = (int)dN - 1;
      dQuot = dN / (double)base;
      if ((dQuot > 1.0000000001) || (dQuot < 0.9999999999))
      {
        base++;
        dQuot = dN / (double)base;
        if ((dQuot > 1.0000000001) || (dQuot < 0.9999999999))
        {
          base++;
          dQuot = dN / (double)base;
          if ((dQuot > 1.0000000001) || (dQuot < 0.9999999999))
          {
            continue;   // Exponent is incorrect. Check next one.
          }
        }
      }
    }
    // If p = prime = k*expon+1, and n = r (mod p), it must be r^k = 1 (mod p)
    k = 1;
    for (;;)
    {
      prime = (k * Exponent) + 1;
      // Test that prime is really prime.
      for (j = 2; (j*j) <= prime; j++)
      {
        if ((prime / j * j) == prime)
        {   // Number is not prime.
          break;
        }
      }
      if ((j*j) > prime)
      {     // Number is prime.
        break;
      }
      k++;
    }
    if (base != prime)
    {     // Prime just generated is different from the number to check primality.
      if (intModPow(getRemainder(pBigNbr, prime), k, prime) != 1)
      {   // Number is not a power of the type a^Exponent.
        continue;
      }
    }
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
      double dLimbRange = (double)LIMB_RANGE;
      dN += 1.0 / dLimbRange;
      ptrLimb->x = (int)trunc(dN);
      dN -= trunc(dN);
      (ptrLimb - 1)->x = (int)trunc(dN * dLimbRange);
    }
    pBase->nbrLimbs = nbrLimbs;
    // Perform Newton iteration for n-th root.
    for (;;)
    {   // Check whether the approximate root is actually exact.
      (void)BigIntPowerIntExp(pBase, Exponent-1, &Temp3); // Temp3 <- x^(e-1)
      (void)BigIntMultiply(&Temp3, pBase, &Temp2);        // Temp2 <- x^e 
      BigIntSubt(pBigNbr, &Temp2, &Temp2);            // Compare to radicand.
      if (BigIntIsZero(&Temp2))
      {                     // Perfect power, so go out.
        return Exponent;
      }
      if (Temp2.sign == SIGN_POSITIVE)
      {                     // x^e > radicand -> not perfect power, so go out.
        break;
      }
      (void)BigIntDivide(pBigNbr, &Temp3, &Temp);         // Temp -> N/x^(e-1)
      BigIntSubt(&Temp, pBase, &Temp2);             // Temp2 -> N/x^(e-1) - x
      if (BigIntIsZero(&Temp2))
      {     // New approximation will be the same as previous. Go out.
        break;
      }
      subtractdivide(&Temp2, Exponent - 1, Exponent);   // Compute (Temp2 - (1 - Exponent)) / Exponent
      BigIntAdd(&Temp2, pBase, pBase);
    }
  }
  pBase->nbrLimbs = pBigNbr->nbrLimbs;
  (void)memcpy(pBase->limbs, pBigNbr->limbs, pBase->nbrLimbs*sizeof(limb));
  return 1;
}

bool checkOne(const limb *value, int nbrLimbs)
{
  const limb *ptrValue = value;
  for (int idx = 0; idx < nbrLimbs; idx++)
  {
    if (ptrValue->x != MontgomeryMultR1[idx].x)
    {
      return false;    // Go out if value is not 1 (mod p)
    }
    ptrValue++;
  }
  return true;
}

bool checkMinusOne(const limb *value, int nbrLimbs)
{
  const limb* limbValue = value;
  unsigned int carry;
  carry = 0;
  for (int idx = 0; idx < nbrLimbs; idx++)
  {
    carry += (unsigned int)limbValue->x + (unsigned int)MontgomeryMultR1[idx].x;
    limbValue++;
    if ((carry & MAX_VALUE_LIMB) != (unsigned int)TestNbr[idx].x)
    {
      return false;    // Go out if value is not -1 (mod p)
    }
    carry >>= BITS_PER_GROUP;
  }
  return true;
}

void BigIntDivideBy2(BigInteger *nbr)
{
  int nbrLimbs = nbr->nbrLimbs;
  int *ptrDest = &nbr->limbs[0].x;
  int curLimb = *ptrDest;
  for (int ctr = 1; ctr < nbrLimbs; ctr++)
  {  // Process starting from least significant limb.
    int nextLimb = *(ptrDest + 1);
    *ptrDest = ((curLimb >> 1) | (nextLimb << (BITS_PER_GROUP - 1))) & MAX_INT_NBR;
    ptrDest++;
    curLimb = nextLimb;
  }
  *ptrDest = (curLimb >> 1) & MAX_INT_NBR;
  if ((nbrLimbs > 1) && (nbr->limbs[nbrLimbs - 1].x == 0))
  {
    nbr->nbrLimbs--;
  }
}

void BigIntMultiplyBy2(BigInteger *nbr)
{
  int nbrLimbs;
  int prevLimb;
  limb *ptrDest = &nbr->limbs[0];
  nbrLimbs = nbr->nbrLimbs;
  prevLimb = 0;
  for (int ctr = 0; ctr < nbrLimbs; ctr++)
  {  // Process starting from least significant limb.
    int curLimb = ptrDest->x;
    ptrDest->x = ((curLimb << 1) | (prevLimb >> (BITS_PER_GROUP - 1))) & MAX_INT_NBR;
    ptrDest++;
    prevLimb = curLimb;
  }
  if ((prevLimb & HALF_INT_RANGE) != 0)
  {
    ptrDest->x = 1;
    nbr->nbrLimbs++;
  }
}

// Find power of 2 that divides the number.
// output: pNbrLimbs = pointer to number of limbs
//         pShRight = pointer to power of 2.
void DivideBigNbrByMaxPowerOf2(int *pShRight, limb *number, int *pNbrLimbs)
{
  int power2 = 0;
  int index;
  int index2;
  int shRg;
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
  if (index == nbrLimbs)
  {   // Input number is zero.
    *pShRight = power2;
    return;
  }
  for (unsigned int mask = 1U; mask <= MAX_VALUE_LIMB; mask *= 2)
  {
    if (((unsigned int)number[index].x & mask) != 0U)
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
  for (index2 = index; index2 < (nbrLimbs-1); index2++)
  {
    number[index2].x = ((number[index2].x >> shRg) |
                        (number[index2+1].x << (BITS_PER_GROUP - shRg))) &
                        MAX_VALUE_LIMB;
  }
  if (index2 < nbrLimbs)
  {
    number[index2].x = (number[index2].x >> shRg) & MAX_VALUE_LIMB;
  }
  if (index > 0)
  {   // Move limbs to final position.
    (void)memmove(number, &number[index], (nbrLimbs - index) * sizeof(limb));
  }
  *pShRight = power2;
}

// Calculate Jacobi symbol by following algorithm 2.3.5 of C&P book.
int JacobiSymbol(int upper, int lower)
{
  int a = upper % lower;
  int m = lower;
  int t = 1;
  while (a != 0)
  {
    int tmp;
    while ((a & 1) == 0)
    {     // a is even.
      a >>= 1;
      if (((m & 7) == 3) || ((m & 7) == 5))
      {   // m = 3 or m = 5 (mod 8)
        t = -t;
      }
    }
    tmp = a; a = m; m = tmp;   // Exchange a and m.
    if ((a & m & 3) == 3)
    {   // a = 3 and m = 3 (mod 4)
      t = -t;
    }
    a = a % m;
  }
  if ((m == 1) || (m == -1))
  {
    return t;
  }
  return 0;
}

int BigIntJacobiSymbol(const BigInteger *upper, const BigInteger *lower)
{
  int t;
  int power2;
  static BigInteger a;
  static BigInteger m;
  static BigInteger tmp;
  CopyBigInt(&m, lower);               // m <- lower
  DivideBigNbrByMaxPowerOf2(&power2, m.limbs, &m.nbrLimbs);
  (void)BigIntRemainder(upper, lower, &a);   // a <- upper % lower
  t = 1;
  if (upper->sign == SIGN_NEGATIVE)
  {
    a.sign = SIGN_POSITIVE;
    if ((m.limbs[0].x & 3) == 3)
    {
      t = -1;
    }
  }
  while (!BigIntIsZero(&a))             // a != 0
  {
    while ((a.limbs[0].x & 1) == 0)
    {     // a is even.
      BigIntDivideBy2(&a);              // a <- a / 2
      if (((m.limbs[0].x & 7) == 3) || ((m.limbs[0].x & 7) == 5))
      {   // m = 3 or m = 5 (mod 8)
        t = -t;
      }
    }
    CopyBigInt(&tmp, &a);               // Exchange a and m.
    CopyBigInt(&a, &m);
    CopyBigInt(&m, &tmp);
    if ((a.limbs[0].x & m.limbs[0].x & 3) == 3)
    {   // a = 3 and m = 3 (mod 4)
      t = -t;
    }
    (void)BigIntRemainder(&a, &m, &tmp);
    CopyBigInt(&a, &tmp);              // a <- a % m  
  }
  if ((m.nbrLimbs == 1) && (m.limbs[0].x == 1))
  {              // Absolute value of m is 1.
    return t;
  }
  return 0;
}

static void Halve(limb *pValue)
{
  if ((pValue[0].x & 1) == 0)
  {    // Number to halve is even. Divide by 2.
    DivBigNbrByInt((int *)pValue, 2, (int *)pValue, NumberLength);
  }
  else
  {    // Number to halve is odd. Add modulus and then divide by 2.
    AddBigNbr((int *)pValue, (int *)TestNbr, (int *)pValue, NumberLength + 1);
    DivBigNbrByInt((int *)pValue, 2, (int *)pValue, NumberLength + 1);
  }
}

void initializeSmallPrimes(int* pSmallPrimes)
{
  int P;
  int* ptrSmallPrimes = pSmallPrimes;
  if (*ptrSmallPrimes != 0)
  {
    return;
  }
  P = 3;
  *ptrSmallPrimes = 2;
  ptrSmallPrimes++;
  for (int ctr = 1; ctr <= SMALL_PRIMES_ARRLEN; ctr++)
  {     // Loop that fills the SmallPrime array.
    int Q;
    *ptrSmallPrimes = P; /* Store prime */
    ptrSmallPrimes++;
    do
    {
      P += 2;
      for (Q = 3; (Q * Q) <= P; Q += 2)
      { /* Check if P is prime */
        if ((P % Q) == 0)
        {
          break;  /* Composite */
        }
      }
    } while ((Q * Q) <= P);
  }
}

// BPSW primality test:
// 1) If the input number is 2-SPRP composite, indicate composite and go out.
// 2) If number is perfect square, indicate it is composite and go out.
// 3) Find the first D in the sequence 5, -7, 9, -11, 13, -15, ...
//    for which the Jacobi symbol (D/n) is −1. Set P = 1 and Q = (1 - D) / 4.
// 4) Perform a strong Lucas probable prime test on n using parameters D, P,
//    and Q. If n is not a strong Lucas probable prime, then n is composite.
//    Otherwise, n is almost certainly prime.
// Output: 0 = probable prime.
//         1 = composite: not 2-Fermat pseudoprime.
//         2 = composite: does not pass 2-SPRP test.
//         3 = composite: does not pass strong Lucas test.
#if FACTORIZATION_APP
int BpswPrimalityTest(const BigInteger* pValue, const struct sFactors* pstFactors)
#else
int BpswPrimalityTest(const BigInteger *pValue)
#endif
{
#if defined(__EMSCRIPTEN__) && !defined(FACTORIZATION_APP)
  char text[200];
#endif
#if !defined(__EMSCRIPTEN__) && defined(FACTORIZATION_APP)
  (void)pstFactors;    // Parameter is not used.
#endif
#ifdef __EMSCRIPTEN__
  char *ptrText;
#endif
  int Mult3Len;
  int ctr;
  int D;
  int absQ;
  int mult;
  int mask;
  int index;
  int signPowQ;
  bool insidePowering = false;
  int nbrLimbs = pValue->nbrLimbs;
  const limb *limbs = pValue->limbs;
  static BigInteger tmp;
  if ((nbrLimbs == 1) && (limbs->x <= 2))
  {
    return 0;    // Indicate prime.
  }
  if ((limbs->x & 1) == 0)
  {
    return 1;    // Number is even and different from 2. Indicate composite.
  }
  if (nbrLimbs > 1)
  {              // Check whether it is divisible by small number.
    initializeSmallPrimes(smallPrimes);
    for (int primeIndex = 0; primeIndex < 180; primeIndex += 3)
    {
      int primeProd = smallPrimes[primeIndex] * smallPrimes[primeIndex+1] * smallPrimes[primeIndex+2];
      int remainder = getRemainder(pValue, primeProd);
      if (((remainder % smallPrimes[primeIndex]) == 0) ||
        ((remainder % smallPrimes[primeIndex + 1]) == 0) ||
        ((remainder % smallPrimes[primeIndex + 2]) == 0))
      {
        return 1;   // Number is divisible by small number. Indicate composite.
      }
    }
  }
#ifdef __EMSCRIPTEN__
#ifdef FACTORIZATION_APP
  StepECM = 3;   // Show progress (in percentage) of BPSW primality test.
  ptrText = ShowFactoredPart(pValue, pstFactors);
  (void)strcpy(ptrText, lang ? "<p>Paso 1 del algoritmo BPSW de primos probables: Miller-Rabin fuerte con base 2.</p>" :
    "<p>Step 1 of BPSW probable prime algorithm: Strong Miller-Rabin with base 2.</p>");
  ShowLowerText();
#else
  databack(lang ? "3<p>Paso 1 del algoritmo BPSW de primos probables: Miller-Rabin fuerte con base 2.</p>" :
    "3<p>Step 1 of BPSW probable prime algorithm: Strong Miller-Rabin with base 2.</p>");
#endif
#endif
  // Perform 2-SPRP test
  (void)memcpy(q, limbs, nbrLimbs * sizeof(limb));
  q[nbrLimbs] = 0;
  q[0]--;                     // q = p - 1 (p is odd, so there is no carry).
  (void)memcpy(Mult3, q, (nbrLimbs + 1) * sizeof(q[0]));
  Mult3Len = nbrLimbs;
  DivideBigNbrByMaxPowerOf2(&ctr, Mult3, &Mult3Len);
  (void)memcpy(TestNbr, limbs, nbrLimbs * sizeof(limb));
  TestNbr[nbrLimbs].x = 0;
  GetMontgomeryParms(nbrLimbs);
  // Find Mult1 = 2^Mult3.
  (void)memcpy(Mult1, MontgomeryMultR1, (NumberLength + 1) * sizeof(limb));  // power <- 1
  for (index = Mult3Len - 1; index >= 0; index--)
  {
    int groupExp = Mult3[index].x;
#ifdef __EMSCRIPTEN__
    percentageBPSW = (Mult3Len - index) * 100 / Mult3Len;
#endif
    for (mask = 1 << (BITS_PER_GROUP - 1); mask > 0; mask >>= 1)
    {
      modmult(Mult1, Mult1, Mult1);
      if ((groupExp & mask) != 0)
      {
        modmultInt(Mult1, 2, Mult1);
      }
    }
  }
                                            // If Mult1 != 1 and Mult1 = TestNbr-1, perform full test.
  if (!checkOne(Mult1, nbrLimbs) && !checkMinusOne(Mult1, nbrLimbs))
  {
    int i;
    for (i = 0; i < ctr; i++)
    {               // Loop that squares number.
      modmult(Mult1, Mult1, Mult4);
      if (checkOne(Mult4, nbrLimbs) != 0)
      {  // Current value is 1 but previous value is not 1 or -1: composite
#if defined(__EMSCRIPTEN__) && defined(FACTORIZATION_APP)
        StepECM = 0;    // Do not show progress.
#endif
        return 2;       // Composite. Not 2-strong probable prime.
      }
      if (checkMinusOne(Mult4, nbrLimbs) != 0)
      {
        i = -1;         // Number is strong pseudoprime.
        break;
      }
      (void)memcpy(Mult1, Mult4, nbrLimbs * sizeof(limb));
    }
    if (i == ctr)
    {
#if defined(__EMSCRIPTEN__) && defined(FACTORIZATION_APP)
      StepECM = 0;      // Do not show progress.
#endif
      return 1;         // Not 2-Fermat probable prime.
    }
    if (i != -1)
    {
#if defined(__EMSCRIPTEN__) && defined(FACTORIZATION_APP)
      StepECM = 0;      // Do not show progress.
#endif
      return 2;         // Composite. Not 2-strong probable prime.
    }
  }
  // At this point, the number is 2-SPRP, so check whether the number is perfect square.
  squareRoot(pValue->limbs, tmp.limbs, pValue->nbrLimbs, &tmp.nbrLimbs);
  tmp.sign = SIGN_POSITIVE;
  (void)BigIntMultiply(&tmp, &tmp, &tmp);
  if (BigIntEqual(pValue, &tmp))
  {                  // Number is perfect square.
#if defined(__EMSCRIPTEN__) && defined(FACTORIZATION_APP)
    StepECM = 0;     // Do not show progress.
#endif
    return 3;        // Indicate number does not pass strong Lucas test.
  }
  // At this point, the number is not perfect square, so find value of D.
  mult = 1;
  D = 5;
  for (;;)
  {
    int rem = getRemainder(pValue, D);
    if (JacobiSymbol(rem, D*mult) == -1)
    {
      break;
    }
    mult = -mult;
    D += 2;
  }
  absQ = (D + 1) >> 2;
#ifdef __EMSCRIPTEN__
  int i;
#ifdef FACTORIZATION_APP
  StepECM = 3;   // Show progress (in percentage) of BPSW primality test.
  ptrText = ShowFactoredPart(pValue, pstFactors);
#else
  ptrText = text;
  *ptrText = '3';
  ptrText++;
#endif
  copyStr(&ptrText, lang ? "<p>Paso 2 del algoritmo BPSW de primos probables: Lucas fuerte con P=1, D=" :
    "<p>Step 2 of BPSW probable prime algorithm: Strong Lucas with P=1, D=");
  int2dec(&ptrText, D);
  copyStr(&ptrText, ", Q=");
  i = -absQ;        // Get value of Q to show it on screen.
  if (i < 0)
  {
    copyStr(&ptrText, "&minus;");
    i = -i;
  }
  int2dec(&ptrText, i);
  copyStr(&ptrText, "</p>");
#ifdef FACTORIZATION_APP
  ShowLowerText();
#else
  databack(text);
#endif
#endif
  // Perform strong Lucas primality test on n with parameters D, P=1, Q just found.
  // Let d*2^s = n+1 where d is odd.
  // Then U_d = 0 or v_{d*2^r} = 0 for some r < s.
  // Use the following recurrences:
  // U_0 = 0, V_0 = 2.
  // U_{2k} = U_k * V_k
  // V_{2k} = (V_k)^2 - 2*Q^K
  // U_{2k+1} = (U_{2k} + V_{2k})/2
  // V_{2k+1} = (D*U_{2k} + V_{2k})/2
  // Use the following temporary variables:
  // Mult1 for Q^n, Mult3 for U, Mult4 for V, Mult2 for temporary.
  (void)memcpy(Mult1, MontgomeryMultR1, (nbrLimbs + 1) * sizeof(limb)); // Q^0 <- 1.
  signPowQ = 1;
  (void)memset(Mult3, 0, (nbrLimbs + 1) * sizeof(limb));                // U_0 <- 0.
  (void)memcpy(Mult4, MontgomeryMultR1, (nbrLimbs + 1) * sizeof(limb)); 
  AddBigNbrMod(Mult4, Mult4, Mult4);                              // V_0 <- 2.
  CopyBigInt(&expon, pValue);
  addbigint(&expon, 1);                            // expon <- n + 1.
  Temp.limbs[nbrLimbs].x = 0;
  Temp2.limbs[nbrLimbs].x = 0;
  expon.limbs[expon.nbrLimbs].x = 0;
  DivideBigNbrByMaxPowerOf2(&ctr, expon.limbs, &expon.nbrLimbs);
  for (index = expon.nbrLimbs - 1; index >= 0; index--)
  {
#ifdef __EMSCRIPTEN__
    percentageBPSW = (expon.nbrLimbs - index) * 100 / expon.nbrLimbs;
#endif
    int groupExp = expon.limbs[index].x;
    for (mask = 1 << (BITS_PER_GROUP - 1); mask > 0; mask >>= 1)
    {
      if (insidePowering)
      {
        // U_{2k} = U_k * V_k
        // V_{2k} = (V_k)^2 - 2*Q^K
        modmult(Mult3, Mult4, Mult3);          // U <- U * V
        modmult(Mult4, Mult4, Mult4);          // V <- V * V
        if (signPowQ > 0)
        {
          SubtBigNbrMod(Mult4, Mult1, Mult4);  // V <- V - Q^k
          SubtBigNbrMod(Mult4, Mult1, Mult4);  // V <- V - Q^k
        }
        else
        {
          AddBigNbrMod(Mult4, Mult1, Mult4);   // V <- V - Q^k
          AddBigNbrMod(Mult4, Mult1, Mult4);   // V <- V - Q^k
        }
        signPowQ = 1;                          // Indicate it is positive. 
        modmult(Mult1, Mult1, Mult1);          // Square power of Q.
      }
      if ((groupExp & mask) != 0)
      {        // Bit of exponent is equal to 1.
        // U_{2k+1} = (U_{2k} + V_{2k})/2
        // V_{2k+1} = (D*U_{2k} + V_{2k})/2
        Mult3[NumberLength].x = 0;
        Mult4[NumberLength].x = 0;
        AddBigNbrMod(Mult3, Mult4, Temp.limbs);
        Halve(Temp.limbs);                     // Temp <- (U + V)/2
        MultBigNbrByIntModN((int *)Mult3, D, (int *)Temp2.limbs, (int *)TestNbr, nbrLimbs);
        if (mult > 0)
        {      // D is positive
          AddBigNbrMod(Mult4, Temp2.limbs, Mult4);
        }
        else
        {      // D is negative.
          SubtBigNbrMod(Mult4, Temp2.limbs, Mult4);
        }
        Halve(Mult4);                       // V <- (V +/- U*D)/2
        (void)memcpy(Mult3, Temp.limbs, NumberLength * sizeof(limb));
        modmultInt(Mult1, absQ, Mult1);     // Multiply power of Q by Q.
        signPowQ = -mult;                   // Attach correct sign to power.
        insidePowering = true;
      }
    }
  }
  // If U is zero, the number passes the BPSW primality test.
  if (BigNbrIsZero(Mult3))
  {
#if defined(__EMSCRIPTEN__) && defined(FACTORIZATION_APP)
    StepECM = 0;      // Do not show progress.
#endif
    return 0;         // Indicate number is probable prime.
  }
  for (index = 0; index < ctr; index++)
  {
    // If V is zero, the number passes the BPSW primality test.
    if (BigNbrIsZero(Mult4))
    {
#if defined(__EMSCRIPTEN__) && defined(FACTORIZATION_APP)
      StepECM = 0;    // Do not show progress.
#endif
      return 0;       // Indicate number is probable prime.
    }
    modmult(Mult4, Mult4, Mult4);          // V <- V * V
    if (signPowQ > 0)
    {
      SubtBigNbrMod(Mult4, Mult1, Mult4);  // V <- V - Q^k
      SubtBigNbrMod(Mult4, Mult1, Mult4);  // V <- V - Q^k
    }
    else
    {
      AddBigNbrMod(Mult4, Mult1, Mult4);   // V <- V - Q^k
      AddBigNbrMod(Mult4, Mult1, Mult4);   // V <- V - Q^k
    }
    modmult(Mult1, Mult1, Mult1);          // Square power of Q.
    signPowQ = 1;                          // Indicate it is positive.
  }
#if defined(__EMSCRIPTEN__) && defined(FACTORIZATION_APP)
  StepECM = 0;     // Do not show progress.
#endif
  return 3;        // Number does not pass strong Lucas test.
}

void NbrToLimbs(int nbr, /*@out@*/limb *limbs, int len)
{
  if (nbr >= MAX_VALUE_LIMB)
  {
    limbs->x = nbr % MAX_VALUE_LIMB;
    (limbs+1)->x = nbr / MAX_VALUE_LIMB;
    (void)memset(limbs + 2, 0, (len-2) * sizeof(limb));
  }
  else
  {
    limbs->x = nbr;
    (void)memset(limbs + 1, 0, (len-1) * sizeof(limb));
  }
}

bool BigNbrIsZero(const limb *value)
{
  const limb* ptrValue = value;
  for (int ctr = 0; ctr < NumberLength; ctr++)
  {
    if (ptrValue->x != 0)
    {
      return false;  // Number is not zero.
    }
    ptrValue++;
  }
  return true;       // Number is zero
}

bool BigIntIsZero(const BigInteger *value)
{
  if ((value->nbrLimbs == 1) && (value->limbs[0].x == 0))
  {
    return true;     // Number is zero.
  }
  return false;      // Number is not zero.
}

bool BigIntIsOne(const BigInteger* value)
{
  if ((value->nbrLimbs == 1) && (value->limbs[0].x == 1) && (value->sign == SIGN_POSITIVE))
  {
    return true;     // Number is zero.
  }
  return false;      // Number is not zero.
}

bool BigIntEqual(const BigInteger *value1, const BigInteger *value2)
{
  int nbrLimbs;
  const limb *ptrValue1;
  const limb *ptrValue2;
  if ((value1->nbrLimbs != value2->nbrLimbs) || (value1->sign != value2->sign))
  {
    return false;    // Numbers are not equal.
  }
  nbrLimbs = value1->nbrLimbs;
  ptrValue1 = value1->limbs;
  ptrValue2 = value2->limbs;
  for (int index = 0; index < nbrLimbs; index++)
  {
    if (ptrValue1->x != ptrValue2->x)
    {
      return false;  // Numbers are not equal.
    }
    ptrValue1++;
    ptrValue2++;
  }
  return true;       // Numbers are equal.
}

double getMantissa(const limb *ptrLimb, int nbrLimbs)
{
  double dN = (double)(ptrLimb - 1)->x;
  double dInvLimb = 1.0 / (double)LIMB_RANGE;
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

void BigIntPowerOf2(BigInteger *pResult, int exponent)
{
  int nbrLimbs = exponent / BITS_PER_GROUP;
  if (nbrLimbs > 0)
  {
    (void)memset(pResult->limbs, 0, nbrLimbs * sizeof(limb));
  }
  pResult->limbs[nbrLimbs].x = 1 << (exponent % BITS_PER_GROUP);
  pResult->nbrLimbs = nbrLimbs + 1;
  pResult->sign = SIGN_POSITIVE;
}

// Find power of 4 that divides the number.
// output: pNbrLimbs = pointer to number of limbs
//         pPower4 = pointer to power of 4.
void DivideBigNbrByMaxPowerOf4(int *pPower4, limb *value, int *pNbrLimbs)
{
  int powerOf4;
  int powerOf2 = 0;
  int numLimbs = *pNbrLimbs;
  int index;
  int power2gr;
  int shRg;
  limb prevLimb;
  limb currLimb;
  // Start from least significant limb (number zero).
  for (index = 0; index < numLimbs; index++)
  {
    if ((value + index)->x != 0)
    {
      break;
    }
    powerOf2 += BITS_PER_GROUP;
  }
  for (int mask = 0x1; (unsigned int)mask <= MAX_VALUE_LIMB; mask *= 2)
  {
    if (((value + index)->x & mask) != 0)
    {
      break;
    }
    powerOf2++;
  }
  powerOf4 = powerOf2 >> 1;
  // Divide value by this power.
  power2gr = powerOf2 % (2 * BITS_PER_GROUP);
  shRg = (power2gr & (-2)) % BITS_PER_GROUP; // Shift right bit counter
  if (power2gr == BITS_PER_GROUP)
  {
    index--;
  }
  prevLimb.x = 0;
  for (int index2 = numLimbs - 1; index2 >= index; index2--)
  {
    currLimb.x = (value + index2)->x;
    (value + index2)->x = ((currLimb.x >> shRg) | (prevLimb.x << (BITS_PER_GROUP - shRg))) & MAX_VALUE_LIMB;
    prevLimb.x = currLimb.x;
  }
  if (index != 0)
  {
    (void)memmove(value, value + index, (numLimbs - index) * sizeof(limb));
  }
  if ((value + numLimbs - 1)->x != 0)
  {
    *pNbrLimbs = numLimbs - index;
  }
  else
  {
    *pNbrLimbs = numLimbs - index - 1;
  }
  *pPower4 = powerOf4;
}

static void InternalBigIntLogical(const BigInteger *firstArg,
  const BigInteger *secondArg, BigInteger *result, enum eOper operation)
{
  int idx;
  int carryFirst = 0;
  int carrySecond = 0;
  int limbFirst;
  int limbSecond;
  const BigInteger *tmpptr;
  if (firstArg->nbrLimbs < secondArg->nbrLimbs)
  {    // After the exchange, firstArg has not fewer limbs than secondArg.
    tmpptr = firstArg;
    firstArg = secondArg;
    secondArg = tmpptr;
  }
  for (idx = 0; idx < secondArg->nbrLimbs; idx++)
  {
    limbFirst = firstArg->limbs[idx].x;
    limbSecond = secondArg->limbs[idx].x;
    if (firstArg->sign == SIGN_NEGATIVE)
    {
      carryFirst -= limbFirst;
      limbFirst = carryFirst & MAX_INT_NBR;
      carryFirst >>= 31;
    }
    if (secondArg->sign == SIGN_NEGATIVE)
    {
      carrySecond -= limbSecond;
      limbSecond = carrySecond & MAX_INT_NBR;
      carrySecond >>= 31;
    }
    if (operation == OPERATION_AND)
    {
      result->limbs[idx].x = limbFirst & limbSecond;
    }
    else if (operation == OPERATION_OR)
    {
      result->limbs[idx].x = limbFirst | limbSecond;
    }
    else
    {
      result->limbs[idx].x = limbFirst ^ limbSecond;
    }
  }
  if (firstArg->sign == SIGN_POSITIVE)
  {
    limbFirst = 0;
  }
  else
  {
    limbFirst = -1;
  }
  for (; idx < firstArg->nbrLimbs; idx++)
  {
    limbSecond = secondArg->limbs[idx].x;
    if (secondArg->sign == SIGN_NEGATIVE)
    {
      carrySecond -= limbSecond;
      limbSecond = carrySecond & MAX_INT_NBR;
      carrySecond >>= 31;
    }
    if (operation == OPERATION_AND)
    {
      result->limbs[idx].x = limbFirst & limbSecond;
    }
    else if (operation == OPERATION_OR)
    {
      result->limbs[idx].x = limbFirst | limbSecond;
    }
    else
    {
      result->limbs[idx].x = limbFirst ^ limbSecond;
    }
  }
  // Generate sign of result according to operation and
  // signs of arguments.
  if (operation == OPERATION_AND)
  {
    if ((firstArg->sign == SIGN_NEGATIVE) && (secondArg->sign == SIGN_NEGATIVE))
    {
      result->sign = SIGN_NEGATIVE;
    }
    else
    {
      result->sign = SIGN_POSITIVE;
    }
  }
  else if (operation == OPERATION_OR)
  {
    if ((firstArg->sign == SIGN_POSITIVE) && (secondArg->sign == SIGN_POSITIVE))
    {
      result->sign = SIGN_POSITIVE;
    }
    else
    {
      result->sign = SIGN_NEGATIVE;
    }
  }
  else     // XOR operation
  {
    if (firstArg->sign == secondArg->sign)
    {
      result->sign = SIGN_POSITIVE;
    }
    else
    {
      result->sign = SIGN_NEGATIVE;
    }
  }
  result->nbrLimbs = firstArg->nbrLimbs;
  ConvertToTwosComplement(result);
}

void BigIntAnd(const BigInteger* firstArg,
  const BigInteger* secondArg, BigInteger* result)
{
  InternalBigIntLogical(firstArg, secondArg, result, OPERATION_AND);
}

void BigIntOr(const BigInteger* firstArg,
  const BigInteger* secondArg, BigInteger* result)
{
  InternalBigIntLogical(firstArg, secondArg, result, OPERATION_OR);
}

void BigIntXor(const BigInteger* firstArg,
  const BigInteger* secondArg, BigInteger* result)
{
  InternalBigIntLogical(firstArg, secondArg, result, OPERATION_XOR);
}

void ConvertToTwosComplement(BigInteger *value)
{
  int idx;
  int nbrLimbs;
  limb *ptrLimb;
  if (value->sign == SIGN_POSITIVE)
  {    // If number is positive, no conversion is needed.
    return;
  }
  nbrLimbs = value->nbrLimbs;
  ptrLimb = &value->limbs[0];
  for (idx = 0; idx < nbrLimbs; idx++)
  {
    if (ptrLimb->x != 0)
    {
      break;
    }
    ptrLimb++;
  }
  if (idx < nbrLimbs)
  {
    ptrLimb->x = 0x80000000 - ptrLimb->x;
    ptrLimb++;
  }
  for (; idx < nbrLimbs; idx++)
  {
    ptrLimb->x = 0x7FFFFFFF - ptrLimb->x;
    ptrLimb++;
  }
}


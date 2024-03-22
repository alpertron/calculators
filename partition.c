//
// This file is part of Alpertron Calculators.
//
// Copyright 2016-2021 Dario Alejandro Alpern
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
#include <stdlib.h>
#include <math.h>
#include "bignbr.h"
#include "expression.h"
#ifdef FACTORIZATION_APP
#include "factor.h"
#endif

#define COMPUTE_NEXT_PRIME_SIEVE_SIZE 2000
#define SQRT_MAX_VALUE_LIMB 46341
static limb partArray[100000 + 1000];
static limb prodModulus[MAX_LEN];
static int prodModulusLimbs;
static BigInteger prod;
static BigInteger theFactor;

static void smallMultiply(int factor1, int factor2, int *product)
{
  int low = UintToInt(((unsigned int)factor1 * (unsigned int)factor2) & MAX_VALUE_LIMB);
  double dAccum = (double)factor1 * (double)factor2;
  *product = low;
  if (low < HALF_INT_RANGE)
  {
    dAccum = ((dAccum + (double)FOURTH_INT_RANGE) / (double)LIMB_RANGE);
  }
  else
  {
    dAccum = ((dAccum - (double)FOURTH_INT_RANGE) / (double)LIMB_RANGE);
  }
  *(product + 1) = (unsigned int)dAccum;
}

// Compute the partitions of an integer p(n)
// The formula p(n) = p(k) = p(k - 1) + p(k - 2) - p(k - 5) -
// - p(k - 7) + p(k - 12) + p(k - 15) - p(k - 22) - ...
// where p(0) = 1, where the numbers have the form n(3n-1)/2 for
// n = 1, -1, 2, -2, 3, -3, ... stopping when k-n(3n-1)/2 is negative.
// The signs in the summation are: +, +, -, -, +, +, -, -, ...
// Modular arithmetic is used modulo the largest primes under 2^31, and
// the result is reconstructed.
// The approximation p(n) = e^(pi*sqrt(2n/3))/(4n*sqrt(3)) is used
// in order to find the number of limbs needed to compute the partition.
void partition(int val, BigInteger *pResult)
{
  int index;
  int currentPrime;
  int Q;
  int k;
  int n;
  int sum;
  int idx;
  // Compute approximate number of limbs: log(p(n))/log(2^31)
  // pi * sqrt(2/3)/log(2^31) < 0.12, so 0.12 is selected.
  int limbs = (int)(0.12*sqrt(val) + 1);

  // Compute the primes which will be used for the modular arithmetic
  // operations. The primes must be ordered in ascending order.
  currentPrime = UintToInt(MAX_VALUE_LIMB + 2U);  // Greatest number representable + 2.
  for (index = limbs + val - 1; index >= val; index--)
  {
    for (;;)
    {         // Loop that computes the previous prime number.
      do
      {
        currentPrime -= 2;
      } while ((currentPrime % 3) == 0);
      for (Q = 5; Q <= SQRT_MAX_VALUE_LIMB; Q += 6)
      { /* Check if Base is prime */
        if ((currentPrime % Q) == 0)
        {
          break;     /* Composite */
        }
        if (currentPrime % (Q + 2) == 0)
        {
          break;    /* Composite */
        }
      }
      if (Q > SQRT_MAX_VALUE_LIMB)
      {
        break; /* Prime found */
      }
    }
    partArray[index].x = currentPrime;
  }
  // Perform modular arithmetic.
  for (index = val; index < (val + limbs); index++)
  {
    currentPrime = partArray[index].x;
    sum = 1;                          // Initialize p(0) mod currentPrime.
    for (k = 1; k <= val; k++)        // Generate all partition numbers
    {                                 // up to the one wanted.
      idx = k;
      partArray[k - 1].x = sum;         // Store p(k-1) mod currentPrime.
      sum = 0;
      n = 1;
      for (;;)                        // Loop for n.
      {
        idx -= n + n - 1;
        if (idx < 0)
        {
          break;                      // Negative index, so go out.
        }
        sum -= currentPrime - partArray[idx].x;
        sum += currentPrime & (sum >> 31);
        idx -= n;
        if (idx < 0)
        {
          break;                      // Negative index, so go out.
        }
        sum -= currentPrime - partArray[idx].x;
        sum += currentPrime & (sum >> BITS_PER_GROUP);
        n++;
        idx -= n + n - 1;
        if (idx < 0)
        {
          break;                      // Negative index, so go out.
        }
        sum -= partArray[idx].x;
        sum += currentPrime & (sum >> 31);
        idx -= n;
        if (idx < 0)
        {
          break;                      // Negative index, so go out.
        }
        sum -= partArray[idx].x;
        sum += currentPrime & (sum >> 31);
        n++;
      }
    }
    partArray[index + limbs].x = sum;
  }
  // Reconstruct the result from p(val) mod all primes.
  // v_1 <- u_1 mod m_1
  // v_2 <- (u_2 - v_1) (m_1)^(-1) mod m_2
  // v_3 <- (u_3 - (v_1+m_1 v_2)) (m_1 m_2)^(-1) mod m_3
  // ...
  // v_r <- (u_n - (v_1+m_1(v_2+m_2(v_3+...+m_(r-2) v_(r-1)...)))*
  //        (m_1 m_2 ... m_(r-1))^(-1) mod m_r
  //
  // u = v_r m_(r-1) ..m_2 m_1 + v_3 m_2 m_1 + v_2 m_1 + v_1
  partArray[0].x = partArray[val + limbs].x;
  for (index = 1; index<limbs; index++)
  {
    limb numerator;
    limb prodmod;
    currentPrime = partArray[val + index].x;
    prodmod.x = 1;
    for (k = index - 1; k >= 0; k--)
    {
      smallmodmult(prodmod.x, partArray[val + k].x, &prodmod,
        currentPrime);
    }
    prodmod.x = modInv(prodmod.x, currentPrime);
    numerator.x = partArray[index - 1].x;
    for (k = index - 2; k >= 0; k--)
    {
      smallmodmult(numerator.x, partArray[val + k].x, &numerator,
        currentPrime);
      numerator.x -= currentPrime - partArray[k].x;
      numerator.x += currentPrime & (numerator.x >> 31);
    }
    sum = partArray[val + limbs + index].x - numerator.x;
    if (sum < 0)
    {
      sum += currentPrime;
    }
    smallmodmult(sum, prodmod.x, &partArray[index], currentPrime);
  }
  // Use Chinese Remainder Theorem to find the partition number from
  // the partition number mod different primes.
  pResult->limbs[0].x = partArray[0].x;
  pResult->nbrLimbs = 1;
  pResult->sign = SIGN_POSITIVE;
  prodModulus[0].x = 1;
  prodModulusLimbs = 1;
  for (index = 1; index<limbs; index++)
  {
    int mult;
    int product[2];
    unsigned int carry1 = 0;
    unsigned int carry2 = 0;
    // Update product of modulus by multiplying by next prime.
    mult = partArray[val + index - 1].x;
    for (idx = 0; idx < prodModulusLimbs; idx++)
    {
      smallMultiply(mult, prodModulus[idx].x, product);
      carry1 += (unsigned int)product[0];
      prodModulus[idx].x = UintToInt(carry1 & MAX_VALUE_LIMB);
      carry1 = (carry1 >> BITS_PER_GROUP) + (unsigned int)product[1];
    }
    if (carry1 != 0U)
    {  // New limb needed.
      prodModulus[idx].x = carry1;
      pResult->limbs[idx].x = 0;
      prodModulusLimbs++;
    }
    // Update result.
    carry1 = 0;
    mult = partArray[index].x;
    for (idx = 0; idx < prodModulusLimbs; idx++)
    {
      smallMultiply(mult, prodModulus[idx].x, product);
      carry1 += (unsigned int)product[0];
      carry2 += (carry1 & MAX_VALUE_LIMB) + (unsigned int)pResult->limbs[idx].x;
      carry1 = (carry1 >> BITS_PER_GROUP) + (unsigned int)product[1];
      pResult->limbs[idx].x = UintToInt(carry2 & MAX_VALUE_LIMB);
      carry2 = (carry2 >> BITS_PER_GROUP);
    }
    if ((carry1 + carry2) != 0U)
    {  // New limb needed.
      prodModulus[idx].x = 0;
      pResult->limbs[idx].x = carry1 + carry2;
      prodModulusLimbs++;
    }
  }
  while (pResult->limbs[prodModulusLimbs - 1].x == 0)
  {
    prodModulusLimbs--;
  }
  pResult->nbrLimbs = prodModulusLimbs;
  return;
}

static int numberofBitsSetToOne(int value)
{
  int bitsSet = 0;
  int shiftedValue = value;
  while (shiftedValue > 0)
  {
    if ((shiftedValue & 1) != 0)
    {
      bitsSet++;
    }
    shiftedValue >>= 1;
  }
  return bitsSet;
}

static enum eExprErr ProcessFactorsFactorial(double factorAccum, 
  int nbrGroupsAccumulated, BigInteger *result)
{
  enum eExprErr rc;
  int index;
  int offset;
  prod.limbs[1].x = (int)(factorAccum / (double)LIMB_RANGE);
  prod.limbs[0].x = (int)(factorAccum - ((double)LIMB_RANGE * (double)prod.limbs[1].x));
  prod.nbrLimbs = ((prod.limbs[1].x == 0) ? 1 : 2);
  prod.sign = SIGN_POSITIVE;
  if (((nbrGroupsAccumulated & 1) == 0) || (result != NULL))
  {     // Even means that k multiplications have to be done, where k is the number of 
        // bits set to zero at the right.
    int nbrGroupsAccum = nbrGroupsAccumulated;
    index = numberofBitsSetToOne(nbrGroupsAccumulated - 1);
    while ((nbrGroupsAccum & 1) == 0)
    {
      index--;
      IntArray2BigInteger(&partArray[partArray[index].x].x, &theFactor);
      rc = BigIntMultiply(&prod, &theFactor, &prod);
      if (rc != EXPR_OK)
      {
        return rc;
      }
      nbrGroupsAccum >>= 1;
    }
    offset = partArray[index].x;
    if (result == NULL)
    {
      NumberLength = prod.nbrLimbs;
      BigInteger2IntArray(&partArray[offset].x, &prod);
    }
    else
    {
      CopyBigInt(result, &prod);
    }
  }
  else
  {     // Odd means that the Big Integer has to be stored into the buffer.
    index = numberofBitsSetToOne(nbrGroupsAccumulated) - 1;
    offset = partArray[index].x;
    NumberLength = prod.nbrLimbs;
    BigInteger2IntArray((int *)&partArray[offset], &prod);
  }
  partArray[index + 1].x = offset + partArray[offset].x + 1;
  return EXPR_OK;
}

// Using partArray as a buffer, accumulate products up to two limbs, then multiply the groups recursively.
enum eExprErr factorial(BigInteger *result, int argument, int multifact)
{
  enum eExprErr rc;
  unsigned int shLeft;
  int power2 = 0;
  int nbrGroupsAccumulated = 1;
  double factorAccum = 1;
  double maxFactorAccum = (double)(1U << 30) * (double)(1U << 23);
  partArray[0].x = 20;     // Index of first big integer.
  if (argument == 0)
  {
    intToBigInteger(result, 1);
    return EXPR_OK;
  }
  for (int ctr = argument; ctr > 0; ctr -= multifact)
  {
    unsigned int multiplier = (unsigned int)ctr;
    while ((multiplier & 1U) == 0U)
    {
      multiplier >>= 1;
      power2++;
    }
    if ((factorAccum * (double)multiplier) > maxFactorAccum)
    {
      rc = ProcessFactorsFactorial(factorAccum, nbrGroupsAccumulated, NULL);
      nbrGroupsAccumulated++;
      if (rc != EXPR_OK)
      {
        return rc;
      }
      factorAccum = multiplier;
    }
    else
    {
      factorAccum *= multiplier;
    }
  }
  shLeft = numberofBitsSetToOne(nbrGroupsAccumulated - 1);
  nbrGroupsAccumulated = UintToInt(1U << shLeft);
  rc = ProcessFactorsFactorial(factorAccum, nbrGroupsAccumulated, result);
  if (rc != EXPR_OK)
  {
    return rc;
  }
  return BigIntMultiplyPower2(result, power2);
}

enum eExprErr primorial(BigInteger *result, int argument)
{
  enum eExprErr rc;
  int j;
  int nbrGroupsAccumulated = 1;
  double factorAccum = 1;
  unsigned int shLeft;
  unsigned int firstFactor = 1U << 30;
  unsigned int secondFactor = 1U << 23;
  double maxFactorAccum = (double)firstFactor * (double)secondFactor;
  partArray[0].x = 20;     // Index of first big integer.
  for (int ctr = 2; ctr <= argument; ctr++)
  {
    for (j = 2; (j*j) <= ctr; j++)
    {
      if (((ctr / j) * j) == ctr)
      {   // Number is not prime.
        break;
      }
    }
    if ((j*j) > ctr)
    {     // Number is prime, perform multiplication.
      if ((factorAccum * (double)ctr) > maxFactorAccum)
      {
        rc = ProcessFactorsFactorial(factorAccum, nbrGroupsAccumulated, NULL);
        if (rc != EXPR_OK)
        {
          return rc;
        }
        nbrGroupsAccumulated++;
        factorAccum = 1;
      }
      factorAccum *= ctr;
    }
  }
  shLeft = numberofBitsSetToOne(nbrGroupsAccumulated - 1);
  nbrGroupsAccumulated = UintToInt(1U << shLeft);
  return ProcessFactorsFactorial(factorAccum, nbrGroupsAccumulated, result);
}

enum eExprErr ComputeFibLucas(int origValue, BigInteger* pArgument)
{
  limb largeVal;
  int val;
  int len;
  int lenBytes;
  limb* pFibonPrev;
  limb* pFibonAct;
  if (pArgument->sign == SIGN_NEGATIVE)
  {
    return EXPR_INVALID_PARAM;
  }
  if (pArgument->nbrLimbs > 1)
  {
    return EXPR_INTERM_TOO_HIGH;
  }
  largeVal.x = pArgument->limbs[0].x;
#ifdef FACTORIZATION_APP
  if (largeVal.x > 956620)
#else
  if (largeVal.x > 95662)
#endif
  {
    return EXPR_INTERM_TOO_HIGH;
  }
  val = largeVal.x;
  pFibonPrev = prod.limbs;
  pFibonAct = theFactor.limbs;
  len = 1;
  if (val == 0)
  {        // F(0) = 0, L(0) = 2
    pFibonAct->x = origValue;
  }
  else
  {
    int j;
    // For Lucas sequences: FibonPrev = 2, FibonAct = 1
    // For Fibonacci sequences: FibonPrev = 0, FibonAct = 1
    pFibonPrev->x = origValue;
    pFibonAct->x = 1;
    for (int i = 1; i < val; i++)
    {
      limb* pTemp;
      unsigned int carry = 0;
      for (j = 0; j < len; j++)
      {
        carry += (unsigned int)(pFibonPrev + j)->x + (unsigned int)(pFibonAct + j)->x;
        (pFibonPrev + j)->x = UintToInt(carry & MAX_VALUE_LIMB);
        carry >>= BITS_PER_GROUP;
      }
      if (carry != 0U)
      {
        (pFibonPrev + j)->x = carry;
        (pFibonAct + j)->x = 0;
        len++;
      }
      pTemp = pFibonAct;
      pFibonAct = pFibonPrev;
      pFibonPrev = pTemp;
    }
  }
  pArgument->sign = SIGN_POSITIVE;
  pArgument->nbrLimbs = len;
  lenBytes = len * (int)sizeof(limb);
  (void)memcpy(pArgument->limbs, pFibonAct, lenBytes);
  return EXPR_OK;
}

static void generateSieve(const int* pSmallPrimes, char* sieve,
  const BigInteger* pArgument, bool isNext)
{
  const int* ptrSmallPrimes = pSmallPrimes;
  int maxPrimeIndex = 1229;
  if (pArgument->nbrLimbs <= 2)
  {
    double value;
    if (pArgument->nbrLimbs == 1)
    {
      value = (double)pArgument->limbs[0].x;
    }
    else
    {
      value = (double)pArgument->limbs[0].x +
        (double)pArgument->limbs[1].x * (double)LIMB_RANGE;
    }
    maxPrimeIndex = (int)sqrt(sqrt(value));
    if (maxPrimeIndex > 1229)
    {
      maxPrimeIndex = 1229;
    }
  }
  // Indicate numbers not divisible by small primes in advance.
  (void)memset(sieve, 0, COMPUTE_NEXT_PRIME_SIEVE_SIZE);
  for (int ctr = 0; ctr < maxPrimeIndex; ctr++)
  {     // For each prime less than 10000...
    int prime = *ptrSmallPrimes;
    ptrSmallPrimes++;
    int remainder = getRemainder(pArgument, prime);
    // Compute first element of sieve to indicate multiple of prime.
    if (isNext)
    {
      if (remainder > 0)
      {
        remainder = prime - remainder;
      }
    }
    else
    {
      remainder = (COMPUTE_NEXT_PRIME_SIEVE_SIZE - remainder) % prime;
    }
    if (remainder < 0)
    {
      remainder += prime;
    }
    for (; remainder < COMPUTE_NEXT_PRIME_SIEVE_SIZE; remainder += prime)
    {        // Indicate numbers are divisible by this prime.
      sieve[remainder] = 1;
    }
  }
}

// Compute previous probable prime
enum eExprErr ComputeBack(BigInteger* pArgument)
{
  char sieve[COMPUTE_NEXT_PRIME_SIEVE_SIZE];
  BigInteger* pResult = pArgument;
  limb* pResultLimbs = pResult->limbs;
  const limb* pArgumentLimbs = pArgument->limbs;
  pResult->sign = SIGN_POSITIVE;
  if (pArgument->sign == SIGN_NEGATIVE)
  {
    return EXPR_INVALID_PARAM;
  }
  if (pArgument->nbrLimbs == 1)
  {
    if (pArgumentLimbs->x < 3)
    {
      return EXPR_INVALID_PARAM;
    }
    if (pArgumentLimbs->x == 3)
    {
      pResult->nbrLimbs = 1;
      pResultLimbs->x = 2;
      return EXPR_OK;
    }
  }
  initializeSmallPrimes(smallPrimes);
  CopyBigInt(pResult, pArgument);
  pResultLimbs->x |= 1;  // If number is even, use next odd number.
  if (pResult->nbrLimbs == 1)
  {
    int rc;
    do
    {        // Loop that searches for previous probable prime.
      addbigint(pResult, -2);
#ifdef FACTORIZATION_APP
      rc = BpswPrimalityTest(pResult, NULL);  // Continue loop if not probable prime.
#else
      rc = BpswPrimalityTest(pResult);        // Continue loop if not probable prime.
#endif
    } while (rc != 0);
  }
  else
  {          // Big number: use sieve.
    for (;;)
    {        // Loop that searches for previous probable prime.
      addbigint(pResult, -COMPUTE_NEXT_PRIME_SIEVE_SIZE);
      generateSieve(smallPrimes, sieve, pResult, true);
      for (int ctr = COMPUTE_NEXT_PRIME_SIEVE_SIZE - 1; ctr >= 0; ctr--)
      {
        if (sieve[ctr] == 0)
        {   // Number is not divisible by primes less than 1000.
          addbigint(pResult, ctr);
#ifdef FACTORIZATION_APP
          if (BpswPrimalityTest(pResult, NULL) == 0)  // End loop if probable prime.
#else
          if (BpswPrimalityTest(pResult) == 0)        // End loop if probable prime.
#endif
          {
            return EXPR_OK;
          }
          addbigint(pResult, -ctr);
        }
      }
    }
  }
  return EXPR_OK;
}

// Compute next probable prime
enum eExprErr ComputeNext(BigInteger* pArgument)
{
  char sieve[COMPUTE_NEXT_PRIME_SIEVE_SIZE];
  BigInteger* pResult = pArgument;
  limb* pResultLimbs = pResult->limbs;
  const limb* pArgumentLimbs = pArgument->limbs;
  if ((pArgument->sign == SIGN_NEGATIVE) ||
    ((pArgument->nbrLimbs == 1) && (pArgumentLimbs->x < 2)))
  {    // If number is less than 2, the next prime is 2.
    pResult->sign = SIGN_POSITIVE;
    pResult->nbrLimbs = 1;
    pResultLimbs->x = 2;
    return EXPR_OK;
  }
  initializeSmallPrimes(smallPrimes);
  CopyBigInt(pResult, pArgument);
  if ((pResult->limbs[0].x & 1) == 0)
  {   // Number is even.
    addbigint(pResult, 1);
  }
  else
  {   // Number is odd.
    addbigint(pResult, 2);
  }
  if (pResult->nbrLimbs == 1)
  {
    for (;;)
    {        // Loop that searches for next probable prime.
#ifdef FACTORIZATION_APP
      if (BpswPrimalityTest(pResult, NULL) == 0)  // Continue loop if not probable prime.
#else
      if (BpswPrimalityTest(pResult) == 0)        // Continue loop if not probable prime.
#endif
      {
        return EXPR_OK;
      }
      addbigint(pResult, 2);
    }
  }
  else
  {          // Big number: use sieve.
    for (;;)
    {        // Loop that searches for next probable prime.
      generateSieve(smallPrimes, sieve, pResult, true);
      for (int ctr = 0; ctr < COMPUTE_NEXT_PRIME_SIEVE_SIZE; ctr++)
      {
        if (sieve[ctr] == 0)
        {   // Number is not divisible by primes less than 1000.
          addbigint(pResult, ctr);
#ifdef FACTORIZATION_APP
          if (BpswPrimalityTest(pResult, NULL) == 0)  // End loop if probable prime.
#else
          if (BpswPrimalityTest(pResult) == 0)        // End loop if probable prime.
#endif
          {
            return EXPR_OK;
          }
          addbigint(pResult, -ctr);
        }
      }
      addbigint(pResult, COMPUTE_NEXT_PRIME_SIEVE_SIZE);
    }
  }
}


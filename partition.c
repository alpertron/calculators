/*
This file is part of Alpertron Calculators.

Copyright 2016 Dario Alejandro Alpern

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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bignbr.h"
#include "expression.h"

#define SQRT_MAX_VALUE_LIMB 46341
static int partArray[100000 + 1000];
static limb prodModulus[MAX_LEN];
static int prodModulusLimbs;

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
  int index, currentPrime, Q, k, n, sum, idx;
  limb carry;
  // Compute approximate number of limbs: log(p(n))/log(2^31)
  // pi * sqrt(2/3)/log(2^31) < 0.12, so 0.12 is selected.
  int limbs = (int)(0.12*sqrt(val) + 1);

  // Compute the primes which will be used for the modular arithmetic
  // operations. The primes must be ordered in ascending order.
  currentPrime = MAX_VALUE_LIMB + 2;    // Greatest number representable + 2.
  for (index = limbs + val - 1; index >= val; index--)
  {
    for (;;)
    {         // Loop that computes the previous prime number.
      do
      {
        currentPrime -= 2;
      } while (currentPrime % 3 == 0);
      for (Q = 5; Q <= SQRT_MAX_VALUE_LIMB; Q += 6)
      { /* Check if Base is prime */
        if (currentPrime % Q == 0)
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
    partArray[index] = currentPrime;
  }
  // Perform modular arithmetic.
  for (index = val; index<val + limbs; index++)
  {
    currentPrime = partArray[index];
    sum = 1;                          // Initialize p(0) mod currentPrime.
    for (k = 1; k <= val; k++)        // Generate all partition numbers
    {                                 // up to the one wanted.
      idx = k;
      partArray[k - 1] = sum;         // Store p(k-1) mod currentPrime.
      sum = 0;
      n = 1;
      for (;;)                        // Loop for n.
      {
        idx -= n + n - 1;
        if (idx < 0)
        {
          break;                      // Negative index, so go out.
        }
        sum -= currentPrime - partArray[idx];
        sum += currentPrime & (sum >> 31);
        idx -= n;
        if (idx < 0)
        {
          break;                      // Negative index, so go out.
        }
        sum -= currentPrime - partArray[idx];
        sum += currentPrime & (sum >> 31);
        n++;
        idx -= n + n - 1;
        if (idx < 0)
        {
          break;                      // Negative index, so go out.
        }
        sum -= partArray[idx];
        sum += currentPrime & (sum >> 31);
        idx -= n;
        if (idx < 0)
        {
          break;                      // Negative index, so go out.
        }
        sum -= partArray[idx];
        sum += currentPrime & (sum >> 31);
        n++;
      }
    }
    partArray[index + limbs] = sum;
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
  partArray[0] = partArray[val + limbs];
  for (index = 1; index<limbs; index++)
  {
    long numerator, prodmod;
    currentPrime = partArray[val + index];
    prodmod = 1;
    for (k = index - 1; k >= 0; k--)
    {
      prodmod = prodmod * partArray[val + k] % currentPrime;
    }
    prodmod = modInv((int)prodmod, currentPrime);
    numerator = partArray[index - 1];
    for (k = index - 2; k >= 0; k--)
    {
      numerator = (numerator*partArray[val + k] + partArray[k]) %
        currentPrime;
    }
    sum = partArray[val + limbs + index] - (int)numerator;
    if (sum<0)
    {
      sum += currentPrime;
    }
    partArray[index] = (int)(sum*prodmod%currentPrime);
  }
  // Use Chinese Remainder Theorem to find the partition number from
  // the partition number mod different primes.
  pResult->limbs[0].x = partArray[0];
  pResult->nbrLimbs = 1;
  pResult->sign = SIGN_POSITIVE;
  prodModulus[0].x = 1;
  prodModulusLimbs = 1;
  for (index = 1; index<limbs; index++)
  {
    int mult;
    // Update product of modulus by multiplying by next prime.
    carry.x = 0;
    mult = partArray[val + index - 1];
    for (idx = 0; idx < prodModulusLimbs; idx++)
    {
      carry.x += mult*prodModulus[idx].x;
      prodModulus[idx].x = carry.x & MAX_VALUE_LIMB;
      carry.x >>= BITS_PER_GROUP;
    }
    if (carry.x != 0)
    {  // New limb needed.
      prodModulus[idx].x = carry.x;
      pResult->limbs[idx].x = 0;
      prodModulusLimbs++;
    }
    // Update result.
    carry.x = 0;
    mult = partArray[index];
    for (idx = 0; idx < prodModulusLimbs; idx++)
    {
      carry.x += mult*prodModulus[idx].x + pResult->limbs[idx].x;
      pResult->limbs[idx].x = carry.x & MAX_VALUE_LIMB;
      carry.x >>= BITS_PER_GROUP;
    }
    if (carry.x != 0)
    {  // New limb needed.
      prodModulus[idx].x = 0;
      pResult->limbs[idx].x = carry.x;
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

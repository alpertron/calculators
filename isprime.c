//
// This file is part of Alpertron Calculators.
//
// Copyright 2021 Dario Alejandro Alpern
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
#include <stdbool.h>
#include "isprime.h"

#define SMALL_NUMBER_BOUND 32768

int TestNbr[NBR_LIMBS];
unsigned int MontgomeryMultN;
int MontgomeryMultR1[NBR_LIMBS+1];  // One more limb required for AdjustModN.
int MontgomeryMultR2[NBR_LIMBS+1];

static char primes[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41,
                         43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97 };
int getPrime(int index)
{
  return primes[index];
}

// Each limb has 31 bits, so with 2 limbs these routines can use numbers up to 2^62
// This is over 10^18.

// Compute Nbr <- Nbr mod TestNbr.
static void AdjustModN(int *Nbr)
{
  int i;
  int carry;
  int TrialQuotient;
  double dNbr;
  double dModulus;
  double dTrialQuotient;
  double dDelta;
  double dVal = 1.0 / (double)LIMB_RANGE;
  double dSquareLimb = (double)LIMB_RANGE * (double)LIMB_RANGE;

  dModulus = (double)TestNbr[1] + (double)TestNbr[0] * dVal;
  dNbr = (double)Nbr[2] * (double)LIMB_RANGE + (double)Nbr[1] + (double)Nbr[0] * dVal;
  TrialQuotient = (int)(unsigned int)floor((dNbr / dModulus) + 0.5);
  if ((unsigned int)TrialQuotient >= LIMB_RANGE)
  {   // Maximum value for limb.
    TrialQuotient = MAX_VALUE_LIMB;
  }
  // Compute Nbr <- Nbr - TrialQuotient * Modulus
  dTrialQuotient = (double)TrialQuotient;
  carry = 0;
  dDelta = 0;
  for (i = 0; i < NBR_LIMBS; i++)
  {
    double dAccumulator;
    int low = (*(Nbr + i) - (TestNbr[i] * TrialQuotient) + carry) & MAX_INT_NBR;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    dAccumulator = *(Nbr+i) - ((double)TestNbr[i] * dTrialQuotient) + (double)carry + dDelta;
    dDelta = 0.0;
    if (dAccumulator < 0.0)
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
    *(Nbr + i) = low;
  }
  *(Nbr + NBR_LIMBS) = (*(Nbr + NBR_LIMBS) + carry) & MAX_INT_NBR;
  if (((unsigned int)*(Nbr+NBR_LIMBS) & MAX_VALUE_LIMB) != 0U)
  {
    unsigned int cy = 0;
    for (i = 0; i < NBR_LIMBS; i++)
    {
      unsigned int tmp;
      cy += (unsigned int)*(Nbr+i) + (unsigned int)TestNbr[i];
      tmp = cy & MAX_VALUE_LIMB;
      *(Nbr+i) = (int)tmp;
      cy >>= BITS_PER_GROUP;
    }
    *(Nbr+NBR_LIMBS) = 0;
  }
}
static void smallmodmult(int factor1, int factor2, int *product, int mod)
{
  if (mod < SMALL_NUMBER_BOUND)
  {
    *product = factor1 * factor2 % mod;
  }
  else
  {   // TestNbr has one limb but it is not small.
#ifdef _USING64BITS_
    *product = (int64_t)factor1 * factor2 % mod;
#else
      // Round up quotient.
    int quotient = (int)floor((double)factor1 * (double)factor2 / (double)mod + 0.5);
    int remainder = (factor1 * factor2) - (quotient * mod);
    if (remainder < 0)
    {    // Quotient was 1 more than expected. Adjust remainder.
      remainder += mod;
    }
    *product = remainder;
#endif
  }
}

static void MontgomeryMult(int *factor1, int *factor2, int *Product)
{
#ifndef _USING64BITS_
  unsigned int carry;
#endif
  if (TestNbr[1] == 0)
  {
    smallmodmult(*factor1, *factor2, Product, TestNbr[0]);
    *(Product+1) = 0;
    return;
  }
  int TestNbr0 = TestNbr[0];
  int TestNbr1 = TestNbr[1];
  uint32_t Prod0;
  uint32_t Prod1;
  int factor2_0 = *factor2;
  int factor2_1 = *(factor2+1);
#ifdef _USING64BITS_
  uint64_t Pr;
  unsigned int Nbr;
  unsigned int MontDig;
  
  Nbr = *factor1;
  Pr = Nbr * (uint64_t)factor2_0;
  MontDig = ((uint32_t)Pr * MontgomeryMultN) & MAX_VALUE_LIMB;
  Pr = ((((uint64_t)MontDig * (uint64_t)TestNbr0) + Pr) >> BITS_PER_GROUP) +
    ((uint64_t)MontDig * (uint64_t)TestNbr1) + ((uint64_t)Nbr * (uint64_t)factor2_1);
  Prod0 = (uint32_t)(Pr & MAX_VALUE_LIMB);
  Prod1 = (uint32_t)(Pr >> BITS_PER_GROUP);
    
  Nbr = *(factor1 + 1);
  Pr = Nbr * (uint64_t)factor2_0 + (uint32_t)Prod0;
  MontDig = ((uint32_t)Pr * MontgomeryMultN) & MAX_VALUE_LIMB;
  Pr = ((((uint64_t)MontDig * (uint64_t)TestNbr0) + Pr) >> BITS_PER_GROUP) +
    ((uint64_t)MontDig * (uint64_t)TestNbr1) + ((uint64_t)Nbr * (uint64_t)factor2_1) + (uint32_t)Prod1;
  Prod0 = (uint32_t)(Pr & MAX_VALUE_LIMB);
  Prod1 = (uint32_t)(Pr >> BITS_PER_GROUP);
    
  if (((Pr >> BITS_PER_GROUP) > (uint64_t)TestNbr1) ||
     ((Prod1 == (uint32_t)TestNbr1) && (Prod0 >= (uint32_t)TestNbr0)))
  {
    uint32_t borrow = Prod0 - (uint32_t)TestNbr0;
    Prod0 = borrow & MAX_VALUE_LIMB;
    Prod1 = (Prod1 - (uint32_t)TestNbr1 - (borrow >> BITS_PER_GROUP)) & MAX_VALUE_LIMB;
  }
#else
  double dInvLimbRange = 1.0 / (double)LIMB_RANGE;
  int Nbr = *(factor1);
  double dNbr = (double)Nbr;
  int low = Nbr * factor2_0;
  double dAccum = dNbr * (double)factor2_0;
  unsigned int tmp = ((unsigned int)low * MontgomeryMultN) & MAX_VALUE_LIMB;
  int MontDig = (int)tmp;
  double dMontDig = (double)MontDig;
  dAccum += dMontDig * (double)TestNbr0;
  // At this moment dAccum is multiple of LIMB_RANGE.
  dAccum = floor((dAccum*dInvLimbRange) + 0.5);
  tmp = ((unsigned int)dAccum + ((unsigned int)MontDig * (unsigned int)TestNbr1) +
               ((unsigned int)Nbr * (unsigned int)factor2_1)) & MAX_VALUE_LIMB;
  low = (int)tmp;
  dAccum += (dMontDig * (double)TestNbr1) + (dNbr * (double)factor2_1);
  Prod0 = low;
  if (low < HALF_INT_RANGE)
  {
    dAccum = ((dAccum + (double)FOURTH_INT_RANGE)*dInvLimbRange);
  }
  else
  {
    dAccum = ((dAccum - (double)FOURTH_INT_RANGE)*dInvLimbRange);
  }
  Prod1 = (unsigned int)dAccum;  // Most significant limb can be greater than LIMB_RANGE
  
  Nbr = *(factor1 + 1);
  dNbr = (double)Nbr;
  low = (Nbr * factor2_0) + (int)Prod0;
  dAccum = (dNbr * (double)factor2_0) + (double)Prod0;
  tmp = ((unsigned int)low * MontgomeryMultN) & MAX_VALUE_LIMB;
  MontDig = (int)tmp;
  dMontDig = (double)MontDig;
  dAccum += dMontDig * (double)TestNbr0;
  // At this moment dAccum is multiple of LIMB_RANGE.
  dAccum = floor((dAccum*dInvLimbRange) + 0.5);
  tmp = ((unsigned int)dAccum + (MontDig * TestNbr1) +
               (Nbr * factor2_1) + Prod1) & MAX_VALUE_LIMB;
  low = (int)tmp;
  dAccum += (dMontDig * (double)TestNbr1) + (dNbr * (double)factor2_1) + (double)Prod1;
  Prod0 = low;
  if (low < HALF_INT_RANGE)
  {
    dAccum = ((dAccum + (double)FOURTH_INT_RANGE)*dInvLimbRange);
  }
  else
  {
    dAccum = ((dAccum - (double)FOURTH_INT_RANGE)*dInvLimbRange);
  }
  Prod1 = (unsigned int)dAccum;  // Most significant limb can be greater than LIMB_RANGE
  
  if ((Prod1 > (unsigned int)TestNbr1) ||
       ((Prod1 == (unsigned int)TestNbr1) && (Prod0 >= (unsigned int)TestNbr0)))
  {        // Prod >= TestNbr, so perform Prod <- Prod - TestNbr
    carry = Prod0 - (unsigned int)TestNbr0;
    tmp = carry & MAX_VALUE_LIMB;
    Prod0 = (int)tmp;
    // On subtraction, carry subtracts also.
    Prod1 = (Prod1 - TestNbr1 - (carry >> BITS_PER_GROUP)) & MAX_VALUE_LIMB;
  }
#endif  
  *Product = Prod0;
  *(Product+1) = Prod1;
}

void AddBigNbr(const int *Nbr1, const int *Nbr2, int *Sum)
{
  unsigned int carry = (unsigned int)*Nbr1 + (unsigned int)*Nbr2;
  *Sum = carry & MAX_VALUE_LIMB;
  carry = (carry >> BITS_PER_GROUP) + *(Nbr1+1) + *(Nbr2+1);
  *(Sum+1) = carry & MAX_VALUE_LIMB;
}

void SubtBigNbr(const int *Nbr1, const int *Nbr2, int *Diff)
{
  unsigned int borrow = (unsigned int)*Nbr1 - (unsigned int)*Nbr2;
  *Diff = borrow & MAX_VALUE_LIMB;
  // On subtraction, borrow subtracts too.
  borrow = *(Nbr1+1) - *(Nbr2+1) - (borrow >> BITS_PER_GROUP);
  *(Diff+1) = borrow & MAX_VALUE_LIMB;
}

static void AddBigNbrModN(const int *Nbr1, const int *Nbr2, int *Sum)
{
  unsigned int Sum0;
  unsigned int Sum1;
  unsigned int TestNbr0 = (unsigned int)TestNbr[0];
  unsigned int TestNbr1 = (unsigned int)TestNbr[1];
  unsigned int carry = (unsigned int)*Nbr1 + (unsigned int)*Nbr2;
  Sum0 = carry & MAX_VALUE_LIMB;
  carry = (carry >> BITS_PER_GROUP) + *(Nbr1 + 1) + *(Nbr2 + 1);
  Sum1 = carry & MAX_VALUE_LIMB;
  if ((carry > TestNbr1) || 
     ((carry == TestNbr1) && (Sum0 >= TestNbr0)))
  {
    unsigned int borrow = Sum0 - TestNbr0;
    Sum0 = borrow & MAX_VALUE_LIMB;
    // On subtraction, borrow subtracts too.
    Sum1 = (Sum1 - TestNbr1 - (borrow >> BITS_PER_GROUP)) & MAX_INT_NBR;
  }
  *Sum = (int)Sum0;
  *(Sum+1) = (int)Sum1;
}

static void GetMontgomeryParms(void)
{
  int N;
  int x;

  MontgomeryMultR1[1] = 0;
  if (TestNbr[1] == 0)
  {
    MontgomeryMultR1[0] = 1;
    return;
  }
  N = (int)TestNbr[0];   // 2 least significant bits of inverse correct.
  x = N;
  x = x * (2 - (N * x)); // 4 least significant bits of inverse correct.
  x = x * (2 - (N * x)); // 8 least significant bits of inverse correct.
  x = x * (2 - (N * x)); // 16 least significant bits of inverse correct.
  x = x * (2 - (N * x)); // 32 least significant bits of inverse correct.
  MontgomeryMultN = (-x) & MAX_INT_NBR;
  MontgomeryMultR1[2] = 1;
  MontgomeryMultR1[0] = 0;
  AdjustModN(MontgomeryMultR1);
}

// Perform Miller-Rabin test of number stored in variable TestNbr.
// The bases to be used are 2, 3, 5, 7, 11, 13, 17, 19 and 23 which
// ensures that any composite less than 3*10^18 is discarded.
bool isPrime(int *value)
{
  static const char bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 0};
  // List of lowest composite numbers that passes Miller-Rabin for above bases (OEIS A014233).
  // The first number is the limit modulo 2^31, the second number is the limit divided 2^31.
  static const int limits[] =
  {
    0, 0,

    // Base 2: limit = 2047
    2047, 0,

    // Bases 2 and 3: limit = 1373653
    1373653, 0,              

    // Bases 2, 3 and 5: limit = 25326001
    25326001, 0,             

    // Bases 2, 3, 5 and 7: limit = 3215031751
    1067548103, 1,
 
    // Bases 2, 3, 5, 7 and 11: limit = 2152302898747
    524283451, 1002,

    // Bases 2, 3, 5, 7, 11 and 13: limit = 3474749660383
    121117919, 1618,

    // Bases 2, 3, 5, 7, 11, 13 and 17: limit = 341550071728321
    1387448513, 159046,

    // Bases 2, 3, 5, 7, 11, 13, 17 and 19: limit = 341550071728321
    1387448513, 159046,

    // Greater than any argument of isPrime()
    0, 2000000000,
  };
  int i;
  int j;
  int index;
  int indexLSB;
  int indexMSB;
  int idxNbrMSB;
  unsigned int mask;
  unsigned int maskMSB;
  unsigned int prevBase = 1;
  unsigned int base;
  int baseInMontRepres[NBR_LIMBS];
  int power[NBR_LIMBS];
  int temp[NBR_LIMBS];
    // Convert parameter to big number (2 limbs of 31 bits each).
  int TestNbr0 = value[0];
  TestNbr[0] = value[0];
  int TestNbr1 = value[1];
  TestNbr[1] = value[1];
  if (TestNbr1 == 0)
  {
    if (TestNbr0 == 1)
    {
      return false;            // 1 is not prime.
    }
    if (TestNbr0 == 2)
    {
      return true;            // 2 is prime.
    }
  }
  if ((TestNbr0 & 1) == 0)
  {
    return false;              // Even numbers different from 2 are not prime.
  }
  for (i=1; i<sizeof(primes); i++)
  {
    base = primes[i];
    if (TestNbr1 == 0)
    {
      if (TestNbr0 == base)
      {
        return true;          // Number is prime.
      }
      if ((TestNbr0 % base) == 0)
      {
        return false;          // Number is multiple of base, so it is composite.
      }
    }
    // Check whether TestNbr is multiple of base. In this case the number would be composite.
    // No overflow possible in next expression.
    else if ((((TestNbr1 % base) * (LIMB_RANGE % base) + (TestNbr0 % base)) % base) == 0)
    {
      return false;            // Number is multiple of base, so it is composite.
    }
    else
    {                          // Nothing to do.
    }
  }
  GetMontgomeryParms();
  baseInMontRepres[0] = MontgomeryMultR1[0];
  baseInMontRepres[1] = MontgomeryMultR1[1];
  
       // Find index of least significant bit set disregarding bit 0.
  mask = TestNbr0 & (MAX_INT_NBR-1);
  indexLSB = 0;
  if (mask == 0U)
  {    // Least significant bit is inside high limb.
    mask = TestNbr1;
    indexLSB = BITS_PER_GROUP;
  }
  if ((mask & 0xFFFFU) == 0U)
  {    // Least significant bit is somewhere between bits 31-16.
    mask >>= 16;
    indexLSB += 16;
  }
  if ((mask & 0xFFU) == 0U)
  {    // Least significant bit is somewhere between bits 15-8.
    mask >>= 8;
    indexLSB += 8;
  }
  if ((mask & 0x0FU) == 0U)
  {    // Least significant bit is somewhere between bits 7-4.
    mask >>= 4;
    indexLSB += 4;
  }
  if ((mask & 0x03U) == 0U)
  {    // Least significant bit is between bits 3-2.
    mask >>= 2;
    indexLSB += 2;
  }
  if ((mask & 0x01U) == 0U)
  {    // Least significant bit is bit 1.
    indexLSB++;
  }
      // Find index of most significant bit set.
  mask = TestNbr1;
  indexMSB = BITS_PER_GROUP;
  idxNbrMSB = 1;
  if (mask == 0U)
  {    // Most significant bit is inside low limb.
    mask = TestNbr0;
    indexMSB = 0;
    idxNbrMSB = 0;
  }
  if ((mask & 0xFFFF0000U) != 0U)
  {    // Most significant bit is somewhere between bits 31-16.
    mask >>= 16;
    indexMSB += 16;
  }
  if ((mask & 0xFF00U) != 0U)
  {    // Most significant bit is somewhere between bits 15-8.
    mask >>= 8;
    indexMSB += 8;
  }
  if ((mask & 0xF0U) != 0U)
  {    // Most significant bit is somewhere between bits 7-4.
    mask >>= 4;
    indexMSB += 4;
  }
  if ((mask & 0x0CU) != 0U)
  {    // Most significant bit is between bits 3-2.
    mask >>= 2;
    indexMSB += 2;
  }
  if ((mask & 0x02U) != 0U)
  {    // Most significant bit is bit 1.
    indexMSB++;
  }
  maskMSB = (1<<(indexMSB % BITS_PER_GROUP));
  i = 0;
  j = 0;
  while ((limits[j+1] < TestNbr1) || ((limits[j+1] == TestNbr1) && (limits[j] < TestNbr0)))
  {
    int idxNbr;
    base = bases[i];
    do                     // Compute next base in Montgomery representation.
    {
      AddBigNbrModN(baseInMontRepres, MontgomeryMultR1, baseInMontRepres);
      prevBase++;
    } while (prevBase < base);
    i++;
    j += 2;
    power[0] = baseInMontRepres[0];
    power[1] = baseInMontRepres[1];
    mask = maskMSB;
    idxNbr = idxNbrMSB;
    for (index = indexMSB-1; index >= indexLSB; index--)
    {
      mask >>= 1;
      if (mask == 0U)
      {
        mask = HALF_INT_RANGE;
        idxNbr--;
      }
      MontgomeryMult(power, power, power);

      if (((unsigned int)TestNbr[idxNbr] & mask) != 0U)
      {
        MontgomeryMult(power, baseInMontRepres, power);
      }
    }
       // If power equals 1 or -1 in Montgomery representation,
       // another base must be tried.
    if ((power[0] == MontgomeryMultR1[0]) && (power[1] == MontgomeryMultR1[1]))
    {
      continue;   // power equals 1, so another base must be tried.
    }
    AddBigNbrModN(power, MontgomeryMultR1, temp);
    if ((temp[0] == 0) && (temp[1] == 0))
    {
      continue;   // power equals -1, so another base must be tried.
    }
    for (; index>=0; index--)
    {
      mask >>= 1;
      if (mask == 0)
      {
        mask = HALF_INT_RANGE;
        idxNbr--;
      }
      MontgomeryMult(power, power, power);
      if ((power[0] == MontgomeryMultR1[0]) && (power[1] == MontgomeryMultR1[1]))
      {
        return false;  // power equals 1, so number is composite.
      }
      AddBigNbrModN(power, MontgomeryMultR1, temp);
      if ((temp[0] == 0) && (temp[1] == 0))
      {                // power equals -1.
        break;
      }
    }
    if (index >= 0)
    {
      continue;        // power equals -1, so another base must be tried.
    }
    return false;      // base^(n-1) != 1 (mod n) is not 1, so number is composite.
  }
  return true;         // All Miller-Rabin tests were passed, so number is prime.
}

void multiply(int factor1, int factor2, int *prod)
{
#ifdef _USING64BITS_
  *prod = (factor1 * factor2) & MAX_VALUE_LIMB;
  *(prod+1) = (int)((long long)factor1 * (long long)factor2 >> BITS_PER_GROUP);
#else
  int low = (factor1 * factor2) & MAX_VALUE_LIMB;
  double dAccum = (double)factor1 * (double)factor2;
  *prod = low;
  if (low < HALF_INT_RANGE)
  {
    dAccum = ((dAccum + (double)FOURTH_INT_RANGE) / LIMB_RANGE);
  }
  else
  {
    dAccum = ((dAccum - (double)FOURTH_INT_RANGE) / LIMB_RANGE);
  }
  *(prod+1) = (unsigned int)dAccum;
#endif 
}

char* appendInt(char* text, int intValue)
{
  char* ptrText = text;
  int value = intValue;
  int div = 1000000000;
  bool zeroIsSignificant = false;
  if (value < 0)
  {
    value = -value;
    *ptrText = '-';
    ptrText++;
  }
  do
  {
    int quot = value / div;
    if ((quot != 0) || zeroIsSignificant)
    {
      zeroIsSignificant = true;
      *ptrText = (char)quot + '0';
      ptrText++;
      value -= quot * div;
    }
    div /= 10;
  } while (div > 1);
  *ptrText = (char)value + '0';
  ptrText++;
  *ptrText = 0;
  return ptrText;
}

//
// This file is part of Alpertron Calculators.
//
// Copyright 2021-2025 Dario Alejandro Alpern
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

int TestNbrGraphic[NBR_LIMBS];
unsigned int MontMultGraphicN;
int MontMultGraphicR1[NBR_LIMBS+1];  // One more limb required for AdjustModN.

static char primes[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41,
                         43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97 };
int getPrime(int index)
{
  return primes[index];
}

// Each limb has 31 bits, so with 2 limbs these routines can use numbers up to 2^62
// This is over 10^18.

// Compute Nbr <- Nbr mod TestNbrGraphic.
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

  dModulus = (double)TestNbrGraphic[1] + (double)TestNbrGraphic[0] * dVal;
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
    int low = (*(Nbr + i) - (TestNbrGraphic[i] * TrialQuotient) + carry) & MAX_INT_NBR;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    dAccumulator = *(Nbr+i) - ((double)TestNbrGraphic[i] * dTrialQuotient) + (double)carry + dDelta;
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
      cy += (unsigned int)*(Nbr+i) + (unsigned int)TestNbrGraphic[i];
      tmp = cy & MAX_VALUE_LIMB;
      *(Nbr+i) = (int)tmp;
      cy >>= BITS_PER_GROUP;
    }
    *(Nbr+NBR_LIMBS) = 0;
  }
}

void AddBigNbrs(const int *Nbr1, const int *Nbr2, int *Sum)
{
  unsigned int carry = (unsigned int)*Nbr1 + (unsigned int)*Nbr2;
  *Sum = carry & MAX_VALUE_LIMB;
  carry = (carry >> BITS_PER_GROUP) + *(Nbr1+1) + *(Nbr2+1);
  *(Sum+1) = carry & MAX_VALUE_LIMB;
}

void SubtBigNbrs(const int *Nbr1, const int *Nbr2, int *Diff)
{
  unsigned int borrow = (unsigned int)*Nbr1 - (unsigned int)*Nbr2;
  *Diff = borrow & MAX_VALUE_LIMB;
  // On subtraction, borrow subtracts too.
  borrow = (unsigned int)*(Nbr1+1) - (unsigned int)*(Nbr2+1) - (borrow >> BITS_PER_GROUP);
  *(Diff+1) = borrow & MAX_VALUE_LIMB;
}

static void AddBigNbrModN(const int *Nbr1, const int *Nbr2, int *Sum)
{
  unsigned int Sum0;
  unsigned int Sum1;
  unsigned int TestNbrGraphic0 = (unsigned int)TestNbrGraphic[0];
  unsigned int TestNbrGraphic1 = (unsigned int)TestNbrGraphic[1];
  unsigned int carry = (unsigned int)*Nbr1 + (unsigned int)*Nbr2;
  Sum0 = carry & MAX_VALUE_LIMB;
  Sum1 = (carry >> BITS_PER_GROUP) + *(Nbr1 + 1) + *(Nbr2 + 1);
  if ((Sum1 > TestNbrGraphic1) || 
     ((Sum1 == TestNbrGraphic1) && (Sum0 >= TestNbrGraphic0)))
  {
    unsigned int borrow = Sum0 - TestNbrGraphic0;
    Sum0 = borrow & MAX_VALUE_LIMB;
    // On subtraction, borrow subtracts too.
    Sum1 = Sum1 - TestNbrGraphic1 - (borrow >> BITS_PER_GROUP);
  }
  *Sum = (int)Sum0;
  *(Sum+1) = (int)Sum1;
}

static void GetMontgomeryParms(void)
{
  int N;
  int x;

  N = TestNbrGraphic[0];   // 2 least significant bits of inverse correct.
  x = N;
  x = x * (2 - (N * x)); // 4 least significant bits of inverse correct.
  x = x * (2 - (N * x)); // 8 least significant bits of inverse correct.
  x = x * (2 - (N * x)); // 16 least significant bits of inverse correct.
  x = x * (2 - (N * x)); // 32 least significant bits of inverse correct.
  MontMultGraphicN = x & MAX_INT_NBR;
  MontMultGraphicR1[1] = 0;
  if (TestNbrGraphic[1] == 0)
  {
#ifdef __ARM_ARCH_7A__
    MontMultGraphicR1[0] = (MAX_INT_NBR - TestNbrGraphic[0] + 1) %
                            TestNbrGraphic[0]; // 2^31 mod TestNbrGraphic.
#else
    MontMultGraphicR1[0] = 1;
#endif
    return;
  }
  MontMultGraphicR1[2] = 1;
  MontMultGraphicR1[0] = 0;
  AdjustModN(MontMultGraphicR1);  // 2^62 mod TestNbrGraphic.
}

// Perform Miller-Rabin test of number stored in variable TestNbrGraphic.
// The bases to be used are 2, 3, 5, 7, 11, 13, 17, 19 and 23 which
// ensures that any composite less than 3*10^18 is discarded.
// On input: value: array of two ints, first low, then high.
bool isPrime(const int *value)
{
  static const char bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 0};
  // List of lowest composite numbers that passes Miller-Rabin
  // for above bases (OEIS A014233).
  // The first number is the limit modulo 2^31,
  // the second number is the limit divided 2^31.
  static const int limits[] =
  {
    0, 0,

    // Base 2: limit = 2047
    2047, 0,              // Least and most 31 significant bits.

    // Bases 2 and 3: limit = 1373653
    1373653, 0,           // Least and most 31 significant bits.

    // Bases 2, 3 and 5: limit = 25326001
    25326001, 0,          // Least and most 31 significant bits.

    // Bases 2, 3, 5 and 7: limit = 3215031751
    1067548103, 1,        // Least and most 31 significant bits.
 
    // Bases 2, 3, 5, 7 and 11: limit = 2152302898747
    524283451, 1002,      // Least and most 31 significant bits.

    // Bases 2, 3, 5, 7, 11 and 13: limit = 3474749660383
    121117919, 1618,      // Least and most 31 significant bits.

    // Bases 2, 3, 5, 7, 11, 13 and 17: limit = 341550071728321
    1387448513, 159046,   // Least and most 31 significant bits.

    // Bases 2, 3, 5, 7, 11, 13, 17 and 19: limit = 341550071728321
    1387448513, 159046,   // Least and most 31 significant bits.

    // Bases 2, ..., 23: limit = 3825123056546413051 = 2^61.73.
    // Greater than any argument of isPrime()
    0, 2000000000,        // Least and most 31 significant bits.
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
  unsigned int tmp;
  int baseInMontRepres[NBR_LIMBS];
  int power[NBR_LIMBS];
  int temp[NBR_LIMBS];
    // Convert parameter to big number (2 limbs of 31 bits each).
  int TestNbrGraphic0 = value[0];
  int TestNbrGraphic1 = value[1];
  if (TestNbrGraphic1 >= HALF_INT_RANGE)
  {    // Number is negative. Change sign.
    if (TestNbrGraphic0 == 0)
    {
      TestNbrGraphic1 = -TestNbrGraphic1 & (int)MAX_VALUE_LIMB;
    }
    else
    {
      TestNbrGraphic0 = -TestNbrGraphic0 & (int)MAX_VALUE_LIMB;
      TestNbrGraphic1 = (-1 - TestNbrGraphic1) & (int)MAX_VALUE_LIMB;
    }
  }
  TestNbrGraphic[0] = TestNbrGraphic0;
  TestNbrGraphic[1] = TestNbrGraphic1;
  if (TestNbrGraphic1 == 0)
  {
    if (TestNbrGraphic0 == 1)
    {
      return false;           // 1 is not prime.
    }
    if (TestNbrGraphic0 == 2)
    {
      return true;            // 2 is prime.
    }
  }
  if ((TestNbrGraphic0 & 1) == 0)
  {
    return false;             // Even numbers different from 2 are not prime.
  }
  for (i=1; i<(int)sizeof(primes); i++)
  {
    base = primes[i];
    if (TestNbrGraphic1 == 0)
    {
      if ((unsigned int)TestNbrGraphic0 == base)
      {
        return true;          // Number is prime.
      }
      if (((unsigned int)TestNbrGraphic0 % base) == 0U)
      {
        return false;         // Number is multiple of base, so it is composite.
      }
    }
    // Check whether TestNbrGraphic is multiple of base. In this case the number would be composite.
    // No overflow possible in next expression.
    else if (((((unsigned int)TestNbrGraphic1 % base) * (LIMB_RANGE % base) +
      (unsigned int)TestNbrGraphic0) % base) == 0U)
    {
      return false;           // Number is multiple of base, so it is composite.
    }
    else
    {                         // Nothing to do.
    }
  }
  GetMontgomeryParms();
  baseInMontRepres[0] = MontMultGraphicR1[0];
  baseInMontRepres[1] = MontMultGraphicR1[1];
  
       // Find index of least significant bit set disregarding bit 0.
  mask = TestNbrGraphic0 & (MAX_INT_NBR-1);
  indexLSB = 0;
  if (mask == 0U)
  {    // Least significant bit is inside high limb.
    mask = TestNbrGraphic1;
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
  mask = TestNbrGraphic1;
  indexMSB = BITS_PER_GROUP;
  idxNbrMSB = 1;
  if (mask == 0U)
  {    // Most significant bit is inside low limb.
    mask = TestNbrGraphic0;
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
  tmp = 1U;
  maskMSB = (tmp << ((unsigned int)indexMSB % (unsigned int)BITS_PER_GROUP));
  i = 0;
  j = 0;
  while ((limits[j+1] < TestNbrGraphic1) || 
         ((limits[j+1] == TestNbrGraphic1) && (limits[j] <= TestNbrGraphic0)))
  {
    int idxNbr;
    base = bases[i];
    do                     // Compute next base in Montgomery representation.
    {
      AddBigNbrModN(baseInMontRepres, MontMultGraphicR1, baseInMontRepres);
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
      MontMultGraphic(power, power, power);

      if (((unsigned int)TestNbrGraphic[idxNbr] & mask) != 0U)
      {
        MontMultGraphic(power, baseInMontRepres, power);
      }
    }
       // If power equals 1 or -1 in Montgomery representation,
       // another base must be tried.
    if ((power[0] == MontMultGraphicR1[0]) && (power[1] == MontMultGraphicR1[1]))
    {
      continue;   // power equals 1, so another base must be tried.
    }
    AddBigNbrModN(power, MontMultGraphicR1, temp);
    if ((temp[0] == 0) && (temp[1] == 0))
    {
      continue;   // power equals -1, so another base must be tried.
    }
    for (; index>=0; index--)
    {
      mask >>= 1;
      if (mask == 0U)
      {
        mask = HALF_INT_RANGE;
        idxNbr--;
      }
      MontMultGraphic(power, power, power);
      if ((power[0] == MontMultGraphicR1[0]) && (power[1] == MontMultGraphicR1[1]))
      {
        return false;  // power equals 1, so number is composite.
      }
      AddBigNbrModN(power, MontMultGraphicR1, temp);
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

// Multiply 32x32 = 64 bits.
void multiplyBigNbrs(int factor1, int factor2, int *prod)
{
  unsigned int tmp = ((unsigned int)factor1 * (unsigned int)factor2) & MAX_VALUE_LIMB;
#ifdef _USING64BITS_
  uint64_t ui64Tmp = ((uint64_t)factor1 * (uint64_t)factor2) >> BITS_PER_GROUP;
  *prod = (int)tmp;
  *(prod + 1) = (int)ui64Tmp;
#else
  int low = (int)tmp;
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

// Generate the ASCII characters from a 32-bit number.
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

// Store in *pNbrLo and *pNbrHi the 64-bit number represented in ASCII.
void getValue64(const char* value, int* pNbrLo, int* pNbrHi)
{
  bool valueIsNegative = false;
  const char* ptrValue = value;
  int nbrLo = 0;
  int nbrHi = 0;
  if (*ptrValue == '-')
  {
    valueIsNegative = true;
    ptrValue++;
  }
  for (int index = 0; index < 19; index++)
  {
    unsigned int tmp;
    int charConverted;
    double dProd;
    if (*ptrValue == 0)
    {      // End of string, so end of conversion from string to number.
      break;
    }
    charConverted = (*ptrValue - '0');
    ptrValue++;
    dProd = ((double)nbrLo * 10.0) + (double)charConverted;
    nbrLo = (nbrLo * 10) + charConverted;
    tmp = (unsigned int)nbrLo & MAX_VALUE_LIMB;
    nbrLo = (int)tmp;
    nbrHi = (nbrHi * 10) + (int)(dProd / (double)LIMB_RANGE);
  }
  if (valueIsNegative)
  {
    if (nbrLo == 0)
    {
      nbrHi = -nbrHi & (int)MAX_VALUE_LIMB;
    }
    else
    {
      nbrLo = -nbrLo & (int)MAX_VALUE_LIMB;
      nbrHi = (-1 - nbrHi) & (int)MAX_VALUE_LIMB;
    }
  }
  *pNbrHi = nbrHi;
  *pNbrLo = nbrLo;
}

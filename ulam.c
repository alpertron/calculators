/*
  ULAM SPIRAL
  Copyright(C) 2018 Dario Alejandro Alpern

  This program is free software : you can redistribute it and / or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warrantgetny of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.If not, see <http://www.gnu.org/licenses/>.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#define SMALL_NUMBER_BOUND 32768
#ifdef __EMSCRIPTEN__
  #define EXTERNALIZE  __attribute__((visibility("default")))
  #define MAX_WIDTH 2048
  #define pixelXY(x, y) &pixels[y * MAX_WIDTH + x];
  unsigned int pixels[2048*MAX_WIDTH];
#else     // Not Emscripten
  #define EXTERNALIZE	
  #define pixelXY(x, y) (Uint32*)doubleBuffer->pixels + y * width + x;
  #include <SDL.h>
  SDL_Surface *screen, *doubleBuffer;
  int oldXCenter;
  int oldYCenter;
  int oldXFraction;
  int oldYFraction;
  int timer;
  int quit;
#endif

#define NBR_LIMBS  2
#define BITS_PER_GROUP 31
#define LIMB_RANGE (1U<<BITS_PER_GROUP)
#define MAX_INT_NBR ((int)((1U << BITS_PER_GROUP)-1))
#define HALF_INT_RANGE (1 << (BITS_PER_GROUP - 1))
#define MAX_VALUE_LIMB (LIMB_RANGE-1)

#define LIMIT(nbr) (int)(nbr%LIMB_RANGE), (int)(nbr/LIMB_RANGE)
#define MAX_LINES  1000
#define MAX_COLUMNS 2000
char infoText[500];
int TestNbr[NBR_LIMBS];
unsigned int MontgomeryMultN;
int MontgomeryMultR1[NBR_LIMBS+1];  // One more limb required for AdjustModN.
int MontgomeryMultR2[NBR_LIMBS+1];
int showAlgebraic;
int thickness = 3;                  // Number of bits to shift.
int xCenter, xFraction;             // Range of fraction: 0 to (1 << thickness) - 1
int yCenter, yFraction;
int width;
int height;
int startNumber[2] = {1, 0};
char initMultipleArrayCalled;
static int multiple[25][97];
static char primes[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41,
                         43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97 };
void initMultipleArray(void)
{
  int i, j;
  for (i = 0; i<25; i++)
  {
    int k = primes[i];
    for (j = k / 2 + 1; j >= 0; j--)
    {
      multiple[i][j*j%k] = 1;
    }
  }
  initMultipleArrayCalled = 1;
}

// Each limb has 31 bits, so with 2 limbs these routines can use numbers up to 2^62
// This is over 10^18.

// Compute Nbr <- Nbr mod TestNbr.
void AdjustModN(int *Nbr)
{
  int i, carry;
  int TrialQuotient;
  double dNbr, dModulus, dTrialQuotient;
  double dDelta;
  double dVal = 1 / (double)LIMB_RANGE;
  double dSquareLimb = (double)LIMB_RANGE * (double)LIMB_RANGE;

  dModulus = (double)TestNbr[1] + (double)TestNbr[0] * dVal;
  dNbr = (double)Nbr[2] * (double)LIMB_RANGE + (double)Nbr[1] + (double)Nbr[0] * dVal;
  TrialQuotient = (int)(unsigned int)floor(dNbr / dModulus + 0.5);
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
    int low = (*(Nbr + i) - TestNbr[i] * TrialQuotient + carry) & MAX_INT_NBR;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    double dAccumulator = *(Nbr+i) - TestNbr[i] * dTrialQuotient + carry + dDelta;
    dDelta = 0;
    if (dAccumulator < 0)
    {
      dAccumulator += dSquareLimb;
      dDelta = -(double)LIMB_RANGE;
    }
    if (low < HALF_INT_RANGE)
    {
      carry = (int)floor((dAccumulator + HALF_INT_RANGE / 2)*dVal);
    }
    else
    {
      carry = (int)floor((dAccumulator - HALF_INT_RANGE / 2)*dVal);
    }
    *(Nbr + i) = low;
  }
  *(Nbr + NBR_LIMBS) = (*(Nbr + NBR_LIMBS) + carry) & MAX_INT_NBR;
  if ((*(Nbr+NBR_LIMBS) & MAX_VALUE_LIMB) != 0)
  {
    unsigned int cy = 0;
    for (i = 0; i < NBR_LIMBS; i++)
    {
      cy += (unsigned int)*(Nbr+i) + (unsigned int)TestNbr[i];
      *(Nbr+i) = (int)(cy & MAX_VALUE_LIMB);
      cy >>= BITS_PER_GROUP;
    }
    *(Nbr+NBR_LIMBS) = 0;
  }
}
// This routine is a lot slower than multiplying by using Montgomery representation,
// so it is used only when converting from normal representation to Montgomery.
void MultBigNbrModN(const int *factor1, const int *factor2, int *Product)
{
  int Prod[4];
  int i;
  double dInvLimbRange = (double)1 / LIMB_RANGE;
  Prod[0] = Prod[1] = Prod[2] = Prod[3] = 0;
  for (i=0; i<NBR_LIMBS; i++)
  {
    unsigned int carry;
    // Multiply least or most significant limb of first factor by least significant limb of second factor.
    int Nbr = *(factor1+i);
    double dNbr = (double)Nbr;
    int low = (Nbr * *factor2 + Prod[i]) & MAX_VALUE_LIMB;
    double dAccum = dNbr * (double)*factor2 + (double)Prod[i];
    Prod[i] = low;
    if (low < HALF_INT_RANGE)
    {
      dAccum = ((dAccum + HALF_INT_RANGE / 2)*dInvLimbRange);
    }
    else
    {
      dAccum = ((dAccum - HALF_INT_RANGE / 2)*dInvLimbRange);
    }
    carry = (unsigned int)dAccum;  // Most significant limb can be greater than LIMB_RANGE

    // Multiply least or most significant limb of first factor by most significant limb of second factor.
    low = (Nbr * *(factor2+1) + Prod[i+1] + carry) & MAX_VALUE_LIMB;
    Prod[i+1] = low;
    dAccum = dNbr * (double)*(factor2 + 1) + (double)Prod[i+1] + (double)carry;
    if (low < HALF_INT_RANGE)
    {
      dAccum = ((dAccum + HALF_INT_RANGE / 2)*dInvLimbRange);
    }
    else
    {
      dAccum = ((dAccum - HALF_INT_RANGE / 2)*dInvLimbRange);
    }
    Prod[i+2] = (unsigned int)dAccum;  // Most significant limb can be greater than LIMB_RANGE
  }
   // Modular adjustment.
  AdjustModN(&Prod[1]);
  AdjustModN(Prod);
  *Product = Prod[0];
  *(Product+1) = Prod[1];
  *(Product+2) = Prod[2];
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
    int remainder = factor1 * factor2 - quotient * mod;
    if (remainder < 0)
    {    // Quotient was 1 more than expected. Adjust remainder.
      remainder += mod;
    }
    *product = remainder;
#endif
  }
}

void MontgomeryMult(const int *factor1, const int *factor2, int *Product)
{
#ifndef _USING64BITS_
  int carry;
#endif  
  if (TestNbr[1] == 0)
  {
    smallmodmult(*factor1, *factor2, Product, TestNbr[0]);
    *(Product+1) = 0;
    return;
  }
  int TestNbr0 = TestNbr[0];
  int TestNbr1 = TestNbr[1];
  uint32_t Prod0, Prod1;
  int factor2_0 = *factor2;
  int factor2_1 = *(factor2+1);
#ifdef _USING64BITS_
  uint64_t Pr;
  unsigned int Nbr, MontDig;
  
  Pr = (Nbr = *factor1) * (uint64_t)factor2_0;
  MontDig = ((uint32_t)Pr * MontgomeryMultN) & MAX_INT_NBR;
  Prod0 = (Pr = (((uint64_t)MontDig * TestNbr0 + Pr) >> BITS_PER_GROUP) +
    (uint64_t)MontDig * TestNbr1 + (uint64_t)Nbr * factor2_1) & MAX_INT_NBR;
  Prod1 = (uint32_t)(Pr >> BITS_PER_GROUP);
    
  Pr = (Nbr = *(factor1 + 1)) * (uint64_t)factor2_0 + (uint32_t)Prod0;
  MontDig = ((uint32_t)Pr * MontgomeryMultN) & MAX_INT_NBR;
  Prod0 = (Pr = (((uint64_t)MontDig * TestNbr0 + Pr) >> BITS_PER_GROUP) +
    (uint64_t)MontDig * TestNbr1 + (uint64_t)Nbr * factor2_1 + (uint32_t)Prod1) & MAX_INT_NBR;
  Prod1 = (uint32_t)(Pr >> BITS_PER_GROUP);
    
  if (Pr >= ((uint64_t)(TestNbr1 + 1) << BITS_PER_GROUP) || (Prod1 == TestNbr1 && Prod0 >= TestNbr0))
  {
    int32_t borrow;
    Prod0 = (borrow = (int32_t)Prod0 - (int32_t)TestNbr0) & MAX_INT_NBR;
    Prod1 = ((borrow >> BITS_PER_GROUP) + (int32_t)Prod1 - (int32_t)TestNbr1) & MAX_INT_NBR;
  }
#else
  double dInvLimbRange = (double)1 / (double)LIMB_RANGE;
  int Nbr = *(factor1);
  double dNbr = (double)Nbr;
  int low = Nbr * factor2_0;
  double dAccum = dNbr * (double)factor2_0;
  int MontDig = (low * MontgomeryMultN) & MAX_VALUE_LIMB;
  double dMontDig = (double)MontDig;
  dAccum += dMontDig * (double)TestNbr0;
  // At this moment dAccum is multiple of LIMB_RANGE.
  dAccum = floor(dAccum*dInvLimbRange + 0.5);
  low = ((unsigned int)dAccum + MontDig * TestNbr1 +
               Nbr * factor2_1) & MAX_VALUE_LIMB;
  dAccum += dMontDig * TestNbr1 + dNbr * factor2_1;
  Prod0 = low;
  if (low < HALF_INT_RANGE)
  {
    dAccum = ((dAccum + HALF_INT_RANGE / 2)*dInvLimbRange);
  }
  else
  {
    dAccum = ((dAccum - HALF_INT_RANGE / 2)*dInvLimbRange);
  }
  Prod1 = (unsigned int)dAccum;  // Most significant limb can be greater than LIMB_RANGE
  
  Nbr = *(factor1 + 1);
  dNbr = (double)Nbr;
  low = Nbr * factor2_0 + Prod0;
  dAccum = dNbr * (double)factor2_0 + (double)Prod0;
  MontDig = (low * MontgomeryMultN) & MAX_VALUE_LIMB;
  dMontDig = (double)MontDig;
  dAccum += dMontDig * (double)TestNbr0;
  // At this moment dAccum is multiple of LIMB_RANGE.
  dAccum = floor(dAccum*dInvLimbRange + 0.5);
  low = ((unsigned int)dAccum + MontDig * TestNbr1 +
               Nbr * factor2_1 + Prod1) & MAX_VALUE_LIMB;
  dAccum += dMontDig * TestNbr1 + dNbr * factor2_1 + (unsigned int)Prod1;
  Prod0 = low;
  if (low < HALF_INT_RANGE)
  {
    dAccum = ((dAccum + HALF_INT_RANGE / 2)*dInvLimbRange);
  }
  else
  {
    dAccum = ((dAccum - HALF_INT_RANGE / 2)*dInvLimbRange);
  }
  Prod1 = (unsigned int)dAccum;  // Most significant limb can be greater than LIMB_RANGE
  
  if ((unsigned int)Prod1 > (unsigned int)TestNbr1 ||
       ((unsigned int)Prod1 == (unsigned int)TestNbr1 && (unsigned int)Prod0 >= (unsigned int)TestNbr0))
  {        // Prod >= TestNbr, so perform Prod <- Prod - TestNbr
    carry = Prod0 - TestNbr0;
    Prod0 = carry & MAX_VALUE_LIMB;
    Prod1 = ((carry >> BITS_PER_GROUP) + Prod1 - TestNbr1) & MAX_VALUE_LIMB;
  }
#endif  
  *Product = Prod0;
  *(Product+1) = Prod1;
}

void AddBigNbr(const int *Nbr1, const int *Nbr2, int *Sum)
{
  unsigned int carry = *(Nbr1) + *(Nbr2);
  *(Sum) = carry & MAX_INT_NBR;
  carry = (carry >> BITS_PER_GROUP) + *(Nbr1+1) + *(Nbr2+1);
  *(Sum+1) = carry & MAX_INT_NBR;
}

void SubtBigNbr(const int *Nbr1, const int *Nbr2, int *Diff)
{
  int borrow = *(Nbr1) - *(Nbr2);
  *(Diff) = borrow & MAX_INT_NBR;
  borrow = (borrow >> BITS_PER_GROUP) + *(Nbr1+1) - *(Nbr2+1);
  *(Diff+1) = borrow & MAX_INT_NBR;
}

void AddBigNbrModN(const int *Nbr1, const int *Nbr2, int *Sum)
{
  int Sum0, Sum1;
  int TestNbr0 = TestNbr[0];
  int TestNbr1 = TestNbr[1];
  unsigned int carry = *Nbr1 + *Nbr2;
  Sum0 = carry & MAX_INT_NBR;
  carry = (carry >> BITS_PER_GROUP) + *(Nbr1 + 1) + *(Nbr2 + 1);
  Sum1 = carry & MAX_INT_NBR;
  if (carry > (unsigned int)TestNbr1 || (carry == (unsigned int)TestNbr1 && Sum0 >= TestNbr0))
  {
    int borrow = Sum0 - TestNbr0;
    Sum0 = borrow & MAX_INT_NBR;
    Sum1 = ((borrow >> BITS_PER_GROUP) + Sum1 - TestNbr1) & MAX_INT_NBR;
  }
  *Sum = Sum0;
  *(Sum+1) = Sum1;
}

void GetMontgomeryParms(void)
{
  int N, x;

  MontgomeryMultR1[1] = 0;
  if (TestNbr[1] == 0)
  {
    MontgomeryMultR1[0] = 1;
    return;
  }
  x = N = (int) TestNbr[0]; // 2 least significant bits of inverse correct.
  x = x * (2 - N * x); // 4 least significant bits of inverse correct.
  x = x * (2 - N * x); // 8 least significant bits of inverse correct.
  x = x * (2 - N * x); // 16 least significant bits of inverse correct.
  x = x * (2 - N * x); // 32 least significant bits of inverse correct.
  MontgomeryMultN = (-x) & MAX_INT_NBR;
  MontgomeryMultR1[2] = 1;
  MontgomeryMultR1[0] = 0;
  AdjustModN(MontgomeryMultR1);
}

// Perform Miller-Rabin test of number stored in variable TestNbr.
// The bases to be used are 2, 3, 5, 7, 11, 13, 17, 19 and 23 which
// ensures that any composite less than 3*10^18 is discarded.
int isPrime(const int *value)
{
  static char bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 0};
  // List of lowest composite numbers that passes Miller-Rabin for above bases (OEIS A014233).
  static int limits[] =
  {
    LIMIT(0),
    LIMIT(2047ll),                 // Base 2
    LIMIT(1373653ll),              // Bases 2 and 3
    LIMIT(25326001ll),             // Bases 2, 3 and 5
    LIMIT(3215031751ll),           // Bases 2, 3, 5 and 7
    LIMIT(2152302898747ll),        // Bases 2, 3, 5, 7 and 11
    LIMIT(3474749660383ll),        // Bases 2, 3, 5, 7, 11 and 13
    LIMIT(341550071728321ll),      // Bases 2, 3, 5, 7, 11, 13 and 17
    LIMIT(341550071728321ll),      // Bases 2, 3, 5, 7, 11, 13, 17 and 19
    LIMIT((1ll << (2*BITS_PER_GROUP)) - 1)
  };
  int i, j;
  int index, indexLSB, indexMSB, idxNbrMSB;
  unsigned int mask, maskMSB;
  unsigned int prevBase = 1;
  unsigned int base;
  int baseInMontRepres[NBR_LIMBS];
  int power[NBR_LIMBS];
  int temp[NBR_LIMBS];
    // Convert parameter to big number (2 limbs of 31 bits each).
  int TestNbr0 = TestNbr[0] = value[0];
  int TestNbr1 = TestNbr[1] = value[1];
  if (TestNbr1 == 0)
  {
    if (TestNbr0 == 1)
    {
      return 0;            // 1 is not prime.
    }
    else if (TestNbr0 == 2)
    {
      return 1;            // 2 is prime.
    }
  }
  if ((TestNbr0 & 1) == 0)
  {
    return 0;              // Even numbers different from 2 are not prime.
  }
  for (i=1; i<sizeof(primes); i++)
  {
    base = primes[i];
    if (TestNbr1 == 0)
    {
      if (TestNbr0 == base)
      {
        return 1;          // Number is prime.
      }
      if (TestNbr0 % base == 0)
      {
        return 0;          // Number is multiple of base, so it is composite.
      }
    }
    // Check whether TestNbr is multiple of base. In this case the number would be composite.
    else if (((TestNbr1 % base) * (LIMB_RANGE % base) + TestNbr0 % base) % base == 0)  // No overflow possible.
    {
      return 0;            // Number is multiple of base, so it is composite.
    }
  }
  GetMontgomeryParms();
  baseInMontRepres[0] = MontgomeryMultR1[0];
  baseInMontRepres[1] = MontgomeryMultR1[1];
  
       // Find index of least significant bit set disregarding bit 0.
  mask = TestNbr0 & (MAX_INT_NBR-1);
  indexLSB = 0;
  if (mask == 0)
  {    // Least significant bit is inside high limb.
    mask = TestNbr1;
    indexLSB = BITS_PER_GROUP;
  }
  if ((mask & 0xFFFF) == 0)
  {    // Least significant bit is somewhere between bits 31-16.
    mask >>= 16;
    indexLSB += 16;
  }
  if ((mask & 0xFF) == 0)
  {    // Least significant bit is somewhere between bits 15-8.
    mask >>= 8;
    indexLSB += 8;
  }
  if ((mask & 0x0F) == 0)
  {    // Least significant bit is somewhere between bits 7-4.
    mask >>= 4;
    indexLSB += 4;
  }
  if ((mask & 0x03) == 0)
  {    // Least significant bit is between bits 3-2.
    mask >>= 2;
    indexLSB += 2;
  }
  if ((mask & 0x01) == 0)
  {    // Least significant bit is bit 1.
    indexLSB++;
  }
      // Find index of most significant bit set.
  mask = TestNbr1;
  indexMSB = BITS_PER_GROUP;
  idxNbrMSB = 1;
  if (mask == 0)
  {    // Most significant bit is inside low limb.
    mask = TestNbr0;
    indexMSB = 0;
    idxNbrMSB = 0;
  }
  if (mask & 0xFFFF0000)
  {    // Most significant bit is somewhere between bits 31-16.
    mask >>= 16;
    indexMSB += 16;
  }
  if (mask & 0xFF00)
  {    // Most significant bit is somewhere between bits 15-8.
    mask >>= 8;
    indexMSB += 8;
  }
  if (mask & 0xF0)
  {    // Most significant bit is somewhere between bits 7-4.
    mask >>= 4;
    indexMSB += 4;
  }
  if (mask & 0x0C)
  {    // Most significant bit is between bits 3-2.
    mask >>= 2;
    indexMSB += 2;
  }
  if (mask & 0x02)
  {    // Most significant bit is bit 1.
    indexMSB++;
  }
  maskMSB = (1<<(indexMSB % BITS_PER_GROUP));
  
  for (i=0, j=0; limits[j+1] < TestNbr1 || (limits[j+1] == TestNbr1 && limits[j] < TestNbr0); i++, j+=2)
  {
    int idxNbr;
    base = bases[i];
    do                     // Compute next base in Montgomery representation.
    {
      AddBigNbrModN(baseInMontRepres, MontgomeryMultR1, baseInMontRepres);
    } while (++prevBase < base);

    power[0] = baseInMontRepres[0];
    power[1] = baseInMontRepres[1];
    mask = maskMSB;
    idxNbr = idxNbrMSB;
    for (index = indexMSB-1; index >= indexLSB; index--)
    {
      mask >>= 1;
      if (mask == 0)
      {
        mask = HALF_INT_RANGE;
        idxNbr--;
      }
      MontgomeryMult(power, power, power);

      if (TestNbr[idxNbr] & mask)
      {
        MontgomeryMult(power, baseInMontRepres, power);
      }
    }
       // If power equals 1 or -1 in Montgomery representation,
       // another base must be tried.
    if (power[0] == MontgomeryMultR1[0] && power[1] == MontgomeryMultR1[1])
    {
      continue;   // power equals 1, so another base must be tried.
    }
    AddBigNbrModN(power, MontgomeryMultR1, temp);
    if (temp[0] == 0 && temp[1] == 0)
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
      if (power[0] == MontgomeryMultR1[0] && power[1] == MontgomeryMultR1[1])
      {
        return 0;  // power equals 1, so number is composite.
      }
      AddBigNbrModN(power, MontgomeryMultR1, temp);
      if (temp[0] == 0 && temp[1] == 0)
      {            // power equals -1.
        break;
      }
    }
    if (index >= 0)
    {
      continue;    // power equals -1, so another base must be tried.
    }
    return 0;      // base^(n-1) != 1 (mod n) is not 1, so number is composite.
  }
  return 1;        // All Miller-Rabin tests were passed, so number is prime.
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
    dAccum = ((dAccum + HALF_INT_RANGE / 2) / LIMB_RANGE);
  }
  else
  {
    dAccum = ((dAccum - HALF_INT_RANGE / 2) / LIMB_RANGE);
  }
  *(prod+1) = (unsigned int)dAccum;
#endif 
}

// Check whether 4x^2 + bx + c has integer roots (b is a very small number).
int algebraicFactor(int linear, int *indep)
{
  int temp[2];
  int iSqDelta, t1, t2;
  double dDelta, dFourAC;
  temp[0] = *indep;
  temp[1] = *(indep+1);
  if ((temp[1] & HALF_INT_RANGE))
  {    // Independent term is negative.
    int carry = -temp[0];
    temp[0] = carry & MAX_INT_NBR;
    carry = (carry >> BITS_PER_GROUP) - temp[1];
    temp[1] = carry & MAX_INT_NBR;
  }
  dFourAC = ((double)temp[1] * MAX_VALUE_LIMB + (double)temp[0]) * 16;
  if (*(indep + 1) & HALF_INT_RANGE)
  {    // Independent term is negative.
    dDelta = linear * linear + dFourAC;
  }
  else
  {    // Independent term is positive.
    dDelta = linear * linear - dFourAC;
  }
  if (dDelta < 0)
  {
    return 0;   // No real roots.
  }
  iSqDelta = (int)(sqrt(dDelta) + 0.5);              // Convert to nearest integer.
  t1 = (int)(((double)(-linear) + iSqDelta) / 4);    // t1 and t2 are less than 2^31 so
  t2 = (int)(((double)(-linear) - iSqDelta) / 4);    // they fit into a double.

  multiply(2*t1 + linear, t1, temp);
  AddBigNbr(temp, indep, temp);
  AddBigNbr(temp, indep, temp);
  if (temp[0] != 0 || temp[1] != 0)
  {
    return 0;
  }
  multiply(2*t2 + linear, t2, temp);
  AddBigNbr(temp, indep, temp);
  AddBigNbr(temp, indep, temp);
  if (temp[0] != 0 || temp[1] != 0)
  {
    return 0;
  }
  return 1;  // Twice the roots are integer numbers
}

void getN(int x, int y, int *value)
{
  int addend[2];
  if (x >= 0 && x >= y && x >= -y)
  {                     // Right quadrant.
    multiply(4 * x + 3, x, value);
    addend[0] = y & MAX_INT_NBR;
    addend[1] = (y >> BITS_PER_GROUP) & MAX_INT_NBR;
  }
  else if (y >= 0 && y>x && y>-x)
  {                     // Top quadrant.
    multiply(4 * y - 3, y, value);
    x = -x;
    addend[0] = x & MAX_INT_NBR;
    addend[1] = (x >> BITS_PER_GROUP) & MAX_INT_NBR;
  }
  else if (x <= 0 && x <= y && x <= -y)
  {                     // Left quadrant.
    multiply(4 * x + 1, x, value);
    y = -y;
    addend[0] = y & MAX_INT_NBR;
    addend[1] = (y >> BITS_PER_GROUP) & MAX_INT_NBR;
  }
  else
  {                     // Bottom quadrant.
    multiply(4 * y - 1, y, value);
    addend[0] = x & MAX_INT_NBR;
    addend[1] = (x >> BITS_PER_GROUP) & MAX_INT_NBR;
  }
  AddBigNbr(value, addend, value);
  AddBigNbr(value, startNumber, value);
}

void AddTwoLimbsPlusOneLimb(const int *addend1, int addend2, int *sum)
{
  int temp[2];
  temp[1] = 0;
  if (addend2 >= 0)
  {
    temp[0] = addend2;
    AddBigNbr(addend1, temp, sum);
  }
  else
  {
    temp[0] = -addend2;
    SubtBigNbr(addend1, temp, sum);
  }
}

void setPoint(int x, int y)
{
  int value[2];
  int xPhysical, yPhysical, absx, absy, currY;
  int row, col, firstRow, firstCol, lastRow, lastCol;
  int t, linear1, linear2;
  int indep1[2], indep2[2];
  int Nminustt[2];
  unsigned int algebraicColor;
  unsigned int color;
#ifdef __EMSCRIPTEN__  
  unsigned int colorBlack = 0xFF000000U;
  unsigned int colorBlue = 0xFFFF0000U;
  unsigned int colorGreen = 0xFF00C000U;
  unsigned int colorWhite = 0xFFC0C0C0U;
#else
  Uint32 colorBlack = SDL_MapRGBA(doubleBuffer->format, 0, 0, 0, 255);
  Uint32 colorBlue = SDL_MapRGBA(doubleBuffer->format, 0, 0, 255, 255);
  Uint32 colorGreen = SDL_MapRGBA(doubleBuffer->format, 0, 192, 0, 255);
  Uint32 colorWhite = SDL_MapRGBA(doubleBuffer->format, 192, 192, 192, 255);
#endif
  unsigned int *ptrPixel;
  getN(x, y, value);
  if (showAlgebraic)
  {                          // Color blue.
    algebraicColor = colorBlue;
  }
  else
  {                          // Color black.
    algebraicColor = colorBlack;
  }
  xPhysical = width / 2 + ((x - xCenter) << thickness) - xFraction;
  yPhysical = height / 2 - ((y - yCenter) << thickness) + yFraction;
  if ((value[1] != 0 || value[0] != 2) && (value[0] & 1) == 0)
  {     // value is not 2 and it is even
    color = algebraicColor;
  }
  else if (isPrime(value))
  {                          // Color green.
    color = colorGreen;
  }
  else
  {
    if (x > y && x > -y)
    {
      t = x;
    }
    else if (x < y && x < -y)
    {
      t = -x;
    }
    else if (y > 0)
    {
      t = y;
    }
    else
    {
      t = -y;
    }
    t = t + t;
    multiply(t, t, Nminustt);
    SubtBigNbr(value, Nminustt, Nminustt);
    if (x + y >= 0)
    {
      if (x - y >= 0)
      {                                 // Right quadrant
        AddTwoLimbsPlusOneLimb(Nminustt, -2*t, indep1);  // n - 4 * t*t - 4 * t
        AddTwoLimbsPlusOneLimb(Nminustt, -t, indep2);    // n - 4 * t*t - 2 * t
        linear1 = 4;
        linear2 = 2;
      }
      else
      {                                 // Upper quadrant
        AddTwoLimbsPlusOneLimb(Nminustt, 2*t, indep1);  // n - 4 * t*t + 4 * t
        AddTwoLimbsPlusOneLimb(Nminustt, t, indep2);    // n - 4 * t*t + 2 * t
        linear1 = -4;
        linear2 = -2;
      }
    }
    else
    {
      indep1[0] = Nminustt[0];                          // n - 4 * t*t
      indep1[1] = Nminustt[1];
      linear1 = 0;
      if (x - y >= 0)
      {                                 // Lower quadrant
        AddTwoLimbsPlusOneLimb(Nminustt, -t, indep2);   // n - 4 * t*t - 2 * t
        linear2 = 2;
      }
      else
      {                             /* Left quadrant */
        AddTwoLimbsPlusOneLimb(Nminustt, t, indep2);    // n - 4 * t*t + 2 * t
        linear2 = -2;
      }
    }
    if (algebraicFactor(linear1, indep1) || algebraicFactor(linear2, indep2))
    {
      color = algebraicColor;
    }
    else
    {                          // Color black.
      color = colorBlack;
    }
  }
  firstCol = (xPhysical<0? 0: xPhysical);
  firstRow = (yPhysical<0? 0: yPhysical);
  lastCol = xPhysical + (1 << thickness);
  if (lastCol > width)
  {
    lastCol = width;
  }
  lastRow = yPhysical + (1 << thickness);
  if (lastRow > height)
  {
    lastRow = height;
  }
  for (row = firstRow; row < lastRow; row++)
  {
    ptrPixel = pixelXY(firstCol, row);
    for (col = firstCol; col < lastCol; col++)
    {
      *ptrPixel++ = color;
    }
  }
  if (thickness >= 2)
  {
    color = colorWhite;
    absx = (x>0 ? x : -x);
    absy = (y>0 ? y : -y);
    if (absx >= absy && !(x == -y && y<0))
    {
      if (xPhysical >= 0 && xPhysical < width)
      {
        ptrPixel = pixelXY(xPhysical, firstRow);
        for (row = firstRow; row < lastRow; row++)
        {
          *ptrPixel = color;
#ifdef __EMSCRIPTEN__
          ptrPixel += MAX_WIDTH;
#else
          ptrPixel += width;
#endif    
        }
      }
    }
    if ((absx < absy && !(x == y - 1 && y>0)) || (absx == absy && y <= 0))
    {
      currY = yPhysical + (1 << thickness) - 1;
      if (currY >= 0 && currY < height)
      {
        ptrPixel = pixelXY(firstCol, currY);
        for (col = firstCol; col < lastCol; col++)
        {
          *ptrPixel++ = color;
        }
      }
    }
  }
}     /* end method setPoint */

EXTERNALIZE void drawPartialUlamSpiral(int xminDisp, int xmaxDisp, int yminDisp, int ymaxDisp)
{
  int x, y;
  if (initMultipleArrayCalled == 0)
  {
    initMultipleArray();
  }
  int xmin = xCenter + ((xFraction + xminDisp) >> thickness);
  int xmax = xCenter + ((xFraction + xmaxDisp) >> thickness);
  int ymin = yCenter - ((-yFraction - yminDisp) >> thickness);
  int ymax = yCenter - ((-yFraction - ymaxDisp) >> thickness);

  for (x = xmin; x <= xmax; x++)
  {
    for (y = ymin; y <= ymax; y++)
    {
      setPoint(x, y);
    }  /* end for y */
  }    /* end for x */
}      /* end method drawUlamSpiral */

#ifdef __EMSCRIPTEN__
char *(void)strcpy(char *dest, const char *src)
{
  char *dst = dest;
  while (*src)
  {
    *dst++ = *src++;
  }
  *dst = 0;
  return dest;
}

size_t strlen(const char *s)
{
  const char *a = s;
  while (*a)
  {
    a++;
  }
  return a - s;
}

unsigned int *getPixels(void)
{
  return pixels;
}
#endif

char *appendInt(char *text, int value)
{
  int div = 1000000000;
  int zeroIsSignificant = 0;
  if (value < 0)
  {
    value = -value;
    *text++ = '-';
  }
  do
  {
    int quot = value / div;
    if (quot != 0 || zeroIsSignificant)
    {
      zeroIsSignificant = 1;
      *text++ = (char)quot + '0';
      value -= quot*div;
    }
    div /= 10;
  } while (div > 1);
  *text++ = (char)value + '0';
  *text = 0;
  return text;
}

// Convert the number from binary to string.
char *appendInt64(char *text, const int *value)
{
  int index, index2;
  int nbr0 = *value;
  int nbr1 = *(value+1);
  // Convert digits from right to left.
  text += 19;
  for (index = 0; index < 19; index++)
  {
    // Divide value by 10 and obtain the remainder.
    double dDivid;
    int rem = nbr1 % 10;
    nbr1 = nbr1 / 10;
    dDivid = (double)rem * (double)LIMB_RANGE + (double)nbr0;
    nbr0 = (int)(dDivid / 10);
    *(--text) = (int)(dDivid - (double)nbr0*10) + '0';
  }
  for (index = 0; index < 18; index++)
  {
    if (*(text+index) != '0')
    {
      break;    // Most significant digit found.
    }
  }
  for (index2 = index; index2 < 19; index2++)
  {
    *(text+index2-index) = *(text+index2);
  }
  text += index2-index;
  *text = 0;    // Set string terminator.
  return text;
}

void ShowLabel(char *text, int b, int *indep)
{
  int p;
  int temp[2];
  int carry;
  char *ptrText = &infoText[strlen(infoText)];
  int firstTime = 1;
  (void)strcpy(ptrText, text);
  ptrText += strlen(ptrText);
  temp[0] = *indep;
  temp[1] = *(indep+1);
  if (temp[1] != 0 || temp[0] != 0)
  {      // Independent term is not zero.
    *ptrText++ = ' ';
    if ((temp[1] & HALF_INT_RANGE) == 0)
    {    // Independent term is positive.
      *ptrText++ = '+';
      *ptrText++ = ' ';
    }
    else
    {    // Independent term is negative.
      *ptrText++ = '-';
      *ptrText++ = ' ';
         // Change its sign.
      carry = -temp[0];
      temp[0] = carry & MAX_INT_NBR;
      carry = (carry >> BITS_PER_GROUP) - temp[1];
      temp[1] = carry & MAX_INT_NBR;
    }
    ptrText = appendInt64(ptrText, temp);
  }
  if (b != 0 || temp[0] != 0 || temp[1] != 0)
  {      // linear and independent terms are not zero.
    if (algebraicFactor(b, indep))
    {
      double dDelta;
      double dFourAC = ((double)temp[1] * MAX_VALUE_LIMB + (double)temp[0]) * 16;
      if (*(indep + 1) & HALF_INT_RANGE)
      {    // Independent term is negative.
        dDelta = b * b + dFourAC;
      }
      else
      {    // Independent term is positive.
        dDelta = b * b - dFourAC;
      }
      int iSqDelta = (int)(sqrt(dDelta) + 0.5);              // Convert to nearest integer.
      int t1 = (int)(((double)(-b) + iSqDelta) / 4);         // t1 and t2 are less than 2^31 so
      int t2 = (int)(((double)(-b) - iSqDelta) / 4);         // they fit into a double.
                                                            /* Twice the roots are integer numbers */
      ptrText += strlen(ptrText);
      if (t1 > 0)
      {
        (void)strcpy(ptrText, " = (2t + ");
        ptrText += strlen(ptrText);
        ptrText = appendInt(ptrText, t1);
        *ptrText++ = ')';
      }
      else if (t1 < 0)
      {
        (void)strcpy(ptrText, " = (2t - ");
        ptrText += strlen(ptrText);
        ptrText = appendInt(ptrText, -t1);
        *ptrText++ = ')';
      }
      else
      {
        (void)strcpy(ptrText, " = 2t");
        ptrText += strlen(ptrText);
      }
      if (t2 > 0)
      {
        (void)strcpy(ptrText, " (2t + ");
        ptrText += strlen(ptrText);
        ptrText = appendInt(ptrText, t2);
        *ptrText++ = ')';
      }
      else if (t2 < 0)
      {
        (void)strcpy(ptrText, " (2t - ");
        ptrText += strlen(ptrText);
        ptrText = appendInt(ptrText, -t2);
        *ptrText++ = ')';
      }
      else
      {
        (void)strcpy(ptrText, " 2t");
        ptrText += strlen(ptrText);
      }
    }
    else
    {
      temp[0] = *indep;
      temp[1] = *(indep+1);
      if ((temp[0] & 1) == 0)
      {           // Independent term is even
        if ((b & 3) == 0 && (temp[0] & 3) == 0)
        {         // Both linear and independent term are multiple of 4.
          (void)strcpy(ptrText, " = 4 (t<sup>2</sup>");
          b /= 4;
          temp[0] = ((temp[0] >> 2) | (temp[1] << (BITS_PER_GROUP-2))) & MAX_INT_NBR;
          temp[1] = (temp[1] << (32-BITS_PER_GROUP)) >> (2 + (32 - BITS_PER_GROUP));   // Divide by 4.
        }
        else
        {
          (void)strcpy(ptrText, " = 2 (2t<sup>2</sup>");
          b /= 2;
          temp[0] = ((temp[0] >> 1) | (temp[1] << (BITS_PER_GROUP-1))) & MAX_INT_NBR;
          temp[1] = (temp[1] << (32-BITS_PER_GROUP)) >> (1 + (32 - BITS_PER_GROUP));   // Divide by 2.
        }
        ptrText += strlen(ptrText);
        if (b != 0)
        {
          if (b < 0)
          {
            (void)strcpy(ptrText, " - ");
            b = -b;
          }
          else
          {
            (void)strcpy(ptrText, " + ");
          }
          ptrText += strlen(ptrText);
          if (b != 1)
          {
            ptrText = appendInt(ptrText, b);
          }
          (void)strcpy(ptrText, "t");
          ptrText++;
        }
        if (temp[1] & HALF_INT_RANGE)
        {     // Independent term is negative.
          *ptrText++ = ' ';
          *ptrText++ = '-';
          *ptrText++ = ' ';
              // Change sign.
          carry = -temp[0];
          temp[0] = carry & MAX_INT_NBR;
          carry = (carry >> BITS_PER_GROUP) - temp[1];
          temp[1] = carry & MAX_INT_NBR;
          ptrText = appendInt64(ptrText, temp);
        }
        else
        {
          *ptrText++ = ' ';
          *ptrText++ = '+';
          *ptrText++ = ' ';
          ptrText = appendInt64(ptrText, temp);
        }
        ptrText += strlen(ptrText);
        (void)strcpy(ptrText, ")");
        ptrText++;
      }
      else
      {
        int i;
        for (i = 0; i<sizeof(primes)/sizeof(primes[0]); i++)
        {
          int deltaModP, indepModP;
          p = primes[i];
          // delta = b*b-16*indep
          // Find delta mod p.
          indepModP = ((*(indep+1)%p)*(LIMB_RANGE%p) + *(indep)%p) % p;
          deltaModP = (((b%p) * (b%p) - 16*indepModP) % p + p) % p;
          if (multiple[i][deltaModP] == 0)
          {
            if (firstTime)
            {
              firstTime = 0;
              *ptrText++ = ' ';
              *ptrText++ = '(';
            }
            else
            {
              *ptrText++ = ',';
              *ptrText++ = ' ';
            }
            ptrText = appendInt(ptrText, p);
          }
        }         /* end for */
        if (firstTime == 0)
        {
          (void)strcpy(ptrText, ")");
          ptrText++;
        }
      }           /* end if */
    }             /* end if */
  }               /* end if */
  (void)strcpy(ptrText, "<br>");
}

EXTERNALIZE char *getInformation(int x, int y)
{
  int value[2];
  int yLogical;
  char *ptrText = infoText;

  infoText[0] = 0;   // Empty string.
  if (x>=0)
  {
    int t;
    int Nminustt[2];
    int indep[2];
    int xLogical = xCenter + ((xFraction + x - width / 2) >> thickness);
    yLogical = yCenter + 1 + ((yFraction - y + height / 2) >> thickness);
    getN(xLogical, yLogical, value);
    if (xLogical > yLogical && xLogical > -yLogical)
    {
      t = xLogical;
    }
    else if (xLogical < yLogical && xLogical < -yLogical)
    {
      t = -xLogical;
    }
    else if (yLogical > 0)
    {
      t = yLogical;
    }
    else
    {
      t = -yLogical;
    }
    t = t + t;
    multiply(t, t, Nminustt);
    SubtBigNbr(value, Nminustt, Nminustt);
    if (xLogical + yLogical >= 0)
    {
      if (xLogical - yLogical >= 0)
      {                              // Right quadrant
        AddTwoLimbsPlusOneLimb(Nminustt, -2*t, indep);  // n - 4 * t*t - 4 * t
        ShowLabel("SW-NE: 4t<sup>2</sup> + 4t", 4, indep);
        AddTwoLimbsPlusOneLimb(Nminustt, -t, indep);    // n - 4 * t*t - 2 * t
        ShowLabel("NW-SE: 4t<sup>2</sup> + 2t", 2, indep);
      }
      else
      {                              // Upper quadrant
        AddTwoLimbsPlusOneLimb(Nminustt, 2*t, indep);  // n - 4 * t*t + 4 * t
        ShowLabel("SW-NE: 4t<sup>2</sup> - 4t", -4, indep);
        AddTwoLimbsPlusOneLimb(Nminustt, t, indep);    // n - 4 * t*t + 2 * t
        ShowLabel("NW-SE: 4t<sup>2</sup> - 2t", -2, indep);
      }
    }
    else
    {
      if (xLogical - yLogical >= 0)
      {                              // Lower quadrant
        ShowLabel("SW-NE: 4t<sup>2</sup>", 0, Nminustt);
        AddTwoLimbsPlusOneLimb(Nminustt, -t, indep);   // n - 4 * t*t - 2 * t
        ShowLabel("NW-SE: 4t<sup>2</sup> + 2t", 2, indep);
      }
      else
      {                             // Left quadrant
        ShowLabel("SW-NE: 4t<sup>2</sup>", 0, Nminustt);
        AddTwoLimbsPlusOneLimb(Nminustt, t, indep);    // n - 4 * t*t + 2 * t
        ShowLabel("NW-SE: 4t<sup>2</sup> - 2t", -2, indep);
      }
    }
    ptrText = &infoText[strlen(infoText)];
    *ptrText++ = 'x';
    *ptrText++ = '=';
    ptrText = appendInt(ptrText, xLogical);
    *ptrText++ = ',';
    *ptrText++ = ' ';
    *ptrText++ = 'y';
    *ptrText++ = '=';
    ptrText = appendInt(ptrText, yLogical);
    *ptrText++ = ',';
    *ptrText++ = ' ';
    *ptrText++ = 'n';
    *ptrText++ = '=';
    *ptrText++ = '<';
    *ptrText++ = 'b';
    *ptrText++ = '>';
    ptrText = appendInt64(ptrText, value);
    *ptrText++ = '<';
    *ptrText++ = '/';
    *ptrText++ = 'b';
    *ptrText++ = '>';
    *ptrText++ = ',';
    *ptrText++ = ' ';
    *ptrText++ = 't';
    *ptrText++ = '=';
    ptrText = appendInt(ptrText, t/2);
  }
  getN(xCenter, yCenter, value);
  *ptrText++ = '^';       // Separator between bottom text and center text.
  ptrText = appendInt64(ptrText, value);
  *ptrText = 0;
  return infoText;
}  

EXTERNALIZE void moveSpiral(int deltaX, int deltaY)
{
  xFraction -= deltaX;
  xCenter += xFraction >> thickness;
  xFraction &= (1 << thickness) - 1;
  yFraction += deltaY;
  yCenter += yFraction >> thickness;
  yFraction &= (1 << thickness) - 1; 
}

// inputBoxNbr = 1 -> changing center
// inputBoxNbr = 2 -> changing start value
EXTERNALIZE char *nbrChanged(char *value, int inputBoxNbr, int newWidth, int newHeight)
{ 
  int temp[2];
  int nbr0, nbr1;
  int a;
  int index;
  width = newWidth;
  height = newHeight;
  nbr0 = nbr1 = 0;
  for (index=0; index<19; index++)
  {
    int charConverted;
    double dProd;
    if (*value == 0)
    {      // End of string, so end of conversion from string to number.
      break;
    }
    charConverted = (*value++ - '0');
    dProd = (double)nbr0 * 10 + charConverted;
    nbr0 = (nbr0 * 10 + charConverted) & MAX_INT_NBR;
    nbr1 = nbr1 * 10 + (int)(dProd / LIMB_RANGE);
  }
  if (inputBoxNbr == 1)
  {           // Changing center
              // nbr <- nbr - startNumber + 1 
    unsigned int carry;
    int borrow = nbr0 - startNumber[0];
    nbr0 = borrow & MAX_INT_NBR;
    nbr1 = ((borrow >> BITS_PER_GROUP) + nbr1 - startNumber[1]) & MAX_INT_NBR;           
    carry = nbr0 + 1;
    nbr0 = carry & MAX_INT_NBR;
    nbr1 = (carry >> BITS_PER_GROUP) + nbr1;
    if (nbr1 > 0 || nbr0 >= 1)
    {       // nbr >= 1
      int diff;
      a = ((int)sqrt(((double)nbr1*LIMB_RANGE + nbr0-1))+1)/2;
      diff = (nbr0 - 4*a*a) & MAX_INT_NBR;
      if (diff & HALF_INT_RANGE)
      {     // Number is negative.
        diff -= MAX_INT_NBR;
      }
      if (diff < -2*a+2)
      {
        xCenter = -3*a + 1 - diff;
        yCenter = a;
      }
      else if (diff <= 0)
      {
        xCenter = -a;
        yCenter = -a + 1 - diff;
      }
      else if (diff <= 2*a)
      {
        xCenter = diff - a - 1;
        yCenter = -a;
      }
      else
      {
        xCenter = a;
        yCenter = diff - 3*a - 1;
      }
    }
    xFraction = 0;
    yFraction = 0;
  }
  else if (inputBoxNbr == 2)
  {           // Changing start number.
    startNumber[0] = nbr0;
    startNumber[1] = nbr1;
  }
  else if (inputBoxNbr == 3)
  {           // Zoom in
    if (thickness == 5)
    {
      return NULL;
    }
    thickness++;
    xFraction <<= 1;
    yFraction <<= 1;
  }
  else
  {           // Zoom out
    if (thickness == 0)
    {
      return NULL;
    }
    thickness--;
    xFraction >>= 1;
    yFraction >>= 1;
  }
  getN(xCenter, yCenter, temp);
  infoText[0] = '^';
  appendInt64(&infoText[1], temp);
  return infoText;
}

#ifndef __EMSCRIPTEN__

void drawUlamSpiral(void)
{
  if (SDL_MUSTLOCK(doubleBuffer))
  {
    SDL_LockSurface(doubleBuffer);
  }
  drawPartialUlamSpiral(-width / 2,
                        width / 2,
                        -height / 2,
                        height / 2);
  if (SDL_MUSTLOCK(doubleBuffer))
  {
    SDL_UnlockSurface(doubleBuffer);
  }
  SDL_BlitSurface(doubleBuffer, NULL, screen, NULL);
  SDL_Flip(screen);
  oldXCenter = xCenter;
  oldYCenter = yCenter;
  oldXFraction = xFraction;
  oldYFraction = yFraction;
  return;
}

void setNewVideoMode(void)
{
  if (doubleBuffer != NULL)
  {
    SDL_FreeSurface(doubleBuffer);
  }
  screen = SDL_SetVideoMode(width, height, 32, SDL_SWSURFACE);
  doubleBuffer = SDL_CreateRGBSurface(0, width, height, 32, 0, 0, 0, 0);
}

void iteration(void)
{
  SDL_Event event;
  SDL_Rect rectSrc, rectDest;
  int xMax, yMin, yMax, yBound, xMove, yMove;
        
  while (SDL_PollEvent(&event))
  {                           // New event arrived.
    if (event.type == SDL_QUIT)
    {                         // Quit event arrived, so exit application.
      quit = 1;
    }
    if (event.type == SDL_MOUSEMOTION)
    {
      if (event.motion.state & SDL_BUTTON_LMASK)
      {                        // Drag operation.
        xFraction -= event.motion.xrel;
        xCenter += xFraction >> thickness;
        xFraction &= (1 << thickness) - 1;
        yFraction += event.motion.yrel;
        yCenter += yFraction >> thickness;
        yFraction &= (1 << thickness) - 1;
      }
      else
      {                        // Show information.
        ShowInformation(event.motion.x, event.motion.y);
      }
    }
    if (event.type == SDL_MOUSEBUTTONDOWN)
    {
      if (event.button.button == SDL_BUTTON_WHEELUP)
      {
        if (thickness < 5)
        {
          thickness++;
          xFraction <<= 1;
          yFraction <<= 1;
          drawUlamSpiral();
        }
      }
      else if (event.button.button == SDL_BUTTON_WHEELDOWN)
      {
        if (thickness > 0)
        {
          thickness--;
          xFraction >>= 1;
          yFraction >>= 1;
          drawUlamSpiral();
        }
      }
      else
      {
        ShowInformation(event.button.x, event.button.y);
      }
    }
  }
  if (++timer == 6)
  {
    timer = 0;
    if (oldXCenter != xCenter || oldYCenter != yCenter ||
        oldXFraction != xFraction || oldYFraction != yFraction)
    {
      int xMin;
          // Move pixels of double buffer according to drag direction.
      xMove = (xCenter << thickness) + xFraction - (oldXCenter << thickness) - oldXFraction;
      yMove = (yCenter << thickness) + yFraction - (oldYCenter << thickness) - oldYFraction;
      if (xMove > 0)
      {           // Move pixels to left.
        rectDest.x = 0;
        rectDest.w = width - xMove;
        rectSrc.x = xMove;
        rectSrc.w = rectDest.w;
      }
      else
      {           // Move pixels to right.
        rectDest.x = -xMove;
        rectDest.w = width + xMove;
        rectSrc.x = 0;
        rectSrc.w = rectDest.w;
      }
      if (yMove > 0)
      {           // Move pixels up.
        rectDest.y = yMove;
        rectDest.h = height - yMove;
        rectSrc.y = 0;
        rectSrc.h = rectDest.h;
      }
      else
      {           // Move pixels down.
        rectDest.y = 0;
        rectDest.h = height + yMove;
        rectSrc.y = -yMove;
        rectSrc.h = rectDest.h;
      }
      SDL_BlitSurface(doubleBuffer, &rectSrc, doubleBuffer, &rectDest);
      if (SDL_MUSTLOCK(doubleBuffer))
      {
        SDL_LockSurface(doubleBuffer);
      }
      xMin = -width / 2;
      xMax = width / 2;
      yMin = -height / 2,
      yMax = height / 2;
      if (yMove > 0)
      {              // Move pixels up.
        yBound = yMax - yMove;
                     // Draw bottom rectangle.
        drawPartialUlamSpiral(xMin, xMax, yBound, yMax);
        if (xMove > 0)
        {            // Move pixels left.
                     // Draw right rectangle.
          drawPartialUlamSpiral(xMax - xMove, xMax, yMin, yBound);
        }
        if (xMove < 0)
        {            // Move pixels right.
                     // Draw left rectangle.
          drawPartialUlamSpiral(xMin, xMin - xMove, yMin, yBound);
        }
      }
      else if (yMove < 0)
      {              // Move pixels up.
        yBound = yMin - yMove;
                     // Draw bottom rectangle.
        drawPartialUlamSpiral(xMin, xMax, yMin, yBound);
        if (xMove > 0)
        {            // Move pixels left.
                     // Draw right rectangle.
          drawPartialUlamSpiral(xMax - xMove, xMax, yBound, yMax);
        }
        if (xMove < 0)
        {            // Move pixels right.
                     // Draw left rectangle.
          drawPartialUlamSpiral(xMin, xMin - xMove, yBound, yMax);
        }
      }
      else
      {
        if (xMove > 0)
        {            // Move pixels left.
                     // Draw right rectangle.
          drawPartialUlamSpiral(xMax - xMove, xMax, yMin, yMax);
        }
        if (xMove < 0)
        {            // Move pixels right.
                     // Draw left rectangle.
          drawPartialUlamSpiral(xMin, xMin - xMove, yMin, yMax);
        }
      }
      if (SDL_MUSTLOCK(doubleBuffer))
      {
        SDL_UnlockSurface(doubleBuffer);
      }
      SDL_BlitSurface(doubleBuffer, NULL, screen, NULL);
      SDL_Flip(screen);
      oldXCenter = xCenter;
      oldYCenter = yCenter;
      oldXFraction = xFraction;
      oldYFraction = yFraction;
    }
  }  
}

int main(int argc, char *argv[])
{
  width = 256;
  height = 256;
  initMultipleArray();
  SDL_Init(SDL_INIT_VIDEO);
  setNewVideoMode();
  drawUlamSpiral();
  while (!quit)
  {
    iteration();
  }
  SDL_Quit();
  return 0;
}
#endif
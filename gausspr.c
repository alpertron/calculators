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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>

#define SMALL_NUMBER_BOUND 32768
#ifdef __EMSCRIPTEN__
  #define MAX_WIDTH 2048
  #define pixelXY(x, y) &pixels[y * MAX_WIDTH + x];
  unsigned int pixels[2048*MAX_WIDTH];
#else
  #include <SDL.h>
  #define pixels (unsigned int *)(4*MAX_WIDTH*32)
  SDL_Surface *screen;
  SDL_Surface *doubleBuffer;
  int oldXCenter;
  int oldYCenter;
  int oldXFraction;
  int oldYFraction;
  int timer;
  int quit;
#endif

char *appendInt(char *text, int value);
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
int xCenter;
int xFraction;                      // Range of fraction: 0 to (1 << thickness) - 1
int yCenter;
int yFraction;
int width;
int height;
char initMultipleArrayCalled;
static int multiple[25][97];
static char primes[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41,
                         43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97 };

static void initMultipleArray(void)
{
  int i;
  int j;
  for (i = 0; i<25; i++)
  {
    int k = primes[i];
    for (j = (k / 2) + 1; j >= 0; j--)
    {
      multiple[i][j*j%k] = 1;
    }
  }
  initMultipleArrayCalled = 1;
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
    dAccumulator = *(Nbr+i) - (TestNbr[i] * dTrialQuotient) + carry + dDelta;
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
  MontDig = ((uint32_t)Pr * MontgomeryMultN) & MAX_INT_NBR;
  Pr = ((((uint64_t)MontDig * TestNbr0) + Pr) >> BITS_PER_GROUP) +
    (uint64_t)MontDig * TestNbr1 + (uint64_t)Nbr * factor2_1;
  Prod0 = Pr & MAX_INT_NBR;
  Prod1 = (uint32_t)(Pr >> BITS_PER_GROUP);
    
  Nbr = *(factor1 + 1);
  Pr = Nbr * (uint64_t)factor2_0 + (uint32_t)Prod0;
  MontDig = ((uint32_t)Pr * MontgomeryMultN) & MAX_INT_NBR;
  Pr = ((((uint64_t)MontDig * TestNbr0) + Pr) >> BITS_PER_GROUP) +
    (uint64_t)MontDig * TestNbr1 + (uint64_t)Nbr * factor2_1 + (uint32_t)Prod1;
  Prod0 = Pr & MAX_INT_NBR;
  Prod1 = (uint32_t)(Pr >> BITS_PER_GROUP);
    
  if ((Pr >= ((uint64_t)(TestNbr1 + 1) << BITS_PER_GROUP)) ||
     ((Prod1 == (uint32_t)TestNbr1) && (Prod0 >= (uint32_t)TestNbr0)))
  {
    int32_t borrow = (int32_t)Prod0 - (int32_t)TestNbr0;
    Prod0 = borrow & MAX_INT_NBR;
    Prod1 = ((borrow >> BITS_PER_GROUP) + (int32_t)Prod1 - (int32_t)TestNbr1) & MAX_INT_NBR;
  }
#else
  double dInvLimbRange = 1.0 / (double)LIMB_RANGE;
  int Nbr = *(factor1);
  double dNbr = (double)Nbr;
  int low = Nbr * factor2_0;
  double dAccum = dNbr * (double)factor2_0;
  int MontDig = (low * MontgomeryMultN) & MAX_VALUE_LIMB;
  double dMontDig = (double)MontDig;
  dAccum += dMontDig * (double)TestNbr0;
  // At this moment dAccum is multiple of LIMB_RANGE.
  dAccum = floor((dAccum*dInvLimbRange) + 0.5);
  low = ((unsigned int)dAccum + (MontDig * TestNbr1) +
               (Nbr * factor2_1)) & MAX_VALUE_LIMB;
  dAccum += (dMontDig * TestNbr1) + (dNbr * factor2_1);
  Prod0 = low;
  if (low < HALF_INT_RANGE)
  {
    dAccum = ((dAccum + (double)(HALF_INT_RANGE / 2U))*dInvLimbRange);
  }
  else
  {
    dAccum = ((dAccum - (double)(HALF_INT_RANGE / 2U))*dInvLimbRange);
  }
  Prod1 = (unsigned int)dAccum;  // Most significant limb can be greater than LIMB_RANGE
  
  Nbr = *(factor1 + 1);
  dNbr = (double)Nbr;
  low = (Nbr * factor2_0) + Prod0;
  dAccum = (dNbr * (double)factor2_0) + (double)Prod0;
  MontDig = (low * MontgomeryMultN) & MAX_VALUE_LIMB;
  dMontDig = (double)MontDig;
  dAccum += dMontDig * (double)TestNbr0;
  // At this moment dAccum is multiple of LIMB_RANGE.
  dAccum = floor((dAccum*dInvLimbRange) + 0.5);
  low = ((unsigned int)dAccum + (MontDig * TestNbr1) +
               (Nbr * factor2_1) + Prod1) & MAX_VALUE_LIMB;
  dAccum += (dMontDig * TestNbr1) + (dNbr * factor2_1) + (unsigned int)Prod1;
  Prod0 = low;
  if (low < HALF_INT_RANGE)
  {
    dAccum = ((dAccum + (double)(HALF_INT_RANGE / 2U))*dInvLimbRange);
  }
  else
  {
    dAccum = ((dAccum - (double)(HALF_INT_RANGE / 2U))*dInvLimbRange);
  }
  Prod1 = (unsigned int)dAccum;  // Most significant limb can be greater than LIMB_RANGE
  
  if (((unsigned int)Prod1 > (unsigned int)TestNbr1) ||
       (((unsigned int)Prod1 == (unsigned int)TestNbr1) && ((unsigned int)Prod0 >= (unsigned int)TestNbr0)))
  {        // Prod >= TestNbr, so perform Prod <- Prod - TestNbr
    carry = Prod0 - TestNbr0;
    Prod0 = carry & MAX_VALUE_LIMB;
    Prod1 = ((carry >> BITS_PER_GROUP) + Prod1 - TestNbr1) & MAX_VALUE_LIMB;
  }
#endif  
  *Product = Prod0;
  *(Product+1) = Prod1;
}

static void AddBigNbr(int *Nbr1, int *Nbr2, int *Sum)
{
  unsigned int carry = *(Nbr1) + *(Nbr2);
  *(Sum) = carry & MAX_INT_NBR;
  carry = (carry >> BITS_PER_GROUP) + *(Nbr1+1) + *(Nbr2+1);
  *(Sum+1) = carry & MAX_INT_NBR;
}

static void AddBigNbrModN(int *Nbr1, int *Nbr2, int *Sum)
{
  int Sum0;
  int Sum1;
  int TestNbr0 = TestNbr[0];
  int TestNbr1 = TestNbr[1];
  unsigned int carry = *Nbr1 + *Nbr2;
  Sum0 = carry & MAX_INT_NBR;
  carry = (carry >> BITS_PER_GROUP) + *(Nbr1 + 1) + *(Nbr2 + 1);
  Sum1 = carry & MAX_INT_NBR;
  if ((carry > (unsigned int)TestNbr1) || 
     ((carry == (unsigned int)TestNbr1) && (Sum0 >= TestNbr0)))
  {
    int borrow = Sum0 - TestNbr0;
    Sum0 = borrow & MAX_INT_NBR;
    Sum1 = ((borrow >> BITS_PER_GROUP) + Sum1 - TestNbr1) & MAX_INT_NBR;
  }
  *Sum = Sum0;
  *(Sum+1) = Sum1;
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
static bool isPrime(int *value)
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
    LIMIT(3000000000000000000ll)   // Greater than any argument of isPrime()
  };
  int i;
  int j;
  int index, indexLSB, indexMSB, idxNbrMSB;
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
    else if (TestNbr0 == 2)
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
      if (TestNbr0 % base == 0)
      {
        return false;          // Number is multiple of base, so it is composite.
      }
    }
    // Check whether TestNbr is multiple of base. In this case the number would be composite.
    else if (((TestNbr1 % base) * (LIMB_RANGE % base) + TestNbr0 % base) % base == 0)  // No overflow possible.
    {
      return false;            // Number is multiple of base, so it is composite.
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
  
  for (i=0, j=0; (limits[j+1] < TestNbr1) || ((limits[j+1] == TestNbr1) && (limits[j] < TestNbr0)); i++, j+=2)
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

static void multiply(int factor1, int factor2, int *prod)
{
#ifdef _USING64BITS_
  int64_t prod64 = (int64_t)factor1 * (int64_t)factor2;
  // Obtain most significant limb.
  *(prod+1) = (int)(prod64 >> BITS_PER_GROUP);
  // Obtain least significant limb.
  *prod = (int)prod64 & MAX_VALUE_LIMB;
#else
  int low = (factor1 * factor2) & MAX_VALUE_LIMB;
  double dAccum = (double)factor1 * (double)factor2;
  *prod = low;
  if (low < HALF_INT_RANGE)
  {
    dAccum = ((dAccum + (double)(HALF_INT_RANGE / 2U)) / LIMB_RANGE);
  }
  else
  {
    dAccum = ((dAccum - (double)(HALF_INT_RANGE / 2U)) / LIMB_RANGE);
  }
  *(prod+1) = (unsigned int)dAccum;
#endif 
}

 // A Gaussian number x+iy is prime if one of these 3 conditions is true:
 // 1) x is not zero and y is not zero and x^2+y^2 is prime
 // 2) x=0, y is prime and y=3 (mod 4)
 // 3) y=0, x is prime and x=3 (mod 4)
static void setPoint(int x, int y)
{
  int xPhysical;
  int yPhysical;
  int row;
  int col;
  int firstRow;
  int firstCol;
  int lastRow;
  int lastCol;
  int firstRow2;
  int lastRow2;
  int firstCol2;
  int lastCol2;
  int xSquared[NBR_LIMBS];
  int ySquared[NBR_LIMBS];
  int check[NBR_LIMBS];
  unsigned int *ptrPixel;
  unsigned int color;
#ifdef __EMSCRIPTEN__  
  unsigned int colorBlack = 0xFF000000U;
  unsigned int colorGreen = 0xFF00C000U;
  unsigned int colorWhite = 0xFFC0C0C0U;
#else
  Uint32 colorBlack = SDL_MapRGBA(doubleBuffer->format, 0, 0, 0, 255);
  Uint32 colorGreen = SDL_MapRGBA(doubleBuffer->format, 0, 192, 0, 255);
  Uint32 colorWhite = SDL_MapRGBA(doubleBuffer->format, 192, 192, 192, 255);
#endif
  xPhysical = width / 2 + ((x - xCenter) << thickness) - xFraction;
  yPhysical = height / 2 - ((y - yCenter) << thickness) + yFraction;
  color = colorBlack;          // Indicate not prime in advance.
  if (x == 0)
  {   // Number is imaginary.
    check[0] = (y >= 0 ? y : -y);
    check[1] = 0;
    if (((check[0] & 0x03) == 3) && isPrime(check))
    {                          // Number is gaussian prime.
      color = colorGreen;
    }
  }
  else if (y == 0)
  {   // Number is real.
    check[0] = (x >= 0 ? x : -x);
    check[1] = 0;
    if (((check[0] & 0x03) == 3) && isPrime(check))
    {                          // Number is gaussian prime.
      color = colorGreen;
    }
  }
  else
  {   // Number is not real or imaginary.
    multiply(x, x, xSquared);
    multiply(y, y, ySquared);
    AddBigNbr(xSquared, ySquared, check);   // Check x^2 + y^2.
    if (isPrime(check))
    {                          // Number is gaussian prime.
      color = colorGreen;
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
    if (x == 0)
    {                     // Draw Y axis if possible.
      col = xPhysical + (1 << (thickness - 1));
      if ((col >= 0) && (col < width))
      {
        ptrPixel = pixelXY(col, firstRow);
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
    else if (x % 10 == 0)
    {
      col = xPhysical + (1 << (thickness - 1));
      if ((col >= 0) && (col < width))
      {
        firstRow2 = yPhysical + (1 << (thickness-2));
        lastRow2 = firstRow2 + (1 << (thickness-1));
        firstRow2 = (firstRow2<0 ? 0 : firstRow2);
        if (lastRow2 > height)
        {
          lastRow2 = height;
        }
        ptrPixel = pixelXY(col, firstRow2);
        for (row = firstRow2; row < lastRow2; row++)
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
    if (y == 0)
    {                     // Draw X axis if possible.
      row = yPhysical + (1 << (thickness - 1));
      if ((row >= 0) && (row < height))
      {
        ptrPixel = pixelXY(firstCol, row);
        for (col = firstCol; col < lastCol; col++)
        {
          *ptrPixel = color;
          ptrPixel++;
        }
      }
    }
    else if (y % 10 == 0)
    {
      row = yPhysical + (1 << (thickness - 1));
      if ((row >= 0) && (row < height))
      {
        firstCol2 = xPhysical + (1 << (thickness - 2));
        lastCol2 = firstCol2 + (1 << (thickness - 1));
        firstCol2 = (firstCol2<0 ? 0 : firstCol2);
        if (lastCol2 > width)
        {
          lastCol2 = width;
        }
        ptrPixel = pixelXY(firstCol2, row);
        for (col = firstCol2; col < lastCol2; col++)
        {
          *ptrPixel = color;
          ptrPixel++;
        }
      }
    }
  }
}     /* end method setPoint */

void drawPartialGraphic(int xminDisp, int xmaxDisp, int yminDisp, int ymaxDisp)
{
  int x;
  int y;
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
}      /* end method drawPartialGraphic */

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
    if ((quot != 0) || zeroIsSignificant)
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

char *getInformation(int x, int y)
{
  int yLogical;
  char *ptrText = infoText;
  
  infoText[0] = 0;   // Empty string.
  if (x >= 0)
  {
    int xLogical = xCenter + ((xFraction + x - width / 2) >> thickness);
    yLogical = yCenter + 1 + ((yFraction - y + height / 2) >> thickness);
    ptrText = appendInt(infoText, xLogical);
    if (yLogical >= 0)
    {
      *ptrText++ = ' ';
      *ptrText++ = '+';
      *ptrText++ = ' ';
      ptrText = appendInt(ptrText, yLogical);
    }
    else
    {
      *ptrText++ = ' ';
      *ptrText++ = '-';
      *ptrText++ = ' ';
      ptrText = appendInt(ptrText, -yLogical);
    }
    *ptrText++ = 'i';
  }
  *ptrText++ = '^';       // Append separator.
  ptrText = appendInt(ptrText, xCenter);
  *ptrText++ = '^';       // Append separator.
  ptrText = appendInt(ptrText, yCenter);
  *ptrText = 0;
  return infoText;
}

void moveGraphic(int deltaX, int deltaY)
{
  xFraction -= deltaX;
  xCenter += xFraction >> thickness;
  xFraction &= (1 << thickness) - 1;
  yFraction += deltaY;
  yCenter += yFraction >> thickness;
  yFraction &= (1 << thickness) - 1; 
}

static int getValue(char *value)
{
  int nbr = 0;
  int sign = 0;
  if (*value == '-')
  {
    sign = 1;
	value++;
  }
  while (*value)
  {
    nbr = nbr*10 + (*value++ - '0');
  }
  return sign? -nbr: nbr;
}

int nbrChanged(char *value, int inputBoxNbr, int newWidth, int newHeight)
{
  width = newWidth;
  height = newHeight;
  if (inputBoxNbr == 1)
  {           
    xCenter = getValue(value);  // Changing center X.
    xFraction = 0;
    value += strlen(value) + 1;
    yCenter = getValue(value);  // Changing center Y.
    yFraction = 0;
  }
  else if (inputBoxNbr == 3)
  {           // Zoom in
    if (thickness == 5)
    {
      return 0;
    }
    thickness++;
    xFraction <<= 1;
    yFraction <<= 1;
  }
  else
  {           // Zoom out
    if (thickness == 0)
    {
      return 0;
    }
    thickness--;
    xFraction >>= 1;
    yFraction >>= 1;
  }
  return 0;
}

#ifdef __EMSCRIPTEN__
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
#else
static void drawGraphic(void)
{
  if (SDL_MUSTLOCK(doubleBuffer))
  {
    SDL_LockSurface(doubleBuffer);
  }
  drawPartialGraphic(-width / 2,
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
  return;
}

static void setNewVideoMode(void)
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
  SDL_Rect rectSrc;
  SDL_Rect rectDest;
  int xMax;
  int yMin;
  int yMax;
  int yBound;
  int xMove;
  int yMove;
        
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
          drawGraphic();
        }
      }
      else if (event.button.button == SDL_BUTTON_WHEELDOWN)
      {
        if (thickness > 0)
        {
          thickness--;
          xFraction >>= 1;
          yFraction >>= 1;
          drawGraphic();
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
    if ((oldXCenter != xCenter) || (oldYCenter != yCenter) ||
        (oldXFraction != xFraction) || (oldYFraction != yFraction))
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
        drawPartialGraphic(xMin, xMax, yBound, yMax);
        if (xMove > 0)
        {            // Move pixels left.
                     // Draw right rectangle.
          drawPartialGraphic(xMax - xMove, xMax, yMin, yBound);
        }
        if (xMove < 0)
        {            // Move pixels right.
                     // Draw left rectangle.
          drawPartialGraphic(xMin, xMin - xMove, yMin, yBound);
        }
      }
      else if (yMove < 0)
      {              // Move pixels up.
        yBound = yMin - yMove;
                     // Draw bottom rectangle.
        drawPartialGraphic(xMin, xMax, yMin, yBound);
        if (xMove > 0)
        {            // Move pixels left.
                     // Draw right rectangle.
          drawPartialGraphic(xMax - xMove, xMax, yBound, yMax);
        }
        if (xMove < 0)
        {            // Move pixels right.
                     // Draw left rectangle.
          drawPartialGraphic(xMin, xMin - xMove, yBound, yMax);
        }
      }
      else
      {
        if (xMove > 0)
        {            // Move pixels left.
                     // Draw right rectangle.
          drawPartialGraphic(xMax - xMove, xMax, yMin, yMax);
        }
        if (xMove < 0)
        {            // Move pixels right.
                     // Draw left rectangle.
          drawPartialGraphic(xMin, xMin - xMove, yMin, yMax);
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
  SDL_Init(SDL_INIT_VIDEO);
  setNewVideoMode();
  drawGraphic();
  while (!quit)
  {
    iteration();
  }
  SDL_Quit();
  return 0;
}
#endif
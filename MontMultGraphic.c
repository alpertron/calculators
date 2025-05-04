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

#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include "isprime.h"
#define SMALL_NUMBER_BOUND 32768

extern int TestNbrGraphic[NBR_LIMBS];
extern uint32_t MontMultGraphicN;

/*
 For large numbers the REDC algorithm is:
 m <- ((T mod R)N') mod R
 t <- (T - mN) / R
 if t < 0 then
   return t + N
 else
   return t
 end if
*/

void MontMultGraphic(const int *factor1, const int *factor2, int *Product)
{
  int TestNbr0 = TestNbrGraphic[0];
  if (TestNbrGraphic[1] == 0)
  {   // Only one limb.
#ifdef __ARM_ARCH_7A__
    int Nbr = *factor1;
    int64_t Pr = ((int64_t)Nbr * (int64_t)*factor2);
    int32_t MontDig = ((int32_t)Pr * MontMultGraphicN) & MAX_VALUE_LIMB;
    Pr -= (int64_t)MontDig * (int64_t)TestNbr0;
    *Product = (int32_t)(Pr >> BITS_PER_GROUP);
    if (*Product < 0)
    {
      *Product += TestNbrGraphic[0];
    }
#else
    int mod = TestNbr0;
    int multiplicand = *factor1;
    int multiplier = *factor2;
    if (mod < SMALL_NUMBER_BOUND)
    {
      *Product = multiplicand * multiplier % mod;
    }
    else
    {   // TestNbr has one limb but it is not small.
#ifdef _USING64BITS_
      *Product = (int64_t)*factor1 * *factor2 % mod;
#else
      // Round up quotient.
      int quotient = (int)floor((double)multiplicand * (double)multiplier / (double)mod + 0.5);
      int remainder = (multiplicand * multiplier) - (quotient * mod);
      if (remainder < 0)
      {    // Quotient was 1 more than expected. Adjust remainder.
        remainder += mod;
      }
      *Product = remainder;
#endif
    }
    *(Product + 1) = 0;
    return;
#endif
  }
  // Two limbs.
  int TestNbr1 = TestNbrGraphic[1];
  int32_t Prod0;
  int32_t Prod1;
  int factor2_0 = *factor2;
  int factor2_1 = *(factor2+1);
#ifdef _USING64BITS_
  int64_t Pr;
  int Nbr;
  int MontDig;
  
  Nbr = *factor1;
  Pr = Nbr * (int64_t)factor2_0;
  MontDig = ((int32_t)Pr * MontMultGraphicN) & MAX_VALUE_LIMB;
  Pr = ((Pr - ((int64_t)MontDig * (int64_t)TestNbr0)) >> BITS_PER_GROUP) -
    ((int64_t)MontDig * (int64_t)TestNbr1) +
    ((int64_t)Nbr * (int64_t)factor2_1);
  Prod0 = (int32_t)(Pr & MAX_VALUE_LIMB);
  Prod1 = (int32_t)(Pr >> BITS_PER_GROUP);
    
  Nbr = *(factor1 + 1);
  Pr = Nbr * (int64_t)factor2_0 + (int64_t)Prod0;
  MontDig = ((int32_t)Pr * MontMultGraphicN) & MAX_VALUE_LIMB;
  Pr = ((Pr - ((int64_t)MontDig * (int64_t)TestNbr0)) >> BITS_PER_GROUP) -
    ((int64_t)MontDig * (int64_t)TestNbr1) +
    ((int64_t)Nbr * (int64_t)factor2_1) + Prod1;
  Prod0 = (int32_t)(Pr & MAX_VALUE_LIMB);
  Prod1 = (int32_t)(Pr >> BITS_PER_GROUP);
    
#else
  double dInvLimbRange = 1.0 / (double)LIMB_RANGE;
  int Nbr = *factor1;
  double dNbr = (double)Nbr;
  int low = Nbr * factor2_0;
  double dAccum = dNbr * (double)factor2_0;
  int MontDig = (low * MontMultGraphicN) & MAX_INT_NBR;
  double dMontDig = (double)MontDig;
  dAccum -= dMontDig * (double)TestNbr0;
  // At this moment dAccum is multiple of LIMB_RANGE.
  dAccum = floor((dAccum*dInvLimbRange) + 0.5);
  low = ((int)dAccum - (MontDig * TestNbr1) + (Nbr * factor2_1)) &
        MAX_INT_NBR;
  dAccum += (dNbr * (double)factor2_1) - (dMontDig * (double)TestNbr1);
  Prod0 = low;
  // Casting from double to int truncates to nearest to zero,
  // so the conversion is done on non-negative numbers.
  if (low < HALF_INT_RANGE)
  {
    dAccum = (dAccum * dInvLimbRange) + 2147483648.25;
  }
  else
  {
    dAccum = (dAccum * dInvLimbRange) + 2147483647.75;
  }
  Prod1 = (unsigned int)dAccum ^ LIMB_RANGE;
  
  Nbr = *(factor1 + 1);
  dNbr = (double)Nbr;
  low = (Nbr * factor2_0) + Prod0;
  dAccum = (dNbr * (double)factor2_0) + (double)Prod0;
  MontDig = (low * MontMultGraphicN) & MAX_VALUE_LIMB;
  dMontDig = (double)MontDig;
  dAccum -= dMontDig * (double)TestNbr0;
  // At this moment dAccum is multiple of LIMB_RANGE.
  dAccum = floor((dAccum*dInvLimbRange) + 0.5);
  low = ((int)dAccum - (MontDig * TestNbr1) + (Nbr * factor2_1) + Prod1) &
        MAX_INT_NBR;
  dAccum += (dNbr * (double)factor2_1) - (dMontDig * (double)TestNbr1) +
            (double)Prod1;
  Prod0 = low;
  // Casting from double to int truncates to nearest to zero,
  // so the conversion is done on non-negative numbers.
  if (low < HALF_INT_RANGE)
  {
    dAccum = (dAccum * dInvLimbRange) + 2147483648.25;
  }
  else
  {
    dAccum = (dAccum * dInvLimbRange) + 2147483647.75;
  }
  Prod1 = (unsigned int)dAccum ^ LIMB_RANGE;
#endif  
  if (Prod1 < 0)
  {
    uint32_t carry = Prod0 + TestNbr0;
    Prod0 = carry & MAX_VALUE_LIMB;
    Prod1 += TestNbr1 + (carry >> BITS_PER_GROUP);
  }
  *Product = Prod0;
  *(Product+1) = Prod1;
}


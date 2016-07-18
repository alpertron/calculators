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
#include "bignbr.h"

#if BITS_PER_GROUP == 15
#define DIGITS_PER_LIMB 4
#define MAX_LIMB_CONVERSION 10000
#else
#define DIGITS_PER_LIMB 9
#define MAX_LIMB_CONVERSION 1000000000
#endif

extern int lang;
static limb power10000[MAX_LEN];
static limb temp[MAX_LEN];
static void add(limb *addend1, limb *addend2, limb *sum, int length);

  // Convert number to little-endian.
void Dec2Bin(char *decimal, limb *binary, int digits, int *bitGroups)
{
  // First step: separate in groups of DIGITS_PER_LIMB digits.
  char *ptrSrc;
  limb *ptrDest, *ptrBinary;
  int digit, groupNbr, multiplier;
  int outerGroup, innerGroup;
  int nbrGroups = 1;
  for (nbrGroups = 1; nbrGroups*DIGITS_PER_LIMB < digits; nbrGroups *= 2)
  {
  }
  memset(binary, 0, nbrGroups * sizeof(limb));
  memset(power10000, 0, nbrGroups * sizeof(limb));
  power10000[0].x = MAX_LIMB_CONVERSION;
  ptrDest = binary;
  for (ptrSrc = decimal + digits - 1; ptrSrc >= decimal + DIGITS_PER_LIMB-1; ptrSrc -= DIGITS_PER_LIMB)
  {
    int limbContents = 0;
    for (digit = DIGITS_PER_LIMB-1; digit >= 0; digit--)
    {
      limbContents = (limbContents * 10) + *(ptrSrc - digit) - '0';
    }
    (ptrDest++)->x = limbContents;
  }
  digit = 0;
  multiplier = 1;
  for (; ptrSrc >= decimal; ptrSrc--)
  {
    digit += multiplier * (*ptrSrc - '0');
    multiplier *= 10;
  }
  (ptrDest++)->x = digit;
  for (outerGroup = 1; outerGroup < nbrGroups; outerGroup += outerGroup)
  {
    for (innerGroup = 0; innerGroup < nbrGroups; innerGroup += 2*outerGroup)
    {
      ptrBinary = binary + innerGroup;
      multiply(power10000, ptrBinary + outerGroup, temp, outerGroup, NULL);
      memset(ptrBinary + outerGroup, 0, outerGroup*sizeof(limb));
      add(temp, ptrBinary, ptrBinary, 2*outerGroup);
    }
    if (outerGroup * 2 < nbrGroups)
    {    // Square power10000.
      multiply(power10000, power10000, temp, outerGroup, NULL);
      memcpy(power10000, temp, (outerGroup * 2)*sizeof(limb));
    }
  }
  // Determine first non-significant group.
  *bitGroups = 1;    // Initialize number of groups in advance.
  for (groupNbr = nbrGroups-1; groupNbr > 0; groupNbr--)
  {
    if ((binary + groupNbr)->x != 0)
    {            // First non-significant group found.
      *bitGroups = groupNbr + 1;
      break;
    }
  }
}

void int2dec(char **pOutput, int nbr)
{
  char *ptrOutput = *pOutput;
  int significantZero = 0;
  int div = 1000000000;
  while (div > 0)
  {
    int digit;

    digit = nbr/div;
    if (digit > 0 || significantZero != 0)
    {
      significantZero = 1;
      *ptrOutput++ = (char)(digit + (int)'0');
    }
    nbr %= div;
    div /= 10;
  }
  if (significantZero == 0)
  {
    *ptrOutput++ = '0';
  }
  *pOutput = ptrOutput;
}

  // Convert little-endian number to a string with space every groupLen digits.
  // In order to perform a faster conversion, use groups of DIGITS_PER_LIMB digits.
void Bin2Dec(limb *binary, char *decimal, int nbrLimbs, int groupLen)
{
  int len, index, index2, count;
  limb carry, carry2, val;
  limb *ptrSrc = binary + nbrLimbs - 1;
  limb *ptrPower;
  char *ptrDest;
  int significantZero = 0;
  int groupCtr, digit[DIGITS_PER_LIMB];
  int digits=0;
  power10000[0].x = ptrSrc->x % MAX_LIMB_CONVERSION;
  power10000[1].x = ptrSrc->x / MAX_LIMB_CONVERSION;
     // Initialize array length.
  if (power10000[1].x == 0)
  {
    len = 1;
  }
  else
  {
    len = 2;
  }
  for (index = nbrLimbs - 2; index >= 0; index--)
  {
    // Multiply by 32768: first multiply by 128 and then by 256.
    ptrSrc--;
    carry.x = ptrSrc->x >> 8;
    carry2.x = ptrSrc->x & 0xFF;
    ptrPower = power10000;
    power10000[len++].x = 0;
    for (index2 = len; index2 > 0; index2--)
    {
      val.x = ptrPower->x;
      carry.x += val.x << (BITS_PER_GROUP/2);
      val.x = carry.x % MAX_LIMB_CONVERSION;
      carry.x /= MAX_LIMB_CONVERSION;
      carry2.x += val.x << (BITS_PER_GROUP - BITS_PER_GROUP / 2);
      (ptrPower++)->x = carry2.x % MAX_LIMB_CONVERSION;
      carry2.x /= MAX_LIMB_CONVERSION;
    }
    if (carry2.x != 0)
    {
      ptrPower->x = carry2.x;
      len++;
    }
  }
  // At this moment the array power10000 has the representation
  // of the number in base 10000 in little-endian. Convert to
  // ASCII separating every groupLen characters.
  ptrDest = decimal;
  ptrSrc = &power10000[len-1];
  groupCtr = len * DIGITS_PER_LIMB % groupLen;
  if (groupCtr == 0)
  {
    groupCtr = groupLen;
  }
  for (index = len; index > 0; index--)
  {
    int value = (int)(ptrSrc--)->x;
    for (count = 0; count < DIGITS_PER_LIMB; count++)
    {
      digit[count] = value % 10;
      value /= 10;
    }
    for (count = DIGITS_PER_LIMB-1; count >= 0; count--)
    {
      if (digit[count] != 0 || significantZero != 0)
      {
        digits++;
        *ptrDest++ = (char)(digit[count] + (int)'0');
        if (groupCtr == 1)
        {
          *ptrDest++ = ' ';
        }
        significantZero = 1;
      }
      if (--groupCtr == 0)
      {
        groupCtr = groupLen;
      }
    }
  }
  if (significantZero == 0)
  {     // Number is zero.
    *ptrDest++ = '0';
    *ptrDest = '\0';
    return;
  }
  if (digits > 30)
  {
    *ptrDest++ = '(';
    int2dec(&ptrDest, digits);
    strcpy(ptrDest, (lang==0?" digits)":" dígitos)"));
    ptrDest += strlen(ptrDest);
  }
  else if (ptrDest > decimal)
  {
    *(ptrDest-1) = '\0';       // Add terminator.
  }
}

static void add(limb *addend1, limb *addend2, limb *sum, int length)
{
  limb carry;
  int i;
  carry.x = 0;
  for (i = 0; i < length; i++)
  {
    carry.x += (addend1++)->x + (addend2++)->x;
    if (carry.x >= LIMB_RANGE)
    {
      (sum++)->x = carry.x - LIMB_RANGE;
      carry.x = 1;
    }
    else
    {
      (sum++)->x = carry.x;
      carry.x = 0;
    }
  }
  return;
}

void BigInteger2Dec(BigInteger *pBigInt, char *decimal, int groupLen)
{
  if (pBigInt->sign == SIGN_NEGATIVE)
  {
    *decimal++ = '-';
  }
  Bin2Dec(pBigInt->limbs, decimal, pBigInt->nbrLimbs, groupLen);
}
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
#include <stdio.h>
#include <math.h>

#define DIGITS_PER_LIMB 9
#define MAX_LIMB_CONVERSION 1000000000
#define FIRST_MULT  (1 << (BITS_PER_GROUP/2))
#define SECOND_MULT (LIMB_RANGE / FIRST_MULT)

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
  unsigned int div = 1000000000;
  unsigned int value = (unsigned int)nbr;
  while (div > 0)
  {
    int digit;

    digit = value/div;
    if (digit > 0 || significantZero != 0)
    {
      significantZero = 1;
      *ptrOutput++ = (char)(digit + (int)'0');
    }
    value %= div;
    div /= 10;
  }
  if (significantZero == 0)
  {
    *ptrOutput++ = '0';
  }
  *pOutput = ptrOutput;
}

// Convert little-endian number to a string with space every groupLen digits.
void Bin2Hex(limb *binary, char *decimal, int nbrLimbs, int groupLength)
{
  int showDigitsText = TRUE;
  int nbrBits, nbrHexDigits, mask, value;
  int currentGroupDigit = -1;

  if (groupLength <= 0)
  {
    groupLength = -groupLength;
    showDigitsText = FALSE;
  }
  strcpy(decimal, "<span class=\"hex\">");
  decimal += strlen(decimal);
  nbrBits = nbrLimbs * BITS_PER_GROUP;
  mask = LIMB_RANGE / 2;
  value = (binary + nbrLimbs - 1)->x;
  if (value == 0)
  {
    *decimal++ = '0';
  }
  else
  {
    while ((value & mask) == 0)
    {   // Loop that finds the most significant bit.
      mask >>= 1;
      nbrBits--;
    }
    nbrHexDigits = (nbrBits + 3) >> 2;
    if (groupLength > 0)
    {
      currentGroupDigit = nbrHexDigits % groupLength;
      if (currentGroupDigit == 0)
      {
        currentGroupDigit = groupLength;
      }
    }
    do
    {  // Get 4 bits.
      int digit = 0;
      do
      {
        digit *= 2;
        if (value & mask)
        {
          digit++;
        }
        mask >>= 1;
        if (mask == 0)
        {
          nbrLimbs--;
          value = (binary + nbrLimbs - 1)->x;
          mask = LIMB_RANGE / 2;
        }
      } while (--nbrBits & 3);
      if (digit < 10)
      {        // Convert 0 - 9 to '0' - '9'.
        *decimal++ = (char)(digit + '0');
      }
      else
      {        // Convert 10 - 15 to 'A' - 'F'.
        *decimal++ = (char)(digit + 'A' - 10);
      }
      if (--currentGroupDigit == 0 && nbrHexDigits != 1)
      {
        *decimal++ = ' ';
        currentGroupDigit = groupLength;
      }
    } while (--nbrHexDigits > 0);
  }
  strcpy(decimal, "</span>");
}

  // Convert little-endian number to a string with space every groupLen digits.
  // In order to perform a faster conversion, use groups of DIGITS_PER_LIMB digits.
void Bin2Dec(limb *binary, char *decimal, int nbrLimbs, int groupLength)
{
  int len, index, index2, count;
  limb *ptrSrc = binary + nbrLimbs - 1;
  char *ptrDest;
  int significantZero = 0;
  int groupCtr;
  int digit[DIGITS_PER_LIMB];
  int digits=0;
  int showDigitsText = TRUE;

  if (groupLength <= 0)
  {
    groupLength = -groupLength;
    showDigitsText = FALSE;
  }
  power10000[0].x = ptrSrc->x % MAX_LIMB_CONVERSION;
  power10000[1].x = ptrSrc->x / MAX_LIMB_CONVERSION;
  len = (power10000[1].x == 0 ? 1 : 2); // Initialize array length.
  for (index = nbrLimbs - 2; index >= 0; index--)
  {
    double dCarry, dQuotient;
    limb *ptrPower;

    // Multiply by FIRST_MULT and then by SECOND_MULT, so there is never
    // more than 53 bits in the product.

    ptrPower = power10000;
    dQuotient = 0;
    for (index2 = 0; index2 < len; index2++)
    {
      dCarry = dQuotient + (double)ptrPower->x * FIRST_MULT;
      dQuotient = floor(dCarry / MAX_LIMB_CONVERSION);
      (ptrPower++)->x = (int)(dCarry - dQuotient * MAX_LIMB_CONVERSION);
    }
    if (dQuotient != 0)
    {
      (ptrPower++)->x = (int)dQuotient;
      len++;
    }
    ptrPower = power10000;
    dQuotient = (--ptrSrc)->x;
    for (index2 = 0; index2 < len; index2++)
    {
      dCarry = dQuotient + (double)ptrPower->x * SECOND_MULT;
      dQuotient = floor(dCarry / MAX_LIMB_CONVERSION);
      (ptrPower++)->x = (int)(dCarry - dQuotient * MAX_LIMB_CONVERSION);
    }
    if (dQuotient != 0)
    {
      (ptrPower++)->x = (int)dQuotient;
      len++;
    }
  }
  // At this moment the array power10000 has the representation
  // of the number in base 10000 in little-endian. Convert to
  // ASCII separating every groupLength characters.
  ptrDest = decimal;
  ptrSrc = &power10000[len-1];
  groupCtr = len * DIGITS_PER_LIMB;
  if (groupLength != 0)
  {
    groupCtr %= groupLength;
    if (groupCtr == 0)
    {
      groupCtr = groupLength;
    }
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
        groupCtr = groupLength;
      }
    }
  }
  if (significantZero == 0)
  {     // Number is zero.
    *ptrDest++ = '0';
    *ptrDest = '\0';
    return;
  }
  if (digits > 30 && showDigitsText)
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
  unsigned int carry;
  int i;
  carry = 0;
  for (i = 0; i < length; i++)
  {
    carry += (unsigned int)(addend1++)->x + (unsigned int)(addend2++)->x;
    if (carry >= LIMB_RANGE)
    {
      (sum++)->x = (int)(carry - LIMB_RANGE);
      carry = 1;
    }
    else
    {
      (sum++)->x = (int)carry;
      carry = 0;
    }
  }
  return;
}

void BigInteger2Dec(BigInteger *pBigInt, char *decimal, int groupLength)
{
  if (pBigInt->sign == SIGN_NEGATIVE)
  {
    *decimal++ = '-';
  }
  Bin2Dec(pBigInt->limbs, decimal, pBigInt->nbrLimbs, groupLength);
}

void BigInteger2Hex(BigInteger *pBigInt, char *decimal, int groupLength)
{
  if (pBigInt->sign == SIGN_NEGATIVE)
  {
    *decimal++ = '-';
  }
  Bin2Hex(pBigInt->limbs, decimal, pBigInt->nbrLimbs, groupLength);
}

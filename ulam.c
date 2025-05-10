//
// This file is part of Alpertron Calculators.
//
// Copyright 2018-2021 Dario Alejandro Alpern
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
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "isprime.h"
#include "graphics.h"
#include "copyStr.h"
#if !defined(__EMSCRIPTEN__) && TEST_GRAPHICS
  #include <SDL.h>
#endif

#define NBR_SMALL_PRIMES   25
#if !defined(__EMSCRIPTEN__) && TEST_GRAPHICS
  SDL_Surface* doubleBuffer;
  int oldXCenter;
  int oldYCenter;
  int oldXFraction;
  int oldYFraction;
#else
  #ifdef __ANDROID__
    extern unsigned int *pixelArrPtr;
  #else          // Emscripten
    unsigned int pixelArray[PIXEL_ARRAY_SIZE];
  #endif
#endif

#ifdef __EMSCRIPTEN__
  extern setPointFunc setPoint;
  extern getInfoFunc getInfo;
  extern nbrChangedFunc nbrChgd;
#endif
#define MAX_LINES  1000
#define MAX_COLUMNS 2000
static char infoText[500];
static bool showAlgebraic;
static int startNumber[2] = {1, 0};
static int multiple[NBR_SMALL_PRIMES][97];
static bool initMultipleArrayCalled;
extern int thickness;             // Number of pixels for each square.
extern int xCenter;
extern int xFraction;             // Range of fraction: 0 to thickness - 1
extern int yCenter;
extern int yFraction;
extern int width;
extern int height;

static void initMultipleArray(void)
{
  if (initMultipleArrayCalled == false)
  {
    for (int i = 0; i<NBR_SMALL_PRIMES; i++)
    {
      int k = getPrime(i);
      for (int j = (k / 2) + 1; j >= 0; j--)
      {
        multiple[i][j*j%k] = 1;
      }
    }
    initMultipleArrayCalled = true;
  }
}

// Check whether 4x^2 + bx + c has integer roots (b is a very small number).
bool algebraicFactor(int linear, int *indep)
{
  int temp[2];
  int iSqDelta;
  int t1;
  int t2;
  double dLinear;
  double dDelta;
  double dFourAC;
  double dSqDelta;
  temp[0] = *indep;
  temp[1] = *(indep+1);
  if (((unsigned int)temp[1] & HALF_INT_RANGE_U) != 0U)
  {    // Independent term is negative.
    int carry = -temp[0];
    unsigned int tmp = (unsigned int)carry & MAX_VALUE_LIMB;
    temp[0] = (int)tmp;
    carry = ((carry >= 0)? -temp[1]: (MAX_INT_NBR - temp[1]));
    tmp = (unsigned int)carry & MAX_VALUE_LIMB;
    temp[1] = (int)tmp;
  }
  dFourAC = (((double)temp[1] * MAX_VALUE_LIMB) + (double)temp[0]) * 16;
  dLinear = (double)linear;
  if (((unsigned int)*(indep + 1) & HALF_INT_RANGE_U) != 0U)
  {    // Independent term is negative.
    dDelta = (dLinear * dLinear) + dFourAC;
  }
  else
  {    // Independent term is positive.
    dDelta = (dLinear * dLinear) - dFourAC;
  }
  if (dDelta < 0.0)
  {
    return false;   // No real roots.
  }
  dSqDelta = sqrt(dDelta) + 0.5;                     // Convert to nearest integer.
  iSqDelta = (int)dSqDelta;
  t1 = (int)(((double)(-linear) + iSqDelta) / 4);    // t1 and t2 are less than 2^31 so
  t2 = (int)(((double)(-linear) - iSqDelta) / 4);    // they fit into a double.

  multiplyBigNbrs((2*t1) + linear, t1, temp);
  AddBigNbrs(temp, indep, temp);
  AddBigNbrs(temp, indep, temp);
  if ((temp[0] != 0) || (temp[1] != 0))
  {
    return false;
  }
  multiplyBigNbrs((2*t2) + linear, t2, temp);
  AddBigNbrs(temp, indep, temp);
  AddBigNbrs(temp, indep, temp);
  if ((temp[0] != 0) || (temp[1] != 0))
  {
    return false;
  }
  return true;  // Twice the roots are integer numbers
}

static void getN(int valX, int valY, int *value)
{
  int x = valX;
  int y = valY;
  int addend[2];
  unsigned int tmp;
  if ((x >= 0) && (x >= y) && (x >= -y))
  {                     // Right quadrant.
    multiplyBigNbrs((4 * x) + 3, x, value);
    tmp = (unsigned int)y & MAX_VALUE_LIMB;
    addend[0] = (int)tmp;
    addend[1] = ((y >= 0)? 0: MAX_INT_NBR);
  }
  else if ((y >= 0) && (y > x) && (y > -x))
  {                     // Top quadrant.
    multiplyBigNbrs((4 * y) - 3, y, value);
    x = -x;
    tmp = (unsigned int)x & MAX_VALUE_LIMB;
    addend[0] = (int)tmp;
    addend[1] = ((x >= 0) ? 0 : MAX_INT_NBR);
  }
  else if ((x <= 0) && (x <= y) && (x <= -y))
  {                     // Left quadrant.
    multiplyBigNbrs((4 * x) + 1, x, value);
    y = -y;
    tmp = (unsigned int)y & MAX_VALUE_LIMB;
    addend[0] = (int)tmp;
    addend[1] = ((y >= 0) ? 0 : MAX_INT_NBR);
  }
  else
  {                     // Bottom quadrant.
    multiplyBigNbrs((4 * y) - 1, y, value);
    tmp = (unsigned int)x & MAX_VALUE_LIMB;
    addend[0] = (int)tmp;
    addend[1] = ((x >= 0) ? 0 : MAX_INT_NBR);
  }
  AddBigNbrs(value, addend, value);
  AddBigNbrs(value, startNumber, value);
}

void AddTwoLimbsPlusOneLimb(const int *addend1, int addend2, int *sum)
{
  int temp[2];
  temp[1] = 0;
  if (addend2 >= 0)
  {
    temp[0] = addend2;
    AddBigNbrs(addend1, temp, sum);
  }
  else
  {
    temp[0] = -addend2;
    SubtBigNbrs(addend1, temp, sum);
  }
}

void setPointUlamSpiral(unsigned int *ptrPixelArr, int x, int y)
{
  int value[2];
  int xPhysical;
  int yPhysical;
  int absx;
  int absy;
  int currY;
  int row;
  int col;
  int firstRow;
  int firstCol;
  int lastRow;
  int lastCol;
  int t;
  int linear1;
  int linear2;
  int indep1[2];
  int indep2[2];
  int Nminustt[2];
  unsigned int algebraicColor;
  unsigned int color;
#if !defined(__EMSCRIPTEN__) && TEST_GRAPHICS
  Uint32 colorBlack = SDL_MapRGBA(doubleBuffer->format, 0, 0, 0, 255);
  Uint32 colorBlue = SDL_MapRGBA(doubleBuffer->format, 0, 0, 255, 255);
  Uint32 colorGreen = SDL_MapRGBA(doubleBuffer->format, 0, 192, 0, 255);
  Uint32 colorWhite = SDL_MapRGBA(doubleBuffer->format, 192, 192, 192, 255);
#else
  unsigned int colorBlack = 0xFF000000U;
  unsigned int colorBlue = 0xFFFF0000U;
  unsigned int colorGreen = 0xFF00C000U;
  unsigned int colorWhite = 0xFFC0C0C0U;
#endif
  unsigned int *ptrPixel;
  initMultipleArray();
  getN(x, y, value);
  if (showAlgebraic)
  {                          // Color blue.
    algebraicColor = colorBlue;
  }
  else
  {                          // Color black.
    algebraicColor = colorBlack;
  }
  xPhysical = (width / 2) + ((x - xCenter) * thickness) - xFraction;
  yPhysical = (height / 2) - ((y - yCenter) * thickness) + yFraction;
  if (((value[1] != 0) || (value[0] != 2)) && ((value[0] & 1) == 0))
  {     // value is not 2 and it is even
    color = algebraicColor;
  }
  else if (isPrime(value))
  {                          // Color green.
    color = colorGreen;
  }
  else
  {
    if ((x > y) && (x > -y))
    {
      t = x;
    }
    else if ((x < y) && (x < -y))
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
    multiplyBigNbrs(t, t, Nminustt);
    SubtBigNbrs(value, Nminustt, Nminustt);
    if ((x + y) >= 0)
    {
      if ((x - y) >= 0)
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
      if ((x - y) >= 0)
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
  firstCol = ((xPhysical < 0)? 0: xPhysical);
  firstRow = ((yPhysical < 0)? 0: yPhysical);
  lastCol = xPhysical + thickness;
  if (lastCol > width)
  {
    lastCol = width;
  }
  lastRow = yPhysical + thickness;
  if (lastRow > height)
  {
    lastRow = height;
  }
  for (row = firstRow; row < lastRow; row++)
  {
    ptrPixel = ptrPixelArr + (row * ROW_WIDTH) + firstCol;
    for (col = firstCol; col < lastCol; col++)
    {
      *ptrPixel = color;
      ptrPixel++;
    }
  }
  if (thickness >= 4)
  {
    color = colorWhite;
    absx = ((x > 0)? x : -x);
    absy = ((y > 0)? y : -y);
    if ((absx >= absy) && !((x == -y) && (y < 0)))
    {
      if ((xPhysical >= 0) && (xPhysical < width))
      {
        ptrPixel = ptrPixelArr + (firstRow * ROW_WIDTH) + xPhysical;
        for (row = firstRow; row < lastRow; row++)
        {
          *ptrPixel = color;
          ptrPixel += ROW_WIDTH;
        }
      }
    }
    if (((absx < absy) && !(x == (y - 1) && (y > 0))) ||
        ((absx == absy) && (y <= 0)))
    {
      currY = yPhysical + thickness - 1;
      if ((currY >= 0) && (currY < height))
      {
        ptrPixel = ptrPixelArr + (currY * ROW_WIDTH) + firstCol;
        for (col = firstCol; col < lastCol; col++)
        {
          *ptrPixel = color;
          ptrPixel++;
        }
      }
    }
  }
}     /* end method setPoint */

// Convert the number from binary to string.
char *appendInt64(char *text, const int *value)
{
  bool valueIsNegative = false;
  int index;
  int index2;
  int nbrLo = *value;
  int nbrHi = *(value+1);
  if (nbrHi >= HALF_INT_RANGE)
  {       // Number is negative.
    if (nbrLo == 0)
    {
      nbrHi = -nbrHi & (int)MAX_VALUE_LIMB;
    }
    else
    {
      nbrLo = -nbrLo & (int)MAX_VALUE_LIMB;
      nbrHi = (-1 - nbrHi) & (int)MAX_VALUE_LIMB;
    }
    valueIsNegative = true;
  }
  // Convert digits from right to left.
  char* ptrText = text + 19;
  for (index = 0; index < 19; index++)
  {
    // Divide value by 10 and obtain the remainder.
    double dDivid;
    int rem = nbrHi % 10;
    nbrHi = nbrHi / 10;
    dDivid = (double)rem * (double)LIMB_RANGE + (double)nbrLo;
    nbrLo = (int)(dDivid / 10.0);
    ptrText--;
    *ptrText = (char)((int)(dDivid - (double)nbrLo * 10.0) + '0');
  }
  for (index = 0; index < 18; index++)
  {
    if (*(ptrText+index) != '0')
    {
      break;    // Most significant digit found.
    }
  }
  for (index2 = index; index2 < 19; index2++)
  {
    *(ptrText+index2-index) = *(ptrText+index2);
  }
  ptrText += index2-index;
  *ptrText = 0;    // Set string terminator.
  if (valueIsNegative)
  {
    char* ptrText2;
    ptrText++;     // Point to new string terminator.
    ptrText2 = ptrText;
    do
    {
      *ptrText2 = *(ptrText2-1);
      ptrText2--;
    } while (ptrText2 > text);
    *text = '-';
  }
  return ptrText;
}

static void ShowLabel(char **pptrText, const char *text, int linear, int *indep)
{
  int b = linear;
  int temp[2];
  int carry;
  char *ptrText = *pptrText;
  copyStr(&ptrText, text);
  temp[0] = *indep;
  temp[1] = *(indep+1);
  if ((temp[1] != 0) || (temp[0] != 0))
  {      // Independent term is not zero.
    *ptrText = ' ';
    ptrText++;
    if (((unsigned int)temp[1] & HALF_INT_RANGE_U) == 0U)
    {    // Independent term is positive.
      *ptrText = '+';
      ptrText++;
    }
    else
    {    // Independent term is negative.
      unsigned int tmp;
      *ptrText = '-';
      ptrText++;
      // Change its sign.
      carry = -temp[0];
      tmp = (unsigned int)carry & MAX_VALUE_LIMB;
      temp[0] = (int)tmp;
      carry = ((carry >= 0)? - temp[1]: (MAX_INT_NBR - temp[1]));
      tmp = (unsigned int)carry & MAX_VALUE_LIMB;
      temp[1] = (int)tmp;
    }
    *ptrText = ' ';
    ptrText++;
    ptrText = appendInt64(ptrText, temp);
  }
  if ((b != 0) || (temp[0] != 0) || (temp[1] != 0))
  {      // linear and independent terms are not zero.
    if (algebraicFactor(b, indep))
    {
      double dDelta;
      double dB = (double)b;
      double dFourAC = (((double)temp[1] * (double)MAX_VALUE_LIMB) + (double)temp[0]) * 16;
      if (((unsigned int)*(indep + 1) & HALF_INT_RANGE_U) != 0U)
      {    // Independent term is negative.
        dDelta = (dB * dB) + dFourAC;
      }
      else
      {    // Independent term is positive.
        dDelta = (dB * dB) - dFourAC;
      }
      double dSqDelta = sqrt(dDelta) + 0.5;      // Convert to nearest integer.
      int iSqDelta = (int)dSqDelta;
      int t1 = (int)((-dB + iSqDelta) / 4);         // t1 and t2 are less than 2^31 so
      int t2 = (int)((-dB - iSqDelta) / 4);         // they fit into a double.
                                                    // Twice the roots are integer numbers
      if (t1 > 0)
      {
        copyStr(&ptrText, " = (2t + ");
        ptrText = appendInt(ptrText, t1);
        *ptrText = ')';
        ptrText++;
      }
      else if (t1 < 0)
      {
        copyStr(&ptrText, " = (2t - ");
        ptrText = appendInt(ptrText, -t1);
        *ptrText = ')';
        ptrText++;
      }
      else
      {
        copyStr(&ptrText, " = 2t");
      }
      if (t2 > 0)
      {
        copyStr(&ptrText, " (2t + ");
        ptrText = appendInt(ptrText, t2);
        *ptrText = ')';
        ptrText++;
      }
      else if (t2 < 0)
      {
        copyStr(&ptrText, " (2t - ");
        ptrText = appendInt(ptrText, -t2);
        *ptrText = ')';
        ptrText++;
      }
      else
      {
        copyStr(&ptrText, " 2t");
      }
    }
    else
    {
      unsigned int tmp;
      char textSign = '+';
      temp[0] = *indep;
      temp[1] = *(indep+1);
      if (((unsigned int)temp[1] & HALF_INT_RANGE_U) != 0U)
      {     // Independent term is negative.
        textSign = '-';
        // Change sign of number.
        carry = -temp[0];
        tmp = (unsigned int)carry & MAX_VALUE_LIMB;
        temp[0] = (int)tmp;
        carry = ((carry >= 0) ? -temp[1] : (MAX_INT_NBR - temp[1]));
        tmp = (unsigned int)carry & MAX_VALUE_LIMB;
        temp[1] = (int)tmp;
      }
      if ((temp[0] & 1) == 0)
      {           // Independent term is even
        if (((b & 3) == 0) && ((temp[0] & 3) == 0))
        {         // Both linear and independent term are multiple of 4.
          copyStr(&ptrText, " = 4 (t<sup>2</sup>");
          b /= 4;
          // Divide by temp[1]: temp[0] by 4.
          tmp = (((unsigned int)temp[0] >> 2) |
            ((unsigned int)temp[1] << BITS_PER_GROUP_MINUS_2)) & MAX_VALUE_LIMB;
          temp[0] = (int)tmp;
          tmp = (unsigned int)temp[1] >> 2;
          temp[1] = (int)tmp;
        }
        else
        {
          copyStr(&ptrText, " = 2 (2t<sup>2</sup>");
          b /= 2;
          // Divide by temp[1]: temp[0] by 2.
          tmp = (((unsigned int)temp[0] >> 1) |
            ((unsigned int)temp[1] << BITS_PER_GROUP_MINUS_1)) & MAX_VALUE_LIMB;
          temp[0] = (int)tmp;
          tmp = (unsigned int)temp[1] >> 1;
          temp[1] = (int)tmp;
        }
        if (b != 0)
        {
          if (b < 0)
          {
            copyStr(&ptrText, " - ");
            b = -b;
          }
          else
          {
            copyStr(&ptrText, " + ");
          }
          if (b != 1)
          {
            ptrText = appendInt(ptrText, b);
          }
          copyStr(&ptrText, "t");
        }
        *ptrText = ' ';
        ptrText++;
        *ptrText = textSign;
        ptrText++;
        *ptrText = ' ';
        ptrText++;
        ptrText = appendInt64(ptrText, temp);
        copyStr(&ptrText, ")");
      }
      else
      {
        bool firstTime = true;
        for (int i = 0; i < NBR_SMALL_PRIMES; i++)
        {
          int deltaModP;
          int indepModP;
          int p = getPrime(i);
          // delta = b*b - 16*indep
          // Find delta mod p.
          indepModP = (((*(indep+1)%p)*((MAX_INT_NBR-p+1)%p)) + *(indep)%p) % p;
          deltaModP = (((((b%p) * (b%p)) - (16*indepModP)) % p) + p) % p;
          if (multiple[i][deltaModP] == 0)
          {
            if (firstTime)
            {
              firstTime = false;
              *ptrText = ' ';
              ptrText++;
              *ptrText = '(';
              ptrText++;
            }
            else
            {
              *ptrText = ',';
              ptrText++;
              *ptrText = ' ';
              ptrText++;
            }
            ptrText = appendInt(ptrText, p);
          }
        }         /* end for */
        if (firstTime == false)
        {
          copyStr(&ptrText, ")");
        }
      }           /* end if */
    }             /* end if */
  }               /* end if */
  copyStr(&ptrText, "<br>");
  *pptrText = ptrText;
}

char *ulamGetInformation(int x, int y)
{
  int value[2];
  char *ptrText = infoText;

  infoText[0] = 0;   // Empty string.
  if (x>=0)
  {
    int t;
    int Nminustt[2];
    int indep[2];
    int adjust = 16384 / thickness;
    int xLogical = xCenter + ((xFraction + x - (width / 2) + 16384) / thickness) - adjust;
    int yLogical = yCenter + 1 + ((yFraction - y + (height / 2) + 16384) / thickness) - adjust;
    getN(xLogical, yLogical, value);
    if ((xLogical > yLogical) && (xLogical > -yLogical))
    {
      t = xLogical;
    }
    else if ((xLogical < yLogical) && (xLogical < -yLogical))
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
    multiplyBigNbrs(t, t, Nminustt);
    SubtBigNbrs(value, Nminustt, Nminustt);
    if ((xLogical + yLogical) >= 0)
    {
      if ((xLogical - yLogical) >= 0)
      {                              // Right quadrant
        AddTwoLimbsPlusOneLimb(Nminustt, -2*t, indep);  // n - 4 * t*t - 4 * t
        ShowLabel(&ptrText, "SW-NE: 4t<sup>2</sup> + 4t", 4, indep);
        AddTwoLimbsPlusOneLimb(Nminustt, -t, indep);    // n - 4 * t*t - 2 * t
        ShowLabel(&ptrText, "NW-SE: 4t<sup>2</sup> + 2t", 2, indep);
      }
      else
      {                              // Upper quadrant
        AddTwoLimbsPlusOneLimb(Nminustt, 2*t, indep);  // n - 4 * t*t + 4 * t
        ShowLabel(&ptrText, "SW-NE: 4t<sup>2</sup> - 4t", -4, indep);
        AddTwoLimbsPlusOneLimb(Nminustt, t, indep);    // n - 4 * t*t + 2 * t
        ShowLabel(&ptrText, "NW-SE: 4t<sup>2</sup> - 2t", -2, indep);
      }
    }
    else
    {
      if ((xLogical - yLogical) >= 0)
      {                              // Lower quadrant
        ShowLabel(&ptrText, "SW-NE: 4t<sup>2</sup>", 0, Nminustt);
        AddTwoLimbsPlusOneLimb(Nminustt, -t, indep);   // n - 4 * t*t - 2 * t
        ShowLabel(&ptrText, "NW-SE: 4t<sup>2</sup> + 2t", 2, indep);
      }
      else
      {                             // Left quadrant
        ShowLabel(&ptrText, "SW-NE: 4t<sup>2</sup>", 0, Nminustt);
        AddTwoLimbsPlusOneLimb(Nminustt, t, indep);    // n - 4 * t*t + 2 * t
        ShowLabel(&ptrText, "NW-SE: 4t<sup>2</sup> - 2t", -2, indep);
      }
    }
    *ptrText = 'x';
    ptrText++;
    *ptrText = '=';
    ptrText++;
    ptrText = appendInt(ptrText, xLogical);
    *ptrText = ',';
    ptrText++;
    *ptrText = ' ';
    ptrText++;
    *ptrText = 'y';
    ptrText++;
    *ptrText = '=';
    ptrText++;
    ptrText = appendInt(ptrText, yLogical);
    *ptrText = ',';
    ptrText++;
    *ptrText = ' ';
    ptrText++;
    *ptrText = 'n';
    ptrText++;
    *ptrText = '=';
    ptrText++;
    *ptrText = '<';
    ptrText++;
    *ptrText = 'b';
    ptrText++;
    *ptrText = '>';
    ptrText++;
    ptrText = appendInt64(ptrText, value);
    *ptrText = '<';
    ptrText++;
    *ptrText = '/';
    ptrText++;
    *ptrText = 'b';
    ptrText++;
    *ptrText = '>';
    ptrText++;
    *ptrText = ',';
    ptrText++;
    *ptrText = ' ';
    ptrText++;
    *ptrText = 't';
    ptrText++;
    *ptrText = '=';
    ptrText++;
    ptrText = appendInt(ptrText, t/2);
  }
  getN(xCenter, yCenter, value);
  *ptrText = '^';       // Separator between bottom text and center text.
  ptrText++;
  ptrText = appendInt64(ptrText, value);
  *ptrText = 0;
  return infoText;
}  

// inputBoxNbr = 1 -> changing center
// inputBoxNbr = 2 -> changing start value
char *ulamNbrChanged(const char* value, int inputBoxNbr, int newWidth, int newHeight)
{
  int temp[2];
  int nbrLo;
  int nbrHi;
  unsigned int tmp;
  width = newWidth;
  height = newHeight;
  getValue64(value, &nbrLo, &nbrHi);
  if (inputBoxNbr == 1)
  {           // Changing center
              // nbr <- nbr - startNumber + 1 
    unsigned int carry;
    int borrow = nbrLo - startNumber[0];
    tmp = (unsigned int)borrow & MAX_VALUE_LIMB;
    nbrLo = (int)tmp;
    nbrHi += ((borrow >= 0) ? -startNumber[1] : (MAX_INT_NBR - startNumber[1]));
    tmp = (unsigned int)nbrHi & MAX_VALUE_LIMB;
    nbrHi = (int)tmp;
    carry = nbrLo + 1;
    tmp = carry & MAX_VALUE_LIMB;
    nbrLo = (int)tmp;
    nbrHi += (int)(carry >> BITS_PER_GROUP);
    if ((nbrHi > 0) || (nbrLo >= 1))
    {       // nbr >= 1
      int diff;
      int a = ((int)sqrt((((double)nbrHi * (double)LIMB_RANGE) + nbrLo - 1))+1)/2;
      diff = nbrLo - (4 * a * a);
      tmp = (unsigned int)diff & MAX_VALUE_LIMB;
      diff = (int)tmp;
      if (((unsigned int)diff & HALF_INT_RANGE_U) != 0U)
      {     // Number is negative.
        diff -= MAX_INT_NBR;
      }
      if (diff < (2 - (2*a)))
      {
        xCenter = (-3*a) + 1 - diff;
        yCenter = a;
      }
      else if (diff <= 0)
      {
        xCenter = -a;
        yCenter = -a + 1 - diff;
      }
      else if (diff <= (2*a))
      {
        xCenter = diff - a - 1;
        yCenter = -a;
      }
      else
      {
        xCenter = a;
        yCenter = diff - (3*a) - 1;
      }
    }
    xFraction = 0;
    yFraction = 0;
  }
  else if (inputBoxNbr == 2)
  {           // Changing start number.
    startNumber[0] = nbrLo;
    startNumber[1] = nbrHi;
  }
  else if (inputBoxNbr == 3)
  {           // Zoom in
    if (thickness == 32)
    {
      return NULL;
    }
    thickness *= 2;
    xFraction <<= 1;
    yFraction <<= 1;
  }
  else
  {           // Zoom out
    if (thickness == 1)
    {
      return NULL;
    }
    thickness /= 2;
    xFraction >>= 1;
    yFraction >>= 1;
  }
  getN(xCenter, yCenter, temp);
  infoText[0] = '^';
  (void)appendInt64(&infoText[1], temp);
  return infoText;
}

void initUlam(void)
{
#ifdef __EMSCRIPTEN__
  setPoint = setPointUlamSpiral;
  getInfo = ulamGetInformation;
  nbrChgd = ulamNbrChanged;
#endif
}

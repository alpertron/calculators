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

#define NBR_SMALL_PRIMES   25
#ifndef __EMSCRIPTEN__
  #define EXTERNALIZE	
  #include <SDL.h>
  SDL_Surface* screen;
  SDL_Surface* doubleBuffer;
  int oldXCenter;
  int oldYCenter;
  int oldXFraction;
  int oldYFraction;
  int timer;
  int quit;
#else     // Emscripten
  #define EXTERNALIZE  __attribute__((visibility("default")))
  #define MAX_WIDTH 2048
  unsigned int pixelArray[2048 * MAX_WIDTH];
#endif

#define MAX_LINES  1000
#define MAX_COLUMNS 2000
static char infoText[500];
static bool showAlgebraic;
static int thickness = 8;         // Number of pixels for each square.
static int xCenter;
static int xFraction;             // Range of fraction: 0 to thickness - 1
static int yCenter;
static int yFraction;
static int width;
static int height;
static int startNumber[2] = {1, 0};
static int multiple[NBR_SMALL_PRIMES][97];
static bool initMultipleArrayCalled;

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

  multiply((2*t1) + linear, t1, temp);
  AddBigNbr(temp, indep, temp);
  AddBigNbr(temp, indep, temp);
  if ((temp[0] != 0) || (temp[1] != 0))
  {
    return false;
  }
  multiply((2*t2) + linear, t2, temp);
  AddBigNbr(temp, indep, temp);
  AddBigNbr(temp, indep, temp);
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
    multiply((4 * x) + 3, x, value);
    tmp = (unsigned int)y & MAX_VALUE_LIMB;
    addend[0] = (int)tmp;
    addend[1] = ((y >= 0)? 0: MAX_INT_NBR);
  }
  else if ((y >= 0) && (y > x) && (y > -x))
  {                     // Top quadrant.
    multiply((4 * y) - 3, y, value);
    x = -x;
    tmp = (unsigned int)x & MAX_VALUE_LIMB;
    addend[0] = (int)tmp;
    addend[1] = ((x >= 0) ? 0 : MAX_INT_NBR);
  }
  else if ((x <= 0) && (x <= y) && (x <= -y))
  {                     // Left quadrant.
    multiply((4 * x) + 1, x, value);
    y = -y;
    tmp = (unsigned int)y & MAX_VALUE_LIMB;
    addend[0] = (int)tmp;
    addend[1] = ((y >= 0) ? 0 : MAX_INT_NBR);
  }
  else
  {                     // Bottom quadrant.
    multiply((4 * y) - 1, y, value);
    tmp = (unsigned int)x & MAX_VALUE_LIMB;
    addend[0] = (int)tmp;
    addend[1] = ((x >= 0) ? 0 : MAX_INT_NBR);
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
    multiply(t, t, Nminustt);
    SubtBigNbr(value, Nminustt, Nminustt);
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
#ifdef EMSCRIPTEN
    ptrPixel = &pixelArray[(row * MAX_WIDTH) + firstCol];
#else
    ptrPixel = (Uint32*)doubleBuffer->pixels + (row * width) + firstCol;
#endif
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
#ifdef EMSCRIPTEN
        ptrPixel = &pixelArray[(firstRow * MAX_WIDTH) + xPhysical];
#else
        ptrPixel = (Uint32*)doubleBuffer->pixels + (firstRow * width) + xPhysical;
#endif
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
    if (((absx < absy) && !(x == (y - 1) && (y > 0))) ||
        ((absx == absy) && (y <= 0)))
    {
      currY = yPhysical + thickness - 1;
      if ((currY >= 0) && (currY < height))
      {
#ifdef EMSCRIPTEN
        ptrPixel = &pixelArray[(currY * MAX_WIDTH) + firstCol];
#else
        ptrPixel = (Uint32*)doubleBuffer->pixels + (currY * width) + firstCol;
#endif
        for (col = firstCol; col < lastCol; col++)
        {
          *ptrPixel = color;
          ptrPixel++;
        }
      }
    }
  }
}     /* end method setPoint */

EXTERNALIZE void drawPartialUlamSpiral(int xminDisp, int xmaxDisp, int yminDisp, int ymaxDisp)
{
  int x;
  int y;
  int adjust = 16384 / thickness;
  int xmin = xCenter + ((xFraction + xminDisp + 16384) / thickness) - adjust;
  int xmax = xCenter + ((xFraction + xmaxDisp + 16384) / thickness) - adjust;
  int ymin = yCenter - ((-yFraction - yminDisp + 16384) / thickness) + adjust;
  int ymax = yCenter - ((-yFraction - ymaxDisp + 16384) / thickness) + adjust;

  initMultipleArray();
  for (x = xmin; x <= xmax; x++)
  {
    for (y = ymin; y <= ymax; y++)
    {
      setPoint(x, y);
    }  /* end for y */
  }    /* end for x */
}      /* end method drawUlamSpiral */

#ifdef __EMSCRIPTEN__
void copyStr(char** pptrString, const char* stringToCopy)
{
  char* ptrString = *pptrString;
  const char* ptrStringToCopy = stringToCopy;
  while (*ptrStringToCopy != '\0')
  {
    *ptrString = *ptrStringToCopy;
    ptrString++;
    ptrStringToCopy++;
  }
  *ptrString = '\0';
  *pptrString = ptrString;
}

unsigned int *getPixels(void)
{
  return pixelArray;
}
#endif

// Convert the number from binary to string.
char *appendInt64(char *text, const int *value)
{
  int index;
  int index2;
  int nbrLo = *value;
  int nbrHi = *(value+1);
  
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
    *ptrText = (int)(dDivid - (double)nbrLo * 10.0) + '0';
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
  return ptrText;
}

void ShowLabel(char **pptrText, char *text, int linear, int *indep)
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
      double dSqDelta = sqrt(dDelta) + 0.5;                 // Convert to nearest integer.
      int iSqDelta = (int)dSqDelta;
      int t1 = (int)((-dB + iSqDelta) / 4);                 // t1 and t2 are less than 2^31 so
      int t2 = (int)((-dB - iSqDelta) / 4);                 // they fit into a double.
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
        int i;
        bool firstTime = true;
        for (i = 0; i < NBR_SMALL_PRIMES; i++)
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

EXTERNALIZE char *getInformation(int x, int y)
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
    multiply(t, t, Nminustt);
    SubtBigNbr(value, Nminustt, Nminustt);
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

EXTERNALIZE void moveSpiral(int deltaX, int deltaY)
{
  int adjust = 16384 / thickness;
  xFraction -= deltaX;
  xCenter += ((xFraction + 16384) / thickness) - adjust;
  xFraction &= thickness - 1;
  yFraction += deltaY;
  yCenter += ((yFraction + 16384) / thickness) - adjust;
  yFraction &= thickness - 1; 
}

// inputBoxNbr = 1 -> changing center
// inputBoxNbr = 2 -> changing start value
EXTERNALIZE char *nbrChanged(char *value, int inputBoxNbr, int newWidth, int newHeight)
{
  char* ptrValue = value;
  int temp[2];
  int nbrLo = 0;
  int nbrHi = 0;
  int index;
  unsigned int tmp;
  width = newWidth;
  height = newHeight;
  for (index=0; index<19; index++)
  {
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

#ifndef __EMSCRIPTEN__

void drawUlamSpiral(void)
{
  if (SDL_MUSTLOCK(doubleBuffer) != 0)
  {
    SDL_LockSurface(doubleBuffer);
  }
  drawPartialUlamSpiral(-width / 2,
                        width / 2,
                        -height / 2,
                        height / 2);
  if (SDL_MUSTLOCK(doubleBuffer) != 0)
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
  SDL_Rect rectSrc;
  SDL_Rect rectDest;
        
  while (SDL_PollEvent(&event) != 0)
  {                           // New event arrived.
    if (event.type == SDL_QUIT)
    {                         // Quit event arrived, so exit application.
      quit = 1;
    }
    if (event.type == SDL_MOUSEMOTION)
    {
      if ((event.motion.state & SDL_BUTTON_LMASK) != 0)
      {                        // Drag operation.
        int adjust = 16384 / thickness;
        xFraction -= event.motion.xrel;
        xCenter += ((xFraction + 16384) / thickness) - adjust;
        xFraction &= thickness - 1;
        yFraction += event.motion.yrel;
        yCenter += ((yFraction + 16384) / thickness) - adjust;
        yFraction &= thickness - 1;
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
        if (thickness < 32)
        {
          thickness *= 2;
          xFraction <<= 1;
          yFraction <<= 1;
          drawUlamSpiral();
        }
      }
      else if (event.button.button == SDL_BUTTON_WHEELDOWN)
      {
        if (thickness > 1)
        {
          thickness /= 2;
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
  timer++;
  if (timer == 6)
  {
    timer = 0;
    if ((oldXCenter != xCenter) || (oldYCenter != yCenter) ||
        (oldXFraction != xFraction) || (oldYFraction != yFraction))
    {
      int xMin;
      int xMax;
      int yMin;
      int yMax;
      // Move pixels of double buffer according to drag direction.
      int xMove = (xCenter * thickness) + xFraction - (oldXCenter * thickness) - oldXFraction;
      int yMove = (yCenter * thickness) + yFraction - (oldYCenter * thickness) - oldYFraction;
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
      if (SDL_MUSTLOCK(doubleBuffer) != 0)
      {
        SDL_LockSurface(doubleBuffer);
      }
      xMin = -width / 2;
      xMax = width / 2;
      yMin = -height / 2;
      yMax = height / 2;
      if (yMove > 0)
      {              // Move pixels up.
        int yBound = yMax - yMove;
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
      if (SDL_MUSTLOCK(doubleBuffer) != 0)
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
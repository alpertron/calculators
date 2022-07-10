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
#include "isprime.h"
#include "graphics.h"
#ifndef __EMSCRIPTEN__
  #define EXTERNALIZE	
  #include <SDL.h>
#endif
static void setPointGaussian(int x, int y);
setPointFunc setPoint = setPointGaussian;

#define SMALL_NUMBER_BOUND 32768
#ifndef __EMSCRIPTEN__
  SDL_Surface* doubleBuffer;
  int oldXCenter;
  int oldYCenter;
  int oldXFraction;
  int oldYFraction;
#else        // Emscripten
  unsigned int pixelArray[PIXEL_ARRAY_SIZE];
#endif

#define MAX_LINES  1000
#define MAX_COLUMNS 2000
static char infoText[500];
int thickness = 8;              // Number of bits to shift.
int xCenter;
int xFraction;                  // Range of fraction: 0 to thickness - 1
int yCenter;
int yFraction;
int width;
int height;

 // A Gaussian number x+iy is prime if one of these 3 conditions is true:
 // 1) x is not zero and y is not zero and x^2+y^2 is prime
 // 2) x=0, y is prime and y=3 (mod 4)
 // 3) y=0, x is prime and x=3 (mod 4)
static void setPointGaussian(int x, int y)
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
  xPhysical = (width / 2) + ((x - xCenter) * thickness) - xFraction;
  yPhysical = (height / 2) - ((y - yCenter) * thickness) + yFraction;
  color = colorBlack;          // Indicate not prime in advance.
  if (x == 0)
  {   // Number is imaginary.
    check[0] = ((y >= 0)? y : -y);
    check[1] = 0;
    if (((check[0] & 0x03) == 3) && isPrime(check))
    {                          // Number is gaussian prime.
      color = colorGreen;
    }
  }
  else if (y == 0)
  {   // Number is real.
    check[0] = ((x >= 0)? x : -x);
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
    if (x == 0)
    {                     // Draw Y axis if possible.
      col = xPhysical + (thickness / 2);
      if ((col >= 0) && (col < width))
      {
#ifdef EMSCRIPTEN
        ptrPixel = &pixelArray[(firstRow * MAX_WIDTH) + col];
#else
        ptrPixel = (Uint32*)doubleBuffer->pixels + (firstRow * width) + col;
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
    else if ((x % 10) == 0)
    {
      col = xPhysical + (thickness / 2);
      if ((col >= 0) && (col < width))
      {
        firstRow2 = yPhysical + (thickness / 4);
        lastRow2 = firstRow2 + (thickness / 2);
        firstRow2 = ((firstRow2 < 0) ? 0 : firstRow2);
        if (lastRow2 > height)
        {
          lastRow2 = height;
        }
#ifdef EMSCRIPTEN
        ptrPixel = &pixelArray[(firstRow2 * MAX_WIDTH) + col];
#else
        ptrPixel = (Uint32*)doubleBuffer->pixels + (firstRow2 * width) + col;
#endif
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
    else
    {                     // Nothing to do.
    }
    if (y == 0)
    {                     // Draw X axis if possible.
      row = yPhysical + (thickness / 2);
      if ((row >= 0) && (row < height))
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
    }
    else if ((y % 10) == 0)
    {
      row = yPhysical + (thickness / 2);
      if ((row >= 0) && (row < height))
      {
        firstCol2 = xPhysical + (thickness / 4);
        lastCol2 = firstCol2 + (thickness / 2);
        firstCol2 = ((firstCol2 < 0)? 0 : firstCol2);
        if (lastCol2 > width)
        {
          lastCol2 = width;
        }
#ifdef EMSCRIPTEN
        ptrPixel = &pixelArray[(row * MAX_WIDTH) + firstCol2];
#else
        ptrPixel = (Uint32*)doubleBuffer->pixels + (row * width) + firstCol2;
#endif
        for (col = firstCol2; col < lastCol2; col++)
        {
          *ptrPixel = color;
          ptrPixel++;
        }
      }
    }
    else
    {              // Nothing to do.
    }
  }
}     /* end method setPoint */

char *getInformation(int x, int y)
{
  char *ptrText = infoText;
  
  infoText[0] = 0;   // Empty string.
  if (x >= 0)
  {
    int adjust = 16384 / thickness;
    int xLogical = xCenter + ((xFraction + x - (width / 2) + 16384) / thickness) - adjust;
    int yLogical = yCenter + 1 + ((yFraction - y + (height / 2) + 16384) / thickness) - adjust;
    ptrText = appendInt(infoText, xLogical);
    if (yLogical >= 0)
    {
      *ptrText = ' ';
      ptrText++;
      *ptrText = '+';
      ptrText++;
      *ptrText = ' ';
      ptrText++;
      ptrText = appendInt(ptrText, yLogical);
    }
    else
    {
      *ptrText = ' ';
      ptrText++;
      *ptrText = '-';
      ptrText++;
      *ptrText = ' ';
      ptrText++;
      ptrText = appendInt(ptrText, -yLogical);
    }
    *ptrText = 'i';
    ptrText++;
  }
  *ptrText = '^';       // Append separator.
  ptrText++;
  ptrText = appendInt(ptrText, xCenter);
  *ptrText = '^';       // Append separator.
  ptrText++;
  ptrText = appendInt(ptrText, yCenter);
  *ptrText = 0;
  return infoText;
}

static int getValue(char **ppValue)
{
  char *ptrValue = *ppValue;
  int nbr = 0;
  int sign = 0;
  if (*ptrValue == '-')
  {
    sign = 1;
    ptrValue++;
  }
  while (*ptrValue != '\0')
  {
    nbr = (nbr * 10) + (*ptrValue - '0');
    ptrValue++;
  }
  *ppValue = ptrValue;
  return sign? -nbr: nbr;
}

int nbrChanged(char *value, int inputBoxNbr, int newWidth, int newHeight)
{
  char* ptrValue = value;
  width = newWidth;
  height = newHeight;
  if (inputBoxNbr == 1)
  {           
    xCenter = getValue(&ptrValue);  // Changing center X.
    xFraction = 0;
    ptrValue++;
    yCenter = getValue(&ptrValue);  // Changing center Y.
    yFraction = 0;
  }
  else if (inputBoxNbr == 3)
  {           // Zoom in
    if (thickness == 32)
    {
      return 0;
    }
    thickness *= 2;
    xFraction <<= 1;
    yFraction <<= 1;
  }
  else
  {           // Zoom out
    if (thickness == 1)
    {
      return 0;
    }
    thickness /= 2;
    xFraction >>= 1;
    yFraction >>= 1;
  }
  return 0;
}

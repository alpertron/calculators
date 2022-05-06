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
#include <stdint.h>
#include "graphics.h"

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
bool quit;
#else     // Emscripten
#define EXTERNALIZE  __attribute__((visibility("default")))
extern unsigned int pixelArray[PIXEL_ARRAY_SIZE];
#endif

extern setPointFunc setPoint;
extern int thickness;                  // Number of bits to shift.
extern int xCenter;
extern int xFraction;                  // Range of fraction: 0 to thickness - 1
extern int yCenter;
extern int yFraction;
extern int width;
extern int height;

EXTERNALIZE void drawPartialGraphic(int xminDisp, int xmaxDisp, int yminDisp, int ymaxDisp)
{
  int x;
  int y;
  int adjust = 16384 / thickness;
  int xmin = xCenter + ((xFraction + xminDisp + 16384) / thickness) - adjust;
  int xmax = xCenter + ((xFraction + xmaxDisp + 16384) / thickness) - adjust;
  int ymin = yCenter - ((-yFraction - yminDisp + 16384) / thickness) + adjust;
  int ymax = yCenter - ((-yFraction - ymaxDisp + 16384) / thickness) + adjust;

  for (x = xmin; x <= xmax; x++)
  {
    for (y = ymin; y <= ymax; y++)
    {
      setPoint(x, y);
    }  /* end for y */
  }    /* end for x */
}      /* end method drawPartialGraphic */

EXTERNALIZE void moveGraphic(int deltaX, int deltaY)
{
  int adjust = 16384 / thickness;
  xFraction -= deltaX;
  xCenter += ((xFraction + 16384) / thickness) - adjust;
  xFraction &= thickness - 1;
  yFraction += deltaY;
  yCenter += ((yFraction + 16384) / thickness) - adjust;
  yFraction &= thickness - 1;
}

#ifdef __EMSCRIPTEN__
unsigned int* getPixels(void)
{
  return pixelArray;
}
#else
void drawGraphic(void)
{
  if (SDL_MUSTLOCK(doubleBuffer) != 0)
  {
    SDL_LockSurface(doubleBuffer);
  }
  drawPartialGraphic(-width / 2,
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
          drawGraphic();
        }
      }
      else if (event.button.button == SDL_BUTTON_WHEELDOWN)
      {
        if (thickness > 1)
        {
          thickness /= 2;
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

int main(int argc, char* argv[])
{
  width = 256;
  height = 256;
  initMultipleArray();
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

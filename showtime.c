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

#include <string.h>
#include <stdlib.h>
#include "bignbr.h"
#include "expression.h"
#include "showtime.h"
#include "output.h"

#ifdef __EMSCRIPTEN__
double originalTenthSecond;
int oldTimeElapsed;
// Convert tenths of seconds to days, hours, minutes and seconds.
void GetDHMS(char **pptrText, int seconds)
{
  char *ptrText = *pptrText;
  int2dec(&ptrText, seconds / 86400);         // Show number of days.
  *ptrText = 'd';
  ptrText++;
  *ptrText = ' ';
  ptrText++;
  int2dec(&ptrText, (seconds / 3600) % 24);   // Show number of hours.
  *ptrText = 'h';
  ptrText++;
  *ptrText = ' ';
  ptrText++;
  int2dec(&ptrText, (seconds / 60) % 60);     // Show number of minutes.
  *ptrText = 'm';
  ptrText++;
  *ptrText = ' ';
  ptrText++;
  int2dec(&ptrText, seconds % 60);            // Show number of seconds.
  *ptrText = 's';
  ptrText++;
  *ptrText = ' ';
  ptrText++;
  *pptrText = ptrText;
}

void GetDHMSt(char **pptrText, int tenths)
{
  char *ptrText;
  GetDHMS(pptrText, tenths / 10);
  ptrText = *pptrText - 2;
  *ptrText = '.';
  ptrText++;
  *ptrText = (char)((tenths % 10) + '0');
  ptrText++;
  *ptrText = 's';
  ptrText++;
  *ptrText = ' ';
  ptrText++;
  *pptrText = ptrText;
}

#endif

void showElapsedTime(char **pptrOutput)
{
  char *ptrOutput = *pptrOutput;
  copyStr(&ptrOutput, lang ? "Tiempo transcurrido: " : "Time elapsed: ");
#ifdef __EMSCRIPTEN__
  GetDHMSt(&ptrOutput, (int)(tenths() - originalTenthSecond));
#endif
  *pptrOutput = ptrOutput;
}

#ifdef __EMSCRIPTEN__
void showElapsedTimeSec(char **pptrOutput)
{
  char *ptrOutput = *pptrOutput;
  copyStr(&ptrOutput, lang ? "<p>Tiempo transcurrido: " : "<p>Time elapsed: ");
  GetDHMS(&ptrOutput, (int)(tenths() - originalTenthSecond) / 10);
  copyStr(&ptrOutput, "</p>");
  *pptrOutput = ptrOutput;
}
#endif

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
#include <stdio.h>
#include <stdlib.h>
#include "bignbr.h"
#include "highlevel.h"
#include "polynomial.h"
#include "showtime.h"

#ifdef __EMSCRIPTEN__
EXTERNALIZE void doWork(void)
{
  int flags;
  int groupLen = 0;
  char* ptrData = inputString;
  while (*ptrData != ',')
  {
    groupLen = groupLen * 10 + (*ptrData++ - '0');
  }
  ptrData++;             // Skip comma.
  flags = *ptrData;
#ifndef lang  
  lang = flags & 1;
#endif
  onlyEvaluate = (unsigned char)(flags & 2);
  pretty = (unsigned char)(flags & 4);
  ptrData += 2;          // Skip flags and comma.
  polyFactText(ptrData, ptrData + strlen(ptrData) + 1, groupLen);
  ptrData += strlen(ptrData);
  databack(output);
}
#endif

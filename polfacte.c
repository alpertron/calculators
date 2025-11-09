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
#include "highlevel.h"
#include "polynomial.h"
#include "showtime.h"

extern bool teach;
#if defined(__EMSCRIPTEN__) && !defined(_MSC_VER)
#ifdef __ANDROID__
EXTERNALIZE void doWorkPolfact(void)
#else
EXTERNALIZE void doWork(void)
#endif
{
  int flags = 0;
  int groupLength = 0;
  char* ptrData = inputString;
  while (*ptrData != ',')
  {
    groupLength = (groupLength * 10) + (*ptrData - '0');
    ptrData++;
  }
  ptrData++;         // Skip comma.
  do
  {
    flags = (flags * 10) + *ptrData - '0';
  } while (*(++ptrData) != ',');
  onlyEvaluate = ((flags & 2) != 0);
  teach = ((flags & 4) != 0);
  switch (flags / 8)
  {
  case 0:
    pretty = PRETTY_PRINT;
    break;
  case 1:
    pretty = TEX;
    break;
  default:
    pretty = PARI_GP;
    break;
  }
  ptrData++;          // Skip comma.
  polyFactText(ptrData, ptrData + (int)strlen(ptrData) + 1, groupLength);
  databack(output);
}
#endif

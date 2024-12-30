//
// This file is part of Alpertron Calculators.
//
// Copyright 2019-2021 Dario Alejandro Alpern
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
#include <stdbool.h>
#include "output.h"
#include "bignbr.h"

char* ptrOutput;
int eqNbr;
#ifndef lang  
extern bool lang;
#endif

void showText(const char* text)
{
  copyStr(&ptrOutput, text);
}

void shownbr(const BigInteger* value)
{
  BigInteger2Dec(&ptrOutput, value, groupLen);
}

void generateEqNbr(void)
{
  showText("<e-q>");
  eqNbr++;
  int2dec(&ptrOutput, eqNbr);
  showText("</e-q>");
}

void showEqNbrs(int eqNbr1, int eqNbr2)
{
  showText("(");
  int2dec(&ptrOutput, eqNbr1);
  showText(lang ? ") y (" : ") and (");
  int2dec(&ptrOutput, eqNbr2);
  showText(")");
}

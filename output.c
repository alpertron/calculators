/*
This file is part of Alpertron Calculators.

Copyright 2019 Dario Alejandro Alpern

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
#include "output.h"

char* ptrOutput;
void showText(char* text)
{
  strcpy(ptrOutput, text);
  ptrOutput += strlen(ptrOutput);
}

void showMinus(void)
{
  showText("&minus;");
}

void shownbr(BigInteger* value)
{
  BigInteger2Dec(value, ptrOutput, groupLen);
  ptrOutput += strlen(ptrOutput);
}

void showInt(int value)
{
  int2dec(&ptrOutput, value);
}

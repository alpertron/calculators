//
// This file is part of Alpertron Calculators.
//
// Copyright 2017-2021 Dario Alejandro Alpern
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
#ifndef _BATCH_H
#define _BATCH_H
extern int valuesProcessed;

enum eExprErr BatchProcessing(char *batchText, BigInteger *valueFound, char **pptrOutput, bool *pIsBatch);
void batchCallback(char **pptrOutput);
char *findChar(char *str, char c);
#ifdef __EMSCRIPTEN__
extern char *ptrInputText;
#endif
#endif

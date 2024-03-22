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
#ifndef _HIGHLEVEL_H
#define _HIGHLEVEL_H
void fsquaresText(char* input, int grpLen);
void fcubesText(char *input, int grpLen);
void tsqcubesText(char* input, int grpLen, int expon);
void contfracText(char *input, int GroupLen, bool hex, bool converg);
void polyFactText(const char *modText, const char *polyText, int groupLen);
void quadmodText(const char *quadrText, const char *linearText, const char *constText,
  const char *modText, int groupLength);
void quadText(char *coefAText, char *coefBText, char *coefCText,
              char *coefDText, char *coefEText, char *coefFText);
#endif

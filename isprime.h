//
// This file is part of Alpertron Calculators.
//
// Copyright 2021 Dario Alejandro Alpern
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
#ifndef _ISPRIME_H
#define _ISPRIME_H

#define NBR_LIMBS        2
#define BITS_PER_GROUP   31
#define LIMB_RANGE       0x80000000U
#define MAX_INT_NBR      0x7FFFFFFF
#define HALF_INT_RANGE   0x40000000
#define FOURTH_INT_RANGE 0x20000000
#define MAX_VALUE_LIMB   0x7FFFFFFFU

bool isPrime(int *value);
void multiply(int factor1, int factor2, int *prod);
void AddBigNbr(const int *Nbr1, const int *Nbr2, int *Sum);
void SubtBigNbr(const int *Nbr1, const int *Nbr2, int *Diff);
int getPrime(int index);
char* appendInt(char* text, int intValue);
#endif
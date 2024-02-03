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
#ifndef _TSQUARES_H
#define _TSQUARES_H
#define MAX_SIEVE 65536

bool FindTwoSquaresNoNumTheory(void);
void FillSieveArray(const limb *value, int *sieveArr);
bool FindSumTwoSquaresNoSmallDivisors(void);
void FindSumOfTwoSquaresUsingCrt(int nbrDivisors);
bool isSumOfTwoSquares(void);
void ShowStatus(void);
void InitSieveArray(void);
extern const char* square;
extern int sieve[MAX_SIEVE];
extern limb Mult1[MAX_LEN];
extern limb Mult2[MAX_LEN];
extern limb Mult3[MAX_LEN];
extern limb Mult4[MAX_LEN];
extern limb valueP[MAX_LEN];
extern limb valueQ[MAX_LEN];
extern limb number[MAX_LEN];
extern int Mult1Len;
extern int Mult2Len;
extern int Mult3Len;
extern int Mult4Len;
extern int nbrLimbsP;
extern int nbrModExp;
extern int nbrLimbs;
extern bool hexadecimal;
extern int attempts;
extern BigInteger SquareMult1;
extern BigInteger SquareMult2;
extern BigInteger SquareMult3;
extern BigInteger SquareMult4;
extern BigInteger biMult1;
extern BigInteger biMult2;
extern BigInteger biMult3;
extern BigInteger biMult4;
#endif

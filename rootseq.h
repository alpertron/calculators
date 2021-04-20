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
#include "output.h"
#include "polynomial.h"

#define TYPE_PM_SPACE_BEFORE 1
#define TYPE_PM_SPACE_AFTER  2

void QuinticEquation(int* ptrPolynomial, int multiplicity);
void startLine(void);
void endLine(void);
void showX(int multiplicity);
void showPlusSignOn(bool condPlus, int type);
void startParen(void);
void endParen(void);
void startSqrt(void);
void endSqrt(void);
void showRatConstants(char* numerator, char* denominator);
void getRootsPolynomial(int nbrFactor, char **pptrOutput, struct sFactorInfo* pstFactorInfo, int groupLength);

extern BigInteger Quintic;
extern BigInteger Quartic;
extern BigInteger Cubic;
extern BigInteger Quadratic;
extern BigInteger Linear;
extern BigInteger Independent;
extern BigInteger discr;
extern BigInteger commonDenom;
extern BigInteger tmp0;
extern BigInteger tmp1;
extern BigInteger tmp2;
extern BigInteger tmp3;
extern BigInteger tmp4;
extern BigInteger tmp5;
extern BigInteger tmp6;
extern BigInteger tmp7;
extern BigRational RatQuartic;
extern BigRational RatCubic;
extern BigRational RatQuadratic;
extern BigRational RatLinear;
extern BigRational RatIndependent;
extern BigRational RatDeprCubic;
extern BigRational RatDeprQuadratic;
extern BigRational RatDeprLinear;
extern BigRational RatDeprIndependent;
extern BigRational RatDiscr;
extern BigRational RatDelta0;
extern BigRational RatDelta1;
extern BigRational RatD;
extern BigRational Rat1;
extern BigRational Rat2;
extern BigRational Rat3;
extern BigRational Rat4;
extern BigRational Rat5;
extern BigRational RatS;
extern int indexRoot;
extern char *ptrMinus;
extern char *ptrTimes;
extern int polyInteger[1000000];

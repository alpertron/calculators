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
#ifndef _ROOTSEQ_H
#define _ROOTSEQ_H
#include "output.h"
#include "polynomial.h"

#define TYPE_PM_SPACE_BEFORE 1
#define TYPE_PM_SPACE_AFTER  2

void LinearEquation(const int* polynomial, int multiplicity);
void QuadraticEquation(const int* polynomial, int multiplicity);
void CubicEquation(const int* polynomial, int multiplicity);
void QuarticEquation(const int* polynomial, int multiplicity);
void QuinticEquation(const int* ptrPolynomial, int multiplicity);
bool isLinearExponential(const int* ptrPolynomial, int polyDegree, int multiplicity);
bool isQuadraticExponential(const int* ptrPolynomial, int polyDegree,
  int multiplicity);
void GenerateRoots(int multiplicity, const char* rationalRoot,
  bool isNegative, int polyDegree);
void StartRadicand(int polyDegree);
void EndRadicand(int polyDegree);
void showTimesPi(char** pptrString);
void showRatString(const char* num, const char* den);
void startLine(void);
void endLine(void);
void showX(int multiplicity);
void endShowX(void);
void showPlusSignOn(bool condPlus, int type);
void startParen(void);
void endParen(void);
void startSqrt(void);
void endSqrt(void);
void startCbrt(void);
void endCbrt(void);
void showRatConstants(const char* numerator, const char* denominator);
void getRootsPolynomial(int nbrFactor, char **pptrOutput, struct sFactorInfo* pstFactorInfo, int groupLength);
void showVariable(char **pptrOutput, char letter);
void stepsForQuadraticEquation(char origVar, char substVar);
void showVarIndex(char letter, int index);
void showCoeffBeforeParen(BigRational* rat);

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
extern const char *ptrMinus;
extern const char *ptrTimes;
extern int polyInteger[1000000];
extern bool teach;
#endif

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
#ifndef __POLYNOMIAL_H
#define __POLYNOMIAL_H
#define MAX_DEGREE          1000
#define COMPRESSED_POLY_MAX_LENGTH 1000000

#define DEBUG_VANHOEIJ        0

enum eOutput
{
  PRETTY_PRINT = 0,
  TEX = 1,
  PARI_GP = 2,
};

struct sFactorInfo
{
  int *ptr;
  int *ptrPolyLifted;
  int degree;
  int multiplicity;
  int expectedDegree;
};

enum eDivType
{
  TYPE_DIVISION,
  TYPE_MODULUS
};

#ifndef lang  
  extern bool lang;
#endif
extern BigInteger Quintic;
extern BigInteger Quartic;
extern BigInteger Cubic;
extern BigInteger Quadratic;
extern BigInteger Linear;
extern BigInteger Independent;
extern struct sFactorInfo factorInfo[MAX_DEGREE];
extern struct sFactorInfo factorInfoInteger[MAX_DEGREE];
extern BigInteger primeMod;              // p
extern int exponentMod;                  // k
extern BigInteger powerMod;              // p^k
extern bool modulusIsZero;
extern BigInteger operand1;
extern BigInteger operand2;
extern BigInteger operand3;
extern BigInteger operand4;
extern BigInteger operand5;
extern int values[1000000];
extern int poly1[1000000];
extern int poly2[1000000];
extern int poly3[1000000];
extern int polyMultTemp[1000000];
extern int valuesPrime[1000000];
extern int nbrFactorsFound;
extern enum eOutput pretty;
extern bool onlyEvaluate;
extern int degree;
extern int *ptrOrigPoly;
extern int degreeOrigPoly;
extern int polyInvCached;
extern int valuesIndex;
extern int poly5[COMPRESSED_POLY_MAX_LENGTH];
extern int polyNonRepeatedFactors[1000000];


typedef void (*powerCback)(int percentage);
void SetNumberToOne(/*@out@*/int *ptrValue1);
void PolyModularGcd(const int *arg1, int degree1, int *arg2, int degree2, int *gcd, int *degreeGcd);
void powerPolynomial(int* polyBase, int* polyMod, int polyDegree, const BigInteger* expon,
  int* polyPower, powerCback callback, int curMultip, int nbrMultip);
int getDegreePoly(const int *poly, int polyDegree);
void DividePolynomial(/*@in@*/int *pDividend, int dividendDegree, /*@in@*/int *pDivisor,
                      int divisorDegree, /*@out@*/int *ptrQuotient);
void multPolynomialModPoly(const int* polyFact1, const int* polyFact2,
  /*@out@*/int* polyProduct, int polyDegree, const int* polyMod);
void MultPolynomial(int degree1, int degree2, /*@in@*/int *factor1, /*@in@*/int *factor2);
void GetPolyInvParm(int polyDegree, /*@in@*/int *polyMod);
int ComputePolynomial(const char *input, int expo);
void OrigPolyFromMontgomeryToStandard(void);
void ConvertToMonic(int *poly, int polyDegree);
void SquareFreeFactorization(int polyDegree, int *poly, int expon);
int HenselLifting(struct sFactorInfo* factorInfo, bool compressPoly);
void polyToStandardNotation(int *nbr, int qtyNbrs);
void textErrorPol(char **pptrOutput, enum eExprErr rc);
void outputOriginalPolynomial(char **pptrOutput, int groupLen);
void outputPolynomialFactor(char** pptrOutput, int groupLength, const struct sFactorInfo* pstFactorInfo);
void DerPolynomial(int *ptrArgument);
void SubtractIntegerPolynomial(const int* minuend, const int* subtrahend, int* difference);
void PolynomialGcd(int *argF, int *argG, int *gcd);
int DivideIntegerPolynomial(int *pDividend, const int *pDivisor, enum eDivType type);
int getModPolynomial(int *polyMod, const int *poly, const BigInteger *content);
int *CopyPolynomial(int *ptrDest, const int *ptrSrc, int polyDegree);
int *CopyPolynomialFixedCoeffSize(int *ptrDest, const int *ptrSrc, int degree, int coeffSize);
void computePower(int expo);
void UncompressBigIntegerB(const int *ptrValues, BigInteger *bigint);
int numLimbs(const int *pLen);
void polyToStandardNotation(int *nbr, int qtyNbrs);
void polyToMontgomeryNotation(int *nbr, int qtyNbrs);
int *getContent(int *poly, BigInteger *content);
int *CopyPolyProduct(const int *ptrSrc, int *ptrDest, int polyDegree);
int FactorPolyOverIntegers(void);
void polyFactText(char* modText, char* polyText, int groupLength);
void showPowerX(char** pptrOutput, int polyDegree);
int FactorModularPolynomial(bool inputMontgomery);
int getNextPrimeNoDuplicatedFactors(int prime);
void FactorPolynomialModPrime(int prime);
void fftPolyMult(const int* factor1, const int* factor2, int* result, int len, int maxLen);
void multUsingInvPolynomial(/*@in@*/int* polyFact1, /*@in@*/int* polyFact2,
  /*@out@*/int* polyProduct,
  int polyDegree, /*@in@*/int* polyMod);
void SameDegreeFactorization(void);
int DivPolynomialExpr(int* ptrArgument1, const int* ptrArgument2, enum eDivType type);
#endif

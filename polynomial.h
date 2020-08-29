/*
This file is part of Alpertron Calculators.

Copyright 2015 Dario Alejandro Alpern

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
#ifndef __POLYNOMIAL_H
#define __POLYNOMIAL_H
#define MAX_DEGREE          1000

#define DEBUG_VANHOEIJ        0
#define DEBUG_HENSEL_LIFTING  0

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
  extern int lang;
#endif
extern BigInteger Quintic, Quartic, Cubic, Quadratic, Linear, Independent;
extern struct sFactorInfo factorInfo[MAX_DEGREE];
extern struct sFactorInfo factorInfoInteger[MAX_DEGREE];
extern BigInteger primeMod;              // p
extern int exponentMod;                  // k
extern BigInteger powerMod;              // p^k
extern int modulusIsZero;
extern BigInteger operand1, operand2, operand3, operand4, operand5;
extern int values[1000000];
extern int poly1[1000000];
extern int poly2[1000000];
extern int poly3[1000000];
extern int polyMultTemp[1000000];
extern int valuesPrime[1000000];
extern int nbrFactorsFound;
extern unsigned char pretty, onlyEvaluate;
extern int degree;
extern int *ptrOrigPoly;
extern int degreeOrigPoly;

typedef void (*powerCback)(int percentage);
void SetNumberToOne(/*@out@*/int *ptrValue1);
void PolyModularGcd(int *arg1, int degree1, int *arg2, int degree2, int *gcd, int *degreeGcd);
void powerPolynomial(int *polyBase, int *polyMod, int polyDegree, BigInteger *expon,
                     int *polyPower, powerCback cback);
int getDegreePoly(int *poly, int polyDegree);
void DividePolynomial(/*@in@*/int *pDividend, int dividendDegree, /*@in@*/int *pDivisor,
                      int divisorDegree, /*@out@*/int *ptrQuotient);
void multPolynomialModPoly(/*@in@*/int *polyFact1, /*@in@*/int *polyFact2, /*@out@*/int *polyProduct,
                    int polyDegree, /*@in@*/int *polyMod);
void MultPolynomial(int degree1, int degree2, /*@in@*/int *factor1, /*@in@*/int *factor2);
void GetPolyInvParm(int polyDegree, /*@in@*/int *polyMod);
int ComputePolynomial(char *input, int expo);
void OrigPolyFromMontgomeryToStandard(void);
void ConvertToMonic(int *poly, int polyDegree);
void SquareFreeFactorization(int polyDegree, int *poly, int expon);
int HenselLifting(struct sFactorInfo* factorInfo);
void polyToStandardNotation(int *nbr, int qtyNbrs);
void textErrorPol(char *output, enum eExprErr rc);
void outputOriginalPolynomial(char *ptrOutput, int groupLen);
void outputPolynomialFactor(char* ptrOutput, int groupLength, struct sFactorInfo* pstFactorInfo);
int DerPolynomial(int *ptrArgument);
void PolynomialGcd(int *argF, int *argG, int *gcd);
int DivideIntegerPolynomial(int *pDividend, int *pDivisor, enum eDivType type);
int getModPolynomial(int *polyMod, int *poly, BigInteger *content);
int *CopyPolynomial(int *ptrDest, int *ptrSrc, int degree);
int *CopyPolynomialFixedCoeffSize(int *ptrDest, int *ptrSrc, int degree, int coeffSize);
void computePower(int expo);
void UncompressBigIntegerB(int *ptrValues, BigInteger *bigint);
int numLimbs(int *pLen);
void polyToStandardNotation(int *nbr, int qtyNbrs);
void polyToMontgomeryNotation(int *nbr, int qtyNbrs);
int *getContent(int *poly, BigInteger *content);
int *CopyPolyProduct(int *ptrSrc, int *ptrDest, int degree);
int FactorPolyOverIntegers(void);
void polyFactText(char* modText, char* polyText, int groupLength);
void showPowerX(char** pptrOutput, int polyDegree);
int FactorModularPolynomial(int inputMontgomery);

#endif

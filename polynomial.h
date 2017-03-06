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

struct sFactorInfo
{
  int *ptr;
  int *ptrPolyLifted;
  int degree;
  int multiplicity;
  int expectedDegree;
};

extern int lang;
extern struct sFactorInfo factorInfo[MAX_DEGREE];
extern BigInteger primeMod;              // p
extern BigInteger powerMod;              // p^k
extern BigInteger operand1, operand2, operand3, operand4, operand5;
extern int values[1000000];
extern int poly1[1000000];
extern int poly2[1000000];
extern int poly3[1000000];
extern int polyMultTemp[1000000];
extern int valuesPrime[1000000];
extern int nbrFactorsFound;
extern unsigned char superscripts, onlyEvaluate;
extern int degree;
extern int *ptrOrigPoly;
extern int degreeOrigPoly;
extern char *output;

void SetNumberToOne(/*@out@*/int *ptrValue1);
void PolynomialGcd(int *arg1, int degree1, int *arg2, int degree2, int *gcd, int *degreeGcd);
void powerPolynomial(int *polyBase, int *polyMod, int polyDegree, BigInteger *expon, int *polyPower);
int getDegreePoly(int *poly, int polyDegree);
void DividePolynomial(/*@in@*/int *pDividend, int dividendDegree, /*@in@*/int *pDivisor, int divisorDegree, /*@out@*/int *ptrQuotient);
void multPolynomial(/*@in@*/int *polyFact1, /*@in@*/int *polyFact2, /*@out@*/int *polyProduct, int polyDegree, /*@in@*/int *polyMod);
void GetPolyInvParm(int polyDegree, /*@in@*/int *polyMod);
int ComputePolynomial(char *input, int expo);
void OrigPolyFromMontgomeryToStandard(void);
void ConvertToMonic(int *poly, int polyDegree);
void SquareFreeFactorization(int polyDegree, int *poly, int expon);
int HenselLifting(void);
void ToStandardNotation(int *nbr, int qtyNbrs);
void textErrorPol(char *output, enum eExprErr rc);
void outputPolynomial(char *ptrOutput, int groupLen);

#endif
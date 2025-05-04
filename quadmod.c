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
#ifndef _QUADMOD_H
#define _QUADMOD_H

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bignbr.h"
#include "expression.h"
#include "factor.h"
#include "quadmodLL.h"

static BigInteger discriminant;
static BigInteger sqrtDiscriminant;
static BigInteger Aux0;
static BigInteger Aux1;
static BigInteger Aux2;
static BigInteger GcdAll;
static BigInteger ValNn;

BigInteger ValA;
BigInteger ValB;
BigInteger ValC;
BigInteger ValN;
static int SolNbr;
#ifdef __ANDROID__
extern char *ptrOutput;
#else
char *ptrOutput;
#endif
extern int factorsMod[20000];

static int Show(const BigInteger *num, const char *str, int t)
{
  if (!BigIntIsZero(num))
  {     // num is not zero.
    if (((t & 1) != 0) && (num->sign == SIGN_POSITIVE))
    {
      copyStr(&ptrOutput, " +");
    }
    if (num->sign == SIGN_NEGATIVE)
    {
      copyStr(&ptrOutput, " &minus;");
    }
    if ((num->nbrLimbs != 1) || (num->limbs[0].x != 1))
    {    // num is not 1 or -1.
      *ptrOutput = ' ';
      ptrOutput++;
      Bin2Dec(&ptrOutput, num->limbs, num->nbrLimbs, groupLen);
    }
    copyStr(&ptrOutput, str);
    return t | 1;
  }
  return t;
}

void Show1(const BigInteger *num, int t)
{
  int u = Show(num, "", t);
  if (((u & 1) == 0) || ((num->nbrLimbs == 1) && (num->limbs[0].x == 1)))
  {
    *ptrOutput = ' ';
    ptrOutput++;
    CopyBigInt(&Aux1, num);
    Aux1.sign = SIGN_POSITIVE;
    BigInteger2Dec(&ptrOutput, &Aux1, groupLen);
  }
}

void Solution(BigInteger *value)
{
  SolNbr++;
  copyStr(&ptrOutput, "<li>x = ");
  BigInteger2Dec(&ptrOutput, value, groupLen);
  copyStr(&ptrOutput, "</li>");
}

// Solve Ax^2 + Bx + C = 0 in integers.
static void SolveIntegerEquation(void)
{
  if (BigIntIsZero(&ValA))
  {      // Not Quadratic Equation
    if (BigIntIsZero(&ValB))
    {    // Constant Equation
      if (BigIntIsZero(&ValC))
      {  // 0 = 0
        copyStr(&ptrOutput, "<p>The equation is satisfied by any integer <var>x</var>.</p>");
      }
      else
      {
        return;
      }
    }
    else
    {        // Linear Equation: ValB * x + ValC = 0.
      (void)BigIntRemainder(&ValC, &ValB, &Aux0);
      if (BigIntIsZero(&Aux0))
      {      // ValC is multiple of ValB: solution is -ValC / ValB.
        (void)BigIntDivide(&ValC, &ValB, &Aux0);
        BigIntNegate(&Aux0, &Aux0);
        Solution(&Aux0);
      }
      else
      {      // No solutions on integers.
        return;
      }
    }
  }
  else
  {                       // Quadratic Equation
    // The discriminant is ValB * ValB - 4 * ValA * ValC.
    (void)BigIntMultiply(&ValB, &ValB, &Aux0);
    (void)BigIntMultiply(&ValA, &ValC, &Aux1);
    multint(&Aux1, &Aux1, 4);
    BigIntSubt(&Aux0, &Aux1, &discriminant);
    if (discriminant.sign == SIGN_NEGATIVE)
    {        // No integer solutions.
      return;
    }
    // Set sqrtDiscriminant to square root of discriminant.
    squareRoot(discriminant.limbs, sqrtDiscriminant.limbs, discriminant.nbrLimbs, &sqrtDiscriminant.nbrLimbs);
    sqrtDiscriminant.sign = SIGN_POSITIVE;
    (void)BigIntMultiply(&sqrtDiscriminant, &sqrtDiscriminant, &Aux0);
    BigIntSubt(&Aux0, &discriminant, &Aux0);
    if (!BigIntIsZero(&Aux0))
    {  // discriminant has no integer square root.
      return;
    }
    multint(&Aux0, &ValA, 2);    // Denominator
    BigIntSubt(&sqrtDiscriminant, &ValB, &Aux1);
    (void)BigIntRemainder(&Aux1, &Aux0, &Aux2);
    if (BigIntIsZero(&Aux2))
    {      // (sqrtDiscriminant-ValB)/(2*ValA) is integer: it is a solution.
      (void)BigIntDivide(&Aux1, &Aux0, &Aux2);
      Solution(&Aux2);
    }
    BigIntNegate(&sqrtDiscriminant, &sqrtDiscriminant);
    BigIntSubt(&sqrtDiscriminant, &ValB, &Aux1);
    (void)BigIntRemainder(&Aux1, &Aux0, &Aux2);
    if (BigIntIsZero(&Aux2))
    {      // (-sqrtDiscriminant-ValB)/(2*ValA) is integer: it is a solution.
      (void)BigIntDivide(&Aux1, &Aux0, &Aux2);
      Solution(&Aux2);
    }
  }
  SolNbr = 1;
}

void textErrorQuadMod(char **pptrOutput, enum eExprErr rc)
{
  if (rc == EXPR_MODULUS_MUST_BE_NONNEGATIVE)
  {
    copyStr(pptrOutput, lang ? "No debe ser negativo" :
      "Must not be negative");
  }
  else
  {
    textError(pptrOutput, rc);
  }
}

static void ModulusIsNotZero(void)
{
  // Modulus is not zero.
  ValN.sign = SIGN_POSITIVE;
  (void)BigIntRemainder(&ValA, &ValN, &ValA);
  (void)BigIntRemainder(&ValB, &ValN, &ValB);
  (void)BigIntRemainder(&ValC, &ValN, &ValC);
  BigIntGcd(&ValA, &ValB, &Aux0);
  Aux0.sign = SIGN_POSITIVE;
  BigIntGcd(&ValN, &Aux0, &GcdAll);
  GcdAll.sign = SIGN_POSITIVE;
  (void)BigIntRemainder(&ValC, &GcdAll, &Aux0);
  if (!BigIntIsZero(&Aux0))
  {  // ValC must be multiple of gcd(ValA, ValB).
     // Otherwise go out because there are no solutions.
    return;
  }
  BigIntGcd(&ValN, &GcdAll, &Aux0);
  CopyBigInt(&GcdAll, &Aux0);
  // Divide all coefficients by gcd(ValA, ValB).
  (void)BigIntDivide(&ValA, &GcdAll, &ValA);
  (void)BigIntDivide(&ValB, &GcdAll, &ValB);
  (void)BigIntDivide(&ValC, &GcdAll, &ValC);
  (void)BigIntDivide(&ValN, &GcdAll, &ValN);
  CopyBigInt(&ValNn, &ValN);
  if ((ValNn.nbrLimbs == 1) && (ValNn.limbs[0].x == 1))
  {     // All values from 0 to GcdAll - 1 are solutions.
    if ((GcdAll.nbrLimbs > 1) || (GcdAll.limbs[0].x > 5))
    {
      copyStr(&ptrOutput, "<p>All values of <var>x</var> between 0 and ");
      addbigint(&GcdAll, -1);
      BigInteger2Dec(&ptrOutput, &GcdAll, groupLen);
      copyStr(&ptrOutput, " are solutions.</p>");
    }
    else
    {
      for (int ctr = 0; ctr < GcdAll.limbs[0].x; ctr++)
      {
        intToBigInteger(&Aux0, ctr);
        Solution(&Aux0);
      }
    }
    return;
  }
  SetCallbacksForSolveEquation(&Solution, NULL,
    NULL);
  SolveEquation(&ValA, &ValB, &ValC, &ValN, &GcdAll, &ValNn);
}

void quadmodText(const char *quadrText, const char *linearText, const char *constText,
  const char *modText, int groupLength)
{
  groupLen = groupLength;
  char *ptrBeginSol;
  enum eExprErr rc;
  ptrOutput = output;
  *ptrOutput = '2';  // Prefix for sending output to screen.
  ptrOutput++;
  rc = ComputeExpression(quadrText, &ValA);
  if (rc != EXPR_OK)
  {
    copyStr(&ptrOutput, lang ? "<p>Coeficiente cuadrático: ": "<p>Quadratic coefficient: ");
    textErrorQuadMod(&ptrOutput, rc);
    copyStr(&ptrOutput, "</p>");
  }
  rc = ComputeExpression(linearText, &ValB);
  if (rc != EXPR_OK)
  {
    copyStr(&ptrOutput, lang ? "<p>Coeficiente lineal: " : "<p>Linear coefficient: ");
    textErrorQuadMod(&ptrOutput, rc);
    copyStr(&ptrOutput, "</p>");
  }
  rc = ComputeExpression(constText, &ValC);
  if (rc != EXPR_OK)
  {
    copyStr(&ptrOutput, lang ? "<p>Término independiente: " : "<p>Constant coefficient: ");
    textErrorQuadMod(&ptrOutput, rc);
    copyStr(&ptrOutput, "</p>");
  }
  rc = ComputeExpression(modText, &ValN);
  if ((rc == EXPR_OK) && (ValN.sign == SIGN_NEGATIVE))
  {
    rc = EXPR_MODULUS_MUST_BE_NONNEGATIVE;
  }
  if (rc != EXPR_OK)
  {
    copyStr(&ptrOutput, lang ? "Módulo: " : "Modulus: ");
    textErrorQuadMod(&ptrOutput, rc);
    copyStr(&ptrOutput, "</p>");
  }
  if (ptrOutput == &output[1])
  {    // No errors found.
    copyStr(&ptrOutput, "<p>");
    int u = Show(&ValA, " x&sup2;", 2);
    u = Show(&ValB, " x", u);
    Show1(&ValC, u);
    copyStr(&ptrOutput, " &equiv; 0 (mod ");
    BigInteger2Dec(&ptrOutput, &ValN, groupLen);
    copyStr(&ptrOutput, ")</p>");
    SolNbr = 0;
    ptrBeginSol = ptrOutput;
    copyStr(&ptrOutput, "<ol>");
    if (BigIntIsZero(&ValN))
    {        // Mod zero => Equation in integer numbers
      SolveIntegerEquation();
    }
    else
    {
      ModulusIsNotZero();
    }
    if (SolNbr == 0)
    {
      ptrOutput = ptrBeginSol;
      copyStr(&ptrOutput, lang? "<p>No hay soluciones.</p>": "<p>There are no solutions.</p>");
    }
    else
    {
      copyStr(&ptrOutput, "</ol>");
    }
  }
  copyStr(&ptrOutput, "<p>");
  copyStr(&ptrOutput, lang ? COPYRIGHT_SPANISH: COPYRIGHT_ENGLISH);
  copyStr(&ptrOutput, "</p>");
}

#if defined(__EMSCRIPTEN__) && !defined(_MSC_VER)
#ifdef __ANDROID__
EXTERNALIZE void doWorkQuadMod(void)
#else
EXTERNALIZE void doWork(void)
#endif
{
  char *ptrData = inputString;
  char* ptrQuadrCoeff;
  char* ptrLinearCoeff;
  char* ptrConstCoeff;
  char* ptrMod;
  int groupLength = 0;
  while (*ptrData != ',')
  {
    groupLength = (groupLength * 10) + (*ptrData - '0');
    ptrData++;
  }
  ptrData++;                    // Skip comma.
#ifndef lang  
  int flags = *ptrData;
  lang = ((flags & 1)? true: false);
#endif
  ptrQuadrCoeff = ptrData + 2;  // Skip flags and comma.
  ptrLinearCoeff = ptrQuadrCoeff + (int)strlen(ptrQuadrCoeff) + 1;
  ptrConstCoeff = ptrLinearCoeff + (int)strlen(ptrLinearCoeff) + 1;
  ptrMod = ptrConstCoeff + (int)strlen(ptrConstCoeff) + 1;
  quadmodText(ptrQuadrCoeff, ptrLinearCoeff, ptrConstCoeff, ptrMod, groupLength);
  databack(output);
}
#endif
#endif

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

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include "bignbr.h"
#include "polynomial.h"
#include "expression.h"

#define TOKEN_GCD    34
#define TOKEN_LCM    35
#define TOKEN_DER    36
static enum eExprErr GcdPolynomialExpr(int* ptrArgument1, int* ptrArgument2);
static enum eExprErr LcmPolynomialExpr(int* ptrArgument1, int* ptrArgument2);
static enum eExprErr RandomPolynomialExpr(const int* pMinDegree, const int* pMaxDegree,
  const int* pMinCoeff, const int* pMaxCoeff, int* randomPoly);
static int* stackValues[STACK_OPER_SIZE];
extern int polyA[1000000];
extern int polyB[1000000];
extern int polyC[1000000];
extern int polyD[1000000];
static int LastAnswerPoly[1000000];
extern bool onlyEvaluate;

struct sFuncOperExpr stFuncOperPolyExpr[] =
{
  // First section: functions
  {"ANS", TOKEN_ANS + NO_PARMS, 0},
  {"GCD", TOKEN_GCD + MANY_PARMS, 0},
  {"LCM", TOKEN_LCM + MANY_PARMS, 0},
  {"DER", TOKEN_DER + ONE_PARM, 0},
  {"RANDOM", TOKEN_RANDOM + FOUR_PARMS, 0},
  {NULL, 0},
  // Second section: functions written at right of argument.
  {NULL, 0},
  // Third section: unary operators.
  {"-", OPER_UNARY_MINUS, 3},
  {NULL, 0},
  // Fourth section: binary operators.
  {"**", OPER_POWER, 1}, // This must be located before multiplication operator.
  {"+", OPER_ADD, 3},
  {"-", OPER_SUBT, 3},
  {"*", OPER_MULTIPLY, 2},
  {"%", OPER_REMAINDER, 2},
  {"/", OPER_DIVIDE, 2},
  {"^", OPER_POWER, 1},
  {"=", OPER_EQUAL, 4},
  {NULL, 0},
};

static int NegatePolynomialExpr(int* ptrArgument)
{
  int* ptrValue1 = ptrArgument;
  if (!modulusIsZero)
  {
    NumberLength = powerMod.nbrLimbs;
  }
  int val = *ptrValue1;
  if (val <= 0)
  {          // Monomial
    if ((*(ptrValue1 + 1) != 1) || (*(ptrValue1 + 2) != 0))
    {        // Coefficient is not zero
      if (modulusIsZero)
      {
        *(ptrValue1 + 1) = -*(ptrValue1 + 1);    // Negate number.
      }
      else
      {
        IntArray2BigInteger(ptrValue1 + 1, &operand1);
        BigIntSubt(&powerMod, &operand1, &operand1);
        BigInteger2IntArray(ptrValue1 + 1, &operand1);
      }
    }
    return EXPR_OK;
  }
  else
  {          // Polynomial. Use poly1 as temporary polynomial.
    int* ptrValue2 = poly1;
    ptrValue1++;
    for (int currentDegree = val; currentDegree >= 0; currentDegree--)
    {
      if ((*ptrValue1 == 1) && (*(ptrValue1 + 1) == 0))
      {          // If value is zero, it does not have to be changed.
        *ptrValue2 = 1;
        ptrValue2++;
        *ptrValue2 = 0;
        ptrValue2++;
        ptrValue1 += 2;              // Point to next coefficient.
      }
      else if (modulusIsZero)
      {                              // Integer polynomial.
        int nbrLimbs = *ptrValue1;
        *ptrValue1 = -nbrLimbs;      // Negate number.
        if (nbrLimbs < 0)
        {
          nbrLimbs = -nbrLimbs;
        }
        ptrValue1 += nbrLimbs + 1;   // Point to next coefficient.
      }
      else
      {                              // Polynomial modulo prime.
        IntArray2BigInteger(ptrValue1, &operand1);
        BigIntSubt(&powerMod, &operand1, &operand1);
        BigInteger2IntArray(ptrValue2, &operand1);
        ptrValue1 += *ptrValue1 + 1; // Point to next coefficient.
        ptrValue2 += *ptrValue2 + 1;
      }
    }
    if (!modulusIsZero)
    {
      if ((ptrArgument + 1 + (ptrValue2 - &poly1[0])) > (values + sizeof(values) / sizeof(values[0])))
      {
        return EXPR_OUT_OF_MEMORY;
      }
      (void)memcpy(ptrArgument + 1, poly1, (ptrValue2 - &poly1[0]) * sizeof(int));
    }
  }
  return EXPR_OK;
}

void UncompressBigIntegerB(const int* ptrValues, BigInteger* bigint)
{
  if (modulusIsZero)
  {
    NumberLength = numLimbs(ptrValues);
  }
  assert(NumberLength >= 1);
  IntArray2BigInteger(ptrValues, bigint);
}

// Add two polynomials starting on buffers pointed by ptrArgument1
// and ptrArgument2. Use poly1 as a temporary buffer to hold the sum.
static int AddPolynomialExpr(int* ptrArgument1, int* ptrArgument2)
{
  int* ptrValue1;
  const int* ptrValue2;
  int currentDegree;
  int degreeMin;
  int degreeMax;
  int degreePoly = 0;
  int degree1 = *ptrArgument1;
  int degree2 = *ptrArgument2;
  size_t diffPtrs;
  if (!modulusIsZero)
  {
    NumberLength = powerMod.nbrLimbs;
  }
  if (degree1 <= 0)
  {
    if (degree1 == degree2)
    {      // Sum of two monomials of same degree.
      UncompressBigIntegerB(ptrArgument1 + 1, &operand1);
      UncompressBigIntegerB(ptrArgument2 + 1, &operand2);
      if (modulusIsZero)
      {
        BigIntAdd(&operand1, &operand2, &operand1);
        NumberLength = operand1.nbrLimbs;
        BigInteger2IntArray(ptrArgument1 + 1, &operand1);
      }
      else
      {
        AddBigNbrModN(operand1.limbs, operand2.limbs, operand1.limbs, powerMod.limbs, powerMod.nbrLimbs);
        BigInteger2IntArray(ptrArgument1 + 1, &operand1);
      }
      if (BigIntIsZero(&operand1))
      {                     // Sum is zero: set degree to zero.
        *ptrArgument1 = 0;
      }
      diffPtrs = ptrArgument1 - &values[0];
      valuesIndex = (int)diffPtrs + 2 + NumberLength;
      return EXPR_OK;
    }
    if (degree2 <= 0)
    {     // Sum of two monomials of different degree. Generate polynomial.
      degree1 = -degree1;   // Convert to positive.
      degree2 = -degree2;
      if (degree1 > degree2)
      {
        degreeMax = degree1;
        degreeMin = degree2;
        UncompressBigIntegerB(ptrArgument1 + 1, &operand2);
        UncompressBigIntegerB(ptrArgument2 + 1, &operand1);
      }
      else
      {
        degreeMax = degree2;
        degreeMin = degree1;
        UncompressBigIntegerB(ptrArgument1 + 1, &operand1);
        UncompressBigIntegerB(ptrArgument2 + 1, &operand2);
      }
      ptrValue1 = ptrArgument1;
      *ptrValue1 = degreeMax;
      ptrValue1++;
      for (currentDegree = 0; currentDegree < degreeMin; currentDegree++)
      {
        *ptrValue1 = 1;  // Number of limbs
        ptrValue1++;
        *ptrValue1 = 0;  // Value = zero.
        ptrValue1++;
      }
      BigInteger2IntArray(ptrValue1, &operand1);
      ptrValue1 += numLimbs(ptrValue1);
      ptrValue1++;
      currentDegree++;
      for (; currentDegree < degreeMax; currentDegree++)
      {
        *ptrValue1 = 1;  // Number of limbs
        ptrValue1++;
        *ptrValue1 = 0;  // Value = zero.
        ptrValue1++;
      }
      NumberLength = operand2.nbrLimbs;
      BigInteger2IntArray(ptrValue1, &operand2);
      diffPtrs = ptrValue1 - &values[0];
      valuesIndex = (int)diffPtrs + 1 + numLimbs(ptrValue1);
      return EXPR_OK;
    }
  }
  if ((degree1 > 0) && (degree2 > 0))
  {           // Sum of two polynomials.
    const int* ptrPolyMin;
    const int* ptrPolyMax;
    if (degree1 > degree2)
    {
      degreeMin = degree2;
      degreeMax = degree1;
      ptrPolyMin = ptrArgument2 + 1;
      ptrPolyMax = ptrArgument1 + 1;
    }
    else
    {
      degreeMin = degree1;
      degreeMax = degree2;
      ptrPolyMin = ptrArgument1 + 1;
      ptrPolyMax = ptrArgument2 + 1;
    }
    ptrValue1 = poly1;    // Point to temporary polynomial sum.
    *ptrValue1 = degreeMax;
    ptrValue1++;
    for (currentDegree = 0; currentDegree <= degreeMin; currentDegree++)
    {
      UncompressBigIntegerB(ptrPolyMin, &operand1);
      UncompressBigIntegerB(ptrPolyMax, &operand2);
      if (modulusIsZero)
      {
        BigIntAdd(&operand1, &operand2, &operand1);
        NumberLength = operand1.nbrLimbs;
      }
      else
      {
        AddBigNbrModN(operand1.limbs, operand2.limbs, operand1.limbs, powerMod.limbs, powerMod.nbrLimbs);
      }
      BigInteger2IntArray(ptrValue1, &operand1);
      ptrPolyMin += numLimbs(ptrPolyMin);
      ptrPolyMin++;
      ptrPolyMax += numLimbs(ptrPolyMax);
      ptrPolyMax++;
      ptrValue1 += numLimbs(ptrValue1);
      ptrValue1++;
    }
    for (; currentDegree <= degreeMax; currentDegree++)
    {
      int lenBytes = (1 + numLimbs(ptrPolyMax)) * (int)sizeof(int);
      (void)memcpy(ptrValue1, ptrPolyMax, lenBytes);
      ptrPolyMax += numLimbs(ptrPolyMax);
      ptrPolyMax++;
      ptrValue1 += numLimbs(ptrValue1);
      ptrValue1++;
    }
    (void)memcpy(ptrArgument1, poly1, (ptrValue1 - &poly1[0]) * sizeof(int));
    degreePoly = degreeMax;   // New degree of polynomial.
  }
  else
  {                 // Sum of polynomial and monomial.
    int degreeMono;
    if (degree1 < 0)
    {
      UncompressBigIntegerB(ptrArgument1 + 1, &operand1);
      (void)memmove(ptrArgument1 + 1, ptrArgument1 + 3 + numLimbs(ptrArgument1 + 1), (&values[valuesIndex] - ptrArgument2) * sizeof(int));
      degreeMono = -degree1;
      degreePoly = degree2;
      ptrValue1 = ptrArgument2 + 1;
    }
    else
    {
      UncompressBigIntegerB(ptrArgument2 + 1, &operand1);
      degreeMono = -degree2;
      degreePoly = degree1;
      ptrValue1 = ptrArgument1 + 1;
    }
    if (degreeMono > degreePoly)
    {
      int* ptrDest = ptrArgument1 + 1;
      *ptrArgument1 = degreeMono;
      for (currentDegree = 0; currentDegree <= degreePoly; currentDegree++)
      {
        int numLen = 1 + numLimbs(ptrValue1);
        int lenBytes = numLen * (int)sizeof(int);
        (void)memcpy(ptrDest, ptrValue1, lenBytes);
        ptrValue1 += numLen;
        ptrDest += numLen;
      }
      for (; currentDegree < degreeMono; currentDegree++)
      {
        *ptrDest = 1;  // Number of limbs
        ptrDest++;
        *ptrDest = 0;  // Value = zero.
        ptrDest++;
      }
      BigInteger2IntArray(ptrDest, &operand1);
      degreePoly = degreeMono;   // Set new degree.
    }
    else
    {        // Degree of polynomial greater than degree of monomial.
      int differenceOfDegrees;
      int nbrLimbsSum;
      int lenBytes;

      *ptrArgument1 = degreePoly;
      ptrValue1 = ptrArgument1 + 1;
      for (currentDegree = 0; currentDegree < degreeMono; currentDegree++)
      {
        ptrValue1 += numLimbs(ptrValue1);
        ptrValue1++;
      }
      UncompressBigIntegerB(ptrValue1, &operand2);
      if (modulusIsZero)
      {
        BigIntAdd(&operand1, &operand2, &operand1);
        NumberLength = operand1.nbrLimbs;
      }
      else
      {
        AddBigNbrModN(operand1.limbs, operand2.limbs, operand1.limbs, powerMod.limbs, powerMod.nbrLimbs);
      }
      BigInteger2IntArray(poly1, &operand1);
      nbrLimbsSum = numLimbs(&poly1[0]);
      differenceOfDegrees = nbrLimbsSum - numLimbs(ptrValue1);
      diffPtrs = ptrArgument2 - ptrValue1;
      if (differenceOfDegrees > 0)
      {    // New coefficient is greater than old coefficient.
        lenBytes = (int)diffPtrs * (int)sizeof(int);
        (void)memmove(ptrValue1 + differenceOfDegrees, ptrValue1, lenBytes);
      }
      if (differenceOfDegrees < 0)
      {
        lenBytes = ((int)diffPtrs - differenceOfDegrees) * (int)sizeof(int);
        (void)memmove(ptrValue1, ptrValue1 - differenceOfDegrees, lenBytes);
      }
      lenBytes = (1 + nbrLimbsSum) * (int)sizeof(int);
      (void)memcpy(ptrValue1, poly1, lenBytes);
    }
  }
  // Reduce degree if leading coefficient is zero.
  degreeMax = 0;
  ptrValue1 = ptrArgument1 + 1;
  ptrValue2 = ptrValue1;
  for (currentDegree = 0; currentDegree <= degreePoly; currentDegree++)
  {
    if ((*ptrValue1 != 1) || (*(ptrValue1 + 1) != 0))
    {                    // Coefficient is not zero
      degreeMax = currentDegree;
      // Store point to coefficient not zero of maximum degree.
      ptrValue2 = ptrValue1;
    }
    ptrValue1 += numLimbs(ptrValue1);
    ptrValue1++;
  }
  *ptrArgument1 = degreeMax;
  diffPtrs = ptrValue2 - &values[0];
  valuesIndex = (int)diffPtrs + 1 + numLimbs(ptrValue2);
  return EXPR_OK;
}

int* CopyPolyProduct(const int* ptrSrc, int* ptrDest, int polyDegree)
{
  const int* ptrValueSrc = ptrSrc;
  int* ptrValueDest = ptrDest;
  for (int currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {
    int nbrLength = 1 + numLimbs(ptrValueSrc);
    int lenBytes = nbrLength * (int)sizeof(int);
    (void)memcpy(ptrValueDest, ptrValueSrc, lenBytes);
    ptrValueSrc += nbrLength;
    ptrValueDest += nbrLength;
  }
  return ptrValueDest;
}

static enum eExprErr MultPolynomialExpr(int* ptrArgument1, const int* ptrArgument2)
{
  int degreeMono;
  int degreePoly;
  int nbrLimbs = powerMod.nbrLimbs + 1;
  int* ptrValue1;
  int* ptrValue2;
  int currentDegree;
  int degree1 = *ptrArgument1;
  int degree2 = *ptrArgument2;
  const int* ptrValueSrc;

  if (!modulusIsZero)
  {
    NumberLength = powerMod.nbrLimbs;
  }
  if ((degree1 <= 0) && (degree2 <= 0))
  {        // Product of two monomials.
    size_t diffPtrs;
    if ((degree1 + degree2) < -MAX_DEGREE)
    {
      return EXPR_DEGREE_TOO_HIGH;
    }
    *ptrArgument1 = degree1 + degree2;
    UncompressBigIntegerB(ptrArgument1 + 1, &operand1);
    UncompressBigIntegerB(ptrArgument2 + 1, &operand2);
    if (modulusIsZero)
    {
      (void)BigIntMultiply(&operand1, &operand2, &operand1);
      NumberLength = operand1.nbrLimbs;
    }
    else
    {
      modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    }
    BigInteger2IntArray(ptrArgument1 + 1, &operand1);
    diffPtrs = ptrArgument1 - &values[0];
    valuesIndex = (int)diffPtrs + 2 + *(ptrArgument1 + 1);
    return EXPR_OK;
  }
  if ((degree1 > 0) && (degree2 > 0))
  {        // Product of two polynomials.
    size_t diffPtrs;
    if ((degree1 + degree2) > MAX_DEGREE)
    {
      return EXPR_DEGREE_TOO_HIGH;
    }
    *ptrArgument1 = degree1 + degree2;
    if (modulusIsZero)
    {
      (void)CopyPolyProduct(ptrArgument1 + 1, poly1, degree1);  // Copy first factor to poly1
      (void)CopyPolyProduct(ptrArgument2 + 1, poly2, degree2);  // Copy second factor to poly2
    }
    else
    {
      // Copy first factor to poly1
      ptrValue1 = ptrArgument1 + 1;
      ptrValue2 = poly1;
      for (currentDegree = 0; currentDegree <= degree1; currentDegree++)
      {
        int lenLimbs = 1 + *ptrValue1;
        int lenBytes = lenLimbs * (int)sizeof(int);
        (void)memcpy(ptrValue2, ptrValue1, lenBytes);
        ptrValue1 += lenLimbs;
        ptrValue2 += nbrLimbs;
      }
      // Copy second factor to poly2
      ptrValueSrc = ptrArgument2 + 1;
      ptrValue2 = poly2;
      for (currentDegree = 0; currentDegree <= degree2; currentDegree++)
      {
        int lenLimbs = 1 + *ptrValueSrc;
        int lenBytes = lenLimbs * (int)sizeof(int);
        (void)memcpy(ptrValue2, ptrValueSrc, lenBytes);
        ptrValueSrc += lenLimbs;
        ptrValue2 += nbrLimbs;
      }
    }
    // Perform multiplication storing the result in polyMultTemp
    MultPolynomial(degree1, degree2, poly1, poly2);
    // Move product back to values stack.
    if (modulusIsZero)
    {
      ptrValue1 = CopyPolyProduct(polyMultTemp, ptrArgument1 + 1, degree1 + degree2);
    }
    else
    {
      ptrValue1 = ptrArgument1 + 1;
      ptrValue2 = polyMultTemp;
      for (currentDegree = 0; currentDegree <= (degree1 + degree2); currentDegree++)
      {
        int lenLimbs = 1 + *ptrValue2;
        int lenBytes = lenLimbs * (int)sizeof(int);
        (void)memcpy(ptrValue1, ptrValue2, lenBytes);
        ptrValue1 += lenLimbs;
        ptrValue2 += nbrLimbs;
      }
    }
    diffPtrs = ptrValue1 - &values[0];
    valuesIndex = (int)diffPtrs;
    return EXPR_OK;
  }
  // Product of monomial and polynomial.
  if (degree1 <= 0)
  {      // First factor is monomial.
    degreeMono = -degree1;
    degreePoly = degree2;
    // Point to first coefficient of polynomial.
    ptrValueSrc = ptrArgument2 + 1;
    // Get coefficient of monomial.
    UncompressBigIntegerB(ptrArgument1 + 1, &operand1);
    if ((degreeMono + degreePoly) > MAX_DEGREE)
    {
      return EXPR_DEGREE_TOO_HIGH;
    }
    // Multiply all coefficients of polynomial by the coefficient
    // of monomial storing the resulting polynomial on poly1.
    ptrValue2 = poly1;    // Initialize pointer to product.
    for (currentDegree = 0; currentDegree < degreeMono; currentDegree++)
    {
      *ptrValue2 = 1;
      ptrValue2++;
      *ptrValue2 = 0;
      ptrValue2++;
    }
    for (currentDegree = 0; currentDegree <= degreePoly; currentDegree++)
    {
      UncompressBigIntegerB(ptrValueSrc, &operand2);
      ptrValueSrc += numLimbs(ptrValueSrc);
      ptrValueSrc++;
      if (modulusIsZero)
      {
        (void)BigIntMultiply(&operand1, &operand2, &operand3);
        NumberLength = operand3.nbrLimbs;
      }
      else
      {
        modmult(operand1.limbs, operand2.limbs, operand3.limbs);
      }
      BigInteger2IntArray(ptrValue2, &operand3);
      ptrValue2 += numLimbs(ptrValue2);
      ptrValue2++;
    }
  }
  else
  {      // Second factor is monomial.
    degreeMono = -degree2;
    degreePoly = degree1;
    // Point to first coefficient of polynomial.
    ptrValue1 = ptrArgument1 + 1;
    // Get coefficient of monomial.
    UncompressBigIntegerB(ptrArgument2 + 1, &operand1);
    if ((degreeMono + degreePoly) > MAX_DEGREE)
    {
      return EXPR_DEGREE_TOO_HIGH;
    }
    // Multiply all coefficients of polynomial by the coefficient
    // of monomial storing the resulting polynomial on poly1.
    ptrValue2 = poly1;    // Initialize pointer to product.
    for (currentDegree = 0; currentDegree < degreeMono; currentDegree++)
    {
      *ptrValue2 = 1;
      ptrValue2++;
      *ptrValue2 = 0;
      ptrValue2++;
    }
    for (currentDegree = 0; currentDegree <= degreePoly; currentDegree++)
    {
      UncompressBigIntegerB(ptrValue1, &operand2);
      ptrValue1 += numLimbs(ptrValue1);
      ptrValue1++;
      if (modulusIsZero)
      {
        (void)BigIntMultiply(&operand1, &operand2, &operand3);
        NumberLength = operand3.nbrLimbs;
      }
      else
      {
        modmult(operand1.limbs, operand2.limbs, operand3.limbs);
      }
      BigInteger2IntArray(ptrValue2, &operand3);
      ptrValue2 += numLimbs(ptrValue2);
      ptrValue2++;
    }
  }
  *ptrArgument1 = degreeMono + degreePoly;
  ptrValue1 = ptrArgument1 + 1;
  (void)memcpy(ptrValue1, poly1, (ptrValue2 - &poly1[0]) * sizeof(int));
  valuesIndex = (int)(ptrValue1 - &values[0] + (ptrValue2 - &poly1[0]));
  return EXPR_OK;
}

void SetNumberToOne(/*@out@*/int* ptrValue1)
{
  int* ptrToValue1 = ptrValue1;
  const limb* destLimb;
  *ptrToValue1 = NumberLengthR1;
  ptrToValue1++;
  destLimb = MontgomeryMultR1;
  for (int ctr = 0; ctr < NumberLengthR1; ctr++)
  {
    *ptrToValue1 = destLimb->x;
    ptrToValue1++;
    destLimb++;
  }
}

int* CopyPolynomial(int* dest, const int* src, int polyDegree)
{
  int* ptrDest = dest;
  const int* ptrSrc = src;
  int currentDegree;
  if (polyDegree >= 0)
  {          // Polynomial
    currentDegree = 0;
  }
  else
  {          // Monomial
    currentDegree = polyDegree;
  }
  for (; currentDegree <= polyDegree; currentDegree++)
  {
    int numLength = numLimbs(ptrSrc) + 1;
    int lenBytes = numLength * (int)sizeof(int);
    (void)memcpy(ptrDest, ptrSrc, lenBytes);
    ptrSrc += numLength;
    ptrDest += numLength;
  }
  return ptrDest;
}

int* CopyPolynomialFixedCoeffSize(int* dest, const int* src, int polyDegree, int coeffSize)
{
  int* ptrDest = dest;
  const int* ptrSrc = src;
  for (int currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {
    int numLength = numLimbs(ptrSrc) + 1;
    int lenBytes = numLength * (int)sizeof(int);
    (void)memcpy(ptrDest, ptrSrc, lenBytes);
    ptrSrc += coeffSize;
    ptrDest += coeffSize;
  }
  return ptrDest;
}

static int PowerPolynomialExpr(int* ptrArgument1, int expon)
{
  limb exponLimb;
  int degreePower;
  int* ptrValue1;
  int* ptrValue2;
  int currentDegree;
  int nbrLimbs = powerMod.nbrLimbs + 1;
  int degreeBase = *ptrArgument1;
  size_t diffPtrs;
  if (!modulusIsZero)
  {
    NumberLength = powerMod.nbrLimbs;
  }
  if (degreeBase <= 0)
  {              // Monomial.
    if ((-degreeBase * expon) > MAX_DEGREE)
    {
      return EXPR_DEGREE_TOO_HIGH;
    }
    *ptrArgument1 = degreeBase * expon;
    UncompressBigIntegerB(ptrArgument1 + 1, &operand1);
    if (expon < 0)
    {
      return EXPR_EXPONENT_NEGATIVE;
    }
    if (modulusIsZero)
    {
      (void)BigIntPowerIntExp(&operand1, expon, &operand2);
      NumberLength = operand2.nbrLimbs;
    }
    else
    {
      exponLimb.x = expon;
      modPowLimb(operand1.limbs, &exponLimb, operand2.limbs);
      operand2.sign = SIGN_POSITIVE;
    }
    BigInteger2IntArray(ptrArgument1 + 1, &operand2);
    diffPtrs = ptrArgument1 - &values[0];
    valuesIndex = (int)diffPtrs + 2 + numLimbs(ptrArgument1 + 1);
    return EXPR_OK;
  }
  // Polynomial.
  ptrValue1 = ptrArgument1 + 1;
  ptrValue2 = poly1;
  if (!modulusIsZero)
  {
    // Copy base to poly1
    for (currentDegree = 0; currentDegree <= degreeBase; currentDegree++)
    {
      int lenBytes = (1 + *ptrValue1) * (int)sizeof(int);
      (void)memcpy(ptrValue2, ptrValue1, lenBytes);
      ptrValue1 += *ptrValue1;
      ptrValue1++;
      ptrValue2 += nbrLimbs;
    }
  }
  else
  {
    (void)CopyPolynomial(ptrValue2, ptrValue1, degreeBase);
  }
  SetNumberToOne(&poly2[0]); // Initialize power with polynomial 1.
  degreePower = 0;
  for (unsigned int mask = HALF_INT_RANGE_U; mask > 0U; mask >>= 1)
  {
    // Square polynomial.
    MultPolynomial(degreePower, degreePower, poly2, poly2);
    degreePower <<= 1;
    if (modulusIsZero)
    {
      (void)CopyPolynomial(poly2, polyMultTemp, degreePower);
    }
    else
    {
      int lenBytes = (degreePower + 1) * nbrLimbs * (int)sizeof(int);
      (void)memcpy(poly2, polyMultTemp, lenBytes);
    }
    if (((unsigned int)expon & mask) != 0U)
    {
      MultPolynomial(degreeBase, degreePower, poly1, poly2);
      degreePower += degreeBase;
      if (modulusIsZero)
      {
        (void)CopyPolynomial(poly2, polyMultTemp, degreePower);
      }
      else
      {
        int lenBytes = (degreePower + 1) * nbrLimbs * (int)sizeof(int);
        (void)memcpy(poly2, polyMultTemp, lenBytes);
      }
    }
  }
  *ptrArgument1 = degreePower;
  ptrValue1 = ptrArgument1 + 1;
  ptrValue2 = polyMultTemp;
  // Move power back to values stack.
  if (modulusIsZero)
  {
    ptrValue1 = CopyPolynomial(ptrValue1, polyMultTemp, degreePower);
  }
  else
  {
    for (currentDegree = 0; currentDegree <= degreePower; currentDegree++)
    {
      int len = 1 + *ptrValue2;
      int lenBytes = len * (int)sizeof(int);
      (void)memcpy(ptrValue1, ptrValue2, lenBytes);
      ptrValue1 += len;
      ptrValue2 += nbrLimbs;
    }
  }
  diffPtrs = ptrValue1 - &values[0];
  valuesIndex = (int)diffPtrs;
  return EXPR_OK;
}

void computePower(int expo)
{
  int lenBytes;
  (void)BigIntPowerIntExp(&primeMod, expo, &powerMod);
  lenBytes = powerMod.nbrLimbs * (int)sizeof(limb);
  (void)memcpy(TestNbr, powerMod.limbs, lenBytes);
  NumberLength = powerMod.nbrLimbs;
  TestNbr[NumberLength].x = 0;
  GetMontgomeryParms(powerMod.nbrLimbs);
}
// Compute polynomial using value array as stack.
// Inside stack:
// Monomial: -degree, nbrLimbs, xxx (limbs)
// Polynomial: degree, coeff deg 0, coeff deg 1, ...
// where coefficient: nbrlimbs, xxx (limbs)
// Exponent: exponent (integer)
int ComputePolynomial(const char* input, int expo)
{
  int stackIndex = 0;
  bool insideExpon = false;
  int* ptrValue1;
  int* ptrValue2;
  int len;
  char* ptrRPNbuffer;
  enum eExprErr rc;
  int val;
  int pwr;
  int expon;
  int nbrEqualSigns = 0;
  bool usingVariables;
  bool randomUsed;
  degree = 1;
  exponentMod = expo;
  if (modulusIsZero)
  {
    intToBigInteger(&primeMod, 1);
  }
  // Use operand1 as temporary variable to store the exponent.
  computePower(expo);
  rc = ConvertToReversePolishNotation(input, &ptrRPNbuffer, stFuncOperPolyExpr,
    PARSE_EXPR_POLYNOMIAL, &usingVariables, &randomUsed);
  if (rc != EXPR_OK)
  {
    return rc;
  }
  if (!usingVariables && !randomUsed)
  {   // Input string has no variables.
#ifdef __EMSCRIPTEN__
    databack(onlyEvaluate ? "N" : "M");  // Use integer factorization calculator.
#endif
    return EXPR_OK;
  }
  stackIndex = 0;
  valuesIndex = 0;
  while (*ptrRPNbuffer != '\0')
  {
    int nbrSizeBytes;
    int nbrParameters;
    switch (*ptrRPNbuffer)
    {
    case TOKEN_NUMBER:
      stackValues[stackIndex] = &values[valuesIndex];
      stackIndex++;
      ptrRPNbuffer++;
      len = 1;
      if (!insideExpon)
      {
        len = ((int)(unsigned char)*ptrRPNbuffer * 256) + (int)(unsigned char)*(ptrRPNbuffer + 1);
        values[valuesIndex] = 0;   // Degree.
        valuesIndex++;
        values[valuesIndex] = len;
        valuesIndex++;
        ptrRPNbuffer += 2;
      }
      nbrSizeBytes = len * (int)sizeof(limb);
      (void)memcpy(&values[valuesIndex], ptrRPNbuffer, nbrSizeBytes);
      ptrRPNbuffer += nbrSizeBytes - 1;
      valuesIndex += len;
      break;
    case TOKEN_VAR:
      stackValues[stackIndex] = &values[valuesIndex];
      stackIndex++;
      values[valuesIndex] = -1;   // Degree of monomial
      valuesIndex++;
         // Set coefficient to 1 in Montgomery notation.
      SetNumberToOne(&values[valuesIndex]);
      valuesIndex += 1 + values[valuesIndex];
      break;
    case TOKEN_START_EXPON:
      insideExpon = true;
      break;
    case TOKEN_END_EXPON:
      insideExpon = false;
      break;
    case OPER_UNARY_MINUS:
      if (stackIndex < 1)
      {
        return EXPR_CANNOT_PARSE_EXPRESSION;
      }
      ptrValue1 = stackValues[stackIndex - 1];
      if (insideExpon)
      {
        *ptrValue1 = -*ptrValue1;
      }
      else
      {
        rc = NegatePolynomialExpr(ptrValue1);
        if (rc != EXPR_OK)
        {
          return rc;
        }
      }
      break;
    case TOKEN_ANS:
      stackValues[stackIndex] = &values[valuesIndex];
      ptrValue1 = stackValues[stackIndex];
      stackIndex++;
      *ptrValue1 = LastAnswerPoly[0];   // Copy degree.
      if (LastAnswerPoly[1] == 0)
      {
        *(ptrValue1 + 1) = 1;
        *(ptrValue1 + 2) = 0;
        valuesIndex += 3;   // This polynomial requires 3 limbs.
      }
      else
      {
        ptrValue2 = CopyPolynomial(ptrValue1+1, &LastAnswerPoly[1],
                                   LastAnswerPoly[0]);
        valuesIndex = (int)(ptrValue2 - &values[0]);
      }
      break;
    case TOKEN_GCD:
      if (stackIndex < 1)
      {
        return EXPR_CANNOT_PARSE_EXPRESSION;
      }
      ptrRPNbuffer++;
      nbrParameters = (int)(unsigned char)*ptrRPNbuffer;
      ptrValue1 = stackValues[stackIndex - 1];
      for (int parmNbr = 1; parmNbr < nbrParameters; parmNbr++)
      {
        stackIndex--;
        if (stackIndex < 1)
        {
          return EXPR_CANNOT_PARSE_EXPRESSION;
        }
        ptrValue2 = stackValues[stackIndex - 1];
        rc = GcdPolynomialExpr(ptrValue2, ptrValue1);
        if (rc != EXPR_OK)
        {
          return rc;
        }
        ptrValue1 = ptrValue2;
      }
      break;
    case TOKEN_LCM:
      ptrRPNbuffer++;
      nbrParameters = (int)(unsigned char)*ptrRPNbuffer;
      if (stackIndex < 1)
      {
        return EXPR_CANNOT_PARSE_EXPRESSION;
      }
      ptrValue1 = stackValues[stackIndex - 1];
      for (int parmNbr = 1; parmNbr < nbrParameters; parmNbr++)
      {
        stackIndex--;
        if (stackIndex < 1)
        {
          return EXPR_CANNOT_PARSE_EXPRESSION;
        }
        ptrValue2 = stackValues[stackIndex - 1];
        rc = LcmPolynomialExpr(ptrValue1, ptrValue2);
        if (rc != EXPR_OK)
        {
          return rc;
        }
        ptrValue1 = ptrValue2;
      }
      break;
    case TOKEN_DER:
      if (stackIndex < 1)
      {
        return EXPR_CANNOT_PARSE_EXPRESSION;
      }
      ptrValue1 = stackValues[stackIndex - 1];
      DerPolynomial(ptrValue1);
      break;
    case TOKEN_RANDOM:
      if (stackIndex < 4)
      {
        return EXPR_CANNOT_PARSE_EXPRESSION;
      }
      rc = RandomPolynomialExpr(stackValues[stackIndex - 4], stackValues[stackIndex - 3],
        stackValues[stackIndex - 2], stackValues[stackIndex - 1], stackValues[stackIndex - 4]);
      stackIndex -= 3;
      if (rc != EXPR_OK)
      {
        return rc;
      }
      break;
    case OPER_ADD:
      stackIndex--;
      if (stackIndex < 1)
      {
        return EXPR_CANNOT_PARSE_EXPRESSION;
      }
      ptrValue2 = stackValues[stackIndex];
      ptrValue1 = stackValues[stackIndex - 1];
      if (insideExpon)
      {
        *ptrValue1 += *ptrValue2;
        if ((unsigned int)*ptrValue1 >= LIMB_RANGE)
        {
          return EXPR_EXPONENT_TOO_LARGE;
        }
      }
      else
      {
        rc = AddPolynomialExpr(ptrValue1, ptrValue2);
        if (rc != EXPR_OK)
        {
          return rc;
        }
      }
      break;
    case OPER_EQUAL:
      stackIndex--;
      if (stackIndex < 1)
      {
        return EXPR_CANNOT_PARSE_EXPRESSION;
      }
      if (nbrEqualSigns > 0)
      {
        return EXPR_MORE_THAN_ONE_EQUAL_SIGN;
      }
      if ((stackIndex != 1) || insideExpon)
      {
        return EXPR_EQUAL_SIGN_INSIDE_PAREN;
      }
      nbrEqualSigns++;
      ptrValue2 = stackValues[stackIndex];
      ptrValue1 = stackValues[stackIndex - 1];
      rc = NegatePolynomialExpr(ptrValue2);
      if (rc != EXPR_OK)
      {
        return rc;
      }
      rc = AddPolynomialExpr(ptrValue1, ptrValue2);
      if (rc != EXPR_OK)
      {
        return rc;
      }
      break;
    case OPER_SUBT:
      stackIndex--;
      if (stackIndex < 1)
      {
        return EXPR_CANNOT_PARSE_EXPRESSION;
      }
      ptrValue2 = stackValues[stackIndex];
      ptrValue1 = stackValues[stackIndex - 1];
      if (insideExpon)
      {
        *ptrValue1 -= *ptrValue2;
        if ((unsigned int)*ptrValue1 >= LIMB_RANGE)
        {
          return EXPR_EXPONENT_TOO_LARGE;
        }
      }
      else
      {
        rc = NegatePolynomialExpr(ptrValue2);
        if (rc != EXPR_OK)
        {
          return rc;
        }
        rc = AddPolynomialExpr(ptrValue1, ptrValue2);
        if (rc != EXPR_OK)
        {
          return rc;
        }
      }
      break;
    case OPER_MULTIPLY:
      stackIndex--;
      if (stackIndex < 1)
      {
        return EXPR_CANNOT_PARSE_EXPRESSION;
      }
      ptrValue2 = stackValues[stackIndex];
      ptrValue1 = stackValues[stackIndex - 1];
      if (insideExpon)
      {
        *ptrValue1 *= *ptrValue2;
        if ((unsigned int)*ptrValue1 >= LIMB_RANGE)
        {
          return EXPR_EXPONENT_TOO_LARGE;
        }
      }
      else
      {
        rc = MultPolynomialExpr(ptrValue1, ptrValue2);
        if (rc != EXPR_OK)
        {
          return rc;
        }
      }
      break;
    case OPER_DIVIDE:
      stackIndex--;
      if (stackIndex < 1)
      {
        return EXPR_CANNOT_PARSE_EXPRESSION;
      }
      ptrValue2 = stackValues[stackIndex];
      ptrValue1 = stackValues[stackIndex - 1];
      if (insideExpon)
      {
        *ptrValue1 /= *ptrValue2;
      }
      else
      {
        rc = DivPolynomialExpr(ptrValue1, ptrValue2, TYPE_DIVISION);
        if (rc != EXPR_OK)
        {
          return rc;
        }
      }
      break;
    case OPER_REMAINDER:
      stackIndex--;
      if (stackIndex < 1)
      {
        return EXPR_CANNOT_PARSE_EXPRESSION;
      }
      ptrValue2 = stackValues[stackIndex];
      ptrValue1 = stackValues[stackIndex - 1];
      if (insideExpon)
      {
        *ptrValue1 %= *ptrValue2;
      }
      else
      {
        rc = DivPolynomialExpr(ptrValue1, ptrValue2, TYPE_MODULUS);
        if (rc != EXPR_OK)
        {
          return rc;
        }
      }
      break;
    case OPER_POWER:
      stackIndex--;
      if (stackIndex < 1)
      {
        return EXPR_CANNOT_PARSE_EXPRESSION;
      }
      expon = *stackValues[stackIndex];
      if (expon < 0)
      {
        return EXPR_EXPONENT_NEGATIVE;
      }
      if (insideExpon)
      {
        val = *stackValues[stackIndex - 1];
        pwr = 1;
        if (pow(val, pwr) >= MAX_VALUE_LIMB)
        {
          return EXPR_EXPONENT_TOO_LARGE;
        }
        for (int count = expon; count > 0; count--)
        {
          pwr *= val;
        }
        *stackValues[stackIndex - 1] = pwr;
      }
      else
      {
        rc = PowerPolynomialExpr(stackValues[stackIndex - 1], expon);
        if (rc != EXPR_OK)
        {
          return rc;
        }
      }
      break;
    default:
      break;
    }
    ptrRPNbuffer++;
  }
  // Adjust degree so the leading coefficient is not zero.
  if (values[0] < 0)
  {  // Monomial
    UncompressBigIntegerB(&values[1], &operand1);
    degree = -values[0];
    values[0] = degree;
    ptrValue1 = &values[1];
    for (int currentDegree = 0; currentDegree < degree; currentDegree++)
    {
      *ptrValue1 = 1;   // Initialize coefficient to zero.
      ptrValue1++;
      *ptrValue1 = 0;
      ptrValue1++;
    }
    BigInteger2IntArray(ptrValue1, &operand1);
  }
  LastAnswerPoly[0] = values[0];   // Copy degree.
  (void)CopyPolynomial(&LastAnswerPoly[1], &values[1], values[0]);
  return EXPR_OK;
}

static enum eExprErr GcdPolynomialExpr(int* ptrArgument1, int* ptrArgument2)
{
  int degreeGcd;
  const int* ptrNextArgument;
  size_t diffPtrs;
  if (modulusIsZero)
  {    // Integer GCD polynomial.
    int polyDegree;
    PolynomialGcd(ptrArgument1, ptrArgument2, ptrArgument1);
    polyDegree = *ptrArgument1;
    ptrNextArgument = ptrArgument1 + 1;
    for (int currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
    {
      int numLen = 1 + numLimbs(ptrNextArgument);
      ptrNextArgument += numLen;
    }
  }
  else
  {    // Modular GCD polynomial.
    int degree1 = *ptrArgument1;
    int degree2 = *ptrArgument2;
    NumberLength = powerMod.nbrLimbs;
    PolyModularGcd(ptrArgument1 + 1, degree1, ptrArgument2 + 1, degree2,
      poly5, &degreeGcd);
    ptrNextArgument = CopyPolynomial(ptrArgument1 + 1, poly5, degreeGcd);
    *ptrArgument1 = degreeGcd;
  }
  diffPtrs = ptrNextArgument - &values[0];
  valuesIndex = (int)diffPtrs;
  return EXPR_OK;
}

static enum eExprErr LcmPolynomialExpr(int* ptrArgument1, int* ptrArgument2)
{
  enum eExprErr retcode;
  const int* ptrNextArgument;
  size_t diffPtrs;
      // Integer LCM polynomial.
  if ((*ptrArgument1 == 0) && (*(ptrArgument1 + 1) == 1) && (*(ptrArgument1 + 2) == 0))
  {         // argument1 equals zero, so the LCM is also zero.
    valuesIndex = 3;
    return EXPR_OK;
  }
  if ((*ptrArgument2 == 0) && (*(ptrArgument2 + 1) == 1) && (*(ptrArgument2 + 2) == 0))
  {         // argument2 equals zero, so the LCM is also zero.
    *ptrArgument1 = 0;
    *(ptrArgument1 + 1) = 1;
    *(ptrArgument1 + 2) = 0;
    valuesIndex = 3;
    return EXPR_OK;
  }
  // Save polynomials before they are overwritten by PolynomialGcd.
  polyA[0] = *ptrArgument1;
  (void)CopyPolynomial(&polyA[1], ptrArgument1 + 1, *ptrArgument1);
  polyB[0] = *ptrArgument2;
  (void)CopyPolynomial(&polyB[1], ptrArgument2 + 1, *ptrArgument2);
  retcode = GcdPolynomialExpr(ptrArgument1, ptrArgument2);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  retcode = DivPolynomialExpr(polyA, ptrArgument1, TYPE_DIVISION);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  retcode = MultPolynomialExpr(polyA, polyB);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  *ptrArgument2 = polyA[0];
  ptrNextArgument = CopyPolynomial(ptrArgument2 + 1, &polyA[1], polyA[0]);
  diffPtrs = ptrNextArgument - ptrArgument2;
  valuesIndex = (int)diffPtrs;
  return EXPR_OK;
}

static enum eExprErr RandomPolynomialExpr(const int* pMinDegree, const int* pMaxDegree,
  const int* pMinCoeff, const int* pMaxCoeff, int *randomPoly)
{
  int minDegree;
  int maxDegree;
  int polyDegree;
  static BigInteger bigMinCoeff;
  static BigInteger bigMaxCoeff;
  static BigInteger bigCurrentCoeff;
  int* ptrRandomPoly = randomPoly;
  if ((*pMinDegree != 0) || (*pMaxDegree != 0) ||
    (*pMinCoeff != 0) || (*pMaxCoeff != 0))
  {
    return EXPR_ARGUMENT_MUST_BE_CONSTANT;
  }
  if ((*(pMinDegree + 1) != 1) || (*(pMaxDegree + 1) != 1))
  {
    return EXPR_DEGREE_TOO_HIGH;
  }
  minDegree = *(pMinDegree + 2);
  maxDegree = *(pMaxDegree + 2);
  if (minDegree > maxDegree)
  {       // Incorrect order. Exchange them.
    int temp = minDegree;
    minDegree = maxDegree;
    maxDegree = temp;
  }
  if (minDegree > MAX_DEGREE)
  {
    return EXPR_DEGREE_TOO_HIGH;
  }
  polyDegree = intRandom(minDegree, maxDegree);
  if (modulusIsZero)
  {
    UncompressBigIntegerB(pMinCoeff + 1, &bigMinCoeff);
    UncompressBigIntegerB(pMaxCoeff + 1, &bigMaxCoeff);
  }
  else
  {
    UncompressLimbsBigInteger((const limb *)pMinCoeff + 1, &bigMinCoeff);
    UncompressLimbsBigInteger((const limb *)pMaxCoeff + 1, &bigMaxCoeff);
  }
  // Generate polynomial of degree "degree" with coefficients
  // between "bigMinCoeff" and "bigMaxCoeff".
  *ptrRandomPoly = polyDegree;
  ptrRandomPoly++;
  for (int currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {
    BigIntRandom(&bigMinCoeff, &bigMaxCoeff, &bigCurrentCoeff);
    if (modulusIsZero)
    {
      NumberLength = bigCurrentCoeff.nbrLimbs;
      BigInteger2IntArray(ptrRandomPoly, &bigCurrentCoeff);
      ptrRandomPoly += NumberLength;
      ptrRandomPoly++;
    }
  }
  return EXPR_OK;
}
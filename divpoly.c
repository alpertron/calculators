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
#include "bignbr.h"
#include "expression.h"
#include "highlevel.h"
#include "polynomial.h"
#include "showtime.h"

extern int valuesIndex;
extern int NumberLength;
extern int poly4[COMPRESSED_POLY_MAX_LENGTH];
static int revDividend[COMPRESSED_POLY_MAX_LENGTH];
static int inverseDivisor[COMPRESSED_POLY_MAX_LENGTH];
static int polyTmp[COMPRESSED_POLY_MAX_LENGTH];

// Decompress polynomial.
static void ToPoly(int polyDegree, const int* polySrc, int* polyDest)
{
  int lenBytes;
  int currentDegree;
  int* ptrPolyDest = polyDest;
  const int *ptrPolySrc = polySrc + 1;
  if (polyDegree < 0)
  {    // Polynomial is a monomial
    for (currentDegree = 0; currentDegree < -polyDegree; currentDegree++)
    {
      *ptrPolyDest = 1;
      *(ptrPolyDest + 1) = 0;
      ptrPolyDest += NumberLength + 1;
    }
    lenBytes = (*ptrPolySrc+1) * (int)sizeof(int);
    (void)memcpy(ptrPolyDest, ptrPolySrc, lenBytes);
  }
  else
  {   // Polynomial
    for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
    {
      int nbrLimbs = *ptrPolySrc+1;
      lenBytes = nbrLimbs * (int)sizeof(int);
      (void)memcpy(ptrPolyDest, ptrPolySrc, lenBytes);
      ptrPolyDest += NumberLength + 1;
      ptrPolySrc += nbrLimbs;
    }
  }
}

// Compress polynomial.
static void FromPoly(int polyDegree, int* polyDest, const int* polySrc)
{
  const int* ptrPolySrc = polySrc;
  int* ptrPolyDest = polyDest;
  *ptrPolyDest = polyDegree;
  ptrPolyDest++;
  for (int currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {
    int nbrLimbs = *ptrPolySrc+1;
    int lenBytes = nbrLimbs * (int)sizeof(int);
    (void)memcpy(ptrPolyDest, ptrPolySrc, lenBytes);
    ptrPolySrc += NumberLength + 1;
    ptrPolyDest += nbrLimbs;
  }
}

static void ReversePolynomial(int* ptrDest, const int* ptrSrc, const int *ptrDenom)
{
  int* pDest = ptrDest;
  int indexes[(2 * MAX_DEGREE) + 1];
  int* ptrIndex;
  int index;
  int numLength;
  int lenBytes;
  int degreePoly = *ptrSrc;
  if (degreePoly < 0)
  {    // Monomial.
    degreePoly = -degreePoly;
    *pDest = degreePoly;
    pDest++;
    numLength = numLimbs(ptrSrc + 1) + 1;
    lenBytes = numLength * (int)sizeof(int);
    (void)memcpy(pDest, ptrSrc + 1, lenBytes);
    pDest += numLength;
    if (ptrDenom != NULL)
    {             // Append denominator to coefficient.
      numLength = numLimbs(ptrDenom);
      *pDest = *ptrDenom;
      pDest++;
      lenBytes = numLength * (int)sizeof(int);
      (void)memcpy(pDest, ptrDenom + 1, lenBytes);
      pDest += numLength;
    }
    for (degree = 0; degree < degreePoly; degree++)
    {
      *pDest = 1;   // Set coefficient to zero.
      pDest++;
      *pDest = 0;
      pDest++;
      if (ptrDenom != NULL)
      {             // Append denominator to coefficient.
        *pDest = 1;   // Set denominator to one.
        pDest++;
        *pDest = 1;
        pDest++;
      }
    }
    return;
  }
  // Fill indexes to start of each coefficient.
  ptrIndex = &indexes[0];
  index = 1;
  for (degree = 0; degree <= degreePoly; degree++)
  {
    *ptrIndex = index;
    ptrIndex++;
    index += numLimbs(ptrSrc + index) + 1;
  }
  // Copy to destination.
  *pDest = degreePoly;
  pDest++;
  for (degree = degreePoly; degree >= 0; degree--)
  {
    const int* ptrSrcCoeff = ptrSrc + indexes[degree];
    numLength = numLimbs(ptrSrcCoeff) + 1;
    lenBytes = numLength * (int)sizeof(int);
    (void)memcpy(pDest, ptrSrcCoeff, lenBytes);
    pDest += numLength;
    if (ptrDenom != NULL)
    {     // Copy denominator to polynomial with rational coefficients.
      numLength = numLimbs(ptrDenom) + 1;
      lenBytes = numLength * (int)sizeof(int);
      (void)memcpy(pDest, ptrDenom, lenBytes);
      pDest += numLength;
    }
  }
}

int DivideIntegerPolynomial(int* pDividend, const int* pDivisor, enum eDivType type)
{
  int* ptrResult;
  int degreeDividend;
  int degreeDivisor;
  int* ptrQuotient;
  // Move arguments to temporary storage with most significant coefficient
  // first.
  ReversePolynomial(poly1, pDividend, NULL);
  ReversePolynomial(poly2, pDivisor, NULL);
  degreeDividend = poly1[0];
  degreeDivisor = poly2[0];
  if (degreeDividend < degreeDivisor)
  {      // Degree of dividend is less than degree of divisor.
    if (type == TYPE_DIVISION)
    {    // Get pointer to quotient.
      poly3[0] = 0;     // Degree of quotient is zero.
      poly3[1] = 1;     // Coefficient is zero.
      poly3[2] = 0;
    }
    else
    {                   // Remainder is equal to dividend.
      (void)CopyPolynomial(poly1, pDividend, *pDividend);
    }
    return EXPR_OK;
  }
  ptrQuotient = poly3;
  *ptrQuotient = degreeDividend - degreeDivisor;
  ptrQuotient++;
  for (int degreeQuotient = degreeDividend - degreeDivisor;
    degreeQuotient >= 0; degreeQuotient--)
  {
    const int* ptrDividend = &poly1[1];
    const int* ptrDivisor = &poly2[1];
    int* ptrRemainder = &poly4[1];
    UncompressBigIntegerB(ptrDividend, &operand1);
    UncompressBigIntegerB(ptrDivisor, &operand2);
    (void)BigIntRemainder(&operand1, &operand2, &operand3);
    if (!BigIntIsZero(&operand3))
    {
      return EXPR_POLYNOMIAL_DIVISION_NOT_INTEGER;
    }
    (void)BigIntDivide(&operand1, &operand2, &operand3);
    NumberLength = operand3.nbrLimbs;
    BigInteger2IntArray(ptrQuotient, &operand3);
    ptrQuotient += NumberLength;
    ptrQuotient++;
    // Calculate remainder.
    if (degreeDivisor == 0)
    {     // Strip leading coefficient of dividend.
      int numLength = 1 + numLimbs(ptrDividend);
      ptrRemainder = &poly1[1];

      for (degree = degreeDividend; degree > 0; degree--)
      {
        int lenBytes;
        ptrDividend += numLength;
        numLength = 1 + numLimbs(ptrDividend);
        lenBytes = numLength * (int)sizeof(int);
        (void)memcpy(ptrRemainder, ptrDividend, lenBytes);
        ptrRemainder += numLength;
      }
    }
    else
    {     // Degree of divisor is greater than zero.
      for (degree = degreeDivisor; degree > 0; degree--)
      {
        ptrDividend += numLimbs(ptrDividend) + 1;
        ptrDivisor += numLimbs(ptrDivisor) + 1;
        UncompressBigIntegerB(ptrDivisor, &operand2);
        // Multiply by quotient (operand3).
        (void)BigIntMultiply(&operand2, &operand3, &operand2);
        UncompressBigIntegerB(ptrDividend, &operand1);
        BigIntSubt(&operand1, &operand2, &operand2);
        NumberLength = operand2.nbrLimbs;
        BigInteger2IntArray(ptrRemainder, &operand2);
        ptrRemainder += NumberLength;
        ptrRemainder++;
      }
      // Copy least significant coefficients of dividend into remainder.
      for (degree = degreeDividend - degreeDivisor; degree > 0; degree--)
      {
        int lenBytes;
        ptrDividend += numLimbs(ptrDividend);
        ptrDividend++;
        int numLength = 1 + numLimbs(ptrDividend);
        lenBytes = numLength * (int)sizeof(int);
        (void)memcpy(ptrRemainder, ptrDividend, lenBytes);
        ptrRemainder += numLength;
      }
      // Copy remainder to dividend.
      (void)memcpy(&poly1[1], &poly4[1], (char*)ptrRemainder - (char*)&poly4[1]);
    }
    degreeDividend--;
  }
  if (type == TYPE_DIVISION)
  {    // Get pointer to quotient.
    ptrResult = poly3;
    degree = poly1[0] - degreeDivisor;
  }
  else
  {    // Get pointer to remainder.
    ptrResult = poly1;
    degree = degreeDivisor - 1;
  }
  // Compute degree discarding leading coefficients set to zero.
  while (degree > 0)
  {
    if ((*(ptrResult + 1) != 1) || (*(ptrResult + 2) != 0))
    {            // Coefficient is not zero.
      break;
    }
    ptrResult += 2;
    degree--;
  }
  *ptrResult = degree;
  // Copy result to first parameter.
  ReversePolynomial(pDividend, ptrResult, NULL);
  // Discard most significant 
  return EXPR_OK;
}

static enum eExprErr divideNumAndDenByGcd(BigInteger* num, BigInteger* den,
  BigInteger* tmp)
{
  BigIntGcd(num, den, tmp);
  enum eExprErr err = BigIntDivide(num, tmp, num);
  if (err != EXPR_OK)
  {
    return err;
  }
  return BigIntDivide(den, tmp, den);
}

// Divide (ptrDividend/ptrDividendDen) / (ptrDivisor/ptrDivisorDen)
// Place the result in operand1/operand2.
static enum eExprErr RationalDivide(const int *ptrDividend, const int *ptrDividendDen,
  const int *ptrDivisor, const int *ptrDivisorDen)
{
  enum eExprErr err;
  // Get numerator of division.
  UncompressBigIntegerB(ptrDividend, &operand1);
  UncompressBigIntegerB(ptrDivisorDen, &operand2);
  err = BigIntMultiply(&operand1, &operand2, &operand1);
  if (err != EXPR_OK)
  {
    return err;
  }
  // Get denominator of division.
  UncompressBigIntegerB(ptrDivisor, &operand2);
  UncompressBigIntegerB(ptrDividendDen, &operand3);
  err = BigIntMultiply(&operand2, &operand3, &operand2);
  if (err != EXPR_OK)
  {
    return err;
  }
  return divideNumAndDenByGcd(&operand1, &operand2, &operand3);
}

// Compute minNum/minDen - subt1Num/subt1Den * subt2Num/subt2Den
// Place the result in operand5/operand4.
static enum eExprErr RationalSubtractProduct(const int *minNum, const int *minDen,
  const int *subt1Num, const int *subt1Den, const int* subt2Num, const int* subt2Den)
{
  enum eExprErr err;
  // Multiply numerators.
  UncompressBigIntegerB(subt1Num, &operand1);
  UncompressBigIntegerB(subt2Num, &operand2);
  err = BigIntMultiply(&operand1, &operand2, &operand3);
  if (err != EXPR_OK)
  {
    return err;
  }
  // Multiply denominators.
  UncompressBigIntegerB(subt1Den, &operand1);
  UncompressBigIntegerB(subt2Den, &operand2);
  err = BigIntMultiply(&operand1, &operand2, &operand4);
  if (err != EXPR_OK)
  {
    return err;
  }
  // Divide them by their GCD.
  err = divideNumAndDenByGcd(&operand3, &operand4, &operand1);
  if (err != EXPR_OK)
  {
    return err;
  }
  // At this moment subtrahend is operand3/operand4.
  // Compute numerator of subtraction.
  UncompressBigIntegerB(minNum, &operand1);
  err = BigIntMultiply(&operand1, &operand4, &operand2);
  if (err != EXPR_OK)
  {
    return err;
  }
  UncompressBigIntegerB(minDen, &operand1);
  err = BigIntMultiply(&operand1, &operand3, &operand5);
  if (err != EXPR_OK)
  {
    return err;
  }
  BigIntSubt(&operand2, &operand5, &operand5);
  // Compute denominator of subtraction.
  err = BigIntMultiply(&operand1, &operand4, &operand4);
  if (err != EXPR_OK)
  {
    return err;
  }
  // Divide them by their GCD.
  return divideNumAndDenByGcd(&operand5, &operand4, &operand3);
}

// polyDest must be different from polySrc.
static enum eExprErr ConvertPolynomialRatCoeffToRatPoly(const int* polySrc, int* polyDest,
  bool reverse, int polyDegree)
{
  int currentDegree;
  const int* indexes[(2 * MAX_DEGREE) + 1];
  const int* ptrSrc = polySrc;
  int* ptrDest;
  // Set operand1 to the LCM of the dividends of all coefficients.
  intToBigInteger(&operand1, 1);
  for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {  // Point to denominator.
    if (reverse)
    {
      indexes[polyDegree - currentDegree] = ptrSrc;
    }
    else
    {
      indexes[currentDegree] = ptrSrc;
    }
    // Point to denominator.
    ptrSrc += numLimbs(ptrSrc) + 1;
    UncompressBigIntegerB(ptrSrc, &operand2);
    // Compute LCM of denominators.
    enum eExprErr err = BigIntLcm(&operand1, &operand2, &operand1);
    if (err != EXPR_OK)
    {
      return err;
    }
    // point to next numerator.
    ptrSrc += numLimbs(ptrSrc) + 1;
  }
  // Generate result.
  *polyDest = degree;
  ptrDest = polyDest + 1;
  for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {  // Point to denominator.
    ptrSrc = indexes[currentDegree];
    UncompressBigIntegerB(ptrSrc, &operand2);
    enum eExprErr err = BigIntMultiply(&operand1, &operand2,
      &operand2);
    if (err != EXPR_OK)
    {
      return err;
    }
    ptrSrc += numLimbs(ptrSrc) + 1;
    UncompressBigIntegerB(ptrSrc, &operand3);
    err = BigIntDivide(&operand2, &operand3, &operand2);
    if (err != EXPR_OK)
    {
      return err;
    }
    NumberLength = operand2.nbrLimbs;
    BigInteger2IntArray(ptrDest, &operand2);
    ptrDest += NumberLength;
    ptrDest++;
  }
  // Copy denominator.
  *ptrDest = 0;   // Initialize degree of denominator
  ptrDest++;
  BigInteger2IntArray(ptrDest, &operand1);
  return EXPR_OK;
}

int DivideRationalPolynomial(int* pDividend, const int* pDivisor, enum eDivType type)
{
  enum eExprErr err;
  const int* ptrResult;
  int degreeDividend;
  int degreeDivisor;
  int* ptrQuotient;
  int* denomin1 = getNextElement(pDividend);
  int* denomin2 = getNextElement(pDivisor);
  if ((*denomin1 != 0) || (*denomin2 != 0))
  {   // Degree of denominator is not zero.
    return EXPR_DENOMINATOR_MUST_BE_CONSTANT;
  }
  if (type == TYPE_MODULUS)
  {
    if (*pDivisor == 0)
    {     // Degree of divisor is zero.
      if ((*(pDivisor + 1) == 1) && (*(pDivisor + 2) == 0))
      {   // Divisor is zero. Return dividend.
        return EXPR_OK;
      }
      else
      {   // Divisor is constant. Return zero.
        *pDividend = 0;       // Degree of numerator.
        *(pDividend + 1) = 1; // Numerator is zero.
        *(pDividend + 2) = 0;
        *(pDividend + 3) = 0; // Degree of denominator.
        *(pDividend + 4) = 1; // Denominator is one.
        *(pDividend + 5) = 1;
        return EXPR_OK;
      }
    }
  }
  // Move arguments to temporary storage with most significant coefficient
  // first. Append the divisor for each coefficient.
  ReversePolynomial(poly1, pDividend, denomin1 + 1);
  ReversePolynomial(poly2, pDivisor, denomin2 + 1);
  degreeDividend = poly1[0];
  degreeDivisor = poly2[0];
  if (degreeDividend < degreeDivisor)
  {      // Degree of dividend is less than degree of divisor.
    if (type == TYPE_DIVISION)
    {    // Get pointer to quotient.
      poly3[0] = 0;     // Degree of quotient is zero.
      poly3[1] = 1;     // Coefficient is zero.
      poly3[2] = 0;
    }
    else
    {                   // Remainder is equal to dividend.
      (void)CopyPolynomial(poly1, pDividend, *pDividend);
    }
    return EXPR_OK;
  }
  ptrQuotient = poly3;
  *ptrQuotient = degreeDividend - degreeDivisor;
  ptrQuotient++;
  for (int degreeQuotient = degreeDividend - degreeDivisor;
    degreeQuotient >= 0; degreeQuotient--)
  {
    const int* ptrDividend = &poly1[1];
    const int* ptrDivisor = &poly2[1];
    int* ptrRemainder;
    const int* ptrDividendDen = ptrDividend + numLimbs(ptrDividend) + 1;
    const int* ptrDivisorDen = ptrDivisor + numLimbs(ptrDivisor) + 1;
    err = RationalDivide(ptrDividend, ptrDividendDen,
      ptrDivisor, ptrDivisorDen);
    if (err != EXPR_OK)
    {
      return err;
    }
    // Save numerator of quotient (operand1).
    NumberLength = operand1.nbrLimbs;
    BigInteger2IntArray(ptrQuotient, &operand1);
    int* ptrQuotientDen = ptrQuotient + numLimbs(ptrQuotient) + 1;
    // Save denominator of quotient (operand2).
    NumberLength = operand2.nbrLimbs;
    BigInteger2IntArray(ptrQuotientDen, &operand2);
    // Calculate remainder.
    if (degreeDivisor == 0)
    {     // Strip leading coefficient of dividend.
      int numLength = 1 + numLimbs(ptrDividend);
      ptrRemainder = &poly1[1];

      for (degree = degreeDividend*2; degree > 0; degree--)
      {   // Degree is multiplied by 2 to move both numerator and denominator.
        int lenBytes;
        ptrDividend += numLength;
        numLength = 1 + numLimbs(ptrDividend);
        lenBytes = numLength * (int)sizeof(int);
        (void)memcpy(ptrRemainder, ptrDividend, lenBytes);
        ptrRemainder += numLength;
      }
    }
    else
    {     // Degree of divisor is greater than zero.
      ptrRemainder = &poly4[1];
          // Skip most significant coefficients of dividend and divisor.
      ptrDividendDen = ptrDividend + numLimbs(ptrDividend) + 1;
      ptrDivisorDen = ptrDivisor + numLimbs(ptrDivisor) + 1;
      ptrDividend = ptrDividendDen + numLimbs(ptrDividendDen) + 1;
      ptrDivisor = ptrDivisorDen + numLimbs(ptrDivisorDen) + 1;
      for (degree = degreeDivisor; degree > 0; degree--)
      {
        ptrDividendDen = ptrDividend + numLimbs(ptrDividend) + 1;
        ptrDivisorDen = ptrDivisor + numLimbs(ptrDivisor) + 1;
        err = RationalSubtractProduct(ptrDividend, ptrDividendDen,
          ptrDivisor, ptrDivisorDen, ptrQuotient, ptrQuotientDen);
        if (err != EXPR_OK)
        {
          return err;
        }
        // Remainder is located in operand5/operand4.
        NumberLength = operand5.nbrLimbs;
        BigInteger2IntArray(ptrRemainder, &operand5);
        ptrRemainder += NumberLength;
        ptrRemainder++;
        NumberLength = operand4.nbrLimbs;
        BigInteger2IntArray(ptrRemainder, &operand4);
        ptrRemainder += NumberLength;
        ptrRemainder++;
        ptrDividend = ptrDividendDen + numLimbs(ptrDividendDen) + 1;
        ptrDivisor = ptrDivisorDen + numLimbs(ptrDivisorDen) + 1;
      }
      // Copy least significant coefficients of dividend into remainder.
      for (degree = (degreeDividend - degreeDivisor)*2; degree > 0; degree--)
      {
        int numLength = 1 + numLimbs(ptrDividend);
        int lenBytes = numLength * (int)sizeof(int);
        (void)memcpy(ptrRemainder, ptrDividend, lenBytes);
        ptrRemainder += numLength;
        ptrDividend += numLength;
      }
      // Copy remainder to dividend.
      (void)memcpy(&poly1[1], &poly4[1], (char*)ptrRemainder - (char*)&poly4[1]);
    }
    ptrQuotient = ptrQuotientDen + numLimbs(ptrQuotientDen) + 1;
    degreeDividend--;
  }
  if (type == TYPE_DIVISION)
  {    // Get pointer to quotient.
    degree = poly1[0] - degreeDivisor;
    return ConvertPolynomialRatCoeffToRatPoly(&poly3[1], pDividend,
      true, degree);
  }
  // Get pointer to remainder.
  ptrResult = &poly1[1];
  degree = degreeDivisor - 1;
  // Compute degree discarding leading coefficients set to zero.
  while (degree > 0)
  {
    if ((*ptrResult != 1) || (*(ptrResult + 1) != 0))
    {            // Coefficient is not zero.
      break;
    }
    ptrResult += 2;
    ptrResult += numLimbs(ptrResult) + 1; // Skip divisor.
    degree--;
  }
  return ConvertPolynomialRatCoeffToRatPoly(ptrResult, pDividend,
    false, degree);
}

// ptrArgument1 is the dividend and ptrArgument2 is the divisor.
// If type equals TYPE_DIVISION, ptrArgument1 is overwritten with the quotient.
// If type equals TYPE_MODULUS, ptrArgument1 is overwritten with the remainder.
enum eExprErr DivPolynomialExpr(int* ptrArgument1, const int* ptrArgument2,
  enum eDivType type)
{
  int currentDegree;
  int degree1 = *ptrArgument1;
  int degree2 = *ptrArgument2;
  if ((*ptrArgument2 == 0) && (*(ptrArgument2 + 1) == 1) && (*(ptrArgument2 + 2) == 0))
  {        // Divisor is zero
    if (type == TYPE_DIVISION)
    {
      return EXPR_DIVIDE_BY_ZERO;
    }
    return EXPR_OK;   // a mod 0 = a.
  }
  if ((degree1 <= 0) && (degree2 <= 0))
  {        // Division of two monomials.
    size_t diffPtrs;
    if (degree1 > degree2)
    {      // Degree of dividend less than degree of divisor.
      if (type == TYPE_DIVISION)
      {       // Result is zero.
        *ptrArgument1 = 0;
        *(ptrArgument1 + 1) = 1;
        *(ptrArgument1 + 2) = 0;
      }
      return EXPR_OK;
    }
    if (type == TYPE_MODULUS)
    {       // Result is zero.
      *ptrArgument1 = 0;
      *(ptrArgument1 + 1) = 1;
      *(ptrArgument1 + 2) = 0;
      return EXPR_OK;
    }
    *ptrArgument1 = degree1 - degree2;
    UncompressBigIntegerB(ptrArgument1 + 1, &operand1);
    UncompressBigIntegerB(ptrArgument2 + 1, &operand2);
    if (modulusIsZero)
    {
      (void)BigIntRemainder(&operand1, &operand2, &operand3);
      if (!BigIntIsZero(&operand3))
      {    // Remainder is not zero.
        return EXPR_POLYNOMIAL_DIVISION_NOT_INTEGER;
      }
      (void)BigIntDivide(&operand1, &operand2, &operand3);
      CopyBigInt(&operand1, &operand3);
      NumberLength = operand1.nbrLimbs;
    }
    else
    {
      (void)ModInvBigNbr(operand2.limbs, operand2.limbs, TestNbr, NumberLength);
      modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    }
    BigInteger2IntArray(ptrArgument1 + 1, &operand1);
    diffPtrs = ptrArgument1 - &values[0];
    valuesIndex = (int)diffPtrs + 2 + *(ptrArgument1 + 1);
    return EXPR_OK;
  }
  if (modulusIsZero)
  {
    return DivideIntegerPolynomial(ptrArgument1, ptrArgument2, type);
  }
  ToPoly(degree1, ptrArgument1, poly1); // Move dividend to poly1.
  ToPoly(degree2, ptrArgument2, poly2); // Move divisor to poly2.
  if (degree1 < 0)
  {
    degree1 = -degree1;
  }
  if (degree2 < 0)
  {
    degree2 = -degree2;
  }
  if (degree1 < degree2)
  {       // Degree of dividend less than degree of divisor.
    if (type == TYPE_DIVISION)
    {     // Result is zero.
      *ptrArgument1 = 0;
      *(ptrArgument1 + 1) = 1;
      *(ptrArgument1 + 2) = 0;
    }
    return EXPR_OK;
  }
  DividePolynomial(poly1, degree1, poly2, degree2, poly3); // Set poly3 to quotient.
  if (type == TYPE_DIVISION)
  {
    currentDegree = degree1 - degree2;
    *ptrArgument1 = currentDegree;
    FromPoly(currentDegree, ptrArgument1, poly3); // Move dividend to poly1.
  }
  else
  {
    currentDegree = getDegreePoly(poly1, degree2 - 1);
    *ptrArgument1 = currentDegree;
    FromPoly(currentDegree, ptrArgument1, poly1); // Move modulus to poly1.
  }
  return EXPR_OK;
}

// Reverse coefficients of polynomials. Both polynomials must be different.
static void ReverseModularPolynomial(const int* ptrSrc, int* ptrRev, int polyDegree)
{
  const int* pSrc = ptrSrc;
  int* pRev = ptrRev;
  int nbrLimbs = NumberLength + 1;
  pRev += polyDegree * nbrLimbs;
  for (int currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {
    int lenBytes = nbrLimbs * (int)sizeof(limb);
    (void)memcpy(pRev, pSrc, lenBytes);
    pSrc += nbrLimbs;
    pRev -= nbrLimbs;
  }
}

// Perform modular division of polynomials pDividend/pDivisor (pDivisor is monic).
static void PolynomialNewtonDivision(/*@in@*/int* pDividend, int dividendDegree,
  const int* pDivisor, int divisorDegree, /*@out@*/int* ptrQuotient)
{
  int quotientDegree = dividendDegree - divisorDegree;
  int oldDegree = 1;
  int currentDegree;
  int* ptrRemainder;
  int *ptrProd;
  int nbrLimbs = NumberLength + 1;
  int degrees[15];
  int nbrDegrees = 0;
  int newtonDegree;

  // Use revDividend as temporary storage for f.
  ReverseModularPolynomial(pDivisor, revDividend, divisorDegree);
  if (divisorDegree < quotientDegree)
  {
    ptrProd = &revDividend[(divisorDegree + 1) * nbrLimbs];
    for (currentDegree = divisorDegree; currentDegree < quotientDegree; currentDegree++)
    {
      *ptrProd = 1;       // Initialize coefficient to zero.
      *(ptrProd + 1) = 0;
      ptrProd += nbrLimbs;
    }
  }

  // Inversion of polynomial f:
  // 1) g <- 1, pow2 <- 2
  // 2) while pow2 < deg f:
  // 2.1)  g <- g * (2 - f*g) (mod x^pow2)
  // 2.2)  pow2 <- 2*pow2
  // 3) g <- g * (2 - f*g) (mod x^deg f)

  // Set polynomial g to 1.
  if (NumberLength == 1)
  {
    inverseDivisor[0] = 1;
    inverseDivisor[1] = 1;
  }
  else
  {
    ArrLimbs2LenAndLimbs(inverseDivisor, MontgomeryMultR1, nbrLimbs);
  }
  newtonDegree = quotientDegree + 1;
  // Compute degrees to use in Newton loop.
  while (newtonDegree > 1)
  {
    degrees[nbrDegrees] = newtonDegree;
    nbrDegrees++;
    newtonDegree = (newtonDegree + 1) / 2;
  }
  // Perform Newton loop.
  while (--nbrDegrees >= 0)
  {
    // Compute g <- g * (2 - f*g) (mod x^degree)
    // f is stored in revDividend.
    // g is stored in inverseDivisor.
    // oldDegree = size of g.
    const int* ptrProduct;
    int *ptrDest;
    int currDegree;
    int lenBytes;
    int newDegree = degrees[nbrDegrees];
    // Compute f*g.
    MultPolynomial(newDegree - 1, oldDegree - 1, revDividend, inverseDivisor);
    // Set operand1.limbs to 2 in Montgomery notation.
    ptrProduct = polyMultTemp; // Point to start of f*g.
    ptrDest = polyTmp;         // Point to start of 2 - f*g.
    if (NumberLength == 1)
    {
      int mod = TestNbr[0].x;
      // Subtract 2 minus the trailing coefficient of f*g.
      polyTmp[0] = 1;
      polyTmp[1] = 2 - polyMultTemp[1];
      if (polyTmp[1] < 0)
      {
        polyTmp[1] += mod;
      }
      for (currDegree = 1; currDegree < newDegree; currDegree++)
      {                    // Get the negative of all coefficients of f*g.
        ptrProduct += 2;   // Point to next coefficient of f*g.
        ptrDest += 2;      // Point to next coefficient of 2 - f*g.
        *ptrDest = 1;
        if (*(ptrProduct + 1) == 0)
        {
          *(ptrDest+1) = 0;
        }
        else
        {
          *(ptrDest+1) = mod - *(ptrProduct + 1);
        }
      }
    }
    else
    {
      AddBigIntModN(MontgomeryMultR1, MontgomeryMultR1,
        operand1.limbs, TestNbr, NumberLength);
      // Subtract 2 minus the trailing coefficient of f*g.
      LenAndLimbs2ArrLimbs(polyMultTemp, operand2.limbs, nbrLimbs);
      SubtBigNbrMod(operand1.limbs, operand2.limbs, operand2.limbs);
      ArrLimbs2LenAndLimbs(polyTmp, operand2.limbs, nbrLimbs);
      lenBytes = nbrLimbs * (int)sizeof(limb);
      (void)memset(operand1.limbs, 0, lenBytes);
      for (currDegree = 1; currDegree < newDegree; currDegree++)
      {                    // Get the negative of all coefficients of f*g.
        ptrProduct += nbrLimbs;  // Point to next coefficient of f*g.
        ptrDest += nbrLimbs;     // Point to next coefficient of 2 - f*g.
        LenAndLimbs2ArrLimbs(ptrProduct, operand2.limbs, nbrLimbs);
        SubtBigNbrMod(operand1.limbs, operand2.limbs, operand2.limbs);
        ArrLimbs2LenAndLimbs(ptrDest, operand2.limbs, nbrLimbs);
      }
    }
    // Compute g * (2 - f*g).
    MultPolynomial(oldDegree - 1, newDegree - 1, inverseDivisor, polyTmp);
    // Store g * (2 - f*g) into g.
    lenBytes = newDegree * nbrLimbs * (int)sizeof(limb);
    (void)memcpy(inverseDivisor, polyMultTemp, lenBytes);
    oldDegree = newDegree;
  }
  ReverseModularPolynomial(pDividend, revDividend, dividendDegree);

  // reverse quotient <- revDividend * inverse mod x^(difDegrees+1)
  MultPolynomial(quotientDegree, quotientDegree, revDividend, inverseDivisor);
  // Reverse coefficients of quotient.
  // Do not store it in ptrQuotient because it may be NULL.
  ReverseModularPolynomial(polyMultTemp, polyTmp, quotientDegree);
  if (ptrQuotient != NULL)
  {     // Copy polyTmp to ptrQuotient.
    int lenBytes = (quotientDegree + 1) * nbrLimbs * (int)sizeof(limb);
    (void)memcpy(ptrQuotient, polyTmp, lenBytes);
  }
  // Compute remainder: pDividend <- pDividend - ptrQuotient * pDivisor.
  // Degree of remainder is at most one less than degree of divisor,
  // so trim quotient polynomial up to degree of divisor.
  if (quotientDegree > divisorDegree)
  {
    quotientDegree = divisorDegree;
  }
  MultPolynomial(quotientDegree, divisorDegree, polyTmp, pDivisor);
  ptrRemainder = pDividend;
  ptrProd = polyMultTemp;
  // Degree of remainder is at most one less than degree of divisor.
  for (currentDegree = 0; currentDegree <= divisorDegree; currentDegree++)
  {
    LenAndLimbs2ArrLimbs(ptrRemainder, operand1.limbs, nbrLimbs);
    LenAndLimbs2ArrLimbs(ptrProd, operand2.limbs, nbrLimbs);
    SubtBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
    ArrLimbs2LenAndLimbs(ptrRemainder, operand1.limbs, nbrLimbs);
    ptrRemainder += nbrLimbs;
    ptrProd += nbrLimbs;
  }
}

// In this routine, the dividend is replaced by the remainder of the division.
// Input and output coefficients are expressed in Montgomery notation.
// If only the remainder is needed, ptrQuotient can be NULL.
void DividePolynomial(/*@in@*/int* pDividend, int dividendDegree,
  /*@in@*/int* pDivisor, int divisorDegree, /*@out@*/int* ptrQuotient)
{
  int currentDegree;
  int nbrLimbs = NumberLength + 1;
  bool divisorIsOne;
  int* ptrDivisor;
  int* ptrDividend;
  int lenBytes;
  if (divisorDegree > dividendDegree)
  {    // Quotient is zero.
    if (ptrQuotient != NULL)
    {
      *ptrQuotient = 1;
      *(ptrQuotient + 1) = 0;
    }
    return;
  }
  IntArray2BigInteger(pDivisor + (divisorDegree * nbrLimbs), &operand1);
  lenBytes = NumberLength * (int)sizeof(int);
  (void)memcpy(operand5.limbs, operand1.limbs, lenBytes);
  if (NumberLength == 1)
  {
    divisorIsOne = operand1.limbs[0].x == 1;
  }
  else
  {
    divisorIsOne = !memcmp(operand1.limbs, MontgomeryMultR1, lenBytes);
  }
  if (!divisorIsOne)
  {        // Leading coefficient is not 1.
    ConvertToMonic(pDivisor, divisorDegree);
    // operand1 holds the inverse of the leading coefficient of divisor.
    // Multiply dividend by this number.
    for (currentDegree = 0; currentDegree <= dividendDegree; currentDegree++)
    {
      IntArray2BigInteger(pDividend + (currentDegree * nbrLimbs), &operand2);
      modmult(operand1.limbs, operand2.limbs, operand2.limbs);
      BigInteger2IntArray(pDividend + (currentDegree * nbrLimbs), &operand2);
    }
  }
  if ((divisorDegree > 16) && (dividendDegree < (4*divisorDegree)))
  {         // Newton division is faster.
    PolynomialNewtonDivision(pDividend, dividendDegree,
      pDivisor, divisorDegree, ptrQuotient);
  }
  else
  {
    int* ptrQuot = NULL;
    if (ptrQuotient != NULL)
    {
      ptrQuot = ptrQuotient + (dividendDegree - divisorDegree) * nbrLimbs;
    }
    for (currentDegree = dividendDegree; currentDegree >= divisorDegree; currentDegree--)
    {
      int index;
      ptrDividend = pDividend + (currentDegree * nbrLimbs);
      IntArray2BigInteger(ptrDividend, &operand1);
      if (ptrQuotient != NULL)
      {
        BigInteger2IntArray(ptrQuot, &operand1);  // Store coefficient of quotient.
      }
      ptrDivisor = pDivisor + (divisorDegree * nbrLimbs);
      if ((NumberLength == 1) && (TestNbr[0].x <= 32768))
      {
        int mod = TestNbr[0].x;
        ptrDividend++;
        ptrDivisor++;
        for (index = 0; index <= divisorDegree; index++)
        {
          *ptrDividend = (*ptrDividend - (*ptrDivisor * operand1.limbs[0].x)) % mod;
          if (*ptrDividend < 0)
          {
            *ptrDividend += mod;
          }
          ptrDividend -= nbrLimbs;
          ptrDivisor -= nbrLimbs;
        }
      }
      else
      {
        for (index = 0; index <= divisorDegree; index++)
        {
          IntArray2BigInteger(ptrDivisor, &operand2);
          modmult(operand1.limbs, operand2.limbs, operand2.limbs);
          IntArray2BigInteger(ptrDividend, &operand3);
          SubtBigNbrMod(operand3.limbs, operand2.limbs, operand3.limbs);
          BigInteger2IntArray(ptrDividend, &operand3);
          ptrDividend -= nbrLimbs;
          ptrDivisor -= nbrLimbs;
        }
      }
      ptrQuot -= nbrLimbs;
    }
  }
  if (!divisorIsOne)
  {        // Leading coefficient is not 1.
           // Adjust remainder by multiplying each coefficient by leading
           // coefficient of divisor.
    ptrDividend = pDividend;
    for (currentDegree = 0; currentDegree <= dividendDegree; currentDegree++)
    {
      IntArray2BigInteger(ptrDividend, &operand2);
      modmult(operand5.limbs, operand2.limbs, operand2.limbs);
      BigInteger2IntArray(ptrDividend, &operand2);
      ptrDividend += nbrLimbs;
    }
    // Restore divisor.
    ptrDivisor = pDivisor;
    for (currentDegree = 0; currentDegree <= divisorDegree; currentDegree++)
    {
      IntArray2BigInteger(ptrDivisor, &operand2);
      modmult(operand5.limbs, operand2.limbs, operand2.limbs);
      BigInteger2IntArray(ptrDivisor, &operand2);
      ptrDivisor += nbrLimbs;
    }
  }
}

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

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "bignbr.h"
#include "highlevel.h"
#include "polynomial.h"
#include "showtime.h"

extern int valuesIndex;
extern int NumberLength;
extern int poly4[COMPRESSED_POLY_MAX_LENGTH];
static int revDividend[COMPRESSED_POLY_MAX_LENGTH];
static int inverseDivisor[COMPRESSED_POLY_MAX_LENGTH];
static int polyTmp[COMPRESSED_POLY_MAX_LENGTH];

static void ToPoly(int polyDegree, int* polySrc, int* polyDest)
{
  int currentDegree;
  polySrc++;
  if (polyDegree < 0)
  {    // Polynomial is a monomial
    for (currentDegree = 0; currentDegree < -polyDegree; currentDegree++)
    {
      *polyDest = 1;
      *(polyDest + 1) = 0;
      polyDest += NumberLength + 1;
    }
    memcpy(polyDest, polySrc, (*(polySrc)+1) * sizeof(int));
  }
  else
  {   // Polynomial
    for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
    {
      int nbrLimbs = *(polySrc)+1;
      memcpy(polyDest, polySrc, nbrLimbs * sizeof(int));
      polyDest += NumberLength + 1;
      polySrc += nbrLimbs;
    }
  }
}

static void FromPoly(int polyDegree, int* polyDest, int* polySrc)
{
  int currentDegree;
  *polyDest++ = polyDegree;
  for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {
    int nbrLimbs = *(polySrc)+1;
    memcpy(polyDest, polySrc, nbrLimbs * sizeof(int));
    polyDest += NumberLength + 1;
    polySrc += nbrLimbs;
  }
}

static void ReversePolynomial(int* ptrDest, int* ptrSrc)
{
  int indexes[2 * MAX_DEGREE + 1];
  int* ptrIndex;
  int index, numLength;
  int degreePoly = *ptrSrc;
  if (degreePoly < 0)
  {    // Monomial.
    degreePoly = -degreePoly;
    *ptrDest++ = degreePoly;
    numLength = numLimbs(ptrSrc + 1) + 1;
    memcpy(ptrDest, ptrSrc + 1, numLength * sizeof(int));
    ptrDest += numLength;
    for (degree = 0; degree < degreePoly; degree++)
    {
      *ptrDest++ = 1;   // Set coefficient to zero.
      *ptrDest++ = 0;
    }
    return;
  }
  // Fill indexes to start of each coefficient.
  ptrIndex = &indexes[0];
  index = 1;
  for (degree = 0; degree <= degreePoly; degree++)
  {
    *ptrIndex++ = index;
    index += numLimbs(ptrSrc + index) + 1;
  }
  // Copy to destination.
  *ptrDest++ = degreePoly;
  for (degree = degreePoly; degree >= 0; degree--)
  {
    int* ptrSrcCoeff = ptrSrc + indexes[degree];
    numLength = numLimbs(ptrSrcCoeff) + 1;
    memcpy(ptrDest, ptrSrcCoeff, numLength * sizeof(int));
    ptrDest += numLength;
  }
}

int DivideIntegerPolynomial(int* pDividend, int* pDivisor, enum eDivType type)
{
  int* ptrResult;
  int degreeDividend;
  int degreeDivisor;
  int degreeQuotient;
  int* ptrQuotient;
  // Move arguments to temporary storage with most significant coefficient
  // first.
  ReversePolynomial(poly1, pDividend);
  ReversePolynomial(poly2, pDivisor);
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
      CopyPolynomial(poly1, pDividend, *pDividend);
    }
    return EXPR_OK;
  }
  ptrQuotient = poly3;
  *ptrQuotient++ = degreeDividend - degreeDivisor;
  for (degreeQuotient = degreeDividend - degreeDivisor;
    degreeQuotient >= 0; degreeQuotient--)
  {
    int* ptrDividend = &poly1[1];
    int* ptrDivisor = &poly2[1];
    int* ptrRemainder = &poly4[1];
    UncompressBigIntegerB(ptrDividend, &operand1);
    UncompressBigIntegerB(ptrDivisor, &operand2);
    BigIntRemainder(&operand1, &operand2, &operand3);
    if (!BigIntIsZero(&operand3))
    {
      return EXPR_POLYNOMIAL_DIVISION_NOT_INTEGER;
    }
    BigIntDivide(&operand1, &operand2, &operand3);
    NumberLength = operand3.nbrLimbs;
    BigInteger2IntArray(ptrQuotient, &operand3);
    ptrQuotient += 1 + NumberLength;
    // Calculate remainder.
    if (degreeDivisor == 0)
    {     // Strip leading coefficient of dividend.
      int numLength = 1 + numLimbs(ptrDividend);
      ptrRemainder = &poly1[1];

      for (degree = degreeDividend; degree > 0; degree--)
      {
        ptrDividend += numLength;
        numLength = 1 + numLimbs(ptrDividend);
        memcpy(ptrRemainder, ptrDividend, numLength * sizeof(int));
        ptrRemainder += numLength;
      }
    }
    else
    {
      for (degree = degreeDivisor; degree > 0; degree--)
      {
        ptrDividend += numLimbs(ptrDividend) + 1;
        ptrDivisor += numLimbs(ptrDivisor) + 1;
        UncompressBigIntegerB(ptrDivisor, &operand2);
        BigIntMultiply(&operand2, &operand3, &operand2);
        UncompressBigIntegerB(ptrDividend, &operand1);
        BigIntSubt(&operand1, &operand2, &operand2);
        NumberLength = operand2.nbrLimbs;
        BigInteger2IntArray(ptrRemainder, &operand2);
        ptrRemainder += 1 + NumberLength;
      }
      // Copy least significant coefficients of dividend into remainder.
      for (degree = degreeDividend - degreeDivisor; degree > 0; degree--)
      {
        ptrDividend += 1 + numLimbs(ptrDividend);
        int numLength = 1 + numLimbs(ptrDividend);
        memcpy(ptrRemainder, ptrDividend, numLength * sizeof(int));
        ptrRemainder += numLength;
      }
      // Copy remainder to dividend.
      memcpy(&poly1[1], &poly4[1], (char*)ptrRemainder - (char*)&poly4[1]);
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
    if (*(ptrResult + 1) != 1 || *(ptrResult + 2) != 0)
    {            // Coefficient is not zero.
      break;
    }
    ptrResult += 2;
    degree--;
  }
  *ptrResult = degree;
  // Copy result to first parameter.
  ReversePolynomial(pDividend, ptrResult);
  // Discard most significant 
  return EXPR_OK;
}

int DivPolynomialExpr(int* ptrArgument1, int* ptrArgument2, enum eDivType type)
{
  int currentDegree;
  int degree1 = *ptrArgument1;
  int degree2 = *ptrArgument2;
  if (*ptrArgument2 == 0 && *(ptrArgument2 + 1) == 1 && *(ptrArgument2 + 2) == 0)
  {        // Divisor is zero
    if (type == TYPE_DIVISION)
    {
      return EXPR_DIVIDE_BY_ZERO;
    }
    return EXPR_OK;   // a mod 0 = a.
  }
  if (degree1 <= 0 && degree2 <= 0)
  {        // Division of two monomials.
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
      BigIntRemainder(&operand1, &operand2, &operand3);
      if (!BigIntIsZero(&operand3))
      {    // Remainder is not zero.
        return EXPR_POLYNOMIAL_DIVISION_NOT_INTEGER;
      }
      BigIntDivide(&operand1, &operand2, &operand3);
      CopyBigInt(&operand1, &operand3);
      NumberLength = operand1.nbrLimbs;
    }
    else
    {
      ModInvBigNbr(operand2.limbs, operand2.limbs, TestNbr, NumberLength);
      modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    }
    BigInteger2IntArray(ptrArgument1 + 1, &operand1);
    valuesIndex = (int)(ptrArgument1 + 2 + *(ptrArgument1 + 1) - &values[0]);
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
static void ReverseModularPolynomial(int *ptrSrc, int *ptrRev, int degree)
{
  int currentDegree;
  int nbrLimbs = NumberLength + 1;
  ptrRev += degree * nbrLimbs;
  for (currentDegree = 0; currentDegree <= degree; currentDegree++)
  {
    memcpy(ptrRev, ptrSrc, nbrLimbs * sizeof(limb));
    ptrSrc += nbrLimbs;
    ptrRev -= nbrLimbs;
  }
}

// Perform modular division of polynomials pDividend/pDivisor (pDivisor is monic).
static void PolynomialNewtonDivision(/*@in@*/int* pDividend, int dividendDegree,
  /*@in@*/int* pDivisor, int divisorDegree, /*@out@*/int* ptrQuotient)
{
  int quotientDegree = dividendDegree - divisorDegree;
  int oldDegree = 1;
  int currentDegree;
  int* ptrRemainder, * ptrProd;
  int nbrLimbs = NumberLength + 1;
  int degrees[15];
  int nbrDegrees = 0;
  int newtonDegree;

  // Use revDividend as temporary storage for f.
  ReverseModularPolynomial(pDivisor, revDividend, divisorDegree);
  if (divisorDegree < quotientDegree)
  {
    ptrProd = &revDividend[(divisorDegree+1) * nbrLimbs];
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
  ArrLimbs2LenAndLimbs(inverseDivisor, MontgomeryMultR1, nbrLimbs);
  newtonDegree = quotientDegree + 1;
  // Compute degrees to use in Newton loop.
  while (newtonDegree > 1)
  {
    degrees[nbrDegrees++] = newtonDegree;
    newtonDegree = (newtonDegree + 1) / 2;
  }
  // Perform Newton loop.
  while (--nbrDegrees >= 0)
  {
    // Compute g <- g * (2 - f*g) (mod x^degree)
    // f is stored in revDividend.
    // g is stored in inverseDivisor.
    // oldDegree = size of g.
    int* ptrProduct, * ptrDest;
    int currentDegree;
    int newDegree = degrees[nbrDegrees];
    // Compute f*g.
    MultPolynomial(newDegree - 1, oldDegree - 1, revDividend, inverseDivisor);
    // Set operand1.limbs to 2 in Montgomery notation.
    AddBigIntModN((int*)MontgomeryMultR1, (int*)MontgomeryMultR1,
      (int*)operand1.limbs, (int*)TestNbr, NumberLength);
    // Subtract 2 minus the trailing coefficient of f*g.
    LenAndLimbs2ArrLimbs(polyMultTemp, operand2.limbs, nbrLimbs);
    SubtBigNbrMod(operand1.limbs, operand2.limbs, operand2.limbs);
    ArrLimbs2LenAndLimbs(polyTmp, operand2.limbs, nbrLimbs);
    ptrProduct = polyMultTemp; // Point to start of f*g.
    ptrDest = polyTmp;         // Point to start of 2 - f*g.
    memset(operand1.limbs, 0, nbrLimbs * sizeof(limb));
    for (currentDegree = 1; currentDegree < newDegree; currentDegree++)
    {                    // Get the negative of all coefficients of f*g.
      ptrProduct += nbrLimbs;  // Point to next coefficient of f*g.
      ptrDest += nbrLimbs;     // Point to next coefficient of 2 - f*g.
      LenAndLimbs2ArrLimbs(ptrProduct, operand2.limbs, nbrLimbs);
      SubtBigNbrMod(operand1.limbs, operand2.limbs, operand2.limbs);
      ArrLimbs2LenAndLimbs(ptrDest, operand2.limbs, nbrLimbs);
    }
    // Compute g * (2 - f*g).
    MultPolynomial(oldDegree - 1, newDegree - 1, inverseDivisor, polyTmp);
    // Store g * (2 - f*g) into g.
    memcpy(inverseDivisor, polyMultTemp, newDegree * nbrLimbs * sizeof(limb));
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
    memcpy(ptrQuotient, polyTmp, (quotientDegree + 1) * nbrLimbs * sizeof(limb));
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
  int currentDegree, index;
  int nbrLimbs = NumberLength + 1;
  int divisorIsOne;
  int* ptrQuot;
  int remainderDegree;
  int* ptrDivisor;
  int* ptrDividend;
  if (divisorDegree > dividendDegree)
  {    // Quotient is zero.
    if (ptrQuotient != NULL)
    {
      *ptrQuotient = 1;
      *(ptrQuotient + 1) = 0;
    }
    return;
  }
  remainderDegree = dividendDegree - divisorDegree;
  IntArray2BigInteger(pDivisor + divisorDegree * nbrLimbs, &operand1);
  memcpy(operand5.limbs, operand1.limbs, NumberLength * sizeof(int));
  divisorIsOne = !memcmp(operand1.limbs, MontgomeryMultR1, NumberLength * sizeof(int));
  if (!divisorIsOne)
  {        // Leading coefficient is not 1.
    ConvertToMonic(pDivisor, divisorDegree);
    // operand1 holds the inverse of the leading coefficient of divisor.
    // Multiply dividend by this number.
    for (currentDegree = 0; currentDegree <= dividendDegree; currentDegree++)
    {
      IntArray2BigInteger(pDividend + currentDegree * nbrLimbs, &operand2);
      modmult(operand1.limbs, operand2.limbs, operand2.limbs);
      BigInteger2IntArray(pDividend + currentDegree * nbrLimbs, &operand2);
    }
  }
  if (divisorDegree > 16 && dividendDegree<4*divisorDegree)
//  if (0)
  {         // Newton division is faster.
    PolynomialNewtonDivision(pDividend, dividendDegree,
      pDivisor, divisorDegree, ptrQuotient);
  }
  else
  {
    ptrQuot = ptrQuotient + (dividendDegree - divisorDegree) * nbrLimbs;
    for (currentDegree = dividendDegree; currentDegree >= divisorDegree; currentDegree--)
    {
      ptrDividend = pDividend + currentDegree * nbrLimbs;
      IntArray2BigInteger(ptrDividend, &operand1);
      if (ptrQuotient != NULL)
      {
        BigInteger2IntArray(ptrQuot, &operand1);  // Store coefficient of quotient.
      }
      ptrDivisor = pDivisor + divisorDegree * nbrLimbs;
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

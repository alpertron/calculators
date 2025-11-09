//
// This file is part of Alpertron Calculators.
//
// Copyright 2019-2024 Dario Alejandro Alpern
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

#include <string.h>
#include <math.h>
#include <assert.h>
#include "string/strings.h"
#include "rootseq.h"
#include "expression.h"
#include "copyStr.h"

extern BigInteger Quadratic;
extern BigInteger Linear;
extern BigInteger Independent;
extern char* ptrOutput;
extern const char* ptrSin;
extern const char* ptrCos;
extern const char* ptrPi;
extern const char* ptrI;

void stepsForQuadraticEquation(char origVar, char substVar)
{
  char currVar = origVar;
  if (!BigIntIsZero(&RatLinear.numerator))
  {
    // To eliminate the linear term, we will perform the following substitution:
    formatString(&ptrOutput, "<p>$1s</p><p>", LITERAL_STEPS_QUADR1);
    showVariable(&ptrOutput, origVar);
    showText(" = ");
    showVariable(&ptrOutput, substVar);
    BigRationalDivideByInt(&RatLinear, -2, &Rat2);
    showPlusMinusRational(&Rat2);
    showText("</p><p>");
    // The constant value in the substitution equals half of the linear coefficient.
    showText(LITERAL_STEPS_QUADR2);
    currVar = substVar;
    // Show a(y-k)^2 + b(y-k) + c = 0.
    showText("</p><p>");
    startParen();
    showVariable(&ptrOutput, substVar);
    showPlusMinusRational(&Rat2);
    endParen();
    showPower(&ptrOutput, 2);
    showCoeffBeforeParen(&RatLinear);
    showText(ptrTimes);
    startParen();
    showVariable(&ptrOutput, substVar);
    showPlusMinusRational(&Rat2);
    endParen();
    showPlusMinusRational(&RatIndependent);
    showText(" = 0</p><p>");
    // Expanding brackets:
    showText(LITERAL_STEPS_QUADR3);
    showText("</p><p>");
    // Expand all terms.
    // (y-B/2)^2
    showRatCoeffAndPowerVar(NULL, -2, substVar);
    BigIntChSign(&RatLinear.numerator);
    showRatCoeffAndPowerVar(&RatLinear, 1, substVar);
    BigIntChSign(&RatLinear.numerator);
    BigRationalMultiply(&Rat2, &Rat2, &Rat4);
    showPlusMinusRational(&Rat4);
    // B(y-B/2)
    showRatCoeffAndPowerVar(&RatLinear, 1, substVar);
    BigRationalMultiplyByInt(&Rat4, -2, &Rat4);
    showPlusMinusRational(&Rat4);
    // C
    showPlusMinusRational(&RatIndependent);
    showText(" = 0</p><p>");
     // Simplifying:
    showText(LITERAL_STEPS_QUADR4);
    showText("</p><p>");
    // Show y^2 - delta = 0.
    BigRationalDivideByInt(&Rat4, 2, &Rat4);
    BigRationalAdd(&Rat4, &RatIndependent, &Rat4);
    showRatCoeffAndPowerVar(NULL, -2, substVar);
    showPlusMinusRational(&Rat4);
    showText(" = 0</p>");
  }
  else
  {         // Copy constant term.
    CopyBigInt(&Rat4.numerator, &RatIndependent.numerator);
    CopyBigInt(&Rat4.denominator, &RatIndependent.denominator);
  }
  BigIntChSign(&Rat4.numerator);
  showText("<p>");
  showVariable(&ptrOutput, currVar);
  if (pretty == TEX)
  {
    showText(" = \\pm");
  }
  else
  {
    showText(" = &pm;");
  }
  showSquareRootOfRational(&Rat4, 2, ptrTimes);
  showText("</p>");
  if (currVar != origVar)
  {
    showText("<p>");
    showVariable(&ptrOutput, origVar);
    BigIntChSign(&Rat2.numerator);
    showPlusMinusRational(&Rat2);
    if (pretty == TEX)
    {
      showText(" = \\pm");
    }
    else
    {
      showText(" = &pm;");
    }
    showSquareRootOfRational(&Rat4, 2, ptrTimes);
    showText("</p>");
  }
}

// From variables named Quadratic, Linear and Independent, compute the
// discriminant and find the rational roots if there are. Otherwise find
// the real or imaginary roots.
// The calculator factors the polynomial before finding the factors, so
// the roots cannot be rational now because in this case it would
// found the linear factors.
void ProcessQuadraticEquation(enum eSign* pSignDescr)
{
  // Compute discriminant (delta = linear^2 - 4*quadratic*independent).
  (void)BigIntMultiply(&Linear, &Linear, &tmp1);
  (void)BigIntMultiply(&Quadratic, &Independent, &tmp2);
  multint(&tmp2, &tmp2, 4);
  BigIntSubt(&tmp1, &tmp2, &discr);
  *pSignDescr = discr.sign;
  discr.sign = SIGN_POSITIVE;
  CopyBigInt(&RatLinear.numerator, &Linear);
  CopyBigInt(&RatLinear.denominator, &Quadratic);
  CopyBigInt(&RatIndependent.numerator, &Independent);
  CopyBigInt(&RatIndependent.denominator, &Quadratic);
  // Reduce numerators and denominators.
  BigRationalMultiplyByInt(&RatLinear, 1, &RatLinear);
  BigRationalMultiplyByInt(&RatIndependent, 1, &RatIndependent);
  if (teach)
  {
    if (!BigIntIsOne(&Quadratic))
    {
      // Dividing the equation by the quadratic coefficient:
      formatString(&ptrOutput, "<p>$1s</p></p>", LITERAL_QUADRATIC_EQ1);
      showRatCoeffAndPowerVar(NULL, -2, 'x');
      showRatCoeffAndPowerVar(&RatLinear, 1, 'x');
      showRatCoeffAndPowerVar(&RatIndependent, 0, 'x');
      showText(" = 0</p>");
    }
    stepsForQuadraticEquation('x', 'y');
  }
  // Compute Rat1 as -linear/(2*quadratic), Rat2 as abs(1/(2*quadratic))
  // and Rat3 as abs(delta).
  intToBigInteger(&Rat2.numerator, 1);
  CopyBigInt(&Rat2.denominator, &Quadratic);
  CopyBigInt(&Rat3.numerator, &discr);
  intToBigInteger(&Rat3.denominator, 1);
  BigRationalDivideByInt(&RatLinear, -2, &Rat1);
  BigRationalDivideByInt(&Rat2, 2, &Rat2);
  ForceDenominatorPositive(&Rat1);
  ForceDenominatorPositive(&Rat2);
  Rat2.numerator.sign = SIGN_POSITIVE;
}

// Compute delta = c_1^2 - 4c_0 * c_2
// If delta > 0 -> x = (-c_1 +/- sqrt(delta))/(2*c_2)
// If delta < 0 -> x = (-c_1 +/- i*sqrt(-delta))/(2*c_2)
// Delta cannot be a perfect square because the polynomial is already factored.
void QuadraticEquation(const int* polynomial, int multiplicity)
{
  const int* ptrPolynomial = polynomial;
  enum eSign signDiscr;
  UncompressBigIntegerB(ptrPolynomial, &Independent);
  ptrPolynomial += numLimbs(ptrPolynomial);
  ptrPolynomial++;
  UncompressBigIntegerB(ptrPolynomial, &Linear);
  ptrPolynomial += numLimbs(ptrPolynomial);
  ptrPolynomial++;
  UncompressBigIntegerB(ptrPolynomial, &Quadratic);
  ProcessQuadraticEquation(&signDiscr);
  // Show irrational / complex roots.
  for (int ctr = 0; ctr < 2; ctr++)
  {
    showX(multiplicity);
    if (!BigIntIsZero(&Linear))
    {
      showRationalNoParen(&Rat1);
    }
    *ptrOutput = ' ';
    ptrOutput++;
    if (ctr == 1)
    {
      if (!BigIntIsZero(&Linear))
      {
        *ptrOutput = '+';
        ptrOutput++;
      }
    }
    else
    {
      showText(ptrMinus);
    }
    *ptrOutput = ' ';
    ptrOutput++;
    MultiplyRationalBySqrtRational(&Rat2, &Rat3);
    ShowRationalAndSqrParts(&Rat2, &Rat3, 2, ptrTimes);
    if (signDiscr == SIGN_NEGATIVE)
    {
      *ptrOutput = ' ';
      ptrOutput++;
      showText(ptrTimes);
      showText(ptrI);
    }
    endShowX();
  }
}

// If polynomial has the form ax^(2n) + bx^n + c, show its roots.
bool isQuadraticExponential(const int* ptrPolynomial, int polyDegree, int multiplicity)
{
  const int* ptrPoly = ptrPolynomial;
  char rootQuadr[30000];
  char* ptrOutputBak;
  enum eSign signDiscr;
  int ctr;
  int currentDegree;
  int halfDegree = polyDegree / 2;
  char degreeStr[20];
  char* ptrDegreeStr;
  if ((polyDegree % 2) != 0)
  {          // If degree is odd, it is not quadratic, so go out.
    return false;
  }
  // Get constant term.
  NumberLength = *ptrPoly;
  IntArray2BigInteger(ptrPoly, &Independent);
  ptrPoly += numLimbs(ptrPoly);
  ptrPoly++;
  for (currentDegree = 1; currentDegree < halfDegree; currentDegree++)
  {
    if ((*ptrPoly != 1) || (*(ptrPoly + 1) != 0))
    {
      return false;    // Polynomial does not have format ax^n + b.
    }
    ptrPoly += numLimbs(ptrPoly);
    ptrPoly++;
  }
  // Get linear coefficient.
  NumberLength = *ptrPoly;
  IntArray2BigInteger(ptrPoly, &Linear);
  ptrPoly += numLimbs(ptrPoly);
  ptrPoly++;
  for (currentDegree = 1; currentDegree < halfDegree; currentDegree++)
  {
    if ((*ptrPoly != 1) || (*(ptrPoly + 1) != 0))
    {
      return false;    // Polynomial does not have format ax^n + b.
    }
    ptrPoly += numLimbs(ptrPoly);
    ptrPoly++;
  }
  // Get leading coefficient.
  NumberLength = *ptrPoly;
  IntArray2BigInteger(ptrPoly, &Quadratic);
  // Discriminant cannot be a perfect square.
  ProcessQuadraticEquation(&signDiscr);
  if (signDiscr == SIGN_POSITIVE)
  {           // Roots of quadratic equation are real.
    for (ctr = 0; ctr < 2; ctr++)
    {         // Show (Rat1 +/- Rat2*sqrt(Rat3))
      enum eSign Rat2SignBak;
      bool isNegative = false;
      if (Rat1.numerator.sign == Rat2.numerator.sign)
      {
        if (Rat1.numerator.sign == SIGN_NEGATIVE)
        {
          isNegative = true;
        }
      }
      else
      {   // Both terms have different signs.
        BigRationalMultiply(&Rat2, &Rat2, &Rat4);
        BigRationalMultiply(&Rat4, &Rat3, &Rat5);
        BigRationalMultiply(&Rat1, &Rat1, &Rat4);
        BigRationalSubt(&Rat5, &Rat4, &Rat4);
        if (Rat1.numerator.sign == Rat4.numerator.sign)
        {
          isNegative = true;
        }
      }
      if (isNegative)
      {     // Change sign of Rat1.
        Rat1.numerator.sign ^= SIGN_POSITIVE ^ SIGN_NEGATIVE;
        Rat2.numerator.sign ^= SIGN_POSITIVE ^ SIGN_NEGATIVE;
      }
      ptrOutputBak = ptrOutput;
      ptrOutput = rootQuadr;
      if (!BigIntIsZero(&Linear))
      {
        showRationalNoParen(&Rat1);
      }
      *ptrOutput = ' ';
      ptrOutput++;
      Rat2SignBak = Rat2.numerator.sign;
      if (Rat2.numerator.sign == SIGN_POSITIVE)
      {
        if (!BigIntIsZero(&Linear))
        {
          *ptrOutput = '+';
          ptrOutput++;
        }
      }
      else
      {
        Rat2.numerator.sign = SIGN_POSITIVE;
        showText(ptrMinus);
      }
      *ptrOutput = ' ';
      ptrOutput++;
      MultiplyRationalBySqrtRational(&Rat2, &Rat3);
      ShowRationalAndSqrParts(&Rat2, &Rat3, 2, ptrTimes);
      *ptrOutput = 0;   // Mark end of string
      Rat2.numerator.sign = Rat2SignBak;
      ptrOutput = ptrOutputBak;
      GenerateRoots(multiplicity, rootQuadr, isNegative, halfDegree);
      if (isNegative)
      {     // Change sign of Rat1.
        Rat1.numerator.sign ^= SIGN_POSITIVE ^ SIGN_NEGATIVE;
        Rat2.numerator.sign ^= SIGN_POSITIVE ^ SIGN_NEGATIVE;
      }
      // Change sign of second term.
      Rat2.numerator.sign ^= SIGN_POSITIVE ^ SIGN_NEGATIVE;
    }
    return true;
  }
  // Roots: (Rat1^2+Rat2^2*Rat3)^(1/deg)*cos/sin((1/(deg/2))*(n*pi+/-arctan(Rat2*sqrt(Rat3)/Rat1))
  BigRationalMultiply(&Rat2, &Rat2, &Rat4);
  BigRationalMultiply(&Rat4, &Rat3, &Rat5);
  BigRationalMultiply(&Rat1, &Rat1, &Rat4);
  BigRationalAdd(&Rat5, &Rat4, &Rat4);  // Square of absolute value.
  BigRationalDivide(&Rat2, &Rat1, &Rat5);
  Rat5.numerator.sign = SIGN_POSITIVE;
  Rat5.denominator.sign = SIGN_POSITIVE;
  MultiplyRationalBySqrtRational(&Rat2, &Rat3); // Argument of arc tangent.
  ptrDegreeStr = degreeStr;
  int2dec(&ptrDegreeStr, halfDegree);
  *ptrDegreeStr = 0;
  for (currentDegree = 0; currentDegree < halfDegree; currentDegree++)
  {
    int multiplicand;
    for (ctr = 0; ctr < 2; ctr++)
    {
      for (int component = 0; component < 2; component++)
      {
        if (component == 0)
        {            // Showing real part.
          showX(multiplicity);
        }
        else
        {            // Showing imaginary part.
          copyStr(&ptrOutput, (pretty == PARI_GP) ? " + I " : " + i ");
          copyStr(&ptrOutput, ptrTimes);
          *ptrOutput = ' ';
          ptrOutput++;
        }
        StartRadicand(polyDegree);
        showRationalNoParen(&Rat4);
        EndRadicand(polyDegree);
        *ptrOutput = ' ';
        ptrOutput++;
        copyStr(&ptrOutput, ptrTimes);
        if (component != 0)
        {
          copyStr(&ptrOutput, ptrSin);
        }
        else
        {
          copyStr(&ptrOutput, ptrCos);
        }
        startParen();
        if (pretty != PARI_GP)
        {
          showRatString("1", degreeStr);
          *ptrOutput = ' ';
          ptrOutput++;
          copyStr(&ptrOutput, ptrTimes);
        }
        else
        {
          *ptrOutput = '1';
          ptrOutput++;
          *ptrOutput = '/';
          ptrOutput++;
          copyStr(&ptrOutput, degreeStr);
          *ptrOutput = '*';
          ptrOutput++;
        }
        startParen();
        multiplicand = 2 * currentDegree;
        if (Rat1.numerator.sign == SIGN_NEGATIVE)
        {
          multiplicand++;
        }
        if (multiplicand != 0)
        {
          if (multiplicand != 1)
          {
            int2dec(&ptrOutput, multiplicand);
            if (pretty != TEX)
            {
              copyStr(&ptrOutput, (pretty != PARI_GP) ? "&#8290; " : "*");
            }
          }
          copyStr(&ptrOutput, ptrPi);
          *ptrOutput = ' ';
          ptrOutput++;
        }
        if (ctr == 1)
        {
          if (multiplicand != 0)
          {
            copyStr(&ptrOutput, "+ ");
          }
        }
        else
        {
          copyStr(&ptrOutput, (pretty == PRETTY_PRINT) ? "&minus; " : "- ");
        }
        if (pretty == PRETTY_PRINT)
        {
          copyStr(&ptrOutput, "arctan");
        }
        else if (pretty == TEX)
        {
          copyStr(&ptrOutput, "\\arctan{");
        }
        else
        {
          copyStr(&ptrOutput, "atan(");
        }
        ShowRationalAndSqrParts(&Rat2, &Rat3, 2, ptrTimes);
        if (pretty == TEX)
        {
          *ptrOutput = '}';  // Close arc tangent.
          ptrOutput++;
        }
        endParen();
        endParen();
        if (pretty == PARI_GP)
        {
          *ptrOutput = ')';  // Close sine or cosine.
          ptrOutput++;
        }
        if (pretty == TEX)
        {
          *ptrOutput = '}';  // Close sine or cosine.
          ptrOutput++;
        }
      }
      endShowX();
    }
  }
  return true;
}


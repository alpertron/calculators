//
// This file is part of Alpertron Calculators.
//
// Copyright 2019-2023 Dario Alejandro Alpern
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
#include "rootseq.h"
#include "expression.h"

extern const char* ptrI;
extern const char* ptrCos;
extern char* ptrOutput;

static void showFirstTermQuarticEq(int ctr)
{
  if (!BigIntIsZero(&RatCubic.numerator))
  {
    showRationalNoParen(&RatCubic);
    showPlusSignOn((ctr == 0) || (ctr == 1), TYPE_PM_SPACE_BEFORE);
  }
  else
  {
    if ((ctr == 2) || (ctr == 3))
    {
      showText(ptrMinus);
    }
  }
}

static void biquadraticEquation(int multiplicity)
{
  int ctr;
  // For biquadratic equation, the depressed independent term equals the independent term (r = e).

  // Solutions are x = +/- sqrt((sqrt(r)-p/2)/2) +/- sqrt((-sqrt(r)-p/2)/2)
  // Test whether e is a perfect square.
  if (BigRationalSquareRoot(&RatDeprIndependent, &Rat3))
  {           // e is a perfect square. Rat3 = sqrt(r).
    BigRationalDivideByInt(&RatDeprQuadratic, 2, &RatDeprQuadratic);
    for (ctr = 0; ctr < 4; ctr++)
    {
      bool isSquareRoot1;
      bool isSquareRoot2;
      showX(multiplicity);
      BigRationalSubt(&Rat3, &RatDeprQuadratic, &Rat1);
      BigRationalAdd(&Rat3, &RatDeprQuadratic, &Rat2);
      BigIntChSign(&Rat2.numerator);
      BigRationalDivideByInt(&Rat1, 2, &Rat1);
      BigRationalDivideByInt(&Rat2, 2, &Rat2);
      ForceDenominatorPositive(&Rat1);
      ForceDenominatorPositive(&Rat2);
      // Rat1 must be greater than Rat2. Exchange them if needed.
      BigRationalSubt(&Rat1, &Rat2, &Rat4);
      if (Rat4.numerator.sign != Rat4.denominator.sign)
      {            // Rat1 < Rat2 -> exchange them.
        CopyBigInt(&tmp1, &Rat1.numerator);
        CopyBigInt(&Rat1.numerator, &Rat2.numerator);
        CopyBigInt(&Rat2.numerator, &tmp1);
        CopyBigInt(&tmp1, &Rat1.denominator);
        CopyBigInt(&Rat1.denominator, &Rat2.denominator);
        CopyBigInt(&Rat2.denominator, &tmp1);
      }
      isSquareRoot1 = BigRationalSquareRoot(&Rat1, &Rat4);
      isSquareRoot2 = BigRationalSquareRoot(&Rat2, &Rat5);
      // They cannot be both perfect squares because in that case the polynomial is reducible.
      showFirstTermQuarticEq(ctr);
      if (Rat1.numerator.sign == SIGN_NEGATIVE)
      {
        showText(ptrI);
        showText(" ");
        showText(ptrTimes);
        showText(" (");
      }
      if (isSquareRoot1)
      {
        showRational(&Rat4);
      }
      else
      {
        Rat1.numerator.sign = SIGN_POSITIVE;
        showSquareRootOfRational(&Rat1, 2, ptrTimes);
      }
      showPlusSignOn((ctr == 0) || (ctr == 2), TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
      if (isSquareRoot2)
      {
        showRational(&Rat5);
      }
      else
      {
        enum eSign sign = Rat2.numerator.sign;  // Back up sign.
        Rat2.numerator.sign = SIGN_POSITIVE;
        showSquareRootOfRational(&Rat2, 2, ptrTimes);
        Rat2.numerator.sign = sign;             // Restore sign.
      }
      if (Rat2.numerator.sign == SIGN_NEGATIVE)
      {
        if (Rat1.numerator.sign == SIGN_NEGATIVE)
        {
          showText(")");
        }
        else
        {
          showText(ptrTimes);
          showText(ptrI);
        }
      }
      if (Rat1.numerator.sign == SIGN_NEGATIVE)
      {
        showText(")");
      }
      endLine();
    }
  }
  else
  {      // e is not a perfect square.
         // Compute discriminant as c^2-4e
    BigRationalMultiply(&RatDeprQuadratic, &RatDeprQuadratic, &RatDiscr);
    BigRationalMultiplyByInt(&RatDeprIndependent, 4, &Rat1);
    BigRationalSubt(&RatDiscr, &Rat1, &RatDiscr);
    ForceDenominatorPositive(&RatDiscr);
    intToBigInteger(&Rat1.numerator, 1);
    intToBigInteger(&Rat1.denominator, 2);
    CopyBigInt(&Rat2.numerator, &RatDeprIndependent.numerator);
    CopyBigInt(&Rat2.denominator, &RatDeprIndependent.denominator);
    MultiplyRationalBySqrtRational(&Rat1, &Rat2);
    if (RatDiscr.numerator.sign == SIGN_NEGATIVE)
    {
      BigRationalDivideByInt(&RatDeprQuadratic, 4, &RatDeprQuadratic);
    }
    else
    {
      BigRationalDivideByInt(&RatDeprQuadratic, 2, &RatDeprQuadratic);
    }
    for (ctr = 0; ctr < 4; ctr++)
    {
      enum eSign sign = RatDeprQuadratic.numerator.sign;
      RatDeprQuadratic.numerator.sign = SIGN_POSITIVE;
      showX(multiplicity);
      showFirstTermQuarticEq(ctr);
      if (RatDiscr.numerator.sign == SIGN_NEGATIVE)
      {  // x = +/- sqrt((sqrt(e) - c) / 2) +/- i * sqrt((sqrt(e) + c) / 2)
        for (int ctr2 = 0; ctr2 < 2; ctr2++)
        {
          if (ctr2 == 1)
          {
            showPlusSignOn((ctr == 0) || (ctr == 2), TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
            showText(ptrI);
            showText(ptrTimes);
            *ptrOutput = ' ';
            ptrOutput++;
          }
          if (BigIntIsZero(&RatDeprQuadratic.numerator))
          {
            ShowRationalAndSqrParts(&Rat1, &Rat2, 4, ptrTimes);
          }
          else
          {
            startSqrt();
            ShowRationalAndSqrParts(&Rat1, &Rat2, 2, ptrTimes);
            if (sign == SIGN_POSITIVE)
            {
              showPlusSignOn(ctr2 != 0, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
            }
            else
            {
              showPlusSignOn(ctr2 == 0, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
            }
            showRational(&RatDeprQuadratic);
            endSqrt();
          }
        }
      }
      else
      {  // x = sqrt(-c/2 +/- sqrt(c^2 - 4e)/2)
        bool isX2Positive;
        // If c >= 0 and positive sqrt => c^2 < c^2 - 4e => e < 0.
        // If c > 0 and negative sqrt = > x^2 is always negative.
        // If c < 0 and positive sqrt = > x^2 is always positive.
        // If c <= 0 and negative sqrt = > c^2 > c^2 - 4e => e > 0
        if ((ctr == 0) || (ctr == 2))
        {    // Positive sqrt
          isX2Positive = ((RatDeprQuadratic.numerator.sign == SIGN_NEGATIVE) ||
            (RatDeprIndependent.numerator.sign == SIGN_NEGATIVE));
        }
        else
        {    // Negative sqrt
          isX2Positive = (((RatDeprQuadratic.numerator.sign == SIGN_NEGATIVE) ||
            BigIntIsZero(&RatDeprQuadratic.numerator)) &&
            (RatDeprIndependent.numerator.sign == SIGN_POSITIVE));
        }
        if (!isX2Positive)
        {
          showText(ptrI);
          showText(" ");
          showText(ptrTimes);
        }
        ForceDenominatorPositive(&RatDeprIndependent);
        if (BigIntIsZero(&RatDeprQuadratic.numerator))
        {
          BigIntChSign(&RatDeprIndependent.numerator);
          showSquareRootOfRational(&RatDeprIndependent, 4, ptrTimes);
          BigIntChSign(&RatDeprIndependent.numerator);
          endLine();
        }
        else
        {
          startSqrt();
          CopyBigInt(&Rat1.numerator, &RatDeprQuadratic.numerator);
          CopyBigInt(&Rat1.denominator, &RatDeprQuadratic.denominator);
          if (isX2Positive)
          {
            BigIntChSign(&Rat1.numerator);
          }
          showRational(&Rat1);
          if ((ctr == 0) || (ctr == 2))
          {    // Positive sqrt
            showPlusSignOn(isX2Positive, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
          }
          else
          {    // Negative sqrt
            showPlusSignOn(!isX2Positive, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
          }
          BigRationalMultiply(&RatDeprQuadratic, &RatDeprQuadratic, &Rat1);
          BigRationalSubt(&Rat1, &RatDeprIndependent, &Rat1);
          showSquareRootOfRational(&Rat1, 2, ptrTimes);
          endSqrt();
          endLine();
        }
      }
    }
  }
}

// Show real or imaginary part of square root of u + q/S
static void showSquareRootOfComplex(const char* plus, const char* minus)
{
  startSqrt();
  // Divide square root part by 4*16.
  BigRationalDivideByInt(&Rat4, 4 * 16, &Rat4);
  // Divide rational part by 2*4.
  BigRationalDivideByInt(&Rat3, 2 * 4, &Rat3);
  showSquareRootOfRational(&Rat4, 2, ptrTimes);
  if (Rat3.numerator.sign == SIGN_POSITIVE)
  {
    showText(plus);
    showRationalNoParen(&Rat3);
  }
  else
  {
    showText(minus);
    BigIntChSign(&Rat3.numerator);
    showRationalNoParen(&Rat3);
    BigIntChSign(&Rat3.numerator);
  }
  endSqrt();
  // Restore square root part.
  BigRationalMultiplyByInt(&Rat4, 4 * 16, &Rat4);
  // Restore rational part.
  BigRationalMultiplyByInt(&Rat3, 2 * 4, &Rat3);
}

static void FerrariResolventHasRationalRoot(int multiplicity)
{
  const int* ptrValues = factorInfoInteger[0].ptrPolyLifted;
  UncompressBigIntegerB(ptrValues, &RatS.numerator);
  ptrValues += numLimbs(ptrValues);
  ptrValues++;
  UncompressBigIntegerB(ptrValues, &RatS.denominator);   // RatS <- -root
  BigRationalDivideByInt(&RatS, -2, &RatS);              // RatS <- S^2 (as on Wikipedia article).
  ForceDenominatorPositive(&RatS);
  BigRationalMultiplyByInt(&RatS, 4, &Rat3);             // 4S^2
  BigRationalAdd(&Rat3, &RatDeprQuadratic, &Rat3);
  BigRationalAdd(&Rat3, &RatDeprQuadratic, &Rat3);       // 4S^2 + 2p
  ForceDenominatorPositive(&Rat3);
  BigIntChSign(&Rat3.numerator);                         // u = -(4S^2 + 2p)
  // According to Quartic function Wikipedia article the roots are:
  // x1,2 = -(b/4a) - S +/- sqrt(u + q/S)
  // x3,4 = -(b/4a) + S +/- sqrt(u - q/S)
  CopyBigInt(&Rat1.numerator, &RatDeprLinear.numerator);
  CopyBigInt(&Rat1.denominator, &RatDeprLinear.denominator);
  CopyBigInt(&Rat2.numerator, &RatS.denominator);
  CopyBigInt(&Rat2.denominator, &RatS.numerator);
  MultiplyRationalBySqrtRational(&Rat1, &Rat2);          // q/S
  RatDeprLinear.numerator.sign = SIGN_POSITIVE;
  for (int ctr = 0; ctr < 4; ctr++)
  {
    showX(multiplicity);
    showFirstTermQuarticEq(ctr);                         // Show -b/4a and next sign.
    *ptrOutput = ' ';
    ptrOutput++;
    if (RatS.numerator.sign == SIGN_POSITIVE)
    {             // S is real.
      bool isImaginary;
      showSquareRootOfRational(&RatS, 2, ptrTimes);
      if ((ctr == 0) || (ctr == 1))
      {   // Change sign so we always consider q/S.
        BigIntChSign(&Rat1.numerator);
      }
      // Determine whether the radicand is positive or not.
      if ((Rat1.numerator.sign == SIGN_POSITIVE) && (Rat3.numerator.sign == SIGN_POSITIVE))
      {
        isImaginary = false;
      }
      else if ((Rat1.numerator.sign == SIGN_NEGATIVE) && (Rat3.numerator.sign == SIGN_NEGATIVE))
      {
        isImaginary = true;
      }
      else
      {
        // Compute Rat5 as Rat1^2*Rat2 - Rat3^2
        BigRationalMultiply(&Rat1, &Rat1, &Rat5);
        BigRationalMultiply(&Rat5, &Rat2, &Rat5);
        BigRationalMultiply(&Rat3, &Rat3, &Rat4);
        BigRationalSubt(&Rat5, &Rat4, &Rat5);
        ForceDenominatorPositive(&Rat5);
        if ((Rat1.numerator.sign == SIGN_POSITIVE) && (Rat3.numerator.sign == SIGN_NEGATIVE))
        {
          isImaginary = (Rat5.numerator.sign != SIGN_POSITIVE);
        }
        else
        {
          isImaginary = (Rat5.numerator.sign == SIGN_POSITIVE);
        }
      }
      showPlusSignOn((ctr == 0) || (ctr == 2), TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
      if (isImaginary)
      {
        showText(ptrI);
        if (pretty != PRETTY_PRINT)
        {
          *ptrOutput = ' ';
          ptrOutput++;
        }
        showText(ptrTimes);
        *ptrOutput = ' ';
        ptrOutput++;
        BigIntChSign(&Rat1.numerator);
        BigIntChSign(&Rat3.numerator);
      }
      startSqrt();
      BigRationalDivideByInt(&Rat1, 4, &Rat1);
      BigRationalDivideByInt(&Rat3, 4, &Rat3);
      ShowRationalAndSqrParts(&Rat1, &Rat2, 2, ptrTimes);
      showPlusSignOn(Rat3.numerator.sign == SIGN_POSITIVE, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
      if (Rat3.numerator.sign == SIGN_POSITIVE)
      {
        showRationalNoParen(&Rat3);
      }
      else
      {
        BigIntChSign(&Rat3.numerator);
        showRationalNoParen(&Rat3);
        BigIntChSign(&Rat3.numerator);
      }
      BigRationalMultiplyByInt(&Rat1, 4, &Rat1);
      BigRationalMultiplyByInt(&Rat3, 4, &Rat3);
      if (ctr <= 1)
      {   // Restore sign of q/S.
        BigIntChSign(&Rat1.numerator);
      }
      endSqrt();
      if (isImaginary)
      {
        BigIntChSign(&Rat1.numerator);
        BigIntChSign(&Rat3.numerator);
      }
    }
    else
    {             // S is imaginary.
      // k^2 = sqrt(u^2 + q^2/S^2)
      // x1,2 = -(b/4a) +/- sqrt((k+u)/2) - i*(-S) +/- sqrt((k-u)/2))
      // x3,4 = -(b/4a) +/- sqrt((k+u)/2) + i*(-S) +/- sqrt((k-u)/2))
      // Get value of k^2.
      char szMinus[10];
      char* pszMinus = &szMinus[1];
      szMinus[0] = ' ';
      copyStr(&pszMinus, ptrMinus);
      *pszMinus = ' ';
      *(pszMinus + 1) = 0;

      BigRationalMultiply(&RatDeprLinear, &RatDeprLinear, &Rat4);
      // Compute q^2 / S^2 (S^2 negative)
      BigRationalDivide(&Rat4, &RatS, &Rat4);
      // Compute u^2
      BigRationalMultiply(&Rat3, &Rat3, &Rat5);
      // Compute k^2 as u^2 + q^2 / |S^2|
      BigRationalSubt(&Rat5, &Rat4, &Rat4);
      showSquareRootOfComplex(" + ", szMinus);
      showText(" + ");
      showText(ptrI);
      showText(ptrTimes);
      startParen();
      if ((ctr == 1) || (ctr == 3))
      {
        showText(ptrMinus);
        *ptrOutput = ' ';
        ptrOutput++;
      }
      BigIntChSign(&RatS.numerator);
      showSquareRootOfRational(&RatS, 2, ptrTimes);
      BigIntChSign(&RatS.numerator);
      showPlusSignOn((ctr == 1) || (ctr == 2), TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
      showSquareRootOfComplex(szMinus, " + ");
      endParen();
    }
  }
}

void QuarticEquation(const int* polynomial, int multiplicity)
{
  int ctr;
  bool isImaginary;
  int* ptrValues;
  enum eSign sign1;
  enum eSign sign2;
  const int* ptrPolynomial = polynomial;
  UncompressBigIntegerB(ptrPolynomial, &Independent);
  ptrPolynomial += numLimbs(ptrPolynomial);
  ptrPolynomial++;
  UncompressBigIntegerB(ptrPolynomial, &Linear);
  ptrPolynomial += numLimbs(ptrPolynomial);
  ptrPolynomial++;
  UncompressBigIntegerB(ptrPolynomial, &Quadratic);
  ptrPolynomial += numLimbs(ptrPolynomial);
  ptrPolynomial++;
  UncompressBigIntegerB(ptrPolynomial, &Cubic);
  ptrPolynomial += numLimbs(ptrPolynomial);
  ptrPolynomial++;
  UncompressBigIntegerB(ptrPolynomial, &Quartic);
  // Get rational coefficients of monic equation.
  CopyBigInt(&RatCubic.numerator, &Cubic);
  CopyBigInt(&RatCubic.denominator, &Quartic);
  CopyBigInt(&RatQuadratic.numerator, &Quadratic);
  CopyBigInt(&RatQuadratic.denominator, &Quartic);
  CopyBigInt(&RatLinear.numerator, &Linear);
  CopyBigInt(&RatLinear.denominator, &Quartic);
  CopyBigInt(&RatIndependent.numerator, &Independent);
  CopyBigInt(&RatIndependent.denominator, &Quartic);
  // Compute coefficients of depressed equation x^4 + px^2 + qx + r
  // where: p = (8c - 3b^2)/8, q = (b^3 - 4bc + 8d)/8,
  // r = (-3b^4 + 256e - 64bd + 16b^2*c)/256
  BigRationalMultiply(&RatCubic, &RatCubic, &RatDeprQuadratic);       // b^2
  BigRationalMultiply(&RatCubic, &RatDeprQuadratic, &RatDeprLinear);  // b^3
  BigRationalMultiply(&RatCubic, &RatQuadratic, &RatDeprIndependent); // bc
  BigRationalMultiplyByInt(&RatDeprIndependent, 16, &RatDeprIndependent);  // 16bc
  BigRationalMultiplyByInt(&RatLinear, -64, &Rat1);                   // -64d
  BigRationalAdd(&RatDeprIndependent, &Rat1, &RatDeprIndependent);    // 16bc - 64d
  BigRationalMultiply(&RatDeprIndependent, &RatCubic, &RatDeprIndependent);   // 16b^2*c - 64bd
  BigRationalMultiply(&RatDeprQuadratic, &RatDeprQuadratic, &Rat1);   // b^4
  BigRationalMultiplyByInt(&Rat1, -3, &Rat1);                         // -3b^4
  BigRationalAdd(&RatDeprIndependent, &Rat1, &RatDeprIndependent);    // -3b^4 + 16b^2*c - 64bd
  BigRationalDivideByInt(&RatDeprIndependent, 256, &RatDeprIndependent); // (-3b^4 + 16b^2*c - 64bd)/256
  BigRationalAdd(&RatDeprIndependent, &RatIndependent, &RatDeprIndependent); // (3b^4 + 16b^2*c - 64bd)/256 + e

  BigRationalMultiplyByInt(&RatDeprQuadratic, 3, &RatDeprQuadratic);  // 3b^2
  BigRationalDivideByInt(&RatDeprQuadratic, 8, &RatDeprQuadratic);    // 3b^2/8
  BigRationalSubt(&RatQuadratic, &RatDeprQuadratic, &RatDeprQuadratic); // c - 3b^2/8
  ForceDenominatorPositive(&RatDeprQuadratic);
  BigRationalMultiply(&RatCubic, &RatQuadratic, &Rat1);               // bc
  BigRationalMultiplyByInt(&Rat1, 4, &Rat1);                          // 4bc
  BigRationalSubt(&RatDeprLinear, &Rat1, &RatDeprLinear);             // b^3 - 4bc
  BigRationalDivideByInt(&RatDeprLinear, 8, &RatDeprLinear);          // (b^3 - 4bc)/8
  BigRationalAdd(&RatDeprLinear, &RatLinear, &RatDeprLinear);         // (b^3 - 4bc)/8 + d
  // Compute delta0 = c^2 - 3bd + 12e and delta1 = 2c^3 - 9bcd + 27b^2*e + 27d^2 - 72ce
  // delta1 = c(2c^2 - 9bd - 72e) + 27(b^2*e + d^2)
  BigRationalMultiply(&RatQuadratic, &RatQuadratic, &RatDelta0);      // c^2
  BigRationalMultiplyByInt(&RatDelta0, 2, &RatDelta1);                // 2c^2
  BigRationalMultiply(&RatLinear, &RatCubic, &Rat1);                  // bd
  BigRationalMultiplyByInt(&Rat1, 3, &Rat1);                          // 3bd
  BigRationalSubt(&RatDelta0, &Rat1, &RatDelta0);                     // c^2 - 3bd
  BigRationalMultiplyByInt(&Rat1, 3, &Rat1);                          // 9bd
  BigRationalSubt(&RatDelta1, &Rat1, &RatDelta1);                     // 2c^2 - 9bd
  BigRationalMultiplyByInt(&RatIndependent, 12, &Rat1);               // 12e
  BigRationalAdd(&RatDelta0, &Rat1, &RatDelta0);                      // c^2 - 3bd + 12e
  BigRationalMultiplyByInt(&RatIndependent, 72, &Rat1);               // 72e
  BigRationalSubt(&RatDelta1, &Rat1, &RatDelta1);                     // 2c^2 - 9bd - 72e
  BigRationalMultiply(&RatDelta1, &RatQuadratic, &RatDelta1);         // c(2c^2 - 9bd - 72e)
  BigRationalMultiply(&RatCubic, &RatCubic, &Rat1);                   // b^2
  BigRationalMultiply(&Rat1, &RatIndependent, &Rat1);                 // b^2*e
  BigRationalMultiply(&RatLinear, &RatLinear, &Rat2);                 // d^2
  BigRationalAdd(&Rat1, &Rat2, &Rat1);                                // b^2*e + d^2
  BigRationalMultiplyByInt(&Rat1, 27, &Rat1);                         // 27(b^2*e + d^2)
  BigRationalAdd(&RatDelta1, &Rat1, &RatDelta1);                      // c(2c^2 - 9bd - 72e) + 27(b^2*e + d^2)
  // Compute discriminant (delta)/(-27) = delta1^2 - 4 delta0^3.
  BigRationalMultiply(&RatDelta1, &RatDelta1, &RatDiscr);             // delta1^2
  BigRationalMultiply(&RatDelta0, &RatDelta0, &Rat1);                 // delta0^2
  BigRationalMultiply(&Rat1, &RatDelta0, &Rat1);                      // delta0^3
  BigRationalMultiplyByInt(&Rat1, 4, &Rat1);                          // 4 delta0^3
  BigRationalSubt(&RatDiscr, &Rat1, &RatDiscr);                       // delta1^2 - 4 delta0^3
  ForceDenominatorPositive(&RatDiscr);
  // Compute D = 64e - 16c^2 + 16b^2*c - 16bd - 3b^4
  BigRationalMultiplyByInt(&RatIndependent, 64, &RatD);               // 64e
  BigRationalMultiply(&RatCubic, &RatCubic, &Rat1);                   // b^2
  BigRationalSubt(&Rat1, &RatQuadratic, &Rat1);                       // b^2 - c
  BigRationalMultiply(&Rat1, &RatQuadratic, &Rat1);                   // b^2*c - c^2
  BigRationalMultiply(&RatCubic, &RatLinear, &Rat2);                  // bd
  BigRationalSubt(&Rat1, &Rat2, &Rat1);                               // b^2*c - c - bd
  BigRationalMultiplyByInt(&Rat1, 16, &Rat1);                         // 16(b^2*c - c - bd)
  BigRationalAdd(&RatD, &Rat1, &RatD);                                // 64e - 16c^2 + 16b^2*c - 16bd
  BigRationalMultiply(&RatCubic, &RatCubic, &Rat2);                   // b^2
  BigRationalMultiply(&Rat2, &Rat2, &Rat2);                           // b^4
  BigRationalMultiplyByInt(&Rat2, 3, &Rat2);                          // 3b^4
  BigRationalSubt(&RatD, &Rat2, &RatD);                               // D
  ForceDenominatorPositive(&RatD);

  BigRationalMultiplyByInt(&Rat1, 16, &Rat1);                         // 16c^2
  BigRationalSubt(&RatD, &Rat1, &RatD);                               // 64e - 16c^2

  BigRationalDivideByInt(&RatCubic, -4, &RatCubic);                   // -b/4
  ForceDenominatorPositive(&RatCubic);
  if (BigIntIsZero(&RatDeprLinear.numerator))
  {             // Biquadratic equation. No cube root needed in this case.
    biquadraticEquation(multiplicity);
    return;
  }
  // If resolvent equation can be factored, no cube roots are needed.
  // Lodovico Ferrari's resolvent equation is:
  // 8x^3 + 8px^2 + (2p^2-8r)x - q^2 = 0
  BigRationalMultiply(&RatDeprQuadratic, &RatDeprQuadratic, &Rat2);
  BigRationalMultiplyByInt(&Rat2, 2, &Rat2);                  // 2p^2
  BigRationalMultiplyByInt(&RatDeprIndependent, -8, &Rat3);   // -8r
  BigRationalAdd(&Rat2, &Rat3, &Rat2);
  BigRationalMultiplyByInt(&RatDeprQuadratic, 8, &Rat1);      // 8p
  BigRationalMultiply(&RatDeprLinear, &RatDeprLinear, &Rat3); // q^2
  BigRationalNegate(&Rat3, &Rat3);                            // -q^2
  // Convert from rational to integer coefficients.
  // Cubic coefficient.
  multint(&tmp3, &Rat1.denominator, 8);
  (void)BigIntMultiply(&tmp3, &Rat2.denominator, &tmp3);
  (void)BigIntMultiply(&tmp3, &Rat3.denominator, &tmp3);
  // Quadratic coefficient.
  (void)BigIntMultiply(&Rat1.numerator, &Rat2.denominator, &tmp2);
  (void)BigIntMultiply(&tmp2, &Rat3.denominator, &tmp2);
  // Linear coefficient.
  (void)BigIntMultiply(&Rat2.numerator, &Rat1.denominator, &tmp1);
  (void)BigIntMultiply(&tmp1, &Rat3.denominator, &tmp1);
  // Independent term.
  (void)BigIntMultiply(&Rat3.numerator, &Rat1.denominator, &tmp0);
  (void)BigIntMultiply(&tmp0, &Rat2.denominator, &tmp0);
  values[0] = 3;
  ptrValues = &values[1];
  NumberLength = tmp0.nbrLimbs;
  BigInteger2IntArray(ptrValues, &tmp0);
  ptrValues += numLimbs(ptrValues);
  ptrValues++;
  NumberLength = tmp1.nbrLimbs;
  BigInteger2IntArray(ptrValues, &tmp1);
  ptrValues += numLimbs(ptrValues);
  ptrValues++;
  NumberLength = tmp2.nbrLimbs;
  BigInteger2IntArray(ptrValues, &tmp2);
  ptrValues += numLimbs(ptrValues);
  ptrValues++;
  NumberLength = tmp3.nbrLimbs;
  BigInteger2IntArray(ptrValues, &tmp3);
  (void)FactorPolyOverIntegers();
  if (factorInfoInteger[0].degree == 1)
  {   // Rational root found. Get root.
    FerrariResolventHasRationalRoot(multiplicity);
  }
  else if (RatDiscr.numerator.sign == SIGN_POSITIVE)
  {
    ForceDenominatorPositive(&RatDelta1);
    if (RatDelta1.numerator.sign == SIGN_POSITIVE)
    {
      intToBigInteger(&Rat1.numerator, 1);
    }
    else
    {
      intToBigInteger(&Rat1.numerator, -1);
    }
    intToBigInteger(&Rat1.denominator, 2);
    CopyBigInt(&Rat2.numerator, &RatDiscr.numerator);
    CopyBigInt(&Rat2.denominator, &RatDiscr.denominator);
    MultiplyRationalBySqrtRational(&Rat1, &Rat2);
    startLine();
    if (pretty == PRETTY_PRINT)
    {
      showText("<var>Q</var> = ");
    }
    else
    {
      showText("Q = ");
    }
    startCbrt();
    if ((Rat2.numerator.nbrLimbs == 1) && (Rat2.numerator.limbs[0].x == 1))
    {
      BigRationalAdd(&RatDelta1, &Rat2, &Rat1);
      BigRationalDivideByInt(&Rat1, 2, &Rat1);
      showRationalNoParen(&Rat1);
    }
    else
    {
      BigRationalDivideByInt(&RatDelta1, 2, &RatDelta1);
      showRationalNoParen(&RatDelta1);
      BigRationalMultiplyByInt(&RatDelta1, 2, &RatDelta1);
      showText(" + ");
      ShowRationalAndSqrParts(&Rat1, &Rat2, 2, ptrTimes);
    }
    endCbrt();
    endLine();
    startLine();
    if (pretty == PRETTY_PRINT)
    {
      showText("<var>S</var> = ");
    }
    else
    {
      showText("S = ");
    }
    if (pretty != PARI_GP)
    {
      showRatConstants("1", "2");
    }
    else
    {
      showText("(1/2)");
    }
    showText(ptrTimes);
    startSqrt();
    intToBigInteger(&Rat1.numerator, -2);
    intToBigInteger(&Rat1.denominator, 3);
    BigRationalMultiply(&RatDeprQuadratic, &Rat1, &Rat1);
    if (!BigIntIsZero(&Rat1.numerator))
    {
      ForceDenominatorPositive(&Rat1);
      showRationalNoParen(&Rat1);
      showText(" + ");
    }
    if (pretty != PARI_GP)
    {
      showRatConstants("<var>Q</var>", "3");
    }
    else
    {
      showText("Q / 3 ");
    }
    BigRationalDivideByInt(&RatDelta0, 3, &Rat1);
    ForceDenominatorPositive(&Rat1);
    showPlusSignOn(Rat1.numerator.sign == SIGN_POSITIVE, TYPE_PM_SPACE_AFTER);
    CopyBigInt(&Rat2.numerator, &Rat1.numerator);
    CopyBigInt(&Rat2.denominator, &Rat1.denominator);
    Rat2.numerator.sign = SIGN_POSITIVE;
    showRationalOverStr(&Rat2, "Q", ptrTimes);
    endSqrt();
    endLine();
    sign1 = RatDeprLinear.numerator.sign;
    RatDeprLinear.numerator.sign = SIGN_POSITIVE;
    for (ctr = 0; ctr < 4; ctr++)
    {
      bool signS = (ctr <= 1);
      bool signNumerator = RatDeprLinear.numerator.sign == SIGN_POSITIVE;
      isImaginary = false;
      if ((signS && signNumerator) || (!signS && !signNumerator))
      {
        isImaginary = true;
      }
      showX(multiplicity);
      showFirstTermQuarticEq(ctr);
      showText((pretty == PRETTY_PRINT) ? " <var>S</var> " : " S ");
      showPlusSignOn((ctr == 0) || (ctr == 2), TYPE_PM_SPACE_AFTER);
      if (pretty != PARI_GP)
      {
        showRatConstants((isImaginary ? "i" : "1"), "2");
      }
      else
      {
        showText(isImaginary ? "(I/2) " : "(1/2) ");
      }
      showText(ptrTimes);
      *ptrOutput = ' ';
      ptrOutput++;
      startSqrt();
      if (!isImaginary)
      {
        showText(ptrMinus);
      }
      showText("4 ");
      showText(ptrTimes);
      if (pretty == PRETTY_PRINT)
      {
        showText("<var>S</var>&sup2; ");
      }
      else
      {
        showText("S^2 ");
      }
      if (!BigIntIsZero(&RatDeprQuadratic.numerator))
      {
        if (RatDeprQuadratic.numerator.sign == SIGN_NEGATIVE)
        {
          showPlusSignOn(!isImaginary, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
        }
        else
        {
          showPlusSignOn(isImaginary, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
        }
        BigRationalMultiplyByInt(&RatDeprQuadratic, 2, &Rat1);  // 2p
        Rat1.numerator.sign = SIGN_POSITIVE;
        showRationalNoParen(&Rat1);
      }
      if (!BigIntIsZero(&RatDeprLinear.numerator))
      {
        if ((signS && (sign1 == SIGN_NEGATIVE)) || (!signS && (sign1 != SIGN_NEGATIVE)))
        {
          showPlusSignOn(!isImaginary, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
        }
        else
        {
          showPlusSignOn(isImaginary, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
        }
        showRationalOverStr(&RatDeprLinear, "S", ptrTimes);
      }
      endSqrt();
    }
  }
  else
  {
    startLine();
    if (pretty == PRETTY_PRINT)
    {
      showText("<var>t</var> = arccos");
    }
    else if (pretty == TEX)
    {
      showText("t = \\arccos");
    }
    else
    {
      showText("t = acos");
    }
    startParen();
    BigRationalDivide(&RatDelta1, &RatDelta0, &Rat1);
    BigRationalDivideByInt(&Rat1, 2, &Rat1);
    CopyBigInt(&Rat2.numerator, &RatDelta0.denominator);
    CopyBigInt(&Rat2.denominator, &RatDelta0.numerator);
    MultiplyRationalBySqrtRational(&Rat1, &Rat2);
    ShowRationalAndSqrParts(&Rat1, &Rat2, 2, ptrTimes);
    endParen();
    endLine();
    startLine();
    if (pretty == PRETTY_PRINT)
    {
      showText("<var>S</var> = ");
    }
    else
    {
      showText("S = ");
    }
    if (pretty != PARI_GP)
    {
      showRatConstants("1", "2");
    }
    else
    {
      showText("(1/2)");
    }
    showText(ptrTimes);
    startSqrt();
    intToBigInteger(&Rat1.numerator, -2);
    intToBigInteger(&Rat1.denominator, 3);
    BigRationalMultiply(&RatDeprQuadratic, &Rat1, &Rat2);
    ForceDenominatorPositive(&Rat2);
    if (!BigIntIsZero(&Rat2.numerator))
    {
      showRationalNoParen(&Rat2);
      showText(" + ");
    }
    BigIntChSign(&Rat1.numerator);
    CopyBigInt(&Rat2.numerator, &RatDelta0.numerator);
    CopyBigInt(&Rat2.denominator, &RatDelta0.denominator);
    MultiplyRationalBySqrtRational(&Rat1, &Rat2);
    ForceDenominatorPositive(&Rat1);
    ShowRationalAndSqrParts(&Rat1, &Rat2, 2, ptrTimes);
    showText(ptrTimes);
    showText(ptrCos);
    if (pretty != PARI_GP)
    {
      showRatConstants("<var>t</var>", "3");
    }
    else
    {
      showText("(<var>t</var> / 3)");
    }
    endSqrt();
    if (pretty == TEX)
    {
      *ptrOutput = '}';
      ptrOutput++;
    }
    endLine();
    BigRationalMultiplyByInt(&RatDeprQuadratic, -2, &Rat1);
    ForceDenominatorPositive(&Rat1);
    Rat1.numerator.sign = SIGN_POSITIVE;
    CopyBigInt(&Rat2.numerator, &RatDeprLinear.numerator);
    CopyBigInt(&Rat2.denominator, &RatDeprLinear.denominator);
    sign2 = Rat2.numerator.sign;
    Rat2.numerator.sign = SIGN_POSITIVE;
    isImaginary = (RatDeprQuadratic.numerator.sign == SIGN_POSITIVE) ||
      (RatD.numerator.sign == SIGN_POSITIVE);
    RatDeprLinear.numerator.sign = SIGN_POSITIVE;
    for (ctr = 0; ctr < 4; ctr++)
    {
      showX(multiplicity);
      showFirstTermQuarticEq(ctr);
      showText(" <var>S</var> ");
      showPlusSignOn((ctr == 0) || (ctr == 2), TYPE_PM_SPACE_AFTER);
      if (pretty != PARI_GP)
      {
        showRatConstants((isImaginary ? "i" : "1"), "2");
      }
      else
      {
        showText(isImaginary ? "(I/2)" : "(1/2)");
      }
      showText(ptrTimes);
      startSqrt();
      if (!isImaginary)
      {
        showText("&minus;");
      }
      showText("4 ");
      showText(ptrTimes);
      showText("<var>S</var>");
      showText((pretty == PRETTY_PRINT) ? "&sup2;" : "^2");
      *ptrOutput = ' ';
      ptrOutput++;
      if (!BigIntIsZero(&RatDeprQuadratic.numerator))
      {
        if (RatDeprQuadratic.numerator.sign == SIGN_NEGATIVE)
        {
          showPlusSignOn(!isImaginary, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
        }
        else
        {
          showPlusSignOn(isImaginary, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
        }
        BigRationalMultiplyByInt(&RatDeprQuadratic, 2, &Rat1);  // 2p
        Rat1.numerator.sign = SIGN_POSITIVE;
        showRationalNoParen(&Rat1);
      }
      if (!BigIntIsZero(&RatDeprLinear.numerator))
      {
        bool signS = (ctr <= 1);
        if ((signS && (sign2 == SIGN_NEGATIVE)) || (!signS && (sign2 != SIGN_NEGATIVE)))
        {
          showPlusSignOn(!isImaginary, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
        }
        else
        {
          showPlusSignOn(isImaginary, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
        }
        showRationalOverStr(&RatDeprLinear, "S", ptrTimes);
      }
      endSqrt();
      endLine();
    }
  }
}


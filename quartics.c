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
#include "string/strings.h"
#include "rootseq.h"
#include "expression.h"
#include "copyStr.h"

#define NBR_COEFFS_BACKUP 7
extern const char* ptrI;
extern const char* ptrCos;
extern char* ptrOutput;
extern int eqNbr;
char currLetter;
static int eqPlusMinusAfterRoot;
static int FerrariEq;
static BigRational coeffBackup[NBR_COEFFS_BACKUP];

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

// Show steps for solving the equation x^4 + mx^2 + n = 0
static void showStepsForRealBiquadratic(BigRational *pQuadratic,
  BigRational *pIndependent, char letter)
{
  BigRationalDivideByInt(pQuadratic, 2, &Rat1);
  BigRationalMultiply(pQuadratic, pQuadratic, &Rat2);
  BigRationalDivideByInt(&Rat2, 4, &Rat2);
  showText("<p>");
  // Since 
  showText(LITERAL_STEPS_BIQUADR1);
  startParen();
  showRatCoeffAndPowerVar(NULL, -2, letter);
  showPlusMinusRational(&Rat1);
  endParen();
  showPower(&ptrOutput, 2);
  showText(" = ");
  showRatCoeffAndPowerVar(NULL, -4, letter);
  showRatCoeffAndPowerVar(pQuadratic, 2, letter);
  showRatCoeffAndPowerVar(&Rat2, 0, letter);
  showText("</p><p>");
  // we get
  showText(LITERAL_STEPS_BIQUADR2);
  startParen();
  showRatCoeffAndPowerVar(NULL, -2, letter);
  showPlusMinusRational(&Rat1);
  endParen();
  showPower(&ptrOutput, 2);
  BigRationalMultiplyByInt(&Rat2, -1, &Rat2);
  showRatCoeffAndPowerVar(&Rat2, 0, letter);
  showRatCoeffAndPowerVar(pIndependent, 0, letter);
  showText(" = 0</p><p>");
  showRatCoeffAndPowerVar(NULL, -2, letter);
  showText(" = ");
  BigRationalMultiplyByInt(&Rat1, -1, &Rat1);
  showRationalNoParen(&Rat1);
  showText(ptrPlusMinus);
  BigRationalMultiply(pQuadratic, pQuadratic, &Rat2);
  BigRationalDivideByInt(&Rat2, 4, &Rat2);
  BigRationalSubt(&Rat2, pIndependent, &Rat2);
  showSquareRootOfRational(&Rat2, 2, ptrTimes);
  showText("</p>");
}

static void showM(void)
{
  showVariable(&ptrOutput, 'm');
}

static void showU(void)
{
  showVariable(&ptrOutput, 'u');
}

static void showAdjustForCubic(void)
{
  if (!BigIntIsZero(&RatCubic.numerator))
  {
    showText("<p>");
    // Since
    showText(LITERAL_STEPS_BIQUADR1);
    showVariable(&ptrOutput, 'x');
    showText(" = ");
    showVariable(&ptrOutput, 'y');
    // RatCubic is already divided by -4.
    showPlusMinusRational(&RatCubic);
    // we get
    showText(LITERAL_STEPS_BIQUADR2);
  }
}

// Show steps to find the square root of the complex number Rat1 + i*sqrt(-Rat2).
static void showStepsForComplexSquareRoot(char letter, const char *pszPlus)
{
  showText("<p>");
  showRatCoeffAndPowerVar(NULL, -2, letter);
  showText(" = ");
  showRationalNoParen(&Rat1);
  showText(pszPlus);
  showText(" ");
  showText(ptrI);
  showText(" ");
  showText(ptrTimes);
  BigRationalMultiplyByInt(&Rat2, -1, &Rat2);
  showSquareRootOfRational(&Rat2, 2, ptrTimes);
  generateEqNbr();     // Equation 1.
  showText("</p><p>");
  // Let
  showText(LITERAL_STEPS_COMPLEX_SQROOT1);
  showRatCoeffAndPowerVar(NULL, -1, letter);
  showText(" = ");
  showU();
  showText(" + ");
  showText(ptrI);
  showRatCoeffAndPowerVar(NULL, -1, 'v');
  showText("</p><p>");
  showText(LITERAL_STEPS_COMPLEX_SQROOT2);
  showRatCoeffAndPowerVar(NULL, -2, letter);
  showText(" = ");
  showRatCoeffAndPowerVar(NULL, -2, 'u');
  showText(ptrMinus);
  showRatCoeffAndPowerVar(NULL, -2, 'v');
  showText(" + 2");
  showText(ptrTimes);
  showText(ptrI);
  showU();
  showText(ptrTimes);
  showRatCoeffAndPowerVar(NULL, -1, 'v');
  generateEqNbr();     // Equation 2.
  showText("</p><p>");
  // Equating the real parts of $1q and $2q:  (equations 1 and 2)
  formatString(&ptrOutput, LITERAL_STEPS_COMPLEX_SQROOT3, eqNbr - 1, eqNbr);
  showText("</p><p>");
  showRatCoeffAndPowerVar(NULL, -2, 'u');
  showText(ptrMinus);
  showRatCoeffAndPowerVar(NULL, -2, 'v');
  showText(" = ");
  showRationalNoParen(&Rat1);
  generateEqNbr();     // Equation 3.
  showText("</p><p>");
  // Multiplying by
  showText(LITERAL_STEPS_COMPLEX_SQROOT4);
  showRatCoeffAndPowerVar(NULL, -2, 'u');
  showText(":</p><p>");
  BigRationalMultiplyByInt(&Rat1, -1, &Rat1);
  showRatCoeffAndPowerVar(NULL, -4, 'u');
  showRatCoeffAndPowerVar(&Rat1, 2, 'u');
  showText(ptrMinus);
  showRatCoeffAndPowerVar(NULL, -2, 'u');
  showText(ptrTimes);
  showRatCoeffAndPowerVar(NULL, -2, 'v');
  showText(" = 0");
  generateEqNbr();      // Equation 4.
  showText("</p><p>");
  // Equating the imaginary parts of $1q and $2q (equations 1 and 2) and dividing by 2:
  formatString(&ptrOutput, LITERAL_STEPS_COMPLEX_SQROOT5, eqNbr - 3, eqNbr - 2);
  showText("</p><p>");
  showU();
  showText(ptrTimes);
  showRatCoeffAndPowerVar(NULL, -1, 'v');
  showText(" = ");
  if (pszPlus == ptrMinus)
  {
    showText(ptrMinus);
  }
  BigRationalDivideByInt(&Rat2, 4, &Rat2);
  showSquareRootOfRational(&Rat2, 2, ptrTimes);
  generateEqNbr();      // Equation 5.
  showText("</p><p>");
  showRatCoeffAndPowerVar(NULL, -2, 'u');
  showText(ptrTimes);
  showRatCoeffAndPowerVar(NULL, -2, 'v');
  showText(" = ");
  showRational(&Rat2);
  generateEqNbr();      // Equation 6.
  showText("</p><p>");
  // From $1q y $2q: (Equations 4 and 6)
  formatString(&ptrOutput, LITERAL_STEPS_COMPLEX_SQROOT6, eqNbr - 2, eqNbr);
  showText("</p><p>");
  showRatCoeffAndPowerVar(NULL, -4, 'u');
  showRatCoeffAndPowerVar(&Rat1, 2, 'u');
  BigRationalMultiplyByInt(&Rat2, -1, &Rat4);
  showRatCoeffAndPowerVar(&Rat4, 0, 'u');
  showText(" = 0</p>");
  CopyBigInt(&Rat3.numerator, &Rat1.numerator);
  CopyBigInt(&Rat3.denominator, &Rat1.denominator);
  showStepsForRealBiquadratic(&Rat3, &Rat4, 'u');
  // We will use the plus sign so the argument of the square root is positive.
  formatString(&ptrOutput, "<p>$1s</p></p>", LITERAL_STEPS_COMPLEX_SQROOT7);
  // From
  formatString(&ptrOutput, LITERAL_STEPS_COMPLEX_SQROOT8, eqNbr - 3);   // Equation 3.
  showText("</p><p>");
  showRatCoeffAndPowerVar(NULL, -2, 'v');
  showText(" = ");
  BigRationalMultiplyByInt(&Rat1, -2, &Rat1);
  showRationalNoParen(&Rat1);
  showText(" + ");
  showRatCoeffAndPowerVar(NULL, -2, 'u');
  showText(" = ");
  BigRationalDivideByInt(&Rat1, 2, &Rat1);
  showRationalNoParen(&Rat1);
  showText(" + ");
  showSquareRootOfRational(&Rat2, 2, ptrTimes);
  showText("</p><p>");
  // Since
  showText(LITERAL_STEPS_BIQUADR1);
  showU();
  // is the real part of
  showText(LITERAL_STEPS_COMPLEX_SQROOT9);
  showRatCoeffAndPowerVar(NULL, -1, letter);
  // and
  showText(LITERAL_STEPS_COMPLEX_SQROOT10);
  showRatCoeffAndPowerVar(NULL, -1, 'v');
  // is the imaginary part of
  showText(LITERAL_STEPS_COMPLEX_SQROOT11);
  showRatCoeffAndPowerVar(NULL, -1, letter);
  // , we get:
  showText(LITERAL_STEPS_COMPLEX_SQROOT12);
  showText("</p>");
}

static void showStepsForBiquadratic(void)
{
  showStepsForRealBiquadratic(&RatDeprQuadratic, &RatDeprIndependent, currLetter);
  if (RatDiscr.numerator.sign == SIGN_NEGATIVE)
  {
    showStepsForComplexSquareRoot(currLetter, ptrPlusMinus);
  }
  showAdjustForCubic();
}

// Solve equation x^4 + cx^2 + e = 0 where e is a perfect square.
// c and e are rational numbers.
static void biquadraticConstantSquare(int multiplicity)
{  // At this moment Rat3 = sqrt(e).
  BigRationalDivideByInt(&RatDeprQuadratic, 2, &RatDeprQuadratic);
  for (int ctr = 0; ctr < 4; ctr++)
  {
    enum eSign sign;
    showX(multiplicity);
    BigRationalSubt(&Rat3, &RatDeprQuadratic, &Rat1);
    BigRationalAdd(&Rat3, &RatDeprQuadratic, &Rat2);
    BigIntChSign(&Rat2.numerator);
    BigRationalDivideByInt(&Rat1, 2, &Rat1);
    BigRationalDivideByInt(&Rat2, 2, &Rat2);
    ForceDenominatorPositive(&Rat1);
    ForceDenominatorPositive(&Rat2);
    // Rat3 is positive, so Rat1 is greater than Rat2.
    // If one of Rat1 or Rat2 were square, the root would have the form
    // a + b sqrt(c), so the quartic polynomial would be reducible.
    // If both Rat1 and Rat2 were squares, the root would be a rational
    // number, so it would be reducible.
    showFirstTermQuarticEq(ctr);
    if (Rat1.numerator.sign == SIGN_NEGATIVE)
    {
      showText(ptrI);
      showText(" ");
      showText(ptrTimes);
      showText(" (");
    }
    sign = Rat1.numerator.sign;             // Back up sign.
    Rat1.numerator.sign = SIGN_POSITIVE;
    showSquareRootOfRational(&Rat1, 2, ptrTimes);
    Rat1.numerator.sign = sign;             // Restore sign.
    showPlusSignOn((ctr == 0) || (ctr == 2), TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
    sign = Rat2.numerator.sign;             // Back up sign.
    Rat2.numerator.sign = SIGN_POSITIVE;
    showSquareRootOfRational(&Rat2, 2, ptrTimes);
    Rat2.numerator.sign = sign;             // Restore sign.
    if (Rat1.numerator.sign == SIGN_NEGATIVE)
    {
      showText(")");
    }
    else
    {
      if (Rat2.numerator.sign == SIGN_NEGATIVE)
      {
        showText(ptrTimes);
        showText(ptrI);
      }
    }
    endShowX();
  }
}

// Solve equation x^4 + cx^2 + e = 0 where e is not a perfect square.
// c and e are rational numbers.
static void biquadraticConstantNotSquare(int multiplicity)
{
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
  enum eSign signQuadr = RatDeprQuadratic.numerator.sign;
  RatDeprQuadratic.numerator.sign = SIGN_POSITIVE;
  for (int ctr = 0; ctr < 4; ctr++)
  {
    showX(multiplicity);
    showFirstTermQuarticEq(ctr);
    if (RatDiscr.numerator.sign == SIGN_NEGATIVE)
    {  // x = +/- sqrt((sqrt(e) - c) / 2) +/- i * sqrt((sqrt(e) + c) / 2)
      for (int ctr2 = 0; ctr2 < 2; ctr2++)
      {   // ctr2 = 0: real part, ctr2 = 1: imaginary part.
        if (ctr2 == 1)
        { // Start of imaginary part.
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
          if (signQuadr == SIGN_POSITIVE)
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
      bool isXSquaredPositive;
      ForceDenominatorPositive(&RatDeprIndependent);
      // If c >= 0 and positive sqrt => c^2 < c^2 - 4e => e < 0.
      // If c > 0 and negative sqrt => x^2 is always negative.
      // If c < 0 and positive sqrt => x^2 is always positive.
      // If c <= 0 and negative sqrt => c^2 > c^2 - 4e => e > 0
      if ((ctr == 0) || (ctr == 2))
      {    // Positive sqrt
        isXSquaredPositive = ((signQuadr == SIGN_NEGATIVE) ||
          (RatDeprIndependent.numerator.sign == SIGN_NEGATIVE));
      }
      else
      {    // Negative sqrt
        isXSquaredPositive = (((signQuadr == SIGN_NEGATIVE) ||
          BigIntIsZero(&RatDeprQuadratic.numerator)) &&
          (RatDeprIndependent.numerator.sign == SIGN_POSITIVE));
      }
      if (!isXSquaredPositive)
      {    // Root is an imaginary number.
        showText(ptrI);
        showText(" ");
        showText(ptrTimes);
      }
      if (BigIntIsZero(&RatDeprQuadratic.numerator))
      {
        BigIntChSign(&RatDeprIndependent.numerator);
        showSquareRootOfRational(&RatDeprIndependent, 4, ptrTimes);
        BigIntChSign(&RatDeprIndependent.numerator);
      }
      else
      {
        startSqrt();
        CopyBigInt(&Rat1.numerator, &RatDeprQuadratic.numerator);
        CopyBigInt(&Rat1.denominator, &RatDeprQuadratic.denominator);
        if (isXSquaredPositive != signQuadr)
        {
          BigIntChSign(&Rat1.numerator);
        }
        showRationalNoParen(&Rat1);
        if ((ctr == 0) || (ctr == 2))
        {    // Positive sqrt
          showPlusSignOn(isXSquaredPositive, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
        }
        else
        {    // Negative sqrt
          showPlusSignOn(!isXSquaredPositive, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
        }
        BigRationalMultiply(&RatDeprQuadratic, &RatDeprQuadratic, &Rat1);
        BigRationalSubt(&Rat1, &RatDeprIndependent, &Rat1);
        showSquareRootOfRational(&Rat1, 2, ptrTimes);
        endSqrt();
      }
    }
    endShowX();
  }
}

// Solve equation x^4 + cx^2 + e = 0 where c and e are rational numbers.
static void biquadraticEquation(int multiplicity)
{
  if (teach)
  {
    // The cubic and linear coefficients are equal to zero,
    // so this is a biquadratic equation.
    formatString(&ptrOutput, "<p>$1s</p>", LITERAL_BIQUADR_EQ1);
  }
  // Compute discriminant as c^2-4e
  BigRationalMultiply(&RatDeprQuadratic, &RatDeprQuadratic, &Rat2);
  BigRationalMultiplyByInt(&RatDeprIndependent, 4, &Rat1);
  BigRationalSubt(&Rat2, &Rat1, &RatDiscr);
  ForceDenominatorPositive(&RatDiscr);
  if (teach)
  {
    showStepsForBiquadratic();
  }
  // For biquadratic equation, the depressed independent term equals the independent term (r = e).

  // Solutions are x = +/- sqrt((sqrt(r)-p/2)/2) +/- sqrt((-sqrt(r)-p/2)/2)
  // Test whether e is a perfect square.
  if (BigRationalSquareRoot(&RatDeprIndependent, &Rat3))
  {      // e is a perfect square. Rat3 = sqrt(r).
    biquadraticConstantSquare(multiplicity);
  }
  else
  {      // e is not a perfect square.
    biquadraticConstantNotSquare(multiplicity);
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

static void showSqrt2M(void)
{
  startSqrt();
  showText("2");
  showText(ptrTimes);
  showM();
  endSqrt();
}

static void showSquareLHS(void)
{
  startParen();
  showRatCoeffAndPowerVar(NULL, -2, currLetter);
  showRatCoeffAndPowerVar(&Rat1, 0, currLetter);
  showText(" + ");
  showM();
  endParen();
  showPower(&ptrOutput, 2);
}

static void showNbrOverSqrt2M(void)
{
  showPlusSignOn(Rat4.numerator.sign == SIGN_POSITIVE,
    TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
  if (Rat4.numerator.sign == SIGN_POSITIVE)
  {
    showRationalNoParenOverGeneric(&Rat4, showSqrt2M);
  }
  else
  {
    Rat4.numerator.sign = SIGN_POSITIVE;
    showRationalNoParenOverGeneric(&Rat4, showSqrt2M);
    Rat4.numerator.sign = SIGN_NEGATIVE;
  }
}

static void showSquareRHS(bool isSquare)
{
  startParen();
  showSqrt2M();
  showText(ptrTimes);
  showRatCoeffAndPowerVar(NULL, -1, currLetter);
  BigRationalDivideByInt(&RatDeprLinear, 2, &Rat4);
  BigRationalMultiplyByInt(&Rat4, -1, &Rat4);
  showNbrOverSqrt2M();
  endParen();
  if (isSquare)
  {
    showPower(&ptrOutput, 2);
  }
}

static void ComputeLeftSqrtFerrariRational(int divisor)
{
  BigRationalDivideByInt(&RatS, divisor, &Rat2);
  BigRationalMultiplyByInt(&Rat2, 2, &Rat2);
  Rat2.numerator.sign = SIGN_POSITIVE;
  showSquareRootOfRational(&Rat2, 2, ptrTimes);
  if (RatS.numerator.sign == SIGN_POSITIVE)
  {
    showText(ptrTimes);
    showText(ptrI);
  }
}

static void ComputeRightSqrtFerrariRational(int divisor, enum eSign sign)
{
  enum eSign signIfPositive = SIGN_POSITIVE;
  if (RatDeprLinear.numerator.sign != RatS.numerator.sign)
  {
    signIfPositive = SIGN_NEGATIVE;
  }
  showPlusSignOn(signIfPositive == sign, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
  BigRationalMultiplyByInt(&RatDeprLinear, -1, &Rat2);
  BigRationalDivideByInt(&Rat2, 2 * divisor, &Rat2);
  Rat2.numerator.sign = SIGN_POSITIVE;
  CopyBigInt(&Rat3.numerator, &RatS.denominator);
  CopyBigInt(&Rat3.denominator, &RatS.numerator);
  Rat3.denominator.sign = SIGN_POSITIVE;
  BigRationalDivideByInt(&Rat3, 2, &Rat3);
  MultiplyRationalBySqrtRational(&Rat2, &Rat3);
  ShowRationalAndSqrParts(&Rat2, &Rat3, 2, ptrTimes);
  if (RatS.numerator.sign == SIGN_POSITIVE)
  {
    showText(ptrTimes);
    showText(" ");
    showText(ptrI);
  }
}

// RatS is the rational root divided by 2.
static void showStepsOnRationalRootOfResolvent(const char *pszMinus, enum eSign sign)
{
  // Change value of RatS and RatDeprLinear to match other debug code.
  BigRationalMultiplyByInt(&RatS, -2, &RatS);
  BigRationalMultiplyByInt(&RatDeprLinear, -1, &RatDeprLinear);
  showText("<p>");
  showRatCoeffAndPowerVar(NULL, -2, currLetter);
  showText(" ");
  showText(pszMinus);
  showText(" ");
  ComputeLeftSqrtFerrariRational(1);
  showText(ptrTimes);
  showRatCoeffAndPowerVar(NULL, -1, currLetter);
  BigRationalDivideByInt(&RatDeprQuadratic, 2, &Rat1);
  BigRationalSubt(&Rat1, &RatS, &Rat1);
  showRatCoeffAndPowerVar(&Rat1, 0, currLetter);
  ComputeRightSqrtFerrariRational(1, sign == SIGN_POSITIVE ? SIGN_NEGATIVE : SIGN_POSITIVE);
  showText(" = 0</p><p>");
  // From the identity:
  showText("</p><p>");
  startParen();
  showRatCoeffAndPowerVar(NULL, -1, currLetter);
  showText(" ");
  showText(pszMinus);
  showText(" ");
  ComputeLeftSqrtFerrariRational(4);
  endParen();
  showPower(&ptrOutput, 2);
  showText(" = ");
  showRatCoeffAndPowerVar(NULL, -2, currLetter);
  showText(" ");
  showText(pszMinus);
  showText(" ");
  ComputeLeftSqrtFerrariRational(1);
  showText(ptrTimes);
  showRatCoeffAndPowerVar(NULL, -1, currLetter);
  CopyBigInt(&Rat3.numerator, &RatS.numerator);
  CopyBigInt(&Rat3.denominator, &RatS.denominator);
  BigRationalMultiplyByInt(&Rat3, -1, &Rat3);
  BigRationalDivideByInt(&Rat3, 2, &Rat3);
  showRatCoeffAndPowerVar(&Rat3, 0, currLetter);
  showText("</p><p>");
  // we get:
  showText(LITERAL_BIQUADR_EQ3);
  showText("</p><p>");
  startParen();
  showRatCoeffAndPowerVar(NULL, -1, currLetter);
  showText(" ");
  showText(pszMinus);
  showText(" ");
  ComputeLeftSqrtFerrariRational(4);
  endParen();
  showPower(&ptrOutput, 2);
  BigRationalDivideByInt(&RatDeprQuadratic, 2, &Rat1);
  BigRationalSubt(&Rat1, &RatS, &Rat1);
  CopyBigInt(&Rat3.numerator, &RatS.numerator);
  CopyBigInt(&Rat3.denominator, &RatS.denominator);
  BigRationalDivideByInt(&Rat3, 2, &Rat3);
  BigRationalAdd(&Rat1, &Rat3, &Rat1);
  showRatCoeffAndPowerVar(&Rat1, 0, currLetter);
  ComputeRightSqrtFerrariRational(1, sign == SIGN_POSITIVE ? SIGN_NEGATIVE : SIGN_POSITIVE);
  showText(" = 0</p><p>");
  showRatCoeffAndPowerVar(NULL, -1, currLetter);
  showText(" = ");
  if (sign == SIGN_NEGATIVE)
  {
    showText(ptrMinus);
    showText(" ");
  }
  ComputeLeftSqrtFerrariRational(4);
  showText(ptrPlusMinus);
  startSqrt();
  BigRationalDivideByInt(&RatDeprQuadratic, 2, &Rat1);
  BigRationalSubt(&RatS, &Rat1, &Rat1);
  CopyBigInt(&Rat3.numerator, &RatS.numerator);
  CopyBigInt(&Rat3.denominator, &RatS.denominator);
  BigRationalDivideByInt(&Rat3, 2, &Rat3);
  BigRationalSubt(&Rat1, &Rat3, &Rat1);
  showRationalNoParen(&Rat1);
  ComputeRightSqrtFerrariRational(1, sign);
  endSqrt();
  BigRationalMultiply(&Rat2, &Rat2, &Rat2);
  BigRationalMultiply(&Rat2, &Rat3, &Rat2);
  showText("</p>");
  if (RatS.numerator.sign == SIGN_POSITIVE)
  {   // Square root of complex number Rat1 +/- i*sqrt(-Rat2).
    // We must compute the square root of a complex number.
    formatString(&ptrOutput, "<p>$1s</p>", LITERAL_BIQUADR_EQ4);
    BigRationalMultiplyByInt(&Rat2, -1, &Rat2);
    if (RatDeprLinear.numerator.sign == sign)
    {
      showStepsForComplexSquareRoot('g', "+");
    }
    else
    {
      showStepsForComplexSquareRoot('g', ptrMinus);
    }
    showText("<p>");
    showVariable(&ptrOutput, currLetter);
    showText(" = ");
    if (sign == SIGN_NEGATIVE)
    {
      showText(ptrMinus);
      showText(" ");
    }
    BigRationalDivideByInt(&RatS, 2, &Rat2);
    showSquareRootOfRational(&Rat2, 2, ptrTimes);
    showText(ptrTimes);
    showText(ptrI);
    showText(" ");
    showText(ptrPlusMinus);
    showText(" ");
    showVariable(&ptrOutput, 'u');
    showText(" ");
    showText(ptrPlusMinus);
    showText(" ");
    showText(ptrI);
    showText(ptrTimes);
    showVariable(&ptrOutput, 'v');
    // From $1q, $2v y $3v have $4?1?the same sign??$4?0?different sign??.
    formatString(&ptrOutput, LITERAL_BIQUADR_EQ5, eqNbr - 1, 'u', 'v',
      (RatDeprLinear.numerator.sign == sign ? 0 : 1));
    showText("</p>");
  }
  showAdjustForCubic();
  // Restore value of RatS and RetDeprLinear.
  BigRationalMultiplyByInt(&RatDeprLinear, -1, &RatDeprLinear);
  BigRationalMultiplyByInt(&RatS, -1, &RatS);
  BigRationalDivideByInt(&RatS, 2, &RatS);
}

static void showFerrariResolventRationalRoot(void)
{
  showText("<p>");
  // This equation has a rational root:
  showText(LITERAL_FERRARI_RAT_ROOT1);
  showText("</p><p>");
  showM();
  showText(" = ");
  BigRationalMultiplyByInt(&RatS, -1, &Rat4);
  showRationalNoParen(&Rat4);
  showText("</p><p>");
  // Replacing this root in $1q:
  formatString(&ptrOutput, LITERAL_FERRARI_RAT_ROOT2, eqNbr);  // Equation 1.
  int2dec(&ptrOutput, eqNbr);    
  showText("</p><p>");
  BigRationalDivideByInt(&RatDeprQuadratic, 2, &Rat1);
  BigRationalAdd(&Rat1, &Rat4, &Rat1);
  startParen();
  showRatCoeffAndPowerVar(NULL, -2, currLetter);
  showRatCoeffAndPowerVar(&Rat1, 0, currLetter);
  endParen();
  showPower(&ptrOutput, 2);
  showText(" = ");
  startParen();
  ComputeLeftSqrtFerrariRational(1);
  showText(ptrTimes);
  showRatCoeffAndPowerVar(NULL, -1, currLetter);
  showPlusSignOn(RatDeprLinear.numerator.sign == SIGN_NEGATIVE,
    TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
  BigRationalDivideByInt(&RatDeprLinear, 2, &Rat2);
  Rat2.numerator.sign = SIGN_POSITIVE;
  CopyBigInt(&Rat3.numerator, &Rat4.denominator);
  CopyBigInt(&Rat3.denominator, &Rat4.numerator);
  BigRationalDivideByInt(&Rat3, 2, &Rat3);
  Rat3.denominator.sign = SIGN_POSITIVE;
  ShowRationalAndSqrParts(&Rat2, &Rat3, 2, ptrTimes);
  if (RatS.numerator.sign == SIGN_POSITIVE)
  {
    showText(ptrTimes);
    showText(" ");
    showRatConstants("1", ptrI);
  }
  showText(ptrTimes);
  endParen();
  showPower(&ptrOutput, 2);
  showText("</p><p>");
  showRatCoeffAndPowerVar(NULL, -2, currLetter);
  showRatCoeffAndPowerVar(&Rat1, 0, currLetter);
  showText(" = ");
  BigRationalMultiplyByInt(&Rat4, 2, &Rat2);
  showText(ptrPlusMinus);
  startParen();
  ComputeLeftSqrtFerrariRational(1);
  showText(ptrTimes);
  showRatCoeffAndPowerVar(NULL, -1, currLetter);
  ComputeRightSqrtFerrariRational(1, SIGN_POSITIVE);
  endParen();
  generateEqNbr();
  eqPlusMinusAfterRoot = eqNbr;
  showText("</p>");
}

static void showSolFerrariResolventRatRoot(int ctr, enum eSign oldLinearSign, int multiplicity)
{
  RatDeprLinear.numerator.sign = oldLinearSign;
  CopyBigInt(&Rat1.numerator, &RatDeprLinear.numerator);
  CopyBigInt(&Rat1.denominator, &RatDeprLinear.denominator);
  CopyBigInt(&Rat2.numerator, &RatS.denominator);
  CopyBigInt(&Rat2.denominator, &RatS.numerator);
  MultiplyRationalBySqrtRational(&Rat1, &Rat2);          // q/S
  BigRationalMultiplyByInt(&RatS, 4, &Rat3);             // 4S^2
  BigRationalAdd(&Rat3, &RatDeprQuadratic, &Rat3);
  BigRationalAdd(&Rat3, &RatDeprQuadratic, &Rat3);       // 4S^2 + 2p
  ForceDenominatorPositive(&Rat3);
  BigIntChSign(&Rat3.numerator);                         // u = -(4S^2 + 2p)
  RatDeprLinear.numerator.sign = SIGN_POSITIVE;
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
  endShowX();
}

static void showEndFerrariResolventNoRationalRoot(const char* minus, int multiplicity)
{
  showRatCoeffAndPowerVar(NULL, -2, currLetter);
  showText(" ");
  showText(minus);
  showText(" ");
  showSqrt2M();
  showText(ptrTimes);
  showText(" ");
  showVariable(&ptrOutput, currLetter);
  BigRationalDivideByInt(&RatDeprQuadratic, 2, &Rat1);
  showRatCoeffAndPowerVar(&Rat1, 0, currLetter);
  showText(" + ");
  showVariable(&ptrOutput, 'm');
  BigRationalDivideByInt(&RatDeprLinear, 2, &Rat4);
  if (minus != ptrMinus)
  {
    BigIntChSign(&Rat4.numerator);
  }
  showNbrOverSqrt2M();
  showText("= 0</p><p>");
  // Since
  showText(LITERAL_STEPS_BIQUADR1);
  intToBigInteger(&Rat2.numerator, (minus == ptrMinus? -1: 1));
  intToBigInteger(&Rat2.denominator, 2);
  startParen();
  showRatCoeffAndPowerVar(NULL, -1, currLetter);
  showRatCoeffAndPowerVar(&Rat2, 0, currLetter);
  showText(ptrTimes);
  showText(" ");
  showSqrt2M();
  endParen();
  showPower(&ptrOutput, 2);
  showText(" = ");
  showRatCoeffAndPowerVar(NULL, -2, currLetter);
  showText(" ");
  showText(minus);
  showText(" ");
  showSqrt2M();
  showText(ptrTimes);
  showText(" ");
  showVariable(&ptrOutput, currLetter);
  Rat2.numerator.sign = SIGN_POSITIVE;
  showRatCoeffAndPowerVar(&Rat2, 1, 'm');
  showText("</p><p>");
  // we get:
  showText(LITERAL_BIQUADR_EQ3);
  showText("</p><p>");
  intToBigInteger(&Rat2.numerator, (minus == ptrMinus ? -1 : 1));
  startParen();
  showRatCoeffAndPowerVar(NULL, -1, currLetter);
  showRatCoeffAndPowerVar(&Rat2, 0, currLetter);
  showText(ptrTimes);
  showText(" ");
  showSqrt2M();
  endParen();
  showPower(&ptrOutput, 2);
  showRatCoeffAndPowerVar(&Rat1, 0, currLetter);
  showRatCoeffAndPowerVar(&Rat2, 1, 'm');
  showNbrOverSqrt2M();
  showText(" = 0</p>");
  for (int ctr = 0; ctr < 3; ctr++)
  {
    if (ctr == 0)
    {
      showText("<p>");
      showRatCoeffAndPowerVar(NULL, -1, currLetter);
      showText(" = ");
    }
    else
    {
      showX(multiplicity);
      showPlusMinusRational(&RatCubic);
    }
    intToBigInteger(&Rat2.numerator, (minus == ptrMinus ? 1 : -1));
    showRatCoeffAndPowerVar(&Rat2, 0, currLetter);
    showText(ptrTimes);
    showText(" ");
    showSqrt2M();
    showText(" ");
    if (ctr == 0)
    {
      showText(ptrPlusMinus);
    }
    else
    {
      showText(ctr == 1 ? "+" : ptrMinus);
    }
    showText(" ");
    startSqrt();
    BigRationalDivideByInt(&RatDeprQuadratic, 2, &Rat1);
    BigIntChSign(&Rat1.numerator);
    showRatCoeffAndPowerVar(&Rat1, 0, currLetter);
    Rat2.numerator.sign = SIGN_NEGATIVE;
    showRatCoeffAndPowerVar(&Rat2, 1, 'm');
    if (ctr == 0)
    {
      BigIntChSign(&Rat4.numerator);
    }
    showNbrOverSqrt2M();
    endSqrt();
    if (ctr == 0)
    {
      showAdjustForCubic();
    }
    else
    {
      endShowX();
    }
  }
}

static void showFerrariResolventNoRationalRoot(int multiplicity)
{
  BigRational* coeffsToBackup[NBR_COEFFS_BACKUP] =
  {
    &RatDeprLinear,
    &RatDeprIndependent,
    &RatCubic,
    &RatQuadratic,
    &RatLinear,
    &RatIndependent,
    &RatDiscr
  };
  BigRational** pRat = coeffsToBackup;
  for (int i = 0; i < NBR_COEFFS_BACKUP; i++)
  {
    CopyBigInt(&coeffBackup[i].numerator, &(*pRat)->numerator);
    CopyBigInt(&coeffBackup[i].denominator, &(*pRat)->denominator);
    pRat++;
  }
  // Lodovico Ferrari's resolvent equation is:
  // 8x^3 + 8px^2 + (2p^2-8r)x - q^2 = 0
  // Initialize RatQuadratic, RatLinear and RatIndependent.
  CopyBigInt(&RatQuadratic.numerator, &RatDeprQuadratic.numerator);
  CopyBigInt(&RatQuadratic.denominator, &RatDeprQuadratic.denominator);
  BigRationalMultiply(&RatDeprQuadratic, &RatDeprQuadratic, &Rat2);
  BigRationalDivideByInt(&Rat2, 4, &Rat2);                       // (1/4)p^2
  BigRationalSubt(&Rat2, &RatDeprIndependent, &RatLinear);       // (1/4)p^2 - r
  BigRationalMultiply(&RatDeprLinear, &RatDeprLinear, &Rat3);    // q^2
  BigRationalDivideByInt(&Rat3, 8, &RatIndependent);             // (1/8)q^2
  BigIntChSign(&RatIndependent.numerator);                       // (-1/8)q^2
  solveCubic(1, true, 'm');
  pRat = coeffsToBackup;
  for (int i = 0; i < NBR_COEFFS_BACKUP; i++)
  {
    CopyBigInt(&(*pRat)->numerator, &coeffBackup[i].numerator);
    CopyBigInt(&(*pRat)->denominator, &coeffBackup[i].denominator);
    pRat++;
  }
  showText("<p>");
  // From $1q y $2q:
  formatString(&ptrOutput, LITERAL_STEPS_COMPLEX_SQROOT6, FerrariEq, FerrariEq + 1);
  showText("</p><p>");
  showRatCoeffAndPowerVar(NULL, -2, currLetter);
  BigRationalDivideByInt(&RatDeprQuadratic, 2, &Rat1);
  showRatCoeffAndPowerVar(&Rat1, 0, currLetter);
  showText(" + ");
  showM();
  showText(" = ");
  showText(ptrPlusMinus);
  showSquareRHS(false);
  generateEqNbr();
  showText("</p><p>");
  // Using the plus sign:
  showText(LITERAL_FERRARI_RAT_ROOT3);
  showText("</p>");
  showEndFerrariResolventNoRationalRoot(ptrMinus, multiplicity);
  showText("<p>");
  // Using the minus sign in $1q:
  formatString(&ptrOutput, LITERAL_FERRARI_RAT_ROOT4, eqNbr);
  showText("</p><p>");
  showEndFerrariResolventNoRationalRoot(" + ", multiplicity);
  showText("</p>");
}

static void FerrariResolventHasRationalRoot(int multiplicity)
{
  const int* ptrValues = factorInfoInteger[0].ptrPolyLifted;
  UncompressBigIntegerB(ptrValues, &RatS.numerator);
  ptrValues += numLimbs(ptrValues);
  ptrValues++;
  UncompressBigIntegerB(ptrValues, &RatS.denominator);   // RatS <- -root
  if (teach)
  {
    showFerrariResolventRationalRoot();
  }
  BigRationalDivideByInt(&RatS, -2, &RatS);              // RatS <- S^2 (as on Wikipedia article).
  ForceDenominatorPositive(&RatS);
  // According to Quartic function Wikipedia article the roots are:
  // x1,2 = -(b/4a) - S +/- sqrt(u + q/S)
  // x3,4 = -(b/4a) + S +/- sqrt(u - q/S)
  enum eSign oldLinearSign = RatDeprLinear.numerator.sign;
  if (teach)
  {
    // Using the plus sign:
    showText(LITERAL_FERRARI_RAT_ROOT3);
    showStepsOnRationalRootOfResolvent(ptrMinus, SIGN_POSITIVE);
  }
  showSolFerrariResolventRatRoot(0, oldLinearSign, multiplicity);
  showSolFerrariResolventRatRoot(2, oldLinearSign, multiplicity);
  if (teach)
  {
    showText("<p>");
    // Using the minus sign in $1q:
    formatString(&ptrOutput, LITERAL_FERRARI_RAT_ROOT4, eqPlusMinusAfterRoot);
    showText("</p>");
    showStepsOnRationalRootOfResolvent("+", SIGN_NEGATIVE);
  }
  showSolFerrariResolventRatRoot(1, oldLinearSign, multiplicity);
  showSolFerrariResolventRatRoot(3, oldLinearSign, multiplicity);
}

static void showDepressedQuartic(void)
{
  if (!BigIntIsOne(&Quartic))
  {
    // Dividing the equation by the quartic coefficient:
    formatString(&ptrOutput, "<p>$1s</p><p>", LITERAL_SHOW_DEPR_QUARTIC1);
    showRatCoeffAndPowerVar(NULL, -4, 'x');
    showRatCoeffAndPowerVar(&RatCubic, 3, 'x');
    showRatCoeffAndPowerVar(&RatQuadratic, 2, 'x');
    showRatCoeffAndPowerVar(&RatLinear, 1, 'x');
    showRatCoeffAndPowerVar(&RatIndependent, 0, 'x');
    showText(" = 0</p>");
  }
  if (!BigIntIsZero(&Cubic))
  {
    // To eliminate the cubic term, we will perform the following substitution:
    formatString(&ptrOutput, "<p>$1s</p>", LITERAL_SHOW_DEPR_QUARTIC2);
    showVariable(&ptrOutput, 'x');
    showText(" = ");
    showVariable(&ptrOutput, 'y');
    BigRationalDivideByInt(&RatCubic, -4, &Rat2);
    showPlusMinusRational(&Rat2);
    // The constant value in the substitution equals the fourth part of the cubic coefficient.
    formatString(&ptrOutput, "</p><p>$1s</p>", LITERAL_SHOW_DEPR_QUARTIC3);
    currLetter = 'y';
    // Show a(y-k)^4 + b(y-k)^3 + c(y-k)^2 + d(y-k) + e = 0.
    showText("<p>");
    startParen();
    showVariable(&ptrOutput, 'y');
    showPlusMinusRational(&Rat2);
    endParen();
    showPower(&ptrOutput, 4);
    // Cubic term is already non-zero.
    showCoeffBeforeParen(&RatCubic);
    startParen();
    showVariable(&ptrOutput, 'y');
    showPlusMinusRational(&Rat2);
    endParen();
    showPower(&ptrOutput, 3);
    if (!BigIntIsZero(&RatQuadratic.numerator))
    {  // Do not show quadratic term if it equals zero.
      showCoeffBeforeParen(&RatQuadratic);
      startParen();
      showVariable(&ptrOutput, 'y');
      showPlusMinusRational(&Rat2);
      endParen();
      showPower(&ptrOutput, 2);
    }
    if (!BigIntIsZero(&RatLinear.numerator))
    {  // Do not show linear term if it equals zero.
      showCoeffBeforeParen(&RatLinear);
      startParen();
      showVariable(&ptrOutput, 'y');
      showPlusMinusRational(&Rat2);
      endParen();
    }
    showPlusMinusRational(&RatIndependent);
    showText(" = 0</p><p>");
    // Expanding brackets:
    showText(LITERAL_SHOW_DEPR_QUARTIC4);
    showText("</p><p>");
    // Expand all terms.
    // Show (y-B/4)^4 = y^4 - B*y^3 + 3/8*B^2*y^2 - 1/16*B^3*y + 1/256*B^4
    showRatCoeffAndPowerVar(NULL, -4, 'y');
    CopyBigInt(&Rat4.numerator, &RatCubic.numerator);
    CopyBigInt(&Rat4.denominator, &RatCubic.denominator);
    BigIntChSign(&Rat4.numerator);
    showRatCoeffAndPowerVar(&Rat4, 3, 'y');
    BigRationalMultiplyByInt(&Rat4, -3, &Rat4);
    BigRationalDivideByInt(&Rat4, 8, &Rat4);
    BigRationalMultiply(&Rat4, &RatCubic, &Rat4);
    showRatCoeffAndPowerVar(&Rat4, 2, 'y');
    BigRationalDivideByInt(&Rat4, -6, &Rat4);
    BigRationalMultiply(&Rat4, &RatCubic, &Rat4);
    showRatCoeffAndPowerVar(&Rat4, 1, 'y');
    BigRationalDivideByInt(&Rat4, -16, &Rat4);
    BigRationalMultiply(&Rat4, &RatCubic, &Rat4);
    showPlusMinusRational(&Rat4);
    // Show B(y-B/4)^3 = B*y^3 - 3/4*B^2*y^2 + 3/16*B^3*y - 1/64*B^4
    showRatCoeffAndPowerVar(&RatCubic, 3, 'y');
    BigRationalMultiplyByInt(&RatCubic, -3, &Rat4);
    BigRationalDivideByInt(&Rat4, 4, &Rat4);
    BigRationalMultiply(&Rat4, &RatCubic, &Rat4);
    showRatCoeffAndPowerVar(&Rat4, 2, 'y');
    BigRationalDivideByInt(&Rat4, -4, &Rat4);
    BigRationalMultiply(&Rat4, &RatCubic, &Rat4);
    showRatCoeffAndPowerVar(&Rat4, 1, 'y');
    BigRationalDivideByInt(&Rat4, -12, &Rat4);
    BigRationalMultiply(&Rat4, &RatCubic, &Rat4);
    showPlusMinusRational(&Rat4);
    // Show C(y-B/4)^2 = C*y^2 - 1/2*C*B*y + 1/16*C*B^2
    showRatCoeffAndPowerVar(&RatQuadratic, 2, 'y');
    BigRationalDivideByInt(&RatQuadratic, -2, &Rat4);
    BigRationalMultiply(&Rat4, &RatCubic, &Rat4);
    showRatCoeffAndPowerVar(&Rat4, 1, 'y');
    BigRationalDivideByInt(&Rat4, -8, &Rat4);
    BigRationalMultiply(&Rat4, &RatCubic, &Rat4);
    showPlusMinusRational(&Rat4);
    // Show D(y-B/4)
    showRatCoeffAndPowerVar(&RatLinear, 1, 'y');
    BigRationalMultiply(&RatLinear, &RatCubic, &Rat4);
    BigRationalDivideByInt(&Rat4, -4, &Rat4);
    showPlusMinusRational(&Rat4);
    // Show E
    showPlusMinusRational(&RatIndependent);
    showText(" = 0</p><p>");
    // Simplifying:
    showText(LITERAL_SHOW_DEPR_QUARTIC5);
    showText("</p><p>");
    // Show y^4 + Py^2 + Qy + R = 0.
    showRatCoeffAndPowerVar(NULL, -4, 'y');
    showRatCoeffAndPowerVar(&RatDeprQuadratic, 2, 'y');
    showRatCoeffAndPowerVar(&RatDeprLinear, 1, 'y');
    showPlusMinusRational(&RatDeprIndependent);
    showText(" = 0</p>");
  }
}

static void showFerrariMethodDerivation(void)
{
  // We will use Lodovico Ferrari's method to solve the quartic equation
  // using an auxiliar cubic equation (named resolvent equation).
  // The right hand side of the following identity includes the first
  // two terms of the previous equation:
  formatString(&ptrOutput, "<p>$1s</p><p>$2s</p><p>", LITERAL_SHOW_FERRARI_METHOD1,
    LITERAL_SHOW_FERRARI_METHOD2);
  BigRationalDivideByInt(&RatDeprQuadratic, 2, &Rat1);
  showSquareLHS();
  showText(" = ");
  showRatCoeffAndPowerVar(NULL, -4, currLetter);
  showRatCoeffAndPowerVar(&RatDeprQuadratic, 2, currLetter);
  intToBigInteger(&Rat2.numerator, 2);
  intToBigInteger(&Rat2.denominator, 1);
  showRatCoeffAndPowerVar(&Rat2, 1, 'm');
  showText(ptrTimes);
  showRatCoeffAndPowerVar(NULL, -2, currLetter);
  showText(" + ");
  showRatCoeffAndPowerVar(NULL, -2, 'm');
  showRatCoeffAndPowerVar(&RatDeprQuadratic, 1, 'm');
  BigRationalMultiply(&RatDeprQuadratic, &RatDeprQuadratic, &Rat3);
  BigRationalDivideByInt(&Rat3, 4, &Rat3);
  showRatCoeffAndPowerVar(&Rat3, 0, 'm');
  showText("</p><p>");
  // Since the equation equals zero,
  // we can subtract it from the right hand side of the identity.
  showText(LITERAL_SHOW_FERRARI_METHOD3);
  showText("</p><p>");
  showSquareLHS();
  showText(" = 2");
  showText(ptrTimes);
  showM();
  showText(ptrTimes);
  showRatCoeffAndPowerVar(NULL, -2, currLetter);
  BigRationalMultiplyByInt(&RatDeprLinear, -1, &Rat4);
  showRatCoeffAndPowerVar(&Rat4, 1, currLetter);
  showText(" + ");
  showRatCoeffAndPowerVar(NULL, -2, 'm');
  showRatCoeffAndPowerVar(&RatDeprQuadratic, 1, 'm');
  BigRationalSubt(&Rat3, &RatDeprIndependent, &Rat3);
  showRatCoeffAndPowerVar(&Rat3, 0, 'm');
  showText("</p><p>");
  // We select the value of $1v such that the right hand side be a perfect square.
  formatString(&ptrOutput, LITERAL_SHOW_FERRARI_METHOD4, 'm');
  showText("</p><p>");
  // From the identity
  showText(LITERAL_SHOW_FERRARI_METHOD5);
  showText("</p><p>");
  showSquareRHS(true);
  showText(" = 2");
  showText(ptrTimes);
  showVariable(&ptrOutput, 'm');
  showText(ptrTimes);
  showRatCoeffAndPowerVar(NULL, -2, currLetter);
  BigRationalMultiplyByInt(&RatDeprLinear, -1, &Rat4);
  showRatCoeffAndPowerVar(&Rat4, 1, currLetter);
  BigRationalMultiply(&RatDeprLinear, &RatDeprLinear, &Rat4);
  BigRationalDivideByInt(&Rat4, 8, &Rat4);
  showText(" + ");
  showRationalNoParenOverGeneric(&Rat4, showM);
  showText("</p><p>");
  // we get:
  showText(LITERAL_BIQUADR_EQ3);
  showText("</p><p>");
  showSquareLHS();
  showText(" = ");
  showSquareRHS(true);
  showText(" + ");
  showRatCoeffAndPowerVar(NULL, -2, 'm');
  showRatCoeffAndPowerVar(&RatDeprQuadratic, 1, 'm');
  showRatCoeffAndPowerVar(&Rat3, 0, 'm');
  BigRationalMultiply(&RatDeprLinear, &RatDeprLinear, &Rat4);
  BigRationalDivideByInt(&Rat4, 8, &Rat4);
  showText(" ");
  showText(ptrMinus);
  showText(" ");
  showRationalNoParenOverGeneric(&Rat4, showM);
  generateEqNbr();    // Equation 1.
  FerrariEq = eqNbr;
  showText("</p><p>");
  // We can make the right hand side a perfect square by setting to zero what is outside
  // the parentheses. In this way we get a cubic equation:
  showText(LITERAL_SHOW_FERRARI_METHOD6);
  showText("</p><p>");
  BigRationalMultiplyByInt(&Rat4, -1, &Rat4);
  showRatCoeffAndPowerVar(NULL, -3, 'm');
  showRatCoeffAndPowerVar(&RatDeprQuadratic, 2, 'm');
  showRatCoeffAndPowerVar(&Rat3, 1, 'm');
  showRatCoeffAndPowerVar(&Rat4, 0, 'm');
  showText(" = 0");
  generateEqNbr();
  showText("</p>");
}

static void FerrariDiscriminantIsPositive(int multiplicity)
{
  bool isImaginary;
  enum eSign sign1;
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
  if (BigIntIsOne(&Rat2.numerator))
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
  for (int ctr = 0; ctr < 4; ctr++)
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
    endShowX();
  }
}

static void FerrariDiscriminantIsNegative(int multiplicity)
{
  bool isImaginary;
  enum eSign sign2;
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
  for (int ctr = 0; ctr < 4; ctr++)
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
    endShowX();
  }
}

void QuarticEquation(const int* polynomial, int multiplicity)
{
  int* ptrValues;
  const int* ptrPolynomial = polynomial;
  currLetter = 'x';
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
  if (teach)
  {
    showDepressedQuartic();
  }
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
  if (teach)
  {
    showFerrariMethodDerivation();
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
  common.poly.values[0] = 3;
  ptrValues = &common.poly.values[1];
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
  else if (teach)
  {
    showFerrariResolventNoRationalRoot(multiplicity);
  }
  else
  {
    if (RatDiscr.numerator.sign == SIGN_POSITIVE)
    {
      FerrariDiscriminantIsPositive(multiplicity);
    }
    else
    {
      FerrariDiscriminantIsNegative(multiplicity);
    }
  }
}

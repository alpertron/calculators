//
// This file is part of Alpertron Calculators.
//
// Copyright 2019-2021 Dario Alejandro Alpern
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

#define NBR_COEFF 6

BigInteger Quintic;
BigInteger Quartic;
BigInteger Cubic;
BigInteger Quadratic;
BigInteger Linear;
BigInteger Independent;
BigInteger discr;
BigInteger commonDenom;
BigInteger tmp0;
BigInteger tmp1;
BigInteger tmp2;
BigInteger tmp3;
BigInteger tmp4;
BigInteger tmp5;
BigInteger tmp6;
BigInteger tmp7;
BigRational RatQuartic;
BigRational RatCubic;
BigRational RatQuadratic;
BigRational RatLinear;
BigRational RatIndependent;
BigRational RatDeprCubic;
BigRational RatDeprQuadratic;
BigRational RatDeprLinear;
BigRational RatDeprIndependent;
BigRational RatDiscr;
BigRational RatDelta0;
BigRational RatDelta1;
BigRational RatD;
BigRational Rat1;
BigRational Rat2;
BigRational Rat3;
BigRational Rat4;
BigRational Rat5;
BigRational RatS;
int indexRoot;
const char *ptrMinus;
const char *ptrTimes;
const char *ptrSin;
const char *ptrCos;
const char* ptrACos;
const char *ptrPi;
const char *ptrI;
static int totients[(2 * MAX_DEGREE) + 1];
extern char* ptrOutput;

enum
{
  SHOW_REAL = 0,
  SHOW_IMAG,
};

struct cosine12
{
  char* radicands;
  int denominator;
};

static struct cosine12 aCosine12[] =
{   // S = sqrt next digit, [ = start sqrt, ] = end sqrt
  {"(S6*[5+S5]+S5-1)", 8},          // cos 12°
  {"(S6*[5-S5]+S5+1)", 8},          // cos 24°
  {"(S5+1)", 4},                    // cos 36°
  {"(S6*[5+S5]-S5+1)", 8},          // cos 48°
  {"1", 2},                         // cos 60°
  {"(S5-1)", 4},                    // cos 72°
  {"(S6*[5-S5]-S5-1)", 8},          // cos 84°
};

static struct cosine12 a2Plus2Cosine12[] =
{   // S = sqrt next digits, [ = start sqrt, ] = end sqrt
  {"7+S6*[5+S5]+S5", 4},          // 2 + 2*cos 12°
  {"9+S6*[5-S5]+S5", 4},          // 2 + 2*cos 24°
  {"10+2*S5", 4},                 // 2 + 2*cos 36°
  {"9+S6*[5+S5]-S5", 4},          // 2 + 2*cos 48°
  {"3", 2},                       // 2 + 2*cos 60°
  {"6-2*S5", 4},                  // 2 + 2*cos 72°
  {"7+S6*[5-S5]-S5", 4},          // 2 + 2*cos 84°
};

static struct cosine12 a2Minus2Cosine12[] =
{   // S = sqrt next digits, [ = start sqrt, ] = end sqrt
  {"9-S6*[5+S5]-S5", 4},          // 2 - 2*cos 12°
  {"7-S6*[5-S5]-S5", 4},          // 2 - 2*cos 24°
  {"6-2*S5", 4},                  // 2 - 2*cos 36°
  {"7-S6*[5+S5]+S5", 4},          // 2 - 2*cos 48°
  {"1", 2},                       // 2 - 2*cos 60°
  {"10-2*S5", 4},                 // 2 - 2*cos 72°
  {"9-S6*[5-S5]+S5", 4},          // 2 - 2*cos 84°
};

#if 0
-cos(Pi/17) + (1/16)*(1-sqrt(17)+sqrt(34-2*sqrt(17))+2*sqrt(17+3*sqrt(17)+sqrt(170+38*sqrt(17))))
-cos(2*Pi/17) + (1/16)*(-1+sqrt(17)+sqrt(34-2*sqrt(17))+2*sqrt(17+3*sqrt(17)-sqrt(170+38*sqrt(17))))
-cos(3*Pi/17) + (1/16)*(1+sqrt(17)+sqrt(34+2*sqrt(17))+2*sqrt(17-3*sqrt(17)-sqrt(170-38*sqrt(17))))
-cos(4*Pi/17) + (1/16)*(-1+sqrt(17)-sqrt(34-2*sqrt(17))+2*sqrt(17+3*sqrt(17)+sqrt(170+38*sqrt(17))))
-cos(5*Pi/17) + (1/16)*(1+sqrt(17)+sqrt(34+2*sqrt(17))-2*sqrt(17-3*sqrt(17)-sqrt(170-38*sqrt(17))))
-cos(6*Pi/17) + (1/16)*(-1-sqrt(17)+sqrt(34+2*sqrt(17))+2*sqrt(17-3*sqrt(17)+sqrt(170-38*sqrt(17))))
-cos(7*Pi/17) + (1/16)*(1+sqrt(17)-sqrt(34+2*sqrt(17))+2*sqrt(17-3*sqrt(17)+sqrt(170-38*sqrt(17))))
-cos(8*Pi/17) + (1/16)*(-1+sqrt(17)+sqrt(34-2*sqrt(17))-2*sqrt(17+3*sqrt(17)-sqrt(170+38*sqrt(17))))
#endif


static struct cosine12 aCosine17[] =
{   // S = sqrt next digits, [ = start sqrt, ] = end sqrt
  {"(1-S17+[34-2*S17]+2*[17+3*S17+[170+38*S17]])", 16},   // cos(pi/17)
  {"(-1+S17+[34-2*S17]+2*[17+3*S17-[170+38*S17]])", 16},  // cos(2*pi/17)
  {"(1+S17+[34+2*S17]+2*[17-3*S17-[170-38*S17]])", 16},   // cos(3*pi/17)
  {"(-1+S17-[34-2*S17]+2*[17+3*S17+[170+38*S17]])", 16},  // cos(4*pi/17)
  {"(1+S17+[34+2*S17]-2*[17-3*S17-[170-38*S17]])", 16},   // cos(5*pi/17)
  {"(-1-S17+[34+2*S17]+2*[17-3*S17+[170-38*S17]])", 16},  // cos(6*pi/17)
  {"(1+S17-[34+2*S17]+2*[17-3*S17+[170-38*S17]])", 16},   // cos(7*pi/17)
  {"(-1+S17+[34-2*S17]-2*[17+3*S17-[170+38*S17]])", 16},  // cos(8*pi/17)
};

static struct cosine12 a2PM2Cosine17[] =
{   // S = sqrt next digits, [ = start sqrt, ] = end sqrt
  {"34-2*S17+2*[34-2*S17]+4*[17+3*S17+[170+38*S17]]", 1},  // 2+2*cos(pi/17)
  {"34-2*S17-2*[34-2*S17]-4*[17+3*S17-[170+38*S17]]", 1},  // 2-2*cos(2*pi/17)
  {"34+2*S17+2*[34+2*S17]+4*[17-3*S17-[170-38*S17]]", 1},  // 2+2*cos(3*pi/17)
  {"34-2*S17+2*[34-2*S17]-4*[17+3*S17+[170+38*S17]]", 1},  // 2-2*cos(4*pi/17)
  {"34+2*S17+2*[34+2*S17]-4*[17-3*S17-[170-38*S17]]", 1},  // 2+2*cos(5*pi/17)
  {"34+2*S17-2*[34+2*S17]-4*[17-3*S17+[170-38*S17]]", 1},  // 2-2*cos(6*pi/17)
  {"34+2*S17-2*[34+2*S17]+4*[17-3*S17+[170-38*S17]]", 1},  // 2+2*cos(7*pi/17)
  {"34-2*S17-2*[34-2*S17]+4*[17+3*S17-[170+38*S17]]", 1},  // 2-2*cos(8*pi/17)
};

void startLine(void)
{
  showText("<li>");
  if (pretty == TEX)
  {
    showText("\\bullet\\,\\,");
  }
}

void endLine(void)
{
  if (pretty == TEX)
  {
    showText("\\\\");
  }
  showText("</li>");
}

static void showXindex(int index)
{
  if (pretty == PRETTY_PRINT)
  {
    showText("<var>x</var><sub>");
    int2dec(&ptrOutput, index);
    showText("</sub>");
  }
  else if (pretty == TEX)
  {
    showText("x_{");
    int2dec(&ptrOutput, index);
    showText("}");
  }
  else
  {
    showText("x");
    int2dec(&ptrOutput, index);
  }
}

void showX(int multiplicity)
{
  startLine();
  if (multiplicity > 2)
  {
    showXindex(indexRoot);
    showText(lang ? " a " : " to ");
    indexRoot += multiplicity;
    showXindex(indexRoot - 1);
    showText(" = ");
  }
  else
  {
    for (int ctr = 0; ctr < multiplicity; ctr++)
    {
      showXindex(indexRoot);
      indexRoot++;
      showText(" = ");
    }
  }
}

void startSqrt(void)
{
  if (pretty == PRETTY_PRINT)
  {
    showText("<span class=\"root\"><span class=\"radicand2\">");
  }
  else if (pretty == TEX)
  {
    showText("\\sqrt{");
  }
  else
  {
    showText("(");
  }
}

void endSqrt(void)
{
  if (pretty == PRETTY_PRINT)
  {
    showText("</span></span>");
  }
  else if (pretty == TEX)
  {
    showText("}");
  }
  else
  {
    showText(")^(1/2)");
  }
}

void startCbrt(void)
{
  if (pretty == PRETTY_PRINT)
  {
    showText("<span class=\"root\"><span class=\"befrad\" aria-hidden=\"true\">3</span><span class=\"radicand3\">");
  }
  else if (pretty == TEX)
  {
    showText("\\sqrt[3]{");
  }
  else
  {
    showText("(");
  }
}

void endCbrt(void)
{
  if (pretty == PRETTY_PRINT)
  {
    showText("</span></span>");
  }
  else if (pretty == TEX)
  {
    showText("}");
  }
  else
  {
    showText(")^(1/3)");
  }
}

void startParen(void)
{
  if (pretty == PRETTY_PRINT)
  {
    if (lang)
    {
      showText("<span class=\"hide\"> paréntesis izquierdo </span>");
    }
    else
    {
      showText("<span class=\"hide\"> left parenthesis </span>");
    }
    showText("<span class=\"paren\">");
  }
  else if (pretty == TEX)
  {
    showText("\\left(");
  }
  else
  {
    *ptrOutput ='(';
    ptrOutput++;
  }
}

void endParen(void)
{
  if (pretty == PRETTY_PRINT)
  {
    showText("</span>");
    if (lang)
    {
      showText("<span class=\"hide\"> paréntesis derecho </span>");
    }
    else
    {
      showText("<span class=\"hide\"> right parenthesis </span>");
    }
  }
  else if (pretty == TEX)
  {
    showText("\\right)");
  }
  else
  {
    *ptrOutput =')';
    ptrOutput++;
  }
}

void showPlusSignOn(bool condPlus, int type)
{
  if ((type & TYPE_PM_SPACE_BEFORE) != 0)
  {
    *ptrOutput = ' ';
    ptrOutput++;
  }
  if (condPlus)
  {
    *ptrOutput = '+';
    ptrOutput++;
  }
  else
  {
    showText(ptrMinus);
  }
  if ((type & TYPE_PM_SPACE_AFTER) != 0)
  {
    *ptrOutput = ' ';
    ptrOutput++;
  }
}
// Compute x = -c_0 / c_1
static void LinearEquation(const int *polynomial, int multiplicity)
{
  const int* ptrPolynomial = polynomial;
  showX(multiplicity);
  UncompressBigIntegerB(ptrPolynomial, &Independent);
  ptrPolynomial += numLimbs(ptrPolynomial);
  ptrPolynomial++;
  UncompressBigIntegerB(ptrPolynomial, &Linear);
  CopyBigInt(&Rat1.numerator, &Independent);
  CopyBigInt(&Rat1.denominator, &Linear);
  BigIntChSign(&Rat1.numerator);
  ForceDenominatorPositive(&Rat1);
  showRationalNoParen(&Rat1);
  endLine();
}

// From variables named Quadratic, Linear and Independent, compute the
// discriminant and find the rational roots if there are. Otherwise find
// the real or imaginary roots.
// The output is zero if the roots are rational, 1 otherwise.
static int ProcessQuadraticEquation(enum eSign* pSignDescr)
{
  // Compute discriminant (delta = linear^2 - 4*quadratic*independent).
  (void)BigIntMultiply(&Linear, &Linear, &tmp1);
  (void)BigIntMultiply(&Quadratic, &Independent, &tmp2);
  multint(&tmp2, &tmp2, 4);
  BigIntSubt(&tmp1, &tmp2, &discr);
  if (discr.sign == SIGN_POSITIVE)
  {
    squareRoot(discr.limbs, tmp4.limbs, discr.nbrLimbs, &tmp4.nbrLimbs);
    tmp4.sign = SIGN_POSITIVE;
    (void)BigIntMultiply(&tmp4, &tmp4, &tmp5);
    BigIntSubt(&tmp5, &discr, &tmp5);
    if (BigIntIsZero(&tmp5))
    {           // Discriminant is perfect square. Roots are rational numbers.
      return 0;
    }
  }
  *pSignDescr = discr.sign;
  discr.sign = SIGN_POSITIVE;
  // Compute Rat1 as -linear/(2*quadratic), Rat2 as abs(1/(2*quadratic))
  // and Rat3 as abs(delta).
  CopyBigInt(&Rat1.numerator, &Linear);
  CopyBigInt(&Rat1.denominator, &Quadratic);
  intToBigInteger(&Rat2.numerator, 1);
  CopyBigInt(&Rat2.denominator, &Quadratic);
  CopyBigInt(&Rat3.numerator, &discr);
  intToBigInteger(&Rat3.denominator, 1);
  BigRationalDivideByInt(&Rat1, -2, &Rat1);
  BigRationalDivideByInt(&Rat2, 2, &Rat2);
  ForceDenominatorPositive(&Rat1);
  ForceDenominatorPositive(&Rat2);
  Rat2.numerator.sign = SIGN_POSITIVE;
  return 1;
}

// Compute delta = c_1^2 - 4c_0 * c_2
// If delta > 0 -> x = (-c_1 +/- sqrt(delta))/(2*c_2)
// If delta < 0 -> x = (-c_1 +/- i*sqrt(-delta))/(2*c_2)
// Delta cannot be zero because the roots are different.
static void QuadraticEquation(const int* polynomial, int multiplicity)
{
  const int* ptrPolynomial = polynomial;
  int ctr;
  enum eSign signDiscr;
  UncompressBigIntegerB(ptrPolynomial, &Independent);
  ptrPolynomial += numLimbs(ptrPolynomial);
  ptrPolynomial++;
  UncompressBigIntegerB(ptrPolynomial, &Linear);
  ptrPolynomial += numLimbs(ptrPolynomial);
  ptrPolynomial++;
  UncompressBigIntegerB(ptrPolynomial, &Quadratic);
  if (ProcessQuadraticEquation(&signDiscr) == 0)
  {           // Discriminant is perfect square. Roots are rational numbers.
    for (ctr = 0; ctr < 2; ctr++)
    {         // Compute (-linear +/- sqrt(discr))/(2*quadratic)
      showX(multiplicity);
      if (ctr == 0)
      {
        BigIntAdd(&Linear, &tmp4, &tmp1);
      }
      else
      {
        BigIntSubt(&Linear, &tmp4, &tmp1);
      }
      CopyBigInt(&Rat1.numerator, &tmp1);        
      CopyBigInt(&Rat1.denominator, &Quadratic);
      BigRationalDivideByInt(&Rat1, -2, &Rat1);
      ForceDenominatorPositive(&Rat1);
      showRationalNoParen(&Rat1);        
      endLine();
    }
    return;
  }
  for (ctr = 0; ctr < 2; ctr++)
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
      *ptrOutput = 'i';  
      ptrOutput++;
    }
    endLine();
  }
}

static void CbrtIndep(void)
{
  if ((pretty != PARI_GP) || !BigIntIsOne(&RatDeprIndependent.denominator))
  {
    startCbrt();
    showRationalNoParen(&RatDeprIndependent);
    endCbrt();
  }
  else
  {
    showRationalNoParen(&RatDeprIndependent);
    showText("^(1/3)");
  }
}

void showRatConstants(const char* numerator, const char* denominator)
{
  if (pretty == PRETTY_PRINT)
  {
    showText("<span class=\"fraction\"><span class=\"numerator\">");
    showText(numerator);
    showText("</span><span class=\"denominator\">");
    showText(denominator);
    showText("</span></span>");
  }
  else if (pretty == TEX)
  {
    showText("\\frac{");
    showText(numerator);
    showText("}{");
    showText(denominator);
    showText("}");
  }
  else
  {
    showText("(");
    showText(numerator);
    showText("/");
    showText(denominator);
    showText(")");
  }
}

static void CubicEquation(const int* polynomial, int multiplicity)
{
  int ctr;
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
  // Get rational coefficients of monic equation.
  CopyBigInt(&RatQuadratic.numerator, &Quadratic);
  CopyBigInt(&RatQuadratic.denominator, &Cubic);
  CopyBigInt(&RatLinear.numerator, &Linear);
  CopyBigInt(&RatLinear.denominator, &Cubic);
  CopyBigInt(&RatIndependent.numerator, &Independent);
  CopyBigInt(&RatIndependent.denominator, &Cubic);
  // Compute coefficients of depressed equation x^3 + px + q
  // where: p = -b^2/3 + c, q = 2b^3/27 - bc/3 + d.
  BigRationalDivideByInt(&RatQuadratic, 3, &Rat1);
  BigRationalMultiply(&Rat1, &RatQuadratic, &RatDeprLinear);   // b^2/3
  BigRationalMultiply(&Rat1, &RatLinear, &RatDeprIndependent); // bc/3
  BigRationalMultiply(&RatDeprLinear, &RatQuadratic, &Rat1);   // b^3/3
  BigRationalMultiplyByInt(&Rat1, 2, &Rat1);                   // 2b^3/3
  BigRationalDivideByInt(&Rat1, 9, &Rat1);                     // 2b^3/27
  BigRationalSubt(&RatLinear, &RatDeprLinear, &RatDeprLinear); // p
  BigRationalSubt(&RatIndependent, &RatDeprIndependent, &RatDeprIndependent); // d - bc/3
  BigRationalAdd(&RatDeprIndependent, &Rat1, &RatDeprIndependent); // q
  // Compute discriminant (delta)/(-27) = q^2 + 4p^3/27.
  BigRationalMultiply(&RatDeprIndependent, &RatDeprIndependent, &RatDiscr);  // q^2
  BigRationalMultiply(&RatDeprLinear, &RatDeprLinear, &Rat1);  // p^2
  BigRationalMultiply(&RatDeprLinear, &Rat1, &Rat1);           // p^3
  BigRationalMultiplyByInt(&Rat1, 4, &Rat1);                   // 4p^3
  BigRationalDivideByInt(&Rat1, 27, &Rat1);                    // 4p^3/27
  BigRationalAdd(&RatDiscr, &Rat1, &RatDiscr);                 // q^2 + 4p^3/27
  BigRationalDivideByInt(&RatQuadratic, -3, &RatQuadratic);    // -b/3
  ForceDenominatorPositive(&RatQuadratic);
  if (RatDiscr.numerator.sign == RatDiscr.denominator.sign)
  {   // Discriminant is negative. Use Cardan's formula.
    if (BigIntIsZero(&RatDeprLinear.numerator))
    {
      // The first root x1 is q ^ (1/3)
      // The second root x2 is (-1/2) * x1 + i * (3 ^ (1/2) * x1) / 2
      // The third root x3 is (-1/2) * x1 - i * (3 ^ (1/2) * x1) / 2
      bool isCbrtNegative;
      BigRationalNegate(&RatDeprIndependent, &RatDeprIndependent);
      ForceDenominatorPositive(&RatDeprIndependent);
      isCbrtNegative = (RatDeprIndependent.numerator.sign == SIGN_NEGATIVE);
      RatDeprIndependent.numerator.sign = SIGN_POSITIVE;
      for (ctr = 0; ctr < 3; ctr++)
      {
        showX(multiplicity);
        if (BigIntIsZero(&Quadratic))
        {
          if ((ctr == 0) && isCbrtNegative)
          {
            showText(ptrMinus);
          }
        }
        else
        {
          showRationalNoParen(&RatQuadratic);
          if ((ctr == 0) && isCbrtNegative)
          {
            *ptrOutput = ' ';
            ptrOutput++;
            showText(ptrMinus);
            *ptrOutput = ' ';
            ptrOutput++;
          }
          else
          {
            showText(" + ");
          }
        }
        if (ctr == 0)
        {
          CbrtIndep();     // q^(1/3)
        }
        else
        {
          if (pretty != PARI_GP)
          {
            if (!isCbrtNegative)
            {
              showText(ptrMinus);
            }
            *ptrOutput = ' ';
            ptrOutput++;
            showRatConstants("1", "2");
            *ptrOutput = ' ';
            ptrOutput++;
            showText(ptrTimes);
          }
          else
          {
            if (isCbrtNegative)
            {
              showText("(1/2) *");
            }
            else
            {
              showText("(-1/2) *");
            }
          }
          CbrtIndep();     // q^(1/3)
          showPlusSignOn(ctr == 1, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
          if (pretty != PARI_GP)
          {                // Show i/2.
            showRatConstants("i", "2");
          }
          else
          {
            showText("I");
          }
          showText(ptrTimes);
          if (pretty != PARI_GP)
          {     // square root of 3.
            startSqrt();
            *ptrOutput = '3';
            ptrOutput++;
            endSqrt();
          }
          else
          {     // square root of 3.
            showText("3^(1/2)");
          }
          *ptrOutput = ' ';
          ptrOutput++;
          showText(ptrTimes);
          CbrtIndep();     // q^(1/3)
          if (pretty == PARI_GP)
          {
            showText("/2");
          }
        }
        endLine();
      }
    }
    else
    {
      BigRationalDivideByInt(&RatDeprIndependent, -2, &Rat3);
      intToBigInteger(&Rat1.numerator, 1);
      intToBigInteger(&Rat1.denominator, 1);
      BigRationalDivideByInt(&RatDiscr, 4, &Rat2);
      MultiplyRationalBySqrtRational(&Rat1, &Rat2);
      for (ctr = 0; ctr < 2; ctr++)
      {
        startLine();
        if (pretty == PRETTY_PRINT)
        {
          showText("<var>");
        }
        *ptrOutput = ((ctr == 0)? 'r' : 's');
        ptrOutput++;
        if (pretty == PRETTY_PRINT)
        {
          showText("</var>");
        }
        *ptrOutput = ' ';
        ptrOutput++;
        *ptrOutput = '=';
        ptrOutput++;
        *ptrOutput = ' ';
        ptrOutput++;
        // Pari-GP requires the argument of the cube root to be positive.
        if ((Rat2.numerator.nbrLimbs == 1) && (Rat2.numerator.limbs[0].x == 1) &&
          (Rat2.denominator.nbrLimbs == 1) && (Rat2.denominator.limbs[0].x == 1))
        {             // Cube root of rational number.
          if (ctr == 0)
          {
            BigRationalAdd(&Rat3, &Rat1, &Rat4);
          }
          else
          {
            BigRationalSubt(&Rat3, &Rat1, &Rat4);
          }
          ForceDenominatorPositive(&Rat4);
          if (Rat4.numerator.sign == SIGN_NEGATIVE)
          {
            showText(ptrMinus);
            Rat4.numerator.sign = SIGN_POSITIVE;
          }
          startCbrt();
          showRationalNoParen(&Rat4);
          endCbrt();
        }
        else
        {           // Cube root of rational plus/minus a square root.
          bool isCbrtNegative;
          ForceDenominatorPositive(&Rat3);
          double dLogRat = logBigNbr(&Rat3.numerator) - logBigNbr(&Rat3.denominator);
          double dLogSqrt = 
            (logBigNbr(&Rat1.numerator) - logBigNbr(&Rat1.denominator)) +
            (logBigNbr(&Rat2.numerator) - logBigNbr(&Rat2.denominator)) / 2.0;
          if (ctr == 0)
          {          // Square root positive
            isCbrtNegative = ((Rat3.numerator.sign == SIGN_NEGATIVE) &&
              (dLogRat > dLogSqrt));
          }
          else
          {
            isCbrtNegative = ((Rat3.numerator.sign == SIGN_NEGATIVE) ||
              (dLogRat < dLogSqrt));
          }
          if (isCbrtNegative)
          {
            showText(ptrMinus);
            BigIntChSign(&Rat3.numerator);
          }
          startCbrt();
          showRationalNoParen(&Rat3);
          if (isCbrtNegative)
          {
            BigIntChSign(&Rat3.numerator);
          }
          showPlusSignOn(ctr == (isCbrtNegative? 1: 0),
            TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
          ShowRationalAndSqrParts(&Rat1, &Rat2, 2, ptrTimes);
          endCbrt();
        }
        endLine();
      }
      for (ctr = 0; ctr < 3; ctr++)
      {
        showX(multiplicity);
        if (!BigIntIsZero(&Quadratic))
        {
          showRationalNoParen(&RatQuadratic);
          showPlusSignOn(ctr == 0, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
        }
        if (ctr == 0)
        {
          showText("<var>r</var> + <var>s</var>");
        }
        else
        {
          if (BigIntIsZero(&RatQuadratic.numerator))
          {
            showText(ptrMinus);
          }
          if (pretty == PRETTY_PRINT)
          {
            showRatConstants("<var>r</var> + <var>s</var>", "2");
          }
          else if (pretty == TEX)
          {
            showRatConstants("r + s", "2");
          }
          else
          {
            showText("(<var>r</var> + <var>s</var>) / 2");
          }
          showPlusSignOn(ctr == 1, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
          if (pretty != PARI_GP)
          {
            showText("i");
            *ptrOutput = ' ';
            ptrOutput++;
            showText(ptrTimes);
            *ptrOutput = ' ';
            ptrOutput++;
            showRatConstants((pretty == TEX)? "r - s" :
              "<var>r</var> &minus; <var>s</var>", "2");
            showText(ptrTimes);
            startSqrt();
            *ptrOutput = '3';
            ptrOutput++;
            endSqrt();
          }
          else
          {
            showText("(I/2) * (<var>r</var> - <var>s</var>) * 3^(1/2)");
          }
        }
        endLine();
      }
    }
  }
  else
  {   // Discriminant is positive. Use Viete's formula.
    enum eSign signRat1;
    startLine();
    if (pretty == PRETTY_PRINT)
    {
      showText("<var>t</var> = ");
    }
    else
    {
      showText("t = ");
    }
    if (pretty != PARI_GP)
    {
      showRatConstants("1", "3");
      showText(ptrTimes);
    }
    else
    {
      showText("(1/3) * ");
    }
    showText(ptrACos);
    startParen();
    BigRationalDivide(&RatDeprIndependent, &RatDeprLinear, &Rat1); // q/p
    BigRationalMultiplyByInt(&Rat1, 3, &Rat1);                     // 3q/p
    BigRationalDivideByInt(&Rat1, 2, &Rat1);                       // 3q/(2p)
    ForceDenominatorPositive(&Rat1);
    CopyBigInt(&Rat2.numerator, &RatDeprLinear.denominator);
    CopyBigInt(&Rat2.denominator, &RatDeprLinear.numerator);       // 1/p
    BigRationalMultiplyByInt(&Rat2, -3, &Rat2);                    // -3/p
    ForceDenominatorPositive(&Rat2);
    MultiplyRationalBySqrtRational(&Rat1, &Rat2);
    ShowRationalAndSqrParts(&Rat1, &Rat2, 2, ptrTimes);
    endParen();
    endLine();
    intToBigInteger(&Rat1.numerator, 2);
    intToBigInteger(&Rat1.denominator, 1);               // 2
    BigRationalDivideByInt(&RatDeprLinear, -3, &Rat2);   // -p/3
    MultiplyRationalBySqrtRational(&Rat1, &Rat2);
    signRat1 = Rat1.numerator.sign;
    Rat1.numerator.sign = SIGN_POSITIVE;
    for (ctr = 0; ctr <= 4; ctr += 2)
    {
      showX(multiplicity);
      if (!BigIntIsZero(&Quadratic))
      {
        showRationalNoParen(&RatQuadratic);
        if (signRat1 == SIGN_POSITIVE)
        {
          showText(" + ");
        }
      }
      if (signRat1 == SIGN_NEGATIVE)
      {
        *ptrOutput = ' ';
        ptrOutput++;
        showText(ptrMinus);
        *ptrOutput = ' ';
        ptrOutput++;
      }
      ShowRationalAndSqrParts(&Rat1, &Rat2, 2, ptrTimes);
      showText(ptrTimes);
      *ptrOutput = ' ';
      ptrOutput++;
      showText(ptrCos);
      if (ctr == 0)
      {
        showText("(<var>t</var>)");
      }
      else
      {  
        startParen();    
        showText("<var>t</var> + ");
        if (ctr == 2)
        {
          if (pretty != PARI_GP)
          {
            char numer[200];
            char* ptrNumer = numer;
            *ptrNumer = '2';
            ptrNumer++;
            copyStr(&ptrNumer, ptrTimes);
            *ptrNumer = ' ';
            ptrNumer++;
            copyStr(&ptrNumer, ptrPi);
            showRatConstants(numer, "3");
          }
          else
          {
            showText("2 * Pi / 3");
          }
        }
        if (ctr == 4)
        {
          if (pretty != PARI_GP)
          {
            char numer[200];
            char* ptrNumer = numer;
            *ptrNumer = '4';
            ptrNumer++;
            copyStr(&ptrNumer, ptrTimes);
            *ptrNumer = ' ';
            ptrNumer++;
            copyStr(&ptrNumer, ptrPi);
            showRatConstants(numer, "3");
          }
          else
          {
            showText("4 * Pi / 3");
          }
        }
        endParen();
        if (pretty == TEX)
        {
          *ptrOutput = '}';
          ptrOutput++;
        }
      }
    }
    endLine();
  }
}

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
  BigRationalDivideByInt(&Rat4, 4*16, &Rat4);
  // Divide rational part by 2*4.
  BigRationalDivideByInt(&Rat3, 2*4, &Rat3);
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
  BigRationalMultiplyByInt(&Rat4, 4*16, &Rat4);
  // Restore rational part.
  BigRationalMultiplyByInt(&Rat3, 2*4, &Rat3);
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

static void QuarticEquation(const int* polynomial, int multiplicity)
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
      showText((pretty == PRETTY_PRINT)?" <var>S</var> ": " S ");
      showPlusSignOn((ctr == 0) || (ctr == 2), TYPE_PM_SPACE_AFTER);
      if (pretty != PARI_GP)
      {
        showRatConstants((isImaginary? "i": "1"), "2");
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
    for (ctr = 0; ctr < 4; ctr++)
    {
      showX(multiplicity);
      showFirstTermQuarticEq(ctr);
      showText(" <var>S</var> ");
      showPlusSignOn((ctr == 0) || (ctr == 2), TYPE_PM_SPACE_AFTER);
      if (pretty != PARI_GP)
      {
        showRatConstants((isImaginary? "i":"1"), "2");
      }
      else
      {
        showText(isImaginary ? "(I/2)": "(1/2)");
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
      showText((pretty == PRETTY_PRINT)?"&sup2;": "^2");
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

static int gcd(int first, int second)
{
  int a = first;
  int b = second;
  while (b != 0)
  {
    int c = a % b;
    a = b;
    b = c;
  }
  return a;
}

// When gcd(first, second) = 1, find Mult1st and Mult2nd such that:
// first * Mult1st + second * Mult2nd = 1.
static void ExtendedGCD(int first, int second, int* pMult1st, int* pMult2nd)
{
  int a = first;
  int b = second;
  int mult1stOld = 1;
  int mult1st = 0;
  int mult2ndOld = 0;
  int mult2nd = 1;

  while (b != 0)
  {
    int quot = a / b;
    int temp = b;
    b = a - (quot * b);
    a = temp;
    temp = mult1st;
    mult1st = mult1stOld - (quot * mult1st);
    mult1stOld = temp;
    temp = mult2nd;
    mult2nd = mult2ndOld - (quot * mult2nd);
    mult2ndOld = temp;
  }
  *pMult1st = mult1stOld;
  *pMult2nd = mult2ndOld;
}

static void showRatString(const char* num, const char* den)
{
  if (pretty != PARI_GP)
  {
    showRatConstants(num, den);
  }
  else
  {
    showText("(");
    showText(num);
    showText("/");
    showText(den);
    showText(")");
  }
}

static void ParseExpression(const char* expr)
{
  const char* ptrExpr = expr;
  while (*ptrExpr != 0)
  {
    if (*ptrExpr == 'S')
    {    // Square root of next digits.
      if (pretty != PARI_GP)
      {
        startSqrt();
      }
      for (;;)
      {
        ptrExpr++;
        char c = *ptrExpr;
        if ((c < '0') || (c > '9'))
        {
          break;
        }
        *ptrOutput = c;
        ptrOutput++;
      }
      ptrExpr--;
      if (pretty != PARI_GP)
      {
        endSqrt();
      }
      else
      {
        showText("^(1/2)");
      }
    }
    else if (*ptrExpr == '[')
    {  // Start of square root.
      startSqrt();
    }
    else if (*ptrExpr == ']')
    {  // End of square root.
      endSqrt();
    }
    else if (*ptrExpr == '-')
    {
      showText(ptrMinus);
    }
    else if (*ptrExpr == '*')
    {
      showText(ptrTimes);
    }
    else if (*ptrExpr == '(')
    {
      startParen();
    }
    else if (*ptrExpr == ')')
    {
      endParen();
    }
    else
    {
      *ptrOutput = *ptrExpr;
      ptrOutput++;
    }
    ptrExpr++;
  }
}

// Show cos(num*pi/den) as radicals.
// The output is the sign and the value of the denominator, or zero
// if the result is zero.
static int showRadicals(int numerator, int denominator, int multipl,
  int powerOf2, const char *times)
{
  int multiple = multipl;
  int power2 = powerOf2;
  int num = numerator;
  int den = denominator;
  int arraySigns[20];
  int indexSigns = 0;
  int exprDen;
  int den2 = 2 * den;
  int angle;
  const char* ptrExpr;
  int sign;
  int mult;
  int result;
  char denom[15];
  num = num % den2;  // Convert to range 0 to 360 degrees.
  if (num < 0)
  {
    num += den2;
  }
  if (num == 0)
  {
    return 1;
  }
  if (num == den)
  {
    return -1;
  }
  if (((num * 2) == den) || ((num*2) == (den*3)))
  {
    return 0;
  }
  while (((num % 2) == 0) && ((den % 2) == 0))
  {
    num /= 2;
    den /= 2;
    power2--;
  }
  if (multiple == 1)
  {
    power2 -= 2;
    multiple = 4;
    angle = num % 8;
    if ((angle == 1) || (angle == 7))
    {    // 45 or 315 degrees.
      sign = 1;
    }
    else
    {
      sign = -1;
    }
    if (power2 > 0)
    {
      if (sign > 0)
      {
        ptrExpr = "2+S2";   // 2+2*cos(Pi/4)
      }
      else
      {
        ptrExpr = "2-S2";   // 2-2*cos(Pi/4)
      }
    }
    else
    {
      ptrExpr = "S2";       // cos(Pi/4)
    }
    exprDen = 2;
  }
  else
  {  // Obtain cosine of num*pi/multiple.
    int indexExpr;
    angle = (2*num * 15 / multiple) % 60;
    if (angle < 15)
    {     // Angle between 0 and 90 degrees.
      indexExpr = (angle - 1) / 2;
      sign = 1;
    }
    else if (angle < 30)
    {     // Angle between 90 and 180 degrees.
      indexExpr = (30 - angle - 1) / 2;
      sign = -1;
    }
    else if (angle < 45)
    {     // Angle between 180 and 270 degrees.
      indexExpr = (angle - 30 - 1) / 2;
      sign = -1;
    }
    else
    {     // Angle between 270 and 360 degrees.
      indexExpr = (60 - angle - 1) / 2;
      sign = 1;
    }
    if (power2 > 0)
    {
      if (sign > 0)
      {
        ptrExpr = a2Plus2Cosine12[indexExpr].radicands;
        exprDen = a2Plus2Cosine12[indexExpr].denominator;
      }
      else
      {
        ptrExpr = a2Minus2Cosine12[indexExpr].radicands;
        exprDen = a2Minus2Cosine12[indexExpr].denominator;
      }
    }
    else
    {
      ptrExpr = aCosine12[indexExpr].radicands;
      exprDen = aCosine12[indexExpr].denominator;
    }
  }
  mult = 1;
  for (indexSigns = power2 - 1; indexSigns >= 0; indexSigns--)
  {
    if (indexSigns >= ((int)sizeof(arraySigns) / (int)sizeof(arraySigns[0])))
    {    // This cannot occur.
      return 0;
    }
    arraySigns[indexSigns] = sign;
    angle = num * 15 / multiple % (60*mult);
    if ((angle < (15*mult)) || (angle > (45*mult)))
    {     // Angle between 0 and 90 degrees.
          // or from 270 to 360 degrees.
      sign = 1;
    }
    else
    {     // Angle between 90 and 270 degrees.
      sign = -1;
    }
    mult *= 2;
  }
  for (indexSigns = 0; indexSigns < power2; indexSigns++)
  {
    if (indexSigns == (power2 - 1))
    {
      if ((indexSigns > 0) && (exprDen != 2))
      {
        char* ptrDenom = denom;
        if (exprDen == 1)
        {
          *ptrOutput = '2';
          ptrOutput++;
        }
        else
        {
          int2dec(&ptrDenom, exprDen / 2);
          showRatString("1", denom);
        }
        showText(ptrTimes);
      }
      startSqrt();
      break;
    }
    startSqrt();
    *ptrOutput = '2';
    ptrOutput++;
    if (arraySigns[indexSigns] < 0)
    {
      *ptrOutput = ' ';
      ptrOutput++;
      showText(ptrMinus);
      *ptrOutput = ' ';
      ptrOutput++;
    }
    else
    {
      showText(" + ");
    }
  }
  // Interpret expression.
  if ((power2 == 0) && (strcmp(ptrExpr, "1") == 0))
  {
    return sign * exprDen;
  }
  ParseExpression(ptrExpr);
  for (indexSigns = 0; indexSigns < power2; indexSigns++)
  {
    endSqrt();
  }
  result = ((power2 > 1)? 2 : exprDen);
  if (strcmp(ptrExpr, "1") && (result != 1))
  {
    showText(times);
  }
  return sign*result;
}

// Show cos(numerator34*Pi/34)
static int showRadicals17(int numerator34)
{
  int index;
  int angle2;
  int angle = numerator34 % 68;        // Convert to range 0 to 360 degrees.
  if (angle < 0)
  {
    angle += 68;
  }
  assert((angle != 17) && (angle != 51)); // It cannot be 90 deg or 270 deg.
  if ((angle % 2) == 0)
  {
    int sign;
    angle /= 2;                       // Show cos(angle*Pi/17).
    if (angle < 9)
    {           // Range 0 to 90 degrees.
      index = angle - 1;
      sign = 1;
    }
    else if (angle < 17)
    {           // Range 90 to 180 degrees.
      index = 17 - angle - 1;
      sign = -1;
    }
    else if (angle < 26)
    {           // Range 180 to 270 degrees.
      index = angle - 17 - 1;
      sign = -1;
    }
    else
    {           // Range 270 to 360 degrees.
      index = 34 - angle - 1;
      sign = 1;
    }
    ParseExpression(aCosine17[index].radicands);
    return aCosine17[index].denominator * sign;
  }
  angle2 = angle % 34;
  if (angle2 < 9)
  {           // Range 0 to 90 degrees. Cosine positive.
    index = angle2 - 1;
  }
  else if (angle2 < 17)
  {           // Range 90 to 180 degrees. Cosine negative.
    index = 17 - angle2 - 1;
  }
  else if (angle2 < 26)
  {           // Range 180 to 270 degrees. Cosine negative.
    index = angle2 - 17 - 1;
  }
  else
  {           // Range 270 to 360 degrees. Cosine positive.
    index = 34 - angle2 - 1;
  }
  startSqrt();
  ParseExpression(a2PM2Cosine17[index].radicands);
  endSqrt();
  if ((angle < 17) || (angle > 51))
  {
    return 8;
  }
  return -8;
}

static void AdjustComponent(int denominator, char* ptrStart, int toShow,
  int isFirst, const char *realRoot)
{
  int denomin = denominator;
  char beginning[500];
  char* ptrBeginning = beginning;
  int lenBeginning;
  *ptrBeginning = 0;
  size_t diffPtrs;
  if (denomin == 0)
  {     // Discard all output if result is zero.
    ptrOutput = ptrStart;
    *ptrOutput = 0;
    return;
  }
  if (denomin < 0)
  {
    copyStr(&ptrBeginning, (pretty == PRETTY_PRINT)? "&minus;" : "-");
    denomin = -denomin;    // Make it positive.
  }
  else if ((toShow == SHOW_IMAG) || (isFirst == 0))
  {
    copyStr(&ptrBeginning, " + ");
  }
  else
  {          // Nothing to do.
  }
  if (denomin == 1)
  {
    if (toShow == SHOW_IMAG)
    {
      if (pretty == PARI_GP)
      {
        *ptrBeginning = 'I';
        ptrBeginning++;
      }
      else
      {
        *ptrBeginning = 'i';
        ptrBeginning++;
        *ptrBeginning = ' ';
        ptrBeginning++;
      }
      copyStr(&ptrBeginning, ptrTimes);
    }
  }
  else if (pretty != PARI_GP)
  {
    copyStr(&ptrBeginning, (pretty == TEX)? "\\frac{":
      "<span class=\"fraction\"><span class=\"numerator\">");
    if (pretty == PARI_GP)
    {
      *ptrBeginning = ((toShow == SHOW_REAL)? '1' : 'I');
    }
    else
    {
      *ptrBeginning = ((toShow == SHOW_REAL)? '1' : 'i');
    }
    ptrBeginning++;
    copyStr(&ptrBeginning, (pretty == TEX)? "}{": "</span><span class=\"denominator\">");
    int2dec(&ptrBeginning, denomin);
    copyStr(&ptrBeginning, (pretty == TEX)? "}": "</span></span> ");
    copyStr(&ptrBeginning, ptrTimes);
  }
  else
  {             // PARI-GP
    *ptrBeginning = '(';
    ptrBeginning++;
    *ptrBeginning = ((toShow == SHOW_REAL)? '1' : 'I');
    ptrBeginning++;
    *ptrBeginning = '/';
    ptrBeginning++;
    int2dec(&ptrBeginning, denomin);
    *ptrBeginning = ')';
    ptrBeginning++;
    *ptrBeginning = '*';
    ptrBeginning++;
    *ptrBeginning = 0;          // Add terminator at end of string.
  }
  copyStr(&ptrBeginning, realRoot);
  diffPtrs = ptrBeginning - &beginning[0];
  lenBeginning = (int)diffPtrs;
  if ((*realRoot != 0) && (lenBeginning != 0) && (*ptrStart != 0))
  {
    copyStr(&ptrBeginning, ptrTimes);
  }
  diffPtrs = ptrBeginning - &beginning[0];
  lenBeginning = (int)diffPtrs;
  (void)memmove(ptrStart + lenBeginning, ptrStart, strlen(ptrStart));
  (void)memcpy(ptrStart, beginning, lenBeginning);
  ptrOutput += lenBeginning;
  *ptrOutput = 0;
}

// If den is multiple of 17, compute the extended GCD of 17 and den/17.
// Let the results be A and B respectively. so num*pi/den = num*A*pi/17 +
// num*B*pi/(den/17). Use the formula cos(a+b) = cos a * cos b - sin a * sin b.
static void showComponent(int num, int den, int multipl, int power2, int toShow, const char *realRoot)
{
  int numer = num;
  int denom = den;
  int multiple = multipl;
  int denominCos;
  int den2 = 2 * denom;
  char * ptrStartRadicals = ptrOutput;

  *ptrOutput = 0;
  numer = numer % den2;  // Convert to range 0 to 360 degrees.
  if (numer < 0)
  {
    numer += den2;
  }
  if ((denom % 17) == 0)
  {           // Denominator is multiple of 17.
    int numerator34;
    int numeratorDen;
    int denominCos17;
    denom /= 17;
    multiple /= 17;
    ExtendedGCD(denom, 17, &numerator34, &numeratorDen);
    numerator34 *= 2*numer;
    numeratorDen *= numer;
    // The original angle is multDen*pi/den + mult17*pi/17 = A + B.
    // Show cos(A) * cos(B) - sin(A) * sin(B).
        // Show cos(A).
    denominCos = showRadicals(numeratorDen, denom, multiple, power2, ptrTimes);
    if (denominCos != 0)
    {   // Show cos(B).
      denominCos17 = showRadicals17(numerator34);
      AdjustComponent(denominCos * denominCos17, ptrStartRadicals, toShow, 1, realRoot);
    }
        // Show sin(A).
    ptrStartRadicals = ptrOutput;
    *ptrOutput = 0;
    if ((denom & 0x01) == 0x01)
    {           // denom is odd.
      denominCos = showRadicals(denom - (numeratorDen * 2), denom * 2, multiple, 1, ptrTimes);
    }
    else
    {           // denom is even.
      denominCos = showRadicals((denom / 2) - numeratorDen, denom, multiple, power2, ptrTimes);
    }
    if (denominCos != 0)
    {    // Show sin(B).
      denominCos17 = showRadicals17(17 - numerator34);
      AdjustComponent(-denominCos * denominCos17, ptrStartRadicals, toShow, 0, realRoot);
    }
    return;
  }
  denominCos = showRadicals(numer, denom, multiple, power2, "");
  AdjustComponent(denominCos, ptrStartRadicals, toShow, 1, realRoot);
}

// Try to find a radical expression for cos(num*pi/den) + 
// i*sin(num*pi/den). If it is not possible, return with no output.
static void outputRadicandsForCosSin(int num, int den, const char *realRoot)
{
  // den must be 3^e3 * 5^e5 * 17^e17 * 2^k, where the exponents can be 0 or 1.
  // At this moment we do not consider other Fermat prime numbers.
  int multiple = 1;
  int power2 = 0;
  int quot = den;
  if ((quot % 3) == 0)
  {
    multiple *= 3;
    quot /= 3;
  }
  if ((quot % 5) == 0)
  {
    multiple *= 5;
    quot /= 5;
  }
  if ((quot % 17) == 0)
  {
    multiple *= 17;
    quot /= 17;
  }
  // At this moment quot should be a power of 2.
  while ((quot % 2) == 0)
  {
    power2++;
    quot /= 2;
  }
  if (quot != 1)
  {                   // Denominator is multiple of another prime.
    return;           // Cannot express by radicals.
  }
  showText(" = ");

  showComponent(num, den, multiple, power2, SHOW_REAL, realRoot);
  if (power2 == 0)
  {
    showComponent(den-(num*2), den*2, multiple, 1, SHOW_IMAG, realRoot);
  }
  else
  {
    showComponent((den/2) - num, den, multiple, power2, SHOW_IMAG, realRoot);
  }
}

// Show multiplicand * cos(M) + i * multiplicand * sin(M)
// where M equals realNum*pi/realDen
static void showTrig(int numerator, int denominator, const char* multiplicand)
{
  char num[300];
  char den[300];
  char* ptrNum = den;
  int2dec(&ptrNum, denominator);
  *ptrNum = 0;  // Include string terminator.
  showText(multiplicand);
  if (*multiplicand != 0)
  {
    if (pretty != PARI_GP)
    {
      *ptrOutput = ' ';
      ptrOutput++;
    }
    showText(ptrTimes);
  }
  showText(ptrCos);
  ptrNum = num;
  if (numerator != 1)
  {
    int2dec(&ptrNum, numerator);
    *ptrOutput = ' ';
    ptrOutput++;
    copyStr(&ptrNum, ptrTimes);
  }
  copyStr(&ptrNum, ptrPi);
  showRatString(num, den);
  if (pretty == TEX)
  {
    *ptrOutput = '}';
    ptrOutput++;
  }
  showText(" + ");
  showText(ptrI);
  showText(" ");
  showText(ptrTimes);
  showText(multiplicand);
  if (*multiplicand != 0)
  {
    *ptrOutput = ' ';
    ptrOutput++;
    showText(ptrTimes);
    *ptrOutput = ' ';
    ptrOutput++;
  }
  showText(ptrSin);
  if (pretty != PARI_GP)
  {
    showRatConstants(num, den);
    if (pretty == TEX)
    {
      *ptrOutput = '}';
      ptrOutput++;
    }
  }
  else
  {
    showText("(");
    showText(num);
    showText("/");
    showText(den);
    showText(")");
  }
}

static bool TestCyclotomic(const int* ptrPolynomial, int multiplicity, int polyDegree)
{
  int index;
  const int* ptrCoeff;
  int currentDegree;
  if (polyDegree <= 0)
  {
    return false;   // Constant polynomials are not cyclotomic.
  }
  // Polynomial must be palindromic of even degree
  // and the absolute value must be less than 10.
  if ((polyDegree & 0x01) == 0x01)
  {      // Polynomial has odd degree so it cannot be cyclotomic.
    return false;
  }
  ptrCoeff = ptrPolynomial;
  for (index = 0; index <= polyDegree; index++)
  {
    if ((*ptrCoeff != 1) && (*ptrCoeff != -1))
    {            // Too low ot too high.
      return false;
    }
    if (*(ptrCoeff + 1) >= 10)
    {            // Absolute value too high.
      return false;
    }
    ptrCoeff += 2;
  }
  if (totients[0] == 0)
  {    // Table not initialized.
    initializeSmallPrimes(smallPrimes);
    totients[0] = 1;          // indicate table initialized
    totients[2] = 1;
    totients[3] = 2;
    // Compute totient[index] as index*(p1-1)/p1*(p2-1)/p2*...
    // where p_n is a prime factor of index.
    for (index = 4; index <= (2*MAX_DEGREE); index++)
    {
      int prime = 2;
      int primeIndex = 0;
      int quotient = index;
      int totient = index;
      while ((quotient != 1) && ((prime*prime) <= index))
      {
        if (((quotient / prime) * prime) == quotient)
        {   // Prime dividing index was found.
          totient = totient * (prime - 1) / prime;
          do
          {
            quotient /= prime;
          } while (((quotient / prime) * prime) == quotient);
        }
        primeIndex++;
        prime = smallPrimes[primeIndex];
      }
      if (quotient > 1)
      {
        totient = totient * (quotient - 1) / quotient;
      }
      totients[index] = totient;
    }
  }
  // Degree of cyclotomic polynomial n is totient(n).
  for (index = polyDegree; index <= (2 * MAX_DEGREE); index++)
  {
    if (polyDegree != totients[index])
    {   // Degree is not correct.
      continue;
    }
    // Test whether x^degree - 1 divides this polynomial.
    int base[MAX_DEGREE];
    int prod[MAX_DEGREE];
    int lenBytes;
    int* ptrBase = base;
    ptrCoeff = ptrPolynomial;
    for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
    {
      *ptrBase = *ptrCoeff* *(ptrCoeff + 1);
      ptrBase++;
      ptrCoeff += 2;
    }
    lenBytes = polyDegree * (int)sizeof(int);
    (void)memset(prod, 0, lenBytes);
    prod[1] = 1;    // Initialize polynomial to x.
    for (int expon = 1; expon < index; expon++)
    {               // Multiply by x.
      int coeffMaxDegree = prod[polyDegree - 1];
      for (currentDegree = polyDegree - 1; currentDegree > 0; currentDegree--)
      {
        prod[currentDegree] = prod[currentDegree-1] - (coeffMaxDegree * base[currentDegree]);
      }
      prod[0] = -coeffMaxDegree * base[0];
    }
    // If result is 1, then the polynomial is cyclotomic.
    prod[0]--;
    lenBytes = polyDegree * (int)sizeof(int);
    (void)memset(base, 0, lenBytes);
    if (memcmp(base, prod, lenBytes) == 0)
    {     // The polynomial is cyclotomic.
      bool denIsOdd = ((index & 0x01) == 0x01)? true: false;
      int realDen = denIsOdd ? index : (index / 2);
      for (int numerator = 1; numerator < index; numerator++)
      {
        int realNum;
        if (gcd(numerator, index) > 1)
        {        // Numbers are not coprime. Next loop.
          continue;
        }
        showX(multiplicity);
        if (denIsOdd)
        {
          realNum = 2 * numerator;
        }
        else
        {
          realNum = numerator;
        }
        // Show cos(M) + i sin(M)
        // where M equals realNum*pi/realDen
        showTrig(realNum, realDen, "");
        outputRadicandsForCosSin(realNum, realDen, "");
        endLine();
      }
      return true;
    }
  }
  return false;         // polynomial is not cyclotomic.
}

static bool isPalindromic(int* ptrPolynomial, int polyDegree)
{
  int currentDegree;
  int *ptrs[MAX_DEGREE+1];
  int* ptrCoeff = ptrPolynomial;
  // Get pointers to coefficients.
  for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {
    ptrs[currentDegree] = ptrCoeff;
    ptrCoeff += numLimbs(ptrCoeff);
    ptrCoeff++;
  }
  for (currentDegree = 0; (currentDegree * 2) < polyDegree; currentDegree++)
  {
    if (*ptrs[currentDegree] != *ptrs[polyDegree - currentDegree])
    {
      return false;    // Polynomial is not palindromic.
    }
    if (memcmp(ptrs[currentDegree] + 1, ptrs[polyDegree - currentDegree] + 1,
      numLimbs(ptrs[currentDegree])) != 0)
    {
      return false;    // Polynomial is not palindromic.
    }
  }
  return true;        // Polynomial is palindromic.
}

static void StartRadicand(int polyDegree)
{
  if (pretty == PRETTY_PRINT)
  {
    copyStr(&ptrOutput, "<span class=\"");
    if (polyDegree < 10)
    {       // degree has 1 digit.
      copyStr(&ptrOutput, "root");
    }
    else if (polyDegree < 100)
    {       // degree has 2 digits.
      copyStr(&ptrOutput, "root2dig");
    }
    else
    {       // degree has 3 digits.
      copyStr(&ptrOutput, "root3dig");
    }
    copyStr(&ptrOutput, "\"><span class=\"befrad\">");
    int2dec(&ptrOutput, polyDegree);
    copyStr(&ptrOutput, "</span><span class=\"radicand2\">");
  }
  if (pretty == TEX)
  {
    copyStr(&ptrOutput, "\\sqrt[");
    int2dec(&ptrOutput, polyDegree);
    copyStr(&ptrOutput, "]{");
  }
}

static void EndRadicand(int polyDegree)
{
  if (pretty == PRETTY_PRINT)
  {
    copyStr(&ptrOutput, "</span></span>");
  }
  else if (pretty == TEX)
  {
    *ptrOutput = '}';
    ptrOutput++;
  }
  else
  {
    *ptrOutput = '^';
    ptrOutput++;
    *ptrOutput = '(';
    ptrOutput++;
    *ptrOutput = '1';
    ptrOutput++;
    *ptrOutput = '/';
    ptrOutput++;
    int2dec(&ptrOutput, polyDegree);
    *ptrOutput = ')';
    ptrOutput++;
  }
  *ptrOutput = 0;
}

static void GenerateRoots(int multiplicity, const char* rationalRoot,
  bool isNegative, int polyDegree)
{
  for (int currentDegree = 0; currentDegree < polyDegree; currentDegree++)
  {
    int realNum;
    int realDen;
    char realRoot[30000];
    char* ptrOutputBak = ptrOutput;
    ptrOutput = realRoot;
    int numer = (2 * currentDegree) + (isNegative ? 1 : 0);
    int gcdNumDen = gcd(numer, polyDegree);
    realNum = numer / gcdNumDen;
    realDen = polyDegree / gcdNumDen;
    StartRadicand(polyDegree);
    if (pretty == PARI_GP)
    {
      *ptrOutput = '(';
      ptrOutput++;
    }
    copyStr(&ptrOutput, rationalRoot);
    if (pretty == PARI_GP)
    {
      *ptrOutput = ')';
      ptrOutput++;
    }
    EndRadicand(polyDegree);
    ptrOutput = ptrOutputBak;
    showX(multiplicity);
    showTrig(realNum, realDen, realRoot);
    outputRadicandsForCosSin(realNum, realDen, realRoot);
    endLine();
  }
}

static void ShowRootsOfRationalNumbers(int polyDegree, int multiplicity)
{
  char rationalRoot[30000];
  char* ptrOutputBak = ptrOutput;
  bool isNegative = false;
  if (Rat1.numerator.sign == Rat1.denominator.sign)
  {
    isNegative = true;
  }
  Rat1.numerator.sign = SIGN_POSITIVE;
  Rat1.denominator.sign = SIGN_POSITIVE;
  ptrOutput = rationalRoot;
  showRational(&Rat1);
  *ptrOutput = 0;
  ptrOutput = ptrOutputBak;
  // Generate roots.
  GenerateRoots(multiplicity, rationalRoot, isNegative, polyDegree);
}

// If polynomial has the form ax^n + b, show its roots.
static bool isLinearExponential(const int* ptrPolynomial, int polyDegree, int multiplicity)
{
  const int* ptrPoly = ptrPolynomial;
  // Get constant term.
  NumberLength = *ptrPoly;
  IntArray2BigInteger(ptrPoly, &Rat1.numerator);
  ptrPoly += numLimbs(ptrPoly);
  ptrPoly++;
  for (int currentDegree=1; currentDegree < polyDegree; currentDegree++)
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
  IntArray2BigInteger(ptrPoly, &Rat1.denominator);
  ShowRootsOfRationalNumbers(polyDegree, multiplicity);
  return true;
}

// If polynomial has the form ax^(2n) + bx^n + c, show its roots.
static bool isQuadraticExponential(const int* ptrPolynomial, int polyDegree,
  int multiplicity)
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
  if (ProcessQuadraticEquation(&signDiscr) == 0)
  {           // Discriminant is perfect square. Roots are rational numbers.
    for (ctr = 0; ctr < 2; ctr++)
    {         // Compute (-linear +/- sqrt(discr))/(2*quadratic)
      if (ctr == 0)
      {
        BigIntAdd(&Linear, &tmp4, &tmp1);
      }
      else
      {
        BigIntSubt(&Linear, &tmp4, &tmp1);
      }
      CopyBigInt(&Rat1.numerator, &tmp1);
      CopyBigInt(&Rat1.denominator, &Quadratic);
      BigRationalDivideByInt(&Rat1, -2, &Rat1);
      ShowRootsOfRationalNumbers(halfDegree, multiplicity);
    }
    return true;
  }
  if (signDiscr == SIGN_POSITIVE)
  {           // Roots of quadratic equation are real.
    enum eSign Rat2SignBak;
    for (ctr = 0; ctr < 2; ctr++)
    {         // Show (Rat1 +/- Rat2*sqrt(Rat3))
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
          copyStr(&ptrOutput, (pretty==PARI_GP)? " + I ": " + i ");
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
          copyStr(&ptrOutput, (pretty == PRETTY_PRINT)? "&minus; " : "- ");
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
      endLine();
    }
  }
  return true;
}

// Save factor degrees sorting by ascending order.
static void SaveFactorDegrees(int prime, int *factors, int nbrFactors)
{
  int* piFactors = factors;
  int* ptrFactors;
  int* ptrOldFactor;
  int tmpDegree;

  *piFactors = prime;
  piFactors++;
  *piFactors = nbrFactors;
  piFactors++;
  ptrFactors = piFactors;
  for (int currentFactor = 0; currentFactor < nbrFactors; currentFactor++)
  {
    int newDegree = factorInfo[currentFactor].degree;
    for (ptrOldFactor = ptrFactors; ptrOldFactor < piFactors; ptrOldFactor++)
    {
      if (newDegree < *ptrOldFactor)
      {
        for (; ptrOldFactor < piFactors; ptrOldFactor++)
        {
          tmpDegree = *ptrOldFactor;
          *ptrOldFactor = newDegree;
          newDegree = tmpDegree;
        }
        break;
      }
    }
    *ptrOldFactor = newDegree;
    piFactors++;
  }
}

static void showDegrees(const int *factors)
{
  const int* piFactors = factors;
  piFactors++;
  int nbrFactors = *piFactors;
  for (int currentFactor = 0; currentFactor < nbrFactors; currentFactor++)
  {
    if (currentFactor != 0)
    {
      if (currentFactor == (nbrFactors - 1))
      {
        showText(lang ? " y " : " and ");
      }
      else
      {
        *ptrOutput = ',';
        ptrOutput++;
        *ptrOutput = ' ';
        ptrOutput++;
      }
    }
    piFactors++;
    int2dec(&ptrOutput, *piFactors);
  }
}

static void ShowNoSolvableSpanish(const int* firstArray, const int* secondArray,
  int numberDifferentX, int nbrFactor, int gcdDegrees)
{
  if (firstArray != secondArray)
  {
    if (secondArray == NULL)
    {
      showX(numberDifferentX);
      *(ptrOutput - 2) = ':';  // Replace equal sign by colon.
      showText("Las raíces del polinomio");
      if (nbrFactor >= 0)
      {
        showText(" número ");
        int2dec(&ptrOutput, nbrFactor + 1);
      }
      showText(" no se pueden expresar mediante expresiones radicales");
      if (gcdDegrees > 1)
      {
        showText(". Reemplazamos <var>y</var> = ");
        showPowerX(&ptrOutput, gcdDegrees);
      }
      showText(". Los");
    }
    else
    {
      showText(" y los");
    }
    showText(" grados de los factores del polinomio módulo ");
    int2dec(&ptrOutput, *firstArray);
    showText(" son ");
    showDegrees(firstArray);
  }
  showText(" (el grupo de Galois contiene un ciclo de longitud ");
}

static void ShowNoSolvableEnglish(const int* firstArray, const int* secondArray,
  int numberDifferentX, int nbrFactor, int gcdDegrees)
{
  if (firstArray != secondArray)
  {
    if (secondArray == NULL)
    {
      showX(numberDifferentX);
      *(ptrOutput - 2) = ':';  // Replace equal sign by colon.
      showText("The roots of the polynomial");
      if (nbrFactor >= 0)
      {
        showText(" number ");
        int2dec(&ptrOutput, nbrFactor + 1);
      }
      showText(" cannot be expressed by radicals");
      if (gcdDegrees > 1)
      {
        showText(". We set <var>y</var> = ");
        showPowerX(&ptrOutput, gcdDegrees);
      }
      showText(". The");
    }
    else
    {
      showText(" and the");
    }
    showText(" degrees of the factors of polynomial modulo ");
    int2dec(&ptrOutput, *firstArray);
    showText(" are ");
    showDegrees(firstArray);
  }
  showText(" (the Galois group contains a cycle of ");
}

static void showExplanation(int left, const char* oper1, int middle, const char* oper2, int right)
{
  *ptrOutput = ' ';
  ptrOutput++;
  *ptrOutput = '(';
  ptrOutput++;
  int2dec(&ptrOutput, left);
  *ptrOutput = ' ';
  ptrOutput++;
  showText(oper1);
  *ptrOutput = ' ';
  ptrOutput++;
  int2dec(&ptrOutput, middle);
  *ptrOutput = ' ';
  ptrOutput++;
  showText(oper2);
  *ptrOutput = ' ';
  ptrOutput++;
  int2dec(&ptrOutput, right);
  *ptrOutput = ')';
  ptrOutput++;
}

static bool isPrime(int value)
{
  int divisor = 3;
  if (value == 2)
  {
    return true;
  }
  if ((value & 1) == 0)
  {
    return false;        // Even value different from 2: composite.
  }
  while ((divisor * divisor) <= value)
  {
    if ((value % divisor) == 0)
    {
      return false;      // Composite.
    }
    divisor += 2;
  }
  return true;           // Prime value.
}

// If polynomial is S_n or A_n, indicate that the roots are not solvable.
// Factor the polynomial modulo different primes less than 100.
// If with the same or different factorizations, there is a factor of
// degree 2 or 3 and there is no factor of degree multiple of 2 or 3 respectively
// and if there is a factor with prime degree greater than half the degree of
// the polynomial, the polynomial is not solvable with radicals.
// Discard the primes such that the factorizations give duplicated roots.
// If the factor of prime degree is n/2 < p < n-2, or if n/3 < p <= n/2 and
// the largest degree is odd and greater than n/2, then the polynomial is not
// solvable with radicals.
static bool isSymmetricOrAlternating(int nbrFactor, const int* ptrPolynomial, 
  int polyDeg, int multiplicity)
{
  int polyDegree = polyDeg;
  const char* ptrStrPrimeHalfSp = "primo mayor que la mitad del grado del polinomio";
  const char* ptrStrPrimeHalfEn = "prime length greater than half the degree of polynomial";
  int factorDegreesCycle2Or3[MAX_DEGREE+2];
  int factorDegreesCycleP[MAX_DEGREE + 2];  // n/2 < p
  int factorDegreesCycleOther[MAX_DEGREE + 2]; // n/3 < p < n/2
  int cycle2Found = 0;
  int cycle3Found = 0;
  int cyclePrGtNOver2Found = 0;
  int cyclePrGtNOver2ToLess2Found = 0;
  int cycleOddGtNOver2Found = 0;
  int cyclePrGtNOver3Found = 0;
  int prime;
  int primeIndex = -1;
  int gcdDegrees = 0;
  int currentDegree;
  const int* ptrCoeff;
  int* ptrCoeffDest;
  factorDegreesCycleP[1] = 0;
  // Compute GCD of all degrees of coefficients different from zero.
  // Generate polynomial polyNonRepeatedFactors stripping all zero
  // coefficients with degree not multiple of this GCD.
  ptrCoeff = ptrPolynomial;
  for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {
    if ((*ptrCoeff != 1) || (*(ptrCoeff + 1) != 0))
    {            // Coefficient is not zero.
      gcdDegrees = gcd(currentDegree, gcdDegrees);
      if (gcdDegrees == 1)
      {          // GCD of degrees is one. Further loops will not change this.
        break;
      }
    }
    ptrCoeff += numLimbs(ptrCoeff);
    ptrCoeff++;
  }
  if (gcdDegrees == 0)
  {    // No gcd computed. Go out.
    return false;
  }
  polyDegree /= gcdDegrees;
  if (polyDegree < 5)
  {     // Polynomial is solvable with radical expressions.
    showX(multiplicity * polyDegree * gcdDegrees);
    *(ptrOutput - 2) = ':';  // Replace equal sign by colon.
    if (lang)
    {
      showText("Las raíces del polinomio");
      if (nbrFactor >= 0)
      {
        showText(" número ");
        int2dec(&ptrOutput, nbrFactor + 1);
      }
      showText(" se pueden expresar mediante expresiones radicales");
      if (gcdDegrees > 1)
      {
        showText(". Reemplazamos <var>y</var> = ");
        showPowerX(&ptrOutput, gcdDegrees);
      }
      showText(". El polinomio es de grado ");
      int2dec(&ptrOutput, polyDegree);
      showText(" que es menor que 5.");
    }
    else
    {
      showText("The roots of the polynomial");
      if (nbrFactor >= 0)
      {
        showText(" number ");
        int2dec(&ptrOutput, nbrFactor + 1);
      }
      showText(" can be expressed by radicals");
      if (gcdDegrees > 1)
      {
        showText(". We set <var>y</var> = ");
        showPowerX(&ptrOutput, gcdDegrees);
      }
      showText(". The polynomial has degree ");
      int2dec(&ptrOutput, polyDegree);
      showText(" which is less than 5.");
    }
    return true;
  }
  ptrCoeff = ptrPolynomial;
  ptrCoeffDest = polyNonRepeatedFactors;
  *ptrCoeffDest = polyDegree;
  ptrCoeffDest++;
  for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {    // Copy coefficient and skip zero coefficients.
    int nbrLen = 1 + numLimbs(ptrCoeff);
    int lenBytes = nbrLen * (int)sizeof(int);
    int offset = nbrLen + (2 * (gcdDegrees - 1));
    (void)memcpy(ptrCoeffDest, ptrCoeff, lenBytes);
    ptrCoeff += offset;
    ptrCoeffDest += nbrLen;
  }
  do
  {
    int nbrFactors;
    int cycle2FoundInThisFactor = 0;
    int cycle3FoundInThisFactor = 0;
    const struct sFactorInfo* pstFactorInfo;
    (void)memset(factorInfo, 0, sizeof(factorInfo));
    primeIndex = getNextPrimeNoDuplicatedFactors(primeIndex);
    prime = smallPrimes[primeIndex];
    FactorPolynomialModPrime(prime);
    pstFactorInfo = factorInfo;
    // Check whether there is a factor of degree 2, 3 or prime
    // greater than half the degree of the polynomial being analyzed.
    for (nbrFactors = 0; nbrFactors < MAX_DEGREE; nbrFactors++)
    {
      if (pstFactorInfo->ptr == NULL)
      {    // No more factors.
        break;
      }
      pstFactorInfo++;
    }
    pstFactorInfo = factorInfo;
    for (int factorNbr = 0; factorNbr < nbrFactors; factorNbr++)
    {
      int currDegree;
      currDegree = pstFactorInfo->degree;
      if ((currDegree % 2) == 0)
      {
        if (currDegree == 2)
        {
          cycle2FoundInThisFactor++;
        }
        else
        {
          cycle2FoundInThisFactor = 2;
        }
      }
      if ((currDegree % 3) == 0)
      {
        if (currDegree == 3)
        {
          cycle3FoundInThisFactor++;
        }
        else
        {
          cycle3FoundInThisFactor = 2;
        }
      }
      if ((currDegree > 3) && (currDegree > (degree / 2)))
      {
        if (isPrime(currDegree))
        {      // Current degree > n/2 and is prime.
          if (currDegree < (degree - 2))
          {
            SaveFactorDegrees(prime, factorDegreesCycleP, nbrFactors);
            cyclePrGtNOver2ToLess2Found = currDegree;
          }
          else if (cyclePrGtNOver2Found == 0)
          {    // If first condition holds, the polynomial is not solvable.
            SaveFactorDegrees(prime, factorDegreesCycleP, nbrFactors);
            cyclePrGtNOver2Found = currDegree;
            cycleOddGtNOver2Found = currDegree;
          }
          else
          {              // Nothing to do.
          }
        }
        else if (((currDegree & 0x01) == 0x01) && (currDegree < degree))
        {      // Current degree > n/2, it is odd and less than the polynomial degree.
          if (cycleOddGtNOver2Found == 0)
          {    // If first condition holds, the polynomial is not solvable.
               // Test that this degree is coprime to all other degrees of
               // factors of this polynomial.
            int factNbr;
            pstFactorInfo = factorInfo;
            for (factNbr = 0; factNbr < nbrFactors; factNbr++)
            {
              if ((pstFactorInfo->degree < currDegree) && 
                (gcd(pstFactorInfo->degree, currDegree) != 1))
              {
                break;
              }
              pstFactorInfo++;
            }
            if (factNbr == nbrFactors)
            {      // All degrees are coprime to current degree.
              SaveFactorDegrees(prime, factorDegreesCycleP, nbrFactors);
              cycleOddGtNOver2Found = currDegree;
            }
          }
        }
        else
        {            // Nothing to do.
        }
      }
      else if ((currDegree > 3) && ((currDegree > degree) / 3))
      {
        if (isPrime(currDegree) && ((degree % currDegree) != 0))
        {      // Current degree > n/3, it is prime, and it does not divide the degree.
               // Ensure that only this degree is multiple of itself.
          int nbrMultiples = 0;
          pstFactorInfo = factorInfo;
          for (int factNbr = 0; factNbr < nbrFactors; factNbr++)
          {
            if ((pstFactorInfo->degree % currDegree) == 0)
            {
              nbrMultiples++;
            }
            pstFactorInfo++;
          }
          if (nbrMultiples == 1)
          {           // Only one multiple of currentDegree expected.
            SaveFactorDegrees(prime, factorDegreesCycleOther, nbrFactors);
            cyclePrGtNOver3Found = currDegree;
            if ((degree & 0x01) == 0x01)
            {         // If degree is odd, the group is very transitive.
              break;
            }
          }
        }
      }
      else
      {                // Nothing to do.
      }
      pstFactorInfo++;
    }
    if ((cycle2Found == 0) && (cycle3Found == 0))
    {
      if (cycle2FoundInThisFactor == 1)
      {
        cycle2Found = 1;
        SaveFactorDegrees(prime, factorDegreesCycle2Or3, nbrFactors);
      }
      else if (cycle3FoundInThisFactor == 1)
      {
        cycle3Found = 1;
        SaveFactorDegrees(prime, factorDegreesCycle2Or3, nbrFactors);
      }
      else
      {   // No more cases.
      }
    }
    if (cyclePrGtNOver2ToLess2Found != 0)
    {           // Group is very transitive.
      break;
    }
    if ((cyclePrGtNOver3Found != 0) && 
       (((degree & 0x01) == 0x01) || (cycleOddGtNOver2Found != 0) ||
         (cyclePrGtNOver2Found != 0)))
    {           // Group is very transitive.
      break;
    }
    if (((cycle2Found != 0) || (cycle3Found != 0)) && (cyclePrGtNOver2Found != 0))
    {           // Polynomial is not solvable with radicals. Exit loop.
      break;
    }
  } while (prime < 100);
  int numberDifferentX = multiplicity * polyDegree * gcdDegrees;
  if (cyclePrGtNOver2ToLess2Found != 0)
  {      // Group is very transitive.
    if (lang)
    {    // Spanish
      ShowNoSolvableSpanish(factorDegreesCycleP, NULL, numberDifferentX,
        nbrFactor, gcdDegrees);
      showText(ptrStrPrimeHalfSp);
      showExplanation(cyclePrGtNOver2ToLess2Found, "&gt;", degree, "&divide;", 2);
      showText(" y menor que el grado menos 2");
    }
    else
    {   // English
      ShowNoSolvableEnglish(factorDegreesCycleP, NULL, numberDifferentX,
        nbrFactor, gcdDegrees);
      showText(ptrStrPrimeHalfEn);
      showExplanation(cyclePrGtNOver2ToLess2Found, "&gt;", degree, "&divide;", 2);
      showText(" and less than the degree minus 2");
    }
    showExplanation(cyclePrGtNOver2ToLess2Found, "&lt;", degree, "&minus;", 2);
    *ptrOutput = ')';
    ptrOutput++;
  }
  else if ((cyclePrGtNOver3Found != 0) && ((degree & 0x01) == 0x01))
  {     // Group is very transitive.
    if (lang)
    {    // Spanish
      ShowNoSolvableSpanish(factorDegreesCycleOther, NULL, numberDifferentX,
        nbrFactor, gcdDegrees);
      showText("primo mayor que la tercera parte del grado del polinomio");
      showExplanation(cyclePrGtNOver3Found, "&gt;", degree, "&divide;", 3);
      showText(" y éste (");
      int2dec(&ptrOutput, degree);
      showText(") es impar)");
    }
    else
    {   // English
      ShowNoSolvableEnglish(factorDegreesCycleOther, NULL, numberDifferentX,
        nbrFactor, gcdDegrees);
      showText("prime length greater than a third of the degree of the polynomial");
      showExplanation(cyclePrGtNOver3Found, "&gt;", degree, "&divide;", 3);
      showText(" and the latter (");
      int2dec(&ptrOutput, degree);
      showText(") is odd)");
    }
  }
  else if ((cyclePrGtNOver3Found != 0) && (cycleOddGtNOver2Found != 0))
  {     // Group is very transitive.
    if (lang)
    {    // Spanish
      ShowNoSolvableSpanish(factorDegreesCycleOther, NULL, numberDifferentX,
        nbrFactor, gcdDegrees);
      showText("primo mayor que la tercera parte del grado del polinomio");
      showExplanation(cyclePrGtNOver3Found, "&gt;", degree, "&divide;", 3);
      *ptrOutput = ')';
      ptrOutput++;
      ShowNoSolvableSpanish(factorDegreesCycleP, factorDegreesCycleOther, numberDifferentX,
        nbrFactor, gcdDegrees);
      showText("impar mayor que la mitad del grado del polinomio");
    }
    else
    {    // English
      ShowNoSolvableEnglish(factorDegreesCycleOther, NULL, numberDifferentX,
        nbrFactor, gcdDegrees);
      showText("prime length greater than a third of the degree of the polynomial");
      showExplanation(cyclePrGtNOver3Found, "&gt;", degree, "&divide;", 3);
      *ptrOutput = ')';
      ptrOutput++;
      ShowNoSolvableEnglish(factorDegreesCycleP, factorDegreesCycleOther, numberDifferentX,
        nbrFactor, gcdDegrees);
      showText("odd length greater than half the degree of the polynomial");
    }
    showExplanation(cycleOddGtNOver2Found, "&gt;", degree, "&divide;", 2);
    *ptrOutput = ')';
    ptrOutput++;
  }
  else if (((cycle2Found != 0) || (cycle3Found != 0)) && (cyclePrGtNOver2Found != 0))
  {
    if (lang)
    {    // Spanish
      ShowNoSolvableSpanish(factorDegreesCycle2Or3, NULL, numberDifferentX,
        nbrFactor, gcdDegrees);
      *ptrOutput = (cycle2Found ? '2' : '3');
      ptrOutput++;
      *ptrOutput = ')';
      ptrOutput++;
      ShowNoSolvableSpanish(factorDegreesCycleP, factorDegreesCycle2Or3, numberDifferentX,
        nbrFactor, gcdDegrees);
      showText(ptrStrPrimeHalfSp);
    }
    else
    {    // English
      ShowNoSolvableEnglish(factorDegreesCycle2Or3, NULL, numberDifferentX,
        nbrFactor, gcdDegrees);
      showText("length ");
      *ptrOutput = (cycle2Found ? '2' : '3');
      ptrOutput++;
      *ptrOutput = ')';
      ptrOutput++;
      ShowNoSolvableEnglish(factorDegreesCycleP, factorDegreesCycle2Or3, numberDifferentX,
        nbrFactor, gcdDegrees);
      showText(ptrStrPrimeHalfEn);
    }
    showExplanation(cyclePrGtNOver2Found, "&gt;", degree, "&divide;", 2);
    *ptrOutput = ')';
    ptrOutput++;
  }
  else
  {
    return false;
  }
  return true;
}

void getRootsPolynomial(int nbrFactor, char **pptrOutput, struct sFactorInfo* pstFactorInfo, int groupLength)
{
  int nbrFactorsFoundBak = nbrFactorsFound;
  static struct sFactorInfo factorInfoIntegerBak[MAX_DEGREE];
  static int polyIntegerBak[1000000];
  int multiplicity = pstFactorInfo->multiplicity;
  if (lang)
  {       // Spanish
    if (pretty == PRETTY_PRINT)
    {
      ptrSin = "<span aria-hidden=\"true\">sen</span><span class=\"hide\"> seno de </span>";
      ptrCos = "<span aria-hidden=\"true\">cos</span><span class=\"hide\"> coseno de </span>";
      ptrACos = "<span aria-hidden=\"true\">arc cos</span><span class=\"hide\"> arco coseno de </span>";
    }
    else if (pretty == TEX)
    {
      ptrSin = "\\sin{";
      ptrCos = "\\cos{";
      ptrACos = "\\arccos";
    }
    else
    {
      ptrSin = "sen";
      ptrCos = "cos";
      ptrACos = "acos";
    }
  }
  else
  {       // English
    if (pretty == PRETTY_PRINT)
    {
      ptrSin = "<span aria-hidden=\"true\">sin</span><span class=\"hide\"> sine of </span>";
      ptrCos = "<span aria-hidden=\"true\">cos</span><span class=\"hide\"> cosine of </span>";
      ptrACos = "<span aria-hidden=\"true\">arccos</span><span class=\"hide\"> arc cosine of </span>";
    }
    else if (pretty == TEX)
    {
      ptrSin = "\\sin{";
      ptrCos = "\\cos{";
      ptrACos = "\\arccos";
    }
    else
    {
      ptrSin = "sin";
      ptrCos = "cos";
      ptrACos = "acos";
    }
  }
  groupLen = groupLength;
  ptrOutput = *pptrOutput;
  if (pretty == PRETTY_PRINT)
  {
    ptrMinus = "&minus;";
    if (lang)
    {
      ptrTimes = "<span class=\"hide\">por</span>";
    }
    else
    {
      ptrTimes = "<span class=\"hide\">times</span>";
    }
    ptrPi = "<span aria-hidden=\"true\">&pi;</span><span class=\"hide\">pi</span>";
    ptrI = "i";
  }
  else if (pretty == TEX)
  {
    ptrMinus = "-";
    ptrTimes = "";  // Do not show multiplication sign in TeX.
    ptrPi = "\\pi ";
    ptrI = "i";
  }
  else
  {
    ptrMinus = "-";
    ptrTimes = "*";
    ptrPi = "Pi";
    ptrI = "I";
  }
  switch (pstFactorInfo->degree)
  {
  case 1:
    LinearEquation(pstFactorInfo->ptrPolyLifted, multiplicity);
    break;
  case 2:
    QuadraticEquation(pstFactorInfo->ptrPolyLifted, multiplicity);
    break;
  case 3:
    CubicEquation(pstFactorInfo->ptrPolyLifted, multiplicity);
    break;
  case 4:
    (void)memcpy(polyIntegerBak, polyInteger, sizeof(polyInteger));
    (void)memcpy(factorInfoIntegerBak, factorInfoInteger, sizeof(factorInfoInteger));
    QuarticEquation(pstFactorInfo->ptrPolyLifted, multiplicity);
    (void)memcpy(polyInteger, polyIntegerBak, sizeof(polyInteger));
    (void)memcpy(factorInfoInteger, factorInfoIntegerBak, sizeof(factorInfoInteger));
    break;
  case 5:
    if (isLinearExponential(pstFactorInfo->ptrPolyLifted, pstFactorInfo->degree,
      multiplicity))
    {          // If polynomial is ax^n+b = 0, show roots.
      break;
    }
    if (isSymmetricOrAlternating(nbrFactor, pstFactorInfo->ptrPolyLifted, pstFactorInfo->degree,
      multiplicity))
    {          // If polynomial is S_n or A_n, indicate that the roots are not solvable.
      break;
    }
    (void)memcpy(polyIntegerBak, polyInteger, sizeof(polyInteger));
    (void)memcpy(factorInfoIntegerBak, factorInfoInteger, sizeof(factorInfoInteger));
    QuinticEquation(pstFactorInfo->ptrPolyLifted, pstFactorInfo->multiplicity);
    (void)memcpy(polyInteger, polyIntegerBak, sizeof(polyInteger));
    (void)memcpy(factorInfoInteger, factorInfoIntegerBak, sizeof(factorInfoInteger));
    break;
  default:
    if (isPalindromic(pstFactorInfo->ptrPolyLifted, pstFactorInfo->degree) &&
        TestCyclotomic(pstFactorInfo->ptrPolyLifted, multiplicity,
        pstFactorInfo->degree))
    {
      break;
    }
    
    if (isLinearExponential(pstFactorInfo->ptrPolyLifted, pstFactorInfo->degree,
      multiplicity))
    {          // If polynomial is ax^n+b = 0, show roots.
      break;
    }
    if (isQuadraticExponential(pstFactorInfo->ptrPolyLifted, pstFactorInfo->degree,
      multiplicity))
    {          // If polynomial is ax^(2n)+bx^n+c = 0, show roots.
      break;
    }
    if (isSymmetricOrAlternating(nbrFactor, pstFactorInfo->ptrPolyLifted, pstFactorInfo->degree,
      multiplicity))
    {          // If polynomial is S_n or A_n, indicate that the roots are not solvable.
      break;
    }
    showX(multiplicity * pstFactorInfo->degree);
    *(ptrOutput - 2) = ':';
    if (lang)
    {
      showText("No puedo determinar si las raíces del polinomio");
      if (nbrFactor >= 0)
      {
        showText(" número ");
        int2dec(&ptrOutput, nbrFactor + 1);
      }
      showText(" se pueden resolver mediante expresiones radicales o no.");
    }
    else
    {
      showText("I cannot determine whether the roots of the polynomial");
      if (nbrFactor >= 0)
      {
        showText(" number ");
        int2dec(&ptrOutput, nbrFactor + 1);
      }
      showText(" can be solved using radical expressions or not.");
    }
    break;
  }
  nbrFactorsFound = nbrFactorsFoundBak;
  *pptrOutput = ptrOutput;
}

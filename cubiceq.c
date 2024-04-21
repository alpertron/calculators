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
#include "rootseq.h"
#include "expression.h"

extern BigInteger Cubic;
extern BigInteger Quadratic;
extern BigInteger Linear;
extern BigInteger Independent;
extern char* ptrOutput;
extern const char* ptrACos;
extern const char* ptrCos;

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

static void showRAndS(void)
{
  BigRationalDivideByInt(&RatDeprIndependent, -2, &Rat3);
  intToBigInteger(&Rat1.numerator, 1);
  intToBigInteger(&Rat1.denominator, 1);
  BigRationalDivideByInt(&RatDiscr, 4, &Rat2);
  MultiplyRationalBySqrtRational(&Rat1, &Rat2);
  for (int ctr = 0; ctr < 2; ctr++)
  {
    if (teach)
    {
      showText("<p>");
    }
    else
    {
      startLine();
    }
    if (pretty == PRETTY_PRINT)
    {
      showText("<var>");
    }
    *ptrOutput = ((ctr == 0) ? 'r' : 's');
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
      showPlusSignOn(ctr == (isCbrtNegative ? 1 : 0),
        TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
      ShowRationalAndSqrParts(&Rat1, &Rat2, 2, ptrTimes);
      endCbrt();
    }
    if (teach)
    {
      showText("</p>");
    }
    else
    {
      endLine();
    }
  }
}

static void showCosT(void)
{
  showText(ptrCos);
  showText(" ");
  showVariable(&ptrOutput, 't');
}

static void showCosTPower3(void)
{
  startParen();
  showCosT();
  endParen();
  showPower(&ptrOutput, 3);
}

static void showArgumentArcCos(void)
{
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
}

static void CasusIrreducibilis(int multiplicity, char currLetter)
{
  enum eSign signRat1;
  if (teach)
  {
    char buf[1000];
    char* ptrBuf;
    showText(lang ? "<p>El discriminante es positivo, lo que indica que las tres raíces son reales. "
      "En este caso no es posible representar las raíces mediante expresiones radicales de números reales, "
      "por lo que se debe utilizar trigonometría.</p>"
      "<p>Comenzando con la fórmula de la triplicación del ángulo:</p><p>4" :
      "<p>The discriminant is positive, which implies that all three roots are real. "
      "In this case the roots cannot be represented by radical expressions of real numbers, "
      "so we must use trigonometry.</p>"
      "<p>Starting with the formula of the triple angle:</p><p>4");
    showText(ptrTimes);
    showCosTPower3();
    showText(" ");
    showText(ptrMinus);
    showText(" 3");
    showText(ptrTimes);
    showCosT();
    showText(" ");
    showText(ptrMinus);
    showText(" ");
    showText(ptrCos);
    startParen();
    showText("3");
    showText(ptrTimes);
    showVariable(&ptrOutput, 't');
    endParen();
    showText(" = 0</p><p>");
    showText(lang ? "Sea " : "Let ");
    showVariable(&ptrOutput, currLetter);
    showText(" = ");
    showVariable(&ptrOutput, 'u');
    showText(" ");
    showText(ptrTimes);
    showCosT();
    showText(lang ? ". De la ecuación anterior a la definición del discriminante:</p><p>":
      ". From the previous equation to the definition of the discriminant:</p><p>");
    showPowerVar(&ptrOutput, 3, 'u');
    showText(ptrTimes);
    showCosTPower3();
    showRatCoeffAndPowerVar(&RatDeprLinear, 1, 'u');
    showText(ptrTimes);
    showText(" ");
    showCosT();
    showPlusMinusRational(&RatDeprIndependent);
    showText(" = 0</p><p>");
    showText(lang ? "El criterio a seguir consiste en igualar los términos de ambas ecuaciones de izquierda a derecha. Dividiendo por " :
      "We will equate the terms of both equations from left to right. Dividing by ");
    showPowerVar(&ptrOutput, 3, 'u');
    showText(" / 4:</p>4");
    showCosTPower3();
    BigRationalMultiplyByInt(&RatDeprLinear, 4, &Rat1);
    showText(ptrTimes);
    ptrBuf = buf;
    showPowerVar(&ptrBuf, 2, 'u');
    *ptrBuf = 0;
    showPlusSignOn(Rat1.numerator.sign == SIGN_POSITIVE, 
      TYPE_PM_SPACE_BEFORE|TYPE_PM_SPACE_AFTER);
    Rat1.numerator.sign = SIGN_POSITIVE;
    showRationalOverStr(&Rat1, buf, ptrTimes);
    showText(ptrTimes);
    showCosT();
    BigRationalMultiplyByInt(&RatDeprIndependent, 4, &Rat1);
    ptrBuf = buf;
    showPowerVar(&ptrBuf, 3, 'u');
    *ptrBuf = 0;
    showPlusSignOn(Rat1.numerator.sign == SIGN_POSITIVE,
      TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
    Rat1.numerator.sign = SIGN_POSITIVE;
    showRationalOverStr(&Rat1, buf, ptrTimes);
    showText(" = 0</p><p>");
    showText(lang ? "El segundo coeficiente debe valer " :
      "The second coefficient must equal ");
    showText(ptrMinus);
    showText(lang ? "3, así que:</p><p>" : "3, so:</p><p>");
    showVariable(&ptrOutput, 'u');
    showText(" = ");
    intToBigInteger(&Rat1.numerator, 2);
    intToBigInteger(&Rat1.denominator, 1);               // 2
    BigRationalDivideByInt(&RatDeprLinear, -3, &Rat2);   // -p/3
    MultiplyRationalBySqrtRational(&Rat1, &Rat2);
    ShowRationalAndSqrParts(&Rat1, &Rat2, 2, ptrTimes);
    showText(lang ? "</p><p>Igualando el último término:</p><p>" : "</p><p>Equating the last term:</p><p>");
    showText(ptrCos);
    startParen();
    showText("3");
    showText(ptrTimes);
    showVariable(&ptrOutput, 't');
    endParen();
    showText(" = ");
    BigRationalMultiplyByInt(&RatDeprIndependent, -4, &Rat1);
    ptrBuf = buf;
    showPowerVar(&ptrBuf, 3, 'u');
    *ptrBuf = 0;
    if (Rat1.numerator.sign == SIGN_NEGATIVE)
    {
      showText(" ");
      showText(ptrMinus);
    }
    showText(" ");
    Rat1.numerator.sign = SIGN_POSITIVE;
    showRationalOverStr(&Rat1, buf, ptrTimes);
    showText(" = ");
    showArgumentArcCos();
    showText("</p>");
  }
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
  showArgumentArcCos();
  endParen();
  if (teach)
  {
    showText("</p>");
  }
  else
  {
    endLine();
  }
  intToBigInteger(&Rat1.numerator, 2);
  intToBigInteger(&Rat1.denominator, 1);               // 2
  BigRationalDivideByInt(&RatDeprLinear, -3, &Rat2);   // -p/3
  MultiplyRationalBySqrtRational(&Rat1, &Rat2);
  signRat1 = Rat1.numerator.sign;
  Rat1.numerator.sign = SIGN_POSITIVE;
  for (int ctr = 0; ctr <= 4; ctr += 2)
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
          showTimesPi(&ptrNumer);
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
          showTimesPi(&ptrNumer);
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
    endShowX();
  }
}

static void showImagPartCubicRoot1(void)
{
  if (pretty != PARI_GP)
  {                // Show i/2.
    showRatConstants("i", "2");
  }
  else
  {
    showText("I/2");
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
}

static void showCompleteCbrt1(const char* ptrPlusMinus)
{
  showText(ptrMinus);
  showText(" ");
  showRatConstants("1", "2");
  showText(" ");
  showText(ptrPlusMinus);
  showText(" ");
  showImagPartCubicRoot1();
}

static void showNonRealCubeRootsOf1(void)
{
  for (char letter = 'e'; letter <= 'f'; letter++)
  {
    if (letter == 'f')
    {
      showText(", ");
    }
    showVariable(&ptrOutput, letter);
    showText(" = ");
    showCompleteCbrt1((letter == 'e') ? "+" : ptrMinus);
  }
}

// The first root x1 is q ^ (1/3)
// The second root x2 is (-1/2) * x1 + i * (3 ^ (1/2) * x1) / 2
// The third root x3 is (-1/2) * x1 - i * (3 ^ (1/2) * x1) / 2
static void bothQuadraticAndLinearCoeffZero(int multiplicity)
{
  bool isCbrtNegative;
  BigRationalNegate(&RatDeprIndependent, &RatDeprIndependent);
  ForceDenominatorPositive(&RatDeprIndependent);
  if (teach)
  {
    showText(lang ? "<p>Las soluciones son la raíz cúbica real de " :
      "<p>The solutions are the real cube root of ");
    showRationalNoParen(&RatDeprIndependent);
    showText(lang ? " y su producto por las dos raíces cúbicas no reales de 1, que son:</p><p>":
      " and the multiplication by both non-real cube roots of 1:</p><p>");
    showNonRealCubeRootsOf1();
    showText("</p>");
  }
  isCbrtNegative = (RatDeprIndependent.numerator.sign == SIGN_NEGATIVE);
  RatDeprIndependent.numerator.sign = SIGN_POSITIVE;
  for (int ctr = 0; ctr < 3; ctr++)
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
      showImagPartCubicRoot1();
      *ptrOutput = ' ';
      ptrOutput++;
      showText(ptrTimes);
      CbrtIndep();     // q^(1/3)
    }
    endShowX();
  }
}

static void showRPlusS(void)
{
  startParen();
  showVariable(&ptrOutput, 'r');
  showText(" + ");
  showVariable(&ptrOutput, 's');
  endParen();
}

static void showN(void)
{
  showCoeffBeforeParen(&RatDeprLinear);
  showRPlusS();
  showPlusMinusRational(&RatDeprIndependent);
  showText(" = 0</p><p>");
}

static void showCardanoMethod(char currLetter)
{
  showText(lang ? "<p>Usando el método de Cardano, asignando " : "<p>Using Cardano's method, setting ");
  showVariable(&ptrOutput, currLetter);
  showText(" = ");
  showVariable(&ptrOutput, 'r');
  showText(" + ");
  showVariable(&ptrOutput, 's');
  showText(":</p><p>");
  // Show (r+s)^3 + N = 0, where N = p(r+s) + q
  showRPlusS();
  showPower(&ptrOutput, 3);
  showN();
  // Show r^3 + 3r^2*s + 3rs^2 + s^3 + N = 0
  showPowerVar(&ptrOutput, 3, 'r');
  intToBigInteger(&Rat1.numerator, 3);
  intToBigInteger(&Rat1.denominator, 1);
  showRatCoeffAndPowerVar(&Rat1, 2, 'r');
  showText(ptrTimes);
  showVariable(&ptrOutput, 's');
  showRatCoeffAndPowerVar(&Rat1, 1, 'r');
  showText(ptrTimes);
  showPowerVar(&ptrOutput, 2, 's');
  showText(" + ");
  showPowerVar(&ptrOutput, 3, 's');
  showN();
  // Show r^3 + s^3 + 3rs(r+s) + N = 0
  showPowerVar(&ptrOutput, 3, 'r');
  showText(" + ");
  showPowerVar(&ptrOutput, 3, 's');
  showRatCoeffAndPowerVar(&Rat1, 1, 'r');
  showText(ptrTimes);
  showVariable(&ptrOutput, 's');
  showText(ptrTimes);
  showRPlusS();
  showN();
  // Show r^3 + s^3 + (3rs + p)(r+s) + q = 0
  showPowerVar(&ptrOutput, 3, 'r');
  showText(" + ");
  showPowerVar(&ptrOutput, 3, 's');
  showText(" + ");
  startParen();
  showRational(&Rat1);
  showText(ptrTimes);
  showVariable(&ptrOutput, 'r');
  showText(ptrTimes);
  showVariable(&ptrOutput, 's');
  showPlusMinusRational(&RatDeprLinear);
  endParen();
  showText(ptrTimes);
  showRPlusS();
  showPlusMinusRational(&RatDeprIndependent);
  showText(" = 0&nbsp;&nbsp;&nbsp;&nbsp;(1)</p><p>");
  showText(lang ? "Como hay una variable adicional, se puede imponer una condición extra. En este caso la elección es:</p><p>" :
    "Since there is an extra variable, we can impose an additional condition. In our case it is:</p><p>");
  // Show 3rs + p = 0
  showRational(&Rat1);
  showText(ptrTimes);
  showVariable(&ptrOutput, 'r');
  showText(ptrTimes);
  showVariable(&ptrOutput, 's');
  showPlusMinusRational(&RatDeprLinear);
  showText(" = 0&nbsp;&nbsp;&nbsp;&nbsp;(2)</p><p>");
  // Show rs = -p/3
  showVariable(&ptrOutput, 'r');
  showText(ptrTimes);
  showVariable(&ptrOutput, 's');
  showText(" = ");
  BigRationalDivideByInt(&RatDeprLinear, -3, &Rat1);
  ForceDenominatorPositive(&Rat1);
  showRational(&Rat1);
  showText("&nbsp;&nbsp;&nbsp;&nbsp;(3)</p>");
  showPowerVar(&ptrOutput, 3, 'r');
  showText(ptrTimes);
  showPowerVar(&ptrOutput, 3, 's');
  showText(" = ");
  BigRationalMultiply(&Rat1, &Rat1, &Rat2);
  BigRationalMultiply(&Rat2, &Rat1, &Rat2);
  showRational(&Rat2);
  showText("&nbsp;&nbsp;&nbsp;&nbsp;(4)</p>");
  showText(lang ? "<p>De (1) y (2):</p><p>" : "<p>From (1) and (2):</p><p>");
  // Show r^3 + s^3 + q = 0
  showPowerVar(&ptrOutput, 3, 'r');
  showText(" + ");
  showPowerVar(&ptrOutput, 3, 's');
  showPlusMinusRational(&RatDeprIndependent);
  showText(" = 0&nbsp;&nbsp;&nbsp;&nbsp;(5)</p><p>");
  showText(lang ? "Multiplicando por " : "Multiplying by ");
  showPowerVar(&ptrOutput, 3, 'r');
  showText(":</p></p>");
  // Show r^6 + r^3*s^3 + q*r^3 = 0
  showPowerVar(&ptrOutput, 6, 'r');
  showText(" + ");
  showPowerVar(&ptrOutput, 3, 'r');
  showText(ptrTimes);
  showPowerVar(&ptrOutput, 3, 's');
  showRatCoeffAndPowerVar(&RatDeprIndependent, 3, 'r');
  showText(" = 0</p><p>");
  showText(lang ? "De (4):</p>" : "From (4):</p>");
  showPowerVar(&ptrOutput, 6, 'r');
  showPlusMinusRational(&Rat2);
  showRatCoeffAndPowerVar(&RatDeprIndependent, 3, 'r');
  showText(" = 0</p><p>");
  showPowerVar(&ptrOutput, 6, 'r');
  showRatCoeffAndPowerVar(&RatDeprIndependent, 3, 'r');
  showPlusMinusRational(&Rat2);
  showText(" = 0</p><p>");
  showText(lang ? "Esta es una ecuación cuadrática en " :
    "This is a quadratic equation in ");
  showPowerVar(&ptrOutput, 3, 'r');
  showText(lang ? ". Al multiplicar (5) por " : ". If we multiplied (5) by ");
  showPowerVar(&ptrOutput, 3, 's');
  showText(lang ? " en vez de " : " instead of ");
  showPowerVar(&ptrOutput, 3, 'r');
  showText(lang ? ", los coeficientes de la ecuación serían los mismos, así que la ecuación cuadrática también permite obtener " :
    ", the equation coefficients would be the same, so the quadratic equation can also give us the value of ");
  showPowerVar(&ptrOutput, 3, 's');
  showText(lang ? ". Sea " : ". Let ");
  showVariable(&ptrOutput, 'w');
  showText(" = ");
  showPowerVar(&ptrOutput, 3, 'r');
  showText(lang ? " o " : " or ");
  showPowerVar(&ptrOutput, 3, 's');
  showText(".</p><p>");
  // Back-up linear and independent coefficients.
  CopyBigInt(&RatQuartic.numerator, &RatLinear.numerator);
  CopyBigInt(&RatQuartic.denominator, &RatLinear.denominator);
  CopyBigInt(&RatCubic.numerator, &RatIndependent.numerator);
  CopyBigInt(&RatCubic.denominator, &RatIndependent.denominator);
  // Populate linear and independent coefficients.
  CopyBigInt(&RatLinear.numerator, &RatDeprIndependent.numerator);
  CopyBigInt(&RatLinear.denominator, &RatDeprIndependent.denominator);
  CopyBigInt(&RatIndependent.numerator, &Rat2.numerator);
  CopyBigInt(&RatIndependent.denominator, &Rat2.denominator);
  showRatCoeffAndPowerVar(NULL, -2, 'w');
  showRatCoeffAndPowerVar(&RatLinear, 1, 'w');
  showRatCoeffAndPowerVar(&RatIndependent, 0, 'w');
  showText(" = 0</p>");
  stepsForQuadraticEquation('w', 'z');
  // Restore linear and independent coefficients.
  CopyBigInt(&RatLinear.numerator, &RatQuartic.numerator);
  CopyBigInt(&RatLinear.denominator, &RatQuartic.denominator);
  CopyBigInt(&RatIndependent.numerator, &RatCubic.numerator);
  CopyBigInt(&RatIndependent.denominator, &RatCubic.denominator);
}

static void linearCoeffNotZero(int multiplicity, char currLetter)
{
  if (teach)
  {
    showText(lang ? "<p>El discriminante es negativo, lo que indica que hay una raíz real y dos raíces complejas conjugadas.</p>" :
      "<p>The discriminant is negative, so there is a real root and two complex conjugate roots.</p>");
    showCardanoMethod(currLetter);
  }
  showRAndS();
  if (teach)
  {
    showText(lang ? "<p>Una raiz cúbica tiene tres soluciones en el campo complejo. Sin embargo, no se puede elegir cualquier valor para " :
      "<p>A cubic root has three solutions in the complex field. But we cannot select any value for ");
    showVariable(&ptrOutput, 'r');
    showText(lang ? " y " : " and ");
    showVariable(&ptrOutput, 's');
    showText(lang ? " porque se debe cumplir la condición (3). Así que el producto debe ser un número real.</p><p>Sea " :
      " because the condition (3) must be true. So their product must be a real number.</p><p>Let ");
    for (char letter = 'r'; letter <= 's'; letter++)
    {
      if (letter == 's')
      {
        showText(", ");
      }
      showVarIndex(letter, 1);
      showText(" = ");
      showVariable(&ptrOutput, letter);
      showText(", ");
      showVarIndex(letter, 2);
      showText(" = ");
      showVariable(&ptrOutput, letter);
      showText(ptrTimes);
      showVariable(&ptrOutput, 'e');
      showText(", ");
      showVarIndex(letter, 3);
      showText(" = ");
      showVariable(&ptrOutput, letter);
      showText(ptrTimes);
      showVariable(&ptrOutput, 'f');
    }
    showText(lang ? ", donde</p><p>" : ", where</p><p>");
    showNonRealCubeRootsOf1();
    showText(lang ? "</p><p>son las raíces cúbicas no reales de 1. Como " :
      "</p><p>are the non-real cubic roots of 1. Since ");
    showVariable(&ptrOutput, 'e');
    showText(ptrTimes);
    showVariable(&ptrOutput, 'f');
    showText(lang ? " es real, pero " : " is real, but ");
    showPowerVar(&ptrOutput, 2, 'e');
    showText(lang ? " y " : " and ");
    showPowerVar(&ptrOutput, 2, 'f');
    showText(lang ? " no lo son, los valores de " : " are not, the values of ");
    showVariable(&ptrOutput, currLetter);
    showText(lang ? " que cumplen la condición (3) son:</p>" : " that follow the condition (3) are:</p>");
    for (int rootNbr = 1; rootNbr <= 3; rootNbr++)
    {
      showText("<p>");
      if (currLetter == 'y')
      {
        showVarIndex('y', rootNbr);
        showText(" = ");
      }
      showVarIndex('x', indexRoot + (rootNbr - 1) * multiplicity);
      if (!BigIntIsZero(&RatQuadratic.numerator))
      {
        BigIntChSign(&RatQuadratic.numerator);
        showPlusMinusRational(&RatQuadratic);
        BigIntChSign(&RatQuadratic.numerator);
      }
      showText(" = ");
      switch (rootNbr)
      {
      case 1:
        showVarIndex('r', 1);
        showText(" + ");
        showVarIndex('s', 1);
        showText(" = ");
        showVariable(&ptrOutput, 'r');
        showText(" + ");
        showVariable(&ptrOutput, 's');
        break;
      default:
        showVarIndex('r', rootNbr);
        showText(" + ");
        showVarIndex('s', 5 - rootNbr);
        showText(" = ");
        for (char letter = 'r'; letter <= 's'; letter++)
        {
          showVariable(&ptrOutput, letter);
          showText(ptrTimes);
          *ptrOutput = ' ';
          ptrOutput++;
          startParen();
          if ((rootNbr == 2) == (letter == 'r'))
          {
            showCompleteCbrt1("+");
          }
          else
          {
            showCompleteCbrt1(ptrMinus);
          }
          endParen();
          if (letter == 'r')
          {
            showText(" + ");
          }
        }
        break;
      }
      showText("</p>");
    }


  }
  for (int ctr = 0; ctr < 3; ctr++)
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
        showRatConstants((pretty == TEX) ? "r - s" :
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
    endShowX();
  }
}

void showCoeffBeforeParen(BigRational* rat)
{
  if (BigIntEqual(&rat->numerator, &rat->denominator))
  {
    showText(" + ");
  }
  else
  {
    bool ratIsMinusOne;
    BigIntChSign(&rat->numerator);
    ratIsMinusOne = BigIntEqual(&rat->numerator, &rat->denominator);
    BigIntChSign(&rat->numerator);
    if (ratIsMinusOne)
    {
      showText(" ");
      showText(ptrMinus);
      showText(" ");
    }
    else
    {
      showPlusMinusRational(rat);
      showText(ptrTimes);
    }
  }
}

void CubicEquation(const int* polynomial, int multiplicity)
{
  const int* ptrPolynomial = polynomial;
  char currLetter = 'x';
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
  // Reduce numerators and denominators.
  BigRationalMultiplyByInt(&RatQuadratic, 1, &RatQuadratic);
  BigRationalMultiplyByInt(&RatLinear, 1, &RatLinear);
  BigRationalMultiplyByInt(&RatIndependent, 1, &RatIndependent);
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
  BigRationalMultiply(&RatDeprLinear, &RatDeprLinear, &Rat1); // p^2
  BigRationalMultiply(&RatDeprLinear, &Rat1, &Rat1);          // p^3
  BigRationalMultiplyByInt(&Rat1, 4, &Rat1);                   // 4p^3
  BigRationalDivideByInt(&Rat1, 27, &Rat1);                  // 4p^3/27
  BigRationalAdd(&RatDiscr, &Rat1, &RatDiscr);                   // q^2 + 4p^3/27
  if (teach && !BigIntIsZero(&RatDeprLinear.numerator))
  {
    const char* ptrDelta;
    if (pretty == PRETTY_PRINT)
    {
      ptrDelta = "&Delta;";
    }
    else if (pretty == TEX)
    {
      ptrDelta = "\\Delta";
    }
    else
    {
      ptrDelta = "D";
    }
    if (!BigIntIsOne(&Cubic))
    {
      showText(lang ? "<p>Dividiendo la ecuación por el coeficiente cúbico:</p><p>" :
        "<p>Dividing the equation by the cubic coefficient:</p><p>");
      showRatCoeffAndPowerVar(NULL, -3, 'x');
      showRatCoeffAndPowerVar(&RatQuadratic, 2, 'x');
      showRatCoeffAndPowerVar(&RatLinear, 1, 'x');
      showRatCoeffAndPowerVar(&RatIndependent, 0, 'x');
      showText(" = 0</p>");
    }
    if (!BigIntIsZero(&Quadratic))
    {
      showText(lang ? "<p>Para eliminar el término cuadrático se debe hacer la sustitución:</p><p>" :
        "<p>To eliminate the quadratic term, we will perform the following substitution:</p><p>");
      showVariable(&ptrOutput, 'x');
      showText(" = ");
      showVariable(&ptrOutput, 'y');
      BigRationalDivideByInt(&RatQuadratic, -3, &Rat2);
      showPlusMinusRational(&Rat2);
      showText("</p><p>");
      showText(lang ? "El valor indicado en la sustitución es la tercera parte del coeficiente cuadrático.</p>" :
        "The constant value in the substitution equals a third of the quadratic coefficient.</p>");
      currLetter = 'y';
      // Show a(y-k)^3 + b(y-k)^2 + c(y-k) + d = 0.
      showText("<p>");
      startParen();
      showVariable(&ptrOutput, 'y');
      showPlusMinusRational(&Rat2);
      endParen();
      showPower(&ptrOutput, 3);
      // Quadratic term is already non-zero.
      showCoeffBeforeParen(&RatQuadratic);
      startParen();
      showVariable(&ptrOutput, 'y');
      showPlusMinusRational(&Rat2);
      endParen();
      showPower(&ptrOutput, 2);
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
      showText(lang ? "Distribuyendo:</p><p>" : "Expanding brackets:</p><p>");
      // Expand all terms.
      // (y-B/3)^3
      showRatCoeffAndPowerVar(NULL, -3, 'y');
      BigIntChSign(&Rat1.numerator);
      CopyBigInt(&Rat4.numerator, &RatQuadratic.numerator);
      CopyBigInt(&Rat4.denominator, &RatQuadratic.denominator);
      BigIntChSign(&Rat4.numerator);
      showRatCoeffAndPowerVar(&Rat4, 2, 'y');
      BigRationalMultiply(&Rat4, &RatQuadratic, &Rat4);
      BigRationalDivideByInt(&Rat4, -3, &Rat4);
      showRatCoeffAndPowerVar(&Rat4, 1, 'y');
      BigRationalMultiply(&Rat4, &RatQuadratic, &Rat4);
      BigRationalDivideByInt(&Rat4, -9, &Rat4);
      showPlusMinusRational(&Rat4);
      // B(y-B/3)^2
      showRatCoeffAndPowerVar(&RatQuadratic, 2, 'y');
      BigRationalMultiply(&RatQuadratic, &RatQuadratic, &Rat4);
      BigRationalDivideByInt(&Rat4, -3, &Rat4);
      BigRationalMultiplyByInt(&Rat4, 2, &Rat3);
      showRatCoeffAndPowerVar(&Rat3, 1, 'y');
      BigRationalMultiply(&Rat4, &RatQuadratic, &Rat4);
      BigRationalDivideByInt(&Rat4, -3, &Rat4);
      showPlusMinusRational(&Rat4);
      // C(y-B/3)
      showRatCoeffAndPowerVar(&RatLinear, 1, 'y');
      BigRationalMultiply(&RatLinear, &RatQuadratic, &Rat4);
      BigRationalDivideByInt(&Rat4, -3, &Rat4);
      showPlusMinusRational(&Rat4);
      // D
      showPlusMinusRational(&RatIndependent);
      showText(" = 0</p><p>");
      showText(lang ? "Simplificando:</p><p>" : "Simplifying:</p><p>");
      // Show y^3 + Py + Q = 0.
      showRatCoeffAndPowerVar(NULL, -3, 'y');
      showRatCoeffAndPowerVar(&RatDeprLinear, 1, 'y');
      showPlusMinusRational(&RatDeprIndependent);
      showText(" = 0</p>");
    }
    showText(lang ? "<p>La naturaleza de las raíces está dada por el valor del discriminante.</p><p>" :
      "<p>The nature of the roots depends on the value of the discriminant.</p><p>");
    showText(ptrDelta);
    showText(" = ");
    intToBigInteger(&Rat1.numerator, -4);
    intToBigInteger(&Rat1.denominator, 1);
    showRatCoeffAndPowerVar(&Rat1, 3, 'p');
    showText(" ");
    intToBigInteger(&Rat1.numerator, -27);
    showRatCoeffAndPowerVar(&Rat1, 2, 'q');
    showText(lang ? "<p>donde " : "<p>where ");
    showVariable(&ptrOutput, 'p');
    showText(lang ? " es el coeficiente lineal y " : " is the linear coefficient and ");
    showVariable(&ptrOutput, 'q');
    showText(lang ? " es el término independiente.</p><p>" : " is the constant term.</p><p>");
    showText(ptrDelta);
    showText(" = ");
    showText(ptrMinus);
    showText("4");
    showText(ptrTimes);
    startParen();
    showRationalNoParen(&RatDeprLinear);
    endParen();
    showPower(&ptrOutput, 3);
    showText(" ");
    showText(ptrMinus);
    showText(" ");
    showText("27");
    showText(ptrTimes);
    startParen();
    showRationalNoParen(&RatDeprIndependent);
    endParen();
    showPower(&ptrOutput, 2);
    showText(" = ");
    BigRationalMultiplyByInt(&RatDiscr, -27, &Rat1);
    ForceDenominatorPositive(&Rat1);
    showRationalNoParen(&Rat1);
    showText("</p>");
  }
  BigRationalDivideByInt(&RatQuadratic, -3, &RatQuadratic);  // -b/3
  ForceDenominatorPositive(&RatQuadratic);
  if (RatDiscr.numerator.sign == RatDiscr.denominator.sign)
  {   // Discriminant is negative.
    if (BigIntIsZero(&RatDeprLinear.numerator))
    {
      bothQuadraticAndLinearCoeffZero(multiplicity);
    }
    else
    {  // Use Cardano's formula.
      linearCoeffNotZero(multiplicity, currLetter);
    }
  }
  else
  {   // Discriminant is positive. Use Viete's formula.
    CasusIrreducibilis(multiplicity, currLetter);
  }
}

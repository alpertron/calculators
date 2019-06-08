/*
This file is part of Alpertron Calculators.

Copyright 2019 Dario Alejandro Alpern

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
#include <string.h>
#include <math.h>
#include "output.h"
#include "polynomial.h"
#define NBR_COEFF 6

static BigInteger Quintic, Quartic, Cubic, Quadratic, Linear, Independent;
static BigInteger discr, commonDenom, tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
static BigRational RatQuartic, RatCubic, RatQuadratic, RatLinear, RatIndependent;
static BigRational RatDeprCubic, RatDeprQuadratic, RatDeprLinear, RatDeprIndependent;
static BigRational RatDiscr, RatDelta0, RatDelta1, RatD;
static BigRational Rat1, Rat2, Rat3, Rat4, Rat5, RatS;
static char teach = 1;
static int indexRoot;

struct stValidateCoeff
{
  char* expression;
  BigInteger* bigint;
  char* textSpanish;
  char* textEnglish;
};

static struct stValidateCoeff astValidateCoeff[NBR_COEFF] =
{
  { NULL, &Quintic, "Coeficiente quíntico <var>a</var>: ", "Quintic coefficient <var>a</var>: " },
  { NULL, &Quartic, "Coeficiente cuártico <var>b</var>: ", "Quartic coefficient <var>b</var>: " },
  { NULL, &Cubic, "Coeficiente cúbico <var>c</var>: ", "Cubic coefficient <var>c</var>: " },
  { NULL, &Quadratic, "Coeficiente cuadrático <var>d</var>: ", "Quadratic coefficient <var>d</var>: " },
  { NULL, &Linear, "Coeficiente lineal <var>e</var>: ", "Linear coefficient <var>e</var>: " },
  { NULL, &Independent, "Término independiente <var>f</var>: ", "Constant term <var>f</var>: " },
};

#define P0 0x0000
#define P1 0x1000
#define P2 0x2000
#define P3 0x3000
#define P4 0x4000
#define P5 0x5000
#define P6 0x6000
#define P7 0x7000

#define Q0 0x0000
#define Q1 0x0100
#define Q2 0x0200
#define Q3 0x0300
#define Q4 0x0400
#define Q5 0x0500
#define Q6 0x0600
#define Q8 0x0800

#define R0 0x0000
#define R1 0x0010
#define R2 0x0020
#define R3 0x0030
#define R4 0x0040
#define R5 0x0050
#define R6 0x0060

#define S0 0x0000
#define S1 0x0001
#define S2 0x0002
#define S3 0x0003
#define S4 0x0004
#define S5 0x0005
#define S6 0x0006

#define END_COEFF 0xFFFF
// Coefficients taken from Dummit's Solving Solvable Quintics article.
static struct stQuinticF20
{
  unsigned short pqrs;
  short coeff;
} astQuinticF20[] =
{
  // Coefficients of independent term
{P0 + Q8 + R0 + S0, 1},      // q^8
{P1 + Q6 + R1 + S0, -13},    // -13pq^6*r
{P5 + Q2 + R2 + S0, 1},      // p^5*q^2*r^2
{P2 + Q4 + R2 + S0, 65},     // 65p^2*q^4*r^2
{P6 + Q0 + R3 + S0, -4},     // -4p^6*r^3
{P3 + Q2 + R3 + S0, -128},   // -128p^3*q^2*r^3
{P0 + Q4 + R3 + S0, 17},     // 17q^4*r^3
{P4 + Q0 + R4 + S0, 48},     // 48p^4*r^4
{P1 + Q2 + R4 + S0, -16},    // -16pq^2*r^4
{P2 + Q0 + R5 + S0, -192},   // -192p^2*r^5
{P0 + Q0 + R6 + S0, 256},    // 256r^6
{P5 + Q3 + R0 + S1, -4},     // -4p^5*q^3*s
{P2 + Q5 + R0 + S1, -12},    // -12p^2*q^5*s
{P6 + Q1 + R1 + S1, 18},     // 18p^6*qrs
{P3 + Q3 + R1 + S1, 12},     // 12p^3*q^3*rs
{P0 + Q5 + R1 + S1, -124},   // -124q^5*r*s
{P4 + Q1 + R2 + S1, 196},    // 196p^4*qr^2*s
{P1 + Q3 + R2 + S1, 590},    // 590pq^3*r^2*s
{P2 + Q1 + R3 + S1, -160},   // -160p^2*qr^3*s
{P0 + Q1 + R4 + S1, -1600},  // -1600qr^4*s
{P7 + Q0 + R0 + S2, -27},    // -27p^7*s^2
{P4 + Q2 + R0 + S2, -150},   // -150p^4*q^2*s^2
{P1 + Q4 + R0 + S2, -125},   // -125pq^4*s^2
{P5 + Q0 + R1 + S2, -99},    // -99p^5*rs^2
{P2 + Q2 + R1 + S2, -725},   // -725p^2*q^2*rs^2
{P3 + Q0 + R2 + S2, 1200},   // 1200p^3*r^2*s^2
{P0 + Q2 + R2 + S2, 3250},   // 3250q^2*r^2*s^2
{P1 + Q0 + R3 + S2, -2000},  // -2000pr^3*s^2
{P1 + Q1 + R1 + S3, -1250},  // -1250pqrs^3
{P2 + Q0 + R0 + S4, 3125},   // 3125p^2*s^4
{P0 + Q0 + R1 + S4, -9375},  // -9375rs^4
{END_COEFF, 0},
  // Coefficients of x
{P1 + Q6 + R0 + S0, -2},     // -2pq^6
{P2 + Q4 + R1 + S0, 19},     // 19p^2*q^4*r
{P3 + Q2 + R2 + S0, -51},    // -51p^3*q^2*r^2
{P0 + Q4 + R2 + S0, 3},      // 3q^4*r^2
{P4 + Q0 + R3 + S0, 32},     // 32p^4*r^3
{P1 + Q2 + R3 + S0, 76},     // 76pq^2*r^3
{P2 + Q0 + R4 + S0, -256},   // -256p^2*r^4
{P0 + Q0 + R5 + S0, 512},    // 512r^5
{P3 + Q3 + R0 + S1, -31},    // -31p^3*q^3*s
{P0 + Q5 + R0 + S1, -58},    // -58q^5*s
{P4 + Q1 + R1 + S1, 117},    // 117p^4*qrs
{P1 + Q3 + R1 + S1, 105},    // 105pq^3*rs
{P2 + Q1 + R2 + S1, 260},    // 260p^2*qr^2*s
{P0 + Q1 + R3 + S1, -2400},  // -2400qr^3*s
{P5 + R0 + R0 + S2, -108},   // -108p^5*s^2
{P2 + Q2 + R0 + S2, -325},   // -325p^2*q^2*s^2
{P3 + Q0 + R1 + S2, 525},    // 525p^3*rs^2
{P0 + Q2 + R1 + S2, 2750},   // 2750q^2*rs^2
{P1 + Q0 + R2 + S2, -500},   // -500pr^2*s^2
{P1 + Q1 + R0 + S3, 625},    // 625pqs^3
{P0 + Q0 + R0 + S4, -3125},  // -3125s^4
{END_COEFF, 0},
  // Coefficients of x^2
{P2 + Q4 + R0 + S0, 1},      // p^2*q^4
{P3 + Q2 + R1 + S0, -6},     // -6p^3*q^2*r
{P0 + Q4 + R1 + S0, -8},     // -8q^4*r
{P4 + Q0 + R2 + S0, 9},      // 9p^4*r^2
{P1 + Q2 + R2 + S0, 76},     // 76pq^2*r^2
{P2 + Q0 + R3 + S0, -136},   // -136p^2*r^3
{P0 + Q0 + R4 + S0, 400},    // 400r^4
{P1 + Q3 + R0 + S1, -50},    // -50pq^3*s
{P2 + Q1 + R1 + S1, 90},     // 90p^2*qrs
{P0 + Q1 + R2 + S1, -1400},  // -1400qr^2*s
{P0 + Q2 + R0 + S2, 625},    // 625q^2*s^2
{P1 + Q0 + R1 + S2, 500},    // 500prs^2
{END_COEFF, 0},
  // Coefficients of x^3
{P0 + Q4 + R0 + S0, -2},     // -2q^4
{P1 + Q2 + R1 + S0, 21},     // 21pq^2*r
{P2 + Q0 + R2 + S0, -40},    // -40p^2*r^2
{P0 + Q0 + R3 + S0, 160},    // 160r^3
{P2 + Q1 + R0 + S1, -15},    // -15p^2*qs
{P0 + Q1 + R1 + S1, -400},   // -400qrs
{P1 + Q0 + R0 + S2, 125},     // 125ps^2
{END_COEFF, 0},
  // Coefficients of x^4
{P1 + Q2 + R0 + S0, 2},      // 2pq^2
{P2 + Q0 + R1 + S0, -6},     // -6p^2*r
{P0 + Q0 + R2 + S0, 40},     // 40r^2
{P0 + Q1 + R0 + S1, -50},    // -50qs
{END_COEFF, 0},
  // Coefficients of x^5
{P0 + Q0 + R1 + S0, 8},      // 8r
{END_COEFF, 0},
  // Coefficients of x^6
{P0 + Q0 + R0 + S0, 1},      // 1
{END_COEFF, -1},
};

static void showX(int multiplicity)
{
  int ctr;
  showText("<li>");
  for (ctr = 0; ctr < multiplicity; ctr++)
  {
    showText("<var>x</var><sub>");
    int2dec(&ptrOutput, indexRoot++);
    showText("</sub> = ");
  }
}

// Compute x = -c_0 / c_1
static void LinearEquation(int *ptrPolynomial, int multiplicity)
{
  showX(multiplicity);
  UncompressBigIntegerB(ptrPolynomial, &Independent);
  ptrPolynomial += 1 + numLimbs(ptrPolynomial);
  UncompressBigIntegerB(ptrPolynomial, &Linear);
  // Negate Independent / Linear but make the denominator positive.
  if (Linear.sign == SIGN_POSITIVE)
  {
    BigIntChSign(&Independent);
  }
  else
  {
    BigIntChSign(&Linear);
  }
  shownbr(&Independent);
  if (Linear.nbrLimbs != 1 || Linear.limbs[0].x != 1)
  {
    showText(" / ");
    shownbr(&Linear);
  }
  showText("</li>");
}

// Compute delta = c_1^2 - 4c_0 * c_2
// If delta > 0 -> x = (-c_1 +/- sqrt(delta))/(2*c_2)
// If delta < 0 -> x = (-c_1 +/- i*sqrt(-delta))/(2*c_2)
// Delta cannot be zero.
static void QuadraticEquation(int* ptrPolynomial, int multiplicity)
{
  int ctr;
  enum eSign signDiscr;
  UncompressBigIntegerB(ptrPolynomial, &Independent);
  ptrPolynomial += 1 + numLimbs(ptrPolynomial);
  UncompressBigIntegerB(ptrPolynomial, &Linear);
  ptrPolynomial += 1 + numLimbs(ptrPolynomial);
  UncompressBigIntegerB(ptrPolynomial, &Quadratic);
  // Compute discriminant (delta).
  BigIntMultiply(&Linear, &Linear, &tmp1);
  BigIntMultiply(&Quadratic, &Independent, &tmp2);
  multint(&tmp2, &tmp2, 4);
  BigIntSubt(&tmp1, &tmp2, &discr);
  // Try to simplify the equation.
  BigIntAdd(&Quadratic, &Quadratic, &Quadratic);
  BigIntGcd(&Linear, &Quadratic, &tmp3);
  BigIntGcd(&tmp3, &discr, &tmp4);
  BigIntDivide(&discr, &tmp4, &tmp5);
  BigIntGcd(&tmp4, &tmp5, &tmp3);
  BigIntDivide(&Linear, &tmp3, &Linear);
  BigIntDivide(&Quadratic, &tmp3, &Quadratic);
  BigIntDivide(&discr, &tmp3, &discr);
  BigIntDivide(&discr, &tmp3, &discr);
  signDiscr = discr.sign;
  discr.sign = SIGN_POSITIVE;
  // If discriminant is perfect square, show its square root.
  squareRoot(discr.limbs, tmp4.limbs, discr.nbrLimbs, &tmp4.nbrLimbs);
  BigIntMultiply(&tmp4, &tmp4, &tmp5);
  BigIntSubt(&tmp5, &discr, &tmp5);
  if (Quadratic.sign == SIGN_NEGATIVE)
  {
    BigIntChSign(&Quadratic);
    BigIntChSign(&Linear);
  }
  BigIntChSign(&Linear);
  for (ctr = 0; ctr < 2; ctr++)
  {
    showX(multiplicity);
    if (!BigIntIsZero(&Linear))
    {
      if (Quadratic.nbrLimbs != 1 || Quadratic.limbs[0].x != 1)
      {
        *ptrOutput++ = '(';
      }
      shownbr(&Linear);
    }
    if (ctr == 1)
    {
      showText(" + ");
    }
    else
    {
      showText(" &minus; ");
    }
    if (BigIntIsZero(&tmp5))
    {
      shownbr(&tmp4);
    }
    else
    {
      shownbr(&discr);
      showText("^(1/2)");
    }
    if (signDiscr == SIGN_NEGATIVE)
    {    // Discriminant is negative.
      showText(" &#8290<var>i</var>");
    }
    if (!BigIntIsZero(&Linear))
    {
      if (Quadratic.nbrLimbs != 1 || Quadratic.limbs[0].x != 1)
      {
        showText(") / ");
        shownbr(&Quadratic);
      }
    }
    showText("</li>");
  }
}

static void CubicEquation(int* ptrPolynomial, int multiplicity)
{
  int ctr;
  UncompressBigIntegerB(ptrPolynomial, &Independent);
  ptrPolynomial += 1 + numLimbs(ptrPolynomial);
  UncompressBigIntegerB(ptrPolynomial, &Linear);
  ptrPolynomial += 1 + numLimbs(ptrPolynomial);
  UncompressBigIntegerB(ptrPolynomial, &Quadratic);
  ptrPolynomial += 1 + numLimbs(ptrPolynomial);
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
  BigRationalAdd(&RatIndependent, &Rat1, &RatDeprIndependent); // q
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
      // The roots are:
      // x1 = q ^ (1/3)
      // x2 = (-1/2) * q ^ (1/3) + i * (3 ^ (1/2) * q ^ (1/3)) / 2
      // x3 = (-1/2) * q ^ (1/3) - i * (3 ^ (1/2) * q ^ (1/3)) / 2
      BigRationalNegate(&RatDeprIndependent, &RatDeprIndependent);
      ForceDenominatorPositive(&RatDeprIndependent);
      for (ctr = 0; ctr < 3; ctr++)
      {
        showX(multiplicity);
        if (ctr == 0)
        {
          showRational(&RatDeprIndependent);
          showText("^(1/3)");
        }
        else
        {
          showText("(-1/2) &#8290;");
          showRational(&RatDeprIndependent);
          showText("^(1/3)");
          if (ctr == 1)
          {
            showText(" +");
          }
          else
          {
            showText(" &minus;");
          }
          showText(" <var>i</var>&#8290;(3^(1/2) &#8290;");
          showRational(&RatDeprIndependent);
          showText("^(1/3)) / 2");
        }
        showText("</li>");
      }
    }
    else
    {
      BigRationalDivideByInt(&RatDeprIndependent, -2, &Rat3);
      intToBigInteger(&Rat1.numerator, 1);
      intToBigInteger(&Rat1.denominator, 2);
      BigRationalDivideByInt(&RatDiscr, 4, &Rat2);
      MultiplyRationalBySqrtRational(&Rat1, &Rat2);
      for (ctr = 0; ctr < 2; ctr++)
      {
        showText("<li><var>");
        *ptrOutput++ = (ctr == 0 ? 'r' : 's');
        showText("</var> = (");
        if (Rat2.numerator.nbrLimbs == 1 && Rat2.numerator.limbs[0].x == 1 &&
          Rat2.denominator.nbrLimbs == 1 && Rat2.denominator.limbs[0].x == 1)
        {
          if (ctr == 0)
          {
            BigRationalAdd(&Rat3, &Rat1, &Rat4);
          }
          else
          {
            BigRationalSubt(&Rat3, &Rat1, &Rat4);
          }
          ForceDenominatorPositive(&Rat4);
          showRational(&Rat4);
        }
        else
        {
          ForceDenominatorPositive(&Rat3);
          showRational(&Rat3);
          if (ctr == 0)
          {
            showText(" + ");
          }
          else
          {
            showText(" &minus; ");
          }
          ShowRationalAndSqrParts(&Rat1, &Rat2, 2);
        }
        showText(")^(1/3)</li>");
      }
      for (ctr = 0; ctr < 3; ctr++)
      {
        showX(multiplicity);
        if (!BigIntIsZero(&Quadratic))
        {
          showRational(&RatQuadratic);
          showText(" + ");
        }
        if (ctr == 0)
        {
          showText("<var>r</var> + <var>s</var>");
        }
        else
        {
          showText("(<var>r</var> + <var>s</var>) / 2");
          if (ctr == 1)
          {
            showText(" + ");
          }
          else
          {
            showText(" &minus; ");
          }
          showText("<var>i</var>&#8290; (<var>r</var> - <var>s</var>)&#8290; 3^(1/2) / 2");
        }
        showText("</li>");
      }
    }
  }
  else
  {   // Discriminant is positive. Use Viete's formula.
    enum eSign signRat1;
    showText("<p><var>&theta;</var> = (1/3) arccos (");
    BigRationalDivide(&RatDeprIndependent, &RatDeprLinear, &Rat1); // q/p
    BigRationalMultiplyByInt(&Rat1, 3, &Rat1);                     // 3q/p
    BigRationalDivideByInt(&Rat1, 2, &Rat1);                       // 3q/(2p)
    ForceDenominatorPositive(&Rat1);
    CopyBigInt(&Rat2.numerator, &RatDeprLinear.denominator);
    CopyBigInt(&Rat2.denominator, &RatDeprLinear.numerator);       // 1/p
    BigRationalMultiplyByInt(&Rat2, -3, &Rat2);                    // -3/p
    ForceDenominatorPositive(&Rat2);
    MultiplyRationalBySqrtRational(&Rat1, &Rat2);
    ShowRationalAndSqrParts(&Rat1, &Rat2, 2);
    showText(")</p><ul>");
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
        showRational(&RatQuadratic);
        if (signRat1 == SIGN_POSITIVE)
        {
          showText(" + ");
        }
      }
      if (signRat1 == SIGN_NEGATIVE)
      {
        showText(" &minus; ");
      }
      ShowRationalAndSqrParts(&Rat1, &Rat2, 2);
      showText("&#8290; cos(<var>&theta;</var>");
      if (ctr == 2)
      {
        showText(" + 2&#8290;&pi; / 3");
      }
      if (ctr == 4)
      {
        showText(" + 4&#8290;&pi; / 3");
      }
      showText(")");
    }
    showText("</li>");
  }
}

static void showFirstTermQuarticEq(int ctr)
{
  if (!BigIntIsZero(&RatCubic.numerator))
  {
    showRational(&RatCubic);
    if (ctr == 0 || ctr == 1)
    {
      showText(" +");
    }
    else
    {
      showText(" &minus;");
    }
  }
  else
  {
    if (ctr == 2 || ctr == 3)
    {
      showText("&minus;");
    }
  }
}

static void biquadraticEquation(int multiplicity)
{
  int isSquareRoot1, isSquareRoot2, ctr;
  // For biquadratic equation, the depressed independent term equals the independent term (r = e).

  // Solutions are x = +/- sqrt((sqrt(r)-p/2)/2) +/- sqrt((-sqrt(r)-p/2)/2)
  // Test whether e is a perfect square.
  BigRationalDivideByInt(&RatDeprQuadratic, 2, &RatDeprQuadratic);
  if (BigRationalSquareRoot(&RatIndependent, &Rat3))
  {           // e is a perfect square. Rat3 = sqrt(r).
    for (ctr = 0; ctr < 4; ctr++)
    {
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
        showText("<var>i</var> &#8290; (");
      }
      if (isSquareRoot1)
      {
        showRational(&Rat4);
      }
      else
      {
        Rat1.numerator.sign = SIGN_POSITIVE;
        showSquareRootOfRational(&Rat1, 2);
      }
      if (ctr == 0 || ctr == 2)
      {
        showText(" + ");
      }
      else
      {
        showText(" &minus; ");
      }
      if (isSquareRoot2)
      {
        showRational(&Rat5);
      }
      else
      {
        enum eSign sign = Rat2.numerator.sign;  // Back up sign.
        Rat2.numerator.sign = SIGN_POSITIVE;
        showSquareRootOfRational(&Rat2, 2);
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
          showText("&#8290;<var>i</var>");
        }
      }
      if (Rat1.numerator.sign == SIGN_NEGATIVE)
      {
        showText(")");
      }
      showText("</li>");
    }
  }
  else
  {      // e is not a perfect square.
         // Compute discriminant as c^2-4e
    BigRationalMultiply(&RatDeprQuadratic, &RatDeprQuadratic, &RatDiscr);
    BigRationalMultiplyByInt(&RatIndependent, 4, &Rat1);
    BigRationalSubt(&RatDiscr, &Rat1, &RatDiscr);
    ForceDenominatorPositive(&RatDiscr);
    intToBigInteger(&Rat1.numerator, 1);
    intToBigInteger(&Rat1.denominator, 2);
    CopyBigInt(&Rat2.numerator, &RatIndependent.numerator);
    CopyBigInt(&Rat2.denominator, &RatIndependent.denominator);
    MultiplyRationalBySqrtRational(&Rat1, &Rat2);
    for (ctr = 0; ctr < 4; ctr++)
    {
      int ctr2;
      int sign = RatDeprQuadratic.numerator.sign;
      RatDeprQuadratic.numerator.sign = SIGN_POSITIVE;
      showX(multiplicity);
      showFirstTermQuarticEq(ctr);
      if (RatDiscr.numerator.sign == SIGN_NEGATIVE)
      {  // x = +/- sqrt((sqrt(e) - c) / 2) +/- i * sqrt((sqrt(e) + c) / 2)
        for (ctr2 = 0; ctr2 < 2; ctr2++)
        {
          if (ctr2 == 1)
          {
            if (ctr == 0 || ctr == 2)
            {
              showText(" + ");
            }
            else
            {
              showText(" &minus; ");
            }
            showText("<var>i</var>&#8290; ");
          }
          if (BigIntIsZero(&RatDeprQuadratic.numerator))
          {
            ShowRationalAndSqrParts(&Rat1, &Rat2, 4);
          }
          else
          {
            showText("(");
            ShowRationalAndSqrParts(&Rat1, &Rat2, 2);
            if (sign == SIGN_POSITIVE)
            {
              showText(ctr2 == 0 ? " &minus; " : " + ");
            }
            else
            {
              showText(ctr2 == 0 ? " + " : " &minus; ");
            }
            showRational(&RatDeprQuadratic);
            showText(")^(1/2)");
          }
        }
      }
      else
      {  // x = sqrt(-c/2 +/- sqrt(c^2 - 4e)/2)
        int isX2Positive;
         // If c >= 0 and positive sqrt => c^2 < c^2 - 4e => e < 0.
         // If c > 0 and negative sqrt = > x^2 is always negative.
         // If c < 0 and positive sqrt = > x^2 is always positive.
         // If c <= 0 and negative sqrt = > c^2 > c^2 - 4e => e > 0
        if (ctr == 0 || ctr == 2)
        {    // Positive sqrt
          isX2Positive = (RatDeprQuadratic.numerator.sign == SIGN_NEGATIVE ||
            RatIndependent.numerator.sign == SIGN_NEGATIVE);
        }
        else
        {    // Negative sqrt
          isX2Positive = ((RatDeprQuadratic.numerator.sign == SIGN_NEGATIVE || 
            BigIntIsZero(&RatDeprQuadratic.numerator)) &&
            RatIndependent.numerator.sign == SIGN_POSITIVE);
        }
        if (!isX2Positive)
        {
          showText("<var>i</var>&#8290; ");
        }
        ForceDenominatorPositive(&RatIndependent);
        if (BigIntIsZero(&RatDeprQuadratic.numerator))
        {
          BigIntChSign(&RatIndependent.numerator);
          showSquareRootOfRational(&RatIndependent, 4);
          BigIntChSign(&RatIndependent.numerator);
          showText("</li>");
        }
        else
        {
          showText("(");
          CopyBigInt(&Rat1.numerator, &RatDeprQuadratic.numerator);
          CopyBigInt(&Rat1.denominator, &RatDeprQuadratic.denominator);
          if (isX2Positive)
          {
            BigIntChSign(&Rat1.numerator);
          }
          showRational(&Rat1);
          if (ctr == 0 || ctr == 2)
          {    // Positive sqrt
            showText(isX2Positive ? " + " : " &minus; ");
          }
          else
          {    // Negative sqrt
            showText(isX2Positive ? " &minus; " : " + ");
          }
          BigRationalMultiply(&RatDeprQuadratic, &RatDeprQuadratic, &Rat1);
          BigRationalSubt(&Rat1, &RatIndependent, &Rat1);
          showSquareRootOfRational(&Rat1, 2);
          showText(")^(1/2)</li>");
        }
      }
    }
  }
}

// Show real or imaginary part of square root of u + q/S
static void showSquareRootOfComplex(char* plus, char* minus)
{
  showText("(");
  // Divide square root part by 4*16.
  BigRationalDivideByInt(&Rat4, 4*16, &Rat4);  
  // Divide rational part by 2*4.
  BigRationalDivideByInt(&Rat3, 2*4, &Rat3);
  showSquareRootOfRational(&Rat4, 2);
  if (Rat3.numerator.sign == SIGN_POSITIVE)
  {
    showText(plus);
    showRational(&Rat3);
  }
  else
  {
    showText(minus);
    BigIntChSign(&Rat3.numerator);
    showRational(&Rat3);
    BigIntChSign(&Rat3.numerator);
  }
  showText(")^(1/2)");
  // Restore square root part.
  BigRationalMultiplyByInt(&Rat4, 4*16, &Rat4);
  // Restore rational part.
  BigRationalMultiplyByInt(&Rat3, 2*4, &Rat3);
}

static void FerrariResolventHasRationalRoot(int multiplicity)
{
  enum eSign sign1;
  int ctr;
  int* ptrValues = factorInfoInteger[0].ptrPolyLifted;
  UncompressBigIntegerB(ptrValues, &RatS.numerator);
  ptrValues += 1 + numLimbs(ptrValues);
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
  sign1 = RatDeprLinear.numerator.sign;
  RatDeprLinear.numerator.sign = SIGN_POSITIVE;
  for (ctr = 0; ctr < 4; ctr++)
  {
    showX(multiplicity);
    showFirstTermQuarticEq(ctr);                         // Show -b/4a and next sign.
    *ptrOutput++ = ' ';
    if (RatS.numerator.sign == SIGN_POSITIVE)
    {             // S is real.
      int isImaginary;
      showSquareRootOfRational(&RatS, 2);
      if (ctr == 0 || ctr == 1)
      {   // Change sign so we always consider q/S.
        BigIntChSign(&Rat1.numerator);
      }
      // Compute k = u^2 - q^2/S^2
      BigRationalMultiply(&RatDeprLinear, &RatDeprLinear, &Rat4);
      BigRationalDivide(&Rat4, &RatS, &Rat4);            // q^2 / S^2
      BigRationalMultiply(&Rat3, &Rat3, &Rat5);          // u^2
      BigRationalSubt(&Rat5, &Rat4, &Rat4);              // k^2 = u^2 - q^2 / S^2
      ForceDenominatorPositive(&Rat4);
      // Determine whether the radicand is positive or not.
      if (Rat1.numerator.sign == SIGN_POSITIVE && Rat3.numerator.sign == SIGN_POSITIVE)
      {
        isImaginary = FALSE;
      }
      else if (Rat1.numerator.sign == SIGN_NEGATIVE && Rat3.numerator.sign == SIGN_NEGATIVE)
      {
        isImaginary = TRUE;
      }
      else if (Rat1.numerator.sign == SIGN_POSITIVE && Rat3.numerator.sign == SIGN_NEGATIVE)
      {
        isImaginary = (Rat5.numerator.sign == SIGN_POSITIVE);
      }
      else
      {
        isImaginary = (Rat5.numerator.sign != SIGN_POSITIVE);
      }
      showText(ctr==0 || ctr==2? " + ": " &minus; ");
      if (isImaginary)
      {
        showText("<var>i</var>&#8290; ");
        BigIntChSign(&Rat1.numerator);
        BigIntChSign(&Rat3.numerator);
      }
      showText("(");
      BigRationalDivideByInt(&Rat1, 4, &Rat1);
      BigRationalDivideByInt(&Rat3, 4, &Rat3);
      ShowRationalAndSqrParts(&Rat1, &Rat2, 2);
      if (Rat3.numerator.sign == SIGN_POSITIVE)
      {
        showText(" + ");
        showRational(&Rat3);
      }
      else
      {
        showText(" &minus; ");
        BigIntChSign(&Rat3.numerator);
        showRational(&Rat3);
        BigIntChSign(&Rat3.numerator);
      }
      if (isImaginary)
      {
        showText(")");
        BigIntChSign(&Rat1.numerator);
        BigIntChSign(&Rat3.numerator);
      }
      BigRationalMultiplyByInt(&Rat1, 4, &Rat1);
      BigRationalMultiplyByInt(&Rat3, 4, &Rat3);
      if (ctr == 0 || ctr == 1)
      {   // Restore sign of q/S.
        BigIntChSign(&Rat1.numerator);
      }
      showText(")^(1/2)");
    }
    else
    {             // S is imaginary.
      // k^2 = sqrt(u^2 + q^2/S^2)
      // x1,2 = -(b/4a) +/- sqrt((k+u)/2) - i*(-S) +/- sqrt((k-u)/2))
      // x3,4 = -(b/4a) +/- sqrt((k+u)/2) + i*(-S) +/- sqrt((k-u)/2))
      // Get value of k^2.
      BigRationalMultiply(&RatDeprLinear, &RatDeprLinear, &Rat4);
      BigRationalDivide(&Rat4, &RatS, &Rat4);            // q^2 / S^2 (S^2 negative)
      BigRationalMultiply(&Rat3, &Rat3, &Rat5);          // u^2
      BigRationalSubt(&Rat5, &Rat4, &Rat4);              // k^2 = u^2 + q^2 / |S^2|
      showSquareRootOfComplex(" + ", " &minus; ");
      showText(" + <var>i</var>&#8290; (");
      if (ctr == 1 || ctr == 3)
      {
        showText("&minus; ");
      }
      BigIntChSign(&RatS.numerator);
      showSquareRootOfRational(&RatS, 2);
      BigIntChSign(&RatS.numerator);
      showText(" &minus; ");
      showSquareRootOfComplex(" &minus; ", " + ");
      showText(")");
    }
  }
}

static void QuarticEquation(int* ptrPolynomial, int multiplicity)
{
  int ctr, isImaginary;
  int* ptrValues;
  enum eSign sign1, sign2;
  UncompressBigIntegerB(ptrPolynomial, &Independent);
  ptrPolynomial += 1 + numLimbs(ptrPolynomial);
  UncompressBigIntegerB(ptrPolynomial, &Linear);
  ptrPolynomial += 1 + numLimbs(ptrPolynomial);
  UncompressBigIntegerB(ptrPolynomial, &Quadratic);
  ptrPolynomial += 1 + numLimbs(ptrPolynomial);
  UncompressBigIntegerB(ptrPolynomial, &Cubic);
  ptrPolynomial += 1 + numLimbs(ptrPolynomial);
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
  BigRationalAdd(&RatDeprIndependent, &Rat1, &RatDeprIndependent);    // 3b^4 + 16b^2*c - 64bd
  BigRationalDivideByInt(&RatDeprIndependent, 256, &RatDeprIndependent); // (3b^4 + 16b^2*c - 64bd)/256
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
  BigIntMultiply(&tmp3, &Rat2.denominator, &tmp3);
  BigIntMultiply(&tmp3, &Rat3.denominator, &tmp3);
  // Quadratic coefficient.
  BigIntMultiply(&Rat1.numerator, &Rat2.denominator, &tmp2);
  BigIntMultiply(&tmp2, &Rat3.denominator, &tmp2);
  // Linear coefficient.
  BigIntMultiply(&Rat2.numerator, &Rat1.denominator, &tmp1);
  BigIntMultiply(&tmp1, &Rat3.denominator, &tmp1);
  // Independent term.
  BigIntMultiply(&Rat3.numerator, &Rat1.denominator, &tmp0);
  BigIntMultiply(&tmp0, &Rat2.denominator, &tmp0);
  values[0] = 3;
  ptrValues = &values[1];
  NumberLength = tmp0.nbrLimbs;
  CompressBigInteger(ptrValues, &tmp0);
  ptrValues += 1 + numLimbs(ptrValues);
  NumberLength = tmp1.nbrLimbs;
  CompressBigInteger(ptrValues, &tmp1);
  ptrValues += 1 + numLimbs(ptrValues);
  NumberLength = tmp2.nbrLimbs;
  CompressBigInteger(ptrValues, &tmp2);
  ptrValues += 1 + numLimbs(ptrValues);
  NumberLength = tmp3.nbrLimbs;
  CompressBigInteger(ptrValues, &tmp3);
  FactorPolyOverIntegers();
  nbrFactorsFound = 0;
  if (factorInfoInteger[0].degree == 1)
  {   // Rational root find. Get root.
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
    showText("<li><var>Q</var> = (");
    if (Rat2.numerator.nbrLimbs == 1 && Rat2.numerator.limbs[0].x == 1)
    {
      BigRationalAdd(&RatDelta1, &Rat2, &Rat1);
      showRational(&Rat1);
    }
    else
    {
      showRational(&RatDelta1);
      showText(" + ");
      ShowRationalAndSqrParts(&Rat1, &Rat2, 2);
    }
    showText(")^(1/3)</li><li><var>S</var> = (1/2)(");
    intToBigInteger(&Rat1.denominator, -2);
    intToBigInteger(&Rat1.numerator, 3);
    BigRationalMultiply(&RatDeprQuadratic, &Rat1, &Rat1);
    showText("<var>Q</var> / 3 ");
    BigRationalDivideByInt(&RatDelta0, 3, &Rat1);
    ForceDenominatorPositive(&Rat1);
    if (Rat1.numerator.sign == SIGN_POSITIVE)
    {
      showText("+ ");
    }
    else
    {
      showText("&minus; ");
    }
    showRational(&Rat1);
    showText("(1/<var>Q</var>))^(1/2)</li>");
    sign1 = RatDeprLinear.numerator.sign;
    RatDeprLinear.numerator.sign = SIGN_POSITIVE;
    for (ctr = 0; ctr < 4; ctr++)
    {
      isImaginary = ((ctr == 0 || ctr == 1) == (RatDeprLinear.numerator.sign == SIGN_POSITIVE));
      showX(multiplicity);
      showFirstTermQuarticEq(ctr);
      showText(" <var>S</var> ");
      if (ctr == 0 || ctr == 2)
      {
        showText("+ ");
      }
      else
      {
        showText("&minus; ");
      }
      if (isImaginary)
      {
        showText("<var>i</var>&#8290; ");
      }
      showText("(1/2)&#8290; (");
      if (!isImaginary)
      {
        showText("&minus;");
      }
      showText("4 &#8290;<var>S</var>&sup2; ");
      if (sign1 == SIGN_POSITIVE)
      {
        showText(isImaginary ? " &minus; " : " + ");
      }
      else
      {
        showText(isImaginary ? " + " : " &minus; ");
      }
      showRational(&Rat1);
      if (!BigIntIsZero(&Rat2.numerator))
      {
        if ((ctr == 0 || ctr == 1) == (sign1 == SIGN_NEGATIVE))
        {
          showText(isImaginary ? " &minus; " : " + ");
        }
        else
        {
          showText(isImaginary ? " + " : " &minus; ");
        }
        showRational(&Rat2);
        showText("&#8290;(1/<var>S</var>)");
      }
      showText(")^(1/2)");
    }
  }
  else
  {
    showText("<li><var>t</var> = arccos(");
    BigRationalDivide(&RatDelta1, &RatDelta0, &Rat1);
    BigRationalDivideByInt(&Rat1, 2, &Rat1);
    CopyBigInt(&Rat2.numerator, &RatDelta0.denominator);
    CopyBigInt(&Rat2.denominator, &RatDelta0.numerator);
    MultiplyRationalBySqrtRational(&Rat1, &Rat2);
    ShowRationalAndSqrParts(&Rat1, &Rat2, 2);
    showText(")</li><li><var>S</var> = (1/2)&#8290;(");
    intToBigInteger(&Rat1.numerator, -2);
    intToBigInteger(&Rat1.denominator, 3);
    BigRationalMultiply(&RatDeprQuadratic, &Rat1, &Rat2);
    ForceDenominatorPositive(&Rat2);
    showRational(&Rat2);
    showText(" + ");
    BigIntChSign(&Rat1.numerator);
    CopyBigInt(&Rat2.numerator, &RatDelta0.numerator);
    CopyBigInt(&Rat2.denominator, &RatDelta0.denominator);
    MultiplyRationalBySqrtRational(&Rat1, &Rat2);
    ForceDenominatorPositive(&Rat1);
    ShowRationalAndSqrParts(&Rat1, &Rat2, 2);
    showText(" &#8290; cos(<var>t</var> / 3))^(1/2)</li>");
    BigRationalMultiplyByInt(&RatDeprQuadratic, -2, &Rat1);
    ForceDenominatorPositive(&Rat1);
    sign1 = Rat1.numerator.sign;
    Rat1.numerator.sign = SIGN_POSITIVE;
    CopyBigInt(&Rat2.numerator, &RatDeprLinear.numerator);
    CopyBigInt(&Rat2.denominator, &RatDeprLinear.denominator);
    sign2 = Rat1.numerator.sign;
    Rat2.numerator.sign = SIGN_POSITIVE;
    isImaginary = RatDeprQuadratic.numerator.sign == SIGN_POSITIVE ||
      RatD.numerator.sign == SIGN_POSITIVE;
    for (ctr = 0; ctr < 4; ctr++)
    {
      showX(multiplicity);
      showFirstTermQuarticEq(ctr);
      showText(" <var>S</var> ");
      if (ctr == 0 || ctr == 2)
      {
        showText("+ ");
      }
      else
      {
        showText("&minus; ");
      }
      if (isImaginary)
      {
        showText("<var>i</var>&#8290; ");
      }
      showText("(1/2)&#8290; (");
      if (!isImaginary)
      {
        showText("&minus;");
      }
      showText("4 &#8290;<var>S</var>&sup2; ");
      if (sign1 == SIGN_POSITIVE)
      {
        showText(isImaginary ? " &minus; " : " + ");
      }
      else
      {
        showText(isImaginary ? " + " : " &minus; ");
      }
      showRational(&Rat1);
      if (!BigIntIsZero(&Rat2.numerator))
      {
        if ((ctr == 0 || ctr == 1) == (sign2 == SIGN_NEGATIVE))
        {
          showText(isImaginary? " &minus; ": " + ");
        }
        else
        {
          showText(isImaginary? " + ": " &minus; ");
        }
        showRational(&Rat2);
        showText("&#8290;(1/<var>S</var>)");
      }
      showText(")^(1/2)");
    }
  }
}

static void QuinticEquation(int* ptrPolynomial, int multiplicity)
{
  int ctr;
  int* ptrValues;
  struct stQuinticF20 *pstQuinticF20;
  UncompressBigIntegerB(ptrPolynomial, &Independent);
  ptrPolynomial += 1 + numLimbs(ptrPolynomial);
  UncompressBigIntegerB(ptrPolynomial, &Linear);
  ptrPolynomial += 1 + numLimbs(ptrPolynomial);
  UncompressBigIntegerB(ptrPolynomial, &Quadratic);
  ptrPolynomial += 1 + numLimbs(ptrPolynomial);
  UncompressBigIntegerB(ptrPolynomial, &Cubic);
  ptrPolynomial += 1 + numLimbs(ptrPolynomial);
  UncompressBigIntegerB(ptrPolynomial, &Quartic);
  ptrPolynomial += 1 + numLimbs(ptrPolynomial);
  UncompressBigIntegerB(ptrPolynomial, &Quintic);
  // Get rational coefficients of monic equation.
  CopyBigInt(&RatQuartic.numerator, &Quartic);
  CopyBigInt(&RatQuartic.denominator, &Quintic);
  CopyBigInt(&RatCubic.numerator, &Cubic);
  CopyBigInt(&RatCubic.denominator, &Quintic);
  CopyBigInt(&RatQuadratic.numerator, &Quadratic);
  CopyBigInt(&RatQuadratic.denominator, &Quintic);
  CopyBigInt(&RatLinear.numerator, &Linear);
  CopyBigInt(&RatLinear.denominator, &Quintic);
  CopyBigInt(&RatIndependent.numerator, &Independent);
  CopyBigInt(&RatIndependent.denominator, &Quintic);
  // Compute coefficients of depressed equation x^5 + px^3 + qx^2 + rx + s
  // where: p = (5c - 2b^2)/5, q = (25d - 15bc + 4b^3)/25, r = (125e - 50bd + 15b^2*c - 3b^4)/125,
  // s = (3125f - 625be + 125b^2*e - 25b^3*c + 4b^5)/3125
  BigRationalMultiply(&RatQuartic, &RatCubic, &Rat2);
  BigRationalMultiplyByInt(&Rat2, -15, &RatDeprQuadratic);         // -15bc
  BigRationalMultiply(&RatQuartic, &RatQuadratic, &Rat2);
  BigRationalMultiplyByInt(&Rat2, -50, &RatDeprLinear);            // -50bd
  BigRationalMultiply(&RatQuartic, &RatLinear, &Rat2);
  BigRationalMultiplyByInt(&Rat2, -625, &RatDeprIndependent);      // -625be
  BigRationalMultiply(&RatQuartic, &RatQuartic, &Rat1);            // b^2
  BigRationalMultiplyByInt(&Rat1, -2, &RatDeprCubic);              // -2b^2
  BigRationalMultiply(&Rat1, &RatCubic, &Rat2);
  BigRationalMultiplyByInt(&Rat2, 15, &Rat2);                      // 15b^2*c
  BigRationalAdd(&RatDeprLinear, &Rat2, &RatDeprLinear);
  BigRationalMultiply(&Rat1, &RatLinear, &Rat2);
  BigRationalMultiplyByInt(&Rat2, 125, &Rat2);                     // 125b^2*e
  BigRationalAdd(&RatDeprIndependent, &Rat2, &RatDeprIndependent);
  BigRationalMultiply(&Rat1, &RatQuartic, &Rat1);                  // b^3
  BigRationalMultiplyByInt(&Rat1, 4, &Rat2);                       // 4b^3
  BigRationalAdd(&RatDeprQuadratic, &Rat2, &RatDeprQuadratic);
  BigRationalMultiply(&Rat1, &RatCubic, &Rat2);
  BigRationalMultiplyByInt(&Rat2, -25, &Rat2);                     // -15b^3*c
  BigRationalAdd(&RatDeprIndependent, &Rat2, &RatDeprIndependent);
  BigRationalMultiply(&Rat1, &RatQuartic, &Rat1);                  // b^4
  BigRationalMultiplyByInt(&Rat1, -3, &Rat2);                      // -3b^4
  BigRationalAdd(&RatDeprLinear, &Rat2, &RatDeprLinear);
  BigRationalMultiply(&Rat1, &RatQuartic, &Rat1);                  // b^5
  BigRationalMultiplyByInt(&Rat1, 4, &Rat2);                       // 4b^5
  BigRationalAdd(&RatDeprIndependent, &Rat2, &RatDeprIndependent);
  BigRationalDivideByInt(&RatDeprCubic, 5, &RatDeprCubic);
  BigRationalAdd(&RatDeprCubic, &RatCubic, &RatDeprCubic);
  BigRationalDivideByInt(&RatDeprQuadratic, 25, &RatDeprQuadratic);
  BigRationalAdd(&RatDeprQuadratic, &RatQuadratic, &RatDeprQuadratic);
  BigRationalDivideByInt(&RatDeprLinear, 125, &RatDeprLinear);
  BigRationalAdd(&RatDeprLinear, &RatLinear, &RatDeprLinear);
  BigRationalDivideByInt(&RatDeprIndependent, 3125, &RatDeprIndependent);
  BigRationalAdd(&RatDeprIndependent, &RatIndependent, &RatDeprIndependent);
  
  // Generate sixth degree polynomial F20 as indicated in
  // Dummit's Solving Solvable Quintics article.
  values[0] = 6;          // Degree of polynomial to factor.
  ptrValues = &values[1];
  
  intToBigInteger(&Rat1.numerator, 0);
  intToBigInteger(&Rat1.denominator, 1);
  intToBigInteger(&commonDenom, 1);             // Common denominator.
  intToBigInteger(&tmp0, 0);
  intToBigInteger(&tmp1, 0);
  intToBigInteger(&tmp2, 0);
  intToBigInteger(&tmp3, 0);
  intToBigInteger(&tmp4, 0);
  intToBigInteger(&tmp5, 0);
  for (pstQuinticF20 = astQuinticF20; ; pstQuinticF20++)
  {
    if (pstQuinticF20->pqrs == END_COEFF)
    {
      BigIntGcd(&commonDenom, &Rat1.denominator, &tmp7);
      BigIntMultiply(&tmp1, &Rat1.denominator, &tmp0);
      BigIntDivide(&tmp0, &tmp7, &tmp0);
      BigIntMultiply(&tmp2, &Rat1.denominator, &tmp1);
      BigIntDivide(&tmp1, &tmp7, &tmp1);
      BigIntMultiply(&tmp3, &Rat1.denominator, &tmp2);
      BigIntDivide(&tmp2, &tmp7, &tmp2);
      BigIntMultiply(&tmp4, &Rat1.denominator, &tmp3);
      BigIntDivide(&tmp3, &tmp7, &tmp3);
      BigIntMultiply(&tmp5, &Rat1.denominator, &tmp4);
      BigIntDivide(&tmp4, &tmp7, &tmp4);
      BigIntMultiply(&tmp6, &Rat1.denominator, &tmp5);
      BigIntDivide(&tmp5, &tmp7, &tmp5);
      BigIntMultiply(&Rat1.numerator, &commonDenom, &tmp6);
      BigIntDivide(&tmp6, &tmp7, &tmp6);
      BigIntMultiply(&commonDenom, &Rat1.denominator, &commonDenom);
      BigIntDivide(&commonDenom, &tmp7, &commonDenom);
      intToBigInteger(&Rat1.numerator, 0);
      intToBigInteger(&Rat1.denominator, 1);
      if (pstQuinticF20->coeff < 0)
      {
        break;
      }
      continue;
    }
    intToBigInteger(&Rat2.numerator, pstQuinticF20->coeff);
    intToBigInteger(&Rat2.denominator, 1);
    // Multiply by power of p.
    for (ctr = pstQuinticF20->pqrs & 0xF000; ctr > 0; ctr -= 0x1000)
    {
      BigRationalMultiply(&Rat2, &RatDeprCubic, &Rat2);
    }
    // Multiply by power of q.
    for (ctr = pstQuinticF20->pqrs & 0x0F00; ctr > 0; ctr -= 0x0100)
    {
      BigRationalMultiply(&Rat2, &RatDeprQuadratic, &Rat2);
    }
    // Multiply by power of r.
    for (ctr = pstQuinticF20->pqrs & 0x00F0; ctr > 0; ctr -= 0x0010)
    {
      BigRationalMultiply(&Rat2, &RatDeprLinear, &Rat2);
    }
    // Multiply by power of s.
    for (ctr = pstQuinticF20->pqrs & 0x000F; ctr > 0; ctr--)
    {
      BigRationalMultiply(&Rat2, &RatDeprIndependent, &Rat2);
    }
    BigRationalAdd(&Rat1, &Rat2, &Rat1);
  }
  NumberLength = tmp0.nbrLimbs;
  CompressBigInteger(ptrValues, &tmp0);
  ptrValues += 1 + numLimbs(ptrValues);
  NumberLength = tmp1.nbrLimbs;
  CompressBigInteger(ptrValues, &tmp1);
  ptrValues += 1 + numLimbs(ptrValues);
  NumberLength = tmp2.nbrLimbs;
  CompressBigInteger(ptrValues, &tmp2);
  ptrValues += 1 + numLimbs(ptrValues);
  NumberLength = tmp3.nbrLimbs;
  CompressBigInteger(ptrValues, &tmp3);
  ptrValues += 1 + numLimbs(ptrValues);
  NumberLength = tmp4.nbrLimbs;
  CompressBigInteger(ptrValues, &tmp4);
  ptrValues += 1 + numLimbs(ptrValues);
  NumberLength = tmp5.nbrLimbs;
  CompressBigInteger(ptrValues, &tmp5);
  ptrValues += 1 + numLimbs(ptrValues);
  NumberLength = tmp6.nbrLimbs;
  CompressBigInteger(ptrValues, &tmp6);
  FactorPolyOverIntegers();
  // If there is a linear factor, that means that the quintic can be expressed with radicands.
  if (factorInfoInteger[0].degree == 1)
  {
    showText(lang ? "<p>La ecuación quíntica se puede expresar mediante radicandos.</p>" :
      "The quintic equation can be expressed with radicands.");
  }
  else
  {
    showText(lang ? "<p>La ecuación quíntica no se puede expresar mediante radicandos.</p>" :
      "The quintic equation cannot be expressed with radicands.");
  }
  nbrFactorsFound = 0;   // Exit loop in caller so flow does not reenter here.
}

void rootsEqText(char* coefAText, char* coefBText, char* coefCText,
  char* coefDText, char* coefEText, char* coefFText)
{
  int coeffNbr, nbrFactor;
  int* ptrValues;
  enum eExprErr rc;
  struct stValidateCoeff* pstValidateCoeff = astValidateCoeff;
  struct sFactorInfo* pstFactorInfo;
  astValidateCoeff[0].expression = coefAText;
  astValidateCoeff[1].expression = coefBText;
  astValidateCoeff[2].expression = coefCText;
  astValidateCoeff[3].expression = coefDText;
  astValidateCoeff[4].expression = coefEText;
  astValidateCoeff[5].expression = coefFText;
  ptrOutput = output;
  showText("2<p>");
  for (coeffNbr = 0; coeffNbr < NBR_COEFF; coeffNbr++)
  {
    rc = ComputeExpression(pstValidateCoeff->expression, 1,
      pstValidateCoeff->bigint);
    if (rc != EXPR_OK)
    {
      strcpy(ptrOutput, lang ? pstValidateCoeff->textSpanish : pstValidateCoeff->textEnglish);
      ptrOutput = output + strlen(output);
      textError(ptrOutput, rc);
      ptrOutput = output + strlen(output);
      strcpy(ptrOutput, "</p>");
      ptrOutput = output + strlen(output);
      break;
    }
    pstValidateCoeff++;
  }
  if (coeffNbr == NBR_COEFF)
  {
    // Send to "values" array all coefficients.
    values[0] = 5;
    if (BigIntIsZero(&Quintic))
    {
      values[0] = 4;
      if (BigIntIsZero(&Quartic))
      {
        values[0] = 3;
        if (BigIntIsZero(&Cubic))
        {
          values[0] = 2;
          if (BigIntIsZero(&Quadratic))
          {
            values[0] = 1;
            if (BigIntIsZero(&Linear))
            {
              values[0] = 0;
            }
          }
        }
      }
    }
    ptrValues = &values[1];
    NumberLength = Independent.nbrLimbs;
    CompressBigInteger(ptrValues, &Independent);
    ptrValues += 1 + numLimbs(ptrValues);
    NumberLength = Linear.nbrLimbs;
    CompressBigInteger(ptrValues, &Linear);
    ptrValues += 1 + numLimbs(ptrValues);
    NumberLength = Quadratic.nbrLimbs;
    CompressBigInteger(ptrValues, &Quadratic);
    ptrValues += 1 + numLimbs(ptrValues);
    NumberLength = Cubic.nbrLimbs;
    CompressBigInteger(ptrValues, &Cubic);
    ptrValues += 1 + numLimbs(ptrValues);
    NumberLength = Quartic.nbrLimbs;
    CompressBigInteger(ptrValues, &Quartic);
    ptrValues += 1 + numLimbs(ptrValues);
    NumberLength = Quintic.nbrLimbs;
    CompressBigInteger(ptrValues, &Quintic);
    FactorPolyOverIntegers();
    ptrOutput = output;
    showText("2<p><h2>");
    outputOriginalPolynomial(ptrOutput, groupLen);
    ptrOutput += strlen(ptrOutput);
    showText(" = 0</h2>");
    indexRoot = 1;
    if (nbrFactorsFound > 1)
    {
      showText("<p>");
      pstFactorInfo = factorInfoInteger;
      for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
      {
        *ptrOutput++ = '(';
        outputPolynomialFactor(ptrOutput, groupLen, pstFactorInfo);
        ptrOutput += strlen(ptrOutput);
        *ptrOutput++ = ')';
        if (nbrFactor < nbrFactorsFound - 1)
        {
          showText(" &#8290;");
        }
        pstFactorInfo++;
      }
      showText(" = 0</p>");
    }
    pstFactorInfo = factorInfoInteger;
    for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
    {
      switch (pstFactorInfo->degree)
      {
      case 1:
        LinearEquation(pstFactorInfo->ptrPolyLifted, pstFactorInfo->multiplicity);
        break;
      case 2:
        QuadraticEquation(pstFactorInfo->ptrPolyLifted, pstFactorInfo->multiplicity);
        break;
      case 3:
        CubicEquation(pstFactorInfo->ptrPolyLifted, pstFactorInfo->multiplicity);
        break;
      case 4:
        QuarticEquation(pstFactorInfo->ptrPolyLifted, pstFactorInfo->multiplicity);
        break;
      case 5:
        QuinticEquation(pstFactorInfo->ptrPolyLifted, pstFactorInfo->multiplicity);
        break;
      }
      pstFactorInfo++;
    }
  }
  showText("</ul>");
  showText(lang ? "<p>" COPYRIGHT_SPANISH "</p>" :
    "<p>" COPYRIGHT_ENGLISH "</p>");
}

#ifdef __EMSCRIPTEN__
EXTERNALIZE void doWork(void)
{
  int flags;
  char* ptrData = inputString;
  char* ptrCoeffA, * ptrCoeffB, * ptrCoeffC;
  char* ptrCoeffD, * ptrCoeffE, * ptrCoeffF;
  groupLen = 0;
  while (*ptrData != ',')
  {
    groupLen = groupLen * 10 + (*ptrData++ - '0');
  }
  ptrData++;                    // Skip comma.
  flags = *ptrData;
  lang = flags & 1;
  teach = flags & 2;
  ptrCoeffA = ptrData + 2;  // Skip flags and comma.
  ptrCoeffB = ptrCoeffA + strlen(ptrCoeffA) + 1;
  ptrCoeffC = ptrCoeffB + strlen(ptrCoeffB) + 1;
  ptrCoeffD = ptrCoeffC + strlen(ptrCoeffC) + 1;
  ptrCoeffE = ptrCoeffD + strlen(ptrCoeffD) + 1;
  ptrCoeffF = ptrCoeffE + strlen(ptrCoeffE) + 1;
  rootsEqText(ptrCoeffA, ptrCoeffB, ptrCoeffC, ptrCoeffD, ptrCoeffE, ptrCoeffF);
  databack(output);
}
#endif

//
// This file is part of Alpertron Calculators.
//
// Copyright 2019-2022 Dario Alejandro Alpern
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

// David Dummit's coefficients for solvable quintics
// See http://www.emba.uvm.edu/~ddummit/quintics/quintics.html

#include <string.h>
#include "rootseq.h"
#include "expression.h"
#include "quintics.h"
#include "musl.h"

extern struct monomial arrayB[];
extern struct monomial arrayD[];
extern struct monomial arrayF[];

enum
{
  index_T1 = 0,
  index_T2,
  index_T3,
  index_T4,
  index_l0,
  index_v,
  index_O,
};

struct stQuinticF20
{
  unsigned short pqrs;
  short coeff;
};

static BigRational RatValues[7];  // T1, T2, T3, T4, l0, v and O.
static BigRational RatRoot;       // Rational root of polynomial F20.
static BigRational RatM;
static BigRational RatN;
static BigRational RatR;
extern char* ptrOutput;
extern const char* ptrCos;
extern const char* ptrACos;
extern const char* ptrPi;
extern BigInteger Quintic;
extern BigInteger Quartic;
extern BigInteger Cubic;
extern BigInteger Quadratic;
extern BigInteger Linear;
extern BigInteger Independent;

// Coefficients taken from Dummit's Solving Solvable Quintics article.
static struct stQuinticF20 astQuinticF20[] =
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

// If pretty == TEX, show text as is.
// If pretty == PARI_GP, show all characters except underscores.
static void showExpr(const char *expr)
{
  const char* ptrExpr = expr;
  while (*ptrExpr != '\0')
  {
    if ((pretty == TEX) || (*ptrExpr != '_'))
    {
      *ptrOutput = *ptrExpr;
      ptrOutput++;
    }
    ptrExpr++;
  }
}

// Generate sixth degree polynomial F20 as indicated in
// Dummit's Solving Solvable Quintics article.
static void FactorPolynomialF20(void)
{
  unsigned int ctr;
  int* ptrValues;
  const struct stQuinticF20* pstQuinticF20;
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
  intToBigInteger(&tmp6, 0);
  pstQuinticF20 = astQuinticF20;
  for (;;)
  {
    if (pstQuinticF20->pqrs == END_COEFF)
    {
      BigIntGcd(&commonDenom, &Rat1.denominator, &tmp7);
      (void)BigIntMultiply(&tmp1, &Rat1.denominator, &tmp0);
      (void)BigIntDivide(&tmp0, &tmp7, &tmp0);
      (void)BigIntMultiply(&tmp2, &Rat1.denominator, &tmp1);
      (void)BigIntDivide(&tmp1, &tmp7, &tmp1);
      (void)BigIntMultiply(&tmp3, &Rat1.denominator, &tmp2);
      (void)BigIntDivide(&tmp2, &tmp7, &tmp2);
      (void)BigIntMultiply(&tmp4, &Rat1.denominator, &tmp3);
      (void)BigIntDivide(&tmp3, &tmp7, &tmp3);
      (void)BigIntMultiply(&tmp5, &Rat1.denominator, &tmp4);
      (void)BigIntDivide(&tmp4, &tmp7, &tmp4);
      (void)BigIntMultiply(&tmp6, &Rat1.denominator, &tmp5);
      (void)BigIntDivide(&tmp5, &tmp7, &tmp5);
      (void)BigIntMultiply(&Rat1.numerator, &commonDenom, &tmp6);
      (void)BigIntDivide(&tmp6, &tmp7, &tmp6);
      (void)BigIntMultiply(&commonDenom, &Rat1.denominator, &commonDenom);
      (void)BigIntDivide(&commonDenom, &tmp7, &commonDenom);
      intToBigInteger(&Rat1.numerator, 0);
      intToBigInteger(&Rat1.denominator, 1);
      if (pstQuinticF20->coeff < 0)
      {
        break;
      }
    }
    else
    {
      intToBigInteger(&Rat2.numerator, pstQuinticF20->coeff);
      intToBigInteger(&Rat2.denominator, 1);
      // Multiply by power of p.
      for (ctr = pstQuinticF20->pqrs & 0xF000U; ctr > 0U; ctr -= 0x1000U)
      {
        BigRationalMultiply(&Rat2, &RatDeprCubic, &Rat2);
      }
      // Multiply by power of q.
      for (ctr = pstQuinticF20->pqrs & 0x0F00U; ctr > 0U; ctr -= 0x0100U)
      {
        BigRationalMultiply(&Rat2, &RatDeprQuadratic, &Rat2);
      }
      // Multiply by power of r.
      for (ctr = pstQuinticF20->pqrs & 0x00F0U; ctr > 0U; ctr -= 0x0010U)
      {
        BigRationalMultiply(&Rat2, &RatDeprLinear, &Rat2);
      }
      // Multiply by power of s.
      for (ctr = pstQuinticF20->pqrs & 0x000FU; ctr > 0U; ctr--)
      {
        BigRationalMultiply(&Rat2, &RatDeprIndependent, &Rat2);
      }
      BigRationalAdd(&Rat1, &Rat2, &Rat1);
    }
    pstQuinticF20++;
  }
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
  ptrValues += numLimbs(ptrValues);
  ptrValues++;
  NumberLength = tmp4.nbrLimbs;
  BigInteger2IntArray(ptrValues, &tmp4);
  ptrValues += numLimbs(ptrValues);
  ptrValues++;
  NumberLength = tmp5.nbrLimbs;
  BigInteger2IntArray(ptrValues, &tmp5);
  ptrValues += numLimbs(ptrValues);
  ptrValues++;
  NumberLength = tmp6.nbrLimbs;
  BigInteger2IntArray(ptrValues, &tmp6);
  (void)FactorPolyOverIntegers();
}

// Compute rational value from array of monomials.
static void computeFormula(struct monomial** ppstMonomial, BigRational *rat)
{
  struct monomial* pstMonomial = *ppstMonomial;
  intToBigInteger(&rat->numerator, 0);
  intToBigInteger(&rat->denominator, 1);
  while (pstMonomial->coefficient != 0)
  {
    unsigned int exponents = (unsigned int)pstMonomial->exponents;
    unsigned int ctr = 4U;
    // Find Rat1 = 5^exp5 * cubic^expP * quadr^expQ * linear^expR * const^expS
    // cubic, quadr, etc. are depressed coefficients and they are rational.
    intToBigInteger(&Rat1.numerator, 1);
    intToBigInteger(&Rat1.denominator, 1);
    for (;;)
    {
      BigRationalMultiply(&Rat1, &Rat1, &Rat1);
      if ((exponents & (FIVEexp << 4)) != 0U)
      {
        BigRationalMultiplyByInt(&Rat1, 5, &Rat1);
      }
      if ((exponents & (Pexp << 4)) != 0U)
      {
        BigRationalMultiply(&Rat1, &RatDeprCubic, &Rat1);
      }
      if ((exponents & (Qexp << 4)) != 0U)
      {
        BigRationalMultiply(&Rat1, &RatDeprQuadratic, &Rat1);
      }
      if ((exponents & (Rexp << 4)) != 0U)
      {
        BigRationalMultiply(&Rat1, &RatDeprLinear, &Rat1);
      }
      if ((exponents & (Sexp << 4)) != 0U)
      {
        BigRationalMultiply(&Rat1, &RatDeprIndependent, &Rat1);
      }
      if (ctr == 0U)
      {
        break;
      }
      exponents <<= 1;
      ctr--;
    }
    BigRationalMultiplyByInt(&Rat1, pstMonomial->coefficient, &Rat1);
    BigRationalAdd(rat, &Rat1, rat);
    pstMonomial++;
  }
  *ppstMonomial = pstMonomial + 1;  // Skip terminator.
}

// Compute T1, T2, T3, T4, l0, v and O as follows:
// T1 = (b10 + b11*t + b12*t^2 + b13*t^3 + b14*t^4 + b15*t^5)/(2F)
// T2 = (b20 + b21*t + b22*t^2 + b23*t^3 + b24*t^4 + b25*t^5)/(2DF)
// T3 = (b30 + b31*t + b32*t^2 + b33*t^3 + b34*t^4 + b35*t^5)/(2F)
// T4 = (b40 + b41*t + b42*t^2 + b43*t^3 + b44*t^4 + b45*t^5)/(2DF)
// l0 = (a0 + a1*t + a2*t^2 + a3*t^3 + a4*t^4 + a5*t^5)/F
// v = (c0 + c1*t + c2*t^2 + c3*t^3 + c4*t^4 + c5*t^5)/(2DF)
// O = (o0 + o1*t + o2*t^2 + o3*t^3 + o4*t^4 + o5*t^5)/(DF)
// where t is the rational root of polynomial F20.
static void computeArrayValues(void)
{
  struct BigRational* pRatValues = &RatValues[0];
  struct monomial* pstMonomial = arrayF;
  size_t RatValuesSize = sizeof(RatValues) / sizeof(RatValues[0]);
  computeFormula(&pstMonomial, &Rat5);    // Compute F.
  pstMonomial = arrayB;
  for (int index = 0; index < (int)RatValuesSize; index++)
  {
    computeFormula(&pstMonomial, pRatValues);
    intToBigInteger(&Rat3.numerator, 1);
    intToBigInteger(&Rat3.denominator, 1);
    for (int ctr = 1; ctr <= 5; ctr++)
    {
      BigRationalMultiply(&Rat3, &RatRoot, &Rat3);
      computeFormula(&pstMonomial, &Rat4);
      BigRationalMultiply(&Rat3, &Rat4, &Rat4);
      BigRationalAdd(pRatValues, &Rat4, pRatValues);
    }
    BigRationalDivide(pRatValues, &Rat5, pRatValues);           // Divide by F.
    ForceDenominatorPositive(pRatValues);
    pRatValues++;
  }
  BigRationalDivide(&RatValues[index_T2], &RatDiscr, &RatValues[index_T2]);   // Divide by discriminant.
  BigRationalDivide(&RatValues[index_T4], &RatDiscr, &RatValues[index_T4]);   // Divide by discriminant.
  BigRationalDivide(&RatValues[index_v], &RatDiscr, &RatValues[index_v]);     // Divide by discriminant.
  BigRationalDivide(&RatValues[index_O], &RatDiscr, &RatValues[index_O]);     // Divide by discriminant.
  BigRationalDivideByInt(&RatValues[index_T1], -2, &RatValues[index_T1]);     // Divide by 2.
  ForceDenominatorPositive(&RatValues[index_T1]);
  BigRationalDivideByInt(&RatValues[index_T2], 2, &RatValues[index_T2]);      // Divide by 2.
  BigRationalDivideByInt(&RatValues[index_T3], 2, &RatValues[index_T3]);      // Divide by 2.
  BigRationalDivideByInt(&RatValues[index_T4], 2, &RatValues[index_T4]);      // Divide by 2.
  BigRationalDivideByInt(&RatValues[index_v], 2, &RatValues[index_v]);        // Divide by 2.
}

void start5thRoot(void)
{
  if (pretty == PRETTY_PRINT)
  {
    showText("<span class=\"root\"><span class=\"befrad\" aria-hidden=\"true\">5</span><span class=\"radicand5\">");
  }
  else if (pretty == TEX)
  {
    showText("\\sqrt[5]{");
  }
  else
  {
    showText("(");
  }
}

void end5thRoot(void)
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
    showText(")^(1/5)");
  }
}

static void showSqrt5(void)
{
  if (pretty == PRETTY_PRINT)
  {
    startSqrt();
    *ptrOutput = '5';
    ptrOutput++;
    endSqrt();
  }
  else if (pretty == TEX)
  {
    showText("\\sqrt{5}");
  }
  else
  {
    showText("5^(1/2)");
  }
}

static void showMinusOnePlusMinusSqrt5(const char *sign)
{
  *ptrOutput = '(';
  ptrOutput++;
  showText(ptrMinus);
  *ptrOutput = '1';
  ptrOutput++;
  showText(sign);
  showSqrt5();
  *ptrOutput = ')';
  ptrOutput++;
}

static void showSqrtTenPlusMinusTwoTimesSqrt5(const char *sign)
{
  startSqrt();
  showText("10 ");
  showText(sign);
  *ptrOutput = ' ';
  ptrOutput++;
  showText(" 2");
  if (pretty == PARI_GP)
  {
    *ptrOutput = ' ';
    ptrOutput++;
  }
  showText(ptrTimes);
  if (pretty == PARI_GP)
  {
    *ptrOutput = ' ';
    ptrOutput++;
  }
  showSqrt5();
  endSqrt();
}

// Show (M +/- N*d)^(1/2)
static bool showSqRoot2(enum eSign signSqrt5, enum eSign signMN)
{
  enum eSign signBak = RatN.numerator.sign;
  RatN.numerator.sign = SIGN_POSITIVE;
  CopyBigInt(&Rat3.numerator, &RatN.numerator);
  CopyBigInt(&Rat3.denominator, &RatN.denominator);
  CopyBigInt(&Rat4.numerator, &RatDiscr.numerator);
  CopyBigInt(&Rat4.denominator, &RatDiscr.denominator);
  MultiplyRationalBySqrtRational(&Rat3, &Rat4);
  if (BigIntIsOne(&Rat4.numerator) && BigIntIsOne(&Rat4.denominator))
  {   // No square root inside square root.
      // Compute rational argument of square root.
    if (signMN == signBak)
    {
      BigRationalAdd(&RatM, &Rat3, &Rat3);
    }
    else
    {
      BigRationalSubt(&RatM, &Rat3, &Rat3);
    }
    if (BigRationalSquareRoot(&Rat3, &Rat4))
    {
      RatN.numerator.sign = signBak;
      return true;
    }
    MultiplyRationalBySqrtRational(&Rat4, &Rat3);
    ShowRationalAndSqrParts(&Rat4, &Rat3, 2, ptrTimes);
    RatN.numerator.sign = signBak;
    return false;
  }
  startSqrt();
  if (!BigIntIsZero(&RatM.numerator))
  {
    showRationalNoParen(&RatM);
    showPlusSignOn(signMN == signBak, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
  }
  ShowRationalAndSqrParts(&Rat3, &Rat4, 2, ptrTimes);
  endSqrt();
  RatN.numerator.sign = signBak;
  return false;
}

static void ShowQuinticsRootsRealR(int multiplicity)
{
  const char* ptrDenominator;
  ptrDenominator = "4";
  BigRationalDivideByInt(&RatQuartic, -5, &RatQuartic);
  ForceDenominatorPositive(&RatQuartic);
  startLine();
  if (pretty == PRETTY_PRINT)
  {
    showText("<var>S</var><sub>1</sub> = ");
  }
  else
  {
    showExpr("S_1 = ");
  }
  showMinusOnePlusMinusSqrt5("+");
  showText(ptrTimes);
  if (pretty == PRETTY_PRINT)
  {
    showText("(<var>R</var><sub>1</sub> + <var>R</var><sub>4</sub>) + ");
  }
  else
  {
    showExpr("(R_1 + R_4) + ");
  }
  showMinusOnePlusMinusSqrt5(ptrMinus);
  showText(ptrTimes);
  if (pretty == PRETTY_PRINT)
  {
    showText("(<var>R</var><sub>2</sub> + <var>R</var><sub>3</sub>)");
    endLine();
    startLine();
    showText("<var>S</var><sub>2</sub> = ");
  }
  else
  {
    showExpr("(R_2 + R_3)");
    endLine();
    startLine();
    showExpr("S_2 = ");
  }
  showMinusOnePlusMinusSqrt5("+");
  showText(ptrTimes);
  if (pretty == PRETTY_PRINT)
  {
    showText("(<var>R</var><sub>2</sub> + <var>R</var><sub>3</sub>) + ");
  }
  else
  {
    showExpr("(R_2 + R_3) + ");
  }
  showMinusOnePlusMinusSqrt5(ptrMinus);
  showText(ptrTimes);
  if (pretty == PRETTY_PRINT)
  {
    showText("(<var>R</var><sub>1</sub> + <var>R</var><sub>4</sub>)");
    endLine();
    startLine();
    showText("<var>T</var><sub>1</sub> = ");
  }
  else
  {
    showExpr("(R_1 + R_4)");
    endLine();
    startLine();
    showExpr("T_1 = ");
  }
  showSqrtTenPlusMinusTwoTimesSqrt5("+");
  showText(ptrTimes);
  if (pretty == PRETTY_PRINT)
  {
    showText("(<var>R</var><sub>4</sub> ");
    showText(ptrMinus);
    showText(" <var>R</var><sub>1</sub>) +");
  }
  else
  {
    showExpr("(R4 - R1) +");
  }
  if (pretty == PARI_GP)
  {
    *ptrOutput = ' ';
    ptrOutput++;
  }
  showSqrtTenPlusMinusTwoTimesSqrt5(ptrMinus);
  showText(ptrTimes);
  if (pretty == PRETTY_PRINT)
  {
    showText("(<var>R</var><sub>3</sub> ");
    showText(ptrMinus);
    showText(" <var>R</var><sub>2</sub>)");
    endLine();
    startLine();
    showText("<var>T</var><sub>2</sub> = ");
  }
  else
  {
    showExpr("(R_3 - R_2)");
    endLine();
    startLine();
    showExpr("T_2 = ");
  }
  showSqrtTenPlusMinusTwoTimesSqrt5("+");
  showText(ptrTimes);
  if (pretty == PRETTY_PRINT)
  {
    showText("(<var>R</var><sub>3</sub> ");
    showText(ptrMinus);
    showText(" <var>R</var><sub>2</sub>) +");
  }
  else
  {
    showExpr("(R_3 - R_2) +");
  }
  if (pretty == PARI_GP)
  {
    *ptrOutput = ' ';
    ptrOutput++;
  }
  showSqrtTenPlusMinusTwoTimesSqrt5(ptrMinus);
  showText(ptrTimes);
  if (pretty == PRETTY_PRINT)
  {
    showText("(<var>R</var><sub>1</sub> ");
    showText(ptrMinus);
    showText(" <var>R</var><sub>4</sub>)");
  }
  else
  {
    showExpr("(R_1 - R_4)");
  }
  endLine();
  for (int ctr = 1; ctr <= 5; ctr++)
  {
    showX(multiplicity);
    if (!BigIntIsZero(&RatQuartic.numerator))
    {
      showRationalNoParen(&RatQuartic);
      showText(" + ");
    }
    switch (ctr)
    {
    case 1:
      if (pretty == PRETTY_PRINT)
      {
        const char* ptrNumerator = "<var>R</var><sub>1</sub> + <var>R</var><sub>2</sub> + <var>R</var><sub>3</sub> + <var>R</var><sub>4</sub>";
        showText(ptrNumerator);
      }
      else if (pretty == TEX)
      {
        const char* ptrNumerator = "R_1 + R_2 + R_3 + R_4";
        showText(ptrNumerator);
      }
      else
      {
        const char* ptrNumerator = "R1 + R2 + R3 + R4";
        showText(ptrNumerator);
      }
      break;
    case 2:
      if (pretty == PRETTY_PRINT)
      {
        showRatConstants("<var>S</var><sub>1</sub>", ptrDenominator);
        showText("+ i ");
        showText(ptrTimes);
        showRatConstants("<var>T</var><sub>1</sub>", ptrDenominator);
      }
      else if (pretty == TEX)
      {
        showRatConstants("S_1", ptrDenominator);
        showText("+ i ");
        showRatConstants("T_1", ptrDenominator);
      }
      else
      {
        showText("(S1 + I * T1) / ");
        showText(ptrDenominator);
      }
      break;
    case 3:
      if (pretty == PRETTY_PRINT)
      {
        showRatConstants("<var>S</var><sub>1</sub>", ptrDenominator);
        showText(ptrMinus);
        showText(" i ");
        showText(ptrTimes);
        showRatConstants("<var>T</var><sub>1</sub>", ptrDenominator);
      }
      else if (pretty == TEX)
      {
        showRatConstants("S_1", ptrDenominator);
        showText("- i ");
        showRatConstants("T_1", ptrDenominator);
      }
      else
      {
        showText("(S1 - I * T1) / ");
        showText(ptrDenominator);
      }
      break;
    case 4:
      if (pretty == PRETTY_PRINT)
      {
        showRatConstants("<var>S</var><sub>2</sub>", ptrDenominator);
        showText("+ i ");
        showText(ptrTimes);
        showRatConstants("<var>T</var><sub>2</sub>", ptrDenominator);
      }
      else if (pretty == TEX)
      {
        showRatConstants("S_2", ptrDenominator);
        showText("+ i ");
        showRatConstants("T_2", ptrDenominator);
      }
      else
      {
        showText("(S2 + I * T2) / ");
        showText(ptrDenominator);
      }
      break;
    case 5:
      if (pretty == PRETTY_PRINT)
      {
        showRatConstants("<var>S</var><sub>2</sub>", ptrDenominator);
        showText(ptrMinus);
        showText(" i ");
        showText(ptrTimes);
        showRatConstants("<var>T</var><sub>2</sub>", ptrDenominator);
      }
      else if (pretty == TEX)
      {
        showRatConstants("S_2", ptrDenominator);
        showText("- i ");
        showText(ptrTimes);
        showRatConstants("T_2", ptrDenominator);
      }
      else
      {
        showText("(S2 - I * T2) / ");
        showText(ptrDenominator);
      }
      break;
    default:
      break;
    }
    endLine();
  }
}

static void getSignsSqrt(int groupOrder, int ctr,
  enum eSign *pFirstSign, enum eSign *pSecondSign)
{
  int ctrOrder;
  if (RatValues[index_O].numerator.sign == SIGN_POSITIVE)
  {   // O > 0.
    ctrOrder = 3;
  }
  else
  {   // O < 0.
    ctrOrder = 2;
  }
  if (groupOrder == 10)
  {
    *pFirstSign = (((ctr == ctrOrder) || (ctr == 4)) ? SIGN_NEGATIVE : SIGN_POSITIVE);
    *pSecondSign = (((ctr == ctrOrder) || (ctr == 1)) ? SIGN_POSITIVE : SIGN_NEGATIVE);
  }
  else
  {
    *pFirstSign = (((ctr == ctrOrder) || (ctr == 4)) ? SIGN_POSITIVE : SIGN_NEGATIVE);
    *pSecondSign = (((ctr == ctrOrder) || (ctr == 1)) ? SIGN_POSITIVE : SIGN_NEGATIVE);
  }
}

// Compute Rn so the argument of the fifth root shown is always positive.
  // R_i = l0 + T1 +/- T2*sqrt(5*d) +/- (1/4)*(sqrt(5)+1)*sqrt(M+N*d)
  //                                +/- (1/4)*(sqrt(5)-1)*sqrt(M-N*d)
  // First  +/-: plus for R_2 and R_3, minus for R_1 and R_4.
  // Second +/-: if O > 0: plus for R_3 and R_4, minus for R_1 and R_2.
  //             if O < 0: plus for R_2 and R_4, minus for R_1 and R_3.
  // Third  +/-: if O > 0: plus for R_1 and R_3, minus for R_2 and R_4.
  //             if O < 0: plus for R_1 and R_2, minus for R_3 and R_4.
static double computeRnGroup10Or20(int groupOrder, int ctr,
  bool changeSignFirstSqrt)
{
  double dTemp;
  double dTemp2;
  double dInsideSqrt;
  enum eSign firstSign;
  enum eSign secondSign;
  double dFifthRootArgument = BigRational2double(&RatValues[index_l0]) +
    BigRational2double(&RatValues[index_T1]) * 0.5;
  dTemp = BigRational2double(&RatValues[index_T2]) *
    sqrt(BigRational2double(&RatDiscr) * 0.25);
  if ((ctr == 2) || (ctr == 3))
  {
    dFifthRootArgument += dTemp;
  }
  else
  {
    dFifthRootArgument -= dTemp;
  }
  getSignsSqrt(groupOrder, ctr, &firstSign, &secondSign);
  dTemp2 = BigRational2double(&RatM);
  dTemp = BigRational2double(&RatN) * sqrt(BigRational2double(&RatDiscr));
  if (((groupOrder == 20) && ((ctr == 2) || (ctr == 3))) ||
    ((groupOrder == 10) && ((ctr == 1) || (ctr == 4))))
  {
    dTemp = -dTemp;
  }
  dInsideSqrt = sqrt(dTemp2 + dTemp);
  if (changeSignFirstSqrt && (firstSign != secondSign) && (dTemp < 0.0))
  {
    if (firstSign == SIGN_NEGATIVE)
    {
      dFifthRootArgument += dInsideSqrt;
    }
    else
    {
      dFifthRootArgument -= dInsideSqrt;
    }
  }
  else
  {
    if (firstSign == SIGN_POSITIVE)
    {
      dFifthRootArgument += dInsideSqrt;
    }
    else
    {
      dFifthRootArgument -= dInsideSqrt;
    }
  }
  return dFifthRootArgument;
}

// Change signs if the value of R_n computed is negative, so 
// the argument of the fifth root is always positive.
static void showRnGroup10Or20(int groupOrder, int ctr, double dFifthRootArgument,
  bool changeSignFirstSqrt)
{
  bool signPlus;
  enum eSign firstSign;
  enum eSign secondSign;
  enum eSign currentSign;
  enum eSign signInsideFirstSquareRoot;
  startLine();
  if (pretty == PRETTY_PRINT)
  {
    showText("<var>R</var><sub>");
    *ptrOutput = (char)(ctr + '0');
    ptrOutput++;
    showText("</sub></var>");
  }
  else
  {
    showExpr("R_");
    *ptrOutput = (char)(ctr + '0');
    ptrOutput++;
  }
  showText(" = ");
  BigRationalDivideByInt(&RatValues[index_T2], 2, &Rat3);
  ForceDenominatorPositive(&Rat3);
  signPlus = ((ctr == 2) || (ctr == 3));
  BigRationalMultiplyByInt(&RatValues[index_l0], 1, &Rat1);
  BigRationalDivideByInt(&RatValues[index_T1], 2, &Rat2);
  BigRationalAdd(&Rat1, &Rat2, &Rat1);   // l0 + t1/2.
  if (dFifthRootArgument < 0.0)
  {                      // Fifth root is negative.
    showText(ptrMinus);  // Prepend minus sign to fifth root.
    BigIntChSign(&Rat1.numerator);
    BigIntChSign(&Rat3.numerator);
  }
  if (Rat3.numerator.sign == SIGN_NEGATIVE)
  {
    signPlus = !signPlus;
    Rat3.numerator.sign = SIGN_POSITIVE;
  }
  start5thRoot();
  CopyBigInt(&Rat4.numerator, &RatDiscr.numerator);
  CopyBigInt(&Rat4.denominator, &RatDiscr.denominator);
  MultiplyRationalBySqrtRational(&Rat3, &Rat4);
  if (BigIntIsOne(&Rat4.numerator) && BigIntIsOne(&Rat4.denominator))
  {
    if (signPlus)
    {
      BigRationalAdd(&Rat1, &Rat3, &Rat1);
    }
    else
    {
      BigRationalSubt(&Rat1, &Rat3, &Rat1);
    }
    if (!BigIntIsZero(&Rat1.numerator))
    {
      showRationalNoParen(&Rat1);        // Show l0 + t1/2 +/- t2/2.
    }
  }
  else
  {
    if (!BigIntIsZero(&Rat1.numerator))
    {
      showRationalNoParen(&Rat1);        // Show l0 + t1/2.
    }
    // The rational part of this term can be added to or subtracted from l0.
    showPlusSignOn(signPlus, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
    ShowRationalAndSqrParts(&Rat3, &Rat4, 2, ptrTimes);  // Show +/- T2*sqrt(discr)
  }
  getSignsSqrt(groupOrder, ctr, &firstSign, &secondSign);
  if (dFifthRootArgument > 0.0)
  {
    currentSign = SIGN_POSITIVE;
  }
  else
  {
    currentSign = SIGN_NEGATIVE;
  }
  if (((groupOrder == 20) && ((ctr == 2) || (ctr == 3))) ||
    ((groupOrder == 10) && ((ctr == 1) || (ctr == 4))))
  {
    signInsideFirstSquareRoot = SIGN_NEGATIVE;
  }
  else
  {
    signInsideFirstSquareRoot = SIGN_POSITIVE;
  }
  if (changeSignFirstSqrt && (firstSign != secondSign) &&
    (signInsideFirstSquareRoot != RatN.numerator.sign))
  {
    showPlusSignOn(firstSign != currentSign,
      TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
  }
  else
  {
    showPlusSignOn(firstSign == currentSign,
      TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
  }
  if (showSqRoot2(SIGN_POSITIVE, signInsideFirstSquareRoot))
  {
    showRationalNoParen(&Rat4);
  }
  end5thRoot();
}

static void GaloisGroupHasOrder10Or20(int groupOrder)
{
  bool changeSignFirstSqrt;
  // Change signs of M and N, so M + Nd and M - Nd are positive.
  BigIntChSign(&RatM.numerator);
  BigIntChSign(&RatN.numerator);
  CopyBigInt(&Rat1.numerator, &RatDiscr.numerator);
  CopyBigInt(&Rat1.denominator, &RatDiscr.denominator);
  MultiplyRationalBySqrtRational(&RatN, &Rat1);
  MultiplyRationalBySqrtRational(&RatValues[index_T2], &RatDiscr);
  BigRationalMultiply(&RatM, &RatM, &RatR);   // RatR <- M^2
  BigRationalDivide(&RatR, &RatDiscr, &RatR); // RatR <- M^2/d^2
  BigRationalMultiply(&RatN, &RatN, &RatS);   // RatS <- N^2
  // Check whether the first sign of square root has to be changed.
  // Compare magnitude of A = (5+sqrt(5))*(M-abs(N)*d) versus
  // B = (5-sqrt(5))*(M+abs(N)*d)
  // A - B = 2*M*sqrt(5) - 2*5*abs(N)*d
  // Compare sqrt(5)*M versus 5*N*d
  // Compare 5*M^2 versus 25*N^2*d^2
  // Compare M^2/(5*N^2*d^2) versus 1 (numerator vs. denominator).
  BigRationalDivide(&RatR, &RatS, &Rat2);    // Rat2 <- M^2/(d^2*N^2)
  BigRationalDivideByInt(&Rat2, 5, &Rat2);    // Rat2 <- M^2/(5*d^2*N^2)
  // Compare numerator vs. denominator.
  BigIntSubt(&Rat2.numerator, &Rat2.denominator, &tmp1);
  changeSignFirstSqrt = (tmp1.sign == SIGN_NEGATIVE);
  BigRationalSubt(&RatR, &RatS, &RatR);
  BigRationalDivideByInt(&RatValues[index_l0], 3125, &RatValues[index_l0]);
  BigRationalDivideByInt(&RatValues[index_T1], 3125, &RatValues[index_T1]);
  BigRationalDivideByInt(&RatValues[index_T2], 3125, &RatValues[index_T2]);
  (void)BigRationalSquareRoot(&RatR, &Rat2);   // Rat2 <- sqrt(M^2/d^2 - N^2).
  BigRationalMultiplyByInt(&RatM, 10, &RatM);
  BigRationalMultiplyByInt(&Rat2, 2, &Rat2);
  if (groupOrder == 10)
  {
    RatN.numerator.sign = SIGN_NEGATIVE;
  }
  BigRationalSubt(&RatN, &Rat2, &RatN);
  BigRationalMultiplyByInt(&RatN, 2, &RatN);
  BigRationalMultiplyByInt(&RatDiscr, 5, &RatDiscr);
  BigRationalDivideByInt(&RatM, 3125 * 3125 * 8, &RatM);
  BigRationalDivideByInt(&RatN, 3125 * 3125 * 8, &RatN);
  for (int ctr = 1; ctr <= 4; ctr++)
  {
    double dRn = computeRnGroup10Or20(groupOrder, ctr, changeSignFirstSqrt);
    showRnGroup10Or20(groupOrder, ctr, dRn, changeSignFirstSqrt);
  }
}

static int getCircleNbr(int circleNbrLeft)
{
  int recordCircleNbrRight = 0;
  double rootRecord = 0.0;
  double firstProd = BigRational2double(&Rat1) * sqrt(BigRational2double(&Rat2));
  double firstTerm = BigRational2double(&Rat3) * sqrt(BigRational2double(&Rat4));
  double secondTerm = BigRational2double(&RatM) * sqrt(BigRational2double(&RatN));
  double firstAcos = acos(firstTerm + secondTerm);
  double secondAcos = acos(firstTerm - secondTerm);
  for (int circleNbrRight = 0; circleNbrRight < 5; circleNbrRight++)
  {
    double RHS;
    // Compute root = RatQuartic + Rat1 * sqrt(Rat2) *
    //     (cos((1/5)(2*circleNbrLeft*Pi + acos(Rat3 * sqrt(Rat4) + RatM * sqrt(RatN)))) +
    //      cos((1/5)(2*circleNbrRight*Pi + acos(Rat3 * sqrt(Rat4) - RatM * sqrt(RatN)))))
    double argFirstCos = (6.2831853071795864769252 * (double)circleNbrLeft + firstAcos)/5.0;
    double argSecondCos = (6.2831853071795864769252 * (double)circleNbrRight + secondAcos)/5.0;
    double root = BigRational2double(&RatQuartic) + (firstProd *
      (cos(argFirstCos) + cos(argSecondCos)));
    // Compute A*root^5 + B*root^4 + C*root^3 + D*root^2 + E*root + F.
    RHS = ((((((((((BigInt2double(&Quintic) * root) + 
      BigInt2double(&Quartic)) * root) +
      BigInt2double(&Cubic)) * root) +
      BigInt2double(&Quadratic))* root) +
      BigInt2double(&Linear)) * root) +
      BigInt2double(&Independent));
    if (RHS < 0.0)
    {   // Compute absolute value of RHS.
      RHS = -RHS;
    }
    if (circleNbrRight == 0)
    {
      rootRecord = RHS;
    }
    else
    {
      if (RHS < rootRecord)
      {
        rootRecord = RHS;
        recordCircleNbrRight = circleNbrRight;
      }
    }
  }
  return recordCircleNbrRight;
}

static void GaloisGroupHasOrder5(int multiplicity)
{
  (void)multiplicity;
  int circleNbr[2] = { 0, 0 };
  BigRational* ptrL;
  BigRational* ptrM;
  // Compute the roots l1, l2, l3 and l4 of the pair of quadratic equations.
  // l1, l4 = (-(T1 + T2*d) +/- sqrt(M + N*d)) / 2
  // l2, l3 = (-(T1 - T2*d) +/- sqrt(M - N*d)) / 2
  // The order of the roots must comply with: (l1 − l4)(l2 − l3) = O*d
  // If O > 0 -> l1 and l2: plus sign, l3 and l4: minus sign
  // If O < 0 -> l1 and l3: plus sign, l2 and l4: minus sign
  // sqrt(M + N*d) is already computed as Rat3 and d is Rat4.
  // Compute sqrt(M - N*d) as Rat5.
  BigRationalMultiply(&Rat4, &RatN, &Rat1);
  BigRationalSubt(&RatM, &Rat1, &Rat1);
  ForceDenominatorPositive(&Rat1);         // M - N*d
  squareRoot(Rat1.numerator.limbs, tmp4.limbs, Rat1.numerator.nbrLimbs, &tmp4.nbrLimbs);
  squareRoot(Rat1.denominator.limbs, tmp1.limbs, Rat1.denominator.nbrLimbs, &tmp1.nbrLimbs);
  CopyBigInt(&Rat5.numerator, &tmp4);      // Rat5 = square root of discriminant
  CopyBigInt(&Rat5.denominator, &tmp1);    // of quadratic equation.
  BigRationalDivideByInt(&Rat5, 2, &Rat5);

  BigRationalDivideByInt(&Rat3, 2, &Rat3);
  BigRationalMultiply(&RatValues[index_T2], &Rat4, &Rat2);    // T2*d
  BigRationalAdd(&RatValues[index_T1], &Rat2, &Rat1);         // T1 + T2*d
  BigRationalDivideByInt(&Rat1, -2, &Rat1);
  ForceDenominatorPositive(&Rat1);                            // -(T1 + T2*d)/2
  BigRationalSubt(&RatValues[index_T1], &Rat2, &Rat2);        // T1 - T2*d
  BigRationalDivideByInt(&Rat2, -2, &Rat2);
  ForceDenominatorPositive(&Rat2);                            // -(T1 - T2*d)/2
  BigRationalAdd(&Rat1, &Rat3, &RatValues[index_T1]);         // l1
  BigRationalSubt(&Rat1, &Rat3, &RatValues[index_T4]);        // l4
  if (RatValues[index_O].numerator.sign == SIGN_POSITIVE)
  {   // O > 0.
    BigRationalAdd(&Rat2, &Rat5, &RatValues[index_T2]);       // l2
    BigRationalSubt(&Rat2, &Rat5, &RatValues[index_T3]);      // l3
  }
  else
  {   // O < 0.
    BigRationalAdd(&Rat2, &Rat5, &RatValues[index_T3]);       // l3
    BigRationalSubt(&Rat2, &Rat5, &RatValues[index_T2]);      // l2
  }
  // l = l0 - (l1 + l2 + l3 + l4)/4
  // m = (l1 - l2 - l3 + l4) / 4
  // U = sqrt(10+2*sqrt(5)) / 4
  // V = sqrt(10-2*sqrt(5)) / 4
  // 
  // Roots of resolvent:
  // R1 = l + m * sqrt(5) + i((l1 - l4)*U + (l2 - l3)*V)
  // R2 = l - m * sqrt(5) + i((l3 - l2)*U + (l1 - l4)*V)
  // R3 = l - m * sqrt(5) + i((l2 - l3)*U + (l4 - l1)*V)
  // R4 = l + m * sqrt(5) + i((l4 - l1)*U + (l3 - l2)*V)
  // Do not compute these roots.
  //
  // Re(R1) = sqrt(r) * ((l + m) * sqrt(5))/K
  // r = l4*(2*l4-l3-l2-l1-l0)/2 + l3*(2*l3-l2-l1-l0)/2 +
  //     l2*(2*l2-l1-l0)/2 + l1*(2*l1-l0)/2 + l0^2
  // r is a fifth power. Let s = r^(1/5)
  // If Q is the quartic coefficient of the original polynomial:
  // x_n = (1/5)(-Q + 2*sqrt(s)*(cos((1/5)*(2*n*Pi+acos((l+m*sqrt(5)/sqrt(r)) +
  //                             cos((1/5)*(2*k_n*Pi+acos((l-m*sqrt(5)/sqrt(r))
  // Compute k_n by inspection minimizing the RHS for the computed value of x.
  // Compute gcd of l, m and sqrt(x) to minimize the numerator and denominator
  // of the argument of acos.
  //
  // Compute ptrL <- l and ptrM <- m.
  ptrL = &RatValues[index_v];
  ptrM = &RatValues[index_O];
  BigRationalAdd(&RatValues[index_T1], &RatValues[index_T4], &Rat2);
  BigRationalAdd(&RatValues[index_T2], &RatValues[index_T3], &Rat3);
  BigRationalAdd(&Rat2, &Rat3, &Rat1);     // l1 + l2 + l3 + l4
  BigRationalSubt(&Rat2, &Rat3, &Rat2);    // l1 - l2 - l3 + l4
  BigRationalDivideByInt(&Rat1, 4, &Rat1); // (l1 + l2 + l3 + l4)/4
  BigRationalDivideByInt(&Rat2, 4, ptrM);  // (l1 - l2 - l3 + l4)/4
  BigRationalSubt(&RatValues[index_l0], &Rat1, ptrL);  // l
  ptrM->numerator.sign = SIGN_POSITIVE;

  // Compute RatR <- r.
  BigRationalAdd(&RatValues[index_T4], &RatValues[index_T4], &Rat2);
  BigRationalSubt(&Rat2, &RatValues[index_T3], &Rat2);
  BigRationalSubt(&Rat2, &RatValues[index_T2], &Rat2);
  BigRationalSubt(&Rat2, &RatValues[index_T1], &Rat2);
  BigRationalSubt(&Rat2, &RatValues[index_l0], &Rat2);
  BigRationalMultiply(&Rat2, &RatValues[index_T4], &RatR);

  BigRationalAdd(&RatValues[index_T3], &RatValues[index_T3], &Rat2);
  BigRationalSubt(&Rat2, &RatValues[index_T2], &Rat2);
  BigRationalSubt(&Rat2, &RatValues[index_T1], &Rat2);
  BigRationalSubt(&Rat2, &RatValues[index_l0], &Rat2);
  BigRationalMultiply(&Rat2, &RatValues[index_T3], &Rat2);
  BigRationalAdd(&RatR, &Rat2, &RatR);

  BigRationalAdd(&RatValues[index_T2], &RatValues[index_T2], &Rat2);
  BigRationalSubt(&Rat2, &RatValues[index_T1], &Rat2);
  BigRationalSubt(&Rat2, &RatValues[index_l0], &Rat2);
  BigRationalMultiply(&Rat2, &RatValues[index_T2], &Rat2);
  BigRationalAdd(&RatR, &Rat2, &RatR);

  BigRationalAdd(&RatValues[index_T1], &RatValues[index_T1], &Rat2);
  BigRationalSubt(&Rat2, &RatValues[index_l0], &Rat2);
  BigRationalMultiply(&Rat2, &RatValues[index_T1], &Rat2);
  BigRationalAdd(&RatR, &Rat2, &RatR);
  BigRationalDivideByInt(&RatR, 2, &RatR);

  BigRationalMultiply(&RatValues[index_l0], &RatValues[index_l0], &Rat2);
  BigRationalAdd(&RatR, &Rat2, &RatR);

  // Compute s <- (2/5)r^(1/5).
  computeRoot(&RatR.numerator, &Rat2.numerator, 5);
  computeRoot(&RatR.denominator, &Rat2.denominator, 5);
  intToBigInteger(&Rat1.numerator, 2);
  intToBigInteger(&Rat1.denominator, 5);
  MultiplyRationalBySqrtRational(&Rat1, &Rat2);

  // Compute L/sqrt(R): result in Rat3, Rat4.
  CopyBigInt(&Rat3.numerator, &ptrL->numerator);
  (void)BigIntMultiply(&ptrL->denominator, &RatR.numerator,
    &Rat3.denominator);
  (void)BigIntMultiply(&RatR.numerator, &RatR.denominator,
    &Rat4.numerator);
  intToBigInteger(&Rat4.denominator, 1);
  MultiplyRationalBySqrtRational(&Rat3, &Rat4);

  // Compute M*sqrt(5)/sqrt(R): result in RatM, RatN.
  CopyBigInt(&RatM.numerator, &ptrM->numerator);
  (void)BigIntMultiply(&ptrM->denominator, &RatR.numerator,
    &RatM.denominator);
  (void)BigIntMultiply(&RatR.numerator, &RatR.denominator,
    &RatN.numerator);
  multint(&RatN.numerator, &RatN.numerator, 5);
  intToBigInteger(&RatN.denominator, 1);
  MultiplyRationalBySqrtRational(&RatM, &RatN);
  BigIntChSign(&RatQuartic.numerator);
  BigRationalDivideByInt(&RatQuartic, 5, &RatQuartic);
  for (circleNbr[0] = 0; circleNbr[0] < 5; circleNbr[0]++)
  {
    showX(multiplicity);
    if (!BigIntIsZero(&RatQuartic.numerator))
    {
      showRationalNoParen(&RatQuartic);
      showText(" + ");
    }
    ShowRationalAndSqrParts(&Rat1, &Rat2, 2, ptrTimes);
    showText(ptrTimes);
    startParen();
    for (int cosNbr = 1; cosNbr <= 2; cosNbr++)
    {
      int circleNumber = circleNbr[cosNbr - 1];
      showText(ptrCos);
      startParen();
      showRatConstants("1", "5");
      showText(ptrTimes);
      if (circleNumber > 0)
      {
        startParen();
        int2dec(&ptrOutput, 2 * circleNumber);
        showText(ptrTimes);
        showText(ptrPi);
        showText(" + ");
      }
      showText(ptrACos);
      startParen();
      ShowRationalAndSqrParts(&Rat3, &Rat4, 2, ptrTimes);
      if (cosNbr == 1)
      {
        showText(" + ");
      }
      else
      {
        showText(" ");
        showText(ptrMinus);
        showText(" ");
      }
      ShowRationalAndSqrParts(&RatM, &RatN, 2, ptrTimes);
      endParen();   // End paren for arc cos.
      if (circleNumber > 0)
      {
        endParen();
      }
      endParen();   // End paren for cos.
      if (pretty == TEX)
      {
        showText("}");
      }
      if (cosNbr == 1)
      {
        showText(" + ");
        circleNbr[1] = getCircleNbr(circleNbr[0]);
      }
    }
    endParen();
    endLine();
  }
}

void QuinticEquation(const int* ptrPoly, int multiplicity)
{
  struct monomial* pstMonomial;
  const int* ptr;
  const int* ptrPolynomial = ptrPoly;
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
  ptrPolynomial += numLimbs(ptrPolynomial);
  ptrPolynomial++;
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
  // s = (3125f - 625be + 125b^2*d - 25b^3*c + 4b^5)/3125
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
  BigRationalMultiply(&Rat1, &RatQuadratic, &Rat2);
  BigRationalMultiplyByInt(&Rat2, 125, &Rat2);                     // 125b^2*d
  BigRationalAdd(&RatDeprIndependent, &Rat2, &RatDeprIndependent);
  BigRationalMultiply(&Rat1, &RatQuartic, &Rat1);                  // b^3
  BigRationalMultiplyByInt(&Rat1, 4, &Rat2);                       // 4b^3
  BigRationalAdd(&RatDeprQuadratic, &Rat2, &RatDeprQuadratic);
  BigRationalMultiply(&Rat1, &RatCubic, &Rat2);
  BigRationalMultiplyByInt(&Rat2, -25, &Rat2);                     // -25b^3*c
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

  FactorPolynomialF20();
  // If there is a linear factor, that means that the quintic can be expressed with radicands.
  if (factorInfoInteger[0].degree != 1)
  {
    showText(lang ? "<p>La ecuación quíntica no se puede expresar mediante radicandos.</p>" :
      "The quintic equation cannot be expressed with radicands.");
    return;
  }
  ptr = factorInfoInteger[0].ptrPolyLifted;
  UncompressBigIntegerB(ptr, &RatRoot.numerator);          // Independent term.
  ptr += numLimbs(ptr);
  ptr++;
  UncompressBigIntegerB(ptr, &RatRoot.denominator);        // Linear coefficient.
  BigIntChSign(&RatRoot.numerator);                        // root = -independent / linear.
  // Compute discriminant
  pstMonomial = arrayD;
  computeFormula(&pstMonomial, &RatDiscr);
  computeArrayValues();
  // Let M = T1^2 + T2^2*D - 4*T3, N = 2*T1*T2 - 4*T4
  // and d = sqrt(D), then:
  // l1, l4 = (-(T1 + T2*d) +/- sqrt(M + N*d)) / 2
  // l2, l3 = (-(T1 - T2*d) +/- sqrt(M - N*d)) / 2
  // The order of the roots must comply with: (l1 − l4)(l2 − l3) = O*d
  // l1 - l4 = +/- sqrt(M + N*d)
  // l2 - l3 = +/- sqrt(M - N*d)
  // Compute M.
  BigRationalMultiply(&RatValues[index_T1], &RatValues[index_T1], &RatM);
  BigRationalMultiply(&RatValues[index_T2], &RatValues[index_T2], &Rat1);
  BigRationalMultiply(&Rat1, &RatDiscr, &Rat1);
  BigRationalAdd(&RatM, &Rat1, &RatM);
  BigRationalMultiplyByInt(&RatValues[index_T3], 4, &Rat1);
  BigRationalSubt(&RatM, &Rat1, &RatM);
  // Compute N.
  BigRationalMultiply(&RatValues[index_T1], &RatValues[index_T2], &RatN);
  BigRationalMultiplyByInt(&RatN, 2, &RatN);
  BigRationalMultiplyByInt(&RatValues[index_T4], 4, &Rat1);
  BigRationalSubt(&RatN, &Rat1, &RatN);
  // The real part of the primitive fifth root of 1 is Re5 = (-1+sqrt(5))/4
  // The real part of its square is: Re5_2 = (-1-sqrt(5))/4
  if (!BigRationalSquareRoot(&RatDiscr, &Rat4))
  {   // Discriminant is not a perfect square.
      // Then the Galois group of f(x) is the Frobenius group of order 20.
    GaloisGroupHasOrder10Or20(20);
    ShowQuinticsRootsRealR(multiplicity);
    return;
  }
  // Test whether M+N*d is a perfect square. In this case the quadratic is
  // reducible. d = Rat4 at this moment.
  BigRationalMultiply(&Rat4, &RatN, &Rat1);
  BigRationalAdd(&Rat1, &RatM, &Rat1);
  ForceDenominatorPositive(&Rat1);         // M + N*d
  if (BigRationalSquareRoot(&Rat1, &Rat3))
  {  // M + N*d is perfect square. its square root is stored in Rat3.
    GaloisGroupHasOrder5(multiplicity);
    return;
  }
  GaloisGroupHasOrder10Or20(10);
  ShowQuinticsRootsRealR(multiplicity);
  return;
}

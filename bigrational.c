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
#include "bignbr.h"
#include "output.h"
#include "polynomial.h"

static BigInteger tmp1;
static BigInteger tmp2;
static BigInteger tmp3;
static BigRational Rat1;
static BigRational Rat2;

// Let A = a_n/a_d, B = b_n/b_d, C = c_n/d_n
// A = B+C means a_n = b_n*c_d + c_n*b_d, a_d = b_d * c_d
// Then divide A_n and A_d by their GCD.
enum eExprErr BigRationalAdd(BigRational* pAddend1, BigRational* pAddend2, BigRational* pSum)
{
  enum eExprErr rc;
  rc = BigIntMultiply(&pAddend1->numerator, &pAddend2->denominator, &tmp1);
  if (rc != EXPR_OK)
  {
    return rc;
  }
  rc = BigIntMultiply(&pAddend2->numerator, &pAddend1->denominator, &tmp2);
  if (rc != EXPR_OK)
  {
    return rc;
  }
  BigIntAdd(&tmp1, &tmp2, &tmp1);
  rc = BigIntMultiply(&pAddend1->denominator, &pAddend2->denominator, &tmp2);
  if (rc != EXPR_OK)
  {
    return rc;
  }
  // At this moment tmp1 is the numerator and tmp2 is the denominator.
  BigIntGcd(&tmp1, &tmp2, &tmp3);
  rc = BigIntDivide(&tmp1, &tmp3, &pSum->numerator);
  if (rc != EXPR_OK)
  {
    return rc;
  }
  return BigIntDivide(&tmp2, &tmp3, &pSum->denominator);
}

// Let A = a_n/a_d, B = b_n/b_d, C = c_n/d_n
// A = B-C means a_n = b_n*c_d - c_n*b_d, a_d = b_d * c_d
// Then divide A_n and A_d by their GCD.
enum eExprErr BigRationalSubt(BigRational* pAddend1, BigRational* pAddend2, BigRational* pSum)
{
  enum eExprErr rc;
  rc = BigIntMultiply(&pAddend1->numerator, &pAddend2->denominator, &tmp1);
  if (rc != EXPR_OK)
  {
    return rc;
  }
  rc = BigIntMultiply(&pAddend2->numerator, &pAddend1->denominator, &tmp2);
  if (rc != EXPR_OK)
  {
    return rc;
  }
  BigIntSubt(&tmp1, &tmp2, &tmp1);
  rc = BigIntMultiply(&pAddend1->denominator, &pAddend2->denominator, &tmp2);
  if (rc != EXPR_OK)
  {
    return rc;
  }
  // At this moment tmp1 is the numerator and tmp2 is the denominator.
  BigIntGcd(&tmp1, &tmp2, &tmp3);
  rc = BigIntDivide(&tmp1, &tmp3, &pSum->numerator);
  if (rc != EXPR_OK)
  {
    return rc;
  }
  return BigIntDivide(&tmp2, &tmp3, &pSum->denominator);
}

void BigRationalNegate(BigRational* pSrc, BigRational* pDest)
{
  (void)pDest;
  BigIntChSign(&pSrc->numerator);
}

enum eExprErr BigRationalDivide(BigRational* pDividend, BigRational* pDivisor, BigRational* pQuotient)
{
  enum eExprErr rc;
  rc = BigIntMultiply(&pDividend->numerator, &pDivisor->denominator, &tmp1);
  if (rc != EXPR_OK)
  {
    return rc;
  }
  rc = BigIntMultiply(&pDivisor->numerator, &pDividend->denominator, &tmp2);
  if (rc != EXPR_OK)
  {
    return rc;
  }
  // At this moment tmp1 is the numerator and tmp2 is the denominator.
  BigIntGcd(&tmp1, &tmp2, &tmp3);
  rc = BigIntDivide(&tmp1, &tmp3, &pQuotient->numerator);
  if (rc != EXPR_OK)
  {
    return rc;
  }
  return BigIntDivide(&tmp2, &tmp3, &pQuotient->denominator);
}

enum eExprErr BigRationalMultiply(BigRational* pFactor1, BigRational* pFactor2, BigRational* pProduct)
{
  enum eExprErr rc;
  rc = BigIntMultiply(&pFactor1->numerator, &pFactor2->numerator, &tmp1);
  if (rc != EXPR_OK)
  {
    return rc;
  }
  rc = BigIntMultiply(&pFactor1->denominator, &pFactor2->denominator, &tmp2);
  if (rc != EXPR_OK)
  {
    return rc;
  }
  // At this moment tmp1 is the numerator and tmp2 is the denominator.
  BigIntGcd(&tmp1, &tmp2, &tmp3);
  rc = BigIntDivide(&tmp1, &tmp3, &pProduct->numerator);
  if (rc != EXPR_OK)
  {
    return rc;
  }
  return BigIntDivide(&tmp2, &tmp3, &pProduct->denominator);
}

enum eExprErr BigRationalMultiplyByInt(BigRational* pFactor1, int factor2, BigRational* pProduct)
{
  multint(&tmp1, &pFactor1->numerator, factor2);
  BigIntGcd(&tmp1, &pFactor1->denominator, &tmp3);
  BigIntDivide(&tmp1, &tmp3, &pProduct->numerator);
  BigIntDivide(&pFactor1->denominator, &tmp3, &pProduct->denominator);
  return EXPR_OK;
}

enum eExprErr BigRationalDivideByInt(BigRational* pDividend, int divisor, BigRational* pQuotient)
{
  multint(&tmp1, &pDividend->denominator, divisor);
  BigIntGcd(&tmp1, &pDividend->numerator, &tmp3);
  BigIntDivide(&tmp1, &tmp3, &pQuotient->denominator);
  BigIntDivide(&pDividend->numerator, &tmp3, &pQuotient->numerator);
  return EXPR_OK;
}

void ForceDenominatorPositive(BigRational* rat)
{
  if (rat->denominator.sign == SIGN_NEGATIVE)
  {
    BigIntChSign(&rat->numerator);
    BigIntChSign(&rat->denominator);
  }
}

static void TryToSimplifyRatSqrt(BigInteger* RatPart, BigInteger* SqrPart, int divisor)
{
  while (getRemainder(SqrPart, divisor * divisor) == 0)
  {
    multint(RatPart, RatPart, divisor);
    subtractdivide(SqrPart, 0, divisor);
    subtractdivide(SqrPart, 0, divisor);
  }
}

void MultiplyRationalBySqrtRational(BigRational* RatPart, BigRational* SqrPart)
{
  int ctr;
  if (BigIntIsZero(&SqrPart->numerator))
  {
    intToBigInteger(&RatPart->numerator, 0);
    intToBigInteger(&RatPart->denominator, 1);
    intToBigInteger(&SqrPart->numerator, 0);
    intToBigInteger(&SqrPart->denominator, 1);
    return;
  }
  // Divide numerator and denominator by their gcd.
  BigIntGcd(&SqrPart->numerator, &SqrPart->denominator, &tmp1);
  BigIntDivide(&SqrPart->numerator, &tmp1, &SqrPart->numerator);
  BigIntDivide(&SqrPart->denominator, &tmp1, &SqrPart->denominator);
  // Get A = gcd(numerator of rational, denominator of sqrt)
  BigIntGcd(&RatPart->numerator, &SqrPart->denominator, &tmp1);
  // Get B = gcd(numerator of rational, denominator of sqrt / A)
  BigIntDivide(&SqrPart->denominator, &tmp1, &tmp2);
  BigIntGcd(&tmp1, &tmp2, &tmp3);
  // Compute numerator of rational <- numerator of rational / B
  BigIntDivide(&RatPart->numerator, &tmp3, &RatPart->numerator);
  // Compute denominator of sqrt <- denominator of sqrt / B^2
  BigIntDivide(&SqrPart->denominator, &tmp3, &SqrPart->denominator);
  BigIntDivide(&SqrPart->denominator, &tmp3, &SqrPart->denominator);
  // Get A = gcd(denominator of rational, numerator of sqrt)
  BigIntGcd(&RatPart->denominator, &SqrPart->numerator, &tmp1);
  // Get B = gcd(denominator of rational, numerator of sqrt / A)
  BigIntDivide(&SqrPart->numerator, &tmp1, &tmp2);
  BigIntGcd(&tmp1, &tmp2, &tmp3);
  // Compute denominator of rational <- denominator of rational / B
  BigIntDivide(&RatPart->denominator, &tmp3, &RatPart->denominator);
  // Compute numerator of sqrt <- numerator of sqrt / B^2
  BigIntDivide(&SqrPart->numerator, &tmp3, &SqrPart->numerator);
  BigIntDivide(&SqrPart->numerator, &tmp3, &SqrPart->numerator);
  // Let A = gcd(numerator of rational, denominator of sqrt)
  BigIntGcd(&RatPart->numerator, &SqrPart->denominator, &tmp1);
  // Divide numerator of rational by this gcd.
  BigIntDivide(&RatPart->numerator, &tmp1, &RatPart->numerator);
  // Divide denominator of sqrt by this gcd.
  BigIntDivide(&SqrPart->denominator, &tmp1, &SqrPart->denominator);
  // Multiply numerator of sqrt by the gcd.
  BigIntMultiply(&SqrPart->numerator, &tmp1, &SqrPart->numerator);
  ForceDenominatorPositive(SqrPart);
  ForceDenominatorPositive(RatPart);
  // 
  //
  // If numerator of SqrPart is a perfect square, multiply the square root by
  // the numerator of RatPart and replace the numerator of SqrPart by 1.
  squareRoot(SqrPart->numerator.limbs, tmp1.limbs, SqrPart->numerator.nbrLimbs, &tmp1.nbrLimbs);
  BigIntMultiply(&tmp1, &tmp1, &tmp2);
  BigIntSubt(&tmp2, &SqrPart->numerator, &tmp2);
  if (BigIntIsZero(&tmp2))
  {
    BigIntMultiply(&RatPart->numerator, &tmp1, &RatPart->numerator);
    intToBigInteger(&SqrPart->numerator, 1);
  }
  // If denominator of SqrPart is a perfect square, multiply the square root by
  // the denominator of RatPart and replace the denominator of SqrPart by 1.
  squareRoot(SqrPart->denominator.limbs, tmp1.limbs, SqrPart->denominator.nbrLimbs, &tmp1.nbrLimbs);
  BigIntMultiply(&tmp1, &tmp1, &tmp2);
  BigIntSubt(&tmp2, &SqrPart->denominator, &tmp2);
  if (BigIntIsZero(&tmp2))
  {
    BigIntMultiply(&RatPart->denominator, &tmp1, &RatPart->denominator);
    intToBigInteger(&SqrPart->denominator, 1);
  }
  // Divide numerator and denominator by small numbers.
  TryToSimplifyRatSqrt(&RatPart->numerator, &SqrPart->numerator, 2);
  TryToSimplifyRatSqrt(&RatPart->denominator, &SqrPart->denominator, 2);
  for (ctr = 3; ctr < 300; ctr += 2)
  {
    TryToSimplifyRatSqrt(&RatPart->numerator, &SqrPart->numerator, ctr);
    TryToSimplifyRatSqrt(&RatPart->denominator, &SqrPart->denominator, ctr);
  }
  // Divide numerator and denominator by their gcd.
  BigIntGcd(&RatPart->numerator, &RatPart->denominator, &tmp1);
  BigIntDivide(&RatPart->numerator, &tmp1, &RatPart->numerator);
  BigIntDivide(&RatPart->denominator, &tmp1, &RatPart->denominator);
}

bool BigRationalSquareRoot(BigRational* RatArgum, BigRational* RatSqRoot)
{
  if (BigIntIsZero(&RatArgum->numerator))
  {
    intToBigInteger(&RatSqRoot->numerator, 0);
    intToBigInteger(&RatSqRoot->denominator, 1);
    return true;
  }
  // Divide numerator and denominator by their gcd.
  BigIntGcd(&RatArgum->numerator, &RatArgum->denominator, &tmp1);
  BigIntDivide(&RatArgum->numerator, &tmp1, &RatArgum->numerator);
  BigIntDivide(&RatArgum->denominator, &tmp1, &RatArgum->denominator);
  // Process numerator.
  squareRoot(RatArgum->numerator.limbs, tmp1.limbs, RatArgum->numerator.nbrLimbs, &tmp1.nbrLimbs);
  BigIntMultiply(&tmp1, &tmp1, &tmp2);
  BigIntSubt(&tmp2, &RatArgum->numerator, &tmp2);
  if (!BigIntIsZero(&tmp2))
  {          // Numerator is not a perfect square.
    return false;
  }
  
  CopyBigInt(&RatSqRoot->numerator, &tmp1);
  // Process denominator.
  squareRoot(RatArgum->denominator.limbs, tmp1.limbs, RatArgum->denominator.nbrLimbs, &tmp1.nbrLimbs);
  BigIntMultiply(&tmp1, &tmp1, &tmp2);
  BigIntSubt(&tmp2, &RatArgum->denominator, &tmp2);
  if (!BigIntIsZero(&tmp2))
  {          // Denominator is not a perfect square.
    return false;
  }
  CopyBigInt(&RatSqRoot->denominator, &tmp1);
  return true;
}

static void showRationalPretty(BigRational* rat)
{
  if (pretty == PRETTY_PRINT)
  {
    showText("<span class=\"fraction\"><span class=\"numerator\">");
    shownbr(&rat->numerator);
    showText("</span><span class=\"denominator\">");
    shownbr(&rat->denominator);
    showText("</span></span>");
  }
  else
  {      // Tex
    showText("\\frac{");
    shownbr(&rat->numerator);
    showText("}{");
    shownbr(&rat->denominator);
    showText("}");
  }
}

void showRationalNoParen(BigRational* rat)
{
  bool denominatorIsNotOne = ((rat->denominator.nbrLimbs != 1) || (rat->denominator.limbs[0].x != 1));
  if ((pretty != PARI_GP) && denominatorIsNotOne)
  {
    showRationalPretty(rat);
    return;
  }
  shownbr(&rat->numerator);
  if (denominatorIsNotOne)
  {
    showText(" / ");
    shownbr(&rat->denominator);
  }
}

void showRationalOverStr(BigRational* rat, char *str, char *ptrTimes)
{
  bool denominatorIsNotOne = ((rat->denominator.nbrLimbs != 1) || (rat->denominator.limbs[0].x != 1));
  if (pretty != PARI_GP)
  {
    showText(pretty == PRETTY_PRINT?
      "<span class=\"fraction\"><span class=\"numerator\">": "\\frac{");
    shownbr(&rat->numerator);
    showText(pretty == PRETTY_PRINT ? "</span><span class=\"denominator\">": "}{");
    if (denominatorIsNotOne)
    {
      shownbr(&rat->denominator);
      showText(ptrTimes);
      showText(" ");
    }
    showText(str);
    showText(pretty == PRETTY_PRINT ? "</span></span>": "}");
    return;
  }
  shownbr(&rat->numerator);
  showText(" / ");
  if (denominatorIsNotOne)
  {
    showText("(");
    shownbr(&rat->denominator);
    showText(" * ");
    showText(str);
    showText(")");
  }
  else
  {
    showText(str);
  }
}

void showRational(BigRational* rat)
{
  bool denominatorIsNotOne = ((rat->denominator.nbrLimbs != 1) || (rat->denominator.limbs[0].x != 1));
  bool showParen;
  if ((pretty != PARI_GP) && denominatorIsNotOne)
  {
    showRationalPretty(rat);
    return;
  }
  showParen = denominatorIsNotOne || (rat->numerator.sign == SIGN_NEGATIVE);
  if (showParen)
  {
    *ptrOutput++ = '(';
  }
  shownbr(&rat->numerator);
  if (denominatorIsNotOne)
  {
    showText(" / ");
    shownbr(&rat->denominator);
  }
  if (showParen)
  {
    *ptrOutput++ = ')';
  }
}

void ShowRationalAndSqrParts(BigRational* RatPart, BigRational* SqrPart, int root,
  char *ptrTimes)
{
  if ((SqrPart->numerator.nbrLimbs != 1) || (SqrPart->numerator.limbs[0].x != 1) ||
    (SqrPart->denominator.nbrLimbs != 1) || (SqrPart->denominator.limbs[0].x != 1))
  {       // Square root part is not 1.
    if ((RatPart->numerator.nbrLimbs != 1) || (RatPart->numerator.limbs[0].x != 1) ||
      (RatPart->denominator.nbrLimbs != 1) || (RatPart->denominator.limbs[0].x != 1))
    {     // Absolute value of rational part is not 1.
      if (root == 2)
      {
        showRational(RatPart);
      }
      else
      {
        if (pretty == PRETTY_PRINT)
        {
          showText("<span class=\"root\"><span class=\"radicand2\">");
          showRational(RatPart);
          showText("</span></span>");
        }
        else if (pretty == TEX)
        {
          showText("\\sqrt{");
          showRational(RatPart);
          showText("}");
        }
        else
        {
          showRational(RatPart);
          showText("^(1/2)");
        }
      }
      if (pretty == 0)
      {
        *ptrOutput++ = ' ';
      }
      showText(ptrTimes);
      if (pretty == 0)
      {
        *ptrOutput++ = ' ';
      }
    }
    else
    {     // Absolute value of rational part is 1.
      if (RatPart->numerator.sign == SIGN_NEGATIVE)
      {   // Rational part is 1. Show negative sign.
        showText(pretty != PARI_GP? "&minus;": "-");
      }
    }
    if (pretty == PRETTY_PRINT)
    {
      if (root == 2)
      {
        showText("<span class=\"root\"><span class=\"radicand2\">");
      }
      else
      {
        showText("<span class=\"root\"><span class=\"befrad\" aria-hidden=\"true\">4</span><span class=\"radicand4\">");
      }
      showRational(SqrPart);
      showText("</span></span>");
    }
    else if (pretty == TEX)
    {
      if (root == 2)
      {
        showText("\\sqrt{");
      }
      else
      {
        showText("\\sqrt[4]{");
      }
      showRational(SqrPart);
      showText("}");
    }
    else
    {
      showRational(SqrPart);
      showText(root == 2 ? "^(1/2)" : "^(1/4)");
    }
  }
  else
  {
    showRational(RatPart);
  }
}

void showSquareRootOfRational(BigRational* rat, int root, char *ptrTimes)
{
  intToBigInteger(&Rat1.numerator, 1);
  intToBigInteger(&Rat1.denominator, 1);
  CopyBigInt(&Rat2.numerator, &rat->numerator);
  CopyBigInt(&Rat2.denominator, &rat->denominator);
  MultiplyRationalBySqrtRational(&Rat1, &Rat2);
  ShowRationalAndSqrParts(&Rat1, &Rat2, root, ptrTimes);
}
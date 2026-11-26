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
#include <string.h>
#include "string/strings.h"
#include "bignbr.h"
#include "expression.h"
#include "highlevel.h"
#include "showtime.h"
#include "copyStr.h"
#include "commonstruc.h"

#ifndef __ANDROID__
union uCommon common;
#endif

extern bool hexadecimal;
static void ShowRational(BigInteger* pNum, BigInteger* pDen);

static void showText(const char *text)
{
  copyStr(&common.misc.ptrOutput, text);
}

static void BigInteger2Out(char** ppDecimal, const BigInteger* pBigInt, int groupLength)
{
  if (hexadecimal)
  {
    BigInteger2Hex(ppDecimal, pBigInt, groupLength);
  }
  else
  {
    BigInteger2Dec(ppDecimal, pBigInt, groupLength);
  }
}

static void ShowConvergents(int coeffIndex, const BigInteger *coeff)
{
  CopyBigInt(&common.misc.U3, &common.misc.U2);                   // U3 <- U2, U2 <- U1, U1 <- a*U2 + U3
  CopyBigInt(&common.misc.U2, &common.misc.U1);
  (void)BigIntMultiply(coeff, &common.misc.U2, &common.misc.U1);
  BigIntAdd(&common.misc.U1, &common.misc.U3, &common.misc.U1);
  CopyBigInt(&common.misc.V3, &common.misc.V2);                   // V3 <- V2, V2 <- V1, V1 <- a*V2 + V3
  CopyBigInt(&common.misc.V2, &common.misc.V1);
  (void)BigIntMultiply(coeff, &common.misc.V2, &common.misc.V1);
  BigIntAdd(&common.misc.V1, &common.misc.V3, &common.misc.V1);
  formatString(&common.misc.ptrOutput, "A[$1d]/B[$1d] = $2b/$3b, a[$1d] = $4b<br>",
    coeffIndex, &common.misc.U1, &common.misc.V1, coeff);
  *common.misc.ptrOutput = 0;
}

// PQa algorithm for 
// (P+G)/Q where G is the square root of discriminant
// If D - U^2 is not multiple of V then 
//   U = U*V
//   V = V*V
//   G = G*V
// U1 <- 1, U2 <- 0
// V1 <- 0, V2 <- 1
// Perform loop:
// a = floor((U + G)/V)
// U3 <- U2, U2 <- U1, U1 <- a*U2 + U3
// V3 <- V2, V2 <- V1, V1 <- a*V2 + V3
// U <- a*V - U
// V <- (D - U^2)/V
// Inside period when: 0 <= G - U < V
static void PeriodicContinuedFraction(bool isContinuation)
{
  bool ended;
  size_t diffPtrs;
  if (!isContinuation)
  {
    // Check whether parameters have to be changed.
    // Initialize variables.
    intToBigInteger(&common.misc.startPeriodNum, -1);     // Less than zero means outside period.
    intToBigInteger(&common.misc.startPeriodDen, -1);
    common.misc.coefIndex = 0;
    common.misc.periodIndex = 0;
  }
  ended = false;
  do
  {
    BigIntAdd(&common.misc.num, &common.misc.intSqrt, &common.misc.bigTmp);
    if (common.misc.den.sign == SIGN_NEGATIVE)
    {   // If denominator is negative, round square root upwards.
      addbigint(&common.misc.bigTmp, 1);
    }
    floordiv(&common.misc.bigTmp, &common.misc.den, &common.misc.Temp);         // Temp = Term of continued fraction.
    if (common.misc.coefIndex == 0)
    {
      CopyBigInt(&common.misc.integerPart, &common.misc.Temp);
    }
    if (common.misc.isShowingConvergents)
    {      // Show convergent checkbox is checked.
      ShowConvergents(common.misc.coefIndex, &common.misc.Temp);
    }
    else
    {      // Show convergent checkbox not checked.
      BigInteger2Out(&common.misc.ptrOutput, &common.misc.Temp, groupLen);  // Show continued fraction coefficient.
      copyStr(&common.misc.ptrOutput, ((common.misc.coefIndex == 0) ? " + //" : ", ")); // Show separator.
    }
    (void)BigIntMultiply(&common.misc.Temp, &common.misc.den, &common.misc.bigTmp);  // U <- a*V - U
    BigIntSubt(&common.misc.bigTmp, &common.misc.num, &common.misc.num);
    (void)BigIntMultiply(&common.misc.num, &common.misc.num, &common.misc.bigTmp);   // V <- (D - U^2)/V
    BigIntSubt(&common.misc.discr, &common.misc.bigTmp, &common.misc.bigTmp);
    (void)BigIntDivide(&common.misc.bigTmp, &common.misc.den, &common.misc.Temp);
    CopyBigInt(&common.misc.den, &common.misc.Temp);
    common.misc.coefIndex++;
    if (common.misc.startPeriodNum.sign == SIGN_POSITIVE)
    {             // Already inside period.
      common.misc.periodIndex++;
      if (BigIntEqual(&common.misc.num, &common.misc.startPeriodNum) &&
        BigIntEqual(&common.misc.den, &common.misc.startPeriodDen))
      {           // New period started.
        ended = true;
        break;    // Go out in this case.
      }
    }
    else
    {             // Check if periodic part of continued fraction has started.
      BigIntSubt(&common.misc.intSqrt, &common.misc.num, &common.misc.bigTmp);  // G - U must be >= 0 and < V.
      if (common.misc.bigTmp.sign == SIGN_POSITIVE)
      {
        BigIntSubt(&common.misc.bigTmp, &common.misc.den, &common.misc.bigTmp);
        if (common.misc.bigTmp.sign == SIGN_NEGATIVE)
        {                                 // Periodic part of continued fraction started.
          CopyBigInt(&common.misc.startPeriodNum, &common.misc.num);       // Save U and V to check period end.
          CopyBigInt(&common.misc.startPeriodDen, &common.misc.den);
          showText("<span class=\"offscr\">");
          // start periodic part
          showText(LITERAL_START_PERIODIC_PART);
          showText("</span><span class=\"bold\">");
        }
      }
    }
    diffPtrs = common.misc.ptrOutput - &output[0];
  } while ((int)diffPtrs < ((int)sizeof(output) - 30000));
  if (!common.misc.isShowingConvergents)
  {        // Show convergent checkbox not checked.
    common.misc.ptrOutput -= 2;                       // Delete extra comma and space at the end.
  }
  if (ended)
  {
    showText("</span>");
    if (!common.misc.isShowingConvergents)
    {        // Show convergent checkbox not checked.
      showText("//");
    }
    showText("<br><span aria-hidden=\"true\">");
    // where the periodic part is marked in bold
    showText(LITERAL_PERIODIC_PART_IN_BOLD);
    showText("</span>");
    if (common.misc.periodIndex > 1)
    {
      // (the period has $1d coefficients)
      formatString(&common.misc.ptrOutput, LITERAL_COEFFICIENTS_IN_PERIOD, common.misc.periodIndex);
    }
  }
  else
  {  // Too many coefficients.
    if (!common.misc.isShowingConvergents)
    {        // Show convergents checkbox not checked.
      showText(", ... </span>//");
    }
    else
    {
      showText("... </span>");
    }
    showText("<br>");
    // where the periodic part (truncated after $1d coefficients) is marked in bold.
    formatString(&common.misc.ptrOutput, LITERAL_TRUNCATED_PERIOD, common.misc.periodIndex);
    output[0] = '6';         // Show Continue button.
  }
  showText("</p>");
}

static void showFooter(void)
{
  showElapsedTime(&common.misc.ptrOutput);
  copyStr(&common.misc.ptrOutput, "<p>");
  showCopyright(&common.misc.ptrOutput);
  copyStr(&common.misc.ptrOutput, "</p>");
}

static void showFormula(void)
{
  common.misc.ptrOutput = output;
  // Show formula.
  // the fraction whose numerator is
  // and the denominator is
  formatString(&common.misc.ptrOutput, "2<p><var>x</var> = <span class=\"fraction\">"
    "<span class=\"offscr\">$1s </span><span class=\"fup\">$2b" // Numerator
    " + <span class=\"sqrtout\"><span class=\"sqrtin\">$3b"     // Radicand
    "</span></span></span><span class=\"bar\"> </span><span class=\"fdn\">"
    "<span class=\"offscr\">$4s</span>$5b</span></span></p>",   // Denominator
    LITERAL_SHOW_FORMULA_NUM, &common.misc.origNum, &common.misc.origDelta, LITERAL_SHOW_FORMULA_DEN,
    &common.misc.origDen);
}

static void showX(void)
{
  showText("<p>");
  if (!common.misc.isShowingConvergents)
  {        // Show convergent checkbox not checked.
    showText("<span class = \"offscr\">");
    // The expansion in continued fraction of
    showText(LITERAL_EXPANSION_CONTFRAC);
    showText(" </span><var>x</var> = ");
  }
}

static void ContFrac(void)
{
  CopyBigInt(&common.misc.origNum, &common.misc.num);
  CopyBigInt(&common.misc.origDen, &common.misc.den);
  CopyBigInt(&common.misc.origDelta, &common.misc.discr);
  showFormula();
  // Validate input.
  if (BigIntIsZero(&common.misc.den))
  {
    // Error: The denominator is zero.
    formatString(&common.misc.ptrOutput, "<p>$1s</p>", LITERAL_ZERO_DENOMINATOR);
    return;
  }
  if (common.misc.discr.sign==SIGN_NEGATIVE)
  {   /* Complex number */
    // The number is not real, so it does not have continued fraction expansion.
    formatString(&common.misc.ptrOutput, "<p>$1s</p>", LITERAL_NUMBER_NOT_REAL);
    return;
  }
  showX();
  intToBigInteger(&common.misc.U1, 1);
  intToBigInteger(&common.misc.U2, 0);
  intToBigInteger(&common.misc.V1, 0);
  intToBigInteger(&common.misc.V2, 1);
  if (BigIntIsZero(&common.misc.discr))
  {   /* Rational number */
    ShowRational(&common.misc.num, &common.misc.den);
    return;
  }
  (void)BigIntMultiply(&common.misc.num, &common.misc.num, &common.misc.Temp);
  BigIntSubt(&common.misc.discr, &common.misc.Temp, &common.misc.Temp);  // Temp <- discr - num^2.
  (void)BigIntRemainder(&common.misc.Temp, &common.misc.den, &common.misc.Temp);
  if (!BigIntIsZero(&common.misc.Temp))
  {            // If discr - num^2 is not multiple of den...
    enum eSign sign;
    (void)BigIntMultiply(&common.misc.discr, &common.misc.den, &common.misc.discr);
    (void)BigIntMultiply(&common.misc.discr, &common.misc.den, &common.misc.discr);   // discr <- discr*den^2
    sign = common.misc.den.sign;
    (void)BigIntMultiply(&common.misc.num, &common.misc.den, &common.misc.num);       // num <- num*den
    (void)BigIntMultiply(&common.misc.den, &common.misc.den, &common.misc.den);       // den <- den*den
    if (sign == SIGN_NEGATIVE)
    {
      BigIntNegate(&common.misc.num, &common.misc.num);             // Change sign to both
      BigIntNegate(&common.misc.den, &common.misc.den);             // numerator and denominator.
    }
  }
  squareRoot(common.misc.discr.limbs, common.misc.intSqrt.limbs,
    common.misc.discr.nbrLimbs, &common.misc.intSqrt.nbrLimbs);
  common.misc.intSqrt.sign = SIGN_POSITIVE;
  (void)BigIntMultiply(&common.misc.intSqrt, &common.misc.intSqrt, &common.misc.Temp);
  if (TestBigNbrEqual(&common.misc.Temp, &common.misc.discr) != 0)
  {     // discr is a perfect square, so number is rational.
    BigIntAdd(&common.misc.num, &common.misc.intSqrt, &common.misc.Temp);
    ShowRational(&common.misc.Temp, &common.misc.den);
  }
  else
  {     // discr is not a perfect square. Periodic continued fraction.
    PeriodicContinuedFraction(false);   // Start showing periodic continued fraction.
  }
  showFooter();
}

static void ShowRational(BigInteger *pNum, BigInteger *pDen)
{
  static BigInteger Tmp;
  const char *sep;
  int coeffIndex = 0;

  BigIntGcd(pNum, pDen, &Tmp);     // Tmp <- GCD of numerator and denominator.
  (void)BigIntDivide(pNum, &Tmp, pNum);
  (void)BigIntDivide(pDen, &Tmp, pDen);
  floordiv(pNum, pDen, &Tmp);
  (void)BigIntMultiply(pDen, &Tmp, &common.misc.bigTmp);
  BigIntSubt(pNum, &common.misc.bigTmp, pNum);   // Reduce numerator.
  if (common.misc.isShowingConvergents)
  {      // Show convergent checkbox is checked.
    ShowConvergents(coeffIndex, &Tmp);
  }
  else
  {
    BigInteger2Out(&common.misc.ptrOutput, &Tmp, groupLen);  // Show convergent.
  }
  (void)BigIntRemainder(pNum, pDen, pNum);
  sep = " + //";
  while (!BigIntIsZero(pNum))
  {      // Numerator greater than zero.
    coeffIndex++;
    if (!BigIntIsZero(pDen))
    {    // Denominator greater than zero.
      floordiv(pDen, pNum, &Tmp);
      if (common.misc.isShowingConvergents)
      {      // Show convergent checkbox is checked.
        ShowConvergents(coeffIndex, &Tmp);
      }
      else
      {
        showText(sep);
        BigInteger2Out(&common.misc.ptrOutput, &Tmp, groupLen);  // Show convergent.
      }
      sep = ", ";
    }
    (void)BigIntMultiply(pNum, &Tmp, &common.misc.bigTmp);
    BigIntSubt(pDen, &common.misc.bigTmp, &Tmp);
    CopyBigInt(pDen, pNum);
    CopyBigInt(pNum, &Tmp);
  }
  if (!common.misc.isShowingConvergents && (sep[0] == ','))
  {         // Inside continued fraction. Close it.
    showText("//");
  }
}

static void getNumber(BigInteger *pNumber, const char *title, char** pptrInput)
{
  enum eExprErr rc;
  rc = ComputeExpression(*pptrInput, pNumber);
  if (rc != EXPR_OK)
  {
    if (common.misc.ptrOutput == output)
    {
      *common.misc.ptrOutput = '2';
      common.misc.ptrOutput++;
    }
    showText("<p>");
    showText(title);
    *common.misc.ptrOutput = ':';
    common.misc.ptrOutput++;
    *common.misc.ptrOutput = ' ';
    common.misc.ptrOutput++;
    textError(&common.misc.ptrOutput, rc);
    showText("</p>");
  }
  *pptrInput += strlen(*pptrInput) + 1U;  // Skip terminator.
}

// input contains three expressions separated by 00h (null character).
void contfracText(char *input, int GroupLen, bool hex, bool converg)
{
  if (input == NULL)
  {
    showFormula();
    if (!common.misc.isShowingConvergents)
    {        // Show convergents checkbox not checked.
      showX();
      BigInteger2Out(&common.misc.ptrOutput, &common.misc.integerPart, groupLen);  // Show integer part.
      copyStr(&common.misc.ptrOutput, " + // ");       // Show separator.
    }
    showText("<span class=\"offscr\">");   // Indicate start of periodic part.
    // start periodic part
    showText(LITERAL_START_PERIODIC_PART);
    showText("</span><span class=\"bold\">");
    copyStr(&common.misc.ptrOutput, "...");            // Show continuation.
    PeriodicContinuedFraction(true);
    showFooter();
    return;
  }
  char *ptrInput = input;
  groupLen = GroupLen;
  common.misc.ptrOutput = output;
  common.misc.isShowingConvergents = converg;
  hexadecimal = hex;
  getNumber(&common.misc.num, LITERAL_NUMERATOR, &ptrInput);    // Numerator
  getNumber(&common.misc.discr, LITERAL_ARGUMENT_SQROOT, &ptrInput);  // Square root argument
  getNumber(&common.misc.den, LITERAL_DENOMINATOR, &ptrInput);  // Denominator
  if (common.misc.ptrOutput != output)
  {       // Error written to output. Do not compute continued fraction.
    return;
  }
  groupLen = GroupLen;
  ContFrac();
}

#if defined(__EMSCRIPTEN__) && !defined(_MSC_VER)
#ifdef __ANDROID__
EXTERNALIZE void doWorkContFrac(void)
#else
EXTERNALIZE void doWork(void)
#endif
{
  int app;
  bool hex;
  bool converg;
  int grpLen = 0;
  char* ptrData = inputString;
  if (*ptrData == 'C')
  {    // User pressed Continue button.
    contfracText(NULL, 0, false, false); // Routine does not use parameters in this case.
    databack(output);
    return;
  }
  originalTenthSecond = tenths();
  while (*ptrData != ',')
  {
    grpLen = (grpLen * 10) + (*ptrData - '0');
    ptrData++;
  }
  ptrData++;             // Skip comma.
  app = *ptrData - '0';
  while (*(ptrData + 1) != ',')
  {
    ptrData++;
    app = (app * 10) + *ptrData - '0';
  }
  hex = ((app & 0x40) != 0);
  groupSeparator = app >> 8;
  converg = ((app & 0x20) != 0);
  contfracText(ptrData + 2, grpLen, hex, converg);
  databack(output);
}
#endif

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
#include "bignbr.h"
#include "expression.h"
#include "highlevel.h"
static BigInteger num;
static BigInteger den;
static BigInteger delta;
static BigInteger Temp;
static BigInteger bigTmp;
static BigInteger U1;
static BigInteger U2;
static BigInteger U3;
static BigInteger V1;
static BigInteger V2;
static BigInteger V3;
static BigInteger startPeriodNum;
static BigInteger startPeriodDen;
static char *ptrOutput;
static void ShowRational(BigInteger *pNum, BigInteger *pDen);
extern bool hexadecimal;
static void showText(const char *text)
{
  copyStr(&ptrOutput, text);
}

static void ShowConvergents(int index, const BigInteger *coeff)
{
  CopyBigInt(&U3, &U2);                   // U3 <- U2, U2 <- U1, U1 <- a*U2 + U3
  CopyBigInt(&U2, &U1);
  (void)BigIntMultiply(coeff, &U2, &U1);
  BigIntAdd(&U1, &U3, &U1);
  CopyBigInt(&V3, &V2);                   // V3 <- V2, V2 <- V1, V1 <- a*V2 + V3
  CopyBigInt(&V2, &V1);
  (void)BigIntMultiply(coeff, &V2, &V1);
  BigIntAdd(&V1, &V3, &V1);
  *ptrOutput = 'A';
  ptrOutput++;
  *ptrOutput = '[';
  ptrOutput++;
  int2dec(&ptrOutput, index);
  *ptrOutput = ']';
  ptrOutput++;
  *ptrOutput = '/';
  ptrOutput++;
  *ptrOutput = 'B';
  ptrOutput++;
  *ptrOutput = '[';
  ptrOutput++;
  int2dec(&ptrOutput, index);
  *ptrOutput = ']';
  ptrOutput++;
  *ptrOutput = ' ';
  ptrOutput++;
  *ptrOutput = '=';
  ptrOutput++;
  *ptrOutput = ' ';
  ptrOutput++;
  BigInteger2Dec(&ptrOutput, &U1, groupLen);  // Show continued fraction convergent.
  *ptrOutput = '/';
  ptrOutput++;
  BigInteger2Dec(&ptrOutput, &V1, groupLen);  // Show continued fraction convergent.
  *ptrOutput = ',';
  ptrOutput++;
  *ptrOutput = ' ';
  ptrOutput++;
  *ptrOutput = 'a';
  ptrOutput++;
  *ptrOutput = '[';
  ptrOutput++;
  int2dec(&ptrOutput, index);
  *ptrOutput = ']';
  ptrOutput++;
  *ptrOutput = ' ';
  ptrOutput++;
  *ptrOutput = '=';
  ptrOutput++;
  *ptrOutput = ' ';
  ptrOutput++;
  BigInteger2Dec(&ptrOutput, coeff, groupLen);  // Show continued fraction coefficient.
  *ptrOutput = '<';
  ptrOutput++;
  *ptrOutput = 'b';
  ptrOutput++;
  *ptrOutput = 'r';
  ptrOutput++;
  *ptrOutput = '>';
  ptrOutput++;
  *ptrOutput = 0;
}

static void ContFrac(void)
{
  static BigInteger intSqrt;

  ptrOutput = output;
  // Show formula.
  showText("2<p><var>x</var> = <span class=\"fraction\"><span class=\"offscr\">");
  copyStr(&ptrOutput, lang ? " la fracción cuyo numerador es </span>" : " the fraction whose numerator is </span>");
  showText("<span class=\"fup\">");
  BigInteger2Dec(&ptrOutput, &num, groupLen);    // Show numerator.
  showText(" + <span class=\"sqrtout\"><span class=\"sqrtin\">");
  BigInteger2Dec(&ptrOutput, &delta, groupLen);  // Show radicand.
  showText("</span></span></span><span class=\"bar\"> </span><span class=\"fdn\"><span class=\"offscr\">");
  copyStr(&ptrOutput, lang ? " y el denominador es </span>" : " and the denominator is </span>");
  BigInteger2Dec(&ptrOutput, &den, groupLen);    // Show denominator.
  showText("</span></span></span></p>");
  // Validate input.
  if (BigIntIsZero(&den))
  {
    showText(lang? "<p>Error: El denominador es cero.</p>": "<p>Error: The denominator is zero.</p>");
    return;
  }
  if (delta.sign==SIGN_NEGATIVE)
  {   /* Complex number */
    showText(lang? "<p>El número no es real, por lo que no tiene desarrollo en fracciones continuas.</p>":
                   "<p>The number is not real, so it does not have continued fraction expansion.</p>");
    return;
  }
  showText("<p>");
  if (!hexadecimal)
  {        // Show convergent checkbox not checked.
    showText("<span class = \"offscr\">");
    showText(lang ? "El desarrollo en fracción continua de" :
      "The expansion in continued fraction of");
    showText(" </span><var>x</var> = ");
  }
  intToBigInteger(&U1, 1);
  intToBigInteger(&U2, 0);
  intToBigInteger(&V1, 0);
  intToBigInteger(&V2, 1);
  if (BigIntIsZero(&delta))
  {   /* Rational number */
    ShowRational(&num, &den);
    return;
  }
  (void)BigIntMultiply(&num, &num, &Temp);
  BigIntSubt(&delta, &Temp, &Temp);  // Temp <- delta - num^2.
  (void)BigIntRemainder(&Temp, &den, &Temp);
  if (!BigIntIsZero(&Temp))
  {            // If delta - num^2 is not multiple of den...
    enum eSign sign;
    (void)BigIntMultiply(&delta, &den, &delta);
    (void)BigIntMultiply(&delta, &den, &delta);   // delta <- delta*den^2
    sign = den.sign;
    (void)BigIntMultiply(&num, &den, &num);       // num <- num*den
    (void)BigIntMultiply(&den, &den, &den);       // den <- den*den
    if (sign == SIGN_NEGATIVE)
    {
      BigIntNegate(&num, &num);             // Change sign to both
      BigIntNegate(&den, &den);             // numerator and denominator.
    }
  }
  squareRoot(delta.limbs, intSqrt.limbs, delta.nbrLimbs, &intSqrt.nbrLimbs);
  intSqrt.sign = SIGN_POSITIVE;
  (void)BigIntMultiply(&intSqrt, &intSqrt, &Temp);
  if (TestBigNbrEqual(&Temp, &delta) != 0)
  {     // delta is a perfect square, so number is rational.
    BigIntAdd(&num, &intSqrt, &Temp);
    ShowRational(&Temp, &den);
  }
  else
  {     // delta is not a perfect square. Periodic continued fraction.
    bool ended;
    int periodIndex;
    int index;
    size_t diffPtrs;
    // PQa algorithm for (P+G)/Q where G = sqrt(discriminant)
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
        // Check whether parameters have to be changed.
        // Initialize variables.
    intToBigInteger(&startPeriodNum, -1);     // Less than zero means outside period.
    intToBigInteger(&startPeriodDen, -1);
    index = 0;
    periodIndex = 0;
    ended = false;
    do
    {
      BigIntAdd(&num, &intSqrt, &bigTmp);
      if (den.sign == SIGN_NEGATIVE)
      {   // If denominator is negative, round square root upwards.
        addbigint(&bigTmp, 1);
      }      
      floordiv(&bigTmp, &den, &Temp);         // Temp = Term of continued fraction.
      if (!hexadecimal)
      {      // Show convergent checkbox not checked.
        BigInteger2Dec(&ptrOutput, &Temp, groupLen);  // Show continued fraction coefficient.
        copyStr(&ptrOutput, ((index == 0)? " + //" : ", ")); // Show separator.
      }
      if (hexadecimal)
      {      // Show convergent checkbox is checked.
        ShowConvergents(index, &Temp);
      }
      (void)BigIntMultiply(&Temp, &den, &bigTmp);  // U <- a*V - U
      BigIntSubt(&bigTmp, &num, &num);
      (void)BigIntMultiply(&num, &num, &bigTmp);   // V <- (D - U^2)/V
      BigIntSubt(&delta, &bigTmp, &bigTmp);
      (void)BigIntDivide(&bigTmp, &den, &Temp);
      CopyBigInt(&den, &Temp);
      index++;
      if (startPeriodNum.sign == SIGN_POSITIVE)
      {             // Already inside period.
        periodIndex++;
        if (BigIntEqual(&num, &startPeriodNum) && BigIntEqual(&den, &startPeriodDen))
        {           // New period started.
          ended = true;
          break;    // Go out in this case.
        }
      }
      else
      {             // Check if periodic part of continued fraction has started.
        BigIntSubt(&intSqrt, &num, &bigTmp);  // G - U must be >= 0 and < V.
        if (bigTmp.sign == SIGN_POSITIVE)
        {
          BigIntSubt(&bigTmp, &den, &bigTmp);
          if (bigTmp.sign == SIGN_NEGATIVE)
          {                                 // Periodic part of continued fraction started.
            CopyBigInt(&startPeriodNum, &num);       // Save U and V to check period end.
            CopyBigInt(&startPeriodDen, &den);
            showText("<span class=\"offscr\">");
            showText(lang ? "inicio del período" : "start periodic part");
            showText("</span><span class=\"bold\">");
          }
        }
      }
      diffPtrs = ptrOutput - &output[0];
    } while ((int)diffPtrs < ((int)sizeof(output) - 30000));
    if (!hexadecimal)
    {        // Show convergent checkbox not checked.
      ptrOutput -= 2;                       // Delete extra comma and space at the end.
    }
    if (ended)
    {
      showText("</span>");
      if (!hexadecimal)
      {        // Show convergent checkbox not checked.
        showText("//");
      }
      showText("<br /><span aria-hidden=\"true\">");
      showText(lang? "donde la parte periódica está señalada en negrita</span>" :
        "where the periodic part is marked in bold</span>");
      if (periodIndex > 1)
      {
        showText(lang? " (el período tiene " : " (the period has ");
        int2dec(&ptrOutput, periodIndex);
        showText(lang? " coeficientes)" : " coefficients)");
      }
    }
    else
    {  // Too many convergents.
      if (!hexadecimal)
      {        // Show convergent checkbox not checked.
        showText(", ... </span>//");
      }
      else
      {
        showText("... </span>");
      }
      showText(lang? "<br />donde la parte periódica (truncada a partir de los " :
        "//<br />where the periodic part (truncated after ");
      int2dec(&ptrOutput, periodIndex);
      showText(lang? " convergentes) está señalada en negrita.</p>" :
        " convergents) is marked in bold.</p>");
    }
  }
  showText(lang? "<p>" COPYRIGHT_SPANISH "</p>" :
                 "<p>" COPYRIGHT_ENGLISH "</p>");
}

static void ShowRational(BigInteger *pNum, BigInteger *pDen)
{
  static BigInteger Tmp;
  const char *sep;
  int index = 0;

  BigIntGcd(pNum, pDen, &Tmp);     // Tmp <- GCD of numerator and denominator.
  (void)BigIntDivide(pNum, &Tmp, pNum);
  (void)BigIntDivide(pDen, &Tmp, pDen);
  floordiv(pNum, pDen, &Tmp);
  (void)BigIntMultiply(pDen, &Tmp, &bigTmp);
  BigIntSubt(pNum, &bigTmp, pNum);   // Reduce numerator.
  if (hexadecimal)
  {      // Show convergent checkbox is checked.
    ShowConvergents(index, &Tmp);
  }
  else
  {
    BigInteger2Dec(&ptrOutput, &Tmp, groupLen);  // Show convergent.
  }
  (void)BigIntRemainder(pNum, pDen, pNum);
  sep = " + //";
  while (!BigIntIsZero(pNum))
  {      // Numerator greater than zero.
    index++;
    if (!BigIntIsZero(pDen))
    {    // Denominator greater than zero.
      floordiv(pDen, pNum, &Tmp);
      if (hexadecimal)
      {      // Show convergent checkbox is checked.
        ShowConvergents(index, &Tmp);
      }
      else
      {
        showText(sep);
        BigInteger2Dec(&ptrOutput, &Tmp, groupLen);  // Show convergent.
      }
      sep = ", ";
    }
    (void)BigIntMultiply(pNum, &Tmp, &bigTmp);
    BigIntSubt(pDen, &bigTmp, &Tmp);
    CopyBigInt(pDen, pNum);
    CopyBigInt(pNum, &Tmp);
  }
  if (!hexadecimal && (sep[0] == ','))
  {         // Inside continued fraction. Close it.
    showText("//");
  }
}

static int getNumber(BigInteger *pNumber, const char *title, char **pptrInput)
{
  enum eExprErr rc;
  rc = ComputeExpression(*pptrInput, pNumber);
  if (rc != EXPR_OK)
  {
    showText(title);
    *ptrOutput = ':';
    ptrOutput++;
    *ptrOutput = ' ';
    ptrOutput++;
    textError(&ptrOutput, rc);
    return 1;
  }
  *pptrInput += strlen(*pptrInput) + 1U;  // Skip terminator.
  return 0;
}

// input contains three expressions separated by 00h (null character).
void contfracText(char *input, int GroupLen)
{
  char *ptrInput = input;
  ptrOutput = output;
  if (getNumber(&num, lang? "Numerador": "Numerator", &ptrInput) != 0)
  {
    return;
  }
  if (getNumber(&delta, lang?"Argumento de la raíz cuadrada": "Square root argument", &ptrInput) != 0)
  {
    return;
  }
  if (getNumber(&den, lang?"Denominador": "Denominator", &ptrInput) != 0)
  {
    return;
  }
  groupLen = GroupLen;
  ContFrac();
}
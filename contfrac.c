/*
This file is part of Alpertron Calculators.

Copyright 2015 Dario Alejandro Alpern

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
#include "bignbr.h"
#include "expression.h"
#include "highlevel.h"
static BigInteger num, den, delta;
static BigInteger Temp, bigTmp;
static BigInteger U1, U2, U3;
static BigInteger V1, V2, V3;
static BigInteger startPeriodNum, startPeriodDen;
static char *ptrOutput;
static void ShowRational(BigInteger *pNum, BigInteger *pDen);
extern char hexadecimal;
static void showText(char *text)
{
  strcpy(ptrOutput, text);
  ptrOutput += strlen(ptrOutput);
}

static void ShowConvergents(int index, BigInteger *coeff)
{
  CopyBigInt(&U3, &U2);                   // U3 <- U2, U2 <- U1, U1 <- a*U2 + U3
  CopyBigInt(&U2, &U1);
  BigIntMultiply(coeff, &U2, &U1);
  BigIntAdd(&U1, &U3, &U1);
  CopyBigInt(&V3, &V2);                   // V3 <- V2, V2 <- V1, V1 <- a*V2 + V3
  CopyBigInt(&V2, &V1);
  BigIntMultiply(coeff, &V2, &V1);
  BigIntAdd(&V1, &V3, &V1);
  *ptrOutput++ = 'A';
  *ptrOutput++ = '[';
  int2dec(&ptrOutput, index);
  *ptrOutput++ = ']';
  *ptrOutput++ = '/';
  *ptrOutput++ = 'B';
  *ptrOutput++ = '[';
  int2dec(&ptrOutput, index);
  *ptrOutput++ = ']';
  *ptrOutput++ = ' ';
  *ptrOutput++ = '=';
  *ptrOutput++ = ' ';
  BigInteger2Dec(&U1, ptrOutput, groupLen);  // Show continued fraction convergent.
  ptrOutput += strlen(ptrOutput);
  *ptrOutput++ = '/';
  BigInteger2Dec(&V1, ptrOutput, groupLen);  // Show continued fraction convergent.
  ptrOutput += strlen(ptrOutput);
  *ptrOutput++ = ',';
  *ptrOutput++ = ' ';
  *ptrOutput++ = 'a';
  *ptrOutput++ = '[';
  int2dec(&ptrOutput, index);
  *ptrOutput++ = ']';
  *ptrOutput++ = ' ';
  *ptrOutput++ = '=';
  *ptrOutput++ = ' ';
  BigInteger2Dec(coeff, ptrOutput, groupLen);  // Show continued fraction coefficient.
  ptrOutput += strlen(ptrOutput);
  *ptrOutput++ = '<';
  *ptrOutput++ = 'b';
  *ptrOutput++ = 'r';
  *ptrOutput++ = '>';
  *ptrOutput = 0;
}

static void ContFrac(void)
{
  BigInteger intSqrt;

  ptrOutput = output;
  // Show formula.
  showText("2<p><var>x</var> = <span class=\"fraction\"><span class=\"offscr\">");
  strcpy(ptrOutput, lang ? " la fracción cuyo numerador es </span>" : " the fraction whose numerator is </span>");
  ptrOutput += strlen(ptrOutput);
  showText("<span class=\"fup\">");
  BigInteger2Dec(&num, ptrOutput, groupLen);    // Show numerator.
  ptrOutput += strlen(ptrOutput);
//  strcpy(ptrOutput, lang ? " más la raíz cuadrada de " : " plus the square root of ");
//  ptrOutput += strlen(ptrOutput);
  showText(" + <span class=\"sqrtout\"><span class=\"sqrtin\">");
  BigInteger2Dec(&delta, ptrOutput, groupLen);  // Show radicand.
  ptrOutput += strlen(ptrOutput);
  showText("</span></span></span><span class=\"bar\"> </span><span class=\"fdn\"><span class=\"offscr\">");
  strcpy(ptrOutput, lang ? " y el denominador es </span>" : " and the denominator is </span>");
  ptrOutput += strlen(ptrOutput);
  BigInteger2Dec(&den, ptrOutput, groupLen);    // Show denominator.
  ptrOutput += strlen(ptrOutput);
  showText("</span></span></span></p>");
  // Validate input.
  if (den.nbrLimbs==1 && den.limbs[0].x==0)
  {
    showText(lang != 0? "<p>Error: El denominador es cero.</p>": "<p>Error: The denominator is zero.</p>");
    return;
  }
  if (delta.sign==SIGN_NEGATIVE)
  {   /* Complex number */
    showText(lang != 0? "<p>El número no es real, por lo que no tiene desarrollo en fracciones continuas.</p>":
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
  if (delta.nbrLimbs == 1 && delta.limbs[0].x == 0)
  {   /* Rational number */
    ShowRational(&num, &den);
    return;
  }
  (void)BigIntMultiply(&num, &num, &Temp);
  BigIntSubt(&delta, &Temp, &Temp);  // Temp <- delta - num^2.
  (void)BigIntRemainder(&Temp, &den, &Temp);
  if (Temp.nbrLimbs != 1 || Temp.limbs[0].x != 0)
  {            // If delta - num^2 is not multiple of den...
    int sign;
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
    char ended;
    int periodIndex, index;
        // PQa algorithm for (P+G)/Q where G = sqrt(discriminant):
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
    ended = 0;
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
        BigInteger2Dec(&Temp, ptrOutput, groupLen);  // Show continued fraction coefficient.
        ptrOutput += strlen(ptrOutput);
        strcpy(ptrOutput, (index == 0 ? " + //" : ", ")); // Show separator.
        ptrOutput += strlen(ptrOutput);
      }
      if (hexadecimal)
      {      // Show convergent checkbox is checked.
        ShowConvergents(index, &Temp);
      }
      BigIntMultiply(&Temp, &den, &bigTmp);  // U <- a*V - U
      BigIntSubt(&bigTmp, &num, &num);
      BigIntMultiply(&num, &num, &bigTmp);   // V <- (D - U^2)/V
      BigIntSubt(&delta, &bigTmp, &bigTmp);
      BigIntDivide(&bigTmp, &den, &Temp);
      CopyBigInt(&den, &Temp);
      index++;
      if (startPeriodNum.sign == SIGN_POSITIVE)
      {             // Already inside period.
        if (++periodIndex == 100000)
        {
          break;    // Too many coefficients. Go out.
        }
        if (BigIntEqual(&num, &startPeriodNum) && BigIntEqual(&den, &startPeriodDen))
        {           // New period started.
          ended = 1;
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
    } while (ptrOutput - &output[0] < sizeof(output) - 30000);
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
      showText(lang != 0 ? "donde la parte periódica está señalada en negrita</span>" :
        "where the periodic part is marked in bold</span>");
      if (periodIndex > 1)
      {
        showText(lang != 0 ? " (el período tiene " : " (the period has ");
        int2dec(&ptrOutput, periodIndex);
        showText(lang != 0 ? " coeficientes)" : " coefficients)");
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
      showText(lang != 0 ? "<br />donde la parte periódica (truncada a partir de los " :
        "//<br />where the periodic part (truncated after ");
      int2dec(&ptrOutput, periodIndex);
      showText(lang != 0 ? " convergentes) está señalada en negrita.</p>" :
        " convergents) is marked in bold.</p>");
    }
  }
  showText(lang? "<p>" COPYRIGHT_SPANISH "</p>" :
                 "<p>" COPYRIGHT_ENGLISH "</p>");
}

static void ShowRational(BigInteger *pNum, BigInteger *pDen)
{
  BigInteger Tmp;
  char *sep;
  int index = 0;

  BigIntGcd(pNum, pDen, &Tmp);     // Tmp <- GCD of numerator and denominator.
  (void)BigIntDivide(pNum, &Tmp, pNum);
  (void)BigIntDivide(pDen, &Tmp, pDen);
  floordiv(pNum, pDen, &Tmp);
  BigIntMultiply(pDen, &Tmp, &bigTmp);
  BigIntSubt(pNum, &bigTmp, pNum);   // Reduce numerator.
  if (hexadecimal)
  {      // Show convergent checkbox is checked.
    ShowConvergents(index, &Tmp);
  }
  else
  {
    BigInteger2Dec(&Tmp, ptrOutput, groupLen);  // Show convergent.
    ptrOutput += strlen(ptrOutput);
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
        BigInteger2Dec(&Tmp, ptrOutput, groupLen);  // Show convergent.
        ptrOutput += strlen(ptrOutput);
      }
      sep = ", ";
    }
    BigIntMultiply(pNum, &Tmp, &bigTmp);
    BigIntSubt(pDen, &bigTmp, &Tmp);
    CopyBigInt(pDen, pNum);
    CopyBigInt(pNum, &Tmp);
  }
  if (!hexadecimal && sep[0] == ',')
  {         // Inside continued fraction. Close it.
    showText("//");
  }
}

static int getNumber(BigInteger *pNumber, char *title, char **pptrInput)
{
  enum eExprErr rc;
  rc = ComputeExpression(*pptrInput, 1, pNumber);
  if (rc != EXPR_OK)
  {
    showText(title);
    *ptrOutput++ = ':';
    *ptrOutput++ = ' ';
    textError(output, rc);
    return 1;
  }
  *pptrInput += strlen(*pptrInput) + 1;  // Skip terminator.
  return 0;
}

// input contains three expressions separated by 00h (null character).
void contfracText(char *input, int GroupLen)
{
  char *ptrInput = input;
  ptrOutput = output;
  if (getNumber(&num, lang != 0? "Numerador": "Numerator", &ptrInput) != 0)
  {
    return;
  }
  if (getNumber(&delta, lang != 0 ?"Argumento de la raíz cuadrada": "Square root argument", &ptrInput) != 0)
  {
    return;
  }
  if (getNumber(&den, lang != 0 ?"Denominador": "Denominator", &ptrInput) != 0)
  {
    return;
  }
  groupLen = GroupLen;
  ContFrac();
}
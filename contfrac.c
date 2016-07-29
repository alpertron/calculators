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
static BigInteger Temp;
extern char *output;
static char *ptrOutput;
static void ShowRational(BigInteger *pNum, BigInteger *pDen);
static int groupLen;
extern int lang;

static void showText(char *text)
{
  strcpy(ptrOutput, text);
  ptrOutput += strlen(ptrOutput);
}

static void ContFrac(void)
{
  BigInteger intSqrt;
  BigInteger biK, biP, biM;
  BigInteger K, L, M, P, Z;

  ptrOutput = output;
  if (den.nbrLimbs==1 && den.limbs[0].x==0)
  {
    showText(lang != 0? "Error: El denominador es cero.": "Error: The denominator is zero.");
    return;
  }
  if (delta.sign==SIGN_NEGATIVE)
  {   /* Complex number */
    showText(lang != 0? "El número no es real, por lo que no tiene desarrollo en fracciones continuas":
                   "The number is not real, so it does not have continued fraction expansion.");
    return;
  }
  showText("x = ");
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
  {
    int cont;
    char *sep;

    CopyBigInt(&biP, &den);
    CopyBigInt(&biK, &intSqrt);
    if (biP.sign == SIGN_NEGATIVE)
    {
      addbigint(&biK, 1);
    }
    BigIntAdd(&biK, &num, &biK);
    if (biK.sign == SIGN_POSITIVE)
    {
      if (den.sign == SIGN_POSITIVE)
      {
        (void)BigIntDivide(&biK, &den, &biM);
      }
      else
      {
        CopyBigInt(&biM, &den);
        addbigint(&biM, 1);
        BigIntSubt(&biM, &biK, &biM);
        BigIntNegate(&den, &den);
        (void)BigIntDivide(&biM, &den, &biM);   // biM <- ((den+1)-biK)/(-den)
        BigIntNegate(&den, &den);
      }
    }
    else
    {
      if (den.sign == SIGN_POSITIVE)
      {
        CopyBigInt(&biM, &biK);
        addbigint(&biM, 1);
        BigIntSubt(&biM, &den, &biM);
        (void)BigIntDivide(&biM, &den, &biM);   // biM <- ((biK+1)-den)/den
      }
      else
      {
        BigIntNegate(&biK, &biM);
        BigIntNegate(&den, &den);
        (void)BigIntDivide(&biM, &den, &biM);   // biM <- (-biK)/(-den)
        BigIntNegate(&den, &den);
      }
    }

    BigInteger2Dec(&biM, ptrOutput, groupLen);  // Show integer part.
    ptrOutput += strlen(ptrOutput);
    (void)BigIntMultiply(&biM, &den, &biM);
    BigIntSubt(&biM, &num, &biM);
    sep = " + //";
    cont = -1;
    K.sign = L.sign = P.sign = M.sign = SIGN_NEGATIVE;
    K.nbrLimbs = L.nbrLimbs = P.nbrLimbs = M.nbrLimbs = 1;
    K.limbs[0].x = L.limbs[0].x = P.limbs[0].x = M.limbs[0].x = 1;   // K <- L <- P <- M <- -1. 
    while ((cont < 0 || TestBigNbrEqual(&K, &P) == 0|| TestBigNbrEqual(&L, &M)==0) && cont < 100000)
    {
      showText(sep);
      sep = ", ";
      if (cont < 0 && biP.sign == SIGN_POSITIVE && biM.sign == SIGN_POSITIVE)
      {
        BigIntSubt(&biP, &intSqrt, &Temp);
        BigIntSubt(&Temp, &biM, &Temp);
        if (Temp.sign == SIGN_NEGATIVE || (Temp.nbrLimbs == 1 && Temp.limbs[0].x == 0))
        {
          BigIntSubt(&biM, &intSqrt, &Temp);
          if (Temp.sign == SIGN_NEGATIVE || (Temp.nbrLimbs == 1 && Temp.limbs[0].x == 0))
          {
            showText("<span class=\"bold\">");
            CopyBigInt(&P, &biP);
            CopyBigInt(&K, &P);
            CopyBigInt(&M, &biM);
            CopyBigInt(&L, &M);
            cont = 0;
          }
        }
      }
      if (cont >= 0)
      {
        /* both numerator and denominator are positive */
        (void)BigIntMultiply(&M, &M, &Temp);
        BigIntSubt(&delta, &Temp, &Temp);
        (void)BigIntDivide(&Temp, &P, &P);     // P <- (delta - M^2) / P
        BigIntAdd(&intSqrt, &M, &Z);
        (void)BigIntDivide(&Z, &P, &Z);        // Z <- (intSqrt + M) / P
        (void)BigIntMultiply(&Z, &P, &Temp);
        BigIntSubt(&Temp, &M, &M);       // M <- Z * P - M
        cont++;
      }
      else
      {
        (void)BigIntMultiply(&biM, &biM, &Temp);
        BigIntSubt(&delta, &Temp, &Temp);
        (void)BigIntDivide(&Temp, &biP, &biP);    // biP <- (delta - biM^2) / biP
        if (biP.sign == SIGN_POSITIVE)
        {
          BigIntAdd(&intSqrt, &biM, &Z);
          (void)BigIntDivide(&Z, &biP, &Z);      // Z <- (intSqrt + biM) / biP
        }
        else
        {
          CopyBigInt(&Z, &intSqrt);
          addbigint(&Z, 1);
          BigIntAdd(&Z, &biM, &Z);
          (void)BigIntDivide(&Z, &biP, &Z);      // Z <- (intSqrt + 1 + biM) / biP
        }
        (void)BigIntMultiply(&Z, &biP, &Temp);
        BigIntSubt(&Temp, &biM, &biM);     // biM <- Z * biP - biM
      }
      BigInteger2Dec(&Z, ptrOutput, groupLen);  // Show convergent.
      ptrOutput += strlen(ptrOutput);
    }
    if (cont >= 100000)
    {  // Too many convergents.
      showText(lang != 0? ", ... </span>//<br />donde la parte periódica (truncada a partir de los 100000 convergentes) está señalada en negrita.</p>" :
        ", ... </span>//<br />where the periodic part (truncated after 100000 convergents) is marked in bold.</p>");
      return;
    }
    else
    {
      showText(lang != 0 ? "</span>//<br />donde la parte periódica está señalada en negrita" :
        "</span>//<br />where the periodic part is marked in bold");
      if (cont > 1)
      {
        showText(lang != 0 ? " (el período tiene " : " (the period has ");
        int2dec(&ptrOutput, cont);
        showText(lang != 0 ? " coeficientes)" : " coefficients)");
      }
    }
  }
  showText(".");
  showText(lang? "<p>" COPYRIGHT_SPANISH "</p>" :
                 "<p>" COPYRIGHT_ENGLISH "</p>");
}

static void ShowRational(BigInteger *pNum, BigInteger *pDen)
{
  BigInteger GcdAll;
  char *sep;

  BigIntGcd(pNum, pDen, &GcdAll);
  (void)BigIntDivide(pNum, &GcdAll, pNum);
  (void)BigIntDivide(pDen, &GcdAll, pDen);
  if (pDen->sign == SIGN_NEGATIVE)
  {
    BigIntNegate(pNum, pNum);
    BigIntNegate(pDen, pDen);
  }
  (void)BigIntDivide(pNum, pDen, &Temp); //  Get integer part.
  (void)BigIntRemainder(pNum, pDen, pNum);
  if (pNum->sign == SIGN_NEGATIVE)
  {
    BigIntAdd(pDen, pNum, pNum);   // Convert numerator to positive.
    addbigint(&Temp, -1);          // Adjust integer part to floor.
  }
  BigInteger2Dec(&Temp, ptrOutput, groupLen);  // Show convergent.
  ptrOutput += strlen(ptrOutput);
  (void)BigIntRemainder(pNum, pDen, pNum);
  sep = " + //";
  while (pNum->nbrLimbs!=1 || pNum->limbs[0].x != 0)
  {      // Numerator greater than zero.
    if (pDen->nbrLimbs!=1 || pDen->limbs[0].x != 0)
    {    // Denominator greater than zero.
      (void)BigIntDivide(pDen, pNum, &Temp);
      showText(sep);
      BigInteger2Dec(&Temp, ptrOutput, groupLen);  // Show convergent.
      ptrOutput += strlen(ptrOutput);
      sep = ", ";
    }
    (void)BigIntRemainder(pDen, pNum, &Temp);
    CopyBigInt(pDen, pNum);
    CopyBigInt(pNum, &Temp);
  }
  if (sep[0] == ',')
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
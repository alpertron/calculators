/*
This file is part of Alpertron Calculators.

Copyright 2016 Dario Alejandro Alpern

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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bignbr.h"
#include "expression.h"
#include "factor.h"
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
extern long long lModularMult;
#endif
extern char *output;
char batch;
extern char verbose, prettyprint, cunningham;
extern BigInteger tofactor;
static int nbrToFactor[MAX_LEN];
static BigInteger Quad1, Quad2, Quad3, Quad4;
int groupLen = 6;
struct sFactors astFactorsMod[1000];
int factorsMod[10000];
extern BigInteger factorValue;
static BigInteger result;
BigInteger valueX;
char outputExpr[100000];
char *ptrInputText;
extern int NextEC;

#ifdef __EMSCRIPTEN__
extern double originalTenthSecond;
void GetDHMSt(char **pptrText, int tenths);
#endif

static void stringToHTML(char **pptrOutput, char *ptrString)
{
  int character;
  char c;
  char *ptrOutput = *pptrOutput;
  for (;;)
  {
    c = *ptrString;
    if (c == 0)
    {
      break;    // End of string, so go out.
    }
    if ((c & 0x80) == 0)
    {           // 1-byte UTF-8 character
      character = c;
      ptrString++;
    }
    else if ((c & 0x60) == 0x40)
    {            // 2-byte UTF-8 character
      character = (int)(c & 0x1F) * 64 + (*(ptrString + 1) & 0x3F);
      ptrString += 2;
    }
    else
    {            // 3-byte UTF-8 character
      character = ((int)(c & 0x1F) << 12) +
        ((int)(*(ptrString + 1) & 0x3F) << 6) +
        (*(ptrString + 2) & 0x3F);
      ptrString += 3;
    }
    if (character >= ' ' && character < 127 && character != '<' &&
      character != '>' && character != '&')
    {  // Safe to copy character.
      *ptrOutput++ = (char)character;
    }
    else
    {  // Insert HTML entity.
      *ptrOutput++ = '&';
      *ptrOutput++ = '#';
      int2dec(&ptrOutput, character);
      *ptrOutput++ = ';';
    }
  }
  *pptrOutput = ptrOutput;
}

static void SkipSpaces(char **pptrText)
{
  char *ptrText = *pptrText;
  while (*ptrText == ' ' || *ptrText == 9)
  {  // Skip spaces or tabs.
    ptrText++;
  }
  *pptrText = ptrText;
}

static char evalExpression(char *expr, int counter, BigInteger *result)
{
  char *ptrInputExpr = expr;
  char *ptrOutputExpr = outputExpr;
  char c;
  while (*ptrInputExpr != 0)
  {
    c = *ptrInputExpr;
    if (c == ';')
    {
      break;
    }
    if (c == 'x' || c == 'X')
    {
      BigInteger2Dec(&valueX, ptrOutputExpr, -100000);
      ptrOutputExpr += strlen(ptrOutputExpr);
    }
    else if (c == 'c' || c == 'C')
    {
      int2dec(&ptrOutputExpr, counter);
    }
    else
    {  // Copy character to output expression.
      *ptrOutputExpr++ = *ptrInputExpr;
    }
    ptrInputExpr++;
  }
  *ptrOutputExpr = 0;   // Append string terminator.
  return ComputeExpression(outputExpr, 1, result);
}

static void BatchError(char **pptrOutput, char *tofactorText, char *errorText)
{
  char *ptrOutput = *pptrOutput;
  stringToHTML(&ptrOutput, tofactorText);
  *ptrOutput++ = ':';
  *ptrOutput++ = ' ';
  strcpy(ptrOutput, errorText);
  ptrOutput += strlen(ptrOutput);
  *pptrOutput = ptrOutput;
}

static void BatchFactorization(char *tofactorText, int doFactorization)
{
  char *ptrNextBatchFactor = tofactorText;
  char *ptrEndBatchFactor = tofactorText + strlen(tofactorText);
  char *ptrCurrBatchFactor;
  int counter = 0;
  char *ptrOutput = output;
  char *NextExpr, *EndExpr, *FactorExpr;
  char c;
  enum eExprErr rc;
  char *ptrCharFound, *ptrSrcString, *ptrStartExpr;
  strcpy(ptrOutput, "2<ul><li>");
  ptrOutput += strlen(ptrOutput);
  while (ptrNextBatchFactor != NULL && ptrNextBatchFactor < ptrEndBatchFactor)
  {  // Get next line.
    counter = 1;
    ptrCurrBatchFactor = ptrNextBatchFactor;
    ptrNextBatchFactor = findChar(ptrCurrBatchFactor, '\n');
    if (ptrNextBatchFactor != NULL)
    {    // End of line found.
      *ptrNextBatchFactor++ = 0;   // Replace it by string terminator.
    }
    // Skip leading spaces.
    ptrSrcString = ptrCurrBatchFactor;
    SkipSpaces(&ptrSrcString);
    c = *ptrSrcString;
    if (c == 0)
    {   // Empty line.
      strcpy(ptrOutput, "<br>");
      ptrOutput += strlen(ptrOutput);
    }
    else if (c == '#')
    {   // Copy comment to output, but convert non-safe characters to entities.
      *ptrOutput++ = '#';
      stringToHTML(&ptrOutput, ptrCurrBatchFactor + 1);
    }
    else if (c == 'x' || c == 'X')
    {   // Loop format: x=<orig expr>; x=<next expr>; <end expr>; <expr to factor>
      ptrCharFound = findChar(ptrSrcString + 1, ';');
      if (ptrCharFound == NULL)
      {
        BatchError(&ptrOutput, tofactorText,
          lang ? "se esperaban tres puntos y comas pero no hay ninguno" :
          "three semicolons expected but none found");
        continue;
      }
      ptrStartExpr = ptrSrcString + 1;
      SkipSpaces(&ptrStartExpr);
      if (*ptrStartExpr != '=')
      {
        BatchError(&ptrOutput, tofactorText,
          lang ? "falta signo igual en la primera expresión" :
          "equal sign missing in first expression");
        continue;
      }
      ptrCharFound = findChar(ptrSrcString + 1, ';');
      *ptrCharFound++ = 0;   // Replace first semicolon by end of string.
      if (ComputeExpression(ptrStartExpr + 1, 1, &valueX) != 0)
      {
        continue;
      }
      ptrStartExpr = ptrCharFound;
      SkipSpaces(&ptrStartExpr);
      if (*ptrStartExpr != 'x' && *ptrStartExpr != 'X')
      {
        BatchError(&ptrOutput, tofactorText,
          lang ? "falta variable x en la segunda expresión" :
          "variable x missing in second expression");
        continue;
      }
      ptrStartExpr++;               // Skip variable 'x'.
      SkipSpaces(&ptrStartExpr);
      if (*ptrStartExpr != '=')
      {
        BatchError(&ptrOutput, tofactorText,
          lang ? "falta signo igual en la segunda expresión" :
          "equal sign missing in second expression");
        continue;
      }
      NextExpr = ptrStartExpr + 1;  // Skip equal sign.
      ptrCharFound = findChar(ptrStartExpr, ';');  // Find second semicolon.
      if (ptrCharFound == NULL)
      {      // Third semicolon not found.
        BatchError(&ptrOutput, tofactorText,
          lang ? "se esperaban tres puntos y comas pero solo hay uno" :
          "three semicolons expected but there are only one");
        continue;
      }
      EndExpr = ptrCharFound + 1;  // Point to end expression.
      ptrCharFound = findChar(EndExpr, ';');  // Find third semicolon.
      if (ptrCharFound == NULL)
      {      // Third semicolon not found.
        BatchError(&ptrOutput, tofactorText,
          lang ? "se esperaban tres puntos y comas pero solo hay dos" :
          "three semicolons expected but there are only two");
        continue;
      }
      FactorExpr = ptrCharFound + 1;
      for (;;)
      {   // Perform loop.
        if (evalExpression(EndExpr, counter, &result) != 0)
        {
          break;
        }
        if (result.sign == SIGN_POSITIVE && (result.nbrLimbs > 1 || result.limbs[0].x > 1))
        {   // result is greater than zero: end of loop
          break;
        }
        rc = evalExpression(FactorExpr, counter, &tofactor);
        if (rc == EXPR_OK && doFactorization)
        {
          NumberLength = tofactor.nbrLimbs;
          CompressBigInteger(nbrToFactor, &tofactor);
          factor(&tofactor, nbrToFactor, factorsMod, astFactorsMod, NULL);
        }
        SendFactorizationToOutput(rc, astFactorsMod, &ptrOutput, doFactorization);
        if (evalExpression(NextExpr, counter, &result) != 0)
        {
          break;
        }
        CopyBigInt(&valueX, &result);
        counter++;
        strcpy(ptrOutput, "</li><li>");
        ptrOutput += strlen(ptrOutput);
      }
    }
    else
    {       // Factor expression (not loop).
      rc = ComputeExpression(ptrSrcString, 1, &tofactor);
      if (rc == EXPR_OK && doFactorization)
      {
        NumberLength = tofactor.nbrLimbs;
        CompressBigInteger(nbrToFactor, &tofactor);
        factor(&tofactor, nbrToFactor, factorsMod, astFactorsMod, NULL);
      }
      SendFactorizationToOutput(rc, astFactorsMod, &ptrOutput, doFactorization);
      counter = 2;
      strcpy(ptrOutput, "</li><li>");
      ptrOutput += strlen(ptrOutput);
    }
    if (counter == 1)
    {
      strcpy(ptrOutput, "</li>");
      ptrOutput += strlen(ptrOutput);
    }
    ptrCurrBatchFactor = ptrNextBatchFactor;
  }
  if (counter > 1)
  {
    ptrOutput -= 4;   // Erase start tag <li> without contents.
  }
  strcpy(ptrOutput, "</ul>");
  ptrOutput += strlen(ptrOutput);
  strcpy(ptrOutput, lang ? "<p>" COPYRIGHT_SPANISH "</p>" :
    "<p>" COPYRIGHT_ENGLISH "</p>");
}

static void ExponentToBigInteger(int exponent, BigInteger *bigint)
{
  if (exponent > MAX_VALUE_LIMB)
  {
    bigint->limbs[0].x = exponent - MAX_VALUE_LIMB;
    bigint->limbs[1].x = 1;
    bigint->nbrLimbs = 2;
  }
  else
  {
    bigint->limbs[0].x = exponent + 1;
    bigint->nbrLimbs = 1;
  }
  bigint->sign = SIGN_POSITIVE;
}

// Find number of divisors as the product of all exponents plus 1.
static void GetNumberOfDivisors(char **pptrOutput)
{
  char *ptrOutput = *pptrOutput;
  int *ptrFactor;
  struct sFactors *pstFactor;
  int factorNumber;
  result.limbs[0].x = 1;    // Set result to 1.
  result.nbrLimbs = 1;
  result.sign = SIGN_POSITIVE;
  pstFactor = &astFactorsMod[1];
  for (factorNumber = 1; factorNumber <= astFactorsMod[0].multiplicity; factorNumber++)
  {
    ptrFactor = pstFactor->ptrFactor;
    if (*ptrFactor == 1 && *(ptrFactor + 1) < 2)
    {                        // Factor is 1.
      break;
    }
    ExponentToBigInteger(pstFactor->multiplicity, &factorValue);
    BigIntMultiply(&factorValue, &result, &result);
    pstFactor++;
  }
  strcpy(ptrOutput, lang ? "<p>Cantidad de divisores: " : "<p>Number of divisors: ");
  ptrOutput += strlen(ptrOutput);
  BigInteger2Dec(&result, ptrOutput, groupLen);
  ptrOutput += strlen(ptrOutput);
  strcpy(ptrOutput, "</p>");
  ptrOutput += strlen(ptrOutput);
  *pptrOutput = ptrOutput;
}

// Find sum of divisors as the product of (p^(e+1)-1)/(p-1) where p=prime and e=exponent.
static void GetSumOfDivisors(char **pptrOutput)
{
  BigInteger Temp1;
  char *ptrOutput = *pptrOutput;
  struct sFactors *pstFactor;
  int factorNumber;
  result.limbs[0].x = 1;    // Set result to 1.
  result.nbrLimbs = 1;
  result.sign = SIGN_POSITIVE;
  pstFactor = &astFactorsMod[1];
  for (factorNumber = 1; factorNumber <= astFactorsMod[0].multiplicity; factorNumber++)
  {
    UncompressBigInteger(pstFactor->ptrFactor, &factorValue);
    BigIntPowerIntExp(&factorValue, pstFactor->multiplicity + 1, &Temp1);   // p^(e+1)
    addbigint(&Temp1, -1);   // p^(e+1)-1
    BigIntMultiply(&result, &Temp1, &result);
    UncompressBigInteger(pstFactor->ptrFactor, &Temp1);
    addbigint(&Temp1, -1);   // p-1
    BigIntDivide(&result, &Temp1, &result);
    pstFactor++;
  }
  strcpy(ptrOutput, lang ? "<p>Suma de divisores: " : "<p>Sum of divisors: ");
  ptrOutput += strlen(ptrOutput);
  BigInteger2Dec(&result, ptrOutput, groupLen);
  ptrOutput += strlen(ptrOutput);
  strcpy(ptrOutput, "</p>");
  ptrOutput += strlen(ptrOutput);
  *pptrOutput = ptrOutput;
}

// Find Euler's Totient as the product of p^(e-1)*(p-1) where p=prime and e=exponent.
static void GetEulerTotient(char **pptrOutput)
{
  BigInteger Temp1;
  char *ptrOutput = *pptrOutput;
  struct sFactors *pstFactor;
  int factorNumber;
  result.limbs[0].x = 1;    // Set result to 1.
  result.nbrLimbs = 1;
  result.sign = SIGN_POSITIVE;
  pstFactor = &astFactorsMod[1];
  for (factorNumber = 1; factorNumber <= astFactorsMod[0].multiplicity; factorNumber++)
  {
    UncompressBigInteger(pstFactor->ptrFactor, &factorValue);
    BigIntPowerIntExp(&factorValue, pstFactor->multiplicity - 1, &Temp1);   // p^(e-1)
    BigIntMultiply(&result, &Temp1, &result);
    UncompressBigInteger(pstFactor->ptrFactor, &Temp1);
    addbigint(&Temp1, -1);   // p-1
    BigIntMultiply(&result, &Temp1, &result);
    pstFactor++;
  }
  strcpy(ptrOutput, lang ? "<p>Phi de Euler: " : "<p>Euler's totient: ");
  ptrOutput += strlen(ptrOutput);
  BigInteger2Dec(&result, ptrOutput, groupLen);
  ptrOutput += strlen(ptrOutput);
  strcpy(ptrOutput, "</p>");
  ptrOutput += strlen(ptrOutput);
  *pptrOutput = ptrOutput;
}

// Find Mobius as zero if some exponent is > 1, 1 if the number of factors is even, -1 if it is odd.
static void GetMobius(char **pptrOutput)
{
  char *ptrOutput = *pptrOutput;
  struct sFactors *pstFactor;
  int factorNumber;
  int mobius = 1;
  pstFactor = &astFactorsMod[1];
  if (astFactorsMod[0].multiplicity > 1 || *pstFactor->ptrFactor != 1 ||
    *(pstFactor->ptrFactor + 1) != 1)
  {                                // Number to factor is not 1.
    for (factorNumber = 1; factorNumber <= astFactorsMod[0].multiplicity; factorNumber++)
    {
      if (pstFactor->multiplicity == 1)
      {
        mobius = -mobius;
      }
      else
      {
        mobius = 0;
      }
      pstFactor++;
    }
  }
  strcpy(ptrOutput, "<p>Möbius: ");
  ptrOutput += strlen(ptrOutput);
  if (mobius < 0)
  {
    mobius = -mobius;
    *ptrOutput++ = '-';
  }
  int2dec(&ptrOutput, mobius);
  strcpy(ptrOutput, "</p>");
  ptrOutput += strlen(ptrOutput);
  *pptrOutput = ptrOutput;
}

static void ComputeFourSquares(struct sFactors *pstFactors)
{
  int indexPrimes;
  BigInteger p, q, K, Mult1, Mult2, Mult3, Mult4;
  BigInteger Tmp, Tmp1, Tmp2, Tmp3, Tmp4, M1, M2, M3, M4;
  struct sFactors *pstFactor;
  static limb minusOneMont[MAX_LEN];

  intToBigInteger(&Quad1, 1);     // 1 = 1^2 + 0^2 + 0^2 + 0^2
  intToBigInteger(&Quad2, 0);
  intToBigInteger(&Quad3, 0);
  intToBigInteger(&Quad4, 0);
  pstFactor = pstFactors + 1;      // Point to first factor in array of factors.
  if (pstFactors->multiplicity == 1 && *pstFactor->ptrFactor == 1)
  {
    if (*(pstFactor->ptrFactor + 1) == 1)
    {                              // Number to factor is 1.
      return;
    }
    if (*(pstFactor->ptrFactor + 1) == 0)
    {                             // Number to factor is 0.
      intToBigInteger(&Quad1, 0);     // 0 = 0^2 + 0^2 + 0^2 + 0^2
      return;
    }
  }
  for (indexPrimes = pstFactors -> multiplicity - 1; indexPrimes >= 0; 
       indexPrimes--, pstFactor++)
  {
    if (pstFactor -> multiplicity % 2 == 0)
    {                              // Prime factor appears twice.
      continue;
    }
    NumberLength = *pstFactor->ptrFactor;
    UncompressBigInteger(pstFactor->ptrFactor, &p);
    p.sign = SIGN_POSITIVE;
    CopyBigInt(&q, &p);
    addbigint(&q, -1);             // q <- p-1
    if (p.nbrLimbs == 1 && p.limbs[0].x == 2)
    {
      intToBigInteger(&Mult1, 1); // 2 = 1^2 + 1^2 + 0^2 + 0^2
      intToBigInteger(&Mult2, 1);
      intToBigInteger(&Mult3, 0);
      intToBigInteger(&Mult4, 0);
    }
    else
    { /* Prime not 2 */
      NumberLength = p.nbrLimbs;
      memcpy(&TestNbr, p.limbs, NumberLength * sizeof(limb));
      TestNbr[NumberLength].x = 0;
      GetMontgomeryParms(NumberLength);
      memset(minusOneMont, 0, NumberLength * sizeof(limb));
      SubtBigNbrModN(minusOneMont, MontgomeryMultR1, minusOneMont, TestNbr, NumberLength);
      memset(K.limbs, 0, NumberLength * sizeof(limb));
      if ((p.limbs[0].x & 3) == 1)
      { /* if p = 1 (mod 4) */
        CopyBigInt(&q, &p);
        subtractdivide(&q, 1, 4);     // q = (prime-1)/4
        K.limbs[0].x = 1;
        do
        {    // Loop that finds mult1 = sqrt(-1) mod prime in Montgomery notation.
          K.limbs[0].x++;
          modPow(K.limbs, q.limbs, q.nbrLimbs, Mult1.limbs);
        } while (!memcmp(Mult1.limbs, MontgomeryMultR1, NumberLength * sizeof(limb)) ||
          !memcmp(Mult1.limbs, minusOneMont, NumberLength * sizeof(limb)));
        Mult1.sign = SIGN_POSITIVE;
        memset(Mult2.limbs, 0, p.nbrLimbs * sizeof(limb));
        Mult2.limbs[0].x = 1;
        Mult2.nbrLimbs = 1;
        Mult2.sign = SIGN_POSITIVE;
        // Convert Mult1 to standard notation by multiplying by 1 in
        // Montgomery notation.
        modmult(Mult1.limbs, Mult2.limbs, Mult3.limbs);
        memcpy(Mult1.limbs, Mult3.limbs, p.nbrLimbs * sizeof(limb));
        for (Mult1.nbrLimbs = p.nbrLimbs; Mult1.nbrLimbs > 1; Mult1.nbrLimbs--)
        {  // Adjust number of limbs so the most significant limb is not zero.
          if (Mult1.limbs[Mult1.nbrLimbs - 1].x != 0)
          {
            break;
          }
        }
        for (;;)
        {
          BigIntMultiply(&Mult1, &Mult1, &Tmp);
          BigIntMultiply(&Mult2, &Mult2, &Tmp1);
          BigIntAdd(&Tmp, &Tmp1, &Tmp);
          BigIntDivide(&Tmp, &p, &K);        // K <- (mult1^2 + mult2^2) / p
          if (K.nbrLimbs == 1 && K.limbs[0].x == 1)
          {    // If K = 1...
            intToBigInteger(&Mult3, 0);
            intToBigInteger(&Mult4, 0);
            break;
          }
          BigIntRemainder(&Mult1, &K, &M1);  // M1 <- Mult1 % K
          if (M1.sign == SIGN_NEGATIVE)
          {
            BigIntAdd(&M1, &K, &M1);
          }
          BigIntRemainder(&Mult2, &K, &M2);  // M2 <- Mult2 % K
          if (M2.sign == SIGN_NEGATIVE)
          {
            BigIntAdd(&M2, &K, &M2);
          }
          CopyBigInt(&Tmp, &K);
          subtractdivide(&Tmp, -1, 2);       // Tmp <- (K+1) / 2
          BigIntSubt(&M1, &Tmp, &Tmp1);      // Tmp1 <- M1 - Tmp
          if (Tmp1.sign == SIGN_POSITIVE)    // If M1 >= K / 2 ... 
          {
            BigIntSubt(&M1, &K, &M1);        // M1 <- M1 - K
          }
          BigIntSubt(&M2, &Tmp, &Tmp1);      // Tmp1 <- M2 - Tmp
          if (Tmp1.sign == SIGN_POSITIVE)    // If M2 >= K / 2 ... 
          {
            BigIntSubt(&M2, &K, &M2);        // M1 <- M1 - K
          }
          BigIntMultiply(&Mult1, &M1, &Tmp);
          BigIntMultiply(&Mult2, &M2, &Tmp1);
          BigIntAdd(&Tmp, &Tmp1, &Tmp);
          BigIntDivide(&Tmp, &K, &Tmp2);     // Tmp2 <- (mult1*m1 + mult2*m2) / K
          BigIntMultiply(&Mult1, &M2, &Tmp);
          BigIntMultiply(&Mult2, &M1, &Tmp1);
          BigIntSubt(&Tmp, &Tmp1, &Tmp);
          BigIntDivide(&Tmp, &K, &Mult2);    // Mult2 <- (mult1*m2 - mult2*m1) / K
          CopyBigInt(&Mult1, &Tmp2);
        } /* end while */
      } /* end p = 1 (mod 4) */
      else
      { /* if p = 3 (mod 4) */
        int mult1 = 0;
        CopyBigInt(&q, &p);
        subtractdivide(&q, 1, 2);     // q = (prime-1)/2
        memcpy(K.limbs, q.limbs, q.nbrLimbs * sizeof(limb));
        // Compute Mult1 and Mult2 so Mult1^2 + Mult2^2 = -1 (mod p)
        memset(Mult1.limbs, 0, p.nbrLimbs*sizeof(limb));
        do
        {
          mult1++;
          // Increment Mult1 by 1 in Montgomery notation.
          AddBigNbrModN(Mult1.limbs, MontgomeryMultR1, Mult1.limbs, p.limbs, p.nbrLimbs);
          modmult(Mult1.limbs, Mult1.limbs, Tmp.limbs);
          SubtBigNbrModN(minusOneMont, Tmp.limbs, Tmp.limbs, p.limbs, p.nbrLimbs);
          modPow(Tmp.limbs, K.limbs, p.nbrLimbs, Tmp1.limbs);
          // At this moment Tmp1 = (-1 - Mult1^2)^((p-1)/2)
          // in Montgomery notation. Continue loop if it is not 1.
        } while (memcmp(Tmp1.limbs, MontgomeryMultR1, p.nbrLimbs));
        // After the loop finishes, Tmp = (-1 - Mult1^2) is a quadratic residue mod p.
        // Convert Mult1 to standard notation by multiplying by 1 in
        // Montgomery notation.
        intToBigInteger(&Mult1, mult1);
        CopyBigInt(&q, &p);
        subtractdivide(&q, -1, 4);  // q <- (p+1)/4.
        // Find Mult2 <- square root of Tmp = Tmp^q (mod p) in Montgomery notation.
        modPow(Tmp.limbs, q.limbs, p.nbrLimbs, Mult2.limbs);
        // Convert Mult2 from Montgomery notation to standard notation.
        memset(Tmp.limbs, 0, p.nbrLimbs * sizeof(limb));
        Tmp.limbs[0].x = 1;
        intToBigInteger(&Mult3, 1);
        intToBigInteger(&Mult4, 0);
        // Convert Mult2 to standard notation by multiplying by 1 in
        // Montgomery notation.
        modmult(Mult2.limbs, Tmp.limbs, Mult2.limbs);
        for (Mult2.nbrLimbs = p.nbrLimbs; Mult2.nbrLimbs > 1; Mult2.nbrLimbs--)
        {  // Adjust number of limbs so the most significant limb is not zero.
          if (Mult2.limbs[Mult2.nbrLimbs - 1].x != 0)
          {
            break;
          }
        }
        Mult2.sign = SIGN_POSITIVE;
        for (;;)
        {
          // Compute K <- (Mult1^2 + Mult2^2 + Mult3^2 + Mult4^2) / p
          BigIntMultiply(&Mult1, &Mult1, &Tmp);
          BigIntMultiply(&Mult2, &Mult2, &Tmp1);
          BigIntAdd(&Tmp, &Tmp1, &Tmp);
          BigIntMultiply(&Mult3, &Mult3, &Tmp1);
          BigIntAdd(&Tmp, &Tmp1, &Tmp);
          BigIntMultiply(&Mult4, &Mult4, &Tmp1);
          BigIntAdd(&Tmp, &Tmp1, &Tmp);
          BigIntDivide(&Tmp, &p, &K);
          if (K.nbrLimbs == 1 && K.limbs[0].x == 1)
          {   // K equals 1
            break;
          }
          if ((K.limbs[0].x & 1) == 0)
          { // If K is even ...
            if ((Mult1.limbs[0].x + Mult2.limbs[0].x) & 1)
            {  // If Mult1 + Mult2 is odd...
              if (((Mult1.limbs[0].x + Mult3.limbs[0].x) & 1) == 0)
              {   // If Mult1 + Mult3 is even...
                CopyBigInt(&Tmp, &Mult2);
                CopyBigInt(&Mult2, &Mult3);
                CopyBigInt(&Mult3, &Tmp);
              }
              else
              {
                CopyBigInt(&Tmp, &Mult2);
                CopyBigInt(&Mult2, &Mult4);
                CopyBigInt(&Mult4, &Tmp);
              }
            } // At this moment Mult1+Mult2 = even, Mult3+Mult4 = even
            BigIntAdd(&Mult1, &Mult2, &Tmp1);
            subtractdivide(&Tmp1, 0, 2);  // Tmp1 <- (Mult1 + Mult2) / 2
            BigIntSubt(&Mult1, &Mult2, &Tmp2);
            subtractdivide(&Tmp2, 0, 2);  // Tmp2 <- (Mult1 - Mult2) / 2
            BigIntAdd(&Mult3, &Mult4, &Tmp3);
            subtractdivide(&Tmp3, 0, 2);  // Tmp3 <- (Mult3 + Mult4) / 2
            BigIntSubt(&Mult3, &Mult4, &Mult4);
            subtractdivide(&Mult4, 0, 2);  // Mult4 <- (Mult3 + Mult4) / 2
            CopyBigInt(&Mult3, &Tmp3);
            CopyBigInt(&Mult2, &Tmp2);
            CopyBigInt(&Mult1, &Tmp1);
            continue;
          } /* end if k is even */
          BigIntRemainder(&Mult1, &K, &M1);    // M1 <- Mult1 % K.
          if (M1.sign == SIGN_NEGATIVE)
          {
            BigIntAdd(&M1, &K, &M1);
          }
          BigIntRemainder(&Mult2, &K, &M2);    // M2 <- Mult2 % K.
          if (M2.sign == SIGN_NEGATIVE)
          {
            BigIntAdd(&M2, &K, &M2);
          }
          BigIntRemainder(&Mult3, &K, &M3);    // M3 <- Mult3 % K.
          if (M3.sign == SIGN_NEGATIVE)
          {
            BigIntAdd(&M3, &K, &M3);
          }
          BigIntRemainder(&Mult4, &K, &M4);    // M4 <- Mult4 % K.
          if (M4.sign == SIGN_NEGATIVE)
          {
            BigIntAdd(&M4, &K, &M4);
          }
          CopyBigInt(&Tmp, &K);
          subtractdivide(&Tmp, -1, 2);       // Tmp <- (K+1) / 2
          BigIntSubt(&M1, &Tmp, &Tmp1);      // Tmp1 <- M1 - Tmp
          if (Tmp1.sign == SIGN_POSITIVE)    // If M1 >= K / 2 ... 
          {
            BigIntSubt(&M1, &K, &M1);        // M1 <- M1 - K
          }
          BigIntSubt(&M2, &Tmp, &Tmp1);      // Tmp1 <- M2 - Tmp
          if (Tmp1.sign == SIGN_POSITIVE)    // If M2 >= K / 2 ... 
          {
            BigIntSubt(&M2, &K, &M2);        // M2 <- M2 - K
          }
          BigIntSubt(&M3, &Tmp, &Tmp1);      // Tmp1 <- M3 - Tmp
          if (Tmp1.sign == SIGN_POSITIVE)    // If M3 >= K / 2 ... 
          {
            BigIntSubt(&M3, &K, &M3);        // M3 <- M3 - K
          }
          BigIntSubt(&M4, &Tmp, &Tmp1);      // Tmp1 <- M4 - Tmp
          if (Tmp1.sign == SIGN_POSITIVE)    // If M4 >= K / 2 ... 
          {
            BigIntSubt(&M4, &K, &M4);        // M4 <- M4 - K
          }
          // Compute Tmp1 <- (Mult1*M1 + Mult2*M2 + Mult3*M3 + Mult4*M4) / K
          BigIntMultiply(&Mult1, &M1, &Tmp);
          BigIntMultiply(&Mult2, &M2, &Tmp4);
          BigIntAdd(&Tmp, &Tmp4, &Tmp);
          BigIntMultiply(&Mult3, &M3, &Tmp4);
          BigIntAdd(&Tmp, &Tmp4, &Tmp);
          BigIntMultiply(&Mult4, &M4, &Tmp4);
          BigIntAdd(&Tmp, &Tmp4, &Tmp);
          BigIntDivide(&Tmp, &K, &Tmp1);

          // Compute Tmp2 <- (Mult1*M2 - Mult2*M1 + Mult3*M4 - Mult4*M3) / K
          BigIntMultiply(&Mult1, &M2, &Tmp);
          BigIntMultiply(&Mult2, &M1, &Tmp4);
          BigIntSubt(&Tmp, &Tmp4, &Tmp);
          BigIntMultiply(&Mult3, &M4, &Tmp4);
          BigIntAdd(&Tmp, &Tmp4, &Tmp);
          BigIntMultiply(&Mult4, &M3, &Tmp4);
          BigIntSubt(&Tmp, &Tmp4, &Tmp);
          BigIntDivide(&Tmp, &K, &Tmp2);

          // Compute Tmp3 <- (Mult1*M3 - Mult3*M1 - Mult2*M4 + Mult4*M2) / K
          BigIntMultiply(&Mult1, &M3, &Tmp);
          BigIntMultiply(&Mult3, &M1, &Tmp4);
          BigIntSubt(&Tmp, &Tmp4, &Tmp);
          BigIntMultiply(&Mult2, &M4, &Tmp4);
          BigIntSubt(&Tmp, &Tmp4, &Tmp);
          BigIntMultiply(&Mult4, &M2, &Tmp4);
          BigIntAdd(&Tmp, &Tmp4, &Tmp);
          BigIntDivide(&Tmp, &K, &Tmp3);

          // Compute Mult4 <- (Mult1*M4 - Mult4*M1 + Mult2*M3 - Mult3*M2) / K
          BigIntMultiply(&Mult1, &M4, &Tmp);
          BigIntMultiply(&Mult4, &M1, &Tmp4);
          BigIntSubt(&Tmp, &Tmp4, &Tmp);
          BigIntMultiply(&Mult2, &M3, &Tmp4);
          BigIntAdd(&Tmp, &Tmp4, &Tmp);
          BigIntMultiply(&Mult3, &M2, &Tmp4);
          BigIntSubt(&Tmp, &Tmp4, &Tmp);
          BigIntDivide(&Tmp, &K, &Mult4);

          CopyBigInt(&Mult3, &Tmp3);
          CopyBigInt(&Mult2, &Tmp2);
          CopyBigInt(&Mult1, &Tmp1);
        } /* end while */
      } /* end if p = 3 (mod 4) */
    } /* end prime not 2 */

    // Compute Tmp1 <- Mult1*Quad1 + Mult2*Quad2 + Mult3*Quad3 + Mult4*Quad4
    BigIntMultiply(&Mult1, &Quad1, &Tmp);
    BigIntMultiply(&Mult2, &Quad2, &Tmp4);
    BigIntAdd(&Tmp, &Tmp4, &Tmp);
    BigIntMultiply(&Mult3, &Quad3, &Tmp4);
    BigIntAdd(&Tmp, &Tmp4, &Tmp);
    BigIntMultiply(&Mult4, &Quad4, &Tmp4);
    BigIntAdd(&Tmp, &Tmp4, &Tmp1);

    // Compute Tmp2 <- Mult1*Quad2 - Mult2*Quad1 + Mult3*Quad4 - Mult4*Quad3
    BigIntMultiply(&Mult1, &Quad2, &Tmp);
    BigIntMultiply(&Mult2, &Quad1, &Tmp4);
    BigIntSubt(&Tmp, &Tmp4, &Tmp);
    BigIntMultiply(&Mult3, &Quad4, &Tmp4);
    BigIntAdd(&Tmp, &Tmp4, &Tmp);
    BigIntMultiply(&Mult4, &Quad3, &Tmp4);
    BigIntSubt(&Tmp, &Tmp4, &Tmp2);

    // Compute Tmp3 <- Mult1*Quad3 - Mult3*Quad1 - Mult2*Quad4 + Mult4*Quad2
    BigIntMultiply(&Mult1, &Quad3, &Tmp);
    BigIntMultiply(&Mult3, &Quad1, &Tmp4);
    BigIntSubt(&Tmp, &Tmp4, &Tmp);
    BigIntMultiply(&Mult2, &Quad4, &Tmp4);
    BigIntSubt(&Tmp, &Tmp4, &Tmp);
    BigIntMultiply(&Mult4, &Quad2, &Tmp4);
    BigIntAdd(&Tmp, &Tmp4, &Tmp3);

    // Compute Quad4 <- Mult1*Quad4 - Mult4*Quad1 + Mult2*Quad3 - Mult3*Quad2
    BigIntMultiply(&Mult1, &Quad4, &Tmp);
    BigIntMultiply(&Mult4, &Quad1, &Tmp4);
    BigIntSubt(&Tmp, &Tmp4, &Tmp);
    BigIntMultiply(&Mult2, &Quad3, &Tmp4);
    BigIntAdd(&Tmp, &Tmp4, &Tmp);
    BigIntMultiply(&Mult3, &Quad2, &Tmp4);
    BigIntSubt(&Tmp, &Tmp4, &Quad4);

    CopyBigInt(&Quad3, &Tmp3);
    CopyBigInt(&Quad2, &Tmp2);
    CopyBigInt(&Quad1, &Tmp1);
  } /* end for indexPrimes */
  pstFactor = pstFactors + 1;      // Point to first factor in array of factors.
  for (indexPrimes = pstFactors->multiplicity - 1; indexPrimes >= 0;
    indexPrimes--, pstFactor++)
  {
    NumberLength = *pstFactor->ptrFactor;
    UncompressBigInteger(pstFactor->ptrFactor, &p);
    BigIntPowerIntExp(&p, pstFactor->multiplicity / 2, &K);
    BigIntMultiply(&Quad1, &K, &Quad1);
    BigIntMultiply(&Quad2, &K, &Quad2);
    BigIntMultiply(&Quad3, &K, &Quad3);
    BigIntMultiply(&Quad4, &K, &Quad4);
  }
  Quad1.sign = SIGN_POSITIVE;
  Quad2.sign = SIGN_POSITIVE;
  Quad3.sign = SIGN_POSITIVE;
  Quad4.sign = SIGN_POSITIVE;
  // Sort squares
  BigIntSubt(&Quad1, &Quad2, &Tmp);
  if (Tmp.sign == SIGN_NEGATIVE)
  {   // Quad1 < Quad2, so exchange them.
    CopyBigInt(&Tmp, &Quad1);
    CopyBigInt(&Quad1, &Quad2);
    CopyBigInt(&Quad2, &Tmp);
  }
  BigIntSubt(&Quad1, &Quad3, &Tmp);
  if (Tmp.sign == SIGN_NEGATIVE)
  {   // Quad1 < Quad3, so exchange them.
    CopyBigInt(&Tmp, &Quad1);
    CopyBigInt(&Quad1, &Quad3);
    CopyBigInt(&Quad3, &Tmp);
  }
  BigIntSubt(&Quad1, &Quad4, &Tmp);
  if (Tmp.sign == SIGN_NEGATIVE)
  {   // Quad1 < Quad4, so exchange them.
    CopyBigInt(&Tmp, &Quad1);
    CopyBigInt(&Quad1, &Quad4);
    CopyBigInt(&Quad4, &Tmp);
  }
  BigIntSubt(&Quad2, &Quad3, &Tmp);
  if (Tmp.sign == SIGN_NEGATIVE)
  {   // Quad2 < Quad3, so exchange them.
    CopyBigInt(&Tmp, &Quad2);
    CopyBigInt(&Quad2, &Quad3);
    CopyBigInt(&Quad3, &Tmp);
  }
  BigIntSubt(&Quad2, &Quad4, &Tmp);
  if (Tmp.sign == SIGN_NEGATIVE)
  {   // Quad2 < Quad4, so exchange them.
    CopyBigInt(&Tmp, &Quad2);
    CopyBigInt(&Quad2, &Quad4);
    CopyBigInt(&Quad4, &Tmp);
  }
  BigIntSubt(&Quad3, &Quad4, &Tmp);
  if (Tmp.sign == SIGN_NEGATIVE)
  {   // Quad3 < Quad4, so exchange them.
    CopyBigInt(&Tmp, &Quad3);
    CopyBigInt(&Quad3, &Quad4);
    CopyBigInt(&Quad4, &Tmp);
  }
}

static void varSquared(char **pptrOutput, char letter, char sign)
{
  char *ptrOutput = *pptrOutput;
  *ptrOutput++ = ' ';
  *ptrOutput++ = letter;
  strcpy(ptrOutput, (prettyprint? "&sup2;": "^2"));
  ptrOutput += strlen(ptrOutput);
  *ptrOutput++ = ' ';
  *ptrOutput++ = sign;
  *pptrOutput = ptrOutput;
}

static void valueVar(char **pptrOutput, char letter, BigInteger *value)
{
  char *ptrOutput = *pptrOutput;
  strcpy(ptrOutput, "<p>");
  ptrOutput += strlen(ptrOutput);
  *ptrOutput++ = letter;
  *ptrOutput++ = ' ';
  *ptrOutput++ = '=';
  *ptrOutput++ = ' ';
  BigInteger2Dec(value, ptrOutput, groupLen);
  ptrOutput += strlen(ptrOutput);
  strcpy(ptrOutput, "</p>");
  ptrOutput += strlen(ptrOutput);
  *pptrOutput = ptrOutput;
}

static void ShowFourSquares(char **pptrOutput)
{
  char *ptrOutput = *pptrOutput;
  strcpy(ptrOutput, "<p>n =");
  ptrOutput += strlen(ptrOutput);
  if (Quad4.nbrLimbs == 1 && Quad4.limbs[0].x == 0)
  {          // Quad4 equals zero.
    if (Quad3.nbrLimbs == 1 && Quad3.limbs[0].x == 0)
    {        // Quad3 and Quad4 equal zero.
      if (Quad2.nbrLimbs == 1 && Quad2.limbs[0].x == 0)
      {      // Quad2, Quad3 and Quad4 equal zero.
        varSquared(&ptrOutput, 'a', ' ');
        strcpy(ptrOutput, "</p>");
        ptrOutput += strlen(ptrOutput);
        valueVar(&ptrOutput, 'a', &Quad1);
        *pptrOutput = ptrOutput;
        return;
      }
      varSquared(&ptrOutput, 'a', '+');
      varSquared(&ptrOutput, 'b', ' ');
      strcpy(ptrOutput, "</p>");
      ptrOutput += strlen(ptrOutput);
      valueVar(&ptrOutput, 'a', &Quad1);
      valueVar(&ptrOutput, 'b', &Quad2);
      *pptrOutput = ptrOutput;
      return;
    }
    varSquared(&ptrOutput, 'a', '+');
    varSquared(&ptrOutput, 'b', '+');
    varSquared(&ptrOutput, 'c', ' ');
    strcpy(ptrOutput, "</p>");
    ptrOutput += strlen(ptrOutput);
    valueVar(&ptrOutput, 'a', &Quad1);
    valueVar(&ptrOutput, 'b', &Quad2);
    valueVar(&ptrOutput, 'c', &Quad3);
    *pptrOutput = ptrOutput;
    return;
  }
  varSquared(&ptrOutput, 'a', '+');
  varSquared(&ptrOutput, 'b', '+');
  varSquared(&ptrOutput, 'c', '+');
  varSquared(&ptrOutput, 'd', ' ');
  strcpy(ptrOutput, "</p>");
  ptrOutput += strlen(ptrOutput);
  valueVar(&ptrOutput, 'a', &Quad1);
  valueVar(&ptrOutput, 'b', &Quad2);
  valueVar(&ptrOutput, 'c', &Quad3);
  valueVar(&ptrOutput, 'd', &Quad4);
  *pptrOutput = ptrOutput;
}

void ecmFrontText(char *tofactorText, int doFactorization, char *knownFactors)
{
  char *ptrOutput;
  enum eExprErr rc;
  ptrInputText = tofactorText;
  if (batch)
  {
    BatchFactorization(tofactorText, doFactorization);
    return;
  }
  rc = ComputeExpression(tofactorText, 1, &tofactor);
  if (output == NULL)
  {
    output = (char *)malloc(1000000);
  }
  if (output == NULL)
  {
    return;   // Go out if cannot generate output string.
  }
  if (rc == EXPR_OK && doFactorization)
  {
    NumberLength = tofactor.nbrLimbs;
    CompressBigInteger(nbrToFactor, &tofactor);
#ifdef __EMSCRIPTEN__
    lModularMult = 0;
#endif
    factor(&tofactor, nbrToFactor, factorsMod, astFactorsMod, knownFactors);
  }
  ptrOutput = output;
  strcpy(output, "2<p>");
  ptrOutput += strlen(output);
  SendFactorizationToOutput(rc, astFactorsMod, &ptrOutput, doFactorization);
  if (rc == EXPR_OK && doFactorization)
  {
    if (tofactor.nbrLimbs > 1 || tofactor.limbs[0].x > 0)
    {      // Number to factor is not zero.
      GetNumberOfDivisors(&ptrOutput);
      GetSumOfDivisors(&ptrOutput);
      GetEulerTotient(&ptrOutput);
      GetMobius(&ptrOutput);
    }
    ComputeFourSquares(astFactorsMod);
    ShowFourSquares(&ptrOutput);
    strcpy(ptrOutput, lang ? "<p>Tiempo transcurrido: " : "<p>Time elapsed: ");
    ptrOutput += strlen(ptrOutput);
#ifdef __EMSCRIPTEN__
    GetDHMSt(&ptrOutput, (int)(tenths() - originalTenthSecond));
    strcpy(ptrOutput, "</p>");
    ptrOutput += strlen(ptrOutput);
#endif
  }
  strcpy(ptrOutput, lang ? "<p>" COPYRIGHT_SPANISH "</p>" :
    "<p>" COPYRIGHT_ENGLISH "</p>");
}

void doWork(char* data)
{
  int flags;
  char *ptrText;
  char *ptrData = data;
  char *ptrPower, *ptrMod;
  char *ptrWebStorage, *ptrKnownFactors;
  if (output == NULL)
  {
    output = malloc(3000000);
  }
  groupLen = 0;
  while (*ptrData != ',')
  {
    groupLen = groupLen * 10 + (*ptrData++ - '0');
  }
  ptrData++;             // Skip comma.
  flags = *ptrData;
  if (flags == '-')
  {
    flags = -*(++ptrData);
  }
  lang = flags & 1;
  ptrData += 2;          // Skip app number and second comma.
  batch = (*ptrData == '1');
  verbose = (*(ptrData + 1) == '1');
  prettyprint = (*(ptrData + 2) == '1');
  cunningham = (*(ptrData + 3) == '1');
  ptrData += 4;
  ptrWebStorage = ptrData + strlen(ptrData) + 1;
  ptrKnownFactors = findChar(ptrWebStorage, '=');
  if (prettyprint == 0)
  {
    groupLen = -groupLen;  // Do not show number of digts.
  }
  if (ptrKnownFactors)
  {
    ptrKnownFactors++;
  }
  if (flags & 0x80)
  {
    if (ptrKnownFactors)
    {
      ptrText = ptrKnownFactors + strlen(ptrKnownFactors) + 1;
      NextEC = 0;
      while (*ptrText != 0)
      {
        NextEC = NextEC * 10 + (*ptrText++ & 0x0F);
      }
      flags = 2;  // do factorization, no batch mode.
    }
  }
  ecmFrontText(ptrData, flags & 2, ptrKnownFactors); // The 3rd parameter includes known factors.
#ifdef __EMSCRIPTEN__
  databack(output);
#endif
}

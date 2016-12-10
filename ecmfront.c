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
static BigInteger tofactor;
static int nbrToFactor[MAX_LEN];
static int NbrFactorsMod = 0;
int groupLen = 6;
struct sFactors astFactorsMod[1000];
int factorsMod[10000];
static BigInteger factorValue, result;
BigInteger valueX;
char outputExpr[100000];
char *ptrInputText;

#ifdef __EMSCRIPTEN__
extern int NextEC;
void databack(char *data);
extern double originalTenthSecond;
int tenths(void);
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

static void SendFactorizationToOutput(enum eExprErr rc, char **pptrOutput, int doFactorization)
{
  char *ptrOutput = *pptrOutput;
  if (rc != EXPR_OK)
  {
    textError(output + 1, rc);
    ptrOutput = output + strlen(output);
  }
  else
  {
    Bin2Dec(tofactor.limbs, ptrOutput, tofactor.nbrLimbs, groupLen);
    ptrOutput += strlen(ptrOutput);
    if (doFactorization)
    {
      struct sFactors *pstFactors;
      int i = 0;
      pstFactors = &astFactorsMod[1];
      if (NbrFactorsMod == 1 && pstFactors->multiplicity == 1 &&
        (*pstFactors->ptrFactor > 1 || *(pstFactors->ptrFactor + 1) > 1))
      {    // Do not show zero or one as prime.
        strcpy(ptrOutput, lang ? " es primo" : " is prime");
        ptrOutput += strlen(ptrOutput);
      }
      else
      {
        strcpy(ptrOutput, " = ");
        ptrOutput += strlen(ptrOutput);
        for (;;)
        {
          UncompressBigInteger(pstFactors->ptrFactor, &factorValue);
          Bin2Dec(factorValue.limbs, ptrOutput, factorValue.nbrLimbs, groupLen);
          ptrOutput += strlen(ptrOutput);
          if (pstFactors->multiplicity > 1)
          {
            strcpy(ptrOutput, "<sup>");
            ptrOutput += strlen(ptrOutput);
            int2dec(&ptrOutput, pstFactors->multiplicity);
            strcpy(ptrOutput, "</sup>");
            ptrOutput += strlen(ptrOutput);
          }
          if (++i == NbrFactorsMod)
          {
            break;
          }
          strcpy(ptrOutput, " &times; ");
          ptrOutput += strlen(ptrOutput);
          pstFactors++;
        }
      }
    }
  }
  *pptrOutput = ptrOutput;
}

static void BatchFactorization(char *tofactorText, int doFactorization)
{
  char *ptrNextBatchFactor = tofactorText;
  char *ptrEndBatchFactor = tofactorText + strlen(tofactorText);
  char *ptrCurrBatchFactor;
  int counter;
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
          factor(nbrToFactor, factorsMod, astFactorsMod, NULL);
          NbrFactorsMod = astFactorsMod[0].multiplicity;
        }
        SendFactorizationToOutput(rc, &ptrOutput, doFactorization);
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
        factor(nbrToFactor, factorsMod, astFactorsMod, NULL);
        NbrFactorsMod = astFactorsMod[0].multiplicity;
      }
      SendFactorizationToOutput(rc, &ptrOutput, doFactorization);
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
  struct sFactors *pstFactor;
  int factorNumber;
  result.limbs[0].x = 1;    // Set result to 1.
  result.nbrLimbs = 1;
  result.sign = SIGN_POSITIVE;
  pstFactor = &astFactorsMod[1];
  for (factorNumber = 1; factorNumber <= astFactorsMod[0].multiplicity; factorNumber++)
  {
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
    addbigint(&Temp1, -1);   // p^(e-1)-1
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
  }
  strcpy(ptrOutput, "<p>Möbius: ");
  ptrOutput += strlen(ptrOutput);
  int2dec(&ptrOutput, mobius);
  strcpy(ptrOutput, "</p>");
  ptrOutput += strlen(ptrOutput);
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
    factor(nbrToFactor, factorsMod, astFactorsMod, knownFactors);
    NbrFactorsMod = astFactorsMod[0].multiplicity;
  }
  ptrOutput = output;
  strcpy(output, "2<p>");
  ptrOutput += strlen(output);
  SendFactorizationToOutput(rc, &ptrOutput, doFactorization);
  if (rc == EXPR_OK && doFactorization)
  {
    GetNumberOfDivisors(&ptrOutput);
    GetSumOfDivisors(&ptrOutput);
    GetEulerTotient(&ptrOutput);
    GetMobius(&ptrOutput);
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

#ifdef __EMSCRIPTEN__
void doWork(char* data, int size)
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
  ptrWebStorage = ptrData + strlen(ptrData) + 1;
  ptrKnownFactors = findChar(ptrWebStorage, '=');
  if (ptrKnownFactors)
  {
    ptrKnownFactors++;
  }
  if (flags & 0x80)
  {
    batch = 0;
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
  else
  {
    batch = flags & 4;
  }
  ecmFrontText(ptrData, flags & 2, ptrKnownFactors); // The 3rd parameter includes known factors.
  databack(output);
}
#endif

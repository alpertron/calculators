//
// This file is part of Alpertron Calculators.
//
// Copyright 2017-2021 Dario Alejandro Alpern
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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bignbr.h"
#include "expression.h"
#include "factor.h"
#include "showtime.h"
#include "batch.h"

static char *ptrEndBatchFactor;
static char *ptrCurrBatchFactor;
static char *ptrNextBatchFactor;
static bool firstExprProcessed;
bool fromFile;
bool lineEndingCRLF;
int valuesProcessed;
char outputExpr[200000];
#ifdef __EMSCRIPTEN__
char *ptrInputText;
char emptyInputText;
#endif
static char* ptrExprToProcess;
static char* ptrOutput;
static char* ptrSrcString;
static const char* ptrConditionExpr;
static bool errorDisplayed;
static int endValuesProcessed;
static pBatchCallback callback;
static const char* ptrStartExprs[MAX_EXPRESSIONS];
static int nbrExpressions;

void beginLine(char** pptrOutput)
{
  if (!fromFile)
  {
    copyStr(pptrOutput, "<p>");
  }
}

void finishLine(char** pptrOutput)
{
  if (fromFile)
  {
    copyStr(pptrOutput, (lineEndingCRLF ? "\r\n\r\n" : "\n\n"));
  }
  else
  {
    copyStr(pptrOutput, "</p>");
  }
}

static void stringToHTML(char **pptrOutput, const char *ptrString)
{
  const char* ptrStr = ptrString;
  int character;
  char *ptrOut = *pptrOutput;
  for (;;)
  {
    char c = *ptrStr;
    if (c == 0)
    {
      break;    // End of string, so go out.
    }
    if (((int)c & 0x80) == 0)
    {           // 1-byte UTF-8 character
      character = c;
      ptrStr++;
    }
    else if (((int)c & 0x60) == 0x40)
    {            // 2-byte UTF-8 character
      character = (((int)c & 0x1F) * 64) + (*(ptrStr + 1) & 0x3F);
      ptrStr += 2;
    }
    else
    {            // 3-byte UTF-8 character
      character = (((int)c & 0x1F) * 0x1000) +
        (((int)*(ptrStr + 1) & 0x3F) * 0x40) +
        ((int)*(ptrStr + 2) & 0x3F);
      ptrStr += 3;
    }
    if ((character >= ' ') && (character < 127) && (character != '<') &&
      (character != '>') && (character != '&'))
    {  // Safe to copy character.
      *ptrOut = (char)character;
      ptrOut++;
    }
    else
    {  // Insert HTML entity.
      *ptrOut = '&';
      ptrOut++;
      *ptrOut = '#';
      ptrOut++;
      int2dec(&ptrOut, character);
      *ptrOut = ';';
      ptrOut++;
    }
  }
  *pptrOutput = ptrOut;
}

static enum eExprErr evalExpression(const char *expr, BigInteger *ptrResult)
{
  const char *ptrInputExpr = expr;
  char *ptrOutputExpr = outputExpr;
  while (*ptrInputExpr != 0)
  {
    char c = *ptrInputExpr;
    if ((c == ';') || (c == ':'))
    { // End of expression characters for batch mode.
      break;
    }
    // Copy character to output expression.
    *ptrOutputExpr = c;
    ptrOutputExpr++;
    ptrInputExpr++;
  }
  *ptrOutputExpr = 0;   // Append string terminator.
  return ComputeExpression(outputExpr, ptrResult, true);
}

static void SkipSpaces(char **pptrText)
{
  char *ptrText = *pptrText;
  while ((*ptrText == ' ') || (*ptrText == 9))
  {  // Skip spaces or tabs.
    ptrText++;
  }
  *pptrText = ptrText;
}

static void BatchError(char **pptrOutput, const char *batchText, const char *errorText)
{
  char *ptrOut = *pptrOutput;
  stringToHTML(&ptrOut, batchText);
  *ptrOut = ':';
  ptrOut++;
  *ptrOut = ' ';
  ptrOut++;
  copyStr(&ptrOut, errorText);
  *pptrOutput = ptrOut;
  counterC = 0;
}

static bool doCallback(const char * ptrExpr, BigInteger *valueFound, int type)
{
  enum eExprErr rcode = evalExpression(ptrExpr, valueFound);
  if (rcode == EXPR_OK)
  {
    bool hexBak = hexadecimal;
    if ((type & BATCH_MASK_HEX) != 0)
    {
      hexadecimal = true;
    }
    else
    {
      hexadecimal = false;
    }
    callback(&ptrOutput, type);
    hexadecimal = hexBak;
  }
  else
  {
    textError(&ptrOutput, rcode);
    return true;
  }
  return false;
}
static bool ProcessLoop(bool* pIsBatch, const char* batchText, BigInteger* valueFound)
{
  char* ptrCharFound;
  const char* NextExpr;
  const char* EndExpr;
  char* ptrStartExpr;
  const char* ptrStartQuote;
  const char* ptrEndQuote;
#ifdef __EMSCRIPTEN__
  ptrInputText = &emptyInputText;
#endif
  ptrCharFound = findChar(ptrSrcString + 1, ';');
  if (ptrCharFound == NULL)
  {
    BatchError(&ptrOutput, batchText,
      lang ? "se esperaban tres o cuatro puntos y comas pero no hay ninguno" :
      "three or four semicolons expected but none found");
    ptrOutput += 4;
    return false;
  }
  ptrStartExpr = ptrSrcString + 1;
  SkipSpaces(&ptrStartExpr);
  if (*ptrStartExpr != '=')
  {
    BatchError(&ptrOutput, batchText,
      lang ? "falta signo igual en la primera expresión" :
      "equal sign missing in first expression");
    ptrOutput += 4;
    return false;
  }
  ptrCharFound = findChar(ptrSrcString + 1, ';');
  expressionNbr = 1;
  if (!firstExprProcessed)
  {
    enum eExprErr rc = evalExpression(ptrStartExpr + 1, valueFound);
    CopyBigInt(&valueX, valueFound);
    if (rc != EXPR_OK)
    {
      textError(&ptrOutput, rc);
      ptrOutput += 4;
      return false;
    }
  }
  ptrStartExpr = ptrCharFound + 1;
  SkipSpaces(&ptrStartExpr);
  if ((*ptrStartExpr != 'x') && (*ptrStartExpr != 'X'))
  {
    BatchError(&ptrOutput, batchText,
      lang ? "falta variable x en la segunda expresión" :
      "variable x missing in second expression");
    ptrOutput += 4;
    return false;
  }
  ptrStartExpr++;               // Skip variable 'x'.
  SkipSpaces(&ptrStartExpr);
  if (*ptrStartExpr != '=')
  {
    BatchError(&ptrOutput, batchText,
      lang ? "falta signo igual en la segunda expresión" :
      "equal sign missing in second expression");
    ptrOutput += 4;
    return false;
  }
  NextExpr = ptrStartExpr + 1;  // Skip equal sign.
  ptrCharFound = findChar(ptrStartExpr, ';');  // Find second semicolon.
  if (ptrCharFound == NULL)
  {      // Third semicolon not found.
    BatchError(&ptrOutput, batchText,
      lang ? "se esperaban tres o cuatro puntos y comas pero solo hay uno" :
      "three or four semicolons expected but there are only one");
    ptrOutput += 4;
    return false;
  }
  EndExpr = ptrCharFound + 1;  // Point to end expression.
  ptrCharFound = findChar(EndExpr, ';');  // Find third semicolon.
  if (ptrCharFound == NULL)
  {      // Third semicolon not found.
    BatchError(&ptrOutput, batchText,
      lang ? "se esperaban tres o cuatro puntos y comas pero solo hay dos" :
      "three or four semicolons expected but there are only two");
    ptrOutput += 4;
    return false;
  }
  ptrCharFound++;
  // Skip spaces.
  while ((*ptrCharFound != 0) && (*ptrCharFound <= ' '))
  {
    ptrCharFound++;
  }
  ptrExprToProcess = ptrCharFound;
  ptrStartQuote = NULL;
  ptrEndQuote = NULL;
  if (*ptrCharFound == '\"')
  {
    const char* ptrColon;
    nbrExpressions = 0;
    ptrStartQuote = ptrCharFound + 1;
    ptrEndQuote = ptrStartQuote;
    // Find closing quote.
    while ((*ptrEndQuote != 0) && (*ptrEndQuote != '\"'))
    {
      if (*ptrEndQuote == '%')
      {
        ptrEndQuote++;  // Discard character after percent sign.
      }
      ptrEndQuote++;
    }
    if (*ptrEndQuote == 0)
    {
      BatchError(&ptrOutput, batchText,
        lang ? "falta comilla de cierre" :
        "missing closing quote");
      ptrOutput += 4;
      return false;
    }
    // Find number of conversion clauses.
    while (ptrCharFound < ptrEndQuote)
    {
      if (*ptrCharFound == '%')
      {
        char upper = *(ptrCharFound + 1) & 0xDF;
        char upper2 = *(ptrCharFound + 2) & 0xDF;
        if ((upper == 'D') || (upper == 'X') || (upper == 'L'))
        {    // Decimal, hexadecimal or logic.
          ptrCharFound += 2;
          if (nbrExpressions == MAX_EXPRESSIONS)
          {
            BatchError(&ptrOutput, batchText,
              lang ? "demasiadas cláusulas de conversion" :
              "too many conversion clauses");
            ptrOutput += 4;
            return false;
          }
          nbrExpressions++;
        }
        else if ((*(ptrCharFound + 1) == '\"') || (*(ptrCharFound + 1) == '%'))
        {    // Quote or percent.
          ptrCharFound += 2;
        }
        else if ((upper == 'F') && ((upper2 == 'D') || (upper2 == 'X')))
        {    // Factoring with decimal or hexadecimal output.
          ptrCharFound += 3;
          if (nbrExpressions == MAX_EXPRESSIONS)
          {
            BatchError(&ptrOutput, batchText,
              lang ? "demasiadas cláusulas de conversion" :
              "too many conversion clauses");
            ptrOutput += 4;
            return false;
          }
          nbrExpressions++;
        }
        else
        {
          BatchError(&ptrOutput, batchText,
            lang ? "carácter extraño después de %" :
            "strange character after %");
          ptrOutput += 4;
          return false;
        }
      }
      else
      {
        ptrCharFound++;
      }
    }
    // Find start of expressions separated by colons.
    ptrColon = ptrEndQuote + 1;
    for (int nbrColons = 0; nbrColons < nbrExpressions; nbrColons++)
    {
      ptrColon = findChar(ptrColon, ':');
      if (ptrColon == NULL)
      {
        BatchError(&ptrOutput, batchText,
          lang ? "la cantidad de clásulas de conversión es mayor que la cantidad de dos puntos" :
          "the number of conversion clauses is greater than the number of colons");
        ptrOutput += 4;
        return false;
      }
      ptrColon++;
      ptrStartExprs[nbrColons] = ptrColon;
    }
    ptrConditionExpr = findChar(ptrColon, ';');  // Find optional fourth semicolon. 
    ptrColon = findChar(ptrColon, ':');
    if (ptrColon != NULL)
    {
      BatchError(&ptrOutput, batchText,
        lang ? "la cantidad de clásulas de conversión es menor que la cantidad de dos puntos" :
        "the number of conversion clauses is less than the number of colons");
      ptrOutput += 4;
      return false;
    }
  }
  else
  {      // No quotes.
    ptrConditionExpr = findChar(ptrExprToProcess, ';');  // Find optional fourth semicolon. 
  }
  firstExprProcessed = true;
  if (ptrConditionExpr != NULL)
  {
    ptrConditionExpr++;
  }
  while (ptrOutput < &output[(int)sizeof(output) - 200000])
  {      // Perform loop while there is space in output buffer.
    bool processExpression = true;
    static enum eExprErr rcode;
    expressionNbr = 3;
    rcode = evalExpression(EndExpr, valueFound);
    if (rcode != EXPR_OK)
    {
      textError(&ptrOutput, rcode);
      break;   // Cannot compute end expression, so go out.
    }
    if (BigIntIsZero(valueFound))
    {          // End expression result is zero: end of loop
      firstExprProcessed = false;
      break;
    }
    if (valuesProcessed >= endValuesProcessed)
    {
      output[0] = (fromFile ? 'A' : '6');  // Show Continue button.
      break;
    }
    if (ptrConditionExpr != NULL)
    {
      expressionNbr = 5;
      rcode = evalExpression(ptrConditionExpr, valueFound);
      if (rcode == EXPR_OK)
      {
        if (BigIntIsZero(valueFound))
        {   // Do not compute factor expression if condition is false.
          processExpression = false;
        }
      }
      else
      {
        textError(&ptrOutput, rcode);
        if ((rcode == EXPR_SYNTAX_ERROR) || (rcode == EXPR_VAR_OR_COUNTER_REQUIRED))
        {   // Do not show multiple errors.
          break;
        }
        processExpression = false;
      }
    }
    if (processExpression)
    {
      expressionNbr = 4;
      if (ptrStartQuote == NULL)
      {
        rcode = evalExpression(ptrExprToProcess, valueFound);
        if (rcode == EXPR_OK)
        {
          callback(&ptrOutput, BATCH_NO_QUOTE);
        }
        else
        {
          textError(&ptrOutput, rcode);
        }
      }
      else
      {
        int colonNbr = 0;
        const char* ptrInsideQuotes = ptrStartQuote;
        while (ptrInsideQuotes < ptrEndQuote)
        {
          if (*ptrInsideQuotes == '%')
          {
            switch (*(ptrInsideQuotes + 1))
            {
            case '%':
            case '\"':
              *ptrOutput = *(ptrInsideQuotes + 1);
              ptrOutput += 2;
              ptrInsideQuotes++;
              break;
            case 'd':
            case 'D':
              if (doCallback(ptrStartExprs[colonNbr], valueFound, BATCH_NO_PROCESS_DEC))
              {
                return false;
              }
              colonNbr++;
              ptrInsideQuotes += 2;
              break;
            case 'x':
            case 'X':
              if (doCallback(ptrStartExprs[colonNbr], valueFound, BATCH_NO_PROCESS_HEX))
              {
                return false;
              }
              colonNbr++;
              ptrInsideQuotes += 2;
              break;
            case 'L':
            case 'l':
              rcode = evalExpression(ptrStartExprs[colonNbr], valueFound);
              if (rcode == EXPR_OK)
              {
                if (BigIntIsZero(valueFound))
                {
                  copyStr(&ptrOutput, "no");
                }
                else
                {
                  copyStr(&ptrOutput, lang ? "sí": "yes");
                }
                colonNbr++;
                ptrInsideQuotes += 2;
              }
              else
              {
                textError(&ptrOutput, rcode);
              }
              break;
            default:
              switch (*(ptrInsideQuotes + 2))
              {
              case 'd':
              case 'D':
                if (doCallback(ptrStartExprs[colonNbr], valueFound, BATCH_PROCESS_DEC))
                {
                  return false;
                }
                colonNbr++;
                ptrInsideQuotes += 3;
                break;
              default:
                if (doCallback(ptrStartExprs[colonNbr], valueFound, BATCH_PROCESS_HEX))
                {
                  return false;
                }
                colonNbr++;
                ptrInsideQuotes += 3;
                break;
              }
            }
          }
          else
          {       // Not %. Copy character to output.
            *ptrOutput = *ptrInsideQuotes;
            ptrOutput++;
            ptrInsideQuotes++;
          }
        }
      }
      if ((rcode == EXPR_SYNTAX_ERROR) || (rcode == EXPR_VAR_OR_COUNTER_REQUIRED) ||
        (rcode == EXPR_VAR_IN_EXPRESSION))
      {   // Do not show multiple errors.
        errorDisplayed = true;
        break;
      }
    }
    expressionNbr = 2;
    rcode = evalExpression(NextExpr, valueFound);
    if (rcode != EXPR_OK)
    {
      textError(&ptrOutput, rcode);
      break;   // Cannot compute next expression, so go out.
    }
    CopyBigInt(&valueX, valueFound);
    counterC++;
    if ((counterC == 2) && (pIsBatch != NULL))
    {         // Not processing first line: Indicate batch processing.
      *pIsBatch = true;
    }
    if (processExpression)
    {
      if (fromFile)
      {
        copyStr(&ptrOutput, (lineEndingCRLF ? "\r\n" : "\n"));
      }
      else
      {
        copyStr(&ptrOutput, "</li><li>");
      }
      valuesProcessed++;
    }
  }
  return true;
}

enum eExprErr BatchProcessing(char *batchText, BigInteger *valueFound, char **pptrOutput,
  bool *pIsBatch, pBatchCallback batchCallback)
{
  enum eExprErr rc = EXPR_OK;
  errorDisplayed = false;
  ptrExprToProcess = NULL;
  ptrOutput = output;
  ptrConditionExpr = NULL;
  callback = batchCallback;
  if (fromFile)
  {
    copyStr(&ptrOutput, "B");
  }
  else
  {
    copyStr(&ptrOutput, "2<ul><li>");
  }
  endValuesProcessed = valuesProcessed + 1000;
  if (pIsBatch != NULL)
  {
    *pIsBatch = false;    // Indicate not batch processing in advance.
  }
  if (valuesProcessed == 0)
  {        // Start batch factorization.
    ptrCurrBatchFactor = batchText;
    ptrEndBatchFactor = batchText + strlen(batchText);
    firstExprProcessed = false;
  }
  for (; ptrCurrBatchFactor < ptrEndBatchFactor;
    ptrCurrBatchFactor += (int)strlen(ptrCurrBatchFactor) + 1)
  {  // Get next line.
    char c;
    if ((ptrCurrBatchFactor != batchText) && (pIsBatch != NULL))
    {         // Not processing first line: Indicate batch processing.
      *pIsBatch = true;
    }
    expressionNbr = 0;
    if (firstExprProcessed == false)
    {
      valueX.nbrLimbs = 0;     // Invalidate variable x and counter c.
      ptrNextBatchFactor = ptrCurrBatchFactor;
      while ((*ptrNextBatchFactor != '\0') && (*ptrNextBatchFactor != '\r') &&
            (*ptrNextBatchFactor != '\n'))
      {
        ptrNextBatchFactor++;
      }
      if ((*ptrNextBatchFactor == '\r') && (*(ptrNextBatchFactor+1) == '\n'))
      {
        lineEndingCRLF = true;
        *ptrNextBatchFactor = ' ';
        ptrNextBatchFactor++;
      }
      *ptrNextBatchFactor = 0;    // Indicate end of line.
      counterC = 1;
    }
    // Skip leading spaces.
    ptrSrcString = ptrCurrBatchFactor;
    SkipSpaces(&ptrSrcString);
    c = *ptrSrcString;
    if (c == 0)
    {   // Empty line.
      copyStr(&ptrOutput, "<br>");
    }
    else if (c == '#')
    {   // Copy comment to output, but convert non-safe characters to entities.
      *ptrOutput = '#';
      ptrOutput++;
      stringToHTML(&ptrOutput, ptrCurrBatchFactor + 1);
    }
    else if ((c == 'x') || (c == 'X'))
    {   // Loop format: x=<orig expr>; x=<next expr>; <end expr>; <expr to factor>[; <factor cond>]
      if (!ProcessLoop(pIsBatch, batchText, valueFound))
      {
        continue;
      }
      if (valuesProcessed >= endValuesProcessed)
      {
        break;
      }
    }
    else
    {       // Factor expression (not loop).
#ifdef __EMSCRIPTEN__
      ptrInputText = ptrSrcString;
#endif
      if (valuesProcessed >= endValuesProcessed)
      {
        output[0] = (fromFile ? 'A' : '6');  // Show Continue button.
      }
      rc = ComputeExpression(ptrSrcString, valueFound, false);
      if (rc == EXPR_OK)
      {
        callback(&ptrOutput, BATCH_NO_QUOTE);
      }
      else
      {
        if (rc == EXPR_VAR_IN_EXPRESSION)
        {
          return rc;
        }
        textError(&ptrOutput, rc);
      }
      counterC = -1;
      if (fromFile)
      {
        copyStr(&ptrOutput, (lineEndingCRLF ? "\r\n" : "\n"));
      }
      else
      {
        copyStr(&ptrOutput, "</li><li>");
      }
      valuesProcessed++;
    }
    if (counterC == 1)
    {
      if (fromFile)
      {
        copyStr(&ptrOutput, (lineEndingCRLF ? "\r\n" : "\n"));
      }
      else
      {
        copyStr(&ptrOutput, "</li>");
      }
    }
    if (ptrOutput >= &output[(int)sizeof(output) - 200000])
    {
      output[0] = (fromFile ? 'A' : '6');  // Show Continue button.
      break;
    }
    if ((*ptrCurrBatchFactor & 0xDF) == 'x')
    {      // Loop mode.
      if (ptrConditionExpr != NULL)
      {
        ptrCurrBatchFactor += strlen(ptrConditionExpr);
        ptrConditionExpr = NULL;
      }
      else
      {
        ptrCurrBatchFactor += strlen(ptrExprToProcess);
      }
      ptrCurrBatchFactor++;
    }
    valueX.nbrLimbs = 0;       // Invalidate variable x and counter c.
  }
  if (!fromFile)
  {
    ptrOutput -= 4;            // Erase start tag <li> without contents.
  }
  if (counterC == 1)
  {
    if (!fromFile)
    {
      ptrOutput--;             // Erase extra character of </li>.
    }
    if (!errorDisplayed)
    {
      copyStr(&ptrOutput, lang ? "No hay valores para la expresión ingresada." :
        "There are no values for the requested expression.");
    }
  }
  if (fromFile)
  {
    copyStr(&ptrOutput, (lineEndingCRLF ? "\r\n" : "\n"));
  }
  else
  {
    copyStr(&ptrOutput, "</ul>");
  }
  *pptrOutput = ptrOutput;
  return rc;
}

char *findChar(const char *str, char c)
{
  const char* ptrStr = str;
  while (*ptrStr != '\0')
  {
    if (*ptrStr == c)
    {
      return (char *)ptrStr;
    }
    ptrStr++;
  }
  return NULL;
}

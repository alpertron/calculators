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
static const char* NextExpr;
static const char* EndExpr;
static const char* ptrStartExpr;
static const char* ptrStartQuote;
static const char* ptrEndQuote;

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

static void showErrorInExpr(char** pptrOutput, int rcode)
{
  copyStr(pptrOutput, lang ? "Error en la expresión " :
    "Error in expression #");
  int2dec(pptrOutput, expressionNbr);
  copyStr(pptrOutput, ": ");
  textError(pptrOutput, rcode);
  counterC = 0;
  ptrOutput += 4;
}

static enum eExprErr generateRPN(const char* expr, char **ptrRPN, bool varsExpected)
{
  enum eExprErr retcode;
  const char* ptrInputExpr = expr;
  char* ptrOutputExpr = outputExpr;
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
  retcode = convertToRPN(expr, ptrRPN, varsExpected);
  if (retcode != EXPR_OK)
  {
    showErrorInExpr(&ptrOutput, retcode);
  }
  return retcode;
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

static void SkipSpacesConst(const char** pptrText)
{
  const char* ptrText = *pptrText;
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
  ptrOut += 4;
  *pptrOutput = ptrOut;
  counterC = 0;
}

static bool doCallback(const char * ptrExpr, BigInteger *valueFound, int type)
{
  enum eExprErr rcode = ComputeExpression(ptrExpr,
    valueFound);
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
    showErrorInExpr(&ptrOutput, rcode);
    return true;
  }
  return false;
}

static const char* missingVariableES[2] =
{
  "falta variable x en la primera expresión",
  "falta variable x en la segunda expresión",
};

static const char* missingVariableEN[2] =
{
  "variable x missing in first expression",
  "variable x missing in second expression",
};

static const char *missingSemicolonES[2] =
{
  "se esperaban tres o cuatro puntos y comas pero no hay ninguno",
  "se esperaban tres o cuatro puntos y comas pero solo hay uno",
};

static const char* missingSemicolonEN[2] =
{
  "three or four semicolons expected but none found",
  "three or four semicolons expected but there are only one",
};

static const char *missingEqualSignES[2] =
{
  "falta signo igual en la primera expresión",
  "falta signo igual en la segunda expresión",
};

static const char* missingEqualSignEN[2] =
{
  "equal sign missing in first expression",
  "equal sign missing in second expression",
};

static bool convertExpressionsToRPN(const char* batchText)
{
  enum eExprErr retcode;
  char* ptrCharFound;
  const char* ptrBatchText = batchText;
  char* ptrRPN;
#ifdef __EMSCRIPTEN__
  ptrInputText = &emptyInputText;
#endif
  // Process first and second expressions "x = expr"
  // separated by semicolons.
  for (int exprNbr = 0; exprNbr < 2; exprNbr++)
  {
    SkipSpacesConst(&ptrBatchText);
    if ((*ptrBatchText != 'x') && (*ptrBatchText != 'X'))
    {
      BatchError(&ptrOutput, batchText,
        lang ? missingVariableES[exprNbr] : missingVariableEN[exprNbr]);
      return false;
    }
    ptrBatchText++;               // Skip variable 'x'.
    SkipSpacesConst(&ptrBatchText);
    if (*ptrBatchText != '=')
    {
      BatchError(&ptrOutput, batchText,
        lang ? missingEqualSignES[exprNbr] : missingEqualSignEN[exprNbr]);
      return false;
    }
    ptrCharFound = findChar(ptrBatchText + 1, ';');
    if (ptrCharFound == NULL)
    {
      BatchError(&ptrOutput, batchText,
        lang ? missingSemicolonES[exprNbr] : missingSemicolonEN[exprNbr]);
      return false;
    }
    ptrBatchText++;   // Point after equal (assignment) sign.
    expressionNbr = 1;
    if (exprNbr == 0)
    {  // Generate RPN for first expression "x = expr".
      expressionNbr = 1;
      retcode = generateRPN(ptrBatchText, &ptrRPN, false);
      ptrStartExpr = (const char*)ptrRPN;
    }
    else
    {  // Generate RPN for second expression "x = expr".
      expressionNbr = 2;
      retcode = generateRPN(ptrBatchText, &ptrRPN, true);
      NextExpr = (const char*)ptrRPN;
    }
    if (retcode != EXPR_OK)
    {
      return false;
    }
    ptrBatchText = ptrCharFound + 1;
  }
  EndExpr = ptrCharFound + 1;  // Point to end expression.
  ptrCharFound = findChar(EndExpr, ';');  // Find third semicolon.
  if (ptrCharFound == NULL)
  {      // Third semicolon not found.
    BatchError(&ptrOutput, batchText,
      lang ? "se esperaban tres o cuatro puntos y comas pero solo hay dos" :
      "three or four semicolons expected but there are only two");
    return false;
  }
  // Generate RPN for third expression (end expression).
  expressionNbr = 3;
  retcode = generateRPN(EndExpr, &ptrRPN, true);
  EndExpr = (const char*)ptrRPN;
  if (retcode != EXPR_OK)
  {
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
        if (*ptrEndQuote == 0)
        {
          break;        // Go out of loop if sequence % NULL is found.
        }
      }
      ptrEndQuote++;
    }
    if (*ptrEndQuote == 0)
    {
      BatchError(&ptrOutput, batchText,
        lang ? "falta comilla de cierre" :
        "missing closing quote");
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
            return false;
          }
          nbrExpressions++;
        }
        else
        {
          BatchError(&ptrOutput, batchText,
            lang ? "carácter extraño después de %" :
            "strange character after %");
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
        return false;
      }
      ptrColon++;
      expressionNbr = 4;
      retcode = generateRPN(ptrColon, &ptrRPN, true);
      ptrStartExprs[nbrColons] = (const char*)ptrRPN;
      if (retcode != EXPR_OK)
      {
        return false;
      }
    }
    ptrConditionExpr = findChar(ptrColon, ';');  // Find optional fourth semicolon. 
    ptrColon = findChar(ptrColon, ':');
    if (ptrColon != NULL)
    {
      BatchError(&ptrOutput, batchText,
        lang ? "la cantidad de clásulas de conversión es menor que la cantidad de dos puntos" :
        "the number of conversion clauses is less than the number of colons");
      return false;
    }
  }
  else
  {      // No quotes.
    ptrConditionExpr = findChar(ptrExprToProcess, ';');  // Find optional fourth semicolon. 
    expressionNbr = 4;
    retcode = generateRPN(ptrExprToProcess, &ptrExprToProcess, true);
    if (retcode != EXPR_OK)
    {
      return false;
    }
  }
  if (ptrConditionExpr != NULL)
  {
    ptrConditionExpr++;
    expressionNbr = 5;
    retcode = generateRPN(ptrConditionExpr, &ptrRPN, true);
    ptrConditionExpr = (const char*)ptrRPN;
    if (retcode != EXPR_OK)
    {
      return false;
    }
  }
  return true;
}

static bool InternalProcessLoop(bool* pIsBatch, const char* batchText,
  BigInteger* valueFound)
{
  bool retcode = convertExpressionsToRPN(batchText);
  if (!retcode)
  {
    return false;
  }
  if (!firstExprProcessed)
  {
    enum eExprErr rc = ComputeExpression(ptrStartExpr,
      valueFound);
    if (rc != EXPR_OK)
    {
      expressionNbr = 1;
      showErrorInExpr(&ptrOutput, rc);
      return false;
    }
    CopyBigInt(&valueX, valueFound);
    firstExprProcessed = true;
  }
  while (ptrOutput < &output[(int)sizeof(output) - 200000])
  {      // Perform loop while there is space in output buffer.
    bool processExpression = true;
    static enum eExprErr rcode;
    expressionNbr = 3;
    rcode = ComputeExpression(EndExpr, valueFound);
    if (rcode != EXPR_OK)
    {
      showErrorInExpr(&ptrOutput, rcode);
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
      rcode = ComputeExpression(ptrConditionExpr,
        valueFound);
      if (rcode == EXPR_OK)
      {
        if (BigIntIsZero(valueFound))
        {   // Do not compute factor expression if condition is false.
          processExpression = false;
        }
      }
      else
      {
        showErrorInExpr(&ptrOutput, rcode);
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
        rcode = ComputeExpression(ptrExprToProcess, valueFound);
        if (rcode == EXPR_OK)
        {
          callback(&ptrOutput, BATCH_NO_QUOTE);
        }
        else
        {
          showErrorInExpr(&ptrOutput, rcode);
          break;
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
              ptrOutput++;
              ptrInsideQuotes += 2;
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
              rcode = ComputeExpression(ptrStartExprs[colonNbr], valueFound);
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
                showErrorInExpr(&ptrOutput, rcode);
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
    rcode = ComputeExpression(NextExpr, valueFound);
    if (rcode != EXPR_OK)
    {
      showErrorInExpr(&ptrOutput, rcode);
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

static bool ProcessLoop(bool* pIsBatch, const char* batchText, BigInteger* valueFound)
{
  bool rc;
  setInsideExpressionLoop(true);
  rc = InternalProcessLoop(pIsBatch, batchText, valueFound);
  setInsideExpressionLoop(false);
  return rc;
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
      if (fromFile)
      {
        copyStr(&ptrOutput, (lineEndingCRLF ? "\r\n" : "\n"));
      }
      else
      {
        copyStr(&ptrOutput, "<li>");
      }
    }
    else if ((c == 'x') || (c == 'X'))
    {   // Loop format: x=<orig expr>; x=<next expr>; <end expr>; <expr to factor>[; <factor cond>]      
      if (!ProcessLoop(pIsBatch, ptrCurrBatchFactor, valueFound))
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
      rc = ComputeExpression(ptrSrcString, valueFound);
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
      return (char *)ptrStr;  // Pointer to character found.
    }
    ptrStr++;
  }
  return NULL;                // Character not found.
}

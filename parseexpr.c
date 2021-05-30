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

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include "bignbr.h"
#include "polynomial.h"

static char RPNbuffer[1000000];
static bool forceMultiplication;
static short stackOper[STACK_OPER_SIZE];
static char stackArgumNbrPriority[STACK_OPER_SIZE];  // For functions: argument number
                                                     // For operators: priority.
static int stackOperIndex;
static int exponOperatorCounter;
static bool prevTokenIsNumber;

static bool isFunc(const char** ppcInput, const struct sFuncOperExpr** ppstFuncOperExpr)
{
  const struct sFuncOperExpr* pstFuncOperExpr = *ppstFuncOperExpr;
  while (pstFuncOperExpr->name != NULL)
  {
    const char* pcInput = *ppcInput;
    const char* ptrFuncName = pstFuncOperExpr->name;
    while (*ptrFuncName != '\0')
    {
      if ((*ptrFuncName >= 'A') && (*ptrFuncName < 'Z'))
      {           // Function name character is a letter.
        if ((*pcInput != *ptrFuncName) && ((*pcInput - 0x20) != *ptrFuncName))
        {         // Insensitive case comparison.
          break;
        }
      }
      else if (*pcInput != *ptrFuncName)
      {           // Not a letter. Perform exact comparison.
        break;    // Function name not found.
      }
      else
      {           // Nothing to do.
      }
      pcInput++;
      ptrFuncName++;
    }
    if (*ptrFuncName == '\0')
    {
      *ppcInput = pcInput;
      *ppstFuncOperExpr = pstFuncOperExpr;
      return true;         // Function or operator name was found.
    }
    pstFuncOperExpr++;
  }
  return false;            // Function or operator name not found.
}

static BigInteger value;

static void getHexValue(const char** pptrInput)
{
  const char* ptrFirstHexDigit = NULL;
  const char* ptrInput = *pptrInput + 2;   // Skip prefix 0x.
  while (((*ptrInput >= '0') && (*ptrInput <= '9')) ||
    ((*ptrInput >= 'A') && (*ptrInput <= 'F')) ||
    ((*ptrInput >= 'a') && (*ptrInput <= 'f')))
  {                // Find end of number.
    if ((ptrFirstHexDigit == NULL) && (*ptrInput != '0'))
    {              // First non-zero digit.
      ptrFirstHexDigit = ptrInput;
    }
    ptrInput++;
  }
  if (ptrFirstHexDigit == NULL)
  {    // Number is zero
    value.limbs[0].x = 0;
    value.nbrLimbs = 1;
  }
  else
  {    // Generate big integer from hexadecimal number from right to left.
    unsigned int unsignedLimb;
    unsigned int carry = 0;
    int shLeft = 0;
    limb* ptrLimb = &value.limbs[0];
    for (const char* ptrHexNbr = ptrInput - 1; ptrHexNbr >= ptrFirstHexDigit; ptrHexNbr--)
    {
      char c = *ptrHexNbr;
      if ((c >= '0') && (c <= '9'))
      {
        c -= '0';
      }
      else if ((c >= 'A') && (c <= 'F'))
      {
        c -= 'A' - 10;
      }
      else
      {
        c -= 'a' - 10;
      }
      carry += (unsigned int)c << shLeft;
      shLeft += 4;   // 4 bits per hex digit.
      if (shLeft >= BITS_PER_GROUP)
      {
        shLeft -= BITS_PER_GROUP;
        unsignedLimb = carry & MAX_VALUE_LIMB;
        ptrLimb->x = (int)unsignedLimb;
        ptrLimb++;
        carry = (unsigned int)c >> (4 - shLeft);
      }
    }
    if ((carry != 0) || (ptrLimb == &value.limbs[0]))
    {
      ptrLimb->x = carry;
      ptrLimb++;
    }
    unsignedLimb = (unsigned int)(ptrLimb - &value.limbs[0]);
    value.nbrLimbs = (int)unsignedLimb;
  }
  *pptrInput = ptrInput;
}

static enum eExprErr parseNumberInsideExpr(const char** ppInput, char** ppOutput)
{
  (void)exponOperatorCounter;
  const char* pInput = *ppInput;
  char* ptrOutput = *ppOutput;
  const char* ptrInput = pInput;
  if ((*ptrInput == '0') && ((*(ptrInput + 1) == 'x') || (*(ptrInput + 1) == 'X')))
  {              // Hexadecimal number.
    getHexValue(&ptrInput);
  }
  else
  {              // Decimal number.
    while ((*ptrInput >= '0') && (*ptrInput <= '9'))
    {                // Find end of number.
      ptrInput++;
    }
    Dec2Bin(pInput - 1, value.limbs, (int)(ptrInput + 1 - pInput), &value.nbrLimbs);
  }
  pInput = ptrInput;
  value.sign = SIGN_POSITIVE;
  *ptrOutput = TOKEN_NUMBER;
  ptrOutput++;
  while (value.nbrLimbs > 1)
  {
    if (value.limbs[value.nbrLimbs - 1].x != 0)
    {
      break;
    }
    value.nbrLimbs--;
  }
  *ptrOutput = (char)((unsigned int)value.nbrLimbs >> 8);
  ptrOutput++;
  *ptrOutput = (char)value.nbrLimbs;
  ptrOutput++;
  for (int index = 0; index < value.nbrLimbs; index++)
  {
    unsigned int limbValue = (unsigned int)value.limbs[index].x;
    for (int bitNbr = BITS_PER_GROUP; bitNbr > 0; bitNbr -= 8)
    {
      *ptrOutput = (char)limbValue;
      ptrOutput++;
      limbValue >>= 8;
    }
  }
  *ppInput = pInput;
  *ppOutput = ptrOutput;
  return EXPR_OK;
}

#ifdef POLYEXPR
static enum eExprErr parseNumberInsidePolyExpr(const char** ppInput, char** ppOutput)
{
  const char* pInput = *ppInput;
  char* ptrOutput = *ppOutput;
  const char* ptrInput = pInput;
  int limb;
  int bitNbr;
  while ((*ptrInput >= '0') && (*ptrInput <= '9'))
  {        // Find end of number.
    ptrInput++;
  }
  Dec2Bin(pInput - 1, value.limbs, (int)(ptrInput + 1 - pInput), &value.nbrLimbs);
  pInput = ptrInput;
  value.sign = SIGN_POSITIVE;
  if (exponOperatorCounter != 0)
  {
    if (value.nbrLimbs != 1)
    {
      return EXPR_EXPONENT_TOO_LARGE;
    }
    *ptrOutput = TOKEN_NUMBER;
    ptrOutput++;
    limb = value.limbs[0].x;
    for (bitNbr = BITS_PER_GROUP; bitNbr > 0; bitNbr -= 8)
    {
      *ptrOutput = (char)limb;
      ptrOutput++;
      limb >>= 8;
    }
  }
  else
  {
    if (!modulusIsZero)
    {
      // Reduce mod power outside exponent.
      (void)BigIntRemainder(&value, &powerMod, &value);
      if (value.nbrLimbs < NumberLength)
      {    // Fill with zeros.
        (void)memset(&value.limbs[value.nbrLimbs], 0, (NumberLength - value.nbrLimbs) * sizeof(int));
      }
      // Convert to Montgomery notation.
      modmult(value.limbs, MontgomeryMultR2, value.limbs);
      value.nbrLimbs = NumberLength;
    }
    *ptrOutput = TOKEN_NUMBER;
    ptrOutput++;
    while (value.nbrLimbs > 1)
    {
      if (value.limbs[value.nbrLimbs - 1].x != 0)
      {
        break;
      }
      value.nbrLimbs--;
    }
    *ptrOutput = (char)(value.nbrLimbs >> 8);
    ptrOutput++;
    *ptrOutput = (char)value.nbrLimbs;
    ptrOutput++;
    for (int index = 0; index < value.nbrLimbs; index++)
    {
      limb = value.limbs[index].x;
      for (bitNbr = BITS_PER_GROUP; bitNbr > 0; bitNbr -= 8)
      {
        *ptrOutput = (char)limb;
        ptrOutput++;
        limb >>= 8;
      }
    }
  }
  *ppInput = pInput;
  *ppOutput = ptrOutput;
  return EXPR_OK;
}
#endif
static enum eExprErr processClosingParenOrComma(char **ppOutput, char c,
  enum eParseExpr eParseExpr)
{
  char* ptrOutput = *ppOutput;
  char s = '\0'; // Assume parenthesis mismatch.
  while (stackOperIndex > 0)
  {      // Send operators to output.
    stackOperIndex--;
    s = (char)stackOper[stackOperIndex];
    if (s == OPER_PAREN)
    {    // Operator on stack has less precedence.
      break;
    }
    if (s == OPER_POWER)
    {
      if (eParseExpr == PARSE_EXPR_POLYNOMIAL)
      {
        exponOperatorCounter--;
      }
      if (exponOperatorCounter == 0)
      {
        *ptrOutput = TOKEN_END_EXPON;
        ptrOutput++;
      }
      if (eParseExpr == PARSE_EXPR_POLYNOMIAL)
      {
        *ptrOutput = s;
        ptrOutput++;
      }
    }
    else
    {
      *ptrOutput = s;
      ptrOutput++;
    }
  }
  if (s != OPER_PAREN)
  {     // Parentheses mismatch.
    return EXPR_PAREN_MISMATCH;
  }
  if (c == ',')
  {
    if ((stackOperIndex <= 0) ||
      (((unsigned short)stackOper[stackOperIndex - 1] & 0xFF00U) == 0U))
    {           // Previous element in stack is not a function token.
      return EXPR_SYNTAX_ERROR;
    }
    // Increment number of arguments given by user.
    stackArgumNbrPriority[stackOperIndex - 1]++;
    // Check whether the user provided extra arguments.
    if (stackArgumNbrPriority[stackOperIndex - 1] ==
      ((unsigned short)stackOper[stackOperIndex - 1] >> 8))
    {
      return EXPR_TOO_MANY_ARGUMENTS;
    }
    stackOper[stackOperIndex] = s;  // Push back paren.
    stackOperIndex++;
    prevTokenIsNumber = false;
  }
  else
  {
    if ((stackOperIndex > 0) &&
      (((unsigned short)stackOper[stackOperIndex - 1] & 0xFF00U) != 0U))
    {           // Previous element in stack is a function token.
      if ((stackArgumNbrPriority[stackOperIndex - 1] + 1) !=
        ((unsigned short)stackOper[stackOperIndex - 1] >> 8))
      {
        return EXPR_TOO_FEW_ARGUMENTS;
      }
      stackOperIndex--;
      *ptrOutput = (char)stackOper[stackOperIndex];
      ptrOutput++;
    }
    prevTokenIsNumber = true;
  }
  *ppOutput = ptrOutput;
  return EXPR_OK;
}

static enum eExprErr parsePrevTokenIsNumber(const char** ppInput, char** ppOutput,
  const struct sFuncOperExpr* binaryOperExpr, const struct sFuncOperExpr* rightFuncExpr,
  enum eParseExpr eParseExpr)
{
  enum eExprErr rc;
  const char* pInput = *ppInput;
  char* ptrOutput = *ppOutput;
  char c = *pInput;
  const struct sFuncOperExpr* pstOperatorExpr = rightFuncExpr;
  if (isFunc(&pInput, &pstOperatorExpr))
  {           // Right function found.
    *ptrOutput = (char)(pstOperatorExpr->token);
    ptrOutput++;
    *ppInput = pInput;
    *ppOutput = ptrOutput;
    return EXPR_OK;
  }
  prevTokenIsNumber = false;
  pstOperatorExpr = binaryOperExpr;
  if (isFunc(&pInput, &pstOperatorExpr) || forceMultiplication)
  {              // Operator name was found.
    short oper;
    char priority;
    char s;
    bool isInfix = false;
    if (forceMultiplication)
    {
      priority = 2;
      oper = OPER_MULTIPLY;
    }
    else
    {
      priority = pstOperatorExpr->priority;
      oper = pstOperatorExpr->token;
      if (priority < 0)
      {
        priority = -priority;
        isInfix = true;
      }
    }
    while (stackOperIndex > 0)
    {      // Send operators that have more precedence to output.
      stackOperIndex--;
      if (priority < stackArgumNbrPriority[stackOperIndex])
      {    // Operator on stack has less precedence.
        stackOperIndex++;
        break;
      }
      s = (char)(stackOper[stackOperIndex]);
      if (s == OPER_POWER)
      {
        if (oper == OPER_POWER)
        {  // Power is right associative and maximum precedence
          stackOperIndex++;
          break;
        }
        if (eParseExpr == PARSE_EXPR_POLYNOMIAL)
        {
          exponOperatorCounter--;
        }
        if (exponOperatorCounter == 0)
        {
          *ptrOutput = TOKEN_END_EXPON;
          ptrOutput++;
        }
        if (eParseExpr == PARSE_EXPR_POLYNOMIAL)
        {
          *ptrOutput = s;
          ptrOutput++;
        }
      }
      else
      {
        *ptrOutput = s;
        ptrOutput++;
      }
    }
    if (isInfix)
    {
      *ptrOutput = (char)oper + 1;       // Insert infix operator.
      ptrOutput++;
    }
    stackOper[stackOperIndex] = oper;    // Push operator onto stack.
    stackArgumNbrPriority[stackOperIndex] = priority;
    stackOperIndex++;
    if (oper == OPER_POWER)
    {
      if (exponOperatorCounter == 0)
      {
        *ptrOutput = TOKEN_START_EXPON;
        ptrOutput++;
      }
      if (eParseExpr == PARSE_EXPR_POLYNOMIAL)
      {
        exponOperatorCounter++;
      }
    }
  }
  else if ((c == ')') || (c == ','))
  {
    pInput++;    // Point to character ater paren or comma.
    rc = processClosingParenOrComma(&ptrOutput, c, eParseExpr);
    if (rc != EXPR_OK)
    {
      return rc;
    }
  }
  else
  {
    return EXPR_SYNTAX_ERROR;
  }
  *ppInput = pInput;
  *ppOutput = ptrOutput;
  return EXPR_OK;
}

// Convert to Reverse Polish Notation using Shunting-yard algorithm
// If the number is part of the exponent it has to be reduced by
// using the modulus phi = (p-1)*p^(k-1), otherwise the modulus will
// be power = p^k.
int ConvertToReversePolishNotation(const char* input, char** pptrOut,
  const struct sFuncOperExpr* funcOperExpr, enum eParseExpr eParseExpr,
  bool *pUsingVariables)
{
  const struct sFuncOperExpr* leftFuncExpr;
  const struct sFuncOperExpr* rightFuncExpr;
  const struct sFuncOperExpr* unaryOperExpr;
  const struct sFuncOperExpr* binaryOperExpr;
  const struct sFuncOperExpr* pstFuncOperExpr;
  char* ptrOutput = RPNbuffer;
  const char* pInput = input;
  if (pUsingVariables != NULL)
  {
    *pUsingVariables = false;
  }
  prevTokenIsNumber = false;
  exponOperatorCounter = 0;
  char s;
  char variableLetter = ' ';  // Indicate variable letter not known yet.
  stackOperIndex = 0;
  *pptrOut = RPNbuffer;
  // Initialize pointers to tokens. Sections are separated by elements
  // whose field "name" is equal to NULL.
  pstFuncOperExpr = funcOperExpr;
  leftFuncExpr = pstFuncOperExpr;
  while (pstFuncOperExpr->name != NULL)
  {
    pstFuncOperExpr++;
  }
  pstFuncOperExpr++;
  rightFuncExpr = pstFuncOperExpr;
  while (pstFuncOperExpr->name != NULL)
  {
    pstFuncOperExpr++;
  }
  pstFuncOperExpr++;
  unaryOperExpr = pstFuncOperExpr;
  while (pstFuncOperExpr->name != NULL)
  {
    pstFuncOperExpr++;
  }
  pstFuncOperExpr++;
  binaryOperExpr = pstFuncOperExpr;

  while ((*pInput != '\0') && (*pInput != ';'))
  {
    const char* inputTemp;
    char c = *pInput;
    char cUppercase;
    char cUppercaseNext;
    forceMultiplication = false;
    if ((c == ' ') || (c == 9))
    {          // Ignore any spaces and tabs.
      pInput++;
      continue;
    }
    inputTemp = pInput;
    pstFuncOperExpr = leftFuncExpr;
    if (isFunc(&inputTemp, &pstFuncOperExpr))
    {              // Function name was found.
      if (prevTokenIsNumber == false)
      {
        pInput = inputTemp;
        stackOper[stackOperIndex] = pstFuncOperExpr->token;  // Push token onto stack.
        stackArgumNbrPriority[stackOperIndex] = 0;    // Indicate no arguments found yet.
        stackOperIndex++;
        continue;
      }
      // No operand between number and function. Force multiplication.
      forceMultiplication = true;
    }
    else
    {           // Nothing to do.
    }
    cUppercase = c & 0xDF;
    cUppercaseNext = *(pInput+1) & 0xDF;
    if ((cUppercase >= 'A') && (cUppercase <= 'Z') &&
      ((cUppercaseNext < 'A') || (cUppercaseNext > 'Z')))
    {           // Unique letter found.
      if ((eParseExpr == PARSE_EXPR_POLYNOMIAL) &&
        (variableLetter != cUppercase) && (variableLetter != ' '))
      {
        return EXPR_MULTIPLE_VARIABLES_NOT_ACCEPTED;
      }
      if ((eParseExpr == PARSE_EXPR_POLYNOMIAL) ||
        ((eParseExpr == PARSE_EXPR_GAUSSIAN) && (cUppercase == 'I')) ||
        ((eParseExpr == PARSE_EXPR_INTEGER) && ((cUppercase == 'X') || (cUppercase == 'C'))))
      {
        if (pUsingVariables != NULL)
        {
          *pUsingVariables = true;
        }
        variableLetter = cUppercase;    // Convert to uppercase.
      }
    }
    if (prevTokenIsNumber)
    {
      if (c == '(')
      {
        forceMultiplication = true;
      }
      if (cUppercase == variableLetter)
      {
        forceMultiplication = true;
      }
      enum eExprErr rc = parsePrevTokenIsNumber(&pInput, &ptrOutput, 
        binaryOperExpr, rightFuncExpr, eParseExpr);
      if (rc != EXPR_OK)
      {
        return rc;
      }
    }
    else
    {    // Not a number before this token.
      pstFuncOperExpr = unaryOperExpr;
      if (isFunc(&pInput, &pstFuncOperExpr))
      {          // Unary operator name was found.
                 // Push operator onto stack.
        stackOper[stackOperIndex] = pstFuncOperExpr->token;
        stackOperIndex++;
        prevTokenIsNumber = false;
        continue;
      }
      pInput++;
      if (c == '+')
      {          // Unary plus, nothing to do.

      }
      else if (c == '(')
      {          // Open parenthesis.
        stackOper[stackOperIndex] = OPER_PAREN;    // Push operator onto stack.
        stackArgumNbrPriority[stackOperIndex] = 0x7F; // Minimum priority
        stackOperIndex++;
      }
      else if (cUppercase == variableLetter)
      {
        if ((eParseExpr == PARSE_EXPR_POLYNOMIAL) && (exponOperatorCounter != 0))
        {
          return EXPR_CANNOT_USE_X_IN_EXPONENT;
        }
        else if ((eParseExpr == PARSE_EXPR_INTEGER) && (cUppercase == 'C'))
        {
          *ptrOutput = TOKEN_COUNTER;
        }
        else
        {
          *ptrOutput = TOKEN_VAR;
        }
        ptrOutput++;
        prevTokenIsNumber = true;
      }
      else if ((c >= '0') && (c <= '9'))
      {          // Number.
        enum eExprErr retcode;
#ifdef POLYEXPR
        if (eParseExpr == PARSE_EXPR_POLYNOMIAL)
        {
          retcode = parseNumberInsidePolyExpr(&pInput, &ptrOutput);
        }
        else
#endif
        {
          retcode = parseNumberInsideExpr(&pInput, &ptrOutput);
        }
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
        prevTokenIsNumber = true;
      }
      else 
      {
        return EXPR_SYNTAX_ERROR;
      }
    }
  }
  if (prevTokenIsNumber == false)
  {
    return EXPR_SYNTAX_ERROR;
  }
  while (stackOperIndex > 0)
  {      // Send operators to output.
    stackOperIndex--;
    s = (char)stackOper[stackOperIndex];
    if (s == OPER_POWER)
    {
      if (eParseExpr == PARSE_EXPR_POLYNOMIAL)
      {
        exponOperatorCounter--;
      }
      if (exponOperatorCounter == 0)
      {
        *ptrOutput = TOKEN_END_EXPON;
        ptrOutput++;
      }
      if (eParseExpr == PARSE_EXPR_POLYNOMIAL)
      {
        *ptrOutput = s;
        ptrOutput++;
      }
    }
    else if (s == OPER_PAREN)
    {    // Operator on stack has less precedence.
      return EXPR_PAREN_MISMATCH;
    }
    else
    {
      *ptrOutput = s;
      ptrOutput++;
    }
  }
  *ptrOutput = '\0';
  return EXPR_OK;
}

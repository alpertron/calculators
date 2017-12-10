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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bignbr.h"
#include "expression.h"
#define PAREN_STACK_SIZE 100
static BigInteger stackValues[PAREN_STACK_SIZE];
static char stackOperators[PAREN_STACK_SIZE];
static limb fibon2[MAX_LEN];
extern limb MontgomeryR1[MAX_LEN];
static int stackIndex, exprIndex;
static int exprLength;
int lang;
char output[3000000];
limb Mult1[MAX_LEN];
limb Mult3[MAX_LEN];
limb Mult4[MAX_LEN];
int q[MAX_LEN];
#define fibon1 MontgomeryR1
static enum eExprErr ComputeSubExpr(void);
// Operators accepted: +, -, *, /, ^, !, F(, L(, P(, B(. N(.
static void SkipSpaces(char *expr);
static int ComputeBack(void);
static int ComputeNext(void);
static enum eExprErr ComputeModInv(void);
static enum eExprErr ComputeFibLucas(int origValue);
static enum eExprErr ComputePartition(void);
static enum eExprErr ComputeExpr(char *expr, BigInteger *ExpressionResult);
static int func(char *expr, BigInteger *ExpressionResult,
  char *funcName, int funcArgs, int leftNumberFlag);
static int type;

enum eExprErr ComputeExpression(char *expr, int typ, BigInteger *ExpressionResult)
{
  int retcode;
  stackIndex = 0;
  exprIndex = 0;
  type = typ;
  retcode = ComputeExpr(expr, ExpressionResult);
  if (retcode != 0) { return retcode; }
  if (ExpressionResult[0].nbrLimbs > 2215 &&    // 10000/log_10(32768)
    ExpressionResult[1].nbrLimbs > 2215)
  {
    return EXPR_NUMBER_TOO_HIGH;
  }
  return 0;
}

static enum eExprErr ComputeExpr(char *expr, BigInteger *ExpressionResult)
{
  int i, j, shLeft;
  int retcode;
  limb carry;
  boolean leftNumberFlag = FALSE;
  int exprIndexAux, offset;
  enum eExprErr SubExprResult;
  int len;
  limb largeLen;
  limb *ptrLimb;
  BigInteger factorial;
  BigInteger *pBigInt;
  int c;
  int startStackIndex = stackIndex;

  exprLength = (int)strlen(expr);
  while (exprIndex < exprLength)
  {
    char charValue;

    charValue = *(expr+exprIndex);
    if (charValue == ' ' || charValue == 9)
    {           // Ignore spaces and horizontal tabs.
      exprIndex++;
      continue;
    }
    if (charValue == '*' && *(expr + exprIndex + 1) == '*')
    {           // Convert double asterisk to exponentiation.
      charValue = '^';
      exprIndex++;
    }
    if (charValue == '!')
    {           // Calculating factorial.
      if (leftNumberFlag == FALSE)
      {
        return EXPR_SYNTAX_ERROR;
      }
      if (stackValues[stackIndex].nbrLimbs > 1)
      {
        return EXPR_INTERM_TOO_HIGH;
      }
      if (stackValues[stackIndex].limbs[0].x < 0 || stackValues[stackIndex].limbs[0].x > 5984)
      {
        return EXPR_INTERM_TOO_HIGH;
      }
      len = (int)stackValues[stackIndex].limbs[0].x;
      factorial.limbs[0].x = 1;
      factorial.nbrLimbs = 1;
      factorial.sign = SIGN_POSITIVE;
      for (i = 2; i <= len; i++)
      {   // Multiply by all integers up to the argument of factorial.
        multint(&factorial, &factorial, i);
      }
      stackValues[stackIndex] = factorial;
      exprIndex++;
      continue;
    }
    if (charValue == '#')
    {           // Calculating primorial.
      if (leftNumberFlag == FALSE)
      {
        return EXPR_SYNTAX_ERROR;
      }
      if (stackValues[stackIndex].nbrLimbs > 2)
      {
        return EXPR_INTERM_TOO_HIGH;
      }
      if (stackValues[stackIndex].nbrLimbs == 2)
      {
        largeLen.x = stackValues[stackIndex].limbs[0].x +
             (stackValues[stackIndex].limbs[1].x << BITS_PER_GROUP);
      }
      else
      {
        largeLen.x = stackValues[stackIndex].limbs[0].x;
      }
      if (largeLen.x < 0 || largeLen.x > 46049)
      {
        return EXPR_INTERM_TOO_HIGH;
      }
      len = (int)largeLen.x;
        // Check if number is prime
      for (i = 2; i*i <= len; i++)
      {
        if (len / i*i == len)
        {   // Number is not prime, so go out.
          return EXPR_INVALID_PARAM;
        }
      }
      factorial.limbs[0].x = 1;
      factorial.nbrLimbs = 1;
      factorial.sign = SIGN_POSITIVE;
      for (i = 2; i <= len; i++)
      {      // Multiply by prime numbers only.
        for (j = 2; j*j <= i; j++)
        {
          if (i / j*j == i)
          {   // Number is not prime.
            break;
          }
        }
        if (j*j > i)
        {     // Number is prime, perform multiplication.
          multint(&factorial, &factorial, i);
        }
      }
      stackValues[stackIndex] = factorial;
      exprIndex++;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "GCD", 2, leftNumberFlag)) <= 0)
    {
      if (retcode != 0) { return retcode; }
      BigIntGcd(&stackValues[stackIndex], &stackValues[stackIndex + 1], &stackValues[stackIndex]);
      leftNumberFlag = 1;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "MODPOW", 3, leftNumberFlag)) <= 0)
    {
      if (retcode != 0) { return retcode; }
      retcode = BigIntGeneralModularPower(&stackValues[stackIndex], &stackValues[stackIndex + 1],
        &stackValues[stackIndex + 2], &stackValues[stackIndex]);
      if (retcode != 0) { return retcode; }
      leftNumberFlag = 1;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "MODINV", 2, leftNumberFlag)) <= 0)
    {
      if (retcode != 0) { return retcode; }
      retcode = ComputeModInv();
      if (retcode != 0) { return retcode; }
      leftNumberFlag = 1;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "F", 1, leftNumberFlag)) <= 0)
    {
      if (retcode != 0) { return retcode; }
      retcode = ComputeFibLucas(0);
      if (retcode != 0) { return retcode; }
      leftNumberFlag = 1;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "L", 1, leftNumberFlag)) <= 0)
    {
      if (retcode != 0) { return retcode; }
      retcode = ComputeFibLucas(2);
      if (retcode != 0) { return retcode; }
      leftNumberFlag = 1;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "P", 1, leftNumberFlag)) <= 0)
    {
      if (retcode != 0) { return retcode; }
      retcode = ComputePartition();
      if (retcode != 0) { return retcode; }
      leftNumberFlag = 1;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "N", 1, leftNumberFlag)) <= 0)
    {
      if (retcode != 0) { return retcode; }
      retcode = ComputeNext();
      if (retcode != 0) { return retcode; }
      leftNumberFlag = 1;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "B", 1, leftNumberFlag)) <= 0)
    {
      if (retcode != 0) { return retcode; }
      retcode = ComputeBack();
      if (retcode != 0) { return retcode; }
      leftNumberFlag = 1;
    }
    else if (charValue == '+' || charValue == '-')
    {
      if (leftNumberFlag == 0)
      {      // Unary plus/minus operator
        exprIndex++;
        if (charValue == '+')
        {
          continue;
        }
        else
        {
          if (stackIndex > startStackIndex && stackOperators[stackIndex - 1] == '_')
          {
            stackIndex--;
            continue;
          }
          if (stackIndex >= PAREN_STACK_SIZE)
          {
            return EXPR_TOO_MANY_PAREN;
          }
          stackOperators[stackIndex++] = '_'; /* Unary minus */
          continue;
        }
      }
      if (stackIndex > startStackIndex &&
        stackOperators[stackIndex - 1] != '(')
      {
        if ((SubExprResult = ComputeSubExpr()) != 0)
        {
          return SubExprResult;
        }
        if (stackIndex > startStackIndex &&
          stackOperators[stackIndex - 1] != '(')
        {
          if ((SubExprResult = ComputeSubExpr()) != 0)
          {
            return SubExprResult;
          }
          if (stackIndex > startStackIndex &&
            stackOperators[stackIndex - 1] != '(')
          {
            if ((SubExprResult = ComputeSubExpr()) != 0)
            {
              return SubExprResult;
            }
          }                         /* end if */
        }                           /* end if */
      }                             /* end if */
      stackOperators[stackIndex++] = charValue;
      leftNumberFlag = 0;
    }                               /* end if */
    else if (charValue == '*' || charValue == '/' || charValue == '%')
    {
      if (leftNumberFlag == 0)
      {
        return EXPR_SYNTAX_ERROR;
      }
      if (stackIndex > startStackIndex && (stackOperators[stackIndex - 1] == '^' ||
        stackOperators[stackIndex - 1] == '*' ||
        stackOperators[stackIndex - 1] == '/'))
      {
        if ((SubExprResult = ComputeSubExpr()) != 0)
        {
          return SubExprResult;
        }
        if (stackIndex > startStackIndex &&
          (stackOperators[stackIndex - 1] == '^' ||
            stackOperators[stackIndex - 1] == '*' ||
            stackOperators[stackIndex - 1] == '/' ||
            stackOperators[stackIndex - 1] == '%'))
        {
          if ((SubExprResult = ComputeSubExpr()) != 0)
          {
            return SubExprResult;
          }
        }                         /* end if */
      }                           /* end if */
      stackOperators[stackIndex++] = charValue;
      leftNumberFlag = 0;
    }
    else if (charValue == '^')
    {
      if (leftNumberFlag == 0)
      {
        return EXPR_SYNTAX_ERROR;
      }
      stackOperators[stackIndex++] = charValue;
      leftNumberFlag = 0;
    }                           /* end if */
    else if (charValue == '(')
    {
      if (leftNumberFlag == 1)
      {
        return EXPR_SYNTAX_ERROR;
      }
      if (stackIndex >= PAREN_STACK_SIZE)
      {
        return EXPR_TOO_MANY_PAREN;
      }
      stackOperators[stackIndex++] = charValue;
    }
    else if (charValue == ')' || charValue == ',')
    {
      if (leftNumberFlag == 0)
      {
        return EXPR_SYNTAX_ERROR;
      }
      if (stackIndex > startStackIndex &&
        stackOperators[stackIndex - 1] != '(')
      {
        if ((SubExprResult = ComputeSubExpr()) != 0)
        {
          return SubExprResult;
        }
        if (stackIndex > startStackIndex &&
          stackOperators[stackIndex - 1] != '(')
        {
          if ((SubExprResult = ComputeSubExpr()) != 0)
          {
            return SubExprResult;
          }
          if (stackIndex > startStackIndex &&
            stackOperators[stackIndex - 1] != '(')
          {
            if ((SubExprResult = ComputeSubExpr()) != 0)
            {
              return SubExprResult;
            }
          }
        }
      }
      if (stackIndex == startStackIndex)
      {
        break;
      }
      if (charValue == ',')
      {
        return EXPR_PAREN_MISMATCH;
      }
      stackIndex--;    /* Discard ')' */
      stackValues[stackIndex] = stackValues[stackIndex + 1];
      leftNumberFlag = 1;
    }
    else if (charValue >= '0' && charValue <= '9')
    {
      exprIndexAux = exprIndex;
      if (charValue == '0' && exprIndexAux < exprLength - 2 &&
          *(expr+exprIndexAux + 1) == 'x')
      {  // hexadecimal
        exprIndexAux += 2;
        while (exprIndexAux < exprLength - 1)
        {
          charValue = *(expr+exprIndexAux + 1);
          if ((charValue >= '0' && charValue <= '9') ||
              (charValue >= 'A' && charValue <= 'F') ||
              (charValue >= 'a' && charValue <= 'f'))
          {
            exprIndexAux++;
          }
          else
          {
            break;
          }
        }
        // Generate big integer from hexadecimal number from right to left.
        carry.x = 0;
        i = 0;  // limb number.
        shLeft = 0;
        offset = exprIndexAux;
        ptrLimb = &stackValues[stackIndex].limbs[0];
        for (; exprIndexAux >= exprIndex + 2; exprIndexAux--)
        {
          c = *(expr + exprIndexAux);
          if (c >= '0' && c <= '9')
          {
            c -= '0';
          }
          else if (c >= 'A' && c <= 'F')
          {
            c -= 'A' - 10;
          }
          else
          {
            c -= 'a' - 10;
          }
          carry.x += c << shLeft;
          shLeft += 4;   // 4 bits per hex digit.
          if (shLeft >= BITS_PER_GROUP)
          {
            shLeft -= BITS_PER_GROUP;
            (ptrLimb++)->x = carry.x & MAX_VALUE_LIMB;
            carry.x = c >> (4-shLeft);
          }
        }
        if (carry.x != 0 || ptrLimb == &stackValues[stackIndex].limbs[0])
        {
          (ptrLimb++)->x = carry.x;
        }
        exprIndex = offset;
        stackValues[stackIndex].nbrLimbs = (int)(ptrLimb - &stackValues[stackIndex].limbs[0]);
        stackValues[stackIndex].sign = SIGN_POSITIVE;
      }
      else
      {                   // Decimal number.
        while (exprIndexAux < exprLength - 1)
        {
          charValue = *(expr+exprIndexAux + 1);
          if (charValue >= '0' && charValue <= '9')
          {
            exprIndexAux++;
          }
          else
          {
            break;
          }
        }
        // Generate big integer from decimal number
        pBigInt = &stackValues[stackIndex];
        Dec2Bin(expr + exprIndex, pBigInt->limbs,
                exprIndexAux + 1 - exprIndex, &pBigInt -> nbrLimbs);
        pBigInt -> sign = SIGN_POSITIVE;
        exprIndex = exprIndexAux;
      }
      leftNumberFlag = TRUE;
    }
    else
    {
      return EXPR_SYNTAX_ERROR;
    }
    exprIndex++;
  }                              /* end while */
  if (leftNumberFlag == FALSE)
  {
    return EXPR_SYNTAX_ERROR;
  }
  if (stackIndex > startStackIndex && stackOperators[stackIndex - 1] != '(')
  {
    if ((SubExprResult = ComputeSubExpr()) != 0)
    {
      return SubExprResult;
    }
    if (stackIndex > startStackIndex && stackOperators[stackIndex - 1] != '(')
    {
      if ((SubExprResult = ComputeSubExpr()) != 0)
      {
        return SubExprResult;
      }
      if (stackIndex > startStackIndex && stackOperators[stackIndex - 1] != '(')
      {
        if ((SubExprResult = ComputeSubExpr()) != 0)
        {
          return SubExprResult;
        }
      }
    }
  }
  if (stackIndex != startStackIndex)
  {
    return EXPR_PAREN_MISMATCH;
  }
  CopyBigInt(ExpressionResult, &stackValues[0]);
  return EXPR_OK;
}

static void SkipSpaces(char *expr)
{
  while (*(expr + exprIndex))
  {
    if (*(expr + exprIndex) > ' ')
    {
      break;
    }
    exprIndex++;
  }
  return;
}

static enum eExprErr ComputeSubExpr(void)
{
  char stackOper = stackOperators[--stackIndex];
  switch (stackOper)
  {
  case '+':
    BigIntAdd(&stackValues[stackIndex], &stackValues[stackIndex + 1], &stackValues[stackIndex]);
    return EXPR_OK;
  case '-':
    BigIntSubt(&stackValues[stackIndex], &stackValues[stackIndex + 1], &stackValues[stackIndex]);
    return EXPR_OK;
  case '_':
    BigIntNegate(&stackValues[stackIndex + 1], &stackValues[stackIndex]);
    return EXPR_OK;
  case '/':
    return BigIntDivide(&stackValues[stackIndex], &stackValues[stackIndex + 1],
      &stackValues[stackIndex]);
  case '*':
    return BigIntMultiply(&stackValues[stackIndex], &stackValues[stackIndex + 1],
      &stackValues[stackIndex]);
  case '%':
    return BigIntRemainder(&stackValues[stackIndex], &stackValues[stackIndex + 1],
      &stackValues[stackIndex]);
  case '^':
    return BigIntPower(&stackValues[stackIndex], &stackValues[stackIndex + 1],
      &stackValues[stackIndex]);
  }
  return EXPR_OK;
}

static int func(char *expr, BigInteger *ExpressionResult,
  char *funcName, int funcArgs, int leftNumberFlag)
{
  int index;
  int funcNameLen = (int)strlen(funcName);
  char *ptrExpr, *ptrFuncName;

  if (exprIndex + funcNameLen > exprLength)
  {
    return 1;
  }
  ptrExpr = expr + exprIndex;
  ptrFuncName = funcName;
  while (*ptrFuncName)
  {
    if ((*ptrExpr & 0xDF) != *ptrFuncName)
    {
      return 1;
    }
    ptrExpr++;
    ptrFuncName++;
  }
  exprIndex += funcNameLen;
  if (leftNumberFlag == 1)
  {
    return EXPR_SYNTAX_ERROR;
  }
  SkipSpaces(expr);
  if (exprIndex == exprLength || *(expr + exprIndex++) != '(')
  {
    return EXPR_SYNTAX_ERROR;
  }
  for (index = 0; index < funcArgs; index++)
  {
    int retcode;
    char compareChar;

    SkipSpaces(expr);
    if (stackIndex >= PAREN_STACK_SIZE)
    {
      return EXPR_TOO_MANY_PAREN;
    }
    retcode = ComputeExpr(expr, ExpressionResult);
    if (retcode != 0) { return retcode; }
    SkipSpaces(expr);
    compareChar = (index == funcArgs - 1 ? ')' : ',');
    if (exprIndex == exprLength || *(expr + exprIndex++) != compareChar)
    {
      return EXPR_SYNTAX_ERROR;
    }
    stackIndex++;
  }
  stackIndex -= funcArgs;
  exprIndex--;
  return 0;
}

// Compute previous probable prime
static int ComputeBack(void)
{
  BigInteger *pArgument = &stackValues[stackIndex];;
  BigInteger *pResult = &stackValues[stackIndex];
  limb *pResultLimbs = pResult->limbs;
  limb *pArgumentLimbs = pArgument->limbs;
  int nbrLimbs = pArgument->nbrLimbs;
  pResult->sign = SIGN_POSITIVE;
  if (pArgument->sign == SIGN_NEGATIVE || (nbrLimbs == 1 && pArgumentLimbs->x < 3))
  {
    return EXPR_INVALID_PARAM;
  }
  if (nbrLimbs == 1 && pArgumentLimbs->x == 3)
  {
    pResult->nbrLimbs = 1;
    pResultLimbs->x = 2;
    return EXPR_OK;
  }
  memcpy(pResultLimbs, pArgumentLimbs, nbrLimbs * sizeof(limb));
  pResultLimbs->x |= 1;  // If number is even, use next odd number.
  pResult->nbrLimbs = nbrLimbs;
  do
  {        // Loop that searches for previous or next probable prime.
    int ctr;
    pResultLimbs->x -= 2;
    if (pResultLimbs->x < 0)
    {
      limb *pTemp = pResultLimbs;
      for (ctr = 1; ctr < nbrLimbs; ctr++)
      {
        (pTemp++)->x = MAX_VALUE_LIMB;
        if (--(pTemp->x) >= 0)
        {
          break;
        }
      }
      if (pTemp->x == 0)
      {  // Most significant limb is zero. Decrement number of limbs.
        nbrLimbs--;
        pResult->nbrLimbs = nbrLimbs;
      }
    }
  } while (BpswPrimalityTest(pResult));  // Continue loop if not probable prime.
  return EXPR_OK;
}

// Compute next probable prime
static int ComputeNext(void)
{
  BigInteger *pArgument = &stackValues[stackIndex];;
  BigInteger *pResult = &stackValues[stackIndex];
  limb *pResultLimbs = pResult->limbs;
  limb *pArgumentLimbs = pArgument->limbs;
  int nbrLimbs = pArgument->nbrLimbs;
  int skipUpdate = 0;    // Do not skip first update in advance.
  pResult->sign = SIGN_POSITIVE;
  if (pArgument->sign == SIGN_NEGATIVE || (nbrLimbs == 1 && pArgumentLimbs->x < 2))
  {
    pResult->nbrLimbs = 1;
    pResultLimbs->x = 2;
    return EXPR_OK;
  }

  memcpy(pResultLimbs, pArgumentLimbs, nbrLimbs * sizeof(limb));
  if ((pResultLimbs->x & 1) == 0)
  {   // Number is even.
    pResultLimbs->x++;
    skipUpdate = 1;   // Skip first update.
  }
  pResult->nbrLimbs = nbrLimbs;
  do
  {        // Loop that searches for previous or next probable prime.
    if (skipUpdate == 0)
    {
      if (pResultLimbs->x < MAX_VALUE_LIMB)
      {                       // No overflow.
        pResultLimbs->x += 2; // Add 2.
      }
      else
      {                       // Overflow.
        int ctr;
        limb *pTemp = pResultLimbs;
        (pTemp++)->x = 1;
        for (ctr = 1; ctr < nbrLimbs; ctr++)
        {
          if (pTemp->x < MAX_VALUE_LIMB)
          {
            pTemp->x++;
            break;
          }
          (pTemp++)->x = 0;
        }
        if (ctr == nbrLimbs)
        {
          pTemp->x = 1;
          nbrLimbs++;
          pResult->nbrLimbs = nbrLimbs;
        }
      }
    }
    else
    {
      skipUpdate = 0;
    }
  } while (BpswPrimalityTest(pResult));  // Continue loop if not probable prime.
  return EXPR_OK;
}

static enum eExprErr ComputeModInv(void)
{
  BigInteger one;
  BigInteger *pDiv = &stackValues[stackIndex + 1];
  if (pDiv->nbrLimbs == 1 && pDiv->limbs[0].x == 0)
  {
    return EXPR_DIVIDE_BY_ZERO;
  }
  // Check that the arguments are relatively prime.
  BigIntGcd(&stackValues[stackIndex], &stackValues[stackIndex + 1], &one);
  if (one.nbrLimbs != 1 || one.limbs[0].x != 1)
  {
    return EXPR_ARGUMENTS_NOT_RELATIVELY_PRIME;
  }
  intToBigInteger(&one, 1);
  BigIntGeneralModularDivision(&one, &stackValues[stackIndex], pDiv, &stackValues[stackIndex]);
  return EXPR_OK;
}

static enum eExprErr ComputeFibLucas(int origValue)
{
  BigInteger *pArgument = &stackValues[stackIndex];
  limb largeVal;
  int val, len;
  limb *pFibonPrev, *pFibonAct;
  if (pArgument->sign == SIGN_NEGATIVE)
  {
    return EXPR_INVALID_PARAM;
  }
  if (pArgument->nbrLimbs > 1)
  {
    return EXPR_INTERM_TOO_HIGH;
  }
  largeVal.x = pArgument->limbs[0].x;
  if (largeVal.x > 95662)
  {
    return EXPR_INTERM_TOO_HIGH;
  }
  val = (int)largeVal.x;
  pFibonPrev = fibon1;
  pFibonAct = fibon2;
  len = 1;
  if (val == 0)
  {        // F(0) = 0, L(0) = 2
    fibon2[0].x = origValue;
  }
  else
  {
    int i, j;
    // For Lucas sequences: FibonPrev = 2, FibonAct = 1
    // For Fibonacci sequences: FibonPrev = 0, FibonAct = 1
    fibon1[0].x = origValue;
    fibon2[0].x = 1;
    for (i = 1; i < val; i++)
    {
      limb *pTemp;
      unsigned int carry = 0;
      for (j = 0; j < len; j++)
      {
        carry += (pFibonPrev + j)->x + (pFibonAct + j)->x;
        (pFibonPrev + j)->x = carry & MAX_VALUE_LIMB;
        carry >>= BITS_PER_GROUP;
      }
      if (carry != 0)
      {
        (pFibonPrev + j)->x = carry;
        (pFibonAct + j)->x = 0;
        len++;
      }
      pTemp = pFibonAct;
      pFibonAct = pFibonPrev;
      pFibonPrev = pTemp;
    }
  }
  pArgument = &stackValues[stackIndex];
  pArgument->sign = SIGN_POSITIVE;
  pArgument->nbrLimbs = len;
  memcpy(pArgument->limbs, pFibonAct, len * sizeof(limb));
  return EXPR_OK;
}

static enum eExprErr ComputePartition(void)
{
  BigInteger *pArgument = &stackValues[stackIndex];
  limb largeVal;
  int val;

  if (pArgument->sign == SIGN_NEGATIVE)
  {
    return EXPR_INVALID_PARAM;
  }
  if (pArgument->nbrLimbs > 2)
  {
    return EXPR_INTERM_TOO_HIGH;
  }
  if (pArgument->nbrLimbs == 2)
  {
    largeVal.x = pArgument->limbs[0].x + (pArgument->limbs[1].x << BITS_PER_GROUP);
  }
  else
  {
    largeVal.x = pArgument->limbs[0].x;
  }
  if (largeVal.x > 3520000)
  {
    return EXPR_INTERM_TOO_HIGH;
  }
  val = (int)largeVal.x;
  if (val > 100000)
  {
    return EXPR_INVALID_PARAM;
  }
  partition(val, &stackValues[stackIndex]);
  return EXPR_OK;
}

#ifndef __EMSCRIPTEN__
void databack(char *data)
{
  (void *)data;
}
#endif
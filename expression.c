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
#include "factor.h"
#define PAREN_STACK_SIZE            100
#define OPER_POWER                    1
#define OPER_MULTIPLY                 2
#define OPER_DIVIDE                   3
#define OPER_REMAINDER                4
#define OPER_UNARY_MINUS              5
#define OPER_PLUS                     6
#define OPER_MINUS                    7
#define OPER_SHR                      8
#define OPER_SHL                      9
#define OPER_NOT_GREATER             10
#define OPER_NOT_LESS                11
#define OPER_NOT_EQUAL               12
#define OPER_EQUAL                   13
#define OPER_GREATER                 14
#define OPER_LESS                    15
#define OPER_NOT                     16
#define OPER_AND                     17
#define OPER_OR                      18
#define OPER_XOR                     19
#define MAXIMUM_OPERATOR             19

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
BigInteger valueX;
int counterC;
#define fibon1 MontgomeryR1
static enum eExprErr ComputeSubExpr(void);
static void SkipSpaces(char *expr);
static int ComputeBack(void);
static int ComputeNext(void);
#ifdef FACTORIZATION_FUNCTIONS
static int ComputeTotient(void);
static int ComputeNumDivs(void);
static int ComputeSumDivs(void);
static int ComputeConcatFact(void);
static char textFactor[30000];
#endif
static int ComputeSumDigits(void);
static int ComputeRevDigits(void);
static int ComputeNumDigits(void);
static enum eExprErr ComputeModInv(void);
static enum eExprErr ComputeFibLucas(int origValue);
static enum eExprErr ComputePartition(void);
static enum eExprErr ComputeExpr(char *expr, BigInteger *ExpressionResult);
static enum eExprErr ShiftLeft(BigInteger* first, BigInteger *second, BigInteger *result);
static int func(char *expr, BigInteger *ExpressionResult,
  char *funcName, int funcArgs, int leftNumberFlag);
static int type;
static int valueXused;

static char priority[] = 
{
  1,                // Power
  2, 2, 2,          // Multiply, divide and remainder.
  3,                // Unary minus.
  4, 4,             // Plus and minus.
  5, 5,             // Shift right and left.
  6, 6, 6, 6, 6, 6, // Six comparison operators (equal, greater, less, etc.)
  7,                // NOT.
  8, 8, 8,          // AND, OR, XOR.
};

enum eExprErr ComputeExpression(char *expr, int typ, BigInteger *ExpressionResult)
{
  int retcode;
  valueXused = FALSE;
  stackIndex = 0;
  exprIndex = 0;
  type = typ;
  retcode = ComputeExpr(expr, ExpressionResult);
  if (retcode != 0) { return retcode; }
  if (ExpressionResult[0].nbrLimbs > 33219 / BITS_PER_GROUP + 1)    // 10000/log_10(2) = 33219
  {
    return EXPR_NUMBER_TOO_HIGH;
  }
  if (valueX.nbrLimbs > 0 && valueXused == FALSE)
  {
    return EXPR_VAR_OR_COUNTER_REQUIRED;
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
    if (charValue == '^')
    {           // Caret is exponentiation operation.
      charValue = OPER_POWER;
      exprIndex++;
    }
    else if (charValue == '*' && *(expr + exprIndex + 1) == '*')
    {           // Double asterisk is exponentiation operation too.
      charValue = OPER_POWER;
      exprIndex += 2;
    }
    else if (charValue == '*')
    {
      charValue = OPER_MULTIPLY;
      exprIndex++;
    }
    else if (charValue == '/')
    {
      charValue = OPER_DIVIDE;
      exprIndex++;
    }
    else if (charValue == '%')
    {
      charValue = OPER_REMAINDER;
      exprIndex++;
    }
    else if (charValue == '+')
    {
      charValue = OPER_PLUS;
      exprIndex++;
    }
    else if (charValue == '-')
    {
      charValue = OPER_MINUS;
      exprIndex++;
    }
    else if (charValue == '<' && *(expr + exprIndex + 1) == '=')
    {
      charValue = OPER_NOT_GREATER;
      exprIndex += 2;
    }
    else if (charValue == '>' && *(expr + exprIndex + 1) == '=')
    {
      charValue = OPER_NOT_LESS;
      exprIndex += 2;
    }
    else if (charValue == '!' && *(expr + exprIndex + 1) == '=')
    {
      charValue = OPER_NOT_EQUAL;
      exprIndex += 2;
    }
    else if (charValue == '=' && *(expr + exprIndex + 1) == '=')
    {
      charValue = OPER_EQUAL;
      exprIndex += 2;
    }
    else if (charValue == '>')
    {
      charValue = OPER_GREATER;
      exprIndex++;
    }
    else if (charValue == '<')
    {
      charValue = OPER_LESS;
      exprIndex++;
    }
    else if ((charValue & 0xDF) == 'N' && (*(expr + exprIndex + 1) & 0xDF) == 'O' &&
      (*(expr + exprIndex + 2) & 0xDF) == 'T')
    {
      charValue = OPER_NOT;
      exprIndex += 3;
    }
    else if ((charValue & 0xDF) == 'A' && (*(expr + exprIndex + 1) & 0xDF) == 'N' &&
      (*(expr + exprIndex + 2) & 0xDF) == 'D')
    {
      charValue = OPER_AND;
      exprIndex += 3;
    }
    else if ((charValue & 0xDF) == 'O' && (*(expr + exprIndex + 1) & 0xDF) == 'R')
    {
      charValue = OPER_OR;
      exprIndex += 2;
    }
    else if ((charValue & 0xDF) == 'X' && (*(expr + exprIndex + 1) & 0xDF) == 'O' &&
      (*(expr + exprIndex + 2) & 0xDF) == 'R')
    {
      charValue = OPER_XOR;
      exprIndex += 3;
    }
    else if ((charValue & 0xDF) == 'S' && (*(expr + exprIndex + 1) & 0xDF) == 'H' &&
      (*(expr + exprIndex + 2) & 0xDF) == 'L')
    {
      charValue = OPER_SHL;
      exprIndex += 3;
    }
    else if ((charValue & 0xDF) == 'S' && (*(expr + exprIndex + 1) & 0xDF) == 'H' &&
      (*(expr + exprIndex + 2) & 0xDF) == 'R')
    {
      charValue = OPER_SHR;
      exprIndex += 3;
    }
    else if (charValue == '!')
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
    else if (charValue == '#')
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
      len = (int)largeLen.x;
      if (len < 0 || len > 46049)
      {
        return EXPR_INTERM_TOO_HIGH;
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
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "MODPOW", 3, leftNumberFlag)) <= 0)
    {
      if (retcode != 0) { return retcode; }
      retcode = BigIntGeneralModularPower(&stackValues[stackIndex], &stackValues[stackIndex + 1],
        &stackValues[stackIndex + 2], &stackValues[stackIndex]);
      if (retcode != 0) { return retcode; }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "MODINV", 2, leftNumberFlag)) <= 0)
    {
      if (retcode != 0) { return retcode; }
      retcode = ComputeModInv();
      if (retcode != 0) { return retcode; }
      leftNumberFlag = 1;
      continue;
    }
#ifdef FACTORIZATION_FUNCTIONS
    else if ((retcode = func(expr, ExpressionResult,
      "TOTIENT", 1, leftNumberFlag)) <= 0)
    {
      if (retcode != 0) { return retcode; }
      retcode = ComputeTotient();
      if (retcode != 0) { return retcode; }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "NUMDIVS", 1, leftNumberFlag)) <= 0)
    {
      if (retcode != 0) { return retcode; }
      retcode = ComputeNumDivs();
      if (retcode != 0) { return retcode; }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "SUMDIVS", 1, leftNumberFlag)) <= 0)
    {
      if (retcode != 0) { return retcode; }
      retcode = ComputeSumDivs();
      if (retcode != 0) { return retcode; }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "CONCATFACT", 2, leftNumberFlag)) <= 0)
    {
      if (retcode != 0) { return retcode; }
      retcode = ComputeConcatFact();
      if (retcode != 0) { return retcode; }
      leftNumberFlag = 1;
      continue;
    }
#endif
    else if ((retcode = func(expr, ExpressionResult,
      "SUMDIGITS", 2, leftNumberFlag)) <= 0)
    {
      if (retcode != 0) { return retcode; }
      retcode = ComputeSumDigits();
      if (retcode != 0) { return retcode; }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "NUMDIGITS", 2, leftNumberFlag)) <= 0)
    {
      if (retcode != 0) { return retcode; }
      retcode = ComputeNumDigits();
      if (retcode != 0) { return retcode; }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "REVDIGITS", 2, leftNumberFlag)) <= 0)
    {
      if (retcode != 0) { return retcode; }
      retcode = ComputeRevDigits();
      if (retcode != 0) { return retcode; }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "ISPRIME", 1, leftNumberFlag)) <= 0)
    {
      if (retcode != 0) { return retcode; }
      if (BpswPrimalityTest(&stackValues[stackIndex]) == 0)
      {    // Argument is a probable prime.
        intToBigInteger(&stackValues[stackIndex], -1);
      }
      else
      {    // Argument is not a probable prime.
        intToBigInteger(&stackValues[stackIndex], 0);
      }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "F", 1, leftNumberFlag)) <= 0)
    {
      if (retcode != 0) { return retcode; }
      retcode = ComputeFibLucas(0);
      if (retcode != 0) { return retcode; }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "L", 1, leftNumberFlag)) <= 0)
    {
      if (retcode != 0) { return retcode; }
      retcode = ComputeFibLucas(2);
      if (retcode != 0) { return retcode; }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "P", 1, leftNumberFlag)) <= 0)
    {
      if (retcode != 0) { return retcode; }
      retcode = ComputePartition();
      if (retcode != 0) { return retcode; }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "N", 1, leftNumberFlag)) <= 0)
    {
      if (retcode != 0) { return retcode; }
      retcode = ComputeNext();
      if (retcode != 0) { return retcode; }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "B", 1, leftNumberFlag)) <= 0)
    {
      if (retcode != 0) { return retcode; }
      retcode = ComputeBack();
      if (retcode != 0) { return retcode; }
      leftNumberFlag = 1;
      continue;
    }
    else if ((charValue & 0xDF) == 'X')
    {
      if (leftNumberFlag || valueX.nbrLimbs == 0)
      {
        return EXPR_SYNTAX_ERROR;
      }
      CopyBigInt(&stackValues[stackIndex], &valueX);
      valueXused = TRUE;
      exprIndex++;
      leftNumberFlag = TRUE;
      continue;
    }
    else if ((charValue & 0xDF) == 'C')
    {
      if (leftNumberFlag || valueX.nbrLimbs == 0)
      {
        return EXPR_SYNTAX_ERROR;
      }
      intToBigInteger(&stackValues[stackIndex], counterC);
      valueXused = TRUE;
      exprIndex++;
      leftNumberFlag = TRUE;
      continue;
    }
    else if (charValue == '(')
    {
      if (leftNumberFlag == TRUE)
      {
        return EXPR_SYNTAX_ERROR;
      }
      if (stackIndex >= PAREN_STACK_SIZE)
      {
        return EXPR_TOO_MANY_PAREN;
      }
      stackOperators[stackIndex++] = charValue;
      exprIndex++;
      continue;
    }
    else if (charValue == ')' || charValue == ',')
    {
      if (leftNumberFlag == 0)
      {
        return EXPR_SYNTAX_ERROR;
      }
      while (stackIndex > startStackIndex &&
        stackOperators[stackIndex - 1] != '(')
      {
        if ((SubExprResult = ComputeSubExpr()) != 0)
        {
          return SubExprResult;
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
      exprIndex++;
      continue;
    }
    else if (charValue >= '0' && charValue <= '9')
    {
      exprIndexAux = exprIndex;
      if (charValue == '0' && exprIndexAux < exprLength - 2 &&
          *(expr+exprIndexAux + 1) == 'x')
      {  // hexadecimal
        exprIndexAux++;
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
        exprIndex = offset+1;
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
        exprIndex = exprIndexAux + 1;
      }
      leftNumberFlag = TRUE;
      continue;
    }
    if (charValue <= MAXIMUM_OPERATOR)
    {
      if ((charValue == OPER_PLUS || charValue == OPER_MINUS) && leftNumberFlag == 0)
      {                    // Unary plus/minus operator
        if (charValue == OPER_PLUS)
        {
          continue;
        }
        else
        {
          if (stackIndex > startStackIndex && stackOperators[stackIndex - 1] == OPER_UNARY_MINUS)
          {
            stackIndex--;
            continue;
          }
          if (stackIndex >= PAREN_STACK_SIZE)
          {
            return EXPR_TOO_MANY_PAREN;
          }
          stackOperators[stackIndex++] = OPER_UNARY_MINUS; /* Unary minus */
          continue;
        }
      }
      if ((leftNumberFlag == 0) != (charValue == OPER_NOT))
      {     // Missing left operator if operator is not NOT or
            // extra left operator if operator is NOT.
        return EXPR_SYNTAX_ERROR;
      }
      if (charValue != OPER_POWER)
      {  // Power operator has right associativity.
        while (stackIndex > startStackIndex &&
          stackOperators[stackIndex - 1] != '(' &&
          priority[(int)stackOperators[stackIndex - 1] - 1] <= priority[(int)charValue] - 1)
        {
          if ((SubExprResult = ComputeSubExpr()) != 0)
          {
            return SubExprResult;
          }
        }
      }
      stackOperators[stackIndex++] = charValue;
      leftNumberFlag = 0;
      continue;
    }
    return EXPR_SYNTAX_ERROR;
  }                              /* end while */
  if (leftNumberFlag == FALSE)
  {
    return EXPR_SYNTAX_ERROR;
  }
  while (stackIndex > startStackIndex && stackOperators[stackIndex - 1] != '(')
  {
    if ((SubExprResult = ComputeSubExpr()) != 0)
    {
      return SubExprResult;
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

static void ConvertToTwosComplement(BigInteger *value)
{
  int idx;
  int nbrLimbs;
  limb *ptrLimb;
  if (value->sign == SIGN_POSITIVE)
  {    // If number is positive, no conversion is needed.
    return;
  }
  nbrLimbs = value->nbrLimbs;
  ptrLimb = &value->limbs[0];
  for (idx = 0; idx < nbrLimbs; idx++)
  {
    if (ptrLimb->x != 0)
    {
      break;
    }
    ptrLimb++;
  }
  if (idx < nbrLimbs)
  {
    ptrLimb->x = 0x80000000 - ptrLimb->x;
    ptrLimb++;
  }
  for (; idx < nbrLimbs; idx++)
  {
    ptrLimb->x = 0x7FFFFFFF - ptrLimb->x;
    ptrLimb++;
  }
}

static enum eExprErr ComputeSubExpr(void)
{
  char stackOper = stackOperators[--stackIndex];
  BigInteger *firstArg = &stackValues[stackIndex];
  BigInteger *secondArg = &stackValues[stackIndex + 1];
  BigInteger *result = &stackValues[stackIndex];
  BigInteger *tmpptr;
  int idx;
  switch (stackOper)
  {
  case OPER_PLUS:
    BigIntAdd(firstArg, secondArg, result);
    return EXPR_OK;
  case OPER_MINUS:
    BigIntSubt(firstArg, secondArg, result);
    return EXPR_OK;
  case OPER_UNARY_MINUS:
    BigIntNegate(secondArg, result);
    return EXPR_OK;
  case OPER_DIVIDE:
    return BigIntDivide(firstArg, secondArg, result);
  case OPER_MULTIPLY:
    return BigIntMultiply(firstArg, secondArg, result);
  case OPER_REMAINDER:
    return BigIntRemainder(firstArg, secondArg, result);
  case OPER_POWER:
    return BigIntPower(firstArg, secondArg, result);
  case OPER_EQUAL:
    BigIntSubt(firstArg, secondArg, result);
    intToBigInteger(result, (result->nbrLimbs == 1 && result->limbs[0].x == 0? -1: 0));
    return EXPR_OK;
  case OPER_NOT_EQUAL:
    BigIntSubt(firstArg, secondArg, result);
    intToBigInteger(result, (result->nbrLimbs == 1 && result->limbs[0].x == 0 ? 0 : -1));
    return EXPR_OK;
  case OPER_GREATER:
    BigIntSubt(secondArg, firstArg, result);
    intToBigInteger(result, result->sign == SIGN_NEGATIVE ? -1 : 0);
    return EXPR_OK;
  case OPER_NOT_GREATER:
    BigIntSubt(secondArg, firstArg, result);
    intToBigInteger(result, result->sign == SIGN_NEGATIVE ? 0 : -1);
    return EXPR_OK;
  case OPER_LESS:
    BigIntSubt(firstArg, secondArg, result);
    intToBigInteger(result, result->sign == SIGN_NEGATIVE ? -1 : 0);
    return EXPR_OK;
  case OPER_NOT_LESS:
    BigIntSubt(firstArg, secondArg, result);
    intToBigInteger(result, result->sign == SIGN_NEGATIVE ? 0 : -1);
    return EXPR_OK;
  case OPER_SHL:
    ShiftLeft(firstArg, secondArg, result);
    return EXPR_OK;
  case OPER_SHR:
    if (secondArg->sign == SIGN_POSITIVE)
    {
      secondArg->sign = SIGN_NEGATIVE;
    }
    else
    {
      secondArg->sign = SIGN_POSITIVE;
    }
    ShiftLeft(firstArg, secondArg, result);
    return EXPR_OK;
  case OPER_NOT:    // Perform binary NOT as result <- -1 - argument.
    intToBigInteger(firstArg, -1);
    BigIntSubt(firstArg, secondArg, result);
    return EXPR_OK;
  case OPER_AND:    // Perform binary AND.
    if (firstArg->nbrLimbs < secondArg->nbrLimbs)
    {    // After the exchange, firstArg has not fewer limbs than secondArg.
      tmpptr = firstArg;
      firstArg = secondArg;
      secondArg = tmpptr;
    }
    ConvertToTwosComplement(firstArg);
    ConvertToTwosComplement(secondArg);
    for (idx = 0; idx < secondArg->nbrLimbs; idx++)
    {
      result->limbs[idx].x = firstArg->limbs[idx].x & secondArg->limbs[idx].x;
    }
    if (secondArg->sign == SIGN_POSITIVE)
    {
      result->nbrLimbs = secondArg->nbrLimbs;
    }
    else
    {
      result->nbrLimbs = firstArg->nbrLimbs;
      for (; idx < firstArg->nbrLimbs; idx++)
      {
        result->limbs[idx].x = firstArg->limbs[idx].x;
      }
    }
    if (firstArg->sign == SIGN_POSITIVE || secondArg->sign == SIGN_POSITIVE)
    {
      result->sign = SIGN_POSITIVE;
    }
    else
    {
      result->sign = SIGN_NEGATIVE;
    }
    ConvertToTwosComplement(result);
    return EXPR_OK;
  case OPER_OR:    // Perform binary OR.
    if (firstArg->nbrLimbs < secondArg->nbrLimbs)
    {    // After the exchange, firstArg has not fewer limbs than secondArg.
      tmpptr = firstArg;
      firstArg = secondArg;
      secondArg = tmpptr;
    }
    ConvertToTwosComplement(firstArg);
    ConvertToTwosComplement(secondArg);
    for (idx = 0; idx < secondArg->nbrLimbs; idx++)
    {
      result->limbs[idx].x = firstArg->limbs[idx].x | secondArg->limbs[idx].x;
    }
    if (secondArg->sign == SIGN_NEGATIVE)
    {
      result->nbrLimbs = secondArg->nbrLimbs;
    }
    else
    {
      result->nbrLimbs = firstArg->nbrLimbs;
      for (; idx < firstArg->nbrLimbs; idx++)
      {
        result->limbs[idx].x = firstArg->limbs[idx].x;
      }
    }
    if (firstArg->sign == SIGN_NEGATIVE || secondArg->sign == SIGN_NEGATIVE)
    {
      result->sign = SIGN_NEGATIVE;
    }
    else
    {
      result->sign = SIGN_POSITIVE;
    }
    ConvertToTwosComplement(result);
    return EXPR_OK;
  case OPER_XOR:    // Perform binary XOR.
    if (firstArg->nbrLimbs < secondArg->nbrLimbs)
    {    // After the exchange, firstArg has not fewer limbs than secondArg.
      tmpptr = firstArg;
      firstArg = secondArg;
      secondArg = tmpptr;
    }
    ConvertToTwosComplement(firstArg);
    ConvertToTwosComplement(secondArg);
    for (idx = 0; idx < secondArg->nbrLimbs; idx++)
    {
      result->limbs[idx].x = firstArg->limbs[idx].x ^ secondArg->limbs[idx].x;
    }
    if (secondArg->sign == SIGN_POSITIVE)
    {
      for (; idx < firstArg->nbrLimbs; idx++)
      {
        result->limbs[idx].x = firstArg->limbs[idx].x;
      }
    }
    else
    {
      for (; idx < firstArg->nbrLimbs; idx++)
      {
        result->limbs[idx].x = firstArg->limbs[idx].x ^ MAX_INT_NBR;
      }
    }
    if ((firstArg->sign == SIGN_NEGATIVE) != (secondArg->sign == SIGN_NEGATIVE))
    {
      result->sign = SIGN_NEGATIVE;
    }
    else
    {
      result->sign = SIGN_POSITIVE;
    }
    result->nbrLimbs = firstArg->nbrLimbs;
    ConvertToTwosComplement(result);
    return EXPR_OK;
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

static void PerformFactorization(BigInteger *tofactor)
{
  NumberLength = tofactor->nbrLimbs;
  CompressBigInteger(nbrToFactor, tofactor);
  if (hexadecimal)
  {
    Bin2Hex(tofactor->limbs, tofactorDec, tofactor->nbrLimbs, groupLen);
  }
  else
  {
    Bin2Dec(tofactor->limbs, tofactorDec, tofactor->nbrLimbs, groupLen);
  }
  factor(tofactor, nbrToFactor, factorsMod, astFactorsMod, NULL);
}

#ifdef FACTORIZATION_FUNCTIONS
static int ComputeTotient(void)
{
  PerformFactorization(&stackValues[stackIndex]);
  Totient(&stackValues[stackIndex]);
  return EXPR_OK;
}

static int ComputeNumDivs(void)
{
  PerformFactorization(&stackValues[stackIndex]);
  NumberOfDivisors(&stackValues[stackIndex]);
  return EXPR_OK;
}

static int ComputeSumDivs(void)
{
  PerformFactorization(&stackValues[stackIndex]);
  SumOfDivisors(&stackValues[stackIndex]);
  return EXPR_OK;
}

static int ComputeConcatFact(void)
{
  BigInteger *mode = &stackValues[stackIndex];
  BigInteger factorValue;
  struct sFactors *pstFactor;
  int factorNumber, nbrFactors;
  int descend = mode->limbs[0].x & 1;
  int repeated = mode->limbs[0].x & 2;
  char *ptrTextFactor = textFactor;
  if (mode->nbrLimbs > 1 || mode->sign == SIGN_NEGATIVE || mode->limbs[0].x > 3)
  {      // The valid modes are 0, 1, 2 and 3.
    return EXPR_INVALID_PARAM;
  }
  PerformFactorization(&stackValues[stackIndex + 1]); // Factor second argument.
  nbrFactors = astFactorsMod[0].multiplicity;
  for (factorNumber = 1; factorNumber <= nbrFactors; factorNumber++)
  {
    int ctr;
    pstFactor = &astFactorsMod[descend ? nbrFactors - factorNumber + 1 : factorNumber];
    UncompressBigInteger(pstFactor->ptrFactor, &factorValue);
    ctr = (repeated ? pstFactor->multiplicity : 1);
    for (; ctr > 0; ctr--)
    {
      BigInteger2Dec(&factorValue, ptrTextFactor, 0);
      ptrTextFactor += strlen(ptrTextFactor);
    }
  }
  if (ptrTextFactor == &textFactor[0])
  {
    intToBigInteger(&stackValues[stackIndex], 0);
  }
  else
  {
    Dec2Bin(textFactor, stackValues[stackIndex].limbs, (int)(ptrTextFactor - &textFactor[0]),
      &stackValues[stackIndex].nbrLimbs);
  }
  return EXPR_OK;
}
#endif

static int ComputeSumDigits(void)
{
  BigInteger argum, Temp;
  BigInteger *result = &stackValues[stackIndex];
  BigInteger *radix = &stackValues[stackIndex + 1];
  CopyBigInt(&argum, &stackValues[stackIndex]);
  intToBigInteger(result, 0);
  while (argum.nbrLimbs > 1 || argum.limbs[0].x > 0)
  {
    BigIntRemainder(&argum, radix, &Temp);
    BigIntAdd(result, &Temp, result);
    BigIntDivide(&argum, radix, &argum);
  }
  return EXPR_OK;
}

static int ComputeNumDigits(void)
{
  BigInteger *result = &stackValues[stackIndex];
  BigInteger *radix = &stackValues[stackIndex + 1];
  int digits = 0;
  while (result->nbrLimbs > 1 || result->limbs[0].x > 0)
  {
    BigIntDivide(result, radix, result);
    digits++;
  }
  intToBigInteger(result, digits);
  return EXPR_OK;
}

static int ComputeRevDigits(void)
{
  BigInteger argum, Temp;
  BigInteger *result = &stackValues[stackIndex];
  BigInteger *radix = &stackValues[stackIndex+1];
  CopyBigInt(&argum, &stackValues[stackIndex]);
  intToBigInteger(result, 0);
  while (argum.nbrLimbs > 1 || argum.limbs[0].x > 0)
  {
    BigIntRemainder(&argum, radix, &Temp);
    BigIntMultiply(result, radix, result);
    BigIntAdd(result, &Temp, result);
    BigIntDivide(&argum, radix, &argum);
  }
  return EXPR_OK;
}

static enum eExprErr ShiftLeft(BigInteger* first, BigInteger *second, BigInteger *result)
{
  int ctr;
  int prevLimb, curLimb;
  int *ptrDest, *ptrSrc;
  int shiftCtr = second->limbs[0].x;
  int delta = shiftCtr / BITS_PER_GROUP;
  int rem = shiftCtr % BITS_PER_GROUP;
  int nbrLimbs = first->nbrLimbs;
  if (second->sign == SIGN_POSITIVE)
  {     // Perform shift left.
    if (second->nbrLimbs > 1)
    {   // Shift too much to the left.
      return EXPR_INTERM_TOO_HIGH;
    }
    if ((unsigned int)(first->nbrLimbs * BITS_PER_GROUP + shiftCtr) > 66438)
    {   // Shift too much to the left.
      return EXPR_INTERM_TOO_HIGH;
    }
    prevLimb = 0;
    ptrSrc = &first->limbs[nbrLimbs - 1].x;
    curLimb = *ptrSrc;
    ptrDest = &first->limbs[nbrLimbs+delta].x;

    for (ctr = nbrLimbs; ctr > 0; ctr--)
    {  // Process starting from most significant limb.
      *ptrDest-- = ((curLimb >> (BITS_PER_GROUP - rem)) | (prevLimb << rem)) & MAX_INT_NBR;
      prevLimb = curLimb;
      curLimb = *--ptrSrc;
    }
    *ptrDest = ((curLimb >> (BITS_PER_GROUP - rem)) | (prevLimb << rem)) & MAX_INT_NBR;
    if (delta > 0)
    {
      memset(first->limbs, 0, delta * sizeof(limb));
    }
    result->nbrLimbs += delta;
    if (result->limbs[result->nbrLimbs].x)
    {
      result->nbrLimbs++;
    }
  }
  else
  {     // Perform shift right.
    int isNegative = 0;
    if (second->nbrLimbs > 1 || shiftCtr > first->nbrLimbs * BITS_PER_GROUP)
    {   // Shift too much to the right. Result is zero or -1.
      if (first->sign == SIGN_POSITIVE)
      {
        intToBigInteger(result, 0);
      }
      else
      {
        intToBigInteger(result, -1);
      }
      return EXPR_OK;
    }
    if (first->sign == SIGN_NEGATIVE)
    {   // If it is negative, add 1, perform shift right, and finally subtract 1 to result.
      isNegative = 1;
      addbigint(first, 1);
    }
    // Shift right the absolute value.
    first->limbs[nbrLimbs].x = 0;
    ptrSrc = &first->limbs[delta+1].x;
    prevLimb = *(ptrSrc-1);
    curLimb = *ptrSrc;
    ptrDest = &first->limbs[0].x;

    for (ctr = delta; ctr < nbrLimbs; ctr++)
    {  // Process starting from least significant limb.
      *ptrDest++ = ((prevLimb >> rem) | (curLimb << (BITS_PER_GROUP - rem))) & MAX_INT_NBR;
      prevLimb = curLimb;
      curLimb = *++ptrSrc;
    }
    *ptrDest = ((prevLimb >> rem) | (curLimb << (BITS_PER_GROUP - rem))) & MAX_INT_NBR;
    result->nbrLimbs -= delta + 1;
    if (result->nbrLimbs == 0 || result->limbs[result->nbrLimbs].x)
    {
      result->nbrLimbs++;
    }
    if (isNegative)
    {    // Adjust negative number.
      addbigint(result, -1);
    }
  }
  return EXPR_OK;
}

#ifndef __EMSCRIPTEN__
void databack(char *data)
{
  (void *)data;
}
#endif
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
#include <stdlib.h>
#include <math.h>
#include "bignbr.h"
#include "expression.h"
#include "factor.h"
#define PAREN_STACK_SIZE           5000
#define COMPR_STACK_SIZE        1000000
#define OPER_POWER                    0
#define OPER_MULTIPLY                 1
#define OPER_DIVIDE                   2
#define OPER_REMAINDER                3
#define OPER_UNARY_MINUS              4
#define OPER_PLUS                     5
#define OPER_MINUS                    6
#define OPER_SHR                      7
#define OPER_SHL                      8
#define OPER_NOT_GREATER              9
#define OPER_NOT_LESS                10
#define OPER_NOT_EQUAL               11
#define OPER_EQUAL                   12
#define OPER_GREATER                 13
#define OPER_LESS                    14
#define OPER_NOT                     15
#define OPER_AND                     16
#define OPER_OR                      17
#define OPER_XOR                     18
#define MAXIMUM_OPERATOR             18

#define COMPUTE_NEXT_PRIME_SIEVE_SIZE 2000

static limb comprStackValues[COMPR_STACK_SIZE];
static int comprStackOffset[PAREN_STACK_SIZE];
static char stackOperators[PAREN_STACK_SIZE];
static limb fibon2[MAX_LEN];
extern limb MontgomeryR1[MAX_LEN];
static int stackIndex;
static int exprIndex;
static int exprLength;
static bool doComputeSubExpression = true;
static int computeSubExprStackThreshold = 0;
#ifndef lang  
  bool lang;
#endif
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
static int ComputeNumFact(void);
static int ComputeMaxFact(void);
static int ComputeMinFact(void);
static int ComputeNumDivs(void);
static int ComputeSumDivs(void);
static int ComputeConcatFact(void);
static char textFactor[MAX_LEN*12];
#endif
static int ComputeSumDigits(void);
static int ComputeRevDigits(void);
static int ComputeNumDigits(void);
static enum eExprErr ComputeModInv(void);
static enum eExprErr ComputeFibLucas(int origValue);
static enum eExprErr ComputePartition(void);
static enum eExprErr ComputeExpr(char *expr, BigInteger *ExpressionResult);
static enum eExprErr ShiftLeft(BigInteger* first, BigInteger *second, BigInteger *result);
static enum eExprErr func(char *expr, BigInteger *ExpressionResult,
  char *funcName, int funcArgs, int leftNumberFlag);
static int type;
static bool valueXused;
static BigInteger curStack;
static BigInteger curStack2;
static BigInteger curStack3;

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
  int nbrParen = 0;
  char* ptrExpr = expr;
  // Check that the parentheses are balanced.
  while ((*ptrExpr != 0) && (*ptrExpr != ';'))
  {
    if (*ptrExpr == '(')
    {
      nbrParen++;
    }
    else if (*ptrExpr == ')')
    {
      nbrParen--;
    }
    if (nbrParen < 0)
    {    // More closing parentheses than opening ones.
      return EXPR_PAREN_MISMATCH;
    }
    ptrExpr++;
  }
  if (nbrParen != 0)
  {      // Number of open and closing parentheses do not match.
    return EXPR_PAREN_MISMATCH;
  }
  doComputeSubExpression = true;
  valueXused = false;
  stackIndex = 0;
  comprStackOffset[0] = 0;
  exprIndex = 0;
  type = typ;
  retcode = ComputeExpr(expr, ExpressionResult);
  if (retcode != 0) { return retcode; }
#ifdef FACTORIZATION_APP
  if (ExpressionResult[0].nbrLimbs > 332192 / BITS_PER_GROUP + 1)   // 100000/log_10(2) = 332192
#else
  if (ExpressionResult[0].nbrLimbs > 33219 / BITS_PER_GROUP + 1)    // 10000/log_10(2) = 33219
#endif
  {
    return EXPR_NUMBER_TOO_HIGH;
  }
  if ((valueX.nbrLimbs > 0) && !valueXused)
  {
    return EXPR_VAR_OR_COUNTER_REQUIRED;
  }
  return 0;
}

static int numLimbs(int* pLen)
{
  int nbrLimbs = *pLen;
  if (nbrLimbs < 0)
  {
    nbrLimbs = -nbrLimbs;
  }
  return nbrLimbs;
}

static void getCurrentStackValue(BigInteger* pValue)
{
  limb* ptrStackValue = &comprStackValues[comprStackOffset[stackIndex]];
  NumberLength = numLimbs((int *)ptrStackValue);
  IntArray2BigInteger((int *)ptrStackValue, pValue);
}

enum eExprErr setStackValue(BigInteger* pValue)
{
  int currentOffset = comprStackOffset[stackIndex];
  if (currentOffset >= COMPR_STACK_SIZE - sizeof(BigInteger) / sizeof(limb))
  {
    return EXPR_OUT_OF_MEMORY;
  }
  NumberLength = pValue->nbrLimbs;
  BigInteger2IntArray((int*)&comprStackValues[currentOffset], pValue);
  return EXPR_OK;
}

static enum eExprErr ComputeExpr(char *expr, BigInteger *ExpressionResult)
{
  int shLeft;
  enum eExprErr retcode;
  limb carry;
  bool leftNumberFlag = false;
  int exprIndexAux;
  int offset;
  enum eExprErr SubExprResult;
  int len;
  limb largeLen;
  limb *ptrLimb;
  BigInteger *pBigInt;
  int c;
  int startStackIndex = stackIndex;
  
  exprLength = (int)strlen(expr);
  while (exprIndex < exprLength)
  {
    char charValue;

    charValue = *(expr+exprIndex);
    if ((charValue == ' ') || (charValue == 9))
    {           // Ignore spaces and horizontal tabs.
      exprIndex++;
      continue;
    }
    if (charValue == '^')
    {           // Caret is exponentiation operation.
      charValue = OPER_POWER;
      exprIndex++;
    }
    else if ((charValue == '*') && (*(expr + exprIndex + 1) == '*'))
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
    else if ((charValue == '<') && (*(expr + exprIndex + 1) == '='))
    {
      charValue = OPER_NOT_GREATER;
      exprIndex += 2;
    }
    else if ((charValue == '>') && (*(expr + exprIndex + 1) == '='))
    {
      charValue = OPER_NOT_LESS;
      exprIndex += 2;
    }
    else if ((charValue == '!') && (*(expr + exprIndex + 1) == '='))
    {
      charValue = OPER_NOT_EQUAL;
      exprIndex += 2;
    }
    else if ((charValue == '=') && (*(expr + exprIndex + 1) == '='))
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
    else if (((charValue & 0xDF) == 'N') && ((*(expr + exprIndex + 1) & 0xDF) == 'O') &&
      ((*(expr + exprIndex + 2) & 0xDF) == 'T'))
    {
      charValue = OPER_NOT;
      exprIndex += 3;
    }
    else if (((charValue & 0xDF) == 'A') && ((*(expr + exprIndex + 1) & 0xDF)) == 'N' &&
      ((*(expr + exprIndex + 2) & 0xDF) == 'D'))
    {
      charValue = OPER_AND;
      exprIndex += 3;
    }
    else if (((charValue & 0xDF) == 'O') && ((*(expr + exprIndex + 1) & 0xDF) == 'R'))
    {
      charValue = OPER_OR;
      exprIndex += 2;
    }
    else if (((charValue & 0xDF) == 'X') && ((*(expr + exprIndex + 1) & 0xDF) == 'O') &&
      ((*(expr + exprIndex + 2) & 0xDF) == 'R'))
    {
      charValue = OPER_XOR;
      exprIndex += 3;
    }
    else if (((charValue & 0xDF) == 'S') && ((*(expr + exprIndex + 1) & 0xDF) == 'H') &&
      ((*(expr + exprIndex + 2) & 0xDF) == 'L'))
    {
      charValue = OPER_SHL;
      exprIndex += 3;
    }
    else if (((charValue & 0xDF) == 'S') && ((*(expr + exprIndex + 1) & 0xDF) == 'H') &&
      ((*(expr + exprIndex + 2) & 0xDF) == 'R'))
    {
      charValue = OPER_SHR;
      exprIndex += 3;
    }
    else if (charValue == '!')
    {           // Calculating factorial.
      if (!leftNumberFlag)
      {
        return EXPR_SYNTAX_ERROR;
      }
      getCurrentStackValue(&curStack);
      if (curStack.nbrLimbs != 1)
      {
        return EXPR_INTERM_TOO_HIGH;
      }
#ifdef FACTORIZATION_APP
      if ((curStack.limbs[0].x < 0) || (curStack.limbs[0].x >= 47177))
#else
      if ((curStack.limbs[0].x < 0) || (curStack.limbs[0].x >= 5984))
#endif
      {
        return EXPR_INTERM_TOO_HIGH;
      }
      factorial(&curStack, (int)curStack.limbs[0].x);
      retcode = setStackValue(&curStack);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      exprIndex++;
      continue;
    }
    else if (charValue == '#')
    {           // Calculating primorial.
      if (!leftNumberFlag)
      {
        return EXPR_SYNTAX_ERROR;
      }
      if (doComputeSubExpression)
      {
        getCurrentStackValue(&curStack);
        if (curStack.nbrLimbs > 2)
        {
          return EXPR_INTERM_TOO_HIGH;
        }
        if (curStack.nbrLimbs == 2)
        {
          largeLen.x = curStack.limbs[0].x +
            (curStack.limbs[1].x << BITS_PER_GROUP);
        }
        else
        {
          largeLen.x = curStack.limbs[0].x;
        }
        len = (int)largeLen.x;
#ifdef FACTORIZATION_APP
        if ((len < 0) || (len > 460490))
#else
        if ((len < 0) || (len > 46049))
#endif
        {
          return EXPR_INTERM_TOO_HIGH;
        }
        primorial(&curStack, (int)curStack.limbs[0].x);
        retcode = setStackValue(&curStack);
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
      }
      exprIndex++;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "GCD", 2, leftNumberFlag)) != EXPR_NOT_FOUND)
    {
      if (retcode != EXPR_OK) { return retcode; }
      if (doComputeSubExpression)
      {
        getCurrentStackValue(&curStack);    // Get first argument of GCD.
        stackIndex++;
        getCurrentStackValue(&curStack2);   // Get second argument of GCD.
        stackIndex--;
        BigIntGcd(&curStack, &curStack2, &curStack);
        retcode = setStackValue(&curStack);
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
      }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "MODPOW", 3, leftNumberFlag)) != EXPR_NOT_FOUND)
    {
      if (retcode != EXPR_OK) { return retcode; }
      if (doComputeSubExpression)
      {
        getCurrentStackValue(&curStack);    // Get first argument of MODPOW.
        stackIndex++;
        getCurrentStackValue(&curStack2);   // Get second argument of MODPOW.
        stackIndex++;
        getCurrentStackValue(&curStack3);   // Get third argument of MODPOW.
        stackIndex -= 2;
        retcode = BigIntGeneralModularPower(&curStack, &curStack2, &curStack3, &curStack);
        if (retcode != EXPR_OK) { return retcode; }
        retcode = setStackValue(&curStack);
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
      }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "MODINV", 2, leftNumberFlag)) != EXPR_NOT_FOUND)
    {
      if (retcode != 0) { return retcode; }
      if (doComputeSubExpression)
      {
        retcode = ComputeModInv();
        if (retcode != 0) { return retcode; }
        retcode = setStackValue(&curStack);
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
      }
      leftNumberFlag = 1;
      continue;
    }
#ifdef FACTORIZATION_FUNCTIONS
    else if ((retcode = func(expr, ExpressionResult,
      "TOTIENT", 1, leftNumberFlag)) != EXPR_NOT_FOUND)
    {
      if (retcode != 0) { return retcode; }
      if (doComputeSubExpression)
      {
        retcode = ComputeTotient();
        if (retcode != 0) { return retcode; }
        retcode = setStackValue(&curStack);
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
      }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "NUMDIVS", 1, leftNumberFlag)) != EXPR_NOT_FOUND)
    {
      if (retcode != 0) { return retcode; }
      if (doComputeSubExpression)
      {
        retcode = ComputeNumDivs();
        if (retcode != 0) { return retcode; }
        retcode = setStackValue(&curStack);
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
      }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "SUMDIVS", 1, leftNumberFlag)) != EXPR_NOT_FOUND)
    {
      if (retcode != 0) { return retcode; }
      if (doComputeSubExpression)
      {
        retcode = ComputeSumDivs();
        if (retcode != 0) { return retcode; }
        retcode = setStackValue(&curStack);
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
        leftNumberFlag = 1;
      }
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "MINFACT", 1, leftNumberFlag)) != EXPR_NOT_FOUND)
    {
      if (retcode != 0) { return retcode; }
      if (doComputeSubExpression)
      {
        retcode = ComputeMinFact();
        if (retcode != 0) { return retcode; }
        retcode = setStackValue(&curStack);
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
      }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "MAXFACT", 1, leftNumberFlag)) != EXPR_NOT_FOUND)
    {
      if (retcode != 0) { return retcode; }
      if (doComputeSubExpression)
      {
        retcode = ComputeMaxFact();
        if (retcode != 0) { return retcode; }
        retcode = setStackValue(&curStack);
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
      }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "NUMFACT", 1, leftNumberFlag)) != EXPR_NOT_FOUND)
    {
      if (retcode != 0) { return retcode; }
      if (doComputeSubExpression)
      {
        retcode = ComputeNumFact();
        if (retcode != 0) { return retcode; }
        retcode = setStackValue(&curStack);
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
      }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "CONCATFACT", 2, leftNumberFlag)) != EXPR_NOT_FOUND)
    {
      if (retcode != 0) { return retcode; }
      if (doComputeSubExpression)
      {
        retcode = ComputeConcatFact();
        if (retcode != 0) { return retcode; }
        retcode = setStackValue(&curStack);
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
      }
      leftNumberFlag = 1;
      continue;
    }
#endif
    else if ((retcode = func(expr, ExpressionResult,
      "SUMDIGITS", 2, leftNumberFlag)) != EXPR_NOT_FOUND)
    {
      if (retcode != 0) { return retcode; }
      if (doComputeSubExpression)
      {
        retcode = ComputeSumDigits();
        if (retcode != 0) { return retcode; }
        retcode = setStackValue(&curStack);
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
      }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "NUMDIGITS", 2, leftNumberFlag)) != EXPR_NOT_FOUND)
    {
      if (retcode != 0) { return retcode; }
      if (doComputeSubExpression)
      {
        retcode = ComputeNumDigits();
        if (retcode != 0) { return retcode; }
        retcode = setStackValue(&curStack);
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
      }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "REVDIGITS", 2, leftNumberFlag)) != EXPR_NOT_FOUND)
    {
      if (retcode != 0) { return retcode; }
      if (doComputeSubExpression)
      {
        retcode = ComputeRevDigits();
        if (retcode != 0) { return retcode; }
        retcode = setStackValue(&curStack);
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
      }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "ISPRIME", 1, leftNumberFlag)) != EXPR_NOT_FOUND)
    {
      if (retcode != 0) { return retcode; }
      if (doComputeSubExpression)
      {
        getCurrentStackValue(&curStack);    // Get argument.
#ifdef FACTORIZATION_APP
        if (BpswPrimalityTest(&curStack, NULL) == 0)
#else
        if (BpswPrimalityTest(&curStack) == 0)
#endif
        {    // Argument is a probable prime.
          intToBigInteger(&curStack, -1);
        }
        else
        {    // Argument is not a probable prime.
          intToBigInteger(&curStack, 0);
        }
        retcode = setStackValue(&curStack);
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
      }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "F", 1, leftNumberFlag)) != EXPR_NOT_FOUND)
    {
      if (retcode != 0) { return retcode; }
      if (doComputeSubExpression)
      {
        getCurrentStackValue(&curStack);    // Get argument.
        retcode = ComputeFibLucas(0);
        if (retcode != 0) { return retcode; }
        retcode = setStackValue(&curStack);
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
      }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "L", 1, leftNumberFlag)) != EXPR_NOT_FOUND)
    {
      if (retcode != 0) { return retcode; }
      if (doComputeSubExpression)
      {
        getCurrentStackValue(&curStack);    // Get argument.
        retcode = ComputeFibLucas(2);
        if (retcode != 0) { return retcode; }
        retcode = setStackValue(&curStack);
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
      }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "P", 1, leftNumberFlag)) != EXPR_NOT_FOUND)
    {
      if (retcode != 0) { return retcode; }
      if (doComputeSubExpression)
      {
        getCurrentStackValue(&curStack);    // Get argument.
        retcode = ComputePartition();
        if (retcode != 0) { return retcode; }
        retcode = setStackValue(&curStack);
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
      }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "N", 1, leftNumberFlag)) != EXPR_NOT_FOUND)
    {
      if (retcode != 0) { return retcode; }
      if (doComputeSubExpression)
      {
        getCurrentStackValue(&curStack);    // Get argument.
        retcode = ComputeNext();
        if (retcode != 0) { return retcode; }
        retcode = setStackValue(&curStack);
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
      }
      leftNumberFlag = 1;
      continue;
    }
    else if ((retcode = func(expr, ExpressionResult,
      "B", 1, leftNumberFlag)) != EXPR_NOT_FOUND)
    {
      if (retcode != 0) { return retcode; }
      if (doComputeSubExpression)
      {
        getCurrentStackValue(&curStack);    // Get argument.
        retcode = ComputeBack();
        if (retcode != 0) { return retcode; }
        retcode = setStackValue(&curStack);
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
      }
      leftNumberFlag = 1;
      continue;
    }
    else if ((charValue & 0xDF) == 'X')
    {
      if (leftNumberFlag || (valueX.nbrLimbs == 0))
      {
        return EXPR_SYNTAX_ERROR;
      }
      valueXused = true;
      exprIndex++;
      if (stackIndex > 0)
      {
        int currentStackOffset = comprStackOffset[stackIndex - 1];
        comprStackOffset[stackIndex] = currentStackOffset +
          numLimbs((int *)&comprStackValues[currentStackOffset]) + 1;
      }
      retcode = setStackValue(&valueX);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      leftNumberFlag = true;
      continue;
    }
    else if ((charValue & 0xDF) == 'C')
    {
      if (leftNumberFlag || (valueX.nbrLimbs == 0))
      {
        return EXPR_SYNTAX_ERROR;
      }
      intToBigInteger(&curStack, counterC);
      valueXused = true;
      exprIndex++;
      if (stackIndex > 0)
      {
        int currentStackOffset = comprStackOffset[stackIndex - 1];
        comprStackOffset[stackIndex] = currentStackOffset +
          numLimbs((int*)&comprStackValues[currentStackOffset]) + 1;
      }
      retcode = setStackValue(&curStack);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      leftNumberFlag = true;
      continue;
    }
    else if (charValue == '(')
    {
      if (leftNumberFlag)
      {
        return EXPR_SYNTAX_ERROR;
      }
      if (stackIndex >= PAREN_STACK_SIZE)
      {
        return EXPR_TOO_MANY_PAREN;
      }
      comprStackOffset[stackIndex + 1] = comprStackOffset[stackIndex];
      stackOperators[stackIndex++] = charValue;
      exprIndex++;
      continue;
    }
    else if ((charValue == ')') || (charValue == ','))
    {
      int curStackOperand;
      if (leftNumberFlag == 0)
      {       // Previous item should be a number or variable.
        return EXPR_SYNTAX_ERROR;
      }
      while ((stackIndex > startStackIndex) &&
        (stackOperators[stackIndex - 1] != '('))
      {
        if ((SubExprResult = ComputeSubExpr()) != 0)
        {
          return SubExprResult;
        }
        if ((stackOperators[stackIndex-1] == OPER_AND) ||
          (stackOperators[stackIndex-1] == OPER_OR))
        {
          if (!doComputeSubExpression &&
              (stackIndex == computeSubExprStackThreshold))
          {
            doComputeSubExpression = true;
          }
        }
      }
      curStackOperand = comprStackOffset[stackIndex];
      if (stackIndex == startStackIndex)
      {
        comprStackOffset[stackIndex+1] = curStackOperand +
          numLimbs(&comprStackValues[curStackOperand].x) + 1;
        break;
      }
      if (charValue == ',')
      {
        return EXPR_PAREN_MISMATCH;
      }
      comprStackOffset[stackIndex - 1] = curStackOperand;
      comprStackOffset[stackIndex] += numLimbs(&comprStackValues[curStackOperand].x) + 1;
      stackIndex--;    /* Discard ')' */
      leftNumberFlag = 1;
      exprIndex++;
      continue;
    }
    else if ((charValue >= '0') && (charValue <= '9'))
    {
      int currentStackOffset;
      exprIndexAux = exprIndex;
      if ((charValue == '0') && (exprIndexAux < exprLength - 2) &&
          (*(expr+exprIndexAux + 1) == 'x'))
      {  // hexadecimal
        int exprIndexFirstHexDigit = -1;
        exprIndexAux++;
        while (exprIndexAux < exprLength - 1)
        {
          charValue = *(expr+exprIndexAux + 1);
          if (((charValue >= '0') && (charValue <= '9')) ||
              ((charValue >= 'A') && (charValue <= 'F')) ||
              ((charValue >= 'a') && (charValue <= 'f')))
          {
            exprIndexAux++;
            if ((charValue != '0') && (exprIndexFirstHexDigit < 0))
            {
              exprIndexFirstHexDigit = exprIndexAux;
            }
          }
          else
          {
            break;
          }
        }
        if (exprIndexFirstHexDigit < 0)
        {    // Number is zero
          intToBigInteger(&curStack, 0);
          exprIndex = exprIndexAux + 1;
        }
        else
        {    // Generate big integer from hexadecimal number from right to left.
          carry.x = 0;
          shLeft = 0;
          offset = exprIndexAux;
          ptrLimb = &curStack.limbs[0];
          for (; exprIndexAux >= exprIndexFirstHexDigit; exprIndexAux--)
          {
            c = *(expr + exprIndexAux);
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
            carry.x += c << shLeft;
            shLeft += 4;   // 4 bits per hex digit.
            if (shLeft >= BITS_PER_GROUP)
            {
              shLeft -= BITS_PER_GROUP;
              (ptrLimb++)->x = carry.x & MAX_VALUE_LIMB;
              carry.x = c >> (4 - shLeft);
            }
          }
          if ((carry.x != 0) || (ptrLimb == &curStack.limbs[0]))
          {
            (ptrLimb++)->x = carry.x;
          }
          exprIndex = offset + 1;
          curStack.nbrLimbs = (int)(ptrLimb - &curStack.limbs[0]);
          curStack.sign = SIGN_POSITIVE;
        }
      }
      else
      {                   // Decimal number.
        while (exprIndexAux < exprLength - 1)
        {
          charValue = *(expr+exprIndexAux + 1);
          if ((charValue >= '0') && (charValue <= '9'))
          {
            exprIndexAux++;
          }
          else
          {
            break;
          }
        }
        // Generate big integer from decimal number
        pBigInt = &curStack;
        Dec2Bin(expr + exprIndex, pBigInt->limbs,
                exprIndexAux + 1 - exprIndex, &pBigInt -> nbrLimbs);
        pBigInt -> sign = SIGN_POSITIVE;
        exprIndex = exprIndexAux + 1;
      }
      retcode = setStackValue(&curStack);   // Push number onto stack.
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      currentStackOffset = comprStackOffset[stackIndex];
      comprStackOffset[stackIndex + 1] = currentStackOffset +
          numLimbs((int*)&comprStackValues[currentStackOffset]) + 1;
      leftNumberFlag = true;
      continue;
    }
    if (charValue <= MAXIMUM_OPERATOR)
    {
      if (((charValue == OPER_PLUS) || (charValue == OPER_MINUS)) && (leftNumberFlag == false))
      {                    // Unary plus/minus operator
        if (charValue == OPER_PLUS)
        {
          continue;
        }
        else
        {
          if ((stackIndex > startStackIndex) && (stackOperators[stackIndex - 1] == OPER_UNARY_MINUS))
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
        while ((stackIndex > startStackIndex) &&
          (stackOperators[stackIndex - 1] != '(') &&
          (priority[(int)stackOperators[stackIndex - 1]] <= priority[(int)charValue]))
        {
          if ((SubExprResult = ComputeSubExpr()) != 0)
          {
            return SubExprResult;
          }
          if ((stackOperators[stackIndex-1] == OPER_AND) ||
            (stackOperators[stackIndex-1] == OPER_OR))
          {
            if (!doComputeSubExpression &&
              (stackIndex == computeSubExprStackThreshold))
            {
              doComputeSubExpression = true;
            }
          }
        }
        if (stackIndex >= 0)
        {
          int currentStackOffset = comprStackOffset[stackIndex];
          comprStackOffset[stackIndex+1] = currentStackOffset +
            numLimbs((int*)&comprStackValues[currentStackOffset]) + 1;
        }
      }
      if (charValue == OPER_AND)
      {
        getCurrentStackValue(&curStack);
        if (BigIntIsZero(&curStack))
        {
          doComputeSubExpression = false;
          computeSubExprStackThreshold = stackIndex;
        }
      }
      else if (charValue == OPER_OR)
      {
        getCurrentStackValue(&curStack);
        if ((curStack.sign == SIGN_NEGATIVE) && (curStack.nbrLimbs == 1) &&
          (curStack.limbs[0].x == 1))
        {       // Number is -1.
          doComputeSubExpression = false;
          computeSubExprStackThreshold = stackIndex;
        }
      }
      stackOperators[stackIndex++] = charValue;
      leftNumberFlag = 0;
      continue;
    }
    return EXPR_SYNTAX_ERROR;
  }                              /* end while */
  if (!leftNumberFlag)
  {
    return EXPR_SYNTAX_ERROR;
  }
  while ((stackIndex > startStackIndex) && (stackOperators[stackIndex - 1] != '('))
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
  getCurrentStackValue(ExpressionResult);
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
  enum eExprErr retcode;
  char stackOper = stackOperators[--stackIndex];
  getCurrentStackValue(&curStack);    // Get first argument of GCD.
  stackIndex++;
  getCurrentStackValue(&curStack2);   // Get second argument of GCD.
  stackIndex--;
  BigInteger *firstArg = &curStack;
  BigInteger *secondArg = &curStack2;
  BigInteger *result = &curStack;
  switch (stackOper)
  {
  case OPER_PLUS:
    BigIntAdd(firstArg, secondArg, result);
    break;
  case OPER_MINUS:
    BigIntSubt(firstArg, secondArg, result);
    break;
  case OPER_UNARY_MINUS:
    BigIntNegate(secondArg, result);
    break;
  case OPER_DIVIDE:
    retcode = BigIntDivide(firstArg, secondArg, result);
    if (retcode != EXPR_OK)
    {
      return retcode;
    }
    break;
  case OPER_MULTIPLY:
    retcode = BigIntMultiply(firstArg, secondArg, result);
    if (retcode != EXPR_OK)
    {
      return retcode;
    }
    break;
  case OPER_REMAINDER:
    retcode = BigIntRemainder(firstArg, secondArg, result);
    if (retcode != EXPR_OK)
    {
      return retcode;
    }
    break;
  case OPER_POWER:
    retcode = BigIntPower(firstArg, secondArg, result);
    if (retcode != EXPR_OK)
    {
      return retcode;
    }
    break;
  case OPER_EQUAL:
    BigIntSubt(firstArg, secondArg, result);
    intToBigInteger(result, ((result->nbrLimbs == 1) && (result->limbs[0].x == 0? -1: 0)));
    break;
  case OPER_NOT_EQUAL:
    BigIntSubt(firstArg, secondArg, result);
    intToBigInteger(result, ((result->nbrLimbs == 1) && (result->limbs[0].x == 0 ? 0 : -1)));
    break;
  case OPER_GREATER:
    BigIntSubt(secondArg, firstArg, result);
    intToBigInteger(result, (result->sign == SIGN_NEGATIVE) ? -1 : 0);
    break;
  case OPER_NOT_GREATER:
    BigIntSubt(secondArg, firstArg, result);
    intToBigInteger(result, (result->sign == SIGN_NEGATIVE) ? 0 : -1);
    break;
  case OPER_LESS:
    BigIntSubt(firstArg, secondArg, result);
    intToBigInteger(result, (result->sign == SIGN_NEGATIVE) ? -1 : 0);
    break;
  case OPER_NOT_LESS:
    BigIntSubt(firstArg, secondArg, result);
    intToBigInteger(result, (result->sign == SIGN_NEGATIVE) ? 0 : -1);
    break;
  case OPER_SHL:
    ShiftLeft(firstArg, secondArg, result);
    break;
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
    break;
  case OPER_NOT:    // Perform binary NOT as result <- -1 - argument.
    intToBigInteger(firstArg, -1);
    BigIntSubt(firstArg, secondArg, result);
    break;
  case OPER_AND:    // Perform binary AND.
    BigIntAnd(firstArg, secondArg, result);
    break;
  case OPER_OR:     // Perform binary OR.
    BigIntOr(firstArg, secondArg, result);
    break;
  case OPER_XOR:    // Perform binary XOR.
    BigIntXor(firstArg, secondArg, result);
    break;
  }
  retcode = setStackValue(&curStack);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  return EXPR_OK;
}

static enum eExprErr func(char *expr, BigInteger *ExpressionResult,
  char *funcName, int funcArgs, int leftNumberFlag)
{
  int index;
  int funcNameLen = (int)strlen(funcName);
  char *ptrExpr;
  char *ptrFuncName;

  if (exprIndex + funcNameLen > exprLength)
  {
    return EXPR_NOT_FOUND;
  }
  ptrExpr = expr + exprIndex;
  ptrFuncName = funcName;
  while (*ptrFuncName != 0)
  {
    if ((*ptrExpr & 0xDF) != *ptrFuncName)
    {
      return EXPR_NOT_FOUND;
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
  if ((exprIndex == exprLength) || (*(expr + exprIndex++) != '('))
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
    if ((exprIndex == exprLength) || (*(expr + exprIndex++) != compareChar))
    {
      return EXPR_SYNTAX_ERROR;
    }
    stackIndex++;
  }
  stackIndex -= funcArgs;
  return EXPR_OK;
}

static void generateSieve(int* pSmallPrimes, char* sieve, BigInteger* pArgument, bool isNext)
{
  int ctr;
  // Indicate numbers not divisible by small primes in advance.
  (void)memset(sieve, 0, COMPUTE_NEXT_PRIME_SIEVE_SIZE);
  for (ctr = 0; ctr < 1229; ctr++)
  {     // For each prime less than 10000...
    int prime = *pSmallPrimes++;
    int remainder = getRemainder(pArgument, prime);
    // Compute first element of sieve to indicate multiple of prime.
    if (isNext)
    {
      if (remainder > 0)
      {
        remainder = prime - remainder;
      }
    }
    else
    {
      remainder = (COMPUTE_NEXT_PRIME_SIEVE_SIZE - remainder) % prime;
    }
    if (remainder < 0)
    {
      remainder += prime;
    }
    for (; remainder < COMPUTE_NEXT_PRIME_SIEVE_SIZE; remainder += prime)
    {        // Indicate numbers are divisible by this prime.
      sieve[remainder] = 1;
    }
  }
}
// Compute previous probable prime
static int ComputeBack(void)
{
  char sieve[COMPUTE_NEXT_PRIME_SIEVE_SIZE];
  BigInteger *pArgument = &curStack;
  BigInteger *pResult = &curStack;
  limb *pResultLimbs = pResult->limbs;
  limb *pArgumentLimbs = pArgument->limbs;
  pResult->sign = SIGN_POSITIVE;
  if (pArgument->sign == SIGN_NEGATIVE)
  {
    return EXPR_INVALID_PARAM;
  }
  if (pArgument->nbrLimbs == 1)
  {
    if (pArgumentLimbs->x < 3)
    {
      return EXPR_INVALID_PARAM;
    }
    if (pArgumentLimbs->x == 3)
    {
      pResult->nbrLimbs = 1;
      pResultLimbs->x = 2;
      return EXPR_OK;
    }
  }
  initializeSmallPrimes(smallPrimes);
  CopyBigInt(pResult, pArgument);
  pResultLimbs->x |= 1;  // If number is even, use next odd number.
  if (pResult->nbrLimbs == 1)
  {
    do
    {        // Loop that searches for previous probable prime.
      addbigint(pResult, -2);
#ifdef FACTORIZATION_APP
    } while (BpswPrimalityTest(pResult, NULL));  // Continue loop if not probable prime.
#else
    } while (BpswPrimalityTest(pResult));        // Continue loop if not probable prime.
#endif
  }
  else
  {          // Big number: use sieve.
    for (;;)
    {        // Loop that searches for previous probable prime.
      int ctr;
      addbigint(pResult, -COMPUTE_NEXT_PRIME_SIEVE_SIZE);
      generateSieve(smallPrimes, sieve, pResult, true);
      for (ctr = COMPUTE_NEXT_PRIME_SIEVE_SIZE-1; ctr >= 0; ctr--)
      {
        if (sieve[ctr] == 0)
        {   // Number is not divisible by primes less than 1000.
          addbigint(pResult, ctr);
#ifdef FACTORIZATION_APP
          if (BpswPrimalityTest(pResult, NULL) == 0)  // End loop if probable prime.
#else
          if (BpswPrimalityTest(pResult) == 0)        // End loop if probable prime.
#endif
          {
            return EXPR_OK;
          }
          addbigint(pResult, -ctr);
        }
      }
    }
  }
  return EXPR_OK;
}

// Compute next probable prime
static int ComputeNext(void)
{
  char sieve[COMPUTE_NEXT_PRIME_SIEVE_SIZE];
  BigInteger *pArgument = &curStack;
  BigInteger *pResult = &curStack;
  limb *pResultLimbs = pResult->limbs;
  limb *pArgumentLimbs = pArgument->limbs;
  pResult->sign = SIGN_POSITIVE;
  if ((pArgument->sign == SIGN_NEGATIVE) ||
    ((pArgument->nbrLimbs == 1) && (pArgumentLimbs->x < 2)))
  {
    pResult->nbrLimbs = 1;
    pResultLimbs->x = 2;
    return EXPR_OK;
  }
  initializeSmallPrimes(smallPrimes);
  CopyBigInt(pResult, pArgument);
  if ((pResult->limbs[0].x & 1) == 0)
  {   // Number is even.
    addbigint(pResult, 1);
  }
  else
  {   // Number is odd.
    addbigint(pResult, 2);
  }
  if (pResult->nbrLimbs == 1)
  {
    for (;;)
    {        // Loop that searches for next probable prime.
#ifdef FACTORIZATION_APP
      if (BpswPrimalityTest(pResult, NULL) == 0)  // Continue loop if not probable prime.
#else
      if (BpswPrimalityTest(pResult) == 0)        // Continue loop if not probable prime.
#endif
      {
        return EXPR_OK;
      }
      addbigint(pResult, 2);
    }
  }
  else
  {          // Big number: use sieve.
    for (;;)
    {        // Loop that searches for next probable prime.
      int ctr;
      generateSieve(smallPrimes, sieve, pResult, true);
      for (ctr = 0; ctr<COMPUTE_NEXT_PRIME_SIEVE_SIZE; ctr++)
      {
        if (sieve[ctr] == 0)
        {   // Number is not divisible by primes less than 1000.
          addbigint(pResult, ctr);
#ifdef FACTORIZATION_APP
          if (BpswPrimalityTest(pResult, NULL) == 0)  // End loop if probable prime.
#else
          if (BpswPrimalityTest(pResult) == 0)        // End loop if probable prime.
#endif
          {
            return EXPR_OK;
          }
          addbigint(pResult, -ctr);
        }
      }
      addbigint(pResult, COMPUTE_NEXT_PRIME_SIEVE_SIZE);
    }
  }
}

static enum eExprErr ComputeModInv(void)
{
  BigInteger one;
  BigInteger* pDiv;
  getCurrentStackValue(&curStack);    // Get first argument of MODINV.
  stackIndex++;
  getCurrentStackValue(&curStack2);   // Get second argument of MODINV.
  stackIndex--;
  pDiv = &curStack2;
  if ((pDiv->nbrLimbs == 1) && (pDiv->limbs[0].x == 0))
  {
    return EXPR_DIVIDE_BY_ZERO;
  }
  // Check that the arguments are relatively prime.
  BigIntGcd(&curStack, &curStack2, &one);
  if ((one.nbrLimbs != 1) || (one.limbs[0].x != 1))
  {
    return EXPR_ARGUMENTS_NOT_RELATIVELY_PRIME;
  }
  intToBigInteger(&one, 1);
  BigIntGeneralModularDivision(&one, &curStack, pDiv, &curStack2);
  CopyBigInt(&curStack, pDiv);
  return EXPR_OK;
}

static enum eExprErr ComputeFibLucas(int origValue)
{
  BigInteger *pArgument = &curStack;
  limb largeVal;
  int val;
  int len;
  limb *pFibonPrev;
  limb *pFibonAct;
  if (pArgument->sign == SIGN_NEGATIVE)
  {
    return EXPR_INVALID_PARAM;
  }
  if (pArgument->nbrLimbs > 1)
  {
    return EXPR_INTERM_TOO_HIGH;
  }
  largeVal.x = pArgument->limbs[0].x;
#ifdef FACTORIZATION_APP
  if (largeVal.x > 956620)
#else
  if (largeVal.x > 95662)
#endif
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
  pArgument = &curStack;
  pArgument->sign = SIGN_POSITIVE;
  pArgument->nbrLimbs = len;
  (void)memcpy(pArgument->limbs, pFibonAct, len * sizeof(limb));
  return EXPR_OK;
}

static enum eExprErr ComputePartition(void)
{
  BigInteger *pArgument = &curStack;
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
  partition(val, &curStack);
  return EXPR_OK;
}

#ifdef FACTORIZATION_FUNCTIONS
static void PerformFactorization(BigInteger *tofactor)
{
  NumberLength = tofactor->nbrLimbs;
  BigInteger2IntArray(nbrToFactor, tofactor);
  if (hexadecimal)
  {
    Bin2Hex(tofactor->limbs, tofactorDec, tofactor->nbrLimbs, groupLen);
  }
  else
  {
    Bin2Dec(tofactor->limbs, tofactorDec, tofactor->nbrLimbs, groupLen);
  }
  factor(tofactor, nbrToFactor, factorsMod, astFactorsMod);
}

static int ComputeTotient(void)
{
  getCurrentStackValue(&curStack);    // Get argument.
  PerformFactorization(&curStack);
  Totient(&curStack);
  return EXPR_OK;
}

static int ComputeNumDivs(void)
{
  getCurrentStackValue(&curStack);    // Get argument.
  PerformFactorization(&curStack);
  NumberOfDivisors(&curStack);
  return EXPR_OK;
}

static int ComputeSumDivs(void)
{
  getCurrentStackValue(&curStack);    // Get argument.
  PerformFactorization(&curStack);
  SumOfDivisors(&curStack);
  return EXPR_OK;
}

static int ComputeNumFact(void)
{
  getCurrentStackValue(&curStack);    // Get argument.
  PerformFactorization(&curStack);
  NumFactors(&curStack);
  return EXPR_OK;
}

static int ComputeMaxFact(void)
{
  getCurrentStackValue(&curStack);    // Get argument.
  PerformFactorization(&curStack);
  MaxFactor(&curStack);
  return EXPR_OK;
}

static int ComputeMinFact(void)
{
  int primeIndex = 0;
  int prime = 2;
  getCurrentStackValue(&curStack);    // Get argument.
  // Try to find a small factor (less than 32767).
  initializeSmallPrimes(smallPrimes);
  do
  {
    if (getRemainder(&curStack, prime) == 0)
    {
      intToBigInteger(&curStack, prime);
      return EXPR_OK;
    }
    prime = smallPrimes[++primeIndex];
  } while (prime < 32767);
  PerformFactorization(&curStack);
  MinFactor(&curStack);
  return EXPR_OK;
}

static int ComputeConcatFact(void)
{
  getCurrentStackValue(&curStack);    // Get first argument.
  stackIndex++;
  getCurrentStackValue(&curStack2);   // Get second argument.
  stackIndex--;
  BigInteger *mode = &curStack;
  BigInteger factorValue;
  int factorNumber;
  int nbrFactors;
  int descend = mode->limbs[0].x & 1;
  int repeated = mode->limbs[0].x & 2;
  char *ptrTextFactor = textFactor;
  if ((mode->nbrLimbs > 1) || (mode->sign == SIGN_NEGATIVE) || (mode->limbs[0].x > 3))
  {      // The valid modes are 0, 1, 2 and 3.
    return EXPR_INVALID_PARAM;
  }
  PerformFactorization(&curStack2);   // Factor second argument.
  nbrFactors = astFactorsMod[0].multiplicity;
  for (factorNumber = 1; factorNumber <= nbrFactors; factorNumber++)
  {
    int ctr;
    struct sFactors* pstFactor = &astFactorsMod[descend ? nbrFactors - factorNumber + 1 : factorNumber];
    NumberLength = *(pstFactor->ptrFactor);
    IntArray2BigInteger(pstFactor->ptrFactor, &factorValue);
    ctr = (repeated ? pstFactor->multiplicity : 1);
    for (; ctr > 0; ctr--)
    {
      BigInteger2Dec(&factorValue, ptrTextFactor, 0);
      ptrTextFactor += strlen(ptrTextFactor);
    }
  }
  if (ptrTextFactor == &textFactor[0])
  {
    intToBigInteger(&curStack, 0);
  }
  else
  {
    Dec2Bin(textFactor, curStack.limbs, (int)(ptrTextFactor - &textFactor[0]),
      &curStack.nbrLimbs);
  }
  return EXPR_OK;
}
#endif

static int ComputeSumDigits(void)
{
  getCurrentStackValue(&curStack);    // Get first argument.
  stackIndex++;
  getCurrentStackValue(&curStack2);   // Get second argument.
  stackIndex--;
  BigInteger argum, Temp;
  BigInteger *result = &curStack;
  BigInteger *radix = &curStack2;
  CopyBigInt(&argum, &curStack);
  intToBigInteger(result, 0);
  while ((argum.nbrLimbs > 1) || (argum.limbs[0].x > 0))
  {
    BigIntRemainder(&argum, radix, &Temp);
    BigIntAdd(result, &Temp, result);
    BigIntDivide(&argum, radix, &argum);
  }
  return EXPR_OK;
}

static int ComputeNumDigits(void)
{
  getCurrentStackValue(&curStack);    // Get first argument.
  stackIndex++;
  getCurrentStackValue(&curStack2);   // Get second argument.
  stackIndex--;
  BigInteger *result = &curStack;
  BigInteger *radix = &curStack2;
  int digits = 0;
  while (!BigIntIsZero(result))
  {
    BigIntDivide(result, radix, result);
    digits++;
  }
  intToBigInteger(result, digits);
  return EXPR_OK;
}

static int ComputeRevDigits(void)
{
  getCurrentStackValue(&curStack);    // Get first argument.
  stackIndex++;
  getCurrentStackValue(&curStack2);   // Get second argument.
  stackIndex--;
  BigInteger argum;
  BigInteger Temp;
  BigInteger *result = &curStack;
  BigInteger *radix = &curStack2;
  CopyBigInt(&argum, &curStack);
  intToBigInteger(result, 0);
  while (!BigIntIsZero(&argum))
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
  int prevLimb;
  int curLimb;
  int *ptrDest;
  int *ptrSrc;
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
#ifdef FACTORIZATION_APP
    if ((unsigned int)(first->nbrLimbs * BITS_PER_GROUP + shiftCtr) > 664380)
#else
    if ((unsigned int)(first->nbrLimbs * BITS_PER_GROUP + shiftCtr) > 66438)
#endif
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
      (void)memset(first->limbs, 0, delta * sizeof(limb));
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
    if ((second->nbrLimbs > 1) || (shiftCtr > first->nbrLimbs * BITS_PER_GROUP))
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
    if ((result->nbrLimbs == 0) || (result->limbs[result->nbrLimbs].x))
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
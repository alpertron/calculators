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
int lang;
char output[3000000];
limb Mult1[MAX_LEN];
limb Mult3[MAX_LEN];
limb Mult4[MAX_LEN];
int q[MAX_LEN];
#define fibon1 MontgomeryR1
static enum eExprErr ComputeSubExpr(int stackIndex);
// Operators accepted: +, -, *, /, ^, !, F(, L(, P(, B(. N(.

enum eExprErr ComputeExpression(char *expr, int type, BigInteger *ExpressionResult)
{
  int stackIndex = 0;
  int exprIndex = 0;
  int exprLength = (int)strlen(expr);
  int i, j, shLeft;
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
    if (charValue == 'B' || charValue == 'b' ||
        charValue == 'N' || charValue == 'n' ||
        charValue == 'F' || charValue == 'f' ||
        charValue == 'P' || charValue == 'p' ||
        charValue == 'L' || charValue == 'l')
    {
      if (leftNumberFlag != 0 || exprIndex == exprLength - 1)
      {
        return EXPR_SYNTAX_ERROR;
      }
      exprIndex++;
      if (*(expr + exprIndex) != '(')
      {
        return EXPR_SYNTAX_ERROR;
      }
      if (stackIndex > PAREN_STACK_SIZE-5)
      {
        return EXPR_TOO_MANY_PAREN;
      }
      stackOperators[stackIndex++] = charValue & 0xDF; /* Convert to uppercase */
      charValue = '(';
    }
    if (charValue == '+' || charValue == '-')
    {
      if (leftNumberFlag == FALSE)
      {      // Unary plus/minus operator
        exprIndex++;
        if (charValue == '+')
        {
          continue;
        }
        else
        {
          if (stackIndex > 0 && stackOperators[stackIndex - 1] == '_')
          {
            stackIndex--;
            continue;
          }
          if (stackIndex > PAREN_STACK_SIZE - 5)
          {
            return EXPR_TOO_MANY_PAREN;
          }
          stackOperators[stackIndex++] = (int)'_'; /* Unitary minus */
          continue;
        }
      }
      if (stackIndex > 0 && stackOperators[stackIndex - 1] != '(')
      {
        if ((SubExprResult = ComputeSubExpr(--stackIndex)) != 0)
        {
          return SubExprResult;
        }
        if (stackIndex > 0 && stackOperators[stackIndex - 1] != '(')
        {
          if ((SubExprResult = ComputeSubExpr(--stackIndex)) != 0)
          {
            return SubExprResult;
          }
          if (stackIndex > 0 && stackOperators[stackIndex - 1] != '(')
          {
            if ((SubExprResult = ComputeSubExpr(--stackIndex)) != 0)
            {
              return SubExprResult;
            }
          }                         /* end if */
        }                           /* end if */
      }                             /* end if */
      stackOperators[stackIndex++] = charValue;
      leftNumberFlag = FALSE;
    }                               /* end if */
    else if (charValue == '*' || charValue == '/' || charValue == '%')
    {
      if (leftNumberFlag == FALSE)
      {
        return EXPR_SYNTAX_ERROR;
      }
      if (stackIndex > 0 && (stackOperators[stackIndex - 1] == '^' ||
          stackOperators[stackIndex - 1] == '*' ||
          stackOperators[stackIndex - 1] == '/' ||
          stackOperators[stackIndex - 1] == '%' ||
          stackOperators[stackIndex - 1] == 'B' ||
          stackOperators[stackIndex - 1] == 'N' ||
          stackOperators[stackIndex - 1] == 'F' ||
          stackOperators[stackIndex - 1] == 'L' ||
          stackOperators[stackIndex - 1] == 'P'))
      {
        if ((SubExprResult = ComputeSubExpr(--stackIndex)) != 0)
        {
          return SubExprResult;
        }
        if (stackIndex > 0 && (stackOperators[stackIndex - 1] == '^' ||
            stackOperators[stackIndex - 1] == '*' ||
            stackOperators[stackIndex - 1] == '/' ||
            stackOperators[stackIndex - 1] == '%' ||
            stackOperators[stackIndex - 1] == 'B' ||
            stackOperators[stackIndex - 1] == 'N' ||
            stackOperators[stackIndex - 1] == 'F' ||
            stackOperators[stackIndex - 1] == 'L' ||
            stackOperators[stackIndex - 1] == 'P'))
        {
          if ((SubExprResult = ComputeSubExpr(--stackIndex)) != 0)
          {
            return SubExprResult;
          }
        }                         /* end if */
      }                           /* end if */
      stackOperators[stackIndex++] = charValue;
      leftNumberFlag = FALSE;
    }
    else if (charValue == '^')
    {
      if (leftNumberFlag == FALSE)
      {
        return EXPR_SYNTAX_ERROR;
      }
      if (stackIndex > 0 && (stackOperators[stackIndex - 1] == '^' ||
          stackOperators[stackIndex - 1] == 'B' ||
          stackOperators[stackIndex - 1] == 'N' ||
          stackOperators[stackIndex - 1] == 'F' ||
          stackOperators[stackIndex - 1] == 'L' ||
          stackOperators[stackIndex - 1] == 'P'))
      {
        if ((SubExprResult = ComputeSubExpr(--stackIndex)) != 0)
        {
          return SubExprResult;
        }
      }                         /* end if */
      stackOperators[stackIndex++] = charValue;
      leftNumberFlag = FALSE;
    }                           /* end if */
    else if (charValue == '(')
    {
      if (leftNumberFlag == TRUE)
      {
        return EXPR_SYNTAX_ERROR;
      }
      if (stackIndex > PAREN_STACK_SIZE - 5)
      {
        return EXPR_TOO_MANY_PAREN;
      }
      stackOperators[stackIndex++] = charValue;
    }
    else if (charValue == ')')
    {
      if (leftNumberFlag == FALSE)
      {
        return EXPR_SYNTAX_ERROR;
      }
      if (stackIndex > 0 && stackOperators[stackIndex - 1] != '(')
      {
        if ((SubExprResult = ComputeSubExpr(--stackIndex)) != 0)
        {
          return SubExprResult;
        }
        if (stackIndex > 0 && stackOperators[stackIndex - 1] != '(')
        {
          if ((SubExprResult = ComputeSubExpr(--stackIndex)) != 0)
          {
            return SubExprResult;
          }
          if (stackIndex > 0 && stackOperators[stackIndex - 1] != '(')
          {
            if ((SubExprResult = ComputeSubExpr(--stackIndex)) != 0)
            {
              return SubExprResult;
            }
          }
        }
      }
      if (stackIndex == 0)
      {
        return EXPR_PAREN_MISMATCH;
      }
      stackIndex--;             /* Discard ")" */
      stackValues[stackIndex] = stackValues[stackIndex + 1];
      leftNumberFlag = TRUE;
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
            c -= 'A'+10;
          }
          else
          {
            c -= 'a' + 10;
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
        (ptrLimb++)->x = carry.x;
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
  if (stackIndex > 0 && stackOperators[stackIndex - 1] != '(')
  {
    if ((SubExprResult = ComputeSubExpr(--stackIndex)) != 0)
    {
      return SubExprResult;
    }
    if (stackIndex > 0 && stackOperators[stackIndex - 1] != '(')
    {
      if ((SubExprResult = ComputeSubExpr(--stackIndex)) != 0)
      {
        return SubExprResult;
      }
      if (stackIndex > 0 && stackOperators[stackIndex - 1] != '(')
      {
        if ((SubExprResult = ComputeSubExpr(--stackIndex)) != 0)
        {
          return SubExprResult;
        }
      }
    }
  }
  if (stackIndex != 0)
  {
    return EXPR_PAREN_MISMATCH;
  }
  if (type == 0)
  {
    if (stackValues[0].nbrLimbs == 1 && stackValues[0].limbs[0].x < 1)
    {
      return EXPR_NUMBER_TOO_LOW;
    }
    if (stackValues[0].sign == SIGN_NEGATIVE)
    {
      return EXPR_NUMBER_TOO_LOW;
    }
  }
  CopyBigInt(ExpressionResult, &stackValues[0]);
  return EXPR_OK;
}

static enum eExprErr ComputeSubExpr(int stackIndex)
{
  int len, val, ctr;
  BigInteger *pArgument, *pResult;
  limb *pResultLimbs, *pArgumentLimbs, *pTemp, *pFibonPrev, *pFibonAct;
  limb largeVal;
  char stackOper;
  int nbrLimbs;
  int FirstTime;

  stackOper = stackOperators[stackIndex];
  switch (stackOper)
  {
  case '+':
    BigIntAdd(&stackValues[stackIndex], &stackValues[stackIndex + 1], &stackValues[stackIndex]);
    return EXPR_OK;
  case '-':
    BigIntSubt(&stackValues[stackIndex], &stackValues[stackIndex + 1], &stackValues[stackIndex]);
    return EXPR_OK;
  case '_':
    BigIntNegate(&stackValues[stackIndex+1], &stackValues[stackIndex]);
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
  case 'F':
  case 'L':
    pArgument = &stackValues[stackIndex + 1];
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
      fibon2[0].x = (stackOper == 'L'? 2: 0);
    }
    else
    {
      int i, j;
      // For Lucas sequences: FibonPrev = 2, FibonAct = 1
      // For Fibonacci sequences: FibonPrev = 0, FibonAct = 1
      fibon1[0].x = (stackOper == 'L' ? 2 : 0);
      fibon2[0].x = 1;
      for (i = 1; i < val; i++)
      {
        unsigned int carry = 0;
        for (j=0; j<len; j++)
        {
          carry += (pFibonPrev+j)->x + (pFibonAct+j)->x;
          (pFibonPrev+j)->x = carry & MAX_VALUE_LIMB;
          carry >>= BITS_PER_GROUP;
        }
        if (carry != 0)
        {
          (pFibonPrev+j)->x = carry;
          (pFibonAct+j)->x = 0;
          len++;
        }
        pTemp = pFibonAct;
        pFibonAct = pFibonPrev;
        pFibonPrev = pTemp;
      }
    }
    pArgument = &stackValues[stackIndex];
    pArgument -> sign = SIGN_POSITIVE;
    pArgument -> nbrLimbs = len;
    memcpy(pArgument->limbs, pFibonAct, len*sizeof(limb));
    return EXPR_OK;

  case 'P':
    pArgument = &stackValues[stackIndex + 1];
    if (pArgument->sign == SIGN_NEGATIVE)
    {
      return EXPR_INVALID_PARAM;
    }
    if (pArgument->nbrLimbs > 2)
    {
      return EXPR_INTERM_TOO_HIGH;
    }
    if (pArgument->nbrLimbs==2)
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
    break;
  case 'B':
  case 'N':
    pArgument = &stackValues[stackIndex + 1];
    pResult = &stackValues[stackIndex];
    pResultLimbs = pResult -> limbs;
    pArgumentLimbs = pArgument -> limbs;
    nbrLimbs = pArgument->nbrLimbs;
    FirstTime = 1;
    if (stackOper == 'B')
    {        // Previous pseudoprime
      if (pArgument->sign == SIGN_NEGATIVE || (nbrLimbs == 1 && pArgumentLimbs->x<3))
      {
        return EXPR_INVALID_PARAM;
      }
      pResult -> sign = SIGN_POSITIVE;
      if (nbrLimbs == 1 && pArgumentLimbs->x == 3)
      {
        pResult -> nbrLimbs = 1;
        pResultLimbs->x = 2;
        return EXPR_OK;
      }
      memcpy(pResultLimbs, pArgumentLimbs, nbrLimbs*sizeof(limb));
      pResultLimbs->x = pResultLimbs->x | 1;
    }
    else
    {        // Next pseudoprime
      pResult -> sign = SIGN_POSITIVE;
      if (pArgument->sign == SIGN_NEGATIVE || (nbrLimbs == 1 && pArgumentLimbs->x < 2))
      {
        pResult -> nbrLimbs = 1;
        pResultLimbs->x = 2;
        return EXPR_OK;
      }
      memcpy(pResultLimbs, pArgumentLimbs, nbrLimbs*sizeof(limb));
    }
    for (;;)
    {        // Search for next pseudoprime.
      if (stackOper == 'B')
      {      // Subtract 2
        pResultLimbs->x = (pResultLimbs->x - 2);
        if (pResultLimbs->x < 0)
        {
          pTemp = pResultLimbs;
          for (ctr = 0; ctr < nbrLimbs; ctr++)
          {
            (pTemp++)->x = MAX_VALUE_LIMB;
            if (--(pTemp->x) >= 0)
            {
              break;
            }
          }
        }
      }
      else
      {      // Next pseudoprime
        if (FirstTime != 0)
        {
          if (nbrLimbs == 1 && pResultLimbs->x == 2)
          {   // Find next pseudoprime to 2.
            pResultLimbs->x = 1;
          }
          else
          {
            pResultLimbs->x |= 1;  // If it is even, change to next odd value.
          }
          FirstTime = 0;
          pResult->nbrLimbs = nbrLimbs;
        }
        pResultLimbs->x += 2;   // Add 2.
        if (pResultLimbs->x > MAX_VALUE_LIMB)
        {
          pTemp = pResultLimbs;
          for (ctr = 0; ctr < nbrLimbs; ctr++)
          {
            (pTemp++)->x -= MAX_VALUE_LIMB + 1;
            if (++(pTemp->x) <= MAX_VALUE_LIMB)
            {
              break;
            }
          }
          if (ctr == nbrLimbs)
          {
            (pTemp++)->x = 0;
            pTemp->x = 1;
            nbrLimbs++;
            pResult->nbrLimbs = nbrLimbs;
          }
        }
      }
      if (BpswPrimalityTest(pResult) == 0)
      {   // Number is probable prime.
        pResult->nbrLimbs = nbrLimbs;
        pResult->sign = SIGN_POSITIVE;
        return EXPR_OK;
      }
    }
  }              /* end switch */
  return EXPR_OK;
}

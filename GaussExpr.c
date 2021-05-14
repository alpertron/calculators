//
// This file is part of Alpertron Calculators.
//
// Copyright 2016-2021 Dario Alejandro Alpern
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
#include "factor.h"
#include "expression.h"

#define PAREN_STACK_SIZE 100

  // Errors for next routine:
  // >0: Ok.
  // -2: Number too high (more than 10000 digits).
  // EXPR_INTERM_TOO_HIGH: Intermediate expression too high (more than 20000 digits).
  // EXPR_PAREN_MISMATCH: Parenthesis mismatch.
  // EXPR_SYNTAX_ERROR: Syntax error
  // EXPR_TOO_MANY_PAREN: Too many parentheses.
  // EXPR_INVALID_PARAM: Invalid parameter.
  // -100: Break.
  // Operators accepted: +, -, *, /, ^, !, F(, L(, P(.

BigInteger stackRealValues[PAREN_STACK_SIZE];
BigInteger stackImagValues[PAREN_STACK_SIZE];
int stackOperators[PAREN_STACK_SIZE];
static int stackIndex;
static int exprIndex;
#ifndef lang  
  bool lang;
#endif
char output[3000000];
limb Mult1[MAX_LEN];
limb Mult3[MAX_LEN];
limb Mult4[MAX_LEN];
int q[MAX_LEN];
static int exprLength;
static int ComputeExpr(char *expr, BigInteger *ExpressionResult);
static int ComputeSubExpr(void);
static int ComputeGCD(void);
static bool func(char* expr, BigInteger* ExpressionResult,
  const char* funcName, int funcArgs, bool leftNumberFlag, int* pResult);
static void SkipSpaces(const char *expr);
static int Modulo(BigInteger* ReNum, BigInteger* ImNum,
  const BigInteger* ReDen, const BigInteger* ImDen,
  BigInteger* Result);
static int ComputeFibonacci(void);
static int ComputeLucas(void);
static int ComputeModPow(void);
static int ComputeModInv(void);
static int ComputePartition(void);
static int ComputePower(const BigInteger* Re1, const BigInteger* Re2,
  const BigInteger* Im1, const BigInteger* Im2);
static int ModInv(const BigInteger* RealNbr, const BigInteger* ImagNbr,
  BigInteger* RealMod, BigInteger* ImagMod,
  BigInteger* Result);

enum eExprErr ComputeGaussianExpression(char *expr, BigInteger *ExpressionResult)
{
  int retcode;
  stackIndex = 0;
  exprIndex = 0;
  retcode = ComputeExpr(expr, ExpressionResult);
  if (retcode != 0) {return retcode;}
  if ((ExpressionResult[0].nbrLimbs > 2215) &&    // 10000/log_10(32768)
      (ExpressionResult[1].nbrLimbs > 2215))
  {
    return EXPR_NUMBER_TOO_HIGH;
  }
  return 0;
}

static int ComputeExpr(char *expr, BigInteger *ExpressionResult)
{
  int c;
  int i;
  int j;
  int exprIndexAux;
  int SubExprResult;
  int len;
  static BigInteger factorialResult;
  static BigInteger Tmp;
  static BigInteger *ptrBigInt;
  const static BigInteger *ptrRe;
  const static BigInteger *ptrIm;
  limb carry;
  limb largeLen;
  limb *ptrLimb;
  int retcode;
  int shLeft;
  int offset;
  bool leftNumberFlag = false;
  int startStackIndex = stackIndex;
  
  exprLength = (int)strlen(expr);
  while (exprIndex < exprLength)
  {
    char charValue = *(expr + exprIndex);
    if (charValue == '!')
    {           // Calculating factorial.
      if (leftNumberFlag == false)
      {
        return EXPR_SYNTAX_ERROR;
      }
      if (!BigIntIsZero(&stackImagValues[stackIndex]))
      {         // Imaginary part must be zero.
        return EXPR_INVALID_PARAM;
      }
      if (stackRealValues[stackIndex].nbrLimbs > 1)
      {
        return EXPR_INTERM_TOO_HIGH;
      }
      if ((stackRealValues[stackIndex].limbs[0].x < 0) ||
          (stackRealValues[stackIndex].limbs[0].x > 5984))
      {
        return EXPR_INTERM_TOO_HIGH;
      }
      len = stackRealValues[stackIndex].limbs[0].x;
      factorialResult.limbs[0].x = 1;
      factorialResult.nbrLimbs = 1;
      factorialResult.sign = SIGN_POSITIVE;
      for (i = 2; i <= len; i++)
      {   // Multiply by all integers up to the argument of factorial.
        carry.x = 0;
        for (j = 0; j < factorialResult.nbrLimbs; j++)
        {
          carry.x += i*factorialResult.limbs[j].x;
          factorialResult.limbs[j].x = carry.x & MAX_VALUE_LIMB;
          carry.x >>= BITS_PER_GROUP;
        }
        if (carry.x != 0)
        {  // New limb needed.
          factorialResult.limbs[j].x = carry.x;
          factorialResult.nbrLimbs++;
        }
      }
      CopyBigInt(&stackRealValues[stackIndex], &factorialResult);
    }
    else if (charValue == '#')
    {           // Calculating primorial.
      if (leftNumberFlag == false)
      {
        return EXPR_SYNTAX_ERROR;
      }
      if (!BigIntIsZero(&stackImagValues[stackIndex]))
      {         // Imaginary part must be zero.
        return EXPR_INVALID_PARAM;
      }
      if (stackRealValues[stackIndex].nbrLimbs > 2)
      {
        return EXPR_INTERM_TOO_HIGH;
      }
      if (stackRealValues[stackIndex].nbrLimbs == 2)
      {
        largeLen.x = stackRealValues[stackIndex].limbs[0].x +
          (stackRealValues[stackIndex].limbs[1].x << BITS_PER_GROUP);
      }
      else
      {
        largeLen.x = stackRealValues[stackIndex].limbs[0].x;
      }
      if ((largeLen.x < 0) || (largeLen.x > 46049))
      {
        return EXPR_INTERM_TOO_HIGH;
      }
      len = largeLen.x;
      factorialResult.limbs[0].x = 1;
      factorialResult.nbrLimbs = 1;
      factorialResult.sign = SIGN_POSITIVE;
      for (i = 2; i <= len; i++)
      {      // Multiply by prime numbers only.
        for (j = 2; (j*j) <= i; j++)
        {
          if ((i / j*j) == i)
          {   // Number is not prime.
            break;
          }
        }
        if ((j*j) > i)
        {     // Number is prime, perform multiplication.
          carry.x = 0;
          for (j = 0; j < factorialResult.nbrLimbs; j++)
          {
            carry.x += i*factorialResult.limbs[j].x;
            factorialResult.limbs[j].x = carry.x & MAX_VALUE_LIMB;
            carry.x >>= BITS_PER_GROUP;
          }
          if (carry.x != 0)
          {  // New limb needed.
            factorialResult.limbs[j].x = carry.x;
            factorialResult.nbrLimbs++;
          }
        }
      }
      CopyBigInt(&stackRealValues[stackIndex], &factorialResult);
    }
    else if (func(expr, ExpressionResult, "GCD", 2, leftNumberFlag, &retcode))
    {
      if (retcode != 0) {return retcode;}
      retcode = ComputeGCD();
      if (retcode != 0) {return retcode;}
      leftNumberFlag = true;
    }
    else if (func(expr, ExpressionResult, "RE", 1, leftNumberFlag, &retcode))
    {
      ptrBigInt = &stackImagValues[stackIndex];
      ptrBigInt->limbs[0].x = 0;
      ptrBigInt->nbrLimbs = 1;
      ptrBigInt->sign = SIGN_POSITIVE;
      if (retcode != 0)
      {
        return retcode;
      }
      leftNumberFlag = true;
    }
    else if (func(expr, ExpressionResult, "NORM", 1, leftNumberFlag, &retcode))
    {
      if (retcode != 0) {return retcode;}
      ptrRe = &stackRealValues[stackIndex];
      ptrIm = &stackImagValues[stackIndex];
      (void)BigIntMultiply(ptrRe, ptrRe, &factorialResult);    // norm <- re2^2 + im2^2.
      (void)BigIntMultiply(ptrIm, ptrIm, &Tmp);
      BigIntAdd(&factorialResult, &Tmp, &stackRealValues[stackIndex]);  
      ptrBigInt = &stackImagValues[stackIndex];
      ptrBigInt->limbs[0].x = 0;
      ptrBigInt->nbrLimbs = 1;
      ptrBigInt->sign = SIGN_POSITIVE;
      leftNumberFlag = true;
    }
    else if (func(expr, ExpressionResult, "IM", 1, leftNumberFlag, &retcode))
    {
      if (retcode != 0) {return retcode;}
      stackRealValues[stackIndex] = stackImagValues[stackIndex];
      ptrBigInt = &stackImagValues[stackIndex];
      ptrBigInt->limbs[0].x = 0;
      ptrBigInt->nbrLimbs = 1;
      ptrBigInt->sign = SIGN_POSITIVE;
      leftNumberFlag = true;
    }
    else if (func(expr, ExpressionResult, "MODPOW", 3, leftNumberFlag, &retcode))
    {
      if (retcode != 0) {return retcode;}
      retcode = ComputeModPow();
      if (retcode != 0) {return retcode;}
      leftNumberFlag = true;
    }
    else if (func(expr, ExpressionResult, "MODINV", 2, leftNumberFlag, &retcode))
    {
      if (retcode != 0) {return retcode;}
      retcode = ComputeModInv();
      if (retcode != 0) {return retcode;}
      leftNumberFlag = true;
    }
    else if (func(expr, ExpressionResult, "F", 1, leftNumberFlag, &retcode))
    {
      if (retcode != 0) {return retcode;}
      retcode = ComputeFibonacci();
      if (retcode != 0) {return retcode;}
      leftNumberFlag = true;
    }
    else if (func(expr, ExpressionResult, "L", 1, leftNumberFlag, &retcode))
    {
      if (retcode != 0) {return retcode;}
      retcode = ComputeLucas();
      if (retcode != 0) {return retcode;}
      leftNumberFlag = true;
    }
    else if (func(expr, ExpressionResult, "P", 1, leftNumberFlag, &retcode))
    {
      if (retcode != 0) {return retcode;}
      retcode = ComputePartition();
      if (retcode != 0) {return retcode;}
      leftNumberFlag = true;
    }
    else if ((charValue == '+') || (charValue == '-'))
    {
      if (!leftNumberFlag)
      {      // Unary plus/minus operator
        exprIndex++;
        if (charValue == '+')
        {
          continue;
        }
        else
        {
          if ((stackIndex > startStackIndex) && (stackOperators[stackIndex-1] == '_'))
          {
            stackIndex--;
            continue;
          }
          if (stackIndex >= PAREN_STACK_SIZE)
          {
            return EXPR_TOO_MANY_PAREN;
          }
          stackOperators[stackIndex] = '_'; /* Unary minus */
          stackIndex++;
          continue;
        }
      }
      if ((stackIndex > startStackIndex) &&
          (stackOperators[stackIndex-1] != '('))
      {
        SubExprResult = ComputeSubExpr();
        if (SubExprResult != 0)
        {
          return SubExprResult;
        }
        if ((stackIndex > startStackIndex) &&
            (stackOperators[stackIndex-1] != '('))
        {
          SubExprResult = ComputeSubExpr();
          if (SubExprResult != 0)
          {
            return SubExprResult;
          }
          if ((stackIndex > startStackIndex) &&
              (stackOperators[stackIndex-1] != '('))
          {
            SubExprResult = ComputeSubExpr();
            if (SubExprResult != 0)
            {
              return SubExprResult;
            }
          }                         /* end if */
        }                           /* end if */
      }                             /* end if */
      stackOperators[stackIndex] = charValue;
      stackIndex++;
      leftNumberFlag = false;
    }                               /* end if */
    else if ((charValue == '*') || (charValue == '/') || (charValue == '%'))
    {
      if (leftNumberFlag == 0)
      {
        return EXPR_SYNTAX_ERROR;
      }
      if ((stackIndex > startStackIndex) && ((stackOperators[stackIndex-1] == '^') ||
          (stackOperators[stackIndex-1] == '*') ||
          (stackOperators[stackIndex-1] == '/')))
      {
        SubExprResult = ComputeSubExpr();
        if (SubExprResult != 0)
        {
          return SubExprResult;
        }
        if ((stackIndex > startStackIndex) &&
             ((stackOperators[stackIndex-1] == '^') ||
              (stackOperators[stackIndex-1] == '*') ||
              (stackOperators[stackIndex-1] == '/') ||
              (stackOperators[stackIndex-1] == '%')))
        {
          SubExprResult = ComputeSubExpr();
          if (SubExprResult != 0)
          {
            return SubExprResult;
          }
        }                         /* end if */
      }                           /* end if */
      stackOperators[stackIndex] = charValue;
      stackIndex++;
      leftNumberFlag = false;
    }                             
    else if (charValue == '^')
    {
      if (leftNumberFlag == 0)
      {
        return EXPR_SYNTAX_ERROR;
      }
      stackOperators[stackIndex] = charValue;
      stackIndex++;
      leftNumberFlag = false;
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
      stackOperators[stackIndex] = charValue;
      stackIndex++;
    }                           
    else if ((charValue == ')') || (charValue == ','))
    {
      if (leftNumberFlag == 0)
      {
        return EXPR_SYNTAX_ERROR;
      }
      if ((stackIndex > startStackIndex) &&
          (stackOperators[stackIndex-1] != '('))
      {
        if ((SubExprResult = ComputeSubExpr()) != 0)
        {
          return SubExprResult;
        }
        if ((stackIndex > startStackIndex) &&
            (stackOperators[stackIndex-1] != '('))
        {
          if ((SubExprResult = ComputeSubExpr()) != 0)
          {
            return SubExprResult;
          }
          if ((stackIndex > startStackIndex) &&
              (stackOperators[stackIndex-1] != '('))
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
      stackRealValues[stackIndex] = stackRealValues[stackIndex+1];
      stackImagValues[stackIndex] = stackImagValues[stackIndex+1];
      leftNumberFlag = true;
    }
    else if ((charValue == 'i') || (charValue == 'I'))
    {
      if (leftNumberFlag == 0)
      {
        // Initialize real part to zero.
        ptrBigInt = &stackRealValues[stackIndex];
        ptrBigInt->limbs[0].x = 0;
        ptrBigInt->nbrLimbs = 1;
        ptrBigInt->sign = SIGN_POSITIVE;
        // Initialize imaginary part to zero.
        ptrBigInt = &stackImagValues[stackIndex];
        ptrBigInt->limbs[0].x = 1;
        ptrBigInt->nbrLimbs = 1;
        ptrBigInt->sign = SIGN_POSITIVE;

        leftNumberFlag = true;
      }
      else
      {     // Multiply by i. First exchange imaginary and real coefficients.
        ptrBigInt = &stackRealValues[stackIndex];
        CopyBigInt(&Tmp, &stackImagValues[stackIndex]);
        CopyBigInt(&stackImagValues[stackIndex], ptrBigInt);
        CopyBigInt(ptrBigInt, &Tmp);
        if (!BigIntIsZero(ptrBigInt))
        {            // Real coefficient is not zero: change sign
          ptrBigInt->sign = ((ptrBigInt->sign == SIGN_POSITIVE)? SIGN_NEGATIVE : SIGN_POSITIVE);
        }
      }
    }
    else if ((charValue >= '0') && (charValue <= '9'))
    {
      exprIndexAux = exprIndex;
      if ((charValue == '0') && (exprIndexAux < (exprLength - 2)) &&
        (*(expr + exprIndexAux + 1) == 'x'))
      {  // hexadecimal
        exprIndexAux += 2;
        while (exprIndexAux < (exprLength - 1))
        {
          charValue = *(expr + exprIndexAux + 1);
          if (((charValue >= '0') && (charValue <= '9')) ||
            ((charValue >= 'A') && (charValue <= 'F')) ||
            ((charValue >= 'a') && (charValue <= 'f')))
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
        ptrLimb = &stackRealValues[stackIndex].limbs[0];
        for (; exprIndexAux >= (exprIndex + 2); exprIndexAux--)
        {
          c = *(expr + exprIndexAux);
          if ((c >= '0') && (c <= '9'))
          {
            c -= '0';
          }
          else if ((c >= 'A') && (c <= 'F'))
          {
            c -= 'A' + 10;
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
            ptrLimb->x = carry.x & MAX_VALUE_LIMB;
            ptrLimb++;
            carry.x >>= BITS_PER_GROUP;
          }
        }
        ptrLimb->x = carry.x;
        ptrLimb++;
        exprIndex = offset;
        stackRealValues[stackIndex].nbrLimbs = (int)(ptrLimb - &stackRealValues[stackIndex].limbs[0]);
        stackRealValues[stackIndex].sign = SIGN_POSITIVE;
      }
      else
      {                   // Decimal number.
        while (exprIndexAux < (exprLength - 1))
        {
          charValue = *(expr + exprIndexAux + 1);
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
        ptrBigInt = &stackRealValues[stackIndex];
        Dec2Bin(expr + exprIndex, ptrBigInt->limbs,
          exprIndexAux + 1 - exprIndex, &ptrBigInt->nbrLimbs);
        ptrBigInt->sign = SIGN_POSITIVE;
        exprIndex = exprIndexAux;
      }
      ptrBigInt = &stackImagValues[stackIndex];
      ptrBigInt->limbs[0].x = 0;     // Initialize imaginary part to zero.
      ptrBigInt->nbrLimbs = 1;
      ptrBigInt->sign = SIGN_POSITIVE;
      leftNumberFlag = true;
    }                            /* end if */
    exprIndex++;
  }                              /* end while */
  if (leftNumberFlag == 0)
  {
    return EXPR_SYNTAX_ERROR;
  }
  if ((stackIndex > startStackIndex) && (stackOperators[stackIndex-1] != '('))
  {
    if ((SubExprResult = ComputeSubExpr()) != 0)
    {
      return SubExprResult;
    }
    if ((stackIndex > startStackIndex) && (stackOperators[stackIndex-1] != '('))
    {
      if ((SubExprResult = ComputeSubExpr()) != 0)
      {
        return SubExprResult;
      }
      if ((stackIndex > startStackIndex) && (stackOperators[stackIndex-1] != '('))
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
  CopyBigInt(&ExpressionResult[0], &stackRealValues[startStackIndex]);
  CopyBigInt(&ExpressionResult[1], &stackImagValues[startStackIndex]);
  return 0;
}

static bool func(char *expr, BigInteger *ExpressionResult,
                const char *funcName, int funcArgs, bool leftNumberFlag, int *pResult)
{
  int funcNameLen = (int)strlen(funcName);
  const char *ptrExpr;
  const char *ptrFuncName;

  if ((exprIndex + funcNameLen) > exprLength)
  {
    return false;           // Next characters are not a function name.
  }
  ptrExpr = expr + exprIndex;
  ptrFuncName = funcName;
  while (*ptrFuncName != '\0')
  {
    if ((*ptrExpr & 0xDF) != *ptrFuncName)
    {
      return false;         // Next characters are not a function name.
    }
    ptrExpr++;
    ptrFuncName++;
  }
  exprIndex += funcNameLen;
  if (leftNumberFlag)
  {
    *pResult = EXPR_SYNTAX_ERROR;
    return true;
  }
  SkipSpaces(expr);
  if ((exprIndex == exprLength) || (*(expr+exprIndex) != '('))
  {
    *pResult = EXPR_SYNTAX_ERROR;
    return true;
  }
  exprIndex++;
  for (int index = 0; index < funcArgs; index++)
  {
    int retcode;
    char compareChar;

    SkipSpaces(expr);
    if (stackIndex >= PAREN_STACK_SIZE)
    {
      *pResult = EXPR_TOO_MANY_PAREN;
      return true;
    }
    retcode = ComputeExpr(expr, ExpressionResult);
    if (retcode != 0) {return retcode;}
    SkipSpaces(expr);
    compareChar = ((index == (funcArgs-1))? ')': ',');
    if ((exprIndex == exprLength) || (*(expr+exprIndex) != compareChar))
    {
      *pResult = EXPR_SYNTAX_ERROR;
      return true;
    }
    exprIndex++;
    stackIndex++;
  }
  stackIndex -= funcArgs;
  exprIndex--;
  *pResult = EXPR_OK;
  return true;
}

static void SkipSpaces(const char *expr)
{
  while (*(expr+exprIndex) != '\0')
  {
    if (*(expr + exprIndex) > ' ')
    {
      break;
    }
    exprIndex++;
  }
  return;
}

static int ComputeSubExpr(void)
{
  int retcode;
  int stackOper;
  static BigInteger Re1;
  static BigInteger Re2;
  static BigInteger Im1;
  static BigInteger Im2;
  static BigInteger Re;
  static BigInteger Im;
  static BigInteger ReTmp;
  static BigInteger ImTmp;
  static BigInteger Result[2];
  static BigInteger norm;
  stackIndex--;

  stackOper = stackOperators[stackIndex];
  Re1 = stackRealValues[stackIndex];
  Re2 = stackRealValues[stackIndex+1];
  Im1 = stackImagValues[stackIndex];
  Im2 = stackImagValues[stackIndex+1];

  switch (stackOper) {
    case '+':
      BigIntAdd(&Re1, &Re2, &stackRealValues[stackIndex]);
      BigIntAdd(&Im1, &Im2, &stackImagValues[stackIndex]);
      return 0;
    case '-':
      BigIntSubt(&Re1, &Re2, &stackRealValues[stackIndex]);
      BigIntSubt(&Im1, &Im2, &stackImagValues[stackIndex]);
      return 0;
    case '_':
      BigIntNegate(&Re2, &stackRealValues[stackIndex]);
      BigIntNegate(&Im2, &stackImagValues[stackIndex]);
      return 0;
    case '/':
      (void)BigIntMultiply(&Re2, &Re2, &ReTmp);
      (void)BigIntMultiply(&Im2, &Im2, &ImTmp);
      BigIntAdd(&ReTmp, &ImTmp, &norm);       // norm <- re2^2 + im2^2.
      if (BigIntIsZero(&norm))
      {              // norm is zero.
        return EXPR_INTERM_TOO_HIGH;
      }
      (void)BigIntMultiply(&Re1, &Re2, &ReTmp);
      (void)BigIntMultiply(&Im1, &Im2, &ImTmp);
      BigIntAdd(&ReTmp, &ImTmp, &Re);         // Re <- re1*re2 + im1*im2.
      (void)BigIntMultiply(&Im1, &Re2, &ReTmp);
      (void)BigIntMultiply(&Re1, &Im2, &ImTmp);
      BigIntSubt(&ReTmp, &ImTmp, &Im);        // Im <- im1*re2 - re1*im2.
      (void)BigIntDivide(&Re, &norm, &stackRealValues[stackIndex]);
      (void)BigIntDivide(&Im, &norm, &stackImagValues[stackIndex]);
      return 0;
    case '%':
      retcode = Modulo(&Re1, &Im1, &Re2, &Im2, Result);
      if (retcode != 0)
      {
        return retcode;
      }
      CopyBigInt(&stackRealValues[stackIndex], &Result[0]);
      CopyBigInt(&stackImagValues[stackIndex], &Result[1]);
      return 0;
    case '*':
      if (((Re1.nbrLimbs + Im1.nbrLimbs) > 4430) || ((Re1.nbrLimbs + Im2.nbrLimbs) > 4430) ||
          ((Re2.nbrLimbs + Im1.nbrLimbs) > 4430) || ((Re2.nbrLimbs + Im2.nbrLimbs) > 4430))
      {           // Result with more than 20000 digits.
        return EXPR_INTERM_TOO_HIGH;
      }
      (void)BigIntMultiply(&Re1, &Re2, &ReTmp);       // Re <- re1*re2 - im1*im2.
      (void)BigIntMultiply(&Im1, &Im2, &ImTmp);
      BigIntSubt(&ReTmp, &ImTmp, &stackRealValues[stackIndex]);
      (void)BigIntMultiply(&Im1, &Re2, &ReTmp);       // Im <- im1*re2 + re1*im2.
      (void)BigIntMultiply(&Re1, &Im2, &ImTmp);
      BigIntAdd(&ReTmp, &ImTmp, &stackImagValues[stackIndex]);
      return 0;
    case '^':
      return ComputePower(&Re1, &Re2, &Im1, &Im2);
    default:
      break;
  }                /* end switch */
  return 0;
}

static int ComputeFibonacci(void)
{
  BigInteger *ptrBigInt;
  BigInteger FibonPrev;
  BigInteger FibonAct;
  BigInteger FibonNext;
  BigInteger Re = stackRealValues[stackIndex];
  int arg;
  if (!BigIntIsZero(&stackImagValues[stackIndex]))
  {     // Imaginary part of argument must be zero.
    return EXPR_INVALID_PARAM;
  }
  if (Re.nbrLimbs > 1)
  {
    return EXPR_INTERM_TOO_HIGH;
  }
  arg = (int)Re.limbs[0].x;
  if (arg > 9571)
  {
    return EXPR_INTERM_TOO_HIGH;
  }
  if (arg < 0)
  {
    return EXPR_INVALID_PARAM;
  }
  FibonPrev.limbs[0].x = 1;
  FibonAct.limbs[0].x = 0;
  FibonPrev.nbrLimbs = 1;
  FibonAct.nbrLimbs = 1;
  FibonPrev.sign = SIGN_POSITIVE; 
  FibonAct.sign = SIGN_POSITIVE;
  for (int i=1; i<=arg; i++)
  {
    BigIntAdd(&FibonPrev, &FibonAct, &FibonNext);
    CopyBigInt(&FibonPrev, &FibonAct);
    CopyBigInt(&FibonAct, &FibonNext);
  }
  CopyBigInt(&stackRealValues[stackIndex], &FibonAct);
  ptrBigInt = &stackImagValues[stackIndex];
  ptrBigInt->limbs[0].x = 0;     // Initialize imaginary part to zero.
  ptrBigInt->nbrLimbs = 1;
  ptrBigInt->sign = SIGN_POSITIVE;
  return 0;
}

static int ComputeLucas(void)
{
  int arg;
  BigInteger FibonPrev;
  BigInteger FibonAct;
  BigInteger FibonNext;
  BigInteger Re = stackRealValues[stackIndex];
  BigInteger *ptrBigInt;

  if (!BigIntIsZero(&stackImagValues[stackIndex]))
  {     // Imaginary part of argument must be zero.
    return EXPR_INVALID_PARAM;
  }
  if (Re.nbrLimbs > 1)
  {
    return EXPR_INTERM_TOO_HIGH;
  }
  arg = Re.limbs[0].x;
  if (arg > 9572)
  {
    return EXPR_INTERM_TOO_HIGH;
  }
  if (arg < 0)
  {
    return EXPR_INVALID_PARAM;
  }
  // Initialize FibonPrev <- -1, FibonAct <-2.
  FibonPrev.limbs[0].x = 1;
  FibonAct.limbs[0].x = 2;
  FibonPrev.nbrLimbs = 1;
  FibonAct.nbrLimbs = 1;
  FibonPrev.sign = SIGN_NEGATIVE;
  FibonAct.sign = SIGN_POSITIVE;
  for (int i = 1; i <= arg; i++)
  {
    BigIntAdd(&FibonPrev, &FibonAct, &FibonNext);
    CopyBigInt(&FibonPrev, &FibonAct);
    CopyBigInt(&FibonAct, &FibonNext);
  }
  CopyBigInt(&stackRealValues[stackIndex], &FibonAct);
  ptrBigInt = &stackImagValues[stackIndex];
  ptrBigInt->limbs[0].x = 0;     // Initialize imaginary part to zero.
  ptrBigInt->nbrLimbs = 1;
  ptrBigInt->sign = SIGN_POSITIVE;
  return 0;
}

static int ComputePartition(void)
{
  BigInteger *pArgument;
  limb largeVal;
  int val;

  if (!BigIntIsZero(&stackImagValues[stackIndex]))
  {         // Imaginary part must be zero.
    return EXPR_INVALID_PARAM;
  }
  pArgument = &stackRealValues[stackIndex];
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
  val = largeVal.x;
  if (val > 100000)
  {
    return EXPR_INVALID_PARAM;
  }
  partition(val, &stackRealValues[stackIndex]);
  pArgument = &stackImagValues[stackIndex];
  pArgument->limbs[0].x = 0;
  pArgument->nbrLimbs = 1;
  pArgument->sign = SIGN_POSITIVE;
  return 0;
}

// Replace dividend by the remainder.
// Back-up sign of dividend and then divide positive by positive (divisor is always positive).
static void GetRemainder(const BigInteger *norm, BigInteger *ReDividend, BigInteger *ImDividend,
  const BigInteger *ReDivisor, const BigInteger *ImDivisor,
  BigInteger *Re, BigInteger *Im)
{
  BigInteger ReTmp;
  BigInteger ImTmp;
  int signBak;
  // Re <- ((Re1*Re2+Im1*Im2)*2/norm+1)/2
  (void)BigIntMultiply(ReDividend, ReDivisor, &ReTmp);
  (void)BigIntMultiply(ImDividend, ImDivisor, &ImTmp);
  BigIntAdd(&ReTmp, &ImTmp, &ReTmp);
  multint(&ReTmp, &ReTmp, 2);
  signBak = ReTmp.sign;
  ReTmp.sign = SIGN_POSITIVE;
  (void)BigIntDivide(&ReTmp, norm, Re);
  subtractdivide(Re, -1, 2);
  if (!BigIntIsZero(Re))
  {
    Re->sign = signBak;
  }
  // Im <- ((Im1*Re2-Re1*Im2)*2/norm+1)/2
  (void)BigIntMultiply(ImDividend, ReDivisor, &ReTmp);
  (void)BigIntMultiply(ReDividend, ImDivisor, &ImTmp);
  BigIntSubt(&ReTmp, &ImTmp, &ImTmp);
  multint(&ImTmp, &ImTmp, 2);
  signBak = ImTmp.sign;
  ImTmp.sign = SIGN_POSITIVE;
  (void)BigIntDivide(&ImTmp, norm, Im);
  subtractdivide(Im, -1, 2);
  if (!BigIntIsZero(Im))
  {
    Im->sign = signBak;
  }
  // Re1 <- Re1 - Re*Re2 + Im*Im2
  (void)BigIntMultiply(Re, ReDivisor, &ReTmp);
  BigIntSubt(ReDividend, &ReTmp, ReDividend);
  (void)BigIntMultiply(Im, ImDivisor, &ReTmp);
  BigIntAdd(ReDividend, &ReTmp, ReDividend);
  // Im1 <- Im1 - Im*Re2 - Re*Im2
  (void)BigIntMultiply(Im, ReDivisor, &ReTmp);
  BigIntSubt(ImDividend, &ReTmp, ImDividend);
  (void)BigIntMultiply(Re, ImDivisor, &ReTmp);
  BigIntSubt(ImDividend, &ReTmp, ImDividend);
}

static int ComputeGCD(void)
{
  BigInteger Re;
  BigInteger Im;
  BigInteger Tmp;
  BigInteger Re1 = stackRealValues[stackIndex];
  BigInteger Re2 = stackRealValues[stackIndex+1];
  BigInteger Im1 = stackImagValues[stackIndex];
  BigInteger Im2 = stackImagValues[stackIndex+1];
  while (!BigIntIsZero(&Re2) || !BigIntIsZero(&Im2))
  {   // Second argument is not zero.
    BigInteger norm;
    (void)BigIntMultiply(&Re2, &Re2, &Re);
    (void)BigIntMultiply(&Im2, &Im2, &Im);
    BigIntAdd(&Re, &Im, &norm);
    // Get remainder of (Re1+i*Im1)/(Re2*i*Im2)
    // Overwrite Re1 and Im1 with that remainder.
    // Re and Im are used as temporary values.
    GetRemainder(&norm, &Re1, &Im1, &Re2, &Im2, &Re, &Im);
    // Exchange Re1 and Re2.
    CopyBigInt(&Tmp, &Re1);
    CopyBigInt(&Re1, &Re2);
    CopyBigInt(&Re2, &Tmp);
    // Exchange Im1 and Im2.
    CopyBigInt(&Tmp, &Im1);
    CopyBigInt(&Im1, &Im2);
    CopyBigInt(&Im2, &Tmp);
    Re2.sign = SIGN_POSITIVE;  // Re2 <- abs(Re2)
    Im2.sign = SIGN_POSITIVE;  // Im2 <- abs(Im2)
  }
  CopyBigInt(&Tmp, &Im1);
  Tmp.sign = SIGN_POSITIVE;  // Tmp <- abs(Tmp)
  BigIntSubt(&Re1, &Tmp, &Tmp);
  while (Tmp.sign == SIGN_NEGATIVE)
  {    // Multiply by i.
    CopyBigInt(&Re, &Re1);
    BigIntNegate(&Im1, &Re1);
    CopyBigInt(&Im1, &Re);
  }
  stackRealValues[stackIndex] = Re1;
  stackImagValues[stackIndex] = Im1;
  return 0;
}

static int ComputePower(const BigInteger *Re1, const BigInteger *Re2,
  const BigInteger *Im1, const BigInteger *Im2)
{
  unsigned int expon;
  double base;
  bool performPower = false;
  BigInteger norm;
  BigInteger ReTmp;
  BigInteger ImTmp;
  BigInteger Re;
  BigInteger Im;

  if (!BigIntIsZero(Im2) || (Re2->sign == SIGN_NEGATIVE))
  {          // Exponent must be positive or zero.
    return EXPR_INVALID_PARAM;
  }
  if (Re2->nbrLimbs > 2)
  {
    return EXPR_INTERM_TOO_HIGH;
  }
  expon = (unsigned int)Re2->limbs[0].x;
  if (Re2->nbrLimbs == 2)
  {
    expon += (unsigned int)Re2->limbs[1].x << BITS_PER_GROUP;
  }
  (void)BigIntMultiply(Re1, Re1, &ReTmp);
  (void)BigIntMultiply(Im1, Im1, &ImTmp);
  BigIntAdd(&ReTmp, &ImTmp, &norm);  // norm <- re1^2 + im1^2.
  if (norm.nbrLimbs > 1)
  {
    base = log((double)(norm.limbs[norm.nbrLimbs - 2].x + (norm.limbs[norm.nbrLimbs - 1].x << BITS_PER_GROUP)) +
      (double)(norm.nbrLimbs - 2) * LOG_2 * (double)BITS_PER_GROUP);
  }
  else
  {
    base = log((double)(norm.limbs[norm.nbrLimbs - 1].x) +
      (double)(norm.nbrLimbs - 1) * LOG_2 * (double)BITS_PER_GROUP);
  }
  if (base*(double)expon > 23026)
  {   // More than 20000 digits. 23026 = log(10^10000) (norm is already squared).
    return EXPR_INTERM_TOO_HIGH;
  }
  /* Compute actual power */
  Re.limbs[0].x = 1;      // Initialize Re <- 1, Im <- 0.
  Im.limbs[0].x = 0;
  Re.nbrLimbs = 1;
  Im.nbrLimbs = 1;
  Re.sign = SIGN_POSITIVE; 
  Im.sign = SIGN_POSITIVE;
  for (unsigned int mask = 1U << 30; mask != 0; mask >>= 1)
  {
    if ((expon & mask) != 0)
    {
      performPower = true;
    }
    if (performPower)
    {
      (void)BigIntMultiply(&Re, &Re, &ReTmp);         // ReTmp <- re*re - im*im.
      (void)BigIntMultiply(&Im, &Im, &ImTmp);
      BigIntSubt(&ReTmp, &ImTmp, &ReTmp);
      (void)BigIntMultiply(&Re, &Im, &Im);            // Im <- 2*re*im
      BigIntAdd(&Im, &Im, &Im);
      CopyBigInt(&Re, &ReTmp);
      if ((expon & mask) != 0)
      {
        (void)BigIntMultiply(Re1, &Re, &ReTmp);       // Re2 <- re1*re - im1*im.
        (void)BigIntMultiply(Im1, &Im, &ImTmp);
        BigIntSubt(&ReTmp, &ImTmp, &stackRealValues[stackIndex]);
        (void)BigIntMultiply(Re1, &Im, &ReTmp);       // Im <- re1*im + im1*re.
        (void)BigIntMultiply(Im1, &Re, &ImTmp);
        BigIntAdd(&ReTmp, &ImTmp, &Im);
        CopyBigInt(&Re, &stackRealValues[stackIndex]);
      }
    }
  }
  CopyBigInt(&stackRealValues[stackIndex], &Re);
  CopyBigInt(&stackImagValues[stackIndex], &Im);
  return EXPR_OK;
}

static int ComputeModPow(void)
{
  static BigInteger Result[2];
  static BigInteger norm;
  static BigInteger ReTmp;
  static BigInteger ImTmp;
  static BigInteger ReBase;
  static BigInteger ImBase;
  static BigInteger ReExp;
  static BigInteger ImExp;
  static BigInteger ReMod;
  static BigInteger ImMod;
  static BigInteger Re;
  static BigInteger Im;
  CopyBigInt(&ReBase, &stackRealValues[stackIndex]);
  CopyBigInt(&ImBase, &stackImagValues[stackIndex]);
  CopyBigInt(&ReExp, &stackRealValues[stackIndex+1]);
  CopyBigInt(&ImExp, &stackImagValues[stackIndex+1]);
  CopyBigInt(&ReMod, &stackRealValues[stackIndex+2]);
  CopyBigInt(&ImMod, &stackImagValues[stackIndex+2]);

  if (!BigIntIsZero(&ImExp))
  {         // Imaginary part must be zero.
    return EXPR_INVALID_PARAM;
  }
  if (ReExp.sign == SIGN_NEGATIVE)
  {
    int retcode;
    BigIntNegate(&ReExp, &ReExp);
    retcode = ModInv(&ReBase, &ImBase, &ReMod, &ImMod, Result);
    if (retcode != 0)
    {
      return retcode;
    }
    ReBase = Result[0];
    ImBase = Result[1];
  }
  // Re <- 1, Im <- 0
  Re.limbs[0].x = 1;
  Im.limbs[0].x = 0;
  Re.nbrLimbs = 1;
  Im.nbrLimbs = 1;
  Re.sign = SIGN_POSITIVE; 
  Im.sign = SIGN_POSITIVE;
  if (BigIntIsZero(&ReMod) && BigIntIsZero(&ImMod))
  {   /* Modulus is zero */
    (void)ComputePower(&ReBase, &ReExp, &ImBase, &ImExp);
  }
  else
  {                            /* Modulus is not zero */
    (void)BigIntMultiply(&ReMod, &ReMod, &ReTmp);
    (void)BigIntMultiply(&ImMod, &ImMod, &ImTmp);
    BigIntAdd(&ReTmp, &ImTmp, &norm);
    for (int index = ReExp.nbrLimbs - 1; index >= 0; index--)
    {
      int groupExp = (int)ReExp.limbs[index].x;
      for (int mask = 1 << (BITS_PER_GROUP - 1); mask > 0; mask >>= 1)
      {
        // Let Re + i*Im <- (Re + i*Im)^2
        (void)BigIntMultiply(&Re, &Re, &ReTmp);
        (void)BigIntMultiply(&Im, &Im, &ImTmp);
        BigIntSubt(&ReTmp, &ImTmp, &ReTmp);
        (void)BigIntMultiply(&Re, &Im, &ImTmp);
        BigIntAdd(&ImTmp, &ImTmp, &Im);
        CopyBigInt(&Re, &ReTmp);
        // Replace (Re + i*Im) by the remainder of
        // (ReTmp + i*ImTmp) / (ReMod + i*ImMod)
        // ReTmp and ImTmp are used as temporary values.
        GetRemainder(&norm, &Re, &Im, &ReMod, &ImMod, &ReTmp, &ImTmp);
        if ((groupExp & mask) != 0)
        {
          // Let Re + i*Im <- (Re + i*Im)*(ReBase + i*ImBase)
          (void)BigIntMultiply(&ReBase, &Re, &ReTmp);
          (void)BigIntMultiply(&ImBase, &Im, &ImTmp);
          BigIntSubt(&ReTmp, &ImTmp, &ReTmp);
          (void)BigIntMultiply(&ImBase, &Re, &ImTmp);
          (void)BigIntMultiply(&ReBase, &Im, &Re);
          BigIntAdd(&Re, &ImTmp, &Im);
          CopyBigInt(&Re, &ReTmp);
          // Replace (Re + i*Im) by the remainder of
          // (ReTmp + i*ImTmp) / (ReMod + i*ImMod)
          // ReTmp and ImTmp are used as temporary values.
          GetRemainder(&norm, &Re, &Im, &ReMod, &ImMod, &ReTmp, &ImTmp);
        }
      }
    }
    (void)Modulo(&Re, &Im, &ReMod, &ImMod, Result);
    stackRealValues[stackIndex] = Result[0];
    stackImagValues[stackIndex] = Result[1];
    return 0;
  }
  return 0;
}

static int ComputeModInv(void)
{
  BigInteger Result[2];
  int retcode;

  retcode = ModInv(&stackRealValues[stackIndex],
    &stackImagValues[stackIndex],
    &stackRealValues[stackIndex + 1],
    &stackImagValues[stackIndex + 1],
    Result);
  if (retcode != 0)
  {
    return retcode;
  }
  stackRealValues[stackIndex] = Result[0];
  stackImagValues[stackIndex] = Result[1];
  return 0;
}

static int ModInv(const BigInteger *RealNbr, const BigInteger *ImagNbr,
                  BigInteger *RealMod, BigInteger *ImagMod,
                  BigInteger *Result)
{
  static BigInteger ReG0;
  static BigInteger ReG1;
  static BigInteger ImG0;
  static BigInteger ImG1;
  static BigInteger ReU0;
  static BigInteger ReU1;
  static BigInteger ImU0;
  static BigInteger ImU1;
  static BigInteger Re;
  static BigInteger Im;
  static BigInteger Tmp;

  if (BigIntIsZero(RealMod) && BigIntIsZero(ImagMod))
  {          // Argument is zero.
    return EXPR_INVALID_PARAM;
  }
  CopyBigInt(&ReG0, RealNbr);
  CopyBigInt(&ImG0, ImagNbr);
  CopyBigInt(&ReG1, RealMod);
  CopyBigInt(&ImG1, ImagMod);
   // Initialize U0 <- 1, U1 <- 0.
  ReU0.limbs[0].x = 1;
  ReU0.nbrLimbs = 1;
  ReU0.sign = SIGN_POSITIVE;
  ReU1.limbs[0].x = 0; 
  ImU0.limbs[0].x = 0; 
  ImU1.limbs[0].x = 0;
  ReU1.nbrLimbs = 1; 
  ImU0.nbrLimbs = 1; 
  ImU1.nbrLimbs = 1;
  ReU1.sign = SIGN_POSITIVE; 
  ImU0.sign = SIGN_POSITIVE; 
  ImU1.sign = SIGN_POSITIVE;
  while (!BigIntIsZero(&ReG1) || !BigIntIsZero(&ImG1))
  {            // G1 is not zero.
    BigInteger norm;
    (void)BigIntMultiply(&ReG1, &ReG1, &Re);
    (void)BigIntMultiply(&ImG1, &ImG1, &Im);
    BigIntAdd(&Re, &Im, &norm);
    // Replace G0 by the remainder of G0/G1.
    GetRemainder(&norm, &ReG0, &ImG0, &ReG1, &ImG1, &Re, &Im);
    // Exchange G0 and G1.
    CopyBigInt(&Tmp, &ReG0);
    CopyBigInt(&ReG0, &ReG1);
    CopyBigInt(&ReG1, &Tmp);
    CopyBigInt(&Tmp, &ImG0);
    CopyBigInt(&ImG0, &ImG1);
    CopyBigInt(&ImG1, &Tmp);
    (void)BigIntMultiply(&Re, &ReU1, &Tmp);
    BigIntSubt(&ReU0, &Tmp, &ReU0);
    (void)BigIntMultiply(&Im, &ImU1, &Tmp);
    BigIntAdd(&ReU0, &Tmp, &ReU0);
    (void)BigIntMultiply(&Im, &ReU1, &Tmp);
    BigIntSubt(&ImU0, &Tmp, &ImU0);
    (void)BigIntMultiply(&Re, &ImU1, &Tmp);
    BigIntSubt(&ImU0, &Tmp, &ImU0);
    // Exchange U0 and U1.
    CopyBigInt(&Tmp, &ReU0);
    CopyBigInt(&ReU0, &ReU1);
    CopyBigInt(&ReU1, &Tmp);
    CopyBigInt(&Tmp, &ImU0);
    CopyBigInt(&ImU0, &ImU1);
    CopyBigInt(&ImU1, &Tmp);
  }

  for (;;)
  {
    CopyBigInt(&Tmp, &ImG0);
    Tmp.sign = SIGN_POSITIVE;   // Tmp <- abs(Tmp)
    BigIntSubt(&ReG0, &Tmp, &Tmp);
    if (Tmp.sign == SIGN_POSITIVE)
    {
      break;
    }
      // Multiply G0 and U0 by i.
    CopyBigInt(&Tmp, &ReG0);
    BigIntNegate(&ImG0, &ReG0);
    CopyBigInt(&ImG0, &Tmp);
    CopyBigInt(&Tmp, &ReU0);
    BigIntNegate(&ImU0, &ReU0);
    CopyBigInt(&ImU0, &Tmp);
  }
  if ((ReG0.nbrLimbs != 1) || (ReG0.limbs[0].x != 1) || !BigIntIsZero(&ImG0))
  {         // G0 is not 1.
    return EXPR_INVALID_PARAM;
  }
  return Modulo(&ReU0, &ImU0, RealMod, ImagMod, Result);
}

static int Modulo(BigInteger *ReNum, BigInteger *ImNum,
                  const BigInteger *ReDen, const BigInteger *ImDen,
                  BigInteger *Result)
{
  BigInteger Re;
  BigInteger Im;
  BigInteger ReMin;
  BigInteger ImMin;
  BigInteger norm;
  BigInteger Tmp;
  BigInteger normmin;
  if (BigIntIsZero(ReDen) && BigIntIsZero(ImDen))
  {      // Denominator is zero.
    CopyBigInt(Result, ReNum);
    CopyBigInt(Result+1, ImNum);
    return 0;
  }
  ReMin.limbs[0].x = 0; 
  ImMin.limbs[0].x = 0;
  ReMin.nbrLimbs = 1; 
  ImMin.nbrLimbs = 1;
  ReMin.sign = SIGN_POSITIVE; 
  ImMin.sign = SIGN_POSITIVE;
  (void)BigIntMultiply(ReDen, ReDen, &norm);
  (void)BigIntMultiply(ImDen, ImDen, &Tmp);
  BigIntAdd(&norm, &Tmp, &norm);
  // Replace Num by the remainder of Num/Den.
  GetRemainder(&norm, ReNum, ImNum, ReDen, ImDen, &Re, &Im);
  // Initialize normmin to -1.
  normmin.limbs[0].x = 1;
  normmin.nbrLimbs = 1;
  normmin.sign = SIGN_NEGATIVE;
  for (int i=0; i<9; i++)
  {
    switch (i)
    {
      case 0:
        CopyBigInt(&Re, ReNum);
        CopyBigInt(&Im, ImNum);
        break;
      case 1:
        BigIntAdd(ReNum, ReDen, &Re);
        BigIntAdd(ImNum, ImDen, &Im);
        break;
      case 2:
        BigIntAdd(ReNum, ReDen, &Re);
        BigIntSubt(&Re, ImDen, &Re);
        BigIntAdd(ImNum, ReDen, &Im);
        BigIntAdd(&Im, ImDen, &Im);
        break;
      case 3:
        BigIntSubt(ReNum, ImDen, &Re);
        BigIntAdd(ImNum, ReDen, &Im);
        break;
      case 4:
        BigIntSubt(ReNum, ReDen, &Re);
        BigIntSubt(&Re, ImDen, &Re);
        BigIntAdd(ImNum, ReDen, &Im);
        BigIntSubt(&Im, ImDen, &Im);
        break;
      case 5:
        BigIntSubt(ReNum, ReDen, &Re);
        BigIntSubt(ImNum, ImDen, &Im);
        break;
      case 6:
        BigIntSubt(ReNum, ReDen, &Re);
        BigIntAdd(&Re, ImDen, &Re);
        BigIntSubt(ImNum, ReDen, &Im);
        BigIntSubt(&Im, ImDen, &Im);
        break;
      case 7:
        BigIntAdd(ReNum, ImDen, &Re);
        BigIntSubt(ImNum, ReDen, &Im);
        break;
      case 8:
        BigIntAdd(ReNum, ReDen, &Re);
        BigIntAdd(&Re, ImDen, &Re);
        BigIntSubt(ImNum, ReDen, &Im);
        BigIntAdd(&Im, ImDen, &Im);
        break;
      default:
        break;
    }
    if (Re.sign == SIGN_POSITIVE)
    {
      (void)BigIntMultiply(&Re, &Re, &norm);
      (void)BigIntMultiply(&Im, &Im, &Tmp);
      BigIntAdd(&norm, &Tmp, &norm);
      BigIntSubt(&norm, &normmin, &Tmp);
      if ((normmin.sign == SIGN_NEGATIVE) || (Tmp.sign == SIGN_NEGATIVE))
      {
        CopyBigInt(&normmin, &norm);
        CopyBigInt(&ReMin, &Re);
        CopyBigInt(&ImMin, &Im);
      }
    }
  }
  CopyBigInt(&Result[0], &ReMin);
  CopyBigInt(&Result[1], &ImMin);
  return 0;
}

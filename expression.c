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
#if BITS_PER_GROUP == 15
  #define SQRT_MAX_VALUE_LIMB 182
#else
  #define SQRT_MAX_VALUE_LIMB 46341
#endif
static BigInteger stackValues[PAREN_STACK_SIZE];
static char stackOperators[PAREN_STACK_SIZE];
static int partArray[100000+1000];
static limb prodModulus[MAX_LEN];
extern limb MontgomeryR1[MAX_LEN];
extern limb TestNbr[MAX_LEN];
extern limb MontgomeryMultN[MAX_LEN];
extern limb MontgomeryMultR1[MAX_LEN];
int lang;
char *output;
limb Mult1[MAX_LEN];
limb Mult3[MAX_LEN];
limb Mult4[MAX_LEN];
int q[MAX_LEN];
#define fibon1 MontgomeryR1
#define fibon2 prodModulus
static int prodModulusLimbs;
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
  int exprIndexAux;
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
        carry.x = 0;
        for (j = 0; j < factorial.nbrLimbs; j++)
        {
          carry.x += i*factorial.limbs[j].x;
          factorial.limbs[j].x = carry.x & MAX_VALUE_LIMB;
          carry.x >>= BITS_PER_GROUP;
        }
        if (carry.x != 0)
        {  // New limb needed.
          factorial.limbs[j].x = carry.x;
          factorial.nbrLimbs++;
        }
      }
      stackValues[stackIndex] = factorial;
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
          carry.x = 0;
          for (j = 0; j < factorial.nbrLimbs; j++)
          {
            carry.x += i*factorial.limbs[j].x;
            factorial.limbs[j].x = carry.x & MAX_VALUE_LIMB;
            carry.x >>= BITS_PER_GROUP;
          }
          if (carry.x != 0)
          {  // New limb needed.
            factorial.limbs[j].x = carry.x;
            factorial.nbrLimbs++;
          }
        }
      }
      stackValues[stackIndex] = factorial;
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
          }
          (ptrLimb++)->x = carry.x & MAX_VALUE_LIMB;
          carry.x >>= BITS_PER_GROUP;
        }
        stackValues[stackIndex].nbrLimbs = (int)(ptrLimb - &stackValues[stackIndex].limbs[0]);
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
      }
      leftNumberFlag = TRUE;
      exprIndex = exprIndexAux;
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
  *ExpressionResult = stackValues[0];
  return EXPR_OK;
}

  // Compute the partitions of an integer p(n)
  // The formula p(n) = p(k) = p(k - 1) + p(k - 2) - p(k - 5) -
  // - p(k - 7) + p(k - 12) + p(k - 15) - p(k - 22) - ...
  // where p(0) = 1, where the numbers have the form n(3n-1)/2 for
  // n = 1, -1, 2, -2, 3, -3, ... stopping when k-n(3n-1)/2 is negative.
  // The signs in the summation are: +, +, -, -, +, +, -, -, ...
  // Modular arithmetic is used modulo the largest primes under 2^31, and
  // the result is reconstructed.
  // The approximation p(n) = e^(pi*sqrt(2n/3))/(4n*sqrt(3)) is used
  // in order to find the number of limbs needed to compute the partition.
static void partition(int val, BigInteger *pResult)
{
  int index, currentPrime, Q, k, n, sum, idx;
  limb carry;
    // Compute approximate number of limbs: log(p(n))/log(2^31)
    // pi * sqrt(2/3)/log(2^31) < 0.12, so 0.12 is selected.
    // for 15 bits is < 0.25
#if BITS_PER_GROUP == 15
  int limbs = (int)(0.25*sqrt((double)val) + 1);
#else
  int limbs = (int)(0.12*sqrt(val) + 1);
#endif

    // Compute the primes which will be used for the modular arithmetic
    // operations. The primes must be ordered in ascending order.
  currentPrime = MAX_VALUE_LIMB+2;    // Greatest number representable + 2.
  for (index = limbs + val - 1; index >= val; index--)
  {
    for (;;)
    {         // Loop that computes the previous prime number.
      do
      {
        currentPrime -= 2;
      } while (currentPrime % 3 == 0);
      for (Q = 5; Q <= SQRT_MAX_VALUE_LIMB; Q += 6)
      { /* Check if Base is prime */
        if (currentPrime % Q == 0)
        {
          break;     /* Composite */
        }
        if (currentPrime % (Q + 2) == 0)
        {
          break;    /* Composite */
        }
      }
      if (Q > SQRT_MAX_VALUE_LIMB)
      {
        break; /* Prime found */
      }
    }
    partArray[index] = currentPrime;
  }
    // Perform modular arithmetic.
  for (index = val; index<val + limbs; index++)
  {
    currentPrime = partArray[index];
    sum = 1;                          // Initialize p(0) mod currentPrime.
    for (k = 1; k <= val; k++)        // Generate all partition numbers
    {                                 // up to the one wanted.
      idx = k;
      partArray[k - 1] = sum;         // Store p(k-1) mod currentPrime.
      sum = 0;
      n = 1;
      for (;;)                        // Loop for n.
      {
        idx -= n + n - 1;
        if (idx < 0)
        {
          break;                      // Negative index, so go out.
        }
        sum -= currentPrime - partArray[idx];
        sum += currentPrime & (sum >> 31);
        idx -= n;
        if (idx < 0)
        {
          break;                      // Negative index, so go out.
        }
        sum -= currentPrime - partArray[idx];
        sum += currentPrime & (sum >> 31);
        n++;
        idx -= n + n - 1;
        if (idx < 0)
        {
          break;                      // Negative index, so go out.
        }
        sum -= partArray[idx];
        sum += currentPrime & (sum >> 31);
        idx -= n;
        if (idx < 0)
        {
          break;                      // Negative index, so go out.
        }
        sum -= partArray[idx];
        sum += currentPrime & (sum >> 31);
        n++;
      }
    }
    partArray[index + limbs] = sum;
  }
    // Reconstruct the result from p(val) mod all primes.
    // v_1 <- u_1 mod m_1
    // v_2 <- (u_2 - v_1) (m_1)^(-1) mod m_2
    // v_3 <- (u_3 - (v_1+m_1 v_2)) (m_1 m_2)^(-1) mod m_3
    // ...
    // v_r <- (u_n - (v_1+m_1(v_2+m_2(v_3+...+m_(r-2) v_(r-1)...)))*
    //        (m_1 m_2 ... m_(r-1))^(-1) mod m_r
    //
    // u = v_r m_(r-1) ..m_2 m_1 + v_3 m_2 m_1 + v_2 m_1 + v_1
  partArray[0] = partArray[val + limbs];
  for (index = 1; index<limbs; index++)
  {
    long numerator, prodmod;
    currentPrime = partArray[val + index];
    prodmod = 1;
    for (k = index - 1; k >= 0; k--)
    {
      prodmod = prodmod * partArray[val + k] % currentPrime;
    }
    prodmod = modInv((int)prodmod, currentPrime);
    numerator = partArray[index - 1];
    for (k = index - 2; k >= 0; k--)
    {
      numerator = (numerator*partArray[val + k] + partArray[k]) %
          currentPrime;
    }
    sum = partArray[val + limbs + index] - (int)numerator;
    if (sum<0)
    {
      sum += currentPrime;
    }
    partArray[index] = (int)(sum*prodmod%currentPrime);
  }
  // Use Chinese Remainder Theorem to find the partition number from
  // the partition number mod different primes.
  pResult->limbs[0].x = partArray[0];
  pResult->nbrLimbs = 1;
  pResult->sign = SIGN_POSITIVE;
  prodModulus[0].x = 1;
  prodModulusLimbs = 1;
  for (index = 1; index<limbs; index++)
  {
    int mult;
    // Update product of modulus by multiplying by next prime.
    carry.x = 0;
    mult = partArray[val + index - 1];
    for (idx = 0; idx < prodModulusLimbs; idx++)
    {
      carry.x += mult*prodModulus[idx].x;
      prodModulus[idx].x = carry.x & MAX_VALUE_LIMB;
      carry.x >>= BITS_PER_GROUP;
    }
    if (carry.x != 0)
    {  // New limb needed.
      prodModulus[idx].x = carry.x;
      pResult->limbs[idx].x = 0;
      prodModulusLimbs++;
    }
    // Update result.
    carry.x = 0;
    mult = partArray[index];
    for (idx = 0; idx < prodModulusLimbs; idx++)
    {
      carry.x += mult*prodModulus[idx].x + pResult->limbs[idx].x;
      pResult->limbs[idx].x = carry.x & MAX_VALUE_LIMB;
      carry.x >>= BITS_PER_GROUP;
    }
    if (carry.x != 0)
    {  // New limb needed.
      prodModulus[idx].x = 0;
      pResult->limbs[idx].x = carry.x;
      prodModulusLimbs++;
    }
  }
  while (pResult->limbs[prodModulusLimbs - 1].x == 0)
  {
    prodModulusLimbs--;
  }
  pResult->nbrLimbs = prodModulusLimbs;
  return;
}

static enum eExprErr ComputeSubExpr(int stackIndex)
{
  int len, val, ctr;
  BigInteger *pArgument, *pResult;
  limb *pResultLimbs, *pArgumentLimbs, *pTemp, *pFibonPrev, *pFibonAct;
  limb carry, largeVal;
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
      fibon1[0].x = (stackOper == 'L'? 2: 0);
      fibon2[0].x = 1;
      for (i = 1; i <= val; i++)
      {
        carry.x = 0;
        for (j=0; j<len; j++)
        {
          carry.x += (pFibonPrev+j)->x + (pFibonAct+j)->x;
          (pFibonPrev+j)->x = carry.x & MAX_VALUE_LIMB;
          carry.x >>= BITS_PER_GROUP;
        }
        if (carry.x != 0)
        {
          (pFibonPrev+j)->x = carry.x;
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
    {
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
    {
      pResult -> sign = SIGN_POSITIVE;
      if (pArgument->sign == SIGN_NEGATIVE || (nbrLimbs == 1 && pArgumentLimbs->x<2))
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
      {
        if (FirstTime != 0)
        {
          pResultLimbs->x |= 1;  // If it is even, change to next odd value.
          FirstTime = 0;
          pResult->nbrLimbs = nbrLimbs;
        }
        else
        {     // Add 2.
          pResultLimbs->x += 2;
        }
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
      if (isPseudoprime(pResult))
      {   // Number is SPRP.
        pResult->nbrLimbs = nbrLimbs;
        pResult->sign = SIGN_POSITIVE;
        return EXPR_OK;
      }
    }
  }              /* end switch */
  return EXPR_OK;
}

void textError(char *output, enum eExprErr rc)
{
  switch (rc)
  {
  case EXPR_NUMBER_TOO_LOW:
    strcpy(output, "Number too low");
    break;
  case EXPR_NUMBER_TOO_HIGH:
    strcpy(output, "Number too high (more than 10000 digits)");
    break;
  case EXPR_INTERM_TOO_HIGH:
    strcpy(output, "Intermediate number too high (more than 20000 digits)");
    break;
  case EXPR_DIVIDE_BY_ZERO:
    strcpy(output, "Division by zero");
    break;
  case EXPR_PAREN_MISMATCH:
    strcpy(output, "Parenthesis mismatch");
    break;
  case EXPR_SYNTAX_ERROR:
    strcpy(output, "Syntax error");
    break;
  case EXPR_TOO_MANY_PAREN:
    strcpy(output, "Too many parenthesis");
    break;
  case EXPR_INVALID_PARAM:
    strcpy(output, "Invalid parameter");
    break;
  case EXPR_BREAK:
    strcpy(output, "Stopped by user");
    break;
  default: 
    break;
  }
}

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

#ifdef FACTORIZATION_FUNCTIONS
#define TOKEN_TOTIENT    32
#define TOKEN_SUMDIVS    33
#define TOKEN_NUMDIVS    34
#define TOKEN_MINFACT    35
#define TOKEN_MAXFACT    36
#define TOKEN_NUMFACT    37
#define TOKEN_CONCATFACT 38
#endif
#define TOKEN_PRIMORIAL  39
#define TOKEN_GCD        40
#define TOKEN_LCM        41
#define TOKEN_MODPOW     42
#define TOKEN_MODINV     43
#define TOKEN_SUMDIGITS  44
#define TOKEN_NUMDIGITS  45
#define TOKEN_REVDIGITS  46
#define TOKEN_ISPRIME    47
#define TOKEN_JACOBI     48
#define TOKEN_SQRT       49
#define TOKEN_F          50
#define TOKEN_L          51
#define TOKEN_P          52
#define TOKEN_N          53
#define TOKEN_B          54

#define PAREN_STACK_SIZE           5000
#define COMPR_STACK_SIZE        1000000

#define DO_NOT_SHORT_CIRCUIT  (COMPR_STACK_SIZE + 100)  // Larger than stack size.

struct sFuncOperExpr stFuncOperIntExpr[] =
{
  // First section: functions
#ifdef FACTORIZATION_FUNCTIONS
  {"TOTIENT", TOKEN_TOTIENT + ONE_PARM, 0},
  {"SUMDIVS", TOKEN_SUMDIVS + ONE_PARM, 0},
  {"NUMDIVS", TOKEN_NUMDIVS + ONE_PARM, 0},
  {"MINFACT", TOKEN_MINFACT + ONE_PARM, 0},
  {"MAXFACT", TOKEN_MAXFACT + ONE_PARM, 0},
  {"NUMFACT", TOKEN_NUMFACT + ONE_PARM, 0},
  {"CONCATFACT", TOKEN_CONCATFACT + TWO_PARMS, 0},
#endif
  {"GCD", TOKEN_GCD + MANY_PARMS, 0},
  {"LCM", TOKEN_LCM + MANY_PARMS, 0},
  {"MODPOW", TOKEN_MODPOW + THREE_PARMS, 0},
  {"MODINV", TOKEN_MODINV + TWO_PARMS, 0},
  {"SUMDIGITS", TOKEN_SUMDIGITS + TWO_PARMS, 0},
  {"NUMDIGITS", TOKEN_NUMDIGITS + TWO_PARMS, 0},
  {"REVDIGITS", TOKEN_REVDIGITS + TWO_PARMS, 0},
  {"ISPRIME", TOKEN_ISPRIME + ONE_PARM, 0},
  {"JACOBI", TOKEN_JACOBI + TWO_PARMS, 0},
  {"SQRT", TOKEN_SQRT + ONE_PARM, 0},
  {"F", TOKEN_F + ONE_PARM, 0},
  {"L", TOKEN_L + ONE_PARM, 0},
  {"P", TOKEN_P + ONE_PARM, 0},
  {"N", TOKEN_N + ONE_PARM, 0},
  {"B", TOKEN_B + ONE_PARM, 0},
  {NULL, 0, 0},
  // Second section: functions written at right of argument.
  {"#", TOKEN_PRIMORIAL, 0},
  {NULL, 0, 0},
  // Third section: unary operators.
  {"-", OPER_UNARY_MINUS, 3},
  {"NOT", OPER_NOT, 7},
  {NULL, 0, 0},
  // Fourth section: binary operators.
  {"**", OPER_POWER, 1}, // This must be located before multiplication operator.
  {"XOR", OPER_XOR, 8},
  {"OR", OPER_OR, -8},
  {"AND", OPER_AND, -8},
  {"==", OPER_EQUAL, 6},
  {"!=", OPER_NOT_EQUAL, 6},
  {"<<", OPER_SHL, 5},
  {">>", OPER_SHR, 5},
  {">=", OPER_NOT_LESS, 6},
  {"<=", OPER_NOT_GREATER, 6},
  {"<", OPER_LESS, 6},
  {">", OPER_GREATER, 6},
  {"SHL", OPER_SHL, 5},
  {"SHR", OPER_SHR, 5},
  {"+", OPER_PLUS, 4},
  {"-", OPER_MINUS, 4},
  {"*", OPER_MULTIPLY, 2},
  {"%", OPER_REMAINDER, 2},
  {"/", OPER_DIVIDE, 2},
  {"^", OPER_POWER, 1},
  {NULL, 0, 0},
};

static int comprStackValues[COMPR_STACK_SIZE];
static int comprStackOffset[PAREN_STACK_SIZE];
extern limb MontgomeryR1[MAX_LEN];
static int stackIndex;
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
#ifdef FACTORIZATION_FUNCTIONS
static enum eExprErr ComputeTotient(void);
static enum eExprErr ComputeNumFact(void);
static enum eExprErr ComputeMaxFact(void);
static enum eExprErr ComputeMinFact(void);
static enum eExprErr ComputeNumDivs(void);
static enum eExprErr ComputeSumDivs(void);
static enum eExprErr ComputeConcatFact(void);
static char textFactor[MAX_LEN*12];
#endif
static int ComputeSumDigits(void);
static int ComputeRevDigits(void);
static int ComputeNumDigits(void);
static enum eExprErr ComputeModInv(void);
static enum eExprErr ComputePartition(void);
static enum eExprErr ShiftLeft(BigInteger* first, const BigInteger* second, BigInteger* result);
static BigInteger curStack;
static BigInteger curStack2;
static BigInteger curStack3;

static int numLimbs(const int* pLen)
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
  const int* ptrStackValue = &comprStackValues[comprStackOffset[stackIndex]];
  NumberLength = numLimbs(ptrStackValue);
  IntArray2BigInteger(ptrStackValue, pValue);
}

static enum eExprErr setStackValue(const BigInteger* pValue)
{
  int currentOffset = comprStackOffset[stackIndex];
  if (currentOffset >= (COMPR_STACK_SIZE - (int)sizeof(BigInteger) / (int)sizeof(limb)))
  {
    return EXPR_OUT_OF_MEMORY;
  }
  NumberLength = pValue->nbrLimbs;
  BigInteger2IntArray(&comprStackValues[currentOffset], pValue);
  comprStackOffset[stackIndex + 1] = currentOffset + pValue->nbrLimbs + 1;
  return EXPR_OK;
}

enum eExprErr ComputeExpression(const char *expr, BigInteger *ExpressionResult)
{
  bool valueXused;
  enum eExprErr retcode;
  char* ptrRPNbuffer;
  int len;
  int stackIndexThreshold = DO_NOT_SHORT_CIRCUIT;
  retcode = ConvertToReversePolishNotation(expr, &ptrRPNbuffer, stFuncOperIntExpr,
    PARSE_EXPR_INTEGER, &valueXused);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  stackIndex = -1;
  comprStackOffset[0] = 0;
  while (*ptrRPNbuffer != '\0')
  {
    char c = *ptrRPNbuffer;
    int currentOffset;
    int nbrLenBytes;
    int jacobi;
    int nbrParameters;
    switch (c)
    {
    case TOKEN_NUMBER:
      ptrRPNbuffer++;           // Skip token.
      stackIndex++;
      len = ((int)(unsigned char)*ptrRPNbuffer * 256) + (int)(unsigned char)*(ptrRPNbuffer + 1);
      nbrLenBytes = len * (int)sizeof(limb);
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        ptrRPNbuffer += nbrLenBytes;
        ptrRPNbuffer += 2;
        break;
      }
      // Move number to compressed stack.
      currentOffset = comprStackOffset[stackIndex];
      if (currentOffset >= (COMPR_STACK_SIZE - 
        ((int)sizeof(BigInteger) / (int)sizeof(limb))))
      {
        return EXPR_OUT_OF_MEMORY;
      }
      comprStackValues[currentOffset] = len;
      ptrRPNbuffer += 2;   // Skip length.
      (void)memcpy(&comprStackValues[currentOffset+1], ptrRPNbuffer, nbrLenBytes);
      ptrRPNbuffer += nbrLenBytes;
      comprStackOffset[stackIndex+1] = currentOffset + 1 + len;
      break;

    case TOKEN_VAR:
      stackIndex++;
      if (stackIndexThreshold >= stackIndex)
      {     // Part of second operand of binary AND/OR not short-circuited.
        CopyBigInt(&curStack, &valueX);
      }
      break;

    case TOKEN_COUNTER:
      stackIndex++;
      if (stackIndexThreshold >= stackIndex)
      {     // Part of second operand of binary AND/OR not short-circuited.
        intToBigInteger(&curStack, counterC);
      }
      break;

    case TOKEN_GCD:
      ptrRPNbuffer++;
      nbrParameters = (int)(unsigned char)*ptrRPNbuffer;
      if (stackIndexThreshold < stackIndex)
      {     // Delete from operand stack when binary AND/OR short-circuited.
        stackIndex -= nbrParameters - 1;
        break;
      }
      getCurrentStackValue(&curStack);
      for (int parmNbr = 1; parmNbr < nbrParameters; parmNbr++)
      {
        stackIndex--;
        getCurrentStackValue(&curStack2);
        BigIntGcd(&curStack, &curStack2, &curStack);
      }
      break;

    case TOKEN_LCM:
      ptrRPNbuffer++;
      nbrParameters = (int)(unsigned char)*ptrRPNbuffer;
      if (stackIndexThreshold < stackIndex)
      {     // Delete from operand stack when binary AND/OR short-circuited.
        stackIndex -= nbrParameters - 1;
        break;
      }
      getCurrentStackValue(&curStack);
      for (int parmNbr = 1; parmNbr < nbrParameters; parmNbr++)
      {
        stackIndex--;
        getCurrentStackValue(&curStack2);
        BigIntLcm(&curStack, &curStack2, &curStack);
      }
      break;

    case TOKEN_JACOBI:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        stackIndex--;
        break;
      }
      getCurrentStackValue(&curStack2);
      stackIndex--;
      getCurrentStackValue(&curStack);
      jacobi = BigIntJacobiSymbol(&curStack, &curStack2);
      intToBigInteger(&curStack, jacobi);
      break;

    case TOKEN_SQRT:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        break;
      }
      getCurrentStackValue(&curStack);
      if (curStack.sign == SIGN_NEGATIVE)
      {
        return EXPR_INVALID_PARAM;
      }
      squareRoot(curStack.limbs, curStack2.limbs, curStack.nbrLimbs, &curStack2.nbrLimbs);
      curStack2.sign = SIGN_POSITIVE;
      CopyBigInt(&curStack, &curStack2);
      break;

    case TOKEN_MODPOW:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        stackIndex -= 2;
        break;
      }
      getCurrentStackValue(&curStack3);
      stackIndex--;
      getCurrentStackValue(&curStack2);
      stackIndex--;
      getCurrentStackValue(&curStack);
      retcode = BigIntGeneralModularPower(&curStack, &curStack2, &curStack3, &curStack);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_MODINV:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        stackIndex--;
        break;
      }
      getCurrentStackValue(&curStack2);
      stackIndex--;
      getCurrentStackValue(&curStack);
      retcode = ComputeModInv();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

#ifdef FACTORIZATION_FUNCTIONS
    case TOKEN_TOTIENT:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        break;
      }
      getCurrentStackValue(&curStack);
      retcode = ComputeTotient();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_NUMDIVS:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        break;
      }
      getCurrentStackValue(&curStack);
      retcode = ComputeNumDivs();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_SUMDIVS:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        break;
      }
      getCurrentStackValue(&curStack);
      retcode = ComputeSumDivs();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_MINFACT:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        break;
      }
      getCurrentStackValue(&curStack);
      retcode = ComputeMinFact();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_MAXFACT:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        break;
      }
      getCurrentStackValue(&curStack);
      retcode = ComputeMaxFact();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_NUMFACT:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        break;
      }
      getCurrentStackValue(&curStack);
      retcode = ComputeNumFact();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_CONCATFACT:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        stackIndex--;
        break;
      }
      getCurrentStackValue(&curStack2);
      stackIndex--;
      getCurrentStackValue(&curStack);
      retcode = ComputeConcatFact();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;
#endif

    case TOKEN_SUMDIGITS:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        stackIndex--;
        break;
      }
      getCurrentStackValue(&curStack2);
      stackIndex--;
      getCurrentStackValue(&curStack);
      retcode = ComputeSumDigits();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_NUMDIGITS:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        stackIndex--;
        break;
      }
      getCurrentStackValue(&curStack2);
      stackIndex--;
      getCurrentStackValue(&curStack);
      retcode = ComputeNumDigits();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_REVDIGITS:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        stackIndex--;
        break;
      }
      getCurrentStackValue(&curStack2);
      stackIndex--;
      getCurrentStackValue(&curStack);
      retcode = ComputeRevDigits();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_ISPRIME:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        break;
      }
      getCurrentStackValue(&curStack);
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
      break;

    case TOKEN_F:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        break;
      }
      getCurrentStackValue(&curStack);
      retcode = ComputeFibLucas(0, &curStack);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_L:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        break;
      }
      getCurrentStackValue(&curStack);
      retcode = ComputeFibLucas(2, &curStack);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_P:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        break;
      }
      getCurrentStackValue(&curStack);
      retcode = ComputePartition();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_N:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        break;
      }
      getCurrentStackValue(&curStack);
      retcode = ComputeNext(&curStack);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_B:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        break;
      }
      getCurrentStackValue(&curStack);
      retcode = ComputeBack(&curStack);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_FACTORIAL:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        stackIndex--;
        break;
      }
      getCurrentStackValue(&curStack2);
      stackIndex--;
      getCurrentStackValue(&curStack);
      if (curStack.sign == SIGN_NEGATIVE)
      {
        return EXPR_INVALID_PARAM;
      }
      retcode = factorial(&curStack, curStack.limbs[0].x, curStack2.limbs[0].x);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_PRIMORIAL:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        break;
      }
      getCurrentStackValue(&curStack);
      if (curStack.sign == SIGN_NEGATIVE)
      {
        return EXPR_INVALID_PARAM;
      }
      if (curStack.nbrLimbs > 1)
      {
        return EXPR_INTERM_TOO_HIGH;
      }
      retcode = primorial(&curStack, curStack.limbs[0].x);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case OPER_PLUS:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        stackIndex--;
        break;
      }
      getCurrentStackValue(&curStack2);
      stackIndex--;
      getCurrentStackValue(&curStack);
      BigIntAdd(&curStack, &curStack2, &curStack);
      break;

    case OPER_MINUS:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        stackIndex--;
        break;
      }
      getCurrentStackValue(&curStack2);
      stackIndex--;
      getCurrentStackValue(&curStack);
      BigIntSubt(&curStack, &curStack2, &curStack);
      break;

    case OPER_UNARY_MINUS:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        break;
      }
      getCurrentStackValue(&curStack);
      BigIntNegate(&curStack, &curStack);
      break;

    case OPER_DIVIDE:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        stackIndex--;
        break;
      }
      getCurrentStackValue(&curStack2);
      stackIndex--;
      getCurrentStackValue(&curStack);
      retcode = BigIntDivide(&curStack, &curStack2, &curStack);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case OPER_MULTIPLY:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        stackIndex--;
        break;
      }
      getCurrentStackValue(&curStack2);
      stackIndex--;
      getCurrentStackValue(&curStack);
      retcode = BigIntMultiply(&curStack, &curStack2, &curStack);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case OPER_REMAINDER:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        stackIndex--;
        break;
      }
      getCurrentStackValue(&curStack2);
      stackIndex--;
      getCurrentStackValue(&curStack);
      retcode = BigIntRemainder(&curStack, &curStack2, &curStack);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;
    case TOKEN_END_EXPON:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        stackIndex--;
        break;
      }
      getCurrentStackValue(&curStack2);
      stackIndex--;
      getCurrentStackValue(&curStack);
      retcode = BigIntPower(&curStack, &curStack2, &curStack3);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      CopyBigInt(&curStack, &curStack3);    // Copy result.
      break;

    case OPER_EQUAL:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        stackIndex--;
        break;
      }
      getCurrentStackValue(&curStack2);
      stackIndex--;
      getCurrentStackValue(&curStack);
      intToBigInteger(&curStack,
        (BigIntEqual(&curStack, &curStack2) ? -1 : 0));
      break;
    case OPER_NOT_EQUAL:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        stackIndex--;
        break;
      }
      getCurrentStackValue(&curStack2);
      stackIndex--;
      getCurrentStackValue(&curStack);
      intToBigInteger(&curStack,
        (BigIntEqual(&curStack, &curStack2) ? 0: -1));
      break;

    case OPER_GREATER:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        stackIndex--;
        break;
      }
      getCurrentStackValue(&curStack2);
      stackIndex--;
      getCurrentStackValue(&curStack);
      BigIntSubt(&curStack2, &curStack, &curStack3);
      intToBigInteger(&curStack, (curStack3.sign == SIGN_NEGATIVE) ? -1 : 0);
      break;

    case OPER_NOT_GREATER:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        stackIndex--;
        break;
      }
      getCurrentStackValue(&curStack2);
      stackIndex--;
      getCurrentStackValue(&curStack);
      BigIntSubt(&curStack2, &curStack, &curStack3);
      intToBigInteger(&curStack, (curStack3.sign == SIGN_NEGATIVE) ? 0 : -1);
      break;

    case OPER_LESS:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        stackIndex--;
        break;
      }
      getCurrentStackValue(&curStack2);
      stackIndex--;
      getCurrentStackValue(&curStack);
      BigIntSubt(&curStack, &curStack2, &curStack3);
      intToBigInteger(&curStack, (curStack3.sign == SIGN_NEGATIVE) ? -1 : 0);
      break;

    case OPER_NOT_LESS:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        stackIndex--;
        break;
      }
      getCurrentStackValue(&curStack2);
      stackIndex--;
      getCurrentStackValue(&curStack);
      BigIntSubt(&curStack, &curStack2, &curStack3);
      intToBigInteger(&curStack, (curStack3.sign == SIGN_NEGATIVE) ? 0 : -1);
      break;

    case OPER_SHL:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        stackIndex--;
        break;
      }
      getCurrentStackValue(&curStack2);
      stackIndex--;
      getCurrentStackValue(&curStack);
      retcode = ShiftLeft(&curStack, &curStack2, &curStack3);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      CopyBigInt(&curStack, &curStack3);    // Copy result.
      break;

    case OPER_SHR:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        stackIndex--;
        break;
      }
      getCurrentStackValue(&curStack2);
      stackIndex--;
      getCurrentStackValue(&curStack);
      BigIntChSign(&curStack2);
      retcode = ShiftLeft(&curStack, &curStack2, &curStack3);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      CopyBigInt(&curStack, &curStack3);    // Copy result.
      break;

    case OPER_NOT:    // Perform binary NOT as result <- -1 - argument.
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        break;
      }
      getCurrentStackValue(&curStack);
      intToBigInteger(&curStack2, -1);
      BigIntSubt(&curStack2, &curStack, &curStack);
      break;

    case OPER_INFIX_AND:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        break;
      }
      // If value in stack is zero, do not compute second operand.
      getCurrentStackValue(&curStack);
      if (BigIntIsZero(&curStack))
      {   // Result of AND is zero. Ignore all component of second operand.
        stackIndexThreshold = stackIndex;
      }
      break;

    case OPER_AND:    // Perform binary AND.
    if ((stackIndexThreshold+1) == stackIndex)
      {
        stackIndex--;
        stackIndexThreshold = DO_NOT_SHORT_CIRCUIT;
        break;
      }
      getCurrentStackValue(&curStack2);
      stackIndex--;
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        break;
      }
      getCurrentStackValue(&curStack);
      BigIntAnd(&curStack, &curStack2, &curStack3);
      CopyBigInt(&curStack, &curStack3);
      break;

    case OPER_INFIX_OR:
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        break;
      }
      // If value in stack is -1, do not compute second operand.
      getCurrentStackValue(&curStack);
      if ((curStack.sign == SIGN_NEGATIVE) && (curStack.nbrLimbs == 1) &&
        (curStack.limbs[0].x == 1))
      {     // Result of OR is -1. Ignore all component of second operand.
        stackIndexThreshold = stackIndex;
      }
      break;

    case OPER_OR:     // Perform binary OR.
      if ((stackIndexThreshold+1) == stackIndex)
      {
        stackIndex--;
        stackIndexThreshold = DO_NOT_SHORT_CIRCUIT;
        break;
      }
      getCurrentStackValue(&curStack2);
      stackIndex--;
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        break;
      }
      getCurrentStackValue(&curStack);
      BigIntOr(&curStack, &curStack2, &curStack3);
      CopyBigInt(&curStack, &curStack3);
      break;
    case OPER_XOR:    // Perform binary XOR.
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        stackIndex--;
        break;
      }
      getCurrentStackValue(&curStack2);
      stackIndex--;
      getCurrentStackValue(&curStack);
      BigIntXor(&curStack, &curStack2, &curStack3);
      CopyBigInt(&curStack, &curStack3);
      break;

    default:
      break;
    }
    if (c != TOKEN_NUMBER)
    {
      if ((stackIndexThreshold >= stackIndex) && (c != TOKEN_START_EXPON))
      {
        retcode = setStackValue(&curStack);
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
      }
      ptrRPNbuffer++;            // Skip token.
    }
  }
  getCurrentStackValue(ExpressionResult);
#ifdef FACTORIZATION_APP
  if (ExpressionResult->nbrLimbs > ((332192 / BITS_PER_GROUP) + 1))   // 100000/log_10(2) = 332192
#else
  if (ExpressionResult->nbrLimbs > ((33219 / BITS_PER_GROUP) + 1))    // 10000/log_10(2) = 33219
#endif
  {
    return EXPR_NUMBER_TOO_HIGH;
  }
  if ((valueX.nbrLimbs > 0) && !valueXused)
  {
    return EXPR_VAR_OR_COUNTER_REQUIRED;
  }
  return EXPR_OK;
}

static enum eExprErr ComputeModInv(void)
{
  static BigInteger one;
  const BigInteger* pDiv;
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

static enum eExprErr ComputePartition(void)
{
  const BigInteger *pArgument = &curStack;
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
    largeVal.x = pArgument->limbs[0].x + 
      (int)((unsigned int)pArgument->limbs[1].x << BITS_PER_GROUP);
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
  partition(val, &curStack);
  return EXPR_OK;
}

#ifdef FACTORIZATION_FUNCTIONS
static void PerformFactorization(const BigInteger* numToFactor)
{
  char* ptrFactorDec = tofactorDec;
  NumberLength = numToFactor->nbrLimbs;
  BigInteger2IntArray(nbrToFactor, numToFactor);
  if (hexadecimal)
  {
    Bin2Hex(&ptrFactorDec, numToFactor->limbs, numToFactor->nbrLimbs, groupLen);
  }
  else
  {
    Bin2Dec(&ptrFactorDec, numToFactor->limbs, numToFactor->nbrLimbs, groupLen);
  }
  factor(numToFactor, nbrToFactor, factorsMod, astFactorsMod);
}

static enum eExprErr ComputeTotient(void)
{
  getCurrentStackValue(&curStack);    // Get argument.
  PerformFactorization(&curStack);
  Totient(&curStack);
  return EXPR_OK;
}

static enum eExprErr ComputeNumDivs(void)
{
  getCurrentStackValue(&curStack);    // Get argument.
  PerformFactorization(&curStack);
  NumberOfDivisors(&curStack);
  return EXPR_OK;
}

static enum eExprErr ComputeSumDivs(void)
{
  getCurrentStackValue(&curStack);    // Get argument.
  PerformFactorization(&curStack);
  SumOfDivisors(&curStack);
  return EXPR_OK;
}

static enum eExprErr ComputeNumFact(void)
{
  getCurrentStackValue(&curStack);    // Get argument.
  PerformFactorization(&curStack);
  NumFactors(&curStack);
  return EXPR_OK;
}

static enum eExprErr ComputeMaxFact(void)
{
  getCurrentStackValue(&curStack);    // Get argument.
  PerformFactorization(&curStack);
  MaxFactor(&curStack);
  return EXPR_OK;
}

static enum eExprErr ComputeMinFact(void)
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
    primeIndex++;
    prime = smallPrimes[primeIndex];
  } while (prime < 32767);
  PerformFactorization(&curStack);
  MinFactor(&curStack);
  return EXPR_OK;
}

static enum eExprErr ComputeConcatFact(void)
{
  getCurrentStackValue(&curStack);    // Get first argument.
  stackIndex++;
  getCurrentStackValue(&curStack2);   // Get second argument.
  stackIndex--;
  const BigInteger *mode = &curStack;
  static BigInteger factorValue;
  int nbrFactors;
  bool descend = ((mode->limbs[0].x & 1)? true: false);
  bool repeated = ((mode->limbs[0].x & 2)? true: false);
  char *ptrTextFactor = textFactor;
  if ((mode->nbrLimbs > 1) || (mode->sign == SIGN_NEGATIVE) || (mode->limbs[0].x > 3))
  {      // The valid modes are 0, 1, 2 and 3.
    return EXPR_INVALID_PARAM;
  }
  PerformFactorization(&curStack2);   // Factor second argument.
  nbrFactors = astFactorsMod[0].multiplicity;
  for (int factorNumber = 1; factorNumber <= nbrFactors; factorNumber++)
  {
    int ctr;
    const struct sFactors* pstFactor;
    if (descend)
    {
      pstFactor = &astFactorsMod[nbrFactors - factorNumber + 1];
    }
    else
    {
      pstFactor = &astFactorsMod[factorNumber];
    }
    NumberLength = *(pstFactor->ptrFactor);
    IntArray2BigInteger(pstFactor->ptrFactor, &factorValue);
    ctr = (repeated ? pstFactor->multiplicity : 1);
    for (; ctr > 0; ctr--)
    {
      BigInteger2Dec(&ptrTextFactor, &factorValue, 0);
    }
  }
  if (ptrTextFactor == &textFactor[0])
  {
    intToBigInteger(&curStack, 0);
  }
  else
  {
    size_t diffPtrs = ptrTextFactor - &textFactor[0];
    int size = (int)diffPtrs;
    Dec2Bin(textFactor, curStack.limbs, size, &curStack.nbrLimbs);
  }
  return EXPR_OK;
}
#endif

static enum eExprErr ComputeSumDigits(void)
{
  getCurrentStackValue(&curStack);    // Get first argument.
  stackIndex++;
  getCurrentStackValue(&curStack2);   // Get second argument.
  stackIndex--;
  static BigInteger argum;
  static BigInteger Temp;
  BigInteger *result = &curStack;
  const BigInteger *radix = &curStack2;
  CopyBigInt(&argum, &curStack);
  intToBigInteger(result, 0);
  while (!BigIntIsZero(&argum))
  {
    (void)BigIntRemainder(&argum, radix, &Temp);
    BigIntAdd(result, &Temp, result);
    (void)BigIntDivide(&argum, radix, &argum);
  }
  return EXPR_OK;
}

static enum eExprErr ComputeNumDigits(void)
{
  getCurrentStackValue(&curStack);    // Get first argument.
  stackIndex++;
  getCurrentStackValue(&curStack2);   // Get second argument.
  stackIndex--;
  BigInteger *result = &curStack;
  const BigInteger *radix = &curStack2;
  int digits = 0;
  while (!BigIntIsZero(result))
  {
    (void)BigIntDivide(result, radix, result);
    digits++;
  }
  intToBigInteger(result, digits);
  return EXPR_OK;
}

static enum eExprErr ComputeRevDigits(void)
{
  getCurrentStackValue(&curStack);    // Get first argument.
  stackIndex++;
  getCurrentStackValue(&curStack2);   // Get second argument.
  stackIndex--;
  static BigInteger argum;
  static BigInteger Temp;
  BigInteger *result = &curStack;
  const BigInteger *radix = &curStack2;
  CopyBigInt(&argum, &curStack);
  intToBigInteger(result, 0);
  while (!BigIntIsZero(&argum))
  {
    (void)BigIntRemainder(&argum, radix, &Temp);
    (void)BigIntMultiply(result, radix, result);
    BigIntAdd(result, &Temp, result);
    (void)BigIntDivide(&argum, radix, &argum);
  }
  return EXPR_OK;
}

static enum eExprErr ShiftLeft(BigInteger* first, const BigInteger *second, BigInteger *result)
{
  int ctr;
  unsigned int prevLimb;
  unsigned int curLimb;
  unsigned int shRight;
  unsigned int shLeft;
  int *ptrDest;
  const int *ptrSrc;
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
    if (((nbrLimbs * BITS_PER_GROUP) + shiftCtr) > 664380)
#else
    if (((nbrLimbs * BITS_PER_GROUP) + shiftCtr) > 66438)
#endif
    {   // Shift too much to the left.
      return EXPR_INTERM_TOO_HIGH;
    }
    prevLimb = 0U;
    ptrSrc = &first->limbs[nbrLimbs - 1].x;
    curLimb = (unsigned int)*ptrSrc;
    ptrDest = &result->limbs[nbrLimbs + delta].x;
    shLeft = (unsigned int)rem;
    shRight = (unsigned int)BITS_PER_GROUP - shLeft;
    for (ctr = nbrLimbs; ctr > 0; ctr--)
    {  // Process starting from most significant limb.
      *ptrDest = UintToInt(((curLimb >> shRight) | (prevLimb << shLeft)) & MAX_VALUE_LIMB);
      ptrDest--;
      prevLimb = curLimb;
      ptrSrc--;
      curLimb = (unsigned int)*ptrSrc;
    }
    *ptrDest = UintToInt(((curLimb >> shRight) | (prevLimb << shLeft)) & MAX_VALUE_LIMB);
    if (delta > 0)
    {
      int lenBytes = delta * (int)sizeof(limb);
      (void)memset(first->limbs, 0, lenBytes);
    }
    result->nbrLimbs = first->nbrLimbs + delta;
    if (result->limbs[result->nbrLimbs].x != 0)
    {
      result->nbrLimbs++;
    }
  }
  else
  {     // Perform shift right.
    bool isNegative = false;
    if ((second->nbrLimbs > 1) || (shiftCtr > (first->nbrLimbs * BITS_PER_GROUP)))
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
      isNegative = true;
      addbigint(first, 1);
    }
    // Shift right the absolute value.
    first->limbs[nbrLimbs].x = 0;
    ptrSrc = &first->limbs[delta+1].x;
    prevLimb = (unsigned int)*(ptrSrc-1);
    curLimb = (unsigned int)*ptrSrc;
    ptrDest = &result->limbs[0].x;
    shLeft = (unsigned int)rem;
    shRight = (unsigned int)BITS_PER_GROUP - shLeft;
    for (ctr = delta; ctr < nbrLimbs; ctr++)
    {  // Process starting from least significant limb.
      *ptrDest = UintToInt(((prevLimb >> shLeft) | (curLimb << shRight)) & MAX_VALUE_LIMB);
      ptrDest++;
      prevLimb = curLimb;
      ptrSrc++;
      curLimb = (unsigned int)*ptrSrc;
    }
    *ptrDest = UintToInt(((prevLimb >> shLeft) | (curLimb << shRight)) & MAX_VALUE_LIMB);
    result->nbrLimbs = first->nbrLimbs - (delta + 1);
    if ((result->nbrLimbs == 0) || (result->limbs[result->nbrLimbs].x))
    {
      result->nbrLimbs++;
    }
    result->sign = first->sign;
    if ((result->nbrLimbs == 1) && (result->limbs[0].x == 0))
    {    // Result is zero.
      result->sign = SIGN_POSITIVE;
    }
    if (isNegative)
    {    // Adjust negative number.
      addbigint(result, -1);
    }
  }
  return EXPR_OK;
}

#ifndef __EMSCRIPTEN__
void databack(const char *data)
{
  (void)data;
}
#endif

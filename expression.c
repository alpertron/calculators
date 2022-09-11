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
#include "exprfact.h"
#include "showtime.h"
#ifdef USING_BLOCKLY
#include "batch.h"
#include "fromBlockly.h"
extern bool doFactorization;
extern char* ptrBlocklyOutput;
extern BigInteger tofactor;
extern int nbrBlocklyOutputLines;
extern bool doShowPrime;
#endif

#define PAREN_STACK_SIZE           5000
#define COMPR_STACK_SIZE        1000000

#define DO_NOT_SHORT_CIRCUIT  (COMPR_STACK_SIZE + 100)  // Larger than stack size.

BigInteger LastAnswer;
const struct sFuncOperExpr stFuncOperIntExpr[] =
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
  {"ANS", TOKEN_ANS + NO_PARMS, 0},
  {"GCD", TOKEN_GCD + MANY_PARMS, 0},
  {"LCM", TOKEN_LCM + MANY_PARMS, 0},
  {"MODPOW", TOKEN_MODPOW + THREE_PARMS, 0},
  {"MODINV", TOKEN_MODINV + TWO_PARMS, 0},
  {"IROOT", TOKEN_IROOT + TWO_PARMS, 0},
  {"SUMDIGITS", TOKEN_SUMDIGITS + TWO_PARMS, 0},
  {"NUMDIGITS", TOKEN_NUMDIGITS + TWO_PARMS, 0},
  {"REVDIGITS", TOKEN_REVDIGITS + TWO_PARMS, 0},
  {"ISPRIME", TOKEN_ISPRIME + ONE_PARM, 0},
  {"JACOBI", TOKEN_JACOBI + TWO_PARMS, 0},
  {"RANDOM", TOKEN_RANDOM + TWO_PARMS, 0},
  {"SQRT", TOKEN_SQRT + ONE_PARM, 0},
  {"ABS", TOKEN_ABS + ONE_PARM, 0},
  {"SGN", TOKEN_SGN + ONE_PARM, 0},
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
  {"+", OPER_ADD, 4},
  {"-", OPER_SUBT, 4},
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
int valueQ[MAX_LEN];
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
static int ComputeRandom(void);
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

static enum eExprErr getCurrentStackValue(BigInteger* pValue)
{
  if ((stackIndex < 0) || (stackIndex >= PAREN_STACK_SIZE))
  {
    return EXPR_CANNOT_PARSE_EXPRESSION;
  }
  const int* ptrStackValue = &comprStackValues[comprStackOffset[stackIndex]];
  NumberLength = numLimbs(ptrStackValue);
  IntArray2BigInteger(ptrStackValue, pValue);
  return EXPR_OK;
}

static enum eExprErr setStackValue(const BigInteger* pValue)
{
  if ((stackIndex < 0) || (stackIndex >= PAREN_STACK_SIZE))
  {
    return EXPR_CANNOT_PARSE_EXPRESSION;
  }
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

static enum eExprErr PerformGCDorLCM(char token, int stackIndexThreshold, 
  int nbrParameters)
{
  if (nbrParameters <= 0)
  {
    return EXPR_CANNOT_PARSE_EXPRESSION;
  }
  if (stackIndexThreshold < stackIndex)
  {     // Delete from operand stack when binary AND/OR short-circuited.
    stackIndex -= nbrParameters - 1;
    return EXPR_OK;
  }
  if (getCurrentStackValue(&curStack) != EXPR_OK)
  {
    return EXPR_CANNOT_PARSE_EXPRESSION;
  }
  for (int parmNbr = 1; parmNbr < nbrParameters; parmNbr++)
  {
    stackIndex--;
    if (getCurrentStackValue(&curStack2) != EXPR_OK)
    {
      return EXPR_CANNOT_PARSE_EXPRESSION;
    }
    if (token == TOKEN_LCM)
    {
      return BigIntLcm(&curStack, &curStack2, &curStack);
    }
    BigIntGcd(&curStack, &curStack2, &curStack);
  }
  return EXPR_OK;
}

static void GetLastAnswer(int stackIndexThreshold)
{
  stackIndex++;
  if (stackIndexThreshold >= stackIndex)
  {     // Part of second operand of binary AND/OR not short-circuited.
    if (LastAnswer.nbrLimbs == 0)
    {
      intToBigInteger(&curStack, 0);
    }
    else
    {
      CopyBigInt(&curStack, &LastAnswer);
    }
  }
}

static enum eExprErr getParms(int nbrParms, int stackIndexThreshold)
{
  if (stackIndexThreshold < stackIndex)
  {     // Part of second operand of binary AND/OR short-circuited.
    stackIndex -= nbrParms - 1;
    return EXPR_SHORT_CIRCUIT;
  }
  if (nbrParms == 3)
  {
    if (getCurrentStackValue(&curStack3) != EXPR_OK)
    {
      return EXPR_CANNOT_PARSE_EXPRESSION;
    }
    stackIndex--;
  }
  if (nbrParms >= 2)
  {
    if (getCurrentStackValue(&curStack2) != EXPR_OK)
    {
      return EXPR_CANNOT_PARSE_EXPRESSION;
    }
    stackIndex--;
  }
  if (getCurrentStackValue(&curStack) != EXPR_OK)
  {
    return EXPR_CANNOT_PARSE_EXPRESSION;
  }
  return EXPR_OK;
}

enum eExprErr ComputeExpression(const char *expr, BigInteger *ExpressionResult,
  bool varsExpected)
{
  bool valueXused;
  bool randomUsed;
  enum eExprErr retcode;
  char* pointerRPNbuffer;
  const char* ptrRPNbuffer;
  int len;
#ifdef USING_BLOCKLY
  const char* ptrRPNstartBuffer;
  bool hexBak;
  bool doFactBak;
  int offset;
#endif
  int stackIndexThreshold = DO_NOT_SHORT_CIRCUIT;
#ifdef USING_BLOCKLY
  if (ExpressionResult != NULL)
  {
#endif
    retcode = ConvertToReversePolishNotation(expr, &pointerRPNbuffer, stFuncOperIntExpr,
      PARSE_EXPR_INTEGER, &valueXused, &randomUsed);
    if (retcode != EXPR_OK)
    {
      return retcode;
    }
    if (varsExpected)
    {      // Variable or random function expected.
      if (randomUsed)
      {
        valueXused = true;
      }
    }
    else
    {
      if (valueXused)
      {    // Variable present is an error here.
        return EXPR_VAR_IN_EXPRESSION;
      }
    }
    ptrRPNbuffer = pointerRPNbuffer;
#ifdef USING_BLOCKLY
  }
  else
  {
    ptrRPNbuffer = expr;
    valueXused = false;
  }
  ptrRPNstartBuffer = ptrRPNbuffer;
#endif
  stackIndex = -1;
  comprStackOffset[0] = 0;
  while (*ptrRPNbuffer != '\0')
  {
    char c = *ptrRPNbuffer;
    int currentOffset;
    int nbrLenBytes;
    int jacobi;
    int Exponent;
    switch (c)
    {
    case TOKEN_NUMBER:
      ptrRPNbuffer++;           // Skip token.
      stackIndex++;
      len = ((int)(unsigned char)*ptrRPNbuffer * 256) +
        (int)(unsigned char)*(ptrRPNbuffer + 1);
      nbrLenBytes = len * (int)sizeof(limb);
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        ptrRPNbuffer += nbrLenBytes;
        ptrRPNbuffer += 2;
        break;
      }
      // Move number to compressed stack.
      if (stackIndex < 0)
      {
        return EXPR_CANNOT_PARSE_EXPRESSION;
      }
      currentOffset = comprStackOffset[stackIndex];
      if (currentOffset >= (COMPR_STACK_SIZE -
        ((int)sizeof(BigInteger) / (int)sizeof(limb))))
      {
        return EXPR_OUT_OF_MEMORY;
      }
      comprStackValues[currentOffset] = len;
      ptrRPNbuffer += 2;   // Skip length.
      (void)memcpy(&comprStackValues[currentOffset + 1], ptrRPNbuffer, nbrLenBytes);
      ptrRPNbuffer += nbrLenBytes;
      comprStackOffset[stackIndex + 1] = currentOffset + 1 + len;
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

    case TOKEN_ANS:
      GetLastAnswer(stackIndexThreshold);
      break;

    case TOKEN_GCD:
    case TOKEN_LCM:
      ptrRPNbuffer++;
      retcode = PerformGCDorLCM(c, stackIndexThreshold,
        (int)(unsigned char)*ptrRPNbuffer);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_JACOBI:
      retcode = getParms(2, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      jacobi = BigIntJacobiSymbol(&curStack, &curStack2);
      intToBigInteger(&curStack, jacobi);
      break;

    case TOKEN_IROOT:
      retcode = getParms(2, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      if ((curStack2.sign == SIGN_NEGATIVE) || BigIntIsZero(&curStack2))
      {
        return EXPR_INVALID_PARAM;
      }
      if (!BigIntIsZero(&curStack))
      {
        Exponent = curStack2.limbs[0].x;
        if (curStack2.nbrLimbs > 1)
        {
          if (curStack.sign == SIGN_POSITIVE)
          {
            intToBigInteger(&curStack, 1);
          }
          else if ((Exponent % 2) == 0)
          {
            return EXPR_BASE_MUST_BE_POSITIVE;
          }
          else if ((curStack.nbrLimbs == 1) && (curStack.limbs[0].x == 1))
          {       // First parameters is -1.
            intToBigInteger(&curStack, -1);
          }
          else
          {
            intToBigInteger(&curStack, -2);
          }
        }
        else
        {
          retcode = BigIntRoot(&curStack, &curStack2, Exponent);
          if (retcode != EXPR_OK)
          {
            return retcode;
          }
          CopyBigInt(&curStack, &curStack2);
        }
      }
      break;

    case TOKEN_ABS:
      retcode = getParms(1, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      curStack.sign = SIGN_POSITIVE;
      break;

    case TOKEN_SGN:
      retcode = getParms(1, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      if (!BigIntIsZero(&curStack))
      {     // If less than zero, return -1. If greater than zero, return 1.
        curStack.nbrLimbs = 1;
        curStack.limbs[0].x = 1;
      }
      break;

    case TOKEN_SQRT:
      retcode = getParms(1, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      if (curStack.sign == SIGN_NEGATIVE)
      {
        return EXPR_INVALID_PARAM;
      }
      squareRoot(curStack.limbs, curStack2.limbs, curStack.nbrLimbs, &curStack2.nbrLimbs);
      curStack2.sign = SIGN_POSITIVE;
      CopyBigInt(&curStack, &curStack2);
      break;

    case TOKEN_MODPOW:
      retcode = getParms(3, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = BigIntGeneralModularPower(&curStack, &curStack2, &curStack3, &curStack);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_MODINV:
      retcode = getParms(2, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = ComputeModInv();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

#ifdef FACTORIZATION_FUNCTIONS
    case TOKEN_TOTIENT:
      retcode = getParms(1, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = ComputeTotient();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_NUMDIVS:
      retcode = getParms(1, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = ComputeNumDivs();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_SUMDIVS:
      retcode = getParms(1, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = ComputeSumDivs();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_MINFACT:
      retcode = getParms(1, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = ComputeMinFact();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_MAXFACT:
      retcode = getParms(1, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = ComputeMaxFact();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_NUMFACT:
      retcode = getParms(1, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = ComputeNumFact();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_CONCATFACT:
      retcode = getParms(2, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = ComputeConcatFact();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;
#endif

    case TOKEN_SUMDIGITS:
      retcode = getParms(2, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = ComputeSumDigits();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_RANDOM:
      retcode = getParms(2, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = ComputeRandom();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_NUMDIGITS:
      retcode = getParms(2, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = ComputeNumDigits();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_REVDIGITS:
      retcode = getParms(2, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = ComputeRevDigits();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_ISPRIME:
      retcode = getParms(1, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
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
      retcode = getParms(1, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = ComputeFibLucas(0, &curStack);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_L:
      retcode = getParms(1, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = ComputeFibLucas(2, &curStack);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_P:
      retcode = getParms(1, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = ComputePartition();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_N:
      retcode = getParms(1, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = ComputeNext(&curStack);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_B:
      retcode = getParms(1, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = ComputeBack(&curStack);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_FACTORIAL:
      retcode = getParms(2, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
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
      retcode = getParms(1, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
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

    case OPER_ADD:
      retcode = getParms(2, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      BigIntAdd(&curStack, &curStack2, &curStack);
      break;

    case OPER_SUBT:
      retcode = getParms(2, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      BigIntSubt(&curStack, &curStack2, &curStack);
      break;

    case OPER_UNARY_MINUS:
      retcode = getParms(1, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      BigIntNegate(&curStack, &curStack);
      break;

    case OPER_DIVIDE:
      retcode = getParms(2, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = BigIntDivide(&curStack, &curStack2, &curStack);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case OPER_MULTIPLY:
      retcode = getParms(2, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = BigIntMultiply(&curStack, &curStack2, &curStack);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case OPER_REMAINDER:
      retcode = getParms(2, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = BigIntRemainder(&curStack, &curStack2, &curStack);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;
    case TOKEN_END_EXPON:
    case OPER_POWER:       // Coming here from Blockly.
      retcode = getParms(2, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = BigIntPower(&curStack, &curStack2, &curStack3);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      CopyBigInt(&curStack, &curStack3);    // Copy result.
      break;

    case OPER_EQUAL:
      retcode = getParms(2, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      intToBigInteger(&curStack,
        (BigIntEqual(&curStack, &curStack2) ? -1 : 0));
      break;
    case OPER_NOT_EQUAL:
      retcode = getParms(2, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      intToBigInteger(&curStack,
        (BigIntEqual(&curStack, &curStack2) ? 0 : -1));
      break;

    case OPER_GREATER:
      retcode = getParms(2, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      BigIntSubt(&curStack2, &curStack, &curStack3);
      intToBigInteger(&curStack, (curStack3.sign == SIGN_NEGATIVE) ? -1 : 0);
      break;

    case OPER_NOT_GREATER:
      retcode = getParms(2, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      BigIntSubt(&curStack2, &curStack, &curStack3);
      intToBigInteger(&curStack, (curStack3.sign == SIGN_NEGATIVE) ? 0 : -1);
      break;

    case OPER_LESS:
      retcode = getParms(2, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      BigIntSubt(&curStack, &curStack2, &curStack3);
      intToBigInteger(&curStack, (curStack3.sign == SIGN_NEGATIVE) ? -1 : 0);
      break;

    case OPER_NOT_LESS:
      retcode = getParms(2, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      BigIntSubt(&curStack, &curStack2, &curStack3);
      intToBigInteger(&curStack, (curStack3.sign == SIGN_NEGATIVE) ? 0 : -1);
      break;

    case OPER_SHL:
      retcode = getParms(2, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = ShiftLeft(&curStack, &curStack2, &curStack3);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      CopyBigInt(&curStack, &curStack3);    // Copy result.
      break;

    case OPER_SHR:
      retcode = getParms(2, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      BigIntChSign(&curStack2);
      retcode = ShiftLeft(&curStack, &curStack2, &curStack3);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      CopyBigInt(&curStack, &curStack3);    // Copy result.
      break;

    case OPER_NOT:    // Perform binary NOT as result <- -1 - argument.
      retcode = getParms(1, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      intToBigInteger(&curStack2, -1);
      BigIntSubt(&curStack2, &curStack, &curStack);
      break;

    case OPER_INFIX_AND:
      retcode = getParms(1, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      // If value in stack is zero, do not compute second operand.
      if (BigIntIsZero(&curStack))
      {   // Result of AND is zero. Ignore all component of second operand.
        stackIndexThreshold = stackIndex;
      }
      break;

    case OPER_AND:    // Perform binary AND.
      if ((stackIndexThreshold + 1) == stackIndex)
      {
        stackIndex--;
        stackIndexThreshold = DO_NOT_SHORT_CIRCUIT;
        break;
      }
      if (getCurrentStackValue(&curStack2) != EXPR_OK)
      {
        return EXPR_CANNOT_PARSE_EXPRESSION;
      }
      stackIndex--;
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        break;
      }
      if (getCurrentStackValue(&curStack) != EXPR_OK)
      {
        return EXPR_CANNOT_PARSE_EXPRESSION;
      }
      BigIntAnd(&curStack, &curStack2, &curStack3);
      CopyBigInt(&curStack, &curStack3);
      break;

    case OPER_INFIX_OR:
      retcode = getParms(1, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      // If value in stack is -1, do not compute second operand.
      if ((curStack.sign == SIGN_NEGATIVE) && (curStack.nbrLimbs == 1) &&
        (curStack.limbs[0].x == 1))
      {     // Result of OR is -1. Ignore all component of second operand.
        stackIndexThreshold = stackIndex;
      }
      break;

    case OPER_OR:     // Perform binary OR.
      if ((stackIndexThreshold + 1) == stackIndex)
      {
        stackIndex--;
        stackIndexThreshold = DO_NOT_SHORT_CIRCUIT;
        break;
      }
      if (getCurrentStackValue(&curStack2) != EXPR_OK)
      {
        return EXPR_CANNOT_PARSE_EXPRESSION;
      }
      stackIndex--;
      if (stackIndexThreshold < stackIndex)
      {     // Part of second operand of binary AND/OR short-circuited.
        break;
      }
      if (getCurrentStackValue(&curStack) != EXPR_OK)
      {
        return EXPR_CANNOT_PARSE_EXPRESSION;
      }
      BigIntOr(&curStack, &curStack2, &curStack3);
      CopyBigInt(&curStack, &curStack3);
      break;
    case OPER_XOR:    // Perform binary XOR.
      retcode = getParms(2, stackIndexThreshold);
      if (retcode == EXPR_SHORT_CIRCUIT)
      {
        break;
      }
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      BigIntXor(&curStack, &curStack2, &curStack3);
      CopyBigInt(&curStack, &curStack3);
      break;

#ifdef USING_BLOCKLY
    case TOKEN_GET_VAR:
      ptrRPNbuffer++;
      stackIndex++;
      getBlocklyVar((unsigned char)*ptrRPNbuffer, &curStack);
      break;

    case TOKEN_SET_VAR:
      if (getCurrentStackValue(&curStack) != EXPR_OK)
      {
        return EXPR_CANNOT_PARSE_EXPRESSION;
      }
      ptrRPNbuffer++;
      setBlocklyVar((unsigned char)*ptrRPNbuffer, &curStack);
      stackIndex--;
      break;

    case TOKEN_JMP:
      offset = ((int)(unsigned char)*(ptrRPNbuffer + 1)) +
        ((int)(unsigned char)*(ptrRPNbuffer + 2)) * 0x100 +
        ((int)(unsigned char)*(ptrRPNbuffer + 3)) * 0x10000 +
        ((int)(unsigned char)*(ptrRPNbuffer + 4)) * 0x1000000;
      ptrRPNbuffer = ptrRPNstartBuffer + offset - 1;
      break;

    case TOKEN_IF:
      if (getCurrentStackValue(&curStack) != EXPR_OK)
      {
        return EXPR_CANNOT_PARSE_EXPRESSION;
      }
      if (BigIntIsZero(&curStack))
      {
        offset = ((int)(unsigned char)*(ptrRPNbuffer + 1)) +
          ((int)(unsigned char)*(ptrRPNbuffer + 2)) * 0x100 +
          ((int)(unsigned char)*(ptrRPNbuffer + 3)) * 0x10000 +
          ((int)(unsigned char)*(ptrRPNbuffer + 4)) * 0x1000000;
        ptrRPNbuffer = ptrRPNstartBuffer + offset - 1;
      }
      else
      {      // Skip jump.
        ptrRPNbuffer += 4;
      }
      stackIndex--;
      break;

    case TOKEN_PRINT:
    case TOKEN_PRINT_HEX:
    case TOKEN_PRINTFACT:
    case TOKEN_PRINTFACT_HEX:
    case TOKEN_PRINTPRIME:
    case TOKEN_PRINTPRIME_HEX:
      hexBak = hexadecimal;
      doFactBak = doFactorization;
      hexadecimal = ((c == TOKEN_PRINT_HEX) || (c == TOKEN_PRINTFACT_HEX) || (c == TOKEN_PRINTPRIME_HEX));
      doFactorization = ((c == TOKEN_PRINTFACT) || (c == TOKEN_PRINTFACT_HEX));
      doShowPrime = ((c == TOKEN_PRINTPRIME) || (c == TOKEN_PRINTPRIME_HEX));
      if (getCurrentStackValue(&tofactor) != EXPR_OK)
      {
        return EXPR_CANNOT_PARSE_EXPRESSION;
      }
      stackIndex--;
      copyStr(&ptrBlocklyOutput, "<li>");
      batchEcmCallback(&ptrBlocklyOutput);
      copyStr(&ptrBlocklyOutput, "</li>");
      nbrBlocklyOutputLines++;
      if (nbrBlocklyOutputLines == 1000)
      {
        copyStr(&ptrBlocklyOutput, "</ul>");
        output[0] = '6';  // Show Continue button.
        databack(output);
      }
      hexadecimal = hexBak;
      doFactorization = doFactBak;
      break;
#endif
    default:
      break;
    }
    if (c != TOKEN_NUMBER)
    {
#ifdef USING_BLOCKLY
      if ((c != TOKEN_SET_VAR) && (c != TOKEN_JMP) && (c != TOKEN_IF) &&
        (c != TOKEN_PRINT) && (c != TOKEN_PRINTFACT) &&
        (c != TOKEN_PRINT_HEX) && (c != TOKEN_PRINTFACT_HEX))
#endif
      {
        if ((stackIndexThreshold >= stackIndex) && (c != TOKEN_START_EXPON))
        {         
          retcode = setStackValue(&curStack);
          if (retcode != EXPR_OK)
          {
            return retcode;
          }
        }
      }
      ptrRPNbuffer++;            // Skip token.
    }
  }
  if (ExpressionResult != NULL)
  {
    retcode = getCurrentStackValue(ExpressionResult);
    if (retcode != EXPR_OK)
    {
      return retcode;
    }
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
    CopyBigInt(&LastAnswer, ExpressionResult);
  }
  else
  {
    return getCurrentStackValue(&LastAnswer);
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
  PerformFactorization(&curStack);
  Totient(&curStack);
  return EXPR_OK;
}

static enum eExprErr ComputeNumDivs(void)
{
  PerformFactorization(&curStack);
  NumberOfDivisors(&curStack);
  return EXPR_OK;
}

static enum eExprErr ComputeSumDivs(void)
{
  PerformFactorization(&curStack);
  SumOfDivisors(&curStack);
  return EXPR_OK;
}

static enum eExprErr ComputeNumFact(void)
{
  PerformFactorization(&curStack);
  NumFactors(&curStack);
  return EXPR_OK;
}

static enum eExprErr ComputeMaxFact(void)
{
  PerformFactorization(&curStack);
  MaxFactor(&curStack);
  return EXPR_OK;
}

static enum eExprErr ComputeMinFact(void)
{
  int primeIndex = 0;
  int prime = 2;
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
  if ((curStack2.nbrLimbs == 1) && (curStack2.limbs[0].x <= 1))
  {   // If number is zero or one, they have no prime factors, so result is zero.
    curStack.limbs[0].x = 0;
    return EXPR_OK;
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

static enum eExprErr ComputeRandom(void)
{
  BigIntRandom(&curStack, &curStack2, &curStack);
  return EXPR_OK;
}

static enum eExprErr ComputeNumDigits(void)
{
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

// Perform arithmetic shift left: sign of result must be the sign of
// value, and mantissa of result must be the mantissa of value
// shifted left shiftCtr bits.
static enum eExprErr PerformShiftLeft(BigInteger* value, int shiftCtr,
  BigInteger* result)
{
  unsigned int prevLimb;
  unsigned int curLimb;
  unsigned int shRight;
  unsigned int shLeft;
  limb* ptrDest;
  const limb* ptrSrc;
  int delta = shiftCtr / BITS_PER_GROUP;
  int rem = shiftCtr % BITS_PER_GROUP;
  int nbrLimbs = value->nbrLimbs;
#ifdef FACTORIZATION_APP
  if (((nbrLimbs * BITS_PER_GROUP) + shiftCtr) > 664380)
#else
  if (((nbrLimbs * BITS_PER_GROUP) + shiftCtr) > 66438)
#endif
  {   // Shift too much to the left.
    return EXPR_INTERM_TOO_HIGH;
  }
  ptrDest = &result->limbs[nbrLimbs + delta];
  shLeft = (unsigned int)rem;
  shRight = (unsigned int)BITS_PER_GROUP - shLeft;
  ptrSrc = &value->limbs[nbrLimbs - 1];
  prevLimb = 0U;
  while (ptrSrc >= &value->limbs[0])
  {  // Process starting from most significant limb.
    curLimb = (unsigned int)ptrSrc->x;
    ptrDest->x = UintToInt(((curLimb >> shRight) | (prevLimb << shLeft)) & MAX_VALUE_LIMB);
    ptrDest--;
    prevLimb = curLimb;
    ptrSrc--;
  }
  ptrDest->x = UintToInt((prevLimb << shLeft) & MAX_VALUE_LIMB);
  if (delta > 0)
  {
    int lenBytes = delta * (int)sizeof(limb);
    (void)memset(value->limbs, 0, lenBytes);
  }
  result->nbrLimbs = value->nbrLimbs + delta;
  if (result->limbs[result->nbrLimbs].x != 0)
  {
    result->nbrLimbs++;
  }
  return EXPR_OK;
}

// Perform arithmetic shift right: sign of result must be the sign of
// value, and mantissa of result depends on sign:
// If sign is positive, the mantissa of result is the mantissa of value
// shifted right shiftCtr bits.
// If sign is negative, add 1, perform shift right shiftCtr bits, and finally
// subtract 1 to result.
static enum eExprErr PerformShiftRight(BigInteger* value, int shiftCtr,
  BigInteger* result)
{
  unsigned int prevLimb;
  unsigned int curLimb;
  unsigned int shRight;
  unsigned int shLeft;
  limb* ptrDest;
  const limb* ptrSrc;
  int delta = shiftCtr / BITS_PER_GROUP;
  int rem = shiftCtr % BITS_PER_GROUP;
  int nbrLimbs = value->nbrLimbs;
  bool isNegative = false;
  if (value->sign == SIGN_NEGATIVE)
  {   // If it is negative, add 1, perform shift right, and finally subtract 1 to result.
    isNegative = true;
    addbigint(value, 1);
  }
  // Shift right the absolute value.
  value->limbs[nbrLimbs].x = 0;
  ptrSrc = &value->limbs[delta + 1];
  prevLimb = (unsigned int)(ptrSrc - 1)->x;
  curLimb = (unsigned int)ptrSrc->x;
  ptrDest = &result->limbs[0];
  shLeft = (unsigned int)rem;
  shRight = (unsigned int)BITS_PER_GROUP - shLeft;
  for (int ctr = delta; ctr < nbrLimbs; ctr++)
  {  // Process starting from least significant limb.
    ptrDest->x = UintToInt(((prevLimb >> shLeft) | (curLimb << shRight)) & MAX_VALUE_LIMB);
    ptrDest++;
    prevLimb = curLimb;
    ptrSrc++;
    curLimb = (unsigned int)ptrSrc->x;
  }
  ptrDest->x = UintToInt(((prevLimb >> shLeft) | (curLimb << shRight)) & MAX_VALUE_LIMB);
  result->nbrLimbs = value->nbrLimbs - (delta + 1);
  if ((result->nbrLimbs == 0) || (result->limbs[result->nbrLimbs].x))
  {
    result->nbrLimbs++;
  }
  result->sign = value->sign;
  if ((result->nbrLimbs == 1) && (result->limbs[0].x == 0))
  {    // Result is zero.
    result->sign = SIGN_POSITIVE;
  }
  if (isNegative)
  {    // Adjust negative number.
    addbigint(result, -1);
  }
  return EXPR_OK;
}

static enum eExprErr ShiftLeft(BigInteger* first, const BigInteger *second,
  BigInteger *result)
{
  int shiftCtr = second->limbs[0].x;
  if (second->sign == SIGN_POSITIVE)
  {     // Perform shift left.
    if (second->nbrLimbs > 1)
    {   // Shift too much to the left.
      return EXPR_INTERM_TOO_HIGH;
    }
    return PerformShiftLeft(first, shiftCtr, result);
  }
  // Perform shift right.
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
  return PerformShiftRight(first, shiftCtr, result);
}

#ifndef __EMSCRIPTEN__
void databack(const char *data)
{
  (void)data;
}
#endif

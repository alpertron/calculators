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
#include <assert.h>
#include "bignbr.h"
#include "factor.h"
#include "expression.h"
#include "commonstruc.h"

#define TOKEN_PRIMORIAL  39
#define TOKEN_GCD        40
#define TOKEN_LCM        41
#define TOKEN_MODPOW     42
#define TOKEN_MODINV     43
#define TOKEN_NORM       44
#define TOKEN_RE         45
#define TOKEN_IM         46
#define TOKEN_ISPRIME    47
#define TOKEN_F          48
#define TOKEN_L          49
#define TOKEN_P          50
#define TOKEN_N          51
#define TOKEN_B          52

struct sFuncOperExpr stFuncOperGaussianExpr[] =
{
  // First section: functions
  {"ANS", TOKEN_ANS + NO_PARMS, 0},
  {"GCD", TOKEN_GCD + MANY_PARMS, 0},
  {"LCM", TOKEN_LCM + MANY_PARMS, 0},
  {"MODPOW", TOKEN_MODPOW + THREE_PARMS, 0},
  {"MODINV", TOKEN_MODINV + TWO_PARMS, 0},
  {"NORM", TOKEN_NORM + ONE_PARM, 0},
  {"RE", TOKEN_RE + ONE_PARM, 0},
  {"IM", TOKEN_IM + ONE_PARM, 0},
  {"ISPRIME", TOKEN_ISPRIME + ONE_PARM, 0},
  {"F", TOKEN_F + ONE_PARM, 0},
  {"L", TOKEN_L + ONE_PARM, 0},
  {"P", TOKEN_P + ONE_PARM, 0},
  {"N", TOKEN_N + ONE_PARM, 0},
  {"B", TOKEN_B + ONE_PARM, 0},
  {"RANDOM", TOKEN_RANDOM + TWO_PARMS, 0},
  {NULL, 0, 0},
  // Second section: functions written at right of argument.
  {"#", TOKEN_PRIMORIAL, 0},
  {NULL, 0, 0},
  // Third section: unary operators.
  {"-", OPER_UNARY_MINUS, 3},
  {NULL, 0, 0},
  // Fourth section: binary operators.
  {"**", OPER_POWER, 1}, // This must be located before multiplication operator.
  {"+", OPER_ADD, 4},
  {"-", OPER_SUBT, 4},
  {"*", OPER_MULTIPLY, 2},
  {"%", OPER_REMAINDER, 2},
  {"/", OPER_DIVIDE, 2},
  {"^", OPER_POWER, 1},
  {NULL, 0, 0},
};

extern int comprStackValues[COMPR_STACK_SIZE];
extern int comprStackOffset[PAREN_STACK_SIZE];
extern int stackIndex;
#ifdef __ANDROID__
  extern limb Mult1[MAX_LEN];
  extern limb Mult3[MAX_LEN];
  extern limb Mult4[MAX_LEN];
  extern int valueQ[MAX_LEN];
#else
  char output[3000000];
  limb Mult1[MAX_LEN];
  limb Mult3[MAX_LEN];
  limb Mult4[MAX_LEN];
  int valueQ[MAX_LEN];
#endif
static void ComputeGCD(void);
static enum eExprErr ComputeLCM(void);
static enum eExprErr GaussianDivide(void);
static enum eExprErr GaussianMultiply(void);
static int Modulo(BigInteger* ReNum, BigInteger* ImNum,
  const BigInteger* ReDen, const BigInteger* ImDen,
  BigInteger* Result);
static enum eExprErr isPrime(void);
static enum eExprErr ComputeModPow(void);
static enum eExprErr ComputeModInv(void);
static enum eExprErr ComputePartition(void);
static enum eExprErr ComputePower(BigInteger* Re1, const BigInteger* Re2,
  BigInteger* Im1, const BigInteger* Im2);
static enum eExprErr ModInv(const BigInteger* RealNbr, const BigInteger* ImagNbr,
  const BigInteger* RealMod, const BigInteger* ImagMod,
  BigInteger* Result);

static int numLimbs(const int* pLen)
{
  int nbrLimbs = *pLen;
  if (nbrLimbs < 0)
  {
    nbrLimbs = -nbrLimbs;
  }
  return nbrLimbs;
}

static void getCurrentStackValue(BigInteger* pValueRe, BigInteger *pValueIm)
{
  assert(stackIndex >= 0);
  const int* ptrStackValue = &comprStackValues[comprStackOffset[2*stackIndex]];
  NumberLength = numLimbs(ptrStackValue);
  IntArray2BigInteger(ptrStackValue, pValueRe);
  ptrStackValue = &comprStackValues[comprStackOffset[(2 * stackIndex) + 1]];
  NumberLength = numLimbs(ptrStackValue);
  IntArray2BigInteger(ptrStackValue, pValueIm);
}

static enum eExprErr setStackValue(const BigInteger* pValueRe, const BigInteger *pValueIm)
{
  assert(stackIndex >= 0);
  int currentOffset = comprStackOffset[2 * stackIndex];
  if (currentOffset >= (COMPR_STACK_SIZE - ((int)sizeof(BigInteger) / (int)sizeof(limb))))
  {
    return EXPR_OUT_OF_MEMORY;
  }
  NumberLength = pValueRe->nbrLimbs;
  BigInteger2IntArray(&comprStackValues[currentOffset], pValueRe);
  currentOffset += pValueRe->nbrLimbs + 1;
  comprStackOffset[(2 * stackIndex) + 1] = currentOffset;
  NumberLength = pValueIm->nbrLimbs;
  BigInteger2IntArray(&comprStackValues[currentOffset], pValueIm);
  currentOffset += pValueIm->nbrLimbs + 1;
  comprStackOffset[(2 * stackIndex) + 2] = currentOffset;
  return EXPR_OK;
}

enum eExprErr ComputeGaussianExpression(const char *expr, BigInteger *ExpressionResult)
{
  char* ptrRPNbuffer;
  enum eExprErr retcode;
  retcode = ConvertToReversePolishNotation(expr, &ptrRPNbuffer, stFuncOperGaussianExpr,
    PARSE_EXPR_GAUSSIAN, NULL, NULL);
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
    int len;
    int nbrParameters;
    switch (c)
    {
    case TOKEN_NUMBER:
      ptrRPNbuffer++;           // Skip token.
      stackIndex++;
      // Move number to compressed stack.
      currentOffset = comprStackOffset[2 * stackIndex];
      if (currentOffset >= (COMPR_STACK_SIZE - ((int)sizeof(BigInteger) / (int)sizeof(limb))))
      {
        return EXPR_OUT_OF_MEMORY;
      }
      len = ((int)(unsigned char)*ptrRPNbuffer * 256) + (int)(unsigned char)*(ptrRPNbuffer + 1);
      comprStackValues[currentOffset] = len;
      ptrRPNbuffer += 2;   // Skip length.
      nbrLenBytes = len * (int)sizeof(limb);
      (void)memcpy(&comprStackValues[currentOffset + 1], ptrRPNbuffer, nbrLenBytes);
      ptrRPNbuffer += nbrLenBytes;
      currentOffset += 1 + len;
      comprStackOffset[(2 * stackIndex) + 1] = currentOffset;
      comprStackValues[currentOffset] = 1;      // Imaginary part is zero.
      comprStackValues[currentOffset + 1] = 0;
      comprStackOffset[(2 * stackIndex) + 2] = currentOffset + 2;
      break;

    case TOKEN_VAR:    // Push number i.
      ptrRPNbuffer++;                             // Skip token.
      stackIndex++;
      currentOffset = comprStackOffset[2 * stackIndex];
      comprStackValues[currentOffset] = 1;      // Real part is zero.
      comprStackValues[currentOffset + 1] = 0;
      currentOffset += 2;
      comprStackOffset[(2 * stackIndex) + 1] = currentOffset;
      comprStackValues[currentOffset] = 1;      // Imaginary part is one.
      comprStackValues[currentOffset + 1] = 1;
      comprStackOffset[(2 * stackIndex) + 2] = currentOffset + 2;
      break;

    case TOKEN_ANS:    // Get last answer.
      if (common.gaussian.LastAnswerRe.nbrLimbs == 0)
      {
        intToBigInteger(&common.gaussian.curStackRe, 0);
      }
      else
      {
        CopyBigInt(&common.gaussian.curStackRe, &common.gaussian.LastAnswerRe);
      }
      if (common.gaussian.LastAnswerIm.nbrLimbs == 0)
      {
        intToBigInteger(&common.gaussian.curStackIm, 0);
      }
      else
      {
        CopyBigInt(&common.gaussian.curStackIm, &common.gaussian.LastAnswerIm);
      }
      stackIndex++;
      break;

    case TOKEN_RE:
      if (stackIndex < 0)
      {
        return EXPR_INTERNAL_ERROR;
      }
      getCurrentStackValue(&common.gaussian.curStackRe, &common.gaussian.curStackIm);
      intToBigInteger(&common.gaussian.curStackIm, 0);
      break;

    case TOKEN_IM:
      if (stackIndex < 0)
      {
        return EXPR_INTERNAL_ERROR;
      }
      getCurrentStackValue(&common.gaussian.curStackIm, &common.gaussian.curStackRe);
      intToBigInteger(&common.gaussian.curStackIm, 0);
      break;

    case TOKEN_NORM:
      if (stackIndex < 0)
      {
        return EXPR_INTERNAL_ERROR;
      }
      // norm = Re^2 + Im^2
      getCurrentStackValue(&common.gaussian.curStackRe, &common.gaussian.curStackIm);
      retcode = BigIntMultiply(&common.gaussian.curStackRe, &common.gaussian.curStackRe,
                               &common.gaussian.curStack2Re);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = BigIntMultiply(&common.gaussian.curStackIm, &common.gaussian.curStackIm,
                               &common.gaussian.curStack2Im);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      BigIntAdd(&common.gaussian.curStack2Re, &common.gaussian.curStack2Im, &common.gaussian.curStackRe);
      intToBigInteger(&common.gaussian.curStackIm, 0);
      break;

    case TOKEN_GCD:
      ptrRPNbuffer++;
      nbrParameters = (int)(unsigned char)*ptrRPNbuffer;
      getCurrentStackValue(&common.gaussian.curStackRe, &common.gaussian.curStackIm);
      for (int parmNbr = 1; parmNbr < nbrParameters; parmNbr++)
      {
        stackIndex--;
        getCurrentStackValue(&common.gaussian.curStack2Re, &common.gaussian.curStack2Im);
        ComputeGCD();
      }
      break;

    case TOKEN_LCM:
      ptrRPNbuffer++;
      nbrParameters = (int)(unsigned char)*ptrRPNbuffer;
      getCurrentStackValue(&common.gaussian.curStackRe, &common.gaussian.curStackIm);
      for (int parmNbr = 1; parmNbr < nbrParameters; parmNbr++)
      {
        stackIndex--;
        getCurrentStackValue(&common.gaussian.curStack2Re, &common.gaussian.curStack2Im);
        retcode = ComputeLCM();
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
      }
      break;

    case TOKEN_MODPOW:
      getCurrentStackValue(&common.gaussian.curStack3Re, &common.gaussian.curStack3Im);
      stackIndex--;
      getCurrentStackValue(&common.gaussian.curStack2Re, &common.gaussian.curStack2Im);
      stackIndex--;
      getCurrentStackValue(&common.gaussian.curStackRe, &common.gaussian.curStackIm);
      retcode = ComputeModPow();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_MODINV:
      getCurrentStackValue(&common.gaussian.curStack2Re, &common.gaussian.curStack2Im);
      stackIndex--;
      getCurrentStackValue(&common.gaussian.curStackRe, &common.gaussian.curStackIm);
      retcode = ComputeModInv();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_ISPRIME:
      getCurrentStackValue(&common.gaussian.curStackRe, &common.gaussian.curStackIm);
      retcode = isPrime();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_F:
      getCurrentStackValue(&common.gaussian.curStackRe, &common.gaussian.curStackIm);
      if (!BigIntIsZero(&common.gaussian.curStackIm))
      {     // Imaginary part of argument must be zero.
        return EXPR_INVALID_PARAM;
      }
      retcode = ComputeFibLucas(0, &common.gaussian.curStackRe);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_L:
      getCurrentStackValue(&common.gaussian.curStackRe, &common.gaussian.curStackIm);
      if (!BigIntIsZero(&common.gaussian.curStackIm))
      {     // Imaginary part of argument must be zero.
        return EXPR_INVALID_PARAM;
      }
      retcode = ComputeFibLucas(2, &common.gaussian.curStackRe);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_N:
      getCurrentStackValue(&common.gaussian.curStackRe, &common.gaussian.curStackIm);
      if (!BigIntIsZero(&common.gaussian.curStackIm))
      {     // Imaginary part of argument must be zero.
        return EXPR_INVALID_PARAM;
      }
      retcode = ComputeNext(&common.gaussian.curStackRe);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_B:
      getCurrentStackValue(&common.gaussian.curStackRe, &common.gaussian.curStackIm);
      if (!BigIntIsZero(&common.gaussian.curStackIm))
      {     // Imaginary part of argument must be zero.
        return EXPR_INVALID_PARAM;
      }
      retcode = ComputeBack(&common.gaussian.curStackRe);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_P:
      getCurrentStackValue(&common.gaussian.curStackRe, &common.gaussian.curStackIm);
      retcode = ComputePartition();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_RANDOM:
      getCurrentStackValue(&common.gaussian.curStack2Re, &common.gaussian.curStack2Im);
      stackIndex--;
      getCurrentStackValue(&common.gaussian.curStackRe, &common.gaussian.curStackIm);
      BigIntRandom(&common.gaussian.curStackRe, &common.gaussian.curStack2Re, &common.gaussian.curStackRe);
      BigIntRandom(&common.gaussian.curStackIm, &common.gaussian.curStack2Im, &common.gaussian.curStackIm);
      break;

    case TOKEN_FACTORIAL:
      getCurrentStackValue(&common.gaussian.curStack2Re, &common.gaussian.curStack2Im);
      stackIndex--;
      getCurrentStackValue(&common.gaussian.curStackRe, &common.gaussian.curStackIm);
      if (!BigIntIsZero(&common.gaussian.curStackIm) || (common.gaussian.curStackRe.sign == SIGN_NEGATIVE))
      {         // Imaginary part must be zero.
        return EXPR_INVALID_PARAM;
      }
      retcode = factorial(&common.gaussian.curStackRe, common.gaussian.curStackRe.limbs[0].x,
                          common.gaussian.curStack2Re.limbs[0].x);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_PRIMORIAL:
      getCurrentStackValue(&common.gaussian.curStackRe, &common.gaussian.curStackIm);
      if (!BigIntIsZero(&common.gaussian.curStackIm) || (common.gaussian.curStackRe.sign == SIGN_NEGATIVE))
      {         // Imaginary part must be zero.
        return EXPR_INVALID_PARAM;
      }
      if ((common.gaussian.curStackRe.nbrLimbs > 1) || (common.gaussian.curStackRe.limbs[0].x > 46049))
      {
        return EXPR_INTERM_TOO_HIGH;
      }
      retcode = primorial(&common.gaussian.curStackRe, common.gaussian.curStackRe.limbs[0].x);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case OPER_ADD:
      getCurrentStackValue(&common.gaussian.curStack2Re, &common.gaussian.curStack2Im);
      stackIndex--;
      getCurrentStackValue(&common.gaussian.curStackRe, &common.gaussian.curStackIm);
      BigIntAdd(&common.gaussian.curStackRe, &common.gaussian.curStack2Re, &common.gaussian.curStackRe);
      BigIntAdd(&common.gaussian.curStackIm, &common.gaussian.curStack2Im, &common.gaussian.curStackIm);
      break;

    case OPER_SUBT:
      getCurrentStackValue(&common.gaussian.curStack2Re, &common.gaussian.curStack2Im);
      stackIndex--;
      getCurrentStackValue(&common.gaussian.curStackRe, &common.gaussian.curStackIm);
      BigIntSubt(&common.gaussian.curStackRe, &common.gaussian.curStack2Re, &common.gaussian.curStackRe);
      BigIntSubt(&common.gaussian.curStackIm, &common.gaussian.curStack2Im, &common.gaussian.curStackIm);
      break;

    case OPER_UNARY_MINUS:
      getCurrentStackValue(&common.gaussian.curStackRe, &common.gaussian.curStackIm);
      BigIntChSign(&common.gaussian.curStackRe);
      BigIntChSign(&common.gaussian.curStackIm);
      break;

    case OPER_MULTIPLY:
      getCurrentStackValue(&common.gaussian.curStack2Re, &common.gaussian.curStack2Im);
      stackIndex--;
      getCurrentStackValue(&common.gaussian.curStackRe, &common.gaussian.curStackIm);
      retcode = GaussianMultiply();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case OPER_DIVIDE:
      getCurrentStackValue(&common.gaussian.curStack2Re, &common.gaussian.curStack2Im);
      stackIndex--;
      getCurrentStackValue(&common.gaussian.curStackRe, &common.gaussian.curStackIm);
      retcode = GaussianDivide();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case OPER_REMAINDER:
      getCurrentStackValue(&common.gaussian.curStack2Re, &common.gaussian.curStack2Im);
      stackIndex--;
      getCurrentStackValue(&common.gaussian.curStackRe, &common.gaussian.curStackIm);
      retcode = Modulo(&common.gaussian.curStackRe, &common.gaussian.curStackIm,
        &common.gaussian.curStack2Re, &common.gaussian.curStack2Im, common.gaussian.Result);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      CopyBigInt(&common.gaussian.curStackRe, &common.gaussian.Result[0]);
      CopyBigInt(&common.gaussian.curStackIm, &common.gaussian.Result[1]);
      break;

    case TOKEN_END_EXPON:
      getCurrentStackValue(&common.gaussian.curStack2Re, &common.gaussian.curStack2Im);
      stackIndex--;
      getCurrentStackValue(&common.gaussian.curStackRe, &common.gaussian.curStackIm);
      retcode = ComputePower(&common.gaussian.curStackRe, &common.gaussian.curStack2Re,
        &common.gaussian.curStackIm, &common.gaussian.curStack2Im);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    default:
      break;
    }
    if ((c != TOKEN_NUMBER) && (c != TOKEN_VAR))
    {
      if (c != TOKEN_START_EXPON)
      {
        retcode = setStackValue(&common.gaussian.curStackRe, &common.gaussian.curStackIm);
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
      }
      ptrRPNbuffer++;            // Skip token.
    }
  }
  getCurrentStackValue(ExpressionResult, ExpressionResult + 1);
  if ((ExpressionResult->nbrLimbs > 2215) &&    // 10000/log_10(32768)
      ((ExpressionResult+1)->nbrLimbs > 2215))
  {
    return EXPR_NUMBER_TOO_HIGH;
  }
  CopyBigInt(&common.gaussian.LastAnswerRe, ExpressionResult);
  CopyBigInt(&common.gaussian.LastAnswerIm, ExpressionResult + 1);
  return EXPR_OK;
}

static enum eExprErr isPrime(void)
{
  if (BigIntIsZero(&common.gaussian.curStackRe))
  {
    if ((common.gaussian.curStackIm.limbs[0].x % 4) != 3)
    {
      intToBigInteger(&common.gaussian.curStackRe, 0);    // Number not prime.
      intToBigInteger(&common.gaussian.curStackIm, 0);
      return EXPR_OK;
    }
    CopyBigInt(&common.gaussian.curStack2Re, &common.gaussian.curStackIm);
  }
  else if (BigIntIsZero(&common.gaussian.curStackIm))
  {
    if ((common.gaussian.curStackRe.limbs[0].x % 4) != 3)
    {
      intToBigInteger(&common.gaussian.curStackRe, 0);    // Number not prime.
      intToBigInteger(&common.gaussian.curStackIm, 0);
      return EXPR_OK;
    }
    CopyBigInt(&common.gaussian.curStack2Re, &common.gaussian.curStackRe);
  }
  else
  {
    enum eExprErr retcode;
    retcode = BigIntMultiply(&common.gaussian.curStackRe, &common.gaussian.curStackRe,
                             &common.gaussian.curStack2Re);
    if (retcode != EXPR_OK)
    {
      return retcode;
    }
    retcode = BigIntMultiply(&common.gaussian.curStackIm, &common.gaussian.curStackIm,
                             &common.gaussian.curStack2Im);
    if (retcode != EXPR_OK)
    {
      return retcode;
    }
    BigIntAdd(&common.gaussian.curStack2Re, &common.gaussian.curStack2Im, &common.gaussian.curStack2Re);
  }
#ifdef FACTORIZATION_APP
  if (BpswPrimalityTest(&common.gaussian.curStack2Re, NULL) == 0)
#else
  if (BpswPrimalityTest(&common.gaussian.curStack2Re) == 0)
#endif
  {    // Argument is a probable prime.
    intToBigInteger(&common.gaussian.curStackRe, -1);
  }
  else
  {    // Argument is not a probable prime.
    intToBigInteger(&common.gaussian.curStackRe, 0);
  }
  intToBigInteger(&common.gaussian.curStackIm, 0);
  return EXPR_OK;
}

static enum eExprErr ComputePartition(void)
{
  int val;

  if (!BigIntIsZero(&common.gaussian.curStackIm))
  {         // Imaginary part must be zero.
    return EXPR_INVALID_PARAM;
  }
  if (common.gaussian.curStackRe.sign == SIGN_NEGATIVE)
  {
    return EXPR_INVALID_PARAM;
  }
  if (common.gaussian.curStackRe.nbrLimbs > 1)
  {
    return EXPR_INTERM_TOO_HIGH;
  }
  val = common.gaussian.curStackRe.limbs[0].x;
  if (val > 100000)
  {
    return EXPR_INVALID_PARAM;
  }
  partition(val, &common.gaussian.curStackRe);
  return EXPR_OK;
}

// Replace dividend by the remainder.
// Back-up sign of dividend and then divide positive by positive (divisor is always positive).
static void GetRemainder(const BigInteger *Norm, BigInteger *ReDividend, BigInteger *ImDividend,
  const BigInteger *ReDivisor, const BigInteger *ImDivisor,
  BigInteger *Re, BigInteger *Im)
{
  BigInteger ReTmp;
  BigInteger ImTmp;
  int signBak;
  // Compute Re as ((Re1*Re2+Im1*Im2)*2/norm+1)/2
  (void)BigIntMultiply(ReDividend, ReDivisor, &ReTmp);
  (void)BigIntMultiply(ImDividend, ImDivisor, &ImTmp);
  BigIntAdd(&ReTmp, &ImTmp, &ReTmp);
  multint(&ReTmp, &ReTmp, 2);
  signBak = ReTmp.sign;
  ReTmp.sign = SIGN_POSITIVE;
  (void)BigIntDivide(&ReTmp, Norm, Re);
  subtractdivide(Re, -1, 2);
  if (!BigIntIsZero(Re))
  {
    Re->sign = signBak;
  }
  // Compute Im as ((Im1*Re2-Re1*Im2)*2/norm+1)/2
  (void)BigIntMultiply(ImDividend, ReDivisor, &ReTmp);
  (void)BigIntMultiply(ReDividend, ImDivisor, &ImTmp);
  BigIntSubt(&ReTmp, &ImTmp, &ImTmp);
  multint(&ImTmp, &ImTmp, 2);
  signBak = ImTmp.sign;
  ImTmp.sign = SIGN_POSITIVE;
  (void)BigIntDivide(&ImTmp, Norm, Im);
  subtractdivide(Im, -1, 2);
  if (!BigIntIsZero(Im))
  {
    Im->sign = signBak;
  }
  // Compute Re1 as Re1 - Re*Re2 + Im*Im2.
  (void)BigIntMultiply(Re, ReDivisor, &ReTmp);
  BigIntSubt(ReDividend, &ReTmp, ReDividend);
  (void)BigIntMultiply(Im, ImDivisor, &ReTmp);
  BigIntAdd(ReDividend, &ReTmp, ReDividend);
  // Compute Im1 as Im1 - Im*Re2 - Re*Im2.
  (void)BigIntMultiply(Im, ReDivisor, &ReTmp);
  BigIntSubt(ImDividend, &ReTmp, ImDividend);
  (void)BigIntMultiply(Re, ImDivisor, &ReTmp);
  BigIntSubt(ImDividend, &ReTmp, ImDividend);
}

// Compute curStack as gcd(curStack, curStack2)
static void ComputeGCD(void)
{
  BigInteger ReTmp;
  BigInteger ImTmp;
  GaussianGCD(&common.gaussian.curStackRe, &common.gaussian.curStackIm,
    &common.gaussian.curStack2Re, &common.gaussian.curStack2Im,
    &common.gaussian.curStack3Re, &common.gaussian.curStack3Im, &ReTmp, &ImTmp);
  CopyBigInt(&common.gaussian.curStackRe, &common.gaussian.curStack3Re);
  CopyBigInt(&common.gaussian.curStackIm, &common.gaussian.curStack3Im);
  while ((common.gaussian.curStackRe.sign == SIGN_NEGATIVE) || (common.gaussian.curStackIm.sign == SIGN_NEGATIVE))
  {    // Multiply by i until number is in first quadrant.
    CopyBigInt(&common.gaussian.curStack3Re, &common.gaussian.curStackRe);
    BigIntNegate(&common.gaussian.curStackIm, &common.gaussian.curStackRe);
    CopyBigInt(&common.gaussian.curStackIm, &common.gaussian.curStack3Re);
  }
}

static enum eExprErr ComputeLCM(void)
{
  enum eExprErr retcode;
  if ((BigIntIsZero(&common.gaussian.curStackRe) && BigIntIsZero(&common.gaussian.curStackIm)) ||
    (BigIntIsZero(&common.gaussian.curStack2Re) && BigIntIsZero(&common.gaussian.curStack2Im)))
  {    // If any of the arguments is zero, the LCM is zero.
    intToBigInteger(&common.gaussian.curStackRe, 0);
    intToBigInteger(&common.gaussian.curStackIm, 0);
    return EXPR_OK;
  }
  // Save both parameters because ComputeGCD will overwrite them.
  CopyBigInt(&common.gaussian.curStackReBak, &common.gaussian.curStackRe);
  CopyBigInt(&common.gaussian.curStackImBak, &common.gaussian.curStackIm);
  CopyBigInt(&common.gaussian.curStack2ReBak, &common.gaussian.curStack2Re);
  CopyBigInt(&common.gaussian.curStack2ImBak, &common.gaussian.curStack2Im);
  ComputeGCD();   // curStackRe + i*curStackIm = GCD.
  // Divide first parameter by GCD.
  CopyBigInt(&common.gaussian.curStack2Re, &common.gaussian.curStackRe);    // Divisor
  CopyBigInt(&common.gaussian.curStack2Im, &common.gaussian.curStackIm);
  CopyBigInt(&common.gaussian.curStackRe, &common.gaussian.curStackReBak);  // Dividend
  CopyBigInt(&common.gaussian.curStackIm, &common.gaussian.curStackImBak);
  retcode = GaussianDivide();
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  // Multiply the quotient by the second parameter.
  CopyBigInt(&common.gaussian.curStack2Re, &common.gaussian.curStack2ReBak);
  CopyBigInt(&common.gaussian.curStack2Im, &common.gaussian.curStack2ImBak);
  return GaussianMultiply();
}

static int ComputePower(BigInteger *Re1, const BigInteger *Re2,
  BigInteger *Im1, const BigInteger *Im2)
{
  unsigned int expon;
  double logNorm;
  bool performPower = false;
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
  BigIntAdd(&ReTmp, &ImTmp, &common.gaussian.norm);  // norm <- re1^2 + im1^2.
  logNorm = logLimbs(common.gaussian.norm.limbs, common.gaussian.norm.nbrLimbs);  // Get logarithm of norm.
  if (logNorm*(double)expon > 23026.0)
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
  for (unsigned int mask = HALF_INT_RANGE_U; mask != 0U; mask >>= 1)
  {
    if ((expon & mask) != 0U)
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
      if ((expon & mask) != 0U)
      {
        (void)BigIntMultiply(Re1, &Re, &ReTmp);       // Re2 <- re1*re - im1*im.
        (void)BigIntMultiply(Im1, &Im, &ImTmp);
        BigIntSubt(&ReTmp, &ImTmp, &common.gaussian.curTmp);
        (void)BigIntMultiply(Re1, &Im, &ReTmp);       // Im <- re1*im + im1*re.
        (void)BigIntMultiply(Im1, &Re, &ImTmp);
        BigIntAdd(&ReTmp, &ImTmp, &Im);
        CopyBigInt(&Re, &common.gaussian.curTmp);
      }
    }
  }
  CopyBigInt(Re1, &Re);
  CopyBigInt(Im1, &Im);
  return EXPR_OK;
}

static enum eExprErr ComputeModPow(void)
{
  CopyBigInt(&common.gaussian.ReBase, &common.gaussian.curStackRe);
  CopyBigInt(&common.gaussian.ImBase, &common.gaussian.curStackIm);
  CopyBigInt(&common.gaussian.ReExp, &common.gaussian.curStack2Re);
  CopyBigInt(&common.gaussian.ImExp, &common.gaussian.curStack2Im);
  CopyBigInt(&common.gaussian.ReMod, &common.gaussian.curStack3Re);
  CopyBigInt(&common.gaussian.ImMod, &common.gaussian.curStack3Im);

  if (!BigIntIsZero(&common.gaussian.ImExp))
  {         // Imaginary part must be zero.
    return EXPR_INVALID_PARAM;
  }
  if (common.gaussian.ReExp.sign == SIGN_NEGATIVE)
  {
    enum eExprErr retcode;
    BigIntNegate(&common.gaussian.ReExp, &common.gaussian.ReExp);
    retcode = ModInv(&common.gaussian.ReBase, &common.gaussian.ImBase,
                     &common.gaussian.ReMod, &common.gaussian.ImMod, common.gaussian.Result);
    if (retcode != EXPR_OK)
    {
      return retcode;
    }
    common.gaussian.ReBase = common.gaussian.Result[0];
    common.gaussian.ImBase = common.gaussian.Result[1];
  }
  // Re <- 1, Im <- 0
  common.gaussian.Re.limbs[0].x = 1;
  common.gaussian.Im.limbs[0].x = 0;
  common.gaussian.Re.nbrLimbs = 1;
  common.gaussian.Im.nbrLimbs = 1;
  common.gaussian.Re.sign = SIGN_POSITIVE;
  common.gaussian.Im.sign = SIGN_POSITIVE;
  if (BigIntIsZero(&common.gaussian.ReMod) && BigIntIsZero(&common.gaussian.ImMod))
  {   /* Modulus is zero */
    (void)ComputePower(&common.gaussian.ReBase, &common.gaussian.ReExp,
                       &common.gaussian.ImBase, &common.gaussian.ImExp);
    CopyBigInt(&common.gaussian.curStackRe, &common.gaussian.ReBase);
    CopyBigInt(&common.gaussian.curStackIm, &common.gaussian.ImBase);
  }
  else
  {                            /* Modulus is not zero */
    (void)BigIntMultiply(&common.gaussian.ReMod, &common.gaussian.ReMod, &common.gaussian.ReTmp);
    (void)BigIntMultiply(&common.gaussian.ImMod, &common.gaussian.ImMod, &common.gaussian.ImTmp);
    BigIntAdd(&common.gaussian.ReTmp, &common.gaussian.ImTmp, &common.gaussian.norm);
    for (int index = common.gaussian.ReExp.nbrLimbs - 1; index >= 0; index--)
    {
      int groupExp = common.gaussian.ReExp.limbs[index].x;
      for (unsigned int mask = HALF_INT_RANGE_U; mask > 0U; mask >>= 1)
      {
        // Compute Re + i*Im as (Re + i*Im)^2
        (void)BigIntMultiply(&common.gaussian.Re, &common.gaussian.Re, &common.gaussian.ReTmp);
        (void)BigIntMultiply(&common.gaussian.Im, &common.gaussian.Im, &common.gaussian.ImTmp);
        BigIntSubt(&common.gaussian.ReTmp, &common.gaussian.ImTmp, &common.gaussian.ReTmp);
        (void)BigIntMultiply(&common.gaussian.Re, &common.gaussian.Im, &common.gaussian.ImTmp);
        BigIntAdd(&common.gaussian.ImTmp, &common.gaussian.ImTmp, &common.gaussian.Im);
        CopyBigInt(&common.gaussian.Re, &common.gaussian.ReTmp);
        // Replace (Re + i*Im) by the remainder of
        // (ReTmp + i*ImTmp) / (ReMod + i*ImMod)
        // ReTmp and ImTmp are used as temporary values.
        GetRemainder(&common.gaussian.norm, &common.gaussian.Re, &common.gaussian.Im, &common.gaussian.ReMod,
                     &common.gaussian.ImMod, &common.gaussian.ReTmp, &common.gaussian.ImTmp);
        if (((unsigned int)groupExp & mask) != 0U)
        {
          // Compute Re + i*Im as (Re + i*Im)*(ReBase + i*ImBase)
          (void)BigIntMultiply(&common.gaussian.ReBase, &common.gaussian.Re, &common.gaussian.ReTmp);
          (void)BigIntMultiply(&common.gaussian.ImBase, &common.gaussian.Im, &common.gaussian.ImTmp);
          BigIntSubt(&common.gaussian.ReTmp, &common.gaussian.ImTmp, &common.gaussian.ReTmp);
          (void)BigIntMultiply(&common.gaussian.ImBase, &common.gaussian.Re, &common.gaussian.ImTmp);
          (void)BigIntMultiply(&common.gaussian.ReBase, &common.gaussian.Im, &common.gaussian.Re);
          BigIntAdd(&common.gaussian.Re, &common.gaussian.ImTmp, &common.gaussian.Im);
          CopyBigInt(&common.gaussian.Re, &common.gaussian.ReTmp);
          // Replace (Re + i*Im) by the remainder of
          // (ReTmp + i*ImTmp) / (ReMod + i*ImMod)
          // ReTmp and ImTmp are used as temporary values.
          GetRemainder(&common.gaussian.norm, &common.gaussian.Re, &common.gaussian.Im,
            &common.gaussian.ReMod, &common.gaussian.ImMod, &common.gaussian.ReTmp, &common.gaussian.ImTmp);
        }
      }
    }
    (void)Modulo(&common.gaussian.Re, &common.gaussian.Im, &common.gaussian.ReMod,
                 &common.gaussian.ImMod, common.gaussian.Result);
    CopyBigInt(&common.gaussian.curStackRe, &common.gaussian.Result[0]);
    CopyBigInt(&common.gaussian.curStackIm, &common.gaussian.Result[1]);
    return EXPR_OK;
  }
  return EXPR_OK;
}

static enum eExprErr ComputeModInv(void)
{
  enum eExprErr retcode = ModInv(&common.gaussian.curStackRe, &common.gaussian.curStackIm,
    &common.gaussian.curStack2Re, &common.gaussian.curStack2Im, common.gaussian.Result);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  CopyBigInt(&common.gaussian.curStackRe, &common.gaussian.Result[0]);
  CopyBigInt(&common.gaussian.curStackIm, &common.gaussian.Result[1]);
  return EXPR_OK;
}

static enum eExprErr ModInv(const BigInteger *RealNbr, const BigInteger *ImagNbr,
                  const BigInteger *RealMod, const BigInteger *ImagMod,
                  BigInteger *result)
{
  if (BigIntIsZero(RealMod) && BigIntIsZero(ImagMod))
  {          // Argument is zero.
    return EXPR_INVALID_PARAM;
  }
  CopyBigInt(&common.gaussian.ReG0, RealNbr);
  CopyBigInt(&common.gaussian.ImG0, ImagNbr);
  CopyBigInt(&common.gaussian.ReG1, RealMod);
  CopyBigInt(&common.gaussian.ImG1, ImagMod);
   // Initialize U0 <- 1, U1 <- 0.
  common.gaussian.ReU0.limbs[0].x = 1;
  common.gaussian.ReU0.nbrLimbs = 1;
  common.gaussian.ReU0.sign = SIGN_POSITIVE;
  common.gaussian.ReU1.limbs[0].x = 0;
  common.gaussian.ImU0.limbs[0].x = 0;
  common.gaussian.ImU1.limbs[0].x = 0;
  common.gaussian.ReU1.nbrLimbs = 1;
  common.gaussian.ImU0.nbrLimbs = 1;
  common.gaussian.ImU1.nbrLimbs = 1;
  common.gaussian.ReU1.sign = SIGN_POSITIVE;
  common.gaussian.ImU0.sign = SIGN_POSITIVE;
  common.gaussian.ImU1.sign = SIGN_POSITIVE;
  while (!BigIntIsZero(&common.gaussian.ReG1) || !BigIntIsZero(&common.gaussian.ImG1))
  {            // G1 is not zero.
    (void)BigIntMultiply(&common.gaussian.ReG1, &common.gaussian.ReG1, &common.gaussian.Re);
    (void)BigIntMultiply(&common.gaussian.ImG1, &common.gaussian.ImG1, &common.gaussian.Im);
    BigIntAdd(&common.gaussian.Re, &common.gaussian.Im, &common.gaussian.norm);
    // Replace G0 by the remainder of G0/G1.
    GetRemainder(&common.gaussian.norm, &common.gaussian.ReG0, &common.gaussian.ImG0, &common.gaussian.ReG1,
      &common.gaussian.ImG1, &common.gaussian.Re, &common.gaussian.Im);
    // Exchange G0 and G1.
    CopyBigInt(&common.gaussian.Tmp, &common.gaussian.ReG0);
    CopyBigInt(&common.gaussian.ReG0, &common.gaussian.ReG1);
    CopyBigInt(&common.gaussian.ReG1, &common.gaussian.Tmp);
    CopyBigInt(&common.gaussian.Tmp, &common.gaussian.ImG0);
    CopyBigInt(&common.gaussian.ImG0, &common.gaussian.ImG1);
    CopyBigInt(&common.gaussian.ImG1, &common.gaussian.Tmp);
    (void)BigIntMultiply(&common.gaussian.Re, &common.gaussian.ReU1, &common.gaussian.Tmp);
    BigIntSubt(&common.gaussian.ReU0, &common.gaussian.Tmp, &common.gaussian.ReU0);
    (void)BigIntMultiply(&common.gaussian.Im, &common.gaussian.ImU1, &common.gaussian.Tmp);
    BigIntAdd(&common.gaussian.ReU0, &common.gaussian.Tmp, &common.gaussian.ReU0);
    (void)BigIntMultiply(&common.gaussian.Im, &common.gaussian.ReU1, &common.gaussian.Tmp);
    BigIntSubt(&common.gaussian.ImU0, &common.gaussian.Tmp, &common.gaussian.ImU0);
    (void)BigIntMultiply(&common.gaussian.Re, &common.gaussian.ImU1, &common.gaussian.Tmp);
    BigIntSubt(&common.gaussian.ImU0, &common.gaussian.Tmp, &common.gaussian.ImU0);
    // Exchange U0 and U1.
    CopyBigInt(&common.gaussian.Tmp, &common.gaussian.ReU0);
    CopyBigInt(&common.gaussian.ReU0, &common.gaussian.ReU1);
    CopyBigInt(&common.gaussian.ReU1, &common.gaussian.Tmp);
    CopyBigInt(&common.gaussian.Tmp, &common.gaussian.ImU0);
    CopyBigInt(&common.gaussian.ImU0, &common.gaussian.ImU1);
    CopyBigInt(&common.gaussian.ImU1, &common.gaussian.Tmp);
  }

  for (;;)
  {
    CopyBigInt(&common.gaussian.Tmp, &common.gaussian.ImG0);
    common.gaussian.Tmp.sign = SIGN_POSITIVE;   // Tmp <- abs(Tmp)
    BigIntSubt(&common.gaussian.ReG0, &common.gaussian.Tmp, &common.gaussian.Tmp);
    if (common.gaussian.Tmp.sign == SIGN_POSITIVE)
    {
      break;
    }
      // Multiply G0 and U0 by i.
    CopyBigInt(&common.gaussian.Tmp, &common.gaussian.ReG0);
    BigIntNegate(&common.gaussian.ImG0, &common.gaussian.ReG0);
    CopyBigInt(&common.gaussian.ImG0, &common.gaussian.Tmp);
    CopyBigInt(&common.gaussian.Tmp, &common.gaussian.ReU0);
    BigIntNegate(&common.gaussian.ImU0, &common.gaussian.ReU0);
    CopyBigInt(&common.gaussian.ImU0, &common.gaussian.Tmp);
  }
  if ((common.gaussian.ReG0.nbrLimbs != 1) || (common.gaussian.ReG0.limbs[0].x != 1) ||
      !BigIntIsZero(&common.gaussian.ImG0))
  {         // G0 is not 1.
    return EXPR_INVALID_PARAM;
  }
  return Modulo(&common.gaussian.ReU0, &common.gaussian.ImU0, RealMod, ImagMod, result);
}

static int Modulo(BigInteger *ReNum, BigInteger *ImNum,
                  const BigInteger *ReDen, const BigInteger *ImDen,
                  BigInteger *result)
{
  BigInteger Re;
  BigInteger Im;
  BigInteger ReMin;
  BigInteger ImMin;
  BigInteger Tmp;
  BigInteger normmin;
  if (BigIntIsZero(ReDen) && BigIntIsZero(ImDen))
  {      // Denominator is zero.
    CopyBigInt(result, ReNum);
    CopyBigInt(result+1, ImNum);
    return 0;
  }
  ReMin.limbs[0].x = 0; 
  ImMin.limbs[0].x = 0;
  ReMin.nbrLimbs = 1; 
  ImMin.nbrLimbs = 1;
  ReMin.sign = SIGN_POSITIVE; 
  ImMin.sign = SIGN_POSITIVE;
  (void)BigIntMultiply(ReDen, ReDen, &common.gaussian.norm);
  (void)BigIntMultiply(ImDen, ImDen, &Tmp);
  BigIntAdd(&common.gaussian.norm, &Tmp, &common.gaussian.norm);
  // Replace Num by the remainder of Num/Den.
  GetRemainder(&common.gaussian.norm, ReNum, ImNum, ReDen, ImDen, &Re, &Im);
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
      (void)BigIntMultiply(&Re, &Re, &common.gaussian.norm);
      (void)BigIntMultiply(&Im, &Im, &Tmp);
      BigIntAdd(&common.gaussian.norm, &Tmp, &common.gaussian.norm);
      BigIntSubt(&common.gaussian.norm, &normmin, &Tmp);
      if ((normmin.sign == SIGN_NEGATIVE) || (Tmp.sign == SIGN_NEGATIVE))
      {
        CopyBigInt(&normmin, &common.gaussian.norm);
        CopyBigInt(&ReMin, &Re);
        CopyBigInt(&ImMin, &Im);
      }
    }
  }
  CopyBigInt(&result[0], &ReMin);
  CopyBigInt(&result[1], &ImMin);
  return 0;
}

// Compute curStackRe + i*curStackIm as the quotient of
// curStackRe + i*curStackIm over curStack2Re + i*curStack2Im
static enum eExprErr GaussianDivide(void)
{
  enum eExprErr retcode = BigIntMultiply(&common.gaussian.curStack2Re, &common.gaussian.curStack2Re, 
    &common.gaussian.curStack3Re);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  retcode = BigIntMultiply(&common.gaussian.curStack2Im, &common.gaussian.curStack2Im,
    &common.gaussian.curStack3Im);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  BigIntAdd(&common.gaussian.curStack3Re, &common.gaussian.curStack3Im, &common.gaussian.norm);       // norm <- re2^2 + im2^2.
  if (BigIntIsZero(&common.gaussian.norm))
  {              // norm is zero.
    return EXPR_INTERM_TOO_HIGH;
  }
  retcode = BigIntMultiply(&common.gaussian.curStackRe, &common.gaussian.curStack2Re, 
    &common.gaussian.curStack3Re);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  retcode = BigIntMultiply(&common.gaussian.curStackIm, &common.gaussian.curStack2Im,
    &common.gaussian.curStack3Im);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  BigIntAdd(&common.gaussian.curStack3Re, &common.gaussian.curStack3Im, &common.gaussian.curTmp);    // Re <- re1*re2 + im1*im2.
  retcode = BigIntMultiply(&common.gaussian.curStackIm, &common.gaussian.curStack2Re,
    &common.gaussian.curStack3Re);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  retcode = BigIntMultiply(&common.gaussian.curStackRe, &common.gaussian.curStack2Im,
    &common.gaussian.curStack3Im);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  BigIntSubt(&common.gaussian.curStack3Re, &common.gaussian.curStack3Im, &common.gaussian.curStack3Im);   // Im <- im1*re2 - re1*im2.
  retcode = BigIntDivide(&common.gaussian.curTmp, &common.gaussian.norm, &common.gaussian.curStackRe);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  retcode = BigIntDivide(&common.gaussian.curStack3Im, &common.gaussian.norm, &common.gaussian.curStackIm);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  return EXPR_OK;
}

// Compute curStackRe + i*curStackIm as the product of
// (curStack2Re + i*curStack2Im) and (curStackRe + i*curStackIm)
static enum eExprErr GaussianMultiply(void)
{
  // Re <- re1*re2 - im1*im2.
  enum eExprErr retcode = BigIntMultiply(&common.gaussian.curStackRe, &common.gaussian.curStack2Re,
    &common.gaussian.curStack3Re);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  retcode = BigIntMultiply(&common.gaussian.curStackIm, &common.gaussian.curStack2Im,
    &common.gaussian.curStack3Im);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  BigIntSubt(&common.gaussian.curStack3Re, &common.gaussian.curStack3Im, &common.gaussian.curTmp);   // Save real part of product in temp var.
  retcode = BigIntMultiply(&common.gaussian.curStackIm, &common.gaussian.curStack2Re, &common.gaussian.curStack3Re); // Im <- im1*re2 + re1*im2.
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  retcode = BigIntMultiply(&common.gaussian.curStackRe, &common.gaussian.curStack2Im, &common.gaussian.curStack3Im);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  BigIntAdd(&common.gaussian.curStack3Re, &common.gaussian.curStack3Im, &common.gaussian.curStackIm);
  CopyBigInt(&common.gaussian.curStackRe, &common.gaussian.curTmp);
  return EXPR_OK;
}

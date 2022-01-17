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

#define PAREN_STACK_SIZE           5000
#define COMPR_STACK_SIZE        1000000

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

static int comprStackValues[COMPR_STACK_SIZE];
static int comprStackOffset[PAREN_STACK_SIZE];
static int stackIndex;
#ifndef lang  
  bool lang;
#endif
char output[3000000];
limb Mult1[MAX_LEN];
limb Mult3[MAX_LEN];
limb Mult4[MAX_LEN];
int q[MAX_LEN];
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
static BigInteger curStackRe;
static BigInteger curStackIm;
static BigInteger curStack2Re;
static BigInteger curStack2Im;
static BigInteger curStack3Re;
static BigInteger curStack3Im;
static BigInteger curStackReBak;
static BigInteger curStackImBak;
static BigInteger curStack2ReBak;
static BigInteger curStack2ImBak;
static BigInteger curTmp;
static BigInteger norm;
static BigInteger Result[2];

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
  const int* ptrStackValue = &comprStackValues[comprStackOffset[2*stackIndex]];
  NumberLength = numLimbs(ptrStackValue);
  IntArray2BigInteger(ptrStackValue, pValueRe);
  ptrStackValue = &comprStackValues[comprStackOffset[(2 * stackIndex) + 1]];
  NumberLength = numLimbs(ptrStackValue);
  IntArray2BigInteger(ptrStackValue, pValueIm);
}

static enum eExprErr setStackValue(const BigInteger* pValueRe, const BigInteger *pValueIm)
{
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
    PARSE_EXPR_GAUSSIAN, NULL);
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

    case TOKEN_RE:
      getCurrentStackValue(&curStackRe, &curStackIm);
      intToBigInteger(&curStackIm, 0);
      break;

    case TOKEN_IM:
      getCurrentStackValue(&curStackIm, &curStackRe);
      intToBigInteger(&curStackIm, 0);
      break;

    case TOKEN_NORM:
      // norm = Re^2 + Im^2
      getCurrentStackValue(&curStackRe, &curStackIm);
      retcode = BigIntMultiply(&curStackRe, &curStackRe, &curStack2Re);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      retcode = BigIntMultiply(&curStackIm, &curStackIm, &curStack2Im);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      BigIntAdd(&curStack2Re, &curStack2Im, &curStackRe);
      intToBigInteger(&curStackIm, 0);
      break;

    case TOKEN_GCD:
      ptrRPNbuffer++;
      nbrParameters = (int)(unsigned char)*ptrRPNbuffer;
      getCurrentStackValue(&curStackRe, &curStackIm);
      for (int parmNbr = 1; parmNbr < nbrParameters; parmNbr++)
      {
        stackIndex--;
        getCurrentStackValue(&curStack2Re, &curStack2Im);
        ComputeGCD();
      }
      break;

    case TOKEN_LCM:
      ptrRPNbuffer++;
      nbrParameters = (int)(unsigned char)*ptrRPNbuffer;
      getCurrentStackValue(&curStackRe, &curStackIm);
      for (int parmNbr = 1; parmNbr < nbrParameters; parmNbr++)
      {
        stackIndex--;
        getCurrentStackValue(&curStack2Re, &curStack2Im);
        retcode = ComputeLCM();
        if (retcode != EXPR_OK)
        {
          return retcode;
        }
      }
      break;

    case TOKEN_MODPOW:
      getCurrentStackValue(&curStack3Re, &curStack3Im);
      stackIndex--;
      getCurrentStackValue(&curStack2Re, &curStack2Im);
      stackIndex--;
      getCurrentStackValue(&curStackRe, &curStackIm);
      retcode = ComputeModPow();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_MODINV:
      getCurrentStackValue(&curStack2Re, &curStack2Im);
      stackIndex--;
      getCurrentStackValue(&curStackRe, &curStackIm);
      retcode = ComputeModInv();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_ISPRIME:
      getCurrentStackValue(&curStackRe, &curStackIm);
      retcode = isPrime();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_F:
      getCurrentStackValue(&curStackRe, &curStackIm);
      if (!BigIntIsZero(&curStackIm))
      {     // Imaginary part of argument must be zero.
        return EXPR_INVALID_PARAM;
      }
      retcode = ComputeFibLucas(0, &curStackRe);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_L:
      getCurrentStackValue(&curStackRe, &curStackIm);
      if (!BigIntIsZero(&curStackIm))
      {     // Imaginary part of argument must be zero.
        return EXPR_INVALID_PARAM;
      }
      retcode = ComputeFibLucas(2, &curStackRe);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_N:
      getCurrentStackValue(&curStackRe, &curStackIm);
      if (!BigIntIsZero(&curStackIm))
      {     // Imaginary part of argument must be zero.
        return EXPR_INVALID_PARAM;
      }
      retcode = ComputeNext(&curStackRe);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_B:
      getCurrentStackValue(&curStackRe, &curStackIm);
      if (!BigIntIsZero(&curStackIm))
      {     // Imaginary part of argument must be zero.
        return EXPR_INVALID_PARAM;
      }
      retcode = ComputeBack(&curStackRe);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_P:
      getCurrentStackValue(&curStackRe, &curStackIm);
      retcode = ComputePartition();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_FACTORIAL:
      getCurrentStackValue(&curStack2Re, &curStack2Im);
      stackIndex--;
      getCurrentStackValue(&curStackRe, &curStackIm);
      if (!BigIntIsZero(&curStackIm) || (curStackRe.sign == SIGN_NEGATIVE))
      {         // Imaginary part must be zero.
        return EXPR_INVALID_PARAM;
      }
      retcode = factorial(&curStackRe, curStackRe.limbs[0].x, curStack2Re.limbs[0].x);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case TOKEN_PRIMORIAL:
      getCurrentStackValue(&curStackRe, &curStackIm);
      if (!BigIntIsZero(&curStackIm) || (curStackRe.sign == SIGN_NEGATIVE))
      {         // Imaginary part must be zero.
        return EXPR_INVALID_PARAM;
      }
      if ((curStackRe.nbrLimbs > 1) || (curStackRe.limbs[0].x > 46049))
      {
        return EXPR_INTERM_TOO_HIGH;
      }
      retcode = primorial(&curStackRe, curStackRe.limbs[0].x);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case OPER_ADD:
      getCurrentStackValue(&curStack2Re, &curStack2Im);
      stackIndex--;
      getCurrentStackValue(&curStackRe, &curStackIm);
      BigIntAdd(&curStackRe, &curStack2Re, &curStackRe);
      BigIntAdd(&curStackIm, &curStack2Im, &curStackIm);
      break;

    case OPER_SUBT:
      getCurrentStackValue(&curStack2Re, &curStack2Im);
      stackIndex--;
      getCurrentStackValue(&curStackRe, &curStackIm);
      BigIntSubt(&curStackRe, &curStack2Re, &curStackRe);
      BigIntSubt(&curStackIm, &curStack2Im, &curStackIm);
      break;

    case OPER_UNARY_MINUS:
      getCurrentStackValue(&curStackRe, &curStackIm);
      BigIntChSign(&curStackRe);
      BigIntChSign(&curStackIm);
      break;

    case OPER_MULTIPLY:
      getCurrentStackValue(&curStack2Re, &curStack2Im);
      stackIndex--;
      getCurrentStackValue(&curStackRe, &curStackIm);
      retcode = GaussianMultiply();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case OPER_DIVIDE:
      getCurrentStackValue(&curStack2Re, &curStack2Im);
      stackIndex--;
      getCurrentStackValue(&curStackRe, &curStackIm);
      retcode = GaussianDivide();
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      break;

    case OPER_REMAINDER:
      getCurrentStackValue(&curStack2Re, &curStack2Im);
      stackIndex--;
      getCurrentStackValue(&curStackRe, &curStackIm);
      retcode = Modulo(&curStackRe, &curStackIm, &curStack2Re, &curStack2Im, Result);
      if (retcode != EXPR_OK)
      {
        return retcode;
      }
      CopyBigInt(&curStackRe, &Result[0]);
      CopyBigInt(&curStackIm, &Result[1]);
      break;

    case TOKEN_END_EXPON:
      getCurrentStackValue(&curStack2Re, &curStack2Im);
      stackIndex--;
      getCurrentStackValue(&curStackRe, &curStackIm);
      retcode = ComputePower(&curStackRe, &curStack2Re, &curStackIm, &curStack2Im);
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
        retcode = setStackValue(&curStackRe, &curStackIm);
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
  return EXPR_OK;
}

static enum eExprErr isPrime(void)
{
  if (BigIntIsZero(&curStackRe))
  {
    if ((curStackIm.limbs[0].x % 4) != 3)
    {
      intToBigInteger(&curStackRe, 0);    // Number not prime.
      intToBigInteger(&curStackIm, 0);
      return EXPR_OK;
    }
    CopyBigInt(&curStack2Re, &curStackIm);
  }
  else if (BigIntIsZero(&curStackIm))
  {
    if ((curStackRe.limbs[0].x % 4) != 3)
    {
      intToBigInteger(&curStackRe, 0);    // Number not prime.
      intToBigInteger(&curStackIm, 0);
      return EXPR_OK;
    }
    CopyBigInt(&curStack2Re, &curStackRe);
  }
  else
  {
    enum eExprErr retcode;
    retcode = BigIntMultiply(&curStackRe, &curStackRe, &curStack2Re);
    if (retcode != EXPR_OK)
    {
      return retcode;
    }
    retcode = BigIntMultiply(&curStackIm, &curStackIm, &curStack2Im);
    if (retcode != EXPR_OK)
    {
      return retcode;
    }
    BigIntAdd(&curStack2Re, &curStack2Im, &curStack2Re);
  }
#ifdef FACTORIZATION_APP
  if (BpswPrimalityTest(&curStack2Re, NULL) == 0)
#else
  if (BpswPrimalityTest(&curStack2Re) == 0)
#endif
  {    // Argument is a probable prime.
    intToBigInteger(&curStackRe, -1);
  }
  else
  {    // Argument is not a probable prime.
    intToBigInteger(&curStackRe, 0);
  }
  intToBigInteger(&curStackIm, 0);
  return EXPR_OK;
}

static enum eExprErr ComputePartition(void)
{
  int val;

  if (!BigIntIsZero(&curStackIm))
  {         // Imaginary part must be zero.
    return EXPR_INVALID_PARAM;
  }
  if (curStackRe.sign == SIGN_NEGATIVE)
  {
    return EXPR_INVALID_PARAM;
  }
  if (curStackRe.nbrLimbs > 1)
  {
    return EXPR_INTERM_TOO_HIGH;
  }
  val = curStackRe.limbs[0].x;
  if (val > 100000)
  {
    return EXPR_INTERM_TOO_HIGH;
  }
  partition(val, &curStackRe);
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
  // Re <- ((Re1*Re2+Im1*Im2)*2/norm+1)/2
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
  // Im <- ((Im1*Re2-Re1*Im2)*2/norm+1)/2
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

// Compute curStack = gcd(curStack, curStack2)
static void ComputeGCD(void)
{
  BigInteger Tmp;
  while (!BigIntIsZero(&curStack2Re) || !BigIntIsZero(&curStack2Im))
  {   // Second argument is not zero.
    (void)BigIntMultiply(&curStack2Re, &curStack2Re, &curStack3Re);
    (void)BigIntMultiply(&curStack2Im, &curStack2Im, &curStack3Im);
    BigIntAdd(&curStack3Re, &curStack3Im, &Tmp);
    // Get remainder of (Re1+i*Im1)/(Re2*i*Im2)
    // Overwrite Re1 and Im1 with that remainder.
    // Re and Im are used as temporary values.
    GetRemainder(&Tmp, &curStackRe, &curStackIm, &curStack2Re, &curStack2Im,
      &curStack3Re, &curStack3Im);
    // Exchange Re1 and Re2.
    CopyBigInt(&Tmp, &curStackRe);
    CopyBigInt(&curStackRe, &curStack2Re);
    CopyBigInt(&curStack2Re, &Tmp);
    // Exchange Im1 and Im2.
    CopyBigInt(&Tmp, &curStackIm);
    CopyBigInt(&curStackIm, &curStack2Im);
    CopyBigInt(&curStack2Im, &Tmp);
    curStack2Re.sign = SIGN_POSITIVE;  // Re2 <- abs(Re2)
    curStack2Im.sign = SIGN_POSITIVE;  // Im2 <- abs(Im2)
  }
  CopyBigInt(&Tmp, &curStackIm);
  Tmp.sign = SIGN_POSITIVE;  // Tmp <- abs(Tmp)
  BigIntSubt(&curStackRe, &Tmp, &Tmp);
  while (Tmp.sign == SIGN_NEGATIVE)
  {    // Multiply by i.
    CopyBigInt(&curStack3Re, &curStackRe);
    BigIntNegate(&curStackIm, &curStackRe);
    CopyBigInt(&curStackIm, &curStack3Re);
  }
}

static enum eExprErr ComputeLCM(void)
{
  enum eExprErr retcode;
  if ((BigIntIsZero(&curStackRe) && BigIntIsZero(&curStackIm)) ||
    (BigIntIsZero(&curStack2Re) && BigIntIsZero(&curStack2Im)))
  {    // If any of the arguments is zero, the LCM is zero.
    intToBigInteger(&curStackRe, 0);
    intToBigInteger(&curStackIm, 0);
    return EXPR_OK;
  }
  // Save both parameters because ComputeGCD will overwrite them.
  CopyBigInt(&curStackReBak, &curStackRe);
  CopyBigInt(&curStackImBak, &curStackIm);
  CopyBigInt(&curStack2ReBak, &curStack2Re);
  CopyBigInt(&curStack2ImBak, &curStack2Im);
  ComputeGCD();   // curStackRe + i*curStackIm = GCD.
  // Divide first parameter by GCD.
  CopyBigInt(&curStack2Re, &curStackRe);    // Divisor
  CopyBigInt(&curStack2Im, &curStackIm);
  CopyBigInt(&curStackRe, &curStackReBak);  // Dividend
  CopyBigInt(&curStackIm, &curStackImBak);
  retcode = GaussianDivide();
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  // Multiply the quotient by the second parameter.
  CopyBigInt(&curStack2Re, &curStack2ReBak);
  CopyBigInt(&curStack2Im, &curStack2ImBak);
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
  BigIntAdd(&ReTmp, &ImTmp, &norm);  // norm <- re1^2 + im1^2.
  logNorm = logLimbs(norm.limbs, norm.nbrLimbs);  // Get logarithm of norm.
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
        BigIntSubt(&ReTmp, &ImTmp, &curTmp);
        (void)BigIntMultiply(Re1, &Im, &ReTmp);       // Im <- re1*im + im1*re.
        (void)BigIntMultiply(Im1, &Re, &ImTmp);
        BigIntAdd(&ReTmp, &ImTmp, &Im);
        CopyBigInt(&Re, &curTmp);
      }
    }
  }
  CopyBigInt(Re1, &Re);
  CopyBigInt(Im1, &Im);
  return EXPR_OK;
}

static enum eExprErr ComputeModPow(void)
{
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
  CopyBigInt(&ReBase, &curStackRe);
  CopyBigInt(&ImBase, &curStackIm);
  CopyBigInt(&ReExp, &curStack2Re);
  CopyBigInt(&ImExp, &curStack2Im);
  CopyBigInt(&ReMod, &curStack3Re);
  CopyBigInt(&ImMod, &curStack3Im);

  if (!BigIntIsZero(&ImExp))
  {         // Imaginary part must be zero.
    return EXPR_INVALID_PARAM;
  }
  if (ReExp.sign == SIGN_NEGATIVE)
  {
    enum eExprErr retcode;
    BigIntNegate(&ReExp, &ReExp);
    retcode = ModInv(&ReBase, &ImBase, &ReMod, &ImMod, Result);
    if (retcode != EXPR_OK)
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
    CopyBigInt(&curStackRe, &ReBase);
    CopyBigInt(&curStackIm, &ImBase);
  }
  else
  {                            /* Modulus is not zero */
    (void)BigIntMultiply(&ReMod, &ReMod, &ReTmp);
    (void)BigIntMultiply(&ImMod, &ImMod, &ImTmp);
    BigIntAdd(&ReTmp, &ImTmp, &norm);
    for (int index = ReExp.nbrLimbs - 1; index >= 0; index--)
    {
      int groupExp = ReExp.limbs[index].x;
      for (unsigned int mask = HALF_INT_RANGE_U; mask > 0U; mask >>= 1)
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
        if (((unsigned int)groupExp & mask) != 0U)
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
    CopyBigInt(&curStackRe, &Result[0]);
    CopyBigInt(&curStackIm, &Result[1]);
    return EXPR_OK;
  }
  return EXPR_OK;
}

static enum eExprErr ComputeModInv(void)
{
  enum eExprErr retcode = ModInv(&curStackRe, &curStackIm,
    &curStack2Re, &curStack2Im, Result);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  CopyBigInt(&curStackRe, &Result[0]);
  CopyBigInt(&curStackIm, &Result[1]);
  return EXPR_OK;
}

static enum eExprErr ModInv(const BigInteger *RealNbr, const BigInteger *ImagNbr,
                  const BigInteger *RealMod, const BigInteger *ImagMod,
                  BigInteger *result)
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
  return Modulo(&ReU0, &ImU0, RealMod, ImagMod, result);
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
  CopyBigInt(&result[0], &ReMin);
  CopyBigInt(&result[1], &ImMin);
  return 0;
}

// (curStackRe + i*curStackIm) <- (curStackRe + i*curStackIm) / (curStack2Re + i*curStack2Im)
static enum eExprErr GaussianDivide(void)
{
  enum eExprErr retcode = BigIntMultiply(&curStack2Re, &curStack2Re, &curStack3Re);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  retcode = BigIntMultiply(&curStack2Im, &curStack2Im, &curStack3Im);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  BigIntAdd(&curStack3Re, &curStack3Im, &norm);       // norm <- re2^2 + im2^2.
  if (BigIntIsZero(&norm))
  {              // norm is zero.
    return EXPR_INTERM_TOO_HIGH;
  }
  retcode = BigIntMultiply(&curStackRe, &curStack2Re, &curStack3Re);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  retcode = BigIntMultiply(&curStackIm, &curStack2Im, &curStack3Im);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  BigIntAdd(&curStack3Re, &curStack3Im, &curTmp);    // Re <- re1*re2 + im1*im2.
  retcode = BigIntMultiply(&curStackIm, &curStack2Re, &curStack3Re);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  retcode = BigIntMultiply(&curStackRe, &curStack2Im, &curStack3Im);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  BigIntSubt(&curStack3Re, &curStack3Im, &curStack3Im);   // Im <- im1*re2 - re1*im2.
  retcode = BigIntDivide(&curTmp, &norm, &curStackRe);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  retcode = BigIntDivide(&curStack3Im, &norm, &curStackIm);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  return EXPR_OK;
}

// (curStackRe + i*curStackIm) <- (curStack2Re + i*curStack2Im) * (curStackRe + i*curStackIm)
static enum eExprErr GaussianMultiply(void)
{
  // Re <- re1*re2 - im1*im2.
  enum eExprErr retcode = BigIntMultiply(&curStackRe, &curStack2Re, &curStack3Re);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  retcode = BigIntMultiply(&curStackIm, &curStack2Im, &curStack3Im);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  BigIntSubt(&curStack3Re, &curStack3Im, &curTmp);   // Save real part of product in temp var.
  retcode = BigIntMultiply(&curStackIm, &curStack2Re, &curStack3Re); // Im <- im1*re2 + re1*im2.
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  retcode = BigIntMultiply(&curStackRe, &curStack2Im, &curStack3Im);
  if (retcode != EXPR_OK)
  {
    return retcode;
  }
  BigIntAdd(&curStack3Re, &curStack3Im, &curStackIm);
  CopyBigInt(&curStackRe, &curTmp);
  return EXPR_OK;
}

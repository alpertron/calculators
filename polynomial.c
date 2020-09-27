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

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "bignbr.h"
#include "highlevel.h"
#include "polynomial.h"
#include "showtime.h"
#define STACK_OPER_SIZE      100
#define TOKEN_NUMBER         '0'
#define TOKEN_START_EXPON    '1'
#define TOKEN_END_EXPON      '2'
#define TOKEN_UNARY_MINUS    '3'
#define TOKEN_GCD            '4'
#define TOKEN_DER            '5'
#define KARATSUBA_POLY_CUTOFF 16
#define COMPRESSED_POLY_MAX_LENGTH 1000000
extern char* ptrOutput2;
static void showPolynomial(char** pptrOutput, int* ptrPoly, int polyDegree, int groupLength);
BigInteger primeMod;              // p
int exponentMod;                  // k
BigInteger powerMod;              // p^k
int modulusIsZero;
int degree;
static BigInteger value;
BigInteger operand1, operand2, operand3, operand4, operand5;
extern limb TestNbr[MAX_LEN];
extern limb MontgomeryMultR2[MAX_LEN];
extern limb MontgomeryMultR1[MAX_LEN];
extern int NumberLength, NumberLengthR1;
static int prime;
static char RPNbuffer[COMPRESSED_POLY_MAX_LENGTH];
int values[COMPRESSED_POLY_MAX_LENGTH];
static int polyInv[COMPRESSED_POLY_MAX_LENGTH];
static int valuesIndex;
int valuesPrime[COMPRESSED_POLY_MAX_LENGTH];
int poly1[COMPRESSED_POLY_MAX_LENGTH];
int poly2[COMPRESSED_POLY_MAX_LENGTH];
int poly3[COMPRESSED_POLY_MAX_LENGTH];
int poly4[COMPRESSED_POLY_MAX_LENGTH];
int poly5[COMPRESSED_POLY_MAX_LENGTH];
int polyS[COMPRESSED_POLY_MAX_LENGTH];
static int polyT[COMPRESSED_POLY_MAX_LENGTH];
static int polyMultM[COMPRESSED_POLY_MAX_LENGTH];
static int polyMultT[COMPRESSED_POLY_MAX_LENGTH];
int polyMultTemp[COMPRESSED_POLY_MAX_LENGTH];
int polyLifted[COMPRESSED_POLY_MAX_LENGTH];
int polyLiftedBak[COMPRESSED_POLY_MAX_LENGTH];
unsigned char pretty, onlyEvaluate = 0;
struct sFactorInfo factorInfo[MAX_DEGREE];
static BigInteger coeff[2 * KARATSUBA_POLY_CUTOFF];
int nbrFactorsFound;
int *ptrOrigPoly;
int degreeOrigPoly;
static int GcdPolynomialExpr(int *ptrArgument1, int *ptrArgument2);

enum eAdjustType
{
  ADJUST_PERFORM_ADDITION = 0,
  ADJUST_PERFORM_SUBTRACTION
};

int numLimbs(int *pLen)
{
  int nbrLimbs = *pLen;
  if (nbrLimbs < 0)
  {
    nbrLimbs = -nbrLimbs;
  }
  return nbrLimbs;
}

static int isFunc(char **ppcInput, char *funcName)
{
  char *pcInput = *ppcInput - 1;
  char *ptrFuncName = funcName;
  while (*ptrFuncName)
  {
    if ((*pcInput & 0xDF) != *ptrFuncName)
    {
      return 0;    // Function name not found.
    }
    pcInput++;
    ptrFuncName++;
  }
  *ppcInput = pcInput;
  return 1;        // Function name was found.
}

// Convert to Reverse Polish Notation using Shunting-yard algorithm
// If the number is part of the exponent it has to be reduced by
// using the modulus phi = (p-1)*p^(k-1), otherwise the modulus will
// be power = p^k.
static int ConvertToReversePolishNotation(char *input, char *ptrOutput)
{
  int exponOperatorCounter = 0;
  char stackOper[STACK_OPER_SIZE];
  int stackOperIndex = 0;
  int prevTokenIsNumber = FALSE;
  char s;
  char *ptrInput;
  int index, limb, bitNbr;
  while (*input != '\0')
  {
    char c = *input++;
    if (c == ' ' || c == 9)
    {          // Ignore any spaces and tabs.
      continue;
    }
    if (c == '*' && *input == '*')
    {          // Convert double asterisk to exponentiation.
      c = '^';
      input++;
    }
    if (prevTokenIsNumber != FALSE)
    {
      prevTokenIsNumber = FALSE;
      if (c == '+' || c == '-')
      {        // Binary plus or minus
        while (stackOperIndex > 0)
        {      // Send operators to output.
          s = stackOper[--stackOperIndex];
          if (s != '+' && s != '-' && s != '*' && s != '^' && s != '/' && s != '%' && s != TOKEN_UNARY_MINUS)
          {    // Operator on stack has less precedence.
            stackOperIndex++;
            break;
          }
          if (s == '^')
          {
            if (--exponOperatorCounter == 0)
            {
              *ptrOutput++ = TOKEN_END_EXPON;
            }
          }
          *ptrOutput++ = s;
        }
        stackOper[stackOperIndex++] = c;  // Push operator onto stack.
      }
      else if (c == '*' || c == '(' || c == 'x' || c == 'X' || c == '/' || c == '%')
      {
        while (stackOperIndex > 0)
        {      // Send operators to output.
          s = stackOper[--stackOperIndex];
          if (s != '^' && s != '*' && s != '/' && s != '%')
          {    // Operator on stack has less precedence.
            stackOperIndex++;
            break;
          }
          if (s == '^')
          {
            if (--exponOperatorCounter == 0)
            {
              *ptrOutput++ = TOKEN_END_EXPON;
            }
          }
          *ptrOutput++ = s;
        }
        if (c == 'x' || c == 'X')
        {
          if (exponOperatorCounter != 0)
          {
            return EXPR_CANNOT_USE_X_IN_EXPONENT;
          }
          stackOper[stackOperIndex++] = '*';    // Push operator onto stack.
          *ptrOutput++ = 'x';
          prevTokenIsNumber = TRUE;
        }
        else
        {
          if (c == '(')
          {
            stackOper[stackOperIndex++] = '*';  // Push operator onto stack.
          }
          stackOper[stackOperIndex++] = c;      // Push operator onto stack.
        }
      }
      else if (c == '^')
      {
        while (stackOperIndex > 0)
        {      // Send operators to output.
          s = stackOper[--stackOperIndex];
          if (s != '^')
          {    // Operator on stack has less precedence.
            stackOperIndex++;
            break;
          }
          if (--exponOperatorCounter == 0)
          {
            *ptrOutput++ = TOKEN_END_EXPON;
          }
          *ptrOutput++ = s;
        }
        stackOper[stackOperIndex++] = c;  // Push operator onto stack.
        if (++exponOperatorCounter == 1)
        {
          *ptrOutput++ = TOKEN_START_EXPON;
        }
      }
      else if (c == ')' || c == ',')
      {
        s = '\0'; // Assume parenthesis mismatch.
        while (stackOperIndex > 0)
        {      // Send operators to output.
          s = stackOper[--stackOperIndex];
          if (s == '(')
          {    // Operator on stack has less precedence.
            break;
          }
          if (s == '^')
          {
            if (--exponOperatorCounter == 0)
            {
              *ptrOutput++ = TOKEN_END_EXPON;
            }
          }
          *ptrOutput++ = s;
        }
        if (s != '(')
        {     // Parentheses mismatch.
          return EXPR_PAREN_MISMATCH;
        }
        if (c == ',')
        {
          stackOper[stackOperIndex++] = s;  // Push back paren.
          prevTokenIsNumber = FALSE;
        }
        else
        {
          prevTokenIsNumber = TRUE;
        }
      }
      else
      {
        return EXPR_SYNTAX_ERROR;
      }
    }
    else
    {    // Not a number before this token.
      if (c == '+')
      {          // Unary plus, nothing to do.

      }
      else if (c == '-')
      {          // Unary minus.
        stackOper[stackOperIndex++] = TOKEN_UNARY_MINUS;  // Push operator onto stack.
      }
      else if (c == '(')
      {          // Open parenthesis.
        stackOper[stackOperIndex++] = '(';  // Push operator onto stack.
      }
      else if (c == 'x' || c == 'X')
      {
        if (exponOperatorCounter != 0)
        {
          return EXPR_CANNOT_USE_X_IN_EXPONENT;
        }
        *ptrOutput++ = 'x';
        prevTokenIsNumber = TRUE;
      }
      else if (c >= '0' && c <= '9')
      {          // Number.
        ptrInput = input;
        while (*ptrInput >= '0' && *ptrInput <= '9')
        {        // Find end of number.
          ptrInput++;
        }
        Dec2Bin(input - 1, value.limbs, (int)(ptrInput + 1 - input), &value.nbrLimbs);
        input = ptrInput;
        value.sign = SIGN_POSITIVE;
        if (exponOperatorCounter != 0)
        {
          if (value.nbrLimbs != 1)
          {
            return EXPR_EXPONENT_TOO_LARGE;
          }
          *ptrOutput++ = TOKEN_NUMBER;
          limb = (int)value.limbs[0].x;
          for (bitNbr = BITS_PER_GROUP; bitNbr > 0; bitNbr -= 8)
          {
            *ptrOutput++ = (char)limb;
            limb >>= 8;
          }
        }
        else
        {
          if (!modulusIsZero)
          {
            // Reduce mod power outside exponent.
            (void)BigIntRemainder(&value, &powerMod, &value);
            if (value.nbrLimbs < NumberLength)
            {    // Fill with zeros.
              memset(&value.limbs[value.nbrLimbs], 0, (NumberLength - value.nbrLimbs) * sizeof(int));
            }
            // Convert to Montgomery notation.
            modmult(value.limbs, MontgomeryMultR2, value.limbs);
            value.nbrLimbs = NumberLength;
          }
          *ptrOutput++ = TOKEN_NUMBER;
          while (value.nbrLimbs > 1)
          {
            if (value.limbs[value.nbrLimbs-1].x != 0)
            {
              break;
            }
            value.nbrLimbs--;
          }
          *ptrOutput++ = (char)(value.nbrLimbs >> 8);
          *ptrOutput++ = (char)value.nbrLimbs;
          for (index = 0; index < value.nbrLimbs; index++)
          {
            limb = (int)value.limbs[index].x;
            for (bitNbr = BITS_PER_GROUP; bitNbr > 0; bitNbr -= 8)
            {
              *ptrOutput++ = (char)limb;
              limb >>= 8;
            }
          }
        }
        prevTokenIsNumber = TRUE;
      }
      else if (isFunc(&input, "GCD"))
      {
        stackOper[stackOperIndex++] = TOKEN_GCD;  // Push token onto stack.
      }
      else if (isFunc(&input, "DER"))
      {
        stackOper[stackOperIndex++] = TOKEN_DER;  // Push token onto stack.
      }
      else
      {
        return EXPR_SYNTAX_ERROR;
      }
    }
  }
  if (prevTokenIsNumber == FALSE)
  {
    return EXPR_SYNTAX_ERROR;
  }
  while (stackOperIndex > 0)
  {      // Send operators to output.
    s = stackOper[--stackOperIndex];
    if (s == '^')
    {
      if (--exponOperatorCounter == 0)
      {
        *ptrOutput++ = TOKEN_END_EXPON;
      }
    }
    if (s == '(')
    {    // Operator on stack has less precedence.
      return EXPR_PAREN_MISMATCH;
    }
    *ptrOutput++ = s;
  }
  *ptrOutput = '\0';
  return EXPR_OK;
}

static int NegatePolynomialExpr(int *ptrArgument)
{
  int *ptrValue1 = ptrArgument;
  int val = *ptrValue1;
  if (val <= 0)
  {          // Monomial
    if (*(ptrValue1 + 1) != 1 || *(ptrValue1 + 2) != 0)
    {        // Coefficient is not zero
      if (modulusIsZero)
      {
        *(ptrValue1 + 1) = -*(ptrValue1 + 1);    // Negate number.
      }
      else
      {
        IntArray2BigInteger(ptrValue1 + 1, &operand1);
        BigIntSubt(&powerMod, &operand1, &operand1);
        BigInteger2IntArray(ptrValue1 + 1, &operand1);
      }
    }
    return EXPR_OK;
  }
  else
  {          // Polynomial. Use poly1 as temporary polynomial.
    int currentDegree;
    int *ptrValue2 = poly1;
    ptrValue1++;
    for (currentDegree = val; currentDegree >= 0; currentDegree--)
    {
      if (*ptrValue1 == 1 && *(ptrValue1 + 1) == 0)
      {          // If value is zero, it does not have to be changed.
        *ptrValue2++ = 1;
        *ptrValue2++ = 0;
        ptrValue1 += 2;              // Point to next coefficient.
      }
      else if (modulusIsZero)
      {                              // Integer polynomial.
        int nbrLimbs = *ptrValue1;
        *ptrValue1 = -nbrLimbs;      // Negate number.
        if (nbrLimbs < 0)
        {
          nbrLimbs = -nbrLimbs;
        }
        ptrValue1 += nbrLimbs + 1;   // Point to next coefficient.
      }
      else
      {                              // Polynomial modulo prime.
        IntArray2BigInteger(ptrValue1, &operand1);
        BigIntSubt(&powerMod, &operand1, &operand1);
        BigInteger2IntArray(ptrValue2, &operand1);
        ptrValue1 += *ptrValue1 + 1; // Point to next coefficient.
        ptrValue2 += *ptrValue2 + 1;
      }
    }
    if (!modulusIsZero)
    {
      if (ptrArgument + 1 + (ptrValue2 - &poly1[0]) > values + sizeof(values) / sizeof(values[0]))
      {
        return EXPR_OUT_OF_MEMORY;
      }
      memcpy(ptrArgument + 1, poly1, (ptrValue2 - &poly1[0]) * sizeof(int));
    }
  }
  return EXPR_OK;
}

void UncompressBigIntegerB(int *ptrValues, BigInteger *bigint)
{
  if (modulusIsZero)
  {
    NumberLength = numLimbs(ptrValues);
  }
  IntArray2BigInteger(ptrValues, bigint);
}

static int AddPolynomialExpr(int *ptrArgument1, int *ptrArgument2)
{
  int *ptrValue1, currentDegree;
  int degreeMin, degreeMax;
  int degreePoly=0, degreeMono;
  int degree1 = *ptrArgument1;
  int degree2 = *ptrArgument2;
  if (degree1 <= 0)
  {
    if (degree1 == degree2)
    {      // Sum of two monomials of same degree.
      UncompressBigIntegerB(ptrArgument1 + 1, &operand1);
      UncompressBigIntegerB(ptrArgument2 + 1, &operand2);
      if (modulusIsZero)
      {
        BigIntAdd(&operand1, &operand2, &operand1);
        NumberLength = operand1.nbrLimbs;
      }
      else
      {
        AddBigNbrModN(operand1.limbs, operand2.limbs, operand1.limbs, powerMod.limbs, powerMod.nbrLimbs);
      }
      BigInteger2IntArray(ptrArgument1 + 1, &operand1);
      if (BigIntIsZero(&operand1))
      {                     // Sum is zero: set degree to zero.
        *ptrArgument1 = 0;
      }
      valuesIndex = (int)(ptrArgument1 + 2 + *(ptrArgument1 + 1) - &values[0]);
      return EXPR_OK;
    }
    if (degree2 <= 0)
    {     // Sum of two monomials of different degree. Generate polynomial.
      degree1 = -degree1;   // Convert to positive.
      degree2 = -degree2;
      if (degree1 > degree2)
      {
        degreeMax = degree1;
        degreeMin = degree2;
        UncompressBigIntegerB(ptrArgument1+1, &operand2);
        UncompressBigIntegerB(ptrArgument2+1, &operand1);
      }
      else
      {
        degreeMax = degree2;
        degreeMin = degree1;
        UncompressBigIntegerB(ptrArgument1+1, &operand1);
        UncompressBigIntegerB(ptrArgument2+1, &operand2);
      }
      ptrValue1 = ptrArgument1;
      *ptrValue1++ = degreeMax;
      for (currentDegree = 0; currentDegree < degreeMin; currentDegree++)
      {
        *ptrValue1++ = 1;  // Number of limbs
        *ptrValue1++ = 0;  // Value = zero.
      }
      BigInteger2IntArray(ptrValue1, &operand1);
      ptrValue1 += 1 + numLimbs(ptrValue1);
      for (currentDegree++; currentDegree < degreeMax; currentDegree++)
      {
        *ptrValue1++ = 1;  // Number of limbs
        *ptrValue1++ = 0;  // Value = zero.
      }
      NumberLength = operand2.nbrLimbs;
      BigInteger2IntArray(ptrValue1, &operand2);
      valuesIndex = (int)(ptrValue1 + 1 + numLimbs(ptrValue1) - &values[0]);
      return EXPR_OK;
    }
  }
  if (degree1 > 0 && degree2 > 0)
  {           // Sum of two polynomials.
    int *ptrPolyMin, *ptrPolyMax;
    if (degree1 > degree2)
    {
      degreeMin = degree2;
      degreeMax = degree1;
      ptrPolyMin = ptrArgument2+1;
      ptrPolyMax = ptrArgument1+1;
    }
    else
    {
      degreeMin = degree1;
      degreeMax = degree2;
      ptrPolyMin = ptrArgument1+1;
      ptrPolyMax = ptrArgument2+1;
    }
    ptrValue1 = poly1;
    *ptrValue1++ = degreeMax;
    for (currentDegree = 0; currentDegree <= degreeMin; currentDegree++)
    {
      UncompressBigIntegerB(ptrPolyMin, &operand1);
      UncompressBigIntegerB(ptrPolyMax, &operand2);
      if (modulusIsZero)
      {
        BigIntAdd(&operand1, &operand2, &operand1);
        NumberLength = operand1.nbrLimbs;
      }
      else
      {
        AddBigNbrModN(operand1.limbs, operand2.limbs, operand1.limbs, powerMod.limbs, powerMod.nbrLimbs);
      }
      BigInteger2IntArray(ptrValue1, &operand1);
      ptrPolyMin += 1 + numLimbs(ptrPolyMin);
      ptrPolyMax += 1 + numLimbs(ptrPolyMax);
      ptrValue1 += 1 + numLimbs(ptrValue1);
    }
    for (; currentDegree <= degreeMax; currentDegree++)
    {
      memcpy(ptrValue1, ptrPolyMax, (1 + numLimbs(ptrPolyMax))*sizeof(int));
      ptrPolyMax += 1 + numLimbs(ptrPolyMax);
      ptrValue1 += 1 + numLimbs(ptrValue1);
    }
    memcpy(ptrArgument1, poly1, (ptrValue1 - &poly1[0])*sizeof(int));
    valuesIndex = (int)(ptrArgument1 - &values[0] + (ptrValue1 - &poly1[0]));
  }
  else
  {                 // Sum of polynomial and monomial.
    if (degree1 < 0)
    {
      UncompressBigIntegerB(ptrArgument1 + 1, &operand1);
      memmove(ptrArgument1 + 1, ptrArgument1 + 3 + numLimbs(ptrArgument1 + 1), (&values[valuesIndex] - ptrArgument2) * sizeof(int));
      degreeMono = -degree1;
      degreePoly = degree2;
      ptrValue1 = ptrArgument2 + 1;
    }
    else
    {
      UncompressBigIntegerB(ptrArgument2 + 1, &operand1);
      degreeMono = -degree2;
      degreePoly = degree1;
      ptrValue1 = ptrArgument1 + 1;
    }
    if (degreeMono > degreePoly)
    {
      int* ptrValue2 = ptrArgument1 + 1;
      *ptrArgument1 = degreeMono;
      for (currentDegree = 0; currentDegree <= degreePoly; currentDegree++)
      {
        UncompressBigIntegerB(ptrValue1, &operand2);
        BigInteger2IntArray(ptrValue2, &operand2);
        ptrValue1 += 1 + numLimbs(ptrValue1);
        ptrValue2 += 1 + numLimbs(ptrValue2);
      }
      for (; currentDegree < degreeMono; currentDegree++)
      {
        *ptrValue1++ = 1;  // Number of limbs
        *ptrValue1++ = 0;  // Value = zero.
      }
      BigInteger2IntArray(ptrValue1, &operand1);
      valuesIndex = (int)(ptrValue1 + 1 + numLimbs(ptrValue1) - &values[0]);
    }
    else
    {        // Degree of polynomial greater than degree of monomial.
      int differenceOfDegrees, nbrLimbsSum;

      *ptrArgument1 = degreePoly;
      ptrValue1 = ptrArgument1 + 1;
      for (currentDegree = 0; currentDegree < degreeMono; currentDegree++)
      {
        ptrValue1 += 1 + numLimbs(ptrValue1);
      }
      UncompressBigIntegerB(ptrValue1, &operand2);
      if (modulusIsZero)
      {
        BigIntAdd(&operand1, &operand2, &operand1);
        NumberLength = operand1.nbrLimbs;
      }
      else
      {
        AddBigNbrModN(operand1.limbs, operand2.limbs, operand1.limbs, powerMod.limbs, powerMod.nbrLimbs);
      }
      BigInteger2IntArray(poly1, &operand1);
      nbrLimbsSum = numLimbs(&poly1[0]);
      differenceOfDegrees = nbrLimbsSum - numLimbs(ptrValue1);
      if (differenceOfDegrees > 0)
      {    // New coefficient is greater than old coefficient.
        memmove(ptrValue1 + differenceOfDegrees, ptrValue1, (ptrArgument2 - ptrValue1) * sizeof(int));
      }
      else if (differenceOfDegrees < 0)
      {
        memmove(ptrValue1, ptrValue1 - differenceOfDegrees, (ptrArgument2 - ptrValue1 - differenceOfDegrees) * sizeof(int));
      }
      memcpy(ptrValue1, poly1, (1 + nbrLimbsSum) * sizeof(int));
      valuesIndex += differenceOfDegrees;
    }
  }
           // Reduce degree if leading coefficient is zero.
  degreeMax = 0;
  ptrValue1 = ptrArgument1 + 1;
  for (currentDegree = 0; currentDegree <= degreePoly; currentDegree++)
  {
    if (*ptrValue1 != 1 || *(ptrValue1 + 1) != 0)
    {                    // Coefficient is not zero
      degreeMax = currentDegree;
    }
    ptrValue1 += 1 + numLimbs(ptrValue1);
  }
  *ptrArgument1 = degreeMax;
  return EXPR_OK;
}

// Multiply two groups of nbrLen coefficients. The first one starts at
// idxFactor1 and the second one at idxFactor2. The 2*nbrLen coefficient
// result is stored starting at idxFactor1. Use arrAux as temporary storage.
// Accumulate products by result coefficient.
static void ClassicalPolyMult(int idxFactor1, int idxFactor2, int coeffLen, int nbrLimbs)
{
  int i, j;
  int *ptrFactor1, *ptrFactor2;
#ifndef _USING64BITS_
  limb result;
#endif
  for (i = 0; i < 2 * coeffLen - 1; i++)
  {    // Process each limb of product (least to most significant limb).
    if (i < coeffLen)
    {   // Processing first half (least significant) of product.
      ptrFactor2 = &polyMultTemp[(idxFactor2 + i)*nbrLimbs];
      ptrFactor1 = &polyMultTemp[idxFactor1*nbrLimbs];
      j = i;
    }
    else
    {  // Processing second half (most significant) of product.
      ptrFactor2 = &polyMultTemp[(idxFactor2 + coeffLen - 1)*nbrLimbs];
      ptrFactor1 = &polyMultTemp[(idxFactor1 + i - coeffLen + 1)*nbrLimbs];
      j = 2 * (coeffLen-1) - i;
    }
    memset(coeff[i].limbs, 0, nbrLimbs * sizeof(limb));
    if (nbrLimbs == 2)
    {    // Optimization for the case when there is only one limb.
      int modulus = TestNbr[0].x;
      int sum = 0;
      ptrFactor1++;
      ptrFactor2++;
      if (modulus < 32768)
      {
#ifdef _USING64BITS_
        uint64_t dSum = 0;
#else
        double dSum = 0;
#endif
        for (; j >= 3; j-=4)
        {
#ifdef _USING64BITS_
          dSum += (uint64_t)(*ptrFactor1 * *ptrFactor2) +
            (uint64_t)(*(ptrFactor1 + 2) * *(ptrFactor2 - 2)) +
            (uint64_t)(*(ptrFactor1 + 4) * *(ptrFactor2 - 4)) +
            (uint64_t)(*(ptrFactor1 + 6) * *(ptrFactor2 - 6));
#else
          dSum += (double)(*ptrFactor1 * *ptrFactor2) +
            (double)(*(ptrFactor1 + 2) * *(ptrFactor2 - 2)) +
            (double)(*(ptrFactor1 + 4) * *(ptrFactor2 - 4)) +
            (double)(*(ptrFactor1 + 6) * *(ptrFactor2 - 6));
#endif
          ptrFactor1 += 8;
          ptrFactor2 -= 8;
        }
        while (j >= 0)
        {
#ifdef _USING64BITS_
          dSum += (uint64_t)(*ptrFactor1 * *ptrFactor2);
#else
          dSum += (double)(*ptrFactor1 * *ptrFactor2);
#endif
          ptrFactor1 += 2;
          ptrFactor2 -= 2;
          j--;
        }
#ifdef _USING64BITS_
        sum = (int)(dSum % modulus);
#else
        sum = (int)(dSum - floor(dSum / modulus) * modulus);
#endif
      }
      else
      {
        for (; j >= 0; j--)
        {
#ifdef _USING64BITS_
          sum += (int)((int64_t)*ptrFactor1 * *ptrFactor2 % modulus) - modulus;
#else
          smallmodmult(*ptrFactor1, *ptrFactor2, &result, modulus);
          sum += result.x - modulus;
#endif
          // If sum < 0 do sum <- sum + modulus else do nothing.
          sum += modulus & (sum >> BITS_PER_GROUP);
          ptrFactor1 += 2;
          ptrFactor2 -= 2;
        }
      }
      coeff[i].limbs[0].x = sum;
    }
    else
    {          // General case.
      for (; j >= 0; j--)
      {
        LenAndLimbs2ArrLimbs(ptrFactor1, operand3.limbs, nbrLimbs);
        LenAndLimbs2ArrLimbs(ptrFactor2, operand2.limbs, nbrLimbs);
        modmult(operand3.limbs, operand2.limbs, operand3.limbs);
        AddBigNbrMod(coeff[i].limbs, operand3.limbs, coeff[i].limbs);
        ptrFactor1 += nbrLimbs;
        ptrFactor2 -= nbrLimbs;
      }
    }
  }
  ptrFactor1 = &polyMultTemp[idxFactor1*nbrLimbs];
  if (nbrLimbs == 2)
  {    // Optimization for the case when there is only one limb.
    for (i = 0; i < 2 * coeffLen - 1; i++)
    {
      *ptrFactor1++ = 1;
      *ptrFactor1++ = coeff[i].limbs[0].x;
    }
  }
  else
  {
    for (i = 0; i < 2 * coeffLen - 1; i++)
    {
      ArrLimbs2LenAndLimbs(ptrFactor1, coeff[i].limbs, nbrLimbs);
      ptrFactor1 += nbrLimbs;
    }
  }
  *ptrFactor1 = 1;
  *(ptrFactor1+1) = 0;
  return;
}

static struct stKaratsubaStack
{
  int idxFactor1;
  int stage;
} astKaratsubaStack[10];

// Recursive Karatsuba function.
static void KaratsubaPoly(int idxFactor1, int nbrLen, int nbrLimbs)
{
  int i, idxFactor2;
  int *ptrResult, *ptrHigh, *ptr1, *ptr2;
  int sum, modulus;
  int halfLength;
  int diffIndex = 2 * nbrLen;
  static struct stKaratsubaStack *pstKaratsubaStack = astKaratsubaStack;
  static int coeff[MAX_LEN];
  int stage = 0;
  // Save current parameters in stack.
  pstKaratsubaStack->idxFactor1 = idxFactor1;
  pstKaratsubaStack->stage = -1;
  pstKaratsubaStack++;
  do
  {
    switch (stage)
    {
    case 0:
      idxFactor2 = idxFactor1 + nbrLen;
      if (nbrLen <= KARATSUBA_POLY_CUTOFF)
      {
        // Check if one of the factors is equal to zero.
        ptrResult = &polyMultTemp[idxFactor1*nbrLimbs];
        for (i = nbrLen; i > 0; i--)
        {
          if (*ptrResult != 1 || *(ptrResult + 1) != 0)
          {      // Coefficient is not zero.
            break;
          }
          ptrResult += nbrLimbs;
        }
        if (i > 0)
        {     // First factor is not zero. Check second.
          ptrResult = &polyMultTemp[idxFactor2*nbrLimbs];
          for (i = nbrLen; i > 0; i--)
          {
            if (*ptrResult != 1 || *(ptrResult + 1) != 0)
            {
              break;
            }
            ptrResult += nbrLimbs;
          }
        }
        if (i == 0)
        {    // One of the factors is equal to zero.
          for (i = nbrLen - 1; i >= 0; i--)
          {
            polyMultTemp[idxFactor1*nbrLimbs] = 1;
            polyMultTemp[idxFactor2*nbrLimbs] = 1;
            polyMultTemp[idxFactor1*nbrLimbs + 1] = 0;
            polyMultTemp[idxFactor2*nbrLimbs + 1] = 0;
            idxFactor1++;
            idxFactor2++;
          }
        }
        else
        {
          // Below cutoff: perform standard classical polynomial multiplcation.
          ClassicalPolyMult(idxFactor1, idxFactor2, nbrLen, nbrLimbs);
        }
        pstKaratsubaStack--;
        idxFactor1 = pstKaratsubaStack->idxFactor1;
        nbrLen *= 2;
        diffIndex -= nbrLen;
        stage = pstKaratsubaStack->stage;
        break;
      }
      // Length > KARATSUBA_CUTOFF: Use Karatsuba multiplication.
      // It uses three half-length multiplications instead of four.
      //  x*y = (xH*b + xL)*(yH*b + yL)
      //  x*y = (b + 1)*(xH*yH*b + xL*yL) + (xH - xL)*(yL - yH)*b
      // The length of b is stored in variable halfLength.

      // At this moment the order is: xL, xH, yL, yH.
      // Exchange high part of first factor with low part of 2nd factor.
      halfLength = nbrLen >> 1;
      for (i = idxFactor1 + halfLength; i < idxFactor2; i++)
      {
        memcpy(coeff, &polyMultTemp[i*nbrLimbs], nbrLimbs * sizeof(int));
        memcpy(&polyMultTemp[i*nbrLimbs], &polyMultTemp[(i + halfLength)*nbrLimbs], nbrLimbs * sizeof(int));
        memcpy(&polyMultTemp[(i + halfLength)*nbrLimbs], coeff, nbrLimbs * sizeof(int));
      }
      // At this moment the order is: xL, yL, xH, yH.
      // Compute (xH-xL) and (yL-yH) and store them starting from index diffIndex.
      ptr1 = &polyMultTemp[idxFactor1*nbrLimbs];
      ptr2 = &polyMultTemp[idxFactor2*nbrLimbs];
      ptrResult = &polyMultTemp[diffIndex*nbrLimbs];
      if (nbrLimbs == 2)
      {    // Small modulus.
        modulus = TestNbr[0].x;
        for (i = 0; i < halfLength; i++)
        {
          sum = *(ptr2 + 1) - *(ptr1 + 1);
          // If sum < 0 do sum <- sum + modulus else do nothing.
          sum += modulus & (sum >> BITS_PER_GROUP);
          *(ptrResult + 1) = sum;
          ptr1 += 2;
          ptr2 += 2;
          ptrResult += 2;
        }
        for (i = 0; i < halfLength; i++)
        {
          sum = *(ptr1 + 1) - *(ptr2 + 1);
          // If sum < 0 do sum <- sum + modulus else do nothing.
          sum += modulus & (sum >> BITS_PER_GROUP);
          *(ptrResult + 1) = sum;
          ptr1 += 2;
          ptr2 += 2;
          ptrResult += 2;
        }
      }
      else
      {    // General case.
        for (i = 0; i < halfLength; i++)
        {
          LenAndLimbs2ArrLimbs(ptr1, operand3.limbs, nbrLimbs);
          LenAndLimbs2ArrLimbs(ptr2, operand2.limbs, nbrLimbs);
          SubtBigNbrMod(operand2.limbs, operand3.limbs, operand3.limbs);
          ArrLimbs2LenAndLimbs(ptrResult, operand3.limbs, nbrLimbs);
          ptr1 += nbrLimbs;
          ptr2 += nbrLimbs;
          ptrResult += nbrLimbs;
        }
        for (i = 0; i < halfLength; i++)
        {
          LenAndLimbs2ArrLimbs(ptr1, operand3.limbs, nbrLimbs);
          LenAndLimbs2ArrLimbs(ptr2, operand2.limbs, nbrLimbs);
          SubtBigNbrMod(operand3.limbs, operand2.limbs, operand3.limbs);
          ArrLimbs2LenAndLimbs(ptrResult, operand3.limbs, nbrLimbs);
          ptr1 += nbrLimbs;
          ptr2 += nbrLimbs;
          ptrResult += nbrLimbs;
        }
      }
      // Save current parameters in stack.
      pstKaratsubaStack->idxFactor1 = idxFactor1;
      pstKaratsubaStack->stage = 1;
      pstKaratsubaStack++;
      // Multiply both low parts.
      diffIndex += nbrLen;
      nbrLen = halfLength;
      break;
    case 1:
      // Multiply both high parts.
      idxFactor1 += nbrLen;
      diffIndex += nbrLen;
      nbrLen >>= 1;
      pstKaratsubaStack->stage = 2;
      pstKaratsubaStack++;
      stage = 0;         // Start new Karatsuba multiplication.
      break;
    case 2:
      // Multiply the differences.
      idxFactor1 = diffIndex;
      diffIndex += nbrLen;
      nbrLen >>= 1;
      pstKaratsubaStack->stage = 3;
      pstKaratsubaStack++;
      stage = 0;         // Start new Karatsuba multiplication.
      break;
    default:
      halfLength = nbrLen >> 1;
         // Obtain (b+1)(xH*yH*b + xL*yL) = xH*yH*b^2 + (xL*yL+xH*yH)*b + xL*yL
         // The first and last terms are already in correct locations.
         // Add (xL*yL+xH*yH)*b.
      ptrResult = &polyMultTemp[(idxFactor1 + halfLength) * nbrLimbs];
      if (nbrLimbs == 2)
      {        // Optimization for small numbers.
        int nbrLen2 = nbrLen * 2;
        modulus = TestNbr[0].x;
        ptrResult++;
        for (i = halfLength; i > 0; i--)
        {
          // First addend is the coefficient from xH*yH*b^2 + xL*yL
          // Second addend is the coefficient from xL*yL
          int coeff = *(ptrResult);
          int coeff1 = *(ptrResult + nbrLen);
          sum = coeff + *(ptrResult - nbrLen) - modulus;
          // If sum < 0 do sum <- sum + modulus else do nothing.
          sum += modulus & (sum >> BITS_PER_GROUP);
          // Addend is the coefficient from xH*yH
          sum += coeff1 - modulus;
          // If sum < 0 do sum <- sum + modulus else do nothing.
          sum += modulus & (sum >> BITS_PER_GROUP);
          // Store coefficient of xH*yH*b^2 + (xL*yL+xH*yH)*b + xL*yL
          *(ptrResult) = sum;

          // First addend is the coefficient from xL*yL
          // Second addend is the coefficient from xH*yH
          sum = coeff + *(ptrResult + nbrLen2) - modulus;
          // If sum < 0 do sum <- sum + modulus else do nothing.
          sum += modulus & (sum >> BITS_PER_GROUP);
          // Addend is the coefficient from xH*yH*b^2 + xL*yL
          sum += coeff1 - modulus;
          // If sum < 0 do sum <- sum + modulus else do nothing.
          sum += modulus & (sum >> BITS_PER_GROUP);
          // Store coefficient of xH*yH*b^2 + (xL*yL+xH*yH)*b + xL*yL
          *(ptrResult + nbrLen) = sum;
          // Point to next address.
          ptrResult += 2;
        }
      }
      else
      {        // General case.
        for (i = halfLength; i > 0; i--)
        {
          // Obtain coefficient from xH*yH*b^2 + xL*yL
          LenAndLimbs2ArrLimbs(ptrResult, operand3.limbs, nbrLimbs);
          // Obtain coefficient from xL*yL
          LenAndLimbs2ArrLimbs(ptrResult - halfLength * nbrLimbs, operand2.limbs, nbrLimbs);
          // Obtain coefficient from xH*yH
          LenAndLimbs2ArrLimbs(ptrResult + halfLength * nbrLimbs, operand1.limbs, nbrLimbs);
          // Add all three coefficients.
          AddBigNbrMod(operand3.limbs, operand2.limbs, operand2.limbs);
          AddBigNbrMod(operand2.limbs, operand1.limbs, operand2.limbs);
          // Store coefficient of xH*yH*b^2 + (xL*yL+xH*yH)*b + xL*yL
          ArrLimbs2LenAndLimbs(ptrResult, operand2.limbs, nbrLimbs);
          // Obtain coefficient from xH*yH
          LenAndLimbs2ArrLimbs(ptrResult + nbrLen * nbrLimbs, operand2.limbs, nbrLimbs);
          // Add coefficient from xL*yL
          AddBigNbrMod(operand3.limbs, operand2.limbs, operand3.limbs);
          // Add coefficient from xH*yH*b^2 + xL*yL
          AddBigNbrMod(operand3.limbs, operand1.limbs, operand3.limbs);
          // Store coefficient of xH*yH*b^2 + (xL*yL+xH*yH)*b + xL*yL
          ArrLimbs2LenAndLimbs(ptrResult + halfLength * nbrLimbs, operand3.limbs, nbrLimbs);
          // Point to next address.
          ptrResult += nbrLimbs;
        }
      }
      // Compute final product by adding (xH - xL)*(yL - yH)*b.
      ptrHigh = &polyMultTemp[diffIndex*nbrLimbs];
      ptrResult = &polyMultTemp[(idxFactor1 + halfLength)*nbrLimbs];
      if (nbrLimbs == 2)
      {        // Optimization for small numbers.
        modulus = TestNbr[0].x;
        for (i = nbrLen; i >= 2; i -= 2)
        {
          sum = *(ptrResult + 1) + *(ptrHigh + 1) - modulus;
          // If sum < 0 do sum <- sum + modulus else do nothing.
          sum += modulus & (sum >> BITS_PER_GROUP);
          *(ptrResult + 1) = sum;
          sum = *(ptrResult + 3) + *(ptrHigh + 3) - modulus;
          // If sum < 0 do sum <- sum + modulus else do nothing.
          sum += modulus & (sum >> BITS_PER_GROUP);
          *(ptrResult + 3) = sum;
          ptrHigh += 4;
          ptrResult += 4;
        }
        if (i > 0)
        {
          sum = *(ptrResult + 1) + *(ptrHigh + 1) - modulus;
          // If sum < 0 do sum <- sum + modulus else do nothing.
          sum += modulus & (sum >> BITS_PER_GROUP);
          *(ptrResult + 1) = sum;
        }
      }
      else
      {        // General case.
        for (i = nbrLen; i > 0; i--)
        {
          LenAndLimbs2ArrLimbs(ptrResult, operand3.limbs, nbrLimbs);
          LenAndLimbs2ArrLimbs(ptrHigh, operand2.limbs, nbrLimbs);
          AddBigNbrMod(operand3.limbs, operand2.limbs, operand3.limbs);
          ArrLimbs2LenAndLimbs(ptrResult, operand3.limbs, nbrLimbs);
          ptrHigh += nbrLimbs;
          ptrResult += nbrLimbs;
        }
      }
      nbrLen *= 2;
      diffIndex -= nbrLen;
      pstKaratsubaStack--;
      idxFactor1 = pstKaratsubaStack->idxFactor1;
      stage = pstKaratsubaStack->stage;
    }     // End switch
  } while (stage >= 0);
}

// Multiply factor1 by factor2.The result will be stored in polyMultTemp.
static void MultIntegerPolynomial(int degree1, int degree2,
                                  /*@in@*/int *factor1, /*@in@*/int *factor2)
{
  int indexes[2][MAX_DEGREE];
  int *ptrIndex, *piTemp;
  int currentDegree, index, degreeF2;
  int *piDest;
  
  if (degree1 < degree2)
  {       // Force degree1 >= degree2.
    int tmp = degree1;
    degree1 = degree2;
    degree2 = tmp;
    piTemp = factor1;
    factor1 = factor2;
    factor2 = piTemp;
  }
  // Fill indexes to start of each coefficient.
  ptrIndex = &indexes[0][0];
  index = 0;
  for (currentDegree = 0; currentDegree <= degree1; currentDegree++)
  {
    *ptrIndex++ = index;
    index += numLimbs(factor1 + index) + 1;
  }
  ptrIndex = &indexes[1][0];
  index = 0;
  for (currentDegree = 0; currentDegree <= degree2; currentDegree++)
  {
    *ptrIndex++ = index;
    index += numLimbs(factor2 + index) + 1;
  }
  piDest = polyMultTemp;
  for (currentDegree = 0; currentDegree < degree2; currentDegree++)
  {
    intToBigInteger(&operand4, 0);
    for (degreeF2 = 0; degreeF2 <= currentDegree; degreeF2++)
    {
      UncompressBigIntegerB(factor2 + indexes[1][degreeF2], &operand1);
      UncompressBigIntegerB(factor1 + indexes[0][currentDegree - degreeF2], &operand2);
      BigIntMultiply(&operand1, &operand2, &operand3);
      BigIntAdd(&operand4, &operand3, &operand4);
      NumberLength = operand4.nbrLimbs;
    }
    BigInteger2IntArray(piDest, &operand4);
    piDest += 1 + numLimbs(piDest);
  }
  for (; currentDegree <= degree1; currentDegree++)
  {
    intToBigInteger(&operand4, 0);
    for (degreeF2 = 0; degreeF2 <= degree2; degreeF2++)
    {
      UncompressBigIntegerB(factor2 + indexes[1][degreeF2], &operand1);
      UncompressBigIntegerB(factor1 + indexes[0][currentDegree - degreeF2], &operand2);
      BigIntMultiply(&operand1, &operand2, &operand3);
      BigIntAdd(&operand4, &operand3, &operand4);
      NumberLength = operand4.nbrLimbs;
    }
    BigInteger2IntArray(piDest, &operand4);
    piDest += 1 + numLimbs(piDest);
  }
  for (; currentDegree <= degree1 + degree2; currentDegree++)
  {
    intToBigInteger(&operand4, 0);
    for (degreeF2 = currentDegree - degree1; degreeF2 <= degree2; degreeF2++)
    {
      UncompressBigIntegerB(factor2 + indexes[1][degreeF2], &operand1);
      UncompressBigIntegerB(factor1 + indexes[0][currentDegree - degreeF2], &operand2);
      BigIntMultiply(&operand1, &operand2, &operand3);
      BigIntAdd(&operand4, &operand3, &operand4);
      NumberLength = operand4.nbrLimbs;
    }
    BigInteger2IntArray(piDest, &operand4);
    piDest += 1 + numLimbs(piDest);
  }
}

// Multiply factor1 by factor2. The result will be stored in polyMultTemp.
void MultPolynomial(int degree1, int degree2, /*@in@*/int *factor1, /*@in@*/int *factor2)
{
  int currentDegree;
  int *ptrValue1;
  if (modulusIsZero)
  {
    MultIntegerPolynomial(degree1, degree2, factor1, factor2);
    return;
  }
  int nbrLimbs = NumberLength + 1;
  // Find the least power of 2 greater or equal than the maximum of factor1 and factor2.
  int karatDegree = (degree1 > degree2? degree1: degree2) + 1;
  // Compute length of numbers for each recursion.
  if (karatDegree > KARATSUBA_POLY_CUTOFF)
  {
    int div = 1;
    while (karatDegree > KARATSUBA_POLY_CUTOFF)
    {
      div *= 2;
      karatDegree = (karatDegree + 1) / 2;
    }
    karatDegree *= div;
  }
  // Initialize Karatsuba polynomial.
  ptrValue1 = polyMultTemp;
  for (currentDegree = 2 * karatDegree; currentDegree > 0; currentDegree--)
  {
    *ptrValue1 = 1;        // Initialize coefficient to zero.
    *(ptrValue1 + 1) = 0;
    ptrValue1 += nbrLimbs;
  }
  memset(polyMultTemp, 0, 3 * karatDegree * nbrLimbs * sizeof(limb));
  memcpy(polyMultTemp, factor1, (degree1+1)*nbrLimbs * sizeof(limb));
  memcpy(&polyMultTemp[karatDegree*nbrLimbs], factor2, (degree2+1)*nbrLimbs * sizeof(limb));
  KaratsubaPoly(0, karatDegree, nbrLimbs);
}

int *CopyPolyProduct(int *ptrSrc, int *ptrDest, int degree)
{
  int currentDegree;
  int *ptrValue1 = ptrSrc;
  int *ptrValue2 = ptrDest;
  for (currentDegree = 0; currentDegree <= degree; currentDegree++)
  {
    int nbrLength = 1 + numLimbs(ptrValue1);
    memcpy(ptrValue2, ptrValue1, nbrLength * sizeof(int));
    ptrValue1 += nbrLength;
    ptrValue2 += nbrLength;
  }
  return ptrValue1;
}

static int MultPolynomialExpr(int *ptrArgument1, int *ptrArgument2)
{
  int degreeMono, degreePoly;
  int nbrLimbs = powerMod.nbrLimbs + 1;
  int *ptrValue1, *ptrValue2;
  int currentDegree;
  int degree1 = *ptrArgument1;
  int degree2 = *ptrArgument2;
  if (degree1 <= 0 && degree2 <= 0)
  {        // Product of two monomials.
    if (degree1 + degree2 < -MAX_DEGREE)
    {
      return EXPR_DEGREE_TOO_HIGH;
    }
    *ptrArgument1 = degree1 + degree2;
    UncompressBigIntegerB(ptrArgument1 + 1, &operand1);
    UncompressBigIntegerB(ptrArgument2 + 1, &operand2);
    if (modulusIsZero)
    {
      BigIntMultiply(&operand1, &operand2, &operand1);
      NumberLength = operand1.nbrLimbs;
    }
    else
    {
      modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    }
    BigInteger2IntArray(ptrArgument1 + 1, &operand1);
    valuesIndex = (int)(ptrArgument1+2+*(ptrArgument1+1) - &values[0]);
    return EXPR_OK;
  }
  if (degree1 > 0 && degree2 > 0)
  {        // Product of two polynomials.
    if (degree1 + degree2 > MAX_DEGREE)
    {
      return EXPR_DEGREE_TOO_HIGH;
    }
    *ptrArgument1 = degree1 + degree2;
    if (modulusIsZero)
    {
      CopyPolyProduct(ptrArgument1 + 1, poly1, degree1);  // Copy first factor to poly1
      CopyPolyProduct(ptrArgument2 + 1, poly2, degree2);  // Copy second factor to poly2
    }
    else
    {
      // Copy first factor to poly1
      ptrValue1 = ptrArgument1 + 1;
      ptrValue2 = poly1;
      for (currentDegree = 0; currentDegree <= degree1; currentDegree++)
      {
        memcpy(ptrValue2, ptrValue1, (1 + *ptrValue1) * sizeof(int));
        ptrValue1 += 1 + *ptrValue1;
        ptrValue2 += nbrLimbs;
      }
      // Copy second factor to poly2
      ptrValue1 = ptrArgument2 + 1;
      ptrValue2 = poly2;
      for (currentDegree = 0; currentDegree <= degree2; currentDegree++)
      {
        memcpy(ptrValue2, ptrValue1, (1 + *ptrValue1) * sizeof(int));
        ptrValue1 += 1 + *ptrValue1;
        ptrValue2 += nbrLimbs;
      }
    }
    // Perform multiplication storing the result in polyMultTemp
    MultPolynomial(degree1, degree2, poly1, poly2);
    // Move product back to values stack.
    if (modulusIsZero)
    {
      ptrValue1 = CopyPolyProduct(polyMultTemp, ptrArgument1 + 1, degree1 + degree2);
    }
    else
    {
      ptrValue1 = ptrArgument1 + 1;
      ptrValue2 = polyMultTemp;
      for (currentDegree = 0; currentDegree <= degree1 + degree2; currentDegree++)
      {
        memcpy(ptrValue1, ptrValue2, (1 + *ptrValue2) * sizeof(int));
        ptrValue1 += 1 + *ptrValue1;
        ptrValue2 += nbrLimbs;
      }
    }
    valuesIndex = (int)(ptrValue1 - &values[0]);
    return EXPR_OK;
  }
  // Product of monomial and polynomial.
  if (degree1 <= 0)
  {      // First factor is monomial.
    degreeMono = -degree1;
    degreePoly = degree2;
         // Point to first coefficient of polynomial.
    ptrValue1 = ptrArgument2 + 1;
         // Get coefficient of monomial.
    UncompressBigIntegerB(ptrArgument1 + 1, &operand1);
  }
  else
  {      // Second factor is monomial.
    degreeMono = -degree2;
    degreePoly = degree1;
         // Point to first coefficient of polynomial.
    ptrValue1 = ptrArgument1 + 1;
         // Get coefficient of monomial.
    UncompressBigIntegerB(ptrArgument2 + 1, &operand1);
  }
  if (degreeMono + degreePoly > MAX_DEGREE)
  {
    return EXPR_DEGREE_TOO_HIGH;
  }
     // Multiply all coefficients of polynomial by the coefficient
     // of monomial storing the resulting polynomial on poly1.
  ptrValue2 = poly1;
  for (currentDegree = 0; currentDegree <= degreePoly; currentDegree++)
  {
    UncompressBigIntegerB(ptrValue1, &operand2);
    if (modulusIsZero)
    {
      BigIntMultiply(&operand1, &operand2, &operand3);
      NumberLength = operand3.nbrLimbs;
    }
    else
    {
      modmult(operand1.limbs, operand2.limbs, operand3.limbs);
    }
    BigInteger2IntArray(ptrValue2, &operand3);
    ptrValue1 += 1 + *ptrValue1;
    ptrValue2 += 1 + *ptrValue2;
  }
  *ptrArgument1 = degreeMono + degreePoly;
  ptrValue1 = ptrArgument1 + 1;
  for (currentDegree = 0; currentDegree < degreeMono; currentDegree++)
  {
    *ptrValue1++ = 1;
    *ptrValue1++ = 0;
  }
  memcpy(ptrValue1, poly1, (ptrValue2 - &poly1[0])*sizeof(int));
  valuesIndex = (int)(ptrValue1 - &values[0] + (ptrValue2 - &poly1[0]));
  return EXPR_OK;
}

static void ToPoly(int polyDegree, int *polySrc, int *polyDest)
{
  int currentDegree;
  polySrc++;
  if (polyDegree < 0)
  {    // Polynomial is a monomial
    for (currentDegree = 0; currentDegree < -polyDegree; currentDegree++)
    {
      *polyDest = 1;
      *(polyDest + 1) = 0;
      polyDest += NumberLength + 1;
    }
    memcpy(polyDest, polySrc, (*(polySrc)+1)*sizeof(int));
  }
  else
  {   // Polynomial
    for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
    {
      int nbrLimbs = *(polySrc)+1;
      memcpy(polyDest, polySrc, nbrLimbs*sizeof(int));
      polyDest += NumberLength + 1;
      polySrc += nbrLimbs;
    }
  }
}

static void FromPoly(int polyDegree, int *polyDest, int *polySrc)
{
  int currentDegree;
  *polyDest++ = polyDegree;
  for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {
    int nbrLimbs = *(polySrc)+1;
    memcpy(polyDest, polySrc, nbrLimbs*sizeof(int));
    polyDest += NumberLength + 1;
    polySrc += nbrLimbs;
  }
}

static void ReversePolynomial(int *ptrDest, int *ptrSrc)
{
  int indexes[MAX_DEGREE];
  int *ptrIndex;
  int index, numLength;
  int degreePoly = *ptrSrc;
  if (degreePoly < 0)
  {    // Monomial.
    degreePoly = -degreePoly;
    *ptrDest++ = degreePoly;
    numLength = numLimbs(ptrSrc+1) + 1;
    memcpy(ptrDest, ptrSrc + 1, numLength * sizeof(int));
    ptrDest += numLength;
    for (degree = 0; degree < degreePoly; degree++)
    {
      *ptrDest++ = 1;   // Set coefficient to zero.
      *ptrDest++ = 0;
    }
    return;
  }
  // Fill indexes to start of each coefficient.
  ptrIndex = &indexes[0];
  index = 1;
  for (degree = 0; degree <= degreePoly; degree++)
  {
    *ptrIndex++ = index;
    index += numLimbs(ptrSrc + index) + 1;
  }
  // Copy to destination.
  *ptrDest++ = degreePoly;
  for (degree = degreePoly; degree >= 0; degree--)
  {
    int *ptrSrcCoeff = ptrSrc + indexes[degree];
    numLength = numLimbs(ptrSrcCoeff) + 1;
    memcpy(ptrDest, ptrSrcCoeff, numLength * sizeof(int));
    ptrDest += numLength;
  }
}

int DivideIntegerPolynomial(int *pDividend, int *pDivisor, enum eDivType type)
{
  int *ptrResult;
  int degreeDividend;
  int degreeDivisor;
  int degreeQuotient;
  int* ptrQuotient;
  if (*pDividend < *pDivisor)
  {      // Degree of dividend is less than degree of divisor.
    if (type == TYPE_DIVISION)
    {    // Get pointer to quotient.
      poly3[0] = 0;     // Degree of quotient is zero.
      poly3[1] = 1;     // Coefficient is zero.
      poly3[2] = 0;
    }
    else
    {                   // Remainder is equal to dividend.
      CopyPolynomial(poly1, pDividend, *pDividend);
    }
    return EXPR_OK;
  }
  // Move arguments to temporary storage with most significant coefficient
  // first.
  ReversePolynomial(poly1, pDividend);
  ReversePolynomial(poly2, pDivisor);
  degreeDividend = poly1[0];
  degreeDivisor = poly2[0];
  ptrQuotient = poly3;
  *ptrQuotient++ = degreeDividend - degreeDivisor;
  for (degreeQuotient = degreeDividend - degreeDivisor;
    degreeQuotient >= 0; degreeQuotient--)
  {
    int *ptrDividend = &poly1[1];
    int *ptrDivisor = &poly2[1];
    int *ptrRemainder = &poly4[1];
    UncompressBigIntegerB(ptrDividend, &operand1);
    UncompressBigIntegerB(ptrDivisor, &operand2);
    BigIntRemainder(&operand1, &operand2, &operand3);
    if (!BigIntIsZero(&operand3))
    {
      return EXPR_POLYNOMIAL_DIVISION_NOT_INTEGER;
    }
    BigIntDivide(&operand1, &operand2, &operand3);
    NumberLength = operand3.nbrLimbs;
    BigInteger2IntArray(ptrQuotient, &operand3);
    ptrQuotient += 1 + NumberLength;
    // Calculate remainder.
    if (degreeDivisor == 0)
    {     // Strip leading coefficient of dividend.
      int numLength = 1 + numLimbs(ptrDividend);
      ptrRemainder = &poly1[1];

      for (degree = degreeDividend; degree > 0; degree--)
      {
        ptrDividend += numLength;
        numLength = 1 + numLimbs(ptrDividend);
        memcpy(ptrRemainder, ptrDividend, numLength * sizeof(int));
        ptrRemainder += numLength;
      }
    }
    else
    {
      for (degree = degreeDivisor; degree > 0; degree--)
      {
        ptrDividend += numLimbs(ptrDividend) + 1;
        ptrDivisor += numLimbs(ptrDivisor) + 1;
        UncompressBigIntegerB(ptrDivisor, &operand2);
        BigIntMultiply(&operand2, &operand3, &operand2);
        UncompressBigIntegerB(ptrDividend, &operand1);
        BigIntSubt(&operand1, &operand2, &operand2);
        NumberLength = operand2.nbrLimbs;
        BigInteger2IntArray(ptrRemainder, &operand2);
        ptrRemainder += 1 + NumberLength;
      }
      // Copy least significant coefficients of dividend into remainder.
      for (degree = degreeDividend - degreeDivisor; degree > 0; degree--)
      {
        ptrDividend += 1 + numLimbs(ptrDividend);
        int numLength = 1 + numLimbs(ptrDividend);
        memcpy(ptrRemainder, ptrDividend, numLength * sizeof(int));
        ptrRemainder += numLength;
      }
      // Copy remainder to dividend.
      memcpy(&poly1[1], &poly4[1], (char *)ptrRemainder - (char *)&poly4[1]);
    }
    degreeDividend--;
  }
  if (type == TYPE_DIVISION)
  {    // Get pointer to quotient.
    ptrResult = poly3;
    degree = poly1[0] - degreeDivisor;
  }
  else
  {    // Get pointer to remainder.
    ptrResult = poly1;
    degree = degreeDivisor - 1;
  }
       // Compute degree discarding leading coefficients set to zero.
  while (degree > 0)
  {
    if (*(ptrResult + 1) != 1 || *(ptrResult + 2) != 0)
    {            // Coefficient is not zero.
      break;
    }
    ptrResult += 2;
    degree--;
  }
  *ptrResult = degree;
       // Copy result to first parameter.
  ReversePolynomial(pDividend, ptrResult);
  // Discard most significant 
  return EXPR_OK;
}

static int DivPolynomialExpr(int *ptrArgument1, int *ptrArgument2, enum eDivType type)
{
  int currentDegree;
  int degree1 = *ptrArgument1;
  int degree2 = *ptrArgument2;
  if (*ptrArgument2 == 0 && *(ptrArgument2 + 1) == 1 && *(ptrArgument2 + 2) == 0)
  {        // Divisor is zero
    if (type == TYPE_DIVISION)
    {
      return EXPR_DIVIDE_BY_ZERO;
    }
    return EXPR_OK;   // a mod 0 = a.
  }
  if (degree1 <= 0 && degree2 <= 0)
  {        // Division of two monomials.
    if (degree1 > degree2)
    {      // Degree of dividend less than degree of divisor.
      if (type == TYPE_DIVISION)
      {       // Result is zero.
        *ptrArgument1 = 0;
        *(ptrArgument1 + 1) = 1;
        *(ptrArgument1 + 2) = 0;
      }
      return EXPR_OK;
    }
    if (type == TYPE_MODULUS)
    {       // Result is zero.
      *ptrArgument1 = 0;
      *(ptrArgument1 + 1) = 1;
      *(ptrArgument1 + 2) = 0;
      return EXPR_OK;
    }
    *ptrArgument1 = degree1 - degree2;
    UncompressBigIntegerB(ptrArgument1 + 1, &operand1);
    UncompressBigIntegerB(ptrArgument2 + 1, &operand2);
    if (modulusIsZero)
    {
      BigIntRemainder(&operand1, &operand2, &operand3);
      if (!BigIntIsZero(&operand3))
      {    // Remainder is not zero.
        return EXPR_POLYNOMIAL_DIVISION_NOT_INTEGER;
      }
      BigIntDivide(&operand1, &operand2, &operand3);
      CopyBigInt(&operand1, &operand3);
      NumberLength = operand1.nbrLimbs;
    }
    else
    {
      ModInvBigNbr(operand2.limbs, operand2.limbs, TestNbr, NumberLength);
      modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    }
    BigInteger2IntArray(ptrArgument1 + 1, &operand1);
    valuesIndex = (int)(ptrArgument1 + 2 + *(ptrArgument1 + 1) - &values[0]);
    return EXPR_OK;
  }
  if (modulusIsZero)
  {
    return DivideIntegerPolynomial(ptrArgument1, ptrArgument2, type);
  }
  ToPoly(degree1, ptrArgument1, poly1); // Move dividend to poly1.
  ToPoly(degree2, ptrArgument2, poly2); // Move divisor to poly2.
  if (degree1 < 0)
  {
    degree1 = -degree1;
  }
  if (degree2 < 0)
  {
    degree2 = -degree2;
  }
  if (degree1 < degree2)
  {       // Degree of dividend less than degree of divisor.
    if (type == TYPE_DIVISION)
    {     // Result is zero.
      *ptrArgument1 = 0;
      *(ptrArgument1 + 1) = 1;
      *(ptrArgument1 + 2) = 0;
    }
    return EXPR_OK;
  }
  DividePolynomial(poly1, degree1, poly2, degree2, poly3); // Set poly3 to quotient.
  if (type == TYPE_DIVISION)
  {
    currentDegree = degree1 - degree2;
    *ptrArgument1 = currentDegree;
    FromPoly(currentDegree, ptrArgument1, poly3); // Move dividend to poly1.
  }
  else
  {
    currentDegree = getDegreePoly(poly1, degree2-1);
    *ptrArgument1 = currentDegree;
    FromPoly(currentDegree, ptrArgument1, poly1); // Move modulus to poly1.
  }
  return EXPR_OK;
}

void SetNumberToOne(/*@out@*/int *ptrValue1)
{
  limb *destLimb;
  int ctr;
  *ptrValue1++ = NumberLengthR1;
  destLimb = MontgomeryMultR1;
  for (ctr = 0; ctr < NumberLengthR1; ctr++)
  {
    *ptrValue1++ = (int)((destLimb++)->x);
  }
}

int *CopyPolynomial(int *ptrDest, int *ptrSrc, int degree)
{
  int currentDegree;

  for (currentDegree = 0; currentDegree <= degree; currentDegree++)
  {
    int numLength = numLimbs(ptrSrc) + 1;
    memcpy(ptrDest, ptrSrc, numLength * sizeof(int));
    ptrSrc += numLength;
    ptrDest += numLength;
  }
  return ptrDest;
}

int *CopyPolynomialFixedCoeffSize(int *ptrDest, int *ptrSrc, int degree, int coeffSize)
{
  int currentDegree;

  for (currentDegree = 0; currentDegree <= degree; currentDegree++)
  {
    int numLength = numLimbs(ptrSrc) + 1;
    memcpy(ptrDest, ptrSrc, numLength * sizeof(int));
    ptrSrc += coeffSize;
    ptrDest += coeffSize;
  }
  return ptrDest;
}

static int PowerPolynomialExpr(int *ptrArgument1, int expon)
{
  limb exponLimb;
  int degreePower;
  int mask;
  int *ptrValue1, *ptrValue2;
  int currentDegree;
  int nbrLimbs = powerMod.nbrLimbs + 1;
  int degreeBase = *ptrArgument1;
  if (degreeBase <= 0)
  {              // Monomial.
    if (-degreeBase*expon > MAX_DEGREE)
    {
      return EXPR_DEGREE_TOO_HIGH;
    }
    *ptrArgument1 = degreeBase*expon;
    UncompressBigIntegerB(ptrArgument1+1, &operand1);
    if (expon < 0)
    {
      return EXPR_EXPONENT_NEGATIVE;
    }
    if (modulusIsZero)
    {
      BigIntPowerIntExp(&operand1, expon, &operand2);
      NumberLength = operand2.nbrLimbs;
    }
    else
    {
      exponLimb.x = expon;
      modPowLimb(operand1.limbs, &exponLimb, operand2.limbs);
    }
    BigInteger2IntArray(ptrArgument1 + 1, &operand2);
    valuesIndex = (int)(ptrArgument1+2+numLimbs(ptrArgument1+1) - &values[0]);
    return EXPR_OK;
  }
                // Polynomial.
  ptrValue1 = ptrArgument1 + 1;
  ptrValue2 = poly1;
  if (!modulusIsZero)
  {
    // Copy base to poly1
    for (currentDegree = 0; currentDegree <= degreeBase; currentDegree++)
    {
      memcpy(ptrValue2, ptrValue1, (1 + *ptrValue1) * sizeof(int));
      ptrValue1 += 1 + *ptrValue1;
      ptrValue2 += nbrLimbs;
    }
  }
  else
  {
    CopyPolynomial(ptrValue2, ptrValue1, degreeBase);
  }
  SetNumberToOne(&poly2[0]); // Initialize power with polynomial 1.
  degreePower = 0;
  for (mask = 1 << (BITS_PER_GROUP - 1); mask > 0; mask >>= 1)
  {
    // Square polynomial.
    MultPolynomial(degreePower, degreePower, poly2, poly2);
    degreePower <<= 1;
    if (modulusIsZero)
    {
      CopyPolynomial(poly2, polyMultTemp, degreePower);
    }
    else
    {
      memcpy(poly2, polyMultTemp, (degreePower + 1)*nbrLimbs * sizeof(int));
    }
    if ((expon & mask) != 0)
    {
      MultPolynomial(degreeBase, degreePower, poly1, poly2);
      degreePower += degreeBase;
      if (modulusIsZero)
      {
        CopyPolynomial(poly2, polyMultTemp, degreePower);
      }
      else
      {
        memcpy(poly2, polyMultTemp, (degreePower + 1)*nbrLimbs * sizeof(int));
      }
    }
  }
  *ptrArgument1 = degreePower;
  ptrValue1 = ptrArgument1 + 1;
  ptrValue2 = polyMultTemp;
  // Move power back to values stack.
  if (modulusIsZero)
  {
    ptrValue1 = CopyPolynomial(ptrValue1, polyMultTemp, degreePower);
  }
  else
  {
    for (currentDegree = 0; currentDegree <= degreePower; currentDegree++)
    {
      int len = 1 + *ptrValue2;
      memcpy(ptrValue1, ptrValue2, len * sizeof(int));
      ptrValue1 += len;
      ptrValue2 += nbrLimbs;
    }
  }
  valuesIndex = (int)(ptrValue1 - &values[0]);
  return EXPR_OK;
}

void computePower(int expo)
{
  BigIntPowerIntExp(&primeMod, expo, &powerMod);
  memcpy(TestNbr, powerMod.limbs, powerMod.nbrLimbs*sizeof(limb));
  NumberLength = powerMod.nbrLimbs;
  TestNbr[NumberLength].x = 0;
  GetMontgomeryParms(powerMod.nbrLimbs);
}
// Compute polynomial using value array as stack.
// Inside stack:
// Monomial: -degree, nbrLimbs, xxx (limbs)
// Polynomial: degree, coeff deg 0, coeff deg 1, ...
// where coefficient: nbrlimbs, xxx (limbs)
// Exponent: exponent (integer)
int ComputePolynomial(char *input, int expo)
{
  int *stackValues[STACK_OPER_SIZE];
  int stackIndex = 0;
  int insideExpon = FALSE;
  int *ptrValue1, *ptrValue2;
  int len, index, bitNbr, bitCtr, valueNbr;
  char *ptrRPNbuffer;
  int rc, count, val, pwr, expon;
  int currentDegree;
  degree = 1;
  exponentMod = expo;
  // Use operand1 as temporary variable to store the exponent.
  computePower(expo);
  rc = ConvertToReversePolishNotation(input, RPNbuffer);
  if (rc != EXPR_OK)
  {
    return rc;
  }
  currentDegree = 0;
  stackIndex = 0;
  valuesIndex = 0;
  ptrRPNbuffer = RPNbuffer;
  while (*ptrRPNbuffer != '\0')
  {
    switch (*ptrRPNbuffer++)
    {
    case TOKEN_NUMBER:
      stackValues[stackIndex++] = &values[valuesIndex];
      if (insideExpon != 0)
      {                    // Inside exponent: only one limb accepted.
        bitCtr = 0;
        valueNbr = 0;
        for (bitNbr = BITS_PER_GROUP; bitNbr > 0; bitNbr -= 8)
        {
          valueNbr += (int)(unsigned char)*ptrRPNbuffer++ << bitCtr;
          bitCtr += 8;
        }
        values[valuesIndex++] = valueNbr;
        break;
      }
      values[valuesIndex++] = 0;   // Degree.
      len = (int)(unsigned char)*(ptrRPNbuffer)* 256 + (unsigned char)*(ptrRPNbuffer + 1);
      values[valuesIndex++] = len;
      ptrRPNbuffer += 2;
      for (index = 0; index < len; index++)
      {
        bitCtr = 0;
        valueNbr = 0;
        for (bitNbr = BITS_PER_GROUP; bitNbr > 0; bitNbr -= 8)
        {
          valueNbr += (int)(unsigned char)*ptrRPNbuffer++ << bitCtr;
          bitCtr += 8;
        }
        values[valuesIndex++] = valueNbr;
      }
      break;
    case 'x':
      stackValues[stackIndex++] = &values[valuesIndex];
      values[valuesIndex++] = -1;   // Degree of monomial
         // Set coefficient to 1 in Montgomery notation.
      SetNumberToOne(&values[valuesIndex]);
      valuesIndex += 1 + values[valuesIndex];
      break;
    case TOKEN_START_EXPON:
      insideExpon = TRUE;
      break;
    case TOKEN_END_EXPON:
      insideExpon = FALSE;
      break;
    case TOKEN_UNARY_MINUS:
      ptrValue1 = stackValues[stackIndex - 1];
      if (insideExpon != 0)
      {
        *ptrValue1 = -*ptrValue1;
      }
      else
      {
        rc = NegatePolynomialExpr(ptrValue1);
        if (rc != EXPR_OK)
        {
          return rc;
        }
      }
      break;
    case TOKEN_GCD:
      ptrValue2 = stackValues[--stackIndex];
      ptrValue1 = stackValues[stackIndex - 1];
      rc = GcdPolynomialExpr(ptrValue1, ptrValue2);
      if (rc != EXPR_OK)
      {
        return rc;
      }
      break;
    case TOKEN_DER:
      ptrValue1 = stackValues[stackIndex - 1];
      rc = DerPolynomial(ptrValue1);
      if (rc != EXPR_OK)
      {
        return rc;
      }
      break;
    case '+':
      ptrValue2 = stackValues[--stackIndex];
      ptrValue1 = stackValues[stackIndex - 1];
      if (insideExpon != 0)
      {
        *ptrValue1 += *ptrValue2;
        if ((unsigned int)*ptrValue1 >= LIMB_RANGE)
        {
          return EXPR_EXPONENT_TOO_LARGE;
        }
      }
      else
      {
        rc = AddPolynomialExpr(ptrValue1, ptrValue2);
        if (rc != EXPR_OK)
        {
          return rc;
        }
      }
      break;
    case '-':
      ptrValue2 = stackValues[--stackIndex];
      ptrValue1 = stackValues[stackIndex - 1];
      if (insideExpon != 0)
      {
        *ptrValue1 -= *ptrValue2;
        if ((unsigned int)*ptrValue1 >= LIMB_RANGE)
        {
          return EXPR_EXPONENT_TOO_LARGE;
        }
      }
      else
      {
        rc = NegatePolynomialExpr(ptrValue2);
        if (rc != EXPR_OK)
        {
          return rc;
        }
        rc = AddPolynomialExpr(ptrValue1, ptrValue2);
        if (rc != EXPR_OK)
        {
          return rc;
        }
      }
      break;
    case '*':
      ptrValue2 = stackValues[--stackIndex];
      ptrValue1 = stackValues[stackIndex - 1];
      if (insideExpon != 0)
      {
        *ptrValue1 *= *ptrValue2;
        if ((unsigned int)*ptrValue1 >= LIMB_RANGE)
        {
          return EXPR_EXPONENT_TOO_LARGE;
        }
      }
      else
      {
        rc = MultPolynomialExpr(ptrValue1, ptrValue2);
        if (rc != EXPR_OK)
        {
          return rc;
        }
      }
      break;
    case '/':
      ptrValue2 = stackValues[--stackIndex];
      ptrValue1 = stackValues[stackIndex - 1];
      if (insideExpon != 0)
      {
        *ptrValue1 /= *ptrValue2;
      }
      else
      {
        rc = DivPolynomialExpr(ptrValue1, ptrValue2, TYPE_DIVISION);
        if (rc != EXPR_OK)
        {
          return rc;
        }
      }
      break;
    case '%':
      ptrValue2 = stackValues[--stackIndex];
      ptrValue1 = stackValues[stackIndex - 1];
      if (insideExpon != 0)
      {
        *ptrValue1 %= *ptrValue2;
      }
      else
      {
        rc = DivPolynomialExpr(ptrValue1, ptrValue2, TYPE_MODULUS);
        if (rc != EXPR_OK)
        {
          return rc;
        }
      }
      break;
    case '^':
      expon = *stackValues[--stackIndex];
      if (expon < 0)
      {
        return EXPR_EXPONENT_NEGATIVE;
      }
      if (insideExpon != 0)
      {
        val = *stackValues[stackIndex - 1];
        pwr = 1;
        if (pow(val, pwr) >= MAX_VALUE_LIMB)
        {
          return EXPR_EXPONENT_TOO_LARGE;
        }
        for (count = expon; count > 0; count--)
        {
          pwr *= val;
        }
        *stackValues[stackIndex - 1] = pwr;
      }
      else
      {
        rc = PowerPolynomialExpr(stackValues[stackIndex - 1], expon);
        if (rc != EXPR_OK)
        {
          return rc;
        }
      }
    }
  }
  // Adjust degree so the leading coefficient is not zero.
  if (values[0] < 0)
  {  // Monomial
    UncompressBigIntegerB(&values[1], &operand1);
    degree = -values[0];
    values[0] = degree;
    ptrValue1 = &values[1];
    for (currentDegree = 0; currentDegree < degree; currentDegree++)
    {
      *ptrValue1++ = 1;   // Initialize coefficient to zero.
      *ptrValue1++ = 0;
    }
    BigInteger2IntArray(ptrValue1, &operand1);
  }
  return EXPR_OK;
}

// Divide the polynomial by its leading coefficient.
// On output, operand1 is the inverse of the leading coefficient.
void ConvertToMonic(int *poly, int polyDegree)
{
  int nbrLimbs = NumberLength + 1;
  int currentDegree;
  if (NumberLength == 1)
  {         // Modulus size is one limb.
    int inverse = modInv(*(poly + polyDegree * 2 + 1), TestNbr[0].x);
    int *ptrPoly = poly + 1;
    intToBigInteger(&operand1, inverse);
    for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
    {
#ifdef _USING64BITS_
      *ptrPoly = (int)((uint64_t)*ptrPoly * inverse % TestNbr[0].x);
#else
      smallmodmult(*ptrPoly, inverse, (limb *)ptrPoly, TestNbr[0].x);
#endif
      ptrPoly += 2;
    }
  }
  else
  {         // General case.
    IntArray2BigInteger(poly + polyDegree * nbrLimbs, &operand1);
    ModInvBigNbr(operand1.limbs, operand1.limbs, TestNbr, NumberLength);
    for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
    {
      IntArray2BigInteger(poly + currentDegree * nbrLimbs, &operand2);
      modmult(operand1.limbs, operand2.limbs, operand2.limbs);
      BigInteger2IntArray(poly + currentDegree * nbrLimbs, &operand2);
    }
  }
}

// The content is the GCD of all coefficients with the sign equal to the sign of the leading coefficient.
int *getContent(int *poly, BigInteger *content)
{
  int currentDegree = *poly++;
  if (currentDegree < 0)
  {
    currentDegree = 0;    // Monomial: process only one coefficient.
  }
  UncompressBigIntegerB(poly, content);
  for (; currentDegree > 0; currentDegree--)
  {
    poly += 1 + numLimbs(poly);
    UncompressBigIntegerB(poly, &operand1);
    BigIntGcd(content, &operand1, &operand2);
    CopyBigInt(content, &operand2);
    content->sign = operand1.sign;
  }
  return poly;
}

// Get polynomial divided content mod prime.
int getModPolynomial(int *polyMod, int *poly, BigInteger *content)
{
  int currentDegree;
  int degreePoly = *poly++;
  if (degreePoly < 0)
  {     // Monomial.
    for (currentDegree = degreePoly; currentDegree < 0; currentDegree++)
    {
      *polyMod++ = 1;       // Initialize coefficient to zero.
      *polyMod++ = 0;
    }
    UncompressBigIntegerB(poly, &operand1);  // Get coefficient.
    BigIntDivide(&operand1, content, &operand2);
    BigIntRemainder(&operand2, &powerMod, &operand1);
    if (operand1.sign == SIGN_NEGATIVE)
    {
      BigIntAdd(&operand1, &powerMod, &operand1);
    }
    NumberLength = operand1.nbrLimbs;
    BigInteger2IntArray(polyMod, &operand1);
    return -degreePoly;
  }
  for (currentDegree = 0; currentDegree <= degreePoly; currentDegree++)
  {
    UncompressBigIntegerB(poly, &operand1);  // Get coefficient.
    BigIntDivide(&operand1, content, &operand2);
    BigIntRemainder(&operand2, &powerMod, &operand1);
    if (operand1.sign == SIGN_NEGATIVE)
    {
      BigIntAdd(&operand1, &powerMod, &operand1);
    }
    NumberLength = operand1.nbrLimbs;
    BigInteger2IntArray(polyMod, &operand1);
    poly += 1 + numLimbs(poly);
    polyMod += 1 + numLimbs(polyMod);
  }
  return degreePoly;
}

void PolynomialGcd(int *argF, int *argG, int *gcd)
{
  static BigInteger contentF, contentG, contentH;
  static BigInteger modulus, gcdLeadingCoeff;
  int *ptrLeadingCoeff;
  int *ptrSrc, *ptrDest, *ptrPrev;
  int potentialDegreeGcd, tmp;
  int rc;
  prime = 65537;
  // Get content (gcd of all coefficients).
  ptrLeadingCoeff = getContent(argF, &contentF);
  UncompressBigIntegerB(ptrLeadingCoeff, &operand3);
  ptrLeadingCoeff = getContent(argG, &contentG);
  UncompressBigIntegerB(ptrLeadingCoeff, &operand2);
  BigIntGcd(&operand2, &operand3, &gcdLeadingCoeff);
  // Get potential degree of GCD as the minimum degree of arguments.
  potentialDegreeGcd = *argF;
  if (potentialDegreeGcd < 0)
  {           // Monomial
    potentialDegreeGcd = -potentialDegreeGcd;
  }
  tmp = *argG;
  if (tmp < 0)
  {           // Monomial
    tmp = -tmp;
  }
  if (potentialDegreeGcd > tmp)
  {
    potentialDegreeGcd = tmp;
  }
  intToBigInteger(&modulus, 1);  // Initialize modulus to 1.
  polyS[0] = 0;                  // Initialize polynomial g_m to zero.
  polyS[1] = 1;
  polyS[2] = 0;
  for (;;)
  {
    int degree1, degree2, degreeGcdMod;
    do
    { // Find next prime.
      prime = nextPrime(prime);
      // Ensure that the prime does not divide the modulus or
      // the gcd of the leading coefficients.
    } while (getRemainder(&modulus, prime) == 0 ||
      getRemainder(&gcdLeadingCoeff, prime) == 0);
    modulusIsZero = 0;
    intToBigInteger(&primeMod, prime);
    computePower(1);
    degree1 = getModPolynomial(poly1, argF, &contentF);
    degree2 = getModPolynomial(poly2, argG, &contentG);
    PolyModularGcd(poly1, degree1, poly2, degree2, poly3, &degreeGcdMod);
    modulusIsZero = 1;
    if (degreeGcdMod == 0)
    {                                // Polynomial GCD is the GCD of both contents.
      BigIntGcd(&contentF, &contentG, &gcdLeadingCoeff);
      NumberLength = gcdLeadingCoeff.nbrLimbs;
      *gcd = 0;                      // Degree of GCD.
      BigInteger2IntArray(gcd + 1, &gcdLeadingCoeff);
      return;
    }
    if (degreeGcdMod > potentialDegreeGcd)
    {                                // Incorrect degree for GCD.
      continue;
    }
    if (degreeGcdMod < potentialDegreeGcd)
    {                                // Potential degree is incorrect.
      intToBigInteger(&modulus, 1);  // Initialize modulus to 1.
      polyS[0] = 0;                  // Initialize polynomial g_m to zero.
      polyS[1] = 1;
      polyS[2] = 0;
      potentialDegreeGcd = degreeGcdMod;
    }
    ptrSrc = &poly3[0];
    if (polyS[0] == 0 && polyS[1] == 1 && polyS[2] == 0)
    {                                // If previous polynomial is zero, g_m <- (b * g_p) mod p.
      ptrDest = &polyS[1];
      for (degree = 0; degree <= potentialDegreeGcd; degree++)
      {
        UncompressBigIntegerB(ptrSrc, &operand1);
        BigIntMultiply(&operand1, &gcdLeadingCoeff, &operand1);
        intToBigInteger(&operand1, getRemainder(&operand1, prime));
        NumberLength = operand1.nbrLimbs;
        BigInteger2IntArray(ptrDest, &operand1);
        ptrSrc += 1 + numLimbs(ptrSrc);
        ptrDest += 1 + numLimbs(ptrDest);
      }
    }
    else
    {                                // Perform Chinese Remainder Theorem between g_p and g_m.
      int gcdLCModPrime;
      // coeff = c = a (mod p), c = b (mod q) => c = jq + b
      // a = jq + b (mod p) => j = (a-b)/q = r (mod p)
      // j = pk + r => c = (pk + r)q + b = pqk + rq + b => c = rq + b (mod pq).
      polyT[0] = potentialDegreeGcd;
      ptrPrev = &polyS[1];
      ptrDest = &polyT[1];
      gcdLCModPrime = getRemainder(&gcdLeadingCoeff, prime);
      for (degree = 0; degree <= potentialDegreeGcd; degree++)
      {
        int valueR;
        int invModPrime;
        int coeffModPrime = *(ptrSrc + 1);
        coeffModPrime = (int)(((int64_t)coeffModPrime*gcdLCModPrime)%prime);
        UncompressBigIntegerB(ptrPrev, &operand1);
        valueR = coeffModPrime - getRemainder(&operand1, prime);
        if (valueR < 0)
        {
          valueR += prime;
        }
        invModPrime = modInv(getRemainder(&modulus, prime), prime);
        valueR = (int)(((int64_t)valueR * (int64_t)invModPrime) % prime);
        multint(&operand2, &modulus, valueR);
        BigIntAdd(&operand2, &operand1, &operand2);
        NumberLength = operand2.nbrLimbs;
        BigInteger2IntArray(ptrDest, &operand2);
        ptrSrc += 2;
        ptrPrev += 1 + numLimbs(ptrPrev);
        ptrDest += 1 + numLimbs(ptrDest);
      }
      CopyPolynomial(&polyS[1], &polyT[1], potentialDegreeGcd);
    }
    polyS[0] = potentialDegreeGcd;
    multint(&modulus, &modulus, prime);   // modulus <- modulus*prime.
    CopyBigInt(&operand4, &modulus);
    BigIntDivideBy2(&operand4);           // operand4 <- modulus / 2.
    // Get the primitive part of g_m in polyT.
    getContent(polyS, &contentH);
    polyT[0] = potentialDegreeGcd;
    ptrSrc = &polyS[1];
    ptrDest = &polyT[1];
    for (degree = 0; degree <= potentialDegreeGcd; degree++)
    {
      UncompressBigIntegerB(ptrSrc, &operand1);
      // If coefficient is greater than modulus/2, 
      // let coefficient <- coefficient - modulus.
      BigIntSubt(&operand1, &operand4, &operand2);
      if (operand2.sign == SIGN_POSITIVE)
      {
        BigIntSubt(&operand1, &modulus, &operand1);
      }
      BigIntDivide(&operand1, &contentH, &operand2);
      NumberLength = operand2.nbrLimbs;
      BigInteger2IntArray(ptrDest, &operand2);
      ptrSrc += 1 + numLimbs(ptrSrc);
      ptrDest += 1 + numLimbs(ptrDest);
    }
    // if h divides f and h divides g, then return h.
    if (argF[0] < 0)
    {
      degree = 0;
    }
    else
    {
      degree = argF[0];
    }
    poly5[0] = argF[0];
    CopyPolynomial(&poly5[1], &argF[1], degree);
    rc = DivideIntegerPolynomial(poly5, polyT, TYPE_MODULUS);
    if (rc == EXPR_OK && poly5[0] == 0 && poly5[1] == 1 && poly5[2] == 0)
    {           // Remainder is zero.
      if (argG[0] < 0)
      {
        degree = 0;
      }
      else
      {
        degree = argG[0];
      }
      poly5[0] = argG[0];
      CopyPolynomial(&poly5[1], &argG[1], degree);
      rc = DivideIntegerPolynomial(poly5, polyT, TYPE_MODULUS);
      if (rc == EXPR_OK && poly5[0] == 0 && poly5[1] == 1 && poly5[2] == 0)
      {           // Remainder is zero.
                  // Return h*gcd(contentF, contentG)
        BigIntGcd(&contentF, &contentG, &operand1);
        *gcd = potentialDegreeGcd;
        ptrSrc = &polyT[1];
        ptrDest = gcd + 1;
        for (degree = 0; degree <= potentialDegreeGcd; degree++)
        {
          UncompressBigIntegerB(ptrSrc, &operand2);
          BigIntMultiply(&operand1, &operand2, &operand3);
          NumberLength = operand3.nbrLimbs;
          BigInteger2IntArray(ptrDest, &operand3);
          ptrSrc += 1 + numLimbs(ptrSrc);
          ptrDest += 1 + numLimbs(ptrDest);
        }
        return;
      }
    }
  }
}

// arg1 must not be changed by this routine.
// arg2 can be changed by this routine.
void PolyModularGcd(int *arg1, int degree1, int *arg2, int degree2, int *gcd, int *degreeGcd)
{
  int nbrLimbs = NumberLength + 1;
  int *ptrArgMax = gcd;
  int *ptrArgMin = arg2;
  int *ptrTemp;
  int degreeMax = degree1;
  int degreeMin = degree2;
  int temp, index;
  int currentDegree;
  if (degree2 == 0)
  {
    if (*arg2 == 1 && *(arg2 + 1) == 0)
    {     // Number is zero. GCD is first argument.
      *degreeGcd = degree1;
      memcpy(gcd, arg1, (degree1+1)*nbrLimbs*sizeof(int));
      return;
    }
    *degreeGcd = 0;
    return;
  }
  memcpy(gcd, arg1, (degree1+1)*nbrLimbs*sizeof(int));
  if (degreeMax < degreeMin)
  {
    temp = degreeMax;
    degreeMax = degreeMin;
    degreeMin = temp;
    ptrTemp = ptrArgMax;
    ptrArgMax = ptrArgMin;
    ptrArgMin = ptrTemp;
  }
  for (;;)
  {
    ConvertToMonic(ptrArgMin, degreeMin);
    for (currentDegree = degreeMax; currentDegree >= degreeMin; currentDegree--)
    {          // Replace ptrArgMax by remainder of long
               // division of ptrArgMax / ptrArgMin.
      ptrTemp = ptrArgMax + (currentDegree - degreeMin)*nbrLimbs;
      if (nbrLimbs == 2)
      {        // Modulus size is one limb.
        int modulus = TestNbr[0].x;
        int value = *(ptrArgMax + currentDegree * 2 + 1);
        int *ptrPoly = ptrArgMin + 1;
        ptrTemp++;
        for (index = 0; index <= degreeMin; index++)
        {
#ifdef _USING64BITS_
          temp = *ptrTemp - (int)((uint64_t)*ptrPoly * value % modulus);
#else
          smallmodmult(*ptrPoly, value, (limb *)&temp, modulus);
          temp = *ptrTemp - temp;
#endif
          temp += modulus & (temp >> BITS_PER_GROUP);
          *ptrTemp = temp;
          ptrTemp += 2;
          ptrPoly += 2;
        }
      }
      else
      {        // General case.
        IntArray2BigInteger(ptrArgMax + currentDegree * nbrLimbs, &operand1);
        for (index = 0; index <= degreeMin; index++)
        {
          IntArray2BigInteger(ptrArgMin + index * nbrLimbs, &operand2);
          modmult(operand1.limbs, operand2.limbs, operand2.limbs);
          IntArray2BigInteger(ptrTemp, &operand3);
          SubtBigNbrMod(operand3.limbs, operand2.limbs, operand3.limbs);
          BigInteger2IntArray(ptrTemp, &operand3);
          ptrTemp += nbrLimbs;
        }
      }
    }
    degreeMax = degreeMin - 1;
               // Discard any leading zero coefficient.
    while (degreeMax >= 0)
    {
      ptrTemp = ptrArgMax + degreeMax*nbrLimbs;
      for (index = *ptrTemp++; index > 0; index--)
      {
        if (*ptrTemp++ != 0)
        {      // Coefficient is not zero.
          break;
        }
      }
      if (index > 0)
      {
        break;
      }
      degreeMax--;
    }
    if (degreeMax < 0)
    {          // Polynomial is zero, so Gcd was found.
      if (ptrArgMin != gcd)
      {        // Move gcd to output buffer.
        memcpy(gcd, ptrArgMin, (degreeMin+1)*nbrLimbs*sizeof(int));
      }
      *degreeGcd = degreeMin;
      return;
    }
    temp = degreeMax;
    degreeMax = degreeMin;
    degreeMin = temp;
    ptrTemp = ptrArgMax;
    ptrArgMax = ptrArgMin;
    ptrArgMin = ptrTemp;
  }
}

static int GcdPolynomialExpr(int *ptrArgument1, int *ptrArgument2)
{
  int degreeGcd;
  if (modulusIsZero)
  {    // Integer GCD polynomial.
    PolynomialGcd(ptrArgument1, ptrArgument2, ptrArgument1);
  }
  else
  {    // Modular GCD polynomial.
    int degree1 = *ptrArgument1;
    int degree2 = *ptrArgument2;
    PolyModularGcd(ptrArgument1+1, degree1, ptrArgument2 + 1, degree2,
      poly5, &degreeGcd);
    CopyPolynomial(ptrArgument1 + 1, poly5, degreeGcd);
    *ptrArgument1 = degreeGcd;
  }
  return EXPR_OK;
}

int DerPolynomial(int *ptrArgument)
{
  int *ptrSrc, *ptrDest;
  int nbrLimbs;
  int degreePoly = *ptrArgument;
  int derivDegreePoly = 0;
  if (degreePoly < 0)
  {           // Monomial.
    *ptrArgument = degreePoly + 1;   // Decrement degree.
    if (modulusIsZero)
    {
      UncompressBigIntegerB(ptrArgument + 1, &operand1);
      multint(&operand1, &operand1, -degreePoly);
      NumberLength = operand1.nbrLimbs;
      BigInteger2IntArray(ptrArgument + 1, &operand1);
    }
    else
    {
      nbrLimbs = NumberLength + 1;
      LenAndLimbs2ArrLimbs(ptrArgument + 1, operand1.limbs, nbrLimbs);
      modmultInt(operand1.limbs, -degreePoly, operand1.limbs);
      ArrLimbs2LenAndLimbs(ptrArgument + 1, operand1.limbs, nbrLimbs);
      if (*(ptrArgument + 1) == 1 && *(ptrArgument + 2) == 0)
      {   // Only coefficient is zero. Set degree to zero.
        *ptrArgument = 0;
      }
    }
    return EXPR_OK;
  }
  ptrSrc = ptrArgument + 1;
  ptrDest = &poly1[1];
  for (degree = 1; degree <= degreePoly; degree++)
  {
    ptrSrc += 1 + numLimbs(ptrSrc);
    if (modulusIsZero)
    {
      UncompressBigIntegerB(ptrSrc, &operand1);
      multint(&operand1, &operand1, degree);
      NumberLength = operand1.nbrLimbs;
      BigInteger2IntArray(ptrDest, &operand1);
      derivDegreePoly = degree - 1;
    }
    else
    {
      nbrLimbs = NumberLength + 1;
      LenAndLimbs2ArrLimbs(ptrSrc, operand1.limbs, nbrLimbs);
      modmultInt(operand1.limbs, degree, operand1.limbs);
      ArrLimbs2LenAndLimbs(ptrDest, operand1.limbs, nbrLimbs);
      if (*ptrDest != 1 || *(ptrDest + 1) != 0)
      {
        derivDegreePoly = degree-1;
      }
    }
    ptrDest += 1 + numLimbs(ptrDest);
  }
  CopyPolynomial(ptrArgument + 1, &poly1[1], derivDegreePoly);
  *ptrArgument = derivDegreePoly;
  return EXPR_OK;
}

// Convert from Montgomery notation to standard notation by multiplying by 1.
void polyToStandardNotation(int *nbr, int qtyNbrs)
{
  int currentNbr;
  // Initialize operand2.limbs to 1 (little-endian).
  operand2.limbs[0].x = 1;
  if (NumberLength > 1)
  {
    memset(&operand2.limbs[1], 0, (NumberLength - 1) * sizeof(limb));
  }
  for (currentNbr = 0; currentNbr < qtyNbrs; currentNbr++)
  {
    IntArray2BigInteger(nbr, &operand1);
    modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    BigInteger2IntArray(nbr, &operand1);
    nbr += NumberLength + 1;
  }
}

void polyToMontgomeryNotation(int *nbr, int qtyNbrs)
{
  int currentNbr;

  for (currentNbr = qtyNbrs - 1; currentNbr >= 0; currentNbr--)
  {
    IntArray2BigInteger(nbr, &operand1);
    modmult(operand1.limbs, MontgomeryMultR2, operand1.limbs);
    BigInteger2IntArray(nbr, &operand1);
    nbr += NumberLength + 1;
  }
}

// In this routine, the dividend is replaced by the remainder of the division.
// Input and output coefficients are expressed in Montgomery notation.
void DividePolynomial(/*@in@*/int *pDividend, int dividendDegree,
 /*@in@*/int *pDivisor, int divisorDegree, /*@out@*/int *ptrQuotient)
{
  int currentDegree, index;
  int nbrLimbs = NumberLength + 1;
  int divisorIsOne;
  int *ptrQuot;
  int remainderDegree;
  if (divisorDegree > dividendDegree)
  {    // Quotient is zero.
    *ptrQuotient = 1;
    *(ptrQuotient+1) = 0;
    return;
  }
  remainderDegree = dividendDegree - divisorDegree;
  IntArray2BigInteger(pDivisor + divisorDegree*nbrLimbs, &operand1);
  memcpy(operand5.limbs, operand1.limbs, NumberLength*sizeof(int));
  divisorIsOne = !memcmp(operand1.limbs, MontgomeryMultR1, NumberLength*sizeof(int));
  if (!divisorIsOne)
  {        // Leading coefficient is not 1.
    ConvertToMonic(pDivisor, divisorDegree);
    // operand1 holds the inverse of the leading coefficient of divisor.
    // Multiply dividend by this number.
    for (currentDegree = 0; currentDegree <= dividendDegree; currentDegree++)
    {
      IntArray2BigInteger(pDividend + currentDegree*nbrLimbs, &operand2);
      modmult(operand1.limbs, operand2.limbs, operand2.limbs);
      BigInteger2IntArray(pDividend + currentDegree*nbrLimbs, &operand2);
    }
  }
  ptrQuot = ptrQuotient + (dividendDegree - divisorDegree)*nbrLimbs;
  for (currentDegree = dividendDegree; currentDegree >= divisorDegree; currentDegree--)
  {
    int *ptrDivisor;
    int *ptrDividend = pDividend + currentDegree*nbrLimbs;
    IntArray2BigInteger(ptrDividend, &operand1);
    BigInteger2IntArray(ptrQuot, &operand1);  // Store coefficient of quotient.
    ptrDivisor = pDivisor + divisorDegree*nbrLimbs;
    for (index = 0; index <= divisorDegree; index++)
    {
      IntArray2BigInteger(ptrDivisor, &operand2);
      modmult(operand1.limbs, operand2.limbs, operand2.limbs);
      IntArray2BigInteger(ptrDividend, &operand3);
      SubtBigNbrMod(operand3.limbs, operand2.limbs, operand3.limbs);
      BigInteger2IntArray(ptrDividend, &operand3);
      ptrDividend -= nbrLimbs;
      ptrDivisor -= nbrLimbs;
    }
    ptrQuot -= nbrLimbs;
  }
  if (!divisorIsOne)
  {        // Leading coefficient is not 1.
           // Adjust remainder by multiplying each coefficient by leading
           // coefficient of divisor.
    for (currentDegree = 0; currentDegree <= dividendDegree; currentDegree++)
    {
      IntArray2BigInteger(pDividend + currentDegree*nbrLimbs, &operand2);
      modmult(operand5.limbs, operand2.limbs, operand2.limbs);
      BigInteger2IntArray(pDividend + currentDegree*nbrLimbs, &operand2);
    }
           // Restore divisor.
    for (currentDegree = 0; currentDegree <= divisorDegree; currentDegree++)
    {
      IntArray2BigInteger(pDivisor + currentDegree*nbrLimbs, &operand2);
      modmult(operand5.limbs, operand2.limbs, operand2.limbs);
      BigInteger2IntArray(pDivisor + currentDegree*nbrLimbs, &operand2);
    }
  }
}

// Compute the polynomial polyInv as:
// polyInv = x^(2*polyDegree) / polymod
void GetPolyInvParm(int polyDegree, /*@in@*/int *polyMod)
{
  int* ptrCoeff;
  int deg, outerDeg;
  int nbrLimbs = NumberLength + 1;
  // Point to leading coefficient of polyMod.
  // Initialize polyInv as x^polyDegree - K*x^(polyDegree-1)
  // where K is the coefficient of x^(polyDegree-1) of polyMod.
  SetNumberToOne(&polyInv[polyDegree*nbrLimbs]);
  LenAndLimbs2ArrLimbs(polyMod + (polyDegree-1)*nbrLimbs, operand1.limbs, nbrLimbs);
  memset(operand2.limbs, 0, nbrLimbs * sizeof(limb));
  SubtBigNbrMod(operand2.limbs, operand1.limbs, operand1.limbs);
  ptrCoeff = &polyInv[(polyDegree - 1)*nbrLimbs];
  ArrLimbs2LenAndLimbs(ptrCoeff, operand1.limbs, nbrLimbs);
  for (deg = polyDegree - 2; deg >= 0; deg--)
  {
    ptrCoeff -= nbrLimbs;
    *ptrCoeff = 1;          // Set coefficient to zero.
    *(ptrCoeff + 1) = 0;
  }
  // Use Newton iteration: F_{2n}(x) = F_n(x)*(2-D(x)*F_n(x))
  memcpy(poly4, polyMod, nbrLimbs*polyDegree * sizeof(limb));
  SetNumberToOne(&poly4[nbrLimbs*polyDegree]);
  for (outerDeg = 1; outerDeg < polyDegree; outerDeg *= 2)
  {  // Use poly5 as temporary polynomial
    int* ptrCoeff2;
    MultPolynomial(polyDegree, polyDegree, polyInv, poly4);
    ptrCoeff = poly5;
    ptrCoeff2 = &polyMultTemp[polyDegree*nbrLimbs];
    memset(operand2.limbs, 0, nbrLimbs * sizeof(limb));
    for (deg = 0; deg < polyDegree; deg++)
    {
      LenAndLimbs2ArrLimbs(ptrCoeff2, operand3.limbs, nbrLimbs);
      SubtBigNbrMod(operand2.limbs, operand3.limbs, operand3.limbs);
      ArrLimbs2LenAndLimbs(ptrCoeff, operand3.limbs, nbrLimbs);
      ptrCoeff += nbrLimbs;
      ptrCoeff2 += nbrLimbs;
    }
    SetNumberToOne(&poly5[polyDegree*nbrLimbs]);
    MultPolynomial(polyDegree, polyDegree, polyInv, poly5);
    memcpy(polyInv, &polyMultTemp[polyDegree*nbrLimbs], polyDegree * nbrLimbs* sizeof(limb));
  }
}

// Multiply two polynomials mod polyMod using inverse polynomial polyInv.
// The algorithm is:
// m <- (T*polyInv)/x^polyDegree
// return T - m*polyMod

void multUsingInvPolynomial(/*@in@*/int *polyFact1, /*@in@*/int *polyFact2,
                           /*@out@*/int *polyProduct,
                           int polyDegree, /*@in@*/int *polyMod)
{
  int currentDegree, index;
  int nbrLimbs = NumberLength + 1;
  // Compute T
  MultPolynomial(polyDegree, polyDegree, polyFact1, polyFact2);
  memcpy(polyMultT, polyMultTemp, (2*polyDegree + 1)*nbrLimbs*sizeof(limb));
  // Compute m
  MultPolynomial(polyDegree, polyDegree, &polyMultT[polyDegree*nbrLimbs], polyInv);
  memcpy(polyMultM, &polyMultTemp[polyDegree*nbrLimbs],
         (polyDegree + 1)*nbrLimbs * sizeof(limb));
  // Compute m*polyMod
  MultPolynomial(polyDegree, polyDegree, polyMultM, polyMod);
  // Compute T - mN.
  index = 0;
  for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {
    IntArray2BigInteger(&polyMultTemp[index], &operand1);
    IntArray2BigInteger(&polyMultT[index], &operand2);
    SubtBigNbrMod(operand2.limbs, operand1.limbs, operand1.limbs);
    BigInteger2IntArray(polyProduct, &operand1);
    index += nbrLimbs;
    polyProduct += nbrLimbs;
  }
}

// Multiply two polynomials mod polyMod.
void multPolynomialModPoly(/*@in@*/int *polyFact1, /*@in@*/int *polyFact2,
                           /*@out@*/int *polyProduct,
                           int polyDegree, /*@in@*/int *polyMod)
{
  int index1, index2;
  int nbrLimbs = NumberLength + 1;
  int *ptrPoly1, *ptrPoly2, *ptrPolyTemp;
  // Initialize polyMultTemp with the most significant half of product.
  ptrPolyTemp = polyMultTemp + (polyDegree-1)*nbrLimbs;
  for (index1 = polyDegree-1; index1 >= 0; index1--)
  {
    ptrPoly1 = polyFact1 + index1*nbrLimbs;
    ptrPoly2 = polyFact2 + (polyDegree-1)*nbrLimbs;
    IntArray2BigInteger(ptrPoly1, &operand1);
    IntArray2BigInteger(ptrPoly2, &operand2);
    modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    for (index2 = polyDegree - 2; index2 >= index1; index2--)
    {
      ptrPoly2 -= nbrLimbs;
      ptrPoly1 += nbrLimbs;
      IntArray2BigInteger(ptrPoly1, &operand2);
      IntArray2BigInteger(ptrPoly2, &operand3);
      modmult(operand2.limbs, operand3.limbs, operand2.limbs);
      AddBigNbrMod(operand2.limbs, operand1.limbs, operand1.limbs);
    }
    BigInteger2IntArray(ptrPolyTemp, &operand1);
    ptrPolyTemp -= nbrLimbs;
  }
  // Get remainder of long division by polyMod and append next limbs of the product.
  for (index1 = polyDegree-2; index1 >= 0; index1--)
  {
    ptrPoly2 = &polyMultTemp[(polyDegree - 1)*nbrLimbs];
    // Back up leading coefficient.
    memcpy(ptrPoly2 + nbrLimbs, ptrPoly2, nbrLimbs*sizeof(int));
    IntArray2BigInteger(ptrPoly2, &operand3);
    ptrPoly1 = polyMod+polyDegree*nbrLimbs;
    for (index2 = polyDegree - 2; index2 >= 0; index2--)
    {
      ptrPoly1 -= nbrLimbs;
      ptrPoly2 -= nbrLimbs;
      IntArray2BigInteger(ptrPoly1, &operand1);
      modmult(operand3.limbs, operand1.limbs, operand1.limbs);
      IntArray2BigInteger(ptrPoly2, &operand2);
      SubtBigNbrMod(operand2.limbs, operand1.limbs, operand1.limbs);
      BigInteger2IntArray(ptrPoly2+nbrLimbs, &operand1);
    }
    ptrPoly1 = polyFact1 + index1*nbrLimbs;
    ptrPoly2 = polyFact2;
    IntArray2BigInteger(ptrPoly1, &operand1);
    IntArray2BigInteger(ptrPoly2, &operand2);
    modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    for (index2 = 1; index2 <= index1; index2++)
    {
      ptrPoly1 -= nbrLimbs;
      ptrPoly2 += nbrLimbs;
      IntArray2BigInteger(ptrPoly1, &operand2);
      IntArray2BigInteger(ptrPoly2, &operand3);
      modmult(operand2.limbs, operand3.limbs, operand2.limbs);
      AddBigNbrMod(operand2.limbs, operand1.limbs, operand1.limbs);
    }
    IntArray2BigInteger(polyMod, &operand2);
    IntArray2BigInteger(&polyMultTemp[polyDegree*nbrLimbs], &operand3);
    modmult(operand3.limbs, operand2.limbs, operand2.limbs);
    SubtBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
    BigInteger2IntArray(polyMultTemp, &operand1);
  }
  memcpy(polyProduct, polyMultTemp, polyDegree*nbrLimbs*sizeof(int));
}

  // Perform polyPower <- polyBase ^ expon (mod polyMod)
void powerPolynomial(int *polyBase, int *polyMod, int polyDegree, BigInteger *expon,
                     int *polyPower, powerCback callback)
{
  int mask, index;
  int nbrLimbs = NumberLength + 1;
  int powerIsOne = TRUE;
  int nbrBits = 0;
  int bitCounter;
  // Initialize polyPower to 1.
  memset(polyPower, 0, (polyDegree + 1)*nbrLimbs * sizeof(int));
  for (index = 0; index <= polyDegree; index++)
  {
    *(polyPower + index*nbrLimbs) = 1;
  }
  SetNumberToOne(polyPower);
  index = expon->nbrLimbs - 1;
  bitCounter = (index+1) * BITS_PER_GROUP;
  for (; index >= 0; index--)
  {
    int groupExp = (int)(expon->limbs[index].x);
    for (mask = 1 << (BITS_PER_GROUP - 1); mask > 0; mask >>= 1)
    {
      if (!powerIsOne)
      {
        if (callback)
        {
          callback(100 - 100 * bitCounter/nbrBits);
        }
        multUsingInvPolynomial(polyPower, polyPower, polyPower, polyDegree, polyMod);
      }
      if ((groupExp & mask) != 0)
      {
        multUsingInvPolynomial(polyPower, polyBase, polyPower, polyDegree, polyMod);
        if (powerIsOne)
        {
          nbrBits = bitCounter;
          powerIsOne = FALSE;
        }
      }
      bitCounter--;
    }
  }
}

int getDegreePoly(int *poly, int polyDegree)
{
  if (polyDegree <= 0)
  {
    return 0;
  }
  while (polyDegree > 0)
  {
    int len, index, *ptrTemp;
    ptrTemp = poly + polyDegree*(NumberLength + 1);
    len = *ptrTemp;
    for (index = len; index > 0; index--)
    {
      if (*++ptrTemp != 0)
      {      // Coefficient is not zero.
        return polyDegree;
      }
    }
    polyDegree--;
  }
  return polyDegree;
}

// Square-free factorization (from Wikipedia).
// i <- 1; g <- f';
// if g != 0 then {
//   c <- gcd(f, g);
//   w <- f / c;
//   while w != 1 do {
//     y <- gcd(w, c); z <- w / y;
//     Output(z^i); i <- i + 1;
//     w <- y; c <- c / y
//   }
//   if c != 1 then {
//     c <- c^(1/p);
//     Output(SFF(c)^p)
//   }
// }
// else {
//    f <- f^(1/p);
//    Output(SFF(f)^p)
//  }
//  end.
void SquareFreeFactorization(int polyDegree, int *poly, int expon)
{
  int currentDegree, degreeC, degreeY;
  int index;
  int primeInt = (int)primeMod.limbs[0].x;
  int nbrLimbs = primeMod.nbrLimbs + 1;
  int *ptrValue1, *ptrValue2;
  // Generate derivative in poly1.
  ptrValue1 = poly;
  ptrValue2 = poly1;
  for (currentDegree = 1; currentDegree <= polyDegree; currentDegree++)
  {
    ptrValue1 += nbrLimbs;
    IntArray2BigInteger(ptrValue1, &operand1);
    modmultInt(operand1.limbs, currentDegree, operand1.limbs);
    BigInteger2IntArray(ptrValue2, &operand1);
    ptrValue2 += nbrLimbs;
  }
  // Check whether the derivative is zero.
  for (currentDegree = polyDegree-1; currentDegree >= 0; currentDegree--)
  {
    ptrValue1 = &poly1[currentDegree*nbrLimbs+1];
    for (index = 1; index < nbrLimbs; index++)
    {
      if (*ptrValue1++ != 0)
      {         // Coefficient of derivative is not zero.
        break;
      }
    }
    if (index < nbrLimbs)
    {           // Coefficient of derivative is not zero.
      break;
    }
  }
  if (currentDegree >= 0)
  {             // Derivative is not zero.
    int degreeW;
    int i = 1;

    PolyModularGcd(poly, polyDegree, poly1, currentDegree, poly2, &degreeC); // poly2 = c
    memcpy(polyMultTemp, poly, (polyDegree + 1)*nbrLimbs*sizeof(int));      // Backup poly
    if (degreeC == 0)
    {
      memcpy(poly1, poly, (polyDegree + 1)*nbrLimbs*sizeof(int));           // poly1 = w
    }
    else
    {
      SetNumberToOne(&poly2[degreeC*nbrLimbs]);
      DividePolynomial(polyMultTemp, polyDegree, poly2, degreeC, poly1);    // poly1 = w
    }
    degreeW = polyDegree - degreeC;
    while (degreeW != 0)
    {
      memcpy(polyMultTemp, poly1, (degreeW+1)*nbrLimbs*sizeof(int));        // Backup w.
      PolyModularGcd(poly2, degreeC, poly1, degreeW, poly3, &degreeY);       // poly3 = y
      SetNumberToOne(&poly3[degreeY*nbrLimbs]);
      if (degreeW != degreeY)
      {
        int degreeZ;
        struct sFactorInfo *pstFactorInfo;

        DividePolynomial(polyMultTemp, degreeW, poly3, degreeY, poly1);     // poly1 = z
        // z^i is divisor of the original polynomial.
        degreeZ = degreeW - degreeY;
        for (currentDegree = 0; currentDegree < i*expon; currentDegree++)
        {
          DividePolynomial(ptrOrigPoly, degreeOrigPoly, poly1, degreeZ, polyMultTemp);
          degreeOrigPoly -= degreeZ;
          memcpy(ptrOrigPoly, polyMultTemp, (degreeOrigPoly+1)*nbrLimbs*sizeof(int));
        }
        pstFactorInfo = &factorInfo[nbrFactorsFound++];
        pstFactorInfo -> ptr = ptrOrigPoly;
        pstFactorInfo -> degree = degreeZ;
        pstFactorInfo -> multiplicity = i*expon;
        pstFactorInfo -> expectedDegree = 0;    // Unknown at this moment.
        memcpy(ptrOrigPoly, poly1, degreeZ*nbrLimbs*sizeof(int));
        ptrOrigPoly += (degreeZ+1)*nbrLimbs;
        memcpy(ptrOrigPoly, polyMultTemp, (degreeOrigPoly + 1)*nbrLimbs*sizeof(int));
      }
      i++;
      memcpy(poly1, poly3, (degreeY+1)*nbrLimbs*sizeof(int));             // Copy y to w.
      degreeW = getDegreePoly(poly1, degreeY);
      DividePolynomial(poly2, degreeC, poly1, degreeW, polyMultTemp); // Compute c.
      degreeC -= degreeW;
      memcpy(poly2, polyMultTemp, (degreeC+1)*nbrLimbs*sizeof(int));
    }
    if (degreeC != 0)
    {      // C is a perfect power.
      ptrValue1 = ptrOrigPoly + (degreeOrigPoly + 1)*nbrLimbs;
      ptrValue2 = poly2;
      for (currentDegree = 0; currentDegree <= degreeC; currentDegree += primeInt)
      {
        memcpy(ptrValue1, ptrValue2, nbrLimbs*sizeof(int));
        ptrValue1 += nbrLimbs;
        ptrValue2 += primeInt*nbrLimbs;
      }
      SquareFreeFactorization(degreeC / primeInt, ptrOrigPoly + (degreeOrigPoly + 1)*nbrLimbs, expon * primeInt);
    }
  }
  else
  {           // Derivative is zero.
    ptrValue1 = ptrOrigPoly + (degreeOrigPoly + 1)*nbrLimbs;
    ptrValue2 = poly;
    for (currentDegree = 0; currentDegree <= polyDegree; currentDegree += primeInt)
    {
      memcpy(ptrValue1, ptrValue2, nbrLimbs*sizeof(int));
      ptrValue1 += nbrLimbs;
      ptrValue2 += primeInt*nbrLimbs;
    }
    SquareFreeFactorization(polyDegree / primeInt, ptrOrigPoly + (degreeOrigPoly + 1)*nbrLimbs, expon * primeInt);
  }
}

// Extended GCD algorithm of polynomials a, b:
// The algorithm finds polynomials u and v such that au + bv = g
// deg(u) < deg(b) - deg(g), deg(v) < deg(a) - deg(g).
//
// r <- b, old_r <- a
// for (i <- 1; r != 0; i <- i+1) do
//   q[i] = quotient(old_r, r)
//   (old_r, r) <- (r, old_r - q[i] * r)
// end do
// r <- 0, old_r <- 1
// for (j <- 1; j < i; j <- j+1) do
//   (old_r, r) <- (r, old_r - q[j] * r)
// end do
// u <- r
// r <- 1, old_r <- 0
// for (j <- 1; j < i; j <- j+1) do
//   (old_r, r) <- (r, old_r - q[j] * r)
// end do
// v <- r

// Routine that performs one of the last two loops. It has to be called
// twice: once to compute U and again to compute V.
static int *ComputeUorV(/*@in@*/int *ptrQuotients[], int nbrQuotients,
  /*@out@*/int *ptrP1, /*@out@*/int *ptrP2, /*@out@*/int *pDegreeR)
{
  int *tmpPtr;
  int *ptrR = ptrP1;
  int *ptrOldR = ptrP2;
  int degreeR = 0;
  int degreeOldR = 0;
  int counter, currentDegree;
  int nbrLimbs = NumberLength + 1;
  for (counter = 0; counter < nbrQuotients; counter++)
  {
    int tmpDegree, offset;
    int degreeQ = (int)((ptrQuotients[counter + 1] - ptrQuotients[counter])/nbrLimbs - 1);
    // Multiply ptrR by ptrQuotients[counter]. The result will be stored in polyMultTemp.
    MultPolynomial(degreeR, degreeQ, ptrR, ptrQuotients[counter]);
    // Compute ptrOldR <- ptrOldR - polyMultTemp.
    degreeQ += degreeR; // Degree of product (always greater than degree of old R).
    offset = 0;
    for (currentDegree = 0; currentDegree <= degreeOldR; currentDegree++)
    {
      IntArray2BigInteger(ptrOldR + offset, &operand1);
      IntArray2BigInteger(&polyMultTemp[offset], &operand2);
      SubtBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
      BigInteger2IntArray(ptrOldR + offset, &operand1);
      offset += nbrLimbs;
    }
    for (; currentDegree <= degreeQ; currentDegree++)
    {
      // Set operand1 to zero.
      memset(operand1.limbs, 0, NumberLength*sizeof(limb));
      IntArray2BigInteger(&polyMultTemp[offset], &operand2);
      SubtBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
      operand1.sign = SIGN_POSITIVE;
      BigInteger2IntArray(ptrOldR + offset, &operand1);
      offset += nbrLimbs;
    }
    degreeOldR = getDegreePoly(ptrOldR, degreeQ);
    // Exchange degrees.
    tmpDegree = degreeR;
    degreeR = degreeOldR;
    degreeOldR = tmpDegree;
    // Exchange pointers.
    tmpPtr = ptrR;
    ptrR = ptrOldR;
    ptrOldR = tmpPtr;
  }
  *pDegreeR = degreeR;
  return ptrR;
}
// Given coprime polynomials A and B, find polynomials U and V such that A*U + B*V = 1 (mod powerMod).
// Use Montgomery notation.
static void ExtendedGcdPolynomial(/*@in@*/int *ptrA, int degreeA, /*@in@*/int *ptrB, int degreeB, /*@out@*/int *ptrQ,
  /*@out@*/int *ptrP1, /*@out@*/int *ptrP2, /*@out@*/int *ptrU, /*@out@*/int *pDegreeU, /*@out@*/int *ptrV, /*@out@*/int *pDegreeV)
{
  int degreeFirst, degreeSecond;
  int *ptrQuotient = ptrQ;
  int *ptrQuotients[MAX_DEGREE];
  int *tmpPtr, *ptrR, *ptrOldR;
  int degreeR, degreeOldR;
  int nbrQuotients = 0;
  int nbrLimbs = NumberLength + 1;
  int polyExchanged = FALSE;
  int currentDegree;
  int* ptrCoeff;
  int degreeU;
  int degreeV;

  // Ensure that degree of A is greater than degree of B.
	if (degreeA < degreeB)
	{
    int tmpDegree = degreeA;
    degreeA = degreeB;
    degreeB = tmpDegree;
    tmpPtr = ptrA;
    ptrA = ptrB;
    ptrB = tmpPtr;
    polyExchanged = TRUE;
	}
  ptrR = ptrP1;
  ptrOldR = ptrP2;
  memcpy(ptrOldR, ptrA, (degreeA + 1)*nbrLimbs*sizeof(int));
  memcpy(ptrR, ptrB, (degreeB + 1)*nbrLimbs*sizeof(int));
  degreeR = degreeB;
  degreeOldR = degreeA;
     // Loop that computes all quotients.
     // Since gcd = 1, the loop will finish when the degree of R is zero.
  while (degreeR != 0)
  {  // If R is not constant...
#if DEBUG_HENSEL_LIFTING
    if (ptrOutput2 != NULL)
    {
      strcpy(ptrOutput2, "OldR = ");
      ptrOutput2 += strlen(ptrOutput2);
      showPolynomial(&ptrOutput2, ptrOldR, degreeOldR + 1, 0);
      strcpy(ptrOutput2, "\nR = ");
      ptrOutput2 += strlen(ptrOutput2);
      showPolynomial(&ptrOutput2, ptrR, degreeR + 1, 0);
      *ptrOutput2++ = '\n';
    }
#endif
    do
    {
      IntArray2BigInteger(ptrR + degreeR * (NumberLength+1), &operand1);
      BigIntRemainder(&operand1, &primeMod, &operand2);
      if (!BigIntIsZero(&operand2))
      {
        break;
      }
    } while (--degreeR != 0);
    if (degreeR == 0)
    {
      break;
    }
    DividePolynomial(ptrOldR, degreeOldR, ptrR, degreeR, ptrQuotient);
    // ptrQuotient points to quotient of division.
    // ptrOldR points to the remainder.
    // Save pointer to the quotient.
#if DEBUG_HENSEL_LIFTING
    if (ptrOutput2 != NULL)
    {
      strcpy(ptrOutput2, "Quotient = ");
      ptrOutput2 += strlen(ptrOutput2);
      showPolynomial(&ptrOutput2, ptrQuotient, degreeOldR - degreeR + 1, 0);
      *ptrOutput2++ = '\n';
    }
#endif
    ptrQuotients[nbrQuotients++] = ptrQuotient;
    // Move pointer to new quotient after this quotient.
    ptrQuotient += (degreeOldR - degreeR + 1)*nbrLimbs;
    // Exchange pointers.
    tmpPtr = ptrR;
    ptrR = ptrOldR;
    ptrOldR = tmpPtr;
    // Get degrees of dividend and divisor of next division.
    degreeOldR = degreeR;
    degreeR = getDegreePoly(ptrR, degreeR-1);
  }
  // Save pointer to last quotient.
  ptrQuotients[nbrQuotients] = ptrQuotient;

  // R <- inverse of gcd, oldR <- 0
  *ptrP1 = NumberLength;
  memcpy(ptrP1+1, MontgomeryMultR1, NumberLength * sizeof(int));
  *ptrP2 = 1;
  *(ptrP2 + 1) = 0; 
  ptrR = ComputeUorV(ptrQuotients, nbrQuotients, ptrP1, ptrP2, &degreeFirst);
  memcpy((polyExchanged? ptrU: ptrV), ptrR, (degreeFirst+1)*nbrLimbs*sizeof(int));
  // R <- 0, oldR <- inverse of gcd.
  *ptrP2 = NumberLength;
  memcpy(ptrP2+1, MontgomeryMultR1, NumberLength * sizeof(int));
  *ptrP1 = 1;
  *(ptrP1 + 1) = 0;
  ptrR = ComputeUorV(ptrQuotients, nbrQuotients, ptrP1, ptrP2, &degreeSecond);
  memcpy((polyExchanged? ptrV : ptrU), ptrR, (degreeSecond + 1)*nbrLimbs*sizeof(int));
  degreeV = (polyExchanged ? degreeSecond : degreeFirst);
  degreeU = (polyExchanged ? degreeFirst : degreeSecond);
  *pDegreeU = degreeU;
  *pDegreeV = degreeV;
  if (polyExchanged)
  {     // If exchanged polynomials, restore pointer to polynomials
    tmpPtr = ptrA;
    ptrA = ptrB;
    ptrB = tmpPtr;
  }

  // Compute A*U+B*V using only the constant terms.
  LenAndLimbs2ArrLimbs(ptrA, operand2.limbs, nbrLimbs);
  LenAndLimbs2ArrLimbs(ptrU, operand3.limbs, nbrLimbs);
  modmult(operand2.limbs, operand3.limbs, operand4.limbs);
  LenAndLimbs2ArrLimbs(ptrB, operand2.limbs, nbrLimbs);
  LenAndLimbs2ArrLimbs(ptrV, operand3.limbs, nbrLimbs);
  modmult(operand2.limbs, operand3.limbs, operand3.limbs);
  AddBigNbrMod(operand4.limbs, operand3.limbs, operand4.limbs);

  // Find the inverse of A*U+B*V
  // Invert gcd mod prime using Montgomery notation.
  ModInvBigNbr(operand4.limbs, operand4.limbs, powerMod.limbs, powerMod.nbrLimbs);
  NumberLength = nbrLimbs - 1;

  // Multiply all coefficients of polynomial U by this inverse.
  ptrCoeff = ptrU;
  for (currentDegree = 0; currentDegree <= degreeU; currentDegree++)
  {
    LenAndLimbs2ArrLimbs(ptrCoeff, operand3.limbs, nbrLimbs);
    modmult(operand3.limbs, operand4.limbs, operand3.limbs);
    ArrLimbs2LenAndLimbs(ptrCoeff, operand3.limbs, nbrLimbs);
    ptrCoeff += nbrLimbs;
  }

  // Multiply all coefficients of polynomial V by this inverse.
  ptrCoeff = ptrV;
  for (currentDegree = 0; currentDegree <= degreeV; currentDegree++)
  {
    LenAndLimbs2ArrLimbs(ptrCoeff, operand3.limbs, nbrLimbs);
    modmult(operand3.limbs, operand4.limbs, operand3.limbs);
    ArrLimbs2LenAndLimbs(ptrCoeff, operand3.limbs, nbrLimbs);
    ptrCoeff += nbrLimbs;
  }
}

static void WidenCoefficients(int *poly, int qtyNbrs, int oldNbrLen, int newNbrLen)
{
  int currentNbr;

  for (currentNbr = qtyNbrs - 1; currentNbr >= 0; currentNbr--)
  {
    memmove(poly + currentNbr*newNbrLen, poly + currentNbr*oldNbrLen, oldNbrLen*sizeof(int));
  }
}

// Lift factorization of coefficients of factors from mod prime to
// mod power = prime^exp.
// Perform lifting mod prime, prime^2, prime^4, etc.
// All calculations are done in Montgomery domain, so the roots must 
// be converted to standard notation and then to Montgomery notation
// with greater modulus.

// Hensel lift from mod m to mod m^2:
// f = g*h (mod m), s*g + t*h = 1 (mod m)
// f = g'*h' (mod m^2), s'*g' + t'*h' = 1 (mod m^2)
// deg(s') < deg(h'), deg(t') < deg(g')
// All calculations below are performed mod m^2.
// e <- f - g*h
// Compute q, r such that s*e = q*h + r
// g' <- g + t*e + q*g
// h' <- h + r
//
// e <- s*g' + t*h' - 1
// Compute q, r such that s*e = q*h + r
// t' <- t - t*e - q*g'
// s' <- s - r
//
// The last four lines do not need to be executed on the last loop
// because the output is g' and h'.

// Memory usage:
// poly1: s', g'
// poly2: t', h'
// poly3: e
// poly4: q, r, f
// poly5: s*e, t*e, t*e + q*g'

static void Adjust(int *polyFirst, int degreeFirst,
  int *polySecond, int degreeSecond,
  int *polyH, int degreeH,
  int degreeS,
  int degreeT,
  int *polyG, int degreeG, int degreeE,
  enum eAdjustType adjustType)
{
  int nbrLimbs = NumberLength + 1;
  int degreeQ = degreeS + degreeE - degreeH;
  int *polyR;
  int currentDegree;
  int* ptrCoeff;
  if (degreeE == 0 && poly3[0] == 1 && poly3[1] == 0)
  {
    return;     // e = 0 so there is nothing to adjust.
  }
  // Compute q, r such that s*e = q*h + r
  if (degreeQ < 0)
  {
    degreeQ = 0;
  }
  polyR = &poly4[(degreeQ + 1)*nbrLimbs];
  // Compute q, r such that s*e = q*h + r
  MultPolynomial(degreeS, degreeE, polyS, poly3);   // polyMultTemp <- s*e
#if DEBUG_HENSEL_LIFTING
  if (ptrOutput2 != NULL)
  {
    strcpy(ptrOutput2, "s*e (dividend) = ");
    ptrOutput2 += strlen(ptrOutput2);
    showPolynomial(&ptrOutput2, polyMultTemp, degreeS + degreeE + 1, 0);
    strcpy(ptrOutput2, "\nh (divisor) = ");
    ptrOutput2 += strlen(ptrOutput2);
    showPolynomial(&ptrOutput2, polyH, degreeH + 1, 0);
  }
#endif
  DividePolynomial(polyMultTemp, degreeS + degreeE, polyH, degreeH, poly4);
  // At this moment poly4 = q (quotient) and polyMultTemp = r (remainder).
  memcpy(polyR, polyMultTemp, (degreeH+1)*nbrLimbs*sizeof(int));
#if DEBUG_HENSEL_LIFTING
  if (ptrOutput2 != NULL)
  {
    strcpy(ptrOutput2, "\nquotient = ");
    ptrOutput2 += strlen(ptrOutput2);
    showPolynomial(&ptrOutput2, poly4, degreeQ + 1, 0);
    strcpy(ptrOutput2, "\nremainder = ");
    ptrOutput2 += strlen(ptrOutput2);
    showPolynomial(&ptrOutput2, polyR, degreeH + 1, 0);
    *ptrOutput2++ = '\n';
  }
#endif
  // Compute t*e + q*g and then add or subtract it from polyFirst
  // according to the last parameter in this function.
  MultPolynomial(degreeT, degreeE, polyT, poly3);   // polyMultTemp <- t*e
  memcpy(poly5, polyMultTemp, (degreeT + degreeE + 1)*nbrLimbs*sizeof(int));
  // At this moment poly5 = t*e
  MultPolynomial(degreeQ, degreeG, poly4, polyG);   // polyMultTemp <- q*g
  ptrCoeff = polyFirst;       // Point to constant coefficient.
  for (currentDegree = 0; currentDegree <= degreeFirst; currentDegree++)
  {                           // For each coefficient...
    IntArray2BigInteger(&poly5[currentDegree*nbrLimbs], &operand1);
    IntArray2BigInteger(&polyMultTemp[currentDegree*nbrLimbs], &operand2);
    AddBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);  // t*e + q*g
    IntArray2BigInteger(ptrCoeff, &operand2);
    if (adjustType == ADJUST_PERFORM_SUBTRACTION)
    {                         // t' <- t - t*e - q*g'
      SubtBigNbrMod(operand2.limbs, operand1.limbs, operand1.limbs);
    }
    else
    {                         // g' <- g + t*e + q*g
      AddBigNbrMod(operand2.limbs, operand1.limbs, operand1.limbs);
    }
    BigInteger2IntArray(ptrCoeff, &operand1);
    ptrCoeff += nbrLimbs;     // Point to next coefficient.
  }
  ptrCoeff = polySecond;      // Point to constant coefficient.
  for (currentDegree = 0; currentDegree <= degreeH && currentDegree <= degreeSecond; currentDegree++)
  {                           // For each coefficient...
    IntArray2BigInteger(ptrCoeff, &operand2);
    IntArray2BigInteger(polyR+currentDegree*nbrLimbs, &operand1);
    if (adjustType == ADJUST_PERFORM_SUBTRACTION)
    {                         // s' <- s - r
      SubtBigNbrMod(operand2.limbs, operand1.limbs, operand1.limbs);
    }
    else
    {                         // h' <- h + r
      AddBigNbrMod(operand2.limbs, operand1.limbs, operand1.limbs);
    }
    BigInteger2IntArray(ptrCoeff, &operand1);
    ptrCoeff += nbrLimbs;     // Point to next coefficient.
  }
}

// Perform Hensel lifting of polynomials in factorInfo.
// On entry:
// values = original polynomial. The product of all lifted factors
// will be equal to this polynomial.
// exponentMod = exponent to perform the lift.
// primeMod = prime as a big integer.
// nbrFactorsFound = number of factors of original polynomial.
// factorInfo = pointer to array of structures holding the polynomial factors.
// Hensel lifting only works when there are no repeated factors mod p,
// so this condition is tested.

int HenselLifting(struct sFactorInfo* factorInfo)
{
  int nbrFactor, nbrFactor2, currentExp, oldNumberLength;
  int degreeFactor, degreeS, degreeT, degreeG, currentDegree;
  int index, nbrLimbs, *ptrValue1;
  int* ptrPolyLifted = polyLifted;
  int* ptrPolyLiftedBak = polyLiftedBak;
  int* ptrPoly2;
  struct sFactorInfo *pstFactorInfo, *pstFactorInfo2;

  if (exponentMod == 1)
  {      // No lift to be done. Copy polynomials to polyLifted.
    int* ptrDest = ptrPolyLifted;
    pstFactorInfo = factorInfo;
    for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
    {
      pstFactorInfo->ptrPolyLifted = ptrDest;
      ptrDest = CopyPolynomial(ptrDest, pstFactorInfo->ptr, pstFactorInfo->degree);
      pstFactorInfo++;
    }
    return EXPR_OK;   // No lift has to be done.
  }
  CopyBigInt(&powerMod, &primeMod);
  // Hensel lifting only works when there are no repeated factors mod p.
  // Go out if the multiplicity of a factor is not 1.
  // Compute in variable degree the degree of the product.
  degree = 0;
  ptrPoly2 = poly2;
  pstFactorInfo = factorInfo;
  for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    if (pstFactorInfo->multiplicity != 1)
    {   // Repeated factor found. Indicate error.
      return EXPR_CANNOT_LIFT;
    }
    // Init poly2 coefficients to zero.
    for (index = 0; index < pstFactorInfo->degree; index++)
    {
      *ptrPoly2++ = 1;
      *ptrPoly2 = 0;
      ptrPoly2 += NumberLength;
    }
    degree += pstFactorInfo->degree;
    pstFactorInfo++;
  }
  // At this moment, variable degree holds the degree of the original polynomial.
  pstFactorInfo = factorInfo;
  for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
#if DEBUG_HENSEL_LIFTING
    if (ptrOutput2 != NULL)
    {
      strcpy(ptrOutput2, "---------------------------------------\n");
      ptrOutput2 += strlen(ptrOutput2);
      sprintf(ptrOutput2, "nbrFactor = %d\n", nbrFactor);
      ptrOutput2 += strlen(ptrOutput2);
    }
#endif
    computePower(1);
    nbrLimbs = NumberLength + 1;
    SetNumberToOne(poly1);
    degreeFactor = 0;
    pstFactorInfo2 = factorInfo;
    for (nbrFactor2 = 0; nbrFactor2 < nbrFactorsFound; nbrFactor2++)
    {    // Loop that multiplies all factors except the current one. poly1 holds the result.
         // Coefficients are in Montgomery format.
      if (nbrFactor != nbrFactor2)
      {
        memcpy(poly2, pstFactorInfo2->ptr, pstFactorInfo2->degree * nbrLimbs * sizeof(int));
        SetNumberToOne(&poly2[pstFactorInfo2->degree*nbrLimbs]);
        MultPolynomial(degreeFactor, pstFactorInfo2->degree, poly1, poly2);
#if DEBUG_HENSEL_LIFTING
        if (ptrOutput2 != NULL)
        {
          strcpy(ptrOutput2, "poly1 = ");
          ptrOutput2 += strlen(ptrOutput2);
          showPolynomial(&ptrOutput2, poly1, degreeFactor + 1, 0);
          strcpy(ptrOutput2, "\npoly2 = ");
          ptrOutput2 += strlen(ptrOutput2);
          showPolynomial(&ptrOutput2, poly2, pstFactorInfo2->degree + 1, 0);
          strcpy(ptrOutput2, "\npolyMultTemp = ");
          ptrOutput2 += strlen(ptrOutput2);
          showPolynomial(&ptrOutput2, polyMultTemp, degreeFactor + pstFactorInfo2->degree + 1, 0);
          *ptrOutput2++ = '\n';
        }
#endif
        degreeFactor += pstFactorInfo2->degree;
        memcpy(poly1, polyMultTemp, (degreeFactor + 1)*nbrLimbs*sizeof(int));
      }
      pstFactorInfo2++;
    }
    memcpy(poly2, pstFactorInfo->ptr, pstFactorInfo->degree * nbrLimbs * sizeof(int));
    degreeFactor = pstFactorInfo->degree;
    SetNumberToOne(&poly2[degreeFactor*nbrLimbs]);
    // At this moment g = poly1, h = poly2.
    degreeG = degree - degreeFactor;

    // From coprime polynomials poly1 and poly2 of degree degreeG and degreeFactor
    // respectively, find polynomials poly3, poly4, poly5 and polyS such that
    // poly1*polyS + poly2*polyT = 1 (mod powerMod) where polyS has degree degreeS
    // and polyT has degree degreeT. Use Montgomery notation.
    ExtendedGcdPolynomial(poly1, degreeG, poly2, degreeFactor, poly3,
      poly4, poly5, polyS, &degreeS, polyT, &degreeT);
    // At this moment, s = polyS, t = polyT.
    // Complete leading coefficients of s and t such that:
    // degree(s) = degree(h) - 1, degree(t) = degree(g) - 1.
    for (currentDegree = degreeS+1; currentDegree < degreeFactor; currentDegree++)
    {
      polyS[currentDegree*nbrLimbs] = 1;   // Set coefficient to zero.
      polyS[currentDegree*nbrLimbs+1] = 0;
    }
    for (currentDegree = degreeT+1; currentDegree < degreeG; currentDegree++)
    {
      polyT[currentDegree*nbrLimbs] = 1;   // Set coefficient to zero.
      polyT[currentDegree*nbrLimbs + 1] = 0;
    }
#if DEBUG_HENSEL_LIFTING
    if (ptrOutput2 != NULL)
    {
      strcpy(ptrOutput2, "poly1 = ");
      ptrOutput2 += strlen(ptrOutput2);
      showPolynomial(&ptrOutput2, poly1, degreeG+1, 0);
      strcpy(ptrOutput2, "\npoly2 = ");
      ptrOutput2 += strlen(ptrOutput2);
      showPolynomial(&ptrOutput2, poly2, degreeFactor+1, 0);
      strcpy(ptrOutput2, "\npolyS = ");
      ptrOutput2 += strlen(ptrOutput2);
      showPolynomial(&ptrOutput2, polyS, degreeS+1, 0);
      strcpy(ptrOutput2, "\npolyT = ");
      ptrOutput2 += strlen(ptrOutput2);
      showPolynomial(&ptrOutput2, polyT, degreeT+1, 0);
      *ptrOutput2++ = '\n';
    }
#endif
    degreeS = degreeFactor - 1;
    degreeT = degreeG - 1;
    if (degreeT < 0)
    {
      degreeT = 0;
    }
    currentExp = 1;
    pstFactorInfo++;
    // Loop that performs the lifting.
    while (currentExp != exponentMod)
    {
#if DEBUG_HENSEL_LIFTING
      if (ptrOutput2 != NULL)
      {
        sprintf(ptrOutput2, "currentExp = %d, ptrPolyLiftedBak = ", currentExp);
        ptrOutput2 += strlen(ptrOutput2);
        showPolynomial(&ptrOutput2, poly2, degreeFactor + 1, 0);
        *ptrOutput2++ = '\n';
      }
#endif
#ifdef __EMSCRIPTEN__
      int elapsedTime = (int)(tenths() - originalTenthSecond);
      if (elapsedTime / 10 != oldTimeElapsed / 10)
      {
        char outputInfo[1000];
        char *ptrOutput = outputInfo;
        if (lang)
        {
          strcpy(ptrOutput, "1<p>Aplicando lema de Hensel en el factor número ");
          ptrOutput += strlen(ptrOutput);
          int2dec(&ptrOutput, nbrFactor + 1);
          strcpy(ptrOutput, " de ");
          ptrOutput += strlen(ptrOutput);
          int2dec(&ptrOutput, nbrFactorsFound);
          strcpy(ptrOutput, " usando el número primo ");
          ptrOutput += strlen(ptrOutput);
          int2dec(&ptrOutput, primeMod.limbs[0].x);
          strcpy(ptrOutput, " procesando exponente ");
          ptrOutput += strlen(ptrOutput);
          int2dec(&ptrOutput, currentExp);
          strcpy(ptrOutput, " de ");
        }
        else
        {
          strcpy(ptrOutput, "1<p>Hensel lifting of factor number ");
          ptrOutput += strlen(ptrOutput);
          int2dec(&ptrOutput, nbrFactor + 1);
          strcpy(ptrOutput, " of ");
          ptrOutput += strlen(ptrOutput);
          int2dec(&ptrOutput, nbrFactorsFound);
          strcpy(ptrOutput, " using prime number ");
          ptrOutput += strlen(ptrOutput);
          int2dec(&ptrOutput, primeMod.limbs[0].x);
          strcpy(ptrOutput, " processing exponent ");
          ptrOutput += strlen(ptrOutput);
          int2dec(&ptrOutput, currentExp);
          strcpy(ptrOutput, " of ");
        }
        ptrOutput += strlen(ptrOutput);
        int2dec(&ptrOutput, exponentMod);
        strcpy(ptrOutput, ".</p>");
        ptrOutput += strlen(ptrOutput);
        showElapsedTime(&ptrOutput);
        databack(outputInfo);
      }
#endif
#if DEBUG_HENSEL_LIFTING
      if (ptrOutput2 != NULL)
      {
        strcpy(ptrOutput2, "\npolyS = ");
        ptrOutput2 += strlen(ptrOutput2);
        showPolynomial(&ptrOutput2, polyS, degreeS + 1, 0);
        strcpy(ptrOutput2, "\npolyT = ");
        ptrOutput2 += strlen(ptrOutput2);
        showPolynomial(&ptrOutput2, polyT, degreeT + 1, 0);
        *ptrOutput2++ = '\n';
      }
#endif
      // Convert g, h, s, t from Montgomery notation to standard notation
      // by multiplying by 1 in Montgomery notation.
      polyToStandardNotation(poly1, degreeG + 1);      // Convert g to standard notation.
      polyToStandardNotation(poly2, degreeFactor + 1); // Convert h to standard notation.
      polyToStandardNotation(polyS, degreeS + 1);      // Convert s to standard notation.
      polyToStandardNotation(polyT, degreeT + 1);      // Convert t to standard notation.
#if DEBUG_HENSEL_LIFTING
      if (ptrOutput2 != NULL)
      {
        strcpy(ptrOutput2, "Standard notation:\npoly1 = ");
        ptrOutput2 += strlen(ptrOutput2);
        showPolynomial(&ptrOutput2, poly1, degreeG + 1, 0);
        *ptrOutput2++ = '\n';
      }
#endif
      oldNumberLength = NumberLength;
      currentExp <<= 1;                                // We can double the exponent in each step.
      if (currentExp > exponentMod)                    // Do not exceed exponentMod.
      {
        currentExp = exponentMod;
      }
      computePower(currentExp);                        // Compute powerMod and init Montgomery parms.
      // Compute polynomial poly4 <- value mod prime^currentExp.
      degree = values[0];                              // Get degree of polynomial
      ptrValue1 = &values[1];                          // Point to constant coefficient.
      for (currentDegree = 0; currentDegree <= degree; currentDegree++)
      {
        NumberLength = numLimbs(ptrValue1);
        IntArray2BigInteger(ptrValue1, &operand1);
        (void)BigIntRemainder(&operand1, &powerMod, &operand1);
        NumberLength = powerMod.nbrLimbs;
        if (operand1.nbrLimbs < NumberLength)
        {
          memset(&operand1.limbs[operand1.nbrLimbs], 0,
            (NumberLength - operand1.nbrLimbs) * sizeof(limb));
        }
        // Convert operand1 from standard to Montgomery notation.
        modmult(operand1.limbs, MontgomeryMultR2, operand1.limbs);
        BigInteger2IntArray(&poly4[currentDegree*(NumberLength + 1)], &operand1);
        ptrValue1 += 1 + numLimbs(ptrValue1);           // Point to next coefficient.
      }
      // Convert polynomial mod prime to monic (leading coefficient must be 1).
      ConvertToMonic(poly4, degree);
      if (oldNumberLength != NumberLength)
      {
        // Now the coefficients require more limbs. Widen them so they can fit.
        // Old size = oldNumberLength, new size = NumberLength.
        WidenCoefficients(poly1, degreeG + 1, oldNumberLength + 1, NumberLength + 1);
        WidenCoefficients(poly2, degreeFactor + 1, oldNumberLength + 1, NumberLength + 1);
        WidenCoefficients(polyS, degreeS + 1, oldNumberLength + 1, NumberLength + 1);
        WidenCoefficients(polyT, degreeT + 1, oldNumberLength + 1, NumberLength + 1);
        nbrLimbs = NumberLength + 1;
      }
        // Convert standard notation to Montgomery notation using new modulus.
      polyToMontgomeryNotation(poly1, degreeG + 1);                // Convert g to Montgomery notation.
      polyToMontgomeryNotation(poly2, degreeFactor + 1);           // Convert h to Montgomery notation.
      polyToMontgomeryNotation(polyS, degreeS + 1);                // Convert s to Montgomery notation.
      polyToMontgomeryNotation(polyT, degreeT + 1);                // Convert t to Montgomery notation.
      // Compute e <- f - g*h and store it into poly3.
      MultPolynomial(degree - degreeFactor, degreeFactor, poly1, poly2);
      for (index = 0; index <= degree; index++)
      {                                                  // For each coefficient...
        IntArray2BigInteger(&poly4[index*nbrLimbs], &operand1);         // f
        IntArray2BigInteger(&polyMultTemp[index*nbrLimbs], &operand2);  // g*h
        SubtBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);  // f - g*h
        BigInteger2IntArray(&poly3[index*nbrLimbs], &operand1);         // Store e <- f - g*h.
      }
#if DEBUG_HENSEL_LIFTING
      if (ptrOutput2 != NULL)
      {
        strcpy(ptrOutput2, "After widening:\npoly1 = ");
        ptrOutput2 += strlen(ptrOutput2);
        showPolynomial(&ptrOutput2, poly1, degreeG + 1, 0);
        strcpy(ptrOutput2, "\npoly2 = ");
        ptrOutput2 += strlen(ptrOutput2);
        showPolynomial(&ptrOutput2, poly2, degreeFactor + 1, 0);
        strcpy(ptrOutput2, "\npoly3 = ");
        ptrOutput2 += strlen(ptrOutput2);
        showPolynomial(&ptrOutput2, poly3, degree + 1, 0);
        strcpy(ptrOutput2, "\npoly4 = ");
        ptrOutput2 += strlen(ptrOutput2);
        showPolynomial(&ptrOutput2, poly4, degree + 1, 0);
        strcpy(ptrOutput2, "\npolyS = ");
        ptrOutput2 += strlen(ptrOutput2);
        showPolynomial(&ptrOutput2, polyS, degreeS + 1, 0);
        strcpy(ptrOutput2, "\npolyT = ");
        ptrOutput2 += strlen(ptrOutput2);
        showPolynomial(&ptrOutput2, polyT, degreeT + 1, 0);
        *ptrOutput2++ = '\n';
      }
#endif
      // Compute q, r such that s*e = q*h + r
      // g' <- g + t*e + q*g
      // h' <- h + r
      Adjust(poly1, degreeG,                             // g
        poly2, degreeFactor,                             // h
        poly2, degreeFactor,                             // h
        degreeS,                                         // s
        degreeT,                                         // t
        poly1, degreeG,                                  // g
        getDegreePoly(poly3, degree),                    // degree of e
        ADJUST_PERFORM_ADDITION);                        // Perform addition.
      if (currentExp == exponentMod)
      {                 // Final values of g' and h' found. Exit loop.
        break;
      }
      // Compute e <- s*g' + t*h' - 1 and store it into poly3.
      MultPolynomial(degreeS, degreeG, polyS, poly1);        // Compute s*g'
      memcpy(poly3, polyMultTemp, (degreeS + degreeG + 1)*nbrLimbs*sizeof(int));
      MultPolynomial(degreeT, degreeFactor, polyT, poly2);   // Compute t*h'
      for (currentDegree = 0; currentDegree <= degreeS+degreeG; currentDegree++)
      {                                                      // Loop that computes s*g' + t*h'
        IntArray2BigInteger(&poly3[currentDegree*nbrLimbs], &operand1);
        IntArray2BigInteger(&polyMultTemp[currentDegree*nbrLimbs], &operand2);
        AddBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
        BigInteger2IntArray(&poly3[currentDegree*nbrLimbs], &operand1);
      }
      IntArray2BigInteger(poly3, &operand1);                             // Get constant coefficient.
      SubtBigNbrMod(operand1.limbs, MontgomeryMultR1, operand1.limbs);   // Subtract 1 in Montgomery notation.
      BigInteger2IntArray(poly3, &operand1);                             // Store e.
      // Compute q, r such that s*e = q*h + r
      // t' <- t - t*e - q*g'
      // s' <- s - r
      Adjust(polyT, degreeT,             // t
        polyS, degreeS,                  // s
        poly2, degreeFactor,             // h
        degreeS,                         // s
        degreeT,                         // t
        poly1, degreeG,                  // g
        getDegreePoly(poly3, degreeS+degreeG),  // degree of e
        ADJUST_PERFORM_SUBTRACTION);     // Perform subtraction.
    }
    // Copy lifted factor to temporary array of polynomials.
    memcpy(ptrPolyLiftedBak, poly2, degreeFactor*nbrLimbs*sizeof(int));
#if DEBUG_HENSEL_LIFTING
    if (ptrOutput2 != NULL)
    {
      sprintf(ptrOutput2, "currentExp = %d, ptrPolyLiftedBak = ", currentExp);
      ptrOutput2 += strlen(ptrOutput2);
      showPolynomial(&ptrOutput2, poly2, degreeFactor + 1, 0);
      *ptrOutput2++ = '\n';
    }
#endif
    // Point to next factor polynomial.
    ptrPolyLifted += degreeFactor * (oldNumberLength+1);
    // Point to next lifted factor polynomial (to fill in next loop).
    ptrPolyLiftedBak += degreeFactor * nbrLimbs;
  }
  // Convert factors to standard notation.
  pstFactorInfo = factorInfo;
  ptrPolyLifted = polyLifted;
  ptrPolyLiftedBak = polyLiftedBak;
  for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {              // For each factor...
    int polyLength;
    degreeFactor = pstFactorInfo->degree;
    polyLength = degreeFactor * nbrLimbs;
    pstFactorInfo->ptrPolyLifted = ptrPolyLifted;
    memcpy(ptrPolyLifted, ptrPolyLiftedBak, polyLength*sizeof(int));
    ptrPolyLifted += polyLength;
    ptrPolyLiftedBak += polyLength;
    polyToStandardNotation(pstFactorInfo->ptrPolyLifted, degreeFactor);
    pstFactorInfo++;
  }
  return EXPR_OK;
}

void OrigPolyFromMontgomeryToStandard(void)
{
  int *ptrValue1, *ptrValue2;
  int currentDegree;
  memcpy(&TestNbr, powerMod.limbs, powerMod.nbrLimbs*sizeof(limb));
  NumberLength = powerMod.nbrLimbs;
  TestNbr[NumberLength].x = 0;
  GetMontgomeryParms(powerMod.nbrLimbs);
  degree = values[0];
  ptrValue1 = &values[1];
  ptrValue2 = &poly4[0];
  for (currentDegree = 0; currentDegree <= degree; currentDegree++)
  {
    IntArray2BigInteger(ptrValue1, &operand1);
    // Convert operand1 from Montgomery to standard notation.
    operand2.limbs[0].x = 1;
    if (NumberLength > 1)
    {
      memset(&operand2.limbs[1], 0, (NumberLength - 1)*sizeof(limb));
    }
    modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    BigInteger2IntArray(ptrValue2, &operand1);
    ptrValue1 += 1 + *ptrValue1;
    ptrValue2 += 1 + *ptrValue2;
  }
  memcpy(&values[1], &poly4[0], (ptrValue2 - &poly4[0])*sizeof(int));
}

static void showPower(char **pptrOutput, int exponent)
{
  char *ptrOutput = *pptrOutput;
  if (pretty)
  {
    *ptrOutput++ = '<';
    *ptrOutput++ = 's';
    *ptrOutput++ = 'u';
    *ptrOutput++ = 'p';
    *ptrOutput++ = '>';
    int2dec(&ptrOutput, exponent);
    *ptrOutput++ = '<';
    *ptrOutput++ = '/';
    *ptrOutput++ = 's';
    *ptrOutput++ = 'u';
    *ptrOutput++ = 'p';
    *ptrOutput++ = '>';
  }
  else
  {
    *ptrOutput++ = '^';
    int2dec(&ptrOutput, exponent);
  }
  *pptrOutput = ptrOutput;
}

void showPowerX(char **pptrOutput, int polyDegree)
{
  char *ptrOutput = *pptrOutput;
  if (polyDegree == 0)
  {
    *ptrOutput++ = '1';
  }
  else
  {
    if (pretty)
    {
      strcpy(ptrOutput, "<var>");
      ptrOutput += strlen(ptrOutput);
    }
    *ptrOutput++ = 'x';
    if (pretty)
    {
      strcpy(ptrOutput, "</var>");
      ptrOutput += strlen(ptrOutput);
    }
    if (polyDegree != 1)
    {
      showPower(&ptrOutput, polyDegree);
    }
  }
  *pptrOutput = ptrOutput;
}

static void showPolynomial(char **pptrOutput, int *ptrPoly, int polyDegree, int groupLength)
{
  int nbrLimbs = NumberLength + 1;
  int currentDegree;
  char *ptrOutput = *pptrOutput;
  int *ptrIndex;
  int indexes[MAX_DEGREE];
  int NumberLengthBak = NumberLength;

  ptrIndex = &indexes[0];
  if (modulusIsZero)
  {
    // Fill indexes to start of each coefficient.
    int index = 0;
    for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
    {
      *ptrIndex++ = index;
      index += numLimbs(ptrPoly + index) + 1;
    }
  }
  else
  {
    NumberLength = powerMod.nbrLimbs;
    nbrLimbs = NumberLength + 1;
    for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
    {
      *ptrIndex++ = currentDegree * nbrLimbs;
    }
  }
  for (currentDegree = polyDegree - 1; currentDegree >= 0; currentDegree--)
  {
    int* ptrValue1 = ptrPoly + indexes[currentDegree];
    int len = numLimbs(ptrValue1);
    if (len != 1 || *(ptrValue1 + 1) != 0)
    {            // Coefficient is not zero.
      *ptrOutput++ = ' ';
      if (*ptrValue1 > 0)
      {
        *ptrOutput++ = '+';
      }
      else
      {
        strcpy(ptrOutput, pretty? "&minus;": "-");
        ptrOutput += strlen(ptrOutput);
      }
      *ptrOutput++ = ' ';
      if (len != 1 || *(ptrValue1 + 1) != 1)
      {            // Absolute value of coefficient is not one.
        NumberLength = numLimbs(ptrValue1);
        IntArray2BigInteger(ptrValue1, &operand1);
        operand1.sign = SIGN_POSITIVE;
        Bin2Dec(operand1.limbs, ptrOutput, operand1.nbrLimbs, groupLength);
        ptrOutput += strlen(ptrOutput);
        if (currentDegree > 0)
        {
          if (pretty)
          {
            strcpy(ptrOutput, "&#8290;");
            ptrOutput += strlen(ptrOutput);
          }
          else
          {
            *ptrOutput++ = '*';
          }
          showPowerX(&ptrOutput, currentDegree);
        }
      }
      else
      {
        showPowerX(&ptrOutput, currentDegree);
      }
    }
    ptrValue1 -= nbrLimbs;
  }
  *pptrOutput = ptrOutput;
  NumberLength = NumberLengthBak;
}

void outputOriginalPolynomial(char* ptrOutput, int groupLength)
{
  int currentDegree;
  int* ptrValue1;
  int nbrLimbs = powerMod.nbrLimbs + 1;
  degree = values[0];
  ptrValue1 = &values[1];
  if (!modulusIsZero)
  {
    int* ptrValue2;
    // Output polynomial to factor. First move polynomial to poly4
    ptrValue2 = &poly4[0];
    for (currentDegree = 0; currentDegree <= degree; currentDegree++)
    {
      memcpy(ptrValue2, ptrValue1, (1 + *ptrValue1) * sizeof(int));
      ptrValue1 += 1 + *ptrValue1;
      ptrValue2 += nbrLimbs;
    }
    // Get leading coefficient.
    IntArray2BigInteger(&poly4[degree * nbrLimbs], &operand1);
  }
  else
  { // Find leading coefficient.
    for (currentDegree = 0; currentDegree < degree; currentDegree++)
    {
      ptrValue1 += 1 + numLimbs(ptrValue1);
    }
    IntArray2BigInteger(ptrValue1, &operand1);
  }
  if (operand1.sign == SIGN_NEGATIVE)
  {
    strcpy(ptrOutput, " &minus;");
    ptrOutput += strlen(ptrOutput);
  }
  if ((operand1.nbrLimbs != 1 || operand1.limbs[0].x != 1) || degree == 0)
  {     // Leading coefficient is not 1 or degree is zero.
    Bin2Dec(operand1.limbs, ptrOutput, operand1.nbrLimbs, groupLength);
    ptrOutput += strlen(ptrOutput);
  }
  if (degree > 0)
  {
    showPowerX(&ptrOutput, degree);
  }
  showPolynomial(&ptrOutput, (modulusIsZero ? &values[1] : poly4), degree, groupLength);
  if (!modulusIsZero)
  {
    strcpy(ptrOutput, " (mod ");
    ptrOutput += strlen(ptrOutput);
    Bin2Dec(primeMod.limbs, ptrOutput, primeMod.nbrLimbs, groupLength);
    ptrOutput += strlen(ptrOutput);
    if (exponentMod != 1)
    {
      showPower(&ptrOutput, exponentMod);
    }
    *ptrOutput++ = ')';
  }
  *ptrOutput = 0;    // Append string terminator.
}

void outputPolynomialFactor(char *ptrOutput, int groupLength, struct sFactorInfo* pstFactorInfo)
{
  int polyDegree = pstFactorInfo->degree;
  int multiplicity = pstFactorInfo->multiplicity;
  int isMonomial = (polyDegree == 1 && *pstFactorInfo->ptrPolyLifted == 1 &&
     *(pstFactorInfo->ptrPolyLifted+1) == 0);
  if (multiplicity > 1 && !isMonomial)
  {
    *ptrOutput++ = '(';
  }
  if (modulusIsZero)
  {
    int currentDegree;
    // Get leading coefficient.
    int *ptrSrc = pstFactorInfo->ptrPolyLifted;
    for (currentDegree = 0; currentDegree < polyDegree; currentDegree++)
    {
      ptrSrc += 1 + numLimbs(ptrSrc);
    }
    UncompressBigIntegerB(ptrSrc, &operand1);
    if (operand1.sign == SIGN_NEGATIVE)
    {
      strcpy(ptrOutput, "&minus;");
      ptrOutput += strlen(ptrOutput);
    }
    if (operand1.nbrLimbs != 1 || operand1.limbs[0].x != 1)
    {     // Absolute value is not 1.
      Bin2Dec(operand1.limbs, ptrOutput, operand1.nbrLimbs, groupLength);
      ptrOutput += strlen(ptrOutput);
    }
  }
  showPowerX(&ptrOutput, polyDegree);
  showPolynomial(&ptrOutput, pstFactorInfo->ptrPolyLifted, polyDegree, groupLength);
  if (multiplicity > 1)
  {
    if (!isMonomial)
    {
      *ptrOutput++ = ')';
    }
    showPower(&ptrOutput, multiplicity);
  }
  *ptrOutput = 0;    // Append string terminator.
}

void textErrorPol(char *ptrOutput, enum eExprErr rc)
{
  char text[150];

  switch (rc)
  {
  case EXPR_CANNOT_USE_X_IN_EXPONENT:
    strcpy(text, lang?"No se puede usar x en el exponente":
                      "Cannot use x in exponent");
    break;
  case EXPR_POLYNOMIAL_DIVISION_NOT_INTEGER:
    strcpy(text, lang ? "La división de polinomios no es entera" :
                        "Polynomial division is not integer");
    break;
  case EXPR_DEGREE_TOO_HIGH:
    strcpy(text, lang?"El grado del polinomio es muy elevado":
                      "Degree is too high");
    break;
  case EXPR_EXPONENT_TOO_LARGE:
    strcpy(text, lang?"Exponente muy grande":"Exponent is too large");
    break;
  case EXPR_EXPONENT_NEGATIVE:
    strcpy(text, lang?"Exponente negativo":"Exponent is negative");
    break;
  case EXPR_LEADING_COFF_MULTIPLE_OF_PRIME:
    strcpy(text, lang?"El primer coeficiente es múltiplo del número primo":
      "Leading coefficient multiple of prime");
    break;
  case EXPR_CANNOT_LIFT:
    strcpy(text, lang?"No se puede elevar porque hay factores duplicados":
      "Cannot lift because of duplicate factors modulo prime");
    break;
  case EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE:
    strcpy(text, lang?"El módulo debe ser mayor que 1":"Modulus must be greater than one");
    break;
  case EXPR_MODULUS_MUST_BE_PRIME_EXP:
    strcpy(text, lang ? "El módulo debe ser un número primo o una potencia de número primo" :
      "Modulus must be a prime number or a power of a prime");
    break;
  default:
    textError(text, rc);
  }
  *ptrOutput++ = '<';
  *ptrOutput++ = 'p';
  *ptrOutput++ = '>';
  strcpy(ptrOutput, text);
  ptrOutput += strlen(ptrOutput);
  *ptrOutput++ = '<';
  *ptrOutput++ = '/';
  *ptrOutput++ = 'p';
  *ptrOutput++ = '>';
  *ptrOutput = 0;    // Add terminator character.
}


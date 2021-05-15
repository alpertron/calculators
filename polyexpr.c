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

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include "bignbr.h"
#include "polynomial.h"

#define STACK_OPER_SIZE            100
#define TOKEN_NUMBER               '0'
#define TOKEN_START_EXPON          '1'
#define TOKEN_END_EXPON            '2'
#define TOKEN_UNARY_MINUS          '3'
#define TOKEN_GCD                  '4'
#define TOKEN_DER                  '5'

static BigInteger value;
static char RPNbuffer[COMPRESSED_POLY_MAX_LENGTH];
static int GcdPolynomialExpr(int* ptrArgument1, int* ptrArgument2);

static bool isFunc(char** ppcInput, const char* funcName)
{
  char* pcInput = *ppcInput;
  const char* ptrFuncName = funcName;
  while (*ptrFuncName != '\0')
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
static int ConvertToReversePolishNotation(char* input, char* ptrOut)
{
  char* ptrOutput = ptrOut;
  char* pInput = input;
  int exponOperatorCounter = 0;
  char stackOper[STACK_OPER_SIZE];
  int stackOperIndex = 0;
  bool prevTokenIsNumber = false;
  char s;
  char* ptrInput;
  int limb;
  int bitNbr;
  char variableLetter = ' ';  // Indicate variable letter not known yet.
  while (*pInput != '\0')
  {
    char* inputTemp;
    char c = *pInput;
    char cUppercase;
    if ((c == ' ') || (c == 9))
    {          // Ignore any spaces and tabs.
      pInput++;
      continue;
    }
    inputTemp = pInput;
    if (isFunc(&inputTemp, "GCD"))
    {
      if (prevTokenIsNumber == false)
      {
        pInput = inputTemp;
        stackOper[stackOperIndex] = TOKEN_GCD;  // Push token onto stack.
        stackOperIndex++;
        continue;
      }
      c = '*';
    }
    else if (isFunc(&inputTemp, "DER"))
    {
      if (prevTokenIsNumber == false)
      {
        pInput = inputTemp;
        stackOper[stackOperIndex] = TOKEN_DER;  // Push token onto stack.
        stackOperIndex++;
        continue;
      }
      c = '*';
    }
    else
    {
      pInput++;
    }
    cUppercase = c & 0xDF;
    if ((cUppercase >= 'A') && (cUppercase <= 'Z'))
    {          // Letter found.
      if ((variableLetter != cUppercase) && (variableLetter != ' '))
      {
        return EXPR_MULTIPLE_VARIABLES_NOT_ACCEPTED;
      }
      variableLetter = cUppercase;
    }
    if ((c == '*') && (*pInput == '*'))
    {          // Convert double asterisk to exponentiation.
      c = '^';
      pInput++;
    }
    if (prevTokenIsNumber)
    {
      prevTokenIsNumber = false;
      if ((c == '+') || (c == '-'))
      {        // Binary plus or minus
        while (stackOperIndex > 0)
        {      // Send operators to output.
          stackOperIndex--;
          s = stackOper[stackOperIndex];
          if ((s != '+') && (s != '-') && (s != '*') && (s != '^') &&
            (s != '/') && (s != '%') && (s != TOKEN_UNARY_MINUS))
          {    // Operator on stack has less precedence.
            stackOperIndex++;
            break;
          }
          if (s == '^')
          {
            exponOperatorCounter--;
            if (exponOperatorCounter == 0)
            {
              *ptrOutput = TOKEN_END_EXPON;
              ptrOutput++;
            }
          }
          *ptrOutput = s;
          ptrOutput++;
        }
        stackOper[stackOperIndex] = c;  // Push operator onto stack.
        stackOperIndex++;
      }
      else if ((c == '*') || (c == '(') || (c == '/') || (c == '%') ||
        (cUppercase == variableLetter))
      {
        while (stackOperIndex > 0)
        {      // Send operators to output.
          stackOperIndex--;
          s = stackOper[stackOperIndex];
          if ((s != '^') && (s != '*') && (s != '/') && (s != '%'))
          {    // Operator on stack has less precedence.
            stackOperIndex++;
            break;
          }
          if (s == '^')
          {
            exponOperatorCounter--;
            if (exponOperatorCounter == 0)
            {
              *ptrOutput = TOKEN_END_EXPON;
              ptrOutput++;
            }
          }
          *ptrOutput = s;
          ptrOutput++;
        }
        if (cUppercase == variableLetter)
        {
          if (exponOperatorCounter != 0)
          {
            return EXPR_CANNOT_USE_X_IN_EXPONENT;
          }
          stackOper[stackOperIndex] = '*';    // Push operator onto stack.
          stackOperIndex++;
          *ptrOutput = 'x';
          ptrOutput++;
          prevTokenIsNumber = true;
        }
        else
        {
          if (c == '(')
          {
            stackOper[stackOperIndex] = '*';  // Push operator onto stack.
            stackOperIndex++;
          }
          stackOper[stackOperIndex] = c;      // Push operator onto stack.
          stackOperIndex++;
        }
      }
      else if (c == '^')
      {
        while (stackOperIndex > 0)
        {      // Send operators to output.
          stackOperIndex--;
          s = stackOper[stackOperIndex];
          if (s != '^')
          {    // Operator on stack has less precedence.
            stackOperIndex++;
            break;
          }
          exponOperatorCounter--;
          if (exponOperatorCounter == 0)
          {
            *ptrOutput = TOKEN_END_EXPON;
            ptrOutput++;
          }
          *ptrOutput = s;
          ptrOutput++;
        }
        stackOper[stackOperIndex] = c;  // Push operator onto stack.
        stackOperIndex++;
        exponOperatorCounter++;
        if (exponOperatorCounter == 1)
        {
          *ptrOutput = TOKEN_START_EXPON;
          ptrOutput++;
        }
      }
      else if ((c == ')') || (c == ','))
      {
        s = '\0'; // Assume parenthesis mismatch.
        while (stackOperIndex > 0)
        {      // Send operators to output.
          stackOperIndex--;
          s = stackOper[stackOperIndex];
          if (s == '(')
          {    // Operator on stack has less precedence.
            break;
          }
          if (s == '^')
          {
            exponOperatorCounter--;
            if (exponOperatorCounter == 0)
            {
              *ptrOutput = TOKEN_END_EXPON;
              ptrOutput++;
            }
          }
          *ptrOutput = s;
          ptrOutput++;
        }
        if (s != '(')
        {     // Parentheses mismatch.
          return EXPR_PAREN_MISMATCH;
        }
        if (c == ',')
        {
          stackOper[stackOperIndex] = s;  // Push back paren.
          stackOperIndex++;
          prevTokenIsNumber = false;
        }
        else
        {
          if (stackOperIndex > 0)
          {
            if ((stackOper[stackOperIndex - 1] == TOKEN_GCD) ||
              (stackOper[stackOperIndex - 1] == TOKEN_DER))
            {
              stackOperIndex--;
              *ptrOutput = stackOper[stackOperIndex];
              ptrOutput++;
            }
          }
          prevTokenIsNumber = true;
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
        stackOper[stackOperIndex] = TOKEN_UNARY_MINUS;  // Push operator onto stack.
        stackOperIndex++;
      }
      else if (c == '(')
      {          // Open parenthesis.
        stackOper[stackOperIndex] = '(';  // Push operator onto stack.
        stackOperIndex++;
      }
      else if (cUppercase == variableLetter)
      {
        if (exponOperatorCounter != 0)
        {
          return EXPR_CANNOT_USE_X_IN_EXPONENT;
        }
        *ptrOutput = 'x';
        ptrOutput++;
        prevTokenIsNumber = true;
      }
      else if ((c >= '0') && (c <= '9'))
      {          // Number.
        ptrInput = pInput;
        while ((*ptrInput >= '0') && (*ptrInput <= '9'))
        {        // Find end of number.
          ptrInput++;
        }
        Dec2Bin(pInput - 1, value.limbs, (int)(ptrInput + 1 - pInput), &value.nbrLimbs);
        pInput = ptrInput;
        value.sign = SIGN_POSITIVE;
        if (exponOperatorCounter != 0)
        {
          if (value.nbrLimbs != 1)
          {
            return EXPR_EXPONENT_TOO_LARGE;
          }
          *ptrOutput = TOKEN_NUMBER;
          ptrOutput++;
          limb = value.limbs[0].x;
          for (bitNbr = BITS_PER_GROUP; bitNbr > 0; bitNbr -= 8)
          {
            *ptrOutput = (char)limb;
            ptrOutput++;
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
              (void)memset(&value.limbs[value.nbrLimbs], 0, (NumberLength - value.nbrLimbs) * sizeof(int));
            }
            // Convert to Montgomery notation.
            modmult(value.limbs, MontgomeryMultR2, value.limbs);
            value.nbrLimbs = NumberLength;
          }
          *ptrOutput = TOKEN_NUMBER;
          ptrOutput++;
          while (value.nbrLimbs > 1)
          {
            if (value.limbs[value.nbrLimbs - 1].x != 0)
            {
              break;
            }
            value.nbrLimbs--;
          }
          *ptrOutput = (char)(value.nbrLimbs >> 8);
          ptrOutput++;
          *ptrOutput = (char)value.nbrLimbs;
          ptrOutput++;
          for (int index = 0; index < value.nbrLimbs; index++)
          {
            limb = value.limbs[index].x;
            for (bitNbr = BITS_PER_GROUP; bitNbr > 0; bitNbr -= 8)
            {
              *ptrOutput = (char)limb;
              ptrOutput++;
              limb >>= 8;
            }
          }
        }
        prevTokenIsNumber = true;
      }
      else
      {
        return EXPR_SYNTAX_ERROR;
      }
    }
  }
  if (prevTokenIsNumber == false)
  {
    return EXPR_SYNTAX_ERROR;
  }
  while (stackOperIndex > 0)
  {      // Send operators to output.
    stackOperIndex--;
    s = stackOper[stackOperIndex];
    if (s == '^')
    {
      exponOperatorCounter--;
      if (exponOperatorCounter == 0)
      {
        *ptrOutput = TOKEN_END_EXPON;
        ptrOutput++;
      }
    }
    if (s == '(')
    {    // Operator on stack has less precedence.
      return EXPR_PAREN_MISMATCH;
    }
    *ptrOutput = s;
    ptrOutput++;
  }
  *ptrOutput = '\0';
  return EXPR_OK;
}

static int NegatePolynomialExpr(int* ptrArgument)
{
  int* ptrValue1 = ptrArgument;
  int val = *ptrValue1;
  if (val <= 0)
  {          // Monomial
    if ((*(ptrValue1 + 1) != 1) || (*(ptrValue1 + 2) != 0))
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
    int* ptrValue2 = poly1;
    ptrValue1++;
    for (int currentDegree = val; currentDegree >= 0; currentDegree--)
    {
      if ((*ptrValue1 == 1) && (*(ptrValue1 + 1) == 0))
      {          // If value is zero, it does not have to be changed.
        *ptrValue2 = 1;
        ptrValue2++;
        *ptrValue2 = 0;
        ptrValue2++;
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
      if ((ptrArgument + 1 + (ptrValue2 - &poly1[0])) > (values + sizeof(values) / sizeof(values[0])))
      {
        return EXPR_OUT_OF_MEMORY;
      }
      (void)memcpy(ptrArgument + 1, poly1, (ptrValue2 - &poly1[0]) * sizeof(int));
    }
  }
  return EXPR_OK;
}

void UncompressBigIntegerB(const int* ptrValues, BigInteger* bigint)
{
  if (modulusIsZero)
  {
    NumberLength = numLimbs(ptrValues);
  }
  IntArray2BigInteger(ptrValues, bigint);
}

// Add two polynomials starting on buffers pointed by ptrArgument1
// and ptrArgument2. Use poly1 as a temporary buffer to hold the sum.
static int AddPolynomialExpr(int* ptrArgument1, int* ptrArgument2)
{
  int* ptrValue1;
  const int* ptrValue2;
  int currentDegree;
  int degreeMin;
  int degreeMax;
  int degreePoly = 0;
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
        BigInteger2IntArray(ptrArgument1 + 1, &operand1);
      }
      else
      {
        AddBigNbrModN(operand1.limbs, operand2.limbs, operand1.limbs, powerMod.limbs, powerMod.nbrLimbs);
        BigInteger2IntArray(ptrArgument1 + 1, &operand1);
      }
      if (BigIntIsZero(&operand1))
      {                     // Sum is zero: set degree to zero.
        *ptrArgument1 = 0;
      }
      valuesIndex = (int)(ptrArgument1 + 2 + NumberLength - &values[0]);
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
        UncompressBigIntegerB(ptrArgument1 + 1, &operand2);
        UncompressBigIntegerB(ptrArgument2 + 1, &operand1);
      }
      else
      {
        degreeMax = degree2;
        degreeMin = degree1;
        UncompressBigIntegerB(ptrArgument1 + 1, &operand1);
        UncompressBigIntegerB(ptrArgument2 + 1, &operand2);
      }
      ptrValue1 = ptrArgument1;
      *ptrValue1 = degreeMax;
      ptrValue1++;
      for (currentDegree = 0; currentDegree < degreeMin; currentDegree++)
      {
        *ptrValue1 = 1;  // Number of limbs
        ptrValue1++;
        *ptrValue1 = 0;  // Value = zero.
        ptrValue1++;
      }
      BigInteger2IntArray(ptrValue1, &operand1);
      ptrValue1 += 1 + numLimbs(ptrValue1);
      for (currentDegree++; currentDegree < degreeMax; currentDegree++)
      {
        *ptrValue1 = 1;  // Number of limbs
        ptrValue1++;
        *ptrValue1 = 0;  // Value = zero.
        ptrValue1++;
      }
      NumberLength = operand2.nbrLimbs;
      BigInteger2IntArray(ptrValue1, &operand2);
      valuesIndex = (int)(ptrValue1 + 1 + numLimbs(ptrValue1) - &values[0]);
      return EXPR_OK;
    }
  }
  if ((degree1 > 0) && (degree2 > 0))
  {           // Sum of two polynomials.
    const int* ptrPolyMin;
    const int* ptrPolyMax;
    if (degree1 > degree2)
    {
      degreeMin = degree2;
      degreeMax = degree1;
      ptrPolyMin = ptrArgument2 + 1;
      ptrPolyMax = ptrArgument1 + 1;
    }
    else
    {
      degreeMin = degree1;
      degreeMax = degree2;
      ptrPolyMin = ptrArgument1 + 1;
      ptrPolyMax = ptrArgument2 + 1;
    }
    ptrValue1 = poly1;    // Point to temporary polynomial sum.
    *ptrValue1 = degreeMax;
    ptrValue1++;
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
      (void)memcpy(ptrValue1, ptrPolyMax, (1 + numLimbs(ptrPolyMax)) * sizeof(int));
      ptrPolyMax += 1 + numLimbs(ptrPolyMax);
      ptrValue1 += 1 + numLimbs(ptrValue1);
    }
    (void)memcpy(ptrArgument1, poly1, (ptrValue1 - &poly1[0]) * sizeof(int));
    degreePoly = degreeMax;   // New degree of polynomial.
  }
  else
  {                 // Sum of polynomial and monomial.
    int degreeMono;
    if (degree1 < 0)
    {
      UncompressBigIntegerB(ptrArgument1 + 1, &operand1);
      (void)memmove(ptrArgument1 + 1, ptrArgument1 + 3 + numLimbs(ptrArgument1 + 1), (&values[valuesIndex] - ptrArgument2) * sizeof(int));
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
      int* ptrDest = ptrArgument1 + 1;
      *ptrArgument1 = degreeMono;
      for (currentDegree = 0; currentDegree <= degreePoly; currentDegree++)
      {
        int numLen = 1 + numLimbs(ptrValue1);
        (void)memcpy(ptrDest, ptrValue1, numLen * sizeof(int));
        ptrValue1 += numLen;
        ptrDest += numLen;
      }
      for (; currentDegree < degreeMono; currentDegree++)
      {
        *ptrDest = 1;  // Number of limbs
        ptrDest++;
        *ptrDest = 0;  // Value = zero.
        ptrDest++;
      }
      BigInteger2IntArray(ptrDest, &operand1);
      degreePoly = degreeMono;   // Set new degree.
    }
    else
    {        // Degree of polynomial greater than degree of monomial.
      int differenceOfDegrees;
      int nbrLimbsSum;

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
        (void)memmove(ptrValue1 + differenceOfDegrees, ptrValue1, (ptrArgument2 - ptrValue1) * sizeof(int));
      }
      if (differenceOfDegrees < 0)
      {
        (void)memmove(ptrValue1, ptrValue1 - differenceOfDegrees, (ptrArgument2 - ptrValue1 - differenceOfDegrees) * sizeof(int));
      }
      (void)memcpy(ptrValue1, poly1, (1 + nbrLimbsSum) * sizeof(int));
    }
  }
  // Reduce degree if leading coefficient is zero.
  degreeMax = 0;
  ptrValue1 = ptrArgument1 + 1;
  ptrValue2 = ptrValue1;
  for (currentDegree = 0; currentDegree <= degreePoly; currentDegree++)
  {
    if ((*ptrValue1 != 1) || (*(ptrValue1 + 1) != 0))
    {                    // Coefficient is not zero
      degreeMax = currentDegree;
      // Store point to coefficient not zero of maximum degree.
      ptrValue2 = ptrValue1;
    }
    ptrValue1 += 1 + numLimbs(ptrValue1);
  }
  *ptrArgument1 = degreeMax;
  valuesIndex = (int)(ptrValue2 + 1 + numLimbs(ptrValue2) - &values[0]);
  return EXPR_OK;
}

int* CopyPolyProduct(const int* ptrSrc, int* ptrDest, int polyDegree)
{
  const int* ptrValueSrc = ptrSrc;
  int* ptrValueDest = ptrDest;
  for (int currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {
    int nbrLength = 1 + numLimbs(ptrValueSrc);
    (void)memcpy(ptrValueDest, ptrValueSrc, nbrLength * sizeof(int));
    ptrValueSrc += nbrLength;
    ptrValueDest += nbrLength;
  }
  return ptrValueDest;
}

static int MultPolynomialExpr(int* ptrArgument1, const int* ptrArgument2)
{
  int degreeMono;
  int degreePoly;
  int nbrLimbs = powerMod.nbrLimbs + 1;
  int* ptrValue1;
  int* ptrValue2;
  int currentDegree;
  int degree1 = *ptrArgument1;
  int degree2 = *ptrArgument2;
  const int* ptrValueSrc;
  if ((degree1 <= 0) && (degree2 <= 0))
  {        // Product of two monomials.
    if ((degree1 + degree2) < -MAX_DEGREE)
    {
      return EXPR_DEGREE_TOO_HIGH;
    }
    *ptrArgument1 = degree1 + degree2;
    UncompressBigIntegerB(ptrArgument1 + 1, &operand1);
    UncompressBigIntegerB(ptrArgument2 + 1, &operand2);
    if (modulusIsZero)
    {
      (void)BigIntMultiply(&operand1, &operand2, &operand1);
      NumberLength = operand1.nbrLimbs;
    }
    else
    {
      modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    }
    BigInteger2IntArray(ptrArgument1 + 1, &operand1);
    valuesIndex = (int)(ptrArgument1 + 2 + *(ptrArgument1 + 1) - &values[0]);
    return EXPR_OK;
  }
  if ((degree1 > 0) && (degree2 > 0))
  {        // Product of two polynomials.
    if ((degree1 + degree2) > MAX_DEGREE)
    {
      return EXPR_DEGREE_TOO_HIGH;
    }
    *ptrArgument1 = degree1 + degree2;
    if (modulusIsZero)
    {
      (void)CopyPolyProduct(ptrArgument1 + 1, poly1, degree1);  // Copy first factor to poly1
      (void)CopyPolyProduct(ptrArgument2 + 1, poly2, degree2);  // Copy second factor to poly2
    }
    else
    {
      // Copy first factor to poly1
      ptrValue1 = ptrArgument1 + 1;
      ptrValue2 = poly1;
      for (currentDegree = 0; currentDegree <= degree1; currentDegree++)
      {
        (void)memcpy(ptrValue2, ptrValue1, (1 + *ptrValue1) * sizeof(int));
        ptrValue1 += 1 + *ptrValue1;
        ptrValue2 += nbrLimbs;
      }
      // Copy second factor to poly2
      ptrValueSrc = ptrArgument2 + 1;
      ptrValue2 = poly2;
      for (currentDegree = 0; currentDegree <= degree2; currentDegree++)
      {
        (void)memcpy(ptrValue2, ptrValueSrc, (1 + *ptrValueSrc) * sizeof(int));
        ptrValueSrc += 1 + *ptrValueSrc;
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
      for (currentDegree = 0; currentDegree <= (degree1 + degree2); currentDegree++)
      {
        (void)memcpy(ptrValue1, ptrValue2, (1 + *ptrValue2) * sizeof(int));
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
    ptrValueSrc = ptrArgument2 + 1;
    // Get coefficient of monomial.
    UncompressBigIntegerB(ptrArgument1 + 1, &operand1);
    if ((degreeMono + degreePoly) > MAX_DEGREE)
    {
      return EXPR_DEGREE_TOO_HIGH;
    }
    // Multiply all coefficients of polynomial by the coefficient
    // of monomial storing the resulting polynomial on poly1.
    ptrValue2 = poly1;    // Initialize pointer to product.
    for (currentDegree = 0; currentDegree < degreeMono; currentDegree++)
    {
      *ptrValue2 = 1;
      ptrValue2++;
      *ptrValue2 = 0;
      ptrValue2++;
    }
    for (currentDegree = 0; currentDegree <= degreePoly; currentDegree++)
    {
      UncompressBigIntegerB(ptrValueSrc, &operand2);
      ptrValueSrc += 1 + numLimbs(ptrValueSrc);
      if (modulusIsZero)
      {
        (void)BigIntMultiply(&operand1, &operand2, &operand3);
        NumberLength = operand3.nbrLimbs;
      }
      else
      {
        modmult(operand1.limbs, operand2.limbs, operand3.limbs);
      }
      BigInteger2IntArray(ptrValue2, &operand3);
      ptrValue2 += 1 + numLimbs(ptrValue2);
    }
  }
  else
  {      // Second factor is monomial.
    degreeMono = -degree2;
    degreePoly = degree1;
    // Point to first coefficient of polynomial.
    ptrValue1 = ptrArgument1 + 1;
    // Get coefficient of monomial.
    UncompressBigIntegerB(ptrArgument2 + 1, &operand1);
    if ((degreeMono + degreePoly) > MAX_DEGREE)
    {
      return EXPR_DEGREE_TOO_HIGH;
    }
    // Multiply all coefficients of polynomial by the coefficient
    // of monomial storing the resulting polynomial on poly1.
    ptrValue2 = poly1;    // Initialize pointer to product.
    for (currentDegree = 0; currentDegree < degreeMono; currentDegree++)
    {
      *ptrValue2 = 1;
      ptrValue2++;
      *ptrValue2 = 0;
      ptrValue2++;
    }
    for (currentDegree = 0; currentDegree <= degreePoly; currentDegree++)
    {
      UncompressBigIntegerB(ptrValue1, &operand2);
      ptrValue1 += 1 + numLimbs(ptrValue1);
      if (modulusIsZero)
      {
        (void)BigIntMultiply(&operand1, &operand2, &operand3);
        NumberLength = operand3.nbrLimbs;
      }
      else
      {
        modmult(operand1.limbs, operand2.limbs, operand3.limbs);
      }
      BigInteger2IntArray(ptrValue2, &operand3);
      ptrValue2 += 1 + numLimbs(ptrValue2);
    }
  }
  *ptrArgument1 = degreeMono + degreePoly;
  ptrValue1 = ptrArgument1 + 1;
  (void)memcpy(ptrValue1, poly1, (ptrValue2 - &poly1[0]) * sizeof(int));
  valuesIndex = (int)(ptrValue1 - &values[0] + (ptrValue2 - &poly1[0]));
  return EXPR_OK;
}

void SetNumberToOne(/*@out@*/int* ptrValue1)
{
  int* ptrToValue1 = ptrValue1;
  const limb* destLimb;
  *ptrToValue1 = NumberLengthR1;
  ptrToValue1++;
  destLimb = MontgomeryMultR1;
  for (int ctr = 0; ctr < NumberLengthR1; ctr++)
  {
    *ptrToValue1 = destLimb->x;
    ptrToValue1++;
    destLimb++;
  }
}

int* CopyPolynomial(int* dest, const int* src, int polyDegree)
{
  int* ptrDest = dest;
  const int* ptrSrc = src;
  for (int currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {
    int numLength = numLimbs(ptrSrc) + 1;
    (void)memcpy(ptrDest, ptrSrc, numLength * sizeof(int));
    ptrSrc += numLength;
    ptrDest += numLength;
  }
  return ptrDest;
}

int* CopyPolynomialFixedCoeffSize(int* dest, const int* src, int polyDegree, int coeffSize)
{
  int* ptrDest = dest;
  const int* ptrSrc = src;
  for (int currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {
    int numLength = numLimbs(ptrSrc) + 1;
    (void)memcpy(ptrDest, ptrSrc, numLength * sizeof(int));
    ptrSrc += coeffSize;
    ptrDest += coeffSize;
  }
  return ptrDest;
}

static int PowerPolynomialExpr(int* ptrArgument1, int expon)
{
  limb exponLimb;
  int degreePower;
  int* ptrValue1;
  int* ptrValue2;
  int currentDegree;
  int nbrLimbs = powerMod.nbrLimbs + 1;
  int degreeBase = *ptrArgument1;
  if (degreeBase <= 0)
  {              // Monomial.
    if ((-degreeBase * expon) > MAX_DEGREE)
    {
      return EXPR_DEGREE_TOO_HIGH;
    }
    *ptrArgument1 = degreeBase * expon;
    UncompressBigIntegerB(ptrArgument1 + 1, &operand1);
    if (expon < 0)
    {
      return EXPR_EXPONENT_NEGATIVE;
    }
    if (modulusIsZero)
    {
      (void)BigIntPowerIntExp(&operand1, expon, &operand2);
      NumberLength = operand2.nbrLimbs;
    }
    else
    {
      exponLimb.x = expon;
      modPowLimb(operand1.limbs, &exponLimb, operand2.limbs);
    }
    BigInteger2IntArray(ptrArgument1 + 1, &operand2);
    valuesIndex = (int)(ptrArgument1 + 2 + numLimbs(ptrArgument1 + 1) - &values[0]);
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
      (void)memcpy(ptrValue2, ptrValue1, (1 + *ptrValue1) * sizeof(int));
      ptrValue1 += 1 + *ptrValue1;
      ptrValue2 += nbrLimbs;
    }
  }
  else
  {
    (void)CopyPolynomial(ptrValue2, ptrValue1, degreeBase);
  }
  SetNumberToOne(&poly2[0]); // Initialize power with polynomial 1.
  degreePower = 0;
  for (int mask = 1 << (BITS_PER_GROUP - 1); mask > 0; mask >>= 1)
  {
    // Square polynomial.
    MultPolynomial(degreePower, degreePower, poly2, poly2);
    degreePower <<= 1;
    if (modulusIsZero)
    {
      (void)CopyPolynomial(poly2, polyMultTemp, degreePower);
    }
    else
    {
      (void)memcpy(poly2, polyMultTemp, (degreePower + 1) * nbrLimbs * sizeof(int));
    }
    if ((expon & mask) != 0)
    {
      MultPolynomial(degreeBase, degreePower, poly1, poly2);
      degreePower += degreeBase;
      if (modulusIsZero)
      {
        (void)CopyPolynomial(poly2, polyMultTemp, degreePower);
      }
      else
      {
        (void)memcpy(poly2, polyMultTemp, (degreePower + 1) * nbrLimbs * sizeof(int));
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
      (void)memcpy(ptrValue1, ptrValue2, len * sizeof(int));
      ptrValue1 += len;
      ptrValue2 += nbrLimbs;
    }
  }
  valuesIndex = (int)(ptrValue1 - &values[0]);
  return EXPR_OK;
}

void computePower(int expo)
{
  (void)BigIntPowerIntExp(&primeMod, expo, &powerMod);
  (void)memcpy(TestNbr, powerMod.limbs, powerMod.nbrLimbs * sizeof(limb));
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
int ComputePolynomial(char* input, int expo)
{
  int* stackValues[STACK_OPER_SIZE];
  int stackIndex = 0;
  bool insideExpon = false;
  int* ptrValue1;
  int* ptrValue2;
  int len;
  int bitNbr;
  int bitCtr;
  int valueNbr;
  const char* ptrRPNbuffer;
  int rc;
  int val;
  int pwr;
  int expon;
  degree = 1;
  exponentMod = expo;
  // Use operand1 as temporary variable to store the exponent.
  computePower(expo);
  rc = ConvertToReversePolishNotation(input, RPNbuffer);
  if (rc != EXPR_OK)
  {
    return rc;
  }
  stackIndex = 0;
  valuesIndex = 0;
  ptrRPNbuffer = RPNbuffer;
  while (*ptrRPNbuffer != '\0')
  {
    switch (*ptrRPNbuffer)
    {
    case TOKEN_NUMBER:
      stackValues[stackIndex] = &values[valuesIndex];
      stackIndex++;
      if (insideExpon != 0)
      {                    // Inside exponent: only one limb accepted.
        bitCtr = 0;
        valueNbr = 0;
        for (bitNbr = BITS_PER_GROUP; bitNbr > 0; bitNbr -= 8)
        {
          ptrRPNbuffer++;
          valueNbr += (int)(unsigned char)*ptrRPNbuffer << bitCtr;
          bitCtr += 8;
        }
        values[valuesIndex] = valueNbr;
        valuesIndex++;
        break;
      }
      ptrRPNbuffer++;
      values[valuesIndex] = 0;   // Degree.
      valuesIndex++;
      len = ((int)(unsigned char)*(ptrRPNbuffer) * 256) + (unsigned char)*(ptrRPNbuffer + 1);
      values[valuesIndex] = len;
      valuesIndex++;
      ptrRPNbuffer++;
      for (int index = 0; index < len; index++)
      {
        bitCtr = 0;
        valueNbr = 0;
        for (bitNbr = BITS_PER_GROUP; bitNbr > 0; bitNbr -= 8)
        {
          ptrRPNbuffer++;
          valueNbr += (int)(unsigned char)*ptrRPNbuffer << bitCtr;
          bitCtr += 8;
        }
        values[valuesIndex] = valueNbr;
        valuesIndex++;
      }
      break;
    case 'x':
      stackValues[stackIndex] = &values[valuesIndex];
      stackIndex++;
      values[valuesIndex] = -1;   // Degree of monomial
      valuesIndex++;
         // Set coefficient to 1 in Montgomery notation.
      SetNumberToOne(&values[valuesIndex]);
      valuesIndex += 1 + values[valuesIndex];
      break;
    case TOKEN_START_EXPON:
      insideExpon = true;
      break;
    case TOKEN_END_EXPON:
      insideExpon = false;
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
      stackIndex--;
      ptrValue2 = stackValues[stackIndex];
      ptrValue1 = stackValues[stackIndex - 1];
      rc = GcdPolynomialExpr(ptrValue1, ptrValue2);
      if (rc != EXPR_OK)
      {
        return rc;
      }
      break;
    case TOKEN_DER:
      ptrValue1 = stackValues[stackIndex - 1];
      DerPolynomial(ptrValue1);
      break;
    case '+':
      stackIndex--;
      ptrValue2 = stackValues[stackIndex];
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
      stackIndex--;
      ptrValue2 = stackValues[stackIndex];
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
      stackIndex--;
      ptrValue2 = stackValues[stackIndex];
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
      stackIndex--;
      ptrValue2 = stackValues[stackIndex];
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
      stackIndex--;
      ptrValue2 = stackValues[stackIndex];
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
      stackIndex--;
      expon = *stackValues[stackIndex];
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
        for (int count = expon; count > 0; count--)
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
      break;
    default:
      break;
    }
    ptrRPNbuffer++;
  }
  // Adjust degree so the leading coefficient is not zero.
  if (values[0] < 0)
  {  // Monomial
    UncompressBigIntegerB(&values[1], &operand1);
    degree = -values[0];
    values[0] = degree;
    ptrValue1 = &values[1];
    for (int currentDegree = 0; currentDegree < degree; currentDegree++)
    {
      *ptrValue1 = 1;   // Initialize coefficient to zero.
      ptrValue1++;
      *ptrValue1 = 0;
      ptrValue1++;
    }
    BigInteger2IntArray(ptrValue1, &operand1);
  }
  return EXPR_OK;
}

static int GcdPolynomialExpr(int* ptrArgument1, int* ptrArgument2)
{
  int degreeGcd;
  const int* ptrNextArgument;
  if (modulusIsZero)
  {    // Integer GCD polynomial.
    int polyDegree;
    PolynomialGcd(ptrArgument1, ptrArgument2, ptrArgument1);
    polyDegree = *ptrArgument1;
    ptrNextArgument = ptrArgument1 + 1;
    for (int currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
    {
      ptrNextArgument += 1 + numLimbs(ptrNextArgument);
    }
  }
  else
  {    // Modular GCD polynomial.
    int degree1 = *ptrArgument1;
    int degree2 = *ptrArgument2;
    PolyModularGcd(ptrArgument1 + 1, degree1, ptrArgument2 + 1, degree2,
      poly5, &degreeGcd);
    ptrNextArgument = CopyPolynomial(ptrArgument1 + 1, poly5, degreeGcd);
    *ptrArgument1 = degreeGcd;
  }
  valuesIndex = (int)(ptrNextArgument - &values[0]);
  return EXPR_OK;
}


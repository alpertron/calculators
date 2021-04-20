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
#include "highlevel.h"
#include "polynomial.h"
#include "showtime.h"
#define STACK_OPER_SIZE            100
#define TOKEN_NUMBER               '0'
#define TOKEN_START_EXPON          '1'
#define TOKEN_END_EXPON            '2'
#define TOKEN_UNARY_MINUS          '3'
#define TOKEN_GCD                  '4'
#define TOKEN_DER                  '5'
extern char* ptrOutput2;
extern int DivPolynomialExpr(int* ptrArgument1, int* ptrArgument2, enum eDivType type);
static void showPolynomial(char** pptrOutput, int* ptrPoly, int polyDegree, int groupLength);
BigInteger primeMod;              // p
int exponentMod;                  // k
BigInteger powerMod;              // p^k
bool modulusIsZero;
int degree;
static BigInteger value;
BigInteger operand1;
BigInteger operand2;
BigInteger operand3;
BigInteger operand4;
BigInteger operand5;
extern limb TestNbr[MAX_LEN];
extern limb MontgomeryMultR2[MAX_LEN];
extern limb MontgomeryMultR1[MAX_LEN];
extern int NumberLength;
extern int NumberLengthR1;
static int prime;
static int primeIndex;
static char RPNbuffer[COMPRESSED_POLY_MAX_LENGTH];
int values[COMPRESSED_POLY_MAX_LENGTH];
int valuesIndex;
int valuesPrime[COMPRESSED_POLY_MAX_LENGTH];
int poly1[COMPRESSED_POLY_MAX_LENGTH];
int poly2[COMPRESSED_POLY_MAX_LENGTH];
int poly3[COMPRESSED_POLY_MAX_LENGTH];
int poly4[COMPRESSED_POLY_MAX_LENGTH];
int poly5[COMPRESSED_POLY_MAX_LENGTH];
int polyS[COMPRESSED_POLY_MAX_LENGTH];
static int polyT[COMPRESSED_POLY_MAX_LENGTH];
int polyMultTemp[COMPRESSED_POLY_MAX_LENGTH];
int polyLifted[COMPRESSED_POLY_MAX_LENGTH];
int polyLiftedNew[COMPRESSED_POLY_MAX_LENGTH];
enum eOutput pretty;
unsigned char onlyEvaluate = 0;
struct sFactorInfo factorInfo[MAX_DEGREE];
int nbrFactorsFound;
int *ptrOrigPoly;
int degreeOrigPoly;
static int *ptrA[MAX_DEGREE];
static int* ptrNewFactors[MAX_DEGREE];
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

static bool isFunc(char **ppcInput, char *funcName)
{
  char *pcInput = *ppcInput;
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
  bool prevTokenIsNumber = false;
  char s;
  char *ptrInput;
  int index;
  int limb;
  int bitNbr;
  char variableLetter = ' ';  // Indicate variable letter not known yet.
  while (*input != '\0')
  {
    char* inputTemp;
    char c = *input;
    char cUppercase;
    if (c == ' ' || c == 9)
    {          // Ignore any spaces and tabs.
      input++;
      continue;
    }
    inputTemp = input;
    if (isFunc(&inputTemp, "GCD"))
    {
      if (prevTokenIsNumber == false)
      {
        input = inputTemp;
        stackOper[stackOperIndex++] = TOKEN_GCD;  // Push token onto stack.
        continue;
      }
      c = '*';
    }
    else if (isFunc(&inputTemp, "DER"))
    {
      if (prevTokenIsNumber == false)
      {
        input = inputTemp;
        stackOper[stackOperIndex++] = TOKEN_DER;  // Push token onto stack.
        continue;
      }
      c = '*';
    }
    else
    {
      input++;
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
    if ((c == '*') && (*input == '*'))
    {          // Convert double asterisk to exponentiation.
      c = '^';
      input++;
    }
    if (prevTokenIsNumber)
    {
      prevTokenIsNumber = false;
      if (c == '+' || c == '-')
      {        // Binary plus or minus
        while (stackOperIndex > 0)
        {      // Send operators to output.
          s = stackOper[--stackOperIndex];
          if ((s != '+') && (s != '-') && (s != '*') && (s != '^') &&
            (s != '/') && (s != '%') && (s != TOKEN_UNARY_MINUS))
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
      else if (c == '*' || c == '(' || c == '/' || c == '%' ||
        cUppercase == variableLetter)
      {
        while (stackOperIndex > 0)
        {      // Send operators to output.
          s = stackOper[--stackOperIndex];
          if ((s != '^') && (s != '*') && (s != '/') && (s != '%'))
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
        if (cUppercase == variableLetter)
        {
          if (exponOperatorCounter != 0)
          {
            return EXPR_CANNOT_USE_X_IN_EXPONENT;
          }
          stackOper[stackOperIndex++] = '*';    // Push operator onto stack.
          *ptrOutput++ = 'x';
          prevTokenIsNumber = true;
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
          prevTokenIsNumber = false;
        }
        else
        {
          if (stackOperIndex > 0)
          {
            if (stackOper[stackOperIndex - 1] == TOKEN_GCD ||
              stackOper[stackOperIndex - 1] == TOKEN_DER)
            {
              *ptrOutput++ = stackOper[--stackOperIndex];
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
        stackOper[stackOperIndex++] = TOKEN_UNARY_MINUS;  // Push operator onto stack.
      }
      else if (c == '(')
      {          // Open parenthesis.
        stackOper[stackOperIndex++] = '(';  // Push operator onto stack.
      }
      else if (cUppercase == variableLetter)
      {
        if (exponOperatorCounter != 0)
        {
          return EXPR_CANNOT_USE_X_IN_EXPONENT;
        }
        *ptrOutput++ = 'x';
        prevTokenIsNumber = true;
      }
      else if ((c >= '0') && (c <= '9'))
      {          // Number.
        ptrInput = input;
        while ((*ptrInput >= '0') && (*ptrInput <= '9'))
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
              (void)memset(&value.limbs[value.nbrLimbs], 0, (NumberLength - value.nbrLimbs) * sizeof(int));
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
      if ((*ptrValue1 == 1) && (*(ptrValue1 + 1) == 0))
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
      (void)memcpy(ptrArgument + 1, poly1, (ptrValue2 - &poly1[0]) * sizeof(int));
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

// Add two polynomials starting on buffers pointed by ptrArgument1
// and ptrArgument2. Use poly1 as a temporary buffer to hold the sum.
static int AddPolynomialExpr(int *ptrArgument1, int *ptrArgument2)
{
  int *ptrValue1;
  int *ptrValue2;
  int currentDegree;
  int degreeMin;
  int degreeMax;
  int degreePoly=0;
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
  if ((degree1 > 0) && (degree2 > 0))
  {           // Sum of two polynomials.
    int *ptrPolyMin;
    int *ptrPolyMax;
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
    ptrValue1 = poly1;    // Point to temporary polynomial sum.
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
      (void)memcpy(ptrValue1, ptrPolyMax, (1 + numLimbs(ptrPolyMax))*sizeof(int));
      ptrPolyMax += 1 + numLimbs(ptrPolyMax);
      ptrValue1 += 1 + numLimbs(ptrValue1);
    }
    (void)memcpy(ptrArgument1, poly1, (ptrValue1 - &poly1[0])*sizeof(int));
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
      int* ptrValue2 = ptrArgument1 + 1;
      *ptrArgument1 = degreeMono;
      for (currentDegree = 0; currentDegree <= degreePoly; currentDegree++)
      {
        int numLen = 1 + numLimbs(ptrValue1);
        (void)memcpy(ptrValue2, ptrValue1, numLen * sizeof(int));
        ptrValue1 += numLen;
        ptrValue2 += numLen;
      }
      for (; currentDegree < degreeMono; currentDegree++)
      {
        *ptrValue2++ = 1;  // Number of limbs
        *ptrValue2++ = 0;  // Value = zero.
      }
      BigInteger2IntArray(ptrValue2, &operand1);
      degreePoly = degreeMono;   // Set new degree.
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
        (void)memmove(ptrValue1 + differenceOfDegrees, ptrValue1, (ptrArgument2 - ptrValue1) * sizeof(int));
      }
      else if (differenceOfDegrees < 0)
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
    if (*ptrValue1 != 1 || *(ptrValue1 + 1) != 0)
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

int *CopyPolyProduct(int *ptrSrc, int *ptrDest, int degree)
{
  int currentDegree;
  int *ptrValue1 = ptrSrc;
  int *ptrValue2 = ptrDest;
  for (currentDegree = 0; currentDegree <= degree; currentDegree++)
  {
    int nbrLength = 1 + numLimbs(ptrValue1);
    (void)memcpy(ptrValue2, ptrValue1, nbrLength * sizeof(int));
    ptrValue1 += nbrLength;
    ptrValue2 += nbrLength;
  }
  return ptrValue1;
}

static int MultPolynomialExpr(int *ptrArgument1, int *ptrArgument2)
{
  int degreeMono;
  int degreePoly;
  int nbrLimbs = powerMod.nbrLimbs + 1;
  int *ptrValue1;
  int *ptrValue2;
  int currentDegree;
  int degree1 = *ptrArgument1;
  int degree2 = *ptrArgument2;
  if ((degree1 <= 0) && (degree2 <= 0))
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
  if ((degree1 > 0) && (degree2 > 0))
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
        (void)memcpy(ptrValue2, ptrValue1, (1 + *ptrValue1) * sizeof(int));
        ptrValue1 += 1 + *ptrValue1;
        ptrValue2 += nbrLimbs;
      }
      // Copy second factor to poly2
      ptrValue1 = ptrArgument2 + 1;
      ptrValue2 = poly2;
      for (currentDegree = 0; currentDegree <= degree2; currentDegree++)
      {
        (void)memcpy(ptrValue2, ptrValue1, (1 + *ptrValue1) * sizeof(int));
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
  ptrValue2 = poly1;    // Initialize pointer to product.
  for (currentDegree = 0; currentDegree < degreeMono; currentDegree++)
  {
    *ptrValue2++ = 1;
    *ptrValue2++ = 0;
  }
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
    ptrValue1 += 1 + numLimbs(ptrValue1);
    ptrValue2 += 1 + numLimbs(ptrValue2);
  }
  *ptrArgument1 = degreeMono + degreePoly;
  ptrValue1 = ptrArgument1 + 1;
  (void)memcpy(ptrValue1, poly1, (ptrValue2 - &poly1[0])*sizeof(int));
  valuesIndex = (int)(ptrValue1 - &values[0] + (ptrValue2 - &poly1[0]));
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
    (void)memcpy(ptrDest, ptrSrc, numLength * sizeof(int));
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
    (void)memcpy(ptrDest, ptrSrc, numLength * sizeof(int));
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
  int *ptrValue1;
  int *ptrValue2;
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
      (void)memcpy(ptrValue2, ptrValue1, (1 + *ptrValue1) * sizeof(int));
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
      (void)memcpy(poly2, polyMultTemp, (degreePower + 1)*nbrLimbs * sizeof(int));
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
        (void)memcpy(poly2, polyMultTemp, (degreePower + 1)*nbrLimbs * sizeof(int));
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
  BigIntPowerIntExp(&primeMod, expo, &powerMod);
  (void)memcpy(TestNbr, powerMod.limbs, powerMod.nbrLimbs*sizeof(limb));
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
  bool insideExpon = false;
  int *ptrValue1;
  int *ptrValue2;
  int len;
  int index;
  int bitNbr;
  int bitCtr;
  int valueNbr;
  char *ptrRPNbuffer;
  int rc;
  int count;
  int val;
  int pwr;
  int expon;
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
    int* ptrPoly = poly + 1;
    intToBigInteger(&operand1, inverse);
    if (TestNbr[0].x <= 32768)
    {
      for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
      {
        *ptrPoly = *ptrPoly * inverse % TestNbr[0].x;
        ptrPoly += 2;
      }
    }
    else
    {
      for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
      {
#ifdef _USING64BITS_
        *ptrPoly = (int)((uint64_t)*ptrPoly * inverse % TestNbr[0].x);
#else
        smallmodmult(*ptrPoly, inverse, (limb*)ptrPoly, TestNbr[0].x);
#endif
        ptrPoly += 2;
      }
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
    NumberLength = numLimbs(poly);
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
    NumberLength = numLimbs(poly);
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
  static BigInteger contentF;
  static BigInteger contentG;
  static BigInteger contentH;
  static BigInteger modulus;
  static BigInteger gcdLeadingCoeff;
  int *ptrLeadingCoeff;
  int *ptrSrc;
  int *ptrDest;
  int *ptrPrev;
  int potentialDegreeGcd;
  int tmp;
  int rc;
  initializeSmallPrimes(smallPrimes);
  prime = 65537;
  primeIndex = 6543;
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
      prime = smallPrimes[++primeIndex];
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
    if ((polyS[0] == 0) && (polyS[1] == 1) && (polyS[2] == 0))
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
    if ((rc == EXPR_OK) && (poly5[0] == 0) && (poly5[1] == 1) && (poly5[2] == 0))
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
      if ((rc == EXPR_OK) && (poly5[0] == 0) && (poly5[1] == 1) && (poly5[2] == 0))
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
  int temp;
  int index;
  int currentDegree;
  if (degree2 == 0)
  {
    if ((*arg2 == 1) && (*(arg2 + 1) == 0))
    {     // Number is zero. GCD is first argument.
      *degreeGcd = degree1;
      (void)memcpy(gcd, arg1, (degree1+1)*nbrLimbs*sizeof(int));
      return;
    }
    *degreeGcd = 0;
    return;
  }
  (void)memcpy(gcd, arg1, (degree1+1)*nbrLimbs*sizeof(int));
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
        int value = modulus - *(ptrArgMax + currentDegree * 2 + 1);
        int *ptrPoly = ptrArgMin + 1;
        ptrTemp++;
        if (modulus < 32767)
        {
          for (index = 0; index <= degreeMin; index++)
          {
            *ptrTemp = (*ptrTemp + *ptrPoly * value) % modulus;
            ptrTemp += 2;
            ptrPoly += 2;
          }
        }
        else
        {
          for (index = 0; index <= degreeMin; index++)
          {
#ifdef _USING64BITS_
            *ptrTemp = (*ptrTemp + (uint64_t)*ptrPoly * value) % modulus;
#else
            smallmodmult(*ptrPoly, value, (limb*)&temp, modulus);
            *ptrTemp = (int)(((unsigned int)*ptrTemp + (unsigned int)temp) %
              (unsigned int)modulus);
#endif
            ptrTemp += 2;
            ptrPoly += 2;
          }
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
        (void)memcpy(gcd, ptrArgMin, (degreeMin+1)*nbrLimbs*sizeof(int));
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
  int* ptrNextArgument;
  if (modulusIsZero)
  {    // Integer GCD polynomial.
    int currentDegree, degree;
    PolynomialGcd(ptrArgument1, ptrArgument2, ptrArgument1);
    degree = *ptrArgument1;
    ptrNextArgument = ptrArgument1 + 1;
    for (currentDegree = 0; currentDegree <= degree; currentDegree++)
    {
      ptrNextArgument += 1 + numLimbs(ptrNextArgument);
    }
  }
  else
  {    // Modular GCD polynomial.
    int degree1 = *ptrArgument1;
    int degree2 = *ptrArgument2;
    PolyModularGcd(ptrArgument1+1, degree1, ptrArgument2 + 1, degree2,
      poly5, &degreeGcd);
    ptrNextArgument = CopyPolynomial(ptrArgument1 + 1, poly5, degreeGcd);
    *ptrArgument1 = degreeGcd;
  }
  valuesIndex = (int)(ptrNextArgument - &values[0]);
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
      if ((*(ptrArgument + 1) == 1) && (*(ptrArgument + 2) == 0))
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
    (void)memset(&operand2.limbs[1], 0, (NumberLength - 1) * sizeof(limb));
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

  // Perform polyPower <- polyBase ^ expon (mod polyMod)
void powerPolynomial(int *polyBase, int *polyMod, int polyDegree, BigInteger *expon,
                     int *polyPower, powerCback callback, int curMultip, int nbrMultip)
{
  int mask;
  int nbrLimbs = NumberLength + 1;
  bool powerIsOne = true;
  int nbrBits = 0;
  int index = expon->nbrLimbs - 1;
  int bitCounter = 0;
  *(polyBase + polyDegree * nbrLimbs) = 1;
  *(polyBase + polyDegree * nbrLimbs + 1) = 0;
  for (; index >= 0; index--)
  {
    int groupExp = (int)(expon->limbs[index].x);
    for (mask = 1 << (BITS_PER_GROUP - 1); mask > 0; mask >>= 1)
    {
      if (!powerIsOne)
      {
        if (callback)
        {
          callback(100 * bitCounter/(nbrBits*nbrMultip));
        }
        multUsingInvPolynomial(polyPower, polyPower, polyPower, polyDegree, polyMod);
      }
      if ((groupExp & mask) != 0)
      {
        if (powerIsOne)
        {
          CopyPolynomialFixedCoeffSize(polyPower, polyBase, polyDegree, nbrLimbs);
          nbrBits = (index + 1) * BITS_PER_GROUP - bitCounter;
          powerIsOne = false;
          bitCounter = curMultip * nbrBits;
        }
        else
        {
          multUsingInvPolynomial(polyPower, polyBase, polyPower, polyDegree, polyMod);
        }
      }
      bitCounter++;
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
    (void)memcpy(poly4, poly, (polyDegree + 1)*nbrLimbs*sizeof(int));      // Backup poly
    if (degreeC == 0)
    {
      (void)memcpy(poly1, poly, (polyDegree + 1)*nbrLimbs*sizeof(int));           // poly1 = w
    }
    else
    {
      SetNumberToOne(&poly2[degreeC*nbrLimbs]);
      DividePolynomial(poly4, polyDegree, poly2, degreeC, poly1);    // poly1 = w
    }
    degreeW = polyDegree - degreeC;
    while (degreeW != 0)
    {
      (void)memcpy(poly4, poly1, (degreeW+1)*nbrLimbs*sizeof(int));        // Backup w.
      PolyModularGcd(poly2, degreeC, poly1, degreeW, poly3, &degreeY);      // poly3 = y
      SetNumberToOne(&poly3[degreeY*nbrLimbs]);
      if (degreeW != degreeY)
      {
        int degreeZ;
        struct sFactorInfo *pstFactorInfo;

        DividePolynomial(poly4, degreeW, poly3, degreeY, poly1);     // poly1 = z
        // z^i is divisor of the original polynomial.
        degreeZ = degreeW - degreeY;
        for (currentDegree = 0; currentDegree < i*expon; currentDegree++)
        {
          DividePolynomial(ptrOrigPoly, degreeOrigPoly, poly1, degreeZ, poly4);
          degreeOrigPoly -= degreeZ;
          (void)memcpy(ptrOrigPoly, poly4, (degreeOrigPoly+1)*nbrLimbs*sizeof(int));
        }
        pstFactorInfo = &factorInfo[nbrFactorsFound++];
        pstFactorInfo -> ptr = ptrOrigPoly;
        pstFactorInfo -> degree = degreeZ;
        pstFactorInfo -> multiplicity = i*expon;
        pstFactorInfo -> expectedDegree = 0;    // Unknown at this moment.
        (void)memcpy(ptrOrigPoly, poly1, degreeZ*nbrLimbs*sizeof(int));
        ptrOrigPoly += (degreeZ+1)*nbrLimbs;
        (void)memcpy(ptrOrigPoly, poly4, (degreeOrigPoly + 1)*nbrLimbs*sizeof(int));
      }
      i++;
      (void)memcpy(poly1, poly3, (degreeY+1)*nbrLimbs*sizeof(int));  // Copy y to w.
      degreeW = getDegreePoly(poly1, degreeY);
      DividePolynomial(poly2, degreeC, poly1, degreeW, poly4); // Compute c.
      degreeC -= degreeW;
      (void)memcpy(poly2, poly4, (degreeC+1)*nbrLimbs*sizeof(int));
    }
    if (degreeC != 0)
    {      // C is a perfect power.
      ptrValue1 = ptrOrigPoly + (degreeOrigPoly + 1)*nbrLimbs;
      ptrValue2 = poly2;
      for (currentDegree = 0; currentDegree <= degreeC; currentDegree += primeInt)
      {
        (void)memcpy(ptrValue1, ptrValue2, nbrLimbs*sizeof(int));
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
      (void)memcpy(ptrValue1, ptrValue2, nbrLimbs*sizeof(int));
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

// Given coprime polynomials A and B, find polynomials U and V such that
// A*U + B*V = 1 (mod powerMod). Polynomial V is not required on output.
// Use Montgomery notation.
static void ExtendedGcdPolynomial(/*@in@*/int *ptrA, int degreeA, /*@in@*/int *ptrB,
  int degreeB, /*@out@*/int *ptrQ,
  /*@out@*/int *ptrP1, /*@out@*/int *ptrP2, /*@out@*/int *ptrU, /*@out@*/int *pDegreeU)
{
  int *ptrQuotient = ptrQ;
  int *ptrQuotients[MAX_DEGREE+1];
  int *tmpPtr, *ptrR, *ptrOldR;
  int degreeR, degreeOldR;
  int nbrQuotients = 0;
  int nbrLimbs = NumberLength + 1;
  bool polyExchanged = false;
  int currentDegree;
  int degreeU;
  int counter;

  // Ensure that degree of A is greater than degree of B.
  // Variable polyExchange indicates whether the polynomials were exchanged or not.
	if (degreeA < degreeB)
	{
    int tmpDegree = degreeA;
    degreeA = degreeB;
    degreeB = tmpDegree;
    tmpPtr = ptrA;
    ptrA = ptrB;
    ptrB = tmpPtr;
    polyExchanged = true;
	}
  ptrR = ptrP1;
  ptrOldR = ptrP2;
  (void)memcpy(ptrOldR, ptrA, (degreeA + 1)*nbrLimbs*sizeof(int));
  (void)memcpy(ptrR, ptrB, (degreeB + 1)*nbrLimbs*sizeof(int));
  degreeR = degreeB;
  degreeOldR = degreeA;
     // Loop that computes all quotients.
     // Since gcd = 1, the loop will finish when the degree of R is zero.
  while (degreeR != 0)
  {  // If R is not constant...
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
  LenAndLimbs2ArrLimbs(ptrR, operand5.limbs, nbrLimbs);
  // Save pointer to last quotient.
  ptrQuotients[nbrQuotients] = ptrQuotient;
  // P1 <- 0
  *ptrP1 = 1;
  *(ptrP1 + 1) = 0;
  // P2 <- modular inverse of gcd (stored in operand5).
  nbrLimbs = NumberLength;
  ModInvBigNbr(operand5.limbs, (limb*)(ptrP2 + 1), powerMod.limbs, powerMod.nbrLimbs);
  while ((nbrLimbs > 1) && (*(ptrP2 + nbrLimbs) == 0))
  {
    nbrLimbs--;
  }
  *ptrP2 = nbrLimbs;
  if (polyExchanged)
  {
    ptrR = ptrP2;
    ptrOldR = ptrP1;
  }
  else
  {
    ptrR = ptrP1;
    ptrOldR = ptrP2;
  }
  degreeR = 0;
  degreeOldR = 0;
  nbrLimbs = NumberLength + 1;
  for (counter = 0; counter < nbrQuotients; counter++)
  {
    int tmpDegree, offset;
    int degreeQ = (int)((ptrQuotients[counter + 1] - ptrQuotients[counter]) / nbrLimbs - 1);
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
      (void)memset(operand1.limbs, 0, NumberLength * sizeof(limb));
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
  degreeU = degreeR;
  (void)memcpy(ptrU, ptrR, (degreeU + 1) * nbrLimbs * sizeof(int));
  *pDegreeU = degreeU;
}

// Lift factorization of coefficients of factors from mod prime to
// mod power = prime^exp.
// Perform lifting mod prime, prime^2, prime^4, etc.
// All calculations are done in Montgomery domain, so the roots must 
// be converted to standard notation and then to Montgomery notation
// with greater modulus.
//
// Let f(x) = f_1(x) * f_2(x) * ... * f_n(x)
// 
// Compute 1 = (f(x)/f_1(x)) * a1(x) + ... + (f(x)/f_n(x)) * a_n(x)
// This is done by computing n extended GCDs between f_i(x) and
// f(x)/f_i(x).
//
// Hensel lift from mod m to mod m^2:
// u(x) <- (1/m)*(f - f_1*f_2*...*f_n)
// For i = 1 to n compute:
//   g(x) <- u(x) * a_i(x) mod f_i(x) (all computations done mod m).
//   f_i(x) <- f_i(x) + m*g(x)
//
// Compute u(x) <- (1/m)*(f(x)/f_1(x)*a_1(x) + ... + f(x)/f_n(x)*a_n(x) - 1).
// For i = 1 to n compute:
//   g(x) <- u(x) * a_i(x) mod f_i(x) (all computations done mod m).
//   a_i(x) <- a_i(x) - m*g(x)
//
// Memory usage:
// poly3: u
// poly4: f, g
// poly5: m*g, f/f_i*a_i
// polyT: a_i
// polyLifted: f_i.
// Change the fields ptrPolyLifted to point to f_i.

// Compute polynomial poly4 <- values divided by the leading coefficient
// mod prime. This is the polynomial f.
static void ComputeF(void)
{
  int degree, currentDegree;
  int* ptrValue1;
  degree = values[0];                              // Get degree of polynomial
  // Find the leading coefficient.
  ptrValue1 = &values[1];                          // Point to constant coefficient.
  for (currentDegree = 0; currentDegree < degree; currentDegree++)
  {
    ptrValue1 += 1 + numLimbs(ptrValue1);
  }
  NumberLength = numLimbs(ptrValue1);
  IntArray2BigInteger(ptrValue1, &operand1);
  (void)BigIntRemainder(&operand1, &powerMod, &operand1);
  NumberLength = powerMod.nbrLimbs;
  if (operand1.nbrLimbs < NumberLength)
  {
    (void)memset(&operand1.limbs[operand1.nbrLimbs], 0,
      (NumberLength - operand1.nbrLimbs) * sizeof(limb));
  }
  // Compute the inverse of leading coefficient.
  ModInvBigNbr(operand1.limbs, operand2.limbs, TestNbr, NumberLength);
  // Convert operand1 from standard to Montgomery notation.
  ptrValue1 = &values[1];                          // Point to constant coefficient.
  for (currentDegree = 0; currentDegree <= degree; currentDegree++)
  {
    NumberLength = numLimbs(ptrValue1);
    IntArray2BigInteger(ptrValue1, &operand1);
    (void)BigIntRemainder(&operand1, &powerMod, &operand1);
    NumberLength = powerMod.nbrLimbs;
    if (operand1.nbrLimbs < NumberLength)
    {
      (void)memset(&operand1.limbs[operand1.nbrLimbs], 0,
        (NumberLength - operand1.nbrLimbs) * sizeof(limb));
    }
    // Multiply by inverse. Result is in Montgomery notation.
    modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    BigInteger2IntArray(&poly4[currentDegree * (NumberLength + 1)], &operand1);
    ptrValue1 += 1 + numLimbs(ptrValue1);           // Point to next coefficient.
  }
}

// Get polynomial a_i. This polynomial is stored with no spaces
// between coefficients. The output is in Montgomery notation in
// polynomial poly2.
static int getAi(int nbrFactor, int degreeA)
{
  int nbrLimbs = NumberLength+1;
  int currentDegree;
  int* ptrSrc = ptrA[nbrFactor];
  int* ptrDest = poly1;
  for (currentDegree = 0; currentDegree <= degreeA; currentDegree++)
  {
    int numLen = *ptrSrc + 1;
    (void)memcpy(ptrDest, ptrSrc, numLen * sizeof(int));
    ptrDest += nbrLimbs;
    ptrSrc += numLen;
  }
  polyToMontgomeryNotation(poly1, degreeA+1);
  return degreeA;
}

// Get polynomial factor f_i. This factor is stored with no spaces
// between coefficients. The output is in Montgomery notation in
// polynomial poly2.
static void getFi(int degreeFactor, int *ptrSrc, int nbrLimbs)
{
  int *ptrDest = poly2;       // Copy f_i to poly2.
  int currentDegree;
  for (currentDegree = 0; currentDegree < degreeFactor; currentDegree++)
  {
    int nbrLen = *ptrSrc + 1;
    (void)memcpy(ptrDest, ptrSrc, nbrLen * sizeof(int));
    ptrSrc += nbrLen;
    ptrDest += nbrLimbs;
  }
  SetNumberToOne(ptrDest);
  polyToMontgomeryNotation(poly2, degreeFactor);
}

// Move factors f_i from polyLiftedNew to polyLifted, then adjust
// factorInfo[...].ptrPolyLifted to point to the new factors.
// There are no spaces between coefficients.
static void MoveFactorsAndFixPointers(struct sFactorInfo* factorInfo, int compressPoly)
{
  int nbrFactor, currentDegree;
  struct sFactorInfo* pstFactorInfo = factorInfo;
  int* ptrSrc = polyLiftedNew;
  int * ptrDest = polyLifted;
  for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    int degreeFactor = pstFactorInfo->degree;
    pstFactorInfo->ptrPolyLifted = ptrDest;
    for (currentDegree = 0; currentDegree < degreeFactor; currentDegree++)
    {
      int numLength = *ptrSrc + 1;
      (void)memcpy(ptrDest, ptrSrc, numLength * sizeof(int));
      ptrSrc += numLength;
      if (compressPoly)
      {
        ptrDest += numLength;
      }
      else
      {
        ptrDest += NumberLength+1;
      }
    }
    pstFactorInfo++;
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

int HenselLifting(struct sFactorInfo* factorInfo, int compressPoly)
{
  int currentExp = 1;
  int nbrFactor;
  int oldNumberLength, newNumberLength;
  int degreeFactor, degreeS, currentDegree;
  int index, nbrLimbs;
  int* ptrPolyLifted = polyLifted;
  int* ptrPoly2;
  int** ptrNewFactor;
  struct sFactorInfo *pstFactorInfo;
  // Copy polynomials f_i(x) to polyLifted.
  int* ptrDest = ptrPolyLifted;
  pstFactorInfo = factorInfo;
  for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    pstFactorInfo->ptrPolyLifted = ptrDest;
    ptrDest = CopyPolynomialFixedCoeffSize(ptrDest, pstFactorInfo->ptr,
      pstFactorInfo->degree, NumberLength+1);
    pstFactorInfo++;
  }
  if (exponentMod == 1)
  {      // No lift to be done. Go out.
    return EXPR_OK;
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
  computePower(1);
  // Compute a_i for each factor f_i by computing the extended gcd between
  // f_i(x) and f(x)/f_i(x).
  ptrDest = polyT;
  pstFactorInfo = factorInfo;
  for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    ComputeF();                       // poly4 <- f mod prime (Montgomery notation).
    nbrLimbs = NumberLength + 1;
    degreeFactor = pstFactorInfo->degree;
    (void)memcpy(poly2, pstFactorInfo->ptr, pstFactorInfo->degree * nbrLimbs * sizeof(int));
    // The factors do not include their leading coefficient. Set it to one.
    SetNumberToOne(&poly2[pstFactorInfo->degree * nbrLimbs]);
    DividePolynomial(poly4, degree,   // Dividend f(x).
      poly2, degreeFactor,            // Divisor f_i(x).
      poly1);                         // Quotient f(x)/f_i(x).
    ExtendedGcdPolynomial(poly1, degree - degreeFactor, poly2, degreeFactor,
      poly3, poly4, poly5, polyS, &degreeS);
    ptrA[nbrFactor] = ptrDest;
    ptrDest = CopyPolynomial(ptrDest, polyS, degreeS);
    for (currentDegree = degreeS + 1; currentDegree < degreeFactor; currentDegree++)
    {
      *ptrDest++ = 1;
      *ptrDest++ = 0;
    }
    pstFactorInfo++;
  }
  // Loop that performs the lifting.
  while (currentExp != exponentMod)
  {      // Hensel lift from mod m to mod m^2:
    int newExponent;
    int* ptrSrc;
    int degreeA;
#ifdef __EMSCRIPTEN__
    int elapsedTime = (int)(tenths() - originalTenthSecond);
    if (elapsedTime / 10 != oldTimeElapsed / 10)
    {
      char outputInfo[1000];
      char* ptrOutput = outputInfo;
      if (lang)
      {
        (void)strcpy(ptrOutput, "1<p>Aplicando lema de Hensel usando el número primo ");
        ptrOutput += strlen(ptrOutput);
        int2dec(&ptrOutput, primeMod.limbs[0].x);
        (void)strcpy(ptrOutput, " procesando exponente ");
        ptrOutput += strlen(ptrOutput);
        int2dec(&ptrOutput, currentExp);
        (void)strcpy(ptrOutput, " de ");
      }
      else
      {
        (void)strcpy(ptrOutput, "1<p>Hensel lifting using prime number ");
        ptrOutput += strlen(ptrOutput);
        int2dec(&ptrOutput, primeMod.limbs[0].x);
        (void)strcpy(ptrOutput, " processing exponent ");
        ptrOutput += strlen(ptrOutput);
        int2dec(&ptrOutput, currentExp);
        (void)strcpy(ptrOutput, " of ");
      }
      ptrOutput += strlen(ptrOutput);
      int2dec(&ptrOutput, exponentMod);
      (void)strcpy(ptrOutput, ".</p>");
      ptrOutput += strlen(ptrOutput);
      showElapsedTimeSec(&ptrOutput);
      databack(outputInfo);
    }
#endif
    // Compute u(x) <- (1/m) * (f - f_1 * f_2 *...* f_n)
    oldNumberLength = NumberLength;
    newExponent = currentExp << 1;     // We can double the exponent in each step.
    if (newExponent > exponentMod)     // Do not exceed exponentMod.
    {
      newExponent = exponentMod;
    }
    BigIntPowerIntExp(&primeMod, currentExp, &operand5);
    computePower(newExponent);         // Compute powerMod and init Montgomery parms.
    newNumberLength = NumberLength;
    nbrLimbs = NumberLength + 1;
    ComputeF();                        // poly4 <- f mod m^2 (Montgomery notation).
    SetNumberToOne(poly1);             // Initialize product of factors.
    pstFactorInfo = factorInfo;
    degree = 0;
    for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
    {    // Loop that multiplies all factors. poly1 holds the result.
      degreeFactor = pstFactorInfo->degree;
      getFi(degreeFactor, pstFactorInfo->ptrPolyLifted, nbrLimbs);   // Copy f_i to poly2.
      MultPolynomial(degree, degreeFactor, poly1, poly2);
      degree += degreeFactor;
      (void)memcpy(poly1, polyMultTemp, (degree + 1) * nbrLimbs * sizeof(int));
      pstFactorInfo++;
    }
    polyToStandardNotation(poly1, degree + 1);
    polyToStandardNotation(poly4, degree + 1);
    for (currentDegree = 0; currentDegree <= degree; currentDegree++)
    {                 // Loop that computes (1/m)*(f - f_1 * f_2 * ... * f_n)
      int nbrLen;
      // Get coefficient of f.
      IntArray2BigInteger(&poly4[currentDegree * nbrLimbs], &operand1);
      // Get coefficient of product of factors.
      IntArray2BigInteger(&poly1[currentDegree * nbrLimbs], &operand2);
      SubtBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
      // Get number of significant limbs before performing division.
      nbrLen = NumberLength;
      while (nbrLen > 1)
      {
        if (operand1.limbs[nbrLen - 1].x != 0)
        {
          break;
        }
        nbrLen--;
      }
      operand1.nbrLimbs = nbrLen;
      BigIntDivide(&operand1, &operand5, &operand1);
      // Store coefficient of subtraction.
      ptrDest = &poly3[currentDegree * (oldNumberLength + 1)];
      *ptrDest++ = operand1.nbrLimbs;
      (void)memcpy(ptrDest, operand1.limbs, operand1.nbrLimbs * sizeof(int));
    }
    computePower(currentExp);
    polyToMontgomeryNotation(poly3, degree+1);
    nbrLimbs = NumberLength + 1;
    pstFactorInfo = factorInfo;
    ptrDest = polyLiftedNew;
    ptrNewFactor = ptrNewFactors;
    for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
    {
      degreeFactor = pstFactorInfo->degree;
      // g(x) <- u(x) * a_i(x) mod f_i(x) (all computations done mod m).
      degreeA = getAi(nbrFactor, degreeFactor - 1); // poly1 <- a_i.
      MultPolynomial(degree, degreeA, poly3, poly1);
      (void)memcpy(poly1, polyMultTemp, (degree + degreeA + 1)* nbrLimbs * sizeof(int));
                                                 // poly1 <- u * a_i
                                                 // Copy f_i to poly2.
      getFi(degreeFactor, pstFactorInfo->ptrPolyLifted, nbrLimbs);
      DividePolynomial(poly1, degree + degreeA,  // Dividend = u * a_i
        poly2, degreeFactor,                     // Divisor = f_i
        NULL);                                   // Quotient: not needed.
      polyToStandardNotation(poly1, degreeFactor);  // Convert g to standard notation.
      // f_i(x) <- f_i(x) + m*g(x)
      int *ptrFi = pstFactorInfo->ptrPolyLifted;
      *(ptrNewFactor++) = ptrDest;                 // Point to f_i(x) mod new modulus.
      for (currentDegree = 0; currentDegree < degreeFactor; currentDegree++)
      {                 // Loop that computes f_i + m*g
        // Get coefficient of g.
        ptrSrc = &poly1[currentDegree * nbrLimbs];
        operand1.nbrLimbs = *ptrSrc++;
        (void)memcpy(operand1.limbs, ptrSrc, operand1.nbrLimbs * sizeof(int));
        if (currentExp * 2 > newExponent)
        {
          BigIntPowerIntExp(&primeMod, newExponent - currentExp, &operand5);
          BigIntRemainder(&operand1, &operand5, &operand1);
        }
        // Get coefficient of f_i.
        operand2.nbrLimbs = *ptrFi++;
        (void)memcpy(operand2.limbs, ptrFi, operand2.nbrLimbs * sizeof(int));
        ptrFi += operand2.nbrLimbs;
        BigIntMultiply(&operand1, &powerMod, &operand1);  // poly1 <- m*g
        BigIntAdd(&operand1, &operand2, &operand1);       // poly1 <- f_i + m*g
        // Store coefficient of new f_i (no spaces between coefficients).
        *ptrDest++ = operand1.nbrLimbs;
        (void)memcpy(ptrDest, &operand1.limbs, operand1.nbrLimbs*sizeof(int));
        ptrDest += operand1.nbrLimbs;    // Point to next coeffficient of f_i.
      }
      pstFactorInfo++;
    }
    computePower(newExponent);         // Compute powerMod and init Montgomery parms.
    if (newExponent == exponentMod)
    {                 // Final values of f_1, f_2,..., f_n found. Exit loop.
      break;
    }
    // Compute u(x) <- (1/m)*(f(x)/f_1(x)*a_1(x) + ... + f(x)/f_n(x)*a_n(x) - 1).
    // Init u(x) to zero.
    ptrDest = poly3;
    nbrLimbs = newNumberLength + 1;
    for (currentDegree = 0; currentDegree <= degree; currentDegree++)
    {
      *ptrDest = 1;
      *(ptrDest + 1) = 0;
      ptrDest += nbrLimbs;
    }
    pstFactorInfo = factorInfo;
    ptrNewFactor = ptrNewFactors;
    ComputeF();              // poly4 <- f (mod m^2) in Montgomery notation.
    for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
    {
      degreeFactor = pstFactorInfo->degree;
      (void)memcpy(poly1, poly4, (degree + 1) * nbrLimbs * sizeof(int));
                                        // poly2 <- f_i(x) in Montgomery notation.
      getFi(degreeFactor, *ptrNewFactor++, nbrLimbs);
      DividePolynomial(poly1, degree,   // Dividend = f.
        poly2, degreeFactor,            // Divisor = f_i.
        poly5);                         // Quotient = f / f_i.
      degreeA = getAi(nbrFactor, degreeFactor - 1);       // poly1 <- a_i.
                                        // polyMultTemp <- (f / f_i) * a_i
      MultPolynomial(degreeA, degree - degreeFactor, poly1, poly5);
      ptrDest = poly3;
      for (currentDegree = 0; currentDegree <= degree - degreeFactor + degreeA; currentDegree++)
      {
        // Get coefficient of sum.
        IntArray2BigInteger(ptrDest, &operand1);
        // Get coefficient of (f / f_i) * a_i.
        IntArray2BigInteger(&polyMultTemp[currentDegree * nbrLimbs], &operand2);
        // Add them.
        AddBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
        // Save coefficient of sum.
        BigInteger2IntArray(ptrDest, &operand1);
        ptrDest += nbrLimbs;
      }
      pstFactorInfo++;
    }
    polyToStandardNotation(poly3, degree + 1);
    // Divide sum by m.
    computePower(currentExp);
    ptrDest = poly3;
    ptrSrc = poly3;
    for (currentDegree = 0; currentDegree <= degree; currentDegree++)
    { // Get coefficient of u.
      operand1.nbrLimbs = *ptrSrc;
      (void)memcpy(operand1.limbs, ptrSrc + 1, operand1.nbrLimbs * sizeof(int));
      // Divide coefficient by m.
      BigIntDivide(&operand1, &powerMod, &operand1);
      // Store coefficient of quotient.
      *ptrDest++ = operand1.nbrLimbs;
      (void)memcpy(ptrDest, operand1.limbs, operand1.nbrLimbs * sizeof(int));
      ptrDest += oldNumberLength;
      ptrSrc += nbrLimbs;
    }
    polyToMontgomeryNotation(poly3, degree);   // poly3 <- u(x)
    nbrLimbs = NumberLength + 1;
    pstFactorInfo = factorInfo;
    // Use polyS as a temporary storage for a_i.
    ptrDest = polyS;
    for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
    {
      degreeFactor = pstFactorInfo->degree;
      // g(x) <- u(x) * a_i(x) mod f_i(x) (all computations done mod m).
      NumberLength = oldNumberLength;
      degreeA = getAi(nbrFactor, degreeFactor - 1);  // poly1 <- a_i
      MultPolynomial(degree, degreeA, poly3, poly1);
      (void)memcpy(poly4, polyMultTemp, (degree + degreeA + 1)* nbrLimbs * sizeof(int));
                                                // poly4 <- u * a_i
                                                // poly2 <- f_i(x) in Montgomery notation. 
      getFi(degreeFactor, pstFactorInfo->ptrPolyLifted, nbrLimbs);
      DividePolynomial(poly4, degree + degreeA, // Dividend = u * a_i
        poly2, degreeFactor,                    // Divisor = f_i
        NULL);                                  // Quotient: not needed.
      polyToStandardNotation(poly1, degreeA+1); // Convert a_i to standard notation.
                                                // Convert g to standard notation.
      polyToStandardNotation(poly4, degree + degreeA - degreeFactor + 1);
      // a_i(x) <- a_i(x) - m*g(x)
      for (currentDegree = 0; currentDegree <= degreeA; currentDegree++)
      {                 // Loop that computes a_i - m*g
        // Get coefficient of g.
        ptrSrc = &poly4[currentDegree * nbrLimbs];
        operand1.nbrLimbs = *ptrSrc++;
        (void)memcpy(operand1.limbs, ptrSrc, operand1.nbrLimbs * sizeof(int));
        if (!BigIntIsZero(&operand1))                     // g <- -g (mod m)
        {
          BigIntSubt(&powerMod, &operand1, &operand1);
        }
        // Get coefficient of a_i.
        ptrSrc = &poly1[currentDegree * nbrLimbs];
        operand2.nbrLimbs = *ptrSrc++;
        (void)memcpy(operand2.limbs, ptrSrc, operand2.nbrLimbs * sizeof(int));
        BigIntMultiply(&operand1, &powerMod, &operand1);  // poly1 <- m*g
        BigIntAdd(&operand1, &operand2, &operand1);       // poly1 <- a_i + m*g
        // Store coefficient of new a_i.
        *ptrDest++ = operand1.nbrLimbs;
        (void)memcpy(ptrDest, operand1.limbs, operand1.nbrLimbs * sizeof(int));
        ptrDest += operand1.nbrLimbs;  // Point to next coefficient of a_i.
      }
      pstFactorInfo++;
    }
    currentExp = newExponent;
    NumberLength = newNumberLength;
    nbrLimbs = NumberLength + 1;
    // Move a_i from polyS to polyT.
    ptrSrc = polyS;
    ptrDest = polyT;
    pstFactorInfo = factorInfo;
    for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
    {
      degreeA = pstFactorInfo -> degree - 1;
      ptrA[nbrFactor] = ptrDest;
      for (currentDegree = 0; currentDegree <= degreeA; currentDegree++)
      {
        nbrLimbs = *ptrSrc++;
        *ptrDest++ = nbrLimbs;
        (void)memcpy(ptrDest, ptrSrc, nbrLimbs * sizeof(int));
        ptrDest += nbrLimbs;
        ptrSrc += nbrLimbs;
      }
      pstFactorInfo++;
    }
    MoveFactorsAndFixPointers(factorInfo, 1);   // Output is compressed.
  }
  MoveFactorsAndFixPointers(factorInfo, compressPoly);
  return EXPR_OK;
}

void OrigPolyFromMontgomeryToStandard(void)
{
  int *ptrValue1, *ptrValue2;
  int currentDegree;
  (void)memcpy(&TestNbr, powerMod.limbs, powerMod.nbrLimbs*sizeof(limb));
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
      (void)memset(&operand2.limbs[1], 0, (NumberLength - 1)*sizeof(limb));
    }
    modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    BigInteger2IntArray(ptrValue2, &operand1);
    ptrValue1 += 1 + *ptrValue1;
    ptrValue2 += 1 + *ptrValue2;
  }
  (void)memcpy(&values[1], &poly4[0], (ptrValue2 - &poly4[0])*sizeof(int));
}

static void showPower(char **pptrOutput, int exponent)
{
  char *ptrOutput = *pptrOutput;
  if (pretty == PRETTY_PRINT)
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
  else if (pretty == TEX)
  {
    *ptrOutput++ = '^';
    *ptrOutput++ = '{';
    int2dec(&ptrOutput, exponent);
    *ptrOutput++ = '}';
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
    if (pretty == PRETTY_PRINT)
    {
      (void)strcpy(ptrOutput, lang?"<span class=\"hide\">equis </span><span aria-hidden=\"true\"><var>x</var></span>":"<var>x</var>");
      ptrOutput += strlen(ptrOutput);
    }
    else
    {
      *ptrOutput++ = 'x';
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
  int indexes[MAX_DEGREE+1];
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
        (void)strcpy(ptrOutput, pretty == PRETTY_PRINT? "&minus;": "-");
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
          if (pretty == PRETTY_PRINT)
          {
            (void)strcpy(ptrOutput, "&#8290;");
            ptrOutput += strlen(ptrOutput);
          }
          else if (pretty == PARI_GP)
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
      (void)memcpy(ptrValue2, ptrValue1, (1 + *ptrValue1) * sizeof(int));
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
    NumberLength = numLimbs(ptrValue1);
    IntArray2BigInteger(ptrValue1, &operand1);
  }
  if (operand1.sign == SIGN_NEGATIVE)
  {
    (void)strcpy(ptrOutput, " &minus;");
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
    (void)strcpy(ptrOutput, " (mod ");
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
  int isMonomial = ((polyDegree == 1) && (*pstFactorInfo->ptrPolyLifted == 1) &&
     (*(pstFactorInfo->ptrPolyLifted+1) == 0));
  if ((multiplicity > 1) && !isMonomial)
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
      (void)strcpy(ptrOutput, "&minus;");
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
    (void)strcpy(text, lang?"No se puede usar variable en el exponente":
                      "Cannot use variable in exponent");
    break;
  case EXPR_POLYNOMIAL_DIVISION_NOT_INTEGER:
    (void)strcpy(text, lang ? "La división de polinomios no es entera" :
                        "Polynomial division is not integer");
    break;
  case EXPR_DEGREE_TOO_HIGH:
    (void)strcpy(text, lang?"El grado del polinomio es muy elevado":
                      "Degree is too high");
    break;
  case EXPR_EXPONENT_TOO_LARGE:
    (void)strcpy(text, lang?"Exponente muy grande":"Exponent is too large");
    break;
  case EXPR_EXPONENT_NEGATIVE:
    (void)strcpy(text, lang?"Exponente negativo":"Exponent is negative");
    break;
  case EXPR_LEADING_COFF_MULTIPLE_OF_PRIME:
    (void)strcpy(text, lang?"El primer coeficiente es múltiplo del número primo":
      "Leading coefficient multiple of prime");
    break;
  case EXPR_CANNOT_LIFT:
    (void)strcpy(text, lang?"No se puede elevar porque hay factores duplicados":
      "Cannot lift because of duplicate factors modulo prime");
    break;
  case EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE:
    (void)strcpy(text, lang?"El módulo debe ser mayor que 1":"Modulus must be greater than one");
    break;
  case EXPR_MODULUS_MUST_BE_PRIME_EXP:
    (void)strcpy(text, lang ? "El módulo debe ser un número primo o una potencia de número primo" :
      "Modulus must be a prime number or a power of a prime");
    break;
  case EXPR_MULTIPLE_VARIABLES_NOT_ACCEPTED:
    (void)strcpy(text, lang ? "No se aceptan múltiples variables" :
      "Multiple variables are not accepted");
    break;
  default:
    textError(text, rc);
  }
  *ptrOutput++ = '<';
  *ptrOutput++ = 'p';
  *ptrOutput++ = '>';
  (void)strcpy(ptrOutput, text);
  ptrOutput += strlen(ptrOutput);
  *ptrOutput++ = '<';
  *ptrOutput++ = '/';
  *ptrOutput++ = 'p';
  *ptrOutput++ = '>';
  *ptrOutput = 0;    // Add terminator character.
}


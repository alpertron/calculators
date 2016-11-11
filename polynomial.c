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
#include "bignbr.h"
#include "highlevel.h"
#include "polynomial.h"
#define STACK_OPER_SIZE      100
#define TOKEN_NUMBER         '0'
#define TOKEN_START_EXPON    '1'
#define TOKEN_END_EXPON      '2'
#define TOKEN_UNARY_MINUS    '3'
#ifdef __EMSCRIPTEN__
void databack(char *data);
int stamp(void);
int newStamp, oldStamp;
extern char *output;
#endif
BigInteger primeMod;              // p
int exponentMod;                  // k
BigInteger powerMod;              // p^k
int degree;
extern int lang;
static BigInteger value;
BigInteger operand1, operand2, operand3, operand4, operand5;
extern limb TestNbr[MAX_LEN];
extern limb MontgomeryMultR2[MAX_LEN];
extern limb MontgomeryMultR1[MAX_LEN];
extern int NumberLength, NumberLengthR1;
static char RPNbuffer[1000000];
int values[1000000];
static int valuesIndex;
int valuesPrime[1000000];
int poly1[1000000];
int poly2[1000000];
int poly3[1000000];
static int poly4[1000000];
static int poly5[1000000];
int polyMultTemp[1000000];
static int polyLifted[1000000];
unsigned char superscripts, onlyEvaluate;
struct sFactorInfo factorInfo[MAX_DEGREE];
int nbrFactorsFound;
int *ptrOrigPoly;
int degreeOrigPoly;
enum eDivType
{
  TYPE_DIVISION,
  TYPE_MODULUS
};
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
          if (s != '+' && s != '-' && s != '*' && s != '^' && s != '/' && s != '%')
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
      else if (c == ')')
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
        prevTokenIsNumber = TRUE;
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
        {                  // Reduce mod power outside exponent.
          (void)BigIntRemainder(&value, &powerMod, &value);
          if (value.nbrLimbs < NumberLength)
          {    // Fill with zeros.
            memset(&value.limbs[value.nbrLimbs], 0, (NumberLength - value.nbrLimbs)*sizeof(int));
          }
          // Convert to Montgomery notation.
          modmult(value.limbs, MontgomeryMultR2, value.limbs);
          *ptrOutput++ = TOKEN_NUMBER;
          value.nbrLimbs = NumberLength;
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
      UncompressBigInteger(ptrValue1 + 1, &operand1);
      BigIntSubt(&powerMod, &operand1, &operand1);
      CompressBigInteger(ptrValue1 + 1, &operand1);
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
      }
      else
      {
        UncompressBigInteger(ptrValue1, &operand1);
        BigIntSubt(&powerMod, &operand1, &operand1);
        CompressBigInteger(ptrValue2, &operand1);
        ptrValue1 += *ptrValue1 + 1;
        ptrValue2 += *ptrValue2 + 1;
      }
    }
    if (ptrArgument + 1 + (ptrValue2 - &poly1[0]) > values + sizeof(values) / sizeof(values[0]))
    {
      return EXPR_OUT_OF_MEMORY;
    }
    memcpy(ptrArgument + 1, poly1, (ptrValue2 - &poly1[0])*sizeof(int));
  }
  return EXPR_OK;
}

static int AddPolynomialExpr(int *ptrArgument1, int *ptrArgument2)
{
  int *ptrValue1, currentDegree;
  int degreeMin, degreeMax;
  int degreePoly, degreeMono;
  int degree1 = *ptrArgument1;
  int degree2 = *ptrArgument2;
  if (degree1 <= 0)
  {
    if (degree1 == degree2)
    {      // Sum of two monomials of same degree.
      UncompressBigInteger(ptrArgument1, &operand1);
      UncompressBigInteger(ptrArgument2, &operand2);
      AddBigNbrModN(operand1.limbs, operand2.limbs, operand1.limbs, powerMod.limbs, powerMod.nbrLimbs);
      CompressBigInteger(ptrArgument1, &operand1);
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
        UncompressBigInteger(ptrArgument1+1, &operand2);
        UncompressBigInteger(ptrArgument2+1, &operand1);
      }
      else
      {
        degreeMax = degree2;
        degreeMin = degree1;
        UncompressBigInteger(ptrArgument1+1, &operand1);
        UncompressBigInteger(ptrArgument2+1, &operand2);
      }
      ptrValue1 = ptrArgument1;
      *ptrValue1++ = degreeMax;
      for (currentDegree = 0; currentDegree < degreeMin; currentDegree++)
      {
        *ptrValue1++ = 1;  // Number of limbs
        *ptrValue1++ = 0;  // Value = zero.
      }
      CompressBigInteger(ptrValue1, &operand1);
      ptrValue1 += 1 + *ptrValue1;
      for (currentDegree++; currentDegree < degreeMax; currentDegree++)
      {
        *ptrValue1++ = 1;  // Number of limbs
        *ptrValue1++ = 0;  // Value = zero.
      }
      CompressBigInteger(ptrValue1, &operand2);
      valuesIndex = (int)(ptrValue1 + 1 + *ptrValue1 - &values[0]);
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
      ptrPolyMin = ptrArgument2;
      ptrPolyMax = ptrArgument1;
    }
    else
    {
      degreeMin = degree1;
      degreeMax = degree2;
      ptrPolyMin = ptrArgument1;
      ptrPolyMax = ptrArgument2;
    }
    ptrValue1 = poly1;
    *ptrValue1++ = degreeMax;
    for (currentDegree = 0; currentDegree < degreeMin; currentDegree++)
    {
      UncompressBigInteger(ptrPolyMin, &operand1);
      UncompressBigInteger(ptrPolyMax, &operand2);
      AddBigNbrModN(operand1.limbs, operand2.limbs, operand1.limbs, powerMod.limbs, powerMod.nbrLimbs);
      CompressBigInteger(ptrValue1, &operand1);
      ptrPolyMin += *ptrPolyMin + 1;
      ptrPolyMax += *ptrPolyMax + 1;
      ptrValue1 += *ptrValue1 + 1;
    }
    for (; currentDegree < degreeMax; currentDegree++)
    {
      memcpy(ptrValue1, ptrPolyMax, (1 + *ptrValue1)*sizeof(int));
      ptrPolyMax += *ptrPolyMax + 1;
      ptrValue1 += *ptrValue1 + 1;
    }
    memcpy(ptrArgument1, poly1, (ptrValue1 - &poly1[0])*sizeof(int));
    valuesIndex = (int)(ptrArgument1 - &values[0] + (ptrValue1 - &poly1[0]));
    return EXPR_OK;
  }
  // Sum of polynomial and monomial.
  if (degree1 < 0)
  {
    UncompressBigInteger(ptrArgument1+1, &operand1);
    memmove(ptrArgument1 + 1, ptrArgument1 + 3 + *(ptrArgument1 + 1), (&values[valuesIndex] - ptrArgument2)*sizeof(int));
    degreeMono = -degree1;
    degreePoly = degree2;
    ptrValue1 = ptrArgument2 + 1;
  }
  else
  {
    UncompressBigInteger(ptrArgument2 + 1, &operand1);
    degreeMono = -degree2;
    degreePoly = degree1;
    ptrValue1 = ptrArgument1 + 1;
  }
  if (degreeMono > degreePoly)
  {
    int *ptrValue2 = ptrArgument1 + 1;
    *ptrArgument1 = degreeMono;
    for (currentDegree = 0; currentDegree <= degreePoly; currentDegree++)
    {
      UncompressBigInteger(ptrValue1, &operand2);
      CompressBigInteger(ptrValue2, &operand2);
      ptrValue1 += 1 + *ptrValue1;
      ptrValue2 += 1 + *ptrValue2;
    }
    for (; currentDegree < degreeMono; currentDegree++)
    {
      *ptrValue1++ = 1;  // Number of limbs
      *ptrValue1++ = 0;  // Value = zero.
    }
    CompressBigInteger(ptrValue1, &operand1);
    valuesIndex = (int)(ptrValue1 + 1 + *ptrValue1 - &values[0]);
  }
  else
  {        // Degree of polynomial greater than degree of monomial.
    int differenceOfDegrees;

    *ptrArgument1 = degreePoly;
    ptrValue1 = ptrArgument1 + 1;
    for (currentDegree = 0; currentDegree < degreeMono; currentDegree++)
    {
      ptrValue1 += 1 + *ptrValue1;
    }
    UncompressBigInteger(ptrValue1, &operand2);
    AddBigNbrModN(operand1.limbs, operand2.limbs, operand1.limbs, powerMod.limbs, powerMod.nbrLimbs);
    CompressBigInteger(poly1, &operand1);
    differenceOfDegrees = poly1[0] - *ptrValue1;
    if (differenceOfDegrees > 0)
    {    // New coefficient is greater than old coefficient.
      memmove(ptrValue1 + differenceOfDegrees, ptrValue1, (ptrArgument2 - ptrValue1)*sizeof(int));
    }
    else if (differenceOfDegrees < 0)
    {
      memmove(ptrValue1, ptrValue1 - differenceOfDegrees, (ptrArgument2 - ptrValue1 - differenceOfDegrees)*sizeof(int));
    }
    memcpy(ptrValue1, poly1, (1 + poly1[0])*sizeof(int));
    valuesIndex += differenceOfDegrees;
  }
  return EXPR_OK;
}

// Multiply poly1 by poly2. The result will be stored in polyMultTemp.
static void MultPolynomialInternal(int degree1, int degree2, /*@in@*/int *factor1, /*@in@*/int *factor2)
{
  int currentDegree, currentDegree2;
  int *ptrValue1 = polyMultTemp;
  int nbrLimbs = powerMod.nbrLimbs + 1;
  for (currentDegree = 0; currentDegree <= degree1 + degree2; currentDegree++)
  {
    *ptrValue1 = 1;        // Initialize number to zero.
    *(ptrValue1+1) = 0;
    ptrValue1 += nbrLimbs;
  }
  ptrValue1 = factor1;
  for (currentDegree = 0; currentDegree <= degree1; currentDegree++)
  {
    int *ptrValue2, *ptrValueProd;
    UncompressBigInteger(ptrValue1, &operand1);
    ptrValue2 = factor2;
    ptrValueProd = &polyMultTemp[currentDegree*nbrLimbs];
    for (currentDegree2 = 0; currentDegree2 <= degree2; currentDegree2++)
    {
      UncompressBigInteger(ptrValue2, &operand2);
      modmult(operand1.limbs, operand2.limbs, operand2.limbs);
      UncompressBigInteger(ptrValueProd, &operand3);
      AddBigNbrMod(operand2.limbs, operand3.limbs, operand3.limbs);
      CompressBigInteger(ptrValueProd, &operand3);
      ptrValue2 += nbrLimbs;
      ptrValueProd += nbrLimbs;
    }
    ptrValue1 += nbrLimbs;
  }
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
    UncompressBigInteger(ptrArgument1 + 1, &operand1);
    UncompressBigInteger(ptrArgument2 + 1, &operand2);
    modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    CompressBigInteger(ptrArgument1 + 1, &operand1);
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
    // Copy first factor to poly1
    ptrValue1 = ptrArgument1 + 1;
    ptrValue2 = poly1;
    for (currentDegree = 0; currentDegree <= degree1; currentDegree++)
    {
      memcpy(ptrValue2, ptrValue1, (1 + *ptrValue1)*sizeof(int));
      ptrValue1 += 1 + *ptrValue1;
      ptrValue2 += nbrLimbs;
    }
    // Copy second factor to poly2
    ptrValue1 = ptrArgument2 + 1;
    ptrValue2 = poly2;
    for (currentDegree = 0; currentDegree <= degree2; currentDegree++)
    {
      memcpy(ptrValue2, ptrValue1, (1 + *ptrValue1)*sizeof(int));
      ptrValue1 += 1 + *ptrValue1;
      ptrValue2 += nbrLimbs;
    }
    // Perform multiplication storing the result in polyMultTemp
    MultPolynomialInternal(degree1, degree2, poly1, poly2);
    // Move product back to values stack.
    ptrValue1 = ptrArgument1 + 1;
    ptrValue2 = polyMultTemp;
    for (currentDegree = 0; currentDegree <= degree1+degree2; currentDegree++)
    {
      memcpy(ptrValue1, ptrValue2, (1 + *ptrValue2)*sizeof(int));
      ptrValue1 += 1 + *ptrValue1;
      ptrValue2 += nbrLimbs;
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
    UncompressBigInteger(ptrArgument1 + 1, &operand1);
  }
  else
  {      // Second factor is monomial.
    degreeMono = -degree2;
    degreePoly = degree1;
         // Point to first coefficient of polynomial.
    ptrValue1 = ptrArgument1 + 1;
         // Get coefficient of monomial.
    UncompressBigInteger(ptrArgument2 + 1, &operand1);
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
    UncompressBigInteger(ptrValue1, &operand2);
    modmult(operand1.limbs, operand2.limbs, operand2.limbs);
    CompressBigInteger(ptrValue2, &operand2);
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
  valuesIndex = (int)(ptrValue1 - &values[0]);
  return EXPR_OK;
}

static void ToPoly(int degree, int *polySrc, int *polyDest)
{
  int currentDegree;
  polySrc++;
  if (degree < 0)
  {    // Polynomial is a monomial
    for (currentDegree = 0; currentDegree < -degree; currentDegree++)
    {
      *polyDest = 1;
      *(polyDest + 1) = 0;
      polyDest += NumberLength + 1;
    }
    memcpy(polyDest, polySrc, (*(polySrc)+1)*sizeof(int));
  }
  else
  {   // Polynomial
    for (currentDegree = 0; currentDegree <= degree; currentDegree++)
    {
      int nbrLimbs = *(polySrc)+1;
      memcpy(polyDest, polySrc, nbrLimbs*sizeof(int));
      polyDest += NumberLength + 1;
      polySrc += nbrLimbs;
    }
  }
}

static void FromPoly(int degree, int *polyDest, int *polySrc)
{
  int currentDegree;
  *polyDest++ = degree;
  for (currentDegree = 0; currentDegree <= degree; currentDegree++)
  {
    int nbrLimbs = *(polySrc)+1;
    memcpy(polyDest, polySrc, nbrLimbs*sizeof(int));
    polyDest += NumberLength + 1;
    polySrc += nbrLimbs;
  }
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
    UncompressBigInteger(ptrArgument1 + 1, &operand1);
    UncompressBigInteger(ptrArgument2 + 1, &operand2);
    ModInvBigNbr(operand2.limbs, operand2.limbs, TestNbr, NumberLength);
    modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    CompressBigInteger(ptrArgument1 + 1, &operand1);
    valuesIndex = (int)(ptrArgument1 + 2 + *(ptrArgument1 + 1) - &values[0]);
    return EXPR_OK;
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
    if (degreeBase*expon > MAX_DEGREE)
    {
      return EXPR_DEGREE_TOO_HIGH;
    }
    *ptrArgument1 = degreeBase*expon;
    UncompressBigInteger(ptrArgument1+1, &operand1);
    if (expon < 0)
    {
      return EXPR_EXPONENT_NEGATIVE;
    }
    exponLimb.x = expon;
    modPowLimb(operand1.limbs, &exponLimb, operand2.limbs);
    CompressBigInteger(ptrArgument1 + 1, &operand2);
    valuesIndex = (int)(ptrArgument1+2+*(ptrArgument1+1) - &values[0]);
    return EXPR_OK;
  }
                // Polynomial.
  // Copy base to poly1
  ptrValue1 = ptrArgument1 + 1;
  ptrValue2 = poly1;
  for (currentDegree = 0; currentDegree <= degreeBase; currentDegree++)
  {
    memcpy(ptrValue2, ptrValue1, (1 + *ptrValue1)*sizeof(int));
    ptrValue1 += 1 + *ptrValue1;
    ptrValue2 += nbrLimbs;
  }
  SetNumberToOne(&poly2[0]); // Initialize power with polynomial 1.
  degreePower = 0;
  for (mask = 1 << (BITS_PER_GROUP - 1); mask > 0; mask >>= 1)
  {
    MultPolynomialInternal(degreePower, degreePower, poly2, poly2);
    degreePower <<= 1;
    memcpy(poly2, polyMultTemp, (degreePower+1)*nbrLimbs*sizeof(int));
    if ((expon & mask) != 0)
    {
      MultPolynomialInternal(degreeBase, degreePower, poly1, poly2);
      degreePower += degreeBase;
      memcpy(poly2, polyMultTemp, (degreePower+1)*nbrLimbs*sizeof(int));
    }
  }
  // Move power back to values stack.
  *ptrArgument1 = degreePower;
  ptrValue1 = ptrArgument1 + 1;
  ptrValue2 = polyMultTemp;
  for (currentDegree = 0; currentDegree <= degreePower; currentDegree++)
  {
    int len = 1 + *ptrValue2;
    memcpy(ptrValue1, ptrValue2, len*sizeof(int));
    ptrValue1 += len;
    ptrValue2 += nbrLimbs;
  }
  valuesIndex = (int)(ptrValue1 - &values[0]);
  return EXPR_OK;
}

static void computePower(int expo)
{
  operand1.limbs[0].x = expo & MAX_VALUE_LIMB;
  operand1.limbs[1].x = expo >> BITS_PER_GROUP;
  if (operand1.limbs[1].x == 0)
  {
    operand1.nbrLimbs = 1;
  }
  else
  {
    operand1.nbrLimbs = 2;
  }
  BigIntPower(&primeMod, &operand1, &powerMod);
  memcpy(&TestNbr, &powerMod, powerMod.nbrLimbs*sizeof(limb));
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
        *ptrValue1 = -val;
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
    case '+':
      ptrValue2 = stackValues[--stackIndex];
      ptrValue1 = stackValues[stackIndex - 1];
      if (insideExpon != 0)
      {
        *ptrValue1 += *ptrValue2;
        if (*ptrValue1 < -MAX_VALUE_LIMB || *ptrValue1 > MAX_VALUE_LIMB)
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
        if (*ptrValue1 < -MAX_VALUE_LIMB || *ptrValue1 > MAX_VALUE_LIMB)
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
        if (*ptrValue1 < -MAX_VALUE_LIMB || *ptrValue1 > MAX_VALUE_LIMB)
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
        for (count = expon; count > 0; count--)
        {
          pwr *= val;
          if (pwr < -MAX_VALUE_LIMB || pwr > MAX_VALUE_LIMB)
          {
            return EXPR_EXPONENT_TOO_LARGE;
          }
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
    UncompressBigInteger(&values[1], &operand1);
    degree = -values[0];
    values[0] = degree;
    ptrValue1 = &values[1];
    for (currentDegree = 0; currentDegree < degree; currentDegree++)
    {
      *ptrValue1++ = 1;   // Initialize coefficient to zero.
      *ptrValue1++ = 0;
    }
    CompressBigInteger(ptrValue1, &operand1);
  }
  return EXPR_OK;
}

void ConvertToMonic(int *poly, int polyDegree)
{
  int nbrLimbs = NumberLength + 1;
  int currentDegree;
  UncompressBigInteger(poly+polyDegree*nbrLimbs, &operand1);
  ModInvBigNbr(operand1.limbs, operand1.limbs, TestNbr, NumberLength);
  for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {
    UncompressBigInteger(poly+currentDegree*nbrLimbs, &operand2);
    modmult(operand1.limbs, operand2.limbs, operand2.limbs);
    CompressBigInteger(poly+currentDegree*nbrLimbs, &operand2);
  }
}

// arg1 must not be changed by this routine.
// arg2 can be changed by this routine.
void PolynomialGcd(int *arg1, int degree1, int *arg2, int degree2, int *gcd, int *degreeGcd)
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
    {          // Get remainder of long division.
      UncompressBigInteger(ptrArgMax + currentDegree*nbrLimbs, &operand1);
      ptrTemp = ptrArgMax + (currentDegree - degreeMin)*nbrLimbs;
      for (index = 0; index <= degreeMin; index++)
      {
        UncompressBigInteger(ptrArgMin + index*nbrLimbs, &operand2);
        modmult(operand1.limbs, operand2.limbs, operand2.limbs);
        UncompressBigInteger(ptrTemp, &operand3);
        SubtBigNbrMod(operand3.limbs, operand2.limbs, operand3.limbs);
        CompressBigInteger(ptrTemp, &operand3);
        ptrTemp += nbrLimbs;
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
        memcpy(gcd, ptrArgMin, degreeMin*nbrLimbs*sizeof(int));
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

// Convert from Montgomery notation to standard notation by multiplying by 1.
void ToStandardNotation(int *nbr, int qtyNbrs)
{
  int currentNbr;
  for (currentNbr = 0; currentNbr < qtyNbrs; currentNbr++)
  {
    UncompressBigInteger(nbr, &operand1);
    operand2.nbrLimbs = NumberLength;
    operand2.limbs[0].x = 1;
    if (NumberLength > 1)
    {
      memset(&operand2.limbs[1], 0, (NumberLength - 1)*sizeof(limb));
    }
    modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    CompressBigInteger(nbr, &operand1);
    nbr += NumberLength + 1;
  }
}

static void ToMontgomeryNotation(int *nbr, int qtyNbrs)
{
  int currentNbr;

  for (currentNbr = qtyNbrs - 1; currentNbr >= 0; currentNbr--)
  {
    UncompressBigInteger(nbr, &operand1);
    modmult(operand1.limbs, MontgomeryMultR2, operand1.limbs);
    CompressBigInteger(nbr, &operand1);
    nbr += NumberLength + 1;
  }
}

// In this routine, the dividend is replaced by the remainder of the division.
void DividePolynomial(/*@in@*/int *pDividend, int dividendDegree, /*@in@*/int *pDivisor, int divisorDegree, /*@out@*/int *ptrQuotient)
{
  int currentDegree, index;
  int nbrLimbs = NumberLength + 1;
  int divisorIsOne;
  int *ptrQuot;
  if (divisorDegree > dividendDegree)
  {    // Quotient is zero.
    *ptrQuotient = 1;
    *(ptrQuotient+1) = 0;
    return;
  }
  UncompressBigInteger(pDivisor + divisorDegree*nbrLimbs, &operand1);
  memcpy(operand5.limbs, operand1.limbs, NumberLength*sizeof(int));
  divisorIsOne = !memcmp(operand1.limbs, MontgomeryMultR1, NumberLength*sizeof(int));
  if (!divisorIsOne)
  {        // Leading coefficient is not 1.
    ConvertToMonic(pDivisor, divisorDegree);
    // operand1 holds the inverse of the leading coefficient of divisor.
    // Multiply dividend by this number.
    for (currentDegree = 0; currentDegree <= dividendDegree; currentDegree++)
    {
      UncompressBigInteger(pDividend + currentDegree*nbrLimbs, &operand2);
      modmult(operand1.limbs, operand2.limbs, operand2.limbs);
      CompressBigInteger(pDividend + currentDegree*nbrLimbs, &operand2);
    }
  }
  ptrQuot = ptrQuotient + (dividendDegree - divisorDegree)*nbrLimbs;
  for (currentDegree = dividendDegree; currentDegree >= divisorDegree; currentDegree--)
  {
    int *ptrDivisor;
    int *ptrDividend = pDividend + currentDegree*nbrLimbs;
    UncompressBigInteger(ptrDividend, &operand1);
    CompressBigInteger(ptrQuot, &operand1);  // Store coefficient of quotient.
    ptrDivisor = pDivisor + divisorDegree*nbrLimbs;
    for (index = 0; index <= divisorDegree; index++)
    {
      UncompressBigInteger(ptrDivisor, &operand2);
      modmult(operand1.limbs, operand2.limbs, operand2.limbs);
      UncompressBigInteger(ptrDividend, &operand3);
      SubtBigNbrMod(operand3.limbs, operand2.limbs, operand3.limbs);
      CompressBigInteger(ptrDividend, &operand3);
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
      UncompressBigInteger(pDividend + currentDegree*nbrLimbs, &operand2);
      modmult(operand5.limbs, operand2.limbs, operand2.limbs);
      CompressBigInteger(pDividend + currentDegree*nbrLimbs, &operand2);
    }
           // Restore divisor.
    for (currentDegree = 0; currentDegree <= divisorDegree; currentDegree++)
    {
      UncompressBigInteger(pDivisor + currentDegree*nbrLimbs, &operand2);
      modmult(operand5.limbs, operand2.limbs, operand2.limbs);
      CompressBigInteger(pDivisor + currentDegree*nbrLimbs, &operand2);
    }
  }
}

// Multiply two polynomials mod polyMod.
void multPolynomial(/*@in@*/int *polyFact1, /*@in@*/int *polyFact2, /*@out@*/int *polyProduct, int polyDegree, /*@in@*/int *polyMod)
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
    UncompressBigInteger(ptrPoly1, &operand1);
    UncompressBigInteger(ptrPoly2, &operand2);
    modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    for (index2 = polyDegree - 2; index2 >= index1; index2--)
    {
      ptrPoly2 -= nbrLimbs;
      ptrPoly1 += nbrLimbs;
      UncompressBigInteger(ptrPoly1, &operand2);
      UncompressBigInteger(ptrPoly2, &operand3);
      modmult(operand2.limbs, operand3.limbs, operand2.limbs);
      AddBigNbrMod(operand2.limbs, operand1.limbs, operand1.limbs);
    }
    CompressBigInteger(ptrPolyTemp, &operand1);
    ptrPolyTemp -= nbrLimbs;
  }
  // Get remainder of long division by polyMod and append next limbs of the product.
  for (index1 = polyDegree-2; index1 >= 0; index1--)
  {
    ptrPoly2 = &polyMultTemp[(polyDegree - 1)*nbrLimbs];
    // Back up leading coefficient.
    memcpy(ptrPoly2 + nbrLimbs, ptrPoly2, nbrLimbs*sizeof(int));
    UncompressBigInteger(ptrPoly2, &operand3);
    ptrPoly1 = polyMod+polyDegree*nbrLimbs;
    for (index2 = polyDegree - 2; index2 >= 0; index2--)
    {
      ptrPoly1 -= nbrLimbs;
      ptrPoly2 -= nbrLimbs;
      UncompressBigInteger(ptrPoly1, &operand1);
      modmult(operand3.limbs, operand1.limbs, operand1.limbs);
      UncompressBigInteger(ptrPoly2, &operand2);
      SubtBigNbrMod(operand2.limbs, operand1.limbs, operand1.limbs);
      CompressBigInteger(ptrPoly2+nbrLimbs, &operand1);
    }
    ptrPoly1 = polyFact1 + index1*nbrLimbs;
    ptrPoly2 = polyFact2;
    UncompressBigInteger(ptrPoly1, &operand1);
    UncompressBigInteger(ptrPoly2, &operand2);
    modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    for (index2 = 1; index2 <= index1; index2++)
    {
      ptrPoly1 -= nbrLimbs;
      ptrPoly2 += nbrLimbs;
      UncompressBigInteger(ptrPoly1, &operand2);
      UncompressBigInteger(ptrPoly2, &operand3);
      modmult(operand2.limbs, operand3.limbs, operand2.limbs);
      AddBigNbrMod(operand2.limbs, operand1.limbs, operand1.limbs);
    }
    UncompressBigInteger(polyMod, &operand2);
    UncompressBigInteger(&polyMultTemp[polyDegree*nbrLimbs], &operand3);
    modmult(operand3.limbs, operand2.limbs, operand2.limbs);
    SubtBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
    CompressBigInteger(polyMultTemp, &operand1);
  }
  memcpy(polyProduct, polyMultTemp, polyDegree*nbrLimbs*sizeof(int));
}
  // Perform polyPower <- polyBase ^ expon (mod polyMod)
void powerPolynomial(int *polyBase, int *polyMod, int polyDegree, BigInteger *expon, int *polyPower)
{
  int mask, index;
  int nbrLimbs = NumberLength + 1;
  int powerIsOne = TRUE;
  // Initialize polyPower to 1.
  memset(polyPower, 0, (polyDegree + 1)*nbrLimbs*sizeof(int));
  for (index = 0; index <= polyDegree; index++)
  {
    *(polyPower+index*nbrLimbs) = 1;
  }
  SetNumberToOne(polyPower);
  for (index = expon->nbrLimbs - 1; index >= 0; index--)
  {
    int groupExp = (int)(expon->limbs[index].x);
    for (mask = 1 << (BITS_PER_GROUP - 1); mask > 0; mask >>= 1)
    {
      if (!powerIsOne)
      {
        multPolynomial(polyPower, polyPower, polyPower, polyDegree, polyMod);
      }
      if ((groupExp & mask) != 0)
      {
        multPolynomial(polyPower, polyBase, polyPower, polyDegree, polyMod);
        powerIsOne = FALSE;
      }
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
    UncompressBigInteger(ptrValue1, &operand1);
    modmultInt(operand1.limbs, currentDegree, operand1.limbs);
    CompressBigInteger(ptrValue2, &operand1);
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

    PolynomialGcd(poly, polyDegree, poly1, currentDegree, poly2, &degreeC); // poly2 = c
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
      PolynomialGcd(poly2, degreeC, poly1, degreeW, poly3, &degreeY);       // poly3 = y
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

static int *ComputeUorV(/*@in@*/int *ptrQuotients[], int nbrQuotients, int startR, /*@out@*/int *ptrP1, /*@out@*/int *ptrP2, /*@out@*/int *pDegreeR)
{
  int *tmpPtr;
  int *ptrR = ptrP1;
  int *ptrOldR = ptrP2;
  int degreeR = 0;
  int degreeOldR = 0;
  int counter, currentDegree;
  int nbrLimbs = NumberLength + 1;
  int zero[2] = { 1, 0 };
  if (startR != 0)
  {    // R <- inverse of gcd, oldR <- 0
    CompressBigInteger(ptrR, &operand4);
    *ptrOldR = 1;
    *(ptrOldR + 1) = 0;
  }
  else
  {    // R <- 0, oldR <- inverse of gcd.
    CompressBigInteger(ptrOldR, &operand4);
    *ptrR = 1;
    *(ptrR + 1) = 0;
  }
  for (counter = 0; counter < nbrQuotients; counter++)
  {
    int tmpDegree, offset;
    int degreeQ = (int)((ptrQuotients[counter + 1] - ptrQuotients[counter])/nbrLimbs - 1);
    // Multiply ptrR by ptrQuotients[counter]. The result will be stored in polyMultTemp.
    MultPolynomialInternal(degreeR, degreeQ, ptrR, ptrQuotients[counter]);
    // Compute ptrOldR <- ptrOldR - polyMultTemp.
    degreeQ += degreeR; // Degree of product (always greater than degree of old R).
    offset = 0;
    for (currentDegree = 0; currentDegree <= degreeOldR; currentDegree++)
    {
      UncompressBigInteger(ptrOldR + offset, &operand1);
      UncompressBigInteger(&polyMultTemp[offset], &operand2);
      SubtBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
      CompressBigInteger(ptrOldR + offset, &operand1);
      offset += nbrLimbs;
    }
    for (; currentDegree <= degreeQ; currentDegree++)
    {
      UncompressBigInteger(zero, &operand1);
      UncompressBigInteger(&polyMultTemp[offset], &operand2);
      SubtBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
      CompressBigInteger(ptrOldR + offset, &operand1);
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
// Given coprime polynomials A and B, find polynomials U and V such that A*U + B*V = 1.
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
    ptrQuotients[nbrQuotients++] = ptrQuotient;
    DividePolynomial(ptrOldR, degreeOldR, ptrR, degreeR, ptrQuotient);
    ptrQuotient += (degreeOldR - degreeR + 1)*nbrLimbs;
    degreeOldR = degreeR;
    degreeR = getDegreePoly(ptrOldR, degreeR-1);
    // Exchange pointers.
    tmpPtr = ptrR;
    ptrR = ptrOldR;
    ptrOldR = tmpPtr;
  }
  ptrQuotients[nbrQuotients] = ptrQuotient;
  UncompressBigInteger(ptrR, &operand4);     // Get gcd (which is a constant).
  ModInvBigNbr(operand4.limbs, operand4.limbs, primeMod.limbs, primeMod.nbrLimbs); // Invert it.
  ptrR = ComputeUorV(ptrQuotients, nbrQuotients, 1, ptrP1, ptrP2, &degreeFirst);
  memcpy((polyExchanged? ptrU: ptrV), ptrR, (degreeFirst+1)*nbrLimbs*sizeof(int));
  ptrR = ComputeUorV(ptrQuotients, nbrQuotients, 0, ptrP1, ptrP2, &degreeSecond);
  memcpy((polyExchanged? ptrV : ptrU), ptrR, (degreeSecond + 1)*nbrLimbs*sizeof(int));
  *pDegreeV = (polyExchanged ? degreeSecond : degreeFirst);
  *pDegreeU = (polyExchanged ? degreeFirst : degreeSecond);
}

static void MoveNumbers(int *poly, int qtyNbrs, int oldNbrLen, int newNbrLen)
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
// s' <- s - r
// t' <- t - t*e - q*g'

// Memory usage:
// poly1: s', g'
// poly2: t', h'
// poly3: e
// poly4: q, r, f
// poly5: s*e, t*e, t*e + q*g'

static void Adjust(int *polyFirst, int degreeFirst,
  int *polySecond,
  int *polyH, int degreeH,
  int *polyS, int degreeS,
  int *polyT, int degreeT,
  int *polyG, int degreeG, int degreeE, int subtract)
{
  int nbrLimbs = NumberLength + 1;
  int degreeQ = degreeS + degreeE - degreeH;
  int *polyR;
  int currentDegree;
  if (degreeE == 0 && poly3[0] == 1 && poly3[1] == 0)
  {
    return;     // e = 0 so there is nothing to adjust.
  }
  if (degreeQ < 0)
  {
    degreeQ = 0;
  }
  polyR = &poly4[(degreeQ + 1)*nbrLimbs];
  // Compute q, r such that s*e = q*h + r
  MultPolynomialInternal(degreeS, degreeE, polyS, poly3);
  DividePolynomial(polyMultTemp, degreeS + degreeE, polyH, degreeH, poly4);
  // At this moment poly4 = Q.
  memcpy(polyR, polyMultTemp, (degreeH+1)*nbrLimbs*sizeof(int));
  // Compute t*e + q*g and then add or subtract it from polyFirst.
  MultPolynomialInternal(degreeT, degreeE, polyT, poly3);
  memcpy(poly5, polyMultTemp, (degreeT + degreeE + 1)*nbrLimbs*sizeof(int));
  MultPolynomialInternal(degreeQ, degreeG, poly4, polyG);
  for (currentDegree = 0; currentDegree <= degreeFirst; currentDegree++)
  {
    UncompressBigInteger(&poly5[currentDegree*nbrLimbs], &operand1);
    UncompressBigInteger(&polyMultTemp[currentDegree*nbrLimbs], &operand2);
    AddBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
    UncompressBigInteger(polyFirst + currentDegree*nbrLimbs, &operand2);
    if (subtract != 0)
    {
      SubtBigNbrMod(operand2.limbs, operand1.limbs, operand1.limbs);
    }
    else
    {
      AddBigNbrMod(operand2.limbs, operand1.limbs, operand1.limbs);
    }
    CompressBigInteger(polyFirst + currentDegree*nbrLimbs, &operand1);
  }
  for (currentDegree = 0; currentDegree < degreeH; currentDegree++)
  {
    UncompressBigInteger(polySecond + currentDegree*nbrLimbs, &operand2);
    UncompressBigInteger(polyR+currentDegree*nbrLimbs, &operand1);
    if (subtract != 0)
    {
      SubtBigNbrMod(operand2.limbs, operand1.limbs, operand1.limbs);
    }
    else
    {
      AddBigNbrMod(operand2.limbs, operand1.limbs, operand1.limbs);
    }
    CompressBigInteger(polySecond + currentDegree*nbrLimbs, &operand1);
  }
}

int HenselLifting(void)
{
  int nbrFactor, nbrFactor2, currentExp, oldNumberLength;
  int degreeFactor, degreeS, degreeT, degreeG, currentDegree;
  int index, nbrLimbs, *ptrValue1, *polyS, *polyT;
  int *ptrPolyLifted = polyLifted;
  struct sFactorInfo *pstFactorInfo, *pstFactorInfo2;
  if (exponentMod == 1)
  {
    pstFactorInfo = factorInfo;
    for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
    {
      pstFactorInfo->ptrPolyLifted = pstFactorInfo->ptr;
      pstFactorInfo++;
    }
    return EXPR_OK;   // No lift has to be done.
  }
  // Hensel lifting only works when there are no repeated factors mod p.
  pstFactorInfo = factorInfo;
  for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    if (pstFactorInfo->multiplicity != 1)
    {   // Repeated factor found.
      return EXPR_CANNOT_LIFT;
    }
    pstFactorInfo++;
  }
  pstFactorInfo = factorInfo;
  for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    if (exponentMod == 1)
    {
      pstFactorInfo->ptrPolyLifted = pstFactorInfo->ptr;
      continue;
    }
    powerMod.nbrLimbs = primeMod.nbrLimbs;
    memcpy(powerMod.limbs, primeMod.limbs, primeMod.nbrLimbs*sizeof(limb));
    memcpy(&TestNbr, powerMod.limbs, primeMod.nbrLimbs*sizeof(limb));
    GetMontgomeryParms(powerMod.nbrLimbs);
    nbrLimbs = NumberLength + 1;
    pstFactorInfo->ptrPolyLifted = ptrPolyLifted;
    SetNumberToOne(poly1);
    degreeFactor = 0;
    pstFactorInfo2 = factorInfo;
    for (nbrFactor2 = 0; nbrFactor2 < nbrFactorsFound; nbrFactor2++)
    {    // Loop that multiplies all factors except the current one. poly1 holds the result.
      if (nbrFactor != nbrFactor2)
      {
        memcpy(poly2, pstFactorInfo2->ptr, pstFactorInfo2->degree*nbrLimbs*sizeof(int));
        SetNumberToOne(&poly2[pstFactorInfo2->degree*nbrLimbs]);
        MultPolynomialInternal(degreeFactor, pstFactorInfo2->degree, poly1, poly2);
        degreeFactor += pstFactorInfo2->degree;
        memcpy(poly1, polyMultTemp, (degreeFactor + 1)*nbrLimbs*sizeof(int));
      }
      pstFactorInfo2++;
    }
    memcpy(poly2, pstFactorInfo->ptr, pstFactorInfo->degree*nbrLimbs*sizeof(int));
    SetNumberToOne(&poly2[pstFactorInfo->degree*nbrLimbs]);
    degreeFactor = pstFactorInfo->degree;
    // At this moment g = poly1, h = poly2.
    degreeG = degree - degreeFactor;
    polyS = &poly1[(degreeG + 1)*nbrLimbs];
    polyT = &poly2[(degreeFactor + 1)*nbrLimbs];
    ExtendedGcdPolynomial(poly1, degree - degreeFactor, poly2, degreeFactor, poly3,
      poly4, poly5, polyS, &degreeS, polyT, &degreeT);
    // At this moment, s = poly1 (after g), t = poly2 (after h).
    currentExp = 1;
    pstFactorInfo++;
    // Loop that performs the lifting.
    while (currentExp != exponentMod)
    {
#ifdef __EMSCRIPTEN__
      newStamp = stamp();
      if (newStamp != oldStamp)
      {
        char *ptrOutput = output;
        oldStamp = newStamp;
        if (lang)
        {
          strcpy(ptrOutput, "1<p>Aplicando lema de Hensel en el factor número ");
          ptrOutput += strlen(ptrOutput);
          int2dec(&ptrOutput, nbrFactor + 1);
          strcpy(ptrOutput, " de ");
          ptrOutput += strlen(ptrOutput);
          int2dec(&ptrOutput, nbrFactorsFound);
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
          strcpy(ptrOutput, " processing exponent ");
          ptrOutput += strlen(ptrOutput);
          int2dec(&ptrOutput, currentExp);
          strcpy(ptrOutput, " of ");
        }
        ptrOutput += strlen(ptrOutput);
        int2dec(&ptrOutput, exponentMod);
        strcpy(ptrOutput, ".</p>");
        databack(output);
      }
#endif
        // Convert g, h, s, t from Montgomery notation to standard notation by multiplying by 1 in Montgomery notation.
      ToStandardNotation(poly1, degreeG+1);                  // Convert g to standard notation.
      ToStandardNotation(poly2, degreeFactor+1);             // Convert h to standard notation.
      ToStandardNotation(polyS, degreeS+1);                  // Convert s to standard notation.
      ToStandardNotation(polyT, degreeT+1);                  // Convert t to standard notation.
      oldNumberLength = NumberLength;
      currentExp <<= 1;
      if (currentExp > exponentMod)
      {
        currentExp = exponentMod;
      }
      computePower(currentExp);
      // Compute polynomial poly4 <- value mod prime^currentExp.
      degree = values[0];
      ptrValue1 = &values[1];
      for (currentDegree = 0; currentDegree <= degree; currentDegree++)
      {
        UncompressBigInteger(ptrValue1, &operand1);
        (void)BigIntRemainder(&operand1, &powerMod, &operand1);
        // Convert operand1 from standard to Montgomery notation.
        modmult(operand1.limbs, MontgomeryMultR2, operand1.limbs);
        CompressBigInteger(&poly4[currentDegree*(NumberLength + 1)], &operand1);
        ptrValue1 += 1 + *ptrValue1;
      }
      // Convert polynomial mod prime to monic (leading coefficient must be 1).
      ConvertToMonic(poly4, degree);
      // Now the numbers are greater. Move them so they can fit. Old size = oldNumberLength, new size = NumberLength
      MoveNumbers(poly1, degree - degreeFactor + degreeS + 2, oldNumberLength+1, NumberLength+1);
      MoveNumbers(poly2, degreeFactor + degreeT + 2, oldNumberLength+1, NumberLength+1);
      nbrLimbs = NumberLength + 1;
      polyS = &poly1[(degree - degreeFactor + 1)*nbrLimbs];
      polyT = &poly2[(degreeFactor + 1)*nbrLimbs];
      // Convert standard notation to Montgomery notation.
      ToMontgomeryNotation(poly1, degreeG+1);                // Convert g to Montgomery notation.
      ToMontgomeryNotation(poly2, degreeFactor+1);           // Convert h to Montgomery notation.
      ToMontgomeryNotation(polyS, degreeS+1);                // Convert s to Montgomery notation.
      ToMontgomeryNotation(polyT, degreeT+1);                // Convert t to Montgomery notation.
      // Compute e <- f - g*h and store it into poly3.
      MultPolynomialInternal(degree - degreeFactor, degreeFactor, poly1, poly2);
      for (index = 0; index <= degree; index++)
      {
        UncompressBigInteger(&poly4[index*nbrLimbs], &operand1);
        UncompressBigInteger(&polyMultTemp[index*nbrLimbs], &operand2);
        SubtBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
        CompressBigInteger(&poly3[index*nbrLimbs], &operand1);
      }
      Adjust(poly1, degreeG,                             // g
        poly2,                                           // h
        poly2, degreeFactor,                             // h
        polyS, degreeS,                                  // s
        polyT, degreeT,                                  // t
        poly1, degreeG,                                  // g
        getDegreePoly(poly3, degree),                    // degree of e
        0);
      if (currentExp == exponentMod)
      {
        break;
      }
      // Compute e <- s*g' + t*h' - 1 and store it into poly3.
      MultPolynomialInternal(degreeS, degreeG, polyS, poly1);
      memcpy(poly3, polyMultTemp, (degreeS + degreeG + 1)*nbrLimbs*sizeof(int));
      MultPolynomialInternal(degreeT, degreeFactor, polyT, poly2);
      for (currentDegree = 0; currentDegree <= degreeS+degreeG; currentDegree++)
      {
        UncompressBigInteger(&poly3[currentDegree*nbrLimbs], &operand1);
        UncompressBigInteger(&polyMultTemp[currentDegree*nbrLimbs], &operand2);
        AddBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
        CompressBigInteger(&poly3[currentDegree*nbrLimbs], &operand1);
      }
      UncompressBigInteger(poly3, &operand1);
      SubtBigNbrMod(operand1.limbs, MontgomeryMultR1, operand1.limbs);
      CompressBigInteger(poly3, &operand1);                               // Store e.
      Adjust(polyT, degreeT,             // t
        polyS,                           // s
        poly2, degreeFactor,             // h
        polyS, degreeS,                  // s
        polyT, degreeT,                  // t
        poly1, degreeG,                  // g
        getDegreePoly(poly3, degreeS+degreeG),  // degree of e
        1);
    }
    memcpy(ptrPolyLifted, poly2, degreeFactor*nbrLimbs*sizeof(int));
    ptrPolyLifted += degreeFactor*nbrLimbs;
  }
  return EXPR_OK;
}

void OrigPolyFromMontgomeryToStandard(void)
{
  int *ptrValue1, *ptrValue2;
  int currentDegree;
  memcpy(&TestNbr, &powerMod, powerMod.nbrLimbs*sizeof(limb));
  GetMontgomeryParms(powerMod.nbrLimbs);
  degree = values[0];
  ptrValue1 = &values[1];
  ptrValue2 = &poly4[0];
  for (currentDegree = 0; currentDegree <= degree; currentDegree++)
  {
    UncompressBigInteger(ptrValue1, &operand1);
    // Convert operand1 from Montgomery to standard notation.
    operand2.limbs[0].x = 1;
    if (NumberLength > 1)
    {
      memset(&operand2.limbs[1], 0, (NumberLength - 1)*sizeof(limb));
    }
    modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    CompressBigInteger(ptrValue2, &operand1);
    ptrValue1 += 1 + *ptrValue1;
    ptrValue2 += 1 + *ptrValue2;
  }
  memcpy(&values[1], &poly4[0], (ptrValue2 - &poly4[0])*sizeof(int));
}

static void showPower(char **pptrOutput, int exponent)
{
  char *ptrOutput = *pptrOutput;
  if (superscripts)
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

static void showPowerX(char **pptrOutput, int polyDegree)
{
  char *ptrOutput = *pptrOutput;
  if (polyDegree == 0)
  {
    *ptrOutput++ = '1';
  }
  else
  {
    *ptrOutput++ = 'x';
    if (polyDegree != 1)
    {
      showPower(&ptrOutput, polyDegree);
    }
  }
  *pptrOutput = ptrOutput;
}

static void showPolynomial(char **pptrOutput, int *ptrPoly, int polyDegree, int groupLen)
{
  int nbrLimbs = NumberLength + 1;
  int *ptrValue1 = ptrPoly + (polyDegree - 1)*nbrLimbs;
  int currentDegree;
  char *ptrOutput = *pptrOutput;

  for (currentDegree = polyDegree - 1; currentDegree >= 0; currentDegree--)
  {
    if (*ptrValue1 != 1 || *(ptrValue1 + 1) != 0)
    {            // Coefficient is not zero.
      *ptrOutput++ = ' ';
      *ptrOutput++ = '+';
      *ptrOutput++ = ' ';
      if (*ptrValue1 != 1 || *(ptrValue1 + 1) != 1)
      {            // Coefficient is not one.
        UncompressBigInteger(ptrValue1, &operand1);
        Bin2Dec(operand1.limbs, ptrOutput, operand1.nbrLimbs, groupLen);
        ptrOutput += strlen(ptrOutput);
        if (currentDegree > 0)
        {
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
}

void outputPolynomial(char *ptrOutput, int groupLen)
{
  struct sFactorInfo *pstFactorInfo = factorInfo;
  int nbrFactor, currentDegree;
  int *ptrValue1, *ptrValue2;
  int nbrLimbs = powerMod.nbrLimbs + 1;
  // Output polynomial to factor. First move polynomial to poly4
  degree = values[0];
  ptrValue1 = &values[1];
  ptrValue2 = &poly4[0];
  for (currentDegree = 0; currentDegree <= degree; currentDegree++)
  {
    memcpy(ptrValue2, ptrValue1, (1 + *ptrValue1)*sizeof(int));
    ptrValue1 += 1 + *ptrValue1;
    ptrValue2 += nbrLimbs;
  }

  if (onlyEvaluate)
  {
    strcpy(ptrOutput, "<p>");
  }
  else
  {
    strcpy(ptrOutput, "<p id=\"pol\">");
  }
  ptrOutput += strlen(ptrOutput);
  UncompressBigInteger(&poly4[degree*nbrLimbs], &operand1);
  if ((operand1.nbrLimbs != 1 || operand1.limbs[0].x != 1) || degree == 0)
  {     // Leading coefficient is not 1 or degree is zero.
    Bin2Dec(operand1.limbs, ptrOutput, operand1.nbrLimbs, groupLen);
    ptrOutput += strlen(ptrOutput);
  }
  if (degree > 0)
  {
    showPowerX(&ptrOutput, degree);
  }
  showPolynomial(&ptrOutput, poly4, degree, groupLen);
  strcpy(ptrOutput, " (mod ");
  ptrOutput += strlen(ptrOutput);
  Bin2Dec(primeMod.limbs, ptrOutput, primeMod.nbrLimbs, groupLen);
  ptrOutput += strlen(ptrOutput);
  if (exponentMod != 1)
  {
    showPower(&ptrOutput, exponentMod);
  }
  strcpy(ptrOutput, ")</p>");
  ptrOutput += strlen(ptrOutput);
  if (onlyEvaluate)
  {
    return;
  }
  strcpy(ptrOutput,"<p>");
  ptrOutput += strlen(ptrOutput);

  // Output factors
  *ptrOutput++ = '<';
  *ptrOutput++ = 'u';
  *ptrOutput++ = 'l';
  *ptrOutput++ = '>';
  UncompressBigInteger(&poly4[degree*nbrLimbs], &operand1);
  if ((operand1.nbrLimbs != 1 || operand1.limbs[0].x != 1) || nbrFactorsFound == 0)
  {     // Leading coefficient is not 1 or degree is zero.
    *ptrOutput++ = '<';
    *ptrOutput++ = 'l';
    *ptrOutput++ = 'i';
    *ptrOutput++ = '>';
    Bin2Dec(operand1.limbs, ptrOutput, operand1.nbrLimbs, groupLen);
    ptrOutput += strlen(ptrOutput);
    *ptrOutput++ = '<';
    *ptrOutput++ = '/';
    *ptrOutput++ = 'l';
    *ptrOutput++ = 'i';
    *ptrOutput++ = '>';
  }
  for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    int polyDegree = pstFactorInfo->degree;
    int multiplicity = pstFactorInfo->multiplicity;
    *ptrOutput++ = '<';
    *ptrOutput++ = 'l';
    *ptrOutput++ = 'i';
    *ptrOutput++ = '>';
    if (multiplicity > 1)
    {
      *ptrOutput++ = '(';
    }
    showPowerX(&ptrOutput, polyDegree);
    showPolynomial(&ptrOutput, pstFactorInfo->ptrPolyLifted, polyDegree, groupLen);
    if (multiplicity > 1)
    {
      *ptrOutput++ = ')';
      showPower(&ptrOutput, multiplicity);
    }
    *ptrOutput++ = '<';
    *ptrOutput++ = '/';
    *ptrOutput++ = 'l';
    *ptrOutput++ = 'i';
    *ptrOutput++ = '>';
    pstFactorInfo++;
  }
  *ptrOutput++ = '<';
  *ptrOutput++ = '/';
  *ptrOutput++ = 'u';
  *ptrOutput++ = 'l';
  *ptrOutput++ = '>';
  *ptrOutput++ = '\0';
}

void textErrorPol(char *output, enum eExprErr rc)
{
  char text[150];

  switch (rc)
  {
  case EXPR_CANNOT_USE_X_IN_EXPONENT:
    strcpy(text, lang?"No se puede usar x en el exponente":
                        "Cannot use x in exponent");
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
  *output++ = '<';
  *output++ = 'p';
  *output++ = '>';
  strcpy(output, text);
  output += strlen(output);
  *output++ = '<';
  *output++ = '/';
  *output++ = 'p';
  *output++ = '>';
  *output = 0;    // Add terminator character.
}


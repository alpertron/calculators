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
#ifndef _EXPRESSION_H
#define _EXPRESSION_H

#define COPYRIGHT_SPANISH "Hecho por Darío Alpern. Actualizado el 18 de abril de 2021."
#define COPYRIGHT_ENGLISH "Written by Dario Alpern. Last updated on 18 April 2021."

#include <stdbool.h>
#ifdef __EMSCRIPTEN__
int stamp(void);
#endif
#if defined _USING64BITS_ && defined __EMSCRIPTEN__
#define EXTERNALIZE  __attribute__((visibility("default")))
#else
#define EXTERNALIZE	
#endif
void databack(char *data);

enum eExprErr
{
  EXPR_NUMBER_TOO_LOW = -100,
  EXPR_NUMBER_TOO_HIGH,
  EXPR_INTERM_TOO_HIGH,
  EXPR_DIVIDE_BY_ZERO,
  EXPR_PAREN_MISMATCH,
  EXPR_SYNTAX_ERROR,
  EXPR_TOO_MANY_PAREN,
  EXPR_INVALID_PARAM,
  EXPR_ARGUMENTS_NOT_RELATIVELY_PRIME,
  EXPR_BREAK,
  EXPR_OUT_OF_MEMORY,
  EXPR_CANNOT_USE_X_IN_EXPONENT,
  EXPR_POLYNOMIAL_DIVISION_NOT_INTEGER,
  EXPR_DEGREE_TOO_HIGH,
  EXPR_EXPONENT_TOO_LARGE,
  EXPR_EXPONENT_NEGATIVE,
  EXPR_LEADING_COFF_MULTIPLE_OF_PRIME,
  EXPR_CANNOT_LIFT,
  EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE,
  EXPR_MODULUS_MUST_BE_PRIME_EXP,
  EXPR_BASE_MUST_BE_POSITIVE,
  EXPR_POWER_MUST_BE_POSITIVE,
  EXPR_MODULUS_MUST_BE_NONNEGATIVE,
  EXPR_VAR_OR_COUNTER_REQUIRED,
  EXPR_MULTIPLE_VARIABLES_NOT_ACCEPTED,
  EXPR_OK = 0,
  EXPR_NOT_FOUND,
};
#ifndef lang  
  extern bool lang;
#endif
extern char inputString[1000000];
extern char output[3000000];
extern BigInteger valueX;
extern int counterC;
extern int expressionNbr;
enum eExprErr ComputeGaussianExpression(char *expr, BigInteger *ExpressionResult);
enum eExprErr ComputeExpression(char *expr, int typ, BigInteger *ExpressionResult);
void partition(int val, BigInteger *pResult);
void factorial(BigInteger *result, int argument);
void primorial(BigInteger *result, int argument);
void textError(char * ptrOutput, enum eExprErr rc);
void initializeSmallPrimes(int* pSmallPrimes);
#else
#endif

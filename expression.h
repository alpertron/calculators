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

#define COPYRIGHT_SPANISH "Hecho por Darío Alpern. Actualizado el 27 de julio de 2023."
#define COPYRIGHT_ENGLISH "Written by Dario Alpern. Last updated on 27 July 2023."

#include "bignbr.h"
#ifdef __EMSCRIPTEN__
int stamp(void);
#endif
#if defined _USING64BITS_ && defined __EMSCRIPTEN__ && !defined _MSC_VER
#define EXTERNALIZE  __attribute__((visibility("default")))
#else
#define EXTERNALIZE	
#endif
void databack(const char *data);

#define STACK_OPER_SIZE             100

#define ONE_PARM                     (1 * 0x100)
#define TWO_PARMS                    (2 * 0x100)
#define THREE_PARMS                  (3 * 0x100)
#define FOUR_PARMS                   (4 * 0x100)
#define MANY_PARMS                   (8 * 0x100)
#define NO_PARMS                     (9 * 0x100)
#define TOKEN_NUMBER                  1
#define TOKEN_START_EXPON             2
#define TOKEN_END_EXPON               3
#define TOKEN_VAR                     4
#define TOKEN_COUNTER                 5
#define TOKEN_FACTORIAL               6
#define OPER_POWER                    7
#define OPER_MULTIPLY                 8
#define OPER_DIVIDE                   9
#define OPER_REMAINDER               10
#define OPER_UNARY_MINUS             11
#define OPER_ADD                     12
#define OPER_SUBT                    13
#define OPER_SHR                     14
#define OPER_SHL                     15
#define OPER_NOT_GREATER             16
#define OPER_NOT_LESS                17
#define OPER_NOT_EQUAL               18
#define OPER_EQUAL                   19
#define OPER_GREATER                 20
#define OPER_LESS                    21
#define OPER_NOT                     22
#define OPER_AND                     23
#define OPER_INFIX_AND               24
#define OPER_OR                      25
#define OPER_INFIX_OR                26
#define OPER_XOR                     27
#define OPER_PAREN                   28
#define TOKEN_RANDOM                 29
#define TOKEN_ABS                    30
#define TOKEN_SGN                    31
#define TOKEN_ANS                    32


enum eExprErr
{
  EXPR_NUMBER_TOO_LOW = -100,
  EXPR_NUMBER_TOO_HIGH,
  EXPR_INTERM_TOO_HIGH,
  EXPR_DIVIDE_BY_ZERO,
  EXPR_PAREN_MISMATCH,
  EXPR_SYNTAX_ERROR,
  EXPR_LITERAL_NOT_INTEGER,
  EXPR_INVALID_PARAM,
  EXPR_ARGUMENTS_NOT_RELATIVELY_PRIME,
  EXPR_OUT_OF_MEMORY,
  EXPR_CANNOT_USE_X_IN_EXPONENT,
  EXPR_POLYNOMIAL_DIVISION_NOT_INTEGER,
  EXPR_DENOMINATOR_MUST_BE_CONSTANT,
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
  EXPR_TOO_MANY_ARGUMENTS,
  EXPR_TOO_FEW_ARGUMENTS,
  EXPR_MORE_THAN_ONE_EQUAL_SIGN,
  EXPR_EQUAL_SIGN_INSIDE_PAREN,
  EXPR_VAR_IN_EXPRESSION,
  EXPR_ARGUMENT_MUST_BE_CONSTANT,
  EXPR_CANNOT_PARSE_EXPRESSION,
  EXPR_SHORT_CIRCUIT,
  EXPR_OK = 0,
  EXPR_NOT_FOUND,
};

struct sFuncOperExpr
{
  const char* name;
  short token;
  char priority;
};

enum eParseExpr
{
  PARSE_EXPR_INTEGER = 0,
  PARSE_EXPR_GAUSSIAN,
  PARSE_EXPR_POLYNOMIAL,
};

#ifndef lang  
  extern bool lang;
#endif
extern char inputString[1000000];
extern char output[3000000];
extern BigInteger valueX;
extern int counterC;
extern int expressionNbr;
enum eExprErr ComputeGaussianExpression(const char *expr, BigInteger *ExpressionResult);
enum eExprErr ComputeExpression(const char *ptrStartRPN, BigInteger *ExpressionResult);
int ConvertToReversePolishNotation(const char* input, char** pptrOut,
  const struct sFuncOperExpr* funcOperExpr, enum eParseExpr eParseExpr,
  bool* pUsingVariables, bool *pUsingRandom);
void partition(int val, BigInteger *pResult);
enum eExprErr factorial(BigInteger *result, int argument, int multifact);
enum eExprErr primorial(BigInteger *result, int argument);
enum eExprErr ComputeFibLucas(int origValue, BigInteger* pArgument);
enum eExprErr ComputeBack(BigInteger* pArgument);
enum eExprErr ComputeNext(BigInteger* pArgument);
void textError(char **pptrOutput, enum eExprErr rc);
void initializeSmallPrimes(int* pSmallPrimes);
enum eExprErr parseNumberInsideExpr(const char** ppInput, char** ppOutput);
void setInsideExpressionLoop(bool inside);
bool isInsideExpressionLoop(void);
enum eExprErr convertToRPN(const char* expr, char** pointerRPNbuffer, bool varsExpected);
#endif

#ifndef _EXPRESSION_H
#define _EXPRESSION_H

#define COPYRIGHT_SPANISH "Hecho por Darío Alpern. Actualizado el 3 de noviembre de 2017."
#define COPYRIGHT_ENGLISH "Written by Dario Alpern. Last updated on 3 November 2017."

#ifdef __EMSCRIPTEN__
void databack(char *data);
int stamp(void);
#endif

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
  EXPR_BREAK,
  EXPR_OUT_OF_MEMORY,
  EXPR_CANNOT_USE_X_IN_EXPONENT,
  EXPR_DEGREE_TOO_HIGH,
  EXPR_EXPONENT_TOO_LARGE,
  EXPR_EXPONENT_NEGATIVE,
  EXPR_LEADING_COFF_MULTIPLE_OF_PRIME,
  EXPR_CANNOT_LIFT,
  EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE,
  EXPR_MODULUS_MUST_BE_PRIME_EXP,
  EXPR_BASE_MUST_BE_POSITIVE,
  EXPR_POWER_MUST_BE_POSITIVE,
  EXPR_MODULUS_BASE_NOT_RELATIVELY_PRIME,
  EXPR_MODULUS_POWER_NOT_RELATIVELY_PRIME,
  EXPR_MODULUS_MUST_BE_NONNEGATIVE,
  EXPR_OK = 0
};
extern int lang;
extern char inputString[1000000];
extern char output[3000000];
enum eExprErr ComputeGaussianExpression(char *expr, BigInteger *ExpressionResult);
enum eExprErr ComputeExpression(char *expr, int type, BigInteger *ExpressionResult);
void partition(int val, BigInteger *pResult);
void textError(char *output, enum eExprErr rc);
#else
#endif

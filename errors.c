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
#include <stdbool.h>
#include "bignbr.h"
#include "expression.h"

int expressionNbr;
void textError(char **pptrOutput, enum eExprErr rc)
{
  char* ptrOut = *pptrOutput;
  switch (rc)
  {
  case EXPR_NUMBER_TOO_LOW:
    copyStr(&ptrOut, lang ? "Número muy pequeño" : "Number too low");
    break;
  case EXPR_NUMBER_TOO_HIGH:
    copyStr(&ptrOut, lang ? "Número muy grande (más de 10000 dígitos)" :
      "Number too high (more than 10000 digits)");
    break;
  case EXPR_INTERM_TOO_HIGH:
    copyStr(&ptrOut, lang ? "Número intermedio muy grande (más de 20000 dígitos" :
      "Intermediate number too high (more than 20000 digits)");
    break;
  case EXPR_DIVIDE_BY_ZERO:
    copyStr(&ptrOut, lang ? "División por cero" : "Division by zero");
    break;
  case EXPR_PAREN_MISMATCH:
    copyStr(&ptrOut, lang ? "Error de paréntesis" : "Parenthesis mismatch");
    break;
  case EXPR_SYNTAX_ERROR:
    if (lang)
    {
      copyStr(&ptrOut, "Error de sintaxis");
      if (expressionNbr > 0)
      {
        copyStr(&ptrOut, " en la expresión ");
        *ptrOut = (char)(expressionNbr + '0');
        ptrOut++;
        *ptrOut = 0;
      }
    }
    else
    {
      copyStr(&ptrOut, "Syntax error");
      if (expressionNbr > 0)
      {
        copyStr(&ptrOut, " in expression #");
        *ptrOut = (char)(expressionNbr + '0');
        ptrOut++;
        *ptrOut = 0;
      }
    }
    break;
  case EXPR_TOO_MANY_PAREN:
    copyStr(&ptrOut, lang ? "Demasiados paréntesis" : "Too many parenthesis");
    break;
  case EXPR_INVALID_PARAM:
    copyStr(&ptrOut, lang ? "Parámetro inválido" : "Invalid parameter");
    break;
  case EXPR_ARGUMENTS_NOT_RELATIVELY_PRIME:
    copyStr(&ptrOut, lang ? "MCD de los argumentos no es 1" : "GCD of arguments is not 1");
    break;
  case EXPR_BREAK:
    copyStr(&ptrOut, lang ? "Detenido por el usuario" : "Stopped by user");
    break;
  case EXPR_VAR_OR_COUNTER_REQUIRED:
    if (lang)
    {
      copyStr(&ptrOut, "La expresión ");
      *ptrOut = (char)(expressionNbr + '0');
      ptrOut++;
      copyStr(&ptrOut, " debe incluir la variable <var>x</var> y/o el contador <var>c</var>");
    }
    else
    {
      copyStr(&ptrOut, "Expression #");
      *ptrOut = (char)(expressionNbr + '0');
      ptrOut++;
      copyStr(&ptrOut, " must include the variable <var>x</var> and/or the counter <var>c</var>");
    }
    break;
  default:
    break;
  }
  *pptrOutput = ptrOut;
}

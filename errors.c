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
#include "string/strings.h"
#include "bignbr.h"
#include "expression.h"
#include "copyStr.h"

int expressionNbr;
void textError(char **pptrOutput, enum eExprErr rc)
{
  char* ptrOut = *pptrOutput;
  switch (rc)
  {
  case EXPR_NUMBER_TOO_LOW:
    copyStr(&ptrOut, LITERAL_TEXT_ERROR1);
    break;
  case EXPR_NUMBER_TOO_HIGH:
    copyStr(&ptrOut, LITERAL_TEXT_ERROR2);
    break;
  case EXPR_INTERM_TOO_HIGH:
    copyStr(&ptrOut, LITERAL_TEXT_ERROR3);
    break;
  case EXPR_DIVIDE_BY_ZERO:
    copyStr(&ptrOut, LITERAL_TEXT_ERROR4);
    break;
  case EXPR_PAREN_MISMATCH:
    copyStr(&ptrOut, LITERAL_TEXT_ERROR5);
    break;
  case EXPR_LITERAL_NOT_INTEGER:
    copyStr(&ptrOut, LITERAL_TEXT_ERROR6);
    break;
  case EXPR_INTERNAL_ERROR:
    copyStr(&ptrOut, LITERAL_TEXT_ERROR7);
    break;
  case EXPR_SYNTAX_ERROR:
    formatString(&ptrOut, LITERAL_TEXT_ERROR8, expressionNbr);
    break;
  case EXPR_INVALID_PARAM:
    copyStr(&ptrOut, LITERAL_TEXT_ERROR9);
    break;
  case EXPR_TOO_FEW_ARGUMENTS:
    copyStr(&ptrOut, LITERAL_TEXT_ERROR10);
    break;
  case EXPR_TOO_MANY_ARGUMENTS:
    copyStr(&ptrOut, LITERAL_TEXT_ERROR11);
    break;
  case EXPR_ARGUMENTS_NOT_RELATIVELY_PRIME:
    copyStr(&ptrOut, LITERAL_TEXT_ERROR12);
    break;
  case EXPR_VAR_OR_COUNTER_REQUIRED:
    // Expression #$1d must include the variable $2v and/or the counter $3v
    formatString(&ptrOut, LITERAL_TEXT_ERROR13, expressionNbr, 'x', 'c');
    break;
  case EXPR_VAR_IN_EXPRESSION:
    copyStr(&ptrOut, LITERAL_TEXT_ERROR14);
    break;
  case EXPR_CANNOT_PARSE_EXPRESSION:
    copyStr(&ptrOut, LITERAL_TEXT_ERROR15);
    break;
  default:
    break;
  }
  *pptrOutput = ptrOut;
}

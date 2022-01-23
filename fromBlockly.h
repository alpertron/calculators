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
#ifndef _FROMBLOCKLY_H
#define _FROMBLOCKLY_H

#define BLOCKLY_INSTR_BUFFER_SIZE    10000000
#define BLOCKLY_INSTR_BUFFER_PTR_SIZE 1000000
#define BLOCKLY_INSTR_STACK_SIZE         5000
#define BLOCKLY_NBR_VARIABLES             200   // It cannot be greater than 255.
#define BLOCKLY_VARIABLE_NAME_LEN          50

enum
{
  BLOCKLY_NO_ERROR = 0,
  BLOCKLY_TOO_MANY_VARIABLES,
  BLOCKLY_VARIABLE_NAME_TOO_LONG,
  BLOCKLY_ONE_TOP_BLOCK,
  BLOCKLY_INVALID_BLOCK_TYPE,
  BLOCKLY_INVALID_NAME_FOR_SHADOW,
  BLOCKLY_TOO_MANY_CLOSING_BLOCK_TAGS,
  BLOCKLY_NO_OUTPUT_EXPECTED
  };

void fromBlockly(const char* ptrXMLFromBlockly);
void setBlocklyVar(int index, const BigInteger* value);
void getBlocklyVar(int index, BigInteger* value);
#endif
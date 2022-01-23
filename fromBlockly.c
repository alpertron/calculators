//
// This file is part of Alpertron Calculators.
//
// Copyright 2019-2021 Dario Alejandro Alpern
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
// Input: XML exported from Blockly.
#include <string.h>
#include "bignbr.h"
#include "fromBlockly.h"
#include "expression.h"
#include "exprfact.h"
#include "linkedbignbr.h"
#include "output.h"
#include <stdio.h>

#ifdef __EMSCRIPTEN__
#define SEND_DATA_TO_OUTPUT(n)   databack(n)
#else
#define SEND_DATA_TO_OUTPUT(n)   printf("%s\n",n)
#endif
char* ptrBlocklyOutput;
int nbrBlocklyOutputLines;
static char bufferInstr[BLOCKLY_INSTR_BUFFER_SIZE];
static char variableNames[BLOCKLY_NBR_VARIABLES * BLOCKLY_VARIABLE_NAME_LEN];
static char blockStack[BLOCKLY_INSTR_STACK_SIZE];
static int stackIf[BLOCKLY_INSTR_STACK_SIZE];
static int nbrVariables;
static int nbrNamedVariables;
static char* ptrInstr;
struct linkedBigInt* pstVariables[BLOCKLY_NBR_VARIABLES];
static char startRepeatCode[] =
{
  TOKEN_NUMBER, 0, 1, 1, 0, 0, 0, TOKEN_SET_VAR
};

// Conversion from block number to operator/function
static const char blockBasicMath[] =
{
  OPER_ADD,
  OPER_SUBT,
  OPER_MULTIPLY,
  OPER_DIVIDE,
  OPER_REMAINDER,
  OPER_POWER,
  TOKEN_SQRT,
  TOKEN_RANDOM,
  TOKEN_ABS,
  TOKEN_SGN,
};

static const char blockComparisons[] =
{
  OPER_EQUAL,
  OPER_NOT_EQUAL,
  OPER_GREATER,
  OPER_NOT_GREATER,
  OPER_LESS,
  OPER_NOT_LESS,
};

static const char blockLogic[] =
{
  OPER_AND,
  OPER_OR,
  OPER_XOR,
  OPER_NOT,
  OPER_SHL,
  OPER_SHR,
};

static const char blockDivisibility[] =
{
  TOKEN_GCD,
  TOKEN_LCM,
  TOKEN_ISPRIME,
  TOKEN_NUMFACT,
  TOKEN_MINFACT,
  TOKEN_MAXFACT,
  TOKEN_NUMDIVS,
  TOKEN_SUMDIVS,
};

static const char blockRecreationalMath[] =
{
  TOKEN_N,
  TOKEN_B,
  TOKEN_NUMDIGITS,
  TOKEN_SUMDIGITS,
  TOKEN_REVDIGITS,
  TOKEN_CONCATFACT,
};

static const char blockNumberTheory[] =
{
  TOKEN_MODINV,
  TOKEN_MODPOW,
  TOKEN_TOTIENT,
  TOKEN_JACOBI,
};

static const char blockOther[] =
{
  TOKEN_FACTORIAL1,
  TOKEN_FACTORIAL,
  TOKEN_PRIMORIAL,
  TOKEN_F,
  TOKEN_L,
  TOKEN_P,
};

static const char blockOutput[] =
{
  TOKEN_PRINT,
  TOKEN_PRINT_HEX,
  TOKEN_PRINTFACT,
  TOKEN_PRINTFACT_HEX,
  TOKEN_PRINTPRIME,
  TOKEN_PRINTPRIME_HEX,
};

static const char *blockToToken[] =
{
  blockBasicMath,
  blockComparisons,
  blockLogic,
  blockDivisibility,
  blockRecreationalMath,
  blockNumberTheory,
  blockOther,
  blockOutput,
};

/*
  Tags XML:
  *variables*: container for list of variables.
  *variable*: represents a variable.
               Format: <variable id="...">var</variable>
  *block*: if it includes fields x and y, it is a top level block.
  Exactly one top level block outside procedures is accepted.
  *value*: it represents an operand inside a block. It can include
  another block, a shadow or a field.
  *shadow*: field that can be overwritten by a block. If a block
  is immediately after a shadow, the shadow is discarded.
  *field*: it indicates a number or a variable:
           <field name="NUM">number</field> or
           <field name="VAR">variable</field> or
  *next*: it indicates the next item in a stack of statements.
          This tag can be ignored.
  Set variable:
  <block type="variables_set" id="...">
    <field name="VAR" id = "...">i</field>
    <value name="VALUE">
      <block type="math_number" id = "...">
        <field name="NUM">123</field>
      </block>
    </value>
  </block>
  
  Get variable:
  <block type="variables_get" id="...">
    <field name="VAR" id="...">i</field>
  </block>
 
  XML is in Polish notation. Use a stack to store block types
  and then the tag </block> to store block type in instruction buffer.
  
  Encoding for FOR loop:
  
  for VAR FROM TO BY
 
  FROM set_var
  to:
  get_var TO <= TOKEN_IF end_loop
  TOKEN_JMP do
  by:
  get_var BY + set_var TOKEN_JMP to
  do:
  ....
  JMP by
  end_loop
*/

static int xmlcmp(const char* source, const char* cmp)
{
  return memcmp(source, cmp, strlen(cmp));
}

static void skipID(const char** pptrXML)
{
  const char* ptrXML = *pptrXML;
  while (xmlcmp(ptrXML, "id=") != 0)
  {
    ptrXML++;
  }
  *pptrXML = strchr(ptrXML + 4, '\"') + 2;         // Skip both quotes.
}

static void setJmpOffset(int offsetJmp, int offsetTarget)
{
  unsigned int shiftedTarget;
  if (offsetJmp >= 0)
  {
    bufferInstr[offsetJmp] = (char)offsetTarget;
    shiftedTarget = (unsigned int)offsetTarget >> 8;
    bufferInstr[offsetJmp + 1] = (char)shiftedTarget;
    shiftedTarget = (unsigned int)offsetTarget >> 16;
    bufferInstr[offsetJmp + 2] = (char)shiftedTarget;
    shiftedTarget = (unsigned int)offsetTarget >> 24;
    bufferInstr[offsetJmp + 3] = (char)shiftedTarget;
  }
}

static int getVariableNbr(const char** ppXML)
{
  int varNameSize;
  const char* ptrVariableNames = variableNames;
  size_t varNameSize_t;
  // Point to end of block tag.
  const char* ptrXML = strchr(*ppXML, '>');
  // Point to start of value tag.
  ptrXML++;
  // Point to end of value tag.
  ptrXML = strchr(ptrXML, '>');
  // Point to start of name of variable.
  ptrXML++;
  // Get length of name of variable.
  varNameSize_t = strchr(ptrXML, '<') - ptrXML;
  varNameSize = (int)varNameSize_t;
  for (int varNbr = 0; varNbr < nbrNamedVariables; varNbr++)
  {
    int index;
    int curVarNameSize = *ptrVariableNames;
    if (curVarNameSize != varNameSize)
    {
      continue;
    }
    ptrVariableNames++;
    for (index = 0; index < curVarNameSize; index++)
    {
      if (*(ptrVariableNames + index) != *(ptrXML + index))
      {
        break;
      }
    }
    if (index == curVarNameSize)
    {        // Variable name found.
      *ppXML = ptrXML;
      return varNbr;
    }
    ptrVariableNames += curVarNameSize; // Point to next name.
  }
  return -1;
}

static void showBlocklyError(int rc)
{
  switch (rc)
  {
  case BLOCKLY_TOO_MANY_VARIABLES:
    SEND_DATA_TO_OUTPUT(lang ? "KDemasiadas variables definidas" :
      "KToo many variables defined");
    return;
  case BLOCKLY_VARIABLE_NAME_TOO_LONG:
    SEND_DATA_TO_OUTPUT(lang ? "KEl nombre de la variable es muy largo" :
      "KVariable name is too long");
    return;
  case BLOCKLY_ONE_TOP_BLOCK:
    SEND_DATA_TO_OUTPUT(lang ? "KDebe haber un solo conjunto de bloques unidos" :
      "KThere must be only one set of joined blocks");
    return;
  case BLOCKLY_INVALID_BLOCK_TYPE:
    SEND_DATA_TO_OUTPUT(lang ? "KTipo de bloque inválido" :
      "KInvalid block type.");
    return;
  case BLOCKLY_INVALID_NAME_FOR_SHADOW:
    SEND_DATA_TO_OUTPUT(lang ? "KNombre inválido para el tag shadow en XML" :
      "KInvalid name for shadow tag in XML");
    return;
  case BLOCKLY_TOO_MANY_CLOSING_BLOCK_TAGS:
    SEND_DATA_TO_OUTPUT(lang ? "KDemasiados tags /block" :
      "KToo many /block tags");
    return;
  case BLOCKLY_NO_OUTPUT_EXPECTED:
    SEND_DATA_TO_OUTPUT(lang ? "KEl bloque suelto debe estar dentro de otro bloque" :
      "KThe detached block must be inside another block");
    return;
  default:
    return;
  }
}

static void EncodeElse(int *ptrStackControl)
{
  // Encode JMP after DO statement.
  size_t offsetEndIf = ptrInstr - &bufferInstr[0];
  *ptrInstr = TOKEN_JMP;
  *(ptrStackControl - 1) = (int)offsetEndIf + 1;
  ptrInstr++;              // Skip token.
  ptrInstr += sizeof(int); // Skip jump offset.
  // Next ELSE IF.
  offsetEndIf = ptrInstr - &bufferInstr[0];
  setJmpOffset(*(ptrStackControl - 2), (int)offsetEndIf);
  *(ptrStackControl - 2) = -1;
}

static int parseBlocklyXml(const char* ptrXMLFromBlockly)
{
  const char* ptrXML;
  int nbrTopBlocks = 0;
  char* ptrLastVariableName;
  char* ptrBlockStack = blockStack;
  char* ptrStartShadow = NULL;
  int* ptrStackControl = stackIf;

  ptrInstr = bufferInstr;
  nbrVariables = 0;
  // Point to first tag after initial XML tag.
  ptrXML = strchr(ptrXMLFromBlockly + 1, '<');
  if (xmlcmp(ptrXML, "<variables") == 0)
  {   // Blockly XML has variables.
    ptrXML = strchr(ptrXML + 1, '<');  // Point to first variable
    ptrLastVariableName = variableNames;
    for (;;)
    {
      size_t variableNameSize;
      // Point to next variable.
      const char* startVariable;
      const char* endVariable;
      if (xmlcmp(ptrXML, "</") == 0)
      {
        ptrXML = strchr(ptrXML, '>');
        ptrXML++;   // Point after tag "</variables>"
        break;
      }
      if (nbrVariables == BLOCKLY_NBR_VARIABLES)
      {    // Too many variables.
        return BLOCKLY_TOO_MANY_VARIABLES;
      }
      startVariable = strchr(ptrXML, '>') + 1;
      endVariable = strchr(startVariable, '<');
      variableNameSize = endVariable - startVariable;
      if ((endVariable - startVariable) > 127)
      {
        return BLOCKLY_VARIABLE_NAME_TOO_LONG;
      }
      if (ptrLastVariableName + variableNameSize >=
          &variableNames[sizeof(variableNames)])
      {    // Too many variables.
        return BLOCKLY_TOO_MANY_VARIABLES;
      }
      *ptrLastVariableName = (char)variableNameSize;
      ptrLastVariableName++;
      (void)memcpy(ptrLastVariableName, startVariable, variableNameSize);
      ptrLastVariableName += variableNameSize;
      nbrVariables++;
      ptrXML = strchr(endVariable, '>') + 1;
    }
  }
  nbrNamedVariables = nbrVariables;
  // Convert blocks expressed in XML to operators and operands in prefix format.
  while (*ptrXML != 0)
  {
    if (xmlcmp(ptrXML, "<block") == 0)
    {                 // Block XML found.
      bool isTopBlock = false;
      const char* ptrXMLbak = ptrXML;
      skipID(&ptrXMLbak);
      if (*ptrXMLbak == 'x')
      {
        isTopBlock = true;
        nbrTopBlocks++;
        if (nbrTopBlocks != 1)
        {
          return BLOCKLY_ONE_TOP_BLOCK;
        }
      }
      ptrXML += 13;   // Point to block type.
      if (*(ptrXML + 2) == '\"')
      {     // Custom Blockly type
        const char* ptrBlockTranslation;
        size_t nbrGroups;
        int groupNbr = *ptrXML - 'A';
        if (isTopBlock && (groupNbr != 7))
        {
          return BLOCKLY_NO_OUTPUT_EXPECTED;
        }
        nbrGroups = sizeof(blockToToken) / sizeof(blockToToken[0]);
        if ((groupNbr < 0) || (groupNbr >= (int)nbrGroups))
        {
          return BLOCKLY_INVALID_BLOCK_TYPE;
        }
        ptrBlockTranslation = blockToToken[groupNbr];
        *ptrBlockStack = *(ptrBlockTranslation + *(ptrXML + 1) - 'A');
        ptrBlockStack++;
      }
      else if (xmlcmp(ptrXML, "M") == 0)
      {
        if (isTopBlock)
        {
          return BLOCKLY_NO_OUTPUT_EXPECTED;
        }
        ptrXML = strchr(ptrXML, '>') + 1;
        ptrXML = strchr(ptrXML, '>') + 2; // Point to byte after start of number.
        if (*(ptrXML - 1) == '-')
        {
          ptrXML++;     // Skip minus sign.
          (void)parseNumberInsideExpr(&ptrXML, &ptrInstr);
          *ptrInstr = OPER_UNARY_MINUS;
          ptrInstr++;
        }
        else
        {
          (void)parseNumberInsideExpr(&ptrXML, &ptrInstr);
        }
      }
      else if (xmlcmp(ptrXML, "controls_if") == 0)
      {            // IF block.
        *ptrStackControl = -1;
        *(ptrStackControl + 1) = -1;
        ptrStackControl += 2;
      }
      else if (xmlcmp(ptrXML, "controls_whileUntil") == 0)
      {      // While or Until block.
        size_t offsetStartWhile = ptrInstr - &bufferInstr[0];
        *ptrStackControl = -1;
        *(ptrStackControl + 1) = (int)offsetStartWhile;
        ptrStackControl += 2;
        ptrXML = strchr(ptrXML, '>') + 1;   // Skip block.
        ptrXML = strchr(ptrXML, '>') + 1;   // Skip field.
        if (xmlcmp(ptrXML, "WHILE") == 0)
        {
          *ptrBlockStack = START_WHILE;
        }
        else
        {
          *ptrBlockStack = START_UNTIL;
        }
        ptrBlockStack++;
        ptrXML = strchr(ptrXML, '>') + 1;   // Skip value.
      }
      else if (xmlcmp(ptrXML, "controls_repeat_ext") == 0)
      {      // Repeat block.
        if (nbrVariables >= BLOCKLY_NBR_VARIABLES)
        {    // Too many variables.
          return BLOCKLY_TOO_MANY_VARIABLES;
        }
        // Store data needed for statements processing in block stack
        *ptrBlockStack = (char)nbrVariables;
        ptrBlockStack++;
        *ptrBlockStack = START_REPEAT;
        ptrBlockStack++;
        // Store instructions to set the variable to zero.
        (void)memcpy(ptrInstr, startRepeatCode, sizeof(startRepeatCode));
        ptrInstr += sizeof(startRepeatCode);
        *ptrInstr = (char)nbrVariables;
        ptrInstr++;
        // Store this address to go back at the end of the loop.
        size_t offsetRepeat = ptrInstr - &bufferInstr[0];
        *(ptrStackControl + 1) = (int)offsetRepeat;
        ptrStackControl += 2;
        // Get count variable.
        *ptrInstr = (char)TOKEN_GET_VAR;
        ptrInstr++;
        *ptrInstr = (char)nbrVariables;
        ptrInstr++;
        nbrVariables++;
      }
      else if (xmlcmp(ptrXML, "controls_for") == 0)
      {      // For block.
        *ptrBlockStack = (char)getVariableNbr(&ptrXML);
        ptrXML = strchr(ptrXML, '>') + 1;
        ptrBlockStack++;
        *ptrBlockStack = (char)nbrVariables;
        ptrBlockStack++;
        *ptrBlockStack = START_FOR;
        ptrBlockStack++;     
      }
      else if (xmlcmp(ptrXML, "variables_set") == 0)
      {
        *ptrBlockStack = (char)getVariableNbr(&ptrXML);
        ptrXML = strchr(ptrXML, '>') + 1;
        ptrBlockStack++;
        *ptrBlockStack = TOKEN_SET_VAR;
        ptrBlockStack++;
      }
      else if (xmlcmp(ptrXML, "variables_get") == 0)
      {
        if (isTopBlock)
        {
          return BLOCKLY_NO_OUTPUT_EXPECTED;
        }
        *ptrBlockStack = (char)getVariableNbr(&ptrXML);
        ptrXML = strchr(ptrXML, '>') + 1;
        ptrBlockStack++;
        *ptrBlockStack = TOKEN_GET_VAR;
        ptrBlockStack++;
      }
      else
      { // No more cases.
      }
    }
    else if (xmlcmp(ptrXML, "<value") == 0)
    {
      ptrXML += 13;  // Point to value name.
      if (xmlcmp(ptrXML, "TIMES") == 0)
      {
        // Nothing to do.
      }
      else if (xmlcmp(ptrXML, "TO") == 0)
      {     // Set counter variable in FOR. 
        *ptrInstr = TOKEN_SET_VAR;
        ptrInstr++;
        *ptrInstr = *(ptrBlockStack - 3);
        ptrInstr++;
      }
      else if (xmlcmp(ptrXML, "BY") == 0)
      {     // Set limit variable in FOR.
        *ptrInstr = TOKEN_SET_VAR;
        ptrInstr++;
        *ptrInstr = *(ptrBlockStack - 2);
        ptrInstr++;
      }
      else if (xmlcmp(ptrXML, "IF") == 0)
      {
        *ptrBlockStack = START_IF;
        ptrBlockStack++;
        if (*(ptrXML + 2) != '0')
        { // Not first IF in IF - ELSE IF - ELSE chain.
          EncodeElse(ptrStackControl);
        }
      }
      else
      {   // No more cases.
      }
    }
    else if (xmlcmp(ptrXML, "<statement") == 0)
    {
      if (xmlcmp(ptrXML + 17, "DO") == 0)
      {
        size_t offsetIf;
        while (ptrBlockStack > blockStack)
        {
          ptrBlockStack--;
          if (*ptrBlockStack == START_FOR)
          {     // Set increment variable in FOR.
                // abs(increment) * sgn(to - counter)
            *ptrInstr = TOKEN_ABS;
            ptrInstr++;
            *ptrInstr = TOKEN_GET_VAR;
            ptrInstr++;
            *ptrInstr = *(ptrBlockStack - 1);          // to
            ptrInstr++;
            *ptrInstr = TOKEN_GET_VAR;
            ptrInstr++;
            *ptrInstr = *(ptrBlockStack - 2);          // counter
            ptrInstr++;
            *ptrInstr = OPER_SUBT;
            ptrInstr++;
            *ptrInstr = TOKEN_SGN;
            ptrInstr++;
            *ptrInstr = OPER_MULTIPLY;
            ptrInstr++;
            *ptrInstr = TOKEN_SET_VAR;
            ptrInstr++;
            *ptrInstr = *(ptrBlockStack - 1) + 1;      // increment
            ptrInstr++;
            ptrStackControl += 2;
            offsetIf = ptrInstr - &bufferInstr[0];
            *(ptrStackControl - 1) = (int)offsetIf;
            // Test limit. If sgn(counter - to) == sgn(increment)
            *ptrInstr = TOKEN_GET_VAR;
            ptrInstr++;
            *ptrInstr = *(ptrBlockStack - 2);          // counter
            ptrInstr++;
            *ptrInstr = TOKEN_GET_VAR;
            ptrInstr++;
            *ptrInstr = *(ptrBlockStack - 1);          // to
            ptrInstr++;
            *ptrInstr = OPER_EQUAL;
            ptrInstr++;
            *ptrInstr = TOKEN_GET_VAR;
            ptrInstr++;
            *ptrInstr = *(ptrBlockStack - 1);          // to
            ptrInstr++;
            *ptrInstr = TOKEN_GET_VAR;
            ptrInstr++;
            *ptrInstr = *(ptrBlockStack - 2);          // counter
            ptrInstr++;
            *ptrInstr = OPER_SUBT;
            ptrInstr++;
            *ptrInstr = TOKEN_SGN;
            ptrInstr++;
            *ptrInstr = TOKEN_GET_VAR;
            ptrInstr++;
            *ptrInstr = *(ptrBlockStack - 1) + 1;      // increment
            ptrInstr++;
            *ptrInstr = TOKEN_SGN;
            ptrInstr++;
            *ptrInstr = OPER_EQUAL;
            ptrInstr++;
            *ptrInstr = OPER_OR;
            ptrInstr++;
            *ptrInstr = TOKEN_IF;
            ptrInstr++;
            offsetIf = ptrInstr - &bufferInstr[0];
            *(ptrStackControl - 2) = (int)offsetIf;
            ptrInstr += 4;
            break;
          }
          if (*ptrBlockStack == START_UNTIL)
          {
            *ptrInstr = OPER_NOT;
            ptrInstr++;
          }
          if (*ptrBlockStack == START_REPEAT)
          {
            *ptrInstr = OPER_NOT_GREATER;
            ptrInstr++;
            *ptrInstr = TOKEN_IF;
            ptrInstr++;
            offsetIf = ptrInstr - &bufferInstr[0];
            *(ptrStackControl - 2) = (int)offsetIf;
            ptrInstr += 4;
            break;
          }
          if ((*ptrBlockStack == START_UNTIL) ||
             (*ptrBlockStack == START_WHILE) ||
             (*ptrBlockStack == START_IF))
          {
            *ptrInstr = TOKEN_IF;
            ptrInstr++;
            offsetIf = ptrInstr - &bufferInstr[0];
            *(ptrStackControl - 2) = (int)offsetIf;
            ptrInstr += 4;
            break;
          }
          if ((*ptrBlockStack == TOKEN_GET_VAR) || (*ptrBlockStack == TOKEN_SET_VAR))
          {
            *ptrInstr = *ptrBlockStack;
            ptrInstr++;
            ptrBlockStack--;
          }
          *ptrInstr = *ptrBlockStack;
          ptrInstr++;
        }
      }
      else if (xmlcmp(ptrXML + 17, "ELSE") == 0)
      {
        while (ptrBlockStack > blockStack)
        {
          ptrBlockStack--;
          if (*ptrBlockStack == START_IF)
          {
            break;
          }
          else if (*ptrBlockStack >= START_FLOW_CONTROL)
          {
            ptrBlockStack++;
            break;
          }
          else if ((*ptrBlockStack == TOKEN_GET_VAR) || (*ptrBlockStack == TOKEN_SET_VAR))
          {
            *ptrInstr = *ptrBlockStack;
            ptrInstr++;
            ptrBlockStack--;
            *ptrInstr = *ptrBlockStack;
            ptrInstr++;
          }
          else
          {
            *ptrInstr = *ptrBlockStack;
            ptrInstr++;
          }
        }
        EncodeElse(ptrStackControl);
        *ptrBlockStack = START_ELSE;
      }
      else
      {
        *ptrBlockStack = START_STATEMENT;
      }
      ptrBlockStack++;
    }
    else if (xmlcmp(ptrXML, "<shadow") == 0)
    {
      ptrStartShadow = ptrInstr;
      ptrXML += 14;  // Point to shadow type.
      if (xmlcmp(ptrXML, "M") == 0)
      {
        ptrXML = strchr(ptrXML, '>') + 1;
        ptrXML = strchr(ptrXML, '>') + 2; // Point to byte after start of number.
        if (*(ptrXML - 1) == '-')
        {
          ptrXML++;     // Skip minus sign.
          (void)parseNumberInsideExpr(&ptrXML, &ptrInstr);
          *ptrInstr = OPER_UNARY_MINUS;
          ptrInstr++;
        }
        else
        {
          (void)parseNumberInsideExpr(&ptrXML, &ptrInstr);
        }
        ptrXML += 17;    // Discard "</field></shadow> tags.
        if ((*ptrXML == '<') && xmlcmp(ptrXML, "<block") == 0)
        {   // Discard number just parsed (block has priority over shadow).
          ptrInstr = ptrStartShadow;
          continue;
        }
      }
      else
      {
        return BLOCKLY_INVALID_NAME_FOR_SHADOW;
      }
    }
    else if (xmlcmp(ptrXML, "<next") == 0)
    {
      if (ptrBlockStack > blockStack)
      {
        ptrBlockStack--;
        if (*ptrBlockStack >= START_FLOW_CONTROL)
        {     // Ignore these entries in stack of blocks.
          ptrBlockStack++;
        }
        else if ((*ptrBlockStack == TOKEN_GET_VAR) || (*ptrBlockStack == TOKEN_SET_VAR))
        {
          *ptrInstr = *ptrBlockStack;
          ptrInstr++;
          ptrBlockStack--;
          *ptrInstr = *ptrBlockStack;
          ptrInstr++;
        }
        else
        {
          *ptrInstr = *ptrBlockStack;
          ptrInstr++;
        }
      }
    }
    else if (xmlcmp(ptrXML, "</statement") == 0)
    {
      while (ptrBlockStack > blockStack)
      {
        ptrBlockStack--;
        if (*ptrBlockStack == START_STATEMENT)
        {
          break;
        }
        if (*ptrBlockStack == START_FOR)
        {
          size_t offsetJmpStartFor;
          size_t offsetAfterForExpression;
          *ptrInstr = TOKEN_GET_VAR;
          ptrInstr++;
          *ptrInstr = *(ptrBlockStack - 2);          // counter
          ptrInstr++;
          *ptrInstr = TOKEN_GET_VAR;
          ptrInstr++;
          *ptrInstr = *(ptrBlockStack - 1) + 1;      // increment
          ptrInstr++;
          *ptrInstr = OPER_ADD;
          ptrInstr++;
          *ptrInstr = TOKEN_SET_VAR;
          ptrInstr++;
          *ptrInstr = *(ptrBlockStack - 2);          // counter
          ptrInstr++;
          ptrBlockStack -= 2;
          // Jump to beginning of loop.
          *ptrInstr = TOKEN_JMP;
          ptrInstr++;
          offsetJmpStartFor = ptrInstr - &bufferInstr[0];
          setJmpOffset((int)offsetJmpStartFor, *(ptrStackControl - 1));
          ptrInstr += 4;
          offsetAfterForExpression = ptrInstr - &bufferInstr[0];
          setJmpOffset(*(ptrStackControl - 2), (int)offsetAfterForExpression);
          ptrStackControl -= 2;
          break;
        }
        if (*ptrBlockStack == START_REPEAT)
        {
          size_t offsetJmpStartRepeat;
          size_t offsetAfterRepeatExpression;
          // Increment repeat variable.
          *ptrInstr = TOKEN_GET_VAR;
          ptrInstr++;
          *ptrInstr = *(ptrBlockStack-1);
          ptrInstr++;
          *ptrInstr = 1;
          ptrInstr++;
          *ptrInstr = 0;
          ptrInstr++;
          *ptrInstr = 1;
          ptrInstr++;
          *ptrInstr = 1;
          ptrInstr++;
          *ptrInstr = 0;
          ptrInstr++;
          *ptrInstr = 0;
          ptrInstr++;
          *ptrInstr = 0;
          ptrInstr++;
          *ptrInstr = OPER_ADD;
          ptrInstr++;
          *ptrInstr = TOKEN_SET_VAR;
          ptrInstr++;
          *ptrInstr = *(ptrBlockStack - 1);
          ptrInstr++;
          // Jump to beginning of loop.
          *ptrInstr = TOKEN_JMP;
          ptrInstr++;
          offsetJmpStartRepeat = ptrInstr - &bufferInstr[0];
          setJmpOffset((int)offsetJmpStartRepeat, *(ptrStackControl - 1));
          ptrInstr += 4;
          offsetAfterRepeatExpression = ptrInstr - &bufferInstr[0];
          setJmpOffset(*(ptrStackControl - 2), (int)offsetAfterRepeatExpression);
          ptrStackControl -= 2;
          ptrBlockStack--;
          break;
        }
        if ((*ptrBlockStack == START_WHILE) ||
          (*ptrBlockStack == START_UNTIL))
        {
          size_t offsetJmpStartWhile;
          size_t offsetAfterWhileExpression;
          *ptrInstr = TOKEN_JMP;
          ptrInstr++;
          offsetJmpStartWhile = ptrInstr - &bufferInstr[0];
          setJmpOffset((int)offsetJmpStartWhile, *(ptrStackControl - 1));
          ptrInstr += 4;
          offsetAfterWhileExpression = ptrInstr - &bufferInstr[0];
          setJmpOffset(*(ptrStackControl - 2), (int)offsetAfterWhileExpression);
          ptrStackControl -= 2;
          break;
        }
        if (*ptrBlockStack == START_IF)
        {
          // *(ptrStackIf-1) = Point to JMP just after DO statement.
          // *(ptrStackIf-2) = Point to next ELSE IF.
          size_t offsetEndIf = ptrInstr - &bufferInstr[0];
          setJmpOffset(*(ptrStackControl - 1), (int)offsetEndIf);
          if ((xmlcmp(ptrXML+12, "</block>") == 0))
          {    // End of block IF.
            setJmpOffset(*(ptrStackControl - 2), (int)offsetEndIf);
            ptrStackControl -= 2;
          }
          break;
        }
        if (*ptrBlockStack == START_ELSE)
        {
          // *(ptrStackIf-1) = Point to JMP just after DO statement.
          // *(ptrStackIf-2) = Point to next ELSE IF.
          size_t offsetEndIf = ptrInstr - &bufferInstr[0];
          setJmpOffset(*(ptrStackControl - 1), (int)offsetEndIf);
          setJmpOffset(*(ptrStackControl - 2), (int)offsetEndIf);
          ptrStackControl -= 2;     // End of IF block.
          break;
        }
        if ((*ptrBlockStack == TOKEN_GET_VAR) || (*ptrBlockStack == TOKEN_SET_VAR))
        {
          *ptrInstr = *ptrBlockStack;
          ptrInstr++;
          ptrBlockStack--;
        }
        *ptrInstr = *ptrBlockStack;
        ptrInstr++;
      }
    }
    else
    {   // No more cases.
    }
    ptrXML = strchr(ptrXML + 1, '>') + 1;  // Point to next XML entry.
  }
  while (ptrBlockStack > blockStack)
  {
    ptrBlockStack--;
    if ((*ptrBlockStack == TOKEN_GET_VAR) || (*ptrBlockStack == TOKEN_SET_VAR))
    {
      *ptrInstr = *ptrBlockStack;
      ptrInstr++;
      ptrBlockStack--;
    }
    *ptrInstr = *ptrBlockStack;
    ptrInstr++;
  }
  if (nbrTopBlocks != 1)
  {
    return BLOCKLY_ONE_TOP_BLOCK;
  }
  *ptrInstr = '\0';    // Add program terminator.
  return BLOCKLY_NO_ERROR;
}

void fromBlockly(const char* ptrXMLFromBlockly)
{
  int rc = parseBlocklyXml(ptrXMLFromBlockly);
  if (rc != BLOCKLY_NO_ERROR)
  {
    showBlocklyError(rc);
    return;
  }
  ptrBlocklyOutput = output;
  copyStr(&ptrBlocklyOutput, "2<ul>");
  nbrBlocklyOutputLines = 0;
  // Use linked big integers for Blockly variables.
  initLinkedBigInt();
  for (int index = 0; index < BLOCKLY_NBR_VARIABLES; index++)
  {   // Initialize all variables to zero.
    intToLinkedBigInt(&pstVariables[index], 0);
  }
#ifdef __EMSCRIPTEN__
  databack("L");  // Exit Blockly mode.
#endif
  (void)ComputeExpression(bufferInstr, NULL);
  if (nbrBlocklyOutputLines == 0)
  {
    copyStr(&ptrBlocklyOutput, lang? "<li>No hay nada para mostrar</li>":
      "<li>There is nothing to print</li>");
  }
  copyStr(&ptrBlocklyOutput, "</ul>");
  beginLine(&ptrBlocklyOutput);
  copyStr(&ptrBlocklyOutput, lang ? COPYRIGHT_SPANISH : COPYRIGHT_ENGLISH );
  finishLine(&ptrBlocklyOutput);
  SEND_DATA_TO_OUTPUT(output);
}

void setBlocklyVar(int index, const BigInteger* value)
{
  setLinkedBigInteger(&pstVariables[index], value);
}

void getBlocklyVar(int index, BigInteger* value)
{
  getBigIntegerFromLinked(pstVariables[index], value);
}

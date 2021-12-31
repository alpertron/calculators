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
#include "fromBlockly.h"
#include "bignbr.h"
#include "expression.h"
#include "exprfact.h"

static char bufferInstr[BLOCKLY_INSTR_BUFFER_SIZE];
static char variableNames[BLOCKLY_NBR_VARIABLES][BLOCKLY_VARIABLE_NAME_LEN];
static char blockStack[BLOCKLY_INSTR_STACK_SIZE];
static int stackIf[BLOCKLY_INSTR_STACK_SIZE];
static int nbrVariables;
static char* ptrInstr;

// Conversion from block number to operator/function
static const char blockBasicMath[] =
{
  OPER_PLUS,
  OPER_MINUS,
  OPER_MULTIPLY,
  OPER_DIVIDE,
  OPER_REMAINDER,
  OPER_POWER,
  TOKEN_SQRT,
  TOKEN_ABS,
  TOKEN_RANDOM,
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
  bufferInstr[offsetJmp] = (char)offsetTarget;
  bufferInstr[offsetJmp+1] = (char)(offsetTarget >> 8);
  bufferInstr[offsetJmp+2] = (char)(offsetTarget >> 16);
  bufferInstr[offsetJmp+3] = (char)(offsetTarget >> 24);
}

void fromBlockly(const char* ptrXMLFromBlockly)
{
  const char* ptrXML;
  int nbrTopBlocks = 0;
  char* ptrBlockStack = blockStack;
  char* ptrStartShadow = NULL;
  int* ptrStackControl = stackIf;

  ptrInstr = bufferInstr;
  nbrVariables = 0;
  // Point to first tag after initial XML tag.
  ptrXML = strchr(ptrXMLFromBlockly + 1, '<');
  if (xmlcmp(ptrXML, "<variables") == 0)
  {   // Blockly XML has variables.
    for (;;)
    {
      // Point to next variable.
      const char* startVariable;
      const char* endVariable;
      ptrXML = strchr(ptrXML + 1, '<');
      if (xmlcmp(ptrXML, "</") == 0)
      {
        ptrXML = strchr(ptrXML, '>');
        ptrXML++;   // Point after tag "</variables>"
        break;
      }
      if (nbrVariables == BLOCKLY_NBR_VARIABLES)
      {    // Too many variables.
#ifdef __EMSCRIPTEN__
        databack("KToo many variables defined");
#endif
        return;
      }
      startVariable = strchr(ptrXML, '>');
      startVariable++;
      endVariable = strchr(startVariable, '<');
      memcpy(variableNames[nbrVariables], startVariable, endVariable - startVariable);
      ptrXML = strchr(endVariable, '>');
      ptrXML++;
    }
  }
  // Convert blocks expressed in XML to operators and operands in prefix format.
  while (*ptrXML != 0)
  {
    if (xmlcmp(ptrXML, "<block") == 0)
    {                 // Block XML found.
      const char* ptrXMLbak = ptrXML;
      skipID(&ptrXMLbak);
      if (*ptrXMLbak == 'x')
      {
        nbrTopBlocks++;
        if (nbrTopBlocks > 1)
        {
#ifdef __EMSCRIPTEN__
          databack("KOnly one top block expected.");
#endif
          return;
        }
      }
      ptrXML += 13;   // Point to block type.
      if (*(ptrXML + 2) == '\"')
      {     // Custom Blockly type
        const char* ptrBlockTranslation;
        int groupNbr = *ptrXML - 'A';
        size_t nbrGroups = sizeof(blockToToken) / sizeof(blockToToken[0]);
        if ((groupNbr < 0) || (groupNbr >= (int)nbrGroups))
        {
#ifdef __EMSCRIPTEN__
          databack("KInvalid block type.");
#endif
          return;
        }
        ptrBlockTranslation = blockToToken[groupNbr];
        *ptrBlockStack = *(ptrBlockTranslation + *(ptrXML + 1) - 'A');
        ptrBlockStack++;
      }
      else if (xmlcmp(ptrXML, "M") == 0)
      {
        ptrXML = strchr(ptrXML, '>') + 1;
        ptrXML = strchr(ptrXML, '>') + 2; // Point to byte after start of number.
        *ptrInstr = TOKEN_NUMBER;
        ptrInstr++;
        (void)parseNumberInsideExpr(&ptrXML, &ptrInstr);
        *ptrBlockStack = NUMBER_FIELD;
        ptrBlockStack++;
      }
      else if (xmlcmp(ptrXML, "controls_if") == 0)
      {            // IF block.
        *ptrStackControl = -1;
        *(ptrStackControl + 1) = -1;
        ptrStackControl += 2;
      }
      else if (xmlcmp(ptrXML, "controls_whileUntil") == 0)
      {
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
      else
      {
      }
    }
    else if (xmlcmp(ptrXML, "<value") == 0)
    {
      ptrXML += 13;  // Point to value name.
      if (xmlcmp(ptrXML, "IF") == 0)
      {
        *ptrBlockStack = START_IF;
        ptrBlockStack++;
      }
      else if (*ptrXML < '1' || *ptrXML > '3')
      {              // Values outside custom block.
#ifdef __EMSCRIPTEN__
        databack("KInvalid name for value in XML");
#endif
        return;
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
          if (*ptrBlockStack == START_UNTIL)
          {
            *ptrInstr = OPER_NOT;
            ptrInstr++;
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
          if (*ptrBlockStack != NUMBER_FIELD)
          {
            *ptrInstr = *ptrBlockStack;
            ptrInstr++;
          }
        }
      }
      else if (xmlcmp(ptrXML + 17, "ELSE") == 0)
      {
        size_t offsetIf;
        while (ptrBlockStack > blockStack)
        {
          ptrBlockStack--;
          if (*ptrBlockStack == START_IF)
          {
            break;
          }
          if (*ptrBlockStack != NUMBER_FIELD)
          {
            *ptrInstr = *ptrBlockStack;
            ptrInstr++;
          }
        }
        *ptrInstr = TOKEN_IF;
        ptrInstr++;
        offsetIf = ptrInstr - &bufferInstr[0];
        *(ptrStackControl - 2) = (int)offsetIf;
        ptrInstr += 4;
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
        *ptrInstr = TOKEN_NUMBER;
        ptrInstr++;
        (void)parseNumberInsideExpr(&ptrXML, &ptrInstr);
        ptrXML += 17;    // Discard "</field></shadow> tags.
        if (*ptrXML == '<' && xmlcmp(ptrXML, "<block") == 0)
        {   // Discard number just parsed (block has priority over shadow).
          ptrInstr = ptrStartShadow;
          continue;
        }
      }
      else
      {
#ifdef __EMSCRIPTEN__
        databack("KInvalid name for shadow in XML");
#endif
        return;
      }
    }
    else if (xmlcmp(ptrXML, "<next") == 0)
    {
      ptrBlockStack--;
      if (ptrBlockStack < blockStack)
      {
#ifdef __EMSCRIPTEN__
        databack("KToo many /block tags");
#endif
        return;
      }
      if ((*ptrBlockStack == START_IF) || (*ptrBlockStack == START_ELSE))
      {     // Ignore these entries in stack of blocks.
        ptrBlockStack++;
      }
      if (*ptrBlockStack != NUMBER_FIELD)
      {
        *ptrInstr = *ptrBlockStack;
        ptrInstr++;
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
        else if (*ptrBlockStack == START_IF)
        {
          // *(ptrStackIf-1) = Point to JMP just after DO statement.
          // *(ptrStackIf-2) = Point to next ELSE IF.
          size_t offsetEndIf = ptrInstr - &bufferInstr[0];
          setJmpOffset(*(ptrStackControl - 1), (int)offsetEndIf);
          // Encode JMP after DO statement.
          *ptrInstr = TOKEN_JMP;
          *(ptrStackControl - 1) = (int)offsetEndIf + 1;
          ptrInstr++;              // Skip token.
          ptrInstr += sizeof(int); // Skip jump offset.
          // Next ELSE IF.
          offsetEndIf = ptrInstr - &bufferInstr[0];
          if (*(ptrStackControl - 2) >= 0)
          {
            setJmpOffset(*(ptrStackControl - 2), (int)offsetEndIf);
          }
          break;
        }
        if (*ptrBlockStack == START_ELSE)
        {
          // *(ptrStackIf-1) = Point to JMP just after DO statement.
          // *(ptrStackIf-2) = Point to next ELSE IF.
          size_t offsetEndIf = ptrInstr - &bufferInstr[0];
          setJmpOffset(*(ptrStackControl - 1), (int)offsetEndIf);
          if (*(ptrStackControl - 2) >= 0)
          {
            setJmpOffset(*(ptrStackControl - 2), (int)offsetEndIf);
          }
          ptrStackControl -= 2;     // End of IF block.
          break;
        }
        if (*ptrBlockStack != NUMBER_FIELD)
        {
          *ptrInstr = *ptrBlockStack;
          ptrInstr++;
        }
      }
    }
    ptrXML = strchr(ptrXML + 1, '>') + 1;  // Point to next XML entry.
  }
  while (ptrBlockStack > blockStack)
  {
    ptrBlockStack--;
    if (*ptrBlockStack != NUMBER_FIELD)
    {
      *ptrInstr = *ptrBlockStack;
      ptrInstr++;
    }
  }
}


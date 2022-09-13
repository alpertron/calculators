//
// This file is part of Alpertron Calculators.
//
// Copyright 2022 Dario Alejandro Alpern
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
#include "bignbr.h"
#include "expression.h"
#include "highlevel.h"
#include "showtime.h"
#include "batch.h"
#include "tsquares.h"
static BigInteger Base1;
static BigInteger Base2;
static BigInteger Base3;
static BigInteger valueN;
static BigInteger powerN;
static BigInteger toProcess;
static int groupLength;
static const char* cube = "<span class=\"bigger\">³</span>";
static const char* fifth = "<sup>5</sup>";
static const char* seventh = "<sup>7</sup>";
static const char* powerStr;
static void batchSqCubesCallback(char** pptrOutput);
static int power4;
static int Exponent = 3;

static bool checkSumOfTwoSquares(const BigInteger *pArgument, int expon)
{
  int lenBytes;
  int numLimbs;
  int mod8;
  power4 = 0;
  (void)BigIntPowerIntExp(&Base3, expon, &powerN);
  BigIntSubt(pArgument, &powerN, &powerN);
  // Test whether powerN is a perfect square.
  numLimbs = powerN.nbrLimbs;
  lenBytes = numLimbs * (int)sizeof(limb);
  (void)memcpy(Base2.limbs, powerN.limbs, lenBytes);
  squareRoot(Base2.limbs, Base1.limbs,
    numLimbs, &Base1.nbrLimbs);
  multiply(Base1.limbs, Base1.limbs,
    Base2.limbs, Base1.nbrLimbs, NULL);
  if (memcmp(Base2.limbs, powerN.limbs, lenBytes) == 0)
  {     // powerN can be expressed as a square of Base1.
    intToBigInteger(&Base2, 0);
    return true;
  }
  DivideBigNbrByMaxPowerOf4(&power4, powerN.limbs, &powerN.nbrLimbs);
  mod8 = powerN.limbs[0].x & 7;
  if ((mod8 > 2) && (mod8 != 5))
  {     // If number mod 8 is not 1, 2 or 5,
        // it cannot be expressed as sum of two squares.
    return false;
  }
  if ((powerN.nbrLimbs == 1) && (powerN.limbs[0].x == 2))
  {
    intToBigInteger(&Base1, 1);
    intToBigInteger(&Base2, 1);
    return true;
  }
  nbrLimbsP = powerN.nbrLimbs;
  nbrLimbs = nbrLimbsP;
  (void)memcpy(valueP, powerN.limbs, nbrLimbsP * sizeof(int));
  sum = 0;
  FillSieveArray(valueP);
  if (!isSumOfTwoSquares())
  {
    return false;
  }
  CopyBigInt(&Base1, &biMult1);
  CopyBigInt(&Base2, &biMult2);
  return true;
}

// Compute Base1, Base2 and Base3 such that 
// Base1^2 + Base2^2 + Base3^exp = *pArgument.
// Compute k = floor(root(pArgument)).
// Then for i = k, k - 1, k - 2, ... try to represent
// pArgument - i^exp as a sum of two squares.
static int tsqcubes(const BigInteger *pArgument, int expon)
{
  // Compute Base3 <- floor(root(pArgument)).
  (void)BigIntRoot(pArgument, &Base3, expon);
  InitSieveArray();
  while (!checkSumOfTwoSquares(pArgument, expon))
  { // If we cannot express powerN as a sum of two squares,
    // try previous value of Base3.
    addbigint(&Base3, -1); // Use previous value of Base3.
  }
  Base1.sign = SIGN_POSITIVE;
  Base2.sign = SIGN_POSITIVE;
  (void)BigIntMultiplyPower2(&Base1, power4);
  (void)BigIntMultiplyPower2(&Base2, power4);

  // Validate

  (void)BigIntMultiply(&Base1, &Base1, &valueN);
  (void)BigIntMultiply(&Base2, &Base2, &powerN);
  BigIntAdd(&powerN, &valueN, &valueN);
  (void)BigIntPowerIntExp(&Base3, expon, &powerN);
  BigIntAdd(&powerN, &valueN, &valueN);
  if (!TestBigNbrEqual(&valueN, pArgument))
  {
    return 1;       // Result does not validate.
  }
  return 0;
}

void tsqcubesText(char *input, int grpLen)
{
  char *ptrOutput;
  if (valuesProcessed == 0)
  {
    groupLength = grpLen;
  }
  (void)BatchProcessing(input, &toProcess, &ptrOutput, NULL, batchSqCubesCallback);
#ifdef __EMSCRIPTEN__
  copyStr(&ptrOutput, lang ? "<p>Transcurrió " : "<p>Time elapsed: ");
  int elapsedTime = (int)(tenths() - originalTenthSecond);
  GetDHMSt(&ptrOutput, elapsedTime);
#endif
  copyStr(&ptrOutput, "</p><p>");
  copyStr(&ptrOutput, (lang ? COPYRIGHT_SPANISH: COPYRIGHT_ENGLISH));
  copyStr(&ptrOutput, "</p>");
}

// Show cube number. Use parentheses for negative numbers.
static void showTerm(char** pptrOutput, const BigInteger* pBase, const char *expon)
{
  char* ptrOutput = *pptrOutput;
  if (pBase->sign == SIGN_NEGATIVE)
  {
    *ptrOutput = '(';
    ptrOutput++;
  }
  if (hexadecimal)
  {
    BigInteger2Hex(&ptrOutput, pBase, groupLength);
  }
  else
  {
    BigInteger2Dec(&ptrOutput, pBase, groupLength);
  }
  if (pBase->sign == SIGN_NEGATIVE)
  {
    *ptrOutput = ')';
    ptrOutput++;
  }
  copyStr(&ptrOutput, expon);
  *pptrOutput = ptrOutput;
}

static void batchSqCubesCallback(char **pptrOutput)
{
  int result;
  char *ptrOutput = *pptrOutput;
  NumberLength = toProcess.nbrLimbs;
  result = tsqcubes(&toProcess, Exponent);
  // Show the number to be decomposed into sum of cubes.
  copyStr(&ptrOutput, "<p>");
  if (hexadecimal)
  {
    BigInteger2Hex(&ptrOutput, &toProcess, groupLength);
  }
  else
  {
    BigInteger2Dec(&ptrOutput, &toProcess, groupLength);
  }
  if (result == 1)
  {
    copyStr(&ptrOutput, (lang?": ¡Error interno! Por favor envíe este número al autor del applet.</p>":
      ": Internal error! Please send the number to the author of the applet.</p>"));
    *pptrOutput = ptrOutput;
    return;
  }
  // Use powerN as a temporary variable.
  BigIntSubt(&Base1, &Base2, &powerN);
  if (powerN.sign == SIGN_NEGATIVE)
  {    // Swap Base1 and Base2 so the first is always greater than the second.
    CopyBigInt(&powerN, &Base1);
    CopyBigInt(&Base1, &Base2);
    CopyBigInt(&Base2, &powerN);
  }
  // Show decomposition in sum of up to two squares and a cube.
  copyStr(&ptrOutput, " = ");
  if (BigIntIsZero(&Base2) && BigIntIsZero(&Base1))
  {
    showTerm(&ptrOutput, &Base3, powerStr);
  }
  else
  {
    showTerm(&ptrOutput, &Base1, square);
    if (!BigIntIsZero(&Base2))
    {
      copyStr(&ptrOutput, " + ");
      showTerm(&ptrOutput, &Base2, square);
    }
    if (!BigIntIsZero(&Base3))
    {
      copyStr(&ptrOutput, " + ");
      showTerm(&ptrOutput, &Base3, powerStr);
    }
  }
  *pptrOutput = ptrOutput;
}

void assignExponent(char c)
{
  switch (c)
  {
  case '0':
    Exponent = 3;
    powerStr = cube;
    break;
  case '1':
    Exponent = 5;
    powerStr = fifth;
    break;
  default:
    Exponent = 7;
    powerStr = seventh;
    break;
  }
}

#if defined __EMSCRIPTEN__ && !defined _MSC_VER
EXTERNALIZE void doWork(void)
{
  int app;
  int grpLen = 0;
  char* ptrData = inputString;
  originalTenthSecond = tenths();
  if (*ptrData == 'C')
  {    // User pressed Continue button.
    tsqcubesText(NULL, 0); // Routine does not use parameters in this case.
    databack(output);
    return;
  }
  assignExponent(*ptrData);
  ptrData += 2;   // Skip from and comma.
  valuesProcessed = 0;
  while (*ptrData != ',')
  {
    grpLen = (grpLen * 10) + (*ptrData - '0');
    ptrData++;
  }
  ptrData++;             // Skip comma.
  app = *ptrData - '0';
  if (*(ptrData + 1) != ',')
  {
    ptrData++;
    app = (app * 10) + *ptrData - '0';
  }
#ifndef lang  
  lang = ((app & 1) ? true : false);
#endif
  if ((app & 0x40) != 0)
  {
    hexadecimal = true;
  }
  else
  {
    hexadecimal = false;
  }
  tsqcubesText(ptrData + 2, grpLen);
  databack(output);
}
#endif

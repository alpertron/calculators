//
// This file is part of Alpertron Calculators.
//
// Copyright 2016-2023 Dario Alejandro Alpern
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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bignbr.h"
#include "expression.h"
#include "factor.h"
#include "showtime.h"
#include "batch.h"
#include "output.h"
#include "commonstruc.h"
#include "fromBlockly.h"

#ifdef FACTORIZATION_APP
#ifdef __EMSCRIPTEN__
extern bool skipPrimality;
extern int64_t lModularMult;
#endif
extern bool fromFile;
extern bool lineEndingCRLF;
extern BigInteger tofactor;
extern BigInteger factorValue;
static BigInteger result;
static BigInteger Tmp;
static void GetEulerTotient(char **pptrOutput);
static void GetMobius(char **pptrOutput);
static void GetNumberOfDivisors(char **pptrOutput);
static void GetSumOfDivisors(char **pptrOutput);
void ComputeFourSquares(const struct sFactors* pstFactors);
void ShowFourSquares(char **pptrOutput);
bool doFactorization;
bool doShowPrime;
static char *knownFactors;
bool useBlockly;
void heapsort(int count, int **ptrBigInts);
void showSumTwoSquares(void);

#ifdef FACTORIZATION_APP
void batchEcmCallback(char **pptrOutput, int type)
{
  bool doFactBak = doFactorization;
  bool doShowPrimeBak = doShowPrime;
  if (type != BATCH_NO_QUOTE)
  {
    doShowPrime = false;
    if ((type & BATCH_MASK_PROCESS) != 0)
    {
      doFactorization = true;
    }
    else
    {
      doFactorization = false;
    }
  }
  char *ptrFactorDec = tofactorDec;
  *ptrFactorDec = 0;
  NumberLength = tofactor.nbrLimbs;
  BigInteger2IntArray(nbrToFactor, &tofactor);
  if (*nbrToFactor < 0)
  {    // If number is negative, make it positive.
    *nbrToFactor = -*nbrToFactor;
  }
  if (type == BATCH_NO_QUOTE)
  {
    if (tofactor.sign == SIGN_NEGATIVE)
    {
      *ptrFactorDec = '-';
      ptrFactorDec++;
    }
    if (hexadecimal)
    {
      Bin2Hex(&ptrFactorDec, tofactor.limbs, tofactor.nbrLimbs, groupLen);
    }
    else
    {
      Bin2Dec(&ptrFactorDec, tofactor.limbs, tofactor.nbrLimbs, groupLen);
    }
  }
  if (doFactorization)
  {
    factorExt(&tofactor, nbrToFactor, factorsMod, astFactorsMod, knownFactors);
    knownFactors = NULL;
  }
  else if (doShowPrime)
  {
    if (tofactor.sign == SIGN_NEGATIVE)
    {
      copyStr(pptrOutput, "-");
    }
    if (hexadecimal)
    {
      Bin2Hex(pptrOutput, tofactor.limbs, tofactor.nbrLimbs, groupLen);
    }
    else
    {
      Bin2Dec(pptrOutput, tofactor.limbs, tofactor.nbrLimbs, groupLen);
    }
    if (BpswPrimalityTest(&tofactor, NULL) == 0)
    {    // Argument is a probable prime.
      copyStr(pptrOutput, lang ? " es primo" : " is prime");
    }
    else
    {
      copyStr(pptrOutput, lang ? " no es primo" : " is not prime");
    }
  }
  else
  {
    if (tofactor.sign == SIGN_NEGATIVE)
    {
      copyStr(pptrOutput, "-");
    }
    if (hexadecimal)
    {
      Bin2Hex(pptrOutput, tofactor.limbs, tofactor.nbrLimbs, groupLen);
    }
    else
    {
      Bin2Dec(pptrOutput, tofactor.limbs, tofactor.nbrLimbs, groupLen);
    }
  }
  SendFactorizationToOutput(astFactorsMod, pptrOutput, doFactorization,
    (type != BATCH_NO_QUOTE));
  doFactorization = doFactBak;
  doShowPrime = doShowPrimeBak;
}
#endif

#ifdef __EMSCRIPTEN__
static void startList(char** pptrOutput)
{
  if (!fromFile)
  {
    copyStr(pptrOutput, "<ul>");
  }
}

static void startListLine(char** pptrOutput)
{
  if (!fromFile)
  {
    copyStr(pptrOutput, "<li>");
  }
}

static void endListLine(char** pptrOutput)
{
  if (fromFile)
  {
    copyStr(pptrOutput, (lineEndingCRLF ? "\r\n" : "\n"));
  }
  else
  {
    copyStr(pptrOutput, "</li>");
  }
}

static void endList(char** pptrOutput)
{
  if (!fromFile)
  {
    copyStr(pptrOutput, "</ul>");
  }
}
#endif

static void ExponentToBigInteger(int exponent, BigInteger *bigint)
{
  if ((unsigned int)exponent > MAX_VALUE_LIMB)
  {
    bigint->limbs[0].x = exponent - (int)MAX_VALUE_LIMB;
    bigint->limbs[1].x = 1;
    bigint->nbrLimbs = 2;
  }
  else
  {
    bigint->limbs[0].x = exponent + 1;
    bigint->nbrLimbs = 1;
  }
  bigint->sign = SIGN_POSITIVE;
}

// Find number of divisors as the product of all exponents plus 1.
static void GetNumberOfDivisors(char **pptrOutput)
{
  char *ptrOutput = *pptrOutput;
  const struct sFactors *pstFactor;
  result.limbs[0].x = 1;    // Set result to 1.
  result.nbrLimbs = 1;
  result.sign = SIGN_POSITIVE;
  pstFactor = &astFactorsMod[1];
  for (int factorNumber = 1; factorNumber <= astFactorsMod[0].multiplicity; factorNumber++)
  {
    const int *ptrFactor = pstFactor->ptrFactor;
    if ((*ptrFactor == 1) && (*(ptrFactor + 1) < 2))
    {                        // Factor is 1.
      break;
    }
    ExponentToBigInteger(pstFactor->multiplicity, &factorValue);
    (void)BigIntMultiply(&factorValue, &result, &result);
    pstFactor++;
  }
  beginLine(&ptrOutput);
  copyStr(&ptrOutput, lang ? "Cantidad de divisores: " : "Number of divisors: ");
  if (hexadecimal)
  {
    BigInteger2Hex(&ptrOutput, &result, groupLen);
  }
  else
  {
    BigInteger2Dec(&ptrOutput, &result, groupLen);
  }
  finishLine(&ptrOutput);
  *pptrOutput = ptrOutput;
}

static void GetSumOfDivisors(char **pptrOutput)
{
  char *ptrOutput = *pptrOutput;
  SumOfDivisors(&result);
  beginLine(&ptrOutput);
  copyStr(&ptrOutput, lang ? "Suma de divisores: " : "Sum of divisors: ");
  if (hexadecimal)
  {
    BigInteger2Hex(&ptrOutput, &result, groupLen);
  }
  else
  {
    BigInteger2Dec(&ptrOutput, &result, groupLen);
  }
  finishLine(&ptrOutput);
  *pptrOutput = ptrOutput;
}

static void GetEulerTotient(char **pptrOutput)
{
  char *ptrOutput = *pptrOutput;
  Totient(&result);
  beginLine(&ptrOutput);
  copyStr(&ptrOutput, lang ? "Phi de Euler: " : "Euler's totient: ");
  if (hexadecimal)
  {
    BigInteger2Hex(&ptrOutput, &result, groupLen);
  }
  else
  {
    BigInteger2Dec(&ptrOutput, &result, groupLen);
  }
  finishLine(&ptrOutput);
  *pptrOutput = ptrOutput;
}

// Find Mobius as zero if some exponent is > 1, 1 if the number of factors is even, -1 if it is odd.
static void GetMobius(char **pptrOutput)
{
  char *ptrOutput = *pptrOutput;
  const struct sFactors *pstFactor;
  int mobius = 1;
  pstFactor = &astFactorsMod[1];
  if ((astFactorsMod[0].multiplicity > 1) || (*pstFactor->ptrFactor != 1) ||
    (*(pstFactor->ptrFactor + 1) != 1))
  {                                // Number to factor is not 1.
    for (int factorNumber = 1; factorNumber <= astFactorsMod[0].multiplicity; factorNumber++)
    {
      if (pstFactor->multiplicity == 1)
      {
        mobius = -mobius;
      }
      else
      {
        mobius = 0;
      }
      pstFactor++;
    }
  }
  beginLine(&ptrOutput);
  copyStr(&ptrOutput, "Möbius: ");
  if (mobius < 0)
  {
    mobius = -mobius;
    *ptrOutput = '-';
    ptrOutput++;
  }
  int2dec(&ptrOutput, mobius);
  finishLine(&ptrOutput);
  *pptrOutput = ptrOutput;
}

void ecmFrontText(char *tofactorText, bool performFactorization, char *factors)
{
  char *ptrOutput;
  bool isBatch;
  knownFactors = factors;
  if (valuesProcessed == 0)
  {
    doFactorization = performFactorization;
  }
  if ((tofactorText != NULL) && ((*tofactorText & 0xDF) == 'X') &&
     (findChar(tofactorText + 1, ';') == NULL))
  {
#ifdef __EMSCRIPTEN__
    databack(doFactorization? "M": "N");    // Use polynomial factorization calculator.
#endif
  }
  enum eExprErr rc = BatchProcessing(tofactorText, &tofactor, &ptrOutput,
       &isBatch, batchEcmCallback);
  if (rc == EXPR_VAR_IN_EXPRESSION)
  {
#ifdef __EMSCRIPTEN__
    databack(doFactorization ? "M" : "N");  // Use polynomial factorization calculator.
#endif
  }
  if (!isBatch && (rc == EXPR_OK) && (counterC == -1) && doFactorization)
  {
#ifdef __EMSCRIPTEN__
    int64_t sumSquaresModMult = 0;
#endif
    if (tofactor.sign == SIGN_POSITIVE)
    {        // Number to factor is non-negative.
#ifdef __EMSCRIPTEN__
      char* ptrText;
#endif
      if (!BigIntIsZero(&tofactor))
      {      // Number to factor is not zero.
        GetNumberOfDivisors(&ptrOutput);
        GetSumOfDivisors(&ptrOutput);
        GetEulerTotient(&ptrOutput);
        GetMobius(&ptrOutput);
      }
#ifdef __EMSCRIPTEN__
      StepECM = 3;   // Show progress (in percentage) of sum of squares.
      ptrText = ShowFactoredPart(&tofactor, astFactorsMod);
      beginLine(&ptrText);
      copyStr(&ptrText, lang ? "<p>Hallando suma de cuadrados.</p>" :
        "<p>Searching for sum of squares.</p>");
      ShowLowerText();
      sumSquaresModMult = lModularMult;
#endif
      ComputeFourSquares(astFactorsMod);
      ShowFourSquares(&ptrOutput);
#ifdef __EMSCRIPTEN__
      sumSquaresModMult = lModularMult - sumSquaresModMult;
      StepECM = 0;   // Do not show progress.
#endif
    }
#ifdef __EMSCRIPTEN__
    copyStr(&ptrOutput, "<div id=\"divisors\"><p><button type=\"button\" id=\"showdiv\">");
    copyStr(&ptrOutput, lang? "Mostrar divisores": "Show divisors");
    copyStr(&ptrOutput, "</button></p></div>");
#endif
    beginLine(&ptrOutput);
    showElapsedTime(&ptrOutput);
    finishLine(&ptrOutput);
#ifdef __EMSCRIPTEN__
    if (lModularMult > 0)
    {
      beginLine(&ptrOutput);
      copyStr(&ptrOutput, lang ? "Multiplicaciones modulares:" :
        "Modular multiplications:");
      finishLine(&ptrOutput);
      startList(&ptrOutput);
      if ((lModularMult - primeModMult - SIQSModMult - sumSquaresModMult) > 0)
      {
        startListLine(&ptrOutput);
        copyStr(&ptrOutput, "ECM: ");
        long2dec(&ptrOutput, lModularMult - primeModMult - SIQSModMult - sumSquaresModMult);
        endListLine(&ptrOutput);
      }
      if (primeModMult > 0)
      {
        startListLine(&ptrOutput);
        copyStr(&ptrOutput, lang ? "Verificación de números primos probables: " :
          "Probable prime checking: ");
        long2dec(&ptrOutput, primeModMult);
        endListLine(&ptrOutput);
      }
      if (SIQSModMult > 0)
      {
        startListLine(&ptrOutput);
        copyStr(&ptrOutput, "SIQS: ");
        long2dec(&ptrOutput, SIQSModMult);
        endListLine(&ptrOutput);
      }
      if (sumSquaresModMult > 0)
      {
        startListLine(&ptrOutput);
        copyStr(&ptrOutput, lang ? "Suma de cuadrados: " : "Sum of squares: ");
        long2dec(&ptrOutput, sumSquaresModMult);
        endListLine(&ptrOutput);
      }
      endList(&ptrOutput);
    }
    if (nbrSIQS > 0)
    {
      uint64_t quotient;
      beginLine(&ptrOutput);
      copyStr(&ptrOutput, "SIQS:");
      startList(&ptrOutput);
      startListLine(&ptrOutput);
      int2dec(&ptrOutput, polynomialsSieved);
      copyStr(&ptrOutput, lang ? " polinomios utilizados" : " polynomials sieved");
      endListLine(&ptrOutput);
      startListLine(&ptrOutput);
      int2dec(&ptrOutput, trialDivisions);
      copyStr(&ptrOutput, lang ? " conjuntos de divisiones de prueba" : " sets of trial divisions");
      endListLine(&ptrOutput);
      startListLine(&ptrOutput);
      int2dec(&ptrOutput, smoothsFound);
      copyStr(&ptrOutput, lang ? " congruencias completas (1 de cada " :
        " smooth congruences found (1 out of every ");
      quotient = ValuesSieved / (uint64_t)smoothsFound;
      int2dec(&ptrOutput, (int)quotient);
      copyStr(&ptrOutput, lang ? " valores)</li><li>" : " values)</li><li>");
      int2dec(&ptrOutput, totalPartials);
      copyStr(&ptrOutput, lang ? " congruencias parciales (1 de cada " :
        " partial congruences found (1 out of every ");
      quotient = ValuesSieved / (uint64_t)totalPartials;
      int2dec(&ptrOutput, (int)quotient);
      copyStr(&ptrOutput, lang ? " valores)" : " values)");
      endListLine(&ptrOutput);
      startListLine(&ptrOutput);
      int2dec(&ptrOutput, partialsFound);
      copyStr(&ptrOutput, lang ? " congruencias parciales útiles</li><li>Tamaño de la matriz binaria: " :
        " useful partial congruences</li><li>Size of binary matrix: ");
      int2dec(&ptrOutput, matrixRows);
      copyStr(&ptrOutput, " &times; ");
      int2dec(&ptrOutput, matrixCols);
      endListLine(&ptrOutput);
      endList(&ptrOutput);
    }
    if ((nbrSIQS > 0) || (nbrECM > 0) || (nbrPrimalityTests > 0))
    {
      beginLine(&ptrOutput);
      copyStr(&ptrOutput, lang ? "Tiempos:" : "Timings:");
      startList(&ptrOutput);
      if (nbrPrimalityTests > 0)
      {
        startListLine(&ptrOutput);
        copyStr(&ptrOutput, lang ? "Test de primo probable de " : "Probable prime test of ");
        int2dec(&ptrOutput, nbrPrimalityTests);
        copyStr(&ptrOutput, lang ? " número" : " number");
        if (nbrPrimalityTests != 1)
        {
          *ptrOutput = 's';
          ptrOutput++;
        }
        *ptrOutput = ':';
        ptrOutput++;
        *ptrOutput = ' ';
        ptrOutput++;
        GetDHMSt(&ptrOutput, timePrimalityTests);
        endListLine(&ptrOutput);
      }
      if (nbrECM > 0)
      {
        startListLine(&ptrOutput);
        copyStr(&ptrOutput, lang ? "Factorización " : "Factoring ");
        int2dec(&ptrOutput, nbrECM);
        copyStr(&ptrOutput, lang ? " número" : " number");
        if (nbrECM != 1)
        {
          *ptrOutput = 's';
          ptrOutput++;
        }
        copyStr(&ptrOutput, lang ? " mediante ECM" : " using ECM:");
        *ptrOutput = ' ';
        ptrOutput++;
        GetDHMSt(&ptrOutput, timeECM - timeSIQS);
        endListLine(&ptrOutput);
      }
      if (nbrSIQS > 0)
      {
        startListLine(&ptrOutput);
        copyStr(&ptrOutput, lang ? "Factorización " : "Factoring ");
        int2dec(&ptrOutput, nbrSIQS);
        copyStr(&ptrOutput, lang ? " número" : " number");
        if (nbrSIQS != 1)
        {
          *ptrOutput = 's';
          ptrOutput++;
        }
        copyStr(&ptrOutput, lang ? " mediante SIQS" : " using SIQS:");
        *ptrOutput = ' ';
        ptrOutput++;
        GetDHMSt(&ptrOutput, timeSIQS);
        endListLine(&ptrOutput);
      }
      endList(&ptrOutput);
    }
    copyStr(&ptrOutput, "</div>");
#endif
  }
  beginLine(&ptrOutput);
  copyStr(&ptrOutput, lang ? COPYRIGHT_SPANISH : COPYRIGHT_ENGLISH );
  finishLine(&ptrOutput);
  // Initialize all exponents to zero.
  (void)memset(common.k.divisors.currentExp, 0, sizeof(common.k.divisors.currentExp));
  (void)memset(common.k.divisors.currentExpGray, 0, sizeof(common.k.divisors.currentExpGray));
  intToBigInteger(&common.k.divisors.divisor, 1);   // First divisor will be 1.
}

// Implementation of Heapsort as in Wikipedia.
static void siftDown(int start, int end, int** pptrBigInts)
{
  int root = start;
  while (((2 * root) + 1) <= end)
  {
    int child = (2 * root) + 1;
    int swap = root;
    // Compare a[swap] vs. a[child]
    if (IntArrayCompare(*(pptrBigInts + swap),
      *(pptrBigInts + child)) < 0)
    {        // a[swap] < a[child]
      swap = child;
    }
    if (((child + 1) <= end) &&
      (IntArrayCompare(*(pptrBigInts + swap),
        *(pptrBigInts + child + 1)) < 0))
    {             // a[swap] < a[child + 1]
      swap = child + 1;
    }
    if (swap == root)
    {
      return;
    }
    int* ptrTemp = *(pptrBigInts + root);
    *(pptrBigInts + root) = *(pptrBigInts + swap);
    *(pptrBigInts + swap) = ptrTemp;
    root = swap;
  }
}

static void heapify(int count, int** pptrBigInts)
{
  int start = (count - 2) / 2;
  while (start >= 0)
  {
    siftDown(start, count - 1, pptrBigInts);
    start--;
  }
}

void heapsort(int count, int** pptrBigInts)
{
  int end = count - 1;
  heapify(count, pptrBigInts);
  while (end > 0)
  {
    int* ptrTemp = *(pptrBigInts + end);
    *(pptrBigInts + end) = *pptrBigInts;
    *pptrBigInts = ptrTemp;
    end--;
    siftDown(0, end, pptrBigInts);
  }
}

// In order to minimize multiplications of big numbers, Gray code is used.
void showDivisors(void)
{
  bool showMoreDivisors = true;
  int nbrDivisors;
  int divisorNbr;
  int nbrExponents = astFactorsMod[0].multiplicity;
  int* ptrFoundDivisors = common.k.divisors.foundDivisors;
  char *ptrOutput = output;
  *ptrOutput = 'D';      // Indicate this output is the list of divisors.
  ptrOutput++;
  copyStr(&ptrOutput, lang ? "<p>Lista de divisores:</p><ul>" :
    "<p>List of divisors:</p><ul>");
  if ((tofactor.nbrLimbs == 1) && (tofactor.limbs[0].x <= 1))
  {
    if (tofactor.limbs[0].x == 0)
    {
      copyStr(&ptrOutput, lang? "<li>Cualquier número natural</li></ul>":
        "<li>Any natural number</li></ul>");
    }
    else
    {
      copyStr(&ptrOutput, "<li>1</li></ul>");
    }
    return;
  }
  if (nbrExponents > 50)
  {                      // Process only the first 50 exponents.
    nbrExponents = 50;
  }
  for (divisorNbr = 0; divisorNbr < 1000; divisorNbr++)
  {
    int exponentNbr;
    int arrLen;
    const struct sFactors* pstFactors;
    if (ptrFoundDivisors > &common.k.divisors.foundDivisors[900000])
    {
      break;             // Divisors are very large.
    }
    common.k.divisors.ptrFoundDivisors[divisorNbr] = ptrFoundDivisors;
    NumberLength = common.k.divisors.divisor.nbrLimbs;
    BigInteger2IntArray(ptrFoundDivisors, &common.k.divisors.divisor);
    arrLen = 1 + common.k.divisors.divisor.nbrLimbs;
    ptrFoundDivisors += arrLen;

    pstFactors = &astFactorsMod[1];
    // Find next divisor.
    for (exponentNbr = 0; exponentNbr < nbrExponents; exponentNbr++)
    {
      if (common.k.divisors.currentExpGray[exponentNbr] == 0)
      {  // Ascending exponents.
        if (common.k.divisors.currentExp[exponentNbr] <
          pstFactors->multiplicity)
        {
          common.k.divisors.currentExp[exponentNbr]++;
          NumberLength = *pstFactors->ptrFactor;
          IntArray2BigInteger(pstFactors->ptrFactor, &Tmp);
          (void)BigIntMultiply(&common.k.divisors.divisor, &Tmp,
            &common.k.divisors.divisor);
          NumberLength = common.k.divisors.divisor.nbrLimbs;
          BigInteger2IntArray(ptrFoundDivisors, &common.k.divisors.divisor);
          break;
        }
      }
      else
      {  // Descending exponents.
        if (common.k.divisors.currentExp[exponentNbr] > 0)
        {
          common.k.divisors.currentExp[exponentNbr]--;
          NumberLength = *pstFactors->ptrFactor;
          IntArray2BigInteger(pstFactors->ptrFactor, &Tmp);
          (void)BigIntDivide(&common.k.divisors.divisor, &Tmp,
            &common.k.divisors.divisor);
          NumberLength = common.k.divisors.divisor.nbrLimbs;
          BigInteger2IntArray(ptrFoundDivisors, &common.k.divisors.divisor);
          break;
        }
      }
      common.k.divisors.currentExpGray[exponentNbr] ^= 1;
      pstFactors++;
    }
    if (exponentNbr == nbrExponents)
    {           // All exponents processed.
      showMoreDivisors = false;
      divisorNbr++;
      break;
    }
  }
  // Sort factors.
  nbrDivisors = divisorNbr;
  heapsort(nbrDivisors, common.k.divisors.ptrFoundDivisors);
  // Generate output from sorted divisors.
  for (divisorNbr = 0; divisorNbr < nbrDivisors; divisorNbr++)
  {
    const int* ptrIntArray = common.k.divisors.ptrFoundDivisors[divisorNbr];
    int arrLen;
    NumberLength = *common.k.divisors.ptrFoundDivisors[divisorNbr];
    IntArray2BigInteger(common.k.divisors.ptrFoundDivisors[divisorNbr], &Tmp);
    copyStr(&ptrOutput, "<li>");
    arrLen = *ptrIntArray;
    if (hexadecimal)
    {
      Bin2Hex(&ptrOutput, (const limb *)ptrIntArray + 1, arrLen, groupLen);
    }
    else
    {
      Bin2Dec(&ptrOutput, (const limb*)ptrIntArray + 1, arrLen, groupLen);
    }
    copyStr(&ptrOutput, "</li>");
  }
  copyStr(&ptrOutput, "</ul>");
  if (showMoreDivisors)
  {
#ifdef __EMSCRIPTEN__
    copyStr(&ptrOutput, "<p><button type=\"button\" id=\"showdiv\">");
    copyStr(&ptrOutput, lang ? "Mostrar más divisores" : "Show more divisors");
    copyStr(&ptrOutput, "</button></p>");
#endif
    output[0] = 'E';    // Indicate button present.
  }
}

#ifndef _MSC_VER
EXTERNALIZE void doWork(void)
{
  int flags;
  char *ptrData = inputString;
  const char *ptrWebStorage;
  char *ptrKnownFactors;
#ifdef __EMSCRIPTEN__
  originalTenthSecond = tenths();
#endif
  if (*ptrData == 'C')
  {    // User pressed Continue button.
    ecmFrontText(NULL, false, NULL); // The 3rd parameter includes known factors.
#ifdef __EMSCRIPTEN__
    databack(output);
#endif
    return;
  }
  if (*ptrData == 'D')
  {    // User pressed Show Divisors button.
    showDivisors();
#ifdef __EMSCRIPTEN__
    databack(output);
#endif
    return;
  }
  if (*ptrData == 'S')
  {    // User pressed Show Sum of squares button.
    showSumTwoSquares();
#ifdef __EMSCRIPTEN__
    databack(output);
#endif
    return;
  }
  valuesProcessed = 0;
  groupLen = 0;
  while (*ptrData != ',')
  {
    groupLen = (groupLen * 10) + (*ptrData - '0');
    ptrData++;
  }
  ptrData++;             // Skip comma.
  flags = *ptrData;
  if (flags == '-')
  {
    ptrData++;
    flags = -*ptrData;
  }
#ifndef lang  
  lang = ((flags & 1)? true: false);
#endif
#ifdef __EMSCRIPTEN__
  useBlockly = false;
  doShowPrime = false;
  if ((flags & (-2)) == '8')
  {
    useBlockly = true;
  }
  else if ((flags & (-2)) == '4')
  {
    skipPrimality = true;
    flags = 2;           // Do factorization.
  }
  else if ((flags & (-2)) == '6')
  {
    doShowPrime = true;
    flags = 0;           // Do not perform factorization.
  }
  else
  {                      // No more cases.
  }
#endif
  ptrData += 2;          // Skip app number and second comma.
  fromFile = (*ptrData == '1');
  ptrData++;
  prettyprint = (*(ptrData + 1) == '1');
  if (fromFile)
  {
    prettyprint = false;
  }
  cunningham = (*(ptrData + 2) == '1');
  hexadecimal = (*(ptrData + 3) == '1');
  while (*ptrData != ',')
  {    // Loop that finds the end of configuration data.
    ptrData++;
  }
  ptrData++;    // Skip comma.
  ptrWebStorage = ptrData + (int)strlen(ptrData) + 1;
  if (useBlockly)
  {
    ptrKnownFactors = NULL;
  }
  else
  {
    ptrKnownFactors = findChar(ptrWebStorage, '=');
  }
  if (!prettyprint)
  {
    groupLen = -groupLen;  // Do not show number of digits.
  }
  if (ptrKnownFactors != NULL)
  {
    ptrKnownFactors++;
  }
  if ((flags & 0x80) && (ptrKnownFactors != NULL))
  {
    flags = 2;    // Do factorization.
  }
  if (useBlockly)
  {
    fromBlockly(ptrData);
  }
  else
  {               // Do factorization if second parameter is true.
                  // The 3rd parameter includes known factors.
    ecmFrontText(ptrData, (flags & 2) != 0, ptrKnownFactors);
  }
#ifdef __EMSCRIPTEN__
  databack(output);
#endif
}
#endif
#endif
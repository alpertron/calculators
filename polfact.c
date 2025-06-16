//
// This file is part of Alpertron Calculators.
//
// Copyright 2015-2024 Dario Alejandro Alpern
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
#include <stdlib.h>
#include "bignbr.h"
#include "expression.h"
#include "highlevel.h"
#include "polynomial.h"
#include "showtime.h"
#include "rootseq.h"
#ifdef FACTORIZATION_APP
#include "factor.h"
#endif
extern int primeEisenstein;
extern char* ptrOutput;
extern int eqNbr;
int grpLen;
bool teachMod;

// Sort factors on ascending degree, and then by coefficient.
void SortFactors(const BigInteger *modulus)
{
  struct sFactorInfo *pstFactorInfo;
  struct sFactorInfo *pstFactorInfo2;
  struct sFactorInfo stFactorInfoTemp;
  int currentDegree;
  int index;
  const int *ptrValue1;
  const int *ptrValue2;
  int nbrLimbs = modulus->nbrLimbs + 1;
  pstFactorInfo = factorInfo;
  for (int nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    pstFactorInfo2 = pstFactorInfo + 1;
    for (int nbrFactor2 = nbrFactor + 1; nbrFactor2 < nbrFactorsFound; nbrFactor2++)
    {
      if (pstFactorInfo->degree > pstFactorInfo2->degree)
      {
        stFactorInfoTemp = *pstFactorInfo;
        *pstFactorInfo = *pstFactorInfo2;
        *pstFactorInfo2 = stFactorInfoTemp;
      }
      else if (pstFactorInfo->degree == pstFactorInfo2->degree)
      {
        index = 0;
        ptrValue1 = pstFactorInfo->ptrPolyLifted + (pstFactorInfo->degree * nbrLimbs);
        ptrValue2 = pstFactorInfo2->ptrPolyLifted + (pstFactorInfo->degree * nbrLimbs);
        for (currentDegree = pstFactorInfo->degree - 1; currentDegree >= 0; currentDegree--)
        {
          ptrValue1 -= nbrLimbs;
          ptrValue2 -= nbrLimbs;
          if (*ptrValue1 > *ptrValue2)
          {        // Coefficient of first polynomial is greater than coeff of second.
            break;
          }
          if (*ptrValue1 == *ptrValue2)
          {        // Number of limbs of coefficients match.
            for (index = *ptrValue1; index > 0; index--)
            {
              if (*(ptrValue1 + index) != *(ptrValue2 + index))
              {
                break;
              }
            }
            if (index > 0)
            {
              break;
            }
          }
        }
        if ((currentDegree >= 0) && (*(ptrValue1 + index) > *(ptrValue2 + index)))
        {
          stFactorInfoTemp = *pstFactorInfo;
          *pstFactorInfo = *pstFactorInfo2;
          *pstFactorInfo2 = stFactorInfoTemp;
        }
      }
      else
      {              // Nothing to do.
      }
      pstFactorInfo2++;
    }
    pstFactorInfo++;
  }
}

static int FactorPolynomial(void)
{
  // At this moment the array "common.poly.values" contains the polynomial.
  if (onlyEvaluate)
  {
    if (!modulusIsZero)
    {
      OrigPolyFromMontgomeryToStandard();
    }
    return EXPR_OK;
  }
  // Generate polynomial mod prime.
#ifdef __EMSCRIPTEN__
  originalTenthSecond = tenths();
#endif
  if (modulusIsZero)
  {
    teachMod = false;
    return FactorPolyOverIntegers();
  }
  // Back up input polynomial.
  common.poly.polyBackup[0] = common.poly.values[0];
  (void)CopyPolynomial(&common.poly.polyBackup[1], &common.poly.values[1], (common.poly.values[0] >= 0) ? common.poly.values[0] : 0);
  // Input is in Montgomery notation.
  teachMod = teach;
  return FactorModularPolynomial(true, false);
}

void polyFactText(const char *modText, const char *polyText, int groupLength)
{
  enum eExprErr rc;
  int expon = 0;
  bool isFraction;
  grpLen = groupLength;
  eqNbr = 0;
  rc = ComputeExpression(modText, &powerMod);
  modulusIsZero = false;
  if (pretty == PRETTY_PRINT)
  {
    ptrTimes = "&#8290;";
    ptrMinus = "&minus;";
    ptrPlusMinus = " &pm; ";
  }
  else
  {
    ptrTimes = (pretty == TEX ? "" : "*");
    ptrMinus = "-";
    ptrPlusMinus = (pretty == TEX ? " \\pm " : " &pm; ");
  }
  if (rc == EXPR_OK)
  {
    if (BigIntIsZero(&powerMod))
    {
      modulusIsZero = true;
    }
    else if (powerMod.sign == SIGN_NEGATIVE)
    {                  // Negative is not greater than 1.
      rc = EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE;
    }
    else if ((powerMod.nbrLimbs == 1) && (powerMod.limbs[0].x < 2))
    {                  // Positive number is less 2.
      rc = EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE;
    }
    else
    {              // Nothing to do.
    }
  }
  if ((rc == EXPR_OK) && (!modulusIsZero))
  {
    expon = PowerCheck(&powerMod, &primeMod);
#if FACTORIZATION_APP
    if (BpswPrimalityTest(&primeMod, NULL) != 0)
#else
    if (BpswPrimalityTest(&primeMod) != 0)
#endif
    {    // Number is composite
      rc = EXPR_MODULUS_MUST_BE_PRIME_EXP;
    }
  }
  output[0] = '2';
  ptrOutput = &output[1];
  if (rc == EXPR_OK)
  {
    rc = ComputePolynomial(polyText, expon);
    if (rc == EXPR_OK)
    {
      isFraction = true;
      if ((common.poly.denom[0] == 0) && (((common.poly.denom[1] == 1) && (common.poly.denom[2] == 1)) || !modulusIsZero))
      {    // If modulus is zero: denominator is one.
           // If modulus is not zero: degree is zero.
        isFraction = false;
        copyStr(&ptrOutput, lang ? "<h2>Polinomio ingresado</h2>" :
          "<h2>Your polynomial</h2>");
      }
      else
      {
        copyStr(&ptrOutput, lang ? "<h2>Fracción de polinomios</h2>" :
          "<h2>Your polynomial fraction</h2>");
      }
      if (onlyEvaluate)
      {
        copyStr(&ptrOutput, "<p>");
      }
      else
      {
        copyStr(&ptrOutput, "<p id=\"pol\">");
      }
      outputOriginalPolynomial(&ptrOutput, groupLength);
      copyStr(&ptrOutput, "</p>");
      if (!onlyEvaluate)
      {
        if (isFraction)
        {
          copyStr(&ptrOutput, lang ? "<h2>Factores irreducibles del polinomio numerador</h2>" :
            "<h2>Irreducible numerator factors</h2>");
        }
        else
        {
          copyStr(&ptrOutput, lang ? "<h2>Factores irreducibles del polinomio</h2>" :
            "<h2>Irreducible polynomial factors</h2>");
        }
        rc = FactorPolynomial();
      }
    }
  }
  if (rc != EXPR_OK)
  {
    textErrorPol(&ptrOutput, rc);
  }
  else
  {
    degree = common.poly.values[0];
    if (!onlyEvaluate)
    {
      int nbrFactor;
      struct sFactorInfo* pstFactorInfo;
      if (modulusIsZero)
      {
        pstFactorInfo = factorInfoInteger;
      }
      else
      {
        pstFactorInfo = factorInfo;
        // Get leading coefficient if using modular arithmetic.
        const int* ptrCoeff = &common.poly.values[1];
        for (int currDegree = 0; currDegree < degree; currDegree++)
        {
          ptrCoeff += numLimbs(ptrCoeff) + 1;
        }
        IntArray2BigInteger(ptrCoeff, &operand5);
      }
      if ((nbrFactorsFound == 0) || ((nbrFactorsFound == 1) &&
        (pstFactorInfo->multiplicity == 1) && BigIntIsOne(&operand5)))
      {
        copyStr(&ptrOutput, lang ? "<p>El polinomio es irreducible" : "<p>The polynomial is irreducible");
        if (modulusIsZero && (primeEisenstein != 0))
        {
          copyStr(&ptrOutput, lang ? " debido al criterio de Eisenstein (primo = " :
            " because of Eisenstein's criterion (prime = ");
          int2dec(&ptrOutput, primeEisenstein);
          *ptrOutput = ')';
          ptrOutput++;
        }
        copyStr(&ptrOutput, "</p>");
      }
      else
      {   // Get number of factors including multiplicity.
        int totalFactors = 0;
        for (int ctr = 0; ctr < nbrFactorsFound; ctr++)
        {
          totalFactors += pstFactorInfo->multiplicity;
          pstFactorInfo++;
        }
        if (!BigIntIsOne(&operand5))
        {    // Add factor of degree zero if it is not one.
          totalFactors++;
        }
        copyStr(&ptrOutput, lang ? "<p>Los " : "<p>The ");
        int2dec(&ptrOutput, totalFactors);
        copyStr(&ptrOutput, lang ? " factores son:</p>" : " factors are:</p>");
        if (modulusIsZero)
        {
          pstFactorInfo = factorInfoInteger;
        }
        else
        {
          pstFactorInfo = factorInfo;
        }
        // Output factors
        showText("<ul>");
        if (pretty == TEX)
        {
          showText("<li>\\begin{array}{l}</li>");
        }
        if (!BigIntIsOne(&operand5) || (nbrFactorsFound == 0))
        {     // Leading coefficient is not 1 or degree is zero.
          showText("<li>");
          if (operand5.sign == SIGN_NEGATIVE)
          {
            copyStr(&ptrOutput, " &minus;");
          }
          if (pretty == PRETTY_PRINT)
          {          // Show number of digits if there are more than 30.
            Bin2Dec(&ptrOutput, operand5.limbs, operand1.nbrLimbs, groupLength);
          }
          else
          {         // Do not show number of digits.
            Bin2Dec(&ptrOutput, operand5.limbs, operand1.nbrLimbs, -groupLength);
          }
          showText("</li>");
        }
        for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
        {
          showText("<li>");
          if (pretty == TEX)
          {
            showText("\\bullet\\,\\,");
          }
          outputPolynomialFactor(&ptrOutput, groupLength, pstFactorInfo);
          if (pretty == TEX)
          {
            showText("\\\\");
          }
          showText("</li>");
          pstFactorInfo++;
        }
        if (pretty == TEX)
        {
          copyStr(&ptrOutput, "<li>\\end{array}</li>");
        }
        showText("</ul>");
      }
      if (modulusIsZero)
      {
        copyStr(&ptrOutput, lang ? "<h2>Raíces</h2>" : "<h2>Roots</h2>");
        if (degree > 1)
        {
          copyStr(&ptrOutput, lang ? "<p>Las " : "<p>The ");
          int2dec(&ptrOutput, degree);
          copyStr(&ptrOutput, lang ? " raíces son:</p>" : " roots are:</p>");
        }
        copyStr(&ptrOutput, "<ul>");
        if (pretty == TEX)
        {
          copyStr(&ptrOutput, "<li>\\begin{array}{l}</li>");
        }
        indexRoot = 1;
        pstFactorInfo = factorInfoInteger;
        for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
        {
          if (nbrFactorsFound == 1)
          {    // Do not show polynomial factor number.
            getRootsPolynomial(-1, &ptrOutput, pstFactorInfo, groupLength);
          }
          else
          {
            getRootsPolynomial(nbrFactor, &ptrOutput, pstFactorInfo, groupLength);
          }
          pstFactorInfo++;
        }
        if (pretty == TEX)
        {
          copyStr(&ptrOutput, "<li>\\end{array}</li>");
        }
        copyStr(&ptrOutput, "</ul>");
      }
      // Show time only when factoring, not when just evaluating polynomial.
      copyStr(&ptrOutput, "<p>");
      showElapsedTime(&ptrOutput);
      copyStr(&ptrOutput, "</p>");
    }
  }
  copyStr(&ptrOutput, "<p>");
  copyStr(&ptrOutput, lang ? COPYRIGHT_SPANISH: COPYRIGHT_ENGLISH);
  copyStr(&ptrOutput, "</p>");
}

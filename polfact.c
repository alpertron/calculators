/*
This file is part of Alpertron Calculators.

Copyright 2015 Dario Alejandro Alpern

Alpertron Calculators is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Alpertron Calculators is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Alpertron Calculators.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "bignbr.h"
#include "highlevel.h"
#include "polynomial.h"
#include "showtime.h"
#include "rootseq.h"

int attemptNbr;
#ifdef __EMSCRIPTEN__
char *ptrPercentageOutput;
#endif

static char outputText[20000];
extern int poly4[1000000];

// Perform distinct degree factorization
static void DistinctDegreeFactorization(int polyDegree)
{
  struct sFactorInfo *pstFactorInfo, *pstNewFactorInfo;
  int nbrLimbs = primeMod.nbrLimbs + 1;
  int currentDegree;
  int nbrFactor, degreeMin, degreeGcd;
  int *ptrPolyToFactor, *ptrValue1;
  // Set poly1 to x.
  memset(poly1, 0, nbrLimbs*(polyDegree + 1)*sizeof(int));
  for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {
    poly1[currentDegree*nbrLimbs] = 1;
  }
  SetNumberToOne(&poly1[nbrLimbs]);
  pstFactorInfo = factorInfo;
  for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    if (pstFactorInfo->degree < 2 || pstFactorInfo->expectedDegree != 0)
    {             // Polynomial is completely factored. Try next one.
      if (pstFactorInfo->expectedDegree == 0)
      {
        pstFactorInfo->expectedDegree = pstFactorInfo->degree;
      }
      pstFactorInfo++;
      continue;
    }
    ptrPolyToFactor = pstFactorInfo->ptr;
    polyDegree = pstFactorInfo->degree;
    GetPolyInvParm(polyDegree, ptrPolyToFactor);
    // For each loop, raise this polynomial to the primeth power and 
    // then compute the gcd between the polynomial to be factored and
    // the computed polynomial less x. If the degree of GCD is > 0, then the
    // GCD is the product of all factors of degree indicated by currentDegree.
    for (currentDegree = 1; currentDegree * 2 <= polyDegree; currentDegree++)
    {
#ifdef __EMSCRIPTEN__
      int elapsedTime = (int)(tenths() - originalTenthSecond);
      if (elapsedTime / 10 != oldTimeElapsed / 10)
      {
        char *ptrOutput = outputText;
        oldTimeElapsed = elapsedTime;
        if (lang)
        {
          strcpy(ptrOutput, "1<p>Factorización de distintos grados: buscando factores de grado ");
          ptrOutput += strlen(ptrOutput);
          int2dec(&ptrOutput, currentDegree);
          strcpy(ptrOutput, " (máx.  ");
          ptrOutput += strlen(ptrOutput);
          int2dec(&ptrOutput, (polyDegree + 1) / 2);
          strcpy(ptrOutput, ") del factor número ");
          ptrOutput += strlen(ptrOutput);
          int2dec(&ptrOutput, nbrFactor + 1);
          strcpy(ptrOutput, " de ");
        }
        else
        {
          strcpy(ptrOutput, "1<p>Distinct degree factorization: searching for factors of degree ");
          ptrOutput += strlen(ptrOutput);
          int2dec(&ptrOutput, currentDegree);
          strcpy(ptrOutput, " (max.  ");
          ptrOutput += strlen(ptrOutput);
          int2dec(&ptrOutput, (polyDegree + 1) / 2);
          strcpy(ptrOutput, ") of factor number ");
          ptrOutput += strlen(ptrOutput);
          int2dec(&ptrOutput, nbrFactor + 1);
          strcpy(ptrOutput, " of ");
        }
        ptrOutput += strlen(ptrOutput);
        int2dec(&ptrOutput, nbrFactorsFound);
        strcpy(ptrOutput, lang ? ".</p><p>Transcurrió " : ".</p><p>Time elapsed: ");
        ptrOutput += strlen(ptrOutput);
        GetDHMS(&ptrOutput, elapsedTime / 10);
        strcpy(ptrOutput, "</p>");
        databack(outputText);
      }
#endif
      // Copy polynomial to factor to poly3 and set leading coefficient to 1.
      ptrValue1 = &poly3[polyDegree*nbrLimbs];
      memcpy(poly3, ptrPolyToFactor, (ptrValue1 - &poly3[0])*sizeof(int));
      SetNumberToOne(ptrValue1);  // Set leading coefficient to 1.
      powerPolynomial(poly1, poly3,  // Base and polynomial modulus.
        polyDegree, &primeMod,       // Degree of polynomials and exponent.
        poly2, NULL,                 // Power and pointer to callback.
        0, 1);
      memcpy(poly1, poly2, polyDegree*nbrLimbs * sizeof(int));
      // Subtract x.
      IntArray2BigInteger(&poly2[nbrLimbs], &operand1);
      memcpy(operand2.limbs, MontgomeryMultR1, NumberLength*sizeof(limb));
      operand2.nbrLimbs = NumberLengthR1;
      SubtBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
      BigInteger2IntArray(&poly2[nbrLimbs], &operand1);
      // Perform Gcd.
      degreeMin = getDegreePoly(poly2, polyDegree - 1);
      PolyModularGcd(poly3, polyDegree, poly2, degreeMin, poly4, &degreeGcd);
      if (degreeGcd == polyDegree)
      {
        pstFactorInfo->expectedDegree = currentDegree;
        polyDegree = 0;
      }
      else if (degreeGcd > 0)
      {         // Non-trivial factor of polynomial has been found.
                // Divide polynomial by GCD. Put the GCD in the first limbs
                // and the quotient in the last limbs.
        ptrValue1 = &poly4[degreeGcd*nbrLimbs];
        SetNumberToOne(ptrValue1);
        DividePolynomial(poly3, polyDegree, poly4, degreeGcd, poly2);
        // Quotient located in poly2.
        pstNewFactorInfo = &factorInfo[nbrFactorsFound++];
        pstNewFactorInfo->ptr = ptrPolyToFactor;
        pstNewFactorInfo->degree = degreeGcd;
        pstNewFactorInfo->multiplicity = pstFactorInfo->multiplicity;
        pstNewFactorInfo->expectedDegree = currentDegree;
        pstFactorInfo->degree = polyDegree - degreeGcd;
        pstFactorInfo->ptr = &ptrPolyToFactor[degreeGcd*nbrLimbs];
        memcpy(ptrPolyToFactor, poly4, degreeGcd*nbrLimbs*sizeof(int));
        memcpy(pstFactorInfo->ptr, poly2, (polyDegree - degreeGcd + 1)*nbrLimbs*sizeof(int));
        polyDegree -= degreeGcd;
        ptrPolyToFactor += degreeGcd*nbrLimbs;
        // Replace poly1 by poly1 mod ptrPolyToFactor
        DividePolynomial(poly1, polyDegree + degreeGcd - 1, ptrPolyToFactor, polyDegree, poly2);
        if (polyDegree > 0)
        {
          GetPolyInvParm(polyDegree, ptrPolyToFactor);
        }
      }
    }
    if (polyDegree > 0)
    {
      pstFactorInfo->expectedDegree = polyDegree;
    }
  }
}

static void percentageCallback(int percentage)
{
#ifdef __EMSCRIPTEN__
  int elapsedTime = (int)(tenths() - originalTenthSecond);
  char *ptrOutput = ptrPercentageOutput;
  if (elapsedTime / 10 != oldTimeElapsed / 10)
  {
    oldTimeElapsed = elapsedTime;
    int2dec(&ptrOutput, percentage);
    if (lang)
    {
      strcpy(ptrOutput, "% del ");
      ptrOutput += strlen(ptrOutput);
      int2dec(&ptrOutput, attemptNbr);
      strcpy(ptrOutput, ".º intento");
      ptrOutput += strlen(ptrOutput);
    }
    else
    {
      strcpy(ptrOutput, "% of attempt #");
      ptrOutput += strlen(ptrOutput);
      int2dec(&ptrOutput, attemptNbr);
    }
    strcpy(ptrOutput, lang ? ".</p><p>Transcurrió " : ".</p><p>Time elapsed: ");
    ptrOutput += strlen(ptrOutput);
    GetDHMS(&ptrOutput, elapsedTime / 10);
    strcpy(ptrOutput, "</p>");
    databack(outputText);
  }
#else
  (void)percentage;
#endif
}


// Perform Cantor-Zassenhaus algorithm to factor polynomials of the same degree.
void SameDegreeFactorization(void)
{
  struct sFactorInfo *pstFactorInfo = factorInfo;
  struct sFactorInfo *pstNewFactorInfo;
  int *ptrValue1, *ptrPolyToFactor;
  int nbrFactor, currentDegree, index, degreeGcd;
  int primeInt = (int)primeMod.limbs[0].x;
  int nbrLimbs = primeMod.nbrLimbs + 1;
  int polyNbr = 1;
  int isCharacteristic2 = primeMod.nbrLimbs == 1 && primeMod.limbs[0].x == 2;
  for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    int polyDegree = pstFactorInfo->degree;
    if (polyDegree < 2 || polyDegree == pstFactorInfo->expectedDegree ||
      pstFactorInfo->expectedDegree == 0)
    {             // Polynomial is completely factored. Try next one.
      pstFactorInfo++;
      continue;
    }
    if (isCharacteristic2 == 0)
    { // If prime is not 2,
      // Calculate operand2 <- (prime^degree-1)/2
      // Use operand1 as temporary variable to store the exponent.
      operand1.limbs[0].x = pstFactorInfo->expectedDegree & MAX_VALUE_LIMB;
      operand1.nbrLimbs = 1;
    }
    ptrPolyToFactor = pstFactorInfo->ptr;
    for (attemptNbr = 1;; attemptNbr++, polyNbr++)
    {
#ifdef __EMSCRIPTEN__
      char *ptrOutput = outputText;
      if (lang)
      {
        strcpy(ptrOutput, "1<p>Factorización del mismo grado: buscando ");
        ptrOutput += strlen(ptrOutput);
        int2dec(&ptrOutput, polyDegree / pstFactorInfo->expectedDegree);
        strcpy(ptrOutput, " factores de grado ");
      }
      else
      {
        strcpy(ptrOutput, "1<p>Equal degree factorization: searching for ");
        ptrOutput += strlen(ptrOutput);
        int2dec(&ptrOutput, polyDegree / pstFactorInfo->expectedDegree);
        strcpy(ptrOutput, " factors of degree ");
      }
      ptrOutput += strlen(ptrOutput);
      int2dec(&ptrOutput, pstFactorInfo->expectedDegree);
      strcpy(ptrOutput, ".</p><p>");
      ptrPercentageOutput = ptrOutput + strlen(ptrOutput);
#endif
      // Copy polynomial to factor to poly3 and set leading coefficient to 1.
      // All operations below will be done modulo this polynomial.
      ptrValue1 = &poly3[pstFactorInfo->degree*nbrLimbs];
      memcpy(poly3, ptrPolyToFactor, (ptrValue1 - &poly3[0])*sizeof(int));
      SetNumberToOne(ptrValue1);  // Set leading coefficient to 1.
      if (attemptNbr == 1)
      {
        GetPolyInvParm(polyDegree, poly3);
      }
      // Initialize polynomial poly1 with different values of coefficients
      // in different iterations.
      ptrValue1 = poly1;
      if (nbrLimbs > 2)
      {    // Coefficient can be any number.
        *ptrValue1 = 1;
        *(ptrValue1 + 1) = polyNbr;
          ptrValue1 += nbrLimbs;
        *ptrValue1 = 1;
        *(ptrValue1 + 1) = 1;
        ptrValue1 += nbrLimbs;
        for (currentDegree = 2; currentDegree < polyDegree; currentDegree++)
        {
          *ptrValue1 = 1;
          *(ptrValue1 + 1) = 0;
          ptrValue1 += nbrLimbs;
        }
      }
      else
      {   // Coefficient range can be from 0 to prime-1.
        int polyCoeff = polyNbr;
        for (currentDegree = 0; polyCoeff != 0; currentDegree++)
        {
          *ptrValue1 = 1;
          *(ptrValue1 + 1) = polyCoeff % primeInt;
          polyCoeff /= primeInt;
          ptrValue1 += nbrLimbs;
        }
        *ptrValue1 = 1;
        *(ptrValue1 + 1) = 1;
        ptrValue1 += nbrLimbs;
        for (currentDegree++; currentDegree < polyDegree; currentDegree++)
        {
          *ptrValue1 = 1;
          *(ptrValue1 + 1) = 0;
          ptrValue1 += nbrLimbs;
        }
      }
      if (isCharacteristic2 == 0)
      { // If prime is not 2: compute base^((p^d-1)/2).
        int degreeFactor = pstFactorInfo->expectedDegree;
        CopyBigInt(&operand4, &primeMod);
        subtractdivide(&operand4, 1, 2);  // operand4 <- exponent = (p-1)/2.
        // Start by raising base to power (p-1)/2.
        powerPolynomial(poly1, poly3,   // Base and polynomial modulus.
          polyDegree, &operand4,        // Degree of polynomials and exponent.
          poly2, percentageCallback,   // Power and pointer to callback.
          0, degreeFactor);
        // Save base^((p-1)/2) on poly4.
        CopyPolynomialFixedCoeffSize(poly4, poly2, polyDegree, nbrLimbs);
        for (currentDegree = 1; currentDegree < degreeFactor; currentDegree++)
        {
          // Square polynomial and multiply by base to get base^(p^(q-1)).
          multUsingInvPolynomial(poly2, poly2, // Multiplicands.
            poly2, polyDegree,                 // Product and degree of poly.
            poly3);                            // Polynomial modulus.
          multUsingInvPolynomial(poly2, poly1, // Multiplicands.
            poly1, polyDegree,                 // Product and degree of poly.
            poly3);                            // Polynomial modulus.
          CopyBigInt(&operand4, &primeMod);
          subtractdivide(&operand4, 1, 2);  // operand4 <- exponent = (p-1)/2.
            // Raise previous power to exponent (p-1)/2 to get
          // base^(p^(q-1)*(p-1)/2))
          powerPolynomial(poly1, poly3,   // Base and polynomial modulus.
            polyDegree, &operand4,        // Degree of polynomials and exponent.
            poly2, percentageCallback,    // Power and pointer to callback.
            currentDegree, degreeFactor);
          // Multiply base^(p^(q-1)*(p-1)/2) * base^(p^(q-2)*(p-1)/2)
          multUsingInvPolynomial(poly2, poly4, // Multiplicands.
            poly4, polyDegree,                 // Product and degree of poly.
            poly3);                            // Polynomial modulus.
        }
        CopyPolynomialFixedCoeffSize(poly2, poly4, polyDegree, nbrLimbs);
        // Subtract 1.
        IntArray2BigInteger(&poly2[0], &operand1);
        SubtBigNbrMod(operand1.limbs, MontgomeryMultR1, operand1.limbs);
        BigInteger2IntArray(&poly2[0], &operand1);
      }
      else
      { // If prime is 2, Compute poly2 = T+T^2+T^4+...+T^2^(d-1) mod f(x)
        // where T is the random polynomial.
        // Z <- T mod F.
        memcpy(poly2, poly1, polyDegree*nbrLimbs*sizeof(int));
        for (currentDegree = 1; currentDegree < pstFactorInfo->expectedDegree; currentDegree++)
        {
          multPolynomialModPoly(poly1, poly1, poly1, polyDegree, poly3);
          for (index = 0; index < polyDegree; index++)
          {
            IntArray2BigInteger(&poly1[index*nbrLimbs], &operand1);
            IntArray2BigInteger(&poly2[index*nbrLimbs], &operand2);
            AddBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
            BigInteger2IntArray(&poly2[index*nbrLimbs], &operand1);
          }
        }
      }
      PolyModularGcd(poly3, polyDegree, poly2, getDegreePoly(poly2, polyDegree - 1), poly4, &degreeGcd);
      if (degreeGcd != 0 && degreeGcd != polyDegree)
      {   // Non-trivial factor found.
        ptrValue1 = &poly4[degreeGcd*nbrLimbs];
        SetNumberToOne(ptrValue1);
        DividePolynomial(poly3, polyDegree, poly4, degreeGcd, poly2);
        // Quotient located in poly2.
        pstNewFactorInfo = &factorInfo[nbrFactorsFound++];
        pstNewFactorInfo->ptr = &ptrPolyToFactor[degreeGcd*nbrLimbs];
        pstNewFactorInfo->degree = polyDegree - degreeGcd;
        pstNewFactorInfo->multiplicity = pstFactorInfo->multiplicity;
        pstNewFactorInfo->expectedDegree = pstFactorInfo->expectedDegree;
        pstFactorInfo->degree = degreeGcd;
        memcpy(ptrPolyToFactor, poly4, degreeGcd*nbrLimbs*sizeof(int));
        memcpy(pstNewFactorInfo->ptr, poly2, (polyDegree - degreeGcd)*nbrLimbs*sizeof(int));
        polyDegree = degreeGcd;
        attemptNbr = 0;
        if (pstFactorInfo->expectedDegree == pstFactorInfo->degree)
        {
          break;
        }
        GetPolyInvParm(polyDegree, ptrPolyToFactor);
      }
    }
    pstFactorInfo++;
  }
}

// Sort factors on ascending degree, and then by coefficient.
static void SortFactors(BigInteger *modulus)
{
  struct sFactorInfo *pstFactorInfo, *pstFactorInfo2;
  struct sFactorInfo stFactorInfoTemp;
  int currentDegree, nbrFactor, nbrFactor2, index;
  int *ptrValue1, *ptrValue2;
  int nbrLimbs = modulus->nbrLimbs + 1;
  pstFactorInfo = factorInfo;
  for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    pstFactorInfo2 = pstFactorInfo + 1;
    for (nbrFactor2 = nbrFactor + 1; nbrFactor2 < nbrFactorsFound; nbrFactor2++)
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
        ptrValue1 = pstFactorInfo->ptrPolyLifted + pstFactorInfo->degree * nbrLimbs;
        ptrValue2 = pstFactorInfo2->ptrPolyLifted + pstFactorInfo->degree * nbrLimbs;
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
        if (currentDegree >= 0 && *(ptrValue1 + index) > *(ptrValue2 + index))
        {
          stFactorInfoTemp = *pstFactorInfo;
          *pstFactorInfo = *pstFactorInfo2;
          *pstFactorInfo2 = stFactorInfoTemp;
        }
      }
      pstFactorInfo2++;
    }
    pstFactorInfo++;
  }
}

// Input: values = degree, coefficient degree 0, coefficient degree 1, etc.
// Output: factorInfo = structure that holds the factors.
int FactorModularPolynomial(int inputMontgomery)
{
  int nbrFactor;
  struct sFactorInfo* ptrFactorInfo;
  int currentDegree, rc;
  int *ptrValue1;
  int nbrLimbsPrime = primeMod.nbrLimbs + 1; // Add 1 for length;
  degree = values[0];
  ptrValue1 = &values[1];
  for (currentDegree = 0; currentDegree <= degree; currentDegree++)
  {
    NumberLength = numLimbs(ptrValue1);
    IntArray2BigInteger(ptrValue1, &operand1);
    NumberLength = powerMod.nbrLimbs;
    if (inputMontgomery)
    {
      // Convert from Montgomery to standard notation.
      operand2.limbs[0].x = 1;
      if (NumberLength > 1)
      {
        memset(&operand2.limbs[1], 0, (NumberLength - 1) * sizeof(limb));
      }
      modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    }
    rc = BigIntRemainder(&operand1, &primeMod, &operand1);
    if (rc != EXPR_OK)
    {
      return rc;
    }
    NumberLength = primeMod.nbrLimbs;
    BigInteger2IntArray(&valuesPrime[currentDegree*nbrLimbsPrime], &operand1);
    ptrValue1 += 1 + numLimbs(ptrValue1);
  }
  if (operand1.nbrLimbs == 1 && operand1.limbs[0].x == 0)
  {
    return EXPR_LEADING_COFF_MULTIPLE_OF_PRIME;
  }
  memcpy(&TestNbr, primeMod.limbs, primeMod.nbrLimbs * sizeof(limb));
  NumberLength = primeMod.nbrLimbs;
  TestNbr[NumberLength].x = 0;
  GetMontgomeryParms(primeMod.nbrLimbs);
  // Convert polynomial mod prime to monic (leading coefficient must be 1).
  ConvertToMonic(valuesPrime, degree);
  // Perform square free factorization.
  ptrOrigPoly = valuesPrime;
  degreeOrigPoly = degree;
  nbrFactorsFound = 0;
  if (degree != 0)
  {
    SquareFreeFactorization(degree, valuesPrime, 1);
    DistinctDegreeFactorization(degree);
    if (inputMontgomery)
    {    // Do not perform same degree factorization if only counting
         // the number of modular factors of input polynomial.
      SameDegreeFactorization();
    }
  }
  if (inputMontgomery == FALSE)
  {
    return EXPR_OK;
  }
  // Convert original polynomial from Montgomery to standard notation.
  ptrFactorInfo = factorInfo;
  for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    polyToStandardNotation(ptrFactorInfo->ptr, ptrFactorInfo->degree);
    ptrFactorInfo++;
  }
  OrigPolyFromMontgomeryToStandard();
  rc = HenselLifting(factorInfo, 0);
  if (rc != EXPR_OK)
  {
    return rc;
  }
  SortFactors(rc != 0 ? &primeMod : &powerMod);
  return rc;
}

static int FactorPolynomial(char *input, int expo)
{
  int rc = ComputePolynomial(input, expo);
  if (rc != 0)
  {
    return rc;
  }
  // At this moment the array "values" contains the polynomial.
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
    return FactorPolyOverIntegers();
  }
  return FactorModularPolynomial(TRUE);   // Input is in Montgomery notation.
}

void polyFactText(char *modText, char *polyText, int groupLength)
{
  char *ptrOutput;
  enum eExprErr rc;
  int expon = 0;
  rc = ComputeExpression(modText, 1, &powerMod);
  modulusIsZero = 0;
  if (rc == EXPR_OK)
  {
    if (powerMod.nbrLimbs == 1 && powerMod.limbs[0].x == 0)
    {
      modulusIsZero = 1;
    }
    else if (powerMod.sign == SIGN_NEGATIVE || (powerMod.nbrLimbs == 1 && powerMod.limbs[0].x < 2))
    {
      rc = EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE;
    }
  }
  if (rc == EXPR_OK && !modulusIsZero)
  {
    expon = PowerCheck(&powerMod, &primeMod);
#if FACTORIZATION_APP
    if (BpswPrimalityTest(&primeMod, NULL))
#else
    if (BpswPrimalityTest(&primeMod))
#endif
    {    // Number is composite
      rc = EXPR_MODULUS_MUST_BE_PRIME_EXP;
    }
  }
  if (rc == EXPR_OK)
  {
    rc = FactorPolynomial(polyText, expon);
  }
  output[0] = '2';
  ptrOutput = &output[1];
  if (rc != EXPR_OK)
  {
    textErrorPol(ptrOutput, rc);
    ptrOutput += strlen(ptrOutput);
  }
  else
  {
    strcpy(ptrOutput, lang ? "<h2>Polinomio ingresado</h2>" : "<h2>Your polynomial</h2>");
    ptrOutput += strlen(ptrOutput);
    if (onlyEvaluate)
    {
      strcpy(ptrOutput, "<p>");
    }
    else
    {
      strcpy(ptrOutput, "<p id=\"pol\">");
    }
    ptrOutput += strlen(ptrOutput);
    outputOriginalPolynomial(ptrOutput, groupLength);
    ptrOutput += strlen(ptrOutput);
    strcpy(ptrOutput, "</p>");
    ptrOutput += strlen(ptrOutput);
    if (onlyEvaluate == 0)
    {
      int nbrFactor;
      int nbrLimbs = powerMod.nbrLimbs + 1;
      struct sFactorInfo* pstFactorInfo;
      if (modulusIsZero)
      {
        pstFactorInfo = factorInfoInteger;
      }
      else
      {
        pstFactorInfo = factorInfo;
      }
      strcpy(ptrOutput, lang? "<h2>Factores irreducibles del polinomio</h2>": "<h2>Irreducible polynomial factors</h2>");
      ptrOutput += strlen(ptrOutput);
      if (nbrFactorsFound == 0 || (nbrFactorsFound == 1 && pstFactorInfo->multiplicity == 1))
      {
        strcpy(ptrOutput, lang ? "<p>El polinomio es irreducible</p>" : "<p>The polynomial is irreducible</p>");
        ptrOutput += strlen(ptrOutput);
      }
      else
      {   // Get number of factors including multiplicity.
        int totalFactors = 0;
        int ctr;
        for (ctr = 0; ctr < nbrFactorsFound; ctr++)
        {
          totalFactors += pstFactorInfo->multiplicity;
          pstFactorInfo++;
        }
        strcpy(ptrOutput, lang ? "Los " : "The ");
        ptrOutput += strlen(ptrOutput);
        int2dec(&ptrOutput, totalFactors);
        strcpy(ptrOutput, lang ? " factores son:</p>" : " factors are:</p>");
        ptrOutput += strlen(ptrOutput);
        if (modulusIsZero)
        {
          pstFactorInfo = factorInfoInteger;
        }
        else
        {
          pstFactorInfo = factorInfo;
        }
        // Output factors
        *ptrOutput++ = '<';
        *ptrOutput++ = 'u';
        *ptrOutput++ = 'l';
        *ptrOutput++ = '>';
        if (pretty == TEX)
        {
          strcpy(ptrOutput, "<li>\\begin{array}{l}</li>");
          ptrOutput += strlen(ptrOutput);
        }
        if (!modulusIsZero)
        {
          IntArray2BigInteger(&poly4[degree * nbrLimbs], &operand5);
        }
        if ((operand5.nbrLimbs != 1 || operand5.limbs[0].x != 1 || operand5.sign == SIGN_NEGATIVE) || nbrFactorsFound == 0)
        {     // Leading coefficient is not 1 or degree is zero.
          *ptrOutput++ = '<';
          *ptrOutput++ = 'l';
          *ptrOutput++ = 'i';
          *ptrOutput++ = '>';
          if (operand5.sign == SIGN_NEGATIVE)
          {
            strcpy(ptrOutput, " &minus;");
            ptrOutput += strlen(ptrOutput);
          }
          Bin2Dec(operand5.limbs, ptrOutput, operand5.nbrLimbs, groupLength);
          ptrOutput += strlen(ptrOutput);
          *ptrOutput++ = '<';
          *ptrOutput++ = '/';
          *ptrOutput++ = 'l';
          *ptrOutput++ = 'i';
          *ptrOutput++ = '>';
        }
        for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
        {
          *ptrOutput++ = '<';
          *ptrOutput++ = 'l';
          *ptrOutput++ = 'i';
          *ptrOutput++ = '>';
          if (pretty == TEX)
          {
            strcpy(ptrOutput, "\\bullet\\,\\,");
            ptrOutput += strlen(ptrOutput);
          }
          outputPolynomialFactor(ptrOutput, groupLength, pstFactorInfo);
          ptrOutput += strlen(ptrOutput);
          if (pretty == TEX)
          {
            strcpy(ptrOutput, "\\\\");
            ptrOutput += strlen(ptrOutput);
          }
          *ptrOutput++ = '<';
          *ptrOutput++ = '/';
          *ptrOutput++ = 'l';
          *ptrOutput++ = 'i';
          *ptrOutput++ = '>';
          pstFactorInfo++;
        }
        if (pretty == TEX)
        {
          strcpy(ptrOutput, "<li>\\end{array}</li>");
          ptrOutput += strlen(ptrOutput);
        }
        *ptrOutput++ = '<';
        *ptrOutput++ = '/';
        *ptrOutput++ = 'u';
        *ptrOutput++ = 'l';
        *ptrOutput++ = '>';
      }
      if (modulusIsZero)
      {
        strcpy(ptrOutput, lang ? "<h2>Raíces</h2>" : "<h2>Roots</h2>");
        ptrOutput += strlen(ptrOutput);
        if (degree > 1)
        {
          strcpy(ptrOutput, lang ? "Las " : "The ");
          ptrOutput += strlen(ptrOutput);
          int2dec(&ptrOutput, degree);
          strcpy(ptrOutput, lang ? " raíces son:</p>" : " roots are:</p>");
          ptrOutput += strlen(ptrOutput);
        }
        strcpy(ptrOutput, "<ul>");
        ptrOutput += strlen(ptrOutput);
        if (pretty == TEX)
        {
          strcpy(ptrOutput, "<li>\\begin{array}{l}</li>");
          ptrOutput += strlen(ptrOutput);
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
          strcpy(ptrOutput, "<li>\\end{array}</li>");
          ptrOutput += strlen(ptrOutput);
        }
      }
      strcpy(ptrOutput, "</ul>");
      ptrOutput += strlen(ptrOutput);
      // Show time only when factoring, not when just evaluating polynomial.
      showElapsedTime(&ptrOutput);
    }
  }
  strcpy(ptrOutput, lang ? "<p>" COPYRIGHT_SPANISH "</p>" :
                           "<p>" COPYRIGHT_ENGLISH "</p>");
}

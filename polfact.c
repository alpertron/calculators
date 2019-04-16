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

#ifdef __EMSCRIPTEN__
int attemptNbr;
char *ptrPercentageOutput;
#endif

extern int poly4[1000000];
extern int poly5[1000000];
extern int polyS[1000000];
int polyNonRepeatedFactors[1000000];
int polyToFactor[1000000];
int origPolyToFactor[1000000];
int tempPoly[1000000];
int polyBackup[1000000];
int polyLiftedRecord[1000000];
extern int polyLifted[1000000];
static int factorX[] = { 1, 0, 1, 1 }; // Polynomial is x.
static BigInteger bound, trailingCoeff, leadingCoeff, halfPowerMod;
static BigInteger contentPolyToFactor;
struct sFactorInfo factorInfoRecord[MAX_DEGREE];
struct sFactorInfo factorInfoInteger[MAX_DEGREE];
int polyInteger[1000000];
static int arrNbrFactors[MAX_DEGREE];
static int FactorModularPolynomial(int inputMontgomery)
;

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
        char *ptrOutput = output;
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
        databack(output);
      }
#endif
      // Copy polynomial to factor to poly3 and set leading coefficient to 1.
      ptrValue1 = &poly3[polyDegree*nbrLimbs];
      memcpy(poly3, ptrPolyToFactor, (ptrValue1 - &poly3[0])*sizeof(int));
      SetNumberToOne(ptrValue1);  // Set leading coefficient to 1.
      powerPolynomial(poly1, poly3, polyDegree, &primeMod, poly2, NULL);
      memcpy(poly1, poly2, polyDegree*nbrLimbs * sizeof(int));
      // Subtract x.
      UncompressBigInteger(&poly2[nbrLimbs], &operand1);
      memcpy(operand2.limbs, MontgomeryMultR1, NumberLength*sizeof(limb));
      operand2.nbrLimbs = NumberLengthR1;
      SubtBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
      CompressBigInteger(&poly2[nbrLimbs], &operand1);
      // Perform Gcd.
      degreeMin = getDegreePoly(poly2, polyDegree - 1);
      PolyModularGcd(poly3, polyDegree, poly2, degreeMin, polyMultTemp, &degreeGcd);
      if (degreeGcd == polyDegree)
      {
        pstFactorInfo->expectedDegree = currentDegree;
        polyDegree = 0;
      }
      else if (degreeGcd > 0)
      {         // Non-trivial factor of polynomial has been found.
                // Divide polynomial by GCD. Put the GCD in the first limbs
                // and the quotient in the last limbs.
        ptrValue1 = polyMultTemp + degreeGcd*nbrLimbs;
        SetNumberToOne(ptrValue1);
        DividePolynomial(poly3, polyDegree, polyMultTemp, degreeGcd, poly2);
        // Quotient located in poly2.
        pstNewFactorInfo = &factorInfo[nbrFactorsFound++];
        pstNewFactorInfo->ptr = ptrPolyToFactor;
        pstNewFactorInfo->degree = degreeGcd;
        pstNewFactorInfo->multiplicity = pstFactorInfo->multiplicity;
        pstNewFactorInfo->expectedDegree = currentDegree;
        pstFactorInfo->degree = polyDegree - degreeGcd;
        pstFactorInfo->ptr = &ptrPolyToFactor[degreeGcd*nbrLimbs];
        memcpy(ptrPolyToFactor, polyMultTemp, degreeGcd*nbrLimbs*sizeof(int));
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
    databack(output);
  }
#else
  (void)percentage;
#endif
}


// Perform Cantor-Zassenhaus algorithm to factor polynomials of the same degree.
static void SameDegreeFactorization(void)
{
  unsigned int seed = 1;  // Initialize pseudorandom sequence.
  struct sFactorInfo *pstFactorInfo = factorInfo;
  struct sFactorInfo *pstNewFactorInfo;
  int *ptrValue1, *ptrPolyToFactor;
  int nbrFactor, currentDegree, index, degreeGcd;
  int primeInt = (int)primeMod.limbs[0].x;
  int nbrLimbs = primeMod.nbrLimbs + 1;
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
#ifdef __EMSCRIPTEN__
    attemptNbr = 1;
#endif
    if (isCharacteristic2 == 0)
    { // If prime is not 2,
      // Calculate operand2 <- (prime^degree-1)/2
      // Use operand1 as temporary variable to store the exponent.
      operand1.limbs[0].x = pstFactorInfo->expectedDegree & MAX_VALUE_LIMB;
      operand1.nbrLimbs = 1;
      BigIntPower(&primeMod, &operand1, &operand4);
      BigIntDivide2(&operand4);
    }
    ptrPolyToFactor = pstFactorInfo->ptr;
    GetPolyInvParm(polyDegree, ptrPolyToFactor);
    for (;;)
    {
#ifdef __EMSCRIPTEN__
      char *ptrOutput = output;
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
      ptrValue1 = &poly3[pstFactorInfo->degree*nbrLimbs];
      memcpy(poly3, ptrPolyToFactor, (ptrValue1 - &poly3[0])*sizeof(int));
      SetNumberToOne(ptrValue1);  // Set leading coefficient to 1.
      ptrValue1 = poly1;
      if (nbrLimbs > 2)
      {    // Pseudorandom number range can be from 0 to MAX_LIMB_VALUE-1.
        for (currentDegree = 0; currentDegree < polyDegree; currentDegree++)
        {
          *ptrValue1 = 1;
          seed = (314159265U)*seed + 123456789;
          *(ptrValue1 + 1) = seed & MAX_VALUE_LIMB;
          ptrValue1 += nbrLimbs;
        }
      }
      else
      {   // Pseudorandom number range can be from 0 to prime-1.
        for (currentDegree = 0; currentDegree < polyDegree; currentDegree++)
        {
          *ptrValue1 = 1;
          seed = (314159265U)*seed + 123456789;
          *(ptrValue1 + 1) = ((seed >> 16)*primeInt) >> 16;
          ptrValue1 += nbrLimbs;
        }
      }
      if (isCharacteristic2 == 0)
      { // If prime is not 2: compute (random poly)^((p^d-1)/2)
        powerPolynomial(poly1, poly3, polyDegree, &operand4, poly2, percentageCallback);
        // Subtract 1.
        UncompressBigInteger(&poly2[0], &operand1);
        SubtBigNbrMod(operand1.limbs, MontgomeryMultR1, operand1.limbs);
        CompressBigInteger(&poly2[0], &operand1);
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
            UncompressBigInteger(&poly1[index*nbrLimbs], &operand1);
            UncompressBigInteger(&poly2[index*nbrLimbs], &operand2);
            AddBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
            CompressBigInteger(&poly2[index*nbrLimbs], &operand1);
          }
        }
      }
      PolyModularGcd(poly3, polyDegree, poly2, getDegreePoly(poly2, polyDegree - 1), polyMultTemp, &degreeGcd);
      if (degreeGcd != 0 && degreeGcd != polyDegree)
      {   // Non-trivial factor found.
        ptrValue1 = polyMultTemp + degreeGcd*nbrLimbs;
        SetNumberToOne(ptrValue1);
        DividePolynomial(poly3, polyDegree, polyMultTemp, degreeGcd, poly2);
        // Quotient located in poly2.
        pstNewFactorInfo = &factorInfo[nbrFactorsFound++];
        pstNewFactorInfo->ptr = &ptrPolyToFactor[degreeGcd*nbrLimbs];
        pstNewFactorInfo->degree = polyDegree - degreeGcd;
        pstNewFactorInfo->multiplicity = pstFactorInfo->multiplicity;
        pstNewFactorInfo->expectedDegree = pstFactorInfo->expectedDegree;
        pstFactorInfo->degree = degreeGcd;
        memcpy(ptrPolyToFactor, polyMultTemp, degreeGcd*nbrLimbs*sizeof(int));
        memcpy(pstNewFactorInfo->ptr, poly2, (polyDegree - degreeGcd)*nbrLimbs*sizeof(int));
        polyDegree = degreeGcd;
#ifdef __EMSCRIPTEN__
        attemptNbr = 0;
#endif
        if (pstFactorInfo->expectedDegree == pstFactorInfo->degree)
        {
          break;
        }
        GetPolyInvParm(polyDegree, ptrPolyToFactor);
      }
#ifdef __EMSCRIPTEN__
      attemptNbr++;
#endif
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
        ptrValue1 = pstFactorInfo->ptr + (pstFactorInfo->degree - 1)*nbrLimbs;
        ptrValue2 = pstFactorInfo2->ptr + (pstFactorInfo->degree - 1)*nbrLimbs;
        for (currentDegree = pstFactorInfo->degree - 1; currentDegree >= 0; currentDegree--)
        {
          if (*ptrValue1 == *ptrValue2)
          {
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

// Generate integer polynomial from modular polynomial.
static void GenerateIntegerPolynomial(int *polyMod, int *polyInt, int degreePoly)
{
  int currentDegree;
  *polyInt = degreePoly;
  int *ptrSrc = polyMod;
  int *ptrDest = polyInt+1;
  for (currentDegree = 0; currentDegree <= degreePoly; currentDegree++)
  {
    NumberLength = powerMod.nbrLimbs;
    UncompressBigInteger(ptrSrc, &operand1);
    ptrSrc += 1 + *ptrSrc;
    // If operand1 >= halfPowerMod, subtract powerMod.
    BigIntSubt(&operand1, &halfPowerMod, &operand2);
    if (operand2.sign == SIGN_POSITIVE)
    {
      BigIntSubt(&operand1, &powerMod, &operand1);
    }
    CompressBigInteger(ptrDest, &operand1);
    ptrDest += 1 + numLimbs(ptrDest);
  }
}

// Insert new factor in factorInfoInteger array sorting this array
// on ascending order of degree, and then on ascending order of
// leading coefficients. Merge equal factors.
static void InsertIntegerPolynomialFactor(int *ptrFactor, int degreePoly)
{
  struct sFactorInfo *pstFactorInfoInteger, *pstFactorInfo;
  int indexNewFactor[MAX_DEGREE];
  int indexOldFactor[MAX_DEGREE];
  int currentDegree, index;
  int *ptrIndex, *ptrOldFactor;

  // Fill indexes to start of each coefficient.
  ptrIndex = &indexNewFactor[0];
  index = 0;
  for (currentDegree = 0; currentDegree <= degreePoly; currentDegree++)
  {
    *ptrIndex++ = index;
    index += numLimbs(ptrFactor + index) + 1;
  }
  for (pstFactorInfoInteger = factorInfoInteger;
    pstFactorInfoInteger->ptrPolyLifted != NULL;
    pstFactorInfoInteger++)
  {
    if (pstFactorInfoInteger->degree < degreePoly)
    {
      continue;
    }
    if (pstFactorInfoInteger->degree > degreePoly)
    {
      break;
    }
    // Degree of the new factor is the same as the one already stored.
    ptrOldFactor = pstFactorInfoInteger->ptrPolyLifted;
    ptrIndex = &indexOldFactor[0];
    index = 0;
    for (currentDegree = 0; currentDegree <= degreePoly; currentDegree++)
    {
      *ptrIndex++ = index;
      index += numLimbs(ptrOldFactor + index) + 1;
    }
    // Compare coefficients.
    for (currentDegree = degreePoly; currentDegree >= 0; currentDegree--)
    {
      UncompressBigIntegerB(ptrFactor + indexNewFactor[currentDegree], &operand1);
      UncompressBigIntegerB(ptrOldFactor + indexOldFactor[currentDegree], &operand2);
      BigIntSubt(&operand1, &operand2, &operand1);
      if (!BigIntIsZero(&operand1))
      {
        break;
      }
    }
    if (currentDegree < 0)
    {   // Both polynomials are the same.
      pstFactorInfoInteger->multiplicity++;
      return;
    }
    if (operand1.sign == SIGN_NEGATIVE)
    {
      break;
    }
  }
  // Move elements of array factorInfoInteger to make room for new factor.
  for (pstFactorInfo = pstFactorInfoInteger;
    pstFactorInfo->ptrPolyLifted != NULL;
    pstFactorInfo++)
  {
  }
  if (pstFactorInfo - pstFactorInfoInteger > 0)
  {
    memmove(pstFactorInfoInteger + 1, pstFactorInfoInteger,
      (pstFactorInfo - pstFactorInfoInteger)*sizeof(*pstFactorInfo));
  }
  pstFactorInfoInteger->multiplicity = 1;
  pstFactorInfoInteger->ptrPolyLifted = ptrFactor;
  pstFactorInfoInteger->degree = degreePoly;
}

// Input: values = degree, coefficient degree 0, coefficient degree 1, etc.
// Output: factorInfo = structure that holds the factors.
static int FactorPolyOverIntegers(void)
{
  int degree = values[0];
  int degree1, degree2;
  int prime;
  int primeRecord = 0;
  int exponRecord = 0;
  int expon, maxDegreeFactor;
  int degreeGcdMod;
  int *ptrSrc, *ptrDest;
  int attemptNbr;
  int nbrFactorsRecord;
  int halfDegree;
  int currentFactor;
  int currentDegree, sumOfDegrees;
  int polXprocessed = FALSE;
  int nbrTmp[1000], nbrTmp2[1000], nbrTmp3[1000];
  int *ptrFactorIntegerBak;
  struct sFactorInfo *pstFactorInfoOrig, *pstFactorInfoRecord;
  struct sFactorInfo *pstFactorInfoInteger = factorInfoInteger;
  int *ptrFactorInteger = polyInteger;
  memset(factorInfoInteger, 0, sizeof(factorInfoInteger));
  getContent(values, &contentPolyToFactor);
  CopyPolynomial(&origPolyToFactor[1], &values[1], degree);
  origPolyToFactor[0] = degree;
  // polyToFactor -> original polynomial / content of polynomial.
  // Let n the degree of the least coefficient different from zero.
  // Then x^n divides the polynomial.
  polyToFactor[0] = degree;
  ptrSrc = &origPolyToFactor[1];
  ptrDest = &polyToFactor[1];
  for (currentDegree = 0; currentDegree <= degree; currentDegree++)
  {
    UncompressBigIntegerB(ptrSrc, &operand1);
    if (polXprocessed == FALSE)
    {
      if (BigIntIsZero(&operand1))
      {
        ptrSrc += 1 + numLimbs(ptrSrc);
        continue;
      }
      // x^currentDegree divides the original polynomial.
      if (currentDegree > 0)
      {
        pstFactorInfoInteger->degree = 1;
        pstFactorInfoInteger->expectedDegree = 1;
        pstFactorInfoInteger->multiplicity = currentDegree;
        pstFactorInfoInteger->ptrPolyLifted = factorX;
        pstFactorInfoInteger++;
        polyToFactor[0] = degree - currentDegree;
      }
      polXprocessed = TRUE;      
    }
    BigIntDivide(&operand1, &contentPolyToFactor, &operand2);
    NumberLength = operand2.nbrLimbs;
    CompressBigInteger(ptrDest, &operand2);
    ptrSrc += 1 + numLimbs(ptrSrc);
    ptrDest += 1 + numLimbs(ptrDest);
  }
  while (polyToFactor[0] != 0)
  {    // At least degree 1.
       // The trailing coefficient of factors must divide the product of the trailing and leading
       // coefficients of the original polynomial.
       // Get trailing coefficient.
    ptrSrc = &values[1];
    UncompressBigIntegerB(ptrSrc, &operand1);
       // Get leading coefficient.
    for (degree1 = 0; degree1 < degree; degree1++)
    {
      ptrSrc += 1 + numLimbs(ptrSrc);
    }
    UncompressBigIntegerB(ptrSrc, &operand2);
    BigIntMultiply(&operand1, &operand2, &trailingCoeff);
    // Compute F/gcd(F, F') where F is the polynomial to factor.
    // In the next loop we will attempt to factor gcd(F, F').
    degree = polyToFactor[0];
    CopyPolynomial(&polyNonRepeatedFactors[1], &polyToFactor[1], degree);
    CopyPolynomial(&tempPoly[1], &polyToFactor[1], degree);
    tempPoly[0] = degree;
    DerPolynomial(tempPoly);
    PolynomialGcd(tempPoly, polyToFactor, polyToFactor);
    polyNonRepeatedFactors[0] = degree;
    DivideIntegerPolynomial(polyNonRepeatedFactors, polyToFactor, TYPE_DIVISION);
    prime = 3;
        // Find Knuth-Cohen bound for coefficients of polynomial factors:
        // If polynomial B divides A we have for all j:
        // |Bj| <= binomial(n-1, j)*SUM(i, |Ai|^2))^(1/2) + binomial(n-1, j-1) * |Am|
        // where m is the degree of A and n is the degree of B.
        // Maximum degree to be considered is n = ceil(m/2).
        // We need to find max(Bj).
    modulusIsZero = 1;
    degree = polyNonRepeatedFactors[0];
    ptrSrc = &polyNonRepeatedFactors[1];
    UncompressBigIntegerB(ptrSrc, &operand1);
    if (degree < 0)
    {      // Monomial.
      maxDegreeFactor = (-degree + 1) / 2;
      operand1.sign = SIGN_POSITIVE;
      CopyBigInt(&operand3, &operand1);         // Get leading coefficient.
    }
    else
    {      // Polynomial.
      maxDegreeFactor = (degree + 1) / 2;
      BigIntMultiply(&operand1, &operand1, &operand1);
      for (degree1 = 1; degree1 <= degree; degree1++)
      {
        ptrSrc += 1 + numLimbs(ptrSrc);
        UncompressBigIntegerB(ptrSrc, &operand3);   // The last loop sets operand3 to the leading coefficient.
        BigIntMultiply(&operand3, &operand3, &operand2);
        BigIntAdd(&operand1, &operand2, &operand1);
      }
      squareRoot(operand1.limbs, operand2.limbs, operand1.nbrLimbs, &operand2.nbrLimbs);
      CopyBigInt(&operand1, &operand2);
    }
          // Loop that finds the maximum value of bound for |Bj|.
    intToBigInteger(&operand2, 1);  // binomial(n-1, 0)
    UncompressBigIntegerB(&polyNonRepeatedFactors[1], &bound);  // bound <- |A0|
    for (degree1 = 1; degree1 <= maxDegreeFactor; degree1++)
    {
      CopyBigInt(&operand4, &operand2);
      multint(&operand2, &operand2, maxDegreeFactor - degree1);
      subtractdivide(&operand2, 0, degree1);
      BigIntMultiply(&operand1, &operand2, &operand5);
      BigIntMultiply(&operand3, &operand4, &operand4);
      BigIntAdd(&operand5, &operand4, &operand5);
        // If operand5 > bound, set bound to operand5.
      BigIntSubt(&operand5, &bound, &operand4);
      if (operand4.sign == SIGN_POSITIVE)
      {
        CopyBigInt(&bound, &operand5);
      }
    }
        // If the modular factorization has k irreducible factors, the
        // number of combinations to try will be 2^k because all
        // factors are different. So up to 5
        // different prime moduli are tested to minimize the number
        // of factors. If the number of factors is less than 10,
        // no more modular factorizations are attempted.

        // The big number ensures that the first copy is done.
    nbrFactorsRecord = 100000;
    for (attemptNbr = 1; attemptNbr < 5; attemptNbr++)
    {
      int factorNbr, *ptrPolyLiftedRecord;
      int nbrFactors;
      do
      {   // Loop that finds a prime modulus such that the factorization
          // has no repeated factors. That means gcd(F, F') = 1 (mod prime).
          // This is required because Hensel lift does not work when
          // repeated factors are present.
        prime = nextPrime(prime);
        modulusIsZero = 0;
        intToBigInteger(&primeMod, prime);
        computePower(1);
        intToBigInteger(&operand5, 1);
        degree2 = getModPolynomial(&poly2[1], polyNonRepeatedFactors, &operand5);
        poly2[0] = degree2;
        DerPolynomial(poly2);   // This function overwrites poly1.
        degree1 = getModPolynomial(poly1, polyNonRepeatedFactors, &operand5);
        PolyModularGcd(poly1, degree1, &poly2[1], poly2[0], poly3, &degreeGcdMod);
      } while (degreeGcdMod > 0);
      // Find expon such that prime^expon >= 2 * bound
      BigIntMultiplyBy2(&bound);
      intToBigInteger(&operand1, prime);
      for (expon = 1; ; expon++)
      {
        // Compare operand1 = prime^expon against 2 * bound.
        BigIntSubt(&operand1, &bound, &operand4);
        if (operand4.sign == SIGN_POSITIVE)
        {
          break;     // prime^expon >= 2 * bound -> go out.
        }
        multint(&operand1, &operand1, prime);
      }
      modulusIsZero = 1;
      intToBigInteger(&primeMod, prime);
      computePower(expon);
      exponentMod = expon;
      intToBigInteger(&operand5, 1);
      degree = getModPolynomial(&poly1[1], polyNonRepeatedFactors, &operand5);
      poly1[0] = degree;
      polyBackup[0] = values[0];
      CopyPolynomial(&polyBackup[1], &values[1], values[0] >= 0 ? values[0] : 0);
      values[0] = degree;
      CopyPolynomial(&values[1], &poly1[1], degree);
      memset(factorInfo, 0, sizeof(factorInfo));
      modulusIsZero = 0;
      FactorModularPolynomial(FALSE);   // Input is not in Montgomery notation.
          // Get number of factors found.
      pstFactorInfoOrig = factorInfo;
      nbrFactors = 0;
      for (factorNbr = 0; factorNbr < MAX_DEGREE; factorNbr++)
      {
        if (pstFactorInfoOrig->ptr == NULL)
        {    // No more factors.
          break;
        }
        nbrFactors++;
        pstFactorInfoOrig++;
      }
      if (nbrFactors < nbrFactorsRecord)
      {    // Copy factors found to records arrays.
        primeRecord = prime;
        exponRecord = expon;
        pstFactorInfoOrig = factorInfo;
        pstFactorInfoRecord = factorInfoRecord;
        ptrPolyLiftedRecord = polyLiftedRecord;
        nbrFactorsRecord = 0;
        for (factorNbr = 0; factorNbr < MAX_DEGREE; factorNbr++)
        {
          if (pstFactorInfoOrig->ptrPolyLifted == NULL)
          {    // No more factors.
            break;
          }
          *pstFactorInfoRecord = *pstFactorInfoOrig;
          pstFactorInfoRecord->ptrPolyLifted = ptrPolyLiftedRecord;
          nbrFactorsRecord += pstFactorInfoOrig->multiplicity;
          ptrPolyLiftedRecord = CopyPolynomialFixedCoeffSize(ptrPolyLiftedRecord,
             pstFactorInfoOrig->ptrPolyLifted,
             pstFactorInfoOrig->degree - 1, powerMod.nbrLimbs + 1);
          *ptrPolyLiftedRecord++ = 1;    // Leading coefficient should be 1.
          *ptrPolyLiftedRecord++ = 1;
          pstFactorInfoOrig++;
          pstFactorInfoRecord++;
        }
        if (factorNbr < MAX_DEGREE)
        {
          pstFactorInfoRecord->ptrPolyLifted = NULL;
        }
        if (nbrFactors < 10)
        {
          break;    // Up to 512 factors: test all of them.
        }
      }
    }
    prime = primeRecord;
    expon = exponRecord;
    intToBigInteger(&primeMod, prime);
    computePower(expon);
    exponentMod = expon;
    // Combine factors. Ensure that the trailing coefficient multiplied by the leading coefficient
    // of the polynomial to factor divides trailingCoeff.
    // The sum of degree of factors must not exceed half the degree of the original polynomial.
    // Save in polyLifted the cumulative products, so we do not have to repeat multiplications.
    // The multiplicity of factors is always 1.
    CopyBigInt(&halfPowerMod, &powerMod);
    subtractdivide(&halfPowerMod, -1, 2); // halfPowerMod <- (powerMod+1)/2
    sumOfDegrees = 0;
    memset(arrNbrFactors, 0, sizeof(arrNbrFactors));
    halfDegree = polyNonRepeatedFactors[0] / 2;
    // Get leading coefficient of polyNonRepeatedFactors.
    ptrSrc = &polyNonRepeatedFactors[1];
    for (degree1 = 0; degree1 < degree; degree1++)
    {
      ptrSrc += 1 + numLimbs(ptrSrc);
    }
    UncompressBigIntegerB(ptrSrc, &leadingCoeff);
    ptrDest = polyLifted;
    for (currentFactor = 0; currentFactor <= nbrFactorsRecord; currentFactor++)
    {   // Initialize products to leading coefficient.
      int nbrLimbs = leadingCoeff.nbrLimbs;
      memcpy(ptrDest, leadingCoeff.limbs, nbrLimbs * sizeof(limb));
      memset(ptrDest + nbrLimbs, 0, (powerMod.nbrLimbs - nbrLimbs) * sizeof(limb));
      ptrDest += powerMod.nbrLimbs;
    }
    for (;;)
    {       // Get next product.
      int degreeProd, degreeFactor;
      int *ptrTrailingCoeff;

      pstFactorInfoRecord = factorInfoRecord;
      for (currentFactor = 0; currentFactor < nbrFactorsRecord; currentFactor++)
      {
        arrNbrFactors[currentFactor]++;
        sumOfDegrees += pstFactorInfoRecord->degree;
        if (arrNbrFactors[currentFactor] <= pstFactorInfoRecord->multiplicity &&
          sumOfDegrees <= halfDegree)
        {
          break;
        }
        sumOfDegrees -= pstFactorInfoRecord->degree * arrNbrFactors[currentFactor];
        arrNbrFactors[currentFactor] = 0;
        memcpy(&polyLifted[currentFactor * powerMod.nbrLimbs],
          &polyLifted[(currentFactor+1) * powerMod.nbrLimbs], powerMod.nbrLimbs*sizeof(int));
        pstFactorInfoRecord++;
      }
      if (currentFactor == nbrFactorsRecord)
      {            // All factors found.
        break;
      }
      ptrTrailingCoeff = &polyLifted[(currentFactor+1) * powerMod.nbrLimbs];
      UncompressIntLimbs(pstFactorInfoRecord->ptrPolyLifted, (limb *)nbrTmp2, powerMod.nbrLimbs);
      MultBigNbrModN(ptrTrailingCoeff, nbrTmp2, nbrTmp3, (int *)powerMod.limbs, powerMod.nbrLimbs);
      do
      {
        ptrTrailingCoeff = &polyLifted[currentFactor * powerMod.nbrLimbs];
        memcpy(ptrTrailingCoeff, nbrTmp3, powerMod.nbrLimbs * sizeof(int));
      } while (--currentFactor >= 0);
      NumberLength = powerMod.nbrLimbs;
      UncompressLimbsBigInteger((limb *)ptrTrailingCoeff, &operand1);
      operand1.sign = SIGN_POSITIVE;
      // If ptrTrailingCoeff >= halfPowerMod, subtract powerMod.
      BigIntSubt(&operand1, &halfPowerMod, &operand2);
      if (operand2.sign == SIGN_POSITIVE)
      {
        BigIntSubt(&operand1, &powerMod, &operand1);
      }
      BigIntRemainder(&trailingCoeff, &operand1, &operand2);
      if (!BigIntIsZero(&operand2))
      {    // Factor not found, Try next one.
        continue;
      }
      // Test whether the factor divides the original polynomial.
      // Initialize factor to 1.
      // Coefficients must be converted to Montgomery notation.
      modulusIsZero = 0;  // Perform modular operations.
      degreeProd = 0;
      poly1[0] = NumberLength;
      memcpy(&poly1[1], MontgomeryMultR1, NumberLength * sizeof(int));
      pstFactorInfoRecord = factorInfoRecord;
      for (currentFactor = 0; currentFactor < nbrFactorsRecord; currentFactor++)
      {
        if (arrNbrFactors[currentFactor] > 0)
        {
          int *ptrValue1 = pstFactorInfoRecord->ptrPolyLifted;    // Source
          int *ptrValue2 = poly2;                                 // Destination
          int nbrLength;
          degreeFactor = pstFactorInfoRecord->degree;
          for (currentDegree = 0; currentDegree <= degreeFactor; currentDegree++)
          {
            nbrLength = 1 + numLimbs(ptrValue1);
            memcpy(ptrValue2, ptrValue1, nbrLength * sizeof(int));
            ptrValue1 += 1 + NumberLength;
            ptrValue2 += 1 + NumberLength;
          }
          // Convert factor to Montgomery notation.
          polyToMontgomeryNotation(poly2, degreeFactor + 1);
          MultPolynomial(degreeProd, degreeFactor, poly1, poly2);
          degreeProd += degreeFactor;
          ptrValue1 = polyMultTemp;              // Source is the product
          ptrValue2 = poly1;                     // Destination
          degreeFactor = pstFactorInfoRecord->degree;
          for (currentDegree = 0; currentDegree <= degreeProd; currentDegree++)
          {
            nbrLength = 1 + numLimbs(ptrValue1);
            memcpy(ptrValue2, ptrValue1, nbrLength * sizeof(int));
            ptrValue1 += 1 + NumberLength;
            ptrValue2 += 1 + NumberLength;
          }
        }
        pstFactorInfoRecord++;
      }
      // Convert from Montgomery to standard notation.
      polyToStandardNotation(poly1, degreeProd + 1);
      // Multiply all coefficients by leadingCoeff and store in poly2.
      ptrSrc = poly1;
      ptrDest = poly2;
      CompressLimbsBigInteger((limb *)nbrTmp, &leadingCoeff);
      for (currentDegree = 0; currentDegree <= degreeProd; currentDegree++)
      {
        UncompressIntLimbs(ptrSrc, (limb *)nbrTmp2, powerMod.nbrLimbs);
        MultBigNbrModN(nbrTmp, nbrTmp2, nbrTmp3, (int *)powerMod.limbs, powerMod.nbrLimbs);
        CompressIntLimbs(ptrDest, (limb *)nbrTmp3, powerMod.nbrLimbs + 1);
        ptrSrc += 1 + NumberLength;
        ptrDest += 1 + numLimbs(ptrDest);
      }
      GenerateIntegerPolynomial(poly2, poly5, degreeProd);
      modulusIsZero = 1;   // Perform integer division.
      // Multiply all coefficients by leadingCoeff and store in polyS.
      polyS[0] = polyNonRepeatedFactors[0];
      ptrSrc = &polyNonRepeatedFactors[1];
      ptrDest = &polyS[1];
      for (currentDegree = 0; currentDegree <= polyS[0]; currentDegree++)
      {
        UncompressBigIntegerB(ptrSrc, &operand1);
        BigIntMultiply(&operand1, &leadingCoeff, &operand2);
        CompressBigInteger(ptrDest, &operand2);
        ptrSrc += 1 + numLimbs(ptrSrc);
        ptrDest += 1 + numLimbs(ptrDest);
      }
      DivideIntegerPolynomial(polyS, poly5, TYPE_MODULUS);
      if (polyS[0] != 0)
      {
        continue;    // Number to factor does not divide this polynomial.
      }
      // Get principal part of poly5 and store it to poly2.
      getContent(poly5, &operand4);   // Content of polynomial.
      poly2[0] = poly5[0];
      ptrSrc = &poly5[1];
      ptrDest = &poly2[1];
      for (currentDegree = 0; currentDegree <= poly5[0]; currentDegree++)
      {
        UncompressBigIntegerB(ptrSrc, &operand2);
        BigIntDivide(&operand2, &operand4, &operand3);
        CompressBigInteger(ptrDest, &operand3);
        ptrSrc += 1 + numLimbs(ptrSrc);
        ptrDest += 1 + numLimbs(ptrDest);
      }
      // Copy this principal part to poly5.
      CopyPolynomial(&poly5[1], &poly2[1], poly5[0]);
      DivideIntegerPolynomial(polyNonRepeatedFactors, poly5, TYPE_DIVISION);
      int degreePoly = poly5[0];
      ptrFactorIntegerBak = CopyPolynomial(ptrFactorInteger, &poly5[1], degreePoly);
      InsertIntegerPolynomialFactor(ptrFactorInteger, degreePoly);
      ptrFactorInteger = ptrFactorIntegerBak;
      pstFactorInfoInteger++;
      modulusIsZero = 0;   // Perform modular operations.
      // Discard factors already used.
      pstFactorInfoRecord = factorInfoRecord;
      for (currentFactor = 0; currentFactor < nbrFactorsRecord; currentFactor++)
      {
        if (arrNbrFactors[currentFactor] > 0)
        {
          pstFactorInfoRecord->multiplicity = 0;
        }
        pstFactorInfoRecord++;
      }
      // Restart finding factors.
      sumOfDegrees -= degreePoly;
      memset(arrNbrFactors, 0, sizeof(arrNbrFactors));
      polyLifted[0] = 1;                    // Initialize first product to 1.
      halfDegree = polyNonRepeatedFactors[0] / 2;
      if (powerMod.nbrLimbs > 1)
      {
        memset(&polyLifted[1], 0, powerMod.nbrLimbs - 1);
      }
    }
    // Polynomial is irreducible.
    if (polyNonRepeatedFactors[0] > 0)
    {    // Degree is greater than zero. Copy it to integer polynomial factor array.
      ptrFactorIntegerBak = CopyPolynomial(ptrFactorInteger, &polyNonRepeatedFactors[1], polyNonRepeatedFactors[0]);
      InsertIntegerPolynomialFactor(ptrFactorInteger, polyNonRepeatedFactors[0]);
      ptrFactorInteger = ptrFactorIntegerBak;
      pstFactorInfoInteger++;
    }
  }
  modulusIsZero = 1;
  CopyPolynomial(&values[1], &origPolyToFactor[1], origPolyToFactor[0]);
  values[0] = origPolyToFactor[0];
  CopyBigInt(&operand5, &contentPolyToFactor);
  // Find number of factors.
  for (nbrFactorsFound = 0; nbrFactorsFound < MAX_DEGREE; nbrFactorsFound++)
  {
    if (factorInfoInteger[nbrFactorsFound].ptrPolyLifted == NULL)
    {
      break;
    }
  }
  return EXPR_OK;
}

// Input: values = degree, coefficient degree 0, coefficient degree 1, etc.
// Output: factorInfo = structure that holds the factors.
static int FactorModularPolynomial(int inputMontgomery)
{
  struct sFactorInfo *pstFactorInfo;
  int currentDegree, nbrFactor, rc;
  int *ptrValue1;
  int nbrLimbsPrime = primeMod.nbrLimbs + 1; // Add 1 for length;
  degree = values[0];
  ptrValue1 = &values[1];
  for (currentDegree = 0; currentDegree <= degree; currentDegree++)
  {
    NumberLength = numLimbs(ptrValue1);
    UncompressBigInteger(ptrValue1, &operand1);
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
    CompressBigInteger(&valuesPrime[currentDegree*nbrLimbsPrime], &operand1);
    ptrValue1 += 1 + numLimbs(ptrValue1);
  }
  if (operand1.nbrLimbs == 1 && operand1.limbs[0].x == 0)
  {
    return EXPR_LEADING_COFF_MULTIPLE_OF_PRIME;
  }
  memcpy(&TestNbr, &primeMod, primeMod.nbrLimbs * sizeof(limb));
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
    SameDegreeFactorization();
  }
  if (inputMontgomery)
  {
    // Convert original polynomial from Montgomery to standard notation.
    OrigPolyFromMontgomeryToStandard();
  }
  rc = HenselLifting();
  if (rc != EXPR_OK)
  {
    return rc;
  }
  // Convert factors to standard notation.
  pstFactorInfo = factorInfo;
  for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    polyToStandardNotation(pstFactorInfo->ptrPolyLifted, pstFactorInfo->degree);
    pstFactorInfo++;
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
    outputPolynomial(ptrOutput, groupLength);
    ptrOutput += strlen(ptrOutput);
    if (onlyEvaluate == 0)
    {        // Show time only when factoring, not when just evaluating polynomial.
      showElapsedTime(&ptrOutput);
    }
  }
  strcpy(ptrOutput, lang ? "<p>" COPYRIGHT_SPANISH "</p>" :
                           "<p>" COPYRIGHT_ENGLISH "</p>");
}

#ifdef __EMSCRIPTEN__
EXTERNALIZE void doWork(void)
{
  int flags;
  int groupLen = 0;
  char *ptrData = inputString;
  while (*ptrData != ',')
  {
    groupLen = groupLen * 10 + (*ptrData++ - '0');
  }
  ptrData++;             // Skip comma.
  flags = *ptrData;
  lang = flags & 1;
  onlyEvaluate = (unsigned char)(flags & 2);
  superscripts = (unsigned char)(flags & 4);
  ptrData += 2;          // Skip flags and comma.
  polyFactText(ptrData, ptrData + strlen(ptrData) + 1, groupLen);
  ptrData += strlen(ptrData);
  databack(output);
}
#endif

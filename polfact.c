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
#ifdef __EMSCRIPTEN__
void databack(char *data);
int stamp(void);
extern int newStamp, oldStamp;
#endif

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
    // For each loop, raise this polynomial to the primeth power and 
    // then compute the gcd between the polynomial to be factored and
    // the computed polynomial less x. If the degree of GCD is > 0, then the
    // GCD is the product of all factors of degree indicated by currentDegree.
    for (currentDegree = 1; currentDegree * 2 <= polyDegree; currentDegree++)
    {
#ifdef __EMSCRIPTEN__
      newStamp = stamp();
      if (newStamp != oldStamp)
      {
        char *ptrOutput = output;
        oldStamp = newStamp;
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
        strcpy(ptrOutput, ".</p>");
        databack(output);
      }
#endif
      // Copy polynomial to factor to poly3 and set leading coefficient to 1.
      ptrValue1 = &poly3[polyDegree*nbrLimbs];
      memcpy(poly3, ptrPolyToFactor, (ptrValue1 - &poly3[0])*sizeof(int));
      SetNumberToOne(ptrValue1);  // Set leading coefficient to 1.
      powerPolynomial(poly1, poly3, polyDegree, &primeMod, poly2);
      memcpy(poly1, poly2, polyDegree*nbrLimbs*sizeof(int));
      // Subtract x.
      UncompressBigInteger(&poly2[nbrLimbs], &operand1);
      memcpy(operand2.limbs, MontgomeryMultR1, NumberLength*sizeof(limb));
      operand2.nbrLimbs = NumberLengthR1;
      SubtBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
      CompressBigInteger(&poly2[nbrLimbs], &operand1);
      // Perform Gcd.
      degreeMin = getDegreePoly(poly2, polyDegree - 1);
      PolynomialGcd(poly3, polyDegree, poly2, degreeMin, polyMultTemp, &degreeGcd);
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
      }
    }
    if (polyDegree > 0)
    {
      pstFactorInfo->expectedDegree = polyDegree;
    }
  }
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
    if (polyDegree < 2 || polyDegree == pstFactorInfo->expectedDegree)
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
      BigIntPower(&primeMod, &operand1, &operand4);
      BigIntDivide2(&operand4);
    }
    ptrPolyToFactor = pstFactorInfo->ptr;
    for (;;)
    {
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
        powerPolynomial(poly1, poly3, polyDegree, &operand4, poly2);
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
          multPolynomial(poly1, poly1, poly1, polyDegree, poly3);
          for (index = 0; index < polyDegree; index++)
          {
            UncompressBigInteger(&poly1[index*nbrLimbs], &operand1);
            UncompressBigInteger(&poly2[index*nbrLimbs], &operand2);
            AddBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
            CompressBigInteger(&poly2[index*nbrLimbs], &operand1);
          }
        }
      }
      PolynomialGcd(poly3, polyDegree, poly2, getDegreePoly(poly2, polyDegree - 1), polyMultTemp, &degreeGcd);
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
        if (pstFactorInfo->expectedDegree == pstFactorInfo->degree)
        {
          break;
        }
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

static int FactorPolynomial(char *input, int expo)
{
  struct sFactorInfo *pstFactorInfo;
  int currentDegree, nbrFactor;
  int *ptrValue1;
  int nbrLimbsPrime = primeMod.nbrLimbs + 1; // Add 1 for length;
  int rc = ComputePolynomial(input, expo);
  if (rc != 0)
  {
    return rc;
  }
  if (onlyEvaluate)
  {
    OrigPolyFromMontgomeryToStandard();
    return EXPR_OK;
  }
  // Generate polynomial mod prime.
  degree = values[0];
  ptrValue1 = &values[1];
  for (currentDegree = 0; currentDegree <= degree; currentDegree++)
  {
    UncompressBigInteger(ptrValue1, &operand1);
    // Convert from Montgomery to standard notation.
    operand2.limbs[0].x = 1;
    if (NumberLength > 1)
    {
      memset(&operand2.limbs[1], 0, (NumberLength - 1)*sizeof(limb));
    }
    modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    rc = BigIntRemainder(&operand1, &primeMod, &operand1);
    if (rc != EXPR_OK)
    {
      return rc;
    }
    CompressBigInteger(&valuesPrime[currentDegree*nbrLimbsPrime], &operand1);
    ptrValue1 += 1 + *ptrValue1;
  }
  if (operand1.nbrLimbs == 1 && operand1.limbs[0].x == 0)
  {
    return EXPR_LEADING_COFF_MULTIPLE_OF_PRIME;
  }
  memcpy(&TestNbr, &primeMod, primeMod.nbrLimbs*sizeof(limb));
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
  // Convert original polynomial from Montgomery to standard notation.
  OrigPolyFromMontgomeryToStandard();
  rc = HenselLifting();
  if (rc != EXPR_OK)
  {
    return rc;
  }
  // Convert factors to standard notation.
  pstFactorInfo = factorInfo;
  for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    ToStandardNotation(pstFactorInfo->ptrPolyLifted, pstFactorInfo->degree);
    pstFactorInfo++;
  }
  SortFactors(rc != 0 ? &primeMod : &powerMod);
  return rc;
}

void polyFactText(char *modText, char *polyText, int groupLen)
{
  enum eExprErr rc;
  int expon = 0;
  rc = ComputeExpression(modText, 1, &powerMod);
  if (output == NULL)
  {
    output = (char *)malloc(1000000);
  }
  if (output == NULL)
  {
    return;   // Go out if cannot generate output string.
  }
  if (rc == EXPR_OK)
  {
    if (powerMod.sign == SIGN_NEGATIVE || (powerMod.nbrLimbs == 1 && powerMod.limbs[0].x < 2))
    {
      rc = EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE;
    }
  }
  if (rc == EXPR_OK)
  {
    expon = PowerCheck(&powerMod, &primeMod);
    if (!isPseudoprime(&primeMod))
    {
      rc = EXPR_MODULUS_MUST_BE_PRIME_EXP;
    }
  }
  if (rc == EXPR_OK)
  {
    rc = FactorPolynomial(polyText, expon);
  }
  output[0] = '2';
  if (rc != EXPR_OK)
  {
    textErrorPol(output + 1, rc);
  }
  else
  {
    outputPolynomial(output + 1, groupLen);
  }
  strcat(output + 1, lang ? "<p>" COPYRIGHT_SPANISH "</p>" :
                            "<p>" COPYRIGHT_ENGLISH "</p>");
}

#ifdef __EMSCRIPTEN__
void doWork(char* data, int size)
{
  int flags;
  int groupLen = 0;
  char *ptrData = data;
  if (output == NULL)
  {
    output = malloc(3000000);
  }
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
  databack(output);
}
#endif

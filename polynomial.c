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

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include "bignbr.h"
#include "highlevel.h"
#include "polynomial.h"
#include "showtime.h"
#define STACK_OPER_SIZE            100
#define TOKEN_NUMBER               '0'
#define TOKEN_START_EXPON          '1'
#define TOKEN_END_EXPON            '2'
#define TOKEN_UNARY_MINUS          '3'
#define TOKEN_GCD                  '4'
#define TOKEN_DER                  '5'
extern char* ptrOutput2;
static void showPolynomial(char** pptrOutput, const int* ptrPoly, int polyDegree, int groupLength);
BigInteger primeMod;              // p
int exponentMod;                  // k
BigInteger powerMod;              // p^k
bool modulusIsZero;
int degree;
BigInteger operand1;
BigInteger operand2;
BigInteger operand3;
BigInteger operand4;
BigInteger operand5;
extern limb TestNbr[MAX_LEN];
extern limb MontgomeryMultR2[MAX_LEN];
extern limb MontgomeryMultR1[MAX_LEN];
extern int NumberLength;
extern int NumberLengthR1;
static int prime;
static int primeIndex;
int values[COMPRESSED_POLY_MAX_LENGTH];
int valuesPrime[COMPRESSED_POLY_MAX_LENGTH];
int poly1[COMPRESSED_POLY_MAX_LENGTH];
int poly2[COMPRESSED_POLY_MAX_LENGTH];
int poly3[COMPRESSED_POLY_MAX_LENGTH];
int poly4[COMPRESSED_POLY_MAX_LENGTH];
int poly5[COMPRESSED_POLY_MAX_LENGTH];
int polyS[COMPRESSED_POLY_MAX_LENGTH];
static int polyT[COMPRESSED_POLY_MAX_LENGTH];
int polyMultTemp[COMPRESSED_POLY_MAX_LENGTH];
int polyLifted[COMPRESSED_POLY_MAX_LENGTH];
int polyLiftedNew[COMPRESSED_POLY_MAX_LENGTH];
enum eOutput pretty;
bool onlyEvaluate = 0;
struct sFactorInfo factorInfo[MAX_DEGREE];
int nbrFactorsFound;
int *ptrOrigPoly;
int degreeOrigPoly;
int valuesIndex;
static int *ptrA[MAX_DEGREE];
static int* ptrNewFactors[MAX_DEGREE];

int numLimbs(const int *pLen)
{
  int nbrLimbs = *pLen;
  if (nbrLimbs < 0)
  {
    nbrLimbs = -nbrLimbs;
  }
  return nbrLimbs;
}

// Divide the polynomial by its leading coefficient.
// On output, operand1 is the inverse of the leading coefficient.
void ConvertToMonic(int *poly, int polyDegree)
{
  int nbrLimbs = NumberLength + 1;
  int currentDegree;
  if (NumberLength == 1)
  {         // Modulus size is one limb.
    int inverse = modInv(*(poly + (polyDegree * 2) + 1), TestNbr[0].x);
    int* ptrPoly = poly + 1;
    intToBigInteger(&operand1, inverse);
    if (TestNbr[0].x <= 32768)
    {
      for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
      {
        *ptrPoly = *ptrPoly * inverse % TestNbr[0].x;
        ptrPoly += 2;
      }
    }
    else
    {
      for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
      {
#ifdef _USING64BITS_
        *ptrPoly = (int)((uint64_t)*ptrPoly * inverse % TestNbr[0].x);
#else
        smallmodmult(*ptrPoly, inverse, (limb*)ptrPoly, TestNbr[0].x);
#endif
        ptrPoly += 2;
      }
    }
  }
  else
  {         // General case.
    IntArray2BigInteger(poly + (polyDegree * nbrLimbs), &operand1);
    ModInvBigNbr(operand1.limbs, operand1.limbs, TestNbr, NumberLength);
    for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
    {
      IntArray2BigInteger(poly + (currentDegree * nbrLimbs), &operand2);
      modmult(operand1.limbs, operand2.limbs, operand2.limbs);
      BigInteger2IntArray(poly + (currentDegree * nbrLimbs), &operand2);
    }
  }
}

// The content is the GCD of all coefficients with the sign equal to the sign of the leading coefficient.
int *getContent(int *poly, BigInteger *content)
{
  int* ptrPoly = poly;
  int currentDegree = *ptrPoly;
  ptrPoly++;
  if (currentDegree < 0)
  {
    currentDegree = 0;    // Monomial: process only one coefficient.
  }
  UncompressBigIntegerB(ptrPoly, content);
  for (; currentDegree > 0; currentDegree--)
  {
    ptrPoly += 1 + numLimbs(ptrPoly);
    UncompressBigIntegerB(ptrPoly, &operand1);
    BigIntGcd(content, &operand1, &operand2);
    CopyBigInt(content, &operand2);
    content->sign = operand1.sign;
  }
  return ptrPoly;
}

// Get polynomial divided content mod prime.
int getModPolynomial(int *polyMod, const int *poly, const BigInteger *content)
{
  int* ptrPolyMod = polyMod;
  const int* ptrPoly = poly;
  int currentDegree;
  int degreePoly = *ptrPoly;
  ptrPoly++;
  if (degreePoly < 0)
  {     // Monomial.
    for (currentDegree = degreePoly; currentDegree < 0; currentDegree++)
    {
      *ptrPolyMod = 1;       // Initialize coefficient to zero.
      ptrPolyMod++;
      *ptrPolyMod = 0;
      ptrPolyMod++;
    }
    NumberLength = numLimbs(ptrPoly);
    UncompressBigIntegerB(ptrPoly, &operand1);  // Get coefficient.
    (void)BigIntDivide(&operand1, content, &operand2);
    (void)BigIntRemainder(&operand2, &powerMod, &operand1);
    if (operand1.sign == SIGN_NEGATIVE)
    {
      BigIntAdd(&operand1, &powerMod, &operand1);
    }
    NumberLength = operand1.nbrLimbs;
    BigInteger2IntArray(ptrPolyMod, &operand1);
    return -degreePoly;
  }
  for (currentDegree = 0; currentDegree <= degreePoly; currentDegree++)
  {
    NumberLength = numLimbs(ptrPoly);
    UncompressBigIntegerB(ptrPoly, &operand1);  // Get coefficient.
    (void)BigIntDivide(&operand1, content, &operand2);
    (void)BigIntRemainder(&operand2, &powerMod, &operand1);
    if (operand1.sign == SIGN_NEGATIVE)
    {
      BigIntAdd(&operand1, &powerMod, &operand1);
    }
    NumberLength = operand1.nbrLimbs;
    BigInteger2IntArray(ptrPolyMod, &operand1);
    ptrPoly += 1 + numLimbs(ptrPoly);
    ptrPolyMod += 1 + numLimbs(ptrPolyMod);
  }
  return degreePoly;
}

void PolynomialGcd(int *argF, int *argG, int *gcd)
{
  static BigInteger contentF;
  static BigInteger contentG;
  static BigInteger contentH;
  static BigInteger modulus;
  static BigInteger gcdLeadingCoeff;
  const int *ptrLeadingCoeff;
  const int *ptrSrc;
  int *ptrDest;
  const int *ptrPrev;
  int potentialDegreeGcd;
  int tmp;
  int rc;
  if ((*argF == 0) && (*(argF + 1) == 1) && (*(argF + 2) == 0))
  {         // F equals zero, so the GCD is G.
    *gcd = *argG;
    (void)CopyPolynomial(gcd+1, argG+1, *argG);
    return;
  }
  if ((*argG == 0) && (*(argG + 1) == 1) && (*(argG + 2) == 0))
  {         // G equals zero, so the GCD is F.
    *gcd = *argF;
    (void)CopyPolynomial(gcd + 1, argF + 1, *argF);
    return;
  }
  initializeSmallPrimes(smallPrimes);
  prime = 65537;
  primeIndex = 6543;
  // Get content (gcd of all coefficients).
  ptrLeadingCoeff = getContent(argF, &contentF);
  UncompressBigIntegerB(ptrLeadingCoeff, &operand3);
  ptrLeadingCoeff = getContent(argG, &contentG);
  UncompressBigIntegerB(ptrLeadingCoeff, &operand2);
  BigIntGcd(&operand2, &operand3, &gcdLeadingCoeff);
  // Get potential degree of GCD as the minimum degree of arguments.
  potentialDegreeGcd = *argF;
  if (potentialDegreeGcd < 0)
  {           // Monomial
    potentialDegreeGcd = -potentialDegreeGcd;
  }
  tmp = *argG;
  if (tmp < 0)
  {           // Monomial
    tmp = -tmp;
  }
  if (potentialDegreeGcd > tmp)
  {
    potentialDegreeGcd = tmp;
  }
  intToBigInteger(&modulus, 1);  // Initialize modulus to 1.
  polyS[0] = 0;                  // Initialize polynomial g_m to zero.
  polyS[1] = 1;
  polyS[2] = 0;
  for (;;)
  {
    int degree1;
    int degree2;
    int degreeGcdMod;
    do
    { // Find next prime.
      primeIndex++;
      prime = smallPrimes[primeIndex];
      // Ensure that the prime does not divide the modulus or
      // the gcd of the leading coefficients.
    } while ((getRemainder(&modulus, prime) == 0) ||
      (getRemainder(&gcdLeadingCoeff, prime) == 0));
    modulusIsZero = false;
    intToBigInteger(&primeMod, prime);
    computePower(1);
    degree1 = getModPolynomial(poly1, argF, &contentF);
    degree2 = getModPolynomial(poly2, argG, &contentG);
    PolyModularGcd(poly1, degree1, poly2, degree2, poly3, &degreeGcdMod);
    modulusIsZero = true;
    if (degreeGcdMod == 0)
    {                                // Polynomial GCD is the GCD of both contents.
      BigIntGcd(&contentF, &contentG, &gcdLeadingCoeff);
      NumberLength = gcdLeadingCoeff.nbrLimbs;
      *gcd = 0;                      // Degree of GCD.
      BigInteger2IntArray(gcd + 1, &gcdLeadingCoeff);
      return;
    }
    if (degreeGcdMod > potentialDegreeGcd)
    {                                // Incorrect degree for GCD.
      continue;
    }
    if (degreeGcdMod < potentialDegreeGcd)
    {                                // Potential degree is incorrect.
      intToBigInteger(&modulus, 1);  // Initialize modulus to 1.
      polyS[0] = 0;                  // Initialize polynomial g_m to zero.
      polyS[1] = 1;
      polyS[2] = 0;
      potentialDegreeGcd = degreeGcdMod;
    }
    ptrSrc = &poly3[0];
    if ((polyS[0] == 0) && (polyS[1] == 1) && (polyS[2] == 0))
    {                                // If previous polynomial is zero, g_m <- (b * g_p) mod p.
      ptrDest = &polyS[1];
      for (degree = 0; degree <= potentialDegreeGcd; degree++)
      {
        UncompressBigIntegerB(ptrSrc, &operand1);
        BigIntMultiply(&operand1, &gcdLeadingCoeff, &operand1);
        intToBigInteger(&operand1, getRemainder(&operand1, prime));
        NumberLength = operand1.nbrLimbs;
        BigInteger2IntArray(ptrDest, &operand1);
        ptrSrc += 1 + numLimbs(ptrSrc);
        ptrDest += 1 + numLimbs(ptrDest);
      }
    }
    else
    {                                // Perform Chinese Remainder Theorem between g_p and g_m.
      int gcdLCModPrime;
      // coeff = c = a (mod p), c = b (mod q) => c = jq + b
      // a = jq + b (mod p) => j = (a-b)/q = r (mod p)
      // j = pk + r => c = (pk + r)q + b = pqk + rq + b => c = rq + b (mod pq).
      polyT[0] = potentialDegreeGcd;
      ptrPrev = &polyS[1];
      ptrDest = &polyT[1];
      gcdLCModPrime = getRemainder(&gcdLeadingCoeff, prime);
      for (degree = 0; degree <= potentialDegreeGcd; degree++)
      {
        int valueR;
        int invModPrime;
        int coeffModPrime = *(ptrSrc + 1);
        coeffModPrime = (int)(((int64_t)coeffModPrime*gcdLCModPrime)%prime);
        UncompressBigIntegerB(ptrPrev, &operand1);
        valueR = coeffModPrime - getRemainder(&operand1, prime);
        if (valueR < 0)
        {
          valueR += prime;
        }
        invModPrime = modInv(getRemainder(&modulus, prime), prime);
        valueR = (int)(((int64_t)valueR * (int64_t)invModPrime) % prime);
        multint(&operand2, &modulus, valueR);
        BigIntAdd(&operand2, &operand1, &operand2);
        NumberLength = operand2.nbrLimbs;
        BigInteger2IntArray(ptrDest, &operand2);
        ptrSrc += 2;
        ptrPrev += 1 + numLimbs(ptrPrev);
        ptrDest += 1 + numLimbs(ptrDest);
      }
      (void)CopyPolynomial(&polyS[1], &polyT[1], potentialDegreeGcd);
    }
    polyS[0] = potentialDegreeGcd;
    multint(&modulus, &modulus, prime);   // modulus <- modulus*prime.
    CopyBigInt(&operand4, &modulus);
    BigIntDivideBy2(&operand4);           // operand4 <- modulus / 2.
    // Get the primitive part of g_m in polyT.
    polyT[0] = potentialDegreeGcd;
    ptrSrc = &polyS[1];
    ptrDest = &polyT[1];
    for (degree = 0; degree <= potentialDegreeGcd; degree++)
    {
      UncompressBigIntegerB(ptrSrc, &operand1);
      // If coefficient is greater than modulus/2, 
      // let coefficient <- coefficient - modulus.
      BigIntSubt(&operand1, &operand4, &operand2);
      if (operand2.sign == SIGN_POSITIVE)
      {
        BigIntSubt(&operand1, &modulus, &operand1);
      }
      NumberLength = operand1.nbrLimbs;
      BigInteger2IntArray(ptrDest, &operand1);
      ptrSrc += 1 + numLimbs(ptrSrc);
      ptrDest += 1 + numLimbs(ptrDest);
    }
    (void)CopyPolynomial(&polyS[1], &polyT[1], potentialDegreeGcd);
    (void)getContent(polyS, &contentH);
    ptrSrc = &polyS[1];
    ptrDest = &polyT[1];
    for (degree = 0; degree <= potentialDegreeGcd; degree++)
    {
      UncompressBigIntegerB(ptrSrc, &operand1);
      (void)BigIntDivide(&operand1, &contentH, &operand2);
      NumberLength = operand2.nbrLimbs;
      BigInteger2IntArray(ptrDest, &operand2);
      ptrSrc += 1 + numLimbs(ptrSrc);
      ptrDest += 1 + numLimbs(ptrDest);
    }
    // if h divides f and h divides g, then return h.
    if (argF[0] < 0)
    {
      degree = 0;
    }
    else
    {
      degree = argF[0];
    }
    poly5[0] = argF[0];
    (void)CopyPolynomial(&poly5[1], &argF[1], degree);
    rc = DivideIntegerPolynomial(poly5, polyT, TYPE_MODULUS);
    if ((rc == EXPR_OK) && (poly5[0] == 0) && (poly5[1] == 1) && (poly5[2] == 0))
    {           // Remainder is zero.
      if (argG[0] < 0)
      {
        degree = 0;
      }
      else
      {
        degree = argG[0];
      }
      poly5[0] = argG[0];
      (void)CopyPolynomial(&poly5[1], &argG[1], degree);
      rc = DivideIntegerPolynomial(poly5, polyT, TYPE_MODULUS);
      if ((rc == EXPR_OK) && (poly5[0] == 0) && (poly5[1] == 1) && (poly5[2] == 0))
      {           // Remainder is zero.
                  // Return h*gcd(contentF, contentG)
        BigIntGcd(&contentF, &contentG, &operand1);
        *gcd = potentialDegreeGcd;
        ptrSrc = &polyT[1];
        ptrDest = gcd + 1;
        for (degree = 0; degree <= potentialDegreeGcd; degree++)
        {
          UncompressBigIntegerB(ptrSrc, &operand2);
          BigIntMultiply(&operand1, &operand2, &operand3);
          NumberLength = operand3.nbrLimbs;
          BigInteger2IntArray(ptrDest, &operand3);
          ptrSrc += 1 + numLimbs(ptrSrc);
          ptrDest += 1 + numLimbs(ptrDest);
        }
        return;
      }
    }
  }
}

// arg1 must not be changed by this routine.
// arg2 can be changed by this routine.
void PolyModularGcd(const int *arg1, int degree1, int *arg2, int degree2, int *gcd, int *degreeGcd)
{
  int nbrLimbs = NumberLength + 1;
  int *ptrArgMax = gcd;
  int *ptrArgMin = arg2;
  int *ptrTemp;
  int degreeMax = degree1;
  int degreeMin = degree2;
  int temp;
  int index;
  if (degree2 == 0)
  {
    if ((*arg2 == 1) && (*(arg2 + 1) == 0))
    {     // Number is zero. GCD is first argument.
      *degreeGcd = degree1;
      (void)memcpy(gcd, arg1, (degree1+1)*nbrLimbs*sizeof(int));
      return;
    }
    *degreeGcd = 0;
    return;
  }
  (void)memcpy(gcd, arg1, (degree1+1)*nbrLimbs*sizeof(int));
  if (degreeMax < degreeMin)
  {
    temp = degreeMax;
    degreeMax = degreeMin;
    degreeMin = temp;
    ptrTemp = ptrArgMax;
    ptrArgMax = ptrArgMin;
    ptrArgMin = ptrTemp;
  }
  for (;;)
  {
    ConvertToMonic(ptrArgMin, degreeMin);
    for (int currentDegree = degreeMax; currentDegree >= degreeMin; currentDegree--)
    {          // Replace ptrArgMax by remainder of long
               // division of ptrArgMax / ptrArgMin.
      ptrTemp = ptrArgMax + (currentDegree - degreeMin)*nbrLimbs;
      if (nbrLimbs == 2)
      {        // Modulus size is one limb.
        int modulus = TestNbr[0].x;
        int value1 = modulus - *(ptrArgMax + (currentDegree * 2) + 1);
        const int *ptrPolynomial = ptrArgMin + 1;
        ptrTemp++;
        if (modulus < 32767)
        {
          for (index = 0; index <= degreeMin; index++)
          {
            *ptrTemp = (*ptrTemp + (*ptrPolynomial * value1)) % modulus;
            ptrTemp += 2;
            ptrPolynomial += 2;
          }
        }
        else
        {
          for (index = 0; index <= degreeMin; index++)
          {
#ifdef _USING64BITS_
            *ptrTemp = (*ptrTemp + (uint64_t)*ptrPolynomial * value1) % modulus;
#else
            smallmodmult(*ptrPolynomial, value1, (limb*)&temp, modulus);
            *ptrTemp = (int)(((unsigned int)*ptrTemp + (unsigned int)temp) %
              (unsigned int)modulus);
#endif
            ptrTemp += 2;
            ptrPolynomial += 2;
          }
        }
      }
      else
      {        // General case.
        IntArray2BigInteger(ptrArgMax + (currentDegree * nbrLimbs), &operand1);
        for (index = 0; index <= degreeMin; index++)
        {
          IntArray2BigInteger(ptrArgMin + (index * nbrLimbs), &operand2);
          modmult(operand1.limbs, operand2.limbs, operand2.limbs);
          IntArray2BigInteger(ptrTemp, &operand3);
          SubtBigNbrMod(operand3.limbs, operand2.limbs, operand3.limbs);
          BigInteger2IntArray(ptrTemp, &operand3);
          ptrTemp += nbrLimbs;
        }
      }
    }
    degreeMax = degreeMin - 1;
               // Discard any leading zero coefficient.
    while (degreeMax >= 0)
    {
      ptrTemp = ptrArgMax + (degreeMax*nbrLimbs);
      for (index = *ptrTemp; index > 0; index--)
      {
        ptrTemp++;
        if (*ptrTemp != 0)
        {      // Coefficient is not zero.
          break;
        }
      }
      if (index > 0)
      {
        break;
      }
      degreeMax--;
    }
    if (degreeMax < 0)
    {          // Polynomial is zero, so Gcd was found.
      if (ptrArgMin != gcd)
      {        // Move gcd to output buffer.
        (void)memcpy(gcd, ptrArgMin, (degreeMin+1)*nbrLimbs*sizeof(int));
      }
      *degreeGcd = degreeMin;
      return;
    }
    temp = degreeMax;
    degreeMax = degreeMin;
    degreeMin = temp;
    ptrTemp = ptrArgMax;
    ptrArgMax = ptrArgMin;
    ptrArgMin = ptrTemp;
  }
}

void DerPolynomial(int *ptrArgument)
{
  const int* ptrSrc;
  int *ptrDest;
  int nbrLimbs;
  int degreePoly = *ptrArgument;
  int derivDegreePoly = 0;
  if (degreePoly < 0)
  {           // Monomial.
    *ptrArgument = degreePoly + 1;   // Decrement degree.
    if (modulusIsZero)
    {
      UncompressBigIntegerB(ptrArgument + 1, &operand1);
      multint(&operand1, &operand1, -degreePoly);
      NumberLength = operand1.nbrLimbs;
      BigInteger2IntArray(ptrArgument + 1, &operand1);
    }
    else
    {
      nbrLimbs = NumberLength + 1;
      LenAndLimbs2ArrLimbs(ptrArgument + 1, operand1.limbs, nbrLimbs);
      modmultInt(operand1.limbs, -degreePoly, operand1.limbs);
      ArrLimbs2LenAndLimbs(ptrArgument + 1, operand1.limbs, nbrLimbs);
      if ((*(ptrArgument + 1) == 1) && (*(ptrArgument + 2) == 0))
      {   // Only coefficient is zero. Set degree to zero.
        *ptrArgument = 0;
      }
    }
    return;
  }
  ptrSrc = ptrArgument + 1;
  ptrDest = &poly1[1];
  for (degree = 1; degree <= degreePoly; degree++)
  {
    ptrSrc += 1 + numLimbs(ptrSrc);
    if (modulusIsZero)
    {
      UncompressBigIntegerB(ptrSrc, &operand1);
      multint(&operand1, &operand1, degree);
      NumberLength = operand1.nbrLimbs;
      BigInteger2IntArray(ptrDest, &operand1);
      derivDegreePoly = degree - 1;
    }
    else
    {
      nbrLimbs = NumberLength + 1;
      LenAndLimbs2ArrLimbs(ptrSrc, operand1.limbs, nbrLimbs);
      modmultInt(operand1.limbs, degree, operand1.limbs);
      ArrLimbs2LenAndLimbs(ptrDest, operand1.limbs, nbrLimbs);
      if ((*ptrDest != 1) || (*(ptrDest + 1) != 0))
      {
        derivDegreePoly = degree-1;
      }
    }
    ptrDest += 1 + numLimbs(ptrDest);
  }
  (void)CopyPolynomial(ptrArgument + 1, &poly1[1], derivDegreePoly);
  *ptrArgument = derivDegreePoly;
  return;
}

// Convert from Montgomery notation to standard notation by multiplying by 1.
void polyToStandardNotation(int *nbr, int qtyNbrs)
{
  int* ptrNbr = nbr;
  // Initialize operand2.limbs to 1 (little-endian).
  operand2.limbs[0].x = 1;
  if (NumberLength > 1)
  {
    (void)memset(&operand2.limbs[1], 0, (NumberLength - 1) * sizeof(limb));
  }
  for (int currentNbr = 0; currentNbr < qtyNbrs; currentNbr++)
  {
    IntArray2BigInteger(ptrNbr, &operand1);
    modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    BigInteger2IntArray(ptrNbr, &operand1);
    ptrNbr += NumberLength + 1;
  }
}

void polyToMontgomeryNotation(int *nbr, int qtyNbrs)
{
  for (int currentNbr = qtyNbrs - 1; currentNbr >= 0; currentNbr--)
  {
    IntArray2BigInteger(nbr, &operand1);
    modmult(operand1.limbs, MontgomeryMultR2, operand1.limbs);
    BigInteger2IntArray(nbr, &operand1);
    nbr += NumberLength + 1;
  }
}

  // Perform polyPower <- polyBase ^ expon (mod polyMod)
void powerPolynomial(int *polyBase, int *polyMod, int polyDegree, const BigInteger *expon,
                     int *polyPower, powerCback callback, int curMultip, int nbrMultip)
{
  int nbrLimbs = NumberLength + 1;
  bool powerIsOne = true;
  int nbrBits = 0;
  int index = expon->nbrLimbs - 1;
  int bitCounter = 0;
  *(polyBase + (polyDegree * nbrLimbs)) = 1;
  *(polyBase + (polyDegree * nbrLimbs) + 1) = 0;
  for (; index >= 0; index--)
  {
    int groupExp = expon->limbs[index].x;
    for (int mask = 1 << (BITS_PER_GROUP - 1); mask > 0; mask >>= 1)
    {
      if (!powerIsOne)
      {
        if (callback != NULL)
        {
          callback(100 * bitCounter/(nbrBits*nbrMultip));
        }
        multUsingInvPolynomial(polyPower, polyPower, polyPower, polyDegree, polyMod);
      }
      if ((groupExp & mask) != 0)
      {
        if (powerIsOne)
        {
          (void)CopyPolynomialFixedCoeffSize(polyPower, polyBase, polyDegree, nbrLimbs);
          nbrBits = ((index + 1) * BITS_PER_GROUP) - bitCounter;
          powerIsOne = false;
          bitCounter = curMultip * nbrBits;
        }
        else
        {
          multUsingInvPolynomial(polyPower, polyBase, polyPower, polyDegree, polyMod);
        }
      }
      bitCounter++;
    }
  }
}

int getDegreePoly(const int *poly, int polyDegree)
{
  int currentDegree;
  if (polyDegree <= 0)
  {
    return 0;
  }
  currentDegree = polyDegree;
  while (currentDegree > 0)
  {
    const int *ptrTemp = poly + (currentDegree*(NumberLength + 1));
    int len = *ptrTemp;
    for (int index = len; index > 0; index--)
    {
      ptrTemp++;
      if (*ptrTemp != 0)
      {      // Coefficient is not zero.
        return currentDegree;
      }
    }
    currentDegree--;
  }
  return currentDegree;
}

#if 0
// Square-free factorization (from Wikipedia).
// i <- 1; g <- f';
// if g != 0 then {
//   c <- gcd(f, g);
//   w <- f / c;
//   while w != 1 do {
//     y <- gcd(w, c); z <- w / y;
//     Output(z^i); i <- i + 1;
//     w <- y; c <- c / y
//   }
//   if c != 1 then {
//     c <- c^(1/p);
//     Output(SFF(c)^p)
//   }
// }
// else {
//    f <- f^(1/p);
//    Output(SFF(f)^p)
//  }
//  end.
#endif
void SquareFreeFactorization(int polyDegree, int *poly, int expon)
{
  int currentDegree;
  int degreeC;
  int degreeY;
  int index;
  int primeInt = primeMod.limbs[0].x;
  int nbrLimbs = primeMod.nbrLimbs + 1;
  int* ptrValue1;
  int *ptrValue2;
  // Generate derivative in poly1.
  ptrValue1 = poly;
  ptrValue2 = poly1;
  for (currentDegree = 1; currentDegree <= polyDegree; currentDegree++)
  {
    ptrValue1 += nbrLimbs;
    IntArray2BigInteger(ptrValue1, &operand1);
    modmultInt(operand1.limbs, currentDegree, operand1.limbs);
    BigInteger2IntArray(ptrValue2, &operand1);
    ptrValue2 += nbrLimbs;
  }
  // Check whether the derivative is zero.
  for (currentDegree = polyDegree-1; currentDegree >= 0; currentDegree--)
  {
    ptrValue1 = &poly1[(currentDegree*nbrLimbs)+1];
    for (index = 1; index < nbrLimbs; index++)
    {
      if (*ptrValue1 != 0)
      {         // Coefficient of derivative is not zero.
        break;
      }
      ptrValue1++;
    }
    if (index < nbrLimbs)
    {           // Coefficient of derivative is not zero.
      break;
    }
  }
  if (currentDegree >= 0)
  {             // Derivative is not zero.
    int degreeW;
    int i = 1;

    PolyModularGcd(poly, polyDegree, poly1, currentDegree, poly2, &degreeC); // poly2 = c
    (void)memcpy(poly4, poly, (polyDegree + 1)*nbrLimbs*sizeof(int));      // Backup poly
    if (degreeC == 0)
    {
      (void)memcpy(poly1, poly, (polyDegree + 1)*nbrLimbs*sizeof(int));           // poly1 = w
    }
    else
    {
      SetNumberToOne(&poly2[degreeC*nbrLimbs]);
      DividePolynomial(poly4, polyDegree, poly2, degreeC, poly1);    // poly1 = w
    }
    degreeW = polyDegree - degreeC;
    while (degreeW != 0)
    {
      (void)memcpy(poly4, poly1, (degreeW+1)*nbrLimbs*sizeof(int));        // Backup w.
      PolyModularGcd(poly2, degreeC, poly1, degreeW, poly3, &degreeY);      // poly3 = y
      SetNumberToOne(&poly3[degreeY*nbrLimbs]);
      if (degreeW != degreeY)
      {
        int degreeZ;
        struct sFactorInfo *pstFactorInfo;

        DividePolynomial(poly4, degreeW, poly3, degreeY, poly1);     // poly1 = z
        // z^i is divisor of the original polynomial.
        degreeZ = degreeW - degreeY;
        for (currentDegree = 0; currentDegree < (i*expon); currentDegree++)
        {
          DividePolynomial(ptrOrigPoly, degreeOrigPoly, poly1, degreeZ, poly4);
          degreeOrigPoly -= degreeZ;
          (void)memcpy(ptrOrigPoly, poly4, (degreeOrigPoly+1)*nbrLimbs*sizeof(int));
        }
        pstFactorInfo = &factorInfo[nbrFactorsFound];
        nbrFactorsFound++;
        pstFactorInfo -> ptr = ptrOrigPoly;
        pstFactorInfo -> degree = degreeZ;
        pstFactorInfo -> multiplicity = i*expon;
        pstFactorInfo -> expectedDegree = 0;    // Unknown at this moment.
        (void)memcpy(ptrOrigPoly, poly1, degreeZ*nbrLimbs*sizeof(int));
        ptrOrigPoly += (degreeZ+1)*nbrLimbs;
        (void)memcpy(ptrOrigPoly, poly4, (degreeOrigPoly + 1)*nbrLimbs*sizeof(int));
      }
      i++;
      (void)memcpy(poly1, poly3, (degreeY+1)*nbrLimbs*sizeof(int));  // Copy y to w.
      degreeW = getDegreePoly(poly1, degreeY);
      DividePolynomial(poly2, degreeC, poly1, degreeW, poly4); // Compute c.
      degreeC -= degreeW;
      (void)memcpy(poly2, poly4, (degreeC+1)*nbrLimbs*sizeof(int));
    }
    if (degreeC != 0)
    {      // C is a perfect power.
      ptrValue1 = ptrOrigPoly + (degreeOrigPoly + 1)*nbrLimbs;
      ptrValue2 = poly2;
      for (currentDegree = 0; currentDegree <= degreeC; currentDegree += primeInt)
      {
        (void)memcpy(ptrValue1, ptrValue2, nbrLimbs*sizeof(int));
        ptrValue1 += nbrLimbs;
        ptrValue2 += primeInt*nbrLimbs;
      }
      SquareFreeFactorization(degreeC / primeInt, ptrOrigPoly + (degreeOrigPoly + 1)*nbrLimbs, expon * primeInt);
    }
  }
  else
  {           // Derivative is zero.
    ptrValue1 = ptrOrigPoly + (degreeOrigPoly + 1)*nbrLimbs;
    ptrValue2 = poly;
    for (currentDegree = 0; currentDegree <= polyDegree; currentDegree += primeInt)
    {
      (void)memcpy(ptrValue1, ptrValue2, nbrLimbs*sizeof(int));
      ptrValue1 += nbrLimbs;
      ptrValue2 += primeInt*nbrLimbs;
    }
    SquareFreeFactorization(polyDegree / primeInt, ptrOrigPoly + (degreeOrigPoly + 1)*nbrLimbs, expon * primeInt);
  }
}

#if 0
// Extended GCD algorithm of polynomials a, b:
// The algorithm finds polynomials u and v such that au + bv = g
// deg(u) < deg(b) - deg(g), deg(v) < deg(a) - deg(g).
//
// r <- b, old_r <- a
// for (i <- 1; r != 0; i <- i+1) do
//   q[i] = quotient(old_r, r)
//   (old_r, r) <- (r, old_r - q[i] * r)
// end do
// r <- 0, old_r <- 1
// for (j <- 1; j < i; j <- j+1) do
//   (old_r, r) <- (r, old_r - q[j] * r)
// end do
// u <- r
// r <- 1, old_r <- 0
// for (j <- 1; j < i; j <- j+1) do
//   (old_r, r) <- (r, old_r - q[j] * r)
// end do
// v <- r
#endif

// Given coprime polynomials A and B, find polynomials U and V such that
// A*U + B*V = 1 (mod powerMod). Polynomial V is not required on output.
// Use Montgomery notation.
static void ExtendedGcdPolynomial(/*@in@*/int *pointerA, int degreeA, /*@in@*/int *pointerB,
  int degreeB, /*@out@*/int *ptrQ,
  /*@out@*/int *ptrP1, /*@out@*/int *ptrP2, /*@out@*/int *ptrU, /*@out@*/int *pDegreeU)
{
  int *ptrQuotient = ptrQ;
  int *ptrQuotients[MAX_DEGREE+1];
  int* tmpPtr;
  int* ptrR;
  int *ptrOldR;
  int degreeR;
  int degreeOldR;
  int nbrQuotients = 0;
  int nbrLimbs = NumberLength + 1;
  bool polyExchanged = false;
  int currentDegree;
  int degreeU;

  // Ensure that degree of A is greater than degree of B.
  // Variable polyExchange indicates whether the polynomials were exchanged or not.
	if (degreeA < degreeB)
	{
    int tmpDegree = degreeA;
    degreeA = degreeB;
    degreeB = tmpDegree;
    tmpPtr = pointerA;
    pointerA = pointerB;
    pointerB = tmpPtr;
    polyExchanged = true;
	}
  ptrR = ptrP1;
  ptrOldR = ptrP2;
  (void)memcpy(ptrOldR, pointerA, (degreeA + 1)*nbrLimbs*sizeof(int));
  (void)memcpy(ptrR, pointerB, (degreeB + 1)*nbrLimbs*sizeof(int));
  degreeR = degreeB;
  degreeOldR = degreeA;
     // Loop that computes all quotients.
     // Since gcd = 1, the loop will finish when the degree of R is zero.
  while (degreeR != 0)
  {  // If R is not constant...
    do
    {
      IntArray2BigInteger(ptrR + (degreeR * (NumberLength+1)), &operand1);
      (void)BigIntRemainder(&operand1, &primeMod, &operand2);
      if (!BigIntIsZero(&operand2))
      {
        break;
      }
      degreeR--;
    } while (degreeR != 0);
    if (degreeR == 0)
    {
      break;
    }
    DividePolynomial(ptrOldR, degreeOldR, ptrR, degreeR, ptrQuotient);
    // ptrQuotient points to quotient of division.
    // ptrOldR points to the remainder.
    // Save pointer to the quotient.
    ptrQuotients[nbrQuotients] = ptrQuotient;
    nbrQuotients++;
    // Move pointer to new quotient after this quotient.
    ptrQuotient += (degreeOldR - degreeR + 1)*nbrLimbs;
    // Exchange pointers.
    tmpPtr = ptrR;
    ptrR = ptrOldR;
    ptrOldR = tmpPtr;
    // Get degrees of dividend and divisor of next division.
    degreeOldR = degreeR;
    degreeR = getDegreePoly(ptrR, degreeR-1);
  }
  LenAndLimbs2ArrLimbs(ptrR, operand5.limbs, nbrLimbs);
  // Save pointer to last quotient.
  ptrQuotients[nbrQuotients] = ptrQuotient;
  // P1 <- 0
  *ptrP1 = 1;
  *(ptrP1 + 1) = 0;
  // P2 <- modular inverse of gcd (stored in operand5).
  nbrLimbs = NumberLength;
  ModInvBigNbr(operand5.limbs, (limb*)(ptrP2 + 1), powerMod.limbs, powerMod.nbrLimbs);
  while ((nbrLimbs > 1) && (*(ptrP2 + nbrLimbs) == 0))
  {
    nbrLimbs--;
  }
  *ptrP2 = nbrLimbs;
  if (polyExchanged)
  {
    ptrR = ptrP2;
    ptrOldR = ptrP1;
  }
  else
  {
    ptrR = ptrP1;
    ptrOldR = ptrP2;
  }
  degreeR = 0;
  degreeOldR = 0;
  nbrLimbs = NumberLength + 1;
  for (int counter = 0; counter < nbrQuotients; counter++)
  {
    int tmpDegree;
    int offset;
    int degreeQ = (int)(((ptrQuotients[counter + 1] - ptrQuotients[counter]) / nbrLimbs) - 1);
    // Multiply ptrR by ptrQuotients[counter]. The result will be stored in polyMultTemp.
    MultPolynomial(degreeR, degreeQ, ptrR, ptrQuotients[counter]);
    // Compute ptrOldR <- ptrOldR - polyMultTemp.
    degreeQ += degreeR; // Degree of product (always greater than degree of old R).
    offset = 0;
    for (currentDegree = 0; currentDegree <= degreeOldR; currentDegree++)
    {
      IntArray2BigInteger(ptrOldR + offset, &operand1);
      IntArray2BigInteger(&polyMultTemp[offset], &operand2);
      SubtBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
      BigInteger2IntArray(ptrOldR + offset, &operand1);
      offset += nbrLimbs;
    }
    for (; currentDegree <= degreeQ; currentDegree++)
    {
      // Set operand1 to zero.
      (void)memset(operand1.limbs, 0, NumberLength * sizeof(limb));
      IntArray2BigInteger(&polyMultTemp[offset], &operand2);
      SubtBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
      operand1.sign = SIGN_POSITIVE;
      BigInteger2IntArray(ptrOldR + offset, &operand1);
      offset += nbrLimbs;
    }
    degreeOldR = getDegreePoly(ptrOldR, degreeQ);
    // Exchange degrees.
    tmpDegree = degreeR;
    degreeR = degreeOldR;
    degreeOldR = tmpDegree;
    // Exchange pointers.
    tmpPtr = ptrR;
    ptrR = ptrOldR;
    ptrOldR = tmpPtr;
  }
  degreeU = degreeR;
  (void)memcpy(ptrU, ptrR, (degreeU + 1) * nbrLimbs * sizeof(int));
  *pDegreeU = degreeU;
}

// Lift factorization of coefficients of factors from mod prime to
// mod power = prime^exp.
// Perform lifting mod prime, prime^2, prime^4, etc.
// All calculations are done in Montgomery domain, so the roots must 
// be converted to standard notation and then to Montgomery notation
// with greater modulus.
//
// Let f(x) = f_1(x) * f_2(x) * ... * f_n(x)
// 
// Compute 1 = (f(x)/f_1(x)) * a1(x) + ... + (f(x)/f_n(x)) * a_n(x)
// This is done by computing n extended GCDs between f_i(x) and
// f(x)/f_i(x).
//
// Hensel lift from mod m to mod m^2:
// u(x) <- (1/m)*(f - f_1*f_2*...*f_n)
// For i = 1 to n compute:
//   g(x) <- u(x) * a_i(x) mod f_i(x) (all computations done mod m).
//   f_i(x) <- f_i(x) + m*g(x)
//
// Compute u(x) <- (1/m)*(f(x)/f_1(x)*a_1(x) + ... + f(x)/f_n(x)*a_n(x) - 1).
// For i = 1 to n compute:
//   g(x) <- u(x) * a_i(x) mod f_i(x) (all computations done mod m).
//   a_i(x) <- a_i(x) - m*g(x)
//
// Memory usage:
// poly3: u
// poly4: f, g
// poly5: m*g, f/f_i*a_i
// polyT: a_i
// polyLifted: f_i.
// Change the fields ptrPolyLifted to point to f_i.

// Compute polynomial poly4 <- values divided by the leading coefficient
// mod prime. This is the polynomial f.
static void ComputeF(void)
{
  int currentDegree;
  const int* ptrValue1;
  int polyDegree = values[0];                      // Get degree of polynomial
  // Find the leading coefficient.
  ptrValue1 = &values[1];                          // Point to constant coefficient.
  for (currentDegree = 0; currentDegree < polyDegree; currentDegree++)
  {
    ptrValue1 += 1 + numLimbs(ptrValue1);
  }
  NumberLength = numLimbs(ptrValue1);
  IntArray2BigInteger(ptrValue1, &operand1);
  (void)BigIntRemainder(&operand1, &powerMod, &operand1);
  NumberLength = powerMod.nbrLimbs;
  if (operand1.nbrLimbs < NumberLength)
  {
    (void)memset(&operand1.limbs[operand1.nbrLimbs], 0,
      (NumberLength - operand1.nbrLimbs) * sizeof(limb));
  }
  // Compute the inverse of leading coefficient.
  ModInvBigNbr(operand1.limbs, operand2.limbs, TestNbr, NumberLength);
  // Convert operand1 from standard to Montgomery notation.
  ptrValue1 = &values[1];                          // Point to constant coefficient.
  for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {
    NumberLength = numLimbs(ptrValue1);
    IntArray2BigInteger(ptrValue1, &operand1);
    (void)BigIntRemainder(&operand1, &powerMod, &operand1);
    NumberLength = powerMod.nbrLimbs;
    if (operand1.nbrLimbs < NumberLength)
    {
      (void)memset(&operand1.limbs[operand1.nbrLimbs], 0,
        (NumberLength - operand1.nbrLimbs) * sizeof(limb));
    }
    // Multiply by inverse. Result is in Montgomery notation.
    modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    BigInteger2IntArray(&poly4[currentDegree * (NumberLength + 1)], &operand1);
    ptrValue1 += 1 + numLimbs(ptrValue1);           // Point to next coefficient.
  }
}

// Get polynomial a_i. This polynomial is stored with no spaces
// between coefficients. The output is in Montgomery notation in
// polynomial poly2.
static int getAi(int nbrFactor, int degreeA)
{
  int nbrLimbs = NumberLength+1;
  const int* ptrSrc = ptrA[nbrFactor];
  int* ptrDest = poly1;
  for (int currentDegree = 0; currentDegree <= degreeA; currentDegree++)
  {
    int numLen = *ptrSrc + 1;
    (void)memcpy(ptrDest, ptrSrc, numLen * sizeof(int));
    ptrDest += nbrLimbs;
    ptrSrc += numLen;
  }
  polyToMontgomeryNotation(poly1, degreeA+1);
  return degreeA;
}

// Get polynomial factor f_i. This factor is stored with no spaces
// between coefficients. The output is in Montgomery notation in
// polynomial poly2.
static void getFi(int degreeFactor, const int *src, int nbrLimbs)
{
  const int* ptrSrc = src;
  int *ptrDest = poly2;       // Copy f_i to poly2.
  for (int currentDegree = 0; currentDegree < degreeFactor; currentDegree++)
  {
    int nbrLen = *ptrSrc + 1;
    (void)memcpy(ptrDest, ptrSrc, nbrLen * sizeof(int));
    ptrSrc += nbrLen;
    ptrDest += nbrLimbs;
  }
  SetNumberToOne(ptrDest);
  polyToMontgomeryNotation(poly2, degreeFactor);
}

// Move factors f_i from polyLiftedNew to polyLifted, then adjust
// factorInfo[...].ptrPolyLifted to point to the new factors.
// There are no spaces between coefficients.
static void MoveFactorsAndFixPointers(struct sFactorInfo* ptrFactorInfo, bool compressPoly)
{
  struct sFactorInfo* pstFactorInfo = ptrFactorInfo;
  const int* ptrSrc = polyLiftedNew;
  int * ptrDest = polyLifted;
  for (int nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    int degreeFactor = pstFactorInfo->degree;
    pstFactorInfo->ptrPolyLifted = ptrDest;
    for (int currentDegree = 0; currentDegree < degreeFactor; currentDegree++)
    {
      int numLength = *ptrSrc + 1;
      (void)memcpy(ptrDest, ptrSrc, numLength * sizeof(int));
      ptrSrc += numLength;
      if (compressPoly)
      {
        ptrDest += numLength;
      }
      else
      {
        ptrDest += NumberLength+1;
      }
    }
    pstFactorInfo++;
  }
}

// Perform Hensel lifting of polynomials in factorInfo.
// On entry:
// values = original polynomial. The product of all lifted factors
// will be equal to this polynomial.
// exponentMod = exponent to perform the lift.
// primeMod = prime as a big integer.
// nbrFactorsFound = number of factors of original polynomial.
// factorInfo = pointer to array of structures holding the polynomial factors.
// Hensel lifting only works when there are no repeated factors mod p,
// so this condition is tested.

int HenselLifting(struct sFactorInfo* ptrFactorInfo, bool compressPoly)
{
  int currentExp = 1;
  int nbrFactor;
  int oldNumberLength;
  int newNumberLength;
  int degreeFactor;
  int degreeS;
  int currentDegree;
  int nbrLimbs;
  int* ptrPolyLifted = polyLifted;
  int* ptrPoly2;
  int** ptrNewFactor;
  struct sFactorInfo *pstFactorInfo;
  // Copy polynomials f_i(x) to polyLifted.
  int* ptrDest = ptrPolyLifted;
  pstFactorInfo = ptrFactorInfo;
  for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    pstFactorInfo->ptrPolyLifted = ptrDest;
    ptrDest = CopyPolynomialFixedCoeffSize(ptrDest, pstFactorInfo->ptr,
      pstFactorInfo->degree, NumberLength+1);
    pstFactorInfo++;
  }
  if (exponentMod == 1)
  {      // No lift to be done. Go out.
    return EXPR_OK;
  }
  CopyBigInt(&powerMod, &primeMod);
  // Hensel lifting only works when there are no repeated factors mod p.
  // Go out if the multiplicity of a factor is not 1.
  // Compute in variable degree the degree of the product.
  degree = 0;
  ptrPoly2 = poly2;
  pstFactorInfo = ptrFactorInfo;
  for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    if (pstFactorInfo->multiplicity != 1)
    {   // Repeated factor found. Indicate error.
      return EXPR_CANNOT_LIFT;
    }
    // Init poly2 coefficients to zero.
    for (int index = 0; index < pstFactorInfo->degree; index++)
    {
      *ptrPoly2++ = 1;
      *ptrPoly2 = 0;
      ptrPoly2 += NumberLength;
    }
    degree += pstFactorInfo->degree;
    pstFactorInfo++;
  }
  computePower(1);
  // Compute a_i for each factor f_i by computing the extended gcd between
  // f_i(x) and f(x)/f_i(x).
  ptrDest = polyT;
  pstFactorInfo = ptrFactorInfo;
  for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    ComputeF();                       // poly4 <- f mod prime (Montgomery notation).
    nbrLimbs = NumberLength + 1;
    degreeFactor = pstFactorInfo->degree;
    (void)memcpy(poly2, pstFactorInfo->ptr, pstFactorInfo->degree * nbrLimbs * sizeof(int));
    // The factors do not include their leading coefficient. Set it to one.
    SetNumberToOne(&poly2[pstFactorInfo->degree * nbrLimbs]);
    DividePolynomial(poly4, degree,   // Dividend f(x).
      poly2, degreeFactor,            // Divisor f_i(x).
      poly1);                         // Quotient f(x)/f_i(x).
    ExtendedGcdPolynomial(poly1, degree - degreeFactor, poly2, degreeFactor,
      poly3, poly4, poly5, polyS, &degreeS);
    ptrA[nbrFactor] = ptrDest;
    ptrDest = CopyPolynomial(ptrDest, polyS, degreeS);
    for (currentDegree = degreeS + 1; currentDegree < degreeFactor; currentDegree++)
    {
      *ptrDest = 1;
      ptrDest++;
      *ptrDest = 0;
      ptrDest++;
    }
    pstFactorInfo++;
  }
  // Loop that performs the lifting.
  while (currentExp != exponentMod)
  {      // Hensel lift from mod m to mod m^2:
    int newExponent;
    const int* ptrSrc;
    int degreeA;
#ifdef __EMSCRIPTEN__
    int elapsedTime = (int)(tenths() - originalTenthSecond);
    if (elapsedTime / 10 != oldTimeElapsed / 10)
    {
      char outputInfo[1000];
      char* ptrOutput = outputInfo;
      if (lang)
      {
        copyStr(&ptrOutput, "1<p>Aplicando lema de Hensel usando el número primo ");
        int2dec(&ptrOutput, primeMod.limbs[0].x);
        copyStr(&ptrOutput, " procesando exponente ");
        int2dec(&ptrOutput, currentExp);
        copyStr(&ptrOutput, " de ");
      }
      else
      {
        copyStr(&ptrOutput, "1<p>Hensel lifting using prime number ");
        int2dec(&ptrOutput, primeMod.limbs[0].x);
        copyStr(&ptrOutput, " processing exponent ");
        int2dec(&ptrOutput, currentExp);
        copyStr(&ptrOutput, " of ");
      }
      int2dec(&ptrOutput, exponentMod);
      copyStr(&ptrOutput, ".</p>");
      showElapsedTimeSec(&ptrOutput);
      databack(outputInfo);
    }
#endif
    // Compute u(x) <- (1/m) * (f - f_1 * f_2 *...* f_n)
    oldNumberLength = NumberLength;
    newExponent = currentExp << 1;     // We can double the exponent in each step.
    if (newExponent > exponentMod)     // Do not exceed exponentMod.
    {
      newExponent = exponentMod;
    }
    (void)BigIntPowerIntExp(&primeMod, currentExp, &operand5);
    computePower(newExponent);         // Compute powerMod and init Montgomery parms.
    newNumberLength = NumberLength;
    nbrLimbs = NumberLength + 1;
    ComputeF();                        // poly4 <- f mod m^2 (Montgomery notation).
    SetNumberToOne(poly1);             // Initialize product of factors.
    pstFactorInfo = ptrFactorInfo;
    degree = 0;
    for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
    {    // Loop that multiplies all factors. poly1 holds the result.
      degreeFactor = pstFactorInfo->degree;
      getFi(degreeFactor, pstFactorInfo->ptrPolyLifted, nbrLimbs);   // Copy f_i to poly2.
      MultPolynomial(degree, degreeFactor, poly1, poly2);
      degree += degreeFactor;
      (void)memcpy(poly1, polyMultTemp, (degree + 1) * nbrLimbs * sizeof(int));
      pstFactorInfo++;
    }
    polyToStandardNotation(poly1, degree + 1);
    polyToStandardNotation(poly4, degree + 1);
    for (currentDegree = 0; currentDegree <= degree; currentDegree++)
    {                 // Loop that computes (1/m)*(f - f_1 * f_2 * ... * f_n)
      int nbrLen;
      // Get coefficient of f.
      IntArray2BigInteger(&poly4[currentDegree * nbrLimbs], &operand1);
      // Get coefficient of product of factors.
      IntArray2BigInteger(&poly1[currentDegree * nbrLimbs], &operand2);
      SubtBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
      // Get number of significant limbs before performing division.
      nbrLen = NumberLength;
      while (nbrLen > 1)
      {
        if (operand1.limbs[nbrLen - 1].x != 0)
        {
          break;
        }
        nbrLen--;
      }
      operand1.nbrLimbs = nbrLen;
      (void)BigIntDivide(&operand1, &operand5, &operand1);
      // Store coefficient of subtraction.
      ptrDest = &poly3[currentDegree * (oldNumberLength + 1)];
      *ptrDest = operand1.nbrLimbs;
      ptrDest++;
      (void)memcpy(ptrDest, operand1.limbs, operand1.nbrLimbs * sizeof(int));
    }
    computePower(currentExp);
    polyToMontgomeryNotation(poly3, degree+1);
    nbrLimbs = NumberLength + 1;
    pstFactorInfo = ptrFactorInfo;
    ptrDest = polyLiftedNew;
    ptrNewFactor = ptrNewFactors;
    for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
    {
      degreeFactor = pstFactorInfo->degree;
      // g(x) <- u(x) * a_i(x) mod f_i(x) (all computations done mod m).
      degreeA = getAi(nbrFactor, degreeFactor - 1); // poly1 <- a_i.
      MultPolynomial(degree, degreeA, poly3, poly1);
      (void)memcpy(poly1, polyMultTemp, (degree + degreeA + 1)* nbrLimbs * sizeof(int));
                                                 // poly1 <- u * a_i
                                                 // Copy f_i to poly2.
      getFi(degreeFactor, pstFactorInfo->ptrPolyLifted, nbrLimbs);
      DividePolynomial(poly1, degree + degreeA,  // Dividend = u * a_i
        poly2, degreeFactor,                     // Divisor = f_i
        NULL);                                   // Quotient: not needed.
      polyToStandardNotation(poly1, degreeFactor);  // Convert g to standard notation.
      // f_i(x) <- f_i(x) + m*g(x)
      const int *ptrFi = pstFactorInfo->ptrPolyLifted;
      *ptrNewFactor = ptrDest;                   // Point to f_i(x) mod new modulus.
      ptrNewFactor++;
      for (currentDegree = 0; currentDegree < degreeFactor; currentDegree++)
      {                 // Loop that computes f_i + m*g
        // Get coefficient of g.
        ptrSrc = &poly1[currentDegree * nbrLimbs];
        operand1.nbrLimbs = *ptrSrc;
        ptrSrc++;
        (void)memcpy(operand1.limbs, ptrSrc, operand1.nbrLimbs * sizeof(int));
        if (currentExp * 2 > newExponent)
        {
          (void)BigIntPowerIntExp(&primeMod, newExponent - currentExp, &operand5);
          BigIntRemainder(&operand1, &operand5, &operand1);
        }
        // Get coefficient of f_i.
        operand2.nbrLimbs = *ptrFi;
        ptrFi++;
        (void)memcpy(operand2.limbs, ptrFi, operand2.nbrLimbs * sizeof(int));
        ptrFi += operand2.nbrLimbs;
        BigIntMultiply(&operand1, &powerMod, &operand1);  // poly1 <- m*g
        BigIntAdd(&operand1, &operand2, &operand1);       // poly1 <- f_i + m*g
        // Store coefficient of new f_i (no spaces between coefficients).
        *ptrDest = operand1.nbrLimbs;
        ptrDest++;
        (void)memcpy(ptrDest, &operand1.limbs, operand1.nbrLimbs*sizeof(int));
        ptrDest += operand1.nbrLimbs;    // Point to next coeffficient of f_i.
      }
      pstFactorInfo++;
    }
    computePower(newExponent);         // Compute powerMod and init Montgomery parms.
    if (newExponent == exponentMod)
    {                 // Final values of f_1, f_2,..., f_n found. Exit loop.
      break;
    }
    // Compute u(x) <- (1/m)*(f(x)/f_1(x)*a_1(x) + ... + f(x)/f_n(x)*a_n(x) - 1).
    // Init u(x) to zero.
    ptrDest = poly3;
    nbrLimbs = newNumberLength + 1;
    for (currentDegree = 0; currentDegree <= degree; currentDegree++)
    {
      *ptrDest = 1;
      *(ptrDest + 1) = 0;
      ptrDest += nbrLimbs;
    }
    pstFactorInfo = ptrFactorInfo;
    ptrNewFactor = ptrNewFactors;
    ComputeF();              // poly4 <- f (mod m^2) in Montgomery notation.
    for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
    {
      degreeFactor = pstFactorInfo->degree;
      (void)memcpy(poly1, poly4, (degree + 1) * nbrLimbs * sizeof(int));
                                        // poly2 <- f_i(x) in Montgomery notation.
      getFi(degreeFactor, *ptrNewFactor, nbrLimbs);
      ptrNewFactor++;
      DividePolynomial(poly1, degree,   // Dividend = f.
        poly2, degreeFactor,            // Divisor = f_i.
        poly5);                         // Quotient = f / f_i.
      degreeA = getAi(nbrFactor, degreeFactor - 1);       // poly1 <- a_i.
                                        // polyMultTemp <- (f / f_i) * a_i
      MultPolynomial(degreeA, degree - degreeFactor, poly1, poly5);
      ptrDest = poly3;
      for (currentDegree = 0; currentDegree <= degree - degreeFactor + degreeA; currentDegree++)
      {
        // Get coefficient of sum.
        IntArray2BigInteger(ptrDest, &operand1);
        // Get coefficient of (f / f_i) * a_i.
        IntArray2BigInteger(&polyMultTemp[currentDegree * nbrLimbs], &operand2);
        // Add them.
        AddBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
        // Save coefficient of sum.
        BigInteger2IntArray(ptrDest, &operand1);
        ptrDest += nbrLimbs;
      }
      pstFactorInfo++;
    }
    polyToStandardNotation(poly3, degree + 1);
    // Divide sum by m.
    computePower(currentExp);
    ptrDest = poly3;
    ptrSrc = poly3;
    for (currentDegree = 0; currentDegree <= degree; currentDegree++)
    { // Get coefficient of u.
      operand1.nbrLimbs = *ptrSrc;
      (void)memcpy(operand1.limbs, ptrSrc + 1, operand1.nbrLimbs * sizeof(int));
      // Divide coefficient by m.
      (void)BigIntDivide(&operand1, &powerMod, &operand1);
      // Store coefficient of quotient.
      *ptrDest = operand1.nbrLimbs;
      ptrDest++;
      (void)memcpy(ptrDest, operand1.limbs, operand1.nbrLimbs * sizeof(int));
      ptrDest += oldNumberLength;
      ptrSrc += nbrLimbs;
    }
    polyToMontgomeryNotation(poly3, degree);   // poly3 <- u(x)
    nbrLimbs = NumberLength + 1;
    pstFactorInfo = ptrFactorInfo;
    // Use polyS as a temporary storage for a_i.
    ptrDest = polyS;
    for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
    {
      degreeFactor = pstFactorInfo->degree;
      // g(x) <- u(x) * a_i(x) mod f_i(x) (all computations done mod m).
      NumberLength = oldNumberLength;
      degreeA = getAi(nbrFactor, degreeFactor - 1);  // poly1 <- a_i
      MultPolynomial(degree, degreeA, poly3, poly1);
      (void)memcpy(poly4, polyMultTemp, (degree + degreeA + 1)* nbrLimbs * sizeof(int));
                                                // poly4 <- u * a_i
                                                // poly2 <- f_i(x) in Montgomery notation. 
      getFi(degreeFactor, pstFactorInfo->ptrPolyLifted, nbrLimbs);
      DividePolynomial(poly4, degree + degreeA, // Dividend = u * a_i
        poly2, degreeFactor,                    // Divisor = f_i
        NULL);                                  // Quotient: not needed.
      polyToStandardNotation(poly1, degreeA+1); // Convert a_i to standard notation.
                                                // Convert g to standard notation.
      polyToStandardNotation(poly4, degree + degreeA - degreeFactor + 1);
      // a_i(x) <- a_i(x) - m*g(x)
      for (currentDegree = 0; currentDegree <= degreeA; currentDegree++)
      {                 // Loop that computes a_i - m*g
        // Get coefficient of g.
        ptrSrc = &poly4[currentDegree * nbrLimbs];
        operand1.nbrLimbs = *ptrSrc;
        ptrSrc++;
        (void)memcpy(operand1.limbs, ptrSrc, operand1.nbrLimbs * sizeof(int));
        if (!BigIntIsZero(&operand1))                     // g <- -g (mod m)
        {
          BigIntSubt(&powerMod, &operand1, &operand1);
        }
        // Get coefficient of a_i.
        ptrSrc = &poly1[currentDegree * nbrLimbs];
        operand2.nbrLimbs = *ptrSrc++;
        (void)memcpy(operand2.limbs, ptrSrc, operand2.nbrLimbs * sizeof(int));
        BigIntMultiply(&operand1, &powerMod, &operand1);  // poly1 <- m*g
        BigIntAdd(&operand1, &operand2, &operand1);       // poly1 <- a_i + m*g
        // Store coefficient of new a_i.
        *ptrDest = operand1.nbrLimbs;
        ptrDest++;
        (void)memcpy(ptrDest, operand1.limbs, operand1.nbrLimbs * sizeof(int));
        ptrDest += operand1.nbrLimbs;  // Point to next coefficient of a_i.
      }
      pstFactorInfo++;
    }
    currentExp = newExponent;
    NumberLength = newNumberLength;
    // Move a_i from polyS to polyT.
    ptrSrc = polyS;
    ptrDest = polyT;
    pstFactorInfo = ptrFactorInfo;
    for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
    {
      degreeA = pstFactorInfo -> degree - 1;
      ptrA[nbrFactor] = ptrDest;
      for (currentDegree = 0; currentDegree <= degreeA; currentDegree++)
      {
        nbrLimbs = *ptrSrc;
        ptrSrc++;
        *ptrDest = nbrLimbs;
        ptrDest++;
        (void)memcpy(ptrDest, ptrSrc, nbrLimbs * sizeof(int));
        ptrDest += nbrLimbs;
        ptrSrc += nbrLimbs;
      }
      pstFactorInfo++;
    }
    MoveFactorsAndFixPointers(ptrFactorInfo, true);   // Output is compressed.
  }
  MoveFactorsAndFixPointers(ptrFactorInfo, compressPoly);
  return EXPR_OK;
}

void OrigPolyFromMontgomeryToStandard(void)
{
  const int* ptrValue1;
  int *ptrValue2;
  (void)memcpy(&TestNbr, powerMod.limbs, powerMod.nbrLimbs*sizeof(limb));
  NumberLength = powerMod.nbrLimbs;
  TestNbr[NumberLength].x = 0;
  GetMontgomeryParms(powerMod.nbrLimbs);
  degree = values[0];
  ptrValue1 = &values[1];
  ptrValue2 = &poly4[0];
  for (int currentDegree = 0; currentDegree <= degree; currentDegree++)
  {
    IntArray2BigInteger(ptrValue1, &operand1);
    // Convert operand1 from Montgomery to standard notation.
    operand2.limbs[0].x = 1;
    if (NumberLength > 1)
    {
      (void)memset(&operand2.limbs[1], 0, (NumberLength - 1)*sizeof(limb));
    }
    modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    BigInteger2IntArray(ptrValue2, &operand1);
    ptrValue1 += 1 + *ptrValue1;
    ptrValue2 += 1 + *ptrValue2;
  }
  (void)memcpy(&values[1], &poly4[0], (ptrValue2 - &poly4[0])*sizeof(int));
}

static void showPower(char **pptrOutput, int exponent)
{
  char *ptrOutput = *pptrOutput;
  if (pretty == PRETTY_PRINT)
  {
    *ptrOutput = '<';
    ptrOutput++;
    *ptrOutput = 's';
    ptrOutput++;
    *ptrOutput = 'u';
    ptrOutput++;
    *ptrOutput = 'p';
    ptrOutput++;
    *ptrOutput = '>';
    ptrOutput++;
    int2dec(&ptrOutput, exponent);
    *ptrOutput = '<';
    ptrOutput++;
    *ptrOutput = '/';
    ptrOutput++;
    *ptrOutput = 's';
    ptrOutput++;
    *ptrOutput = 'u';
    ptrOutput++;
    *ptrOutput = 'p';
    ptrOutput++;
    *ptrOutput = '>';
    ptrOutput++;
  }
  else if (pretty == TEX)
  {
    *ptrOutput = '^';
    ptrOutput++;
    *ptrOutput = '{';
    ptrOutput++;
    int2dec(&ptrOutput, exponent);
    *ptrOutput = '}';
    ptrOutput++;
  }
  else
  {
    *ptrOutput = '^';
    ptrOutput++;
    int2dec(&ptrOutput, exponent);
  }
  *pptrOutput = ptrOutput;
}

void showPowerX(char **pptrOutput, int polyDegree)
{
  char *ptrOutput = *pptrOutput;
  if (polyDegree == 0)
  {
    *ptrOutput = '1';
    ptrOutput++;
  }
  else
  {
    if (pretty == PRETTY_PRINT)
    {
      copyStr(&ptrOutput, lang?"<span class=\"hide\">equis </span><span aria-hidden=\"true\"><var>x</var></span>":"<var>x</var>");
    }
    else
    {
      *ptrOutput = 'x';
      ptrOutput++;
    }
    if (polyDegree != 1)
    {
      showPower(&ptrOutput, polyDegree);
    }
  }
  *pptrOutput = ptrOutput;
}

static void showPolynomial(char **pptrOutput, const int *ptrPoly, int polyDegree, int groupLength)
{
  int nbrLimbs = NumberLength + 1;
  int currentDegree;
  char *ptrOutput = *pptrOutput;
  int *ptrIndex;
  int indexes[MAX_DEGREE+1];
  int NumberLengthBak = NumberLength;

  ptrIndex = &indexes[0];
  if (modulusIsZero)
  {
    // Fill indexes to start of each coefficient.
    int index = 0;
    for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
    {
      *ptrIndex = index;
      ptrIndex++;
      index += numLimbs(ptrPoly + index) + 1;
    }
  }
  else
  {
    NumberLength = powerMod.nbrLimbs;
    nbrLimbs = NumberLength + 1;
    for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
    {
      *ptrIndex = currentDegree * nbrLimbs;
      ptrIndex++;
    }
  }
  for (currentDegree = polyDegree - 1; currentDegree >= 0; currentDegree--)
  {
    const int* ptrValue1 = ptrPoly + indexes[currentDegree];
    int len = numLimbs(ptrValue1);
    if ((len != 1) || (*(ptrValue1 + 1) != 0))
    {            // Coefficient is not zero.
      *ptrOutput = ' ';
      ptrOutput++;
      if (*ptrValue1 > 0)
      {
        *ptrOutput = '+';
        ptrOutput++;
      }
      else
      {
        copyStr(&ptrOutput, pretty == PRETTY_PRINT? "&minus;": "-");
      }
      *ptrOutput = ' ';
      ptrOutput++;
      if ((len != 1) || (*(ptrValue1 + 1) != 1))
      {            // Absolute value of coefficient is not one.
        NumberLength = numLimbs(ptrValue1);
        IntArray2BigInteger(ptrValue1, &operand1);
        operand1.sign = SIGN_POSITIVE;
        Bin2Dec(operand1.limbs, ptrOutput, operand1.nbrLimbs, groupLength);
        ptrOutput += strlen(ptrOutput);
        if (currentDegree > 0)
        {
          if (pretty == PRETTY_PRINT)
          {
            copyStr(&ptrOutput, "&#8290;");
          }
          else if (pretty == PARI_GP)
          {
            *ptrOutput = '*';
            ptrOutput++;
          }
          showPowerX(&ptrOutput, currentDegree);
        }
      }
      else
      {
        showPowerX(&ptrOutput, currentDegree);
      }
    }
  }
  *pptrOutput = ptrOutput;
  NumberLength = NumberLengthBak;
}

void outputOriginalPolynomial(char* ptrOutput, int groupLength)
{
  int currentDegree;
  const int* ptrValue1;
  int nbrLimbs = powerMod.nbrLimbs + 1;
  degree = values[0];
  ptrValue1 = &values[1];
  if (!modulusIsZero)
  {
    int* ptrValue2;
    // Output polynomial to factor. First move polynomial to poly4
    ptrValue2 = &poly4[0];
    for (currentDegree = 0; currentDegree <= degree; currentDegree++)
    {
      (void)memcpy(ptrValue2, ptrValue1, (1 + *ptrValue1) * sizeof(int));
      ptrValue1 += 1 + *ptrValue1;
      ptrValue2 += nbrLimbs;
    }
    // Get leading coefficient.
    IntArray2BigInteger(&poly4[degree * nbrLimbs], &operand1);
  }
  else
  { // Find leading coefficient.
    for (currentDegree = 0; currentDegree < degree; currentDegree++)
    {
      ptrValue1 += 1 + numLimbs(ptrValue1);
    }
    NumberLength = numLimbs(ptrValue1);
    IntArray2BigInteger(ptrValue1, &operand1);
  }
  if (operand1.sign == SIGN_NEGATIVE)
  {
    copyStr(&ptrOutput, " &minus;");
  }
  if ((operand1.nbrLimbs != 1) || (operand1.limbs[0].x != 1) || (degree == 0))
  {     // Leading coefficient is not 1 or degree is zero.
    Bin2Dec(operand1.limbs, ptrOutput, operand1.nbrLimbs, groupLength);
    ptrOutput += strlen(ptrOutput);
  }
  if (degree > 0)
  {
    showPowerX(&ptrOutput, degree);
  }
  showPolynomial(&ptrOutput, (modulusIsZero ? &values[1] : poly4), degree, groupLength);
  if (!modulusIsZero)
  {
    copyStr(&ptrOutput, " (mod ");
    Bin2Dec(primeMod.limbs, ptrOutput, primeMod.nbrLimbs, groupLength);
    ptrOutput += strlen(ptrOutput);
    if (exponentMod != 1)
    {
      showPower(&ptrOutput, exponentMod);
    }
    *ptrOutput = ')';
    ptrOutput++;
  }
  *ptrOutput = 0;    // Append string terminator.
}

void outputPolynomialFactor(char *ptrOutput, int groupLength, const struct sFactorInfo* pstFactorInfo)
{
  int polyDegree = pstFactorInfo->degree;
  int multiplicity = pstFactorInfo->multiplicity;
  int isMonomial = ((polyDegree == 1) && (*pstFactorInfo->ptrPolyLifted == 1) &&
     (*(pstFactorInfo->ptrPolyLifted+1) == 0));
  if ((multiplicity > 1) && !isMonomial)
  {
    *ptrOutput = '(';
    ptrOutput++;
  }
  if (modulusIsZero)
  {
    // Get leading coefficient.
    const int *ptrSrc = pstFactorInfo->ptrPolyLifted;
    for (int currentDegree = 0; currentDegree < polyDegree; currentDegree++)
    {
      ptrSrc += 1 + numLimbs(ptrSrc);
    }
    UncompressBigIntegerB(ptrSrc, &operand1);
    if (operand1.sign == SIGN_NEGATIVE)
    {
      copyStr(&ptrOutput, "&minus;");
    }
    if ((operand1.nbrLimbs != 1) || (operand1.limbs[0].x != 1))
    {     // Absolute value is not 1.
      Bin2Dec(operand1.limbs, ptrOutput, operand1.nbrLimbs, groupLength);
      ptrOutput += strlen(ptrOutput);
    }
  }
  showPowerX(&ptrOutput, polyDegree);
  showPolynomial(&ptrOutput, pstFactorInfo->ptrPolyLifted, polyDegree, groupLength);
  if (multiplicity > 1)
  {
    if (!isMonomial)
    {
      *ptrOutput = ')';
      ptrOutput++;
    }
    showPower(&ptrOutput, multiplicity);
  }
  *ptrOutput = 0;    // Append string terminator.
}

void textErrorPol(char *ptrOutput, enum eExprErr rc)
{
  char text[150];

  switch (rc)
  {
  case EXPR_CANNOT_USE_X_IN_EXPONENT:
    (void)strcpy(text, lang?"No se puede usar variable en el exponente":
                      "Cannot use variable in exponent");
    break;
  case EXPR_POLYNOMIAL_DIVISION_NOT_INTEGER:
    (void)strcpy(text, lang ? "La división de polinomios no es entera" :
                        "Polynomial division is not integer");
    break;
  case EXPR_DEGREE_TOO_HIGH:
    (void)strcpy(text, lang?"El grado del polinomio es muy elevado":
                      "Degree is too high");
    break;
  case EXPR_EXPONENT_TOO_LARGE:
    (void)strcpy(text, lang?"Exponente muy grande":"Exponent is too large");
    break;
  case EXPR_EXPONENT_NEGATIVE:
    (void)strcpy(text, lang?"Exponente negativo":"Exponent is negative");
    break;
  case EXPR_LEADING_COFF_MULTIPLE_OF_PRIME:
    (void)strcpy(text, lang?"El primer coeficiente es múltiplo del número primo":
      "Leading coefficient multiple of prime");
    break;
  case EXPR_CANNOT_LIFT:
    (void)strcpy(text, lang?"No se puede elevar porque hay factores duplicados":
      "Cannot lift because of duplicate factors modulo prime");
    break;
  case EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE:
    (void)strcpy(text, lang?"El módulo debe ser mayor que 1":"Modulus must be greater than one");
    break;
  case EXPR_MODULUS_MUST_BE_PRIME_EXP:
    (void)strcpy(text, lang ? "El módulo debe ser un número primo o una potencia de número primo" :
      "Modulus must be a prime number or a power of a prime");
    break;
  case EXPR_MULTIPLE_VARIABLES_NOT_ACCEPTED:
    (void)strcpy(text, lang ? "No se aceptan múltiples variables" :
      "Multiple variables are not accepted");
    break;
  default:
    textError(text, rc);
  }
  *ptrOutput = '<';
  ptrOutput++;
  *ptrOutput = 'p';
  ptrOutput++;
  *ptrOutput = '>';
  ptrOutput++;
  copyStr(&ptrOutput, text);
  *ptrOutput = '<';
  ptrOutput++;
  *ptrOutput = '/';
  ptrOutput++;
  *ptrOutput = 'p';
  ptrOutput++;
  *ptrOutput = '>';
  ptrOutput++;
  *ptrOutput = 0;    // Add terminator character.
}

void SubtractIntegerPolynomial(const int* minuend, const int* subtrahend, int* difference)
{
  int degreeMinuend = *minuend;
  int degreeSubtrahend = *subtrahend;
  const int* ptrMinuend = minuend + 1;
  const int* ptrSubtrahend = subtrahend + 1;
  int* ptrDifference = difference + 1;
  int minDegree = degreeMinuend;
  int currentDegree;
  int degreeDifference;
  if (minDegree > degreeSubtrahend)
  {
    minDegree = degreeSubtrahend;
  }
  for (currentDegree = 0; currentDegree <= minDegree; currentDegree++)
  {
    UncompressBigIntegerB(ptrMinuend, &operand1);
    UncompressBigIntegerB(ptrSubtrahend, &operand2);
    BigIntSubt(&operand1, &operand2, &operand3);
    NumberLength = operand3.nbrLimbs;
    BigInteger2IntArray(ptrDifference, &operand3);
    ptrMinuend += 1 + numLimbs(ptrMinuend);
    ptrSubtrahend += 1 + numLimbs(ptrSubtrahend);
    ptrDifference += 1 + numLimbs(ptrDifference);
  }
  if (degreeMinuend > minDegree)
  {
    (void)CopyPolynomial(ptrDifference, ptrMinuend, degreeMinuend - minDegree);
    *difference = degreeMinuend;
  }
  else if (degreeSubtrahend > minDegree)
  {
    for (; currentDegree <= degreeSubtrahend; currentDegree++)
    {
      *ptrDifference = -*ptrSubtrahend;
      memcpy(ptrDifference + 1, ptrSubtrahend + 1, numLimbs(ptrSubtrahend));
      ptrSubtrahend += 1 + numLimbs(ptrSubtrahend);
      ptrDifference += 1 + numLimbs(ptrDifference);
    }
  }
  else
  {  // Both minuend and subtrahend have the same degree.
     // Erase the most significant coefficients which are zero.
    degreeDifference = 0;
    ptrDifference = difference + 1;
    for (currentDegree = 0; currentDegree <= degreeMinuend; currentDegree++)
    {
      if ((*ptrDifference != 1) || (*(ptrDifference + 1) != 0))
      {
        degreeDifference = currentDegree;
      }
      ptrDifference += numLimbs(ptrDifference) + 1;
    }
    *difference = degreeDifference;
  }
}

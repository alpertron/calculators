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
#include <assert.h>
#include "bignbr.h"
#include "expression.h"
#include "highlevel.h"
#include "polynomial.h"
#include "rootseq.h"
#include "showtime.h"
#include "output.h"

#define MAX_MODULUS 100

#ifndef __ANDROID__
union uCommon common;
#endif

extern bool teachMod;
extern int grpLen;
extern char* ptrOutput; 
extern char* ptrOutput2;
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

static void showExponent(char **pptrOutput, int expon)
{
  if (expon == 1)
  {    // Do not show exponent if it equals 1.
    return;
  }
  if (pretty == PRETTY_PRINT)
  {
    copyStr(pptrOutput, "<sup>");
    int2dec(pptrOutput, expon);
    copyStr(pptrOutput, "</sup>");
  }
  else
  {
    char* ptrOut = *pptrOutput;
    *ptrOut++ = '^';
    if ((pretty == TEX) && (expon >= 10))
    {
      *ptrOut++ = '{';
      int2dec(&ptrOut, expon);
      *ptrOut++ = '}';
    }
    else
    {
      int2dec(&ptrOut, expon);
    }
    *pptrOutput = ptrOut;
  }
}

static void show1overM(int expon)
{
  if (pretty == PRETTY_PRINT)
  {
    copyStr(&ptrOutput, " = <f-f><f-n>1</f-n><f-d>p");
    showExponent(&ptrOutput, expon);
    copyStr(&ptrOutput, "</f-d></f-f> &#8290; (");
  }
  else
  {
    copyStr(&ptrOutput, " = (1/p");
    showExponent(&ptrOutput, expon);
    copyStr(&ptrOutput, ") * (");
  }
}

static void showSumOfProducts(void)
{
  if (pretty == PRETTY_PRINT)
  {
    copyStr(&ptrOutput,
      "<f-f><f-n>f<sub>1</sub>(x)</f-n><f-d>f<sub>1, 1</sub>(x)</f-d></f-f> &#8290; a<sub>1, 1</sub>(x) + ... + "
      "<f-f><f-n>f<sub>1</sub>(x)</f-n><f-d>f<sub>1, n</sub>(x)</f-d></f-f> &#8290; a<sub>1, n</sub>(x)");
  }
  else
  {
    copyStr(&ptrOutput, "(f_1(x)/f_{1, 1}(x)) * a_{1, 1}(x) + ... + (f_1(x)/f_{1, n}(x)) * a_{1, n}(x)");
  }
}

// Divide the polynomial by its leading coefficient.
// On output, operand1 is the inverse of the leading coefficient.
void ConvertToMonic(int *poly, int polyDegree)
{
  int nbrLimbs = NumberLength + 1;
  int currentDegree;
  if (NumberLength == 1)
  {         // Modulus size is one limb.
    int lenLimbs = (polyDegree * 2) + 1;
    int inverse = modInv(*(poly + lenLimbs), TestNbr[0].x);
    int* ptrPoly = poly + 1;
    intToBigInteger(&operand1, inverse);
    if (TestNbr[0].x <= 32768)
    {
      currentDegree = polyDegree;
      if ((currentDegree & 0x01) == 0x00)
      {
        *ptrPoly = *ptrPoly * inverse % TestNbr[0].x;
        ptrPoly += 2;
        currentDegree--;
      }
      for (; currentDegree >= 0; currentDegree -= 2)
      {
        *ptrPoly = *ptrPoly * inverse % TestNbr[0].x;
        *(ptrPoly + 2) = *(ptrPoly + 2) * inverse % TestNbr[0].x;
        ptrPoly += 4;
      }
    }
    else
    {
      for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
      {
#ifdef _USING64BITS_
        uint64_t u64Limb = (uint64_t)*ptrPoly * (uint64_t)inverse % (uint64_t)TestNbr[0].x;
        *ptrPoly = (int)u64Limb;
#else
        smallmodmult(*ptrPoly, inverse, (limb*)ptrPoly, TestNbr[0].x);
#endif
        ptrPoly += 2;
      }
    }
  }
  else
  {         // General case.
    int lenLimbs = polyDegree * nbrLimbs;
    int* ptrPoly = poly;
    IntArray2BigInteger(poly + lenLimbs, &operand1);
    (void)ModInvBigNbr(operand1.limbs, operand1.limbs, TestNbr, NumberLength);
    for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
    {
      IntArray2BigInteger(ptrPoly, &operand2);
      modmult(operand1.limbs, operand2.limbs, operand2.limbs);
      BigInteger2IntArray(ptrPoly, &operand2);
      ptrPoly += nbrLimbs;
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
    int lenLimbs = 1 + numLimbs(ptrPoly);
    ptrPoly += lenLimbs;
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
    int lenLimbs;
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
    lenLimbs = 1 + numLimbs(ptrPoly);
    ptrPoly += lenLimbs;
    lenLimbs = 1 + numLimbs(ptrPolyMod);
    ptrPolyMod += lenLimbs;
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
  enum eExprErr rc;
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
  common.poly.polyS[0] = 0;                  // Initialize polynomial g_m to zero.
  common.poly.polyS[1] = 1;
  common.poly.polyS[2] = 0;
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
    degree1 = getModPolynomial(common.poly.poly1, argF, &contentF);
    degree2 = getModPolynomial(common.poly.poly2, argG, &contentG);
    PolyModularGcd(common.poly.poly1, degree1, common.poly.poly2, degree2, common.poly.poly3, &degreeGcdMod);
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
      common.poly.polyS[0] = 0;                  // Initialize polynomial g_m to zero.
      common.poly.polyS[1] = 1;
      common.poly.polyS[2] = 0;
      potentialDegreeGcd = degreeGcdMod;
    }
    ptrSrc = &common.poly.poly3[0];
    if ((common.poly.polyS[0] == 0) && (common.poly.polyS[1] == 1) && (common.poly.polyS[2] == 0))
    {                                // If previous polynomial is zero, g_m <- (b * g_p) mod p.
      ptrDest = &common.poly.polyS[1];
      for (degree = 0; degree <= potentialDegreeGcd; degree++)
      {
        int lenLimbs;
        UncompressBigIntegerB(ptrSrc, &operand1);
        (void)BigIntMultiply(&operand1, &gcdLeadingCoeff, &operand1);
        intToBigInteger(&operand1, getRemainder(&operand1, prime));
        NumberLength = operand1.nbrLimbs;
        BigInteger2IntArray(ptrDest, &operand1);
        lenLimbs = 1 + numLimbs(ptrSrc);
        ptrSrc += lenLimbs;
        lenLimbs = 1 + numLimbs(ptrDest);
        ptrDest += lenLimbs;
      }
    }
    else
    {                                // Perform Chinese Remainder Theorem between g_p and g_m.
      int gcdLCModPrime;
      // coeff = c = a (mod p), c = b (mod q) => c = jq + b
      // a = jq + b (mod p) => j = (a-b)/q = r (mod p)
      // j = pk + r => c = (pk + r)q + b = pqk + rq + b => c = rq + b (mod pq).
      common.poly.polyT[0] = potentialDegreeGcd;
      ptrPrev = &common.poly.polyS[1];
      ptrDest = &common.poly.polyT[1];
      gcdLCModPrime = getRemainder(&gcdLeadingCoeff, prime);
      for (degree = 0; degree <= potentialDegreeGcd; degree++)
      {
        int lenLimbs;
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
        lenLimbs = 1 + numLimbs(ptrPrev);
        ptrPrev += lenLimbs;
        lenLimbs = 1 + numLimbs(ptrDest);
        ptrDest += lenLimbs;
      }
      (void)CopyPolynomial(&common.poly.polyS[1], &common.poly.polyT[1], potentialDegreeGcd);
    }
    common.poly.polyS[0] = potentialDegreeGcd;
    multint(&modulus, &modulus, prime);   // modulus <- modulus*prime.
    CopyBigInt(&operand4, &modulus);
    BigIntDivideBy2(&operand4);           // operand4 <- modulus / 2.
    // Get the primitive part of g_m in common.poly.polyT.
    common.poly.polyT[0] = potentialDegreeGcd;
    ptrSrc = &common.poly.polyS[1];
    ptrDest = &common.poly.polyT[1];
    for (degree = 0; degree <= potentialDegreeGcd; degree++)
    {
      int lenLimbs;
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
      lenLimbs = 1 + numLimbs(ptrSrc);
      ptrSrc += lenLimbs;
      lenLimbs = 1 + numLimbs(ptrDest);
      ptrDest += lenLimbs;
    }
    (void)CopyPolynomial(&common.poly.polyS[1], &common.poly.polyT[1], potentialDegreeGcd);
    (void)getContent(common.poly.polyS, &contentH);
    ptrSrc = &common.poly.polyS[1];
    ptrDest = &common.poly.polyT[1];
    for (degree = 0; degree <= potentialDegreeGcd; degree++)
    {
      int lenLimbs;
      UncompressBigIntegerB(ptrSrc, &operand1);
      (void)BigIntDivide(&operand1, &contentH, &operand2);
      NumberLength = operand2.nbrLimbs;
      BigInteger2IntArray(ptrDest, &operand2);
      lenLimbs = 1 + numLimbs(ptrSrc);
      ptrSrc += lenLimbs;
      lenLimbs = 1 + numLimbs(ptrDest);
      ptrDest += lenLimbs;
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
    common.poly.poly5[0] = argF[0];
    (void)CopyPolynomial(&common.poly.poly5[1], &argF[1], degree);
    rc = DivideIntegerPolynomial(common.poly.poly5, common.poly.polyT, TYPE_MODULUS);
    if ((rc == EXPR_OK) && (common.poly.poly5[0] == 0) && (common.poly.poly5[1] == 1) && (common.poly.poly5[2] == 0))
    {           // Remainder is zero.
      if (argG[0] < 0)
      {
        degree = 0;
      }
      else
      {
        degree = argG[0];
      }
      common.poly.poly5[0] = argG[0];
      (void)CopyPolynomial(&common.poly.poly5[1], &argG[1], degree);
      rc = DivideIntegerPolynomial(common.poly.poly5, common.poly.polyT, TYPE_MODULUS);
      if ((rc == EXPR_OK) && (common.poly.poly5[0] == 0) && (common.poly.poly5[1] == 1) && (common.poly.poly5[2] == 0))
      {           // Remainder is zero.
                  // Return h*gcd(contentF, contentG)
        BigIntGcd(&contentF, &contentG, &operand1);
        *gcd = potentialDegreeGcd;
        ptrSrc = &common.poly.polyT[1];
        ptrDest = gcd + 1;
        for (degree = 0; degree <= potentialDegreeGcd; degree++)
        {
          int lenLimbs;
          UncompressBigIntegerB(ptrSrc, &operand2);
          (void)BigIntMultiply(&operand1, &operand2, &operand3);
          NumberLength = operand3.nbrLimbs;
          BigInteger2IntArray(ptrDest, &operand3);
          lenLimbs = 1 + numLimbs(ptrSrc);
          ptrSrc += lenLimbs;
          lenLimbs = 1 + numLimbs(ptrDest);
          ptrDest += lenLimbs;
        }
        return;
      }
    }
  }
}

static void GcdModSizeOneLimb(int currentDegree, int degreeMin,
  const unsigned char * reduceModulus, const int* ptrArgMin, int* ptrArgMax)
{
  int mod = TestNbr[0].x;
  int index = (currentDegree * 2) + 1;
  int* ptrTemp = ptrArgMax + (currentDegree - degreeMin) * 2;
  int value1 = mod - *(ptrArgMax + index);
  const int* ptrPolynomial = ptrArgMin + 1;
  ptrTemp++;
  if ((NumberLength == 1) && (mod < 100))
  {
    for (index = degreeMin; index >= 3; index -= 4)
    {
      *ptrTemp = (int)reduceModulus[*ptrTemp + (*ptrPolynomial * value1)];
      *(ptrTemp + 2) = (int)reduceModulus[*(ptrTemp + 2) +
        (*(ptrPolynomial + 2) * value1)];
      *(ptrTemp + 4) = (int)reduceModulus[*(ptrTemp + 4) +
        (*(ptrPolynomial + 4) * value1)];
      *(ptrTemp + 6) = (int)reduceModulus[*(ptrTemp + 6) +
        (*(ptrPolynomial + 6) * value1)];
      ptrTemp += 8;
      ptrPolynomial += 8;
    }
    for (; index >= 1; index -= 2)
    {
      *ptrTemp = reduceModulus[*ptrTemp + (*ptrPolynomial * value1)];
      *(ptrTemp + 2) = reduceModulus[*(ptrTemp + 2) + (*(ptrPolynomial + 2) * value1)];
      ptrTemp += 4;
      ptrPolynomial += 4;
    }
    if (index == 0)
    {
      *ptrTemp = reduceModulus[*ptrTemp + (*ptrPolynomial * value1)];
    }
  }
  else if (mod < 32767)
  {
    for (index = 0; index <= degreeMin; index++)
    {
      *ptrTemp = (*ptrTemp + (*ptrPolynomial * value1)) % mod;
      ptrTemp += 2;
      ptrPolynomial += 2;
    }
  }
  else
  {
    for (index = 0; index <= degreeMin; index++)
    {
#ifdef _USING64BITS_
      uint64_t u64Limb = ((uint64_t)*ptrTemp + (uint64_t)*ptrPolynomial * (uint64_t)value1) %
        (uint64_t)mod;
      *ptrTemp = (int)u64Limb;
#else
      int temp;
      smallmodmult(*ptrPolynomial, value1, (limb*)&temp, mod);
      *ptrTemp = UintToInt(((unsigned int)*ptrTemp + (unsigned int)temp) %
        (unsigned int)mod);
#endif
      ptrTemp += 2;
      ptrPolynomial += 2;
    }
  }
}

static void ConvertToMonicSmallMod(int degreeMin, int* ptrArgMin, int modulus,
  const unsigned char* reduceModulus)
{
  int lenLimbs = (degreeMin * 2) + 1;
  int inverse = modInv(*(ptrArgMin + lenLimbs), modulus);
  if (inverse != 1)
  {
    int index;
    int* ptrPoly = ptrArgMin + 1;
    intToBigInteger(&operand1, inverse);
    for (index = degreeMin; index >= 3; index -= 4)
    {
      *ptrPoly = (int)reduceModulus[*ptrPoly * inverse];
      *(ptrPoly + 2) = (int)reduceModulus[*(ptrPoly + 2) * inverse];
      *(ptrPoly + 4) = (int)reduceModulus[*(ptrPoly + 4) * inverse];
      *(ptrPoly + 6) = (int)reduceModulus[*(ptrPoly + 6) * inverse];
      ptrPoly += 8;
    }
    for (; index >= 1; index -= 2)
    {
      *ptrPoly = (int)reduceModulus[*ptrPoly * inverse];
      *(ptrPoly + 2) = (int)reduceModulus[*(ptrPoly + 2) * inverse];
      ptrPoly += 4;
    }
    if (index == 0)
    {
      *ptrPoly = (int)reduceModulus[*ptrPoly * inverse];
    }
  }
}

// arg1 must not be changed by this routine.
// arg2 can be changed by this routine.
void PolyModularGcd(const int *arg1, int degree1, int *arg2, int degree2,
  int *gcd, int *degreeGcd)
{
  int nbrLimbs = NumberLength + 1;
  int *ptrArgMax = gcd;
  int *ptrArgMin = arg2;
  int *ptrTemp;
  int degreeMax = degree1;
  int degreeMin = degree2;
  int temp;
  int index;
  int lenBytes;
  unsigned char reduceModulus[MAX_MODULUS * MAX_MODULUS];
  int modulus = TestNbr[0].x;
  if (NumberLength > 1)
  {
    modulus = MAX_MODULUS;
  }
  if (degree1 == 0)
  {
    if ((*arg1 == 1) && (*(arg1 + 1) == 0))
    {     // First polynomial is zero. GCD is second argument.
      *degreeGcd = degree2;
      lenBytes = (degree2 + 1) * nbrLimbs * (int)sizeof(int);
      (void)memcpy(gcd, arg2, lenBytes);
      return;
    }
    *degreeGcd = 0;
    SetNumberToOne(common.poly.poly5);
    return;
  }
  if (degree2 == 0)
  {
    if ((*arg2 == 1) && (*(arg2 + 1) == 0))
    {     // Second polynomial is zero. GCD is first argument.
      *degreeGcd = degree1;
      lenBytes = (degree1 + 1) * nbrLimbs * (int)sizeof(int);
      (void)memcpy(gcd, arg1, lenBytes);
      return;
    }
    *degreeGcd = 0;
    SetNumberToOne(common.poly.poly5);
    return;
  }
  if ((degree1 < 0) && (degree2 < 0))
  {  // Both arguments are monomials. gcd is x^{min(|degree1|, |degree2|)}
    *degreeGcd = degree1;
    if (degree1 < degree2)
    {
      *degreeGcd = degree2;
    }
    *gcd = 1;         // Set coefficient to 1.
    *(gcd + 1) = 1;
    return;
  }
  if ((degree1 < 0) || (degree2 < 0))
  {  // One of the arguments is a monomial.
    int monomialDegree;
    int polynomialDegree;
    int currCoeff;
    const int* polynomialCoeffs;
    if (degree1 < 0)
    {
      monomialDegree = -degree1;
      polynomialDegree = degree2;
      polynomialCoeffs = arg2;
    }
    else
    {
      monomialDegree = -degree2;
      polynomialDegree = degree1;
      polynomialCoeffs = arg1;
    }
    for (currCoeff = 0;
      (currCoeff < polynomialDegree) && (currCoeff < monomialDegree);
      currCoeff++)
    {
      int coeffLength;
      if ((*polynomialCoeffs != 1) || (*(polynomialCoeffs+1) != 0))
      {
        break;    // Coefficient is not zero. Exit loop.
      }
      coeffLength = 1 + *polynomialCoeffs;
      polynomialCoeffs += coeffLength;
    }
    *degreeGcd = -currCoeff;
    *gcd = 1;         // Set coefficient to 1.
    *(gcd + 1) = 1;
    return;
  }
  if (modulus < MAX_MODULUS)
  {  // Initialize array reduceModulus when modulus is less than 100.
    unsigned char* ptrReduce = reduceModulus;
    for (int row = 0; row < modulus; row++)
    {
      for (int col = 0; col < modulus; col++)
      {
        *ptrReduce = (unsigned char)col;
        ptrReduce++;
      }
    }
  }
  // gcd of two polynomials.
  lenBytes = (degree1 + 1) * nbrLimbs * (int)sizeof(int);
  (void)memcpy(gcd, arg1, lenBytes);
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
    if (modulus < 100)
    {           // Modulus is less than 100.
      ConvertToMonicSmallMod(degreeMin, ptrArgMin, modulus, reduceModulus);
    }
    else
    {
      ConvertToMonic(ptrArgMin, degreeMin);
    }
    for (int currentDegree = degreeMax; currentDegree >= degreeMin; currentDegree--)
    {          // Replace ptrArgMax by remainder of long
               // division of ptrArgMax / ptrArgMin.
      index = (currentDegree - degreeMin) * nbrLimbs;
      if (nbrLimbs == 2)
      {        // Modulus size is one limb.
        GcdModSizeOneLimb(currentDegree, degreeMin, reduceModulus,
          ptrArgMin, ptrArgMax);
      }
      else
      {        // General case.
        int lenLimbs = currentDegree * nbrLimbs;
        ptrTemp = ptrArgMax + index;
        IntArray2BigInteger(ptrArgMax + lenLimbs, &operand1);
        for (index = 0; index <= degreeMin; index++)
        {
          lenLimbs = index * nbrLimbs;
          IntArray2BigInteger(ptrArgMin + lenLimbs, &operand2);
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
      int lenLimbs = degreeMax * nbrLimbs;
      ptrTemp = ptrArgMax + lenLimbs;
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
        lenBytes = (degreeMin + 1) * nbrLimbs * (int)sizeof(int);
        (void)memcpy(gcd, ptrArgMin, lenBytes);
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
  if (degreePoly == 0)
  {      // Degree of polynomial is zero: result polynomial must be zero.
    *ptrArgument = 0;
    *(ptrArgument+1) = 1;
    *(ptrArgument+2) = 0;
    return;
  }
  ptrSrc = ptrArgument + 1;
  ptrDest = &common.poly.poly1[1];
  for (degree = 1; degree <= degreePoly; degree++)
  {
    int lenLimbs = 1 + numLimbs(ptrSrc);
    ptrSrc += lenLimbs;
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
    lenLimbs = 1 + numLimbs(ptrDest);
    ptrDest += lenLimbs;
  }
  (void)CopyPolynomial(ptrArgument + 1, &common.poly.poly1[1], derivDegreePoly);
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
    int lenBytes = (NumberLength - 1) * (int)sizeof(limb);
    (void)memset(&operand2.limbs[1], 0, lenBytes);
  }
  for (int currentNbr = 0; currentNbr < qtyNbrs; currentNbr++)
  {
    int lenLimbs;
    IntArray2BigInteger(ptrNbr, &operand1);
    modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    BigInteger2IntArray(ptrNbr, &operand1);
    lenLimbs = NumberLength + 1;
    ptrNbr += lenLimbs;
  }
}

void polyToMontgomeryNotation(int *nbr, int qtyNbrs)
{
  int* ptrNbr = nbr;
  for (int currentNbr = qtyNbrs - 1; currentNbr >= 0; currentNbr--)
  {
    int lenLimbs;
    IntArray2BigInteger(ptrNbr, &operand1);
    modmult(operand1.limbs, MontgomeryMultR2, operand1.limbs);
    BigInteger2IntArray(ptrNbr, &operand1);
    lenLimbs = NumberLength + 1;
    ptrNbr += lenLimbs;
  }
}

  // Perform polyPower <- polyBase ^ expon (mod polyMod)
void powerPolynomial(int *polyBase, const int *polyMod, int polyDegree, const BigInteger *expon,
                     int *polyPower, powerCback callback, int curMultip, int nbrMultip)
{
  int lenLimbs;
  int nbrLimbs = NumberLength + 1;
  bool powerIsOne = true;
  int nbrBits = 0;
  int index = expon->nbrLimbs - 1;
  int bitCounter = 0;
  lenLimbs = polyDegree * nbrLimbs;
  *(polyBase + lenLimbs) = 1;
  *(polyBase + lenLimbs + 1) = 0;
  for (; index >= 0; index--)
  {
    int groupExp = expon->limbs[index].x;
    for (unsigned int mask = HALF_INT_RANGE_U; mask > 0U; mask >>= 1)
    {
      if (!powerIsOne)
      {
        if (callback != NULL)
        {
          callback(100 * bitCounter/(nbrBits*nbrMultip));
        }
        multUsingInvPolynomial(polyPower, polyPower, polyPower, polyDegree, polyMod);
      }
      if (((unsigned int)groupExp & mask) != 0U)
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
    int lenLimbs = (currentDegree * (NumberLength + 1));
    const int *ptrTemp = poly + lenLimbs;
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
static void ExtendedGcdPolynomial(/*@in@*/int *pointrA, int degA, /*@in@*/int * pointrB,
  int degB, /*@out@*/int *ptrQ,
  /*@out@*/int *ptrP1, /*@out@*/int *ptrP2, /*@out@*/int *ptrU, /*@out@*/int *pDegreeU)
{
  int* pointerA = pointrA;
  int* pointerB = pointrB;
  int degreeA = degA;
  int degreeB = degB;
  int *ptrQuotient = ptrQ;
  int *ptrQuotients[MAX_DEGREE+1];
  int* tmpPtr;
  int* ptrR;
  int *ptrOldR;
  int degreeR;
  int degreeOldR;
  int nbrQuotients = 0;
  int nbrLimbs = NumberLength + 1;
  int currentDegree;
  int degreeU;
  int lenBytes;
  bool polyExchanged = false;

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
  lenBytes = (degreeA + 1) * nbrLimbs * (int)sizeof(int);
  (void)memcpy(ptrOldR, pointerA, lenBytes);
  lenBytes = (degreeB + 1) * nbrLimbs * (int)sizeof(int);
  (void)memcpy(ptrR, pointerB, lenBytes);
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
  (void)ModInvBigNbr(operand5.limbs, (limb*)(ptrP2 + 1), powerMod.limbs, powerMod.nbrLimbs);
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
    // Multiply ptrR by ptrQuotients[counter]. The result will be stored in common.poly.polyMultTemp.
    MultPolynomial(degreeR, degreeQ, ptrR, ptrQuotients[counter]);
    // Compute ptrOldR <- ptrOldR - common.poly.polyMultTemp.
    degreeQ += degreeR; // Degree of product (always greater than degree of old R).
    offset = 0;
    for (currentDegree = 0; currentDegree <= degreeOldR; currentDegree++)
    {
      IntArray2BigInteger(ptrOldR + offset, &operand1);
      IntArray2BigInteger(&common.poly.polyMultTemp[offset], &operand2);
      SubtBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
      BigInteger2IntArray(ptrOldR + offset, &operand1);
      offset += nbrLimbs;
    }
    for (; currentDegree <= degreeQ; currentDegree++)
    {
      // Set operand1 to zero.
      lenBytes = NumberLength * (int)sizeof(limb);
      (void)memset(operand1.limbs, 0, lenBytes);
      IntArray2BigInteger(&common.poly.polyMultTemp[offset], &operand2);
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
  lenBytes = (degreeU + 1) * nbrLimbs * (int)sizeof(int);
  (void)memcpy(ptrU, ptrR, lenBytes);
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
// Compute 1 = (f(x)/f_1(x)) * a_1(x) + ... + (f(x)/f_n(x)) * a_n(x)
// This is done by computing n extended GCDs between f_i(x) and
// f(x)/f_i(x).
//
// Hensel lift from mod m to mod m^2:
// u(x) <- (1/m)*(f - f_1*f_2*...*f_n)  (mod m).
// For i = 1 to n compute:
//   g(x) <- u(x) * a_i(x) mod f_i(x) (all computations done mod m).
//   f_i(x) <- f_i(x) + m*g(x)
//
// Compute v(x) <- (1/m)*(f(x)/f_1(x)*a_1(x) + ... + f(x)/f_n(x)*a_n(x) - 1).
// For i = 1 to n compute:
//   h(x) <- v(x) * a_i(x) mod f_i(x) (all computations done mod m).
//   a_i(x) <- a_i(x) - m*h(x)
//
// Memory usage:
// common.poly.poly3: u, v
// common.poly.poly4: f, g, h
// common.poly.poly5: m*g, f/f_i*a_i
// common.poly.polyT: a_i
// common.poly.polyLifted: f_i.
// Change the fields ptrPolyLifted to point to f_i.

// Compute polynomial common.poly.poly4 <- common.poly.values divided by the leading coefficient
// mod prime. This is the polynomial f. Polynomial common.poly.values is in Montgomery notation.
static void ComputeF(bool isMontg)
{
  int currentDegree;
  const int* ptrValue1;
  int polyDegree = common.poly.values[0];                      // Get degree of polynomial
  // Find the leading coefficient.
  ptrValue1 = &common.poly.values[1];                          // Point to constant coefficient.
  for (currentDegree = 0; currentDegree < polyDegree; currentDegree++)
  {
    int lenLimbs = 1 + numLimbs(ptrValue1);
    ptrValue1 += lenLimbs;
  }
  NumberLength = numLimbs(ptrValue1);
  IntArray2BigInteger(ptrValue1, &operand1);
  (void)BigIntRemainder(&operand1, &powerMod, &operand1);
  NumberLength = powerMod.nbrLimbs;
  if (operand1.nbrLimbs < NumberLength)
  {
    int lenBytes = (NumberLength - operand1.nbrLimbs) * (int)sizeof(limb);
    (void)memset(&operand1.limbs[operand1.nbrLimbs], 0, lenBytes);
  }
  ptrValue1 = &common.poly.values[1];                          // Point to constant coefficient.
  // Compute the inverse of leading coefficient.
  (void)ModInvBigNbr(operand1.limbs, operand2.limbs, TestNbr, NumberLength);
  // Convert operand1 from standard to Montgomery notation.
  for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {
    int lenLimbs;
    NumberLength = numLimbs(ptrValue1);
    IntArray2BigInteger(ptrValue1, &operand1);
    (void)BigIntRemainder(&operand1, &powerMod, &operand1);
    NumberLength = powerMod.nbrLimbs;
    if (operand1.nbrLimbs < NumberLength)
    {
      int lenBytes = (NumberLength - operand1.nbrLimbs) * (int)sizeof(limb);
      (void)memset(&operand1.limbs[operand1.nbrLimbs], 0, lenBytes);
    }
    // Multiply by inverse. Result is in Montgomery notation.
    modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    if (!isMontg)
    {     // Convert to normal notation by multiplying by 1.
      int lenBytes = NumberLength * (int)sizeof(limb);
      memset(operand3.limbs, 0, lenBytes);
      operand3.limbs[0].x = 1;
      modmult(operand1.limbs, operand3.limbs, operand1.limbs);
    }
    BigInteger2IntArray(&common.poly.poly4[currentDegree * (NumberLength + 1)], &operand1);
    lenLimbs = 1 + numLimbs(ptrValue1);
    ptrValue1 += lenLimbs;               // Point to next coefficient.
  }
}

// Get polynomial a_i. This polynomial is stored with no spaces
// between coefficients. The output is in Montgomery notation in
// polynomial common.poly.poly1.
static int getAi(int nbrFactor, int degreeA)
{
  int nbrLimbs = NumberLength+1;
  const int* ptrSrc = ptrA[nbrFactor];
  int* ptrDest = common.poly.poly1;
  for (int currentDegree = 0; currentDegree <= degreeA; currentDegree++)
  {
    int numLen = *ptrSrc + 1;
    int lenBytes = numLen * (int)sizeof(int);
    (void)memcpy(ptrDest, ptrSrc, lenBytes);
    ptrDest += nbrLimbs;
    ptrSrc += numLen;
  }
  polyToMontgomeryNotation(common.poly.poly1, degreeA+1);
  return degreeA;
}

// Get polynomial factor f_i. This factor is stored with no spaces
// between coefficients. The output is in Montgomery notation in
// polynomial common.poly.poly2.
static void getFi(int degreeFactor, const int *src, int nbrLimbsDest, int nbrLimbsSrc)
{
  const int* ptrSrc = src;
  int NumberLengthBak = NumberLength;
  int *ptrDest = common.poly.poly2;       // Copy f_i to common.poly.poly2.
  for (int currentDegree = 0; currentDegree < degreeFactor; currentDegree++)
  {
    int nbrLen = *ptrSrc + 1;
    assert((nbrLen >= 1) && (nbrLen < MAX_LEN));
    int lenBytes = nbrLen * (int)sizeof(int);
    (void)memcpy(ptrDest, ptrSrc, lenBytes);
    if (nbrLimbsSrc == 0)
    {
      ptrSrc += nbrLen;
    }
    else
    {
      ptrSrc += nbrLimbsSrc;
    }
    ptrDest += nbrLimbsDest;
  }
  NumberLength = nbrLimbsDest - 1;
  SetNumberToOne(ptrDest);
  polyToMontgomeryNotation(common.poly.poly2, degreeFactor);
  NumberLength = NumberLengthBak;
}

// Move factors f_i from common.poly.polyLiftedNew to common.poly.polyLifted, then adjust
// factorInfo[...].ptrPolyLifted to point to the new factors.
// There are no spaces between coefficients.
static void MoveFactorsAndFixPointers(struct sFactorInfo* ptrFactorInfo, bool compressPoly)
{
  struct sFactorInfo* pstFactorInfo = ptrFactorInfo;
  const int* ptrSrc = common.poly.polyLiftedNew;
  int * ptrDest = common.poly.polyLifted;
  for (int nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    int degreeFactor = pstFactorInfo->degree;
    pstFactorInfo->ptrPolyLifted = ptrDest;
    for (int currentDegree = 0; currentDegree < degreeFactor; currentDegree++)
    {
      int numLength = *ptrSrc + 1;
      int lenBytes = numLength * (int)sizeof(int);
      (void)memcpy(ptrDest, ptrSrc, lenBytes);
      ptrSrc += numLength;
      if (compressPoly)
      {
        ptrDest += numLength;
      }
      else
      {
        int lenLimbs = NumberLength + 1;
        ptrDest += lenLimbs;
      }
    }
    pstFactorInfo++;
  }
}

static void showOneSubindex(char poly, int index)
{
  *ptrOutput++ = poly;
  if (pretty == PRETTY_PRINT)
  {
    copyStr(&ptrOutput, "<sub>");
    int2dec(&ptrOutput, index);
    copyStr(&ptrOutput, "</sub>");
  }
  else
  {
    copyStr(&ptrOutput, "_");
    int2dec(&ptrOutput, index);
  }
}

static void showTwoSubindexes(char poly, int firstIndex, int secondIndex)
{
  *ptrOutput++ = poly;
  if (pretty == PRETTY_PRINT)
  {
    copyStr(&ptrOutput, "<sub>");
    int2dec(&ptrOutput, firstIndex);
    copyStr(&ptrOutput, ", ");
    int2dec(&ptrOutput, secondIndex);
    copyStr(&ptrOutput, "</sub>");
  }
  else
  {
    copyStr(&ptrOutput, "_{");
    int2dec(&ptrOutput, firstIndex);
    copyStr(&ptrOutput, ", ");
    int2dec(&ptrOutput, secondIndex);
    copyStr(&ptrOutput, "}");
  }
}

// Perform Hensel lifting of polynomials in factorInfo.
// On entry:
// common.poly.values = original polynomial. The product of all lifted factors
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
  int* ptrPolyLifted = common.poly.polyLifted;
  int* ptrPoly2;
  int** ptrNewFactor;
  struct sFactorInfo *pstFactorInfo;
  // Copy polynomials f_i(x) to common.poly.polyLifted.
  int* ptrDest = ptrPolyLifted;
  pstFactorInfo = ptrFactorInfo;
  NumberLength = primeMod.nbrLimbs;
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
  ptrPoly2 = common.poly.poly2;
  pstFactorInfo = ptrFactorInfo;
  for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    if (pstFactorInfo->multiplicity != 1)
    {   // Repeated factor found. Indicate error.
      return EXPR_CANNOT_LIFT;
    }
    // Init common.poly.poly2 coefficients to zero.
    for (int index = 0; index < pstFactorInfo->degree; index++)
    {
      *ptrPoly2 = 1;
      ptrPoly2++;
      *ptrPoly2 = 0;
      ptrPoly2 += NumberLength;
    }
    degree += pstFactorInfo->degree;
    pstFactorInfo++;
  }
  if (teachMod)
  {
    showText(lang ? "<h3>Lema de Hensel</h3>" :
      "<h3>Hensel Lifting</h3>");
  }
  computePower(1);
  // Compute a_i for each factor f_i by computing the extended gcd between
  // f_i(x) and f(x)/f_i(x).
  ptrDest = common.poly.polyT;
  pstFactorInfo = ptrFactorInfo;
  for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    int lenBytes;
    ComputeF(true);                   // common.poly.poly4 <- f mod prime (Montgomery notation).
    nbrLimbs = NumberLength + 1;
    degreeFactor = pstFactorInfo->degree;
    lenBytes = pstFactorInfo->degree * nbrLimbs * (int)sizeof(int);
    (void)memcpy(common.poly.poly2, pstFactorInfo->ptr, lenBytes);
    polyToMontgomeryNotation(common.poly.poly2, degreeFactor + 1);
    // The factors do not include their leading coefficient. Set it to one.
    SetNumberToOne(&common.poly.poly2[pstFactorInfo->degree * nbrLimbs]);
    DividePolynomial(common.poly.poly4, degree,   // Dividend f(x).
      common.poly.poly2, degreeFactor,            // Divisor f_i(x).
      common.poly.poly1);                         // Quotient f(x)/f_i(x).
    ExtendedGcdPolynomial(common.poly.poly1, degree - degreeFactor, common.poly.poly2, degreeFactor,
      common.poly.poly3, common.poly.poly4, common.poly.poly5, common.poly.polyS, &degreeS);
    polyToStandardNotation(common.poly.polyS, degreeS + 1);
    // At this moment: common.poly.polyS = a_i on Montgomery notation, degreeS = degree of a_i.
    ptrA[nbrFactor] = ptrDest;
    ptrDest = CopyPolynomial(ptrDest, common.poly.polyS, degreeS);
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
    if ((elapsedTime / 10) != (oldTimeElapsed / 10))
    {
      char outputInfo[1000];
      char* ptrOut = outputInfo;
      if (lang)
      {
        copyStr(&ptrOut, "1<p>Aplicando lema de Hensel usando el número primo ");
        int2dec(&ptrOut, primeMod.limbs[0].x);
        copyStr(&ptrOut, " procesando exponente ");
        int2dec(&ptrOut, currentExp);
        copyStr(&ptrOut, " de ");
      }
      else
      {
        copyStr(&ptrOut, "1<p>Hensel lifting using prime number ");
        int2dec(&ptrOut, primeMod.limbs[0].x);
        copyStr(&ptrOut, " processing exponent ");
        int2dec(&ptrOut, currentExp);
        copyStr(&ptrOut, " of ");
      }
      int2dec(&ptrOut, exponentMod);
      copyStr(&ptrOut, ".</p>");
      showElapsedTimeSec(&ptrOut);
      databack(outputInfo);
    }
#endif
    // Compute u(x) as (1/m) * (f - f_1 * f_2 *...* f_n)
    oldNumberLength = NumberLength;
    newExponent = currentExp * 2;    // We can double the exponent in each step.
    if (newExponent > exponentMod)   // Do not exceed exponentMod.
    {
      newExponent = exponentMod;
    }
    (void)BigIntPowerIntExp(&primeMod, currentExp, &operand5);
    computePower(newExponent);       // Compute powerMod and init Montgomery parms.
    newNumberLength = NumberLength;
    nbrLimbs = NumberLength + 1;
    if (teachMod)
    {
      ComputeF(false);               // common.poly.poly4 <- f mod m^newExponent (normal notation).
      copyStr(&ptrOutput, "<p>");
      showOneSubindex('f', currentExp);
      copyStr(&ptrOutput, " = x");
      showExponent(&ptrOutput, common.poly.values[0]);
      *ptrOutput++ = ' ';
      showPolynomial(&ptrOutput, common.poly.poly4, degree, grpLen);
      copyStr(&ptrOutput, "</p>");
      if (currentExp == 1)
      {   // If exponent is not 1, f_{m, n} is already shown.
        computePower(currentExp);
        pstFactorInfo = ptrFactorInfo;
        for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
        {    // Loop that shows all factors.
          degreeFactor = pstFactorInfo->degree;
          copyStr(&ptrOutput, "<p>");
          showTwoSubindexes('f', 1, nbrFactor + 1);
          copyStr(&ptrOutput, " = x");
          showExponent(&ptrOutput, degreeFactor);
          *ptrOutput++ = ' ';
          computePower(newExponent);
          if (currentExp == 1)
          {
            getFi(degreeFactor, pstFactorInfo->ptrPolyLifted, nbrLimbs, oldNumberLength + 1);
          }
          else
          {     // Copy f_i to common.poly.poly2.
            getFi(degreeFactor, pstFactorInfo->ptrPolyLifted, nbrLimbs, 0);
          }
          showMontPolynomial(&ptrOutput, common.poly.poly2, degreeFactor, grpLen);
          copyStr(&ptrOutput, "</p>");
          pstFactorInfo++;
        }
        copyStr(&ptrOutput, lang? "<p>Usando el algoritmo de MCD extendido, calcular ":
          "<p>Using the extended GCD algorithm, compute ");
        if (pretty == PRETTY_PRINT)
        {
          copyStr(&ptrOutput, "a<sub>1, 1</sub>, ..., a<sub>1, n</sub>");
        }
        else
        {
          copyStr(&ptrOutput, "a_{1, 1}, ..., a_{1, n}");
        }
        copyStr(&ptrOutput, lang ? " tales que</p><p>" : " such that</p><p>");
        copyStr(&ptrOutput, "1 = ");
        showSumOfProducts();
        if (pretty != PRETTY_PRINT)
        {
          *ptrOutput++ = ')';
        }
        copyStr(&ptrOutput, "</p>");
        pstFactorInfo = ptrFactorInfo;
        for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
        {    // Loop that shows all a_i.
          degreeFactor = pstFactorInfo->degree;
          computePower(currentExp);
          degreeA = getAi(nbrFactor, degreeFactor - 1);
          copyStr(&ptrOutput, "<p>");
          showTwoSubindexes('a', 1, nbrFactor + 1);
          copyStr(&ptrOutput, " = ");
          showMontPolynomial(&ptrOutput, common.poly.poly1, degreeA, grpLen);
          copyStr(&ptrOutput, "</p>");
          pstFactorInfo++;
        }
      }
      computePower(newExponent);
    }
    ComputeF(true);                 // common.poly.poly4 <- f mod m^newExponent (Montgomery notation).
    SetNumberToOne(common.poly.poly1);          // Initialize product of factors.
    pstFactorInfo = ptrFactorInfo;
    degree = 0;
    for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
    {    // Loop that multiplies all factors. common.poly.poly1 holds the result.
      int lenBytes;
      degreeFactor = pstFactorInfo->degree;
      if (currentExp == 1)
      {
        getFi(degreeFactor, pstFactorInfo->ptrPolyLifted, nbrLimbs, oldNumberLength + 1);
      }
      else
      {     // Copy f_i to common.poly.poly2.
        getFi(degreeFactor, pstFactorInfo->ptrPolyLifted, nbrLimbs, 0);
      }
      MultPolynomial(degree, degreeFactor, common.poly.poly1, common.poly.poly2);
      degree += degreeFactor;
      lenBytes = (degree + 1) * nbrLimbs * (int)sizeof(int);
      (void)memcpy(common.poly.poly1, common.poly.polyMultTemp, lenBytes);              // common.poly.poly1 <- product.
      pstFactorInfo++;
    }
    polyToStandardNotation(common.poly.poly1, degree + 1);   // f_1 * f_2 *...* f_n
    polyToStandardNotation(common.poly.poly4, degree + 1);   // f
    for (currentDegree = 0; currentDegree <= degree; currentDegree++)
    {                 // Loop that computes (1/m)*(f - f_1 * f_2 * ... * f_n)
      int nbrLen;
      int lenBytes;
      // Get coefficient of f.
      IntArray2BigInteger(&common.poly.poly4[currentDegree * nbrLimbs], &operand1);
      // Get coefficient of product of factors.
      IntArray2BigInteger(&common.poly.poly1[currentDegree * nbrLimbs], &operand2);
      SubtBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
      // At this moment: operand1.limbs = coeff(f - f_1 * f_2 * ... * f_n)
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
      ptrDest = &common.poly.poly3[currentDegree * (oldNumberLength + 1)];
      *ptrDest = operand1.nbrLimbs;
      ptrDest++;
      lenBytes = operand1.nbrLimbs * (int)sizeof(int);
      (void)memcpy(ptrDest, operand1.limbs, lenBytes);
    }
    computePower(currentExp);
    if (teachMod)
    {
      copyStr(&ptrOutput, "<p>");
      showOneSubindex('u', currentExp);
      show1overM(currentExp);
      showOneSubindex('f', currentExp);
      *ptrOutput++ = ' ';
      copyStr(&ptrOutput, ptrMinus);
      *ptrOutput++ = ' ';
      showTwoSubindexes('f', currentExp, 1);
      for (nbrFactor = 2; nbrFactor <= nbrFactorsFound; nbrFactor++)
      {
        *ptrOutput++ = ' ';
        copyStr(&ptrOutput, ptrTimes);
        *ptrOutput++ = ' ';
        showTwoSubindexes('f', currentExp, nbrFactor);
      }
      copyStr(&ptrOutput, ") = ");
      showPolynomial(&ptrOutput, common.poly.poly3, degree, grpLen);
      copyStr(&ptrOutput, "</p>");
    }
    polyToMontgomeryNotation(common.poly.poly3, degree + 1);
    nbrLimbs = NumberLength + 1;
    pstFactorInfo = ptrFactorInfo;
    ptrDest = common.poly.polyLiftedNew;
    ptrNewFactor = ptrNewFactors;
    computePower(currentExp);
    for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
    {
      int lenBytes;
      degreeFactor = pstFactorInfo->degree;
      // g(x) <- u(x) * a_i(x) mod f_i(x) (all computations done mod p^currentExp).
      degreeA = getAi(nbrFactor, degreeFactor - 1); // common.poly.poly1 <- a_i.
      MultPolynomial(degree, degreeA, common.poly.poly3, common.poly.poly1);
      lenBytes = (degree + degreeA + 1) * nbrLimbs * (int)sizeof(int);
      (void)memcpy(common.poly.poly1, common.poly.polyMultTemp, lenBytes);
                                                 // common.poly.poly1 <- u * a_i
                                                 // Copy f_i to common.poly.poly2.
      if (currentExp == 1)
      {
        getFi(degreeFactor, pstFactorInfo->ptrPolyLifted, nbrLimbs, oldNumberLength + 1);
      }
      else
      {     // Copy f_i to common.poly.poly2.
        getFi(degreeFactor, pstFactorInfo->ptrPolyLifted, nbrLimbs, 0);
      }
      DividePolynomial(common.poly.poly1, degree + degreeA,  // Dividend = u * a_i
        common.poly.poly2, degreeFactor,                     // Divisor = f_i
        NULL);                                   // Quotient: not needed.
                                                 // Remainder will be stored in common.poly.poly1.
      polyToStandardNotation(common.poly.poly1, degreeFactor);  // Convert g to standard notation.
      if (teachMod)
      {
        copyStr(&ptrOutput, "<p>");
        showTwoSubindexes('g', currentExp, nbrFactor + 1);
        copyStr(&ptrOutput, " = ");
        showOneSubindex('u', currentExp);
        *ptrOutput++ = ' ';
        copyStr(&ptrOutput, ptrTimes);
        *ptrOutput++ = ' ';
        showTwoSubindexes('a', currentExp, nbrFactor + 1);
        copyStr(&ptrOutput, (pretty == PRETTY_PRINT? " mod ": " % "));
        showTwoSubindexes('f', currentExp, nbrFactor + 1);
        copyStr(&ptrOutput, " = ");
        showPolynomial(&ptrOutput, common.poly.poly1, degreeFactor, grpLen);
        copyStr(&ptrOutput, "</p>");
      }
      // Compute new value of f_i(x) as f_i(x) + m*g(x)
      const int *ptrFi = pstFactorInfo->ptrPolyLifted;
      const int *ptrStartPolyF = ptrDest;
      *ptrNewFactor = ptrDest;                   // Point to f_i(x) mod new modulus.
      ptrNewFactor++;
      for (currentDegree = 0; currentDegree < degreeFactor; currentDegree++)
      {                 // Loop that computes f_i + m*g
        // Get coefficient of g.
        ptrSrc = &common.poly.poly1[currentDegree * nbrLimbs];
        operand1.nbrLimbs = *ptrSrc;
        ptrSrc++;
        lenBytes = operand1.nbrLimbs * (int)sizeof(int);
        (void)memcpy(operand1.limbs, ptrSrc, lenBytes);
        if ((currentExp * 2) > newExponent)
        {
          (void)BigIntPowerIntExp(&primeMod, newExponent - currentExp, &operand5);
          (void)BigIntRemainder(&operand1, &operand5, &operand1);
        }
        // Get coefficient of f_i.
        operand2.nbrLimbs = *ptrFi;
        ptrFi++;
        lenBytes = operand2.nbrLimbs * (int)sizeof(int);
        (void)memcpy(operand2.limbs, ptrFi, lenBytes);
        if (currentExp == 1)
        {
          ptrFi += NumberLength;
        }
        else
        {
          ptrFi += operand2.nbrLimbs;
        }
        (void)BigIntMultiply(&operand1, &powerMod, &operand1);  // common.poly.poly1 <- m*g
        BigIntAdd(&operand1, &operand2, &operand1);       // common.poly.poly1 <- f_i + m*g
        // Store coefficient of new f_i (no spaces between coefficients).
        *ptrDest = operand1.nbrLimbs;
        ptrDest++;
        lenBytes = operand1.nbrLimbs * (int)sizeof(int);
        (void)memcpy(ptrDest, &operand1.limbs, lenBytes);
        ptrDest += operand1.nbrLimbs;    // Point to next coeffficient of f_i.
      }
      if (teachMod)
      {
        computePower(newExponent);         // Compute powerMod and init Montgomery parms.
        copyStr(&ptrOutput, "<p>");
        showTwoSubindexes('f', newExponent, nbrFactor + 1);
        copyStr(&ptrOutput, " = ");
        showTwoSubindexes('f', currentExp, nbrFactor + 1);
        copyStr(&ptrOutput, " + ");
        *ptrOutput++ = 'p';
        showExponent(&ptrOutput, currentExp);
        *ptrOutput++ = ' ';
        copyStr(&ptrOutput, ptrTimes);
        *ptrOutput++ = ' ';
        showTwoSubindexes('g', currentExp, nbrFactor + 1);
        copyStr(&ptrOutput, " = x");
        computePower(newExponent);
        showExponent(&ptrOutput, degreeFactor);
        showPolynomial(&ptrOutput, ptrStartPolyF, degreeFactor, grpLen);
        copyStr(&ptrOutput, "</p>");
        computePower(currentExp);
      }
      pstFactorInfo++;
    }
    computePower(newExponent);         // Compute powerMod and init Montgomery parms.
    if (newExponent == exponentMod)
    {                 // Final common.poly.values of f_1, f_2,..., f_n found. Exit loop.
      break;
    }
    // Compute v(x) <- (1/m)*(f(x)/f_1(x)*a_1(x) + ... + f(x)/f_n(x)*a_n(x) - 1).
    // Init v(x) to zero.
    ptrDest = common.poly.poly3;
    nbrLimbs = newNumberLength + 1;
    for (currentDegree = 0; currentDegree <= degree; currentDegree++)
    {
      *ptrDest = 1;
      *(ptrDest + 1) = 0;
      ptrDest += nbrLimbs;
    }
    pstFactorInfo = ptrFactorInfo;
    ptrNewFactor = ptrNewFactors;
    ComputeF(true);              // common.poly.poly4 <- f (mod m^2) in Montgomery notation.
    for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
    {
      int lenBytes;
      degreeFactor = pstFactorInfo->degree;
      lenBytes = (degree + 1) * nbrLimbs * (int)sizeof(int);
      (void)memcpy(common.poly.poly1, common.poly.poly4, lenBytes);
                                        // common.poly.poly2 <- f_i(x) in Montgomery notation.
      getFi(degreeFactor, *ptrNewFactor, nbrLimbs, 0);
      ptrNewFactor++;
      DividePolynomial(common.poly.poly1, degree,   // Dividend = f.
        common.poly.poly2, degreeFactor,            // Divisor = f_i.
        common.poly.poly5);                         // Quotient = f / f_i.
      degreeA = getAi(nbrFactor, degreeFactor - 1);       // common.poly.poly1 <- a_i.
                                        // common.poly.polyMultTemp <- (f / f_i) * a_i
      MultPolynomial(degreeA, degree - degreeFactor, common.poly.poly1, common.poly.poly5);
      ptrDest = common.poly.poly3;
      for (currentDegree = 0; currentDegree <= (degree - degreeFactor + degreeA); currentDegree++)
      {
        // Get coefficient of sum.
        IntArray2BigInteger(ptrDest, &operand1);
        // Get coefficient of (f / f_i) * a_i.
        IntArray2BigInteger(&common.poly.polyMultTemp[currentDegree * nbrLimbs], &operand2);
        // Add them.
        AddBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
        // Save coefficient of sum.
        BigInteger2IntArray(ptrDest, &operand1);
        ptrDest += nbrLimbs;
      }
      pstFactorInfo++;
    }
    polyToStandardNotation(common.poly.poly3, degree + 1);
    // Divide sum by m.
    computePower(currentExp);
    ptrDest = common.poly.poly3;
    ptrSrc = common.poly.poly3;
    for (currentDegree = 0; currentDegree <= degree; currentDegree++)
    { // Get coefficient of u.
      int lenBytes;
      operand1.nbrLimbs = *ptrSrc;
      lenBytes = operand1.nbrLimbs * (int)sizeof(int);
      (void)memcpy(operand1.limbs, ptrSrc + 1, lenBytes);
      // Divide coefficient by m.
      (void)BigIntDivide(&operand1, &powerMod, &operand1);
      // Store coefficient of quotient.
      *ptrDest = operand1.nbrLimbs;
      ptrDest++;
      lenBytes = operand1.nbrLimbs * (int)sizeof(int);
      (void)memcpy(ptrDest, operand1.limbs, lenBytes);
      ptrDest += oldNumberLength;
      ptrSrc += nbrLimbs;
    }
    polyToMontgomeryNotation(common.poly.poly3, degree);   // common.poly.poly3 <- v(x)
    if (teachMod)
    {
      copyStr(&ptrOutput, "<p>");
      showOneSubindex('v', currentExp);
      show1overM(currentExp);
      showSumOfProducts();
      *ptrOutput++ = ' ';
      copyStr(&ptrOutput, ptrMinus);
      copyStr(&ptrOutput, " 1) = ");
      showMontPolynomial(&ptrOutput, common.poly.poly3, degree, grpLen);
      copyStr(&ptrOutput, "</p>");
    }
    nbrLimbs = NumberLength + 1;
    pstFactorInfo = ptrFactorInfo;
    // Use common.poly.polyS as a temporary storage for a_i.
    ptrDest = common.poly.polyS;
    for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
    {
      int lenBytes;
      degreeFactor = pstFactorInfo->degree;
      // h(x) <- v(x) * a_i(x) mod f_i(x) (all computations done mod m).
      NumberLength = oldNumberLength;
      degreeA = getAi(nbrFactor, degreeFactor - 1);  // common.poly.poly1 <- a_i
      MultPolynomial(degree, degreeA, common.poly.poly3, common.poly.poly1);
      lenBytes = (degree + degreeA + 1) * nbrLimbs * (int)sizeof(int);
      (void)memcpy(common.poly.poly4, common.poly.polyMultTemp, lenBytes);
                                                // common.poly.poly4 <- u * a_i
                                                // common.poly.poly2 <- f_i(x) in Montgomery notation. 
      getFi(degreeFactor, pstFactorInfo->ptrPolyLifted, nbrLimbs, 0);
      DividePolynomial(common.poly.poly4, degree + degreeA, // Dividend = u * a_i
        common.poly.poly2, degreeFactor,                    // Divisor = f_i
        NULL);                                  // Quotient: not needed.
      polyToStandardNotation(common.poly.poly1, degreeA+1); // Convert a_i to standard notation.
                                                // Convert g to standard notation.
      polyToStandardNotation(common.poly.poly4, degree + degreeA - degreeFactor + 1);
      // Subtract m*g(x) from a_i(x).
      const int* ptrStartPolyG = ptrDest;
      for (currentDegree = 0; currentDegree <= degreeA; currentDegree++)
      {                 // Loop that computes a_i - m*g
        // Get coefficient of g.
        ptrSrc = &common.poly.poly4[currentDegree * nbrLimbs];
        operand1.nbrLimbs = *ptrSrc;
        ptrSrc++;
        lenBytes = operand1.nbrLimbs * (int)sizeof(int);
        (void)memcpy(operand1.limbs, ptrSrc, lenBytes);
        if (!BigIntIsZero(&operand1))                     // g <- -g (mod m)
        {
          BigIntSubt(&powerMod, &operand1, &operand1);
        }
        // Get coefficient of a_i.
        ptrSrc = &common.poly.poly1[currentDegree * nbrLimbs];
        operand2.nbrLimbs = *ptrSrc;
        ptrSrc++;
        lenBytes = operand2.nbrLimbs * (int)sizeof(int);
        (void)memcpy(operand2.limbs, ptrSrc, lenBytes);
        (void)BigIntMultiply(&operand1, &powerMod, &operand1);  // common.poly.poly1 <- m*g
        BigIntAdd(&operand1, &operand2, &operand1);       // common.poly.poly1 <- a_i + m*g
        // Store coefficient of new a_i.
        *ptrDest = operand1.nbrLimbs;
        ptrDest++;
        lenBytes = operand1.nbrLimbs * (int)sizeof(int);
        (void)memcpy(ptrDest, operand1.limbs, lenBytes);
        ptrDest += operand1.nbrLimbs;  // Point to next coefficient of a_i.
      }
      if (teachMod)
      {
        copyStr(&ptrOutput, "<p>");
        showTwoSubindexes('h', currentExp, nbrFactor + 1);
        copyStr(&ptrOutput, " = ");
        showOneSubindex('v', currentExp);
        *ptrOutput++ = ' ';
        copyStr(&ptrOutput, ptrTimes);
        *ptrOutput++ = ' ';
        showTwoSubindexes('a', currentExp, nbrFactor + 1);
        copyStr(&ptrOutput, (pretty == PRETTY_PRINT ? " mod " : " % "));
        showTwoSubindexes('f', currentExp, nbrFactor + 1);
        copyStr(&ptrOutput, " = ");
        showMontPolynomial(&ptrOutput, common.poly.poly3, degreeFactor, grpLen);
        copyStr(&ptrOutput, "</p>");
        computePower(newExponent);         // Compute powerMod and init Montgomery parms.
        copyStr(&ptrOutput, "<p>");
        showTwoSubindexes('a', newExponent, nbrFactor + 1);
        copyStr(&ptrOutput, " = ");
        showTwoSubindexes('a', currentExp, nbrFactor + 1);
        copyStr(&ptrOutput, " + ");
        *ptrOutput++ = 'p';
        showExponent(&ptrOutput, currentExp);
        *ptrOutput++ = ' ';
        copyStr(&ptrOutput, ptrTimes);
        *ptrOutput++ = ' ';
        showTwoSubindexes('g', currentExp, nbrFactor + 1);
        copyStr(&ptrOutput, " = ");
        computePower(newExponent);
        showPolynomial(&ptrOutput, ptrStartPolyG, degreeFactor, grpLen);
        copyStr(&ptrOutput, "</p>");
        computePower(currentExp);
      }
      pstFactorInfo++;
    }
    currentExp = newExponent;
    NumberLength = newNumberLength;
    // Move a_i from common.poly.polyS to common.poly.polyT.
    ptrSrc = common.poly.polyS;
    ptrDest = common.poly.polyT;
    pstFactorInfo = ptrFactorInfo;
    for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
    {
      degreeA = pstFactorInfo -> degree - 1;
      ptrA[nbrFactor] = ptrDest;
      for (currentDegree = 0; currentDegree <= degreeA; currentDegree++)
      {
        int lenBytes;
        nbrLimbs = *ptrSrc;
        ptrSrc++;
        *ptrDest = nbrLimbs;
        ptrDest++;
        lenBytes = nbrLimbs * (int)sizeof(int);
        (void)memcpy(ptrDest, ptrSrc, lenBytes);
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
  int lenBytes = powerMod.nbrLimbs * (int)sizeof(limb);
  (void)memcpy(&TestNbr, powerMod.limbs, lenBytes);
  NumberLength = powerMod.nbrLimbs;
  TestNbr[NumberLength].x = 0;
  GetMontgomeryParms(powerMod.nbrLimbs);
  degree = common.poly.values[0];
  ptrValue1 = &common.poly.values[1];
  ptrValue2 = &common.poly.poly4[0];
  for (int currentDegree = 0; currentDegree <= degree; currentDegree++)
  {
    int lenLimbs;
    IntArray2BigInteger(ptrValue1, &operand1);
    // Convert operand1 from Montgomery to standard notation.
    operand2.limbs[0].x = 1;
    if (NumberLength > 1)
    {
      lenBytes = (NumberLength - 1) * (int)sizeof(limb);
      (void)memset(&operand2.limbs[1], 0, lenBytes);
    }
    modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    BigInteger2IntArray(ptrValue2, &operand1);
    lenLimbs = 1 + *ptrValue1;
    ptrValue1 += lenLimbs;
    lenLimbs = 1 + *ptrValue2;
    ptrValue2 += lenLimbs;
  }
  (void)memcpy(&common.poly.values[1], &common.poly.poly4[0], (ptrValue2 - &common.poly.poly4[0])*sizeof(int));
}

void showPower(char **pptrOutput, int exponent)
{
  char *ptrOut = *pptrOutput;
  if (pretty == PRETTY_PRINT)
  {
    *ptrOut = '<';
    ptrOut++;
    *ptrOut = 's';
    ptrOut++;
    *ptrOut = 'u';
    ptrOut++;
    *ptrOut = 'p';
    ptrOut++;
    *ptrOut = '>';
    ptrOut++;
    int2dec(&ptrOut, exponent);
    *ptrOut = '<';
    ptrOut++;
    *ptrOut = '/';
    ptrOut++;
    *ptrOut = 's';
    ptrOut++;
    *ptrOut = 'u';
    ptrOut++;
    *ptrOut = 'p';
    ptrOut++;
    *ptrOut = '>';
    ptrOut++;
  }
  else if (pretty == TEX)
  {
    *ptrOut = '^';
    ptrOut++;
    *ptrOut = '{';
    ptrOut++;
    int2dec(&ptrOut, exponent);
    *ptrOut = '}';
    ptrOut++;
  }
  else
  {
    *ptrOut = '^';
    ptrOut++;
    int2dec(&ptrOut, exponent);
  }
  *pptrOutput = ptrOut;
}

void showPowerVar(char** pptrOutput, int polyDegree, char letter)
{
  char* ptrOut = *pptrOutput;
  if (polyDegree == 0)
  {
    *ptrOut = '1';
    ptrOut++;
  }
  else
  {
    if (pretty == PRETTY_PRINT)
    {
      if (letter == 'x')
      {
        copyStr(&ptrOut, lang ? "<var role=\"img\" aria-label=\"equis\">x</var>" : "<var>x</var>");
      }
      else
      {
        showVariable(&ptrOut, letter);
      }
    }
    else
    {
      *ptrOut = letter;
      ptrOut++;
    }
    if (polyDegree != 1)
    {
      showPower(&ptrOut, polyDegree);
    }
  }
  *pptrOutput = ptrOut;
}

void showPowerX(char **pptrOutput, int polyDegree)
{
  showPowerVar(pptrOutput, polyDegree, 'x');
}

static void getMontCoeff(const int* ptrCoeff)
{
  NumberLength = powerMod.nbrLimbs;
  IntArray2BigInteger(ptrCoeff, &operand2);
  memset(operand1.limbs, 0, NumberLength * sizeof(limb));
  operand1.limbs[0].x = 1;
  modmult(operand1.limbs, operand2.limbs, operand1.limbs);
  operand1.nbrLimbs = NumberLength;
  while ((operand1.nbrLimbs > 1) && (operand1.limbs[operand1.nbrLimbs - 1].x == 0))
  {   // Loop that discards non-significant limbs.
    operand1.nbrLimbs--;
  }
}

static void showPolynomialMontOrNorm(char **pptrOutput, const int *ptrPoly,
      int polyDegree, int groupLength, bool isMontNotation)
{
  int currentDegree;
  char *ptrOut = *pptrOutput;
  int *ptrIndex;
  int indexes[MAX_DEGREE+1];
  int NumberLengthBak = NumberLength;

  assert(polyDegree >= 0);
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
    int nbrLimbs;
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
      *ptrOut = ' ';
      ptrOut++;
      if (*ptrValue1 > 0)
      {
        *ptrOut = '+';
        ptrOut++;
      }
      else
      {
        copyStr(&ptrOut, (pretty == PRETTY_PRINT)? "&minus;": "-");
      }
      *ptrOut = ' ';
      ptrOut++;
      if (isMontNotation)
      {
        getMontCoeff(ptrValue1);
      }
      else
      {
        NumberLength = numLimbs(ptrValue1);
        IntArray2BigInteger(ptrValue1, &operand1);
        NumberLength = powerMod.nbrLimbs;
      }
      if ((operand1.nbrLimbs > 1) || (operand1.limbs[0].x != 1))
      {            // Absolute value of coefficient is not one.
        operand1.sign = SIGN_POSITIVE;
        if (pretty == PRETTY_PRINT)
        {          // Show number of digits if there are more than 30.
          Bin2Dec(&ptrOut, operand1.limbs, operand1.nbrLimbs, groupLength);
        }
        else
        {         // Do not show number of digits.
          Bin2Dec(&ptrOut, operand1.limbs, operand1.nbrLimbs, -groupLength);
        }
        if (currentDegree > 0)
        {
          copyStr(&ptrOut, ptrTimes);
          showPowerX(&ptrOut, currentDegree);
        }
      }
      else
      {
        showPowerX(&ptrOut, currentDegree);
      }
    }
  }
  *pptrOutput = ptrOut;
  NumberLength = NumberLengthBak;
}

void showPolynomial(char** pptrOutput, const int* ptrPoly, int polyDegree, int groupLength)
{
  showPolynomialMontOrNorm(pptrOutput, ptrPoly, polyDegree, groupLength, false);
}

void showMontPolynomial(char** pptrOutput, const int* ptrPoly, int polyDegree, int groupLength)
{
  showPolynomialMontOrNorm(pptrOutput, ptrPoly, polyDegree, groupLength, true);
}

static void outputOriginalPolynomialElem(char** pptrOutput, const int* ptrPoly, int groupLength)
{
  int currentDegree;
  const int* ptrValue1;
  int nbrLimbs = powerMod.nbrLimbs + 1;
  char* ptrOut = *pptrOutput;
  degree = *ptrPoly;
  ptrValue1 = ptrPoly + 1;
  if (!modulusIsZero)
  {
    int* ptrValue2;
    if (pretty == PARI_GP)
    {
      copyStr(&ptrOut, "Mod(");
    }
    // Output polynomial to factor. First move polynomial to common.poly.poly4
    ptrValue2 = &common.poly.poly4[0];
    for (currentDegree = 0; currentDegree <= degree; currentDegree++)
    {
      int lenLimbs = 1 + *ptrValue1;
      int lenBytes = lenLimbs * (int)sizeof(int);
      (void)memcpy(ptrValue2, ptrValue1, lenBytes);
      ptrValue1 += lenLimbs;
      ptrValue2 += nbrLimbs;
    }
    // Get leading coefficient.
    getMontCoeff(&common.poly.poly4[degree * nbrLimbs]);
  }
  else
  { // Find leading coefficient.
    for (currentDegree = 0; currentDegree < degree; currentDegree++)
    {
      int lenLimbs = 1 + numLimbs(ptrValue1);
      ptrValue1 += lenLimbs;
    }
    NumberLength = numLimbs(ptrValue1);
    IntArray2BigInteger(ptrValue1, &operand1);
  }
  if (operand1.sign == SIGN_NEGATIVE)
  {
    copyStr(&ptrOut, " &minus;");
  }
  if ((operand1.nbrLimbs != 1) || (operand1.limbs[0].x != 1) || (degree == 0))
  {     // Leading coefficient is not 1 or degree is zero.
    if (pretty == PRETTY_PRINT)
    {          // Show number of digits if there are more than 30.
      Bin2Dec(&ptrOut, operand1.limbs, operand1.nbrLimbs, groupLength);
    }
    else
    {         // Do not show number of digits.
      Bin2Dec(&ptrOut, operand1.limbs, operand1.nbrLimbs, -groupLength);
    }
    if ((pretty == PARI_GP) && (degree > 0))
    {
      *ptrOut = '*';
      ptrOut++;
    }
  }
  if (degree < 0)
  {
    showPowerX(&ptrOut, -degree);
  }
  else
  {
    if (degree > 0)
    {
      showPowerX(&ptrOut, degree);
    }
    if (modulusIsZero)
    {
      showPolynomial(&ptrOut, ptrPoly + 1, degree, groupLength);
    }
    else
    {
      showMontPolynomial(&ptrOut, common.poly.poly4, degree, groupLength);
    }
  }
  if (!modulusIsZero)
  {
    if (pretty == PRETTY_PRINT)
    {
      copyStr(&ptrOut, " (mod ");
    }
    else if (pretty == TEX)
    {
      copyStr(&ptrOut, " (\\pmod ");
    }
    else
    {     // Pari-GP
      copyStr(&ptrOut, ", ");
    }
    if (pretty == PRETTY_PRINT)
    {          // Show number of digits if there are more than 30.
      Bin2Dec(&ptrOut, primeMod.limbs, primeMod.nbrLimbs, groupLength);
    }
    else
    {         // Do not show number of digits.
      Bin2Dec(&ptrOut, primeMod.limbs, primeMod.nbrLimbs, -groupLength);
    }
    if (exponentMod != 1)
    {
      showPower(&ptrOut, exponentMod);
    }
    *ptrOut = ')';
    ptrOut++;
  }
  *ptrOut = 0;    // Append string terminator.
  *pptrOutput = ptrOut;
}

void outputOriginalPolynomial(char** pptrOutput, int groupLength)
{
  bool isFraction = false;
  if ((modulusIsZero && ((common.poly.denom[0] != 0) || (common.poly.denom[1] != 1) || (common.poly.denom[2] != 1))) ||
    (!modulusIsZero && (common.poly.denom[0] != 0)))
  {    // If modulus is zero: denominator is not one.
       // If modulus is not zero: degree is not zero.
    isFraction = true;
    if (pretty == PRETTY_PRINT)
    {
      copyStr(pptrOutput, "<f-f><f-n>");
    }
    else if (pretty == PARI_GP)
    {
      copyStr(pptrOutput, "(");
    }
    else
    {
      copyStr(pptrOutput, "\\frac{");
    }
  }
  outputOriginalPolynomialElem(pptrOutput, common.poly.values, groupLength);
  if (!isFraction)
  {
    return;
  }
  if (pretty == PRETTY_PRINT)
  {
    copyStr(pptrOutput, "</f-n><f-d>");
  }
  else if (pretty == PARI_GP)
  {
    copyStr(pptrOutput, ") / (");
  }
  else
  {
    copyStr(pptrOutput, "}{");
  }
  outputOriginalPolynomialElem(pptrOutput, common.poly.denom, groupLength);
  if (pretty == PRETTY_PRINT)
  {
    copyStr(pptrOutput, "</f-d></f-f>");
  }
  else if (pretty == PARI_GP)
  {
    copyStr(pptrOutput, ")");
  }
  else
  {
    copyStr(pptrOutput, "}");
  }
}

void outputPolynomialFactor(char **pptrOutput, int groupLength, const struct sFactorInfo* pstFactorInfo)
{
  int polyDegree = pstFactorInfo->degree;
  int multiplicity = pstFactorInfo->multiplicity;
  int isMonomial = ((polyDegree == 1) && (*pstFactorInfo->ptrPolyLifted == 1) &&
     (*(pstFactorInfo->ptrPolyLifted+1) == 0));
  char* ptrOut = *pptrOutput;
  if ((multiplicity > 1) && !isMonomial)
  {
    *ptrOut = '(';
    ptrOut++;
  }
  if (modulusIsZero)
  {
    // Get leading coefficient.
    const int *ptrSrc = pstFactorInfo->ptrPolyLifted;
    for (int currentDegree = 0; currentDegree < polyDegree; currentDegree++)
    {
      int lenLimbs = 1 + numLimbs(ptrSrc);
      ptrSrc += lenLimbs;
    }
    UncompressBigIntegerB(ptrSrc, &operand1);
    if (operand1.sign == SIGN_NEGATIVE)
    {
      copyStr(&ptrOut, "&minus;");
    }
    if ((operand1.nbrLimbs != 1) || (operand1.limbs[0].x != 1))
    {     // Absolute value is not 1.
      Bin2Dec(&ptrOut, operand1.limbs, operand1.nbrLimbs, groupLength);
    }
  }
  showPowerX(&ptrOut, polyDegree);
  showPolynomial(&ptrOut, pstFactorInfo->ptrPolyLifted, polyDegree, groupLength);
  if (multiplicity > 1)
  {
    if (!isMonomial)
    {
      *ptrOut = ')';
      ptrOut++;
    }
    showPower(&ptrOut, multiplicity);
  }
  *ptrOut = 0;    // Append string terminator.
  *pptrOutput = ptrOut;
}

void textErrorPol(char **pptrOutput, enum eExprErr rc)
{
  char *ptrOut = *pptrOutput;
  *ptrOut = '<';
  ptrOut++;
  *ptrOut = 'p';
  ptrOut++;
  *ptrOut = '>';
  ptrOut++;
  switch (rc)
  {
  case EXPR_CANNOT_USE_X_IN_EXPONENT:
    copyStr(&ptrOut, lang ? "No se puede usar variable en el exponente" :
      "Cannot use variable in exponent");
    break;
  case EXPR_POLYNOMIAL_DIVISION_NOT_INTEGER:
    copyStr(&ptrOut, lang ? "La división de polinomios no es entera" :
      "Polynomial division is not integer");
    break;
  case EXPR_DENOMINATOR_MUST_BE_CONSTANT:
    copyStr(&ptrOut, lang ? "El denominador debe ser constante" :
      "Denominator must be constant");
    break;
  case EXPR_DEGREE_TOO_HIGH:
    copyStr(&ptrOut, lang ? "El grado del polinomio es muy elevado" :
      "Degree is too high");
    break;
  case EXPR_EXPONENT_TOO_LARGE:
    copyStr(&ptrOut, lang ? "Exponente muy grande" : "Exponent is too large");
    break;
  case EXPR_EXPONENT_NEGATIVE:
    copyStr(&ptrOut, lang ? "Exponente negativo" : "Exponent is negative");
    break;
  case EXPR_LEADING_COFF_MULTIPLE_OF_PRIME:
    copyStr(&ptrOut, lang ? "El primer coeficiente es múltiplo del número primo" :
      "Leading coefficient multiple of prime");
    break;
  case EXPR_CANNOT_LIFT:
    copyStr(&ptrOut, lang ? "No se puede elevar porque hay factores duplicados" :
      "Cannot lift because of duplicate factors modulo prime");
    break;
  case EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE:
    copyStr(&ptrOut, lang ? "El módulo debe ser mayor que 1" : "Modulus must be greater than one");
    break;
  case EXPR_MODULUS_MUST_BE_PRIME_EXP:
    copyStr(&ptrOut, lang ? "El módulo debe ser un número primo o una potencia de número primo" :
      "Modulus must be a prime number or a power of a prime");
    break;
  case EXPR_MULTIPLE_VARIABLES_NOT_ACCEPTED:
    copyStr(&ptrOut, lang ? "No se aceptan múltiples variables" :
      "Multiple variables are not accepted");
    break;
  case EXPR_MORE_THAN_ONE_EQUAL_SIGN:
    copyStr(&ptrOut, lang ? "No se acepta más de un signo igual" :
      "More than one equal sign is not accepted");
    break;
  case EXPR_EQUAL_SIGN_INSIDE_PAREN:
    copyStr(&ptrOut, lang ? "Hay un signo igual dentro de paréntesis" :
      "Equal sign inside parentheses");
    break;
  default:
    textError(&ptrOut, rc);
    break;
  }
  *ptrOut = '<';
  ptrOut++;
  *ptrOut = '/';
  ptrOut++;
  *ptrOut = 'p';
  ptrOut++;
  *ptrOut = '>';
  ptrOut++;
  *ptrOut = 0;    // Add terminator character.
  *pptrOutput = ptrOut;
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
    int lenLimbs;
    UncompressBigIntegerB(ptrMinuend, &operand1);
    UncompressBigIntegerB(ptrSubtrahend, &operand2);
    BigIntSubt(&operand1, &operand2, &operand3);
    NumberLength = operand3.nbrLimbs;
    BigInteger2IntArray(ptrDifference, &operand3);
    lenLimbs = 1 + numLimbs(ptrMinuend);
    ptrMinuend += lenLimbs;
    lenLimbs = 1 + numLimbs(ptrSubtrahend);
    ptrSubtrahend += lenLimbs;
    lenLimbs = 1 + numLimbs(ptrDifference);
    ptrDifference += lenLimbs;
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
      int lenLimbs;
      *ptrDifference = -*ptrSubtrahend;
      (void)memcpy(ptrDifference + 1, ptrSubtrahend + 1, numLimbs(ptrSubtrahend));
      lenLimbs = 1 + numLimbs(ptrSubtrahend);
      ptrSubtrahend += lenLimbs;
      lenLimbs = 1 + numLimbs(ptrDifference);
      ptrDifference += lenLimbs;
    }
  }
  else
  {  // Both minuend and subtrahend have the same degree.
     // Erase the most significant coefficients which are zero.
    degreeDifference = 0;
    ptrDifference = difference + 1;
    for (currentDegree = 0; currentDegree <= degreeMinuend; currentDegree++)
    {
      int lenLimbs;
      if ((*ptrDifference != 1) || (*(ptrDifference + 1) != 0))
      {
        degreeDifference = currentDegree;
      }
      lenLimbs = numLimbs(ptrDifference) + 1;
      ptrDifference += lenLimbs;
    }
    *difference = degreeDifference;
  }
}

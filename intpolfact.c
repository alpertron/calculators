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

#define MAX_MATRIX_SIZE  15

BigInteger basis[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE];
BigInteger basisStar[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE];
BigInteger lambda[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE];
BigInteger prodB[MAX_MATRIX_SIZE];
BigInteger detProdB[MAX_MATRIX_SIZE];

static BigInteger contentPolyToFactor, halfPowerMod;
static int origPolyToFactor[1000000];
static int polyToFactor[1000000];
static int factorX[] = { 1, 0, 1, 1 }; // Polynomial is x.
static BigInteger bound, trailingCoeff, leadingCoeff;
static int polyNonRepeatedFactors[1000000];
static int tempPoly[1000000];
static int polyBackup[1000000];
static int polyLiftedRecord[1000000];
static struct sFactorInfo factorInfoRecord[MAX_DEGREE];
struct sFactorInfo factorInfoInteger[MAX_DEGREE];
static int arrNbrFactors[MAX_DEGREE];
extern int polyLifted[1000000];
extern int polyS[1000000];
extern int poly5[1000000];
int polyInteger[1000000];

// Compute the orthogonal basis of E (basisStar), lambda, prodB
// and the determinant of products of B.
static void GramSchmidtOrthogonalization(int size)
{
  BigInteger* pDest;
  int colI, colJ, row;
  for (colI = 0; colI < size; colI++)
  {
    // Compute b*_i = b_i - sum_j(u_{i, j} b*_j /d_j) (u = lambda)
    for (row = 0; row < size; row++)
    {
      CopyBigInt(&basisStar[row][colI], &basis[row][colI]);
    }
    for (colJ = 0; colJ < colI; colJ++)
    {
      for (row = 0; row < size; row++)
      {
        BigIntMultiply(&lambda[colI][colJ], &basisStar[colJ][row], &tmp0);
        BigIntDivide(&tmp0, &detProdB[colJ], &tmp1);
        BigIntSubt(&basisStar[row][colJ], &tmp1, &basisStar[row][colJ]);
      }
    }
    // Compute B_i = b*_i * b*_i
    intToBigInteger(&prodB[colI], 0);
    for (row = 0; row < size; row++)
    {
      BigIntMultiply(&basisStar[row][colI], &basisStar[row][colI], &tmp0);
      BigIntAdd(&prodB[colI], &tmp0, &prodB[colI]);
    }
    if (colI == 0)
    {
      CopyBigInt(&detProdB[0], &prodB[0]);
    }
    else
    {
      BigIntMultiply(&detProdB[colI-1], &prodB[colI], &detProdB[colI]);
    }
    // Compute lambda_{i, j} = b_i * b*_j * d_j / B_j = b_i * b*_j * d_{j-1}
    for (colJ = 0; colJ < colI; colJ++)
    {
      pDest = &lambda[colI][colJ];
      intToBigInteger(pDest, 0);
      for (row = 0; row < size; row++)
      {
        BigIntMultiply(&basis[row][colI], &basisStar[row][colJ], &tmp0);
        BigIntAdd(pDest, &tmp0, pDest);
      }
      if (colI > 0)
      {
        BigIntMultiply(pDest, &detProdB[colJ - 1], pDest);
      }
    }
  }
}

static void PerformREDI(int k, int l, int size)
{
  int i, row;
  // If |2 lambda_{k. l}| <= d_l, go out.
  BigInteger* pLambda = &lambda[k][l];
  BigInteger* pDet = &detProdB[l];
  enum eSign signDivision;
  multint(&tmp0, pLambda, 2);
  tmp0.sign = SIGN_POSITIVE;
  BigIntSubt(pDet, &tmp0, &tmp1);
  if (tmp1.sign == SIGN_NEGATIVE)
  {
    return;
  }
  // Compute q = nearest integer to lambda{k, l} / d_l
  if (pLambda->sign == pDet->sign)
  {
    signDivision = SIGN_POSITIVE;
  }
  else
  {
    signDivision = SIGN_NEGATIVE;
  }
  CopyBigInt(&tmp1, pDet);
  tmp1.sign = SIGN_POSITIVE;
  BigIntAdd(&tmp0, &tmp1, &tmp0);
  multint(&tmp1, &tmp1, 2);       // Compute 2 d_l.
  BigIntDivide(&tmp0, &tmp1, &tmp2);
  tmp2.sign = signDivision;       // tmp2 <- q.
  for (row = 0; row < size; row++)
  {  // b_k <- b_k - q*b_l
    BigIntMultiply(&tmp2, &basis[row][l], &tmp0);
    BigIntSubt(&basis[row][k], &tmp0, &basis[row][k]);
  }
    // lambda_{k, l} <- lambda_{k, l} - q*d_l
  BigIntMultiply(&tmp2, &detProdB[l], &tmp0);
  BigIntSubt(&lambda[k][l], &tmp0, &lambda[k][l]);
  for (i = 0; i < l; i++)
  { // lambda_{k, i} <- lambda_{k, i} - q*lambda_{l, i}
    BigIntMultiply(&tmp2, &lambda[l][i], &tmp0);
    BigIntSubt(&lambda[k][i], &tmp0, &lambda[k][i]);
  }
}

static void PerformSWAPI(int k, int kMax, int size)
{
  int row, i, j;
  // Exchange b_k with b_{k-1}
  for (row = 0; row < size; row++)
  {
    CopyBigInt(&tmp0, &basis[row][k]);
    CopyBigInt(&basis[row][k], &basis[row][k-1]);
    CopyBigInt(&basis[row][k-1], &tmp0);
  }
  for (j = 0; j < k - 1; j++)
  {  // Exchange lambda_{k, j} with lambda_{k-1, j}
    CopyBigInt(&tmp0, &lambda[k][j]);
    CopyBigInt(&lambda[k][j], &lambda[k-1][j]);
    CopyBigInt(&lambda[k-1][j], &tmp0);
  }
    // Set lambda <- lambda_{k, k-1}
  CopyBigInt(&tmp0, &lambda[k][k - 1]);    // tmp0 <- lambda.
    // Set B <- (d_{k-2}*d_k + lambda^2)/d_{k-1}
  if (k > 1)
  {
    BigIntMultiply(&detProdB[k - 2], &detProdB[k], &tmp1);
  }
  else
  {
    CopyBigInt(&tmp1, &detProdB[k]);
  }
  BigIntMultiply(&tmp0, &tmp0, &tmp2);
  BigIntAdd(&tmp1, &tmp2, &tmp1);
  BigIntDivide(&tmp1, &detProdB[k - 1], &tmp1); // tmp1 <- B.
  for (i = k; i <= kMax; i++)
  {
    // t <- lambda_{i, k}
    CopyBigInt(&tmp2, &lambda[i][k]);
    // lambda_{i, k} <- (d_k*lambda_{i, k-1} - lambda * t)/d_{k-1}
    BigIntMultiply(&detProdB[k], &lambda[i][k-1], &tmp3);
    BigIntMultiply(&tmp0, &tmp2, &tmp4);
    BigIntSubt(&tmp3, &tmp4, &tmp3);
    BigIntDivide(&tmp3, &detProdB[k - 1], &lambda[i][k]);
    // lambda_{i, k-1} <- (B*t + lambda * lambda_{i, k})/d_k
    BigIntMultiply(&tmp0, &tmp1, &tmp3);
    BigIntMultiply(&tmp2, &lambda[i][k], &tmp4);
    BigIntAdd(&tmp3, &tmp4, &tmp3);
    BigIntDivide(&tmp3, &detProdB[k], &lambda[i][k - 1]);
  }
  CopyBigInt(&detProdB[k - 1], &tmp1);  // d_{k-1} <- B.
}

void integralLLL(int size)
{
  int k, kMax, row;
  int colI, colJ;
  GramSchmidtOrthogonalization(size);
  k = 1;
  kMax = 0;
  do
  {
    // Perform incremental Gram-Schmidt
    if (k > kMax)
    {
      kMax = k;
      for (colJ = 0; colJ <= k; colJ++)
      {
        // Compute u = b_k * b_j (use tmp2 for u)
        intToBigInteger(&tmp2, 0);
        for (row = 0; row < size; row++)
        {
          BigIntMultiply(&basis[row][k], &basisStar[row][colJ], &tmp0);
          BigIntAdd(&tmp2, &tmp0, &tmp2);
        }
        for (colI = 0; colI < colJ; colI++)
        {    // Set u = (d_i * u - lambda_{k, i} * lambda_{j, i}) / d_{i-1}
          BigIntMultiply(&detProdB[colI], &tmp2, &tmp0);
          BigIntMultiply(&lambda[k][colI], &lambda[colJ][colI], &tmp1);
          if (colI == 0)
          {
            BigIntSubt(&tmp0, &tmp1, &tmp2);
          }
          else
          {
            BigIntSubt(&tmp0, &tmp1, &tmp0);
            BigIntDivide(&tmp0, &detProdB[colI - 1], &tmp2);
          }
        }
        if (colJ < k)
        {
          CopyBigInt(&lambda[k][colJ], &tmp2);
        }
        else
        {
          CopyBigInt(&detProdB[k], &tmp2);
        }
      }
    }
    // Test LLL condition
    for (;;)
    {
      int l;
      PerformREDI(k, k - 1, size);
      // Check whether 4 * d_k * d_{k-2} < 3 (d_{k-1})^2 - 4*(lambda_{k, k-1})^2
      if (k > 1)
      {
        BigIntMultiply(&detProdB[k], &detProdB[k - 2], &tmp0);
        multint(&tmp0, &tmp0, 4);
        BigIntMultiply(&detProdB[k - 1], &detProdB[k - 1], &tmp1);
        multint(&tmp1, &tmp1, 3);
        BigIntMultiply(&lambda[k][k-1], &lambda[k][k - 1], &tmp2);
        multint(&tmp2, &tmp2, 4);
        BigIntSubt(&tmp1, &tmp2, &tmp1);    // tmp1 = Right Hand Side.
        BigIntSubt(&tmp0, &tmp1, &tmp0);
        if (tmp0.sign == SIGN_POSITIVE)
        {
          break;
        }
        PerformSWAPI(k, kMax, size);
        if (--k < 1)
        {
          k = 1;
        }
      }
      for (l = k - 2; l >= 0; l--)
      {
        PerformREDI(k, l, size);
      }
      k++;
    }
  } while (k < size);
}

// Generate integer polynomial from modular polynomial.
static void GenerateIntegerPolynomial(int* polyMod, int* polyInt, int degreePoly)
{
  int currentDegree;
  *polyInt = degreePoly;
  int* ptrSrc = polyMod;
  int* ptrDest = polyInt + 1;
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
static void InsertIntegerPolynomialFactor(int* ptrFactor, int degreePoly)
{
  struct sFactorInfo* pstFactorInfoInteger, * pstFactorInfo;
  int indexNewFactor[MAX_DEGREE];
  int indexOldFactor[MAX_DEGREE];
  int currentDegree, index;
  int* ptrIndex, * ptrOldFactor;

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
      (pstFactorInfo - pstFactorInfoInteger) * sizeof(*pstFactorInfo));
  }
  pstFactorInfoInteger->multiplicity = 1;
  pstFactorInfoInteger->ptrPolyLifted = ptrFactor;
  pstFactorInfoInteger->degree = degreePoly;
}

// Input: values = degree, coefficient degree 0, coefficient degree 1, etc.
// Output: factorInfo = structure that holds the factors.
int FactorPolyOverIntegers(void)
{
  int degreePolyToFactor = values[0];
  int degree1, degree2, rc;
  int primeRecord = 0;
  int exponRecord = 0;
  int expon, maxDegreeFactor;
  int degreeGcdMod;
  int* ptrSrc, * ptrDest;
  int attemptNbr;
  int currentFactor;
  int currentDegree;
  int polXprocessed = FALSE;
  int nbrTmp[1000], nbrTmp2[1000], nbrTmp3[1000];
  int* ptrFactorIntegerBak;
  struct sFactorInfo* pstFactorInfoOrig, * pstFactorInfoRecord;
  struct sFactorInfo* pstFactorInfoInteger = factorInfoInteger;
  int* ptrFactorInteger = polyInteger;
  modulusIsZero = 1;
  memset(factorInfoInteger, 0, sizeof(factorInfoInteger));
  getContent(values, &contentPolyToFactor);
  CopyPolynomial(&origPolyToFactor[1], &values[1], degreePolyToFactor);
  origPolyToFactor[0] = degreePolyToFactor;
  // polyToFactor -> original polynomial / content of polynomial.
  // Let n the degree of the least coefficient different from zero.
  // Then x^n divides the polynomial.
  polyToFactor[0] = degreePolyToFactor;
  ptrSrc = &origPolyToFactor[1];
  ptrDest = &polyToFactor[1];
  for (currentDegree = 0; currentDegree <= degreePolyToFactor; currentDegree++)
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
        polyToFactor[0] = degreePolyToFactor - currentDegree;
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
    int halfDegree;
    int nbrFactorsRecord;
    int prime;
    int sumOfDegrees;
    modulusIsZero = 1;
    // Get trailing coefficient.
    ptrSrc = &values[1];
    UncompressBigIntegerB(ptrSrc, &operand1);
    // Get leading coefficient.
    for (degree1 = 0; degree1 < degreePolyToFactor; degree1++)
    {
      ptrSrc += 1 + numLimbs(ptrSrc);
    }
    UncompressBigIntegerB(ptrSrc, &operand2);
    BigIntMultiply(&operand1, &operand2, &trailingCoeff);
    // Compute F/gcd(F, F') where F is the polynomial to factor.
    // In the next loop we will attempt to factor gcd(F, F').
    degreePolyToFactor = polyToFactor[0];
    CopyPolynomial(&polyNonRepeatedFactors[1], &polyToFactor[1], degreePolyToFactor);
    CopyPolynomial(&tempPoly[1], &polyToFactor[1], degreePolyToFactor);
    tempPoly[0] = degreePolyToFactor;
    DerPolynomial(tempPoly);
    PolynomialGcd(tempPoly, polyToFactor, polyToFactor);
    polyNonRepeatedFactors[0] = degreePolyToFactor;
    DivideIntegerPolynomial(polyNonRepeatedFactors, polyToFactor, TYPE_DIVISION);
    prime = 3;
    // Find Knuth-Cohen bound for coefficients of polynomial factors:
    // If polynomial B divides A we have for all j:
    // |Bj| <= binomial(n-1, j)*SUM(i, |Ai|^2))^(1/2) + binomial(n-1, j-1) * |Am|
    // where m is the degree of A and n is the degree of B.
    // Maximum degree to be considered is n = ceil(m/2).
    // We need to find max(Bj).
    modulusIsZero = 1;
    degreePolyToFactor = polyNonRepeatedFactors[0];
    ptrSrc = &polyNonRepeatedFactors[1];
    UncompressBigIntegerB(ptrSrc, &operand1);
    if (degreePolyToFactor < 0)
    {      // Monomial.
      maxDegreeFactor = (-degreePolyToFactor + 1) / 2;
      operand1.sign = SIGN_POSITIVE;
      CopyBigInt(&operand3, &operand1);         // Get leading coefficient.
    }
    else
    {      // Polynomial.
      maxDegreeFactor = (degreePolyToFactor + 1) / 2;
      BigIntMultiply(&operand1, &operand1, &operand1);
      for (degree1 = 1; degree1 <= degreePolyToFactor; degree1++)
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
    bound.sign = SIGN_POSITIVE;
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
      int factorNbr, * ptrPolyLiftedRecord;
      int nbrFactors;
      // Get leading coefficient of polyNonRepeatedFactors.
      ptrSrc = &polyNonRepeatedFactors[1];
      for (degree1 = 0; degree1 < degreePolyToFactor; degree1++)
      {
        ptrSrc += 1 + numLimbs(ptrSrc);
      }
      UncompressBigIntegerB(ptrSrc, &leadingCoeff);
      do
      {   // Loop that finds a prime modulus such that the factorization
          // has no repeated factors. That means gcd(F, F') = 1 (mod prime).
          // This is required because Hensel lift does not work when
          // repeated factors are present.
          // If the leading coefficient is multiple of prime, the prime cannot
          // be used.
        do
        {      // Loop while the leading coefficient is multiple of prime.
          prime = nextPrime(prime);
        } while (getRemainder(&leadingCoeff, prime) == 0);
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
      degreePolyToFactor = getModPolynomial(&poly1[1], polyNonRepeatedFactors, &operand5);
      poly1[0] = degreePolyToFactor;
      polyBackup[0] = values[0];
      CopyPolynomial(&polyBackup[1], &values[1], values[0] >= 0 ? values[0] : 0);
      values[0] = degreePolyToFactor;
      CopyPolynomial(&values[1], &poly1[1], degreePolyToFactor);
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
    for (degree1 = 0; degree1 < degreePolyToFactor; degree1++)
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
      int* ptrTrailingCoeff;

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
          &polyLifted[(currentFactor + 1) * powerMod.nbrLimbs], powerMod.nbrLimbs * sizeof(int));
        pstFactorInfoRecord++;
      }
      if (currentFactor == nbrFactorsRecord)
      {            // All factors found.
        break;
      }
      ptrTrailingCoeff = &polyLifted[(currentFactor + 1) * powerMod.nbrLimbs];
      UncompressIntLimbs(pstFactorInfoRecord->ptrPolyLifted, (limb*)nbrTmp2, powerMod.nbrLimbs);
      MultBigNbrModN(ptrTrailingCoeff, nbrTmp2, nbrTmp3, (int*)powerMod.limbs, powerMod.nbrLimbs);
      do
      {
        ptrTrailingCoeff = &polyLifted[currentFactor * powerMod.nbrLimbs];
        memcpy(ptrTrailingCoeff, nbrTmp3, powerMod.nbrLimbs * sizeof(int));
      } while (--currentFactor >= 0);
      NumberLength = powerMod.nbrLimbs;
      UncompressLimbsBigInteger((limb*)ptrTrailingCoeff, &operand1);
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
          int* ptrValue1 = pstFactorInfoRecord->ptrPolyLifted;    // Source
          int* ptrValue2 = poly2;                                 // Destination
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
      CompressLimbsBigInteger((limb*)nbrTmp, &leadingCoeff);
      for (currentDegree = 0; currentDegree <= degreeProd; currentDegree++)
      {
        UncompressIntLimbs(ptrSrc, (limb*)nbrTmp2, powerMod.nbrLimbs);
        MultBigNbrModN(nbrTmp, nbrTmp2, nbrTmp3, (int*)powerMod.limbs, powerMod.nbrLimbs);
        CompressIntLimbs(ptrDest, (limb*)nbrTmp3, powerMod.nbrLimbs + 1);
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
        NumberLength = operand2.nbrLimbs;
        CompressBigInteger(ptrDest, &operand2);
        ptrSrc += 1 + numLimbs(ptrSrc);
        ptrDest += 1 + numLimbs(ptrDest);
      }
      rc = DivideIntegerPolynomial(polyS, poly5, TYPE_MODULUS);
      if (rc == EXPR_POLYNOMIAL_DIVISION_NOT_INTEGER)
      {
        continue;    // Cannot perform the division.
      }
      if (polyS[0] != 0 || polyS[1] != 1 || polyS[2] != 0)
      {              // Remainder is not zero. Number to factor does not divide this polynomial.
        continue;
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
        NumberLength = operand3.nbrLimbs;
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


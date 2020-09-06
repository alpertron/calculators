/*
This file is part of Alpertron Calculators.

Copyright 2019 Dario Alejandro Alpern

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
#include "musl.h"

#define MAX_MATRIX_SIZE  100
#define MAX_LEN_MINI     700   // 5000 digits

char* ptrOutput2;
#if DEBUG_VANHOEIJ
char debugOutput[20000000];
char* ptrDebugOutput = debugOutput;
#endif

#ifdef __EMSCRIPTEN__
#define LF "<br>"
#else
#define LF "\r\n"
#endif

typedef struct MiniBigInteger
{
  int nbrLimbs;
  enum eSign sign;
  limb limbs[MAX_LEN_MINI];
} MiniBigInteger;

MiniBigInteger basis[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE];
MiniBigInteger basisStar[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE];
MiniBigInteger lambda[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE];
MiniBigInteger prodB[MAX_MATRIX_SIZE];
MiniBigInteger detProdB[MAX_MATRIX_SIZE];
MiniBigInteger matrixBL[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE];
BigInteger traces[MAX_MATRIX_SIZE*10];
BigInteger powerBoundA, powerExtraBits;
int *ptrCoeffs[MAX_MATRIX_SIZE];
static int* ptrFactorInteger;

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
static void GenerateIntegerPolynomial(int* polyMod, int* polyInt, int degreePoly);
static void InsertIntegerPolynomialFactor(int* ptrFactor, int degreePoly);

// Generate row echelon form from matrixBL.
// The output will be located in matrix lambda.
static int gauss(int nbrCols, int nbrRows)
{
  int k, l;
  int row, col;
  int pos;

#ifdef __EMSCRIPTEN__
  char outputInfo[1000];
  char* ptrOutput = outputInfo;
  if (lang)
  {
    strcpy(ptrOutput, "1<p>Resolviendo matriz de");
  }
  else
  {
    strcpy(ptrOutput, "1<p>Solving matrix of");
  }
  ptrOutput += strlen(ptrOutput);
  *ptrOutput++ = ' ';
  int2dec(&ptrOutput, nbrRows);
  strcpy(ptrOutput, " &times; ");
  ptrOutput += strlen(ptrOutput);
  int2dec(&ptrOutput, nbrCols);
  strcpy(ptrOutput, ".</p>");
  ptrOutput += strlen(ptrOutput);
  databack(outputInfo);
#endif
  for (row = 0; row < nbrRows; row++)
  {
    for (col = 0; col < nbrCols; col++)
    {
      CopyBigInt((BigInteger*)&lambda[row][col], (BigInteger*)&matrixBL[col][row]);
    }
  }

  l = 0;
  for (k = 0; k < nbrCols && l < nbrRows; k++)
  {
    pos = -1;
    /* Look for a pivot under the diagonal. */
    row = l;
    while (row < nbrRows && BigIntIsZero((BigInteger*)&lambda[row][k]))
    {
      row++;
    }

    if (row == nbrRows)
      /* No pivot found; try to find a unique 1 above */
    {
      row = 0;
      while (row < l && BigIntIsZero((BigInteger *)&lambda[row][k]))
      {
        row++;
      }
      if (row == l || !BigIntIsOne((BigInteger*)&lambda[row][k]))
      {
        return -1;           // No pivot found, go out.
      }
      row++;
      for (; row < l; row++)
      {
        if (!BigIntIsZero((BigInteger*)&lambda[row][k]))
        {
          return -1;
        }
      }
      continue;
    }

    // We have found a non-zero element on the k-th column
    for (; row < nbrRows && pos == -1; row++)
    {
      BigInteger* pLambda = (BigInteger*)&lambda[row][k];
      if (pLambda->nbrLimbs == 1 && pLambda->sign == SIGN_NEGATIVE &&
        pLambda->limbs[0].x == 1)   // Value is -1.
      {   // invert all elements on i-th row
        for (col = k; col < nbrCols; col++)
        {
          BigIntChSign((BigInteger*)&lambda[row][col]);
        }
        pos = row;
      }
      if (BigIntIsOne(pLambda))
      {
        pos = row;
      }
    }

    if (pos != -1)
    { /* necessarily M[pos][k]=1 */
      for (col = 0; col < nbrCols; col++)
      {
        CopyBigInt(&tmp0, (BigInteger *)&lambda[pos][col]);
        CopyBigInt((BigInteger*)&lambda[pos][col], (BigInteger*)&lambda[l][col]);
        CopyBigInt((BigInteger*)&lambda[l][col], &tmp0);
      }
      for (row = 0; row < nbrRows; row++)
      {
        // M[i] = M[i] - M[l]
        CopyBigInt(&tmp1, (BigInteger*)&lambda[row][k]);
        if (row != l && !BigIntIsZero(&tmp1))
        {
          intToBigInteger((BigInteger*)&lambda[row][k], 0);
          for (col = k + 1; col < nbrCols; col++)
          {
            // *x = *x + (*y)*t1
            BigIntMultiply((BigInteger*)&lambda[l][col], &tmp1, &tmp0);
            BigIntSubt((BigInteger*)&lambda[row][col], &tmp0, (BigInteger*)&lambda[row][col]);
          }
        }
      }
      l++;
    }
    else
    {
      return -1;    // No 1 found as pivot.
    }
  }
  return l;
}


// Compute the matrix lambda of Gram-Schimdt orthogonalization.
// Use algorithm 2 of "Symplectic Lattice Reduction and NTRU"
// (Gama1, Howgrave - Graham, Nguyen)
static void GramSchmidtOrthogonalization(int nbrRows, int nbrCols)
{
  int colI, colJ, k, row;
#ifdef __EMSCRIPTEN__
  char outputInfo[1000];
  char* ptrOutput = outputInfo;
  if (lang)
  {
    strcpy(ptrOutput, "1<p>Calculando Gram-Schmidt en matriz de");
  }
  else
  {
    strcpy(ptrOutput, "1<p>Computing Gram-Schmidt in matrix of");
  }
  ptrOutput += strlen(ptrOutput);
  *ptrOutput++ = ' ';
  int2dec(&ptrOutput, nbrRows);
  strcpy(ptrOutput, " &times; ");
  ptrOutput += strlen(ptrOutput);
  int2dec(&ptrOutput, nbrCols);
  strcpy(ptrOutput, ".</p>");
  ptrOutput += strlen(ptrOutput);
  databack(outputInfo);
#endif
  for (colI = 0; colI < nbrCols; colI++)
  {
    if (colI == 0)
    {             // U_0 <- 1
      intToBigInteger((BigInteger*)&detProdB[0], 1);
    }
    else
    {             // U_i <- lambda_{i-1, i-1}; U_{i-1} <- -lambda_{i, i-1}
      CopyBigInt((BigInteger*)&detProdB[colI], (BigInteger*)&lambda[colI - 1][colI - 1]);
      CopyBigInt((BigInteger*)&detProdB[colI - 1], (BigInteger*)&lambda[colI][colI - 1]);
      BigIntChSign((BigInteger*)&detProdB[colI - 1]);
    }
    for (colJ = colI - 2; colJ >= 0; colJ--)
    {             // U_j <- 0
      intToBigInteger(&tmp1, 0);
      for (k = colJ + 1; k <= colI; k++)
      {           // U_j <- U_j + lambda_{k,j} * U_k
        BigIntMultiply((BigInteger*)&lambda[k][colJ], (BigInteger*)&detProdB[k], &tmp2);
        BigIntSubt(&tmp1, &tmp2, &tmp1);
      }
      // U_j <- U_j / lambda_{j,j}
      BigIntDivide(&tmp1, (BigInteger*)&lambda[colJ][colJ], (BigInteger*)&detProdB[colJ]);
    }
    for (colJ = colI; colJ < nbrCols; colJ++)
    {           // lambda_{j,i} <- 0
      BigInteger* plambdaJI = (BigInteger*)&lambda[colJ][colI];
      intToBigInteger(plambdaJI, 0);
      for (k = 0; k <= colI; k++)
      {         // tmp2 <- scalar product b_j * b_k 
        intToBigInteger(&tmp2, 0);
        for (row = 0; row < nbrRows; row++)
        {
          BigIntMultiply((BigInteger*)&basisStar[row][colJ], (BigInteger*)&basisStar[row][k], &tmp3);
          BigIntAdd(&tmp2, &tmp3, &tmp2);
        }
        // lambda_{j,i} <- lambda_{j,i} + tmp2 * U_k
        BigIntMultiply(&tmp2, (BigInteger*)&detProdB[k], &tmp2);
        BigIntAdd(plambdaJI, &tmp2, plambdaJI);
      }
    }
  }
}

static void PerformREDI(int k, int l, int size)
{
  int i, row;
  // If |2 lambda_{k, l}| <= d_l, go out.
  BigInteger* pLambda = (BigInteger *)&lambda[k][l];
  BigInteger* pDet = (BigInteger *)&detProdB[l];
  multint(&tmp0, pLambda, 2);        // tmp0 <- 2 lambda_{k, l}
  tmp0.sign = SIGN_POSITIVE;         // tmp0 <- |2 lambda_{k, l}|
  BigIntSubt(pDet, &tmp0, &tmp1);    // tmp1 <- d_l - |2 lambda_{k, l}|
  if (tmp1.sign == SIGN_POSITIVE)
  {
    return;                          // Go out if d_l >= |2 lambda_{k, l}|
  }
  // Compute q = nearest integer to lambda_{k, l} / d_l.
  // d_l is always positive.
  // If lambda is positive, compute (2 lambda_{k, l} + d_l) / (2 d_l)
  // If lambda is negative, compute (2 lambda_{k, l} - d_l) / (2 d_l)
  multint(&tmp0, pLambda, 2);        // tmp0 <- 2 lambda_{k, l}
  if (pLambda->sign == SIGN_POSITIVE)
  {
    BigIntAdd(&tmp0, pDet, &tmp0);   // tmp0 <- 2 lambda_{k, l} + d_l
  }
  else
  {
    BigIntSubt(&tmp0, pDet, &tmp0);  // tmp0 <- 2 lambda_{k, l} - d_l
  }
  multint(&tmp1, pDet, 2);           // tmp1 <- 2 d_l
  BigIntDivide(&tmp0, &tmp1, &tmp2); // tmp2 <- q.
  for (row = 0; row < size; row++)
  {  // b_k <- b_k - q*b_l
    BigIntMultiply(&tmp2, (BigInteger *)&basis[row][l], &tmp0);  // tmp0 <- q*b_l
    BigIntSubt((BigInteger *)&basis[row][k], &tmp0, (BigInteger *)&basis[row][k]);
  }
    // lambda_{k, l} <- lambda_{k, l} - q*d_l
  BigIntMultiply(&tmp2, (BigInteger *)&detProdB[l], &tmp0);      // tmp0 <- q*d_l
  BigIntSubt((BigInteger *)&lambda[k][l], &tmp0, (BigInteger *)&lambda[k][l]);
  for (i = 0; i < l; i++)
  { // lambda_{k, i} <- lambda_{k, i} - q*lambda_{l, i}
    BigIntMultiply(&tmp2, (BigInteger *)&lambda[l][i], &tmp0);   // tmp0 <- q*lambda_{l, i}
    BigIntSubt((BigInteger *)&lambda[k][i], &tmp0, (BigInteger *)&lambda[k][i]);
  }
}

static void PerformSWAPI(int k, int kMax, int size)
{                   // On entry: k >= 1
  int row, i, j;
  // Exchange b_k with b_{k-1}
  for (row = 0; row < size; row++)
  {
    CopyBigInt(&tmp0, (BigInteger *)&basis[row][k]);
    CopyBigInt((BigInteger *)&basis[row][k], (BigInteger *)&basis[row][k-1]);
    CopyBigInt((BigInteger *)&basis[row][k-1], &tmp0);
  }
  if (k > 1)
  {
    for (j = 0; j < k - 1; j++)
    {  // Exchange lambda_{k, j} with lambda_{k-1, j}
      CopyBigInt(&tmp0, (BigInteger*)&lambda[k][j]);
      CopyBigInt((BigInteger*)&lambda[k][j], (BigInteger*)&lambda[k - 1][j]);
      CopyBigInt((BigInteger*)&lambda[k - 1][j], &tmp0);
    }
  }
    // Set lambda <- lambda_{k, k-1}
  CopyBigInt(&tmp0, (BigInteger *)&lambda[k][k - 1]);    // tmp0 <- lambda.
    // Set B <- (d_{k-2}*d_k + lambda^2)/d_{k-1}
  if (k > 1)
  {
    BigIntMultiply((BigInteger *)&detProdB[k - 2], (BigInteger *)&detProdB[k], &tmp1);
  }
  else
  {                                     // detProdB[k] = 1 for k < 0.
    CopyBigInt(&tmp1, (BigInteger *)&detProdB[k]);
  }
  BigIntMultiply(&tmp0, &tmp0, &tmp2);  // tmp2 <- lambda^2
  BigIntAdd(&tmp1, &tmp2, &tmp1);       // tmp1 <- d_{k-2}*d_k + lambda^2
  BigIntDivide(&tmp1, (BigInteger*)&detProdB[k - 1], &tmp1); // tmp1 <- B.
  for (i = k+1; i <= kMax; i++)
  {
    // t <- lambda_{i, k}
    CopyBigInt(&tmp2, (BigInteger *)&lambda[i][k]);   // tmp2 <- t
    // lambda_{i, k} <- (d_k*lambda_{i, k-1} - lambda * t)/d_{k-1}
    BigIntMultiply((BigInteger *)&detProdB[k], (BigInteger *)&lambda[i][k-1], &tmp3);
    BigIntMultiply(&tmp0, &tmp2, &tmp4);       // tmp4 <- lambda * t
    BigIntSubt(&tmp3, &tmp4, &tmp3);
    BigIntDivide(&tmp3, (BigInteger*)&detProdB[k - 1], (BigInteger*)&lambda[i][k]);
    // lambda_{i, k-1} <- (B*t + lambda * lambda_{i, k})/d_k
    BigIntMultiply(&tmp1, &tmp2, &tmp3);       // tmp3 <- B * t
    BigIntMultiply(&tmp0, (BigInteger *)&lambda[i][k], &tmp4);
    BigIntAdd(&tmp3, &tmp4, &tmp3);
    BigIntDivide(&tmp3, (BigInteger *)&detProdB[k], (BigInteger *)&lambda[i][k - 1]);
  }
  CopyBigInt((BigInteger *)&detProdB[k - 1], &tmp1);  // d_{k-1} <- B.
}

// Use algorithm 2.6.7 of Henri Cohen's book A Course in Computational Algebraic Number Theory
// B_0, B_1, ... are the different columns of matrix basis.
// d_k = 1 when k < 0.
void integralLLL(int size)
{
  int k, kMax, row;
  int colI, colJ;
  int l;
#ifdef __EMSCRIPTEN__
  char outputInfo[1000];
  char* ptrOutput = outputInfo;
  if (lang)
  {
    strcpy(ptrOutput, "1<p>Calculando LLL en matriz de");
  }
  else
  {
    strcpy(ptrOutput, "1<p>Computing LLL in matrix of");
  }
  ptrOutput += strlen(ptrOutput);
  *ptrOutput++ = ' ';
  int2dec(&ptrOutput, size);
  strcpy(ptrOutput, " &times; ");
  ptrOutput += strlen(ptrOutput);
  int2dec(&ptrOutput, size);
  strcpy(ptrOutput, ".</p>");
  ptrOutput += strlen(ptrOutput);
  databack(outputInfo);
#endif
  k = 1;
  kMax = 0;
  // Set d_0 to scalar product B_0 * B_0
  intToBigInteger((BigInteger *)&detProdB[0], 0);
  for (row = 0; row < size; row++)
  {
    BigIntMultiply((BigInteger*)&basis[row][0], (BigInteger*)&basis[row][0], &tmp0);
    BigIntAdd((BigInteger*)&detProdB[0], &tmp0, (BigInteger*)&detProdB[0]);
  }
  do
  {
    // Perform incremental Gram-Schmidt
    if (k > kMax)
    {
      kMax = k;
      for (colJ = 0; colJ <= k; colJ++)
      {
        // Compute u <- scalar product of b_k and b_j (use tmp2 for u)
        intToBigInteger(&tmp2, 0);
        for (row = 0; row < size; row++)
        {
          BigIntMultiply((BigInteger *)&basis[row][k], (BigInteger *)&basis[row][colJ], &tmp0);
          BigIntAdd(&tmp2, &tmp0, &tmp2);
        }
        for (colI = 0; colI < colJ; colI++)
        {    // Set u <- (d_i * u - lambda_{k, i} * lambda_{j, i}) / d_{i-1}
          BigIntMultiply((BigInteger *)&detProdB[colI], &tmp2, &tmp0);
          BigIntMultiply((BigInteger *)&lambda[k][colI], (BigInteger *)&lambda[colJ][colI], &tmp1);
          if (colI == 0)
          {
            BigIntSubt(&tmp0, &tmp1, &tmp2);
          }
          else
          {
            BigIntSubt(&tmp0, &tmp1, &tmp0);
            BigIntDivide(&tmp0, (BigInteger *)&detProdB[colI - 1], &tmp2);
          }
        }
        if (colJ < k)
        {
          CopyBigInt((BigInteger *)&lambda[k][colJ], &tmp2);
        }
        else
        {
          CopyBigInt((BigInteger *)&detProdB[k], &tmp2);
        }
      }
    }
    // Test LLL condition
    for (;;)
    {
      PerformREDI(k, k - 1, size);
      // Check whether 4 * d_k * d_{k-2} < 3 (d_{k-1})^2 - 4*(lambda_{k, k-1})^2
      // that means 4 (d_k * d_{k-2} + lambda_{k, k-1})^2) < 3 (d_{k-1})^2 
      if (k > 1)
      {
        BigIntMultiply((BigInteger*)&detProdB[k], (BigInteger*)&detProdB[k - 2], &tmp0);
      }
      else
      {                // detProdB[k] = 1 for k < 0.
        CopyBigInt(&tmp0, (BigInteger*)&detProdB[k]);
      }
      BigIntMultiply((BigInteger*)&lambda[k][k - 1], (BigInteger*)&lambda[k][k - 1], &tmp1);
      BigIntAdd(&tmp0, &tmp1, &tmp0);
      multint(&tmp0, &tmp0, 4);           // tmp0 = Left Hand Side.
      BigIntMultiply((BigInteger *)&detProdB[k - 1], (BigInteger *)&detProdB[k - 1], &tmp1);
      multint(&tmp1, &tmp1, 3);
      BigIntSubt(&tmp0, &tmp1, &tmp0);
      if (tmp0.sign == SIGN_POSITIVE)
      {
        break;                            // Exit loop if condition above does not hold.
      }
      PerformSWAPI(k, kMax, size);
      if (k > 1)                          // k <- max(1,k-1)
      {
        k--;
      }
    }
    for (l = k - 2; l >= 0; l--)
    {
      PerformREDI(k, l, size);
    }
  } while (++k < size);
}

// Compute remainder such that the result is in range -divisor/2 to divisor/2.
static void BigIntSymmetricRemainder(BigInteger* dividend, BigInteger* divisor, BigInteger* result)
{
  BigIntRemainder(dividend, divisor, result);
  // Convert trace to range -powerExtraBits/2 to powerExtraBits/2
  if (result->sign == SIGN_NEGATIVE)
  {
    BigIntAdd(result, divisor, result);
  }
  multint(&operand1, result, 2);
  BigIntSubt(&operand1, divisor, &operand1);
  if (operand1.sign == SIGN_POSITIVE)
  {
    BigIntSubt(result, divisor, result);
  }
}

// Compute the traces of all modular factors in a matrix.
// traces[row] = trace_row(poly_col). Column number is a parameter.
// Let the polynomial factor be:
// P(x) = x^d + E_1*x^(d-1) + E_2*x^(d-2) + ... + E_n
// Tr_n(P) = -n*E_n - E_{n-1}*Tr_1(P) - E_{n-2}*Tr_2(P) - ...
// The coefficients have to be processed in the order inverted,
// so the first step is to get the pointers to the coefficients.
//
static void ComputeTraces(int nbrTraces, int nbrCol)
{
  int traceNbr;
  int* ptrCoeff, *ptrCoeffG;
  int polyDegree, degree;
  int* ptrTemp = tempPoly;
  BigInteger* ptrTrace;
  struct sFactorInfo* pstFactorInfo = &factorInfoRecord[nbrCol];
  // Get pointers to the coefficients.
  polyDegree = pstFactorInfo->degree;          // Degree of polynomial factor
  ptrCoeff = pstFactorInfo->ptrPolyLifted;     // coefficients of polynomial factor.
  for (traceNbr = polyDegree - 1; traceNbr >= 0; traceNbr--)
  {
    ptrCoeffs[traceNbr] = ptrTemp;
    operand4.nbrLimbs = *ptrCoeff;
    operand4.sign = SIGN_POSITIVE;
    memcpy(operand4.limbs, ptrCoeff+1, *ptrCoeff *sizeof(int));
    BigIntRemainder(&operand4, &powerMod, &operand3);
    *ptrTemp = operand3.nbrLimbs;
    memcpy(ptrTemp+1, operand3.limbs, *ptrTemp * sizeof(int));
    ptrTemp += 1 + numLimbs(ptrTemp);
    ptrCoeff += 1 + numLimbs(ptrCoeff);
  }
  // Compute coefficients of P.
  ptrCoeffG = (int*)basis;
  intToBigInteger(&operand1, 1);
#if DEBUG_VANHOEIJ
  strcpy(ptrDebugOutput, "Coefficients: ");
  ptrDebugOutput += strlen(ptrDebugOutput);
#endif
  for (traceNbr = 0; traceNbr < polyDegree; traceNbr++)
  {
    UncompressBigIntegerB(ptrCoeffs[traceNbr], &operand3);
#if DEBUG_VANHOEIJ
    BigInteger2Dec(&operand3, ptrDebugOutput, 0);
    ptrDebugOutput += strlen(ptrDebugOutput);
    *ptrDebugOutput++ = ',';
    *ptrDebugOutput++ = ' ';
#endif
    BigIntMultiply(&operand3, &operand1, &operand3);
    BigInteger2IntArray(ptrCoeffG, &operand3);
    ptrCoeffs[traceNbr] = ptrCoeffG;
    ptrCoeffG += 1 + numLimbs(ptrCoeffG);
  }
#if DEBUG_VANHOEIJ
  strcpy(ptrDebugOutput, LF);
  ptrDebugOutput += strlen(ptrDebugOutput);
#endif
  // Store traces of polynomial P in matrix traces.
  ptrTrace = &traces[0];
  intToBigInteger(ptrTrace, polyDegree);
  BigIntRemainder(ptrTrace, &powerMod, &traces[0]);
  for (traceNbr = 1; traceNbr < nbrTraces; traceNbr++)
  {
    ptrTrace++;
    // Initialize trace to -n*E_n.
    if (traceNbr <= polyDegree)
    {
      UncompressBigIntegerB(ptrCoeffs[traceNbr - 1], &operand1);
      multint(ptrTrace, &operand1, -traceNbr);
    }
    else
    {
      intToBigInteger(ptrTrace, 0);
    }
    for (degree = 1; degree < traceNbr; degree++)
    {   // Subtract - E_{n-deg}*Tr_deg(P)
      if (traceNbr - degree <= polyDegree)
      {
        UncompressBigIntegerB(ptrCoeffs[traceNbr - degree - 1], &operand3);
        BigIntMultiply(&traces[degree], &operand3, &operand1);
        BigIntSubt(ptrTrace, &operand1, ptrTrace);
        // Get remainder of quotient of result by p^a.
        BigIntRemainder(ptrTrace, &powerMod, ptrTrace);
      }
    }
  }
  ptrTrace = &traces[0];
  for (traceNbr = 1; traceNbr < nbrTraces; traceNbr++)
  {       // Loop that divides all traces by p^b and get modulus p^(a-b).
    ptrTrace++;
    // Get remainder of quotient of trace divided by p^b
    BigIntSymmetricRemainder(ptrTrace, &powerBoundA, &operand2);
    // Subtract the remainder from the computed trace.
    BigIntSubt(ptrTrace, &operand2, ptrTrace);
    // Divide the result by p^b
    BigIntDivide(ptrTrace, &powerBoundA, ptrTrace);
    // Get remainder of quotient of result divided by p^(a-b)
    BigIntSymmetricRemainder(ptrTrace, &powerExtraBits, ptrTrace);
  }
}

static int AttemptToFactor(int nbrVectors, int nbrFactors, int *pNbrFactors)
{
  int currentDegree, currentFactor;
  int nbrVector;
  int stepNbr;
  int degreeProd, degreeFactor;
  int rc;
  struct sFactorInfo* pstFactorInfo;
  int* ptrSrc, * ptrDest;
  int nbrTmp[1000], nbrTmp2[1000], nbrTmp3[1000];
  // In step 1 we check that all factors can be found. If this succeeds,
  // in step 2 we insert factors in final array of polynomial factors.
  for (stepNbr = 1; stepNbr <= 2; stepNbr++)
  {
#if DEBUG_VANHOEIJ
    sprintf(ptrDebugOutput, "stepNbr = %d\n", stepNbr);
    ptrDebugOutput += strlen(ptrDebugOutput);
#endif
    NumberLength = powerMod.nbrLimbs;
    for (nbrVector = 0; nbrVector < nbrVectors; nbrVector++)
    {
      // Test whether the product of the polynomial factors
      // defined by this vector divides the original polynomial.
      // First reduce coefficients mod powerMod and then
      // convert them to Montgomery notation.
      modulusIsZero = 0;  // Perform modular operations.
      NumberLength = powerMod.nbrLimbs;
      degreeProd = 0;     // Set product to 1.
      pstFactorInfo = factorInfoRecord;
      for (currentFactor = 0; currentFactor < nbrFactors; currentFactor++)
      {
        if (!BigIntIsZero((BigInteger*)&lambda[nbrVector][currentFactor]))
        {                 // Multiply this polynomial.
          int* ptrCoeffSrc = pstFactorInfo->ptrPolyLifted;    // Source
          int* ptrCoeffDest = poly2;                          // Destination
          int nbrLength;
          degreeFactor = pstFactorInfo->degree;
          // Reduce coefficients mod powerMod and store them on poly2.
          for (currentDegree = 0; currentDegree < degreeFactor; currentDegree++)
          {
            nbrLength = numLimbs(ptrCoeffSrc);
            if (nbrVectors == 1)
            {              // Coefficient is already reduced.
              *ptrCoeffDest = nbrLength;
              memcpy(ptrCoeffDest + 1, ptrCoeffSrc + 1, nbrLength * sizeof(int));
            }
            else
            {              // Reduce the coefficient mod powerMod.
              operand1.nbrLimbs = nbrLength;
              operand1.sign = SIGN_POSITIVE;
              memcpy(operand1.limbs, ptrCoeffSrc + 1, nbrLength * sizeof(int));
              BigIntRemainder(&operand1, &powerMod, &operand2);
              *ptrCoeffDest = operand2.nbrLimbs;
              memcpy(ptrCoeffDest + 1, operand2.limbs, operand2.nbrLimbs * sizeof(int));
            }
            ptrCoeffSrc += 1 + nbrLength;
            ptrCoeffDest += 1 + NumberLength;
          }
          *ptrCoeffDest = 1;            // Store 1 as the leading coefficient.
          *(ptrCoeffDest + 1) = 1;
          // Convert factor to Montgomery notation.
          polyToMontgomeryNotation(poly2, degreeFactor + 1);
          if (degreeProd == 0)
          {
            ptrCoeffSrc = poly2;              // Source is the new factor.
          }
          else
          {
            MultPolynomial(degreeProd, degreeFactor, poly1, poly2);
            ptrCoeffSrc = polyMultTemp;       // Source is the product
          }
          degreeProd += degreeFactor;
          ptrCoeffDest = poly1;                     // Destination
          degreeFactor = pstFactorInfo->degree;
          for (currentDegree = 0; currentDegree <= degreeProd; currentDegree++)
          {
            nbrLength = 1 + numLimbs(ptrCoeffSrc);
            memcpy(ptrCoeffDest, ptrCoeffSrc, nbrLength * sizeof(int));
            ptrCoeffSrc += 1 + NumberLength;
            ptrCoeffDest += 1 + NumberLength;
          }
          if (stepNbr == 2)
          {
            pstFactorInfo->multiplicity = 0;
            (*pNbrFactors)--;
          }
        }
        pstFactorInfo++;
      }
      // Multiply all coefficients by leadingCoeff by using modmult
      // (this converts from Montgomery to standard notation)
      // and store them in poly2.
      ptrSrc = poly1;
      ptrDest = poly2;
      CompressLimbsBigInteger((limb*)nbrTmp, &leadingCoeff);
      for (currentDegree = 0; currentDegree <= degreeProd; currentDegree++)
      {
        LenAndLimbs2ArrLimbs(ptrSrc, (limb*)nbrTmp2, NumberLength);
        modmult((limb*)nbrTmp, (limb*)nbrTmp2, (limb*)nbrTmp3);
        ArrLimbs2LenAndLimbs(ptrDest, (limb*)nbrTmp3, NumberLength + 1);
        ptrSrc += 1 + NumberLength;
        ptrDest += 1 + numLimbs(ptrDest);
      }
      GenerateIntegerPolynomial(poly2, poly5, degreeProd);
      // Ensure that the absolute value of the coefficients 
      // are less than the bound.
      ptrSrc = &poly5[1];
      for (currentDegree = 0; currentDegree <= poly5[0]; currentDegree++)
      {
        UncompressBigIntegerB(ptrSrc, &operand1);
        BigIntSubt(&operand1, &bound, &operand2);
        if (operand2.sign == SIGN_POSITIVE)
        {
          return 0;     // Coefficient is too high.
        }
        BigIntAdd(&operand1, &bound, &operand2);
        if (operand2.sign == SIGN_NEGATIVE)
        {
          return 0;     // Coefficient is too low.
        }
        ptrSrc += 1 + numLimbs(ptrSrc);
      }
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
        BigInteger2IntArray(ptrDest, &operand2);
        ptrSrc += 1 + numLimbs(ptrSrc);
        ptrDest += 1 + numLimbs(ptrDest);
      }
      if (stepNbr == 1)
      {   // Test whether the product divides the original polynomial.
        UncompressBigIntegerB(&polyS[1], &operand1);
        UncompressBigIntegerB(&poly5[1], &operand2);
        BigIntRemainder(&operand1, &operand2, &operand3);
        if (!BigIntIsZero(&operand3))
        {    // Constant term of product does not divide constant term of
             // original polynomial, so go out.
          return 0;
        }
        rc = DivideIntegerPolynomial(polyS, poly5, TYPE_MODULUS);
        if (rc == EXPR_POLYNOMIAL_DIVISION_NOT_INTEGER)
        {               // Cannot perform the division.
          return 0;
        }
        if (polyS[0] != 0 || polyS[1] != 1 || polyS[2] != 0)
        {              // Remainder is not zero.
          return 0;    // Number to factor does not divide this polynomial.
        }
      }
      else
      {
        int degreePoly;
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
          BigInteger2IntArray(ptrDest, &operand3);
          ptrSrc += 1 + numLimbs(ptrSrc);
          ptrDest += 1 + numLimbs(ptrDest);
        }
        // Copy this principal part to poly5.
        CopyPolynomial(&poly5[1], &poly2[1], poly5[0]);
        DivideIntegerPolynomial(polyNonRepeatedFactors, poly5, TYPE_DIVISION);
        degreePoly = poly5[0];
        CopyPolynomial(ptrFactorInteger, &poly5[1], degreePoly);
        InsertIntegerPolynomialFactor(ptrFactorInteger, degreePoly);
      }
      modulusIsZero = 0;   // Perform modular operations.
      // Restart finding factors.
      memset(arrNbrFactors, 0, sizeof(arrNbrFactors));
    }
  }
  return 1;              // Factorization is complete.
}

static void CompressPolynomial(int nbrFactorsFound)
{
  int nbrFactor;
  struct sFactorInfo* pstFactorInfo = factorInfoRecord;
  for (nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    int* ptrDest = pstFactorInfo->ptrPolyLifted;
    int* ptrSrc = ptrDest;
    for (degree = 0; degree < pstFactorInfo->degree; degree++)
    {
      int coeffLength = 1 + numLimbs(ptrSrc);
      memmove(ptrDest, ptrSrc, coeffLength * sizeof(int));
      ptrDest += coeffLength;
      ptrSrc += NumberLength + 1;
    }
    pstFactorInfo++;
  }
}

// Perform Van Hoeij algorithm 
// On input: values = integer polynomial.
// prime: integer prime used for factoring the polynomial
// primeMod: BigInteger representing prime.
// exponentMod: Exponent used for Knuth-Cohen bound.
// powerMod: BigInteger representing Knuth-Cohen bound.
//
// Use algorithm in paper "Tuning and generalizing van Hoeij's algorithm"
// by Karim Belabas, Gauillaume Hanrot and Paul Zimmermann (2001)
//
// Relation from variables in paper to variables here:
// n_0 = nbrFactors
// n = nbrColsTraces
// r = nbrVectors (number of columns of BL)
// s = nbrRequiredTraces
static void vanHoeij(int prime, int nbrFactors)
{
#if DEBUG_VANHOEIJ
  int nbrStepsDone = 0;
#endif
  int nbrVectors, nbrColsTraces, r1;
  int nbrRequiredTraces, firstTrace;
  int nbrRow, nbrCol;
  int exponDifference;
  int delta;
  // Store into polyLifted the polynomial values / lc(values) (mod powerMod).
  double logPrime = log(prime);
  double b0 = log(2 * values[0]) / logPrime;
  double log_rootbound = logBigNbr(&bound) / logPrime;
  int a0;
  int b, C, j, rank;
  int degree1;
  int* ptrSrc;
  struct sFactorInfo* pstFactorInfo;
  int degreePolyToFactor = polyNonRepeatedFactors[0];
  int nbrFactor;
  int newNumberLength, newNbrFactors;
  int ctr1, ctr2, ctr3;

  exponentMod = (int)ceil(log_rootbound) + 2;
#if DEBUG_VANHOEIJ
  sprintf(ptrDebugOutput, LF "====================================================="
    LF "prime = %d, root bound = ", prime);
  ptrDebugOutput += strlen(ptrDebugOutput);
  BigInteger2Dec(&bound, ptrDebugOutput, 0);
  ptrDebugOutput += strlen(ptrDebugOutput);
  *ptrDebugOutput++ = '\n';
#endif
  memset(arrNbrFactors, 0, sizeof(arrNbrFactors));
  // Get leading coefficient of polyNonRepeatedFactors.
  ptrSrc = &polyNonRepeatedFactors[1];
  for (degree1 = 0; degree1 < degreePolyToFactor; degree1++)
  {
    ptrSrc += 1 + numLimbs(ptrSrc);
  }
  computePower(exponentMod);
  modulusIsZero = 0;    // Use modular arithmetic for polynomials.
  intToBigInteger(&operand5, 1);
  values[0] = polyNonRepeatedFactors[0];
  getModPolynomial(&values[1], polyNonRepeatedFactors, &operand5);
#if DEBUG_HENSEL_LIFTING
  ptrOutput2 = ptrDebugOutput;
#endif
  HenselLifting(factorInfoRecord);
#if DEBUG_HENSEL_LIFTING
  ptrDebugOutput = ptrOutput2; ptrOutput2 = NULL;
#endif
  // Compress polynomials so there are no spaces between coefficients.
  CompressPolynomial(nbrFactorsFound);
  // Check whether the product of up to 3 factors mod p^b
  // is a divisor of the original polynomial.
  CopyBigInt(&halfPowerMod, &powerMod);
  subtractdivide(&halfPowerMod, -1, 2); // halfPowerMod <- (powerMod+1)/2
  GetMontgomeryParms(powerMod.nbrLimbs);
  for (ctr1 = 0; ctr1 < nbrFactors; ctr1++)
  {
    intToBigInteger((BigInteger*)&lambda[0][ctr1], 0);
  }
  newNbrFactors = nbrFactors;
#if 1
  for (ctr1 = 0; ctr1 < nbrFactors; ctr1++)
  {
    intToBigInteger((BigInteger*)&lambda[0][ctr1], 1);
    AttemptToFactor(1, nbrFactors, &newNbrFactors);
    intToBigInteger((BigInteger*)&lambda[0][ctr1], 0);
  }
  for (ctr1 = 0; ctr1 < nbrFactors; ctr1++)
  {
    intToBigInteger((BigInteger*)&lambda[0][ctr1], 1);
    for (ctr2 = ctr1 + 1; ctr2 < nbrFactors; ctr2++)
    {
      intToBigInteger((BigInteger*)&lambda[0][ctr2], 1);
      AttemptToFactor(1, nbrFactors, &newNbrFactors);
      intToBigInteger((BigInteger*)&lambda[0][ctr2], 0);
    }
    intToBigInteger((BigInteger*)&lambda[0][ctr1], 0);
  }
  for (ctr1 = 0; ctr1 < nbrFactors; ctr1++)
  {
    intToBigInteger((BigInteger*)&lambda[0][ctr1], 1);
    for (ctr2 = ctr1 + 1; ctr2 < nbrFactors; ctr2++)
    {
      intToBigInteger((BigInteger*)&lambda[0][ctr2], 1);
      for (ctr3 = ctr2 + 1; ctr3 < nbrFactors; ctr3++)
      {
        intToBigInteger((BigInteger*)&lambda[0][ctr3], 1);
        AttemptToFactor(1, nbrFactors, &newNbrFactors);
        intToBigInteger((BigInteger*)&lambda[0][ctr3], 0);
      }
      intToBigInteger((BigInteger*)&lambda[0][ctr2], 0);
    }
    intToBigInteger((BigInteger*)&lambda[0][ctr1], 0);
  }
#endif
  if (newNbrFactors <= 1)
  {   // Zero or 1 factor left. Polynomial completely factored, so go out.
    return;
  }
  // Discard the entries with multiplicity zero in the stored factors array.
  pstFactorInfo = factorInfoRecord;
  delta = 0;
  for (nbrFactor = 0; nbrFactor < nbrFactors; nbrFactor++)
  {
    if (pstFactorInfo->multiplicity == 0)
    {
      delta++;
    }
    else if (delta > 0)
    {
      memcpy(pstFactorInfo - delta, pstFactorInfo, sizeof(*pstFactorInfo));
    }
    pstFactorInfo++;
  }
  nbrFactors = newNbrFactors;
  nbrFactorsFound = nbrFactors;
  // Compute powerExtraBits as a power of prime with 3 times
  // the number of bits of nbrFactors (the number of polynomial factors).
  // Compute powerBoundA as powerExtraBits * powerMod.
  exponDifference = (int)(4.0 * (double)nbrFactors * log(2.0) / log(prime));
  intToBigInteger(&operand1, prime);
  BigIntPowerIntExp(&operand1, exponDifference, &powerExtraBits);
  CopyBigInt(&powerBoundA, &powerMod);
  BigIntMultiply(&powerExtraBits, &powerMod, &powerBoundA);
  b = (int)(b0 + (ceil(log(3) * log_rootbound) / logPrime) + 3);
  a0 = b + exponDifference;
  exponentMod = a0 * 2 + 50;
  computePower(exponentMod);
  modulusIsZero = 0;    // Use modular arithmetic for polynomials.
  intToBigInteger(&operand5, 1);
  values[0] = polyNonRepeatedFactors[0];
  getModPolynomial(&values[1], polyNonRepeatedFactors, &operand5);
#if DEBUG_HENSEL_LIFTING
  ptrOutput2 = ptrDebugOutput;
#endif
  HenselLifting(factorInfoRecord);
#if DEBUG_HENSEL_LIFTING
  ptrDebugOutput = ptrOutput2; ptrOutput2 = NULL;
#endif
  // Compress polynomials so there are no spaces between coefficients.
  CompressPolynomial(nbrFactorsFound);
  // Initialize matrix BL with identity matrix (length nbrFactors).
  for (nbrRow = 0; nbrRow < nbrFactors; nbrRow++)
  {
    for (nbrCol = 0; nbrCol < nbrFactors; nbrCol++)
    {
      intToBigInteger((BigInteger*)&matrixBL[nbrRow][nbrCol],
        (nbrRow == nbrCol ? 1 : 0));
    }
  }
  // Initialize variables for van Hoeij algorithm.
  nbrRequiredTraces = 1;        // Initialize nbrRequiredTraces.
  firstTrace = 1;               // First trace to use for matrix M.
  nbrVectors = nbrFactors;      // Initialize nbrVectors.
  nbrColsTraces = nbrFactors;
  for (;;)
  {
    int squareFormula;
    modulusIsZero = 0;    // Use modular arithmetic for polynomials.
    firstTrace += nbrRequiredTraces - 1;
    if (nbrVectors >= nbrFactors)
    {
      nbrColsTraces = nbrFactors;
      nbrVectors = nbrFactors;
      // Initialize matrix BL with identity matrix (length nbrFactors).
      for (nbrRow = 0; nbrRow < nbrFactors; nbrRow++)
      {
        for (nbrCol = 0; nbrCol < nbrFactors; nbrCol++)
        {
          intToBigInteger((BigInteger*)&matrixBL[nbrRow][nbrCol],
            (nbrRow == nbrCol ? 1 : 0));
        }
      }
    }
    // Get b such that p^b is greater than the bounds on the cofficients.
    b = (int)(b0 + (ceil(log(firstTrace + nbrRequiredTraces + 1) * log_rootbound) / logPrime) + 3);
    a0 = b + exponDifference;
    computePower(a0);
    newNumberLength = NumberLength;
    exponentMod = a0;
    intToBigInteger(&operand1, prime);
    BigIntPowerIntExp(&operand1, b, &powerBoundA);
#if DEBUG_VANHOEIJ
    sprintf(ptrDebugOutput, LF "====================================================="
      LF "prime = %d, a0 = %d, b0 = %d, exponDifference = %d" LF, prime, a0, b, exponDifference);
    ptrDebugOutput += strlen(ptrDebugOutput);
#endif
    intToBigInteger(&operand5, 1);
    values[0] = polyNonRepeatedFactors[0];
    getModPolynomial(&values[1], polyNonRepeatedFactors, &operand5);
#if DEBUG_HENSEL_LIFTING
    ptrOutput2 = ptrDebugOutput;
#endif
#if DEBUG_HENSEL_LIFTING
    ptrDebugOutput = ptrOutput2; ptrOutput2 = NULL;
#endif
    NumberLength = newNumberLength;
    // use exponDifference additional bits instead of a fixed number

    intToBigInteger(&operand1, prime);
    BigIntPowerIntExp(&operand1, exponDifference, &powerExtraBits);

    C = (int)(ceil(sqrt(nbrRequiredTraces * nbrFactors) / 2));

    // Step 1: Construction of the LLL input matrix.
    // Generate matrix M in matrix basis.
    // The matrix to generate is:
    // ( C * I_r      0     )
    // (   S      p^l * I_s )
    //
    //     ( S_1(G_1)   ...   S_1(G_n)  )
    // S = (   ...      ...     ...     ) * BL
    //     ( S_k(G_1)   ...   S_k(G_n)  )
    // where G_1, G_2, ..., G_n are the p-adic factors of P.
    // I_r is the r*r identity matrix.
    // S_r(G_k) = (T_r(G_i) - T_r(G_i) mod p^b) / p^b (mod p^(a-b))
    // where T_r(G_i) is the trace (Newton sum of powers of roots of
    // polynomial G_i).
    // r = nbrVectors
    // s = nbrRequiredTraces
    // Use basisStar for matrix m.
    for (nbrRow = 0; nbrRow < nbrVectors; nbrRow++)   // Loop for upper half
    {
      for (nbrCol = 0; nbrCol < nbrVectors + nbrRequiredTraces; nbrCol++)
      {
        intToBigInteger((BigInteger*)&basisStar[nbrRow][nbrCol],
          (nbrRow == nbrCol ? C : 0));
      }
    }

    // Use matrix lambda to hold the traces.
    for (int factorNbr = 0; factorNbr < nbrFactors; factorNbr++)
    {
      ComputeTraces(firstTrace + nbrRequiredTraces, factorNbr);
#if DEBUG_VANHOEIJ
      strcpy(ptrDebugOutput, "Traces: ");
      ptrDebugOutput += strlen(ptrDebugOutput);
      for (nbrRow = firstTrace; nbrRow < firstTrace + nbrRequiredTraces; nbrRow++)
      {
        BigInteger2Dec(&traces[nbrRow], ptrDebugOutput, 0);
        ptrDebugOutput += strlen(ptrDebugOutput);
        if (nbrRow < firstTrace + nbrRequiredTraces - 1)
        {
          *ptrDebugOutput++ = ',';
          *ptrDebugOutput++ = ' ';
        }
      }
      strcpy(ptrDebugOutput, LF);
      ptrDebugOutput += strlen(ptrDebugOutput);
#endif
      // Use traces 1 to RequiredTraces
      // (i.e. Tr[1] to Tr[nbrRequiredTraces])
      // Convert traces matrix to range -powerExtraBits/2 to powerExtraBits/2.
      for (j = 0; j < nbrRequiredTraces; j++)
      {
        CopyBigInt(&operand1, &traces[firstTrace + j]);
        BigIntMultiplyBy2(&operand1);
        BigIntSubt(&operand1, &powerExtraBits, &operand1); // operand1 <- 2*traces[j] - powerExtraBits
        if (operand1.sign == SIGN_POSITIVE)
        {
          BigIntSubt(&traces[firstTrace + j], &powerExtraBits, (BigInteger*)&lambda[j][factorNbr]);
        }
        else
        {
          CopyBigInt((BigInteger*)&lambda[j][factorNbr], (BigInteger*)&traces[firstTrace + j]);
        }
      }
    }

    // Multiply traces (in lambda matrix) * BL.
    // Size of matrix trace: nbrRequiredTraces rows and nbrFactors columns.
    // Size of matrix BL: nbrFactors rows and nbrVectors columns.
    // Size of product matrix: nbrRequiredTraces rows and nbrVectors columns.
    // Store product in lower-left part of matrix basisStar.

    for (nbrRow = 0; nbrRow < nbrRequiredTraces; nbrRow++)
    {
      for (nbrCol = 0; nbrCol < nbrVectors; nbrCol++)
      {
        intToBigInteger(&operand1, 0);     // Initialize sum.
        for (j = 0; j < nbrFactors; j++)
        {
          BigIntMultiply((BigInteger*)&lambda[nbrRow][j], (BigInteger*)&matrixBL[j][nbrCol], &operand2);
          BigIntAdd(&operand1, &operand2, &operand1);
        }
        CopyBigInt((BigInteger*)&basisStar[nbrRow + nbrVectors][nbrCol], &operand1);
      }
    }

    // Initialize lower-right half of matrix M.

    for (nbrRow = nbrVectors; nbrRow < nbrVectors + nbrRequiredTraces; nbrRow++)
    {
      for (nbrCol = nbrVectors; nbrCol < nbrVectors + nbrRequiredTraces; nbrCol++)
      {
        if (nbrRow == nbrCol)
        {
          CopyBigInt((BigInteger*)&basisStar[nbrRow][nbrCol], &powerExtraBits);
        }
        else
        {
          intToBigInteger((BigInteger*)&basisStar[nbrRow][nbrCol], 0);
        }
      }
    }

    // Copy matrix M to basis.
    for (nbrRow = 0; nbrRow < nbrVectors + nbrRequiredTraces; nbrRow++)
    {
      for (nbrCol = 0; nbrCol < nbrVectors + nbrRequiredTraces; nbrCol++)
      {
        CopyBigInt((BigInteger*)&basis[nbrRow][nbrCol],
          (BigInteger*)&basisStar[nbrRow][nbrCol]);
      }
    }
#if DEBUG_VANHOEIJ
    strcpy(ptrDebugOutput, LF "Matrix M before LLL: ");
    ptrDebugOutput += strlen(ptrDebugOutput);
    for (nbrRow = 0; nbrRow < nbrVectors + nbrRequiredTraces; nbrRow++)
    {
      for (nbrCol = 0; nbrCol < nbrVectors + nbrRequiredTraces; nbrCol++)
      {
        BigInteger2Dec((BigInteger*)&basisStar[nbrRow][nbrCol], ptrDebugOutput, 0);
        ptrDebugOutput += strlen(ptrDebugOutput);
        if (nbrCol < nbrVectors + nbrRequiredTraces - 1)
        {
          *ptrDebugOutput++ = ',';
        }
      }
      *ptrDebugOutput++ = ';';
    }
    *ptrDebugOutput++ = '.';
#endif
#if DEBUG_VANHOEIJ
    if (++nbrStepsDone == 2)
    {
#ifdef __EMSCRIPTEN__
      output[0] = '1';
      strcpy(&output[1], debugOutput);
      databack(output);
#else
      printf("%s", debugOutput);
#endif
    }
#endif
    // Step 2: LLL-reduce the (r+s)*(r+s) matrix M (of rank r+s).
    integralLLL(nbrVectors + nbrRequiredTraces);

    // Copy matrix M to basisStar.
    for (nbrRow = 0; nbrRow < nbrVectors + nbrRequiredTraces; nbrRow++)
    {
      for (nbrCol = 0; nbrCol < nbrVectors + nbrRequiredTraces; nbrCol++)
      {
        CopyBigInt((BigInteger*)&basisStar[nbrRow][nbrCol],
          (BigInteger*)&basis[nbrRow][nbrCol]);
      }
    }
#if DEBUG_VANHOEIJ
    strcpy(ptrDebugOutput, LF "Matrix M after LLL: " LF);
    ptrDebugOutput += strlen(ptrDebugOutput);
    for (nbrRow = 0; nbrRow < nbrVectors + nbrRequiredTraces; nbrRow++)
    {
      for (nbrCol = 0; nbrCol < nbrVectors + nbrRequiredTraces; nbrCol++)
      {
        BigInteger2Dec((BigInteger*)&basisStar[nbrRow][nbrCol], ptrDebugOutput, 0);
        ptrDebugOutput += strlen(ptrDebugOutput);
        *ptrDebugOutput++ = ' ';
      }
      *ptrDebugOutput++ = ';';
    }
    *ptrDebugOutput++ = '.';
    *ptrDebugOutput = 0;
#endif
    // Step 3: Replace the upper r*(r+s) submatrix L of M by BL*L;
    // Then set BL <- BL*L of dimension n0*(r+s)
    // M is now of dimension (n0+s)*(r+s) of rank r+s.
    // First move down the rows of M.
    for (nbrRow = nbrVectors + nbrRequiredTraces - 1; nbrRow >= nbrVectors; nbrRow--)
    {
      MiniBigInteger* ptrMBISrc = &basisStar[nbrRow][0];
      MiniBigInteger* ptrMBIDest = &basisStar[nbrRow - nbrVectors + nbrFactors][0];
      for (nbrCol = 0; nbrCol < nbrVectors + nbrRequiredTraces; nbrCol++)
      {
        CopyBigInt((BigInteger*)ptrMBIDest++, (BigInteger*)ptrMBISrc++);
      }
    }
    // Copy L to array lambda.
    for (; nbrRow >= 0; nbrRow--)
    {
      MiniBigInteger* ptrMBISrc = &basisStar[nbrRow][0];
      MiniBigInteger* ptrMBIDest = &lambda[nbrRow][0];
      for (nbrCol = 0; nbrCol < nbrVectors + nbrRequiredTraces; nbrCol++)
      {
        CopyBigInt((BigInteger*)ptrMBIDest++, (BigInteger*)ptrMBISrc++);
      }
    }
    // Multiply matrixBL (nbrFactors rows and nbrVectors columns) by
    // L (nbrVectors rows and nbrVectors + nbrRequiredTraces columns).
    // Send the product of dimension nbrFactors by nbrVectors + nbrRequiredTraces
    // to top of array M.
    for (nbrRow = 0; nbrRow < nbrFactors; nbrRow++)
    {
      for (nbrCol = 0; nbrCol < nbrVectors + nbrRequiredTraces; nbrCol++)
      {
        intToBigInteger(&operand1, 0);     // Initialize sum.
        for (j = 0; j < nbrVectors; j++)
        {
          BigIntMultiply((BigInteger*)&matrixBL[nbrRow][j], (BigInteger*)&lambda[j][nbrCol], &operand2);
          BigIntAdd(&operand1, &operand2, &operand1);
        }
        CopyBigInt((BigInteger*)&basisStar[nbrRow][nbrCol], &operand1);
      }
    }

    // Copy BL*L to array matrixBL.
    for (nbrRow = 0; nbrRow < nbrFactors; nbrRow++)
    {
      MiniBigInteger* ptrMBISrc = &basisStar[nbrRow][0];
      MiniBigInteger* ptrMBIDest = &matrixBL[nbrRow][0];
      for (nbrCol = 0; nbrCol < nbrVectors; nbrCol++)
      {
        CopyBigInt((BigInteger*)ptrMBIDest++, (BigInteger*)ptrMBISrc++);
      }
    }

#if DEBUG_VANHOEIJ
    strcpy(ptrDebugOutput, LF "Matrix M before Gram-Schmidt" LF);
    ptrDebugOutput += strlen(ptrDebugOutput);
    for (nbrRow = 0; nbrRow < nbrVectors + nbrRequiredTraces; nbrRow++)
    {
      for (nbrCol = 0; nbrCol < nbrVectors + nbrRequiredTraces; nbrCol++)
      {
        BigInteger2Dec((BigInteger*)&basisStar[nbrRow][nbrCol], ptrDebugOutput, 0);
        ptrDebugOutput += strlen(ptrDebugOutput);
        *ptrDebugOutput++ = ' ';
      }
      *ptrDebugOutput++ = ';';
    }
    *ptrDebugOutput++ = '.';
    *ptrDebugOutput = 0;
#endif
    // Step 4: Perform Gram-Schmidt orthogonalization on M.
    // Let c_i, 1 <= i <= r+s be the norm of the ith 
    // orthogonal vector.
    GramSchmidtOrthogonalization(nbrFactors + nbrRequiredTraces, nbrVectors + nbrRequiredTraces);

#if DEBUG_VANHOEIJ
    strcpy(ptrDebugOutput, LF "Norms after Gram-Schmidt: ");
    ptrDebugOutput += strlen(ptrDebugOutput);
    for (r1 = 2; r1 <= nbrVectors + nbrRequiredTraces; r1++)
    {
      BigIntDivide((BigInteger*)&lambda[r1 - 1][r1 - 1], (BigInteger*)&lambda[r1 - 2][r1 - 2], &operand1);
      BigInteger2Dec(&operand1, ptrDebugOutput, 0);
      ptrDebugOutput += strlen(ptrDebugOutput);
      *ptrDebugOutput++ = ' ';
    }
    *ptrDebugOutput++ = '.';
    *ptrDebugOutput = 0;
#endif
    // Step 5: Let r' the largest value such that all c_i of
    // index larger than r' are of norm greater than
    // sqrt(c^2*n_0 + s*n_0^2/4)

    // Compute square of formula.
    squareFormula = nbrFactors * C * C + nbrRequiredTraces * nbrFactors * nbrFactors / 4;
    for (r1 = nbrVectors + nbrRequiredTraces; r1 > 1; r1--)
    {
      // The norm of B*[r1] is lambda_{r1, r1} / lambda_{r1-1, r1-1}.
      BigIntDivide((BigInteger*)&lambda[r1 - 1][r1 - 1], (BigInteger*)&lambda[r1-2][r1-2], &operand1);
      if (operand1.nbrLimbs == 1 && operand1.limbs[0].x <= squareFormula)
      {
        break;        // r' was found.
      }
    }
#if DEBUG_VANHOEIJ
    sprintf(ptrDebugOutput, LF "squareFormula = %d", squareFormula);
    ptrDebugOutput += strlen(ptrDebugOutput);
#endif

    // Step 6: If r' = 1, return "irreducible".
    if (r1 == 1)
    {                 // Polynomial is irreducible.
#if DEBUG_VANHOEIJ
      sprintf(ptrDebugOutput, LF "r' = 1 -> irreducible." LF);
#ifdef __EMSCRIPTEN__
      output[0] = '1';
      strcpy(&output[1], debugOutput);
      databack(output);
#else
      printf("%s", debugOutput);
#endif
#endif
      return;
    }

#ifdef __EMSCRIPTEN__
    {
      char outputInfo[1000];
      char* ptrOutput = outputInfo;
      if (lang)
      {
        strcpy(ptrOutput, "1<p>Paso 7</p>");
      }
      else
      {
        strcpy(ptrOutput, "1<p>Step 7</p>");
      }
      databack(outputInfo);
    }
#endif
    // Step 7: BL <- BL / C. Now BL has dimension 
    // nbrFactors rows by r'
    for (nbrRow = 0; nbrRow < nbrFactors; nbrRow++)
    {
      MiniBigInteger* ptrMBIDest = &matrixBL[nbrRow][0];
      for (nbrCol = 0; nbrCol < r1; nbrCol++)
      {
        subtractdivide((BigInteger*)ptrMBIDest++, 0, C);
      }
    }

#if DEBUG_VANHOEIJ
    strcpy(ptrDebugOutput, LF "Matrix BL before Gauss: ");
    ptrDebugOutput += strlen(ptrDebugOutput);
    for (nbrRow = 0; nbrRow < nbrFactors; nbrRow++)
    {
      for (nbrCol = 0; nbrCol < r1; nbrCol++)
      {
        BigInteger2Dec((BigInteger*)&basisStar[nbrRow][nbrCol], ptrDebugOutput, 0);
        ptrDebugOutput += strlen(ptrDebugOutput);
        *ptrDebugOutput++ = ' ';
      }
      *ptrDebugOutput++ = ';';
    }
    *ptrDebugOutput++ = '.';
    *ptrDebugOutput = 0;
#endif
#ifdef __EMSCRIPTEN__
    {
      char outputInfo[1000];
      char* ptrOutput = outputInfo;
      if (lang)
      {
        strcpy(ptrOutput, "1<p>Paso 8</p>");
      }
      else
      {
        strcpy(ptrOutput, "1<p>Step 8</p>");
      }
      databack(outputInfo);
    }
#endif
    // Step 8: Put BL in Gauss-Jordan form. If possibly
    // after division by a constant factor, each column
    // contains only 0 and 1, and each row contains 
    // exactly one 1, we might have a valid factorization,
    // done <- true. Otherwise, increase s, replace n
    // by r and r by r'
    rank = gauss(nbrFactors, r1);  // Dimension of matrix BL.
#if DEBUG_VANHOEIJ
    sprintf(ptrDebugOutput, LF "squareFormula = %d, rank = %d" LF, squareFormula, rank);
    ptrDebugOutput += strlen(ptrDebugOutput);
#endif
    nbrRequiredTraces++;
    nbrVectors = r1;
#if 0
    if (++nbrStepsDone == 10)
    {
#ifdef __EMSCRIPTEN__
      output[0] = '1';
      strcpy(&output[1], debugOutput);
      databack(output);
#else
      printf("%s", debugOutput);
#endif
    }
#endif
    if (rank < 0)
    {            // Cannot convert to row echelon form. Try again with bigger matrix.
      continue;
    }
#if DEBUG_VANHOEIJ
    strcpy(ptrDebugOutput, LF "Matrix R after Gauss: ");
    ptrDebugOutput += strlen(ptrDebugOutput);
    for (nbrRow = 0; nbrRow < nbrFactors; nbrRow++)
    {
      for (nbrCol = 0; nbrCol < r1; nbrCol++)
      {
        BigInteger2Dec((BigInteger*)&lambda[nbrRow][nbrCol], ptrDebugOutput, 0);
        ptrDebugOutput += strlen(ptrDebugOutput);
        *ptrDebugOutput++ = ' ';
      }
      *ptrDebugOutput++ = ';';
    }
    *ptrDebugOutput++ = '.';
    *ptrDebugOutput = 0;
#endif
    // There must be only one number different from zero in each column.
    for (nbrCol = 0; nbrCol < nbrFactors; nbrCol++)
    {
      int differentFromZero = 0;
      for (nbrRow = 0; nbrRow < nbrVectors; nbrRow++)
      {
        if (!BigIntIsZero((BigInteger *)&lambda[nbrRow][nbrCol]))
        {
          differentFromZero++;
        }
      }
      if (differentFromZero != 1)
      {                 // Column should have only one 1.
        break;
      }
    }
    if (nbrCol < nbrFactors)
    {                   // Invalid matrix. Try again with bigger matrix.
#if DEBUG_VANHOEIJ
      sprintf(ptrDebugOutput, "nbrCol = %d, nbrFactors = %d. They should have been equal." LF, nbrCol, nbrFactors);
      ptrDebugOutput += strlen(ptrDebugOutput);
#endif
      continue;
    }
    CopyBigInt(&halfPowerMod, &powerMod);
    subtractdivide(&halfPowerMod, -1, 2); // halfPowerMod <- (powerMod+1)/2
    newNbrFactors = nbrFactors;
    // Combine factors.
    GetMontgomeryParms(powerMod.nbrLimbs);
    if (AttemptToFactor(nbrVectors, nbrFactors, &newNbrFactors))
    {
#if DEBUG_VANHOEIJ
      sprintf(ptrDebugOutput, "nbrVector == nbrVectors" LF);
#ifdef __EMSCRIPTEN__
      output[0] = '1';
      strcpy(&output[1], debugOutput);
      databack(output);
#else
      printf("%s", debugOutput);
#endif
#endif
      return;               // Factorization found.
    }
  }
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
    IntArray2BigInteger(ptrSrc, &operand1);
    ptrSrc += 1 + *ptrSrc;
    // If operand1 >= halfPowerMod, subtract powerMod.
    BigIntSubt(&operand1, &halfPowerMod, &operand2);
    if (operand2.sign == SIGN_POSITIVE)
    {
      BigIntSubt(&operand1, &powerMod, &operand1);
    }
    BigInteger2IntArray(ptrDest, &operand1);
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
  for (currentDegree = 0; currentDegree <= degreePoly+1; currentDegree++)
  {
    *ptrIndex++ = index;
    index += numLimbs(ptrFactor + index) + 1;
  }
  for (pstFactorInfoInteger = factorInfoInteger;
    pstFactorInfoInteger->ptrPolyLifted != NULL;
    pstFactorInfoInteger++)
  {      // Polynomial factors are sorted by degree and then
         // by coefficients.
    if (degreePoly > pstFactorInfoInteger->degree)
    {    // New factor degree is greater than stored factor degree.
      continue;  // Check next stored factor.
    }
    if (degreePoly < pstFactorInfoInteger->degree)
    {    // New factor degree is less than stored factor degree.
      break;     // Exit loop.
    }
    // Degree of the new factor is the same as the one already stored.
    ptrOldFactor = pstFactorInfoInteger->ptrPolyLifted;
    ptrIndex = &indexOldFactor[0];
    index = 0;
    for (currentDegree = 0; currentDegree <= degreePoly; currentDegree++)
    {   // Loop that gets the pointers to coefficients of stored factor.
      *ptrIndex++ = index;
      index += numLimbs(ptrOldFactor + index) + 1;
    }
    // Compare coefficients starting from higher degrees.
    for (currentDegree = degreePoly; currentDegree >= 0; currentDegree--)
    {
      UncompressBigIntegerB(ptrFactor + indexNewFactor[currentDegree], &operand1);
      UncompressBigIntegerB(ptrOldFactor + indexOldFactor[currentDegree], &operand2);
      BigIntSubt(&operand1, &operand2, &operand1);
      if (!BigIntIsZero(&operand1))
      {       // Polynomials are different.
        break;
      }
    }
    if (currentDegree < 0)
    {         // Both polynomials are the same.
      pstFactorInfoInteger->multiplicity++;
#if DEBUG_VANHOEIJ
#ifndef __EMSCRIPTEN__
      printf("%s", debugOutput);
#endif
#endif
      return;
    }
    if (operand1.sign == SIGN_NEGATIVE)
    {   // Coefficient of new factor is less than
        // coefficient of stored factor. Exit loop.
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
  ptrFactorInteger += indexNewFactor[degreePoly + 1];
}

// Input: values = degree, coefficient degree 0, coefficient degree 1, etc.
// Output: factorInfo = structure that holds the factors.
int FactorPolyOverIntegers(void)
{
  int degreePolyToFactor = values[0];
  int degree1, degree2;
  int primeRecord = 0;
  int expon, maxDegreeFactor;
  int degreeGcdMod;
  int* ptrSrc, * ptrDest;
  int attemptNbr;
  int currentDegree;
  int polXprocessed = FALSE;
  int* ptrFactorIntegerBak;
  struct sFactorInfo* pstFactorInfoOrig, * pstFactorInfoRecord;
  struct sFactorInfo* pstFactorInfoInteger = factorInfoInteger;
  ptrFactorInteger = polyInteger;
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
    BigInteger2IntArray(ptrDest, &operand2);
    ptrSrc += 1 + numLimbs(ptrSrc);
    ptrDest += 1 + numLimbs(ptrDest);
  }
  while (polyToFactor[0] != 0)
  {    // At least degree 1.
       // The trailing coefficient of factors must divide the product of the trailing and leading
       // coefficients of the original polynomial.
    int nbrFactorsRecord;
    int prime;
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
    // Up to 5 different prime moduli are tested to minimize the number
    // of factors because the Van Hoeij algorithm speed depends on
    // this number. If the number of factors is less than 10,
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
        ptrSrc = &polyNonRepeatedFactors[1];
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
      computePower(1);
      exponentMod = 1;
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
        pstFactorInfoOrig = factorInfo;
        pstFactorInfoRecord = factorInfoRecord;
        ptrPolyLiftedRecord = polyLiftedRecord;
        nbrFactorsRecord = 0;
        for (factorNbr = 0; factorNbr < MAX_DEGREE; factorNbr++)
        {
          if (pstFactorInfoOrig->ptr == NULL)
          {    // No more factors.
            break;
          }
          *pstFactorInfoRecord = *pstFactorInfoOrig;
          pstFactorInfoRecord->ptr = ptrPolyLiftedRecord;
          nbrFactorsRecord += pstFactorInfoOrig->multiplicity;
          ptrPolyLiftedRecord = CopyPolynomialFixedCoeffSize(ptrPolyLiftedRecord,
            pstFactorInfoOrig->ptr,
            pstFactorInfoOrig->degree - 1, primeMod.nbrLimbs + 1);
          *ptrPolyLiftedRecord++ = 1;    // Leading coefficient should be 1.
          *ptrPolyLiftedRecord++ = 1;
          pstFactorInfoOrig++;
          pstFactorInfoRecord++;
        }
        if (factorNbr < MAX_DEGREE)
        {
          pstFactorInfoRecord->ptr = NULL;
        }
        if (nbrFactors < 10)
        {
          break;    // Small enough number of factors. Go out of loop.
        }
      }
    }
    prime = primeRecord;
    nbrFactorsFound = nbrFactorsRecord;
    intToBigInteger(&primeMod, prime);
    vanHoeij(prime, nbrFactorsRecord);
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

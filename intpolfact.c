//
// This file is part of Alpertron Calculators.
//
// Copyright 2019-2021 Dario Alejandro Alpern
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
#include "linkedbignbr.h"
#include "highlevel.h"
#include "polynomial.h"
#include "showtime.h"
#include "rootseq.h"
#include "musl.h"

#define MAX_MATRIX_SIZE  200

char* ptrOutput2;
#if DEBUG_VANHOEIJ
char debugOutput[20000000];
char* ptrDebugOutput = debugOutput;
static int z;
#endif

#ifdef __EMSCRIPTEN__
#define LF "<br>"
#else
#define LF "\r\n"
#endif

struct linkedBigInt* basis[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE];
struct linkedBigInt* basisStar[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE];
struct linkedBigInt* lambda[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE];
struct linkedBigInt* matrixBL[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE];
struct linkedBigInt* detProdB[MAX_MATRIX_SIZE];
struct linkedBigInt* traces[MAX_MATRIX_SIZE * 10];
struct linkedBigInt* ptrCoeffs[MAX_MATRIX_SIZE * 10];
BigInteger powerBoundA;
BigInteger powerExtraBits;
static int* ptrFactorInteger;

static BigInteger contentPolyToFactor;
static BigInteger halfPowerMod;
static int origPolyToFactor[1000000];
static int polyToFactor[1000000];
static int factorX[] = { 1, 0, 1, 1 }; // Polynomial is x.
static BigInteger bound;
static BigInteger trailingCoeff;
static BigInteger leadingCoeff;
int polyNonRepeatedFactors[1000000];
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
static void GenerateIntegerPolynomial(const int* polyMod, int* polyInt, int degreePoly);
static void InsertIntegerPolynomialFactor(int* ptrFactor, int degreePoly);
static int numberLLL;

// Generate row echelon form from matrixBL.
// The output will be located in matrix lambda.
static int gauss(int nbrCols, int nbrRows)
{
  int l;
  int row;
  int col;

  for (row = 0; row < nbrRows; row++)
  {
    for (col = 0; col < nbrCols; col++)
    {
      getBigIntegerFromLinked(matrixBL[col][row], &tmp5);
      setLinkedBigInteger(&lambda[row][col], &tmp5);
    }
  }

  l = 0;
  for (int k = 0; (k < nbrCols) && (l < nbrRows); k++)
  {
    int pos = -1;
    /* Look for a pivot under the diagonal. */
    row = l;
    while ((row < nbrRows) && (linkedBigIntIsZero(lambda[row][k])))
    {
      row++;
    }

    if (row == nbrRows)
      /* No pivot found; try to find a unique 1 above */
    {
      row = 0;
      while ((row < l) && linkedBigIntIsZero(lambda[row][k]))
      {
        row++;
      }
      if ((row == l) || !linkedBigIntIsOne(lambda[row][k]))
      {
        return -1;           // No pivot found, go out.
      }
      row++;
      for (; row < l; row++)
      {
        if (!linkedBigIntIsZero(lambda[row][k]))
        {
          return -1;
        }
      }
      continue;
    }

    // We have found a non-zero element on the k-th column
    for (; (row < nbrRows) && (pos == -1); row++)
    {
      if (linkedBigIntIsMinusOne(lambda[row][k]))   // Value is -1.
      {   // invert all elements on i-th row
        for (col = k; col < nbrCols; col++)
        {
          linkedBigIntChSign(lambda[row][col]);
        }
        pos = row;
      }
      if (linkedBigIntIsOne(lambda[row][k]))
      {
        pos = row;
      }
    }

    if (pos != -1)
    { /* necessarily M[pos][k]=1 */
      for (col = 0; col < nbrCols; col++)
      {
        struct linkedBigInt* pstLinkedBigInt;
        pstLinkedBigInt = lambda[pos][col];
        lambda[pos][col] = lambda[l][col];
        lambda[l][col] = pstLinkedBigInt;
      }
      for (row = 0; row < nbrRows; row++)
      {
        // M[i] = M[i] - M[l]
        getBigIntegerFromLinked(lambda[row][k], &tmp1);
        if ((row != l) && !BigIntIsZero(&tmp1))
        {
          intToLinkedBigInt(&lambda[row][k], 0);
          for (col = k + 1; col < nbrCols; col++)
          {
            // *x = *x + (*y)*t1
            getBigIntegerFromLinked(lambda[l][col], &tmp5);
            (void)BigIntMultiply(&tmp5, &tmp1, &tmp0);
            getBigIntegerFromLinked(lambda[row][col], &tmp5);
            BigIntSubt(&tmp5, &tmp0, &tmp5);
            setLinkedBigInteger(&lambda[row][col], &tmp5);
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
  int colJ;
  int k;
  for (int colI = 0; colI < nbrCols; colI++)
  {
    if (colI == 0)
    {             // U_0 <- 1
      intToLinkedBigInt(&detProdB[0], 1);
    }
    else
    {             // U_i <- lambda_{i-1, i-1}
      getBigIntegerFromLinked(lambda[colI - 1][colI - 1], &tmp5);
      setLinkedBigInteger(&detProdB[colI], &tmp5);

                  // U_{i-1} <- -lambda_{i, i-1}
      getBigIntegerFromLinked(lambda[colI][colI - 1], &tmp5);
      BigIntChSign(&tmp5);
      setLinkedBigInteger(&detProdB[colI - 1], &tmp5);
    }
    for (colJ = colI - 2; colJ >= 0; colJ--)
    {             // U_j <- 0
      intToBigInteger(&tmp1, 0);
      for (k = colJ + 1; k <= colI; k++)
      {           // U_j <- U_j + lambda_{k,j} * U_k
        getBigIntegerFromLinked(lambda[k][colJ], &tmp4);
        getBigIntegerFromLinked(detProdB[k], &tmp5);
        (void)BigIntMultiply(&tmp4, &tmp5, &tmp2);
        BigIntSubt(&tmp1, &tmp2, &tmp1);
      }
      // U_j <- U_j / lambda_{j,j}
      getBigIntegerFromLinked(lambda[colJ][colJ], &tmp4);
      (void)BigIntDivide(&tmp1, &tmp4, &tmp5);
      setLinkedBigInteger(&detProdB[colJ], &tmp5);
    }
    for (colJ = colI; colJ < nbrCols; colJ++)
    {           // lambda_{j,i} <- 0. Use tmp1 for lambda_{j,i}
      intToBigInteger(&tmp1, 0);
      for (k = 0; k <= colI; k++)
      {         // tmp2 <- scalar product b_j * b_k 
        intToBigInteger(&tmp2, 0);
        for (int row = 0; row < nbrRows; row++)
        {
          getBigIntegerFromLinked(basisStar[row][colJ], &tmp4);
          getBigIntegerFromLinked(basisStar[row][k], &tmp5);
          (void)BigIntMultiply(&tmp4, &tmp5, &tmp3);
          BigIntAdd(&tmp2, &tmp3, &tmp2);
        }
        // lambda_{j,i} <- lambda_{j,i} + tmp2 * U_k
        getBigIntegerFromLinked(detProdB[k], &tmp5);
        (void)BigIntMultiply(&tmp2, &tmp5, &tmp2);
        BigIntAdd(&tmp1, &tmp2, &tmp1);
      }
      setLinkedBigInteger(&lambda[colJ][colI], &tmp1);
    }
  }
}

static void PerformREDI(int k, int l, int size)
{
  // If |2 lambda_{k, l}| <= d_l, go out.
  getBigIntegerFromLinked(lambda[k][l], &tmp4);
  multint(&tmp0, &tmp4, 2);          // tmp0 <- 2 lambda_{k, l}
  tmp0.sign = SIGN_POSITIVE;         // tmp0 <- |2 lambda_{k, l}|
  getBigIntegerFromLinked(detProdB[l], &tmp5);
  BigIntSubt(&tmp5, &tmp0, &tmp1);   // tmp1 <- d_l - |2 lambda_{k, l}|
  if (tmp1.sign == SIGN_POSITIVE)
  {
    return;                          // Go out if d_l >= |2 lambda_{k, l}|
  }
  // Compute q = nearest integer to lambda_{k, l} / d_l.
  // d_l is always positive.
  // If lambda is positive, compute (2 lambda_{k, l} + d_l) / (2 d_l)
  // If lambda is negative, compute (2 lambda_{k, l} - d_l) / (2 d_l)
  // the division, change sign of quotient.
  multint(&tmp0, &tmp4, 2);          // tmp0 <- 2 lambda_{k, l}
  if (tmp4.sign == SIGN_POSITIVE)
  {
    BigIntAdd(&tmp0, &tmp5, &tmp0);   // tmp0 <- 2 lambda_{k, l} + d_l
  }
  else
  {
    BigIntSubt(&tmp0, &tmp5, &tmp0);  // tmp0 <- 2 lambda_{k, l} - d_l
  }
  multint(&tmp1, &tmp5, 2);           // tmp1 <- 2 d_l
  (void)BigIntDivide(&tmp0, &tmp1, &tmp2); // tmp2 <- q.
  for (int row = 0; row < size; row++)
  {  // Loop that computes b_k <- b_k - q*b_l
    getBigIntegerFromLinked(basis[row][l], &tmp5);
    (void)BigIntMultiply(&tmp2, &tmp5, &tmp0);  // tmp0 <- q*b_l
    getBigIntegerFromLinked(basis[row][k], &tmp5);
    BigIntSubt(&tmp5, &tmp0, &tmp5);
    setLinkedBigInteger(&basis[row][k], &tmp5);
  }
    // lambda_{k, l} <- lambda_{k, l} - q*d_l
  getBigIntegerFromLinked(detProdB[l], &tmp5);
  (void)BigIntMultiply(&tmp2, &tmp5, &tmp0);    // tmp0 <- q*d_l
  getBigIntegerFromLinked(lambda[k][l], &tmp5);
  BigIntSubt(&tmp5, &tmp0, &tmp5);
  setLinkedBigInteger(&lambda[k][l], &tmp5);
  for (int i = 0; i < l; i++)
  { // lambda_{k, i} <- lambda_{k, i} - q*lambda_{l, i}
    getBigIntegerFromLinked(lambda[l][i], &tmp5);
    (void)BigIntMultiply(&tmp2, &tmp5, &tmp0);   // tmp0 <- q*lambda_{l, i}
    getBigIntegerFromLinked(lambda[k][i], &tmp5);
    BigIntSubt(&tmp5, &tmp0, &tmp5);
    setLinkedBigInteger(&lambda[k][i], &tmp5);
  }
}

static void PerformSWAPI(int k, int kMax, int size)
{                   // On entry: k >= 1
  struct linkedBigInt* pstLinkedBigInt;
  // Exchange b_k with b_{k-1}
  for (int row = 0; row < size; row++)
  {
    pstLinkedBigInt = basis[row][k];
    basis[row][k] = basis[row][k - 1];
    basis[row][k - 1] = pstLinkedBigInt;
  }
  for (int j = 0; j < k - 1; j++)
  {  // Exchange lambda_{k, j} with lambda_{k-1, j}
    pstLinkedBigInt = lambda[k][j];
    lambda[k][j] = lambda[k - 1][j];
    lambda[k - 1][j] = pstLinkedBigInt;
  }
    // Set lambda <- lambda_{k, k-1}
  getBigIntegerFromLinked(lambda[k][k - 1], &tmp0);    // tmp0 <- lambda.
    // Set B <- (d_{k-2}*d_k + lambda^2)/d_{k-1}
    // d_{k-2}*d_k + lambda^2 is already in tmp3.
  getBigIntegerFromLinked(detProdB[k - 1], &tmp5);
  (void)BigIntDivide(&tmp3, &tmp5, &tmp1);                   // tmp1 <- B
#if DEBUG_VANHOEIJ
  if ((size == 16) && (z < 1000))
  {
    z++;
    (void)strcpy(ptrDebugOutput, "k = ");
    ptrDebugOutput += strlen(ptrDebugOutput);
    int2dec(&ptrDebugOutput, k);
    (void)strcpy(ptrDebugOutput, ", z = ");
    ptrDebugOutput += strlen(ptrDebugOutput);
    int2dec(&ptrDebugOutput, z);
    *ptrDebugOutput++ = '\r';
    *ptrDebugOutput++ = '\n';
    (void)strcpy(ptrDebugOutput, "tmp3 = ");
    ptrDebugOutput += strlen(ptrDebugOutput);
    BigInteger2Dec(&tmp3, ptrDebugOutput, 0);
    ptrDebugOutput += strlen(ptrDebugOutput);
    *ptrDebugOutput++ = '\r';
    *ptrDebugOutput++ = '\n';
    (void)strcpy(ptrDebugOutput, "detProdB[k-1] = ");
    ptrDebugOutput += strlen(ptrDebugOutput);
    getBigIntegerFromLinked(detProdB[k - 1], &tmp5);
    BigInteger2Dec(&tmp5, ptrDebugOutput, 0);
    ptrDebugOutput += strlen(ptrDebugOutput);
    *ptrDebugOutput++ = '\r';
    *ptrDebugOutput++ = '\n';
    (void)strcpy(ptrDebugOutput, "B = ");
    ptrDebugOutput += strlen(ptrDebugOutput);
    BigInteger2Dec(&tmp1, ptrDebugOutput, 0);
    ptrDebugOutput += strlen(ptrDebugOutput);
    *ptrDebugOutput++ = '\r';
    *ptrDebugOutput++ = '\n';
  }
#endif
  for (int i = k+1; i <= kMax; i++)
  {
    // t <- lambda_{i, k}
    getBigIntegerFromLinked(lambda[i][k], &tmp2);     // tmp2 <- t
    // lambda_{i, k} <- (d_k*lambda_{i, k-1} - lambda * t)/d_{k-1}
    getBigIntegerFromLinked(detProdB[k], &tmp4);
    getBigIntegerFromLinked(lambda[i][k - 1], &tmp5);
    (void)BigIntMultiply(&tmp4, &tmp5, &tmp3);
    (void)BigIntMultiply(&tmp0, &tmp2, &tmp4);       // tmp4 <- lambda * t
    BigIntSubt(&tmp3, &tmp4, &tmp3);
    getBigIntegerFromLinked(detProdB[k - 1], &tmp5);
    (void)BigIntDivide(&tmp3, &tmp5, &tmp4);
    setLinkedBigInteger(&lambda[i][k], &tmp4);
    // lambda_{i, k-1} <- (B*t + lambda * lambda_{i, k})/d_k
    (void)BigIntMultiply(&tmp1, &tmp2, &tmp3);       // tmp3 <- B * t
    (void)BigIntMultiply(&tmp0, &tmp4, &tmp4);
    BigIntAdd(&tmp3, &tmp4, &tmp3);
    getBigIntegerFromLinked(detProdB[k], &tmp4);
    (void)BigIntDivide(&tmp3, &tmp4, &tmp5);
    setLinkedBigInteger(&lambda[i][k - 1], &tmp5);
  }
  setLinkedBigInteger(&detProdB[k - 1], &tmp1);     // d_{k-1} <- B.
}

// Use algorithm 2.6.7 of Henri Cohen's book
// A Course in Computational Algebraic Number Theory.
// Here the indexes start at zero instead of one as in the book.
// B_0, B_1, ... are the different columns of matrix basis.
// d_k = 1 when k < 0.
void integralLLL(int size)
{
  int colK;
  int colKMax;
  int row;
#ifdef __EMSCRIPTEN__
  char outputInfo[1000];
  char* ptrOutput = outputInfo;
  if (lang)
  {
    (void)strcpy(ptrOutput, "1<p>Calculando LLL número ");
    ptrOutput += strlen(ptrOutput);
    int2dec(&ptrOutput, ++numberLLL);
    (void)strcpy(ptrOutput, " en matriz de");
  }
  else
  {
    (void)strcpy(ptrOutput, "1<p>Computing LLL #");
    ptrOutput += strlen(ptrOutput);
    int2dec(&ptrOutput, ++numberLLL);
    (void)strcpy(ptrOutput, " in matrix of");
  }
  ptrOutput += strlen(ptrOutput);
  *ptrOutput = ' ';
  ptrOutput++;
  int2dec(&ptrOutput, size);
  (void)strcpy(ptrOutput, " &times; ");
  ptrOutput += strlen(ptrOutput);
  int2dec(&ptrOutput, size);
  (void)strcpy(ptrOutput, ".</p>");
  ptrOutput += strlen(ptrOutput);
  showElapsedTimeSec(&ptrOutput);
  databack(outputInfo);
#endif
  colK = 1;
  colKMax = 0;  // Vectors b_0 to b_kMax are LLL-reduced.
  // Set d_0 to scalar product B_0 * B_0
  // Array detProdB holds the array d in the algorithm.
  intToLinkedBigInt(&detProdB[0], 0);
  for (row = 0; row < size; row++)
  {       // Loop that generates the scalar product B_0 * B_0.
    getBigIntegerFromLinked(basis[row][0], &tmp5);
    (void)BigIntMultiply(&tmp5, &tmp5, &tmp0);
    getBigIntegerFromLinked(detProdB[0], &tmp5);
    BigIntAdd(&tmp5, &tmp0, &tmp5);
    setLinkedBigInteger(&detProdB[0], &tmp5);     // d_{k-1} <- B.
  }
  do
  {    
    if (colK > colKMax)
    { // Perform incremental Gram-Schmidt
      colKMax = colK;
      for (int colJ = 0; colJ <= colK; colJ++)
      {
        // Compute u <- scalar product of b_k and b_j (use tmp2 for u)
        intToBigInteger(&tmp2, 0);
        for (row = 0; row < size; row++)
        {
          getBigIntegerFromLinked(basis[row][colK], &tmp4);
          getBigIntegerFromLinked(basis[row][colJ], &tmp5);
          (void)BigIntMultiply(&tmp4, &tmp5, &tmp0);
          BigIntAdd(&tmp2, &tmp0, &tmp2);
        }
        for (int colI = 0; colI < colJ; colI++)
        {    // Set u <- (d_i * u - lambda_{k, i} * lambda_{j, i}) / d_{i-1}
          getBigIntegerFromLinked(detProdB[colI], &tmp5);
          (void)BigIntMultiply(&tmp5, &tmp2, &tmp0);  // d_i * u
          getBigIntegerFromLinked(lambda[colK][colI], &tmp4);
          getBigIntegerFromLinked(lambda[colJ][colI], &tmp5);
          (void)BigIntMultiply(&tmp4, &tmp5, &tmp1);  // lambda_{k, i} * lambda_{j, i}
          if (colI == 0)
          {
            BigIntSubt(&tmp0, &tmp1, &tmp2);           // d_{i-1} = 0
          }
          else
          {
            BigIntSubt(&tmp0, &tmp1, &tmp0);
            getBigIntegerFromLinked(detProdB[colI - 1], &tmp5);
            (void)BigIntDivide(&tmp0, &tmp5, &tmp2);
          }
        }
        if (colJ < colK)
        {      // Set lambda_{j, k} <- u
          setLinkedBigInteger(&lambda[colK][colJ], &tmp2);
        }
        else
        {      // Set d_k <- u
          setLinkedBigInteger(&detProdB[colK], &tmp2);
        }
      }
    }
#if DEBUG_VANHOEIJ
    if (size == 16)
    {
      (void)strcpy(ptrDebugOutput, "lambda: ");
      ptrDebugOutput += strlen(ptrDebugOutput);
      for (row = 0; row < size; row++)
      {
        (void)strcpy(ptrDebugOutput, "Row #");
        ptrDebugOutput += strlen(ptrDebugOutput);
        int2dec(&ptrDebugOutput, row);
        *ptrDebugOutput++ = ':';
        *ptrDebugOutput++ = ' ';
        for (colI = 0; colI < size; colI++)
        {
          getBigIntegerFromLinked(lambda[row][colI], &tmp5);
          BigInteger2Dec(&tmp5, ptrDebugOutput, 0);
          ptrDebugOutput += strlen(ptrDebugOutput);
          *ptrDebugOutput++ = ',';
          *ptrDebugOutput++ = ' ';
        }
        *(ptrDebugOutput - 2) = '\r';
        *(ptrDebugOutput - 1) = '\n';
      }
      (void)strcpy(ptrDebugOutput, "detProdB: ");
      ptrDebugOutput += strlen(ptrDebugOutput);
      for (colI = 0; colI < size; colI++)
      {
        BigInteger2Dec(&detProdB[colI], ptrDebugOutput, 0);
        ptrDebugOutput += strlen(ptrDebugOutput);
        *ptrDebugOutput++ = ',';
        *ptrDebugOutput++ = ' ';
      }
      *(ptrDebugOutput - 2) = '\r';
      *(ptrDebugOutput - 1) = '\n';
    }
#endif
    // Test LLL condition
    for (;;)
    {
      PerformREDI(colK, colK - 1, size);
      // Check whether 4 * d_k * d_{k-2} < 3 (d_{k-1})^2 - 4*(lambda_{k, k-1})^2
      // that means 4 (d_k * d_{k-2} + lambda_{k, k-1})^2) < 3 (d_{k-1})^2 
                       // Get d_k
      getBigIntegerFromLinked(detProdB[colK], &tmp0);
      if (colK > 1)
      {                // Multiply by d_{k-2} if k>=0
        getBigIntegerFromLinked(detProdB[colK - 2], &tmp5);
        (void)BigIntMultiply(&tmp0, &tmp5, &tmp0);
      }
                       // Compute lambda_{k, k-1})^2
      getBigIntegerFromLinked(lambda[colK][colK - 1], &tmp5);
      (void)BigIntMultiply(&tmp5, &tmp5, &tmp1);
                       // Compute d_k * d_{k-2} + lambda_{k, k-1})^2
      BigIntAdd(&tmp0, &tmp1, &tmp3);
      multint(&tmp0, &tmp3, 4);           // tmp0 = Left Hand Side.
                       // Compute (d_{k-1})^2 
      getBigIntegerFromLinked(detProdB[colK - 1], &tmp5);
      (void)BigIntMultiply(&tmp5, &tmp5, &tmp1);
                       // Compute 3 * (d_{k-1})^2 = Right hand side.
      multint(&tmp1, &tmp1, 3);
      BigIntSubt(&tmp0, &tmp1, &tmp0);    // tmp0 = LHS - RHS.
      if (tmp0.sign == SIGN_POSITIVE)
      {
        break;                            // Exit loop if LHS - RHS >= 0.
      }
      PerformSWAPI(colK, colKMax, size);
      if (colK > 1)                       // k <- max(1,k-1)
      {
        colK--;
      }
    }
    for (int l = colK - 2; l >= 0; l--)
    {
      PerformREDI(colK, l, size);
    }
#if DEBUG_VANHOEIJ
    if (size == 16)
    {
      (void)strcpy(ptrDebugOutput, "basis: ");
      ptrDebugOutput += strlen(ptrDebugOutput);
      for (row = 0; row < size; row++)
      {
        (void)strcpy(ptrDebugOutput, "Row #");
        ptrDebugOutput += strlen(ptrDebugOutput);
        int2dec(&ptrDebugOutput, row);
        *ptrDebugOutput++ = ':';
        *ptrDebugOutput++ = ' ';
        for (colI = 0; colI < size; colI++)
        {
          getBigIntegerFromLinked(basis[row][colI], &tmp5);
          BigInteger2Dec(&tmp5, ptrDebugOutput, 0);
          ptrDebugOutput += strlen(ptrDebugOutput);
          *ptrDebugOutput++ = ',';
          *ptrDebugOutput++ = ' ';
        }
        *(ptrDebugOutput - 2) = '\r';
        *(ptrDebugOutput - 1) = '\n';
      }
    }
#endif
  } while (++colK < size);
}

// Compute remainder such that the result is in range -divisor/2 to divisor/2.
static void BigIntSymmetricRemainder(BigInteger* dividend, BigInteger* divisor, BigInteger* result)
{
  (void)BigIntRemainder(dividend, divisor, result);
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
  const int* ptrCoeff;
  int polyDegree;
  struct linkedBigInt** ptrTrace;
  const struct sFactorInfo* pstFactorInfo = &factorInfoRecord[nbrCol];
  // Get pointers to the coefficients.
  polyDegree = pstFactorInfo->degree;          // Degree of polynomial factor
  ptrCoeff = pstFactorInfo->ptrPolyLifted;     // coefficients of polynomial factor.
  for (traceNbr = polyDegree - 1; traceNbr >= 0; traceNbr--)
  {
    operand4.nbrLimbs = *ptrCoeff;
    operand4.sign = SIGN_POSITIVE;
    (void)memcpy(operand4.limbs, ptrCoeff+1, *ptrCoeff *sizeof(int));
    (void)BigIntRemainder(&operand4, &powerMod, &operand3);
    setLinkedBigInteger(&ptrCoeffs[traceNbr], &operand3);
    ptrCoeff += 1 + numLimbs(ptrCoeff);
  }
  // Compute coefficients of P.
  intToBigInteger(&operand1, 1);
#if DEBUG_VANHOEIJ
  (void)strcpy(ptrDebugOutput, "Coefficients: ");
  ptrDebugOutput += strlen(ptrDebugOutput);
#endif
  for (traceNbr = 0; traceNbr < polyDegree; traceNbr++)
  {
    getBigIntegerFromLinked(ptrCoeffs[traceNbr], &operand3);
#if DEBUG_VANHOEIJ
    BigInteger2Dec(&operand3, ptrDebugOutput, 0);
    ptrDebugOutput += strlen(ptrDebugOutput);
    *ptrDebugOutput++ = ',';
    *ptrDebugOutput++ = ' ';
#endif
    (void)BigIntMultiply(&operand3, &operand1, &operand3);
    setLinkedBigInteger(&ptrCoeffs[traceNbr], &operand3);
  }
#if DEBUG_VANHOEIJ
  (void)strcpy(ptrDebugOutput, LF);
  ptrDebugOutput += strlen(ptrDebugOutput);
#endif
  // Store traces of polynomial P in matrix traces.
  ptrTrace = &traces[0];
  intToLinkedBigInt(ptrTrace, polyDegree);
  intToBigInteger(&operand1, polyDegree);
  (void)BigIntRemainder(&operand1, &powerMod, &tmp5);
  setLinkedBigInteger(&traces[0], &tmp5);
  for (traceNbr = 1; traceNbr < nbrTraces; traceNbr++)
  {
    ptrTrace++;
    // Initialize trace to -n*E_n.
    if (traceNbr <= polyDegree)
    {
      getBigIntegerFromLinked(ptrCoeffs[traceNbr - 1], &operand1);
      multint(&tmp5, &operand1, -traceNbr);
    }
    else
    {
      intToBigInteger(&tmp5, 0);
    }
    for (int currDegree = 1; currDegree < traceNbr; currDegree++)
    {   // Subtract - E_{n-deg}*Tr_deg(P)
      if (traceNbr - currDegree <= polyDegree)
      {
        getBigIntegerFromLinked(ptrCoeffs[traceNbr - currDegree - 1], &operand3);
        getBigIntegerFromLinked(traces[currDegree], &tmp4);
        (void)BigIntMultiply(&tmp4, &operand3, &operand1);
        BigIntSubt(&tmp5, &operand1, &tmp5);
        // Get remainder of quotient of result by p^a.
        (void)BigIntRemainder(&tmp5, &powerMod, &tmp5);
      }
    }
    setLinkedBigInteger(ptrTrace, &tmp5);
  }
  ptrTrace = &traces[0];
  for (traceNbr = 1; traceNbr < nbrTraces; traceNbr++)
  {       // Loop that divides all traces by p^b and get modulus p^(a-b).
    ptrTrace++;
    getBigIntegerFromLinked(*ptrTrace, &tmp5);
    // Get remainder of quotient of trace divided by p^b
    BigIntSymmetricRemainder(&tmp5, &powerBoundA, &operand2);
    // Subtract the remainder from the computed trace.
    BigIntSubt(&tmp5, &operand2, &tmp5);
    // Divide the result by p^b
    (void)BigIntDivide(&tmp5, &powerBoundA, &tmp5);
    // Get remainder of quotient of result divided by p^(a-b)
    BigIntSymmetricRemainder(&tmp5, &powerExtraBits, &tmp5);
    setLinkedBigInteger(ptrTrace, &tmp5);
  }
}

static bool AttemptToFactor(int nbrVectors, int nbrFactors, int *pNbrFactors)
{
  int currentDegree;
  int degreeProd;
  int degreeFactor;
  int rc;
  struct sFactorInfo* pstFactorInfo;
  int* ptrSrc;
  int* ptrDest;
  int nbrTmp[1000];
  int nbrTmp2[1000];
  int nbrTmp3[1000];
  int TestNbr0Bak;
  int dividendMod32768[2 * MAX_DEGREE + 1];
  int divisorMod32768[2 * MAX_DEGREE + 1];
  // In step 1 we check that all factors can be found. If this succeeds,
  // in step 2 we insert factors in final array of polynomial factors.
  for (int stepNbr = 1; stepNbr <= 2; stepNbr++)
  {
    int* ptrMod32768;
#if DEBUG_VANHOEIJ
    (void)strcpy(ptrDebugOutput, "stepNbr = ");
    ptrDebugOutput += strlen(ptrDebugOutput);
    int2dec(&ptrDebugOutput, stepNbr);
#endif
    NumberLength = powerMod.nbrLimbs;
    for (int nbrVector = 0; nbrVector < nbrVectors; nbrVector++)
    {
      // Test whether the product of the polynomial factors
      // defined by this vector divides the original polynomial.
      // First reduce coefficients mod powerMod and then
      // convert them to Montgomery notation.
      modulusIsZero = false;  // Perform modular operations.
      NumberLength = powerMod.nbrLimbs;
      degreeProd = 0;     // Set product to 1.
      pstFactorInfo = factorInfoRecord;
      for (int currentFactor = 0; currentFactor < nbrFactors; currentFactor++)
      {
        if (!linkedBigIntIsZero(lambda[nbrVector][currentFactor]))
        {                 // Multiply this polynomial.
          const int* ptrCoeffSrc = pstFactorInfo->ptrPolyLifted;    // Source
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
              (void)memcpy(ptrCoeffDest + 1, ptrCoeffSrc + 1, nbrLength * sizeof(int));
            }
            else
            {              // Reduce the coefficient mod powerMod.
              operand1.nbrLimbs = nbrLength;
              operand1.sign = SIGN_POSITIVE;
              (void)memcpy(operand1.limbs, ptrCoeffSrc + 1, nbrLength * sizeof(int));
              (void)BigIntRemainder(&operand1, &powerMod, &operand2);
              *ptrCoeffDest = operand2.nbrLimbs;
              (void)memcpy(ptrCoeffDest + 1, operand2.limbs, operand2.nbrLimbs * sizeof(int));
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
            (void)memcpy(ptrCoeffDest, ptrCoeffSrc, nbrLength * sizeof(int));
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
      // In the same loop get the coefficients of the divisor polynomial
      // by computing the coefficient mod 32768.
      ptrSrc = &poly5[1];
      ptrMod32768 = divisorMod32768;
      for (currentDegree = 0; currentDegree <= poly5[0]; currentDegree++)
      {
        UncompressBigIntegerB(ptrSrc, &operand1);
        BigIntSubt(&operand1, &bound, &operand2);
        if (operand2.sign == SIGN_POSITIVE)
        {
          return false;     // Coefficient is too high.
        }
        BigIntAdd(&operand1, &bound, &operand2);
        if (operand2.sign == SIGN_NEGATIVE)
        {
          return false;     // Coefficient is too low.
        }
        *ptrMod32768++ = 1;
        if (*ptrSrc >= 0)
        {
          *ptrMod32768++ = *(ptrSrc+1) & 32767;
        }
        else
        {
          *ptrMod32768++ = (-*(ptrSrc + 1)) & 32767;
        }
        ptrSrc += 1 + numLimbs(ptrSrc);
      }
      modulusIsZero = true;   // Perform integer division.
      // Multiply all coefficients by leadingCoeff and store in polyS.
      // In the same loop get the polynomial mod 32768.
      ptrMod32768 = dividendMod32768;
      polyS[0] = polyNonRepeatedFactors[0];
      ptrSrc = &polyNonRepeatedFactors[1];
      ptrDest = &polyS[1];
      for (currentDegree = 0; currentDegree <= polyS[0]; currentDegree++)
      {
        UncompressBigIntegerB(ptrSrc, &operand1);
        *ptrMod32768++ = 1;
        if (operand1.sign == SIGN_POSITIVE)
        {
          *ptrMod32768++ = operand1.limbs[0].x & 32767;
        }
        else
        {
          *ptrMod32768++ = (-operand1.limbs[0].x) & 32767;
        }
        (void)BigIntMultiply(&operand1, &leadingCoeff, &operand2);
        NumberLength = operand2.nbrLimbs;
        BigInteger2IntArray(ptrDest, &operand2);
        ptrSrc += 1 + numLimbs(ptrSrc);
        ptrDest += 1 + numLimbs(ptrDest);
      }
      if (stepNbr == 1)
      {   // Test whether the product divides the original polynomial.
          // This means that f(0) divides F(0) where f(x) is the
          // polynomial just computed and F(x) is the original
          // polynomial times the leading coefficient.
        UncompressBigIntegerB(&polyS[1], &operand1);
        UncompressBigIntegerB(&poly5[1], &operand2);
        (void)BigIntRemainder(&operand1, &operand2, &operand3);
        if (!BigIntIsZero(&operand3))
        {    // Constant term of product does not divide constant term of
             // original polynomial, so go out.
          return false;
        }
          // Check that f(1) divides F(1) and f(-1) divides F(-1).
          // Accumulate even coefficients of F in operand1 and odd
          // coefficients in operand2.
        intToBigInteger(&operand1, 0);
        intToBigInteger(&operand2, 0);
        ptrDest = &polyS[1];
        for (currentDegree = 0; currentDegree <= polyS[0]; currentDegree++)
        {
          UncompressBigIntegerB(ptrDest, &tmp1);
          ptrDest += 1 + numLimbs(ptrDest);
          if (currentDegree & 1)
          {
            BigIntAdd(&operand2, &tmp1, &operand2);
          }
          else
          {
            BigIntAdd(&operand1, &tmp1, &operand1);
          }
        }
        intToBigInteger(&operand3, 0);
        intToBigInteger(&operand4, 0);
        ptrDest = &poly5[1];
        for (currentDegree = 0; currentDegree <= poly5[0]; currentDegree++)
        {
          UncompressBigIntegerB(ptrDest, &tmp1);
          ptrDest += 1 + numLimbs(ptrDest);
          if (currentDegree & 1)
          {
            BigIntAdd(&operand4, &tmp1, &operand4);
          }
          else
          {
            BigIntAdd(&operand3, &tmp1, &operand3);
          }
        }
        BigIntAdd(&operand3, &operand4, &tmp1);     // Get f(1).
        BigIntAdd(&operand1, &operand2, &tmp2);     // Get F(1).
        if (BigIntIsZero(&tmp1))
        {
          if (!BigIntIsZero(&tmp2))
          {
            return false;   // f(1) is zero and F(1) is not. Go out.
          }
        }
        else
        {
          BigIntAdd(&operand1, &operand2, &tmp2);   // Get F(1).
          (void)BigIntRemainder(&tmp2, &tmp1, &tmp3);
          if (!BigIntIsZero(&tmp3))
          {    // f(1) does not divide F(1), so go out.
            return false;
          }
        }
        BigIntSubt(&operand3, &operand4, &tmp1);    // Get f(-1).
        BigIntSubt(&operand1, &operand2, &tmp2);    // Get F(-1).
        if (BigIntIsZero(&tmp1))
        {
          if (!BigIntIsZero(&tmp2))
          {
            return false;   // f(-1) is zero and F(-1) is not. Go out.
          }
        }
        else
        {
          (void)BigIntRemainder(&tmp2, &tmp1, &tmp3);
          if (!BigIntIsZero(&tmp3))
          {    // f(-1) does not divide F(-1), so go out.
            return false;
          }
        }
        // Divide both polynomials mod 32768. If the remainder is
        // not zero, the integer division will not be performed.
        TestNbr0Bak = TestNbr[0].x;
        TestNbr[0].x = 32768;
        modulusIsZero = false;
        DividePolynomial(dividendMod32768, polyS[0], divisorMod32768, poly5[0], NULL);
        modulusIsZero = true;
        TestNbr[0].x = TestNbr0Bak;
        // Test whether the remainder is zero.
        ptrMod32768 = dividendMod32768;
        for (currentDegree = 0; currentDegree < poly5[0]; currentDegree++)
        {
          if (*(ptrMod32768+1) != 0)
          {
            return false;         // Coefficient of remainder is not zero. Go out.
          }
          ptrMod32768 += 2;   // Point to next coefficient.
        }
        rc = DivideIntegerPolynomial(polyS, poly5, TYPE_MODULUS);
        if (rc == EXPR_POLYNOMIAL_DIVISION_NOT_INTEGER)
        {               // Cannot perform the division.
          return false;
        }
        if ((polyS[0] != 0) || (polyS[1] != 1) || (polyS[2] != 0))
        {              // Remainder is not zero.
          return false;    // Number to factor does not divide this polynomial.
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
          (void)BigIntDivide(&operand2, &operand4, &operand3);
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
      modulusIsZero = false;   // Perform modular operations.
      // Restart finding factors.
      (void)memset(arrNbrFactors, 0, sizeof(arrNbrFactors));
    }
  }
  return true;              // Factorization is complete.
}

// Find Knuth-Cohen bound for coefficients of polynomial factors:
// If polynomial B divides A we have for all j:
// |Bj| <= binomial(n-1, j)*SUM(i, |Ai|^2))^(1/2) + binomial(n-1, j-1) * Q
// where m is the degree of A and n is the degree of B.
// Q is the minimum of |A0| and |Am|.
// Maximum degree to be considered is n = ceil(m/2).
// We need to find max(Bj).
static void ComputeCoeffBounds(void)
{
  int degreePolyToFactor;
  const int* ptrSrc;
  int maxDegreeFactor;
  int degree1;

  modulusIsZero = true;
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
    (void)BigIntMultiply(&operand1, &operand1, &operand1);
    for (degree1 = 1; degree1 <= degreePolyToFactor; degree1++)
    {
      ptrSrc += 1 + numLimbs(ptrSrc);
      UncompressBigIntegerB(ptrSrc, &operand3);   // The last loop sets operand3 to the leading coefficient.
      (void)BigIntMultiply(&operand3, &operand3, &operand2);
      BigIntAdd(&operand1, &operand2, &operand1);
    }
    squareRoot(operand1.limbs, operand2.limbs, operand1.nbrLimbs, &operand2.nbrLimbs);
    CopyBigInt(&operand1, &operand2);  // operand1 <- SUM(i, |Ai|^2))^(1/2)
    operand3.sign = SIGN_POSITIVE;
  }
  // If the absolute value of trailing coefficient is less than the
  // absolute value of leading coefficient (stored in operand3), replace it.
  UncompressBigIntegerB(&polyNonRepeatedFactors[1], &operand2);
  operand2.sign = SIGN_POSITIVE;
  BigIntSubt(&operand3, &operand2, &operand4);
  if (operand4.sign == SIGN_POSITIVE)
  {      // |Am| >= |A0|
    CopyBigInt(&operand3, &operand2);
  }

  // Loop that finds the maximum value of bound for |Bj|.
  intToBigInteger(&operand2, 1);  // binomial(n-1, 0)
  UncompressBigIntegerB(&polyNonRepeatedFactors[1], &bound);  // bound <- |A0|
  bound.sign = SIGN_POSITIVE;
  for (degree1 = 1; degree1 <= maxDegreeFactor; degree1++)
  {
    CopyBigInt(&operand4, &operand2);
    multint(&operand2, &operand2, maxDegreeFactor - degree1);
    subtractdivide(&operand2, 0, degree1);  // operand2 <- binomial(n-1, j)
    (void)BigIntMultiply(&operand1, &operand2, &operand5);
    (void)BigIntMultiply(&operand3, &operand4, &operand4);
    BigIntAdd(&operand5, &operand4, &operand5);
    // If operand5 > bound, set bound to operand5.
    BigIntSubt(&operand5, &bound, &operand4);
    if (operand4.sign == SIGN_POSITIVE)
    {
      CopyBigInt(&bound, &operand5);
    }
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
  int nbrVectors;
  int r1;
  int nbrRequiredTraces;
  int firstTrace;
  int nbrRow;
  int nbrCol;
  int exponDifference;
  int delta;
  // Store into polyLifted the polynomial values / lc(values) (mod powerMod).
  double logPrime = log(prime);
  double b0 = log(2 * values[0]) / logPrime;
  double log_rootbound = logBigNbr(&bound) / logPrime;
  int a0;
  int b;
  int C;
  int j;
  int rank;
  const int* ptrSrc;
  struct sFactorInfo* pstFactorInfo;
  int degreePolyToFactor = polyNonRepeatedFactors[0];
  int newNumberLength;
  int newNbrFactors;
  int ctr1;
  int currentAttempts;
#ifdef __EMSCRIPTEN__
  int maxAttempts;
#endif
  numberLLL = 0;
  exponDifference = (int)(6.0 * (double)nbrFactors * log(2.0) / log(prime));
  b = (int)(b0 + (ceil(log(3) * log_rootbound) / logPrime) + 3);
  a0 = b + exponDifference;
  exponentMod = a0;
#if DEBUG_VANHOEIJ
  (void)strcpy(ptrDebugOutput, LF "====================================================="
    LF "prime = ");
  ptrDebugOutput += strlen(ptrDebugOutput);
  int2dec(&ptrDebugOutput, prime);
  (void)strcpy(ptrDebugOutput, "root bound = ");
  ptrDebugOutput += strlen(ptrDebugOutput);
  BigInteger2Dec(&bound, ptrDebugOutput, 0);
  ptrDebugOutput += strlen(ptrDebugOutput);
  *ptrDebugOutput++ = '\n';
#endif
  (void)memset(arrNbrFactors, 0, sizeof(arrNbrFactors));
  // Get leading coefficient of polyNonRepeatedFactors.
  ptrSrc = &polyNonRepeatedFactors[1];
  for (int degree1 = 0; degree1 < degreePolyToFactor; degree1++)
  {
    ptrSrc += 1 + numLimbs(ptrSrc);
  }
  computePower(exponentMod);
  modulusIsZero = false;    // Use modular arithmetic for polynomials.
  intToBigInteger(&operand5, 1);
  values[0] = polyNonRepeatedFactors[0];
  getModPolynomial(&values[1], polyNonRepeatedFactors, &operand5);
#if DEBUG_HENSEL_LIFTING
  ptrOutput2 = ptrDebugOutput;
#endif
  HenselLifting(factorInfoRecord, 1);
#if DEBUG_HENSEL_LIFTING
  ptrDebugOutput = ptrOutput2; ptrOutput2 = NULL;
#endif
  // Check whether the product of up to 3 factors mod p^b
  // is a divisor of the original polynomial.
  CopyBigInt(&halfPowerMod, &powerMod);
  subtractdivide(&halfPowerMod, -1, 2); // halfPowerMod <- (powerMod+1)/2
  GetMontgomeryParms(powerMod.nbrLimbs);
  // Polynomials in AttemptToFactor are multiplied by the leading
  // coefficient, so multiply the coefficient bound as well.
  (void)BigIntMultiply(&bound, &leadingCoeff, &bound);
  for (ctr1 = 0; ctr1 < nbrFactors; ctr1++)
  {
    intToLinkedBigInt(&lambda[0][ctr1], 0);
  }
  newNbrFactors = nbrFactors;
  for (ctr1 = 0; ctr1 < nbrFactors; ctr1++)
  {
    intToLinkedBigInt(&lambda[0][ctr1], 1);
    (void)AttemptToFactor(1, nbrFactors, &newNbrFactors);
    intToLinkedBigInt(&lambda[0][ctr1], 0);
  }
  currentAttempts = 0;
#ifdef __EMSCRIPTEN__
  maxAttempts = nbrFactors * (nbrFactors - 1) / 2;
#endif
  for (ctr1 = 0; ctr1 < nbrFactors; ctr1++)
  {
    intToLinkedBigInt(&lambda[0][ctr1], 1);
    for (int ctr2 = ctr1 + 1; ctr2 < nbrFactors; ctr2++)
    {
      intToLinkedBigInt(&lambda[0][ctr2], 1);
      (void)AttemptToFactor(1, nbrFactors, &newNbrFactors);
      currentAttempts++;
#ifdef __EMSCRIPTEN__
      int elapsedTime = (int)(tenths() - originalTenthSecond);
      if (elapsedTime / 10 != oldTimeElapsed / 10)
      {
        oldTimeElapsed = elapsedTime;
        ptrOutput = output;
        (void)strcpy(ptrOutput, lang? "1<p>Obteniendo factores de dos factores modulares: prueba ":
                                "1<p>Finding factors from two modular factors: attempt ");
        ptrOutput += strlen(ptrOutput);
        int2dec(&ptrOutput, currentAttempts);
        (void)strcpy(ptrOutput, lang ? " de " : " of ");
        ptrOutput += strlen(ptrOutput);
        int2dec(&ptrOutput, maxAttempts);
        (void)strcpy(ptrOutput, "</p>");
        ptrOutput += strlen(ptrOutput);
        showElapsedTimeSec(&ptrOutput);
        databack(output);
      }
#endif
      intToLinkedBigInt(&lambda[0][ctr2], 0);
    }
    intToLinkedBigInt(&lambda[0][ctr1], 0);
  }
  // Restore value of bound.
  (void)BigIntDivide(&bound, &leadingCoeff, &bound);
  if (newNbrFactors <= 1)
  {   // Zero or 1 factor left. Polynomial completely factored, so go out.
    return;
  }
  // Discard the entries with multiplicity zero in the stored factors array.
  pstFactorInfo = factorInfoRecord;
  delta = 0;
  for (int nbrFactor = 0; nbrFactor < nbrFactors; nbrFactor++)
  {
    if (pstFactorInfo->multiplicity == 0)
    {
      delta++;
    }
    else if (delta > 0)
    {
      (void)memcpy(pstFactorInfo - delta, pstFactorInfo, sizeof(*pstFactorInfo));
    }
    pstFactorInfo++;
  }
  nbrFactors = newNbrFactors;
  nbrFactorsFound = nbrFactors;
  // Compute powerExtraBits as a power of prime with 3 times
  // the number of bits of nbrFactors (the number of polynomial factors).
  // Compute powerBoundA as powerExtraBits * powerMod.
  exponDifference = (int)(6.0 * (double)nbrFactors * log(2.0) / log(prime));
  intToBigInteger(&operand1, prime);
  BigIntPowerIntExp(&operand1, exponDifference, &powerExtraBits);
  CopyBigInt(&powerBoundA, &powerMod);
  (void)BigIntMultiply(&powerExtraBits, &powerMod, &powerBoundA);
  ComputeCoeffBounds();     // bound = Bound of coefficient of factors.
  exponentMod = a0;
  computePower(exponentMod);
  modulusIsZero = false;    // Use modular arithmetic for polynomials.
  // Initialize matrix BL with identity matrix (length nbrFactors).
  for (nbrRow = 0; nbrRow < nbrFactors; nbrRow++)
  {
    for (nbrCol = 0; nbrCol < nbrFactors; nbrCol++)
    {
      intToLinkedBigInt(&matrixBL[nbrRow][nbrCol],
        (nbrRow == nbrCol ? 1 : 0));
    }
  }
  // Initialize variables for van Hoeij algorithm.
  nbrRequiredTraces = 1;        // Initialize nbrRequiredTraces.
  firstTrace = 1;               // First trace to use for matrix M.
  nbrVectors = nbrFactors;      // Initialize nbrVectors.
  for (;;)
  {
    int squareFormula;
    modulusIsZero = false;    // Use modular arithmetic for polynomials.
    firstTrace += nbrRequiredTraces - 1;
    if (nbrVectors >= nbrFactors)
    {
      nbrVectors = nbrFactors;
      // Initialize matrix BL with identity matrix (length nbrFactors).
      for (nbrRow = 0; nbrRow < nbrFactors; nbrRow++)
      {
        for (nbrCol = 0; nbrCol < nbrFactors; nbrCol++)
        {
          intToLinkedBigInt(&matrixBL[nbrRow][nbrCol],
            (nbrRow == nbrCol ? 1 : 0));
        }
      }
    }
// Get b such that p^b is greater than the bounds on the cofficients.
    computePower(a0);
    newNumberLength = NumberLength;
    exponentMod = a0;
    intToBigInteger(&operand1, prime);
    BigIntPowerIntExp(&operand1, b, &powerBoundA);
#if DEBUG_VANHOEIJ
    (void)strcpy(ptrDebugOutput, LF "====================================================="
      LF "prime = ");
    ptrDebugOutput += strlen(ptrDebugOutput);
    int2dec(&ptrDebugOutput, prime);
    (void)strcpy(ptrDebugOutput, ", a0 = ");
    ptrDebugOutput += strlen(ptrDebugOutput);
    int2dec(&ptrDebugOutput, a0);
    (void)strcpy(ptrDebugOutput, ", b0 = ");
    ptrDebugOutput += strlen(ptrDebugOutput);
    int2dec(&ptrDebugOutput, b);
    (void)strcpy(ptrDebugOutput, ", exponDifference = ");
    ptrDebugOutput += strlen(ptrDebugOutput);
    int2dec(&ptrDebugOutput, exponDifference);
    (void)strcpy(ptrDebugOutput, LF);
    ptrDebugOutput += strlen(ptrDebugOutput);
#endif
    intToBigInteger(&operand5, 1);
    values[0] = polyNonRepeatedFactors[0];
    getModPolynomial(&values[1], polyNonRepeatedFactors, &operand5);
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
        intToLinkedBigInt(&basisStar[nbrRow][nbrCol],
          (nbrRow == nbrCol ? C : 0));
      }
    }

    // Use matrix lambda to hold the traces.
    for (int factorNbr = 0; factorNbr < nbrFactors; factorNbr++)
    {
      ComputeTraces(firstTrace + nbrRequiredTraces, factorNbr);
#if DEBUG_VANHOEIJ
      (void)strcpy(ptrDebugOutput, "Traces: ");
      ptrDebugOutput += strlen(ptrDebugOutput);
      for (nbrRow = firstTrace; nbrRow < firstTrace + nbrRequiredTraces; nbrRow++)
      {
        getBigIntegerFromLinked(traces[nbrRow], &tmp5);
        BigInteger2Dec(&tmp5, ptrDebugOutput, 0);
        ptrDebugOutput += strlen(ptrDebugOutput);
        if (nbrRow < firstTrace + nbrRequiredTraces - 1)
        {
          *ptrDebugOutput++ = ',';
          *ptrDebugOutput++ = ' ';
        }
      }
      (void)strcpy(ptrDebugOutput, LF);
      ptrDebugOutput += strlen(ptrDebugOutput);
#endif
      // Use traces 1 to RequiredTraces
      // (i.e. Tr[1] to Tr[nbrRequiredTraces])
      // Convert traces matrix to range -powerExtraBits/2 to powerExtraBits/2.
      for (j = 0; j < nbrRequiredTraces; j++)
      {
        getBigIntegerFromLinked(traces[firstTrace + j], &operand1);
        BigIntMultiplyBy2(&operand1);
        BigIntSubt(&operand1, &powerExtraBits, &operand1); // operand1 <- 2*traces[j] - powerExtraBits
        getBigIntegerFromLinked(traces[firstTrace + j], &tmp5);
        if (operand1.sign == SIGN_POSITIVE)
        {
          BigIntSubt(&tmp5, &powerExtraBits, &tmp5);
        }
        setLinkedBigInteger(&lambda[j][factorNbr], &tmp5);
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
          getBigIntegerFromLinked(lambda[nbrRow][j], &tmp4);
          getBigIntegerFromLinked(matrixBL[j][nbrCol], &tmp5);
          (void)BigIntMultiply(&tmp4, &tmp5, &operand2);
          BigIntAdd(&operand1, &operand2, &operand1);
        }
        setLinkedBigInteger(&basisStar[nbrRow + nbrVectors][nbrCol], &operand1);
      }
    }

    // Initialize lower-right half of matrix M.

    for (nbrRow = nbrVectors; nbrRow < nbrVectors + nbrRequiredTraces; nbrRow++)
    {
      for (nbrCol = nbrVectors; nbrCol < nbrVectors + nbrRequiredTraces; nbrCol++)
      {
        if (nbrRow == nbrCol)
        {
          setLinkedBigInteger(&basisStar[nbrRow][nbrCol], &powerExtraBits);
        }
        else
        {
          intToLinkedBigInt(&basisStar[nbrRow][nbrCol], 0);
        }
      }
    }

    // Copy matrix M to basis.
    for (nbrRow = 0; nbrRow < nbrVectors + nbrRequiredTraces; nbrRow++)
    {
      for (nbrCol = 0; nbrCol < nbrVectors + nbrRequiredTraces; nbrCol++)
      {
        getBigIntegerFromLinked(basisStar[nbrRow][nbrCol], &tmp5);
        setLinkedBigInteger(&basis[nbrRow][nbrCol], &tmp5);
      }
    }
#if DEBUG_VANHOEIJ
    (void)strcpy(ptrDebugOutput, LF "Matrix M before LLL (");
    ptrDebugOutput += strlen(ptrDebugOutput);
    int2dec(&ptrDebugOutput, nbrVectors + nbrRequiredTraces);
    *ptrDebugOutput++ = '*';
    int2dec(&ptrDebugOutput, nbrVectors + nbrRequiredTraces);
    (void)strcpy(ptrDebugOutput, "): ");
    ptrDebugOutput += strlen(ptrDebugOutput);
    for (nbrRow = 0; nbrRow < nbrVectors + nbrRequiredTraces; nbrRow++)
    {
      for (nbrCol = 0; nbrCol < nbrVectors + nbrRequiredTraces; nbrCol++)
      {
        getBigIntegerFromLinked(basisStar[nbrRow][nbrCol], &tmp5);
        BigInteger2Dec(&tmp5, ptrDebugOutput, 0);
        ptrDebugOutput += strlen(ptrDebugOutput);
        if (nbrCol < nbrVectors + nbrRequiredTraces - 1)
        {
          *ptrDebugOutput++ = ',';
        }
        intToLinkedBigInt(&lambda[nbrRow][nbrCol], 0);  // DEBUG BORRAR
      }
      *ptrDebugOutput++ = ';';
    }
    *ptrDebugOutput++ = '.';

#endif
    // Step 2: LLL-reduce the (r+s)*(r+s) matrix M (of rank r+s).
    integralLLL(nbrVectors + nbrRequiredTraces);

    // Copy matrix M to basisStar.
    for (nbrRow = 0; nbrRow < nbrVectors + nbrRequiredTraces; nbrRow++)
    {
      for (nbrCol = 0; nbrCol < nbrVectors + nbrRequiredTraces; nbrCol++)
      {
        getBigIntegerFromLinked(basis[nbrRow][nbrCol], &tmp5);
        setLinkedBigInteger(&basisStar[nbrRow][nbrCol], &tmp5);
      }
    }
#if DEBUG_VANHOEIJ
    (void)strcpy(ptrDebugOutput, LF "Matrix M after LLL: " LF);
    ptrDebugOutput += strlen(ptrDebugOutput);
    for (nbrRow = 0; nbrRow < nbrVectors + nbrRequiredTraces; nbrRow++)
    {
      for (nbrCol = 0; nbrCol < nbrVectors + nbrRequiredTraces; nbrCol++)
      {
        getBigIntegerFromLinked(basisStar[nbrRow][nbrCol], &tmp5);
        BigInteger2Dec(&tmp5, ptrDebugOutput, 0);
        ptrDebugOutput += strlen(ptrDebugOutput);
        *ptrDebugOutput++ = ' ';
      }
      *ptrDebugOutput++ = ';';
    }
    *ptrDebugOutput++ = '.';
    *ptrDebugOutput = 0;
#endif
    // Step 3: Replace the upper r*(r+s) submatrix L of M by BL*L.
    // Then set BL <- BL*L of dimension n0*(r+s)
    // M is now of dimension (n0+s)*(r+s) of rank r+s.
    // First move down the rows of M.
    for (nbrRow = nbrVectors + nbrRequiredTraces - 1; nbrRow >= nbrVectors; nbrRow--)
    {
      struct linkedBigInt** ptrLBISrc = &basisStar[nbrRow][0];
      struct linkedBigInt** ptrLBIDest = &basisStar[nbrRow - nbrVectors + nbrFactors][0];
      for (nbrCol = 0; nbrCol < nbrVectors + nbrRequiredTraces; nbrCol++)
      {
        getBigIntegerFromLinked(*ptrLBISrc, &tmp5);
        ptrLBISrc++;
        setLinkedBigInteger(ptrLBIDest, &tmp5);
        ptrLBIDest++;
      }
    }
    // Copy L to array lambda.
    for (; nbrRow >= 0; nbrRow--)
    {
      struct linkedBigInt** ptrLBISrc = &basisStar[nbrRow][0];
      struct linkedBigInt** ptrLBIDest = &lambda[nbrRow][0];
      for (nbrCol = 0; nbrCol < nbrVectors + nbrRequiredTraces; nbrCol++)
      {
        getBigIntegerFromLinked(*ptrLBISrc, &tmp5);
        ptrLBISrc++;
        setLinkedBigInteger(ptrLBIDest, &tmp5);
        ptrLBIDest++;
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
          getBigIntegerFromLinked(matrixBL[nbrRow][j], &tmp4);
          getBigIntegerFromLinked(lambda[j][nbrCol], &tmp5);
          (void)BigIntMultiply(&tmp4, &tmp5, &operand2);
          BigIntAdd(&operand1, &operand2, &operand1);
        }
        setLinkedBigInteger(&basisStar[nbrRow][nbrCol], &operand1);
      }
    }

    // Copy BL*L to array matrixBL.
    for (nbrRow = 0; nbrRow < nbrFactors; nbrRow++)
    {
      struct linkedBigInt** ptrLBISrc = &basisStar[nbrRow][0];
      struct linkedBigInt** ptrLBIDest = &matrixBL[nbrRow][0];
      for (nbrCol = 0; nbrCol < nbrVectors; nbrCol++)
      {
        getBigIntegerFromLinked(*ptrLBISrc, &tmp5);
        ptrLBISrc++;
        setLinkedBigInteger(ptrLBIDest, &tmp5);
        ptrLBIDest++;
      }
    }

#if DEBUG_VANHOEIJ
    (void)strcpy(ptrDebugOutput, LF "Matrix M before Gram-Schmidt" LF);
    ptrDebugOutput += strlen(ptrDebugOutput);
    for (nbrRow = 0; nbrRow < nbrVectors + nbrRequiredTraces; nbrRow++)
    {
      for (nbrCol = 0; nbrCol < nbrVectors + nbrRequiredTraces; nbrCol++)
      {
        getBigIntegerFromLinked(basisStar[nbrRow][nbrCol], &tmp5);
        BigInteger2Dec(&tmp5, ptrDebugOutput, 0);
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
    (void)strcpy(ptrDebugOutput, LF "Norms after Gram-Schmidt: ");
    ptrDebugOutput += strlen(ptrDebugOutput);
    for (r1 = 2; r1 <= nbrVectors + nbrRequiredTraces; r1++)
    {
      getBigIntegerFromLinked(lambda[r1 - 1][r1 - 1], &tmp4);
      getBigIntegerFromLinked(lambda[r1 - 2][r1 - 2], &tmp5);
      (void)BigIntDivide(&tmp4, &tmp5, &operand1);
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
      getBigIntegerFromLinked(lambda[r1 - 1][r1 - 1], &tmp4);
      getBigIntegerFromLinked(lambda[r1 - 2][r1 - 2], &tmp5);
      (void)BigIntDivide(&tmp4, &tmp5, &operand1);
      if ((operand1.nbrLimbs == 1) && (operand1.limbs[0].x <= squareFormula))
      {
        break;        // r' was found.
      }
    }
#if DEBUG_VANHOEIJ
    (void)strcpy(ptrDebugOutput, LF "squareFormula = ");
    ptrDebugOutput += strlen(ptrDebugOutput);
    int2dec(&ptrDebugOutput, squareFormula);
#endif

    // Step 6: If r' = 1, return "irreducible".
    if (r1 == 1)
    {                 // Polynomial is irreducible.
#if DEBUG_VANHOEIJ
      (void)strcpy(ptrDebugOutput, LF "r' = 1 -> irreducible." LF);
#ifdef __EMSCRIPTEN__
      output[0] = '1';
      (void)strcpy(&output[1], debugOutput);
      databack(output);
#else
      printf("%s", debugOutput);
#endif
#endif
      return;
    }

#if DEBUG_VANHOEIJ
#ifdef __EMSCRIPTEN__
    {
      char outputInfo[1000];
      char* ptrOutput = outputInfo;
      if (lang)
      {
        (void)strcpy(ptrOutput, "1<p>Paso 7</p>");
      }
      else
      {
        (void)strcpy(ptrOutput, "1<p>Step 7</p>");
      }
      databack(outputInfo);
    }
#endif
#endif
    // Step 7: BL <- BL / C. Now BL has dimension 
    // nbrFactors rows by r'
    for (nbrRow = 0; nbrRow < nbrFactors; nbrRow++)
    {
      struct linkedBigInt** ptrLBIDest = &matrixBL[nbrRow][0];
      for (nbrCol = 0; nbrCol < r1; nbrCol++)
      {
        getBigIntegerFromLinked(*ptrLBIDest, &tmp5);
        subtractdivide(&tmp5, 0, C);
        setLinkedBigInteger(ptrLBIDest, &tmp5);
        ptrLBIDest++;
      }
    }

#if DEBUG_VANHOEIJ
    (void)strcpy(ptrDebugOutput, LF "Matrix BL before Gauss: ");
    ptrDebugOutput += strlen(ptrDebugOutput);
    for (nbrRow = 0; nbrRow < nbrFactors; nbrRow++)
    {
      for (nbrCol = 0; nbrCol < r1; nbrCol++)
      {
        getBigIntegerFromLinked(basisStar[nbrRow][nbrCol], &tmp5);
        BigInteger2Dec(&tmp5, ptrDebugOutput, 0);
        ptrDebugOutput += strlen(ptrDebugOutput);
        *ptrDebugOutput++ = ' ';
      }
      *ptrDebugOutput++ = ';';
    }
    *ptrDebugOutput++ = '.';
    *ptrDebugOutput = 0;
#ifdef __EMSCRIPTEN__
    {
      char outputInfo[1000];
      char* ptrOutput = outputInfo;
      if (lang)
      {
        (void)strcpy(ptrOutput, "1<p>Paso 8</p>");
      }
      else
      {
        (void)strcpy(ptrOutput, "1<p>Step 8</p>");
      }
      databack(outputInfo);
    }
#endif
#endif
    // Step 8: Put BL in Gauss-Jordan form. If possibly
    // after division by a constant factor, each column
    // contains only 0 and 1, and each row contains 
    // exactly one 1, we might have a valid factorization,
    // done <- true. Otherwise, increase s, replace n
    // by r and r by r'
    rank = gauss(nbrFactors, r1);  // Dimension of matrix BL.
#if DEBUG_VANHOEIJ
    (void)strcpy(ptrDebugOutput, LF "squareFormula = ");
    ptrDebugOutput += strlen(ptrDebugOutput);
    int2dec(&ptrDebugOutput, squareFormula);
    (void)strcpy(ptrDebugOutput, ", rank = ");
    ptrDebugOutput += strlen(ptrDebugOutput);
    int2dec(&ptrDebugOutput, rank);
    (void)strcpy(ptrDebugOutput, LF);
    ptrDebugOutput += strlen(ptrDebugOutput);
#endif
    nbrRequiredTraces++;
    nbrVectors = r1;
#if 0
    if (++nbrStepsDone == 10)
    {
#ifdef __EMSCRIPTEN__
      output[0] = '1';
      (void)strcpy(&output[1], debugOutput);
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
    (void)strcpy(ptrDebugOutput, LF "Matrix R after Gauss: ");
    ptrDebugOutput += strlen(ptrDebugOutput);
    for (nbrRow = 0; nbrRow < nbrFactors; nbrRow++)
    {
      for (nbrCol = 0; nbrCol < r1; nbrCol++)
      {
        getBigIntegerFromLinked(lambda[nbrRow][nbrCol], &tmp5);
        BigInteger2Dec(&tmp5, ptrDebugOutput, 0);
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
        if (!linkedBigIntIsZero(lambda[nbrRow][nbrCol]))
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
      (void)strcpy(ptrDebugOutput, "nbrCol = ");
      ptrDebugOutput += strlen(ptrDebugOutput);
      int2dec(&ptrDebugOutput, nbrCol);
      (void)strcpy(ptrDebugOutput, ", nbrFactors = ");
      ptrDebugOutput += strlen(ptrDebugOutput);
      int2dec(&ptrDebugOutput, nbrFactors);
      (void)strcpy(ptrDebugOutput, "They should have been equal." LF);
      ptrDebugOutput += strlen(ptrDebugOutput);
#endif
      continue;
    }
    CopyBigInt(&halfPowerMod, &powerMod);
    subtractdivide(&halfPowerMod, -1, 2); // halfPowerMod <- (powerMod+1)/2
    newNbrFactors = nbrFactors;
    // Combine factors.
    GetMontgomeryParms(powerMod.nbrLimbs);
#ifdef __EMSCRIPTEN__
    char outputInfo[1000];
    char* ptrOutput = outputInfo;
    if (lang)
    {
      (void)strcpy(ptrOutput, "1<p>Verificando si los polinomios hallados son factores irreducibles.</p>");
    }
    else
    {
      (void)strcpy(ptrOutput, "1<p>Testing whether the polynomials found are irreducible factors.</p>");
    }
    ptrOutput += strlen(ptrOutput);
    showElapsedTimeSec(&ptrOutput);
    databack(outputInfo);
#endif
    if (AttemptToFactor(nbrVectors, nbrFactors, &newNbrFactors))
    {
#if DEBUG_VANHOEIJ
      (void)strcpy(ptrDebugOutput, "nbrVector == nbrVectors" LF);
#ifdef __EMSCRIPTEN__
      output[0] = '1';
      (void)strcpy(&output[1], debugOutput);
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
static void GenerateIntegerPolynomial(const int* polyMod, int* polyInt, int degreePoly)
{
  *polyInt = degreePoly;
  const int* ptrSrc = polyMod;
  int* ptrDest = polyInt + 1;
  for (int currentDegree = 0; currentDegree <= degreePoly; currentDegree++)
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
  struct sFactorInfo* pstFactorInfoInteger;
  const struct sFactorInfo* pstFactorInfo;
  int indexNewFactor[MAX_DEGREE];
  int indexOldFactor[MAX_DEGREE];
  int currentDegree;
  int index;
  int* ptrIndex;
  const int* ptrOldFactor;

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
  pstFactorInfo = pstFactorInfoInteger;
  while (pstFactorInfo->ptrPolyLifted != NULL)
  {
    pstFactorInfo++;
  }
  if (pstFactorInfo - pstFactorInfoInteger > 0)
  {
    (void)memmove(pstFactorInfoInteger + 1, pstFactorInfoInteger,
      (pstFactorInfo - pstFactorInfoInteger) * sizeof(*pstFactorInfo));
  }
  pstFactorInfoInteger->multiplicity = 1;
  pstFactorInfoInteger->ptrPolyLifted = ptrFactor;
  pstFactorInfoInteger->degree = degreePoly;
  ptrFactorInteger += indexNewFactor[degreePoly + 1];
}

int getNextPrimeNoDuplicatedFactors(int primeIndex)
{
  const int* ptrSrc;
  int prime;
  int degreeGcdMod;
  int degree2;
  int polyDegree = polyNonRepeatedFactors[0];
  initializeSmallPrimes(smallPrimes);
  // Get leading coefficient of polyNonRepeatedFactors.
  ptrSrc = &polyNonRepeatedFactors[1];
  for (int currentDegree = 0; currentDegree < polyDegree; currentDegree++)
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
    int degree1;
    do
    {      // Loop while the leading coefficient is multiple of prime.
      prime = smallPrimes[++primeIndex];
    } while (getRemainder(&leadingCoeff, prime) == 0);
    modulusIsZero = false;
    intToBigInteger(&primeMod, prime);
    computePower(1);
    intToBigInteger(&operand5, 1);
    degree2 = getModPolynomial(&poly2[1], polyNonRepeatedFactors, &operand5);
    poly2[0] = degree2;
    DerPolynomial(poly2);   // This function overwrites poly1.
    degree1 = getModPolynomial(poly1, polyNonRepeatedFactors, &operand5);
    PolyModularGcd(poly1, degree1, &poly2[1], poly2[0], poly3, &degreeGcdMod);
  } while (degreeGcdMod > 0);
  return primeIndex;
}

static void initFactorModularPoly(int prime)
{
  int degreePolyToFactor;
  modulusIsZero = true;
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
  modulusIsZero = false;
}

void PerformSameDegreeFactorization(int prime)
{
  initFactorModularPoly(prime);
  SameDegreeFactorization();
}

void FactorPolynomialModPrime(int prime)
{
  (void)memset(factorInfo, 0, sizeof(factorInfo));
  initFactorModularPoly(prime);
  FactorModularPolynomial(false);   // Input is not in Montgomery notation.
}

static void CopyFactorsFoundToRecord(void)
{
  const struct sFactorInfo* pstFactorInfoOrig = factorInfo;
  struct sFactorInfo* pstFactorInfoRecord = factorInfoRecord;
  int *ptrPolyLiftedRecord = polyLiftedRecord;
  for (int factorNbr = 0; factorNbr < MAX_DEGREE; factorNbr++)
  {
    if (pstFactorInfoOrig->ptr == NULL)
    {    // No more factors.
      break;
    }
    *pstFactorInfoRecord = *pstFactorInfoOrig;
    pstFactorInfoRecord->ptr = ptrPolyLiftedRecord;
    ptrPolyLiftedRecord = CopyPolynomialFixedCoeffSize(ptrPolyLiftedRecord,
      pstFactorInfoOrig->ptr,
      pstFactorInfoOrig->degree - 1, primeMod.nbrLimbs + 1);
    *ptrPolyLiftedRecord++ = 1;    // Leading coefficient should be 1.
    *ptrPolyLiftedRecord++ = 1;
    pstFactorInfoOrig++;
    pstFactorInfoRecord++;
  }
  pstFactorInfoRecord->ptr = NULL; // Indicate that there are no more factors.
}

// Input: values = degree, coefficient degree 0, coefficient degree 1, etc.
// Output: factorInfo = structure that holds the factors.
int FactorPolyOverIntegers(void)
{
  int degreePolyToFactor = values[0];
  int primeRecord = 0;
  int expon;
  const int* ptrSrc;
  int* ptrDest;
  int factorNbr;
  bool polXprocessed = false;
  int* ptrFactorIntegerBak;
  int* ptrPolyLiftedOrig;
  struct sFactorInfo* pstFactorInfoOrig;
  struct sFactorInfo* pstFactorInfoInteger = factorInfoInteger;
  initLinkedBigInt();
  ptrFactorInteger = polyInteger;
  modulusIsZero = true;
  (void)memset(factorInfoInteger, 0, sizeof(factorInfoInteger));
  getContent(values, &contentPolyToFactor);
  CopyPolynomial(&origPolyToFactor[1], &values[1], degreePolyToFactor);
  origPolyToFactor[0] = degreePolyToFactor;
  // polyToFactor -> original polynomial / content of polynomial.
  // Let n the degree of the least coefficient different from zero.
  // Then x^n divides the polynomial.
  polyToFactor[0] = degreePolyToFactor;
  ptrSrc = &origPolyToFactor[1];
  ptrDest = &polyToFactor[1];
  for (int currentDegree = 0; currentDegree <= degreePolyToFactor; currentDegree++)
  {
    UncompressBigIntegerB(ptrSrc, &operand1);
    if (polXprocessed == false)
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
      polXprocessed = true;
    }
    (void)BigIntDivide(&operand1, &contentPolyToFactor, &operand2);
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
    int primeIndex;
    const struct sFactorInfo* pstFactorInfoRecord;
    modulusIsZero = true;
    // Get trailing coefficient.
    ptrSrc = &values[1];
    UncompressBigIntegerB(ptrSrc, &operand1);
    // Get leading coefficient.
    for (int degree1 = 0; degree1 < degreePolyToFactor; degree1++)
    {
      ptrSrc += 1 + numLimbs(ptrSrc);
    }
    UncompressBigIntegerB(ptrSrc, &operand2);
    (void)BigIntMultiply(&operand1, &operand2, &trailingCoeff);
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
    primeIndex = 0;
    ComputeCoeffBounds();   // bound = Bound of coefficient of factors.
    // Up to 5 different prime moduli are tested to minimize the number
    // of factors because the Van Hoeij algorithm speed depends on
    // this number. If the number of factors is less than 10,
    // no more modular factorizations are attempted.

    // The big number ensures that the first copy is done.
    nbrFactorsRecord = 100000;
    for (int attemptNbr = 1; attemptNbr < 5; attemptNbr++)
    {
      int nbrFactors;
      primeIndex = getNextPrimeNoDuplicatedFactors(primeIndex);
      prime = smallPrimes[primeIndex];
      // Find expon such that prime^expon >= 2 * bound
      BigIntMultiplyBy2(&bound);
      intToBigInteger(&operand1, prime);
      expon = 1;
      for (;;)
      {
        // Compare operand1 = prime^expon against 2 * bound.
        BigIntSubt(&operand1, &bound, &operand4);
        if (operand4.sign == SIGN_POSITIVE)
        {
          break;     // prime^expon >= 2 * bound -> go out.
        }
        multint(&operand1, &operand1, prime);
        expon++;
      }
      FactorPolynomialModPrime(prime);
          // Get number of factors found.
      pstFactorInfoOrig = factorInfo;
      nbrFactors = 0;
      for (factorNbr = 0; factorNbr < MAX_DEGREE; factorNbr++)
      {
        if (pstFactorInfoOrig->ptr == NULL)
        {    // No more factors.
          break;
        }
        nbrFactors += pstFactorInfoOrig->degree / pstFactorInfoOrig->expectedDegree;        
        pstFactorInfoOrig++;
      }
      if (nbrFactors < nbrFactorsRecord)
      {    // Copy factors found to records arrays.
        primeRecord = prime;
        nbrFactorsRecord = nbrFactors;
        CopyFactorsFoundToRecord();
        if (nbrFactors < 10)
        {
          break;    // Small enough number of factors. Go out of loop.
        }
      }
    }
    // Modulus that generate the lowest number of polynomials factors found.
    // Perform same degree factorization.
    // Copy back the record factorization to the work area.
    pstFactorInfoOrig = factorInfo;
    pstFactorInfoRecord = factorInfoRecord;
    ptrPolyLiftedOrig = polyLifted;
    for (factorNbr = 0; factorNbr < MAX_DEGREE; factorNbr++)
    {
      if (pstFactorInfoRecord->ptr == NULL)
      {    // No more factors.
        break;
      }
      *pstFactorInfoOrig = *pstFactorInfoRecord;
      pstFactorInfoOrig->ptr = ptrPolyLiftedOrig;
      ptrPolyLiftedOrig = CopyPolynomialFixedCoeffSize(ptrPolyLiftedOrig,
        pstFactorInfoRecord->ptr,
        pstFactorInfoRecord->degree - 1, primeMod.nbrLimbs + 1);
      pstFactorInfoOrig++;
      pstFactorInfoRecord++;
    }
    prime = primeRecord;
    nbrFactorsFound = factorNbr;
    PerformSameDegreeFactorization(prime);
    nbrFactorsRecord = nbrFactorsFound;
    CopyFactorsFoundToRecord();
    if (nbrFactorsRecord > 1)
    {
      vanHoeij(prime, nbrFactorsRecord);
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
  modulusIsZero = true;
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

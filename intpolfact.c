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
#include "expression.h"
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
extern char* ptrOutput;
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
static int intPolyMultiplicity;
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
extern int poly4[1000000];
extern int poly5[1000000];
int polyA[1000000];
int polyB[1000000];
int polyC[1000000];
int polyD[1000000];
int polySqFreeFact[1000000];
int polyInteger[1000000];
int primeEisenstein;
static void GenerateIntegerPolynomial(const int* polyMod, int* polyInt, int degreePoly);
static void InsertIntegerPolynomialFactor(int* ptrFactor, int degreePoly);
static bool EisensteinCriterion(const int* poly);
static bool checkEisenstein(const BigInteger* gcdAll,
  const int* ptrLeading, const int* ptrTrailing);
static int numberLLL;
static unsigned int validDegrees[((MAX_DEGREE / (int)sizeof(int)) / 8) + 2];
static unsigned int validDegreesRecord[((MAX_DEGREE / (int)sizeof(int)) / 8) + 2];

// Generate row echelon form from matrixBL.
// The output will be located in matrix lambda.
static int gauss(int nbrCols, int nbrRows)
{
  int k = 0;
  int l = 0;
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
  while ((k < nbrCols) && (l < nbrRows))
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
      k++;
      continue;
    }

    // We have found a non-zero element on the k-th column
    while ((row < nbrRows) && (pos == -1))
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
      row++;
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
        // Compute new value of M[i] as M[i] - M[l]
        getBigIntegerFromLinked(lambda[row][k], &tmp1);
        if ((row != l) && !BigIntIsZero(&tmp1))
        {
          intToLinkedBigInt(&lambda[row][k], 0);
          for (col = k + 1; col < nbrCols; col++)
          {
            // Compute new value of *x as *x + (*y)*t1
            getBigIntegerFromLinked(lambda[l][col], &tmp5);
            (void)BigIntMultiply(&tmp5, &tmp1, &tmp0);
            getBigIntegerFromLinked(lambda[row][col], &tmp5);
            BigIntSubt(&tmp5, &tmp0, &tmp5);
            setLinkedBigInteger(&lambda[row][col], &tmp5);
          }
        }
      }
      l++;
      k++;
    }
    else
    {
      return -1;    // No 1 found as pivot.
    }
  }
  return l;
}


// Compute the matrix lambda of Gram-Schmidt orthogonalization.
// Use algorithm 2 of "Symplectic Lattice Reduction and NTRU"
// (Gama1, Howgrave - Graham, Nguyen)
static void GramSchmidtOrthogonalization(int nbrRows, int nbrCols)
{
  int colJ;
  int k;
  for (int colI = 0; colI < nbrCols; colI++)
  {
    if (colI == 0)
    {             // Set element U_0 to 1
      intToLinkedBigInt(&detProdB[0], 1);
    }
    else
    {             // Set element U_i to lambda_{i-1, i-1}
      getBigIntegerFromLinked(lambda[colI - 1][colI - 1], &tmp5);
      setLinkedBigInteger(&detProdB[colI], &tmp5);

                  // Set element U_{i-1} to -lambda_{i, i-1}
      getBigIntegerFromLinked(lambda[colI][colI - 1], &tmp5);
      BigIntChSign(&tmp5);
      setLinkedBigInteger(&detProdB[colI - 1], &tmp5);
    }
    for (colJ = colI - 2; colJ >= 0; colJ--)
    {             // Set U_j to 0
      intToBigInteger(&tmp1, 0);
      for (k = colJ + 1; k <= colI; k++)
      {           // Set U_j to U_j + lambda_{k,j} * U_k
        getBigIntegerFromLinked(lambda[k][colJ], &tmp4);
        getBigIntegerFromLinked(detProdB[k], &tmp5);
        (void)BigIntMultiply(&tmp4, &tmp5, &tmp2);
        BigIntSubt(&tmp1, &tmp2, &tmp1);
      }
      // Set U_j to U_j / lambda_{j,j}
      getBigIntegerFromLinked(lambda[colJ][colJ], &tmp4);
      (void)BigIntDivide(&tmp1, &tmp4, &tmp5);
      setLinkedBigInteger(&detProdB[colJ], &tmp5);
    }
    for (colJ = colI; colJ < nbrCols; colJ++)
    {           // Set lambda_{j,i} to 0. Use tmp1 for lambda_{j,i}
      intToBigInteger(&tmp1, 0);
      for (k = 0; k <= colI; k++)
      {         // Set tmp2 to scalar product b_j * b_k 
        intToBigInteger(&tmp2, 0);
        for (int row = 0; row < nbrRows; row++)
        {
          getBigIntegerFromLinked(basisStar[row][colJ], &tmp4);
          getBigIntegerFromLinked(basisStar[row][k], &tmp5);
          (void)BigIntMultiply(&tmp4, &tmp5, &tmp3);
          BigIntAdd(&tmp2, &tmp3, &tmp2);
        }
        // Set matrix cell lambda_{j,i} to lambda_{j,i} + tmp2 * U_k
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
  // Set tmp0 to 2 lambda_{k, l}
  multint(&tmp0, &tmp4, 2);
  // Set tmp0 to |2 lambda_{k, l}|
  tmp0.sign = SIGN_POSITIVE;
  getBigIntegerFromLinked(detProdB[l], &tmp5);
  // Compute tmp1 as d_l - |2 lambda_{k, l}|
  BigIntSubt(&tmp5, &tmp0, &tmp1);
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
    // Subtract q*d_l from lambda_{k, l}.
  getBigIntegerFromLinked(detProdB[l], &tmp5);
  (void)BigIntMultiply(&tmp2, &tmp5, &tmp0);    // tmp0 <- q*d_l
  getBigIntegerFromLinked(lambda[k][l], &tmp5);
  BigIntSubt(&tmp5, &tmp0, &tmp5);
  setLinkedBigInteger(&lambda[k][l], &tmp5);
  for (int i = 0; i < l; i++)
  { // Subtract q*lambda_{l, i} from lambda_{k, i}.
    getBigIntegerFromLinked(lambda[l][i], &tmp5);
    // Set tmp0 to q*lambda_{l, i}
    (void)BigIntMultiply(&tmp2, &tmp5, &tmp0);
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
  for (int j = 0; j < (k - 1); j++)
  {  // Exchange lambda_{k, j} with lambda_{k-1, j}
    pstLinkedBigInt = lambda[k][j];
    lambda[k][j] = lambda[k - 1][j];
    lambda[k - 1][j] = pstLinkedBigInt;
  }
    // Set lambda to lambda_{k, k-1}
  getBigIntegerFromLinked(lambda[k][k - 1], &tmp0);    // tmp0 <- lambda.
    // Set B to (d_{k-2}*d_k + lambda^2)/d_{k-1}
    // d_{k-2}*d_k + lambda^2 is already in tmp3.
  getBigIntegerFromLinked(detProdB[k - 1], &tmp5);
  (void)BigIntDivide(&tmp3, &tmp5, &tmp1);                   // tmp1 <- B
#if DEBUG_VANHOEIJ
  if ((size == 16) && (z < 1000))
  {
    z++;
    copyStr(&ptrDebugOutput, "k = ");
    int2dec(&ptrDebugOutput, k);
    copyStr(&ptrDebugOutput, ", z = ");
    int2dec(&ptrDebugOutput, z);
    *ptrDebugOutput = '\r';
    ptrDebugOutput++;
    *ptrDebugOutput = '\n';
    ptrDebugOutput++;
    copyStr(&ptrDebugOutput, "tmp3 = ");
    BigInteger2Dec(&ptrDebugOutput, &tmp3, 0);
    *ptrDebugOutput = '\r';
    ptrDebugOutput++;
    *ptrDebugOutput = '\n';
    ptrDebugOutput++;
    copyStr(&ptrDebugOutput, "detProdB[k-1] = ");
    getBigIntegerFromLinked(detProdB[k - 1], &tmp5);
    BigInteger2Dec(&ptrDebugOutput, &tmp5, 0);
    *ptrDebugOutput = '\r';
    ptrDebugOutput++;
    *ptrDebugOutput = '\n';
    ptrDebugOutput++;
    copyStr(&ptrDebugOutput, "B = ");
    BigInteger2Dec(&ptrDebugOutput, &tmp1, 0);
    *ptrDebugOutput = '\r';
    ptrDebugOutput++;
    *ptrDebugOutput = '\n';
    ptrDebugOutput++;
  }
#endif
  for (int i = k+1; i <= kMax; i++)
  {
    // Set t to lambda_{i, k}
    // Set tmp2 to t
    getBigIntegerFromLinked(lambda[i][k], &tmp2);
    // Set matrix cell lambda_{i, k} to (d_k*lambda_{i, k-1} minus lambda * t)/d_{k-1}
    getBigIntegerFromLinked(detProdB[k], &tmp4);
    getBigIntegerFromLinked(lambda[i][k - 1], &tmp5);
    (void)BigIntMultiply(&tmp4, &tmp5, &tmp3);
    (void)BigIntMultiply(&tmp0, &tmp2, &tmp4);  // tmp4 <- lambda * t
    BigIntSubt(&tmp3, &tmp4, &tmp3);
    getBigIntegerFromLinked(detProdB[k - 1], &tmp5);
    (void)BigIntDivide(&tmp3, &tmp5, &tmp4);
    setLinkedBigInteger(&lambda[i][k], &tmp4);
    // Set matrix cell lambda_{i, k-1} to (B*t + lambda * lambda_{i, k})/d_k
    (void)BigIntMultiply(&tmp1, &tmp2, &tmp3);  // tmp3 <- B * t
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
  char* ptrOut = outputInfo;
  numberLLL++;
  if (lang)
  {
    copyStr(&ptrOut, "1<p>Calculando LLL número ");
    int2dec(&ptrOut, numberLLL);
    copyStr(&ptrOut, " en matriz de");
  }
  else
  {
    copyStr(&ptrOut, "1<p>Computing LLL #");
    int2dec(&ptrOut, numberLLL);
    copyStr(&ptrOut, " in matrix of");
  }
  *ptrOut = ' ';
  ptrOut++;
  int2dec(&ptrOut, size);
  copyStr(&ptrOut, " &times; ");
  int2dec(&ptrOut, size);
  copyStr(&ptrOut, ".</p>");
  showElapsedTimeSec(&ptrOut);
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
        { // On each loop, set u to (d_i * u - lambda_{k, i} * lambda_{j, i}) / d_{i-1}
          getBigIntegerFromLinked(detProdB[colI], &tmp5);
          (void)BigIntMultiply(&tmp5, &tmp2, &tmp0);  // d_i * u
          getBigIntegerFromLinked(lambda[colK][colI], &tmp4);
          getBigIntegerFromLinked(lambda[colJ][colI], &tmp5);
          // Compute lambda_{k, i} times lambda_{j, i}.
          (void)BigIntMultiply(&tmp4, &tmp5, &tmp1);
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
      copyStr(&ptrDebugOutput, "lambda: ");
      for (row = 0; row < size; row++)
      {
        copyStr(&ptrDebugOutput, "Row #");
        int2dec(&ptrDebugOutput, row);
        *ptrDebugOutput = ':';
        ptrDebugOutput++;
        *ptrDebugOutput = ' ';
        ptrDebugOutput++;
        for (colI = 0; colI < size; colI++)
        {
          getBigIntegerFromLinked(lambda[row][colI], &tmp5);
          BigInteger2Dec(&ptrDebugOutput, &tmp5, 0);
          *ptrDebugOutput = ',';
          ptrDebugOutput++;
          *ptrDebugOutput = ' ';
          ptrDebugOutput++;
        }
        *(ptrDebugOutput - 2) = '\r';
        *(ptrDebugOutput - 1) = '\n';
      }
      copyStr(&ptrDebugOutput, "detProdB: ");
      for (colI = 0; colI < size; colI++)
      {
        BigInteger2Dec(&ptrDebugOutput, &detProdB[colI], 0);
        *ptrDebugOutput = ',';
        ptrDebugOutput++;
        *ptrDebugOutput = ' ';
        ptrDebugOutput++;
      }
      *(ptrDebugOutput - 2) = '\r';
      *(ptrDebugOutput - 1) = '\n';
    }
#endif
    // Test LLL condition
    for (;;)
    {
      PerformREDI(colK, colK - 1, size);
      // Check whether 4 * d_k * d_{k-2}
      // is less than 3 (d_{k-1})^2 minus 4*(lambda_{k, k-1})^2.
      // That means 4 (d_k * d_{k-2} + lambda_{k, k-1})^2)
      // is less than 3 (d_{k-1})^2.
                       // Get d_k
      getBigIntegerFromLinked(detProdB[colK], &tmp0);
      if (colK > 1)
      {                // Multiply by d_{k-2} if k>=0
        getBigIntegerFromLinked(detProdB[colK - 2], &tmp5);
        (void)BigIntMultiply(&tmp0, &tmp5, &tmp0);
      }
                       // Compute the square of matrix cell lambda_{k, k-1}.
      getBigIntegerFromLinked(lambda[colK][colK - 1], &tmp5);
      (void)BigIntMultiply(&tmp5, &tmp5, &tmp1);
                       // Compute d_k * d_{k-2} plus the square of lambda_{k, k-1}.
      BigIntAdd(&tmp0, &tmp1, &tmp3);
      multint(&tmp0, &tmp3, 4);           // tmp0 is the Left Hand Side.
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
      copyStr(&ptrDebugOutput, "basis: ");
      for (row = 0; row < size; row++)
      {
        copyStr(&ptrDebugOutput, "Row #");
        int2dec(&ptrDebugOutput, row);
        *ptrDebugOutput = ':';
        ptrDebugOutput++;
        *ptrDebugOutput = ' ';
        ptrDebugOutput++;
        for (colI = 0; colI < size; colI++)
        {
          getBigIntegerFromLinked(basis[row][colI], &tmp5);
          BigInteger2Dec(&ptrDebugOutput, &tmp5, 0);
          *ptrDebugOutput = ',';
          ptrDebugOutput++;
          *ptrDebugOutput = ' ';
          ptrDebugOutput++;
        }
        *(ptrDebugOutput - 2) = '\r';
        *(ptrDebugOutput - 1) = '\n';
      }
    }
#endif
    colK++;
  } while (colK < size);
}

// Compute remainder such that the result is in range -divisor/2 to divisor/2.
static void BigIntSymmetricRemainder(const BigInteger* dividend,
  const BigInteger* divisor, BigInteger* result)
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
    ptrCoeff += numLimbs(ptrCoeff);
    ptrCoeff++;
  }
  // Compute coefficients of P.
  intToBigInteger(&operand1, 1);
#if DEBUG_VANHOEIJ
  copyStr(&ptrDebugOutput, "Coefficients: ");
#endif
  for (traceNbr = 0; traceNbr < polyDegree; traceNbr++)
  {
    getBigIntegerFromLinked(ptrCoeffs[traceNbr], &operand3);
#if DEBUG_VANHOEIJ
    BigInteger2Dec(&ptrDebugOutput, &operand3, 0);
    *ptrDebugOutput = ',';
    ptrDebugOutput++;
    *ptrDebugOutput = ' ';
    ptrDebugOutput++;
#endif
    (void)BigIntMultiply(&operand3, &operand1, &operand3);
    setLinkedBigInteger(&ptrCoeffs[traceNbr], &operand3);
  }
#if DEBUG_VANHOEIJ
  copyStr(&ptrDebugOutput, LF);
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
      if ((traceNbr - currDegree) <= polyDegree)
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

// Test whether the product divides the original polynomial.
// This means that f(0) divides F(0) where f(x) is the
// polynomial just computed and F(x) is the original
// polynomial times the leading coefficient.
static bool doStep1(int dividendMod32768[], int divisorMod32768[])
{
  int currentDegree;
  const int* ptrDest;

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
    ptrDest += numLimbs(ptrDest);
    ptrDest++;
    if ((currentDegree & 1) != 0)
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
    ptrDest += numLimbs(ptrDest);
    ptrDest++;
    if ((currentDegree & 1) != 0)
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
  // If the leading coefficient of the divisor is odd,
  // divide both polynomials mod 32768. If the remainder is
  // not zero, the integer division will not be performed.
  if ((divisorMod32768[(2 * poly5[0]) + 1] % 2) != 0)
  {
    enum eExprErr rc;
    const int* ptrMod32768;
    int TestNbr0Bak = TestNbr[0].x;
    TestNbr[0].x = 32768;
    modulusIsZero = false;
    NumberLength = 1;
    DividePolynomial(dividendMod32768, polyS[0], divisorMod32768, poly5[0], NULL);
    modulusIsZero = true;
    TestNbr[0].x = TestNbr0Bak;
    // Test whether the remainder is zero.
    ptrMod32768 = dividendMod32768;
    for (currentDegree = 0; currentDegree < poly5[0]; currentDegree++)
    {
      if (*(ptrMod32768 + 1) != 0)
      {
        return false;       // Coefficient of remainder is not zero. Go out.
      }
      ptrMod32768 += 2;     // Point to next coefficient.
    }
    rc = DivideIntegerPolynomial(polyS, poly5, TYPE_MODULUS);
    if (rc == EXPR_POLYNOMIAL_DIVISION_NOT_INTEGER)
    {                       // Cannot perform the division.
      return false;
    }
    if ((polyS[0] != 0) || (polyS[1] != 1) || (polyS[2] != 0))
    {                       // Remainder is not zero.
      return false;         // Number to factor does not divide this polynomial.
    }
  }
  return true;
}

static void doStep2(void)
{
  int degreePoly;
  // Get principal part of poly5 and store it to poly2.
  (void)getContent(poly5, &operand4);   // Content of polynomial.
  poly2[0] = poly5[0];
  const int *ptrSrc = &poly5[1];
  int *ptrDest = &poly2[1];
  for (int currentDegree = 0; currentDegree <= poly5[0]; currentDegree++)
  {
    UncompressBigIntegerB(ptrSrc, &operand2);
    (void)BigIntDivide(&operand2, &operand4, &operand3);
    NumberLength = operand3.nbrLimbs;
    BigInteger2IntArray(ptrDest, &operand3);
    ptrSrc += numLimbs(ptrSrc);
    ptrSrc++;
    ptrDest += numLimbs(ptrDest);
    ptrDest++;
  }
  // Copy this principal part to poly5.
  (void)CopyPolynomial(&poly5[1], &poly2[1], poly5[0]);
  (void)DivideIntegerPolynomial(polyNonRepeatedFactors, poly5, TYPE_DIVISION);
  degreePoly = poly5[0];
  (void)CopyPolynomial(ptrFactorInteger, &poly5[1], degreePoly);
  InsertIntegerPolynomialFactor(ptrFactorInteger, degreePoly);
}

static void forEachCurrentFactor(int stepNbr, int* pNbrFactors,
  int nbrVector, int currentFactor, int nbrVectors, int *pDegreeProd,
  struct sFactorInfo* pstFactorInfo)
{
  int degreeProd = *pDegreeProd;
  if (!linkedBigIntIsZero(lambda[nbrVector][currentFactor]))
  {                 // Multiply this polynomial.
    const int* ptrCoeffSrc = pstFactorInfo->ptrPolyLifted;    // Source
    int* ptrCoeffDest = poly2;                          // Destination
    int nbrLength;
    int degreeFactor = pstFactorInfo->degree;
    int currentDegree;
    // Reduce coefficients mod powerMod and store them on poly2.
    for (currentDegree = 0; currentDegree < degreeFactor; currentDegree++)
    {
      int lenBytes;
      nbrLength = numLimbs(ptrCoeffSrc);
      if (nbrVectors == 1)
      {              // Coefficient is already reduced.
        *ptrCoeffDest = nbrLength;
        lenBytes = nbrLength * (int)sizeof(int);
        (void)memcpy(ptrCoeffDest + 1, ptrCoeffSrc + 1, lenBytes);
      }
      else
      {              // Reduce the coefficient mod powerMod.
        operand1.nbrLimbs = nbrLength;
        operand1.sign = SIGN_POSITIVE;
        lenBytes = nbrLength * (int)sizeof(int);
        (void)memcpy(operand1.limbs, ptrCoeffSrc + 1, lenBytes);
        (void)BigIntRemainder(&operand1, &powerMod, &operand2);
        *ptrCoeffDest = operand2.nbrLimbs;
        lenBytes = operand2.nbrLimbs * (int)sizeof(int);
        (void)memcpy(ptrCoeffDest + 1, operand2.limbs, lenBytes);
      }
      ptrCoeffSrc += nbrLength;
      ptrCoeffSrc++;
      ptrCoeffDest += NumberLength;
      ptrCoeffDest++;
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
    ptrCoeffDest = poly1;               // Destination
    for (currentDegree = 0; currentDegree <= degreeProd; currentDegree++)
    {
      int lenBytes;
      nbrLength = 1 + numLimbs(ptrCoeffSrc);
      lenBytes = nbrLength * (int)sizeof(int);
      (void)memcpy(ptrCoeffDest, ptrCoeffSrc, lenBytes);
      ptrCoeffSrc += NumberLength;
      ptrCoeffSrc++;
      ptrCoeffDest += NumberLength;
      ptrCoeffDest++;
    }
    if (stepNbr == 2)
    {
      pstFactorInfo->multiplicity = 0;
      (*pNbrFactors)--;
    }
  }
  *pDegreeProd = degreeProd;
}

static bool AttemptToFactorStepNbr(int stepNbr, int nbrVectors,
                                   int nbrFactors, int* pNbrFactors)
{
  int currentDegree;
  int degreeProd;
  int* ptrMod32768;
  const int* ptrSrc;
  int* ptrDest;
  struct sFactorInfo* pstFactorInfo;
  int nbrTmp[1000];
  int dividendMod32768[(2 * MAX_DEGREE) + 1];
  int divisorMod32768[(2 * MAX_DEGREE) + 1];
#if DEBUG_VANHOEIJ
  copyStr(&ptrDebugOutput, "stepNbr = ");
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
      forEachCurrentFactor(stepNbr, pNbrFactors, nbrVector, currentFactor,
        nbrVectors, &degreeProd, pstFactorInfo);
      pstFactorInfo++;
    }
    // Multiply all coefficients by leadingCoeff by using modmult
    // (this converts from Montgomery to standard notation)
    // and store them in poly2.
    ptrSrc = poly1;
    ptrDest = poly2;
    CopyBigInt(&operand3, &leadingCoeff);
    if (operand3.nbrLimbs < NumberLength)
    {
      int lenBytes = (NumberLength - operand3.nbrLimbs) * (int)sizeof(limb);
      (void)memset(&operand3.limbs[operand3.nbrLimbs], 0, lenBytes);
    }
    CompressLimbsBigInteger((limb*)nbrTmp, &leadingCoeff);
    for (currentDegree = 0; currentDegree <= degreeProd; currentDegree++)
    {
      LenAndLimbs2ArrLimbs(ptrSrc, operand2.limbs, NumberLength);
      modmult(operand3.limbs, operand2.limbs, operand1.limbs);
      ArrLimbs2LenAndLimbs(ptrDest, operand1.limbs, NumberLength + 1);
      ptrSrc += NumberLength;
      ptrSrc++;
      ptrDest += numLimbs(ptrDest);
      ptrDest++;
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
      *ptrMod32768 = 1;
      ptrMod32768++;
      if (*ptrSrc >= 0)
      {
        *ptrMod32768 = *(ptrSrc + 1) & 32767;
      }
      else
      {
        *ptrMod32768 = (-*(ptrSrc + 1)) & 32767;
      }
      ptrMod32768++;
      ptrSrc += numLimbs(ptrSrc);
      ptrSrc++;
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
      *ptrMod32768 = 1;
      ptrMod32768++;
      if (operand1.sign == SIGN_POSITIVE)
      {
        *ptrMod32768 = operand1.limbs[0].x & 32767;
      }
      else
      {
        *ptrMod32768 = (-operand1.limbs[0].x) & 32767;
      }
      ptrMod32768++;
      (void)BigIntMultiply(&operand1, &leadingCoeff, &operand2);
      NumberLength = operand2.nbrLimbs;
      BigInteger2IntArray(ptrDest, &operand2);
      ptrSrc += numLimbs(ptrSrc);
      ptrSrc++;
      ptrDest += numLimbs(ptrDest);
      ptrDest++;
    }
    if (stepNbr == 1)
    {
      if (!doStep1(dividendMod32768, divisorMod32768))
      {
        return false;
      }
    }
    else
    {
      doStep2();
    }
    modulusIsZero = false;   // Perform modular operations.
    // Restart finding factors.
    (void)memset(arrNbrFactors, 0, sizeof(arrNbrFactors));
  }
  return true;
}

static bool AttemptToFactor(int nbrVectors, int nbrFactors, int *pNbrFactors)
{
  // In step 1 we check that all factors can be found. If this succeeds,
  // in step 2 we insert factors in final array of polynomial factors.
  if (!AttemptToFactorStepNbr(1, nbrVectors, nbrFactors, pNbrFactors))
  {
    return false;           // Cannot factor.
  }
  return AttemptToFactorStepNbr(2, nbrVectors, nbrFactors, pNbrFactors);
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
      ptrSrc += numLimbs(ptrSrc);
      ptrSrc++;
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
  {      // At this moment, |Am| >= |A0|
    CopyBigInt(&operand3, &operand2);
  }

  // Loop that finds the maximum value of bound for |Bj|.
  intToBigInteger(&operand2, 1);  // binomial(n-1, 0)
  // Set bound to |A0|
  UncompressBigIntegerB(&polyNonRepeatedFactors[1], &bound);
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
static void vanHoeij(int prime, int numFactors)
{
#if DEBUG_VANHOEIJ
  int nbrStepsDone = 0;
#endif
  int nbrFactors = numFactors;
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
  double b0 = log(2.0 * (double)values[0]) / logPrime;
  double log_rootbound = logBigNbr(&bound) / logPrime;
  double log_leadingcoeff = logBigNbr(&leadingCoeff) / logPrime;
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
  double dExponDifference = 6.0 * (double)nbrFactors * LOG_2 / log((double)prime);
  exponDifference = (int)dExponDifference;
  numberLLL = 0;
  b = (int)(b0 + (ceil(LOG_3 * log_rootbound) / logPrime) + 3.0);
  a0 = b + exponDifference;
  if (a0 <= (int)log_leadingcoeff + 2)
  {
    exponentMod = (int)log_leadingcoeff + 2;
  }
  else
  {
    exponentMod = a0;
  }
#if DEBUG_VANHOEIJ
  copyStr(&ptrDebugOutput, LF "====================================================="
    LF "prime = ");
  int2dec(&ptrDebugOutput, prime);
  copyStr(&ptrDebugOutput, "root bound = ");
  BigInteger2Dec(&ptrDebugOutput, &bound, 0);
  *ptrDebugOutput = '\n';
  ptrDebugOutput++;
#endif
  (void)memset(arrNbrFactors, 0, sizeof(arrNbrFactors));
  // Get leading coefficient of polyNonRepeatedFactors.
  ptrSrc = &polyNonRepeatedFactors[1];
  for (int degree1 = 0; degree1 < degreePolyToFactor; degree1++)
  {
    ptrSrc += numLimbs(ptrSrc);
    ptrSrc++;
  }
  computePower(exponentMod);
  modulusIsZero = false;    // Use modular arithmetic for polynomials.
  intToBigInteger(&operand5, 1);
  values[0] = polyNonRepeatedFactors[0];
  (void)getModPolynomial(&values[1], polyNonRepeatedFactors, &operand5);
#if DEBUG_HENSEL_LIFTING
  ptrOutput2 = ptrDebugOutput;
#endif
  (void)HenselLifting(factorInfoRecord, true);
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
      if ((elapsedTime / 10) != (oldTimeElapsed / 10))
      {
        char outputInfo[1000];
        char *ptrOut = outputInfo;

        oldTimeElapsed = elapsedTime;
        copyStr(&ptrOut, lang? "1<p>Obteniendo factores de dos factores modulares: prueba ":
                                "1<p>Finding factors from two modular factors: attempt ");
        int2dec(&ptrOut, currentAttempts);
        copyStr(&ptrOut, lang ? " de " : " of ");
        int2dec(&ptrOut, maxAttempts);
        copyStr(&ptrOut, "</p>");
        showElapsedTimeSec(&ptrOut);
        databack(outputInfo);
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
    else
    {           // Nothing to do.
    }
    pstFactorInfo++;
  }
  nbrFactors = newNbrFactors;
  nbrFactorsFound = nbrFactors;
  // Compute powerExtraBits as a power of prime with 3 times
  // the number of bits of nbrFactors (the number of polynomial factors).
  // Compute powerBoundA as powerExtraBits * powerMod.
  dExponDifference = 6.0 * (double)nbrFactors * LOG_2 / log((double)prime);
  exponDifference = (int)dExponDifference;
  intToBigInteger(&operand1, prime);
  (void)BigIntPowerIntExp(&operand1, exponDifference, &powerExtraBits);
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
        ((nbrRow == nbrCol)? 1 : 0));
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
            ((nbrRow == nbrCol)? 1 : 0));
        }
      }
    }
// Get b such that p^b is greater than the bounds on the cofficients.
    computePower(a0);
    newNumberLength = NumberLength;
    exponentMod = a0;
    intToBigInteger(&operand1, prime);
    (void)BigIntPowerIntExp(&operand1, b, &powerBoundA);
#if DEBUG_VANHOEIJ
    copyStr(&ptrDebugOutput, LF "====================================================="
      LF "prime = ");
    int2dec(&ptrDebugOutput, prime);
    copyStr(&ptrDebugOutput, ", a0 = ");
    int2dec(&ptrDebugOutput, a0);
    copyStr(&ptrDebugOutput, ", b0 = ");
    int2dec(&ptrDebugOutput, b);
    copyStr(&ptrDebugOutput, ", exponDifference = ");
    int2dec(&ptrDebugOutput, exponDifference);
    copyStr(&ptrDebugOutput, LF);
#endif
    intToBigInteger(&operand5, 1);
    values[0] = polyNonRepeatedFactors[0];
    (void)getModPolynomial(&values[1], polyNonRepeatedFactors, &operand5);
    NumberLength = newNumberLength;
    // use exponDifference additional bits instead of a fixed number

    intToBigInteger(&operand1, prime);
    (void)BigIntPowerIntExp(&operand1, exponDifference, &powerExtraBits);

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
      for (nbrCol = 0; nbrCol < (nbrVectors + nbrRequiredTraces); nbrCol++)
      {
        intToLinkedBigInt(&basisStar[nbrRow][nbrCol],
          ((nbrRow == nbrCol)? C : 0));
      }
    }

    // Use matrix lambda to hold the traces.
    for (int factorNbr = 0; factorNbr < nbrFactors; factorNbr++)
    {
      ComputeTraces(firstTrace + nbrRequiredTraces, factorNbr);
#if DEBUG_VANHOEIJ
      copyStr(&ptrDebugOutput, "Traces: ");
      for (nbrRow = firstTrace; nbrRow < firstTrace + nbrRequiredTraces; nbrRow++)
      {
        getBigIntegerFromLinked(traces[nbrRow], &tmp5);
        BigInteger2Dec(&ptrDebugOutput, &tmp5, 0);
        if (nbrRow < firstTrace + nbrRequiredTraces - 1)
        {
          *ptrDebugOutput = ',';
          ptrDebugOutput++;
          *ptrDebugOutput = ' ';
          ptrDebugOutput++;
        }
      }
      copyStr(&ptrDebugOutput, LF);
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

    for (nbrRow = nbrVectors; nbrRow < (nbrVectors + nbrRequiredTraces); nbrRow++)
    {
      for (nbrCol = nbrVectors; nbrCol < (nbrVectors + nbrRequiredTraces); nbrCol++)
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
    for (nbrRow = 0; nbrRow < (nbrVectors + nbrRequiredTraces); nbrRow++)
    {
      for (nbrCol = 0; nbrCol < (nbrVectors + nbrRequiredTraces); nbrCol++)
      {
        getBigIntegerFromLinked(basisStar[nbrRow][nbrCol], &tmp5);
        setLinkedBigInteger(&basis[nbrRow][nbrCol], &tmp5);
      }
    }
#if DEBUG_VANHOEIJ
    copyStr(&ptrDebugOutput, LF "Matrix M before LLL (");
    int2dec(&ptrDebugOutput, nbrVectors + nbrRequiredTraces);
    *ptrDebugOutput = '*';
    ptrDebugOutput++;
    int2dec(&ptrDebugOutput, nbrVectors + nbrRequiredTraces);
    copyStr(&ptrDebugOutput, "): ");
    for (nbrRow = 0; nbrRow < nbrVectors + nbrRequiredTraces; nbrRow++)
    {
      for (nbrCol = 0; nbrCol < nbrVectors + nbrRequiredTraces; nbrCol++)
      {
        getBigIntegerFromLinked(basisStar[nbrRow][nbrCol], &tmp5);
        BigInteger2Dec(&ptrDebugOutput, &tmp5, 0);
        if (nbrCol < nbrVectors + nbrRequiredTraces - 1)
        {
          *ptrDebugOutput = ',';
          ptrDebugOutput++;
        }
        intToLinkedBigInt(&lambda[nbrRow][nbrCol], 0);  // DEBUG BORRAR
      }
      *ptrDebugOutput = ';';
      ptrDebugOutput++;
    }
    *ptrDebugOutput = '.';
    ptrDebugOutput++;

#endif
    // Step 2: LLL-reduce the (r+s)*(r+s) matrix M (of rank r+s).
    integralLLL(nbrVectors + nbrRequiredTraces);

    // Copy matrix M to basisStar.
    for (nbrRow = 0; nbrRow < (nbrVectors + nbrRequiredTraces); nbrRow++)
    {
      for (nbrCol = 0; nbrCol < (nbrVectors + nbrRequiredTraces); nbrCol++)
      {
        getBigIntegerFromLinked(basis[nbrRow][nbrCol], &tmp5);
        setLinkedBigInteger(&basisStar[nbrRow][nbrCol], &tmp5);
      }
    }
#if DEBUG_VANHOEIJ
    copyStr(&ptrDebugOutput, LF "Matrix M after LLL: " LF);
    for (nbrRow = 0; nbrRow < nbrVectors + nbrRequiredTraces; nbrRow++)
    {
      for (nbrCol = 0; nbrCol < nbrVectors + nbrRequiredTraces; nbrCol++)
      {
        getBigIntegerFromLinked(basisStar[nbrRow][nbrCol], &tmp5);
        BigInteger2Dec(&ptrDebugOutput, &tmp5, 0);
        *ptrDebugOutput = ' ';
        ptrDebugOutput++;
      }
      *ptrDebugOutput = ';';
      ptrDebugOutput++;
    }
    *ptrDebugOutput = '.';
    ptrDebugOutput++;
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
      for (nbrCol = 0; nbrCol < (nbrVectors + nbrRequiredTraces); nbrCol++)
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
      for (nbrCol = 0; nbrCol < (nbrVectors + nbrRequiredTraces); nbrCol++)
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
      for (nbrCol = 0; nbrCol < (nbrVectors + nbrRequiredTraces); nbrCol++)
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
    copyStr(&ptrDebugOutput, LF "Matrix M before Gram-Schmidt" LF);
    for (nbrRow = 0; nbrRow < nbrVectors + nbrRequiredTraces; nbrRow++)
    {
      for (nbrCol = 0; nbrCol < nbrVectors + nbrRequiredTraces; nbrCol++)
      {
        getBigIntegerFromLinked(basisStar[nbrRow][nbrCol], &tmp5);
        BigInteger2Dec(&ptrDebugOutput, &tmp5, 0);
        *ptrDebugOutput = ' ';
        ptrDebugOutput++;
      }
      *ptrDebugOutput = ';';
      ptrDebugOutput++;
    }
    *ptrDebugOutput = '.';
    ptrDebugOutput++;
    *ptrDebugOutput = 0;
#endif
    // Step 4: Perform Gram-Schmidt orthogonalization on M.
    // Let c_i, 1 <= i <= r+s be the norm of the ith 
    // orthogonal vector.
    GramSchmidtOrthogonalization(nbrFactors + nbrRequiredTraces, nbrVectors + nbrRequiredTraces);

#if DEBUG_VANHOEIJ
    copyStr(&ptrDebugOutput, LF "Norms after Gram-Schmidt: ");
    for (r1 = 2; r1 <= nbrVectors + nbrRequiredTraces; r1++)
    {
      getBigIntegerFromLinked(lambda[r1 - 1][r1 - 1], &tmp4);
      getBigIntegerFromLinked(lambda[r1 - 2][r1 - 2], &tmp5);
      (void)BigIntDivide(&tmp4, &tmp5, &operand1);
      BigInteger2Dec(&ptrDebugOutput, &operand1, 0);
      *ptrDebugOutput = ' ';
      ptrDebugOutput++;
    }
    *ptrDebugOutput = '.';
    ptrDebugOutput++;
    *ptrDebugOutput = 0;
#endif
    // Step 5: Let r' the largest value such that all c_i of
    // index larger than r' are of norm greater than
    // sqrt(c^2*n_0 + s*n_0^2/4)

    // Compute square of formula.
    squareFormula = (nbrFactors * C * C) + (nbrRequiredTraces * nbrFactors * nbrFactors / 4);
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
    copyStr(&ptrDebugOutput, LF "squareFormula = ");
    int2dec(&ptrDebugOutput, squareFormula);
#endif

    // Step 6: If r' = 1, return "irreducible".
    if (r1 == 1)
    {                 // Polynomial is irreducible.
#if DEBUG_VANHOEIJ
      copyStr(&ptrDebugOutput, LF "r' = 1 -> irreducible." LF);
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
      char* ptrOutputInfo = outputInfo;
      if (lang)
      {
        copyStr(&ptrOutputInfo, "1<p>Paso 7</p>");
      }
      else
      {
        copyStr(&ptrOutputInfo, "1<p>Step 7</p>");
      }
      databack(outputInfo);
    }
#endif
#endif
    // Step 7: BL <- BL / C. Now matrix BL has dimension 
    // nbrFactors rows by r' columns.
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
    copyStr(&ptrDebugOutput, LF "Matrix BL before Gauss: ");
    for (nbrRow = 0; nbrRow < nbrFactors; nbrRow++)
    {
      for (nbrCol = 0; nbrCol < r1; nbrCol++)
      {
        getBigIntegerFromLinked(basisStar[nbrRow][nbrCol], &tmp5);
        BigInteger2Dec(&ptrDebugOutput, &tmp5, 0);
        *ptrDebugOutput = ' ';
        ptrDebugOutput++;
      }
      *ptrDebugOutput = ';';
      ptrDebugOutput++;
    }
    *ptrDebugOutput = '.';
    ptrDebugOutput++;
    *ptrDebugOutput = 0;
#ifdef __EMSCRIPTEN__
    {
      char outputInfo[1000];
      char* ptrOutputInfo = outputInfo;
      if (lang)
      {
        copyStr(&ptrOutputInfo, "1<p>Paso 8</p>");
      }
      else
      {
        copyStr(&ptrOutputInfo, "1<p>Step 8</p>");
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
    copyStr(&ptrDebugOutput, LF "squareFormula = ");
    int2dec(&ptrDebugOutput, squareFormula);
    copyStr(&ptrDebugOutput, ", rank = ");
    int2dec(&ptrDebugOutput, rank);
    copyStr(&ptrDebugOutput, LF);
#endif
    nbrRequiredTraces++;
    nbrVectors = r1;
#if 0
    nbsStepsDone++;
    if (nbrStepsDone == 10)
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
    copyStr(&ptrDebugOutput, LF "Matrix R after Gauss: ");
    for (nbrRow = 0; nbrRow < nbrFactors; nbrRow++)
    {
      for (nbrCol = 0; nbrCol < r1; nbrCol++)
      {
        getBigIntegerFromLinked(lambda[nbrRow][nbrCol], &tmp5);
        BigInteger2Dec(&ptrDebugOutput, &tmp5, 0);
        *ptrDebugOutput = ' ';
        ptrDebugOutput++;
      }
      *ptrDebugOutput = ';';
      ptrDebugOutput++;
    }
    *ptrDebugOutput = '.';
    ptrDebugOutput++;
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
      copyStr(&ptrDebugOutput, "nbrCol = ");
      int2dec(&ptrDebugOutput, nbrCol);
      copyStr(&ptrDebugOutput, ", nbrFactors = ");
      int2dec(&ptrDebugOutput, nbrFactors);
      copyStr(&ptrDebugOutput, "They should have been equal." LF);
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
    char* ptrOutputInfo = outputInfo;
    if (lang)
    {
      copyStr(&ptrOutputInfo, "1<p>Verificando si los polinomios hallados son factores irreducibles.</p>");
    }
    else
    {
      copyStr(&ptrOutputInfo, "1<p>Testing whether the polynomials found are irreducible factors.</p>");
    }
    showElapsedTimeSec(&ptrOutputInfo);
    databack(outputInfo);
#endif
    if (AttemptToFactor(nbrVectors, nbrFactors, &newNbrFactors))
    {
#if DEBUG_VANHOEIJ
      copyStr(&ptrDebugOutput, "nbrVector == nbrVectors" LF);
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
    ptrSrc += *ptrSrc;
    ptrSrc++;
    // If operand1 >= halfPowerMod, subtract powerMod.
    BigIntSubt(&operand1, &halfPowerMod, &operand2);
    if (operand2.sign == SIGN_POSITIVE)
    {
      BigIntSubt(&operand1, &powerMod, &operand1);
    }
    BigInteger2IntArray(ptrDest, &operand1);
    ptrDest += numLimbs(ptrDest);
    ptrDest++;
  }
}

// Insert new factor in factorInfoInteger array sorting this array
// on ascending order of degree, and then on ascending order of
// leading coefficients. Merge equal factors.
static void InsertIntegerPolynomialFactor(int* ptrFactor, int degreePoly)
{
  struct sFactorInfo* pstFactorInfoInteger;
  const struct sFactorInfo* pstFactorInfo;
  int indexNewFactor[MAX_DEGREE+2];
  int indexOldFactor[MAX_DEGREE+2];
  int currentDegree;
  int index;
  int* ptrIndex;
  const int* ptrOldFactor;

  if ((degreePoly < 0) || (degreePoly > MAX_DEGREE))
  {    // Invalid values of degreePoly.
    return;
  }
  // Fill indexes to start of each coefficient.
  ptrIndex = &indexNewFactor[0];
  index = 0;
  for (currentDegree = 0; currentDegree <= (degreePoly+1); currentDegree++)
  {
    *ptrIndex = index;
    ptrIndex++;
    index += numLimbs(ptrFactor + index) + 1;
  }
  pstFactorInfoInteger = factorInfoInteger;
  while ((pstFactorInfoInteger->ptrPolyLifted != NULL) &&
         (degreePoly >= pstFactorInfoInteger->degree))
  {      // Polynomial factors are sorted by degree and then
         // by coefficients.
    if (degreePoly > pstFactorInfoInteger->degree)
    {    // New factor degree is greater than stored factor degree.
      pstFactorInfoInteger++;
      continue;  // Check next stored factor.
    }
    // Degree of the new factor is the same as the one already stored.
    ptrOldFactor = pstFactorInfoInteger->ptrPolyLifted;
    ptrIndex = &indexOldFactor[0];
    index = 0;
    for (currentDegree = 0; currentDegree <= degreePoly; currentDegree++)
    {   // Loop that gets the pointers to coefficients of stored factor.
      *ptrIndex = index;
      ptrIndex++;
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
    pstFactorInfoInteger++;
  }
  // Move elements of array factorInfoInteger to make room for new factor.
  pstFactorInfo = pstFactorInfoInteger;
  while (pstFactorInfo->ptrPolyLifted != NULL)
  {
    pstFactorInfo++;
  }
  if (pstFactorInfo > pstFactorInfoInteger)
  {
    (void)memmove(pstFactorInfoInteger + 1, pstFactorInfoInteger,
      (pstFactorInfo - pstFactorInfoInteger) * sizeof(*pstFactorInfo));
  }
  pstFactorInfoInteger->multiplicity = intPolyMultiplicity;
  pstFactorInfoInteger->ptrPolyLifted = ptrFactor;
  pstFactorInfoInteger->degree = degreePoly;
  ptrFactorInteger += indexNewFactor[degreePoly + 1];
}

int getNextPrimeNoDuplicatedFactors(int primeIndex)
{
  const int* ptrSrc;
  int primeIdx = primeIndex;
  int degreeGcdMod;
  int polyDegree = polyNonRepeatedFactors[0];
  initializeSmallPrimes(smallPrimes);
  // Get leading coefficient of polyNonRepeatedFactors.
  ptrSrc = &polyNonRepeatedFactors[1];
  for (int currentDegree = 0; currentDegree < polyDegree; currentDegree++)
  {
    ptrSrc += numLimbs(ptrSrc);
    ptrSrc++;
  }
  modulusIsZero = true;
  UncompressBigIntegerB(ptrSrc, &leadingCoeff);
  modulusIsZero = false;
  do
  {   // Loop that finds a prime modulus such that the factorization
      // has no repeated factors. That means gcd(F, F') = 1 (mod prime).
      // This is required because Hensel lift does not work when
      // repeated factors are present.
      // If the leading coefficient is multiple of prime, the prime cannot
      // be used.
    int degree1;
    int degree2;
    int prime;
    do
    {      // Loop while the leading coefficient is multiple of prime.
      primeIdx++;
      prime = smallPrimes[primeIdx];
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
  return primeIdx;
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
  (void)CopyPolynomial(&polyBackup[1], &values[1], (values[0] >= 0)? values[0] : 0);
  values[0] = degreePolyToFactor;
  (void)CopyPolynomial(&values[1], &poly1[1], degreePolyToFactor);
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
  // Input is not in Montgomery notation.
  (void)FactorModularPolynomial(false);
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
    *ptrPolyLiftedRecord = 1;    // Leading coefficient should be 1.
    ptrPolyLiftedRecord++;
    *ptrPolyLiftedRecord = 1;
    ptrPolyLiftedRecord++;
    pstFactorInfoOrig++;
    pstFactorInfoRecord++;
  }
  pstFactorInfoRecord->ptr = NULL; // Indicate that there are no more factors.
}

// Apply Yun's algorithm for integer squarefree factorization.
// Let f = (a_1) * (a_2)^2 * (a_3)^3 * ... * (a_k)^k (all polynomials).
//
// b <- f, d <- f', a <- gcd(b, d)
// repeat
//   b <- b/a, c <- d/a, d <- c - b', a <- gcd(b, d)
//   output a
// until a = b.
// f = PolyToFactor.
// Store polynomial factors in array pointed by pstFactorInfoInteger.
// poly1 = a
// poly2 = b
// poly3 = c
// poly4 = d

static int IntegerSquarefreeFactorization(void)
{
  int* ptrPolySqFreeFact = polySqFreeFact;
  int multiplicity = 1;
  int nbrFactors = 0;
  polyB[0] = polyToFactor[0];
  (void)CopyPolynomial(&polyB[1], &polyToFactor[1], polyToFactor[0]); // b <- f
  polyD[0] = polyToFactor[0];
  (void)CopyPolynomial(&polyD[1], &polyToFactor[1], polyToFactor[0]); // d <- f
  DerPolynomial(polyD);                                            // d <- f'
  PolynomialGcd(polyB, polyD, polyA);                              // a <- gcd(b, d)
  (void)DivideIntegerPolynomial(polyB, polyA, TYPE_DIVISION);      // b <- b/a
  do
  {
    polyC[0] = polyD[0];
    (void)CopyPolynomial(&polyC[1], &polyD[1], polyD[0]);          // c <- d
    (void)DivideIntegerPolynomial(polyC, polyA, TYPE_DIVISION);    // c <- d/a
    tempPoly[0] = polyB[0];
    (void)CopyPolynomial(&tempPoly[1], &polyB[1], polyB[0]);       // temp <- b
    DerPolynomial(tempPoly);                                       // temp <- b'
    SubtractIntegerPolynomial(polyC, tempPoly, polyD);             // d <- c - b'
    PolynomialGcd(polyB, polyD, polyA);                            // a <- gcd(b, d)
    if (polyA[0] != 0)
    {
      // Copy polynomial polyA to factor array if its degree is not zero.
      *ptrPolySqFreeFact = multiplicity;
      ptrPolySqFreeFact++;
      *ptrPolySqFreeFact = polyA[0];
      ptrPolySqFreeFact++;
      ptrPolySqFreeFact = CopyPolynomial(ptrPolySqFreeFact, &polyA[1], polyA[0]);
      nbrFactors++;
    }
    multiplicity++;
    (void)DivideIntegerPolynomial(polyB, polyA, TYPE_DIVISION);    // b <- b/a
    // Continue loop if b does not equal 1.
  } while ((polyB[0] != 0) || (polyB[1] != 1) || (polyB[2] != 1));
  return nbrFactors;
}

// Up to 5 different prime moduli are tested to minimize the number
// of factors because the Van Hoeij algorithm speed depends on
// this number. If the number of factors is less than 10,
// no more modular factorizations are attempted.
static void factorDifferentModuli(int *pNbrFactorsRecord, int *pPrimeRecord)
{
  int primeRecord = 0;
  const struct sFactorInfo* pstFactorInfoOrig;
  // The big number ensures that the first copy is done.
  int nbrFactorsRecord = 100000;
  int primeIndex = 0;
  for (int attemptNbr = 1; attemptNbr < 5; attemptNbr++)
  {
    int expon = 1;
    int prime;
    int nbrFactors;
    bool isIrreducible;
    primeIndex = getNextPrimeNoDuplicatedFactors(primeIndex);
    prime = smallPrimes[primeIndex];
    // Find expon such that prime^expon >= 2 * bound
    BigIntMultiplyBy2(&bound);
    intToBigInteger(&operand1, prime);
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
    (void)memset(validDegrees, 0x00, sizeof(validDegrees));
    validDegrees[0] = 1;
    int curDegree = 0;
    for (int factorNbr = 0; factorNbr < MAX_DEGREE; factorNbr++)
    {
      if (pstFactorInfoOrig->ptr == NULL)
      {    // No more factors.
        break;
      }
      // Compute D <- D OR (D << d) where D is validDegrees and
      // d is the degree of the factor found.
      int nbrFactorsSameDegree = pstFactorInfoOrig->degree / pstFactorInfoOrig->expectedDegree;
      int shiftLeftCtr = pstFactorInfoOrig->expectedDegree;
      int shiftLeftLimbs = shiftLeftCtr / 32;
      unsigned int shiftLeftRem = (unsigned int)shiftLeftCtr & 0x1FU;
      if ((shiftLeftCtr & 0x1F) == 0)
      {
        for (int curFactor = 0; curFactor < nbrFactorsSameDegree; curFactor++)
        {
          for (int limbIndex = curDegree / 32; limbIndex >= 0; limbIndex--)
          {
            validDegrees[limbIndex + shiftLeftLimbs] |= validDegrees[limbIndex];
          }
          curDegree += shiftLeftCtr;
        }
      }
      else
      {
        for (int curFactor = 0; curFactor < nbrFactorsSameDegree; curFactor++)
        {
          int limbIndex = curDegree / 32;
          unsigned int shiftRightRem = 32U - shiftLeftRem;
          validDegrees[limbIndex + shiftLeftLimbs + 1] |=
            (validDegrees[limbIndex] >> shiftRightRem);
          for (; limbIndex > 0; limbIndex--)
          {
            validDegrees[limbIndex + shiftLeftLimbs] |=
              ((validDegrees[limbIndex] << shiftLeftRem) |
                (validDegrees[limbIndex - 1] >> shiftRightRem));
          }
          validDegrees[shiftLeftLimbs] |= (validDegrees[0] << shiftLeftRem);
          curDegree += shiftLeftCtr;
        }
      }
      nbrFactors += nbrFactorsSameDegree;
      pstFactorInfoOrig++;
    }
    // Update valid degrees record by performing a logic AND with valid degrees.
    for (int limbIndex = curDegree / 32; limbIndex >= 0; limbIndex--)
    {
      validDegreesRecord[limbIndex] &= validDegrees[limbIndex];
    }
    // Test whether the polynomial is irreducible.
    isIrreducible = true;
    if (curDegree < 32)
    {
      if (validDegreesRecord[0] != (1U + (1U << (unsigned int)curDegree)))
      {
        isIrreducible = false;
      }
    }
    else
    {
      if ((validDegreesRecord[0] != 1U) ||
        (validDegreesRecord[curDegree >> 5] != (1U << (curDegree & 0x1F))))
      {
        isIrreducible = false;
      }
    }
    if (isIrreducible)
    {
      for (int limbIndex = (curDegree / 32) - 1; limbIndex > 0; limbIndex--)
      {
        if (validDegreesRecord[limbIndex] != 0U)
        {
          isIrreducible = false;
          break;
        }
      }
    }
    if (isIrreducible)
    {
      nbrFactorsRecord = 1;
      break;
    }
    if (nbrFactors < nbrFactorsRecord)
    {    // Copy factors found to records arrays.
      primeRecord = prime;
      nbrFactorsRecord = nbrFactors;
      CopyFactorsFoundToRecord();
    }
  }
  *pNbrFactorsRecord = nbrFactorsRecord;
  *pPrimeRecord = primeRecord;
}

// Input: values = degree, coefficient degree 0, coefficient degree 1, etc.
// Output: factorInfo = structure that holds the factors.
int FactorPolyOverIntegers(void)
{
  const int* ptrPolySqFreeFact = polySqFreeFact;
  int degreePolyToFactor = values[0];
  int primeRecord = 0;
  const int* ptrSrc;
  int* ptrDest;
  int factorNbr;
  bool polXprocessed = false;
  int* ptrFactorIntegerBak;
  int* ptrPolyLiftedOrig;
  int nbrSquareFreeFactors;
  struct sFactorInfo* pstFactorInfoOrig;
  struct sFactorInfo* pstFactorInfoInteger = factorInfoInteger;
  nbrFactorsFound = 0;
  initLinkedBigInt();
  ptrFactorInteger = polyInteger;
  modulusIsZero = true;
  (void)memset(factorInfoInteger, 0, sizeof(factorInfoInteger));
  (void)getContent(values, &contentPolyToFactor);
  (void)CopyPolynomial(&origPolyToFactor[1], &values[1], degreePolyToFactor);
  origPolyToFactor[0] = degreePolyToFactor;
  // polyToFactor -> original polynomial / content of polynomial.
  // Let n be the degree of the least coefficient different from zero.
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
        ptrSrc += numLimbs(ptrSrc);
        ptrSrc++;
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
        nbrFactorsFound = 1;
        polyToFactor[0] = degreePolyToFactor - currentDegree;
      }
      polXprocessed = true;
    }
    (void)BigIntDivide(&operand1, &contentPolyToFactor, &operand2);
    NumberLength = operand2.nbrLimbs;
    BigInteger2IntArray(ptrDest, &operand2);
    ptrSrc += numLimbs(ptrSrc);
    ptrSrc++;
    ptrDest += numLimbs(ptrDest);
    ptrDest++;
  }
  CopyBigInt(&operand5, &contentPolyToFactor);
  if (polyToFactor[0] == 0)
  { // Degree of polynomial is zero.
    nbrSquareFreeFactors = 0;
  }
  else
  {
    // Use Eisenstein criterion to detect some irreducible polynomials.
    if ((polyToFactor[0] > 1) && EisensteinCriterion(polyToFactor))
    { // Polynomial of degree > 1 is irreducible. Do not attempt to factor it.
      nbrFactorsFound++;
      (void)CopyPolynomial(ptrFactorInteger, &polyToFactor[1],
        polyToFactor[0]);
      pstFactorInfoInteger->degree = polyToFactor[0];
      pstFactorInfoInteger->expectedDegree = polyToFactor[0];
      pstFactorInfoInteger->multiplicity = 1;
      pstFactorInfoInteger->ptrPolyLifted = ptrFactorInteger;
      if ((pstFactorInfoInteger != &factorInfoInteger[0]) ||
        !BigIntIsOne(&contentPolyToFactor))
      {
        primeEisenstein = 0;   // Do not show irreducible polynomial.
      }
      return EXPR_OK;
    }
    nbrSquareFreeFactors = IntegerSquarefreeFactorization();
  }
  (void)memset(validDegreesRecord, 0xFF, sizeof(validDegreesRecord));
  for (int squareFreeFactor = 0; squareFreeFactor < nbrSquareFreeFactors; squareFreeFactor++)
  {    // At least degree 1.
       // The trailing coefficient of factors must divide the product of the trailing and leading
       // coefficients of the original polynomial.
    int nbrFactorsRecord;
    modulusIsZero = true;
    intPolyMultiplicity = *ptrPolySqFreeFact;
    ptrPolySqFreeFact++;
    polyNonRepeatedFactors[0] = *ptrPolySqFreeFact;
    ptrPolySqFreeFact++;
    ptrDest = &polyNonRepeatedFactors[1];
    for (int currentDegree = 0; currentDegree <= polyNonRepeatedFactors[0]; currentDegree++)
    {         // Copy polynomial.
      int numLength = numLimbs(ptrPolySqFreeFact) + 1;
      int lenBytes = numLength * (int)sizeof(int);
      (void)memcpy(ptrDest, ptrPolySqFreeFact, lenBytes);
      ptrPolySqFreeFact += numLength;
      ptrDest += numLength;
    }
    // Get trailing coefficient.
    ptrSrc = &polyNonRepeatedFactors[1];
    UncompressBigIntegerB(ptrSrc, &operand1);
    // Get leading coefficient.
    for (int degree1 = 0; degree1 < polyNonRepeatedFactors[0]; degree1++)
    {
      ptrSrc += numLimbs(ptrSrc);
      ptrSrc++;
    }
    UncompressBigIntegerB(ptrSrc, &operand2);
    (void)BigIntMultiply(&operand1, &operand2, &trailingCoeff);
    ComputeCoeffBounds();   // bound = Bound of coefficient of factors.
    factorDifferentModuli(&nbrFactorsRecord, &primeRecord);
    // Modulus that generate the lowest number of polynomials factors found.
    // Perform same degree factorization.
    // Copy back the record factorization to the work area.
    if (nbrFactorsRecord > 1)
    {
      int prime;
      const struct sFactorInfo* pstFactorInfoRecord = factorInfoRecord;
      pstFactorInfoOrig = factorInfo;
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
      vanHoeij(prime, nbrFactorsRecord);
    }
    // Polynomial is irreducible.
    if (polyNonRepeatedFactors[0] > 0)
    {    // Degree is greater than zero. Copy it to integer polynomial factor array.
      ptrFactorIntegerBak = CopyPolynomial(ptrFactorInteger, &polyNonRepeatedFactors[1],
        polyNonRepeatedFactors[0]);
      InsertIntegerPolynomialFactor(ptrFactorInteger, polyNonRepeatedFactors[0]);
      ptrFactorInteger = ptrFactorIntegerBak;
      pstFactorInfoInteger++;
    }
  }
  modulusIsZero = true;
  (void)CopyPolynomial(&values[1], &origPolyToFactor[1], origPolyToFactor[0]);
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

// Use Eisenstein criterion to determine whether the polynomial poly
// is irreducible or not. It returns true if it is irreducible.
// If this function returns false, it could be irreducible or not.
// Let f the leading coefficient of the polynomial.
// Let t the trailing coefficient of the polynomial.
// Let a = gcd of all coefficients except the leading and trailing.
// Let c = gcd(a, t)
// Let b = c / gcd(c, f)
// For each prime factor of b (name it b_i): if t mod b_i^2 != 0,
// the polynomial is irreducible, so go out.
// Let c = gcd(a, f)
// Let b = c / gcd(c, t)
// For each prime factor of b (name it b_i): if f mod b_i^2 != 0,
// the polynomial is irreducible, so go out.
// Do not process B if it does not fit in a limb.
static bool EisensteinCriterion(const int* poly)
{
  int polyDegree = *poly;
  const int *ptrPoly = poly + 1;
  const int* ptrLeading;
  const int* ptrTrailing = ptrPoly;
  primeEisenstein = 0;
  ptrPoly += numLimbs(ptrPoly);
  ptrPoly++;
  intToBigInteger(&operand1, 0);  // Initialize a.
  for (int currentDegree = 1; currentDegree < polyDegree; currentDegree++)
  {                               // Loop that computes a.
    UncompressBigIntegerB(ptrPoly, &operand2);
    BigIntGcd(&operand1, &operand2, &operand3);
    CopyBigInt(&operand1, &operand3);
    ptrPoly += numLimbs(ptrPoly);
    ptrPoly++;
  }
  ptrLeading = ptrPoly;
  if (checkEisenstein(&operand1, ptrLeading, ptrTrailing))
  {
    return true;
  }
  return checkEisenstein(&operand1, ptrTrailing, ptrLeading);
}

static bool checkEisenstein(const BigInteger *gcdAll, 
  const int *ptrLeading, const int *ptrTrailing)
{
  int B;
  int prime;
  UncompressBigIntegerB(ptrTrailing, &operand4);
  BigIntGcd(gcdAll, &operand4, &operand2);        // c <- operand2
  UncompressBigIntegerB(ptrLeading, &operand4);
  BigIntGcd(&operand2, &operand4, &operand3);
  (void)BigIntDivide(&operand2, &operand3, &operand4);  // b <- operand4
  if (operand4.nbrLimbs > 1)
  {
    return false;    // Cannot process. Number too big.
  }
  UncompressBigIntegerB(ptrLeading, &operand3);
  B = operand4.limbs[0].x;
  // Find prime factors of B.
  if ((B % 2) == 0)
  {
    if (((*(ptrLeading + 1) % 4) != 0) && ((*(ptrTrailing + 1) % 4) != 0))
    {
      primeEisenstein = 2;
      return true;   // Polynomial is irreducible.
    }
    while ((B % 2) == 0)
    {
      B /= 2;
    }
  }
  prime = 3;
  while ((prime * prime) <= B)
  {
    if ((B % prime) != 0)
    {
      prime += 2;
      continue;
    }
    if (getRemainder(&operand3, prime) != 0)
    {
      UncompressBigIntegerB(ptrTrailing, &operand2);
      subtractdivide(&operand2, 0, prime);
      if (getRemainder(&operand2, prime) != 0)
      {
        primeEisenstein = prime;
        return true;   // Polynomial is irreducible.
      }
    }
    while ((B % prime) == 0)
    {
      B /= prime;
    }
    prime += 2;
  }
  if ((B > 1) &&    // At this moment B is prime.
     getRemainder(&operand3, B) != 0)
  {
    UncompressBigIntegerB(ptrTrailing, &operand2);
    subtractdivide(&operand2, 0, B);
    if (getRemainder(&operand2, B) != 0)
    {
      primeEisenstein = B;
      return true;   // Polynomial is irreducible.
    }
  }
  return false;
}

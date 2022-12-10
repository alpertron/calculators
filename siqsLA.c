//
// This file is part of Alpertron Calculators.
//
// Copyright 2016-2021 Dario Alejandro Alpern
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
#include <math.h>
#include <stdint.h>
#include "bignbr.h"
#include "expression.h"
#include "factor.h"
#include "commonstruc.h"
#if (DEBUG_SIQS == 2) && !defined(__EMSCRIPTEN__)
#include <stdio.h>
#endif

#define processCol(col) (*(RightMatr + col) & ((leftMatr << col) >> 31))
#ifdef __EMSCRIPTEN__
extern char lowerText[MAX_LEN * 16];
extern char* ptrLowerText;
unsigned char SIQSInfoText[300];
#endif

#ifdef __EMSCRIPTEN__
static void showMatrixSize(const char* SIQSInformationText, int rows, int cols)
{
  (void)SIQSInformationText;
  char* ptrText = ptrLowerText;  // Point after number that is being factored.
  copyStr(&ptrText, lang ? "<p>Resolviendo la matriz de congruencias de " : "<p>Solving ");
  int2dec(&ptrText, rows);   // Show number of rows.
  copyStr(&ptrText, " &times; ");
  int2dec(&ptrText, cols);   // Show number of columns.
  copyStr(&ptrText, lang ? " usando el algoritmo de Lanczos en bloques.</p>" :
    " congruence matrix using Block Lanczos algorithm.</p>");
  databack(lowerText);
}
#endif
/* Multiply binary matrices of length m x 32 by 32 x 32 */
/* The product matrix has size m x 32. Then add it to a m x 32 matrix. */
static void MatrixMultAdd(const int *LeftMatr, const int *RightMatr, int *ProdMatr)
{
  int matrLength = common.siqs.matrixBLength;
  for (int row = 0; row < matrLength; row++)
  {
    int leftMatr = *(LeftMatr+row);
    // For each column...
    *(ProdMatr + row) ^=
    processCol(0) ^
    processCol(1) ^
    processCol(2) ^
    processCol(3) ^
    processCol(4) ^
    processCol(5) ^
    processCol(6) ^
    processCol(7) ^
    processCol(8) ^
    processCol(9) ^
    processCol(10) ^
    processCol(11) ^
    processCol(12) ^
    processCol(13) ^
    processCol(14) ^
    processCol(15) ^
    processCol(16) ^
    processCol(17) ^
    processCol(18) ^
    processCol(19) ^
    processCol(20) ^
    processCol(21) ^
    processCol(22) ^
    processCol(23) ^
    processCol(24) ^
    processCol(25) ^
    processCol(26) ^
    processCol(27) ^
    processCol(28) ^
    processCol(29) ^
    processCol(30) ^
    processCol(31);
  }
}
/* Multiply binary matrices of length m x 32 by 32 x 32 */
/* The product matrix has size m x 32 */
static void MatrixMultiplication(const int *LeftMatr, const int *RightMatr, int *ProdMatr)
{
  int matrLength = 32;
  for (int row = 0; row < matrLength; row++)
  {
    int prodMatr = 0;
    int leftMatr = *(LeftMatr+row);
    for (int col=0; col < 32; col++)
    {
      if (leftMatr < 0)
      {
        prodMatr ^= *(RightMatr+col);
      }
      leftMatr *= 2;
    }
    *(ProdMatr+row) = prodMatr;
  }
}

/* Multiply the transpose of a binary matrix of length n x 32 by */
/* another binary matrix of length n x 32 */
/* The product matrix has size 32 x 32 */
static void MatrTranspMult(int matrLength, const int *LeftMatr, const int *RightMatr, int *ProdMatr)
{
  for (int col = 31; col >= 0; col--)
  {
    int prodMatr = 0;
    for (int row = 0; row < matrLength; row++)
    {
      prodMatr ^= *(RightMatr + row) & ((*(LeftMatr + row) << col) >> 31);
    }
    *(ProdMatr+col) = prodMatr;
  }
}

static void MatrixAddition(const int *leftMatr, const int *rightMatr, int *sumMatr)
{
  for (int row = 32 - 1; row >= 0; row--)
  {
    *(sumMatr+row) = *(leftMatr+row) ^ *(rightMatr+row);
  }
}

static void MatrMultBySSt(int length, const int *Matr, int diagS, int *Prod)
{
  for (int row = length - 1; row >= 0; row--)
  {
    *(Prod+row) = diagS & *(Matr+row);
  }
}

/* Compute Bt * B * input matrix where B is the matrix that holds the */
/* factorization relations */
static void MultiplyAByMatrix(const int *Matr, int *TempMatr, int *ProdMatr)
{
  int index;
  int row;
  const int *rowMatrixB;

  /* Compute TempMatr = B * Matr */
  (void)memset(TempMatr, 0, common.siqs.matrixBLength*sizeof(int));
  for (row = common.siqs.matrixBLength - 1; row >= 0; row--)
  {
    int rowValue;
    rowMatrixB = common.siqs.matrixB[row];
    rowValue = *(Matr+row);
    for (index = *(rowMatrixB+LENGTH_OFFSET)-1; index >= 1; index--)
    {
      *(TempMatr+*(rowMatrixB+index)) ^= rowValue;
    }
  }

  /* Compute ProdMatr = Bt * TempMatr */
  for (row = common.siqs.matrixBLength - 1; row >= 0; row--)
  {
    int prodMatr = 0;
    rowMatrixB = common.siqs.matrixB[row];
    for (index = *(rowMatrixB+LENGTH_OFFSET)-1; index >= 1; index--)
    {
      prodMatr ^= *(TempMatr+*(rowMatrixB+index));
    }
    *(ProdMatr+row) = prodMatr;
  }
}

static void colexchange(int* XmY, int* V, int* V1, int* V2,
  int col1, int col2)
{
  int* matr1;
  int* matr2;
  int* matr3;
  int* matr4;

  if (col1 == col2)
  {          // Cannot exchange the same column.
    return;
  }
  unsigned int c1 = col1 & 31;
  unsigned int c2 = col2 & 31;
  unsigned int mask1 = (int)(0x80000000U >> c1);
  unsigned int mask2 = (int)(0x80000000U >> c2);
  unsigned int notMask1 = ~mask1;
  unsigned int notMask2 = ~mask2;
  if (col1 >= 32)
  {
    matr1 = V1;
    matr3 = XmY;
  }
  else
  {
    matr1 = V2;
    matr3 = V;
  }
  if (col2 >= 32)
  {
    matr2 = V1;
    matr4 = XmY;
  }
  else
  {
    matr2 = V2;
    matr4 = V;
  }
  for (int row = common.siqs.matrixBLength - 1; row >= 0; row--)
  {
    // Exchange columns col1 and col2 of V1:V2
    unsigned int m1 = *matr1;
    unsigned int m2 = *matr2;
    *matr1 = (int)((m1 & notMask1) | (((m2 << c2) & 0x80000000U) >> c1));
    *matr2 = (int)((m2 & notMask2) | (((m1 << c1) & 0x80000000U) >> c2));
    matr1++;
    matr2++;
    // Exchange columns col1 and col2 of XmY:V
    unsigned int m3 = *matr3;
    unsigned int m4 = *matr4;
    *matr3 = (int)((m3 & notMask1) | (((m4 << c2) & 0x80000000U) >> c1));
    *matr4 = (int)((m4 & notMask2) | (((m3 << c1) & 0x80000000U) >> c2));
    matr3++;
    matr4++;
  }
}

static void coladd(int *XmY, int *V, int *V1, int *V2,
  int col1, int col2)
{
  const int* matr1;
  int* matr2;
  const int* matr3;
  int* matr4;

  if (col1 == col2)
  {          // Nothing to do: go out.
    return;
  }
  unsigned int c1 = col1 & 31;
  unsigned int c2 = col2 & 31;
  if (col1 >= 32)
  {
    matr1 = V1;
    matr3 = XmY;
  }
  else
  {
    matr1 = V2;
    matr3 = V;
  }
  if (col2 >= 32)
  {
    matr2 = V1;
    matr4 = XmY;
  }
  else
  {
    matr2 = V2;
    matr4 = V;
  }
  for (int row = common.siqs.matrixBLength - 1; row >= 0; row--)
  {
    // Add column col1 to column col2 of V1:V2
    unsigned int m1 = *matr1;
    unsigned int m2 = *matr2;
    *matr2 = (int)(m2 ^ (((m1 << c1) & 0x80000000U) >> c2));
    matr1++;
    matr2++;
    // Add column col1 to column col2 of XmY:V
    unsigned int m3 = *matr3;
    unsigned int m4 = *matr4;
    *matr4 = (int)(m4 ^ (((m3 << c1) & 0x80000000U) >> c2));
    matr3++;
    matr4++;
  }
}

static bool BlockLanczos(int seed)
{
  int i;
  int j;
  int k;
  int index;
  int mask;
  int matrixD[32];
  int matrixE[32];
  int matrixF[32];
  int matrixWinv[32];
  int matrixWinv1[32];
  int matrixWinv2[32];
  int matrixVtV0[32];
  int matrixVt1V0[32];
  int matrixVt2V0[32];
  int matrixVtAV[32];
  int matrixVt1AV1[32];
  int matrixCalcParenD[32];
  int vectorIndex[64];
  int matrixTemp[32];
  int matrixCalc1[32]; // Matrix that holds temporary data
  int matrixCalc2[32]; // Matrix that holds temporary data
  const int *matr;
  double dSeed;
  double dMult;
  double dDivisor;
  double dAdd;
  int Temp;
  int Temp1;
  int stepNbr = 0;
  int currentOrder;
  int currentMask;
  int row;
  int col;
  int leftCol;
  int rightCol;
  int minind;
  int min;
  int minanswer;
  const int* rowMatrixB;
  int *ptrMatrixXmY;
  int newDiagonalSSt = -1;

  (void)memset(matrixWinv, 0, sizeof(matrixWinv));
  (void)memset(matrixWinv1, 0, sizeof(matrixWinv1));
  (void)memset(matrixWinv2, 0, sizeof(matrixWinv2));
  (void)memset(matrixVtV0, 0, sizeof(matrixVtV0));
  (void)memset(matrixVt1V0, 0, sizeof(matrixVt1V0));
  (void)memset(matrixVt2V0, 0, sizeof(matrixVt2V0));
  (void)memset(matrixVt1AV1, 0, sizeof(matrixVt1AV1));

  /* Initialize matrix X-Y and matrix V_0 with pseudorandom data */
  dSeed = (double)seed;
  dMult = (double)62089911;
  dAdd = (double)54325442;
  dDivisor = (double)0x7FFFFFFF;
  ptrMatrixXmY = &common.siqs.matrixXmY[common.siqs.matrixBLength - 1];
  for (int *ptrMatrixV = &common.siqs.matrixV[common.siqs.matrixBLength - 1];
    ptrMatrixV >= common.siqs.matrixV; ptrMatrixV--)
  {
    double dSeed2 = (dSeed * dMult) + dAdd;
    dSeed2 -= floor(dSeed2 / dDivisor) * dDivisor;
    *ptrMatrixXmY = (int)dSeed + (int)dSeed2;
    ptrMatrixXmY--;
    dSeed = (dSeed2 * dMult) + dAdd;
    dSeed -= floor(dSeed / dDivisor) * dDivisor;
    dSeed2 = (dSeed * dMult) + dAdd;
    dSeed2 -= floor(dSeed2 / dDivisor) * dDivisor;
    *ptrMatrixV = (int)dSeed + (int)dSeed2;
    dSeed = (dSeed2 * dMult) + dAdd;
    dSeed -= floor(dSeed / dDivisor) * dDivisor;
  }
  // Compute matrix Vt(0) * V(0)
  MatrTranspMult(common.siqs.matrixBLength, common.siqs.matrixV, common.siqs.matrixV, matrixVtV0);
#ifdef __EMSCRIPTEN__
  showMatrixSize((char *)SIQSInfoText, matrixRows, matrixCols);
#endif
#if DEBUG_SIQS == 2
  char *ptrOutput = output;
#ifdef __EMSCRIPTEN__
  *ptrOutput = '9';
  ptrOutput++;
#endif
  copyStr(&ptrOutput, "MatrixBLength = ");
  int2dec(&ptrOutput, common.siqs.matrixBLength);
  copyStr(&ptrOutput, ", matrix = ");
  int2dec(&ptrOutput, matrixRows);
  copyStr(&ptrOutput, " * ");
  int2dec(&ptrOutput, matrixCols);
  *ptrOutput = 0;
#ifdef __EMSCRIPTEN__
  databack(output);
#else
  printf("%s\n", output);
#endif
#endif
  for (;;)
  {
    int indexC;
    int oldDiagonalSSt;
    if ((stepNbr * 3200 / matrixRows) > 105)
    {
      return false;
    }
#ifdef __EMSCRIPTEN__
    int elapsedTime = (int)(tenths() - originalTenthSecond);
    if ((elapsedTime / 10) != (oldTimeElapsed / 10))
    {
      char SIQSInfo[200];
      char *ptrText = SIQSInfo;
      oldTimeElapsed = elapsedTime;
      copyStr(&ptrText, "4<p>");
      GetDHMS(&ptrText, elapsedTime/10);
      copyStr(&ptrText, "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;");
      copyStr(&ptrText, lang? "Progreso del álgebra lineal: ": "Linear algebra progress: ");
      int2dec(&ptrText, stepNbr * 3200 / matrixRows);
      copyStr(&ptrText, "%</p>");
      databack(SIQSInfo);
    }  
#endif
#if 0
    //if (getTerminateThread())
    //{
    //  throw new ArithmeticException();
    // }
#endif
    oldDiagonalSSt = newDiagonalSSt;
    stepNbr++;
    // Compute matrix A * V(i)
    MultiplyAByMatrix(common.siqs.matrixV, common.siqs.matrixCalc3, common.siqs.matrixAV);
    // Compute matrix Vt(i) * A * V(i)
    MatrTranspMult(common.siqs.matrixBLength, common.siqs.matrixV, common.siqs.matrixAV, matrixVtAV);
    /* If Vt(i) * A * V(i) = 0, end of loop */
    for (i = sizeof(matrixVtAV)/sizeof(matrixVtAV[0]) - 1; i >= 0; i--)
    {
      if (matrixVtAV[i] != 0)
      {
        break;
      }
    }
    if (i < 0)
    {
      break;
    } /* End X-Y calculation loop */

      /* Selection of S(i) and W(i) */

    (void)memcpy(matrixTemp, matrixWinv2, sizeof(matrixTemp));
    (void)memcpy(matrixWinv2, matrixWinv1, sizeof(matrixTemp));
    (void)memcpy(matrixWinv1, matrixWinv, sizeof(matrixTemp));
    (void)memcpy(matrixWinv, matrixTemp, sizeof(matrixTemp));

    mask = 1;
    for (j = 31; j >= 0; j--)
    {
      matrixD[j] = matrixVtAV[j]; /*  D <- VtAV    */
      matrixWinv[j] = mask;       /*  Winv <- I    */
      mask *= 2;
    }

    index = 31;
    mask = 1;
    for (indexC = 31; indexC >= 0; indexC--)
    {
      if ((oldDiagonalSSt & mask) != 0)
      {
        matrixE[index] = indexC;
        matrixF[index] = mask;
        index--;
      }
      mask *= 2;
    }
    mask = 1;
    for (indexC = 31; indexC >= 0; indexC--)
    {
      if ((oldDiagonalSSt & mask) == 0)
      {
        matrixE[index] = indexC;
        matrixF[index] = mask;
        index--;
      }
      mask *= 2;
    }
    newDiagonalSSt = 0;
    for (j = 0; j < 32; j++)
    {
      currentOrder = matrixE[j];
      currentMask = matrixF[j];
      for (k = j; k < 32; k++)
      {
        if ((matrixD[matrixE[k]] & currentMask) != 0)
        {
          break;
        }
      }
      if (k < 32)
      {
        i = matrixE[k];
        Temp = matrixWinv[i];
        matrixWinv[i] = matrixWinv[currentOrder];
        matrixWinv[currentOrder] = Temp;
        Temp1 = matrixD[i];
        matrixD[i] = matrixD[currentOrder];
        matrixD[currentOrder] = Temp1;
        newDiagonalSSt |= currentMask;
        for (k = 31; k >= 0; k--)
        {
          if ((k != currentOrder) && ((matrixD[k] & currentMask) != 0))
          {
            matrixWinv[k] ^= Temp;
            matrixD[k] ^= Temp1;
          }
        } /* end for k */
      }
      else
      {
        for (k = j; k < 31; k++)
        {
          if ((matrixWinv[matrixE[k]] & currentMask) != 0)
          {
            break;
          }
        }
        i = matrixE[k];
        Temp = matrixWinv[i];
        matrixWinv[i] = matrixWinv[currentOrder];
        matrixWinv[currentOrder] = Temp;
        Temp1 = matrixD[i];
        matrixD[i] = matrixD[currentOrder];
        matrixD[currentOrder] = Temp1;
        for (k = 31; k >= 0; k--)
        {
          if ((matrixWinv[k] & currentMask) != 0)
          {
            matrixWinv[k] ^= Temp;
            matrixD[k] ^= Temp1;
          }
        } /* end for k */
      } /* end if */
    } /* end for j */
      /* Compute D(i), E(i) and F(i) */

    if (stepNbr >= 3)
    {
      // F = -Winv(i-2) * (I - Vt(i-1)*A*V(i-1)*Winv(i-1)) * ParenD * S*St
      MatrixMultiplication(matrixVt1AV1, matrixWinv1, matrixCalc2);
      mask = 1; /* Add identity matrix */
      for (index = 31; index >= 0; index--)
      {
        matrixCalc2[index] ^= mask;
        mask *= 2;
      }
      MatrixMultiplication(matrixWinv2, matrixCalc2, matrixCalc1);
      MatrixMultiplication(matrixCalc1, matrixCalcParenD, matrixF);
      MatrMultBySSt(32, matrixF, newDiagonalSSt, matrixF);
    }
    // E = -Winv(i-1) * Vt(i)*A*V(i) * S*St
    if (stepNbr >= 2)
    {
      MatrixMultiplication(matrixWinv1, matrixVtAV, matrixE);
      MatrMultBySSt(32, matrixE, newDiagonalSSt, matrixE);
    }
    // ParenD = Vt(i)*A*A*V(i) * S*St + Vt(i)*A*V(i)
    // D = I - Winv(i) * ParenD
    MatrTranspMult(common.siqs.matrixBLength, common.siqs.matrixAV, common.siqs.matrixAV, matrixCalc1); // Vt(i)*A*A*V(i)
    MatrMultBySSt(32, matrixCalc1, newDiagonalSSt, matrixCalc1);
    MatrixAddition(matrixCalc1, matrixVtAV, matrixCalcParenD);
    MatrixMultiplication(matrixWinv, matrixCalcParenD, matrixD);
    mask = 1; /* Add identity matrix */
    for (index = 31; index >= 0; index--)
    {
      matrixD[index] ^= mask;
      mask*=2;
    }

    /* Update value of X - Y */
    MatrixMultiplication(matrixWinv, matrixVtV0, matrixCalc1);
    MatrixMultAdd(common.siqs.matrixV, matrixCalc1, common.siqs.matrixXmY);

    /* Compute value of new matrix V(i) */
    // V(i+1) = A * V(i) * S * St + V(i) * D + V(i-1) * E + V(i-2) * F
    MatrMultBySSt(common.siqs.matrixBLength, common.siqs.matrixAV, newDiagonalSSt, common.siqs.matrixCalc3);
    MatrixMultAdd(common.siqs.matrixV, matrixD, common.siqs.matrixCalc3);
#if DEBUG_SIQS == 2
    ptrOutput = output;
    if (stepNbr < 200)
    {
      int sum;
      int ctr;
#ifdef __EMSCRIPTEN__
      * ptrOutput = '9';
      ptrOutput++;
#endif
      copyStr(&ptrOutput, "Step #");
      int2dec(&ptrOutput, stepNbr);
      copyStr(&ptrOutput, ": matrixWinv1[0] = ");
      int2dec(&ptrOutput, matrixWinv1[0]);
      sum = 0;
      for (ctr = 0; ctr < common.siqs.matrixBLength; ctr++)
      {
        sum += common.siqs.matrixV[ctr];
      }
      copyStr(&ptrOutput, ", sum matrixV = ");
      int2dec(&ptrOutput, sum);
      sum = 0;
      for (ctr = 0; ctr < 32; ctr++)
      {
        sum += matrixD[ctr];
      }
      copyStr(&ptrOutput, ", sum matrixD = ");
      int2dec(&ptrOutput, sum);
      copyStr(&ptrOutput, ", newDiagonalSSt = ");
      int2dec(&ptrOutput, newDiagonalSSt);
      *ptrOutput = 0;
#ifdef __EMSCRIPTEN__
      databack(output);
#else
      printf("%s\n", output);
#endif
    }
#endif
    if (stepNbr >= 2)
    {
      MatrixMultAdd(common.siqs.matrixV1, matrixE, common.siqs.matrixCalc3);
      if (stepNbr >= 3)
      {
        MatrixMultAdd(common.siqs.matrixV2, matrixF, common.siqs.matrixCalc3);
      }
    }
    /* Compute value of new matrix Vt(i)V0 */
    // Vt(i+1)V(0) = Dt * Vt(i)V(0) + Et * Vt(i-1)V(0) + Ft * Vt(i-2)V(0)
    MatrTranspMult(32, matrixD, matrixVtV0, matrixCalc2);
    if (stepNbr >= 2)
    {
      MatrTranspMult(32, matrixE, matrixVt1V0, matrixCalc1);
      MatrixAddition(matrixCalc1, matrixCalc2, matrixCalc2);
      if (stepNbr >= 3)
      {
        MatrTranspMult(32, matrixF, matrixVt2V0, matrixCalc1);
        MatrixAddition(matrixCalc1, matrixCalc2, matrixCalc2);
      }
    }
    (void)memcpy(common.siqs.matrixTemp2, common.siqs.matrixV2, sizeof(common.siqs.matrixTemp2));
    (void)memcpy(common.siqs.matrixV2, common.siqs.matrixV1, sizeof(common.siqs.matrixV2));
    (void)memcpy(common.siqs.matrixV1, common.siqs.matrixV, sizeof(common.siqs.matrixV1));
    (void)memcpy(common.siqs.matrixV, common.siqs.matrixCalc3, sizeof(common.siqs.matrixV));
    (void)memcpy(common.siqs.matrixCalc3, common.siqs.matrixTemp2, sizeof(common.siqs.matrixCalc3));
    (void)memcpy(matrixTemp, matrixVt2V0, sizeof(matrixTemp));
    (void)memcpy(matrixVt2V0, matrixVt1V0, sizeof(matrixVt2V0));
    (void)memcpy(matrixVt1V0, matrixVtV0, sizeof(matrixVt1V0));
    (void)memcpy(matrixVtV0, matrixCalc2, sizeof(matrixVtV0));
    (void)memcpy(matrixCalc2, matrixTemp, sizeof(matrixCalc2));
    (void)memcpy(matrixTemp, matrixVt1AV1, sizeof(matrixTemp));
    (void)memcpy(matrixVt1AV1, matrixVtAV, sizeof(matrixVt1AV1));
    (void)memcpy(matrixVtAV, matrixTemp, sizeof(matrixVtAV));
  } /* end while */

    /* Find matrix V1:V2 = B * (X-Y:V) */
  for (row = common.siqs.matrixBLength - 1; row >= 0; row--)
  {
    common.siqs.matrixV1[row] = 0;
    common.siqs.matrixV2[row] = 0;
  }
  for (row = common.siqs.matrixBLength - 1; row >= 0; row--)
  {
    int rowMatrixXmY;
    int rowMatrixV;

    rowMatrixB = common.siqs.matrixB[row];
    rowMatrixXmY = common.siqs.matrixXmY[row];
    rowMatrixV = common.siqs.matrixV[row];
    // The vector rowMatrixB includes the indexes of the columns set to '1'.
    for (index = rowMatrixB[LENGTH_OFFSET]-1; index >= 1; index--)
    {
      col = rowMatrixB[index];
      common.siqs.matrixV1[col] ^= rowMatrixXmY;
      common.siqs.matrixV2[col] ^= rowMatrixV;
    }
  }
  rightCol = 64;
  leftCol = 0;
  while (leftCol < rightCol)
  {
    for (col = leftCol; col < rightCol; col++)
    {       // For each column find the first row which has a '1'.
            // Columns outside this range must have '0' in all rows.
      matr = ((col >= 32) ? common.siqs.matrixV1 : common.siqs.matrixV2);
      mask = 0x80000000U >> (col & 31);
      vectorIndex[col] = -1;    // indicate all rows in zero in advance.
      for (row = 0; row < common.siqs.matrixBLength; row++)
      {
        if ((matr[row] & mask) != 0)
        {               // First row for this mask is found. Store it.
          vectorIndex[col] = row;
          break;
        }
      }
    }
    for (col = leftCol; col < rightCol; col++)
    {
      if (vectorIndex[col] < 0)
      {  // If all zeros in col 'col', exchange it with first column with
         // data different from zero (leftCol).
        colexchange(common.siqs.matrixXmY, common.siqs.matrixV, common.siqs.matrixV1, common.siqs.matrixV2, leftCol, col);
        vectorIndex[col] = vectorIndex[leftCol];
        vectorIndex[leftCol] = -1;  // This column now has zeros.
        leftCol++;                  // Update leftCol to exclude that column.
      }
    }
    if (leftCol == rightCol)
    {
      break;
    }
    // At this moment all columns from leftCol to rightCol are non-zero.
    // Get the first row that includes a '1'.
    min = vectorIndex[leftCol];
    minind = leftCol;
    for (col = leftCol + 1; col < rightCol; col++)
    {
      if (vectorIndex[col] < min)
      {
        min = vectorIndex[col];
        minind = col;
      }
    }
    minanswer = 0;
    for (col = leftCol; col < rightCol; col++)
    {
      if (vectorIndex[col] == min)
      {
        minanswer++;
      }
    }
    if (minanswer > 1)
    {            // Two columns with the same first row to '1'.
      for (col = minind + 1; col < rightCol; col++)
      {
        if (vectorIndex[col] == min)
        {        // Add first column which has '1' in the same row to
                 // the other columns so they have '0' in this row after
                 // this operation.
          coladd(common.siqs.matrixXmY, common.siqs.matrixV, common.siqs.matrixV1, common.siqs.matrixV2, minind, col);
        }
      }
    }
    else
    {
      rightCol--;
      colexchange(common.siqs.matrixXmY, common.siqs.matrixV, common.siqs.matrixV1, common.siqs.matrixV2, minind, rightCol);
    }
  }
  leftCol = 0; /* find linear independent solutions */
  while (leftCol < rightCol)
  {
    for (col = leftCol; col < rightCol; col++)
    {         // For each column find the first row which has a '1'.
      matr = ((col >= 32) ? common.siqs.matrixXmY : common.siqs.matrixV);
      mask = 0x80000000U >> (col & 31);
      vectorIndex[col] = -1;    // indicate all rows in zero in advance.
      for (row = 0; row < common.siqs.matrixBLength; row++)
      {
        if ((matr[row] & mask) != 0)
        {         // First row for this mask is found. Store it.
          vectorIndex[col] = row;
          break;
        }
      }
    }
    col = leftCol;
    while (col < rightCol)
    {  // If all zeros in col 'col', exchange it with last column with
       // data different from zero (rightCol).
      if (vectorIndex[col] < 0)
      {
        rightCol--;                 // Update rightCol to exclude that column.
        colexchange(common.siqs.matrixXmY, common.siqs.matrixV, common.siqs.matrixV1, common.siqs.matrixV2, rightCol, col);
        vectorIndex[col] = vectorIndex[rightCol];
        vectorIndex[rightCol] = -1; // This column now has zeros.
      }
      col++;
    }
    if (leftCol == rightCol)
    {
      break;
    }
    // At this moment all columns from leftCol to rightCol are non-zero.
    // Get the first row that includes a '1'.
    min = vectorIndex[leftCol];
    minind = leftCol;
    for (col = leftCol + 1; col < rightCol; col++)
    {
      if (vectorIndex[col] < min)
      {
        min = vectorIndex[col];
        minind = col;
      }
    }
    minanswer = 0;
    for (col = leftCol; col < rightCol; col++)
    {
      if (vectorIndex[col] == min)
      {
        minanswer++;
      }
    }
    if (minanswer > 1)
    {            // At least two columns with the same first row to '1'.
      for (col = minind + 1; col < rightCol; col++)
      {
        if (vectorIndex[col] == min)
        {        // Add first column which has '1' in the same row to
                 // the other columns so they have '0' in this row after
                 // this operation.
          coladd(common.siqs.matrixXmY, common.siqs.matrixV, common.siqs.matrixV1, common.siqs.matrixV2, minind, col);
        }
      }
    }
    else
    {
      colexchange(common.siqs.matrixXmY, common.siqs.matrixV, common.siqs.matrixV1, common.siqs.matrixV2, minind, leftCol);
      leftCol++;
    }
  }
  return true;
}

static int EraseSingletons(int nbrFactorBasePrimes)
{
  int row;
  int column;
  int delta;
  int* rowMatrixB;
  int matrixBlength = common.siqs.matrixBLength;

  {
    int nbrBytes = matrixBlength * (int)sizeof(int);
    (void)memset(common.siqs.newColumns, 0, nbrBytes);
  }
  // Find singletons in matrixB storing in array vectExpParity the number
  // of primes in each column.
  do
  {   // The singleton removal phase must run until there are no more
      // singletons to erase.
    {
      int nbrBytes = common.siqs.matrixBLength * (int)sizeof(limb);
      (void)memset(common.siqs.vectExpParity, 0, nbrBytes);
    }
    for (row = matrixBlength - 1; row >= 0; row--)
    {                  // Traverse all rows of the matrix.
      rowMatrixB = common.siqs.matrixB[row];
      for (column = rowMatrixB[LENGTH_OFFSET] - 1; column >= 1; column--)
      {                // A prime appeared in that column.
        common.siqs.vectExpParity[rowMatrixB[column]]++;
      }
    }
    row = 0;
    for (column = 0; column < nbrFactorBasePrimes; column++)
    {
      if (common.siqs.vectExpParity[column] > 1)
      {                // Useful column found with at least 2 primes.
#if DEBUG_SIQS == 2
        common.siqs.primeSieveData[row] = common.siqs.primeSieveData[column];
#endif
        common.siqs.newColumns[column] = row;
        common.siqs.primeTrialDivisionData[row].value =
          common.siqs.primeTrialDivisionData[column].value;
        row++;
      }
    }
    nbrFactorBasePrimes = row;
    delta = 0;
    // Erase singletons from matrixB. The rows to be erased are those where the
    // the corresponding element of the array vectExpParity equals 1.
    for (row = 0; row < matrixBlength; row++)
    {                  // Traverse all rows of the matrix.
      rowMatrixB = common.siqs.matrixB[row];
      for (column = rowMatrixB[LENGTH_OFFSET] - 1; column >= 1; column--)
      {                // Traverse all columns.
        if (common.siqs.vectExpParity[rowMatrixB[column]] == 1)
        {              // Singleton found: erase this row.
          delta++;
          break;
        }
      }
      if ((column == 0) && (delta != 0))
      {                // Singleton not found: move row upwards.
        (void)memcpy(common.siqs.matrixB[row - delta],
          common.siqs.matrixB[row], sizeof(common.siqs.matrixB[0]));
        (void)memcpy(common.siqs.vectLeftHandSide[row - delta],
          common.siqs.vectLeftHandSide[row], sizeof(common.siqs.vectLeftHandSide[0]));
      }
    }
    matrixBlength -= delta;      // Update number of rows of the matrix.
    for (row = 0; row < matrixBlength; row++)
    {                  // Traverse all rows of the matrix.
      rowMatrixB = common.siqs.matrixB[row];
      for (column = rowMatrixB[LENGTH_OFFSET]; column >= 1; column--)
      {                // Change all column indexes in this row.
        rowMatrixB[column] = common.siqs.newColumns[rowMatrixB[column]];
      }
    }
  } while (delta > 0);           // End loop if number of rows did not
                                 // change.
  common.siqs.primeTrialDivisionData[0].exp2 = nbrFactorBasePrimes;
  return matrixBlength;
}

#if DEBUG_SIQS == 2
static void showRelations(void)
{
#ifdef __EMSCRIPTEN__
  databack("9******* START LINEAR ALGEBRA *******\n");
#else
  printf("******* START LINEAR ALGEBRA *******\n");
#endif
  for (int j = 0; j < common.siqs.matrixBLength; j++)
  {
    char* ptrOutput = output;
#ifdef __EMSCRIPTEN__
    * ptrOutput = '9';
    ptrOutput++;
#endif
    copyStr(&ptrOutput, "Mod(");
    for (int i = 1; i < common.siqs.matrixB[j][LENGTH_OFFSET]; i++)
    {
      if (i != 1)
      {
        *ptrOutput++ = '*';
      }
      if (common.siqs.matrixB[j][i] == 0)
      {
        copyStr(&ptrOutput, "(-1)");
      }
      else
      {
        int2dec(&ptrOutput, common.siqs.primeSieveData[common.siqs.matrixB[j][i]].value);
      }
    }
    copyStr(&ptrOutput, " - ");
    static BigInteger k1;
    (void)memcpy(k1.limbs, common.siqs.vectLeftHandSide[j],
      NumberLength * sizeof(limb));
    k1.nbrLimbs = NumberLength;
    k1.sign = SIGN_POSITIVE;
    BigInteger2Dec(&ptrOutput, &k1, 0);
    ShowSquareModP(ptrOutput);
  }
}
#endif

/************************/
/* Linear algebra phase */
/************************/
bool LinearAlgebraPhase(limb* biT, limb* biR, limb* biU, int nbrLength)
{
  int seed = 123456789;
  int mask;
  const int* rowMatrixB;
  int primeIndex;
  // Get new number of rows after erasing singletons.
  int matrixBlength = EraseSingletons(common.siqs.nbrFactorBasePrimes);
  common.siqs.matrixBLength = matrixBlength;
  matrixRows = matrixBlength;
  matrixCols = common.siqs.primeTrialDivisionData[0].exp2;
  common.siqs.primeTrialDivisionData[0].exp2 = 0;         // Restore correct value.
#if DEBUG_SIQS == 2
  showRelations();
#endif
  for (int ctr = 0; ctr < 10; ctr++)
  {
    while (BlockLanczos(seed + 1) == false)
    {   // Block Lanczos does not work with this seed. Try another one.
      seed++;
    }
    // The rows of matrixV indicate which rows must be multiplied so no
    // primes are multiplied an odd number of times.
    mask = 1;
    for (int col = 31; col >= 0; col--)
    {
      int NumberLengthBak;
      int index;

      IntToBigNbr(1, biT, nbrLength + 1);
      IntToBigNbr(1, biR, nbrLength + 1);
      {
        int nbrBytes = matrixBlength * (int)sizeof(common.siqs.vectExpParity[0]);
        (void)memset(common.siqs.vectExpParity, 0, nbrBytes);
      }
      NumberLengthBak = nbrLength;
      if (common.siqs.Modulus[nbrLength - 1].x == 0)
      {
        nbrLength--;
      }
      for (int row = matrixBlength - 1; row >= 0; row--)
      {
        if ((common.siqs.matrixV[row] & mask) != 0)
        {
          MultBigNbrModN(common.siqs.vectLeftHandSide[row], biR, biU, common.siqs.Modulus,
            nbrLength);
          {
            int nbrBytes = (nbrLength + 1) * (int)sizeof(biR[0]);
            (void)memcpy(biR, biU, nbrBytes);
          }
          rowMatrixB = common.siqs.matrixB[row];
          for (int j = rowMatrixB[LENGTH_OFFSET] - 1; j >= 1; j--)
          {
            primeIndex = rowMatrixB[j];
            common.siqs.vectExpParity[primeIndex] ^= 1;
            if (common.siqs.vectExpParity[primeIndex] == 0)
            {
              if (primeIndex == 0)
              {
                SubtractBigNbr(common.siqs.Modulus, biT, biT, nbrLength); // Multiply biT by -1.
              }
              else
              {
                MultBigNbrByIntModN(biT,
                  common.siqs.primeTrialDivisionData[primeIndex].value, biT,
                  common.siqs.Modulus, nbrLength);
              }
            }
          }
        }
      }
      nbrLength = NumberLengthBak;
      SubtractBigNbrModN(biR, biT, biR, common.siqs.Modulus, nbrLength);
      GcdBigNbr(biR, common.siqs.TestNbr2, biT, nbrLength);
      for (index = 1; index < nbrLength; index++)
      {
        if (biT[index].x != 0)
        {
          break;
        }
      }
      if ((index < nbrLength) || (biT[0].x > 1))
      {   // GCD is not zero or 1.
        for (index = 0; index < nbrLength; index++)
        {
          if (biT[index].x != common.siqs.TestNbr2[index].x)
          {
            break;
          }
        }
        if (index < nbrLength)
        { /* GCD is not 1 */
          return true;
        }
      }
      mask *= 2;
    }
  }
  return false;
}



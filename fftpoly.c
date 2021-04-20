//
// This file is part of Alpertron Calculators.
//
// Copyright 2020-2021 Dario Alejandro Alpern
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

#include "bignbr.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "polynomial.h"
#include "fft.h"

#define FFT_LIMB_SIZE   22
#define MAX_FFT_LEN     2048  //  Power of 2 greater than 2*MAX_DEGREE
#define POWERS_2        13
// In the next array, all numbers are represented by two elements,
// first the least significant limb, then the most significant limb.

static struct sCosSin cossin[4 << (POWERS_2 - 2)];
static double Cosine[5 * QUARTER_CIRCLE + 1];
static complex firstFactor[MAX_FFT_LEN];
static complex secondFactor[MAX_FFT_LEN];
static complex transf[MAX_FFT_LEN];
static complex product[MAX_FFT_LEN];
static complex finalProduct[MAX_FFT_LEN];
static complex tempFFT[MAX_FFT_LEN];
static complex polyInvTransf[MAX_FFT_LEN];
extern int polyInv[COMPRESSED_POLY_MAX_LENGTH];

// Use formulas sin(A+B) = sin A cos B + cos A sin B
// and cos(A+B) = cos A cos B - sin A sin B
static void initCosinesArray(void)
{
  struct sCosSin* ptrCosSin;
  const struct sCosSin *ptrOldCosSin;
  const struct sCosSin *ptrCosSinDelta;
  double invLimb = 1 / (double)LIMB_RANGE;
  double invSqLimb = invLimb * invLimb;
  int index;
  cossin[0].Cos[0] = MAX_VALUE_LIMB;                       // cos(0) = 1
  cossin[0].Cos[1] = MAX_VALUE_LIMB;
  cossin[0].Sin[0] = 0;                                    // sin(0) = 0
  cossin[0].Sin[1] = 0;
  ptrCosSin = &cossin[1];
  for (index = 1; ; index++)
  {
    // Get order of least significant non-zero bit.
    int bitNbr;
    int mask = 1;
    for (bitNbr = 0; ; bitNbr++)
    {
      if (index & mask)
      {
        break;
      }
      mask *= 2;
    }
    if (bitNbr == POWERS_2 - 2)
    {
      break;
    }
    ptrCosSinDelta = &cossinPowerOneHalf[(POWERS_2 - 3 - bitNbr)];  // Pointer to cos/sin B.
    if (index == mask)
    {
      ptrCosSin->Cos[0] = ptrCosSinDelta->Cos[0];
      ptrCosSin->Cos[1] = ptrCosSinDelta->Cos[1];
      ptrCosSin->Sin[0] = ptrCosSinDelta->Sin[0];
      ptrCosSin->Sin[1] = ptrCosSinDelta->Sin[1];
    }
    else
    {
      int firstProd[6];
      int secondProd[6];
      // Compute cos(A+B) = cos A cos B - sin A sin B.
      ptrOldCosSin = ptrCosSin - mask;   // Pointer to cos/sin A.
      MultBigNbrComplete(ptrOldCosSin->Cos, ptrCosSinDelta->Cos, firstProd, 2);
      MultBigNbrComplete(ptrOldCosSin->Sin, ptrCosSinDelta->Sin, secondProd, 2);
      SubtractBigNbr(firstProd, secondProd, firstProd, 4);
      ptrCosSin->Cos[0] = *(firstProd + 2);
      ptrCosSin->Cos[1] = *(firstProd + 3);
      // Compute sin(A+B) = sin A cos B + cos A sin B.
      MultBigNbrComplete(ptrOldCosSin->Sin, ptrCosSinDelta->Cos, firstProd, 2);
      MultBigNbrComplete(ptrOldCosSin->Cos, ptrCosSinDelta->Sin, secondProd, 2);
      AddBigNbr(firstProd, secondProd, firstProd, 4);
      ptrCosSin->Sin[0] = *(firstProd + 2);
      ptrCosSin->Sin[1] = *(firstProd + 3);
    }
    ptrCosSin++;
  }
  // Convert from integers to doubles and send the results to the final array.
  ptrCosSin = cossin;
  for (index = 0; index < QUARTER_CIRCLE; index++)
  {
    double cosine = (double)ptrCosSin->Cos[0] * invSqLimb + (double)ptrCosSin->Cos[1] * invLimb;
    Cosine[index] = cosine;
    Cosine[HALF_CIRCLE - index] = -cosine;
    Cosine[HALF_CIRCLE + index] = -cosine;
    Cosine[2 * HALF_CIRCLE - index] = cosine;
    Cosine[2 * HALF_CIRCLE + index] = cosine;
    ptrCosSin++;
  }
  Cosine[QUARTER_CIRCLE] = 0;
  Cosine[3 * QUARTER_CIRCLE] = 0;
}

/*
  Algorithm 9.5.6 of Crandall and Pomerance book Prime Numbers:
  X, Y: Pointers to complex numbers.

  J = 1; X = x; Y = y;
  for (d >= i > 0)
  {
    m = 0;
    while (m < D/2)
    {
      a = exp(-2*pi*i*m/D);
      for (J >= j > 0)
      {
        Y[0] = X[0] + X[D/2];
        Y[J] = a(X[0] - X[D/2]);
        X = X + 1;
        Y = Y + 1;
      }
      Y = Y + J;
      m = m + J;
    }
    J = 2 * J;
    X = X - D/2;
    Y = Y - D;
    (X,Y) = (Y,X);
  }
  if (d even) return complex data at X.
  return complex data at Y.
*/

// length is power of 2.
static void complexPolyFFT(complex* x, complex* y, int length)
{
  int j;
  int J;
  int halfLength = length / 2;
  int step = (1 << POWERS_2) / length;
  int exponentOdd = 0;
  complex* ptrX = x;
  complex* ptrY = y;
  complex* ptrZ;
  complex* ptrTemp;
  int angle;
  if (Cosine[0] == 0)
  {    // Cosines array not initialized yet. Initialize array.
    initCosinesArray();
  }
  ptrZ = ptrX + halfLength;
  for (angle = 0; angle < HALF_CIRCLE; angle += step)
  {
    double rootReal = Cosine[angle];
    double rootImag = Cosine[angle + QUARTER_CIRCLE];
    double tempReal = ptrX->real;
    double tempImag = ptrX->imaginary;
    double Zreal = ptrZ->real;
    double Zimag = ptrZ->imaginary;
    ptrY->real = tempReal + Zreal;
    ptrY->imaginary = tempImag + Zimag;
    tempReal -= Zreal;
    tempImag -= Zimag;
    ptrY++;
    ptrY->real = rootReal * tempReal - rootImag * tempImag;
    ptrY->imaginary = rootReal * tempImag + rootImag * tempReal;
    ptrX++;
    ptrY++;
    ptrZ++;
  }
  for (J = 2; J < length; J *= 2)
  {
    step *= 2;
    ptrTemp = ptrX - halfLength;
    ptrX = ptrY - length;
    ptrY = ptrTemp;
    ptrZ = ptrX + halfLength;
    exponentOdd = 1 - exponentOdd;
    for (angle = 0; angle < HALF_CIRCLE; angle += step)
    {
      double rootReal = Cosine[angle];
      double rootImag = Cosine[angle + QUARTER_CIRCLE];
      complex* ptrW = ptrY + J;
      for (j = J; j > 0; j--)
      {
        double tempReal = ptrX->real;
        double tempImag = ptrX->imaginary;
        double Zreal = ptrZ->real;
        double Zimag = ptrZ->imaginary;
        ptrY->real = tempReal + Zreal;
        ptrY->imaginary = tempImag + Zimag;
        tempReal -= Zreal;
        tempImag -= Zimag;
        ptrW->real = rootReal * tempReal - rootImag * tempImag;
        ptrW->imaginary = rootReal * tempImag + rootImag * tempReal;
        ptrX++;
        ptrY++;
        ptrZ++;
        ptrW++;
      }
      ptrY += J;
    }
  }
  if (exponentOdd)
  {     // Move data from x to y.
    (void)memcpy(y, x, length * sizeof(complex));
  }
}

// Formulas to use:
// Gr(k) = Xr(k)Ar(k) – Xi(k)Ai(k) + Xr(N/2–k)Br(k) + Xi(N/2–k)Bi(k)
// Gi(k) = Xi(k)Ar(k) + Xr(k)Ai(k) + Xr(N/2–k)Bi(k) – Xi(N/2–k)Br(k)
// for k = 0, 1, ..., N/2–1 and X(N/2) = X(0)
// Ar(k) = 1 – sin( PI k / N)
// Ai(k) = –cos( PI k / N)
// Br(k) = 1 + sin( PI k / N)
// Bi(k) = cos( PI k / N)
static void ConvertHalfToFullSizeFFT(complex* halfSizeFFT, complex* fullSizeFFT, int power2)
{
  int k;
  int step = (1 << (POWERS_2 - 1)) / power2;
  complex* ptrFullSizeFFT = fullSizeFFT;
  complex* ptrHalfSizeFFT = halfSizeFFT;
  complex* ptrHalfSizeFFTRev = halfSizeFFT + power2;
  ptrHalfSizeFFTRev->real = halfSizeFFT->real;
  ptrHalfSizeFFTRev->imaginary = halfSizeFFT->imaginary;
  for (k = 0; k < power2; k++)
  {
    int angle = k * step;
    double diffReal = ptrHalfSizeFFT->real - ptrHalfSizeFFTRev->real;
    double sumImag = ptrHalfSizeFFT->imaginary + ptrHalfSizeFFTRev->imaginary;
    double negativeSine = Cosine[angle + QUARTER_CIRCLE];
    double cosine = Cosine[angle];
    ptrFullSizeFFT->real = ptrHalfSizeFFT->real + ptrHalfSizeFFTRev->real +
      diffReal * negativeSine +
      sumImag * cosine;
    ptrFullSizeFFT->imaginary = ptrHalfSizeFFT->imaginary - ptrHalfSizeFFTRev->imaginary +
      sumImag * negativeSine -
      diffReal * cosine;
    ptrHalfSizeFFT++;
    ptrHalfSizeFFTRev--;
    ptrFullSizeFFT++;
  }
  ptrFullSizeFFT->real = 2 * (halfSizeFFT->real - halfSizeFFT->imaginary);
  ptrFullSizeFFT->imaginary = 0;
}

// Formulas to use:
// Xr(k) = Gr(k)IAr(k) – Gi(k)IAi(k) + Gr(N/2–k)IBr(k) + Gi(N/2–k)IBi(k)
// Xi(k) = Gi(k)IAr(k) + Gr(k)IAi(k) + Gr(N/2–k)IBi(k) – Gi(N/2–k)IBr(k)
// for k = 0, 1, ..., N/2–1
// IAr(k) = 1 – sin( PI k / N)
// IAi(k) = cos( PI k / N)
// IBr(k) = 1 + sin( PI k / N)
// IBi(k) = -cos( PI k / N)
static void ConvertFullToHalfSizeFFT(complex* fullSizeFFT, complex* halfSizeFFT, int power2)
{
  int k;
  int step = (1 << (POWERS_2 - 1)) / power2;
  complex* ptrFullSizeFFT = fullSizeFFT;
  complex* ptrFullSizeFFTRev = fullSizeFFT + power2;
  complex* ptrHalfSizeFFT = halfSizeFFT;
  for (k = 0; k < power2; k++)
  {
    int angle = k * step;
    double diffReal = ptrFullSizeFFT->real - ptrFullSizeFFTRev->real;
    double sumImag = ptrFullSizeFFT->imaginary + ptrFullSizeFFTRev->imaginary;
    double negativeSine = Cosine[angle + QUARTER_CIRCLE];
    double cosine = Cosine[angle];
    ptrHalfSizeFFT->real = ptrFullSizeFFT->real + ptrFullSizeFFTRev->real +
      diffReal * negativeSine -
      sumImag * cosine;
    // Negative sign for imaginary part required for inverse FFT.
    ptrHalfSizeFFT->imaginary = -(ptrFullSizeFFT->imaginary - ptrFullSizeFFTRev->imaginary +
      sumImag * negativeSine +
      diffReal * cosine);
    ptrHalfSizeFFT++;
    ptrFullSizeFFT++;
    ptrFullSizeFFTRev--;
  }
}

static void ConvertFactorToInternal(int* factor, complex* fftFactor, int len, int maxLen)
{
  int ctr = 0;
  int* ptrFactor = factor+1;  // Point to constant coefficient.
  complex* ptrInternalFactor = fftFactor;
  for (ctr = 2; ctr <= len; ctr += 2)
  {
    ptrInternalFactor->real = *ptrFactor;
    ptrFactor += 2;            // Point to next coefficient.
    ptrInternalFactor++->imaginary = *ptrFactor;
    ptrFactor += 2;            // Point to next coefficient.
  }
  if (len & 1)
  {
    ctr += 2;
    ptrInternalFactor->real = *ptrFactor;
    ptrInternalFactor++->imaginary = 0;
  }
  maxLen++;
  for (; ctr <= maxLen; ctr += 2)
  {
    ptrInternalFactor->real = 0;
    ptrInternalFactor++->imaginary = 0;
  }
}

/*
   Algorithm 9.5.12 of Crandall and Pomerance book Prime Numbers:

   Zero-pad x and y until each has length 2D.
   X = DFT(x);
   Y = DFT(y);
   Z = X * Y;      // Using convolution
   z = DFT^(-1)(Z)
   z = round(z)    // Round elementwise
   Delete leading zeros.
*/
// If degree of first polynomial is greater than the degree of the second polynomial,
// subdivide the coefficients of the first polynomial in groups of K, where K is the
// lowest power of 2 greater or equal than the length of the second polynomial.
void fftPolyMult(int *factor1, int* factor2, int* result, int len1, int len2)
{
  complex *ptrFirst;
  complex *ptrProduct;
  double invPower2;
  int power2plus1;
  int* ptrResult;
  int chunkLen;
  int index;
  int nbrLimbs = NumberLength + 1;
  int factor1DegreesProcessed = 0;
  int power2SecondFactor = 0;
  int power2;
  int modulus = TestNbr[0].x;
  complex* ptrFinalProduct;
  if (len1 > len2)
  { // Degree of first polynomial is greater than degree of second polynomial.
    // Set results to polynomial zero.
    int ctr;
    ptrFinalProduct = finalProduct;
    chunkLen = (len1 + len2 + 1) / 2;
    for (ctr = 0; ctr <= chunkLen; ctr++)
    {
      ptrFinalProduct->real = 0;         // Initialize coefficient to zero.
      ptrFinalProduct->imaginary = 0;
      ptrFinalProduct++;
    }
  }
  // Get least power of 2 greater or equal than the degree of second factor.
  for (power2SecondFactor = 1; power2SecondFactor < len2; power2SecondFactor *= 2)
  {
  }
  // Get transform of second factor outside multiplication loop,
  // because it has to be computed only once.
  power2 = power2SecondFactor;
  power2plus1 = power2 + 1;
  if (factor1 != factor2)
  {
    if ((polyInvCached == NBR_CACHED) && (factor2 == polyInv))
    {   // Get transform of inverse of polynomial from cache.
      (void)memcpy(transf, polyInvTransf, power2plus1 * sizeof(transf[0]));
    }
    else
    {   // Second factor is not cached. Compute transform.
      ConvertFactorToInternal(factor2, secondFactor, len2, 2 * power2);
      complexPolyFFT(secondFactor, tempFFT, power2);
      ConvertHalfToFullSizeFFT(tempFFT, transf, power2);  // transf <- DFT(secondFactor)
    }
    if ((polyInvCached == NBR_READY_TO_BE_CACHED) && (factor2 == polyInv))
    {   // Save transform of inverse of polynomial to cache.
      (void)memcpy(polyInvTransf, transf, power2plus1 * sizeof(transf[0]));
      polyInvCached = NBR_CACHED;
    }
  }
  for (factor1DegreesProcessed = 0; factor1DegreesProcessed < len1;
    factor1DegreesProcessed += power2SecondFactor)
  {
    int index;
    complex * ptrSecond;
    int lenFirstFactor = power2SecondFactor;
    if (lenFirstFactor > len1 - factor1DegreesProcessed)
    {
      lenFirstFactor = len1 - factor1DegreesProcessed;
    }
    // Get transform of first polynomial.
    ConvertFactorToInternal(factor1 + factor1DegreesProcessed*nbrLimbs,
      firstFactor, lenFirstFactor, 2*power2);
    complexPolyFFT(firstFactor, tempFFT, power2);
    ConvertHalfToFullSizeFFT(tempFFT, product, power2);   // product <- DFT(firstFactor)

    // If first factor is equal to second factor and this is the first loop,
    // use transform of second factor as the transform of first factor.
    if ((factor1DegreesProcessed == 0) && (factor1 == factor2))
    {
      (void)memcpy(transf, product, power2plus1 * sizeof(product[0]));   // transf <- DFT(secondFactor)
    }
    // Perform convolution.
    // Overwrite transform of first factor with transform of product.
    ptrFirst = product;
    ptrSecond = transf;
    for (index = 0; index <= power2; index++)
    {   // Perform complex multiplication componentwise.
      double real = ptrFirst->real * ptrSecond->real - ptrFirst->imaginary * ptrSecond->imaginary;
      ptrFirst->imaginary = ptrFirst->real * ptrSecond->imaginary + ptrFirst->imaginary * ptrSecond->real;
      ptrFirst->real = real;
      ptrFirst++;
      ptrSecond++;
    }
    ConvertFullToHalfSizeFFT(product, tempFFT, power2);
    // Recover product from the transform using inverse DFT.
    complexPolyFFT(tempFFT, product, power2);
    // Divide by 8*power2 rounding to nearest integer to get
    // coefficients of product.
    chunkLen = (len1 + len2 - factor1DegreesProcessed + 1)/2;
    if (chunkLen > power2)
    {
      chunkLen = power2;
    }
    if (len1 > len2)
    {
      ptrFinalProduct = &finalProduct[factor1DegreesProcessed / 2];
      ptrProduct = product;
      for (index = 0; index < chunkLen; index++)
      {
        ptrFinalProduct->real += ptrProduct->real;
        ptrFinalProduct->imaginary += ptrProduct->imaginary;
        ptrProduct++;
        ptrFinalProduct++;
      }
    }
  }
  if (len1 > len2)
  {
    ptrProduct = finalProduct;
  }
  else
  {
    ptrProduct = product;
  }
  invPower2 = (double)1 / ((double)(power2 * 8));
  ptrResult = result;
  chunkLen = (len1 + len2 + 1) / 2;
  for (index = 0; index < chunkLen; index++)
  {
    int coeff = (int)floor(ptrProduct->real * invPower2 + 0.5);
    *ptrResult++ = 1;
    *ptrResult++ = coeff % modulus;
    // Imaginary part. Use negative value for inverse FFT.
    coeff = (int)floor(-ptrProduct->imaginary * invPower2 + 0.5);
    *ptrResult++ = 1;
    *ptrResult++ = coeff % modulus;
    ptrProduct++;
  }
}

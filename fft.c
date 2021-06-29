//
// This file is part of Alpertron Calculators.
//
// Copyright 2018-2021 Dario Alejandro Alpern
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
#include "fft.h"

#define FFT_LIMB_SIZE   18
#define FFT_LIMB_RANGE  0x00040000     // 2^18
#define MAX_FFT_LEN     (((MAX_LEN * BITS_PER_GROUP) / FFT_LIMB_SIZE) + 10)
#define POWERS_2        17
#define FULL_CIRCLE     0x00020000     // 2^17
// In the next array, all numbers are represented by two elements,
// first the least significant limb, then the most significant limb.
const struct sCosSin cossinPowerOneHalf[] =
{  // cos(pi/2^n), then sin(pi/2^n)
  {{{2121767201}, {1518500249}}, {{2121767201}, {1518500249}}},  // n = 2
  {{{1696238673}, {1984016188}}, {{782852818}, {821806413}}},    // n = 3
  {{{1857642581}, {2106220351}}, {{886244699}, {418953276}}},    // n = 4
  {{{575294268}, {2137142927}}, {{174918392}, {210490206}}},     // n = 5
  {{{1926953927}, {2144896909}}, {{565903997}, {105372028}}},    // n = 6
  {{{161094006}, {2146836866}}, {{2050385888}, {52701886}}},     // n = 7
  {{{925218479}, {2147321946}}, {{1727413283}, {26352927}}},     // n = 8
  {{{487924891}, {2147443222}}, {{2040267204}, {13176711}}},     // n = 9
  {{{1144652709}, {2147473541}}, {{2107197813}, {6588386}}},     // n = 10
  {{{819842189}, {2147481121}}, {{786843438}, {3294197}}},       // n = 11
  {{{741631966}, {2147483016}}, {{360077744}, {1647099}}},       // n = 12
  {{{185395523}, {2147483490}}, {{1383830442}, {823549}}},       // n = 13
  {{{1120089925}, {2147483608}}, {{1781913263}, {411774}}},      // n = 14
  {{{280022433}, {2147483638}}, {{892988659}, {205887}}},        // n = 15
  {{{1143747429}, {2147483645}}, {{1520490157}, {102943}}},      // n = 16
};

static struct sCosSin cossin[FULL_CIRCLE];
static double Cosine[(5 * QUARTER_CIRCLE) + 1];
static struct sComplex firstFactor[MAX_FFT_LEN];
static struct sComplex secondFactor[MAX_FFT_LEN];
static struct sComplex transf[MAX_FFT_LEN];
static struct sComplex product[MAX_FFT_LEN];
static struct sComplex tempFFT[MAX_FFT_LEN];
static struct sComplex MontgomeryMultNTransf[MAX_FFT_LEN];
static struct sComplex TestNbrTransf[MAX_FFT_LEN];
// Use formulas sin(A+B) = sin A cos B + cos A sin B
// and cos(A+B) = cos A cos B - sin A sin B
static void initCosinesArray(void)
{
  struct sCosSin* ptrCosSin;
  const struct sCosSin* ptrOldCosSin;
  const struct sCosSin* ptrCosSinDelta;
  double invLimb = 1.0 / (double)LIMB_RANGE;
  double invSqLimb = invLimb * invLimb;
  int index = 1;
  cossin[0].Cos[0].x = (int)MAX_VALUE_LIMB;                  // cos(0) = 1
  cossin[0].Cos[1].x = (int)MAX_VALUE_LIMB;
  cossin[0].Sin[0].x = 0;                                    // sin(0) = 0
  cossin[0].Sin[1].x = 0;
  ptrCosSin = &cossin[1];
  for (;;)
  {
    // Get order of least significant non-zero bit.
    int bitNbr = 0;
    int mask = 1;
    for (;;)
    {
      if ((index & mask) != 0)
      {
        break;
      }
      mask *= 2;
      bitNbr++;
    }
    if (bitNbr == (POWERS_2 - 2))
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
      limb firstProd[6];
      limb secondProd[6];
      // Compute cos(A+B) = cos A cos B - sin A sin B.
      ptrOldCosSin = ptrCosSin - mask;   // Pointer to cos/sin A.
      MultBigNbrComplete(ptrOldCosSin->Cos, ptrCosSinDelta->Cos, firstProd, 2);
      MultBigNbrComplete(ptrOldCosSin->Sin, ptrCosSinDelta->Sin, secondProd, 2);
      SubtractBigNbr(firstProd, secondProd, firstProd, 4);
      ptrCosSin->Cos[0].x = firstProd[2].x;
      ptrCosSin->Cos[1].x = firstProd[3].x;
      // Compute sin(A+B) = sin A cos B + cos A sin B.
      MultBigNbrComplete(ptrOldCosSin->Sin, ptrCosSinDelta->Cos, firstProd, 2);
      MultBigNbrComplete(ptrOldCosSin->Cos, ptrCosSinDelta->Sin, secondProd, 2);
      AddBigNbr(firstProd, secondProd, firstProd, 4);
      ptrCosSin->Sin[0].x = firstProd[2].x;
      ptrCosSin->Sin[1].x = firstProd[3].x;
    }
    ptrCosSin++;
    index++;
  }
  // Convert from integers to doubles and send the results to the final array.
  ptrCosSin = cossin;
  for (index = 0; index < QUARTER_CIRCLE; index++)
  {
    double cosine = ((double)ptrCosSin->Cos[0].x * invSqLimb) +
      ((double)ptrCosSin->Cos[1].x * invLimb);
    Cosine[index] = cosine;
    Cosine[HALF_CIRCLE - index] = -cosine;
    Cosine[HALF_CIRCLE + index] = -cosine;
    Cosine[(2 * HALF_CIRCLE) - index] = cosine;
    Cosine[(2 * HALF_CIRCLE) + index] = cosine;
    ptrCosSin++;
  }
  Cosine[QUARTER_CIRCLE] = 0;
  Cosine[3*QUARTER_CIRCLE] = 0;
}
#if 0
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
  if (d even)
  {
    return complex data at X.
  }
  return complex data at Y.
#endif

// length is power of 2.
static void complexFFT(struct sComplex *x, struct sComplex *y, int length)
{
  int halfLength = length / 2;
  int step = FULL_CIRCLE / length;
  bool exponentOdd = false;
  struct sComplex *ptrX = x;
  struct sComplex *ptrY = y;
  const struct sComplex *ptrZ;
  struct sComplex *ptrTemp;
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
    ptrY->real = (rootReal * tempReal) - (rootImag * tempImag);
    ptrY->imaginary = (rootReal * tempImag) + (rootImag * tempReal);
    ptrX++;
    ptrY++;
    ptrZ++;
  }
  for (int J = 2; J < length; J *= 2)
  {
    step *= 2;
    ptrTemp = ptrX - halfLength;
    ptrX = ptrY - length;
    ptrY = ptrTemp;
    ptrZ = ptrX + halfLength;
    exponentOdd = !exponentOdd;
    for (angle = 0; angle < HALF_CIRCLE; angle += step)
    {
      double rootReal = Cosine[angle];
      double rootImag = Cosine[angle + QUARTER_CIRCLE];
      struct sComplex *ptrW = ptrY + J;
      for (int j = J; j > 0; j--)
      {
        double tempReal = ptrX->real;
        double tempImag = ptrX->imaginary;
        double Zreal = ptrZ->real;
        double Zimag = ptrZ->imaginary;
        ptrY->real = tempReal + Zreal;
        ptrY->imaginary = tempImag + Zimag;
        tempReal -= Zreal;
        tempImag -= Zimag;
        ptrW->real = (rootReal * tempReal) - (rootImag * tempImag);
        ptrW->imaginary = (rootReal * tempImag) + (rootImag * tempReal);
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
    int lengthBytes = length * (int)sizeof(struct sComplex);
    (void)memcpy(y, x, lengthBytes);
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
static void ConvertHalfToFullSizeFFT(struct sComplex *halfSizeFFT, 
  struct sComplex *fullSizeFFT, int power2)
{
  int step = HALF_CIRCLE / power2;
  struct sComplex *ptrFullSizeFFT = fullSizeFFT;
  const struct sComplex *ptrHalfSizeFFT = halfSizeFFT;
  struct sComplex *ptrHalfSizeFFTRev = halfSizeFFT + power2;
  ptrHalfSizeFFTRev->real = halfSizeFFT->real;
  ptrHalfSizeFFTRev->imaginary = halfSizeFFT->imaginary;
  for (int k = 0; k < power2; k++)
  {
    int angle = k * step;
    double diffReal = ptrHalfSizeFFT->real - ptrHalfSizeFFTRev->real;
    double sumImag = ptrHalfSizeFFT->imaginary + ptrHalfSizeFFTRev->imaginary;
    double negativeSine = Cosine[angle + QUARTER_CIRCLE];
    double cosine = Cosine[angle];
    ptrFullSizeFFT->real = ptrHalfSizeFFT->real + ptrHalfSizeFFTRev->real +
      (diffReal * negativeSine) +
      (sumImag * cosine);
    ptrFullSizeFFT->imaginary = ptrHalfSizeFFT->imaginary - ptrHalfSizeFFTRev->imaginary +
      (sumImag * negativeSine) -
      (diffReal * cosine);
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
static void ConvertFullToHalfSizeFFT(const struct sComplex *fullSizeFFT,
  struct sComplex *halfSizeFFT, int power2)
{
  int step = HALF_CIRCLE / power2;
  const struct sComplex *ptrFullSizeFFT = fullSizeFFT;
  const struct sComplex *ptrFullSizeFFTRev = fullSizeFFT + power2;
  struct sComplex *ptrHalfSizeFFT = halfSizeFFT;
  for (int k = 0; k < power2; k++)
  {
    int angle = k * step;
    double diffReal = ptrFullSizeFFT->real - ptrFullSizeFFTRev->real;
    double sumImag = ptrFullSizeFFT->imaginary + ptrFullSizeFFTRev->imaginary;
    double negativeSine = Cosine[angle + QUARTER_CIRCLE];
    double cosine = Cosine[angle];
    ptrHalfSizeFFT->real = ptrFullSizeFFT->real + ptrFullSizeFFTRev->real +
      (diffReal * negativeSine) -
      (sumImag * cosine);
    // Negative sign for imaginary part required for inverse FFT.
    ptrHalfSizeFFT->imaginary = -(ptrFullSizeFFT->imaginary - ptrFullSizeFFTRev->imaginary +
      (sumImag * negativeSine) +
      (diffReal * cosine));
    ptrHalfSizeFFT++;
    ptrFullSizeFFT++;
    ptrFullSizeFFTRev--;
  }
}

// Read limbs of input numbers with length BITS_PER_GROUP and convert them
// to internal FFT limbs with length FFT_LIMB_SIZE.
// The output is the number of internal FFT limbs generated.
// Even limbs correspond to real component and odd limbs
// correspond to imaginary component.
// This code requires that FFT_LIMB_SIZE > BITS_PER_GROUP/2.
// At this moment FFT_LIMB_SIZE = 19 and BITS_PER_GROUP = 31.
static int ReduceLimbs(const limb *factor, struct sComplex *fftFactor, int len)
{
  size_t diffPtrs;
  int bitExternal = 0;  // Least significant bit of current external limb
                        // that corresponds to bit zero of FFT internal limb.
  const limb *ptrFactor = factor;
  struct sComplex *ptrInternalFactor = fftFactor;
  for (;;)
  {
    unsigned int uBitExternal = (unsigned int)bitExternal;
    int maxValueFFTLimb;
    int real = UintToInt((unsigned int)ptrFactor->x >> uBitExternal);
    int imaginary;
    diffPtrs = ptrFactor - factor;
    if ((int)diffPtrs < (len - 1))
    {                   // Do not read outside input buffer.
      unsigned int complementBitExternal = (unsigned int)BITS_PER_GROUP - uBitExternal;
      real += UintToInt((unsigned int)(ptrFactor + 1)->x << complementBitExternal);
    }
    maxValueFFTLimb = MAX_VALUE_FFT_LIMB;
    real &= (unsigned int)maxValueFFTLimb;
    ptrInternalFactor->real = (double)real;
    bitExternal += FFT_LIMB_SIZE;
    if (bitExternal >= BITS_PER_GROUP)
    {                   // All bits of input limb have been used.
      bitExternal -= BITS_PER_GROUP;
      ptrFactor++;
      if ((ptrFactor - factor) == len)
      {
        ptrInternalFactor->imaginary = 0.0;
        ptrInternalFactor++;
        break;
      }
    }
    uBitExternal = (unsigned int)bitExternal;
    imaginary = UintToInt((unsigned int)ptrFactor->x >> uBitExternal);
    diffPtrs = ptrFactor - factor;
    if ((int)diffPtrs < (len - 1))
    {                   // Do not read outside input buffer.
      unsigned int complementBitExternal = (unsigned int)BITS_PER_GROUP - uBitExternal;
      imaginary += UintToInt((unsigned int)(ptrFactor + 1)->x << complementBitExternal);
    }
    imaginary &= MAX_VALUE_FFT_LIMB;
    ptrInternalFactor->imaginary = (double)imaginary;
    ptrInternalFactor++;
    bitExternal += FFT_LIMB_SIZE;
    if (bitExternal >= BITS_PER_GROUP)
    {                   // All bits of input limb have been used.
      bitExternal -= BITS_PER_GROUP;
      ptrFactor++;
      diffPtrs = ptrFactor - factor;
      if ((int)diffPtrs == len)
      {
        break;
      }
    }
  }
  diffPtrs = ptrInternalFactor - fftFactor;
  return (int)diffPtrs;
}

#if 0
   Algorithm 9.5.12 of Crandall and Pomerance book Prime Numbers:

   Zero-pad x and y until each has length 2D.
   X = DFT(x);
   Y = DFT(y);
   Z = X * Y;      // Using convolution
   z = DFT^(-1)(Z)
   z = round(z)    // Round elementwise
   carry = 0;
   for (0 <= n < 2D)
   {
     v = z_n + carry;
     z_n = v mod B
     carry = floor(v/B)
   }
   Delete leading zeros.
#endif

void fftMultiplication(const limb *factor1, const limb *factor2, limb *result,
  int len1, int len2, int *pResultLen)
{
  const struct sComplex *ptrFirst;
  const struct sComplex *ptrSecond;
  struct sComplex *ptrProduct;
  double invPower2;
  double dCarry;
  int fftLen;
  int fftLen1;
  int fftLen2 = 0;
  unsigned int bitExternal;
  int power2plus1;
  int power2plus1Bytes;
  limb *ptrResult;
  int power2 = 1;
  int index;
  int sumLen;
  int sumLenBytes;
  fftLen1 = ReduceLimbs(factor1, firstFactor, len1);
  if ((factor1 != factor2) && !((TestNbrCached == NBR_CACHED) && (factor2 == TestNbr)) &&
    !((MontgomeryMultNCached == NBR_CACHED) && (factor2 == MontgomeryMultN)))
  {
    fftLen2 = ReduceLimbs(factor2, secondFactor, len2);
  }
  // Get next power of 2 to len (the greatest of the two lengths).
  fftLen = fftLen1;
  if (fftLen < fftLen2)
  {
    fftLen = fftLen2;
  }
  while (power2 < fftLen)
  {
    power2 *= 2;
  }
  power2 += power2;
  for (index = fftLen1; index < power2; index++)
  {
    firstFactor[index].real = 0;
    firstFactor[index].imaginary = 0;
  }
  for (index = fftLen2; index < power2; index++)
  {
    secondFactor[index].real = 0;
    secondFactor[index].imaginary = 0;
  }
  complexFFT(firstFactor, tempFFT, power2);
  ConvertHalfToFullSizeFFT(tempFFT, product, power2);   // product <- DFT(firstFactor)
  power2plus1 = power2 + 1;
  power2plus1Bytes = power2plus1 * (int)sizeof(transf[0]);
  if (factor1 != factor2)
  {
    if ((TestNbrCached == NBR_CACHED) && (factor2 == TestNbr))
    {
      (void)memcpy(transf, TestNbrTransf, power2plus1Bytes);
    }
    else if ((MontgomeryMultNCached == NBR_CACHED) && (factor2 == MontgomeryMultN))
    {
      (void)memcpy(transf, MontgomeryMultNTransf, power2plus1Bytes);
    }
    else
    {
      complexFFT(secondFactor, tempFFT, power2);
      ConvertHalfToFullSizeFFT(tempFFT, transf, power2);  // transf <- DFT(secondFactor)
    }
    if ((TestNbrCached == NBR_READY_TO_BE_CACHED) && (factor2 == TestNbr))
    {
      (void)memcpy(TestNbrTransf, transf, power2plus1Bytes);
      TestNbrCached = NBR_CACHED;
    }
    else if ((MontgomeryMultNCached == NBR_READY_TO_BE_CACHED) && (factor2 == MontgomeryMultN))
    {
      (void)memcpy(MontgomeryMultNTransf, transf, power2plus1Bytes);
      MontgomeryMultNCached = NBR_CACHED;
    }
    else
    {  // No more conditions.
    }
  }
  else
  {
    (void)memcpy(transf, product, power2plus1Bytes);   // transf <- DFT(secondFactor)
  }

    // Perform convolution.
  ptrFirst = product;
  ptrSecond = transf;
  ptrProduct = product;
  for (index = 0; index <= power2; index++)
  {
    double real = (ptrFirst->real*ptrSecond->real) - 
      (ptrFirst->imaginary*ptrSecond->imaginary);
    ptrProduct->imaginary = (ptrFirst->real*ptrSecond->imaginary) +
      (ptrFirst->imaginary*ptrSecond->real);
    ptrProduct->real = real;
    ptrFirst++;
    ptrSecond++;
    ptrProduct++;
  }
  ConvertFullToHalfSizeFFT(product, tempFFT, power2);
  // Perform inverse DFT of product.
  complexFFT(tempFFT, transf, power2);
  ptrProduct = transf;
  invPower2 = 0.125 / (double)power2;
  dCarry = 0;
  sumLen = len1 + len2;
  sumLenBytes = sumLen * (int)sizeof(limb);
  (void)memset(result, 0, sumLenBytes);
  bitExternal = 0U;
  ptrResult = result;
  for (index = 0; index < power2; index++)
  {
    double dQuot;
    int fftResult;

    // Real part.
    dCarry += floor((ptrProduct->real * invPower2) + 0.5);
    dQuot = floor(dCarry / (double)FFT_LIMB_RANGE);
    fftResult = (int)(dCarry - (dQuot * (double)FFT_LIMB_RANGE));
    ptrResult->x |= UintToInt(((unsigned int)fftResult << bitExternal) & MAX_VALUE_LIMB);
    if ((int)bitExternal > (BITS_PER_GROUP - FFT_LIMB_SIZE))
    {
      unsigned int shiftRight = (unsigned int)BITS_PER_GROUP - bitExternal;
      (ptrResult+1)->x |= UintToInt(((unsigned int)fftResult >> shiftRight) & MAX_VALUE_LIMB);
    }
    bitExternal += (unsigned int)FFT_LIMB_SIZE;
    if (bitExternal >= (unsigned int)BITS_PER_GROUP)
    {
      bitExternal -= (unsigned int)BITS_PER_GROUP;
      ptrResult++;
      if ((ptrResult - result) == sumLen)
      {
        break;
      }
    }
    dCarry = dQuot;

    // Imaginary part. Use negative value for inverse FFT.
    dCarry += floor((-ptrProduct->imaginary * invPower2) + 0.5);
    ptrProduct++;
    dQuot = floor(dCarry / (double)FFT_LIMB_RANGE);
    fftResult = (int)(dCarry - (dQuot * (double)FFT_LIMB_RANGE));
    ptrResult->x |= UintToInt(((unsigned int)fftResult << bitExternal) & MAX_VALUE_LIMB);
    if ((int)bitExternal > (BITS_PER_GROUP - FFT_LIMB_SIZE))
    {
      unsigned int shRight = (unsigned int)BITS_PER_GROUP - bitExternal;
      (ptrResult + 1)->x |= UintToInt(((unsigned int)fftResult >> shRight) & MAX_VALUE_LIMB);
    }
    bitExternal += (unsigned int)FFT_LIMB_SIZE;
    if (bitExternal >= (unsigned int)BITS_PER_GROUP)
    {
      bitExternal -= (unsigned int)BITS_PER_GROUP;
      ptrResult++;
      if ((ptrResult - result) == sumLen)
      {
        break;
      }
    }
    dCarry = dQuot;
  }
  if (pResultLen != NULL)
  {
    while ((sumLen > 1) && ((result + sumLen - 1)->x == 0))
    {
      sumLen--;
    }
    *pResultLen = sumLen;
  }
}

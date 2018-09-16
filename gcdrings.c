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

#include "bignbr.h"

struct approx
{
  int expon;
  double mantissa;
};

static struct approx appMin;

// Operations with complex numbers
static void MultiplyComplexBy2(BigInteger *real, BigInteger *imag)
{
  BigIntMultiplyBy2(real);
  BigIntMultiplyBy2(imag);
}

static void DivideComplexBy2(BigInteger *real, BigInteger *imag)
{
  BigIntDivideBy2(real);
  BigIntDivideBy2(imag);
}
static void DivideComplexBy1PlusI(BigInteger *real, BigInteger *imag, BigInteger *temp)
{        // Multiply by 1-i and then divide by 2.
  BigIntAdd(real, imag, temp);
  BigIntSubt(imag, real, imag);
  CopyBigInt(real, temp);
  DivideComplexBy2(real, imag);
}

static void MultiplyComplexBy1PlusI(BigInteger *real, BigInteger *imag, BigInteger *temp)
{        // Multiply by 1+i.
  BigIntSubt(real, imag, temp);
  BigIntAdd(real, imag, imag);
  CopyBigInt(real, temp);
}

// Use algorithm "Variants of an algorithm of J. Stein" of Sándor Szabó.
void GaussianGCD(BigInteger *realA, BigInteger *imagA, BigInteger *realB, BigInteger *imagB,
  BigInteger *realGcd, BigInteger *imagGcd, BigInteger *temp1, BigInteger *temp2)
{
  intToBigInteger(realGcd, 1);     // Initialize d to 1.
  intToBigInteger(imagGcd, 0);
  for (;;)
  {
    int GcdIsAtimesD = FALSE;
    int type;
    if (BigIntEqual(realA, realB) && BigIntEqual(imagA, imagB))
    {         // A and B are associates.
      GcdIsAtimesD = TRUE;
    }
    BigIntChSign(realB);
    BigIntChSign(imagB);
    if (BigIntEqual(realA, realB) && BigIntEqual(imagA, imagB))
    {         // A and B are associates.
      GcdIsAtimesD = TRUE;
    }
    BigIntChSign(realB);
    BigIntChSign(imagB);
    if (BigIntIsZero(realB) && BigIntIsZero(imagB))
    {
      GcdIsAtimesD = TRUE;
    }
    if (GcdIsAtimesD)
    {        // The GCD is d*a
      BigIntMultiply(realA, realGcd, temp1);
      BigIntMultiply(imagA, imagGcd, temp2);
      BigIntSubt(temp1, temp2, realB);
      BigIntMultiply(realA, imagGcd, temp1);
      BigIntMultiply(imagA, realGcd, temp2);
      BigIntAdd(temp1, temp2, imagGcd);
      CopyBigInt(realGcd, realB);
      return;
    }
    if (realB->nbrLimbs == 1 && imagB->nbrLimbs == 1 && realB->limbs[0].x + imagB->limbs[0].x == 1)
    {        // B is a unit. GCD already found in d.
      return;
    }
    type = (realA->limbs[0].x & 1 ? 8 : 0) + (imagA->limbs[0].x & 1 ? 4 : 0) +
      (realB->limbs[0].x & 1 ? 2 : 0) + (imagB->limbs[0].x & 1 ? 1 : 0);
    switch (type)
    {
    case 0:      // all components even.
      DivideComplexBy2(realA, imagA);
      DivideComplexBy2(realB, imagB);
      MultiplyComplexBy2(realGcd, imagGcd);
      break;
    case 1:      // real A even, imag A even, real B even, imag B odd.
    case 2:      // real A even, imag A even, real B odd, imag B even.
      DivideComplexBy2(realA, imagA);
      break;
    case 3:      // real A even, imag A even, real B odd, imag B odd.
      DivideComplexBy2(realA, imagA);
      DivideComplexBy1PlusI(realB, imagB, temp1);
      MultiplyComplexBy1PlusI(realGcd, imagGcd, temp1);
      break;
    case 4:      // real A even, imag A odd, real B even, imag B even.
    case 8:      // real A odd, imag A even, real B even, imag B even.
      DivideComplexBy2(realB, imagB);
      break;
    case 5:      // real A even, imag A odd, real B even, imag B odd.
    case 10:     // real A odd, imag A even, real B odd, imag B even.
      BigIntAdd(realA, realB, temp1);       // A_{n+1} <- (A_n + B_n)/2
      BigIntSubt(realA, realB, realB);      // B_{n+1} <- (A_n - B_n)/2
      CopyBigInt(realA, temp1);
      BigIntAdd(imagA, imagB, temp1);
      BigIntSubt(imagA, imagB, imagB);
      CopyBigInt(imagA, temp1);
      DivideComplexBy2(realA, imagA);
      DivideComplexBy2(realB, imagB);
      break;
    case 6:      // real A even, imag A odd, real B odd, imag B even.
    case 9:      // real A odd, imag A even, real B even, imag B odd.
      CopyBigInt(temp1, realB);              // B_{n+1} <- B_n * i.
      CopyBigInt(realB, imagB);
      CopyBigInt(imagB, temp1);
      BigIntChSign(realB);
      break;
    case 7:      // real A even, imag A odd, real B odd, imag B odd.
    case 11:     // real A odd, imag A even, real B odd, imag B odd.
      DivideComplexBy1PlusI(realB, imagB, temp1);
      break;
    case 12:     // real A odd, imag A odd, real B even, imag B even.
      DivideComplexBy1PlusI(realA, imagA, temp1);
      DivideComplexBy2(realB, imagB);
      MultiplyComplexBy1PlusI(realGcd, imagGcd, temp1);
      break;
    case 13:     // real A odd, imag A odd, real B even, imag B odd.
    case 14:     // real A odd, imag A odd, real B odd, imag B even.
      DivideComplexBy1PlusI(realA, imagA, temp1);
      break;
    default:     // real A odd, imag A odd, real B odd, imag B odd.
      DivideComplexBy1PlusI(realA, imagA, temp1);
      DivideComplexBy1PlusI(realB, imagB, temp1);
      MultiplyComplexBy1PlusI(realGcd, imagGcd, temp1);
      break;
    }
  }
}

// Operations with quaternions
void MultiplyQuaternionBy2(BigInteger *scalar, BigInteger *vecI, BigInteger *vecJ, BigInteger *vecK)
{
  BigIntMultiplyBy2(scalar);
  BigIntMultiplyBy2(vecI);
  BigIntMultiplyBy2(vecJ);
  BigIntMultiplyBy2(vecK);
}

void DivideQuaternionBy2(BigInteger *scalar, BigInteger *vecI, BigInteger *vecJ, BigInteger *vecK)
{
  BigIntDivideBy2(scalar);
  BigIntDivideBy2(vecI);
  BigIntDivideBy2(vecJ);
  BigIntDivideBy2(vecK);
}
static void DivideQuaternionBy1PlusI(BigInteger *scalar, BigInteger *vecI, BigInteger *vecJ, BigInteger *vecK, BigInteger *temp)
{        // Multiply by 1-i and then divide by 2.
  BigIntAdd(vecI, scalar, temp);
  BigIntSubt(vecI, scalar, vecI);
  CopyBigInt(scalar, temp);
  BigIntAdd(vecK, vecJ, temp);
  BigIntSubt(vecK, vecJ, vecK);
  CopyBigInt(vecJ, temp);
  DivideQuaternionBy2(scalar, vecI, vecJ, vecK);
}

static void MultiplyQuaternionBy1PlusI(BigInteger *scalar, BigInteger *vecI, BigInteger *vecJ, BigInteger *vecK, BigInteger *temp)
{        // Multiply by 1+i.
  BigIntSubt(scalar, vecI, temp);
  BigIntAdd(scalar, vecI, vecI);
  CopyBigInt(scalar, temp);
  BigIntSubt(vecJ, vecK, temp);
  BigIntAdd(vecJ, vecK, vecK);
  CopyBigInt(vecJ, temp);
}

static int DivideQuaternionByPowerOf1PlusI(BigInteger *scalar, BigInteger *vecI, BigInteger *vecJ, BigInteger *vecK, BigInteger *temp)
{
  int expon = 0;
  while ((scalar->limbs[0].x & 1) == 0)
  {
    int sumCoeff = (scalar->limbs[0].x & 3) + (vecI->limbs[0].x & 3) +
      (vecJ->limbs[0].x & 3) + (vecK->limbs[0].x & 3);
    if (sumCoeff == 0)
    {     // Number is multiple of 2.
      expon += 2;
      DivideQuaternionBy2(scalar, vecI, vecJ, vecK);
    }
    else if (sumCoeff == 4)
    {     // Number is multiple of 1+i.
      expon++;
      DivideQuaternionBy1PlusI(scalar, vecI, vecJ, vecK, temp);
    }
    else
    {     // Number not multiple of 2 or 1+i. Go out.
      break;
    }
  }
  return expon;
}

#define Mod4(bignbr) ((bignbr->limbs[0].x & 3) ^ (bignbr->sign << 1))
// If type of quaternion is odd, if the number of '1' is even, multiply it by (1/2)(1+i+j+k).
// If the number of '3' is odd, multiply it by (1/2)(-1+i+j+k).
static void ConvertQuaternionToEvenType(BigInteger *scalar, BigInteger *vecI, BigInteger *vecJ, BigInteger *vecK,
  BigInteger *temp1, BigInteger *temp2, BigInteger *temp3)
{
  int sumCoeff;
  if ((scalar->limbs[0].x & 1) == 0)
  {
    return;    // Type is even. No need to change quaternion.
  }

  sumCoeff = Mod4(scalar) + Mod4(vecI) + Mod4(vecJ) + Mod4(vecK);
  if (sumCoeff == 6 || sumCoeff == 10)
  {     // Odd number of '3' in coefficients of A. Multiply by (1/2)(1+i+j+k)
    // Compute scalar part of the product: d_1*a_1 - d_i*a_i - d_j*a_j - d_k*a_k.
    BigIntSubt(scalar, vecI, temp1);
    BigIntSubt(temp1, vecJ, temp1);
    BigIntSubt(temp1, vecK, temp1);
    // Compute coefficient of "i" of the product: d_1*a_i + d_i*a_1 + d_j*a_k - d_k*a_j.
    BigIntAdd(scalar, vecI, temp2);
    BigIntSubt(temp2, vecJ, temp2);
    BigIntAdd(temp2, vecK, temp2);
    // Compute coefficient of "j" of the product: d_1*a_j - d_i*a_k + d_j*a_1 + d_k*a_i.
    BigIntAdd(scalar, vecI, temp3);
    BigIntAdd(temp3, vecJ, temp3);
    BigIntSubt(temp3, vecK, temp3);
    // Compute coefficient of "k" of the product: d_1*a_k + d_i*a_j - d_j*a_i + d_k*a_1.
    BigIntAdd(vecK, scalar, vecK);
    BigIntSubt(vecK, vecI, vecK);
    BigIntAdd(vecK, vecJ, vecK);
  }
  else
  {     // Even number of '3' in coefficients of A. Multiply by (1/2)(-1+i+j+k)
    // Compute scalar part of the product: d_1*a_1 - d_i*a_i - d_j*a_j - d_k*a_k.
    BigIntAdd(scalar, vecI, temp1);
    BigIntAdd(temp1, vecJ, temp1);
    BigIntAdd(temp1, vecK, temp1);
    BigIntChSign(temp1);
    // Compute coefficient of "i" of the product: d_1*a_i + d_i*a_1 + d_j*a_k - d_k*a_j.
    BigIntSubt(scalar, vecI, temp2);
    BigIntSubt(temp2, vecJ, temp2);
    BigIntAdd(temp2, vecK, temp2);
    // Compute coefficient of "j" of the product: d_1*a_j - d_i*a_k + d_j*a_1 + d_k*a_i.
    BigIntAdd(scalar, vecI, temp3);
    BigIntSubt(temp3, vecJ, temp3);
    BigIntSubt(temp3, vecK, temp3);
    // Compute coefficient of "k" of the product: d_1*a_k + d_i*a_j - d_j*a_i + d_k*a_1.
    BigIntSubt(scalar, vecK, vecK);
    BigIntSubt(vecK, vecI, vecK);
    BigIntAdd(vecK, vecJ, vecK);
  }
  CopyBigInt(scalar, temp1);
  CopyBigInt(vecI, temp2);
  CopyBigInt(vecJ, temp3);
  DivideQuaternionBy2(scalar, vecI, vecJ, vecK);
}

static void approximate(BigInteger *nbr, struct approx *appNbr)
{
  int nbrLimbs = nbr->nbrLimbs;
  appNbr->expon = nbrLimbs;
  appNbr->mantissa = (double)nbr->limbs[nbrLimbs - 1].x * (double)LIMB_RANGE;
  if (nbrLimbs > 1)
  {
    appNbr->mantissa += (double)nbr->limbs[nbrLimbs - 2].x;
  }
  if (nbr->sign == SIGN_NEGATIVE)
  {
    appNbr->mantissa = -appNbr->mantissa;
  }
}

static void addApprox(struct approx *addend1, struct approx *addend2, struct approx *sum)
{
  int nbrLimbs1 = addend1->expon;
  int nbrLimbs2 = addend2->expon;
  int iTemp;
  double mantissa1 = addend1->mantissa;
  double mantissa2 = addend2->mantissa;
  double dTemp;
  if (nbrLimbs1 < nbrLimbs2)
  {   // Exchange number of limbs and mantissas.
    iTemp = nbrLimbs1;
    nbrLimbs1 = nbrLimbs2;
    nbrLimbs2 = iTemp;
    dTemp = mantissa1;
    mantissa1 = mantissa2;
    mantissa2 = dTemp;
  }
  if (nbrLimbs1 == nbrLimbs2)
  {
    mantissa1 += mantissa2;
  }
  else if (nbrLimbs1 == nbrLimbs2 + 1)
  {
    mantissa1 += mantissa2 / (double)LIMB_RANGE;
  }
  if (mantissa1 >= (double)LIMB_RANGE * (double)LIMB_RANGE ||
    mantissa1 <= -(double)LIMB_RANGE * (double)LIMB_RANGE)
  {
    mantissa1 /= (double)LIMB_RANGE;
    nbrLimbs1++;
  }
  if (nbrLimbs1 > 1 && mantissa1 < (double)LIMB_RANGE && mantissa1 > -(double)LIMB_RANGE)
  {
    nbrLimbs1--;
    mantissa1 *= (double)LIMB_RANGE;
  }
  sum->expon = nbrLimbs1;
  sum->mantissa = mantissa1;
}

static void squareApprox(struct approx *value, struct approx *square)
{
  int nbrLimbs = value->expon + value->expon - 1;
  double mantissa = value->mantissa * value->mantissa / (double)LIMB_RANGE;
  while (mantissa >= (double)LIMB_RANGE * (double)LIMB_RANGE)
  {                // mantissa is positive after squaring it.
    nbrLimbs++;
    mantissa /= (double)LIMB_RANGE;
  }
  square->expon = nbrLimbs;
  square->mantissa = mantissa;
}

// let sum = |scalarA + scalarB| + |vecIA + vecIB| + |vecJA + vecJB| + |vecKA + vecKB|
// if sum < min, set *pCombination to nbr and set min to sum.
// let sum = |scalarA - scalarB| + |vecIA - vecIB| + |vecJA - vecJB| + |vecKA - vecKB|
// if sum < min, set *pCombination to nbr+1 and set min to sum.
static void TestCombination(BigInteger *scalarA, BigInteger *vecIA, BigInteger *vecJA, BigInteger *vecKA,
  BigInteger *scalarB, BigInteger *vecIB, BigInteger *vecJB, BigInteger *vecKB,
  int nbr, int *pCombination)
{
  int ctr;
  struct approx appScalarA, appVecIA, appVecJA, appVecKA;
  struct approx appScalarB, appVecIB, appVecJB, appVecKB;
  struct approx sum, temp;
  approximate(scalarA, &appScalarA);
  approximate(vecIA, &appVecIA);
  approximate(vecJA, &appVecJA);
  approximate(vecKA, &appVecKA);
  approximate(scalarB, &appScalarB);
  approximate(vecIB, &appVecIB);
  approximate(vecJB, &appVecJB);
  approximate(vecKB, &appVecKB);
  for (ctr = 0; ctr < 2; ctr++)
  {
    addApprox(&appScalarA, &appScalarB, &sum);
    squareApprox(&sum, &sum);
    addApprox(&appVecIA, &appVecIB, &temp);
    squareApprox(&temp, &temp);
    addApprox(&sum, &temp, &sum);
    addApprox(&appVecJA, &appVecJB, &temp);
    squareApprox(&temp, &temp);
    addApprox(&sum, &temp, &sum);
    addApprox(&appVecKA, &appVecKB, &temp);
    squareApprox(&temp, &temp);
    addApprox(&sum, &temp, &sum);
    if ((nbr == 0 && ctr == 0) ||
      (sum.expon < appMin.expon || (sum.expon == appMin.expon && sum.mantissa < appMin.mantissa)))
    {
      appMin.expon = sum.expon;
      appMin.mantissa = sum.mantissa;
      *pCombination = nbr + ctr;
    }
    appScalarB.mantissa = -appScalarB.mantissa;
    appVecIB.mantissa = -appVecIB.mantissa;
    appVecJB.mantissa = -appVecJB.mantissa;
    appVecKB.mantissa = -appVecKB.mantissa;
  }
}

// If *pCombination > 1, subtract 2 to it and go out.
// If *pCombination = 0, perform A <- A + B.
// If *pCombination = 1, perform A <- A - B.
static void AddCombination(BigInteger *scalarA, BigInteger *vecIA, BigInteger *vecJA, BigInteger *vecKA,
  BigInteger *scalarB, BigInteger *vecIB, BigInteger *vecJB, BigInteger *vecKB,
  int *pCombination)
{
  int combination = *pCombination;
  *pCombination -= 2;
  if (combination == 0)
  {
    BigIntAdd(scalarA, scalarB, scalarA);
    BigIntAdd(vecIA, vecIB, vecIA);
    BigIntAdd(vecJA, vecJB, vecJA);
    BigIntAdd(vecKA, vecKB, vecKA);
  }
  else if (combination == 1)
  {
    BigIntSubt(scalarA, scalarB, scalarA);
    BigIntSubt(vecIA, vecIB, vecIA);
    BigIntSubt(vecJA, vecJB, vecJA);
    BigIntSubt(vecKA, vecKB, vecKA);
  }
}

// Compute sum as the sum of the squares of all components.
static void getApproxHeight(BigInteger *scalar, BigInteger *vecI, BigInteger *vecJ, BigInteger *vecK,
  struct approx *sum)
{
  struct approx appTemp;
  approximate(scalar, sum);
  squareApprox(sum, sum);
  approximate(vecI, &appTemp);
  squareApprox(&appTemp, &appTemp);
  addApprox(sum, &appTemp, sum);
  approximate(vecJ, &appTemp);
  squareApprox(&appTemp, &appTemp);
  addApprox(sum, &appTemp, sum);
  approximate(vecK, &appTemp);
  squareApprox(&appTemp, &appTemp);
  addApprox(sum, &appTemp, sum);
}

// Use algorithm "Variants of an algorithm of J. Stein" of Sándor Szabó.
// Since integer quaternions can have coefficients that are half of an integer, all
// cofficients in this routine are multiplied by 2.
void QuaternionGCD(BigInteger *scalarA, BigInteger *vecIA, BigInteger *vecJA, BigInteger *vecKA,
  BigInteger *scalarB, BigInteger *vecIB, BigInteger *vecJB, BigInteger *vecKB,
  BigInteger *scalarGcd, BigInteger *vecIGcd, BigInteger *vecJGcd, BigInteger *vecKGcd,
  BigInteger *temp1, BigInteger *temp2, BigInteger *temp3, BigInteger *temp4)
{
  int combination;
  struct approx sumA, sumB;
  intToBigInteger(scalarGcd, 1);     // Initialize d to 1.
  intToBigInteger(vecIGcd, 0);
  intToBigInteger(vecJGcd, 0);
  intToBigInteger(vecKGcd, 0);
  for (;;)
  {
    // if type of quaternions are odd, convert them to even by multiplying by units.
    ConvertQuaternionToEvenType(scalarA, vecIA, vecJA, vecKA, temp1, temp2, temp3);
    ConvertQuaternionToEvenType(scalarB, vecIB, vecJB, vecKB, temp1, temp2, temp3);
    int GcdIsAtimesD = FALSE;
    int exponA, exponB;
    if (BigIntEqual(scalarA, scalarB) && BigIntEqual(vecIA, vecIB) &&
      BigIntEqual(vecJA, vecJB) && BigIntEqual(vecKA, vecKB))
    {         // A and B are associates.
      GcdIsAtimesD = TRUE;
    }
    BigIntChSign(scalarB);
    BigIntChSign(vecIB);
    BigIntChSign(vecJB);
    BigIntChSign(vecKB);
    if (BigIntEqual(scalarA, scalarB) && BigIntEqual(vecIA, vecIB) &&
      BigIntEqual(vecJA, vecJB) && BigIntEqual(vecKA, vecKB))
    {         // A and B are associates.
      GcdIsAtimesD = TRUE;
    }
    BigIntChSign(scalarB);
    BigIntChSign(vecIB);
    BigIntChSign(vecJB);
    BigIntChSign(vecKB);
    if (BigIntIsZero(scalarB) && BigIntIsZero(vecIB) && BigIntIsZero(vecJB) && BigIntIsZero(vecKB))
    {
      GcdIsAtimesD = TRUE;
    }
    if (BigIntIsZero(scalarA) && BigIntIsZero(vecIA) && BigIntIsZero(vecJA) && BigIntIsZero(vecKA))
    {
      CopyBigInt(scalarA, scalarB);
      CopyBigInt(vecIA, vecIB);
      CopyBigInt(vecJA, vecJB);
      CopyBigInt(vecKA, vecKB);
      GcdIsAtimesD = TRUE;
    }
    if (GcdIsAtimesD)
    {        // The GCD is d*a
      // Compute scalar part of the product: d_1*a_1 - d_i*a_i - d_j*a_j - d_k*a_k.
      BigIntMultiply(scalarGcd, scalarA, temp1);
      BigIntMultiply(vecIGcd, vecIA, temp4);
      BigIntSubt(temp1, temp4, temp1);
      BigIntMultiply(vecJGcd, vecJA, temp4);
      BigIntSubt(temp1, temp4, temp1);
      BigIntMultiply(vecKGcd, vecKA, temp4);
      BigIntSubt(temp1, temp4, temp1);

      // Compute coefficient of "i" of the product: d_1*a_i + d_i*a_1 + d_j*a_k - d_k*a_j.
      BigIntMultiply(scalarGcd, vecIA, temp2);
      BigIntMultiply(vecIGcd, scalarA, temp4);
      BigIntAdd(temp2, temp4, temp2);
      BigIntMultiply(vecJGcd, vecKA, temp4);
      BigIntAdd(temp2, temp4, temp2);
      BigIntMultiply(vecKGcd, vecJA, temp4);
      BigIntSubt(temp2, temp4, temp2);

      // Compute coefficient of "j" of the product: d_1*a_j - d_i*a_k + d_j*a_1 + d_k*a_i.
      BigIntMultiply(scalarGcd, vecJA, temp3);
      BigIntMultiply(vecIGcd, vecKA, temp4);
      BigIntSubt(temp3, temp4, temp3);
      BigIntMultiply(vecJGcd, scalarA, temp4);
      BigIntAdd(temp3, temp4, temp3);
      BigIntMultiply(vecKGcd, vecIA, temp4);
      BigIntAdd(temp3, temp4, temp3);

      // Compute coefficient of "k" of the product: d_1*a_k + d_i*a_j - d_j*a_i + d_k*a_1.
      BigIntMultiply(vecKGcd, scalarA, vecKGcd);
      BigIntMultiply(scalarGcd, vecKA, temp4);
      BigIntAdd(vecKGcd, temp4, vecKGcd);
      BigIntMultiply(vecIGcd, vecJA, temp4);
      BigIntSubt(vecKGcd, temp4, vecKGcd);
      BigIntMultiply(vecJGcd, vecIA, temp4);
      BigIntAdd(vecKGcd, temp4, vecKGcd);

      CopyBigInt(scalarGcd, temp1);
      CopyBigInt(vecIGcd, temp2);
      CopyBigInt(vecJGcd, temp3);
      ConvertQuaternionToEvenType(scalarGcd, vecIGcd, vecJGcd, vecKGcd, temp1, temp2, temp3);
      return;
    }
    if (scalarB->nbrLimbs == 1 && vecIB->nbrLimbs == 1 && vecJB->nbrLimbs == 1 && vecKB->nbrLimbs == 1)
    {
      int magnitude = scalarB->limbs[0].x | vecIB->limbs[0].x |
        vecJB->limbs[0].x | vecKB->limbs[0].x;
      if (scalarB->limbs[0].x & 1)
      {          // If the absolute values of all four coefficients equal 1, B is a unit.
        if (magnitude == 1)
        {        // B is a unit. GCD already found in d.
          ConvertQuaternionToEvenType(scalarGcd, vecIGcd, vecJGcd, vecKGcd, temp1, temp2, temp3);
          return;
        }
      }
      else
      {          // If the absolute value of one coefficient is 2 and all other coefficients
                 // equal zero, B is a unit.
        if (magnitude == 2)
        {
          if ((scalarB->limbs[0].x + vecIB->limbs[0].x +
            vecJB->limbs[0].x + vecKB->limbs[0].x) == 2)
          {        // B is a unit. GCD already found in d.
            ConvertQuaternionToEvenType(scalarGcd, vecIGcd, vecJGcd, vecKGcd, temp1, temp2, temp3);
            return;
          }
        }
      }
    }
    // Check whether the type of A is even.
    exponA = DivideQuaternionByPowerOf1PlusI(scalarA, vecIA, vecJA, vecKA, temp1);
    exponB = DivideQuaternionByPowerOf1PlusI(scalarB, vecIB, vecJB, vecKB, temp1);
    if (exponA && exponB)
    {
      // Set in exponA the power of 1+i to be multiplied to the GCD.
      if (exponA > exponB)
      {
        exponA = exponB;
      }
      for (exponB = 0; exponB < exponA; exponB += 2)
      {
        MultiplyQuaternionBy2(scalarGcd, vecIGcd, vecJGcd, vecKGcd);
      }
      if (exponA & 1)
      {
        MultiplyQuaternionBy1PlusI(scalarGcd, vecIGcd, vecJGcd, vecKGcd, temp1);
      }
      continue;
    }
    if (exponA || exponB)
    {
      continue;
    }
    // Let q(a) be the sum of the absolute values of the components of the quaternion a.
    // If q(a) < q(b), exchange a and b.
    // x = minimum of q(a + u*b) where u is one of +/-1, +/-i, +/-j, +/-k (8 combinations).
    // Replace a by x.
    getApproxHeight(scalarA, vecIA, vecJA, vecKA, &sumA);
    getApproxHeight(scalarB, vecIB, vecJB, vecKB, &sumB);
    if (sumA.expon < sumB.expon || (sumA.expon == sumB.expon && sumA.mantissa < sumB.mantissa))
    {    // Exchange A and B.
      CopyBigInt(temp1, scalarA);
      CopyBigInt(scalarA, scalarB);
      CopyBigInt(scalarB, temp1);
      CopyBigInt(temp1, vecIA);
      CopyBigInt(vecIA, vecIB);
      CopyBigInt(vecIB, temp1);
      CopyBigInt(temp1, vecJA);
      CopyBigInt(vecJA, vecJB);
      CopyBigInt(vecJB, temp1);
      CopyBigInt(temp1, vecKA);
      CopyBigInt(vecKA, vecKB);
      CopyBigInt(vecKB, temp1);
    }
    combination = -1;
    // Test a +/- b
    TestCombination(scalarA, vecIA, vecJA, vecKA, scalarB, vecIB, vecJB, vecKB, 0, &combination);
    // Test a +/- i*b  i and k coefficients must be negative.
    BigIntChSign(vecIB);
    BigIntChSign(vecKB);
    TestCombination(scalarA, vecIA, vecJA, vecKA, vecIB, scalarB, vecKB, vecJB, 2, &combination);
    BigIntChSign(vecKB);
    // Test a +/- j*b  i and j coefficients must be negative.
    BigIntChSign(vecJB);
    TestCombination(scalarA, vecIA, vecJA, vecKA, vecJB, vecKB, scalarB, vecIB, 4, &combination);
    BigIntChSign(vecIB);
    // Test a +/- k*b  j and k coefficients must be negative.
    BigIntChSign(vecKB);
    TestCombination(scalarA, vecIA, vecJA, vecKA, vecKB, vecJB, vecIB, scalarB, 6, &combination);
    BigIntChSign(vecKB);
    BigIntChSign(vecJB);

    // Perform a +/- b if minimum combination
    AddCombination(scalarA, vecIA, vecJA, vecKA, scalarB, vecIB, vecJB, vecKB, &combination);
    // Perform a +/- i*b if minimum combination
    BigIntChSign(vecIB);
    BigIntChSign(vecKB);
    AddCombination(scalarA, vecIA, vecJA, vecKA, vecIB, scalarB, vecKB, vecJB, &combination);
    BigIntChSign(vecKB);
    // Perform a +/- j*b if minimum combination
    BigIntChSign(vecJB);
    AddCombination(scalarA, vecIA, vecJA, vecKA, vecJB, vecKB, scalarB, vecIB, &combination);
    BigIntChSign(vecIB);
    // Perform a +/- k*b if minimum combination
    BigIntChSign(vecKB);
    AddCombination(scalarA, vecIA, vecJA, vecKA, vecKB, vecJB, vecIB, scalarB, &combination);
    BigIntChSign(vecKB);
    BigIntChSign(vecJB);
  }
}


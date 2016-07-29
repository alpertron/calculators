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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bignbr.h"
limb MontgomeryR1[MAX_LEN];
limb TestNbr[MAX_LEN];
limb MontgomeryMultN[MAX_LEN];
limb MontgomeryMultR1[MAX_LEN];
limb MontgomeryMultR2[MAX_LEN];
static limb aux[MAX_LEN], aux2[MAX_LEN];
static limb aux3[MAX_LEN], aux4[MAX_LEN];
static int NumberLength2;
int NumberLength, NumberLengthR1;
static limb U[MAX_LEN], V[MAX_LEN], R[MAX_LEN], S[MAX_LEN];
static BigInteger tmpDen, tmpNum;

int getNbrLimbs(limb *bigNbr)
{
  limb *ptrLimb = bigNbr + NumberLength;
  while (ptrLimb > bigNbr)
  {
    if ((--ptrLimb)->x != 0)
    {
      return (int)(ptrLimb - &bigNbr[0] + 1);
    }
  }
  return 1;
}

// Find the inverse of value mod 2^(NumberLength*BITS_PER_GROUP)
void ComputeInversePower2(limb *value, limb *result, limb *aux)
{
  int N, x, j;
  limb Cy;
  int currLen;
  x = N = (int)value->x;       // 2 least significant bits of inverse correct.
  x = x * (2 - N * x);         // 4 least significant bits of inverse correct.
  x = x * (2 - N * x);         // 8 least significant bits of inverse correct.
  x = x * (2 - N * x);         // 16 least significant bits of inverse correct.
  result->x = x & MAX_VALUE_LIMB;
  for (currLen = 2; currLen <= NumberLength * 2; currLen <<= 1)
  {
    multiply(value, result, aux, currLen, NULL);    // aux <- N * x
    Cy.x = 2 - aux[0].x;
    aux[0].x = Cy.x & MAX_VALUE_LIMB;
    for (j = 1; j < currLen; j++)
    {
      Cy.x = (Cy.x >> BITS_PER_GROUP) - aux[j].x;
      aux[j].x = Cy.x & MAX_VALUE_LIMB;
    }                                               // aux <- 2 - N * x
    multiply(result, aux, result, currLen, NULL);   // aux <- x * (2 - N * x)
  }
}

// Let R be a power of 2 of at least len limbs.
// Compute R1 = MontgomeryR1 and N = MontgomeryN using the formulas:
// R1 = R mod M
// N = -M^(-1) mod R

void GetMontgomeryParms(int len)
{
  int j;
  limb Cy;
  TestNbr[len].x = 0;
  NumberLength = len;
  NumberLength2 = len + len;
  if (NumberLength == 1)
  {
    MontgomeryMultR1[0].x = 1;
    MontgomeryMultR2[0].x = 1;
    NumberLengthR1 = 1;
    return;
  }
  // Compute MontgomeryMultN as -1/TestNbr (mod 2^k) using Newton method,
  // which doubles the precision for each iteration.
  // In the formula above: k = BITS_PER_GROUP * NumberLength.
  if (NumberLength >= 8)
  {
    limb *ptrResult;
    ComputeInversePower2(TestNbr, MontgomeryMultN, aux);
    ptrResult = &MontgomeryMultN[0];
    Cy.x = 0;          // Change sign.
    for (j = 0; j < NumberLength; j++)
    {
      Cy.x = (Cy.x >> BITS_PER_GROUP) - ptrResult->x;
      ptrResult->x = Cy.x & MAX_VALUE_LIMB;
      ptrResult++;
    }
    ptrResult->x = 0;
  }
  else
  {
    int x, N;
    x = N = (int)TestNbr[0].x;   // 2 least significant bits of inverse correct.
    x = x * (2 - N * x);         // 4 least significant bits of inverse correct.
    x = x * (2 - N * x);         // 8 least significant bits of inverse correct.
    x = x * (2 - N * x);         // 16 least significant bits of inverse correct.
    MontgomeryMultN[0].x = (-x) & MAX_VALUE_LIMB;    // Change sign
  }
  // Compute MontgomeryMultR1 as 1 in Montgomery notation.
  j = NumberLength;
  MontgomeryMultR1[j].x = 1;
  do
  {
    MontgomeryMultR1[--j].x = 0;
  } while (j > 0);
  AdjustModN(MontgomeryMultR1, TestNbr, len);
  MontgomeryMultR1[NumberLength].x = 0;
  memcpy(MontgomeryMultR2, MontgomeryMultR1, (NumberLength+1)*sizeof(limb));
  for (NumberLengthR1 = NumberLength; NumberLengthR1 > 0; NumberLengthR1--)
  {
    if (MontgomeryMultR1[NumberLengthR1 - 1].x != 0)
    {
      break;
    }
  }
  for (j = NumberLength; j > 0; j--)
  {
    memmove(&MontgomeryMultR2[1], &MontgomeryMultR2[0], NumberLength*sizeof(limb));
    MontgomeryMultR2[0].x = 0;
    AdjustModN(MontgomeryMultR2, TestNbr, len);
  }
}

// Compute Nbr <- Nbr mod TestNbr.
// TestNbr has NumberLength limbs.
void AdjustModN(limb Nbr[], limb Modulus[], int NumberLength)
{
  limb TrialQuotient, carry;
  int i;
  double dAux, dN;

  dN = (double) Modulus[NumberLength - 1].x;
  if (NumberLength > 1)
  {
    dN += (double) Modulus[NumberLength - 2].x / LIMB_RANGE;
  }
  if (NumberLength > 2)
  {
    dN += (double) Modulus[NumberLength - 3].x / (LIMB_RANGE*LIMB_RANGE);
  }
  dAux = (double) Nbr[NumberLength].x * LIMB_RANGE + (double) Nbr[NumberLength - 1].x;
  if (NumberLength > 1)
  {
    dAux += (double) Nbr[NumberLength - 2].x / LIMB_RANGE;
  }
#if BITS_PER_GROUP == 15
  TrialQuotient.x = (int)(dAux / dN) + 3;
#else
  TrialQuotient.x = (long long)(dAux / dN) + 3;
#endif
  if (TrialQuotient.x >= LIMB_RANGE)
  {
    carry.x = 0;
    for (i = 0; i < NumberLength; i++)
    {
      carry.x = Nbr[i + 1].x - (TrialQuotient.x >> BITS_PER_GROUP) * Modulus[i].x - carry.x;
      Nbr[i + 1].x = (int)(carry.x & MAX_VALUE_LIMB);
      carry.x = (MAX_VALUE_LIMB - carry.x) >> BITS_PER_GROUP;
    }
    TrialQuotient.x &= MAX_VALUE_LIMB;
  }
  carry.x = 0;
  for (i = 0; i < NumberLength; i++)
  {
    carry.x = Nbr[i].x - TrialQuotient.x * Modulus[i].x - carry.x;
    Nbr[i].x = (int)(carry.x & MAX_VALUE_LIMB);
    carry.x = (MAX_VALUE_LIMB - carry.x) >> BITS_PER_GROUP;
  }
  Nbr[NumberLength].x -= (int)carry.x;
  while ((Nbr[NumberLength].x & MAX_VALUE_LIMB) != 0)
  {
    carry.x = 0;
    for (i = 0; i < NumberLength; i++)
    {
      carry.x += Nbr[i].x + Modulus[i].x;
      Nbr[i].x = (int)(carry.x & MAX_VALUE_LIMB);
      carry.x >>= BITS_PER_GROUP;
    }
    Nbr[NumberLength].x += carry.x;
  }
}

void AddBigNbrModN(limb *Nbr1, limb *Nbr2, limb *Sum,
                   limb *TestNbr, int NumberLength)
{
  limb carry;
  int i;

  carry.x = 0;
  for (i = 0; i < NumberLength; i++)
  {
    carry.x = (carry.x >> BITS_PER_GROUP) + (Nbr1+i)->x + (Nbr2+i)->x - (TestNbr+i)->x;
    Sum[i].x = (int)(carry.x & MAX_VALUE_LIMB);
  }
  if (carry.x < 0)
  {
    carry.x = 0;
    for (i = 0; i < NumberLength; i++)
    {
      carry.x = (carry.x >> BITS_PER_GROUP) + (Sum+i)->x + (TestNbr+i)->x;
      Sum[i].x = (int)(carry.x & MAX_VALUE_LIMB);
    }
  }
}

void AddBigNbrMod(limb *Nbr1, limb *Nbr2, limb *Sum)
{
  AddBigNbrModN(Nbr1, Nbr2, Sum, TestNbr, NumberLength);
}

void SubtBigNbrModN(limb *Nbr1, limb *Nbr2, limb *Diff,
                    limb *TestNbr, int NumberLength)
{
  limb carry;
  int i;

  carry.x = 0;
  for (i = 0; i < NumberLength; i++)
  {
    carry.x = (carry.x >> BITS_PER_GROUP) + (Nbr1 + i)->x - (Nbr2 + i)->x;
    Diff[i].x = (int)(carry.x & MAX_VALUE_LIMB);
  }
  if (carry.x < 0)
  {
    carry.x = 0;
    for (i = 0; i < NumberLength; i++)
    {
      carry.x = (carry.x >> BITS_PER_GROUP) + (Diff + i)->x + (TestNbr + i)->x;
      Diff[i].x = (int)(carry.x & MAX_VALUE_LIMB);
    }
  }
}

void SubtBigNbrMod(limb *Nbr1, limb *Nbr2, limb *Sum)
{
  SubtBigNbrModN(Nbr1, Nbr2, Sum, TestNbr, NumberLength);
}

// Multiply two numbers in Montgomery notation.
//
// For large numbers the REDC algorithm is:
// m <- ((T mod R)N') mod R
// t <- (T + mN) / R
// if t = N then
//   return t - N
// else
//   return t
// end if

void modmult(limb *factor1, limb *factor2, limb *product)
{
  limb carry;
  int count;
  limb Pr, Nbr, MontDig;
  limb Prod[10];
  if (NumberLength == 1)
  {
    product->x = factor1->x * factor2->x % TestNbr[0].x;
    return;
  }
  if (NumberLength < 8)
  {     // Small numbers.
    int i, j;
    memset(Prod, 0, NumberLength*sizeof(limb));
    for (i = 0; i < NumberLength; i++)
    {
      Pr.x = (Nbr.x = (factor1 + i)->x) * factor2->x + Prod[0].x;
      MontDig.x = (Pr.x * MontgomeryMultN[0].x) & MAX_VALUE_LIMB;
      Prod[0].x = (Pr.x = ((MontDig.x * TestNbr[0].x + Pr.x) >> BITS_PER_GROUP) +
        MontDig.x * TestNbr[1].x + Nbr.x * (factor2+1)->x + Prod[1].x) & MAX_VALUE_LIMB;
      for (j = 2; j < NumberLength; j++)
      {
        Prod[j - 1].x = ((Pr.x = (Pr.x >> BITS_PER_GROUP) +
          MontDig.x * TestNbr[j].x + Nbr.x * (factor2 + j)->x + Prod[j].x) & MAX_VALUE_LIMB);
      }
      Prod[j - 1].x = Pr.x >> BITS_PER_GROUP;
    }
    for (j = NumberLength - 1; j >= 0; j--)
    {
      if (Prod[j].x != TestNbr[j].x)
      {
        break;
      }
    }
    if (j<0 || Prod[j].x >= TestNbr[j].x)
    {        // Prod >= TestNbr, so perform Prod <- Prod - TestNbr
      carry.x = 0;
      for (count = 0; count < NumberLength; count++)
      {
        carry.x += Prod[count].x - TestNbr[count].x;
        Prod[count].x = carry.x & MAX_VALUE_LIMB;
        carry.x >>= BITS_PER_GROUP;
      }
    }
    memcpy(product, Prod, NumberLength*sizeof(limb));
    return;
  }
  // Compute T
  multiply(factor1, factor2, product, NumberLength, NULL);
  // Compute m
  multiply(product, MontgomeryMultN, aux, NumberLength, NULL);
  // Compute mN
  multiply(TestNbr, aux, aux2, NumberLength, NULL);
  // Check if lowest half of mN is not zero
  carry.x = 0;
  for (count = 0; count < NumberLength; count++)
  {
    carry.x += aux2[count].x;
  }
  // If lowest half of mN is zero, compute hi(T) + hi(mN) - N
  // else compute hi(T) + hi(mN) + 1 - N
  // Where hi(number) is the high half of number.
  carry.x = (carry.x > 0 ? 1 : 0);
  for (count = 0; count < NumberLength; count++)
  {
    carry.x += (product+count+NumberLength)->x + aux2[count+NumberLength].x - TestNbr[count].x;
    (product+count)->x = carry.x & MAX_VALUE_LIMB;
    carry.x >>= BITS_PER_GROUP;
  }
  // If result is less than zero, add N to it.
  if (carry.x < 0)
  {
    carry.x = 0;
    for (count = 0; count < NumberLength; count++)
    {
      carry.x += (product+count)->x + TestNbr[count].x;
      (product+count)->x = carry.x & MAX_VALUE_LIMB;
      carry.x >>= BITS_PER_GROUP;
    }
  }
}

// Multiply big number in Montgomery notation by integer.
void modmultInt(limb *factorBig, int factorInt, limb *result)
{
  int index;
  limb carry;
  limb *ptrSrc = factorBig;
  limb *ptrDest = result;
  if (NumberLength == 1)
  {
    result->x = factorBig->x * factorInt % TestNbr[0].x;
    return;
  }
  carry.x = 0;
  for (index = 0; index < NumberLength; index++)
  {
    carry.x += factorInt * (ptrSrc++)->x;
    (ptrDest++)->x = carry.x & MAX_VALUE_LIMB;
    carry.x >>= BITS_PER_GROUP;
  }
  ptrDest->x = carry.x;
  AdjustModN(result, TestNbr, NumberLength);
}

// Input: base = base in Montgomery notation.
//        exp  = exponent.
//        nbrGroupsExp = number of limbs of exponent.
// Output: power = power in Montgomery notation.
void modPow(limb *base, limb *exp, int nbrGroupsExp, limb *power)
{
  int mask, index;
  memcpy(power, MontgomeryMultR1, (NumberLength + 1)*sizeof(*power));  // power <- 1
  for (index = nbrGroupsExp - 1; index >= 0; index--)
  {
    int groupExp = (int)(exp + index)->x;
    for (mask = 1 << (BITS_PER_GROUP - 1); mask > 0; mask >>= 1)
    {
      modmult(power, power, power);
      if ((groupExp & mask) != 0)
      {
        modmult(power, base, power);
      }
    }
  }
}

// Input: base = base in Montgomery notation.
//        exp  = exponent.
// Output: power = power in Montgomery notation.
void modPowLimb(limb *base, limb *exp, limb *power)
{
  int mask, groupExp;
  memcpy(power, MontgomeryMultR1, (NumberLength + 1)*sizeof(*power));  // power <- 1
  groupExp = (int)exp->x;
  for (mask = 1 << (BITS_PER_GROUP - 1); mask > 0; mask >>= 1)
  {
    modmult(power, power, power);
    if ((groupExp & mask) != 0)
    {
      modmult(power, base, power);
    }
  }
}

void modPowBaseInt(int base, limb *exp, int nbrGroupsExp, limb *power)
{
  int mask, index;
  memcpy(power, MontgomeryMultR1, (NumberLength+1)*sizeof(*power));  // power <- 1
  for (index = nbrGroupsExp-1; index>=0; index--)
  {
    int groupExp = (int)(exp+index)->x;
    for (mask = 1<<(BITS_PER_GROUP-1); mask > 0; mask >>= 1)
    {
      modmult(power, power, power);
      if ((groupExp & mask) != 0)
      {
        modmultInt(power, base, power);
      }
    }
  }
}

/* U' <- eU + fV, V' <- gU + hV                                        */
/* U <- U', V <- V'                                                    */
static void AddMult(limb *U, int e, int f, limb *V, int g, int h, int NumberLength)
{
  limb carryU, carryV;
  int i;
  carryU.x = carryV.x = 0;
  for (i = 0; i <= NumberLength; i++)
  {
    carryU.x = (carryU.x >> BITS_PER_GROUP) + U->x * e + V->x * f;
    carryV.x = (carryV.x >> BITS_PER_GROUP) + U->x * g + V->x * h;
    U->x = carryU.x & MAX_VALUE_LIMB;
    V->x = carryV.x & MAX_VALUE_LIMB;
    U++;
    V++;
  }
}

// Perform U <- (U - V) / 2
// U must be greater than V.
static int HalveDifference(limb *U, limb *V, int len)
{
  int i;
  limb carry;
    // Perform U <- U - V.
  carry.x = 0;
  for (i = 0; i < len; i++)
  {
    carry.x += U[i].x - V[i].x;
    U[i].x = carry.x & MAX_VALUE_LIMB;
    carry.x >>= BITS_PER_GROUP;
  }
    // Perform U <- U / 2.
  for (i = len - 1; i >= 0; i--)
  {
    carry.x += U[i].x;
    U[i].x = carry.x >> 1;
    carry.x = (carry.x & 1) << BITS_PER_GROUP;
  }
    // Get length of result.
  for (len--; len >= 0; len--)
  {
    if (U[len].x != 0)
    {
      break;
    }
  }
  return len+1;
}

int modInv(int NbrMod, int currentPrime)
{
  int QQ, T1, T3;
  int V1 = 1;
  int V3 = NbrMod;
  int U1 = 0;
  int U3 = currentPrime;
  while (V3 != 0)
  {
    if (U3 < V3 + V3)
    {               // QQ = 1
      T1 = U1 - V1;
      T3 = U3 - V3;
    }
    else
    {
      QQ = U3 / V3;
      T1 = U1 - V1 * QQ;
      T3 = U3 - V3 * QQ;
    }
    U1 = V1;
    U3 = V3;
    V1 = T1;
    V3 = T3;
  }
  return U1 + (currentPrime & (U1 >> 31));
}

/***********************************************************************/
/* NAME: ModInvBigNbr                                                  */
/*                                                                     */
/* PURPOSE: Find the inverse multiplicative modulo M.                  */
/* The algorithm terminates with inv = X^(-1) mod M.                   */
/*                                                                     */
/* This routine uses Kaliski Montgomery inverse algorithm              */
/* with changes by E. Savas and C. K. Koc.                             */
/*  1. U <- M, V <- X, R <- 0, S <- 1, k <- 0                          */
/*  2. while V > 0 do                                                  */
/*  3.   if U even then U <- U / 2, S <- 2S                            */
/*  4.   elsif V even then V <- V / 2, R <- 2R                         */
/*  5.   elsif U > V  then U <- (U - V) / 2, R <- R + S, S <- 2S       */
/*  6.   else V <- (V - U) / 2, S <- S + R, R <- 2R                    */
/*  7.   k <- k + 1                                                    */
/*  8. if R >= M then R <- R - M                                       */
/*  9. R <- M - R                                                      */
/* 10. R <- MonPro(R, R2)                                              */
/* 11. return MonPro(R, 2^(m-k))                                       */
/*                                                                     */
/*  In order to reduce the calculations, several single precision      */
/*  variables are added:                                               */
/*                                                                     */
/* R' <- aR + bS, S' <- cR + dS                                        */
/***********************************************************************/
void ModInvBigNbr(limb *num, limb *inv, limb *mod, int NumberLength)
{
  int len;
  int k;
  int a, b, c, d;
  int size, i;
  int bitCount = 0;
  int lenU, lenV;
  limb carry;
  if (NumberLength == 1)
  {
    inv->x = modInv(num->x, mod->x);
    return;
  }
  //  1. U <- M, V <- X, R <- 0, S <- 1, k <- 0
  size = (NumberLength+1)*sizeof(limb);
  memcpy(U, mod, size);
  memcpy(V, num, size);
  U[NumberLength].x = 0;
  V[NumberLength].x = 0;
    // Maximum value of R and S can be up to 2*M, so one more limb is needed.
  memset(R, 0, size);   // R <- 0
  memset(S, 0, size);   // S <- 1
  S[0].x = 1;
  k = 0;
  // R' <- aR + bS, S' <- cR + dS
  a = d = 1;  // R' = R, S' = S.
  b = c = 0;
  len = NumberLength;
  // Find length of U.
  for (lenU = NumberLength - 1; lenU > 0; lenU--)
  {
    if (U[lenU].x != 0)
    {
      break;
    }
  }
  lenU++;
  // Find length of V.
  for (lenV = NumberLength - 1; lenV > 0; lenV--)
  {
    if (V[lenV].x != 0)
    {
      break;
    }
  }
  lenV++;
  //  2. while V > 0 do
  while (lenV > 1 || V[0].x > 0)
  {
    carry.x = 0;
    //  3.   if U even then U <- U / 2, S <- 2S
    if ((U[0].x & 1) == 0)
    {     // U is even.
      for (i = lenU - 1; i >= 0; i--)
      {   // Loop that divides U by 2.
        carry.x += U[i].x;
        U[i].x = carry.x >> 1;
        carry.x = (carry.x & 1) << BITS_PER_GROUP;
      }
      if (U[lenU - 1].x == 0)
      {
        lenU--;
      }
      // R' <- aR + bS, S' <-cR + dS
      c *= 2; d *= 2;  // Multiply S by 2.
    }
    //  4.   elsif V even then V <- V / 2, R <- 2R
    else if ((V[0].x & 1) == 0)
    {    // V is even.
      for (i = lenV - 1; i >= 0; i--)
      {  // Loop that divides V by 2.
        carry.x += V[i].x;
        V[i].x = carry.x >> 1;
        carry.x = (carry.x & 1) << BITS_PER_GROUP;
      }
      if (V[lenV - 1].x == 0)
      {
        lenV--;
      }
      // R' <- aR + bS, S' <- cR + dS
      a *= 2; b *= 2;  // Multiply R by 2.
    }
    //  5.   elsif U >= V  then U <- (U - V) / 2, R <- R + S, S <- 2S
    else
    {
      len = (lenU > lenV ? lenU : lenV);
      for (i = len - 1; i > 0; i--)
      {
        if (U[i].x != V[i].x)
        {
          break;
        }
      }
      if (U[i].x > V[i].x)
      {
        lenU = HalveDifference(U, V, len);
        // R' <- aR + bS, S' <- cR + dS
        a += c; b += d;  // R <- R + S
        c *= 2; d *= 2;  // S <- 2S
      }
      //  6.   elsif V > U then V <- (V - U) / 2, S <- S + R, R <- 2R
      else
      {
        lenV = HalveDifference(V, U, len);
        // R' <- aR + bS, S' <- cR + dS
        c += a; d += b;  // S <- S + R
        a *= 2; b *= 2;  // R <- 2R
      }
    }
    //  7.   k <- k + 1
    k++;
    // Adjust variables.
    if (++bitCount == BITS_PER_GROUP-1)
    {
      AddMult(R, a, b, S, c, d, NumberLength+1);
      // R' <- aR + bS, S' <- cR + dS
      a = d = 1;  // R' = R, S' = S.
      b = c = 0;
      bitCount = 0;
    }
  }
  AddMult(R, a, b, S, c, d, NumberLength+1);
  //  8. if R >= M then R <- R - M
  carry.x = 0;
  for (i = 0; i <= NumberLength; i++)
  {
    carry.x += R[i].x - (mod + i)->x;
    R[i].x = carry.x & MAX_VALUE_LIMB;
    carry.x >>= BITS_PER_GROUP;
  }
  if (carry.x < 0)
  {
    carry.x = 0;
    for (i = 0; i <= NumberLength; i++)
    {
      carry.x += R[i].x + (mod + i)->x;
      R[i].x = carry.x & MAX_VALUE_LIMB;
      carry.x >>= BITS_PER_GROUP;
    }
  }
  //  9. R <- M - R
  carry.x = 0;
  for (i = 0; i <= NumberLength; i++)
  {
    carry.x += (mod + i)->x - R[i].x;
    R[i].x = carry.x & MAX_VALUE_LIMB;
    carry.x >>= BITS_PER_GROUP;
  }
  // At this moment R = x^(-1)*2^k
  // 10. R <- MonPro(R, R2)
  modmult(R, MontgomeryMultR2, R);
  // At this moment R = x^(-1)*2^(k+m)
  // 11. return MonPro(R, 2^(m-k))
  memset(S, 0, size);
  bitCount = NumberLength*BITS_PER_GROUP - k;
  if (bitCount < 0)
  {
    bitCount += NumberLength*BITS_PER_GROUP;
    S[bitCount / BITS_PER_GROUP].x = 1 << (bitCount % BITS_PER_GROUP);
    modmult(R, S, inv);
  }
  else
  {
    S[bitCount / BITS_PER_GROUP].x = 1 << (bitCount % BITS_PER_GROUP);
    modmult(R, S, inv);
    modmult(inv, MontgomeryMultR2, inv);
  }
}

void BigIntModularDivision(BigInteger *Num, BigInteger *Den, BigInteger *mod, BigInteger *quotient)
{
  NumberLength = mod -> nbrLimbs;
  memcpy(TestNbr, mod -> limbs, NumberLength * sizeof(limb));
  TestNbr[NumberLength].x = 0;
  GetMontgomeryParms(NumberLength);
  BigIntRemainder(Num, mod, &tmpNum);
  BigIntRemainder(Den, mod, &tmpDen);
  if (tmpNum.sign == SIGN_NEGATIVE)
  {
    BigIntAdd(&tmpNum, mod, &tmpNum);
  }
  if (tmpDen.sign == SIGN_NEGATIVE)
  {
    BigIntAdd(&tmpDen, mod, &tmpDen);
  }
  CompressLimbsBigInteger(aux3, &tmpDen);
  modmult(aux3, MontgomeryMultR2, aux3);  // aux3 <- Den in Montgomery notation
  ModInvBigNbr(aux3, aux3, TestNbr, NumberLength); // aux3 <- 1 / Den in Montg notation.
  CompressLimbsBigInteger(aux4, &tmpNum);
  modmult(aux3, aux4, aux3);              // aux3 <- Num / Dev in standard notation.
  UncompressLimbsBigInteger(aux3, quotient);  // Get Num/Den
}

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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include "bignbr.h"

// Multiply two numbers in Montgomery notation.
//
// For large numbers the REDC algorithm is:
// m <- ((T mod R)N') mod R
// t <- (T - mN) / R
// if t < 0 then
//   return t + N
// else
//   return t
// end if

#define MONTGOMERY_MULT_THRESHOLD 13
#define MONTMULT_LIMB_START   \
    int32_t Nbr = (pNbr1 + i)->x;   \
    int64_t Pr = ((int64_t)Nbr * (int64_t)Nbr2_0) + (int64_t)Prod0; \
    int32_t MontDig = ((int32_t)Pr * (int32_t)MontgomeryMultN[0].x) & MAX_INT_NBR_U; \
    Pr -= (int64_t)MontDig * (int64_t)TestNbr0
#define MONTMULT_LIMB(curr, next)  \
    Pr = ((int64_t)Nbr * (int64_t)Nbr2_##next) + (int64_t)Prod##next + \
         (Pr >> BITS_PER_GROUP) -  \
         ((int64_t)MontDig * (int64_t)TestNbr##next); \
    Prod##curr = (int32_t)Pr & MAX_INT_NBR_U
#define MONTMULT_LIMB_END(curr)   \
    Prod##curr = (int32_t)(Pr >> BITS_PER_GROUP)

#ifdef _USING64BITS_
static void MontgomeryMult2(const limb* pNbr1, const limb* pNbr2, limb* pProd)
{
  int32_t Prod0 = 0;
  int32_t Prod1 = 0;
  uint32_t TestNbr0 = TestNbr[0].x;
  uint32_t TestNbr1 = TestNbr[1].x;
  uint32_t Nbr2_0 = pNbr2->x;
  uint32_t Nbr2_1 = (pNbr2 + 1)->x;
  for (int i = 0; i < 2; i++)
  {
    MONTMULT_LIMB_START;
    MONTMULT_LIMB(0, 1);
    MONTMULT_LIMB_END(1);
  }
  if (Prod1 < 0)
  {
    uint32_t carry = Prod0 + TestNbr0;
    Prod0 = carry & MAX_INT_NBR_U;
    Prod1 += TestNbr1 + (carry >> BITS_PER_GROUP);
  }
  pProd->x = Prod0;
  (pProd + 1)->x = Prod1;
}

static void MontgomeryMult3(const limb* pNbr1, const limb* pNbr2, limb* pProd)
{
  int32_t Prod0 = 0;
  int32_t Prod1 = 0;
  int32_t Prod2 = 0;
  uint32_t TestNbr0 = TestNbr[0].x;
  uint32_t TestNbr1 = TestNbr[1].x;
  uint32_t TestNbr2 = TestNbr[2].x;
  uint32_t Nbr2_0 = pNbr2->x;
  uint32_t Nbr2_1 = (pNbr2 + 1)->x;
  uint32_t Nbr2_2 = (pNbr2 + 2)->x;
  for (int i = 0; i < 3; i++)
  {
    MONTMULT_LIMB_START;
    MONTMULT_LIMB(0, 1);
    MONTMULT_LIMB(1, 2);
    MONTMULT_LIMB_END(2);
  }
  if (Prod2 < 0)
  {
    uint32_t carry = Prod0 + TestNbr0;
    Prod0 = carry & MAX_INT_NBR_U;
    carry = Prod1 + TestNbr1 + (carry >> BITS_PER_GROUP);
    Prod1 = carry & MAX_INT_NBR_U;
    Prod2 += TestNbr2 + (carry >> BITS_PER_GROUP);
  }
  pProd->x = Prod0;
  (pProd + 1)->x = Prod1;
  (pProd + 2)->x = Prod2;
}

static void MontgomeryMult4(const limb* pNbr1, const limb* pNbr2, limb* pProd)
{
  int32_t Prod0 = 0;
  int32_t Prod1 = 0;
  int32_t Prod2 = 0;
  int32_t Prod3 = 0;
  uint32_t TestNbr0 = TestNbr[0].x;
  uint32_t TestNbr1 = TestNbr[1].x;
  uint32_t TestNbr2 = TestNbr[2].x;
  uint32_t TestNbr3 = TestNbr[3].x;
  uint32_t Nbr2_0 = pNbr2->x;
  uint32_t Nbr2_1 = (pNbr2 + 1)->x;
  uint32_t Nbr2_2 = (pNbr2 + 2)->x;
  uint32_t Nbr2_3 = (pNbr2 + 3)->x;
  for (int i = 0; i < 4; i++)
  {
    MONTMULT_LIMB_START;
    MONTMULT_LIMB(0, 1);
    MONTMULT_LIMB(1, 2);
    MONTMULT_LIMB(2, 3);
    MONTMULT_LIMB_END(3);
  }
  if (Prod3 < 0)
  {
    uint32_t carry = Prod0 + TestNbr0;
    Prod0 = carry & MAX_INT_NBR_U;
    carry = Prod1 + TestNbr1 + (carry >> BITS_PER_GROUP);
    Prod1 = carry & MAX_INT_NBR_U;
    carry = Prod2 + TestNbr2 + (carry >> BITS_PER_GROUP);
    Prod2 = carry & MAX_INT_NBR_U;
    Prod3 += TestNbr3 + (carry >> BITS_PER_GROUP);
  }
  pProd->x = Prod0;
  (pProd + 1)->x = Prod1;
  (pProd + 2)->x = Prod2;
  (pProd + 3)->x = Prod3;
}

static void MontgomeryMult5(const limb* pNbr1, const limb* pNbr2, limb* pProd)
{
  int32_t Prod0 = 0;
  int32_t Prod1 = 0;
  int32_t Prod2 = 0;
  int32_t Prod3 = 0;
  int32_t Prod4 = 0;
  uint32_t TestNbr0 = TestNbr[0].x;
  uint32_t TestNbr1 = TestNbr[1].x;
  uint32_t TestNbr2 = TestNbr[2].x;
  uint32_t TestNbr3 = TestNbr[3].x;
  uint32_t TestNbr4 = TestNbr[4].x;
  uint32_t Nbr2_0 = pNbr2->x;
  uint32_t Nbr2_1 = (pNbr2 + 1)->x;
  uint32_t Nbr2_2 = (pNbr2 + 2)->x;
  uint32_t Nbr2_3 = (pNbr2 + 3)->x;
  uint32_t Nbr2_4 = (pNbr2 + 4)->x;
  for (int i = 0; i < 5; i++)
  {
    MONTMULT_LIMB_START;
    MONTMULT_LIMB(0, 1);
    MONTMULT_LIMB(1, 2);
    MONTMULT_LIMB(2, 3);
    MONTMULT_LIMB(3, 4);
    MONTMULT_LIMB_END(4);
  }
  if (Prod4 < 0)
  {
    uint32_t carry = Prod0 + TestNbr0;
    Prod0 = carry & MAX_INT_NBR_U;
    carry = Prod1 + TestNbr1 + (carry >> BITS_PER_GROUP);
    Prod1 = carry & MAX_INT_NBR_U;
    carry = Prod2 + TestNbr2 + (carry >> BITS_PER_GROUP);
    Prod2 = carry & MAX_INT_NBR_U;
    carry = Prod3 + TestNbr3 + (carry >> BITS_PER_GROUP);
    Prod3 = carry & MAX_INT_NBR_U;
    Prod4 += TestNbr4 + (carry >> BITS_PER_GROUP);
  }
  pProd->x = Prod0;
  (pProd + 1)->x = Prod1;
  (pProd + 2)->x = Prod2;
  (pProd + 3)->x = Prod3;
  (pProd + 4)->x = Prod4;
}

static void MontgomeryMult6(const limb* pNbr1, const limb* pNbr2, limb* pProd)
{
  int32_t Prod0 = 0;
  int32_t Prod1 = 0;
  int32_t Prod2 = 0;
  int32_t Prod3 = 0;
  int32_t Prod4 = 0;
  int32_t Prod5 = 0;
  uint32_t TestNbr0 = TestNbr[0].x;
  uint32_t TestNbr1 = TestNbr[1].x;
  uint32_t TestNbr2 = TestNbr[2].x;
  uint32_t TestNbr3 = TestNbr[3].x;
  uint32_t TestNbr4 = TestNbr[4].x;
  uint32_t TestNbr5 = TestNbr[5].x;
  uint32_t Nbr2_0 = pNbr2->x;
  uint32_t Nbr2_1 = (pNbr2 + 1)->x;
  uint32_t Nbr2_2 = (pNbr2 + 2)->x;
  uint32_t Nbr2_3 = (pNbr2 + 3)->x;
  uint32_t Nbr2_4 = (pNbr2 + 4)->x;
  uint32_t Nbr2_5 = (pNbr2 + 5)->x;
  for (int i = 0; i < 6; i++)
  {
    MONTMULT_LIMB_START;
    MONTMULT_LIMB(0, 1);
    MONTMULT_LIMB(1, 2);
    MONTMULT_LIMB(2, 3);
    MONTMULT_LIMB(3, 4);
    MONTMULT_LIMB(4, 5);
    MONTMULT_LIMB_END(5);
  }
  if (Prod5 < 0)
  {
    uint32_t carry = Prod0 + TestNbr0;
    Prod0 = carry & MAX_INT_NBR_U;
    carry = Prod1 + TestNbr1 + (carry >> BITS_PER_GROUP);
    Prod1 = carry & MAX_INT_NBR_U;
    carry = Prod2 + TestNbr2 + (carry >> BITS_PER_GROUP);
    Prod2 = carry & MAX_INT_NBR_U;
    carry = Prod3 + TestNbr3 + (carry >> BITS_PER_GROUP);
    Prod3 = carry & MAX_INT_NBR_U;
    carry = Prod4 + TestNbr4 + (carry >> BITS_PER_GROUP);
    Prod4 = carry & MAX_INT_NBR_U;
    Prod5 += TestNbr5 + (carry >> BITS_PER_GROUP);
  }
  pProd->x = Prod0;
  (pProd + 1)->x = Prod1;
  (pProd + 2)->x = Prod2;
  (pProd + 3)->x = Prod3;
  (pProd + 4)->x = Prod4;
  (pProd + 5)->x = Prod5;
}

static void MontgomeryMult7(const limb* pNbr1, const limb* pNbr2, limb* pProd)
{
  int32_t Prod0 = 0;
  int32_t Prod1 = 0;
  int32_t Prod2 = 0;
  int32_t Prod3 = 0;
  int32_t Prod4 = 0;
  int32_t Prod5 = 0;
  int32_t Prod6 = 0;
  uint32_t TestNbr0 = TestNbr[0].x;
  uint32_t TestNbr1 = TestNbr[1].x;
  uint32_t TestNbr2 = TestNbr[2].x;
  uint32_t TestNbr3 = TestNbr[3].x;
  uint32_t TestNbr4 = TestNbr[4].x;
  uint32_t TestNbr5 = TestNbr[5].x;
  uint32_t TestNbr6 = TestNbr[6].x;
  int Nbr2_0 = pNbr2->x;
  int Nbr2_1 = (pNbr2 + 1)->x;
  int Nbr2_2 = (pNbr2 + 2)->x;
  int Nbr2_3 = (pNbr2 + 3)->x;
  int Nbr2_4 = (pNbr2 + 4)->x;
  int Nbr2_5 = (pNbr2 + 5)->x;
  int Nbr2_6 = (pNbr2 + 6)->x;
  for (int i = 0; i < 7; i++)
  {
    MONTMULT_LIMB_START;
    MONTMULT_LIMB(0, 1);
    MONTMULT_LIMB(1, 2);
    MONTMULT_LIMB(2, 3);
    MONTMULT_LIMB(3, 4);
    MONTMULT_LIMB(4, 5);
    MONTMULT_LIMB(5, 6);
    MONTMULT_LIMB_END(6);
  }
  if (Prod6 < 0)
  {
    uint32_t carry = Prod0 + TestNbr0;
    Prod0 = carry & MAX_INT_NBR_U;
    carry = Prod1 + TestNbr1 + (carry >> BITS_PER_GROUP);
    Prod1 = carry & MAX_INT_NBR_U;
    carry = Prod2 + TestNbr2 + (carry >> BITS_PER_GROUP);
    Prod2 = carry & MAX_INT_NBR_U;
    carry = Prod3 + TestNbr3 + (carry >> BITS_PER_GROUP);
    Prod3 = carry & MAX_INT_NBR_U;
    carry = Prod4 + TestNbr4 + (carry >> BITS_PER_GROUP);
    Prod4 = carry & MAX_INT_NBR_U;
    carry = Prod5 + TestNbr5 + (carry >> BITS_PER_GROUP);
    Prod5 = carry & MAX_INT_NBR_U;
    Prod6 += TestNbr6 + (carry >> BITS_PER_GROUP);
  }
  pProd->x = Prod0;
  (pProd + 1)->x = Prod1;
  (pProd + 2)->x = Prod2;
  (pProd + 3)->x = Prod3;
  (pProd + 4)->x = Prod4;
  (pProd + 5)->x = Prod5;
  (pProd + 6)->x = Prod6;
}

static void MontgomeryMult8(const limb* pNbr1, const limb* pNbr2, limb* pProd)
{
  int32_t Prod0 = 0;
  int32_t Prod1 = 0;
  int32_t Prod2 = 0;
  int32_t Prod3 = 0;
  int32_t Prod4 = 0;
  int32_t Prod5 = 0;
  int32_t Prod6 = 0;
  int32_t Prod7 = 0;
  uint32_t TestNbr0 = TestNbr[0].x;
  uint32_t TestNbr1 = TestNbr[1].x;
  uint32_t TestNbr2 = TestNbr[2].x;
  uint32_t TestNbr3 = TestNbr[3].x;
  uint32_t TestNbr4 = TestNbr[4].x;
  uint32_t TestNbr5 = TestNbr[5].x;
  uint32_t TestNbr6 = TestNbr[6].x;
  uint32_t TestNbr7 = TestNbr[7].x;
  uint32_t Nbr2_0 = pNbr2->x;
  uint32_t Nbr2_1 = (pNbr2 + 1)->x;
  uint32_t Nbr2_2 = (pNbr2 + 2)->x;
  uint32_t Nbr2_3 = (pNbr2 + 3)->x;
  uint32_t Nbr2_4 = (pNbr2 + 4)->x;
  uint32_t Nbr2_5 = (pNbr2 + 5)->x;
  uint32_t Nbr2_6 = (pNbr2 + 6)->x;
  uint32_t Nbr2_7 = (pNbr2 + 7)->x;
  for (int i = 0; i < 8; i++)
  {
    MONTMULT_LIMB_START;
    MONTMULT_LIMB(0, 1);
    MONTMULT_LIMB(1, 2);
    MONTMULT_LIMB(2, 3);
    MONTMULT_LIMB(3, 4);
    MONTMULT_LIMB(4, 5);
    MONTMULT_LIMB(5, 6);
    MONTMULT_LIMB(6, 7);
    MONTMULT_LIMB_END(7);
  }
  if (Prod7 < 0)
  {
    uint32_t carry = Prod0 + TestNbr0;
    Prod0 = carry & MAX_INT_NBR_U;
    carry = Prod1 + TestNbr1 + (carry >> BITS_PER_GROUP);
    Prod1 = carry & MAX_INT_NBR_U;
    carry = Prod2 + TestNbr2 + (carry >> BITS_PER_GROUP);
    Prod2 = carry & MAX_INT_NBR_U;
    carry = Prod3 + TestNbr3 + (carry >> BITS_PER_GROUP);
    Prod3 = carry & MAX_INT_NBR_U;
    carry = Prod4 + TestNbr4 + (carry >> BITS_PER_GROUP);
    Prod4 = carry & MAX_INT_NBR_U;
    carry = Prod5 + TestNbr5 + (carry >> BITS_PER_GROUP);
    Prod5 = carry & MAX_INT_NBR_U;
    carry = Prod6 + TestNbr6 + (carry >> BITS_PER_GROUP);
    Prod6 = carry & MAX_INT_NBR_U;
    Prod7 += TestNbr7 + (carry >> BITS_PER_GROUP);
  }
  pProd->x = Prod0;
  (pProd + 1)->x = Prod1;
  (pProd + 2)->x = Prod2;
  (pProd + 3)->x = Prod3;
  (pProd + 4)->x = Prod4;
  (pProd + 5)->x = Prod5;
  (pProd + 6)->x = Prod6;
  (pProd + 7)->x = Prod7;
}
static void MontgomeryMult9(const limb* pNbr1, const limb* pNbr2, limb* pProd)
{
  int32_t Prod0 = 0;
  int32_t Prod1 = 0;
  int32_t Prod2 = 0;
  int32_t Prod3 = 0;
  int32_t Prod4 = 0;
  int32_t Prod5 = 0;
  int32_t Prod6 = 0;
  int32_t Prod7 = 0;
  int32_t Prod8 = 0;
  uint32_t TestNbr0 = TestNbr[0].x;
  uint32_t TestNbr1 = TestNbr[1].x;
  uint32_t TestNbr2 = TestNbr[2].x;
  uint32_t TestNbr3 = TestNbr[3].x;
  uint32_t TestNbr4 = TestNbr[4].x;
  uint32_t TestNbr5 = TestNbr[5].x;
  uint32_t TestNbr6 = TestNbr[6].x;
  uint32_t TestNbr7 = TestNbr[7].x;
  uint32_t TestNbr8 = TestNbr[8].x;
  uint32_t Nbr2_0 = pNbr2->x;
  uint32_t Nbr2_1 = (pNbr2 + 1)->x;
  uint32_t Nbr2_2 = (pNbr2 + 2)->x;
  uint32_t Nbr2_3 = (pNbr2 + 3)->x;
  uint32_t Nbr2_4 = (pNbr2 + 4)->x;
  uint32_t Nbr2_5 = (pNbr2 + 5)->x;
  uint32_t Nbr2_6 = (pNbr2 + 6)->x;
  uint32_t Nbr2_7 = (pNbr2 + 7)->x;
  uint32_t Nbr2_8 = (pNbr2 + 8)->x;
  for (int i = 0; i < 9; i++)
  {
    MONTMULT_LIMB_START;
    MONTMULT_LIMB(0, 1);
    MONTMULT_LIMB(1, 2);
    MONTMULT_LIMB(2, 3);
    MONTMULT_LIMB(3, 4);
    MONTMULT_LIMB(4, 5);
    MONTMULT_LIMB(5, 6);
    MONTMULT_LIMB(6, 7);
    MONTMULT_LIMB(7, 8);
    MONTMULT_LIMB_END(8);
  }
  if (Prod8 < 0)
  {
    uint32_t carry = Prod0 + TestNbr0;
    Prod0 = carry & MAX_INT_NBR_U;
    carry = Prod1 + TestNbr1 + (carry >> BITS_PER_GROUP);
    Prod1 = carry & MAX_INT_NBR_U;
    carry = Prod2 + TestNbr2 + (carry >> BITS_PER_GROUP);
    Prod2 = carry & MAX_INT_NBR_U;
    carry = Prod3 + TestNbr3 + (carry >> BITS_PER_GROUP);
    Prod3 = carry & MAX_INT_NBR_U;
    carry = Prod4 + TestNbr4 + (carry >> BITS_PER_GROUP);
    Prod4 = carry & MAX_INT_NBR_U;
    carry = Prod5 + TestNbr5 + (carry >> BITS_PER_GROUP);
    Prod5 = carry & MAX_INT_NBR_U;
    carry = Prod6 + TestNbr6 + (carry >> BITS_PER_GROUP);
    Prod6 = carry & MAX_INT_NBR_U;
    carry = Prod7 + TestNbr7 + (carry >> BITS_PER_GROUP);
    Prod7 = carry & MAX_INT_NBR_U;
    Prod8 += TestNbr8 + (carry >> BITS_PER_GROUP);
  }
  pProd->x = Prod0;
  (pProd + 1)->x = Prod1;
  (pProd + 2)->x = Prod2;
  (pProd + 3)->x = Prod3;
  (pProd + 4)->x = Prod4;
  (pProd + 5)->x = Prod5;
  (pProd + 6)->x = Prod6;
  (pProd + 7)->x = Prod7;
  (pProd + 8)->x = Prod8;
}

static void MontgomeryMult10(const limb* pNbr1, const limb* pNbr2, limb* pProd)
{
  int32_t Prod0 = 0;
  int32_t Prod1 = 0;
  int32_t Prod2 = 0;
  int32_t Prod3 = 0;
  int32_t Prod4 = 0;
  int32_t Prod5 = 0;
  int32_t Prod6 = 0;
  int32_t Prod7 = 0;
  int32_t Prod8 = 0;
  int32_t Prod9 = 0;
  uint32_t TestNbr0 = TestNbr[0].x;
  uint32_t TestNbr1 = TestNbr[1].x;
  uint32_t TestNbr2 = TestNbr[2].x;
  uint32_t TestNbr3 = TestNbr[3].x;
  uint32_t TestNbr4 = TestNbr[4].x;
  uint32_t TestNbr5 = TestNbr[5].x;
  uint32_t TestNbr6 = TestNbr[6].x;
  uint32_t TestNbr7 = TestNbr[7].x;
  uint32_t TestNbr8 = TestNbr[8].x;
  uint32_t TestNbr9 = TestNbr[9].x;
  uint32_t Nbr2_0 = pNbr2->x;
  uint32_t Nbr2_1 = (pNbr2 + 1)->x;
  uint32_t Nbr2_2 = (pNbr2 + 2)->x;
  uint32_t Nbr2_3 = (pNbr2 + 3)->x;
  uint32_t Nbr2_4 = (pNbr2 + 4)->x;
  uint32_t Nbr2_5 = (pNbr2 + 5)->x;
  uint32_t Nbr2_6 = (pNbr2 + 6)->x;
  uint32_t Nbr2_7 = (pNbr2 + 7)->x;
  uint32_t Nbr2_8 = (pNbr2 + 8)->x;
  uint32_t Nbr2_9 = (pNbr2 + 9)->x;
  for (int i = 0; i < 10; i++)
  {
    MONTMULT_LIMB_START;
    MONTMULT_LIMB(0, 1);
    MONTMULT_LIMB(1, 2);
    MONTMULT_LIMB(2, 3);
    MONTMULT_LIMB(3, 4);
    MONTMULT_LIMB(4, 5);
    MONTMULT_LIMB(5, 6);
    MONTMULT_LIMB(6, 7);
    MONTMULT_LIMB(7, 8);
    MONTMULT_LIMB(8, 9);
    MONTMULT_LIMB_END(9);
  }
  if (Prod9 < 0)
  {
    uint32_t carry = Prod0 + TestNbr0;
    Prod0 = carry & MAX_INT_NBR_U;
    carry = Prod1 + TestNbr1 + (carry >> BITS_PER_GROUP);
    Prod1 = carry & MAX_INT_NBR_U;
    carry = Prod2 + TestNbr2 + (carry >> BITS_PER_GROUP);
    Prod2 = carry & MAX_INT_NBR_U;
    carry = Prod3 + TestNbr3 + (carry >> BITS_PER_GROUP);
    Prod3 = carry & MAX_INT_NBR_U;
    carry = Prod4 + TestNbr4 + (carry >> BITS_PER_GROUP);
    Prod4 = carry & MAX_INT_NBR_U;
    carry = Prod5 + TestNbr5 + (carry >> BITS_PER_GROUP);
    Prod5 = carry & MAX_INT_NBR_U;
    carry = Prod6 + TestNbr6 + (carry >> BITS_PER_GROUP);
    Prod6 = carry & MAX_INT_NBR_U;
    carry = Prod7 + TestNbr7 + (carry >> BITS_PER_GROUP);
    Prod7 = carry & MAX_INT_NBR_U;
    carry = Prod8 + TestNbr8 + (carry >> BITS_PER_GROUP);
    Prod8 = carry & MAX_INT_NBR_U;
    Prod9 += TestNbr9 + (carry >> BITS_PER_GROUP);
  }
  pProd->x = Prod0;
  (pProd + 1)->x = Prod1;
  (pProd + 2)->x = Prod2;
  (pProd + 3)->x = Prod3;
  (pProd + 4)->x = Prod4;
  (pProd + 5)->x = Prod5;
  (pProd + 6)->x = Prod6;
  (pProd + 7)->x = Prod7;
  (pProd + 8)->x = Prod8;
  (pProd + 9)->x = Prod9;
}

static void MontgomeryMult11(const limb* pNbr1, const limb* pNbr2, limb* pProd)
{
  int32_t Prod0 = 0;
  int32_t Prod1 = 0;
  int32_t Prod2 = 0;
  int32_t Prod3 = 0;
  int32_t Prod4 = 0;
  int32_t Prod5 = 0;
  int32_t Prod6 = 0;
  int32_t Prod7 = 0;
  int32_t Prod8 = 0;
  int32_t Prod9 = 0;
  int32_t Prod10 = 0;
  uint32_t TestNbr0 = TestNbr[0].x;
  uint32_t TestNbr1 = TestNbr[1].x;
  uint32_t TestNbr2 = TestNbr[2].x;
  uint32_t TestNbr3 = TestNbr[3].x;
  uint32_t TestNbr4 = TestNbr[4].x;
  uint32_t TestNbr5 = TestNbr[5].x;
  uint32_t TestNbr6 = TestNbr[6].x;
  uint32_t TestNbr7 = TestNbr[7].x;
  uint32_t TestNbr8 = TestNbr[8].x;
  uint32_t TestNbr9 = TestNbr[9].x;
  uint32_t TestNbr10 = TestNbr[10].x;
  uint32_t Nbr2_0 = pNbr2->x;
  uint32_t Nbr2_1 = (pNbr2 + 1)->x;
  uint32_t Nbr2_2 = (pNbr2 + 2)->x;
  uint32_t Nbr2_3 = (pNbr2 + 3)->x;
  uint32_t Nbr2_4 = (pNbr2 + 4)->x;
  uint32_t Nbr2_5 = (pNbr2 + 5)->x;
  uint32_t Nbr2_6 = (pNbr2 + 6)->x;
  uint32_t Nbr2_7 = (pNbr2 + 7)->x;
  uint32_t Nbr2_8 = (pNbr2 + 8)->x;
  uint32_t Nbr2_9 = (pNbr2 + 9)->x;
  uint32_t Nbr2_10 = (pNbr2 + 10)->x;
  for (int i = 0; i < 11; i++)
  {
    MONTMULT_LIMB_START;
    MONTMULT_LIMB(0, 1);
    MONTMULT_LIMB(1, 2);
    MONTMULT_LIMB(2, 3);
    MONTMULT_LIMB(3, 4);
    MONTMULT_LIMB(4, 5);
    MONTMULT_LIMB(5, 6);
    MONTMULT_LIMB(6, 7);
    MONTMULT_LIMB(7, 8);
    MONTMULT_LIMB(8, 9);
    MONTMULT_LIMB(9, 10);
    MONTMULT_LIMB_END(10);
  }
  if (Prod10 < 0)
  {
    uint32_t carry = Prod0 + TestNbr0;
    Prod0 = carry & MAX_INT_NBR_U;
    carry = Prod1 + TestNbr1 + (carry >> BITS_PER_GROUP);
    Prod1 = carry & MAX_INT_NBR_U;
    carry = Prod2 + TestNbr2 + (carry >> BITS_PER_GROUP);
    Prod2 = carry & MAX_INT_NBR_U;
    carry = Prod3 + TestNbr3 + (carry >> BITS_PER_GROUP);
    Prod3 = carry & MAX_INT_NBR_U;
    carry = Prod4 + TestNbr4 + (carry >> BITS_PER_GROUP);
    Prod4 = carry & MAX_INT_NBR_U;
    carry = Prod5 + TestNbr5 + (carry >> BITS_PER_GROUP);
    Prod5 = carry & MAX_INT_NBR_U;
    carry = Prod6 + TestNbr6 + (carry >> BITS_PER_GROUP);
    Prod6 = carry & MAX_INT_NBR_U;
    carry = Prod7 + TestNbr7 + (carry >> BITS_PER_GROUP);
    Prod7 = carry & MAX_INT_NBR_U;
    carry = Prod8 + TestNbr8 + (carry >> BITS_PER_GROUP);
    Prod8 = carry & MAX_INT_NBR_U;
    carry = Prod9 + TestNbr9 + (carry >> BITS_PER_GROUP);
    Prod9 = carry & MAX_INT_NBR_U;
    Prod10 += TestNbr10 + (carry >> BITS_PER_GROUP);
  }
  pProd->x = Prod0;
  (pProd + 1)->x = Prod1;
  (pProd + 2)->x = Prod2;
  (pProd + 3)->x = Prod3;
  (pProd + 4)->x = Prod4;
  (pProd + 5)->x = Prod5;
  (pProd + 6)->x = Prod6;
  (pProd + 7)->x = Prod7;
  (pProd + 8)->x = Prod8;
  (pProd + 9)->x = Prod9;
  (pProd + 10)->x = Prod10;
}
#endif

void MontgomeryMult(const limb* factor1, const limb* factor2, limb* product)
{
  assert(NumberLength > 1);
  int32_t Prod[MONTGOMERY_MULT_THRESHOLD + 1];
  int NumberLengthBytes;

#ifdef _USING64BITS_
  switch (NumberLength)
#endif  
  {
#ifdef _USING64BITS_
  case 2:
    MontgomeryMult2(factor1, factor2, product);
    return;
  case 3:
    MontgomeryMult3(factor1, factor2, product);
    return;
  case 4:
    MontgomeryMult4(factor1, factor2, product);
    return;
  case 5:
    MontgomeryMult5(factor1, factor2, product);
    return;
  case 6:
    MontgomeryMult6(factor1, factor2, product);
    return;
  case 7:
    MontgomeryMult7(factor1, factor2, product);
    return;
  case 8:
    MontgomeryMult8(factor1, factor2, product);
    return;
  case 9:
    MontgomeryMult9(factor1, factor2, product);
    return;
  case 10:
    MontgomeryMult10(factor1, factor2, product);
    return;
  case 11:
    MontgomeryMult11(factor1, factor2, product);
    return;
  default:    // 12 or more limbs.
    NumberLengthBytes = (NumberLength + 1) * (int)sizeof(limb);
    (void)memset(Prod, 0, NumberLengthBytes);
    for (int i = 0; i < NumberLength; i++)
    {
      int32_t Nbr = (factor1 + i)->x;
      int64_t Pr = ((int64_t)Nbr * (int64_t)factor2->x) + (int64_t)Prod[0];
      int32_t MontDig = ((int32_t)Pr * MontgomeryMultN[0].x) & MAX_INT_NBR;
      Pr -= (int64_t)MontDig * (int64_t)TestNbr[0].x;
      for (int j = 1; j < NumberLength; j++)
      {
        Pr = ((int64_t)Nbr * (int64_t)(factor2 + j)->x) + (int64_t)Prod[j] +
          (Pr >> BITS_PER_GROUP) -
          ((int64_t)MontDig * (int64_t)TestNbr[j].x);
        Prod[j - 1] = (int)Pr & (int)MAX_INT_NBR_U;
      }
      Prod[NumberLength - 1] = (int)(Pr >> BITS_PER_GROUP);
    }
#else
    double dLimbRange = (double)LIMB_RANGE;
    double dInvLimbRange = 1.0 / dLimbRange;
    NumberLengthBytes = NumberLength * (int)sizeof(limb);
    (void)memset(Prod, 0, NumberLengthBytes);
    for (int i = 0; i < NumberLength; i++)
    {
      int Nbr = (factor1 + i)->x;
      double dNbr = (double)Nbr;
      int low = (Nbr * factor2->x) + Prod[0];
      double dAccum = dNbr * (double)factor2->x + (double)Prod[0];
      int MontDig = (low * MontgomeryMultN[0].x) & MAX_INT_NBR;
      double dMontDig = (double)MontDig;
      dAccum -= dMontDig * (double)TestNbr[0].x;
      // At this moment dAccum is multiple of LIMB_RANGE.
      dAccum = floor((dAccum * dInvLimbRange) + 0.5);
      low = ((int)dAccum - (MontDig * TestNbr[1].x) + (Nbr * (factor2 + 1)->x) + Prod[1]) &
        MAX_INT_NBR;
      dAccum += (dNbr * (double)(factor2 + 1)->x) - (dMontDig * (double)TestNbr[1].x) +
        (double)Prod[1];
      Prod[0] = low;
      for (int j = 2; j < NumberLength; j++)
      {
        // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
        // In that case, there would be an error of +/- 1.
        if (low < HALF_INT_RANGE)
        {
          dAccum = ((dAccum + (double)FOURTH_INT_RANGE) * dInvLimbRange);
        }
        else
        {
          dAccum = ((dAccum - (double)FOURTH_INT_RANGE) * dInvLimbRange);
        }
        low = (int)(dAccum - (floor(dAccum * dInvLimbRange) * dLimbRange));
        dAccum += (dNbr * (double)(factor2 + j)->x) - (dMontDig * (double)TestNbr[j].x) +
          (double)Prod[j];
        low = (low - (MontDig * TestNbr[j].x) + (Nbr * (factor2 + j)->x) + Prod[j]) &
          MAX_INT_NBR;
        Prod[j - 1] = low;
      }
      // Casting from double to int truncates to nearest to zero,
      // so the conversion is done on non-negative numbers.
      if (low < HALF_INT_RANGE)
      {
        dAccum = (dAccum * dInvLimbRange) + 2147483648.25;
      }
      else
      {
        dAccum = (dAccum * dInvLimbRange) + 2147483647.75;
      }
      Prod[NumberLength - 1] = (unsigned int)dAccum ^ LIMB_RANGE;
    }
#endif
    if (Prod[NumberLength - 1] < 0)
    {        // Prod < 0, so perform Prod <- Prod + TestNbr.
      unsigned int carry = 0U;
      for (int idx = 0; idx < NumberLength; idx++)
      {
        carry = (unsigned int)Prod[idx] + (unsigned int)TestNbr[idx].x +
          (carry >> BITS_PER_GROUP);
        Prod[idx] = carry & MAX_VALUE_LIMB;
      }
    }
    NumberLengthBytes = NumberLength * (int)sizeof(limb);
    (void)memcpy(product, Prod, NumberLengthBytes);
    return;
  }
}

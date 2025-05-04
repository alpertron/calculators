//
// This file is part of Alpertron Calculators.
//
// Copyright 2015-2025 Dario Alejandro Alpern
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
#include <stdint.h>
#include <math.h>
#include "bignbr.h"

#define PROLOG_MULTIPLICATION_DOUBLE                                    \
  factor2_i = arr[idxFactor2 + i].x;                                    \
  factor2_iPlus1 = arr[idxFactor2 + i + 1].x;                           \
  Pr = prod_iPlus0 + ((uint64_t)factor2_i * factor1_0);                 \
  arrayAux[i].x = UintToInt((uint32_t)Pr & MAX_INT_NBR_U);              \
  Pr = prod_iPlus1 + ((uint64_t)factor2_i * factor1_1) +                \
       ((uint64_t)factor2_iPlus1 * factor1_0) + (Pr >> BITS_PER_GROUP); \
  arrayAux[i + 1].x = UintToInt((uint32_t)Pr & MAX_INT_NBR_U)           \

#define MULT_MACRO_DOUBLE(m, n, p)                                      \
  Pr = prod_iPlus##p + ((uint64_t)factor2_i * factor1_##p) +            \
    (uint64_t)factor2_iPlus1 * factor1_##n + (Pr >> BITS_PER_GROUP);    \
  prod_iPlus##m = UintToInt((uint32_t)Pr & MAX_INT_NBR_U)               \

#define EPILOG_MULTIPLICATION_DOUBLE(m, n)                              \
  Pr = ((uint64_t)factor2_iPlus1 * (uint64_t)factor1_##n) + (Pr >> BITS_PER_GROUP); \
  prod_iPlus##m = UintToInt((uint32_t)Pr & MAX_INT_NBR_U);              \
  prod_iPlus##n = Uint64ToInt(Pr >> BITS_PER_GROUP)
 
#define PROLOG_MULTIPLICATION_SINGLE(m)                                 \
  factor2_i = arr[idxFactor2 + m].x;                                    \
  Pr = prod_iPlus0 + ((uint64_t)factor2_i * (uint64_t)factor1_0);       \
  arrayAux[m].x = (uint32_t)Pr & MAX_INT_NBR_U 

#define MULT_MACRO_SINGLE(m, n)                                         \
  Pr = prod_iPlus##m + ((uint64_t)factor2_i * (uint64_t)factor1_##m) +  \
       (Pr >> BITS_PER_GROUP);                                          \
  arrayAux[n].x = UintToInt((uint32_t)Pr & MAX_INT_NBR_U)

#define EPILOG_MULTIPLICATION_SINGLE(m)                                 \
  arrayAux[m].x = Uint64ToInt(Pr >> BITS_PER_GROUP)

#define M(n)                                                            \
  uint32_t factor1_##n = arr[idxFactor1 + n].x;                         \
  uint32_t prod_iPlus##n = 0

extern limb arr[MAX_LEN/*3* FFT_THRESHOLD*/];
extern limb arrayAux[MAX_LEN/*3* FFT_THRESHOLD*/];

// Multiply two groups of nbrLen limbs. The first one starts at idxFactor1
// and the second one at idxFactor2. The 2*nbrLen limb result is stored
// starting at idxFactor1. Use arrayAux as temporary storage.
// Accumulate products by result limb.
#ifdef _USING64BITS_
static void ClassicalMult2Limbs(int idxFactor1, int idxFactor2)
{
  int i=0;
  uint32_t factor2_i;
  uint32_t factor2_iPlus1;
  uint64_t Pr;
  M(0); M(1);
  PROLOG_MULTIPLICATION_DOUBLE;
  EPILOG_MULTIPLICATION_DOUBLE(0, 1);
  arrayAux[2].x = prod_iPlus0;
  arrayAux[3].x = prod_iPlus1;
}

static void ClassicalMult3Limbs(int idxFactor1, int idxFactor2)
{
  int i=0;
  uint32_t factor2_i;
  uint32_t factor2_iPlus1;
  uint64_t Pr;
  M(0); M(1); M(2);
  PROLOG_MULTIPLICATION_DOUBLE;
  MULT_MACRO_DOUBLE(0, 1, 2);
  EPILOG_MULTIPLICATION_DOUBLE(1, 2);
  PROLOG_MULTIPLICATION_SINGLE(2);
  MULT_MACRO_SINGLE(1, 3);
  MULT_MACRO_SINGLE(2, 4);
  EPILOG_MULTIPLICATION_SINGLE(5);
}

static void ClassicalMult4Limbs(int idxFactor1, int idxFactor2)
{
  M(0); M(1); M(2); M(3);
  for (int i = 0; i < 4; i += 2)
  {
    uint32_t factor2_i;
    uint32_t factor2_iPlus1;
    uint64_t Pr;
    PROLOG_MULTIPLICATION_DOUBLE;
    MULT_MACRO_DOUBLE(0, 1, 2);
    MULT_MACRO_DOUBLE(1, 2, 3);
    EPILOG_MULTIPLICATION_DOUBLE(2, 3);
  }
  arrayAux[4].x = prod_iPlus0;
  arrayAux[5].x = prod_iPlus1;
  arrayAux[6].x = prod_iPlus2;
  arrayAux[7].x = prod_iPlus3;
}

static void ClassicalMult5Limbs(int idxFactor1, int idxFactor2)
{
  uint32_t factor2_i;
  uint64_t Pr;
  M(0); M(1); M(2); M(3); M(4);
  for (int i = 0; i < 4; i+=2)
  {
    uint32_t factor2_iPlus1;
    PROLOG_MULTIPLICATION_DOUBLE;
    MULT_MACRO_DOUBLE(0, 1, 2);
    MULT_MACRO_DOUBLE(1, 2, 3);
    MULT_MACRO_DOUBLE(2, 3, 4);
    EPILOG_MULTIPLICATION_DOUBLE(3, 4);
  }
  PROLOG_MULTIPLICATION_SINGLE(4);
  MULT_MACRO_SINGLE(1, 5);
  MULT_MACRO_SINGLE(2, 6);
  MULT_MACRO_SINGLE(3, 7);
  MULT_MACRO_SINGLE(4, 8);
  EPILOG_MULTIPLICATION_SINGLE(9);
}

static void ClassicalMult6Limbs(int idxFactor1, int idxFactor2)
{
  M(0); M(1); M(2); M(3); M(4); M(5);
  for (int i = 0; i < 6; i+=2)
  {
    uint32_t factor2_i;
    uint32_t factor2_iPlus1;
    uint64_t Pr;
    PROLOG_MULTIPLICATION_DOUBLE;
    MULT_MACRO_DOUBLE(0, 1, 2);
    MULT_MACRO_DOUBLE(1, 2, 3);
    MULT_MACRO_DOUBLE(2, 3, 4);
    MULT_MACRO_DOUBLE(3, 4, 5);
    EPILOG_MULTIPLICATION_DOUBLE(4, 5);
  }
  arrayAux[6].x = prod_iPlus0;
  arrayAux[7].x = prod_iPlus1;
  arrayAux[8].x = prod_iPlus2;
  arrayAux[9].x = prod_iPlus3;
  arrayAux[10].x = prod_iPlus4;
  arrayAux[11].x = prod_iPlus5;
}

static void ClassicalMult7Limbs(int idxFactor1, int idxFactor2)
{
  uint32_t factor2_i;
  uint64_t Pr;
  M(0); M(1); M(2); M(3); M(4); M(5); M(6);
  for (int i = 0; i < 6; i += 2)
  {
    uint32_t factor2_iPlus1;
    PROLOG_MULTIPLICATION_DOUBLE;
    MULT_MACRO_DOUBLE(0, 1, 2);
    MULT_MACRO_DOUBLE(1, 2, 3);
    MULT_MACRO_DOUBLE(2, 3, 4);
    MULT_MACRO_DOUBLE(3, 4, 5);
    MULT_MACRO_DOUBLE(4, 5, 6);
    EPILOG_MULTIPLICATION_DOUBLE(5, 6);
  }
  PROLOG_MULTIPLICATION_SINGLE(6);
  MULT_MACRO_SINGLE(1, 7);
  MULT_MACRO_SINGLE(2, 8);
  MULT_MACRO_SINGLE(3, 9);
  MULT_MACRO_SINGLE(4, 10);
  MULT_MACRO_SINGLE(5, 11);
  MULT_MACRO_SINGLE(6, 12);
  EPILOG_MULTIPLICATION_SINGLE(13);
}

static void ClassicalMult8Limbs(int idxFactor1, int idxFactor2)
{
  M(0); M(1); M(2); M(3); M(4); M(5); M(6); M(7);
  for (int i = 0; i < 8; i += 2)
  {
    uint32_t factor2_i;
    uint32_t factor2_iPlus1;
    uint64_t Pr;
    PROLOG_MULTIPLICATION_DOUBLE;
    MULT_MACRO_DOUBLE(0, 1, 2);
    MULT_MACRO_DOUBLE(1, 2, 3);
    MULT_MACRO_DOUBLE(2, 3, 4);
    MULT_MACRO_DOUBLE(3, 4, 5);
    MULT_MACRO_DOUBLE(4, 5, 6);
    MULT_MACRO_DOUBLE(5, 6, 7);
    EPILOG_MULTIPLICATION_DOUBLE(6, 7);
  }
  arrayAux[8].x = prod_iPlus0;
  arrayAux[9].x = prod_iPlus1;
  arrayAux[10].x = prod_iPlus2;
  arrayAux[11].x = prod_iPlus3;
  arrayAux[12].x = prod_iPlus4;
  arrayAux[13].x = prod_iPlus5;
  arrayAux[14].x = prod_iPlus6;
  arrayAux[15].x = prod_iPlus7;
}

static void ClassicalMult9Limbs(int idxFactor1, int idxFactor2)
{
  uint32_t factor2_i;
  uint64_t Pr;
  M(0); M(1); M(2); M(3); M(4); M(5); M(6); M(7); M(8);
  for (int i = 0; i < 8; i += 2)
  {
    uint32_t factor2_iPlus1;
    PROLOG_MULTIPLICATION_DOUBLE;
    MULT_MACRO_DOUBLE(0, 1, 2);
    MULT_MACRO_DOUBLE(1, 2, 3);
    MULT_MACRO_DOUBLE(2, 3, 4);
    MULT_MACRO_DOUBLE(3, 4, 5);
    MULT_MACRO_DOUBLE(4, 5, 6);
    MULT_MACRO_DOUBLE(5, 6, 7);
    MULT_MACRO_DOUBLE(6, 7, 8);
    EPILOG_MULTIPLICATION_DOUBLE(7, 8);
  }
  PROLOG_MULTIPLICATION_SINGLE(8);
  MULT_MACRO_SINGLE(1, 9);
  MULT_MACRO_SINGLE(2, 10);
  MULT_MACRO_SINGLE(3, 11);
  MULT_MACRO_SINGLE(4, 12);
  MULT_MACRO_SINGLE(5, 13);
  MULT_MACRO_SINGLE(6, 14);
  MULT_MACRO_SINGLE(7, 15);
  MULT_MACRO_SINGLE(8, 16);
  EPILOG_MULTIPLICATION_SINGLE(17);
}

static void ClassicalMult10Limbs(int idxFactor1, int idxFactor2)
{
  M(0); M(1); M(2); M(3); M(4); M(5); M(6); M(7); M(8); M(9);
  for (int i = 0; i < 10; i+=2)
  {
    uint32_t factor2_i;
    uint32_t factor2_iPlus1;
    uint64_t Pr;
    PROLOG_MULTIPLICATION_DOUBLE;
    MULT_MACRO_DOUBLE(0, 1, 2);
    MULT_MACRO_DOUBLE(1, 2, 3);
    MULT_MACRO_DOUBLE(2, 3, 4);
    MULT_MACRO_DOUBLE(3, 4, 5);
    MULT_MACRO_DOUBLE(4, 5, 6);
    MULT_MACRO_DOUBLE(5, 6, 7);
    MULT_MACRO_DOUBLE(6, 7, 8);
    MULT_MACRO_DOUBLE(7, 8, 9);
    EPILOG_MULTIPLICATION_DOUBLE(8, 9);
  }
  arrayAux[10].x = prod_iPlus0;
  arrayAux[11].x = prod_iPlus1;
  arrayAux[12].x = prod_iPlus2;
  arrayAux[13].x = prod_iPlus3;
  arrayAux[14].x = prod_iPlus4;
  arrayAux[15].x = prod_iPlus5;
  arrayAux[16].x = prod_iPlus6;
  arrayAux[17].x = prod_iPlus7;
  arrayAux[18].x = prod_iPlus8;
  arrayAux[19].x = prod_iPlus9;
}

static void ClassicalMult11Limbs(int idxFactor1, int idxFactor2)
{
  uint32_t factor2_i;
  uint64_t Pr;
  M(0); M(1); M(2); M(3); M(4); M(5); M(6); M(7); M(8); M(9); M(10);
  for (int i = 0; i < 10; i += 2)
  {
    uint32_t factor2_iPlus1;
    PROLOG_MULTIPLICATION_DOUBLE;
    MULT_MACRO_DOUBLE(0, 1, 2);
    MULT_MACRO_DOUBLE(1, 2, 3);
    MULT_MACRO_DOUBLE(2, 3, 4);
    MULT_MACRO_DOUBLE(3, 4, 5);
    MULT_MACRO_DOUBLE(4, 5, 6);
    MULT_MACRO_DOUBLE(5, 6, 7);
    MULT_MACRO_DOUBLE(6, 7, 8);
    MULT_MACRO_DOUBLE(7, 8, 9);
    MULT_MACRO_DOUBLE(8, 9, 10);
    EPILOG_MULTIPLICATION_DOUBLE(9, 10);
  }
  PROLOG_MULTIPLICATION_SINGLE(10);
  MULT_MACRO_SINGLE(1, 11);
  MULT_MACRO_SINGLE(2, 12);
  MULT_MACRO_SINGLE(3, 13);
  MULT_MACRO_SINGLE(4, 14);
  MULT_MACRO_SINGLE(5, 15);
  MULT_MACRO_SINGLE(6, 16);
  MULT_MACRO_SINGLE(7, 17);
  MULT_MACRO_SINGLE(8, 18);
  MULT_MACRO_SINGLE(9, 19);
  MULT_MACRO_SINGLE(10, 20);
  EPILOG_MULTIPLICATION_SINGLE(21);
}

static void ClassicalMult12Limbs(int idxFactor1, int idxFactor2)
{
  M(0); M(1); M(2); M(3); M(4); M(5); M(6); M(7); M(8); M(9);
  M(10); M(11);
  for (int i = 0; i < 12; i += 2)
  {
    uint32_t factor2_i;
    uint32_t factor2_iPlus1;
    uint64_t Pr;
    PROLOG_MULTIPLICATION_DOUBLE;
    MULT_MACRO_DOUBLE(0, 1, 2);
    MULT_MACRO_DOUBLE(1, 2, 3);
    MULT_MACRO_DOUBLE(2, 3, 4);
    MULT_MACRO_DOUBLE(3, 4, 5);
    MULT_MACRO_DOUBLE(4, 5, 6);
    MULT_MACRO_DOUBLE(5, 6, 7);
    MULT_MACRO_DOUBLE(6, 7, 8);
    MULT_MACRO_DOUBLE(7, 8, 9);
    MULT_MACRO_DOUBLE(8, 9, 10);
    MULT_MACRO_DOUBLE(9, 10, 11);
    EPILOG_MULTIPLICATION_DOUBLE(10, 11);
  }
  arrayAux[12].x = prod_iPlus0;
  arrayAux[13].x = prod_iPlus1;
  arrayAux[14].x = prod_iPlus2;
  arrayAux[15].x = prod_iPlus3;
  arrayAux[16].x = prod_iPlus4;
  arrayAux[17].x = prod_iPlus5;
  arrayAux[18].x = prod_iPlus6;
  arrayAux[19].x = prod_iPlus7;
  arrayAux[20].x = prod_iPlus8;
  arrayAux[21].x = prod_iPlus9;
  arrayAux[22].x = prod_iPlus10;
  arrayAux[23].x = prod_iPlus11;
}

static void ClassicalMult13Limbs(int idxFactor1, int idxFactor2)
{
  uint32_t factor2_i;
  uint64_t Pr;
  M(0); M(1); M(2); M(3); M(4); M(5); M(6); M(7); M(8); M(9);
  M(10); M(11); M(12);
  for (int i = 0; i < 12; i += 2)
  {
    uint32_t factor2_iPlus1;
    PROLOG_MULTIPLICATION_DOUBLE;
    MULT_MACRO_DOUBLE(0, 1, 2);
    MULT_MACRO_DOUBLE(1, 2, 3);
    MULT_MACRO_DOUBLE(2, 3, 4);
    MULT_MACRO_DOUBLE(3, 4, 5);
    MULT_MACRO_DOUBLE(4, 5, 6);
    MULT_MACRO_DOUBLE(5, 6, 7);
    MULT_MACRO_DOUBLE(6, 7, 8);
    MULT_MACRO_DOUBLE(7, 8, 9);
    MULT_MACRO_DOUBLE(8, 9, 10);
    MULT_MACRO_DOUBLE(9, 10, 11);
    MULT_MACRO_DOUBLE(10, 11, 12);
    EPILOG_MULTIPLICATION_DOUBLE(11, 12);
  }
  PROLOG_MULTIPLICATION_SINGLE(12);
  MULT_MACRO_SINGLE(1, 13);
  MULT_MACRO_SINGLE(2, 14);
  MULT_MACRO_SINGLE(3, 15);
  MULT_MACRO_SINGLE(4, 16);
  MULT_MACRO_SINGLE(5, 17);
  MULT_MACRO_SINGLE(6, 18);
  MULT_MACRO_SINGLE(7, 19);
  MULT_MACRO_SINGLE(8, 20);
  MULT_MACRO_SINGLE(9, 21);
  MULT_MACRO_SINGLE(10, 22);
  MULT_MACRO_SINGLE(11, 23);
  MULT_MACRO_SINGLE(12, 24);
  EPILOG_MULTIPLICATION_SINGLE(25);
}

static void ClassicalMult14Limbs(int idxFactor1, int idxFactor2)
{
  M(0); M(1); M(2); M(3); M(4); M(5); M(6); M(7); M(8); M(9);
  M(10); M(11); M(12); M(13);
  for (int i = 0; i < 14; i += 2)
  {
    uint32_t factor2_i;
    uint32_t factor2_iPlus1;
    uint64_t Pr;
    PROLOG_MULTIPLICATION_DOUBLE;
    MULT_MACRO_DOUBLE(0, 1, 2);
    MULT_MACRO_DOUBLE(1, 2, 3);
    MULT_MACRO_DOUBLE(2, 3, 4);
    MULT_MACRO_DOUBLE(3, 4, 5);
    MULT_MACRO_DOUBLE(4, 5, 6);
    MULT_MACRO_DOUBLE(5, 6, 7);
    MULT_MACRO_DOUBLE(6, 7, 8);
    MULT_MACRO_DOUBLE(7, 8, 9);
    MULT_MACRO_DOUBLE(8, 9, 10);
    MULT_MACRO_DOUBLE(9, 10, 11);
    MULT_MACRO_DOUBLE(10, 11, 12);
    MULT_MACRO_DOUBLE(11, 12, 13);
    EPILOG_MULTIPLICATION_DOUBLE(12, 13);
  }
  arrayAux[14].x = prod_iPlus0;
  arrayAux[15].x = prod_iPlus1;
  arrayAux[16].x = prod_iPlus2;
  arrayAux[17].x = prod_iPlus3;
  arrayAux[18].x = prod_iPlus4;
  arrayAux[19].x = prod_iPlus5;
  arrayAux[20].x = prod_iPlus6;
  arrayAux[21].x = prod_iPlus7;
  arrayAux[22].x = prod_iPlus8;
  arrayAux[23].x = prod_iPlus9;
  arrayAux[24].x = prod_iPlus10;
  arrayAux[25].x = prod_iPlus11;
  arrayAux[26].x = prod_iPlus12;
  arrayAux[27].x = prod_iPlus13;
}

static void ClassicalMult15Limbs(int idxFactor1, int idxFactor2)
{
  uint32_t factor2_i;
  uint64_t Pr;
  M(0); M(1); M(2); M(3); M(4); M(5); M(6); M(7); M(8); M(9);
  M(10); M(11); M(12); M(13); M(14);
  for (int i = 0; i < 14; i += 2)
  {
    uint32_t factor2_iPlus1;
    PROLOG_MULTIPLICATION_DOUBLE;
    MULT_MACRO_DOUBLE(0, 1, 2);
    MULT_MACRO_DOUBLE(1, 2, 3);
    MULT_MACRO_DOUBLE(2, 3, 4);
    MULT_MACRO_DOUBLE(3, 4, 5);
    MULT_MACRO_DOUBLE(4, 5, 6);
    MULT_MACRO_DOUBLE(5, 6, 7);
    MULT_MACRO_DOUBLE(6, 7, 8);
    MULT_MACRO_DOUBLE(7, 8, 9);
    MULT_MACRO_DOUBLE(8, 9, 10);
    MULT_MACRO_DOUBLE(9, 10, 11);
    MULT_MACRO_DOUBLE(10, 11, 12);
    MULT_MACRO_DOUBLE(11, 12, 13);
    MULT_MACRO_DOUBLE(12, 13, 14);
    EPILOG_MULTIPLICATION_DOUBLE(13, 14);
  }
  PROLOG_MULTIPLICATION_SINGLE(14);
  MULT_MACRO_SINGLE(1, 15);
  MULT_MACRO_SINGLE(2, 16);
  MULT_MACRO_SINGLE(3, 17);
  MULT_MACRO_SINGLE(4, 18);
  MULT_MACRO_SINGLE(5, 19);
  MULT_MACRO_SINGLE(6, 20);
  MULT_MACRO_SINGLE(7, 21);
  MULT_MACRO_SINGLE(8, 22);
  MULT_MACRO_SINGLE(9, 23);
  MULT_MACRO_SINGLE(10, 24);
  MULT_MACRO_SINGLE(11, 25);
  MULT_MACRO_SINGLE(12, 26);
  MULT_MACRO_SINGLE(13, 27);
  MULT_MACRO_SINGLE(14, 28);
  EPILOG_MULTIPLICATION_SINGLE(29);
}

static void ClassicalMult16Limbs(int idxFactor1, int idxFactor2)
{
  M(0); M(1); M(2); M(3); M(4); M(5); M(6); M(7); M(8); M(9);
  M(10); M(11); M(12); M(13); M(14); M(15);
  for (int i = 0; i < 16; i += 2)
  {
    uint32_t factor2_i;
    uint32_t factor2_iPlus1;
    uint64_t Pr;
    PROLOG_MULTIPLICATION_DOUBLE;
    MULT_MACRO_DOUBLE(0, 1, 2);
    MULT_MACRO_DOUBLE(1, 2, 3);
    MULT_MACRO_DOUBLE(2, 3, 4);
    MULT_MACRO_DOUBLE(3, 4, 5);
    MULT_MACRO_DOUBLE(4, 5, 6);
    MULT_MACRO_DOUBLE(5, 6, 7);
    MULT_MACRO_DOUBLE(6, 7, 8);
    MULT_MACRO_DOUBLE(7, 8, 9);
    MULT_MACRO_DOUBLE(8, 9, 10);
    MULT_MACRO_DOUBLE(9, 10, 11);
    MULT_MACRO_DOUBLE(10, 11, 12);
    MULT_MACRO_DOUBLE(11, 12, 13);
    MULT_MACRO_DOUBLE(12, 13, 14);
    MULT_MACRO_DOUBLE(13, 14, 15);
    EPILOG_MULTIPLICATION_DOUBLE(14, 15);
  }
  arrayAux[16].x = prod_iPlus0;
  arrayAux[17].x = prod_iPlus1;
  arrayAux[18].x = prod_iPlus2;
  arrayAux[19].x = prod_iPlus3;
  arrayAux[20].x = prod_iPlus4;
  arrayAux[21].x = prod_iPlus5;
  arrayAux[22].x = prod_iPlus6;
  arrayAux[23].x = prod_iPlus7;
  arrayAux[24].x = prod_iPlus8;
  arrayAux[25].x = prod_iPlus9;
  arrayAux[26].x = prod_iPlus10;
  arrayAux[27].x = prod_iPlus11;
  arrayAux[28].x = prod_iPlus12;
  arrayAux[29].x = prod_iPlus13;
  arrayAux[30].x = prod_iPlus14;
  arrayAux[31].x = prod_iPlus15;
}
#endif

void ClassicalMult(int idxFactor1, int idxFactor2, int nbrLen)
{
  int lenBytes;
#ifdef _USING64BITS_
  uint64_t product;
  switch (nbrLen)
  {
  case 1:
    product = (int64_t)arr[idxFactor1].x * arr[idxFactor2].x;
    // Store least significant limb.
    arr[idxFactor1].x = (int32_t)product & MAX_INT_NBR;
    // Store most significant limb.
    arr[idxFactor2].x = (int32_t)(product >> BITS_PER_GROUP) & MAX_INT_NBR;
    return;
  case 2:
    ClassicalMult2Limbs(idxFactor1, idxFactor2);
    break;
  case 3:
    ClassicalMult3Limbs(idxFactor1, idxFactor2);
    break;
  case 4:
    ClassicalMult4Limbs(idxFactor1, idxFactor2);
    break;
  case 5:
    ClassicalMult5Limbs(idxFactor1, idxFactor2);
    break;
  case 6:
    ClassicalMult6Limbs(idxFactor1, idxFactor2);
    break;
  case 7:
    ClassicalMult7Limbs(idxFactor1, idxFactor2);
    break;
  case 8:
    ClassicalMult8Limbs(idxFactor1, idxFactor2);
    break;
  case 9:
    ClassicalMult9Limbs(idxFactor1, idxFactor2);
    break;
  case 10:
    ClassicalMult10Limbs(idxFactor1, idxFactor2);
    break;
  case 11:
    ClassicalMult11Limbs(idxFactor1, idxFactor2);
    break;
  case 12:
    ClassicalMult12Limbs(idxFactor1, idxFactor2);
    break;
  case 13:
    ClassicalMult13Limbs(idxFactor1, idxFactor2);
    break;
  case 14:
    ClassicalMult14Limbs(idxFactor1, idxFactor2);
    break;
  case 15:
    ClassicalMult15Limbs(idxFactor1, idxFactor2);
    break;
  case 16:
    ClassicalMult16Limbs(idxFactor1, idxFactor2);
    break;
  default:
    break;
  }
#else
  const limb *ptrFactor1;
  const limb *ptrFactor2;
  int prodCol;
  int fact1Col;
  double dRangeLimb = (double)LIMB_RANGE;
  double dInvRangeLimb = 1.0 / dRangeLimb;
  int low = 0;              // Low limb of sums of multiplications.
  double dAccumulator = 0;  // Approximation to the sum of multiplications.
  int factor1;
  int factor2;
  for (prodCol = 0; prodCol < (2 * nbrLen) - 1; prodCol++)
  {    // Process each limb of product (least to most significant limb).
    if (prodCol < nbrLen)
    {   // Processing first half (least significant) of product.
      ptrFactor2 = &arr[idxFactor2 + prodCol];
      ptrFactor1 = &arr[idxFactor1];
      fact1Col = prodCol;
    }
    else
    {  // Processing second half (most significant) of product.
      ptrFactor2 = &arr[idxFactor2 + nbrLen - 1];
      ptrFactor1 = &arr[idxFactor1 + prodCol - nbrLen + 1];
      fact1Col = 2 * (nbrLen - 1) - prodCol;
    }
    for (; fact1Col>=0; fact1Col--)
    {
      factor1 = ptrFactor1->x;
      ptrFactor1++;
      factor2 = ptrFactor2->x;
      ptrFactor2--;
      low += factor1 * factor2;
      dAccumulator += (double)factor1 * (double)factor2;
    }
    low &= MAX_VALUE_LIMB;    // Trim extra bits.
    arrayAux[prodCol].x = low;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    if (low < HALF_INT_RANGE)
    {
      dAccumulator = floor((dAccumulator * dInvRangeLimb) + 0.25);
    }
    else
    {
      dAccumulator = floor((dAccumulator * dInvRangeLimb) - 0.25);
    }
    low = (unsigned int)(dAccumulator - floor(dAccumulator * dInvRangeLimb) * dRangeLimb);
  }
  arrayAux[prodCol].x = low;
#endif
  lenBytes = 2 * nbrLen * (int)sizeof(limb);
  (void)memcpy(&arr[idxFactor1], &arrayAux[0], lenBytes);
  return;
}

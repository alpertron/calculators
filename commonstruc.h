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
#ifndef COMMONSTRUCT_H
#define COMMONSTRUCT_H

// These defines are valid for factoring up to 10^110.
#define BITS_SIQS_HASH          14
#define MAX_NBR_FACTORS         13
#define MAX_PRIMES          150000
#define PRIME_SIEVE_DATA_LEN (MAX_PRIMES + 3)
#define MAX_LIMBS_SIQS          15
#define MAX_FACTORS_RELATION    80
#define LENGTH_OFFSET            0
#define MAX_SIEVE_LIMIT     200000
#define DEBUG_SIQS               0

typedef struct
{
  int value;
  int modsqrt;
  int Bainv2_0;
  int soln1;
  int difsoln;
} PrimeSieveData;

typedef struct
{
  int value;
  int exp1;
  int exp2;
  int exp3;
  int exp4;
  int exp5;
  int exp6;
} PrimeTrialDivisionData;

struct stSiqs
{
  int matrixBLength;
  int nbrFactorBasePrimes;
  int multiplier;
  int nbrFactorsA;
  int afact[MAX_NBR_FACTORS];
  limb Modulus[MAX_LIMBS_SIQS];
  limb TestNbr2[MAX_LIMBS_SIQS];
  limb biQuadrCoeff[MAX_LIMBS_SIQS];
  limb biLinearDelta[MAX_LIMBS_SIQS][MAX_LIMBS_SIQS];
  int Bainv2[MAX_NBR_FACTORS][PRIME_SIEVE_DATA_LEN];
  long largePrimeUpperBound;
  int aindex[MAX_NBR_FACTORS];
  PrimeSieveData primeSieveData[PRIME_SIEVE_DATA_LEN];
  PrimeTrialDivisionData primeTrialDivisionData[MAX_PRIMES + 50];
  int span;
  int indexMinFactorA;
  int threadNumber;
  int nbrThreadFinishedPolySet;
  unsigned int oldSeed;
  unsigned int newSeed;
  int NbrPolynomials;
  int SieveLimit;
  int matrixPartial[MAX_PRIMES * 8][(MAX_LIMBS_SIQS/2) + 4];
  limb vectLeftHandSide[(MAX_PRIMES*33/32) + 50][MAX_LIMBS_SIQS + 4];
  int matrixPartialHashIndex[1 << BITS_SIQS_HASH];
  int matrixB[(MAX_PRIMES*33/32) + 50][MAX_FACTORS_RELATION];
  int amodq[MAX_NBR_FACTORS];
  int smallPrimeUpperLimit;
  int firstLimit;
  int secondLimit;
  int thirdLimit;
  int vectExpParity[MAX_PRIMES + 50];
  int matrixAV[MAX_PRIMES];
  int matrixV[MAX_PRIMES];
  int matrixV1[MAX_PRIMES];
  int matrixV2[MAX_PRIMES];
  int matrixXmY[MAX_PRIMES];
  int newColumns[MAX_PRIMES];
 // Matrix that holds temporary data
  int matrixCalc3[MAX_PRIMES];
  int matrixTemp2[MAX_PRIMES];
  char primesUsed[MAX_PRIMES];
  int nbrPrimes2;
  int nbrPrimesUsed;
  BigInteger factorSiqs;
  PrimeSieveData *firstPrimeSieveData;
  BigInteger TempResult;
  unsigned char logar2;
  char threshold;
};

#define SHORT_SIEVE_SIZE (2*3*5*7)
#define SHORT_GROUP_SIZE ((2-1)*(3-1)*(5-1)*(7-1))
#define SIEVE_SIZE (2*3*5*7*11)
#define GROUP_SIZE ((2-1)*(3-1)*(5-1)*(7-1)*(11-1))
#define MAX_SIEVE_PRIME (10*SIEVE_SIZE)

struct stEcm
{
  limb *AA;
  limb *UX;
  limb *UZ;
  limb *W1;
  limb *W2;
  limb* W3;
  limb *W4;
  limb *X;
  limb *Z;
  limb *Aux1;
  limb *Aux2;
  limb *Aux3;
  limb *Aux4;
  limb *Aux5;
  limb *Aux6;
  limb *Aux7;
  limb *Aux8;
  limb *Xbak;
  limb *Zbak;
  limb *root;
  int sieveidx[GROUP_SIZE];
  limb GcdAccumulated[MAX_LEN];
  unsigned char sieve[MAX_SIEVE_PRIME];
  unsigned char isCoprime210_2310[SIEVE_SIZE];
  unsigned char ProcessExpon[(332199 + 7) / 8];
  unsigned char primes[((2 * 332199) + 3 + 7) / 8];
  limb buffer[(20+SHORT_GROUP_SIZE) * MAX_LEN];
  BigInteger Temp1;
  BigInteger Temp2;
  BigInteger Temp3;
};

struct stBigComplex
{
  BigInteger real;
  BigInteger imag;
};

#define MAX_NBR_PRIMES_4K1 31
struct stSumSquares
{
  struct stBigComplex primeDivisors[MAX_NBR_PRIMES_4K1];
  struct stBigComplex divisors[MAX_NBR_PRIMES_4K1];
  struct stBigComplex currentValue;
  int GrayCode;
  int indexes[MAX_NBR_PRIMES_4K1];
  int indexes2toIndexes[MAX_NBR_PRIMES_4K1];
  int currentExp[MAX_NBR_PRIMES_4K1];
  int currentExpGray[MAX_NBR_PRIMES_4K1];
  int foundSumSquares[1000000];
  int* ptrFoundSumSquares[1000];
  int nbrIndexes;
  int nbrIndexes2;
  struct BigInteger productOtherDivisors;
  struct BigInteger temp1;
  struct BigInteger temp2;
  struct BigInteger temp3;
  bool initPending;
  bool manyDecompositions;
  bool twoOddMultiplicity;
};

struct stPoly
{
  int polyLifted[1000000];
  int polyNonRepeatedFactors[1000000];
  int values[1000000];
  int denom[1000000];
  int valuesPrime[1000000];
  int poly1[1000000];
  int poly2[1000000];
  int poly3[1000000];
  int poly4[1000000];
  int poly5[1000000];
  int polyS[1000000];
  int polyT[1000000];
  int polyMultTemp[1000000];
  int polyLiftedNew[1000000];
  int polyBackup[1000000];
  int polyA[1000000];
  int polyB[1000000];
  int polyC[1000000];
  int polyD[1000000];
  int polySqFreeFact[1000000];
  int polyInteger[1000000];
  int LastAnswerPoly[1000000];
  int polyMultM[1000000];
  int polyMultT[1000000];
  int revDividend[1000000];
  int inverseDivisor[1000000];
  int polyTmp[1000000];
  int polyInv[1000000];
};

struct stDivisors
{
  int currentExp[50];
  int currentExpGray[50];
  BigInteger divisor;
  int foundDivisors[1000000];
  int* ptrFoundDivisors[1000];
};

struct stTrialDivision
{
  BigInteger power[20];
  BigInteger cofactor;
  BigInteger quotient;
  BigInteger temp;
};

struct stQuad
{
  BigInteger Solution1[400];
  BigInteger Solution2[400];
  BigInteger Increment[400];
};

// Sum of squares and divisors are used at the same time.
struct stK
{
  struct stSumSquares sumSquares;
  struct stDivisors divisors;
};

union uCommon
{
  struct stSiqs siqs;
  struct stEcm ecm;
  struct stTrialDivision trialDiv;
  struct stQuad quad;
  struct stK k;
  struct stPoly poly;
};

extern union uCommon common;
#endif

void ShowSquareModP(char* pOutput);

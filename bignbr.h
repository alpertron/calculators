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
#ifndef _BIGNBR_H
#define _BIGNBR_H
#include <stdint.h>
#include <stdbool.h>
#define MAX_LEN_MULT  25000  // 200000 digits
#ifdef FACTORIZATION_APP
#define MAX_LEN       25000  // 200000 digits
#else
#define MAX_LEN        2500  // 20000 digits
#endif
#define BITS_PER_GROUP         31
#define BITS_PER_GROUP_MINUS_1 30
#define HALF_INT_RANGE        0x40000000
#define HALF_INT_RANGE_U      0x40000000U
#define FOURTH_INT_RANGE      0x20000000
#define FOURTH_INT_RANGE_U    0x20000000U
#define MAX_INT_NBR           0x7FFFFFFF
#define MAX_INT_NBR_U         0x7FFFFFFFU
#define LIMB_RANGE            0x80000000U
#define SMALL_NUMBER_BOUND    32768
#define SMALL_PRIMES_ARRLEN    9592  // Number of primes less than 100000.
#define LOG_2            0.69314718055994531
#define LOG_3            1.09861228866810969

struct mylimb
{
  int x;
};
#define MAX_VALUE_LIMB   0x7FFFFFFFU
typedef struct mylimb limb;

enum eSign
{
  SIGN_POSITIVE = 0,
  SIGN_NEGATIVE,
};

enum eNbrCached
{
  NBR_NOT_CACHED = 0,
  NBR_READY_TO_BE_CACHED,
  NBR_CACHED,
};

typedef struct BigInteger
{
  int nbrLimbs;
  enum eSign sign;
  limb limbs[MAX_LEN];
} BigInteger;

typedef struct BigRational
{
  BigInteger numerator;
  BigInteger denominator;
} BigRational;

extern limb TestNbr[MAX_LEN];
extern limb MontgomeryMultR2[MAX_LEN];
extern limb MontgomeryMultR1[MAX_LEN];
extern int NumberLength;
extern int NumberLengthR1;
extern int groupLen;
extern enum eNbrCached MontgomeryMultNCached;
extern enum eNbrCached TestNbrCached;
extern enum eNbrCached CustomNbrCached;
extern const limb *CustomNbrAddr;
extern limb MontgomeryMultN[MAX_LEN];
extern int smallPrimes[SMALL_PRIMES_ARRLEN+1];
#ifdef __EMSCRIPTEN__
extern int percentageBPSW;
#endif

void multiply(const limb* factor1, const limb* factor2, limb* result,
  int len, int* pResultLen);
void multiplyWithBothLen(const limb* factor1, const limb* factor2, limb* result,
  int len1, int len2, int* pResultLen);
void squareRoot(const limb *argument, limb *sqRoot, int len, int *pLenSqRoot);
void Dec2Bin(const char *decimal, limb *binary, int digits, int *bitGroups);
void Bin2Dec(char **ppDecimal, const limb *binary, int nbrLimbs, int groupLength);
void Bin2Hex(char** ppDecimal, const limb *binary, int nbrLimbs, int groupLength);
void long2dec(char **pOutput, uint64_t nbr);
void int2dec(char **pOutput, int nbr);
void int2hex(char **pOutput, int nbr);
void GetMontgomeryParms(int len);
void GetMontgomeryParmsPowerOf2(int powerOf2);
void AddBigNbrModN(const limb *Nbr1, const limb *Nbr2, limb *Sum, const limb *TestNbr, int NumberLength);
void SubtBigNbrModN(const limb *Nbr1, const limb *Nbr2, limb *Sum, const limb *TestNbr, int NumberLength);
#define SubtBigNbrMod(Nbr1, Nbr2, Sum) SubtBigNbrModN(Nbr1, Nbr2, Sum, TestNbr, NumberLength) 
void modmult(const limb *factor1, const limb *factor2, limb *product);
void modmultInt(limb *factorBig, int factorInt, limb *result);
void modmultIntExtended(limb* factorBig, int factorInt, limb* result, const limb* pTestNbr, int nbrLen);
void endBigModmult(const limb* prodNotAdjusted, const limb *origProd, limb* product);
#define AddBigNbrMod(Nbr1, Nbr2, Sum) AddBigNbrModN(Nbr1, Nbr2, Sum, TestNbr, NumberLength) 
void modPowBaseInt(int base, const limb *exp, int nbrGroupsExp, limb *power);
void modPow(const limb *base, const limb *exp, int nbrGroupsExp, limb *power);
void modPowLimb(const limb *base, const limb *exp, limb *power);
void AdjustModN(limb *Nbr, const limb *TestNbr, int NumberLength);
int fsquares(void);
void AddBigInt(const limb *pAddend1, const limb *pAddend2, limb *pSum, int nbrLimbs);
bool BigIntIsZero(const BigInteger *value);
bool BigIntIsOne(const BigInteger* value);
bool BigIntEqual(const BigInteger *value1, const BigInteger *value2);
void BigIntChSign(BigInteger *value);
void BigIntAdd(const BigInteger *pAddend1, const BigInteger *pAddend2, BigInteger *pSum);
void BigIntSubt(const BigInteger* pMinuend, const BigInteger* pSubtrahend, BigInteger* pDifference);
void BigIntNegate(const BigInteger *pSrc, BigInteger *pDest);
enum eExprErr BigIntDivide(const BigInteger *pDividend, const BigInteger *pDivisor, BigInteger *pQuotient);
void classicalDivision(const BigInteger* pDividend, const BigInteger* pDivisor,
  BigInteger* pQuotient, BigInteger* pRemainder);
enum eExprErr BigIntMultiply(const BigInteger *pFactor1, const BigInteger *pFactor2, BigInteger *pProduct);
enum eExprErr BigIntRemainder(const BigInteger* pDividend,
  const BigInteger* pDivisor, BigInteger* pRemainder);
enum eExprErr BigIntPower(const BigInteger *pBase, const BigInteger *pExponent, BigInteger *pPower);
enum eExprErr BigIntPowerIntExp(const BigInteger *pBase, int exponent, BigInteger *pPower);
enum eExprErr BigIntRoot(const BigInteger* argument, BigInteger* nthRoot, int Exponent);
void BigIntRandom(BigInteger* pLowerLimit, BigInteger* pUpperLimit, BigInteger* pRandom);
int intRandom(int firstLimit, int secondLimit);
void floordiv(const BigInteger *num, const BigInteger *den, BigInteger *result);
enum eExprErr BigIntMultiplyPower2(BigInteger* pArg, int powerOf2);
void BigIntMultiplyBy2(BigInteger *nbr);
void BigIntDivideBy2(BigInteger *nbr);
void BigInteger2Dec(char **ppDecimal, const BigInteger *pBigInt, int groupLength);
void BigInteger2Hex(char** ppDecimal, const BigInteger *pBigInt, int groupLength);
void BigIntGcd(const BigInteger *pArg1, const BigInteger *pArg2, BigInteger *pResult);
enum eExprErr BigIntLcm(const BigInteger* pArg1, const BigInteger* pArg2,
  BigInteger* pResult);
enum eExprErr BigIntMod(const BigInteger *pVal, const BigInteger *pMod, BigInteger *pRes);
void BigIntGeneralModularDivision(const BigInteger *Num, const BigInteger *Den,
  const BigInteger *mod, BigInteger *quotient);
void BigIntModularDivision(const BigInteger* Num, const BigInteger* Den,
  const BigInteger* mod, BigInteger* quotient);
void BigIntModularDivisionPower2(const BigInteger* Num, const BigInteger* Den,
  const BigInteger* mod, BigInteger* quotient);
void BigIntModularDivisionSaveTestNbr(const BigInteger* Num, const BigInteger* Den,
  const BigInteger* mod, BigInteger* quotient);
void multint(BigInteger *pResult, const BigInteger *pMult, int factor);
void multadd(BigInteger *pResult, int iMult, const BigInteger *pMult, int addend);
void addmult(BigInteger* pResult, const BigInteger* pMult1, int iMult1,
  const BigInteger* pMult2, int iMult2);
void BigIntPowerOf2(BigInteger *pResult, int expon);
int getRemainder(const BigInteger *pBigInt, int divisor);
int getQuotientAndRemainder(const BigInteger* pDividend, int divisor,
  BigInteger* pQuotient);
void subtractdivide(BigInteger *pBigInt, int subt, int divisor);
void addbigint(BigInteger *pResult, int addend);
bool TestBigNbrEqual(const BigInteger *pNbr1, const BigInteger *pNbr2);
void CopyBigInt(BigInteger *pDest, const BigInteger *pSrc);
bool ModInvBigNbr(limb *num, limb *inv, limb *mod, int NumberLength);
int modInv(int NbrMod, int currentPrime);
int getNbrLimbs(const limb *bigNbr);
void BigIntDivide2(BigInteger *pArg);
int PowerCheck(const BigInteger *pBigNbr, BigInteger *pBase);
void BigIntAnd(const BigInteger *firstArg, const BigInteger *secondArg, BigInteger *result);
void BigIntOr(const BigInteger *firstArg, const BigInteger *secondArg, BigInteger *result);
void BigIntXor(const BigInteger *firstArg, const BigInteger *secondArg, BigInteger *result);
void ConvertToTwosComplement(BigInteger *value);
#ifndef FACTORIZATION_APP
int BpswPrimalityTest(const BigInteger *pValue);
#endif
void IntArray2BigInteger(const int *ptrValues, /*@out@*/BigInteger *bigint);
void BigInteger2IntArray(/*@out@*/int *ptrValues, const BigInteger *bigint);
int IntArrayCompare(const int* ptrFirst, const int* ptrSecond);
void UncompressLimbsBigInteger(const limb *ptrValues, /*@out@*/BigInteger *bigint);
void CompressLimbsBigInteger(/*@out@*/limb *ptrValues, const BigInteger *bigint);
void ComputeInversePower2(const limb *value, /*@out@*/limb *result, /*@out@*/limb *aux);
bool BigNbrIsZero(const limb *value);
void intToBigInteger(BigInteger *bigint, int value);
void longToBigInteger(BigInteger *bigint, int64_t value);
void expBigNbr(BigInteger *bignbr, double logar);
double logBigNbr(const BigInteger *pBigNbr);
double logLimbs(const limb *pBigNbr, int nbrLimbs);
double getMantissa(const limb *ptrLimb, int nbrLimbs);
double BigInt2double(const BigInteger* value);
void LenAndLimbs2ArrLimbs(const int *ptrValues, /*@out@*/limb *bigint, int nbrLen);
void ArrLimbs2LenAndLimbs(/*@out@*/int *ptrValues, const limb *bigint, int nbrLen);
bool checkOne(const limb *value, int nbrLimbs);
bool checkMinusOne(const limb *value, int nbrLimbs);
void DivideBigNbrByMaxPowerOf2(int *pShRight, limb *number, int *pNbrLimbs);
void BigIntModularPower(const BigInteger *base, const BigInteger *exponent, BigInteger *power);
enum eExprErr BigIntGeneralModularPower(const BigInteger* base, const BigInteger* exponent,
  const BigInteger* mod, BigInteger* power);

void ChSignBigNbr(limb *nbr, int length);
void ChSignBigNbrB(limb *nbr, int length);
void AddBigNbr(const limb *pNbr1, const limb *pNbr2, limb *pSum, int nbrLen);
void SubtractBigNbr(const limb *pNbr1, const limb *pNbr2, limb *pDiff, int nbrLen);
void AddBigNbrB(const limb *pNbr1, const limb *pNbr2, limb *pSum, int nbrLen);
void SubtractBigNbrB(const limb *pNbr1, const limb *pNbr2, limb *pDiff, int nbrLen);
void AddBigIntModN(const limb *pNbr1, const limb *pNbr2, limb *pSum, const limb *pMod,
  int nbrLen);
void SubtractBigNbrModN(const limb *pNbr1, const limb *pNbr2, limb *pDiff,
  const limb *pMod, int nbrLen);
void MultBigNbrByInt(const limb *bigFactor, int factor, limb *bigProd, int nbrLen);
void MultBigNbrByIntB(const limb *bigFactor, int factor, limb *bigProd, int nbrLen);
void DivBigNbrByInt(const limb *pDividend, int divisor, limb *pQuotient, int nbrLen);
int RemDivBigNbrByInt(const limb *pDividend, int divisor, int nbrLen);
void MultBigNbr(const limb *pFactor1, const limb *pFactor2, limb *pProd, int nbrLen);
void MultBigNbrComplete(const limb *pFactor1, const limb *pFactor2, limb *pProd,
  int nbrLen);
void IntToBigNbr(int value, limb *bigNbr, int nbrLength);
int BigNbrToBigInt(const BigInteger *pBigNbr, limb *pBigInt);
void BigIntToBigNbr(BigInteger *pBigNbr, const limb *pBigInt, int nbrLenBigInt);
void GcdBigNbr(const limb *pNbr1, const limb *pNbr2, limb *pGcd, int nbrLen);
void AdjustBigIntModN(limb *Nbr, const limb *Mod, int nbrLen);
void MultBigNbrModN(const limb *Nbr1, limb *Nbr2, limb *Prod, const limb *Mod, int nbrLen);
void MultBigNbrByIntModN(limb *Nbr1, int Nbr2, limb *Prod, const limb *Mod, int nbrLen);
int intDoubleModPow(int NbrMod, int Expon, int currentPrime);
void ModInvBigInt(const limb *num, limb *inv, const limb *mod, int nbrLenBigInt);
void IntToBigNbr(int value, limb *bigNbr, int nbrLength);
int JacobiSymbol(int upper, int lower);
int BigIntJacobiSymbol(const BigInteger *upper, const BigInteger *lower);
void DivideBigNbrByMaxPowerOf4(int *pPower4, limb *value, int *pNbrLimbs);
void smallmodmult(int factor1, int factor2, limb *product, int mod);
void MontgomeryMult(const limb* factor1, const limb* factor2, limb* product);
void fftMultiplication(const limb *factor1, const limb *factor2, limb *result,
  int len1, int len2, int *pResultLen);

typedef void(*mmCback)(void);
void GaussianGCD(BigInteger *realA, BigInteger *imagA, BigInteger *realB, BigInteger *imagB,
  BigInteger *realGcd, BigInteger *imagGcd, BigInteger *temp1, BigInteger *temp2);
void QuaternionGCD(BigInteger *scalarA, BigInteger *vecIA, BigInteger *vecJA, BigInteger *vecKA,
  BigInteger *scalarB, BigInteger *vecIB, BigInteger *vecJB, BigInteger *vecKB,
  BigInteger *scalarGcd, BigInteger *vecIGcd, BigInteger *vecJGcd, BigInteger *vecKGcd,
  BigInteger *temp1, BigInteger *temp2, BigInteger *temp3, BigInteger *temp4);
void MultiplyQuaternionBy2(BigInteger *scalar, BigInteger *vecI, BigInteger *vecJ, BigInteger *vecK);
void DivideQuaternionBy2(BigInteger *scalar, BigInteger *vecI, BigInteger *vecJ, BigInteger *vecK);

typedef void (*showGenericFunc)(void);
void showRationalNoParenOverGeneric(const BigRational* rat,
  showGenericFunc showGeneric);
double BigRational2double(const BigRational* value);
void BigRationalAdd(const BigRational* pAddend1,
  const BigRational* pAddend2, BigRational* pSum);
void BigRationalSubt(const BigRational* pAddend1,
  const BigRational* pAddend2, BigRational* pSum);
void BigRationalNegate(BigRational* pSrc, const BigRational* pDest);
void BigRationalMultiply(const BigRational* pFactor1,
  const BigRational* pFactor2, BigRational* pProduct);
void BigRationalDivide(const BigRational* pDividend,
  const BigRational* pDivisor, BigRational* pQuotient);
void BigRationalMultiplyByInt(const BigRational* pFactor1,
  int factor2, BigRational* pProduct);
void BigRationalDivideByInt(const BigRational* pDividend,
  int divisor, BigRational* pQuotient);
void MultiplyRationalBySqrtRational(BigRational* RatPart, BigRational* SqrPart);
bool BigRationalSquareRoot(BigRational* RatArgum, BigRational* RatSqRoot);
void ForceDenominatorPositive(BigRational* rat);
void showRational(const BigRational* rat);
void showRationalNoParen(const BigRational* rat);
void showRationalOverStr(const BigRational* rat, const char* str, const char *ptrTimes);
void ShowRationalAndSqrParts(const BigRational* RatPart, const BigRational* SqrPart,
  int root, const char *ptrTimes);
void showSquareRootOfRational(const BigRational* rat, int root, const char *ptrTimes);
void showPlusMinusRational(BigRational* rat);
void showRatCoeffAndPowerVar(BigRational* rat, int expon, char letter);
void copyStr(char** pptrString, const char* stringToCopy);
void computeRoot(const BigInteger* argument, BigInteger* nthRoot, int Exponent);

static inline int UintToInt(unsigned int value)
{
  return (int)value;
}
static inline int Uint64ToInt(uint64_t value)
{
  return (int)value;
}
#endif

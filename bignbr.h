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
#ifndef _BIGNBR_H
#define _BIGNBR_H
#ifdef FACTORIZATION_APP
#define MAX_LEN 25000        // 200000 digits
#else
#define MAX_LEN 2500         // 20000 digits
#endif
#define BITS_PER_GROUP 31
#define BITS_PER_INT_GROUP 31
#define HALF_INT_RANGE (1 << (BITS_PER_INT_GROUP - 1))
#define MAX_INT_NBR ((int)((1U << BITS_PER_INT_GROUP)-1))
#define LIMB_RANGE (1U<<BITS_PER_GROUP)
#define SMALL_NUMBER_BOUND 32768
#define SMALL_PRIMES_ARRLEN  9592  // Number of primes less than 100000.

struct mylimb
{
  int x;
};
#define MAX_VALUE_LIMB (LIMB_RANGE-1)
typedef struct mylimb limb;
typedef enum eBoolean
{
  FALSE = 0,
  TRUE,
} boolean;

enum eSign
{
  SIGN_POSITIVE = 0,
  SIGN_NEGATIVE,
};

enum
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
extern char MontgomeryMultNCached;
extern char TestNbrCached;
extern limb MontgomeryMultN[MAX_LEN];
extern int smallPrimes[SMALL_PRIMES_ARRLEN+1];
#ifdef __EMSCRIPTEN__
extern int percentageBPSW;
#endif

#include "expression.h"
#include <stdint.h>
void multiply(limb *factor1, limb *factor2, limb *result, int len, int *pResultLen);
void squareRoot(limb *argument, limb *sqRoot, int len, int *pLenSqRoot);
void Dec2Bin(char *decimal, limb *binary, int digits, int *bitGroups);
void Bin2Dec(limb *binary, char *decimal, int nbrLimbs, int groupLength);
void Bin2Hex(limb *binary, char *decimal, int nbrLimbs, int groupLength);
void long2dec(char **pOutput, uint64_t nbr);
void int2dec(char **pOutput, int nbr);
void int2hex(char **pOutput, int nbr);
void GetMontgomeryParms(int len);
void AddBigNbrModN(limb *Nbr1, limb *Nbr2, limb *Sum, limb *TestNbr, int NumberLength);
void SubtBigNbrModN(limb *Nbr1, limb *Nbr2, limb *Sum, limb *TestNbr, int NumberLength);
#define SubtBigNbrMod(Nbr1, Nbr2, Sum) SubtBigNbrModN(Nbr1, Nbr2, Sum, TestNbr, NumberLength) 
void modmult(limb *factor1, limb *factor2, limb *product);
void modmultInt(limb *factorBig, int factorInt, limb *result);
void modmultIntExtended(limb *factorBig, int factorInt, limb *result, limb *pTestNbr, int nbrLen);
#define AddBigNbrMod(Nbr1, Nbr2, Sum) AddBigNbrModN(Nbr1, Nbr2, Sum, TestNbr, NumberLength) 
void modPowBaseInt(int base, limb *exp, int nbrGroupsExp, limb *power);
void modPow(limb *base, limb *exp, int nbrGroupsExp, limb *power);
void modPowLimb(limb *base, limb *exp, limb *power);
void AdjustModN(limb *Nbr, limb *TestNbr, int NumberLength);
int fsquares(void);
void AddBigInt(limb *pAddend1, limb *pAddend2, limb *pSum, int nbrLimbs);
void SubtractBigInt(limb* pMinuend, limb* pSubtrahend, limb *pDiff, int nbrLimbs);
bool BigIntIsZero(BigInteger *value);
bool BigIntIsOne(BigInteger* value);
bool BigIntEqual(BigInteger *value1, BigInteger *value2);
void BigIntChSign(BigInteger *value);
void BigIntAdd(BigInteger *pAddend1, BigInteger *pAddend2, BigInteger *pSum);
void BigIntSubt(BigInteger* pMinuend, BigInteger* pSubtrahend, BigInteger* pDifference);
void BigIntNegate(BigInteger *pSrc, BigInteger *pDest);
enum eExprErr BigIntDivide(BigInteger *pDividend, BigInteger *pDivisor, BigInteger *pQuotient);
enum eExprErr BigIntMultiply(BigInteger *pFactor1, BigInteger *pFactor2, BigInteger *pProduct);
enum eExprErr BigIntRemainder(BigInteger *pDividend, BigInteger *pDivisor, BigInteger *pRemainder);
enum eExprErr BigIntPower(BigInteger *pBase, BigInteger *pExponent, BigInteger *pPower);
enum eExprErr BigIntPowerIntExp(BigInteger *pBase, int exponent, BigInteger *pPower);
void floordiv(BigInteger *num, BigInteger *den, BigInteger *result);
void ceildiv(BigInteger *num, BigInteger *den, BigInteger *result);
void BigIntMultiplyBy2(BigInteger *nbr);
void BigIntDivideBy2(BigInteger *nbr);
void BigInteger2Dec(BigInteger *pBigInt, char *decimal, int groupLength);
void BigInteger2Hex(BigInteger *pBigInt, char *decimal, int groupLength);
void BigIntGcd(BigInteger *pArg1, BigInteger *pArg2, BigInteger *pResult);
void BigIntGeneralModularDivision(BigInteger *Num, BigInteger *Den, BigInteger *mod, BigInteger *quotient);
void BigIntModularDivision(BigInteger *Num, BigInteger *Den, BigInteger *mod, BigInteger *quotient);
void BigIntModularDivisionPower2(BigInteger *Num, BigInteger *Den, BigInteger *mod, BigInteger *quotient);
void BigIntModularDivisionSaveTestNbr(BigInteger *Num, BigInteger *Den, BigInteger *mod, BigInteger *quotient);
void multint(BigInteger *pResult, BigInteger *pMult, int factor);
void multadd(BigInteger *pResult, int iMult, BigInteger *pMult, int addend);
void addmult(BigInteger *pResult, BigInteger *pMult1, int iMult1, BigInteger *pMult2, int iMult2);
void BigIntPowerOf2(BigInteger *pResult, int expon);
int getRemainder(BigInteger *pBigInt, int divisor);
void subtractdivide(BigInteger *pBigInt, int subt, int divisor);
void addbigint(BigInteger *pResult, int addend);
boolean TestBigNbrEqual(BigInteger *pNbr1, BigInteger *pNbr2);
void CopyBigInt(BigInteger *pDest, BigInteger *pSrc);
void ModInvBigNbr(limb *num, limb *inv, limb *mod, int NumberLength);
int modInv(int NbrMod, int currentPrime);
int getNbrLimbs(limb *bigNbr);
void BigIntDivide2(BigInteger *pArg);
int PowerCheck(BigInteger *pBigNbr, BigInteger *pBase);
void BigIntAnd(BigInteger *firstArg, BigInteger *secondArg, BigInteger *result);
void BigIntOr(BigInteger *firstArg, BigInteger *secondArg, BigInteger *result);
void BigIntXor(BigInteger *firstArg, BigInteger *secondArg, BigInteger *result);
void ConvertToTwosComplement(BigInteger *value);
#ifdef FACTORIZATION_APP
bool BpswPrimalityTest(/*@in@*/BigInteger *pValue, void *vFactors);
#else
bool BpswPrimalityTest(/*@in@*/BigInteger *pValue);
#endif
void IntArray2BigInteger(/*@in@*/int *ptrValues, /*@out@*/BigInteger *bigint);
void BigInteger2IntArray(/*@out@*/int *ptrValues, /*@in@*/BigInteger *bigint);
void UncompressLimbsBigInteger(/*@in@*/limb *ptrValues, /*@out@*/BigInteger *bigint);
void CompressLimbsBigInteger(/*@out@*/limb *ptrValues, /*@in@*/BigInteger *bigint);
void NbrToLimbs(int nbr, /*@out@*/limb *limbs, int len);
void ComputeInversePower2(/*@in@*/limb *value, /*@out@*/limb *result, /*@out@*/limb *aux);
bool BigNbrIsZero(limb *value);
void intToBigInteger(BigInteger *bigint, int value);
void longToBigInteger(BigInteger *bigint, long long value);
void expBigNbr(BigInteger *bignbr, double logar);
double logBigNbr(BigInteger *pBigNbr);
double logLimbs(limb *pBigNbr, int nbrLimbs);
double getMantissa(limb *ptrLimb, int nbrLimbs);
void LenAndLimbs2ArrLimbs(/*@in@*/int *ptrValues, /*@out@*/limb *bigint, int nbrLen);
void ArrLimbs2LenAndLimbs(/*@out@*/int *ptrValues, /*@in@*/limb *bigint, int nbrLen);
bool checkOne(limb *value, int nbrLimbs);
bool checkMinusOne(limb *value, int nbrLimbs);
void DivideBigNbrByMaxPowerOf2(int *pShRight, limb *number, int *pNbrLimbs);
void BigIntModularPower(BigInteger *base, BigInteger *exponent, BigInteger *power);
enum eExprErr BigIntGeneralModularPower(BigInteger *base, BigInteger *exponent, BigInteger *mod, BigInteger *power);

void ChSignBigNbr(int *nbr, int length);
void ChSignBigNbrB(int *nbr, int length);
void AddBigNbr(int *pNbr1, int *pNbr2, int *pSum, int nbrLen);
void SubtractBigNbr(int *pNbr1, int *pNbr2, int *pDiff, int nbrLen);
void AddBigNbrB(int *pNbr1, int *pNbr2, int *pSum, int nbrLen);
void SubtractBigNbrB(int *pNbr1, int *pNbr2, int *pDiff, int nbrLen);
void AddBigIntModN(int *pNbr1, int *pNbr2, int *pSum, int *pMod, int nbrLen);
void SubtractBigNbrModN(int *pNbr1, int *pNbr2, int *pDiff, int *pMod, int nbrLen);
void MultBigNbrByInt(int *bigFactor, int factor, int *bigProd, int nbrLen);
void MultBigNbrByIntB(int *bigFactor, int factor, int *bigProd, int nbrLen);
void DivBigNbrByInt(int *pDividend, int divisor, int *pQuotient, int nbrLen);
int RemDivBigNbrByInt(int *pDividend, int divisor, int nbrLen);
void MultBigNbr(int *pFactor1, int *pFactor2, int *pProd, int nbrLen);
void MultBigNbrComplete(const int *pFactor1, const int *pFactor2, int *pProd, int nbrLen);
void IntToBigNbr(int value, int *bigNbr, int nbrLength);
int BigNbrToBigInt(BigInteger *pBigNbr, int *pBigInt);
void BigIntToBigNbr(BigInteger *pBigNbr, int *pBigInt, int nbrLenBigInt);
bool BigNbrIsZero(limb *pNbr);
void GcdBigNbr(int *pNbr1, int *pNbr2, int *pGcd, int nbrLen);
void AdjustBigIntModN(int *Nbr, int *Mod, int nbrLen);
void MultBigNbrModN(int *Nbr1, int *Nbr2, int *Prod, int *Mod, int nbrLen);
void MultBigNbrByIntModN(int *Nbr1, int Nbr2, int *Prod, int *Mod, int nbrLen);
int intDoubleModPow(int NbrMod, int Expon, int currentPrime);
void ModInvBigInt(int *num, int *inv, int *mod, int nbrLenBigInt);
void IntToBigNbr(int value, int *bigNbr, int nbrLength);
int JacobiSymbol(int upper, int lower);
int BigIntJacobiSymbol(BigInteger *upper, BigInteger *lower);
void DivideBigNbrByMaxPowerOf4(int *pPower4, limb *value, int *pNbrLimbs);
void smallmodmult(int factor1, int factor2, limb *product, int mod);
void fftMultiplication(limb *factor1, limb *factor2, limb *result, int len, int *pResultLen);

typedef void(*mmCback)(void);
void GaussianGCD(BigInteger *realA, BigInteger *imagA, BigInteger *realB, BigInteger *imagB,
  BigInteger *realGcd, BigInteger *imagGcd, BigInteger *temp1, BigInteger *temp2);
void QuaternionGCD(BigInteger *scalarA, BigInteger *vecIA, BigInteger *vecJA, BigInteger *vecKA,
  BigInteger *scalarB, BigInteger *vecIB, BigInteger *vecJB, BigInteger *vecKB,
  BigInteger *scalarGcd, BigInteger *vecIGcd, BigInteger *vecJGcd, BigInteger *vecKGcd,
  BigInteger *temp1, BigInteger *temp2, BigInteger *temp3, BigInteger *temp4);
void MultiplyQuaternionBy2(BigInteger *scalar, BigInteger *vecI, BigInteger *vecJ, BigInteger *vecK);
void DivideQuaternionBy2(BigInteger *scalar, BigInteger *vecI, BigInteger *vecJ, BigInteger *vecK);

enum eExprErr BigRationalAdd(BigRational* pAddend1, BigRational* pAddend2, BigRational* pSum);
enum eExprErr BigRationalSubt(BigRational* pAddend1, BigRational* pAddend2, BigRational* pSum);
void BigRationalNegate(BigRational* pSrc, BigRational* pDest);
enum eExprErr BigRationalMultiply(BigRational* pFactor1, BigRational* pFactor2, BigRational* pProduct);
enum eExprErr BigRationalDivide(BigRational* pDividend, BigRational* pDivisor, BigRational* pQuotient);
enum eExprErr BigRationalMultiplyByInt(BigRational* pFactor1, int factor2, BigRational* pProduct);
enum eExprErr BigRationalDivideByInt(BigRational* pDividend, int divisor, BigRational* pQuotient);
void MultiplyRationalBySqrtRational(BigRational* RatPart, BigRational* SqrPart);
int BigRationalSquareRoot(BigRational* RatArgum, BigRational* RatSqRoot);
void ForceDenominatorPositive(BigRational* rat);
void showRational(BigRational* rat);
void showRationalNoParen(BigRational* rat);
void showRationalOverStr(BigRational* rat, char* str, char *ptrTimes);
void ShowRationalAndSqrParts(BigRational* RatPart, BigRational* SqrPart, int root, char *ptrTimes);
void showSquareRootOfRational(BigRational* rat, int root, char *ptrTimes);
#endif
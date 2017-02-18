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

#define MAX_LEN 2500        // 20000 digits
#define BITS_PER_GROUP 31
#define BITS_PER_INT_GROUP 31
#define HALF_INT_RANGE (1 << (BITS_PER_INT_GROUP - 1))
#define MAX_INT_NBR ((int)((1U << BITS_PER_INT_GROUP)-1))
#define LIMB_RANGE (1U<<BITS_PER_GROUP)
#define SMALL_NUMBER_BOUND 32768
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

typedef struct BigInteger
{
  limb limbs[MAX_LEN];
  int nbrLimbs;
  enum eSign sign;
} BigInteger;

extern limb TestNbr[MAX_LEN];
extern limb MontgomeryMultR2[MAX_LEN];
extern limb MontgomeryMultR1[MAX_LEN];
extern int NumberLength, NumberLengthR1;

#include "expression.h"
void multiply(limb *factor1, limb *factor2, limb *result, int len, int *pResultLen);
void squareRoot(limb *argument, limb *sqRoot, int len, int *pLenSqRoot);
void Dec2Bin(char *decimal, limb *binary, int digits, int *bitGroups);
void Bin2Dec(limb *binary, char *decimal, int nbrLimbs, int groupLen);
void int2dec(char **pOutput, int nbr);
void GetMontgomeryParms(int len);
void AddBigNbrModN(limb *Nbr1, limb *Nbr2, limb *Sum, limb *TestNbr, int NumberLength);
void SubtBigNbrModN(limb *Nbr1, limb *Nbr2, limb *Sum, limb *TestNbr, int NumberLength);
void SubtBigNbrMod(limb *Nbr1, limb *Nbr2, limb *Sum);
void modmult(limb *factor1, limb *factor2, limb *product);
void modmultInt(limb *factorBig, int factorInt, limb *result);
void modmultIntExtended(limb *factorBig, int factorInt, limb *result, limb *pTestNbr, int nbrLen);
void AddBigNbrMod(limb *Nbr1, limb *Nbr2, limb *Sum);
void modPowBaseInt(int base, limb *exp, int nbrGroupsExp, limb *power);
void modPow(limb *base, limb *exp, int nbrGroupsExp, limb *power);
void modPowLimb(limb *base, limb *exp, limb *power);
void AdjustModN(limb *Nbr, limb *TestNbr, int NumberLength);
int fsquares(void);
void AddBigInt(limb *pAddend1, limb *pAddend2, limb *pSum, int nbrLimbs);
void SubtractBigInt(limb *pAddend1, limb *pAddend2, limb *pSum, int nbrLimbs);
void BigIntAdd(BigInteger *pAddend1, BigInteger *pAddend2, BigInteger *pSum);
void BigIntSubt(BigInteger *pAddend1, BigInteger *pAddend2, BigInteger *pSum);
void BigIntNegate(BigInteger *pSrc, BigInteger *pDest);
enum eExprErr BigIntDivide(BigInteger *pDividend, BigInteger *pDivisor, BigInteger *pQuotient);
enum eExprErr BigIntMultiply(BigInteger *pFactor1, BigInteger *pFactor2, BigInteger *pProduct);
enum eExprErr BigIntRemainder(BigInteger *pDividend, BigInteger *pDivisor, BigInteger *pRemainder);
enum eExprErr BigIntPower(BigInteger *pBase, BigInteger *pExponent, BigInteger *pPower);
enum eExprErr BigIntPowerIntExp(BigInteger *pBase, int expon, BigInteger *pPower);
void BigInteger2Dec(BigInteger *pBigInt, char *decimal, int groupLen);
void BigIntGcd(BigInteger *pArg1, BigInteger *pArg2, BigInteger *pResult);
void BigIntModularDivision(BigInteger *Num, BigInteger *Den, BigInteger *mod, BigInteger *quotient);
void multint(BigInteger *pResult, BigInteger *pMult, int iMult);
void multadd(BigInteger *pResult, int iMult, BigInteger *pMult, int addend);
void addmult(BigInteger *pResult, BigInteger *pMult1, int iMult1, BigInteger *pMult2, int iMult2);
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
int isPseudoprime(BigInteger *pBigNbr);
void UncompressBigInteger(/*@in@*/int *ptrValues, /*@out@*/BigInteger *bigint);
void CompressBigInteger(/*@out@*/int *ptrValues, /*@in@*/BigInteger *bigint);
void UncompressLimbsBigInteger(/*@in@*/limb *ptrValues, /*@out@*/BigInteger *bigint);
void CompressLimbsBigInteger(/*@out@*/limb *ptrValues, /*@in@*/BigInteger *bigint);
void NbrToLimbs(int nbr, /*@out@*/limb *limbs, int len);
void ComputeInversePower2(/*@in@*/limb *value, /*@out@*/limb *result, /*@out@*/limb *aux);
int BigNbrIsZero(limb *value);
void intToBigInteger(BigInteger *bigint, int value);
void longToBigInteger(BigInteger *bigint, long long value);
void expBigNbr(BigInteger *pBigNbr, double logar);
double logBigNbr(BigInteger *pBigNbr);
double logLimbs(limb *pBigNbr, int nbrLimbs);
double getMantissa(limb *ptrLimb, int nbrLimbs);

void ChSignBigNbr(int *nbr, int length);
void AddBigNbr(int *pNbr1, int *pNbr2, int *pSum, int nbrLen);
void SubtractBigNbr(int *pNbr1, int *pNbr2, int *pDiff, int nbrLen);
void AddBigIntModN(int *pNbr1, int *pNbr2, int *pSum, int *pMod, int nbrLen);
void SubtractBigNbrModN(int *pNbr1, int *pNbr2, int *pDiff, int *pMod, int nbrLen);
void MultBigNbrByInt(int *bigFactor, int factor, int *bigProduct, int nbrLen);
void DivBigNbrByInt(int *pDividend, int divisor, int *pQuotient, int nbrLen);
int RemDivBigNbrByInt(int *pDividend, int divisor, int nbrLen);
void MultBigNbr(int *pFactor1, int *pFactor2, int *pProd, int nbrLen);
void IntToBigNbr(int value, int *bigNbr, int nbrLength);
int BigNbrToBigInt(BigInteger *pBigNbr, int *pBigInt);
void BigIntToBigNbr(BigInteger *pBigNbr, int *pBigInt, int nbrLenBigInt);
void GcdBigNbr(int *pNbr1, int *pNbr2, int *pGcd, int nbrLen);
void AdjustBigIntModN(int *Nbr, int *Mod, int nbrLen);
void MultBigNbrModN(int Nbr1[], int Nbr2[], int Prod[], int Mod[], int nbrLen);
void MultBigNbrByIntModN(int Nbr1[], int Nbr2, int Prod[], int Mod[], int nbrLen);
int intDoubleModPow(int NbrMod, int Expon, int currentPrime);
void ModInvBigInt(int *num, int *inv, int *mod, int NumberLength);
void IntToBigNbr(int value, int *bigNbr, int nbrLength);

typedef void(*mmCback)(void);
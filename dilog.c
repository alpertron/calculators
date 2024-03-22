//
//
// This file is part of Alpertron Calculators.
//
// Copyright 2016-2021 Dario Alejandro Alpern
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
#include "bignbr.h"
#include "expression.h"
#include "factor.h"

#ifdef __EMSCRIPTEN__
extern int64_t lModularMult;
#endif

static BigInteger groupOrder;
static BigInteger subGroupOrder;
static BigInteger powSubGroupOrder;
static BigInteger powSubGroupOrder2;
static BigInteger Exponent;
static BigInteger runningExp;
static BigInteger baseExp;
static BigInteger mod;
static BigInteger logar;
static BigInteger logarMult;
static BigInteger currentExp;
static BigInteger DiscreteLog;
static BigInteger DiscreteLogPeriod;
static BigInteger base;
static BigInteger power;
static BigInteger modulus;
static BigInteger tmpBase;
static BigInteger tmp2;
static BigInteger tmp3;
static BigInteger baseModGO;
static BigInteger bigNbrA;
static BigInteger bigNbrB;
static BigInteger LastModulus;
static int ExponentsGOComputed[400];
static BigInteger nbrV[400];
static int NbrFactorsMod = 0;
static limb nbrA[MAX_LEN];
static limb nbrA2[MAX_LEN];
static limb nbrB[MAX_LEN];
static limb nbrB2[MAX_LEN];
static limb nbrR[MAX_LEN];
static limb nbrR2[MAX_LEN];
static limb nbrPower[MAX_LEN];
static limb nbrBase[MAX_LEN];
static limb baseMontg[MAX_LEN];
static limb basePHMontg[MAX_LEN];
static limb powerMontg[MAX_LEN];
static limb powerPHMontg[MAX_LEN];
static limb currPowerMontg[MAX_LEN];
static limb primRoot[MAX_LEN];
static limb primRootPwr[MAX_LEN];
static limb TestNbrOther[MAX_LEN];
static limb MontgomeryMultR1Other[MAX_LEN];
static int NumberLengthOther;
static double dN;
static char textExp[1000];
struct sFactors astFactorsGO[1000];
int factorsGO[10000];
extern int NumberLength;
static void AdjustExponent(limb *nbr, limb mult, limb add, const BigInteger *subGroupOrder);
static void ExchangeMods(void);

static void showText(const char *text)
{
  (void)strcpy(output, text);
}

static void noDiscreteLog(void)
{
  showText(lang? "No existe el logaritmo discreto": "There is no discrete logarithm");
  DiscreteLogPeriod.sign = SIGN_NEGATIVE;
}

void textErrorDilog(char **pptrOutput, enum eExprErr rc)
{
  char* pOutput = *pptrOutput;

  *pOutput = '<';
  pOutput++;
  *pOutput = 'p';
  pOutput++;
  *pOutput = '>';
  pOutput++;
  switch (rc)
  {
  case EXPR_BASE_MUST_BE_POSITIVE:
    copyStr(&pOutput, lang ? "La base debe ser mayor que cero" :
      "Base must be greater than zero");
    break;
  case EXPR_POWER_MUST_BE_POSITIVE:
    copyStr(&pOutput, lang ? "La potencia debe ser mayor que cero" :
      "Power must be greater than zero");
    break;
  case EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE:
    copyStr(&pOutput, lang ? "El módulo debe ser mayor que 1" : "Modulus must be greater than one");
    break;
  default:
    textError(&pOutput, rc);
    break;
  }
  *pOutput = '<';
  pOutput++;
  *pOutput = '/';
  pOutput++;
  *pOutput = 'p';
  pOutput++;
  *pOutput = '>';
  pOutput++;
  *pOutput = 0;    // Add terminator character.
  *pptrOutput = pOutput;
}

static void Bin2Out(char** ppDecimal, const limb* binary, int nbrLimbs, int groupLength)
{
  if (hexadecimal)
  {
    Bin2Hex(ppDecimal, binary, nbrLimbs, groupLength);
  }
  else
  {
    Bin2Dec(ppDecimal, binary, nbrLimbs, groupLength);
  }
}

static void indicateCannotComputeLog(int indexBase, int indexExp)
{
  char *ptrText = textExp;
  const struct sFactors *pstFactors = &astFactorsGO[indexBase + 1];
  copyStr(&ptrText, "Cannot compute discrete logarithm: subgroup=");
  IntArray2BigInteger(pstFactors->ptrFactor, &tmpBase);
  Bin2Out(&ptrText, tmpBase.limbs, tmpBase.nbrLimbs, groupLen);
  copyStr(&ptrText, ", exponent=");
  int2dec(&ptrText, indexExp);
  DiscreteLogPeriod.sign = SIGN_NEGATIVE;
}

// primRootPwr = base
// powerPHmont = power
static bool ComputeDLogModSubGroupOrder(int indexBase, int indexExp, 
  BigInteger *bigExp, const BigInteger *bigSubGroupOrder)
{
  // Set tmpBase to 1 in Montgomery notation.
  int lenBytes = NumberLength * (int)sizeof(limb);
  (void)memcpy(tmpBase.limbs, MontgomeryMultR1, lenBytes);
  // Set Exponent to zero.
  intToBigInteger(bigExp, 0);
  while (!TestBigNbrEqual(bigExp, bigSubGroupOrder))
  {
    if (!memcmp(tmpBase.limbs, powerPHMontg, lenBytes))
    {    // Logarithm for this subgroup has been found. Go out.
      return true;
    }
         // Set tmpBase to next power.
    modmult(tmpBase.limbs, primRootPwr, tmpBase.limbs);
         // Set next exponent.
    addbigint(bigExp, 1);
  }
  // All exponents have been tried and logarithm has not been found, so go out.
  indicateCannotComputeLog(indexBase, indexExp);
  return false;
}

// On input: 
// - primRootPwr = base
// - powerPHmont = power
// - NumberSizeBytes: number of bytes of subgroup (multiple of 4).
// - firstLimit: limit between the sets S_0 and S_1.
// - secondLimit: limit between the sets S_1 and S_2.
// - leastSignificantLimb: index to next to most significant limb.
// - mostSignificantLimb: index to next most significant limb.
// On output: Exponent: result.
static bool PollardRho(int NumberSizeBytes, double firstLimit, double secondLimit,
  int leastSignificantLimb, int mostSignificantLimb, int indexBase, int indexExp)
{
  limb addA;
  limb addB;
  limb mult1;
  // Use intermediate values of addA2, addB2 and mult2 so
  // we do not have to compute modular arithmetic for every cycle.
  limb addA2 = { 0 };
  limb addB2 = { 0 };
  limb mult2 = { 1 };
  double magnitude;
  int64_t brentK;
  int64_t brentR;
  bool EndPollardBrentRho;
  (void)memcpy(nbrPower, powerPHMontg, NumberSizeBytes);
  // Copy base ^ (groupOrder/subGroupOrder)
  (void)memcpy(nbrBase, primRootPwr, NumberSizeBytes);
  (void)memcpy(nbrR2, nbrBase, NumberSizeBytes);
  (void)memset(nbrA2, 0, NumberSizeBytes);
  (void)memset(nbrB2, 0, NumberSizeBytes);
  nbrB2[0].x = 1;
  brentR = 1;
  brentK = 0;
  EndPollardBrentRho = false;
  // Perform Brent's cycle finding loop.
  do
  {
    (void)memcpy(nbrR, nbrR2, NumberSizeBytes);
    (void)memcpy(nbrA, nbrA2, NumberSizeBytes);
    (void)memcpy(nbrB, nbrB2, NumberSizeBytes);
    addA = addA2;
    addB = addB2;
    mult1 = mult2;
    brentR *= 2;
    // Loop inside power of 2.
    do
    {
      brentK++;
      if (NumberLength == 1)
      {
        magnitude = (double)nbrR2[leastSignificantLimb].x;
      }
      else
      {
        magnitude = ((double)nbrR2[mostSignificantLimb].x * (double)LIMB_RANGE) +
          nbrR2[leastSignificantLimb].x;
      }
      if (magnitude < firstLimit)
      {    // nbrR2 is in set S_0.
        modmult(nbrR2, nbrPower, nbrR2);
        addA2.x++;
      }
      else if (magnitude < secondLimit)
      {    // nbrR2 is in set S_1.
        modmult(nbrR2, nbrR2, nbrR2);
        mult2.x *= 2;
        addA2.x *= 2;
        addB2.x *= 2;
      }
      else
      {    // nbrR2 is in set S_2.
        modmult(nbrR2, nbrBase, nbrR2);
        addB2.x++;
      }
      if ((addA2.x >= HALF_INT_RANGE) || (addB2.x >= HALF_INT_RANGE) ||
        (mult2.x >= HALF_INT_RANGE))
      {
        // Get next value of A2 as nbrA2 <- (nbrA2 * mult2 + addA2) % subGroupOrder
        AdjustExponent(nbrA2, mult2, addA2, &subGroupOrder);
        // Get next value of B2 as nbrB2 <- (nbrB2 * mult2 + addB2) % subGroupOrder
        AdjustExponent(nbrB2, mult2, addB2, &subGroupOrder);
        // Reinitialize temporary multiplier and adders.
        mult2.x = 1;
        addA2.x = 0;
        addB2.x = 0;
      }
      if (memcmp(nbrR, nbrR2, NumberSizeBytes) == 0)
      {
        EndPollardBrentRho = true;
        break;
      }
    } while (brentK < brentR);
  } while (!EndPollardBrentRho);
  ExchangeMods();                  // TestNbr <- subGroupOrder
  // Get next value of A as nbrA <- (nbrA * mult1 + addA) % subGroupOrder
  AdjustExponent(nbrA, mult1, addA, &subGroupOrder);
  // Get next value of B as nbrB <- (nbrB * mult1 + addB) % subGroupOrder
  AdjustExponent(nbrB, mult1, addB, &subGroupOrder);
  // Get next value of A2 as nbrA2 <- (nbrA * mult2 + addA2) % subGroupOrder
  AdjustExponent(nbrA2, mult2, addA2, &subGroupOrder);
  // Get next value of B2 as nbrB2 <- (nbrA * mult2 + addB2) % subGroupOrder
  AdjustExponent(nbrB2, mult2, addB2, &subGroupOrder);
  // Get next value of B as nbrB <- (nbrB2 - nbrB) % subGroupOrder
  SubtBigNbrMod(nbrB2, nbrB, nbrB);
  SubtBigNbrMod(nbrA, nbrA2, nbrA);
  if (BigNbrIsZero(nbrA))
  {     // Denominator is zero, so rho does not work.
    ExchangeMods();           // TestNbr <- modulus
    if (!ComputeDLogModSubGroupOrder(indexBase, indexExp, &Exponent, &subGroupOrder))
    {
      return false;   // Cannot compute discrete logarithm.
    }
  }
  else
  {
    // Exponent <- (nbrB / nbrA) (mod subGroupOrder)
    UncompressLimbsBigInteger(nbrA, &bigNbrA);
    UncompressLimbsBigInteger(nbrB, &bigNbrB);
    BigIntModularDivisionSaveTestNbr(&bigNbrB, &bigNbrA, &subGroupOrder, &Exponent);
    Exponent.sign = SIGN_POSITIVE;
    ExchangeMods();           // TestNbr <- modulus
  }
  return true;
}

int pepe;
// On output: logar: Computed logarithm.
static bool ComputeDiscrLogInPrimeSubgroup(int indexBase,
  double firstLimit, double secondLimit,
  int leastSignificantLimb, int mostSignificantLimb)
{
  char* ptr;
  int nbrLimbs;
  int NumberSizeBytes;
  int indexExp;
  int lenBytes;
  int maxSGExp;
  int subGroupMultiplicity = astFactorsGO[indexBase + 1].multiplicity;

  NumberLength = *astFactorsGO[indexBase + 1].ptrFactor;
  IntArray2BigInteger(astFactorsGO[indexBase + 1].ptrFactor, &subGroupOrder);
  subGroupOrder.limbs[subGroupOrder.nbrLimbs].x = 0;
  ptr = textExp;
  copyStr(&ptr, lang? "Calculando el logaritmo discreto en subgrupo de ":
    "Computing discrete logarithm in subgroup of ");
  Bin2Out(&ptr, subGroupOrder.limbs, subGroupOrder.nbrLimbs, groupLen);
  if (astFactorsGO[indexBase + 1].multiplicity > 1)
  {
    *ptr = '<';
    ptr++;
    *ptr = 's';
    ptr++;
    *ptr = 'u';
    ptr++;
    *ptr = 'p';
    ptr++;
    *ptr = '>';
    ptr++;
    int2dec(&ptr, subGroupMultiplicity);
    *ptr = '<';
    ptr++;
    *ptr = '/';
    ptr++;
    *ptr = 's';
    ptr++;
    *ptr = 'u';
    ptr++;
    *ptr = 'p';
    ptr++;
    *ptr = '>';
    ptr++;
  }
  copyStr(&ptr, lang? " elementos": " elements");
  if (NbrFactorsMod > 1)
  {
    copyStr(&ptr, lang ? " módulo " : " modulo ");
    Bin2Out(&ptr, mod.limbs, mod.nbrLimbs, groupLen);
  }
  *ptr = '.';
  ptr++;
  *ptr = 0;
  if (mod.limbs[0].x == 23)
  {
    pepe = 1;
  }
  showText(textExp);
  NumberLength = mod.nbrLimbs;
  lenBytes = NumberLength * (int)sizeof(limb);
  (void)memcpy(TestNbr, mod.limbs, lenBytes);
  NumberLengthOther = subGroupOrder.nbrLimbs;
  lenBytes = NumberLengthOther * (int)sizeof(limb);
  (void)memcpy(TestNbrOther, subGroupOrder.limbs, lenBytes);
  TestNbr[NumberLength].x = 0;
  GetMontgomeryParms(NumberLength);
  NumberSizeBytes = NumberLength * (int)sizeof(limb);
  nbrLimbs = subGroupOrder.nbrLimbs;
  dN = (double)subGroupOrder.limbs[nbrLimbs - 1].x;
  if (nbrLimbs > 1)
  {
    dN += (double)subGroupOrder.limbs[nbrLimbs - 2].x / LIMB_RANGE;
    if (nbrLimbs > 2)
    {
      dN += (double)subGroupOrder.limbs[nbrLimbs - 3].x / LIMB_RANGE / LIMB_RANGE;
    }
  }
  CopyBigInt(&baseExp, &groupOrder);
  (void)BigIntPowerIntExp(&subGroupOrder, subGroupMultiplicity, &powSubGroupOrder2);
  (void)BigIntDivide(&groupOrder, &powSubGroupOrder2, &tmpBase);
  modPow(baseMontg, tmpBase.limbs, tmpBase.nbrLimbs, primRoot);
  (void)memcpy(primRootPwr, primRoot, NumberSizeBytes);
  ExponentsGOComputed[indexBase] = 0;
  if (!memcmp(primRoot, MontgomeryMultR1, NumberSizeBytes))
  {          // Power is one, check power of power.
    intToBigInteger(&nbrV[indexBase], 0);
    modPow(powerMontg, tmpBase.limbs, tmpBase.nbrLimbs, primRoot);
    return memcmp(primRoot, MontgomeryMultR1, NumberSizeBytes) == 0;
  }
  // Get the maximum value of k such that:
  // base^((p-1)/s^k), base^((p-1)/s^(k-1)), ... , base^((p-1)/s^(k-1))
  // is congruent to 1. s = subgroup.
  for (maxSGExp = 1; maxSGExp <= subGroupMultiplicity; maxSGExp++)
  {
    ExponentsGOComputed[indexBase]++;
    modPow(primRootPwr, subGroupOrder.limbs,
      subGroupOrder.nbrLimbs, primRoot);
    if (memcmp(primRoot, MontgomeryMultR1, NumberSizeBytes) == 0)
    {          // Power is one, so check power of power.
      (void)BigIntPowerIntExp(&subGroupOrder, subGroupMultiplicity - maxSGExp,
        &powSubGroupOrder2);
      (void)BigIntDivide(&groupOrder, &powSubGroupOrder2, &tmpBase);
      modPow(powerMontg, tmpBase.limbs, tmpBase.nbrLimbs, primRoot);
      if (memcmp(primRoot, MontgomeryMultR1, NumberSizeBytes) != 0)
      {        // Power is not one, so there is no logarithm.
        return false;
      }
      break;
    }
    (void)memcpy(primRootPwr, primRoot, NumberLength * 4);
  }
  (void)memcpy(primRoot, baseMontg, NumberSizeBytes);

  intToBigInteger(&runningExp, 0);       // runningExp <- 0
  intToBigInteger(&powSubGroupOrder, 1); // powSubGroupOrder <- 1
  CopyBigInt(&currentExp, &groupOrder);
  (void)memcpy(basePHMontg, baseMontg, NumberSizeBytes);
  (void)memcpy(currPowerMontg, powerMontg, NumberSizeBytes);
  BigIntPowerIntExp(&subGroupOrder, subGroupMultiplicity - maxSGExp,
    &powSubGroupOrder2);
  (void)BigIntDivide(&currentExp, &powSubGroupOrder2, &currentExp);
  for (indexExp = 0; indexExp < maxSGExp; indexExp++)
  {
    /* PH below comes from Pohlig-Hellman algorithm */
    (void)BigIntDivide(&currentExp, &subGroupOrder, &currentExp);
    modPow(currPowerMontg, currentExp.limbs, currentExp.nbrLimbs, powerPHMontg);
    (void)BigIntDivide(&baseExp, &subGroupOrder, &baseExp);
    if ((subGroupOrder.nbrLimbs == 1) && (subGroupOrder.limbs[0].x < 20))
    {        // subGroupOrder less than 20. Get Exponent.
      if (!ComputeDLogModSubGroupOrder(indexBase, indexExp, &Exponent, &subGroupOrder))
      {      // Cannot find logarithm, so go out.
        return false;
      }
    }
    else
    {        // Use Pollard's rho method with Brent's modification. Get Exponent.
      if (!PollardRho(NumberSizeBytes, firstLimit, secondLimit,
        leastSignificantLimb, mostSignificantLimb, indexBase, indexExp))
      {
        return false;
      }
    }
    modPow(primRoot, Exponent.limbs, Exponent.nbrLimbs, tmpBase.limbs);
    (void)ModInvBigNbr(tmpBase.limbs, tmpBase.limbs, TestNbr, NumberLength);
    modmult(tmpBase.limbs, currPowerMontg, currPowerMontg);
    (void)BigIntMultiply(&Exponent, &powSubGroupOrder, &tmpBase);
    BigIntAdd(&runningExp, &tmpBase, &runningExp);
    (void)BigIntMultiply(&powSubGroupOrder, &subGroupOrder, &powSubGroupOrder);
    modPow(primRoot, subGroupOrder.limbs, subGroupOrder.nbrLimbs, tmpBase.limbs);
    (void)memcpy(primRoot, tmpBase.limbs, NumberSizeBytes);
  }
  
  CopyBigInt(&nbrV[indexBase], &runningExp);
  NumberLength = powSubGroupOrder.nbrLimbs;
  lenBytes = NumberLength * (int)sizeof(limb);
  (void)memcpy(TestNbr, powSubGroupOrder.limbs, lenBytes);
  TestNbr[NumberLength].x = 0;
  GetMontgomeryParms(NumberLength);
  for (indexExp = 0; indexExp < indexBase; indexExp++)
  {
    const struct sFactors* pstFactors;

    // Compute nbrV[indexBase] as (nbrV[indexBase] - nbrV[indexExp]) times
    // the modular inverse of PrimesGO[indexExp] raised to ExponentsGO[indexExp]
    // mod powSubGroupOrder.
    NumberLength = mod.nbrLimbs;
    BigIntSubt(&nbrV[indexBase], &nbrV[indexExp], &nbrV[indexBase]);
    (void)BigIntRemainder(&nbrV[indexBase], &powSubGroupOrder, &nbrV[indexBase]);
    if (nbrV[indexBase].sign == SIGN_NEGATIVE)
    {
      BigIntAdd(&nbrV[indexBase], &powSubGroupOrder, &nbrV[indexBase]);
    }
    pstFactors = &astFactorsGO[indexExp + 1];
    IntArray2BigInteger(pstFactors->ptrFactor, &tmpBase);
    (void)BigIntPowerIntExp(&tmpBase, ExponentsGOComputed[indexExp], &tmpBase);
    (void)BigIntRemainder(&tmpBase, &powSubGroupOrder, &tmpBase);
    NumberLength = powSubGroupOrder.nbrLimbs;
    CompressLimbsBigInteger(tmp2.limbs, &tmpBase);
    modmult(tmp2.limbs, MontgomeryMultR2, tmp2.limbs);
    if ((NumberLength > 1) || (TestNbr[0].x != 1))
    {           // If TestNbr != 1 ...
      (void)ModInvBigNbr(tmp2.limbs, tmp2.limbs, TestNbr, NumberLength);
    }
    tmpBase.limbs[0].x = 1;
    lenBytes = (NumberLength - 1) * (int)sizeof(limb);
    (void)memset(&tmpBase.limbs[1], 0, lenBytes);
    modmult(tmpBase.limbs, tmp2.limbs, tmp2.limbs);
    UncompressLimbsBigInteger(tmp2.limbs, &tmpBase);
    (void)BigIntMultiply(&tmpBase, &nbrV[indexBase], &nbrV[indexBase]);
  }
  (void)BigIntRemainder(&nbrV[indexBase], &powSubGroupOrder, &nbrV[indexBase]);
  (void)BigIntMultiply(&nbrV[indexBase], &logarMult, &tmpBase);
  BigIntAdd(&logar, &tmpBase, &logar);
  (void)BigIntMultiply(&logarMult, &powSubGroupOrder, &logarMult);
  return true;
}

static bool DiscrLogPowerPrimeSubgroup(int multiplicity, const int *ptrPrime)
{
  int expon = 1;
  IntArray2BigInteger(ptrPrime, &bigNbrB);
  if ((bigNbrB.nbrLimbs == 1) && (bigNbrB.limbs[0].x == 2))
  {            // Prime factor is 2. Base and power are odd at this moment.
    int lsbBase = base.limbs[0].x;
    int lsbPower = power.limbs[0].x;
    if (multiplicity > 1)
    {
      int mask = ((multiplicity == 2)? 3 : 7);
      expon = ((multiplicity == 2)? 2 : 3);
      if ((lsbPower & mask) == 1)
      {
        intToBigInteger(&logar, 0);
        intToBigInteger(&logarMult, ((lsbBase == 1) ? 1 : 2));
      }
      else if (((lsbPower - lsbBase) & mask) == 0)
      {
        intToBigInteger(&logar, 1);
        intToBigInteger(&logarMult, 2);
      }
      else
      {
        noDiscreteLog();
        return false;
      }
    }
  }
  for (; expon < multiplicity; expon++)
  {    // Repeated factor.
    // L = logar, LM = logarMult
    // B = base, P = power, p = prime

    // B^n = P (mod p^(k+1)) -> n = L + m*LM   m = ?
    // B^(L + m*LM) = P
    // (B^LM) ^ m = P*B^(-L)
    // B^LM = r*p^k + 1, P*B^(-L) = s*p^k + 1
    // (r*p^k + 1)^m = s*p^k + 1
    // From binomial theorem: m = s / r (mod p)
    // If r = 0 and s != 0 there is no solution.
    // If r = 0 and s = 0 do not change LM.
    int lenBytes;
    (void)BigIntPowerIntExp(&bigNbrB, expon + 1, &bigNbrA);
    NumberLength = bigNbrA.nbrLimbs;
    lenBytes = NumberLength * (int)sizeof(limb);
    (void)memcpy(TestNbr, bigNbrA.limbs, lenBytes);
    GetMontgomeryParms(NumberLength);
    (void)BigIntRemainder(&base, &bigNbrA, &tmpBase);
    CompressLimbsBigInteger(baseMontg, &tmpBase);
    modmult(baseMontg, MontgomeryMultR2, baseMontg);
    modPow(baseMontg, logarMult.limbs, logarMult.nbrLimbs, primRootPwr); // B^LM
    tmpBase.limbs[0].x = 1;   // Convert from Montgomery to standard notation.
    lenBytes = (NumberLength - 1) * (int)sizeof(limb);
    (void)memset(&tmpBase.limbs[1], 0, lenBytes);
    modmult(primRootPwr, tmpBase.limbs, primRootPwr);                    // B^LM
    (void)ModInvBigNbr(baseMontg, tmpBase.limbs, TestNbr, NumberLength);       // B^(-1)
    modPow(tmpBase.limbs, logar.limbs, logar.nbrLimbs, primRoot);        // B^(-L)
    (void)BigIntRemainder(&power, &bigNbrA, &tmpBase);
    CompressLimbsBigInteger(tmp2.limbs, &tmpBase);
    modmult(primRoot, tmp2.limbs, primRoot);                             // P*B^(-L)
    (void)BigIntDivide(&bigNbrA, &bigNbrB, &tmpBase);
    UncompressLimbsBigInteger(primRootPwr, &tmp2);
    (void)BigIntDivide(&tmp2, &tmpBase, &bigNbrA);                             // r
    UncompressLimbsBigInteger(primRoot, &baseModGO);   // Use baseMontGO as temp var.
    (void)BigIntDivide(&baseModGO, &tmpBase, &tmp2);                           // s
    if (BigIntIsZero(&bigNbrA))
    {            // r equals zero.
      if (!BigIntIsZero(&tmp2))
      {          // s does not equal zero.
        noDiscreteLog();
        return false;
      }
    }
    else
    {            // r does not equal zero.
      BigIntModularDivisionSaveTestNbr(&tmp2, &bigNbrA, &bigNbrB, &tmpBase);          // m
      (void)BigIntMultiply(&tmpBase, &logarMult, &tmp2);
      BigIntAdd(&logar, &tmp2, &logar);
      (void)BigIntMultiply(&logarMult, &bigNbrB, &logarMult);
    }
  }
  return true;
}

void DiscreteLogarithm(void)
{
  double firstLimit;
  double secondLimit;

#ifdef __EMSCRIPTEN__
  lModularMult = 0;
#endif
  intToBigInteger(&nbrV[0], 0);
  NumberLength = modulus.nbrLimbs;
  if ((LastModulus.nbrLimbs == 0) || !TestBigNbrEqual(&LastModulus, &modulus))
  {
    char* ptrToFactorDec;
    BigInteger2IntArray(nbrToFactor, &modulus);
    ptrToFactorDec = tofactorDec;
    Bin2Out(&ptrToFactorDec, modulus.limbs, modulus.nbrLimbs, groupLen);
    factor(&modulus, nbrToFactor, factorsMod, astFactorsMod);
    NbrFactorsMod = astFactorsMod[0].multiplicity;
  }
  intToBigInteger(&DiscreteLog, 0);       // DiscreteLog <- 0
  intToBigInteger(&DiscreteLogPeriod, 1); // DiscreteLogPeriod <- 1
  for (int index = 1; index <= NbrFactorsMod; index++)
  {
    int mostSignificantLimb;
    int leastSignificantLimb;
    int NbrFactors;
    const int *ptrPrime;
    int multiplicity;
    int lenBytes;
    bool rc;

    ptrPrime = astFactorsMod[index].ptrFactor;
    NumberLength = *ptrPrime;
    IntArray2BigInteger(ptrPrime, &groupOrder);
    (void)BigIntRemainder(&base, &groupOrder, &tmpBase);
    if (BigIntIsZero(&tmpBase))
    {     // modulus and base are not relatively prime.
      int ctr;
      multiplicity = astFactorsMod[index].multiplicity;
      CopyBigInt(&bigNbrA, &power);
      for (ctr = multiplicity; ctr > 0; ctr--)
      {
        (void)BigIntRemainder(&bigNbrA, &groupOrder, &bigNbrB);
        if (!BigIntIsZero(&bigNbrB))
        {    // Exit loop if integer division cannot be performed
          break;
        }
        (void)BigIntDivide(&bigNbrA, &groupOrder, &bigNbrB);
        CopyBigInt(&bigNbrA, &bigNbrB);
      }
      if (ctr == 0)
      {  // Power is multiple of prime^exp.
        continue;
      }
      // Compute prime^multiplicity.
      (void)BigIntPowerIntExp(&groupOrder, multiplicity, &tmp2);
      (void)BigIntRemainder(&base, &tmp2, &tmpBase);
      // Get tentative exponent.
      ctr = multiplicity - ctr;
      intToBigInteger(&bigNbrB, ctr);   // Convert exponent to big integer.
      NumberLength = tmp2.nbrLimbs;
      lenBytes = (NumberLength + 1) * (int)sizeof(limb);
      (void)memcpy(TestNbr, tmp2.limbs, lenBytes);
      GetMontgomeryParms(NumberLength);
      BigIntModularPower(&tmpBase, &bigNbrB, &bigNbrA);
      (void)BigIntRemainder(&power, &tmp2, &bigNbrB);
      BigIntSubt(&bigNbrA, &bigNbrB, &bigNbrA);
      if ((bigNbrA.nbrLimbs == 1) && (bigNbrA.limbs[0].x == 0))
      {
        intToBigInteger(&DiscreteLog, ctr);     // DiscreteLog <- exponent
        intToBigInteger(&DiscreteLogPeriod, 0); // DiscreteLogPeriod <- 0
        continue;
      }
      noDiscreteLog();
      return;
    }
    else
    {     // modulus and base are relatively prime.
      (void)BigIntRemainder(&power, &groupOrder, &bigNbrB);
      if (BigIntIsZero(&bigNbrB))
      {   // power is multiple of prime. Error.
        noDiscreteLog();
        return;
      }
    }
    CompressLimbsBigInteger(baseMontg, &tmpBase);
    (void)BigIntRemainder(&power, &groupOrder, &tmpBase);
    CompressLimbsBigInteger(powerMontg, &tmpBase);
    // Compute group order as the prime minus 1.
    groupOrder.limbs[0].x--;
    showText(lang? "Calculando el logaritmo discreto...":
      "Computing discrete logarithm...");
    BigInteger2IntArray(nbrToFactor, &groupOrder);
    factor(&groupOrder, nbrToFactor, factorsGO, astFactorsGO);  // factor groupOrder.
    NbrFactors = astFactorsGO[0].multiplicity;
    NumberLength = *ptrPrime;
    IntArray2BigInteger(ptrPrime, &mod);
    intToBigInteger(&logar, 0);     // logar <- 0
    intToBigInteger(&logarMult, 1); // logarMult <- 1
    NumberLength = mod.nbrLimbs;
    lenBytes = NumberLength * (int)sizeof(limb);
    (void)memcpy(TestNbr, mod.limbs, lenBytes);
    TestNbr[NumberLength].x = 0;
    GetMontgomeryParms(NumberLength);
#if 0
    char *ptrText = textExp;
    copyStr(&ptrText, "<p>NumberLength (2) = ");
    int2dec(&ptrText, NumberLength);
    copyStr(&ptrText, "</p>");
    DiscreteLogPeriod.sign = SIGN_NEGATIVE;
    return;
#endif
    // Convert base and power to Montgomery notation.
    modmult(baseMontg, MontgomeryMultR2, baseMontg);
    modmult(powerMontg, MontgomeryMultR2, powerMontg);
    // Compute firstLimit as TestNbr / 3.
    mostSignificantLimb = NumberLength - 1;
    if (NumberLength == 1)
    {
      leastSignificantLimb = NumberLength - 1;
      firstLimit = (double)TestNbr[leastSignificantLimb].x / 3;
    }
    else
    {
      leastSignificantLimb = NumberLength - 2;
      firstLimit = ((double)TestNbr[mostSignificantLimb].x * (double)LIMB_RANGE +
        (double)TestNbr[leastSignificantLimb].x) / 3.0;
    }
    // Compute secondLimit as 2 * TestNbr / 3.
    secondLimit = firstLimit * 2.0;
    for (int indexBase = 0; indexBase < NbrFactors; indexBase++)
    {
      rc = ComputeDiscrLogInPrimeSubgroup(indexBase, firstLimit, secondLimit,
        leastSignificantLimb, mostSignificantLimb);
      if (!rc)
      {
        DiscreteLogPeriod.sign = SIGN_NEGATIVE;    // There is no logarithm.
        return;
      }
    }
    rc = DiscrLogPowerPrimeSubgroup(astFactorsMod[index].multiplicity, ptrPrime);
    if (!rc)
    {
      DiscreteLogPeriod.sign = SIGN_NEGATIVE;    // There is no logarithm.
      return;
    }
    // Based on logar and logarMult, compute DiscreteLog and DiscreteLogPeriod
    // using the following formulas, that can be deduced from the Chinese
    // Remainder Theorem:
    // L = logar, LM = logarMult, DL = DiscreteLog, DLP = DiscreteLogPeriod.
    // The modular implementation does not allow operating with even moduli.
    //
    // g <- gcd(LM, DLP)
    // if L%g != DL%g there is no discrete logarithm, so go out.
    // h <- LM / g
    // if h is odd:
    //   t <- (L - DL) / DLP (mod h)
    //   t <- DLP * t + DL
    // else
    //   i <- DLP / g
    //   t <- (DL - L) / LM (mod i)
    //   t <- LM * t + L
    // endif
    //   DLP <- DLP * h
    //   DL <- t % DLP

    BigIntGcd(&logarMult, &DiscreteLogPeriod, &tmpBase);
    (void)BigIntRemainder(&logar, &tmpBase, &bigNbrA);
    (void)BigIntRemainder(&DiscreteLog, &tmpBase, &bigNbrB);
    if (!TestBigNbrEqual(&bigNbrA, &bigNbrB))
    {
      noDiscreteLog();
      return;
    }
    (void)BigIntDivide(&logarMult, &tmpBase, &tmp2);
    if ((tmp2.limbs[0].x & 1) != 0)
    {     // h is odd.
      BigIntSubt(&logar, &DiscreteLog, &tmpBase);
      BigIntGcd(&tmpBase, &DiscreteLogPeriod, &tmp3);
      (void)BigIntDivide(&tmpBase, &tmp3, &tmpBase);
      (void)BigIntDivide(&DiscreteLogPeriod, &tmp3, &DiscreteLogPeriod);
      BigIntModularDivisionSaveTestNbr(&tmpBase, &DiscreteLogPeriod, &tmp2, &bigNbrA);
      (void)BigIntMultiply(&DiscreteLogPeriod, &tmp3, &DiscreteLogPeriod);
      (void)BigIntMultiply(&DiscreteLogPeriod, &bigNbrA, &tmpBase);
      BigIntAdd(&tmpBase, &DiscreteLog, &tmpBase);
    }
    else
    {     // h is even.
      (void)BigIntDivide(&DiscreteLogPeriod, &tmpBase, &bigNbrB);
      BigIntSubt(&DiscreteLog, &logar, &tmpBase);
      BigIntGcd(&tmpBase, &logarMult, &tmp3);
      (void)BigIntDivide(&tmpBase, &tmp3, &tmpBase);
      (void)BigIntDivide(&logarMult, &tmp3, &logarMult);
      BigIntModularDivisionSaveTestNbr(&tmpBase, &logarMult, &bigNbrB, &bigNbrA);
      (void)BigIntMultiply(&logarMult, &tmp3, &logarMult);
      (void)BigIntMultiply(&logarMult, &bigNbrA, &tmpBase);
      BigIntAdd(&tmpBase, &logar, &tmpBase);
    }
    (void)BigIntMultiply(&DiscreteLogPeriod, &tmp2, &DiscreteLogPeriod);
    (void)BigIntRemainder(&tmpBase, &DiscreteLogPeriod, &DiscreteLog);
  }
}

// Exchange TestMod and TestModOther.
static void ExchangeMods(void)
{
  int count;
  limb aux = { 0 };
  limb *ptrTestNbr;
  limb *ptrTestNbrOther;
  limb *ptrMontgomeryMultR1;
  limb *ptrMontgomeryMultR1Other;
  int maxNbrLength = NumberLength;
  if (NumberLengthOther > NumberLength)
  {
    maxNbrLength = NumberLengthOther;
  }
  ptrTestNbr = TestNbr;
  ptrTestNbrOther = TestNbrOther;
  ptrMontgomeryMultR1 = MontgomeryMultR1;
  ptrMontgomeryMultR1Other = MontgomeryMultR1Other;
  for (count = 0; count < maxNbrLength; count++)
  {
    aux.x = ptrTestNbr->x;
    ptrTestNbr->x = ptrTestNbrOther->x;
    ptrTestNbrOther->x = aux.x;
    ptrTestNbr++;
    ptrTestNbrOther++;
    aux.x = ptrMontgomeryMultR1->x;
    ptrMontgomeryMultR1->x = ptrMontgomeryMultR1Other->x;
    ptrMontgomeryMultR1Other->x = aux.x;
    ptrMontgomeryMultR1++;
    ptrMontgomeryMultR1Other++;
  }
  count = NumberLength;
  NumberLength = NumberLengthOther;
  NumberLengthOther = count;
  TestNbr[NumberLength].x = 0;
  MontgomeryMultR1[NumberLength].x = 0;
}

// nbr = (nbr * mult + add) % subGroupOrder
static void AdjustExponent(limb *nbr, limb mult, limb add, const BigInteger *bigSubGroupOrder)
{
  unsigned int carry;
  int nbrLimbs = bigSubGroupOrder->nbrLimbs;
  (nbr + nbrLimbs)->x = 0;
  MultBigNbrByInt(nbr, mult.x, nbr, nbrLimbs+1);
  carry = add.x;
  for (int j = 0; j<=nbrLimbs; j++)
  {
    carry += (unsigned int)nbr[j].x;
    nbr[j].x = UintToInt(carry & MAX_VALUE_LIMB);
    carry >>= BITS_PER_GROUP;
  }
  AdjustModN(nbr, bigSubGroupOrder->limbs, nbrLimbs);
}

static void generateOutput(enum eExprErr rc, int groupLength)
{
  char* ptrOutput;
  output[0] = '2';
  ptrOutput = &output[1];
  if (rc != EXPR_OK)
  {
    textErrorDilog(&ptrOutput, rc);
  }
  else
  {
    copyStr(&ptrOutput, lang ? "<p>Hallar <var>exp</var> tal que " :
      "<p>Find <var>exp</var> such that ");
    Bin2Out(&ptrOutput, base.limbs, base.nbrLimbs, groupLength);
    copyStr(&ptrOutput, "<sup><var>exp</var></sup> &equiv; ");
    Bin2Out(&ptrOutput, power.limbs, power.nbrLimbs, groupLength);
    copyStr(&ptrOutput, " (mod ");
    Bin2Out(&ptrOutput, modulus.limbs, modulus.nbrLimbs, groupLength);
    copyStr(&ptrOutput, ")</p><p>");
    if (DiscreteLogPeriod.sign == SIGN_NEGATIVE)
    {
      copyStr(&ptrOutput, lang ? "Ningún valor de <var>exp</var> satisface la congruencia.</p>" :
        "There is no such value of <var>exp</var>.</p>");
      copyStr(&ptrOutput, textExp);
    }
    else
    {
      copyStr(&ptrOutput, "<var>exp</var> = ");
      Bin2Out(&ptrOutput, DiscreteLog.limbs, DiscreteLog.nbrLimbs, groupLength);
      if (!BigIntIsZero(&DiscreteLogPeriod))
      {   // Discrete log period is not zero.
        copyStr(&ptrOutput, " + ");
        Bin2Out(&ptrOutput, DiscreteLogPeriod.limbs, DiscreteLogPeriod.nbrLimbs, groupLength);
        copyStr(&ptrOutput, "<var>k</var>");
      }
      copyStr(&ptrOutput, "</p>");
    }
  }
  copyStr(&ptrOutput, "<p>");
  showElapsedTime(&ptrOutput);
  copyStr(&ptrOutput, "</p><p>");
  copyStr(&ptrOutput, lang ? COPYRIGHT_SPANISH: COPYRIGHT_ENGLISH);
  copyStr(&ptrOutput, "</p>");
}

void dilogText(const char* baseText, const char* powerText, const char* modText, int groupLength)
{
  enum eExprErr rc;
  rc = ComputeExpression(baseText, &base);
  if ((rc == EXPR_OK) &&
    ((base.sign == SIGN_NEGATIVE) || BigIntIsZero(&base)))
  {
    rc = EXPR_BASE_MUST_BE_POSITIVE;
  }
  if (rc == EXPR_OK)
  {
    rc = ComputeExpression(powerText, &power);
  }
  if ((rc == EXPR_OK) &&
    ((power.sign == SIGN_NEGATIVE) || BigIntIsZero(&power)))
  {
    rc = EXPR_POWER_MUST_BE_POSITIVE;
  }
  if (rc == EXPR_OK)
  {
    rc = ComputeExpression(modText, &modulus);
  }
  if ((rc == EXPR_OK) &&
    ((modulus.sign == SIGN_NEGATIVE) || ((modulus.nbrLimbs == 1) && (modulus.limbs[0].x < 2))))
  {
    rc = EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE;
  }
  if (rc == EXPR_OK)
  {
    DiscreteLogarithm();
  }
  generateOutput(rc, groupLength);
}

#ifdef __EMSCRIPTEN__
EXTERNALIZE void doWork(void)
{
  int flags;
  int groupLength = 0;
  char *ptrData = inputString;
  char *ptrPower;
  char *ptrMod;
  groupLength = 0;
  originalTenthSecond = tenths();
  while (*ptrData != ',')
  {
    groupLength = (groupLength * 10) + (*ptrData - '0');
    ptrData++;
  }
  ptrData++;             // Skip comma.
  flags = *ptrData - '0';
  if (*(ptrData+1) != ',')
  {
    ptrData++;
    flags = (flags * 10) + *ptrData - '0';
  }
#ifndef lang
  lang = ((flags & 1)? true: false);
#endif
  hexadecimal = ((flags & 0x10) ? true : false);
  ptrData += 2;          // Skip flags and comma.
  ptrPower = ptrData + (int)strlen(ptrData) + 1;
  ptrMod = ptrPower + (int)strlen(ptrPower) + 1;
  dilogText(ptrData, ptrPower, ptrMod, groupLength);
  databack(output);
}
#endif

/*
This file is part of Alpertron Calculators.

Copyright 2016 Dario Alejandro Alpern

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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bignbr.h"
#include "expression.h"
#include "factor.h"

#ifdef __EMSCRIPTEN__
extern int64_t lModularMult;
#endif

static BigInteger groupOrder, subGroupOrder, powSubGroupOrder, powSubGroupOrderBak;
static BigInteger Exponent, runningExp, baseExp, mod;
static BigInteger logar, logarMult, runningExpBase;
static BigInteger currentExp;
static BigInteger DiscreteLog, DiscreteLogPeriod;
static BigInteger base, power, modulus, tmpBase, tmp2, baseModGO;
static BigInteger bigNbrA, bigNbrB;
static BigInteger LastModulus;
static int ExponentsGOComputed[400];
static BigInteger nbrV[400];
static int NbrFactorsMod = 0;
static limb nbrA[MAX_LEN];
static limb nbrA2[MAX_LEN];
static limb nbrB[MAX_LEN];
static limb nbrB2[MAX_LEN];
static limb nbrR[MAX_LEN];
static limb nbrROther[MAX_LEN];
static limb nbrR2[MAX_LEN];
static limb nbrPower[MAX_LEN];
static limb nbrBase[MAX_LEN];
static limb nbrTemp[MAX_LEN];
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
static void AdjustExponent(limb *nbr, limb mult, limb add, BigInteger *subGroupOrder);
static void ExchangeMods(void);

enum eLogMachineState
{
  BASE_PRIMITIVE_ROOT = 0,
  CALC_LOG_BASE,
  CALC_LOG_POWER,
};
static void showText(char *text)
{
  (void)strcpy(output, text);
}

void textErrorDilog(char *ptrOutput, enum eExprErr rc)
{
  char text[150];

  switch (rc)
  {
  case EXPR_BASE_MUST_BE_POSITIVE:
    (void)strcpy(text, lang ? "La base debe ser mayor que cero" :
      "Base must be greater than zero");
    break;
  case EXPR_POWER_MUST_BE_POSITIVE:
    (void)strcpy(text, lang ? "La potencia debe ser mayor que cero" :
      "Power must be greater than zero");
    break;
  case EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE:
    (void)strcpy(text, lang ? "El módulo debe ser mayor que 1" : "Modulus must be greater than one");
    break;
  default:
    textError(text, rc);
  }
  *ptrOutput++ = '<';
  *ptrOutput++ = 'p';
  *ptrOutput++ = '>';
  (void)strcpy(ptrOutput, text);
  ptrOutput += strlen(ptrOutput);
  *ptrOutput++ = '<';
  *ptrOutput++ = '/';
  *ptrOutput++ = 'p';
  *ptrOutput++ = '>';
  *ptrOutput = 0;    // Add terminator character.
}

static void indicateCannotComputeLog(int indexBase, int indexExp)
{
  char *ptrText;
  struct sFactors *pstFactors = &astFactorsGO[indexBase + 1];
  (void)strcpy(textExp, "Cannot compute discrete logarithm: subgroup=");
  IntArray2BigInteger(pstFactors->ptrFactor, &tmpBase);
  Bin2Dec(tmpBase.limbs, textExp + strlen(textExp), tmpBase.nbrLimbs, groupLen);
  (void)strcpy(textExp + strlen(textExp), ", exponent=");
  ptrText = textExp + strlen(textExp);
  int2dec(&ptrText, indexExp);
  DiscreteLogPeriod.sign = SIGN_NEGATIVE;
}

static int ComputeDLogModSubGroupOrder(int indexBase, int indexExp, BigInteger *Exponent, BigInteger *subGroupOrder)
{
  // Set tmpBase to 1 in Montgomery notation.
  (void)memcpy(tmpBase.limbs, MontgomeryMultR1, NumberLength * sizeof(limb));
  // Set Exponent to zero.
  Exponent->limbs[0].x = 0;
  Exponent->nbrLimbs = 1;
  Exponent->sign = SIGN_POSITIVE;
  for (;;)
  {
    if (TestBigNbrEqual(Exponent, subGroupOrder))
    {    // All exponents have been tried and logarithm has not been found, so go out.
      indicateCannotComputeLog(indexBase, indexExp);
      return 0;
    }
    if (!memcmp(tmpBase.limbs, powerPHMontg, NumberLength * sizeof(limb)))
    {    // Logarithm for this subgroup has been found. Go out.
      return 1;
    }
         // Set tmpBase to next power.
    modmult(tmpBase.limbs, primRootPwr, tmpBase.limbs);
         // Set next exponent.
    addbigint(Exponent, 1);
  }
}

void DiscreteLogarithm(void)
{
  int indexBase, indexExp;
  int index, expon;
  limb addA, addB, addA2, addB2;
  limb mult1, mult2;
  double magnitude, firstLimit, secondLimit;
  long long brentK, brentR;
  unsigned char EndPollardBrentRho;
  int nbrLimbs;
  struct sFactors *pstFactors;
  enum eLogMachineState logMachineState;
  char *ptr;

#ifdef __EMSCRIPTEN__
  lModularMult = 0;
#endif
  NumberLength = modulus.nbrLimbs;
  if (!TestBigNbrEqual(&LastModulus, &modulus))
  {
    BigInteger2IntArray(nbrToFactor, &modulus);
    Bin2Dec(modulus.limbs, tofactorDec, modulus.nbrLimbs, groupLen);
    factor(&modulus, nbrToFactor, factorsMod, astFactorsMod);
    NbrFactorsMod = astFactorsMod[0].multiplicity;
  }
  intToBigInteger(&DiscreteLog, 0);       // DiscreteLog <- 0
  intToBigInteger(&DiscreteLogPeriod, 1); // DiscreteLogPeriod <- 1
  for (index = 1; index <= NbrFactorsMod; index++)
  {
    int mostSignificantDword, leastSignificantDword;
    int NbrFactors;
    int *ptrPrime;
    int multiplicity;

    ptrPrime = astFactorsMod[index].ptrFactor;
    NumberLength = *ptrPrime;
    IntArray2BigInteger(ptrPrime, &groupOrder);
    BigIntRemainder(&base, &groupOrder, &tmpBase);
    if (tmpBase.nbrLimbs == 1 && tmpBase.limbs[0].x == 0)
    {     // modulus and base are not relatively prime.
      int ctr;
      multiplicity = astFactorsMod[index].multiplicity;
      CopyBigInt(&bigNbrA, &power);
      for (ctr = multiplicity; ctr > 0; ctr--)
      {
        BigIntRemainder(&bigNbrA, &groupOrder, &bigNbrB);
        if (bigNbrB.nbrLimbs != 1 || bigNbrB.limbs[0].x != 0)
        {    // Exit loop if integer division cannot be performed
          break;
        }
        BigIntDivide(&bigNbrA, &groupOrder, &bigNbrB);
        CopyBigInt(&bigNbrA, &bigNbrB);
      }
      if (ctr == 0)
      {  // Power is multiple of prime^exp.
        continue;
      }
      // Compute prime^mutliplicity.
      BigIntPowerIntExp(&groupOrder, multiplicity, &tmp2);
      BigIntRemainder(&base, &tmp2, &tmpBase);
      // Get tentative exponent.
      ctr = multiplicity - ctr;
      intToBigInteger(&bigNbrB, ctr);   // Convert exponent to big integer.
      NumberLength = tmp2.nbrLimbs;
      (void)memcpy(TestNbr, tmp2.limbs, (NumberLength + 1) * sizeof(limb));
      GetMontgomeryParms(NumberLength);
      BigIntModularPower(&tmpBase, &bigNbrB, &bigNbrA);
      BigIntRemainder(&power, &tmp2, &bigNbrB);
      BigIntSubt(&bigNbrA, &bigNbrB, &bigNbrA);
      if (bigNbrA.nbrLimbs == 1 && bigNbrA.limbs[0].x == 0)
      {
        intToBigInteger(&DiscreteLog, ctr);     // DiscreteLog <- exponent
        intToBigInteger(&DiscreteLogPeriod, 0); // DiscreteLogPeriod <- 0
        break;
      }
      showText("There is no discrete logarithm");
      DiscreteLogPeriod.sign = SIGN_NEGATIVE;
      return;
    }
    else
    {     // modulus and base are relatively prime.
      BigIntRemainder(&power, &groupOrder, &bigNbrB);
      if (bigNbrB.nbrLimbs == 1 && bigNbrB.limbs[0].x == 0)
      {   // power is multiple of prime. Error.
        showText("There is no discrete logarithm");
        DiscreteLogPeriod.sign = SIGN_NEGATIVE;
        return;
      }
    }
    CompressLimbsBigInteger(baseMontg, &tmpBase);
    BigIntRemainder(&power, &groupOrder, &tmpBase);
    CompressLimbsBigInteger(powerMontg, &tmpBase);
    // Compute group order as the prime minus 1.
    groupOrder.limbs[0].x--;
    showText("Computing discrete logarithm...");
    BigInteger2IntArray(nbrToFactor, &groupOrder);
    factor(&groupOrder, nbrToFactor, factorsGO, astFactorsGO);  // factor groupOrder.
    NbrFactors = astFactorsGO[0].multiplicity;
    NumberLength = *ptrPrime;
    IntArray2BigInteger(ptrPrime, &mod);
    intToBigInteger(&logar, 0);     // logar <- 0
    intToBigInteger(&logarMult, 1); // logarMult <- 1
    NumberLength = mod.nbrLimbs;
    (void)memcpy(TestNbr, mod.limbs, NumberLength * sizeof(limb));
    TestNbr[NumberLength].x = 0;
    //    yieldFreq = 1000000 / (NumberLength*NumberLength);
    GetMontgomeryParms(NumberLength);
#if 0
    char *ptrText = textExp;
    (void)strcpy(ptrText, "<p>NumberLength (2) = ");
    ptrText = ptrText + strlen(ptrText);
    int2dec(&ptrText, NumberLength);
    (void)strcpy(ptrText, "</p>");
    DiscreteLogPeriod.sign = SIGN_NEGATIVE;
    return;
#endif
    // Convert base and power to Montgomery notation.
    modmult(baseMontg, MontgomeryMultR2, baseMontg);
    modmult(powerMontg, MontgomeryMultR2, powerMontg);
    mostSignificantDword = NumberLength - 1;
    if (NumberLength == 1)
    {
      leastSignificantDword = NumberLength - 1;
      firstLimit = (double)TestNbr[leastSignificantDword].x / 3;
    }
    else
    {
      leastSignificantDword = NumberLength - 2;
      firstLimit = ((double)TestNbr[mostSignificantDword].x * LIMB_RANGE +
        TestNbr[leastSignificantDword].x) / 3;
    }
    secondLimit = firstLimit * 2;
    for (indexBase = 0; indexBase < NbrFactors; indexBase++)
    {
      NumberLength = *astFactorsGO[indexBase + 1].ptrFactor;
      IntArray2BigInteger(astFactorsGO[indexBase + 1].ptrFactor, &subGroupOrder);
      subGroupOrder.limbs[subGroupOrder.nbrLimbs].x = 0;
      (void)strcpy(textExp, "Computing discrete logarithm in subgroup of ");
      Bin2Dec(subGroupOrder.limbs, textExp + strlen(textExp), subGroupOrder.nbrLimbs, groupLen);
      ptr = textExp + strlen(textExp);
      if (astFactorsGO[indexBase + 1].multiplicity > 1)
      {
        *ptr++ = '<';
        *ptr++ = 's';
        *ptr++ = 'u';
        *ptr++ = 'p';
        *ptr++ = '>';
        int2dec(&ptr, astFactorsGO[indexBase + 1].multiplicity);
        *ptr++ = '<';
        *ptr++ = '/';
        *ptr++ = 's';
        *ptr++ = 'u';
        *ptr++ = 'p';
        *ptr++ = '>';
      }
      (void)strcpy(ptr, " elements.");
      showText(textExp);
      NumberLength = mod.nbrLimbs;
      (void)memcpy(TestNbr, mod.limbs, NumberLength * sizeof(limb));
      NumberLengthOther = subGroupOrder.nbrLimbs;
      (void)memcpy(TestNbrOther, subGroupOrder.limbs, NumberLengthOther * sizeof(limb));
      TestNbr[NumberLength].x = 0;
      GetMontgomeryParms(NumberLength);
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
      // Check whether base is primitive root.
      BigIntDivide(&groupOrder, &subGroupOrder, &tmpBase);
      modPow(baseMontg, tmpBase.limbs, tmpBase.nbrLimbs, primRootPwr);
      if (!memcmp(primRootPwr, MontgomeryMultR1, NumberLength * sizeof(limb)))
      {        // Power is one, so it is not a primitive root.
        logMachineState = CALC_LOG_BASE;
        // Find primitive root
        primRoot[0].x = 1;
        if (NumberLength > 1)
        {
          (void)memset(&primRoot[1], 0, (NumberLength - 1) * sizeof(limb));
        }
        do
        {
          primRoot[0].x++;
          modPow(primRoot, tmpBase.limbs, tmpBase.nbrLimbs, primRootPwr);
        } while (!memcmp(primRootPwr, MontgomeryMultR1, NumberLength * sizeof(limb)));
      }
      else
      {           // Power is not 1, so the base is a primitive root.
        logMachineState = BASE_PRIMITIVE_ROOT;
        (void)memcpy(primRoot, baseMontg, NumberLength * sizeof(limb));
      }
      for (;;)
      {                  // Calculate discrete logarithm in subgroup.
        runningExp.nbrLimbs = 1;     // runningExp <- 0
        runningExp.limbs[0].x = 0;
        runningExp.sign = SIGN_POSITIVE;
        powSubGroupOrder.nbrLimbs = 1;     // powSubGroupOrder <- 1
        powSubGroupOrder.limbs[0].x = 1;
        powSubGroupOrder.sign = SIGN_POSITIVE;
        CopyBigInt(&currentExp, &groupOrder);
        if (logMachineState == BASE_PRIMITIVE_ROOT)
        {
          (void)memcpy(basePHMontg, baseMontg, NumberLength * sizeof(limb));
          (void)memcpy(currPowerMontg, powerMontg, NumberLength * sizeof(limb));
        }
        else if (logMachineState == CALC_LOG_BASE)
        {
          (void)memcpy(basePHMontg, primRoot, NumberLength * sizeof(limb));
          (void)memcpy(currPowerMontg, baseMontg, NumberLength * sizeof(limb));
        }
        else
        {           // logMachineState == CALC_LOG_POWER
          (void)memcpy(primRoot, basePHMontg, NumberLength * sizeof(limb));
          (void)memcpy(currPowerMontg, powerMontg, NumberLength * sizeof(limb));
        }
        for (indexExp = 0; indexExp < astFactorsGO[indexBase + 1].multiplicity; indexExp++)
        {
          /* PH below comes from Pohlig-Hellman algorithm */
          BigIntDivide(&currentExp, &subGroupOrder, &currentExp);
          modPow(currPowerMontg, currentExp.limbs, currentExp.nbrLimbs, powerPHMontg);
          BigIntDivide(&baseExp, &subGroupOrder, &baseExp);
          if (subGroupOrder.nbrLimbs == 1 && subGroupOrder.limbs[0].x < 20)
          {       // subGroupOrder less than 20.
            if (!ComputeDLogModSubGroupOrder(indexBase, indexExp, &Exponent, &subGroupOrder))
            {
              return;
            }
          }
          else
          {        // Use Pollard's rho method with Brent's modification
            (void)memcpy(nbrPower, powerPHMontg, NumberLength * sizeof(limb));
            (void)memcpy(nbrBase, primRootPwr, NumberLength * sizeof(limb));
            (void)memcpy(nbrR2, nbrBase, NumberLength * sizeof(limb));
            (void)memset(nbrA2, 0, NumberLength * sizeof(limb));
            (void)memset(nbrB2, 0, NumberLength * sizeof(limb));
            nbrB2[0].x = 1;
            addA2.x = addB2.x = 0;
            mult2.x = 1;
            brentR = 1;
            brentK = 0;
            EndPollardBrentRho = FALSE;
            do
            {
              (void)memcpy(nbrR, nbrR2, NumberLength * sizeof(limb));
              (void)memcpy(nbrA, nbrA2, NumberLength * sizeof(limb));
              (void)memcpy(nbrB, nbrB2, NumberLength * sizeof(limb));
              addA = addA2;
              addB = addB2;
              mult1 = mult2;
              brentR *= 2;
              do
              {
                brentK++;
                if (NumberLength == 1)
                {
                  magnitude = (double)nbrR2[leastSignificantDword].x;
                }
                else
                {
                  magnitude = (double)nbrR2[mostSignificantDword].x * LIMB_RANGE +
                    nbrR2[leastSignificantDword].x;
                }
                if (magnitude < firstLimit)
                {
                  modmult(nbrR2, nbrPower, nbrROther);
                  addA2.x++;
                }
                else if (magnitude < secondLimit)
                {
                  modmult(nbrR2, nbrR2, nbrROther);
                  mult2.x *= 2;
                  addA2.x *= 2;
                  addB2.x *= 2;
                }
                else
                {
                  modmult(nbrR2, nbrBase, nbrROther);
                  addB2.x++;
                }
                // Exchange nbrR2 and nbrROther
                (void)memcpy(nbrTemp, nbrR2, NumberLength * sizeof(limb));
                (void)memcpy(nbrR2, nbrROther, NumberLength * sizeof(limb));
                (void)memcpy(nbrROther, nbrTemp, NumberLength * sizeof(limb));
                if (addA2.x >= (int)(LIMB_RANGE / 2) || addB2.x >= (int)(LIMB_RANGE / 2) ||
                    mult2.x >= (int)(LIMB_RANGE / 2))
                {
                  // nbrA2 <- (nbrA2 * mult2 + addA2) % subGroupOrder
                  AdjustExponent(nbrA2, mult2, addA2, &subGroupOrder);
                  // nbrB2 <- (nbrB2 * mult2 + addB2) % subGroupOrder
                  AdjustExponent(nbrB2, mult2, addB2, &subGroupOrder);
                  mult2.x = 1;
                  addA2.x = addB2.x = 0;
                }
                if (!memcmp(nbrR, nbrR2, NumberLength * sizeof(limb)))
                {
                  EndPollardBrentRho = TRUE;
                  break;
                }
              } while (brentK < brentR);
            } while (EndPollardBrentRho == FALSE);
            ExchangeMods();                  // TestNbr <- subGroupOrder
            // nbrA <- (nbrA * mult1 + addA) % subGroupOrder
            AdjustExponent(nbrA, mult1, addA, &subGroupOrder);
            // nbrB <- (nbrB * mult1 + addB) % subGroupOrder
            AdjustExponent(nbrB, mult1, addB, &subGroupOrder);
            // nbrA2 <- (nbrA * mult2 + addA2) % subGroupOrder
            AdjustExponent(nbrA2, mult2, addA2, &subGroupOrder);
            // nbrB2 <- (nbrA * mult2 + addB2) % subGroupOrder
            AdjustExponent(nbrB2, mult2, addB2, &subGroupOrder);
            // nbrB <- (nbrB2 - nbrB) % subGroupOrder
            SubtBigNbrMod(nbrB2, nbrB, nbrB);
            SubtBigNbrMod(nbrA, nbrA2, nbrA);
            if (BigNbrIsZero(nbrA))
            {     // Denominator is zero, so rho does not work.
              ExchangeMods();           // TestNbr <- modulus
              if (!ComputeDLogModSubGroupOrder(indexBase, indexExp, &Exponent, &subGroupOrder))
              {
                return;   // Cannot compute discrete logarithm.
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
          }
          modPow(primRoot, Exponent.limbs, Exponent.nbrLimbs, tmpBase.limbs);
          ModInvBigNbr(tmpBase.limbs, tmpBase.limbs, TestNbr, NumberLength);
          modmult(tmpBase.limbs, currPowerMontg, currPowerMontg);
          BigIntMultiply(&Exponent, &powSubGroupOrder, &tmpBase);
          BigIntAdd(&runningExp, &tmpBase, &runningExp);
          BigIntMultiply(&powSubGroupOrder, &subGroupOrder, &powSubGroupOrder);
          modPow(primRoot, subGroupOrder.limbs, subGroupOrder.nbrLimbs, tmpBase.limbs);
          (void)memcpy(primRoot, tmpBase.limbs, NumberLength * sizeof(limb));
        }
        if (logMachineState == BASE_PRIMITIVE_ROOT)
        {         // Discrete logarithm was determined for this subgroup.
          ExponentsGOComputed[indexBase] = astFactorsGO[indexBase + 1].multiplicity;
          break;
        }
        if (logMachineState == CALC_LOG_BASE)
        {
          CopyBigInt(&runningExpBase, &runningExp);
          logMachineState = CALC_LOG_POWER;
        }
        else
        {  // Set powSubGroupOrderBak to powSubGroupOrder.
           // if runningExpBase is not multiple of subGroupOrder,
           // discrete logarithm is runningExp/runningExpBase mod powSubGroupOrderBak.
           // Otherwise if runningExp is multiple of subGroupOrder, there is no logarithm.
           // Otherwise, divide runningExp, runnignExpBase and powSubGroupOrderBak by subGroupOrder and repeat.
          ExponentsGOComputed[indexBase] = astFactorsGO[indexBase + 1].multiplicity;
          CopyBigInt(&powSubGroupOrderBak, &powSubGroupOrder);
          do
          {
            BigIntRemainder(&runningExpBase, &subGroupOrder, &tmpBase);
            if (tmpBase.nbrLimbs > 1 || tmpBase.limbs[0].x != 0)
            {    // runningExpBase is not multiple of subGroupOrder
              BigIntModularDivisionSaveTestNbr(&runningExp, &runningExpBase, &powSubGroupOrderBak, &tmpBase);
              CopyBigInt(&runningExp, &tmpBase);
              break;
            }
            BigIntRemainder(&runningExp, &subGroupOrder, &tmpBase);
            if (tmpBase.nbrLimbs > 1 || tmpBase.limbs[0].x != 0)
            {    // runningExpBase is not multiple of subGroupOrder
              showText("There is no discrete logarithm");
              DiscreteLogPeriod.sign = SIGN_NEGATIVE;
              return;
            }
            BigIntDivide(&runningExp, &subGroupOrder, &tmpBase);
            CopyBigInt(&runningExp, &tmpBase);
            BigIntDivide(&runningExpBase, &subGroupOrder, &tmpBase);
            CopyBigInt(&runningExpBase, &tmpBase);
            BigIntDivide(&powSubGroupOrderBak, &subGroupOrder, &tmpBase);
            CopyBigInt(&powSubGroupOrderBak, &tmpBase);
            ExponentsGOComputed[indexBase]--;
            if (tmpBase.nbrLimbs == 1 && tmpBase.limbs[0].x == 1)
            {
              break;
            }
            BigIntRemainder(&runningExpBase, &subGroupOrder, &tmpBase);
          } while (tmpBase.nbrLimbs == 1 && tmpBase.limbs[0].x == 0);
          CopyBigInt(&powSubGroupOrder, &powSubGroupOrderBak);
          // The logarithm is runningExp / runningExpBase mod powSubGroupOrder
          // When powSubGroupOrder is even, we cannot use Montgomery.
          if (powSubGroupOrder.limbs[0].x & 1)
          {          // powSubGroupOrder is odd.
            BigIntModularDivisionSaveTestNbr(&runningExp, &runningExpBase, &powSubGroupOrder, &tmpBase);
            CopyBigInt(&runningExp, &tmpBase);
          }
          else
          {          // powSubGroupOrder is even (power of 2).
            NumberLength = powSubGroupOrder.nbrLimbs;
            CompressLimbsBigInteger(nbrB, &runningExpBase);
            ComputeInversePower2(nbrB, nbrA, nbrB2);  // nbrB2 is auxiliary var.
            CompressLimbsBigInteger(nbrB, &runningExp);
            multiply(nbrA, nbrB, nbrA, NumberLength, NULL);   // nbrA <- quotient.
            UncompressLimbsBigInteger(nbrA, &runningExp);
          }
          break;
        }
      }
      CopyBigInt(&nbrV[indexBase], &runningExp);
      NumberLength = powSubGroupOrder.nbrLimbs;
      (void)memcpy(TestNbr, powSubGroupOrder.limbs, NumberLength * sizeof(limb));
      TestNbr[NumberLength].x = 0;
      GetMontgomeryParms(NumberLength);
      for (indexExp = 0; indexExp < indexBase; indexExp++)
      {
        // nbrV[indexBase] <- (nbrV[indexBase] - nbrV[indexExp])*
        //                     modinv(PrimesGO[indexExp]^(ExponentsGO[indexExp]),
        //                     powSubGroupOrder)
        NumberLength = mod.nbrLimbs;
        BigIntSubt(&nbrV[indexBase], &nbrV[indexExp], &nbrV[indexBase]);
        BigIntRemainder(&nbrV[indexBase], &powSubGroupOrder, &nbrV[indexBase]);
        if (nbrV[indexBase].sign == SIGN_NEGATIVE)
        {
          BigIntAdd(&nbrV[indexBase], &powSubGroupOrder, &nbrV[indexBase]);
        }
        pstFactors = &astFactorsGO[indexExp + 1];
        IntArray2BigInteger(pstFactors->ptrFactor, &tmpBase);
        BigIntPowerIntExp(&tmpBase, ExponentsGOComputed[indexExp], &tmpBase);
        BigIntRemainder(&tmpBase, &powSubGroupOrder, &tmpBase);
        NumberLength = powSubGroupOrder.nbrLimbs;
        CompressLimbsBigInteger(tmp2.limbs, &tmpBase);
        modmult(tmp2.limbs, MontgomeryMultR2, tmp2.limbs);
        if (NumberLength > 1 || TestNbr[0].x != 1)
        {           // If TestNbr != 1 ...
          ModInvBigNbr(tmp2.limbs, tmp2.limbs, TestNbr, NumberLength);
        }
        tmpBase.limbs[0].x = 1;
        (void)memset(&tmpBase.limbs[1], 0, (NumberLength - 1) * sizeof(limb));
        modmult(tmpBase.limbs, tmp2.limbs, tmp2.limbs);
        UncompressLimbsBigInteger(tmp2.limbs, &tmpBase);
        BigIntMultiply(&tmpBase, &nbrV[indexBase], &nbrV[indexBase]);
      }
      BigIntRemainder(&nbrV[indexBase], &powSubGroupOrder, &nbrV[indexBase]);
      BigIntMultiply(&nbrV[indexBase], &logarMult, &tmpBase);
      BigIntAdd(&logar, &tmpBase, &logar);
      BigIntMultiply(&logarMult, &powSubGroupOrder, &logarMult);
    }
    multiplicity = astFactorsMod[index].multiplicity;
    IntArray2BigInteger(ptrPrime, &bigNbrB);
    expon = 1;
    if (bigNbrB.nbrLimbs == 1 && bigNbrB.limbs[0].x == 2)
    {            // Prime factor is 2. Base and power are odd at this moment.
      int lsbBase = base.limbs[0].x;
      int lsbPower = power.limbs[0].x;
      if (multiplicity > 1)
      {
        int mask = (multiplicity == 2? 3 : 7);
        expon = (multiplicity == 2 ? 2 : 3);
        if ((lsbPower & mask) == 1)
        {
          intToBigInteger(&logar, 0);
          intToBigInteger(&logarMult, (lsbBase == 1 ? 1 : 2));
        }
        else if (((lsbPower - lsbBase) & mask) == 0)
        {
          intToBigInteger(&logar, 1);
          intToBigInteger(&logarMult, 2);
        }
        else
        {
          showText("There is no discrete logarithm");
          DiscreteLogPeriod.sign = SIGN_NEGATIVE;
          return;
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
      BigIntPowerIntExp(&bigNbrB, expon + 1, &bigNbrA);
      NumberLength = bigNbrA.nbrLimbs;
      (void)memcpy(TestNbr, bigNbrA.limbs, NumberLength * sizeof(limb));
      GetMontgomeryParms(NumberLength);
      BigIntRemainder(&base, &bigNbrA, &tmpBase);
      CompressLimbsBigInteger(baseMontg, &tmpBase);
      modmult(baseMontg, MontgomeryMultR2, baseMontg);
      modPow(baseMontg, logarMult.limbs, logarMult.nbrLimbs, primRootPwr); // B^LM
      tmpBase.limbs[0].x = 1;   // Convert from Montgomery to standard notation.
      (void)memset(&tmpBase.limbs[1], 0, (NumberLength - 1) * sizeof(limb));
      modmult(primRootPwr, tmpBase.limbs, primRootPwr);                    // B^LM
      ModInvBigNbr(baseMontg, tmpBase.limbs, TestNbr, NumberLength);       // B^(-1)
      modPow(tmpBase.limbs, logar.limbs, logar.nbrLimbs, primRoot);        // B^(-L)
      BigIntRemainder(&power, &bigNbrA, &tmpBase);
      CompressLimbsBigInteger(tmp2.limbs, &tmpBase);
      modmult(primRoot, tmp2.limbs, primRoot);                             // P*B^(-L)
      BigIntDivide(&bigNbrA, &bigNbrB, &tmpBase);
      UncompressLimbsBigInteger(primRootPwr, &tmp2);
      BigIntDivide(&tmp2, &tmpBase, &bigNbrA);                             // s
      UncompressLimbsBigInteger(primRoot, &baseModGO);   // Use baseMontGO as temp var.
      BigIntDivide(&baseModGO, &tmpBase, &tmp2);                           // r
      if (bigNbrA.nbrLimbs == 1 && bigNbrA.limbs[0].x == 0)
      {            // r equals zero.
        if (tmp2.nbrLimbs != 1 || tmp2.limbs[0].x != 0)
        {          // s does not equal zero.
          showText("There is no discrete logarithm");
          DiscreteLogPeriod.sign = SIGN_NEGATIVE;
          return;
        }
      }
      else
      {            // r does not equal zero.
        BigIntModularDivisionSaveTestNbr(&tmp2, &bigNbrA, &bigNbrB, &tmpBase);          // m
        BigIntMultiply(&tmpBase, &logarMult, &tmp2);
        BigIntAdd(&logar, &tmp2, &logar);
        BigIntMultiply(&logarMult, &bigNbrB, &logarMult);
      }
    }
    // Based on logar and logarMult, compute DiscreteLog and DiscreteLogPeriod
    // using the following formulas, that can be deduced from the Chinese
    // Remainder Theorem:
    // L = logar, LM = logarMult, DL = DiscreteLog, DLP = DiscreteLogPeriod.
    // The modular implementation does not allow operating with even moduli.
    //
    // g <- gcd(LM, DLP)
    // if (L%g != DL%g) there is no discrete logarithm, so go out.
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
    BigIntRemainder(&logar, &tmpBase, &bigNbrA);
    BigIntRemainder(&DiscreteLog, &tmpBase, &bigNbrB);
    if (!TestBigNbrEqual(&bigNbrA, &bigNbrB))
    {
      showText("There is no discrete logarithm");
      DiscreteLogPeriod.sign = SIGN_NEGATIVE;
      return;
    }
    BigIntDivide(&logarMult, &tmpBase, &tmp2);
    if (tmp2.limbs[0].x & 1)
    {     // h is odd.
      BigIntSubt(&logar, &DiscreteLog, &tmpBase);
      BigIntModularDivisionSaveTestNbr(&tmpBase, &DiscreteLogPeriod, &tmp2, &bigNbrA);
      BigIntMultiply(&DiscreteLogPeriod, &bigNbrA, &tmpBase);
      BigIntAdd(&tmpBase, &DiscreteLog, &tmpBase);
    }
    else
    {     // h is even.
      BigIntDivide(&DiscreteLogPeriod, &tmpBase, &bigNbrB);
      BigIntSubt(&DiscreteLog, &logar, &tmpBase);
      BigIntModularDivisionSaveTestNbr(&tmpBase, &logarMult, &bigNbrB, &bigNbrA);
      BigIntMultiply(&logarMult, &bigNbrA, &tmpBase);
      BigIntAdd(&tmpBase, &logar, &tmpBase);
    }
    BigIntMultiply(&DiscreteLogPeriod, &tmp2, &DiscreteLogPeriod);
    BigIntRemainder(&tmpBase, &DiscreteLogPeriod, &DiscreteLog);
  }
#if 0
  textExp.setText(DiscreteLog.toString());
  textPeriod.setText(DiscreteLogPeriod.toString());
  long t = OldTimeElapsed / 1000;
  labelStatus.setText("Time elapsed: " +
    t / 86400 + "d " + (t % 86400) / 3600 + "h " + ((t % 3600) / 60) + "m " + (t % 60) +
    "s    mod mult: " + lModularMult);
#endif
}

// Exchange TestMod and TestModOther.
static void ExchangeMods(void)
{
  int count;
  limb aux;
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
static void AdjustExponent(limb *nbr, limb mult, limb add, BigInteger *subGroupOrder)
{
  unsigned int carry;
  int j;
  int nbrLimbs = subGroupOrder->nbrLimbs;
  (nbr + nbrLimbs)->x = 0;
  MultBigNbrByInt((int *)nbr, mult.x, (int *)nbr, nbrLimbs+1);
  carry = add.x;
  for (j = 0; j<=nbrLimbs; j++)
  {
    carry += nbr[j].x;
    nbr[j].x = (int)(carry & MAX_VALUE_LIMB);
    carry >>= BITS_PER_GROUP;
  }
  AdjustModN(nbr, subGroupOrder->limbs, nbrLimbs);
}

void dilogText(char *baseText, char *powerText, char *modText, int groupLength)
{
  char *ptrOutput;
  enum eExprErr rc;
  rc = ComputeExpression(baseText, 1, &base);
  if (rc == EXPR_OK)
  {
    if (base.sign == SIGN_NEGATIVE || (base.nbrLimbs == 1 && base.limbs[0].x == 0))
    {
      rc = EXPR_BASE_MUST_BE_POSITIVE;
    }
  }
  if (rc == EXPR_OK)
  {
    rc = ComputeExpression(powerText, 1, &power);
  }
  if (rc == EXPR_OK)
  {
    if (power.sign == SIGN_NEGATIVE || (power.nbrLimbs == 1 && base.limbs[0].x == 0))
    {
      rc = EXPR_POWER_MUST_BE_POSITIVE;
    }
  }
  if (rc == EXPR_OK)
  {
    rc = ComputeExpression(modText, 1, &modulus);
  }
  if (rc == EXPR_OK)
  {
    if (modulus.sign == SIGN_NEGATIVE || (modulus.nbrLimbs == 1 && modulus.limbs[0].x < 2))
    {
      rc = EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE;
    }
  }
  if (rc == EXPR_OK)
  {
    DiscreteLogarithm();
  }
  output[0] = '2';
  ptrOutput = &output[1];
  if (rc != EXPR_OK)
  {
    textErrorDilog(output + 1, rc);
    ptrOutput = output + strlen(output);
  }
  else
  {
    (void)strcpy(ptrOutput, lang?"<p>Hallar <var>exp</var> tal que ": 
                           "<p>Find <var>exp</var> such that ");
    ptrOutput += strlen(ptrOutput);
    Bin2Dec(base.limbs, ptrOutput, base.nbrLimbs, groupLength);
    ptrOutput += strlen(ptrOutput);
    strcat(ptrOutput, "<sup><var>exp</var></sup> &equiv; ");
    ptrOutput += strlen(ptrOutput);
    Bin2Dec(power.limbs, ptrOutput, power.nbrLimbs, groupLength);
    ptrOutput += strlen(ptrOutput);
    strcat(ptrOutput, " (mod ");
    ptrOutput += strlen(ptrOutput);
    Bin2Dec(modulus.limbs, ptrOutput, modulus.nbrLimbs, groupLength);
    ptrOutput += strlen(ptrOutput);
    strcat(ptrOutput, ")</p><p>");
    ptrOutput += strlen(ptrOutput);
    if (DiscreteLogPeriod.sign == SIGN_NEGATIVE)
    {
      strcat(ptrOutput, lang? "Ningún valor de <var>exp</var> satisface la congruencia.</p>":
                              "There is no such value of <var>exp</var>.</p>");
      ptrOutput += strlen(ptrOutput);
      (void)strcpy(ptrOutput, textExp);
    }
    else
    {
      strcat(ptrOutput, "<var>exp</var> = ");
      ptrOutput += strlen(ptrOutput);
      Bin2Dec(DiscreteLog.limbs, ptrOutput, DiscreteLog.nbrLimbs, groupLength);
      ptrOutput += strlen(ptrOutput);
      if (DiscreteLogPeriod.nbrLimbs != 1 || DiscreteLogPeriod.limbs[0].x != 0)
      {   // Discrete log period is not zero.
        strcat(ptrOutput, " + ");
        ptrOutput += strlen(ptrOutput);
        Bin2Dec(DiscreteLogPeriod.limbs, ptrOutput, DiscreteLogPeriod.nbrLimbs, groupLength);
        ptrOutput += strlen(ptrOutput);
        strcat(ptrOutput, "<var>k</var>");
      }
      strcat(ptrOutput, "</p>");
    }
  }
  strcat(ptrOutput, lang ? "<p>" COPYRIGHT_SPANISH "</p>" :
                           "<p>" COPYRIGHT_ENGLISH "</p>");
}

#ifdef __EMSCRIPTEN__
EXTERNALIZE void doWork(void)
{
  int flags;
  int groupLen = 0;
  char *ptrData = inputString;
  char *ptrPower, *ptrMod;
  groupLen = 0;
  while (*ptrData != ',')
  {
    groupLen = groupLen * 10 + (*ptrData++ - '0');
  }
  ptrData++;             // Skip comma.
  flags = *ptrData;
#ifndef lang  
  lang = flags & 1;
#endif
  ptrData += 2;          // Skip flags and comma.
  ptrPower = ptrData + strlen(ptrData) + 1;
  ptrMod = ptrPower + strlen(ptrPower) + 1;
  dilogText(ptrData, ptrPower, ptrMod, groupLen);
  databack(output);
}
#endif

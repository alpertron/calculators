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
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "bignbr.h"
#include "expression.h"
#include "factor.h"
#include "commonstruc.h"
#include "skiptest.h"

int yieldFreq;
int maxIndexM;
int indexM;
static int oldNbrFactors;
#ifdef __EMSCRIPTEN__
char* ptrLowerText;
char upperText[MAX_LEN*16];
char lowerText[MAX_LEN*16];
extern mmCback modmultCallback;
extern int64_t lModularMult;
extern char *ptrInputText;
int nbrPrimes;
int indexPrimes;
#endif

union uCommon common;
int64_t primeModMult;
int StepECM;
bool skipPrimality;
int nbrECM;
int nbrPrimalityTests;
int nbrSIQS;
int timeECM;
int timePrimalityTests;
int timeSIQS;
static int DegreeAurif;
static int NextEC;
static BigInteger power;
static BigInteger prime;
int *factorArr[FACTOR_ARRSIZE];
static bool foundByLehman;
static int EC;
static BigInteger Temp1;
static BigInteger Temp2;
static BigInteger Temp3;
static BigInteger Temp4;
BigInteger factorValue;
BigInteger tofactor;
bool prettyprint;
bool cunningham;
bool hexadecimal;
long long Gamma[386];
long long Delta[386];
long long AurifQ[386];
char tofactorDec[MAX_LEN*12];
extern int q[MAX_LEN];
int nbrToFactor[MAX_LEN];
struct sFactors astFactorsMod[5000];
int factorsMod[20000];
static void insertBigFactor(struct sFactors *pstFactors, const BigInteger *divisor,
  int type);
static char *findChar(char *str, char c);

static void GetYieldFrequency(void)
{
  yieldFreq = 1000000 / (NumberLength * NumberLength) + 1;
  if (yieldFreq > 100000)
  {
    yieldFreq = yieldFreq / 100000 * 100000;
  }
  else if (yieldFreq > 10000)
  {
    yieldFreq = yieldFreq / 10000 * 10000;
  }
  else if (yieldFreq > 1000)
  {
    yieldFreq = yieldFreq / 1000 * 1000;
  }
  else if (yieldFreq > 100)
  {
    yieldFreq = yieldFreq / 100 * 100;
  }
  else
  {      // Nothing to do.
  }
}

static int Cos(int N)
{
  int NMod8 = N % 8;
  if (NMod8 == 0)
  {
    return 1;
  }
  if (NMod8 == 4)
  {
    return -1;
  }
  return 0;
}

static int intTotient(int argument)
{
  int trialDivisor;
  int argumentDivisor = argument;
  int totient = argumentDivisor;
  if ((argumentDivisor % 2) == 0)
  {
    totient /= 2;
    do
    {
      argumentDivisor /= 2;
    } while ((argumentDivisor % 2) == 0);
  }
  if ((argumentDivisor % 3) == 0)
  {
    totient = totient * 2 / 3;
    do
    {
      argumentDivisor /= 3;
    } while ((argumentDivisor % 3) == 0);
  }
  trialDivisor = 5;
  while ((trialDivisor * trialDivisor) <= argumentDivisor)
  {
    if (((trialDivisor % 3) != 0) && ((argumentDivisor % trialDivisor) == 0))
    {
      totient = totient * (trialDivisor - 1) / trialDivisor;
      do
      {
        argumentDivisor /= trialDivisor;
      } while ((argumentDivisor % trialDivisor) == 0);
    }
    trialDivisor += 2;
  }
  if (argumentDivisor > 1)
  {
    totient = totient * (argumentDivisor - 1) / argumentDivisor;
  }
  return totient;
}

int Moebius(int argument)
{
  int moebius;
  int argumentDivisor;
  int trialDivisor;

  moebius = 1;
  argumentDivisor = argument;
  if ((argumentDivisor % 2) == 0)
  {
    moebius = -moebius;
    argumentDivisor /= 2;
    if ((argumentDivisor % 2) == 0)
    {
      return 0;
    }
  }
  if ((argumentDivisor % 3) == 0)
  {
    moebius = -moebius;
    argumentDivisor /= 3;
    if ((argumentDivisor % 3) == 0)
    {
      return 0;
    }
  }
  trialDivisor = 5;
  while ((trialDivisor * trialDivisor) <= argumentDivisor)
  {
    if ((trialDivisor % 3) != 0)
    {
      while ((argumentDivisor % trialDivisor) == 0)
      {
        moebius = -moebius;
        argumentDivisor /= trialDivisor;
        if ((argumentDivisor % trialDivisor) == 0)
        {
          return 0;
        }
      }
    }
    trialDivisor += 2;
  }
  if (argumentDivisor > 1)
  {
    moebius = -moebius;
  }
  return moebius;
}

void GetAurifeuilleFactor(struct sFactors *pstFactors, int L, const BigInteger *BigBase)
{
  static BigInteger x;
  static BigInteger Csal;
  static BigInteger Dsal;
  static BigInteger Nbr1;
  int k;

  (void)BigIntPowerIntExp(BigBase, L, &x);   // x <- BigBase^L.
  intToBigInteger(&Csal, 1);
  intToBigInteger(&Dsal, 1);
  for (k = 1; k < DegreeAurif; k++)
  {
    longToBigInteger(&Nbr1, Gamma[k]);
    (void)BigIntMultiply(&Csal, &x, &Csal);
    BigIntAdd(&Csal, &Nbr1, &Csal);      // Csal <- Csal * x + Gamma[k]
    longToBigInteger(&Nbr1, Delta[k]);
    (void)BigIntMultiply(&Dsal, &x, &Dsal);
    BigIntAdd(&Dsal, &Nbr1, &Dsal);      // Dsal <- Dsal * x + Gamma[k]
  }
  longToBigInteger(&Nbr1, Gamma[k]);
  (void)BigIntMultiply(&Csal, &x, &Csal);
  BigIntAdd(&Csal, &Nbr1, &Csal);        // Csal <- Csal * x + Gamma[k]
  (void)BigIntPowerIntExp(BigBase, (L + 1) / 2, &Nbr1);   // Nbr1 <- Dsal * base^((L+1)/2)
  (void)BigIntMultiply(&Dsal, &Nbr1, &Nbr1);
  BigIntAdd(&Csal, &Nbr1, &Dsal);
  insertBigFactor(pstFactors, &Dsal, TYP_AURIF);
  BigIntSubt(&Csal, &Nbr1, &Dsal);
  insertBigFactor(pstFactors, &Dsal, TYP_AURIF);
}

// Get Aurifeuille factors.
void InsertAurifFactors(struct sFactors *pstFactors, const BigInteger *BigBase,
  int exponent, int increment)
{
  int Incre = increment;
  int Expon = exponent;
  int Base = BigBase->limbs[0].x;
  if ((BigBase->nbrLimbs != 1) || (Base >= 386))
  {
    return;    // Base is very big, so go out.
  }
  if (((Expon % 2) == 0) && (Incre == -1))
  {
    do
    {
      Expon /= 2;
    } while ((Expon % 2) == 0);
    Incre = (Base % 4) - 2;
  }
  if (((Expon % Base) == 0)
    && (((Expon / Base) % 2) != 0)
    && ((((Base % 4) != 1) && (Incre == 1)) || (((Base % 4) == 1) && (Incre == -1))))
  {
    int N1;
    int q1;
    int L;
    int k;
    int N = Base;
    if ((N % 4) == 1)
    {
      N1 = N;
    }
    else
    {
      N1 = 2 * N;
    }
    DegreeAurif = intTotient(N1) / 2;
    for (k = 1; k <= DegreeAurif; k += 2)
    {
      AurifQ[k] = JacobiSymbol(N, k);
    }
    for (k = 2; k <= DegreeAurif; k += 2)
    {
      int t1 = k; // Calculate t2 = gcd(k, N1)
      int t2 = N1;
      while (t1 != 0)
      {
        int t3 = t2 % t1;
        t2 = t1;
        t1 = t3;
      }
      t1 = Moebius(N1 / t2) * intTotient(t2) * Cos((N - 1) * k);
      AurifQ[k] = t1;
    }
    Gamma[0] = 1;
    Delta[0] = 1;
    for (k = 1; k <= (DegreeAurif / 2); k++)
    {
      Gamma[k] = 0;
      Delta[k] = 0;
      for (int j = 0; j < k; j++)
      {
        int m = 2 * (k - j);
        Gamma[k] = Gamma[k] + (N * AurifQ[m - 1] * Delta[j]) - (AurifQ[m] * Gamma[j]);
        Delta[k] = Delta[k] + (AurifQ[m + 1] * Gamma[j]) - (AurifQ[m] * Delta[j]);
      }
      Gamma[k] /= 2 * k;
      Delta[k] = (Delta[k] + Gamma[k]) / ((2 * k) + 1);
    }
    for (k = (DegreeAurif / 2) + 1; k <= DegreeAurif; k++)
    {
      Gamma[k] = Gamma[DegreeAurif - k];
    }
    for (k = (DegreeAurif + 1) / 2; k < DegreeAurif; k++)
    {
      Delta[k] = Delta[DegreeAurif - k - 1];
    }
    q1 = Expon / Base;
    L = 1;
    while ((L * L) <= q1)
    {
      if ((q1 % L) == 0)
      {
        GetAurifeuilleFactor(pstFactors,L, BigBase);
        if (q1 != (L * L))
        {
          GetAurifeuilleFactor(pstFactors, q1 / L, BigBase);
        }
      }
      L += 2;
    }
  }
}

void copyString(const char *textFromServer)
{
  (void)strcpy(common.saveFactors.text, textFromServer);
}

static void Cunningham(struct sFactors *pstFactors, const BigInteger *BigBase, int Expon,
                int increment, const BigInteger *BigOriginal)
{
#ifdef __EMSCRIPTEN__
  char url[200];
  char *ptrUrl;
#endif
  int Expon2;
  int k;
  char *ptrFactorsAscii;
  static BigInteger Nbr1;
  static BigInteger Nbr2;

  common.saveFactors.text[0] = 0;    // Indicate no new factor found in advance.
  Expon2 = Expon;
  if (cunningham && (BigOriginal->nbrLimbs > 4))
  {   // Enter here on numbers of more than 40 digits if the user selected
      // get Cunningham factors from server.
#ifdef __EMSCRIPTEN__
    databack(lang ? "4<p>Obteniendo los factores primitivos conocidos del servidor Web.</p>":
                    "4<p>Requesting known primitive factors from Web server.</p>");
    // Format URL.
    ptrUrl = url;
    copyStr(&ptrUrl, "factors.pl?base=");
    int2dec(&ptrUrl, BigBase->limbs[0].x);
    copyStr(&ptrUrl, "&expon=");
    int2dec(&ptrUrl, Expon);
    copyStr(&ptrUrl, "&type=");
    *ptrUrl = ((increment > 0)? 'p': 'm');
    ptrUrl++;
    *ptrUrl = 0;
    getCunn(url, common.saveFactors.text);
#endif
  }
  ptrFactorsAscii = common.saveFactors.text;
  while (*ptrFactorsAscii > ' ')
  { // Loop through factors found in server.
    int nbrDigits;
    char *ptrEndFactor = findChar(ptrFactorsAscii, '*');
    if (ptrEndFactor == NULL)
    {
      nbrDigits = (int)strlen(ptrFactorsAscii);
      ptrEndFactor = ptrFactorsAscii + nbrDigits;
    }
    else
    {
      nbrDigits = (int)(ptrEndFactor - ptrFactorsAscii);
      ptrEndFactor++;
    }
    Dec2Bin(ptrFactorsAscii, Nbr1.limbs, nbrDigits, &Nbr1.nbrLimbs);
    Nbr1.sign = SIGN_POSITIVE;
    ptrFactorsAscii = ptrEndFactor;
    insertBigFactor(pstFactors, &Nbr1, TYP_TABLE);
  }
  while (((Expon2 % 2) == 0) && (increment == -1))
  {
    Expon2 /= 2;
    (void)BigIntPowerIntExp(BigBase, Expon2, &Nbr1);
    addbigint(&Nbr1, increment);
    insertBigFactor(pstFactors, &Nbr1, TYP_TABLE);
    InsertAurifFactors(pstFactors,BigBase, Expon2, 1);
  }
  k = 1;
  while ((k * k) <= Expon)
  {
    if ((Expon % k) == 0)
    {
      if ((k % 2) != 0)
      { /* Only for odd exponent */
        (void)BigIntPowerIntExp(BigBase, Expon / k, &Nbr1);
        addbigint(&Nbr1, increment);
        BigIntGcd(&Nbr1, BigOriginal, &Nbr2);   // Nbr2 <- gcd(Base^(Expon/k)+incre, original)
        insertBigFactor(pstFactors, &Nbr2, TYP_TABLE);
        CopyBigInt(&Temp1, BigOriginal);
        (void)BigIntDivide(&Temp1, &Nbr2, &Nbr1);
        insertBigFactor(pstFactors, &Nbr1, TYP_TABLE);
        InsertAurifFactors(pstFactors, BigBase, Expon / k, increment);
      }
      if (((Expon / k) % 2) != 0)
      { /* Only for odd exponent */
        (void)BigIntPowerIntExp(BigBase, k, &Nbr1);
        addbigint(&Nbr1, increment); 
        BigIntGcd(&Nbr1, BigOriginal, &Nbr2);   // Nbr2 <- gcd(Base^k+incre, original)
        insertBigFactor(pstFactors, &Nbr2, TYP_TABLE);
        CopyBigInt(&Temp1, BigOriginal);
        (void)BigIntDivide(&Temp1, &Nbr2, &Nbr1);
        insertBigFactor(pstFactors, &Nbr1, TYP_TABLE);
        InsertAurifFactors(pstFactors, BigBase, k, increment);
      }
    }
    k++;
  }
}

static bool ProcessExponent(struct sFactors *pstFactors, const BigInteger *numToFactor, int Exponent)
{
#ifdef __EMSCRIPTEN__
  char status[200];
  char *ptrStatus;
#endif
  static BigInteger NFp1;
  static BigInteger NFm1;
  static BigInteger nthRoot;
  static BigInteger rootN1;
  static BigInteger rootN;
  static BigInteger rootbak;
  static BigInteger nextroot;
  static BigInteger dif;
  double log2N;
#ifdef __EMSCRIPTEN__
  int elapsedTime = (int)(tenths() - originalTenthSecond);
  if ((elapsedTime / 10) != (oldTimeElapsed / 10))
  {
    oldTimeElapsed = elapsedTime;
    ptrStatus = status;
    copyStr(&ptrStatus, lang ? "4<p>Transcurrió " : "4<p>Time elapsed: ");
    GetDHMS(&ptrStatus, elapsedTime / 10);
    copyStr(&ptrStatus, lang ? "&nbsp;&nbsp;&nbsp;Exponente potencia +/- 1: " :
      "&nbsp;&nbsp;&nbsp;Power +/- 1 exponent: ");
    int2dec(&ptrStatus, Exponent);
    *ptrStatus = 0;                       // Add string terminator.
    databack(status);
}
#endif
  CopyBigInt(&NFp1, numToFactor);
  addbigint(&NFp1, 1);                    // NFp1 <- NumberToFactor + 1
  CopyBigInt(&NFm1, numToFactor);
  addbigint(&NFm1, -1);                   // NFm1 <- NumberToFactor - 1
  log2N = logBigNbr(&NFp1) / Exponent;    // Find nth root of number to factor.
  expBigNbr(&nthRoot, log2N);
  rootbak = nthRoot;
  for (;;)
  {
    (void)BigIntPowerIntExp(&nthRoot, Exponent - 1, &rootN1); // rootN1 <- nthRoot ^ (Exponent-1)
    (void)BigIntMultiply(&nthRoot, &rootN1, &rootN);  // rootN <- nthRoot ^ Exponent
    BigIntSubt(&NFp1, &rootN, &dif);            // dif <- NFp1 - rootN
    if (BigIntIsZero(&dif))
    { // Perfect power
      Cunningham(pstFactors, &nthRoot, Exponent, -1, numToFactor);
      return true;
    }
    addbigint(&dif, 1);                         // dif <- dif + 1
    (void)BigIntDivide(&dif, &rootN1, &Temp1);  // Temp1 <- dif / rootN1
    subtractdivide(&Temp1, 0, Exponent);        // Temp1 <- Temp1 / Exponent
    BigIntAdd(&Temp1, &nthRoot, &nextroot);     // nextroot <- Temp1 + nthRoot
    addbigint(&nextroot, -1);                   // nextroot <- nextroot - 1
    BigIntSubt(&nextroot, &nthRoot, &nthRoot);  // nthRoot <- nextroot - nthRoot
    if (nthRoot.sign == SIGN_POSITIVE)
    {
      break; // Not a perfect power
    }
    CopyBigInt(&nthRoot, &nextroot);
  }
  nthRoot = rootbak;
  for (;;)
  {
    (void)BigIntPowerIntExp(&nthRoot, Exponent - 1, &rootN1); // rootN1 <- nthRoot ^ (Exponent-1)
    (void)BigIntMultiply(&nthRoot, &rootN1, &rootN);     // rootN <- nthRoot ^ Exponent
    BigIntSubt(&NFm1, &rootN, &dif);            // dif <- NFm1 - rootN
    if (BigIntIsZero(&dif))
    { // Perfect power
      Cunningham(pstFactors, &nthRoot, Exponent, 1, numToFactor);
      return true;
    }
    addbigint(&dif, 1);                         // dif <- dif + 1
    (void)BigIntDivide(&dif, &rootN1, &Temp1);  // Temp1 <- dif / rootN1
    subtractdivide(&Temp1, 0, Exponent);        // Temp1 <- Temp1 / Exponent
    BigIntAdd(&Temp1, &nthRoot, &nextroot);     // nextroot <- Temp1 + nthRoot
    addbigint(&nextroot, -1);                   // nextroot <- nextroot - 1
    BigIntSubt(&nextroot, &nthRoot, &nthRoot);  // nthRoot <- nextroot - nthRoot
    if (nthRoot.sign == SIGN_POSITIVE)
    {
      break;                                    // Not a perfect power
    }
    CopyBigInt(&nthRoot, &nextroot);
  }
  return false;
}

static void PowerPM1Check(struct sFactors *pstFactors, const BigInteger *numToFactor)
{
  bool plus1 = false;
  bool minus1 = false;
  int Exponent = 0;
  int i;
  int j;
  int modulus;
  int mod9 = getRemainder(numToFactor, 9);
  int maxExpon = numToFactor->nbrLimbs * BITS_PER_GROUP;
  int numPrimes = (2 * maxExpon) + 3;
  double logar = logBigNbr(numToFactor);
  // 33219 = logarithm base 2 of max number supported = 10^10000.
  // Let n = a^b +/- 1 (n = number to factor).
  // If n!=1 or n!=0 or n!=7 (mod 8), then b cannot be even.
  modulus = numToFactor->limbs[0].x & 7;
  if ((modulus == 0) || (modulus == 1) || (modulus == 7))
  {       // b can be even
    (void)memset(common.ecm.ProcessExpon, 0xFF, sizeof(common.ecm.ProcessExpon));
  }
  else
  {       // b cannot be even
    (void)memset(common.ecm.ProcessExpon, 0xAA, sizeof(common.ecm.ProcessExpon));
  }
  (void)memset(common.ecm.primes, 0xFF, sizeof(common.ecm.primes));
  for (i = 2; (i * i) < numPrimes; i++)
  {       // Generation of primes using sieve of Eratosthenes.
    if ((common.ecm.primes[i >> 3] & (1 << (i & 7))) != 0)
    {     // Number i is prime.
      for (j = i * i; j < numPrimes; j += i)
      {   // Mark multiple of i as composite.
        common.ecm.primes[j >> 3] &= ~(1 << (j & 7));
      }
    }
  }
  // Let n = a^b +/- 1 (n = number to factor).
  // If -1<=n<=2 (mod p) does not hold, b cannot be multiple of p-1.
  // If -2<=n<=2 (mod p) does not hold, b cannot be multiple of (p-1)/2.
  for (i = 2; i < numPrimes; i++)
  {
    if ((common.ecm.primes[i>>3] & (1 << (i & 7))) != 0)
    {      // i is prime according to sieve.
           // If n+/-1 is multiple of p, then it must be multiple
           // of p^2, otherwise it cannot be a perfect power.
      uint64_t remainder;
      int index;
      int minRemainder;
      int maxRemainder;
      int rem = getRemainder(numToFactor, i);
      longToBigInteger(&Temp1, (uint64_t)i*(uint64_t)i);
      (void)BigIntRemainder(numToFactor, &Temp1, &Temp2);     // Temp2 <- nbrToFactor % (i*i)
      remainder = (uint64_t)Temp2.limbs[0].x;
      if ((rem == 1) || (rem == (i - 1)))
      {
        if (Temp2.nbrLimbs > 1)
        {
          remainder += (uint64_t)Temp2.limbs[1].x << BITS_PER_GROUP;
        }
        // NumberFactor cannot be a power + 1 if condition holds.
        plus1 = ((rem == 1) && (remainder != 1));
        // NumberFactor cannot be a power - 1 if condition holds.
        minus1 = ((rem == (i - 1)) && (remainder != (((uint64_t)i*(uint64_t)i) - 1)));
      }
      index = i / 2;
      if (!(common.ecm.ProcessExpon[index >> 3] & (1<<(index&7))))
      {
        continue;
      }
      modulus = remainder % i;
      if (plus1)
      {
        minRemainder = 1;
      }
      else
      {
        minRemainder = 2;
      }
      if (minus1)
      {
        maxRemainder = i - 1;
      }
      else
      {
        maxRemainder = i - 2;
      }
      if ((modulus > minRemainder) && (modulus < maxRemainder))
      {
        for (j = index; j <= maxExpon; j += index)
        {
          common.ecm.ProcessExpon[j >> 3] &= ~(1 << (j&7));
        }
      }
      else
      {
        if (modulus == (i - 2))
        {
          for (j = i - 1; j <= maxExpon; j += i - 1)
          {
            common.ecm.ProcessExpon[j >> 3] &= ~(1 << (j & 7));
          }
        }
      }
    }
  }
  for (j = 2; j < 100; j++)
  {
    double u = logar / log(j) + .000005;
    Exponent = (int)floor(u);
    if ((u - Exponent) > .00001)
    {
      continue;
    }
    if (((Exponent % 3) == 0) && (mod9 > 2) && (mod9 < 7))
    {
      continue;
    }
    if (!(common.ecm.ProcessExpon[Exponent >> 3] & (1 << (Exponent & 7))))
    {
      continue;
    }
    if (ProcessExponent(pstFactors, numToFactor, Exponent))
    {
      return;
    }
  }
  for (; Exponent >= 2; Exponent--)
  {
    if (((Exponent % 3) == 0) && (mod9 > 2) && (mod9 < 7))
    {
      continue;
    }
    if (!(common.ecm.ProcessExpon[Exponent >> 3] & (1 << (Exponent & 7))))
    {
      continue;
    }
    if (ProcessExponent(pstFactors, numToFactor, Exponent))
    {
      return;
    }
  }
}

// Perform Lehman algorithm
static void Lehman(const BigInteger *nbr, int k, BigInteger *factor)
{
  const int bitsSqrLow[] =
  {
    0x00000003, // 3
    0x00000013, // 5
    0x00000017, // 7
    0x0000023B, // 11
    0x0000161B, // 13
    0x0001A317, // 17
    0x00030AF3, // 19
    0x0005335F, // 23
    0x13D122F3, // 29
    0x121D47B7, // 31
    0x5E211E9B, // 37
    0x82B50737, // 41
    0x83A3EE53, // 43
    0x1B2753DF, // 47
    0x3303AED3, // 53
    0x3E7B92BB, // 59
    0x0A59F23B, // 61
  };
  const int bitsSqrHigh[] =
  {
    0x00000000, // 3
    0x00000000, // 5
    0x00000000, // 7
    0x00000000, // 11
    0x00000000, // 13
    0x00000000, // 17
    0x00000000, // 19
    0x00000000, // 23
    0x00000000, // 29
    0x00000000, // 31
    0x00000016, // 37
    0x000001B3, // 41
    0x00000358, // 43
    0x00000435, // 47
    0x0012DD70, // 53
    0x022B6218, // 59
    0x1713E694, // 61
  };
  const int primes[] = { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61 };
  int nbrs[17];
  int diffs[17];
  int i;
  int m;
  int r;
  int nbrLimbs;
  int nbrIterations;
  static BigInteger sqrRoot;
  static BigInteger nextroot;
  static BigInteger a;
  static BigInteger c;
  static BigInteger sqr;
  static BigInteger val;
  if ((nbr->limbs[0].x & 1) == 0)
  { // nbr Even
    r = 0;
    m = 1;
  }
  else
  {
    if ((k % 2) == 0)
    { // k Even
      r = 1;
      m = 2;
    }
    else
    { // k Odd
      r = (k + nbr->limbs[0].x) & 3;
      m = 4;
    }
  }
  intToBigInteger(&sqr, k<<2);
  (void)BigIntMultiply(&sqr, nbr, &sqr);
  squareRoot(sqr.limbs, sqrRoot.limbs, sqr.nbrLimbs, &sqrRoot.nbrLimbs);
  sqrRoot.sign = SIGN_POSITIVE;
  CopyBigInt(&a, &sqrRoot);
  for (;;)
  {
    if ((a.limbs[0].x & (m-1)) == r)
    {
      (void)BigIntMultiply(&a, &a, &nextroot);
      BigIntSubt(&nextroot, &sqr, &nextroot);
      if (nextroot.sign == SIGN_POSITIVE)
      {
        break;
      }
    }
    addbigint(&a, 1);                         // a <- a + 1
  }
  (void)BigIntMultiply(&a, &a, &nextroot);
  BigIntSubt(&nextroot, &sqr, &c);
  for (i = 0; i < 17; i++)
  {
    int pr = primes[i];
    nbrs[i] = getRemainder(&c, pr);    // nbrs[i] <- c % primes[i]
    diffs[i] = m * ((getRemainder(&a, pr) * 2) + m) % pr;
  }
  nbrLimbs = factor->nbrLimbs;
  if (nbrLimbs > 10)
  {
    nbrIterations = 10000;
  }
  else
  {
    nbrIterations = 1000 * factor->nbrLimbs;
  }
  for (int j = 0; j < nbrIterations; j++)
  {
    for (i = 0; i < 17; i++)
    {
      int shiftBits = nbrs[i];
      if (shiftBits < 32)
      {
        if ((bitsSqrLow[i] & (1 << shiftBits)) == 0)
        { // Not a perfect square
          break;
        }
      }
      else if ((bitsSqrHigh[i] & (1 << (shiftBits-32))) == 0)
      { // Not a perfect square
        break;
      }
    }
    if (i == 17)
    { // Test for perfect square
      intToBigInteger(&c, m * j);           // c <- m * j
      BigIntAdd(&a, &c, &val);
      (void)BigIntMultiply(&val, &val, &c);       // c <- val * val
      BigIntSubt(&c, &sqr, &c);             // c <- val * val - sqr
      squareRoot(c.limbs, sqrRoot.limbs, c.nbrLimbs, &sqrRoot.nbrLimbs);
      sqrRoot.sign = SIGN_POSITIVE;         // sqrRoot <- sqrt(c)
      BigIntAdd(&sqrRoot, &val, &sqrRoot);
      BigIntGcd(&sqrRoot, nbr, &c);         // Get GCD(sqrRoot + val, nbr)
      if (c.nbrLimbs > 1)
      {    // Non-trivial factor has been found.
        CopyBigInt(factor, &c);
        return;
      }
    }
    for (i = 0; i < 17; i++)
    {
      nbrs[i] = (nbrs[i] + diffs[i]) % primes[i];
      diffs[i] = (diffs[i] + (2 * m * m)) % primes[i];
    }
  }
  intToBigInteger(factor, 1);   // Factor not found.
}

#ifdef __EMSCRIPTEN__
char *ShowFactoredPart(const BigInteger *pNbr, const struct sFactors* pstFactors)
{
  ptrLowerText = lowerText;
  *ptrLowerText = '3';
  ptrLowerText++;
  if ((pstFactors != NULL) && (pstFactors->multiplicity > 1))
  {    // Some factorization known.
    int NumberLengthBak = NumberLength;
    copyStr(&ptrLowerText, "<p class=\"blue\">");
    SendFactorizationToOutput(pstFactors, &ptrLowerText, true);
    copyStr(&ptrLowerText, "</p>");
    NumberLength = NumberLengthBak;
  }
  if (StepECM == 3)
  {
    copyStr(&ptrLowerText, lang ? "<p>Comprobando si es primo " : "<p>Testing primality of ");
  }
  else
  {
    copyStr(&ptrLowerText, lang ? "<p>Factorizando " : "<p>Factoring ");
  }
  if (hexadecimal)
  {
    Bin2Hex(&ptrLowerText, pNbr->limbs, pNbr->nbrLimbs, groupLen);
  }
  else
  {
    Bin2Dec(&ptrLowerText, pNbr->limbs, pNbr->nbrLimbs, groupLen);
  }
  copyStr(&ptrLowerText, "</p>");
  return ptrLowerText;
}

void ShowLowerText(void)
{
  databack(lowerText);
}
#endif

static bool isOne(const limb* nbr, int length)
{
  if (nbr->x != 1)
  {
    return false;
  }
  for (int ctr = 1; ctr < length; ctr++)
  {
    if ((nbr + ctr)->x != 0)
    {
      return false;
    }
  }
  return true;
}

static void performFactorization(const BigInteger *numToFactor, const struct sFactors *pstFactors)
{
#ifndef __EMSCRIPTEN__
  (void)pstFactors;     // Ignore parameter.
#endif
  static BigInteger potentialFactor;
  common.ecm.fieldTX = common.ecm.TX;
  common.ecm.fieldTZ = common.ecm.TZ;
  common.ecm.fieldUX = common.ecm.UX;
  common.ecm.fieldUZ = common.ecm.UZ;
  //  int Prob
  //  BigInteger NN

  common.ecm.fieldAA = common.ecm.AA;
  NumberLength = numToFactor->nbrLimbs;
  (void)memcpy(TestNbr, numToFactor->limbs, NumberLength * sizeof(limb));
  GetYieldFrequency();
  GetMontgomeryParms(NumberLength);
  (void)memset(common.ecm.M, 0, NumberLength * sizeof(limb));
  (void)memset(common.ecm.DX, 0, NumberLength * sizeof(limb));
  (void)memset(common.ecm.DZ, 0, NumberLength * sizeof(limb));
  (void)memset(common.ecm.W3, 0, NumberLength * sizeof(limb));
  (void)memset(common.ecm.W4, 0, NumberLength * sizeof(limb));
  (void)memset(common.ecm.GD, 0, NumberLength * sizeof(limb));
#ifdef __EMSCRIPTEN__
  ptrLowerText = ShowFactoredPart(numToFactor, pstFactors);
#endif
  EC--;
  foundByLehman = false;
  do
  {
    enum eEcmResult ecmResp;
    // Try to factor BigInteger N using Lehman algorithm. Result in potentialFactor.
    Lehman(numToFactor, EC % 50000000, &potentialFactor);
    if (potentialFactor.nbrLimbs > 1)
    {                // Factor found.
      (void)memcpy(common.ecm.GD, potentialFactor.limbs, NumberLength * sizeof(limb));
      (void)memset(&common.ecm.GD[potentialFactor.nbrLimbs], 0,
        (NumberLength - potentialFactor.nbrLimbs) * sizeof(limb));
      foundByLehman = true;
      break;
    }
    ecmResp = ecmCurve(&EC, &NextEC);
    if (ecmResp == CHANGE_TO_SIQS)
    {    // Perform SIQS
#ifdef __EMSCRIPTEN__
      double originalTenths = tenths();
      uint64_t oldModularMult = lModularMult;
#endif
      FactoringSIQS(TestNbr, common.ecm.GD);
#ifdef __EMSCRIPTEN__
      SIQSModMult += lModularMult - oldModularMult;
      timeSIQS += (int)(tenths() - originalTenths);
      nbrSIQS++;
#endif

      break;
    }
    else if (ecmResp == FACTOR_FOUND)
    {
      break;
    }
  } while (!memcmp(common.ecm.GD, TestNbr, NumberLength*sizeof(limb)) ||
           isOne(common.ecm.GD, NumberLength));
  StepECM = 0; /* do not show pass number on screen */
}

void SendFactorizationToOutput(const struct sFactors *pstFactors, char **pptrOutput, bool doFactorization)
{
  char *ptrOutput = *pptrOutput;
  copyStr(&ptrOutput, tofactorDec);
  if (doFactorization)
  {
    const struct sFactors *pstFactor;
    pstFactor = pstFactors+1;
    if ((tofactor.sign == SIGN_POSITIVE) && (pstFactors->multiplicity == 1) && (pstFactor->multiplicity == 1) &&
      ((*pstFactor->ptrFactor > 1) || (*(pstFactor->ptrFactor + 1) > 1)))
    {    // Do not show zero or one as prime.
      copyStr(&ptrOutput, lang ? " es primo" : " is prime");
    }
    else
    {
      int i = 0;
      copyStr(&ptrOutput, " = ");
      if (tofactor.sign == SIGN_NEGATIVE)
      {
        *ptrOutput = '-';
        ptrOutput++;
        if ((tofactor.nbrLimbs > 1) || (tofactor.limbs[0].x > 1))
        {
          if (prettyprint)
          {
            copyStr(&ptrOutput, "1 &times; ");
          }
          else
          {
            copyStr(&ptrOutput, "1 * ");
          }
        }
      }
      for (;;)
      {
        NumberLength = *pstFactor->ptrFactor;
        IntArray2BigInteger(pstFactor->ptrFactor, &factorValue);
        if (hexadecimal)
        {
          Bin2Hex(&ptrOutput, factorValue.limbs, factorValue.nbrLimbs, groupLen);
        }
        else
        {
          Bin2Dec(&ptrOutput, factorValue.limbs, factorValue.nbrLimbs, groupLen);
        }
        if (pstFactor->multiplicity > 1)
        {
          if (prettyprint)
          {
            copyStr(&ptrOutput, "<sup>");
            int2dec(&ptrOutput, pstFactor->multiplicity);
            copyStr(&ptrOutput, "</sup>");
          }
          else
          {
            *ptrOutput = '^';
            ptrOutput++;
            int2dec(&ptrOutput, pstFactor->multiplicity);
          }
        }
#ifdef ENABLE_VERBOSE
        int type = pstFactor->type;
        bool isPrime = (pstFactor->upperBound == 0)? true: false;
        if (type > 0)
        {
          int compositeType = type / 50000000 * 50000000;
          copyStr(&ptrOutput, " <span class=\"verbose\">(");
          if (compositeType == TYP_AURIF)
          {
            copyStr(&ptrOutput, "Aurifeuille");
            if (!isPrime)
            {
              copyStr(&ptrOutput, lang ? " - Compuesto" : " - Composite");
            }
          }
          else if (compositeType == TYP_TABLE)
          {
            copyStr(&ptrOutput, lang ? "Tabla" : "Table");
            if (!isPrime)
            {
              copyStr(&ptrOutput, lang ? " - Compuesto" : " - Composite");
            }
          }
          else if (compositeType == TYP_SIQS)
          {
            copyStr(&ptrOutput, lang? "<abbr title=\"Criba cuadrática autoinicializada\">SIQS</abbr>":
                                "<abbr title=\"Self-Initializing Quadratic Sieve\">SIQS</abbr>");
            if (!isPrime)
            {
              copyStr(&ptrOutput, lang ? " - Compuesto" : " - Composite");
            }
          }
          else if (compositeType == TYP_LEHMAN)
          {
            copyStr(&ptrOutput, "Lehman");
            if (!isPrime)
            {
              copyStr(&ptrOutput, lang ? " - Compuesto" : " - Composite");
            }
          }
          else if (compositeType == TYP_RABIN)
          {
            copyStr(&ptrOutput, lang ? "Miller y Rabin" : "Miller &amp; Rabin");
            if (!isPrime)
            {
              copyStr(&ptrOutput, lang ? " - Compuesto" : " - Composite");
            }
          }
          else if (compositeType == TYP_DIVISION)
          {
            copyStr(&ptrOutput, lang ? "División" : "Division");
            if (!isPrime)
            {
              copyStr(&ptrOutput, lang ? " - Compuesto" : " - Composite");
            }
          }
          else if (type > TYP_EC)
          {
            if (isPrime)
            {
              copyStr(&ptrOutput, lang ? "<abbr title=\"Método de curvas elípticas\">ECM</abbr>, curva " :
                "<abbr title=\"Elliptic curve method\">ECM</abbr>, curve ");
              int2dec(&ptrOutput, type - TYP_EC);
              *ptrOutput = 0;     // Add string terminator.
            }
            else
            {
              copyStr(&ptrOutput, lang ? "Compuesto" : "Composite");
            }
          }
          else if (!isPrime)
          {
            copyStr(&ptrOutput, lang ? "Compuesto": "Composite");
          }
          copyStr(&ptrOutput, ")</span>");
        }
        else if (!isPrime)
        {
          copyStr(&ptrOutput, "<span class=\"terse\"> (");
          copyStr(&ptrOutput, lang ? "Compuesto" : "Composite");
          copyStr(&ptrOutput, ")</span>");
        }
        if (type < 0)
        {
          copyStr(&ptrOutput, " (Unknown)");
        }
#endif
        i++;
        if (i == pstFactors->multiplicity)
        {
          break;
        }
        if (prettyprint)
        {
          copyStr(&ptrOutput, " &times; ");
        }
        else
        {
          copyStr(&ptrOutput, " * ");
        }
        pstFactor++;
      }
    }
  }
  *pptrOutput = ptrOutput;
}

static void SortFactors(struct sFactors *pstFactors)
{
  int factorNumber;
  int ctr;
  struct sFactors *pstCurFactor = pstFactors + 1;
  struct sFactors stTempFactor;
  int *ptrNewFactor;
  struct sFactors *pstNewFactor;
  for (factorNumber = 1; factorNumber <= pstFactors->multiplicity; factorNumber++)
  {
    pstNewFactor = pstCurFactor + 1;
    for (int factorNumber2 = factorNumber + 1; factorNumber2 <= pstFactors->multiplicity; factorNumber2++)
    {
      const int *ptrFactor = pstCurFactor->ptrFactor;
      const int *ptrFactor2 = pstNewFactor->ptrFactor;
      if (*ptrFactor < *ptrFactor2)
      {     // Factors already in correct order.
        pstNewFactor++;
        continue;
      }
      if (*ptrFactor == *ptrFactor2)
      {
        for (ctr = *ptrFactor; ctr > 1; ctr--)
        {
          if (*(ptrFactor + ctr) != *(ptrFactor2 + ctr))
          {
            break;
          }
        }
        if (*(ptrFactor + ctr) < *(ptrFactor2 + ctr))
        {     // Factors already in correct order.
          pstNewFactor++;
          continue;
        }
        if (*(ptrFactor + ctr) == *(ptrFactor2 + ctr))
        {     // Factors are the same.
          pstCurFactor->multiplicity += pstNewFactor->multiplicity;
          ctr = pstFactors->multiplicity - factorNumber2;
          if (ctr > 0)
          {
            (void)memmove(pstNewFactor, pstNewFactor + 1, ctr * sizeof(struct sFactors));
          }
          pstFactors->multiplicity--;   // Indicate one less known factor.
          pstNewFactor++;
          continue;
        }
      }
      // Exchange both factors.
      (void)memcpy(&stTempFactor, pstCurFactor, sizeof(struct sFactors));
      (void)memcpy(pstCurFactor, pstNewFactor, sizeof(struct sFactors));
      (void)memcpy(pstNewFactor, &stTempFactor, sizeof(struct sFactors));
      pstNewFactor++;
    }
    pstCurFactor++;
  }
  // Find location for new factors.
  ptrNewFactor = 0;
  pstCurFactor = pstFactors + 1;
  for (factorNumber = 1; factorNumber <= pstFactors->multiplicity; factorNumber++)
  {
    int *ptrPotentialNewFactor = pstCurFactor->ptrFactor + *(pstCurFactor->ptrFactor) + 1;
    if (ptrPotentialNewFactor > ptrNewFactor)
    {
      ptrNewFactor = ptrPotentialNewFactor;
    }
    pstCurFactor++;
  }
  pstFactors->ptrFactor = ptrNewFactor;
}

// Insert new factor found into factor array. This factor array must be sorted.
static void insertIntFactor(struct sFactors *pstFactors, struct sFactors *pstFactorDividend,
  int divisor, int expon, const BigInteger *cofactor)
{
  struct sFactors *pstCurFactor;
  int multiplicity;
  int factorNumber;
  int *ptrFactor = pstFactorDividend->ptrFactor;
  int nbrLimbs = *ptrFactor;
  int *ptrValue;
  pstFactorDividend->upperBound = divisor;
  // Divide number by factor just found.
  if (cofactor == NULL)
  {        // Find cofactor.
    DivBigNbrByInt(ptrFactor + 1, divisor, ptrFactor + 1, nbrLimbs);
    if (*(ptrFactor + nbrLimbs) == 0)
    {
      (*ptrFactor)--;
    }
  }
  else
  {        // Cofactor given as a parameter.
    NumberLength = cofactor->nbrLimbs;
    BigInteger2IntArray(ptrFactor, cofactor);
  }
  // Check whether prime is already in factor list.
  pstCurFactor = pstFactors+1;
  for (factorNumber = 1; factorNumber <= pstFactors->multiplicity; factorNumber++, pstCurFactor++)
  {
    ptrValue = pstCurFactor->ptrFactor;  // Point to factor in factor array.
    if ((*ptrValue == 1) && (*(ptrValue+1) == divisor))
    {  // Prime already found: increment multiplicity and go out.
      pstCurFactor->multiplicity += pstFactorDividend->multiplicity * expon;
      ptrValue = pstFactorDividend->ptrFactor;
      if ((*ptrValue == 1) && (*(ptrValue + 1) == 1))
      {    // Dividend is 1 now, so discard it.
        *pstFactorDividend = *(pstFactors + pstFactors->multiplicity);
        pstFactors->multiplicity--;
      }
      SortFactors(pstFactors);
      return;
    }
    if ((*ptrValue > 1) || (*(ptrValue + 1) > divisor))
    {   // Factor in factor list is greater than factor to insert. Exit loop.
      break;
    }
  }
  pstFactors->multiplicity++; // Indicate new known factor.
  multiplicity = pstFactorDividend->multiplicity;
  // Move all elements.
  ptrValue = pstFactorDividend->ptrFactor;
  if (pstFactors->multiplicity > factorNumber)
  {
    (void)memmove(pstCurFactor + 1, pstCurFactor,
      (pstFactors->multiplicity - factorNumber) * sizeof(struct sFactors));
  }
  if ((*ptrValue == 1) && (*(ptrValue + 1) == 1))
  {
    pstCurFactor = pstFactorDividend;
  }
  else
  {
    ptrValue = pstFactors->ptrFactor;
    pstCurFactor->ptrFactor = ptrValue;
    pstCurFactor->multiplicity = multiplicity * expon;
    pstFactors->ptrFactor += 2;  // Next free memory.
  }
  pstCurFactor->upperBound = 0;
  *ptrValue = 1;  // Number of limbs.
  *(ptrValue + 1) = divisor;
  SortFactors(pstFactors);
}

// Insert new factor found into factor array. This factor array must be sorted.
// The divisor must be also sorted.
static void insertBigFactor(struct sFactors *pstFactors, const BigInteger *divisor, int type)
{
  int typeFactor = type;
  struct sFactors *pstCurFactor;
  int lastFactorNumber = pstFactors->multiplicity;
  struct sFactors *pstNewFactor = pstFactors + lastFactorNumber + 1;
  int *ptrNewFactorLimbs = pstFactors->ptrFactor;
  pstCurFactor = pstFactors + 1;
  for (int factorNumber = 1; factorNumber <= lastFactorNumber; factorNumber++)
  {     // For each known factor...
    int *ptrFactor = pstCurFactor->ptrFactor;
    NumberLength = *ptrFactor;
    IntArray2BigInteger(ptrFactor, &Temp2);    // Convert known factor to Big Integer.
    BigIntGcd(divisor, &Temp2, &Temp3);         // Temp3 is the GCD between known factor and divisor.
    if ((Temp3.nbrLimbs == 1) && (Temp3.limbs[0].x < 2))
    {                                           // divisor is not a new factor (GCD = 0 or 1).
      pstCurFactor++;
      continue;
    }
    if (TestBigNbrEqual(&Temp2, &Temp3))
    {                                           // GCD is equal to known factor.
      pstCurFactor++;
      continue;
    }
    // At this moment both GCD and known factor / GCD are new known factors. Replace the known factor by
    // known factor / GCD and generate a new known factor entry.
    NumberLength = Temp3.nbrLimbs;
    BigInteger2IntArray(ptrNewFactorLimbs, &Temp3);      // Append new known factor.
    (void)BigIntDivide(&Temp2, &Temp3, &Temp4);          // Divide by this factor.
    NumberLength = Temp4.nbrLimbs;
    BigInteger2IntArray(ptrFactor, &Temp4);              // Overwrite old known factor.
    pstNewFactor->multiplicity = pstCurFactor->multiplicity;
    pstNewFactor->ptrFactor = ptrNewFactorLimbs;
    pstNewFactor->upperBound = pstCurFactor->upperBound;
    if (typeFactor < 50000000)
    {          // Factor found using ECM.
      pstNewFactor->type = TYP_EC + EC;
      typeFactor = pstCurFactor->type / 50000000 * 50000000;
      if (typeFactor == 0)
      {
        pstCurFactor->type = TYP_DIVISION + EC;
      }
      else
      {
        pstCurFactor->type = typeFactor + EC;
      }
    }
    else
    {          // Found otherwise.
      pstNewFactor->type = typeFactor;
    }
    pstNewFactor++;
    pstFactors->multiplicity++;
    ptrNewFactorLimbs += 1 + Temp3.nbrLimbs;
    pstCurFactor++;
  }
  // Sort factors in ascending order. If two factors are equal, coalesce them.
  // Divide number by factor just found.
  SortFactors(pstFactors);
}
#ifdef __EMSCRIPTEN__
static void showECMStatus(void)
{
  char status[200];
  int elapsedTime;
  char *ptrStatus;
  if ((lModularMult % yieldFreq) != 0)
  {
    return;
  }
  elapsedTime = (int)(tenths() - originalTenthSecond);
  if ((elapsedTime / 10) == (oldTimeElapsed / 10))
  {
    return;
  }
  oldTimeElapsed = elapsedTime;
  ptrStatus = status;
  copyStr(&ptrStatus, lang ? "4<p>Transcurrió " : "4<p>Time elapsed: ");
  GetDHMS(&ptrStatus, elapsedTime / 10);
  copyStr(&ptrStatus, "&nbsp;&nbsp;&nbsp;");  // Separate with three spaces.
  switch (StepECM)
  {
  case 1:
    copyStr(&ptrStatus, lang ? "Paso 1: " : "Step 1: ");
    int2dec(&ptrStatus, indexPrimes / (nbrPrimes / 100));
    *ptrStatus = '%';
    ptrStatus++;
    break;
  case 2:
    copyStr(&ptrStatus, lang ? "Paso 2: " : "Step 2: ");
    int2dec(&ptrStatus, (maxIndexM == 0)? 0 : (indexM / (maxIndexM / 100)));
    *ptrStatus = '%';
    ptrStatus++;
    break;
  case 3:
    copyStr(&ptrStatus, lang ? "Progreso: " : "Progress: ");
    int2dec(&ptrStatus, percentageBPSW);
    *ptrStatus = '%';
    ptrStatus++;
    break;
  }
  copyStr(&ptrStatus, "</p>");
  databack(status);
}

// Save factors in Web Storage so factorization can continue the next time the application runs.
static void SaveFactors(struct sFactors *pstFactors)
{
#ifdef FACTORIZATION_APP
  struct sFactors *pstCurFactor = pstFactors + 1;
  int factorNbr;
  BigInteger bigint;
  char *ptrText = common.saveFactors.text;
  if (oldNbrFactors == pstFactors->multiplicity)
  {
    return;
  }
  oldNbrFactors = pstFactors->multiplicity;
  *ptrText = '8';
  ptrText++;
  copyStr(&ptrText, ptrInputText);
  *ptrText = '=';
  ptrText++;
  for (factorNbr = 1; factorNbr <= pstFactors->multiplicity; factorNbr++)
  {
    if (factorNbr > 1)
    {
      *ptrText = '*';
      ptrText++;
    }
    NumberLength = *pstCurFactor->ptrFactor;
    IntArray2BigInteger(pstCurFactor->ptrFactor, &bigint);
    BigInteger2Dec(&ptrText, &bigint, -100000);   // Factors are saved in decimal.
    *ptrText = '^';
    ptrText++;
    int2dec(&ptrText, pstCurFactor->multiplicity);
    *ptrText = '(';
    ptrText++;
    int2dec(&ptrText, pstCurFactor->upperBound);
    *ptrText = ',';
    ptrText++;
    int2dec(&ptrText, pstCurFactor->type);
    *ptrText = ')';
    ptrText++;
    pstCurFactor++;
  }
  *ptrText = 0;
  databack(common.saveFactors.text);
#endif
}

#endif

static char *findChar(char *str, char c)
{
  char* ptrStr = str;
  while (*ptrStr != '\0')
  {
    if (*ptrStr == c)
    {
      return ptrStr;
    }
    ptrStr++;
  }
  return NULL;
}

static bool getNextInteger(char **ppcFactors, int *result, char delimiter)
{
  char *pcFactors = *ppcFactors;
  char *ptrCharFound = findChar(pcFactors, delimiter);
  if (ptrCharFound == NULL)
  {
    return true;
  }
  *ptrCharFound = 0;
  Dec2Bin(pcFactors, Temp1.limbs, (int)(ptrCharFound - pcFactors), &Temp1.nbrLimbs);
  if (Temp1.nbrLimbs != 1)
  {   // Exponent is too large.
    return true;
  }
  *result = Temp1.limbs[0].x;
  *ppcFactors = ptrCharFound+1;
  return false;
}

// Return: 0 = No factors found.
//         1 = Factors found.
// Use: Xaux for square root of -1.
//      Zaux for square root of 1.
static int factorCarmichael(BigInteger *pValue, struct sFactors *pstFactors)
{
  int randomBase = 0;
#ifdef __EMSCRIPTEN__
  timePrimalityTests = 0;  // Reset to zero all timings.
  timeSIQS = 0;
  timeECM = 0;
  nbrPrimalityTests = 0;   // Reset to zero all counters.
  nbrSIQS = 0;
  nbrECM = 0;
  SIQSModMult = 0;
#endif
  bool factorsFound = false;
  int nbrLimbsQ;
  int ctr;
  int nbrLimbs = pValue->nbrLimbs;
  bool sqrtOneFound = false;
  bool sqrtMinusOneFound = false;
  int Aux1Len;
  limb *pValueLimbs = pValue->limbs;
  (pValueLimbs + nbrLimbs)->x = 0;
  (void)memcpy(q, pValueLimbs, (nbrLimbs + 1) * sizeof(limb));
  nbrLimbsQ = nbrLimbs;
  q[0]--;                     // q = p - 1 (p is odd, so there is no carry).
  (void)memcpy(common.ecm.Aux1, q, (nbrLimbsQ + 1) * sizeof(q[0]));
  Aux1Len = nbrLimbs;
  DivideBigNbrByMaxPowerOf2(&ctr, common.ecm.Aux1, &Aux1Len);
  (void)memcpy(TestNbr, pValueLimbs, nbrLimbs * sizeof(limb));
  TestNbr[nbrLimbs].x = 0;
  GetMontgomeryParms(nbrLimbs);
  for (int countdown = 20; countdown > 0; countdown--)
  {
    NumberLength = nbrLimbs;
    randomBase = (int)(((uint64_t)randomBase * 89547121) + 1762281733) & MAX_INT_NBR;
    modPowBaseInt(randomBase, common.ecm.Aux1, Aux1Len, common.ecm.Aux2); // Aux2 = base^Aux1.
                                                 // If Mult1 = 1 or Mult1 = TestNbr-1, then try next base.
    if (checkOne(common.ecm.Aux2, nbrLimbs) || checkMinusOne(common.ecm.Aux2, nbrLimbs))
    {
      continue;    // This base cannot find a factor. Try another one.
    }
    for (int i = 0; i < ctr; i++)
    {              // Loop that squares number.
      modmult(common.ecm.Aux2, common.ecm.Aux2, common.ecm.Aux3);
      if (checkOne(common.ecm.Aux3, nbrLimbs) != 0)
      {            // Non-trivial square root of 1 found.
        if (!sqrtOneFound)
        {          // Save it to perform GCD later.
          (void)memcpy(common.ecm.Zaux, common.ecm.Aux2, nbrLimbs*sizeof(limb));
          sqrtOneFound = true;
        }
        else
        {          // Try to find non-trivial factor by doing GCD.
          SubtBigNbrMod(common.ecm.Aux2, common.ecm.Zaux, common.ecm.Aux4);
          UncompressLimbsBigInteger(common.ecm.Aux4, &Temp2);
          BigIntGcd(pValue, &Temp2, &Temp4);
          if (((Temp4.nbrLimbs != 1) || (Temp4.limbs[0].x > 1)) &&
            ((Temp4.nbrLimbs != NumberLength) ||
              memcmp(pValue->limbs, Temp4.limbs, NumberLength * sizeof(limb))))
          {          // Non-trivial factor found.
            insertBigFactor(pstFactors, &Temp4, TYP_RABIN);
            factorsFound = true;
          }
        }
                   // Try to find non-trivial factor by doing GCD.
        NumberLength = nbrLimbs;
        AddBigNbrMod(common.ecm.Aux2, MontgomeryMultR1, common.ecm.Aux4);
        UncompressLimbsBigInteger(common.ecm.Aux4, &Temp2);
        BigIntGcd(pValue, &Temp2, &Temp4);
        if (((Temp4.nbrLimbs != 1) || (Temp4.limbs[0].x > 1)) &&
          ((Temp4.nbrLimbs != NumberLength) ||
            memcmp(pValue->limbs, Temp4.limbs, NumberLength * sizeof(limb))))
        {          // Non-trivial factor found.
          insertBigFactor(pstFactors, &Temp4, TYP_RABIN);
          factorsFound = true;
        }
        i = ctr;
        continue;  // Find more factors.
      }
      if (checkMinusOne(common.ecm.Aux3, nbrLimbs) != 0)
      {            // Square root of 1 found.
        if (!sqrtMinusOneFound)
        {          // Save it to perform GCD later.
          (void)memcpy(common.ecm.Xaux, common.ecm.Aux2, nbrLimbs * sizeof(limb));
          sqrtOneFound = true;
        }
        else
        {          // Try to find non-trivial factor by doing GCD.
          SubtBigNbrMod(common.ecm.Aux3, common.ecm.Xaux, common.ecm.Aux4);
          UncompressLimbsBigInteger(common.ecm.Aux4, &Temp2);
          BigIntGcd(pValue, &Temp2, &Temp4);
          if (((Temp4.nbrLimbs != 1) || (Temp4.limbs[0].x > 1)) &&
            ((Temp4.nbrLimbs != NumberLength) ||
              memcmp(pValue->limbs, Temp4.limbs, NumberLength * sizeof(limb))))
          {          // Non-trivial factor found.
            insertBigFactor(pstFactors, &Temp4, TYP_RABIN);
            factorsFound = true;
          }
        }
        i = ctr;
        continue;  // Find more factors.
      }
      (void)memcpy(common.ecm.Aux2, common.ecm.Aux3, nbrLimbs * sizeof(limb));
    }
  }
  return factorsFound;
}

static void factorSmallInt(int intToFactor, int* factors, struct sFactors* pstFactors)
{
  int toFactor = intToFactor;
  int factorsFound = 0;
  int primeFactor;
  int multiplicity;
  struct sFactors* ptrFactor = pstFactors + 1;
  int* ptrFactorLimbs = factors;
  if (toFactor <= 3)
  {     // Only one factor.
    ptrFactor->ptrFactor = ptrFactorLimbs;
    ptrFactor->multiplicity = 1;
    ptrFactor->type = 0;
    ptrFactor->upperBound = 0;
    *ptrFactorLimbs = 1;
    ptrFactorLimbs++;
    *ptrFactorLimbs = toFactor;
    pstFactors->multiplicity = 1;
    return;
  }
  // Divide by 2 and 3.
  for (primeFactor = 2; primeFactor <= 3; primeFactor++)
  {
    multiplicity = 0;
    while ((toFactor % primeFactor) == 0)
    {
      toFactor /= primeFactor;
      multiplicity++;
    }
    if (multiplicity > 0)
    {
      factorsFound++;
      ptrFactor->ptrFactor = ptrFactorLimbs;
      ptrFactor->multiplicity = multiplicity;
      ptrFactor->type = 0;
      ptrFactor->upperBound = 0;
      *ptrFactorLimbs = 1;
      ptrFactorLimbs++;
      *ptrFactorLimbs = primeFactor;
      ptrFactorLimbs++;
      ptrFactor++;
    }
  }
  for (primeFactor = 5;
    (unsigned int)primeFactor * (unsigned int)primeFactor <= (unsigned int)toFactor;
    primeFactor += 2)
  {
    if ((primeFactor % 3) == 0)
    {
      continue;
    }
    multiplicity = 0;
    while ((toFactor % primeFactor) == 0)
    {
      toFactor /= primeFactor;
      multiplicity++;
    }
    if (multiplicity > 0)
    {
      factorsFound++;
      ptrFactor->ptrFactor = ptrFactorLimbs;
      ptrFactor->multiplicity = multiplicity;
      ptrFactor->type = 0;
      ptrFactor->upperBound = 0;
      *ptrFactorLimbs = 1;
      ptrFactorLimbs++;
      *ptrFactorLimbs = primeFactor;
      ptrFactorLimbs++;
      ptrFactor++;
    }
  }
  if (toFactor > 1)
  {
    factorsFound++;
    ptrFactor->ptrFactor = ptrFactorLimbs;
    ptrFactor->multiplicity = 1;
    ptrFactor->type = 0;
    ptrFactor->upperBound = 0;
    *ptrFactorLimbs = 1;
    ptrFactorLimbs++;
    *ptrFactorLimbs = toFactor;
  }
  pstFactors->multiplicity = factorsFound;
}

void factor(const BigInteger* toFactor, const int* number, int* factors, struct sFactors* pstFactors)
{
  factorExt(toFactor, number, factors, pstFactors, NULL);
}

// pstFactors -> ptrFactor points to end of factors.
// pstFactors -> multiplicity indicates the number of different factors.
void factorExt(const BigInteger *toFactor, const int *number, 
  int *factors, struct sFactors *pstFactors, char *pcKnownFactors)
{
  char* ptrKnownFactors = pcKnownFactors;
  struct sFactors *pstCurFactor;
  int expon;
  int remainder;
  int nbrLimbs;
  int ctr;
  const int *ptrFactor;
  int dividend;
  char *ptrCharFound;
  int result;
  initializeSmallPrimes(smallPrimes);
  if (toFactor->nbrLimbs == 1)
  {
    factorSmallInt(toFactor->limbs[0].x, factors, pstFactors);
    return;
  }
  oldNbrFactors = 0;
  NextEC = -1;
  EC = 1;
  NumberLength = toFactor->nbrLimbs;
  GetYieldFrequency();
#ifdef __EMSCRIPTEN__
  oldTimeElapsed = 0;
  lModularMult = 0;
  SIQSModMult = 0;
  primeModMult = 0;
  nbrECM = 0;
  nbrSIQS = 0;
  nbrPrimalityTests = 0;
  timeECM = 0;
  timePrimalityTests = 0;
  timeSIQS = 0;
  originalTenthSecond = tenths();
  modmultCallback = showECMStatus;   // Set callback.
#endif
  pstCurFactor = pstFactors + 1;
  if (ptrKnownFactors == NULL)
  {   // No factors known.
    (void)memcpy(factors, number, (1 + *number) * sizeof(int));
    pstFactors->multiplicity = 1;
    pstFactors->ptrFactor = factors + 1 + *factors;
    pstFactors->upperBound = 0;
    pstCurFactor->multiplicity = 1;
    pstCurFactor->ptrFactor = factors;
    pstCurFactor->upperBound = 2;
#ifdef __EMSCRIPTEN__
    SaveFactors(pstFactors);
#endif
  }
  else
  {   // Insert factors saved on Web Storage.
    pstFactors->multiplicity = 0;
    pstFactors->ptrFactor = factors;
    while (*ptrKnownFactors != 0)
    {
      ptrCharFound = findChar(ptrKnownFactors, '^');
      if (ptrCharFound == NULL)
      {
        break;
      }
      *ptrCharFound = 0;
      Dec2Bin(ptrKnownFactors, prime.limbs, (int)(ptrCharFound - ptrKnownFactors), &prime.nbrLimbs);
      BigInteger2IntArray(pstFactors->ptrFactor, &prime);
      ptrKnownFactors = ptrCharFound + 1;
      if (getNextInteger(&ptrKnownFactors, &pstCurFactor->multiplicity, '('))
      {     // Error on processing exponent.
        break;
      }
      ptrCharFound = findChar(ptrKnownFactors, ',');
      if (ptrCharFound != NULL)
      {
        if (getNextInteger(&ptrKnownFactors, &pstCurFactor->upperBound, ','))
        {     // Error on processing upper bound.
          break;
        }
        if (getNextInteger(&ptrKnownFactors, &pstCurFactor->type, ')'))
        {     // Error on processing upper bound.
          break;
        }
      }
      else
      {
        if (getNextInteger(&ptrKnownFactors, &pstCurFactor->upperBound, ')'))
        {     // Error on processing upper bound.
          break;
        }
        pstCurFactor->type = 0;
      }
      pstFactors->multiplicity++;
      pstCurFactor->ptrFactor = pstFactors->ptrFactor;
      pstFactors->ptrFactor += 1 + *pstFactors->ptrFactor;
      pstCurFactor++;
      if (*ptrKnownFactors == '*')
      {
        ptrKnownFactors++;  // Skip multiplication sign.
      }
      if (*ptrKnownFactors == ';')
      {
        ptrKnownFactors++;  // Skip separation between known factors and factor entered by user.
        Dec2Bin(ptrKnownFactors, prime.limbs, (int)strlen(ptrKnownFactors), &prime.nbrLimbs);
        if (foundByLehman)
        {
          insertBigFactor(pstFactors, &prime, TYP_LEHMAN);
        }
        else if (EC >= TYP_SIQS)
        {
          insertBigFactor(pstFactors, &prime, TYP_SIQS);
        }
        else
        {
          insertBigFactor(pstFactors, &prime, 0);
        }
#ifdef __EMSCRIPTEN__
        SaveFactors(pstFactors);
#endif
        break;
      }
      if (*ptrKnownFactors == ',')
      {
        NextEC = 0;        // Curve number
        ptrKnownFactors++;
        while (*ptrKnownFactors != '\0')
        {
          NextEC = (NextEC * 10) + (*ptrKnownFactors & 0x0F);
          ptrKnownFactors++;
        }
        oldNbrFactors = pstFactors->multiplicity;
        break;
      }
    }
  }
  if (tofactor.nbrLimbs > 1)
  {
    PowerPM1Check(pstFactors, toFactor);
  }
  pstCurFactor = pstFactors;
  for (int factorNbr = 1; factorNbr <= pstFactors->multiplicity; factorNbr++)
  {
    int upperBoundIndex;
    int upperBound;
    bool restartFactoring = false;
    pstCurFactor++;
    upperBound = pstCurFactor->upperBound;
    // If number is prime, do not process it.
    if (upperBound == 0)
    {     // Factor is prime.
      continue;
    }
    // Get upperBoundIndex from upperBound.
    upperBoundIndex = 0;
    if (upperBound < 100000)
    {
      for (int delta = 8192; delta > 0; delta >>= 1)
      {
        int tempUpperBoundIndex = upperBoundIndex + delta;
        if (tempUpperBoundIndex >= SMALL_PRIMES_ARRLEN)
        {           // Index too large.
          continue;
        }
        if (upperBound == smallPrimes[tempUpperBoundIndex])
        {
          upperBoundIndex = tempUpperBoundIndex;
          break;
        }
        if (upperBound > smallPrimes[tempUpperBoundIndex])
        {
          upperBoundIndex = tempUpperBoundIndex;
        }
      }
    }
    ptrFactor = pstCurFactor->ptrFactor;
    nbrLimbs = *ptrFactor;
    NumberLength = *pstCurFactor->ptrFactor;
    IntArray2BigInteger(pstCurFactor->ptrFactor, &power);
    NumberLength = power.nbrLimbs;
#ifdef __EMSCRIPTEN__
    char *ptrText = ShowFactoredPart(&prime, pstFactors);
    if (!skipPrimality)
    {
#endif
      expon = PowerCheck(&power, &prime);
      if (expon > 1)
      {
        NumberLength = prime.nbrLimbs;
        BigInteger2IntArray(pstCurFactor->ptrFactor, &prime);
        pstCurFactor->multiplicity *= expon;
        nbrLimbs = *pstCurFactor->ptrFactor;
      }
#ifdef __EMSCRIPTEN__
      copyStr(&ptrText, lang ? "<p>División por primos menores que 100000.</p>" :
        "<p>Trial division by primes less than 100000.</p>");
      ShowLowerText();
#endif
      while ((upperBound < 100000) && (nbrLimbs > 1))
      {        // Number has at least 2 limbs: Trial division by small numbers.
        while (pstCurFactor->upperBound != 0)
        {            // Factor found.
          int exponent;
          int index;
          int deltaIndex;
          ptrFactor = pstCurFactor->ptrFactor;
          remainder = RemDivBigNbrByInt(ptrFactor + 1, upperBound, nbrLimbs);
          if (remainder != 0)
          {    // Factor not found. Use new divisor.
            break;
          }
          // Small factor found. Find the exponent.
          exponent = 1;
          index = 0;
          deltaIndex = 1;
          CopyBigInt(&common.trialDiv.cofactor, &prime);
          subtractdivide(&common.trialDiv.cofactor, 0, upperBound);
          intToBigInteger(&common.trialDiv.power[0], upperBound);
          for (;;)
          {      // Test whether the cofactor is multiple of power.
            (void)BigIntDivide(&common.trialDiv.cofactor, &common.trialDiv.power[index], &common.trialDiv.quotient);
            (void)BigIntMultiply(&common.trialDiv.quotient, &common.trialDiv.power[index], &common.trialDiv.temp);
            if (!BigIntEqual(&common.trialDiv.temp, &common.trialDiv.cofactor))
            {    // Not a multiple, so exit loop.
              break;
            }
            CopyBigInt(&common.trialDiv.cofactor, &common.trialDiv.quotient);
            (void)BigIntMultiply(&common.trialDiv.power[index], &common.trialDiv.power[index], &common.trialDiv.power[index+1]);
            exponent += deltaIndex;
            deltaIndex <<= 1;
            index++;
          }
          for (index--; index >= 0; index--)
          {
            deltaIndex >>= 1;
            (void)BigIntDivide(&common.trialDiv.cofactor, &common.trialDiv.power[index], &common.trialDiv.quotient);
            (void)BigIntMultiply(&common.trialDiv.quotient, &common.trialDiv.power[index], &common.trialDiv.temp);
            if (BigIntEqual(&common.trialDiv.temp, &common.trialDiv.cofactor))
            {    // It is a multiple.
              CopyBigInt(&common.trialDiv.cofactor, &common.trialDiv.quotient);
              exponent += deltaIndex;
            }
          }
          insertIntFactor(pstFactors, pstCurFactor, upperBound, exponent, &common.trialDiv.cofactor);
          restartFactoring = true;
          break;
        }
        if (restartFactoring)
        {
          factorNbr = 0;
          pstCurFactor = pstFactors;
          break;
        }
        upperBoundIndex++;
        upperBound = smallPrimes[upperBoundIndex];
      }
      if (restartFactoring)
      {
        continue;
      }
      if (nbrLimbs == 1)
      {
        dividend = *(ptrFactor + 1);
        while ((upperBound < 65535) &&
                (((unsigned int)upperBound * (unsigned int)upperBound) <= (unsigned int)dividend))
        {              // Trial division by small numbers.
          if ((dividend % upperBound) == 0)
          {            // Factor found.
            insertIntFactor(pstFactors, pstCurFactor, upperBound, 1, NULL);
            restartFactoring = true;
            break;
          }
          upperBoundIndex++;
          upperBound = smallPrimes[upperBoundIndex];
        }
        if (restartFactoring)
        {
          factorNbr = 0;
          pstCurFactor = pstFactors;
          continue;
        }
        pstCurFactor->upperBound = 0;   // Number is prime.
        continue;
      }
      // No small factor. Check whether the number is prime or prime power.
#ifdef __EMSCRIPTEN__
#ifdef FACTORIZATION_APP
      StepECM = 0;
      copyStr(&ptrText, lang ? "<p>Verificando si el número es potencia perfecta.<p>" :
        "<p>Testing whether the number is perfect power or not.</p>");
      ShowLowerText();
#else
      databack(lang ? "3<p>Verificando si el número es potencia perfecta.</p>" :
        "3<p>Testing whether the number is perfect power or not.</p>");
#endif
#endif
#ifdef __EMSCRIPTEN__
    }
#endif
#ifdef __EMSCRIPTEN__
    SaveFactors(pstFactors);
#endif
#ifdef FACTORIZATION_APP
    if (skipPrimality)
    {
      skipPrimality = false;
      result = 1;
    }
    else
    {
#ifdef __EMSCRIPTEN__
      if (prime.nbrLimbs > (3322 / BITS_PER_GROUP))  // 3322 = 1000*log_2(10) -> 1000 digits
      {
        startSkipTest();
      }
      uint64_t oldModularMult = lModularMult;
#endif
      result = BpswPrimalityTest(&prime, pstFactors);
#ifdef __EMSCRIPTEN__
      nbrPrimalityTests++;
      primeModMult += lModularMult - oldModularMult;
      if (prime.nbrLimbs > (3322 / BITS_PER_GROUP))  // 3322 = 1000*log_2(10) -> 1000 digits
      {
        endSkipTest();
      }
#endif
    }
#else
    result = BpswPrimalityTest(&prime);
#endif
    if (result == 0)
    {   // Number is prime power.
      pstCurFactor->upperBound = 0;   // Indicate that number is prime.
      continue;                       // Check next factor.
    }
    if (result > 1)
    {    // Number is 2-Fermat probable prime. Try to factor it.
      if (factorCarmichael(&prime, pstFactors))
      {                               // Factors found.
        factorNbr--;                  // Test whether factor found is prime.
        pstCurFactor--;
        continue;
      }
    }
#ifdef __EMSCRIPTEN__
    double originalTenths = tenths();
#endif
    performFactorization(&prime, pstFactors);          // Factor number.
#ifdef __EMSCRIPTEN__
    nbrECM++;
    timeECM += (int)(tenths() - originalTenths);
#endif
    // Check whether GD is not one. In this case we found a proper factor.
    for (ctr = 1; ctr < NumberLength; ctr++)
    {
      if (common.ecm.GD[ctr].x != 0)
      {
        break;
      }
    }
    if ((ctr != NumberLength) || (common.ecm.GD[0].x != 1))
    {
      int numLimbs;
      Temp1.sign = SIGN_POSITIVE;
      numLimbs = NumberLength;
      while (numLimbs > 1)
      {
        if (common.ecm.GD[numLimbs-1].x != 0)
        {
          break;
        }
        numLimbs--;
      }
      (void)memcpy(Temp1.limbs, common.ecm.GD, numLimbs * sizeof(limb));
      Temp1.nbrLimbs = numLimbs;
      if (foundByLehman)
      {
        insertBigFactor(pstFactors, &Temp1, TYP_LEHMAN + EC);
      }
      else
      {
        insertBigFactor(pstFactors, &Temp1, EC);
      }
#ifdef __EMSCRIPTEN__
      SaveFactors(pstFactors);
#endif
      factorNbr = 0;
      pstCurFactor = pstFactors;
    }    // End if
  }      // End for
#ifdef __EMSCRIPTEN__
  SaveFactors(pstFactors);
#endif
}

EXTERNALIZE char *getFactorsAsciiPtr(void)
{
  return common.saveFactors.text;
}

static void intArrayToBigInteger(const int *ptrValues, BigInteger *bigint)
{
  bigint->sign = SIGN_POSITIVE;
  bigint->nbrLimbs = *ptrValues;
  (void)memcpy(bigint->limbs, ptrValues + 1, *ptrValues * sizeof(int));
}

// Find Euler's Totient as the product of p^(e-1)*(p-1) where p=prime and e=exponent.
void Totient(BigInteger *result)
{
  static BigInteger TempVar;
  const struct sFactors *pstFactor;
  intToBigInteger(result, 1);  // Set result to 1.
  pstFactor = &astFactorsMod[1];
  for (int factorNumber = 1; factorNumber <= astFactorsMod[0].multiplicity; factorNumber++)
  {
    intArrayToBigInteger(pstFactor->ptrFactor, &factorValue);
    if ((factorValue.nbrLimbs == 1) && (factorValue.limbs[0].x == 1))
    {   // If factor is 1 do not do anything.
      continue;
    }
    (void)BigIntPowerIntExp(&factorValue, pstFactor->multiplicity - 1, &TempVar);   // p^(e-1)
    (void)BigIntMultiply(result, &TempVar, result);
    intArrayToBigInteger(pstFactor->ptrFactor, &TempVar);
    addbigint(&TempVar, -1);   // p-1
    (void)BigIntMultiply(result, &TempVar, result);
    pstFactor++;
  }
}

void NumFactors(BigInteger* result)
{
  intToBigInteger(result, astFactorsMod[0].multiplicity);
}

void MinFactor(BigInteger* result)
{
  const struct sFactors* pstFactor = &astFactorsMod[1];
  intArrayToBigInteger(pstFactor->ptrFactor, result);
  for (int factorNumber = 2; factorNumber <= astFactorsMod[0].multiplicity; factorNumber++)
  {
    pstFactor++;
    intArrayToBigInteger(pstFactor->ptrFactor, &factorValue);
    BigIntSubt(&factorValue, result, &factorValue);
    if (factorValue.sign == SIGN_NEGATIVE)
    {
      intArrayToBigInteger(pstFactor->ptrFactor, result);
    }
  }
}

void MaxFactor(BigInteger* result)
{
  const struct sFactors* pstFactor = &astFactorsMod[1];
  intArrayToBigInteger(pstFactor->ptrFactor, result);
  for (int factorNumber = 2; factorNumber <= astFactorsMod[0].multiplicity; factorNumber++)
  {
    pstFactor++;
    intArrayToBigInteger(pstFactor->ptrFactor, &factorValue);
    BigIntSubt(&factorValue, result, &factorValue);
    if (factorValue.sign == SIGN_POSITIVE)
    {
      intArrayToBigInteger(pstFactor->ptrFactor, result);
    }
  }
}

// Find sum of divisors as the product of (p^(e+1)-1)/(p-1) where p=prime and e=exponent.
void SumOfDivisors(BigInteger *result)
{
  const struct sFactors *pstFactor;
  intToBigInteger(result, 1);  // Set result to 1.
  pstFactor = &astFactorsMod[1];
  for (int factorNumber = 1; factorNumber <= astFactorsMod[0].multiplicity; factorNumber++)
  {
    intArrayToBigInteger(pstFactor->ptrFactor, &factorValue);
    if ((factorValue.nbrLimbs == 1) && (factorValue.limbs[0].x == 1))
    {   // If factor is 1 do not do anything.
      continue;
    }
    (void)BigIntPowerIntExp(&factorValue, pstFactor->multiplicity + 1, &Temp1);   // p^(e+1)
    addbigint(&Temp1, -1);   // p^(e+1)-1
    (void)BigIntMultiply(result, &Temp1, &Temp2);
    intArrayToBigInteger(pstFactor->ptrFactor, &Temp1);
    addbigint(&Temp1, -1);   // p-1
    (void)BigIntDivide(&Temp2, &Temp1, result);
    pstFactor++;
  }
}

// Find number of divisors as the product of e+1 where p=prime and e=exponent.
void NumberOfDivisors(BigInteger *result)
{
  const struct sFactors *pstFactor;
  intToBigInteger(result, 1);  // Set result to 1.
  pstFactor = &astFactorsMod[1];
  for (int factorNumber = 1; factorNumber <= astFactorsMod[0].multiplicity; factorNumber++)
  {
    intArrayToBigInteger(pstFactor->ptrFactor, &factorValue);
    if ((factorValue.nbrLimbs == 1) && (factorValue.limbs[0].x == 1))
    {   // If factor is 1 do not do anything.
      continue;
    }
    multint(result, result, pstFactor->multiplicity + 1);
    pstFactor++;
  }
}

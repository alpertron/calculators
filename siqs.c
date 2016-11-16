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
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "bignbr.h"
#include "expression.h"
#include "factor.h"
// These defines are valid for factoring up to 10^110.
#define MAX_NBR_FACTORS        13
#define MAX_PRIMES         150000
#define MAX_LIMBS_SIQS         15
#define MAX_FACTORS_RELATION   50
#define LENGTH_OFFSET           0

#ifdef __EMSCRIPTEN__
void databack(char *data);
int stamp(void);
int tenths(void);
extern int newStamp, oldStamp;
extern char lowerText[], *ptrLowerText;
char *ptrSIQSStrings;
int oldTimeElapsed, startSieveTenths;
double originalTenthSecond;
#endif

typedef struct
{
  int value;
  int modsqrt;
  int Bainv2[MAX_NBR_FACTORS];
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

unsigned char SIQSInfoText[300];
int numberThreads = 1;
extern int NumberLength;
extern int groupLen;
int matrixBLength;
long trialDivisions;
long smoothsFound;
long totalPartials;
long partialsFound;
long ValuesSieved;
int nbrPrimes;
int congruencesFound;
long polynomialsSieved;
int nbrPartials;
int multiplier;
int nbrFactorsA;
int afact[MAX_NBR_FACTORS];
long startTime;
int Modulus[MAX_LIMBS_SIQS];
int TestNbr2[MAX_LIMBS_SIQS];
int biQuadrCoeff[MAX_LIMBS_SIQS];
int biLinearDelta[MAX_LIMBS_SIQS][MAX_LIMBS_SIQS];
static long largePrimeUpperBound;
static unsigned char logar2;
static int aindex[MAX_NBR_FACTORS];
//  static Thread threadArray[];
static PrimeSieveData primeSieveData[MAX_PRIMES+3];
static PrimeTrialDivisionData primeTrialDivisionData[MAX_PRIMES];
static int span;
static int indexMinFactorA;
static int threadNumber;
static int nbrThreadFinishedPolySet;
static unsigned int oldSeed;
static unsigned int newSeed;
static int NbrPolynomials;
static int SieveLimit;
static int matrixPartial[MAX_PRIMES * 8][MAX_LIMBS_SIQS/2+4];
static int matrixPartialLength;
static int vectLeftHandSide[MAX_PRIMES+50][MAX_FACTORS_RELATION];
static int matrixPartialHashIndex[2048];
static int matrixB[MAX_PRIMES + 50][MAX_FACTORS_RELATION];
static int amodq[MAX_NBR_FACTORS];
static int tmodqq[MAX_NBR_FACTORS];
static char threshold;
static int smallPrimeUpperLimit;
static int firstLimit;
static int secondLimit;
static int thirdLimit;
static int vectExpParity[MAX_PRIMES + 50];
static int matrixAV[MAX_PRIMES];
static int matrixV[MAX_PRIMES];
static int matrixV1[MAX_PRIMES];
static int matrixV2[MAX_PRIMES];
static int matrixXmY[MAX_PRIMES];
static int newColumns[MAX_PRIMES];
// Matrix that holds temporary data
static int matrixCalc3[MAX_PRIMES];
static int matrixTemp2[MAX_PRIMES];
static int nbrPrimes2;
static BigInteger factorSiqs;
static unsigned char onlyFactoring;
PrimeSieveData *firstPrimeSieveData;
static unsigned char InsertNewRelation(
  int *rowMatrixB,
  int *biT, int *biU, int *biR,
  int NumberLength);
static void BlockLanczos(void);
void ShowSIQSStatus(void);
static unsigned int getFactorsOfA(unsigned int seed, int *aindex);
static void sieveThread(BigInteger *result);

#ifdef __EMSCRIPTEN__
// Convert tenths of seconds to days, hours, minutes and seconds.
void GetDHMS(char **pptrText, int seconds)
{
  char *ptrText = *pptrText;
  int2dec(&ptrText, seconds / 86400);         // Show number of days.
  *ptrText++ = 'd';
  *ptrText++ = ' ';
  int2dec(&ptrText, (seconds / 3600) % 24);   // Show number of hours.
  *ptrText++ = 'h';
  *ptrText++ = ' ';
  int2dec(&ptrText, (seconds / 60) % 60);     // Show number of minutes.
  *ptrText++ = 'm';
  *ptrText++ = ' ';
  int2dec(&ptrText, seconds % 60);            // Show number of seconds.
  *ptrText++ = 's';
  *ptrText++ = ' ';
  *pptrText = ptrText;
}

static void showMatrixSize(char *SIQSInfoText, int rows, int cols)
{
  char *ptrText = ptrLowerText;  // Point after number that is being factored.
  strcpy(ptrText, lang ? "<p>Resolviendo la matriz de congruencias de " : "<p>Solving ");
  ptrText += strlen(ptrText);
  int2dec(&ptrText, rows);   // Show number of rows.
  strcpy(ptrText, " &times; ");
  ptrText += strlen(ptrText);
  int2dec(&ptrText, cols);   // Show number of columns.
  strcpy(ptrText, lang ? "usando el algoritmo de Lanczos en bloques.</p>" :
                         "congruence matrix using Block Lanczos algorithm.</p>");
  databack(lowerText);
}

static void InitSIQSStrings(int SieveLimit)
{
  char *ptrText = ptrLowerText;  // Point after number that is being factored.
  strcpy(ptrText, lang ? "<p>Parámetros de SIQS: " : "<p>SIQS parameters: ");
  ptrText += strlen(ptrText);
  int2dec(&ptrText, nbrPrimes);   // Show number of primes in factor base.
  strcpy(ptrText, lang ? " primos, límite de la criba: " : " primes, sieve limit: ");
  ptrText += strlen(ptrText);
  int2dec(&ptrText, SieveLimit);  // Show sieve limit.
  strcpy(ptrText, "</p>");
  ptrText += strlen(ptrText);
  ptrSIQSStrings = ptrText;
  strcpy(ptrText, lang ? "<p>Buscando el mejor multiplicador de Knuth-Schroeppel...</p>" :
                         "<p>Searching for Knuth-Schroeppel multiplier...</p>");
  databack(lowerText);
}

// Append multiplier and factor base to SIQS string.
static void getMultAndFactorBase(int multiplier, int FactorBase)
{
  char *ptrText = ptrSIQSStrings;
  strcpy(ptrText, lang ? "<p>Multiplicador: " : "<p>Multiplier: ");
  ptrText += strlen(ptrText);
  int2dec(&ptrText, multiplier);  // Show Knuth-Schroeppel multiplier.
  strcpy(ptrText, lang ? ", base de factores: " : ", factor base: ");
  ptrText += strlen(ptrText);
  int2dec(&ptrText, FactorBase);  // Show factor base.
  strcpy(ptrText, "</p>");
  ptrText += strlen(ptrText);
  ptrSIQSStrings = ptrText;
}

static void ShowSIQSInfo(int timeSieve, int congruencesFound, int matrixBLength, int elapsedTime)
{
  char SIQSInfo[100];
  int percentage = (int)((float)(congruencesFound * 100) / (float)matrixBLength);
  int u = (int)((double)timeSieve * (double)(matrixBLength - congruencesFound) / (double)congruencesFound);
  char *ptrText = SIQSInfo;
  strcpy(ptrText, "4<p>");
  ptrText += strlen(ptrText);
  GetDHMS(&ptrText, elapsedTime);
  strcpy(ptrText, "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;");
  ptrText += strlen(ptrText);
  int2dec(&ptrText, congruencesFound);  // Show number of coungruences found.
  strcpy(ptrText, lang ? " congruencias halladas (" : " congruences found (");
  ptrText += strlen(ptrText);
  int2dec(&ptrText, percentage);  // Show number of coungruences found.
  if (timeSieve > 1 && congruencesFound > 10)
  {
    strcpy(ptrText, lang ? "%). Fin de la criba en " : "%). End sieve in ");
    ptrText += strlen(ptrText);
    GetDHMS(&ptrText, u/2);
    strcpy(ptrText, "</p>");
  }
  else
  {
    strcpy(ptrText, "%)");
  }
  databack(SIQSInfo);
}

#endif

static void PerformSiqsSieveStage(PrimeSieveData *primeSieveData,
  short SieveArray[],
  int PolynomialIndex,
  int *biLinearCoeff,
  int NumberLength)
{
  short logPrimeEvenPoly, logPrimeOddPoly;
  int currentPrime, F1, F2, F3, F4, X1, X2;
  int index, index2, indexFactorA;
  int mask;
  unsigned char polyadd;
  int S1, G0, G1, G2, G3;
  int H0, H1, H2, H3, I0, I1, I2, I3;
  PrimeSieveData *rowPrimeSieveData;

  F1 = PolynomialIndex;
  indexFactorA = 0;
  while ((F1 & 1) == 0)
  {
    F1 >>= 1;
    indexFactorA++;
  }
  polyadd = (F1 & 2) != 0;
  if (polyadd)   // Adjust value of B as appropriate
  {                                // according to the Gray code.
    AddBigNbr(biLinearCoeff, biLinearDelta[indexFactorA], biLinearCoeff,
      NumberLength);
    AddBigNbr(biLinearCoeff, biLinearDelta[indexFactorA], biLinearCoeff,
      NumberLength);
  }
  else
  {
    SubtractBigNbr(biLinearCoeff, biLinearDelta[indexFactorA],
      biLinearCoeff, NumberLength);
    SubtractBigNbr(biLinearCoeff, biLinearDelta[indexFactorA],
      biLinearCoeff, NumberLength);
  }
  indexFactorA--;
  X1 = SieveLimit << 1;
  rowPrimeSieveData = primeSieveData + 1;
  F1 = polyadd ? -rowPrimeSieveData -> Bainv2[indexFactorA] :
    rowPrimeSieveData->Bainv2[indexFactorA];
  if (((rowPrimeSieveData->soln1 += F1) & 1) == 0)
  {
    SieveArray[0] = (short)(logar2 - threshold);
    SieveArray[1] = (short)(-threshold);
  }
  else
  {
    SieveArray[0] = (short)(-threshold);
    SieveArray[1] = (short)(logar2 - threshold);
  }
  if (((rowPrimeSieveData->soln1 + rowPrimeSieveData->Bainv2_0) & 1) == 0)
  {
    SieveArray[0] += (short)((logar2 - threshold) << 8);
    SieveArray[1] += (short)((-threshold) << 8);
  }
  else
  {
    SieveArray[0] += (short)((-threshold) << 8);
    SieveArray[1] += (short)((logar2 - threshold) << 8);
  }
  F2 = 2;
  index = 2;
  for (;;)
  {
    rowPrimeSieveData = primeSieveData + index;
    currentPrime = rowPrimeSieveData->value;
    F3 = F2 * currentPrime;
    if (X1 + 1 < F3)
    {
      F3 = X1 + 1;
    }
    F4 = F2;
    while (F4 * 2 <= F3)
    {
      memcpy(&SieveArray[F4], &SieveArray[0], F4*sizeof(SieveArray[0]));
      F4 *= 2;
    }
    memcpy(&SieveArray[F4], &SieveArray[0], (F3-F4) * sizeof(SieveArray[0]));
    if (F3 == X1 + 1)
    {
      break;
    }
    F1 = currentPrime;
    logPrimeEvenPoly = 1;
    while (F1 >= 5)
    {
      F1 /= 3;
      logPrimeEvenPoly++;
    }
    logPrimeOddPoly = (short)(logPrimeEvenPoly << 8);
    F1 = polyadd ? -rowPrimeSieveData->Bainv2[indexFactorA] :
      rowPrimeSieveData->Bainv2[indexFactorA];
    index2 = (rowPrimeSieveData->soln1 + F1) % currentPrime;
    rowPrimeSieveData->soln1 = index2 += currentPrime & (index2 >> 31);
    for (; index2 < F3; index2 += currentPrime)
    {
      SieveArray[index2] += logPrimeEvenPoly;
    }
    for (index2 = (rowPrimeSieveData->soln1 + currentPrime -
      rowPrimeSieveData->Bainv2_0) % currentPrime;
      index2 < F3;
      index2 += currentPrime)
    {
      SieveArray[index2] += logPrimeOddPoly;
    }
    if (currentPrime != multiplier)
    {
      for (F1 = index2 = (rowPrimeSieveData->soln1 + currentPrime -
        rowPrimeSieveData->difsoln) % currentPrime;
        index2 < F3;
        index2 += currentPrime)
      {
        SieveArray[index2] += logPrimeEvenPoly;
      }
      for (index2 = (F1 + currentPrime -
        rowPrimeSieveData->Bainv2_0) % currentPrime;
        index2 < F3;
        index2 += currentPrime)
      {
        SieveArray[index2] += logPrimeOddPoly;
      }
    }
    index++;
    F2 *= currentPrime;
  }

  F1 = (primeSieveData + smallPrimeUpperLimit) -> value;
  logPrimeEvenPoly = 1;
  logPrimeOddPoly = 0x100;
  mask = 5;
  while (F1 >= 5)
  {
    F1 /= 3;
    logPrimeEvenPoly++;
    logPrimeOddPoly += 0x100;
    mask *= 3;
  }
  if (polyadd)
  {
    for (; index < smallPrimeUpperLimit; index++)
    {
      rowPrimeSieveData = primeSieveData + index;
      currentPrime = rowPrimeSieveData->value;
      if ((S1 = rowPrimeSieveData->soln1 -
        rowPrimeSieveData->Bainv2[indexFactorA]) < 0)
      {
        S1 += currentPrime;
      }
      rowPrimeSieveData->soln1 = S1;
    }
    for (index = smallPrimeUpperLimit; index < firstLimit; index++)
    {
      rowPrimeSieveData = primeSieveData + index;
      currentPrime = rowPrimeSieveData->value;
      if (currentPrime >= mask)
      {
        mask *= 3;
        logPrimeEvenPoly++;
        logPrimeOddPoly += 0x100;
      }
      F2 = currentPrime + currentPrime;
      F3 = F2 + currentPrime;
      F4 = F3 + currentPrime;
      S1 = rowPrimeSieveData->soln1 -
        rowPrimeSieveData->Bainv2[indexFactorA];
      rowPrimeSieveData->soln1 = S1 += (S1 >> 31) & currentPrime;
      index2 = X1 / F4 * F4 + S1;
      G0 = -rowPrimeSieveData->difsoln;
      if (S1 + G0 < 0)
      {
        G0 += currentPrime;
      }
      G1 = G0 + currentPrime;
      G2 = G1 + currentPrime;
      G3 = G2 + currentPrime;
      H0 = -rowPrimeSieveData->Bainv2_0;
      if (S1 + H0 < 0)
      {
        H0 += currentPrime;
      }
      H1 = H0 + currentPrime;
      H2 = H1 + currentPrime;
      H3 = H2 + currentPrime;
      I0 = H0 - rowPrimeSieveData->difsoln;
      if (S1 + I0 < 0)
      {
        I0 += currentPrime;
      }
      I1 = I0 + currentPrime;
      I2 = I1 + currentPrime;
      I3 = I2 + currentPrime;
      do
      {
        SieveArray[index2] += logPrimeEvenPoly;
        SieveArray[index2 + currentPrime] += logPrimeEvenPoly;
        SieveArray[index2 + F2] += logPrimeEvenPoly;
        SieveArray[index2 + F3] += logPrimeEvenPoly;
        SieveArray[index2 + G0] += logPrimeEvenPoly;
        SieveArray[index2 + G1] += logPrimeEvenPoly;
        SieveArray[index2 + G2] += logPrimeEvenPoly;
        SieveArray[index2 + G3] += logPrimeEvenPoly;
        SieveArray[index2 + H0] += logPrimeOddPoly;
        SieveArray[index2 + H1] += logPrimeOddPoly;
        SieveArray[index2 + H2] += logPrimeOddPoly;
        SieveArray[index2 + H3] += logPrimeOddPoly;
        SieveArray[index2 + I0] += logPrimeOddPoly;
        SieveArray[index2 + I1] += logPrimeOddPoly;
        SieveArray[index2 + I2] += logPrimeOddPoly;
        SieveArray[index2 + I3] += logPrimeOddPoly;
      } while ((index2 -= F4) >= 0);
    }
    for (; index < secondLimit; index++)
    {
      rowPrimeSieveData = primeSieveData + index;
      currentPrime = rowPrimeSieveData->value;
      F2 = currentPrime + currentPrime;
      F3 = F2 + currentPrime;
      F4 = F2 + F2;
      X2 = X1 - F4;
      if (currentPrime >= mask)
      {
        mask *= 3;
        logPrimeEvenPoly++;
        logPrimeOddPoly += 0x100;
      }
      if (rowPrimeSieveData->difsoln >= 0)
      {
        F1 = rowPrimeSieveData->soln1 -
          rowPrimeSieveData->Bainv2[indexFactorA];
        F1 += (F1 >> 31) & currentPrime;
        index2 = (rowPrimeSieveData->soln1 = F1);
        do
        {
          SieveArray[index2] += logPrimeEvenPoly;
          SieveArray[index2 + currentPrime] += logPrimeEvenPoly;
          SieveArray[index2 + F2] += logPrimeEvenPoly;
          SieveArray[index2 + F3] += logPrimeEvenPoly;
        } while ((index2 += F4) <= X2);
        for (; index2 <= X1; index2 += currentPrime)
        {
          SieveArray[index2] += logPrimeEvenPoly;
        }
        index2 = F1 - rowPrimeSieveData->Bainv2_0;
        index2 += (index2 >> 31) & currentPrime;
        do
        {
          SieveArray[index2] += logPrimeOddPoly;
          SieveArray[index2 + currentPrime] += logPrimeOddPoly;
          SieveArray[index2 + F2] += logPrimeOddPoly;
          SieveArray[index2 + F3] += logPrimeOddPoly;
        } while ((index2 += F4) <= X2);
        for (; index2 <= X1; index2 += currentPrime)
        {
          SieveArray[index2] += logPrimeOddPoly;
        }
        F1 -= rowPrimeSieveData->difsoln;
        F1 += (F1 >> 31) & currentPrime;
        index2 = F1;
        do
        {
          SieveArray[index2] += logPrimeEvenPoly;
          SieveArray[index2 + currentPrime] += logPrimeEvenPoly;
          SieveArray[index2 + F2] += logPrimeEvenPoly;
          SieveArray[index2 + F3] += logPrimeEvenPoly;
        } while ((index2 += F4) <= X2);
        for (; index2 <= X1; index2 += currentPrime)
        {
          SieveArray[index2] += logPrimeEvenPoly;
        }
        index2 = F1 - rowPrimeSieveData->Bainv2_0;
        index2 += (index2 >> 31) & currentPrime;
        do
        {
          SieveArray[index2] += logPrimeOddPoly;
          SieveArray[index2 + currentPrime] += logPrimeOddPoly;
          SieveArray[index2 + F2] += logPrimeOddPoly;
          SieveArray[index2 + F3] += logPrimeOddPoly;
        } while ((index2 += F4) <= X2);
        for (; index2 <= X1; index2 += currentPrime)
        {
          SieveArray[index2] += logPrimeOddPoly;
        }
      }
    }
    for (; index < thirdLimit; index++)
    {
      rowPrimeSieveData = primeSieveData + index;
      currentPrime = rowPrimeSieveData->value;
      if (currentPrime >= mask)
      {
        mask *= 3;
        logPrimeEvenPoly++;
        logPrimeOddPoly += 0x100;
      }
      F2 = rowPrimeSieveData->soln1 - rowPrimeSieveData->Bainv2[indexFactorA];
      F2 += currentPrime & (F2 >> 31);
      index2 = (rowPrimeSieveData->soln1 = F2);
      do
      {
        SieveArray[index2] += logPrimeEvenPoly;
      } while ((index2 += currentPrime) <= X1);
      F1 = F2 - (F3 = rowPrimeSieveData->Bainv2_0);
      F1 += currentPrime & (F1 >> 31);
      do
      {
        SieveArray[F1] += logPrimeOddPoly;
      } while ((F1 += currentPrime) <= X1);
      F2 -= rowPrimeSieveData->difsoln;
      index2 = F2 += currentPrime & (F2 >> 31);
      do
      {
        SieveArray[index2] += logPrimeEvenPoly;
      } while ((index2 += currentPrime) <= X1);
      F2 += (currentPrime & ((F2 - F3) >> 31)) - F3;
      do
      {
        SieveArray[F2] += logPrimeOddPoly;
      } while ((F2 += currentPrime) <= X1);
    }
    for (; index < nbrPrimes2; index++)
    {
      rowPrimeSieveData = primeSieveData + index;
      currentPrime = rowPrimeSieveData->value;
      if (currentPrime >= mask)
      {
        mask *= 3;
        logPrimeEvenPoly++;
        logPrimeOddPoly += 0x100;
      }
      F2 = rowPrimeSieveData->soln1 - rowPrimeSieveData->Bainv2[indexFactorA];
      if ((rowPrimeSieveData->soln1 = (F2 += currentPrime & (F2 >> 31))) < X1)
      {
        SieveArray[F2] += logPrimeEvenPoly;
      }
      F1 = F2 - (F3 = rowPrimeSieveData->Bainv2_0);
      if ((F1 += currentPrime & (F1 >> 31)) < X1)
      {
        SieveArray[F1] += logPrimeOddPoly;
      }
      F2 -= rowPrimeSieveData->difsoln;
      if ((F2 += currentPrime & (F2 >> 31)) < X1)
      {
        SieveArray[F2] += logPrimeEvenPoly;
      }
      if ((F2 += (currentPrime & ((F2 - F3) >> 31)) - F3) < X1)
      {
        SieveArray[F2] += logPrimeOddPoly;
      }
      rowPrimeSieveData = primeSieveData + ++index;
      currentPrime = rowPrimeSieveData->value;
      F2 = rowPrimeSieveData->soln1 - rowPrimeSieveData->Bainv2[indexFactorA];
      if ((rowPrimeSieveData->soln1 = (F2 += currentPrime & (F2 >> 31))) < X1)
      {
        SieveArray[F2] += logPrimeEvenPoly;
      }
      F1 = F2 - (F3 = rowPrimeSieveData->Bainv2_0);
      if ((F1 += currentPrime & (F1 >> 31)) < X1)
      {
        SieveArray[F1] += logPrimeOddPoly;
      }
      F2 -= rowPrimeSieveData->difsoln;
      if ((F2 += currentPrime & (F2 >> 31)) < X1)
      {
        SieveArray[F2] += logPrimeEvenPoly;
      }
      if ((F2 += (currentPrime & ((F2 - F3) >> 31)) - F3) < X1)
      {
        SieveArray[F2] += logPrimeOddPoly;
      }
      rowPrimeSieveData = primeSieveData + ++index;
      currentPrime = rowPrimeSieveData->value;
      F2 = rowPrimeSieveData->soln1 - rowPrimeSieveData->Bainv2[indexFactorA];
      if ((rowPrimeSieveData->soln1 = (F2 += currentPrime & (F2 >> 31))) < X1)
      {
        SieveArray[F2] += logPrimeEvenPoly;
      }
      F1 = F2 - (F3 = rowPrimeSieveData->Bainv2_0);
      if ((F1 += currentPrime & (F1 >> 31)) < X1)
      {
        SieveArray[F1] += logPrimeOddPoly;
      }
      F2 -= rowPrimeSieveData->difsoln;
      if ((F2 += currentPrime & (F2 >> 31)) < X1)
      {
        SieveArray[F2] += logPrimeEvenPoly;
      }
      F2 -= F3;
      if ((F2 += currentPrime & (F2 >> 31)) < X1)
      {
        SieveArray[F2] += logPrimeOddPoly;
      }
      rowPrimeSieveData = primeSieveData + ++index;
      currentPrime = rowPrimeSieveData->value;
      F2 = rowPrimeSieveData->soln1 - rowPrimeSieveData->Bainv2[indexFactorA];
      if ((rowPrimeSieveData->soln1 = (F2 += currentPrime & (F2 >> 31))) < X1)
      {
        SieveArray[F2] += logPrimeEvenPoly;
      }
      F1 = F2 - (F3 = rowPrimeSieveData->Bainv2_0);
      if ((F1 += currentPrime & (F1 >> 31)) < X1)
      {
        SieveArray[F1] += logPrimeOddPoly;
      }
      F2 -= rowPrimeSieveData->difsoln;
      if ((F2 += currentPrime & (F2 >> 31)) < X1)
      {
        SieveArray[F2] += logPrimeEvenPoly;
      }
      F2 -= F3;
      if ((F2 += currentPrime & (F2 >> 31)) < X1)
      {
        SieveArray[F2] += logPrimeOddPoly;
      }
    }
    for (; index < nbrPrimes; index++)
    {
      rowPrimeSieveData = primeSieveData + index;
      currentPrime = rowPrimeSieveData->value;
      if (currentPrime >= mask)
      {
        mask *= 3;
        logPrimeEvenPoly++;
        logPrimeOddPoly += 0x100;
      }
      F2 = rowPrimeSieveData->soln1 - rowPrimeSieveData->Bainv2[indexFactorA];
      if ((rowPrimeSieveData->soln1 = (F2 += currentPrime & (F2 >> 31))) < X1)
      {
        SieveArray[F2] += logPrimeEvenPoly;
      }
      F1 = F2 - (F3 = rowPrimeSieveData->Bainv2_0);
      if ((F1 += currentPrime & (F1 >> 31)) < X1)
      {
        SieveArray[F1] += logPrimeOddPoly;
      }
      F2 -= rowPrimeSieveData->difsoln;
      if ((F2 += currentPrime & (F2 >> 31)) < X1)
      {
        SieveArray[F2] += logPrimeEvenPoly;
      }
      F2 -= F3;
      if ((F2 += currentPrime & (F2 >> 31)) < X1)
      {
        SieveArray[F2] += logPrimeOddPoly;
      }
    }
  }
  else
  {
    for (; index < smallPrimeUpperLimit; index++)
    {
      rowPrimeSieveData = primeSieveData + index;
      currentPrime = rowPrimeSieveData->value;
      S1 = rowPrimeSieveData->soln1 +
        rowPrimeSieveData->Bainv2[indexFactorA] - currentPrime;
      S1 += currentPrime & (S1 >> 31);
      rowPrimeSieveData->soln1 = S1;
    }
    for (index = smallPrimeUpperLimit; index < firstLimit; index++)
    {
      rowPrimeSieveData = primeSieveData + index;
      currentPrime = rowPrimeSieveData->value;
      if (currentPrime >= mask)
      {
        mask *= 3;
        logPrimeEvenPoly++;
        logPrimeOddPoly += 0x100;
      }
      F2 = currentPrime + currentPrime;
      F3 = F2 + currentPrime;
      F4 = F3 + currentPrime;
      S1 = rowPrimeSieveData->soln1 +
        rowPrimeSieveData->Bainv2[indexFactorA] - currentPrime;
      rowPrimeSieveData->soln1 = S1 += (S1 >> 31) & currentPrime;
      index2 = X1 / F4 * F4 + S1;
      G0 = -rowPrimeSieveData->difsoln;
      if (S1 + G0 < 0)
      {
        G0 += currentPrime;
      }
      G1 = G0 + currentPrime;
      G2 = G1 + currentPrime;
      G3 = G2 + currentPrime;
      H0 = -rowPrimeSieveData->Bainv2_0;
      if (S1 + H0 < 0)
      {
        H0 += currentPrime;
      }
      H1 = H0 + currentPrime;
      H2 = H1 + currentPrime;
      H3 = H2 + currentPrime;
      I0 = H0 - rowPrimeSieveData->difsoln;
      if (S1 + I0 < 0)
      {
        I0 += currentPrime;
      }
      I1 = I0 + currentPrime;
      I2 = I1 + currentPrime;
      I3 = I2 + currentPrime;
      do
      {
        SieveArray[index2] += logPrimeEvenPoly;
        SieveArray[index2 + currentPrime] += logPrimeEvenPoly;
        SieveArray[index2 + F2] += logPrimeEvenPoly;
        SieveArray[index2 + F3] += logPrimeEvenPoly;
        SieveArray[index2 + G0] += logPrimeEvenPoly;
        SieveArray[index2 + G1] += logPrimeEvenPoly;
        SieveArray[index2 + G2] += logPrimeEvenPoly;
        SieveArray[index2 + G3] += logPrimeEvenPoly;
        SieveArray[index2 + H0] += logPrimeOddPoly;
        SieveArray[index2 + H1] += logPrimeOddPoly;
        SieveArray[index2 + H2] += logPrimeOddPoly;
        SieveArray[index2 + H3] += logPrimeOddPoly;
        SieveArray[index2 + I0] += logPrimeOddPoly;
        SieveArray[index2 + I1] += logPrimeOddPoly;
        SieveArray[index2 + I2] += logPrimeOddPoly;
        SieveArray[index2 + I3] += logPrimeOddPoly;
      } while ((index2 -= F4) >= 0);
    }
    for (; index < secondLimit; index++)
    {
      rowPrimeSieveData = primeSieveData + index;
      currentPrime = rowPrimeSieveData->value;
      F2 = currentPrime + currentPrime;
      F3 = F2 + currentPrime;
      F4 = F2 + F2;
      X2 = X1 - F4;
      if (currentPrime >= mask)
      {
        mask *= 3;
        logPrimeEvenPoly++;
        logPrimeOddPoly += 0x100;
      }
      if (rowPrimeSieveData->difsoln >= 0)
      {
        F1 = rowPrimeSieveData->soln1 +
          rowPrimeSieveData->Bainv2[indexFactorA] - currentPrime;
        F1 += currentPrime & (F1 >> 31);
        index2 = (rowPrimeSieveData->soln1 = F1);
        do
        {
          SieveArray[index2] += logPrimeEvenPoly;
          SieveArray[index2 + currentPrime] += logPrimeEvenPoly;
          SieveArray[index2 + F2] += logPrimeEvenPoly;
          SieveArray[index2 + F3] += logPrimeEvenPoly;
        } while ((index2 += F4) <= X2);
        for (; index2 <= X1; index2 += currentPrime)
        {
          SieveArray[index2] += logPrimeEvenPoly;
        }
        index2 = F1 - rowPrimeSieveData->Bainv2_0;
        index2 += (index2 >> 31) & currentPrime;
        do
        {
          SieveArray[index2] += logPrimeOddPoly;
          SieveArray[index2 + currentPrime] += logPrimeOddPoly;
          SieveArray[index2 + F2] += logPrimeOddPoly;
          SieveArray[index2 + F3] += logPrimeOddPoly;
        } while ((index2 += F4) <= X2);
        for (; index2 <= X1; index2 += currentPrime)
        {
          SieveArray[index2] += logPrimeOddPoly;
        }
        F1 -= rowPrimeSieveData->difsoln;
        F1 += (F1 >> 31) & currentPrime;
        index2 = F1;
        do
        {
          SieveArray[index2] += logPrimeEvenPoly;
          SieveArray[index2 + currentPrime] += logPrimeEvenPoly;
          SieveArray[index2 + F2] += logPrimeEvenPoly;
          SieveArray[index2 + F3] += logPrimeEvenPoly;
        } while ((index2 += F4) <= X2);
        for (; index2 <= X1; index2 += currentPrime)
        {
          SieveArray[index2] += logPrimeEvenPoly;
        }
        index2 = F1 - rowPrimeSieveData->Bainv2_0;
        index2 += (index2 >> 31) & currentPrime;
        do
        {
          SieveArray[index2] += logPrimeOddPoly;
          SieveArray[index2 + currentPrime] += logPrimeOddPoly;
          SieveArray[index2 + F2] += logPrimeOddPoly;
          SieveArray[index2 + F3] += logPrimeOddPoly;
        } while ((index2 += F4) <= X2);
        for (; index2 <= X1; index2 += currentPrime)
        {
          SieveArray[index2] += logPrimeOddPoly;
        }
      }
    }
    for (; index < thirdLimit; index++)
    {
      rowPrimeSieveData = primeSieveData + index;
      currentPrime = rowPrimeSieveData->value;
      if (currentPrime >= mask)
      {
        mask *= 3;
        logPrimeEvenPoly++;
        logPrimeOddPoly += 0x100;
      }
      F2 = rowPrimeSieveData->soln1 +
        rowPrimeSieveData->Bainv2[indexFactorA] - currentPrime;
      index2 = (rowPrimeSieveData->soln1 = (F2 += currentPrime & (F2 >> 31)));
      do
      {
        SieveArray[index2] += logPrimeEvenPoly;
      } while ((index2 += currentPrime) <= X1);
      F1 = F2 - (F3 = rowPrimeSieveData->Bainv2_0);
      F1 += currentPrime & (F1 >> 31);
      do
      {
        SieveArray[F1] += logPrimeOddPoly;
      } while ((F1 += currentPrime) <= X1);
      F2 -= rowPrimeSieveData->difsoln;
      F1 = F2 += currentPrime & (F2 >> 31);
      do
      {
        SieveArray[F2] += logPrimeEvenPoly;
      } while ((F2 += currentPrime) <= X1);
      F1 -= F3;
      F1 += currentPrime & (F1 >> 31);
      do
      {
        SieveArray[F1] += logPrimeOddPoly;
      } while ((F1 += currentPrime) <= X1);
    }
    for (; index < nbrPrimes2; index++)
    {
      rowPrimeSieveData = primeSieveData + index;
      currentPrime = rowPrimeSieveData->value;
      if (currentPrime >= mask)
      {
        mask *= 3;
        logPrimeEvenPoly++;
        logPrimeOddPoly += 0x100;
      }
      F2 = rowPrimeSieveData->soln1 +
        rowPrimeSieveData->Bainv2[indexFactorA] - currentPrime;
      if ((rowPrimeSieveData->soln1 = (F2 += currentPrime & (F2 >> 31))) < X1)
      {
        SieveArray[F2] += logPrimeEvenPoly;
      }
      F1 = F2 - (F3 = rowPrimeSieveData->Bainv2_0);
      if ((F1 += currentPrime & (F1 >> 31)) < X1)
      {
        SieveArray[F1] += logPrimeOddPoly;
      }
      F2 -= rowPrimeSieveData->difsoln;
      if ((F2 += currentPrime & (F2 >> 31)) < X1)
      {
        SieveArray[F2] += logPrimeEvenPoly;
      }
      F2 -= F3;
      if ((F2 += currentPrime & (F2 >> 31)) < X1)
      {
        SieveArray[F2] += logPrimeOddPoly;
      }
      rowPrimeSieveData = primeSieveData + ++index;
      currentPrime = rowPrimeSieveData->value;
      F2 = rowPrimeSieveData->soln1 +
        rowPrimeSieveData->Bainv2[indexFactorA] - currentPrime;
      if ((rowPrimeSieveData->soln1 = (F2 += currentPrime & (F2 >> 31))) < X1)
      {
        SieveArray[F2] += logPrimeEvenPoly;
      }
      F1 = F2 - (F3 = rowPrimeSieveData->Bainv2_0);
      if ((F1 += currentPrime & (F1 >> 31)) < X1)
      {
        SieveArray[F1] += logPrimeOddPoly;
      }
      F2 -= rowPrimeSieveData->difsoln;
      if ((F2 += currentPrime & (F2 >> 31)) < X1)
      {
        SieveArray[F2] += logPrimeEvenPoly;
      }
      F2 -= F3;
      if ((F2 += currentPrime & (F2 >> 31)) < X1)
      {
        SieveArray[F2] += logPrimeOddPoly;
      }
      rowPrimeSieveData = primeSieveData + ++index;
      currentPrime = rowPrimeSieveData->value;
      F2 = rowPrimeSieveData->soln1 +
        rowPrimeSieveData->Bainv2[indexFactorA] - currentPrime;
      if ((rowPrimeSieveData->soln1 = (F2 += currentPrime & (F2 >> 31))) < X1)
      {
        SieveArray[F2] += logPrimeEvenPoly;
      }
      F1 = F2 - (F3 = rowPrimeSieveData->Bainv2_0);
      if ((F1 += currentPrime & (F1 >> 31)) < X1)
      {
        SieveArray[F1] += logPrimeOddPoly;
      }
      F2 -= rowPrimeSieveData->difsoln;
      if ((F2 += currentPrime & (F2 >> 31)) < X1)
      {
        SieveArray[F2] += logPrimeEvenPoly;
      }
      F2 -= F3;
      if ((F2 += currentPrime & (F2 >> 31)) < X1)
      {
        SieveArray[F2] += logPrimeOddPoly;
      }
      rowPrimeSieveData = primeSieveData + ++index;
      currentPrime = rowPrimeSieveData->value;
      F2 = rowPrimeSieveData->soln1 +
        rowPrimeSieveData->Bainv2[indexFactorA] - currentPrime;
      if ((rowPrimeSieveData->soln1 = (F2 += currentPrime & (F2 >> 31))) < X1)
      {
        SieveArray[F2] += logPrimeEvenPoly;
      }
      F1 = F2 - (F3 = rowPrimeSieveData->Bainv2_0);
      if ((F1 += currentPrime & (F1 >> 31)) < X1)
      {
        SieveArray[F1] += logPrimeOddPoly;
      }
      F2 -= rowPrimeSieveData->difsoln;
      if ((F2 += currentPrime & (F2 >> 31)) < X1)
      {
        SieveArray[F2] += logPrimeEvenPoly;
      }
      F2 -= F3;
      if ((F2 += currentPrime & (F2 >> 31)) < X1)
      {
        SieveArray[F2] += logPrimeOddPoly;
      }
    }
    for (; index < nbrPrimes; index++)
    {
      rowPrimeSieveData = primeSieveData + index;
      currentPrime = rowPrimeSieveData->value;
      if (currentPrime >= mask)
      {
        mask *= 3;
        logPrimeEvenPoly++;
        logPrimeOddPoly += 0x100;
      }
      F2 = rowPrimeSieveData->soln1 +
        rowPrimeSieveData->Bainv2[indexFactorA] - currentPrime;
      if ((rowPrimeSieveData->soln1 = (F2 += currentPrime & (F2 >> 31))) < X1)
      {
        SieveArray[F2] += logPrimeEvenPoly;
      }
      F1 = F2 - (F3 = rowPrimeSieveData->Bainv2_0);
      if ((F1 += currentPrime & (F1 >> 31)) < X1)
      {
        SieveArray[F1] += logPrimeOddPoly;
      }
      F2 -= rowPrimeSieveData->difsoln;
      if ((F2 += currentPrime & (F2 >> 31)) < X1)
      {
        SieveArray[F2] += logPrimeEvenPoly;
      }
      F2 -= F3;
      if ((F2 += currentPrime & (F2 >> 31)) < X1)
      {
        SieveArray[F2] += logPrimeOddPoly;
      }
    }
  }
}

static int PerformTrialDivision(PrimeSieveData *primeSieveData,
  int rowMatrixBbeforeMerge[],
  int index2,
  int biDividend[], int rowSquares[],
  int NumberLengthDividend,
  unsigned char oddPolynomial)
{
  int biR0 = 0, biR1 = 0, biR2 = 0, biR3 = 0, biR4 = 0, biR5 = 0;
  int biR6 = 0;
  int nbrSquares = rowSquares[0];
  int Divisor;
  int divis;
  int index;
  int expParity;
  int nbrColumns = rowMatrixBbeforeMerge[0];
  PrimeSieveData *rowPrimeSieveData;
  PrimeTrialDivisionData *rowPrimeTrialDivisionData;
  switch (NumberLengthDividend)
  {
  case 7:
    biR6 = biDividend[6];
    /* no break */
  case 6:
    biR5 = biDividend[5];
    /* no break */
  case 5:
    biR4 = biDividend[4];
    /* no break */
  case 4:
    biR3 = biDividend[3];
    /* no break */
  case 3:
    biR2 = biDividend[2];
    /* no break */
  case 1:
  case 2:
    biR1 = biDividend[1];
    biR0 = biDividend[0];
  }
  expParity = 0;
  if (NumberLengthDividend <= 1)
  {
    for (index = 1; index < nbrPrimes; index++)
    {
      Divisor = primeTrialDivisionData[index].value;
      while (biR0 % Divisor == 0)
      {
        biR0 /= Divisor;
        expParity = 1 - expParity;
        if (expParity == 0)
        {
          rowSquares[nbrSquares++] = (int)Divisor;
        }
      }
      if (expParity != 0)
      {
        rowMatrixBbeforeMerge[nbrColumns++] = index;
        expParity = 0;
      }
    }
  }
  else
  {
    unsigned char mostSignificantLimbZero = FALSE;
    double dRem, dCurrentPrime;
    double dDivid, dLimbMult;
    int LimbMult = 1 << BITS_PER_INT_GROUP;
    int quot, Dividend;
    int nbr, iRem, Rem;
    int left, right, median;
    unsigned char fullRemainder;
    int indexFactorA = 0;
    int newFactorAIndex;
    unsigned char testFactorA = TRUE;
    newFactorAIndex = aindex[0];
    for (index = 1; testFactorA; index++)
    {
      fullRemainder = FALSE;
      if (index < 3)
      {
        fullRemainder = TRUE;
      }
      else if (index == newFactorAIndex)
      {
        fullRemainder = TRUE;
        if (++indexFactorA == nbrFactorsA)
        {
          testFactorA = FALSE;   // All factors of A were tested.
        }
        else
        {
          newFactorAIndex = aindex[indexFactorA];
        }
      }
      for (;;)
      {
        rowPrimeTrialDivisionData = &primeTrialDivisionData[index];
        if (fullRemainder == FALSE)
        {
          rowPrimeSieveData = primeSieveData + index;
          Divisor = rowPrimeSieveData->value;
          divis = (int)Divisor;
          if (oddPolynomial)
          {
            iRem = index2 - rowPrimeSieveData->soln1 +
              rowPrimeSieveData->Bainv2_0;
          }
          else
          {
            iRem = index2 - rowPrimeSieveData->soln1;
          }
          if (iRem >= divis)
          {
            if ((iRem -= divis) >= divis)
            {
              if ((iRem -= divis) >= divis)
              {
                iRem %= divis;
              }
            }
          }
          else
          {
            iRem += (iRem >> 31) & divis;
          }
          if (iRem != 0 && iRem != divis - rowPrimeSieveData->difsoln)
          {
            if (expParity != 0)
            {
              rowMatrixBbeforeMerge[nbrColumns++] = index;
              expParity = 0;
            }
            break;              // Process next prime.
          }
          fullRemainder = TRUE;
        }
        else
        {
          Divisor = rowPrimeTrialDivisionData->value;
          divis = (int)Divisor;
          switch (NumberLengthDividend)
          {
          case 7:
            dRem = (double)biR6*(double)rowPrimeTrialDivisionData->exp6 +
              (double)biR5*(double)rowPrimeTrialDivisionData->exp5 +
              (double)biR4*(double)rowPrimeTrialDivisionData->exp4 +
              (double)biR3*(double)rowPrimeTrialDivisionData->exp3 +
              (double)biR2*(double)rowPrimeTrialDivisionData->exp2 +
              (double)biR1*(double)rowPrimeTrialDivisionData->exp1 + biR0;
            break;
          case 6:
            dRem = (double)biR5*(double)rowPrimeTrialDivisionData->exp5 +
              (double)biR4*(double)rowPrimeTrialDivisionData->exp4 +
              (double)biR3*(double)rowPrimeTrialDivisionData->exp3 +
              (double)biR2*(double)rowPrimeTrialDivisionData->exp2 +
              (double)biR1*(double)rowPrimeTrialDivisionData->exp1 + biR0;
            break;
          case 5:
            dRem = (double)biR4*(double)rowPrimeTrialDivisionData->exp4 +
              (double)biR3*(double)rowPrimeTrialDivisionData->exp3 +
              (double)biR2*(double)rowPrimeTrialDivisionData->exp2 +
              (double)biR1*(double)rowPrimeTrialDivisionData->exp1 + biR0;
            break;
          case 4:
            dRem = (double)biR3*(double)rowPrimeTrialDivisionData->exp3 +
              (double)biR2*(double)rowPrimeTrialDivisionData->exp2 +
              (double)biR1*(double)rowPrimeTrialDivisionData->exp1 + biR0;
            break;
          case 3:
            dRem = (double)biR2*(double)rowPrimeTrialDivisionData->exp2 +
              (double)biR1*(double)rowPrimeTrialDivisionData->exp1 + biR0;
            break;
          default:
            dRem = (double)biR1*(double)rowPrimeTrialDivisionData->exp1 + biR0;
            break;
          }
          dCurrentPrime = (double)divis;
          Rem = (int)(dRem - floor(dRem / dCurrentPrime)*dCurrentPrime);
          if (Rem%divis != 0)
          {                     // Number is not a multiple of prime.
            if (expParity != 0)
            {
              rowMatrixBbeforeMerge[nbrColumns++] = index;
              expParity = 0;
            }
            break;              // Process next prime.
          }
        }
        expParity = 1 - expParity;
        if (expParity == 0)
        {
          rowSquares[nbrSquares++] = (int)Divisor;
        }
        Rem = 0;
        // Perform division
        dLimbMult = (double)(1U << BITS_PER_INT_GROUP);
        dCurrentPrime = (double)divis;
        switch (NumberLengthDividend)
        {
        case 7:     // {biR6 - biR0} <- {biR6 - biR0} / divis
          Dividend = biR6;
          Rem = Dividend - (biR6 = Dividend / Divisor) * Divisor;
          /* no break */
        case 6:     // {biR5 - biR0} <- {biR5 - biR0} / divis
          dDivid = (double)biR5 + (double)Rem*dLimbMult;
          quot = (int)(dDivid / dCurrentPrime);
          Rem = biR5 + Rem*LimbMult - quot * divis;
          biR5 = quot;
          /* no break */
        case 5:     // {biR4 - biR0} <- {biR4 - biR0} / divis
          dDivid = (double)biR4 + (double)Rem*dLimbMult;
          quot = (int)(dDivid / dCurrentPrime);
          Rem = biR4 + Rem*LimbMult - quot * divis;
          biR4 = quot;
          /* no break */
        case 4:     // {biR3 - biR0} <- {biR3 - biR0} / divis
          dDivid = (double)biR3 + (double)Rem*dLimbMult;
          quot = (int)(dDivid / dCurrentPrime);
          Rem = biR3 + Rem*LimbMult - quot * divis;
          biR3 = quot;
          /* no break */
        case 3:     // {biR2 - biR0} <- {biR2 - biR0} / divis
          dDivid = (double)biR2 + (double)Rem*dLimbMult;
          quot = (int)(dDivid / dCurrentPrime);
          Rem = biR2 + Rem*LimbMult - quot * divis;
          biR2 = quot;
          /* no break */
        case 2:     // {biR1 - biR0} <- {biR1 - biR0} / divis
          dDivid = (double)biR1 + (double)Rem*dLimbMult;
          quot = (int)(dDivid / dCurrentPrime);
          Rem = biR1 + Rem*LimbMult - quot * divis;
          biR1 = quot;
          dDivid = (double)biR0 + (double)Rem*dLimbMult;
          biR0 = (int)(dDivid / dCurrentPrime);
        }
        switch (NumberLengthDividend)
        {
        case 7:
          mostSignificantLimbZero = (biR6 == 0 && biR5 < HALF_INT_RANGE);
          break;
        case 6:
          mostSignificantLimbZero = (biR5 == 0 && biR4 < HALF_INT_RANGE);
          break;
        case 5:
          mostSignificantLimbZero = (biR4 == 0 && biR3 < HALF_INT_RANGE);
          break;
        case 4:
          mostSignificantLimbZero = (biR3 == 0 && biR2 < HALF_INT_RANGE);
          break;
        case 3:
          mostSignificantLimbZero = (biR2 == 0 && biR1 < HALF_INT_RANGE);
          break;
        case 2:    // Criteria is to fit in a double (52 bits).
          mostSignificantLimbZero = (biR1 < (1 << (52 - BITS_PER_INT_GROUP)));
          break;
        }
        if (mostSignificantLimbZero)
        {
          NumberLengthDividend--;
          if (NumberLengthDividend == 1)
          {     // Number fits in a double
            double dDivid = (double)biR1 * (double)(1U << BITS_PER_INT_GROUP) + (double)biR0;
            int sqrtDivid = (int)(floor(sqrt(dDivid)));
            if (trialDivisions == 12135)
            {
              Divisor = 0;
            }
            fullRemainder = TRUE;
            for (; index < nbrPrimes; index++)
            {
              rowPrimeSieveData = primeSieveData + index;
              Divisor = rowPrimeSieveData->value;
              if (testFactorA && index == newFactorAIndex)
              {
                fullRemainder = TRUE;
                if (++indexFactorA == nbrFactorsA)
                {
                  testFactorA = FALSE;   // All factors of A were tested.
                }
                else
                {
                  newFactorAIndex = aindex[indexFactorA];
                }
              }
              for (;;)
              {
                if (fullRemainder == FALSE)
                {
                  divis = (int)Divisor;
                  if (oddPolynomial)
                  {
                    iRem = index2 - rowPrimeSieveData->soln1 +
                      rowPrimeSieveData->Bainv2_0;
                  }
                  else
                  {
                    iRem = index2 - rowPrimeSieveData->soln1;
                  }
                  if (iRem >= divis)
                  {
                    if ((iRem -= divis) >= divis)
                    {
                      if ((iRem -= divis) >= divis)
                      {
                        iRem %= divis;
                      }
                    }
                  }
                  else
                  {
                    iRem += (iRem >> 31) & divis;
                  }
                  if (iRem != 0 && iRem != divis - rowPrimeSieveData->difsoln)
                  {
                    break;
                  }
                  fullRemainder = TRUE;
                }
                else
                {
                  double dQuot = floor(dDivid / Divisor);
                  if (dDivid - dQuot * Divisor != 0)
                  {
                    break;
                  }
                }
                dDivid /= Divisor;
                sqrtDivid = (int)(floor(sqrt(dDivid))) + 1;
                expParity = 1 - expParity;
                if (expParity == 0)
                {
                  rowSquares[nbrSquares++] = (int)Divisor;
                }
              }
              if (expParity != 0)
              {
                rowMatrixBbeforeMerge[nbrColumns++] = index;
                expParity = 0;
              }
              if (Divisor > sqrtDivid)
              {                     // End of trial division.
                rowSquares[0] = nbrSquares;
                index = nbrPrimes - 1;
                if (dDivid <= primeTrialDivisionData[index].value &&
                  dDivid > 1)
                {          // Perform binary search to find the index.
                  left = -1;
                  median = right = nbrPrimes;
                  while (left != right)
                  {
                    median = ((right - left) >> 1) + left;
                    nbr = primeTrialDivisionData[median].value;
                    if (nbr < dDivid)
                    {
                      if (median == left &&
                        congruencesFound >= matrixBLength)
                      {
                        return 0;
                      }
                      left = median;
                    }
                    else if (nbr > dDivid)
                    {
                      right = median;
                    }
                    else
                    {
                      break;
                    }
                  }
                  rowMatrixBbeforeMerge[nbrColumns++] = median;
                  rowMatrixBbeforeMerge[0] = nbrColumns;
                  return 1;
                }
                rowMatrixBbeforeMerge[0] = nbrColumns;
                return (int)dDivid;
              }
              fullRemainder = FALSE;
            }
            break;
          }
        }
      }             /* end while */
    }               /* end for */
    for (; index < nbrPrimes; index++)
    {
      fullRemainder = FALSE;
      for (;;)
      {
        rowPrimeSieveData = primeSieveData + index;
        rowPrimeTrialDivisionData = &primeTrialDivisionData[index];
        if (fullRemainder == FALSE)
        {
          Divisor = rowPrimeSieveData->value;
          divis = (int)Divisor;
          if (trialDivisions == 12135 && Divisor == 2857)
          {
            Rem = 0;
          }
          if (oddPolynomial)
          {
            iRem = index2 - rowPrimeSieveData->soln1 +
              rowPrimeSieveData->Bainv2_0;
          }
          else
          {
            iRem = index2 - rowPrimeSieveData->soln1;
          }
          if (iRem >= divis)
          {
            if ((iRem -= divis) >= divis)
            {
              if ((iRem -= divis) >= divis)
              {
                iRem %= divis;
              }
            }
          }
          else
          {
            iRem += (iRem >> 31) & divis;
          }
          if (iRem != 0 && iRem != divis - rowPrimeSieveData->difsoln)
          {
            if (expParity != 0)
            {
              rowMatrixBbeforeMerge[nbrColumns++] = index;
              expParity = 0;
            }
            break;              // Process next prime.
          }
          fullRemainder = TRUE;
        }
        else
        {
          Divisor = rowPrimeTrialDivisionData->value;
          divis = (int)Divisor;
          switch (NumberLengthDividend)
          {
          case 7:
            dRem = (double)biR6*(double)rowPrimeTrialDivisionData->exp6 +
              (double)biR5*(double)rowPrimeTrialDivisionData->exp5 +
              (double)biR4*(double)rowPrimeTrialDivisionData->exp4 +
              (double)biR3*(double)rowPrimeTrialDivisionData->exp3 +
              (double)biR2*(double)rowPrimeTrialDivisionData->exp2 +
              (double)biR1*(double)rowPrimeTrialDivisionData->exp1 + biR0;
            break;
          case 6:
            dRem = (double)biR5*(double)rowPrimeTrialDivisionData->exp5 +
              (double)biR4*(double)rowPrimeTrialDivisionData->exp4 +
              (double)biR3*(double)rowPrimeTrialDivisionData->exp3 +
              (double)biR2*(double)rowPrimeTrialDivisionData->exp2 +
              (double)biR1*(double)rowPrimeTrialDivisionData->exp1 + biR0;
            break;
          case 5:
            dRem = (double)biR4*(double)rowPrimeTrialDivisionData->exp4 +
              (double)biR3*(double)rowPrimeTrialDivisionData->exp3 +
              (double)biR2*(double)rowPrimeTrialDivisionData->exp2 +
              (double)biR1*(double)rowPrimeTrialDivisionData->exp1 + biR0;
            break;
          case 4:
            dRem = (double)biR3*(double)rowPrimeTrialDivisionData->exp3 +
              (double)biR2*(double)rowPrimeTrialDivisionData->exp2 +
              (double)biR1*(double)rowPrimeTrialDivisionData->exp1 + biR0;
            break;
          case 3:
            dRem = (double)biR2*(double)rowPrimeTrialDivisionData->exp2 +
              (double)biR1*(double)rowPrimeTrialDivisionData->exp1 + biR0;
            break;
          default:
            dRem = (double)biR1*(double)rowPrimeTrialDivisionData->exp1 + biR0;
            break;
          }
          dCurrentPrime = (double)divis;
          Rem = (int)(dRem - floor(dRem / dCurrentPrime)*dCurrentPrime);
          if (Rem%divis != 0)
          {                     // Number is not a multiple of prime.
            if (expParity != 0)
            {
              rowMatrixBbeforeMerge[nbrColumns++] = index;
              expParity = 0;
            }
            break;              // Process next prime.
          }
        }
        expParity = 1 - expParity;
        if (expParity == 0)
        {
          rowSquares[nbrSquares++] = (int)Divisor;
        }
        Rem = 0;
        // Perform division
        dLimbMult = (double)(1U << BITS_PER_INT_GROUP);
        dCurrentPrime = (double)divis;
        switch (NumberLengthDividend)
        {
        case 7:     // {biR6 - biR0} <- {biR6 - biR0} / divis
          Dividend = biR6;
          Rem = Dividend - (biR6 = Dividend / Divisor) * Divisor;
          /* no break */
        case 6:     // {biR5 - biR0} <- {biR5 - biR0} / divis
          dDivid = (double)biR5 + (double)Rem*dLimbMult;
          quot = (int)(dDivid / dCurrentPrime);
          Rem = biR5 + Rem*LimbMult - quot * divis;
          biR5 = quot;
          /* no break */
        case 5:     // {biR4 - biR0} <- {biR4 - biR0} / divis
          dDivid = (double)biR4 + (double)Rem*dLimbMult;
          quot = (int)(dDivid / dCurrentPrime);
          Rem = biR4 + Rem*LimbMult - quot * divis;
          biR4 = quot;
          /* no break */
        case 4:     // {biR3 - biR0} <- {biR3 - biR0} / divis
          dDivid = (double)biR3 + (double)Rem*dLimbMult;
          quot = (int)(dDivid / dCurrentPrime);
          Rem = biR3 + Rem*LimbMult - quot * divis;
          biR3 = quot;
          /* no break */
        case 3:     // {biR2 - biR0} <- {biR2 - biR0} / divis
          dDivid = (double)biR2 + (double)Rem*dLimbMult;
          quot = (int)(dDivid / dCurrentPrime);
          Rem = biR2 + Rem*LimbMult - quot * divis;
          biR2 = quot;
          /* no break */
        case 2:     // {biR1 - biR0} <- {biR1 - biR0} / divis
          dDivid = (double)biR1 + (double)Rem*dLimbMult;
          quot = (int)(dDivid / dCurrentPrime);
          Rem = biR1 + Rem*LimbMult - quot * divis;
          biR1 = quot;
          dDivid = (double)biR0 + (double)Rem*dLimbMult;
          biR0 = (int)(dDivid / dCurrentPrime);
        }
        switch (NumberLengthDividend)
        {
        case 7:
          mostSignificantLimbZero = (biR6 == 0 && biR5 < HALF_INT_RANGE);
          break;
        case 6:
          mostSignificantLimbZero = (biR5 == 0 && biR4 < HALF_INT_RANGE);
          break;
        case 5:
          mostSignificantLimbZero = (biR4 == 0 && biR3 < HALF_INT_RANGE);
          break;
        case 4:
          mostSignificantLimbZero = (biR3 == 0 && biR2 < HALF_INT_RANGE);
          break;
        case 3:
          mostSignificantLimbZero = (biR2 == 0 && biR1 < HALF_INT_RANGE);
          break;
        case 2:    // Criteria is to fit in a double (52 bits).
          mostSignificantLimbZero = (biR1 < (1 << (52 - BITS_PER_INT_GROUP)));
          break;
        }
        if (mostSignificantLimbZero)
        {
          NumberLengthDividend--;
          if (NumberLengthDividend == 1)
          {
            double dDivid = (double)biR1 * (double)(1U << BITS_PER_INT_GROUP) + (double)biR0;
            int sqrtDivid = (int)(floor(sqrt(dDivid)));
            fullRemainder = TRUE;
            for (; index < nbrPrimes; index++)
            {
              rowPrimeSieveData = primeSieveData + index;
              if (rowPrimeSieveData->value == 41893)
              {
                divis = 5;
              }
              Divisor = rowPrimeSieveData->value;
              if (testFactorA && index == newFactorAIndex)
              {
                fullRemainder = TRUE;
                if (++indexFactorA == nbrFactorsA)
                {
                  testFactorA = FALSE;   // All factors of A were tested.
                }
                else
                {
                  newFactorAIndex = aindex[indexFactorA];
                }
              }
              for (;;)
              {
                if (fullRemainder == FALSE)
                {
                  divis = (int)Divisor;
                  if (oddPolynomial)
                  {
                    iRem = index2 - rowPrimeSieveData->soln1 +
                      rowPrimeSieveData->Bainv2_0;
                  }
                  else
                  {
                    iRem = index2 - rowPrimeSieveData->soln1;
                  }
                  if (iRem >= divis)
                  {
                    if ((iRem -= divis) >= divis)
                    {
                      if ((iRem -= divis) >= divis)
                      {
                        iRem %= divis;
                      }
                    }
                  }
                  else
                  {
                    iRem += (iRem >> 31) & divis;
                  }
                  if (iRem != 0 && iRem != divis - rowPrimeSieveData->difsoln)
                  {
                    break;
                  }
                  fullRemainder = TRUE;
                }
                else
                {
                  double dQuot = floor(dDivid / Divisor);
                  if (dDivid - dQuot * Divisor != 0)
                  {
                    break;
                  }
                }
                dDivid /= Divisor;
                sqrtDivid = (int)floor(sqrt(dDivid));
                expParity = 1 - expParity;
                if (expParity == 0)
                {
                  rowSquares[nbrSquares++] = (int)Divisor;
                }
              }
              if (expParity != 0)
              {
                rowMatrixBbeforeMerge[nbrColumns++] = index;
                expParity = 0;
              }
              if (Divisor > sqrtDivid)
              {                     // End of trial division.
                if (dDivid >= (double)(1U << BITS_PER_INT_GROUP))
                {                   // Dividend is too large.
                  return 0;
                }
                Divisor = (int)dDivid;
                rowSquares[0] = nbrSquares;
                index = nbrPrimes - 1;
                if (Divisor <= primeTrialDivisionData[index].value &&
                  Divisor > 1)
                {          // Perform binary search to find the index.
                  left = -1;
                  median = right = nbrPrimes;
                  while (left != right)
                  {
                    median = ((right - left) >> 1) + left;
                    nbr = primeTrialDivisionData[median].value;
                    if (nbr < Divisor)
                    {
                      if (median == left &&
                        congruencesFound >= matrixBLength)
                      {
                        return 0;
                      }
                      left = median;
                    }
                    else if (nbr > Divisor)
                    {
                      right = median;
                    }
                    else
                    {
                      break;
                    }
                  }
                  rowMatrixBbeforeMerge[nbrColumns++] = median;
                  rowMatrixBbeforeMerge[0] = nbrColumns;
                  return 1;
                }
                rowMatrixBbeforeMerge[0] = nbrColumns;
                return Divisor;
              }
              fullRemainder = FALSE;
            }
            if (dDivid >= (double)(1U << BITS_PER_INT_GROUP))
            {                   // Dividend is too large.
              return 0;
            }
            biR0 = (int)dDivid;
            break;
          }
        }
      }             /* end while */
    }               /* end for */
  }
  rowSquares[0] = nbrSquares;
  rowMatrixBbeforeMerge[0] = nbrColumns;
  if (NumberLengthDividend > 1)
  {
    return 0;           // Very large quotient.
  }
  return biR0;
}

static void mergeArrays(int aindex[], int nbrFactorsA, int rowMatrixB[], int rowMatrixBeforeMerge[],
  int rowSquares[])
{
  int indexAindex = 0;
  int indexRMBBM = 1;
  int indexRMB = 1;
  int nbrColumns = rowMatrixBeforeMerge[0];

  while (indexAindex < nbrFactorsA && indexRMBBM < nbrColumns)
  {
    if (aindex[indexAindex] < rowMatrixBeforeMerge[indexRMBBM])
    {
      rowMatrixB[indexRMB++] = aindex[indexAindex++];
    }
    else if (aindex[indexAindex] > rowMatrixBeforeMerge[indexRMBBM])
    {
      rowMatrixB[indexRMB++] = rowMatrixBeforeMerge[indexRMBBM++];
    }
    else
    {
      rowSquares[rowSquares[0]++] =
        primeTrialDivisionData[aindex[indexAindex++]].value;
      indexRMBBM++;
    }
  }
  while (indexAindex < nbrFactorsA)
  {
    rowMatrixB[indexRMB++] = aindex[indexAindex++];
  }
  while (indexRMBBM < nbrColumns)
  {
    rowMatrixB[indexRMB++] = rowMatrixBeforeMerge[indexRMBBM++];
  }
  rowMatrixB[LENGTH_OFFSET] = indexRMB;
}

static void SmoothRelationFound(
  unsigned char positive,
  int *rowMatrixB, int *rowMatrixBbeforeMerge,
  int index2,
  int *rowSquares, int *biLinearCoeff,
  int NumberLength, int *biT, int *biU,
  int *biR, unsigned char oddPolynomial)
{
  int index;
  int nbrSquares;
  if (congruencesFound == matrixBLength)
  {
    return;            // All congruences already found.
  }
  // Add all elements of aindex array to the rowMatrixB array discarding
  // duplicates.
  mergeArrays(aindex, nbrFactorsA, rowMatrixB, rowMatrixBbeforeMerge, rowSquares);
  nbrSquares = rowSquares[0];
  IntToBigNbr(1, biR, NumberLength);
  IntToBigNbr(positive ? 1 : -1, biT, NumberLength);
  MultBigNbrByInt(biQuadrCoeff, index2 - SieveLimit, biU,
    NumberLength);                     // Ax
  AddBigNbr(biU, biLinearCoeff, biU, NumberLength);         // Ax+B
  if (oddPolynomial)
  {
    SubtractBigNbr(biU, biLinearDelta[0], biU, NumberLength);// Ax+B (odd)
    SubtractBigNbr(biU, biLinearDelta[0], biU, NumberLength);// Ax+B (odd)
  }
  if ((biU[NumberLength - 1] & HALF_INT_RANGE) != 0)
  {                                        // If number is negative
    ChSignBigNbr(biU, NumberLength);   // make it positive.
  }
  for (index = 1; index < nbrSquares; index++)
  {
    int D = rowSquares[index];
    if (D == multiplier)
    {
      AddBigNbr(biU, Modulus, biU, NumberLength);
      DivBigNbrByInt(biU, D, biU, NumberLength);
    }
    else
    {
      MultBigNbrByInt(biR, D, biR, NumberLength);
    }
  }
  if (InsertNewRelation(rowMatrixB, biT, biU, biR, NumberLength))
  {
    smoothsFound++;
    ShowSIQSStatus();
  }
  return;
}

static void PartialRelationFound(
  unsigned char positive,
  int *rowMatrixB, int *rowMatrixBbeforeMerge,
  int index2,
  int Divid, int *rowPartials,
  int *rowSquares,
  int *biLinearCoeff, int NumberLength, int *biT,
  int *biR, int *biU, int *biV,
  int *indexFactorsA, unsigned char oddPolynomial)
{
  int index;
  int expParity;
  int D, Divisor;
  int nbrFactorsPartial;
  int prev;
  unsigned int seed;
  int hashIndex;
  int *rowPartial;
  int newDivid = (int)Divid;    // This number is greater than zero.
  int indexFactorA = 0;
  int nbrSquares;
  int NumberLengthDivid;
  int biT0, biT1, biT2, biT3, biT4, biT5, biT6;
  int squareRootSize = NumberLength / 2 + 1;
  int nbrColumns;
  PrimeTrialDivisionData *rowPrimeTrialDivisionData;

  if (congruencesFound == matrixBLength)
  {
    return;
  }
  // Partial relation found.
  totalPartials++;
  // Check if there is already another relation with the same
  // factor outside the prime base.
  // Calculate hash index
  hashIndex = matrixPartialHashIndex[(int)(Divid & 0xFFE) >> 1];
  prev = -1;
  while (hashIndex >= 0)
  {
    int oldDivid;

    rowPartial = matrixPartial[hashIndex];
    oldDivid = rowPartial[0];
    if (newDivid == oldDivid || newDivid == -oldDivid)
    {   // Match of partials.
      double dRem, dDivisor;
      for (index = 0; index < squareRootSize; index++)
      {
        biV[index] = rowPartial[index + 2];
      }                           // biV = Old positive square root (Ax+B).
      for (; index < NumberLength; index++)
      {
        biV[index] = 0;
      }
      seed = rowPartial[squareRootSize + 2];
      getFactorsOfA(seed, indexFactorsA);
      IntToBigNbr(newDivid, biR, NumberLength);
      nbrFactorsPartial = 0;
      // biT = old (Ax+B)^2.
      MultBigNbr(biV, biV, biT, NumberLength);
      // biT = old (Ax+B)^2 - N.
      SubtractBigNbr(biT, Modulus, biT, NumberLength);
      if (oldDivid < 0)
      {
        rowPartials[nbrFactorsPartial++] = 0; // Insert -1 as a factor.
      }
      if (biT[NumberLength - 1] >= HALF_INT_RANGE)
      {
        ChSignBigNbr(biT, NumberLength);   // Make it positive.
      }
      NumberLengthDivid = NumberLength;
      // The number is multiple of the big prime, so divide by it.
      DivBigNbrByInt(biT, newDivid, biT, NumberLengthDivid);
      if (biT[NumberLengthDivid - 1] == 0 &&
        biT[NumberLengthDivid - 2] < HALF_INT_RANGE)
      {
        NumberLengthDivid--;
      }
      for (index = 0; index < nbrFactorsA; index++)
      {
        DivBigNbrByInt(biT,
          primeTrialDivisionData[indexFactorsA[index]].value, biT,
          NumberLengthDivid);
        if (biT[NumberLengthDivid - 1] == 0 &&
          biT[NumberLengthDivid - 2] < HALF_INT_RANGE)
        {
          NumberLengthDivid--;
        }
      }
      biT0 = biT[0];
      biT1 = biT[1];
      biT2 = biT[2];
      biT3 = biT[3];
      biT4 = biT[4];
      biT5 = biT[5];
      biT6 = biT[6];
      for (index = 1; index < nbrPrimes; index++)
      {
        expParity = 0;
        if (index >= indexMinFactorA && indexFactorA < nbrFactorsA)
        {
          if (index == indexFactorsA[indexFactorA])
          {
            expParity = 1;
            indexFactorA++;
          }
        }
        rowPrimeTrialDivisionData = &primeTrialDivisionData[index];
        Divisor = rowPrimeTrialDivisionData->value;
        for (;;)
        {
          switch (NumberLengthDivid)
          {
          case 7:
            dRem = (double)biT6*(double)rowPrimeTrialDivisionData->exp6 +
              (double)biT5*(double)rowPrimeTrialDivisionData->exp5 +
              (double)biT4*(double)rowPrimeTrialDivisionData->exp4 +
              (double)biT3*(double)rowPrimeTrialDivisionData->exp3 +
              (double)biT2*(double)rowPrimeTrialDivisionData->exp2;
            break;
          case 6:
            dRem = (double)biT5*(double)rowPrimeTrialDivisionData->exp5 +
              (double)biT4*(double)rowPrimeTrialDivisionData->exp4 +
              (double)biT3*(double)rowPrimeTrialDivisionData->exp3 +
              (double)biT2*(double)rowPrimeTrialDivisionData->exp2;
            break;
          case 5:
            dRem = (double)biT4*(double)rowPrimeTrialDivisionData->exp4 +
              (double)biT3*(double)rowPrimeTrialDivisionData->exp3 +
              (double)biT2*(double)rowPrimeTrialDivisionData->exp2;
            break;
          case 4:
            dRem = (double)biT3*(double)rowPrimeTrialDivisionData->exp3 +
              (double)biT2*(double)rowPrimeTrialDivisionData->exp2;
            break;
          case 3:
            dRem = (double)biT2*(double)rowPrimeTrialDivisionData->exp2;
            break;
          default:
            dRem = 0;
            break;
          }
          dRem += (double)biT1*(double)rowPrimeTrialDivisionData->exp1 + biT0;
          dDivisor = (double)Divisor;
          dRem -= floor(dRem / dDivisor)*dDivisor;
          if (dRem != 0)
          {
            break;
          }
          expParity = 1 - expParity;
          DivBigNbrByInt(biT, Divisor, biT, NumberLengthDivid);
          biT0 = biT[0];
          biT1 = biT[1];
          biT2 = biT[2];
          biT3 = biT[3];
          biT4 = biT[4];
          biT5 = biT[5];
          biT6 = biT[6];
          if (expParity == 0)
          {
            rowSquares[rowSquares[0]++] = (int)Divisor;
          }
          if (NumberLengthDivid <= 2)
          {
            if (biT0 == 1 && biT1 == 0)
            {               // biT = 1, so division has ended.
              break;
            }
          }
          else if (biT[NumberLengthDivid - 1] == 0 &&
            biT[NumberLengthDivid - 2] < HALF_INT_RANGE)
          {
            NumberLengthDivid--;
          }
        }
        if (expParity != 0)
        {
          rowPartials[nbrFactorsPartial++] = index;
        }
      }
      MultBigNbrByInt(biQuadrCoeff, index2 - SieveLimit, biT,
        NumberLength);
      AddBigNbr(biT, biLinearCoeff, biT, NumberLength); // biT = Ax+B
      if (oddPolynomial)
      {                                                     // Ax+B (odd)
        SubtractBigNbr(biT, biLinearDelta[0], biT, NumberLength);
        SubtractBigNbr(biT, biLinearDelta[0], biT, NumberLength);
      }
      if ((biT[NumberLength - 1] & HALF_INT_RANGE) != 0)
      {                                        // If number is negative
        ChSignBigNbr(biT, NumberLength);   // make it positive.
      }
      // biU = Product of old Ax+B times new Ax+B
      MultBigNbrModN(biV, biT, biU, Modulus, NumberLength);
      // Add all elements of aindex array to the rowMatrixB array discarding
      // duplicates.
      mergeArrays(aindex, nbrFactorsA, rowMatrixB, rowMatrixBbeforeMerge, rowSquares);
      rowMatrixBbeforeMerge[0] = nbrColumns = rowMatrixB[LENGTH_OFFSET];
      memcpy(&rowMatrixBbeforeMerge[1], &rowMatrixB[1], nbrColumns*sizeof(int));
      mergeArrays(rowPartials, nbrFactorsPartial, rowMatrixB, rowMatrixBbeforeMerge, rowSquares);
      nbrSquares = rowSquares[0];
      for (index = 1; index < nbrSquares; index++)
      {
        D = rowSquares[index];
        if (D != multiplier)
        {
          MultBigNbrByInt(biR, D, biR, NumberLength);
        }
        else
        {
          AddBigNbr(biU, Modulus, biU, NumberLength);
          DivBigNbrByInt(biU, multiplier, biU, NumberLength);
        }
      }
      if (rowMatrixB[0] > 1 &&
        InsertNewRelation(rowMatrixB, biT, biU, biR, NumberLength))
      {
        partialsFound++;
        ShowSIQSStatus();
      }
      return;
    }
    else
    {
      prev = hashIndex;
      hashIndex = rowPartial[1]; // Get next index for same hash.
    }
  } /* end while */
//  synchronized(firstPrimeSieveData)
  {
    if (hashIndex == -1 && nbrPartials < matrixPartialLength)
    { // No match and partials table is not full.
      // Add partial to table of partials.
      if (prev >= 0)
      {
        matrixPartial[prev][1] = nbrPartials;
      }
      else
      {
        matrixPartialHashIndex[(newDivid & 0xFFE) >> 1] = nbrPartials;
      }
      rowPartial = matrixPartial[nbrPartials];
      // Add all elements of aindex array to the rowMatrixB array discarding
      // duplicates.
      mergeArrays(aindex, nbrFactorsA, rowMatrixB, rowMatrixBbeforeMerge, rowSquares);
      IntToBigNbr(Divid, biR, NumberLength);
      nbrSquares = rowSquares[0];
      for (index = 1; index < nbrSquares; index++)
      {
        D = rowSquares[index];
        MultBigNbrByIntModN(biR, D, biR, Modulus, NumberLength);
        if (D == multiplier)
        {
          DivBigNbrByInt(biU, D, biU, NumberLength);
        }
      }
      rowPartial[0] = (positive ? newDivid : -newDivid);
      // Indicate last index with this hash.
      rowPartial[1] = -1;
      MultBigNbrByInt(biQuadrCoeff, index2 - SieveLimit, biT,
        NumberLength);
      AddBigNbr(biT, biLinearCoeff, biT, NumberLength); // biT = Ax+B
      if (oddPolynomial)
      {                                                     // Ax+B (odd)
        SubtractBigNbr(biT, biLinearDelta[0], biT, NumberLength);
        SubtractBigNbr(biT, biLinearDelta[0], biT, NumberLength);
      }
      if ((biT[NumberLength - 1] & HALF_INT_RANGE) != 0)
      {                      // If square root is negative convert to positive.
        ChSignBigNbr(biT, NumberLength);
      }
      for (index = 0; index < squareRootSize; index++)
      {
        rowPartial[index + 2] = (int)biT[index];
      }
      rowPartial[squareRootSize + 2] = (int)oldSeed;
      nbrPartials++;
    }
  }               // End synchronized block.
  return;
}

static void SieveLocationHit(int rowMatrixB[], int rowMatrixBbeforeMerge[],
  int index2,
  PrimeSieveData *primeSieveData,
  int rowPartials[],
  int rowSquares[], int biDividend[],
  int NumberLength, int biT[], int *biLinearCoeff,
  int biR[], int biU[], int biV[],
  int indexFactorsA[], unsigned char oddPolynomial)
{
  unsigned char positive;
  int NumberLengthDivid;
  int index;
  int Divid;
  int nbrColumns;

  trialDivisions++;
  MultBigNbrByInt(biQuadrCoeff, index2 - SieveLimit, biT,
    NumberLength);                                      // Ax
  AddBigNbr(biT, biLinearCoeff, biT, NumberLength);     // Ax+B
  if (oddPolynomial)
  {                                                     // Ax+B (odd)
    SubtractBigNbr(biT, biLinearDelta[0], biT, NumberLength);
    SubtractBigNbr(biT, biLinearDelta[0], biT, NumberLength);
  }
  MultBigNbr(biT, biT, biDividend, NumberLength);       // (Ax+B)^2
                                                        // To factor: (Ax+B)^2-N
  SubtractBigNbr(biDividend, Modulus, biDividend, NumberLength);
  /* factor biDividend */

  NumberLengthDivid = NumberLength; /* Number length for dividend */
  positive = TRUE;
  if (biDividend[NumberLengthDivid - 1] >= HALF_INT_RANGE)
  { /* Negative */
    positive = FALSE;
    ChSignBigNbr(biDividend, NumberLengthDivid); // Convert to positive
  }
  rowSquares[0] = 1;
  for (index = 0; index < nbrFactorsA; index++)
  {
    DivBigNbrByInt(biDividend, afact[index], biDividend,
      NumberLengthDivid);
    if ((biDividend[NumberLengthDivid - 1] == 0
      && biDividend[NumberLengthDivid - 2] < HALF_INT_RANGE))
    {
      NumberLengthDivid--;
    }
  }
  nbrColumns = 1;
  if (!positive)
  {                                  // Insert -1 as a factor.
    rowMatrixBbeforeMerge[nbrColumns++] = 0;
  }
  rowMatrixBbeforeMerge[0] = nbrColumns;
  Divid = PerformTrialDivision(primeSieveData,
    rowMatrixBbeforeMerge,
    index2, biDividend,
    rowSquares, NumberLengthDivid,
    oddPolynomial);
  if (Divid == 1)
  { // Smooth relation found.
    SmoothRelationFound(positive, rowMatrixB,
      rowMatrixBbeforeMerge,
      index2,
      rowSquares,
      biLinearCoeff, NumberLength, biT, biU, biR,
      oddPolynomial);
  }
  else
  {
    if (Divid > 0 && Divid < largePrimeUpperBound)
    {
      PartialRelationFound(positive, rowMatrixB,
        rowMatrixBbeforeMerge,
        index2,
        Divid, rowPartials,
        rowSquares, biLinearCoeff,
        NumberLength, biT, biR, biU, biV,
        indexFactorsA, oddPolynomial);
    }
  }
  return;
}

static unsigned int getFactorsOfA(unsigned int seed, int *indexA)
{
  int index, index2, i, tmp;
  for (index = 0; index < nbrFactorsA; index++)
  {
    do
    {
      seed = 1141592621 * seed + 321435;
      i = (int)(((double)seed * (double)span)/(double)0x100000000ll + indexMinFactorA);
      for (index2 = 0; index2 < index; index2++)
      {
        if (indexA[index2] == i || indexA[index2] == i + 1)
        {
          break;
        }
      }
    } while (index2 < index);
    indexA[index] = i;
  }
  for (index = 0; index<nbrFactorsA; index++)    // Sort factors of A.
  {
    for (index2 = index + 1; index2<nbrFactorsA; index2++)
    {
      if (indexA[index] > indexA[index2])
      {
        tmp = indexA[index];
        indexA[index] = indexA[index2];
        indexA[index2] = tmp;
      }
    }
  }
  return seed;
}

/************************************************************************/
/* Multithread procedure:                                               */
/*                                                                      */
/* 1) Main thread generates factor base and other parameters.           */
/* 2) Start M threads where the number M is specified by the user in a  */
/*    box beneath the applet.                                           */
/* 3) For each polynomial:                                              */
/*    3a) Main thread generates the data for the set of 2^n polynomials.*/
/*    3b) Each child thread computes a range of polynomials             */
/*        (u*2^n/M to (u+1)*2^n/M exclusive).                           */
/* Partial and full relation routines must be synchronized.             */
/************************************************************************/
void FactoringSIQS(limb *pNbrToFactor, limb *pFactor)
{
  int origNumberLength;
  int FactorBase;
  int currentPrime;
  int NbrMod;
  BigInteger TempResult;
  PrimeSieveData *rowPrimeSieveData;
  PrimeTrialDivisionData *rowPrimeTrialDivisionData;
  int Power2, SqrRootMod, fact;
  int D, E, Q, V, W, X, Y, Z, T1, V1, W1, Y1;
  double Temp, Prod;
  double bestadjust;
  int i, j;
  int arrmult[] = { 1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
    47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97 };
  double adjustment[sizeof(arrmult)/sizeof(arrmult[0])];
  double dNumberToFactor, dlogNumberToFactor;
  origNumberLength = NumberLength;
  nbrThreadFinishedPolySet = 0;
  trialDivisions = 0;
  smoothsFound = 0;
  totalPartials = 0;
  partialsFound = 0;
  ValuesSieved = 0;
  congruencesFound = 0;
  polynomialsSieved = 0;
  nbrPartials = 0;
  matrixPartialLength = 0;
  newSeed = 0;

//  threadArray = new Thread[numberThreads];
  Temp = logLimbs(pNbrToFactor, origNumberLength);
  nbrPrimes = (int)exp(sqrt(Temp * log(Temp)) * 0.318);
  SieveLimit = (int)exp(8.5 + 0.015 * Temp) & 0xFFFFFFF8;
  nbrFactorsA = (int)(Temp*0.051 + 1);
  NbrPolynomials = (1 << (nbrFactorsA - 1)) - 1;
  factorSiqs.nbrLimbs = NumberLength;
  factorSiqs.sign = SIGN_POSITIVE;
  memcpy(factorSiqs.limbs, pNbrToFactor, NumberLength * sizeof(limb));
  NumberLength = BigNbrToBigInt(&factorSiqs, Modulus);
  Modulus[NumberLength++] = 0;
  memcpy(TestNbr2, Modulus, NumberLength*sizeof(int));
  memset(matrixPartialHashIndex, 0xFF, sizeof(matrixPartialHashIndex));
#ifdef __EMSCRIPTEN__
  InitSIQSStrings(SieveLimit);
  startSieveTenths = (int)(tenths() - originalTenthSecond);
#endif
  /************************/
  /* Compute startup data */
  /************************/

  /* search for best Knuth-Schroeppel multiplier */
  bestadjust = -10.0e0;
  primeSieveData[0].value = 1;
  primeTrialDivisionData[0].value = 1;
  rowPrimeSieveData = &primeSieveData[1];
  rowPrimeTrialDivisionData = &primeTrialDivisionData[1];
  rowPrimeSieveData->value = 2;
  rowPrimeTrialDivisionData->value = 2;
  // (2^31)^(j+1) mod 2
  rowPrimeTrialDivisionData->exp1 = rowPrimeTrialDivisionData->exp2 =
    rowPrimeTrialDivisionData->exp3 = rowPrimeTrialDivisionData->exp4 =
    rowPrimeTrialDivisionData->exp5 = rowPrimeTrialDivisionData->exp6 = 0;

  NbrMod = pNbrToFactor->x & 7;
  for (j = 0; j<sizeof(arrmult)/sizeof(arrmult[0]); j++)
  {
    int mod = (NbrMod * arrmult[j]) & 7;
    adjustment[j] = 0.34657359; /*  (ln 2)/2  */
    if (mod == 1)
      adjustment[j] *= (4.0e0);
    if (mod == 5)
      adjustment[j] *= (2.0e0);
    adjustment[j] -= log((double)arrmult[j]) / (2.0e0);
  }
  currentPrime = 3;
  while (currentPrime < 10000)
  {
    int halfCurrentPrime;

    NbrMod = (int)RemDivBigNbrByInt(Modulus, currentPrime, NumberLength);
    halfCurrentPrime = (currentPrime - 1) / 2;
    int jacobi = intDoubleModPow(NbrMod, halfCurrentPrime, currentPrime);
    double dp = (double)currentPrime;
    double logp = log(dp) / dp;
    for (j = 0; j<sizeof(arrmult)/sizeof(arrmult[0]); j++)
    {
      if (arrmult[j] == currentPrime)
      {
        adjustment[j] += logp;
      }
      else if (jacobi * intDoubleModPow(arrmult[j], halfCurrentPrime,
        currentPrime) % currentPrime == 1)
      {
        adjustment[j] += 2 * logp;
      }
    }
    do
    {
      currentPrime += 2;
      for (Q = 3; Q * Q <= currentPrime; Q += 2)
      { /* Check if currentPrime is prime */
        if (currentPrime % Q == 0)
        {
          break;  /* Composite */
        }
      }
    } while (Q * Q <= currentPrime);
  }  /* end while */
  multiplier = 1;
  for (j = 0; j<sizeof(arrmult)/sizeof(arrmult[0]); j++)
  {
    if (adjustment[j] > bestadjust)
    { /* find biggest adjustment */
      bestadjust = adjustment[j];
      multiplier = arrmult[j];
    }
  } /* end while */
  MultBigNbrByInt(TestNbr2, multiplier, Modulus, NumberLength);
  matrixPartialLength = nbrPrimes * 8;
  FactorBase = currentPrime;
  matrixBLength = nbrPrimes + 50;
  rowPrimeSieveData->modsqrt = (pNbrToFactor->x & 1) ? 1 : 0;
  switch ((int)Modulus[0] & 0x07)
  {
  case 1:
    logar2 = (unsigned char)3;
    break;
  case 5:
    logar2 = (unsigned char)1;
    break;
  default:
    logar2 = (unsigned char)1;
    break;
  }
  if (multiplier != 1 && multiplier != 2)
  {
    rowPrimeSieveData = &primeSieveData[2];
    rowPrimeTrialDivisionData = &primeTrialDivisionData[2];
    rowPrimeSieveData->value = multiplier;
    rowPrimeTrialDivisionData->value = multiplier;
    rowPrimeSieveData->modsqrt = 0;
    // The following works because multiplier has less than 16 significant bits.
    E = (int)((1U << BITS_PER_INT_GROUP) % multiplier);
    rowPrimeTrialDivisionData->exp1 = E;  // (2^31) mod multiplier
    D = E * E % multiplier;
    rowPrimeTrialDivisionData->exp2 = D;  // (2^31)^2 mod multiplier
    D = D * E % multiplier;
    rowPrimeTrialDivisionData->exp3 = D;  // (2^31)^3 mod multiplier
    D = D * E % multiplier;
    rowPrimeTrialDivisionData->exp4 = D;  // (2^31)^4 mod multiplier
    D = D * E % multiplier;
    rowPrimeTrialDivisionData->exp5 = D;  // (2^31)^5 mod multiplier
    D = D * E % multiplier;
    rowPrimeTrialDivisionData->exp6 = D;  // (2^31)^6 mod multiplier
    j = 3;
  }
  else
  {
    j = 2;
  }
  currentPrime = 3;
  while (j < nbrPrimes)
  { /* select small primes */
    NbrMod = (int)RemDivBigNbrByInt(Modulus, currentPrime,
      NumberLength);
    if (currentPrime != multiplier &&
      intDoubleModPow(NbrMod, (currentPrime - 1) / 2, currentPrime) == 1)
    {
      double dBase, dPower, dCurrentPrime, dRem;
      /* use only if Jacobi symbol = 0 or 1 */
      rowPrimeSieveData = &primeSieveData[j];
      rowPrimeTrialDivisionData = &primeTrialDivisionData[j];
      rowPrimeSieveData->value = (int)currentPrime;
      rowPrimeTrialDivisionData->value = (int)currentPrime;
      // The following works because multiplier has less than 26 significant bits.
      dBase = (double)((1U << BITS_PER_INT_GROUP) % currentPrime);
      rowPrimeTrialDivisionData->exp1 = (int)dBase;  // (2^31) mod currentPrime
      dCurrentPrime = (double)currentPrime;
      dPower = dBase * dBase;
      dPower -= floor(dPower / dCurrentPrime)*dCurrentPrime;
      rowPrimeTrialDivisionData->exp2 = (int)dPower; // (2^31)^2 mod currentPrime
      dPower *= dBase;
      dPower -= floor(dPower / dCurrentPrime)*dCurrentPrime;
      rowPrimeTrialDivisionData->exp3 = (int)dPower; // (2^31)^3 mod currentPrime
      dPower *= dBase;
      dPower -= floor(dPower / dCurrentPrime)*dCurrentPrime;
      rowPrimeTrialDivisionData->exp4 = (int)dPower; // (2^31)^4 mod currentPrime
      dPower *= dBase;
      dPower -= floor(dPower / dCurrentPrime)*dCurrentPrime;
      rowPrimeTrialDivisionData->exp5 = (int)dPower; // (2^31)^5 mod currentPrime
      dPower *= dBase;
      dPower -= floor(dPower / dCurrentPrime)*dCurrentPrime;
      rowPrimeTrialDivisionData->exp6 = (int)dPower; // (2^31)^6 mod currentPrime
      if ((currentPrime & 3) == 3)
      {
        SqrRootMod = intDoubleModPow(NbrMod, (currentPrime + 1) / 4, currentPrime);
      }
      else if ((currentPrime & 7) == 5)    // currentPrime = 5 (mod 8)
      {
        SqrRootMod =
          intDoubleModPow(NbrMod * 2, (currentPrime - 5) / 8, currentPrime);
        dRem = (double)2 * NbrMod * (double)SqrRootMod;
        dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
        dRem = dRem * (double)SqrRootMod - 1;
        dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
        dRem = dRem * (double)NbrMod;
        dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
        dRem = dRem * (double)SqrRootMod;
        dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
        SqrRootMod = (int)dRem;
      }
      else
      {
        Q = currentPrime - 1;
        E = 0;
        Power2 = 1;
        do
        {
          E++;
          Q /= 2;
          Power2 *= 2;
        } while ((Q & 1) == 0); /* E >= 3 */
        Power2 /= 2;
        X = 1;
        do
        {
          X++;
          Z = intDoubleModPow(X, Q, currentPrime);
        } while (intDoubleModPow(Z, Power2, currentPrime) == 1);
        Y = Z;
        X = intDoubleModPow(NbrMod, (Q - 1) / 2, currentPrime);
        dBase = (double)NbrMod * (double)X;
        V = (int)(dBase - floor(dBase / dCurrentPrime)*dCurrentPrime);
        dBase = (double)V * (double)X;
        W = (int)(dBase - floor(dBase / dCurrentPrime)*dCurrentPrime);
        while (W != 1)
        {
          T1 = 0;
          D = W;
          while (D != 1)
          {
            dBase = (double)D * (double)D;
            D = (int)(dBase - floor(dBase / dCurrentPrime)*dCurrentPrime);
            T1++;
          }
          D = intDoubleModPow(Y, 1 << (E - T1 - 1), currentPrime);
          dBase = (double)D * (double)D;
          Y1 = (int)(dBase - floor(dBase / dCurrentPrime)*dCurrentPrime);
          E = T1;
          dBase = (double)V * (double)D;
          V1 = (int)(dBase - floor(dBase / dCurrentPrime)*dCurrentPrime);
          dBase = (double)W * (double)Y1;
          W1 = (int)(dBase - floor(dBase / dCurrentPrime)*dCurrentPrime);
          Y = Y1;
          V = V1;
          W = W1;
        } /* end while */
        SqrRootMod = V;
      } /* end if */
      rowPrimeSieveData->modsqrt = (int)SqrRootMod;
      j++;
    } /* end while */
    do
    {
      currentPrime += 2;
      for (Q = 3; Q * Q <= currentPrime; Q += 2)
      { /* Check if currentPrime is prime */
        if (currentPrime % Q == 0)
        {
          break;  /* Composite */
        }
      }
    } while (Q * Q <= currentPrime);
  } /* End while */

  FactorBase = currentPrime;
  largePrimeUpperBound = 100 * FactorBase;
  // Convert array of limbs to BigInteger and find its logarithm.
  dlogNumberToFactor = logLimbs(pNbrToFactor, origNumberLength);
  dNumberToFactor = exp(dlogNumberToFactor);
#ifdef __EMSCRIPTEN__
  getMultAndFactorBase(multiplier, FactorBase);
  databack(lowerText);
#endif
  firstLimit = 2;
  for (j = 2; j < nbrPrimes; j++)
  {
    firstLimit *= (int)(primeSieveData[j].value);
    if (firstLimit > 2 * SieveLimit)
    {
      break;
    }
  }
  dNumberToFactor *= multiplier;
  smallPrimeUpperLimit = j + 1;
  threshold =
    (unsigned char)(log(
      sqrt(dNumberToFactor) * SieveLimit /
      (FactorBase * 64) /
      primeSieveData[j + 1].value)
      / log(3) + 0x81);
  firstLimit = (int)(log(dNumberToFactor) / 3);
  for (secondLimit = firstLimit; secondLimit < nbrPrimes; secondLimit++)
  {
    if (primeSieveData[secondLimit].value * 2 > SieveLimit)
    {
      break;
    }
  }
  for (thirdLimit = secondLimit; thirdLimit < nbrPrimes; thirdLimit++)
  {
    if (primeSieveData[thirdLimit].value > 2 * SieveLimit)
    {
      break;
    }
  }
  nbrPrimes2 = nbrPrimes - 4;
#ifdef __EMSCRIPTEN__
  newStamp = stamp();
  // Sieve start time in milliseconds.
  //  startTime = System.currentTimeMillis();
#endif
  Prod = sqrt(2 * dNumberToFactor) / (double)SieveLimit;
  fact = (int)pow(Prod, 1 / (float)nbrFactorsA);
  for (i = 2;; i++)
  {
    if (primeSieveData[i].value > fact)
    {
      break;
    }
  }
  span = nbrPrimes / (2 * nbrFactorsA*nbrFactorsA);
  if (nbrPrimes < 500)
  {
    span *= 2;
  }
  indexMinFactorA = i - span / 2;
  /*********************************************/
  /* Generate sieve threads                    */
  /*********************************************/
  sieveThread(&TempResult);
  NumberLength = origNumberLength;
  memcpy(pFactor, TempResult.limbs, TempResult.nbrLimbs * sizeof(limb));
  memset(pFactor + TempResult.nbrLimbs, 0, (NumberLength - TempResult.nbrLimbs) * sizeof(limb));

#if 0
  for (threadNumber = 0; threadNumber<numberThreads; threadNumber++)
  {
    //new Thread(this).start();                // Start new thread.
    //synchronized(amodq)
    {
//      while (threadArray[threadNumber] == null &&
//        getTerminateThread() == FALSE)
      {
        //try
        //{
        //  amodq.wait();
        //}
        //catch (InterruptedException ie) {}
      }
    }
  }
  //synchronized(matrixB)
  {
    while (factorSiqs == null && getTerminateThread() == FALSE)
    {
      try
      {
        matrixB.wait();
      }
      catch (InterruptedException ie) {}
    }
  }
#endif
  if (/*getTerminateThread() ||*/ (TempResult.nbrLimbs == 1 && TempResult.limbs[0].x == 0))
  {
    //throw new ArithmeticException();
  }
#if 0
  for (threadNumber = 0; threadNumber<numberThreads; threadNumber++)
  {                 // Wake up all sieve threads so they can terminate.
    if (threadArray[threadNumber].isAlive())
    {
      //try
      {
        threadArray[threadNumber].interrupt();
      }
      //catch (Exception e) {}
    }
  }
#endif
#if 0
  synchronized(this)
  {
    saveSIQSStatistics(polynomialsSieved, trialDivisions,
      smoothsFound, totalPartials,
      partialsFound, ValuesSieved);
  }
  return factorSiqs;
#endif
}

void ShowSIQSStatus(void)
{
#ifdef __EMSCRIPTEN__
  int elapsedTime = (int)(tenths() - originalTenthSecond);
  if (elapsedTime / 10 != oldTimeElapsed / 10)
  {
    oldTimeElapsed = elapsedTime;
    ShowSIQSInfo((elapsedTime - startSieveTenths)/10, congruencesFound, matrixBLength,
                 elapsedTime / 10);
  }
#endif
}

static int EraseSingletons(int nbrPrimes)
{
  int row, column, delta;
  int *rowMatrixB;
  int matrixBlength = matrixBLength;
  memset(newColumns, 0, matrixBlength*sizeof(int));
  // Find singletons in matrixB storing in array vectExpParity the number
  // of primes in each column.
  do
  {   // The singleton removal phase must run until there are no more
      // singletons to erase.
    for (column = nbrPrimes - 1; column >= 0; column--)
    {                  // Initialize number of primes per column to zero.
      vectExpParity[column] = 0;
    }
    for (row = matrixBlength - 1; row >= 0; row--)
    {                  // Traverse all rows of the matrix.
      rowMatrixB = matrixB[row];
      for (column = rowMatrixB[LENGTH_OFFSET]; column >= 1; column--)
      {                // A prime appeared in that column.
        vectExpParity[rowMatrixB[column]]++;
      }
    }
    row = 0;
    for (column = 0; column<nbrPrimes; column++)
    {
      if (vectExpParity[column] > 1)
      {                // Useful column found with at least 2 primes.
        newColumns[column] = row;
        primeTrialDivisionData[row++].value =
          primeTrialDivisionData[column].value;
      }
    }
    nbrPrimes = row;
    delta = 0;
    // Erase singletons from matrixB. The rows to be erased are those where the
    // the corresponding element of the array vectExpParity equals 1.
    for (row = 0; row < matrixBlength; row++)
    {                  // Traverse all rows of the matrix.
      rowMatrixB = matrixB[row];
      for (column = rowMatrixB[LENGTH_OFFSET]; column >= 1; column--)
      {                // Traverse all columns.
        if (vectExpParity[rowMatrixB[column]] == 1)
        {              // Singleton found: erase this row.
          delta++;
          break;
        }
      }
      if (column == 0)
      {                // Singleton not found: move row upwards.
        memcpy(matrixB[row - delta], matrixB[row], sizeof(matrixB[0]));
        memcpy(vectLeftHandSide[row - delta], vectLeftHandSide[row], sizeof(vectLeftHandSide[0]));
      }
    }
    matrixBlength -= delta;      // Update number of rows of the matrix.
    for (row = 0; row < matrixBlength; row++)
    {                  // Traverse all rows of the matrix.
      rowMatrixB = matrixB[row];
      for (column = rowMatrixB[LENGTH_OFFSET]; column >= 1; column--)
      {                // Change all column indexes in this row.
        rowMatrixB[column] = newColumns[rowMatrixB[column]];
      }
    }
  } while (delta > 0);           // End loop if number of rows did not
                                 // change.
  primeTrialDivisionData[0].exp2 = nbrPrimes;
  return matrixBlength;
}

/************************/
/* Linear algebra phase */
/************************/
static unsigned char LinearAlgebraPhase(
  int *biT, int *biR, int *biU,
  int NumberLength)
{
  int mask, row, col, j;
  int *rowMatrixB;
  int primeIndex;
  // Get new number of rows after erasing singletons.
  int matrixBlength = EraseSingletons(nbrPrimes);
  matrixBLength = matrixBlength;
#ifdef __EMSCRIPTEN__
  showMatrixSize((char *)SIQSInfoText, matrixBlength, primeTrialDivisionData[0].exp2);
#endif
  primeTrialDivisionData[0].exp2 = 0;         // Restore correct value.
  BlockLanczos();
  // The rows of matrixV indicate which rows must be multiplied so no
  // primes are multiplied an odd number of times.
  mask = 1;
  for (col = 31; col >= 0; col--)
  {
    int NumberLengthBak, index;

    IntToBigNbr(1, biT, NumberLength+1);
    IntToBigNbr(1, biR, NumberLength+1);
    memset(vectExpParity, 0, matrixBlength * sizeof(vectExpParity[0]));
    NumberLengthBak = NumberLength;
    if (Modulus[NumberLength - 1] == 0 && Modulus[NumberLength - 2] < HALF_INT_RANGE)
    {
      NumberLength--;
    }
    for (row = matrixBlength - 1; row >= 0; row--)
    {
      if ((matrixV[row] & mask) != 0)
      {
        MultBigNbrModN(vectLeftHandSide[row], biR, biU, Modulus,
          NumberLength);
        memcpy(biR, biU, (NumberLength + 1) * sizeof(biR[0]));
        rowMatrixB = matrixB[row];
        for (j = rowMatrixB[LENGTH_OFFSET]-1; j >= 1; j--)
        {
          primeIndex = rowMatrixB[j];
          vectExpParity[primeIndex] ^= 1;
          if (vectExpParity[primeIndex] == 0)
          {
            if (primeIndex == 0)
            {
              SubtractBigNbr(Modulus, biT, biT, NumberLength); // Multiply biT by -1.
            }
            else
            {
              MultBigNbrByIntModN(biT,
                primeTrialDivisionData[primeIndex].value, biT,
                Modulus, NumberLength);
            }
          }
        }
      }
    }
    NumberLength = NumberLengthBak;
    SubtractBigNbrModN(biR, biT, biR, Modulus, NumberLength);
    GcdBigNbr(biR, TestNbr2, biT, NumberLength);
    for (index = 1; index < NumberLength; index++)
    {
      if (biT[index] != 0)
      {
        break;
      }
    }
    if (index < NumberLength || biT[0] > 1)
    {   // GCD is not zero or 1.
      for (index = 0; index < NumberLength; index++)
      {
        if (biT[index] != TestNbr2[index])
        {
          break;
        }
      }
      if (index < NumberLength)
      { /* GCD is not 1 */
        return TRUE;
      }
    }
    mask *= 2;
  }
  return FALSE;
}

static unsigned char InsertNewRelation(
  int *rowMatrixB,
  int *biT, int *biU, int *biR,
  int NumberLengthMod)
{
  int i, k;
  int nbrColumns = rowMatrixB[LENGTH_OFFSET];
  // Insert it only if it is different from previous relations.
  if (congruencesFound >= matrixBLength)
  {                   // Discard excess congruences.
    return TRUE;
  }
  for (i = 0; i < congruencesFound; i++)
  {
    int *curRowMatrixB = matrixB[i];
    if (nbrColumns > curRowMatrixB[LENGTH_OFFSET])
    {
      continue;
    }
    if (nbrColumns == curRowMatrixB[LENGTH_OFFSET])
    {
      for (k = 1; k < nbrColumns; k++)
      {
        if (rowMatrixB[k] != curRowMatrixB[k])
        {
          break;
        }
      }
      if (k == nbrColumns)
      {
        return FALSE; // Do not insert same relation.
      }
      if (rowMatrixB[k] > curRowMatrixB[k])
      {
        continue;
      }
    }
    for (k = congruencesFound - 1; k >= i; k--)
    {
      memcpy(matrixB[k + 1], matrixB[k], matrixB[k][LENGTH_OFFSET]*sizeof(int));
      memcpy(vectLeftHandSide[k + 1], vectLeftHandSide[k], sizeof(vectLeftHandSide[0]));
    }
    break;
  }
  /* Convert negative numbers to the range 0 <= n < Modulus */
  if ((Modulus[0] & 1) == 0)
  {
    DivBigNbrByInt(Modulus, 2, TestNbr2, NumberLengthMod);
    // If biR >= Modulus perform biR = biR - Modulus.
    for (k = 0; k < NumberLengthMod; k++)
    {
      biT[k] = 0;
    }
    AddBigIntModN(biR, biT, biR, TestNbr2, NumberLengthMod);
    ModInvBigInt(biR, biT, TestNbr2, NumberLengthMod);
  }
  else
  {
    ModInvBigInt(biR, biT, Modulus, NumberLengthMod);
  }
  if ((biU[NumberLengthMod - 1] & HALF_INT_RANGE) != 0)
  {
    AddBigNbr(biU, Modulus, biU, NumberLengthMod);
  }

  // Compute biU / biR  (mod Modulus)
  MultBigNbrModN(biU, biT, biR, Modulus, NumberLengthMod);

  // Add relation to matrix B.
  memcpy(matrixB[i], &rowMatrixB[0], nbrColumns * sizeof(int));
  memcpy(vectLeftHandSide[i], biR, NumberLengthMod * sizeof(int));
  congruencesFound++;
  return TRUE;
}

static int intModInv(int NbrMod, int currentPrime)
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

/* Multiply binary matrices of length m x 32 by 32 x 32 */
/* The product matrix has size m x 32. Then add it to a m x 32 matrix. */
static void MatrixMultAdd(int *LeftMatr, int *RightMatr, int *ProdMatr)
{
  int matrLength = matrixBLength;
  int row;
  for (row = 0; row < matrLength; row++)
  {
    int col;
    int prodMatr = *(ProdMatr+row);
    int leftMatr = *(LeftMatr+row);
    for (col=0; col<32; col++)
    {
      if (leftMatr < 0)
      {
        prodMatr ^= *(RightMatr+col);
      }
      leftMatr *= 2;
    }
    *(ProdMatr+row) = prodMatr;
  }
}
/* Multiply binary matrices of length m x 32 by 32 x 32 */
/* The product matrix has size m x 32 */
static void MatrixMultiplication(int *LeftMatr, int *RightMatr, int *ProdMatr)
{
  int matrLength = 32;
  int row;
  for (row = 0; row < matrLength; row++)
  {
    int col;
    int prodMatr = 0;
    int leftMatr = *(LeftMatr+row);
    for (col=0; col < 32; col++)
    {
      if (leftMatr < 0)
      {
        prodMatr ^= *(RightMatr+col);
      }
      leftMatr *= 2;
    }
    *(ProdMatr+row) = prodMatr;
  }
}

/* Multiply the transpose of a binary matrix of length n x 32 by */
/* another binary matrix of length n x 32 */
/* The product matrix has size 32 x 32 */
static void MatrTranspMult(int matrLength, int *LeftMatr, int *RightMatr, int *ProdMatr)
{
  int row, col;
  int iMask = 1;
  for (col = 31; col >= 0; col--)
  {
    int prodMatr = 0;
    for (row = 0; row < matrLength; row++)
    {
      if ((*(LeftMatr+row) & iMask) != 0)
      {
        prodMatr ^= *(RightMatr+row);
      }
    }
    *(ProdMatr+col) = prodMatr;
    iMask *= 2;
  }
}

static void MatrixAddition(int *leftMatr, int *rightMatr, int *sumMatr)
{
  int row;
  for (row = 32 - 1; row >= 0; row--)
  {
    *(sumMatr+row) = *(leftMatr+row) ^ *(rightMatr+row);
  }
}

static void MatrMultBySSt(int length, int *Matr, int diagS, int *Prod)
{
  int row;
  for (row = length - 1; row >= 0; row--)
  {
    *(Prod+row) = diagS & *(Matr+row);
  }
}

/* Compute Bt * B * input matrix where B is the matrix that holds the */
/* factorization relations */
static void MultiplyAByMatrix(int *Matr, int *TempMatr, int *ProdMatr)
{
  int index;
  int row, rowValue;
  int *rowMatrixB;

  /* Compute TempMatr = B * Matr */
  memset(TempMatr, 0, matrixBLength*sizeof(int));
  for (row = matrixBLength - 1; row >= 0; row--)
  {
    rowMatrixB = matrixB[row];
    rowValue = *(Matr+row);
    for (index = *(rowMatrixB+LENGTH_OFFSET)-1; index >= 1; index--)
    {
      *(TempMatr+*(rowMatrixB+index)) ^= rowValue;
    }
  }

  /* Compute ProdMatr = Bt * TempMatr */
  for (row = matrixBLength - 1; row >= 0; row--)
  {
    int prodMatr = 0;
    rowMatrixB = matrixB[row];
    for (index = *(rowMatrixB+LENGTH_OFFSET)-1; index >= 1; index--)
    {
      prodMatr ^= *(TempMatr+*(rowMatrixB+index));
    }
    *(ProdMatr+row) = prodMatr;
  }
}

static void colexchange(int *XmY, int *V, int *V1, int *V2,
  int col1, int col2)
{
  int row;
  int mask1, mask2;
  int *matr1, *matr2;

  if (col1 == col2)
  {          // Cannot exchange the same column.
    return;
  }          // Exchange columns col1 and col2 of V1:V2
  mask1 = 0x80000000 >> (col1 & 31);
  mask2 = 0x80000000 >> (col2 & 31);
  matr1 = (col1 >= 32 ? V1 : V2);
  matr2 = (col2 >= 32 ? V1 : V2);
  for (row = matrixBLength - 1; row >= 0; row--)
  {             // If both bits are different toggle them.
    if (((matr1[row] & mask1) == 0) != ((matr2[row] & mask2) == 0))
    {           // If both bits are different toggle them.
      matr1[row] ^= mask1;
      matr2[row] ^= mask2;
    }
  }
  // Exchange columns col1 and col2 of XmY:V
  matr1 = (col1 >= 32 ? XmY : V);
  matr2 = (col2 >= 32 ? XmY : V);
  for (row = matrixBLength - 1; row >= 0; row--)
  {             // If both bits are different toggle them.
    if (((matr1[row] & mask1) == 0) != ((matr2[row] & mask2) == 0))
    {
      matr1[row] ^= mask1;
      matr2[row] ^= mask2;
    }
  }
}

static void coladd(int *XmY, int *V, int *V1, int *V2,
  int col1, int col2)
{
  int row;
  int mask1, mask2;
  int *matr1, *matr2;

  if (col1 == col2)
  {
    return;
  }               // Add column col1 to column col2 of V1:V2
  mask1 = 0x80000000 >> (col1 & 31);
  mask2 = 0x80000000 >> (col2 & 31);
  matr1 = (col1 >= 32 ? V1 : V2);
  matr2 = (col2 >= 32 ? V1 : V2);
  for (row = matrixBLength - 1; row >= 0; row--)
  {              // If bit to add is '1'...
    if ((matr1[row] & mask1) != 0)
    {            // Toggle bit in destination.
      matr2[row] ^= mask2;
    }
  }
  // Add column col1 to column col2 of XmY:V
  matr1 = (col1 >= 32 ? XmY : V);
  matr2 = (col2 >= 32 ? XmY : V);
  for (row = matrixBLength - 1; row >= 0; row--)
  {              // If bit to add is '1'...
    if ((matr1[row] & mask1) != 0)
    {            // Toggle bit in destination.
      matr2[row] ^= mask2;
    }
  }
}

static void BlockLanczos(void)
{
  int i, j, k;
  int oldDiagonalSSt, newDiagonalSSt;
  int index, mask;
  int matrixD[32];
  int matrixE[32];
  int matrixF[32];
  int matrixWinv[32];
  int matrixWinv1[32];
  int matrixWinv2[32];
  int matrixVtV0[32];
  int matrixVt1V0[32];
  int matrixVt2V0[32];
  int matrixVtAV[32];
  int matrixVt1AV1[32];
  int matrixCalcParenD[32];
  int vectorIndex[64];
  int matrixTemp[32];
  int matrixCalc1[32]; // Matrix that holds temporary data
  int matrixCalc2[32]; // Matrix that holds temporary data
  int *matr;
  double dSeed, dSeed2, dMult, dDivisor, dAdd;
  int Temp, Temp1;
  int stepNbr = 0;
  int currentOrder, currentMask;
  int row, col;
  int leftCol, rightCol;
  int minind, min, minanswer;
  int *rowMatrixB, *ptrMatrixV, *ptrMatrixXmY;

  newDiagonalSSt = oldDiagonalSSt = -1;
  memset(matrixWinv, 0, sizeof(matrixWinv));
  memset(matrixWinv1, 0, sizeof(matrixWinv1));
  memset(matrixWinv2, 0, sizeof(matrixWinv2));
  memset(matrixVtV0, 0, sizeof(matrixVtV0));
  memset(matrixVt1V0, 0, sizeof(matrixVt1V0));
  memset(matrixVt2V0, 0, sizeof(matrixVt2V0));
  memset(matrixVt1AV1, 0, sizeof(matrixVt1AV1));

  /* Initialize matrix X-Y and matrix V_0 with random data */
  dSeed = (double)123456789;
  dMult = (double)62089911;
  dAdd = (double)54325442;
  dDivisor = (double)0x7FFFFFFF;
  ptrMatrixXmY = &matrixXmY[matrixBLength - 1];
  for (ptrMatrixV = &matrixV[matrixBLength - 1]; ptrMatrixV >= matrixV; ptrMatrixV--)
  {
    dSeed2 = (dSeed * dMult + dAdd);
    dSeed2 -= floor(dSeed2 / dDivisor) * dDivisor;
    *ptrMatrixXmY-- = (int)dSeed + (int)dSeed2;
    dSeed = (dSeed2 * dMult + dAdd);
    dSeed -= floor(dSeed / dDivisor) * dDivisor;
    dSeed2 = (dSeed * dMult + dAdd);
    dSeed2 -= floor(dSeed2 / dDivisor) * dDivisor;
    *ptrMatrixV = (int)dSeed + (int)dSeed2;
    dSeed = (dSeed2 * dMult + dAdd);
    dSeed -= floor(dSeed / dDivisor) * dDivisor;
  }
  // Compute matrix Vt(0) * V(0)
  MatrTranspMult(matrixBLength, matrixV, matrixV, matrixVtV0);
  for (;;)
  {
    int indexC;
    //if (getTerminateThread())
    //{
    //  throw new ArithmeticException();
    // }
    oldDiagonalSSt = newDiagonalSSt;
    stepNbr++;
    // Compute matrix A * V(i)
    MultiplyAByMatrix(matrixV, matrixCalc3, matrixAV);
    // Compute matrix Vt(i) * A * V(i)
    MatrTranspMult(matrixBLength, matrixV, matrixAV, matrixVtAV);

    /* If Vt(i) * A * V(i) = 0, end of loop */
    for (i = sizeof(matrixVtAV)/sizeof(matrixVtAV[0]) - 1; i >= 0; i--)
    {
      if (matrixVtAV[i] != 0)
      {
        break;
      }
    }
    if (i < 0)
    {
      break;
    } /* End X-Y calculation loop */

      /* Selection of S(i) and W(i) */

    memcpy(matrixTemp, matrixWinv2, sizeof(matrixTemp));
    memcpy(matrixWinv2, matrixWinv1, sizeof(matrixTemp));
    memcpy(matrixWinv1, matrixWinv, sizeof(matrixTemp));
    memcpy(matrixWinv, matrixTemp, sizeof(matrixTemp));

    mask = 1;
    for (j = 31; j >= 0; j--)
    {
      matrixD[j] = matrixVtAV[j]; /*  D = VtAV    */
      matrixWinv[j] = mask; /*  Winv = I    */
      mask *= 2;
    }

    index = 31;
    mask = 1;
    for (indexC = 31; indexC >= 0; indexC--)
    {
      if ((oldDiagonalSSt & mask) != 0)
      {
        matrixE[index] = indexC;
        matrixF[index] = mask;
        index--;
      }
      mask *= 2;
    }
    mask = 1;
    for (indexC = 31; indexC >= 0; indexC--)
    {
      if ((oldDiagonalSSt & mask) == 0)
      {
        matrixE[index] = indexC;
        matrixF[index] = mask;
        index--;
      }
      mask *= 2;
    }
    newDiagonalSSt = 0;
    for (j = 0; j < 32; j++)
    {
      currentOrder = matrixE[j];
      currentMask = matrixF[j];
      for (k = j; k < 32; k++)
      {
        if ((matrixD[matrixE[k]] & currentMask) != 0)
        {
          break;
        }
      }
      if (k < 32)
      {
        i = matrixE[k];
        Temp = matrixWinv[i];
        matrixWinv[i] = matrixWinv[currentOrder];
        matrixWinv[currentOrder] = Temp;
        Temp1 = matrixD[i];
        matrixD[i] = matrixD[currentOrder];
        matrixD[currentOrder] = Temp1;
        newDiagonalSSt |= currentMask;
        for (k = 31; k >= 0; k--)
        {
          if (k != currentOrder && ((matrixD[k] & currentMask) != 0))
          {
            matrixWinv[k] ^= Temp;
            matrixD[k] ^= Temp1;
          }
        } /* end for k */
      }
      else
      {
        for (k = j; k < 31; k++)
        {
          if ((matrixWinv[matrixE[k]] & currentMask) != 0)
          {
            break;
          }
        }
        i = matrixE[k];
        Temp = matrixWinv[i];
        matrixWinv[i] = matrixWinv[currentOrder];
        matrixWinv[currentOrder] = Temp;
        Temp1 = matrixD[i];
        matrixD[i] = matrixD[currentOrder];
        matrixD[currentOrder] = Temp1;
        for (k = 31; k >= 0; k--)
        {
          if ((matrixWinv[k] & currentMask) != 0)
          {
            matrixWinv[k] ^= Temp;
            matrixD[k] ^= Temp1;
          }
        } /* end for k */
      } /* end if */
    } /* end for j */
      /* Compute D(i), E(i) and F(i) */
    if (stepNbr >= 3)
    {
      // F = -Winv(i-2) * (I - Vt(i-1)*A*V(i-1)*Winv(i-1)) * ParenD * S*St
      MatrixMultiplication(matrixVt1AV1, matrixWinv1, matrixCalc2);
      mask = 1; /* Add identity matrix */
      for (index = 31; index >= 0; index--)
      {
        matrixCalc2[index] ^= mask;
        mask*=2;
      }
      MatrixMultiplication(matrixWinv2, matrixCalc2, matrixCalc1);
      MatrixMultiplication(matrixCalc1, matrixCalcParenD, matrixF);
      MatrMultBySSt(32, matrixF, newDiagonalSSt, matrixF);
    }
    // E = -Winv(i-1) * Vt(i)*A*V(i) * S*St
    if (stepNbr >= 2)
    {
      MatrixMultiplication(matrixWinv1, matrixVtAV, matrixE);
      MatrMultBySSt(32, matrixE, newDiagonalSSt, matrixE);
    }
    // ParenD = Vt(i)*A*A*V(i) * S*St + Vt(i)*A*V(i)
    // D = I - Winv(i) * ParenD
    MatrTranspMult(matrixBLength, matrixAV, matrixAV, matrixCalc1); // Vt(i)*A*A*V(i)
    MatrMultBySSt(32, matrixCalc1, newDiagonalSSt, matrixCalc1);
    MatrixAddition(matrixCalc1, matrixVtAV, matrixCalcParenD);
    MatrixMultiplication(matrixWinv, matrixCalcParenD, matrixD);
    mask = 1; /* Add identity matrix */
    for (index = 31; index >= 0; index--)
    {
      matrixD[index] ^= mask;
      mask*=2;
    }

    /* Update value of X - Y */
    MatrixMultiplication(matrixWinv, matrixVtV0, matrixCalc1);
    MatrixMultAdd(matrixV, matrixCalc1, matrixXmY);

    /* Compute value of new matrix V(i) */
    // V(i+1) = A * V(i) * S * St + V(i) * D + V(i-1) * E + V(i-2) * F
    MatrMultBySSt(matrixBLength, matrixAV, newDiagonalSSt, matrixCalc3);
    MatrixMultAdd(matrixV, matrixD, matrixCalc3);
    if (stepNbr >= 2)
    {
      MatrixMultAdd(matrixV1, matrixE, matrixCalc3);
      if (stepNbr >= 3)
      {
        MatrixMultAdd(matrixV2, matrixF, matrixCalc3);
      }
    }
    /* Compute value of new matrix Vt(i)V0 */
    // Vt(i+1)V(0) = Dt * Vt(i)V(0) + Et * Vt(i-1)V(0) + Ft * Vt(i-2)V(0)
    MatrTranspMult(32, matrixD, matrixVtV0, matrixCalc2);
    if (stepNbr >= 2)
    {
      MatrTranspMult(32, matrixE, matrixVt1V0, matrixCalc1);
      MatrixAddition(matrixCalc1, matrixCalc2, matrixCalc2);
      if (stepNbr >= 3)
      {
        MatrTranspMult(32, matrixF, matrixVt2V0, matrixCalc1);
        MatrixAddition(matrixCalc1, matrixCalc2, matrixCalc2);
      }
    }
    memcpy(matrixTemp2, matrixV2, sizeof(matrixTemp2));
    memcpy(matrixV2, matrixV1, sizeof(matrixV2));
    memcpy(matrixV1, matrixV, sizeof(matrixV1));
    memcpy(matrixV, matrixCalc3, sizeof(matrixV));
    memcpy(matrixCalc3, matrixTemp2, sizeof(matrixCalc3));
    memcpy(matrixTemp, matrixVt2V0, sizeof(matrixTemp));
    memcpy(matrixVt2V0, matrixVt1V0, sizeof(matrixVt2V0));
    memcpy(matrixVt1V0, matrixVtV0, sizeof(matrixVt1V0));
    memcpy(matrixVtV0, matrixCalc2, sizeof(matrixVtV0));
    memcpy(matrixCalc2, matrixTemp, sizeof(matrixCalc2));
    memcpy(matrixTemp, matrixVt1AV1, sizeof(matrixTemp));
    memcpy(matrixVt1AV1, matrixVtAV, sizeof(matrixVt1AV1));
    memcpy(matrixVtAV, matrixTemp, sizeof(matrixVtAV));
  } /* end while */

    /* Find matrix V1:V2 = B * (X-Y:V) */
  for (row = matrixBLength - 1; row >= 0; row--)
  {
    matrixV1[row] = matrixV2[row] = 0;
  }
  for (row = matrixBLength - 1; row >= 0; row--)
  {
    int rowMatrixXmY, rowMatrixV;

    rowMatrixB = matrixB[row];
    rowMatrixXmY = matrixXmY[row];
    rowMatrixV = matrixV[row];
    // The vector rowMatrixB includes the indexes of the columns set to '1'.
    for (index = rowMatrixB[LENGTH_OFFSET]-1; index >= 1; index--)
    {
      col = rowMatrixB[index];
      matrixV1[col] ^= rowMatrixXmY;
      matrixV2[col] ^= rowMatrixV;
    }
  }
  rightCol = 64;
  leftCol = 0;
  while (leftCol < rightCol)
  {
    for (col = leftCol; col < rightCol; col++)
    {       // For each column find the first row which has a '1'.
            // Columns outside this range must have '0' in all rows.
      matr = (col >= 32 ? matrixV1 : matrixV2);
      mask = 0x80000000 >> (col & 31);
      vectorIndex[col] = -1;    // indicate all rows in zero in advance.
      for (row = 0; row < matrixBLength; row++)
      {
        if ((matr[row] & mask) != 0)
        {               // First row for this mask is found. Store it.
          vectorIndex[col] = row;
          break;
        }
      }
    }
    for (col = leftCol; col < rightCol; col++)
    {
      if (vectorIndex[col] < 0)
      {  // If all zeros in col 'col', exchange it with first column with
         // data different from zero (leftCol).
        colexchange(matrixXmY, matrixV, matrixV1, matrixV2, leftCol, col);
        vectorIndex[col] = vectorIndex[leftCol];
        vectorIndex[leftCol] = -1;  // This column now has zeros.
        leftCol++;                  // Update leftCol to exclude that column.
      }
    }
    if (leftCol == rightCol)
    {
      break;
    }
    // At this moment all columns from leftCol to rightCol are non-zero.
    // Get the first row that includes a '1'.
    min = vectorIndex[leftCol];
    minind = leftCol;
    for (col = leftCol + 1; col < rightCol; col++)
    {
      if (vectorIndex[col] < min)
      {
        min = vectorIndex[col];
        minind = col;
      }
    }
    minanswer = 0;
    for (col = leftCol; col < rightCol; col++)
    {
      if (vectorIndex[col] == min)
      {
        minanswer++;
      }
    }
    if (minanswer > 1)
    {            // Two columns with the same first row to '1'.
      for (col = minind + 1; col < rightCol; col++)
      {
        if (vectorIndex[col] == min)
        {        // Add first column which has '1' in the same row to
                 // the other columns so they have '0' in this row after
                 // this operation.
          coladd(matrixXmY, matrixV, matrixV1, matrixV2, minind, col);
        }
      }
    }
    else
    {
      rightCol--;
      colexchange(matrixXmY, matrixV, matrixV1, matrixV2, minind, rightCol);
    }
  }
  leftCol = 0; /* find linear independent solutions */
  while (leftCol < rightCol)
  {
    for (col = leftCol; col < rightCol; col++)
    {         // For each column find the first row which has a '1'.
      matr = (col >= 32 ? matrixXmY : matrixV);
      mask = 0x80000000 >> (col & 31);
      vectorIndex[col] = -1;    // indicate all rows in zero in advance.
      for (row = 0; row < matrixBLength; row++)
      {
        if ((matr[row] & mask) != 0)
        {         // First row for this mask is found. Store it.
          vectorIndex[col] = row;
          break;
        }
      }
    }
    for (col = leftCol; col < rightCol; col++)
    {  // If all zeros in col 'col', exchange it with last column with
       // data different from zero (rightCol).
      if (vectorIndex[col] < 0)
      {
        rightCol--;                 // Update rightCol to exclude that column.
        colexchange(matrixXmY, matrixV, matrixV1, matrixV2, rightCol, col);
        vectorIndex[col] = vectorIndex[rightCol];
        vectorIndex[rightCol] = -1; // This column now has zeros.
      }
    }
    if (leftCol == rightCol)
    {
      break;
    }
    // At this moment all columns from leftCol to rightCol are non-zero.
    // Get the first row that includes a '1'.
    min = vectorIndex[leftCol];
    minind = leftCol;
    for (col = leftCol + 1; col < rightCol; col++)
    {
      if (vectorIndex[col] < min)
      {
        min = vectorIndex[col];
        minind = col;
      }
    }
    minanswer = 0;
    for (col = leftCol; col < rightCol; col++)
    {
      if (vectorIndex[col] == min)
      {
        minanswer++;
      }
    }
    if (minanswer > 1)
    {            // At least two columns with the same first row to '1'.
      for (col = minind + 1; col < rightCol; col++)
      {
        if (vectorIndex[col] == min)
        {        // Add first column which has '1' in the same row to
                 // the other columns so they have '0' in this row after
                 // this operation.
          coladd(matrixXmY, matrixV, matrixV1, matrixV2, minind, col);
        }
      }
    }
    else
    {
      colexchange(matrixXmY, matrixV, matrixV1, matrixV2, minind, leftCol);
      leftCol++;
    }
  }
}

/****************/
/* Sieve thread */
/****************/
static void sieveThread(BigInteger *result)
{
  int polySet;
  int biT[MAX_LIMBS_SIQS];
  int biU[MAX_LIMBS_SIQS];
  int biV[MAX_LIMBS_SIQS];
  int biR[MAX_LIMBS_SIQS];
//  int multiplier = this.multiplier;
  PrimeSieveData *rowPrimeSieveData;
  PrimeSieveData *rowPrimeSieveData0;
  PrimeTrialDivisionData *rowPrimeTrialDivisionData;
  short SieveArray[100000];
  int rowPartials[200];
  int biLinearCoeff[MAX_LIMBS_SIQS];
//  int threadNumber = this.threadNumber;
  int biDividend[MAX_LIMBS_SIQS];
  int biAbsLinearCoeff[MAX_LIMBS_SIQS];
  int indexFactorsA[50];
  int rowSquares[200];
  int polynomialsPerThread = ((NbrPolynomials - 1) / numberThreads) & 0xFFFFFFFE;
  int firstPolynomial = threadNumber*polynomialsPerThread;
  int lastPolynomial = firstPolynomial + polynomialsPerThread;
  firstPolynomial |= 1;
  int grayCode = firstPolynomial ^ (firstPolynomial >> 1);
  firstPolynomial++;
  int i, PolynomialIndex, index, index2;
  int currentPrime;
  int RemB, D, Q;
  int Dividend[MAX_LIMBS_SIQS];
  int rowMatrixBbeforeMerge[200];
  int rowMatrixB[200];
  unsigned char positive;
  int inverseA, twiceInverseA;
  int NumberLengthA, NumberLengthB;

  memset(biLinearCoeff, 0, sizeof(biLinearCoeff));
//  synchronized(amodq)
  {
    if (threadNumber == 0)
    {
//      primeSieveData = this.primeSieveData;
    }
    else
    {
      for (i = 0; i<nbrPrimes; i++)
      {
        rowPrimeSieveData = &primeSieveData[i];
        rowPrimeSieveData0 = &primeSieveData[i];
        rowPrimeSieveData->value = rowPrimeSieveData0->value;
        rowPrimeSieveData->modsqrt = rowPrimeSieveData0->modsqrt;
        memcpy(rowPrimeSieveData->Bainv2, rowPrimeSieveData0->Bainv2, sizeof(rowPrimeSieveData0->Bainv2));
      }
    }
//    threadArray[threadNumber] = Thread.currentThread();
//    amodq.notifyAll();
  }               // End synchronized block.
    for (polySet = 1;; polySet++)
    {                         // For each polynomial set...
      //if (getTerminateThread())
      //{
      //  throw new ArithmeticException();
      //}
      //synchronized(amodq)
      {
        nbrThreadFinishedPolySet++;
        if (congruencesFound >= matrixBLength/* || factorSiqs != null*/)
        {
          if (nbrThreadFinishedPolySet < polySet * numberThreads)
          {
            return;
          }
          if (1/*factorSiqs == null*/)
          {
            while (!LinearAlgebraPhase(biT, biR, biU, NumberLength));
            BigIntToBigNbr(result, biT, NumberLength);  // Factor found.
#if 0
            synchronized(matrixB)
            {
              factorSiqs = result;
              matrixB.notify();
            }
#endif
          }
          else
          {
#if 0
            synchronized(matrixB)
            {
              matrixB.notify();
            }
#endif
          }
          return;
        }
        if (nbrThreadFinishedPolySet == polySet * numberThreads)
        {
          /*********************************************/
          /* Initialization stage for first polynomial */
          /*********************************************/
          firstPrimeSieveData = primeSieveData;
          oldSeed = newSeed;
          newSeed = getFactorsOfA(oldSeed, aindex);
          for (index = 0; index<nbrFactorsA; index++)
          {                        // Get the values of the factors of A.
            afact[index] = primeSieveData[aindex[index]].value;
          }
          // Compute the leading coefficient in biQuadrCoeff.

          IntToBigNbr(afact[0], biQuadrCoeff, NumberLength);
          for (index = 1; index < nbrFactorsA; index++)
          {
            MultBigNbrByInt(biQuadrCoeff, afact[index], biQuadrCoeff,
              NumberLength);
          }
          for (NumberLengthA = NumberLength; NumberLengthA >= 2; NumberLengthA--)
          {
            if (biQuadrCoeff[NumberLengthA - 1] != 0 ||
              biQuadrCoeff[NumberLengthA - 2] >= HALF_INT_RANGE)
            {
              break;
            }
          }
          for (index = 0; index < nbrFactorsA; index++)
          {
            currentPrime = afact[index];
            D = RemDivBigNbrByInt(biQuadrCoeff,
              currentPrime*currentPrime, NumberLengthA) / currentPrime;
            Q = primeSieveData[aindex[index]].modsqrt *
              intModInv(D, currentPrime) % currentPrime;
            amodq[index] = D << 1;
            tmodqq[index] = RemDivBigNbrByInt(Modulus,
              currentPrime*currentPrime, NumberLength);
            if (Q + Q > currentPrime)
            {
              Q = currentPrime - Q;
            }
            DivBigNbrByInt(biQuadrCoeff, currentPrime, biDividend,
              NumberLengthA);
            MultBigNbrByInt(biDividend, Q, biLinearDelta[index],
              NumberLengthA);
            for (index2 = NumberLengthA; index2 < NumberLength; index2++)
            {
              biLinearDelta[index][index2] = 0;
            }
          }
          for (index = 1; index < nbrPrimes; index++)
          {
            double dRem, dCurrentPrime;
            rowPrimeTrialDivisionData = &primeTrialDivisionData[index];
            rowPrimeSieveData = &primeSieveData[index];
            // Get current prime.
            currentPrime = rowPrimeTrialDivisionData->value;
            memcpy(Dividend, biQuadrCoeff, sizeof(Dividend));   // Get A mod current prime.
            switch (NumberLengthA)
            {
            case 7:
              dRem = (double)Dividend[6] * (double)rowPrimeTrialDivisionData->exp6 +
                (double)Dividend[5] * (double)rowPrimeTrialDivisionData->exp5 +
                (double)Dividend[4] * (double)rowPrimeTrialDivisionData->exp4 +
                (double)Dividend[3] * (double)rowPrimeTrialDivisionData->exp3 +
                (double)Dividend[2] * (double)rowPrimeTrialDivisionData->exp2 +
                (double)Dividend[1] * (double)rowPrimeTrialDivisionData->exp1;
              break;
            case 6:
              dRem = (double)Dividend[5] * (double)rowPrimeTrialDivisionData->exp5 +
                (double)Dividend[4] * (double)rowPrimeTrialDivisionData->exp4 +
                (double)Dividend[3] * (double)rowPrimeTrialDivisionData->exp3 +
                (double)Dividend[2] * (double)rowPrimeTrialDivisionData->exp2 +
                (double)Dividend[1] * (double)rowPrimeTrialDivisionData->exp1;
              break;
            case 5:
              dRem = (double)Dividend[4] * (double)rowPrimeTrialDivisionData->exp4 +
                (double)Dividend[3] * (double)rowPrimeTrialDivisionData->exp3 +
                (double)Dividend[2] * (double)rowPrimeTrialDivisionData->exp2 +
                (double)Dividend[1] * (double)rowPrimeTrialDivisionData->exp1;
              break;
            case 4:
              dRem = (double)Dividend[3] * (double)rowPrimeTrialDivisionData->exp3 +
                (double)Dividend[2] * (double)rowPrimeTrialDivisionData->exp2 +
                (double)Dividend[1] * (double)rowPrimeTrialDivisionData->exp1;
              break;
            default:
              dRem = (double)Dividend[2] * (double)rowPrimeTrialDivisionData->exp2 +
                (double)Dividend[1] * (double)rowPrimeTrialDivisionData->exp1;
              break;
            }
            dRem += (double)Dividend[0];
            dCurrentPrime = (double)currentPrime;
            dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
            // Get its inverse
            inverseA = intModInv((int)dRem, currentPrime);
            twiceInverseA = inverseA << 1;       // and twice this value.
            dRem = (double)twiceInverseA * (double)rowPrimeSieveData->modsqrt;
            dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
            rowPrimeSieveData->difsoln = (int)dRem;
            switch (NumberLengthA)
            {
            case 7:
              for (index2 = nbrFactorsA - 1; index2 > 0; index2--)
              {
                memcpy(Dividend, biLinearDelta[index2], sizeof(Dividend));
                dRem = (double)Dividend[6] * (double)rowPrimeTrialDivisionData->exp6 +
                  (double)Dividend[5] * (double)rowPrimeTrialDivisionData->exp5 +
                  (double)Dividend[4] * (double)rowPrimeTrialDivisionData->exp4 +
                  (double)Dividend[3] * (double)rowPrimeTrialDivisionData->exp3 +
                  (double)Dividend[2] * (double)rowPrimeTrialDivisionData->exp2 +
                  (double)Dividend[1] * (double)rowPrimeTrialDivisionData->exp1 +
                  (double)Dividend[0];
                dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
                dRem *= twiceInverseA;
                dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
                rowPrimeSieveData->Bainv2[index2 - 1] = (int)dRem;
              }
              memcpy(Dividend, biLinearDelta[0], sizeof(Dividend));
              dRem = (double)Dividend[6] * (double)rowPrimeTrialDivisionData->exp6 +
                (double)Dividend[5] * (double)rowPrimeTrialDivisionData->exp5 +
                (double)Dividend[4] * (double)rowPrimeTrialDivisionData->exp4 +
                (double)Dividend[3] * (double)rowPrimeTrialDivisionData->exp3 +
                (double)Dividend[2] * (double)rowPrimeTrialDivisionData->exp2 +
                (double)Dividend[1] * (double)rowPrimeTrialDivisionData->exp1 +
                (double)Dividend[0];
              dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
              break;
            case 6:
              for (index2 = nbrFactorsA - 1; index2 > 0; index2--)
              {
                memcpy(Dividend, biLinearDelta[index2], sizeof(Dividend));
                dRem = (double)Dividend[5] * (double)rowPrimeTrialDivisionData->exp5 +
                  (double)Dividend[4] * (double)rowPrimeTrialDivisionData->exp4 +
                  (double)Dividend[3] * (double)rowPrimeTrialDivisionData->exp3 +
                  (double)Dividend[2] * (double)rowPrimeTrialDivisionData->exp2 +
                  (double)Dividend[1] * (double)rowPrimeTrialDivisionData->exp1 +
                  (double)Dividend[0];
                dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
                dRem *= twiceInverseA;
                dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
                rowPrimeSieveData->Bainv2[index2 - 1] = (int)dRem;
              }
              memcpy(Dividend, biLinearDelta[0], sizeof(Dividend));
              dRem = (double)Dividend[5] * (double)rowPrimeTrialDivisionData->exp5 +
                (double)Dividend[4] * (double)rowPrimeTrialDivisionData->exp4 +
                (double)Dividend[3] * (double)rowPrimeTrialDivisionData->exp3 +
                (double)Dividend[2] * (double)rowPrimeTrialDivisionData->exp2 +
                (double)Dividend[1] * (double)rowPrimeTrialDivisionData->exp1 +
                (double)Dividend[0];
              dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
              break;
            case 5:
              for (index2 = nbrFactorsA - 1; index2 > 0; index2--)
              {
                memcpy(Dividend, biLinearDelta[index2], sizeof(Dividend));
                dRem = (double)Dividend[4] * (double)rowPrimeTrialDivisionData->exp4 +
                  (double)Dividend[3] * (double)rowPrimeTrialDivisionData->exp3 +
                  (double)Dividend[2] * (double)rowPrimeTrialDivisionData->exp2 +
                  (double)Dividend[1] * (double)rowPrimeTrialDivisionData->exp1 +
                  (double)Dividend[0];
                dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
                dRem *= twiceInverseA;
                dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
                rowPrimeSieveData->Bainv2[index2 - 1] = (int)dRem;
              }
              memcpy(Dividend, biLinearDelta[0], sizeof(Dividend));
              dRem = (double)Dividend[4] * (double)rowPrimeTrialDivisionData->exp4 +
                (double)Dividend[3] * (double)rowPrimeTrialDivisionData->exp3 +
                (double)Dividend[2] * (double)rowPrimeTrialDivisionData->exp2 +
                (double)Dividend[1] * (double)rowPrimeTrialDivisionData->exp1 +
                (double)Dividend[0];
              dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
              break;
            case 4:
              for (index2 = nbrFactorsA - 1; index2 > 0; index2--)
              {
                memcpy(Dividend, biLinearDelta[index2], sizeof(Dividend));
                dRem = (double)Dividend[3] * (double)rowPrimeTrialDivisionData->exp3 +
                  (double)Dividend[2] * (double)rowPrimeTrialDivisionData->exp2 +
                  (double)Dividend[1] * (double)rowPrimeTrialDivisionData->exp1 +
                  (double)Dividend[0];
                dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
                dRem *= twiceInverseA;
                dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
                rowPrimeSieveData->Bainv2[index2 - 1] = (int)dRem;
              }
              memcpy(Dividend, biLinearDelta[0], sizeof(Dividend));
              dRem = (double)Dividend[3] * (double)rowPrimeTrialDivisionData->exp3 +
                (double)Dividend[2] * (double)rowPrimeTrialDivisionData->exp2 +
                (double)Dividend[1] * (double)rowPrimeTrialDivisionData->exp1 +
                (double)Dividend[0];
              dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
              break;
            default:
              for (index2 = nbrFactorsA - 1; index2 > 0; index2--)
              {
                memcpy(Dividend, biLinearDelta[index2], sizeof(Dividend));
                dRem = (double)Dividend[2] * (double)rowPrimeTrialDivisionData->exp2 +
                  (double)Dividend[1] * (double)rowPrimeTrialDivisionData->exp1 +
                  (double)Dividend[0];
                dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
                dRem *= twiceInverseA;
                dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
                rowPrimeSieveData->Bainv2[index2 - 1] = (int)dRem;
              }
              memcpy(Dividend, biLinearDelta[0], sizeof(Dividend));
              dRem = (double)Dividend[2] * (double)rowPrimeTrialDivisionData->exp2 +
                (double)Dividend[1] * (double)rowPrimeTrialDivisionData->exp1 +
                (double)Dividend[0];
              dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
              break;
            }
            dRem *= twiceInverseA;
            dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
            rowPrimeSieveData->Bainv2_0 = (int)dRem;
            if (rowPrimeSieveData->Bainv2_0 != 0)
            {
              rowPrimeSieveData->Bainv2_0 =
                currentPrime - rowPrimeSieveData->Bainv2_0;
            }
          }
          for (index2 = 0; index2 < nbrFactorsA; index2++)
          {
            primeSieveData[aindex[index2]].difsoln = -1; // Do not sieve.
          }
          //synchronized(TestNbr2)
          {
            //TestNbr2.notifyAll();
          }
        }           // End initializing first polynomial
      }             // End synchronized
#if 0
      synchronized(TestNbr2)
      {
        while (nbrThreadFinishedPolySet < polySet * numberThreads)
        {
          try
          {
            TestNbr2.wait();
          }
          catch (InterruptedException ie) {}
        }
      }
#endif
      if (/*factorSiqs != null ||*/ congruencesFound >= matrixBLength)
      {
        if (nbrThreadFinishedPolySet > numberThreads*polySet)
        {
          continue;
        }
        //synchronized(amodq)
        {
          nbrThreadFinishedPolySet++;
        }
        return;
      }
      PolynomialIndex = firstPolynomial;
      // Compute first polynomial parameters.
      for (i = 0; i<NumberLength; i++)
      {
        biLinearCoeff[i] = biLinearDelta[0][i];
      }
      for (i = 1; i<nbrFactorsA; i++)
      {
        if ((grayCode & (1 << i)) == 0)
        {
          AddBigNbr(biLinearCoeff, biLinearDelta[i], biLinearCoeff,
            NumberLength);
        }
        else
        {
          SubtractBigNbr(biLinearCoeff, biLinearDelta[i], biLinearCoeff,
            NumberLength);
        }
      }
      for (NumberLengthA = NumberLength; NumberLengthA >= 2; NumberLengthA--)
      {
        if (biQuadrCoeff[NumberLengthA - 1] != 0 ||
          biQuadrCoeff[NumberLengthA - 2] >= HALF_INT_RANGE)
        {                             // Go out if significant limb.
          break;
        }
      }
      if (biLinearCoeff[NumberLength - 1] >= HALF_INT_RANGE)
      {                               // Number is negative.
        positive = FALSE;
        memcpy(biT, biLinearCoeff, NumberLength * sizeof(biT[0]));
        ChSignBigNbr(biT, NumberLength);   // Make it positive.
        memcpy(biAbsLinearCoeff, biT, sizeof(biT));
      }
      else
      {
        positive = TRUE;                       // B is positive.
                                               // Get B mod current prime. 
        memcpy(biAbsLinearCoeff, biLinearCoeff, sizeof(biLinearCoeff));
      }
      for (NumberLengthB = NumberLength; NumberLengthB >= 2; NumberLengthB--)
      {
        if (biAbsLinearCoeff[NumberLengthB - 1] != 0 ||
          biAbsLinearCoeff[NumberLengthB - 2] >= HALF_INT_RANGE)
        {                                // Go out if significant limb.
          break;
        }
      }
      for (i = nbrPrimes - 1; i>0; i--)
      {
        double dRem, dCurrentPrime;
        rowPrimeSieveData = &primeSieveData[i];
        rowPrimeSieveData0 = firstPrimeSieveData+i;
        rowPrimeSieveData->difsoln = rowPrimeSieveData0->difsoln;
        rowPrimeSieveData->Bainv2_0 = rowPrimeSieveData0->Bainv2_0;
        rowPrimeTrialDivisionData = &primeTrialDivisionData[i];
        currentPrime = rowPrimeTrialDivisionData->value;     // Get current prime.
        memcpy(Dividend, biQuadrCoeff,sizeof(biQuadrCoeff)); // Get A mod current prime.
        switch (NumberLengthA)
        {
        case 7:
          dRem = (double)Dividend[6] * (double)rowPrimeTrialDivisionData->exp6 +
            (double)Dividend[5] * (double)rowPrimeTrialDivisionData->exp5 +
            (double)Dividend[4] * (double)rowPrimeTrialDivisionData->exp4 +
            (double)Dividend[3] * (double)rowPrimeTrialDivisionData->exp3 +
            (double)Dividend[2] * (double)rowPrimeTrialDivisionData->exp2 +
            (double)Dividend[1] * (double)rowPrimeTrialDivisionData->exp1;
          break;
        case 6:
          dRem = (double)Dividend[5] * (double)rowPrimeTrialDivisionData->exp5 +
            (double)Dividend[4] * (double)rowPrimeTrialDivisionData->exp4 +
            (double)Dividend[3] * (double)rowPrimeTrialDivisionData->exp3 +
            (double)Dividend[2] * (double)rowPrimeTrialDivisionData->exp2 +
            (double)Dividend[1] * (double)rowPrimeTrialDivisionData->exp1;
          break;
        case 5:
          dRem = (double)Dividend[4] * (double)rowPrimeTrialDivisionData->exp4 +
            (double)Dividend[3] * (double)rowPrimeTrialDivisionData->exp3 +
            (double)Dividend[2] * (double)rowPrimeTrialDivisionData->exp2 +
            (double)Dividend[1] * (double)rowPrimeTrialDivisionData->exp1;
          break;
        case 4:
          dRem = (double)Dividend[3] * (double)rowPrimeTrialDivisionData->exp3 +
            (double)Dividend[2] * (double)rowPrimeTrialDivisionData->exp2 +
            (double)Dividend[1] * (double)rowPrimeTrialDivisionData->exp1;
          break;
        case 3:
          dRem = (double)Dividend[2] * (double)rowPrimeTrialDivisionData->exp2 +
            (double)Dividend[1] * (double)rowPrimeTrialDivisionData->exp1;
          break;
        default:
          dRem = (double)Dividend[1] * (double)rowPrimeTrialDivisionData->exp1;
          break;
        }
        dCurrentPrime = (double)currentPrime;
        dRem += (double)Dividend[0];
        dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
        // Get its inverse
        inverseA = intModInv((int)dRem, currentPrime);
        switch (NumberLengthB)
        {
        case 7:
          dRem = (double)biAbsLinearCoeff[6] * (double)rowPrimeTrialDivisionData->exp6 +
            (double)biAbsLinearCoeff[5] * (double)rowPrimeTrialDivisionData->exp5 +
            (double)biAbsLinearCoeff[4] * (double)rowPrimeTrialDivisionData->exp4 +
            (double)biAbsLinearCoeff[3] * (double)rowPrimeTrialDivisionData->exp3 +
            (double)biAbsLinearCoeff[2] * (double)rowPrimeTrialDivisionData->exp2 +
            (double)biAbsLinearCoeff[1] * (double)rowPrimeTrialDivisionData->exp1;
          break;
        case 6:
          dRem = (double)biAbsLinearCoeff[5] * (double)rowPrimeTrialDivisionData->exp5 +
            (double)biAbsLinearCoeff[4] * (double)rowPrimeTrialDivisionData->exp4 +
            (double)biAbsLinearCoeff[3] * (double)rowPrimeTrialDivisionData->exp3 +
            (double)biAbsLinearCoeff[2] * (double)rowPrimeTrialDivisionData->exp2 +
            (double)biAbsLinearCoeff[1] * (double)rowPrimeTrialDivisionData->exp1;
          break;
        case 5:
          dRem = (double)biAbsLinearCoeff[4] * (double)rowPrimeTrialDivisionData->exp4 +
            (double)biAbsLinearCoeff[3] * (double)rowPrimeTrialDivisionData->exp3 +
            (double)biAbsLinearCoeff[2] * (double)rowPrimeTrialDivisionData->exp2 +
            (double)biAbsLinearCoeff[1] * (double)rowPrimeTrialDivisionData->exp1;
          break;
        case 4:
          dRem = (double)biAbsLinearCoeff[3] * (double)rowPrimeTrialDivisionData->exp3 +
            (double)biAbsLinearCoeff[2] * (double)rowPrimeTrialDivisionData->exp2 +
            (double)biAbsLinearCoeff[1] * (double)rowPrimeTrialDivisionData->exp1;
          break;
        case 3:
          dRem = (double)biAbsLinearCoeff[2] * (double)rowPrimeTrialDivisionData->exp2 +
            (double)biAbsLinearCoeff[1] * (double)rowPrimeTrialDivisionData->exp1;
          break;
        default:
          dRem = (double)biAbsLinearCoeff[1] * (double)rowPrimeTrialDivisionData->exp1;
          break;
        }
        dRem += (double)biAbsLinearCoeff[0];
        dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
        RemB = (int)dRem;
        if (positive)
        {
          RemB = currentPrime - RemB;
        }
        dRem = (double)inverseA * (double)(rowPrimeSieveData0->modsqrt + RemB) +
               (double)SieveLimit;
        dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
        rowPrimeSieveData->soln1 = (int)dRem;
      }
      do
      {                       // For each polynomial...
        if (congruencesFound >= matrixBLength /*|| factorSiqs != null*/)
        {
          if (nbrThreadFinishedPolySet > numberThreads*polySet)
          {
            break;
          }
          //synchronized(amodq)
          {
            nbrThreadFinishedPolySet++;
          }
          return;             // Another thread finished factoring.
        }
        if (onlyFactoring)
        {
          polynomialsSieved += 2;
        }
        /***************/
        /* Sieve stage */
        /***************/
        PerformSiqsSieveStage(primeSieveData, SieveArray,
          PolynomialIndex,
          biLinearCoeff,
          NumberLength);
        ValuesSieved += 2 * SieveLimit;
        /************************/
        /* Trial division stage */
        /************************/
        index2 = 2 * SieveLimit - 1;
        do
        {
          index2 -= 16;
          if (((SieveArray[index2+1] | SieveArray[index2+2] |
            SieveArray[index2+3] | SieveArray[index2+4] |
            SieveArray[index2+5] | SieveArray[index2+6] |
            SieveArray[index2+7] | SieveArray[index2+8] |
            SieveArray[index2+9] | SieveArray[index2+10] |
            SieveArray[index2+11] | SieveArray[index2+12] |
            SieveArray[index2+13] | SieveArray[index2+14] |
            SieveArray[index2+15] | SieveArray[index2+16]) & 0x8080) != 0)
          {
            for (i = 16; i>0; i--)
            {
              if ((SieveArray[index2 + i] & 0x80) != 0)
              {
                if (congruencesFound >= matrixBLength)
                {       // All congruences were found: stop sieving.
                  index2 = 0;
                  break;
                }
                SieveLocationHit(rowMatrixB,
                  rowMatrixBbeforeMerge,
                  index2 + i, primeSieveData,
                  rowPartials,
                  rowSquares,
                  biDividend, NumberLength, biT,
                  biLinearCoeff, biR, biU, biV,
                  indexFactorsA, FALSE);
                if (congruencesFound >= matrixBLength)
                {               // All congruences were found: stop sieving.
                  index2 = 0;
                  break;
                }
              }
              if (SieveArray[index2 + i] < 0)
              {
                if (congruencesFound >= matrixBLength)
                {       // All congruences were found: stop sieving.
                  index2 = 0;
                  break;
                }
                SieveLocationHit(rowMatrixB,
                  rowMatrixBbeforeMerge,
                  index2 + i, primeSieveData,
                  rowPartials,
                  rowSquares,
                  biDividend, NumberLength, biT,
                  biLinearCoeff, biR, biU, biV,
                  indexFactorsA, TRUE);
                if (congruencesFound >= matrixBLength)
                {               // All congruences were found: stop sieving.
                  index2 = 0;
                  break;
                }
              }
            }
          }
        } while (index2 > 0);
        /*******************/
        /* Next polynomial */
        /*******************/
        PolynomialIndex += 2;
      } while (PolynomialIndex <= lastPolynomial &&
        congruencesFound < matrixBLength);
    }
#if 0
  }
  catch (ArithmeticException ae)
  {
    synchronized(matrixB)
    {
      factorSiqs = null;
      matrixB.notify();
    }
  }
#endif
}
#if 0
/* Implementation of algorithm explained in Gower and Wagstaff paper */
/* The variables with suffix 3 correspond to multiplier = 3 */
static int SQUFOF(long N, int queue[])
{
  double sqrt;
  int Q, Q1, P, P1, L, S;
  int i, j, r, s, t, q;
  int queueHead, queueTail, queueIndex;
  long N3;
  int Q3, Q13, P3, P13, L3, S3;
  int r3, s3, t3, q3;
  int queueHead3, queueTail3, queueIndex3;
  int QRev, Q1Rev, PRev, P1Rev;
  int tRev, qRev, uRev;
  /* Step 1: Initialize */
  N3 = 3 * N;
  if ((N & 3) == 1)
  {
    N <<= 1;
  }
  if ((N3 & 3) == 1)
  {
    N3 <<= 1;
  }
  sqrt = sqrt(N);
  S = (int)sqrt;
  if ((long)(S + 1)*(long)(S + 1) <= N)
  {
    S++;
  }
  if ((long)S*(long)S > N)
  {
    S--;
  }
  if ((long)S*(long)S == N)
  {
    return S;
  }
  Q1 = 1;
  P = S;
  Q = (int)N - P*P;
  L = (int)(2 * sqrt(2 * sqrt));
  queueHead = 0;
  queueTail = 0;

  sqrt = sqrt(N3);
  S3 = (int)sqrt;
  if ((long)(S3 + 1)*(long)(S3 + 1) <= N3)
  {
    S3++;
  }
  if ((long)S3*(long)S3 > N3)
  {
    S3--;
  }
  if ((long)S3*(long)S3 == N3)
  {
    return S3;
  }
  Q13 = 1;
  P3 = S3;
  Q3 = (int)N3 - P3*P3;
  L3 = (int)(2 * sqrt(2 * sqrt));
  queueHead3 = 100;
  queueTail3 = 100;

  /* Step 2: Cycle forward to find a proper square form */
  for (i = 0; i <= L; i++)
  {
    /* Multiplier == 1 */
    q = (S + P) / Q;
    P1 = q*Q - P;
    if (Q <= L)
    {
      if ((Q & 1) == 0)
      {
        queue[queueHead++] = Q >> 1;
        queue[queueHead++] = P % (Q >> 1);
        if (queueHead == 100)
        {
          queueHead = 0;
        }
      }
      else if (Q + Q <= L)
      {
        queue[queueHead++] = Q;
        queue[queueHead++] = P % Q;
        if (queueHead == 100)
        {
          queueHead = 0;
        }
      }
    }
    t = Q1 + q*(P - P1);
    Q1 = Q;
    Q = t;
    P = P1;
    {
      r = (int)sqrt(Q);
      if (r*r == Q)
      {
        queueIndex = queueTail;
        for (;;)
        {
          if (queueIndex == queueHead)
          {
            /* Step 3: Compute inverse square root of the square form */
            PRev = P;
            Q1Rev = r;
            uRev = (S - PRev) % r;
            uRev += (uRev >> 31) & r;
            PRev = S - uRev;
            QRev = (int)((N - (long)PRev*(long)PRev) / Q1Rev);
            /* Step 4: Cycle in the reverse direction to find a factor of N */
            for (j = i; j >= 0; j--)
            {
              qRev = (S + PRev) / QRev;
              P1Rev = qRev*QRev - PRev;
              if (PRev == P1Rev)
              {
                /* Step 5: Get the factor of N */
                if ((QRev & 1) == 0)
                {
                  return QRev >> 1;
                }
                return QRev;
              }
              tRev = Q1Rev + qRev*(PRev - P1Rev);
              Q1Rev = QRev;
              QRev = tRev;
              PRev = P1Rev;
              qRev = (S + PRev) / QRev;
              P1Rev = qRev*QRev - PRev;
              if (PRev == P1Rev)
              {
                /* Step 5: Get the factor of N */
                if ((QRev & 1) == 0)
                {
                  return QRev >> 1;
                }
                return QRev;
              }
              tRev = Q1Rev + qRev*(PRev - P1Rev);
              Q1Rev = QRev;
              QRev = tRev;
              PRev = P1Rev;
            }
            break;
          }
          s = queue[queueIndex++];
          t = queue[queueIndex++];
          if (queueIndex == 100)
          {
            queueIndex = 0;
          }
          if ((P - t) % s == 0)
          {
            break;
          }
        }
        if (r > 1)
        {
          queueTail = queueIndex;
        }
        if (r == 1)
        {
          queueIndex = queueTail;
          for (;;)
          {
            if (queueIndex == queueHead)
            {
              break;
            }
            if (queue[queueIndex] == 1)
            {
              return 0;
            }
            queueIndex += 2;
            if (queueIndex == 100)
            {
              queueIndex = 0;
            }
          }
        }
      }
    }
    q = (S + P) / Q;
    P1 = q*Q - P;
    if (Q <= L)
    {
      if ((Q & 1) == 0)
      {
        queue[queueHead++] = Q >> 1;
        queue[queueHead++] = P % (Q >> 1);
        if (queueHead == 100)
        {
          queueHead = 0;
        }
      }
      else if (Q + Q <= L)
      {
        queue[queueHead++] = Q;
        queue[queueHead++] = P % Q;
        if (queueHead == 100)
        {
          queueHead = 0;
        }
      }
    }
    t = Q1 + q*(P - P1);
    Q1 = Q;
    Q = t;
    P = P1;

    /* Multiplier == 3 */
    q3 = (S3 + P3) / Q3;
    P13 = q3*Q3 - P3;
    if (Q3 <= L3)
    {
      if ((Q3 & 1) == 0)
      {
        queue[queueHead3++] = Q3 >> 1;
        queue[queueHead3++] = P3 % (Q3 >> 1);
        if (queueHead3 == 200)
        {
          queueHead3 = 100;
        }
      }
      else if (Q3 + Q3 <= L3)
      {
        queue[queueHead3++] = Q3;
        queue[queueHead3++] = P3 % Q3;
        if (queueHead3 == 200)
        {
          queueHead3 = 100;
        }
      }
    }
    t3 = Q13 + q3*(P3 - P13);
    Q13 = Q3;
    Q3 = t3;
    P3 = P13;
    {
      r3 = (int)sqrt(Q3);
      if (r3*r3 == Q3)
      {
        queueIndex3 = queueTail3;
        for (;;)
        {
          if (queueIndex3 == queueHead3)
          {
            /* Step 3: Compute inverse square root of the square form */
            PRev = P3;
            Q1Rev = r3;
            uRev = (S3 - PRev) % r3;
            uRev += (uRev >> 31) & r3;
            PRev = S3 - uRev;
            QRev = (int)((N3 - (long)PRev*(long)PRev) / Q1Rev);
            /* Step 4: Cycle in the reverse direction to find a factor of N */
            for (j = i; j >= 0; j--)
            {
              qRev = (S3 + PRev) / QRev;
              P1Rev = qRev*QRev - PRev;
              if (PRev == P1Rev)
              {
                /* Step 5: Get the factor of N */
                if ((QRev & 1) == 0)
                {
                  return QRev >> 1;
                }
                return QRev;
              }
              tRev = Q1Rev + qRev*(PRev - P1Rev);
              Q1Rev = QRev;
              QRev = tRev;
              PRev = P1Rev;
              qRev = (S3 + PRev) / QRev;
              P1Rev = qRev*QRev - PRev;
              if (PRev == P1Rev)
              {
                /* Step 5: Get the factor of N */
                if ((QRev & 1) == 0)
                {
                  return QRev >> 1;
                }
                return QRev;
              }
              tRev = Q1Rev + qRev*(PRev - P1Rev);
              Q1Rev = QRev;
              QRev = tRev;
              PRev = P1Rev;
            }
            break;
          }
          s3 = queue[queueIndex3++];
          t3 = queue[queueIndex3++];
          if (queueIndex3 == 200)
          {
            queueIndex3 = 100;
          }
          if ((P3 - t3) % s3 == 0)
          {
            break;
          }
        }
        if (r3 > 1)
        {
          queueTail3 = queueIndex3;
        }
        if (r3 == 1)
        {
          queueIndex3 = queueTail3;
          for (;;)
          {
            if (queueIndex3 == queueHead3)
            {
              break;
            }
            if (queue[queueIndex3] == 1)
            {
              return 0;
            }
            queueIndex3 += 2;
            if (queueIndex3 == 200)
            {
              queueIndex3 = 100;
            }
          }
        }
      }
    }
    q3 = (S3 + P3) / Q3;
    P13 = q3*Q3 - P3;
    if (Q3 <= L3)
    {
      if ((Q3 & 1) == 0)
      {
        queue[queueHead3++] = Q3 >> 1;
        queue[queueHead3++] = P3 % (Q3 >> 1);
        if (queueHead3 == 200)
        {
          queueHead3 = 100;
        }
      }
      else if (Q3 + Q3 <= L3)
      {
        queue[queueHead3++] = Q3;
        queue[queueHead3++] = P3 % Q3;
        if (queueHead3 == 200)
        {
          queueHead3 = 100;
        }
      }
    }
    t3 = Q13 + q3*(P3 - P13);
    Q13 = Q3;
    Q3 = t3;
    P3 = P13;
  }
  return 0;
}
/* If 2^value mod value = 2, then the value is a probable prime (value odd) */
static unsigned char isProbablePrime(long value)
{
  long mask, montgomery2, Pr, Prod0, Prod1, MontDig;
  int N, MontgomeryMultN;
  long BaseHI, BaseLO, valueHI, valueLO;
  int x = N = (int)value; // 2 least significant bits of inverse correct.
  x = x * (2 - N * x);     // 4 least significant bits of inverse correct.
  x = x * (2 - N * x);     // 8 least significant bits of inverse correct.
  x = x * (2 - N * x);     // 16 least significant bits of inverse correct.
  x = x * (2 - N * x);     // 32 least significant bits of inverse correct.
  MontgomeryMultN = (-x) & 0x7FFFFFFF;
  mask = 1L << 62;
  montgomery2 = 2 * (mask % value);
  if (montgomery2 >= value)
  {
    montgomery2 -= value;
  }
  BaseHI = (int)(montgomery2 >> 31);
  BaseLO = (int)montgomery2 & 0x7FFFFFFF;
  valueHI = (int)(value >> 31);
  valueLO = (int)value & 0x7FFFFFFF;
  while ((mask & value) == 0)
  {
    mask >>= 1;
  }
  mask >>= 1;
  while (mask > 0)
  {
    /* Square the base */
    Pr = BaseLO * BaseLO;
    MontDig = ((int)Pr * MontgomeryMultN) & 0x7FFFFFFFL;
    Prod0 = (Pr = ((MontDig * valueLO + Pr) >> 31) +
      MontDig * valueHI + BaseLO * BaseHI) & 0x7FFFFFFFL;
    Prod1 = Pr >> 31;
    Pr = BaseHI * BaseLO + Prod0;
    MontDig = ((int)Pr * MontgomeryMultN) & 0x7FFFFFFFL;
    Prod0 = (Pr = ((MontDig * valueLO + Pr) >> 31) +
      MontDig * valueHI + BaseHI * BaseHI + Prod1) & 0x7FFFFFFFL;
    Prod1 = Pr >> 31;
    if (Prod1 > valueHI || (Prod1 == valueHI && Prod0 >= valueLO))
    {
      Prod0 = (Pr = Prod0 - valueLO) & 0x7FFFFFFFL;
      Prod1 = ((Pr >> 31) + Prod1 - valueHI) & 0x7FFFFFFFL;
    }
    BaseLO = Prod0;
    BaseHI = Prod1;

    if ((mask & value) != 0)
    {
      /* Multiply by 2 */
      Pr = 2 * ((BaseHI << 31) + BaseLO);
      if (Pr >= value)
      {
        Pr -= value;
      }
      BaseHI = (int)(Pr >> 31);
      BaseLO = (int)Pr & 0x7FFFFFFF;
    }
    mask >>= 1;
  }
  Pr = (BaseHI << 31) + BaseLO;
  return Pr == montgomery2;
}

#endif

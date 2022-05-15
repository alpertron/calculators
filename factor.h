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
#ifndef _FACTOR_H
#define _FACTOR_H
#define MAX_FACTORS 5000
#include "showtime.h"
#ifdef __EMSCRIPTEN__
void getCunn(const char *url, char *factorsFromServer);
#endif

#define TYP_AURIF    100000000
#define TYP_TABLE    150000000
#define TYP_SIQS     200000000
#define TYP_LEHMAN   250000000
#define TYP_RABIN    300000000
#define TYP_DIVISION 350000000
#define TYP_EC       400000000

enum eEcmResult
{
  FACTOR_NOT_FOUND = 0,
  FACTOR_FOUND,
  FACTOR_NOT_FOUND_GCD,
  CHANGE_TO_SIQS,
};

struct sFactors
{
  int *ptrFactor;
  int multiplicity;
  int upperBound;
  int type;
};
#ifdef FACTORIZATION_APP
extern int StepECM;
#endif

extern int64_t primeModMult;
extern int64_t SIQSModMult;
extern int timePrimalityTests;
extern int timeSIQS;
extern int timeECM;
extern int nbrPrimalityTests;
extern int nbrSIQS;
extern int nbrECM;

extern int trialDivisions;
extern int smoothsFound;
extern int totalPartials;
extern int partialsFound;
extern uint64_t ValuesSieved;
extern int congruencesFound;
extern int polynomialsSieved;
extern int matrixRows;
extern int matrixCols;
extern int nbrPartials;

extern char tofactorDec[MAX_LEN*12];
extern bool prettyprint;
extern bool cunningham;
extern bool hexadecimal;
extern struct sFactors stFactors[MAX_FACTORS];
void factor(const BigInteger *toFactor, const int *number, int *factors, struct sFactors *pstFactors);
void factorExt(const BigInteger* toFactor, const int* number,
  int* factors, struct sFactors* pstFactors, char* pcKnownFactors);
void FactoringSIQS(const limb *pNbrToFactor, limb *pFactor);
#ifndef lang  
  extern bool lang;
#endif
extern int nbrToFactor[MAX_LEN];
extern struct sFactors astFactorsMod[MAX_FACTORS];
extern int factorsMod[20000];
void SendFactorizationToOutput(const struct sFactors *pstFactors, char **pptrOutput, bool doFactorization);
void Totient(BigInteger *result);
void NumFactors(BigInteger* result);
void MinFactor(BigInteger* result);
void MaxFactor(BigInteger* result);
void SumOfDivisors(BigInteger *result);
void NumberOfDivisors(BigInteger *result);
enum eEcmResult ecmCurve(int* pEC, int* pNextEC);
#ifdef FACTORIZATION_APP
char* ShowFactoredPart(const BigInteger* pNbr, const struct sFactors* pstFactors);
void ShowLowerText(void);
int BpswPrimalityTest(const BigInteger* pValue, const struct sFactors* pstFactors);
void batchEcmCallback(char** pptrOutput);
#endif
#endif

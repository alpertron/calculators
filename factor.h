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
#ifndef _FACTOR_H
#define _FACTOR_H
#define MAX_FACTORS 5000
#define FACTOR_ARRSIZE (2*MAX_FACTORS)
#include "showtime.h"
#ifdef __EMSCRIPTEN__
void getCunn(char *url, char *factorsFromServer);
#endif
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

extern char tofactorDec[MAX_LEN*12];
extern char prettyprint, cunningham, hexadecimal;
extern struct sFactors stFactors[MAX_FACTORS];
extern int *factorArr[FACTOR_ARRSIZE];
void factor(BigInteger *nbrToFactor, int *number, int *factors, struct sFactors *pstFactors);
void factorExt(BigInteger *nbrToFactor, int *number, int *factors, struct sFactors *pstFactors, char *pcKnownFactors);
void FactoringSIQS(limb *pNbrToFactor, limb *pFactor);
extern int lang;
extern int nbrToFactor[MAX_LEN];
extern struct sFactors astFactorsMod[MAX_FACTORS];
extern int factorsMod[20000];
void SendFactorizationToOutput(struct sFactors *pstFactors, char **pptrOutput, int doFactorization);
void Totient(BigInteger *result);
void SumOfDivisors(BigInteger *result);
void NumberOfDivisors(BigInteger *result);
#ifdef FACTORIZATION_APP
char *ShowFactoredPart(BigInteger *nbr, void *vFactors);
void ShowLowerText(void);
#endif
#endif

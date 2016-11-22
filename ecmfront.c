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
#include <emscripten.h>
extern long long lModularMult;
#endif
extern char *output;
static BigInteger tofactor;
static int nbrToFactor[MAX_LEN];
static int NbrFactorsMod = 0;
int groupLen = 6;
struct sFactors astFactorsMod[1000];
int factorsMod[10000];
static BigInteger factorValue;

#ifdef __EMSCRIPTEN__
void databack(char *data);
extern double originalTenthSecond;
int tenths(void);
void GetDHMSt(char **pptrText, int tenths);

#endif

void ecmFrontText(char *tofactorText, int doFactorization)
{
  char *ptrOutput;
  enum eExprErr rc;
  rc = ComputeExpression(tofactorText, 1, &tofactor);
  if (output == NULL)
  {
    output = (char *)malloc(1000000);
  }
  if (output == NULL)
  {
    return;   // Go out if cannot generate output string.
  }
  if (rc == EXPR_OK && doFactorization)
  {
    NumberLength = tofactor.nbrLimbs;
    CompressBigInteger(nbrToFactor, &tofactor);
#ifdef __EMSCRIPTEN__
    lModularMult = 0;
#endif
    factor(nbrToFactor, factorsMod, astFactorsMod);
    NbrFactorsMod = astFactorsMod[0].multiplicity;
  }
  ptrOutput = output;
  strcpy(output, "2<p>");
  ptrOutput += strlen(output);
  if (rc != EXPR_OK)
  {
    textError(output + 1, rc);
    ptrOutput = output + strlen(output);
  }
  else
  {
    Bin2Dec(tofactor.limbs, ptrOutput, tofactor.nbrLimbs, groupLen);
    ptrOutput += strlen(ptrOutput);
    if (doFactorization)
    {
      struct sFactors *pstFactors;
      int i = 0;
      pstFactors = &astFactorsMod[1];
      if (NbrFactorsMod == 1 && pstFactors->multiplicity == 1 &&
          (*pstFactors->ptrFactor > 1 || *(pstFactors->ptrFactor+1) > 1))
      {    // Do not show zero or one as prime.
        strcpy(ptrOutput, lang ? " es primo" : " is prime");
      }
      else
      {
        strcpy(ptrOutput, " = ");
        ptrOutput += strlen(ptrOutput);
        for (;;)
        {
          UncompressBigInteger(pstFactors->ptrFactor, &factorValue);
          Bin2Dec(factorValue.limbs, ptrOutput, factorValue.nbrLimbs, groupLen);
          ptrOutput += strlen(ptrOutput);
          if (pstFactors->multiplicity > 1)
          {
            strcpy(ptrOutput, "<sup>");
            ptrOutput += strlen(ptrOutput);
            int2dec(&ptrOutput, pstFactors->multiplicity);
            strcpy(ptrOutput, "</sup>");
            ptrOutput += strlen(ptrOutput);
          }
          if (++i == NbrFactorsMod)
          {
            break;
          }
          strcpy(ptrOutput, " &times; ");
          ptrOutput += strlen(ptrOutput);
          pstFactors++;
        }
      }
    }
  }
  strcpy(ptrOutput, lang ? "<p>Tiempo transcurrido: " : "<p>Time elapsed: ");
  ptrOutput += strlen(ptrOutput);
#ifdef __EMSCRIPTEN__
  GetDHMSt(&ptrOutput, (int)(tenths() - originalTenthSecond));
  strcpy(ptrOutput, "</p>");
  ptrOutput += strlen(ptrOutput);
#endif
  strcpy(ptrOutput, lang ? "<p>" COPYRIGHT_SPANISH "</p>" :
    "<p>" COPYRIGHT_ENGLISH "</p>");
}

#ifdef __EMSCRIPTEN__
void doWork(char* data, int size)
{
  int flags;
  char *ptrData = data;
  char *ptrPower, *ptrMod;
  if (output == NULL)
  {
    output = malloc(3000000);
  }
  groupLen = 0;
  while (*ptrData != ',')
  {
    groupLen = groupLen * 10 + (*ptrData++ - '0');
  }
  ptrData++;             // Skip comma.
  flags = *ptrData;
  lang = flags & 1;
  ecmFrontText(ptrData+2, flags & 2); // Skip app number and second comma.
  databack(output);
}
#endif

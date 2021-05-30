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
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "bignbr.h"
#include "highlevel.h"
#include "factor.h"
#include "batch.h"
#include "polynomial.h"
#ifndef DEBUG_CODE
#define DEBUG_CODE 13
#endif
void dilogText(char *baseText, char *powerText, char *modText, int groupLen);
void gaussianText(char *valueText, int doFactorization);
void ecmFrontText(char *tofactorText, bool doFactorization, char *knownFactors);
int Factor1[] = { 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x00 };
int Factor2[] = { 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x00 };
int Factor3[] = { 29504, 29490, 19798, 633, 181, 0, 0, 0, 0, 0 };
int Factor4[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
int Factor5[] = { 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0, 0, 0, 0 };
int Factor6[] = { 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0, 0, 0, 0 };
int factor7[2002];
int factor8[2002];
int factors[5000];
int quotientPoly[20];
struct sFactors astFactors[1000];
extern int number[MAX_LEN];
extern int nbrLimbs;
extern int groupLen;
extern limb TestNbr[MAX_LEN];
char expr[] = "123456789012345";
int Product[32];
char input[MAX_LEN*4];
extern char tofactorDec[MAX_LEN*12];
extern int app;
extern bool skipPrimality;
BigInteger dividend;
BigInteger divisor;
BigInteger quotient;
int main(int argc, char* argv[])
{
  (void)argc;  // Parameter is not used. 
#if DEBUG_CODE == 1
//  fsquaresText("n(10^32)", 6);
  fsquaresText(argv[1], 6);
  printf("%s\n", output);
#elif DEBUG_CODE == 2
  app = 1;
  fcubesText(argv[1], 6);
  printf("%s\n", output);
#elif DEBUG_CODE == 3
  int i, len;
  for (i = 0; i < 20; i++)
  {
    Factor1[i] = 7 - i;
    Factor2[i] = 1 + i;
    Product[i] = 21 + i;
  }
  multiply(Factor5, Factor6, Product, 4, NULL);
  for (i = 0; i < 32; i++)
  {
    printf("%04X ", Product[i]);
  }
  printf("\n");
  Factor1[0] = 0;
  Factor1[1] = 0;
  Factor1[2] = 0;
  Factor1[3] = 0;
  Factor1[4] = 728;
  Factor1[5] = 32767;
  Factor1[6] = 32767;
  squareRoot(&Factor1[2], Factor2, 4, &len);
  //squareRoot(Factor1, Factor2, 7, &len);
  while (1);
#elif DEBUG_CODE == 4
  char* ptrOutput;
  if (argv[1][0] == '-')
  {
    Dec2Bin(&argv[1][1], dividend.limbs, (int)strlen(&argv[1][1]), &dividend.nbrLimbs);
    dividend.sign = SIGN_NEGATIVE;
  }
  else
  {
    Dec2Bin(argv[1], dividend.limbs, (int)strlen(argv[1]), &dividend.nbrLimbs);
    dividend.sign = SIGN_POSITIVE;
  }
  if (argv[2][0] == '-')
  {
    Dec2Bin(&argv[2][1], divisor.limbs, (int)strlen(&argv[2][1]), &divisor.nbrLimbs);
    divisor.sign = SIGN_NEGATIVE;
  }
  else
  {
    Dec2Bin(argv[2], divisor.limbs, (int)strlen(argv[2]), &divisor.nbrLimbs);
    divisor.sign = SIGN_POSITIVE;
  }
  // Insert garbage after dividend and divisor.
  (void)memset(&dividend.limbs[dividend.nbrLimbs], 0x45, 12);
  (void)memset(&divisor.limbs[divisor.nbrLimbs], 0x25, 12);
  BigIntDivide(&dividend, &divisor, &quotient);
  ptrOutput = output;
  Bin2Dec(&ptrOutput, quotient.limbs, quotient.nbrLimbs, 0);
  if (quotient.sign == SIGN_POSITIVE)
  {
    printf("Quotient = %s\n", output);
  }
  else
  {
    printf("Quotient = -%s\n", output);
  }
#elif DEBUG_CODE == 5
  int i;
  for (i = sizeof(expr) - 1; i >= 0; i -= 3)
  {
    expr[i] = 0;
    fcubesText(expr, 6);
    printf("%s\n\n", output);
  }
#elif DEBUG_CODE == 6
  char* ptrInput;
  if (argc != 4)
  {
    printf("num delta den\n");
    return 0;
  }
  ptrInput = input;
  copyStr(&ptrInput, argv[1]);
  ptrInput++;
  copyStr(&ptrInput, argv[2]);
  ptrInput++;
  copyStr(&ptrInput, argv[3]);
  ptrInput++;
  contfracText(input, 20000);
  printf("%s\n", output);
#elif DEBUG_CODE == 7
  char* ptrInput;
  if (argc != 2)
  {
    printf("num\n");
    return 0;
  }
  ptrInput = input;
  (void)strcpy(ptrInput, argv[1]);
  fsquaresText(input, 6);
  printf("%s\n", output);
#elif DEBUG_CODE == 8
  char* ptrInput;
  char* ptrOutput;
  enum eExprErr rc;
  int NumberLength;
  BigInteger num, mod, inv;

  if (argc != 3)
  {
    printf("num mod\n");
    return 0;
  }
  rc = ComputeExpression(argv[1], &num);
  if (rc != EXPR_OK)
  {
    ptrOutput = output;
    textError(&ptrOutput, rc);
  }
  else
  {
    rc = ComputeExpression(argv[2], &mod);
    if (rc != EXPR_OK)
    {
      ptrOutput = output;
      textError(&ptrOutput, rc);
    }
    else
    {
      NumberLength = mod.nbrLimbs;
      if (num.nbrLimbs < mod.nbrLimbs)
      {
        (void)memset(&num.limbs[num.nbrLimbs], 0, (mod.nbrLimbs - num.nbrLimbs) * sizeof(limb));
      }
      (void)memcpy(TestNbr, mod.limbs, NumberLength * sizeof(limb));
      TestNbr[NumberLength].x = 0;
      GetMontgomeryParms(NumberLength);
      ModInvBigNbr(num.limbs, inv.limbs, mod.limbs, NumberLength);
      ptrOutput = output;
      Bin2Dec(&ptrOutput, inv.limbs, NumberLength, 200);
    }
  }
  printf("%s", output);
#elif DEBUG_CODE == 9
  if (argc != 3)
  {
    (void)printf("modulus polynomial\n");
    return 0;
  }
  pretty = PRETTY_PRINT;
  for (int ctr = 0; ctr < 1; ctr++)
  {
    polyFactText(argv[1], argv[2], 7);
  }
  (void)printf("%s\n", output);
#elif DEBUG_CODE == 11
  if (argc != 4)
  {
    printf("base power modulus\n");
    return 0;
  }
  dilogText(argv[1], argv[2], argv[3], 6);
  dilogText(argv[1], argv[2], argv[3], 6);
  printf("%s\n", output);
#elif DEBUG_CODE == 12
  if (argc != 3)
  {
    printf("value factorize\n");
    return 0;
  }
  gaussianText(argv[1], argv[2][0]);
  printf("%s\n", output);
#elif DEBUG_CODE == 13
  skipPrimality = false;
  lang = false;
  hexadecimal = false;
  char text[40000];
  char* ptrText = text;
  copyStr(&ptrText, argv[1]);
  copyStr(&ptrText, "\n");
  ecmFrontText(text, true, NULL);
  (void)printf("%s\n", output);
#if 0
  if (argc == 3)
  {
    ecmFrontText(argv[1], true, argv[2]);
    printf("%s\n", output);
  }
  else if (argc == 2)
  {
    char *ptrKnownFactors = strchr(argv[1], '=');
#if 0
    (void)strcpy(text, "x=10**45+572;x=x+1;c<1000;x");
    ecmFrontText(text, true, ptrKnownFactors);
    printf("%s\n", output);
    ecmFrontText(NULL, false, NULL);
    printf("%s\n", output);
    valuesProcessed = 0;
    (void)strcpy(text, "10**45+573");
    ecmFrontText(text, true, NULL);
    printf("%s\n", output);
    ecmFrontText(NULL, false, NULL);
    printf("%s\n", output);
    return 0;
#endif
    if (ptrKnownFactors != NULL)
    {                          // There is equal sign.
      *ptrKnownFactors = 0;    // Replace equal sign by string terminator.
      ptrKnownFactors++;
    }
    (void)sprintf(text, "%s\n", argv[1]);
    ecmFrontText(text, true, ptrKnownFactors);
    printf("%s\n", output);
  }
  else
  {
    printf("value [known factors]\n");
    return 0;
  }
#endif
#elif DEBUG_CODE == 14
  {
    limb internalNotation[100];
    static int bitGroups;
    char textInput[500];
    char* ptrTextInput = textInput;
    textInput[0] = '1';
    (void)memset(&textInput[1], '0', 150);
    textInput[151] = 0;
    Dec2Bin(textInput, internalNotation, 151, &bitGroups);
    Bin2Dec(&ptrTextInput, internalNotation, 17, 6);
    textInput[200] = 0;
  }
#elif DEBUG_CODE == 15
  (void)memcpy(inputString, "6,-2,00102^1042+1""\0""2^1042+1=5^1(0)*16673^1(0)*627186185377^1(16673)*131294792925870751515684960383613518415615538737991528767912593379854404518341858118366491474959205710499826133822402120149306175263402700301^1(16673)*6864797660130609714981900799081393217269435300143305409394463459185543183397652346775704046543201000705776033378429553397612687501667381169885775070966579201^1(2)""\0\0""222""\0", 10001314 - 10000928 + 1);
  doWork();
#elif DEBUG_CODE == 16
  quadmodText(argv[1], argv[2], argv[3], argv[4], 6);
  printf("%s\n", output);
#elif DEBUG_CODE == 17
  if (argc != 7)
  {
    printf("6 arguments needed\n");
    return 1;
  }
  quadText(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);
  printf("%s\n", output);
#elif DEBUG_CODE == 18
  int resultLen, k;
  (void)memset(factors, 0xEE, 1000 * sizeof(limb));
  for (k = 0; k < 500; k++)
  {
    factor7[k] = 0x7FFFFFFF;
  }
  (void)memset(factors, 0x00, 2000 * sizeof(limb));
  fftMultiplication((limb *)factor7, (limb *)factor7, (limb *)factors, 4, &resultLen);
#elif DEBUG_CODE == 19
  limb tempVal[4];
  limb tempRes[4];
  char* ptrOutput;
  // TestNbr = 5^18
  TestNbr[0].x = 1128244537;
  TestNbr[1].x = 441554605;
  TestNbr[2].x = 1;
  TestNbr[3].x = 0;
  NumberLength = 3;
  GetMontgomeryParms(NumberLength);
  tempVal[0].x = 312571964;
  tempVal[1].x = 1720066743;
  tempVal[2].x = 0;
  tempVal[3].x = 0;
  ptrOutput = output;
  Bin2Dec(&ptrOutput, tempVal, 3, 0);
  printf("Value = %s\n", output);
  ptrOutput = output;
  Bin2Dec(&ptrOutput, MontgomeryMultR1, NumberLength, 0);
  printf("MontgomeryMultR1 = %s\n", output);
  ptrOutput = output;
  Bin2Dec(&ptrOutput, TestNbr, NumberLength, 0);
  printf("TestNbr = %s\n", output);
  ModInvBigNbr(tempVal, tempRes, TestNbr, 3);
  ptrOutput = output;
  Bin2Dec(&ptrOutput, tempRes, 3, 0);
  printf("Inverse = %s", output);
#elif DEBUG_CODE == 20
  int ctr;
  for (ctr = 0; ctr < 901; ctr++)
  {
    factor7[2 * ctr] = 1;
    factor7[2 * ctr + 1] = 10;
    factor8[2 * ctr] = 1;
    factor8[2 * ctr + 1] = 10;
  }
  fftPolyMult(factor7, factor8, polyMultTemp, 127, 64);
  for (ctr = 0; ctr < 501; ctr++)
  {
    factor7[2 * ctr] = 1;
    factor7[2 * ctr + 1] = 10;
    factor8[2 * ctr] = 1;
    factor8[2 * ctr + 1] = 10;
  }
  fftPolyMult(factor7, factor8, polyMultTemp, 501, 501);
#elif DEBUG_CODE == 21
//  int dividendPoly[] = {1, 15, 1, 14, 1, 13, 1, 12, 1, 11, 1, 10, 1, 9, 1, 8,
//    1, 7, 1, 6, 1, 5, 1, 4, 1, 3, 1, 2, 1, 1};
  int dividendPoly[] = { 1, 7, 1, 6, 1, 5, 1, 4, 1, 3, 1, 2, 1, 1 };
  int divisorPoly[] = { 1, 11, 1, 9, 1, 7, 1, 5, 1, 3 };
  int dividendDegree = 6;
  int divisorDegree = 4;
  NumberLength = 1;
  TestNbr[0].x = 23;
  GetMontgomeryParms(1);
  DividePolynomial(dividendPoly, dividendDegree, divisorPoly, divisorDegree, quotientPoly);
#endif
  return 0;
}
#ifdef __EMSCRIPTEN__
void databack(const char* data)
{
  (void)printf("%s\n", data);
}

double tenths(void)
{
  return 0;
}

void startSkipTest(void)
{
  // Nothing to do in this stub implementation.
}

void endSkipTest(void)
{
  // Nothing to do in this stub implementation.
}

void getCunn(const char* url, char* factorsFromServer)
{
  (void)url;
  (void)factorsFromServer;
  // Nothing to do in this stub implementation.
}

#endif


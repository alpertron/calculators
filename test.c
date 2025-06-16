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
#include "expression.h"
#include "highlevel.h"
#include "factor.h"
#include "batch.h"
#include "polynomial.h"
#include "fromBlockly.h"
#include "isprime.h"
#ifndef DEBUG_CODE
#define DEBUG_CODE 13
#endif
#ifdef __EMSCRIPTEN__
extern char inputString[];
#endif
#if (DEBUG_CODE == 9) || (DEBUG_CODE == 17)
extern bool teach;
#endif
#if DEBUG_CODE == 13
extern bool doShowPrime;
#endif
#if DEBUG_CODE == 22
static char bufferXML[1000000];
#endif
void showSumTwoSquares(void);
void showDivisors(void);
void dilogText(const char *baseText, const char *powerText, const char *modText, int groupLen);
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
char results[300000];
struct sFactors astFactors[1000];
extern int number[MAX_LEN];
extern int nbrLimbs;
extern int groupLen;
extern limb TestNbr[MAX_LEN];
extern limb MontgomeryMultN[MAX_LEN];
char expr[] = "123456789012345";
int Product[32];
char input[MAX_LEN*4];
extern char tofactorDec[MAX_LEN*12];
extern bool skipPrimality;
BigInteger dividend;
BigInteger divisor;
BigInteger quotient;
int main(int argc, char* argv[])
{
  (void)argc;  // Parameter is not used. 
#if DEBUG_CODE == 1
  fsquaresText(argv[1], 6);
  (void)printf("%s\n", output);
#elif DEBUG_CODE == 2
  fcubesText(argv[1], 6);
  (void)printf("%s\n", output);
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
    (void)printf("%04X ", Product[i]);
  }
  (void)printf("\n");
  Factor1[0] = 0;
  Factor1[1] = 0;
  Factor1[2] = 0;
  Factor1[3] = 0;
  Factor1[4] = 728;
  Factor1[5] = 32767;
  Factor1[6] = 32767;
  squareRoot(&Factor1[2], Factor2, 4, &len);
  return 0;
#elif DEBUG_CODE == 4
  char* ptrOutput;
  if (argc != 3)
  {
    printf("Arguments: dividend divisor");
  }
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
    (void)printf("Quotient = %s\n", output);
  }
  else
  {
    (void)printf("Quotient = -%s\n", output);
  }
  BigIntRemainder(&dividend, &divisor, &quotient);
  ptrOutput = output;
  Bin2Dec(&ptrOutput, quotient.limbs, quotient.nbrLimbs, 0);
  if (quotient.sign == SIGN_POSITIVE)
  {
    (void)printf("Remainder = %s\n", output);
  }
  else
  {
    (void)printf("Remainder = -%s\n", output);
  }
#elif DEBUG_CODE == 5
  int i;
  for (i = sizeof(expr) - 1; i >= 0; i -= 3)
  {
    expr[i] = 0;
    fcubesText(expr, 6);
    (void)printf("%s\n\n", output);
  }
#elif DEBUG_CODE == 6
  char* ptrInput;
  if (argc != 6)
  {
    (void)printf("num delta den show_hex show_convergents\n");
    return 0;
  }
  ptrInput = input;
  copyStr(&ptrInput, argv[1]);
  ptrInput++;
  copyStr(&ptrInput, argv[2]);
  ptrInput++;
  copyStr(&ptrInput, argv[3]);
  ptrInput++;
  bool hex = (argv[4][0] != '0');
  bool converg = (argv[5][0] != '0');
  contfracText(input, 6, hex, converg);
  (void)printf("%s\n", output);
#elif DEBUG_CODE == 7
  char* ptrInput;
  if (argc != 2)
  {
    (void)printf("num\n");
    return 0;
  }
  ptrInput = input;
  (void)strcpy(ptrInput, argv[1]);
  fsquaresText(input, 6);
  (void)printf("%s\n", output);
#elif DEBUG_CODE == 8
  char* ptrInput;
  char* ptrOutput;
  enum eExprErr rc;
  int NumberLength;
  BigInteger num, mod, inv;

  if (argc != 3)
  {
    (void)printf("num mod\n");
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
      (void)ModInvBigNbr(num.limbs, inv.limbs, mod.limbs, NumberLength);
      ptrOutput = output;
      Bin2Dec(&ptrOutput, inv.limbs, NumberLength, 200);
    }
  }
  (void)printf("%s", output);
#elif DEBUG_CODE == 9
  if (argc != 5)
  {
    (void)printf("modulus polynomial teach pretty (0: pretty, 1: TEX, 2: Pari-GP)\n");
    return 0;
  }
  teach = (argv[3][0] != '0');
  switch (argv[4][0] - '0')
  {
  case 0:
    pretty = PRETTY_PRINT;
    break;
  case 1:
    pretty = TEX;
    break;
  default:
    pretty = PARI_GP;
    break;
  }
  for (int ctr = 0; ctr < 2; ctr++)
  {
    polyFactText(argv[1], argv[2], 7);
  }
  (void)printf("%s\n", output);
#elif DEBUG_CODE == 11
  if (argc != 4)
  {
    (void)printf("base power modulus\n");
    return 0;
  }
  dilogText(argv[1], argv[2], argv[3], 6);
  dilogText(argv[1], argv[2], argv[3], 6);
  (void)printf("%s\n", output);
#elif DEBUG_CODE == 12
  if (argc != 3)
  {
    (void)printf("value factorize\n");
    return 0;
  }
  gaussianText(argv[1], argv[2][0]);
  (void)printf("%s\n", output);
#elif DEBUG_CODE == 13
  skipPrimality = false;
  lang = false;
  hexadecimal = false;
  char text[40000];
  char* ptrText = text;
  char* ptrCurveNbr = NULL;
  bool doFactorization = true;
  if ((argc < 2) || (argc > 4))
  {
    (void)printf("nbr_to_factor [flags] [dofact]\n"
      "flags = 1: hexadecimal output, do not use SIQS\n"
      "flags = 2: decimal output, use SIQS\n"
      "flags = 3: hexadecimal output, use SIQS\n"
      "flags = 4: decimal output, do not use SIQS\n"
      "dofact = 0: do not factor, just show expression result\n"
      "dofact = 1: factor expression (default)\n"
      "dofact = 2: show whether the expression is prime\n");
    return 0;
  }
  if (argc == 3)
  {
    if ((argv[2][0] == '1') || (argv[2][0] == '3'))
    {
      hexadecimal = true;
    }
    if ((argv[2][0] == '2') || (argv[2][0] == '3'))
    {
      ptrCurveNbr = "!";    // Use SIQS.
    }
  }
  if (argc == 4)
  {
    if (argv[3][0] == '0')
    {
      doFactorization = false;
    }
    if (argv[3][0] == '2')
    {
      doShowPrime = true;
      doFactorization = false;
    }
  }
  copyStr(&ptrText, argv[1]);
  copyStr(&ptrText, "\n");
  ecmFrontText(text, doFactorization, ptrCurveNbr);
  (void)printf("%s\n", output);
#if 0
  if (argc == 3)
  {
    ecmFrontText(argv[1], true, argv[2]);
    (void)printf("%s\n", output);
  }
  else if (argc == 2)
  {
    char* ptrKnownFactors = strchr(argv[1], '=');
    if (ptrKnownFactors != NULL)
    {                          // There is equal sign.
      *ptrKnownFactors = 0;    // Replace equal sign by string terminator.
      ptrKnownFactors++;
    }
    (void)sprintf(text, "%s\n", argv[1]);
    ecmFrontText(text, true, ptrKnownFactors);
    (void)printf("%s\n", output);
  }
  else
  {
    (void)printf("value [known factors]\n");
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
  if (argc != 5)
  {
    (void)printf("Enter: quadr linear const modulus\n");
    return 1;
  }
  quadmodText(argv[1], argv[2], argv[3], argv[4], 6);
  quadmodText(argv[1], argv[2], argv[3], argv[4], 6);
  (void)printf("%s\n", output);
#elif DEBUG_CODE == 17
#if 1
  if (argc != 8)
  {
    (void)printf("Enter 6 coefficients and teach flag (0 or 1):\n"
      "   x^2, xy, y^2, x, y, const teach.\n");
    return 1;
  }
  if (argv[7][0] == '1')
  {
    teach = true;
  }
  else
  {
    teach = false;
  }
  quadText(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);
  (void)printf("%s\n", output);
#else
  for (int k = 0; k < 100; k++)
  {
    char buf1[100];
    char buf2[100];
    char buf3[100];
    char buf4[100];
    char buf5[100];
    char buf6[100];
    sprintf(buf1, "%d", 543253452 + 7 * k);
    sprintf(buf2, "%d", 1547654759 + 11 * k);
    sprintf(buf3, "%d", 945432534 + 13 * k);
    sprintf(buf4, "%d", 743134318 + 17 * k);
    sprintf(buf5, "%d", 142435467 + 19 * k);
    sprintf(buf6, "%d", 854325454 + 23 * k);
    teach = false;
    quadText(buf1, buf2, buf3, buf4, buf5, buf6);
  }
#endif
#elif DEBUG_CODE == 18
  int resultLen, k;
  (void)memset(factors, 0xEE, 1000 * sizeof(limb));
  for (k = 0; k < 500; k++)
  {
    factor7[k] = 0x7FFFFFFF;
  }
  (void)memset(factors, 0x00, 2000 * sizeof(limb));
  fftMultiplication((limb*)factor7, (limb*)factor7, (limb*)factors, 4, &resultLen);
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
  (void)printf("Value = %s\n", output);
  ptrOutput = output;
  Bin2Dec(&ptrOutput, MontgomeryMultR1, NumberLength, 0);
  (void)printf("MontgomeryMultR1 = %s\n", output);
  ptrOutput = output;
  Bin2Dec(&ptrOutput, TestNbr, NumberLength, 0);
  (void)printf("TestNbr = %s\n", output);
  (void)ModInvBigNbr(tempVal, tempRes, TestNbr, 3);
  ptrOutput = output;
  Bin2Dec(&ptrOutput, tempRes, 3, 0);
  (void)printf("Inverse = %s", output);
#elif DEBUG_CODE == 20
  int ctr;
  for (ctr = 0; ctr < 901; ctr++)
  {
    factor7[2 * ctr] = 1;
    factor7[2 * ctr + 1] = 10;
    factor8[2 * ctr] = 1;
    factor8[2 * ctr + 1] = 10;
  }
  fftPolyMult(factor7, factor8, common.poly.polyMultTemp, 127, 64);
  for (ctr = 0; ctr < 501; ctr++)
  {
    factor7[2 * ctr] = 1;
    factor7[2 * ctr + 1] = 10;
    factor8[2 * ctr] = 1;
    factor8[2 * ctr + 1] = 10;
  }
  fftPolyMult(factor7, factor8, common.poly.polyMultTemp, 501, 501);
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
#elif DEBUG_CODE == 22
  if (argc == 3)
  {
    int lineNbr;
    FILE* fpFile = fopen(argv[1], "rb");
    if (fpFile == NULL)
    {
      printf("Cannot open file.\n");
      return 1;
    }
    if (sscanf(argv[2], "%d", &lineNbr) != 1)
    {
      printf("Second parameter must be an integer.\n");
      fclose(fpFile);
      return 1;
    }
    for (int ctr = 0; ctr <= lineNbr; ctr++)
    {
      const char* rc = fgets(bufferXML, sizeof(bufferXML), fpFile);
      if (rc == NULL)
      {
        printf("Line not found in XML file.\n");
        fclose(fpFile);
        return 1;
      }
    }
    for (;;)
    {
      int len = (int)strlen(bufferXML) - 1;
      if ((bufferXML[len] == 0) || (bufferXML[len] & 0xE0))
      {
        break;
      }
      bufferXML[len] = 0;     // Delete line feed.
      len--;
    }
    fclose(fpFile);
    fromBlockly(bufferXML);
  }
  else
  {
    printf("Command line: calculator <XML file> <line number>.\n");
  }
#elif DEBUG_CODE == 23
  if (argc != 3)
  {
    printf("Command line: calculator number expon (expon must be odd)");
  }
  else
  {
    tsqcubesText(argv[1], 6, atoi(argv[2]));
    (void)printf("%s\n", output);
  }
#elif DEBUG_CODE == 24
  limb product[501];
  bool linePrinted = false;
  for (int limbNbr = 0; limbNbr < 500; limbNbr++)
  {
    TestNbr[limbNbr].x = 0;
    product[limbNbr + 1].x = 0x32423424;
    for (int ctr = 0; ctr < 7; ctr++)
    {
      TestNbr[limbNbr].x = TestNbr[limbNbr].x * 16 + 5;
      GetMontgomeryParms(limbNbr + 1);
      modmult(MontgomeryMultR1, MontgomeryMultR1, product);
      int nbrBytes = NumberLength * (int)sizeof(limb);
      if (memcmp(MontgomeryMultR1, product, nbrBytes) != 0)
      {
        (void)printf("Product is not equal to MontgomeryMultR1 in limbNbr = %d, ctr = %d\n", limbNbr, ctr);
        linePrinted = true;
        break;
      }
      if (product[limbNbr + 1].x != 0x32423424)
      {
        (void)printf("Overwriting after product when limbNbr = %d, ctr = %d\n", limbNbr, ctr);
        linePrinted = true;
        break;
      }
    }
    if (linePrinted)
    {
      break;
    }
  }
  if (!linePrinted)
  {
    (void)printf("No errors found in modmult test\n");
  }
#elif DEBUG_CODE == 25
  if (argc != 2)
  {
    printf("sumSquares {number}\n");
    return 1;
  }
  ecmFrontText(argv[1], true, NULL);
  showSumTwoSquares();
  (void)printf("%s\n", output);
#elif DEBUG_CODE == 26
  char buffer[100];
  strcpy(buffer, "422500");
  ecmFrontText(buffer, true, NULL);
  showSumTwoSquares();
  strcpy(buffer, "422500*17^2");
  ecmFrontText(buffer, true, NULL);
  showSumTwoSquares();
  (void)printf("%s\n", output);
#elif DEBUG_CODE == 27
  if (argc != 2)
  {
    printf("divisors {number}\n");
    return 1;
  }
  ecmFrontText(argv[1], true, NULL);
  showDivisors();
  (void)printf("%s\n", output);
#elif DEBUG_CODE == 28
  char* ptrResults = results;
  for (int cont = 1; cont < argc; cont++)
  {
    int value[2];
    getValue64(argv[cont], &value[0], &value[1]);
    ptrResults = appendInt(ptrResults, (isPrime(value) ? -109 : 0));
  }
  *ptrResults = 0;
  (void)printf("%s\n", results);
#elif DEBUG_CODE == 29
  static BigInteger multiplier;
  static BigInteger multiplicand;
  static BigInteger modulus;
  if (argc != 4)
  {
    printf("modmult {multiplicand} {multiplier} {modulus}\n");
    return 1;
  }
  int rc = ComputeExpression(argv[1], &multiplicand);
  if (rc != EXPR_OK)
  {
    printf("Error in the multiplicand: %s\n", argv[1]);
    return 1;
  }
  rc = ComputeExpression(argv[2], &multiplier);
  if (rc != EXPR_OK)
  {
    printf("Error in the multiplier: %s\n", argv[2]);
    return 1;
  }
  rc = ComputeExpression(argv[3], &modulus);
  if (rc != EXPR_OK)
  {
    printf("Error in the modulus: %s.\n", argv[3]);
    return 1;
  }
  int nbrBytes = (int)sizeof(limb) * modulus.nbrLimbs;
  memcpy(TestNbr, modulus.limbs, nbrBytes);
  GetMontgomeryParms(modulus.nbrLimbs);
  if (multiplicand.nbrLimbs < NumberLength)
  {
    nbrBytes = (NumberLength - multiplicand.nbrLimbs) * (int)sizeof(limb);
    memset(&multiplier.limbs[multiplicand.nbrLimbs], 0, nbrBytes);
  }
  if (multiplier.nbrLimbs < NumberLength)
  {
    nbrBytes = (NumberLength - multiplier.nbrLimbs) * (int)sizeof(limb);
    memset(&multiplier.limbs[multiplier.nbrLimbs], 0, nbrBytes);
  }
  // Convert to Montgomery notation.
  modmult(multiplicand.limbs, MontgomeryMultR2, multiplicand.limbs);
  modmult(multiplier.limbs, MontgomeryMultR2, multiplier.limbs);
  // Multiply in Montgomery notation.
  modmult(multiplicand.limbs, multiplier.limbs, multiplier.limbs);
  // Convert to normal notation.
  nbrBytes = NumberLength * (int)sizeof(limb);
  memset(multiplicand.limbs, 0, nbrBytes);
  multiplicand.limbs[0].x = 1;
  modmult(multiplicand.limbs, multiplier.limbs, multiplier.limbs);
  while (NumberLength > 1)
  {
    if (&multiplier.limbs[NumberLength - 1].x != 0)
    {
      break;
    }
    NumberLength--;
  }
  char* ptrResults = results;
  Bin2Dec(&ptrResults, multiplier.limbs, NumberLength, 0);
  *ptrResults = 0;     // Add string terminator.
  (void)printf("%s\n", results);
#elif DEBUG_CODE == 30
  static BigInteger multiplier;
  static BigInteger multiplicand;
  static BigInteger modulus;
  static BigInteger power;
  intToBigInteger(&multiplier, 10);
  BigIntPowerIntExp(&multiplier, 100, &power);
  for (int exponent = 100; exponent < 100000; exponent += 10)
  {
    if (exponent % 100 == 0)
    {
      fprintf(stderr, "%d\n", exponent);
    }
    intToBigInteger(&multiplicand, 3);
    intToBigInteger(&multiplier, 6);
    intToBigInteger(&modulus, 18745);
    BigIntAdd(&power, &modulus, &modulus);
    multint(&power, &power, 100000);
    multint(&power, &power, 100000);
    int nbrBytes = (int)sizeof(limb) * modulus.nbrLimbs;
    memcpy(TestNbr, modulus.limbs, nbrBytes);
    GetMontgomeryParms(modulus.nbrLimbs);
    if (multiplicand.nbrLimbs < NumberLength)
    {
      nbrBytes = (NumberLength - multiplicand.nbrLimbs) * (int)sizeof(limb);
      memset(&multiplier.limbs[multiplicand.nbrLimbs], 0, nbrBytes);
    }
    if (multiplier.nbrLimbs < NumberLength)
    {
      nbrBytes = (NumberLength - multiplier.nbrLimbs) * (int)sizeof(limb);
      memset(&multiplier.limbs[multiplier.nbrLimbs], 0, nbrBytes);
    }
    // Convert to Montgomery notation.
    modmult(multiplicand.limbs, MontgomeryMultR2, multiplicand.limbs);
    modmult(multiplier.limbs, MontgomeryMultR2, multiplier.limbs);
    // Multiply in Montgomery notation.
    modmult(multiplicand.limbs, multiplier.limbs, multiplier.limbs);
    // Convert to normal notation.
    nbrBytes = NumberLength * (int)sizeof(limb);
    memset(multiplicand.limbs, 0, nbrBytes);
    multiplicand.limbs[0].x = 1;
    modmult(multiplicand.limbs, multiplier.limbs, multiplier.limbs);
    while (NumberLength > 1)
    {
      if (&multiplier.limbs[NumberLength - 1].x != 0)
      {
        break;
      }
      NumberLength--;
    }
    char* ptrResults = results;
    Bin2Dec(&ptrResults, multiplier.limbs, NumberLength, 0);
    *ptrResults = 0;     // Add string terminator.
    if (strcmp(results, "18") != 0)
    {                    // Modular multiplication is incorrect.
      (void)printf("Exponent = %d, Result = %s\n", exponent, results);
    }
  }
#elif DEBUG_CODE == 31
  static BigInteger multiplier;
  static BigInteger multiplicand;
  if (argc != 3)
  {
    printf("mult {multiplicand} {multiplier}\n");
    return 1;
  }
  int rc = ComputeExpression(argv[1], &multiplicand);
  if (rc != EXPR_OK)
  {
    printf("Error in the multiplicand: %s\n", argv[1]);
    return 1;
  }
  rc = ComputeExpression(argv[2], &multiplier);
  if (rc != EXPR_OK)
  {
    printf("Error in the multiplier: %s\n", argv[2]);
    return 1;
  }
  BigIntMultiply(&multiplicand, &multiplier, &multiplier);
  char* ptrResults = results;
  Bin2Dec(&ptrResults, multiplier.limbs, multiplier.nbrLimbs, 0);
  *ptrResults = 0;     // Add string terminator.
  (void)printf("%s\n", results);
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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bignbr.h"
#include "highlevel.h"
#include "factor.h"
#define DEBUG_CODE  9
void dilogText(char *baseText, char *powerText, char *modText, int groupLen);
void gaussianText(char *valueText, int doFactorization);
void ecmFrontText(char *tofactorText, int doFactorization, char *knownFactors);
int Factor1[] = { 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x00, 0x00, 0x00 };
int Factor2[] = { 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x00, 0x00, 0x00 };
int Factor3[] = { 29504, 29490, 19798, 633, 181, 0, 0, 0, 0, 0 };
int Factor4[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
int Factor5[] = { 32767, 32767, 32767, 32767, 0, 0, 0, 0 };
int Factor6[] = { 32767, 32767, 32767, 32767, 0, 0, 0, 0 };
int factors[5000];
struct sFactors astFactors[1000];
extern int number[MAX_LEN];
extern int nbrLimbs;
extern int lang, groupLen;
extern limb TestNbr[MAX_LEN];
char expr[] = "123456789012345";
extern char batch;
int Product[32];
char input[10000];
extern char tofactorDec[30000];
BigInteger dividend, divisor, quotient;
int main(int argc, char *argv[])
{
  int len, i;
#if DEBUG_CODE == 1
  fsquaresText(argv[1], 6);
  printf("%s\n", output);
#elif DEBUG_CODE == 2
  fcubesText(argv[1], 6);
  printf("%s\n", output);
#elif DEBUG_CODE == 3
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
  dividend.limbs[0] = 0;
  dividend.limbs[1] = 0;
  dividend.limbs[2] = 0;
  dividend.limbs[3] = 0x20;
  dividend.nbrLimbs = 4;
  dividend.sign = SIGN_POSITIVE;
  divisor.limbs[0] = 0;
  divisor.limbs[1] = 0x05;
  divisor.nbrLimbs = 2;
  dividend.sign = SIGN_NEGATIVE;
  BigIntDivide(&dividend, &divisor, &quotient);
#elif DEBUG_CODE == 5
  for (i = sizeof(expr) - 1; i >= 0; i -= 3)
  {
    expr[i] = 0;
    fcubesText(expr, 6);
    printf("%s\n\n", output);
  }
#elif DEBUG_CODE == 6
  char *ptrInput;
  if (argc != 4)
  {
    printf("num delta den\n");
    return 0;
  }
  ptrInput = input;
  strcpy(ptrInput, argv[1]);
  ptrInput += strlen(ptrInput) + 1;
  strcpy(ptrInput, argv[2]);
  ptrInput += strlen(ptrInput) + 1;
  strcpy(ptrInput, argv[3]);
  ptrInput += strlen(ptrInput) + 1;
  contfracText(input, 20000);
  printf("%s\n", output);
#elif DEBUG_CODE == 7
  char *ptrInput;
  if (argc != 2)
  {
    printf("num\n");
    return 0;
  }
  ptrInput = input;
  strcpy(ptrInput, argv[1]);
  fsquaresText(input, 6);
  printf("%s\n", output);
#elif DEBUG_CODE == 8
  char *ptrInput;
  enum eExprErr rc;
  int NumberLength;
  BigInteger num, mod, inv;

  if (argc != 3)
  {
    printf("num mod\n");
    return 0;
  }
  rc = ComputeExpression(argv[1], 1, &num);
  if (rc != EXPR_OK)
  {
    textError(output, rc);
  }
  else
  {
    rc = ComputeExpression(argv[2], 1, &mod);
    if (rc != EXPR_OK)
    {
      textError(output, rc);
    }
    else
    {
      NumberLength = mod.nbrLimbs;
      if (num.nbrLimbs < mod.nbrLimbs)
      {
        memset(&num.limbs[num.nbrLimbs], 0, (mod.nbrLimbs - num.nbrLimbs) * sizeof(limb));
      }
      memcpy(TestNbr, mod.limbs, NumberLength * sizeof(limb));
      TestNbr[NumberLength].x = 0;
      GetMontgomeryParms(NumberLength);
      ModInvBigNbr(num.limbs, inv.limbs, mod.limbs, NumberLength);
      Bin2Dec(inv.limbs, output, NumberLength, 200);
    }
  }
  printf("%s", output);
#elif DEBUG_CODE == 9
  if (argc != 3)
  {
    printf("modulus polynomial\n");
    return 0;
  }
  polyFactText(argv[1], argv[2], 6);
  printf("%s\n", output);
#elif DEBUG_CODE == 11
  if (argc != 4)
  {
    printf("base power modulus\n");
    return 0;
  }
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
  batch = 0;
  lang = 0;
  if (argc == 3)
  {
    ecmFrontText(argv[1], 1, argv[2]);
    printf("%s\n", output);
  }
  else if (argc == 2)
  {
    char *ptrKnownFactors = strchr(argv[1], '=');
    if (ptrKnownFactors)
    {                          // There is equal sign.
      *ptrKnownFactors = 0;    // Replace equal sign by string terminator.
      ptrKnownFactors++;
    }
    ecmFrontText(argv[1], 1, ptrKnownFactors);
    printf("%s\n", output);
  }
  else
  {
    printf("value [known factors]\n");
    return 0;
  }
#elif DEBUG_CODE == 14
  {
    limb internalNotation[100];
    static int bitGroups;
    char textInput[500];
    textInput[0] = '1';
    memset(&textInput[1], '0', 150);
    textInput[151] = 0;
    Dec2Bin(textInput, internalNotation, 151, &bitGroups);
    Bin2Dec(internalNotation, textInput, 17, 6);
    textInput[200] = 0;
  }
#elif DEBUG_CODE == 15
  memcpy(inputString, "6,-2,00102^1042+1""\0""2^1042+1=5^1(0)*16673^1(0)*627186185377^1(16673)*131294792925870751515684960383613518415615538737991528767912593379854404518341858118366491474959205710499826133822402120149306175263402700301^1(16673)*6864797660130609714981900799081393217269435300143305409394463459185543183397652346775704046543201000705776033378429553397612687501667381169885775070966579201^1(2)""\0\0""222""\0", 10001314 - 10000928 + 1);
  doWork();
#elif DEBUG_CODE == 16
  quadmodText("1", "1", "0", "112856782", 6);
//  quadmodText("1", "1", "0", "56428391", 6);
//  quadmodText("1", "1", "0", "2", 6);
  printf("%s\n", output);
#endif
  return 0;
}

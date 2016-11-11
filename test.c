#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bignbr.h"
#include "highlevel.h"
#include "factor.h"
#define DEBUG_CODE 13
void dilogText(char *baseText, char *powerText, char *modText, int groupLen);
void gaussianText(char *valueText, int doFactorization);
void ecmFrontText(char *tofactorText, int doFactorization);
int Factor1[] = { 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x00, 0x00, 0x00 };
int Factor2[] = { 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x00, 0x00, 0x00 };
int Factor3[] = { 29504, 29490, 19798, 633, 181, 0, 0, 0, 0, 0 };
int Factor4[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
int Factor5[] = { 32767, 32767, 32767, 32767, 0, 0, 0, 0 };
int Factor6[] = { 32767, 32767, 32767, 32767, 0, 0, 0, 0 };
int factors[5000];
struct sFactors astFactors[1000];
extern int karatCtr, multCtr;
extern int number[MAX_LEN];
extern int nbrLimbs;
extern int lang, groupLen;
extern limb TestNbr[MAX_LEN];
char expr[] = "123456789012345";
extern char *output;
int Product[32];
char input[10000];
BigInteger dividend, divisor, quotient;
int main(int argc, char *argv[])
{
  int len, i;
  output = (char *)malloc(10000000);
#if DEBUG_CODE == 1
  fsquaresText(argv[1], 6);
  printf("%s\n", output);
  printf("multiplication count: %d, Karatsuba count: %d", multCtr, karatCtr);
#elif DEBUG_CODE == 2
  fcubesText(argv[1], 6);
  printf("%s\n", output);
  printf("multiplication count: %d, Karatsuba count: %d", multCtr, karatCtr);
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
#elif DEBUG_CODE == 10
  int index;
  limb number[5000];
  BigInteger value;
  if (argc != 2)
  {
    printf("number_to_factor\n");
    return 0;
  }
  Dec2Bin(argv[1], &number[1], (int)strlen(argv[1]), &number[0]);
  factor(number, factors, astFactors);
  for (index = 1; index <= astFactors[0].multiplicity; index++)
  {
    UncompressBigInteger(astFactors[index].ptrFactor, &value);
    Bin2Dec(value.limbs, output, value.nbrLimbs, 6);
    printf("%s^%d\n", output, astFactors[index].multiplicity);
  }
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
  if (argc == 3)
  {
    ecmFrontText(argv[1], 1);
    printf("%s\n", output);
    ecmFrontText(argv[2], 1);
    printf("%s\n", output);
  }
  else if (argc == 2)
  {
    ecmFrontText(argv[1], 1);
    printf("%s\n", output);
  }
  else
  {
    printf("value [value]\n");
    return 0;
  }
#endif
  return 0;
}

#include "output.h"
#include "polynomial.h"

#define TYPE_PM_SPACE_BEFORE 1
#define TYPE_PM_SPACE_AFTER  2

void QuinticEquation(int* ptrPolynomial, int multiplicity);
void showX(int multiplicity);
void showPlusSignOn(int condPlus, int type);
void startSqrt(void);
void endSqrt(void);
void showRatConstants(char* numerator, char* denominator);
void getRootsPolynomial(char **pptrOutput, struct sFactorInfo* pstFactorInfo, int groupLength);

extern BigInteger Quintic, Quartic, Cubic, Quadratic, Linear, Independent;
extern BigInteger discr, commonDenom, tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
extern BigRational RatQuartic, RatCubic, RatQuadratic, RatLinear, RatIndependent;
extern BigRational RatDeprCubic, RatDeprQuadratic, RatDeprLinear, RatDeprIndependent;
extern BigRational RatDiscr, RatDelta0, RatDelta1, RatD;
extern BigRational Rat1, Rat2, Rat3, Rat4, Rat5, RatS;
extern int indexRoot;
extern char *ptrMinus, *ptrTimes;
extern int polyInteger[1000000];

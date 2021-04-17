#include "output.h"
#include "polynomial.h"

#define TYPE_PM_SPACE_BEFORE 1
#define TYPE_PM_SPACE_AFTER  2

void QuinticEquation(int* ptrPolynomial, int multiplicity);
void startLine(void);
void endLine(void);
void showX(int multiplicity);
void showPlusSignOn(int condPlus, int type);
void startParen(void);
void endParen(void);
void startSqrt(void);
void endSqrt(void);
void showRatConstants(char* numerator, char* denominator);
void getRootsPolynomial(int nbrFactor, char **pptrOutput, struct sFactorInfo* pstFactorInfo, int groupLength);

extern BigInteger Quintic;
extern BigInteger Quartic;
extern BigInteger Cubic;
extern BigInteger Quadratic;
extern BigInteger Linear;
extern BigInteger Independent;
extern BigInteger discr;
extern BigInteger commonDenom;
extern BigInteger tmp0;
extern BigInteger tmp1;
extern BigInteger tmp2;
extern BigInteger tmp3;
extern BigInteger tmp4;
extern BigInteger tmp5;
extern BigInteger tmp6;
extern BigInteger tmp7;
extern BigRational RatQuartic;
extern BigRational RatCubic;
extern BigRational RatQuadratic;
extern BigRational RatLinear;
extern BigRational RatIndependent;
extern BigRational RatDeprCubic;
extern BigRational RatDeprQuadratic;
extern BigRational RatDeprLinear;
extern BigRational RatDeprIndependent;
extern BigRational RatDiscr;
extern BigRational RatDelta0;
extern BigRational RatDelta1;
extern BigRational RatD;
extern BigRational Rat1;
extern BigRational Rat2;
extern BigRational Rat3;
extern BigRational Rat4;
extern BigRational Rat5;
extern BigRational RatS;
extern int indexRoot;
extern char *ptrMinus;
extern char *ptrTimes;
extern int polyInteger[1000000];

/* Code from MUSL library */

#include <stdint.h>
#include <string.h>
union fi { double f; uint64_t i; };
#define asuint64(f) (((union{double _f; uint64_t _i;}){f})._i)
#define asdouble(i) (((union{uint64_t _i; double _f;}){i})._f)
#define GET_HIGH_WORD(hi,d)                       \
do {                                              \
  (hi) = asuint64(d) >> 32;                       \
} while (0)

#define GET_LOW_WORD(lo,d)                        \
do {                                              \
  (lo) = (uint32_t)asuint64(d);                   \
} while (0)

#define INSERT_WORDS(d,hi,lo)                     \
do {                                              \
  (d) = asdouble(((uint64_t)(hi)<<32) | (uint32_t)(lo)); \
} while (0)

#define SET_LOW_WORD(d,lo)                        \
  INSERT_WORDS(d, asuint64(d)>>32, lo)

static inline void fp_force_evalf(float x)
{
  volatile float y;
  y = x;
}

static inline void fp_force_eval(double x)
{
  volatile double y;
  y = x;
}

static inline void fp_force_evall(long double x)
{
  volatile long double y;
  y = x;
}

#define FORCE_EVAL(x) do {                        \
	if (sizeof(x) == sizeof(float)) {         \
		fp_force_evalf(x);                \
	} else if (sizeof(x) == sizeof(double)) { \
		fp_force_eval(x);                 \
	} else {                                  \
		fp_force_evall(x);                \
	}                                         \
} while(0)

#define EPS 2.2204460492503131e-16
static const double toint = 1.5 / EPS;
static const double pio4 = 0x1.921fb54442d18p-1;
static const double invpio2 = 6.36619772367581382433e-01; /* 0x3FE45F30, 0x6DC9C883 */
static const double pio2_1 = 1.57079632673412561417e+00; /* 0x3FF921FB, 0x54400000 */
static const double pio2_1t = 6.07710050650619224932e-11; /* 0x3DD0B461, 0x1A626331 */
static const double pio2_2 = 6.07710050630396597660e-11; /* 0x3DD0B461, 0x1A600000 */
static const double pio2_2t = 2.02226624879595063154e-21; /* 0x3BA3198A, 0x2E037073 */
static const double pio2_3 = 2.02226624871116645580e-21; /* 0x3BA3198A, 0x2E000000 */
static const double pio2_3t = 8.47842766036889956997e-32; /* 0x397B839A, 0x252049C1 */

/* caller must handle the case when reduction is not needed: |x| ~<= pi/4 */
int __rem_pio2(double x, double* y)
{
  union { double f; uint64_t i; } u = { x };
  double z;
  int sign;

  sign = u.i >> 63;
  if ((x <= 3.926990816987241548078) && (x >= -3.926990816987241548078))
  {  /* |x| ~<= 5pi/4 */
    if ((x <= 2.356194490192344928846) && (x >= -2.356194490192344928846))
    {  /* |x| ~<= 3pi/4 */
      if (!sign) {
        z = x - pio2_1;  /* one round good to 85 bits */
        y[0] = z - pio2_1t;
        y[1] = (z - y[0]) - pio2_1t;
        return 1;
      }
      else {
        z = x + pio2_1;
        y[0] = z + pio2_1t;
        y[1] = (z - y[0]) + pio2_1t;
        return -1;
      }
    }
    else {
      if (!sign) {
        z = x - (2.0 * pio2_1);
        y[0] = z - (2.0 * pio2_1t);
        y[1] = (z - y[0]) - (2.0 * pio2_1t);
        return 2;
      }
      else {
        z = x + (2.0 * pio2_1);
        y[0] = z + (2.0 * pio2_1t);
        y[1] = (z - y[0]) + (2.0 * pio2_1t);
        return -2;
      }
    }
  }
  if ((x <= 7.068583470577034786540) && (x >= -7.068583470577034786540))
  {  /* |x| ~<= 9pi/4 */
    if ((x <= 5.497787143782138167309) && (x >= -5.497787143782138167309))
    {  /* |x| ~<= 7pi/4 */
      if (!sign)
      {
        z = x - (3.0 * pio2_1);
        y[0] = z - (3.0 * pio2_1t);
        y[1] = (z - y[0]) - (3.0 * pio2_1t);
        return 3;
      }
      else
      {
        z = x + (3.0 * pio2_1);
        y[0] = z + (3.0 * pio2_1t);
        y[1] = (z - y[0]) + (3.0 * pio2_1t);
        return -3;
      }
    }
    else
    {
      if (!sign)
      {
        z = x - (4.0 * pio2_1);
        y[0] = z - (4.0 * pio2_1t);
        y[1] = (z - y[0]) - (4.0 * pio2_1t);
        return 4;
      }
      else
      {
        z = x + (4.0 * pio2_1);
        y[0] = z + (4.0 * pio2_1t);
        y[1] = (z - y[0]) + (4.0 * pio2_1t);
        return -4;
      }
    }
  }
  y[0] = 0;
  y[1] = 0;
  return 0;   // It should never come here.
}

int memcmp(const void *vl, const void *vr, size_t count)
{
  size_t n = count;
  const unsigned char *l=(const unsigned char *)vl;
  const unsigned char *r=(const unsigned char *)vr;
  for (; (n != 0U) && (*l == *r); n--)
  {
    l++;
    r++;
  }
  if (n != 0U)
  {
    char diff = *l - *r;
    int rc = (int)diff;
    return rc;
  }
  return 0;
}

void *memcpy(void *dest, const void *src, size_t count)
{
  size_t n = count;
  unsigned char *d = (unsigned char *)dest;
  const unsigned char *s = (const unsigned char *)src;

#define LS >>
#define RS <<

  typedef uint32_t u32;
  if (n == 0U)
  {
    return dest;
  }
  for (; (((uintptr_t)s % 4U) != 0U) && (n != 0U); n--)
  {
    *d = *s;
    d++;
    s++;
  }

  if (((uintptr_t)d % 4U) == 0U)
  {
    for (; n>=16U; n-=16U)
    {
      *(u32 *)(d+0) = *(const u32 *)(s+0);
      *(u32 *)(d+4) = *(const u32 *)(s+4);
      *(u32 *)(d+8) = *(const u32 *)(s+8);
      *(u32 *)(d+12) = *(const u32 *)(s+12);
      s+=16;
      d+=16;
    }
    if ((n & 8U) != 0U)
    {
      *(u32 *)(d+0) = *(const u32 *)(s+0);
      *(u32 *)(d+4) = *(const u32 *)(s+4);
      d += 8;
      s += 8;
    }
    if ((n & 4U) != 0U)
    {
      *(u32 *)(d+0) = *(const u32 *)(s+0);
      d += 4;
      s += 4;
    }
    if ((n & 2U) != 0U)
    {
      *d = *s;
      *(d + 1) = *(s + 1);
      d += 2;
      s += 2;
    }
    if ((n & 1U) != 0U)
    {
      *d = *s;
    }
    return dest;
  }

  if (n >= 32U)
  {
    uint32_t x;
    uint32_t w;
    switch ((uintptr_t)d % 4U)
    {
      case 1:
        w = *(const u32 *)s;
        *d = *s;
        *(d + 1) = *(s + 1);
        *(d + 2) = *(s + 2);
        d += 3;
        s += 3;
        n -= 3U;
        for (; n>=17U; n-=16U)
        {
          const u32* ptrSrc = (const u32*)(s + 1);
          u32* ptrDest = (u32*)d;
          x = *ptrSrc;
          *ptrDest = (w LS 24) | (x RS 8);
          w = *(ptrSrc + 1);
          *(ptrDest + 1) = (x LS 24) | (w RS 8);
          x = *(ptrSrc + 2);
          *(ptrDest + 2) = (w LS 24) | (x RS 8);
          w = *(ptrSrc + 3);
          *(ptrDest + 3) = (x LS 24) | (w RS 8);
          s += 16;
          d += 16;
        }
        break;
      case 2:
        w = *(const u32 *)s;
        *d = *s;
        *(d + 1) = *(s + 1);
        d += 2;
        s += 2;
        n -= 2U;
        for (; n>=18U; n-=16U)
        {
          const u32* ptrSrc = (const u32*)(s + 2);
          u32* ptrDest = (u32*)d;
          x = *ptrSrc;
          *ptrDest = (w LS 16) | (x RS 16);
          w = *(ptrSrc + 1);
          *(ptrDest + 1) = (x LS 16) | (w RS 16);
          x = *(ptrSrc + 2);
          *(ptrDest + 2) = (w LS 16) | (x RS 16);
          w = *(ptrSrc + 3);
          *(ptrDest + 3) = (x LS 16) | (w RS 16);
          s += 16;
          d += 16;
        }
        break;
      case 3:
        w = *(const u32 *)s;
        *d = *s;
        d++;
        s++;
        n -= 1U;
        for (; n>=19U; n-=16U)
        {
          const u32* ptrSrc = (const u32*)(s + 3);
          u32* ptrDest = (u32*)d;
          x = *ptrSrc;
          *ptrDest = (w LS 8) | (x RS 24);
          w = *(ptrSrc + 1);
          *(ptrDest + 1) = (x LS 8) | (w RS 24);
          x = *(ptrSrc + 2);
          *(ptrDest + 2) = (w LS 8) | (x RS 24);
          w = *(ptrSrc + 3);
          *(ptrDest + 3) = (x LS 8) | (w RS 24);
          s += 16;
          d += 16;
        }
        break;
      default:
        break;
    }
  }
  if ((n & 16U) != 0U)
  {
    *d = *s;
    *(d + 1) = *(s + 1);
    *(d + 2) = *(s + 2);
    *(d + 3) = *(s + 3);
    *(d + 4) = *(s + 4);
    *(d + 5) = *(s + 5);
    *(d + 6) = *(s + 6);
    *(d + 7) = *(s + 7);
    *(d + 8) = *(s + 8);
    *(d + 9) = *(s + 9);
    *(d + 10) = *(s + 10);
    *(d + 11) = *(s + 11);
    *(d + 12) = *(s + 12);
    *(d + 13) = *(s + 13);
    *(d + 14) = *(s + 14);
    *(d + 15) = *(s + 15);
    d += 16;
    s += 16;
  }
  if ((n & 8U) != 0U)
  {
    *d = *s;
    *(d + 1) = *(s + 1);
    *(d + 2) = *(s + 2);
    *(d + 3) = *(s + 3);
    *(d + 4) = *(s + 4);
    *(d + 5) = *(s + 5);
    *(d + 6) = *(s + 6);
    *(d + 7) = *(s + 7);
    d += 8;
    s += 8;
  }
  if ((n & 4U) != 0U)
  {
    *d = *s;
    *(d + 1) = *(s + 1);
    *(d + 2) = *(s + 2);
    *(d + 3) = *(s + 3);
    d += 4;
    s += 4;
  }
  if ((n & 2U) != 0U)
  {
    *d = *s;
    *(d + 1) = *(s + 1);
    d += 2;
    s += 2;
  }
  if ((n & 1U) != 0U)
  {
    *d = *s;
  }
  return dest;
}
#ifdef __clang__
#define WT size_t
#define WS (sizeof(WT))
void *memmove(void *dest, const void *src, size_t count)
{
  char *d = dest;
  const char *s = src;
  size_t n = count;

  if (d == s)
  {
    return d;
  }
  if (((s + n) <= d) || ((d + n) <= s))
  {
    return memcpy(d, s, n);
  }
  if (d<s)
  {
    if (((uintptr_t)s % WS) == ((uintptr_t)d % WS))
    {
      while (((uintptr_t)d % WS) != 0U)
      {
        if (n == 0U)
        {
          return dest;
        }
        n--;
        *d = *s;
        d++;
        s++;
      }
      for (; n>=WS; n-=WS)
      {
        *(WT *)d = *(const WT *)s;
        d+=WS;
        s+=WS;
      }
    }
    for (; n != 0U; n--)
    {
      *d = *s;
      d++;
      s++;
    }
  }
  else
  {
    if (((uintptr_t)s % WS) == ((uintptr_t)d % WS))
    {
      while (((uintptr_t)(d+n) % WS) != 0U)
      {
        if (n == 0U)
        {
          return dest;
        }
        n--;
        d[n] = s[n];
      }
      while (n>=WS)
      {
        n-=WS;
        *(WT *)(d+n) = *(const WT *)(s+n);
      }
    }
    while (n != 0U)
    {
      n--;
      d[n] = s[n];
    }
  }

  return dest;
}
#endif
void *memset(void *dest, int c, size_t count)
{
  unsigned char *s = (unsigned char *)dest;
  int k;
  int n = (int)count;

  /* Fill head and tail with minimal branching. Each
   * conditional ensures that all the subsequently used
   * offsets are well-defined and in the dest region. */

  if (n == 0)
  {
    return dest;
  }
  s[0] = (unsigned char)c; 
  s[n-1] = (unsigned char)c;
  if (n <= 2)
  {
    return dest;
  }
  s[1] = (unsigned char)c; 
  s[n-2] = (unsigned char)c;
  s[2] = (unsigned char)c; 
  s[n-3] = (unsigned char)c;
  if (n <= 6)
  { 
    return dest; 
  }
  s[3] = (unsigned char)c;
  s[n-4] = (unsigned char)c;
  if (n <= 8)
  {
    return dest;
  }

  /* Advance pointer to align it at a 4-byte boundary,
   * and truncate n to a multiple of 4. The previous code
   * already took care of any head/tail that get cut off
   * by the alignment. */

  k = (0x03 - (int)(uintptr_t)s+1) & 3;
  s += k;
  n -= k;
  n &= 0xFFFFFFFCU;

  typedef uint32_t  u32;
  typedef uint64_t  u64;

  u32 c32 = 0x01010101U * (unsigned char)c;
  int words = n / 4;
  u32* ptrStart = (u32 *)s;
  u32* ptrEnd = ptrStart + words;

  /* In preparation to copy 32 bytes at a time, aligned on
   * an 8-byte bounary, fill head/tail up to 28 bytes each.
   * As in the initial byte-based head/tail fill, each
   * conditional below ensures that the subsequent offsets
   * are valid (e.g. !(n<=24) implies n>=28). */

  *ptrStart = c32;
  *(ptrEnd - 1) = c32;
  if (n <= 8)
  {
    return dest;
  }
  *(ptrStart + 1) = c32;
  *(ptrStart + 2) = c32;
  *(ptrEnd - 3) = c32;
  *(ptrEnd - 2) = c32;
  if (n <= 24)
  {
    return dest;
  }
  *(ptrStart + 3) = c32;
  *(ptrStart + 4) = c32;
  *(ptrStart + 5) = c32;
  *(ptrStart + 6) = c32;
  *(ptrEnd - 7) = c32;
  *(ptrEnd - 6) = c32;
  *(ptrEnd - 5) = c32;
  *(ptrEnd - 4) = c32;

  /* Align to a multiple of 8 so we can fill 64 bits at a time,
   * and avoid writing the same bytes twice as much as is
   * practical without introducing additional branching. */

  k = 24 + ((int)(uintptr_t)s & 4);
  s += k;
  n -= k;

  /* If this loop is reached, 28 tail bytes have already been
   * filled, so any remainder when n drops below 32 can be
   * safely ignored. */

  u64 c64 = c32 | ((u64)c32 << 32);
  for (; n >= 32; n-=32)
  {
    u64* ptrStart64 = (u64*)s;
    *ptrStart64 = c64;
    *(ptrStart64 + 1) = c64;
    *(ptrStart64 + 2) = c64;
    *(ptrStart64 + 3) = c64;
    s += 32;
  }
  return dest;
}

char *strcpy(char *dest, const char *src)
{
  char *ptrDest = dest;
  const char* ptrSrc = src;
  while (*ptrSrc != '\0')
  {
    *ptrDest = *ptrSrc;
    ptrDest++;
    ptrSrc++;
  }
  *ptrDest = '\0';
  return dest;
}

size_t strlen(const char *s)
{
  const char *a = s;
  while (*a != 0)
  {
    a++;
  }
  return a - s;
}

double scalbn(double x, int n)
{
  union fi u;
  double y = x;
  int value = 0x3ff + n;

  u.i = (uint64_t)value << 52;
  return y * u.f;
}

/*
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunSoft, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */

static const double ln2_hi = 6.93147180369123816490e-01;  /* 3fe62e42 fee00000 */
static const double ln2_lo = 1.90821492927058770002e-10;  /* 3dea39ef 35793c76 */
static const double Lg1 = 6.666666666666735130e-01;  /* 3FE55555 55555593 */
static const double Lg2 = 3.999999999940941908e-01;  /* 3FD99999 9997FA04 */
static const double Lg3 = 2.857142874366239149e-01;  /* 3FD24924 94229359 */
static const double Lg4 = 2.222219843214978396e-01;  /* 3FCC71C5 1D8E78AF */
static const double Lg5 = 1.818357216161805012e-01;  /* 3FC74664 96CB03DE */
static const double Lg6 = 1.531383769920937332e-01;  /* 3FC39A09 D078C69F */
static const double Lg7 = 1.479819860511658591e-01;  /* 3FC2F112 DF3E5244 */

double log(double argum)
{
  double x = argum;
  union fi u;
  u.f = x;
  double hfsq;
  double f;
  double s;
  double z;
  double w;
  double t1;
  double t2;
  double dk;
  uint32_t hx;
  int k;

  hx = u.i>>32;
  k = 0;

  /* reduce x into [sqrt(2)/2, sqrt(2)] */
  hx += 0x3ff00000U - 0x3fe6a09eU;
  k += (int)(hx>>20) - 0x3ff;
  hx = (hx&0x000fffffU) + 0x3fe6a09eU;
  u.i = ((uint64_t)hx<<32) | (u.i&0xffffffff);
  x = u.f;

  f = x - 1.0;
  hfsq = 0.5*f*f;
  s = f/(2.0+f);
  z = s*s;
  w = z*z;
  t1 = w*(Lg2+(w*(Lg4+(w*Lg6))));
  t2 = z*(Lg1+(w*(Lg3+(w*(Lg5+(w*Lg7))))));
  dk = k;
  return (s*(hfsq+ t2 + t1)) + (dk*ln2_lo) - hfsq + f + (dk*ln2_hi);
}

static const double half[2] = {0.5,-0.5};
static const double ln2hi = 6.93147180369123816490e-01; /* 0x3fe62e42, 0xfee00000 */
static const double ln2lo = 1.90821492927058770002e-10; /* 0x3dea39ef, 0x35793c76 */
static const double invln2 = 1.44269504088896338700e+00; /* 0x3ff71547, 0x652b82fe */
static const double P1   =  1.66666666666666019037e-01; /* 0x3FC55555, 0x5555553E */
static const double P2   = -2.77777777770155933842e-03; /* 0xBF66C16C, 0x16BEBD93 */
static const double P3   =  6.61375632143793436117e-05; /* 0x3F11566A, 0xAF25DE2C */
static const double P4   = -1.65339022054652515390e-06; /* 0xBEBBBD41, 0xC5D26BF1 */
static const double P5   =  4.13813679705723846039e-08; /* 0x3E663769, 0x72BEA4D0 */

double exp(double argum)
{
  double x = argum;
  double hi;
  double lo;
  double c;
  double xx;
  double y;
  int k;
  int sign;
  uint32_t hx;
  union fi __u;
  __u.f = x;
  hx = __u.i >> 32;
  sign = hx>>31;
  hx &= 0x7fffffff;  /* high word of |x| */

  /* argument reduction */
  if (hx > 0x3fd62e42U)
  {  /* if |x| > 0.5 ln2 */
    if (hx >= 0x3ff0a2b2U)
    {/* if |x| >= 1.5 ln2 */
      k = (int)((invln2 * x) + half[sign]);
    }
    else
    {
      k = 1 - sign - sign;
    }
    hi = x - ((double)k*ln2hi);  /* k*ln2hi is exact here */
    lo = (double)k*ln2lo;
    x = hi - lo;
  }
  else if (hx > 0x3e300000U)
  {  /* if |x| > 2**-28 */
    k = 0;
    hi = x;
    lo = 0;
  }
  else
  {
    /* inexact if x!=0 */
    return 1.0 + x;
  }

  /* x is now in primary range */
  xx = x*x;
  c = x - (xx * (P1 + (xx * (P2 + (xx * (P3 + (xx * (P4 + (xx * P5)))))))));
  y = 1.0 + (x*c/(2.0-c)) - lo + hi;
  if (k == 0)
  {
    return y;
  }
  return scalbn(y, k);
}

double pow(double x, double y)
{
  return exp(y*log(x));
}

// Square root.
static const uint16_t __rsqrt_tab[128] = 
{
  0xb451,0xb2f0,0xb196,0xb044,0xaef9,0xadb6,0xac79,0xab43,
  0xaa14,0xa8eb,0xa7c8,0xa6aa,0xa592,0xa480,0xa373,0xa26b,
  0xa168,0xa06a,0x9f70,0x9e7b,0x9d8a,0x9c9d,0x9bb5,0x9ad1,
  0x99f0,0x9913,0x983a,0x9765,0x9693,0x95c4,0x94f8,0x9430,
  0x936b,0x92a9,0x91ea,0x912e,0x9075,0x8fbe,0x8f0a,0x8e59,
  0x8daa,0x8cfe,0x8c54,0x8bac,0x8b07,0x8a64,0x89c4,0x8925,
  0x8889,0x87ee,0x8756,0x86c0,0x862b,0x8599,0x8508,0x8479,
  0x83ec,0x8361,0x82d8,0x8250,0x81c9,0x8145,0x80c2,0x8040,
  0xff02,0xfd0e,0xfb25,0xf947,0xf773,0xf5aa,0xf3ea,0xf234,
  0xf087,0xeee3,0xed47,0xebb3,0xea27,0xe8a3,0xe727,0xe5b2,
  0xe443,0xe2dc,0xe17a,0xe020,0xdecb,0xdd7d,0xdc34,0xdaf1,
  0xd9b3,0xd87b,0xd748,0xd61a,0xd4f1,0xd3cd,0xd2ad,0xd192,
  0xd07b,0xcf69,0xce5b,0xcd51,0xcc4a,0xcb48,0xca4a,0xc94f,
  0xc858,0xc764,0xc674,0xc587,0xc49d,0xc3b7,0xc2d4,0xc1f4,
  0xc116,0xc03c,0xbf65,0xbe90,0xbdbe,0xbcef,0xbc23,0xbb59,
  0xba91,0xb9cc,0xb90a,0xb84a,0xb78c,0xb6d0,0xb617,0xb560,
};

/* returns a*b*2^-32 - e, with error 0 <= e < 1.  */
static inline uint32_t mul32(uint32_t a, uint32_t b)
{
  return ((uint64_t)a * b) >> 32;
}

/* returns a*b*2^-64 - e, with error 0 <= e < 3.  */
static inline uint64_t mul64(uint64_t a, uint64_t b)
{
  uint64_t ahi = a >> 32;
  uint64_t alo = a & 0xffffffff;
  uint64_t bhi = b >> 32;
  uint64_t blo = b & 0xffffffff;
  return (ahi * bhi) + ((ahi * blo) >> 32) + ((alo * bhi) >> 32);
}

double sqrt(double x)
{
  uint64_t ix;
  uint64_t top;
  uint64_t m;

  /* special case handling.  */
  ix = asuint64(x);
  top = ix >> 52;
  if (x == 0.0)
  {
    return 0;
  }

  int even = (int)top & 1;
  m = (ix << 11) | 0x8000000000000000ULL;
  if (even != 0)
  {
    m >>= 1;
  }
  top = (top + 0x3ffULL) >> 1;

  static const uint64_t three = 0xc0000000ULL;
  uint64_t r;
  uint64_t s;
  uint64_t d;
  uint64_t u;
  uint64_t i;

  i = (ix >> 46) % 128;
  r = (uint32_t)__rsqrt_tab[i] << 16;
  /* |r sqrt(m) - 1| < 0x1.fdp-9 */
  s = mul32(m >> 32, r);
  /* |s/sqrt(m) - 1| < 0x1.fdp-9 */
  d = mul32(s, r);
  u = three - d;
  r = mul32(r, u) << 1;
  /* |r sqrt(m) - 1| < 0x1.7bp-16 */
  s = mul32(s, u) << 1;
  /* |s/sqrt(m) - 1| < 0x1.7bp-16 */
  d = mul32(s, r);
  u = three - d;
  r = mul32(r, u) << 1;
  /* |r sqrt(m) - 1| < 0x1.3704p-29 (measured worst-case) */
  r = r << 32;
  s = mul64(m, r);
  d = mul64(s, r);
  u = (three << 32) - d;
  s = mul64(s, u);  /* repr: 3.61 */
  /* -0x1p-57 < s - sqrt(m) < 0x1.8001p-61 */
  s = (s - 2) >> 9; /* repr: 12.52 */
  /* -0x1.09p-52 < s - sqrt(m) < -0x1.fffcp-63 */

  /* s < sqrt(m) < s + 0x1.09p-52,
     compute nearest rounded result:
     the nearest result to 52 bits is either s or s+0x1p-52,
     we can decide by comparing (2^52 s + 0.5)^2 to 2^104 m.  */
  uint64_t d0;
  uint64_t d1;
  double y;
  d0 = (m << 42) - (s * s);
  d1 = s - d0;
  s += d1 >> 63;
  s &= 0x000fffffffffffff;
  s |= top << 52;
  y = asdouble(s);
  return y;
}

// Trigonometric functions

static const double S1 = -1.66666666666666324348e-01; /* 0xBFC55555, 0x55555549 */
static const double S2 = 8.33333333332248946124e-03; /* 0x3F811111, 0x1110F8A6 */
static const double S3 = -1.98412698298579493134e-04; /* 0xBF2A01A0, 0x19C161D5 */
static const double S4 = 2.75573137070700676789e-06; /* 0x3EC71DE3, 0x57B1FE7D */
static const double S5 = -2.50507602534068634195e-08; /* 0xBE5AE5E6, 0x8A2B9CEB */
static const double S6 = 1.58969099521155010221e-10; /* 0x3DE5D93A, 0x5ACFD57C */

double __sin(double x, double y, int iy)
{
  double z;
  double r;
  double v;
  double w;

  z = x * x;
  w = z * z;
  r = S2 + (z * (S3 + (z * S4))) + (z * w * (S5 + (z * S6)));
  v = z * x;
  if (iy == 0)
  {
    return x + (v * (S1 + (z * r)));
  }
  return x - ((z * ((0.5 * y) - (v * r)) - y) - (v * S1));
}


static const double C1 = 4.16666666666666019037e-02; /* 0x3FA55555, 0x5555554C */
static const double C2 = -1.38888888888741095749e-03; /* 0xBF56C16C, 0x16C15177 */
static const double C3 = 2.48015872894767294178e-05; /* 0x3EFA01A0, 0x19CB1590 */
static const double C4 = -2.75573143513906633035e-07; /* 0xBE927E4F, 0x809C52AD */
static const double C5 = 2.08757232129817482790e-09; /* 0x3E21EE9E, 0xBDB4B1C4 */
static const double C6 = -1.13596475577881948265e-11; /* 0xBDA8FAE9, 0xBE8838D4 */

double __cos(double x, double y)
{
  double hz;
  double z;
  double r;
  double w;

  z = x * x;
  w = z * z;
  r = z * (C1 + (z * (C2 + (z * C3)))) + (w * w * (C4 + z * (C5 + (z * C6))));
  hz = 0.5 * z;
  w = 1.0 - hz;
  return w + (((1.0 - w) - hz) + ((z * r) - (x * y)));
}

double cos(double x)
{
  double y[2];
  uint32_t ix;
  unsigned n;

  GET_HIGH_WORD(ix, x);
  ix &= 0x7fffffff;

  /* |x| ~< pi/4 */
  if (ix <= 0x3fe921fbU)
  {
    if (ix < 0x3e46a09eU)
    {  /* |x| < 2**-27 * sqrt(2) */
      /* raise inexact if x!=0 */
      FORCE_EVAL(x + 0x1p120f);
      return 1.0;
    }
    return __cos(x, 0);
  }

  /* cos(Inf or NaN) is NaN */
  if (ix >= 0x7ff00000U)
  {
    return asdouble(0x7FF8000000000000LL);
  }
  /* argument reduction */
  n = __rem_pio2(x, y);
  switch (n & 3U) {
  case 0: return  __cos(y[0], y[1]);
  case 1: return -__sin(y[0], y[1], 1);
  case 2: return -__cos(y[0], y[1]);
  default:
    return  __sin(y[0], y[1], 1);
  }
}

static const double pio2_hi = 1.57079632679489655800e+00; /* 0x3FF921FB, 0x54442D18 */
static const double pio2_lo = 6.12323399573676603587e-17; /* 0x3C91A626, 0x33145C07 */
static const double pS0 = 1.66666666666666657415e-01; /* 0x3FC55555, 0x55555555 */
static const double pS1 = -3.25565818622400915405e-01; /* 0xBFD4D612, 0x03EB6F7D */
static const double pS2 = 2.01212532134862925881e-01; /* 0x3FC9C155, 0x0E884455 */
static const double pS3 = -4.00555345006794114027e-02; /* 0xBFA48228, 0xB5688F3B */
static const double pS4 = 7.91534994289814532176e-04; /* 0x3F49EFE0, 0x7501B288 */
static const double pS5 = 3.47933107596021167570e-05; /* 0x3F023DE1, 0x0DFDF709 */
static const double qS1 = -2.40339491173441421878e+00; /* 0xC0033A27, 0x1C8A2D4B */
static const double qS2 = 2.02094576023350569471e+00; /* 0x40002AE5, 0x9C598AC8 */
static const double qS3 = -6.88283971605453293030e-01; /* 0xBFE6066C, 0x1B8D0159 */
static const double qS4 = 7.70381505559019352791e-02; /* 0x3FB3B8C5, 0xB12E9282 */

static double R(double z)
{
  double p;
  double q;
  p = z * (pS0 + (z * (pS1 + (z * (pS2 + (z * (pS3 + (z * (pS4 + (z * pS5))))))))));
  q = 1.0 + (z * (qS1 + (z * (qS2 + (z * (qS3 + (z * qS4)))))));
  return p / q;
}

double acos(double x)
{
  double z;
  double w;
  double s;
  double c;
  double df;
  uint32_t hx;
  uint32_t ix;

  GET_HIGH_WORD(hx, x);
  ix = hx & 0x7fffffffU;
  /* |x| >= 1 or nan */
  if (ix >= 0x3ff00000U)
  {
    uint32_t lx;

    GET_LOW_WORD(lx, x);
    if ((ix - 0x3ff00000 | lx) == 0) {
      /* acos(1)=0, acos(-1)=pi */
      if ((hx >> 31) != 0U)
      {
        return (2.0 * pio2_hi) + 0x1p-120f;
      }
      return 0.0;
    }
    return asdouble(0x7FF0000000000000LL);
  }
  /* |x| < 0.5 */
  if (ix < 0x3fe00000U)
  {
    if (ix <= 0x3c600000U)
    {  /* |x| < 2**-57 */
      return pio2_hi + 0x1p-120f;
    }
    return pio2_hi - (x - (pio2_lo - (x * R(x * x))));
  }
  /* x < -0.5 */
  if (hx >> 31)
  {
    z = (1.0 + x) * 0.5;
    s = sqrt(z);
    w = (R(z) * s) - pio2_lo;
    return 2.0 * (pio2_hi - (s + w));
  }
  /* x > 0.5 */
  z = (1.0 - x) * 0.5;
  s = sqrt(z);
  df = s;
  SET_LOW_WORD(df, 0);
  c = (z - (df * df)) / (s + df);
  w = (R(z) * s) + c;
  return 2.0 * (df + w);
}


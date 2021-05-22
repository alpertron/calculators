/* Code from MUSL library */

#include <stdint.h>
#include <string.h>
union fi { double f; uint64_t i; };

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
  return (n != 0U)? (*l-*r) : 0;
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
          x = *(const u32 *)(s+1);
          *(u32 *)(d+0) = (w LS 24) | (x RS 8);
          w = *(const u32 *)(s+5);
          *(u32 *)(d+4) = (x LS 24) | (w RS 8);
          x = *(const u32 *)(s+9);
          *(u32 *)(d+8) = (w LS 24) | (x RS 8);
          w = *(const u32 *)(s+13);
          *(u32 *)(d+12) = (x LS 24) | (w RS 8);
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
          x = *(const u32 *)(s+2);
          *(u32 *)(d+0) = (w LS 16) | (x RS 16);
          w = *(const u32 *)(s+6);
          *(u32 *)(d+4) = (x LS 16) | (w RS 16);
          x = *(const u32 *)(s+10);
          *(u32 *)(d+8) = (w LS 16) | (x RS 16);
          w = *(const u32 *)(s+14);
          *(u32 *)(d+12) = (x LS 16) | (w RS 16);
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
          x = *(const u32 *)(s+3);
          *(u32 *)(d+0) = (w LS 8) | (x RS 24);
          w = *(const u32 *)(s+7);
          *(u32 *)(d+4) = (x LS 8) | (w RS 24);
          x = *(const u32 *)(s+11);
          *(u32 *)(d+8) = (w LS 8) | (x RS 24);
          w = *(const u32 *)(s+15);
          *(u32 *)(d+12) = (x LS 8) | (w RS 24);
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
void *memmove(void *dest, const void *src, size_t n)
{
  char *d = dest;
  const char *s = src;

  if (d == s)
  {
    return d;
  }
  if (((s + n) <= d) || ((d + n) <= s))
  {
    return (void)memcpy(d, s, n);
  }
  if (d<s)
  {
    if (((uintptr_t)s % WS) == ((uintptr_t)d % WS))
    {
      while (((uintptr_t)d % WS) != 0)
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
        *(WT *)d = *(WT *)s;
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
      while (((uintptr_t)(d+n) % WS) != 0)
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
        *(WT *)(d+n) = *(WT *)(s+n);
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
void *memset(void *dest, int c, size_t n)
{
  unsigned char *s = (unsigned char *)dest;
  size_t k;

  /* Fill head and tail with minimal branching. Each
   * conditional ensures that all the subsequently used
   * offsets are well-defined and in the dest region. */

  if (n == 0U)
  {
    return dest;
  }
  s[0] = s[n-1U] = (unsigned char)c;
  if (n <= 2U)
  {
    return dest;
  }
  s[1] = s[n-2U] = (unsigned char)c;
  s[2] = s[n-3U] = (unsigned char)c;
  if (n <= 6U)
  { 
    return dest; 
  }
  s[3] = s[n-4U] = (unsigned char)c;
  if (n <= 8U)
  {
    return dest;
  }

  /* Advance pointer to align it at a 4-byte boundary,
   * and truncate n to a multiple of 4. The previous code
   * already took care of any head/tail that get cut off
   * by the alignment. */

  k = (~(uintptr_t)s+1) & 3U;
  s += k;
  n -= k;
  n &= -4U;

  typedef uint32_t  u32;
  typedef uint64_t  u64;

  u32 c32 = ((u32)-1)/255 * (unsigned char)c;

  /* In preparation to copy 32 bytes at a time, aligned on
   * an 8-byte bounary, fill head/tail up to 28 bytes each.
   * As in the initial byte-based head/tail fill, each
   * conditional below ensures that the subsequent offsets
   * are valid (e.g. !(n<=24) implies n>=28). */

  *(u32 *)(s+0) = c32;
  *(u32 *)(s+n-4) = c32;
  if (n <= 8)
  {
    return dest;
  }
  *(u32 *)(s+4) = c32;
  *(u32 *)(s+8) = c32;
  *(u32 *)(s+n-12) = c32;
  *(u32 *)(s+n-8) = c32;
  if (n <= 24)
  {
    return dest;
  }
  *(u32 *)(s+12) = c32;
  *(u32 *)(s+16) = c32;
  *(u32 *)(s+20) = c32;
  *(u32 *)(s+24) = c32;
  *(u32 *)(s+n-28) = c32;
  *(u32 *)(s+n-24) = c32;
  *(u32 *)(s+n-20) = c32;
  *(u32 *)(s+n-16) = c32;

  /* Align to a multiple of 8 so we can fill 64 bits at a time,
   * and avoid writing the same bytes twice as much as is
   * practical without introducing additional branching. */

  k = 24 + ((uintptr_t)s & 4);
  s += k;
  n -= k;

  /* If this loop is reached, 28 tail bytes have already been
   * filled, so any remainder when n drops below 32 can be
   * safely ignored. */

  u64 c64 = c32 | ((u64)c32 << 32);
  for (; n >= 32U; n-=32U)
  {
    *(u64 *)(s+0) = c64;
    *(u64 *)(s+8) = c64;
    *(u64 *)(s+16) = c64;
    *(u64 *)(s+24) = c64;
    s += 32;
  }
  return dest;
}

char *strcpy(char *dest, const char *src)
{
  char *dst = dest;
  while (*src != '\0')
  {
    *dst = *src;
    dst++;
    src++;
  }
  *dst = '\0';
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

  u.i = (uint64_t)(0x3ff + n) << 52;
  x = y * u.f;
  return x;
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

static const double
ln2_hi = 6.93147180369123816490e-01,  /* 3fe62e42 fee00000 */
ln2_lo = 1.90821492927058770002e-10,  /* 3dea39ef 35793c76 */
Lg1 = 6.666666666666735130e-01,  /* 3FE55555 55555593 */
Lg2 = 3.999999999940941908e-01,  /* 3FD99999 9997FA04 */
Lg3 = 2.857142874366239149e-01,  /* 3FD24924 94229359 */
Lg4 = 2.222219843214978396e-01,  /* 3FCC71C5 1D8E78AF */
Lg5 = 1.818357216161805012e-01,  /* 3FC74664 96CB03DE */
Lg6 = 1.531383769920937332e-01,  /* 3FC39A09 D078C69F */
Lg7 = 1.479819860511658591e-01;  /* 3FC2F112 DF3E5244 */

double log(double x)
{
  union fi u;
  u.f = x;
  double hfsq;
  double f;
  double s;
  double z;
  double R;
  double w;
  double t1;
  double t2;
  double dk;
  uint32_t hx;
  int k;

  hx = u.i>>32;
  k = 0;

  /* reduce x into [sqrt(2)/2, sqrt(2)] */
  hx += 0x3ff00000 - 0x3fe6a09e;
  k += (int)(hx>>20) - 0x3ff;
  hx = (hx&0x000fffff) + 0x3fe6a09e;
  u.i = (uint64_t)hx<<32 | (u.i&0xffffffff);
  x = u.f;

  f = x - 1.0;
  hfsq = 0.5*f*f;
  s = f/(2.0+f);
  z = s*s;
  w = z*z;
  t1 = w*(Lg2+w*(Lg4+w*Lg6));
  t2 = z*(Lg1+w*(Lg3+w*(Lg5+w*Lg7)));
  R = t2 + t1;
  dk = k;
  return s*(hfsq+R) + dk*ln2_lo - hfsq + f + dk*ln2_hi;
}

static const double
half[2] = {0.5,-0.5},
ln2hi = 6.93147180369123816490e-01, /* 0x3fe62e42, 0xfee00000 */
ln2lo = 1.90821492927058770002e-10, /* 0x3dea39ef, 0x35793c76 */
invln2 = 1.44269504088896338700e+00, /* 0x3ff71547, 0x652b82fe */
P1   =  1.66666666666666019037e-01, /* 0x3FC55555, 0x5555553E */
P2   = -2.77777777770155933842e-03, /* 0xBF66C16C, 0x16BEBD93 */
P3   =  6.61375632143793436117e-05, /* 0x3F11566A, 0xAF25DE2C */
P4   = -1.65339022054652515390e-06, /* 0xBEBBBD41, 0xC5D26BF1 */
P5   =  4.13813679705723846039e-08; /* 0x3E663769, 0x72BEA4D0 */

double exp(double x)
{
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
  if (hx > 0x3fd62e42)
  {  /* if |x| > 0.5 ln2 */
    if (hx >= 0x3ff0a2b2)
    {/* if |x| >= 1.5 ln2 */
      k = (int)(invln2 * x + half[sign]);
    }
    else
    {
      k = 1 - sign - sign;
    }
    hi = x - k*ln2hi;  /* k*ln2hi is exact here */
    lo = k*ln2lo;
    x = hi - lo;
  }
  else if (hx > 0x3e300000)
  {  /* if |x| > 2**-28 */
    k = 0;
    hi = x;
    lo = 0;
  }
  else
  {
    /* inexact if x!=0 */
    return 1 + x;
  }

  /* x is now in primary range */
  xx = x*x;
  c = x - xx*(P1+xx*(P2+xx*(P3+xx*(P4+xx*P5))));
  y = 1 + (x*c/(2-c) - lo + hi);
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

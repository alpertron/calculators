//
// This file is part of Alpertron Calculators.
//
// Copyright 2019-2021 Dario Alejandro Alpern
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
#include <stdarg.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include "string/strings.h"
#include "bignbr.h"
#ifdef POLYEXPR
#include "polynomial.h"
#include "rootseq.h"
#endif
extern bool hexadecimal;
void copyStr(char** pptrString, const char* stringToCopy)
{
  char* ptrString = *pptrString;
  const char* ptrStringToCopy = stringToCopy;
  while (*ptrStringToCopy != '\0')
  {
    *ptrString = *ptrStringToCopy;
    ptrString++;
    ptrStringToCopy++;
  }
  *ptrString = '\0';
  *pptrString = ptrString;
}

static int getRangeNumber(const char** pptrOutput)
{
  const char* ptrOutput = *pptrOutput;
  int number = 0;
  while ((*ptrOutput >= '0') && (*ptrOutput <= '9'))
  {
    number = number * 10 + (*ptrOutput - '0');
    ptrOutput++;
  }
  *pptrOutput = ptrOutput;
  return number;
}

// Format a string with variable number of arguments.
// Format string (second argument) can contain format specifiers:
// $nf where:
//  n - argument number
//  f - format specifier:
//      d - integer in decimal format
//      u - unsigned integer in decimal format
//      l - long integer in decimal format
//      e - exponent
//      s - string
//      b - BigInteger number
//      a - array of NumberLength limbs
//      v - variable
//      p - power of variable
//      q - equation number
//      ? - conditional (used for plurals)
//          followed by:
//          n?string?? - if argument is n, string is included, else ignored
//          m-n?string?? - if argument is in range m to n, string is included, else ignored
//          n+?string?? - if argument is >= n, string is included, else ignored

typedef union {
  int i;
  int64_t l;
  const char* s;
  char c;
  const BigInteger* b;
  const limb* a;
} ArgValue;

void formatString(char** pptrString, const char* format, ...)
{
  va_list args;
  va_start(args, format);
  ArgValue argv[10]; // supports up to 9 arguments
  char argType[10];  // store 'd' or 's' per index
  memset(argType, 0, sizeof(argType));

  // First pass: find maximum argument index and collect arguments in array
  int maxIndex = 0;
  const char* ptrFormat = format;
  while (*ptrFormat)
  {
    if (*ptrFormat == '$')
    {
      char c = *(ptrFormat+1);
      assert((c >= '0') && (c <= '9')); // ensure valid format
      int idx = c - '0';
      if (idx <= maxIndex)
      {
        ptrFormat++; // Skip '$'.
        continue;    // Argument already processed
      }
      assert(idx == maxIndex + 1);  // It must be the next argument.
      maxIndex = idx;
      char type = *(ptrFormat + 2);
      if (!argType[idx])
      {
        argType[idx] = type;
        if ((type == 'd') || (type == 'u') || (type == 'e') ||
          (type == '?') || (type == 'q'))
        {
          argv[idx].i = va_arg(args, int);
        }
        else if (type == 'l')
        {
          argv[idx].l = va_arg(args, int64_t);
        }
#ifdef POLYEXPR
        else if (type == 'p')
        {
          argv[idx].s = va_arg(args, const char *);
        }
#endif
        else if (type == 'v')
        {
          argv[idx].c = va_arg(args, int);
        }
        else if (type == 's')
        {
          argv[idx].s = va_arg(args, const char*);
        }
        else if (type == 'b')
        {
          argv[idx].b = va_arg(args, const BigInteger*);
        }
        else if (type == 'a')
        {
          argv[idx].a = va_arg(args, const limb*);
        }
        else
        {
          assert(false); // unsupported type
        }
      }
      ptrFormat += 3;
    }
    else
    {
      ptrFormat++;
    }
  }
  va_end(args);

  // Second pass: build output
  char* ptrString = *pptrString;
  ptrFormat = format;
  while (*ptrFormat)
  {
    if (*ptrFormat == '$')
    {
      int idx = *(ptrFormat+1) - '0';
      assert(idx >= 0 && idx <= maxIndex);
      char type = *(ptrFormat+2);
      if (type == 'd')
      {
        if (argv[idx].i < 0)
        {
          copyStr(&ptrString, "&minus;");
          int2dec(&ptrString, -argv[idx].i);
        }
        else
        {
          int2dec(&ptrString, argv[idx].i);
        }
      }
      else if (type == 'e')
      {
        if (argv[idx].i != 1)
        {
          copyStr(&ptrString, "<sup>");
          int2dec(&ptrString, argv[idx].i);
          copyStr(&ptrString, "</sup>");
        }
      }
      else if (type == 'u')
      {
        int2dec(&ptrString, argv[idx].i);
      }
      else if (type == 'q')
      {
        copyStr(&ptrString, "<span class=\"eq\">(");
        int2dec(&ptrString, argv[idx].i);
        copyStr(&ptrString, ")</span>");
      }
      else if (type == 'l')
      {
        long2dec(&ptrString, argv[idx].l);
      }
      else if (type == 's')
      {
        const char* ptrSrcString = argv[idx].s;
        while (*ptrSrcString != '\0')
        {
          *ptrString = *ptrSrcString;
          ptrString++;
          ptrSrcString++;
        }
      }
      else if (type == 'a')
      {
        if (hexadecimal)
        {
          Bin2Hex(&ptrString, argv[idx].a, NumberLength, groupLen);
        }
        else
        {
          Bin2Dec(&ptrString, argv[idx].a, NumberLength, groupLen);
        }
      }
      else if (type == 'b')
      {
        if (hexadecimal)
        {
          BigInteger2Hex(&ptrString, argv[idx].b, groupLen);
        }
        else
        {
          BigInteger2Dec(&ptrString, argv[idx].b, groupLen);
        }
      }
      else if (type == 'v')
      {
#ifdef POLYEXPR
        showVariable(&ptrString, argv[idx].c);
#else
        copyStr(&ptrString, "<var>");
        *ptrString = argv[idx].c;
        ptrString++;
        copyStr(&ptrString, "</var>");
#endif
      }
#ifdef POLYEXPR
      else if (type == 'p')
      {
        showPowerVar(&ptrString, argv[idx].s[1] - '0', argv[idx].s[0]);
      }
#endif
      else if (type == '?')
      {
        // Test range.
        ptrFormat += 3;  // Point to first character after '?'
        int minimum = getRangeNumber(&ptrFormat);
        int maximum;
        if (*ptrFormat == '-')
        {   // Range separator.
          ptrFormat++;  // Skip range separator.
          maximum = getRangeNumber(&ptrFormat);
        }
        else if (*ptrFormat == '+')
        {  // The previous number or greater.
          ptrFormat++;
          maximum = 2000000000; // effectively infinite
        }
        else
        {
          maximum = minimum;
        }
        assert(*ptrFormat == '?');
        ptrFormat++;  // Point to first character after second '?'
        if ((argv[idx].i >= minimum) && (argv[idx].i <= maximum))
        {       // include string
          while ((*ptrFormat != '?') || (*(ptrFormat+1) != '?'))
          {
            assert(*ptrFormat != '\0'); // ensure valid format
            *ptrString = *ptrFormat;
            ptrString++;
            ptrFormat++;
          }
        }
        else
        {       // skip string
          while ((*ptrFormat != '?') || (*(ptrFormat+1) != '?'))
          {
            assert(*ptrFormat != '\0'); // ensure valid format
            ptrFormat++;
          }
        }
        ptrFormat--;  // Adjust for the increment at the end of the loop.
      }
      ptrFormat += 3;
    }
    else
    {
      *ptrString++ = *ptrFormat++;
    }
  }

  *ptrString = '\0';
  *pptrString = ptrString;
}

void showCopyright(char** pptrOutput)
{
  strcpy(*pptrOutput, LITERAL_COPYRIGHT);
  *pptrOutput += strlen(LITERAL_COPYRIGHT);
}

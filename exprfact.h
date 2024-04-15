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
#ifndef _EXPRFACT_H
#define _EXPRFACT_H

#ifdef FACTORIZATION_FUNCTIONS
#define TOKEN_TOTIENT        34
#define TOKEN_SUMDIVS        35
#define TOKEN_NUMDIVS        36
#define TOKEN_MINFACT        37
#define TOKEN_MAXFACT        38
#define TOKEN_NUMFACT        39
#define TOKEN_CONCATFACT     40
#endif
#define TOKEN_PRIMORIAL      41
#define TOKEN_GCD            42
#define TOKEN_LCM            43
#define TOKEN_MODPOW         44
#define TOKEN_MODINV         45
#define TOKEN_MODDIV         46
#define TOKEN_MOD            47
#define TOKEN_FLOORDIV       48
#define TOKEN_SUMDIGITS      49
#define TOKEN_NUMDIGITS      50
#define TOKEN_REVDIGITS      51
#define TOKEN_ISPRIME        52
#define TOKEN_JACOBI         53
#define TOKEN_SQRT           54
#define TOKEN_IROOT          55
#define TOKEN_F              56
#define TOKEN_L              57
#define TOKEN_P              58
#define TOKEN_N              59
#define TOKEN_B              60
#define TOKEN_FACTORIAL1     61
#define TOKEN_PRINT          62
#define TOKEN_PRINT_HEX      63
#define TOKEN_PRINTFACT      64
#define TOKEN_PRINTFACT_HEX  65
#define TOKEN_PRINTPRIME     66
#define TOKEN_PRINTPRIME_HEX 67
#define TOKEN_GET_VAR        68
#define TOKEN_SET_VAR        69
#define TOKEN_IF             70
#define TOKEN_JMP            71
#define TOKEN_SGN1           72
#define START_FLOW_CONTROL  121
#define START_FOR           121
#define START_REPEAT        122
#define START_WHILE         123
#define START_UNTIL         124
#define START_IF            125
#define START_ELSE          126
#define START_STATEMENT     127

#endif

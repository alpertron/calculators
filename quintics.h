//
// This file is part of Alpertron Calculators.
//
// Copyright 2019-2022 Dario Alejandro Alpern
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
#ifndef _QUINTICS_H
#define _QUINTICS_H

#define FIVEexp (1 << 20)
#define Pexp (1 << 15)
#define Qexp (1 << 10)
#define Rexp (1 << 5)
#define Sexp (1 << 0)
#define END_COEFF 0xFFFF
#define END_ARRAY 0

#define P0 0x0000
#define P1 0x1000
#define P2 0x2000
#define P3 0x3000
#define P4 0x4000
#define P5 0x5000
#define P6 0x6000
#define P7 0x7000

#define Q0 0x0000
#define Q1 0x0100
#define Q2 0x0200
#define Q3 0x0300
#define Q4 0x0400
#define Q5 0x0500
#define Q6 0x0600
#define Q8 0x0800

#define R0 0x0000
#define R1 0x0010
#define R2 0x0020
#define R3 0x0030
#define R4 0x0040
#define R5 0x0050
#define R6 0x0060

#define S0 0x0000
#define S1 0x0001
#define S2 0x0002
#define S3 0x0003
#define S4 0x0004
#define S5 0x0005
#define S6 0x0006

struct monomial
{
  int coefficient;
  int exponents;
};

#endif
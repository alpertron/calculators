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

// David Dummit's coefficients for solvable quintics
// See http://www.emba.uvm.edu/~ddummit/quintics/quintics.html

#include "rootseq.h"

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

enum
{
  index_T1 = 0,
  index_T2,
  index_T3,
  index_T4,
  index_l0,
  index_v,
  index_O,
};

struct stQuinticF20
{
  unsigned short pqrs;
  short coeff;
};

static BigRational RatValues[7];  // T1, T2, T3, T4, l0, v and O.
static BigRational RatRoot;       // Rational root of polynomial F20.
static BigRational RatM;
static BigRational RatN;
static BigRational RatR;
static BigRational RatR2;
static BigRational RatS2;
static bool firstNumberShown;

// Coefficients taken from Dummit's Solving Solvable Quintics article.
static struct stQuinticF20 astQuinticF20[] =
{
  // Coefficients of independent term
{P0 + Q8 + R0 + S0, 1},      // q^8
{P1 + Q6 + R1 + S0, -13},    // -13pq^6*r
{P5 + Q2 + R2 + S0, 1},      // p^5*q^2*r^2
{P2 + Q4 + R2 + S0, 65},     // 65p^2*q^4*r^2
{P6 + Q0 + R3 + S0, -4},     // -4p^6*r^3
{P3 + Q2 + R3 + S0, -128},   // -128p^3*q^2*r^3
{P0 + Q4 + R3 + S0, 17},     // 17q^4*r^3
{P4 + Q0 + R4 + S0, 48},     // 48p^4*r^4
{P1 + Q2 + R4 + S0, -16},    // -16pq^2*r^4
{P2 + Q0 + R5 + S0, -192},   // -192p^2*r^5
{P0 + Q0 + R6 + S0, 256},    // 256r^6
{P5 + Q3 + R0 + S1, -4},     // -4p^5*q^3*s
{P2 + Q5 + R0 + S1, -12},    // -12p^2*q^5*s
{P6 + Q1 + R1 + S1, 18},     // 18p^6*qrs
{P3 + Q3 + R1 + S1, 12},     // 12p^3*q^3*rs
{P0 + Q5 + R1 + S1, -124},   // -124q^5*r*s
{P4 + Q1 + R2 + S1, 196},    // 196p^4*qr^2*s
{P1 + Q3 + R2 + S1, 590},    // 590pq^3*r^2*s
{P2 + Q1 + R3 + S1, -160},   // -160p^2*qr^3*s
{P0 + Q1 + R4 + S1, -1600},  // -1600qr^4*s
{P7 + Q0 + R0 + S2, -27},    // -27p^7*s^2
{P4 + Q2 + R0 + S2, -150},   // -150p^4*q^2*s^2
{P1 + Q4 + R0 + S2, -125},   // -125pq^4*s^2
{P5 + Q0 + R1 + S2, -99},    // -99p^5*rs^2
{P2 + Q2 + R1 + S2, -725},   // -725p^2*q^2*rs^2
{P3 + Q0 + R2 + S2, 1200},   // 1200p^3*r^2*s^2
{P0 + Q2 + R2 + S2, 3250},   // 3250q^2*r^2*s^2
{P1 + Q0 + R3 + S2, -2000},  // -2000pr^3*s^2
{P1 + Q1 + R1 + S3, -1250},  // -1250pqrs^3
{P2 + Q0 + R0 + S4, 3125},   // 3125p^2*s^4
{P0 + Q0 + R1 + S4, -9375},  // -9375rs^4
{END_COEFF, 0},
// Coefficients of x
{P1 + Q6 + R0 + S0, -2},     // -2pq^6
{P2 + Q4 + R1 + S0, 19},     // 19p^2*q^4*r
{P3 + Q2 + R2 + S0, -51},    // -51p^3*q^2*r^2
{P0 + Q4 + R2 + S0, 3},      // 3q^4*r^2
{P4 + Q0 + R3 + S0, 32},     // 32p^4*r^3
{P1 + Q2 + R3 + S0, 76},     // 76pq^2*r^3
{P2 + Q0 + R4 + S0, -256},   // -256p^2*r^4
{P0 + Q0 + R5 + S0, 512},    // 512r^5
{P3 + Q3 + R0 + S1, -31},    // -31p^3*q^3*s
{P0 + Q5 + R0 + S1, -58},    // -58q^5*s
{P4 + Q1 + R1 + S1, 117},    // 117p^4*qrs
{P1 + Q3 + R1 + S1, 105},    // 105pq^3*rs
{P2 + Q1 + R2 + S1, 260},    // 260p^2*qr^2*s
{P0 + Q1 + R3 + S1, -2400},  // -2400qr^3*s
{P5 + R0 + R0 + S2, -108},   // -108p^5*s^2
{P2 + Q2 + R0 + S2, -325},   // -325p^2*q^2*s^2
{P3 + Q0 + R1 + S2, 525},    // 525p^3*rs^2
{P0 + Q2 + R1 + S2, 2750},   // 2750q^2*rs^2
{P1 + Q0 + R2 + S2, -500},   // -500pr^2*s^2
{P1 + Q1 + R0 + S3, 625},    // 625pqs^3
{P0 + Q0 + R0 + S4, -3125},  // -3125s^4
{END_COEFF, 0},
// Coefficients of x^2
{P2 + Q4 + R0 + S0, 1},      // p^2*q^4
{P3 + Q2 + R1 + S0, -6},     // -6p^3*q^2*r
{P0 + Q4 + R1 + S0, -8},     // -8q^4*r
{P4 + Q0 + R2 + S0, 9},      // 9p^4*r^2
{P1 + Q2 + R2 + S0, 76},     // 76pq^2*r^2
{P2 + Q0 + R3 + S0, -136},   // -136p^2*r^3
{P0 + Q0 + R4 + S0, 400},    // 400r^4
{P1 + Q3 + R0 + S1, -50},    // -50pq^3*s
{P2 + Q1 + R1 + S1, 90},     // 90p^2*qrs
{P0 + Q1 + R2 + S1, -1400},  // -1400qr^2*s
{P0 + Q2 + R0 + S2, 625},    // 625q^2*s^2
{P1 + Q0 + R1 + S2, 500},    // 500prs^2
{END_COEFF, 0},
// Coefficients of x^3
{P0 + Q4 + R0 + S0, -2},     // -2q^4
{P1 + Q2 + R1 + S0, 21},     // 21pq^2*r
{P2 + Q0 + R2 + S0, -40},    // -40p^2*r^2
{P0 + Q0 + R3 + S0, 160},    // 160r^3
{P2 + Q1 + R0 + S1, -15},    // -15p^2*qs
{P0 + Q1 + R1 + S1, -400},   // -400qrs
{P1 + Q0 + R0 + S2, 125},     // 125ps^2
{END_COEFF, 0},
// Coefficients of x^4
{P1 + Q2 + R0 + S0, 2},      // 2pq^2
{P2 + Q0 + R1 + S0, -6},     // -6p^2*r
{P0 + Q0 + R2 + S0, 40},     // 40r^2
{P0 + Q1 + R0 + S1, -50},    // -50qs
{END_COEFF, 0},
// Coefficients of x^5
{P0 + Q0 + R1 + S0, 8},      // 8r
{END_COEFF, 0},
// Coefficients of x^6
{P0 + Q0 + R0 + S0, 1},      // 1
{END_COEFF, -1},
};

struct monomial
{
  int coefficient;
	int exponents;
};

struct monomial arrayF[] =
{
  {4, FIVEexp*0 + Pexp * 6 + Qexp * 6},
  {59, FIVEexp*0 + Pexp * 3 + Qexp * 8},
  {216, FIVEexp*0 + Qexp * 10},
  {-36, FIVEexp*0 + Pexp * 7 + Qexp * 4 + Rexp},
  {-623, FIVEexp*0 + Pexp * 4 + Qexp * 6 + Rexp},
  {-522, FIVEexp*1 + Pexp + Qexp * 8 + Rexp},
  {81, FIVEexp*0 + Pexp * 8 + Qexp * 2 + Rexp * 2},
  {403, FIVEexp*1 + Pexp * 5 + Qexp * 4 + Rexp * 2},
  {433, FIVEexp*2 + Pexp * 2 + Qexp * 6 + Rexp * 2},
  {-72, FIVEexp*2 + Pexp * 6 + Qexp * 2 + Rexp * 3},
  {-28, FIVEexp*4 + Pexp * 3 + Qexp * 4 + Rexp * 3},
  {1, FIVEexp*4 + Qexp * 6 + Rexp * 3},
  {16, FIVEexp*4 + Pexp * 4 + Qexp * 2 + Rexp * 4},
  {108, FIVEexp*0 + Pexp * 8 + Qexp * 3 + Sexp},
  {1584, FIVEexp*0 + Pexp * 5 + Qexp * 5 + Sexp},
  {228, FIVEexp*2 + Pexp * 2 + Qexp * 7 + Sexp},
  {-486, FIVEexp*0 + Pexp * 9 + Qexp + Rexp + Sexp},
  {-1944, FIVEexp*1 + Pexp * 6 + Qexp * 3 + Rexp + Sexp},
  {-1802, FIVEexp*2 + Pexp * 3 + Qexp * 5 + Rexp + Sexp},
  {-72, FIVEexp*3 + Qexp * 7 + Rexp + Sexp},
  {432, FIVEexp*2 + Pexp * 7 + Qexp + Rexp * 2 + Sexp},
  {148, FIVEexp*4 + Pexp * 4 + Qexp * 3 + Rexp * 2 + Sexp},
  {52, FIVEexp*4 + Pexp + Qexp * 5 + Rexp * 2 + Sexp},
  {-96, FIVEexp*4 + Pexp * 5 + Qexp + Rexp * 3 + Sexp},
  {-16, FIVEexp*5 + Pexp * 2 + Qexp * 3 + Rexp * 3 + Sexp},
  {729, FIVEexp*0 + Pexp * 10 + Sexp * 2},
  {486, FIVEexp*2 + Pexp * 7 + Qexp * 2 + Sexp * 2},
  {96, FIVEexp*4 + Pexp * 4 + Qexp * 4 + Sexp * 2},
  {6, FIVEexp*6 + Pexp + Qexp * 6 + Sexp * 2},
  {-729, FIVEexp*2 + Pexp * 8 + Rexp + Sexp * 2},
  {-1404, FIVEexp*3 + Pexp * 5 + Qexp * 2 + Rexp + Sexp * 2},
  {-153, FIVEexp*5 + Pexp * 2 + Qexp * 4 + Rexp + Sexp * 2},
  {216, FIVEexp*4 + Pexp * 6 + Rexp * 2 + Sexp * 2},
  {272, FIVEexp*5 + Pexp * 3 + Qexp * 2 + Rexp * 2 + Sexp * 2},
  {1, FIVEexp*6 + Qexp * 4 + Rexp * 2 + Sexp * 2},
  {-16, FIVEexp*6 + Pexp * 4 + Rexp * 3 + Sexp * 2},
  {72, FIVEexp*5 + Pexp * 3 + Qexp * 3 + Sexp * 3},
  {56, FIVEexp*5 + Qexp * 5 + Sexp * 3},
  {-324, FIVEexp*5 + Pexp * 4 + Qexp + Rexp + Sexp * 3},
  {-76, FIVEexp*6 + Pexp + Qexp * 3 + Rexp + Sexp * 3},
  {16, FIVEexp*7 + Pexp * 2 + Qexp + Rexp * 2 + Sexp * 3},
  {297, FIVEexp*5 + Pexp * 5 + Sexp * 4},
  {24, FIVEexp*7 + Pexp * 2 + Qexp * 2 + Sexp * 4},
  {-36, FIVEexp*7 + Pexp * 3 + Rexp + Sexp * 4},
  {-1, FIVEexp*8 + Qexp * 2 + Rexp + Sexp * 4},
  {-1, FIVEexp*10 + Sexp * 6},
  {0, END_ARRAY},
};

struct monomial arrayD[] =
{
  {-4, FIVEexp * 0 + Pexp * 3 + Qexp * 2 + Rexp * 2},
  {-27, FIVEexp * 0 + Qexp * 4 + Rexp * 2},
  {16, FIVEexp * 0 + Pexp * 4 + Rexp * 3},
  {144, FIVEexp*0 + Pexp + Qexp * 2 + Rexp * 3},
  {-128, FIVEexp*0 + Pexp * 2 + Rexp * 4},
  {256, FIVEexp*0 + Rexp * 5},
  {16, FIVEexp*0 + Pexp * 3 + Qexp * 3 + Sexp},
  {108, FIVEexp*0 + Qexp * 5 + Sexp},
  {-72, FIVEexp*0 + Pexp * 4 + Qexp + Rexp + Sexp},
  {-126, FIVEexp*1 + Pexp + Qexp * 3 + Rexp + Sexp},
  {112, FIVEexp*1 + Pexp * 2 + Qexp + Rexp * 2 + Sexp},
  {-64, FIVEexp*2 + Qexp + Rexp * 3 + Sexp},
  {108, FIVEexp*0 + Pexp * 5 + Sexp * 2},
  {33, FIVEexp*2 + Pexp * 2 + Qexp * 2 + Sexp * 2},
  {-36, FIVEexp*2 + Pexp * 3 + Rexp + Sexp * 2},
  {18, FIVEexp*3 + Qexp * 2 + Rexp + Sexp * 2},
  {16, FIVEexp*3 + Pexp + Rexp * 2 + Sexp * 2},
  {-6, FIVEexp*4 + Pexp + Qexp + Sexp * 3},
  {1, FIVEexp*5 + Sexp * 4},
  {0, END_ARRAY},
};

struct monomial arrayB[] =
{
  // b[1, 0]
  {4, FIVEexp*2 + Pexp * 7 + Qexp * 7},
  {87, FIVEexp*2 + Pexp * 4 + Qexp * 9},
  {84, FIVEexp*3 + Pexp + Qexp * 11},
  {-44, FIVEexp*2 + Pexp * 8 + Qexp * 5 + Rexp},
  {-1119, FIVEexp*2 + Pexp * 5 + Qexp * 7 + Rexp},
  {-6118, FIVEexp*2 + Pexp * 2 + Qexp * 9 + Rexp},
  {33, FIVEexp*3 + Pexp * 9 + Qexp * 3 + Rexp * 2},
  {1031, FIVEexp*3 + Pexp * 6 + Qexp * 5 + Rexp * 2},
  {33221, FIVEexp*2 + Pexp * 3 + Qexp * 7 + Rexp * 2},
  {-2378, FIVEexp*2 + Qexp * 9 + Rexp * 2},
  {-216, FIVEexp*2 + Pexp * 10 + Qexp + Rexp * 3},
  {-9752, FIVEexp*2 + Pexp * 7 + Qexp * 3 + Rexp * 3},
  {-83306, FIVEexp*2 + Pexp * 4 + Qexp * 5 + Rexp * 3},
  {13357, FIVEexp*2 + Pexp + Qexp * 7 + Rexp * 3},
  {5568, FIVEexp*2 + Pexp * 8 + Qexp + Rexp * 4},
  {19248, FIVEexp*3 + Pexp * 5 + Qexp * 3 + Rexp * 4},
  {4904, FIVEexp*2 + Pexp * 2 + Qexp * 5 + Rexp * 4},
  {-50176, FIVEexp*2 + Pexp * 6 + Qexp + Rexp * 5},
  {-30208, FIVEexp*3 + Pexp * 3 + Qexp * 3 + Rexp * 5},
  {-14656, FIVEexp*3 + Qexp * 5 + Rexp * 5},
  {37888, FIVEexp*3 + Pexp * 4 + Qexp + Rexp * 6},
  {10752, FIVEexp*4 + Pexp + Qexp * 3 + Rexp * 6},
  {-2048, FIVEexp*5 + Pexp * 2 + Qexp + Rexp * 7},
  {36, FIVEexp*2 + Pexp * 9 + Qexp * 4 + Sexp},
  {1496, FIVEexp*2 + Pexp * 6 + Qexp * 6 + Sexp},
  {2253, FIVEexp*3 + Pexp * 3 + Qexp * 8 + Sexp},
  {696, FIVEexp*4 + Qexp * 10 + Sexp},
  {-54, FIVEexp*3 + Pexp * 10 + Qexp * 2 + Rexp + Sexp},
  {-12892, FIVEexp*2 + Pexp * 7 + Qexp * 4 + Rexp + Sexp},
  {-108743, FIVEexp*2 + Pexp * 4 + Qexp * 6 + Rexp + Sexp},
  {-33714, FIVEexp*3 + Pexp + Qexp * 8 + Rexp + Sexp},
  {648, FIVEexp*2 + Pexp * 11 + Rexp * 2 + Sexp},
  {34371, FIVEexp*2 + Pexp * 8 + Qexp * 2 + Rexp * 2 + Sexp},
  {357019, FIVEexp*2 + Pexp * 5 + Qexp * 4 + Rexp * 2 + Sexp},
  {115423, FIVEexp*3 + Pexp * 2 + Qexp * 6 + Rexp * 2 + Sexp},
  {-18144, FIVEexp*2 + Pexp * 9 + Rexp * 3 + Sexp},
  {-401536, FIVEexp*2 + Pexp * 6 + Qexp * 2 + Rexp * 3 + Sexp},
  {-27836, FIVEexp*4 + Pexp * 3 + Qexp * 4 + Rexp * 3 + Sexp},
  {18133, FIVEexp*4 + Qexp * 6 + Rexp * 3 + Sexp},
  {178048, FIVEexp*2 + Pexp * 7 + Rexp * 4 + Sexp},
  {5072, FIVEexp*5 + Pexp * 4 + Qexp * 2 + Rexp * 4 + Sexp},
  {-2176, FIVEexp*6 + Pexp + Qexp * 4 + Rexp * 4 + Sexp},
  {-143872, FIVEexp*3 + Pexp * 5 + Rexp * 5 + Sexp},
  {128, FIVEexp*7 + Pexp * 2 + Qexp * 2 + Rexp * 5 + Sexp},
  {8192, FIVEexp*5 + Pexp * 3 + Rexp * 6 + Sexp},
  {512, FIVEexp*6 + Qexp * 2 + Rexp * 6 + Sexp},
  {-243, FIVEexp*2 + Pexp * 11 + Qexp + Sexp * 2},
  {666, FIVEexp*3 + Pexp * 8 + Qexp * 3 + Sexp * 2},
  {2052, FIVEexp*4 + Pexp * 5 + Qexp * 5 + Sexp * 2},
  {916, FIVEexp*5 + Pexp * 2 + Qexp * 7 + Sexp * 2},
  {-28971, FIVEexp*2 + Pexp * 9 + Qexp + Rexp + Sexp * 2},
  {-78458, FIVEexp*3 + Pexp * 6 + Qexp * 3 + Rexp + Sexp * 2},
  {-45399, FIVEexp*4 + Pexp * 3 + Qexp * 5 + Rexp + Sexp * 2},
  {-7108, FIVEexp*5 + Qexp * 7 + Rexp + Sexp * 2},
  {71856, FIVEexp*3 + Pexp * 7 + Qexp + Rexp * 2 + Sexp * 2},
  {12672, FIVEexp*5 + Pexp * 4 + Qexp * 3 + Rexp * 2 + Sexp * 2},
  {19759, FIVEexp*5 + Pexp + Qexp * 5 + Rexp * 2 + Sexp * 2},
  {1616, FIVEexp*4 + Pexp * 5 + Qexp + Rexp * 3 + Sexp * 2},
  {64, FIVEexp*6 + Pexp * 2 + Qexp * 3 + Rexp * 3 + Sexp * 2},
  {-4992, FIVEexp*6 + Pexp * 3 + Qexp + Rexp * 4 + Sexp * 2},
  {-384, FIVEexp*7 + Qexp * 3 + Rexp * 4 + Sexp * 2},
  {-1024, FIVEexp*7 + Pexp + Qexp + Rexp * 5 + Sexp * 2},
  {243, FIVEexp*5 + Pexp * 10 + Sexp * 3},
  {3132, FIVEexp*5 + Pexp * 7 + Qexp * 2 + Sexp * 3},
  {4, FIVEexp*10 + Pexp * 4 + Qexp * 4 + Sexp * 3},
  {134, FIVEexp*8 + Pexp + Qexp * 6 + Sexp * 3},
  {-19683, FIVEexp*4 + Pexp * 8 + Rexp + Sexp * 3},
  {-31416, FIVEexp*5 + Pexp * 5 + Qexp * 2 + Rexp + Sexp * 3},
  {-2881, FIVEexp*7 + Pexp * 2 + Qexp * 4 + Rexp + Sexp * 3},
  {17568, FIVEexp*5 + Pexp * 6 + Rexp * 2 + Sexp * 3},
  {3968, FIVEexp*7 + Pexp * 3 + Qexp * 2 + Rexp * 2 + Sexp * 3},
  {101, FIVEexp*7 + Qexp * 4 + Rexp * 2 + Sexp * 3},
  {-656, FIVEexp*7 + Pexp * 4 + Rexp * 3 + Sexp * 3},
  {5376, FIVEexp*7 + Pexp + Qexp * 2 + Rexp * 3 + Sexp * 3},
  {-1408, FIVEexp*7 + Pexp * 2 + Rexp * 4 + Sexp * 3},
  {512, FIVEexp*8 + Rexp * 5 + Sexp * 3},
  {-27, FIVEexp*7 + Pexp * 6 + Qexp + Sexp * 4},
  {54, FIVEexp*8 + Pexp * 3 + Qexp * 3 + Sexp * 4},
  {46, FIVEexp*9 + Qexp * 5 + Sexp * 4},
  {-2334, FIVEexp*7 + Pexp * 4 + Qexp + Rexp + Sexp * 4},
  {-1877, FIVEexp*8 + Pexp + Qexp * 3 + Rexp + Sexp * 4},
  {504, FIVEexp*8 + Pexp * 2 + Qexp + Rexp * 2 + Sexp * 4},
  {-576, FIVEexp*9 + Qexp + Rexp * 3 + Sexp * 4},
  {81, FIVEexp*9 + Pexp * 5 + Sexp * 5},
  {58, FIVEexp*10 + Pexp * 2 + Qexp * 2 + Sexp * 5},
  {-52, FIVEexp*9 + Pexp * 3 + Rexp + Sexp * 5},
  {171, FIVEexp*10 + Qexp * 2 + Rexp + Sexp * 5},
  {-128, FIVEexp*10 + Pexp + Rexp * 2 + Sexp * 5},
  {1, FIVEexp*13 + Pexp + Qexp + Sexp * 6},
  {-1, FIVEexp*14 + Sexp * 7},
  {0, END_COEFF},

  // b[1, 1]
  {-8, FIVEexp*3 + Pexp * 5 + Qexp * 7},
  {-58, FIVEexp*3 + Pexp * 2 + Qexp * 9},
  {432, FIVEexp*2 + Pexp * 6 + Qexp * 5 + Rexp},
  {3876, FIVEexp*2 + Pexp * 3 + Qexp * 7 + Rexp},
  {84, FIVEexp*4 + Qexp * 9 + Rexp},
  {-1496, FIVEexp*2 + Pexp * 7 + Qexp * 3 + Rexp * 2},
  {-18834, FIVEexp*2 + Pexp * 4 + Qexp * 5 + Rexp * 2},
  {-25624, FIVEexp*2 + Pexp + Qexp * 7 + Rexp * 2},
  {1584, FIVEexp*2 + Pexp * 8 + Qexp + Rexp * 3},
  {39344, FIVEexp*2 + Pexp * 5 + Qexp * 3 + Rexp * 3},
  {113924, FIVEexp*2 + Pexp * 2 + Qexp * 5 + Rexp * 3},
  {-32576, FIVEexp*2 + Pexp * 6 + Qexp + Rexp * 4},
  {-48608, FIVEexp*3 + Pexp * 3 + Qexp * 3 + Rexp * 4},
  {-18464, FIVEexp*3 + Qexp * 5 + Rexp * 4},
  {40192, FIVEexp*3 + Pexp * 4 + Qexp + Rexp * 5},
  {15488, FIVEexp*4 + Pexp + Qexp * 3 + Rexp * 5},
  {-3072, FIVEexp*5 + Pexp * 2 + Qexp + Rexp * 6},
  {-552, FIVEexp*2 + Pexp * 7 + Qexp * 4 + Sexp},
  {-3786, FIVEexp*2 + Pexp * 4 + Qexp * 6 + Sexp},
  {212, FIVEexp*3 + Pexp + Qexp * 8 + Sexp},
  {3456, FIVEexp*2 + Pexp * 8 + Qexp * 2 + Rexp + Sexp},
  {6532, FIVEexp*3 + Pexp * 5 + Qexp * 4 + Rexp + Sexp},
  {412, FIVEexp*4 + Pexp * 2 + Qexp * 6 + Rexp + Sexp},
  {-3672, FIVEexp*2 + Pexp * 9 + Rexp * 2 + Sexp},
  {-74148, FIVEexp*2 + Pexp * 6 + Qexp * 2 + Rexp * 2 + Sexp},
  {-1008, FIVEexp*4 + Pexp * 3 + Qexp * 4 + Rexp * 2 + Sexp},
  {14354, FIVEexp*4 + Qexp * 6 + Rexp * 2 + Sexp},
  {82848, FIVEexp*2 + Pexp * 7 + Rexp * 3 + Sexp},
  {11584, FIVEexp*4 + Pexp * 4 + Qexp * 2 + Rexp * 3 + Sexp},
  {-376, FIVEexp*7 + Pexp + Qexp * 4 + Rexp * 3 + Sexp},
  {-115328, FIVEexp*3 + Pexp * 5 + Rexp * 4 + Sexp},
  {1664, FIVEexp*5 + Pexp * 2 + Qexp * 2 + Rexp * 4 + Sexp},
  {9728, FIVEexp*5 + Pexp * 3 + Rexp * 5 + Sexp},
  {768, FIVEexp*6 + Qexp * 2 + Rexp * 5 + Sexp},
  {-2592, FIVEexp*2 + Pexp * 9 + Qexp + Sexp * 2},
  {-4536, FIVEexp*3 + Pexp * 6 + Qexp * 3 + Sexp * 2},
  {-2648, FIVEexp*4 + Pexp * 3 + Qexp * 5 + Sexp * 2},
  {-2236, FIVEexp*5 + Qexp * 7 + Sexp * 2},
  {-108, FIVEexp*5 + Pexp * 7 + Qexp + Rexp + Sexp * 2},
  {-2708, FIVEexp*5 + Pexp * 4 + Qexp * 3 + Rexp + Sexp * 2},
  {372, FIVEexp*6 + Pexp + Qexp * 5 + Rexp + Sexp * 2},
  {39888, FIVEexp*4 + Pexp * 5 + Qexp + Rexp * 2 + Sexp * 2},
  {4424, FIVEexp*6 + Pexp * 2 + Qexp * 3 + Rexp * 2 + Sexp * 2},
  {-6624, FIVEexp*6 + Pexp * 3 + Qexp + Rexp * 3 + Sexp * 2},
  {-384, FIVEexp*7 + Qexp * 3 + Rexp * 3 + Sexp * 2},
  {-1152, FIVEexp*7 + Pexp + Qexp + Rexp * 4 + Sexp * 2},
  {1134, FIVEexp*4 + Pexp * 8 + Sexp * 3},
  {1728, FIVEexp*5 + Pexp * 5 + Qexp * 2 + Sexp * 3},
  {-114, FIVEexp*7 + Pexp * 2 + Qexp * 4 + Sexp * 3},
  {-1188, FIVEexp*6 + Pexp * 6 + Rexp + Sexp * 3},
  {8, FIVEexp*7 + Pexp * 3 + Qexp * 2 + Rexp + Sexp * 3},
  {-76, FIVEexp*8 + Qexp * 4 + Rexp + Sexp * 3},
  {192, FIVEexp*8 + Pexp * 4 + Rexp * 2 + Sexp * 3},
  {5328, FIVEexp*7 + Pexp + Qexp * 2 + Rexp * 2 + Sexp * 3},
  {-768, FIVEexp*7 + Pexp * 2 + Rexp * 3 + Sexp * 3},
  {768, FIVEexp*8 + Rexp * 4 + Sexp * 3},
  {-918, FIVEexp*7 + Pexp * 4 + Qexp + Sexp * 4},
  {-484, FIVEexp*8 + Pexp + Qexp * 3 + Sexp * 4},
  {-108, FIVEexp*9 + Pexp * 2 + Qexp + Rexp + Sexp * 4},
  {-608, FIVEexp*9 + Qexp + Rexp * 2 + Sexp * 4},
  {96, FIVEexp*9 + Pexp * 3 + Sexp * 5},
  {82, FIVEexp*10 + Qexp * 2 + Sexp * 5},
  {8, FIVEexp*11 + Pexp + Rexp + Sexp * 5},
  {0, END_COEFF},

  // b[1, 2]
  {4, FIVEexp*3 + Pexp * 6 + Qexp * 5},
  {254, FIVEexp*2 + Pexp * 3 + Qexp * 7},
  {792, FIVEexp*2 + Qexp * 9},
  {-6, FIVEexp*4 + Pexp * 7 + Qexp * 3 + Rexp},
  {-2604, FIVEexp*2 + Pexp * 4 + Qexp * 5 + Rexp},
  {-10598, FIVEexp*2 + Pexp + Qexp * 7 + Rexp},
  {54, FIVEexp*3 + Pexp * 8 + Qexp + Rexp * 2},
  {8362, FIVEexp*2 + Pexp * 5 + Qexp * 3 + Rexp * 2},
  {9738, FIVEexp*3 + Pexp * 2 + Qexp * 5 + Rexp * 2},
  {-1752, FIVEexp*3 + Pexp * 6 + Qexp + Rexp * 3},
  {-4016, FIVEexp*4 + Pexp * 3 + Qexp * 3 + Rexp * 3},
  {-8788, FIVEexp*3 + Qexp * 5 + Rexp * 3},
  {16544, FIVEexp*3 + Pexp * 4 + Qexp + Rexp * 4},
  {8096, FIVEexp*4 + Pexp + Qexp * 3 + Rexp * 4},
  {-1664, FIVEexp*5 + Pexp * 2 + Qexp + Rexp * 5},
  {54, FIVEexp*3 + Pexp * 8 + Qexp * 2 + Sexp},
  {3854, FIVEexp*2 + Pexp * 5 + Qexp * 4 + Sexp},
  {2768, FIVEexp*3 + Pexp * 2 + Qexp * 6 + Sexp},
  {-162, FIVEexp*3 + Pexp * 9 + Rexp + Sexp},
  {-18396, FIVEexp*2 + Pexp * 6 + Qexp * 2 + Rexp + Sexp},
  {-2926, FIVEexp*4 + Pexp * 3 + Qexp * 4 + Rexp + Sexp},
  {4688, FIVEexp*4 + Qexp * 6 + Rexp + Sexp},
  {4752, FIVEexp*3 + Pexp * 7 + Rexp * 2 + Sexp},
  {6882, FIVEexp*4 + Pexp * 4 + Qexp * 2 + Rexp * 2 + Sexp},
  {-698, FIVEexp*6 + Pexp + Qexp * 4 + Rexp * 2 + Sexp},
  {-42016, FIVEexp*3 + Pexp * 5 + Rexp * 3 + Sexp},
  {464, FIVEexp*5 + Pexp * 2 + Qexp * 2 + Rexp * 3 + Sexp},
  {4096, FIVEexp*5 + Pexp * 3 + Rexp * 4 + Sexp},
  {416, FIVEexp*6 + Qexp * 2 + Rexp * 4 + Sexp},
  {-594, FIVEexp*3 + Pexp * 7 + Qexp + Sexp * 2},
  {-454, FIVEexp*5 + Pexp * 4 + Qexp * 3 + Sexp * 2},
  {-1906, FIVEexp*5 + Pexp + Qexp * 5 + Sexp * 2},
  {6876, FIVEexp*4 + Pexp * 5 + Qexp + Rexp + Sexp * 2},
  {1914, FIVEexp*6 + Pexp * 2 + Qexp * 3 + Rexp + Sexp * 2},
  {-2016, FIVEexp*6 + Pexp * 3 + Qexp + Rexp * 2 + Sexp * 2},
  {-32, FIVEexp*8 + Qexp * 3 + Rexp * 2 + Sexp * 2},
  {-448, FIVEexp*7 + Pexp + Qexp + Rexp * 3 + Sexp * 2},
  {-432, FIVEexp*5 + Pexp * 6 + Sexp * 3},
  {-78, FIVEexp*7 + Pexp * 3 + Qexp * 2 + Sexp * 3},
  {-224, FIVEexp*7 + Qexp * 4 + Sexp * 3},
  {18, FIVEexp*8 + Pexp * 4 + Rexp + Sexp * 3},
  {1636, FIVEexp*7 + Pexp + Qexp * 2 + Rexp + Sexp * 3},
  {-48, FIVEexp*8 + Pexp * 2 + Rexp * 2 + Sexp * 3},
  {416, FIVEexp*8 + Rexp * 3 + Sexp * 3},
  {-276, FIVEexp*8 + Pexp * 2 + Qexp + Sexp * 4},
  {-236, FIVEexp*9 + Qexp + Rexp + Sexp * 4},
  {22, FIVEexp*10 + Pexp + Sexp * 5},
  {0, END_COEFF},

  // b[1, 3]
  {-78, FIVEexp*2 + Pexp * 4 + Qexp * 5},
  {-564, FIVEexp*2 + Pexp + Qexp * 7},
  {574, FIVEexp*2 + Pexp * 5 + Qexp * 3 + Rexp},
  {5024, FIVEexp*2 + Pexp * 2 + Qexp * 5 + Rexp},
  {-1116, FIVEexp*2 + Pexp * 6 + Qexp + Rexp * 2},
  {-3218, FIVEexp*3 + Pexp * 3 + Qexp * 3 + Rexp * 2},
  {-2306, FIVEexp*3 + Qexp * 5 + Rexp * 2},
  {3488, FIVEexp*3 + Pexp * 4 + Qexp + Rexp * 3},
  {2152, FIVEexp*4 + Pexp + Qexp * 3 + Rexp * 3},
  {-448, FIVEexp*5 + Pexp * 2 + Qexp + Rexp * 4},
  {-378, FIVEexp*2 + Pexp * 6 + Qexp * 2 + Sexp},
  {2, FIVEexp*4 + Pexp * 3 + Qexp * 4 + Sexp},
  {744, FIVEexp*4 + Qexp * 6 + Sexp},
  {1998, FIVEexp*2 + Pexp * 7 + Rexp + Sexp},
  {484, FIVEexp*4 + Pexp * 4 + Qexp * 2 + Rexp + Sexp},
  {-22, FIVEexp*7 + Pexp + Qexp * 4 + Rexp + Sexp},
  {-6672, FIVEexp*3 + Pexp * 5 + Rexp * 2 + Sexp},
  {-28, FIVEexp*6 + Pexp * 2 + Qexp * 2 + Rexp * 2 + Sexp},
  {992, FIVEexp*5 + Pexp * 3 + Rexp * 3 + Sexp},
  {112, FIVEexp*6 + Qexp * 2 + Rexp * 3 + Sexp},
  {468, FIVEexp*4 + Pexp * 5 + Qexp + Sexp * 2},
  {124, FIVEexp*6 + Pexp * 2 + Qexp * 3 + Sexp * 2},
  {-214, FIVEexp*6 + Pexp * 3 + Qexp + Rexp + Sexp * 2},
  {-24, FIVEexp*7 + Qexp * 3 + Rexp + Sexp * 2},
  {-104, FIVEexp*7 + Pexp + Qexp + Rexp * 2 + Sexp * 2},
  {18, FIVEexp*7 + Pexp * 4 + Sexp * 3},
  {158, FIVEexp*7 + Pexp + Qexp * 2 + Sexp * 3},
  {-68, FIVEexp*7 + Pexp * 2 + Rexp + Sexp * 3},
  {112, FIVEexp*8 + Rexp * 2 + Sexp * 3},
  {-38, FIVEexp*9 + Qexp + Sexp * 4},
  {0, END_COEFF},

  // b[1, 4]
  {12, FIVEexp*2 + Pexp * 5 + Qexp * 3},
  {86, FIVEexp*2 + Pexp * 2 + Qexp * 5},
  {-54, FIVEexp*2 + Pexp * 6 + Qexp + Rexp},
  {-172, FIVEexp*3 + Pexp * 3 + Qexp * 3 + Rexp},
  {-492, FIVEexp*3 + Qexp * 5 + Rexp},
  {336, FIVEexp*3 + Pexp * 4 + Qexp + Rexp * 2},
  {464, FIVEexp*4 + Pexp + Qexp * 3 + Rexp * 2},
  {-96, FIVEexp*5 + Pexp * 2 + Qexp + Rexp * 3},
  {162, FIVEexp*2 + Pexp * 7 + Sexp},
  {72, FIVEexp*4 + Pexp * 4 + Qexp * 2 + Sexp},
  {8, FIVEexp*6 + Pexp + Qexp * 4 + Sexp},
  {-864, FIVEexp*3 + Pexp * 5 + Rexp + Sexp},
  {-206, FIVEexp*5 + Pexp * 2 + Qexp * 2 + Rexp + Sexp},
  {224, FIVEexp*5 + Pexp * 3 + Rexp * 2 + Sexp},
  {24, FIVEexp*6 + Qexp * 2 + Rexp * 2 + Sexp},
  {6, FIVEexp*6 + Pexp * 3 + Qexp + Sexp * 2},
  {4, FIVEexp*7 + Qexp * 3 + Sexp * 2},
  {-24, FIVEexp*7 + Pexp + Qexp + Rexp + Sexp * 2},
  {18, FIVEexp*7 + Pexp * 2 + Sexp * 3},
  {24, FIVEexp*8 + Rexp + Sexp * 3},
  {0, END_COEFF},

  // b[1, 5]
  {-2, FIVEexp*4 + Pexp * 3 + Qexp * 3},
  {-72, FIVEexp*3 + Qexp * 5},
  {36, FIVEexp*3 + Pexp * 4 + Qexp + Rexp},
  {74, FIVEexp*4 + Pexp + Qexp * 3 + Rexp},
  {-16, FIVEexp*5 + Pexp * 2 + Qexp + Rexp * 2},
  {-54, FIVEexp*3 + Pexp * 5 + Sexp},
  {-14, FIVEexp*5 + Pexp * 2 + Qexp * 2 + Sexp},
  {24, FIVEexp*5 + Pexp * 3 + Rexp + Sexp},
  {4, FIVEexp*6 + Qexp * 2 + Rexp + Sexp},
  {-2, FIVEexp*7 + Pexp + Qexp + Sexp * 2},
  {4, FIVEexp*8 + Sexp * 3},
  {0, END_COEFF},

  // b[2, 0]
  {8, FIVEexp*2 + Pexp * 6 + Qexp * 11},
  {-2, FIVEexp*3 + Pexp * 3 + Qexp * 13},
  {-432, FIVEexp*2 + Qexp * 15},
  {-156, FIVEexp*2 + Pexp * 7 + Qexp * 9 + Rexp},
  {-133, FIVEexp*2 + Pexp * 4 + Qexp * 11 + Rexp},
  {7272, FIVEexp*2 + Pexp + Qexp * 13 + Rexp},
  {1078, FIVEexp*2 + Pexp * 8 + Qexp * 7 + Rexp * 2},
  {557, FIVEexp*3 + Pexp * 5 + Qexp * 9 + Rexp * 2},
  {-48578, FIVEexp*2 + Pexp * 2 + Qexp * 11 + Rexp * 2},
  {-3149, FIVEexp*2 + Pexp * 9 + Qexp * 5 + Rexp * 3},
  {-14747, FIVEexp*2 + Pexp * 6 + Qexp * 7 + Rexp * 3},
  {166653, FIVEexp*2 + Pexp * 3 + Qexp * 9 + Rexp * 3},
  {45244, FIVEexp*2 + Qexp * 11 + Rexp * 3},
  {2936, FIVEexp*2 + Pexp * 10 + Qexp * 3 + Rexp * 4},
  {26478, FIVEexp*2 + Pexp * 7 + Qexp * 5 + Rexp * 4},
  {-366078, FIVEexp*2 + Pexp * 4 + Qexp * 7 + Rexp * 4},
  {-665323, FIVEexp*2 + Pexp + Qexp * 9 + Rexp * 4},
  {288, FIVEexp*3 + Pexp * 11 + Qexp + Rexp * 5},
  {5424, FIVEexp*2 + Pexp * 8 + Qexp * 3 + Rexp * 5},
  {692856, FIVEexp*2 + Pexp * 5 + Qexp * 5 + Rexp * 5},
  {3413532, FIVEexp*2 + Pexp * 2 + Qexp * 7 + Rexp * 5},
  {-6656, FIVEexp*3 + Pexp * 9 + Qexp + Rexp * 6},
  {-855168, FIVEexp*2 + Pexp * 6 + Qexp * 3 + Rexp * 6},
  {-1408352, FIVEexp*3 + Pexp * 3 + Qexp * 5 + Rexp * 6},
  {-2256, FIVEexp*4 + Qexp * 7 + Rexp * 6},
  {52224, FIVEexp*3 + Pexp * 7 + Qexp + Rexp * 7},
  {1037312, FIVEexp*3 + Pexp * 4 + Qexp * 3 + Rexp * 7},
  {378752, FIVEexp*3 + Pexp + Qexp * 5 + Rexp * 7},
  {-172032, FIVEexp*3 + Pexp * 5 + Qexp + Rexp * 8},
  {-36864, FIVEexp*5 + Pexp * 2 + Qexp * 3 + Rexp * 8},
  {8192, FIVEexp*5 + Pexp * 3 + Qexp + Rexp * 9},
  {4096, FIVEexp*6 + Qexp * 3 + Rexp * 9},
  {628, FIVEexp*2 + Pexp * 8 + Qexp * 8 + Sexp},
  {4821, FIVEexp*2 + Pexp * 5 + Qexp * 10 + Sexp},
  {906, FIVEexp*3 + Pexp * 2 + Qexp * 12 + Sexp},
  {-7876, FIVEexp*2 + Pexp * 9 + Qexp * 6 + Rexp + Sexp},
  {-71077, FIVEexp*2 + Pexp * 6 + Qexp * 8 + Rexp + Sexp},
  {-122499, FIVEexp*2 + Pexp * 3 + Qexp * 10 + Rexp + Sexp},
  {-33228, FIVEexp*3 + Qexp * 12 + Rexp + Sexp},
  {34317, FIVEexp*2 + Pexp * 10 + Qexp * 4 + Rexp * 2 + Sexp},
  {422511, FIVEexp*2 + Pexp * 7 + Qexp * 6 + Rexp * 2 + Sexp},
  {55786, FIVEexp*4 + Pexp * 4 + Qexp * 8 + Rexp * 2 + Sexp},
  {587894, FIVEexp*3 + Pexp + Qexp * 10 + Rexp * 2 + Sexp},
  {-53352, FIVEexp*2 + Pexp * 11 + Qexp * 2 + Rexp * 3 + Sexp},
  {-233701, FIVEexp*3 + Pexp * 8 + Qexp * 4 + Rexp * 3 + Sexp},
  {-6749187, FIVEexp*2 + Pexp * 5 + Qexp * 6 + Rexp * 3 + Sexp},
  {-3417846, FIVEexp*3 + Pexp * 2 + Qexp * 8 + Rexp * 3 + Sexp},
  {864, FIVEexp*3 + Pexp * 12 + Rexp * 4 + Sexp},
  {1215368, FIVEexp*2 + Pexp * 9 + Qexp * 2 + Rexp * 4 + Sexp},
  {12981404, FIVEexp*2 + Pexp * 6 + Qexp * 4 + Rexp * 4 + Sexp},
  {7621334, FIVEexp*3 + Pexp * 3 + Qexp * 6 + Rexp * 4 + Sexp},
  {-60923, FIVEexp*4 + Qexp * 8 + Rexp * 4 + Sexp},
  {-34368, FIVEexp*3 + Pexp * 10 + Rexp * 5 + Sexp},
  {-8544256, FIVEexp*2 + Pexp * 7 + Qexp * 2 + Rexp * 5 + Sexp},
  {-1347296, FIVEexp*4 + Pexp * 4 + Qexp * 4 + Rexp * 5 + Sexp},
  {-152456, FIVEexp*4 + Pexp + Qexp * 6 + Rexp * 5 + Sexp},
  {489472, FIVEexp*3 + Pexp * 8 + Rexp * 6 + Sexp},
  {908032, FIVEexp*4 + Pexp * 5 + Qexp * 2 + Rexp * 6 + Sexp},
  {35008, FIVEexp*6 + Pexp * 2 + Qexp * 4 + Rexp * 6 + Sexp},
  {-3127296, FIVEexp*3 + Pexp * 6 + Rexp * 7 + Sexp},
  {-260096, FIVEexp*5 + Pexp * 3 + Qexp * 2 + Rexp * 7 + Sexp},
  {-59136, FIVEexp*6 + Qexp * 4 + Rexp * 7 + Sexp},
  {73728, FIVEexp*6 + Pexp * 4 + Rexp * 8 + Sexp},
  {2048, FIVEexp*8 + Pexp + Qexp * 2 + Rexp * 8 + Sexp},
  {-16384, FIVEexp*7 + Pexp * 2 + Rexp * 9 + Sexp},
  {1134, FIVEexp*3 + Pexp * 10 + Qexp * 5 + Sexp * 2},
  {-252, FIVEexp*3 + Pexp * 7 + Qexp * 7 + Sexp * 2},
  {-3624, FIVEexp*5 + Pexp * 4 + Qexp * 9 + Sexp * 2},
  {-2028, FIVEexp*6 + Pexp + Qexp * 11 + Sexp * 2},
  {-51759, FIVEexp*2 + Pexp * 11 + Qexp * 3 + Rexp + Sexp * 2},
  {-192152, FIVEexp*2 + Pexp * 8 + Qexp * 5 + Rexp + Sexp * 2},
  {571186, FIVEexp*3 + Pexp * 5 + Qexp * 7 + Rexp + Sexp * 2},
  {14568, FIVEexp*6 + Pexp * 2 + Qexp * 9 + Rexp + Sexp * 2},
  {130248, FIVEexp*2 + Pexp * 12 + Qexp + Rexp * 2 + Sexp * 2},
  {350961, FIVEexp*3 + Pexp * 9 + Qexp * 3 + Rexp * 2 + Sexp * 2},
  {516652, FIVEexp*3 + Pexp * 6 + Qexp * 5 + Rexp * 2 + Sexp * 2},
  {17949, FIVEexp*5 + Pexp * 3 + Qexp * 7 + Rexp * 2 + Sexp * 2},
  {16654, FIVEexp*6 + Qexp * 9 + Rexp * 2 + Sexp * 2},
  {-119376, FIVEexp*4 + Pexp * 10 + Qexp + Rexp * 3 + Sexp * 2},
  {-5297492, FIVEexp*3 + Pexp * 7 + Qexp * 3 + Rexp * 3 + Sexp * 2},
  {-3180396, FIVEexp*4 + Pexp * 4 + Qexp * 5 + Rexp * 3 + Sexp * 2},
  {-259817, FIVEexp*5 + Pexp + Qexp * 7 + Rexp * 3 + Sexp * 2},
  {3770288, FIVEexp*3 + Pexp * 8 + Qexp + Rexp * 4 + Sexp * 2},
  {3369664, FIVEexp*4 + Pexp * 5 + Qexp * 3 + Rexp * 4 + Sexp * 2},
  {50732, FIVEexp*6 + Pexp * 2 + Qexp * 5 + Rexp * 4 + Sexp * 2},
  {-216192, FIVEexp*4 + Pexp * 6 + Qexp + Rexp * 5 + Sexp * 2},
  {158656, FIVEexp*6 + Pexp * 3 + Qexp * 3 + Rexp * 5 + Sexp * 2},
  {335504, FIVEexp*6 + Qexp * 5 + Rexp * 5 + Sexp * 2},
  {-16384, FIVEexp*8 + Pexp * 4 + Qexp + Rexp * 6 + Sexp * 2},
  {-110336, FIVEexp*7 + Pexp + Qexp * 3 + Rexp * 6 + Sexp * 2},
  {169984, FIVEexp*7 + Pexp * 2 + Qexp + Rexp * 7 + Sexp * 2},
  {4096, FIVEexp*8 + Qexp + Rexp * 8 + Sexp * 2},
  {2187, FIVEexp*3 + Pexp * 12 + Qexp * 2 + Sexp * 3},
  {-4356, FIVEexp*5 + Pexp * 9 + Qexp * 4 + Sexp * 3},
  {-11344, FIVEexp*6 + Pexp * 6 + Qexp * 6 + Sexp * 3},
  {-6541, FIVEexp*7 + Pexp * 3 + Qexp * 8 + Sexp * 3},
  {-4108, FIVEexp*7 + Qexp * 10 + Sexp * 3},
  {-110808, FIVEexp*2 + Pexp * 13 + Rexp + Sexp * 3},
  {100764, FIVEexp*3 + Pexp * 10 + Qexp * 2 + Rexp + Sexp * 3},
  {174064, FIVEexp*5 + Pexp * 7 + Qexp * 4 + Rexp + Sexp * 3},
  {103186, FIVEexp*6 + Pexp * 4 + Qexp * 6 + Rexp + Sexp * 3},
  {12392, FIVEexp*7 + Pexp + Qexp * 8 + Rexp + Sexp * 3},
  {616248, FIVEexp*3 + Pexp * 11 + Rexp * 2 + Sexp * 3},
  {23886, FIVEexp*6 + Pexp * 8 + Qexp * 2 + Rexp * 2 + Sexp * 3},
  {23549, FIVEexp*7 + Pexp * 5 + Qexp * 4 + Rexp * 2 + Sexp * 3},
  {23277, FIVEexp*7 + Pexp * 2 + Qexp * 6 + Rexp * 2 + Sexp * 3},
  {-1242792, FIVEexp*4 + Pexp * 9 + Rexp * 3 + Sexp * 3},
  {-2195544, FIVEexp*5 + Pexp * 6 + Qexp * 2 + Rexp * 3 + Sexp * 3},
  {-256186, FIVEexp*7 + Pexp * 3 + Qexp * 4 + Rexp * 3 + Sexp * 3},
  {-175972, FIVEexp*7 + Qexp * 6 + Rexp * 3 + Sexp * 3},
  {241152, FIVEexp*6 + Pexp * 7 + Rexp * 4 + Sexp * 3},
  {452672, FIVEexp*7 + Pexp * 4 + Qexp * 2 + Rexp * 4 + Sexp * 3},
  {440856, FIVEexp*7 + Pexp + Qexp * 4 + Rexp * 4 + Sexp * 3},
  {-616192, FIVEexp*6 + Pexp * 5 + Rexp * 5 + Sexp * 3},
  {-809344, FIVEexp*7 + Pexp * 2 + Qexp * 2 + Rexp * 5 + Sexp * 3},
  {34816, FIVEexp*8 + Pexp * 3 + Rexp * 6 + Sexp * 3},
  {-1536, FIVEexp*10 + Qexp * 2 + Rexp * 6 + Sexp * 3},
  {-26624, FIVEexp*8 + Pexp + Rexp * 7 + Sexp * 3},
  {-2916, FIVEexp*6 + Pexp * 11 + Qexp + Sexp * 4},
  {-6732, FIVEexp*7 + Pexp * 8 + Qexp * 3 + Sexp * 4},
  {-17462, FIVEexp*7 + Pexp * 5 + Qexp * 5 + Sexp * 4},
  {-708, FIVEexp*9 + Pexp * 2 + Qexp * 7 + Sexp * 4},
  {36612, FIVEexp*6 + Pexp * 9 + Qexp + Rexp + Sexp * 4},
  {31661, FIVEexp*7 + Pexp * 6 + Qexp * 3 + Rexp + Sexp * 4},
  {33772, FIVEexp*8 + Pexp * 3 + Qexp * 5 + Rexp + Sexp * 4},
  {6504, FIVEexp*9 + Qexp * 7 + Rexp + Sexp * 4},
  {-5778, FIVEexp*7 + Pexp * 7 + Qexp + Rexp * 2 + Sexp * 4},
  {-46474, FIVEexp*8 + Pexp * 4 + Qexp * 3 + Rexp * 2 + Sexp * 4},
  {-18897, FIVEexp*9 + Pexp + Qexp * 5 + Rexp * 2 + Sexp * 4},
  {-116088, FIVEexp*7 + Pexp * 5 + Qexp + Rexp * 3 + Sexp * 4},
  {204692, FIVEexp*8 + Pexp * 2 + Qexp * 3 + Rexp * 3 + Sexp * 4},
  {14112, FIVEexp*8 + Pexp * 3 + Qexp + Rexp * 4 + Sexp * 4},
  {25936, FIVEexp*9 + Qexp * 3 + Rexp * 4 + Sexp * 4},
  {38272, FIVEexp*9 + Pexp + Qexp + Rexp * 5 + Sexp * 4},
  {729, FIVEexp*7 + Pexp * 10 + Sexp * 5},
  {3537, FIVEexp*8 + Pexp * 7 + Qexp * 2 + Sexp * 5},
  {-16, FIVEexp*11 + Pexp * 4 + Qexp * 4 + Sexp * 5},
  {18, FIVEexp*11 + Pexp + Qexp * 6 + Sexp * 5},
  {-33993, FIVEexp*7 + Pexp * 8 + Rexp + Sexp * 5},
  {-8251, FIVEexp*8 + Pexp * 5 + Qexp * 2 + Rexp + Sexp * 5},
  {-3596, FIVEexp*10 + Pexp * 2 + Qexp * 4 + Rexp + Sexp * 5},
  {68028, FIVEexp*8 + Pexp * 6 + Rexp * 2 + Sexp * 5},
  {1514, FIVEexp*10 + Pexp * 3 + Qexp * 2 + Rexp * 2 + Sexp * 5},
  {-1066, FIVEexp*11 + Qexp * 4 + Rexp * 2 + Sexp * 5},
  {-52768, FIVEexp*9 + Pexp * 4 + Rexp * 3 + Sexp * 5},
  {-28856, FIVEexp*10 + Pexp + Qexp * 2 + Rexp * 3 + Sexp * 5},
  {3008, FIVEexp*11 + Pexp * 2 + Rexp * 4 + Sexp * 5},
  {-768, FIVEexp*11 + Rexp * 5 + Sexp * 5},
  {-36, FIVEexp*12 + Pexp * 6 + Qexp + Sexp * 6},
  {-16, FIVEexp*12 + Pexp * 3 + Qexp * 3 + Sexp * 6},
  {6, FIVEexp*12 + Qexp * 5 + Sexp * 6},
  {2108, FIVEexp*11 + Pexp * 4 + Qexp + Rexp + Sexp * 6},
  {1217, FIVEexp*12 + Pexp + Qexp * 3 + Rexp + Sexp * 6},
  {-892, FIVEexp*12 + Pexp * 2 + Qexp + Rexp * 2 + Sexp * 6},
  {688, FIVEexp*12 + Qexp + Rexp * 3 + Sexp * 6},
  {9, FIVEexp*13 + Pexp * 5 + Sexp * 7},
  {81, FIVEexp*13 + Pexp * 2 + Qexp * 2 + Sexp * 7},
  {-774, FIVEexp*12 + Pexp * 3 + Rexp + Sexp * 7},
  {-228, FIVEexp*13 + Qexp * 2 + Rexp + Sexp * 7},
  {424, FIVEexp*13 + Pexp + Rexp * 2 + Sexp * 7},
  {-4, FIVEexp*16 + Pexp + Qexp + Sexp * 8},
  {1, FIVEexp*17 + Sexp * 9},
  {0, END_COEFF},

  // b[2, 1]
  {-8, FIVEexp*2 + Pexp * 7 + Qexp * 9},
  {74, FIVEexp*2 + Pexp * 4 + Qexp * 11},
  {864, FIVEexp*2 + Pexp + Qexp * 13},
  {128, FIVEexp*2 + Pexp * 8 + Qexp * 7 + Rexp},
  {-768, FIVEexp*2 + Pexp * 5 + Qexp * 9 + Rexp},
  {-12654, FIVEexp*2 + Pexp * 2 + Qexp * 11 + Rexp},
  {-762, FIVEexp*2 + Pexp * 9 + Qexp * 5 + Rexp * 2},
  {1496, FIVEexp*2 + Pexp * 6 + Qexp * 7 + Rexp * 2},
  {14074, FIVEexp*3 + Pexp * 3 + Qexp * 9 + Rexp * 2},
  {17604, FIVEexp*2 + Qexp * 11 + Rexp * 2},
  {78, FIVEexp*4 + Pexp * 10 + Qexp * 3 + Rexp * 3},
  {7608, FIVEexp*2 + Pexp * 7 + Qexp * 5 + Rexp * 3},
  {-184168, FIVEexp*2 + Pexp * 4 + Qexp * 7 + Rexp * 3},
  {-242912, FIVEexp*2 + Pexp + Qexp * 9 + Rexp * 3},
  {-1728, FIVEexp*2 + Pexp * 11 + Qexp + Rexp * 4},
  {-6676, FIVEexp*3 + Pexp * 8 + Qexp * 3 + Rexp * 4},
  {39328, FIVEexp*3 + Pexp * 5 + Qexp * 5 + Rexp * 4},
  {1117074, FIVEexp*2 + Pexp * 2 + Qexp * 7 + Rexp * 4},
  {38784, FIVEexp*2 + Pexp * 9 + Qexp + Rexp * 5},
  {98688, FIVEexp*2 + Pexp * 6 + Qexp * 3 + Rexp * 5},
  {-1815728, FIVEexp*2 + Pexp * 3 + Qexp * 5 + Rexp * 5},
  {-43196, FIVEexp*3 + Qexp * 7 + Rexp * 5},
  {-291328, FIVEexp*2 + Pexp * 7 + Qexp + Rexp * 6},
  {84288, FIVEexp*3 + Pexp * 4 + Qexp * 3 + Rexp * 6},
  {333248, FIVEexp*3 + Pexp + Qexp * 5 + Rexp * 6},
  {182272, FIVEexp*3 + Pexp * 5 + Qexp + Rexp * 7},
  {-11264, FIVEexp*5 + Pexp * 2 + Qexp * 3 + Rexp * 7},
  {-8192, FIVEexp*5 + Pexp * 3 + Qexp + Rexp * 8},
  {6144, FIVEexp*6 + Qexp * 3 + Rexp * 8},
  {-24, FIVEexp*3 + Pexp * 9 + Qexp * 6 + Sexp},
  {1616, FIVEexp*2 + Pexp * 6 + Qexp * 8 + Sexp},
  {5462, FIVEexp*2 + Pexp * 3 + Qexp * 10 + Sexp},
  {-13176, FIVEexp*3 + Qexp * 12 + Sexp},
  {324, FIVEexp*3 + Pexp * 10 + Qexp * 4 + Rexp + Sexp},
  {-6944, FIVEexp*2 + Pexp * 7 + Qexp * 6 + Rexp + Sexp},
  {-1012, FIVEexp*3 + Pexp * 4 + Qexp * 8 + Rexp + Sexp},
  {191754, FIVEexp*3 + Pexp + Qexp * 10 + Rexp + Sexp},
  {-6156, FIVEexp*2 + Pexp * 11 + Qexp * 2 + Rexp * 2 + Sexp},
  {-19446, FIVEexp*2 + Pexp * 8 + Qexp * 4 + Rexp * 2 + Sexp},
  {-164632, FIVEexp*2 + Pexp * 5 + Qexp * 6 + Rexp * 2 + Sexp},
  {-901226, FIVEexp*3 + Pexp * 2 + Qexp * 8 + Rexp * 2 + Sexp},
  {5184, FIVEexp*2 + Pexp * 12 + Rexp * 3 + Sexp},
  {107334, FIVEexp*2 + Pexp * 9 + Qexp * 2 + Rexp * 3 + Sexp},
  {436266, FIVEexp*2 + Pexp * 6 + Qexp * 4 + Rexp * 3 + Sexp},
  {1498316, FIVEexp*3 + Pexp * 3 + Qexp * 6 + Rexp * 3 + Sexp},
  {70558, FIVEexp*4 + Qexp * 8 + Rexp * 3 + Sexp},
  {-175392, FIVEexp*2 + Pexp * 10 + Rexp * 4 + Sexp},
  {-1426432, FIVEexp*2 + Pexp * 7 + Qexp * 2 + Rexp * 4 + Sexp},
  {-280672, FIVEexp*4 + Pexp * 4 + Qexp * 4 + Rexp * 4 + Sexp},
  {-682462, FIVEexp*4 + Pexp + Qexp * 6 + Rexp * 4 + Sexp},
  {2434304, FIVEexp*2 + Pexp * 8 + Rexp * 5 + Sexp},
  {2795488, FIVEexp*3 + Pexp * 5 + Qexp * 2 + Rexp * 5 + Sexp},
  {288192, FIVEexp*5 + Pexp * 2 + Qexp * 4 + Rexp * 5 + Sexp},
  {-3436544, FIVEexp*3 + Pexp * 6 + Rexp * 6 + Sexp},
  {-483584, FIVEexp*5 + Pexp * 3 + Qexp * 2 + Rexp * 6 + Sexp},
  {-82304, FIVEexp*6 + Qexp * 4 + Rexp * 6 + Sexp},
  {94208, FIVEexp*6 + Pexp * 4 + Rexp * 7 + Sexp},
  {18432, FIVEexp*7 + Pexp + Qexp * 2 + Rexp * 7 + Sexp},
  {-24576, FIVEexp*7 + Pexp * 2 + Rexp * 8 + Sexp},
  {-1458, FIVEexp*2 + Pexp * 11 + Qexp * 3 + Sexp * 2},
  {-7524, FIVEexp*2 + Pexp * 8 + Qexp * 5 + Sexp * 2},
  {-44038, FIVEexp*3 + Pexp * 5 + Qexp * 7 + Sexp * 2},
  {-486, FIVEexp*7 + Pexp * 2 + Qexp * 9 + Sexp * 2},
  {10206, FIVEexp*2 + Pexp * 12 + Qexp + Rexp + Sexp * 2},
  {22032, FIVEexp*3 + Pexp * 9 + Qexp * 3 + Rexp + Sexp * 2},
  {393572, FIVEexp*3 + Pexp * 6 + Qexp * 5 + Rexp + Sexp * 2},
  {103548, FIVEexp*5 + Pexp * 3 + Qexp * 7 + Rexp + Sexp * 2},
  {-5328, FIVEexp*6 + Qexp * 9 + Rexp + Sexp * 2},
  {-3726, FIVEexp*3 + Pexp * 10 + Qexp + Rexp * 2 + Sexp * 2},
  {-255052, FIVEexp*3 + Pexp * 7 + Qexp * 3 + Rexp * 2 + Sexp * 2},
  {-664936, FIVEexp*4 + Pexp * 4 + Qexp * 5 + Rexp * 2 + Sexp * 2},
  {337528, FIVEexp*5 + Pexp + Qexp * 7 + Rexp * 2 + Sexp * 2},
  {-774588, FIVEexp*3 + Pexp * 8 + Qexp + Rexp * 3 + Sexp * 2},
  {-1122384, FIVEexp*4 + Pexp * 5 + Qexp * 3 + Rexp * 3 + Sexp * 2},
  {-189026, FIVEexp*6 + Pexp * 2 + Qexp * 5 + Rexp * 3 + Sexp * 2},
  {2327296, FIVEexp*4 + Pexp * 6 + Qexp + Rexp * 4 + Sexp * 2},
  {490912, FIVEexp*6 + Pexp * 3 + Qexp * 3 + Rexp * 4 + Sexp * 2},
  {362308, FIVEexp*6 + Qexp * 5 + Rexp * 4 + Sexp * 2},
  {-99648, FIVEexp*7 + Pexp * 4 + Qexp + Rexp * 5 + Sexp * 2},
  {-24192, FIVEexp*8 + Pexp + Qexp * 3 + Rexp * 5 + Sexp * 2},
  {7168, FIVEexp*9 + Pexp * 2 + Qexp + Rexp * 6 + Sexp * 2},
  {6144, FIVEexp*8 + Qexp + Rexp * 7 + Sexp * 2},
  {-17496, FIVEexp*2 + Pexp * 13 + Sexp * 3},
  {-81162, FIVEexp*3 + Pexp * 10 + Qexp * 2 + Sexp * 3},
  {-39012, FIVEexp*5 + Pexp * 7 + Qexp * 4 + Sexp * 3},
  {-36898, FIVEexp*6 + Pexp * 4 + Qexp * 6 + Sexp * 3},
  {-6766, FIVEexp*7 + Pexp + Qexp * 8 + Sexp * 3},
  {103518, FIVEexp*3 + Pexp * 11 + Rexp + Sexp * 3},
  {100278, FIVEexp*5 + Pexp * 8 + Qexp * 2 + Rexp + Sexp * 3},
  {138996, FIVEexp*6 + Pexp * 5 + Qexp * 4 + Rexp + Sexp * 3},
  {30486, FIVEexp*7 + Pexp * 2 + Qexp * 6 + Rexp + Sexp * 3},
  {-199422, FIVEexp*4 + Pexp * 9 + Rexp * 2 + Sexp * 3},
  {-960504, FIVEexp*5 + Pexp * 6 + Qexp * 2 + Rexp * 2 + Sexp * 3},
  {-31384, FIVEexp*8 + Pexp * 3 + Qexp * 4 + Rexp * 2 + Sexp * 3},
  {-127812, FIVEexp*7 + Qexp * 6 + Rexp * 2 + Sexp * 3},
  {24576, FIVEexp*6 + Pexp * 7 + Rexp * 3 + Sexp * 3},
  {179168, FIVEexp*7 + Pexp * 4 + Qexp * 2 + Rexp * 3 + Sexp * 3},
  {265594, FIVEexp*7 + Pexp + Qexp * 4 + Rexp * 3 + Sexp * 3},
  {-35424, FIVEexp*6 + Pexp * 5 + Rexp * 4 + Sexp * 3},
  {-535488, FIVEexp*7 + Pexp * 2 + Qexp * 2 + Rexp * 4 + Sexp * 3},
  {69376, FIVEexp*7 + Pexp * 3 + Rexp * 5 + Sexp * 3},
  {-41728, FIVEexp*8 + Qexp * 2 + Rexp * 5 + Sexp * 3},
  {-45056, FIVEexp*8 + Pexp + Rexp * 6 + Sexp * 3},
  {-486, FIVEexp*6 + Pexp * 9 + Qexp + Sexp * 4},
  {3702, FIVEexp*7 + Pexp * 6 + Qexp * 3 + Sexp * 4},
  {9194, FIVEexp*8 + Pexp * 3 + Qexp * 5 + Sexp * 4},
  {3068, FIVEexp*9 + Qexp * 7 + Sexp * 4},
  {8424, FIVEexp*7 + Pexp * 7 + Qexp + Rexp + Sexp * 4},
  {-138, FIVEexp*9 + Pexp * 4 + Qexp * 3 + Rexp + Sexp * 4},
  {-8132, FIVEexp*9 + Pexp + Qexp * 5 + Rexp + Sexp * 4},
  {-61248, FIVEexp*7 + Pexp * 5 + Qexp + Rexp * 2 + Sexp * 4},
  {139202, FIVEexp*8 + Pexp * 2 + Qexp * 3 + Rexp * 2 + Sexp * 4},
  {-14512, FIVEexp*8 + Pexp * 3 + Qexp + Rexp * 3 + Sexp * 4},
  {18364, FIVEexp*9 + Qexp * 3 + Rexp * 3 + Sexp * 4},
  {58304, FIVEexp*9 + Pexp + Qexp + Rexp * 4 + Sexp * 4},
  {-6966, FIVEexp*7 + Pexp * 8 + Sexp * 5},
  {-13842, FIVEexp*8 + Pexp * 5 + Qexp * 2 + Sexp * 5},
  {-1462, FIVEexp*10 + Pexp * 2 + Qexp * 4 + Sexp * 5},
  {13878, FIVEexp*8 + Pexp * 6 + Rexp + Sexp * 5},
  {-2554, FIVEexp*10 + Pexp * 3 + Qexp * 2 + Rexp + Sexp * 5},
  {-522, FIVEexp*11 + Qexp * 4 + Rexp + Sexp * 5},
  {-5888, FIVEexp*9 + Pexp * 4 + Rexp * 2 + Sexp * 5},
  {-23746, FIVEexp*10 + Pexp + Qexp * 2 + Rexp * 2 + Sexp * 5},
  {-128, FIVEexp*11 + Pexp * 2 + Rexp * 3 + Sexp * 5},
  {-896, FIVEexp*11 + Rexp * 4 + Sexp * 5},
  {726, FIVEexp*11 + Pexp * 4 + Qexp + Sexp * 6},
  {564, FIVEexp*12 + Pexp + Qexp * 3 + Sexp * 6},
  {142, FIVEexp*12 + Pexp * 2 + Qexp + Rexp + Sexp * 6},
  {828, FIVEexp*12 + Qexp + Rexp * 2 + Sexp * 6},
  {-138, FIVEexp*12 + Pexp * 3 + Sexp * 7},
  {-176, FIVEexp*13 + Qexp * 2 + Sexp * 7},
  {-26, FIVEexp*13 + Pexp + Rexp + Sexp * 7},
  {0, END_COEFF},

  // b[2, 2]
  {-32, FIVEexp*2 + Pexp * 5 + Qexp * 9},
  {-216, FIVEexp*2 + Pexp * 2 + Qexp * 11},
  {232, FIVEexp*2 + Pexp * 6 + Qexp * 7 + Rexp},
  {78, FIVEexp*4 + Pexp * 3 + Qexp * 9 + Rexp},
  {648, FIVEexp*2 + Qexp * 11 + Rexp},
  {-24, FIVEexp*3 + Pexp * 7 + Qexp * 5 + Rexp * 2},
  {-4334, FIVEexp*2 + Pexp * 4 + Qexp * 7 + Rexp * 2},
  {-2106, FIVEexp*3 + Pexp + Qexp * 9 + Rexp * 2},
  {-2428, FIVEexp*2 + Pexp * 8 + Qexp * 3 + Rexp * 3},
  {-618, FIVEexp*4 + Pexp * 5 + Qexp * 5 + Rexp * 3},
  {10124, FIVEexp*2 + Pexp * 2 + Qexp * 7 + Rexp * 3},
  {5112, FIVEexp*2 + Pexp * 9 + Qexp + Rexp * 4},
  {93068, FIVEexp*2 + Pexp * 6 + Qexp * 3 + Rexp * 4},
  {262622, FIVEexp*2 + Pexp * 3 + Qexp * 5 + Rexp * 4},
  {-5646, FIVEexp*3 + Qexp * 7 + Rexp * 4},
  {-116128, FIVEexp*2 + Pexp * 7 + Qexp + Rexp * 5},
  {-169744, FIVEexp*3 + Pexp * 4 + Qexp * 3 + Rexp * 5},
  {8456, FIVEexp*3 + Pexp + Qexp * 5 + Rexp * 5},
  {162944, FIVEexp*3 + Pexp * 5 + Qexp + Rexp * 6},
  {2112, FIVEexp*6 + Pexp * 2 + Qexp * 3 + Rexp * 6},
  {-13824, FIVEexp*5 + Pexp * 3 + Qexp + Rexp * 7},
  {3328, FIVEexp*6 + Qexp * 3 + Rexp * 7},
  {248, FIVEexp*2 + Pexp * 7 + Qexp * 6 + Sexp},
  {1506, FIVEexp*3 + Pexp * 4 + Qexp * 8 + Sexp},
  {7452, FIVEexp*3 + Pexp + Qexp * 10 + Sexp},
  {-2952, FIVEexp*2 + Pexp * 8 + Qexp * 4 + Rexp + Sexp},
  {-58674, FIVEexp*2 + Pexp * 5 + Qexp * 6 + Rexp + Sexp},
  {-55152, FIVEexp*3 + Pexp * 2 + Qexp * 8 + Rexp + Sexp},
  {12636, FIVEexp*2 + Pexp * 9 + Qexp * 2 + Rexp * 2 + Sexp},
  {36376, FIVEexp*3 + Pexp * 6 + Qexp * 4 + Rexp * 2 + Sexp},
  {6516, FIVEexp*5 + Pexp * 3 + Qexp * 6 + Rexp * 2 + Sexp},
  {4806, FIVEexp*5 + Qexp * 8 + Rexp * 2 + Sexp},
  {-26136, FIVEexp*2 + Pexp * 10 + Rexp * 3 + Sexp},
  {-555902, FIVEexp*2 + Pexp * 7 + Qexp * 2 + Rexp * 3 + Sexp},
  {-122812, FIVEexp*4 + Pexp * 4 + Qexp * 4 + Rexp * 3 + Sexp},
  {-198732, FIVEexp*4 + Pexp + Qexp * 6 + Rexp * 3 + Sexp},
  {742704, FIVEexp*2 + Pexp * 8 + Rexp * 4 + Sexp},
  {1407288, FIVEexp*3 + Pexp * 5 + Qexp * 2 + Rexp * 4 + Sexp},
  {81212, FIVEexp*5 + Pexp * 2 + Qexp * 4 + Rexp * 4 + Sexp},
  {-1470528, FIVEexp*3 + Pexp * 6 + Rexp * 5 + Sexp},
  {-226208, FIVEexp*5 + Pexp * 3 + Qexp * 2 + Rexp * 5 + Sexp},
  {-42608, FIVEexp*6 + Qexp * 4 + Rexp * 5 + Sexp},
  {9472, FIVEexp*7 + Pexp * 4 + Rexp * 6 + Sexp},
  {11392, FIVEexp*7 + Pexp + Qexp * 2 + Rexp * 6 + Sexp},
  {-13312, FIVEexp*7 + Pexp * 2 + Rexp * 7 + Sexp},
  {-6104, FIVEexp*3 + Pexp * 6 + Qexp * 5 + Sexp * 2},
  {-792, FIVEexp*6 + Pexp * 3 + Qexp * 7 + Sexp * 2},
  {-2592, FIVEexp*6 + Qexp * 9 + Sexp * 2},
  {2916, FIVEexp*3 + Pexp * 10 + Qexp + Rexp + Sexp * 2},
  {124296, FIVEexp*3 + Pexp * 7 + Qexp * 3 + Rexp + Sexp * 2},
  {247028, FIVEexp*4 + Pexp * 4 + Qexp * 5 + Rexp + Sexp * 2},
  {119106, FIVEexp*5 + Pexp + Qexp * 7 + Rexp + Sexp * 2},
  {-8154, FIVEexp*5 + Pexp * 8 + Qexp + Rexp * 2 + Sexp * 2},
  {-123616, FIVEexp*5 + Pexp * 5 + Qexp * 3 + Rexp * 2 + Sexp * 2},
  {-63766, FIVEexp*6 + Pexp * 2 + Qexp * 5 + Rexp * 2 + Sexp * 2},
  {735796, FIVEexp*4 + Pexp * 6 + Qexp + Rexp * 3 + Sexp * 2},
  {188412, FIVEexp*6 + Pexp * 3 + Qexp * 3 + Rexp * 3 + Sexp * 2},
  {156018, FIVEexp*6 + Qexp * 5 + Rexp * 3 + Sexp * 2},
  {-36912, FIVEexp*7 + Pexp * 4 + Qexp + Rexp * 4 + Sexp * 2},
  {-2224, FIVEexp*9 + Pexp + Qexp * 3 + Rexp * 4 + Sexp * 2},
  {70336, FIVEexp*7 + Pexp * 2 + Qexp + Rexp * 5 + Sexp * 2},
  {3328, FIVEexp*8 + Qexp + Rexp * 6 + Sexp * 2},
  {-2916, FIVEexp*3 + Pexp * 11 + Sexp * 3},
  {-4374, FIVEexp*5 + Pexp * 8 + Qexp * 2 + Sexp * 3},
  {-7258, FIVEexp*6 + Pexp * 5 + Qexp * 4 + Sexp * 3},
  {-2036, FIVEexp*7 + Pexp * 2 + Qexp * 6 + Sexp * 3},
  {22356, FIVEexp*4 + Pexp * 9 + Rexp + Sexp * 3},
  {19692, FIVEexp*5 + Pexp * 6 + Qexp * 2 + Rexp + Sexp * 3},
  {-4154, FIVEexp*8 + Pexp * 3 + Qexp * 4 + Rexp + Sexp * 3},
  {-34824, FIVEexp*7 + Qexp * 6 + Rexp + Sexp * 3},
  {-12906, FIVEexp*6 + Pexp * 7 + Rexp * 2 + Sexp * 3},
  {24956, FIVEexp*7 + Pexp * 4 + Qexp * 2 + Rexp * 2 + Sexp * 3},
  {12748, FIVEexp*8 + Pexp + Qexp * 4 + Rexp * 2 + Sexp * 3},
  {31816, FIVEexp*6 + Pexp * 5 + Rexp * 3 + Sexp * 3},
  {-142728, FIVEexp*7 + Pexp * 2 + Qexp * 2 + Rexp * 3 + Sexp * 3},
  {38176, FIVEexp*7 + Pexp * 3 + Rexp * 4 + Sexp * 3},
  {-16928, FIVEexp*8 + Qexp * 2 + Rexp * 4 + Sexp * 3},
  {-26752, FIVEexp*8 + Pexp + Rexp * 5 + Sexp * 3},
  {324, FIVEexp*8 + Pexp * 7 + Qexp + Sexp * 4},
  {3696, FIVEexp*8 + Pexp * 4 + Qexp * 3 + Sexp * 4},
  {144, FIVEexp*9 + Pexp + Qexp * 5 + Sexp * 4},
  {-21096, FIVEexp*7 + Pexp * 5 + Qexp + Rexp + Sexp * 4},
  {28854, FIVEexp*8 + Pexp * 2 + Qexp * 3 + Rexp + Sexp * 4},
  {-98, FIVEexp*11 + Pexp * 3 + Qexp + Rexp * 2 + Sexp * 4},
  {902, FIVEexp*10 + Qexp * 3 + Rexp * 2 + Sexp * 4},
  {26824, FIVEexp*9 + Pexp + Qexp + Rexp * 3 + Sexp * 4},
  {-1566, FIVEexp*8 + Pexp * 6 + Sexp * 5},
  {-1334, FIVEexp*10 + Pexp * 3 + Qexp * 2 + Sexp * 5},
  {-6, FIVEexp*12 + Qexp * 4 + Sexp * 5},
  {3324, FIVEexp*9 + Pexp * 4 + Rexp + Sexp * 5},
  {-6692, FIVEexp*10 + Pexp + Qexp * 2 + Rexp + Sexp * 5},
  {-276, FIVEexp*11 + Pexp * 2 + Rexp * 2 + Sexp * 5},
  {-496, FIVEexp*11 + Rexp * 3 + Sexp * 5},
  {172, FIVEexp*12 + Pexp * 2 + Qexp + Sexp * 6},
  {286, FIVEexp*12 + Qexp + Rexp + Sexp * 6},
  {-28, FIVEexp*13 + Pexp + Sexp * 7},
  {0, END_COEFF},

  // b[2, 3]
  {-8, FIVEexp*3 + Pexp * 6 + Qexp * 7},
  {-206, FIVEexp*2 + Pexp * 3 + Qexp * 9},
  {432, FIVEexp*2 + Qexp * 11},
  {88, FIVEexp*3 + Pexp * 7 + Qexp * 5 + Rexp},
  {2658, FIVEexp*2 + Pexp * 4 + Qexp * 7 + Rexp},
  {-5112, FIVEexp*2 + Pexp + Qexp * 9 + Rexp},
  {-66, FIVEexp*4 + Pexp * 8 + Qexp * 3 + Rexp * 2},
  {-14736, FIVEexp*2 + Pexp * 5 + Qexp * 5 + Rexp * 2},
  {8168, FIVEexp*2 + Pexp * 2 + Qexp * 7 + Rexp * 2},
  {432, FIVEexp*3 + Pexp * 9 + Qexp + Rexp * 3},
  {41638, FIVEexp*2 + Pexp * 6 + Qexp * 3 + Rexp * 3},
  {16772, FIVEexp*3 + Pexp * 3 + Qexp * 5 + Rexp * 3},
  {64, FIVEexp*5 + Qexp * 7 + Rexp * 3},
  {-1824, FIVEexp*4 + Pexp * 7 + Qexp + Rexp * 4},
  {-61528, FIVEexp*3 + Pexp * 4 + Qexp * 3 + Rexp * 4},
  {-18248, FIVEexp*3 + Pexp + Qexp * 5 + Rexp * 4},
  {58368, FIVEexp*3 + Pexp * 5 + Qexp + Rexp * 5},
  {4256, FIVEexp*5 + Pexp * 2 + Qexp * 3 + Rexp * 5},
  {-4608, FIVEexp*5 + Pexp * 3 + Qexp + Rexp * 6},
  {896, FIVEexp*6 + Qexp * 3 + Rexp * 6},
  {-72, FIVEexp*3 + Pexp * 8 + Qexp * 4 + Sexp},
  {2084, FIVEexp*2 + Pexp * 5 + Qexp * 6 + Sexp},
  {5682, FIVEexp*3 + Pexp * 2 + Qexp * 8 + Sexp},
  {108, FIVEexp*4 + Pexp * 9 + Qexp * 2 + Rexp + Sexp},
  {-10244, FIVEexp*2 + Pexp * 6 + Qexp * 4 + Rexp + Sexp},
  {-46024, FIVEexp*3 + Pexp * 3 + Qexp * 6 + Rexp + Sexp},
  {468, FIVEexp*4 + Qexp * 8 + Rexp + Sexp},
  {-1296, FIVEexp*3 + Pexp * 10 + Rexp * 2 + Sexp},
  {-57294, FIVEexp*2 + Pexp * 7 + Qexp * 2 + Rexp * 2 + Sexp},
  {8656, FIVEexp*4 + Pexp * 4 + Qexp * 4 + Rexp * 2 + Sexp},
  {-11854, FIVEexp*4 + Pexp + Qexp * 6 + Rexp * 2 + Sexp},
  {35208, FIVEexp*3 + Pexp * 8 + Rexp * 3 + Sexp},
  {38696, FIVEexp*4 + Pexp * 5 + Qexp * 2 + Rexp * 3 + Sexp},
  {266, FIVEexp*7 + Pexp * 2 + Qexp * 4 + Rexp * 3 + Sexp},
  {-344096, FIVEexp*3 + Pexp * 6 + Rexp * 4 + Sexp},
  {-46816, FIVEexp*5 + Pexp * 3 + Qexp * 2 + Rexp * 4 + Sexp},
  {-10616, FIVEexp*6 + Qexp * 4 + Rexp * 4 + Sexp},
  {11648, FIVEexp*6 + Pexp * 4 + Rexp * 5 + Sexp},
  {128, FIVEexp*9 + Pexp + Qexp * 2 + Rexp * 5 + Sexp},
  {-3584, FIVEexp*7 + Pexp * 2 + Rexp * 6 + Sexp},
  {486, FIVEexp*3 + Pexp * 10 + Qexp + Sexp * 2},
  {19314, FIVEexp*3 + Pexp * 7 + Qexp * 3 + Sexp * 2},
  {25232, FIVEexp*4 + Pexp * 4 + Qexp * 5 + Sexp * 2},
  {5064, FIVEexp*5 + Pexp + Qexp * 7 + Sexp * 2},
  {-48168, FIVEexp*3 + Pexp * 8 + Qexp + Rexp + Sexp * 2},
  {-99604, FIVEexp*4 + Pexp * 5 + Qexp * 3 + Rexp + Sexp * 2},
  {-4782, FIVEexp*6 + Pexp * 2 + Qexp * 5 + Rexp + Sexp * 2},
  {145422, FIVEexp*4 + Pexp * 6 + Qexp + Rexp * 2 + Sexp * 2},
  {30164, FIVEexp*6 + Pexp * 3 + Qexp * 3 + Rexp * 2 + Sexp * 2},
  {33656, FIVEexp*6 + Qexp * 5 + Rexp * 2 + Sexp * 2},
  {-6904, FIVEexp*7 + Pexp * 4 + Qexp + Rexp * 3 + Sexp * 2},
  {-13184, FIVEexp*7 + Pexp + Qexp * 3 + Rexp * 3 + Sexp * 2},
  {14624, FIVEexp*7 + Pexp * 2 + Qexp + Rexp * 4 + Sexp * 2},
  {896, FIVEexp*8 + Qexp + Rexp * 5 + Sexp * 2},
  {-486, FIVEexp*4 + Pexp * 9 + Sexp * 3},
  {-11502, FIVEexp*5 + Pexp * 6 + Qexp * 2 + Sexp * 3},
  {-4248, FIVEexp*7 + Pexp * 3 + Qexp * 4 + Sexp * 3},
  {-6476, FIVEexp*7 + Qexp * 6 + Sexp * 3},
  {108, FIVEexp*7 + Pexp * 7 + Rexp + Sexp * 3},
  {6794, FIVEexp*7 + Pexp * 4 + Qexp * 2 + Rexp + Sexp * 3},
  {14724, FIVEexp*7 + Pexp + Qexp * 4 + Rexp + Sexp * 3},
  {-9888, FIVEexp*6 + Pexp * 5 + Rexp * 2 + Sexp * 3},
  {-26356, FIVEexp*7 + Pexp * 2 + Qexp * 2 + Rexp * 2 + Sexp * 3},
  {2944, FIVEexp*8 + Pexp * 3 + Rexp * 3 + Sexp * 3},
  {-688, FIVEexp*9 + Qexp * 2 + Rexp * 3 + Sexp * 3},
  {-7424, FIVEexp*8 + Pexp + Rexp * 4 + Sexp * 3},
  {396, FIVEexp*7 + Pexp * 5 + Qexp + Sexp * 4},
  {2986, FIVEexp*8 + Pexp * 2 + Qexp * 3 + Sexp * 4},
  {-3832, FIVEexp*8 + Pexp * 3 + Qexp + Rexp + Sexp * 4},
  {664, FIVEexp*9 + Qexp * 3 + Rexp + Sexp * 4},
  {5448, FIVEexp*9 + Pexp + Qexp + Rexp * 2 + Sexp * 4},
  {-144, FIVEexp*9 + Pexp * 4 + Sexp * 5},
  {-998, FIVEexp*10 + Pexp + Qexp * 2 + Sexp * 5},
  {26, FIVEexp*11 + Pexp * 2 + Rexp + Sexp * 5},
  {-152, FIVEexp*11 + Rexp * 2 + Sexp * 5},
  {64, FIVEexp*12 + Qexp + Sexp * 6},
  {0, END_COEFF},

  // b[2, 4]
  {-64, FIVEexp*2 + Pexp * 4 + Qexp * 7},
  {-432, FIVEexp*2 + Pexp + Qexp * 9},
  {392, FIVEexp*2 + Pexp * 5 + Qexp * 5 + Rexp},
  {3222, FIVEexp*2 + Pexp * 2 + Qexp * 7 + Rexp},
  {-184, FIVEexp*2 + Pexp * 6 + Qexp * 3 + Rexp * 2},
  {-4508, FIVEexp*2 + Pexp * 3 + Qexp * 5 + Rexp * 2},
  {324, FIVEexp*3 + Qexp * 7 + Rexp * 2},
  {-1368, FIVEexp*2 + Pexp * 7 + Qexp + Rexp * 3},
  {-2236, FIVEexp*3 + Pexp * 4 + Qexp * 3 + Rexp * 3},
  {-5326, FIVEexp*3 + Pexp + Qexp * 5 + Rexp * 3},
  {5056, FIVEexp*3 + Pexp * 5 + Qexp + Rexp * 4},
  {1024, FIVEexp*5 + Pexp * 2 + Qexp * 3 + Rexp * 4},
  {-896, FIVEexp*5 + Pexp * 3 + Qexp + Rexp * 5},
  {192, FIVEexp*6 + Qexp * 3 + Rexp * 5},
  {-744, FIVEexp*2 + Pexp * 6 + Qexp * 4 + Sexp},
  {-414, FIVEexp*3 + Pexp * 3 + Qexp * 6 + Sexp},
  {648, FIVEexp*4 + Qexp * 8 + Sexp},
  {864, FIVEexp*2 + Pexp * 7 + Qexp * 2 + Rexp + Sexp},
  {-196, FIVEexp*4 + Pexp * 4 + Qexp * 4 + Rexp + Sexp},
  {-4626, FIVEexp*4 + Pexp + Qexp * 6 + Rexp + Sexp},
  {6264, FIVEexp*2 + Pexp * 8 + Rexp * 2 + Sexp},
  {12558, FIVEexp*3 + Pexp * 5 + Qexp * 2 + Rexp * 2 + Sexp},
  {2222, FIVEexp*5 + Pexp * 2 + Qexp * 4 + Rexp * 2 + Sexp},
  {-30192, FIVEexp*3 + Pexp * 6 + Rexp * 3 + Sexp},
  {-8672, FIVEexp*5 + Pexp * 3 + Qexp * 2 + Rexp * 3 + Sexp},
  {-1932, FIVEexp*6 + Qexp * 4 + Rexp * 3 + Sexp},
  {1792, FIVEexp*6 + Pexp * 4 + Rexp * 4 + Sexp},
  {672, FIVEexp*7 + Pexp + Qexp * 2 + Rexp * 4 + Sexp},
  {-768, FIVEexp*7 + Pexp * 2 + Rexp * 5 + Sexp},
  {-648, FIVEexp*3 + Pexp * 8 + Qexp + Sexp * 2},
  {-384, FIVEexp*4 + Pexp * 5 + Qexp * 3 + Sexp * 2},
  {12, FIVEexp*7 + Pexp * 2 + Qexp * 5 + Sexp * 2},
  {5238, FIVEexp*4 + Pexp * 6 + Qexp + Rexp + Sexp * 2},
  {1946, FIVEexp*6 + Pexp * 3 + Qexp * 3 + Rexp + Sexp * 2},
  {3564, FIVEexp*6 + Qexp * 5 + Rexp + Sexp * 2},
  {-108, FIVEexp*8 + Pexp * 4 + Qexp + Rexp * 2 + Sexp * 2},
  {-1444, FIVEexp*7 + Pexp + Qexp * 3 + Rexp * 2 + Sexp * 2},
  {1952, FIVEexp*7 + Pexp * 2 + Qexp + Rexp * 3 + Sexp * 2},
  {192, FIVEexp*8 + Qexp + Rexp * 4 + Sexp * 2},
  {-54, FIVEexp*7 + Pexp * 4 + Qexp * 2 + Sexp * 3},
  {194, FIVEexp*7 + Pexp + Qexp * 4 + Sexp * 3},
  {378, FIVEexp*6 + Pexp * 5 + Rexp + Sexp * 3},
  {-2644, FIVEexp*7 + Pexp * 2 + Qexp * 2 + Rexp + Sexp * 3},
  {1376, FIVEexp*7 + Pexp * 3 + Rexp * 2 + Sexp * 3},
  {-408, FIVEexp*8 + Qexp * 2 + Rexp * 2 + Sexp * 3},
  {-1568, FIVEexp*8 + Pexp + Rexp * 3 + Sexp * 3},
  {348, FIVEexp*8 + Pexp * 3 + Qexp + Sexp * 4},
  {24, FIVEexp*9 + Qexp * 3 + Sexp * 4},
  {602, FIVEexp*9 + Pexp + Qexp + Rexp + Sexp * 4},
  {-6, FIVEexp*11 + Pexp * 2 + Sexp * 5},
  {-28, FIVEexp*11 + Rexp + Sexp * 5},
  {0, END_COEFF},

  // b[2, 5]
  {-32, FIVEexp*2 + Pexp * 5 + Qexp * 5},
  {-216, FIVEexp*2 + Pexp * 2 + Qexp * 7},
  {48, FIVEexp*3 + Pexp * 6 + Qexp * 3 + Rexp},
  {2068, FIVEexp*2 + Pexp * 3 + Qexp * 5 + Rexp},
  {216, FIVEexp*3 + Qexp * 7 + Rexp},
  {-432, FIVEexp*2 + Pexp * 7 + Qexp + Rexp * 2},
  {-1306, FIVEexp*3 + Pexp * 4 + Qexp * 3 + Rexp * 2},
  {-2286, FIVEexp*3 + Pexp + Qexp * 5 + Rexp * 2},
  {1536, FIVEexp*3 + Pexp * 5 + Qexp + Rexp * 3},
  {64, FIVEexp*6 + Pexp * 2 + Qexp * 3 + Rexp * 3},
  {-256, FIVEexp*5 + Pexp * 3 + Qexp + Rexp * 4},
  {32, FIVEexp*6 + Qexp * 3 + Rexp * 4},
  {-432, FIVEexp*2 + Pexp * 7 + Qexp * 2 + Sexp},
  {-92, FIVEexp*4 + Pexp * 4 + Qexp * 4 + Sexp},
  {108, FIVEexp*4 + Pexp + Qexp * 6 + Sexp},
  {1296, FIVEexp*2 + Pexp * 8 + Rexp + Sexp},
  {2232, FIVEexp*3 + Pexp * 5 + Qexp * 2 + Rexp + Sexp},
  {-42, FIVEexp*5 + Pexp * 2 + Qexp * 4 + Rexp + Sexp},
  {-5832, FIVEexp*3 + Pexp * 6 + Rexp * 2 + Sexp},
  {-1312, FIVEexp*5 + Pexp * 3 + Qexp * 2 + Rexp * 2 + Sexp},
  {-342, FIVEexp*6 + Qexp * 4 + Rexp * 2 + Sexp},
  {64, FIVEexp*7 + Pexp * 4 + Rexp * 3 + Sexp},
  {128, FIVEexp*7 + Pexp + Qexp * 2 + Rexp * 3 + Sexp},
  {-128, FIVEexp*7 + Pexp * 2 + Rexp * 4 + Sexp},
  {1026, FIVEexp*4 + Pexp * 6 + Qexp + Sexp * 2},
  {372, FIVEexp*6 + Pexp * 3 + Qexp * 3 + Sexp * 2},
  {648, FIVEexp*6 + Qexp * 5 + Sexp * 2},
  {-18, FIVEexp*8 + Pexp * 4 + Qexp + Rexp + Sexp * 2},
  {-264, FIVEexp*7 + Pexp + Qexp * 3 + Rexp + Sexp * 2},
  {224, FIVEexp*7 + Pexp * 2 + Qexp + Rexp * 2 + Sexp * 2},
  {32, FIVEexp*8 + Qexp + Rexp * 3 + Sexp * 2},
  {-54, FIVEexp*6 + Pexp * 5 + Sexp * 3},
  {-248, FIVEexp*7 + Pexp * 2 + Qexp * 2 + Sexp * 3},
  {384, FIVEexp*7 + Pexp * 3 + Rexp + Sexp * 3},
  {-52, FIVEexp*8 + Qexp * 2 + Rexp + Sexp * 3},
  {-288, FIVEexp*8 + Pexp + Rexp * 2 + Sexp * 3},
  {94, FIVEexp*9 + Pexp + Qexp + Sexp * 4},
  {-6, FIVEexp*11 + Sexp * 5},
  {0, END_COEFF},

  // b[3, 0]
  {4, FIVEexp*3 + Pexp * 11 + Qexp * 6},
  {79, FIVEexp*3 + Pexp * 8 + Qexp * 8},
  {341, FIVEexp*3 + Pexp * 5 + Qexp * 10},
  {-56, FIVEexp*4 + Pexp * 2 + Qexp * 12},
  {-36, FIVEexp*3 + Pexp * 12 + Qexp * 4 + Rexp},
  {-867, FIVEexp*3 + Pexp * 9 + Qexp * 6 + Rexp},
  {-4134, FIVEexp*3 + Pexp * 6 + Qexp * 8 + Rexp},
  {8884, FIVEexp*3 + Pexp * 3 + Qexp * 10 + Rexp},
  {4368, FIVEexp*4 + Qexp * 12 + Rexp},
  {81, FIVEexp*3 + Pexp * 13 + Qexp * 2 + Rexp * 2},
  {2866, FIVEexp*3 + Pexp * 10 + Qexp * 4 + Rexp * 2},
  {15269, FIVEexp*3 + Pexp * 7 + Qexp * 6 + Rexp * 2},
  {-93954, FIVEexp*3 + Pexp * 4 + Qexp * 8 + Rexp * 2},
  {-347066, FIVEexp*3 + Pexp + Qexp * 10 + Rexp * 2},
  {-2511, FIVEexp*3 + Pexp * 11 + Qexp * 2 + Rexp * 3},
  {-16599, FIVEexp*3 + Pexp * 8 + Qexp * 4 + Rexp * 3},
  {416758, FIVEexp*3 + Pexp * 5 + Qexp * 6 + Rexp * 3},
  {423308, FIVEexp*4 + Pexp * 2 + Qexp * 8 + Rexp * 3},
  {6369, FIVEexp*3 + Pexp * 9 + Qexp * 2 + Rexp * 4},
  {-147978, FIVEexp*4 + Pexp * 6 + Qexp * 4 + Rexp * 4},
  {-1212732, FIVEexp*4 + Pexp * 3 + Qexp * 6 + Rexp * 4},
  {-46967, FIVEexp*4 + Qexp * 8 + Rexp * 4},
  {97552, FIVEexp*4 + Pexp * 7 + Qexp * 2 + Rexp * 5},
  {355988, FIVEexp*5 + Pexp * 4 + Qexp * 4 + Rexp * 5},
  {913751, FIVEexp*4 + Pexp + Qexp * 6 + Rexp * 5},
  {-1096464, FIVEexp*4 + Pexp * 5 + Qexp * 2 + Rexp * 6},
  {-652096, FIVEexp*5 + Pexp * 2 + Qexp * 4 + Rexp * 6},
  {729472, FIVEexp*5 + Pexp * 3 + Qexp * 2 + Rexp * 7},
  {54336, FIVEexp*6 + Qexp * 4 + Rexp * 7},
  {-18944, FIVEexp*7 + Pexp + Qexp * 2 + Rexp * 8},
  {108, FIVEexp*3 + Pexp * 13 + Qexp * 3 + Sexp},
  {2904, FIVEexp*3 + Pexp * 10 + Qexp * 5 + Sexp},
  {4578, FIVEexp*4 + Pexp * 7 + Qexp * 7 + Sexp},
  {2718, FIVEexp*5 + Pexp * 4 + Qexp * 9 + Sexp},
  {218, FIVEexp*7 + Pexp + Qexp * 11 + Sexp},
  {-486, FIVEexp*3 + Pexp * 14 + Qexp + Rexp + Sexp},
  {-18558, FIVEexp*3 + Pexp * 11 + Qexp * 3 + Rexp + Sexp},
  {-181394, FIVEexp*3 + Pexp * 8 + Qexp * 5 + Rexp + Sexp},
  {-23798, FIVEexp*5 + Pexp * 5 + Qexp * 7 + Rexp + Sexp},
  {-10917, FIVEexp*6 + Pexp * 2 + Qexp * 9 + Rexp + Sexp},
  {22086, FIVEexp*3 + Pexp * 12 + Qexp + Rexp * 2 + Sexp},
  {373752, FIVEexp*3 + Pexp * 9 + Qexp * 3 + Rexp * 2 + Sexp},
  {1306851, FIVEexp*3 + Pexp * 6 + Qexp * 5 + Rexp * 2 + Sexp},
  {576473, FIVEexp*4 + Pexp * 3 + Qexp * 7 + Rexp * 2 + Sexp},
  {-62717, FIVEexp*5 + Qexp * 9 + Rexp * 2 + Sexp},
  {-298734, FIVEexp*3 + Pexp * 10 + Qexp + Rexp * 3 + Sexp},
  {-311583, FIVEexp*4 + Pexp * 7 + Qexp * 3 + Rexp * 3 + Sexp},
  {-168741, FIVEexp*4 + Pexp * 4 + Qexp * 5 + Rexp * 3 + Sexp},
  {-133064, FIVEexp*5 + Pexp + Qexp * 7 + Rexp * 3 + Sexp},
  {144288, FIVEexp*4 + Pexp * 8 + Qexp + Rexp * 4 + Sexp},
  {-1584884, FIVEexp*4 + Pexp * 5 + Qexp * 3 + Rexp * 4 + Sexp},
  {1126148, FIVEexp*5 + Pexp * 2 + Qexp * 5 + Rexp * 4 + Sexp},
  {1779552, FIVEexp*4 + Pexp * 6 + Qexp + Rexp * 5 + Sexp},
  {-57712, FIVEexp*7 + Pexp * 3 + Qexp * 3 + Rexp * 5 + Sexp},
  {-522208, FIVEexp*6 + Qexp * 5 + Rexp * 5 + Sexp},
  {-278784, FIVEexp*6 + Pexp * 4 + Qexp + Rexp * 6 + Sexp},
  {187072, FIVEexp*7 + Pexp + Qexp * 3 + Rexp * 6 + Sexp},
  {-27648, FIVEexp*7 + Pexp * 2 + Qexp + Rexp * 7 + Sexp},
  {729, FIVEexp*3 + Pexp * 15 + Sexp * 2},
  {1053, FIVEexp*5 + Pexp * 12 + Qexp * 2 + Sexp * 2},
  {11232, FIVEexp*5 + Pexp * 9 + Qexp * 4 + Sexp * 2},
  {11226, FIVEexp*6 + Pexp * 6 + Qexp * 6 + Sexp * 2},
  {8052, FIVEexp*7 + Pexp * 3 + Qexp * 8 + Sexp * 2},
  {2332, FIVEexp*8 + Qexp * 10 + Sexp * 2},
  {-45684, FIVEexp*3 + Pexp * 13 + Rexp + Sexp * 2},
  {-160677, FIVEexp*4 + Pexp * 10 + Qexp * 2 + Rexp + Sexp * 2},
  {-193518, FIVEexp*5 + Pexp * 7 + Qexp * 4 + Rexp + Sexp * 2},
  {-37822, FIVEexp*7 + Pexp * 4 + Qexp * 6 + Rexp + Sexp * 2},
  {-58721, FIVEexp*7 + Pexp + Qexp * 8 + Rexp + Sexp * 2},
  {929556, FIVEexp*3 + Pexp * 11 + Rexp * 2 + Sexp * 2},
  {2049146, FIVEexp*4 + Pexp * 8 + Qexp * 2 + Rexp * 2 + Sexp * 2},
  {2368381, FIVEexp*5 + Pexp * 5 + Qexp * 4 + Rexp * 2 + Sexp * 2},
  {743676, FIVEexp*6 + Pexp * 2 + Qexp * 6 + Rexp * 2 + Sexp * 2},
  {-1523477, FIVEexp*4 + Pexp * 9 + Rexp * 3 + Sexp * 2},
  {-2086308, FIVEexp*5 + Pexp * 6 + Qexp * 2 + Rexp * 3 + Sexp * 2},
  {-131054, FIVEexp*7 + Pexp * 3 + Qexp * 4 + Rexp * 3 + Sexp * 2},
  {383803, FIVEexp*7 + Qexp * 6 + Rexp * 3 + Sexp * 2},
  {818656, FIVEexp*5 + Pexp * 7 + Rexp * 4 + Sexp * 2},
  {41536, FIVEexp*8 + Pexp * 4 + Qexp * 2 + Rexp * 4 + Sexp * 2},
  {-33277, FIVEexp*9 + Pexp + Qexp * 4 + Rexp * 4 + Sexp * 2},
  {268944, FIVEexp*6 + Pexp * 5 + Rexp * 5 + Sexp * 2},
  {120128, FIVEexp*8 + Pexp * 2 + Qexp * 2 + Rexp * 5 + Sexp * 2},
  {-74112, FIVEexp*8 + Pexp * 3 + Rexp * 6 + Sexp * 2},
  {-512, FIVEexp*9 + Qexp * 2 + Rexp * 6 + Sexp * 2},
  {18944, FIVEexp*9 + Pexp + Rexp * 7 + Sexp * 2},
  {-3078, FIVEexp*6 + Pexp * 11 + Qexp + Sexp * 3},
  {-8619, FIVEexp*7 + Pexp * 8 + Qexp * 3 + Sexp * 3},
  {-5556, FIVEexp*8 + Pexp * 5 + Qexp * 5 + Sexp * 3},
  {-1263, FIVEexp*9 + Pexp * 2 + Qexp * 7 + Sexp * 3},
  {8289, FIVEexp*7 + Pexp * 9 + Qexp + Rexp + Sexp * 3},
  {7642, FIVEexp*7 + Pexp * 6 + Qexp * 3 + Rexp + Sexp * 3},
  {-19309, FIVEexp*8 + Pexp * 3 + Qexp * 5 + Rexp + Sexp * 3},
  {-21056, FIVEexp*9 + Qexp * 7 + Rexp + Sexp * 3},
  {-139253, FIVEexp*6 + Pexp * 7 + Qexp + Rexp * 2 + Sexp * 3},
  {-3636, FIVEexp*9 + Pexp * 4 + Qexp * 3 + Rexp * 2 + Sexp * 3},
  {257528, FIVEexp*8 + Pexp + Qexp * 5 + Rexp * 2 + Sexp * 3},
  {-115004, FIVEexp*7 + Pexp * 5 + Qexp + Rexp * 3 + Sexp * 3},
  {-61476, FIVEexp*9 + Pexp * 2 + Qexp * 3 + Rexp * 3 + Sexp * 3},
  {5872, FIVEexp*10 + Pexp * 3 + Qexp + Rexp * 4 + Sexp * 3},
  {4864, FIVEexp*9 + Qexp * 3 + Rexp * 4 + Sexp * 3},
  {-175552, FIVEexp*9 + Pexp + Qexp + Rexp * 5 + Sexp * 3},
  {1026, FIVEexp*8 + Pexp * 10 + Sexp * 4},
  {4368, FIVEexp*9 + Pexp * 7 + Qexp * 2 + Sexp * 4},
  {697, FIVEexp*11 + Pexp * 4 + Qexp * 4 + Sexp * 4},
  {7, FIVEexp*14 + Pexp + Qexp * 6 + Sexp * 4},
  {-16101, FIVEexp*8 + Pexp * 8 + Rexp + Sexp * 4},
  {-12307, FIVEexp*9 + Pexp * 5 + Qexp * 2 + Rexp + Sexp * 4},
  {-6413, FIVEexp*10 + Pexp * 2 + Qexp * 4 + Rexp + Sexp * 4},
  {44286, FIVEexp*8 + Pexp * 6 + Rexp * 2 + Sexp * 4},
  {55479, FIVEexp*9 + Pexp * 3 + Qexp * 2 + Rexp * 2 + Sexp * 4},
  {-5671, FIVEexp*10 + Qexp * 4 + Rexp * 2 + Sexp * 4},
  {54092, FIVEexp*9 + Pexp * 4 + Rexp * 3 + Sexp * 4},
  {123313, FIVEexp*10 + Pexp + Qexp * 2 + Rexp * 3 + Sexp * 4},
  {-37376, FIVEexp*10 + Pexp * 2 + Rexp * 4 + Sexp * 4},
  {3776, FIVEexp*11 + Rexp * 5 + Sexp * 4},
  {-657, FIVEexp*11 + Pexp * 6 + Qexp + Sexp * 5},
  {-402, FIVEexp*12 + Pexp * 3 + Qexp * 3 + Sexp * 5},
  {77, FIVEexp*13 + Qexp * 5 + Sexp * 5},
  {-733, FIVEexp*12 + Pexp * 4 + Qexp + Rexp + Sexp * 5},
  {-5324, FIVEexp*12 + Pexp + Qexp * 3 + Rexp + Sexp * 5},
  {6808, FIVEexp*11 + Pexp * 2 + Qexp + Rexp * 2 + Sexp * 5},
  {-4896, FIVEexp*12 + Qexp + Rexp * 3 + Sexp * 5},
  {59, FIVEexp*13 + Pexp * 5 + Sexp * 6},
  {53, FIVEexp*14 + Pexp * 2 + Qexp * 2 + Sexp * 6},
  {559, FIVEexp*13 + Pexp * 3 + Rexp + Sexp * 6},
  {401, FIVEexp*14 + Qexp * 2 + Rexp + Sexp * 6},
  {-2467, FIVEexp*13 + Pexp + Rexp * 2 + Sexp * 6},
  {4, FIVEexp*17 + Pexp + Qexp + Sexp * 7},
  {-2, FIVEexp*18 + Sexp * 8},
  {0, END_COEFF},

  // b[3, 1]
  {12, FIVEexp*3 + Pexp * 9 + Qexp * 6},
  {557, FIVEexp*3 + Pexp * 6 + Qexp * 8},
  {4723, FIVEexp*3 + Pexp * 3 + Qexp * 10},
  {1656, FIVEexp*4 + Qexp * 12},
  {-108, FIVEexp*3 + Pexp * 10 + Qexp * 4 + Rexp},
  {-1217, FIVEexp*4 + Pexp * 7 + Qexp * 6 + Rexp},
  {-63236, FIVEexp*3 + Pexp * 4 + Qexp * 8 + Rexp},
  {-145354, FIVEexp*3 + Pexp + Qexp * 10 + Rexp},
  {243, FIVEexp*3 + Pexp * 11 + Qexp * 2 + Rexp * 2},
  {21029, FIVEexp*3 + Pexp * 8 + Qexp * 4 + Rexp * 2},
  {303032, FIVEexp*3 + Pexp * 5 + Qexp * 6 + Rexp * 2},
  {194188, FIVEexp*4 + Pexp * 2 + Qexp * 8 + Rexp * 2},
  {-21594, FIVEexp*3 + Pexp * 9 + Qexp * 2 + Rexp * 3},
  {-122843, FIVEexp*4 + Pexp * 6 + Qexp * 4 + Rexp * 3},
  {-645733, FIVEexp*4 + Pexp * 3 + Qexp * 6 + Rexp * 3},
  {-25237, FIVEexp*5 + Qexp * 8 + Rexp * 3},
  {97452, FIVEexp*4 + Pexp * 7 + Qexp * 2 + Rexp * 4},
  {1176466, FIVEexp*4 + Pexp * 4 + Qexp * 4 + Rexp * 4},
  {1249828, FIVEexp*4 + Pexp + Qexp * 6 + Rexp * 4},
  {-893232, FIVEexp*4 + Pexp * 5 + Qexp * 2 + Rexp * 5},
  {-688232, FIVEexp*5 + Pexp * 2 + Qexp * 4 + Rexp * 5},
  {644928, FIVEexp*5 + Pexp * 3 + Qexp * 2 + Rexp * 6},
  {75584, FIVEexp*6 + Qexp * 4 + Rexp * 6},
  {-28416, FIVEexp*7 + Pexp + Qexp * 2 + Rexp * 7},
  {324, FIVEexp*3 + Pexp * 11 + Qexp * 3 + Sexp},
  {11012, FIVEexp*3 + Pexp * 8 + Qexp * 5 + Sexp},
  {637, FIVEexp*6 + Pexp * 5 + Qexp * 7 + Sexp},
  {1, FIVEexp*10 + Pexp * 2 + Qexp * 9 + Sexp},
  {-1458, FIVEexp*3 + Pexp * 12 + Qexp + Rexp + Sexp},
  {-70872, FIVEexp*3 + Pexp * 9 + Qexp * 3 + Rexp + Sexp},
  {-662836, FIVEexp*3 + Pexp * 6 + Qexp * 5 + Rexp + Sexp},
  {-114418, FIVEexp*4 + Pexp * 3 + Qexp * 7 + Rexp + Sexp},
  {87562, FIVEexp*5 + Qexp * 9 + Rexp + Sexp},
  {81864, FIVEexp*3 + Pexp * 10 + Qexp + Rexp * 2 + Sexp},
  {287404, FIVEexp*4 + Pexp * 7 + Qexp * 3 + Rexp * 2 + Sexp},
  {40263, FIVEexp*4 + Pexp * 4 + Qexp * 5 + Rexp * 2 + Sexp},
  {-936733, FIVEexp*5 + Pexp + Qexp * 7 + Rexp * 2 + Sexp},
  {-274088, FIVEexp*4 + Pexp * 8 + Qexp + Rexp * 3 + Sexp},
  {-174326, FIVEexp*5 + Pexp * 5 + Qexp * 3 + Rexp * 3 + Sexp},
  {97068, FIVEexp*7 + Pexp * 2 + Qexp * 5 + Rexp * 3 + Sexp},
  {1823776, FIVEexp*4 + Pexp * 6 + Qexp + Rexp * 4 + Sexp},
  {-83064, FIVEexp*7 + Pexp * 3 + Qexp * 3 + Rexp * 4 + Sexp},
  {-616024, FIVEexp*6 + Qexp * 5 + Rexp * 4 + Sexp},
  {-117632, FIVEexp*6 + Pexp * 4 + Qexp + Rexp * 5 + Sexp},
  {253888, FIVEexp*7 + Pexp + Qexp * 3 + Rexp * 5 + Sexp},
  {-41472, FIVEexp*7 + Pexp * 2 + Qexp + Rexp * 6 + Sexp},
  {2187, FIVEexp*3 + Pexp * 13 + Sexp * 2},
  {15606, FIVEexp*4 + Pexp * 10 + Qexp * 2 + Sexp * 2},
  {26424, FIVEexp*5 + Pexp * 7 + Qexp * 4 + Sexp * 2},
  {2589, FIVEexp*7 + Pexp * 4 + Qexp * 6 + Sexp * 2},
  {7118, FIVEexp*7 + Pexp + Qexp * 8 + Sexp * 2},
  {-92421, FIVEexp*3 + Pexp * 11 + Rexp + Sexp * 2},
  {-185301, FIVEexp*4 + Pexp * 8 + Qexp * 2 + Rexp + Sexp * 2},
  {201789, FIVEexp*5 + Pexp * 5 + Qexp * 4 + Rexp + Sexp * 2},
  {86209, FIVEexp*6 + Pexp * 2 + Qexp * 6 + Rexp + Sexp * 2},
  {252126, FIVEexp*4 + Pexp * 9 + Rexp * 2 + Sexp * 2},
  {-220546, FIVEexp*5 + Pexp * 6 + Qexp * 2 + Rexp * 2 + Sexp * 2},
  {-11007, FIVEexp*8 + Pexp * 3 + Qexp * 4 + Rexp * 2 + Sexp * 2},
  {305871, FIVEexp*7 + Qexp * 6 + Rexp * 2 + Sexp * 2},
  {-65436, FIVEexp*6 + Pexp * 7 + Rexp * 3 + Sexp * 2},
  {85096, FIVEexp*7 + Pexp * 4 + Qexp * 2 + Rexp * 3 + Sexp * 2},
  {-135412, FIVEexp*8 + Pexp + Qexp * 4 + Rexp * 3 + Sexp * 2},
  {281712, FIVEexp*6 + Pexp * 5 + Rexp * 4 + Sexp * 2},
  {13568, FIVEexp*9 + Pexp * 2 + Qexp * 2 + Rexp * 4 + Sexp * 2},
  {-11328, FIVEexp*9 + Pexp * 3 + Rexp * 5 + Sexp * 2},
  {-768, FIVEexp*9 + Qexp * 2 + Rexp * 5 + Sexp * 2},
  {28416, FIVEexp*9 + Pexp + Rexp * 6 + Sexp * 2},
  {-351, FIVEexp*8 + Pexp * 9 + Qexp + Sexp * 3},
  {-25036, FIVEexp*7 + Pexp * 6 + Qexp * 3 + Sexp * 3},
  {-17383, FIVEexp*8 + Pexp * 3 + Qexp * 5 + Sexp * 3},
  {-8702, FIVEexp*9 + Qexp * 7 + Sexp * 3},
  {137358, FIVEexp*6 + Pexp * 7 + Qexp + Rexp + Sexp * 3},
  {3364, FIVEexp*9 + Pexp * 4 + Qexp * 3 + Rexp + Sexp * 3},
  {102192, FIVEexp*8 + Pexp + Qexp * 5 + Rexp + Sexp * 3},
  {-98218, FIVEexp*7 + Pexp * 5 + Qexp + Rexp * 2 + Sexp * 3},
  {-26624, FIVEexp*9 + Pexp * 2 + Qexp * 3 + Rexp * 2 + Sexp * 3},
  {26008, FIVEexp*9 + Pexp * 3 + Qexp + Rexp * 3 + Sexp * 3},
  {384, FIVEexp*11 + Qexp * 3 + Rexp * 3 + Sexp * 3},
  {-204736, FIVEexp*9 + Pexp + Qexp + Rexp * 4 + Sexp * 3},
  {1998, FIVEexp*8 + Pexp * 8 + Sexp * 4},
  {7401, FIVEexp*9 + Pexp * 5 + Qexp * 2 + Sexp * 4},
  {1029, FIVEexp*10 + Pexp * 2 + Qexp * 4 + Sexp * 4},
  {-38691, FIVEexp*8 + Pexp * 6 + Rexp + Sexp * 4},
  {20121, FIVEexp*9 + Pexp * 3 + Qexp * 2 + Rexp + Sexp * 4},
  {-5339, FIVEexp*10 + Qexp * 4 + Rexp + Sexp * 4},
  {23534, FIVEexp*9 + Pexp * 4 + Rexp * 2 + Sexp * 4},
  {93636, FIVEexp*10 + Pexp + Qexp * 2 + Rexp * 2 + Sexp * 4},
  {8, FIVEexp*12 + Pexp * 2 + Rexp * 3 + Sexp * 4},
  {6848, FIVEexp*11 + Rexp * 4 + Sexp * 4},
  {-611, FIVEexp*12 + Pexp * 4 + Qexp + Sexp * 5},
  {-1883, FIVEexp*12 + Pexp + Qexp * 3 + Sexp * 5},
  {-6668, FIVEexp*11 + Pexp * 2 + Qexp + Rexp + Sexp * 5},
  {-5992, FIVEexp*12 + Qexp + Rexp * 2 + Sexp * 5},
  {243, FIVEexp*13 + Pexp * 3 + Sexp * 6},
  {217, FIVEexp*14 + Qexp * 2 + Sexp * 6},
  {252, FIVEexp*13 + Pexp + Rexp + Sexp * 6},
  {0, END_COEFF},

  // b[3, 2]
  {-166, FIVEexp*3 + Pexp * 7 + Qexp * 6},
  {-2321, FIVEexp*3 + Pexp * 4 + Qexp * 8},
  {-7944, FIVEexp*3 + Pexp + Qexp * 10},
  {1169, FIVEexp*3 + Pexp * 8 + Qexp * 4 + Rexp},
  {21772, FIVEexp*3 + Pexp * 5 + Qexp * 6 + Rexp},
  {18934, FIVEexp*4 + Pexp * 2 + Qexp * 8 + Rexp},
  {-1899, FIVEexp*3 + Pexp * 9 + Qexp * 2 + Rexp * 2},
  {-13068, FIVEexp*4 + Pexp * 6 + Qexp * 4 + Rexp * 2},
  {-87369, FIVEexp*4 + Pexp * 3 + Qexp * 6 + Rexp * 2},
  {-38084, FIVEexp*4 + Qexp * 8 + Rexp * 2},
  {14284, FIVEexp*4 + Pexp * 7 + Qexp * 2 + Rexp * 3},
  {41979, FIVEexp*5 + Pexp * 4 + Qexp * 4 + Rexp * 3},
  {407512, FIVEexp*4 + Pexp + Qexp * 6 + Rexp * 3},
  {-194498, FIVEexp*4 + Pexp * 5 + Qexp * 2 + Rexp * 4},
  {-232961, FIVEexp*5 + Pexp * 2 + Qexp * 4 + Rexp * 4},
  {224816, FIVEexp*5 + Pexp * 3 + Qexp * 2 + Rexp * 5},
  {38228, FIVEexp*6 + Qexp * 4 + Rexp * 5},
  {-15392, FIVEexp*7 + Pexp + Qexp * 2 + Rexp * 6},
  {-1557, FIVEexp*3 + Pexp * 9 + Qexp * 3 + Sexp},
  {-12551, FIVEexp*3 + Pexp * 6 + Qexp * 5 + Sexp},
  {15497, FIVEexp*4 + Pexp * 3 + Qexp * 7 + Sexp},
  {23892, FIVEexp*5 + Qexp * 9 + Sexp},
  {2619, FIVEexp*3 + Pexp * 10 + Qexp + Rexp + Sexp},
  {2048, FIVEexp*4 + Pexp * 7 + Qexp * 3 + Rexp + Sexp},
  {-197926, FIVEexp*4 + Pexp * 4 + Qexp * 5 + Rexp + Sexp},
  {-272039, FIVEexp*5 + Pexp + Qexp * 7 + Rexp + Sexp},
  {-11898, FIVEexp*4 + Pexp * 8 + Qexp + Rexp * 2 + Sexp},
  {263712, FIVEexp*4 + Pexp * 5 + Qexp * 3 + Rexp * 2 + Sexp},
  {747731, FIVEexp*5 + Pexp * 2 + Qexp * 5 + Rexp * 2 + Sexp},
  {51524, FIVEexp*4 + Pexp * 6 + Qexp + Rexp * 3 + Sexp},
  {-31101, FIVEexp*7 + Pexp * 3 + Qexp * 3 + Rexp * 3 + Sexp},
  {-276391, FIVEexp*6 + Qexp * 5 + Rexp * 3 + Sexp},
  {9472, FIVEexp*6 + Pexp * 4 + Qexp + Rexp * 4 + Sexp},
  {126756, FIVEexp*7 + Pexp + Qexp * 3 + Rexp * 4 + Sexp},
  {-22464, FIVEexp*7 + Pexp * 2 + Qexp + Rexp * 5 + Sexp},
  {9234, FIVEexp*3 + Pexp * 11 + Sexp * 2},
  {58914, FIVEexp*4 + Pexp * 8 + Qexp * 2 + Sexp * 2},
  {119304, FIVEexp*5 + Pexp * 5 + Qexp * 4 + Sexp * 2},
  {74939, FIVEexp*6 + Pexp * 2 + Qexp * 6 + Sexp * 2},
  {-62262, FIVEexp*4 + Pexp * 9 + Rexp + Sexp * 2},
  {-255063, FIVEexp*5 + Pexp * 6 + Qexp * 2 + Rexp + Sexp * 2},
  {-36467, FIVEexp*7 + Pexp * 3 + Qexp * 4 + Rexp + Sexp * 2},
  {97938, FIVEexp*7 + Qexp * 6 + Rexp + Sexp * 2},
  {132822, FIVEexp*5 + Pexp * 7 + Rexp * 2 + Sexp * 2},
  {40337, FIVEexp*7 + Pexp * 4 + Qexp * 2 + Rexp * 2 + Sexp * 2},
  {-45304, FIVEexp*8 + Pexp + Qexp * 4 + Rexp * 2 + Sexp * 2},
  {-46462, FIVEexp*6 + Pexp * 5 + Rexp * 3 + Sexp * 2},
  {13556, FIVEexp*8 + Pexp * 2 + Qexp * 2 + Rexp * 3 + Sexp * 2},
  {-21968, FIVEexp*8 + Pexp * 3 + Rexp * 4 + Sexp * 2},
  {-416, FIVEexp*9 + Qexp * 2 + Rexp * 4 + Sexp * 2},
  {15392, FIVEexp*9 + Pexp + Rexp * 5 + Sexp * 2},
  {-7047, FIVEexp*6 + Pexp * 7 + Qexp + Sexp * 3},
  {-1012, FIVEexp*9 + Pexp * 4 + Qexp * 3 + Sexp * 3},
  {-16203, FIVEexp*8 + Pexp + Qexp * 5 + Sexp * 3},
  {28881, FIVEexp*7 + Pexp * 5 + Qexp + Rexp + Sexp * 3},
  {876, FIVEexp*10 + Pexp * 2 + Qexp * 3 + Rexp + Sexp * 3},
  {6629, FIVEexp*9 + Pexp * 3 + Qexp + Rexp * 2 + Sexp * 3},
  {4088, FIVEexp*9 + Qexp * 3 + Rexp * 2 + Sexp * 3},
  {-85524, FIVEexp*9 + Pexp + Qexp + Rexp * 3 + Sexp * 3},
  {3024, FIVEexp*8 + Pexp * 6 + Sexp * 4},
  {9151, FIVEexp*9 + Pexp * 3 + Qexp * 2 + Sexp * 4},
  {-2094, FIVEexp*10 + Qexp * 4 + Sexp * 4},
  {-8328, FIVEexp*9 + Pexp * 4 + Rexp + Sexp * 4},
  {26713, FIVEexp*10 + Pexp + Qexp * 2 + Rexp + Sexp * 4},
  {2153, FIVEexp*10 + Pexp * 2 + Rexp * 2 + Sexp * 4},
  {4252, FIVEexp*11 + Rexp * 3 + Sexp * 4},
  {-4243, FIVEexp*11 + Pexp * 2 + Qexp + Sexp * 5},
  {-2521, FIVEexp*12 + Qexp + Rexp + Sexp * 5},
  {247, FIVEexp*13 + Pexp + Sexp * 6},
  {0, END_COEFF},

  // b[3, 3]
  {17, FIVEexp*5 + Pexp * 5 + Qexp * 6},
  {136, FIVEexp*5 + Pexp * 2 + Qexp * 8},
  {-631, FIVEexp*4 + Pexp * 6 + Qexp * 4 + Rexp},
  {-6883, FIVEexp*4 + Pexp * 3 + Qexp * 6 + Rexp},
  {-1032, FIVEexp*5 + Qexp * 8 + Rexp},
  {1362, FIVEexp*4 + Pexp * 7 + Qexp * 2 + Rexp * 2},
  {27057, FIVEexp*4 + Pexp * 4 + Qexp * 4 + Rexp * 2},
  {70736, FIVEexp*4 + Pexp + Qexp * 6 + Rexp * 2},
  {-32759, FIVEexp*4 + Pexp * 5 + Qexp * 2 + Rexp * 3},
  {-47101, FIVEexp*5 + Pexp * 2 + Qexp * 4 + Rexp * 3},
  {52472, FIVEexp*5 + Pexp * 3 + Qexp * 2 + Rexp * 4},
  {10036, FIVEexp*6 + Qexp * 4 + Rexp * 4},
  {-4144, FIVEexp*7 + Pexp + Qexp * 2 + Rexp * 5},
  {-159, FIVEexp*4 + Pexp * 7 + Qexp * 3 + Sexp},
  {-10232, FIVEexp*4 + Pexp * 4 + Qexp * 5 + Sexp},
  {-15758, FIVEexp*5 + Pexp + Qexp * 7 + Sexp},
  {-1863, FIVEexp*4 + Pexp * 8 + Qexp + Rexp + Sexp},
  {1429, FIVEexp*5 + Pexp * 5 + Qexp * 3 + Rexp + Sexp},
  {2631, FIVEexp*7 + Pexp * 2 + Qexp * 5 + Rexp + Sexp},
  {19462, FIVEexp*4 + Pexp * 6 + Qexp + Rexp * 2 + Sexp},
  {-4039, FIVEexp*7 + Pexp * 3 + Qexp * 3 + Rexp * 2 + Sexp},
  {-60573, FIVEexp*6 + Qexp * 5 + Rexp * 2 + Sexp},
  {-1504, FIVEexp*6 + Pexp * 4 + Qexp + Rexp * 3 + Sexp},
  {29612, FIVEexp*7 + Pexp + Qexp * 3 + Rexp * 3 + Sexp},
  {-6048, FIVEexp*7 + Pexp * 2 + Qexp + Rexp * 4 + Sexp},
  {2106, FIVEexp*4 + Pexp * 9 + Sexp * 2},
  {7269, FIVEexp*5 + Pexp * 6 + Qexp * 2 + Sexp * 2},
  {2649, FIVEexp*7 + Pexp * 3 + Qexp * 4 + Sexp * 2},
  {15616, FIVEexp*7 + Qexp * 6 + Sexp * 2},
  {-1341, FIVEexp*6 + Pexp * 7 + Rexp + Sexp * 2},
  {-3554, FIVEexp*7 + Pexp * 4 + Qexp * 2 + Rexp + Sexp * 2},
  {-8492, FIVEexp*8 + Pexp + Qexp * 4 + Rexp + Sexp * 2},
  {18799, FIVEexp*6 + Pexp * 5 + Rexp * 2 + Sexp * 2},
  {692, FIVEexp*9 + Pexp * 2 + Qexp * 2 + Rexp * 2 + Sexp * 2},
  {-5832, FIVEexp*8 + Pexp * 3 + Rexp * 3 + Sexp * 2},
  {-112, FIVEexp*9 + Qexp * 2 + Rexp * 3 + Sexp * 2},
  {4144, FIVEexp*9 + Pexp + Rexp * 4 + Sexp * 2},
  {-123, FIVEexp*7 + Pexp * 5 + Qexp + Sexp * 3},
  {123, FIVEexp*9 + Pexp * 2 + Qexp * 3 + Sexp * 3},
  {1183, FIVEexp*9 + Pexp * 3 + Qexp + Rexp + Sexp * 3},
  {24, FIVEexp*11 + Qexp * 3 + Rexp + Sexp * 3},
  {-17132, FIVEexp*9 + Pexp + Qexp + Rexp * 2 + Sexp * 3},
  {1119, FIVEexp*9 + Pexp * 4 + Sexp * 4},
  {3336, FIVEexp*10 + Pexp + Qexp * 2 + Sexp * 4},
  {-7, FIVEexp*13 + Pexp * 2 + Rexp + Sexp * 4},
  {1196, FIVEexp*11 + Rexp * 2 + Sexp * 4},
  {-467, FIVEexp*12 + Qexp + Sexp * 5},
  {0, END_COEFF},

  // b[3, 4]
  {-26, FIVEexp*4 + Pexp * 6 + Qexp * 4},
  {-307, FIVEexp*4 + Pexp * 3 + Qexp * 6},
  {-792, FIVEexp*4 + Qexp * 8},
  {117, FIVEexp*4 + Pexp * 7 + Qexp * 2 + Rexp},
  {92, FIVEexp*6 + Pexp * 4 + Qexp * 4 + Rexp},
  {9386, FIVEexp*4 + Pexp + Qexp * 6 + Rexp},
  {-3269, FIVEexp*4 + Pexp * 5 + Qexp * 2 + Rexp * 2},
  {-1102, FIVEexp*6 + Pexp * 2 + Qexp * 4 + Rexp * 2},
  {6114, FIVEexp*5 + Pexp * 3 + Qexp * 2 + Rexp * 3},
  {2177, FIVEexp*6 + Qexp * 4 + Rexp * 3},
  {-888, FIVEexp*7 + Pexp + Qexp * 2 + Rexp * 4},
  {-351, FIVEexp*4 + Pexp * 8 + Qexp + Sexp},
  {-4554, FIVEexp*4 + Pexp * 5 + Qexp * 3 + Sexp},
  {-2567, FIVEexp*5 + Pexp * 2 + Qexp * 5 + Sexp},
  {5472, FIVEexp*4 + Pexp * 6 + Qexp + Rexp + Sexp},
  {-21, FIVEexp*7 + Pexp * 3 + Qexp * 3 + Rexp + Sexp},
  {-9758, FIVEexp*6 + Qexp * 5 + Rexp + Sexp},
  {196, FIVEexp*6 + Pexp * 4 + Qexp + Rexp * 2 + Sexp},
  {4879, FIVEexp*7 + Pexp + Qexp * 3 + Rexp * 2 + Sexp},
  {-1296, FIVEexp*7 + Pexp * 2 + Qexp + Rexp * 3 + Sexp},
  {891, FIVEexp*5 + Pexp * 7 + Sexp * 2},
  {557, FIVEexp*7 + Pexp * 4 + Qexp * 2 + Sexp * 2},
  {296, FIVEexp*8 + Pexp + Qexp * 4 + Sexp * 2},
  {-3081, FIVEexp*6 + Pexp * 5 + Rexp + Sexp * 2},
  {-787, FIVEexp*8 + Pexp * 2 + Qexp * 2 + Rexp + Sexp * 2},
  {-66, FIVEexp*8 + Pexp * 3 + Rexp * 2 + Sexp * 2},
  {-24, FIVEexp*9 + Qexp * 2 + Rexp * 2 + Sexp * 2},
  {888, FIVEexp*9 + Pexp + Rexp * 3 + Sexp * 2},
  {-66, FIVEexp*9 + Pexp * 3 + Qexp + Sexp * 3},
  {174, FIVEexp*9 + Qexp * 3 + Sexp * 3},
  {-2347, FIVEexp*9 + Pexp + Qexp + Rexp + Sexp * 3},
  {229, FIVEexp*10 + Pexp * 2 + Sexp * 4},
  {251, FIVEexp*11 + Rexp + Sexp * 4},
  {0, END_COEFF},

  // b[3, 5]
  {51, FIVEexp*4 + Pexp * 4 + Qexp * 4},
  {408, FIVEexp*4 + Pexp + Qexp * 6},
  {-132, FIVEexp*4 + Pexp * 5 + Qexp * 2 + Rexp},
  {-354, FIVEexp*5 + Pexp * 2 + Qexp * 4 + Rexp},
  {529, FIVEexp*5 + Pexp * 3 + Qexp * 2 + Rexp * 2},
  {332, FIVEexp*6 + Qexp * 4 + Rexp * 2},
  {-148, FIVEexp*7 + Pexp + Qexp * 2 + Rexp * 3},
  {-189, FIVEexp*4 + Pexp * 6 + Qexp + Sexp},
  {-46, FIVEexp*7 + Pexp * 3 + Qexp * 3 + Sexp},
  {-1524, FIVEexp*6 + Qexp * 5 + Sexp},
  {298, FIVEexp*6 + Pexp * 4 + Qexp + Rexp + Sexp},
  {859, FIVEexp*7 + Pexp + Qexp * 3 + Rexp + Sexp},
  {-216, FIVEexp*7 + Pexp * 2 + Qexp + Rexp * 2 + Sexp},
  {-63, FIVEexp*6 + Pexp * 5 + Sexp * 2},
  {-2, FIVEexp*10 + Pexp * 2 + Qexp * 2 + Sexp * 2},
  {-97, FIVEexp*8 + Pexp * 3 + Rexp + Sexp * 2},
  {-4, FIVEexp*9 + Qexp * 2 + Rexp + Sexp * 2},
  {148, FIVEexp*9 + Pexp + Rexp * 2 + Sexp * 2},
  {-271, FIVEexp*9 + Pexp + Qexp + Sexp * 3},
  {48, FIVEexp*11 + Sexp * 4},
  {0, END_COEFF},

  // b[4, 0]
  {24, FIVEexp*2 + Pexp * 10 + Qexp * 10},
  {554, FIVEexp*2 + Pexp * 7 + Qexp * 12},
  {4246, FIVEexp*2 + Pexp * 4 + Qexp * 14},
  {432, FIVEexp*4 + Pexp + Qexp * 16},
  {-372, FIVEexp*2 + Pexp * 11 + Qexp * 8 + Rexp},
  {-9363, FIVEexp*2 + Pexp * 8 + Qexp * 10 + Rexp},
  {-77713, FIVEexp*2 + Pexp * 5 + Qexp * 12 + Rexp},
  {-212796, FIVEexp*2 + Pexp * 2 + Qexp * 14 + Rexp},
  {2082, FIVEexp*2 + Pexp * 12 + Qexp * 6 + Rexp * 2},
  {59241, FIVEexp*2 + Pexp * 9 + Qexp * 8 + Rexp * 2},
  {543778, FIVEexp*2 + Pexp * 6 + Qexp * 10 + Rexp * 2},
  {320502, FIVEexp*3 + Pexp * 3 + Qexp * 12 + Rexp * 2},
  {-142776, FIVEexp*2 + Qexp * 14 + Rexp * 2},
  {-4887, FIVEexp*2 + Pexp * 13 + Qexp * 4 + Rexp * 3},
  {-170414, FIVEexp*2 + Pexp * 10 + Qexp * 6 + Rexp * 3},
  {-360419, FIVEexp*3 + Pexp * 7 + Qexp * 8 + Rexp * 3},
  {-5705396, FIVEexp*2 + Pexp * 4 + Qexp * 10 + Rexp * 3},
  {2167454, FIVEexp*2 + Pexp + Qexp * 12 + Rexp * 3},
  {3888, FIVEexp*2 + Pexp * 14 + Qexp * 2 + Rexp * 4},
  {211369, FIVEexp*2 + Pexp * 11 + Qexp * 4 + Rexp * 4},
  {2815581, FIVEexp*2 + Pexp * 8 + Qexp * 6 + Rexp * 4},
  {9309314, FIVEexp*2 + Pexp * 5 + Qexp * 8 + Rexp * 4},
  {-12753976, FIVEexp*2 + Pexp * 2 + Qexp * 10 + Rexp * 4},
  {-16368, FIVEexp*3 + Pexp * 12 + Qexp * 2 + Rexp * 5},
  {-350993, FIVEexp*3 + Pexp * 9 + Qexp * 4 + Rexp * 5},
  {-4296474, FIVEexp*2 + Pexp * 6 + Qexp * 6 + Rexp * 5},
  {37932428, FIVEexp*2 + Pexp * 3 + Qexp * 8 + Rexp * 5},
  {-1390543, FIVEexp*2 + Qexp * 10 + Rexp * 5},
  {236624, FIVEexp*2 + Pexp * 10 + Qexp * 2 + Rexp * 6},
  {-4635512, FIVEexp*2 + Pexp * 7 + Qexp * 4 + Rexp * 6},
  {-65981696, FIVEexp*2 + Pexp * 4 + Qexp * 6 + Rexp * 6},
  {1795751, FIVEexp*3 + Pexp + Qexp * 8 + Rexp * 6},
  {4810112, FIVEexp*2 + Pexp * 8 + Qexp * 2 + Rexp * 7},
  {14239216, FIVEexp*3 + Pexp * 5 + Qexp * 4 + Rexp * 7},
  {-18448, FIVEexp*6 + Pexp * 2 + Qexp * 6 + Rexp * 7},
  {-292864, FIVEexp*5 + Pexp * 6 + Qexp * 2 + Rexp * 8},
  {-74496, FIVEexp*6 + Pexp * 3 + Qexp * 4 + Rexp * 8},
  {-142144, FIVEexp*5 + Qexp * 6 + Rexp * 8},
  {800768, FIVEexp*5 + Pexp * 4 + Qexp * 2 + Rexp * 9},
  {126976, FIVEexp*6 + Pexp + Qexp * 4 + Rexp * 9},
  {-36864, FIVEexp*7 + Pexp * 2 + Qexp * 2 + Rexp * 10},
  {828, FIVEexp*2 + Pexp * 12 + Qexp * 7 + Sexp},
  {22059, FIVEexp*2 + Pexp * 9 + Qexp * 9 + Sexp},
  {41559, FIVEexp*3 + Pexp * 6 + Qexp * 11 + Sexp},
  {30376, FIVEexp*4 + Pexp * 3 + Qexp * 13 + Sexp},
  {216, FIVEexp*7 + Qexp * 15 + Sexp},
  {-8748, FIVEexp*2 + Pexp * 13 + Qexp * 5 + Rexp + Sexp},
  {-264259, FIVEexp*2 + Pexp * 10 + Qexp * 7 + Rexp + Sexp},
  {-2790834, FIVEexp*2 + Pexp * 7 + Qexp * 9 + Rexp + Sexp},
  {-2282604, FIVEexp*3 + Pexp * 4 + Qexp * 11 + Rexp + Sexp},
  {-467208, FIVEexp*4 + Pexp + Qexp * 13 + Rexp + Sexp},
  {27783, FIVEexp*2 + Pexp * 14 + Qexp * 3 + Rexp * 2 + Sexp},
  {209502, FIVEexp*3 + Pexp * 11 + Qexp * 5 + Rexp * 2 + Sexp},
  {13159713, FIVEexp*2 + Pexp * 8 + Qexp * 7 + Rexp * 2 + Sexp},
  {62931696, FIVEexp*2 + Pexp * 5 + Qexp * 9 + Rexp * 2 + Sexp},
  {15440351, FIVEexp*3 + Pexp * 2 + Qexp * 11 + Rexp * 2 + Sexp},
  {-23328, FIVEexp*2 + Pexp * 15 + Qexp + Rexp * 3 + Sexp},
  {-1490529, FIVEexp*2 + Pexp * 12 + Qexp * 3 + Rexp * 3 + Sexp},
  {-25543177, FIVEexp*2 + Pexp * 9 + Qexp * 5 + Rexp * 3 + Sexp},
  {-156808489, FIVEexp*2 + Pexp * 6 + Qexp * 7 + Rexp * 3 + Sexp},
  {-10123739, FIVEexp*4 + Pexp * 3 + Qexp * 9 + Rexp * 3 + Sexp},
  {705551, FIVEexp*4 + Qexp * 11 + Rexp * 3 + Sexp},
  {107568, FIVEexp*3 + Pexp * 13 + Qexp + Rexp * 4 + Sexp},
  {18493213, FIVEexp*2 + Pexp * 10 + Qexp * 3 + Rexp * 4 + Sexp},
  {180363531, FIVEexp*2 + Pexp * 7 + Qexp * 5 + Rexp * 4 + Sexp},
  {18735673, FIVEexp*4 + Pexp * 4 + Qexp * 7 + Rexp * 4 + Sexp},
  {-1145461, FIVEexp*5 + Pexp + Qexp * 9 + Rexp * 4 + Sexp},
  {-3401344, FIVEexp*2 + Pexp * 11 + Qexp + Rexp * 5 + Sexp},
  {-85472064, FIVEexp*2 + Pexp * 8 + Qexp * 3 + Rexp * 5 + Sexp},
  {-488863032, FIVEexp*2 + Pexp * 5 + Qexp * 5 + Rexp * 5 + Sexp},
  {15089671, FIVEexp*4 + Pexp * 2 + Qexp * 7 + Rexp * 5 + Sexp},
  {425728, FIVEexp*2 + Pexp * 9 + Qexp + Rexp * 6 + Sexp},
  {36524752, FIVEexp*3 + Pexp * 6 + Qexp * 3 + Rexp * 6 + Sexp},
  {-2845744, FIVEexp*4 + Pexp * 3 + Qexp * 5 + Rexp * 6 + Sexp},
  {1549496, FIVEexp*5 + Qexp * 7 + Rexp * 6 + Sexp},
  {1139712, FIVEexp*4 + Pexp * 7 + Qexp + Rexp * 7 + Sexp},
  {-1035648, FIVEexp*6 + Pexp * 4 + Qexp * 3 + Rexp * 7 + Sexp},
  {-1402752, FIVEexp*6 + Pexp + Qexp * 5 + Rexp * 7 + Sexp},
  {-237568, FIVEexp*5 + Pexp * 5 + Qexp + Rexp * 8 + Sexp},
  {397312, FIVEexp*7 + Pexp * 2 + Qexp * 3 + Rexp * 8 + Sexp},
  {16384, FIVEexp*7 + Pexp * 3 + Qexp + Rexp * 9 + Sexp},
  {12288, FIVEexp*8 + Qexp * 3 + Rexp * 9 + Sexp},
  {9234, FIVEexp*2 + Pexp * 14 + Qexp * 4 + Sexp * 2},
  {58986, FIVEexp*3 + Pexp * 11 + Qexp * 6 + Sexp * 2},
  {136073, FIVEexp*4 + Pexp * 8 + Qexp * 8 + Sexp * 2},
  {5109, FIVEexp*7 + Pexp * 5 + Qexp * 10 + Sexp * 2},
  {36162, FIVEexp*6 + Pexp * 2 + Qexp * 12 + Sexp * 2},
  {-50301, FIVEexp*2 + Pexp * 15 + Qexp * 2 + Rexp + Sexp * 2},
  {-2109159, FIVEexp*2 + Pexp * 12 + Qexp * 4 + Rexp + Sexp * 2},
  {-5947731, FIVEexp*3 + Pexp * 9 + Qexp * 6 + Rexp + Sexp * 2},
  {-6631864, FIVEexp*4 + Pexp * 6 + Qexp * 8 + Rexp + Sexp * 2},
  {-2272861, FIVEexp*5 + Pexp * 3 + Qexp * 10 + Rexp + Sexp * 2},
  {-88956, FIVEexp*6 + Qexp * 12 + Rexp + Sexp * 2},
  {34992, FIVEexp*2 + Pexp * 16 + Rexp * 2 + Sexp * 2},
  {3594051, FIVEexp*2 + Pexp * 13 + Qexp * 2 + Rexp * 2 + Sexp * 2},
  {75889471, FIVEexp*2 + Pexp * 10 + Qexp * 4 + Rexp * 2 + Sexp * 2},
  {113153304, FIVEexp*3 + Pexp * 7 + Qexp * 6 + Rexp * 2 + Sexp * 2},
  {2044347, FIVEexp*6 + Pexp * 4 + Qexp * 8 + Rexp * 2 + Sexp * 2},
  {4257638, FIVEexp*5 + Pexp + Qexp * 10 + Rexp * 2 + Sexp * 2},
  {-184032, FIVEexp*3 + Pexp * 14 + Rexp * 3 + Sexp * 2},
  {-58028619, FIVEexp*2 + Pexp * 11 + Qexp * 2 + Rexp * 3 + Sexp * 2},
  {-31083368, FIVEexp*4 + Pexp * 8 + Qexp * 4 + Rexp * 3 + Sexp * 2},
  {-113014446, FIVEexp*4 + Pexp * 5 + Qexp * 6 + Rexp * 3 + Sexp * 2},
  {-612534, FIVEexp*7 + Pexp * 2 + Qexp * 8 + Rexp * 3 + Sexp * 2},
  {8188416, FIVEexp*2 + Pexp * 12 + Rexp * 4 + Sexp * 2},
  {17400216, FIVEexp*4 + Pexp * 9 + Qexp * 2 + Rexp * 4 + Sexp * 2},
  {26758018, FIVEexp*5 + Pexp * 6 + Qexp * 4 + Rexp * 4 + Sexp * 2},
  {4015648, FIVEexp*6 + Pexp * 3 + Qexp * 6 + Rexp * 4 + Sexp * 2},
  {-253526, FIVEexp*7 + Qexp * 8 + Rexp * 4 + Sexp * 2},
  {-30283792, FIVEexp*2 + Pexp * 10 + Rexp * 5 + Sexp * 2},
  {-60348848, FIVEexp*4 + Pexp * 7 + Qexp * 2 + Rexp * 5 + Sexp * 2},
  {-1438688, FIVEexp*6 + Pexp * 4 + Qexp * 4 + Rexp * 5 + Sexp * 2},
  {1171918, FIVEexp*7 + Pexp + Qexp * 6 + Rexp * 5 + Sexp * 2},
  {3690112, FIVEexp*4 + Pexp * 8 + Rexp * 6 + Sexp * 2},
  {17772672, FIVEexp*5 + Pexp * 5 + Qexp * 2 + Rexp * 6 + Sexp * 2},
  {-1444496, FIVEexp*7 + Pexp * 2 + Qexp * 4 + Rexp * 6 + Sexp * 2},
  {-137216, FIVEexp*7 + Pexp * 6 + Rexp * 7 + Sexp * 2},
  {-828416, FIVEexp*7 + Pexp * 3 + Qexp * 2 + Rexp * 7 + Sexp * 2},
  {-153408, FIVEexp*8 + Qexp * 4 + Rexp * 7 + Sexp * 2},
  {14336, FIVEexp*9 + Pexp * 4 + Rexp * 8 + Sexp * 2},
  {14336, FIVEexp*9 + Pexp + Qexp * 2 + Rexp * 8 + Sexp * 2},
  {-12288, FIVEexp*9 + Pexp * 2 + Rexp * 9 + Sexp * 2},
  {6561, FIVEexp*3 + Pexp * 16 + Qexp + Sexp * 3},
  {58887, FIVEexp*4 + Pexp * 13 + Qexp * 3 + Sexp * 3},
  {176712, FIVEexp*5 + Pexp * 10 + Qexp * 5 + Sexp * 3},
  {214886, FIVEexp*6 + Pexp * 7 + Qexp * 7 + Sexp * 3},
  {91476, FIVEexp*7 + Pexp * 4 + Qexp * 9 + Sexp * 3},
  {1972, FIVEexp*9 + Pexp + Qexp * 11 + Sexp * 3},
  {-739206, FIVEexp*3 + Pexp * 14 + Qexp + Rexp + Sexp * 3},
  {-3760281, FIVEexp*4 + Pexp * 11 + Qexp * 3 + Rexp + Sexp * 3},
  {-1246121, FIVEexp*6 + Pexp * 8 + Qexp * 5 + Rexp + Sexp * 3},
  {-3380518, FIVEexp*6 + Pexp * 5 + Qexp * 7 + Rexp + Sexp * 3},
  {-581731, FIVEexp*7 + Pexp * 2 + Qexp * 9 + Rexp + Sexp * 3},
  {14594904, FIVEexp*3 + Pexp * 12 + Qexp + Rexp * 2 + Sexp * 3},
  {402371, FIVEexp*7 + Pexp * 9 + Qexp * 3 + Rexp * 2 + Sexp * 3},
  {45349612, FIVEexp*5 + Pexp * 6 + Qexp * 5 + Rexp * 2 + Sexp * 3},
  {2923146, FIVEexp*7 + Pexp * 3 + Qexp * 7 + Rexp * 2 + Sexp * 3},
  {443009, FIVEexp*7 + Qexp * 9 + Rexp * 2 + Sexp * 3},
  {-28146921, FIVEexp*4 + Pexp * 10 + Qexp + Rexp * 3 + Sexp * 3},
  {-2418268, FIVEexp*7 + Pexp * 7 + Qexp * 3 + Rexp * 3 + Sexp * 3},
  {-6426739, FIVEexp*7 + Pexp * 4 + Qexp * 5 + Rexp * 3 + Sexp * 3},
  {-2404471, FIVEexp*7 + Pexp + Qexp * 7 + Rexp * 3 + Sexp * 3},
  {4836976, FIVEexp*6 + Pexp * 8 + Qexp + Rexp * 4 + Sexp * 3},
  {877568, FIVEexp*8 + Pexp * 5 + Qexp * 3 + Rexp * 4 + Sexp * 3},
  {3780919, FIVEexp*7 + Pexp * 2 + Qexp * 5 + Rexp * 4 + Sexp * 3},
  {-275664, FIVEexp*8 + Pexp * 6 + Qexp + Rexp * 5 + Sexp * 3},
  {682624, FIVEexp*7 + Pexp * 3 + Qexp * 3 + Rexp * 5 + Sexp * 3},
  {3476312, FIVEexp*7 + Qexp * 5 + Rexp * 5 + Sexp * 3},
  {-338048, FIVEexp*7 + Pexp * 4 + Qexp + Rexp * 6 + Sexp * 3},
  {-481792, FIVEexp*8 + Pexp + Qexp * 3 + Rexp * 6 + Sexp * 3},
  {47104, FIVEexp*9 + Pexp * 2 + Qexp + Rexp * 7 + Sexp * 3},
  {12288, FIVEexp*10 + Qexp + Rexp * 8 + Sexp * 3},
  {15309, FIVEexp*5 + Pexp * 15 + Sexp * 4},
  {73629, FIVEexp*6 + Pexp * 12 + Qexp * 2 + Sexp * 4},
  {118137, FIVEexp*7 + Pexp * 9 + Qexp * 4 + Sexp * 4},
  {63884, FIVEexp*8 + Pexp * 6 + Qexp * 6 + Sexp * 4},
  {2353, FIVEexp*10 + Pexp * 3 + Qexp * 8 + Sexp * 4},
  {28, FIVEexp*11 + Qexp * 10 + Sexp * 4},
  {-381996, FIVEexp*5 + Pexp * 13 + Rexp + Sexp * 4},
  {-1332342, FIVEexp*6 + Pexp * 10 + Qexp * 2 + Rexp + Sexp * 4},
  {-1262592, FIVEexp*7 + Pexp * 7 + Qexp * 4 + Rexp + Sexp * 4},
  {-577964, FIVEexp*8 + Pexp * 4 + Qexp * 6 + Rexp + Sexp * 4},
  {-306, FIVEexp*12 + Pexp + Qexp * 8 + Rexp + Sexp * 4},
  {4298022, FIVEexp*5 + Pexp * 11 + Rexp * 2 + Sexp * 4},
  {2415486, FIVEexp*7 + Pexp * 8 + Qexp * 2 + Rexp * 2 + Sexp * 4},
  {11202018, FIVEexp*7 + Pexp * 5 + Qexp * 4 + Rexp * 2 + Sexp * 4},
  {1520161, FIVEexp*8 + Pexp * 2 + Qexp * 6 + Rexp * 2 + Sexp * 4},
  {-894327, FIVEexp*7 + Pexp * 9 + Rexp * 3 + Sexp * 4},
  {-10941582, FIVEexp*7 + Pexp * 6 + Qexp * 2 + Rexp * 3 + Sexp * 4},
  {-886097, FIVEexp*9 + Pexp * 3 + Qexp * 4 + Rexp * 3 + Sexp * 4},
  {-292284, FIVEexp*9 + Qexp * 6 + Rexp * 3 + Sexp * 4},
  {2074568, FIVEexp*7 + Pexp * 7 + Rexp * 4 + Sexp * 4},
  {786624, FIVEexp*9 + Pexp * 4 + Qexp * 2 + Rexp * 4 + Sexp * 4},
  {78352, FIVEexp*10 + Pexp + Qexp * 4 + Rexp * 4 + Sexp * 4},
  {-2124656, FIVEexp*7 + Pexp * 5 + Rexp * 5 + Sexp * 4},
  {-373232, FIVEexp*9 + Pexp * 2 + Qexp * 2 + Rexp * 5 + Sexp * 4},
  {62208, FIVEexp*9 + Pexp * 3 + Rexp * 6 + Sexp * 4},
  {-106432, FIVEexp*10 + Qexp * 2 + Rexp * 6 + Sexp * 4},
  {-2048, FIVEexp*11 + Pexp + Rexp * 7 + Sexp * 4},
  {-972, FIVEexp*8 + Pexp * 11 + Qexp + Sexp * 5},
  {-5943, FIVEexp*9 + Pexp * 8 + Qexp * 3 + Sexp * 5},
  {-2133, FIVEexp*10 + Pexp * 5 + Qexp * 5 + Sexp * 5},
  {-138, FIVEexp*12 + Pexp * 2 + Qexp * 7 + Sexp * 5},
  {-106218, FIVEexp*8 + Pexp * 9 + Qexp + Rexp + Sexp * 5},
  {-214556, FIVEexp*9 + Pexp * 6 + Qexp * 3 + Rexp + Sexp * 5},
  {-13262, FIVEexp*10 + Pexp * 3 + Qexp * 5 + Rexp + Sexp * 5},
  {1277, FIVEexp*12 + Qexp * 7 + Rexp + Sexp * 5},
  {1587497, FIVEexp*8 + Pexp * 7 + Qexp + Rexp * 2 + Sexp * 5},
  {23648, FIVEexp*11 + Pexp * 4 + Qexp * 3 + Rexp * 2 + Sexp * 5},
  {3733, FIVEexp*10 + Pexp + Qexp * 5 + Rexp * 2 + Sexp * 5},
  {-1160168, FIVEexp*9 + Pexp * 5 + Qexp + Rexp * 3 + Sexp * 5},
  {-6183, FIVEexp*12 + Pexp * 2 + Qexp * 3 + Rexp * 3 + Sexp * 5},
  {51824, FIVEexp*11 + Pexp * 3 + Qexp + Rexp * 4 + Sexp * 5},
  {13352, FIVEexp*12 + Qexp * 3 + Rexp * 4 + Sexp * 5},
  {384, FIVEexp*12 + Pexp + Qexp + Rexp * 5 + Sexp * 5},
  {486, FIVEexp*11 + Pexp * 10 + Sexp * 6},
  {4989, FIVEexp*11 + Pexp * 7 + Qexp * 2 + Sexp * 6},
  {1399, FIVEexp*12 + Pexp * 4 + Qexp * 4 + Sexp * 6},
  {-2, FIVEexp*14 + Pexp + Qexp * 6 + Sexp * 6},
  {-36567, FIVEexp*10 + Pexp * 8 + Rexp + Sexp * 6},
  {-6107, FIVEexp*11 + Pexp * 5 + Qexp * 2 + Rexp + Sexp * 6},
  {1666, FIVEexp*12 + Pexp * 2 + Qexp * 4 + Rexp + Sexp * 6},
  {165446, FIVEexp*10 + Pexp * 6 + Rexp * 2 + Sexp * 6},
  {2288, FIVEexp*12 + Pexp * 3 + Qexp * 2 + Rexp * 2 + Sexp * 6},
  {-11514, FIVEexp*12 + Qexp * 4 + Rexp * 2 + Sexp * 6},
  {-12128, FIVEexp*12 + Pexp * 4 + Rexp * 3 + Sexp * 6},
  {-15574, FIVEexp*12 + Pexp + Qexp * 2 + Rexp * 3 + Sexp * 6},
  {9616, FIVEexp*12 + Pexp * 2 + Rexp * 4 + Sexp * 6},
  {-2752, FIVEexp*12 + Rexp * 5 + Sexp * 6},
  {-534, FIVEexp*13 + Pexp * 6 + Qexp + Sexp * 7},
  {-239, FIVEexp*14 + Pexp * 3 + Qexp * 3 + Sexp * 7},
  {-1, FIVEexp*16 + Qexp * 5 + Sexp * 7},
  {1334, FIVEexp*13 + Pexp * 4 + Qexp + Rexp + Sexp * 7},
  {647, FIVEexp*14 + Pexp + Qexp * 3 + Rexp + Sexp * 7},
  {-751, FIVEexp*13 + Pexp * 2 + Qexp + Rexp * 2 + Sexp * 7},
  {264, FIVEexp*14 + Qexp + Rexp * 3 + Sexp * 7},
  {21, FIVEexp*15 + Pexp * 5 + Sexp * 8},
  {7, FIVEexp*16 + Pexp * 2 + Qexp * 2 + Sexp * 8},
  {-67, FIVEexp*15 + Pexp * 3 + Rexp + Sexp * 8},
  {-21, FIVEexp*16 + Qexp * 2 + Rexp + Sexp * 8},
  {57, FIVEexp*15 + Pexp + Rexp * 2 + Sexp * 8},
  {0, END_COEFF},

  // b[4, 1]
  {-24, FIVEexp*2 + Pexp * 11 + Qexp * 8},
  {-562, FIVEexp*2 + Pexp * 8 + Qexp * 10},
  {-4364, FIVEexp*2 + Pexp * 5 + Qexp * 12},
  {-11232, FIVEexp*2 + Pexp * 2 + Qexp * 14},
  {288, FIVEexp*2 + Pexp * 12 + Qexp * 6 + Rexp},
  {7548, FIVEexp*2 + Pexp * 9 + Qexp * 8 + Rexp},
  {64869, FIVEexp*2 + Pexp * 6 + Qexp * 10 + Rexp},
  {183083, FIVEexp*2 + Pexp * 3 + Qexp * 12 + Rexp},
  {216, FIVEexp*2 + Qexp * 14 + Rexp},
  {-1134, FIVEexp*2 + Pexp * 13 + Qexp * 4 + Rexp * 2},
  {-36424, FIVEexp*2 + Pexp * 10 + Qexp * 6 + Rexp * 2},
  {-369519, FIVEexp*2 + Pexp * 7 + Qexp * 8 + Rexp * 2},
  {-1228756, FIVEexp*2 + Pexp * 4 + Qexp * 10 + Rexp * 2},
  {-223038, FIVEexp*2 + Pexp + Qexp * 12 + Rexp * 2},
  {1458, FIVEexp*2 + Pexp * 14 + Qexp * 2 + Rexp * 3},
  {2957, FIVEexp*4 + Pexp * 11 + Qexp * 4 + Rexp * 3},
  {1005511, FIVEexp*2 + Pexp * 8 + Qexp * 6 + Rexp * 3},
  {4383658, FIVEexp*2 + Pexp * 5 + Qexp * 8 + Rexp * 3},
  {2825106, FIVEexp*2 + Pexp * 2 + Qexp * 10 + Rexp * 3},
  {-52686, FIVEexp*2 + Pexp * 12 + Qexp * 2 + Rexp * 4},
  {-1314284, FIVEexp*2 + Pexp * 9 + Qexp * 4 + Rexp * 4},
  {-8765023, FIVEexp*2 + Pexp * 6 + Qexp * 6 + Rexp * 4},
  {-2620527, FIVEexp*3 + Pexp * 3 + Qexp * 8 + Rexp * 4},
  {-20819, FIVEexp*4 + Qexp * 10 + Rexp * 4},
  {659366, FIVEexp*2 + Pexp * 10 + Qexp * 2 + Rexp * 5},
  {1777938, FIVEexp*3 + Pexp * 7 + Qexp * 4 + Rexp * 5},
  {1027478, FIVEexp*4 + Pexp * 4 + Qexp * 6 + Rexp * 5},
  {162022, FIVEexp*4 + Pexp + Qexp * 8 + Rexp * 5},
  {-126952, FIVEexp*4 + Pexp * 8 + Qexp * 2 + Rexp * 6},
  {-693088, FIVEexp*4 + Pexp * 5 + Qexp * 4 + Rexp * 6},
  {-29994, FIVEexp*5 + Pexp * 2 + Qexp * 6 + Rexp * 6},
  {-23776, FIVEexp*5 + Pexp * 6 + Qexp * 2 + Rexp * 7},
  {-338528, FIVEexp*5 + Pexp * 3 + Qexp * 4 + Rexp * 7},
  {-189176, FIVEexp*5 + Qexp * 6 + Rexp * 7},
  {605312, FIVEexp*5 + Pexp * 4 + Qexp * 2 + Rexp * 8},
  {178944, FIVEexp*6 + Pexp + Qexp * 4 + Rexp * 8},
  {-55296, FIVEexp*7 + Pexp * 2 + Qexp * 2 + Rexp * 9},
  {-648, FIVEexp*2 + Pexp * 13 + Qexp * 5 + Sexp},
  {-2876, FIVEexp*3 + Pexp * 10 + Qexp * 7 + Sexp},
  {-104153, FIVEexp*2 + Pexp * 7 + Qexp * 9 + Sexp},
  {-36723, FIVEexp*3 + Pexp * 4 + Qexp * 11 + Sexp},
  {19764, FIVEexp*4 + Pexp + Qexp * 13 + Sexp},
  {972, FIVEexp*3 + Pexp * 14 + Qexp * 3 + Rexp + Sexp},
  {129096, FIVEexp*2 + Pexp * 11 + Qexp * 5 + Rexp + Sexp},
  {1092069, FIVEexp*2 + Pexp * 8 + Qexp * 7 + Rexp + Sexp},
  {2379239, FIVEexp*2 + Pexp * 5 + Qexp * 9 + Rexp + Sexp},
  {-1098471, FIVEexp*3 + Pexp * 2 + Qexp * 11 + Rexp + Sexp},
  {-8748, FIVEexp*2 + Pexp * 15 + Qexp + Rexp * 2 + Sexp},
  {-356157, FIVEexp*2 + Pexp * 12 + Qexp * 3 + Rexp * 2 + Sexp},
  {-4036729, FIVEexp*2 + Pexp * 9 + Qexp * 5 + Rexp * 2 + Sexp},
  {-13011652, FIVEexp*2 + Pexp * 6 + Qexp * 7 + Rexp * 2 + Sexp},
  {585128, FIVEexp*4 + Pexp * 3 + Qexp * 9 + Rexp * 2 + Sexp},
  {150948, FIVEexp*4 + Qexp * 11 + Rexp * 2 + Sexp},
  {305316, FIVEexp*2 + Pexp * 13 + Qexp + Rexp * 3 + Sexp},
  {6519816, FIVEexp*2 + Pexp * 10 + Qexp * 3 + Rexp * 3 + Sexp},
  {38982359, FIVEexp*2 + Pexp * 7 + Qexp * 5 + Rexp * 3 + Sexp},
  {1489586, FIVEexp*4 + Pexp * 4 + Qexp * 7 + Rexp * 3 + Sexp},
  {-158518, FIVEexp*5 + Pexp + Qexp * 9 + Rexp * 3 + Sexp},
  {-3893796, FIVEexp*2 + Pexp * 11 + Qexp + Rexp * 4 + Sexp},
  {-11253914, FIVEexp*3 + Pexp * 8 + Qexp * 3 + Rexp * 4 + Sexp},
  {-8916042, FIVEexp*4 + Pexp * 5 + Qexp * 5 + Rexp * 4 + Sexp},
  {-609276, FIVEexp*5 + Pexp * 2 + Qexp * 7 + Rexp * 4 + Sexp},
  {217136, FIVEexp*5 + Pexp * 9 + Qexp + Rexp * 5 + Sexp},
  {13081944, FIVEexp*4 + Pexp * 6 + Qexp * 3 + Rexp * 5 + Sexp},
  {5786256, FIVEexp*5 + Pexp * 3 + Qexp * 5 + Rexp * 5 + Sexp},
  {347894, FIVEexp*6 + Qexp * 7 + Rexp * 5 + Sexp},
  {-953536, FIVEexp*5 + Pexp * 7 + Qexp + Rexp * 6 + Sexp},
  {-1866464, FIVEexp*6 + Pexp * 4 + Qexp * 3 + Rexp * 6 + Sexp},
  {-350944, FIVEexp*7 + Pexp + Qexp * 5 + Rexp * 6 + Sexp},
  {2010368, FIVEexp*5 + Pexp * 5 + Qexp + Rexp * 7 + Sexp},
  {623232, FIVEexp*7 + Pexp * 2 + Qexp * 3 + Rexp * 7 + Sexp},
  {-36864, FIVEexp*7 + Pexp * 3 + Qexp + Rexp * 8 + Sexp},
  {18432, FIVEexp*8 + Qexp * 3 + Rexp * 8 + Sexp},
  {-4374, FIVEexp*2 + Pexp * 15 + Qexp * 2 + Sexp * 2},
  {-96228, FIVEexp*2 + Pexp * 12 + Qexp * 4 + Sexp * 2},
  {-129002, FIVEexp*3 + Pexp * 9 + Qexp * 6 + Sexp * 2},
  {-7888, FIVEexp*4 + Pexp * 6 + Qexp * 8 + Sexp * 2},
  {64368, FIVEexp*5 + Pexp * 3 + Qexp * 10 + Sexp * 2},
  {-15552, FIVEexp*6 + Qexp * 12 + Sexp * 2},
  {13122, FIVEexp*2 + Pexp * 16 + Rexp + Sexp * 2},
  {422091, FIVEexp*2 + Pexp * 13 + Qexp * 2 + Rexp + Sexp * 2},
  {3520764, FIVEexp*2 + Pexp * 10 + Qexp * 4 + Rexp + Sexp * 2},
  {-33669, FIVEexp*3 + Pexp * 7 + Qexp * 6 + Rexp + Sexp * 2},
  {-24581, FIVEexp*7 + Pexp * 4 + Qexp * 8 + Rexp + Sexp * 2},
  {563052, FIVEexp*5 + Pexp + Qexp * 10 + Rexp + Sexp * 2},
  {-478224, FIVEexp*2 + Pexp * 14 + Rexp * 2 + Sexp * 2},
  {-7855002, FIVEexp*2 + Pexp * 11 + Qexp * 2 + Rexp * 2 + Sexp * 2},
  {-888314, FIVEexp*4 + Pexp * 8 + Qexp * 4 + Rexp * 2 + Sexp * 2},
  {6741232, FIVEexp*4 + Pexp * 5 + Qexp * 6 + Rexp * 2 + Sexp * 2},
  {-10078, FIVEexp*6 + Pexp * 2 + Qexp * 8 + Rexp * 2 + Sexp * 2},
  {6506244, FIVEexp*2 + Pexp * 12 + Rexp * 3 + Sexp * 2},
  {3009392, FIVEexp*4 + Pexp * 9 + Qexp * 2 + Rexp * 3 + Sexp * 2},
  {241179, FIVEexp*5 + Pexp * 6 + Qexp * 4 + Rexp * 3 + Sexp * 2},
  {-1627076, FIVEexp*6 + Pexp * 3 + Qexp * 6 + Rexp * 3 + Sexp * 2},
  {-181018, FIVEexp*7 + Qexp * 8 + Rexp * 3 + Sexp * 2},
  {-2003118, FIVEexp*4 + Pexp * 10 + Rexp * 4 + Sexp * 2},
  {-4007832, FIVEexp*5 + Pexp * 7 + Qexp * 2 + Rexp * 4 + Sexp * 2},
  {231266, FIVEexp*7 + Pexp * 4 + Qexp * 4 + Rexp * 4 + Sexp * 2},
  {774808, FIVEexp*7 + Pexp + Qexp * 6 + Rexp * 4 + Sexp * 2},
  {2184872, FIVEexp*5 + Pexp * 8 + Rexp * 5 + Sexp * 2},
  {3658048, FIVEexp*6 + Pexp * 5 + Qexp * 2 + Rexp * 5 + Sexp * 2},
  {-194162, FIVEexp*8 + Pexp * 2 + Qexp * 4 + Rexp * 5 + Sexp * 2},
  {-312992, FIVEexp*7 + Pexp * 6 + Rexp * 6 + Sexp * 2},
  {-371072, FIVEexp*8 + Pexp * 3 + Qexp * 2 + Rexp * 6 + Sexp * 2},
  {-8408, FIVEexp*10 + Qexp * 4 + Rexp * 6 + Sexp * 2},
  {120192, FIVEexp*8 + Pexp * 4 + Rexp * 7 + Sexp * 2},
  {6144, FIVEexp*10 + Pexp + Qexp * 2 + Rexp * 7 + Sexp * 2},
  {-18432, FIVEexp*9 + Pexp * 2 + Rexp * 8 + Sexp * 2},
  {12393, FIVEexp*3 + Pexp * 14 + Qexp + Sexp * 3},
  {82998, FIVEexp*4 + Pexp * 11 + Qexp * 3 + Sexp * 3},
  {38386, FIVEexp*6 + Pexp * 8 + Qexp * 5 + Sexp * 3},
  {154954, FIVEexp*6 + Pexp * 5 + Qexp * 7 + Sexp * 3},
  {-21677, FIVEexp*7 + Pexp * 2 + Qexp * 9 + Sexp * 3},
  {-839079, FIVEexp*3 + Pexp * 12 + Qexp + Rexp + Sexp * 3},
  {-123996, FIVEexp*6 + Pexp * 9 + Qexp * 3 + Rexp + Sexp * 3},
  {-3667537, FIVEexp*5 + Pexp * 6 + Qexp * 5 + Rexp + Sexp * 3},
  {131832, FIVEexp*7 + Pexp * 3 + Qexp * 7 + Rexp + Sexp * 3},
  {135056, FIVEexp*7 + Qexp * 9 + Rexp + Sexp * 3},
  {2138022, FIVEexp*4 + Pexp * 10 + Qexp + Rexp * 2 + Sexp * 3},
  {1519184, FIVEexp*6 + Pexp * 7 + Qexp * 3 + Rexp * 2 + Sexp * 3},
  {29254, FIVEexp*9 + Pexp * 4 + Qexp * 5 + Rexp * 2 + Sexp * 3},
  {-106098, FIVEexp*7 + Pexp + Qexp * 7 + Rexp * 2 + Sexp * 3},
  {-698062, FIVEexp*6 + Pexp * 8 + Qexp + Rexp * 3 + Sexp * 3},
  {-2055872, FIVEexp*7 + Pexp * 5 + Qexp * 3 + Rexp * 3 + Sexp * 3},
  {-1422524, FIVEexp*7 + Pexp * 2 + Qexp * 5 + Rexp * 3 + Sexp * 3},
  {598312, FIVEexp*7 + Pexp * 6 + Qexp + Rexp * 4 + Sexp * 3},
  {1460504, FIVEexp*8 + Pexp * 3 + Qexp * 3 + Rexp * 4 + Sexp * 3},
  {140718, FIVEexp*9 + Qexp * 5 + Rexp * 4 + Sexp * 3},
  {-37536, FIVEexp*9 + Pexp * 4 + Qexp + Rexp * 5 + Sexp * 3},
  {-74464, FIVEexp*9 + Pexp + Qexp * 3 + Rexp * 5 + Sexp * 3},
  {896, FIVEexp*10 + Pexp * 2 + Qexp + Rexp * 6 + Sexp * 3},
  {18432, FIVEexp*10 + Qexp + Rexp * 7 + Sexp * 3},
  {5103, FIVEexp*5 + Pexp * 13 + Sexp * 4},
  {81, FIVEexp*6 + Pexp * 10 + Qexp * 2 + Sexp * 4},
  {-42014, FIVEexp*7 + Pexp * 7 + Qexp * 4 + Sexp * 4},
  {-97898, FIVEexp*8 + Pexp * 4 + Qexp * 6 + Sexp * 4},
  {-822, FIVEexp*11 + Pexp + Qexp * 8 + Sexp * 4},
  {-36207, FIVEexp*5 + Pexp * 11 + Rexp + Sexp * 4},
  {-1179, FIVEexp*9 + Pexp * 8 + Qexp * 2 + Rexp + Sexp * 4},
  {876162, FIVEexp*7 + Pexp * 5 + Qexp * 4 + Rexp + Sexp * 4},
  {454084, FIVEexp*8 + Pexp * 2 + Qexp * 6 + Rexp + Sexp * 4},
  {324, FIVEexp*10 + Pexp * 9 + Rexp * 2 + Sexp * 4},
  {184229, FIVEexp*7 + Pexp * 6 + Qexp * 2 + Rexp * 2 + Sexp * 4},
  {-278416, FIVEexp*9 + Pexp * 3 + Qexp * 4 + Rexp * 2 + Sexp * 4},
  {-163722, FIVEexp*9 + Qexp * 6 + Rexp * 2 + Sexp * 4},
  {-269418, FIVEexp*7 + Pexp * 7 + Rexp * 3 + Sexp * 4},
  {-123232, FIVEexp*9 + Pexp * 4 + Qexp * 2 + Rexp * 3 + Sexp * 4},
  {-23364, FIVEexp*10 + Pexp + Qexp * 4 + Rexp * 3 + Sexp * 4},
  {11808, FIVEexp*9 + Pexp * 5 + Rexp * 4 + Sexp * 4},
  {30762, FIVEexp*10 + Pexp * 2 + Qexp * 2 + Rexp * 4 + Sexp * 4},
  {9568, FIVEexp*10 + Pexp * 3 + Rexp * 5 + Sexp * 4},
  {-23368, FIVEexp*11 + Qexp * 2 + Rexp * 5 + Sexp * 4},
  {-768, FIVEexp*12 + Pexp + Rexp * 6 + Sexp * 4},
  {4509, FIVEexp*8 + Pexp * 9 + Qexp + Sexp * 5},
  {-1627, FIVEexp*9 + Pexp * 6 + Qexp * 3 + Sexp * 5},
  {6146, FIVEexp*10 + Pexp * 3 + Qexp * 5 + Sexp * 5},
  {444, FIVEexp*12 + Qexp * 7 + Sexp * 5},
  {61713, FIVEexp*8 + Pexp * 7 + Qexp + Rexp + Sexp * 5},
  {8798, FIVEexp*11 + Pexp * 4 + Qexp * 3 + Rexp + Sexp * 5},
  {42022, FIVEexp*10 + Pexp + Qexp * 5 + Rexp + Sexp * 5},
  {-12074, FIVEexp*9 + Pexp * 5 + Qexp + Rexp * 2 + Sexp * 5},
  {4146, FIVEexp*11 + Pexp * 2 + Qexp * 3 + Rexp * 2 + Sexp * 5},
  {-7848, FIVEexp*11 + Pexp * 3 + Qexp + Rexp * 3 + Sexp * 5},
  {9146, FIVEexp*12 + Qexp * 3 + Rexp * 3 + Sexp * 5},
  {6144, FIVEexp*12 + Pexp + Qexp + Rexp * 4 + Sexp * 5},
  {-1404, FIVEexp*10 + Pexp * 8 + Sexp * 6},
  {-4154, FIVEexp*11 + Pexp * 5 + Qexp * 2 + Sexp * 6},
  {-2178, FIVEexp*12 + Pexp * 2 + Qexp * 4 + Sexp * 6},
  {7509, FIVEexp*10 + Pexp * 6 + Rexp + Sexp * 6},
  {-16, FIVEexp*12 + Pexp * 3 + Qexp * 2 + Rexp + Sexp * 6},
  {-5226, FIVEexp*12 + Qexp * 4 + Rexp + Sexp * 6},
  {-496, FIVEexp*12 + Pexp * 4 + Rexp * 2 + Sexp * 6},
  {-13552, FIVEexp*12 + Pexp + Qexp * 2 + Rexp * 2 + Sexp * 6},
  {74, FIVEexp*12 + Pexp * 2 + Rexp * 3 + Sexp * 6},
  {-8, FIVEexp*15 + Rexp * 4 + Sexp * 6},
  {268, FIVEexp*13 + Pexp * 4 + Qexp + Sexp * 7},
  {274, FIVEexp*14 + Pexp + Qexp * 3 + Sexp * 7},
  {366, FIVEexp*13 + Pexp * 2 + Qexp + Rexp + Sexp * 7},
  {202, FIVEexp*14 + Qexp + Rexp * 2 + Sexp * 7},
  {-9, FIVEexp*15 + Pexp * 3 + Sexp * 8},
  {-7, FIVEexp*16 + Qexp * 2 + Sexp * 8},
  {-2, FIVEexp*15 + Pexp + Rexp + Sexp * 8},
  {0, END_COEFF},

  // b[4, 2]
  {8, FIVEexp*2 + Pexp * 9 + Qexp * 8},
  {302, FIVEexp*2 + Pexp * 6 + Qexp * 10},
  {3146, FIVEexp*2 + Pexp * 3 + Qexp * 12},
  {9936, FIVEexp*2 + Qexp * 14},
  {-192, FIVEexp*2 + Pexp * 10 + Qexp * 6 + Rexp},
  {-6572, FIVEexp*2 + Pexp * 7 + Qexp * 8 + Rexp},
  {-68383, FIVEexp*2 + Pexp * 4 + Qexp * 10 + Rexp},
  {-44532, FIVEexp*3 + Pexp + Qexp * 12 + Rexp},
  {1242, FIVEexp*2 + Pexp * 11 + Qexp * 4 + Rexp * 2},
  {44647, FIVEexp*2 + Pexp * 8 + Qexp * 6 + Rexp * 2},
  {506986, FIVEexp*2 + Pexp * 5 + Qexp * 8 + Rexp * 2},
  {1813354, FIVEexp*2 + Pexp * 2 + Qexp * 10 + Rexp * 2},
  {-486, FIVEexp*3 + Pexp * 12 + Qexp * 2 + Rexp * 3},
  {-114909, FIVEexp*2 + Pexp * 9 + Qexp * 4 + Rexp * 3},
  {-1616122, FIVEexp*2 + Pexp * 6 + Qexp * 6 + Rexp * 3},
  {-277703, FIVEexp*4 + Pexp * 3 + Qexp * 8 + Rexp * 3},
  {-89938, FIVEexp*3 + Qexp * 10 + Rexp * 3},
  {86964, FIVEexp*2 + Pexp * 10 + Qexp * 2 + Rexp * 4},
  {86416, FIVEexp*4 + Pexp * 7 + Qexp * 4 + Rexp * 4},
  {2648599, FIVEexp*3 + Pexp * 4 + Qexp * 6 + Rexp * 4},
  {182678, FIVEexp*4 + Pexp + Qexp * 8 + Rexp * 4},
  {-198868, FIVEexp*3 + Pexp * 8 + Qexp * 2 + Rexp * 5},
  {-19256, FIVEexp*6 + Pexp * 5 + Qexp * 4 + Rexp * 5},
  {-511089, FIVEexp*4 + Pexp * 2 + Qexp * 6 + Rexp * 5},
  {111696, FIVEexp*4 + Pexp * 6 + Qexp * 2 + Rexp * 6},
  {-7648, FIVEexp*5 + Pexp * 3 + Qexp * 4 + Rexp * 6},
  {-94292, FIVEexp*5 + Qexp * 6 + Rexp * 6},
  {167744, FIVEexp*5 + Pexp * 4 + Qexp * 2 + Rexp * 7},
  {91648, FIVEexp*6 + Pexp + Qexp * 4 + Rexp * 7},
  {-29952, FIVEexp*7 + Pexp * 2 + Qexp * 2 + Rexp * 8},
  {216, FIVEexp*2 + Pexp * 11 + Qexp * 5 + Sexp},
  {12416, FIVEexp*2 + Pexp * 8 + Qexp * 7 + Sexp},
  {143669, FIVEexp*2 + Pexp * 5 + Qexp * 9 + Sexp},
  {92454, FIVEexp*3 + Pexp * 2 + Qexp * 11 + Sexp},
  {-4212, FIVEexp*2 + Pexp * 12 + Qexp * 3 + Rexp + Sexp},
  {-169386, FIVEexp*2 + Pexp * 9 + Qexp * 5 + Rexp + Sexp},
  {-399431, FIVEexp*3 + Pexp * 6 + Qexp * 7 + Rexp + Sexp},
  {-11141, FIVEexp*6 + Pexp * 3 + Qexp * 9 + Rexp + Sexp},
  {1152, FIVEexp*6 + Qexp * 11 + Rexp + Sexp},
  {2916, FIVEexp*3 + Pexp * 13 + Qexp + Rexp * 2 + Sexp},
  {630522, FIVEexp*2 + Pexp * 10 + Qexp * 3 + Rexp * 2 + Sexp},
  {8807496, FIVEexp*2 + Pexp * 7 + Qexp * 5 + Rexp * 2 + Sexp},
  {11899, FIVEexp*7 + Pexp * 4 + Qexp * 7 + Rexp * 2 + Sexp},
  {-13969, FIVEexp*5 + Pexp + Qexp * 9 + Rexp * 2 + Sexp},
  {-537084, FIVEexp*2 + Pexp * 11 + Qexp + Rexp * 3 + Sexp},
  {-553706, FIVEexp*4 + Pexp * 8 + Qexp * 3 + Rexp * 3 + Sexp},
  {-18301387, FIVEexp*3 + Pexp * 5 + Qexp * 5 + Rexp * 3 + Sexp},
  {-449249, FIVEexp*5 + Pexp * 2 + Qexp * 7 + Rexp * 3 + Sexp},
  {1476688, FIVEexp*3 + Pexp * 9 + Qexp + Rexp * 4 + Sexp},
  {4773536, FIVEexp*4 + Pexp * 6 + Qexp * 3 + Rexp * 4 + Sexp},
  {2792262, FIVEexp*5 + Pexp * 3 + Qexp * 5 + Rexp * 4 + Sexp},
  {6471, FIVEexp*8 + Qexp * 7 + Rexp * 4 + Sexp},
  {-2055168, FIVEexp*4 + Pexp * 7 + Qexp + Rexp * 5 + Sexp},
  {-840848, FIVEexp*6 + Pexp * 4 + Qexp * 3 + Rexp * 5 + Sexp},
  {-896104, FIVEexp*6 + Pexp + Qexp * 5 + Rexp * 5 + Sexp},
  {1351936, FIVEexp*5 + Pexp * 5 + Qexp + Rexp * 6 + Sexp},
  {351488, FIVEexp*7 + Pexp * 2 + Qexp * 3 + Rexp * 6 + Sexp},
  {-48128, FIVEexp*7 + Pexp * 3 + Qexp + Rexp * 7 + Sexp},
  {9984, FIVEexp*8 + Qexp * 3 + Rexp * 7 + Sexp},
  {1458, FIVEexp*2 + Pexp * 13 + Qexp * 2 + Sexp * 2},
  {110619, FIVEexp*2 + Pexp * 10 + Qexp * 4 + Sexp * 2},
  {272221, FIVEexp*3 + Pexp * 7 + Qexp * 6 + Sexp * 2},
  {6232, FIVEexp*6 + Pexp * 4 + Qexp * 8 + Sexp * 2},
  {-28248, FIVEexp*5 + Pexp + Qexp * 10 + Sexp * 2},
  {-4374, FIVEexp*3 + Pexp * 14 + Rexp + Sexp * 2},
  {-175689, FIVEexp*3 + Pexp * 11 + Qexp * 2 + Rexp + Sexp * 2},
  {-87379, FIVEexp*5 + Pexp * 8 + Qexp * 4 + Rexp + Sexp * 2},
  {-9748, FIVEexp*7 + Pexp * 5 + Qexp * 6 + Rexp + Sexp * 2},
  {119658, FIVEexp*6 + Pexp * 2 + Qexp * 8 + Rexp + Sexp * 2},
  {821826, FIVEexp*2 + Pexp * 12 + Rexp * 2 + Sexp * 2},
  {151659, FIVEexp*5 + Pexp * 9 + Qexp * 2 + Rexp * 2 + Sexp * 2},
  {562577, FIVEexp*5 + Pexp * 6 + Qexp * 4 + Rexp * 2 + Sexp * 2},
  {-22383, FIVEexp*8 + Pexp * 3 + Qexp * 6 + Rexp * 2 + Sexp * 2},
  {-53176, FIVEexp*7 + Qexp * 8 + Rexp * 2 + Sexp * 2},
  {-2372472, FIVEexp*3 + Pexp * 10 + Rexp * 3 + Sexp * 2},
  {-1300818, FIVEexp*5 + Pexp * 7 + Qexp * 2 + Rexp * 3 + Sexp * 2},
  {-477, FIVEexp*8 + Pexp * 4 + Qexp * 4 + Rexp * 3 + Sexp * 2},
  {248569, FIVEexp*7 + Pexp + Qexp * 6 + Rexp * 3 + Sexp * 2},
  {744404, FIVEexp*5 + Pexp * 8 + Rexp * 4 + Sexp * 2},
  {271224, FIVEexp*7 + Pexp * 5 + Qexp * 2 + Rexp * 4 + Sexp * 2},
  {-67333, FIVEexp*8 + Pexp * 2 + Qexp * 4 + Rexp * 4 + Sexp * 2},
  {-672848, FIVEexp*6 + Pexp * 6 + Rexp * 5 + Sexp * 2},
  {-35808, FIVEexp*9 + Pexp * 3 + Qexp * 2 + Rexp * 5 + Sexp * 2},
  {-21636, FIVEexp*9 + Qexp * 4 + Rexp * 5 + Sexp * 2},
  {11968, FIVEexp*9 + Pexp * 4 + Rexp * 6 + Sexp * 2},
  {20864, FIVEexp*9 + Pexp + Qexp * 2 + Rexp * 6 + Sexp * 2},
  {-9984, FIVEexp*9 + Pexp * 2 + Rexp * 7 + Sexp * 2},
  {32076, FIVEexp*3 + Pexp * 12 + Qexp + Sexp * 3},
  {2313, FIVEexp*6 + Pexp * 9 + Qexp * 3 + Sexp * 3},
  {-107347, FIVEexp*5 + Pexp * 6 + Qexp * 5 + Sexp * 3},
  {-276, FIVEexp*10 + Pexp * 3 + Qexp * 7 + Sexp * 3},
  {-19024, FIVEexp*7 + Qexp * 9 + Sexp * 3},
  {6561, FIVEexp*6 + Pexp * 10 + Qexp + Rexp + Sexp * 3},
  {256436, FIVEexp*6 + Pexp * 7 + Qexp * 3 + Rexp + Sexp * 3},
  {353147, FIVEexp*7 + Pexp * 4 + Qexp * 5 + Rexp + Sexp * 3},
  {2068, FIVEexp*10 + Pexp + Qexp * 7 + Rexp + Sexp * 3},
  {-35748, FIVEexp*7 + Pexp * 8 + Qexp + Rexp * 2 + Sexp * 3},
  {-564674, FIVEexp*7 + Pexp * 5 + Qexp * 3 + Rexp * 2 + Sexp * 3},
  {-929401, FIVEexp*7 + Pexp * 2 + Qexp * 5 + Rexp * 2 + Sexp * 3},
  {240032, FIVEexp*7 + Pexp * 6 + Qexp + Rexp * 3 + Sexp * 3},
  {111688, FIVEexp*9 + Pexp * 3 + Qexp * 3 + Rexp * 3 + Sexp * 3},
  {281447, FIVEexp*8 + Qexp * 5 + Rexp * 3 + Sexp * 3},
  {-56048, FIVEexp*8 + Pexp * 4 + Qexp + Rexp * 4 + Sexp * 3},
  {-33376, FIVEexp*9 + Pexp + Qexp * 3 + Rexp * 4 + Sexp * 3},
  {-15872, FIVEexp*9 + Pexp * 2 + Qexp + Rexp * 5 + Sexp * 3},
  {9984, FIVEexp*10 + Qexp + Rexp * 6 + Sexp * 3},
  {-27702, FIVEexp*5 + Pexp * 11 + Sexp * 4},
  {-25029, FIVEexp*7 + Pexp * 8 + Qexp * 2 + Sexp * 4},
  {-114698, FIVEexp*7 + Pexp * 5 + Qexp * 4 + Sexp * 4},
  {-3476, FIVEexp*8 + Pexp * 2 + Qexp * 6 + Sexp * 4},
  {21222, FIVEexp*7 + Pexp * 9 + Rexp + Sexp * 4},
  {26796, FIVEexp*8 + Pexp * 6 + Qexp * 2 + Rexp + Sexp * 4},
  {-5662, FIVEexp*10 + Pexp * 3 + Qexp * 4 + Rexp + Sexp * 4},
  {-6214, FIVEexp*10 + Qexp * 6 + Rexp + Sexp * 4},
  {-146862, FIVEexp*7 + Pexp * 7 + Rexp * 2 + Sexp * 4},
  {-12624, FIVEexp*10 + Pexp * 4 + Qexp * 2 + Rexp * 2 + Sexp * 4},
  {-16887, FIVEexp*10 + Pexp + Qexp * 4 + Rexp * 2 + Sexp * 4},
  {33664, FIVEexp*8 + Pexp * 5 + Rexp * 3 + Sexp * 4},
  {19533, FIVEexp*10 + Pexp * 2 + Qexp * 2 + Rexp * 3 + Sexp * 4},
  {6336, FIVEexp*10 + Pexp * 3 + Rexp * 4 + Sexp * 4},
  {-1916, FIVEexp*12 + Qexp * 2 + Rexp * 4 + Sexp * 4},
  {-2432, FIVEexp*11 + Pexp + Rexp * 5 + Sexp * 4},
  {19413, FIVEexp*8 + Pexp * 7 + Qexp + Sexp * 5},
  {1119, FIVEexp*11 + Pexp * 4 + Qexp * 3 + Sexp * 5},
  {4037, FIVEexp*10 + Pexp + Qexp * 5 + Sexp * 5},
  {3317, FIVEexp*10 + Pexp * 5 + Qexp + Rexp + Sexp * 5},
  {5704, FIVEexp*11 + Pexp * 2 + Qexp * 3 + Rexp + Sexp * 5},
  {-1222, FIVEexp*12 + Pexp * 3 + Qexp + Rexp * 2 + Sexp * 5},
  {2297, FIVEexp*12 + Qexp * 3 + Rexp * 2 + Sexp * 5},
  {2952, FIVEexp*12 + Pexp + Qexp + Rexp * 3 + Sexp * 5},
  {-2016, FIVEexp*10 + Pexp * 6 + Sexp * 6},
  {-654, FIVEexp*12 + Pexp * 3 + Qexp * 2 + Sexp * 6},
  {-296, FIVEexp*12 + Qexp * 4 + Sexp * 6},
  {478, FIVEexp*12 + Pexp * 4 + Rexp + Sexp * 6},
  {-757, FIVEexp*13 + Pexp + Qexp * 2 + Rexp + Sexp * 6},
  {-639, FIVEexp*12 + Pexp * 2 + Rexp * 2 + Sexp * 6},
  {-92, FIVEexp*13 + Rexp * 3 + Sexp * 6},
  {286, FIVEexp*13 + Pexp * 2 + Qexp + Sexp * 7},
  {13, FIVEexp*15 + Qexp + Rexp + Sexp * 7},
  {-7, FIVEexp*15 + Pexp + Sexp * 8},
  {0, END_COEFF},

  // b[4, 3]
  {-24, FIVEexp*2 + Pexp * 10 + Qexp * 6},
  {-738, FIVEexp*2 + Pexp * 7 + Qexp * 8},
  {-1392, FIVEexp*3 + Pexp * 4 + Qexp * 10},
  {-20736, FIVEexp*2 + Pexp + Qexp * 12},
  {216, FIVEexp*2 + Pexp * 11 + Qexp * 4 + Rexp},
  {7902, FIVEexp*2 + Pexp * 8 + Qexp * 6 + Rexp},
  {85911, FIVEexp*2 + Pexp * 5 + Qexp * 8 + Rexp},
  {288792, FIVEexp*2 + Pexp * 2 + Qexp * 10 + Rexp},
  {-486, FIVEexp*2 + Pexp * 12 + Qexp * 2 + Rexp * 2},
  {-26488, FIVEexp*2 + Pexp * 9 + Qexp * 4 + Rexp * 2},
  {-370991, FIVEexp*2 + Pexp * 6 + Qexp * 6 + Rexp * 2},
  {-61329, FIVEexp*4 + Pexp * 3 + Qexp * 8 + Rexp * 2},
  {-44064, FIVEexp*3 + Qexp * 10 + Rexp * 2},
  {26262, FIVEexp*2 + Pexp * 10 + Qexp * 2 + Rexp * 3},
  {25974, FIVEexp*4 + Pexp * 7 + Qexp * 4 + Rexp * 3},
  {778687, FIVEexp*3 + Pexp * 4 + Qexp * 6 + Rexp * 3},
  {93234, FIVEexp*4 + Pexp + Qexp * 8 + Rexp * 3},
  {-78764, FIVEexp*3 + Pexp * 8 + Qexp * 2 + Rexp * 4},
  {-191143, FIVEexp*4 + Pexp * 5 + Qexp * 4 + Rexp * 4},
  {-311091, FIVEexp*4 + Pexp * 2 + Qexp * 6 + Rexp * 4},
  {79144, FIVEexp*4 + Pexp * 6 + Qexp * 2 + Rexp * 5},
  {10624, FIVEexp*6 + Pexp * 3 + Qexp * 4 + Rexp * 5},
  {-25854, FIVEexp*5 + Qexp * 6 + Rexp * 5},
  {17408, FIVEexp*5 + Pexp * 4 + Qexp * 2 + Rexp * 6},
  {24176, FIVEexp*6 + Pexp + Qexp * 4 + Rexp * 6},
  {-8064, FIVEexp*7 + Pexp * 2 + Qexp * 2 + Rexp * 7},
  {-648, FIVEexp*2 + Pexp * 12 + Qexp * 3 + Sexp},
  {-18372, FIVEexp*2 + Pexp * 9 + Qexp * 5 + Sexp},
  {-168289, FIVEexp*2 + Pexp * 6 + Qexp * 7 + Sexp},
  {-17324, FIVEexp*4 + Pexp * 3 + Qexp * 9 + Sexp},
  {21816, FIVEexp*4 + Qexp * 11 + Sexp},
  {2916, FIVEexp*2 + Pexp * 13 + Qexp + Rexp + Sexp},
  {115092, FIVEexp*2 + Pexp * 10 + Qexp * 3 + Rexp + Sexp},
  {1329588, FIVEexp*2 + Pexp * 7 + Qexp * 5 + Rexp + Sexp},
  {171329, FIVEexp*4 + Pexp * 4 + Qexp * 7 + Rexp + Sexp},
  {-36792, FIVEexp*5 + Pexp + Qexp * 9 + Rexp + Sexp},
  {-144072, FIVEexp*2 + Pexp * 11 + Qexp + Rexp * 2 + Sexp},
  {-120343, FIVEexp*4 + Pexp * 8 + Qexp * 3 + Rexp * 2 + Sexp},
  {-3096586, FIVEexp*3 + Pexp * 5 + Qexp * 5 + Rexp * 2 + Sexp},
  {17773, FIVEexp*5 + Pexp * 2 + Qexp * 7 + Rexp * 2 + Sexp},
  {430344, FIVEexp*3 + Pexp * 9 + Qexp + Rexp * 3 + Sexp},
  {1099483, FIVEexp*4 + Pexp * 6 + Qexp * 3 + Rexp * 3 + Sexp},
  {534406, FIVEexp*5 + Pexp * 3 + Qexp * 5 + Rexp * 3 + Sexp},
  {1703, FIVEexp*8 + Qexp * 7 + Rexp * 3 + Sexp},
  {-626512, FIVEexp*4 + Pexp * 7 + Qexp + Rexp * 4 + Sexp},
  {-218888, FIVEexp*6 + Pexp * 4 + Qexp * 3 + Rexp * 4 + Sexp},
  {-231016, FIVEexp*6 + Pexp + Qexp * 5 + Rexp * 4 + Sexp},
  {436992, FIVEexp*5 + Pexp * 5 + Qexp + Rexp * 5 + Sexp},
  {92624, FIVEexp*7 + Pexp * 2 + Qexp * 3 + Rexp * 5 + Sexp},
  {-15616, FIVEexp*7 + Pexp * 3 + Qexp + Rexp * 6 + Sexp},
  {2688, FIVEexp*8 + Qexp * 3 + Rexp * 6 + Sexp},
  {-4374, FIVEexp*2 + Pexp * 14 + Sexp * 2},
  {-122634, FIVEexp*2 + Pexp * 11 + Qexp * 2 + Sexp * 2},
  {-43053, FIVEexp*4 + Pexp * 8 + Qexp * 4 + Sexp * 2},
  {-71371, FIVEexp*4 + Pexp * 5 + Qexp * 6 + Sexp * 2},
  {17268, FIVEexp*6 + Pexp * 2 + Qexp * 8 + Sexp * 2},
  {208008, FIVEexp*2 + Pexp * 12 + Rexp + Sexp * 2},
  {130923, FIVEexp*4 + Pexp * 9 + Qexp * 2 + Rexp + Sexp * 2},
  {49721, FIVEexp*5 + Pexp * 6 + Qexp * 4 + Rexp + Sexp * 2},
  {-123917, FIVEexp*6 + Pexp * 3 + Qexp * 6 + Rexp + Sexp * 2},
  {-2556, FIVEexp*8 + Qexp * 8 + Rexp + Sexp * 2},
  {-617166, FIVEexp*3 + Pexp * 10 + Rexp * 2 + Sexp * 2},
  {-238569, FIVEexp*5 + Pexp * 7 + Qexp * 2 + Rexp * 2 + Sexp * 2},
  {8748, FIVEexp*7 + Pexp * 4 + Qexp * 4 + Rexp * 2 + Sexp * 2},
  {52267, FIVEexp*7 + Pexp + Qexp * 6 + Rexp * 2 + Sexp * 2},
  {189852, FIVEexp*5 + Pexp * 8 + Rexp * 3 + Sexp * 2},
  {61432, FIVEexp*7 + Pexp * 5 + Qexp * 2 + Rexp * 3 + Sexp * 2},
  {-10669, FIVEexp*8 + Pexp * 2 + Qexp * 4 + Rexp * 3 + Sexp * 2},
  {-174792, FIVEexp*6 + Pexp * 6 + Rexp * 4 + Sexp * 2},
  {-47792, FIVEexp*8 + Pexp * 3 + Qexp * 2 + Rexp * 4 + Sexp * 2},
  {-5366, FIVEexp*9 + Qexp * 4 + Rexp * 4 + Sexp * 2},
  {15872, FIVEexp*8 + Pexp * 4 + Rexp * 5 + Sexp * 2},
  {6016, FIVEexp*9 + Pexp + Qexp * 2 + Rexp * 5 + Sexp * 2},
  {-2688, FIVEexp*9 + Pexp * 2 + Rexp * 6 + Sexp * 2},
  {42444, FIVEexp*4 + Pexp * 10 + Qexp + Sexp * 3},
  {33666, FIVEexp*6 + Pexp * 7 + Qexp * 3 + Sexp * 3},
  {40457, FIVEexp*7 + Pexp * 4 + Qexp * 5 + Sexp * 3},
  {33924, FIVEexp*7 + Pexp + Qexp * 7 + Sexp * 3},
  {-28674, FIVEexp*6 + Pexp * 8 + Qexp + Rexp + Sexp * 3},
  {-85542, FIVEexp*7 + Pexp * 5 + Qexp * 3 + Rexp + Sexp * 3},
  {-149028, FIVEexp*7 + Pexp * 2 + Qexp * 5 + Rexp + Sexp * 3},
  {32681, FIVEexp*7 + Pexp * 6 + Qexp + Rexp * 2 + Sexp * 3},
  {19064, FIVEexp*9 + Pexp * 3 + Qexp * 3 + Rexp * 2 + Sexp * 3},
  {55991, FIVEexp*8 + Qexp * 5 + Rexp * 2 + Sexp * 3},
  {7176, FIVEexp*8 + Pexp * 4 + Qexp + Rexp * 3 + Sexp * 3},
  {-5528, FIVEexp*9 + Pexp + Qexp * 3 + Rexp * 3 + Sexp * 3},
  {-8528, FIVEexp*9 + Pexp * 2 + Qexp + Rexp * 4 + Sexp * 3},
  {2688, FIVEexp*10 + Qexp + Rexp * 5 + Sexp * 3},
  {-972, FIVEexp*7 + Pexp * 9 + Sexp * 4},
  {-9012, FIVEexp*7 + Pexp * 6 + Qexp * 2 + Sexp * 4},
  {-4282, FIVEexp*9 + Pexp * 3 + Qexp * 4 + Sexp * 4},
  {-5324, FIVEexp*9 + Qexp * 6 + Sexp * 4},
  {2529, FIVEexp*7 + Pexp * 7 + Rexp + Sexp * 4},
  {-8418, FIVEexp*9 + Pexp * 4 + Qexp * 2 + Rexp + Sexp * 4},
  {-2486, FIVEexp*10 + Pexp + Qexp * 4 + Rexp + Sexp * 4},
  {-14633, FIVEexp*8 + Pexp * 5 + Rexp * 2 + Sexp * 4},
  {3249, FIVEexp*10 + Pexp * 2 + Qexp * 2 + Rexp * 2 + Sexp * 4},
  {2768, FIVEexp*10 + Pexp * 3 + Rexp * 3 + Sexp * 4},
  {-378, FIVEexp*12 + Qexp * 2 + Rexp * 3 + Sexp * 4},
  {-688, FIVEexp*11 + Pexp + Rexp * 4 + Sexp * 4},
  {5322, FIVEexp*9 + Pexp * 5 + Qexp + Sexp * 5},
  {779, FIVEexp*11 + Pexp * 2 + Qexp * 3 + Sexp * 5},
  {-986, FIVEexp*11 + Pexp * 3 + Qexp + Rexp + Sexp * 5},
  {61, FIVEexp*13 + Qexp * 3 + Rexp + Sexp * 5},
  {496, FIVEexp*12 + Pexp + Qexp + Rexp * 2 + Sexp * 5},
  {-9, FIVEexp*12 + Pexp * 4 + Sexp * 6},
  {-379, FIVEexp*12 + Pexp + Qexp * 2 + Sexp * 6},
  {63, FIVEexp*12 + Pexp * 2 + Rexp + Sexp * 6},
  {-26, FIVEexp*13 + Rexp * 2 + Sexp * 6},
  {9, FIVEexp*14 + Qexp + Sexp * 7},
  {0, END_COEFF},

  // b[4, 4]
  {8, FIVEexp*2 + Pexp * 8 + Qexp * 6},
  {118, FIVEexp*2 + Pexp * 5 + Qexp * 8},
  {432, FIVEexp*2 + Pexp * 2 + Qexp * 10},
  {-72, FIVEexp*2 + Pexp * 9 + Qexp * 4 + Rexp},
  {-1986, FIVEexp*2 + Pexp * 6 + Qexp * 6 + Rexp},
  {-3227, FIVEexp*3 + Pexp * 3 + Qexp * 8 + Rexp},
  {-7992, FIVEexp*3 + Qexp * 10 + Rexp},
  {162, FIVEexp*2 + Pexp * 10 + Qexp * 2 + Rexp * 2},
  {1893, FIVEexp*3 + Pexp * 7 + Qexp * 4 + Rexp * 2},
  {24876, FIVEexp*3 + Pexp * 4 + Qexp * 6 + Rexp * 2},
  {18342, FIVEexp*4 + Pexp + Qexp * 8 + Rexp * 2},
  {-2652, FIVEexp*3 + Pexp * 8 + Qexp * 2 + Rexp * 3},
  {-12509, FIVEexp*4 + Pexp * 5 + Qexp * 4 + Rexp * 3},
  {-66258, FIVEexp*4 + Pexp * 2 + Qexp * 6 + Rexp * 3},
  {7652, FIVEexp*4 + Pexp * 6 + Qexp * 2 + Rexp * 4},
  {15192, FIVEexp*5 + Pexp * 3 + Qexp * 4 + Rexp * 4},
  {-5353, FIVEexp*5 + Qexp * 6 + Rexp * 4},
  {-2784, FIVEexp*5 + Pexp * 4 + Qexp * 2 + Rexp * 5},
  {5232, FIVEexp*6 + Pexp + Qexp * 4 + Rexp * 5},
  {-1728, FIVEexp*7 + Pexp * 2 + Qexp * 2 + Rexp * 6},
  {216, FIVEexp*2 + Pexp * 10 + Qexp * 3 + Sexp},
  {5768, FIVEexp*2 + Pexp * 7 + Qexp * 5 + Sexp},
  {1503, FIVEexp*4 + Pexp * 4 + Qexp * 7 + Sexp},
  {324, FIVEexp*5 + Pexp + Qexp * 9 + Sexp},
  {-972, FIVEexp*2 + Pexp * 11 + Qexp + Rexp + Sexp},
  {-9354, FIVEexp*3 + Pexp * 8 + Qexp * 3 + Rexp + Sexp},
  {-112218, FIVEexp*3 + Pexp * 5 + Qexp * 5 + Rexp + Sexp},
  {-14223, FIVEexp*5 + Pexp * 2 + Qexp * 7 + Rexp + Sexp},
  {16092, FIVEexp*3 + Pexp * 9 + Qexp + Rexp * 2 + Sexp},
  {78929, FIVEexp*4 + Pexp * 6 + Qexp * 3 + Rexp * 2 + Sexp},
  {87043, FIVEexp*5 + Pexp * 3 + Qexp * 5 + Rexp * 2 + Sexp},
  {1332, FIVEexp*7 + Qexp * 7 + Rexp * 2 + Sexp},
  {-55456, FIVEexp*4 + Pexp * 7 + Qexp + Rexp * 3 + Sexp},
  {-29124, FIVEexp*6 + Pexp * 4 + Qexp * 3 + Rexp * 3 + Sexp},
  {-40058, FIVEexp*6 + Pexp + Qexp * 5 + Rexp * 3 + Sexp},
  {67264, FIVEexp*5 + Pexp * 5 + Qexp + Rexp * 4 + Sexp},
  {16624, FIVEexp*7 + Pexp * 2 + Qexp * 3 + Rexp * 4 + Sexp},
  {-3072, FIVEexp*7 + Pexp * 3 + Qexp + Rexp * 5 + Sexp},
  {576, FIVEexp*8 + Qexp * 3 + Rexp * 5 + Sexp},
  {1458, FIVEexp*2 + Pexp * 12 + Sexp * 2},
  {1971, FIVEexp*4 + Pexp * 9 + Qexp * 2 + Sexp * 2},
  {3428, FIVEexp*5 + Pexp * 6 + Qexp * 4 + Sexp * 2},
  {278, FIVEexp*7 + Pexp * 3 + Qexp * 6 + Sexp * 2},
  {216, FIVEexp*7 + Qexp * 8 + Sexp * 2},
  {-22518, FIVEexp*3 + Pexp * 10 + Rexp + Sexp * 2},
  {-21636, FIVEexp*5 + Pexp * 7 + Qexp * 2 + Rexp + Sexp * 2},
  {-4418, FIVEexp*7 + Pexp * 4 + Qexp * 4 + Rexp + Sexp * 2},
  {-3624, FIVEexp*7 + Pexp + Qexp * 6 + Rexp + Sexp * 2},
  {16416, FIVEexp*5 + Pexp * 8 + Rexp * 2 + Sexp * 2},
  {9402, FIVEexp*7 + Pexp * 5 + Qexp * 2 + Rexp * 2 + Sexp * 2},
  {3244, FIVEexp*8 + Pexp * 2 + Qexp * 4 + Rexp * 2 + Sexp * 2},
  {-24596, FIVEexp*6 + Pexp * 6 + Rexp * 3 + Sexp * 2},
  {-10016, FIVEexp*8 + Pexp * 3 + Qexp * 2 + Rexp * 3 + Sexp * 2},
  {-933, FIVEexp*9 + Qexp * 4 + Rexp * 3 + Sexp * 2},
  {2848, FIVEexp*8 + Pexp * 4 + Rexp * 4 + Sexp * 2},
  {1248, FIVEexp*9 + Pexp + Qexp * 2 + Rexp * 4 + Sexp * 2},
  {-576, FIVEexp*9 + Pexp * 2 + Rexp * 5 + Sexp * 2},
  {-929, FIVEexp*7 + Pexp * 5 + Qexp * 3 + Sexp * 3},
  {-2423, FIVEexp*7 + Pexp * 2 + Qexp * 5 + Sexp * 3},
  {1629, FIVEexp*7 + Pexp * 6 + Qexp + Rexp + Sexp * 3},
  {3624, FIVEexp*8 + Pexp * 3 + Qexp * 3 + Rexp + Sexp * 3},
  {3148, FIVEexp*8 + Qexp * 5 + Rexp + Sexp * 3},
  {3708, FIVEexp*8 + Pexp * 4 + Qexp + Rexp * 2 + Sexp * 3},
  {1136, FIVEexp*9 + Pexp + Qexp * 3 + Rexp * 2 + Sexp * 3},
  {-2064, FIVEexp*9 + Pexp * 2 + Qexp + Rexp * 3 + Sexp * 3},
  {576, FIVEexp*10 + Qexp + Rexp * 4 + Sexp * 3},
  {-1701, FIVEexp*7 + Pexp * 7 + Sexp * 4},
  {-271, FIVEexp*9 + Pexp * 4 + Qexp * 2 + Sexp * 4},
  {-18, FIVEexp*10 + Pexp + Qexp * 4 + Sexp * 4},
  {-1029, FIVEexp*8 + Pexp * 5 + Rexp + Sexp * 4},
  {-459, FIVEexp*10 + Pexp * 2 + Qexp * 2 + Rexp + Sexp * 4},
  {424, FIVEexp*10 + Pexp * 3 + Rexp * 2 + Sexp * 4},
  {-43, FIVEexp*12 + Qexp * 2 + Rexp * 2 + Sexp * 4},
  {-144, FIVEexp*11 + Pexp + Rexp * 3 + Sexp * 4},
  {1, FIVEexp*13 + Pexp * 3 + Qexp + Sexp * 5},
  {8, FIVEexp*12 + Qexp * 3 + Sexp * 5},
  {58, FIVEexp*12 + Pexp + Qexp + Rexp + Sexp * 5},
  {-7, FIVEexp*12 + Pexp * 2 + Sexp * 6},
  {-3, FIVEexp*13 + Rexp + Sexp * 6},
  {0, END_COEFF},

  // b[4, 5]
  {-184, FIVEexp*2 + Pexp * 6 + Qexp * 6},
  {-2714, FIVEexp*2 + Pexp * 3 + Qexp * 8},
  {-9936, FIVEexp*2 + Qexp * 10},
  {1556, FIVEexp*2 + Pexp * 7 + Qexp * 4 + Rexp},
  {27183, FIVEexp*2 + Pexp * 4 + Qexp * 6 + Rexp},
  {22932, FIVEexp*3 + Pexp + Qexp * 8 + Rexp},
  {-3276, FIVEexp*2 + Pexp * 8 + Qexp * 2 + Rexp * 2},
  {-16078, FIVEexp*3 + Pexp * 5 + Qexp * 4 + Rexp * 2},
  {-17254, FIVEexp*4 + Pexp * 2 + Qexp * 6 + Rexp * 2},
  {2366, FIVEexp*4 + Pexp * 6 + Qexp * 2 + Rexp * 3},
  {4533, FIVEexp*5 + Pexp * 3 + Qexp * 4 + Rexp * 3},
  {-878, FIVEexp*5 + Qexp * 6 + Rexp * 3},
  {-1744, FIVEexp*5 + Pexp * 4 + Qexp * 2 + Rexp * 4},
  {812, FIVEexp*6 + Pexp + Qexp * 4 + Rexp * 4},
  {-288, FIVEexp*7 + Pexp * 2 + Qexp * 2 + Rexp * 5},
  {-4068, FIVEexp*2 + Pexp * 8 + Qexp * 3 + Sexp},
  {-68039, FIVEexp*2 + Pexp * 5 + Qexp * 5 + Sexp},
  {-11298, FIVEexp*4 + Pexp * 2 + Qexp * 7 + Sexp},
  {16956, FIVEexp*2 + Pexp * 9 + Qexp + Rexp + Sexp},
  {74339, FIVEexp*3 + Pexp * 6 + Qexp * 3 + Rexp + Sexp},
  {80702, FIVEexp*4 + Pexp * 3 + Qexp * 5 + Rexp + Sexp},
  {6552, FIVEexp*5 + Qexp * 7 + Rexp + Sexp},
  {-12564, FIVEexp*4 + Pexp * 7 + Qexp + Rexp * 2 + Sexp},
  {-5617, FIVEexp*6 + Pexp * 4 + Qexp * 3 + Rexp * 2 + Sexp},
  {-7783, FIVEexp*6 + Pexp + Qexp * 5 + Rexp * 2 + Sexp},
  {15264, FIVEexp*5 + Pexp * 5 + Qexp + Rexp * 3 + Sexp},
  {3388, FIVEexp*7 + Pexp * 2 + Qexp * 3 + Rexp * 3 + Sexp},
  {-832, FIVEexp*7 + Pexp * 3 + Qexp + Rexp * 4 + Sexp},
  {96, FIVEexp*8 + Qexp * 3 + Rexp * 4 + Sexp},
  {-21384, FIVEexp*2 + Pexp * 10 + Sexp * 2},
  {-16551, FIVEexp*4 + Pexp * 7 + Qexp * 2 + Sexp * 2},
  {-3511, FIVEexp*6 + Pexp * 4 + Qexp * 4 + Sexp * 2},
  {-516, FIVEexp*7 + Pexp + Qexp * 6 + Sexp * 2},
  {16254, FIVEexp*4 + Pexp * 8 + Rexp + Sexp * 2},
  {37689, FIVEexp*5 + Pexp * 5 + Qexp * 2 + Rexp + Sexp * 2},
  {2463, FIVEexp*7 + Pexp * 2 + Qexp * 4 + Rexp + Sexp * 2},
  {-4518, FIVEexp*6 + Pexp * 6 + Rexp * 2 + Sexp * 2},
  {-8772, FIVEexp*7 + Pexp * 3 + Qexp * 2 + Rexp * 2 + Sexp * 2},
  {-856, FIVEexp*8 + Qexp * 4 + Rexp * 2 + Sexp * 2},
  {496, FIVEexp*8 + Pexp * 4 + Rexp * 3 + Sexp * 2},
  {256, FIVEexp*9 + Pexp + Qexp * 2 + Rexp * 3 + Sexp * 2},
  {-96, FIVEexp*9 + Pexp * 2 + Rexp * 4 + Sexp * 2},
  {108, FIVEexp*7 + Pexp * 6 + Qexp + Sexp * 3},
  {2038, FIVEexp*7 + Pexp * 3 + Qexp * 3 + Sexp * 3},
  {2824, FIVEexp*7 + Qexp * 5 + Sexp * 3},
  {4529, FIVEexp*7 + Pexp * 4 + Qexp + Rexp + Sexp * 3},
  {1056, FIVEexp*8 + Pexp + Qexp * 3 + Rexp + Sexp * 3},
  {-524, FIVEexp*9 + Pexp * 2 + Qexp + Rexp * 2 + Sexp * 3},
  {96, FIVEexp*10 + Qexp + Rexp * 3 + Sexp * 3},
  {-2637, FIVEexp*7 + Pexp * 5 + Sexp * 4},
  {-359, FIVEexp*9 + Pexp * 2 + Qexp * 2 + Sexp * 4},
  {511, FIVEexp*9 + Pexp * 3 + Rexp + Sexp * 4},
  {-134, FIVEexp*10 + Qexp * 2 + Rexp + Sexp * 4},
  {-28, FIVEexp*11 + Pexp + Rexp * 2 + Sexp * 4},
  {7, FIVEexp*12 + Pexp + Qexp + Sexp * 5},
  {-4, FIVEexp*12 + Sexp * 6},
  {0, END_ARRAY},

  // a[0]
  {-4, FIVEexp*2 + Pexp * 7 + Qexp * 7},
  {-87, FIVEexp*2 + Pexp * 4 + Qexp * 9},
  {-84, FIVEexp*3 + Pexp + Qexp * 11},
  {44, FIVEexp*2 + Pexp * 8 + Qexp * 5 + Rexp},
  {1119, FIVEexp*2 + Pexp * 5 + Qexp * 7 + Rexp},
  {6118, FIVEexp*2 + Pexp * 2 + Qexp * 9 + Rexp},
  {-33, FIVEexp*3 + Pexp * 9 + Qexp * 3 + Rexp * 2},
  {-1031, FIVEexp*3 + Pexp * 6 + Qexp * 5 + Rexp * 2},
  {-33221, FIVEexp*2 + Pexp * 3 + Qexp * 7 + Rexp * 2},
  {2378, FIVEexp*2 + Qexp * 9 + Rexp * 2},
  {216, FIVEexp*2 + Pexp * 10 + Qexp + Rexp * 3},
  {9752, FIVEexp*2 + Pexp * 7 + Qexp * 3 + Rexp * 3},
  {83306, FIVEexp*2 + Pexp * 4 + Qexp * 5 + Rexp * 3},
  {-13357, FIVEexp*2 + Pexp + Qexp * 7 + Rexp * 3},
  {-5568, FIVEexp*2 + Pexp * 8 + Qexp + Rexp * 4},
  {-19248, FIVEexp*3 + Pexp * 5 + Qexp * 3 + Rexp * 4},
  {-4904, FIVEexp*2 + Pexp * 2 + Qexp * 5 + Rexp * 4},
  {50176, FIVEexp*2 + Pexp * 6 + Qexp + Rexp * 5},
  {30208, FIVEexp*3 + Pexp * 3 + Qexp * 3 + Rexp * 5},
  {14656, FIVEexp*3 + Qexp * 5 + Rexp * 5},
  {-37888, FIVEexp*3 + Pexp * 4 + Qexp + Rexp * 6},
  {-10752, FIVEexp*4 + Pexp + Qexp * 3 + Rexp * 6},
  {2048, FIVEexp*5 + Pexp * 2 + Qexp + Rexp * 7},
  {-36, FIVEexp*2 + Pexp * 9 + Qexp * 4 + Sexp},
  {-1496, FIVEexp*2 + Pexp * 6 + Qexp * 6 + Sexp},
  {-2253, FIVEexp*3 + Pexp * 3 + Qexp * 8 + Sexp},
  {-696, FIVEexp*4 + Qexp * 10 + Sexp},
  {54, FIVEexp*3 + Pexp * 10 + Qexp * 2 + Rexp + Sexp},
  {12892, FIVEexp*2 + Pexp * 7 + Qexp * 4 + Rexp + Sexp},
  {108743, FIVEexp*2 + Pexp * 4 + Qexp * 6 + Rexp + Sexp},
  {33714, FIVEexp*3 + Pexp + Qexp * 8 + Rexp + Sexp},
  {-648, FIVEexp*2 + Pexp * 11 + Rexp * 2 + Sexp},
  {-34371, FIVEexp*2 + Pexp * 8 + Qexp * 2 + Rexp * 2 + Sexp},
  {-357019, FIVEexp*2 + Pexp * 5 + Qexp * 4 + Rexp * 2 + Sexp},
  {-115423, FIVEexp*3 + Pexp * 2 + Qexp * 6 + Rexp * 2 + Sexp},
  {18144, FIVEexp*2 + Pexp * 9 + Rexp * 3 + Sexp},
  {401536, FIVEexp*2 + Pexp * 6 + Qexp * 2 + Rexp * 3 + Sexp},
  {27836, FIVEexp*4 + Pexp * 3 + Qexp * 4 + Rexp * 3 + Sexp},
  {-18133, FIVEexp*4 + Qexp * 6 + Rexp * 3 + Sexp},
  {-178048, FIVEexp*2 + Pexp * 7 + Rexp * 4 + Sexp},
  {-5072, FIVEexp*5 + Pexp * 4 + Qexp * 2 + Rexp * 4 + Sexp},
  {2176, FIVEexp*6 + Pexp + Qexp * 4 + Rexp * 4 + Sexp},
  {143872, FIVEexp*3 + Pexp * 5 + Rexp * 5 + Sexp},
  {-128, FIVEexp*7 + Pexp * 2 + Qexp * 2 + Rexp * 5 + Sexp},
  {-8192, FIVEexp*5 + Pexp * 3 + Rexp * 6 + Sexp},
  {-512, FIVEexp*6 + Qexp * 2 + Rexp * 6 + Sexp},
  {243, FIVEexp*2 + Pexp * 11 + Qexp + Sexp * 2},
  {-666, FIVEexp*3 + Pexp * 8 + Qexp * 3 + Sexp * 2},
  {-2052, FIVEexp*4 + Pexp * 5 + Qexp * 5 + Sexp * 2},
  {-916, FIVEexp*5 + Pexp * 2 + Qexp * 7 + Sexp * 2},
  {28971, FIVEexp*2 + Pexp * 9 + Qexp + Rexp + Sexp * 2},
  {78458, FIVEexp*3 + Pexp * 6 + Qexp * 3 + Rexp + Sexp * 2},
  {45399, FIVEexp*4 + Pexp * 3 + Qexp * 5 + Rexp + Sexp * 2},
  {7108, FIVEexp*5 + Qexp * 7 + Rexp + Sexp * 2},
  {-71856, FIVEexp*3 + Pexp * 7 + Qexp + Rexp * 2 + Sexp * 2},
  {-12672, FIVEexp*5 + Pexp * 4 + Qexp * 3 + Rexp * 2 + Sexp * 2},
  {-19759, FIVEexp*5 + Pexp + Qexp * 5 + Rexp * 2 + Sexp * 2},
  {-1616, FIVEexp*4 + Pexp * 5 + Qexp + Rexp * 3 + Sexp * 2},
  {-64, FIVEexp*6 + Pexp * 2 + Qexp * 3 + Rexp * 3 + Sexp * 2},
  {4992, FIVEexp*6 + Pexp * 3 + Qexp + Rexp * 4 + Sexp * 2},
  {384, FIVEexp*7 + Qexp * 3 + Rexp * 4 + Sexp * 2},
  {1024, FIVEexp*7 + Pexp + Qexp + Rexp * 5 + Sexp * 2},
  {-243, FIVEexp*5 + Pexp * 10 + Sexp * 3},
  {-3132, FIVEexp*5 + Pexp * 7 + Qexp * 2 + Sexp * 3},
  {-4, FIVEexp*10 + Pexp * 4 + Qexp * 4 + Sexp * 3},
  {-134, FIVEexp*8 + Pexp + Qexp * 6 + Sexp * 3},
  {19683, FIVEexp*4 + Pexp * 8 + Rexp + Sexp * 3},
  {31416, FIVEexp*5 + Pexp * 5 + Qexp * 2 + Rexp + Sexp * 3},
  {2881, FIVEexp*7 + Pexp * 2 + Qexp * 4 + Rexp + Sexp * 3},
  {-17568, FIVEexp*5 + Pexp * 6 + Rexp * 2 + Sexp * 3},
  {-3968, FIVEexp*7 + Pexp * 3 + Qexp * 2 + Rexp * 2 + Sexp * 3},
  {-101, FIVEexp*7 + Qexp * 4 + Rexp * 2 + Sexp * 3},
  {656, FIVEexp*7 + Pexp * 4 + Rexp * 3 + Sexp * 3},
  {-5376, FIVEexp*7 + Pexp + Qexp * 2 + Rexp * 3 + Sexp * 3},
  {1408, FIVEexp*7 + Pexp * 2 + Rexp * 4 + Sexp * 3},
  {-512, FIVEexp*8 + Rexp * 5 + Sexp * 3},
  {27, FIVEexp*7 + Pexp * 6 + Qexp + Sexp * 4},
  {-54, FIVEexp*8 + Pexp * 3 + Qexp * 3 + Sexp * 4},
  {-46, FIVEexp*9 + Qexp * 5 + Sexp * 4},
  {2334, FIVEexp*7 + Pexp * 4 + Qexp + Rexp + Sexp * 4},
  {1877, FIVEexp*8 + Pexp + Qexp * 3 + Rexp + Sexp * 4},
  {-504, FIVEexp*8 + Pexp * 2 + Qexp + Rexp * 2 + Sexp * 4},
  {576, FIVEexp*9 + Qexp + Rexp * 3 + Sexp * 4},
  {-81, FIVEexp*9 + Pexp * 5 + Sexp * 5},
  {-58, FIVEexp*10 + Pexp * 2 + Qexp * 2 + Sexp * 5},
  {52, FIVEexp*9 + Pexp * 3 + Rexp + Sexp * 5},
  {-171, FIVEexp*10 + Qexp * 2 + Rexp + Sexp * 5},
  {128, FIVEexp*10 + Pexp + Rexp * 2 + Sexp * 5},
  {-1, FIVEexp*13 + Pexp + Qexp + Sexp * 6},
  {1, FIVEexp*14 + Sexp * 7},
  {0, END_COEFF},

  // a[1]
  {8, FIVEexp*3 + Pexp * 5 + Qexp * 7},
  {58, FIVEexp*3 + Pexp * 2 + Qexp * 9},
  {-432, FIVEexp*2 + Pexp * 6 + Qexp * 5 + Rexp},
  {-3876, FIVEexp*2 + Pexp * 3 + Qexp * 7 + Rexp},
  {-84, FIVEexp*4 + Qexp * 9 + Rexp},
  {1496, FIVEexp*2 + Pexp * 7 + Qexp * 3 + Rexp * 2},
  {18834, FIVEexp*2 + Pexp * 4 + Qexp * 5 + Rexp * 2},
  {25624, FIVEexp*2 + Pexp + Qexp * 7 + Rexp * 2},
  {-1584, FIVEexp*2 + Pexp * 8 + Qexp + Rexp * 3},
  {-39344, FIVEexp*2 + Pexp * 5 + Qexp * 3 + Rexp * 3},
  {-113924, FIVEexp*2 + Pexp * 2 + Qexp * 5 + Rexp * 3},
  {32576, FIVEexp*2 + Pexp * 6 + Qexp + Rexp * 4},
  {48608, FIVEexp*3 + Pexp * 3 + Qexp * 3 + Rexp * 4},
  {18464, FIVEexp*3 + Qexp * 5 + Rexp * 4},
  {-40192, FIVEexp*3 + Pexp * 4 + Qexp + Rexp * 5},
  {-15488, FIVEexp*4 + Pexp + Qexp * 3 + Rexp * 5},
  {3072, FIVEexp*5 + Pexp * 2 + Qexp + Rexp * 6},
  {552, FIVEexp*2 + Pexp * 7 + Qexp * 4 + Sexp},
  {3786, FIVEexp*2 + Pexp * 4 + Qexp * 6 + Sexp},
  {-212, FIVEexp*3 + Pexp + Qexp * 8 + Sexp},
  {-3456, FIVEexp*2 + Pexp * 8 + Qexp * 2 + Rexp + Sexp},
  {-6532, FIVEexp*3 + Pexp * 5 + Qexp * 4 + Rexp + Sexp},
  {-412, FIVEexp*4 + Pexp * 2 + Qexp * 6 + Rexp + Sexp},
  {3672, FIVEexp*2 + Pexp * 9 + Rexp * 2 + Sexp},
  {74148, FIVEexp*2 + Pexp * 6 + Qexp * 2 + Rexp * 2 + Sexp},
  {1008, FIVEexp*4 + Pexp * 3 + Qexp * 4 + Rexp * 2 + Sexp},
  {-14354, FIVEexp*4 + Qexp * 6 + Rexp * 2 + Sexp},
  {-82848, FIVEexp*2 + Pexp * 7 + Rexp * 3 + Sexp},
  {-11584, FIVEexp*4 + Pexp * 4 + Qexp * 2 + Rexp * 3 + Sexp},
  {376, FIVEexp*7 + Pexp + Qexp * 4 + Rexp * 3 + Sexp},
  {115328, FIVEexp*3 + Pexp * 5 + Rexp * 4 + Sexp},
  {-1664, FIVEexp*5 + Pexp * 2 + Qexp * 2 + Rexp * 4 + Sexp},
  {-9728, FIVEexp*5 + Pexp * 3 + Rexp * 5 + Sexp},
  {-768, FIVEexp*6 + Qexp * 2 + Rexp * 5 + Sexp},
  {2592, FIVEexp*2 + Pexp * 9 + Qexp + Sexp * 2},
  {4536, FIVEexp*3 + Pexp * 6 + Qexp * 3 + Sexp * 2},
  {2648, FIVEexp*4 + Pexp * 3 + Qexp * 5 + Sexp * 2},
  {2236, FIVEexp*5 + Qexp * 7 + Sexp * 2},
  {108, FIVEexp*5 + Pexp * 7 + Qexp + Rexp + Sexp * 2},
  {2708, FIVEexp*5 + Pexp * 4 + Qexp * 3 + Rexp + Sexp * 2},
  {-372, FIVEexp*6 + Pexp + Qexp * 5 + Rexp + Sexp * 2},
  {-39888, FIVEexp*4 + Pexp * 5 + Qexp + Rexp * 2 + Sexp * 2},
  {-4424, FIVEexp*6 + Pexp * 2 + Qexp * 3 + Rexp * 2 + Sexp * 2},
  {6624, FIVEexp*6 + Pexp * 3 + Qexp + Rexp * 3 + Sexp * 2},
  {384, FIVEexp*7 + Qexp * 3 + Rexp * 3 + Sexp * 2},
  {1152, FIVEexp*7 + Pexp + Qexp + Rexp * 4 + Sexp * 2},
  {-1134, FIVEexp*4 + Pexp * 8 + Sexp * 3},
  {-1728, FIVEexp*5 + Pexp * 5 + Qexp * 2 + Sexp * 3},
  {114, FIVEexp*7 + Pexp * 2 + Qexp * 4 + Sexp * 3},
  {1188, FIVEexp*6 + Pexp * 6 + Rexp + Sexp * 3},
  {-8, FIVEexp*7 + Pexp * 3 + Qexp * 2 + Rexp + Sexp * 3},
  {76, FIVEexp*8 + Qexp * 4 + Rexp + Sexp * 3},
  {-192, FIVEexp*8 + Pexp * 4 + Rexp * 2 + Sexp * 3},
  {-5328, FIVEexp*7 + Pexp + Qexp * 2 + Rexp * 2 + Sexp * 3},
  {768, FIVEexp*7 + Pexp * 2 + Rexp * 3 + Sexp * 3},
  {-768, FIVEexp*8 + Rexp * 4 + Sexp * 3},
  {918, FIVEexp*7 + Pexp * 4 + Qexp + Sexp * 4},
  {484, FIVEexp*8 + Pexp + Qexp * 3 + Sexp * 4},
  {108, FIVEexp*9 + Pexp * 2 + Qexp + Rexp + Sexp * 4},
  {608, FIVEexp*9 + Qexp + Rexp * 2 + Sexp * 4},
  {-96, FIVEexp*9 + Pexp * 3 + Sexp * 5},
  {-82, FIVEexp*10 + Qexp * 2 + Sexp * 5},
  {-8, FIVEexp*11 + Pexp + Rexp + Sexp * 5},
  {0, END_COEFF},

  // a[2]
  {-4, FIVEexp*3 + Pexp * 6 + Qexp * 5},
  {-254, FIVEexp*2 + Pexp * 3 + Qexp * 7},
  {-792, FIVEexp*2 + Qexp * 9},
  {6, FIVEexp*4 + Pexp * 7 + Qexp * 3 + Rexp},
  {2604, FIVEexp*2 + Pexp * 4 + Qexp * 5 + Rexp},
  {10598, FIVEexp*2 + Pexp + Qexp * 7 + Rexp},
  {-54, FIVEexp*3 + Pexp * 8 + Qexp + Rexp * 2},
  {-8362, FIVEexp*2 + Pexp * 5 + Qexp * 3 + Rexp * 2},
  {-9738, FIVEexp*3 + Pexp * 2 + Qexp * 5 + Rexp * 2},
  {1752, FIVEexp*3 + Pexp * 6 + Qexp + Rexp * 3},
  {4016, FIVEexp*4 + Pexp * 3 + Qexp * 3 + Rexp * 3},
  {8788, FIVEexp*3 + Qexp * 5 + Rexp * 3},
  {-16544, FIVEexp*3 + Pexp * 4 + Qexp + Rexp * 4},
  {-8096, FIVEexp*4 + Pexp + Qexp * 3 + Rexp * 4},
  {1664, FIVEexp*5 + Pexp * 2 + Qexp + Rexp * 5},
  {-54, FIVEexp*3 + Pexp * 8 + Qexp * 2 + Sexp},
  {-3854, FIVEexp*2 + Pexp * 5 + Qexp * 4 + Sexp},
  {-2768, FIVEexp*3 + Pexp * 2 + Qexp * 6 + Sexp},
  {162, FIVEexp*3 + Pexp * 9 + Rexp + Sexp},
  {18396, FIVEexp*2 + Pexp * 6 + Qexp * 2 + Rexp + Sexp},
  {2926, FIVEexp*4 + Pexp * 3 + Qexp * 4 + Rexp + Sexp},
  {-4688, FIVEexp*4 + Qexp * 6 + Rexp + Sexp},
  {-4752, FIVEexp*3 + Pexp * 7 + Rexp * 2 + Sexp},
  {-6882, FIVEexp*4 + Pexp * 4 + Qexp * 2 + Rexp * 2 + Sexp},
  {698, FIVEexp*6 + Pexp + Qexp * 4 + Rexp * 2 + Sexp},
  {42016, FIVEexp*3 + Pexp * 5 + Rexp * 3 + Sexp},
  {-464, FIVEexp*5 + Pexp * 2 + Qexp * 2 + Rexp * 3 + Sexp},
  {-4096, FIVEexp*5 + Pexp * 3 + Rexp * 4 + Sexp},
  {-416, FIVEexp*6 + Qexp * 2 + Rexp * 4 + Sexp},
  {594, FIVEexp*3 + Pexp * 7 + Qexp + Sexp * 2},
  {454, FIVEexp*5 + Pexp * 4 + Qexp * 3 + Sexp * 2},
  {1906, FIVEexp*5 + Pexp + Qexp * 5 + Sexp * 2},
  {-6876, FIVEexp*4 + Pexp * 5 + Qexp + Rexp + Sexp * 2},
  {-1914, FIVEexp*6 + Pexp * 2 + Qexp * 3 + Rexp + Sexp * 2},
  {2016, FIVEexp*6 + Pexp * 3 + Qexp + Rexp * 2 + Sexp * 2},
  {32, FIVEexp*8 + Qexp * 3 + Rexp * 2 + Sexp * 2},
  {448, FIVEexp*7 + Pexp + Qexp + Rexp * 3 + Sexp * 2},
  {432, FIVEexp*5 + Pexp * 6 + Sexp * 3},
  {78, FIVEexp*7 + Pexp * 3 + Qexp * 2 + Sexp * 3},
  {224, FIVEexp*7 + Qexp * 4 + Sexp * 3},
  {-18, FIVEexp*8 + Pexp * 4 + Rexp + Sexp * 3},
  {-1636, FIVEexp*7 + Pexp + Qexp * 2 + Rexp + Sexp * 3},
  {48, FIVEexp*8 + Pexp * 2 + Rexp * 2 + Sexp * 3},
  {-416, FIVEexp*8 + Rexp * 3 + Sexp * 3},
  {276, FIVEexp*8 + Pexp * 2 + Qexp + Sexp * 4},
  {236, FIVEexp*9 + Qexp + Rexp + Sexp * 4},
  {-22, FIVEexp*10 + Pexp + Sexp * 5},
  {0, END_COEFF},

  // a[3]
  {78, FIVEexp*2 + Pexp * 4 + Qexp * 5},
  {564, FIVEexp*2 + Pexp + Qexp * 7},
  {-574, FIVEexp*2 + Pexp * 5 + Qexp * 3 + Rexp},
  {-5024, FIVEexp*2 + Pexp * 2 + Qexp * 5 + Rexp},
  {1116, FIVEexp*2 + Pexp * 6 + Qexp + Rexp * 2},
  {3218, FIVEexp*3 + Pexp * 3 + Qexp * 3 + Rexp * 2},
  {2306, FIVEexp*3 + Qexp * 5 + Rexp * 2},
  {-3488, FIVEexp*3 + Pexp * 4 + Qexp + Rexp * 3},
  {-2152, FIVEexp*4 + Pexp + Qexp * 3 + Rexp * 3},
  {448, FIVEexp*5 + Pexp * 2 + Qexp + Rexp * 4},
  {378, FIVEexp*2 + Pexp * 6 + Qexp * 2 + Sexp},
  {-2, FIVEexp*4 + Pexp * 3 + Qexp * 4 + Sexp},
  {-744, FIVEexp*4 + Qexp * 6 + Sexp},
  {-1998, FIVEexp*2 + Pexp * 7 + Rexp + Sexp},
  {-484, FIVEexp*4 + Pexp * 4 + Qexp * 2 + Rexp + Sexp},
  {22, FIVEexp*7 + Pexp + Qexp * 4 + Rexp + Sexp},
  {6672, FIVEexp*3 + Pexp * 5 + Rexp * 2 + Sexp},
  {28, FIVEexp*6 + Pexp * 2 + Qexp * 2 + Rexp * 2 + Sexp},
  {-992, FIVEexp*5 + Pexp * 3 + Rexp * 3 + Sexp},
  {-112, FIVEexp*6 + Qexp * 2 + Rexp * 3 + Sexp},
  {-468, FIVEexp*4 + Pexp * 5 + Qexp + Sexp * 2},
  {-124, FIVEexp*6 + Pexp * 2 + Qexp * 3 + Sexp * 2},
  {214, FIVEexp*6 + Pexp * 3 + Qexp + Rexp + Sexp * 2},
  {24, FIVEexp*7 + Qexp * 3 + Rexp + Sexp * 2},
  {104, FIVEexp*7 + Pexp + Qexp + Rexp * 2 + Sexp * 2},
  {-18, FIVEexp*7 + Pexp * 4 + Sexp * 3},
  {-158, FIVEexp*7 + Pexp + Qexp * 2 + Sexp * 3},
  {68, FIVEexp*7 + Pexp * 2 + Rexp + Sexp * 3},
  {-112, FIVEexp*8 + Rexp * 2 + Sexp * 3},
  {38, FIVEexp*9 + Qexp + Sexp * 4},
  {0, END_COEFF},

  // a[4]
  {-12, FIVEexp*2 + Pexp * 5 + Qexp * 3},
  {-86, FIVEexp*2 + Pexp * 2 + Qexp * 5},
  {54, FIVEexp*2 + Pexp * 6 + Qexp + Rexp},
  {172, FIVEexp*3 + Pexp * 3 + Qexp * 3 + Rexp},
  {492, FIVEexp*3 + Qexp * 5 + Rexp},
  {-336, FIVEexp*3 + Pexp * 4 + Qexp + Rexp * 2},
  {-464, FIVEexp*4 + Pexp + Qexp * 3 + Rexp * 2},
  {96, FIVEexp*5 + Pexp * 2 + Qexp + Rexp * 3},
  {-162, FIVEexp*2 + Pexp * 7 + Sexp},
  {-72, FIVEexp*4 + Pexp * 4 + Qexp * 2 + Sexp},
  {-8, FIVEexp*6 + Pexp + Qexp * 4 + Sexp},
  {864, FIVEexp*3 + Pexp * 5 + Rexp + Sexp},
  {206, FIVEexp*5 + Pexp * 2 + Qexp * 2 + Rexp + Sexp},
  {-224, FIVEexp*5 + Pexp * 3 + Rexp * 2 + Sexp},
  {-24, FIVEexp*6 + Qexp * 2 + Rexp * 2 + Sexp},
  {-6, FIVEexp*6 + Pexp * 3 + Qexp + Sexp * 2},
  {-4, FIVEexp*7 + Qexp * 3 + Sexp * 2},
  {24, FIVEexp*7 + Pexp + Qexp + Rexp + Sexp * 2},
  {-18, FIVEexp*7 + Pexp * 2 + Sexp * 3},
  {-24, FIVEexp*8 + Rexp + Sexp * 3},
  {0, END_COEFF},

  // a[5]
  {2, FIVEexp*4 + Pexp * 3 + Qexp * 3},
  {72, FIVEexp*3 + Qexp * 5},
  {-36, FIVEexp*3 + Pexp * 4 + Qexp + Rexp},
  {-74, FIVEexp*4 + Pexp + Qexp * 3 + Rexp},
  {16, FIVEexp*5 + Pexp * 2 + Qexp + Rexp * 2},
  {54, FIVEexp*3 + Pexp * 5 + Sexp},
  {14, FIVEexp*5 + Pexp * 2 + Qexp * 2 + Sexp},
  {-24, FIVEexp*5 + Pexp * 3 + Rexp + Sexp},
  {-4, FIVEexp*6 + Qexp * 2 + Rexp + Sexp},
  {2, FIVEexp*7 + Pexp + Qexp + Sexp * 2},
  {-4, FIVEexp*8 + Sexp * 3},
  {0, END_ARRAY},

  // c[0]
  {-8, FIVEexp*1 + Pexp * 5 + Qexp * 11},
  {-54, FIVEexp*1 + Pexp * 2 + Qexp * 13},
  {28, FIVEexp*2 + Pexp * 6 + Qexp * 9 + Rexp},
  {1033, FIVEexp*1 + Pexp * 3 + Qexp * 11 + Rexp},
  {108, FIVEexp*1 + Qexp * 13 + Rexp},
  {-846, FIVEexp*1 + Pexp * 7 + Qexp * 7 + Rexp * 2},
  {-6369, FIVEexp*1 + Pexp * 4 + Qexp * 9 + Rexp * 2},
  {4176, FIVEexp*1 + Pexp + Qexp * 11 + Rexp * 2},
  {1929, FIVEexp*1 + Pexp * 8 + Qexp * 5 + Rexp * 3},
  {11523, FIVEexp*1 + Pexp * 5 + Qexp * 7 + Rexp * 3},
  {-71651, FIVEexp*1 + Pexp * 2 + Qexp * 9 + Rexp * 3},
  {-376, FIVEexp*1 + Pexp * 9 + Qexp * 3 + Rexp * 4},
  {22804, FIVEexp*1 + Pexp * 6 + Qexp * 5 + Rexp * 4},
  {402438, FIVEexp*1 + Pexp * 3 + Qexp * 7 + Rexp * 4},
  {-5371, FIVEexp*1 + Qexp * 9 + Rexp * 4},
  {-576, FIVEexp*2 + Pexp * 10 + Qexp + Rexp * 5},
  {-18816, FIVEexp*2 + Pexp * 7 + Qexp * 3 + Rexp * 5},
  {-1017728, FIVEexp*1 + Pexp * 4 + Qexp * 5 + Rexp * 5},
  {184, FIVEexp*1 + Pexp + Qexp * 7 + Rexp * 5},
  {13312, FIVEexp*2 + Pexp * 8 + Qexp + Rexp * 6},
  {1159424, FIVEexp*1 + Pexp * 5 + Qexp * 3 + Rexp * 6},
  {12864, FIVEexp*3 + Pexp * 2 + Qexp * 5 + Rexp * 6},
  {-104448, FIVEexp*2 + Pexp * 6 + Qexp + Rexp * 7},
  {-59392, FIVEexp*3 + Pexp * 3 + Qexp * 3 + Rexp * 7},
  {-92928, FIVEexp*2 + Qexp * 5 + Rexp * 7},
  {344064, FIVEexp*2 + Pexp * 4 + Qexp + Rexp * 8},
  {75776, FIVEexp*3 + Pexp + Qexp * 3 + Rexp * 8},
  {-16384, FIVEexp*4 + Pexp * 2 + Qexp + Rexp * 9},
  {-612, FIVEexp*1 + Pexp * 7 + Qexp * 8 + Sexp},
  {-7817, FIVEexp*1 + Pexp * 4 + Qexp * 10 + Sexp},
  {-5292, FIVEexp*2 + Pexp + Qexp * 12 + Sexp},
  {7316, FIVEexp*1 + Pexp * 8 + Qexp * 6 + Rexp + Sexp},
  {104037, FIVEexp*1 + Pexp * 5 + Qexp * 8 + Rexp + Sexp},
  {393972, FIVEexp*1 + Pexp * 2 + Qexp * 10 + Rexp + Sexp},
  {-28809, FIVEexp*1 + Pexp * 9 + Qexp * 4 + Rexp * 2 + Sexp},
  {-97537, FIVEexp*2 + Pexp * 6 + Qexp * 6 + Rexp * 2 + Sexp},
  {-432379, FIVEexp*2 + Pexp * 3 + Qexp * 8 + Rexp * 2 + Sexp},
  {20754, FIVEexp*2 + Qexp * 10 + Rexp * 2 + Sexp},
  {36504, FIVEexp*1 + Pexp * 10 + Qexp * 2 + Rexp * 3 + Sexp},
  {906786, FIVEexp*1 + Pexp * 7 + Qexp * 4 + Rexp * 3 + Sexp},
  {5239354, FIVEexp*1 + Pexp * 4 + Qexp * 6 + Rexp * 3 + Sexp},
  {-181693, FIVEexp*2 + Pexp + Qexp * 8 + Rexp * 3 + Sexp},
  {864, FIVEexp*2 + Pexp * 11 + Rexp * 4 + Sexp},
  {-441616, FIVEexp*1 + Pexp * 8 + Qexp * 2 + Rexp * 4 + Sexp},
  {-4957592, FIVEexp*1 + Pexp * 5 + Qexp * 4 + Rexp * 4 + Sexp},
  {432556, FIVEexp*2 + Pexp * 2 + Qexp * 6 + Rexp * 4 + Sexp},
  {-19968, FIVEexp*2 + Pexp * 9 + Rexp * 5 + Sexp},
  {765568, FIVEexp*1 + Pexp * 6 + Qexp * 2 + Rexp * 5 + Sexp},
  {76768, FIVEexp*3 + Pexp * 3 + Qexp * 4 + Rexp * 5 + Sexp},
  {181296, FIVEexp*3 + Qexp * 6 + Rexp * 5 + Sexp},
  {156672, FIVEexp*2 + Pexp * 7 + Rexp * 6 + Sexp},
  {-239616, FIVEexp*3 + Pexp * 4 + Qexp * 2 + Rexp * 6 + Sexp},
  {-29056, FIVEexp*5 + Pexp + Qexp * 4 + Rexp * 6 + Sexp},
  {-516096, FIVEexp*2 + Pexp * 5 + Rexp * 7 + Sexp},
  {139264, FIVEexp*4 + Pexp * 2 + Qexp * 2 + Rexp * 7 + Sexp},
  {24576, FIVEexp*4 + Pexp * 3 + Rexp * 8 + Sexp},
  {4096, FIVEexp*5 + Qexp * 2 + Rexp * 8 + Sexp},
  {-7614, FIVEexp*1 + Pexp * 9 + Qexp * 5 + Sexp * 2},
  {-22668, FIVEexp*2 + Pexp * 6 + Qexp * 7 + Sexp * 2},
  {-4119, FIVEexp*4 + Pexp * 3 + Qexp * 9 + Sexp * 2},
  {-2916, FIVEexp*4 + Qexp * 11 + Sexp * 2},
  {58563, FIVEexp*1 + Pexp * 10 + Qexp * 3 + Rexp + Sexp * 2},
  {1034056, FIVEexp*1 + Pexp * 7 + Qexp * 5 + Rexp + Sexp * 2},
  {44669, FIVEexp*4 + Pexp * 4 + Qexp * 7 + Rexp + Sexp * 2},
  {35196, FIVEexp*4 + Pexp + Qexp * 9 + Rexp + Sexp * 2},
  {-114696, FIVEexp*1 + Pexp * 11 + Qexp + Rexp * 2 + Sexp * 2},
  {-582654, FIVEexp*2 + Pexp * 8 + Qexp * 3 + Rexp * 2 + Sexp * 2},
  {-4194063, FIVEexp*2 + Pexp * 5 + Qexp * 5 + Rexp * 2 + Sexp * 2},
  {-154318, FIVEexp*4 + Pexp * 2 + Qexp * 7 + Rexp * 2 + Sexp * 2},
  {2200248, FIVEexp*1 + Pexp * 9 + Qexp + Rexp * 3 + Sexp * 2},
  {5911944, FIVEexp*2 + Pexp * 6 + Qexp * 3 + Rexp * 3 + Sexp * 2},
  {1269062, FIVEexp*3 + Pexp * 3 + Qexp * 5 + Rexp * 3 + Sexp * 2},
  {-125156, FIVEexp*4 + Qexp * 7 + Rexp * 3 + Sexp * 2},
  {-2512768, FIVEexp*2 + Pexp * 7 + Qexp + Rexp * 4 + Sexp * 2},
  {-217856, FIVEexp*4 + Pexp * 4 + Qexp * 3 + Rexp * 4 + Sexp * 2},
  {508088, FIVEexp*4 + Pexp + Qexp * 5 + Rexp * 4 + Sexp * 2},
  {1281792, FIVEexp*3 + Pexp * 5 + Qexp + Rexp * 5 + Sexp * 2},
  {-85632, FIVEexp*5 + Pexp * 2 + Qexp * 3 + Rexp * 5 + Sexp * 2},
  {-49152, FIVEexp*5 + Pexp * 3 + Qexp + Rexp * 6 + Sexp * 2},
  {-1536, FIVEexp*7 + Qexp * 3 + Rexp * 6 + Sexp * 2},
  {-2048, FIVEexp*6 + Pexp + Qexp + Rexp * 7 + Sexp * 2},
  {-5103, FIVEexp*2 + Pexp * 11 + Qexp * 2 + Sexp * 3},
  {-3438, FIVEexp*4 + Pexp * 8 + Qexp * 4 + Sexp * 3},
  {-21844, FIVEexp*4 + Pexp * 5 + Qexp * 6 + Sexp * 3},
  {-2, FIVEexp*10 + Pexp * 2 + Qexp * 8 + Sexp * 3},
  {99144, FIVEexp*1 + Pexp * 12 + Rexp + Sexp * 3},
  {94851, FIVEexp*3 + Pexp * 9 + Qexp * 2 + Rexp + Sexp * 3},
  {172492, FIVEexp*4 + Pexp * 6 + Qexp * 4 + Rexp + Sexp * 3},
  {71147, FIVEexp*5 + Pexp * 3 + Qexp * 6 + Rexp + Sexp * 3},
  {6724, FIVEexp*6 + Qexp * 8 + Rexp + Sexp * 3},
  {-462672, FIVEexp*2 + Pexp * 10 + Rexp * 2 + Sexp * 3},
  {-69372, FIVEexp*5 + Pexp * 7 + Qexp * 2 + Rexp * 2 + Sexp * 3},
  {-8108, FIVEexp*7 + Pexp * 4 + Qexp * 4 + Rexp * 2 + Sexp * 3},
  {-32287, FIVEexp*6 + Pexp + Qexp * 6 + Rexp * 2 + Sexp * 3},
  {727344, FIVEexp*3 + Pexp * 8 + Rexp * 3 + Sexp * 3},
  {907328, FIVEexp*4 + Pexp * 5 + Qexp * 2 + Rexp * 3 + Sexp * 3},
  {44348, FIVEexp*6 + Pexp * 2 + Qexp * 4 + Rexp * 3 + Sexp * 3},
  {-522624, FIVEexp*4 + Pexp * 6 + Rexp * 4 + Sexp * 3},
  {-21696, FIVEexp*6 + Pexp * 3 + Qexp * 2 + Rexp * 4 + Sexp * 3},
  {23632, FIVEexp*6 + Qexp * 4 + Rexp * 4 + Sexp * 3},
  {7168, FIVEexp*7 + Pexp * 4 + Rexp * 5 + Sexp * 3},
  {32512, FIVEexp*6 + Pexp + Qexp * 2 + Rexp * 5 + Sexp * 3},
  {-6144, FIVEexp*7 + Pexp * 2 + Rexp * 6 + Sexp * 3},
  {4096, FIVEexp*7 + Rexp * 7 + Sexp * 3},
  {-729, FIVEexp*4 + Pexp * 10 + Qexp + Sexp * 4},
  {-1764, FIVEexp*6 + Pexp * 7 + Qexp * 3 + Sexp * 4},
  {-7718, FIVEexp*6 + Pexp * 4 + Qexp * 5 + Sexp * 4},
  {-772, FIVEexp*7 + Pexp + Qexp * 7 + Sexp * 4},
  {35397, FIVEexp*5 + Pexp * 8 + Qexp + Rexp + Sexp * 4},
  {42431, FIVEexp*6 + Pexp * 5 + Qexp * 3 + Rexp + Sexp * 4},
  {1353, FIVEexp*8 + Pexp * 2 + Qexp * 5 + Rexp + Sexp * 4},
  {-34668, FIVEexp*6 + Pexp * 6 + Qexp + Rexp * 2 + Sexp * 4},
  {-16158, FIVEexp*7 + Pexp * 3 + Qexp * 3 + Rexp * 2 + Sexp * 4},
  {-1194, FIVEexp*8 + Qexp * 5 + Rexp * 2 + Sexp * 4},
  {40512, FIVEexp*6 + Pexp * 4 + Qexp + Rexp * 3 + Sexp * 4},
  {-16184, FIVEexp*7 + Pexp + Qexp * 3 + Rexp * 3 + Sexp * 4},
  {13888, FIVEexp*7 + Pexp * 2 + Qexp + Rexp * 4 + Sexp * 4},
  {-6912, FIVEexp*8 + Qexp + Rexp * 5 + Sexp * 4},
  {-4374, FIVEexp*6 + Pexp * 9 + Sexp * 5},
  {-6129, FIVEexp*7 + Pexp * 6 + Qexp * 2 + Sexp * 5},
  {-182, FIVEexp*9 + Pexp * 3 + Qexp * 4 + Sexp * 5},
  {-6, FIVEexp*9 + Qexp * 6 + Sexp * 5},
  {45954, FIVEexp*6 + Pexp * 7 + Rexp + Sexp * 5},
  {4246, FIVEexp*8 + Pexp * 4 + Qexp * 2 + Rexp + Sexp * 5},
  {1133, FIVEexp*9 + Pexp + Qexp * 4 + Rexp + Sexp * 5},
  {-36552, FIVEexp*7 + Pexp * 5 + Rexp * 2 + Sexp * 5},
  {-2188, FIVEexp*9 + Pexp * 2 + Qexp * 2 + Rexp * 2 + Sexp * 5},
  {11616, FIVEexp*8 + Pexp * 3 + Rexp * 3 + Sexp * 5},
  {4112, FIVEexp*9 + Qexp * 2 + Rexp * 3 + Sexp * 5},
  {-896, FIVEexp*9 + Pexp + Rexp * 4 + Sexp * 5},
  {693, FIVEexp*9 + Pexp * 5 + Qexp + Sexp * 6},
  {32, FIVEexp*11 + Pexp * 2 + Qexp * 3 + Sexp * 6},
  {-406, FIVEexp*10 + Pexp * 3 + Qexp + Rexp + Sexp * 6},
  {-148, FIVEexp*11 + Qexp * 3 + Rexp + Sexp * 6},
  {8, FIVEexp*12 + Pexp + Qexp + Rexp * 2 + Sexp * 6},
  {-36, FIVEexp*11 + Pexp * 4 + Sexp * 7},
  {-13, FIVEexp*12 + Pexp + Qexp * 2 + Sexp * 7},
  {132, FIVEexp*11 + Pexp * 2 + Rexp + Sexp * 7},
  {-16, FIVEexp*12 + Rexp * 2 + Sexp * 7},
  {1, FIVEexp*14 + Qexp + Sexp * 8},
  {0, END_COEFF},

  // c[1]
  {8, FIVEexp*1 + Pexp * 6 + Qexp * 9},
  {22, FIVEexp*1 + Pexp * 3 + Qexp * 11},
  {-216, FIVEexp*1 + Qexp * 13},
  {-112, FIVEexp*1 + Pexp * 7 + Qexp * 7 + Rexp},
  {-356, FIVEexp*1 + Pexp * 4 + Qexp * 9 + Rexp},
  {3474, FIVEexp*1 + Pexp + Qexp * 11 + Rexp},
  {114, FIVEexp*2 + Pexp * 8 + Qexp * 5 + Rexp * 2},
  {2104, FIVEexp*1 + Pexp * 5 + Qexp * 7 + Rexp * 2},
  {-23182, FIVEexp*1 + Pexp * 2 + Qexp * 9 + Rexp * 2},
  {-1218, FIVEexp*1 + Pexp * 9 + Qexp * 3 + Rexp * 3},
  {-5066, FIVEexp*1 + Pexp * 6 + Qexp * 5 + Rexp * 3},
  {89748, FIVEexp*1 + Pexp * 3 + Qexp * 7 + Rexp * 3},
  {25646, FIVEexp*1 + Qexp * 9 + Rexp * 3},
  {864, FIVEexp*1 + Pexp * 10 + Qexp + Rexp * 4},
  {3392, FIVEexp*1 + Pexp * 7 + Qexp * 3 + Rexp * 4},
  {-45744, FIVEexp*2 + Pexp * 4 + Qexp * 5 + Rexp * 4},
  {-282062, FIVEexp*1 + Pexp + Qexp * 7 + Rexp * 4},
  {768, FIVEexp*1 + Pexp * 8 + Qexp + Rexp * 5},
  {348896, FIVEexp*1 + Pexp * 5 + Qexp * 3 + Rexp * 5},
  {1123904, FIVEexp*1 + Pexp * 2 + Qexp * 5 + Rexp * 5},
  {-239616, FIVEexp*1 + Pexp * 6 + Qexp + Rexp * 6},
  {-423168, FIVEexp*2 + Pexp * 3 + Qexp * 3 + Rexp * 6},
  {-117632, FIVEexp*2 + Qexp * 5 + Rexp * 6},
  {331776, FIVEexp*2 + Pexp * 4 + Qexp + Rexp * 7},
  {108544, FIVEexp*3 + Pexp + Qexp * 3 + Rexp * 7},
  {-24576, FIVEexp*4 + Pexp * 2 + Qexp + Rexp * 8},
  {168, FIVEexp*1 + Pexp * 8 + Qexp * 6 + Sexp},
  {1516, FIVEexp*1 + Pexp * 5 + Qexp * 8 + Sexp},
  {4884, FIVEexp*1 + Pexp * 2 + Qexp * 10 + Sexp},
  {-324, FIVEexp*2 + Pexp * 9 + Qexp * 4 + Rexp + Sexp},
  {-3764, FIVEexp*2 + Pexp * 6 + Qexp * 6 + Rexp + Sexp},
  {-3784, FIVEexp*3 + Pexp * 3 + Qexp * 8 + Rexp + Sexp},
  {-18936, FIVEexp*2 + Qexp * 10 + Rexp + Sexp},
  {4536, FIVEexp*1 + Pexp * 10 + Qexp * 2 + Rexp * 2 + Sexp},
  {74874, FIVEexp*1 + Pexp * 7 + Qexp * 4 + Rexp * 2 + Sexp},
  {577604, FIVEexp*1 + Pexp * 4 + Qexp * 6 + Rexp * 2 + Sexp},
  {222442, FIVEexp*2 + Pexp + Qexp * 8 + Rexp * 2 + Sexp},
  {-2592, FIVEexp*1 + Pexp * 11 + Rexp * 3 + Sexp},
  {-97164, FIVEexp*1 + Pexp * 8 + Qexp * 2 + Rexp * 3 + Sexp},
  {-1344688, FIVEexp*1 + Pexp * 5 + Qexp * 4 + Rexp * 3 + Sexp},
  {-942456, FIVEexp*2 + Pexp * 2 + Qexp * 6 + Rexp * 3 + Sexp},
  {38016, FIVEexp*1 + Pexp * 9 + Rexp * 4 + Sexp},
  {1178976, FIVEexp*1 + Pexp * 6 + Qexp * 2 + Rexp * 4 + Sexp},
  {407056, FIVEexp*3 + Pexp * 3 + Qexp * 4 + Rexp * 4 + Sexp},
  {179292, FIVEexp*3 + Qexp * 6 + Rexp * 4 + Sexp},
  {-51712, FIVEexp*1 + Pexp * 7 + Rexp * 5 + Sexp},
  {-369984, FIVEexp*3 + Pexp * 4 + Qexp * 2 + Rexp * 5 + Sexp},
  {-33216, FIVEexp*5 + Pexp + Qexp * 4 + Rexp * 5 + Sexp},
  {-149504, FIVEexp*2 + Pexp * 5 + Rexp * 6 + Sexp},
  {191488, FIVEexp*4 + Pexp * 2 + Qexp * 2 + Rexp * 6 + Sexp},
  {16384, FIVEexp*4 + Pexp * 3 + Rexp * 7 + Sexp},
  {6144, FIVEexp*5 + Qexp * 2 + Rexp * 7 + Sexp},
  {1458, FIVEexp*1 + Pexp * 10 + Qexp * 3 + Sexp * 2},
  {23472, FIVEexp*1 + Pexp * 7 + Qexp * 5 + Sexp * 2},
  {1106, FIVEexp*4 + Pexp * 4 + Qexp * 7 + Sexp * 2},
  {-318, FIVEexp*4 + Pexp + Qexp * 9 + Sexp * 2},
  {-1458, FIVEexp*2 + Pexp * 11 + Qexp + Rexp + Sexp * 2},
  {-34182, FIVEexp*2 + Pexp * 8 + Qexp * 3 + Rexp + Sexp * 2},
  {-293628, FIVEexp*2 + Pexp * 5 + Qexp * 5 + Rexp + Sexp * 2},
  {-3246, FIVEexp*4 + Pexp * 2 + Qexp * 7 + Rexp + Sexp * 2},
  {199098, FIVEexp*1 + Pexp * 9 + Qexp + Rexp * 2 + Sexp * 2},
  {755864, FIVEexp*2 + Pexp * 6 + Qexp * 3 + Rexp * 2 + Sexp * 2},
  {40212, FIVEexp*3 + Pexp * 3 + Qexp * 5 + Rexp * 2 + Sexp * 2},
  {-83636, FIVEexp*4 + Qexp * 7 + Rexp * 2 + Sexp * 2},
  {-665472, FIVEexp*2 + Pexp * 7 + Qexp + Rexp * 3 + Sexp * 2},
  {-13824, FIVEexp*5 + Pexp * 4 + Qexp * 3 + Rexp * 3 + Sexp * 2},
  {357482, FIVEexp*4 + Pexp + Qexp * 5 + Rexp * 3 + Sexp * 2},
  {896544, FIVEexp*3 + Pexp * 5 + Qexp + Rexp * 4 + Sexp * 2},
  {-11328, FIVEexp*6 + Pexp * 2 + Qexp * 3 + Rexp * 4 + Sexp * 2},
  {-15616, FIVEexp*6 + Pexp * 3 + Qexp + Rexp * 5 + Sexp * 2},
  {-9984, FIVEexp*6 + Qexp * 3 + Rexp * 5 + Sexp * 2},
  {8748, FIVEexp*1 + Pexp * 12 + Sexp * 3},
  {8262, FIVEexp*3 + Pexp * 9 + Qexp * 2 + Sexp * 3},
  {13764, FIVEexp*4 + Pexp * 6 + Qexp * 4 + Sexp * 3},
  {4994, FIVEexp*5 + Pexp * 3 + Qexp * 6 + Sexp * 3},
  {2536, FIVEexp*6 + Qexp * 8 + Sexp * 3},
  {-64152, FIVEexp*2 + Pexp * 10 + Rexp + Sexp * 3},
  {-43092, FIVEexp*4 + Pexp * 7 + Qexp * 2 + Rexp + Sexp * 3},
  {-1252, FIVEexp*6 + Pexp * 4 + Qexp * 4 + Rexp + Sexp * 3},
  {-9728, FIVEexp*6 + Pexp + Qexp * 6 + Rexp + Sexp * 3},
  {204444, FIVEexp*3 + Pexp * 8 + Rexp * 2 + Sexp * 3},
  {25968, FIVEexp*4 + Pexp * 5 + Qexp * 2 + Rexp * 2 + Sexp * 3},
  {108, FIVEexp*8 + Pexp * 2 + Qexp * 4 + Rexp * 2 + Sexp * 3},
  {-265056, FIVEexp*4 + Pexp * 6 + Rexp * 3 + Sexp * 3},
  {23904, FIVEexp*6 + Pexp * 3 + Qexp * 2 + Rexp * 3 + Sexp * 3},
  {21308, FIVEexp*6 + Qexp * 4 + Rexp * 3 + Sexp * 3},
  {5952, FIVEexp*7 + Pexp * 4 + Rexp * 4 + Sexp * 3},
  {37504, FIVEexp*6 + Pexp + Qexp * 2 + Rexp * 4 + Sexp * 3},
  {-37888, FIVEexp*6 + Pexp * 2 + Rexp * 5 + Sexp * 3},
  {6144, FIVEexp*7 + Rexp * 6 + Sexp * 3},
  {-486, FIVEexp*5 + Pexp * 8 + Qexp + Sexp * 4},
  {-4002, FIVEexp*6 + Pexp * 5 + Qexp * 3 + Sexp * 4},
  {98, FIVEexp*7 + Pexp * 2 + Qexp * 5 + Sexp * 4},
  {11826, FIVEexp*6 + Pexp * 6 + Qexp + Rexp + Sexp * 4},
  {-202, FIVEexp*7 + Pexp * 3 + Qexp * 3 + Rexp + Sexp * 4},
  {-346, FIVEexp*8 + Qexp * 5 + Rexp + Sexp * 4},
  {-73488, FIVEexp*6 + Pexp * 4 + Qexp + Rexp * 2 + Sexp * 4},
  {-27154, FIVEexp*7 + Pexp + Qexp * 3 + Rexp * 2 + Sexp * 4},
  {25472, FIVEexp*7 + Pexp * 2 + Qexp + Rexp * 3 + Sexp * 4},
  {-8064, FIVEexp*8 + Qexp + Rexp * 4 + Sexp * 4},
  {-162, FIVEexp*6 + Pexp * 7 + Sexp * 5},
  {1692, FIVEexp*8 + Pexp * 4 + Qexp * 2 + Sexp * 5},
  {686, FIVEexp*9 + Pexp + Qexp * 4 + Sexp * 5},
  {-432, FIVEexp*7 + Pexp * 5 + Rexp + Sexp * 5},
  {-348, FIVEexp*9 + Pexp * 2 + Qexp * 2 + Rexp + Sexp * 5},
  {16, FIVEexp*8 + Pexp * 3 + Rexp * 2 + Sexp * 5},
  {3172, FIVEexp*9 + Qexp * 2 + Rexp * 2 + Sexp * 5},
  {576, FIVEexp*9 + Pexp + Rexp * 3 + Sexp * 5},
  {-102, FIVEexp*10 + Pexp * 3 + Qexp + Sexp * 6},
  {-64, FIVEexp*11 + Qexp * 3 + Sexp * 6},
  {-66, FIVEexp*11 + Pexp + Qexp + Rexp + Sexp * 6},
  {24, FIVEexp*11 + Pexp * 2 + Sexp * 7},
  {4, FIVEexp*12 + Rexp + Sexp * 7},
  {0, END_COEFF},

  // c[2]
  {16, FIVEexp*1 + Pexp * 4 + Qexp * 9},
  {108, FIVEexp*1 + Pexp + Qexp * 11},
  {-24, FIVEexp*2 + Pexp * 5 + Qexp * 7 + Rexp},
  {-954, FIVEexp*1 + Pexp * 2 + Qexp * 9 + Rexp},
  {246, FIVEexp*1 + Pexp * 6 + Qexp * 5 + Rexp * 2},
  {836, FIVEexp*2 + Pexp * 3 + Qexp * 7 + Rexp * 2},
  {378, FIVEexp*3 + Qexp * 9 + Rexp * 2},
  {-142, FIVEexp*1 + Pexp * 7 + Qexp * 3 + Rexp * 3},
  {-3398, FIVEexp*2 + Pexp * 4 + Qexp * 5 + Rexp * 3},
  {-105262, FIVEexp*1 + Pexp + Qexp * 7 + Rexp * 3},
  {144, FIVEexp*1 + Pexp * 8 + Qexp + Rexp * 4},
  {43256, FIVEexp*1 + Pexp * 5 + Qexp * 3 + Rexp * 4},
  {413604, FIVEexp*1 + Pexp * 2 + Qexp * 5 + Rexp * 4},
  {-39616, FIVEexp*1 + Pexp * 6 + Qexp + Rexp * 5},
  {-148128, FIVEexp*2 + Pexp * 3 + Qexp * 3 + Rexp * 5},
  {-56944, FIVEexp*2 + Qexp * 5 + Rexp * 5},
  {114432, FIVEexp*2 + Pexp * 4 + Qexp + Rexp * 6},
  {56448, FIVEexp*3 + Pexp + Qexp * 3 + Rexp * 6},
  {-13312, FIVEexp*4 + Pexp * 2 + Qexp + Rexp * 7},
  {-544, FIVEexp*1 + Pexp * 6 + Qexp * 6 + Sexp},
  {-1854, FIVEexp*2 + Pexp * 3 + Qexp * 8 + Sexp},
  {-7128, FIVEexp*2 + Qexp * 10 + Sexp},
  {5148, FIVEexp*1 + Pexp * 7 + Qexp * 4 + Rexp + Sexp},
  {97898, FIVEexp*1 + Pexp * 4 + Qexp * 6 + Rexp + Sexp},
  {86094, FIVEexp*2 + Pexp + Qexp * 8 + Rexp + Sexp},
  {-12312, FIVEexp*1 + Pexp * 8 + Qexp * 2 + Rexp * 2 + Sexp},
  {-62726, FIVEexp*2 + Pexp * 5 + Qexp * 4 + Rexp * 2 + Sexp},
  {-72484, FIVEexp*3 + Pexp * 2 + Qexp * 6 + Rexp * 2 + Sexp},
  {4968, FIVEexp*1 + Pexp * 9 + Rexp * 3 + Sexp},
  {338476, FIVEexp*1 + Pexp * 6 + Qexp * 2 + Rexp * 3 + Sexp},
  {144786, FIVEexp*3 + Pexp * 3 + Qexp * 4 + Rexp * 3 + Sexp},
  {75102, FIVEexp*3 + Qexp * 6 + Rexp * 3 + Sexp},
  {-76512, FIVEexp*1 + Pexp * 7 + Rexp * 4 + Sexp},
  {-134544, FIVEexp*3 + Pexp * 4 + Qexp * 2 + Rexp * 4 + Sexp},
  {-15784, FIVEexp*5 + Pexp + Qexp * 4 + Rexp * 4 + Sexp},
  {48512, FIVEexp*2 + Pexp * 5 + Rexp * 5 + Sexp},
  {103744, FIVEexp*4 + Pexp * 2 + Qexp * 2 + Rexp * 5 + Sexp},
  {-512, FIVEexp*4 + Pexp * 3 + Rexp * 6 + Sexp},
  {3328, FIVEexp*5 + Qexp * 2 + Rexp * 6 + Sexp},
  {-1458, FIVEexp*2 + Pexp * 8 + Qexp * 3 + Sexp * 2},
  {-23534, FIVEexp*2 + Pexp * 5 + Qexp * 5 + Sexp * 2},
  {-138, FIVEexp*6 + Pexp * 2 + Qexp * 7 + Sexp * 2},
  {24786, FIVEexp*1 + Pexp * 9 + Qexp + Rexp + Sexp * 2},
  {115188, FIVEexp*2 + Pexp * 6 + Qexp * 3 + Rexp + Sexp * 2},
  {100384, FIVEexp*3 + Pexp * 3 + Qexp * 5 + Rexp + Sexp * 2},
  {-23112, FIVEexp*4 + Qexp * 7 + Rexp + Sexp * 2},
  {-25866, FIVEexp*3 + Pexp * 7 + Qexp + Rexp * 2 + Sexp * 2},
  {-45576, FIVEexp*4 + Pexp * 4 + Qexp * 3 + Rexp * 2 + Sexp * 2},
  {23114, FIVEexp*5 + Pexp + Qexp * 5 + Rexp * 2 + Sexp * 2},
  {256744, FIVEexp*3 + Pexp * 5 + Qexp + Rexp * 3 + Sexp * 2},
  {-19608, FIVEexp*5 + Pexp * 2 + Qexp * 3 + Rexp * 3 + Sexp * 2},
  {-1376, FIVEexp*7 + Pexp * 3 + Qexp + Rexp * 4 + Sexp * 2},
  {-5024, FIVEexp*6 + Qexp * 3 + Rexp * 4 + Sexp * 2},
  {1408, FIVEexp*6 + Pexp + Qexp + Rexp * 5 + Sexp * 2},
  {-2916, FIVEexp*2 + Pexp * 10 + Sexp * 3},
  {-1944, FIVEexp*4 + Pexp * 7 + Qexp * 2 + Sexp * 3},
  {-188, FIVEexp*6 + Pexp * 4 + Qexp * 4 + Sexp * 3},
  {586, FIVEexp*6 + Pexp + Qexp * 6 + Sexp * 3},
  {20898, FIVEexp*3 + Pexp * 8 + Rexp + Sexp * 3},
  {26496, FIVEexp*4 + Pexp * 5 + Qexp * 2 + Rexp + Sexp * 3},
  {-966, FIVEexp*7 + Pexp * 2 + Qexp * 4 + Rexp + Sexp * 3},
  {-10476, FIVEexp*5 + Pexp * 6 + Rexp * 2 + Sexp * 3},
  {2164, FIVEexp*7 + Pexp * 3 + Qexp * 2 + Rexp * 2 + Sexp * 3},
  {1558, FIVEexp*7 + Qexp * 4 + Rexp * 2 + Sexp * 3},
  {10256, FIVEexp*6 + Pexp * 4 + Rexp * 3 + Sexp * 3},
  {14064, FIVEexp*6 + Pexp + Qexp * 2 + Rexp * 3 + Sexp * 3},
  {-20288, FIVEexp*6 + Pexp * 2 + Rexp * 4 + Sexp * 3},
  {3328, FIVEexp*7 + Rexp * 5 + Sexp * 3},
  {162, FIVEexp*6 + Pexp * 6 + Qexp + Sexp * 4},
  {288, FIVEexp*7 + Pexp * 3 + Qexp * 3 + Sexp * 4},
  {102, FIVEexp*8 + Qexp * 5 + Sexp * 4},
  {-17046, FIVEexp*6 + Pexp * 4 + Qexp + Rexp + Sexp * 4},
  {-9938, FIVEexp*7 + Pexp + Qexp * 3 + Rexp + Sexp * 4},
  {404, FIVEexp*9 + Pexp * 2 + Qexp + Rexp * 2 + Sexp * 4},
  {-3504, FIVEexp*8 + Qexp + Rexp * 3 + Sexp * 4},
  {864, FIVEexp*7 + Pexp * 5 + Sexp * 5},
  {226, FIVEexp*9 + Pexp * 2 + Qexp * 2 + Sexp * 5},
  {-798, FIVEexp*8 + Pexp * 3 + Rexp + Sexp * 5},
  {914, FIVEexp*9 + Qexp * 2 + Rexp + Sexp * 5},
  {56, FIVEexp*10 + Pexp + Rexp * 2 + Sexp * 5},
  {-26, FIVEexp*11 + Pexp + Qexp + Sexp * 6},
  {2, FIVEexp*12 + Sexp * 7},
  {0, END_COEFF},

  // c[3]
  {24, FIVEexp*1 + Pexp * 5 + Qexp * 7},
  {162, FIVEexp*1 + Pexp * 2 + Qexp * 9},
  {-256, FIVEexp*1 + Pexp * 6 + Qexp * 5 + Rexp},
  {-1832, FIVEexp*1 + Pexp * 3 + Qexp * 7 + Rexp},
  {756, FIVEexp*1 + Qexp * 9 + Rexp},
  {906, FIVEexp*1 + Pexp * 7 + Qexp * 3 + Rexp * 2},
  {7328, FIVEexp*1 + Pexp * 4 + Qexp * 5 + Rexp * 2},
  {-9054, FIVEexp*1 + Pexp + Qexp * 7 + Rexp * 2},
  {-216, FIVEexp*2 + Pexp * 8 + Qexp + Rexp * 3},
  {-12184, FIVEexp*1 + Pexp * 5 + Qexp * 3 + Rexp * 3},
  {8002, FIVEexp*2 + Pexp * 2 + Qexp * 5 + Rexp * 3},
  {1248, FIVEexp*2 + Pexp * 6 + Qexp + Rexp * 4},
  {-3808, FIVEexp*3 + Pexp * 3 + Qexp * 3 + Rexp * 4},
  {-15128, FIVEexp*2 + Qexp * 5 + Rexp * 4},
  {20864, FIVEexp*2 + Pexp * 4 + Qexp + Rexp * 5},
  {14976, FIVEexp*3 + Pexp + Qexp * 3 + Rexp * 5},
  {-3584, FIVEexp*4 + Pexp * 2 + Qexp + Rexp * 6},
  {288, FIVEexp*1 + Pexp * 7 + Qexp * 4 + Sexp},
  {3062, FIVEexp*1 + Pexp * 4 + Qexp * 6 + Sexp},
  {2376, FIVEexp*2 + Pexp + Qexp * 8 + Sexp},
  {-1836, FIVEexp*1 + Pexp * 8 + Qexp * 2 + Rexp + Sexp},
  {-23048, FIVEexp*1 + Pexp * 5 + Qexp * 4 + Rexp + Sexp},
  {-23586, FIVEexp*2 + Pexp * 2 + Qexp * 6 + Rexp + Sexp},
  {648, FIVEexp*2 + Pexp * 9 + Rexp * 2 + Sexp},
  {63342, FIVEexp*1 + Pexp * 6 + Qexp * 2 + Rexp * 2 + Sexp},
  {20382, FIVEexp*3 + Pexp * 3 + Qexp * 4 + Rexp * 2 + Sexp},
  {17424, FIVEexp*3 + Qexp * 6 + Rexp * 2 + Sexp},
  {-10368, FIVEexp*2 + Pexp * 7 + Rexp * 3 + Sexp},
  {-32984, FIVEexp*3 + Pexp * 4 + Qexp * 2 + Rexp * 3 + Sexp},
  {-3744, FIVEexp*5 + Pexp + Qexp * 4 + Rexp * 3 + Sexp},
  {37504, FIVEexp*2 + Pexp * 5 + Rexp * 4 + Sexp},
  {26144, FIVEexp*4 + Pexp * 2 + Qexp * 2 + Rexp * 4 + Sexp},
  {-1024, FIVEexp*4 + Pexp * 3 + Rexp * 5 + Sexp},
  {896, FIVEexp*5 + Qexp * 2 + Rexp * 5 + Sexp},
  {-486, FIVEexp*1 + Pexp * 9 + Qexp + Sexp * 2},
  {-2178, FIVEexp*2 + Pexp * 6 + Qexp * 3 + Sexp * 2},
  {-2284, FIVEexp*3 + Pexp * 3 + Qexp * 5 + Sexp * 2},
  {-4428, FIVEexp*4 + Qexp * 7 + Sexp * 2},
  {1728, FIVEexp*2 + Pexp * 7 + Qexp + Rexp + Sexp * 2},
  {-1466, FIVEexp*4 + Pexp * 4 + Qexp * 3 + Rexp + Sexp * 2},
  {23172, FIVEexp*4 + Pexp + Qexp * 5 + Rexp + Sexp * 2},
  {38448, FIVEexp*3 + Pexp * 5 + Qexp + Rexp * 2 + Sexp * 2},
  {-4228, FIVEexp*5 + Pexp * 2 + Qexp * 3 + Rexp * 2 + Sexp * 2},
  {-8128, FIVEexp*5 + Pexp * 3 + Qexp + Rexp * 3 + Sexp * 2},
  {-48, FIVEexp*8 + Qexp * 3 + Rexp * 3 + Sexp * 2},
  {512, FIVEexp*6 + Pexp + Qexp + Rexp * 4 + Sexp * 2},
  {972, FIVEexp*3 + Pexp * 8 + Sexp * 3},
  {3294, FIVEexp*4 + Pexp * 5 + Qexp * 2 + Sexp * 3},
  {-426, FIVEexp*6 + Pexp * 2 + Qexp * 4 + Sexp * 3},
  {-10746, FIVEexp*4 + Pexp * 6 + Rexp + Sexp * 3},
  {1544, FIVEexp*6 + Pexp * 3 + Qexp * 2 + Rexp + Sexp * 3},
  {1528, FIVEexp*6 + Qexp * 4 + Rexp + Sexp * 3},
  {552, FIVEexp*7 + Pexp * 4 + Rexp * 2 + Sexp * 3},
  {2928, FIVEexp*6 + Pexp + Qexp * 2 + Rexp * 2 + Sexp * 3},
  {-224, FIVEexp*8 + Pexp * 2 + Rexp * 3 + Sexp * 3},
  {896, FIVEexp*7 + Rexp * 4 + Sexp * 3},
  {-2844, FIVEexp*6 + Pexp * 4 + Qexp + Sexp * 4},
  {-1382, FIVEexp*7 + Pexp + Qexp * 3 + Sexp * 4},
  {2042, FIVEexp*7 + Pexp * 2 + Qexp + Rexp + Sexp * 4},
  {-728, FIVEexp*8 + Qexp + Rexp * 2 + Sexp * 4},
  {18, FIVEexp*8 + Pexp * 3 + Sexp * 5},
  {136, FIVEexp*9 + Qexp * 2 + Sexp * 5},
  {16, FIVEexp*9 + Pexp + Rexp + Sexp * 5},
  {0, END_COEFF},

  // c[4]
  {32, FIVEexp*1 + Pexp * 3 + Qexp * 7},
  {216, FIVEexp*1 + Qexp * 9},
  {-216, FIVEexp*1 + Pexp * 4 + Qexp * 5 + Rexp},
  {-1746, FIVEexp*1 + Pexp + Qexp * 7 + Rexp},
  {302, FIVEexp*1 + Pexp * 5 + Qexp * 3 + Rexp * 2},
  {4084, FIVEexp*1 + Pexp * 2 + Qexp * 5 + Rexp * 2},
  {144, FIVEexp*1 + Pexp * 6 + Qexp + Rexp * 3},
  {-928, FIVEexp*2 + Pexp * 3 + Qexp * 3 + Rexp * 3},
  {-3196, FIVEexp*2 + Qexp * 5 + Rexp * 3},
  {1408, FIVEexp*2 + Pexp * 4 + Qexp + Rexp * 4},
  {3232, FIVEexp*3 + Pexp + Qexp * 3 + Rexp * 4},
  {-768, FIVEexp*4 + Pexp * 2 + Qexp + Rexp * 5},
  {192, FIVEexp*1 + Pexp * 5 + Qexp * 4 + Sexp},
  {114, FIVEexp*2 + Pexp * 2 + Qexp * 6 + Sexp},
  {108, FIVEexp*1 + Pexp * 6 + Qexp * 2 + Rexp + Sexp},
  {508, FIVEexp*3 + Pexp * 3 + Qexp * 4 + Rexp + Sexp},
  {2556, FIVEexp*3 + Qexp * 6 + Rexp + Sexp},
  {-1512, FIVEexp*1 + Pexp * 7 + Rexp * 2 + Sexp},
  {-2028, FIVEexp*3 + Pexp * 4 + Qexp * 2 + Rexp * 2 + Sexp},
  {-578, FIVEexp*5 + Pexp + Qexp * 4 + Rexp * 2 + Sexp},
  {3648, FIVEexp*2 + Pexp * 5 + Rexp * 3 + Sexp},
  {832, FIVEexp*5 + Pexp * 2 + Qexp * 2 + Rexp * 3 + Sexp},
  {-128, FIVEexp*4 + Pexp * 3 + Rexp * 4 + Sexp},
  {192, FIVEexp*5 + Qexp * 2 + Rexp * 4 + Sexp},
  {-162, FIVEexp*2 + Pexp * 7 + Qexp + Sexp * 2},
  {-192, FIVEexp*4 + Pexp * 4 + Qexp * 3 + Sexp * 2},
  {-438, FIVEexp*4 + Pexp + Qexp * 5 + Sexp * 2},
  {3402, FIVEexp*3 + Pexp * 5 + Qexp + Rexp + Sexp * 2},
  {744, FIVEexp*5 + Pexp * 2 + Qexp * 3 + Rexp + Sexp * 2},
  {-1728, FIVEexp*5 + Pexp * 3 + Qexp + Rexp * 2 + Sexp * 2},
  {-184, FIVEexp*6 + Qexp * 3 + Rexp * 2 + Sexp * 2},
  {96, FIVEexp*6 + Pexp + Qexp + Rexp * 3 + Sexp * 2},
  {-486, FIVEexp*4 + Pexp * 6 + Sexp * 3},
  {-54, FIVEexp*6 + Pexp * 3 + Qexp * 2 + Sexp * 3},
  {-52, FIVEexp*6 + Qexp * 4 + Sexp * 3},
  {324, FIVEexp*6 + Pexp * 4 + Rexp + Sexp * 3},
  {852, FIVEexp*6 + Pexp + Qexp * 2 + Rexp + Sexp * 3},
  {-928, FIVEexp*6 + Pexp * 2 + Rexp * 2 + Sexp * 3},
  {192, FIVEexp*7 + Rexp * 3 + Sexp * 3},
  {-48, FIVEexp*7 + Pexp * 2 + Qexp + Sexp * 4},
  {-92, FIVEexp*8 + Qexp + Rexp + Sexp * 4},
  {6, FIVEexp*9 + Pexp + Sexp * 5},
  {0, END_COEFF},

  // c[5]
  {16, FIVEexp*1 + Pexp * 4 + Qexp * 5},
  {108, FIVEexp*1 + Pexp + Qexp * 7},
  {-24, FIVEexp*2 + Pexp * 5 + Qexp * 3 + Rexp},
  {-954, FIVEexp*1 + Pexp * 2 + Qexp * 5 + Rexp},
  {216, FIVEexp*1 + Pexp * 6 + Qexp + Rexp * 2},
  {448, FIVEexp*2 + Pexp * 3 + Qexp * 3 + Rexp * 2},
  {-486, FIVEexp*2 + Qexp * 5 + Rexp * 2},
  {-192, FIVEexp*2 + Pexp * 4 + Qexp + Rexp * 3},
  {512, FIVEexp*3 + Pexp + Qexp * 3 + Rexp * 3},
  {-128, FIVEexp*4 + Pexp * 2 + Qexp + Rexp * 4},
  {216, FIVEexp*1 + Pexp * 6 + Qexp * 2 + Sexp},
  {106, FIVEexp*3 + Pexp * 3 + Qexp * 4 + Sexp},
  {432, FIVEexp*3 + Qexp * 6 + Sexp},
  {-648, FIVEexp*1 + Pexp * 7 + Rexp + Sexp},
  {-18, FIVEexp*5 + Pexp * 4 + Qexp * 2 + Rexp + Sexp},
  {-108, FIVEexp*5 + Pexp + Qexp * 4 + Rexp + Sexp},
  {1728, FIVEexp*2 + Pexp * 5 + Rexp * 2 + Sexp},
  {896, FIVEexp*4 + Pexp * 2 + Qexp * 2 + Rexp * 2 + Sexp},
  {-128, FIVEexp*4 + Pexp * 3 + Rexp * 3 + Sexp},
  {32, FIVEexp*5 + Qexp * 2 + Rexp * 3 + Sexp},
  {54, FIVEexp*3 + Pexp * 5 + Qexp + Sexp * 2},
  {72, FIVEexp*5 + Pexp * 2 + Qexp * 3 + Sexp * 2},
  {-288, FIVEexp*5 + Pexp * 3 + Qexp + Rexp + Sexp * 2},
  {-36, FIVEexp*6 + Qexp * 3 + Rexp + Sexp * 2},
  {32, FIVEexp*6 + Pexp + Qexp + Rexp * 2 + Sexp * 2},
  {54, FIVEexp*6 + Pexp * 4 + Sexp * 3},
  {124, FIVEexp*6 + Pexp + Qexp * 2 + Sexp * 3},
  {-192, FIVEexp*6 + Pexp * 2 + Rexp + Sexp * 3},
  {32, FIVEexp*7 + Rexp * 2 + Sexp * 3},
  {-14, FIVEexp*8 + Qexp + Sexp * 4},
  {0, END_COEFF},

  // o[0]
  {-64, FIVEexp*2 + Pexp * 10 + Qexp * 10},
  {-944, FIVEexp*2 + Pexp * 7 + Qexp * 12},
  {-3456, FIVEexp*2 + Pexp * 4 + Qexp * 14},
  {992, FIVEexp*2 + Pexp * 11 + Qexp * 8 + Rexp},
  {16768, FIVEexp*2 + Pexp * 8 + Qexp * 10 + Rexp},
  {74018, FIVEexp*2 + Pexp * 5 + Qexp * 12 + Rexp},
  {35856, FIVEexp*2 + Pexp * 2 + Qexp * 14 + Rexp},
  {-5552, FIVEexp*2 + Pexp * 12 + Qexp * 6 + Rexp * 2},
  {-116876, FIVEexp*2 + Pexp * 9 + Qexp * 8 + Rexp * 2},
  {-691808, FIVEexp*2 + Pexp * 6 + Qexp * 10 + Rexp * 2},
  {-217022, FIVEexp*3 + Pexp * 3 + Qexp * 12 + Rexp * 2},
  {-1043064, FIVEexp*2 + Qexp * 14 + Rexp * 2},
  {13032, FIVEexp*2 + Pexp * 13 + Qexp * 4 + Rexp * 3},
  {399754, FIVEexp*2 + Pexp * 10 + Qexp * 6 + Rexp * 3},
  {704084, FIVEexp*3 + Pexp * 7 + Qexp * 8 + Rexp * 3},
  {10961906, FIVEexp*2 + Pexp * 4 + Qexp * 10 + Rexp * 3},
  {16406856, FIVEexp*2 + Pexp + Qexp * 12 + Rexp * 3},
  {-10368, FIVEexp*2 + Pexp * 14 + Qexp * 2 + Rexp * 4},
  {-685884, FIVEexp*2 + Pexp * 11 + Qexp * 4 + Rexp * 4},
  {-10171566, FIVEexp*2 + Pexp * 8 + Qexp * 6 + Rexp * 4},
  {-52741929, FIVEexp*2 + Pexp * 5 + Qexp * 8 + Rexp * 4},
  {-105343939, FIVEexp*2 + Pexp * 2 + Qexp * 10 + Rexp * 4},
  {101088, FIVEexp*3 + Pexp * 12 + Qexp * 2 + Rexp * 5},
  {3111288, FIVEexp*3 + Pexp * 9 + Qexp * 4 + Rexp * 5},
  {130788189, FIVEexp*2 + Pexp * 6 + Qexp * 6 + Rexp * 5},
  {351647692, FIVEexp*2 + Pexp * 3 + Qexp * 8 + Rexp * 5},
  {3742423, FIVEexp*2 + Qexp * 10 + Rexp * 5},
  {-9134464, FIVEexp*2 + Pexp * 10 + Qexp * 2 + Rexp * 6},
  {-158047968, FIVEexp*2 + Pexp * 7 + Qexp * 4 + Rexp * 6},
  {-651079244, FIVEexp*2 + Pexp * 4 + Qexp * 6 + Rexp * 6},
  {-12777816, FIVEexp*3 + Pexp + Qexp * 8 + Rexp * 6},
  {77915968, FIVEexp*2 + Pexp * 8 + Qexp * 2 + Rexp * 7},
  {136301184, FIVEexp*3 + Pexp * 5 + Qexp * 4 + Rexp * 7},
  {14271584, FIVEexp*4 + Pexp * 2 + Qexp * 6 + Rexp * 7},
  {-12275456, FIVEexp*4 + Pexp * 6 + Qexp * 2 + Rexp * 8},
  {-991744, FIVEexp*6 + Pexp * 3 + Qexp * 4 + Rexp * 8},
  {270336, FIVEexp*6 + Qexp * 6 + Rexp * 8},
  {637952, FIVEexp*6 + Pexp * 4 + Qexp * 2 + Rexp * 9},
  {-110592, FIVEexp*7 + Pexp + Qexp * 4 + Rexp * 9},
  {12288, FIVEexp*8 + Pexp * 2 + Qexp * 2 + Rexp * 10},
  {-2208, FIVEexp*2 + Pexp * 12 + Qexp * 7 + Sexp},
  {-27424, FIVEexp*2 + Pexp * 9 + Qexp * 9 + Sexp},
  {8226, FIVEexp*3 + Pexp * 6 + Qexp * 11 + Sexp},
  {12048, FIVEexp*5 + Pexp * 3 + Qexp * 13 + Sexp},
  {7128, FIVEexp*6 + Qexp * 15 + Sexp},
  {23328, FIVEexp*2 + Pexp * 13 + Qexp * 5 + Rexp + Sexp},
  {363024, FIVEexp*2 + Pexp * 10 + Qexp * 7 + Rexp + Sexp},
  {-35326, FIVEexp*2 + Pexp * 7 + Qexp * 9 + Rexp + Sexp},
  {-4054646, FIVEexp*3 + Pexp * 4 + Qexp * 11 + Rexp + Sexp},
  {-573804, FIVEexp*5 + Pexp + Qexp * 13 + Rexp + Sexp},
  {-74088, FIVEexp*2 + Pexp * 14 + Qexp * 3 + Rexp * 2 + Sexp},
  {-331482, FIVEexp*3 + Pexp * 11 + Qexp * 5 + Rexp * 2 + Sexp},
  {-3222668, FIVEexp*2 + Pexp * 8 + Qexp * 7 + Rexp * 2 + Sexp},
  {99426944, FIVEexp*2 + Pexp * 5 + Qexp * 9 + Rexp * 2 + Sexp},
  {91538289, FIVEexp*3 + Pexp * 2 + Qexp * 11 + Rexp * 2 + Sexp},
  {62208, FIVEexp*2 + Pexp * 15 + Qexp + Rexp * 3 + Sexp},
  {3233844, FIVEexp*2 + Pexp * 12 + Qexp * 3 + Rexp * 3 + Sexp},
  {22596272, FIVEexp*2 + Pexp * 9 + Qexp * 5 + Rexp * 3 + Sexp},
  {-179720496, FIVEexp*2 + Pexp * 6 + Qexp * 7 + Rexp * 3 + Sexp},
  {-56787826, FIVEexp*4 + Pexp * 3 + Qexp * 9 + Rexp * 3 + Sexp},
  {-1263891, FIVEexp*4 + Qexp * 11 + Rexp * 3 + Sexp},
  {-574128, FIVEexp*3 + Pexp * 13 + Qexp + Rexp * 4 + Sexp},
  {-62045968, FIVEexp*2 + Pexp * 10 + Qexp * 3 + Rexp * 4 + Sexp},
  {-70937516, FIVEexp*2 + Pexp * 7 + Qexp * 5 + Rexp * 4 + Sexp},
  {83132149, FIVEexp*4 + Pexp * 4 + Qexp * 7 + Rexp * 4 + Sexp},
  {191549, FIVEexp*7 + Pexp + Qexp * 9 + Rexp * 4 + Sexp},
  {49262784, FIVEexp*2 + Pexp * 11 + Qexp + Rexp * 5 + Sexp},
  {481719104, FIVEexp*2 + Pexp * 8 + Qexp * 3 + Rexp * 5 + Sexp},
  {-1086040448, FIVEexp*2 + Pexp * 5 + Qexp * 5 + Rexp * 5 + Sexp},
  {-140928976, FIVEexp*4 + Pexp * 2 + Qexp * 7 + Rexp * 5 + Sexp},
  {-396519808, FIVEexp*2 + Pexp * 9 + Qexp + Rexp * 6 + Sexp},
  {-235584832, FIVEexp*3 + Pexp * 6 + Qexp * 3 + Rexp * 6 + Sexp},
  {231927744, FIVEexp*4 + Pexp * 3 + Qexp * 5 + Rexp * 6 + Sexp},
  {-14272416, FIVEexp*5 + Qexp * 7 + Rexp * 6 + Sexp},
  {56726016, FIVEexp*4 + Pexp * 7 + Qexp + Rexp * 7 + Sexp},
  {-4048896, FIVEexp*6 + Pexp * 4 + Qexp * 3 + Rexp * 7 + Sexp},
  {3874816, FIVEexp*6 + Pexp + Qexp * 5 + Rexp * 7 + Sexp},
  {-1923072, FIVEexp*6 + Pexp * 5 + Qexp + Rexp * 8 + Sexp},
  {474112, FIVEexp*7 + Pexp * 2 + Qexp * 3 + Rexp * 8 + Sexp},
  {-155648, FIVEexp*8 + Pexp * 3 + Qexp + Rexp * 9 + Sexp},
  {-24576, FIVEexp*9 + Qexp * 3 + Rexp * 9 + Sexp},
  {-24624, FIVEexp*2 + Pexp * 14 + Qexp * 4 + Sexp * 2},
  {-84196, FIVEexp*3 + Pexp * 11 + Qexp * 6 + Sexp * 2},
  {-10826, FIVEexp*5 + Pexp * 8 + Qexp * 8 + Sexp * 2},
  {71298, FIVEexp*5 + Pexp * 5 + Qexp * 10 + Sexp * 2},
  {14076, FIVEexp*7 + Pexp * 2 + Qexp * 12 + Sexp * 2},
  {134136, FIVEexp*2 + Pexp * 15 + Qexp * 2 + Rexp + Sexp * 2},
  {2970774, FIVEexp*2 + Pexp * 12 + Qexp * 4 + Rexp + Sexp * 2},
  {2211566, FIVEexp*3 + Pexp * 9 + Qexp * 6 + Rexp + Sexp * 2},
  {-837952, FIVEexp*5 + Pexp * 6 + Qexp * 8 + Rexp + Sexp * 2},
  {-4631438, FIVEexp*5 + Pexp * 3 + Qexp * 10 + Rexp + Sexp * 2},
  {17712, FIVEexp*7 + Qexp * 12 + Rexp + Sexp * 2},
  {-93312, FIVEexp*2 + Pexp * 16 + Rexp * 2 + Sexp * 2},
  {-5310036, FIVEexp*2 + Pexp * 13 + Qexp * 2 + Rexp * 2 + Sexp * 2},
  {-36031006, FIVEexp*2 + Pexp * 10 + Qexp * 4 + Rexp * 2 + Sexp * 2},
  {65993956, FIVEexp*3 + Pexp * 7 + Qexp * 6 + Rexp * 2 + Sexp * 2},
  {3809651, FIVEexp*6 + Pexp * 4 + Qexp * 8 + Rexp * 2 + Sexp * 2},
  {-12893718, FIVEexp*5 + Pexp + Qexp * 10 + Rexp * 2 + Sexp * 2},
  {1026432, FIVEexp*3 + Pexp * 14 + Rexp * 3 + Sexp * 2},
  {126409284, FIVEexp*2 + Pexp * 11 + Qexp * 2 + Rexp * 3 + Sexp * 2},
  {13327328, FIVEexp*4 + Pexp * 8 + Qexp * 4 + Rexp * 3 + Sexp * 2},
  {-72893534, FIVEexp*4 + Pexp * 5 + Qexp * 6 + Rexp * 3 + Sexp * 2},
  {3804837, FIVEexp*7 + Pexp * 2 + Qexp * 8 + Rexp * 3 + Sexp * 2},
  {-110794176, FIVEexp*2 + Pexp * 12 + Rexp * 4 + Sexp * 2},
  {-59305552, FIVEexp*4 + Pexp * 9 + Qexp * 2 + Rexp * 4 + Sexp * 2},
  {-5812003, FIVEexp*6 + Pexp * 6 + Qexp * 4 + Rexp * 4 + Sexp * 2},
  {-40185728, FIVEexp*6 + Pexp * 3 + Qexp * 6 + Rexp * 4 + Sexp * 2},
  {2321247, FIVEexp*7 + Qexp * 8 + Rexp * 4 + Sexp * 2},
  {1237877312, FIVEexp*2 + Pexp * 10 + Rexp * 5 + Sexp * 2},
  {399927488, FIVEexp*4 + Pexp * 7 + Qexp * 2 + Rexp * 5 + Sexp * 2},
  {51389108, FIVEexp*6 + Pexp * 4 + Qexp * 4 + Rexp * 5 + Sexp * 2},
  {-207104, FIVEexp*8 + Pexp + Qexp * 6 + Rexp * 5 + Sexp * 2},
  {-307419712, FIVEexp*4 + Pexp * 8 + Rexp * 6 + Sexp * 2},
  {-298445312, FIVEexp*5 + Pexp * 5 + Qexp * 2 + Rexp * 6 + Sexp * 2},
  {-12073504, FIVEexp*7 + Pexp * 2 + Qexp * 4 + Rexp * 6 + Sexp * 2},
  {42138368, FIVEexp*6 + Pexp * 6 + Rexp * 7 + Sexp * 2},
  {24713216, FIVEexp*7 + Pexp * 3 + Qexp * 2 + Rexp * 7 + Sexp * 2},
  {1520128, FIVEexp*8 + Qexp * 4 + Rexp * 7 + Sexp * 2},
  {-2976768, FIVEexp*8 + Pexp * 4 + Rexp * 8 + Sexp * 2},
  {-28672, FIVEexp*10 + Pexp + Qexp * 2 + Rexp * 8 + Sexp * 2},
  {86016, FIVEexp*10 + Pexp * 2 + Rexp * 9 + Sexp * 2},
  {-17496, FIVEexp*3 + Pexp * 16 + Qexp + Sexp * 3},
  {-15174, FIVEexp*5 + Pexp * 13 + Qexp * 3 + Sexp * 3},
  {-57798, FIVEexp*5 + Pexp * 10 + Qexp * 5 + Sexp * 3},
  {5712, FIVEexp*8 + Pexp * 7 + Qexp * 7 + Sexp * 3},
  {228579, FIVEexp*7 + Pexp * 4 + Qexp * 9 + Sexp * 3},
  {612, FIVEexp*11 + Pexp + Qexp * 11 + Sexp * 3},
  {197316, FIVEexp*3 + Pexp * 14 + Qexp + Rexp + Sexp * 3},
  {-273078, FIVEexp*5 + Pexp * 11 + Qexp * 3 + Rexp + Sexp * 3},
  {-8300702, FIVEexp*5 + Pexp * 8 + Qexp * 5 + Rexp + Sexp * 3},
  {-2272532, FIVEexp*7 + Pexp * 5 + Qexp * 7 + Rexp + Sexp * 3},
  {-3814129, FIVEexp*7 + Pexp * 2 + Qexp * 9 + Rexp + Sexp * 3},
  {-1227744, FIVEexp*3 + Pexp * 12 + Qexp + Rexp * 2 + Sexp * 3},
  {386408, FIVEexp*7 + Pexp * 9 + Qexp * 3 + Rexp * 2 + Sexp * 3},
  {110095943, FIVEexp*5 + Pexp * 6 + Qexp * 5 + Rexp * 2 + Sexp * 3},
  {6835544, FIVEexp*7 + Pexp * 3 + Qexp * 7 + Rexp * 2 + Sexp * 3},
  {-4361299, FIVEexp*7 + Qexp * 9 + Rexp * 2 + Sexp * 3},
  {-14489904, FIVEexp*4 + Pexp * 10 + Qexp + Rexp * 3 + Sexp * 3},
  {-4125628, FIVEexp*7 + Pexp * 7 + Qexp * 3 + Rexp * 3 + Sexp * 3},
  {-9851291, FIVEexp*7 + Pexp * 4 + Qexp * 5 + Rexp * 3 + Sexp * 3},
  {-1065944, FIVEexp*7 + Pexp + Qexp * 7 + Rexp * 3 + Sexp * 3},
  {10499904, FIVEexp*6 + Pexp * 8 + Qexp + Rexp * 4 + Sexp * 3},
  {17681392, FIVEexp*7 + Pexp * 5 + Qexp * 3 + Rexp * 4 + Sexp * 3},
  {39526656, FIVEexp*7 + Pexp * 2 + Qexp * 5 + Rexp * 4 + Sexp * 3},
  {-16225984, FIVEexp*7 + Pexp * 6 + Qexp + Rexp * 5 + Sexp * 3},
  {-97825664, FIVEexp*7 + Pexp * 3 + Qexp * 3 + Rexp * 5 + Sexp * 3},
  {-35326432, FIVEexp*7 + Qexp * 5 + Rexp * 5 + Sexp * 3},
  {58851328, FIVEexp*7 + Pexp * 4 + Qexp + Rexp * 6 + Sexp * 3},
  {4726272, FIVEexp*8 + Pexp + Qexp * 3 + Rexp * 6 + Sexp * 3},
  {-3587072, FIVEexp*9 + Pexp * 2 + Qexp + Rexp * 7 + Sexp * 3},
  {-24576, FIVEexp*11 + Qexp + Rexp * 8 + Sexp * 3},
  {5832, FIVEexp*5 + Pexp * 15 + Sexp * 4},
  {3402, FIVEexp*8 + Pexp * 12 + Qexp * 2 + Sexp * 4},
  {316533, FIVEexp*7 + Pexp * 9 + Qexp * 4 + Sexp * 4},
  {433836, FIVEexp*8 + Pexp * 6 + Qexp * 6 + Sexp * 4},
  {8464, FIVEexp*11 + Pexp * 3 + Qexp * 8 + Sexp * 4},
  {4584, FIVEexp*11 + Qexp * 10 + Sexp * 4},
  {227448, FIVEexp*5 + Pexp * 13 + Rexp + Sexp * 4},
  {-238221, FIVEexp*7 + Pexp * 10 + Qexp * 2 + Rexp + Sexp * 4},
  {-4023608, FIVEexp*7 + Pexp * 7 + Qexp * 4 + Rexp + Sexp * 4},
  {-2120806, FIVEexp*8 + Pexp * 4 + Qexp * 6 + Rexp + Sexp * 4},
  {1888, FIVEexp*12 + Pexp + Qexp * 8 + Rexp + Sexp * 4},
  {-8215992, FIVEexp*5 + Pexp * 11 + Rexp * 2 + Sexp * 4},
  {-668461, FIVEexp*7 + Pexp * 8 + Qexp * 2 + Rexp * 2 + Sexp * 4},
  {-4955648, FIVEexp*7 + Pexp * 5 + Qexp * 4 + Rexp * 2 + Sexp * 4},
  {-12118221, FIVEexp*8 + Pexp * 2 + Qexp * 6 + Rexp * 2 + Sexp * 4},
  {4402072, FIVEexp*7 + Pexp * 9 + Rexp * 3 + Sexp * 4},
  {42606827, FIVEexp*7 + Pexp * 6 + Qexp * 2 + Rexp * 3 + Sexp * 4},
  {8464592, FIVEexp*9 + Pexp * 3 + Qexp * 4 + Rexp * 3 + Sexp * 4},
  {3062074, FIVEexp*9 + Qexp * 6 + Rexp * 3 + Sexp * 4},
  {-29382368, FIVEexp*7 + Pexp * 7 + Rexp * 4 + Sexp * 4},
  {-7582628, FIVEexp*9 + Pexp * 4 + Qexp * 2 + Rexp * 4 + Sexp * 4},
  {-132944, FIVEexp*11 + Pexp + Qexp * 4 + Rexp * 4 + Sexp * 4},
  {99436416, FIVEexp*7 + Pexp * 5 + Rexp * 5 + Sexp * 4},
  {17506592, FIVEexp*9 + Pexp * 2 + Qexp * 2 + Rexp * 5 + Sexp * 4},
  {-6071808, FIVEexp*9 + Pexp * 3 + Rexp * 6 + Sexp * 4},
  {1078272, FIVEexp*10 + Qexp * 2 + Rexp * 6 + Sexp * 4},
  {24576, FIVEexp*12 + Pexp + Rexp * 7 + Sexp * 4},
  {29322, FIVEexp*8 + Pexp * 11 + Qexp + Sexp * 5},
  {15336, FIVEexp*10 + Pexp * 8 + Qexp * 3 + Sexp * 5},
  {55863, FIVEexp*10 + Pexp * 5 + Qexp * 5 + Sexp * 5},
  {432, FIVEexp*13 + Pexp * 2 + Qexp * 7 + Sexp * 5},
  {-950562, FIVEexp*8 + Pexp * 9 + Qexp + Rexp + Sexp * 5},
  {-354453, FIVEexp*10 + Pexp * 6 + Qexp * 3 + Rexp + Sexp * 5},
  {-811098, FIVEexp*10 + Pexp * 3 + Qexp * 5 + Rexp + Sexp * 5},
  {-3889, FIVEexp*13 + Qexp * 7 + Rexp + Sexp * 5},
  {6189608, FIVEexp*8 + Pexp * 7 + Qexp + Rexp * 2 + Sexp * 5},
  {22387, FIVEexp*12 + Pexp * 4 + Qexp * 3 + Rexp * 2 + Sexp * 5},
  {607962, FIVEexp*10 + Pexp + Qexp * 5 + Rexp * 2 + Sexp * 5},
  {-5458512, FIVEexp*9 + Pexp * 5 + Qexp + Rexp * 3 + Sexp * 5},
  {-10416, FIVEexp*14 + Pexp * 2 + Qexp * 3 + Rexp * 3 + Sexp * 5},
  {492736, FIVEexp*11 + Pexp * 3 + Qexp + Rexp * 4 + Sexp * 5},
  {-135264, FIVEexp*12 + Qexp * 3 + Rexp * 4 + Sexp * 5},
  {-35328, FIVEexp*13 + Pexp + Qexp + Rexp * 5 + Sexp * 5},
  {9666, FIVEexp*10 + Pexp * 10 + Sexp * 6},
  {4494, FIVEexp*12 + Pexp * 7 + Qexp * 2 + Sexp * 6},
  {11591, FIVEexp*12 + Pexp * 4 + Qexp * 4 + Sexp * 6},
  {-2, FIVEexp*17 + Pexp + Qexp * 6 + Sexp * 6},
  {279, FIVEexp*10 + Pexp * 8 + Rexp + Sexp * 6},
  {57338, FIVEexp*12 + Pexp * 5 + Qexp * 2 + Rexp + Sexp * 6},
  {233859, FIVEexp*12 + Pexp * 2 + Qexp * 4 + Rexp + Sexp * 6},
  {-871681, FIVEexp*10 + Pexp * 6 + Rexp * 2 + Sexp * 6},
  {-155408, FIVEexp*12 + Pexp * 3 + Qexp * 2 + Rexp * 2 + Sexp * 6},
  {135954, FIVEexp*12 + Qexp * 4 + Rexp * 2 + Sexp * 6},
  {206988, FIVEexp*12 + Pexp * 4 + Rexp * 3 + Sexp * 6},
  {641664, FIVEexp*12 + Pexp + Qexp * 2 + Rexp * 3 + Sexp * 6},
  {-348256, FIVEexp*12 + Pexp * 2 + Rexp * 4 + Sexp * 6},
  {41472, FIVEexp*12 + Rexp * 5 + Sexp * 6},
  {-2226, FIVEexp*13 + Pexp * 6 + Qexp + Sexp * 7},
  {-372, FIVEexp*15 + Pexp * 3 + Qexp * 3 + Sexp * 7},
  {-17, FIVEexp*16 + Qexp * 5 + Sexp * 7},
  {-38629, FIVEexp*13 + Pexp * 4 + Qexp + Rexp + Sexp * 7},
  {-5264, FIVEexp*15 + Pexp + Qexp * 3 + Rexp + Sexp * 7},
  {116736, FIVEexp*13 + Pexp * 2 + Qexp + Rexp * 2 + Sexp * 7},
  {-4384, FIVEexp*14 + Qexp + Rexp * 3 + Sexp * 7},
  {-143, FIVEexp*15 + Pexp * 5 + Sexp * 8},
  {-124, FIVEexp*17 + Pexp * 2 + Qexp * 2 + Sexp * 8},
  {3748, FIVEexp*15 + Pexp * 3 + Rexp + Sexp * 8},
  {67, FIVEexp*17 + Qexp * 2 + Rexp + Sexp * 8},
  {-11352, FIVEexp*15 + Pexp + Rexp * 2 + Sexp * 8},
  {1, FIVEexp*21 + Pexp + Qexp + Sexp * 9},
  {-1, FIVEexp*21 + Sexp * 10},
  {0, END_COEFF},

  // o[1]
  {64, FIVEexp*2 + Pexp * 11 + Qexp * 8},
  {832, FIVEexp*2 + Pexp * 8 + Qexp * 10},
  {1804, FIVEexp*2 + Pexp * 5 + Qexp * 12},
  {-6048, FIVEexp*2 + Pexp * 2 + Qexp * 14},
  {-768, FIVEexp*2 + Pexp * 12 + Qexp * 6 + Rexp},
  {-11728, FIVEexp*2 + Pexp * 9 + Qexp * 8 + Rexp},
  {-31784, FIVEexp*2 + Pexp * 6 + Qexp * 10 + Rexp},
  {105387, FIVEexp*2 + Pexp * 3 + Qexp * 12 + Rexp},
  {105624, FIVEexp*2 + Qexp * 14 + Rexp},
  {3024, FIVEexp*2 + Pexp * 13 + Qexp * 4 + Rexp * 2},
  {61164, FIVEexp*2 + Pexp * 10 + Qexp * 6 + Rexp * 2},
  {249334, FIVEexp*2 + Pexp * 7 + Qexp * 8 + Rexp * 2},
  {-480534, FIVEexp*2 + Pexp * 4 + Qexp * 10 + Rexp * 2},
  {-1162782, FIVEexp*2 + Pexp + Qexp * 12 + Rexp * 2},
  {-3888, FIVEexp*2 + Pexp * 14 + Qexp * 2 + Rexp * 3},
  {-228, FIVEexp*6 + Pexp * 11 + Qexp * 4 + Rexp * 3},
  {-1079396, FIVEexp*2 + Pexp * 8 + Qexp * 6 + Rexp * 3},
  {-636013, FIVEexp*2 + Pexp * 5 + Qexp * 8 + Rexp * 3},
  {3050684, FIVEexp*2 + Pexp * 2 + Qexp * 10 + Rexp * 3},
  {130896, FIVEexp*2 + Pexp * 12 + Qexp * 2 + Rexp * 4},
  {2379474, FIVEexp*2 + Pexp * 9 + Qexp * 4 + Rexp * 4},
  {8850803, FIVEexp*2 + Pexp * 6 + Qexp * 6 + Rexp * 4},
  {592522, FIVEexp*3 + Pexp * 3 + Qexp * 8 + Rexp * 4},
  {-480903, FIVEexp*4 + Qexp * 10 + Rexp * 4},
  {-1822776, FIVEexp*2 + Pexp * 10 + Qexp * 2 + Rexp * 5},
  {-3509328, FIVEexp*3 + Pexp * 7 + Qexp * 4 + Rexp * 5},
  {-711714, FIVEexp*4 + Pexp * 4 + Qexp * 6 + Rexp * 5},
  {783442, FIVEexp*5 + Pexp + Qexp * 8 + Rexp * 5},
  {465024, FIVEexp*4 + Pexp * 8 + Qexp * 2 + Rexp * 6},
  {273872, FIVEexp*5 + Pexp * 5 + Qexp * 4 + Rexp * 6},
  {-73496, FIVEexp*7 + Pexp * 2 + Qexp * 6 + Rexp * 6},
  {-41216, FIVEexp*6 + Pexp * 6 + Qexp * 2 + Rexp * 7},
  {356736, FIVEexp*6 + Pexp * 3 + Qexp * 4 + Rexp * 7},
  {297152, FIVEexp*6 + Qexp * 6 + Rexp * 7},
  {-108544, FIVEexp*6 + Pexp * 4 + Qexp * 2 + Rexp * 8},
  {-162048, FIVEexp*7 + Pexp + Qexp * 4 + Rexp * 8},
  {18432, FIVEexp*8 + Pexp * 2 + Qexp * 2 + Rexp * 9},
  {1728, FIVEexp*2 + Pexp * 13 + Qexp * 5 + Sexp},
  {4576, FIVEexp*3 + Pexp * 10 + Qexp * 7 + Sexp},
  {-2392, FIVEexp*2 + Pexp * 7 + Qexp * 9 + Sexp},
  {-193397, FIVEexp*3 + Pexp * 4 + Qexp * 11 + Sexp},
  {-23868, FIVEexp*5 + Pexp + Qexp * 13 + Sexp},
  {-2592, FIVEexp*3 + Pexp * 14 + Qexp * 3 + Rexp + Sexp},
  {-221256, FIVEexp*2 + Pexp * 11 + Qexp * 5 + Rexp + Sexp},
  {-148484, FIVEexp*2 + Pexp * 8 + Qexp * 7 + Rexp + Sexp},
  {11720371, FIVEexp*2 + Pexp * 5 + Qexp * 9 + Rexp + Sexp},
  {8924391, FIVEexp*3 + Pexp * 2 + Qexp * 11 + Rexp + Sexp},
  {23328, FIVEexp*2 + Pexp * 15 + Qexp + Rexp * 2 + Sexp},
  {733752, FIVEexp*2 + Pexp * 12 + Qexp * 3 + Rexp * 2 + Sexp},
  {3116444, FIVEexp*2 + Pexp * 9 + Qexp * 5 + Rexp * 2 + Sexp},
  {-38299553, FIVEexp*2 + Pexp * 6 + Qexp * 7 + Rexp * 2 + Sexp},
  {-8719458, FIVEexp*4 + Pexp * 3 + Qexp * 9 + Rexp * 2 + Sexp},
  {1536192, FIVEexp*4 + Qexp * 11 + Rexp * 2 + Sexp},
  {-947376, FIVEexp*2 + Pexp * 13 + Qexp + Rexp * 3 + Sexp},
  {-14950476, FIVEexp*2 + Pexp * 10 + Qexp * 3 + Rexp * 3 + Sexp},
  {-1117799, FIVEexp*2 + Pexp * 7 + Qexp * 5 + Rexp * 3 + Sexp},
  {3320237, FIVEexp*5 + Pexp * 4 + Qexp * 7 + Rexp * 3 + Sexp},
  {-297542, FIVEexp*6 + Pexp + Qexp * 9 + Rexp * 3 + Sexp},
  {15832656, FIVEexp*2 + Pexp * 11 + Qexp + Rexp * 4 + Sexp},
  {23287744, FIVEexp*3 + Pexp * 8 + Qexp * 3 + Rexp * 4 + Sexp},
  {-2920372, FIVEexp*5 + Pexp * 5 + Qexp * 5 + Rexp * 4 + Sexp},
  {-3742758, FIVEexp*5 + Pexp * 2 + Qexp * 7 + Rexp * 4 + Sexp},
  {-4845824, FIVEexp*4 + Pexp * 9 + Qexp + Rexp * 5 + Sexp},
  {-1040496, FIVEexp*5 + Pexp * 6 + Qexp * 3 + Rexp * 5 + Sexp},
  {3066512, FIVEexp*6 + Pexp * 3 + Qexp * 5 + Rexp * 5 + Sexp},
  {-387272, FIVEexp*7 + Qexp * 7 + Rexp * 5 + Sexp},
  {595456, FIVEexp*6 + Pexp * 7 + Qexp + Rexp * 6 + Sexp},
  {-550016, FIVEexp*7 + Pexp * 4 + Qexp * 3 + Rexp * 6 + Sexp},
  {80576, FIVEexp*8 + Pexp + Qexp * 5 + Rexp * 6 + Sexp},
  {139264, FIVEexp*6 + Pexp * 5 + Qexp + Rexp * 7 + Sexp},
  {158976, FIVEexp*8 + Pexp * 2 + Qexp * 3 + Rexp * 7 + Sexp},
  {-110592, FIVEexp*8 + Pexp * 3 + Qexp + Rexp * 8 + Sexp},
  {-36864, FIVEexp*9 + Qexp * 3 + Rexp * 8 + Sexp},
  {11664, FIVEexp*2 + Pexp * 15 + Qexp * 2 + Sexp * 2},
  {108108, FIVEexp*2 + Pexp * 12 + Qexp * 4 + Sexp * 2},
  {-309538, FIVEexp*3 + Pexp * 9 + Qexp * 6 + Sexp * 2},
  {-172449, FIVEexp*5 + Pexp * 6 + Qexp * 8 + Sexp * 2},
  {-516196, FIVEexp*5 + Pexp * 3 + Qexp * 10 + Sexp * 2},
  {4104, FIVEexp*7 + Qexp * 12 + Sexp * 2},
  {-34992, FIVEexp*2 + Pexp * 16 + Rexp + Sexp * 2},
  {-566676, FIVEexp*2 + Pexp * 13 + Qexp * 2 + Rexp + Sexp * 2},
  {7731396, FIVEexp*2 + Pexp * 10 + Qexp * 4 + Rexp + Sexp * 2},
  {29508164, FIVEexp*3 + Pexp * 7 + Qexp * 6 + Rexp + Sexp * 2},
  {743257, FIVEexp*6 + Pexp * 4 + Qexp * 8 + Rexp + Sexp * 2},
  {-4995162, FIVEexp*5 + Pexp + Qexp * 10 + Rexp + Sexp * 2},
  {1761264, FIVEexp*2 + Pexp * 14 + Rexp * 2 + Sexp * 2},
  {19293822, FIVEexp*2 + Pexp * 11 + Qexp * 2 + Rexp * 2 + Sexp * 2},
  {-3233411, FIVEexp*4 + Pexp * 8 + Qexp * 4 + Rexp * 2 + Sexp * 2},
  {-11851242, FIVEexp*4 + Pexp * 5 + Qexp * 6 + Rexp * 2 + Sexp * 2},
  {8715248, FIVEexp*6 + Pexp * 2 + Qexp * 8 + Rexp * 2 + Sexp * 2},
  {-40035384, FIVEexp*2 + Pexp * 12 + Rexp * 3 + Sexp * 2},
  {-2844464, FIVEexp*5 + Pexp * 9 + Qexp * 2 + Rexp * 3 + Sexp * 2},
  {-1926509, FIVEexp*6 + Pexp * 6 + Qexp * 4 + Rexp * 3 + Sexp * 2},
  {-20464704, FIVEexp*6 + Pexp * 3 + Qexp * 6 + Rexp * 3 + Sexp * 2},
  {659446, FIVEexp*7 + Qexp * 8 + Rexp * 3 + Sexp * 2},
  {20195832, FIVEexp*4 + Pexp * 10 + Rexp * 4 + Sexp * 2},
  {31050384, FIVEexp*5 + Pexp * 7 + Qexp * 2 + Rexp * 4 + Sexp * 2},
  {6015094, FIVEexp*7 + Pexp * 4 + Qexp * 4 + Rexp * 4 + Sexp * 2},
  {745836, FIVEexp*8 + Pexp + Qexp * 6 + Rexp * 4 + Sexp * 2},
  {-6168448, FIVEexp*6 + Pexp * 8 + Rexp * 5 + Sexp * 2},
  {-7655296, FIVEexp*7 + Pexp * 5 + Qexp * 2 + Rexp * 5 + Sexp * 2},
  {-2982456, FIVEexp*8 + Pexp * 2 + Qexp * 4 + Rexp * 5 + Sexp * 2},
  {228608, FIVEexp*9 + Pexp * 6 + Rexp * 6 + Sexp * 2},
  {845568, FIVEexp*9 + Pexp * 3 + Qexp * 2 + Rexp * 6 + Sexp * 2},
  {80832, FIVEexp*10 + Qexp * 4 + Rexp * 6 + Sexp * 2},
  {-589824, FIVEexp*9 + Pexp * 4 + Rexp * 7 + Sexp * 2},
  {-12288, FIVEexp*11 + Pexp + Qexp * 2 + Rexp * 7 + Sexp * 2},
  {129024, FIVEexp*10 + Pexp * 2 + Rexp * 8 + Sexp * 2},
  {-198288, FIVEexp*3 + Pexp * 14 + Qexp + Sexp * 3},
  {-211356, FIVEexp*5 + Pexp * 11 + Qexp * 3 + Sexp * 3},
  {-1883634, FIVEexp*5 + Pexp * 8 + Qexp * 5 + Sexp * 3},
  {-232156, FIVEexp*7 + Pexp * 5 + Qexp * 7 + Sexp * 3},
  {-65543, FIVEexp*7 + Pexp * 2 + Qexp * 9 + Sexp * 3},
  {6621264, FIVEexp*3 + Pexp * 12 + Qexp + Rexp + Sexp * 3},
  {853956, FIVEexp*6 + Pexp * 9 + Qexp * 3 + Rexp + Sexp * 3},
  {12730262, FIVEexp*5 + Pexp * 6 + Qexp * 5 + Rexp + Sexp * 3},
  {-1424492, FIVEexp*7 + Pexp * 3 + Qexp * 7 + Rexp + Sexp * 3},
  {-197616, FIVEexp*7 + Qexp * 9 + Rexp + Sexp * 3},
  {-23265252, FIVEexp*4 + Pexp * 10 + Qexp + Rexp * 2 + Sexp * 3},
  {-8696944, FIVEexp*6 + Pexp * 7 + Qexp * 3 + Rexp * 2 + Sexp * 3},
  {-90569, FIVEexp*9 + Pexp * 4 + Qexp * 5 + Rexp * 2 + Sexp * 3},
  {-11961122, FIVEexp*7 + Pexp + Qexp * 7 + Rexp * 2 + Sexp * 3},
  {10986832, FIVEexp*6 + Pexp * 8 + Qexp + Rexp * 3 + Sexp * 3},
  {14906088, FIVEexp*7 + Pexp * 5 + Qexp * 3 + Rexp * 3 + Sexp * 3},
  {40862034, FIVEexp*7 + Pexp * 2 + Qexp * 5 + Rexp * 3 + Sexp * 3},
  {-16582096, FIVEexp*7 + Pexp * 6 + Qexp + Rexp * 4 + Sexp * 3},
  {-16744224, FIVEexp*8 + Pexp * 3 + Qexp * 3 + Rexp * 4 + Sexp * 3},
  {-1514408, FIVEexp*9 + Qexp * 5 + Rexp * 4 + Sexp * 3},
  {2747776, FIVEexp*9 + Pexp * 4 + Qexp + Rexp * 5 + Sexp * 3},
  {324096, FIVEexp*10 + Pexp + Qexp * 3 + Rexp * 5 + Sexp * 3},
  {-181504, FIVEexp*11 + Pexp * 2 + Qexp + Rexp * 6 + Sexp * 3},
  {-36864, FIVEexp*11 + Qexp + Rexp * 7 + Sexp * 3},
  {75816, FIVEexp*5 + Pexp * 13 + Sexp * 4},
  {113859, FIVEexp*7 + Pexp * 10 + Qexp * 2 + Sexp * 4},
  {1359114, FIVEexp*7 + Pexp * 7 + Qexp * 4 + Sexp * 4},
  {1214008, FIVEexp*8 + Pexp * 4 + Qexp * 6 + Sexp * 4},
  {2524, FIVEexp*12 + Pexp + Qexp * 8 + Sexp * 4},
  {-2238678, FIVEexp*5 + Pexp * 11 + Rexp + Sexp * 4},
  {-471528, FIVEexp*8 + Pexp * 8 + Qexp * 2 + Rexp + Sexp * 4},
  {-16630457, FIVEexp*7 + Pexp * 5 + Qexp * 4 + Rexp + Sexp * 4},
  {-6337174, FIVEexp*8 + Pexp * 2 + Qexp * 6 + Rexp + Sexp * 4},
  {161406, FIVEexp*8 + Pexp * 9 + Rexp * 2 + Sexp * 4},
  {21078931, FIVEexp*7 + Pexp * 6 + Qexp * 2 + Rexp * 2 + Sexp * 4},
  {4652486, FIVEexp*9 + Pexp * 3 + Qexp * 4 + Rexp * 2 + Sexp * 4},
  {2392792, FIVEexp*9 + Qexp * 6 + Rexp * 2 + Sexp * 4},
  {-2760512, FIVEexp*7 + Pexp * 7 + Rexp * 3 + Sexp * 4},
  {-805334, FIVEexp*10 + Pexp * 4 + Qexp * 2 + Rexp * 3 + Sexp * 4},
  {-122632, FIVEexp*11 + Pexp + Qexp * 4 + Rexp * 3 + Sexp * 4},
  {60912, FIVEexp*10 + Pexp * 5 + Rexp * 4 + Sexp * 4},
  {2845816, FIVEexp*10 + Pexp * 2 + Qexp * 2 + Rexp * 4 + Sexp * 4},
  {-79488, FIVEexp*11 + Pexp * 3 + Rexp * 5 + Sexp * 4},
  {49984, FIVEexp*12 + Qexp * 2 + Rexp * 5 + Sexp * 4},
  {8448, FIVEexp*13 + Pexp + Rexp * 6 + Sexp * 4},
  {-88884, FIVEexp*8 + Pexp * 9 + Qexp + Sexp * 5},
  {-55821, FIVEexp*10 + Pexp * 6 + Qexp * 3 + Sexp * 5},
  {-222866, FIVEexp*10 + Pexp * 3 + Qexp * 5 + Sexp * 5},
  {-2288, FIVEexp*13 + Qexp * 7 + Sexp * 5},
  {-3573, FIVEexp*8 + Pexp * 7 + Qexp + Rexp + Sexp * 5},
  {-40091, FIVEexp*11 + Pexp * 4 + Qexp * 3 + Rexp + Sexp * 5},
  {524558, FIVEexp*10 + Pexp + Qexp * 5 + Rexp + Sexp * 5},
  {439404, FIVEexp*9 + Pexp * 5 + Qexp + Rexp * 2 + Sexp * 5},
  {-861186, FIVEexp*11 + Pexp * 2 + Qexp * 3 + Rexp * 2 + Sexp * 5},
  {145168, FIVEexp*11 + Pexp * 3 + Qexp + Rexp * 3 + Sexp * 5},
  {-106392, FIVEexp*12 + Qexp * 3 + Rexp * 3 + Sexp * 5},
  {-56768, FIVEexp*13 + Pexp + Qexp + Rexp * 4 + Sexp * 5},
  {39393, FIVEexp*10 + Pexp * 8 + Sexp * 6},
  {26228, FIVEexp*12 + Pexp * 5 + Qexp * 2 + Sexp * 6},
  {68578, FIVEexp*12 + Pexp * 2 + Qexp * 4 + Sexp * 6},
  {-355029, FIVEexp*10 + Pexp * 6 + Rexp + Sexp * 6},
  {45096, FIVEexp*12 + Pexp * 3 + Qexp * 2 + Rexp + Sexp * 6},
  {67386, FIVEexp*12 + Qexp * 4 + Rexp + Sexp * 6},
  {28586, FIVEexp*12 + Pexp * 4 + Rexp * 2 + Sexp * 6},
  {595572, FIVEexp*12 + Pexp + Qexp * 2 + Rexp * 2 + Sexp * 6},
  {16696, FIVEexp*12 + Pexp * 2 + Rexp * 3 + Sexp * 6},
  {3648, FIVEexp*14 + Rexp * 4 + Sexp * 6},
  {-17953, FIVEexp*13 + Pexp * 4 + Qexp + Sexp * 7},
  {-2838, FIVEexp*15 + Pexp + Qexp * 3 + Sexp * 7},
  {-18666, FIVEexp*13 + Pexp * 2 + Qexp + Rexp + Sexp * 7},
  {-16952, FIVEexp*14 + Qexp + Rexp * 2 + Sexp * 7},
  {616, FIVEexp*15 + Pexp * 3 + Sexp * 8},
  {139, FIVEexp*17 + Qexp * 2 + Sexp * 8},
  {582, FIVEexp*15 + Pexp + Rexp + Sexp * 8},
  {0, END_COEFF},

  // o[2]
  {112, FIVEexp*2 + Pexp * 9 + Qexp * 8},
  {2228, FIVEexp*2 + Pexp * 6 + Qexp * 10},
  {14544, FIVEexp*2 + Pexp * 3 + Qexp * 12},
  {31104, FIVEexp*2 + Qexp * 14},
  {-1088, FIVEexp*2 + Pexp * 10 + Qexp * 6 + Rexp},
  {-28008, FIVEexp*2 + Pexp * 7 + Qexp * 8 + Rexp},
  {-229062, FIVEexp*2 + Pexp * 4 + Qexp * 10 + Rexp},
  {-120528, FIVEexp*3 + Pexp + Qexp * 12 + Rexp},
  {2988, FIVEexp*2 + Pexp * 11 + Qexp * 4 + Rexp * 2},
  {114383, FIVEexp*2 + Pexp * 8 + Qexp * 6 + Rexp * 2},
  {1247029, FIVEexp*2 + Pexp * 5 + Qexp * 8 + Rexp * 2},
  {4125906, FIVEexp*2 + Pexp * 2 + Qexp * 10 + Rexp * 2},
  {-324, FIVEexp*3 + Pexp * 12 + Qexp * 2 + Rexp * 3},
  {-170976, FIVEexp*2 + Pexp * 9 + Qexp * 4 + Rexp * 3},
  {-3042633, FIVEexp*2 + Pexp * 6 + Qexp * 6 + Rexp * 3},
  {-584998, FIVEexp*4 + Pexp * 3 + Qexp * 8 + Rexp * 3},
  {-1058112, FIVEexp*3 + Qexp * 10 + Rexp * 3},
  {87696, FIVEexp*2 + Pexp * 10 + Qexp * 2 + Rexp * 4},
  {5924, FIVEexp*6 + Pexp * 7 + Qexp * 4 + Rexp * 4},
  {6393551, FIVEexp*3 + Pexp * 4 + Qexp * 6 + Rexp * 4},
  {1901109, FIVEexp*4 + Pexp + Qexp * 8 + Rexp * 4},
  {-329852, FIVEexp*3 + Pexp * 8 + Qexp * 2 + Rexp * 5},
  {-1462736, FIVEexp*4 + Pexp * 5 + Qexp * 4 + Rexp * 5},
  {-1062033, FIVEexp*5 + Pexp * 2 + Qexp * 6 + Rexp * 5},
  {127632, FIVEexp*5 + Pexp * 6 + Qexp * 2 + Rexp * 6},
  {252416, FIVEexp*6 + Pexp * 3 + Qexp * 4 + Rexp * 6},
  {141524, FIVEexp*6 + Qexp * 6 + Rexp * 6},
  {-116288, FIVEexp*6 + Pexp * 4 + Qexp * 2 + Rexp * 7},
  {-86016, FIVEexp*7 + Pexp + Qexp * 4 + Rexp * 7},
  {9984, FIVEexp*8 + Pexp * 2 + Qexp * 2 + Rexp * 8},
  {3024, FIVEexp*2 + Pexp * 11 + Qexp * 5 + Sexp},
  {72924, FIVEexp*2 + Pexp * 8 + Qexp * 7 + Sexp},
  {581366, FIVEexp*2 + Pexp * 5 + Qexp * 9 + Sexp},
  {306126, FIVEexp*3 + Pexp * 2 + Qexp * 11 + Sexp},
  {-15768, FIVEexp*2 + Pexp * 12 + Qexp * 3 + Rexp + Sexp},
  {-458154, FIVEexp*2 + Pexp * 9 + Qexp * 5 + Rexp + Sexp},
  {-809704, FIVEexp*3 + Pexp * 6 + Qexp * 7 + Rexp + Sexp},
  {-71541, FIVEexp*5 + Pexp * 3 + Qexp * 9 + Rexp + Sexp},
  {26568, FIVEexp*6 + Qexp * 11 + Rexp + Sexp},
  {1944, FIVEexp*3 + Pexp * 13 + Qexp + Rexp * 2 + Sexp},
  {546183, FIVEexp*2 + Pexp * 10 + Qexp * 3 + Rexp * 2 + Sexp},
  {6552469, FIVEexp*2 + Pexp * 7 + Qexp * 5 + Rexp * 2 + Sexp},
  {277538, FIVEexp*4 + Pexp * 4 + Qexp * 7 + Rexp * 2 + Sexp},
  {-38511, FIVEexp*7 + Pexp + Qexp * 9 + Rexp * 2 + Sexp},
  {-80676, FIVEexp*2 + Pexp * 11 + Qexp + Rexp * 3 + Sexp},
  {-138522, FIVEexp*4 + Pexp * 8 + Qexp * 3 + Rexp * 3 + Sexp},
  {-2593173, FIVEexp*3 + Pexp * 5 + Qexp * 5 + Rexp * 3 + Sexp},
  {1081122, FIVEexp*5 + Pexp * 2 + Qexp * 7 + Rexp * 3 + Sexp},
  {-713688, FIVEexp*3 + Pexp * 9 + Qexp + Rexp * 4 + Sexp},
  {-88332, FIVEexp*4 + Pexp * 6 + Qexp * 3 + Rexp * 4 + Sexp},
  {478064, FIVEexp*5 + Pexp * 3 + Qexp * 5 + Rexp * 4 + Sexp},
  {-804711, FIVEexp*6 + Qexp * 7 + Rexp * 4 + Sexp},
  {364832, FIVEexp*5 + Pexp * 7 + Qexp + Rexp * 5 + Sexp},
  {544, FIVEexp*7 + Pexp * 4 + Qexp * 3 + Rexp * 5 + Sexp},
  {275344, FIVEexp*7 + Pexp + Qexp * 5 + Rexp * 5 + Sexp},
  {-259712, FIVEexp*6 + Pexp * 5 + Qexp + Rexp * 6 + Sexp},
  {17216, FIVEexp*8 + Pexp * 2 + Qexp * 3 + Rexp * 6 + Sexp},
  {-3584, FIVEexp*8 + Pexp * 3 + Qexp + Rexp * 7 + Sexp},
  {-19968, FIVEexp*9 + Qexp * 3 + Rexp * 7 + Sexp},
  {20412, FIVEexp*2 + Pexp * 13 + Qexp * 2 + Sexp * 2},
  {192591, FIVEexp*2 + Pexp * 10 + Qexp * 4 + Sexp * 2},
  {-562121, FIVEexp*3 + Pexp * 7 + Qexp * 6 + Sexp * 2},
  {-13011, FIVEexp*7 + Pexp * 4 + Qexp * 8 + Sexp * 2},
  {-1030752, FIVEexp*5 + Pexp + Qexp * 10 + Sexp * 2},
  {-2916, FIVEexp*3 + Pexp * 14 + Rexp + Sexp * 2},
  {242514, FIVEexp*3 + Pexp * 11 + Qexp * 2 + Rexp + Sexp * 2},
  {70833, FIVEexp*6 + Pexp * 8 + Qexp * 4 + Rexp + Sexp * 2},
  {702989, FIVEexp*6 + Pexp * 5 + Qexp * 6 + Rexp + Sexp * 2},
  {2169972, FIVEexp*6 + Pexp * 2 + Qexp * 8 + Rexp + Sexp * 2},
  {-1519236, FIVEexp*2 + Pexp * 12 + Rexp * 2 + Sexp * 2},
  {-3428649, FIVEexp*4 + Pexp * 9 + Qexp * 2 + Rexp * 2 + Sexp * 2},
  {-1721352, FIVEexp*6 + Pexp * 6 + Qexp * 4 + Rexp * 2 + Sexp * 2},
  {-1223057, FIVEexp*7 + Pexp * 3 + Qexp * 6 + Rexp * 2 + Sexp * 2},
  {5796, FIVEexp*9 + Qexp * 8 + Rexp * 2 + Sexp * 2},
  {13950252, FIVEexp*3 + Pexp * 10 + Rexp * 3 + Sexp * 2},
  {11338324, FIVEexp*5 + Pexp * 7 + Qexp * 2 + Rexp * 3 + Sexp * 2},
  {455311, FIVEexp*8 + Pexp * 4 + Qexp * 4 + Rexp * 3 + Sexp * 2},
  {310427, FIVEexp*8 + Pexp + Qexp * 6 + Rexp * 3 + Sexp * 2},
  {-8301812, FIVEexp*5 + Pexp * 8 + Rexp * 4 + Sexp * 2},
  {-15946592, FIVEexp*6 + Pexp * 5 + Qexp * 2 + Rexp * 4 + Sexp * 2},
  {-1182053, FIVEexp*8 + Pexp * 2 + Qexp * 4 + Rexp * 4 + Sexp * 2},
  {2276144, FIVEexp*7 + Pexp * 6 + Rexp * 5 + Sexp * 2},
  {1859136, FIVEexp*8 + Pexp * 3 + Qexp * 2 + Rexp * 5 + Sexp * 2},
  {206948, FIVEexp*9 + Qexp * 4 + Rexp * 5 + Sexp * 2},
  {-289728, FIVEexp*9 + Pexp * 4 + Rexp * 6 + Sexp * 2},
  {-41728, FIVEexp*10 + Pexp + Qexp * 2 + Rexp * 6 + Sexp * 2},
  {69888, FIVEexp*10 + Pexp * 2 + Rexp * 7 + Sexp * 2},
  {-473121, FIVEexp*3 + Pexp * 12 + Qexp + Sexp * 3},
  {-82593, FIVEexp*6 + Pexp * 9 + Qexp * 3 + Sexp * 3},
  {-2811223, FIVEexp*5 + Pexp * 6 + Qexp * 5 + Sexp * 3},
  {-39906, FIVEexp*8 + Pexp * 3 + Qexp * 7 + Sexp * 3},
  {224064, FIVEexp*7 + Qexp * 9 + Sexp * 3},
  {57348, FIVEexp*6 + Pexp * 10 + Qexp + Rexp + Sexp * 3},
  {165744, FIVEexp*6 + Pexp * 7 + Qexp * 3 + Rexp + Sexp * 3},
  {-1058562, FIVEexp*7 + Pexp * 4 + Qexp * 5 + Rexp + Sexp * 3},
  {-178014, FIVEexp*9 + Pexp + Qexp * 7 + Rexp + Sexp * 3},
  {43578, FIVEexp*8 + Pexp * 8 + Qexp + Rexp * 2 + Sexp * 3},
  {4107822, FIVEexp*7 + Pexp * 5 + Qexp * 3 + Rexp * 2 + Sexp * 3},
  {13339096, FIVEexp*7 + Pexp * 2 + Qexp * 5 + Rexp * 2 + Sexp * 3},
  {-4521764, FIVEexp*7 + Pexp * 6 + Qexp + Rexp * 3 + Sexp * 3},
  {-1132884, FIVEexp*9 + Pexp * 3 + Qexp * 3 + Rexp * 3 + Sexp * 3},
  {-3206167, FIVEexp*8 + Qexp * 5 + Rexp * 3 + Sexp * 3},
  {5036448, FIVEexp*8 + Pexp * 4 + Qexp + Rexp * 4 + Sexp * 3},
  {810672, FIVEexp*9 + Pexp + Qexp * 3 + Rexp * 4 + Sexp * 3},
  {-364608, FIVEexp*10 + Pexp * 2 + Qexp + Rexp * 5 + Sexp * 3},
  {-19968, FIVEexp*11 + Qexp + Rexp * 6 + Sexp * 3},
  {147987, FIVEexp*5 + Pexp * 11 + Sexp * 4},
  {181899, FIVEexp*7 + Pexp * 8 + Qexp * 2 + Sexp * 4},
  {1273878, FIVEexp*7 + Pexp * 5 + Qexp * 4 + Sexp * 4},
  {294286, FIVEexp*8 + Pexp * 2 + Qexp * 6 + Sexp * 4},
  {-226827, FIVEexp*7 + Pexp * 9 + Rexp + Sexp * 4},
  {-256821, FIVEexp*8 + Pexp * 6 + Qexp * 2 + Rexp + Sexp * 4},
  {104647, FIVEexp*10 + Pexp * 3 + Qexp * 4 + Rexp + Sexp * 4},
  {122184, FIVEexp*10 + Qexp * 6 + Rexp + Sexp * 4},
  {2193552, FIVEexp*7 + Pexp * 7 + Rexp * 2 + Sexp * 4},
  {-570056, FIVEexp*9 + Pexp * 4 + Qexp * 2 + Rexp * 2 + Sexp * 4},
  {-24811, FIVEexp*11 + Pexp + Qexp * 4 + Rexp * 2 + Sexp * 4},
  {-701584, FIVEexp*8 + Pexp * 5 + Rexp * 3 + Sexp * 4},
  {848401, FIVEexp*10 + Pexp * 2 + Qexp * 2 + Rexp * 3 + Sexp * 4},
  {-212928, FIVEexp*10 + Pexp * 3 + Rexp * 4 + Sexp * 4},
  {105212, FIVEexp*11 + Qexp * 2 + Rexp * 4 + Sexp * 4},
  {25344, FIVEexp*12 + Pexp + Rexp * 5 + Sexp * 4},
  {-347328, FIVEexp*8 + Pexp * 7 + Qexp + Sexp * 5},
  {-26117, FIVEexp*11 + Pexp * 4 + Qexp * 3 + Sexp * 5},
  {-2907, FIVEexp*10 + Pexp + Qexp * 5 + Sexp * 5},
  {103953, FIVEexp*10 + Pexp * 5 + Qexp + Rexp + Sexp * 5},
  {-216769, FIVEexp*11 + Pexp * 2 + Qexp * 3 + Rexp + Sexp * 5},
  {2916, FIVEexp*13 + Pexp * 3 + Qexp + Rexp * 2 + Sexp * 5},
  {-5677, FIVEexp*13 + Qexp * 3 + Rexp * 2 + Sexp * 5},
  {-26464, FIVEexp*13 + Pexp + Qexp + Rexp * 3 + Sexp * 5},
  {44091, FIVEexp*10 + Pexp * 6 + Sexp * 6},
  {38494, FIVEexp*12 + Pexp * 3 + Qexp * 2 + Sexp * 6},
  {2356, FIVEexp*12 + Qexp * 4 + Sexp * 6},
  {-16738, FIVEexp*12 + Pexp * 4 + Rexp + Sexp * 6},
  {34557, FIVEexp*13 + Pexp + Qexp * 2 + Rexp + Sexp * 6},
  {33409, FIVEexp*12 + Pexp * 2 + Rexp * 2 + Sexp * 6},
  {10332, FIVEexp*13 + Rexp * 3 + Sexp * 6},
  {-20901, FIVEexp*13 + Pexp * 2 + Qexp + Sexp * 7},
  {-1213, FIVEexp*15 + Qexp + Rexp + Sexp * 7},
  {622, FIVEexp*15 + Pexp + Sexp * 8},
  {0, END_COEFF},

  // o[3]
  {64, FIVEexp*2 + Pexp * 10 + Qexp * 6},
  {368, FIVEexp*2 + Pexp * 7 + Qexp * 8},
  {-1008, FIVEexp*3 + Pexp * 4 + Qexp * 10},
  {-31104, FIVEexp*2 + Pexp + Qexp * 12},
  {-576, FIVEexp*2 + Pexp * 11 + Qexp * 4 + Rexp},
  {-4772, FIVEexp*2 + Pexp * 8 + Qexp * 6 + Rexp},
  {48129, FIVEexp*2 + Pexp * 5 + Qexp * 8 + Rexp},
  {376488, FIVEexp*2 + Pexp * 2 + Qexp * 10 + Rexp},
  {1296, FIVEexp*2 + Pexp * 12 + Qexp * 2 + Rexp * 2},
  {16718, FIVEexp*2 + Pexp * 9 + Qexp * 4 + Rexp * 2},
  {-181749, FIVEexp*2 + Pexp * 6 + Qexp * 6 + Rexp * 2},
  {-78413, FIVEexp*4 + Pexp * 3 + Qexp * 8 + Rexp * 2},
  {-193536, FIVEexp*3 + Qexp * 10 + Rexp * 2},
  {-11682, FIVEexp*2 + Pexp * 10 + Qexp * 2 + Rexp * 3},
  {14016, FIVEexp*4 + Pexp * 7 + Qexp * 4 + Rexp * 3},
  {1100053, FIVEexp*3 + Pexp * 4 + Qexp * 6 + Rexp * 3},
  {360702, FIVEexp*4 + Pexp + Qexp * 8 + Rexp * 3},
  {-33706, FIVEexp*3 + Pexp * 8 + Qexp * 2 + Rexp * 4},
  {-277753, FIVEexp*4 + Pexp * 5 + Qexp * 4 + Rexp * 4},
  {-208961, FIVEexp*5 + Pexp * 2 + Qexp * 6 + Rexp * 4},
  {26424, FIVEexp*5 + Pexp * 6 + Qexp * 2 + Rexp * 5},
  {10728, FIVEexp*7 + Pexp * 3 + Qexp * 4 + Rexp * 5},
  {37028, FIVEexp*6 + Qexp * 6 + Rexp * 5},
  {-26976, FIVEexp*6 + Pexp * 4 + Qexp * 2 + Rexp * 6},
  {-22992, FIVEexp*7 + Pexp + Qexp * 4 + Rexp * 6},
  {2688, FIVEexp*8 + Pexp * 2 + Qexp * 2 + Rexp * 7},
  {1728, FIVEexp*2 + Pexp * 12 + Qexp * 3 + Sexp},
  {32292, FIVEexp*2 + Pexp * 9 + Qexp * 5 + Sexp},
  {213129, FIVEexp*2 + Pexp * 6 + Qexp * 7 + Sexp},
  {27114, FIVEexp*4 + Pexp * 3 + Qexp * 9 + Sexp},
  {47304, FIVEexp*4 + Qexp * 11 + Sexp},
  {-7776, FIVEexp*2 + Pexp * 13 + Qexp + Rexp + Sexp},
  {-220212, FIVEexp*2 + Pexp * 10 + Qexp * 3 + Rexp + Sexp},
  {-1995468, FIVEexp*2 + Pexp * 7 + Qexp * 5 + Rexp + Sexp},
  {-57223, FIVEexp*5 + Pexp * 4 + Qexp * 7 + Rexp + Sexp},
  {-14256, FIVEexp*6 + Pexp + Qexp * 9 + Rexp + Sexp},
  {272592, FIVEexp*2 + Pexp * 11 + Qexp + Rexp * 2 + Sexp},
  {192841, FIVEexp*4 + Pexp * 8 + Qexp * 3 + Rexp * 2 + Sexp},
  {4213556, FIVEexp*3 + Pexp * 5 + Qexp * 5 + Rexp * 2 + Sexp},
  {26901, FIVEexp*5 + Pexp * 2 + Qexp * 7 + Rexp * 2 + Sexp},
  {-989364, FIVEexp*3 + Pexp * 9 + Qexp + Rexp * 3 + Sexp},
  {-1770771, FIVEexp*4 + Pexp * 6 + Qexp * 3 + Rexp * 3 + Sexp},
  {-214258, FIVEexp*5 + Pexp * 3 + Qexp * 5 + Rexp * 3 + Sexp},
  {-183633, FIVEexp*6 + Qexp * 7 + Rexp * 3 + Sexp},
  {321392, FIVEexp*5 + Pexp * 7 + Qexp + Rexp * 4 + Sexp},
  {8664, FIVEexp*8 + Pexp * 4 + Qexp * 3 + Rexp * 4 + Sexp},
  {14504, FIVEexp*8 + Pexp + Qexp * 5 + Rexp * 4 + Sexp},
  {-218304, FIVEexp*6 + Pexp * 5 + Qexp + Rexp * 5 + Sexp},
  {-1072, FIVEexp*8 + Pexp * 2 + Qexp * 3 + Rexp * 5 + Sexp},
  {4352, FIVEexp*8 + Pexp * 3 + Qexp + Rexp * 6 + Sexp},
  {-5376, FIVEexp*9 + Qexp * 3 + Rexp * 6 + Sexp},
  {11664, FIVEexp*2 + Pexp * 14 + Sexp * 2},
  {393174, FIVEexp*2 + Pexp * 11 + Qexp * 2 + Sexp * 2},
  {182643, FIVEexp*4 + Pexp * 8 + Qexp * 4 + Sexp * 2},
  {835471, FIVEexp*4 + Pexp * 5 + Qexp * 6 + Sexp * 2},
  {45882, FIVEexp*6 + Pexp * 2 + Qexp * 8 + Sexp * 2},
  {-745038, FIVEexp*2 + Pexp * 12 + Rexp + Sexp * 2},
  {-24813, FIVEexp*6 + Pexp * 9 + Qexp * 2 + Rexp + Sexp * 2},
  {-131591, FIVEexp*6 + Pexp * 6 + Qexp * 4 + Rexp + Sexp * 2},
  {-48653, FIVEexp*6 + Pexp * 3 + Qexp * 6 + Rexp + Sexp * 2},
  {39312, FIVEexp*7 + Qexp * 8 + Rexp + Sexp * 2},
  {4099356, FIVEexp*3 + Pexp * 10 + Rexp * 2 + Sexp * 2},
  {1874897, FIVEexp*5 + Pexp * 7 + Qexp * 2 + Rexp * 2 + Sexp * 2},
  {155572, FIVEexp*7 + Pexp * 4 + Qexp * 4 + Rexp * 2 + Sexp * 2},
  {8157, FIVEexp*9 + Pexp + Qexp * 6 + Rexp * 2 + Sexp * 2},
  {-2111486, FIVEexp*5 + Pexp * 8 + Rexp * 3 + Sexp * 2},
  {-2787176, FIVEexp*6 + Pexp * 5 + Qexp * 2 + Rexp * 3 + Sexp * 2},
  {-209949, FIVEexp*8 + Pexp * 2 + Qexp * 4 + Rexp * 3 + Sexp * 2},
  {557288, FIVEexp*7 + Pexp * 6 + Rexp * 4 + Sexp * 2},
  {411728, FIVEexp*8 + Pexp * 3 + Qexp * 2 + Rexp * 4 + Sexp * 2},
  {50724, FIVEexp*9 + Qexp * 4 + Rexp * 4 + Sexp * 2},
  {-72608, FIVEexp*9 + Pexp * 4 + Rexp * 5 + Sexp * 2},
  {-12032, FIVEexp*10 + Pexp + Qexp * 2 + Rexp * 5 + Sexp * 2},
  {18816, FIVEexp*10 + Pexp * 2 + Rexp * 6 + Sexp * 2},
  {-247374, FIVEexp*4 + Pexp * 10 + Qexp + Sexp * 3},
  {-211806, FIVEexp*6 + Pexp * 7 + Qexp * 3 + Sexp * 3},
  {-266677, FIVEexp*7 + Pexp * 4 + Qexp * 5 + Sexp * 3},
  {-444564, FIVEexp*7 + Pexp + Qexp * 7 + Sexp * 3},
  {382104, FIVEexp*6 + Pexp * 8 + Qexp + Rexp + Sexp * 3},
  {877226, FIVEexp*7 + Pexp * 5 + Qexp * 3 + Rexp + Sexp * 3},
  {1826538, FIVEexp*7 + Pexp * 2 + Qexp * 5 + Rexp + Sexp * 3},
  {-953717, FIVEexp*7 + Pexp * 6 + Qexp + Rexp * 2 + Sexp * 3},
  {-199244, FIVEexp*9 + Pexp * 3 + Qexp * 3 + Rexp * 2 + Sexp * 3},
  {-693601, FIVEexp*8 + Qexp * 5 + Rexp * 2 + Sexp * 3},
  {938744, FIVEexp*8 + Pexp * 4 + Qexp + Rexp * 3 + Sexp * 3},
  {204816, FIVEexp*9 + Pexp + Qexp * 3 + Rexp * 3 + Sexp * 3},
  {-15056, FIVEexp*11 + Pexp * 2 + Qexp + Rexp * 4 + Sexp * 3},
  {-5376, FIVEexp*11 + Qexp + Rexp * 5 + Sexp * 3},
  {4752, FIVEexp*7 + Pexp * 9 + Sexp * 4},
  {272832, FIVEexp*7 + Pexp * 6 + Qexp * 2 + Sexp * 4},
  {91912, FIVEexp*9 + Pexp * 3 + Qexp * 4 + Sexp * 4},
  {122064, FIVEexp*9 + Qexp * 6 + Sexp * 4},
  {-116109, FIVEexp*7 + Pexp * 7 + Rexp + Sexp * 4},
  {-5508, FIVEexp*11 + Pexp * 4 + Qexp * 2 + Rexp + Sexp * 4},
  {-7778, FIVEexp*11 + Pexp + Qexp * 4 + Rexp + Sexp * 4},
  {359073, FIVEexp*8 + Pexp * 5 + Rexp * 2 + Sexp * 4},
  {152193, FIVEexp*10 + Pexp * 2 + Qexp * 2 + Rexp * 2 + Sexp * 4},
  {-82584, FIVEexp*10 + Pexp * 3 + Rexp * 3 + Sexp * 4},
  {21836, FIVEexp*11 + Qexp * 2 + Rexp * 3 + Sexp * 4},
  {7056, FIVEexp*12 + Pexp + Rexp * 4 + Sexp * 4},
  {-64332, FIVEexp*9 + Pexp * 5 + Qexp + Sexp * 5},
  {-25319, FIVEexp*11 + Pexp * 2 + Qexp * 3 + Sexp * 5},
  {25966, FIVEexp*11 + Pexp * 3 + Qexp + Rexp + Sexp * 5},
  {-4279, FIVEexp*12 + Qexp * 3 + Rexp + Sexp * 5},
  {-5432, FIVEexp*13 + Pexp + Qexp + Rexp * 2 + Sexp * 5},
  {659, FIVEexp*12 + Pexp * 4 + Sexp * 6},
  {26019, FIVEexp*12 + Pexp + Qexp * 2 + Sexp * 6},
  {-2903, FIVEexp*12 + Pexp * 2 + Rexp + Sexp * 6},
  {3196, FIVEexp*13 + Rexp * 2 + Sexp * 6},
  {-1319, FIVEexp*14 + Qexp + Sexp * 7},
  {0, END_COEFF},

  // o[4]
  {112, FIVEexp*2 + Pexp * 8 + Qexp * 6},
  {1652, FIVEexp*2 + Pexp * 5 + Qexp * 8},
  {6048, FIVEexp*2 + Pexp * 2 + Qexp * 10},
  {-1008, FIVEexp*2 + Pexp * 9 + Qexp * 4 + Rexp},
  {-21704, FIVEexp*2 + Pexp * 6 + Qexp * 6 + Rexp},
  {-27183, FIVEexp*3 + Pexp * 3 + Qexp * 8 + Rexp},
  {-46008, FIVEexp*3 + Qexp * 10 + Rexp},
  {2268, FIVEexp*2 + Pexp * 10 + Qexp * 2 + Rexp * 2},
  {15777, FIVEexp*3 + Pexp * 7 + Qexp * 4 + Rexp * 2},
  {148994, FIVEexp*3 + Pexp * 4 + Qexp * 6 + Rexp * 2},
  {80406, FIVEexp*4 + Pexp + Qexp * 8 + Rexp * 2},
  {-13608, FIVEexp*3 + Pexp * 8 + Qexp * 2 + Rexp * 3},
  {-52209, FIVEexp*4 + Pexp * 5 + Qexp * 4 + Rexp * 3},
  {-44758, FIVEexp*5 + Pexp * 2 + Qexp * 6 + Rexp * 3},
  {5812, FIVEexp*5 + Pexp * 6 + Qexp * 2 + Rexp * 4},
  {11336, FIVEexp*6 + Pexp * 3 + Qexp * 4 + Rexp * 4},
  {7791, FIVEexp*6 + Qexp * 6 + Rexp * 4},
  {-6432, FIVEexp*6 + Pexp * 4 + Qexp * 2 + Rexp * 5},
  {-4944, FIVEexp*7 + Pexp + Qexp * 4 + Rexp * 5},
  {576, FIVEexp*8 + Pexp * 2 + Qexp * 2 + Rexp * 6},
  {3024, FIVEexp*2 + Pexp * 10 + Qexp * 3 + Sexp},
  {68352, FIVEexp*2 + Pexp * 7 + Qexp * 5 + Sexp},
  {20539, FIVEexp*4 + Pexp * 4 + Qexp * 7 + Sexp},
  {2052, FIVEexp*6 + Pexp + Qexp * 9 + Sexp},
  {-13608, FIVEexp*2 + Pexp * 11 + Qexp + Rexp + Sexp},
  {-81486, FIVEexp*3 + Pexp * 8 + Qexp * 3 + Rexp + Sexp},
  {-780022, FIVEexp*3 + Pexp * 5 + Qexp * 5 + Rexp + Sexp},
  {-3861, FIVEexp*7 + Pexp * 2 + Qexp * 7 + Rexp + Sexp},
  {57348, FIVEexp*3 + Pexp * 9 + Qexp + Rexp * 2 + Sexp},
  {217537, FIVEexp*4 + Pexp * 6 + Qexp * 3 + Rexp * 2 + Sexp},
  {187991, FIVEexp*5 + Pexp * 3 + Qexp * 5 + Rexp * 2 + Sexp},
  {-24624, FIVEexp*6 + Qexp * 7 + Rexp * 2 + Sexp},
  {-9384, FIVEexp*5 + Pexp * 7 + Qexp + Rexp * 3 + Sexp},
  {-4104, FIVEexp*7 + Pexp * 4 + Qexp * 3 + Rexp * 3 + Sexp},
  {268, FIVEexp*9 + Pexp + Qexp * 5 + Rexp * 3 + Sexp},
  {-2688, FIVEexp*6 + Pexp * 5 + Qexp + Rexp * 4 + Sexp},
  {176, FIVEexp*9 + Pexp * 2 + Qexp * 3 + Rexp * 4 + Sexp},
  {384, FIVEexp*8 + Pexp * 3 + Qexp + Rexp * 5 + Sexp},
  {-1152, FIVEexp*9 + Qexp * 3 + Rexp * 5 + Sexp},
  {20412, FIVEexp*2 + Pexp * 12 + Sexp * 2},
  {20493, FIVEexp*4 + Pexp * 9 + Qexp * 2 + Sexp * 2},
  {6852, FIVEexp*6 + Pexp * 6 + Qexp * 4 + Sexp * 2},
  {692, FIVEexp*8 + Pexp * 3 + Qexp * 6 + Sexp * 2},
  {-432, FIVEexp*8 + Qexp * 8 + Sexp * 2},
  {-20412, FIVEexp*3 + Pexp * 10 + Rexp + Sexp * 2},
  {-324, FIVEexp*6 + Pexp * 7 + Qexp * 2 + Rexp + Sexp * 2},
  {9118, FIVEexp*7 + Pexp * 4 + Qexp * 4 + Rexp + Sexp * 2},
  {12258, FIVEexp*8 + Pexp + Qexp * 6 + Rexp + Sexp * 2},
  {-82188, FIVEexp*5 + Pexp * 8 + Rexp * 2 + Sexp * 2},
  {-228788, FIVEexp*6 + Pexp * 5 + Qexp * 2 + Rexp * 2 + Sexp * 2},
  {-38318, FIVEexp*8 + Pexp * 2 + Qexp * 4 + Rexp * 2 + Sexp * 2},
  {51724, FIVEexp*7 + Pexp * 6 + Rexp * 3 + Sexp * 2},
  {67584, FIVEexp*8 + Pexp * 3 + Qexp * 2 + Rexp * 3 + Sexp * 2},
  {8747, FIVEexp*9 + Qexp * 4 + Rexp * 3 + Sexp * 2},
  {-11168, FIVEexp*9 + Pexp * 4 + Rexp * 4 + Sexp * 2},
  {-2496, FIVEexp*10 + Pexp + Qexp * 2 + Rexp * 4 + Sexp * 2},
  {4032, FIVEexp*10 + Pexp * 2 + Rexp * 5 + Sexp * 2},
  {-1629, FIVEexp*7 + Pexp * 5 + Qexp * 3 + Sexp * 3},
  {-8707, FIVEexp*7 + Pexp * 2 + Qexp * 5 + Sexp * 3},
  {-26217, FIVEexp*7 + Pexp * 6 + Qexp + Rexp + Sexp * 3},
  {-48114, FIVEexp*8 + Pexp * 3 + Qexp * 3 + Rexp + Sexp * 3},
  {-64128, FIVEexp*8 + Qexp * 5 + Rexp + Sexp * 3},
  {68152, FIVEexp*8 + Pexp * 4 + Qexp + Rexp * 2 + Sexp * 3},
  {18948, FIVEexp*9 + Pexp + Qexp * 3 + Rexp * 2 + Sexp * 3},
  {-432, FIVEexp*12 + Pexp * 2 + Qexp + Rexp * 3 + Sexp * 3},
  {-1152, FIVEexp*11 + Qexp + Rexp * 4 + Sexp * 3},
  {14391, FIVEexp*7 + Pexp * 7 + Sexp * 4},
  {4737, FIVEexp*9 + Pexp * 4 + Qexp * 2 + Sexp * 4},
  {-164, FIVEexp*11 + Pexp + Qexp * 4 + Sexp * 4},
  {-10251, FIVEexp*8 + Pexp * 5 + Rexp + Sexp * 4},
  {737, FIVEexp*12 + Pexp * 2 + Qexp * 2 + Rexp + Sexp * 4},
  {-7752, FIVEexp*10 + Pexp * 3 + Rexp * 2 + Sexp * 4},
  {2733, FIVEexp*11 + Qexp * 2 + Rexp * 2 + Sexp * 4},
  {1488, FIVEexp*12 + Pexp + Rexp * 3 + Sexp * 4},
  {-3, FIVEexp*15 + Pexp * 3 + Qexp + Sexp * 5},
  {-16, FIVEexp*13 + Qexp * 3 + Sexp * 5},
  {-616, FIVEexp*13 + Pexp + Qexp + Rexp + Sexp * 5},
  {647, FIVEexp*12 + Pexp * 2 + Sexp * 6},
  {613, FIVEexp*13 + Rexp + Sexp * 6},
  {0, END_COEFF},

  // o[5]
  {-576, FIVEexp*2 + Pexp * 6 + Qexp * 6},
  {-8496, FIVEexp*2 + Pexp * 3 + Qexp * 8},
  {-31104, FIVEexp*2 + Qexp * 10},
  {3684, FIVEexp*2 + Pexp * 7 + Qexp * 4 + Rexp},
  {67587, FIVEexp*2 + Pexp * 4 + Qexp * 6 + Rexp},
  {58968, FIVEexp*3 + Pexp + Qexp * 8 + Rexp},
  {-4914, FIVEexp*2 + Pexp * 8 + Qexp * 2 + Rexp * 2},
  {-29882, FIVEexp*3 + Pexp * 5 + Qexp * 4 + Rexp * 2},
  {-35892, FIVEexp*4 + Pexp * 2 + Qexp * 6 + Rexp * 2},
  {3678, FIVEexp*4 + Pexp * 6 + Qexp * 2 + Rexp * 3},
  {1881, FIVEexp*6 + Pexp * 3 + Qexp * 4 + Rexp * 3},
  {1152, FIVEexp*6 + Qexp * 6 + Rexp * 3},
  {-1136, FIVEexp*6 + Pexp * 4 + Qexp * 2 + Rexp * 4},
  {-804, FIVEexp*7 + Pexp + Qexp * 4 + Rexp * 4},
  {96, FIVEexp*8 + Pexp * 2 + Qexp * 2 + Rexp * 5},
  {-2052, FIVEexp*2 + Pexp * 8 + Qexp * 3 + Sexp},
  {-30721, FIVEexp*2 + Pexp * 5 + Qexp * 5 + Sexp},
  {-4482, FIVEexp*4 + Pexp * 2 + Qexp * 7 + Sexp},
  {-11016, FIVEexp*2 + Pexp * 9 + Qexp + Rexp + Sexp},
  {-43839, FIVEexp*3 + Pexp * 6 + Qexp * 3 + Rexp + Sexp},
  {-56862, FIVEexp*4 + Pexp * 3 + Qexp * 5 + Rexp + Sexp},
  {-22032, FIVEexp*5 + Qexp * 7 + Rexp + Sexp},
  {20412, FIVEexp*4 + Pexp * 7 + Qexp + Rexp * 2 + Sexp},
  {8553, FIVEexp*6 + Pexp * 4 + Qexp * 3 + Rexp * 2 + Sexp},
  {14247, FIVEexp*6 + Pexp + Qexp * 5 + Rexp * 2 + Sexp},
  {-6944, FIVEexp*6 + Pexp * 5 + Qexp + Rexp * 3 + Sexp},
  {-3716, FIVEexp*7 + Pexp * 2 + Qexp * 3 + Rexp * 3 + Sexp},
  {704, FIVEexp*8 + Pexp * 3 + Qexp + Rexp * 4 + Sexp},
  {-192, FIVEexp*9 + Qexp * 3 + Rexp * 4 + Sexp},
  {77274, FIVEexp*2 + Pexp * 10 + Sexp * 2},
  {64341, FIVEexp*4 + Pexp * 7 + Qexp * 2 + Sexp * 2},
  {16251, FIVEexp*6 + Pexp * 4 + Qexp * 4 + Sexp * 2},
  {1188, FIVEexp*8 + Pexp + Qexp * 6 + Sexp * 2},
  {-113724, FIVEexp*4 + Pexp * 8 + Rexp + Sexp * 2},
  {-261999, FIVEexp*5 + Pexp * 5 + Qexp * 2 + Rexp + Sexp * 2},
  {-24093, FIVEexp*7 + Pexp * 2 + Qexp * 4 + Rexp + Sexp * 2},
  {52866, FIVEexp*6 + Pexp * 6 + Rexp * 2 + Sexp * 2},
  {55932, FIVEexp*7 + Pexp * 3 + Qexp * 2 + Rexp * 2 + Sexp * 2},
  {7956, FIVEexp*8 + Qexp * 4 + Rexp * 2 + Sexp * 2},
  {-10096, FIVEexp*8 + Pexp * 4 + Rexp * 3 + Sexp * 2},
  {-512, FIVEexp*10 + Pexp + Qexp * 2 + Rexp * 3 + Sexp * 2},
  {672, FIVEexp*10 + Pexp * 2 + Rexp * 4 + Sexp * 2},
  {-756, FIVEexp*8 + Pexp * 6 + Qexp + Sexp * 3},
  {-37618, FIVEexp*7 + Pexp * 3 + Qexp * 3 + Sexp * 3},
  {-62064, FIVEexp*7 + Qexp * 5 + Sexp * 3},
  {48531, FIVEexp*7 + Pexp * 4 + Qexp + Rexp + Sexp * 3},
  {19344, FIVEexp*8 + Pexp + Qexp * 3 + Rexp + Sexp * 3},
  {-6124, FIVEexp*9 + Pexp * 2 + Qexp + Rexp * 2 + Sexp * 3},
  {-192, FIVEexp*11 + Qexp + Rexp * 3 + Sexp * 3},
  {21357, FIVEexp*7 + Pexp * 5 + Sexp * 4},
  {10469, FIVEexp*9 + Pexp * 2 + Qexp * 2 + Sexp * 4},
  {-11241, FIVEexp*9 + Pexp * 3 + Rexp + Sexp * 4},
  {1944, FIVEexp*10 + Qexp * 2 + Rexp + Sexp * 4},
  {276, FIVEexp*12 + Pexp + Rexp * 2 + Sexp * 4},
  {-99, FIVEexp*13 + Pexp + Qexp + Sexp * 5},
  {644, FIVEexp*12 + Sexp * 6},
  {0, END_ARRAY},
};

// Generate sixth degree polynomial F20 as indicated in
// Dummit's Solving Solvable Quintics article.
static void FactorPolynomialF20(void)
{
  int ctr;
  int* ptrValues;
  const struct stQuinticF20* pstQuinticF20;
  values[0] = 6;          // Degree of polynomial to factor.
  ptrValues = &values[1];

  intToBigInteger(&Rat1.numerator, 0);
  intToBigInteger(&Rat1.denominator, 1);
  intToBigInteger(&commonDenom, 1);             // Common denominator.
  intToBigInteger(&tmp0, 0);
  intToBigInteger(&tmp1, 0);
  intToBigInteger(&tmp2, 0);
  intToBigInteger(&tmp3, 0);
  intToBigInteger(&tmp4, 0);
  intToBigInteger(&tmp5, 0);
  pstQuinticF20 = astQuinticF20;
  for (;;)
  {
    if (pstQuinticF20->pqrs == END_COEFF)
    {
      BigIntGcd(&commonDenom, &Rat1.denominator, &tmp7);
      BigIntMultiply(&tmp1, &Rat1.denominator, &tmp0);
      BigIntDivide(&tmp0, &tmp7, &tmp0);
      BigIntMultiply(&tmp2, &Rat1.denominator, &tmp1);
      BigIntDivide(&tmp1, &tmp7, &tmp1);
      BigIntMultiply(&tmp3, &Rat1.denominator, &tmp2);
      BigIntDivide(&tmp2, &tmp7, &tmp2);
      BigIntMultiply(&tmp4, &Rat1.denominator, &tmp3);
      BigIntDivide(&tmp3, &tmp7, &tmp3);
      BigIntMultiply(&tmp5, &Rat1.denominator, &tmp4);
      BigIntDivide(&tmp4, &tmp7, &tmp4);
      BigIntMultiply(&tmp6, &Rat1.denominator, &tmp5);
      BigIntDivide(&tmp5, &tmp7, &tmp5);
      BigIntMultiply(&Rat1.numerator, &commonDenom, &tmp6);
      BigIntDivide(&tmp6, &tmp7, &tmp6);
      BigIntMultiply(&commonDenom, &Rat1.denominator, &commonDenom);
      BigIntDivide(&commonDenom, &tmp7, &commonDenom);
      intToBigInteger(&Rat1.numerator, 0);
      intToBigInteger(&Rat1.denominator, 1);
      if (pstQuinticF20->coeff < 0)
      {
        break;
      }
    }
    else
    {
      intToBigInteger(&Rat2.numerator, pstQuinticF20->coeff);
      intToBigInteger(&Rat2.denominator, 1);
      // Multiply by power of p.
      for (ctr = pstQuinticF20->pqrs & 0xF000; ctr > 0; ctr -= 0x1000)
      {
        BigRationalMultiply(&Rat2, &RatDeprCubic, &Rat2);
      }
      // Multiply by power of q.
      for (ctr = pstQuinticF20->pqrs & 0x0F00; ctr > 0; ctr -= 0x0100)
      {
        BigRationalMultiply(&Rat2, &RatDeprQuadratic, &Rat2);
      }
      // Multiply by power of r.
      for (ctr = pstQuinticF20->pqrs & 0x00F0; ctr > 0; ctr -= 0x0010)
      {
        BigRationalMultiply(&Rat2, &RatDeprLinear, &Rat2);
      }
      // Multiply by power of s.
      for (ctr = pstQuinticF20->pqrs & 0x000F; ctr > 0; ctr--)
      {
        BigRationalMultiply(&Rat2, &RatDeprIndependent, &Rat2);
      }
      BigRationalAdd(&Rat1, &Rat2, &Rat1);
    }
    pstQuinticF20++;
  }
  NumberLength = tmp0.nbrLimbs;
  BigInteger2IntArray(ptrValues, &tmp0);
  ptrValues += 1 + numLimbs(ptrValues);
  NumberLength = tmp1.nbrLimbs;
  BigInteger2IntArray(ptrValues, &tmp1);
  ptrValues += 1 + numLimbs(ptrValues);
  NumberLength = tmp2.nbrLimbs;
  BigInteger2IntArray(ptrValues, &tmp2);
  ptrValues += 1 + numLimbs(ptrValues);
  NumberLength = tmp3.nbrLimbs;
  BigInteger2IntArray(ptrValues, &tmp3);
  ptrValues += 1 + numLimbs(ptrValues);
  NumberLength = tmp4.nbrLimbs;
  BigInteger2IntArray(ptrValues, &tmp4);
  ptrValues += 1 + numLimbs(ptrValues);
  NumberLength = tmp5.nbrLimbs;
  BigInteger2IntArray(ptrValues, &tmp5);
  ptrValues += 1 + numLimbs(ptrValues);
  NumberLength = tmp6.nbrLimbs;
  BigInteger2IntArray(ptrValues, &tmp6);
  FactorPolyOverIntegers();
}

// Compute rational value from array of monomials.
static void computeFormula(struct monomial** ppstMonomial, BigRational *rat)
{
  struct monomial* pstMonomial = *ppstMonomial;
  intToBigInteger(&rat->numerator, 0);
  intToBigInteger(&rat->denominator, 1);
  while (pstMonomial->coefficient != 0)
  {
    int exponents = pstMonomial->exponents;
    intToBigInteger(&Rat1.numerator, 1);
    intToBigInteger(&Rat1.denominator, 1);
    for (int ctr = 4; ctr >= 0; ctr--)
    {
      BigRationalMultiply(&Rat1, &Rat1, &Rat1);
      if (exponents & (FIVEexp << ctr))
      {
        BigRationalMultiplyByInt(&Rat1, 5, &Rat1);
      }
      if (exponents & (Pexp << ctr))
      {
        BigRationalMultiply(&Rat1, &RatDeprCubic, &Rat1);
      }
      if (exponents & (Qexp << ctr))
      {
        BigRationalMultiply(&Rat1, &RatDeprQuadratic, &Rat1);
      }
      if (exponents & (Rexp << ctr))
      {
        BigRationalMultiply(&Rat1, &RatDeprLinear, &Rat1);
      }
      if (exponents & (Sexp << ctr))
      {
        BigRationalMultiply(&Rat1, &RatDeprIndependent, &Rat1);
      }
    }
    BigRationalMultiplyByInt(&Rat1, pstMonomial->coefficient, &Rat1);
    BigRationalAdd(rat, &Rat1, rat);
    pstMonomial++;
  }
  *ppstMonomial = pstMonomial + 1;  // Skip terminator.
}

// Compute T1, T2, T3, T4, l0, v and O as follows:
// T1 = (b10 + b11*t + b12*t^2 + b13*t^3 + b14*t^4 + b15*t^5)/(2F)
// T2 = (b20 + b21*t + b22*t^2 + b23*t^3 + b24*t^4 + b25*t^5)/(2DF)
// T3 = (b30 + b31*t + b32*t^2 + b33*t^3 + b34*t^4 + b35*t^5)/(2F)
// T4 = (b40 + b41*t + b42*t^2 + b43*t^3 + b44*t^4 + b45*t^5)/(2DF)
// l0 = (a0 + a1*t + a2*t^2 + a3*t^3 + a4*t^4 + a5*t^5)/F
// v = (c0 + c1*t + c2*t^2 + c3*t^3 + c4*t^4 + c5*t^5)/(2DF)
// O = (o0 + o1*t + o2*t^2 + o3*t^3 + o4*t^4 + o5*t^5)/(DF)
// where t is the rational root of polynomial F20.
static void computeArrayValues(void)
{
  struct BigRational* pRatValues = &RatValues[0];
  struct monomial* pstMonomial = arrayF;

  computeFormula(&pstMonomial, &Rat5);    // Compute F.
  pstMonomial = arrayB;
  for (int index = 0; index < sizeof(RatValues)/sizeof(RatValues[0]); index++)
  {
    computeFormula(&pstMonomial, pRatValues);
    intToBigInteger(&Rat3.numerator, 1);
    intToBigInteger(&Rat3.denominator, 1);
    for (int ctr = 1; ctr <= 5; ctr++)
    {
      BigRationalMultiply(&Rat3, &RatRoot, &Rat3);
      computeFormula(&pstMonomial, &Rat4);
      BigRationalMultiply(&Rat3, &Rat4, &Rat4);
      BigRationalAdd(pRatValues, &Rat4, pRatValues);
    }
    BigRationalDivide(pRatValues, &Rat5, pRatValues);           // Divide by F.
    ForceDenominatorPositive(pRatValues);
    pRatValues++;
  }
  BigRationalDivide(&RatValues[index_T2], &RatDiscr, &RatValues[index_T2]);   // Divide by discriminant.
  BigRationalDivide(&RatValues[index_T4], &RatDiscr, &RatValues[index_T4]);   // Divide by discriminant.
  BigRationalDivide(&RatValues[index_v], &RatDiscr, &RatValues[index_v]);     // Divide by discriminant.
  BigRationalDivide(&RatValues[index_O], &RatDiscr, &RatValues[index_O]);     // Divide by discriminant.
  BigRationalDivideByInt(&RatValues[index_T1], -2, &RatValues[index_T1]);     // Divide by 2.
  ForceDenominatorPositive(&RatValues[index_T1]);
  BigRationalDivideByInt(&RatValues[index_T2], 2, &RatValues[index_T2]);      // Divide by 2.
  BigRationalDivideByInt(&RatValues[index_T3], 2, &RatValues[index_T3]);      // Divide by 2.
  BigRationalDivideByInt(&RatValues[index_T4], 2, &RatValues[index_T4]);      // Divide by 2.
  BigRationalDivideByInt(&RatValues[index_v], 2, &RatValues[index_v]);        // Divide by 2.
}

void start5thRoot(void)
{
  if (pretty == PRETTY_PRINT)
  {
    showText("<span class=\"root\"><span class=\"befrad\" aria-hidden=\"true\">5</span><span class=\"radicand5\">");
  }
  else if (pretty == TEX)
  {
    showText("\\sqrt[5]{");
  }
  else
  {
    showText("(");
  }
}

void end5thRoot(void)
{
  if (pretty == PRETTY_PRINT)
  {
    showText("</span></span>");
  }
  else if (pretty == TEX)
  {
    showText("}");
  }
  else
  {
    showText(")^(1/5)");
  }
}

static void showSqrt5(void)
{
  if (pretty == PRETTY_PRINT)
  {
    startSqrt();
    *ptrOutput++ = '5';
    endSqrt();
  }
  else if (pretty == TEX)
  {
    showText("\\sqrt{5}");
  }
  else
  {
    showText("5^(1/2)");
  }
}

static void showMinusOnePlusMinusSqrt5(const char *sign)
{
  *ptrOutput++ = '(';
  showText(ptrMinus);
  *ptrOutput++ = '1';
  showText(sign);
  showSqrt5();
  *ptrOutput++ = ')';
}

static void showSqrtTenPlusMinusTwoTimesSqrt5(const char *sign)
{
  startSqrt();
  showText("10 ");
  showText(sign);
  *ptrOutput++ = ' ';
  showText(" 2");
  if (pretty == PARI_GP)
  {
    *ptrOutput++ = ' ';
  }
  showText(ptrTimes);
  if (pretty == PARI_GP)
  {
    *ptrOutput++ = ' ';
  }
  showSqrt5();
  endSqrt();
}

// Show (R +/- S*d)^(1/2). At this moment S is not zero.
static void showSqRoot1(enum eSign sign, BigRational *ptrRatR, BigRational *ptrRatS)
{
  enum eSign signBak = ptrRatS->numerator.sign;
  startSqrt();
  if (!BigIntIsZero(&ptrRatR->numerator))
  {
    ptrRatS->numerator.sign = SIGN_POSITIVE;
    showRational(ptrRatR);
    showPlusSignOn(sign == signBak, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
  }
  CopyBigInt(&Rat3.numerator, &RatDiscr.numerator);
  CopyBigInt(&Rat3.denominator, &RatDiscr.denominator);
  ShowRationalAndSqrParts(ptrRatS, &Rat3, 2, ptrTimes);
  endSqrt();
  ptrRatS->numerator.sign = signBak;
}

// Show ((-1 +/- sqrt(5))(M +/- N*d))^(1/2)
static void showSqRoot2(enum eSign sign)
{
  enum eSign signBak = RatS.numerator.sign;
  startSqrt();
  *ptrOutput++ = '(';
  showText(ptrMinus);
  *ptrOutput++ = '1';
  showPlusSignOn(sign == SIGN_POSITIVE, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
  showSqrt5();
  *ptrOutput++ = ')';
  showText(ptrTimes);
  *ptrOutput++ = ' ';
  if (!BigIntIsZero(&RatM.numerator))
  {
    RatS.numerator.sign = SIGN_POSITIVE;
    showRationalNoParen(&RatM);
    showPlusSignOn(sign == signBak, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
  }
  CopyBigInt(&Rat4.numerator, &RatDiscr.numerator);
  CopyBigInt(&Rat4.denominator, &RatDiscr.denominator);
  ShowRationalAndSqrParts(&RatN, &Rat4, 2, ptrTimes);
  endSqrt();
  RatS.numerator.sign = signBak;
}

static void ShowQuinticsRootsRealR(int multiplicity)
{
  BigRationalDivideByInt(&RatQuartic, -5, &RatQuartic);
  ForceDenominatorPositive(&RatQuartic);
  startLine();
  if (pretty == TEX)
  {
    showText("S_1 = ");
  }
  else
  {
    showText("<var>S</var><sub>1</sub> = ");
  }
  showMinusOnePlusMinusSqrt5("+");
  showText(ptrTimes);
  if (pretty == TEX)
  {
    showText("(R_1 + R_4)");
  }
  else
  {
    showText("(<var>R</var><sub>1</sub> + <var>R</var><sub>4</sub>) + ");
  }
  showMinusOnePlusMinusSqrt5(ptrMinus);
  showText(ptrTimes);
  if (pretty == TEX)
  {
    showText("(R_2 + R_3)");
    endLine();
    startLine();
    showText("S_2 = ");
  }
  else
  {
    showText("(<var>R</var><sub>2</sub> + <var>R</var><sub>3</sub>)");
    endLine();
    startLine();
    showText("<var>S</var><sub>2</sub> = ");
  }
  showMinusOnePlusMinusSqrt5("+");
  showText(ptrTimes);
  if (pretty == TEX)
  {
    showText("(R_2 + R_3) + ");
  }
  else
  {
    showText("(<var>R</var><sub>2</sub> + <var>R</var><sub>3</sub>) + ");
  }
  showMinusOnePlusMinusSqrt5(ptrMinus);
  showText(ptrTimes);
  if (pretty == TEX)
  {
    showText("(R_1 + R_4)");
    endLine();
    startLine();
    showText("T_1 = ");
  }
  else
  {
    showText("(<var>R</var><sub>1</sub> + <var>R</var><sub>4</sub>)");
    endLine();
    startLine();
    showText("<var>T</var><sub>1</sub> = ");
  }
  showSqrtTenPlusMinusTwoTimesSqrt5("+");
  showText(ptrTimes);
  if (pretty == TEX)
  {
    showText("(R_4 - R_1) +");
  }
  else
  {
    showText("(<var>R</var><sub>4</sub> ");
    showText(ptrMinus);
    showText(" <var>R</var><sub>1</sub>) +");
  }
  if (pretty == PARI_GP)
  {
    *ptrOutput++ = ' ';
  }
  showSqrtTenPlusMinusTwoTimesSqrt5(ptrMinus);
  showText(ptrTimes);
  if (pretty == TEX)
  {
    showText("(R_3 - R2)");
    endLine();
    startLine();
    showText("T_2 = ");
  }
  else
  {
    showText("(<var>R</var><sub>3</sub> ");
    showText(ptrMinus);
    showText(" <var>R</var><sub>2</sub>)");
    endLine();
    startLine();
    showText("<var>T</var><sub>2</sub> = ");
  }
  showSqrtTenPlusMinusTwoTimesSqrt5("+");
  showText(ptrTimes);
  if (pretty == TEX)
  {
    showText("(R_3 - R_2) +");
  }
  else
  {
    showText("(<var>R</var><sub>3</sub> ");
    showText(ptrMinus);
    showText(" <var>R</var><sub>2</sub>) +");
  }
  if (pretty == PARI_GP)
  {
    *ptrOutput++ = ' ';
  }
  showSqrtTenPlusMinusTwoTimesSqrt5(ptrMinus);
  showText(ptrTimes);
  if (pretty == TEX)
  {
    showText("(R_4 - R_1)");
  }
  else
  {
    showText("(<var>R</var><sub>4</sub> ");
    showText(ptrMinus);
    showText(" <var>R</var><sub>1</sub>)");
  }
  endLine();
  for (int ctr = 1; ctr <= 5; ctr++)
  {
    showX(multiplicity);
    if (!BigIntIsZero(&RatQuartic.numerator))
    {
      showRationalNoParen(&RatQuartic);
      showText(" + ");
    }
    switch (ctr)
    {
    case 1:
      if (pretty == PRETTY_PRINT)
      {
        showRatConstants("<var>R</var><sub>1</sub> + <var>R</var><sub>2</sub> + <var>R</var><sub>3</sub> + <var>R</var><sub>4</sub>", "5");
      }
      else if (pretty == TEX)
      {
        showRatConstants("R_1 + R_2 + R_3 + R_4", "5");
      }
      else
      {
        showText("(<var>R</var><sub>1</sub> + <var>R</var><sub>2</sub> + <var>R</var><sub>3</sub> + <var>R</var><sub>4</sub>) / 5");
      }
      break;
    case 2:
      if (pretty == PRETTY_PRINT)
      {
        showRatConstants("<var>S</var><sub>1</sub>", "20");
        showText("+ i ");
        showText(ptrTimes);
        showRatConstants("<var>T</var><sub>1</sub>", "20");
      }
      else if (pretty == TEX)
      {
        showRatConstants("S_1", "20");
        showText("+ i ");
        showRatConstants("T_1", "20");
      }
      else
      {
        showText("(<var>S</var><sub>1</sub> + I * <var>T</var><sub>1</sub>) / 20");
      }
      break;
    case 3:
      if (pretty == PRETTY_PRINT)
      {
        showRatConstants("<var>S</var><sub>1</sub>", "20");
        showText(ptrMinus);
        showText(" i ");
        showText(ptrTimes);
        showRatConstants("<var>T</var><sub>1</sub>", "20");
      }
      else if (pretty == TEX)
      {
        showRatConstants("S_1", "20");
        showText("- i ");
        showRatConstants("T_1", "20");
      }
      else
      {
        showText("(<var>S</var><sub>1</sub> - I * <var>T</var><sub>1</sub>) / 20");
      }
      break;
    case 4:
      if (pretty == PRETTY_PRINT)
      {
        showRatConstants("<var>S</var><sub>2</sub>", "20");
        showText("+ i ");
        showText(ptrTimes);
        showRatConstants("<var>T</var><sub>2</sub>", "20");
      }
      else if (pretty == TEX)
      {
        showRatConstants("S_2", "20");
        showText("+ i ");
        showRatConstants("T_2", "20");
      }
      else
      {
        showText("(<var>S</var><sub>2</sub> + I * <var>T</var><sub>2</sub>) / 20");
      }
      break;
    case 5:
      if (pretty == PRETTY_PRINT)
      {
        showRatConstants("<var>S</var><sub>2</sub>", "20");
        showText(ptrMinus);
        showText(" i ");
        showText(ptrTimes);
        showRatConstants("<var>T</var><sub>2</sub>", "20");
      }
      else if (pretty == TEX)
      {
        showRatConstants("S_2", "20");
        showText("- i ");
        showText(ptrTimes);
        showRatConstants("T_2", "20");
      }
      else
      {
        showText("(<var>S</var><sub>2</sub> - I * <var>T</var><sub>2</sub>) / 20");
      }
      break;
    default:
      break;
    }
    endLine();
  }
}

static void NumberIsNotRational(enum eSign sign)
{
  if (firstNumberShown)
  {
    showPlusSignOn(sign == SIGN_POSITIVE, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
  }
  else
  {
    if (BigIntIsZero(&Rat3.numerator))
    {
      start5thRoot();
      if (sign == SIGN_NEGATIVE)
      {
        showText(ptrMinus);
        *ptrOutput++ = ' ';
      }
    }
    else
    {
      start5thRoot();
      showRationalNoParen(&Rat3);
      showPlusSignOn(sign == SIGN_POSITIVE, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
    }
    firstNumberShown = true;
  }
}

static void showRn(int groupOrder)
{
  // If group order is 20:
  // R_i = l0 - T1 +/- 5*T2 +/- sqrt(r+s*d) +/- sqrt(r-s*d)
  // Otherwise:
  // R_i = l0 - T1 +/- sqrt(5)*T2*d +/- sqrt(r+s*d) +/- sqrt(r-s*d)
  // First  +/-: plus for R_2 and R_3, minus for R_1 and R_4.
  // Second +/-: if O > 0: plus for R_3 and R_4, minus for R_1 and R_2.
  //             if O < 0: plus for R_2 and R_4, minus for R_1 and R_3.
  // Third  +/-: if O > 0: plus for R_1 and R_3, minus for R_2 and R_4.
  //             if O < 0: plus for R_1 and R_2, minus for R_3 and R_4.
  for (int ctr = 1; ctr <= 4; ctr++)
  {
    firstNumberShown = false;
    BigRational* ptrRatR = (groupOrder == 10 || ctr == 1 || ctr == 4 ? &RatR : &RatR2);
    BigRational* ptrRatS = (groupOrder == 10 || ctr == 1 || ctr == 4 ? &RatS : &RatS2);
    enum eSign firstSign, secondSign;
    startLine();
    if (pretty == TEX)
    {
      showText("R_");
      *ptrOutput++ = (char)(ctr + '0');
      showText(" = ");
    }
    else
    {
      showText("<var>R</var><sub>");
      *ptrOutput++ = (char)(ctr + '0');
      showText("</sub> = ");
    }
    BigRationalDivideByInt(&RatValues[index_T1], 2, &Rat3);
    BigRationalAdd(&RatValues[index_l0], &Rat3, &Rat3);
    if (groupOrder == 20)
    {         // Group order is 20.
      BigRationalMultiplyByInt(&RatValues[index_T2], 5, &Rat2);
      BigRationalDivideByInt(&Rat2, 2, &Rat2);
      if (ctr == 2 || ctr == 3)
      {
        BigRationalAdd(&Rat3, &Rat2, &Rat3);
      }
      else
      {
        BigRationalSubt(&Rat3, &Rat2, &Rat3);
      }
      ForceDenominatorPositive(&Rat3);
    }
    else
    {        // Group order is 10.
      BigRationalMultiply(&RatValues[index_T2], &Rat4, &Rat2);
      BigRationalDivideByInt(&Rat2, 2, &Rat2);
      ForceDenominatorPositive(&Rat3);
      ForceDenominatorPositive(&Rat2);
      if (!BigIntIsZero(&Rat2.numerator))
      {
        firstNumberShown = true;
        start5thRoot();
        if (!BigIntIsZero(&Rat1.numerator))
        {
          showRational(&Rat3);
        }
        showPlusSignOn(ctr == 2 || ctr == 3, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
        showRational(&Rat2);
        if (pretty == PARI_GP)
        {
          *ptrOutput++ = ' ';
        }
        showText(ptrTimes);
        *ptrOutput++ = ' ';
        showSqrt5();
      }
    }
    if (RatValues[index_O].numerator.sign == SIGN_POSITIVE)
    {   // O > 0.
      firstSign = ((ctr == 3 || ctr == 4) ? SIGN_POSITIVE : SIGN_NEGATIVE);
    }
    else
    {   // O < 0.
      firstSign = ((ctr == 2 || ctr == 4) ? SIGN_POSITIVE : SIGN_NEGATIVE);
    }
    if (RatValues[index_O].numerator.sign == SIGN_POSITIVE)
    {   // O > 0.
      secondSign = ((ctr == 1 || ctr == 3) ? SIGN_POSITIVE : SIGN_NEGATIVE);
    }
    else
    {   // O < 0.
      secondSign = ((ctr == 1 || ctr == 2) ? SIGN_POSITIVE : SIGN_NEGATIVE);
    }
    if (BigIntIsZero(&ptrRatS->numerator))
    {
      if (firstSign == secondSign)
      {
        NumberIsNotRational(firstSign);
        BigRationalMultiplyByInt(ptrRatR, 4, &Rat2);
        showSquareRootOfRational(&Rat2, 2, ptrTimes);
      }
      else
      {
        showText("0");
      }
    }
    else
    {
      // Test whether the sum of two square roots can be expressed with
      // only one square root.
      // (p+q)^2 = p^2 + 2pq + q^2
      // (p^(1/2) + q^(1/2))^2 = (p+q) + 2p^(1/2)*q^(1/2)
      // p^(1/2) + q^(1/2) = ((p+q) + 2(pq)^(1/2))^(1/2)
      // p^(1/2) - q^(1/2) = ((p+q) - 2(pq)^(1/2))^(1/2)
      // (r+s*d)^(1/2) + (r-s*d)^(1/2) = (2*r + 2(r^2 - 5s^2)^(1/2))^(1/2)
      // (r+s*d)^(1/2) - (r-s*d)^(1/2) = (2*r - 2(r^2 - 5s^2)^(1/2))^(1/2)
      // So we have to check whether (r^2 - 5s^2) is a perfect square or not.
      BigRationalMultiply(ptrRatR, ptrRatR, &Rat1);
      BigRationalMultiply(ptrRatS, ptrRatS, &Rat2);
      BigRationalMultiplyByInt(&Rat2, 5, &Rat2);
      BigRationalSubt(&Rat1, &Rat2, &Rat1);
      ForceDenominatorPositive(&Rat1);
      intToBigInteger(&tmp5, 1);    // Indicate not perfect square.
      if (Rat1.numerator.sign == SIGN_POSITIVE)
      {
        squareRoot(Rat1.numerator.limbs, tmp4.limbs, Rat1.numerator.nbrLimbs, &tmp4.nbrLimbs);
        tmp4.sign = SIGN_POSITIVE;
        BigIntMultiply(&tmp4, &tmp4, &tmp5);
        BigIntSubt(&tmp5, &Rat1.numerator, &tmp5);
        if (BigIntIsZero(&tmp5))
        {   // Numerator is perfect square.
          squareRoot(Rat1.denominator.limbs, tmp1.limbs, Rat1.denominator.nbrLimbs, &tmp1.nbrLimbs);
          tmp1.sign = SIGN_POSITIVE;
          BigIntMultiply(&tmp1, &tmp1, &tmp5);
          BigIntSubt(&tmp5, &Rat1.denominator, &tmp5);
        }
      }
      if (BigIntIsZero(&tmp5))
      {   // Number is perfect square.
        CopyBigInt(&Rat2.numerator, &tmp4);
        CopyBigInt(&Rat2.denominator, &tmp1);
          // At this moment Rat2 = (r^2 - 5s^2)^(1/2)
        if (firstSign == secondSign)
        {
          BigRationalAdd(ptrRatR, &Rat2, &Rat1);
        }
        else
        {
          BigRationalSubt(ptrRatR, &Rat2, &Rat1);
        }
        BigRationalMultiplyByInt(&Rat1, 2, &Rat1);
          // At this moment Rat1 = 2*r - 2(r^2 - 5s^2)^(1/2)
          // If this number is a perfect square, it can be reduced.
        squareRoot(Rat1.numerator.limbs, tmp4.limbs, Rat1.numerator.nbrLimbs, &tmp4.nbrLimbs);
        tmp4.sign = SIGN_POSITIVE;
        BigIntMultiply(&tmp4, &tmp4, &tmp5);
        BigIntSubt(&tmp5, &Rat1.numerator, &tmp5);
        if (BigIntIsZero(&tmp5))
        {   // Numerator is perfect square.
          squareRoot(Rat1.denominator.limbs, tmp1.limbs, Rat1.denominator.nbrLimbs, &tmp1.nbrLimbs);
          tmp1.sign = SIGN_POSITIVE;
          BigIntMultiply(&tmp1, &tmp1, &tmp5);
          BigIntSubt(&tmp5, &Rat1.denominator, &tmp5);
        }
        if (BigIntIsZero(&tmp5))
        {   // Number is perfect square.
          CopyBigInt(&Rat2.numerator, &tmp4);
          CopyBigInt(&Rat2.denominator, &tmp1);
          if (firstSign)
          {
            BigRationalAdd(&Rat3, &Rat2, &Rat2);
          }
          else
          {
            BigRationalSubt(&Rat3, &Rat2, &Rat2);
          }
          if (BigIntIsZero(&Rat2.numerator))
          {
            showText("0");
          }
          else
          {
            start5thRoot();
            firstNumberShown = true;
            showRational(&Rat2);
          }
        }
        else
        {
          NumberIsNotRational(firstSign);
          showSquareRootOfRational(&Rat1, 2, ptrTimes);
        }
      }
      else
      {
        NumberIsNotRational(firstSign);
        showSqRoot1(SIGN_POSITIVE, ptrRatR, ptrRatS);
        showPlusSignOn(secondSign == SIGN_POSITIVE, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
        showSqRoot1(SIGN_NEGATIVE, ptrRatR, ptrRatS);
      }
    }
    if (firstNumberShown)
    {
      end5thRoot();
    }
    endLine();
  }
}

static void GaloisGroupHasOrder20(int multiplicity)
{
  (void)multiplicity;
  CopyBigInt(&Rat1.numerator, &RatDiscr.numerator);
  CopyBigInt(&Rat1.denominator, &RatDiscr.denominator);
  MultiplyRationalBySqrtRational(&RatN, &Rat1);
  MultiplyRationalBySqrtRational(&RatValues[index_T2], &RatDiscr);
  if (RatDiscr.denominator.nbrLimbs == 1 && RatDiscr.denominator.limbs[0].x == 1 &&
    getRemainder(&RatDiscr.numerator, 5) == 0 && getRemainder(&RatDiscr.numerator, 25) != 0)
  {      // The discriminant has the form 5*n^2.
         // Compute r = 5(m+n)/8, s = (m+5n)/8.
    BigRationalAdd(&RatM, &RatN, &RatR);         // R <- m + n
    BigRationalMultiplyByInt(&RatR, 5, &RatR);   // R <- 5(m + n)
    BigRationalDivideByInt(&RatR, 8, &RatR);     // R <- 5(m + n)/8
    BigRationalMultiplyByInt(&RatN, 5, &RatS);   // S <- 5n
    BigRationalAdd(&RatS, &RatM, &RatS);         // S <- m + 5n
    BigRationalDivideByInt(&RatS, 8, &RatS);     // S <- 5(m + 5n)/8
    BigIntChSign(&RatR.numerator);
    BigIntChSign(&RatS.numerator);
    // Compute r2 = 5(m-n)/8, s2 = (m-5n)/8.
    BigRationalSubt(&RatM, &RatN, &RatR2);       // R2 <- m - n
    BigRationalMultiplyByInt(&RatR2, 5, &RatR2); // R2 <- 5(m - n)
    BigRationalDivideByInt(&RatR2, 8, &RatR2);   // R2 <- 5(m - n)/8
    BigRationalMultiplyByInt(&RatN, 5, &RatS2);  // S2 <- 5n
    BigRationalSubt(&RatS2, &RatM, &RatS2);      // S2 <- m - 5n
    BigRationalDivideByInt(&RatS2, 8, &RatS2);   // S2 <- 5(m - 5n)/8
    BigIntChSign(&RatR2.numerator);
    BigIntChSign(&RatS2.numerator);
    showRn(20);
  }
  else
  {
    // R_i = l0 - T1 +/- T2*sqrt(5*d) +/- (1/4)*(sqrt(5)+1)*sqrt(M+N*d)
    //                                +/- (1/4)*(sqrt(5)-1)*sqrt(M-N*d)
    // First  +/-: plus for R_2 and R_3, minus for R_1 and R_4.
    // Second +/-: if O > 0: plus for R_3 and R_4, minus for R_1 and R_2.
    //             if O < 0: plus for R_2 and R_4, minus for R_1 and R_3.
    // Third  +/-: if O > 0: plus for R_1 and R_3, minus for R_2 and R_4.
    //             if O < 0: plus for R_1 and R_2, minus for R_3 and R_4.

    for (int ctr = 1; ctr <= 4; ctr++)
    {
      startLine();
      if (pretty == TEX)
      {
        showText("R_");
        *ptrOutput++ = (char)(ctr + '0');
        showText(" = (");
      }
      else
      {
        showText("<var>R</var><sub>");
        *ptrOutput++ = (char)(ctr + '0');
        showText("</sub></var> = (");
      }
      BigRationalDivideByInt(&RatValues[index_T1], 2, &Rat2);
      BigRationalAdd(&RatValues[index_l0], &Rat2, &Rat1);
      ForceDenominatorPositive(&Rat1);
      if (!BigIntIsZero(&Rat1.numerator))
      {
        showRational(&Rat1);
      }
      CopyBigInt(&Rat1.numerator, &RatValues[index_T2].numerator);
      CopyBigInt(&Rat1.denominator, &RatValues[index_T2].denominator);
      showPlusSignOn((ctr == 2 || ctr == 3) == (Rat1.numerator.sign == SIGN_POSITIVE), TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
      Rat1.numerator.sign = SIGN_POSITIVE;
      BigRationalDivideByInt(&Rat1, 2, &Rat2);
      BigRationalMultiplyByInt(&RatDiscr, 5, &Rat2);
      ShowRationalAndSqrParts(&Rat1, &Rat2, 2, ptrTimes);
      if (RatValues[index_O].numerator.sign == SIGN_POSITIVE)
      {   // O > 0.
        showPlusSignOn(ctr == 3 || ctr == 4, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
      }
      else
      {   // O < 0.
        showPlusSignOn(ctr == 2 || ctr == 4, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
      }
      showSqRoot2(SIGN_POSITIVE);
      if (RatValues[index_O].numerator.sign == SIGN_POSITIVE)
      {   // O > 0.
        showPlusSignOn(ctr == 1 || ctr == 3, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
      }
      else
      {   // O < 0.
        showPlusSignOn(ctr == 1 || ctr == 2, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
      }
      showSqRoot2(SIGN_NEGATIVE);
      showText(")^(1/5)");
    }
  }
}

static void GaloisGroupHasOrder10(int multiplicity)
{
  (void)multiplicity;
  // Set Y <- M + N*d. This value is already in Rat1.
  // Compute r = 5m/4, s = m/4.
  BigRationalDivideByInt(&Rat1, -4, &RatS);
  BigRationalMultiplyByInt(&RatS, 5, &RatR);
  RatR.numerator.sign = SIGN_POSITIVE;
  RatR.denominator.sign = SIGN_POSITIVE;
  intToBigInteger(&RatDiscr.numerator, 5);
  intToBigInteger(&RatDiscr.denominator, 1);
  showRn(10);    // Show all four values of R.
}

static void GaloisGroupHasOrder5(int multiplicity)
{
  (void)multiplicity;
  // Compute the roots l1, l2, l3 and l4 of the pair of quadratic equations.
  // l1, l4 = (-(T1 + T2*d) +/- sqrt(M + N*d)) / 2
  // l2, l3 = (-(T1 - T2*d) +/- sqrt(M - N*d)) / 2
  // The order of the roots must comply with: (l1 − l4)(l2 − l3) = O*d
  // If O > 0 -> l1 and l2: plus sign, l3 and l4: minus sign
  // If O < 0 -> l1 and l3: plus sign, l2 and l4: minus sign
  // sqrt(M + N*d) is already computed as Rat3 and d is Rat4.
  // Compute sqrt(M - N*d) as Rat5.
  BigRationalMultiply(&Rat4, &RatN, &Rat1);
  BigRationalSubt(&RatM, &Rat1, &Rat1);
  ForceDenominatorPositive(&Rat1);         // M - N*d
  squareRoot(Rat1.numerator.limbs, tmp4.limbs, Rat1.numerator.nbrLimbs, &tmp4.nbrLimbs);
  squareRoot(Rat1.denominator.limbs, tmp1.limbs, Rat1.denominator.nbrLimbs, &tmp1.nbrLimbs);
  CopyBigInt(&Rat5.numerator, &tmp4);      // Rat5 = square root of discriminant
  CopyBigInt(&Rat5.denominator, &tmp1);    // of quadratic equation.
  BigRationalDivideByInt(&Rat5, 2, &Rat5);

  BigRationalDivideByInt(&Rat3, 2, &Rat3);
  BigRationalMultiply(&RatValues[index_T2], &Rat4, &Rat2);    // T2*d
  BigRationalAdd(&RatValues[index_T1], &Rat2, &Rat1);         // T1 + T2*d
  BigRationalDivideByInt(&Rat1, -2, &Rat1);
  ForceDenominatorPositive(&Rat1);                            // -(T1 + T2*d)/2
  BigRationalSubt(&RatValues[index_T1], &Rat2, &Rat2);        // T1 - T2*d
  BigRationalDivideByInt(&Rat2, -2, &Rat2);
  ForceDenominatorPositive(&Rat2);                            // -(T1 - T2*d)/2
  BigRationalAdd(&Rat1, &Rat3, &RatValues[index_T1]);         // l1
  BigRationalSubt(&Rat1, &Rat3, &RatValues[index_T4]);        // l4
  if (RatValues[index_O].numerator.sign == SIGN_POSITIVE)
  {   // O > 0.
    BigRationalAdd(&Rat2, &Rat5, &RatValues[index_T2]);       // l2
    BigRationalSubt(&Rat2, &Rat5, &RatValues[index_T3]);      // l3
  }
  else
  {   // O < 0.
    BigRationalAdd(&Rat2, &Rat5, &RatValues[index_T3]);       // l3
    BigRationalSubt(&Rat2, &Rat5, &RatValues[index_T2]);      // l2
  }
  // l = l0 - (l1 + l2 + l3 + l4)/4
  // m = (l1 - l2 - l3 + l4) / 4
  // U = sqrt(10+2*sqrt(5)) / 4
  // V = sqrt(10-2*sqrt(5)) / 4
  // K = (l1 - l4)*U + (l2 - l3)*V
  // R1 = l + m * sqrt(5) + i((l1 - l4)*U + (l2 - l3)*V)
  // R2 = l - m * sqrt(5) + i((l3 - l2)*U + (l1 - l4)*V)
  // R3 = l - m * sqrt(5) + i((l2 - l3)*U + (l4 - l1)*V)
  // R4 = l + m * sqrt(5) + i((l4 - l1)*U + (l3 - l2)*V)
  BigRationalAdd(&RatValues[index_T1], &RatValues[index_T4], &Rat2);
  BigRationalAdd(&RatValues[index_T2], &RatValues[index_T3], &Rat3);
  BigRationalAdd(&Rat2, &Rat3, &Rat1);     // l1 + l2 + l3 + l4
  BigRationalSubt(&Rat2, &Rat3, &Rat2);    // l1 - l2 - l3 + l4
  BigRationalDivideByInt(&Rat1, 4, &Rat1); // (l1 + l2 + l3 + l4)/4
  BigRationalDivideByInt(&Rat2, 4, &Rat2); // (l1 - l2 - l3 + l4)/4
  BigRationalSubt(&RatValues[index_l0], &Rat1, &Rat1);  // l
  enum eSign signRat2 = Rat2.numerator.sign;
  Rat2.numerator.sign = SIGN_POSITIVE;
  for (int ctr = 1; ctr <= 4; ctr++)
  {
    startLine();
    if (pretty == TEX)
    {
      showText("R_");
      *ptrOutput++ = (char)(ctr + '0');
      showText(" = ");
    }
    else
    {
      showText("<var>R</var><sub>");
      *ptrOutput++ = (char)(ctr + '0');
      showText("</sub> = ");
    }
    start5thRoot();
    showRationalNoParen(&Rat1);
    showPlusSignOn((signRat2 == SIGN_POSITIVE) == (ctr == 1 || ctr == 4), TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
    showRational(&Rat2);
    if (pretty == PARI_GP)
    {
      *ptrOutput++ = ' ';
    }
    showText(ptrTimes);
    if (pretty == PARI_GP)
    {
      *ptrOutput++ = ' ';
    }
    showSqrt5();
    *ptrOutput++ = ' ';
    if (ctr >= 3)
    {
      showText(ptrMinus);
    }
    else
    {
      *ptrOutput++ = '+';
    }
    showText(pretty == PARI_GP? " I ": " i ");
    showText(ptrTimes);
    showText("&nbsp;");
    startParen();
    BigRationalSubt(&RatValues[index_T1], &RatValues[index_T4], &Rat3);  // l1 - l4
    BigRationalDivideByInt(&Rat3, 4, &Rat3);
    ForceDenominatorPositive(&Rat3);
    showRational(&Rat3);
    if (pretty == PARI_GP)
    {
      *ptrOutput++ = ' ';
    }
    showText(ptrTimes);
    if ((ctr == 1) || (ctr == 3))
    {
      showSqrtTenPlusMinusTwoTimesSqrt5("+");
    }
    else
    {
      showSqrtTenPlusMinusTwoTimesSqrt5(ptrMinus);
    }
    BigRationalSubt(&RatValues[index_T3], &RatValues[index_T2], &Rat3);  // l3 - l2
    if ((ctr == 2) || (ctr == 3))
    {
      BigRationalDivideByInt(&Rat3, 4, &Rat3);
    }
    else
    {
      BigRationalDivideByInt(&Rat3, -4, &Rat3);
    }
    ForceDenominatorPositive(&Rat3);
    showPlusSignOn(Rat3.numerator.sign == SIGN_POSITIVE, TYPE_PM_SPACE_BEFORE | TYPE_PM_SPACE_AFTER);
    Rat3.numerator.sign = SIGN_POSITIVE;
    showRational(&Rat3);
    *ptrOutput++ = ' ';
    showText(ptrTimes);
    if ((ctr == 2) || (ctr == 4))
    {
      showSqrtTenPlusMinusTwoTimesSqrt5("+");
    }
    else
    {
      showSqrtTenPlusMinusTwoTimesSqrt5(ptrMinus);
    }
    endParen();
    end5thRoot();
    endLine();
  }
}

void QuinticEquation(const int* ptrPolynomial, int multiplicity)
{
  struct monomial* pstMonomial;
  const int* ptr;
  UncompressBigIntegerB(ptrPolynomial, &Independent);
  ptrPolynomial += 1 + numLimbs(ptrPolynomial);
  UncompressBigIntegerB(ptrPolynomial, &Linear);
  ptrPolynomial += 1 + numLimbs(ptrPolynomial);
  UncompressBigIntegerB(ptrPolynomial, &Quadratic);
  ptrPolynomial += 1 + numLimbs(ptrPolynomial);
  UncompressBigIntegerB(ptrPolynomial, &Cubic);
  ptrPolynomial += 1 + numLimbs(ptrPolynomial);
  UncompressBigIntegerB(ptrPolynomial, &Quartic);
  ptrPolynomial += 1 + numLimbs(ptrPolynomial);
  UncompressBigIntegerB(ptrPolynomial, &Quintic);
  // Get rational coefficients of monic equation.
  CopyBigInt(&RatQuartic.numerator, &Quartic);
  CopyBigInt(&RatQuartic.denominator, &Quintic);
  CopyBigInt(&RatCubic.numerator, &Cubic);
  CopyBigInt(&RatCubic.denominator, &Quintic);
  CopyBigInt(&RatQuadratic.numerator, &Quadratic);
  CopyBigInt(&RatQuadratic.denominator, &Quintic);
  CopyBigInt(&RatLinear.numerator, &Linear);
  CopyBigInt(&RatLinear.denominator, &Quintic);
  CopyBigInt(&RatIndependent.numerator, &Independent);
  CopyBigInt(&RatIndependent.denominator, &Quintic);
  // Compute coefficients of depressed equation x^5 + px^3 + qx^2 + rx + s
  // where: p = (5c - 2b^2)/5, q = (25d - 15bc + 4b^3)/25, r = (125e - 50bd + 15b^2*c - 3b^4)/125,
  // s = (3125f - 625be + 125b^2*d - 25b^3*c + 4b^5)/3125
  BigRationalMultiply(&RatQuartic, &RatCubic, &Rat2);
  BigRationalMultiplyByInt(&Rat2, -15, &RatDeprQuadratic);         // -15bc
  BigRationalMultiply(&RatQuartic, &RatQuadratic, &Rat2);
  BigRationalMultiplyByInt(&Rat2, -50, &RatDeprLinear);            // -50bd
  BigRationalMultiply(&RatQuartic, &RatLinear, &Rat2);
  BigRationalMultiplyByInt(&Rat2, -625, &RatDeprIndependent);      // -625be
  BigRationalMultiply(&RatQuartic, &RatQuartic, &Rat1);            // b^2
  BigRationalMultiplyByInt(&Rat1, -2, &RatDeprCubic);              // -2b^2
  BigRationalMultiply(&Rat1, &RatCubic, &Rat2);
  BigRationalMultiplyByInt(&Rat2, 15, &Rat2);                      // 15b^2*c
  BigRationalAdd(&RatDeprLinear, &Rat2, &RatDeprLinear);
  BigRationalMultiply(&Rat1, &RatQuadratic, &Rat2);
  BigRationalMultiplyByInt(&Rat2, 125, &Rat2);                     // 125b^2*d
  BigRationalAdd(&RatDeprIndependent, &Rat2, &RatDeprIndependent);
  BigRationalMultiply(&Rat1, &RatQuartic, &Rat1);                  // b^3
  BigRationalMultiplyByInt(&Rat1, 4, &Rat2);                       // 4b^3
  BigRationalAdd(&RatDeprQuadratic, &Rat2, &RatDeprQuadratic);
  BigRationalMultiply(&Rat1, &RatCubic, &Rat2);
  BigRationalMultiplyByInt(&Rat2, -25, &Rat2);                     // -25b^3*c
  BigRationalAdd(&RatDeprIndependent, &Rat2, &RatDeprIndependent);
  BigRationalMultiply(&Rat1, &RatQuartic, &Rat1);                  // b^4
  BigRationalMultiplyByInt(&Rat1, -3, &Rat2);                      // -3b^4
  BigRationalAdd(&RatDeprLinear, &Rat2, &RatDeprLinear);
  BigRationalMultiply(&Rat1, &RatQuartic, &Rat1);                  // b^5
  BigRationalMultiplyByInt(&Rat1, 4, &Rat2);                       // 4b^5
  BigRationalAdd(&RatDeprIndependent, &Rat2, &RatDeprIndependent);
  BigRationalDivideByInt(&RatDeprCubic, 5, &RatDeprCubic);
  BigRationalAdd(&RatDeprCubic, &RatCubic, &RatDeprCubic);
  BigRationalDivideByInt(&RatDeprQuadratic, 25, &RatDeprQuadratic);
  BigRationalAdd(&RatDeprQuadratic, &RatQuadratic, &RatDeprQuadratic);
  BigRationalDivideByInt(&RatDeprLinear, 125, &RatDeprLinear);
  BigRationalAdd(&RatDeprLinear, &RatLinear, &RatDeprLinear);
  BigRationalDivideByInt(&RatDeprIndependent, 3125, &RatDeprIndependent);
  BigRationalAdd(&RatDeprIndependent, &RatIndependent, &RatDeprIndependent);

  FactorPolynomialF20();
  // If there is a linear factor, that means that the quintic can be expressed with radicands.
  if (factorInfoInteger[0].degree != 1)
  {
    showText(lang ? "<p>La ecuación quíntica no se puede expresar mediante radicandos.</p>" :
      "The quintic equation cannot be expressed with radicands.");
    return;
  }
  ptr = factorInfoInteger[0].ptrPolyLifted;
  UncompressBigIntegerB(ptr, &RatRoot.numerator);          // Independent term.
  ptr += 1 + numLimbs(ptr);
  UncompressBigIntegerB(ptr, &RatRoot.denominator);        // Linear coefficient.
  BigIntChSign(&RatRoot.numerator);                        // root = -independent / linear.
  // Compute discriminant
  pstMonomial = arrayD;
  computeFormula(&pstMonomial, &RatDiscr);
  computeArrayValues();
  // Test whether the discriminant is a perfect square.
  squareRoot(RatDiscr.numerator.limbs, tmp4.limbs, RatDiscr.numerator.nbrLimbs, &tmp4.nbrLimbs);
  tmp4.sign = SIGN_POSITIVE;
  BigIntMultiply(&tmp4, &tmp4, &tmp5);
  BigIntSubt(&tmp5, &RatDiscr.numerator, &tmp5);
  if (BigIntIsZero(&tmp5))
  {   // Numerator of discriminant is perfect square.
    squareRoot(RatDiscr.denominator.limbs, tmp1.limbs, RatDiscr.denominator.nbrLimbs, &tmp1.nbrLimbs);
    tmp1.sign = SIGN_POSITIVE;
    BigIntMultiply(&tmp1, &tmp1, &tmp5);
    BigIntSubt(&tmp5, &RatDiscr.denominator, &tmp5);
  }
  // Let M = T1^2 + T2^2*D - 4*T3, N = 2*T1*T2 - 4*T4
  // and d = sqrt(D), then:
  // l1, l4 = (-(T1 + T2*d) +/- sqrt(M + N*d)) / 2
  // l2, l3 = (-(T1 - T2*d) +/- sqrt(M - N*d)) / 2
  // The order of the roots must comply with: (l1 − l4)(l2 − l3) = O*d
  // l1 - l4 = +/- sqrt(M + N*d)
  // l2 - l3 = +/- sqrt(M - N*d)
  // Compute M.
  BigRationalMultiply(&RatValues[index_T1], &RatValues[index_T1], &RatM);
  BigRationalMultiply(&RatValues[index_T2], &RatValues[index_T2], &Rat1);
  BigRationalMultiply(&Rat1, &RatDiscr, &Rat1);
  BigRationalAdd(&RatM, &Rat1, &RatM);
  BigRationalMultiplyByInt(&RatValues[index_T3], 4, &Rat1);
  BigRationalSubt(&RatM, &Rat1, &RatM);
  // Compute N.
  BigRationalMultiply(&RatValues[index_T1], &RatValues[index_T2], &RatN);
  BigRationalMultiplyByInt(&RatN, 2, &RatN);
  BigRationalMultiplyByInt(&RatValues[index_T4], 4, &Rat1);
  BigRationalSubt(&RatN, &Rat1, &RatN);
  // The real part of the primitive fifth root of 1 is Re5 = (-1+sqrt(5))/4
  // The real part of its square is: Re5_2 = (-1-sqrt(5))/4
  if (!BigIntIsZero(&tmp5))
  {   // Discriminant is not a perfect square.
      // Then the Galois group of f(x) is the Frobenius group of order 20.
    GaloisGroupHasOrder20(multiplicity);
    ShowQuinticsRootsRealR(multiplicity);
    return;
  }
  // Test whether M+N*d is a perfect square. In this case the quadratic is
  // reducible. d = Rat4 at this moment.
  CopyBigInt(&Rat4.numerator, &tmp4);
  CopyBigInt(&Rat4.denominator, &tmp1);    // Rat4 = square root of discriminant.
  BigRationalMultiply(&Rat4, &RatN, &Rat1);
  BigRationalAdd(&Rat1, &RatM, &Rat1);
  ForceDenominatorPositive(&Rat1);         // M + N*d
  if (BigIntIsZero(&Rat1.numerator) || Rat1.numerator.sign == SIGN_POSITIVE)
  {
    squareRoot(Rat1.numerator.limbs, tmp4.limbs, Rat1.numerator.nbrLimbs, &tmp4.nbrLimbs);
    tmp4.sign = SIGN_POSITIVE;
    BigIntMultiply(&tmp4, &tmp4, &tmp5);
    BigIntSubt(&tmp5, &Rat1.numerator, &tmp5);
    if (BigIntIsZero(&tmp5))
    {     // Numerator is perfect square.
      squareRoot(Rat1.denominator.limbs, tmp1.limbs, Rat1.denominator.nbrLimbs, &tmp1.nbrLimbs);
      tmp1.sign = SIGN_POSITIVE;
      BigIntMultiply(&tmp1, &tmp1, &tmp5);
      BigIntSubt(&tmp5, &Rat1.denominator, &tmp5);
      if (BigIntIsZero(&tmp5))
      {   // Quadratic has rational solutions.
        CopyBigInt(&Rat3.numerator, &tmp4);    // Rat3 = square root of discriminant
        CopyBigInt(&Rat3.denominator, &tmp1);  // of quadratic equation.
        GaloisGroupHasOrder5(multiplicity);
        return;
      }
    }
  }
  GaloisGroupHasOrder10(multiplicity);
  ShowQuinticsRootsRealR(multiplicity);
  return;
}

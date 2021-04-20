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
#ifndef __FFT_H
#define __FFT_H

typedef struct sComplex
{
  double real;
  double imaginary;
} complex;

struct sCosSin
{
  int Cos[2];
  int Sin[2];
};

#define FFT_LIMB_RANGE  (1 << FFT_LIMB_SIZE)
#define MAX_VALUE_FFT_LIMB (FFT_LIMB_RANGE - 1)
#define QUARTER_CIRCLE (1 << (POWERS_2 - 2))
#define HALF_CIRCLE    (1 << (POWERS_2 - 1))
#define CIRCLE_MASK    ((1 << POWERS_2) - 1)

extern const struct sCosSin cossinPowerOneHalf[15];

#endif


//
// This file is part of Alpertron Calculators.
//
// Copyright 2021 Dario Alejandro Alpern
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
#ifndef _LINKED_BIGNBR_H
#define _LINKED_BIGNBR_H

#define NBR_LINKED_NODES      3000000
#define LIMBS_PER_LINKED_NODE       7  // It must be greater than 2.

struct linkedBigInt
{
  struct linkedBigInt* pstNext;
  int node[LIMBS_PER_LINKED_NODE];
};

void initLinkedBigInt(void);
void getBigIntegerFromLinked(struct linkedBigInt* pstLinkedBigInt, BigInteger* pBigInt);
void setLinkedBigInteger(struct linkedBigInt** pstLinkedBigInt, BigInteger* pBigInt);
bool linkedBigIntIsZero(struct linkedBigInt* pstLinkedBigInt);
bool linkedBigIntIsOne(struct linkedBigInt* pstLinkedBigInt);
bool linkedBigIntIsMinusOne(struct linkedBigInt* pstLinkedBigInt);
void linkedBigIntChSign(struct linkedBigInt* pstLinkedBigInt);
void intToLinkedBigInt(struct linkedBigInt** ppstLinkedBigInt, int value);

#endif
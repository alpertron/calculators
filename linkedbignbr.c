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

#include <stdlib.h>
#include <string.h>
#include "bignbr.h"
#include "linkedbignbr.h"

struct linkedBigInt astLinkedBigInt[NBR_LINKED_NODES];
struct linkedBigInt* pstFirstFree;
bool initLinkedBigIntDone = false;

void initLinkedBigInt(void)
{
  struct linkedBigInt* pstLinkedBigInt = astLinkedBigInt;
  if (initLinkedBigIntDone)
  {            // Linked list already init. Go out.
    return;
  }
  initLinkedBigIntDone = true;
  pstFirstFree = &astLinkedBigInt[0];
  for (int ctr = 0; ctr < (NBR_LINKED_NODES - 1); ctr++)
  {
    pstLinkedBigInt->pstNext = pstLinkedBigInt + 1;
    pstLinkedBigInt++;
  }
  pstLinkedBigInt->pstNext = NULL;
}

void getBigIntegerFromLinked(struct linkedBigInt* pstLinkedBigInt, BigInteger* pBigInt)
{
  struct linkedBigInt* pstCurrentBigInt;
  int* pLimb = (int *)&pBigInt->limbs[LIMBS_PER_LINKED_NODE-2];
  pBigInt->nbrLimbs = pstLinkedBigInt->node[0];
  pBigInt->sign = pstLinkedBigInt->node[1];
  (void)memcpy(pBigInt->limbs, &pstLinkedBigInt->node[2], (LIMBS_PER_LINKED_NODE - 2)*sizeof(int));
  pstCurrentBigInt = pstLinkedBigInt;
  while (pstCurrentBigInt->pstNext != NULL)
  {
    pstCurrentBigInt = pstCurrentBigInt->pstNext;
    (void)memcpy(pLimb, &pstCurrentBigInt->node, LIMBS_PER_LINKED_NODE * sizeof(int));
    pLimb += LIMBS_PER_LINKED_NODE;
  }
}

static void deleteLinkedBigInteger(struct linkedBigInt* pstLinkedBigInt)
{
  if (pstLinkedBigInt != NULL)
  {    // Send this linked big number to linked free area.
    struct linkedBigInt* pstOldFirstFree = pstFirstFree;
    struct linkedBigInt* pstDeletedBigInt = pstLinkedBigInt;
    pstFirstFree = pstLinkedBigInt;
    while (pstDeletedBigInt->pstNext != NULL)
    {   // Loop to find last element of linked big number.
      pstDeletedBigInt = pstDeletedBigInt->pstNext;
    }
    // Make the old linked big number to be the first part of linked free area.
    pstDeletedBigInt->pstNext = pstOldFirstFree;
  }
}

void setLinkedBigInteger(struct linkedBigInt** ppstLinkedBigInt, const BigInteger* pBigInt)
{
  int ctrLimbs;
  struct linkedBigInt* pstLinkedBigInt = *ppstLinkedBigInt;
  deleteLinkedBigInteger(pstLinkedBigInt); // Delete old number.
  *ppstLinkedBigInt = pstFirstFree;
  pstFirstFree->node[0] = pBigInt->nbrLimbs;
  pstFirstFree->node[1] = pBigInt->sign;
  (void)memcpy(&pstFirstFree->node[2], pBigInt->limbs, (LIMBS_PER_LINKED_NODE - 2)*sizeof(int));
  ctrLimbs = LIMBS_PER_LINKED_NODE - 2;
  while (ctrLimbs < pBigInt->nbrLimbs)
  {
    pstFirstFree = pstFirstFree->pstNext;
    (void)memcpy(&pstFirstFree->node, &pBigInt->limbs[ctrLimbs], LIMBS_PER_LINKED_NODE*sizeof(int));
    ctrLimbs += LIMBS_PER_LINKED_NODE;
  }
  pstLinkedBigInt = pstFirstFree->pstNext;
  pstFirstFree->pstNext = NULL;    // End of link of big number.
  pstFirstFree = pstLinkedBigInt;  // First link of first area.
}

bool linkedBigIntIsZero(const struct linkedBigInt* pstLinkedBigInt)
{
  if ((pstLinkedBigInt->node[0] == 1) && (pstLinkedBigInt->node[2] == 0))
  {
    return true;    // Number is zero.
  }
  return false;     // Number is not zero.
}

bool linkedBigIntIsOne(const struct linkedBigInt* pstLinkedBigInt)
{
  if ((pstLinkedBigInt->node[0] == 1) && (pstLinkedBigInt->node[2] == 1) &&
    (pstLinkedBigInt->node[1] == SIGN_POSITIVE))
  {
    return true;    // Number is one.
  }
  return false;     // Number is not one.
}

bool linkedBigIntIsMinusOne(const struct linkedBigInt* pstLinkedBigInt)
{
  if ((pstLinkedBigInt->node[0] == 1) && (pstLinkedBigInt->node[2] == 1) &&
    (pstLinkedBigInt->node[1] == SIGN_NEGATIVE))
  {
    return true;    // Number is minus one.
  }
  return false;     // Number is not minus one.
}

void linkedBigIntChSign(struct linkedBigInt* pstLinkedBigInt)
{
  if ((pstLinkedBigInt->node[0] == 1) && (pstLinkedBigInt->node[2] == 0))
  {    // Value is zero. Do not change sign.
    return;
  }
  if (pstLinkedBigInt->node[1] == SIGN_POSITIVE)
  {
    pstLinkedBigInt->node[1] = SIGN_NEGATIVE;
  }
  else
  {
    pstLinkedBigInt->node[1] = SIGN_POSITIVE;
  }
}

void intToLinkedBigInt(struct linkedBigInt** ppstLinkedBigInt, int value)
{
  int* pNode;
  struct linkedBigInt* pstLinkedBigInt = *ppstLinkedBigInt;
  deleteLinkedBigInteger(pstLinkedBigInt);  // Delete old number.
  pNode = pstFirstFree->node;
  *ppstLinkedBigInt = pstFirstFree;
  if (value >= 0)
  {
    *(pNode+2) = value;
    *(pNode+1) = SIGN_POSITIVE;
  }
  else
  {
    *(pNode+2) = -value;
    *(pNode+1) = SIGN_NEGATIVE;
  }
  *pNode = 1;    // Number of limbs.
  pstLinkedBigInt = pstFirstFree->pstNext;
  pstFirstFree->pstNext = NULL;    // End of link of big number.
  pstFirstFree = pstLinkedBigInt;  // First link of first area.
}

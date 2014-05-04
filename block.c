/* block.c */

/* Copyrighted 24/6/05, Joe Luquette */

#include <stdlib.h>
#include <string.h>

#include "block.h"

BLOCK *
mkblock(void *data)
{
  BLOCK *n;

  n = malloc(sizeof(BLOCK));
  LEFT(n) = NULL;
  RIGHT(n) = NULL;
  PARENT(n) = NULL;
  DATA(n) = data;

  return n;
}

void
rmblock(BLOCK *b)
{
  free(b);
}

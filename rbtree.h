/* rbtree.h */

/* Copyrighted 7/11/2005, Joe Luquette */

#ifndef RBTREE_H
#define RBTREE_H

#include "block.h"

#define RED 0
#define BLACK 1

#define ROOT(_t)             ((_t)->root)
#define BH(_t)               ((_t)->bh)
#define NODENO(_t)           ((_t)->nodeno)
#define BCOMP(_t, _b1, _b2)  ((_t)->bcomp(_b1, _b2))
#define KCOMP(_t, _b, _k)    ((_t)->kcomp(_b, _k))
typedef struct rbtree {
  char *name;
  int bh;       /* The black height of the tree */
  int nodeno;   /* Number of nodes in the tree */
  BLOCK *root;
  int (*bcomp)(BLOCK *, BLOCK *);
  int (*kcomp)(BLOCK *, void *);
} RBTREE;

RBTREE *rbcreate(char *, int (*)(BLOCK *, BLOCK *),
		 int (*)(BLOCK *, void *));
RBTREE *rbclose(RBTREE *);

BLOCK *rbsearch(RBTREE *, void *);
BLOCK *rbmin(RBTREE *);

BLOCK *rbdelete(RBTREE *, BLOCK *);
int rbinsert(RBTREE *, BLOCK *);

#endif /* RBTREE_H */

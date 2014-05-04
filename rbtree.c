/* rbtree.c */

/* Copyrighted 7/11/2005, Joe Luquette */

/*
 * According to Cormen, Leiserson, Rivest and Stein's
 * "Introduction to Algorithms," a red-black tree is defined as
 * follows:
 *    A red-black tree is a binary search tree with one extra
 *    bit* of information per node: its color.  A binary search
 *    tree is considered a red-black tree if the following five
 *    properties hold:
 *        1. Every node in the tree has a color that is either
 *           red or black;
 *        2. The root of the tree is always black;
 *        3. Every external node (leaf node^) of the tree is 
 *           black;
 *        4. Given any red-colored node in the tree, both of its
 *           children must be black;
 *        5. Given any node in the tree, any path from that node
 *           to any leaf which is a descendant of that node con-
 *           tains the same number of black nodes.
 *
 *    * Unfortunately, in practice, we generally gain no benefit
 *      from using only a single bit of information to represent
 *      the color of a node; this is generally due to architecture-
 *      based alignment.
 *    ^ In the case of a red-black tree, every leaf node is a point-
 *      to the NIL node.  Thus, we have a somewhat more specific
 *      meaning in mind for "internal" and "external" nodes.  An ex-
 *      ternal node is any leaf node, and is always a reference to
 *      the special NIL pointer--these nodes do not bear keys and
 *      are not to be considered as data-bearing nodes.  An internal
 *      node is any non-NIL node bearing a key and information.
 *
 * For an in-depth analysis of red-black trees, the proof of the
 * red-black tree height bounds and general commentary on the alg-
 * orithms implemented here, see Cormen et al's "Introduction to
 * Algorithms," Chapter 13.
 *
 * On the implementation side, our red-black trees are made up of
 * the special, generic block structure as defined in block.h.  The
 * only red-black tree-specific structure in our implementation is
 * the RBTREE structure which contains a pointer to the root node,
 * key-based and block-based comparison functions, and diagnostic
 * data such as the tree's keyword, the number of internal nodes,
 * and the black height of the tree as seen from the root node.
 */

#include <stdlib.h>
#include <string.h> /* XXX: For strdup, remove after KEYWORD implementation */
#include "rbtree.h"

/*
 * The special NIL pointer.  This MUST be updated according to
 * the declaration of the BLOCK structure.  NIL is a pointer so
 * that the general macros for working on any node in the tree
 * also work on NIL.
 */
BLOCK _NIL = {
  NULL,
  0,
  BLACK,
  NULL,
  NULL,
  NULL
}, *NIL = &_NIL;

/*
 * Performs a counter-clockwise rotation on 'x' and its right
 * child.
 *
 * Note that since we use a real node as the NIL value, we'll
 * never segfault here due to bad dereferencing; however, it
 * still makes no sense to rotate certain nodes--namely any
 * node whose left child is NIL.  The effect of such a rotat-
 * ion would be moving NIL into the internal area of the tree
 * and rotating the key-bearing node 'x' into the exterior.
 * A symmetrical warning applies for the right_rotate function
 * as well.
 */
void
left_rotate(RBTREE *t, BLOCK *x)
{
  BLOCK *y;
  
  y = RIGHT(x);
  
  PARENT(y) = PARENT(x);
  if (PARENT(x) == NIL)
    ROOT(t) = y;
  else {
    if (x == LEFT(PARENT(x)))
      LEFT(PARENT(x)) = y;
    else
      RIGHT(PARENT(x)) = y;
  }

  RIGHT(x) = LEFT(y);
  if (LEFT(y) != NIL)
    PARENT(LEFT(y)) = x;

  LEFT(y) = x;
  PARENT(x) = y;
}


/*
 * Performs a clock-wise rotation on 'x' and its left child.
 * The comments for left_rotate are worth reading if you think
 * that having a non-NULL NIL sentinel means that you can rot-
 * ate any node you choose.
 */
void
right_rotate(RBTREE *t, BLOCK *x)
{
  BLOCK *y;

  y = LEFT(x);

  PARENT(y) = PARENT(x);
  if (PARENT(x) == NIL)
    ROOT(t) = y;
  else {
    if (x == LEFT(PARENT(x)))
      LEFT(PARENT(x)) = y;
    else
      RIGHT(PARENT(x)) = y;
  }

  LEFT(x) = RIGHT(y);
  if (RIGHT(y) != NIL)
    PARENT(RIGHT(y)) = x;

  RIGHT(y) = x;
  PARENT(x) = y;
}



/*
 * Create a new RBTREE and return it.  Given two arguments, 'a' and
 * 'b', where 'a' is first and 'b' is second, the comparison func-
 * tions should operate as follows:
 *   - If 'a' is less than 'b', then a negative value should be re-
 *     turned;
 *   - If 'a' is equal to 'b', then a value of zero should be re-
 *     turned;
 *   - If 'a' is greater than 'b', then a positive value should be
 *     returned.
 * In this implementation, there are two different types of compar-
 * ison functions: a "block comparator" and a "key comparator".
 * The block comparator should satify the above requirements when
 * the arguments 'a' and 'b' are both BLOCKs; the key comparator
 * should satisfy the above requirements when 'a' is a BLOCK and
 * 'b' is a void pointer, and 'b' holds the same type of data as
 * the DATA field of 'a'.
 *
 * The macro for the comparison functions aren't used here because
 * they are written to call the function, not just reference the
 * structure member.  I think this is more pleasing on the eyes
 * everywhere else, and we only pay with this one inconsistency.
 * It should be readily apparent to anyone who reads rbtree.h that
 * we have no structure member that references the comparison func-
 * tion members directly, and thus no one should be particularly sur-
 * prised if they change the structure member nomenclature and notice
 * some "undeclared member" errors in their C compiler's output.
 */
RBTREE *
rbcreate(char *name, int (*bcomp)(BLOCK *, BLOCK *),
	 int (*kcomp)(BLOCK *, void *))
{
  RBTREE *nt;

  nt = malloc(sizeof(RBTREE));
  NAME(nt) = strdup(name);
  BH(nt) = 0;
  NODENO(nt) = 0;
  ROOT(nt) = NIL;
  nt->bcomp = bcomp;
  nt->kcomp = kcomp;

  return nt;
} 



/* Assistant recursive function for rbclose() */
void
rbclose_recurse(BLOCK *n)
{
  if (n == NIL)
    return;

  rbclose_recurse(LEFT(n));
  rbclose_recurse(RIGHT(n));

  rmblock(n);
}



/*
 * Recursively deletes every node in the tree, frees the top
 * level RBTREE structure, and returns NULL so that the caller
 * may destroy his own reference to the tree using a call like:
 *   RBTREE *t;
 *
 *   t = rbcreate(strcmp, free);
 *   . . .
 *   t = rbclose(t);
 */
RBTREE *
rbclose(RBTREE *t)
{
  rbclose_recurse(ROOT(t));
  free(t);

  return NULL;
}



/*
 * Perform a red-black tree balancing after inserting a new
 * node.  Note that since both children of any red node must
 * be black and that the root must be black; that we can
 * assume that if 'n' has a red parent, then n must have a
 * non-NIL grandparent.  Note that if the uncle of 'n' is not
 * red, then we color the parent of 'n' black, and thus cause
 * the loop to terminate.  Thus, the majority of balance op-
 * erations will be coloring two nodes and setting 'n' to its
 * grandparent; the loop will terminate after performing either
 * one or two rotates.
 *
 * Also note the final lines: in "Introduction to Algorithms",
 * the red-black tree fix function colored the root to black
 * unconditionally.  This is, of course, more simple and some-
 * times faster; but we are interested in keeping track of the
 * black height of our tree.  If you pay special attention to
 * the while loop, you will notice that it is written so that
 * it constantly preserves black heights--that is, the while
 * loop itself never alters the black height of the tree.  The
 * only time the black height of the tree is altered is when
 * the while loop exits with 'n' equal to the root of 't', and
 * with 'n' colored red.  'n' is then forced to black and the
 * black height is incremented.  The black height plays no fun-
 * ctional role in our red-black tree implementation, and is
 * only used for pretty output.
 */ 
void
rbinsert_fix(RBTREE *t, BLOCK *n)
{
  BLOCK *uncle;

  while (COLOR(PARENT(n)) == RED) {
    if (PARENT(n) == LEFT(PARENT(PARENT(n)))) {
      uncle = RIGHT(PARENT(PARENT(n)));
      if (COLOR(uncle) == RED) {           /* 1. */
	COLOR(uncle) = BLACK;
	COLOR(PARENT(n)) = BLACK;
	n = PARENT(PARENT(n));
	COLOR(n) = RED;
      } else {                             /* 2. */
	if (n == RIGHT(PARENT(n))) {
	  n = PARENT(n);
	  left_rotate(t, n);
	}
	COLOR(PARENT(n)) = BLACK;          /* 3. */
	COLOR(PARENT(PARENT(n))) = RED;
	right_rotate(t, PARENT(PARENT(n)));
      }
    } else {
      uncle = LEFT(PARENT(PARENT(n)));
      if (COLOR(uncle) == RED) {           /* 1 (Symmetric) */
	COLOR(PARENT(n)) = BLACK;
	COLOR(uncle) = BLACK;
	n = PARENT(PARENT(n));
	COLOR(n) = RED;
      } else {                             /* 2 (Symmetric) */
	if (n == LEFT(PARENT(n))) {
	  n = PARENT(n);
	  right_rotate(t, n);
	}
	COLOR(PARENT(n)) = BLACK;          /* 3 (Symmetric) */
	COLOR(PARENT(PARENT(n))) = RED;
	left_rotate(t, PARENT(PARENT(n)));
      }
    }
  }

  if (COLOR(ROOT(t)) == RED) {
    COLOR(ROOT(t)) = BLACK;
    ++BH(t);
  }
}


 
/*
 * Insert the block 'n' into the red-black tree 't'.  This
 * algorithm is extremely similar to the algorithm given in
 * the venerable CLRS, Chapter 13, page 280.  Our rbinsert,
 * however, explicitly disallows duplicate entries.
 *
 * Returns non-zero when the block could not be inserted (curr-
 * ently this only happens when a block with this key already
 * exists in the tree); otherwise returns zero.
 */
int
rbinsert(RBTREE *t, BLOCK *n)
{
  int res;
  BLOCK *b, *p;

  p = NIL;
  b = ROOT(t);
  while (b != NIL) {
    p = b;
    res = BCOMP(t, n, b);
    if (res < 0)
      b = LEFT(b);
    else if (res > 0)
      b = RIGHT(b);
    else
      return -1;
  }

  PARENT(n) = p;
  if (p == NIL)
    ROOT(t) = n;
  else {
    if (BCOMP(t, n, p) < 0)
      LEFT(p) = n;
    else
      RIGHT(p) = n;
  }

  LEFT(n) = NIL;
  RIGHT(n) = NIL;
  COLOR(n) = RED;

  rbinsert_fix(t, n);
  
  ++NODENO(t);

  return 0;
}



/*
 * Given a key 'k', search the red-black tree 't' for a block
 * with key 'k' and return it.  If no such block exists in 't',
 * return NULL.
 */
BLOCK *
rbsearch(RBTREE *t, void *k)
{
  int res;
  BLOCK *n;

  n = ROOT(t);
  while (n != NIL) {
    res = KCOMP(t, n, k);
    if (res < 0)
      n = LEFT(n);
    else if (res > 0)
      n = RIGHT(n);
    else
      return n;
  }

  return NULL;
}
 
/* Return the block with minimum key in t and NULL if t is empty */
BLOCK *
rbmin(RBTREE *t)
{
  BLOCK *n;

  n = ROOT(t);
  if (n == NIL)
    return NULL;

  while (LEFT(n) != NIL)
    n = LEFT(n);

  return n;
}


/*
 * Work upwards from 'n' dealing with all 8 cases (4 of which are
 * symmetric) that may occur to violate the properties of a red-
 * black tree.  Note that rbdelete_fix assumes that 'n' is colored
 * black.  There is no reason to rebalance a red-black tree if the
 * spliced node was red as it will cause no violations of the red-
 * black properties!
 *
 * Note that variable names are chosen as: 'n' for "node"; 's' for
 * "sibling."
 *
 * The following comments apply to the first of the two major cases:
 * where 'n' is the left child of its parent, and its sibling 's' is
 * the right child.  Using symmetry, one can use these explanations
 * to explain the symmetric case whre 'n' is the right and 's' the
 * left child.
 * 1. Condition 1 falls through to the following conditions after
 *    performing recolorings and a rotation that forces 's' to be
 *    black.  Since 's' was red at the beginning of 1, by property
 *    4 both of its children must be black.  During the left rot-
 *    ation, 'n's parent remains its parent, and it acquires (as its
 *    right child) 's's left child; which we have just shown must be
 *    black.  's' is updated after the rotation to point to 'n's sib-
 *    ling as required.
 * 2. If case 1 falls through to case 2, then the loop is ended,
 *    since case 1 colors 'n's parent red, which is preserved after
 *    the rotation.
 * 3. Implies LEFT(s) is red since case 2 did not execute, and case
 *    1 assures us that 's' is black.  The aim of case 3 is to morph
 *    itself into case 4, so that case 4 can make assumptions about
 *    the location of the guaranteed red pointer.
 * 4. Now we're guaranteed to be able to finish our balancing proce-
 *    dure.  We color the sibling 's' red; color its red child (which
 *    is guaranteed through case 3 to be its right child) black; then
 *    color the parent of 'n' and 's' black.  After this, we perform
 *    a left rotation causing 's' to become the new parent of the
 *    subtree.  Note that if the parent of 'n' and 's' (prior to the
 *    rotation) is already black, then we're effectively decreasing
 *    the black height of the entire subtree.
 */
void
rbdelete_fix(RBTREE *t, BLOCK *n)
{
  BLOCK *s;

  while (n != ROOT(t) && COLOR(n) == BLACK) {
    if (n == LEFT(PARENT(n))) {
      s = RIGHT(PARENT(n));
      if (COLOR(s) == RED) {             /* 1. */
	COLOR(s) = BLACK;
	COLOR(PARENT(n)) = RED;
	left_rotate(t, PARENT(n));
	s = RIGHT(PARENT(n));
      }
      if (COLOR(LEFT(s)) == BLACK &&
	  COLOR(RIGHT(s)) == BLACK) {    /* 2. */
	COLOR(s) = RED;
	n = PARENT(n);
      } else {
	if (COLOR(RIGHT(s)) == BLACK) {  /* 3. */
	  COLOR(LEFT(s)) = BLACK;
	  COLOR(s) = RED;
	  right_rotate(t, s);
	  s = RIGHT(PARENT(n));
	}
	COLOR(s) = COLOR(PARENT(n));     /* 4. */
	COLOR(PARENT(n)) = BLACK;
	COLOR(RIGHT(s)) = BLACK;
	left_rotate(t, PARENT(n));
	break;
      }
    } else {
      s = LEFT(PARENT(n));
      if (COLOR(s) == RED) {             /* 1. (Symmetric) */
	COLOR(s) = BLACK;
	COLOR(PARENT(n)) = RED;
	right_rotate(t, PARENT(n));
	s = LEFT(PARENT(n));
      }
      if (COLOR(LEFT(s)) == BLACK &&
	  COLOR(RIGHT(s)) == BLACK) {    /* 2. (Symmetric) */
	COLOR(s) = RED;
	n = PARENT(n);
      } else {
	if (COLOR(LEFT(s)) == BLACK) {   /* 3. (Symmetric) */
	  COLOR(RIGHT(s)) = BLACK;
	  COLOR(s) = RED;
	  left_rotate(t, s);
	  s = LEFT(PARENT(n));
	}
	COLOR(s) = COLOR(PARENT(n));     /* 4. (Symmetric) */
	COLOR(PARENT(n)) = BLACK;
	COLOR(LEFT(s)) = BLACK;
	right_rotate(t, PARENT(n));
	break;
      }
    }
  }

  if (n == ROOT(t))
    --BH(t);

  COLOR(n) = BLACK;
}



/*
 * This function ASSUMES that 'n' is a valid BLOCK that is already
 * in the tree 't'!  If you want a safe version of rbdelete() (or
 * perhaps an rbdelete-by-key instead of rbdelete-by-block), write
 * a wrapper using rbsearch().
 *
 * Remove the node 'n' from the red-black tree 't'.  The rbdelete
 * algorithm begins by selecting a node 's' (for "spliced") to be
 * spliced out of the red-black tree.  In the case where 'n' has at
 * least one NIL child, 's' and 'n' are synonymous; otherwise 'n'
 * has a non-NIL right child, and we set 's' to the successor of
 * 'n'.  It is worthwhile to point out here that we deviate from
 * Cormen et. al's algorithm since we only care about a successor
 * of 'n' in the case where RIGHT(n) is not NIL.  Thus, we do not
 * need a full implementation of Tree-Successor as in "Introduction
 * to Algorithms"--we need only the Tree-Minimum portion of Tree-
 * Successor, which is implemented here in a single for loop.
 * 
 * The node 'c' is the only child (and may in fact be NIL in the
 * case where 'n' has no children) of the to-be-spliced node 's'.
 *
 * Note that in the case where 's' is not 'n', 's' is spliced out
 * for the sole purpose of replacing 'n'; thus once 's' has been
 * removed from the tree, we must redirect all 6 links (3 leaving
 * and 3 entering) around 'n' to instead point to 's'.
 *
 * Note that if the node we spliced out was red, then we cannot
 * have violated any properties of the red-black tree:
 *   1. All nodes are still either red or black;
 *   2. The root is still black (the root must have been black to
 *      begin with if 't' was a valid red-black tree, thus if we
 *      removed a red node, it could not have been the root).
 *   3. All leaves are still black (NIL, to be more precise);
 *   4. All red nodes still have only black children;
 *   5. The black height from any node to any leaf is unaltered,
 *      since the node deleted was a red node.
 */
BLOCK *
rbdelete(RBTREE *t, BLOCK *n)
{
  unsigned char color;
  BLOCK *s, *c;

  if (LEFT(n) == NIL || RIGHT(n) == NIL)
    s = n;
  else
    /* Find the successor assuming a non-NIL RIGHT(n) */
    for (s = RIGHT(n); LEFT(s) != NIL; s = LEFT(s))
      continue;

  color = COLOR(s);

  /* Select the only child of 's' (which may be NIL!) */
  if (LEFT(s) != NIL)
    c = LEFT(s);
  else
    c = RIGHT(s);

  /* Now remove 's' from the tree */
  PARENT(c) = PARENT(s);
  if (PARENT(s) == NIL)
    ROOT(t) = c;
  else {
    if (s == LEFT(PARENT(s)))
      LEFT(PARENT(s)) = c;
    else
      RIGHT(PARENT(s)) = c;
  }

  /* If 's' isn't 'n', then 's' was spliced to put in 'n's place */
  if (s != n) {
    COLOR(s) = COLOR(n);

    PARENT(s) = PARENT(n);
    if (PARENT(n) == NIL)
      ROOT(t) = s;
    else {
      if (n == LEFT(PARENT(n)))
	LEFT(PARENT(n)) = s;
      else
	RIGHT(PARENT(n)) = s;
    }
    
    LEFT(s) = LEFT(n);
    PARENT(LEFT(n)) = s;
    
    RIGHT(s) = RIGHT(n);
    PARENT(RIGHT(n)) = s;
  }
  
  if (color == BLACK)
    rbdelete_fix(t, c);

  --NODENO(t);

  return n;
}

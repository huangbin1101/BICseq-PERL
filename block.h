/* block.h */

/* Copyrighted 24/6/05, Joe Luquette */

#ifndef BLOCK_H
#define BLOCK_H

#define NAME(_x)   ((_x)->name)
#define DATA(_b)   ((_b)->data)
#define SIZE(_b)   ((_b)->size)
#define PREV(_b)   ((_b)->prev)
#define NEXT(_x)   ((_x)->next)
#define PARENT(_b) ((_b)->parent)
#define LEFT(_b)   ((_b)->prev)
#define RIGHT(_b)  ((_b)->next)
#define COLOR(_b)  ((_b)->color)
typedef struct block {
  void *data;
  size_t size;
  unsigned char color;
  struct block *parent;
  struct block *prev;
  struct block *next;
} BLOCK;

BLOCK *mkblock(void *);
void rmblock(BLOCK *);

#endif /* BLOCK_H */

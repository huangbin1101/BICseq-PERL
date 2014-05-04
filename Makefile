CC = gcc
CFLAGS = -g -O2 -lm -Wall


BIC-seq: BICseq.o pos_cnt_lst.o read_v0.1.o bic-seq.o block.o rbtree.o gamma.o
	$(CC) $(CFLAGS)  BICseq.o pos_cnt_lst.o read_v0.1.o bic-seq.o block.o rbtree.o gamma.o -o BICseq/BIC-seq



BICseq.o: BICseq/pos_cnt_lst.h BICseq/read.h BICseq/bin.h BICseq/ll.h  BICseq/block.h BICseq/rbtree.h BICseq/gamma.h
	$(CC) $(CFLAGS) -c BICseq/BICseq.c BICseq/pos_cnt_lst.c BICseq/read_v0.1.c BICseq/bic-seq.c BICseq/block.c BICseq/rbtree.c BICseq/gamma.c

bic-seq.o: BICseq/bin.h BICseq/ll.h BICseq/block.h BICseq/rbtree.h
	$(CC) $(CFLAGS) -c BICseq/bic-seq.c BICseq/block.c BICseq/rbtree.c BICseq/gamma.c

block.o: BICseq/block.h
	$(CC) $(CFLAGS) -c BICseq/block.c

rbtree.o: BICseq/rbtree.h BICseq/block.h
	$(CC) $(CFLAGS) -c BICseq/rbtree.c BICseq/block.c

gamma.o: BICseq/gamma.h
	$(CC) $(CFLAGS) -c BICseq/gamma.c

pos_cnt_lst.o: BICseq/pos_cnt_lst.h
	$(CC) $(CFLAGS) -c BICseq/pos_cnt_lst.c

read_v0.1.o: BICseq/read.h
	$(CC) $(CFLAGS) -c BICseq/read_v0.1.c  


clean:
	rm -rf *.o


#ifndef BIN_H
#define BIN_H

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ll.h"
#include "rbtree.h"
#include "gamma.h"

/* The smallest data type the computation handles */
typedef struct bin {
	int tumor;
	int total;
	double freq;
	int start;
	int end;
	struct bin *prev;
	struct bin *next;
} BIN;

typedef LL(BIN) BIN_LIST;

typedef struct SEG_PERMUTE_t{
        BIN_LIST *bins_perm;
        int size;
        } *SEG_PERMUTE;



BIN * bin_new(int tumor, int total, double freq, int start, int end);
void bic_seq(int paired);
void bic_seq_auto(int nbins,double FP, int paired);
void set_lambda(double lambda_in); /*set the penalty parameter*/
void set_BinList(BIN_LIST bins_in); /*set the global variable bins in bic-seq.c as bins_in*/
BIN_LIST get_BinList(); /*get the global variable bins*/

void set_totalreadcount(double tumor, double normal); /*set the total number of reads as tumor and normal (also the total number of reads)*/
void BIN_LIST_print(BIN_LIST bs, FILE *out);


SEG_PERMUTE SEG_PERMUTE_create(int size);
void SEG_PERMUTE_destroy(SEG_PERMUTE segs);
void print_SEG_PERMUTE(SEG_PERMUTE segs, FILE *out);



SEG_PERMUTE bic_seq_perm(BIN_LIST bins_in, double p, double FP,int num_perm ,int auto_select); /*use permutation to get the FDR instead of resampling*/

#endif

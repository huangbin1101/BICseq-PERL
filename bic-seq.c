/*
 * bic-seq.c
 *
 * Concept by Ruibin Xi, algorithm by Joe Luquette.
 * Added two new functions: 1. perform permutation for FDR estimate. 2. choose lambda based on the allowed number of type I errors in the merging process
 */

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <sys/time.h>


#include "ll.h"
#include "rbtree.h"
#include "bin.h"
#include "gamma.h"


SEG_PERMUTE SEG_PERMUTE_create(int size)
{	BIN_LIST bins_tmp = ll_new();
	int i;
	SEG_PERMUTE segs = malloc(sizeof(struct SEG_PERMUTE_t));
	if(segs == NULL){fprintf(stderr,"Error, memory allocation failed in SEG_PERMUTE_create\n"); exit(1);}

	segs->size = (size>=1 ? size: 1);
	segs->bins_perm = malloc(sizeof(BIN_LIST)*segs->size);

	if(segs->bins_perm == NULL){fprintf(stderr,"Error, memory allocation failed in SEG_PERMUTE_create\n"); exit(1);}

	for(i=0;i<size;i++){
		segs->bins_perm[i] = bins_tmp;
		}

	return segs;
}


void SEG_PERMUTE_destroy(SEG_PERMUTE segs)
{	int i;
	if(segs==NULL) return;
	if(segs->bins_perm!=NULL&&segs->size>0){
		for(i=0;i<segs->size;i++)
			ll_dealloc(segs->bins_perm[i]);

		free(segs->bins_perm);
		}
	segs->bins_perm = NULL;
	segs->size = 0;
	free(segs);
	return;
}

void print_SEG_PERMUTE(SEG_PERMUTE segs, FILE *out)
{	int i;
	BIN_LIST bins_tmp;
	BIN *b;
	if(segs==NULL||segs->size<=0||segs->bins_perm==NULL) return;

	fprintf(out,"resample\tcaseRead\ttotalRead\tcontrolFreq\tstart\tend\n");
	for(i=0;i<segs->size;i++){
		bins_tmp = segs->bins_perm[i];
		b = ll_head(bins_tmp);
		while(b!=NULL){
			fprintf(out,"%d\t%d\t%d\t%lf\t%d\t%d\n",i+1,b->tumor, b->total, b->freq, b->start, b->end);
			b = ll_next(b);
			}
		}
	return;
}


static double N = 0,N_tumor_X=0, N_normal_X=0;
static double lambda = 2.0;

static BIN_LIST bins = ll_new();

void set_lambda(double lambda_in)
{	lambda = lambda_in;
	return;
}

void set_BinList(BIN_LIST bins_in)
{	bins = bins_in;
}

BIN_LIST get_BinList()
{	return bins;
}

void set_totalreadcount(double tumor, double normal)
{	N_tumor_X = tumor;
	N_normal_X = normal;
	N = tumor+normal;
	return;
}

BIN *
bin_new(int tumor, int total, double freq, int start, int end)
{
	BIN *nb = malloc(sizeof(BIN));

	if (nb == NULL) {
		fprintf(stderr,"bin_new: malloc: %s\n", strerror(errno));
		exit(1);
	}

	nb->tumor = tumor;
	nb->total = total;
	nb->freq = freq;
	nb->start = start;
	nb->end = end;
	nb->prev = NULL;
	nb->next = NULL;

	return nb;
}

void
bin_print(BIN *b, FILE *out)
{
	fprintf(out,"%d\t%d\t%lf\t%d\t%d\n",
		b->tumor, b->total, b->freq, b->start, b->end);
}

void BIN_LIST_print(BIN_LIST bs, FILE *out)
{	BIN *b;
	b = ll_head(bs);
	while(b!=NULL){
		bin_print(b, out);
		b = ll_next(b);
		}
	return;
}

/*
 * Windows are a doubly linked list for fast next() and prev()
 * and are stored in a red-black tree keyed by BIC difference.
 *
 * Some operations on windows:
 *     - bic_diff(window): compute the difference between the
 *       current BIC of the window and the BIC after running
 *       merge() on it.
 *     - merge(window): merge window.  Merging a window requires
 *       extensive modification of the (b - 1) windows to the left
 *       and right of the merged window.  See merge() for details.
 *     - rbinsert(index, window): insert window into a fast index
 *       keyed by the BIC difference.
 *     - rbdelete(index, window): removes window from the fast index
 *     - ll_next(window): next window in the global window list.  The
 *       global list is ordered by increasing genomic coordinates.
 *     - ll_prev(window): previous window in the global list.
 *     - window_print(window)
 *     - window_new(b, start): a new window with b bins, starting at
 *       bin 'start'.
 */
typedef struct window {
	int id;               /* XXX: hack to allow duplicate keys */
	double bic_diff;
	int b;                /* Number of bins in window */
	BIN *bin_start;       /* The first bin in this window */
	BLOCK *idx_entry;     /* The entry stored in the BIC idx */
	struct window *next;
	struct window *prev;
} WINDOW;

typedef LL(WINDOW) WINDOW_LIST;
WINDOW_LIST windows = ll_new();

void
window_print(WINDOW *w)
{
	int i, start, end;
	BIN *b;

	/* Find start/end */
	b = w->bin_start;
	start = b->start;
	end = b->end;
	for ( ; ll_next(b); b = ll_next(b))
		if (b->end > 0)
			end = b->end;

	fprintf(stderr,"WINDOW (b=%d): (%d,%d), merge gain: %lf\n",
		w->b, start, end, w->bic_diff);
	b = w->bin_start;
	for (i = 0; i < w->b; ++i) {
		bin_print(b,stdout);
		b = ll_next(b);
	}
}

/*
 * Compute w's BIC.  Let:
 *    k be the number of parameters in w's model,
 *    N be the number of Bernoulli trials observed in w,
 *    L be the value of the likelihood function evaluated at MLE
 * then BIC = k*ln(N) - 2*ln(L).
 *
 * bins is a pointer to the first bin in a window, b is the number
 * of bins to be included in this window.
 */
double
compute_bic(BIN *bins, int b)
{
	int i;
	double logL, logN;
	BIN *x;

	logL = 0;  /* ln(P(mle)) */
	x = bins;
	for (i = 0; i < b; ++i) {
		if (x->tumor != 0 && x->tumor != x->total) {
			/* Compute log(L(mle)), L=likelihood function */
			logL += x->tumor * log(x->freq) +
				(x->total - x->tumor) * log(1.0 - x->freq);
		} /* else logL is 0 */

		x = x->next;
	}

	/* Some windows are empty.  In that case, log(N) is set to 0 */
	logN = N == 0 ? 0 : log(N);

	return b*logN*lambda - 2*logL;
}

/* Compute the difference in BIC between w as is and merge(w) */
double
bic_diff(WINDOW *w)
{
	int i;
	BIN *x, merged_bin;
	WINDOW merged_w;

	/* Create a fake bin */
	merged_bin.tumor = 0;
	merged_bin.total = 0;
	merged_bin.next = NULL;
	merged_bin.prev = NULL;
	x = w->bin_start;
	for (i = 0; i < w->b; ++i) {
		merged_bin.tumor += x->tumor;
		merged_bin.total += x->total;
		x = ll_next(x);
	}
	merged_bin.freq = merged_bin.total > 0 ?
		(double)merged_bin.tumor / (double)merged_bin.total : 0;

	merged_w.b = 1;  /* Fake window representing merge(w) */

	return compute_bic(&merged_bin, 1) - compute_bic(w->bin_start, w->b);
}

WINDOW *
window_new(int b, BIN *bin_start)
{
	static long last_id=0;
	WINDOW *nw = malloc(sizeof(WINDOW));

	if (nw == NULL) {
		fprintf(stderr,"window_new: malloc: %s\n", strerror(errno));
		exit(1);
	}

	nw->id = last_id++;
	nw->b = b;
	nw->bin_start = bin_start;
	nw->idx_entry = mkblock(nw);
	nw->prev = NULL;
	nw->next = NULL;

	nw->bic_diff = bic_diff(nw);

	return nw;
}

WINDOW *
window_free(WINDOW *w)
{
	rmblock(w->idx_entry);
	free(w);
	return NULL;
}

/* Must return < 0 if a < b; > 0 if a > b; and =0 if a == b */
/* XXX: hack for equal keys: break tie by a unique ID assigned to */
/* each window. */
int
bcomp(BLOCK *a, BLOCK *b)
{
	double d;
	WINDOW *x = (WINDOW *)DATA(a);
	WINDOW *y = (WINDOW *)DATA(b);

	/* break ties on bic_diff using an ID */
	d = x->bic_diff - y->bic_diff;
	d = d == 0 ? x->id - y->id : d;
	return d < 0 ? -1 : (d > 0 ? 1 : 0);
}

/*
 * After deleting bins in merge(), a window on the right-most edge
 * of the chromosome may no longer be valid.  This function checks
 * the bin list to ensure that from w's start position, there are
 * at least w->b bins for w to represent.
 */
int
window_can_exist(WINDOW *w)
{
	int i;
	BIN *bin;

	bin = w->bin_start;
	for (i = 0; i < w->b && bin != NULL; ++i)
		bin = ll_next(bin);

	return i == w->b;
}

void
merge(RBTREE *idx, WINDOW *w)
{
	int i;
	BIN *b, *new_bin, *tmpb;
	WINDOW *z, *tmpw;

//if (w->bin_start->start > 191550000 && w->bin_start->end < 191750000) window_print(w);

	/*
	 * Step 1: merge the (b-1) bins right of the first bin in window
	 * w into the first bin.
	 * NOTE: Ruibin's R script outputs start=0 and end=0 for bins
	 * with 0 total reads.  So we only extend the bin's end point
	 * when the end is not 0.
	 */
	new_bin = w->bin_start;
	b = ll_next(w->bin_start);
	for (i = 1; i < w->b && b != NULL; ++i) {
		if (new_bin->start == 0)
			new_bin->start = b->start;
		new_bin->tumor += b->tumor;
		new_bin->total += b->total;
		if (b->end > 0)
			new_bin->end = b->end;
		tmpb = ll_next(b);
		ll_delete(bins, b);
		free(b);
		b = tmpb;
	}
	new_bin->freq = new_bin->total > 0 ?
		(double)new_bin->tumor / (double)new_bin->total : 0;

	/*
	 * Step 2: delete the (b-1) windows right of w.  They no longer
	 * exist since the (b-1) starting bins they correspond to no
	 * longer exist.  These windows are not merged; they're simply
	 * deleted.
	 * NOTE: there may be fewer than (b-1) windows right of w.
	 */
	z = ll_next(w);
	for (i = 1; i < w->b && z != NULL; ++i) {
		rbdelete(idx, z->idx_entry);
		tmpw = ll_next(z);
		ll_delete(windows, z);
		z = window_free(z);
		z = tmpw;
	}

	/*
	 * Step 3: recompute the BIC for w and the (b-1) windows left
	 * of w.  Recomputing the BIC implies adjusting the BIC index:
	 * delete the entry in the index for the old BIC and insert a
	 * new entry for the newly computed BIC.
	 * NOTE: there may be fewer than (b-1) windows left of w.
	 * MOREOVER: w and some of the (b-1) windows left of it may not
	 * have reason to exist after the bin merge.
	 * XXX: for right now, we'll use a non-optimal O(b) solution to
	 * determine which windows can stick around.  we'll fix this later
	 * if it's a performance problem.
	 */
	z = w;
	for (i = 0; i < w->b && z != NULL; ++i) {
		tmpw = z;
		rbdelete(idx, z->idx_entry);

		if (window_can_exist(z)) {
			z->bic_diff = bic_diff(z);
			if (rbinsert(idx, z->idx_entry) < 0) {
				fprintf(stderr,"merge: key already exists in index\n");
				exit(1);
			}
		} else {
			ll_delete(windows, z);
			z = window_free(z);
		}

		z = ll_prev(tmpw);
	}
}

BLOCK *
rbmin_exhaustive(RBTREE *t, BLOCK *x)
{
	BLOCK *m, *z;
	extern BLOCK *NIL;

	m = x;
	if (LEFT(x) != NIL) {
		z = rbmin_exhaustive(t, LEFT(x));
		m = BCOMP(t, m, z) < 0 ? m : z;
	}

	if (RIGHT(x) != NIL) {
		z = rbmin_exhaustive(t, RIGHT(x));
		m = BCOMP(t, m, z) < 0 ? m : z;
	}

	return m;
}

WINDOW *
window_min()
{
	WINDOW *w, *m;

	m = ll_head(windows);
	for (w = ll_head(windows); w != NULL; w = ll_next(w))
		if (w->bic_diff < m->bic_diff)
			m = w;

	return m;
}
 
/*
 * Run bic_seq on the current list of bins, attempting to merge
 * b bins at a time.  Continues merging until no bin merge will
 * yield a lower BIC.
 *
 * Returns the number of successful merges.  If 0 merges occur,
 * no further levels should be attempted.
 */
int
bic_seq_level(int b)
{
	int i, p, q;
	int nmerged;
	BIN *x;
	BLOCK *m;
	RBTREE *idx;
	WINDOW *w;

	/* Compute the current BIC (this is just for display) */
	//fprintf(stderr,"\t%d-bins: total bins: %d\n", b, ll_length(bins));
	//fprintf(stderr,"\t%d-bins: initial BIC: %lf\n",
	//	b, compute_bic(ll_head(bins), ll_length(bins)));

	/*
	 * Compute a list of windows (based on b) and index them.  To
	 * index, we compute the BIC of the window before and after a
	 * a simulated merge() and record the difference.  The diff-
	 * erence keys the index.  The most negative value is the best.
	 */
	idx = rbcreate("BIC tree", bcomp, NULL);  /* No block-key comparator */
	x = ll_head(bins);
	p = q = 0;
	for (i = 0; i < ll_length(bins) - b+1; ++i) {
		w = window_new(b, x);
		w->bic_diff < 0 ? ++p : ++q;      /* Just for display */
		ll_append(windows, w);

		if (rbinsert(idx, w->idx_entry) < 0) {
			fprintf(stderr,"rbinsert: key already exists in tree\n");
			exit(1);
		}
		x = ll_next(x);
	}
	//printf("\t%d-bins: %d windows, %d with BIC improvement, %d without\n"
	//	"\t%d-bins: index on BIC improvement: %d nodes, bl-height=%d\n",
	//	b, ll_length(windows), p, q,
	//	b, NODENO(idx), BH(idx));

	/* Merge until no window can lower the BIC of the model */
	nmerged = 0;
	for ( ; ; ) {
		/* m == NULL means idx is empty */
		m = rbmin(idx);
		if (m == NULL)
			break;

		w = (WINDOW *)DATA(m);
		if (w != NULL && w->bic_diff < 0) {
			merge(idx, w);
			++nmerged;
		} else
			break;
	}

	ll_dealloc(windows);
	idx = rbclose(idx);
	return nmerged;
}


BIN_LIST dup_bin_lst(BIN_LIST bins_in) /*duplicate a bin list*/
{	BIN_LIST bins_new = ll_new();
	BIN *b;

	b = ll_head(bins_in);
	while(b!=NULL){
		ll_append(bins_new, bin_new(b->tumor, b->total, b->freq, b->start, b->end));
		b = ll_next(b);
		}

	return bins_new;
}

BIN_LIST random_assign(BIN_LIST bins_in, double p) /*random assign the reads as tumor or normal, p is the probability that a read is tumor read*/
{       BIN_LIST bins_new = ll_new();
        BIN *b;
	int total, tumor;
	double freq;
        struct timeval tv;
        int seed;
	
        gettimeofday(&tv, NULL);
        seed = tv.tv_sec * 1000000 + tv.tv_usec;
	seed_set(seed);

        b = ll_head(bins_in);
        while(b!=NULL){
		total = b->total;
		tumor = (int) rbinom(p, (double) total);
		freq = ((double) tumor)/((double) total);
		
                ll_append(bins_new, bin_new(tumor, total, freq, b->start, b->end));
                b = ll_next(b);
                }
	return bins_new;
}



int merge_refine(int num_ini_bins,double FP, double siglevel)
/*siglevel: the significance level for Bonferroni corrected pvalue, siglevel should be between 0 and 1, though this function does not check for this
  FP: number of false positve allowed; this number should be postive;
*/
{	BIN *b, *btmp;
	double pvalue1, pvalue2, tmp1, tmp2;
	int nmerged=0;

	b = ll_head(bins);

	if(FP<=0.0||siglevel<=0.0||siglevel>1.0){fprintf(stderr,"Invalide parameters in merge_refine\n");exit(1);}

	while(b!=NULL){
		if(b->prev==NULL) pvalue1 = -1.0;
		else {
			tmp1 = 2*pbinom(b->tumor,b->total,b->prev->freq,0)*num_ini_bins/FP;
			tmp2 = 2*pbinom(b->tumor,b->total,b->prev->freq,1)*num_ini_bins/FP;
			pvalue1 = tmp1<tmp2? tmp1 : tmp2;			
			}

		if(b->next==NULL) pvalue2 = -1.0;
		else {
			tmp1 = 2*pbinom(b->tumor,b->total,b->next->freq,0)*num_ini_bins/FP;
			tmp2 = 2*pbinom(b->tumor,b->total,b->next->freq,1)*num_ini_bins/FP;
			pvalue2 = tmp1<tmp2? tmp1 : tmp2;
			}
		
		btmp = ll_next(b);
		if(pvalue1>=pvalue2 && pvalue1>siglevel &&b->prev!=NULL){ //merge b with the previous node;
			b->prev->total += b->total;
			b->prev->tumor += b->tumor;
			b->prev->end = b->end;
			b->prev->freq = ((double) b->prev->tumor)/((double) b->prev->total);
			ll_delete(bins, b);
			free(b); b= NULL;
			nmerged++;	
			}
		else if(pvalue2>pvalue1 && pvalue2 > siglevel && b->next!=NULL){//merge b with the next node;
			b->next->total += b->total;
			b->next->tumor += b->tumor;
			b->next->start = b->start;
			b->next->freq = ((double) b->next->tumor)/((double) b->next->total);
			ll_delete(bins, b);
                        free(b); b= NULL;
			nmerged++;
			}
		b = btmp;
		}

	return nmerged;
}


void
bic_seq(int paired) 
{
	int nmerged = 0,total_merged=0;
	int level = 2;
        int flag = 0;
	if(paired==1){N=N/2; lambda=2*lambda;}


	//fprintf(stderr,"ll_length(bins)=%d\n",ll_length(bins));
	while(flag==0){
		do {
			nmerged = bic_seq_level(level); total_merged += nmerged*(level-1);
			//fprintf(stderr,"Merged %d bins of size %d; %d remaining.\n",
			//	nmerged, level, ll_length(bins));
			++level;
		} while (nmerged>0);

		level = level -1;
		if(level<=2&&nmerged==0) flag = 1; else level = level-1;
	}
	fprintf(stderr,"Merged %d bins; %d remaining\n",total_merged,ll_length(bins));
	return;
	
}

void bic_seq_auto(int nbins,double FP, int paired){ /*nbins is the initial number of bins, and FP is the expected number of type I errors allowed in merging process*/
	int nmerged = 0;
	int level = 2;
	//int flag = 0;
	double siglevel = 0.0001;
	
	if(nbins<=0){fprintf(stderr,"Error in bic_seq_auto: number of bins must be positive\n");exit(1);}
	siglevel = qchisq(1,FP/((double) nbins),0);
	//fprintf(stderr,"nbins=%d\nsiglevel = %g\n,N=%d\n",nbins,siglevel,(int) N);
	if(paired !=1) lambda = siglevel/log(N); else {N=N/2; lambda = 2*siglevel/log(N);}
	//lambda = siglevel/log(N);
	

	nmerged = bic_seq_level(level);
	//fprintf(stderr,"Merged %d bins of size %d; %d remaining.\n",nmerged, level, ll_length(bins));
	nmerged = merge_refine(nbins,FP,0.1);
	//fprintf(stderr,"Merged %d bins of size %d in the refine step; %d remaining\n",nmerged, level, ll_length(bins));
        nmerged = bic_seq_level(level);
        //fprintf(stderr,"Merged %d bins of size %d; %d remaining.\n",nmerged, level, ll_length(bins));	
	return;
}

SEG_PERMUTE bic_seq_perm(BIN_LIST bins_in, double p, double FP,int num_perm ,int auto_select) /*p is the probability that a read is assigned as a tumor read*/
{	SEG_PERMUTE segs = NULL;
	int i;
	if(num_perm<=0) return segs;

	segs = SEG_PERMUTE_create(num_perm);

	for(i=0;i<segs->size;i++){
		bins = random_assign(bins_in, p);
		if(auto_select==1){
			bic_seq_auto(ll_length(bins_in),FP,0);
			}else{
			bic_seq(0);
			}
		segs->bins_perm[i] = bins;
		}

	return segs;
}


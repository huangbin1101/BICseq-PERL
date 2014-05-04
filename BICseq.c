/*This version added the option to report the outlier position to a file*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <getopt.h>
#include "read.h"
#include "pos_cnt_lst.h"
#include "bin.h"
#include "gamma.h"

typedef struct SRM_binning_t{
	char *tumor_file, *normal_file,*output_file, *inbin_file, *binout_file;
	char *outlier_file;
	FILE *in_tmor, *in_nml, *output, *inbin, *outbin, *outlier;
	int bin_size;
	int win_size; /*the window size for removing the singular reads*/
	double quantile, multple;
	int autoselect_lambda; /*if or choose lambda automatically*/
	double FP;
	double tumor_freq; /*parameter for resampling*/
	int B; /*number of resamplings or permutations;*/
	int fdr; /*0 no FDR estimation, 1: estimate FDR*/
	int resampling; /*1 use resampling, 0 use permutation*/ 
	int paired; /*if the data is paired-end or not*/
	double lambda;

	int insert, sd; /*insert size and standard deviation of insert sizes*/
	
	} SRM_binning; /*the argument passed to the program*/


static int cmpdouble(const void *a, const void *b);

static int *aggregate(double *reads, int n_reads, int *num_pos);
/* aggregate the reads located at the same position together.
 * Returned value is a two column matrix looking like
 * 
 * <position>	<number of read at this position>
 *
 * The argument num_pos will record the size of this matrix (number of rows)
 *
 * reads: the read positions which should be ordered nondecreasingly
 * n_reads: total number of reads (length of the vector reads)
 * */
BIN_LIST  binning(int *tumor_1bp_bin, int n_tumor, int *normal_1bp_bin, int n_normal,int bin_size);
/* tumor_1bp_bin: output of the function aggregate for tumor reads
 * n_tumor: nrows of tumor_1bp_bin
 * similar for normal_1bp_bin and n_normal
 * bin_size: the bin size
 * nbins: total number of bins created. 
 *
 * */

BIN_LIST sort_rms_binning(double *tumor,int n_tmor,double *normal, int n_nml,int bin_size, int *nbins,int w, double quantile, double multple,FILE *outlier,char *tumor_file,char *normal_file);
/* sort the reads, remove the singular reads and bin the data.
   arguments:
	nbins: the total number of bins obtained.
	w: the window size used to determine the singular genomic positions.
	quantile: quantile (e.g. 0.95) used to determine the singular genomic positions.
	multple: if a genomic position has more than multple * quantile, it will be identified as an outlier.
 */

static void random_sample(double *tumor, int n_tumor, double *normal, int n_nml, int n,  int bin_size ,int *bin, int n_bins, int paired, int insert, int sd);
/* random sample n reads from tumor and normal reads and bin the sampled reads to bin_size bins,
   the binned data will be stored at *bin which should have enough space for storing the binned data
   paired: if the reads are paired data;
*/
static void bic_seq_resample(double *tumor, int n_tumor, double *normal, int n_nml, SRM_binning args);/*tumor and normal should have been ordered nondecreasingly*/

static void check_read(double *rds, int nrds); /*check if the reads are all nonnegative; rds should have been ordered nondecreasingly*/

static void explain_command(char **argv);
static SRM_binning option_assign(int argc, char **argv);

int main(int argc, char **argv)
{
	int n_tmor,n_nml,ncol,nbins;
	double  *tumor, *normal;
	SRM_binning args; 
	BIN_LIST bins = ll_new();
	
	args = option_assign(argc,argv);
	
	if(args.inbin_file==NULL){ /*input is not binned data*/
		tumor = read_table(args.in_tmor,&n_tmor,&ncol,-1,0);
		if(tumor==NULL) { fprintf(stderr,"Warning: the file %s is empty.\n",args.tumor_file);exit(1);}
		if(ncol!=1) {fprintf(stderr,"Error: tumor file has multiple columns.\n"); exit(1);}
		fprintf(stderr,"%d tumor reads loaded\n",n_tmor);

		normal = read_table(args.in_nml,&n_nml,&ncol,-1,0);
		if(normal==NULL)  { fprintf(stderr,"Warning: the file %s is empty.\n",args.normal_file);exit(1);}
		if(ncol!=1) { fprintf(stderr,"Error: normal file has multiple columns.\n"); exit(1);}
		fprintf(stderr,"%d normal reads loaded\n",n_nml);


		if(args.fdr==1&&args.resampling==1){
                        qsort(tumor,n_tmor,sizeof(double),cmpdouble);
                        fprintf(stderr,"sorted %d case reads\n",n_tmor);
                        qsort(normal,n_nml,sizeof(double),cmpdouble);
                        fprintf(stderr,"sorted %d control reads.\n",n_nml);
                        check_read(tumor,n_tmor);
                        check_read(normal,n_nml);

                        fprintf(stderr,"Start resampling\n");
                        bic_seq_resample(tumor,n_tmor,normal,n_nml,args);
                        free(tumor);tumor=NULL;
                        free(normal);normal = NULL;		
			}

		else{
			fprintf(stderr,"Binning\n");
			bins = sort_rms_binning(tumor,n_tmor,normal,n_nml,args.bin_size,&nbins,args.win_size,args.quantile,args.multple,args.outlier,args.tumor_file,args.normal_file);
			/*sort, remove singular positions and bin*/

			if(args.binout_file!=NULL){
				fprintf(stderr,"Reporting Binned data\n");
				BIN_LIST_print(bins, args.outbin);
				}

			if(args.fdr!=1){
				set_BinList(bins);
				fprintf(stderr,"Merging\n");
				if(args.autoselect_lambda!=1){
					bic_seq(args.paired);
					}else{
					bic_seq_auto(ll_length(bins),args.FP, args.paired);
					}
				bins = get_BinList();
				BIN_LIST_print(bins, args.output);
				ll_dealloc(bins);
				}else{
				SEG_PERMUTE segs = NULL;
				segs = bic_seq_perm(bins, args.tumor_freq, args.FP,args.B ,args.autoselect_lambda);
			        print_SEG_PERMUTE(segs,args.output);
			        SEG_PERMUTE_destroy(segs); segs = NULL;
				}
			}

		}else{
		double N_total=0.0, N_tumor=0.0,freq;
		int start, end, tumor, total;
	        while (fscanf(args.inbin, "%d %d %lf %d %d",
	                &tumor, &total, &freq, &start, &end) != EOF) {
	                ll_append(bins, bin_new(tumor, total, freq, start, end));
	                N_total += total; N_tumor += tumor;
                	}
		 set_totalreadcount(N_tumor,N_total-N_tumor);
		 
		if(args.fdr!=1){
			set_BinList(bins);
			fprintf(stderr,"Merging\n");
			if(args.autoselect_lambda!=1){
				bic_seq(args.paired);
				//bic_seq(0);
				}else{
				bic_seq_auto(ll_length(bins),args.FP, args.paired);
				//bic_seq_auto(ll_length(bins),args.FP, 0);
				}
			bins = get_BinList();
			BIN_LIST_print(bins, args.output);
			ll_dealloc(bins);
			}else{
			SEG_PERMUTE segs = NULL;
			segs = bic_seq_perm(bins, args.tumor_freq, args.FP,args.B ,args.autoselect_lambda);
			print_SEG_PERMUTE(segs,args.output);
			SEG_PERMUTE_destroy(segs); segs = NULL;
			}
		}

	return 0;	

}



static void explain_command(char **argv)
	{ fprintf(stderr,"Usage: %s [options] <case seq file> <control seq file>\n",argv[0]);
	  fprintf(stderr,"<case seq file>: a one column file only containing the positions of case reads.\n");
	  fprintf(stderr,"<control seq file>: a one column file only containing the positions of control reads.\n");
	  fprintf(stderr,"Options:\n");
	  fprintf(stderr,"        -h: print this message.\n");
	  fprintf(stderr,"        -l: the the penalty parameter lambda; default is 2.\n");
	  fprintf(stderr,"        -b <int>: the bin size; default is 100\n");
	  fprintf(stderr,"        -q <float>: the quantile used for identification of the singular genomic positions; default is 0.95\n");
	  fprintf(stderr,"        -w <int>: the window size for calculating the quantiles; default is 200\n");
	  fprintf(stderr,"        -o <string>: the segmentated output; if unspecified, print to the stdout.\n");
	  fprintf(stderr,"        -R <string>: Report the outlier position, count and the quantile in its local window to a file.\n");
	  fprintf(stderr,"        -f <float>: expected number of type I errors in merging process\n");
	  fprintf(stderr,"        -2: the data is paried-end data\n");
	  fprintf(stderr,"        -p <folat>: the overall frequency of tumor reads; if specified use resampling procedure for FDR esimate\n");
	  fprintf(stderr,"        -B <int>: how many resampling should be performed; Default is 0.\n");
	  fprintf(stderr,"        -I <InsertSize,SDofInsertSize>: if the data is paried-end, specify its mean insert size and standard deviation; Useful for resampling procedure; Default <200,20>\n");
	  fprintf(stderr,"        --multiplicity <float>: If a genomic position has more than multiplicity*quantile number of reads,\n");
	  fprintf(stderr,"                        it will be viewed as an outlier\n");
	  fprintf(stderr,"                        and the number of reads at this position will be set as multiplicity*quantile;\n");
	  fprintf(stderr,"                        default is 2.0\n");
	  fprintf(stderr,"       --binOut <string>: if specifed, print the binned data to this file\n");
	  fprintf(stderr,"       --inb <string>: the input is binned data; if this is specifed, <case seq file> and <control seq file> will be ignored (if also provided)\n");
	  fprintf(stderr,"       --perm: if specified,use permutation method for FDR estimate(default: use resampling method); Only valid if case/control reads are provided\n");
	  return;
	}


static SRM_binning option_assign(int argc, char **argv)
	{ SRM_binning args;
	  int c,option_index;
	  static struct option long_options[] =
			{ {"multiplicity", required_argument ,0,0},
			  {"binOut",required_argument,0,2},
			  {"inb",required_argument,0,3},
			  {"perm",no_argument,0,4},
			  {0,0,0,0} 
			};
	  char *insert=NULL;

	  args.multple = 2.0;
	  args.quantile = 0.95;
	  args.bin_size = 100;
	  args.win_size = 200;
	  args.output_file = NULL;
	  args.binout_file = NULL;
	  args.inbin_file = NULL;
	  args.outlier_file=NULL;
	  args.autoselect_lambda = 0;
	  args.FP = 1;
	  args.B = 10;
	  args.fdr = 0;
	  args.lambda = 2;
	  args.paired = 0;
	  args.tumor_file = NULL;
	  args.normal_file = NULL;
	  args.insert = 220;
	  args.sd = 20;
	  args.tumor_freq = 0.5;
	  args.resampling = 1;
	
	  args.in_tmor = NULL;
	  args.in_nml = NULL;
	  args.output = NULL;
	  args.inbin = NULL;
	  args.outbin = NULL;
	  args.outlier = NULL;

	  while((c=getopt_long(argc,argv,"2q:b:w:o:R:hf:p:B:l:",long_options,&option_index))!=-1){
		switch(c){
			case 0:
				args.multple = atoi(optarg);
				if(args.multple<1.0) { fprintf(stderr,"Error,multiplicity must be larger than or equal to 1.0\n"); exit(1);}
				break;
			case 2:
				args.binout_file = strdup(optarg);
				break;
			case 3: 
				args.inbin_file = strdup(optarg);
				break;
			case 4: args.resampling = 0;
				break;
			case 'l': 
				args.lambda = atof(optarg);
				if(args.lambda<=0.0) { fprintf(stderr,"Parameter misspecification: Lambda must be positive number.\n"); exit(1);}
				break;
                        case 'f':
                                args.FP = atof(optarg);
                                args.autoselect_lambda = 1;
                                if(args.FP<=0.0) {fprintf(stderr,"Expected number of type I error must be postive\n");exit(1);}
                                break;

			case 'I':
				insert = strdup(optarg);
				args.paired = 1;
				break;
			case 'q':
				args.quantile = atof(optarg);
				if(args.quantile<=0.0 || args.quantile>1.0)  { fprintf(stderr,"Quantile should be between 0 and 1.\n"); exit(1);}
				break;
			case 'b':
				args.bin_size = atoi(optarg);
				if(args.bin_size<1) {fprintf(stderr,"Bin size should be at least 1.\n");exit(1);}
				break;
			case 'w':
				args.win_size = atoi(optarg);
				if(args.win_size<1) {fprintf(stderr,"The window size must be positive.\n");exit(1);}
				break;
			case '2':
				args.paired = 1;
				break;
			case 'o':
				args.output_file = strdup(optarg);
				break;
			case 'R':
				args.outlier_file = strdup(optarg);
                        case 'p':
                                args.tumor_freq = atof(optarg);
                                if(args.tumor_freq<0.0||args.tumor_freq>=1.0){fprintf(stderr,"The probabilty of a read being a tumor read must be between 0 and 1\n");exit(1);}
                                args.fdr = 1;
                                break;
                        case 'B':
                                args.B = atoi(optarg);
                                if(args.B<=0) {fprintf(stderr,"The option -B must be a positve number\n");exit(1);}
                                break;
			case 'h': 
				explain_command(argv);
				exit(0);
			case '?': /* getopt_long already printed an error message. */
				exit(1);
				break;
			default:
				abort (); 
			}
		}

        if (argc - optind!=2&&args.inbin_file==NULL){
           	explain_command(argv); exit(1);
         	}

	if(args.inbin_file==NULL){
		args.tumor_file = strdup(argv[optind]);
		args.normal_file = strdup(argv[optind+1]);

	        args.in_tmor = fopen(args.tumor_file,"r");
	        if(args.in_tmor == NULL)
	                { fprintf(stderr,"fopen %s: %s\n", args.tumor_file,strerror(errno));
        	          exit(1);
	                }

	        args.in_nml = fopen(args.normal_file,"r");
	        if(args.in_nml == NULL)
        	        { fprintf(stderr,"fopen %s: %s\n", args.normal_file,strerror(errno));
	                  exit(1);
	                }

		}
	else{
		args.inbin =  fopen(args.inbin_file,"r");
		if(args.inbin==NULL){
			fprintf(stderr,"fopen %s: %s\n", args.inbin_file,strerror(errno));
                        exit(1);
			}
		}
	

        if(args.outlier_file!=NULL) {
                args.outlier = fopen(args.outlier_file,"w");
                if(args.outlier==NULL) {fprintf(stderr,"fopen %s: %s\n",args.outlier_file,strerror(errno));exit(1);}
                }

        if(args.output_file!=NULL){
                args.output = fopen(args.output_file,"w");
                if(args.output==NULL) { fprintf(stderr,"fopen %s: cannot create the file.\n",args.output_file); exit(1);}
                }else{
		args.output = stdout;
		}

        if(args.binout_file!=NULL){
                args.outbin = fopen(args.binout_file,"w");
                if(args.outbin==NULL) { fprintf(stderr,"fopen %s: cannot create the file.\n",args.binout_file); exit(1);}
                }

	if(insert!=NULL){
		char *substr;
		substr = strtok(insert,",");
		args.insert = atoi(substr);
		substr=strtok(NULL,",");
		if(substr==NULL){fprintf(stderr,"Incorrect format for option -I\n");exit(1);}
		args.sd = atoi(substr);
		if(args.insert<0||args.sd<0) {fprintf(stderr,"Error, insert size and its standard deviaiton must be positive\n");}
		}

	set_lambda(args.lambda);

	return args;

	}


static int cmpdouble(const void *a, const void *b)
	{ double tmp1,tmp2;
	  tmp1 = *((const double *) a);
	  tmp2 = *((const double *) b);
	  if(tmp1<tmp2) return -1;
	  else if(tmp1>tmp2) return 1;
	  else return 0;
	}




static int *aggregate(double *reads, int n_reads, int *num_pos)
	{ /*The first column of read_dist is the position of the read*/
	  /*The second column of read_list is the number of reads at this position.*/
	  /*If certain position does not have any reads, this position will be ignored*/
	  int *reads_dist, pos_min, pos_max,size;
	  int i,j;
	
	  pos_min = (int) reads[0];
	  pos_max = (int) reads[n_reads-1];

	  size = (n_reads<pos_max-pos_min+1)? n_reads : pos_max-pos_min+1;
	  reads_dist = (int *)  malloc(sizeof(int)*(size*2+10));

	  j = 0;
	  i =0;
	  reads_dist[2*j] = (int) reads[i]; /*The first column of read_dist is the position of the read*/
	  reads_dist[2*j+1] = 1; /*The second column of read_list is the number of reads at this position.*/

	  for(i=1;i<n_reads;i++)
		{ if((int) reads[i]==reads_dist[2*j]) 
			{ reads_dist[2*j+1]++;} /*one more read at position reads_list[2*j]*/
		  else /*get to a new position*/
			{ j++;
			  reads_dist[2*j] = (int) reads[i];
			  reads_dist[2*j+1] = 1;
			}
		}

	 *num_pos = j+1;
	
	 return reads_dist;

	}


BIN_LIST binning(int *tumor_1bp_bin, int n_tumor, int *normal_1bp_bin, int n_normal,int bin_size)
	{ int i,nrow,max_pos,min_pos;
	  int start,end ,i_tum, i_norm;
	  double tumorcnt, normalcnt,freq,total, N_tumor, N_normal;
	  BIN_LIST bins = ll_new();

	  max_pos = (tumor_1bp_bin[2*(n_tumor-1)]>normal_1bp_bin[2*(n_normal-1)])? tumor_1bp_bin[2*(n_tumor-1)]: normal_1bp_bin[2*(n_normal-1)];
	  min_pos = (tumor_1bp_bin[0]<normal_1bp_bin[0])?tumor_1bp_bin[0]:normal_1bp_bin[0];
	  nrow = (max_pos-min_pos+1)/bin_size+10;

	  N_tumor = 0.0; N_normal = 0.0;
	  i_tum = 0;
	  i_norm = 0;
	  start = min_pos - (min_pos-1)%bin_size;/*the start position of the left most bin*/
	  end = start + bin_size -1;
	  for(i=0;i<nrow;i++)
		{ tumorcnt = 0.0;
		  normalcnt = 0.0;
		  while(tumor_1bp_bin[2*i_tum]<=end&&i_tum<n_tumor)
			{ tumorcnt += tumor_1bp_bin[2*i_tum+1];
			  i_tum++;
			}
		  while(normal_1bp_bin[2*i_norm]<=end&&i_norm<n_normal)
			{ normalcnt += normal_1bp_bin[2*i_norm+1];
			  i_norm++;
			}

		  total = tumorcnt + normalcnt;
		  freq = tumorcnt/total;
		  if(total>0.0) ll_append(bins, bin_new(tumorcnt, total, freq, start, end));
		  N_tumor += tumorcnt;
		  N_normal += normalcnt;

		  start += bin_size;
		  end = start + bin_size -1;
		}

	   set_totalreadcount(N_tumor,N_normal);
	  return bins;
	}


BIN_LIST sort_rms_binning(double *tumor, int n_tmor,double *normal, int n_nml,int bin_size, int *nbins,int w, double quantile, double multple,FILE *outlier,char *tumor_file,char *normal_file)
	{ int *tumor_1bpbin,num_tum_1bp, *normal_1bpbin,num_norm_1bp;
	  BIN_LIST bins;

	  qsort(tumor,n_tmor,sizeof(double),cmpdouble);
          fprintf(stderr,"sorted %d case reads\n",n_tmor);
       	  qsort(normal,n_nml,sizeof(double),cmpdouble);
      	  fprintf(stderr,"sorted %d control reads.\n",n_nml);
 
	  /*aggregate the reads; */
	  tumor_1bpbin = aggregate(tumor,n_tmor,&num_tum_1bp);
	  free(tumor); /*If invoke in R, it's maybe better not free tumor here*/
	  normal_1bpbin = aggregate(normal,n_nml,&num_norm_1bp);
	  free(normal);	/*If invoke in R, it's maybe better not free normal here*/

	  /*remove the singular points*/
	  if(outlier!=NULL) {fprintf(outlier,"#case\n#%s\npos\tcount\tquantile_in_windonw\n",tumor_file);}
	  singularity_rm(tumor_1bpbin, num_tum_1bp, w, quantile, multple, outlier);
	  if(outlier!=NULL) fprintf(outlier,"#control\n#%s\npos\tcount\tquantile_in_windonw\n",normal_file);
	  singularity_rm(normal_1bpbin, num_norm_1bp, w, quantile, multple,outlier);

	  /*bin the processed reads*/
	  bins = binning(tumor_1bpbin,num_tum_1bp,normal_1bpbin, num_norm_1bp,bin_size);
	  free(tumor_1bpbin);
	  free(normal_1bpbin);

	 return bins;
	}

/*tumor and normal should have been sorted before call the following function*/
static void random_sample(double *tumor, int n_tumor, double *normal, int n_nml, int n,  int bin_size ,int *bin, int n_bins, int paired, int insert, int sd)
{	int i,k,pos;
	double U, rn_insert, U1;
	//int n_sampled=0, t_sampled=0;
	int N = n_tumor + n_nml, max_pos, min_pos;
	double tmp, ave_dist, offset;

	tmp = (tumor[n_tumor-1]>normal[n_nml-1]? tumor[n_tumor-1]: normal[n_nml-1]);
	max_pos = (int) tmp;
	if(max_pos/bin_size > n_bins-1) {fprintf(stderr,"Error in random_sample: not enough space for binned data\n");exit(1);}
        tmp = (tumor[0]<normal[0]? tumor[0]: normal[0]);
        min_pos = (int) tmp;
	
	ave_dist = ((double) max_pos - min_pos)/((double) n_nml+n_tumor) + 1;  /*plus 1 to make sure ave_dist at least 2*/
	//ave_dist = 100;

	/*fist initialize n_bins*/
	for(i=0;i<n_bins;i++){
		bin[i] = 0;
		}

	/*random sample and binning*/
	if(paired!=1){
		int pos1, pos2;
		for(i=0;i<n;i++){
			//U = N*rand_lp();
			//U1 = 2*rand_lp()-1;
			U = N*drand48();
			U1 = 2*drand48() -1;
			offset = U1*ave_dist;
			if(U > n_tumor){
				k = ((int) U - n_tumor)%n_nml;
				pos = (int) normal[k];
				pos = pos > 0 ? pos : 1;
				pos = pos < max_pos ? pos : max_pos;
				//bin[pos/bin_size] ++;
				//n_sampled++;
				}else{
				k = ((int) U)%n_tumor;
				pos = (int) tumor[k];
                                pos = pos > 0 ? pos : 1;
                                pos = pos < max_pos ? pos : max_pos;
				//bin[pos/bin_size] ++;
				//t_sampled++;
				}
			pos1 = pos; 
			pos = pos + offset;
                        pos = pos>0 ? pos:1;
                        pos = pos < max_pos ? pos : max_pos;
                        //bin[pos/bin_size] ++;
                        pos2 = pos;
			if(drand48()<0.5){ /*add some noise into the simulated data*/
	                        if(drand48()<0.95) bin[pos1/bin_size] ++;
	                        if(drand48()<0.05) bin[pos2/bin_size] ++;						
				}else{
                                if(drand48()<0.05) bin[pos1/bin_size] ++;
                                if(drand48()<0.95) bin[pos2/bin_size] ++;
				}
			}
		}else{
		//int pos1, pos2;
                for(i=0;i<n/2;i++){
                        U = N*rand_lp();
                        if(U > n_tumor){
                                k = ((int) U - n_tumor)%n_nml;
                                pos = (int) normal[k];
                                bin[pos/bin_size] ++;
				//n_sampled++;
                                }else{
                                k = ((int) U)%n_tumor;
                                pos = (int) tumor[k];
                                bin[pos/bin_size] ++;
				//t_sampled++;
                                }
			//pos1 = pos;
			rn_insert = rnorm((double)insert, (double) sd);
			pos = pos+rn_insert;
			pos = pos>0 ? pos:1;
			pos = pos < max_pos ? pos : max_pos;
			bin[pos/bin_size] ++;
			//pos2 = pos;
			//if(rand_lp()<0.5) bin[pos1/bin_size] --;
			//if(rand_lp()<0.5) bin[pos2/bin_size] --;
                        }

		}
	//fprintf(stderr,"sampled %d from normal\n",n_sampled);
	//fprintf(stderr,"sampled %d from tumor\n", t_sampled);
	return;
}


static void check_read(double *rds, int nrds)
{	if(rds[0] < 0.0) {fprintf(stderr,"Error, read position must be nonnegative\n"); exit(1);}
}


static void  bic_seq_resample(double *tumor, int n_tumor, double *normal, int n_nml, SRM_binning args)
{	SEG_PERMUTE segs = NULL;
	int *tumor_bin, *normal_bin, nbins;
	int n_tumor_sample, n_normal_sample,i,k, total,start,end, kmin;
	double tmp, freq, N_tumor, N_normal;
        struct timeval tv;
        int seed;

        gettimeofday(&tv, NULL);
        seed = tv.tv_sec * 1000000 + tv.tv_usec;
        seed_set(seed);
	srand48(seed);
	
	segs = SEG_PERMUTE_create(args.B);

	tmp = tumor[n_tumor-1] > normal[n_nml-1] ? tumor[n_tumor-1]:normal[n_nml-1];
	nbins = floor(tmp/args.bin_size)+10;
	nbins = nbins>10?nbins:10;
	tumor_bin = (int *) malloc(sizeof(int)*nbins);
	normal_bin = (int *)malloc(sizeof(int)*nbins);
	if(tumor_bin==NULL||normal_bin==NULL){
		fprintf(stderr,"Error in bic_seq_resample: memory allocation failed\n");
		exit(1);
		}

        tmp = tumor[0] < normal[0] ? tumor[0]:normal[0];
        kmin = (int) floor(tmp/args.bin_size)-1;
        kmin = (kmin>0? kmin:0);

	for(i=0;i<segs->size;i++){
		n_tumor_sample = rbinom(args.tumor_freq,n_tumor+n_nml);
		n_normal_sample = rbinom(1-args.tumor_freq,n_tumor+n_nml);
		random_sample(tumor, n_tumor, normal, n_nml, n_tumor_sample,  args.bin_size ,tumor_bin, nbins, args.paired, args.insert, args.sd);
		random_sample(tumor, n_tumor, normal, n_nml, n_normal_sample, args.bin_size ,normal_bin,nbins, args.paired, args.insert, args.sd);


		N_tumor=0.0; N_normal = 0.0;
		for(k=kmin;k<nbins;k++){
			start = k*args.bin_size+1;
			end = start+args.bin_size;
			total = tumor_bin[k] + normal_bin[k];
			freq = ((double) tumor_bin[k])/((double) total);
			if(total>0) ll_append(segs->bins_perm[i], bin_new(tumor_bin[k], total, freq, start, end));
			N_tumor += tumor_bin[k];
			N_normal += normal_bin[k];
			}
		set_BinList(segs->bins_perm[i]);
		set_totalreadcount(N_tumor,N_normal);

                if(args.autoselect_lambda!=1){
                        bic_seq(args.paired);
			//bic_seq(0);
                        }else{
                        bic_seq_auto(ll_length(segs->bins_perm[i]),args.FP,args.paired);
			//bic_seq_auto(ll_length(segs->bins_perm[i]),args.FP,0);
                        }
		segs->bins_perm[i] = get_BinList();
		}

	print_SEG_PERMUTE(segs,args.output);
	SEG_PERMUTE_destroy(segs); segs = NULL;
	free(tumor_bin); tumor_bin = NULL;
	free(normal_bin);normal_bin = NULL;

	return;
}

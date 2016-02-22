/* 

The code itself is Copyright (C) 2016, by David J. Cutler.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h> 
#include <ctype.h>
#include <zlib.h>
#include <time.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>


#define FALSE 0
#define TRUE 1

#define SQR(a) ((a)*(a))
#define NO_ALLELES 6
#define MAX_GENOTYPES 14


/* row 0   =  A
   1   =  C
   2   =  G
   3   =  T
   4   =  Del
   5   =  Ins
   6   =  M  AC
   7   =  R  AG
   8   =  W  AT
   9   =  S  CG
   10  =  Y  CT
   11  =  K  GT
   12  = Del Het
   13  = Ins Het
*/
typedef struct sample_node
{
  double post_prob[MAX_GENOTYPES+1];
  double final_p;
  double coef;
  int reads[NO_ALLELES];
  double frac[NO_ALLELES];
  int tot;
  char final_call;
  char initial_call;
  char family[64];
  char indiv[64];
  short sex;
  struct sample_node *mom;
  struct sample_node *dad;
  int which;
} SAMNODE;

typedef struct config_node
{
  int genotype_count[MAX_GENOTYPES];
  double prior;
  double like;
  double post;
  char *sample_calls;
  float avg_depth;
  short no_alleles;
  int allele_count[NO_ALLELES];
  int no_denovo;
  int hets;
  int homs;
} CNODE;

short dyad_denovo[MAX_GENOTYPES+1][MAX_GENOTYPES+1];
short trio_denovo[MAX_GENOTYPES+1][MAX_GENOTYPES+1][MAX_GENOTYPES+1];
short allele_counts[4][MAX_GENOTYPES][NO_ALLELES];

#define LOG10 (double)2.302585

#define minim(atesta,btestb) ((atesta<btestb)?atesta:btestb)
#define maxim(atesta,btestb) ((atesta>btestb)?atesta:btestb)

#define AUTO 0
#define CHRX 1
#define CHRY 2
#define CHRMT 3

#define REF_ALLELE 0
#define SNP 1
#define DELETION 2
#define INSERTION 3
// #define LOW 4
// #define HIGH 5
#define MESS 4

#define DATA_EMPTY 0
#define DATA_LOADED 1
#define DATA_RUNNING 2
#define DATA_ALL_DONE 3 

typedef struct pthread_data_node 
{
  pthread_mutex_t mutex;
  int tid;
  int status;
  int chrom;
  int pos;
  long tot_pos;
  char dom;
  int dom_int;
  char fragment[1024];
  SAMNODE **samples;
  int HAPLOID;
} PTHREAD_DATA_NODE;

int no_threads;
pthread_mutex_t outfile_write_mutex;
pthread_mutex_t snpfile_write_mutex;
pthread_t *threads;
double ln_denovo;
double ***ln_HW;
double config_threshold = 2.3;
int use_ped;

int **genotype_order;
int INDIV;


int gen_to_int(char c);
char int_to_gen(int c);
int clean_config_probs(CNODE **cn, SAMNODE **sn, int n, int max,int max_gen,int **alpha,double *coef,int indiv,int depth,int ref,int HAPLOID);
int sort_configs(const void *a,const void *b);
int fill_config_probs(CNODE **cn, int n, SAMNODE *sn, int max,int **alpha,double *coef,int indiv,int this_depth,int ref,int chrom,int HAPLOID);
SAMNODE *sample_alloc();
CNODE *config_alloc(int N);
void config_free(CNODE *tn,int N);
void fill_hardy_weinberg(double **exact_HW,int asize,int n);
void fill_alpha_prior(int **alpha,int max_gen,int hom,int het,int ref);
void fill_alpha_coef(int **alpha,double *coef,int max_gen);
void fill_hom_like(CNODE *cn, SAMNODE **sn, int n,int **alpha,double *coef);
int add_denovo(int kid, int dad, int mom,int sex,int chrom);
double gammln(double xx);
double exactfactln(int n);
double factln(int n);
void check_alpha_sanity(int **alpha_prior,int max_gen,int **first_alpha_prior,int ref,double **weight);
void get_het_alleles(int i, int *a, int *b,int ref);
int get_snp_type(SAMNODE *sm,int ref,char *minor);
void read_var(char *line,char *result);
int *ivector(int nl,int nh);
unsigned int *uvector(int nl,int nh);
char *cvector(int nl,int nh);
double *dvector(int nl, int nh);
double **dmatrix(int nrl,int nrh,int ncl,int nch);
unsigned short **usmatrix(int nrl,int nrh,int ncl,int nch);
int **imatrix(int nrl,int nrh,int ncl,int nch);
char **cmatrix(int nrl,int nrh,int ncl,int nch);
void free_cvector(char *v, int nl, int nh);
void free_ivector(int *v, int nl, int nh);
void free_uvector(unsigned int *v, int nl, int nh);
void free_dvector(double *v, int nl,int nh);
void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch);
unsigned int **umatrix(int nrl,int nrh,int ncl,int nch);
void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch);
void free_ucmatrix(unsigned char **m,int nrl,int nrh,int ncl,int nch);
void free_cmatrix(char **m,int nrl,int nrh,int ncl,int nch);
unsigned char **ucmatrix(int nrl,int nrh,int ncl,int nch);
void dump_error(char *error_text);
double chidist(double x, double lambda, int df);
double gammq(double a, double x);
double gammp(double a,double x);
void gser(double *gamser,double a,double x,double *gln);
void gcf(double *gammcf,double a,double x,double *gln);
double gammln(double xx);
double brent(double ax, double bx, double cx, double (*f)(double), double tol,double *xmin);
double calc_chi(double gamma);
int sort_compare(const void *a,const void *b);
double exactfactln(int n);
double factln(int n);
int find_chrom(unsigned int *pos,int first, int last,int try,unsigned this);
unsigned int find_lowest(unsigned *list,int n);
void init_genome_buffer(gzFile *mfile);
void free_genome_buffer(void);
char get_genome(unsigned int which,gzFile *mfile);
void *call_single_base(void *threadid);

int max_gen,min_depth_needed;
int ALL_FINISHED = FALSE;

double EPS = 1e-8;
#define MAX_DIST 501
gzFile *outfile,*pilefile;
double ***HW_exact;
double THRESHOLD,theta,ln_theta;
char snp_type[MESS+1][80];
FILE *snpfile;

int main(int argc, char *argv[])
{
  FILE *sfile,*guide_file;
  gzFile *reffile;
  char ss[4096],sdxname[4096],sss[4096];
  char **filenames;
  DIR *pDIR;
  struct dirent *pDirEnt;
  int no_files,i,j,k;
  int max_contigs;
  int use_guide = FALSE;
  int which;
  unsigned int *frag_pos;
  char **contig_names;
  int no_contigs;
  unsigned int tot_bases = 0;
  int HAPLOID;
  char linkage_filename[512];
  int *chrom_type;
  double denovo_rate;


  FILE *distfile;
  gzFile **pileupfile;  

  if(argc < 10 || argc > 13)
  {
      printf("\n Usage %s pileup_extension sdx_file no_files outfile Prob_to_call theta haploid[y,n] no_threads use_pedfile[y,n] [pedfilename] [denovo_mutation_rate] [guide_file_bed_format]\n",argv[0]);
      exit(1);
  }

  no_threads = (int)atoi(argv[8]);
  if(no_threads < 2 || no_threads > 200)
  {
      printf("\n Number of threads is limited to 2 to 200.   You entered %d \n\n",no_threads);
      exit(1);
  }
  no_threads--;

  sprintf(ss,"%s.base.gz",argv[4]);
  if((outfile=gzopen(ss,"w"))==(gzFile *)NULL)
  {
      printf("\n Can not open file %s\n",ss);
      exit(1);
  }

  sprintf(ss,"%s.snp",argv[4]);
  if((snpfile=fopen(ss,"w"))==(FILE *)NULL)
  {
      printf("\n Can not open file %s for writing\n",ss);
      exit(1);
  }

  sprintf(ss,"%s.dist",argv[4]);
  if((distfile=fopen(ss,"w"))==(FILE *)NULL)
  {
      printf("\n Can not open file %s for writing\n",ss);
      exit(1);
  }

  sprintf(ss,"%s.piles.gz",argv[4]);
  if((pilefile=gzopen(ss,"w"))==(gzFile *)NULL)
  {
      printf("\n Can not open file %s for writing\n",ss);
      exit(1);
  }

  guide_file = stdout;


  THRESHOLD = (double)atof(argv[5]);
  theta = (double)atof(argv[6]);
  if( (theta < 1e-10) || (theta > 0.5) )
    {
      printf("\n Encountered impossible value for theta = %g \n",theta);
      exit(1);
    }
  ln_theta = log(theta);

  HAPLOID = FALSE;
  max_gen = min_depth_needed = 0;

  for(i=0;i<=MAX_GENOTYPES;i++)
    for(j=0;j<=MAX_GENOTYPES;j++)
    {
	dyad_denovo[i][j] = 0;
	for(k=0;k<=MAX_GENOTYPES;k++)
	  trio_denovo[i][j][k] = 0;
    }
  strcpy(ss,argv[7]);
  if( (strchr(ss,'Y')) || (strchr(ss,'y')) )
  {
      HAPLOID = TRUE;
      max_gen = 6;
      min_depth_needed = 1;
      for(i=0;i<max_gen;i++)
	for(j=0;j<max_gen;j++)
	  if(i!=j)
	    dyad_denovo[i][j] = 1;
  }
  else
  {
      max_gen = MAX_GENOTYPES;
      min_depth_needed = 2;
      int da, db, ma, mb, ka, kb;
      for(i=0;i<max_gen;i++)
      {
	  get_het_alleles(i,&da,&db,0);
	  for(j=0;j<max_gen;j++)
	  {
	      get_het_alleles(j,&ka,&kb,0);
	      if((ka != da) && (ka != db) && (kb != da) && (kb != db))
		dyad_denovo[i][j] = 1;
	      for(k=0;k<max_gen;k++)
	      {
		  get_het_alleles(k,&ma,&mb,0);
		  if( ((ka == ma) && (kb == da)) || ((ka == ma) && (kb == db)) || ((ka == mb) && (kb == da)) || ((ka == mb) && (kb == db)) ||
		      ((kb == ma) && (ka == da)) || ((kb == ma) && (ka == db)) || ((kb == mb) && (ka == da)) || ((kb == mb) && (ka == db)))
		    trio_denovo[i][k][j] = 0;
		  else
		    if( ((ka != ma) && (kb != db)) && ((kb != ma) && (ka != db)) && ((ka != mb) && (kb != da)) && ((kb != mb) && (ka != da)))
		      trio_denovo[i][k][j] = 2;
		    else
		      trio_denovo[i][k][j] = 1;
	      }
	  }
      }
  }
  strcpy(ss,argv[9]);
  if( (strchr(ss,'Y')) || (strchr(ss,'y')) )
  {
      strcpy(linkage_filename,argv[10]);
      strcpy(ss,argv[11]);
      denovo_rate = (double)atof(ss);
      if( (denovo_rate < 1e-30) || (denovo_rate > theta) )
	{
	  printf("\n Encounted impossible denovo mutation rate of %g with a theta of %g",denovo_rate,theta);
	  exit(1);
	}
      ln_denovo = log(denovo_rate);
      use_ped = TRUE;
      if(argc == 13)
      {
	  use_guide = TRUE;
	  sprintf(ss,"%s",argv[12]);
	  if((guide_file=fopen(ss,"r"))==(FILE *)NULL)
	  {
	      printf("\n Can not open file %s for writing which should contain the guide_file\n",ss);
	      exit(1);
	  }
      }

  }
  else
  {
      use_ped = FALSE;
      ln_denovo = 0;
      if(argc == 11)
      {
	  use_guide = TRUE;
	  sprintf(ss,"%s",argv[10]);
	  if((guide_file=fopen(ss,"r"))==(FILE *)NULL)
	  {
	      printf("\n Can not open file %s for writing which should contain the guide_file\n",ss);
	      exit(1);
	  }
      }
  }
  

  pDIR = opendir(".");

  if ( pDIR == NULL ) 
  {
      fprintf(stderr, "%s %d: opendir() failed (%s)\n",
	      __FILE__, __LINE__, strerror( errno ));
      exit( -1 );
  }

  strcpy(sdxname,argv[2]);
  if((sfile=fopen(sdxname,"r"))==(FILE *)NULL)
  {
      printf("\n Can not open file %s\n",sdxname);
      exit(1);
  }
  
  if(strstr(sdxname,".sdx") != NULL)
  {
      for(i=strlen(sdxname)-1;i>0;i--)
	if(sdxname[i] == '.')
	  {
	    sdxname[i] = '\0';
	    i = 0;
	  }
  }

  sprintf(sss,"%s.seq",sdxname);

  if((reffile=gzopen(sss,"r"))==(gzFile *)NULL)
  {
      printf("\n Can not open file %s for reading\n",sss);
      exit(1);
  }
  printf("\n About to initialize the genome buffer \n\n");
  init_genome_buffer(reffile);

  printf("\n Finished genome buffer initialization \n\n");
  fgets(sss,256,sfile);
  no_contigs = atoi(sss);
  max_contigs = no_contigs;
  frag_pos = uvector(-1,no_contigs);
  contig_names = cmatrix(0,no_contigs,0,256);
  chrom_type  = ivector(0,no_contigs);
  frag_pos[-1] = 0;
  char *token;
  for(i=0;i<no_contigs;i++)
  {
      // printf("\n About to read line %d \n\n",i);
      fgets(sss,1024,sfile);
      token = strtok(sss,"\t \n");
      frag_pos[i] = atoi(token) + 15;
      /* sprintf(sss,"%u\t>%s",&offsets[i],s2); */
      frag_pos[i] += frag_pos[i-1];
      token = strtok(NULL,"\t \n");
      strcpy(contig_names[i],token);
      strcpy(sss,contig_names[i]);
      chrom_type[i] = AUTO;
      char *token;
      token = strtok(sss,":_- \n\0");
      for(j=0;j<strlen(token);j++)
	token[j] = tolower(token[j]);
      if(strcmp(token,"chrx") == 0)
	chrom_type[i] = CHRX;
      else
	if(strcmp(token,"chry") == 0)
	  chrom_type[i] = CHRY;
	else
	  if(strcmp(token,"chrmt") == 0)
	    chrom_type[i] = CHRMT;
      // printf("\nFor contig %d we have offset %u\n\n",i,frag_pos[i]); 
  }
  fclose(sfile);

  printf("\n Finished reading the sdx file \n\n");

  no_files = atoi(argv[3]);  
  filenames = cmatrix(0,no_files,0,128);
  pileupfile = malloc(sizeof(gzFile *)*(no_files+1));
  if(!pileupfile)
    dump_error("\n Error allocating pileupfile array \n");
  
  pDirEnt = readdir( pDIR );
  printf("\n Just read dir \n\n");
  i = 0;
  while (pDirEnt != NULL && i <= no_files ) 
  {
	if( (strstr(pDirEnt->d_name,argv[1]) != NULL) )
	{
    	    if((pileupfile[i]=gzopen(pDirEnt->d_name,"rb"))==(gzFile *)NULL)
	    {
		printf("\n Can not open file %s which should contain pedigree information\n",pDirEnt->d_name);
		exit(1);
	    }
	    strcpy(sss,pDirEnt->d_name);
	    token = strtok(sss,"\n.\t \0");
	    strcpy(filenames[i++],token);
	}
	
	pDirEnt = readdir(pDIR);
  }
  closedir(pDIR);

  if(i > no_files)
    dump_error("\n Found more files than you specified \n");

  no_files = i;
  
  INDIV = no_files;

  sprintf(snp_type[SNP],"SNP");
  sprintf(snp_type[DELETION],"DEL");
  sprintf(snp_type[INSERTION],"INS");
  sprintf(snp_type[MESS],"MESS");



  printf("\n Found a total of %d individuals\n\n",INDIV);
  SAMNODE **samples;

  samples = (SAMNODE **)malloc((unsigned)((INDIV+1)*sizeof(SAMNODE *)));
  if(!samples) dump_error("Allocation failure in samples\n");

  int *sorted_samples;
  sorted_samples = ivector(0,INDIV);

  for(i=0;i<INDIV;i++)
  {
      samples[i] = sample_alloc();
      samples[i]->which = i;
      strcpy(samples[i]->indiv,filenames[i]);
      sorted_samples[i] = i;
  }

  if(use_ped)
  {
      FILE *pedfile;
      printf("\n Reading Ped file \n");
      if( (pedfile = fopen(linkage_filename,"r")) == NULL)
      {
		printf("\n Could Not open %s",linkage_filename);
		exit(1);
      }
      fgets(sss,81919,pedfile);	    
      while(!feof(pedfile) && strlen(sss) > 5)
      {
	  char *fam,*ind,*token;
	  fam = strtok(sss,"\n\t ");
	  ind = strtok(NULL,"\n\t ");
	  // printf("\n Read Pedigree data for fam=%s ind=%s",fam,ind);
	  for(i=0;i<INDIV;i++)
	    if(strcmp(ind,samples[i]->indiv) == 0)
	    {
		printf("\n\t\tFOUND %s ind",ind);
		strcpy(samples[i]->family,fam);
		token = strtok(NULL,"\n\t ");
		if(strcmp(token,"0") != 0)
		  for(j=0;j<INDIV;j++)
		    if(strcmp(token,samples[j]->indiv) == 0)
		      {
			samples[i]->dad = samples[j];
			j = INDIV;
		      }
		token = strtok(NULL,"\n\t ");
		if(strcmp(token,"0") != 0)
		  for(j=0;j<INDIV;j++)
		    if(strcmp(token,samples[j]->indiv) == 0)
		      {
			samples[i]->mom = samples[j];
			j = INDIV;
		      }
		token = strtok(NULL,"\n\t ");
		samples[i]->sex = (short)atoi(token);
	    }
	  sss[0] = '\0';
	  if(!feof(pedfile))
	    fgets(sss,81919,pedfile);
      }
      fclose(pedfile);
      printf("\n Done reading pedfile \n\n");
  }

  printf("\n Before Sorting Samples \n\n");
  //for(i=0;i<INDIV;i++)
  //{
  //	  printf("\n %d is sorted to %d which is .%s. .%s.",i,sorted_samples[i],samples[sorted_samples[i]]->family,samples[sorted_samples[i]]->indiv);
  // }
	
  for(i=0;i<INDIV;i++)
  {
      int k = 0;
      SAMNODE *bigger = NULL;
      if(samples[i]->dad)
      {
	  if(samples[i]->mom)
	  {
	      if(samples[i]->dad->which >= samples[i]->mom->which)
	      {
		bigger = samples[i]->dad;
		k = samples[i]->dad->which;
	      }
	      else
	      {
		bigger = samples[i]->mom;
		k = samples[i]->mom->which;
              }
          }
	  else
	  {
		bigger = samples[i]->dad;
		k = samples[i]->dad->which;
	  }
      }
      else
	if(samples[i]->mom)
        {
		bigger = samples[i]->mom;
		k = samples[i]->mom->which;
        }
      if(k > samples[i]->which)
      {
		int j = samples[i]->which;
		samples[i]->which = k;
		bigger->which = j;
		i = -1;
      }
	
  }
  printf("\n Done Sorting Samples \n\n");
  for(i=0;i<INDIV;i++)
  {
	  sorted_samples[i] = samples[i]->which;
  	  printf("\n %d is sorted to %d which is .%s. .%s. and labeled as number %d",i,sorted_samples[i],samples[i]->family,samples[i]->indiv,samples[i]->which);
  }

  HW_exact = NULL;
  int asize;
  if(!HAPLOID)
  {
      HW_exact = (double ***)malloc((unsigned) (INDIV+1)*sizeof(double **));
      if(!HW_exact) dump_error("Allocation failure in configs\n");
      for(i=1;i<=INDIV;i++)
      {
	  asize = 2*(i);
	  HW_exact[i] = dmatrix(0,asize,0,i);
	  fill_hardy_weinberg(HW_exact[i],asize,i);
      }
  }
  printf("\n Done filling the HW matrix \n\n");

  genotype_order = imatrix(0,NO_ALLELES-1,0,13);
  genotype_order[0][0] = 0;
  genotype_order[0][1] = 7;
  genotype_order[0][2] = 6;
  genotype_order[0][3] = 8;
  genotype_order[0][4] = 12;
  genotype_order[0][5] = 13;
  genotype_order[0][6] = 1;
  genotype_order[0][7] = 2;
  genotype_order[0][8] = 3;
  genotype_order[0][9] = 4;
  genotype_order[0][10] = 5;
  genotype_order[0][11] = 9;
  genotype_order[0][12] = 10;
  genotype_order[0][13] = 11;
  
  genotype_order[1][0] = 1;
  genotype_order[1][1] = 10;
  genotype_order[1][2] = 6;
  genotype_order[1][3] = 9;
  genotype_order[1][4] = 12;
  genotype_order[1][5] = 13;
  genotype_order[1][6] = 0;
  genotype_order[1][7] = 2;
  genotype_order[1][8] = 3;
  genotype_order[1][9] = 4;
  genotype_order[1][10] = 5;
  genotype_order[1][11] = 7;
  genotype_order[1][12] = 8;
  genotype_order[1][13] = 11;

  genotype_order[2][0] = 2;
  genotype_order[2][1] = 7;
  genotype_order[2][2] = 9;
  genotype_order[2][3] = 11;
  genotype_order[2][4] = 12;
  genotype_order[2][5] = 13;
  genotype_order[2][6] = 0;
  genotype_order[2][7] = 1;
  genotype_order[2][8] = 3;
  genotype_order[2][9] = 4;
  genotype_order[2][10] = 5;
  genotype_order[2][11] = 6;
  genotype_order[2][12] = 8;
  genotype_order[2][13] = 10;
  
  genotype_order[3][0] = 3;
  genotype_order[3][1] = 10;
  genotype_order[3][2] = 8;
  genotype_order[3][3] = 11;
  genotype_order[3][4] = 12;
  genotype_order[3][5] = 13;
  genotype_order[3][6] = 1;
  genotype_order[3][7] = 0;
  genotype_order[3][8] = 2;
  genotype_order[3][9] = 4;
  genotype_order[3][10] = 5;
  genotype_order[3][11] = 6;
  genotype_order[3][12] = 7;
  genotype_order[3][13] = 9;

  for(i=0;i<MAX_GENOTYPES;i++)
  {
      int a, b;
      for(j=0;j<4;j++)
      {
	  for(k=0;k<NO_ALLELES;k++)
	    allele_counts[j][i][k] = 0;
	  get_het_alleles(i,&a,&b,j);
	  allele_counts[j][i][a]++;
	  if(!HAPLOID)
	    allele_counts[j][i][b]++;		
      }
  }
  
  unsigned int *current_base;
  unsigned short **data;
  unsigned int *base_count;
  double *mean;
  int *median;
  int *max_coverage;
  unsigned int *tot_1x;
  unsigned int *tot_8x;
  unsigned int **counts;
  data = usmatrix(0,no_files,0,5);
  counts = umatrix(0,no_files,0,MAX_DIST);
  current_base = uvector(0,no_files);
  mean = dvector(0,no_files);
  median = ivector(0,no_files);
  base_count = uvector(0,no_files);
  tot_1x = uvector(0,no_files);
  tot_8x = uvector(0,no_files);
  max_coverage = ivector(0,no_files);
  
  for(i=0;i<no_files;i++)
  {
      mean[i] = 0.0;
      base_count[i] = 0;
      current_base[i] = 0;
      tot_1x[i] = 0;
      tot_8x[i] = 0;
      median[i] = 0;
      max_coverage[i] = 0;
      for(j=0;j<MAX_DIST;j++)
	counts[i][j] = 0;
  }

  int running_files = no_files;
  int this_count;

  // Init threads 
  printf("\n About to initialize threads and mutexes \n\n");

  pthread_mutex_init(&outfile_write_mutex, NULL);
  pthread_mutex_init(&snpfile_write_mutex, NULL);

  PTHREAD_DATA_NODE **thread_data;
  thread_data = (PTHREAD_DATA_NODE **)malloc(sizeof(PTHREAD_DATA_NODE *)*no_threads);
  if(!thread_data)
    dump_error("\n Can not allocate thread_data \n");

  int rc;
  threads = (pthread_t *) malloc(sizeof(pthread_t)*no_threads);
  if(!threads)
    dump_error("\n Could not allocate ram for threads \n");
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);

  for(i=0;i<no_threads;i++)
  {
    // printf("\n Working on thread %d \n\n",i);
      thread_data[i] = (PTHREAD_DATA_NODE *)malloc(sizeof(PTHREAD_DATA_NODE));
      if(!thread_data[i])
	dump_error("\n Can not allocate thread_data 2 \n");
      pthread_mutex_init(&(thread_data[i]->mutex),NULL);
      // pthread_mutex_lock(&(thread_data[i]->mutex));
      thread_data[i]->tid = i;
      thread_data[i]->status = DATA_EMPTY;
      thread_data[i]->chrom = AUTO;
      thread_data[i]->pos = 0;
      thread_data[i]->tot_pos = 0;
      thread_data[i]->samples = (SAMNODE **)malloc(sizeof(SAMNODE *)*INDIV);
      thread_data[i]->HAPLOID = HAPLOID;
      for(j=0;j<INDIV;j++)
      {
	  int kk = sorted_samples[j];
	  thread_data[i]->samples[kk] = sample_alloc();
	  thread_data[i]->samples[kk]->which = samples[j]->which;
	  thread_data[i]->samples[kk]->sex = samples[j]->sex;
	  strcpy(thread_data[i]->samples[kk]->indiv,samples[j]->indiv);
	  strcpy(thread_data[i]->samples[kk]->family,samples[j]->family);
      }
      for(j=0;j<INDIV;j++)
      {
	  int kk = sorted_samples[j];
	  if(samples[j]->dad)
		thread_data[i]->samples[kk]->dad = thread_data[i]->samples[samples[j]->dad->which];
	  else
		thread_data[i]->samples[kk]->dad = NULL;
	  if(samples[j]->mom)
		thread_data[i]->samples[kk]->mom = thread_data[i]->samples[samples[j]->mom->which];
	  else
		thread_data[i]->samples[kk]->mom = NULL;
      }
      /* for(j=0;j<INDIV;j++)
      {
	  printf("\n j = %d  which = %d is %s %s",j,thread_data[i]->samples[j]->which,
		thread_data[i]->samples[j]->family,thread_data[i]->samples[j]->indiv);
	  if(thread_data[i]->samples[j]->dad)
	    	printf("\n\tDad is %d",thread_data[i]->samples[j]->dad->which);
	  if(thread_data[i]->samples[j]->mom)
	    	printf("\n\tMom is %d",thread_data[i]->samples[j]->mom->which);
      }*/
      rc = pthread_create(&threads[i],&attr,call_single_base, (void *)thread_data[i]);
      if (rc)
      {
	  printf("ERROR; return code from pthread_create() is %d\n", rc);
	  exit(-1);
      }
  }

  printf("\n Finished initializing threads and mutex \n\n");
  for(i=0;i<no_files;i++)
  {
      while( (current_base[i] == 0) && !gzeof(pileupfile[i]))
      {
	  this_count = gzread(pileupfile[i],(void *)&current_base[i],sizeof(unsigned int));
	  if(this_count != 0)
	    gzread(pileupfile[i],(void *)data[i],sizeof(unsigned short)*6);
	  else
	      current_base[i] = 0;
	      
      }
      if(current_base[i] == 0)
	  running_files--;
  }

  int which_thread = 0;
  fprintf(snpfile,"Fragment\tPosition\tReference\tMinor_allele\tType");
  gzprintf(outfile,"Fragment\tPosition\tReference");
  gzprintf(pilefile,"Fragment\tPosition\tReference");
  for(i=0;i<INDIV;i++)
  {
      fprintf(snpfile,"\t%s\t",filenames[sorted_samples[i]]);
      gzprintf(outfile,"\t%s\t",filenames[sorted_samples[i]]);
      gzprintf(pilefile,"\t%s\t\t\t\t\t",filenames[sorted_samples[i]]);
  }
  if(!use_guide)
  {
      while(running_files > 0)
      {
	  while(thread_data[which_thread]->status != DATA_EMPTY)
	  {
	      which_thread++;
	      which_thread %= no_threads;
	  }
	  // pthread_mutex_lock(&(thread_data[which_thread]->mutex));

	  unsigned int lowest = find_lowest(current_base,no_files);
	  which = find_chrom(frag_pos,0,no_contigs,7,lowest);
	  char ref = get_genome(lowest,reffile);

	  thread_data[which_thread]->dom = ref;
	  thread_data[which_thread]->dom_int = gen_to_int(ref);
	  strcpy(thread_data[which_thread]->fragment,contig_names[which]);
	  thread_data[which_thread]->pos = 1 + lowest-frag_pos[which-1];
	  thread_data[which_thread]->chrom = chrom_type[which];
	  tot_bases++;
	  thread_data[which_thread]->tot_pos = tot_bases;
	  for(i=0;i<no_files;i++)
	    if(current_base[i] == lowest)
	    {
		for(j=0;j<6;j++)
		  thread_data[which_thread]->samples[sorted_samples[i]]->reads[j] = data[i][j];
		int tot_coverage = (data[i][0] + data[i][1] + data[i][2] + data[i][3] +data[i][4] + data[i][5]);
		mean[i] += (double)(tot_coverage);
		max_coverage[i] = maxim(max_coverage[i],tot_coverage);
		j = minim(tot_coverage,MAX_DIST-1);
		counts[i][j]++;
		base_count[i]++;

		if(!gzeof(pileupfile[i]))
		{
		    this_count = gzread(pileupfile[i],(void *)&current_base[i],sizeof(unsigned int));
		    if(this_count != 0)
		      gzread(pileupfile[i],(void *)data[i],sizeof(unsigned short)*6);
		    else
		    {
			current_base[i] = 0;
			running_files--;
		    }
		}
		else
		{
		    current_base[i] = 0;
		    running_files--;
		}
	    }
	    else
	      for(j=0;j<6;j++)
		thread_data[which_thread]->samples[sorted_samples[i]]->reads[j] = 0;
	  thread_data[which_thread]->status = DATA_LOADED;
	  // pthread_mutex_unlock(&(thread_data[which_thread]->mutex));
	  which_thread++;
	  which_thread %= no_threads;
      }
  }
  else
  {
      char line[4096];
      fgets(line,4095,guide_file);
      char *token;
      token = strtok(line,"\t \n");
      which = -1;
      for(i=0;i<no_contigs;i++)
	if(strcmp(token,contig_names[i]) == 0)
	{
	    which = i;
	    i = no_contigs;
	}
      if(which < 0)
      {
	  printf("\n For line chrom %s \n",token);
	  exit(1);
      }
      unsigned int start;
      unsigned int end;
      token = strtok(NULL,"\t \n");
      start = frag_pos[which-1] + atoi(token) - 1;
      token = strtok(NULL,"\t \n");
      end = frag_pos[which-1] + atoi(token) - 1;
      unsigned int lowest = start;
      // printf("\n About to start circling our threads \n\n");

      while(running_files > 0 )
      {
	  char ref = get_genome(lowest,reffile);

	  while(thread_data[which_thread]->status != DATA_EMPTY)
	  {
	      which_thread++;
	      which_thread %= no_threads;
	  }

	  thread_data[which_thread]->dom = ref;
	  thread_data[which_thread]->dom_int = gen_to_int(ref);
	  strcpy(thread_data[which_thread]->fragment,contig_names[which]);
	  thread_data[which_thread]->pos = 1 + lowest-frag_pos[which-1];
	  thread_data[which_thread]->chrom = chrom_type[which];
	  thread_data[which_thread]->HAPLOID = HAPLOID;
	  if(thread_data[which_thread]->chrom == CHRY || thread_data[which_thread]->chrom == CHRMT )
	    thread_data[which_thread]->HAPLOID = TRUE;
	  tot_bases++;
	  thread_data[which_thread]->tot_pos = tot_bases;

	  // printf("\n Looking for contig = %d name = %s start = %u end = %u lowest = %u\n\n",which,contig_names[which],start,end,lowest);
	  
	  for(i=0;i<no_files;i++)
	  {
	      while(current_base[i] < lowest && current_base[i] > 0)
	      {
		  if(!gzeof(pileupfile[i]))
		  {
		      this_count = gzread(pileupfile[i],(void *)&current_base[i],sizeof(unsigned int));
		      if(this_count != 0)
			gzread(pileupfile[i],(void *)data[i],sizeof(unsigned short)*6);
		      else
		      {
			  current_base[i] = 0;
			  running_files--;
		      }
		  }
		  else
		  {
		      current_base[i] = 0;
		      running_files--;
		  }
	      }
	      if(current_base[i] == lowest)
	      {
		  for(j=0;j<6;j++)
		    thread_data[which_thread]->samples[sorted_samples[i]]->reads[j] = data[i][j];
		  int tot_coverage = (data[i][0] + data[i][1] + data[i][2] + data[i][3] +data[i][4] + data[i][5]);
		  mean[i] += (double)(tot_coverage);
		  max_coverage[i] = maxim(max_coverage[i],tot_coverage);
		  j = minim(tot_coverage,MAX_DIST-1);
		  counts[i][j]++;
		  base_count[i]++;

		  if(!gzeof(pileupfile[i]))
		  {
		      this_count = gzread(pileupfile[i],(void *)&current_base[i],sizeof(unsigned int));
		      if(this_count != 0)
			gzread(pileupfile[i],(void *)data[i],sizeof(unsigned short)*6);
		      else
		      {
			current_base[i] = 0;
			running_files--;
		      }
		  }
		  else
		  {
		      current_base[i] = 0;
		      running_files--;
		  }
	      }
	      else
	      {
		  for(j=0;j<6;j++)
		    thread_data[which_thread]->samples[sorted_samples[i]]->reads[j] = 0;
		  base_count[i]++;
	      }
	  }
	  thread_data[which_thread]->status = DATA_LOADED;
	  // pthread_mutex_unlock(&(thread_data[which_thread]->mutex));
	  which_thread++;
	  which_thread %= no_threads;
	  lowest++;
	  if(lowest > end)
	  {
	      line[0] = '\0';
	      if(!feof(guide_file))
		fgets(line,4095,guide_file);

	      if(strlen(line) < 5)
		running_files = 0;
	      else
	      {
		  char *token;
		  token = strtok(line,"\t \n");
		  which = -1;
		  for(i=0;i<no_contigs;i++)
		    if(strcmp(token,contig_names[i]) == 0)
		    {
			which = i;
			i = no_contigs;
		    }
		  if(which < 0)
		  {
		      printf("\n For line chrom %s \n",token);
		      exit(1);
		  }
		  token = strtok(NULL,"\t \n");
		  start = frag_pos[which-1] + atoi(token) - 1;
		  token = strtok(NULL,"\t \n");
		  end = frag_pos[which-1] + atoi(token) - 1;
		  lowest = start;
	      }
	  }
	  // printf("\n About to switch status of thread %d to DATA_LOADED \n\n",which_thread);
      }
  }

  ALL_FINISHED = TRUE;

  for(i=0;i<no_files;i++)
    if(base_count[i] > 0)
      mean[i] /= (double)base_count[i];

  for(i=0;i<no_files;i++)
  {   
      for(j=8;j<MAX_DIST;j++)
	tot_8x[i] += counts[i][j];
      tot_1x[i] = tot_8x[i] + counts[i][1] + counts[i][2] + counts[i][3] + counts[i][4] + counts[i][5] + counts[i][6] + counts[i][7];

      counts[i][0] = tot_bases-tot_1x[i];
      long median_count = counts[i][0];
      median[i] = 0;
      long stop = tot_bases / 2;
      for(j=1;j<MAX_DIST;j++)
      {
	  if(median_count > stop)
	    j = MAX_DIST;
	  else
	    median_count += counts[i][++median[i]];
      }
  }
  
  fprintf(distfile,"Category");
  for(i=0;i<no_files;i++)
    fprintf(distfile,"\t%s",filenames[i]);

  fprintf(distfile,"\nTotal Number of bases in target");
  for(i=0;i<no_files;i++)
    fprintf(distfile,"\t%u",tot_bases);
  fprintf(distfile,"\nTotal Number of bases with at least 1x coverage");
  for(i=0;i<no_files;i++)
    fprintf(distfile,"\t%u",tot_1x[i]);
  fprintf(distfile,"\nTotal Number of bases with at least 8x coverage");
  for(i=0;i<no_files;i++)
    fprintf(distfile,"\t%u",tot_8x[i]);
  fprintf(distfile,"\nMean depth of coverage");
  for(i=0;i<no_files;i++)
    fprintf(distfile,"\t%g",mean[i]);
  fprintf(distfile,"\nMedian depth of coverage");
  for(i=0;i<no_files;i++)
    fprintf(distfile,"\t%d",median[i]);
  fprintf(distfile,"\nMaximum depth of coverage");
  for(i=0;i<no_files;i++)
    fprintf(distfile,"\t%d",max_coverage[i]);
  fprintf(distfile,"\n\nDepth");
  for(j=0;j<MAX_DIST-1;j++)
  {
      fprintf(distfile,"\n%d",j);
      for(i=0;i<no_files;i++)
	fprintf(distfile,"\t%u",counts[i][j]);
  }
  fprintf(distfile,"\n%d+",MAX_DIST-1);
  for(i=0;i<no_files;i++)
    fprintf(distfile,"\t%u",counts[i][MAX_DIST-1]);
  fprintf(distfile,"\n");
  fclose(distfile);

  j = 0;
  for(i=0;i<no_threads;i++)
	while(thread_data[i]->status != DATA_ALL_DONE)
	{
		j++;
		if(j%1000000 == 0)
			printf("\n Currently stuck waiting for thread %d whose status is %d \n",i,thread_data[i]->status);
	} 
  gzclose(outfile);
  gzclose(pilefile);

  return 0;
  
}

/*-------------------------------------------------------------------------------------------------------------------------------------- */
void *call_single_base(void *threadid)
{
  PTHREAD_DATA_NODE *td;
  int i,j,ii,jj;
  int tid;
  int **alpha_prior,**first_alpha_prior,**new_alpha;
  int normal_factor;
  int last_pass = 5;
  char minor[80],this_minor[80];
  double **d_alpha_mean,**d_alpha_var,**d_alpha_weight;
  double coef_prior[MAX_GENOTYPES];
  int ind,total_configs,issnp;
  char **outline;
  char **snpline,**pline;
  char sss[4096];
  CNODE **configs;
  int max_outlines = 1000;
  int max_snplines = 10;



  outline = cmatrix(0,max_outlines,0,99999);
  snpline = cmatrix(0,max_snplines,0,99999);
  pline = cmatrix(0,max_snplines,0,99999);

  int s_line_c = 0;
  int o_line_c = 0;

  int max_configs = 500;


  td = (PTHREAD_DATA_NODE *)threadid;
  tid = td->tid;

  // printf("\n Thread %d is launched with status %d\n\n",tid,td->status);

  alpha_prior = imatrix(0,MAX_GENOTYPES-1,0,NO_ALLELES-1);
  first_alpha_prior = imatrix(0,MAX_GENOTYPES-1,0,NO_ALLELES-1);
  new_alpha = imatrix(0,MAX_GENOTYPES-1,0,NO_ALLELES-1);
  d_alpha_mean = dmatrix(0,MAX_GENOTYPES-1,0,NO_ALLELES-1);
  d_alpha_var = dmatrix(0,MAX_GENOTYPES-1,0,NO_ALLELES-1);
  d_alpha_weight = dmatrix(0,MAX_GENOTYPES-1,0,NO_ALLELES-1);
  


  configs = (CNODE **)malloc((unsigned) (max_gen*(max_configs+1)+1)*sizeof(CNODE *));
  if(!configs) dump_error("Allocation failure in configs\n");
    	
  ln_HW = HW_exact;
  // long iii = 0;

  while(!ALL_FINISHED)
    if(td->status == DATA_LOADED)
    {
      // pthread_mutex_lock(&(td->mutex));
      td->status = DATA_RUNNING;
      // printf("\n Thread %d has moved into running data status \n\n",tid);
      int tot_depth = 0;
      int dom_int = td->dom_int;
      char dom = td->dom;
      int expos = td->pos;
      int chrom = td->chrom;
      SAMNODE **samples = td->samples;
      int HAPLOID = td->HAPLOID;
      // printf("\n%s\t%d\t%c\n\n",td->fragment,expos,dom);

      for(ind=0;ind<INDIV;ind++)
      {		
	  samples[ind]->tot = samples[ind]->reads[0];
	  for(i=1;i<NO_ALLELES;i++)
	    samples[ind]->tot += samples[ind]->reads[i];
	  tot_depth += samples[ind]->tot;
	  if(samples[ind]->tot > 0)
	    for(i=0;i<NO_ALLELES;i++)
	      samples[ind]->frac[i] = (double)samples[ind]->reads[i] / (double)samples[ind]->tot;
	  samples[ind]->final_p = 0.0;

	  // printf("\n Just sucked in data for individual %d  tot = %d %d %d\n\n",ind,samples[ind]->tot,samples[ind]->reads[0],samples[ind]->reads[1]);
	  for(i=0;i<max_gen;i++)
	    samples[ind]->post_prob[i] = 0.0;
	  samples[ind]->coef = factln(samples[ind]->tot);
	  for(i=0;i<NO_ALLELES;i++)
	    samples[ind]->coef -= factln(samples[ind]->reads[i]);
	  // printf("\n About to make call \n\n");
	  samples[ind]->initial_call = dom_int;
	  // printf("\n Called %c \n\n",int_to_gen(samples[ind]->initial_call));
      }
      // printf("\n About to fill prior \n\n");
      // normal_factor = 20;
      // normal_factor = 5000;
      if(dom != 'N')
      {
	      normal_factor = maxim(10,tot_depth / INDIV);
	      fill_alpha_prior(alpha_prior,max_gen,normal_factor,normal_factor/2,dom_int);
	      for(ii=0;ii<max_gen;ii++)
		for(jj=0;jj<NO_ALLELES;jj++)
		  first_alpha_prior[ii][jj] = alpha_prior[ii][jj];
					
	      // printf("\n About to fill coef matrix \n\n");
	      fill_alpha_coef(alpha_prior,coef_prior,max_gen);
	      // printf("\n Back from coef fill \n\n");
	}
	      int calls_changed = TRUE;
      int pass = 0;
      last_pass = 0;
      while(calls_changed && pass < last_pass)
      {
	  pass++;
	  configs[0] = config_alloc(INDIV);
	  total_configs = 1;
	  
	  /* printf("\n Alpha Matrix\n");
	     for(ii=0;ii<max_gen;ii++)
	     {
	     for(jj=0;jj<NO_ALLELES;jj++)
	     printf("\t%d",alpha_prior[ii][jj]);
	     printf("\t\t%g\n",coef_prior[ii]);
	     } */ 

	  for(ind=0;ind<INDIV;ind++)
	    if(samples[ind]->tot > min_depth_needed && dom != 'N')
	    {
		// printf("\n About to fill config probs \n\n");
		total_configs = fill_config_probs(configs,total_configs,samples[ind],max_gen,alpha_prior,coef_prior,INDIV,ind,dom_int,chrom,HAPLOID);
		// printf("\n About to clean config probs with total configs = %d\n\n",total_configs);
		total_configs = clean_config_probs(configs,samples,total_configs,max_configs,max_gen,alpha_prior,coef_prior,INDIV,ind,dom_int,HAPLOID);
		// printf("\n Out of clean config probs with total configs = %d\n\n",total_configs);
	    }
	    else
	    {
		samples[ind]->final_call = MAX_GENOTYPES;
		for(i=0;i<max_gen;i++)
		  samples[ind]->post_prob[i] = 0.0;
		for(i=0;i<total_configs;i++)
		  configs[i]->sample_calls[ind] = MAX_GENOTYPES;
	        samples[ind]->final_p = 0.0;
	    }
	  // printf("\n About to fill prior with total configs = %d\n\n",total_configs);

	  double max_post = configs[0]->post;
	  double tot_post = 0;
	  for(i=0;i<total_configs;i++)
	  {
	      configs[i]->post -= max_post;
	      configs[i]->post = exp(configs[i]->post);
	      tot_post += configs[i]->post;
	  }
	  for(i=0;i<total_configs;i++)
	    configs[i]->post /= tot_post;

	  for(ind=0;ind<INDIV;ind++)
	    for(i=0;i<max_gen;i++)
	      samples[ind]->post_prob[i] = 0;

	  // printf("\n About to do the sample probs \n\n");
	  for(ind=0;ind<INDIV;ind++)
	    if(samples[ind]->tot > min_depth_needed && dom != 'N')
	      for(i=0;i<total_configs;i++)
		  samples[ind]->post_prob[(int)configs[i]->sample_calls[ind]] += configs[i]->post;
	  
	  // printf("\n About to check if the calls changed \n\n");
	  calls_changed = FALSE;
	  for(ind=0;ind<INDIV;ind++)
	    if(samples[ind]->tot > min_depth_needed && dom != 'N')
	    {
		int besti = 0;
		for(i=1;i<max_gen;i++)
		  if(samples[ind]->post_prob[i] > samples[ind]->post_prob[besti])
		    besti = i;
		samples[ind]->final_p = samples[ind]->post_prob[besti];
		samples[ind]->final_call = besti;
		if(samples[ind]->final_call != samples[ind]->initial_call)
		  calls_changed = TRUE;
		// printf("\nMaking call %c at p = %g",int_to_gen(besti),samples[ind]->post_prob[besti]); 
	    }
	  
	  if(INDIV < 4 || pass == last_pass)
	    calls_changed = FALSE;

	  if(calls_changed && dom != 'N')
	  {
	      for(ii=0;ii<max_gen;ii++)
		for(jj=0;jj<NO_ALLELES;jj++)
		  d_alpha_weight[ii][jj] = d_alpha_mean[ii][jj] = d_alpha_var[ii][jj] = (double)0.0;
	      
	      // printf("\n Calls have changed.  Recalculating alpha matrix  with total_configs = %d\n\n",total_configs);
	      for(i=0;i<total_configs;i++)
		for(ind=0;ind<INDIV;ind++)
		  if(samples[ind]->tot > min_depth_needed)
		    for(j=0;j<NO_ALLELES;j++)
		    {
			// printf("\n Working on individual %d of %d and allele %d of %d which is called %d\n\n",ind,INDIV,j,NO_ALLELES,configs[i]->sample_calls[ind]);
			d_alpha_mean[(int)configs[i]->sample_calls[ind]][j] += samples[ind]->frac[j]*configs[i]->post;
			d_alpha_var[(int)configs[i]->sample_calls[ind]][j] += (samples[ind]->frac[j]*samples[ind]->frac[j])*configs[i]->post;
			d_alpha_weight[(int)configs[i]->sample_calls[ind]][j] += configs[i]->post;
			// d_alpha_mean[(int)configs[i]->sample_calls[ind]][j] += samples[ind]->frac[j];
			// d_alpha_var[(int)configs[i]->sample_calls[ind]][j] += (samples[ind]->frac[j]*samples[ind]->frac[j]);
			// d_alpha_weight[(int)configs[i]->sample_calls[ind]][j] += 1.0;
		    }	      
	      // printf("\n About to finish average calcs \n\n");
	      for(ii=0;ii<max_gen;ii++)
		for(jj=0;jj<NO_ALLELES;jj++)
		  if(d_alpha_weight[ii][jj] > 1e-9)
		  {
		      d_alpha_mean[ii][jj] /= d_alpha_weight[ii][jj];
		      d_alpha_var[ii][jj] /= d_alpha_weight[ii][jj];
		      d_alpha_var[ii][jj] -= d_alpha_mean[ii][jj]*d_alpha_mean[ii][jj];
		  }

	      double var_eps = 1e-6;	
	      for(ii=0;ii<max_gen;ii++)
	      {
		  int non_zero_var = 0;
		  int this_min = 0;
		  int little_up = 0;
		  for(jj=1;jj<NO_ALLELES;jj++)
		    if(d_alpha_mean[ii][jj] > d_alpha_mean[ii][little_up])
		      little_up = jj;
		  
		  for(jj=0;jj<NO_ALLELES;jj++)
		  {
		      if(d_alpha_weight[ii][jj] >= 1.5 && d_alpha_var[ii][jj] > var_eps*d_alpha_mean[ii][jj])
			non_zero_var++;
		      if(d_alpha_mean[ii][jj] < d_alpha_mean[ii][this_min])
			this_min = jj;
		      if(d_alpha_mean[ii][jj] > var_eps && d_alpha_mean[ii][jj] < d_alpha_mean[ii][little_up])
			little_up = jj;
		      // printf("\n For ii = %d jj = %d  mean = %g  var = %g  weight = %g",ii,jj,d_alpha_mean[ii][jj],d_alpha_var[ii][jj],d_alpha_weight[ii][jj]);
		  }
		  if(non_zero_var > 1)
		  {
		      double s0 = 1.0;
		      for(jj=0;jj<NO_ALLELES;jj++)
			if(jj != this_min && d_alpha_var[ii][jj] > var_eps*d_alpha_mean[ii][jj])
			  s0 *= d_alpha_mean[ii][jj]*(1.0 - d_alpha_mean[ii][jj]) / d_alpha_var[ii][jj];
		      s0 = pow(s0-1.0,(double)1.0/(double)(non_zero_var-1.0));
		      s0 = maxim(s0,1.0/d_alpha_mean[ii][little_up]);
		      if(s0>3.0)
			for(jj=0;jj<NO_ALLELES;jj++)
			  alpha_prior[ii][jj] = maxim(1,(int)ceil(d_alpha_mean[ii][jj]*s0));
		      else
			for(jj=0;jj<NO_ALLELES;jj++)
			  alpha_prior[ii][jj] = first_alpha_prior[ii][jj];
		  }
		  else
		    for(jj=0;jj<NO_ALLELES;jj++)
		      alpha_prior[ii][jj] = first_alpha_prior[ii][jj];
	      }
	      
	      /* printf("\n Alpha Matrix\n");
		 for(ii=0;ii<max_gen;ii++)
		 {
		 for(jj=0;jj<NO_ALLELES;jj++)
		 printf("\t%d",alpha_prior[ii][jj]);
		 printf("\t\t\n");
		 } */
	      // printf("\n About to calculate coefficients \n\n");
		 check_alpha_sanity(alpha_prior,max_gen,first_alpha_prior,dom_int,d_alpha_weight);
	      fill_alpha_coef(alpha_prior,coef_prior,max_gen);
	      
	      for(ind=0;ind<INDIV;ind++)
		samples[ind]->initial_call = samples[ind]->final_call;
	  }
	  
	  // printf("\n About to free stuff\n\n");
	  for(i=0;i<total_configs;i++)
	    config_free(configs[i],INDIV);
      }
      sprintf(outline[o_line_c],"\n%s\t%d\t%c",td->fragment,expos,dom);
      issnp = REF_ALLELE;
      sprintf(minor,"N");
      for(ind=0;ind<INDIV;ind++)
	if(samples[ind]->tot > min_depth_needed && dom != 'N')
	{
	    sprintf(sss,"\t%c\t%g",int_to_gen(samples[ind]->final_call),samples[ind]->final_p);
	    strcat(outline[o_line_c],sss);
	    if(samples[ind]->final_p >= THRESHOLD)
	    {
		if(samples[ind]->final_call != dom_int && issnp != MESS)
		{
		    
		    int this_site = get_snp_type(samples[ind],dom_int,this_minor);
		    if(issnp == REF_ALLELE)
		      issnp = this_site;
		    else
		      if(issnp != this_site)
		      {
			  if(this_site <= INSERTION && issnp <= INSERTION)
			    issnp = MESS;
			  else
			    if(issnp > INSERTION && this_site <= INSERTION)
			      issnp = this_site;
			    else
			      if(issnp > INSERTION && this_site > INSERTION)
				issnp = MESS;
		      }
		    if(issnp < MESS && issnp == this_site)
		      strcpy(minor,this_minor);
		}
	    }

	}
	else
	{
	    strcat(outline[o_line_c],"\tN\t1");
	    samples[ind]->final_call = MAX_GENOTYPES;
	}

      //if(issnp)
      {
	  
	  if(use_ped)
	  {
	      int d_count = 0;
	      for(i=0;i<INDIV;i++)
		if(samples[i]->final_p >= THRESHOLD)
		{
		    int dad_called = MAX_GENOTYPES;
		    int mom_called = MAX_GENOTYPES;
		    int kid_called = samples[i]->final_call;
		    if(samples[i]->dad)
		      if(samples[i]->dad->final_p >= THRESHOLD)
			dad_called = samples[i]->dad->final_call;
		    if(samples[i]->mom)
		      if(samples[i]->mom->final_p >= THRESHOLD)
			mom_called = samples[i]->mom->final_call;
		    d_count += add_denovo(kid_called,dad_called,mom_called,samples[i]->sex,chrom);
		}
	      if(d_count > 0)
		sprintf(sss,"DENOVO_%s",snp_type[issnp]);
	      else
		sprintf(sss,"%s",snp_type[issnp]);
	  }
	  else
	    sprintf(sss,"%s",snp_type[issnp]);
	  
	  sprintf(snpline[s_line_c],"\n%s\t%d\t%c\t%s\t%s",td->fragment,td->pos,td->dom,minor,sss);
	  sprintf(pline[s_line_c],"\n%s\t%d\t%c",td->fragment,td->pos,td->dom);       
	  for(ind=0;ind<INDIV;ind++)
	  {
	      sprintf(sss,"\t%c\t%g",int_to_gen(samples[ind]->final_call),samples[ind]->final_p);
	      strcat(snpline[s_line_c],sss);
	      for(i=0;i<6;i++)
	      {
		  sprintf(sss,"\t%d",samples[ind]->reads[i]);
		  strcat(pline[s_line_c],sss);
	      }
	  }
	 // printf("\n In thread = %d and found snp %d",tid,s_line_c);
	  s_line_c++;

      }
      td->status = DATA_EMPTY;

      o_line_c++;
      if(o_line_c == max_outlines)
      {
	  pthread_mutex_lock(&(outfile_write_mutex));
	  for(i=0;i<o_line_c;i++)
	    gzprintf(outfile,"%s",outline[i]);
	  pthread_mutex_unlock(&(outfile_write_mutex));
	  o_line_c = 0;
      }
      if(s_line_c == max_snplines)
      {
	  pthread_mutex_lock(&(snpfile_write_mutex));
	  // printf("\n Just got my mutex.  I am now writing the snp stuff \n");
	  for(i=0;i<s_line_c;i++)
	  {
	      gzprintf(pilefile,"%s",pline[i]);
	      // fprintf(snpfile,"%s",snpline[i]);
	      // printf("\n%s\n%s",pline[i],snpline[i]);
	  }
	  pthread_mutex_unlock(&(snpfile_write_mutex));
	  s_line_c = 0;
      }
  }
  // printf("\n Got here 0 in thread %d \n\n",tid);

  if(o_line_c > 0)
  {
      pthread_mutex_lock(&(outfile_write_mutex));
      for(i=0;i<o_line_c;i++)
	gzprintf(outfile,"%s",outline[i]);
      pthread_mutex_unlock(&(outfile_write_mutex));
  }
  // printf("\n Got here 1 in thread %d \n\n",tid);

  if(s_line_c > 0)
  {
      pthread_mutex_lock(&(snpfile_write_mutex));
      for(i=0;i<s_line_c;i++)
      {
	  gzprintf(pilefile,"%s",pline[i]);
	  // fprintf(snpfile,"%s",snpline[i]);
      }
      pthread_mutex_unlock(&(snpfile_write_mutex));
  }
  td->status = DATA_ALL_DONE;
  // printf("\n Exiting all done in thread %d \n\n",tid);
  pthread_exit(NULL);
     
}
/*-------------------------------------------------------------------------------------------------------------------------------------- */

#define buffer_chunk 50000000
static long current_start;
static long current_end;
static char *genome_buffer;

void init_genome_buffer(gzFile *mfile)
{
  current_start = 0;
  current_end = buffer_chunk;
  genome_buffer = cvector(0,buffer_chunk);
  gzread(mfile,(void *)genome_buffer,sizeof(char)*buffer_chunk);
}
/*-------------------------------------------------------------------------------------------------------------------------------------- */
void free_genome_buffer(void)
{
  free_cvector(genome_buffer,0,buffer_chunk);
}
				     
/*-------------------------------------------------------------------------------------------------------------------------------------- */
char get_genome(unsigned int which,gzFile *mfile)
{
  long current_pos;
  if(which < current_start || which >= current_end)
  {
      current_start = maxim((long)which,0);
      current_end = current_start + buffer_chunk;
      gzseek(mfile,current_start,SEEK_SET);
      gzread(mfile,(void *)genome_buffer,sizeof(char)*buffer_chunk);
  }
 
  current_pos = which - current_start;
  return genome_buffer[current_pos];
}			     

/*-------------------------------------------------------------------------------------------------------------------------------------- */

int find_chrom(unsigned int *pos,int first, int last,int try,unsigned this)
{
  // printf("\n first = %d last = %d try  = %d this = %u pos[try] = %u",first,last,try,this,pos[try]);
  if(first >= try)
  {
      if(this>pos[first])
	return first+1;
      else
	return first;
  }
  if(last <= try)
    return last;

  if(pos[try] < this)
    return find_chrom(pos,try,last,(last+try)/2,this);

  if(pos[try] > this)
    return find_chrom(pos,first,try,(try+first)/2,this);
  
  return try+1;
}


/*---------------------------------------------------------------------*/
unsigned int find_lowest(unsigned *list,int n)
{
  unsigned int low = 0;
  int i = 0;
  while(low < 1 && i < n)
    low = list[i++];

  for(;i<n;i++)
    if(list[i] > 0 && list[i] < low)
	low = list[i];

  return low;
}
/*---------------------------------------------------------------------*/


char *cvector(int nl,int nh)
{
	char *v;
 
	v=(char *)malloc((unsigned) (nh-nl+1)*sizeof(char));
	if (!v) dump_error("allocation failure in cvector()");
	return v-nl;
}

int *ivector(int nl,int nh)
{
	int *v;
 
	v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
	if (!v) dump_error("allocation failure in ivector()");
	return v-nl;
}
unsigned int *uvector(int nl,int nh)
{
	unsigned int *v;
 
	v=(unsigned int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
	if (!v) dump_error("allocation failure in uvector()");
	return v-nl;
}

double *dvector(int nl, int nh)
{
	double *v;
 
	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) dump_error("allocation failure in dvector()");
	return v-nl;
}
 
 
int **imatrix(int nrl,int nrh,int ncl,int nch)
{
	int i,**m;

	m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
	if (!m) dump_error("allocation failure 1 in imatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
		if (!m[i]) dump_error("allocation failure 2 in imatrix()");
		m[i] -= ncl;
	}
	return m;
}

void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

 
double **dmatrix(int nrl,int nrh,int ncl,int nch)
{
	int i;
	double **m;
 
	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) dump_error("allocation failure 1 in dmatrix()");
	m -= nrl;
 
	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) dump_error("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}

unsigned short **usmatrix(int nrl,int nrh,int ncl,int nch)
{
	int i;
	unsigned short **m;

	m=(unsigned short **)malloc((unsigned) (nrh-nrl+1)*sizeof(unsigned short*));
	if (!m) dump_error("allocation failure 1 in cmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(unsigned short *)malloc((unsigned) (nch-ncl+1)*sizeof(unsigned short));
		if (!m[i]) dump_error("allocation failure 2 in cmatrix()");
		m[i] -= ncl;
	}
	return m;
}
unsigned int **umatrix(int nrl,int nrh,int ncl,int nch)
{
	int i;
	unsigned int **m;

	m=(unsigned int **)malloc((unsigned) (nrh-nrl+1)*sizeof(unsigned int*));
	if (!m) dump_error("allocation failure 1 in cmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(unsigned int *)malloc((unsigned) (nch-ncl+1)*sizeof(unsigned int));
		if (!m[i]) dump_error("allocation failure 2 in cmatrix()");
		m[i] -= ncl;
	}
	return m;
}

unsigned char **ucmatrix(int nrl,int nrh,int ncl,int nch)
{
	int i;
	unsigned char **m;

	m=(unsigned char **)malloc((unsigned) (nrh-nrl+1)*sizeof(unsigned char*));
	if (!m) dump_error("allocation failure 1 in cmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(unsigned char *)malloc((unsigned) (nch-ncl+1)*sizeof(unsigned char));
		if (!m[i]) dump_error("allocation failure 2 in cmatrix()");
		m[i] -= ncl;
	}
	return m;
}
 
void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch)
{
	int i;
 
	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((double*) (m+nrl));
}

char **cmatrix(int nrl,int nrh,int ncl,int nch)
{
	int i;
	char **m;
 
	m=(char **) malloc((unsigned) (nrh-nrl+1)*sizeof(char*));
	if (!m) dump_error("allocation failure 1 in cmatrix()");
	m -= nrl;
 
	for(i=nrl;i<=nrh;i++) {
		m[i]=(char *) malloc((unsigned) (nch-ncl+1)*sizeof(char));
		if (!m[i]) dump_error("allocation failure 2 in cmatrix()");
		m[i] -= ncl;
	}
	return m;
}
 
void free_cmatrix(char **m,int nrl,int nrh,int ncl,int nch)
{
	int i;
 
	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}
void free_ucmatrix(unsigned char **m,int nrl,int nrh,int ncl,int nch)
{
	int i;
 
	for(i=nrh;i>=nrl;i--) free((unsigned char*) (m[i]+ncl));
	free((unsigned char*) (m+nrl));
}
 
void free_cvector(char *v, int nl, int nh)
{
	free((char*) (v+nl));
}
 
void free_ivector(int *v, int nl, int nh)
{
	free((int*) (v+nl));
}
void free_uvector(unsigned int *v, int nl, int nh)
{
	free((int*) (v+nl));
}
 
void free_dvector(double *v, int nl,int nh)
{
	free((double*) (v+nl));
}
 
/*---------------------------------------------------------------------*/

/*---------------------------------------------------------------------*/
int get_snp_type(SAMNODE *sm,int ref,char *minor)
{
  int a,b;

  get_het_alleles(sm->final_call,&a,&b,ref);
  if(a < 4 && b < 4)
  {
      if(a != ref)
	sprintf(minor,"%c",int_to_gen(a));
      else
	sprintf(minor,"%c",int_to_gen(b));
 
      return SNP; 
  }

  sprintf(minor,"-");

  if(b == 4)
      return DELETION;

  if(b == 5)
      return INSERTION;

  return MESS;
  
}
/*---------------------------------------------------------------------*/
void check_alpha_sanity(int **alpha_prior,int max_gen,int **first_alpha_prior,int ref,double **weight)
{
  int i,j;
  int max;
  double scale_factor[MAX_GENOTYPES];


  // Homozygote call must have most reads of its own base
  for(i=0;i<4;i++)
  {
      max = 0;
      for(j=1;j<NO_ALLELES;j++)
	if(alpha_prior[i][j] > alpha_prior[i][max])
	  max = j;
      if(max != i)
	for(j=0;j<NO_ALLELES;j++)
	  alpha_prior[i][j] = first_alpha_prior[i][j];
  }

  double frac[MAX_GENOTYPES][NO_ALLELES];
  double temp;
  for(i=0;i<max_gen;i++)
  {
      int tot = alpha_prior[i][0];
      for(j=1;j<NO_ALLELES;j++)
	tot += alpha_prior[i][j];
      for(j=0;j<NO_ALLELES;j++)
	frac[i][j] = (double)alpha_prior[i][j] / (double)tot;
  }

  // Del Homozygote
  i = 4;
  temp = frac[i][i] - frac[ref][i];

  if(temp < 0.5)
    for(j=0;j<NO_ALLELES;j++)
      alpha_prior[i][j] = first_alpha_prior[i][j];
  
  // Ins Homozygote
  i = 5;
  temp = frac[i][i] - frac[ref][i];

  if(temp < 0.7)
    for(j=0;j<NO_ALLELES;j++)
      alpha_prior[i][j] = first_alpha_prior[i][j];

  // printf("\n In check sanity with max = %d \n",max_gen);
  for(i = NO_ALLELES;i<max_gen;i++)
  {
      // HET
      int a,b;
      get_het_alleles(i,&a,&b,ref);
      if(b == ref)
      {
	  int ti;
	  ti = a;	  a = b;	  b = ti;
      }
      if(frac[i][b] - frac[ref][b] < 0.15)
	for(j=0;j<NO_ALLELES;j++)
	  alpha_prior[i][j] = first_alpha_prior[i][j];
      else
      {
	  if(ref == a)
	    frac[i][a] -= 0.05;
	  else
	    frac[i][a] -= maxim(frac[ref][a],0.05);

	  frac[i][b] -= maxim(0.05,frac[ref][b]);


	  int bad = FALSE;
	  
	  for(j=0;j<NO_ALLELES;j++)
	    if(j != a && j != b)
	      if(frac[i][j] > frac[i][a] || frac[i][j] > frac[i][b])
	      {
		  bad = TRUE;
		  j = NO_ALLELES;
	      }
	  
	  if(bad)
	    for(j=0;j<NO_ALLELES;j++)
	      alpha_prior[i][j] = first_alpha_prior[i][j];
      }
  }

  int max_weight = 0;
  for(i=0;i<max_gen;i++)
  {  
      scale_factor[i] = alpha_prior[i][0];
      for(j=1;j<NO_ALLELES;j++)
	scale_factor[i] += alpha_prior[i][j];
      if(weight[i][0] > weight[max_weight][0])
	max_weight = i;
  }

  for(i=0;i<max_gen;i++)
    if(i != max_weight)
      scale_factor[i] = scale_factor[max_weight] / scale_factor[i];

  for(i=0;i<max_gen;i++)
    if(i != max_weight && weight[i][0] < 4.0)
      for(j=0;j<NO_ALLELES;j++)
      alpha_prior[i][j] = maxim(1,(int)ceil(scale_factor[i]*(double)alpha_prior[i][j])); 

}
/*---------------------------------------------------------------------*/
void get_het_alleles(int i, int *a, int *b,int ref)
{
  if(i < NO_ALLELES)
  {
      *a = *b = i;
  }
  else
  if(i == 6)
  {
      *a = 0;
      *b = 1;
  }
  else
  if(i == 7)
  {
      *a = 0;
      *b = 2;
  }
  else
  if(i == 8)
  {
      *a = 0;
      *b = 3;
  }
  else
  if(i == 9)
  {
      *a = 1;
      *b = 2;
  }
  else
  if(i == 10)
  {
      *a = 1;
      *b = 3;
  }
  else
  if(i == 11)
  {
      *a = 2;
      *b = 3;
  }
  else
  if(i == 12)
  {
      *a = ref;
      *b = 4;
  }
  else
  if(i == 13)
  {
      *a = ref;
      *b = 5;
  }
  else
  {
      printf("\n This is impossible in get_het_alleles.  i = %d\n\n",i);
      exit(1);
  }

  return;
}
/*---------------------------------------------------------------------*/
int clean_config_probs(CNODE **cn, SAMNODE **sn, int n, int max,int max_gen,int **alpha,double *coef,int indiv,int depth,int ref,int HAPLOID)
{
  int i;
  
  qsort(cn,n,sizeof(CNODE *),sort_configs);
  /* for(i=0;i<n;i++)
  {
      printf("\n For configuration %d like = %g  prior = %g  post = %g",i,cn[i]->like,cn[i]->prior,cn[i]->post);
      int j;
      for(j=0;j<max_gen;j++)
	printf(" %d",cn[i]->genotype_count[j]); 
	} */
  max = minim(max,n);
  for(i=1;i<max;i++)
    if(cn[0]->post > cn[i]->post + config_threshold)
      max = i;

  for(i=max;i<n;i++)
    config_free(cn[i],indiv);

  int found_hom = FALSE;
  for(i=0;i<max;i++)
    if(cn[i]->no_alleles == 1)
    {
	found_hom = TRUE;
	i = max;
    }


  // printf("\n Finished first round of cleaning, with max = %d and found_hom = %d \n",max,found_hom);
  if(!found_hom)
  {
      int best_hom = 0;
      for(i=1;i<NO_ALLELES;i++)
	if(cn[0]->allele_count[i] > cn[0]->allele_count[best_hom])
	  best_hom = i;
      if(best_hom > 3)
	best_hom = ref;
      cn[max] = config_alloc(indiv);
      cn[max]->genotype_count[best_hom] = depth+1;
      for(i=0;i<=depth;i++)
	cn[max]->sample_calls[i] = best_hom;
      if(HAPLOID)
	cn[max]->allele_count[best_hom] = depth+1;
      else
	cn[max]->allele_count[best_hom] = 2*(depth+1);
      cn[max]->no_alleles = 1;

      cn[max]->prior = 0;
      fill_hom_like(cn[max],sn,depth,alpha,coef);
      cn[max]->post = cn[max]->prior + cn[max]->like;
      max++;
  }
  // printf("\n Leaving the clean \n\n");
  return max;
}
/*---------------------------------------------------------------------*/
void fill_hom_like(CNODE *cn, SAMNODE **sn, int n,int **alpha,double *coef)
{
  int i, j, ii;
  double prob;

  cn->like = 0;
  for(i=0;i<=n;i++)
    if(sn[i]->tot >= min_depth_needed)
    {
	j = cn->sample_calls[i];
	prob = sn[i]->coef + coef[j];
	int tot = 0;
	for(ii=0;ii<NO_ALLELES;ii++)
	  {
	    prob += factln(alpha[j][ii]+sn[i]->reads[ii]-1);
	    tot += alpha[j][ii] + sn[i]->reads[ii];
	  }
	prob -= factln(tot-1);
	cn->like += prob;
    } 
}
/*---------------------------------------------------------------------*/
int sort_configs(const void *a,const void *b)
{
  CNODE *fa, *fb;
  
  fa = *((CNODE **)a);
  fb = *((CNODE **)b);

  if(fa->post > fb->post)
    return -1;
  else
    if(fa->post < fb->post)
      return 1;
    else
      return 0;
}
/*---------------------------------------------------------------------*/
int add_denovo(int kid, int dad, int mom,int sex,int chrom)
{
  if(dad < MAX_GENOTYPES)
  {
      if(mom < MAX_GENOTYPES) // BOTH
      {
	  if(chrom == AUTO)
	    return trio_denovo[dad][mom][kid];

	  if(chrom == CHRX)
	  {
	      if(sex == 1)
		return dyad_denovo[mom][kid];
	      else
		return trio_denovo[dad][mom][kid];
	  }

	  if(chrom == CHRY)
	  {
	      if(sex == 1)
		return dyad_denovo[dad][kid];
	      else
		return 0;
	  }

	  if(chrom == CHRMT)
	    return dyad_denovo[mom][kid];

	  return 0;
      }
      else // DAD ONLY
      {
	  if(chrom == AUTO)
	    return dyad_denovo[dad][kid];
	  if( (chrom == CHRX) && (sex == 2))
	    return dyad_denovo[dad][kid];
	  if( (chrom == CHRY) && (sex == 1) )
	    return dyad_denovo[dad][kid];
	  
	  return 0;
      }
  }
  
  if(mom < MAX_GENOTYPES)  // MOM ONLY
    if(chrom != CHRY )
      return dyad_denovo[mom][kid];
  
  return 0;
}
/*---------------------------------------------------------------------*/

int fill_config_probs(CNODE **cn, int n, SAMNODE *sn, int max,int **alpha,double *coef,int indiv,int this_depth,int ref,int chrom,int HAPLOID)
{
  int i,j,jj,k,newcount,ii;
  CNODE *temp,*old,**new;
  double prob;
  double best = -1e100;
  
  newcount = 0;
  new = (CNODE **)malloc((unsigned) ((max+1)*(n+1)*sizeof(CNODE *)));
  if(!new) dump_error("Allocation failure in fill_config_probs\n");
  // printf("\n Entering fill_config_probs with n = %d \n\n",n);
  for(i=0;i<n;i++)
  {
      // printf("\n At top of the loop with i=%d and cn[i] = %d\n\n",i,cn[i]);
      old = cn[i];
      for(jj=0;jj<max;jj++)
      {
	  j = genotype_order[ref][jj];
	  // printf("\n ABout to allocate temp with ref=%d jj =%d and j = %d \n\n",ref,jj,j);
	  temp = config_alloc(indiv);

	  // printf("\n Copying over \n\n");
	  for(k=0;k<max;k++)
	      temp->genotype_count[k] = old->genotype_count[k];

	  for(k=0;k<this_depth;k++)
	    temp->sample_calls[k] = old->sample_calls[k];

	  temp->genotype_count[j]++;
	  temp->sample_calls[this_depth] = j;
	  prob = sn->coef + coef[j];
	  // printf("\n Calcing Prob \n\n");
	  int tot = 0;
	  for(ii=0;ii<NO_ALLELES;ii++)
	  {
	      prob += factln(alpha[j][ii]+sn->reads[ii]-1);
	      tot += alpha[j][ii] + sn->reads[ii];
	  }
	  prob -= factln(tot-1);
	  // printf("\nSample = %s Config = %d Genotype = %d  Prob = %g old = %g best = %g ",sn->indiv,i,jj,prob,old->like,best);
	  temp->like = old->like + prob;
	  if(temp->like + config_threshold < best)
	    config_free(temp,indiv);
	  else
	  {
	      temp->prior = 0.0;
	      // printf("\n Made it in here to store the config \n");
	      temp->hets = old->hets;
	      temp->homs = old->homs;
	      for(k=0;k<NO_ALLELES;k++)
		temp->allele_count[k] = old->allele_count[k];
	      temp->no_denovo = old->no_denovo;
	      // printf("\n no_denovo = %d for %s and call %d\n\n",temp->no_denovo,sn->indiv,j);
	      if(sn->dad)
		if(sn->mom) // Both
		  temp->no_denovo += add_denovo(j,(int)temp->sample_calls[sn->dad->which],(int)temp->sample_calls[sn->mom->which],sn->sex,chrom);
		else // Dad Only
		  temp->no_denovo += add_denovo(j,(int)temp->sample_calls[sn->dad->which],MAX_GENOTYPES,sn->sex,chrom);
	      else
		if(sn->mom) // Mom Only
		  temp->no_denovo += add_denovo(j,MAX_GENOTYPES,(int)temp->sample_calls[sn->mom->which],sn->sex,chrom);
		else // NEITHER
		{
		    for(k=0;k<NO_ALLELES;k++)
		      temp->allele_count[k] += allele_counts[ref][j][k];
		    if(j >= NO_ALLELES)
		      temp->hets++;
		    else
		      temp->homs++;
		}
	      // printf("\n\t\tAfterwards the number of denovo is %d\n\n",temp->no_denovo);
	      temp->no_alleles = 0;
	      for(k=0;k<NO_ALLELES;k++)
		if(temp->allele_count[k] > 0)
		  temp->no_alleles++;
	      
	      temp->prior = 0;
	      if(temp->no_alleles > 1)
		temp->prior = (temp->no_alleles-1) * ln_theta;
	      
	      // printf("\nWith Theta Prior is %g",temp->prior); 
	      if(temp->no_denovo > 0)
		temp->prior += temp->no_denovo * ln_denovo;
	      // printf("\nWith denovo Prior is %g",temp->prior); 
	      
	      if(!HAPLOID && temp->no_alleles > 1)
	      {
		  int major = 0;
		  int minor = 0;
		  for(k=0;k<NO_ALLELES;k++)
		    if(temp->allele_count[k] > major)
		    {
			minor += major;
			major = temp->allele_count[k];
		    }
		    else
		      minor += temp->allele_count[k];
		  
		  minor = minim(minor,major);
		  int hets = minim(minor,temp->hets);
		  int tot_n = (minor + major) / 2;
		  if( (minor - hets) % 2 == 1)
		  {
		      minor++;
		      major++;
		      // Bad mojo
		  }
		  // printf("\n Using minor = %d major = %d hets = %d tot_n = %d",minor,major,hets,tot_n);
		  temp->prior += ln_HW[tot_n][minor][hets];
	      }
	      
	      /* printf("\nWith HW Prior is %g",temp->prior); 
	      
	      printf("\n j = %d with n = %d newcount = %d Configuration:",j,n,newcount);
	      for(ii=0;ii<max;ii++)
		printf(" %d",temp->genotype_count[ii]);
	      printf("\n\tTotal Like = %g  This Prob = %g",temp->like,exp(prob));
	      for(ii=0;ii<NO_ALLELES;ii++)
		printf(" %d",sn->reads[ii]);  
		printf("\n About to store \n\n"); */
	      temp->post = temp->prior + temp->like;
	      if(temp->post + config_threshold > best)
	      {
		  best = maxim(best,temp->post);
		  new[newcount] = temp;
		  newcount++;
	      }
	      else
		config_free(temp,indiv);
	      
	  }
	  // printf("\n i = %d  about to free old = %d\n\n",i,old);
      }
      config_free(old,indiv);
  }
  for(i=0;i<newcount;i++)
	cn[i] = new[i];
  free(new);

  // printf("\n About to leave with newcount = %d\n\n",newcount);

  return newcount;
    
}
/*---------------------------------------------------------------------*/
void fill_hardy_weinberg(double **exact_HW,int asize,int n)
{
  double **marg,sum,p;
  int i,j,naa,nab,nbb,Na,Nb,start,expect;
  
  // printf("\n Entering fill_hardy_weinberg with asize = %d\n\n",asize);
  marg = dmatrix(0,asize,0,n);
  for(i=0;i<=asize;i++)
    for(j=0;j<=n;j++)
      exact_HW[i][j] = marg[i][j] = 0.0;


  for(i=1;i<=asize;i++)
  {
      Na = 2*n - i;
      Nb = i;
      p = (double)i/(double)(Na+Nb);
      expect = ceil(i*(1.0-p));
      
      if(i%2 == 0)
      {
	  if(expect%2 == 1)
	    start = expect -1;
	  else
	    start = expect;
      }
      else
      {
	  if(expect%2 == 1)
	    start = expect;
	  else
	    start = expect-1;
      }
      // printf("\ni=%d expect = %d  start = %d p = %g ",i,expect,start,p);
      sum = marg[i][start] = 1.0;

      nbb = ((Nb-start)/2);
      naa = ((Na-start)/2);

      // printf("\n naa = %d nbb = %d nab = %d  Nb = %d  Na = %d",naa,nbb,start+2,Nb,Na);
      for(nab=start+2;naa > 0 && nbb > 0;nab+=2,naa--,nbb--)
      {
	  marg[i][nab] = marg[i][nab-2] * 4.0 * ((double)naa*(double)nbb) / ((double)(nab) * (double)(nab-1.0) );
	  // printf("\n nab = %d last = %g  current = %g",nab,marg[i][nab-2],marg[i][nab]);
	  sum += marg[i][nab];
      }

      nbb = ((Nb-start)/2);
      naa = ((Na-start)/2);

      for(nab=start-2;nab>=0;nab-=2,naa++,nbb++)
      {
	  marg[i][nab] = marg[i][nab+2] * ((double)(nab+2.0) * (double)(nab+1.0) ) / ( (double)4.0*((double)(naa+1.0)*(nbb+1.0)));
	  sum += marg[i][nab];
      }
      
      for(j=0;j<=n;j++)
	marg[i][j] /= sum;
  }
  for(i=0;i<=asize;i++)
      for(j=0;j<=n;j++)
      {
	  exact_HW[i][j] = log(marg[i][j]);
	  // printf("\n i = %d  j = %d  p = %g",i,j,exact_HW[i][j]);
      }
  // exit(1);   
  free_dmatrix(marg,0,asize,0,n);
}

/*---------------------------------------------------------------------*/
int gen_to_int(char c)
{
  if(c == 'A') return 0;
  if(c == 'C') return 1;
  if(c == 'G') return 2;
  if(c == 'T') return 3;
  if(c == 'D') return 4;
  if(c == 'I') return 5;
  if(c == 'M') return 6;
  if(c == 'R') return 7;
  if(c == 'W') return 8;
  if(c == 'S') return 9;
  if(c == 'Y') return 10;
  if(c == 'K') return 11;
  if(c == 'E') return 12;
  if(c == 'H') return 13;
  if(c == 'N') return 14;

  printf("\n This is impossible\n Illegal character in gen_to_int %c\n\n",c);
  exit(1);
  return -1;

}
/*---------------------------------------------------------------------*/
char int_to_gen(int c)
{
  if(c == 0) return 'A';
  if(c == 1) return 'C';
  if(c == 2) return 'G';
  if(c == 3) return 'T';
  if(c == 4) return 'D';
  if(c == 5) return 'I';
  if(c == 6) return 'M';
  if(c == 7) return 'R';
  if(c == 8) return 'W';
  if(c == 9) return 'S';
  if(c == 10) return 'Y';
  if(c == 11) return 'K';
  if(c == 12) return 'E';
  if(c == 13) return 'H';

  return 'N';
}

/*---------------------------------------------------------------------*/

SAMNODE *sample_alloc()
{
	SAMNODE *tn;
	int i;

	tn = (SAMNODE *)malloc((unsigned) sizeof(struct sample_node));
	if (!tn) dump_error("allocation failure in sample_alloc()");

	for(i=0;i<=MAX_GENOTYPES;i++)
	  tn->post_prob[i] = 0.0;

	for(i=0;i<NO_ALLELES;i++)
	  tn->reads[i] = 0;
	tn->final_call = MAX_GENOTYPES;
	tn->initial_call = MAX_GENOTYPES;
	tn->final_p = 0.0;
	tn->coef = 0.0;

	tn->family[0] = '\0';
	tn->indiv[0] = '\0';
	tn->mom = NULL;
	tn->dad = NULL;
	tn->which = -1;
	tn->sex = 0;

	return tn;
}

/*---------------------------------------------------------------------*/

CNODE *config_alloc(int N)
{
	CNODE *tn;
	int i;

	tn = (CNODE *)malloc((unsigned) sizeof(struct config_node));
	if (!tn) dump_error("allocation failure in config_alloc()");

	for(i=0;i<MAX_GENOTYPES;i++)
	    tn->genotype_count[i] = 0;

	tn->sample_calls = cvector(0,N-1);
	for(i=0;i<N;i++)
	    tn->sample_calls[i] = (char)MAX_GENOTYPES;

	tn->like = 0;
	tn->prior = 0;
	tn->post = 0;

	tn->no_alleles = 0;
	for(i=0;i<NO_ALLELES;i++)
	  tn->allele_count[i] = 0;
	
	tn->no_denovo = 0;
	tn->hets = 0;
	tn->homs = 0;

	return tn;
}
/*---------------------------------------------------------------------*/
void config_free(CNODE *tn,int N)
{
  // printf("\n In free with tn = %d\n\n",tn);
  // printf("\n\t tn->sample_calls =  %d  tn->alpha = %d\n\n",tn->sample_calls,tn->alpha);
  if(!tn)
    return;

  free_cvector(tn->sample_calls,0,N);
  free(tn);
}
/*---------------------------------------------------------------------*/
void fill_alpha_prior(int **alpha,int max_gen,int hom,int het,int ref)
{
  int i,j,k,hom_err,err;


  hom_err = maxim(1,hom/300);
  err = maxim(1,(2*het)/300);

  for(i=0;i<max_gen;i++)
  {
        if(i<NO_ALLELES-2)
	  for(j=0;j<NO_ALLELES;j++)
	    if(i == j)
	      alpha[i][j] = hom;
	    else
	      alpha[i][j] = hom_err;
	else
	{
	    j = i;
	    if(j == NO_ALLELES - 2)
	    {
		// Del homozygote
		for(k=0;k<4;k++)
		  if(k == ref)
		    alpha[j][k] = hom/5;
		  else
		    alpha[j][k] = err;
		alpha[j][4] = (4*hom)/5;
		alpha[j][5] = err;
	    }
	    else
	    if(j == NO_ALLELES - 1)
	    {
		// Ins homozygote
		for(k=0;k<4;k++)
		  if(k == ref)
		    alpha[j][k] = hom;
		  else
		    alpha[j][k] = err;
		alpha[j][4] = err;
		alpha[j][5] = (4*hom)/5;
	    } 
	    else
	      if(j < NO_ALLELES+6)
	      {
		  int a,b;
		  get_het_alleles(j,&a,&b,ref);
		  if(a == ref)
		  {
		      alpha[j][a] = (51*het) / 50;
		      alpha[j][b] = (49*het) / 50;
		      alpha[j][4] = maxim(1,het/20);
		      alpha[j][5] = err;
		      for(k=0;k<4;k++)
			if(k != a && k != b)
			  alpha[j][k] = err;
		  }
		  else
		    if(b == ref)
		    {
			alpha[j][b] = (51*het) / 50;
			alpha[j][a] = (49*het) / 50;
			alpha[j][4] = maxim(1,het/20);
			alpha[j][5] = err;
			for(k=0;k<4;k++)
			  if(k != a && k != b)
			    alpha[j][k] = err;
		    }
		    else
		    {
			alpha[j][a] = het;
			alpha[j][b] = het;
			for(k=0;k<NO_ALLELES;k++)
			  if(k != a && k != b)
			    alpha[j][k] = err;
		    }
	      }
	      else
		if(j == NO_ALLELES + 6)
		  {
		    // Del - ref
		    alpha[j][4] = (4*het)/5;
		    alpha[j][ref] = (6*het)/5;
		    for(k=0;k<4;k++)
		      if(k != ref)
			alpha[j][k] = err;
		    alpha[j++][5] = err;
		  }
		else
		  {
		    // Ins - ref
		    alpha[j][5] = (2*het)/5;
		    alpha[j][ref] = (8*het)/5;
		    for(k=0;k<5;k++)
		      if(k != ref)
			alpha[j][k] = err;
		  }
	}
  }
}

/*---------------------------------------------------------------------*/
void fill_alpha_coef(int **alpha,double *coef,int max_gen)
{
  int i,j,tot;

  for(j=0;j<max_gen;j++)
  {
      coef[j] = 0;
      tot  = 0;
      for(i=0;i<NO_ALLELES;i++)
      {
	  tot += alpha[j][i];
	  coef[j] -= factln(alpha[j][i]-1);
      }
      coef[j] += factln(tot-1);
  }

}

/*---------------------------------------------------------------------*/

double gammln(double xx)
{
	double x,tmp,ser;
	static double cof[6]={76.18009173,-86.50532033,24.01409822,
		-1.231739516,0.120858003e-2,-0.536382e-5};
	int j;

	x=xx-1.0;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.0;
	for (j=0;j<=5;j++) {
		x += 1.0;
		ser += cof[j]/x;
	}
	return -tmp+log(2.50662827465*ser);
}

/*---------------------------------------------------------------------*/
double exactfactln(int n)
{
  int i;
  double x=1.0;

  for(i=2;i<=n;i++)
    x*=(double)i;

  return log(x);
}
/*---------------------------------------------------------------------*/

double factln(int n)
{
	static double a[10001];

	if (n < 0) dump_error("Negative factorial in routine FACTLN");
	if (n <= 1) return 0.0;
	if (n <= 40) return a[n] ? a[n] : (a[n]=exactfactln(n));
	if (n <= 10000) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
	else return gammln(n+1.0);
}

/* -------------------------------------------------------------- */
void dump_error(char *ss)
{
  printf("\n Unrecoverable error\n %s \n  Exiting now \n",ss);
  exit(1);
}
/* -------------------------------------------------------------- */

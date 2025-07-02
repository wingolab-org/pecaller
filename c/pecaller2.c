/*

The code itself is Copyright (C) 2018, by David J. Cutler.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
This library is distributed in the hope that it will be useful,freed
but WITHOUT ANY WARRANTY; without even the implied warranty offs
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public
License along with library; if not, write to the Free Software
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

#define SQR(a) ((a) * (a))
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
	double post_prob[MAX_GENOTYPES + 1];
	double like[MAX_GENOTYPES + 1];
	double final_p;
	double initial_p;
	double coef;
	int reads[NO_ALLELES];
	double frac[NO_ALLELES];
	int tot;
	char final_call;
	char initial_call;
	char family[1024];
	char indiv[1024];
	short sex;
	struct sample_node *mom;
	struct sample_node *dad;
	int no_kids;
	struct sample_node **kids;
	int which;
} SAMNODE;

typedef struct config_node
{
	int genotype_count[MAX_GENOTYPES + 1];
	double prior;
	double like;
	double post;
	char *sample_calls;
	short no_alleles;
	int allele_count[NO_ALLELES];
	int no_denovo;
	int hets;
	int homs;
} CNODE;

short dyad_denovo[NO_ALLELES][MAX_GENOTYPES + 1][MAX_GENOTYPES + 1];
short trio_denovo[NO_ALLELES][MAX_GENOTYPES + 1][MAX_GENOTYPES + 1][MAX_GENOTYPES + 1];
short allele_counts[NO_ALLELES][MAX_GENOTYPES][NO_ALLELES];

double LOW_BASE;

#define LOG10 (double)2.302585
#define MAX_MODELS 5000000
#define MAX_PRIOR_CALLS 65000

#define minim(atesta, btestb) ((atesta < btestb) ? atesta : btestb)
#define maxim(atesta, btestb) ((atesta > btestb) ? atesta : btestb)

#define AUTO 0
#define CHRX 1
#define CHRY 2
#define CHRMT 3

#define REF_ALLELE 0
#define SNP 1
#define DELETION 2
#define INSERTION 3
#define LOW 4
// #define HIGH 5
#define MULTI 5
#define MESS 6

#define DATA_EMPTY 0
#define DATA_LOADED 1
#define DATA_RUNNING 2
#define DATA_ALL_DONE 3

#define pileup_header_size 240

typedef struct pthread_data_node
{
	pthread_mutex_t mutex;
	int tid;
	int status;
	int chrom;
	int pos;
	char dom;
	int dom_int;
	char fragment[1024];
	unsigned short **reads;
	int HAPLOID;
	void *model_bucket;
	long times_used;
} PTHREAD_DATA_NODE;

int no_threads;
pthread_mutex_t outfile_write_mutex;
pthread_mutex_t snpfile_write_mutex;
pthread_t *threads;
double ln_denovo;
double starting_threshold = 9.2;
int use_ped;
char allele_char[NO_ALLELES + 1];
int **genotype_order;
int INDIV;

int gen_to_int(char c);
char int_to_gen(int c);
int clean_config_probs(CNODE **cn, SAMNODE **sn, int n, int max, int max_gen, int indiv, int depth, int ref, int HAPLOID, double thres, int *starting);
int sort_configs(const void *a, const void *b);
int fill_config_probs(CNODE **cn, int n, SAMNODE **samples, int max, int indiv, int this_depth, int ref, int chrom, int HAPLOID, double thres, int *starting);
SAMNODE *sample_alloc(int kids);
CNODE *config_alloc(int N, int dom, SAMNODE **sn, int is_haploid, int *starting, int chrom, int first_config);
void config_free(CNODE *tn, int N);
void fill_hardy_weinberg(double **exact_HW, int asize, int n);
void fill_alpha_prior(double **alpha, int max_gen, int dom_int);
void fill_alpha_coef(double **alpha, double *coef, int max_gen);
void fill_config_like(CNODE *cn, SAMNODE **sn, int n);
void fill_sample_like(SAMNODE **samples, double **alpha, int max_gen, int indiv, int dom_int, int pass);
int add_denovo(int kid, int dad, int mom, int sex, int chrom, int ref);
double gammln(double xx);
double exactfactln(int n);
double factln(int n);
int check_alpha_sanity(double **alpha_prior, double **weight, int max_gen, int ref);
void get_het_alleles(int i, int *a, int *b, int ref);
void read_var(char *line, char *result);
int *ivector(int nl, int nh);
unsigned int *uvector(int nl, int nh);
char *cvector(int nl, int nh);
double *dvector(int nl, int nh);
double **dmatrix(int nrl, int nrh, int ncl, int nch);
unsigned short **usmatrix(int nrl, int nrh, int ncl, int nch);
int **imatrix(int nrl, int nrh, int ncl, int nch);
char **cmatrix(int nrl, int nrh, int ncl, int nch);
void free_cvector(char *v, int nl, int nh);
void free_ivector(int *v, int nl, int nh);
void free_uvector(unsigned int *v, int nl, int nh);
void free_dvector(double *v, int nl, int nh);
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
unsigned int **umatrix(int nrl, int nrh, int ncl, int nch);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
void free_ucmatrix(unsigned char **m, int nrl, int nrh, int ncl, int nch);
void free_cmatrix(char **m, int nrl, int nrh, int ncl, int nch);
unsigned char **ucmatrix(int nrl, int nrh, int ncl, int nch);
void dump_error(char *error_text);
double chidist(double x, double lambda, int df);
double gammq(double a, double x);
double gammp(double a, double x);
void gser(double *gamser, double a, double x, double *gln);
void gcf(double *gammcf, double a, double x, double *gln);
double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin);
double calc_chi(double gamma);
int sort_compare(const void *a, const void *b);
int find_chrom(unsigned int *pos, int first, int last, int try, unsigned this);
unsigned int find_lowest(unsigned *list, int n);
void init_genome_buffer(gzFile mfile);
void free_genome_buffer(void);
char get_genome(unsigned int which, gzFile mfile);
void *call_single_base(void *threadid);
void fill_prior(CNODE *temp, int HAPLOID);
double get_HW_exact(int i, int j, int k);
void fill_first_alpha_prior(void);
double calc_angle(double *vec1, double mag1, double *vec2, double mag2, int n);
int count_denovos(SAMNODE **samples, int INDIV, int chrom, int dom_int);

int min_depth_needed;
int ALL_FINISHED = FALSE;
int dump_me = FALSE;
SAMNODE **all_samples;

double EPS = 1e-8;
#define MAX_DIST 501
gzFile outfile, pilefile;
double ***HW_exact;
double THRESHOLD, theta, ln_theta;
char snp_type[MESS + 1][80];
FILE *snpfile;
int use_saved_models;
int rewrite_models;
double ***default_alpha_prior;
double ***default_alpha_frac;
double **default_alpha_mag;

int main(int argc, char *argv[])
{
	FILE *sfile, *guide_file;
	gzFile reffile;
	char ss[4096], sdxname[4096], sss[4096];
	char guide_name[4096];
	char **filenames;
	DIR *pDIR;
	struct dirent *pDirEnt;
	int no_files, i, j, k;
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
	long current_byte_pos = pileup_header_size;
	int model_bucket_size, single_model_size;
	void *model_bucket;
	gzFile modelfile_read, modelfile_write;

	modelfile_read = modelfile_write = NULL;
	model_bucket = NULL;

	use_saved_models = FALSE;
	rewrite_models = FALSE;

	FILE *distfile;
	gzFile *pileupfile;

	if (argc < 11 || argc > 14)
	{
		printf("\n Usage %s pileup_extension sdx_file no_files outfile Prob_to_call theta haploid[y,n] no_threads rewrite_models[y,n] use_pedfile[y,n]  [pedfilename] [denovo_mutation_rate] [guide_file_bed_format]\n", argv[0]);
		exit(1);
	}

	no_threads = (int)atoi(argv[8]);
	if (no_threads < 2 || no_threads > 200)
	{
		printf("\n Number of threads is limited to 2 to 200.   You entered %d \n\n", no_threads);
		exit(1);
	}
	no_threads--;

	sprintf(ss, "%s.base.gz", argv[4]);
	if ((outfile = gzopen(ss, "w")) == (gzFile)NULL)
	{
		printf("\n Can not open file %s\n", ss);
		exit(1);
	}
	gzbuffer(outfile, 131072);

	sprintf(ss, "%s.snp", argv[4]);
	if ((snpfile = fopen(ss, "w")) == (FILE *)NULL)
	{
		printf("\n Can not open file %s for writing\n", ss);
		exit(1);
	}

	sprintf(ss, "%s.dist", argv[4]);
	if ((distfile = fopen(ss, "w")) == (FILE *)NULL)
	{
		printf("\n Can not open file %s for writing\n", ss);
		exit(1);
	}

	sprintf(ss, "%s.piles.gz", argv[4]);
	if ((pilefile = gzopen(ss, "w")) == (gzFile)NULL)
	{
		printf("\n Can not open file %s for writing\n", ss);
		exit(1);
	}
	gzbuffer(pilefile, 131072);

	guide_file = stdout;

	THRESHOLD = (double)atof(argv[5]);
	theta = (double)atof(argv[6]);
	if ((theta < 1e-10) || (theta > 0.5))
	{
		printf("\n Encountered impossible value for theta = %g \n", theta);
		exit(1);
	}
	ln_theta = log(theta);

	HAPLOID = FALSE;
	int max_gen;
	max_gen = min_depth_needed = 0;

	int ii;
	for (ii = 0; ii < 4; ii++)
		for (i = 0; i <= MAX_GENOTYPES; i++)
			for (j = 0; j <= MAX_GENOTYPES; j++)
			{
				dyad_denovo[ii][i][j] = 0;
				for (k = 0; k <= MAX_GENOTYPES; k++)
					trio_denovo[ii][i][j][k] = 0;
			}
	strcpy(ss, argv[7]);
	if ((strchr(ss, 'Y')) || (strchr(ss, 'y')))
	{
		HAPLOID = TRUE;
		max_gen = NO_ALLELES;
		min_depth_needed = 1;
		for (ii = 0; ii < 4; ii++)
			for (i = 0; i < max_gen; i++)
				for (j = 0; j < max_gen; j++)
					if (i != j)
						dyad_denovo[ii][i][j] = 1;
	}
	else
	{
		max_gen = MAX_GENOTYPES;
		min_depth_needed = 2;
		int da, db, ma, mb, ka, kb;
		for (ii = 0; ii < 4; ii++)
			for (i = 0; i < max_gen; i++)
			{
				get_het_alleles(i, &da, &db, ii);
				for (j = 0; j < max_gen; j++)
				{
					get_het_alleles(j, &ka, &kb, ii);
					if ((ka != da) && (ka != db) && (kb != da) && (kb != db))
						dyad_denovo[ii][i][j] = 1;
					for (k = 0; k < max_gen; k++)
					{
						get_het_alleles(k, &ma, &mb, ii);
						if (((ka == ma) && (kb == da)) || ((ka == ma) && (kb == db)) || ((ka == mb) && (kb == da)) || ((ka == mb) && (kb == db)) ||
							((kb == ma) && (ka == da)) || ((kb == ma) && (ka == db)) || ((kb == mb) && (ka == da)) || ((kb == mb) && (ka == db)))
							trio_denovo[ii][i][k][j] = 0;
						else if (((ka != ma) && (kb != db)) && ((kb != ma) && (ka != db)) && ((ka != mb) && (kb != da)) && ((kb != mb) && (ka != da)))
							trio_denovo[ii][i][k][j] = 2;
						else
							trio_denovo[ii][i][k][j] = 1;
					}
				}
			}
	}

	strcpy(ss, argv[9]);
	if ((strchr(ss, 'Y')) || (strchr(ss, 'y')))
		rewrite_models = TRUE;

	strcpy(ss, argv[10]);
	if ((strchr(ss, 'Y')) || (strchr(ss, 'y')))
	{
		strcpy(linkage_filename, argv[11]);
		strcpy(ss, argv[12]);
		denovo_rate = (double)atof(ss);
		if ((denovo_rate < 1e-30) || (denovo_rate > theta))
		{
			printf("\n Encounted impossible denovo mutation rate of %g with a theta of %g", denovo_rate, theta);
			exit(1);
		}
		ln_denovo = log(denovo_rate);
		use_ped = TRUE;
		if (argc == 14)
		{
			use_guide = TRUE;
			strcpy(guide_name, argv[13]);
			sprintf(ss, "%s", argv[13]);
			if ((guide_file = fopen(ss, "r")) == (FILE *)NULL)
			{
				printf("\n Can not open file %s for writing which should contain the guide_file\n", ss);
				exit(1);
			}
		}
	}
	else
	{
		use_ped = FALSE;
		ln_denovo = 0;
		if (argc == 12)
		{
			use_guide = TRUE;
			sprintf(ss, "%s", argv[11]);
			strcpy(guide_name, argv[11]);
			if ((guide_file = fopen(ss, "r")) == (FILE *)NULL)
			{
				printf("\n Can not open file %s for writing which should contain the guide_file\n", ss);
				exit(1);
			}
		}
	}

	allele_char[0] = 'A';
	allele_char[1] = 'C';
	allele_char[2] = 'G';
	allele_char[3] = 'T';
	allele_char[4] = 'D';
	allele_char[5] = 'I';
	allele_char[6] = 'N';

	pDIR = opendir(".");

	if (pDIR == NULL)
	{
		fprintf(stderr, "%s %d: opendir() failed (%s)\n",
				__FILE__, __LINE__, strerror(errno));
		exit(-1);
	}

	strcpy(sdxname, argv[2]);
	if ((sfile = fopen(sdxname, "r")) == (FILE *)NULL)
	{
		printf("\n Can not open file %s\n", sdxname);
		exit(1);
	}

	if (strstr(sdxname, ".sdx") != NULL)
	{
		for (i = strlen(sdxname) - 1; i > 0; i--)
			if (sdxname[i] == '.')
			{
				sdxname[i] = '\0';
				i = 0;
			}
	}

	if (use_guide)
		sprintf(ss, "%s.models.gz", guide_name);
	else
		sprintf(ss, "%s.models.gz", sdxname);
	if ((modelfile_read = gzopen(ss, "r")) == (gzFile)NULL)
		printf("\n Not Using Any Saved Genotype Models\n");
	else
	{
		use_saved_models = TRUE;
		printf("\n We are using the genotype models saved in %s \n", ss);
		gzbuffer(modelfile_read, 13107200);
	}

	if (rewrite_models)
	{
		if (use_saved_models)
			sprintf(ss, "%s.models.gz", argv[4]);
		else if (use_guide)
			sprintf(ss, "%s.models.gz", guide_name);
		else
			sprintf(ss, "%s.models.gz", sdxname);

		if ((modelfile_write = gzopen(ss, "w")) == (gzFile)NULL)
		{
			printf("\n Can not open %s for writing the models\n", ss);
			exit(1);
		}
	}
	if (sizeof(float) != sizeof(int))
		dump_error("\n All sorts of hell is going to break loose becasue an int and float are different sizes \n");

	if (use_saved_models || rewrite_models)
	{
		single_model_size = (MAX_GENOTYPES) * (1 + NO_ALLELES);
		model_bucket_size = sizeof(float) * single_model_size * MAX_MODELS;
		printf("\n Model bucket size is %ld \n\n", (long)model_bucket_size);
		model_bucket = (void *)malloc(model_bucket_size);
		if (!model_bucket)
			dump_error("\n Failed to allocated space for the model_buckets \n");
		memset(model_bucket, 0, model_bucket_size);
	}

	sprintf(sss, "%s.seq", sdxname);

	if ((reffile = gzopen(sss, "r")) == (gzFile)NULL)
	{
		printf("\n Can not open file %s for reading\n", sss);
		exit(1);
	}
	gzbuffer(reffile, 131072);

	printf("\n About to initialize the genome buffer \n\n");
	init_genome_buffer(reffile);

	printf("\n Finished genome buffer initialization \n\n");
	fgets(sss, 256, sfile);
	no_contigs = atoi(sss);
	frag_pos = uvector(-1, no_contigs);
	contig_names = cmatrix(0, no_contigs, 0, 256);
	chrom_type = ivector(0, no_contigs);
	frag_pos[-1] = 0;
	char *token;
	for (i = 0; i < no_contigs; i++)
	{
		// printf("\n About to read line %d \n\n",i);
		fgets(sss, 1024, sfile);
		token = strtok(sss, "\t \n");
		frag_pos[i] = atoi(token) + 15;
		frag_pos[i] += frag_pos[i - 1];
		token = strtok(NULL, "\t \n");
		strcpy(contig_names[i], token);
		strcpy(sss, contig_names[i]);
		chrom_type[i] = AUTO;
		char *token;
		token = strtok(sss, ":_- \n\0");
		char lett;
		lett = tolower(token[3]);
		if (lett == 'x')
			chrom_type[i] = CHRX;
		else if (lett == 'y')
			chrom_type[i] = CHRY;
		else if (lett == 'm')
			chrom_type[i] = CHRMT;

		// printf("\nFor contig %d %s we have offset %u\n\n",i,contig_names[i],frag_pos[i]);
	}
	fclose(sfile);

	printf("\n Finished reading the sdx file \n\n");

	no_files = atoi(argv[3]);
	filenames = cmatrix(0, no_files, 0, 128);
	pileupfile = malloc(sizeof(gzFile) * (no_files + 1));
	if (!pileupfile)
		dump_error("\n Error allocating pileupfile array \n");

	pDirEnt = readdir(pDIR);
	printf("\n Just read dir \n\n");
	i = 0;
	while (pDirEnt != NULL && i <= no_files)
	{
		if ((strstr(pDirEnt->d_name, argv[1]) != NULL))
		{
			if ((pileupfile[i] = gzopen(pDirEnt->d_name, "rb")) == (gzFile)NULL)
			{
				printf("\n Can not open file %s which should contain pileup information\n", pDirEnt->d_name);
				exit(1);
			}
			gzbuffer(pileupfile[i], 1310720);
			strcpy(sss, pDirEnt->d_name);
			token = strtok(sss, "\n.\t \0");
			strcpy(filenames[i++], token);
		}

		pDirEnt = readdir(pDIR);
	}
	closedir(pDIR);

	if (i > no_files)
		dump_error("\n Found more files than you specified \n");

	no_files = i;

	INDIV = no_files;

	sprintf(snp_type[SNP], "SNP");
	sprintf(snp_type[DELETION], "DEL");
	sprintf(snp_type[INSERTION], "INS");
	sprintf(snp_type[LOW], "LOW");
	sprintf(snp_type[MULTI], "MULTIALLELIC");
	sprintf(snp_type[MESS], "MESS");

	printf("\n Found a total of %d individuals\n\n", INDIV);

	all_samples = (SAMNODE **)malloc((unsigned)((INDIV) * sizeof(SAMNODE *)));
	if (!all_samples)
		dump_error("Allocation failure in samples\n");

	for (i = 0; i < INDIV; i++)
	{
		all_samples[i] = sample_alloc(500);
		all_samples[i]->which = i;
		all_samples[i]->no_kids = 0;
		strcpy(all_samples[i]->indiv, filenames[i]);
	}

	if (use_ped)
	{
		FILE *pedfile;
		printf("\n Reading Ped file \n");
		if ((pedfile = fopen(linkage_filename, "r")) == NULL)
		{
			printf("\n Could Not open %s", linkage_filename);
			exit(1);
		}
		fgets(sss, 81919, pedfile);
		while (!feof(pedfile) && strlen(sss) > 5)
		{
			char *fam, *ind, *token;
			fam = strtok(sss, "\n\t ");
			ind = strtok(NULL, "\n\t ");
			// printf("\n Read Pedigree data for fam=%s ind=%s",fam,ind);
			for (i = 0; i < INDIV; i++)
				if (strcmp(ind, all_samples[i]->indiv) == 0)
				{
					printf("\n\t\tFOUND %s ind", ind);
					strcpy(all_samples[i]->family, fam);
					token = strtok(NULL, "\n\t ");
					if (strcmp(token, "0") != 0)
						for (j = 0; j < INDIV; j++)
							if (strcmp(token, all_samples[j]->indiv) == 0)
							{
								all_samples[i]->dad = all_samples[j];
								all_samples[j]->kids[all_samples[j]->no_kids] = all_samples[i];
								all_samples[j]->no_kids++;
								j = INDIV;
							}
					token = strtok(NULL, "\n\t ");
					if (strcmp(token, "0") != 0)
						for (j = 0; j < INDIV; j++)
							if (strcmp(token, all_samples[j]->indiv) == 0)
							{
								all_samples[i]->mom = all_samples[j];
								all_samples[j]->kids[all_samples[j]->no_kids] = all_samples[i];
								all_samples[j]->no_kids++;
								j = INDIV;
							}
					token = strtok(NULL, "\n\t ");
					all_samples[i]->sex = (short)atoi(token);
				}
			sss[0] = '\0';
			if (!feof(pedfile))
				fgets(sss, 81919, pedfile);
		}
		fclose(pedfile);
		printf("\n Done reading pedfile \n\n");
	}

	HW_exact = NULL;
	genotype_order = NULL;
	if (!HAPLOID)
	{
		HW_exact = (double ***)malloc((unsigned)(INDIV + MAX_PRIOR_CALLS + 1) * sizeof(double **));
		if (!HW_exact)
			dump_error("Allocation failure in configs\n");
		for (i = 0; i < INDIV + MAX_PRIOR_CALLS; i++)
			HW_exact[i] = NULL;

		// printf("\n Done filling the HW matrix \n\n");

		genotype_order = imatrix(0, NO_ALLELES - 1, 0, 13);
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
	}
	else
	{
		genotype_order = imatrix(0, NO_ALLELES - 1, 0, max_gen);
		genotype_order[0][0] = 0;
		genotype_order[0][1] = 2;
		genotype_order[0][2] = 1;
		genotype_order[0][3] = 3;
		genotype_order[0][4] = 4;
		genotype_order[0][5] = 5;

		genotype_order[1][0] = 1;
		genotype_order[1][1] = 3;
		genotype_order[1][2] = 0;
		genotype_order[1][3] = 2;
		genotype_order[1][4] = 4;
		genotype_order[1][5] = 5;

		genotype_order[2][0] = 2;
		genotype_order[2][1] = 0;
		genotype_order[2][2] = 1;
		genotype_order[2][3] = 3;
		genotype_order[2][4] = 4;
		genotype_order[2][5] = 5;

		genotype_order[3][0] = 3;
		genotype_order[3][1] = 1;
		genotype_order[3][2] = 0;
		genotype_order[3][3] = 2;
		genotype_order[3][4] = 4;
		genotype_order[3][5] = 5;
	}

	for (i = 0; i < MAX_GENOTYPES; i++)
	{
		int a, b;
		for (j = 0; j < 4; j++)
		{
			for (k = 0; k < NO_ALLELES; k++)
				allele_counts[j][i][k] = 0;
			get_het_alleles(i, &a, &b, j);
			allele_counts[j][i][a]++;
			if (!HAPLOID)
				allele_counts[j][i][b]++;
		}
	}

	unsigned int *base_count;
	double *mean;
	int *median;
	int *max_coverage;
	unsigned int *tot_1x;
	unsigned int *tot_8x;
	unsigned int **counts;
	counts = umatrix(0, no_files, 0, MAX_DIST);
	mean = dvector(0, no_files);
	median = ivector(0, no_files);
	base_count = uvector(0, no_files);
	tot_1x = uvector(0, no_files);
	tot_8x = uvector(0, no_files);
	max_coverage = ivector(0, no_files);

	for (i = 0; i < no_files; i++)
	{
		mean[i] = 0.0;
		base_count[i] = 0;
		tot_1x[i] = 0;
		tot_8x[i] = 0;
		median[i] = 0;
		max_coverage[i] = 0;
		for (j = 0; j < MAX_DIST; j++)
			counts[i][j] = 0;
	}

	// Init threads
	fill_first_alpha_prior();
	printf("\n About to initialize threads and mutexes \n\n");

	pthread_mutex_init(&outfile_write_mutex, NULL);
	pthread_mutex_init(&snpfile_write_mutex, NULL);

	PTHREAD_DATA_NODE **thread_data;
	thread_data = (PTHREAD_DATA_NODE **)malloc(sizeof(PTHREAD_DATA_NODE *) * no_threads);
	if (!thread_data)
		dump_error("\n Can not allocate thread_data \n");

	int rc;
	threads = (pthread_t *)malloc(sizeof(pthread_t) * no_threads);
	if (!threads)
		dump_error("\n Could not allocate ram for threads \n");
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);

	for (i = 0; i < no_threads; i++)
	{
		// printf("\n Working on thread %d \n\n",i);
		thread_data[i] = (PTHREAD_DATA_NODE *)malloc(sizeof(PTHREAD_DATA_NODE));
		if (!thread_data[i])
			dump_error("\n Can not allocate thread_data 2 \n");
		pthread_mutex_init(&(thread_data[i]->mutex), NULL);
		pthread_mutex_lock(&(thread_data[i]->mutex));
		thread_data[i]->tid = i;
		thread_data[i]->status = DATA_EMPTY;
		thread_data[i]->chrom = AUTO;
		thread_data[i]->pos = 0;
		thread_data[i]->dom_int = NO_ALLELES;
		thread_data[i]->dom = 'N';
		thread_data[i]->reads = usmatrix(0, INDIV - 1, 0, NO_ALLELES - 1);
		thread_data[i]->HAPLOID = HAPLOID;
		thread_data[i]->times_used = 0;
		rc = pthread_create(&threads[i], &attr, call_single_base, (void *)thread_data[i]);
		if (rc)
		{
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
	}

	printf("\n Finished initializing threads and mutex \n\n");
	int data_size = (NO_ALLELES) * sizeof(unsigned short);
	for (i = 0; i < no_files; i++)
	{
		gzread(pileupfile[i], (void *)ss, pileup_header_size);
		if (strcmp(ss, argv[2]) != 0)
			printf("\n Warning %s was mapped relative to ,%s, and being called relative to %s.\nIf these aren't the same genome all sorts of awful will happen\n",
				   filenames[i], ss, argv[2]);
	}

	int which_thread = 0;
	fprintf(snpfile, "Fragment\tPosition\tReference\tAlleles\tAllele_Counts\tType");
	gzprintf(outfile, "Fragment\tPosition\tReference");
	gzprintf(pilefile, "Fragment\tPosition\tReference");
	for (i = 0; i < INDIV; i++)
	{
		fprintf(snpfile, "\t%s\t", filenames[i]);
		gzprintf(outfile, "\t%s\t", filenames[i]);
		gzprintf(pilefile, "\t%s\t\t\t\t\t", filenames[i]);
	}
	unsigned int current_pos = 0;
	which = 0;
	unsigned int last = frag_pos[no_contigs - 1];
	char *last_chrom;
	last_chrom = cvector(0, 1024);
	sprintf(last_chrom, "!!!!!!!");
	char line[4096];
	long start_bed_line, end_bed_line;
	start_bed_line = 1000;
	end_bed_line = 0;
	int max_seek = 2000000000;
	int this_model_bucket = 0;

	while (current_pos < last)
	{
		if ((use_guide) && (start_bed_line > end_bed_line))
		{
			line[0] = '\0';
			// printf("\n Got in the use_guide section \n");
			if (!feof(guide_file))
				fgets(line, 4095, guide_file);
			if (strlen(line) > 2)
			{
				char *token;
				token = strtok(line, "\t \n");
				if (strcmp(token, last_chrom) != 0)
				{
					which = -1;
					for (i = 0; i < no_contigs; i++)
						if (strcmp(token, contig_names[i]) == 0)
						{
							which = i;
							i = no_contigs;
						}
					if (which < 0)
					{
						printf("\n For line chrom %s which is not in our genome \n", token);
						exit(1);
					}
				}
				strcpy(last_chrom, token);
				token = strtok(NULL, "\t \n");
				start_bed_line = frag_pos[which - 1] + atoi(token) - 1;
				token = strtok(NULL, "\t \n");
				end_bed_line = frag_pos[which - 1] + atoi(token) - 1;
				current_pos = start_bed_line;
				long seek_pos = (long)data_size * (long)current_pos + (long)pileup_header_size;
				// long seek_pos = (long)data_size*(long)current_pos;
				long this_diff = seek_pos - current_byte_pos;
				if (this_diff < 0)
				{
					for (i = 0; i < no_files; i++)
					{
						long so_far = 0;
						long temp_long = minim((long)max_seek, seek_pos);
						unsigned int this_seek = temp_long;
						gzseek(pileupfile[i], this_seek, SEEK_SET);
						so_far += this_seek;
						while (so_far < seek_pos)
						{
							temp_long = minim((long)max_seek, seek_pos - so_far);
							this_seek = temp_long;
							gzseek(pileupfile[i], this_seek, SEEK_CUR);
							so_far += this_seek;
							// printf("\n Just moved another %d bytes for a total of %ld out of %ld",this_seek,so_far,seek_pos);
						}
					}
				}
				else
				{
					for (i = 0; i < no_files; i++)
					{
						long so_far = 0;

						while (so_far < this_diff)
						{
							long temp_long = minim((long)max_seek, this_diff - so_far);
							unsigned int this_seek = temp_long;
							gzseek(pileupfile[i], this_seek, SEEK_CUR);
							// printf("\n Just moved %d bytes",this_seek);
							so_far += this_seek;
							// printf("\n Just moved another %d bytes for a total of %ld out of %ld",this_seek,so_far,seek_pos);
						}
					}
				}
				current_byte_pos = seek_pos;
				// printf("\n About to start circling our threads with start = %u end = %u current_pos = %u last = %u \n\n",start_bed_line,end_bed_line,current_pos,last);
			}
			else
				current_pos = last;
		}
		else
			while ((current_pos > frag_pos[which]) && (which < no_contigs))
				which++;

		if (this_model_bucket >= MAX_MODELS)
		{
			if (rewrite_models)
			{
				printf("\n Saving Models \n\n");

				for (i = 0; i < no_threads; i++)
				{
					long j = 0;
					while ((thread_data[i]->status != DATA_EMPTY) && (thread_data[i]->status != DATA_ALL_DONE))
					{
						j++;
						if (j % 1000000000 == 0)
							printf("\n Waiting for Thread %d to finish before writing \n", i);
					}
					pthread_mutex_lock(&(thread_data[i]->mutex));
				}
				/* int kk = 0;
				int jj = 0;
				for(i=0;i<5;i++)
				{
				printf("\n Dumping bucket %d \n",i);
				for(ii=0;ii<MAX_GENOTYPES;ii++)
				{
					for(jj=0;jj<=NO_ALLELES;jj++)
				  printf("\t%d",(int)model_bucket[kk++]);
					printf("\n");
				}
				} */
				gzwrite(modelfile_write, model_bucket, model_bucket_size);
				for (i = 0; i < no_threads; i++)
					pthread_mutex_unlock(&(thread_data[i]->mutex));
			}
			this_model_bucket = 0;
		}
		if (current_pos < last)
		{
			// printf("\n About to gzread with no_files = %d %d bytes\n\n",no_files,data_size);
			current_byte_pos += data_size;
			start_bed_line++;

			// printf("\n About to test thread %d \n\n",which_thread);
			while (thread_data[which_thread]->status != DATA_EMPTY)
			{
				which_thread++;
				which_thread %= no_threads;
			}

			if (thread_data[which_thread]->times_used != 0)
				pthread_mutex_lock(&(thread_data[which_thread]->mutex));

			thread_data[which_thread]->times_used++;

			// printf("\n About to gzread depth data \n\n");
			for (i = 0; i < INDIV; i++)
				gzread(pileupfile[i], (void *)thread_data[which_thread]->reads[i], data_size);

			int this_b = sizeof(float) * single_model_size * this_model_bucket;
			if (use_saved_models)
			{
				gzread(modelfile_read, (model_bucket + this_b), single_model_size * sizeof(float));
				thread_data[which_thread]->model_bucket = (model_bucket + this_b);
			}
			else if (rewrite_models)
				thread_data[which_thread]->model_bucket = (model_bucket + this_b);

			char ref = get_genome(current_pos, reffile);
			thread_data[which_thread]->dom = ref;
			thread_data[which_thread]->dom_int = gen_to_int(ref);
			strcpy(thread_data[which_thread]->fragment, contig_names[which]);
			thread_data[which_thread]->pos = 1 + current_pos - frag_pos[which - 1];
			thread_data[which_thread]->chrom = chrom_type[which];
			thread_data[which_thread]->HAPLOID = HAPLOID;
			if (thread_data[which_thread]->chrom == CHRY || thread_data[which_thread]->chrom == CHRMT)
				thread_data[which_thread]->HAPLOID = TRUE;

			// printf("\n\tcontig_name = %s  pos = %d chrom = %d",contig_names[which],thread_data[which_thread]->pos,thread_data[which_thread]->chrom);

			thread_data[which_thread]->status = DATA_LOADED;
			pthread_mutex_unlock(&(thread_data[which_thread]->mutex));

			tot_bases++;
			for (i = 0; i < INDIV; i++)
			{
				unsigned short *data;
				data = thread_data[which_thread]->reads[i];
				int tot_coverage = (data[0] + data[1] + data[2] + data[3] + data[4] + data[5]);
				mean[i] += (double)(tot_coverage);
				max_coverage[i] = maxim(max_coverage[i], tot_coverage);
				j = minim(tot_coverage, MAX_DIST - 1);
				counts[i][j]++;
				base_count[i]++;
			}

			which_thread++;
			which_thread %= no_threads;
			this_model_bucket++;
			current_pos++;
		}
	}
	ALL_FINISHED = TRUE;

	long jj = 0;
	for (i = 0; i < no_threads; i++)
		while (thread_data[i]->status != DATA_ALL_DONE)
		{
			jj++;
			if (jj % 1000000000 == 0)
				printf("\n Waiting for Thread %d to finish before doing final cleanup \n", i);
		}

	for (i = 0; i < no_files; i++)
		if (base_count[i] > 0)
			mean[i] /= (double)base_count[i];

	if (rewrite_models && (this_model_bucket > 0))
	{
		printf("\n Saving Models \n\n");
		int s_size = this_model_bucket * single_model_size * sizeof(float);
		gzwrite(modelfile_write, (void *)model_bucket, s_size);
	}
	gzclose(modelfile_write);

	for (i = 0; i < no_files; i++)
	{
		for (j = 8; j < MAX_DIST; j++)
			tot_8x[i] += counts[i][j];
		tot_1x[i] = tot_8x[i] + counts[i][1] + counts[i][2] + counts[i][3] + counts[i][4] + counts[i][5] + counts[i][6] + counts[i][7];

		counts[i][0] = tot_bases - tot_1x[i];
		long median_count = counts[i][0];
		median[i] = 0;
		long stop = tot_bases / 2;
		for (j = 1; j < MAX_DIST; j++)
		{
			if (median_count > stop)
				j = MAX_DIST;
			else
				median_count += counts[i][++median[i]];
		}
	}

	fprintf(distfile, "Category");
	for (i = 0; i < no_files; i++)
		fprintf(distfile, "\t%s", filenames[i]);

	fprintf(distfile, "\nTotal Number of bases in target");
	for (i = 0; i < no_files; i++)
		fprintf(distfile, "\t%u", tot_bases);
	fprintf(distfile, "\nTotal Number of bases with at least 1x coverage");
	for (i = 0; i < no_files; i++)
		fprintf(distfile, "\t%u", tot_1x[i]);
	fprintf(distfile, "\nTotal Number of bases with at least 8x coverage");
	for (i = 0; i < no_files; i++)
		fprintf(distfile, "\t%u", tot_8x[i]);
	fprintf(distfile, "\nMean depth of coverage");
	for (i = 0; i < no_files; i++)
		fprintf(distfile, "\t%g", mean[i]);
	fprintf(distfile, "\nMedian depth of coverage");
	for (i = 0; i < no_files; i++)
		fprintf(distfile, "\t%d", median[i]);
	fprintf(distfile, "\nMaximum depth of coverage");
	for (i = 0; i < no_files; i++)
		fprintf(distfile, "\t%d", max_coverage[i]);
	fprintf(distfile, "\n\nDepth");
	for (j = 0; j < MAX_DIST - 1; j++)
	{
		fprintf(distfile, "\n%d", j);
		for (i = 0; i < no_files; i++)
			fprintf(distfile, "\t%u", counts[i][j]);
	}
	fprintf(distfile, "\n%d+", MAX_DIST - 1);
	for (i = 0; i < no_files; i++)
		fprintf(distfile, "\t%u", counts[i][MAX_DIST - 1]);
	fprintf(distfile, "\n");
	fclose(distfile);
	gzclose(outfile);
	gzclose(pilefile);

	return 0;
}

/*-------------------------------------------------------------------------------------------------------------------------------------- */
void *call_single_base(void *threadid)
{
	PTHREAD_DATA_NODE *td;
	int i, j, k, ii, jj;
	double **alpha_prior, **first_alpha_prior;
	int last_pass = 5;
	char minor[80];
	char am_count[80];
	double **d_alpha_mean, **d_alpha_var, **d_alpha_weight;
	double **saved_mean, **saved_var;
	double coef_prior[MAX_GENOTYPES];
	int ind, total_configs, issnp;
	char **outline;
	char **snpline, **pline;
	char fragment[4096];
	char sss[4096];
	CNODE **configs;
	SAMNODE **my_samples;
	int max_outlines = 1000;
	int max_snplines = 10;
	int *starting_counts;
	void *this_bucket = NULL;

	int line_length = 4096 + 20 * INDIV;

	outline = cmatrix(0, max_outlines, 0, line_length);
	snpline = cmatrix(0, max_snplines, 0, line_length);
	pline = cmatrix(0, max_snplines, 0, line_length);

	int s_line_c = 0;
	int o_line_c = 0;

	int max_configs = 512;

	k = 0;
	td = (PTHREAD_DATA_NODE *)threadid;

	// printf("\n Thread %d is launched with status %d\n\n",tid,td->status);

	saved_mean = saved_var = NULL;
	alpha_prior = dmatrix(0, MAX_GENOTYPES - 1, 0, NO_ALLELES - 1);
	first_alpha_prior = dmatrix(0, MAX_GENOTYPES - 1, 0, NO_ALLELES - 1);
	d_alpha_mean = dmatrix(0, MAX_GENOTYPES - 1, 0, NO_ALLELES - 1);
	d_alpha_var = dmatrix(0, MAX_GENOTYPES - 1, 0, NO_ALLELES - 1);
	d_alpha_weight = dmatrix(0, MAX_GENOTYPES - 1, 0, NO_ALLELES - 1);
	if (use_saved_models)
	{
		saved_mean = dmatrix(0, MAX_GENOTYPES - 1, 0, NO_ALLELES - 1);
		saved_var = dmatrix(0, MAX_GENOTYPES - 1, 0, NO_ALLELES - 1);
	}

	configs = (CNODE **)malloc((unsigned)(MAX_GENOTYPES * (max_configs + 1) + 1) * sizeof(CNODE *));
	if (!configs)
		dump_error("Allocation failure in configs\n");

	my_samples = (SAMNODE **)malloc((unsigned)(INDIV * sizeof(SAMNODE *)));
	for (j = 0; j < INDIV; j++)
	{
		my_samples[j] = sample_alloc(all_samples[j]->no_kids);
		my_samples[j]->which = all_samples[j]->which;
		my_samples[j]->sex = all_samples[j]->sex;
		my_samples[j]->no_kids = all_samples[j]->no_kids;
		strcpy(my_samples[j]->indiv, all_samples[j]->indiv);
		strcpy(my_samples[j]->family, all_samples[j]->family);
	}
	for (j = 0; j < INDIV; j++)
	{
		if (all_samples[j]->dad)
			my_samples[j]->dad = my_samples[all_samples[j]->dad->which];
		else
			my_samples[j]->dad = NULL;
		if (all_samples[j]->mom)
			my_samples[j]->mom = my_samples[all_samples[j]->mom->which];
		else
			my_samples[j]->mom = NULL;

		if (all_samples[j]->no_kids > 0)
			for (k = 0; k < all_samples[j]->no_kids; k++)
				my_samples[j]->kids[k] = my_samples[all_samples[j]->kids[k]->which];
	}

	starting_counts = ivector(0, MAX_GENOTYPES - 1);
	for (i = 0; i < MAX_GENOTYPES; i++)
		starting_counts[i] = 0;

	// long iii = 0;
	// dump_me = FALSE;
	int max_gen;
	while (!ALL_FINISHED)
	{
		pthread_mutex_lock(&(td->mutex));

		if ((td->dom_int < NO_ALLELES) && (td->status == DATA_LOADED))
		{
			td->status = DATA_RUNNING;
			for (ind = 0; ind < INDIV; ind++)
				for (j = 0; j < NO_ALLELES; j++)
					my_samples[ind]->reads[j] = td->reads[ind][j];
			if (use_saved_models || rewrite_models)
				this_bucket = td->model_bucket;
			// printf("\n Thread %d has moved into running data status \n\n",td->tid);
			int dom_int = td->dom_int;
			char dom = td->dom;
			int expos = td->pos;
			int bad_base = FALSE;
			int starting_indiv = 0;
			if (td->HAPLOID)
				max_gen = NO_ALLELES;
			else
				max_gen = MAX_GENOTYPES;
			strcpy(fragment, td->fragment);

			/* if(expos == 10517654)
			   dump_me = TRUE;
			else
			{
			   bad_base = TRUE;
			   dump_me = FALSE;
			   if(expos > 10517660)
				  exit(1);
			} */
			int chrom = td->chrom;
			SAMNODE **samples = my_samples;
			int haploid = td->HAPLOID;
			// if(samples[0]->reads[5] > 10)
			// dump_me = TRUE;
			// else
			//	dump_me = FALSE;

			td->status = DATA_EMPTY;
			pthread_mutex_unlock(&(td->mutex));

			int calls_changed = TRUE;
			int pass = 0;
			double average_depth = 0;
			double total_depth = 0;

			if (!use_saved_models)
			{
				fill_alpha_prior(alpha_prior, MAX_GENOTYPES, dom_int);
				for (ii = 0; ii < MAX_GENOTYPES; ii++)
					for (jj = 0; jj < NO_ALLELES; jj++)
						d_alpha_weight[ii][jj] = 0;
			}
			else
			{
				int kk = 0;
				// printf("\n Got here with kk = %d\n\n",kk);
				for (ii = 0; ii < MAX_GENOTYPES; ii++)
				{
					starting_counts[ii] = *((int *)this_bucket);
					kk += sizeof(float);
					starting_indiv += starting_counts[ii];
					double a0 = 0;
					for (jj = 0; jj < NO_ALLELES; jj++)
					{
						alpha_prior[ii][jj] = (double)*((float *)(this_bucket + kk));
						kk += sizeof(float);
						a0 += alpha_prior[ii][jj];
					}
					if (dump_me)
						printf("\n Got here with count %d and a0 = %g \n\n", starting_counts[ii], a0);
					double den = a0 * a0 * (a0 + 1);
					for (jj = 0; jj < NO_ALLELES; jj++)
					{
						saved_mean[ii][jj] /= a0;
						saved_var[ii][jj] = alpha_prior[ii][jj] * (a0 - alpha_prior[ii][jj]) / den;
					}
				}
			}

			int possible_calls = 0;
			for (ind = 0; ind < INDIV; ind++)
			{
				samples[ind]->tot = samples[ind]->reads[0];
				for (i = 1; i < NO_ALLELES; i++)
					samples[ind]->tot += samples[ind]->reads[i];
				if (samples[ind]->tot > 0)
				{
					for (i = 0; i < NO_ALLELES; i++)
					{
						samples[ind]->frac[i] = (double)samples[ind]->reads[i] / (double)samples[ind]->tot;
						samples[ind]->coef -= factln(samples[ind]->reads[i]);
					}
					samples[ind]->coef = factln(samples[ind]->tot);
				}
				// printf("\n Just sucked in data for individual %d  tot = %d %d %d\n\n",ind,samples[ind]->tot,samples[ind]->reads[0],samples[ind]->reads[1]);
				for (i = 0; i < max_gen; i++)
					samples[ind]->post_prob[i] = 0.0;

				if (samples[ind]->tot > min_depth_needed)
				{
					samples[ind]->initial_call = dom_int;
					samples[ind]->final_call = dom_int;
					possible_calls++;
				}
				else
				{
					samples[ind]->initial_call = MAX_GENOTYPES;
					samples[ind]->final_call = MAX_GENOTYPES;
				}
				samples[ind]->final_p = 1.0;
				// printf("\n About to make call \n\n");
				// printf("\n Called %c \n\n",int_to_gen(samples[ind]->initial_call));
			}

			// printf("\n About to allocate first config \n\n");
			// printf("\n Back from for first allocation \n\n");
			double config_threshold = 0;

			int sample_count = 0;
			for (i = 0; i < INDIV; i++)
			{
				total_depth += samples[i]->tot;
				if (samples[i]->tot >= 8)
					sample_count++;
			}

			average_depth = total_depth / (double)INDIV;
			if (average_depth < 8)
			{
				if (INDIV > 4)
					bad_base = TRUE;
				else if (!(use_saved_models))
					bad_base = TRUE;
			}

			double my_stop = maxim(1, floor(0.8 * (double)INDIV));
			if (chrom == CHRY)
				my_stop = 1;
			if (sample_count < my_stop)
				bad_base = TRUE;

			if (dump_me)
				printf("\n%s\t%d\t%c\tAverage Depth = %g\n\n", fragment, expos, dom, average_depth);

			// printf("\n About to fill prior with bad_base = %d\n\n",bad_base);

			total_configs = 0;
			if (bad_base)
			{
				for (ind = 0; ind < INDIV; ind++)
					samples[ind]->tot = 0;
				calls_changed = FALSE;
			}

			while (calls_changed && (pass < last_pass))
			{
				pass++;
				// config_threshold = minim(6.9,starting_threshold*pass);

				config_threshold = starting_threshold;

				if (dump_me)
				{
					printf("\n Alpha Matrix\n");
					for (ii = 0; ii < max_gen; ii++)
					{
						if (use_saved_models)
							printf("%d", starting_counts[ii]);
						else
							printf("%lg", d_alpha_weight[ii][0]);
						for (jj = 0; jj < NO_ALLELES; jj++)
							printf("\t%lg", alpha_prior[ii][jj]);
						printf("\n");
					}
				}

				for (ii = 0; ii < MAX_GENOTYPES; ii++)
					for (jj = 0; jj < NO_ALLELES; jj++)
						first_alpha_prior[ii][jj] = alpha_prior[ii][jj];

				fill_alpha_coef(alpha_prior, coef_prior, max_gen);
				fill_sample_like(samples, alpha_prior, max_gen, INDIV, dom_int, pass);
				if (pass == 1)
				{
					if (dump_me)
						printf("\n About to enter main while loop \n\n");
					total_configs = 1;
					configs[0] = config_alloc(INDIV, dom_int, samples, haploid, starting_counts, chrom, TRUE);
				}

				for (i = 0; i < total_configs; i++)
					fill_config_like(configs[i], samples, INDIV);
				total_configs = clean_config_probs(configs, samples, total_configs, max_configs, max_gen, INDIV, INDIV - 1, dom_int, haploid, config_threshold, starting_counts);

				if (dump_me)
					printf("\n About to loop the individuals \n\n");

				for (ind = 0; ind < INDIV; ind++)
				{
					if (samples[ind]->tot > min_depth_needed)
					{
						// printf("\n About to clean config probs with total configs = %d\n\n",total_configs);
						/* if(total_configs > max_configs)
					   total_configs = clean_config_probs(configs,samples,total_configs,max_configs,max_gen,alpha_prior,coef_prior,INDIV,ind,dom_int,haploid,config_threshold);
					   else
					   if(total_configs > 1)
					   qsort(configs,total_configs,sizeof(CNODE *),sort_configs); */
						if (dump_me)
							printf("\n About to fill config probs for pass %d indiv = %d and total_Configs = %d \n\n", pass, ind, total_configs);
						total_configs = fill_config_probs(configs, total_configs, samples, max_gen, INDIV, ind, dom_int, chrom, haploid, config_threshold, starting_counts);
						total_configs = clean_config_probs(configs, samples, total_configs, max_configs, max_gen, INDIV, INDIV - 1, dom_int, haploid, config_threshold, starting_counts);
					}
					else
					{
						samples[ind]->final_call = MAX_GENOTYPES;
						for (i = 0; i < max_gen; i++)
							samples[ind]->post_prob[i] = 0.0;
						samples[ind]->post_prob[MAX_GENOTYPES] = 1.0;
						for (i = 0; i < total_configs; i++)
							configs[i]->sample_calls[ind] = MAX_GENOTYPES;
						samples[ind]->final_p = 1.0;
					}
				}

				// total_configs = clean_config_probs(configs,samples,total_configs,max_configs,max_gen,INDIV,INDIV-1,dom_int,haploid,config_threshold,starting_counts);
				if (dump_me)
					printf("\n 2 Out of clean config probs with total configs = %d\n\n", total_configs);
				double max_post = configs[0]->post;
				double tot_post = 0;
				for (i = 0; i < total_configs; i++)
				{
					configs[i]->post -= max_post;
					if (dump_me)
						printf("\n For config %d current post = %g with like =%g and prior =%g\n", i, configs[i]->post, configs[i]->like, configs[i]->prior);
					if (configs[i]->post > -40)
						configs[i]->post = exp(configs[i]->post);
					else
						configs[i]->post = 0;
					tot_post += configs[i]->post;
					if (dump_me)
						for (ind = 0; ind < INDIV; ind++)
							printf("\n\tConfig %d is calling individual %d a %d", i, ind, configs[i]->sample_calls[ind]);
				}
				for (i = 0; i < total_configs; i++)
					configs[i]->post /= tot_post;

				for (ind = 0; ind < INDIV; ind++)
					for (i = 0; i < max_gen; i++)
						samples[ind]->post_prob[i] = 0;

				// printf("\n About to do the sample probs \n\n");
				for (ind = 0; ind < INDIV; ind++)
					if (samples[ind]->tot > min_depth_needed)
						for (i = 0; i < total_configs; i++)
							samples[ind]->post_prob[(int)configs[i]->sample_calls[ind]] += configs[i]->post;

				// printf("\n About to check if the calls changed with indiv = %d\n\n",INDIV);
				int is_variation = FALSE;
				calls_changed = FALSE;
				for (ind = 0; ind < INDIV; ind++)
					if (samples[ind]->tot > min_depth_needed)
					{
						int besti = 0;
						for (i = 1; i < max_gen; i++)
							if (samples[ind]->post_prob[i] > samples[ind]->post_prob[besti])
								besti = i;
						samples[ind]->final_p = samples[ind]->post_prob[besti];
						samples[ind]->final_call = besti;
						if (besti != dom_int)
							is_variation = TRUE;
						// if((samples[ind]->final_call != samples[ind]->initial_call) || (samples[ind]->final_p < THRESHOLD))
						if (samples[ind]->final_call != samples[ind]->initial_call)
							calls_changed = TRUE;
						if (dump_me)
							printf("\nind = %d Making call  %c at p = %g with initial call %c",
								   ind, int_to_gen(besti), samples[ind]->post_prob[besti], int_to_gen(samples[ind]->initial_call));
					}

				if (INDIV + starting_indiv < 4 || pass == last_pass)
					calls_changed = FALSE;

				if (calls_changed || rewrite_models || is_variation)
				{
					for (ii = 0; ii < max_gen; ii++)
						for (jj = 0; jj < NO_ALLELES; jj++)
							d_alpha_weight[ii][jj] = d_alpha_mean[ii][jj] = d_alpha_var[ii][jj] = (double)0.0;

					// printf("\n Calls have changed.  Recalculating alpha matrix  with total_configs = %d\n\n",total_configs);
					for (i = 0; i < total_configs; i++)
						for (ind = 0; ind < INDIV; ind++)
							if (samples[ind]->tot > min_depth_needed)
								for (j = 0; j < NO_ALLELES; j++)
								{
									// printf("\n Working on individual %d of %d and allele %d of %d which is called %d\n\n",ind,INDIV,j,NO_ALLELES,configs[i]->sample_calls[ind]);
									d_alpha_mean[(int)configs[i]->sample_calls[ind]][j] += samples[ind]->frac[j] * configs[i]->post;
									d_alpha_var[(int)configs[i]->sample_calls[ind]][j] += (samples[ind]->frac[j] * samples[ind]->frac[j]) * configs[i]->post;
									d_alpha_weight[(int)configs[i]->sample_calls[ind]][j] += configs[i]->post;
									// d_alpha_mean[(int)configs[i]->sample_calls[ind]][j] += samples[ind]->frac[j];
									// d_alpha_var[(int)configs[i]->sample_calls[ind]][j] += (samples[ind]->frac[j]*samples[ind]->frac[j]);
									// d_alpha_weight[(int)configs[i]->sample_calls[ind]][j] += 1.0;
								}
					for (ii = 0; ii < max_gen; ii++)
						for (jj = 0; jj < NO_ALLELES; jj++)
							if (d_alpha_weight[ii][jj] > 1e-9)
							{
								d_alpha_mean[ii][jj] /= d_alpha_weight[ii][jj];
								d_alpha_var[ii][jj] /= d_alpha_weight[ii][jj];
								d_alpha_var[ii][jj] -= d_alpha_mean[ii][jj] * d_alpha_mean[ii][jj];
							}
					if (use_saved_models)
					{
						for (ii = 0; ii < max_gen; ii++)
						{
							double sc = starting_counts[ii];
							if (starting_counts[ii] > 0)
								for (jj = 0; jj < NO_ALLELES; jj++)
								{
									double new_weight = d_alpha_weight[ii][jj] + sc;
									double new_mean = (d_alpha_weight[ii][jj] * d_alpha_mean[ii][jj] + sc * saved_mean[ii][jj]) / new_weight;
									double mean_var = (d_alpha_weight[ii][jj] * d_alpha_var[ii][jj] + sc * saved_var[ii][jj]) / new_weight;
									double var_mean = (d_alpha_weight[ii][jj] * d_alpha_mean[ii][jj] * d_alpha_mean[ii][jj] + sc * saved_mean[ii][jj] * saved_mean[ii][jj]) / new_weight;
									var_mean -= new_mean * new_mean;
									d_alpha_mean[ii][jj] = new_mean;
									d_alpha_weight[ii][jj] = new_weight;
									d_alpha_var[ii][jj] = maxim(var_mean, 0) + maxim(mean_var, 0);
								}
						}
					}
					// printf("\n About to finish average calcs \n\n");

					double var_eps = 1e-6;
					for (ii = 0; ii < max_gen; ii++)
					{
						int non_zero_var = 0;
						int this_min = 0;

						for (jj = 0; jj < NO_ALLELES; jj++)
						{
							if (d_alpha_weight[ii][jj] >= 1.5 && d_alpha_var[ii][jj] > var_eps * d_alpha_mean[ii][jj])
								non_zero_var++;
							if (d_alpha_mean[ii][jj] < d_alpha_mean[ii][this_min])
								this_min = jj;
							// printf("\n For ii = %d jj = %d  mean = %g  var = %g  weight = %g",ii,jj,d_alpha_mean[ii][jj],d_alpha_var[ii][jj],d_alpha_weight[ii][jj]);
						}
						if (non_zero_var > 1)
						{
							double s0 = 1.0;
							for (jj = 0; jj < NO_ALLELES; jj++)
								if (jj != this_min && d_alpha_var[ii][jj] > var_eps * d_alpha_mean[ii][jj])
									s0 *= d_alpha_mean[ii][jj] * (1.0 - d_alpha_mean[ii][jj]) / d_alpha_var[ii][jj];
							s0 = pow(s0 - 1.0, (double)1.0 / (double)(non_zero_var - 1.0));
							if (s0 > 3.0)
								for (jj = 0; jj < NO_ALLELES; jj++)
									alpha_prior[ii][jj] = maxim(0.001, d_alpha_mean[ii][jj] * s0);
							else
								for (jj = 0; jj < NO_ALLELES; jj++)
									alpha_prior[ii][jj] = first_alpha_prior[ii][jj];
						}
						else
							for (jj = 0; jj < NO_ALLELES; jj++)
								alpha_prior[ii][jj] = first_alpha_prior[ii][jj];
					}

					if (dump_me)
					{

						printf("\n At the bottom Alpha Matrix and Weights \n");
						for (ii = 0; ii < max_gen; ii++)
						{
							printf("%g", d_alpha_weight[ii][0]);
							for (jj = 0; jj < NO_ALLELES; jj++)
								printf("\t%lg", alpha_prior[ii][jj]);
							printf("\t\t\n");
						}
					}
					// printf("\n About to calculate coefficients \n\n");
					if (check_alpha_sanity(alpha_prior, d_alpha_weight, max_gen, dom_int))
					{
						for (ii = 0; ii < NO_ALLELES; ii++)
						{
							if (d_alpha_weight[ii][0] > (double)possible_calls - 0.1)
							{
								for (ind = 0; ind < INDIV; ind++)
									if (samples[ind]->tot > min_depth_needed)
									{
										samples[ind]->initial_call = ii;
										samples[ind]->final_call = ii;
										samples[ind]->final_p = 1.0;
									}
								for (jj = 0; jj < total_configs; jj++)
									config_free(configs[jj], INDIV);
								total_configs = 1;
								configs[0] = config_alloc(INDIV, dom_int, samples, haploid, starting_counts, chrom, TRUE);
								calls_changed = FALSE;
							}
						}
						if (calls_changed)
						{
							int need_to_clean = FALSE;
							for (ii = 0; ii < total_configs; ii++)
								for (jj = 0; jj < max_gen; jj++)
									if ((configs[ii]->genotype_count[jj] > 0) && (d_alpha_weight[jj][0] < 1e-6))
									{
										configs[ii]->post = -1e10;
										need_to_clean = TRUE;
									}
							if (need_to_clean)
							{
								int good_contigs = 0;
								for (ii = 0; ii < total_configs; ii++)
									if (configs[ii]->post > -1e6)
										good_contigs++;
								if (good_contigs)
									total_configs = clean_config_probs(configs, samples, total_configs, max_configs, max_gen, INDIV, INDIV - 1, dom_int, haploid, config_threshold, starting_counts);
								else
								{
									for (jj = 0; jj < total_configs; jj++)
										config_free(configs[jj], INDIV);
									fill_alpha_coef(alpha_prior, coef_prior, max_gen);
									fill_sample_like(samples, alpha_prior, max_gen, INDIV, dom_int, 0);
									total_configs = 1;
									configs[0] = config_alloc(INDIV, dom_int, samples, haploid, starting_counts, chrom, TRUE);
								}
							}
						}
					}

					for (ind = 0; ind < INDIV; ind++)
						samples[ind]->initial_call = samples[ind]->final_call;
				}
			}
			// printf("\n About to free stuff\n\n");
			if (rewrite_models)
			{
				int k = 0;
				if (dump_me)
					printf("\n Storing models into bucket position at position %ld \n", (long)td->model_bucket);

				for (i = 0; i < MAX_GENOTYPES; i++)
				{
					*((int *)(this_bucket + k)) = (int)ceil(d_alpha_weight[i][0]);
					k += sizeof(float);
					for (j = 0; j < NO_ALLELES; j++)
					{
						*((float *)(this_bucket + k)) = (float)alpha_prior[i][j];
						k += sizeof(float);
					}
				}
			}

			for (i = 0; i < total_configs; i++)
				config_free(configs[i], INDIV);
			sprintf(outline[o_line_c], "\n%s\t%d\t%c", fragment, expos, dom);
			issnp = REF_ALLELE;
			minor[0] = '\0';
			am_count[0] = '\0';
			int not_low = 0;
			int this_allele_count[NO_ALLELES];
			for (i = 0; i < NO_ALLELES; i++)
				this_allele_count[i] = 0;
			LOW_BASE = maxim(8, 0.4 * average_depth);
			int on_target = 0;
			int off_target = 0;
			for (ind = 0; ind < INDIV; ind++)
				if (samples[ind]->tot > min_depth_needed)
				{
					sprintf(sss, "\t%c\t%g", int_to_gen(samples[ind]->final_call), samples[ind]->final_p);
					strcat(outline[o_line_c], sss);
					if (samples[ind]->final_p >= THRESHOLD)
					{
						for (i = 0; i < NO_ALLELES; i++)
							if (allele_counts[dom_int][(int)samples[ind]->final_call][i])
							{
								this_allele_count[i] += allele_counts[dom_int][(int)samples[ind]->final_call][i];
								on_target += samples[ind]->reads[i];
							}
							else
							{
								if ((i != dom_int) || (samples[ind]->final_call != (NO_ALLELES - 1)))
									off_target += samples[ind]->reads[i];
							}
						if ((samples[ind]->tot > LOW_BASE) && (samples[ind]->final_call != dom_int))
							not_low++;
					}
				}
				else
				{
					strcat(outline[o_line_c], "\tN\t1");
					samples[ind]->final_call = MAX_GENOTYPES;
				}
			int this_no_alleles = 0;
			int isdel = FALSE;
			int isins = FALSE;
			// printf("\n Still working on freeing \n\n");
			for (i = 0; i < NO_ALLELES; i++)
				if (this_allele_count[i] > 0)
				{
					this_no_alleles++;
					sprintf(sss, "%c,", allele_char[i]);
					strcat(minor, sss);
					sprintf(sss, "%d,", this_allele_count[i]);
					strcat(am_count, sss);
					if (i == 4)
						isdel = TRUE;
					else if (i == 5)
						isins = TRUE;
					else if (i != dom_int)
						issnp = SNP;
				}
			if (this_no_alleles > 1 || ((this_no_alleles > 0) && (this_allele_count[dom_int] < 1)))
			{
				if ((double)off_target / (double)(on_target + off_target) > 0.15)
					issnp = MESS;
				else if (this_no_alleles > 2)
					issnp = MULTI;
				else if (not_low > 0)
					if (isdel)
						issnp = DELETION;
					else if (isins)
						issnp = INSERTION;
					else
						issnp = SNP;
				else
					issnp = LOW;
			}

			/* if(dump_me && (expos == 64078 ))
			   {
			   printf("\n Shutting off dump_me \n\n");
			   dump_me = FALSE;
			   } */

			// printf("\n About to check issnp \n\n");
			if (issnp)
			{
				minor[strlen(minor) - 1] = '\0';
				am_count[strlen(am_count) - 1] = '\0';

				if (use_ped)
				{
					int d_count = count_denovos(samples, INDIV, chrom, dom_int);
					if (d_count > 0)
						sprintf(sss, "DENOVO_%s", snp_type[issnp]);
					else
						sprintf(sss, "%s", snp_type[issnp]);
				}
				else
					sprintf(sss, "%s", snp_type[issnp]);

				sprintf(snpline[s_line_c], "\n%s\t%d\t%c\t%s\t%s\t%s", fragment, expos, dom, minor, am_count, sss);
				sprintf(pline[s_line_c], "\n%s\t%d\t%c", fragment, expos, dom);
				for (ind = 0; ind < INDIV; ind++)
				{
					sprintf(sss, "\t%c\t%g", int_to_gen(samples[ind]->final_call), samples[ind]->final_p);
					strcat(snpline[s_line_c], sss);
					for (i = 0; i < 6; i++)
					{
						sprintf(sss, "\t%d", samples[ind]->reads[i]);
						strcat(pline[s_line_c], sss);
					}
				}
				// printf("\n In thread = %d and found snp %d",tid,s_line_c);
				s_line_c++;
			}
			k = 0;
			o_line_c++;
			if (o_line_c == max_outlines)
			{
				pthread_mutex_lock(&(outfile_write_mutex));
				for (i = 0; i < o_line_c; i++)
					gzprintf(outfile, "%s", outline[i]);
				pthread_mutex_unlock(&(outfile_write_mutex));
				o_line_c = 0;
			}
			if (s_line_c == max_snplines)
			{
				pthread_mutex_lock(&(snpfile_write_mutex));
				// printf("\n Just got my mutex.  I am now writing the snp stuff \n");
				for (i = 0; i < s_line_c; i++)
				{
					gzprintf(pilefile, "%s", pline[i]);
					fprintf(snpfile, "%s", snpline[i]);
					// printf("\n%s\n%s",pline[i],snpline[i]);
				}
				pthread_mutex_unlock(&(snpfile_write_mutex));
				s_line_c = 0;
			}
		}
		else // if dom_int < NO_ALLELES
		{
			td->status = DATA_EMPTY;
			pthread_mutex_unlock(&(td->mutex));
		}
	} // While

	// printf("\n Got here 0 in thread %d \n\n",tid);

	if (o_line_c > 0)
	{
		pthread_mutex_lock(&(outfile_write_mutex));
		for (i = 0; i < o_line_c; i++)
			gzprintf(outfile, "%s", outline[i]);
		pthread_mutex_unlock(&(outfile_write_mutex));
	}
	// printf("\n Got here 1 in thread %d \n\n",tid);

	if (s_line_c > 0)
	{
		pthread_mutex_lock(&(snpfile_write_mutex));
		for (i = 0; i < s_line_c; i++)
		{
			gzprintf(pilefile, "%s", pline[i]);
			fprintf(snpfile, "%s", snpline[i]);
		}
		pthread_mutex_unlock(&(snpfile_write_mutex));
	}
	td->status = DATA_ALL_DONE;
	// printf("\n Exiting all done in thread %d \n\n",tid);
	pthread_exit(NULL);
}
/*-------------------------------------------------------------------------------------------------------------------------------------- */
int count_denovos(SAMNODE **samples, int INDIV, int chrom, int dom_int)
{
	int i;
	int d_count = 0;
	for (i = 0; i < INDIV; i++)
		if (samples[i]->final_p >= THRESHOLD)
		{
			int dad_called = MAX_GENOTYPES;
			int mom_called = MAX_GENOTYPES;
			int kid_called = samples[i]->final_call;
			if (samples[i]->dad)
				if (samples[i]->dad->final_p >= THRESHOLD)
					dad_called = samples[i]->dad->final_call;
			if (samples[i]->mom)
				if (samples[i]->mom->final_p >= THRESHOLD)
					mom_called = samples[i]->mom->final_call;
			d_count += add_denovo(kid_called, dad_called, mom_called, samples[i]->sex, chrom, dom_int);
		}
	return d_count;
}

/*-------------------------------------------------------------------------------------------------------------------------------------- */

#define buffer_chunk 50000000
static long current_start;
static long current_end;
static char *genome_buffer;

void init_genome_buffer(gzFile mfile)
{
	current_start = 0;
	current_end = buffer_chunk;
	genome_buffer = cvector(0, buffer_chunk);
	gzread(mfile, (void *)genome_buffer, sizeof(char) * buffer_chunk);
}
/*-------------------------------------------------------------------------------------------------------------------------------------- */
void free_genome_buffer(void)
{
	free_cvector(genome_buffer, 0, buffer_chunk);
}

/*-------------------------------------------------------------------------------------------------------------------------------------- */
char get_genome(unsigned int which, gzFile mfile)
{
	long current_pos;
	if (which < current_start || which >= current_end)
	{
		current_start = maxim((long)which, 0);
		current_end = current_start + buffer_chunk;
		gzseek(mfile, current_start, SEEK_SET);
		gzread(mfile, (void *)genome_buffer, sizeof(char) * buffer_chunk);
	}

	current_pos = which - current_start;
	return genome_buffer[current_pos];
}

/*-------------------------------------------------------------------------------------------------------------------------------------- */

int find_chrom(unsigned int *pos, int first, int last, int try, unsigned this)
{
	// printf("\n first = %d last = %d try  = %d this = %u pos[try] = %u",first,last,try,this,pos[try]);
	if (first == last)
		return first;
	if (first >= try)
	{
		if (this > pos[first])
			return first + 1;
		else
			return first;
	}
	if (last <= try)
		return last;

	if (pos[try] < this)
		return find_chrom(pos, try, last, (last + try) / 2, this);

	if (pos[try] > this)
		return find_chrom(pos, first, try, (try + first) / 2, this);

	return try + 1;
}

/*---------------------------------------------------------------------*/
unsigned int find_lowest(unsigned *list, int n)
{
	unsigned int low = 0;
	int i = 0;
	while (low < 1 && i < n)
		low = list[i++];

	for (; i < n; i++)
		if (list[i] > 0 && list[i] < low)
			low = list[i];

	return low;
}
/*---------------------------------------------------------------------*/

char *cvector(int nl, int nh)
{
	char *v;

	v = (char *)malloc((unsigned)(nh - nl + 1) * sizeof(char));
	if (!v)
		dump_error("allocation failure in cvector()");
	return v - nl;
}

int *ivector(int nl, int nh)
{
	int *v;

	v = (int *)malloc((unsigned)(nh - nl + 1) * sizeof(int));
	if (!v)
		dump_error("allocation failure in ivector()");
	return v - nl;
}
unsigned int *uvector(int nl, int nh)
{
	unsigned int *v;

	v = (unsigned int *)malloc((unsigned)(nh - nl + 1) * sizeof(int));
	if (!v)
		dump_error("allocation failure in uvector()");
	return v - nl;
}

double *dvector(int nl, int nh)
{
	double *v;

	v = (double *)malloc((unsigned)(nh - nl + 1) * sizeof(double));
	if (!v)
		dump_error("allocation failure in dvector()");
	return v - nl;
}

int **imatrix(int nrl, int nrh, int ncl, int nch)
{
	int i, **m;

	m = (int **)malloc((unsigned)(nrh - nrl + 1) * sizeof(int *));
	if (!m)
		dump_error("allocation failure 1 in imatrix()");
	m -= nrl;

	for (i = nrl; i <= nrh; i++)
	{
		m[i] = (int *)malloc((unsigned)(nch - ncl + 1) * sizeof(int));
		if (!m[i])
			dump_error("allocation failure 2 in imatrix()");
		m[i] -= ncl;
	}
	return m;
}

void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for (i = nrh; i >= nrl; i--)
		free((char *)(m[i] + ncl));
	free((char *)(m + nrl));
}

double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	double **m;

	m = (double **)malloc((unsigned)(nrh - nrl + 1) * sizeof(double *));
	if (!m)
		dump_error("allocation failure 1 in dmatrix()");
	m -= nrl;

	for (i = nrl; i <= nrh; i++)
	{
		m[i] = (double *)malloc((unsigned)(nch - ncl + 1) * sizeof(double));
		if (!m[i])
			dump_error("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}

unsigned short **usmatrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	unsigned short **m;

	m = (unsigned short **)malloc((unsigned)(nrh - nrl + 1) * sizeof(unsigned short *));
	if (!m)
		dump_error("allocation failure 1 in cmatrix()");
	m -= nrl;

	for (i = nrl; i <= nrh; i++)
	{
		m[i] = (unsigned short *)malloc((unsigned)(nch - ncl + 1) * sizeof(unsigned short));
		if (!m[i])
			dump_error("allocation failure 2 in cmatrix()");
		m[i] -= ncl;
	}
	return m;
}
unsigned int **umatrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	unsigned int **m;

	m = (unsigned int **)malloc((unsigned)(nrh - nrl + 1) * sizeof(unsigned int *));
	if (!m)
		dump_error("allocation failure 1 in cmatrix()");
	m -= nrl;

	for (i = nrl; i <= nrh; i++)
	{
		m[i] = (unsigned int *)malloc((unsigned)(nch - ncl + 1) * sizeof(unsigned int));
		if (!m[i])
			dump_error("allocation failure 2 in cmatrix()");
		m[i] -= ncl;
	}
	return m;
}

unsigned char **ucmatrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	unsigned char **m;

	m = (unsigned char **)malloc((unsigned)(nrh - nrl + 1) * sizeof(unsigned char *));
	if (!m)
		dump_error("allocation failure 1 in cmatrix()");
	m -= nrl;

	for (i = nrl; i <= nrh; i++)
	{
		m[i] = (unsigned char *)malloc((unsigned)(nch - ncl + 1) * sizeof(unsigned char));
		if (!m[i])
			dump_error("allocation failure 2 in cmatrix()");
		m[i] -= ncl;
	}
	return m;
}

void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for (i = nrh; i >= nrl; i--)
		free((char *)(m[i] + ncl));
	free((double *)(m + nrl));
}

char **cmatrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	char **m;

	m = (char **)malloc((unsigned)(nrh - nrl + 1) * sizeof(char *));
	if (!m)
		dump_error("allocation failure 1 in cmatrix()");
	m -= nrl;

	for (i = nrl; i <= nrh; i++)
	{
		m[i] = (char *)malloc((unsigned)(nch - ncl + 1) * sizeof(char));
		if (!m[i])
			dump_error("allocation failure 2 in cmatrix()");
		m[i] -= ncl;
	}
	return m;
}

void free_cmatrix(char **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for (i = nrh; i >= nrl; i--)
		free((char *)(m[i] + ncl));
	free((char *)(m + nrl));
}
void free_ucmatrix(unsigned char **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for (i = nrh; i >= nrl; i--)
		free((unsigned char *)(m[i] + ncl));
	free((unsigned char *)(m + nrl));
}

void free_cvector(char *v, int nl, int nh)
{
	free((char *)(v + nl));
}

void free_ivector(int *v, int nl, int nh)
{
	free((int *)(v + nl));
}
void free_uvector(unsigned int *v, int nl, int nh)
{
	free((int *)(v + nl));
}

void free_dvector(double *v, int nl, int nh)
{
	free((double *)(v + nl));
}

/*---------------------------------------------------------------------*/

int check_alpha_sanity(double **alpha_prior, double **weight, int max_gen, int ref)
{
	int i, j, fiddle;
	double frac[MAX_GENOTYPES][NO_ALLELES];
	double mag[MAX_GENOTYPES];
	double var[MAX_GENOTYPES][NO_ALLELES];
	double these_angles[MAX_GENOTYPES][MAX_GENOTYPES];
	int closest[MAX_GENOTYPES];

	fiddle = FALSE;
	for (i = 0; i < max_gen; i++)
		if (weight[i][0] > 0.001)
		{
			double tot = alpha_prior[i][0];
			mag[i] = 0.0;
			for (j = 1; j < NO_ALLELES; j++)
				tot += alpha_prior[i][j];
			double den = tot * tot * (tot + 1.0);
			for (j = 0; j < NO_ALLELES; j++)
			{
				frac[i][j] = alpha_prior[i][j] / tot;
				var[i][j] = alpha_prior[i][j] * (tot - alpha_prior[i][j]) / den;
				mag[i] += frac[i][j] * frac[i][j];
			}
			mag[i] = sqrt(mag[i]);
			double this_close = 400;
			closest[i] = 0;
			for (j = 0; j < max_gen; j++)
			{
				these_angles[i][j] = calc_angle(frac[i], mag[i], default_alpha_frac[ref][j], default_alpha_mag[ref][j], NO_ALLELES);
				if (these_angles[i][j] < this_close)
				{
					closest[i] = j;
					this_close = these_angles[i][j];
				}
			}
		}
	for (i = 0; i < NO_ALLELES; i++)
		if (weight[i][0] > 0.001)
		{
			if ((closest[i] != i) && (closest[i] < NO_ALLELES))
			{
				for (j = 0; j < NO_ALLELES; j++)
				{
					alpha_prior[i][j] = default_alpha_prior[ref][i][j];
					weight[i][j] = 0.0;
				}
			}
		}
	for (i = NO_ALLELES; i < max_gen; i++)
		if (weight[i][0] > 0.001)
			if ((closest[i] != i) || (these_angles[i][i] > 20))
			{
				// Merge these Heterozygotes into the nearest homozygote
				int closest_hom = 0;
				for (j = 0; j < NO_ALLELES; j++)
				{
					alpha_prior[i][j] = default_alpha_prior[ref][i][j];
					if (these_angles[i][j] < these_angles[i][closest_hom])
						closest_hom = j;
				}
				double new_mean[NO_ALLELES];
				double new_var[NO_ALLELES];
				double new_weight = weight[i][0] + weight[closest_hom][0];
				if (new_weight > 0.5)
				{
					for (j = 0; j < NO_ALLELES; j++)
					{
						double nm = (weight[i][0] * frac[i][j] + weight[closest_hom][0] * frac[closest_hom][j]) / new_weight;
						double vm = (weight[i][0] * frac[i][j] * frac[i][j] + weight[closest_hom][0] * frac[closest_hom][j] * frac[closest_hom][j]) / new_weight;
						vm -= nm * nm;
						double mv = (weight[i][0] * var[i][j] + weight[closest_hom][0] * var[closest_hom][j]) / new_weight;
						new_mean[j] = nm;
						new_var[j] = maxim(vm, 0.0) + maxim(mv, 0.0);
					}

					double var_eps = 1e-6;
					int non_zero_var = 0;
					int this_min = 0;
					int little_up = 0;
					for (j = 1; j < NO_ALLELES; j++)
						if (new_mean[j] > new_mean[little_up])
							little_up = j;

					for (j = 0; j < NO_ALLELES; j++)
					{
						if (new_var[j] > var_eps * new_mean[j])
							non_zero_var++;
						if (new_mean[j] < new_mean[this_min])
							this_min = j;
						if (new_mean[j] > var_eps && new_mean[j] < new_mean[little_up])
							little_up = j;
					}
					if (non_zero_var > 1)
					{
						double s0 = 1.0;
						for (j = 0; j < NO_ALLELES; j++)
							if (j != this_min)
								s0 *= new_mean[j] * (1.0 - new_mean[j]) / new_var[j];
						s0 = pow(s0 - 1.0, (double)1.0 / (double)(non_zero_var - 1.0));
						for (j = 0; j < NO_ALLELES; j++)
						{
							alpha_prior[closest_hom][j] = new_mean[j] * s0;
							weight[closest_hom][j] = new_weight;
						}
					}
				}
				weight[i][0] = 0;
				fiddle = TRUE;
			}
	return fiddle;
}
/*---------------------------------------------------------------------*/
void get_het_alleles(int i, int *a, int *b, int ref)
{
	if (i < NO_ALLELES)
	{
		*a = *b = i;
	}
	else if (i == 6)
	{
		*a = 0;
		*b = 1;
	}
	else if (i == 7)
	{
		*a = 0;
		*b = 2;
	}
	else if (i == 8)
	{
		*a = 0;
		*b = 3;
	}
	else if (i == 9)
	{
		*a = 1;
		*b = 2;
	}
	else if (i == 10)
	{
		*a = 1;
		*b = 3;
	}
	else if (i == 11)
	{
		*a = 2;
		*b = 3;
	}
	else if (i == 12)
	{
		*a = ref;
		*b = 4;
	}
	else if (i == 13)
	{
		*a = ref;
		*b = 5;
	}
	else
	{
		printf("\n This is impossible in get_het_alleles.  i = %d\n\n", i);
		exit(1);
	}

	return;
}
/*---------------------------------------------------------------------*/
int clean_config_probs(CNODE **cn, SAMNODE **sn, int n, int max, int max_gen, int indiv, int depth, int ref, int HAPLOID, double ct, int *starting)
{
	int i;

	qsort(cn, n, sizeof(CNODE *), sort_configs);
	if (dump_me)
		for (i = 0; i < n; i++)
		{
			printf("\n For configuration %d like = %g  prior = %g  post = %g", i, cn[i]->like, cn[i]->prior, cn[i]->post);
			int j;
			for (j = 0; j < max_gen; j++)
				printf(" %d", cn[i]->genotype_count[j]);
		}
	max = minim(max, n);
	for (i = 1; i < max; i++)
		if (cn[0]->post > cn[i]->post + ct)
			max = i;

	for (i = max; i < n; i++)
		config_free(cn[i], indiv);

	return max;
}
/*---------------------------------------------------------------------*/
void fill_config_like(CNODE *cn, SAMNODE **sn, int n)
{
	int i;

	cn->like = 0;

	// printf("\n IN fill hom like wiht n = %d \n\n",n);
	for (i = 0; i < n; i++)
		if (sn[i]->tot > min_depth_needed)
			cn->like += sn[i]->like[(int)cn->sample_calls[i]];
	cn->post = cn->like + cn->prior;
	// printf("\n About to leave fill_config_like %g\n\n",cn->like);
}
/*---------------------------------------------------------------------*/
int sort_configs(const void *a, const void *b)
{
	CNODE *fa, *fb;

	fa = *((CNODE **)a);
	fb = *((CNODE **)b);

	if (fa->post > fb->post)
		return -1;
	else if (fa->post < fb->post)
		return 1;
	else
		return 0;
}
/*---------------------------------------------------------------------*/
int add_denovo(int kid, int dad, int mom, int sex, int chrom, int ref)
{
	if (dad < MAX_GENOTYPES)
	{
		if (mom < MAX_GENOTYPES) // BOTH
		{
			if (chrom == AUTO)
				return trio_denovo[ref][dad][mom][kid];

			if (chrom == CHRX)
			{
				if (sex == 1)
					return dyad_denovo[ref][mom][kid];
				else
					return trio_denovo[ref][dad][mom][kid];
			}

			if (chrom == CHRY)
			{
				if (sex == 1)
					return dyad_denovo[ref][dad][kid];
				else
					return 0;
			}

			if (chrom == CHRMT)
				return dyad_denovo[ref][mom][kid];

			return 0;
		}
		else // DAD ONLY
		{
			if (chrom == AUTO)
				return dyad_denovo[ref][dad][kid];
			if ((chrom == CHRX) && (sex == 2))
				return dyad_denovo[ref][dad][kid];
			if ((chrom == CHRY) && (sex == 1))
				return dyad_denovo[ref][dad][kid];

			return 0;
		}
	}

	if (mom < MAX_GENOTYPES) // MOM ONLY
		if (chrom != CHRY)
			return dyad_denovo[ref][mom][kid];

	return 0;
}
/*---------------------------------------------------------------------*/
void fill_sample_like(SAMNODE **samples, double **alpha, int max_gen, int indiv, int dom_int, int pass)
{
	int i, j, ii, best_call, second_call, best_hom;
	double coef;
	SAMNODE *sn;
	for (i = 0; i < indiv; i++)
	{

		sn = samples[i];
		//	printf("\n Working on Sample %d",i);

		if (sn->tot > min_depth_needed)
		{
			best_call = second_call = 0;
			best_hom = 0;
			for (j = 0; j < max_gen; j++)
			{
				double tot_a = 0.0;
				sn->like[j] = 0.0;
				coef = sn->coef;
				for (ii = 0; ii < NO_ALLELES; ii++)
				{
					tot_a += alpha[j][ii];
					coef -= gammln(alpha[j][ii]);
					sn->like[j] += gammln(alpha[j][ii] + sn->reads[ii]);
				}
				coef += gammln(tot_a);
				sn->like[j] += coef;
				sn->like[j] -= gammln(sn->tot + tot_a);
				if (pass < 2)
				{
					if (sn->like[j] > sn->like[best_call])
					{
						second_call = best_call;
						best_call = j;
					}
					else
					{
						if (best_call == second_call)
							second_call = j;
						else if (sn->like[j] > sn->like[second_call])
							second_call = j;
					}
					if (j < NO_ALLELES)
						if (sn->like[j] > sn->like[best_hom])
							best_hom = j;
					// printf(" %lg",sn->like[j]);
				}
			}
			if (pass < 2)
			{
				sn->final_call = best_hom;
				if (sn->like[best_call] - sn->like[second_call] >= 2 * starting_threshold)
					sn->final_p = 0.999;
				else
					sn->initial_call = best_hom;
			}
		}
		else
		{
			sn->initial_call = MAX_GENOTYPES;
			sn->final_p = 1.0;
		}
	}
}
/*---------------------------------------------------------------------*/

int fill_config_probs(CNODE **cn, int n, SAMNODE **samples, int max_gen, int indiv, int this_depth, int ref, int chrom, int HAPLOID, double thres, int *starting)
{
	int i, j, jj, k, newcount, ii;
	CNODE *temp, *old, **new;
	double best_post = cn[0]->post;
	double best_like = cn[0]->like;
	SAMNODE *sn;
	sn = samples[this_depth];

	newcount = 0;
	new = (CNODE **)malloc((unsigned)((max_gen + 1) * (n + 1) * sizeof(CNODE *)));
	if (!new)
		dump_error("Allocation failure in fill_config_probs\n");
	if (dump_me)
		for (i = 0; i < n; i++)
		{
			printf("\nIn fill_config_probs For configuration %d like = %g  prior = %g  post = %g", i, cn[i]->like, cn[i]->prior, cn[i]->post);
			int j;
			for (j = 0; j < max_gen; j++)
				printf(" %d", cn[i]->genotype_count[j]);
		}
	// max_gen = minim(max_gen,n);
	// printf("\n Entering fill_config_probs with n = %d depth = %d\n\n",n,this_depth);
	for (i = 0; i < n; i++)
	{
		int done_it = FALSE;
		// printf("\n About to check i = %d address = %ld\n\n",i,(long)cn[i]);

		for (ii = 0; ii < i; ii++)
		{
			done_it = TRUE;
			// printf("\n Checking i = %d ii = %d n = %d \n\n",i,ii,n);
			for (jj = 0; (jj < indiv) && done_it; jj++)
			{
				// printf("\n Checking i = %d ii = %d n = %d jj = %d  this_depth = %d address of cn[i] = %ld cn[ii] = %ld",
				//	i,ii,n,jj,this_depth,(long)cn[i],(long)cn[ii]);
				// printf("\n\tcn[i]->samples = %ld cn[ii]->samples = %ld \n\n",
				//	(long)cn[i]->sample_calls,(long)cn[ii]->sample_calls);
				if ((jj != this_depth) && (cn[i]->sample_calls[jj] != cn[ii]->sample_calls[jj]))
					done_it = FALSE;
				// printf("\n Done check \n\n");
			}
			if (done_it)
				ii = i;
		}
		// printf("\n Just checked i = %d and found done_it = %d \n\n",i,done_it);
		if (!done_it)
		{

			old = cn[i];
			j = (int)old->sample_calls[this_depth];

			if (dump_me)
				printf("\n Last call was a %d \n", j);
			if (j < MAX_GENOTYPES)
			{
				// if(dump_me)
				//	printf("\n removing stuff with j = %d hets = %d homs = %d\n\n",j,old->hets,old->homs);
				for (k = 0; k < NO_ALLELES; k++)
					old->allele_count[k] -= allele_counts[ref][j][k];
				if (j >= NO_ALLELES)
					old->hets--;
				else
					old->homs--;
				if (sn->dad)
					if (sn->mom) // Both
						old->no_denovo -= add_denovo(j, (int)old->sample_calls[sn->dad->which], (int)old->sample_calls[sn->mom->which], sn->sex, chrom, ref);
					else // Dad Only
						old->no_denovo -= add_denovo(j, (int)old->sample_calls[sn->dad->which], MAX_GENOTYPES, sn->sex, chrom, ref);
				else if (sn->mom) // Mom Only
					old->no_denovo -= add_denovo(j, MAX_GENOTYPES, (int)old->sample_calls[sn->mom->which], sn->sex, chrom, ref);
				if (sn->no_kids > 0)
				{
					int kg, dg, mg;
					kg = dg = mg = MAX_GENOTYPES;
					// printf("\n Old number of denovos is %d",old->no_denovo);
					for (k = 0; k < sn->no_kids; k++)
					{
						kg = (int)old->sample_calls[sn->kids[k]->which];
						if (sn->kids[k]->dad)
							dg = (int)old->sample_calls[sn->kids[k]->dad->which];
						if (sn->kids[k]->mom)
							mg = (int)old->sample_calls[sn->kids[k]->mom->which];
						old->no_denovo -= add_denovo(kg, dg, mg, sn->kids[k]->sex, chrom, ref);
					}
					// printf("\n New number of denovos is %d",old->no_denovo);
				}
				old->like -= sn->like[j];
				old->genotype_count[j]--;
			}
			for (jj = 0; jj < max_gen; jj++)
			{
				if (!HAPLOID)
					j = genotype_order[ref][jj];
				else
					j = jj;

				// printf("\n ABout to allocate temp with ref=%d jj =%d and j = %d max = %d\n\n",ref,jj,j,max);

				// printf("\n Copying over \n\n");
				// if(dump_me)
				//{
				// printf("\nSample = %s Config = %d Genotype = %d jj = %d ref = %d this_like = %g old = %g best_like = %g best_post = %g \n\n",
				//	sn->indiv,i,j,jj,ref,sn->like[j],old->like,best_like,best_post);
				// for(k=0;k<NO_ALLELES;k++)
				//	printf("\n Allele %d has count %d",k,old->allele_count[k]);
				// }
				double templ = old->like + sn->like[j];
				// insertion / deletion weirdness.
				if (((j == 4) || (j == 12)) && (sn->reads[4] < 1))
					templ -= 1e10;
				if (((j == 13) || (j == 5)) && (sn->reads[5] < 1))
					templ -= 1e10;
				if (dump_me)
					printf("\n Just calculated genotype %d with a likelihood of %g jj = %d max = %d\n", j, templ, jj, max_gen);
				if ((templ + thres > best_post) || (templ + 0.01 > best_like))
				{
					temp = config_alloc(indiv, ref, samples, HAPLOID, starting, chrom, FALSE);
					for (k = 0; k < max_gen; k++)
						temp->genotype_count[k] = old->genotype_count[k];
					temp->like = templ;

					for (k = 0; k < indiv; k++)
						temp->sample_calls[k] = old->sample_calls[k];

					temp->genotype_count[j]++;
					temp->sample_calls[this_depth] = j;

					// printf("\n Made it in here to store the config \n");
					temp->hets = old->hets;
					temp->homs = old->homs;
					temp->no_alleles = 0;
					for (k = 0; k < NO_ALLELES; k++)
					{
						temp->allele_count[k] = old->allele_count[k];
						temp->allele_count[k] += allele_counts[ref][j][k];

						if (temp->allele_count[k] > 0)
							temp->no_alleles++;
					}
					temp->no_denovo = old->no_denovo;
					if (dump_me)
					{
						printf("\nOn Individual %s and call %d", sn->indiv, j);
						printf("\n\tBefore the number of denovo is %d hets = %d  homs = %d\n\n", temp->no_denovo, temp->hets, temp->homs);
					}
					if (j >= NO_ALLELES)
						temp->hets++;
					else
						temp->homs++;
					if (use_ped)
					{
						if (sn->dad)
							if (sn->mom) // Both
								temp->no_denovo += add_denovo(j, (int)temp->sample_calls[sn->dad->which], (int)temp->sample_calls[sn->mom->which], sn->sex, chrom, ref);
							else // Dad Only
								temp->no_denovo += add_denovo(j, (int)temp->sample_calls[sn->dad->which], MAX_GENOTYPES, sn->sex, chrom, ref);
						else if (sn->mom) // Mom Only
							temp->no_denovo += add_denovo(j, MAX_GENOTYPES, (int)temp->sample_calls[sn->mom->which], sn->sex, chrom, ref);
						if (sn->no_kids > 0)
						{
							int kg, dg, mg;
							kg = dg = mg = MAX_GENOTYPES;
							for (k = 0; k < sn->no_kids; k++)
							{
								kg = (int)temp->sample_calls[sn->kids[k]->which];
								if (sn->kids[k]->dad)
									dg = (int)temp->sample_calls[sn->kids[k]->dad->which];
								if (sn->kids[k]->mom)
									mg = (int)temp->sample_calls[sn->kids[k]->mom->which];
								temp->no_denovo += add_denovo(kg, dg, mg, sn->kids[k]->sex, chrom, ref);
							}
						}
					}
					// printf("\n\tAfterwards the number of denovo is %d hets = %d  homs = %d\n\n",temp->no_denovo,temp->hets,temp->homs);
					if (dump_me)
					{
						printf("\n\tAfterwards the number of denovo is %d hets = %d  homs = %d alleles = %d\n\n", temp->no_denovo, temp->hets, temp->homs, temp->no_alleles);
						for (k = 0; k < NO_ALLELES; k++)
							printf("\n Allele %d has count %d", k, old->allele_count[k]);
					}

					if (temp->no_alleles > 1)
						fill_prior(temp, HAPLOID);
					else
						temp->prior = 0;
					if (dump_me)
						printf("\n Back from fill prior with prior = %g \n\n", temp->prior);
					temp->post = temp->prior + temp->like;
					best_like = maxim(temp->like, best_like);
					best_post = maxim(temp->post, best_post);
					if (dump_me)
					{
						printf("\nSample data is");
						for (ii = 0; ii < NO_ALLELES; ii++)
							printf(" %d", sn->reads[ii]);
						printf("\nWith Prior is %g", temp->prior);
						printf("\n\tTotal Like = %g  This Prob = %g  Total Post = %g", temp->like, sn->like[j], temp->post);
						printf("\n call = %d  max = %d \n", j, max_gen);
					}
					if (temp->post + thres > best_post)
					{
						// printf("\n About to store \n\n");
						// printf("\n j = %d with n = %d newcount = %d Configuration:",j,n,newcount);
						// for(ii=0;ii<max_gen;ii++)
						//	printf(" %d",temp->genotype_count[ii]);
						new[newcount] = temp;
						newcount++;
					}
					else
						/* if(temp->like + 0.01 > best_like)
						{
						  // printf("\n About to store \n\n");
							  // printf("\n j = %d with n = %d newcount = %d Configuration:",j,n,newcount);
							  // for(ii=0;ii<max_gen;ii++)
							// printf(" %d",temp->genotype_count[ii]);
						  new[newcount] = temp;
						  newcount++;
						}
						else */
						config_free(temp, indiv);
				}
			}
			// printf("\n i = %d  about to free old = %d\n\n",i,old);
		}
	}
	for (i = 0; i < n; i++)
		config_free(cn[i], indiv);
	// printf("\n About to leave with newcount = %d\n\n",newcount);
	for (i = 0; i < newcount; i++)
		cn[i] = new[i];
	free(new);

	// printf("\n About to leave with newcount = %d\n\n",newcount);

	return newcount;
}
/*---------------------------------------------------------------------*/
double get_HW_exact(int i, int j, int k)
{
	int asize;
	if (!HW_exact[i])
	{
		asize = 2 * (i);
		HW_exact[i] = dmatrix(0, asize, 0, i);
		fill_hardy_weinberg(HW_exact[i], asize, i);
	}
	return HW_exact[i][j][k];
}
/*---------------------------------------------------------------------*/
void fill_prior(CNODE *temp, int HAPLOID)
{
	temp->prior = 0;
	int k;
	if (temp->no_alleles > 1)
		temp->prior = (temp->no_alleles - 1) * ln_theta;

	// printf("\nWith Theta Prior is %g",temp->prior);
	if (temp->no_denovo > 0)
		temp->prior += temp->no_denovo * ln_denovo;
	// printf("\nWith denovo Prior is %g",temp->prior);

	if (!HAPLOID && temp->no_alleles > 1)
	{

		int major = 0;
		int minor = 0;
		for (k = 1; k < NO_ALLELES; k++)
			if (temp->allele_count[k] > temp->allele_count[major])
				major = k;

		for (k = 0; k < NO_ALLELES; k++)
			if (k != major)
				minor += temp->allele_count[k];
		major = temp->allele_count[major];
		// printf("\n Major = %d \n\n",major);
		if (minor > major)
		{
			int ii = major;
			major = minor;
			minor = ii;
		}
		int hets = minim(minor, temp->hets);
		int tot_n = (minor + major) / 2;
		if ((minor - hets) % 2 == 1)
		{
			minor++;
			major++;
			// Bad mojo
		}
		// printf("\n Using minor = %d major = %d hets = %d tot_n = %d\n\n",minor,major,hets,tot_n);
		temp->prior += get_HW_exact(tot_n, minor, hets);
	}
}
/*---------------------------------------------------------------------*/
void fill_hardy_weinberg(double **exact_HW, int asize, int n)
{
	double **marg, sum, p;
	int i, j, naa, nab, nbb, Na, Nb, start, expect;

	// printf("\n Entering fill_hardy_weinberg with asize = %d\n\n",asize);
	marg = dmatrix(0, asize, 0, n);
	for (i = 0; i <= asize; i++)
		for (j = 0; j <= n; j++)
			exact_HW[i][j] = marg[i][j] = 0.0;

	for (i = 1; i <= asize; i++)
	{
		Na = 2 * n - i;
		Nb = i;
		p = (double)i / (double)(Na + Nb);
		expect = ceil(i * (1.0 - p));

		if (i % 2 == 0)
		{
			if (expect % 2 == 1)
				start = expect - 1;
			else
				start = expect;
		}
		else
		{
			if (expect % 2 == 1)
				start = expect;
			else
				start = expect - 1;
		}
		// printf("\ni=%d expect = %d  start = %d p = %g ",i,expect,start,p);
		sum = marg[i][start] = 1.0;

		nbb = ((Nb - start) / 2);
		naa = ((Na - start) / 2);

		// printf("\n naa = %d nbb = %d nab = %d  Nb = %d  Na = %d",naa,nbb,start+2,Nb,Na);
		for (nab = start + 2; naa > 0 && nbb > 0; nab += 2, naa--, nbb--)
		{
			marg[i][nab] = marg[i][nab - 2] * 4.0 * ((double)naa * (double)nbb) / ((double)(nab) * (double)(nab - 1.0));
			// printf("\n nab = %d last = %g  current = %g",nab,marg[i][nab-2],marg[i][nab]);
			sum += marg[i][nab];
		}

		nbb = ((Nb - start) / 2);
		naa = ((Na - start) / 2);

		for (nab = start - 2; nab >= 0; nab -= 2, naa++, nbb++)
		{
			marg[i][nab] = marg[i][nab + 2] * ((double)(nab + 2.0) * (double)(nab + 1.0)) / ((double)4.0 * ((double)(naa + 1.0) * (nbb + 1.0)));
			sum += marg[i][nab];
		}

		for (j = 0; j <= n; j++)
			marg[i][j] /= sum;
	}
	for (i = 0; i <= asize; i++)
		for (j = 0; j <= n; j++)
		{
			if (marg[i][j] > 1e-50)
				exact_HW[i][j] = log(marg[i][j]);
			else
				exact_HW[i][j] = -5000;

			// printf("\n i = %d  j = %d  p = %g",i,j,exact_HW[i][j]);
		}
	// exit(1);
	free_dmatrix(marg, 0, asize, 0, n);
}

/*---------------------------------------------------------------------*/
int gen_to_int(char c)
{
	if (c == 'A')
		return 0;
	if (c == 'C')
		return 1;
	if (c == 'G')
		return 2;
	if (c == 'T')
		return 3;
	if (c == 'D')
		return 4;
	if (c == 'I')
		return 5;
	if (c == 'M')
		return 6;
	if (c == 'R')
		return 7;
	if (c == 'W')
		return 8;
	if (c == 'S')
		return 9;
	if (c == 'Y')
		return 10;
	if (c == 'K')
		return 11;
	if (c == 'E')
		return 12;
	if (c == 'H')
		return 13;
	if (c == 'N')
		return 14;

	printf("\n This is impossible\n Illegal character in gen_to_int %c\n\n", c);
	exit(1);
	return -1;
}
/*---------------------------------------------------------------------*/
char int_to_gen(int c)
{
	if (c == 0)
		return 'A';
	if (c == 1)
		return 'C';
	if (c == 2)
		return 'G';
	if (c == 3)
		return 'T';
	if (c == 4)
		return 'D';
	if (c == 5)
		return 'I';
	if (c == 6)
		return 'M';
	if (c == 7)
		return 'R';
	if (c == 8)
		return 'W';
	if (c == 9)
		return 'S';
	if (c == 10)
		return 'Y';
	if (c == 11)
		return 'K';
	if (c == 12)
		return 'E';
	if (c == 13)
		return 'H';

	return 'N';
}

/*---------------------------------------------------------------------*/

SAMNODE *sample_alloc(int kids)
{
	SAMNODE *tn;
	int i;

	tn = (SAMNODE *)malloc((unsigned)sizeof(struct sample_node));
	if (!tn)
		dump_error("allocation failure in sample_alloc()");

	for (i = 0; i <= MAX_GENOTYPES; i++)
	{
		tn->post_prob[i] = 0.0;
		tn->like[i] = 0.0;
	}

	for (i = 0; i < NO_ALLELES; i++)
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
	tn->no_kids = kids;
	if (kids > 0)
	{
		tn->kids = (SAMNODE **)malloc((unsigned)sizeof(SAMNODE *) * kids);
		if (!tn->kids)
			dump_error("Can't allocate space for kids");
	}
	else
		tn->kids = NULL;

	return tn;
}

/*---------------------------------------------------------------------*/

CNODE *config_alloc(int N, int dom, SAMNODE **sn, int is_haploid, int *starting, int chrom, int first_config)
{
	CNODE *tn;
	int i, j, max_gen;

	tn = (CNODE *)malloc((unsigned)sizeof(struct config_node));
	if (!tn)
		dump_error("allocation failure in config_alloc()");
	tn->homs = 0;
	tn->hets = 0;
	if (is_haploid)
		max_gen = NO_ALLELES;
	else
		max_gen = MAX_GENOTYPES;
	tn->sample_calls = cvector(0, N - 1);
	tn->like = 0;
	tn->post = 1;
	tn->no_denovo = 0;
	tn->homs = 0;
	tn->hets = 0;
	tn->no_alleles = 0;

	if (first_config)
	{
		for (j = 0; j < NO_ALLELES; j++)
			tn->allele_count[j] = 0;
		for (i = 0; i < max_gen; i++)
		{
			tn->genotype_count[i] = starting[i];
			for (j = 0; j < NO_ALLELES; j++)
				tn->allele_count[j] += starting[i] * allele_counts[dom][i][j];
			if (i < NO_ALLELES)
				tn->homs += starting[i];
			else
				tn->hets += starting[i];
		}
		for (i = 0; i < N; i++)
			if (sn[i]->tot > min_depth_needed)
			{
				tn->sample_calls[i] = sn[i]->final_call;
				tn->genotype_count[(int)sn[i]->final_call]++;
				for (j = 0; j < NO_ALLELES; j++)
					tn->allele_count[j] += allele_counts[dom][(int)sn[i]->final_call][j];
				if (sn[i]->final_call < MAX_GENOTYPES)
				{
					if (sn[i]->final_call < NO_ALLELES)
						tn->homs++;
					else
						tn->hets++;
				}
			}
			else
				tn->sample_calls[i] = (char)MAX_GENOTYPES;
		tn->no_alleles = 0;
		for (i = 0; i < NO_ALLELES; i++)
			if (tn->allele_count[i] > 0)
				tn->no_alleles++;
		if (use_ped && (tn->no_alleles > 1))
			tn->no_denovo = count_denovos(sn, N, chrom, dom);
	}
	if (dump_me)
		printf("\n Found %d homs %d hets %d alleles\n\n", tn->homs, tn->hets, tn->no_alleles);

	if (tn->no_alleles > 1)
		fill_prior(tn, is_haploid);
	else
		tn->prior = 0;

	if (dump_me)
		printf("\n And prior = %g \n\n", tn->prior);

	return tn;
}
/*---------------------------------------------------------------------*/
void config_free(CNODE *tn, int N)
{
	// printf("\n In free with tn = %d\n\n",tn);
	// printf("\n\t tn->sample_calls =  %d  tn->alpha = %d\n\n",tn->sample_calls,tn->alpha);
	if (!tn)
		return;

	free_cvector(tn->sample_calls, 0, N);
	free(tn);
}
/*---------------------------------------------------------------------*/
void fill_first_alpha_prior(void)
{
	int i, j, k;
	// printf("\n About to fill the first prior \n\n");
	default_alpha_prior = (double ***)malloc((sizeof(double **)) * NO_ALLELES);
	if (!default_alpha_prior)
		dump_error("\n Failure to allocate default_alpha_prior \n");
	default_alpha_frac = (double ***)malloc((sizeof(double **)) * NO_ALLELES);
	if (!default_alpha_frac)
		dump_error("\n Failure to allocate default_alpha_frac \n");
	default_alpha_mag = dmatrix(0, NO_ALLELES - 1, 0, MAX_GENOTYPES - 1);
	for (i = 0; i < NO_ALLELES; i++)
	{
		default_alpha_prior[i] = dmatrix(0, MAX_GENOTYPES - 1, 0, NO_ALLELES - 1);
		default_alpha_frac[i] = dmatrix(0, MAX_GENOTYPES - 1, 0, NO_ALLELES - 1);
		for (j = 0; j < MAX_GENOTYPES; j++)
			for (k = 0; k < NO_ALLELES; k++)
				default_alpha_prior[i][j][k] = 0.0;
	}
	default_alpha_prior[0][0][0] = 1398.75;
	default_alpha_prior[0][0][1] = 2.28014;
	default_alpha_prior[0][0][2] = 1.61355;
	default_alpha_prior[0][0][3] = 1.13666;
	default_alpha_prior[0][0][4] = 3.1617;
	default_alpha_prior[0][0][5] = 1.27181;
	default_alpha_prior[0][1][0] = 0.6675;
	default_alpha_prior[0][1][1] = 212.678;
	default_alpha_prior[0][1][2] = 0.329486;
	default_alpha_prior[0][1][3] = 0.305412;
	default_alpha_prior[0][1][4] = 0.497083;
	default_alpha_prior[0][1][5] = 0.26123;
	default_alpha_prior[0][2][0] = 0.905415;
	default_alpha_prior[0][2][1] = 0.441742;
	default_alpha_prior[0][2][2] = 413.915;
	default_alpha_prior[0][2][3] = 0.62365;
	default_alpha_prior[0][2][4] = 0.660203;
	default_alpha_prior[0][2][5] = 0.359421;
	default_alpha_prior[0][3][0] = 0.758859;
	default_alpha_prior[0][3][1] = 0.471962;
	default_alpha_prior[0][3][2] = 0.706239;
	default_alpha_prior[0][3][3] = 210.849;
	default_alpha_prior[0][3][4] = 1.05241;
	default_alpha_prior[0][3][5] = 0.276823;
	default_alpha_prior[0][4][0] = 13.2677;
	default_alpha_prior[0][4][1] = 0.680801;
	default_alpha_prior[0][4][2] = 0.867212;
	default_alpha_prior[0][4][3] = 0.982356;
	default_alpha_prior[0][4][4] = 111.32;
	default_alpha_prior[0][4][5] = 0.28733;
	default_alpha_prior[0][5][0] = 13.516;
	default_alpha_prior[0][5][1] = 0.230312;
	default_alpha_prior[0][5][2] = 0.25831;
	default_alpha_prior[0][5][3] = 0.272581;
	default_alpha_prior[0][5][4] = 0.153863;
	default_alpha_prior[0][5][5] = 10.0295;
	default_alpha_prior[0][6][0] = 436.259;
	default_alpha_prior[0][6][1] = 422.251;
	default_alpha_prior[0][6][2] = 1.36162;
	default_alpha_prior[0][6][3] = 1.35467;
	default_alpha_prior[0][6][4] = 2.34361;
	default_alpha_prior[0][6][5] = 1.65746;
	default_alpha_prior[0][7][0] = 999.478;
	default_alpha_prior[0][7][1] = 3.351;
	default_alpha_prior[0][7][2] = 971.641;
	default_alpha_prior[0][7][3] = 2.8878;
	default_alpha_prior[0][7][4] = 4.04716;
	default_alpha_prior[0][7][5] = 2.90267;
	default_alpha_prior[0][8][0] = 401.747;
	default_alpha_prior[0][8][1] = 1.66616;
	default_alpha_prior[0][8][2] = 1.6901;
	default_alpha_prior[0][8][3] = 387.222;
	default_alpha_prior[0][8][4] = 3.43448;
	default_alpha_prior[0][8][5] = 2.23022;
	default_alpha_prior[0][9][0] = 0.497367;
	default_alpha_prior[0][9][1] = 78.6425;
	default_alpha_prior[0][9][2] = 75.3493;
	default_alpha_prior[0][9][3] = 0.508448;
	default_alpha_prior[0][9][4] = 0.726492;
	default_alpha_prior[0][9][5] = 0.387046;
	default_alpha_prior[0][10][0] = 0.54991;
	default_alpha_prior[0][10][1] = 93.6725;
	default_alpha_prior[0][10][2] = 0.54011;
	default_alpha_prior[0][10][3] = 94.9624;
	default_alpha_prior[0][10][4] = 0.956855;
	default_alpha_prior[0][10][5] = 0.664605;
	default_alpha_prior[0][11][0] = 0.573645;
	default_alpha_prior[0][11][1] = 0.656561;
	default_alpha_prior[0][11][2] = 70.4953;
	default_alpha_prior[0][11][3] = 77.2136;
	default_alpha_prior[0][11][4] = 1.80456;
	default_alpha_prior[0][11][5] = 0.473358;
	default_alpha_prior[0][12][0] = 48.8631;
	default_alpha_prior[0][12][1] = 0.366134;
	default_alpha_prior[0][12][2] = 0.486298;
	default_alpha_prior[0][12][3] = 0.608978;
	default_alpha_prior[0][12][4] = 37.6816;
	default_alpha_prior[0][12][5] = 0.177597;
	default_alpha_prior[0][13][0] = 19.7941;
	default_alpha_prior[0][13][1] = 0.184073;
	default_alpha_prior[0][13][2] = 0.172082;
	default_alpha_prior[0][13][3] = 0.167972;
	default_alpha_prior[0][13][4] = 0.0978209;
	default_alpha_prior[0][13][5] = 6.72615;
	default_alpha_prior[1][0][0] = 141.841;
	default_alpha_prior[1][0][1] = 0.5191;
	default_alpha_prior[1][0][2] = 0.340512;
	default_alpha_prior[1][0][3] = 0.291512;
	default_alpha_prior[1][0][4] = 0.500353;
	default_alpha_prior[1][0][5] = 0.148894;
	default_alpha_prior[1][1][0] = 1.57741;
	default_alpha_prior[1][1][1] = 1360.43;
	default_alpha_prior[1][1][2] = 0.914483;
	default_alpha_prior[1][1][3] = 1.0835;
	default_alpha_prior[1][1][4] = 1.08013;
	default_alpha_prior[1][1][5] = 1.94094;
	default_alpha_prior[1][2][0] = 0.363248;
	default_alpha_prior[1][2][1] = 0.49266;
	default_alpha_prior[1][2][2] = 216.89;
	default_alpha_prior[1][2][3] = 0.545921;
	default_alpha_prior[1][2][4] = 0.481726;
	default_alpha_prior[1][2][5] = 0.229392;
	default_alpha_prior[1][3][0] = 0.70512;
	default_alpha_prior[1][3][1] = 1.55639;
	default_alpha_prior[1][3][2] = 1.266;
	default_alpha_prior[1][3][3] = 583.963;
	default_alpha_prior[1][3][4] = 1.1176;
	default_alpha_prior[1][3][5] = 0.53722;
	default_alpha_prior[1][4][0] = 0.615346;
	default_alpha_prior[1][4][1] = 5.9153;
	default_alpha_prior[1][4][2] = 0.436753;
	default_alpha_prior[1][4][3] = 0.960804;
	default_alpha_prior[1][4][4] = 65.3034;
	default_alpha_prior[1][4][5] = 0.166339;
	default_alpha_prior[1][5][0] = 0.801274;
	default_alpha_prior[1][5][1] = 31.1124;
	default_alpha_prior[1][5][2] = 0.276082;
	default_alpha_prior[1][5][3] = 0.57533;
	default_alpha_prior[1][5][4] = 0.356165;
	default_alpha_prior[1][5][5] = 21.4554;
	default_alpha_prior[1][6][0] = 468.53;
	default_alpha_prior[1][6][1] = 482.655;
	default_alpha_prior[1][6][2] = 1.58249;
	default_alpha_prior[1][6][3] = 1.55114;
	default_alpha_prior[1][6][4] = 2.78027;
	default_alpha_prior[1][6][5] = 1.72228;
	default_alpha_prior[1][7][0] = 61.8246;
	default_alpha_prior[1][7][1] = 0.35774;
	default_alpha_prior[1][7][2] = 63.7382;
	default_alpha_prior[1][7][3] = 0.259355;
	default_alpha_prior[1][7][4] = 0.49329;
	default_alpha_prior[1][7][5] = 0.287305;
	default_alpha_prior[1][8][0] = 57.5114;
	default_alpha_prior[1][8][1] = 0.568707;
	default_alpha_prior[1][8][2] = 0.408681;
	default_alpha_prior[1][8][3] = 54.8673;
	default_alpha_prior[1][8][4] = 1.19697;
	default_alpha_prior[1][8][5] = 0.500171;
	default_alpha_prior[1][9][0] = 1.97583;
	default_alpha_prior[1][9][1] = 646.431;
	default_alpha_prior[1][9][2] = 630.806;
	default_alpha_prior[1][9][3] = 1.97612;
	default_alpha_prior[1][9][4] = 2.73376;
	default_alpha_prior[1][9][5] = 1.9608;
	default_alpha_prior[1][10][0] = 3.91379;
	default_alpha_prior[1][10][1] = 1264.44;
	default_alpha_prior[1][10][2] = 4.6037;
	default_alpha_prior[1][10][3] = 1231.98;
	default_alpha_prior[1][10][4] = 5.71963;
	default_alpha_prior[1][10][5] = 3.80692;
	default_alpha_prior[1][11][0] = 0.54446;
	default_alpha_prior[1][11][1] = 0.962782;
	default_alpha_prior[1][11][2] = 96.8566;
	default_alpha_prior[1][11][3] = 98.935;
	default_alpha_prior[1][11][4] = 1.70665;
	default_alpha_prior[1][11][5] = 0.604796;
	default_alpha_prior[1][12][0] = 0.523245;
	default_alpha_prior[1][12][1] = 30.9318;
	default_alpha_prior[1][12][2] = 0.341623;
	default_alpha_prior[1][12][3] = 0.766954;
	default_alpha_prior[1][12][4] = 24.6313;
	default_alpha_prior[1][12][5] = 0.120901;
	default_alpha_prior[1][13][0] = 0.217485;
	default_alpha_prior[1][13][1] = 28.7504;
	default_alpha_prior[1][13][2] = 0.175558;
	default_alpha_prior[1][13][3] = 0.234377;
	default_alpha_prior[1][13][4] = 0.141493;
	default_alpha_prior[1][13][5] = 9.0125;
	default_alpha_prior[2][0][0] = 468.581;
	default_alpha_prior[2][0][1] = 1.02699;
	default_alpha_prior[2][0][2] = 1.35759;
	default_alpha_prior[2][0][3] = 0.590319;
	default_alpha_prior[2][0][4] = 0.861131;
	default_alpha_prior[2][0][5] = 0.432728;
	default_alpha_prior[2][1][0] = 0.390813;
	default_alpha_prior[2][1][1] = 151.609;
	default_alpha_prior[2][1][2] = 0.302841;
	default_alpha_prior[2][1][3] = 0.241061;
	default_alpha_prior[2][1][4] = 0.29978;
	default_alpha_prior[2][1][5] = 0.210987;
	default_alpha_prior[2][2][0] = 1.07296;
	default_alpha_prior[2][2][1] = 0.856751;
	default_alpha_prior[2][2][2] = 1300.75;
	default_alpha_prior[2][2][3] = 1.15811;
	default_alpha_prior[2][2][4] = 0.920434;
	default_alpha_prior[2][2][5] = 1.06854;
	default_alpha_prior[2][3][0] = 0.386628;
	default_alpha_prior[2][3][1] = 0.372895;
	default_alpha_prior[2][3][2] = 0.719567;
	default_alpha_prior[2][3][3] = 194.764;
	default_alpha_prior[2][3][4] = 0.682238;
	default_alpha_prior[2][3][5] = 0.192583;
	default_alpha_prior[2][4][0] = 1.039;
	default_alpha_prior[2][4][1] = 0.40001;
	default_alpha_prior[2][4][2] = 5.8913;
	default_alpha_prior[2][4][3] = 0.540887;
	default_alpha_prior[2][4][4] = 64.4763;
	default_alpha_prior[2][4][5] = 0.157423;
	default_alpha_prior[2][5][0] = 0.400904;
	default_alpha_prior[2][5][1] = 0.192724;
	default_alpha_prior[2][5][2] = 12.1472;
	default_alpha_prior[2][5][3] = 0.21622;
	default_alpha_prior[2][5][4] = 0.130075;
	default_alpha_prior[2][5][5] = 8.90533;
	default_alpha_prior[2][6][0] = 64.5065;
	default_alpha_prior[2][6][1] = 65.0027;
	default_alpha_prior[2][6][2] = 0.629015;
	default_alpha_prior[2][6][3] = 0.385327;
	default_alpha_prior[2][6][4] = 0.910418;
	default_alpha_prior[2][6][5] = 0.403687;
	default_alpha_prior[2][7][0] = 1337.99;
	default_alpha_prior[2][7][1] = 4.9897;
	default_alpha_prior[2][7][2] = 1372.09;
	default_alpha_prior[2][7][3] = 4.29543;
	default_alpha_prior[2][7][4] = 6.20526;
	default_alpha_prior[2][7][5] = 3.9931;
	default_alpha_prior[2][8][0] = 56.2898;
	default_alpha_prior[2][8][1] = 0.372861;
	default_alpha_prior[2][8][2] = 0.474429;
	default_alpha_prior[2][8][3] = 56.2369;
	default_alpha_prior[2][8][4] = 1.2012;
	default_alpha_prior[2][8][5] = 0.416838;
	default_alpha_prior[2][9][0] = 1.98084;
	default_alpha_prior[2][9][1] = 653.983;
	default_alpha_prior[2][9][2] = 670.932;
	default_alpha_prior[2][9][3] = 1.97637;
	default_alpha_prior[2][9][4] = 2.78502;
	default_alpha_prior[2][9][5] = 1.84321;
	default_alpha_prior[2][10][0] = 0.543888;
	default_alpha_prior[2][10][1] = 96.925;
	default_alpha_prior[2][10][2] = 0.696742;
	default_alpha_prior[2][10][3] = 95.2395;
	default_alpha_prior[2][10][4] = 1.26314;
	default_alpha_prior[2][10][5] = 0.509389;
	default_alpha_prior[2][11][0] = 1.3535;
	default_alpha_prior[2][11][1] = 1.37881;
	default_alpha_prior[2][11][2] = 427.984;
	default_alpha_prior[2][11][3] = 415.578;
	default_alpha_prior[2][11][4] = 2.40342;
	default_alpha_prior[2][11][5] = 1.57879;
	default_alpha_prior[2][12][0] = 0.94697;
	default_alpha_prior[2][12][1] = 0.302139;
	default_alpha_prior[2][12][2] = 35.8007;
	default_alpha_prior[2][12][3] = 0.533744;
	default_alpha_prior[2][12][4] = 28.3782;
	default_alpha_prior[2][12][5] = 0.146126;
	default_alpha_prior[2][13][0] = 0.181626;
	default_alpha_prior[2][13][1] = 0.165676;
	default_alpha_prior[2][13][2] = 17.8087;
	default_alpha_prior[2][13][3] = 0.153894;
	default_alpha_prior[2][13][4] = 0.0839989;
	default_alpha_prior[2][13][5] = 6.03531;
	default_alpha_prior[3][0][0] = 240.883;
	default_alpha_prior[3][0][1] = 0.735554;
	default_alpha_prior[3][0][2] = 0.523268;
	default_alpha_prior[3][0][3] = 0.945254;
	default_alpha_prior[3][0][4] = 1.16698;
	default_alpha_prior[3][0][5] = 0.29715;
	default_alpha_prior[3][1][0] = 0.64427;
	default_alpha_prior[3][1][1] = 465.94;
	default_alpha_prior[3][1][2] = 0.486377;
	default_alpha_prior[3][1][3] = 1.00094;
	default_alpha_prior[3][1][4] = 0.695201;
	default_alpha_prior[3][1][5] = 0.466236;
	default_alpha_prior[3][2][0] = 0.27082;
	default_alpha_prior[3][2][1] = 0.29516;
	default_alpha_prior[3][2][2] = 199.88;
	default_alpha_prior[3][2][3] = 0.689022;
	default_alpha_prior[3][2][4] = 0.465843;
	default_alpha_prior[3][2][5] = 0.193624;
	default_alpha_prior[3][3][0] = 1.24395;
	default_alpha_prior[3][3][1] = 1.66739;
	default_alpha_prior[3][3][2] = 2.33255;
	default_alpha_prior[3][3][3] = 1432.05;
	default_alpha_prior[3][3][4] = 3.34616;
	default_alpha_prior[3][3][5] = 1.26537;
	default_alpha_prior[3][4][0] = 0.75205;
	default_alpha_prior[3][4][1] = 0.773375;
	default_alpha_prior[3][4][2] = 0.502719;
	default_alpha_prior[3][4][3] = 15.7907;
	default_alpha_prior[3][4][4] = 113.265;
	default_alpha_prior[3][4][5] = 0.248075;
	default_alpha_prior[3][5][0] = 0.405766;
	default_alpha_prior[3][5][1] = 0.42827;
	default_alpha_prior[3][5][2] = 0.187361;
	default_alpha_prior[3][5][3] = 15.1324;
	default_alpha_prior[3][5][4] = 0.185238;
	default_alpha_prior[3][5][5] = 11.3996;
	default_alpha_prior[3][6][0] = 54.045;
	default_alpha_prior[3][6][1] = 50.9537;
	default_alpha_prior[3][6][2] = 0.395894;
	default_alpha_prior[3][6][3] = 0.358199;
	default_alpha_prior[3][6][4] = 1.08529;
	default_alpha_prior[3][6][5] = 0.338531;
	default_alpha_prior[3][7][0] = 66.0728;
	default_alpha_prior[3][7][1] = 0.316336;
	default_alpha_prior[3][7][2] = 69.0806;
	default_alpha_prior[3][7][3] = 0.306993;
	default_alpha_prior[3][7][4] = 0.809972;
	default_alpha_prior[3][7][5] = 0.302452;
	default_alpha_prior[3][8][0] = 484.474;
	default_alpha_prior[3][8][1] = 1.92852;
	default_alpha_prior[3][8][2] = 1.93772;
	default_alpha_prior[3][8][3] = 505.099;
	default_alpha_prior[3][8][4] = 3.86887;
	default_alpha_prior[3][8][5] = 2.35758;
	default_alpha_prior[3][9][0] = 0.299784;
	default_alpha_prior[3][9][1] = 57.6077;
	default_alpha_prior[3][9][2] = 63.2437;
	default_alpha_prior[3][9][3] = 0.39431;
	default_alpha_prior[3][9][4] = 0.569814;
	default_alpha_prior[3][9][5] = 0.308767;
	default_alpha_prior[3][10][0] = 3.25275;
	default_alpha_prior[3][10][1] = 1118.47;
	default_alpha_prior[3][10][2] = 3.79864;
	default_alpha_prior[3][10][3] = 1149.4;
	default_alpha_prior[3][10][4] = 4.57076;
	default_alpha_prior[3][10][5] = 3.06625;
	default_alpha_prior[3][11][0] = 1.05466;
	default_alpha_prior[3][11][1] = 1.11556;
	default_alpha_prior[3][11][2] = 335.433;
	default_alpha_prior[3][11][3] = 346.991;
	default_alpha_prior[3][11][4] = 1.841;
	default_alpha_prior[3][11][5] = 1.14127;
	default_alpha_prior[3][12][0] = 0.501687;
	default_alpha_prior[3][12][1] = 0.506732;
	default_alpha_prior[3][12][2] = 0.367711;
	default_alpha_prior[3][12][3] = 52.0458;
	default_alpha_prior[3][12][4] = 38.2164;
	default_alpha_prior[3][12][5] = 0.169153;
	default_alpha_prior[3][13][0] = 0.202104;
	default_alpha_prior[3][13][1] = 0.330625;
	default_alpha_prior[3][13][2] = 0.175502;
	default_alpha_prior[3][13][3] = 22.2998;
	default_alpha_prior[3][13][4] = 0.113824;
	default_alpha_prior[3][13][5] = 7.8681;

	// printf("\n About the do some angle stuff \n\n");
	for (i = 0; i < NO_ALLELES; i++)
	{
		for (j = 0; j < MAX_GENOTYPES; j++)
		{
			double tot = 0.0;
			default_alpha_mag[i][j] = 0.0;
			for (k = 0; k < NO_ALLELES; k++)
				tot += default_alpha_prior[i][j][k];
			for (k = 0; k < NO_ALLELES; k++)
			{
				default_alpha_frac[i][j][k] = default_alpha_prior[i][j][k] / tot;
				default_alpha_mag[i][j] += default_alpha_frac[i][j][k] * default_alpha_frac[i][j][k];
			}
			default_alpha_mag[i][j] = sqrt(default_alpha_mag[i][j]);
		}
	}
	// printf("\n Finished setting up the first prior \n\n");
}
/*---------------------------------------------------------------------*/
#define pi2 6.283185307
double calc_angle(double *vec1, double mag1, double *vec2, double mag2, int n)
{
	int i;
	double tot = 0.0;
	for (i = 0; i < n; i++)
		tot += vec1[i] * vec2[i];
	double angle = tot / (mag1 * mag2);
	angle = acos(minim(0.999999999999, angle));
	return 360 * angle / pi2;
}
/*---------------------------------------------------------------------*/
void fill_alpha_prior(double **alpha, int max_gen, int dom_int)
{
	int i, j;

	for (i = 0; i < max_gen; i++)
		for (j = 0; j < NO_ALLELES; j++)
			alpha[i][j] = default_alpha_prior[dom_int][i][j];
}

/*---------------------------------------------------------------------*/
void fill_alpha_coef(double **alpha, double *coef, int max_gen)
{
	int i, j;
	double tot;

	for (j = 0; j < max_gen; j++)
	{
		coef[j] = 0;
		tot = 0.0;
		for (i = 0; i < NO_ALLELES; i++)
		{
			// alpha[j][i] = maxim(1,alpha[j][i]);
			tot += alpha[j][i];
			coef[j] -= gammln(alpha[j][i]);
		}
		coef[j] += gammln(tot);
	}
}

/*---------------------------------------------------------------------*/

double gammln(double xx)
{
	double x, tmp, ser;
	static double cof[6] = {76.18009173, -86.50532033, 24.01409822,
							-1.231739516, 0.120858003e-2, -0.536382e-5};
	int j;

	x = xx - 1.0;
	tmp = x + 5.5;
	tmp -= (x + 0.5) * log(tmp);
	ser = 1.0;
	for (j = 0; j <= 5; j++)
	{
		x += 1.0;
		ser += cof[j] / x;
	}
	return -tmp + log(2.50662827465 * ser);
}

/*---------------------------------------------------------------------*/
double exactfactln(int n)
{
	int i;
	double x = 1.0;

	for (i = 2; i <= n; i++)
		x *= (double)i;

	return log(x);
}
/*---------------------------------------------------------------------*/

double factln(int n)
{
	static double a[70001];

	if (n < 0)
	{
		printf("\n n = %d\n", n);
		dump_error("Negative factorial in routine FACTLN");
	}
	if (n <= 1)
		return 0.0;
	if (n <= 40)
		return a[n] ? a[n] : (a[n] = exactfactln(n));
	if (n <= 70000)
		return a[n] ? a[n] : (a[n] = gammln(n + 1.0));
	else
		return gammln(n + 1.0);
}

/* -------------------------------------------------------------- */
void dump_error(char *ss)
{
	printf("\n Unrecoverable error\n %s \n  Exiting now \n", ss);
	exit(1);
}
/* -------------------------------------------------------------- */

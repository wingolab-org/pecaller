#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h> 
#include <ctype.h>
#include <zlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define TRUE 1
#define FALSE 0
#define stdprn  ((FILE *)NULL)
#include <time.h>
 
#define minim(a,b) ((a<b)?a:b)
#define maxim(a,b) ((a>b)?a:b)
#define UNIQUE_MATE 0
#define UNIQUE_SLIP 1
#define UNIQUE_SINGLE 2 
#define UNIQUE_MIS 3
#define NON_MATE 4
#define NON_MIS 5
#define FRAG_MIS 6
#define NON_NO 7
#define NEITHER_MAP 8

#define MISALIGN_SLOP 10

typedef struct base_node 
{
  char ref;
  unsigned int pos;
  unsigned short As;
  unsigned short Cs;
  unsigned short Gs;
  unsigned short Ts;
  unsigned short Dels;
  unsigned short no_ins;
  char **ins;
} BASE_NODE;

typedef struct pthread_data_node 
{
  pthread_mutex_t mutex;
  long tid;
  char **read1;
  int *len1;
  char **read2;
  int *len2;
  long *read_no;
  int this_tot;
  int min_dist;
  int max_dist;
  int idepth;
  int *mapping_type;
  unsigned int *m1;
  unsigned int *m2;
} PTHREAD_DATA_NODE;

pthread_t *threads;
PTHREAD_DATA_NODE *pd_node_alloc(int min_dist,int max_dist,int idepth,int max);
int *ivector(int nl,int nh);
char *cvector(int nl,int nh);
unsigned char *ucvector(int nl,int nh);
unsigned long *ulvector(int nl,int nh);
long *lvector(int nl,int nh);
unsigned int *uvector(int nl,int nh);
unsigned long long *ullvector(int nl,int nh);
double *dvector(int nl, int nh);
double **dmatrix(int nrl,int nrh,int ncl,int nch);
char **cmatrix(int nrl,int nrh,int ncl,int nch);
unsigned char **ucmatrix(int nrl,int nrh,int ncl,int nch);
unsigned long **ulmatrix(int nrl,int nrh,int ncl,int nch);
unsigned int **umatrix(int nrl,int nrh,int ncl,int nch);
void free_cvector(char *v, int nl, int nh);
void free_ivector(int *v, int nl, int nh);
void free_ucvector(unsigned char *v, int nl, int nh);
void free_ulvector(unsigned long *v, int nl, int nh);
void free_uvector(unsigned int *v, int nl, int nh);
void free_dvector(double *v, int nl,int nh);
void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch);
void free_cmatrix(char **m,int nrl,int nrh,int ncl,int nch);
void free_ucmatrix(unsigned char **m,int nrl,int nrh,int ncl,int nch);
void free_ulmatrix(unsigned long **m,int nrl,int nrh,int ncl,int nch);
void free_umatrix(unsigned int **m,int nrl,int nrh,int ncl,int nch);
void dump_error(char *error_text);
int **imatrix(int nrl,int nrh,int ncl,int nch);
void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch);
void init_index_buffer(FILE *mfile,gzFile ifile);
static inline unsigned int *get_mers(unsigned int which,unsigned int *decode_length);
void init_penalty_matrices(double match_bonus,double ****penalty_matrix1,double ****penalty_matrix2,int ms);
void *map_everything(void *threadid);
int find_chrom(unsigned int *pos,int first, int last,int try,unsigned this);
int initial_map(char **read1,int seq_len,char *orient_return,unsigned int *match_returns,int idepth);
void reverse_transcribe(char *contig, char *s, int n);
static inline void convert_ct(char *s,int n);
void find_matches(unsigned int **mers,int max_depth,int *offsets, int idepth,
		  int *min_match,unsigned int *hits,int *hits_off,int *seg_matches,
		  int *tot_hits,char *orient,char or);
static inline void fill_mers(char *seq,unsigned int **forward_mers,int total_cuts,int *offsets, int **mismatch);
void fill_cv_mat(int *cv,unsigned ****mat);
static inline unsigned int convert_seq_int(char *seq,unsigned int ****mat,int *cv);
int sort_unsigned_int(const void *a,const void *b);
static inline double smith_waterman_align(BASE_NODE *base,int nn,char *seq,int mm, double ***S, double **match, double **gap_open, double **gap_extend,int *start);
int smith_waterman_backtrack(BASE_NODE *base,int nn,char *seq,int mm, double ***S, double **gap_open, double **gap_extend,int *start);
void reverse_compliment(char *contig, int n);
int find_mate_pairs(char **seq1, char **seq2,unsigned int *pos1, char *lor1,int rl1,int n1,unsigned int *pos2, char *lor2,int rl2,int n2,int *plen2,int *plen4,
		    int max_dist,int min_dist,char **read1,char **read2,int *len,unsigned int *start_match1,unsigned int *start_match2,
		    int *start1,int *start2,double ****p1_save,double ****p2_save,
		    double match_bonus,double **forward_mat, double **reverse_mat,double **gap_open,double **gap_extend,
		    BASE_NODE **base1,BASE_NODE **base2);
static inline void init_bonus_matrices(double match_bonus,double **forward_bonus_matrix,
			 double **reverse_bonus_matrix,double **gap_open,
			 double **gap_extend,int ms);
long init_gzfile_buffer(char *buf,long *buf_end,gzFile file);
static inline char *my_gzgets(gzFile file,char *buf,long *buf_end,long *buf_current);

static int MAX_THREADS;
static FILE *outfile;
static long genome_size;
static int pair_flag;
static long total_bases;
static long total_reads;
static long total_dist;
static long mate_counts[NEITHER_MAP+1];
static unsigned int *maps1, *maps2;
static long no_dists;

static double MIN_ALIGN = 0.9;
static int IS_BISULFITE = FALSE;
static unsigned int *contig_starts;
static int max_contigs,no_contigs;
#define MAX_READ_LENGTH 300 
static BASE_NODE *all_base_list;
static PTHREAD_DATA_NODE **all_thread_data;
static int reads_per_thread = 20000;

#define MAX_FILE_BUFFER 1200000000

static int max_hits = 200;
static int too_many_spots = 100;
static int **mismatch;
static unsigned int *pos_index;
static unsigned int *mers;
static pthread_mutex_t *all_base_mutex;
static unsigned int genome_chunk = 100;
static int genome_mutex_size;
static double match_bonus = 1.0;
static double **forward_bonus_matrix;
static double **reverse_bonus_matrix;
static double **gap_open;
static double **gap_extend;
static unsigned int ****seq_int_mat;
static int *seq_int_vector;
 
int main(int argc,char *argv[])
{
	char sss[4196],basename[1024];
	char sdxname[1024];
	char *sss1,*sss2;
	char **contig_names,*token;
	char *genome_buffer;
	long current_read = 0;
	int i,j,k;
	int file_num = 1;
	int min_dist,max_dist;

	long max_reads;
	long buffer1_current,buffer2_current;
	long buffer1_end,buffer2_end;
	char *file1_buf,*file2_buf;
	char **file1_names;
	char **file2_names;

	FILE *sfile;
	gzFile reffile,ifile;
	FILE *mfile1,*mfile2;
	FILE *mmfile;
	FILE *afile1,*afile2;
	FILE *summaryfile;
	gzFile indelfile,pileupfile;

	char mate_names[NEITHER_MAP+1][80];
	
  	gzFile innfile1,innfile2;
	sss1 = sss2 = NULL;
	file1_names = file2_names = NULL;
	afile1 = afile2 = NULL;

	char c_end = toupper(argv[3][0]);
	char c_array = toupper(argv[3][1]);
	min_dist = max_dist = 0;
	innfile1 = innfile2 = NULL;
        maps2 = NULL;
        mfile2 = NULL;

	file1_buf = cvector(0,MAX_FILE_BUFFER+100);
	file2_buf = NULL;
	buffer1_current = buffer2_current = 0;
	buffer1_end = buffer2_end = 0;

	max_reads = 30000;
	if(c_end != 'P' && c_end != 'S')
	{
	    printf("\nUsage: %s out_file sdx_file paired_or_single_or_array[p,s,pa,ps] file1 [file2] [max_dist] [min_dist] is_bisulfite[y,n] min_match_percentage max_threads max_reads\n",argv[0]);
	    exit(1);
	}
	if(c_end == 'S')
	{
	    pair_flag = FALSE;
	    if(argc != 9)
	    {
		printf("\nUsage: %s out_file sdx_file [s,sa] file1 is_bisulfite[y,n] min_match_percentage max_threads max_reads\n",argv[0]);
		exit(1);
	    }
	    MAX_THREADS = (int) atoi(argv[7]);
	    MIN_ALIGN = (double)atof(argv[6]);
	    max_reads = (long) atoi(argv[8]);
	    maps1 = malloc( sizeof(unsigned int)*max_reads);
	    if(c_array == 'A')
	    {
		file1_names = cmatrix(0,2000,0,512);
		if((afile1=fopen(argv[4],"r"))==(FILE *)NULL)
		{
		    printf("\n Can not open file %s for reading\n",argv[4]);
		    exit(1);
		}
		i = 0;
		fgets(file1_names[i],511,afile1);
		token = strtok(file1_names[i],"\t \n");
		while( (i < 2000) && (strlen(file1_names[i]) > 2) && (!feof(afile1)))
		{
		    i++;
		    file1_names[i][0] = '\0';
		    fgets(file1_names[i],511,afile1);
		    token = strtok(file1_names[i],"\t \n");

		}
		file_num = i;
		fclose(afile1);
	    }
	    else
	    {
		file1_names = cmatrix(0,1,0,512);
		strcpy(file1_names[0],argv[4]);
	    }
	    if(file_num >= 2000)
	      dump_error("\n Too many files found in array... 2000 is the limit \n");
	    if(!maps1)
	      dump_error("\n Could not allocate space for mapping position of reads1 \n");
	    if( (strchr(argv[5],'Y')) || (strchr(argv[5],'y')) )
	      IS_BISULFITE = TRUE;
	}
	else
	{
	    pair_flag = TRUE;
	    if(argc != 12) 
	    {
		printf("\nUsage: %s out_file sdx_file [p,pa] file1 file2 max_dist min_dist is_bisulfite[y,n] min_match_percentage max_threads max_reads\n",argv[0]);
		exit(1);
	    }
	    MAX_THREADS = (int) atoi(argv[10]);
	    MIN_ALIGN = (double)atof(argv[9]);
	    max_reads = (long) atol(argv[11]);
	    max_dist = (int) atoi(argv[6]);
	    min_dist = (int) atoi(argv[7]);
	    if( (strchr(argv[8],'Y')) || (strchr(argv[8],'y')) )
	      IS_BISULFITE = TRUE;

	    file2_buf = cvector(0,MAX_FILE_BUFFER+100);

	    maps1 = malloc( sizeof(unsigned int)*max_reads);
	    if(!maps1)
	      dump_error("\n Could not allocate space for mapping position of reads1 \n");
	    maps2 = malloc( sizeof(unsigned int)*max_reads);
	    if(!maps2)
	      dump_error("\n Could not allocate space for mapping position of reads2 \n");
	    if(c_array == 'A')
	    {
		file1_names = cmatrix(0,2000,0,512);
		file2_names = cmatrix(0,2000,0,512);
		if((afile1=fopen(argv[4],"r"))==(FILE *)NULL)
		{
		    printf("\n Can not open file %s for reading\n",argv[4]);
		    exit(1);
		}
		if((afile2=fopen(argv[5],"r"))==(FILE *)NULL)
		{
		    printf("\n Can not open file %s for reading\n",argv[5]);
		    exit(1);
		}
		i = 0;
		fgets(file1_names[i],511,afile1);
		token = strtok(file1_names[i],"\t \n");		    
		while( (i < 2000) && (strlen(file1_names[i]) > 2) && (!feof(afile1)))
		{
		    i++;
		    file1_names[i][0] = '\0';
		    fgets(file1_names[i],511,afile1);
		    token = strtok(file1_names[i],"\t \n");

		}
		file_num = i;
		fclose(afile1);
		i = 0;
		fgets(file2_names[i],511,afile2);
		token = strtok(file2_names[i],"\t \n");
		while( (i < 2000) && (strlen(file2_names[i]) > 2) && (!feof(afile2)))
		{
		    i++;
		    file2_names[i][0] = '\0';
		    fgets(file2_names[i],511,afile2);
		    token = strtok(file2_names[i],"\t \n");
		}
		if(file_num != i)
		  dump_error("\n Mismatch in number of files in the two arrays \n");
		fclose(afile2);

	    }
	    else
	    {
		file1_names = cmatrix(0,1,0,512);
		file2_names = cmatrix(0,1,0,512);
		strcpy(file1_names[0],argv[4]);
		strcpy(file2_names[0],argv[5]);
	    }
	    if(file_num >= 2000)
	      dump_error("\n Too many files found in array... 2000 is the limit \n");
	}

	if( (MAX_THREADS < 2) || (MAX_THREADS > 10000) )
	{
	    printf("\n Max_threads is not a sensible number (2,10000).  You gave %d \n",MAX_THREADS);
	    exit(1);
	}
	MAX_THREADS--;

	all_thread_data = malloc( sizeof(PTHREAD_DATA_NODE *) * MAX_THREADS);
	if(!all_thread_data)
	  dump_error("\n Could not allocate space for the PTHREAD_DATA_NODE vector \n");
	
	strcpy(basename,argv[1]);
	outfile = stdout;

	sprintf(sss,"%s.pileup.gz",basename);
	if((pileupfile=gzopen(sss,"wb"))==(gzFile)NULL)
	{
	    printf("\n Can not open file %s for writing\n",sss);
	    exit(1);
	}
	sprintf(sss,"%s.indel.txt.gz",basename);
	if((indelfile=gzopen(sss,"w"))==(gzFile)NULL)
	{
	    printf("\n Can not open file %s for writing\n",sss);
	    exit(1);
	}
	sprintf(sss,"%s.summary.txt",basename);
	if((summaryfile=fopen(sss,"w"))==(FILE *)NULL)
	{
	    printf("\n Can not open file %s for writing\n",sss);
	    exit(1);
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


	fgets(sss,256,sfile);
	max_contigs = atoi(sss);
	no_contigs = max_contigs;
	contig_starts = uvector(0,max_contigs);
	contig_names = cmatrix(0,max_contigs,0,256);
	seq_int_vector = ivector(0,256);
	if( (seq_int_mat = (unsigned int ****) malloc(sizeof(unsigned int ***)*4)) == NULL)
	  dump_error("\n Can not allocate seq_int mat \n");
	for(i=0;i<4;i++)
	{
	    if( (seq_int_mat[i] = (unsigned int ***) malloc(sizeof(unsigned int **)*4)) == NULL)
	      dump_error("\n Can not allocate seq_int mat i\n");
	    for(j=0;j<4;j++)
	    {
		if( (seq_int_mat[i][j] = (unsigned int **) malloc(sizeof(unsigned int *)*4)) == NULL)
		  dump_error("\n Can not allocate seq_int mat i,j\n");
		for(k=0;k<4;k++)		  
		  if( (seq_int_mat[i][j][k] = (unsigned int *) malloc(sizeof(unsigned int)*4)) == NULL)
		    dump_error("\n Can not allocate seq_int mat i,j,k\n");
	    }
	}
	fill_cv_mat(seq_int_vector,seq_int_mat);	

	contig_starts[0] = 0;
	for(i=0;i<max_contigs;i++)
	{
	    fgets(sss,1024,sfile);
	    token = strtok(sss,"\t \n");
	    contig_starts[i+1] = atoi(token);
	    token = strtok(NULL,"\t \n");
	    strcpy(contig_names[i],token);
	    // printf("\nFor contig %d we have offset %d",i,contig_starts[i]); 
	}
	fgets(sss,1024,sfile);
	int idepth = atoi(sss);
	fclose(sfile);
	for(i=1;i<=max_contigs;i++)
	  contig_starts[i] += contig_starts[i-1];

	// for(i=0;i<=max_contigs;i++)
	// printf("\n Contig %d starts at position %u \n",i,contig_starts[i]);

	genome_size = (long)contig_starts[max_contigs]+(long)15*max_contigs;
	genome_buffer = (char *)malloc(sizeof(char)*(genome_size+1));
	if(!genome_buffer)
	  dump_error("\n Failed to allocate memory for the Genome Buffer \n");
        else
	   printf("\n Genome size is %ld \n\n",genome_size);
	sprintf(sss,"%s.seq",sdxname);
	if((reffile=gzopen(sss,"r"))==(gzFile)NULL)
	{
	    printf("\n Can not open file %s for reading\n",sss);
	    exit(1);
	}
	sprintf(sss,"%s.mdx",sdxname);
	if((mmfile=fopen(sss,"r"))==(FILE *)NULL)
	{
	    printf("\nCould Not Open file %s\n",sss);
	    exit(1);
	}

	sprintf(sss,"%s.idx",sdxname);
	if((ifile=gzopen(sss,"r"))==(gzFile)NULL)
	{
	    printf("\nCould Not Open file %s\n",sss);
	    exit(1);
	}
	long g_temp = 0;
	printf("\n About to read genome \n\n");
	if(genome_size < MAX_FILE_BUFFER)	
		g_temp = gzread(reffile,(void *)genome_buffer,(unsigned int)sizeof(char)*genome_size);
	else
	{
		long count = 0;
		while(count < genome_size)
		{
			int ttemp = minim((long)genome_size - count,MAX_FILE_BUFFER);
			g_temp += gzread(reffile,(void *)&genome_buffer[count],ttemp);
			count += ttemp;
		}
	}
	gzclose(reffile);

	printf("\n About to allocate space for genome structure having read %ld\n\n",g_temp);
	long ts = genome_size*sizeof(BASE_NODE);
	ts++;
	all_base_list = (BASE_NODE *) malloc(ts);
	all_base_list += sizeof(BASE_NODE);
	if(!all_base_list) 
	  dump_error("\n Error allocating space for all_base_list");
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
		
	unsigned int pos;
	BASE_NODE *node;
	printf("\n About to initialize genome structure\n\n");
	for(pos=0;pos<genome_size;pos++)
	{
	  	node = &(all_base_list[pos]);
		node->ref = genome_buffer[pos]; 
		node->As = 0;
		node->Cs = 0;
		node->Gs = 0;
		node->Ts = 0;
		node->Dels = 0;
		node->no_ins = 0;
		node->pos = pos;
		node->ins = (char **)NULL;
	}

	free(genome_buffer);

	printf("\n About to allocate and initialize genome mutexes \n\n");
	genome_mutex_size = 1 + genome_size / genome_chunk;
	all_base_mutex = (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t)*genome_mutex_size);
	if(!all_base_mutex)
	  dump_error("\n Can not allocate space for all_base_mutex \n\n");

	for(i=0;i<genome_mutex_size;i++)
	  pthread_mutex_init(&(all_base_mutex[i]), NULL);

	printf("\n About to read kmers index \n\n");
	init_index_buffer(mmfile,ifile);	    
	fclose(mmfile);
	gzclose(ifile);
	printf("\n About to initialize all the thread mutexes \n\n");
	

	for(i=0;i<MAX_THREADS;i++)
	  all_thread_data[i] = pd_node_alloc(min_dist,max_dist,idepth,reads_per_thread);


	mismatch = imatrix(0,255,0,11);
	for(i=0;i<256;i++)
	{
		int a[4];
		int which = 0;
		unsigned int mask = 3;
		unsigned int shift = 0;
		for(j=0;j<4;j++,shift+=2)
		  a[j] = i & (mask << shift);
		//printf("\n%i a[0] = %d  a[1] = %d  a[2] = %d a[3] = %d",i,a[0],a[1],a[2],a[3]);
		shift = 0;
		for(j=0;j<4;j++,shift+=2)
		{
		    int start = i - a[j];
		    int this = (a[j] >> shift) & mask;
		    for(k=0;k<4;k++)
		      if(k!=this)		
			mismatch[i][which++] = start+(k << shift);
		}
	}

	if(pair_flag)
	{
	    sprintf(mate_names[UNIQUE_MATE],"Unique Mate-Paired");
	    sprintf(mate_names[UNIQUE_SLIP],"Unique Mate-Paired with slip");
	    sprintf(mate_names[UNIQUE_SINGLE],"Unique Single End");
	    sprintf(mate_names[UNIQUE_MIS],"Unique Mis-size");
	    sprintf(mate_names[NON_MATE],"Non-Unique Mate-Paired");
	    sprintf(mate_names[NON_MIS],"Non-Unique Mis-size");
	    sprintf(mate_names[FRAG_MIS],"Fragment Mismatch");
	    sprintf(mate_names[NON_NO],"Non-unique with no map");
	    sprintf(mate_names[NEITHER_MAP],"Neither Map");
	}
	else
	{
	    sprintf(mate_names[UNIQUE_MATE],"Not Used");
	    sprintf(mate_names[UNIQUE_SLIP],"Not Used");
	    sprintf(mate_names[UNIQUE_SINGLE],"Unique Mapping");
	    sprintf(mate_names[UNIQUE_MIS],"Not Used");
	    sprintf(mate_names[NON_MATE],"Not Used");
	    sprintf(mate_names[NON_MIS],"Not Used");
	    sprintf(mate_names[FRAG_MIS],"Not Used");
	    sprintf(mate_names[NON_NO],"Non-Unique Mapping, discarded");
	    sprintf(mate_names[NEITHER_MAP],"No mapping reaches threshold");
	}
	for(i=0;i<=NEITHER_MAP;i++)
		mate_counts[i] = 0;

	int ms = MAX_READ_LENGTH;
	gap_open = dmatrix(0,ms,0,ms);
	gap_extend = dmatrix(0,ms,0,ms);
	forward_bonus_matrix = dmatrix(0,ms,0,ms);
	reverse_bonus_matrix = dmatrix(0,ms,0,ms);

	init_bonus_matrices(match_bonus,forward_bonus_matrix,reverse_bonus_matrix,gap_open,gap_extend,MAX_READ_LENGTH);
	ms--;

	//pthread_t threads[MAX_THREADS];
	threads = (pthread_t *) malloc(sizeof(pthread_t)*MAX_THREADS);	


	    
	total_bases = total_reads = total_dist = no_dists = 0;
	printf("\n About to start mapping everything \n\n");
	int no_thread = -1;
	int this_read = 0;

	int iter;
	long tot_pairs = 0;

	for(iter = 0;iter<file_num;iter++)
	{
	    printf("\n About to open new set of files \n");
	    if((innfile1=gzopen(file1_names[iter],"r"))==(gzFile)NULL)
	    {
		printf("\n Can not open file %s for reading\n",file1_names[iter]);
		exit(1);
	    }
	    buffer1_current = init_gzfile_buffer(file1_buf,&buffer1_end,innfile1);
	    sprintf(sss,"%s.mfile",file1_names[iter]);
	    if((mfile1=fopen(sss,"wb"))==(FILE *)NULL)
	    {
		printf("\n Can not open file %s\n",sss);
		exit(1);
	    }
	    int seq_len2 = 0;
	    if(pair_flag)
	    {
		if((innfile2=gzopen(file2_names[iter],"r"))==(gzFile)NULL)
		{
		    printf("\n Can not open file %s for reading\n",file2_names[iter]);
		    exit(1);
		}
		buffer2_current = init_gzfile_buffer(file2_buf,&buffer2_end,innfile2);
		sprintf(sss,"%s.mfile",file2_names[iter]);
		if((mfile2=fopen(sss,"wb"))==(FILE *)NULL)
		{
		    printf("\n Can not open file %s\n",sss);
		    exit(1);
		}
		sss2 = my_gzgets(innfile2,file2_buf,&buffer2_end,&buffer2_current);
		sss2 = my_gzgets(innfile2,file2_buf,&buffer2_end,&buffer2_current);
		seq_len2 = strlen(sss2);
	    }
	    sss1 = my_gzgets(innfile1,file1_buf,&buffer1_end,&buffer1_current);
	    sss1 = my_gzgets(innfile1,file1_buf,&buffer1_end,&buffer1_current);

	    int seq_len = strlen(sss1);
	    int not_done = TRUE;
	    current_read = 0;
	    printf("\n Ready to map \n");
	    no_thread = 0;
	    pthread_mutex_lock(&(all_thread_data[no_thread]->mutex));
	    this_read = 0;
	    while(not_done && seq_len > 12)
	    {
		// if(current_read > 3600000)
		{
		all_thread_data[no_thread]->read_no[this_read] = current_read;
		memcpy(all_thread_data[no_thread]->read1[this_read],sss1,seq_len+1);
		all_thread_data[no_thread]->len1[this_read] = seq_len;
		if(pair_flag)
		{
		    memcpy(all_thread_data[no_thread]->read2[this_read],sss2,seq_len2+1);
		    all_thread_data[no_thread]->len2[this_read] = seq_len2;
		}
		this_read++;

		if(this_read >= reads_per_thread)
		{
		    long tid = no_thread;
		    int rc; 
		    all_thread_data[no_thread]->this_tot = this_read;
		    // printf("\t\tAbout to toss off a thread %ld\n",tid);
		    all_thread_data[no_thread]->tid = tid;
		    rc = pthread_create(&threads[tid],&attr,map_everything, (void *)all_thread_data[no_thread]);
		    if (rc)
		    {
			printf("ERROR; return code from pthread_create() is %d on read %ld\n", rc,current_read);
			exit(-1);
		    }
		    no_thread = -1;
		    while(no_thread < 0)
		      for(i=0;i<MAX_THREADS;i++)
			if(pthread_mutex_trylock(&(all_thread_data[i]->mutex)) == 0)
			{
			    // printf("I have found that thread %d is available\n",i);

			    all_thread_data[i]->this_tot = 0;
			    no_thread = i;
			    i = MAX_THREADS;
			}
		    this_read = 0;	
		}
	    	}

		if(current_read % 100000 == 0)
		  printf("\n We have read %ld reads and %ld have come back from successful mapping\n\n",current_read,total_reads);
		current_read++;
		if(current_read >=max_reads)
		  not_done = FALSE;
		else
		{
		    sss1 = my_gzgets(innfile1,file1_buf,&buffer1_end,&buffer1_current);
		    sss1 = my_gzgets(innfile1,file1_buf,&buffer1_end,&buffer1_current);
		    sss1 = my_gzgets(innfile1,file1_buf,&buffer1_end,&buffer1_current);
		    int not_there = TRUE;
		    seq_len = 0;
		    seq_len2 = 0;
		    while((sss1 != NULL) && not_there)
		    {
			if(sss1[0] == '@')
			  not_there = FALSE;
			sss1 = my_gzgets(innfile1,file1_buf,&buffer1_end,&buffer1_current);		    
		    }
		    if(not_there)
		      not_done = FALSE;
		    if(pair_flag && not_done)
		    {
			not_there = TRUE;
			sss2 = my_gzgets(innfile2,file2_buf,&buffer2_end,&buffer2_current);
			sss2 = my_gzgets(innfile2,file2_buf,&buffer2_end,&buffer2_current);
			sss2 = my_gzgets(innfile2,file2_buf,&buffer2_end,&buffer2_current);
			while((sss2 != NULL) && not_there)
			  {
			    if(sss2[0] == '@')
			      not_there = FALSE;
			    sss2 = my_gzgets(innfile2,file2_buf,&buffer2_end,&buffer2_current);
			  }
		    }
		    if(not_there)
		      not_done = FALSE;
		    else
		    {
			seq_len = strlen(sss1);
			if(pair_flag)
			  seq_len2 = strlen(sss2);
			  
		    }
		}
	    }
	    printf("\n Made it out alive, and have started cleanup \n\n");
	    if(this_read > 0)
	    {	
		long tid = no_thread;
		int rc; 
		// printf("\t\tAbout to toss off a thread %ld\n",tid);
		all_thread_data[no_thread]->tid = tid;
		all_thread_data[no_thread]->this_tot = this_read;
		rc = pthread_create(&threads[tid],&attr,map_everything, (void *)all_thread_data[no_thread]);
		if (rc)
		{
		    printf("ERROR; return code from pthread_create() is %d on read %ld\n", rc,current_read);
		    exit(-1);
		}	
		this_read = 0;
	    }

	    for(i=0;i<MAX_THREADS;i++)
	    {
		if(all_thread_data[i]->this_tot > 0)
			pthread_mutex_lock(&(all_thread_data[i]->mutex));
		all_thread_data[i]->this_tot = 0;
		pthread_mutex_unlock(&(all_thread_data[i]->mutex));
	    }
	    fwrite(maps1,sizeof(unsigned int),current_read,mfile1);
	    fclose(mfile1);
	    if(pair_flag)
	    {
		fwrite(maps2,sizeof(unsigned int),current_read,mfile2);
		fclose(mfile2);
	    }
	    gzclose(innfile1);
	    gzclose(innfile2);
	    tot_pairs += current_read;

	} // End For loop

	long no_bases = total_bases;

	if(no_bases <= 0)
	{
	    fprintf(summaryfile,"\n================================================================");
	    fprintf(summaryfile,"\n================= Summary ======================================");
	    fprintf(summaryfile,"\n================================================================");
	    fprintf(summaryfile,"\n================================================================");
	    fprintf(summaryfile,"\n\nTotal Number of Mapping reads of Any Kind\t0\tWith average Length\t0\tAverage Depth\t0\tAverage Insert Size\t0");
	    fprintf(summaryfile,"\n\nMapping Type\tCount\tFraction");
	    fprintf(summaryfile,"\nAll\t%ld\t1",tot_pairs);
	    for(i=0;i<=NEITHER_MAP;i++)
	      if(strstr(mate_names[i],"Not Used") == NULL)
		fprintf(summaryfile,"\n%s\t%ld\t%g",mate_names[i],mate_counts[i],(double)mate_counts[i]/(double)tot_pairs);
	    fprintf(summaryfile,"\n");
	    
	    fclose(summaryfile);
	    exit(1);
	}

	
	double avg_readlen = (double)total_bases;
	if(total_reads > 0)
	  avg_readlen /= (double) total_reads;
	
	double avg_dist = (double)total_dist;
	if(no_dists > 0)
	  avg_dist /= (double)no_dists;

	gzprintf(indelfile,"Fragment\tPositions\tReference Base\tTotal Coverage\tReference Reads\tNo Deletions\tNo Insertions\tInsertion Sequence");
	
	BASE_NODE *ref_seq;
	int ss = sizeof(unsigned short)*6;
	int s1 = sizeof(unsigned int);
	int tot_c;
	for(pos=0;pos<genome_size;pos++)
	  if( (tot_c = all_base_list[pos].As + all_base_list[pos].Cs + all_base_list[pos].Gs + all_base_list[pos].Ts + all_base_list[pos].Dels + all_base_list[pos].no_ins) > 0)
	  {
	      ref_seq = &(all_base_list[pos]);
	      unsigned short read_counts[6];
	      read_counts[0] = (unsigned short) ref_seq->As;
	      read_counts[1] = (unsigned short) ref_seq->Cs;
	      read_counts[2] = (unsigned short) ref_seq->Gs;
	      read_counts[3] = (unsigned short) ref_seq->Ts;
	      read_counts[4] = (unsigned short) ref_seq->Dels;
	      read_counts[5] = (unsigned short) ref_seq->no_ins;
	      gzwrite(pileupfile,(void *)&pos,s1);
	      gzwrite(pileupfile,(void *)read_counts,ss);
	    
	      if(ref_seq->no_ins > 0)
	      {
		  int ref_reads;
		  if(ref_seq->ref == 'A')
		    ref_reads = ref_seq->As;
		  else
		    if(ref_seq->ref == 'C')
		      ref_reads = ref_seq->Cs;
		    else
		      if(ref_seq->ref == 'G')
			ref_reads = ref_seq->Gs;
		      else
			ref_reads = ref_seq->Ts;
		  int which;
		  which = find_chrom(contig_starts,0,no_contigs-1,7,pos);
		  int contig_pos = 1 + pos - contig_starts[which];
		  gzprintf(indelfile,"\n%s\t%d\t%c\t%d\t%d\t%d\t%d",contig_names[which],
			   contig_pos,ref_seq->ref,tot_c,ref_reads,ref_seq->Dels,ref_seq->no_ins);
		  for(j=0;j<ref_seq->no_ins;j++)
		    gzprintf(indelfile,"\t%s",ref_seq->ins[j]);
	      }
	    
	  }
	gzclose(pileupfile);
	gzclose(indelfile);

	double avg_reads = (double)total_bases / (double) genome_size;

	printf("\n================================================================");
	printf("\n================= Summary ======================================");
	printf("\n================================================================");
	printf("\n================================================================");
	printf("\n\nTotal Number of Mapping reads of Any Kind\t%ld\tWith average Length\t%g\tAverage Depth\t%g\tAverage Insert Size\t%g",
	       total_reads,avg_readlen,avg_reads,avg_dist);
	printf("\n\nMapping Type\tCount\tFraction");
	printf("\nAll\t%ld\t1",tot_pairs);
	for(i=0;i<=NEITHER_MAP;i++)
	  if(strstr(mate_names[i],"Not Used") == NULL)
	    printf("\n%s\t%ld\t%g",mate_names[i],mate_counts[i],(double)mate_counts[i]/(double)tot_pairs);
	printf("\n");

	fprintf(summaryfile,"\n================================================================");
	fprintf(summaryfile,"\n================= Summary ======================================");
	fprintf(summaryfile,"\n================================================================");
	fprintf(summaryfile,"\n================================================================");
	fprintf(summaryfile,"\n\nTotal Number of Mapping reads of Any Kind\t%ld\tWith average Length\t%g\tAverage Depth\t%g\tAverage Insert Size\t%g",
			total_reads,avg_readlen,avg_reads,avg_dist);
	fprintf(summaryfile,"\n\nMapping Type\tCount\tFraction");
	fprintf(summaryfile,"\nAll\t%ld\t1",tot_pairs);
	for(i=0;i<=NEITHER_MAP;i++)
	  if(strstr(mate_names[i],"Not Used") == NULL)
		fprintf(summaryfile,"\n%s\t%ld\t%g",mate_names[i],mate_counts[i],(double)mate_counts[i]/(double)tot_pairs);

	fprintf(summaryfile,"\n");
	
	fclose(summaryfile);

	pthread_exit(NULL);

}
/*-------------------------------------------------------------------------------------------------------------------------------------- */
void *map_everything(void *threadid)
{
   long tid;
   PTHREAD_DATA_NODE *td;

   td = (PTHREAD_DATA_NODE *)threadid;
   tid = td->tid;
   int iter;
   int ms = MAX_READ_LENGTH;
   int i,j,k;
   long pool_size,pool_spot,psize,pline;

   // printf("\nIn map_everything with with thread = %ld\n\n",tid);   

   double *penalty_pool1,*penalty_pool2;
   double ****penalty_matrix1,****penalty_matrix2;
   penalty_pool1 = penalty_pool2 = NULL;
   penalty_matrix1 = penalty_matrix2 = NULL;
   char *orients1,*orients2;
   unsigned int *spots1,*spots2;

   spots1 = spots2 = NULL;
   orients1 = orients2 = NULL;

   pline = sizeof(double)*ms;
   psize = pline*ms;
   pool_size = max_hits * 3 * psize;
   penalty_pool1 = malloc(pool_size);
   if(!penalty_pool1)
     dump_error("\n Could not open space for penalty_pool1 \n\n"); 

   if( (penalty_matrix1 = (double ****) malloc(sizeof(double ***) * max_hits)) == NULL)
     dump_error("\n Error allocating space for the penalty matrix 1\n");
      
   pool_spot = 0;

   for(i=0;i<max_hits;i++)
   {
       if( (penalty_matrix1[i] = (double ***) malloc(sizeof(double **) * 3)) == NULL)
	 dump_error("\n Error allocating space for the penalty matrix 1,i\n");
       
       for(j=0;j<3;j++)
       {
	   if( (penalty_matrix1[i][j] = (double **) malloc(sizeof(double *) *ms)) == NULL)
	     dump_error("\n Error allocating space for penalty matrix 1, i,j \n");
	   for(k=0;k<ms;k++,pool_spot+=ms)
	     penalty_matrix1[i][j][k] = &penalty_pool1[pool_spot];
       }
   }

   // printf("\n Finished allocating penalty pool 1 \n\n");
   BASE_NODE **bases1,**bases2;
   int *pass_len1,*pass_len2;
   bases1 = bases2 = NULL;
   pass_len1 = pass_len2 = NULL;

   if( (bases1 = (BASE_NODE **) malloc(sizeof(BASE_NODE *) * (max_hits))) == NULL)
     dump_error("\n Error allocating space for BASES1 \n");
   pass_len1 = ivector(0,max_hits);

   if(pair_flag)
   {
       penalty_pool2 = malloc(pool_size);
       if(!penalty_pool2)
	 dump_error("\n Could not open space for penalty_pool2 \n\n"); 
       if( (penalty_matrix2 = (double ****) malloc(sizeof(double ***) * max_hits)) == NULL)
	 dump_error("\n Error allocating space for the penalty matrix 2\n");
       pool_spot = 0;
       for(i=0;i<max_hits;i++)
       {
	   if( (penalty_matrix2[i] = (double ***) malloc(sizeof(double **) * 3)) == NULL)
	     dump_error("\n Error allocating space for the penalty matrix 2\n");
	   for(j=0;j<3;j++)
	   {
	       if( (penalty_matrix2[i][j] = (double **) malloc(sizeof(double *) *ms)) == NULL)
		 dump_error("\n Error allocating space for penalty matrix 2, i,j \n");
	       for(k=0;k<ms;k++,pool_spot+=ms)
		 penalty_matrix2[i][j][k] = &penalty_pool2[pool_spot];
	   }

       }      

       if( (bases2 = (BASE_NODE **) malloc(sizeof(BASE_NODE *) * (max_hits))) == NULL)
	 dump_error("\n Error allocating space for BASES1 \n");
       pass_len2 = ivector(0,max_hits);
       orients2 = cvector(0,max_hits);
       spots2 = uvector(0,max_hits);
   }
   // printf("\n About to initialize penalty matrix \n\n");
   init_penalty_matrices(match_bonus,penalty_matrix1,penalty_matrix2,ms);
  
   // printf("\n Back from penalty matrix init wiht tid= %ld \n\n",tid);   
   char **iread1,**iread2;
   iread1 = cmatrix(0,1,0,MAX_READ_LENGTH);
   iread2 = cmatrix(0,1,0,MAX_READ_LENGTH);

   orients1 = cvector(0,max_hits);
   spots1 = uvector(0,max_hits);
   int idepth = td->idepth;
   int min_dist, max_dist;
   min_dist = td->min_dist;
   max_dist = td->max_dist;   
   for(iter=0;iter<td->this_tot;iter++)
   {
       char *read1,*read2;
       int hits1,hits2,len[5],m;
       int start1[3],start2[3],start_temp[3];
       // printf("\n About to start iteration %d\n\n",iter);       
       hits2 = 0;
       len[0] = len[1] = len[2] = len[3] = len[4] = 0;

       strcpy(iread1[0],td->read1[iter]);
       len[1] = td->len1[iter];
       reverse_transcribe(iread1[0],iread1[1],len[1]);

       // printf("\n In map_everything iter = %d with read1 = %s\n\n",iter,td->read1[iter]);
             
       hits1 = initial_map(iread1,len[1],orients1,spots1,idepth);
       if(pair_flag)
       {
	   // printf("\n In map_everything with read1 = %s and read2 = %s\n\n",td->read1[iter],td->read2[iter]);
	   strcpy(iread2[0],td->read2[iter]);
	   len[3] = td->len2[iter];
	   reverse_transcribe(iread2[0],iread2[1],len[3]);

	   len[3] = td->len2[iter];
	   // printf("\n about to do initial mapping with read 2 \n\n");
	   hits2 = initial_map(iread2,len[3],orients2,spots2,idepth);
       }
       // printf("\n Back from initial mapping with hits1 = %d  hits2 = %d for thread %ld\n\n",hits1,hits2,tid);     
       read1 = read2 = NULL; 
       unsigned int start_match1,start_match2;
       start_match1 = start_match2 = 0;
       unsigned int end_match;
       
       unsigned int bsm = 0;
       int first_call = 0;
 
  
       if(hits1 > 0)
       {
	   // printf("\n About to set up local base_list \n\n");
	   for(i=0;i<hits1;i++)
	   {
	       int this_chrom = find_chrom(contig_starts,0,no_contigs-1,7,spots1[i]);
	       unsigned int extra = 15*this_chrom;
	       long ttemp = maxim(0,(long)extra+(long)spots1[i]-(long)MISALIGN_SLOP);
	       start_match1 = maxim(contig_starts[this_chrom]+extra,ttemp);
	       end_match = minim(contig_starts[this_chrom+1]+extra,extra+spots1[i] + len[1]+MISALIGN_SLOP);
	       
	       int blen = 1 + end_match - start_match1;
	       pass_len1[i] = blen;
	       bases1[i] = &(all_base_list[start_match1]);
	       // printf("\n blen = %d this_chrom = %d start = %u end = %u bases1[i] = %ld spots1[i] = %u\n\n",blen,this_chrom,start_match1,end_match,(long)bases1[i],spots1[i]);

	   }
       }
       if(hits2 > 0)
       {
	   //	printf("\n In here with hits2 = %d \n\n",hits2);
	   for(i=0;i<hits2;i++)
	   {
	       int this_chrom = find_chrom(contig_starts,0,no_contigs-1,7,spots2[i]);
	       unsigned int extra = 15*this_chrom;
	       long ttemp = maxim(0,(long)extra+(long)spots2[i]-(long)MISALIGN_SLOP);
	       start_match2 = maxim(contig_starts[this_chrom]+extra,ttemp);
	       end_match = minim(contig_starts[this_chrom+1]+extra,extra+spots2[i] + len[3]+MISALIGN_SLOP);
	       
	       int blen = 1 + end_match - start_match2;
	       pass_len2[i] = blen;
	       bases2[i] = &(all_base_list[start_match2]); 
	       // printf("\n blen = %d this_chrom = %d start = %u end = %u bases2[i] = %ld spots2[i] = %u\n\n",blen,this_chrom,start_match2,end_match,(long)bases2[i],spots2[i]);
	   }
       }
       // printf("\n Made the base list \n\n");

       if(hits1 > 0 && hits2 == 0)    // First is present; second doesn't map
       {
	   double good_score = len[1] * MIN_ALIGN * match_bonus;
	   double top_score = -gap_open[0][0]*len[1];
	   int top_score_count = 0;
	   double this_score = 0;
	   read2 = NULL;
	   for(i=0;i<hits1;i++)
	   {
	       double **this_bonus = forward_bonus_matrix;
	       if(orients1[i])
		   this_bonus = reverse_bonus_matrix;
	       // printf("\n About to Smith_Waterman_align \n\n");
	       this_score = smith_waterman_align(bases1[i],pass_len1[i],iread1[(int)orients1[i]],len[1],penalty_matrix1[i],this_bonus,gap_open,gap_extend,start_temp);
	       // printf("\n Back from SM align with score %g \n\n",this_score);
	       if(this_score > top_score && this_score > good_score)
	       {
		   top_score = this_score;
		   top_score_count = 1;
		   bsm = i;
		   for(m=0;m<3;m++)
		     start1[m] = start_temp[m];
	       }
	       else
		 if(fabs(this_score - top_score) < 0.0001)
		   top_score_count++;
	   }
	   if(top_score_count == 0)
	   {
	       read1 = NULL;
	       first_call = NEITHER_MAP;
	   }
	   else
	     if(top_score_count == 1)
	     {
		 first_call = UNIQUE_SINGLE;
		 read1 = iread1[(int)orients1[bsm]];
		 start_match1 = bsm;
	     }
	     else
	     {
		 read1 = NULL;
		 first_call = NON_NO;
	     }
       
       }
       else
	 if(hits2 > 0 && hits1 == 0)    // Second is present; first doesn't map
	 {
	     double good_score = len[3] * MIN_ALIGN * match_bonus;
	     double top_score = -gap_open[0][0]*len[3];
	     int top_score_count = 0;
	     double this_score = 0;
	     read1 = NULL;
	 
	     for(i=0;i<hits2;i++)
	     {
		 double **this_bonus = forward_bonus_matrix;
		 if(orients2[i])
		   this_bonus = reverse_bonus_matrix;
		 
		 this_score = smith_waterman_align(bases2[i],pass_len2[i],iread2[(int)orients2[i]],len[3],penalty_matrix2[i],this_bonus,gap_open,gap_extend,start_temp);
		 if(this_score > top_score && this_score > good_score)
		 {
		     top_score = this_score;
		     top_score_count = 1;
		     bsm = i;
		     for(m=0;m<3;m++)
		       start2[m] = start_temp[m];
		 }
		 else
		   if(fabs(this_score - top_score) < 0.0001)
		     top_score_count++;
	     }
	     if(top_score_count == 0)
	     {
		 read2 = NULL;
		 first_call = NEITHER_MAP;
	     }
	     else
	       if(top_score_count == 1)
	       {
		   first_call = UNIQUE_SINGLE;
		   read2 = iread1[(int)orients2[bsm]];
		   start_match2 = bsm;
	       }
	       else
	       {
		   read2 = NULL;
		   first_call = NON_NO;
	       }
	 }
	 else
	   if(hits1 >0 && hits2 > 0)    // Both Present
	   {	
	       int this_match = find_mate_pairs(iread1,iread2,spots1,orients1,len[1],(int)hits1,spots2,
						orients2,len[3],(int)hits2,pass_len1,pass_len2,
						max_dist,min_dist,&read1,&read2,len,&start_match1,&start_match2,
						start1,start2,penalty_matrix1,penalty_matrix2,
						match_bonus,forward_bonus_matrix,reverse_bonus_matrix,gap_open,gap_extend,
						bases1,bases2);
	       
	       first_call = this_match;
	   }
	   else
	   {
	       if(hits1 <= 0 && hits2 <= 0)
		 first_call = NEITHER_MAP;
	       else
		 first_call = NON_NO;
	   }
       
       unsigned int temp_spot = 0;
       if(read1)
       {

	   // printf("About to backtrack with read1 %s for thread %ld with start_match1 = %d pass_len = %d  len = %d\n\n",
	  //	read1,tid,start_match1,pass_len1[start_match1],len[1]);

	   if(smith_waterman_backtrack(bases1[start_match1],pass_len1[start_match1],read1,len[1],penalty_matrix1[start_match1],gap_open,gap_extend,start1))
	   {
	       // printf("\n About to store Bn with start_match1 = %d \n\n",start_match1);
	       BASE_NODE *bn = bases1[start_match1];
	       temp_spot = bn[start1[1]].pos+1;
	   }		    
       }
       // printf("\n about to store temp_spot in iter = %d \n\n",iter);
       td->m1[iter] = temp_spot;
       // printf("\n Finished backtrack \n\n");
       temp_spot = 0;
       
       if(read2)
       {
	   // printf("\n In here and about to print something, but I think bad things might happen with start_match2 = %d pass_len2 = %d len3 = %d read2 = %ld \n\n",
	//	start_match2,pass_len2[start_match2],len[3],(long)read2);
	  // printf("About to backtrack with read2 %s for thread %ld with start_match2 = %d pass_len = %d  len = %d\n\n",
	  //	read2,tid,start_match2,pass_len2[start_match2],len[3]);
	   if(smith_waterman_backtrack(bases2[start_match2],pass_len2[start_match2],read2,len[3],penalty_matrix2[start_match2],gap_open,gap_extend,start2))
	   {
	       BASE_NODE *bn = bases2[start_match2];
	       temp_spot = bn[start2[1]].pos+1;
	   }
       }
       td->m2[iter] = temp_spot;  
       td->mapping_type[iter] = first_call;
       
       // printf("\n About to exit map_everything with tid = %ld\n\n",tid);
   } // ITER

 
   for(j=0;j<td->this_tot;j++)
   {
       mate_counts[td->mapping_type[j]]++;
       maps1[td->read_no[j]] = td->m1[j];
       if(td->m1[j])
       {
	   total_reads++;
	   total_bases += td->len1[j];
	   if(td->m2[j])
	   {
	       total_reads++;
	       total_bases += td->len2[j];
	       long test = labs(td->m1[j] - td->m2[j]);
	       if(test < max_dist*4)
	       {
		   total_dist += (long)test;
		   no_dists++;
	       }
	   }
       }
       else
	 if(td->m2[j])
	 {
	     total_reads++;
	     total_bases += td->len2[j];
	     maps2[td->read_no[j]] = td->m2[j];
	 }
   }
   td->this_tot = 0;
   free(bases1);
   free_ivector(pass_len1,0,max_hits);
   free_cvector(orients1,0,max_hits);
   free_uvector(spots1,0,max_hits);
          
   if(pair_flag)
   {
       free(bases2);
       free_ivector(pass_len2,0,max_hits);
       free_cvector(orients2,0,max_hits);
       free_uvector(spots2,0,max_hits);
   }
   

   free_cmatrix(iread1,0,1,0,MAX_READ_LENGTH);
   free_cmatrix(iread2,0,1,0,MAX_READ_LENGTH);
   for(i=0;i<max_hits;i++)
   {
       for(j=0;j<3;j++)
 	   free(penalty_matrix1[i][j]);
       free(penalty_matrix1[i]);
   }
   free(penalty_matrix1);
   free(penalty_pool1);

   if(pair_flag)
   {
       for(i=0;i<max_hits;i++)
       {
	   for(j=0;j<3;j++)
	     free(penalty_matrix2[i][j]);
	   free(penalty_matrix2[i]);
       }
       free(penalty_matrix2);
       free(penalty_pool2);
   }

 

   // printf("Exiting map everything with thread %ld\n\n",tid);
   pthread_mutex_unlock(&(td->mutex));
   pthread_exit(NULL);
}
/*-------------------------------------------------------------------------------------------------------------------------------------- */
/*-------------------------------------------------------------------------------------------------------------------------------------- */
int find_mate_pairs(char **seq1, char **seq2,unsigned int *pos1, char *lor1,int rl1,int n1,unsigned int *pos2, char *lor2,int rl2,int n2,int *plen2,int *plen4,
		    int max_dist,int min_dist,char **read1,char **read2,int *len,unsigned int *start_match1,unsigned int *start_match2,
		    int *start1,int *start2,double ****p1_save,double ****p2_save,
		    double match_bonus,double **forward_mat, double **reverse_mat,double **gap_open,double **gap_extend,
		    BASE_NODE **base1,BASE_NODE **base2)
{
	int k,i;
	int start_temp1[3];
	int start_temp2[3];
        int perfect = 0;
        int size_mismatch = 0;
        int broken = 0;
        int frag_count = 0;
	int l1,l3;
	double *smax1,*smax2;
	double this_1,this_2;
	double **this_mat;
	int best_len[5];
	long temp_dist = 0;
	int which1,which2;
	double tot_best = -1e5;
	read1[0] = NULL;
	read2[0] = NULL;
	
	start_temp1[0] = start_temp1[1] = start_temp1[2] = 0;
	start_temp2[0] = start_temp2[1] = start_temp2[2] = 0;

	if(n1 > 12000 || n2 > 12000)
	  return NON_MIS;

	smax1 = dvector(0,max_hits);
	smax2 = dvector(0,max_hits);
	for(i=0;i<=max_hits;i++)
	    smax1[i] = smax2[i] = -1.0;

	best_len[1] = best_len[2] = best_len[3] = best_len[4] = 0;

	best_len[1] = len[1] = l1 = rl1;
	double good_score1 = l1 * MIN_ALIGN * match_bonus;
	best_len[3] = len[3] = l3 = rl2;
	double good_score2 = l3 * MIN_ALIGN * match_bonus;

	// printf("\n\n Entered find_mate_pairs and looking for .%s. and .%s. n1 = %d  n2 = %d  l1 = %d  l3 = %d\n\n",seq1[0],seq2[0],n1,n2,l1,l3);
	for(which1=0;which1<n1;which1++)
	{
	    for(which2=0;which2<n2;which2++)
	    {
		unsigned int p1 = pos1[which1];
		unsigned int p2 = pos2[which2];
		temp_dist = labs((long)p1-(long)p2);
		// printf("\n p1 = %u  p2 = %u dist = %ld",p1,p2,temp_dist);
		if(temp_dist < 1e6)
		{
		    int or1 = lor1[which1];
		    int or2 = lor2[which2];
		    int is_perfect = ( (temp_dist >= min_dist) && (temp_dist <= max_dist) && (or1 != or2));
		    // printf("\n for i = %d   j = %d p1=%u  p2 = %u is_perfect = %d\n\n",which1,which2,p1,p2,is_perfect);
 		    if(is_perfect || (perfect < 1) )
		    {
			if(smax1[which1] < 0)
			{
			    this_mat = forward_mat;
			    if(or1)
			      this_mat = reverse_mat;
			    this_1 = smith_waterman_align(base1[which1],plen2[which1],seq1[or1],l1,p1_save[which1],this_mat,gap_open,gap_extend,start_temp1);
			    smax1[which1] = this_1;
			}
			else
			  this_1 = smax1[which1];
	
			if(smax2[which2] < 0)
			{
			    this_mat = forward_mat;
			    if(or2)
				this_mat = reverse_mat;

			    this_2 = smith_waterman_align(base2[which2],plen4[which2],seq2[or2],l3,p2_save[which2],this_mat,gap_open,gap_extend,start_temp2);
			    smax2[which2] = this_2;
			}
			else
			    this_2 = smax2[which2];

			double inc = this_1 + this_2 - tot_best;
			// printf("\n Checking whether or not to store %g and this_1 = %g  this_2 = %g\n\n",inc,this_1,this_2);
			if(inc > 0.001)
			{
			    // printf("\n About to store everything i=%d j = %d l1 = %d  l2 = %d  l3 = %d  l4 = %d this1 = %g  this2 = %g  good = %g",
			    //	 which1,which2,l1,l2,l3,l4,this_1,this_2,good_score1);
			    if(is_perfect && this_1 > good_score1 && this_2 > good_score2)
			    {
				perfect = 1;
				size_mismatch = 0;
			    }
			    else
			    {
			      	size_mismatch=1;
				// *dist = temp_dist;
			    }
			    frag_count = 0;
			    tot_best = this_1 + this_2;
			    best_len[2] = plen2[which1];
			    best_len[4] = plen4[which2];
			    *start_match1 = which1;
			    *start_match2 = which2;
			    good_score1 = this_1;
			    good_score2 = this_2;
			    for(k=0;k<3;k++)
			    {
				start1[k] = start_temp1[k];
				start2[k] = start_temp2[k];
			    }			    
			}
			else
			  if(inc > -0.001)
			  {
			      if(is_perfect)
				perfect++;
			      else
				size_mismatch++;
			  }

		    }
		    // printf("\n Made it to the bottom \n\n");
		}
		else
		  broken++;
	    }
	    
	}

	// printf("\n\tLeaving with smax1 = %g  smax2 = %g  perfect = %d l2=%d l4=%d \n\n",smax1[1],smax2[1],perfect,best_len[2],best_len[4]);

	int exit_code = NEITHER_MAP;
	if(frag_count > 0)
	{
	    if(perfect > 0)
	      exit_code = NON_MATE;
	    else
	      exit_code = NON_MIS;
	}
	else
	if(perfect < 1)
	{
	    int best1 = 0;
	    int best2 = 0;
	    for(i=1;i<n1;i++)
	      if(smax1[i] > smax1[best1])
		best1 = i;

	    for(i=1;i<n1;i++)
	      if(smax2[i] > smax2[best2])
		best2 = i;

	    if(smax1[best1] > good_score1)
	      if(smax2[best2] > good_score2)
	      {
		read1[0] = seq1[(int)lor1[best1]];
		read2[0] = seq2[(int)lor2[best2]];
		exit_code = UNIQUE_MIS;
	      }
	      else
	      {
		  read2[0] = NULL;
		  read1[0] = seq1[(int)lor1[best1]];
		  exit_code = UNIQUE_SINGLE;
	      }
	    else
	      if(smax2[best2] > good_score2)
	      {
		  read1[0] = NULL;
		  read2[0] = seq2[(int)lor2[best2]];
		  exit_code = UNIQUE_SINGLE;
	      }
	}
	else
	if( (perfect > 0) || (size_mismatch > 0) )
	{
	    len[2] = best_len[2];
	    len[4] = best_len[4];
	    read1[0] = seq1[(int)lor1[*start_match1]];
	    read2[0] = seq2[(int)lor2[*start_match2]];
		
	    if(perfect == 1)
	      exit_code = UNIQUE_MATE;
	    else
	      if(perfect > 1)
		exit_code = UNIQUE_SLIP;
	      else
		if(size_mismatch == 1)
		  exit_code = UNIQUE_MIS;
		else
		{
			read1[0] = NULL;
			read2[0] = NULL;
			exit_code = NON_MIS;
		}
	}
	else
	  if(broken > 0)
		exit_code = FRAG_MIS;

	free_dvector(smax1,0,max_hits);
	free_dvector(smax2,0,max_hits);
	return exit_code;
}

/*-------------------------------------------------------------------------------------------------------------------------------------- */
int initial_map(char **read1,int seq_len,char *orient_return,unsigned int *match_returns,int idepth)
{
  char *forward_seq,*reverse_seq;
  int total_cuts,*offsets,*hits_off;
  int max_mers = too_many_spots * 50;
  unsigned int **forward_mers, **reverse_mers,**temp_mers;
  unsigned int *hits;
  char *orient;
  int i,j;

  // printf("\n Just entered initial map with %s \n\n",read1);

  int N_limit = 1 + seq_len / 10;
  int n_count = 0;
  for(i=0;i<seq_len;i++)
    if(read1[0][i] == 'N')
      n_count++;

  if(n_count >= N_limit)
    return 0;
 
  forward_seq = cvector(0,seq_len);
  strcpy(forward_seq,read1[0]);
  reverse_seq = cvector(0,seq_len);
  strcpy(reverse_seq,read1[1]);

  if(IS_BISULFITE)
  {
      convert_ct(forward_seq,seq_len);
      convert_ct(reverse_seq,seq_len);
  }


  total_cuts = seq_len / idepth;
  if(seq_len % idepth == 0)
    total_cuts--;

  offsets = ivector(0,total_cuts);
  // printf("\n Made it here with %d total_cuts for seq %s\n\n",total_cuts,forward_seq);	    
  offsets[0] = 0;
  i = 1;
  while(i < total_cuts)
  {
      offsets[i] = offsets[i-1]+idepth;
      i++;
  }
  if(i == total_cuts)
    offsets[i] = seq_len - idepth;

  forward_mers = umatrix(0,total_cuts,0,max_mers);
  reverse_mers = umatrix(0,total_cuts,0,max_mers);
  temp_mers = umatrix(0,total_cuts,0,48);
  
  fill_mers(forward_seq,temp_mers,total_cuts,offsets,mismatch);
  for(i=0;i<=total_cuts;i++)
  {
      unsigned int this_tot = 0;
      unsigned int *tmer;
      forward_mers[i][0] = 0;
      for(j=0;j<=48;j++)
      {
	  tmer = get_mers(temp_mers[i][j],&this_tot);
	  if(this_tot >= too_many_spots)
	  {
	      forward_mers[i][0] = 0;
	      j = 48;
	  }
	  else
	  {
	      memcpy(&forward_mers[i][forward_mers[i][0]+1],tmer,this_tot*sizeof(unsigned int));
	      forward_mers[i][0] += this_tot;
	  }
      }
      if(forward_mers[i][0] > 1)
	qsort(&forward_mers[i][1],forward_mers[i][0],sizeof(unsigned int),sort_unsigned_int);
	
  }
  // printf("\n Made it here with %d total_cuts for reverse seq %s\n\n",total_cuts,reverse_seq);	    
  fill_mers(reverse_seq,temp_mers,total_cuts,offsets,mismatch);
  for(i=0;i<=total_cuts;i++)
  {
      unsigned int this_tot = 0;
      unsigned int *tmer;
      reverse_mers[i][0] = 0;
      for(j=0;j<=48;j++)
      {
	  tmer = get_mers(temp_mers[i][j],&this_tot);
	  if(this_tot >= too_many_spots)
	  {
	      reverse_mers[i][0] = 0;
	      j = 48;
	  }
	  else
	  {
	      memcpy(&reverse_mers[i][reverse_mers[i][0]+1],tmer,this_tot*sizeof(unsigned int));
	      reverse_mers[i][0] += this_tot;
	  }
      }
      if(reverse_mers[i][0] > 1)
	qsort(&reverse_mers[i][1],reverse_mers[i][0],sizeof(unsigned int),sort_unsigned_int);	
  }

  int min_match = maxim(1,total_cuts);
  if(total_cuts > 4)
    min_match = (4*(total_cuts))/5;
	       
  int *segs_matched;
  hits = uvector(0,max_hits);
  hits_off = ivector(0,max_hits);
  segs_matched = ivector(0,max_hits); 
  orient = cvector(0,max_hits);
  int tot1, tot2;
  tot1 = tot2 = 0;
  // printf("\n About to do the real work of initial mapping \n\n");

  find_matches(forward_mers,total_cuts,offsets,idepth,&min_match,hits,hits_off,segs_matched,&tot1,orient,0);
  tot2 = tot1;
  if(tot2 < max_hits)
    find_matches(reverse_mers,total_cuts,offsets,idepth,&min_match,hits,hits_off,segs_matched,&tot2,orient,1);
 


  for(i=0;i<tot2;i++)
  {
      orient_return[i] = orient[i];
      long temp = (long)hits[i] - (long)hits_off[i];
      match_returns[i] = maxim(0,temp);
  }

  free_cvector(forward_seq,0,seq_len);
  free_cvector(reverse_seq,0,seq_len);

  free_ivector(offsets,0,total_cuts);
  free_umatrix(forward_mers,0,total_cuts,0,max_mers);
  free_umatrix(reverse_mers,0,total_cuts,0,max_mers);
  free_umatrix(temp_mers,0,total_cuts,0,48);  

  free_uvector(hits,0,max_hits);
  free_ivector(hits_off,0,max_hits);
  free_ivector(segs_matched,0,max_hits);
  free_cvector(orient,0,max_hits);

  // printf("\n Found a total of %d hits  for %s \n\n",tot2,read1[0]);
  // printf("\n Exiting initial map \n\n");
  

  return tot2;
	    
}

/*-------------------------------------------------------------------------------------------------------------------------------------- */
/*-------------------------------------------------------------------------------------------------------------------------------------- */
inline double smith_waterman_align(BASE_NODE *base,int nn,char *seq,int mm, double ***S, double **match, double **gap_open, double **gap_extend,int *start)
{
  int maxi,maxj,maxk,i,j,i1,j1;
  double bump;

  maxk = 0;
  maxi = 0;
  maxj = mm;

 for(i=1;i<=nn;i++)
    for(j=1;j<mm;j++)
    {
	i1 = i-1;
	j1 = j-1;
	S[2][i][j] = maxim(S[0][i][j1]-gap_open[i1][j1],S[2][i][j1]-gap_extend[i1][j1]);
	S[1][i][j] = maxim(S[0][i1][j]-gap_open[i1][j],S[1][i1][j]-gap_extend[i1][j]);
	bump = match[(int)base[i1].ref][(int)seq[j1]];
	S[0][i][j] = maxim(maxim(S[0][i1][j1]+bump,S[1][i1][j1]+bump),S[2][i1][j1]+bump);
    }
  j = mm;
  j1 = j-1;
  for(i=1;i<=nn;i++)
  {
	i1 = i-1;
	S[2][i][j] = maxim(S[0][i][j1]-gap_open[i1][j1],S[2][i][j1]-gap_extend[i1][j1]);
	S[1][i][j] = maxim(S[0][i1][j]-gap_open[i1][j],S[1][i1][j]-gap_extend[i1][j]);
	bump = match[(int)base[i1].ref][(int)seq[j1]];
	S[0][i][j] = maxim(maxim(S[0][i1][j1]+bump,S[1][i1][j1]+bump),S[2][i1][j1]+bump);
	if(S[0][i][j] > S[maxk][maxi][maxj])
	{
	    maxk = 0; maxi = i; maxj = j; 
	}
	if(S[1][i][j] > S[maxk][maxi][maxj])
	{
	    maxk = 1; maxi = i; maxj = j; 
	}
	if(S[2][i][j] > S[maxk][maxi][maxj])
	{
	    maxk = 2; maxi = i; maxj = j; 
	}
   }
		
  start[0] = maxk;
  start[1] = maxi;
  start[2] = maxj;
  return S[maxk][maxi][maxj];
}

/*-------------------------------------------------------------------------------------------------------------------------------------- */
/*-------------------------------------------------------------------------------------------------------------------------------------- */
int smith_waterman_backtrack(BASE_NODE *base,int nn,char *seq,int mm, double ***S, double **gap_open, double **gap_extend,int *start)
{  
  double smax = 0.0;
  char ins_string[MAX_READ_LENGTH],**temp_ins;
  int ins_len;
  int maxi,maxj,maxk,i,j,k,i1,j1,m;
  int tg,last_lock = -1;

  for(i=0;i<MAX_READ_LENGTH;i++)
	ins_string[i] = '\0';

  k = start[0];
  i = start[1];
  j = start[2];
  
  i1 = j1 = 0;

  /* if(strstr(seq,"GACTGTGTTTATGCTTAGACATTAGGGCAGCCTGCCCCTGATGCTTCATTGAGCCTACCTGCTGCTTGCCTGCTATATTCGTCTTTGATCTTCTGGCTGCCTTCTT") != NULL)
  {
  printf("\n\t In backtrack with k=%d i = %d  j = %d with max = %g\n\n",k,i,j,S[k][i][j]);
  int kk,ii,jj;
   for(kk=0;kk<3;kk++)
  {
      printf("k=%d 0",kk);
      for(ii=0;ii<nn;ii++)
	printf(" %c",base[ii].ref);
      for(jj=0;jj<=mm;jj++)
      {
	  if(jj > 0)
	    printf("\n%c",seq[jj-1]);
	   else
	    printf("\n0");
	  for(ii=0;ii<=nn;ii++)
	    printf(" %g",S[kk][ii][jj]);
      }
      printf("\n");
  } 
  } */

  ins_len = 0;
  while(i > 0 && j > 0)
  {
      i1 = i-1;
      j1 = j-1;
      if(k == 0)
      {
	  maxi = i1; maxj = j1; maxk = 0;
	  smax = S[0][i1][j1];
	  if(S[1][maxi][maxj] > smax)
	  {
	      maxk = 1; smax = S[maxk][maxi][maxj];
	  }
	  if(S[2][maxi][maxj] > smax)
	    maxk = 2; 

      }
      else
      if(k == 2)
      {
	  maxk = 0; maxi = i; maxj = j1; smax = S[maxk][maxi][maxj] - gap_open[i1][maxj];
	  if(S[2][maxi][maxj] - gap_extend[i1][maxj] > smax)
	      maxk = 2; 
      }
      else 
      {
	  maxk = 0; maxi = i1; maxj = j; smax = S[maxk][maxi][maxj] - gap_open[maxi][maxj];
	  if(S[1][maxi][maxj] - gap_extend[maxi][maxj] > smax)
	      maxk = 1; 
      }

      tg = base[i1].pos / genome_chunk;
      if(tg != last_lock)
      {
	  if(last_lock >= 0)
	    pthread_mutex_unlock(&(all_base_mutex[last_lock]));
	  last_lock = tg;
	  pthread_mutex_lock(&(all_base_mutex[last_lock]));
      }
		

      if( maxi != i )
      {

	  if(maxj != j) 
	  {
	      if(seq[j1] == 'A')
		base[i1].As++;
	      else
		if(seq[j1] == 'T')
		  base[i1].Ts++;
		else
		  if(seq[j1] == 'G')
		    base[i1].Gs++;
		  else
		    if(seq[j1] == 'C')
		      base[i1].Cs++;
	  }
	  else
	  {
	    /* int nnn;
            char sss[1024];
            for(nnn=0;nnn<nn;nnn++)
                sss[nnn] = base[nnn]->ref;
            sss[nn] = '\0';
            fprintf(outfile,"\n Deletion at base %d in reference sequence %s comapred to %s  which is a %c Top score = %g\n",
	    i1,sss,seq,base[i1]->ref,S[start[0]][start[1]][start[2]]); */ 
	    base[i1].Dels++;
	  }

	  if(ins_len > 0)
	  {
	      // printf("\n About to do the deed with an existing number of insertions = %d and an insertion of len = %d \n\n",base[i1]->no_ins,ins_len);
	      if(base[i1].no_ins > 65000)
	      {
			printf("\n At base %u we found out %d insertions, and that is too many.\n Exiting now \n\n",base[i1].pos,base[i1].no_ins);
              }
	      if( (temp_ins = (char **) malloc(sizeof(char *) * (base[i1].no_ins+1))) == NULL)
		dump_error("\n Error allocating space for the insertion list \n");
	      if(base[i1].no_ins > 0)
	      {
		  for(m=0;m<base[i1].no_ins;m++)
		    temp_ins[m] = base[i1].ins[m];
		  
		 free(base[i1].ins);
	      }
	      base[i1].ins = temp_ins;
		  
	      base[i1].ins[base[i1].no_ins] = cvector(0,ins_len);
	      base[i1].ins[base[i1].no_ins][ins_len] = '\0';
	      for(m=0;m<ins_len;m++)
		base[i1].ins[base[i1].no_ins][m] = ins_string[ins_len - (m+1)];
	      /* printf("\n The deed is done \n\n");
              int nnn;
              char sss[1024];
              for(nnn=0;nnn<nn;nnn++)
                sss[nnn] = base[nnn]->ref;
              sss[nnn] = '\0';
              fprintf(outfile,"\n Insertion at base %d in reference sequence %s comapred to %s.  The insertion is %s Top score is %g\n\n",
	      i1,sss,seq,base[i1].ins[base[i1].no_ins],S[start[0]][start[1]][start[2]]);  */

	      base[i1].no_ins++;
	  }
	  ins_len = 0; 

      }
      else
      {
	  ins_string[ins_len] = seq[j1];
	  ins_len++;
      }

      i = maxi; j = maxj;  k = maxk;
  }
  if(ins_len > 0 && i >= 1)
  {
              tg = base[i1].pos / genome_chunk;
	      if(tg != last_lock)
	      {
			if(last_lock >= 0)
				pthread_mutex_unlock(&(all_base_mutex[last_lock]));
			last_lock = tg;
			pthread_mutex_lock(&(all_base_mutex[last_lock]));
	      }
	      // printf("\n\n Got Here \n\n");
              if( (temp_ins = (char **) malloc(sizeof(char *) * (base[i1].no_ins+1))) == NULL)
                dump_error("\n Error allocating space for the insertion list \n");
              if(base[i1].no_ins > 0)
              {
                  for(m=0;m<base[i1].no_ins;m++)
                    temp_ins[m] = base[i1].ins[m];

                 free(base[i1].ins);
              }
              base[i1].ins = temp_ins;

              base[i1].ins[base[i1].no_ins] = cvector(0,ins_len);
              base[i1].ins[base[i1].no_ins][ins_len] = '\0';
              for(m=0;m<ins_len;m++)
                base[i1].ins[base[i1].no_ins][m] = ins_string[ins_len - (m+1)];
              /* int nnn;
              char sss[1024];
              for(nnn=0;nnn<nn;nnn++)
                sss[nnn] = base[nnn]->ref;
              sss[nnn] = '\0';
              fprintf(outfile,"\n Insertion at base %d in reference sequence %s comapred to %s.  The insertion is %s \n",
                        i1,sss,seq,base[i1]->ins[base[i1]->no_ins]); */

              base[i1].no_ins++;
	      if(base[i1].no_ins > 65000)
	      {
			printf("\n At base %u we found out %d insertions, and that is too many.\n Exiting now \n\n",base[i1].pos,base[i1].no_ins);
              }
  }
  if(last_lock >= 0)
	pthread_mutex_unlock(&(all_base_mutex[last_lock]));
  // printf("\n Leaving BackTrack \n\n");

  return TRUE;
  
}

/*-------------------------------------------------------------------------------------------------------------------------------------- */

inline void fill_mers(char *seq,unsigned int **forward_mers,int total_cuts,int *offsets, int **mismatch)
{
  int i,j,k,m;
  unsigned int start;
  unsigned int a[4];
  unsigned int mask = 255;
  unsigned int shift[4];

  shift[0] = 0;
  shift[1] = 8;
  shift[2] = 16;
  shift[3] = 24;

  // printf("\n The sequence is %s and offset[0] is %d\n\n",seq,offsets[0]);
  for(i=0;i<=total_cuts;i++)
  {
      start = convert_seq_int(&seq[offsets[i]],seq_int_mat,seq_int_vector);
      forward_mers[i][0] = start;
      a[0] = start & mask;
      a[1] = start & (mask << 8);
      a[2] = start & (mask << 16);
      a[3] = start & (mask << 24);
      // printf("\n Start = %u  a[0] = %u  a[1] = %u a[2] = %u  a[3] = %u\n\n",start,a[0],a[1],a[2],a[3]);
      m = 1;
      for(j=0;j<4;j++)
      {
	  unsigned int this = start - a[j];
	  int which = a[j] >> shift[j];
	  // printf("\n In loop with j = %d  a[j] = %u  which = %d\n\n",j,a[j],which);
	  for(k=0;k<12;k++)
	      forward_mers[i][m++] = this+(mismatch[which][k] << shift[j]);
      }
  }
}

/*---------------------------------------------------------------------*/
inline void init_bonus_matrices(double match_bonus,double **forward_bonus_matrix,
			 double **reverse_bonus_matrix,double **gap_open,
			 double **gap_extend,int ms)
{
        int i,j;
	double mtemp = -1.0 / ((double)3.0 * match_bonus);

	for(i=0;i<=ms;i++)
	{
	    for(j=0;j<=ms;j++)
	      if(i==j)
		forward_bonus_matrix[i][j] = reverse_bonus_matrix[i][j] = match_bonus;
	      else
		forward_bonus_matrix[i][j] = reverse_bonus_matrix[i][j] = mtemp;
	    forward_bonus_matrix[i][(int)'N'] = forward_bonus_matrix[(int)'N'][i] = forward_bonus_matrix[i][(int)'n'] = forward_bonus_matrix[(int)'n'][i] = match_bonus;
	    reverse_bonus_matrix[i][(int)'N'] = reverse_bonus_matrix[(int)'N'][i] = reverse_bonus_matrix[i][(int)'n'] = reverse_bonus_matrix[(int)'n'][i] = match_bonus;
	    if(IS_BISULFITE)
	    {
		forward_bonus_matrix[(int)'C'][(int)'T'] = match_bonus;
		forward_bonus_matrix[(int)'C'][(int)'t'] = match_bonus;
		forward_bonus_matrix[(int)'c'][(int)'T'] = match_bonus;
		forward_bonus_matrix[(int)'c'][(int)'t'] = match_bonus;
		reverse_bonus_matrix[(int)'C'][(int)'T'] = match_bonus;
		reverse_bonus_matrix[(int)'C'][(int)'t'] = match_bonus;
		reverse_bonus_matrix[(int)'c'][(int)'T'] = match_bonus;
		reverse_bonus_matrix[(int)'c'][(int)'t'] = match_bonus;
	    }
	}
	// printf("\n About to read the sequence \n\n");
	

	double gopen_penalty = 2.0*match_bonus;
	double ge = match_bonus / 36.0;
	for(i=0;i<=ms;i++)
	  for(j=0;j<=ms;j++)
	  {
	      gap_open[i][j] = gopen_penalty;
	      gap_extend[i][j] = ge;
	  }

}
/*---------------------------------------------------------------------*/
void init_penalty_matrices(double match_bonus,double ****penalty_matrix1,double ****penalty_matrix2,int ms)
{
	int i;
	int jj;
	double gopen_penalty = 2.0*match_bonus;
	double ge = match_bonus / 36.0;

	for(jj=0;jj<max_hits;jj++)
	{
	    
	    penalty_matrix1[jj][0][0][0] = 0.0;
	    penalty_matrix1[jj][1][0][0] = 0.0; 
	    penalty_matrix1[jj][2][0][0] = -1.0*gopen_penalty;
	    if(pair_flag)
	    {
		penalty_matrix2[jj][0][0][0] = 0.0;
		penalty_matrix2[jj][1][0][0] = 0.0; 
		penalty_matrix2[jj][2][0][0] = -1.0*gopen_penalty;
	    }


	    for(i=1;i<ms;i++)
	    {
		/* penalty_matrix[0][0][i] = penalty_matrix[1][0][i] = penalty_matrix[0][i][0] = penalty_matrix[2][i][0] 
		   = -(gap_open + (double)(i-1)*gap_extend);  */	
		penalty_matrix1[jj][0][0][i] = penalty_matrix1[jj][1][0][i] = penalty_matrix1[jj][2][0][i] = -(gopen_penalty + (double)(i-1)*ge);
		penalty_matrix1[jj][0][i][0] = penalty_matrix1[jj][0][0][0];
		penalty_matrix1[jj][1][i][0] = penalty_matrix1[jj][1][0][0];
		penalty_matrix1[jj][2][i][0] = penalty_matrix1[jj][2][0][0]; 

		if(pair_flag)
		{
		    penalty_matrix2[jj][0][0][i] = penalty_matrix2[jj][1][0][i] = penalty_matrix2[jj][2][0][i] = -(gopen_penalty + (double)(i-1)*ge);
		    penalty_matrix2[jj][0][i][0] = penalty_matrix2[jj][0][0][0];
		    penalty_matrix2[jj][1][i][0] = penalty_matrix2[jj][1][0][0];
		    penalty_matrix2[jj][2][i][0] = penalty_matrix2[jj][2][0][0]; 
		}
	    }
	    

	}
}
/*---------------------------------------------------------------------*/

void reverse_compliment(char *contig, int n)
{
  int i,j=0;
  char c,s[4128];

  for(i=n-1;i>=0;i--)
  {
      if(*(contig+i) == 'A')
	c = 'T';
      else
      if(*(contig+i) == 'C')
	c = 'G';
      else
      if(*(contig+i) == 'G')
	c = 'C';
      else
      if(*(contig+i) == 'T')
	c = 'A';
      else
	c = 'N';

      s[j++] = c;
  }
  for(i=0;i<n;i++)
    contig[i] = s[i];
  contig[i] = '\0';

}
/*---------------------------------------------------------------------*/

/*---------------------------------------------------------------------*/
void init_index_buffer(FILE *mfile,gzFile ifile)
{
  long decode_length;
  long buffer_chunk = 4294967296;
  int i;
  long bpie = buffer_chunk / 16;
  pos_index = (unsigned int *)calloc((long)(buffer_chunk+1),sizeof(unsigned int));
  if(!pos_index)
    dump_error("\n Can not allocate space for the position index \n");
  unsigned int *fake_index;
  fake_index = pos_index;
  unsigned l = (unsigned int) bpie * sizeof(unsigned int);
  for(i=0;i<15;i++)
  {
	// printf("\n Reading bit i = %i \n\n",i);
  	gzread(ifile,(void *)fake_index,l);
	fake_index += bpie;
  }
  gzread(ifile,(void *)fake_index,l+sizeof(unsigned int));
  decode_length = pos_index[buffer_chunk];
  // printf("\n Just read the last bit with the file number being %u %ld\n\n",pos_index[buffer_chunk-1],decode_length);
  mers = (unsigned int *)malloc((long)(decode_length+1)*sizeof(unsigned int));
  if(!mers)
    dump_error("\n Can not allocate space for the mers index\n"); 
  fread(mers,sizeof(unsigned int),pos_index[buffer_chunk],mfile);
}
				     
/*---------------------------------------------------------------------*/
inline unsigned int *get_mers(unsigned int which,unsigned int *decode_length)
{
  //if(start_mer > 4541000)
  //  printf("\n Here with start_mer = %u which = %u current_start = %ld",start_mer,which,current_start);
  *decode_length = pos_index[(long)(which+1)] - pos_index[which];
  return &mers[pos_index[which]];
}			     

/*-------------------------------------------------------------------------------------------------------------------------------------- */
int find_chrom(unsigned int *pos,int first, int last,int try,unsigned this)
{
  // printf("\n first = %d last = %d try  = %d this = %u pos[try] = %u",first,last,try,this,pos[try]);
  if(first == last)
    return first;

  if(pos[try] <= this && pos[try+1] >= this)
    return try;

  if(pos[try] > this)
    last = try-1;
  else
    first = try + 1;

  try = (last + first) / 2;

  return find_chrom(pos,first,last,try,this);
}
/*-------------------------------------------------------------------------------------------------------------------------------------- */
void find_matches(unsigned int **mers,int max_depth,int *offsets, int idepth,
		  int *min_match,unsigned int *hits,int *hits_off,int *seg_matches,
		  int *tot_hits,char *orient,char or)
{

	unsigned int tot_found;
	unsigned int mer_pos[10000];
	unsigned int i,j,k,loop,max_off=maxim(2,idepth-4);
	int start, end;
	unsigned int min_spots = 10000;

	for(i=0;i<=max_depth;i++)
	  min_spots = minim(min_spots,mers[i][0]);
	
	if(min_spots > max_hits)
	{
	    *tot_hits = 0;
	    return;
	}
			    
	/* for(i=0;i<=max_depth;i++)
	{
	    printf("\n Piece %u has %u options",i,mers[i][0]);
	    for(j=1;j<=mers[i][0];j++)
	      printf(" %u",mers[i][j]);
	      } */

	for(loop=0;loop <= 1 + max_depth-(*min_match);loop++)
	{
		start = -(offsets[loop] + max_off);
		end = max_off;
		for(j=loop+1;j<=max_depth;j++)
			end = maxim(end,max_off+offsets[j]-offsets[loop]);
		
		for(i=loop;i<=max_depth;i++)
		  mer_pos[i] = 1;

		// printf("\n loop =%d offset = %d Start = %d   end = %d max_depth = %d",loop,offsets[loop],start,end,max_depth);

		for(i=1;i<=mers[loop][0];i++)
		{
		    long this_start = (long)mers[loop][i] + (long)start;
		    long this_end = (long)mers[loop][i] + (long)end;
		    this_start = maxim(this_start,0);
		    this_end = maxim(this_end,0);
		    for(j=loop+1;j<=max_depth;j++)
			while( (mer_pos[j] < mers[j][0]) && (mers[j][mer_pos[j]] < this_start))
			  mer_pos[j]++;

		    // for(j=loop+1;j<=max_depth;j++)
		    //  printf("\n Starting element for list %d is %d out of %d",j,mer_pos[j],mers[j][0]);

		    tot_found = 1;
		    for(j=loop+1;j<=max_depth;j++)
		      for(k=mer_pos[j];(k<=mers[j][0] && mers[j][k]<= this_end) ; k++)
			if( abs( (mers[loop][i]-mers[j][k]) - (offsets[loop]-offsets[j])) < max_off)
			{
			    tot_found++;
			    // printf("\n Found a part at j = %d k = %d",j,k);
			    k = mers[j][0]+1;
			}

		    if(tot_found > *min_match)
		    {
			*min_match = tot_found;
			*tot_hits = 0;
			seg_matches[*tot_hits] = tot_found;
			hits[*tot_hits] = mers[loop][i];
			hits_off[*tot_hits] = offsets[loop];
			orient[*tot_hits] = or;
			(*tot_hits)++;
			// printf("\n Nailed it tot_found = %d loop = %d i = %d pos = %d offset = %d total_hits = %d",tot_found,loop,i,mers[loop][i],offsets[loop],*tot_hits);
		    }
		    else
		      if(tot_found == *min_match )
		      {
			if(*tot_hits < max_hits)
			{
			    int new = 1;
			    for(k=0;k<*tot_hits;k++)
			      if(hits[k] - hits_off[k] == mers[loop][i] - offsets[loop])
			      {
				  k = *tot_hits;
				  new = 0;
			      }
			    if(new)
			    {
				seg_matches[*tot_hits] = tot_found;
				hits[*tot_hits] = mers[loop][i];
				hits_off[*tot_hits] = offsets[loop];
				orient[*tot_hits] = or;
				(*tot_hits)++;
				// printf("\n Nailed it tot_found = %d loop = %d i = %d pos = %d offset = %d total_hits = %d",tot_found,loop,i,mers[loop][i],offsets[loop],*tot_hits);
			    }
			}
			else
				return;
		      }
		}
	}

}
/*-------------------------------------------------------------------------------------------------------------------------------------- */
inline void convert_ct(char *s,int n)
{
  int i;
  for(i=0;i<n;i++)
    if(s[i] == 'C')
      s[i] = 'T';
  
}
/*-------------------------------------------------------------------------------------------------------------------------------------- */
void reverse_transcribe(char *contig, char *s, int n)
{
  int i;
  char c;
  
  for(i=n-1;i>-1;i--)
  {
      if(*(contig+i) == 'A')
	c = 'T';
      else
      if(*(contig+i) == 'C')
	c = 'G';
      else
      if(*(contig+i) == 'G')
	c = 'C';
      else
      if(*(contig+i) == 'T')
	c = 'A';
      else
      if(*(contig+i) == 'W')
	c = 'W';
      else
      if(*(contig+i) == 'S')
	c = 'S';
      else
      if(*(contig+i) == 'K')
	c = 'M';
      else
      if(*(contig+i) == 'M')
	c = 'K';
      else
      if(*(contig+i) == 'Y')
	c = 'R';
      else
      if(*(contig+i) == 'R')
	c = 'Y';
      else
	c = 'N';
      *s++ = c;
  }
  *s = '\0';

}

/*---------------------------------------------------------------------*/
PTHREAD_DATA_NODE *pd_node_alloc(int min_dist,int max_dist,int idepth,int max)
{
  PTHREAD_DATA_NODE *pn;
  pn = (PTHREAD_DATA_NODE *)malloc(sizeof(PTHREAD_DATA_NODE));
  if(!pn)
    dump_error("\n Can not allocate space for a pthread_data_node \n");

  pn->this_tot = 0;

  pn->read1 = cmatrix(0,max,0,MAX_READ_LENGTH);
  pn->read2 = cmatrix(0,max,0,MAX_READ_LENGTH);
  pn->len1 = ivector(0,max);
  pn->len2 = ivector(0,max);
  pn->read_no = lvector(0,max);
  pn->m1 = uvector(0,max);
  pn->m2 = uvector(0,max);
  pn->mapping_type = ivector(0,max);
  int i;
  for(i=0;i<=max;i++)
  {
	pn->read1[i][0] = pn->read2[i][0] = '\0';
        pn->len1[i] = pn->len2[i] = 0;
	pn->read_no[i] = 0;
	pn->m1[i] = pn->m2[i] = 0;
  }
  pn->idepth = idepth;
  pn->min_dist = min_dist;
  pn->max_dist = max_dist;
  pthread_mutex_init(&(pn->mutex), NULL);
  return pn;
		     
}
/*---------------------------------------------------------------------*/
void fill_cv_mat(int *cv,unsigned ****mat)
{
  unsigned int i,j,k,m;
  for(i=0;i<256;i++)
    cv[i] = 0;
  cv['c'] = cv['C'] = 1;
  cv['g'] = cv['G'] = 2;
  cv['t'] = cv['T'] = 3;

  unsigned int up[3];

  for(i=0;i<4;i++)
  {
      up[0] = i;
      up[0] = up[0] << 2;
      for(j=0;j<4;j++)
      {
	  up[1] = up[0] + j;
	  up[1] = up[1] << 2;
	  for(k=0;k<4;k++)
	  {
	      up[2] = up[1] + k;
	      up[2] = up[2] << 2;
	      for(m=0;m<4;m++)
		  mat[i][j][k][m] = up[2] + m;
	  }
      }
  }

}
/*---------------------------------------------------------------------*/
inline unsigned int convert_seq_int(char *seq,unsigned int ****mat,int *cv)
{
  unsigned int this_mer;

  this_mer = mat[cv[(int)seq[0]]][cv[(int)seq[1]]][cv[(int)seq[2]]][cv[(int)seq[3]]];
  this_mer = this_mer << 8;
  this_mer += mat[cv[(int)seq[4]]][cv[(int)seq[5]]][cv[(int)seq[6]]][cv[(int)seq[7]]];
  this_mer = this_mer << 8;
  this_mer += mat[cv[(int)seq[8]]][cv[(int)seq[9]]][cv[(int)seq[10]]][cv[(int)seq[11]]];
  this_mer = this_mer << 8;
  this_mer += mat[cv[(int)seq[12]]][cv[(int)seq[13]]][cv[(int)seq[14]]][cv[(int)seq[15]]];
  
  return this_mer;

}

/*---------------------------------------------------------------------*/
int sort_unsigned_int(const void *a,const void *b)
{
        if(*((unsigned int *)a)<*((unsigned int *)b))
                return -1;
        else
        if(*((unsigned int *)a)>*((unsigned int *)b))
                return 1;
        else
                return 0;
}

/*-------------------------------------------------------------------------------------------------------------------------------------- */
long init_gzfile_buffer(char *buf,long *buf_end,gzFile file)
{
  (*buf_end) = gzread(file,buf,MAX_FILE_BUFFER);
  // printf("\n Leaving init_gzbuffer with and end of %ld \n",*buf_end);
  return 0;
}

/*-------------------------------------------------------------------------------------------------------------------------------------- */
static inline char *my_gzgets(gzFile file,char *buf,long *buf_end,long *buf_current)
{
  long start = *buf_current;
  // printf("\n Buf_current = %ld  end = %ld c = %c",*buf_current,*buf_end,buf[*buf_current]);
  while((*buf_current) < *buf_end )
  {
      if(buf[*buf_current] == '\n')
      {
	  buf[*buf_current] = '\0';
	  (*buf_current)++;
	  long temp = *buf_current - start;
	  if(temp >= MAX_READ_LENGTH)
	    buf[temp-1] = '\0';
	  return &buf[start];
      }
      (*buf_current)++;
  }
 
  if(gzeof(file))
    return NULL;
  long i = 0;
  long j;
  for(j=start;j<*buf_end;j++,i++)
    buf[i] = buf[j];

  (*buf_end) = gzread(file,&buf[i],MAX_FILE_BUFFER-i);
  *buf_current = 0;
  if( (*buf_end) > 50)
  {
      (*buf_end) += i;
      return my_gzgets(file,buf,buf_end,buf_current);
  }
  else
    return NULL;

}

/*---------------------------------------------------------------------*/

char *cvector(int nl,int nh)
{
	char *v;
 
	v=(char *)malloc((unsigned) (nh-nl+1)*sizeof(char));
	if (!v) dump_error("allocation failure in cvector()");
	return v-nl;
}

unsigned char *ucvector(int nl,int nh)
{
	unsigned char *v;
 
	v=(unsigned char *)malloc((unsigned) (nh-nl+1)*sizeof(unsigned char));
	if (!v) dump_error("allocation failure in ucvector()");
	return v-nl;
}

unsigned long long *ullvector(int nl,int nh)
{
	unsigned long long *v;
 
	v=(unsigned long long *)malloc((unsigned) (nh-nl+1)*sizeof(unsigned long long));
	if (!v) dump_error("allocation failure in ullvector()");
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
	if (!v) dump_error("allocation failure in ulvector()");
	return v-nl;
}
long *lvector(int nl,int nh)
{
	long *v;
 
	v=(long *)malloc((nh-nl+1)*sizeof(long));
	if (!v) dump_error("allocation failure in lvector()");
	return v-nl;
}

unsigned long *ulvector(int nl,int nh)
{
	unsigned long *v;
 
	v=(unsigned long *)malloc((unsigned) (nh-nl+1)*sizeof(long));
	if (!v) dump_error("allocation failure in ulvector()");
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

unsigned int **umatrix(int nrl,int nrh,int ncl,int nch)
{
        int i;
	unsigned int **m;

	m=(unsigned int **)malloc((unsigned) (nrh-nrl+1)*sizeof(unsigned int*));
	if (!m) dump_error("allocation failure 1 in ulmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(unsigned int *)malloc((unsigned) (nch-ncl+1)*sizeof(unsigned int));
		if (!m[i]) dump_error("allocation failure 2 in ulmatrix()");
		m[i] -= ncl;
	}
	return m;
}

unsigned long **ulmatrix(int nrl,int nrh,int ncl,int nch)
{
        int i;
	unsigned long **m;

	m=(unsigned long **)malloc((unsigned) (nrh-nrl+1)*sizeof(unsigned long*));
	if (!m) dump_error("allocation failure 1 in ulmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(unsigned long *)malloc((unsigned) (nch-ncl+1)*sizeof(unsigned long));
		if (!m[i]) dump_error("allocation failure 2 in ulmatrix()");
		m[i] -= ncl;
	}
	return m;
}

void free_umatrix(unsigned int **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((void *) (m[i]+ncl));
	free((void *) (m+nrl));
}


void free_ulmatrix(unsigned long **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((void *) (m[i]+ncl));
	free((void *) (m+nrl));
}


void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((void *) (m[i]+ncl));
	free((void *) (m+nrl));
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
 
void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch)
{
	int i;
 
	for(i=nrh;i>=nrl;i--) free((void *) (m[i]+ncl));
	free((void *) (m+nrl));
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


unsigned char **ucmatrix(int nrl,int nrh,int ncl,int nch)
{
	int i;
	unsigned char **m;
 
	m=(unsigned char **) malloc((unsigned) (nrh-nrl+1)*sizeof(unsigned char*));
	if (!m) dump_error("allocation failure 1 in cmatrix()");
	m -= nrl;
 
	for(i=nrl;i<=nrh;i++) {
		m[i]=(unsigned char *) malloc((unsigned) (nch-ncl+1)*sizeof(unsigned char));
		if (!m[i]) dump_error("allocation failure 2 in cmatrix()");
		m[i] -= ncl;
	}
	return m;
}


 
void free_cmatrix(char **m,int nrl,int nrh,int ncl,int nch)
{
	int i;
 
	for(i=nrh;i>=nrl;i--) free((void *) (m[i]+ncl));
	free((void *) (m+nrl));
}

 
void free_ucmatrix(unsigned char **m,int nrl,int nrh,int ncl,int nch)
{
	int i;
 
	for(i=nrh;i>=nrl;i--) free((void *) (m[i]+ncl));
	free((void *) (m+nrl));
}
 
void free_cvector(char *v, int nl, int nh)
{
	free((void *) (v+nl));
}

void free_ucvector(unsigned char *v, int nl, int nh)
{
	free((void *) (v+nl));
}

 
void free_ivector(int *v, int nl, int nh)
{
	free((void *) (v+nl));
}

void free_ulvector(unsigned long *v, int nl, int nh)
{
	free((void *) (v+nl));
}

void free_uvector(unsigned int *v, int nl, int nh)
{
	free((void *) (v+nl));
}
 
void free_dvector(double *v, int nl,int nh)
{
	free((void *) (v+nl));
}
 
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/

void dump_error(char *error_text)
{
 
	fprintf(outfile,"PEmapper error...\n");
	fprintf(outfile,"%s\n",error_text);
	fprintf(outfile,"...now exiting to system...\n");
	exit(1);
}

/*---------------------------------------------------------------------*/



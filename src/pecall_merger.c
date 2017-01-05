/* 

The code itself is Copyright (C) 2017, by David J. Cutler.

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

typedef struct base_node
{
  short chrom;
  unsigned int pos;
  char *calls;
  unsigned char *quality;
  char ref;
  short no_calls;
} BNODE;

typedef struct sample_node
{
  char name[128];
  int which;
} SNODE;



void read_var (char *line, char *result);
int *ivector (int nl, int nh);
unsigned int *uvector (int nl, int nh);
char *cvector (int nl, int nh);
unsigned char *ucvector (int nl, int nh);
double *dvector (int nl, int nh);
double **dmatrix (int nrl, int nrh, int ncl, int nch);
unsigned short **usmatrix (int nrl, int nrh, int ncl, int nch);
int **imatrix (int nrl, int nrh, int ncl, int nch);
char **cmatrix (int nrl, int nrh, int ncl, int nch);
void free_cvector (char *v, int nl, int nh);
void free_ucvector (unsigned char *v, int nl, int nh);
void free_ivector (int *v, int nl, int nh);
void free_uvector (unsigned int *v, int nl, int nh);
void free_dvector (double *v, int nl, int nh);
void free_dmatrix (double **m, int nrl, int nrh, int ncl, int nch);
unsigned int **umatrix (int nrl, int nrh, int ncl, int nch);
void free_imatrix (int **m, int nrl, int nrh, int ncl, int nch);
void free_ucmatrix (unsigned char **m, int nrl, int nrh, int ncl, int nch);
void free_cmatrix (char **m, int nrl, int nrh, int ncl, int nch);
unsigned char **ucmatrix (int nrl, int nrh, int ncl, int nch);
void dump_error (char *error_text);
SNODE *sample_alloc ();
BNODE *base_alloc ();
int find_chrom (char **list, int n, char *s);
int compare_base (const void *a, const void *b);
void get_het_alleles (int i, int *a, int *b, int ref);

#define HARD_N 15
#define SOFT_N 14

#define minim(atesta,btestb) ((atesta<btestb)?atesta:btestb)
#define maxim(atesta,btestb) ((atesta>btestb)?atesta:btestb)


int
main (int argc, char *argv[])
{
  FILE *outfile, *bedfile, *sfile;
  gzFile basefile;
  int maxsnps;
  int i, j;
  int maxsamples;
  SNODE **samples;
  BNODE **bases;
  int last_base;
  int last_sample;
  char *temp_string;
  int char_to_int[256];
  char sss[4096];
  char **contig_names;
  int no_contigs;
  char int_to_char[256];
  char merge_ints[256][256];
  DIR *pDIR;
  struct dirent *pDirEnt;


  if (argc != 7)
  {
    printf ("\n Usage %s maxsnps maxsamples bedfile outfile sdxfile ishaploid \n", argv[0]);
    exit (1);
  }

  int IS_HAPLOID = FALSE;

  if ((strchr (argv[6], 'Y')) || (strchr (argv[6], 'y')))
    IS_HAPLOID = TRUE;
  if ((sfile = fopen (argv[5], "r")) == (FILE *) NULL)
  {
    printf ("\n Can not open file %s\n", argv[5]);
    exit (1);
  }
  if ((outfile = fopen (argv[4], "w")) == (FILE *) NULL)
  {
    printf ("\n Can not open file %s for writing\n", argv[4]);
    exit (1);
  }


  fgets (sss, 256, sfile);
  no_contigs = atoi (sss);
  contig_names = cmatrix (0, no_contigs, 0, 256);
  char *token;
  printf ("\n Reading Genome index \n");
  for (i = 0; i < no_contigs; i++)
  {
    // printf("\n About to read line %d \n\n",i);
    fgets (sss, 1024, sfile);
    token = strtok (sss, "\t \n");
    token = strtok (NULL, "\t \n");
    strcpy (contig_names[i], token);
  }
  fclose (sfile);
  maxsnps = atoi (argv[1]);
  maxsamples = atoi (argv[2]);

  printf ("\n Done reading Genome Index \n");
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
     14 = Soft N
     15 = Hard N
   */

  for (i = 0; i < 256; i++)
  {
    int_to_char[i] = 'N';
    char_to_int[i] = HARD_N;
    if (i == 14)		// Soft N
      for (j = 0; j < 256; j++)
	merge_ints[i][j] = j;
    else
    {
      for (j = 0; j < 256; j++)
	if (i != j)
	  merge_ints[i][j] = HARD_N;
	else
	  merge_ints[i][j] = i;

      merge_ints[i][14] = i;
    }
  }
  int_to_char[0] = 'A';
  int_to_char[1] = 'C';
  int_to_char[2] = 'G';
  int_to_char[3] = 'T';
  int_to_char[4] = 'D';
  int_to_char[5] = 'I';
  int_to_char[6] = 'M';
  int_to_char[7] = 'R';
  int_to_char[8] = 'W';
  int_to_char[9] = 'S';
  int_to_char[10] = 'Y';
  int_to_char[11] = 'K';
  int_to_char[12] = 'E';
  int_to_char[13] = 'H';

  char_to_int['A'] = char_to_int['a'] = 0;
  char_to_int['C'] = char_to_int['c'] = 1;
  char_to_int['G'] = char_to_int['g'] = 2;
  char_to_int['T'] = char_to_int['t'] = 3;
  char_to_int['D'] = char_to_int['d'] = 4;
  char_to_int['I'] = char_to_int['i'] = 5;
  char_to_int['M'] = char_to_int['m'] = 6;
  char_to_int['R'] = char_to_int['r'] = 7;
  char_to_int['W'] = char_to_int['w'] = 8;
  char_to_int['S'] = char_to_int['s'] = 9;
  char_to_int['Y'] = char_to_int['y'] = 10;
  char_to_int['K'] = char_to_int['k'] = 11;
  char_to_int['E'] = char_to_int['e'] = 12;
  char_to_int['H'] = char_to_int['h'] = 13;
  char_to_int['N'] = char_to_int['n'] = 14;


  samples = (SNODE **) malloc (sizeof (SNODE *) * maxsamples);
  if (!samples)
    dump_error ("\n Not enough memory for samples array \n");

  bases = (BNODE **) malloc (sizeof (BNODE *) * maxsnps);
  if (!bases)
    dump_error ("\n Not enough memory for bases array \n");

  last_base = 0;
  last_sample = 0;
  temp_string = cvector (0, 99999999);

  if ((bedfile = fopen (argv[3], "r")) == (FILE *) NULL)
  {
    printf ("\n Can not open %s which should be the bedfile of positions to report \n", argv[3]);
    exit (1);
  }
  int not_done = TRUE;
  fgets (temp_string, 9999999, bedfile);
  char last_chrom_name[64];
  strcpy (last_chrom_name, "!!!!!!");
  int last_chrom_number = -1;
  printf ("\n Reading Bed file \n");
  while (not_done)
  {
    char *token;
    int i;
    int start, stop;
    int chrom_no;
    token = strtok (temp_string, "\t \n");
    // printf("\n Working on chromosome %s \n",token);
    if (strcmp (last_chrom_name, token) == 0)
      chrom_no = last_chrom_number;
    else
    {
      chrom_no = find_chrom (contig_names, no_contigs, token);
      strcpy (last_chrom_name, token);
      last_chrom_number = chrom_no;
    }
    token = strtok (NULL, "\t \n");
    start = atoi (token);
    token = strtok (NULL, "\t \n");
    stop = atoi (token);

    int j = last_base + 1 + (stop - start);
    if ((j >= maxsnps) || (j < 0))
    {
      printf ("\n the line which is %s %d %d will cause us to have too many bases\n", last_chrom_name, start, stop);
      exit (1);
    }
    for (i = start; i <= stop; i++)
    {
      bases[last_base] = base_alloc ();
      bases[last_base]->chrom = chrom_no;
      bases[last_base]->pos = i;
      last_base++;
    }
    not_done = FALSE;
    temp_string[0] = '\0';
    if (!feof (bedfile))
    {
      fgets (temp_string, 9999999, bedfile);
      if (strlen (temp_string) > 3)
	not_done = TRUE;
    }
  }
  fclose (bedfile);
  printf ("\n Done reading bedfile \n");

  pDIR = opendir (".");

  if (pDIR == NULL)
  {
    fprintf (stderr, "%s %d: opendir() failed (%s)\n", __FILE__, __LINE__, strerror (errno));
    exit (-1);
  }

  temp_string = cvector (0, 99999999);

  pDirEnt = readdir (pDIR);
  printf ("\n Just read dir \n\n");
  BNODE *this_base;
  this_base = base_alloc ();
  int *sample_map;
  sample_map = ivector (0, maxsamples);
  double thres = 0.95;
  unsigned char char_thres = (unsigned char) round (thres * 255);
  printf ("\n The char_thres = %d \n", (int) char_thres);
  while (pDirEnt != NULL)
  {
    // printf("\n About to look at %s \n",pDirEnt->d_name);
    if ((strstr (pDirEnt->d_name, "base.gz") != NULL))
    {

      if ((basefile = gzopen (pDirEnt->d_name, "r")) == (gzFile) NULL)
      {
	printf ("\n Can not open file %s which is the basefile\n", pDirEnt->d_name);
	exit (1);
      }
      printf ("\n Working on file %s \n\n", pDirEnt->d_name);
      gzbuffer (basefile, 33554432);

      int not_done = TRUE;
      gzgets (basefile, temp_string, 99999999);
      char *token;
      // printf("\n Just grabbed %s \n\n",temp_string);
      token = strtok (temp_string, "\t \n");
      token = strtok (NULL, "\t \n");
      token = strtok (NULL, "\t \n");
      token = strtok (NULL, "\t \n");

      int this_sample = 0;
      while (token)
      {
	// printf("\n I've got sample %s \n\n",token);
	if (strlen (token) > 2)
	{
	  int old = FALSE;
	  for (i = 0; i < last_sample; i++)
	    if (strcmp (samples[i]->name, token) == 0)
	    {
	      printf ("\n We matched %s with %s \n\n", samples[i]->name, token);
	      sample_map[this_sample] = i;
	      i = last_sample;
	      old = TRUE;
	    }
	  if (!old)
	  {
	    samples[last_sample] = sample_alloc ();
	    strcpy (samples[last_sample]->name, token);
	    samples[last_sample]->which = last_sample;
	    sample_map[this_sample] = last_sample;
	    last_sample++;
	  }
	  this_sample++;
	}
	token = strtok (NULL, "\t \n");
      }
      printf ("\n Found a total of %d samples this time and a total of %d \n\n", this_sample, last_sample);
      char these_calls;
      float these_quals;

      gzgets (basefile, temp_string, 99999999);
      strcpy (last_chrom_name, "!!!!!!");

      last_chrom_number = -1;

      not_done = TRUE;
      while (not_done)
      {
	// printf("\n Line is %s \n\n",temp_string);
	token = strtok (temp_string, "\t \n");
	// printf("\n ABout to look for %s \n\n",token);
	if (strcmp (token, last_chrom_name) == 0)
	  this_base->chrom = last_chrom_number;
	else
	  this_base->chrom = find_chrom (contig_names, no_contigs, token);
	last_chrom_number = this_base->chrom;
	// printf("\n Chromosome number is %d\n\n",this_base->chrom);
	token = strtok (NULL, "\t \n");
	this_base->pos = atoi (token);
	token = strtok (NULL, "\t \n");
	this_base->ref = token[0];
	BNODE **look_base;
	// printf("\n About to bsearch \n\n");
	look_base = bsearch (&this_base, bases, last_base, sizeof (BNODE *), compare_base);
	// printf("\n Back from bsearch \n\n");
	if (look_base)
	{
	  BNODE *need_base;
	  need_base = *look_base;
	  // if(last_sample > 400)
	  //printf("\n Found something to store which goes in %d %d at memory spot %ld with start %ld so this element %ld  last_sample = %d \n\n",
	  //       need_base->chrom,need_base->pos,(long)look_base,(long)bases,(long)((long)look_base-(long)bases)/sizeof(BNODE *),last_sample);
	  // printf("\n testing 89 chrom = %d  pos = %d \n\n",bases[89]->chrom,bases[89]->pos);
	  need_base->chrom = this_base->chrom;
	  need_base->ref = this_base->ref;
	  need_base->pos = this_base->pos;

	  // printf("\n Stored some stuff \n\n");

	  char *temp_calls;
	  unsigned char *temp_quality;
	  temp_calls = cvector (0, last_sample);
	  temp_quality = ucvector (0, last_sample);
	  for (i = 0; i < last_sample; i++)
	  {
	    temp_calls[i] = 14;
	    temp_quality[i] = 0;
	  }
	  if (need_base->no_calls > 0)
	  {
	    for (i = 0; i < need_base->no_calls; i++)
	    {
	      temp_calls[i] = need_base->calls[i];
	      temp_quality[i] = need_base->quality[i];
	    }
	    if (need_base->calls)
	      free_cvector (need_base->calls, 0, need_base->no_calls);
	    if (need_base->quality)
	      free_ucvector (need_base->quality, 0, need_base->no_calls);
	  }
	  need_base->calls = temp_calls;
	  need_base->quality = temp_quality;
	  need_base->no_calls = last_sample;
	  // printf("\n About to loop through samples with max = %d \n\n",this_sample);
	  for (i = 0; i < this_sample; i++)
	  {
	    token = strtok (NULL, "\t \n");
	    these_calls = char_to_int[(int) token[0]];
	    token = strtok (NULL, "\t \n");
	    these_quals = atof (token);

	    unsigned char qual_char = (unsigned char) round (these_quals * 255);
	    // if(last_sample > 400)
	    // printf("\n i=%d call = %c qual = %g qual char = %d which will go in position %d \n\n",i,these_calls,these_quals,(int)qual_char,sample_map[i]);
	    if ((qual_char >= char_thres))
	    {
	      if (need_base->quality[sample_map[i]] >= char_thres)
	      {
		need_base->quality[sample_map[i]] = maxim (need_base->quality[sample_map[i]], qual_char);
		need_base->calls[sample_map[i]] = merge_ints[(int) need_base->calls[sample_map[i]]][(int) these_calls];
	      }
	      else
	      {
		need_base->quality[sample_map[i]] = qual_char;
		need_base->calls[sample_map[i]] = these_calls;
	      }
	    }
	    else if (qual_char > need_base->quality[sample_map[i]])
	    {
	      need_base->quality[sample_map[i]] = qual_char;
	      need_base->calls[sample_map[i]] = these_calls;
	    }
	  }
	}
	// printf("\n About to grab new line \n\n");
	not_done = FALSE;
	temp_string[0] = '\0';
	if (!gzeof (basefile))
	{
	  gzgets (basefile, temp_string, 99999999);
	  if (strlen (temp_string) > 10)
	    not_done = TRUE;
	}
      }
    }
    pDirEnt = readdir (pDIR);
  }
  printf ("\n About to write output file \n\n");
  fprintf (outfile, "Fragment\tPosition\tReference\tAlleles\tAllele_Counts\tType");
  for (i = 0; i < last_sample; i++)
    fprintf (outfile, "\t%s\t", samples[i]->name);

  int NO_ALLELES = 6;
  short allele_counts[4][HARD_N + 1][HARD_N + 1];
  int k;

  for (i = 0; i <= HARD_N; i++)
  {
    int a, b;
    for (j = 0; j < 4; j++)
    {
      for (k = 0; k < NO_ALLELES; k++)
	allele_counts[j][i][k] = 0;
      get_het_alleles (i, &a, &b, j);
      allele_counts[j][i][a]++;
      if (!IS_HAPLOID)
	allele_counts[j][i][b]++;
    }
  }

  char snptype[128];
  char count_string[128];
  char allele_string[128];
  char temps[128];
  for (i = 0; i < last_base; i++)
  {
    int this_allele_counts[NO_ALLELES];
    int no_alleles;
    int ref = char_to_int[(int) bases[i]->ref];
    for (k = 0; k < NO_ALLELES; k++)
      this_allele_counts[k] = 0;
    int this_s = minim (last_sample, bases[i]->no_calls);
    for (j = 0; j < this_s; j++)
      // if(bases[i]->quality[samples[j]->which] >= char_thres)
      for (k = 0; k < NO_ALLELES; k++)
	this_allele_counts[k] += allele_counts[ref][(int) bases[i]->calls[samples[j]->which]][k];
    no_alleles = 0;
    for (j = 0; j < NO_ALLELES; j++)
      if (this_allele_counts[j] > 0)
	no_alleles++;
    if (no_alleles > 0)
    {
      if (no_alleles > 2)
	sprintf (snptype, "MULTIALLELIC");
      else if (this_allele_counts[4] > 0)
	sprintf (snptype, "DEL");
      else if (this_allele_counts[5] > 0)
	sprintf (snptype, "INS");
      else
	sprintf (snptype, "SNP");
      int thisa = 0;

      for (j = 0; j < NO_ALLELES; j++)
	if (this_allele_counts[j] > 0)
	{
	  if (thisa == 0)
	  {
	    sprintf (allele_string, "%c", int_to_char[j]);
	    sprintf (count_string, "%d", this_allele_counts[j]);
	  }
	  else
	  {
	    sprintf (temps, ",%c", int_to_char[j]);
	    strcat (allele_string, temps);
	    sprintf (temps, ",%d", this_allele_counts[j]);
	    strcat (count_string, temps);
	  }
	  thisa++;
	}
      fprintf (outfile, "\n%s\t%d\t%c\t%s\t%s\t%s",
	       contig_names[bases[i]->chrom], bases[i]->pos, bases[i]->ref, allele_string, count_string, snptype);
      for (j = 0; j < this_s; j++)
	fprintf (outfile, "\t%c\t%g",
		 int_to_char[(int) bases[i]->calls[samples[j]->which]],
		 ((float) bases[i]->quality[samples[j]->which]) / 255.0);
      for (j = this_s; j < last_sample; j++)
	fprintf (outfile, "\tN\t1");
    }
  }
  fprintf (outfile, "\n");


  fclose (outfile);
  return 0;

}

/*---------------------------------------------------------------------*/
void
get_het_alleles (int i, int *a, int *b, int ref)
{
  if (i < 6)
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
    *a = 14;
    *b = 14;
  }

  return;
}

/*---------------------------------------------------------------------*/
int
compare_base (const void *a, const void *b)
{
  BNODE *A, *B;
  A = *((BNODE **) a);
  B = *((BNODE **) b);

  if (A->chrom < B->chrom)
    return -1;
  if (A->chrom > B->chrom)
    return 1;

  if (A->pos < B->pos)
    return -1;
  if (A->pos > B->pos)
    return 1;

  return 0;
}

/*---------------------------------------------------------------------*/
int
find_chrom (char **list, int n, char *s)
{
  int i;
  for (i = 0; i < n; i++)
    if (strcmp (list[i], s) == 0)
      return i;
  printf ("\n We could not find chromosome %s \n", s);
  exit (1);
  return 0;
}

/*---------------------------------------------------------------------*/

SNODE *
sample_alloc ()
{
  SNODE *tn;

  tn = (SNODE *) malloc ((unsigned) sizeof (struct sample_node));
  if (!tn)
    dump_error ("allocation failure in sample_alloc()");

  tn->name[0] = '\0';
  tn->which = -1;
  return tn;
}

/*---------------------------------------------------------------------*/

BNODE *
base_alloc ()
{
  BNODE *tn;
  tn = (BNODE *) malloc ((unsigned) sizeof (struct base_node));
  if (!tn)
    dump_error ("allocation failure in base_alloc()");

  tn->chrom = 0;
  tn->pos = 0;
  tn->calls = NULL;
  tn->quality = NULL;
  tn->ref = 'N';
  tn->no_calls = 0;

  return tn;
}

/*---------------------------------------------------------------------*/


char *
cvector (int nl, int nh)
{
  char *v;

  v = (char *) malloc ((unsigned) (nh - nl + 1) * sizeof (char));
  if (!v)
    dump_error ("allocation failure in cvector()");
  return v - nl;
}

unsigned char *
ucvector (int nl, int nh)
{
  unsigned char *v;

  v = (unsigned char *) malloc ((unsigned) (nh - nl + 1) * sizeof (unsigned char));
  if (!v)
    dump_error ("allocation failure in ucvector()");
  return v - nl;
}

int *
ivector (int nl, int nh)
{
  int *v;

  v = (int *) malloc ((unsigned) (nh - nl + 1) * sizeof (int));
  if (!v)
    dump_error ("allocation failure in ivector()");
  return v - nl;
}

unsigned int *
uvector (int nl, int nh)
{
  unsigned int *v;

  v = (unsigned int *) malloc ((unsigned) (nh - nl + 1) * sizeof (int));
  if (!v)
    dump_error ("allocation failure in uvector()");
  return v - nl;
}

double *
dvector (int nl, int nh)
{
  double *v;

  v = (double *) malloc ((unsigned) (nh - nl + 1) * sizeof (double));
  if (!v)
    dump_error ("allocation failure in dvector()");
  return v - nl;
}


int **
imatrix (int nrl, int nrh, int ncl, int nch)
{
  int i, **m;

  m = (int **) malloc ((unsigned) (nrh - nrl + 1) * sizeof (int *));
  if (!m)
    dump_error ("allocation failure 1 in imatrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++)
  {
    m[i] = (int *) malloc ((unsigned) (nch - ncl + 1) * sizeof (int));
    if (!m[i])
      dump_error ("allocation failure 2 in imatrix()");
    m[i] -= ncl;
  }
  return m;
}

void
free_imatrix (int **m, int nrl, int nrh, int ncl, int nch)
{
  int i;

  for (i = nrh; i >= nrl; i--)
    free ((char *) (m[i] + ncl));
  free ((char *) (m + nrl));
}


double **
dmatrix (int nrl, int nrh, int ncl, int nch)
{
  int i;
  double **m;

  m = (double **) malloc ((unsigned) (nrh - nrl + 1) * sizeof (double *));
  if (!m)
    dump_error ("allocation failure 1 in dmatrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++)
  {
    m[i] = (double *) malloc ((unsigned) (nch - ncl + 1) * sizeof (double));
    if (!m[i])
      dump_error ("allocation failure 2 in dmatrix()");
    m[i] -= ncl;
  }
  return m;
}

unsigned short **
usmatrix (int nrl, int nrh, int ncl, int nch)
{
  int i;
  unsigned short **m;

  m = (unsigned short **) malloc ((unsigned) (nrh - nrl + 1) * sizeof (unsigned short *));
  if (!m)
    dump_error ("allocation failure 1 in cmatrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++)
  {
    m[i] = (unsigned short *) malloc ((unsigned) (nch - ncl + 1) * sizeof (unsigned short));
    if (!m[i])
      dump_error ("allocation failure 2 in cmatrix()");
    m[i] -= ncl;
  }
  return m;
}

unsigned int **
umatrix (int nrl, int nrh, int ncl, int nch)
{
  int i;
  unsigned int **m;

  m = (unsigned int **) malloc ((unsigned) (nrh - nrl + 1) * sizeof (unsigned int *));
  if (!m)
    dump_error ("allocation failure 1 in cmatrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++)
  {
    m[i] = (unsigned int *) malloc ((unsigned) (nch - ncl + 1) * sizeof (unsigned int));
    if (!m[i])
      dump_error ("allocation failure 2 in cmatrix()");
    m[i] -= ncl;
  }
  return m;
}

unsigned char **
ucmatrix (int nrl, int nrh, int ncl, int nch)
{
  int i;
  unsigned char **m;

  m = (unsigned char **) malloc ((unsigned) (nrh - nrl + 1) * sizeof (unsigned char *));
  if (!m)
    dump_error ("allocation failure 1 in cmatrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++)
  {
    m[i] = (unsigned char *) malloc ((unsigned) (nch - ncl + 1) * sizeof (unsigned char));
    if (!m[i])
      dump_error ("allocation failure 2 in cmatrix()");
    m[i] -= ncl;
  }
  return m;
}

void
free_dmatrix (double **m, int nrl, int nrh, int ncl, int nch)
{
  int i;

  for (i = nrh; i >= nrl; i--)
    free ((char *) (m[i] + ncl));
  free ((double *) (m + nrl));
}

char **
cmatrix (int nrl, int nrh, int ncl, int nch)
{
  int i;
  char **m;

  m = (char **) malloc ((unsigned) (nrh - nrl + 1) * sizeof (char *));
  if (!m)
    dump_error ("allocation failure 1 in cmatrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++)
  {
    m[i] = (char *) malloc ((unsigned) (nch - ncl + 1) * sizeof (char));
    if (!m[i])
      dump_error ("allocation failure 2 in cmatrix()");
    m[i] -= ncl;
  }
  return m;
}

void
free_cmatrix (char **m, int nrl, int nrh, int ncl, int nch)
{
  int i;

  for (i = nrh; i >= nrl; i--)
    free ((char *) (m[i] + ncl));
  free ((char *) (m + nrl));
}

void
free_ucmatrix (unsigned char **m, int nrl, int nrh, int ncl, int nch)
{
  int i;

  for (i = nrh; i >= nrl; i--)
    free ((unsigned char *) (m[i] + ncl));
  free ((unsigned char *) (m + nrl));
}

void
free_cvector (char *v, int nl, int nh)
{
  free ((char *) (v + nl));
}

void
free_ucvector (unsigned char *v, int nl, int nh)
{
  free ((unsigned char *) (v + nl));
}

void
free_ivector (int *v, int nl, int nh)
{
  free ((int *) (v + nl));
}

void
free_uvector (unsigned int *v, int nl, int nh)
{
  free ((int *) (v + nl));
}

void
free_dvector (double *v, int nl, int nh)
{
  free ((double *) (v + nl));
}

/*---------------------------------------------------------------------*/

void
dump_error (char *ss)
{
  printf ("\n Unrecoverable error\n %s \n  Exiting now \n", ss);
  exit (1);
}

/* -------------------------------------------------------------- */

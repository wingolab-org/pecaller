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

int *ivector (int nl, int nh);
char *cvector (int nl, int nh);
unsigned char *ucvector (int nl, int nh);
unsigned long *ulvector (int nl, int nh);
long *lvector (int nl, int nh);
unsigned int *uvector (int nl, int nh);
unsigned long long *ullvector (int nl, int nh);
double *dvector (int nl, int nh);
double **dmatrix (int nrl, int nrh, int ncl, int nch);
char **cmatrix (int nrl, int nrh, int ncl, int nch);
unsigned char **ucmatrix (int nrl, int nrh, int ncl, int nch);
unsigned long **ulmatrix (int nrl, int nrh, int ncl, int nch);
unsigned int **umatrix (int nrl, int nrh, int ncl, int nch);
void free_cvector (char *v, int nl, int nh);
void free_ivector (int *v, int nl, int nh);
void free_ucvector (unsigned char *v, int nl, int nh);
void free_ulvector (unsigned long *v, int nl, int nh);
void free_uvector (unsigned int *v, int nl, int nh);
void free_dvector (double *v, int nl, int nh);
void free_dmatrix (double **m, int nrl, int nrh, int ncl, int nch);
void free_cmatrix (char **m, int nrl, int nrh, int ncl, int nch);
void free_ucmatrix (unsigned char **m, int nrl, int nrh, int ncl, int nch);
void free_ulmatrix (unsigned long **m, int nrl, int nrh, int ncl, int nch);
void free_umatrix (unsigned int **m, int nrl, int nrh, int ncl, int nch);
void dump_error (char *error_text);
int **imatrix (int nrl, int nrh, int ncl, int nch);
void free_imatrix (int **m, int nrl, int nrh, int ncl, int nch);
int find_chrom (unsigned int *pos, int first, int last, int try, unsigned this);


static FILE *outfile;
static long genome_size;
static int MAX_FILE_BUFFER = 512000000;

int
main (int argc, char *argv[])
{
  char sss[4196];
  char sdxname[1024];
  char **contig_names, *token;
  unsigned int *contig_starts;
  int max_contigs, no_contigs;
  char *genome_buffer;
  double min_prob;
  int i, j;


  FILE *sfile;
  gzFile reffile, snpfile;


  if ((argc != 4) && (argc != 3))
  {
    printf ("\nUsage: %s sdx_file snpfile [min_prob_to_make_call] \n", argv[0]);
    exit (1);
  }
  min_prob = 0.0;
  double tp = (double) atof (argv[3]);
  if ((tp >= 0.0) && (tp <= 1.0))
    min_prob = tp;

  strcpy (sdxname, argv[1]);
  if ((sfile = fopen (sdxname, "r")) == (FILE *) NULL)
  {
    printf ("\n Can not open file %s\n", sdxname);
    exit (1);
  }
  if (strstr (sdxname, ".sdx") != NULL)
  {
    for (i = strlen (sdxname) - 1; i > 0; i--)
      if (sdxname[i] == '.')
      {
	sdxname[i] = '\0';
	i = 0;
      }
  }


  fgets (sss, 256, sfile);
  max_contigs = atoi (sss);
  no_contigs = max_contigs;
  contig_starts = uvector (0, max_contigs);
  contig_names = cmatrix (0, max_contigs, 0, 256);

  contig_starts[0] = 0;
  for (i = 0; i < max_contigs; i++)
  {
    fgets (sss, 1024, sfile);
    token = strtok (sss, "\t \n");
    contig_starts[i + 1] = atoi (token);
    token = strtok (NULL, "\t \n");
    strcpy (contig_names[i], token);
    // printf("\nFor contig %d we have offset %d",i,contig_starts[i]); 
  }
  fgets (sss, 1024, sfile);
  fclose (sfile);
  for (i = 1; i <= max_contigs; i++)
    contig_starts[i] += contig_starts[i - 1];

  // for(i=0;i<=max_contigs;i++)
  // printf("\n Contig %d starts at position %u \n",i,contig_starts[i]);

  genome_size = (long) contig_starts[max_contigs] + (long) 15 *max_contigs;
  genome_buffer = (char *) malloc (sizeof (char) * (genome_size + 1));
  if (!genome_buffer)
    dump_error ("\n Failed to allocate memory for the Genome Buffer \n");
  // else
  // printf("\n Genome size is %ld \n\n",genome_size);
  sprintf (sss, "%s.seq", sdxname);
  if ((reffile = gzopen (sss, "r")) == (gzFile) NULL)
  {
    printf ("\n Can not open file %s for reading\n", sss);
    exit (1);
  }
  gzbuffer (reffile, 33554432);
  long g_temp = 0;
  // printf("\n About to read genome \n\n");
  if (genome_size < MAX_FILE_BUFFER)
    g_temp = gzread (reffile, (void *) genome_buffer, (unsigned int) sizeof (char) * genome_size);
  else
  {
    long count = 0;
    while (count < genome_size)
    {
      int ttemp = minim ((long) genome_size - count, MAX_FILE_BUFFER);
      g_temp += gzread (reffile, (void *) &genome_buffer[count], ttemp);
      count += ttemp;
    }
  }
  gzclose (reffile);

  for (i = 1, j = 15; i <= no_contigs; i++, j += 15)
    contig_starts[i] += j;

  printf ("##fileformat=VCFv4.0\n");
  time_t ttt = time (NULL);
  struct tm tm = *localtime (&ttt);
  printf ("##fileDate=%d%d%d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday);
  printf ("##reference=%s\n", argv[1]);
  printf ("##phasing=none");
  printf ("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n");
  printf ("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
  printf ("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
  printf ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
  if ((snpfile = gzopen (argv[2], "r")) == (gzFile) NULL)
  {
    printf ("\n Can not open file %s for reading\n", argv[2]);
    exit (1);
  }
  char *buffer;
  int buff_len = 41959999;
  buffer = cvector (0, buff_len + 1);
  gzbuffer (snpfile, 33554432);
  gzgets (snpfile, buffer, buff_len);
  int tot_samples = 0;

  token = strtok (buffer, "\n\t ");
  for (i = 0; i < 6; i++)
    token = strtok (NULL, "\n\t ");
  while (token)
  {
    printf ("\t%s", token);
    tot_samples++;
    token = strtok (NULL, "\n\t ");
    // token = strtok(NULL,"\n\t ");
  }
  gzgets (snpfile, buffer, buff_len);
  int len = strlen (buffer);
  char chrom[128];
  char ref;
  char ref_string[8196];
  char alt_a_temp[8196];
  char alt_a_final[8196];
  int pos;
  int drop_it;
  char **call_map;
  char **het_map;
  call_map = cmatrix (0, 255, 0, 10);
  for (j = 0; j < 256; j++)
    strcpy (call_map[j], "./.");
  char slabel[255];
  het_map = cmatrix (0, 255, 0, 255);
  for (i = 0; i < 256; i++)
    for (j = 0; j < 256; j++)
      het_map[i][j] = 'N';
  het_map[(int) 'A'][(int) 'C'] = 'M';
  het_map[(int) 'A'][(int) 'G'] = 'R';
  het_map[(int) 'A'][(int) 'T'] = 'W';
  het_map[(int) 'A'][(int) 'D'] = 'E';
  het_map[(int) 'A'][(int) 'I'] = 'H';

  het_map[(int) 'C'][(int) 'G'] = 'S';
  het_map[(int) 'C'][(int) 'T'] = 'Y';
  het_map[(int) 'C'][(int) 'A'] = 'M';
  het_map[(int) 'C'][(int) 'D'] = 'E';
  het_map[(int) 'C'][(int) 'I'] = 'H';

  het_map[(int) 'G'][(int) 'T'] = 'K';
  het_map[(int) 'G'][(int) 'A'] = 'R';
  het_map[(int) 'G'][(int) 'C'] = 'S';
  het_map[(int) 'G'][(int) 'D'] = 'E';
  het_map[(int) 'G'][(int) 'I'] = 'H';

  het_map[(int) 'T'][(int) 'A'] = 'W';
  het_map[(int) 'T'][(int) 'C'] = 'Y';
  het_map[(int) 'T'][(int) 'G'] = 'K';
  het_map[(int) 'T'][(int) 'D'] = 'E';
  het_map[(int) 'T'][(int) 'I'] = 'H';

  het_map[(int) 'D'][(int) 'A'] = 'E';
  het_map[(int) 'D'][(int) 'C'] = 'E';
  het_map[(int) 'D'][(int) 'G'] = 'E';
  het_map[(int) 'D'][(int) 'T'] = 'E';
  het_map[(int) 'D'][(int) 'I'] = 'E';

  het_map[(int) 'I'][(int) 'A'] = 'H';
  het_map[(int) 'I'][(int) 'C'] = 'H';
  het_map[(int) 'I'][(int) 'G'] = 'H';
  het_map[(int) 'I'][(int) 'T'] = 'H';
  het_map[(int) 'I'][(int) 'D'] = 'H';


  char allele_char[30];
  for (i = 0; i < 30; i++)
    allele_char[i] = 'N';
  char last_chr[128];
  int last_chr_no = 0;
  sprintf (last_chr, "!!!!!!");

  while (len > 5)
  {
    token = strtok (buffer, "\n\t ");
    strcpy (chrom, token);

    token = strtok (NULL, "\n\t ");
    pos = atoi (token);
    token = strtok (NULL, "\n\t ");
    ref = token[0];
    token = strtok (NULL, "\n\t ");
    strcpy (alt_a_temp, token);
    token = strtok (NULL, "\n\t ");
    token = strtok (NULL, "\n\t ");
    drop_it = FALSE;
    if (strcmp (token, "LOW") == 0)
      drop_it = TRUE;
    else if (strcmp (token, "MESS") == 0)
      drop_it = TRUE;


    if (!drop_it)
    {
      strcpy (call_map[(int) 'A'], "1/1");
      strcpy (call_map[(int) 'C'], "1/1");
      strcpy (call_map[(int) 'G'], "1/1");
      strcpy (call_map[(int) 'T'], "1/1");
      strcpy (call_map[(int) 'I'], "1/1");
      strcpy (call_map[(int) 'D'], "1/1");
      strcpy (call_map[(int) 'Y'], "0/1");
      strcpy (call_map[(int) 'R'], "0/1");
      strcpy (call_map[(int) 'S'], "0/1");
      strcpy (call_map[(int) 'W'], "0/1");
      strcpy (call_map[(int) 'K'], "0/1");
      strcpy (call_map[(int) 'M'], "0/1");
      strcpy (call_map[(int) 'E'], "0/1");
      strcpy (call_map[(int) 'H'], "0/1");
      strcpy (call_map[(int) ref], "0/0");
      strcpy (slabel, "PASS");
      sprintf (ref_string, "%c", ref);
      allele_char[0] = ref;

      if (strcmp (token, "SNP") == 0)
      {
	if (alt_a_temp[0] == ref)
	{
	  sprintf (alt_a_final, "%c", alt_a_temp[2]);
	  strcpy (call_map[(int) alt_a_temp[2]], "1/1");
	  allele_char[1] = alt_a_temp[2];
	}
	else
	{
	  sprintf (alt_a_final, "%c", alt_a_temp[0]);
	  strcpy (call_map[(int) alt_a_temp[0]], "1/1");
	  allele_char[1] = alt_a_temp[0];
	}
      }
      else if (strcmp (token, "MULTIALLELIC") == 0)
      {
	int this_a = 1;
	int this_a_pos = 0;
	int has_del = FALSE;
	char alt_a = 'N';
	int this_stop = strlen (alt_a_temp);

	while (this_a_pos < this_stop)
	{
	  if (alt_a_temp[this_a_pos] == ref)
	    this_a_pos += 2;
	  else if (alt_a_temp[this_a_pos] == '+')
	  {
	    alt_a = 'I';
	    allele_char[this_a] = alt_a;
	    sprintf (call_map[(int) alt_a], "%d/%d", this_a, this_a);
	    sprintf (call_map[(int) 'H'], "0/%d", this_a);
	    if (!has_del)
	    {
	      if (this_a == 1)
		sprintf (alt_a_final, "%c", ref);
	      else
		sprintf (alt_a_final, "%s,%c", alt_a_final, ref);
	    }
	    else
	      sprintf (alt_a_final, "%s,%s", alt_a_final, ref_string);

	    this_a_pos++;
	    while ((this_a_pos < this_stop) && (alt_a_temp[this_a_pos] != ','))
	    {
	      if (isalpha (alt_a_temp[this_a_pos]))
		sprintf (alt_a_final, "%s%c", alt_a_final, alt_a_temp[this_a_pos]);
	      this_a_pos++;
	    }

	    this_a_pos++;
	    this_a++;
	    sprintf (slabel, ".");
	  }
	  else if (alt_a_temp[this_a_pos] == '-')
	  {
	    alt_a = 'D';
	    allele_char[this_a] = alt_a;
	    sprintf (call_map[(int) alt_a], "%d/%d", this_a, this_a);
	    sprintf (call_map[(int) 'E'], "0/%d", this_a);

	    if (strcmp (chrom, last_chr) != 0)
	    {
	      last_chr_no = -1;
	      for (i = 0; i < no_contigs; i++)
		if (strcmp (chrom, contig_names[i]) == 0)
		{
		  strcpy (last_chr, chrom);
		  last_chr_no = i;
		  i = no_contigs;
		}
	      if (last_chr_no < 0)
	      {
		printf ("\n Failed to find chrom = %s \n", chrom);
		exit (1);
	      }
	    }
	    pos--;

	    long this_offset = pos + contig_starts[last_chr_no] - 1;
	    char sn[128];
	    has_del = TRUE;
	    ref = genome_buffer[this_offset];
	    this_a_pos++;
	    i = 0;
	    while ((this_a_pos < this_stop) && (alt_a_temp[this_a_pos] != ','))
	      sn[i++] = alt_a_temp[this_a_pos++];
	    sn[i] = '\0';
	    int del_len = atoi (sn);
	    del_len++;

	    char gb[4196];
	    strncpy (gb, &genome_buffer[this_offset], del_len);
	    gb[del_len] = '\0';
	    sprintf (ref_string, "%s", gb);
	    sprintf (slabel, ".");
	    if (this_a == 1)
	      sprintf (alt_a_final, "%c", ref);
	    else
	    {
	      strcpy (gb, alt_a_final);
	      strcpy (sn, ref_string);
	      sn[1] = gb[0];
	      sprintf (alt_a_final, "%s", sn);
	      for (i = 2, j = 2; i < this_a; i++, j += 2)
	      {
		strcpy (sn, ref_string);
		sn[1] = gb[j];
		sprintf (alt_a_final, "%s,%s", alt_a_final, sn);
	      }
	      sprintf (alt_a_final, "%s,%c", alt_a_final, ref);
	    }
	    this_a_pos++;
	    this_a++;
	    sprintf (slabel, ".");
	  }
	  else
	  {
	    alt_a = alt_a_temp[this_a_pos];
	    sprintf (call_map[(int) alt_a], "%d/%d", this_a, this_a);
	    allele_char[this_a] = alt_a;
	    for (i = 0; i <= this_a; i++)
	      for (j = i + 1; j <= this_a; j++)
		if (i != j)
		  sprintf (call_map[(int) het_map[(int) allele_char[i]][(int) allele_char[j]]], "%d/%d", i, j);
	    if (this_a == 1)
	      sprintf (alt_a_final, "%c", alt_a);
	    else
	      sprintf (alt_a_final, "%s,%c", alt_a_final, alt_a);
	    this_a++;
	    this_a_pos += 2;
	  }
	}
      }
      else if (strcmp (token, "INS") == 0)
      {
	int mono_allelic = TRUE;
	int al = strlen (alt_a_temp);
	for (i = 1; i <= al; i++)
	  if (alt_a_temp[i] == ',')
	    mono_allelic = FALSE;
	if (!mono_allelic)
	  sprintf (alt_a_final, "%c%s", ref, &alt_a_temp[3]);
	else
	  sprintf (alt_a_final, "%c%s", ref, &alt_a_temp[1]);
	sprintf (slabel, ".");
      }
      else			// Deletion
      {
	if (strcmp (chrom, last_chr) != 0)
	{
	  last_chr_no = -1;
	  for (i = 0; i < no_contigs; i++)
	    if (strcmp (chrom, contig_names[i]) == 0)
	    {
	      strcpy (last_chr, chrom);
	      last_chr_no = i;
	      i = no_contigs;
	    }
	  if (last_chr_no < 0)
	  {
	    printf ("\n Failed to find chrom = %s \n", chrom);
	    exit (1);
	  }
	}
	pos--;
	long this_offset = pos + contig_starts[last_chr_no] - 1;
	ref = genome_buffer[this_offset];
	char sn[128];
	int mono_allelic = TRUE;
	int al = strlen (alt_a_temp);
	for (i = 1; i <= al; i++)
	  if (alt_a_temp[i] == ',')
	    mono_allelic = FALSE;
	if (!mono_allelic)
	  sprintf (sn, "%s", &alt_a_temp[3]);
	else
	  sprintf (sn, "%s", &alt_a_temp[1]);

	int del_len = atoi (sn);
	del_len++;
	strncpy (ref_string, &genome_buffer[this_offset], del_len);
	ref_string[del_len] = '\0';
	sprintf (slabel, ".");
	sprintf (alt_a_final, "%c", ref);
      }
      printf ("\n%s\t%d\t.\t%s\t%s\t.\t%s\tNS=%d\tGT:GQ", chrom, pos, ref_string, alt_a_final, slabel, tot_samples);
      char *token2;
      double op;
      for (i = 0; i < tot_samples; i++)
      {
	token = strtok (NULL, "\n\t ");
	token2 = strtok (NULL, "\n\t ");
	op = (double) atof (token2);
	if (op >= min_prob)
	  printf ("\t%s", call_map[(int) token[0]]);
	else
	  printf ("\t./.");
	printf (":%s", token2);
      }
    }

    buffer[0] = '\0';
    if (!gzeof (snpfile))
      gzgets (snpfile, buffer, buff_len);
    len = strlen (buffer);
  }
  printf ("\n");

  exit (0);

}

/*-------------------------------------------------------------------------------------------------------------------------------------- */

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

unsigned long long *
ullvector (int nl, int nh)
{
  unsigned long long *v;

  v = (unsigned long long *) malloc ((unsigned) (nh - nl + 1) * sizeof (unsigned long long));
  if (!v)
    dump_error ("allocation failure in ullvector()");
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
    dump_error ("allocation failure in ulvector()");
  return v - nl;
}

long *
lvector (int nl, int nh)
{
  long *v;

  v = (long *) malloc ((nh - nl + 1) * sizeof (long));
  if (!v)
    dump_error ("allocation failure in lvector()");
  return v - nl;
}

unsigned long *
ulvector (int nl, int nh)
{
  unsigned long *v;

  v = (unsigned long *) malloc ((unsigned) (nh - nl + 1) * sizeof (long));
  if (!v)
    dump_error ("allocation failure in ulvector()");
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

unsigned int **
umatrix (int nrl, int nrh, int ncl, int nch)
{
  int i;
  unsigned int **m;

  m = (unsigned int **) malloc ((unsigned) (nrh - nrl + 1) * sizeof (unsigned int *));
  if (!m)
    dump_error ("allocation failure 1 in ulmatrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++)
  {
    m[i] = (unsigned int *) malloc ((unsigned) (nch - ncl + 1) * sizeof (unsigned int));
    if (!m[i])
      dump_error ("allocation failure 2 in ulmatrix()");
    m[i] -= ncl;
  }
  return m;
}

unsigned long **
ulmatrix (int nrl, int nrh, int ncl, int nch)
{
  int i;
  unsigned long **m;

  m = (unsigned long **) malloc ((unsigned) (nrh - nrl + 1) * sizeof (unsigned long *));
  if (!m)
    dump_error ("allocation failure 1 in ulmatrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++)
  {
    m[i] = (unsigned long *) malloc ((unsigned) (nch - ncl + 1) * sizeof (unsigned long));
    if (!m[i])
      dump_error ("allocation failure 2 in ulmatrix()");
    m[i] -= ncl;
  }
  return m;
}

void
free_umatrix (unsigned int **m, int nrl, int nrh, int ncl, int nch)
{
  int i;

  for (i = nrh; i >= nrl; i--)
    free ((void *) (m[i] + ncl));
  free ((void *) (m + nrl));
}


void
free_ulmatrix (unsigned long **m, int nrl, int nrh, int ncl, int nch)
{
  int i;

  for (i = nrh; i >= nrl; i--)
    free ((void *) (m[i] + ncl));
  free ((void *) (m + nrl));
}


void
free_imatrix (int **m, int nrl, int nrh, int ncl, int nch)
{
  int i;

  for (i = nrh; i >= nrl; i--)
    free ((void *) (m[i] + ncl));
  free ((void *) (m + nrl));
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

void
free_dmatrix (double **m, int nrl, int nrh, int ncl, int nch)
{
  int i;

  for (i = nrh; i >= nrl; i--)
    free ((void *) (m[i] + ncl));
  free ((void *) (m + nrl));
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
free_cmatrix (char **m, int nrl, int nrh, int ncl, int nch)
{
  int i;

  for (i = nrh; i >= nrl; i--)
    free ((void *) (m[i] + ncl));
  free ((void *) (m + nrl));
}


void
free_ucmatrix (unsigned char **m, int nrl, int nrh, int ncl, int nch)
{
  int i;

  for (i = nrh; i >= nrl; i--)
    free ((void *) (m[i] + ncl));
  free ((void *) (m + nrl));
}

void
free_cvector (char *v, int nl, int nh)
{
  free ((void *) (v + nl));
}

void
free_ucvector (unsigned char *v, int nl, int nh)
{
  free ((void *) (v + nl));
}


void
free_ivector (int *v, int nl, int nh)
{
  free ((void *) (v + nl));
}

void
free_ulvector (unsigned long *v, int nl, int nh)
{
  free ((void *) (v + nl));
}

void
free_uvector (unsigned int *v, int nl, int nh)
{
  free ((void *) (v + nl));
}

void
free_dvector (double *v, int nl, int nh)
{
  free ((void *) (v + nl));
}

/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/

void
dump_error (char *error_text)
{

  fprintf (outfile, "PEmapper error...\n");
  fprintf (outfile, "%s\n", error_text);
  fprintf (outfile, "...now exiting to system...\n");
  exit (1);
}

/*---------------------------------------------------------------------*/

/*

The code itself is Copyright (C) 2018, by David J. Cutler.

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
#include <fcntl.h>
#include <sys/mman.h>

#define FALSE 0
#define TRUE 1
char **cmatrix(int nrl, int nrh, int ncl, int nch);
unsigned int **umatrix(int nrl, int nrh, int ncl, int nch);
unsigned *uvector(int nl, int nh);
char *cvector(int nl, int nh);
void dump_error(char *error_text);
int find_chrom(char *this, char **contigs, int n);
FILE *outfile;

#define minim(atesta, btestb) ((atesta < btestb) ? atesta : btestb)
#define maxim(atesta, btestb) ((atesta > btestb) ? atesta : btestb)
#define pileup_header_size 240

int main(int argc, char *argv[])
{
  FILE *sfile;
  gzFile pileupfile;
  gzFile indelfile;
  int i, ii;
  int use_stdin;
  char sss[4096];
  char **contig_names;
  int no_contigs;
  char *pileup_header;
  outfile = stdout;

  pileupfile = NULL;
  indelfile = NULL;

  if (argc != 4)
  {
    printf("\n Usage %s reference_genome_sdxfile input_pileup[stdin,or filename] output_pileup \n", argv[0]);
    exit(1);
  }

  unsigned int *frag_pos;

  if ((sfile = fopen(argv[1], "r")) == (FILE *)NULL)
  {
    printf("\n Can not open file %s\n", argv[1]);
    exit(1);
  }

  fgets(sss, 256, sfile);
  no_contigs = atoi(sss);
  contig_names = cmatrix(0, no_contigs, 0, 256);
  frag_pos = uvector(-1, no_contigs);
  frag_pos[-1] = 0;
  char *token;
  printf("\n Reading Genome index \n");
  for (i = 0; i < no_contigs; i++)
  {
    // printf("\n About to read line %d \n\n",i);
    fgets(sss, 1024, sfile);
    token = strtok(sss, "\t \n");
    frag_pos[i] = atoi(token) + 15;
    frag_pos[i] += frag_pos[i - 1];
    token = strtok(NULL, "\t \n");
    strcpy(contig_names[i], token);
  }
  unsigned int genome_size = frag_pos[no_contigs - 1];
  fclose(sfile);

  printf("\n About to calculate some sizes for the bucket \n\n");

  unsigned int record = 6 * sizeof(unsigned short);
  unsigned short *bucket;
  size_t bucket_size = (size_t)record * (size_t)genome_size;
  long tot_recs = (long)6 * (long)genome_size;
  pileup_header = cvector(0, pileup_header_size + 1);
  for (i = 0; i <= pileup_header_size; i++)
    pileup_header[i] = '\0';

  strcpy(pileup_header, argv[1]);

  // exit(15);
  int fd;
  int result;
  size_t final_size = bucket_size + pileup_header_size;
  printf("\n About to open the memory map of size %ld \n\n", final_size);

  /* Open a file for writing.
   *      *  - Creating the file if it doesn't exist.
   *           *  - Truncating it to 0 size if it already exists. (not really needed)
   *                *
   *                     * Note: "O_WRONLY" mode is not sufficient when mmaping.
   *                          */

  fd = open(argv[3], O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);
  if (fd == -1)
  {
    printf("Failure to open %s for writing", argv[3]);
    exit(2);
  }
  write(fd, pileup_header, pileup_header_size);
  result = lseek(fd, final_size - 1, SEEK_SET);
  if (result == -1)
  {
    close(fd);
    printf("Error calling lseek() to 'stretch' the file");
    exit(3);
  }
  result = write(fd, "", 1);
  if (result != 1)
  {
    close(fd);
    printf("Error writing last byte of the file");
    exit(4);
  }

  void *map;
  printf("\n About to attach the memory map \n");

  map = mmap(0, final_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
  if (map == MAP_FAILED)
  {
    close(fd);
    printf("Error mmapping the file");
    exit(5);
  }
  bucket = (map + pileup_header_size);
  long li;
  printf("\n About the zero out the file \n");
  for (li = 0; li < tot_recs; li++)
    bucket[li] = 0;
  printf("\n Bucket zero'd.  Reading pileup \n");

  sprintf(sss, "%s", argv[2]);
  if (strcmp(sss, "stdin") == 0)
  {
    use_stdin = TRUE;
  }
  else
  {
    use_stdin = FALSE;
    if ((pileupfile = gzopen(sss, "r")) == (gzFile)NULL)
    {
      printf("\n Can not open file %s for reading which is the mpileup file\n", argv[2]);
      exit(1);
    }
    gzbuffer(pileupfile, 33554432);
  }
  strcpy(sss, argv[3]);
  token = strtok(sss, ".pile");
  char ssp[4196];
  sprintf(ssp, "%s.indel.txt.gz", token);
  if ((indelfile = gzopen(ssp, "wb")) == (gzFile)NULL)
  {
    printf("\n Can not open file %s which is the indel file for writing\n", ssp);
    exit(1);
  }

  printf("\n Working on file %s \n About to read a big bucket\n\n", sss);
  gzbuffer(indelfile, 2000000);
  long which_base;
  int char_map[256];
  int ref_map[256];
  for (i = 0; i < 256; i++)
  {
    char_map[i] = -1;
    ref_map[i] = -1;
  }
  char_map[(int)'a'] = ref_map[(int)'a'] = 0;
  char_map[(int)'c'] = ref_map[(int)'c'] = 1;
  char_map[(int)'g'] = ref_map[(int)'g'] = 2;
  char_map[(int)'t'] = ref_map[(int)'t'] = 3;
  char_map[(int)'A'] = ref_map[(int)'A'] = 0;
  char_map[(int)'C'] = ref_map[(int)'C'] = 1;
  char_map[(int)'G'] = ref_map[(int)'G'] = 2;
  char_map[(int)'T'] = ref_map[(int)'T'] = 3;
  char *line;
  line = cvector(0, 1000000);
  if (use_stdin)
    fgets(line, 999999, stdin);
  else
    gzgets(pileupfile, line, 999999);
  printf("\n Just sucked line %s \n\n", line);
  char last_chr[1000];
  sprintf(last_chr, "!!!!!!");
  token = strtok(line, "\t \n");
  char this_chr[1000];
  strcpy(this_chr, token);
  int chr_no = find_chrom(this_chr, contig_names, no_contigs);
  if (chr_no < 0)
  {
    printf("\n Right off the we don't find the chromosome name %s \n\n", this_chr);
    exit(1);
  }
  token = strtok(NULL, "\t \n");
  int expos = atoi(token) - 1;
  which_base = frag_pos[chr_no - 1] + expos;
  strcpy(last_chr, this_chr);
  gzprintf(indelfile, "Fragment\tPositions\tReference Base\tTotal Coverage\tReference Reads\tNo Deletions\tNo Insertions\tInsertion Sequence\n");
  printf("\nAbout to enter the big loop with %s %d %ld \n\n", this_chr, expos + 1, which_base);
  char *indelseq;
  indelseq = cvector(0, 1000000);

  while (which_base > 0)
  {
    long which_b = (long)which_base * 6;
    if (which_base % 10000000 == 0)
      printf("\n Read  %ld bases \n", which_base);
    token = strtok(NULL, "\t \n");
    char ref_a = token[0];
    char_map[(int)','] = ref_map[(int)ref_a];
    char_map[(int)'.'] = ref_map[(int)ref_a];
    token = strtok(NULL, "\t \n");
    token = strtok(NULL, "\t \n");
    // printf("\n%s %d %c %s\n\n",this_chr,expos,ref_a,token);
    indelseq[0] = '\0';
    int indel_len = 0;
    int no_ins = 0;
    ;
    while (*token)
    {
      if (char_map[(int)*token] >= 0)
        bucket[which_b + (int)char_map[(int)*token]]++;
      else if (*token == '^')
        token++;
      else if (*token == '-')
      {
        // printf("\n About to deletion \n\n");
        token++;
        int this_del = atoi(token);
        // printf("\n\t\t This deletion is length %d \n\n",this_del);
        long start = which_b + 10;
        for (ii = 0; ii < this_del; ii++, start += 6)
        {
          token++;
          bucket[start]++;
        }
      }
      else if (*token == '+')
      {
        // printf("\n About to Insertion \n\n");
        if (no_ins)
          indelseq[indel_len++] = ',';
        no_ins++;
        token++;
        bucket[which_b + 5]++;
        bucket[which_b + (int)char_map[(int)',']]++;
        int this_ins = atoi(token);
        // printf("\n\t\t This insertion is length %d \n\n",this_ins);
        while (isdigit(*token))
          token++;
        for (ii = 0; ii < this_ins; ii++)
          indelseq[indel_len++] = toupper(*token++);
        indelseq[indel_len] = '\0';
      }
      token++;
    }
    if (no_ins > 0)
    {
      int tot_depth, ref_depth, ii;
      tot_depth = 0;
      ref_depth = bucket[which_b + (int)char_map[(int)',']];
      for (ii = 0; ii < 6; ii++)
        tot_depth += bucket[which_b + ii];
      // printf("\n%s\t%d\t%c\t%d\t%d\t%d\t%d\t%s\n\n",this_chr,expos,ref_a,tot_depth,ref_depth,bucket[which_b+4],no_ins,indelseq);

      gzprintf(indelfile, "%s\t%d\t%c\t%d\t%d\t%d\t%d\t%s\n", this_chr, expos + 1, ref_a, tot_depth, ref_depth, bucket[which_b + 4], no_ins, indelseq);
    }

    // printf("\n%ld %ld %d %d %d %d %d %d",which_base,which_b,(int)bucket[which_b],(int)bucket[which_b+1],(int)bucket[which_b+2],
    //	 (int)bucket[which_b+3],(int)bucket[which_b+4],(int)bucket[which_b+5]);
    which_base = 0;
    while ((!gzeof(pileupfile)) && (which_base == 0))
    {
      line[0] = '\0';
      if (use_stdin)
        fgets(line, 999999, stdin);
      else
        gzgets(pileupfile, line, 999999);
      if (strlen(line) < 2)
        which_base = -1;
      else
      {
        token = strtok(line, "\t \n");
        strcpy(this_chr, token);
        if (strcmp(this_chr, last_chr) != 0)
        {
          chr_no = find_chrom(this_chr, contig_names, no_contigs);
          if (chr_no < 0)
            printf("\n Could not find chromosome %s ... skipping \n", this_chr);
          else
          {
            token = strtok(NULL, "\t \n");
            expos = atoi(token) - 1;
            which_base = frag_pos[chr_no - 1] + expos;
          }
          strcpy(last_chr, this_chr);
        }
        else
        {
          if (chr_no >= 0)
          {
            token = strtok(NULL, "\t \n");
            expos = atoi(token) - 1;
            which_base = frag_pos[chr_no - 1] + expos;
            if (which_base > frag_pos[chr_no])
              which_base = 0;
          }
        }
      }
    }
  }
  if (!use_stdin)
    gzclose(pileupfile);
  gzclose(indelfile);
  close(fd);

  return 0;
}
/*---------------------------------------------------------------------*/
int find_chrom(char *this, char **contigs, int n)
{
  int i;
  for (i = 0; i < n; i++)
    if (strcmp(this, contigs[i]) == 0)
      return i;
  return -1;
}

/*---------------------------------------------------------------------*/
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
/*---------------------------------------------------------------------*/
unsigned int **umatrix(int nrl, int nrh, int ncl, int nch)
{
  int i;
  unsigned int **m;

  m = (unsigned int **)malloc((unsigned)(nrh - nrl + 1) * sizeof(unsigned int *));
  if (!m)
    dump_error("allocation failure 1 in ulmatrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++)
  {
    m[i] = (unsigned int *)malloc((unsigned)(nch - ncl + 1) * sizeof(unsigned int));
    if (!m[i])
      dump_error("allocation failure 2 in ulmatrix()");
    m[i] -= ncl;
  }
  return m;
}

/*---------------------------------------------------------------------*/
unsigned int *uvector(int nl, int nh)
{
  unsigned int *v;

  v = (unsigned int *)malloc((unsigned)(nh - nl + 1) * sizeof(int));
  if (!v)
    dump_error("allocation failure in ulvector()");
  return v - nl;
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
/*---------------------------------------------------------------------*/

void dump_error(char *error_text)
{

  fprintf(outfile, "PEmapper error...\n");
  fprintf(outfile, "%s\n", error_text);
  fprintf(outfile, "...now exiting to system...\n");
  exit(1);
}

/*---------------------------------------------------------------------*/

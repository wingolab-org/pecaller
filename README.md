# PEMapper/PECaller

**Citation**:

Johnston HR, Chopra P, Wingo TS, Patel V, Epstein MP, Mulle JG, Warren ST,
Zwick ME, Cutler DJ. PEMapper and PECaller provide a simplified approach to
whole-genome sequencing. **Proc Natl Acad Sci U S A** **2017**:114(10):E1923-E1932
DOI: [10.1073/pnas.1618065114](https://www.doi.org/10.1073/pnas.1618065114)
PMID: [28223510](https://pubmed.ncbi.nlm.nih.gov/28223510/)

## Program Versions

This project includes two generations of tools:

**Version 1 (Original)**: `pecaller`, `pemapper`, `snp_to_vcf`

- Original implementation (2014)
- Uses text-based pileup format
- Kept for informational and compatibility purposes.

**Version 2 (Current/Recommended)**: `pecaller2`, `pemapper2`, `snp_to_vcf2`

- Current implementation (2017-2018)
- Uses optimized binary pileup format for improved performance
- Includes `mpileup_to_pileup2` converter for format transformation
- Both versions use the same genome index created by `index_genome_whole`
- **This is the primary version described in this documentation**

**All new projects should use Version 2 tools**.

## Installation

Install using `make` and further instructions are in [Build.md](./BUILD.md).
During the build process, you will see warnings about unused parameters, which is fine.

```shell
# Standard system-wide install
sudo make install

# User-specific install
make install PREFIX=~/.local

# Custom location
make install PREFIX=/opt/mytools
```

## Workflow Overview

The PECaller2 pipeline processes sequencing data through the following stages:

### 1. **Genome Indexing**

Create a searchable index of your reference genome:

```
FASTA files → index_genome_whole → reference.sdx + reference.seq + reference.idx
```

### 2. **Read Mapping & Format Conversion**

Convert sequencing data to binary pileup format:

**Option A - Direct mapping from FASTQ:**

```
FASTQ files → pemapper2 → sample.indel.txt + sample.pileup2
```

**Option B - Convert from existing alignments:**

```
CRAM/BAM → samtools mpileup → mpileup_to_pileup2 → sample.indel.txt + sample.pileup2
```

_Note: For large studies, group samples into batches (B1, B2, etc.) of ~100 samples each._

### 3. **Variant Calling**

Call variants from pileup data:

```
*.indel.txt + *.pileup2 + reference.sdx + regions.bed → pecaller2 → *.base + *.snp + *.dist + *.piles
```

### 4. **Quality Control & Site Selection**

Identify high-quality variant sites:

```
*.snp + reference.sdx → make_snplist_formerge.pl → good.bed + bad.bed
```

Then split by chromosome: `good.bed → chr[1-22,X,Y,M].good.bed`

### 5. **Cross-Batch Merging**

Combine variant calls across all batches:

```
*.base + chr?.good.bed + reference.sdx → pecall_merger → merged.chr?.snp
```

### 6. **Final Integration**

Integrate indels with SNP calls:

```
*.indel.txt + merged.chr?.snp + reference.sdx → merge_indel_snp.pl → final.chr?.snp
```

## Build genome index

### Choosing an Indexing Tool

This project provides two genome indexing programs:

**`index_genome_whole`** (Recommended):

- Designed for large genomes (human, mouse, etc.)
- Fixed 16-mer indexing depth optimized for performance
- Uses compressed storage (gzipped files) to save disk space
- No genome size limitations
- Uses ~108 Gb to encode human genome.

**`index_genome`** (For development and experimental use):

- Configurable indexing depth and genome size limits
- Memory-efficient for smaller genomes
- Uncompressed output files
- Intended for experimental use

### Building the Genomic Index

- Download the genome fa files.

For example, to download all of mm10 to the current directory from UCSC:

```bash
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/ ./
```

- Merge all of the fa files into one large fa file.

Use `merge_dir_fa.pl` to merge the files into a single fa file with chromosomes
sorted in a natural order (e.g., 1-19, M, X, and Y) with unmapped chromosome
after chrY.

```bash
merge_dir_fa.pl -d ../mm10_2015-01-25/ -c '1-19,M,X,Y' -o mm10
```

- Index the merged fa file.

**Using `index_genome_whole` (recommended):**
Run interactively and provide:

- Output destination (Screen or Disk)
- Input FASTA file path
- Output basename for index files
- Whether DNA will be bisulfite converted (y/n)

**Using `index_genome` (for smaller genomes):**
Run interactively or redirect input from a file.
For example, `index_genome < in.cmds` where `in.cmds` contains:

```text
    d
    1000
    mm10.fa
    mm10
    n
```

## Read Mapping

**Using `pemapper2` (recommended):**

- Single-ended or paired-ended mapping of FASTQ files is supported
- Outputs binary pileup2 format for optimal performance with `pecaller2`
- Automatically generates indel.txt files for insertion/deletion tracking

**Legacy `pemapper` support:**

- Available for compatibility with existing workflows
- Outputs standard pileup format for use with original `pecaller`

**File naming conventions:**

- Paired-end files should use `_1_` vs `_2_` or `_R1_` vs `_R2_` naming

## Variant Calling

**Using `pecaller2` (recommended):**

- Place all pileup2 files into a single directory and launch `pecaller2`
- Supports BED file regions to restrict calling to specific genomic intervals
- **Important:** Chromosome order in BED file must match exactly the order in the genome `.sdx` file
- Outputs: `.base`, `.snp`, `.dist`, and `.piles` files for comprehensive variant analysis

**Legacy `pecaller` support:**

- Available for compatibility with text-based pileup inputs
- Use for workflows requiring the original format

## Cross-Sample Merging

**Step 1: Quality Control & Site Selection**

- Run `make_snplist_formerge.pl` on SNP files to identify high-quality variant sites
- This creates "good.bed" and "bad.bed" files for quality-based filtering
- Split good.bed by chromosome for parallel processing: `chr[1-22,X,Y,M].good.bed`

**Step 2: Merge Base Files**

- Use `pecall_merger` to combine `.base` files across samples at high-quality sites
- Input: `.base` files + chromosome-specific `.good.bed` files + genome `.sdx` file
- Output: Merged SNP files per chromosome

## Final Variant Integration

**Combining SNPs and Indels:**

- Both `pecaller2` and `pecaller` call each base independently
- Insertions are stored separately in `.indel.txt` files during mapping/calling
- Use `merge_indel_snp.pl` to integrate indels into final SNP files
- This step sorts variants and places all indels into the final output

## VCF Output

For compatibility with standard genomics tools:

**Using `snp_to_vcf2` (recommended):**

- Converts the optimized SNP format to standard VCF format
- Supports both SNPs and indels in VCF output
- Improved handling of complex variants and metadata

**Legacy `snp_to_vcf` support:**

- Available for compatibility with Version 1 workflows
- Converts original text-based SNP format to VCF

## Quality Control

- `snp_tran_counter.pl` and `snp_tran_silent_rep.pl` give transition to
  transversion counts for variants obsered. The latter provides counts specific
  to different kinds of classes of variants (i.e., replacement sites) and expects
  the annotation to come from the [Bystro](https://github.com/bystrogenomics/bystro)
  genomic annotator.

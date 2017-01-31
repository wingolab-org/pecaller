# PEMapper/PECaller

Current Citation (Pre-submission):
[Johnston HR, Chopra P, Wingo T, et al. PEMapper / PECaller: A simplified approach to whole-genome sequencing. bioRxiv. 2016](http://biorxiv.org/content/early/2016/09/22/076968)

## Build genome index

* Download the genome fa files.

For example, to download all of mm10 to the current directory from UCSC:

		rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/ ./

* Merge all of the fa files into one large fa file.

Use `merge_dir_fa.pl` to merge the files into a single fa file with chromosomes
sorted in a natural order (e.g.,  1-19, M, X, and Y) with unmapped chromosome
after chrY.

		./merge_dir_fa.pl -d ../mm10_2015-01-25/ -c '1-19,M,X,Y' -o mm10

* Index the merged fa file.

Use `index_genome` either interactively redirect input from a file.
For example, `index_genome < in.cmds` where `in.cmds` contains:

    d
    1000
    mm10.fa
    mm10
    n

## Map `Fastq` files

- Single-ended or paired-ended mapping of fastq files is supported.
- Mapping a collection of files is also supported and made easier using
`map_directory_array.pl`, which does make some assumptions about the naming of
the files and that paired-ended files either have `_1_` vs `_2_` or `_R1_` vs
`_R2_`.

## Basecalling

- place all pileup files into a single directory and call using `pecaller`
launched from that directory. Note: if you experience a segmentation fault
just after running the command it is likely that you have not supplied the
command correctly. An `argtable3` or `getopt` interface would be a nice
addition.
- `pecaller` will make a base and snp file for all sites in the user supplied
region (or every site covered in the pileup file if one is not provided).

## Merging basecall files together

- You may have multiple base and snp files from different basecall runs that
would be nice to be merged together. To do this, place all of the called files
(i.e., snp, base, and indel files) into the same directory (or symlink them) and
use `make_snplist_formerge.pl` to create a "good" and "bad" list of sites, based
on variants sites in the snpfiles. Then provide `pecall_merger` with the base
files and sites that should be merged.

## Adding Indels to the snpfiles

- After `pecaller` writes the snpfile it will have each indel site listed but
indels might not be listed in a useful way since all sites are called
independently. To merge all contiguous deletions together and insert the most
common insertion sequence use `merge_indel_snp.pl`.

## Q/C

- `snp_tran_counter.pl` and `snp_tran_silent_rep.pl` give Transition to
Transversion counts for different sites and kinds of changes within the genome.
The `snp_tran_silent_rep.pl` expects the annotation to come from
[SeqAnt](https://seqant.emory.edu/).

## Contributing / Folding in Dave's changes

1. Change EOL characters to unix `set ff=unix` in vim, for example.
2. indent things in a consistent way: `indent -bli0 -l120`.

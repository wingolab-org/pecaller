# pecaller

This software maps and calls next generation genomic sequencing. It was written by David J. Cutler, Ph.D.

## Folding in Dave's code

1. Change EOL characters to unix `set ff=unix` in vim.
2. indent things in a consistent way: `indent -bli0 -l120`.

## Build genome index

* Download the genome fa files. 

For instance, to download all of mm10 into the current working directory from UCSC try the following:

	rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/ ./

* Merge all of the fa files into one large fa file.

Use the Perl script `merge_dir_fa.pl` to merge the files. The intention of the script is to order the index so that the named chromosomes (here, chrs 1-19, M, X, and Y) are in the natural order and after those are placed any unmapped or alternative chromosome.

	perl merge_dir_fa.pl -d ../mm10_2015-01-25/ -c '1-19,M,X,Y' -o mm10

* Index the merged fa file.

Use `index_genome` either interactively to index the merged fa file. You can redirect the expected input from a file. For example, `index_genome < in.cmds` where `in.cmds` is a file with the following information:

    d
    1000
    mm10.fa
    mm10
    n

## Map `Fastq` files

You can either map a single fastq file, paired-ended fastq files, or collection of paired-ended fastq files. 


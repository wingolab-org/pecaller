# Build Instructions

This project uses a Makefile to build C programs and install Perl scripts.

During the build process, you will see warnings about unused parameters, which is fine.

## Quick Start

```bash
# Build all programs
make

# Build with development flags (extra warnings)
make dev

# Install to /usr/local/bin (requires sudo)
sudo make install

# Install to custom location
make install PREFIX=/opt/local

# Clean build artifacts
make clean
```

## Available Targets

- `make` or `make all` - Build all C programs
- `make dev` - Development build with extra compiler warnings
- `make clean` - Remove all build artifacts
- `make install` - Install binaries and Perl scripts to `$(PREFIX)/bin`
- `make uninstall` - Remove installed files
- `make help` - Show available targets and options
- `make debug` - Show detected source files and build targets

## Programs Built

The following C programs are built from the `c/` directory:

- `index_genome_whole` - Genome indexing utility
- `index_genome` - Genome indexing utility
- `mpileup_to_pileup2` - Pileup format converter
- `pecall_merger` - PE call merger utility
- `pecaller` - PE variant caller
- `pecaller2` - PE variant caller (version 2)
- `pemapper` - PE read mapper
- `pemapper2` - PE read mapper (version 2)
- `snp_to_vcf` - SNP to VCF converter
- `snp_to_vcf2` - SNP to VCF converter (version 2)

## Perl Scripts Installed

The following Perl scripts from the `perl/` directory are also installed:

- `make_snplist_formerge.pl`
- `merge_indel_compressed_snp.pl`
- `merge_indel_snp.pl`
- `snp_filter_call_rate.pl`
- `snp_tran_silent_rep_site.pl`

## Build Directory

All compiled binaries are placed in the `build/` directory to keep the source tree clean.

## Dependencies

- GCC compiler
- zlib development libraries (`libz-dev` or `zlib-devel`)
- Math library (usually included with GCC)
- pthread library (usually included with GCC)

## Installation

By default, programs are installed to `/usr/local/bin`. You can change this with:

```bash
make install PREFIX=/your/custom/path
```

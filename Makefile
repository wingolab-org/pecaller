# Makefile for PE Mapper / Caller Tools

CC         = gcc
CC_options = -Wall -Wextra -g -v -O3 -std=gnu11
INCLUDES   =
CFLAGS     = $(CC_OPTIONS) $(INCLUDES)
LIBS       = -lm -lz -lpthread

DUMP_PILEUP_EXE = ./bin/dump_pileups
DUMP_PILEUP_SRC = dump_pileups.c

PEMAPPER_TSW_EXE = ./bin/pemapper_tsw
PEMAPPER_TSW_SRC = pemapper_tsw.c

PECALLER_EXE = ./bin/pecaller
PECALLER_SRC = pecaller.c

PEMAPPER_EXE = ./bin/pemapper
PEMAPPER_SRC = pemapper.c

PROGS = $(DUMP_PILEUP_EXE) $(PEMAPPER_TSW_EXE) $(PECALLER_EXE) $(PEMAPPER_EXE)

$(DUMP_PILEUP_EXE):
	$(CC) $(DUMP_PILEUP_SRC) -o $@ $(CFLAGS) $(LIBS)

$(PEMAPPER_TSW_EXE): 
	$(CC) $(CFLAGS) $(LIBS) $(PEMAPPER_TSW_SRC) -o $@

$(PECALLER_EXE):
	$(CC) $(CFLAGS) $(LIBS) $(PECALLER_SRC) -o $@

$(PEMAPPER_EXE):
	$(CC) $(CFLAGS) $(LIBS) $(PEMAPPER_SRC) -o $@

all: introduce $(PROGS)
	@echo done.

introduce: 
	@echo "building: dump_pileup pecaller pemapper pemapper_tsw"

clean:
	rm -f *.o

distclean:
	rm -f $(PROGS) *.o

## end of Makefile
# DO NOT DELETE THIS LINE -- make depend depends on it.

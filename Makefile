# Makefile for PE Mapper / Caller Tools

CC         = gcc 
CC_OPTIONS = -Wall -Wextra
INCLUDES   =
CFLAGS     = $(CC_OPTIONS) $(INCLUDES)
LIBS       = -lm -lz -lpthread

SRC=$(wildcard src/*.c)
EXE=$(patsubst %.c,%,$(SRC))

all: build $(EXE) slink_pl

$(EXE):
	$(CC) $(CFLAGS) $(LIBS) -o $@ $@.c
	mv $@ bin/

clean:
	rm -rf bin
  
build:
	@mkdir -p bin 

slink_pl:
	cp ./src/*.pl bin/
	chmod 755 bin/*

## end of Makefile
# DO NOT DELETE THIS LINE -- make depend depends on it.


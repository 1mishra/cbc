# Makefile for building C programs

SRC=main.c sam_models.c sam_file_allocation.c qv_codebook.c qv_compression.c read_compression.c read_decompression.c util.c well.c Arithmetic_stream.c stream_model.c pmf.c distortion.c quantizer.c sfio_ssh_session.c compression.c io_functions.c id_compression.c aux_compression.c

OBJ=$(SRC:.c=.o)

CC=gcc
RM=rm -f

CFLAGS=-ffast-math -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -Wall -I../include -I/usr/local/include  -D__APPLE__
LDFLAGS=-L/usr/local/lib -lc -lm -lssh

%.o : %.c
	$(CC) -gdwarf-2 -g $(CFLAGS) -c $< #he add -gdwarf-2 -g para poder debugear bien en XCode http://stackoverflow.com/questions/21132194/is-it-possible-to-debug-a-gcc-compiled-program-using-lldb-or-debug-a-clang-comp

all : program

program : $(OBJ)
	$(CC) $(OBJ) -o program $(LDFLAGS)

debug : CFLAGS += -DDEBUG -ggdb
debug : program

clean :
	$(RM) *.o program
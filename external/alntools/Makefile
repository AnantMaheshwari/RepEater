# makefile for alntools

CFLAGS = -O3
#CFLAGS = -g	# for debugging

LIBS = -lpthread -lz

ALL = tanbed gdbmask svfind taco ONEview

DESTDIR = ~/bin

all: $(ALL)

install:
	cp $(ALL) $(DESTDIR)

clean:
	$(RM) *.o *~ $(ALL)
	$(RM) -r *.dSYM

### object files

UTILS_OBJS = utils.o array.o dict.o
UTILS_HEADERS = utils.h array.h dict.h
$(UTILS_OBJS): $(UTILS_HEADERS)

ONElib.o: ONElib.h 

tanbed.o: alntools.h ONElib.h $(UTILS_HEADERS)

gdb.o: alntools.h ONElib.h $(UTILS_HEADERS)

SEQIO_OPTS = -DONEIO
seqio.o: seqio.c seqio.h ONElib.h $(UTILS_HEADERS)
	$(CC) $(CFLAGS) $(SEQIO_OPTS) -c $^

alnseq.o: alnseq.h ONElib.h

alncode.o: alncode.h align.h

### programs

tanbed: tanbed.o gdb.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

gdbmask: gdb.c ONElib.o $(UTILS_OBJS)
	$(CC) -D GDB_MASK $(CFLAGS) -o $@ $^ $(LIBS)

salsa.o: salsa.h ONElib.h $(UTILS_HEADERS)

taco: taco.c salsa.o gdb.o seqio.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

svfind: svfind.c alnseq.o alncode.o seqio.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

ONEview: ONEview.c ONElib.o
	$(CC) $(CFLAGS) -o $@ $^

### test

### end of file

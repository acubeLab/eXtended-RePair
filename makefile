# Compilation flags
CC=gcc 
CFLAGS= -std=c11 -Wall -g -m64 -O2

# executables in this directory
EXECS = xrepair.x despair.x
OBJS = array.o hash.o heap.o records.o basics.o
HEADERS = hash.h heap.h records.h basics.h

# comment out this definition to get rid of malloc_count 
# MALLOC_FLAGS=tools/malloc_count.c -DMALLOC_COUNT -ldl

# malloc_count dependencies
ifdef MALLOC_FLAGS
MALLOC_FILES=tools/malloc_count.c tools/malloc_count.h
endif


all: $(EXECS)

# targets not producing a file declared phony
.PHONY: all clean


xrepair.x: xrepair.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

despair.x: despair.o basics.o 
	$(CC) $(CFLAGS) -o $@ $^

xrepair.o: xrepair.c $(HEADERS)
	$(CC) $(CFLAGS) -c xrepair.c

despair.o: despair.c basics.h
	$(CC) $(CFLAGS) -c despair.c

irepair: irepair.o array.o hash.o heap.o records.o basics.o makefile
	$(CC) $(CFLAGS) -o irepair irepair.o array.o hash.o heap.o records.o basics.o

irepair.o: irepair.c array.h hash.h heap.h records.h basics.h makefile
	$(CC) $(CFLAGS) -c irepair.c

idespair: idespair.o basics.o makefile
	$(CC) $(CFLAGS) -o idespair idespair.o basics.o

idespair.o: idespair.c basics.h makefile
	$(CC) $(CFLAGS) -c idespair.c

%.o: %.c $(HEADERS)



#~ array.o: array.c array.h hash.h heap.h records.h basics.h makefile
#~ 	$(CC) $(CFLAGS) -c array.c

#~ hash.o: hash.c array.h hash.h heap.h records.h basics.h makefile
#~ 	$(CC) $(CFLAGS) -c hash.c

#~ heap.o: heap.c array.h hash.h heap.h records.h basics.h makefile
#~ 	$(CC) $(CFLAGS) -c heap.c

#~ records.o: records.c array.h hash.h heap.h records.h basics.h makefile
#~ 	$(CC) $(CFLAGS) -c records.c

#~ basics.o: basics.c basics.h makefile
#~ 	$(CC) $(CFLAGS) -c basics.c


test: xrepair.x despair.x
	xrepair.x README 
	despair.x README
	cmp README README.out

clean:
	rm -f $(EXECS) *.o


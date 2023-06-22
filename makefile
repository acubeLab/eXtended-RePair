EXECS = xrepair despair
CFLAGS= -std=c11 -Wall -g -m64 -O2

all: $(EXECS)

# targets not producing a file declared phony
.PHONY: all clean


xrepair: xrepair.o array.o hash.o heap.o records.o basics.o makefile
	gcc $(CFLAGS) -o xrepair xrepair.o array.o hash.o heap.o records.o basics.o

xrepair.o: xrepair.c array.h hash.h heap.h records.h basics.h makefile
	gcc $(CFLAGS) -c xrepair.c

despair: despair.o basics.o makefile
	gcc $(CFLAGS) -o despair despair.o basics.o

despair.o: despair.c basics.h makefile
	gcc $(CFLAGS) -c despair.c

irepair: irepair.o array.o hash.o heap.o records.o basics.o makefile
	gcc $(CFLAGS) -o irepair irepair.o array.o hash.o heap.o records.o basics.o

irepair.o: irepair.c array.h hash.h heap.h records.h basics.h makefile
	gcc $(CFLAGS) -c irepair.c

idespair: idespair.o basics.o makefile
	gcc $(CFLAGS) -o idespair idespair.o basics.o

idespair.o: idespair.c basics.h makefile
	gcc $(CFLAGS) -c idespair.c

array.o: array.c array.h hash.h heap.h records.h basics.h makefile
	gcc $(CFLAGS) -c array.c

hash.o: hash.c array.h hash.h heap.h records.h basics.h makefile
	gcc $(CFLAGS) -c hash.c

heap.o: heap.c array.h hash.h heap.h records.h basics.h makefile
	gcc $(CFLAGS) -c heap.c

records.o: records.c array.h hash.h heap.h records.h basics.h makefile
	gcc $(CFLAGS) -c records.c

basics.o: basics.c basics.h makefile
	gcc $(CFLAGS) -c basics.c

toint: toint.c makefile
	gcc $(CFLAGS) -o toint toint.c


test: xrepair despair
	xrepair README 
	despair README
	cmp README README.out

clean:
	rm -f $(EXECS) *.o


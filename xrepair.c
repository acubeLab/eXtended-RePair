/* repairX.c
   Giovanni Manzini    June 20st, 2023

  Implementation of repair algorithm for characters with bounded memory
  by Gonzalo Navarro with the following additional options:
    -x C    [do not generate rules involving chars <= C]
    -r R    [do not generate more than R rules]
    -s      [try to save space computing the first 2^16-2^8 rules
             using a quadratic algorithm]
  The input file name is a required input parameter. 
  The program produces two files:
    - a file with extension .R containing the rules
    - a file with extension .C containing the compressed text
  The option -s only affects the running time, not the output.


  This code is derived by the large and balanced repair implementation,
  which used the following approach: since the working memory is dominated
  by the term tlen * 3 * sizeof(long long) [=24*tlen] when the option -s
  is specified, we compute the first 2^16-2^8 rules storing the input in an
  array of shorts and using a quadratic algorithm (one a rule is found is 
  applied with a complete scan on the input). This will hopefully significantly 
  reduce the length of the working sequence so that the complete machinery
  for the linear time algoritm fits in the available memory.

*/

/*

Original notice by Gonzalo Navarro:

Repair -- an implementation of Larsson and Moffat's compression and
decompression algorithms.
Copyright (C) 2010-current_year Gonzalo Navarro

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

Author's contact: Gonzalo Navarro, Dept. of Computer Science, University of
Chile. Blanco Encalada 2120, Santiago, Chile. gnavarro@dcc.uchile.cl

*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdint.h>
#include <assert.h>
#include <errno.h>

#include "basics.h"
#include "records.h"
#include "hash.h"
#include "heap.h"

// debug control
int PRNC = 0;  // print current sequence C (verbose!)
int PRNR = 0;  // print active pairs in the heap (verbose!)
int PRNCf = 0; // print final sequence C
int PRNP =  0;  // print forming pairs
int PRNL =  0;  // print progress on text scan
void prnSym(int c);
void prnsC(unsigned short *, relong len);
void prnC(void);
void prnRec(void);
static void quit(const char *s);
static void *mymalloc(size_t size, int line, const char *file);
static void *myrealloc(void *ptr, size_t size, int line, const char *file);
static void usage_and_exit(char *name);


float factor = 0.75; // 1/extra space overhead; set closer to 1 for smaller
                     // and slower execution
int minsize = 256;   // to avoid many reallocs at small sizes, should be ok as is

relong *C; // compressed text, pero tiene -ptrs al mismo texto

relong u; // |text| and later current |C| with gaps
relong c; // real |C| (without gaps)

int alph; // alphabet size and first non terminal symbol (is written to .R file)
int n;    // alph + |R| = next id for a non terminal symbol

Tlist *L; // |L| = c; // next and previous entry in gapped |C|

Thash Hash; // hash table of pairs
Theap Heap; // special heap of pairs
Trarray Rec; // records


// globals controlled by the input parameters
static int32_t maxForbiddenChar; // do not generate rules involving chars <= this
// default is -1, i.e. no restriction
static int64_t maxRules;        // max number of rules to generate
// default -1, i.e. no restriction
static int Verbose;


// return true if this pair of symbols should never appear
// as the left-hand side of a rule
int forbidden_pair(relong left, relong right)
{
  assert(left>=0 && right >=0);
  if (maxForbiddenChar < 0)
    return 0;
  if (left <= maxForbiddenChar)
    return 1;
  if (right <= maxForbiddenChar)
    return 1;
  return 0;
}



// init compute alphabet size and init Rec, Hash, and Heap
void init(uint8_t *text, relong len, FILE *R)
{
  relong i;
  int id;
  Tpair pair;
  // compute largest input symbol
  alph = 0;
  for (i = 0; i < len; i++) {
    if (text[i] > alph)
      alph = text[i];
  }
  // init alphabet size and n as first code usable for non terminal and save it to file
  n = ++alph;
  if (fwrite(&alph, sizeof(int), 1, R) != 1)
    quit("Error writing alphabet size to .R file");
  // init hash and heap of pairs 
  Rec = createRecords(factor, minsize);
  Heap = createHeap(len, &Rec, factor, minsize);
  Hash = createHash(256 * 256, &Rec);
  assocRecords(&Rec, &Hash, &Heap, NULL);
  // init all valid pairs
  for (i = 0; i < len - 1; i++) {
    pair.left = text[i];
    pair.right = text[i + 1];
    if(!forbidden_pair(pair.left,pair.right)) {
      id = searchHash(Hash, pair);
      if (id == -1) { // new pair, insert
        id = insertRecord(&Rec, pair);
      }
      else {
        incFreq(&Heap, id);
      }
      if (PRNL && (i % 1000000 == 0))
        printf("Processed %lli chars\n", i);
    }
  }
  purgeHeap(&Heap); // remove pairs with freq==1
}



// third pass, use the full machinery
void prepare(uint8_t *text, uint16_t *text16, relong len)
{
  relong i;
  int id;
  Tpair pair;
  c = u = len;
  C = mymalloc(u * sizeof(relong),__LINE__,__FILE__);
  if (text16 != NULL) {
    assert(text==NULL);
    for (i = 0; i < u; i++)
      C[i] = text16[i];  // copy from short array (repair1)
    free(text16);
  }
  else {
    assert(text16==NULL);
    for (i = 0; i < u; i++)
      C[i] = text[i]; // copy from text array (repair0)
    free(text);
  }
  // init prev/next list L
  L = mymalloc(u * sizeof(Tlist),__LINE__,__FILE__);
  assocRecords(&Rec, &Hash, &Heap, L);
  for (i = 0; i < c - 1; i++) {
    pair.left = C[i];
    pair.right = C[i + 1];
    if(forbidden_pair(pair.left,pair.right)) {
      L[i].prev = NullFreq;
      L[i].next = -1;
      continue;
    }
    id = searchHash(Hash, pair);
    if (id == -1) { // non existing => occurs only once, don't create
      L[i].prev = NullFreq;
      L[i].next = -1;
    }
    else {
      if (Rec.records[id].cpos == -1) { // first time I see it this pass
        L[i].next = -1;
      }
      else {
        L[i].next = Rec.records[id].cpos;
        L[L[i].next].prev = i;
      }
      L[i].prev = -id - 1;
      Rec.records[id].cpos = i;
    }
    if (PRNL && (i % 1000000 == 0))
      printf("Processed %lli chars\n", i);
  }
  L[i].prev = NullFreq;
  L[i].next = -1;
  purgeHeap(&Heap); // probably not necessary since we did not modify the heap
}



// second pass, work on a version of shorts until shorts are insufficient
// TODO: until the structure fits in MB, or 
relong repair1(uint16_t *sC, relong len, FILE *R)
{
  int oid, id;
  relong cpos, pos;
  Trecord *orec;
  Tpair pair;
  int left, right;
  if (PRNL)
    printf("--- second stage, n=%lli\n", len);
  if (PRNC)
    prnsC(sC,len);
  while (n < 1 << 16 && (n-alph) < maxRules) {
    // if ((len / 1024 / 1024) * 3 * sizeof(relong) <= MB)
    //  return len;
    if (PRNR)
      prnRec();
    oid = extractMax(&Heap);
    if (oid == -1)
      break; // the end!!
    orec = &Rec.records[oid];
    if (fwrite(&orec->pair, sizeof(Tpair), 1, R) != 1)
      quit("Error writing to .R file in repair1()");
    left = orec->pair.left;
    right = orec->pair.right;
    // ofreq = orec->freq;
    if (PRNP) {
      printf("Chosen pair %i = (", n);
      prnSym(orec->pair.left);
      printf(",");
      prnSym(orec->pair.right);
      printf(") (%lli occs)\n", orec->freq);
    }
    pos = 0;
    for (cpos = 0; cpos < len - 1; cpos++) {
      if ((sC[cpos] != left) || (sC[cpos + 1] != right))
        sC[pos] = sC[cpos];
      else { // occurrence of the pair to modify
        // decr freqs of pairs that disappear
        if (pos > 0) {
          pair.left = sC[pos - 1];
          pair.right = sC[cpos];
          id = searchHash(Hash, pair);
          if (id != -1) { // may not exist if purgeHeap'd
            if (id != oid)
              decFreq(&Heap, id); // not to my pair!
          }
        }
        if (cpos < len - 2) {
          pair.left = sC[cpos + 1];
          pair.right = sC[cpos + 2];
          id = searchHash(Hash, pair);
          if (id != -1) { // may not exist if purgeHeap'd
            if (id != oid)
              decFreq(&Heap, id); // not to my pair!
          }
        }
        // create or incr freqs of pairs that appear
        if (pos > 0) {
          pair.left = sC[pos - 1];
          pair.right = n;
          id = searchHash(Hash, pair);
          if (id == -1) { // new pair, insert
            id = insertRecord(&Rec, pair);
          }
          else {
            incFreq(&Heap, id);
          }
        }
        if (cpos < len - 2) {
          pair.left = n;
          pair.right = sC[cpos + 2];
          id = searchHash(Hash, pair);
          if (id == -1) { // new pair, insert
            id = insertRecord(&Rec, pair);
          }
          else {
            incFreq(&Heap, id);
          }
        }
        // actually do the substitution
        sC[pos] = n;
        cpos++;
      }
      pos++;
    }
    if (cpos == len - 1)
      sC[pos++] = sC[cpos];
    len = pos;
    if (PRNC)
      prnsC(sC,len);
    removeRecord(&Rec, oid);
    n++;
    purgeHeap(&Heap);                      // remove freq 1 from heap
  }
  purgeHeap(&Heap); // remove freq 1 from heap, if it exited for oid=-1
  return len;
}

// to avoid another list to access the sparse C we thread it using the
// empty space. if next cell of an active cell is negative, it is
// (minus) a ptr to the next occ. idem previous cell to previous occ,
// except that next ptr dominates over prev ptr if they must be in
// the same cell. but in this case one can find prev in O(1) anyway.
relong repair(FILE *R)
{
  int oid, id;
  relong cpos;
  Trecord *rec, *orec;
  Tpair pair;
  if (PRNL)
    printf("--- third stage, n=%lli\n", c);
  if (PRNC)
    prnC();
  // printf("repair: rules:%lld maxRules=%lld\n",n-alph,maxRules);
  while (n + 1 > 0 && (n-alph) < maxRules) {
    if (PRNR)
      prnRec();
    oid = extractMax(&Heap);
    if (oid == -1)
      break; // the end!!
    orec = &Rec.records[oid];
    cpos = orec->cpos;
    if (fwrite(&orec->pair, sizeof(Tpair), 1, R) != 1)
      quit("Error writing to .R file in repair()");
    if (PRNP) {
      printf("Chosen pair %i = (", n);
      prnSym(orec->pair.left);
      printf(",");
      prnSym(orec->pair.right);
      printf(") (%lli occs)\n", orec->freq);
    }
    while (cpos != -1) {
      relong ant, sgte, ssgte;
      // replacing bc->e in abcd, b = cpos, c = sgte, d = ssgte
      if (C[cpos + 1] < 0)
        sgte = -C[cpos + 1] - 1;
      else
        sgte = cpos + 1;
      if ((sgte + 1 < u) && (C[sgte + 1] < 0))
        ssgte = -C[sgte + 1] - 1;
      else
        ssgte = sgte + 1;
      // remove bc from L
      if (L[cpos].next != -1)
        L[L[cpos].next].prev = -oid - 1;
      orec->cpos = L[cpos].next;
      if (ssgte != u) { // there is ssgte
        // remove occ of cd
        pair.left = C[sgte];
        pair.right = C[ssgte];
        id = searchHash(Hash, pair);
        if (id != -1) { // may not exist if purgeHeap'd
          if (id != oid)
            decFreq(&Heap, id);         // not to my pair!
          if (L[sgte].prev != NullFreq) { // still exists(not removed)
            rec = &Rec.records[id];
            if (L[sgte].prev < 0) // this cd is head of its list
              rec->cpos = L[sgte].next;
            else
              L[L[sgte].prev].next = L[sgte].next;
            if (L[sgte].next != -1) // not tail of its list
              L[L[sgte].next].prev = L[sgte].prev;
          }
        }
        // create occ of ed
        pair.left = n;
        if(forbidden_pair(pair.left,pair.right)) {
          L[cpos].prev = NullFreq;
          L[cpos].next = -1;
        }
        else {
          id = searchHash(Hash, pair);
          if (id == -1) { // new pair, insert
            id = insertRecord(&Rec, pair);
            rec = &Rec.records[id];
            L[cpos].next = -1;
          }
          else {
            incFreq(&Heap, id);
            rec = &Rec.records[id];
            L[cpos].next = rec->cpos;
            L[L[cpos].next].prev = cpos;
          }
          L[cpos].prev = -id - 1;
          rec->cpos = cpos;
        }
      }
      if (cpos != 0) { // there is ant
        // remove occ of ab
        if (C[cpos - 1] < 0) {
          ant = -C[cpos - 1] - 1;
          if (ant == cpos) // sgte and ant clashed -> 1 hole
            ant = cpos - 2;
        }
        else
          ant = cpos - 1;
        pair.left = C[ant];
        pair.right = C[cpos];
        id = searchHash(Hash, pair);
        if (id != -1) { // may not exist if purgeHeap'd
          if (id != oid)
            decFreq(&Heap, id);        // not to my pair!
          if (L[ant].prev != NullFreq) { // still exists (not removed)
            rec = &Rec.records[id];
            if (L[ant].prev < 0) // this ab is head of its list
              rec->cpos = L[ant].next;
            else
              L[L[ant].prev].next = L[ant].next;
            if (L[ant].next != -1) // it is not tail of its list
              L[L[ant].next].prev = L[ant].prev;
          }
        }
        // create occ of ae
        pair.right = n;
        if(forbidden_pair(pair.left,pair.right)) {
          L[ant].prev = NullFreq;
          L[ant].next = -1;
        }
        else {
          id = searchHash(Hash, pair);
          if (id == -1) { // new pair, insert
            id = insertRecord(&Rec, pair);
            rec = &Rec.records[id];
            L[ant].next = -1;
          }
          else {
            incFreq(&Heap, id);
            rec = &Rec.records[id];
            L[ant].next = rec->cpos;
            L[L[ant].next].prev = ant;
          }
          L[ant].prev = -id - 1;
          rec->cpos = ant;
        }
      }
      C[cpos] = n;
      if (ssgte != u)
        C[ssgte - 1] = -cpos - 1;
      C[cpos + 1] = -ssgte - 1;
      c--;
      orec = &Rec.records[oid]; // just in case of Rec.records realloc'd
      cpos = orec->cpos;
    }
    if (PRNC)
      prnC();
    removeRecord(&Rec, oid);
    n++;
    purgeHeap(&Heap);   // remove freq 1 from heap
    if (c < factor * u) { // compact C
      relong i, ni;
      i = 0;
      for (ni = 0; ni < c - 1; ni++) {
        C[ni] = C[i];
        L[ni] = L[i];
        if (L[ni].prev < 0) {
          if (L[ni].prev != NullFreq) // real ptr
            Rec.records[-L[ni].prev - 1].cpos = ni;
        }
        else
          L[L[ni].prev].next = ni;
        if (L[ni].next != -1)
          L[L[ni].next].prev = ni;
        i++;
        if (C[i] < 0)
          i = -C[i] - 1;
      }
      C[ni] = C[i];
      u = c;
      C = realloc(C, c * sizeof(relong));
      L = realloc(L, c * sizeof(Tlist));
      assocRecords(&Rec, &Hash, &Heap, L);
    }
  }
  return 0;
}


int main(int argc, char **argv)
{
  extern char *optarg;
  extern int optind, opterr, optopt;
  char fname[1024];
  FILE *Tf, *Rf, *Cf;
  relong i, olen, len;
  struct stat s;

  /* ----- -------- read options from command line ----------- */
  int save_space=0;
  maxRules= INT64_MAX;
  maxForbiddenChar= -1;
  Verbose=0;
  opterr = 0;
  while ((c=getopt(argc, argv, "r:x:vs")) != -1) {
    switch (c) {
    case 'v':
      Verbose++;
      break;
    case 's':
      save_space=1;
      break;
    case 'r':
      maxRules=atol(optarg);
      break;
    case 'x':
      maxForbiddenChar=atoi(optarg);
      if(maxForbiddenChar<0 || maxForbiddenChar>255) {
        fprintf(stderr,"-x option must be in [0-255]\n");
        exit(1);
      }
      break;
    case '?':
      fprintf(stderr,"Unknown option: %c\n", optopt);
      exit(1);
    }
  }
  if(Verbose>0) {
    fputs("==== Command line:\n",stderr);
    for(int i=0; i<argc; i++)
      fprintf(stderr," %s",argv[i]);
    fputs("\n",stderr);
  }
  if(maxRules<0) {
    fprintf(stderr,"-r option must be > 0\n");
    exit(1);
  }

  // virtually get rid of options from the command line
  optind -=1;
  if (argc-optind != 2) usage_and_exit(argv[0]);
  argv += optind;
  argc -= optind;
  
  // ----- read input
  if (stat(argv[1], &s) != 0) {
    fprintf(stderr, "Error: cannot stat file %s\n", argv[1]);
    exit(1);
  }
  olen = len = s.st_size;
  Tf = fopen(argv[1], "r");
  if (Tf == NULL) {
    fprintf(stderr, "Error: cannot open file %s for reading\n", argv[1]);
    exit(1);
  }
  uint8_t *text8 = mymalloc(len * sizeof(char),__LINE__,__FILE__);
  if (fread(text8, 1, len, Tf) != len) {
    fprintf(stderr, "Error: cannot read file %s\n", argv[1]);
    exit(1);
  }
  fclose(Tf);
  // ---------- open .R file
  strcpy(fname, argv[1]);
  strcat(fname, ".R");
  Rf = fopen(fname, "w");
  if (Rf == NULL) {
    fprintf(stderr, "Error: cannot open file %s for writing\n", fname);
    exit(1);
  }


  // ------------- stage 0  (init)
  init(text8, len, Rf);
  // ------------- stage 1  (if option -s is specified)
  if (save_space>0) {
    // do the quadratic phase using uint16_t
    uint16_t *text16 = myrealloc(text8, len * sizeof(uint16_t),__LINE__,__FILE__);
    text8 = (uint8_t *) text16;
    // convert inplace from char to short
    for (long i = len-1; i >=0; i--) {
      uint16_t c = text8[i];
      assert(c<256);
      text16[i] = c;
    }
    // work on the short version as much as possible 
    len = repair1(text16, len, Rf);
    if (PRNL) printf("---end of stage 1, size=%lld\n",len);
    assert(len>0);
    text16 = myrealloc(text16, len * sizeof(uint16_t),__LINE__,__FILE__);
    prepare(NULL, text16, len); // text16 is freed here
  }
  else {
    // go directly from char to relong 
    if (PRNL) printf("---skipping stage 1\n");
    prepare(text8, NULL, len); // text8 is freed here
  }
  // ------------ final stage
  repair(Rf);
  if (fclose(Rf) != 0) {
    fprintf(stderr, "Error: cannot close file %s\n", fname);
    exit(1);
  }
  strcpy(fname, argv[1]);
  strcat(fname, ".C");
  Cf = fopen(fname, "w");
  if (Cf == NULL) {
    fprintf(stderr, "Error: cannot open file %s for writing\n", fname);
    exit(1);
  }
  i = 0;
  while (i < u) {
    int cc = (int)C[i];
    if (fwrite(&cc, sizeof(int), 1, Cf) != 1) {
      fprintf(stderr, "Error: cannot write file %s\n", fname);
      exit(1);
    }
    i++;
    if ((i < u) && (C[i] < 0))
      i = -C[i] - 1;
  }
  if (fclose(Cf) != 0) {
    fprintf(stderr, "Error: cannot close file %s\n", fname);
    exit(1);
  }
  if (PRNCf)
    prnC();
  // free what is possible
  free(C); free(L);  

  // ------------- report compression statistics

  // n is the highest numbered symbol, since each rule introduces a new symbol
  // starting with alph, n-alph is the number of rules.
  // output size: 2.0*rules for the grammar tree shape
  //              rules+c*log(n-1) for the encoding of the tree leaves
  //                               and the sequence C
  long est_size = (long)((2.0 * (n - alph) + ((n - alph) + c) * (float)blog(n - 1)) / 8) + 1;
  fprintf(stderr, "RePair succeeded\n");
  fprintf(stderr, "   Original chars: %lli\n", olen);
  fprintf(stderr, "   Number of rules: %i\n", n - alph);
  fprintf(stderr, "   Final sequence length: %lli (integers)\n", c);
  fprintf(stderr, "   Estimated output size (bytes): %ld\n", est_size);
  fprintf(stderr, "   Estimated compression ratio: %0.2f%%\n", (100.0 * est_size) / olen);
  // original estimate (4.0*(n-alph)+((n-alph)+c)*(float)blog(n-1))/(olen*8.0)*100.0);
  exit(0);
}



// functions printing debug information 

void prnSym(int c)
{
  if (c < alph)
    printf("%c", c);
  else
    printf("%i", c);
}


void prnsC(unsigned short *sC, relong len)
{
  relong i = 0;
  printf("C[1..%lli] = ", len);
  while (i < len) {
    prnSym(sC[i]);
    printf(" ");
    i++;
  }
  printf("\n\n");
}


void prnC(void)
{
  relong i = 0;
  printf("C[1..%lli] = ", c);
  while (i < u) {
    prnSym((int)C[i]);
    printf(" ");
    i++;
    if ((i < u) && (C[i] < 0))
      i = -C[i] - 1;
  }
  printf("\n\n");
}

void prnRec(void)
{
  int i;
  printf("Active pairs:\n");
  for (i = 0; i < Rec.size; i++) {
    printf("\t(");
    prnSym(Rec.records[i].pair.left);
    printf(",");
    prnSym(Rec.records[i].pair.right);
    printf("), %lli occs\n", Rec.records[i].freq);
  }
  printf("\n");
}

static void usage_and_exit(char *name)
{
  fprintf(stderr,"Usage:\n\t  %s [options] infile MB\n",name);
  fprintf(stderr,"\t\t-v             verbose\n");
  fprintf(stderr,"\t\t-s             slower but save some space \n");
  fprintf(stderr,"\t\t-r num         max number of rules (def. no limit)\n");
  fprintf(stderr,"\t\t-x mul         largest char ignored by rules def. none\n\n");
  exit(1);
}

// write error message and exit
static void quit(const char *s)
{
  if(errno==0) fprintf(stderr,"%s\n",s);
  else  perror(s);
  exit(1);
}

// malloc and exit if out of memory
static void *mymalloc(size_t size, int line, const char *file)
{
  void *v=malloc(size);
  if(v==NULL) {
    fprintf(stderr,"Out of memory allocating %zu bytes" 
                 "at line %d of file %s\n",size,line,file);
    exit(3);
  }
  return v;
}

static void *myrealloc(void *ptr, size_t size, int line, const char *file)
{
  void *v=realloc(ptr,size);
  if(v==NULL) {
    fprintf(stderr,"Out of memory allocating %zu bytes" 
                 "at line %d of file %s\n",size,line,file);
    exit(3);
  }
  return v;
}

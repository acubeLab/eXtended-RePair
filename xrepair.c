/* repairX.c
   Giovanni Manzini    June 20st, 2023

  Implementation of repair algorithm for characters with bounded memory
  by Gonzalo Navarro with the following additional options:
    -x C    [do not generate rules involving chars <= C]
    -r R    [do not generate more than R rules]
    -m M    [use a quadratic algorithm until the full linear time
             machinery fits in M MBs ]
  The input file name is a required input parameter. 
  The program produces two files:
    - a file with extension .R containing the rules
    - a file with extension .C containing the compressed text
  The option -m can greatly affect the running time, it should not 
  the affect the output.

  This code is derived by the large and balanced repair implementation,
  which used the following approach: since the working memory 
  of the linear time algorithm is dominated by the term 
     tlen * 3 * sizeof(long long) [=24*tlen]   
  until that size fits in the available memory we use a quadratic update 
  (once a rule is found it is applied with a complete scan on the input). 
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
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include <stdbool.h>

#include "basics.h"
#include "records.h"
#include "hash.h"
#include "heap.h"


// debug control
#define PRNC  (Verbose>3)  // print current sequence C (verbose!)
#define PRNR  (Verbose>3)  // print active pairs in the heap (verbose!) 
#define PRNCf (Verbose>2)  // print final sequence C
#define PRNP  (Verbose>1)  // print forming pairs
#define PRNL  (Verbose>1)  // print progress on text scan
static void prnSym(uint32_t c);
static void prnsC(uint32_t *, reIdx len);
static void prnC(void);
static void prnRec(void);
static void usage_and_exit(char *name);
void writeAlpha(size_t a, FILE *f, int line, char *file);
void writeRule(Tpair p, FILE *f, int line, char *file);


// ugly globals
float factor = 0.75; // 1/extra space overhead; set closer to 1 for smaller
                     // and slower execution
int minsize = 256;   // to avoid many reallocs at small sizes, should be ok as is

reIdx *C; // compressed text, pero tiene -ptrs al mismo texto

static reIdx uLen; // |text| and later current |C| with gaps
static reIdx cSize; // real |C| (without gaps)

ssize_t alph; // alphabet size and first non terminal symbol 
              // it is the first item written to the .R file as a uint32
ssize_t n;    // alph + |R| = next id for a non terminal symbol

Tlist *L; // |L| = c; // next and previous entry in gapped |C|

Thash Hash;  // hash table of pairs
Theap Heap;  // special heap of pairs
Trarray Rec; // records


// globals controlled by the input parameters
static int64_t maxForbiddenChar;  // do not generate rules involving chars <= this
                                  // default is -1 which means no restriction
static int64_t maxRules;          // max number of rules to generate
                                  // default INT64_MAX, i.e. no restriction
static int32_t maxMB;             // available memory in MBs 
                                  // default INT32_MAX, i.e. no restriction
static int Verbose;



// return true if forbidden chars are enabled 
// and this pair of symbols should never appear
// as the left-hand side of a rule
int forbidden_pair(reSym left, reSym right)
{
  assert(left>=0 && right >=0); // useless as long as reSym is unsigned 
  return (maxForbiddenChar>=0) && 
         ((left <= maxForbiddenChar) || (right <= maxForbiddenChar));
}


// init32: given the input inside a uint32 sequence compute alphabet size 
// and init Rec, Hash, and Heap. If expand==true the input is in
// the first len bytes and must be expanded
void init32(reSym *text, bool expand, reIdx len, FILE *R)
{
  reIdx i;
  int id;
  Tpair pair;
  // compute largest input symbol: if expand==true we must 
  // also move the uint8 in the corresponding uint32 position
  assert(len>0);
  alph = 0;
  if(expand) {
    uint8_t *text8 = (uint8_t *) text;
    for(i=len-1;i>=0;i--) {
      uint32_t ctmp = text8[i];
      text[i]=ctmp; 
      if(ctmp>alph) alph =ctmp;
    }
  }
  else for (i = 0; i < len; i++) {
    if (text[i] > alph)
      alph = text[i];
  }
  
  // init n as first code usable for non terminals and 
  // save alphabet size to the R file
  n = ++alph;
  writeAlpha(alph,R,__LINE__,__FILE__);

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
        printf("Processed %zd chars\n", i);
    }
  }
  purgeHeap(&Heap); // remove pairs with freq==1
}


// prepare for the use of the full machinery:
//   copy text[] to C[], free text[] 
//   init the prev/next array L[]
void prepare(uint32_t *text, reIdx len)
{
  reIdx i;
  int id;
  Tpair pair;
  cSize = uLen = len;
  C = mymalloc(uLen * sizeof(reIdx),__LINE__,__FILE__);
  assert(text!=NULL);
  for (i = 0; i < uLen; i++)
    C[i] = text[i]; // copy from text array
  free(text); 

  // init prev/next list L
  L = mymalloc(uLen * sizeof(Tlist),__LINE__,__FILE__);
  assocRecords(&Rec, &Hash, &Heap, L);
  for (i = 0; i < cSize - 1; i++) {
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
      printf("Processed %zi chars\n", i);
  }
  L[i].prev = NullFreq;
  L[i].next = -1;
  purgeHeap(&Heap); // probably not necessary since we did not modify the heap
}


// first pass consisting in the repair algorithm done using just an
// uint32_t array and a quadratic update algorithm
// called if -m option was used and until the the full 24*len structure
// fits in the available memory
reIdx repair32q(reSym *sC, reIdx len, FILE *R)
{
  int oid, id;
  reIdx cpos, pos;
  Trecord *orec;
  Tpair pair;
  int left, right;
  if (PRNC)
    prnsC(sC,len);
  // n is the id of the next rule, n-alpha the number of rules so far
  // stop if rule id becomes invalid or we reached the desired number of rules
  while (n<=reSym_MAX && (n-alph) < maxRules) {
    if ((len / 1024 / 1024) * 3 * sizeof(reIdx) < maxMB)
      return len;
    if (PRNR)
      prnRec();
    oid = extractMax(&Heap);
    if (oid == -1)
      break; // the end: no more pairs appearing more than once
    orec = &Rec.records[oid];
    writeRule(orec->pair,R,__LINE__,__FILE__); // write rule to .R file 
    left = orec->pair.left;
    right = orec->pair.right;
    if (PRNP) {
      printf("Chosen pair %zi = (", n);
      prnSym(orec->pair.left);
      printf(",");
      prnSym(orec->pair.right);
      printf(") (%zi occs)\n", orec->freq);
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


// linear time repair algorithm. 
// to avoid another list to access the sparse C we thread it using the
// empty space. if next cell of an active cell is negative, it is
// (minus) a ptr to the next occ. idem previous cell to previous occ,
// except that next ptr dominates over prev ptr if they must be in
// the same cell. but in this case one can find prev in O(1) anyway.
reIdx repair(FILE *R)
{
  int oid, id;
  reIdx cpos;
  Trecord *rec, *orec;
  Tpair pair;
  if (Verbose>0)
    printf("--- final stage, n=%zi\n", cSize);
  if (PRNC)
    prnC();
  // n is the id of the next rule, n-alpha the number of rules so far
  // stop if rule id becomes invalid or we reached the desired number of rules
  while (n<=reSym_MAX && (n-alph) < maxRules) {
    if (PRNR)
      prnRec();
    oid = extractMax(&Heap);
    if (oid == -1)
      break; // the end!!
    orec = &Rec.records[oid];
    cpos = orec->cpos;
    writeRule(orec->pair,R,__LINE__,__FILE__);
    if (PRNP) {
      printf("Chosen pair %zi = (", n);
      prnSym(orec->pair.left);
      printf(",");
      prnSym(orec->pair.right);
      printf(") (%zi occs)\n", orec->freq);
    }
    while (cpos != -1) {
      reIdx ant, sgte, ssgte;
      // replacing bc->e in abcd, b = cpos, c = sgte, d = ssgte
      if (C[cpos + 1] < 0)
        sgte = -C[cpos + 1] - 1;
      else
        sgte = cpos + 1;
      if ((sgte + 1 < uLen) && (C[sgte + 1] < 0))
        ssgte = -C[sgte + 1] - 1;
      else
        ssgte = sgte + 1;
      // remove bc from L
      if (L[cpos].next != -1)
        L[L[cpos].next].prev = -oid - 1;
      orec->cpos = L[cpos].next;
      if (ssgte != uLen) { // there is ssgte
        // remove occ of cd
        assert(C[sgte]>=0 && C[sgte]<=reSym_MAX);
        assert(C[ssgte]>=0 && C[ssgte]<=reSym_MAX);
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
        assert(C[ant]>=0 && C[ant]<=reSym_MAX);
        assert(C[cpos]>=0 && C[cpos]<=reSym_MAX);
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
      if (ssgte != uLen)
        C[ssgte - 1] = -cpos - 1;
      C[cpos + 1] = -ssgte - 1;
      cSize--;
      orec = &Rec.records[oid]; // just in case of Rec.records realloc'd
      cpos = orec->cpos;
    }
    if (PRNC)
      prnC();
    removeRecord(&Rec, oid);
    n++;
    purgeHeap(&Heap);   // remove freq 1 from heap
    if (cSize < factor * uLen) { // compact C
      reIdx i, ni;
      i = 0;
      for (ni = 0; ni < cSize - 1; ni++) {
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
      uLen = cSize;
      C = myrealloc(C, cSize * sizeof(reIdx),__LINE__,__FILE__);
      L = myrealloc(L, cSize * sizeof(Tlist),__LINE__,__FILE__);
      assocRecords(&Rec, &Hash, &Heap, L);
    }
  }
  return 0;
}


int main(int argc, char **argv)
{
  extern char *optarg;
  extern int optind, opterr, optopt;
  char fname[PATH_MAX];
  FILE *Tf, *Rf, *Cf;
  reIdx i, olen, len;
  struct stat s;
  int o;

  /* ----- -------- read options from command line ----------- */
  maxRules= INT64_MAX;
  maxForbiddenChar= -1;
  maxMB = INT32_MAX;
  Verbose=0;
  opterr = 0;
  while ((o=getopt(argc, argv, "m:r:x:v")) != -1) {
    switch (o) {
    case 'v':
      Verbose++;
      break;
    case 'm':
      maxMB=atoi(optarg);
      if(maxMB<0) {
        fprintf(stderr,"-m option must be non negative\n"); exit(1);
      }
      break;
    case 'r':
      maxRules=atol(optarg);
      if(maxRules<0) {
        fprintf(stderr,"-r option must be non negative\n"); exit(1);
      }
      break;
    case 'x': {
      long maxfc =atol(optarg);
      if(maxfc<0 || maxfc> reSym_MAX) {
        fprintf(stderr,"-x option must be between 0 and 2^32-1\n");
        exit(1);
      } else maxForbiddenChar = maxfc;
      break; }
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
  // store input in an array of int32
  reSym *text32 = mymalloc(len * sizeof(uint32_t),__LINE__,__FILE__);
  if (fread(text32, 1, len, Tf) != len) {
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
  init32(text32, true, len, Rf);
  // ------------- stage 1  (if necessary)
  if ((len/1024/1024)*3*sizeof(reIdx)>= maxMB) {
    // work on the short version as much as possible 
    len = repair32q(text32, len, Rf);
    if (Verbose>0) fprintf(stderr,"--- end of stage 1, size=%zi\n",len);
    assert(len>0);
  }
  else if (Verbose>0) fprintf(stderr,"--- skipping stage 1\n");

  // --- final stage (linear time but using 24*len)
  prepare(text32, len); // init C[] and L[], text32 is freed here
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
  while (i < uLen) {
    int cc = (int)C[i];
    if (fwrite(&cc, sizeof(int), 1, Cf) != 1) {
      fprintf(stderr, "Error: cannot write file %s\n", fname);
      exit(1);
    }
    i++;
    if ((i < uLen) && (C[i] < 0))
      i = -C[i] - 1;
  }
  if (fclose(Cf) != 0) {
    fprintf(stderr, "Error: cannot close file %s\n", fname);
    exit(1);
  }
  if (PRNCf) prnC();
  // free what is possible
  free(C); free(L);  


  // ------------- report compression statistics

  // n is the highest numbered symbol, since each rule introduces a new symbol
  // starting with alph, n-alph is the number of rules.
  // output size: 2.0*rules for the grammar tree shape
  //              rules+c*log(n-1) for the encoding of the tree leaves
  //                               and the sequence C
  long est_size = (long)((2.0 * (n - alph) + ((n - alph) + cSize) * (float)blog(n - 1)) / 8) + 1;
  fprintf(stderr, "RePair succeeded\n");
  fprintf(stderr, "   Original chars: %zi\n", olen);
  fprintf(stderr, "   Number of rules: %zi\n",  (n - alph));
  fprintf(stderr, "   Final sequence length: %zi (integers)\n", cSize);
  fprintf(stderr, "   Estimated output size (bytes): %ld\n", est_size);
  fprintf(stderr, "   Estimated compression ratio: %0.2f%%\n", (100.0 * est_size) / olen);
  // original estimate (4.0*(n-alph)+((n-alph)+c)*(float)blog(n-1))/(olen*8.0)*100.0);
  exit(0);
}



// functions printing debug information 
static void prnSym(uint32_t c)
{
  if (c < 256)
    printf("%c", c);
  else
    printf("%u", c);
}


static void prnsC(uint32_t *sC, reIdx len)
{
  reIdx i = 0;
  printf("C[1..%zi] = ", len);
  while (i < len) {
    prnSym(sC[i]);
    printf(" ");
    i++;
  }
  printf("\n\n");
}


static void prnC(void)
{
  reIdx i = 0;
  printf("C[1..%zi] = ", cSize);
  while (i < uLen) {
    prnSym((uint32_t)C[i]);
    printf(" ");
    i++;
    if ((i < uLen) && (C[i] < 0))
      i = -C[i] - 1;
  }
  printf("\n\n");
}

static void prnRec(void)
{
  int i;
  printf("Active pairs:\n");
  for (i = 0; i < Rec.size; i++) {
    printf("\t(");
    prnSym(Rec.records[i].pair.left);
    printf(",");
    prnSym(Rec.records[i].pair.right);
    printf("), %zi occs\n", Rec.records[i].freq);
  }
  printf("\n");
}


// --------- output functions

// write (input) alphabet size (usually to R) as an uint32_t 
// shall we change it if we use 5nytes per symbol? probably not...
void writeAlpha(size_t a, FILE *f, int line, char *file) {
  if(a<0 || a>UINT32_MAX) 
    quit("Input alphabet negative or larger than 2^32-1",line,file);
  uint32_t t = (uint32_t) a;
  if(fwrite(&t,sizeof(t),1,f)!=1)
    quit("Error writing alphabet size to file",line,file);
}

// write a rule (a pair of nonterminals) currently as a 
// pair of uint32, maybe in the future using 5+5 bytes
void writeRule(Tpair p, FILE *f, int line, char *file) {
  // !!! to be changed if reSym changes 
  static_assert(sizeof(p)==8,"Tpair should consist of two uint32s");
  if(fwrite(&p,sizeof(p),1,f)!=1)
    quit("Error writing rule to file", line, file); 
}

// write a single symbol (usually to the .C file) currently as a
// uint32, maybe in the future using 5 bytes 
void writeSymbol(reIdx s , FILE *f, int line, char *file) {
  if(s<0 || s>reSym_MAX) 
    quit("Symbol negative or larger than reSym_MAX",line,file);
  // !!! to be changed if reSym changes   
  static_assert(sizeof(reSym)==4,"reSym should consist of an uint32");
  uint32_t t = (uint32_t) s;
  if(fwrite(&t,sizeof(t),1,f)!=1)
    quit("Error writing symbol to file", line, file); 
}



static void usage_and_exit(char *name)
{
  fprintf(stderr,"Usage:\n\t  %s [options] infile\n",name);
  fprintf(stderr,"\t\t-v             verbose\n");
  fprintf(stderr,"\t\t-m maxMB       max memory to use in MB (def. no limit)\n");
  fprintf(stderr,"\t\t-r num         max number of rules (def. no limit)\n");
  fprintf(stderr,"\t\t-x mul         largest char ignored by rules (def. none)\n\n");
  exit(1);
}

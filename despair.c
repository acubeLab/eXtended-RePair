// decompress a char-based .R .C pair

/*

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
// type definitions anche checks for type consistency are in basics.h
#include "basics.h"

static void usage_and_exit(const char *name);

// static reIdx u;    // |text| and later current |C| with gaps
static reIdx alph; // size of terminal alphabet and smallest non terminal symbol
Tpair *R;          // rules (from .R file)
static reIdx n; // |R|

static FILE *outfile;
static size_t maxdepth = 0;

// globals controlled by the input parameters
static int Verbose;



reIdx expand (reIdx i, size_t d)
{
  size_t ret = 1;
  while (i >= alph) { // while i is not a terminal expand recursively
    ret += expand(R[i-alph].left,d+1);
    i = R[i-alph].right;
    d++;  // expansion on the right branch is replaced by iteration
  }
  // here we are using that the output is a uint8
  assert(i<256);
  char c = i;
  if (fwrite(&c,sizeof(char),1,outfile) != 1) 
    quit("Cannot write to output file", __LINE__, __FILE__);
  if (d > maxdepth) maxdepth = d;// keep track of max depth
  return ret;
}

int main (int argc, char **argv)
{
  extern char *optarg;
  extern int optind, opterr, optopt;
  char fname[PATH_MAX];
  char outname[PATH_MAX];
  //char *text;
  FILE *Rf,*Cf;
  unsigned int i;
  reIdx len,c,u;
  struct stat s;
  int o;
  
  Verbose=0;
    opterr = 0;
  while ((o=getopt(argc, argv, "v")) != -1) {
    switch (o) {
    case 'v':
      Verbose++;
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
  // virtually get rid of options from the command line
  optind -=1;
  if (argc-optind != 2) usage_and_exit(argv[0]);
  argv += optind;
  argc -= optind;
  
  // read .R file, store data in alph and R[]
  strcpy(fname,argv[1]);
  strcat(fname,".R");
  if (stat (fname,&s) != 0) {
    fprintf (stderr,"Error: cannot stat file %s\n",fname);
    exit(1);
  }
  len = s.st_size;
  Rf = fopen (fname,"r");
  if (Rf == NULL) {
    fprintf (stderr,"Error: cannot open file %s for reading\n",fname);
    exit(1);
  }
  // currently the alphabet size is stored into an int32_t 
  uint32_t tmp;
  if (fread(&tmp,sizeof(tmp),1,Rf) != 1) {
    fprintf (stderr,"Error: cannot read file %s\n",fname);
    exit(1);
  }
  else {
    alph = tmp;
    fprintf (stderr,"Alphabet size: %zd\n",alph);
  }
  // note that in the original char-based repair the R file contains also
  // a map between the 0...alph-1 and the actual symbols in the input file
  // here there is no such map
  
  // read and store .R file
  // n is the number of rules, sizeof(tmp) accounts for alpha
  if( (len-sizeof(tmp)) % sizeof(Tpair) !=0) 
    quit("Invalid format for the .R file",__LINE__,__FILE__);
  n = (len-sizeof(tmp))/sizeof(Tpair);
  // allocate and read array of rules stored as pairs
  R = mymalloc(n*sizeof(Tpair),__LINE__,__FILE__);
  if (fread(R,sizeof(Tpair),n,Rf) != n) {
    fprintf (stderr,"Error: cannot read file %s\n",fname);
    exit(1);
  }
  fclose(Rf);

  // open C file and get the number of symbols in it
  strcpy(fname,argv[1]);
  strcat(fname,".C");
  if (stat (fname,&s) != 0) {
    fprintf (stderr,"Error: cannot stat file %s\n",fname);
    exit(1);
  }
  // from now on the size of a symbol in the C file
  // is assumed to be 4 bytes, it can change in the future...
  if(s.st_size%4!=0)
    quit("Invalid format for the .R file",__LINE__,__FILE__); 
  c = len = s.st_size/4;
  Cf = fopen (fname,"r");
  if (Cf == NULL) {
    fprintf (stderr,"Error: cannot open file %s for reading\n",fname);
    exit(1);
  }

  // open output file
  strcpy(outname,argv[1]);
  strcat(outname,".out");
  outfile = fopen (outname,"w");
  if (outfile == NULL) {
    fprintf (stderr,"Error: cannot open file %s for writing\n",outname);
    exit(1);
  }

  // actual decompression
  u = 0;
  for (; len>0; len--) {
    if (fread(&i,sizeof(unsigned int),1,Cf) != 1)
      quit("Canot read from .C file",__LINE__,__FILE__);
    u += expand(i,0); // expand non terminal i, 0 is initial depth
  }
  if(fclose(Cf)!=0)
    quit("Cannot close ,C file",__LINE__,__FILE__);
  if (fclose(outfile) != 0) 
    quit("Cannot close output file",__LINE__,__FILE__);

  // here n is the number of rules, n+alpha the effective alphabet in C
  // output size: 2.0*rules for the grammar tree shape
  //              (rules+c)*log(n+alpha) for the encoding of the tree leaves
  //                                     and the sequence C
  long est_size = (long) ( (2.0*n+(n+c)*(double)blog(n+alph))/8) + 1;
  fprintf (stderr,"DesPair succeeded\n");
  fprintf (stderr,"   Original chars: %ld\n",u);
  fprintf (stderr,"   Number of rules: %ld\n",n);
  fprintf (stderr,"   Compressed sequence length: %zd (integers)\n",c);
  fprintf (stderr,"   Maximum rule depth: %zd\n",maxdepth);
  fprintf (stderr,"   Estimated output size (bytes): %ld\n",est_size);
  fprintf (stderr,"   Compression ratio: %0.2f%%\n", (100.0* est_size)/u);
  return 0;
}


static void usage_and_exit(const char *name)
{
  fprintf(stderr,"Usage:\n\t  %s [options] name\n\n",name);
  fprintf(stderr,"\tRepair grammar decompression algorithm\n");
  fprintf(stderr,"\tRead name.C and name.R and decompress them to name.out\n\n");
  fprintf(stderr,"\t\t-v             verbose\n\n");
  exit(1);
}


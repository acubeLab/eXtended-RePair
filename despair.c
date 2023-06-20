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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

typedef struct {
  unsigned int left,right;
} Tpair;

long u; // |text| and later current |C| with gaps


unsigned int alph; // size of terminal alphabet, or smallest non terminal symbol

Tpair *R; // rules

size_t n; // |R|
// int n; // |R|

char *ff;
FILE *f;
size_t maxdepth = 0;

size_t expand (unsigned int i, size_t d)

{
  size_t ret = 1;
  char c;
  while (i >= alph) { // while i is not a terminal expand recursively
    ret += expand(R[i-alph].left,d+1);
    i = R[i-alph].right;
    d++;  // expansion on the right branch is replaced by iteration
  }
  c = i;
  if (fwrite(&c,sizeof(char),1,f) != 1) {
    fprintf (stderr,"Error: cannot write file %s\n",ff);
    exit(1);
  }
  if (d > maxdepth) maxdepth = d;// keep track of max depth
  return ret;
}

static int bits (size_t x)
{
  int l=0;
  while (x) {
    x>>=1;
    l++;
  }
  return l;
}

int main (int argc, char **argv)

{
  char fname[1024];
  char outname[1024];
  //char *text;
  FILE *Tf,*Rf,*Cf;
  unsigned int i;
  size_t len,c,u;
  struct stat s;
  fputs("==== Command line:\n",stderr);
  for(int i=0; i<argc; i++)
    fprintf(stderr," %s",argv[i]);
  fputs("\n",stderr);
  if (argc != 2) {
    fprintf (stderr,"Usage: %s <filename>\n"
             "Decompresses <filename> from its .C and .R "
             "extensions.\n Decompressed file is <filename>.out\n"
             "This is a version for prefix-free parsing\n",argv[0]);
    exit(1);
  }

  // read .R file, store data in alpha and R[]
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
  if (fread(&alph,sizeof(int),1,Rf) != 1) {
    fprintf (stderr,"Error: cannot read file %s\n",fname);
    exit(1);
  }
  else {
    fprintf (stderr,"Alphabet size: %d\n",alph);
  }
  // note that in the original char-based repair the R file contains also
  // a map between the 0...alph-1 and the actual symbols in the input file
  // here alph is 256 and there is no such map (this is why the two .R .C
  // formats are not compatible).

  // n is the number of rules, sizeof(int) accounts for alpha
  n = (len-sizeof(int))/sizeof(Tpair);
  // allocate and read array of rules stored as pairs
  R = (void*)malloc(n*sizeof(Tpair));
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
  c = len = s.st_size/sizeof(unsigned int);
  Cf = fopen (fname,"r");
  if (Cf == NULL) {
    fprintf (stderr,"Error: cannot open file %s for reading\n",fname);
    exit(1);
  }

  // open output file
  strcpy(outname,argv[1]);
  strcat(outname,".out");
  Tf = fopen (outname,"w");
  if (Tf == NULL) {
    fprintf (stderr,"Error: cannot open file %s for writing\n",outname);
    exit(1);
  }

  // actual decompression
  u = 0;
  f = Tf;
  for (; len>0; len--) {
    if (fread(&i,sizeof(unsigned int),1,Cf) != 1) {
      fprintf (stderr,"Error: cannot read file %s\n",fname);
      exit(1);
    }
    u += expand(i,0); // expand non terminal i, 0 is initial depth
  }
  fclose(Cf);
  if (fclose(Tf) != 0) {
    fprintf (stderr,"Error: cannot close file %s\n",outname);
    exit(1);
  }
  // here n is the number of rules, n+alpha the effective alphabet in C
  long est_size = (long) ( (2.0*n+(n+c)*(double)bits(n+256))/8) + 1;
  fprintf (stderr,"DesPair succeeded\n");
  fprintf (stderr,"   Original chars: %ld\n",u);
  fprintf (stderr,"   Number of rules: %ld\n",n);
  fprintf (stderr,"   Compressed sequence length: %zd (integers)\n",c);
  fprintf (stderr,"   Maximum rule depth: %zd\n",maxdepth);
  fprintf (stderr,"   Estimated output size (bytes): %ld\n",est_size);
  fprintf (stderr,"   Compression ratio: %0.2f%%\n", (100.0* est_size)/u);
  return 0;
}


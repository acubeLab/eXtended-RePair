
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

#include "basics.h"


relong NullFreq = ((relong)1) << (8*sizeof(relong)-1);


// malloc and exit if out of memory
void *mymalloc(size_t size, int line, const char *file)
{
  void *v=malloc(size);
  if(v==NULL) {
    fprintf(stderr,"Out of memory allocating %zu bytes" 
                 "at line %d of file %s\n",size,line,file);
    exit(3);
  }
  return v;
}

void *myrealloc(void *ptr, size_t size, int line, const char *file)
{
  void *v=realloc(ptr,size);
  if(v==NULL) {
    fprintf(stderr,"Out of memory allocating %zu bytes" 
                 "at line %d of file %s\n",size,line,file);
    exit(3);
  }
  return v;
}

int blog (int x) { int l=0;
     while (x) { x>>=1; l++; }
     return l;
}

// write error message and exit
void quit(const char *msg, int line, char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);

  exit(1);
}






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
#ifndef BASICSINCLUDED
#define BASICSINCLUDED
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <assert.h>

// type used to represent an index in a sequence, size_t seems a safe choice
typedef ssize_t reIdx;
typedef size_t ureIdx;
#define reIdx_MAX INT64_MAX // replacement fo SSIZE_MAX which is not defined
// type used to represent a single symbol (terminal or non-terminal)
// currently uint32 but we may later go to 5 byte symbols 
typedef uint32_t reSym; 
#define reSym_MAX UINT32_MAX
// check type consistency: 
static_assert (sizeof(reIdx) >=8, "reIdx type must be at least 64 bits");
static_assert (sizeof(ureIdx) >=8, "ureIdx type must be at least 64 bits");
static_assert (sizeof(reSym) >=4, "reSym type must be at least 32 bits");
// reIdx is also used as a signed superset of reSym so we have extra constraints:
// reIdx must be signed and contain values larger than reSym_MAX  
static_assert (((reIdx) -1)  < 0, "reIdx must be a signed type");
static_assert (reIdx_MAX  > reSym_MAX, "reIdx must strictly contain than reSym");



void *mymalloc(size_t size, int line, const char *file);
void *myrealloc(void *ptr, size_t size, int line, const char *file);
void quit(const char *s, int, char *);

typedef struct { 
  reSym left,right;
} Tpair;

extern reIdx NullFreq;

int blog (int x); // bits to represent x

#endif

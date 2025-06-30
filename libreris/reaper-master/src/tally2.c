
/*
 * STALLED DEVELOPMENT.
 * tally2 was a stage in the separation of tally into two different programs,
 * namely (1) core tally functionality (tallying) (2) re-matching read pairs.
*/

/*
 * (C) Copyright 2011, 2012, 2013, 2014 European Molecular Biology Laboratory.
 * Author: Stijn van Dongen <stijn@ebi.ac.uk>.
 * Contact: <kraken@ebi.ac.uk>
 *
 * This file is part of Reaper.   Reaper is free software: you can redistribute
 * it  and/or modify it under the terms of the  GNU  General  Public License as
 * published by the Free Software Foundation;  either version 3 of the License,
 * or  (at your option)  any later version.  This program is distributed in the
 * hope that it will be useful,  but  WITHOUT  ANY  WARRANTY;  without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the  GNU  General  Public  License  for  more  details.  You should have
 * received a copy of the  GNU  General Public License along with this program.
 * If not, see http://www.gnu.org/licenses/.
*/


/*
   README

   If hash gets too full, a new hash is *added* rather
   than rehashing everything.  This works fine for up to 10 or even more hashes,
   whilst noting that such a scenario is highly unlikely to begin with.  The
   reasoning behind it is that rehashing will increase memory pressure, as
   bigger hash and old hash need to coexist at the same time.  This in the
   scenario that we do not really know how many unique reads there will be, and
   noting that the aim of this program is to work on lower-end machines.
      By and large these considerations are moot, as the hash size is settable
      from the command line.
         Even mooter, now by default tally will try to be smart and estimate
         file size and read content and set parameters accordingly.
         This can be turned off with --no-auto.

   -  Reads are stored entirely contiguously, in large chunks of memory.
      This helps efficient usage of memory (no fragmentation).
      Currently pointers to those chunks are discarded, as pointers to
      individual reads (kept in the hash) are all that we need.

   -  Reads are by default stored in compressed format - this leads to a space
      reduction of approximately 3-fold. Compressed format is described
      near compress_sequence.

   -  With --no-tally output is processed by writerecord rather than writeline,
      and output directives are different/more limited.
      The code logic depends on some continue statements.
      Tread with care.

   TODO

   -  test new inheritance code in cases of
         -i -j -o
         -i -o -p

   -  allow sequences of arbitrary length.

   ?  quality information; e.g. use unsigned and accumulate quality for 4 bases,
      average at the end, or simply keep unsigned for each base.
      60M unique reads of length 60 with 4 bytes -> 14.6G.

   ?  does hash function work well on uncompressed data?

   -  called-from-R has been rotting a bit. For example paired end not supported.
*/

#define SEP_PAIR_WIDTH 1

#ifndef TRACING_ON
#define TRACING_ON 0
#endif

#include "version.h"
#include "tally.h"
#include "slib.h"
#include "trint.h"
#include "dna.h"

#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <stdint.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>



#define DATA_SIZE_LOG_MAX 31
#define DATA_SIZE_LOG_MIN 14
#define HASH_SIZE_LOG_MAX 28
#define HASH_SIZE_LOG_MIN 14

#define DATA_SIZE_LOG_DEFAULT 28   /* If not cl-supplied and !auto 256M ish */
#define HASH_SIZE_LOG_DEFAULT 25   /* If not cl-supplied and !auto 32M  ish -- units, with unit size 12/16/20 bytes */



unsigned long g_memusage = 0;
unsigned long g_tally_fieldtosmall = 0;
unsigned g_stingy = 1;
unsigned g_verbose = 0;
unsigned long g_gauge_memory = 0;
unsigned g_have_paired = 0;
unsigned g_ntest = 100000;

unsigned n_bogon = 0;
long g_nt_left_bogon_diff = 0;    /* number of nt lost or gained due to bogons */
size_t g_nt_left_in     = 0;
size_t g_nt_right_in    = 0;
size_t g_nt_left_out    = 0;
size_t g_nt_right_out   = 0;




   /* encoded format is this:    concatenated segments, each of the form

               cookie {nt-cookie: content | nn-cookie: empty}

      *  total length of all segments is not stored in the format itself.
      *  there are two types of cookies: nt (nucleotide) and nn (sequence of Ns).
      *  a cookie is encoded in a single byte.
      *  the nt-cookie has bit 8 set.
      *  the nn-cookie has bit 7 set but bit 8 unset.
      *  the nt-cookie encodes a length L in the lower bits 1-7.
      *  the nn-cookie encodes a length L in the lower bits 1-6.
      *  the nt-cookie implies that L nucleotides are encoded next;
            L/4 bytes are consumed (integer division).
      *  the nn-cookie implies that a sequence of L 'N' follows; no bytes are consumed.
      *  a word of 4 nucletodides ACGT is encoded in bits as follows: TTGGCCAA,
         so in reverse order. This makes it easy to shift them out when decoding.

      Note that the format is not self-synchronising; the content
      following an nt-cookie is length encoded (by the cookie) and may
      contain both 0-bytes and cookie-bytes.
   */


int compress_sequence
(  char* seq
,  unsigned length
,  char* dest
)
   {  char* a = seq
   ;  char* z = seq + length;
   ;  char* d = dest, *cookie = d
   ;  unsigned runlength = 0, enclength = 0        /* runlength is length of current run of N or nt in input (uncompressed) space */
   ;  unsigned state = 0                           /* important (really ?) that it's not 'B' and it's not 'N' */

   ;  if (length)
         cookie[0] = 1 << (BASEMAP((unsigned char) a[0]) < 4 ? 7 : 6)
      ,  enclength++
      ,  state = 0

   ;  for (a=seq; a<z; a++)
      {  unsigned char c = a[0]
      ;  unsigned bm = BASEMAP(c)
      ;  if
         (  (bm < 4 && (state == 'N' || runlength == 128))
         || (bm >=4 && (state == 'B' || runlength == 64))
         )
         {  cookie[0] |= (runlength - 1)           /* write runlength to previous lookahead cookie */
         ;  cookie = ++d                           /* start new lookahead cookie */
         ;  cookie[0] = 1 << (bm < 4 ? 7 : 6)      /* bit 8 indicates nucleotide run, bit 7 indicates run of Ns */
         ;  enclength++                            /* encoding length increased by one (duh) */
         ;  runlength = 0                          /* started new run */
      ;  }
         if (bm < 4)
         {  if (!(runlength & 3))                  /* encoding length increases every 4 nucleotides */
               d++
            ,  d[0] = 0
            ,  enclength++
         ;  d[0] |= bm << (2 * (runlength & 3))    /* offset 0 -> bits {1,2}  1 -> bits { 3,4}  etc */
         ;  state = 'B'
      ;  }
         else
         state = 'N'
      ;  runlength++
   ;  }
      if (runlength)
      cookie[0] |= (runlength-1)
   ;  return enclength
;  }


int uncompress_sequence
(  const char* cpr
,  unsigned length
,  char* dest
)
   {  const char* a = cpr
   ;  const char* z = cpr + length
   ;  char* d = dest
   ;  unsigned runlength = 0, seqlength = 0
   ;  unsigned char pamesab[4] = { 'A', 'C', 'G', 'T' }

   ;  for (a=cpr; a<z; a++)
      {  unsigned char c = a[0]
;if(0)fprintf(stderr, "byte %X\n", (unsigned) c)
      ;  if (!runlength)                              /* read cookie */
         {  unsigned have_Ns = 0

         ;  if (a[0] & (1 << 7))                      /* run of nucleotides */
            runlength = 1 + (a[0] & 0x7F)
         ;  else if (a[0] & (1 << 6))                 /* run of Ns */
               have_Ns = 1
            ,  runlength = 1 + (a[0] & 0x3F)
         ;  else
            X_FINISHED_VAL(-1, "uncompress_sequence corruption at positions %u/%u", (unsigned) (a-cpr), seqlength)
;if(0)fprintf(stderr, "runlength %d\n", (int) runlength)
         ;  if (have_Ns)
            {  unsigned i = 0
            ;  for (i=0; i<runlength; i++)
               d++[0] = 'N'
            ;  seqlength += runlength
            ;  runlength = 0
         ;  }
      ;  }
         else
         {  unsigned n_decode = runlength > 3 ? 4 : runlength
         ;  switch(n_decode)
            {  case 4: d++[0] = pamesab[c & 3]; c >>= 2     /* use switch fall-through */
            ;  case 3: d++[0] = pamesab[c & 3]; c >>= 2
            ;  case 2: d++[0] = pamesab[c & 3]; c >>= 2
            ;  case 1: d++[0] = pamesab[c & 3]
         ;  }
            seqlength += n_decode
         ;  runlength -= n_decode
      ;  }
      }
      if (runlength)
      X_FINISHED_VAL(-1, "uncompress_sequence corruption at end")
   ;  return seqlength
;  }



#if 1

unsigned hash_sequence
(  char* s
,  unsigned length
)
   {  unsigned i, h=0
   ;  for (i=0;i<length;i++)
      h = 137 + ((h << 6) ^ (h << 13) ^ (h >> 21) ^ (s[i] * 71523) ^ (11 << BASEMAP((unsigned char) s[i])))
      // h = 11 + ((h << 6) ^ (h << 13) ^ (h >> 21) ^ (((uint32_t*) (void*) (s+4*i))[0] * 71523))
                  // (523 << BASEMAP((unsigned char) s[i])))  // (s[i] * 71523))   // BASEMAP((unsigned char) s[i]))  /* (s[i] * 71523)) */
   ;  return h
;  }

#else

unsigned hash_sequence
(  char* s
,  unsigned length
)
   {  unsigned i, h=0
   ;  for (i=0;i<length;i++)
      h = 0 + ((h << 6) ^ (h << 21) ^ (h >> 13) ^ (s[i] * 71523))
                  // (523 << BASEMAP((unsigned char) s[i])))  // (s[i] * 71523))   // BASEMAP((unsigned char) s[i]))  /* (s[i] * 71523)) */
   ;  return h >> 6
;  }

#endif

struct unit1
{  const char* nt       /* initialised to NULL */
;  unsigned length      /* if reads are compressed, this encodes compressed length */
;  unsigned count
;
}  ;


struct unit2
{  struct unit1 u1
;  unsigned left_n      /* for paired reads, denotes uncompressed length of left read */
;
}  ;


struct unit3
{  struct unit2 u2
;  char* qp             /* pointer to quality */
;
}  ;



int unit_cmp_count
(  const void* a
,  const void* b
)
   {  return
         ((struct unit1*) a)->count < ((struct unit1*) b)->count
      ?  1
      :     ((struct unit1*) a)->count > ((struct unit1*) b)->count
         ?  -1
         :  0
;  }



struct assoc
{  void* ls
;  unsigned usize
;  unsigned capacity
;  int is_full
;  int n_noconflict
;  int n_used
;  int n_gtone             /* for an optimisation during sorting; bit hairy implications */
;
}  ;

static struct assoc *current_hash = NULL;
static long unsigned old_hash_num = 0;


int assoc_init
(  struct assoc* a
,  unsigned usize          /* unit size */
,  unsigned hsize          /* hash size */
)
   {  a->capacity =  hsize
   ;  a->usize    =  usize
   ;  if (!(a->ls =  calloc(hsize, usize)))
      X_FINISHED("no memory available")
   ;  a->is_full = 0
   ;  a->n_noconflict = 0
   ;  a->n_used = 0
   ;  a->n_gtone = 0
   ;  g_memusage += hsize * usize
   ;  return 0
;  }



#define getcharp(ls, i, usize)   ((char*) ls + ((jlu i) * (jlu usize))) 
#define getunit1(ls, i, usize)   ((struct unit1*) (void*) getcharp(ls, i, usize))
#define getunit2(ls, i, usize)   ((struct unit2*) (void*) getcharp(ls, i, usize))
#define getunit3(ls, i, usize)   ((struct unit3*) (void*) getcharp(ls, i, usize))
#define getcount(ls, i, usize)   (getunit1(ls, i, usize)->count)


                           /* a[i].n_used denotes nr of elements (set by caller).
                            * a[i].n_gtone denotes nr of elements with count > 1 (these will
                            * be at the front because elements are ordered by count).

                            * Either can be used; We now use the latter to
                            * prevent excessive shifting down.  This creates
                            * dependencies with caller follow-up code, as
                            * elements with count == 1 need to be output
                            * separately (necessary in case more than one hash
                            * is used). If a[i].n_used is reinstated then
                            * caller output code needs to change.

                            * The code keep chunks sorted after picking off
                            * next largest element from all of the chunks.

                            * To do this it may need to push down a unit that
                            * is relocated to make way for the currently
                            * largest element that is not yet processed.
                            * Finding the offset where to push down this unit
                            * to is done using a bsearch type loop.

                            * The code operates on an abstract data type, so needs to do char*
                            * manipulation.
                           */
void stack_sortall_bycount
(  struct assoc* a
,  int n_assoc
)
   {  int pi = 0, pj = 0            /* P for pivot (and P-funk) */
   ;  double totalshift = 0.0
   ;  unsigned long stackshift = 0, usize = a->usize

   ;  while (pi+1 < n_assoc)
      {  int qi = pi, i
      ;  char scratch[512]          /* should accomodate the largest of the units; union might be slightly neater */
      ;  unsigned movevalue = getcount(a[pi].ls, pj, usize)
      ;  unsigned newmax    = movevalue

      ;  memcpy(scratch, getcharp(a[pi].ls, pj, usize), usize)

                                    /* pick off largest element from chunk qi */
                                    /* arrays kept sorted, so zero element is the largest */
      ;  for (i=pi+1; i<n_assoc;i++)
         {  unsigned u = getcount(a[i].ls, 0, usize)
         ;  if (u > newmax)
               qi = i
            ,  newmax = u
      ;  }

                                    /* keep chunk qi sorted */
         if (qi > pi)
         {  int j = 1
         ;  memcpy(getcharp(a[pi].ls, pj, usize), getcharp(a[qi].ls, 0, usize), usize)

                                    /* if clause is true, at least one such element exists,
                                     * and the bsearch below is garantueed to work.
                                    */
         ;  if (getcount(a[qi].ls, 1, usize) > movevalue)
            {  int l = 1, r = a[qi].n_gtone
            ;  while (l+1 < r)
               {  int m = l + (r-l) / 2
               ;  if (getcount(a[qi].ls, m, usize) > movevalue)
                  l = m
               ;  else
                  r = m
            ;  }
               j = l+1
         ;  }

         ;  if (j>1)
            memmove(getcharp(a[qi].ls, 0, usize), getcharp(a[qi].ls, 1, usize), (j-1) * usize)

         ;  memcpy(getcharp(a[qi].ls, j-1, usize), scratch, usize)
         ;  stackshift += j-1
      ;  }

         pj++
      ;  if (pj >= a[pi].n_gtone)
         {  pi++
         ;  pj = 0
         ;  totalshift += stackshift
         ;  if (g_verbose)
            argh("hash-merge", "stack shift %lu (average %.1f) new leader %u", stackshift, stackshift * 1.0 / a[pi].n_used,  getcount(a[pi].ls, 0, usize))
         ;  stackshift = 0
      ;  }
      }
   }


            /* Use two pointers moving in different directions.
             * The first (moving forward) stops once it has a singleton
             * -- it needs to move towards the back
             * The second (moving backward) stops once it has a multiplon
             * -- it needs to move forward.
             * At that point the two are swapped and the process continues.
            */
void stack_separate_singletons
(  struct assoc* a
,  int n_assoc
)
   {  unsigned long usize = a->usize
   ;  unsigned i

   ;  for (i=0; i < n_assoc; i++)
      {  struct assoc * b = a+i
      ;  unsigned x = 0, y = b->n_used-1
      ;  char scratch[512]

;if(0)fprintf(stderr, "have used %d gtone %d\n", (int) b->n_used, b->n_gtone)

            /* Need this; otherwise y will run down to zero and wrap back round */
      ;  if (!b->n_gtone)
         continue

            /* In the loop below we may do the same if clause comparison multiple
             * times. The simplicity of the loop is appealing though, and
             * caching hopefully helps.
            */
      ;  while (x <= y)
         {  if (getcount(b->ls, x, usize) > 1)
            x++
         ;  else if (getcount(b->ls, y, usize) == 1)
            y--
         ;  else
            {  memcpy(scratch, getcharp(b->ls, x, usize), usize)     /* count(x) == 1, count(y) > 1 */
            ;  memcpy(getcharp(b->ls, x, usize), getcharp(b->ls, y, usize), usize)
            ;  memcpy(getcharp(b->ls, y, usize), scratch, usize)
            ;  x++
            ;  y--
         ;  }
         }
         if (x != b->n_gtone)
         argh("_____ singleton code", "wobbly boundary in hash %u (boundary=%u gtone=%u) -- file bug please", ju i, ju x, ju b->n_gtone)
   ;  }
   }


#if TRACING_ON
static unsigned long g_cfl = 0;
static double g_hist[256] = { 0.0 };
#endif

void* hash_search
(  struct assoc* a
,  char* s
,  unsigned length
,  unsigned* np_try
)
#  define LINEAR_PROBE 13337
   {  char* h        =  a->ls
   ;  unsigned hsize =  a->capacity
   ;  unsigned usize =  a->usize
   ;  unsigned pos   =  (hsize-1) & hash_sequence(s, length)
   ;  unsigned n_try =  0
   ;  struct unit1* u=  getunit1(h, pos, usize)

#if TRACING_ON
;  g_hist[(int) (pos * 255.99 / hsize)]++
#endif

   ;  while (u->nt && (u->length != length || memcmp(u->nt, s, length)))
      {  pos = (hsize-1) & (pos+LINEAR_PROBE)
      ;  u = getunit1(h, pos, usize)
      ;  n_try++
   ;  }
      if (np_try)
      np_try[0] = n_try
#if TRACING_ON
;g_cfl += n_try
#endif
   ;  return u
;  }


#define HCL_NHASH 50


#if 0
#define MAXFIELDSIZE ((1 << 10)-1)
#else
#define MAXFIELDSIZE ((1 << 15)-1)
#endif

struct record
{  char seq       [MAXFIELDSIZE+1]
;  char q         [MAXFIELDSIZE+1]
;  char discard   [MAXFIELDSIZE+1]
;  char id        [MAXFIELDSIZE+1]

;  unsigned seq_n
;  unsigned q_n
;  unsigned discard_n
;  unsigned id_n

;  unsigned bytesize
;  unsigned count
;  unsigned ID                /* e.g. reaper-created */
;  unsigned roffset
;  unsigned long roffset_counted    /* account for alread-present counts, if present */
;  unsigned error_quality     /* length(q) != length(seq) */
;  unsigned discarded_alien
;  unsigned discarded_unmatched
;  unsigned discarded_length
;  unsigned discarded_trint
;  unsigned discarded_byother
;  unsigned n_paired   /* paired end */
;  unsigned skip
;
}  ;



unsigned int stack_insert
(  struct assoc* a0
,  char* s
,  unsigned length
,  struct record * r1
,  unsigned left_n
,  void** uptr
)
   {  int i
   ;  unsigned usize = a0->usize
   ;  unsigned count = r1->count
   ;  unsigned n_try = 0
   ;  for (i=0;i<HCL_NHASH;i++)
      {  struct unit1* u
;if(0)fprintf(stderr, "searching stack %d for %.*s\n", (int) i, (int) length, s)
      ;  u = hash_search(a0+i, s, length, &n_try)
      ;  uptr[0] = u
;if(0)fprintf(stderr, "[%.*s] %d\n", (int) u->length, u->nt, (int) u->count)
      ;  if (u->nt)
         {  u->count += count

         ;  if (usize > sizeof u[0] && ((struct unit2*) u)->left_n != left_n)
            {  argh
               (  "hiccups"
               , "ambiguous hash [%s] offsets %lu and %lu"
               ,  r1->seq
               ,  jlu ((struct unit2*) u)->left_n
               ,  jlu left_n
               )
            ;  n_bogon++
            ;  g_nt_left_bogon_diff += (int) ((struct unit2*) u)->left_n - (int) left_n
         ;  }
            return 'f'
      ;  }
         else if (!a0[i].is_full)
         {  unsigned captest = a0[i].n_used + (g_stingy ? a0[i].n_used >> g_stingy : a0[i].n_used << 1)
         ;  u->nt       = s
         ;  u->length   = length
         ;  u->count    = count
         ;  if (usize > sizeof u[0])         /* crude test */
            ((struct unit2*) u)->left_n = left_n
         ;  a0[i].n_used++
         ;  if (!n_try)
            a0[i].n_noconflict++
         ;  if (captest >= a0[i].capacity)
            {  a0[i].is_full = 1
            ;  if (i+1<HCL_NHASH)
               {  X_CASCADE(assoc_init(a0+i+1, a0[i].usize, a0[i].capacity))
               ;  current_hash = a0+i+1
               ;  old_hash_num += a0[i].n_used
               ;  if (g_verbose)
                  arrr
                  (  "hash %d created, previous stack used/noconflict/size %u/%u/%u"
                  ,  (int) (i+1)
                  ,  (int) a0[i].n_used
                  ,  (int) a0[i].n_noconflict
                  ,  (int) a0[i].capacity
                  )
            ;  }
            }
            return 'i'
      ;  }
      }
   ;  return 0
;  }



int move_sequences
(  struct record* r1
,  struct record* r2
,  int delta
)
   {  int nospace = SEP_PAIR_WIDTH + r1->seq_n + r2->seq_n > MAXFIELDSIZE + delta
                                 /* below should garantuee we can handle q in the same ways seq;
                                  * any sequence insert check that seq passed successfully automatically
                                  * makes q pass it as well.
                                 */
   ;  int carry_quality = r1->q_n == r1->seq_n
   ;  if (nospace)    /* this loses (does not utilise all) space in the paired case */
      return 1

;if(0)fprintf(stderr, "in s1 %d s2 %d   q1 %d q2 %d\n", (int) r1->seq_n, (int) r2->seq_n, (int) r1->q_n, (int) r2->q_n)
   ;  if (delta)
      {  memmove(r1->seq, r1->seq+delta, r1->seq_n - delta)
      ;  r1->seq_n -= delta
      ;  if (carry_quality)
         {  memmove(r1->q, r1->q+delta, r1->q_n - delta)
         ;  r1->q_n -= delta
      ;  }
      }

      if (g_have_paired)
      {  memcpy(r1->seq+r1->seq_n + SEP_PAIR_WIDTH, r2->seq + delta, r2->seq_n - delta)
      ;  if (SEP_PAIR_WIDTH)
         r1->seq[r1->seq_n] = 'N'
      ;  r1->seq_n += SEP_PAIR_WIDTH + (r2->seq_n - delta)

      ;  if (carry_quality)
         {  memcpy(    r1->q+r1->q_n + SEP_PAIR_WIDTH,   r2->q + delta,   r2->q_n - delta)
         ;  if (SEP_PAIR_WIDTH)
            r1->q[r1->q_n] = ' '
         ;  r1->q_n += SEP_PAIR_WIDTH + (r2->q_n - delta)
      ;  }
      }
      r1->q[r1->q_n] = '\0'
   ;  r1->seq[r1->seq_n] = '\0'
;if(0)fprintf(stderr, "out s1 %d    q1 %d [%s]\n", (int) r1->seq_n, (int) r1->q_n, r1->q)
   ;  return 0
;  }



#define TALLY_ERROR 3
#define TALLY_NOMEM 2
#define TALLY_DONE  1
#define TALLY_OK    0

#define izblank(c) ((unsigned char) (c) == ' ' || (unsigned char) (c) == '\t')


int tally_cpytofield
(  char* dest
,  char* src
,  int   n
)
   {  if (n > MAXFIELDSIZE)
         n = MAXFIELDSIZE
      ,  g_tally_fieldtosmall++
   ;  memcpy(dest, src, n)
   ;  dest[n] = '\0'
   ;  return n
;  }


int read_record3
(  ZFILE  ip
,  struct file_buffer* fb
,  const char* format
,  struct record* rec
)
   {  const char* fmtp = format, *fmtz = format + strlen(format)
   ;  unsigned rlstat = 0
#if 0
#define LINE_LIMIT ((1<<13)-1)
#else
#define LINE_LIMIT ((1<<15)-1)
#endif
   ;  char buf[LINE_LIMIT+1]
   ;  unsigned n_received = 0
   ;  unsigned n_lines = 0
   ;  unsigned truncated = 0
   ;  char* bufp = NULL, *curp = NULL, *bufz = NULL
   ;  int esc = 0

   ;  rec->bytesize  =  0
   ;  rec->seq_n     =  0
   ;  rec->discard_n =  0
   ;  rec->q_n       =  0
   ;  rec->id_n      =  0
   ;  rec->count     =  1
   ;  rec->ID        =  0
   ;  rec->seq[0]    =  '\0'
   ;  rec->discard[0]=  '\0'
   ;  rec->q[0]      =  '\0'
   ;  rec->id[0]     =  '\0'

   ;  while (fmtp < fmtz)
      {  if (!bufp)
         {  rlstat |= kraken_hookline(ip, fb, buf, LINE_LIMIT, &n_received, &truncated)
         ;  if (truncated || rlstat & (RL_DONE | RL_ERROR))
            break
         ;  bufp = buf
         ;  bufz = buf + n_received
         ;  esc  = 0
         ;  rec->bytesize += n_received
         ;  n_lines++
      ;  }
         curp = bufp

      ;  if (esc)
         {  switch((unsigned char) fmtp[0])
            {  case '#':
                  fmtp++
               ;  bufp = NULL
               ;  continue
               ;
            case 'n':
                  if (bufp != bufz)
                  goto DONE
               ;  fmtp++
               ;  bufp = NULL
               ;  continue
               ;
            case '%':
                  if (bufp[0] != '%') goto DONE
               ;  bufp++
               ;  break
               ;
            case 't':
                  if (bufp[0] != '\t') goto DONE
               ;  bufp++
               ;  break
               ;
            case 's':
                  if (bufp[0] != ' ') goto DONE
               ;  bufp++
               ;  break
               ;
            case '.':
                  bufp++
               ;  break
               ;
            case 'b':
                  while (bufp < bufz && izblank(bufp[0]))
                  bufp++
               ;  break
               ;
            case 'X': case 'C':
                  {  int n_scanned = 0
                  ;  if (sscanf(bufp, "%u%n", &rec->count, &n_scanned) < 1)
                     goto DONE
                  ;  bufp += n_scanned
               ;  }
                  break
               ;
            case 'J':
                  {  int n_scanned = 0
                  ;  if (sscanf(bufp, "%u%n", &rec->ID, &n_scanned) < 1)
                     goto DONE
                  ;  bufp += n_scanned
               ;  }
                  break
               ;
            case 'G':
                  while (bufp < bufz && !izblank((unsigned char) bufp[0]))
                  bufp++
               ;  break
               ;
            case 'F':
                  while (bufp < bufz && '\t' != (unsigned char) bufp[0])
                  bufp++
               ;  break
               ;
            case 'H':
                  while (bufp < bufz && (unsigned char) bufp[0] != (unsigned char) fmtp[1])
                  bufp++
               ;  if (bufp == bufz)
                  goto DONE
               ;  fmtp++
               ;  bufp++
               ;  break
               ;
            case 'Q':
                  while (bufp < bufz && !izblank((unsigned char) bufp[0]))    /* need izblank or sth similar to stop on newline */
                  bufp++
               ;  rec->q_n = tally_cpytofield(rec->q, curp, bufp - curp)
               ;  break
               ;
            case 'I':
                  while (bufp < bufz && !izblank((unsigned char) bufp[0]))
                  bufp++
               ;  rec->id_n = tally_cpytofield(rec->id, curp, bufp - curp)
               ;  break
               ;
            case 'R':
                  while (bufp < bufz && !izblank((unsigned char) bufp[0]))
                     bufp[0] = toupper((unsigned char) bufp[0])
                  ,  bufp++
               ;  rec->seq_n = tally_cpytofield(rec->seq, curp, bufp - curp)
               ;  break
               ;
            default:
                  goto DONE
         ;  }
            esc = 0
      ;  }
         else
         {  if (fmtp[0] == '%')
            esc = 1
         ;  else if (bufp[0] != fmtp[0])
            goto DONE
         ;  else
            bufp++
      ;  }
         fmtp++
   ;  }

      DONE :

      if (truncated)
         argh("reaper", "line too long")
      ,  rlstat |= RL_ERROR

   ;  if (n_lines && (fmtp[0] || (bufp && bufp - buf != n_received)))
      {  argh("tally", "parse error at line %lu (remaining format string [%s], buffer [%s])", jlu n_lines, fmtp, bufp)
      ;  return TALLY_ERROR
   ;  }

      if (rlstat & RL_DONE)
      return TALLY_DONE

   ;  if (rlstat & RL_ERROR)
      return TALLY_ERROR

   ;  if (rlstat & RL_NOMEM)           /* fixme not checked yet */
      return TALLY_NOMEM

   ;  rec->roffset++
   ;  rec->roffset_counted += rec->count
   ;  return TALLY_OK
;  }



void writerecord
(  struct record* rec
,  const char* format
,  void* fpo
,  int zippit
)
   {  char output[MAXFIELDSIZE+1]
   ;  const char* f = format
   ;  const char* z = format + strlen(format)
   ;  int escape = 0

   ;  while (f < z)
      {  unsigned c = (unsigned char) f[0]
      ;  const char* p = output
      ;  int n = 0
      ;  if (!escape)
         {  if (c == '%')
            escape = 1
         ;  else
               output[0] = c
            ,  n = 1
      ;  }
         else
         {  switch(c)
            {  case 'R': p = rec->seq ;  n = rec->seq_n ;  break
            ;  case 'Q': p = rec->q ;    n = rec->q_n ;  break
            ;  case 'T': { int trint = trintscore(rec->seq, rec->seq_n);
                           n = snprintf(output, MAXFIELDSIZE, "%u", trint); }
                           break
            ;  case 'q':
                         { int a = (unsigned char) f[1], i;
                           for (i=0;i<rec->seq_n && i < MAXFIELDSIZE; i++) output[i] = a;
                           output[i] = '\0';
                           n = i;
                           f++;
                           break;
                         }
            ;  case 'I': p = rec->id ;   n = rec->id_n ;  break
            ;  case 'J': n = snprintf(output, MAXFIELDSIZE, "%u", rec->ID) ;  break
            ;  case 'n':
               case 't':
               case 's':
               case '%':   output[0] = c == 'n' ? '\n' : c == 't' ? '\t' : c == 's' ? ' ' : '%'
                        ;  n = 1
                        ;  break
         ;  }
            escape = 0
      ;  }
         if (n > 0)
         {  
#if WE_USE_ZLIB
         if (zippit)
            gzwrite(fpo, p, n)
         ;  else
#endif
            fwrite(p, n, 1, fpo)
      ;  }
         f++
   ;  }
   }



      /* if delta > 0, read layout is
       *       iiiii***********************                                 in single file mode
       *       iiiii***********************Niiiii***********************    in paired end mode
       *       <-------------------------->                                 left_n
      */
void writeline
(  const char* read
,  unsigned n_read
,  unsigned delta
,  unsigned count
,  unsigned left_n
,  unsigned trint
,  unsigned output_id
,  const char* format
,  void* fpo
,  int zippit
,  unsigned lr             /* left/right bits */
,  const char* q
)
   {  char output[MAXFIELDSIZE+1]
   ;  const char* f = format
   ;  const char* z = format + strlen(format)
   ;  int escape = 0

   ;  if (delta > left_n || (lr == 2 && (delta > n_read -left_n)))
      {  argh("tally", "delta error")
      ;  return
   ;  }

#define READ_OFS  (delta)
#define READ_LEN  (n_read - delta)
#define LFT_LEN   (left_n - delta)
#define LFT_OFS   (delta)
#define RGT_LEN   (n_read - SEP_PAIR_WIDTH - left_n - delta)
#define RGT_OFS   (left_n + delta + SEP_PAIR_WIDTH)

      while (f < z)
      {  unsigned c = (unsigned char) f[0]
      ;  const char* p = output
      ;  int n = 0
      ;  if (!escape)
         {  if (c == '%')
            escape = 1
         ;  else
               output[0] = c
            ,  n = 1
      ;  }
         else
         {  switch(c)
            {  case 'R':
                  p = read + READ_OFS
               ;  n = READ_LEN

               ;  if (!lr)
                  g_nt_left_out += n * count
               ;  else if (lr == 1)
                     n = LFT_LEN
                  ,  g_nt_left_out += n * count
               ;  else if (lr == 2)
                     p = read + RGT_OFS
                  ,  n = RGT_LEN
                  ,  g_nt_right_out += n * count
               ;  break

            ;  case 'Q':
                  if (!q)
                     p = "noqual"            /* fixme: fake quality would be better; also tally this error */
                  ,  n = 6
               ;  else
                  {  p = q + READ_OFS
                  ;  n = READ_LEN
                  ;  if (lr == 1)
                     n = LFT_LEN
                  ;  else if (lr == 2)
                        p = q + RGT_OFS
                     ,  n = RGT_LEN          /* noteme */
               ;  }
                  break

            ;  case 'A': p = read + LFT_OFS; n = LFT_LEN; break;
            ;  case 'B': p = read + RGT_OFS; n = RGT_LEN; break;

            ;  case 'U': p = q + LFT_OFS; n = LFT_LEN; break;
            ;  case 'V': p = q + RGT_OFS; n = RGT_LEN; break;

            ;  case 'T': n = snprintf(output, MAXFIELDSIZE, "%u", trint) ;  break
            ;  case 'X' : case 'C': n = snprintf(output, MAXFIELDSIZE, "%u", count) ;  break
            ;  case 'K': n = snprintf(output, MAXFIELDSIZE, "%u", left_n) ;  break
            ;  case 'L': n = snprintf(
                                 output, MAXFIELDSIZE, "%u",
                                 (unsigned) (!lr ? READ_LEN : lr == 1 ? LFT_LEN : RGT_LEN)
                             )
                  ;  break
            ;  case 'I': n = snprintf(output, MAXFIELDSIZE, "%u", output_id) ;  break
            ;  case 'n':
               case 't':
               case 's':
               case '%':   output[0] = c == 'n' ? '\n' : c == 't' ? '\t' : c == 's' ? ' ' : '%'
                        ;  n = 1
                        ;  break
            ;  case 'q':
                           {  int a = (unsigned char) f[1], i;
                              n  =
                                    lr == 1
                                 ?  LFT_LEN
                                 :     lr == 2
                                    ?  RGT_LEN
                                    :  READ_LEN

                           ;  for (i=0;i<n && i < MAXFIELDSIZE; i++)
                              output[i] = a
                           ;  output[i] = '\0'
                           ;  n = i
                           ;  f++
                           ;  break
                        ;  }
         ;  }
            escape = 0
      ;  }
         if (n > 0)
         {  
#if WE_USE_ZLIB
            if (zippit)
            gzwrite(fpo, p, n)
         ;  else
#endif
            fwrite(p, n, 1, fpo)
      ;  }
         f++
   ;  }
   }


int gauge_mem_pams
(  ZFILE input
,  const char* format_in
,  unsigned usize
,  unsigned long fsize
,  double gzip_factor
,  int with_quality
,  int logds_delta
,  int loghs_delta
,  unsigned* logdsp
,  unsigned* loghsp
)
   {  int ntest = 0
   ;  struct record rc = { { 0 }, { 0 }, { 0 }, { 0 }, 0 }
   ;  struct file_buffer fb = { { 0 }, 0 }
   ;  unsigned sumrecsize = 0
   ;  unsigned sumreadsize = 0
   ;  unsigned truncated = 0
   ;  unsigned st = 0

   ;  while (!(st = read_record3(input, &fb, format_in, &rc)))
      {  sumrecsize +=  rc.bytesize
      ;  sumreadsize += rc.seq_n
      ;  if (++ntest == g_ntest)                      /* magic parameter */
         break
   ;  }

      kraken_readline(NULL, NULL, 0, NULL, &truncated)           /* resets static buffers */
   ;  kraken_hookline(NULL, &fb, NULL, 0, NULL, &truncated)      /* resets buffers */

   ;  if (sumreadsize && ntest)
      {  double est_readsize = 0.0, est_readnum = 0.0, est_qsize = 0.0
      ;  double redundancy_reduction = 0.8
      ;  double my_compression = 3.4
      ;  double hash_expansion = 1.6               /* 1.5 in reality */
      ;  double paired_expansion = g_have_paired ? 2.0 : 1.0
      ;  double logdsf, loghsf
      ;  int logds = 0, loghs = 0, memmegs

      ;  est_readsize = fsize * gzip_factor * paired_expansion * redundancy_reduction * (sumreadsize * 1.0 / sumrecsize)
                                                                  /* ^ does not take into account 2-bit compression */
      ;  if (with_quality)
         est_qsize = est_readsize + est_readnum                   /* bit pedantic; est_readnum is the extra '\0' terminator */

      ;  est_readnum  = est_readsize / (sumreadsize * 1.0 / ntest)
      ;  logdsf = log(est_qsize + est_readsize/my_compression) / log(2.0)
      ;  loghsf = log(hash_expansion * est_readnum) / log(2.0)
      ;  logds = logdsf - 2.5                                     /* round down; chunks are cheap       */
      ;  loghs = loghsf - 1.5                                     /* round down; 2 hashes is acceptable */
      ;  memmegs = ((1UL << logds) + (1UL << loghs) * usize) >> 20

      ;  if (g_verbose)
         {  argh("peekaboo", "average read length: %.0f", 0.5 + sumreadsize * 1.0 / ntest)
         ;  argh("peekaboo", "read/record ratio: %.2f", sumreadsize * 1.0 / sumrecsize)
         ;  argh("peekaboo", "estimated read count: %.1fM", est_readnum * 1.0 / (redundancy_reduction * (1L << 20)))
         ;  argh("peekaboo", "nonredundant read count: %.1fM", est_readnum * 1.0 / (1L << 20))
         ;  argh("peekaboo", "data size log: %.1f (rounded %d)", logdsf, logds)
         ;  argh("peekaboo", "hash size log: %.1f (rounded %d)", loghsf, loghs)
         ;  argh("peekaboo", "associated memory usage: %dM (unit size %d)", memmegs, (int) usize)
         ;  argh("peekaboo", "expect %dM, use %dM to be very safe", (int) (1.0 + memmegs * 1.1), (int) (1.0 + memmegs * 2.2)) 
      ;  }

         if (!logdsp[0])
         {  logds += logds_delta
         ;  if (logds > DATA_SIZE_LOG_MAX)
            logds = DATA_SIZE_LOG_MAX
         ;  if (logds < DATA_SIZE_LOG_MIN)
            logds = DATA_SIZE_LOG_MIN
         ;  logdsp[0] = logds
      ;  }

         if (!loghsp[0])
         {  loghs += loghs_delta
         ;  if (loghs > HASH_SIZE_LOG_MAX)
            loghs = HASH_SIZE_LOG_MAX
         ;  if (loghs < HASH_SIZE_LOG_MIN)
            loghs = HASH_SIZE_LOG_MIN
         ;  loghsp[0] = loghs
      ;  }
      }
      return 0
;  }


clock_t clockit
(  clock_t t1
,  double* task_time
,  double* sum_time
)
   {  clock_t t2 = clock()
   ;  task_time[0] = (t2 - t1) * 1.0 / CLOCKS_PER_SEC
   ;  sum_time[0] += task_time[0]
   ;  return t2
;  }


#define TALLY_VERBOSE_CMP  1 << 0
#define TALLY_NO_TALLY     1 << 1
#define TALLY_FORCE_PAIRS  1 << 2
#define TALLY_CLEVER       1 << 3
#define TALLY_PEEK         1 << 4


                                                 /* 27 is 134M-ish */
                                                 /* 26 is  67M-ish */
                                                 /* 25 is  33M-ish */
                                                 /* 24 is  17M-ish */
                                                 /* 23 is   8M-ish */
                                                 /* 22 is   4M-ish */
                                        /* hsize should be multipled by 4 to obtain memory usage */
#define STATUS_CLINE_PROCESSING 22
#define STATUS_DATA_PROCESSING  99


int tally_main
(  int argc
,  const char* argv[]
)
   {  unsigned long storage_size = 0
   ;  unsigned long hsize = 0
   ;  unsigned logds = 0
   ;  unsigned loghs = 0
   ;  int      loghs_delta = 0
   ;  int      logds_delta = 0

   ;  int status = STATUS_CLINE_PROCESSING
   ;  int called_from_R = 0
   ;  int changed_name = 0
   ;  int lower = 0, upper = 0, trintco = 0
   ;  unsigned long len_compressed = 0, len_uncompressed = 0
   ;  unsigned n_pass = 0, n_pass_unique = 0
   ;  unsigned n_chunk = 1
   ;  clock_t t1 = clock()
   ;  double time_iput = 0.0
   ;  double time_oput = 0.0
   ;  double time_process = 0.0
   ;  double time_task = 0.0

   ;  struct assoc lsohash[HCL_NHASH] 
   ;  char* storage, *cp

   ;  unsigned long storage_used = 0
   ;  const char* gzopenmode = ZWMODE

   ;  ZFILE fpo = NULL, input = NULL, input2 = NULL, fpo2 = NULL
   ;  FILE* fpoh = NULL
   ;  const char* g_fnin = "-", *g_fnin2 = NULL
   ;  const char* g_fnout = "out.tally.gz", *g_fnout2 = NULL
   ;  const char* g_fnouth = NULL

   ;  const char* fastqformatin     =  "@%I%#%R%n+%#%Q%n"
   ;  const char* fastqxformatin    =  "@%I%brecno=%J%#%R%n+%#%Q%n"

   ;  const char* fastaformatin     =  ">%I%#%R%n"
   ;  const char* fastaxformatin     =  ">%I%brecno=%J%#%R%n"

   ;  const char* fastqformatout    =  "@" "trn_%I%s%C" "%n%R%n+%n%Q%n"
   ;  const char* fastqformatoutnoq =  "@" "trn_%I%s%C" "%n%R%n+%n%q~%n"
   ;  const char* fastaformatout    =  ">" "trn_%I%s%C" "%n%R%n"

/* no-tally (nt) variants */
   ;  const char* fastqformatoutnt    =  "@" "%I" "%n%R%n+%n%Q%n"
   ;  const char* fastqformatoutntnoq =  "@" "%I" "%n%R%n+%n%q~%n"
   ;  const char* fastaformatoutnt    =  ">" "%I" "%n%R%n"

   ;  const char* g_format_in       =  fastqformatin
   ;  const char* g_format_out      =  fastqformatoutnoq

   ;  char* my_fnout   = NULL
   ;  unsigned limit   = 0
   ;  unsigned limit2  = 0
   ;  int use_compression = 1, do_output = 1, do_sort = 'C'
   ;  unsigned sinsert = 0
   ;  int zippit = 1
   ;  int early = 0
   ;  int with_quality = 0
   ;  unsigned truncated = 0
   ;  unsigned n_pair_error = 0, n_pair_id = 0
   ;  unsigned tally_modes = TALLY_CLEVER
   ;  unsigned long inputfilesize = 0
   ;  double gzip_factor = 3.2
   ;  int n_hash_used = 0

   ;  struct record whoopsie[2] = { { { 0 }, { 0 }, { 0 }, { 0 }, 0 }, { { 0 }, { 0 }, { 0 }, { 0 }, 0 } }
   ;  struct record* r1 = whoopsie+0
   ;  struct record* r2 = whoopsie+1

   ;  struct file_buffer fb1 = { { 0 }, 0 }
   ;  struct file_buffer fb2 = { { 0 }, 0 }

   ;  unsigned id_prev = 0
   ;  unsigned id_prev2 = 0

   ;  unsigned usize = sizeof(struct unit1)

#ifdef WINDOWS_BUILD
   ;  set_argh("reaperlog.txt", 1)
#endif

   ;  kraken_readline(NULL, NULL, 0, NULL, &truncated)            /* resets static buffers */
   ;  kraken_hookline(NULL, &fb1, NULL, 0, NULL, &truncated)      /* resets buffers */
   ;  kraken_hookline(NULL, &fb2, NULL, 0, NULL, &truncated)      /* resets buffers */

   ;  memset(lsohash, 0, sizeof(lsohash))
   ;  current_hash = lsohash

   ;  themap['A'] = 0
   ;  themap['C'] = 1
   ;  themap['G'] = 2
   ;  themap['T'] = 3
   ;  themap['N'] = 4
   ;  themap['X'] = 4
   ;
      arg_switch()

      optarg("-o")   g_fnout  = thearg(); endarg()
      optarg("-p")   g_fnout2 = thearg(); endarg()
      optarg("-i")   g_fnin   = thearg(); endarg()
      optarg("-j")   g_fnin2  = thearg();  g_have_paired = 1; endarg()
      optarg("-si")  sinsert  = atoi(thearg()); early = 1; endarg()
      optarg("-dsi") sinsert  = atoi(thearg()); endarg()       /* early == 0 by default */
      optarg("-do")  limit    = atoi(thearg()); endarg()
      optarg("-gzwrite")  gzopenmode = thearg(); endarg()
      optarg("-sumstat") g_fnouth = thearg(); endarg()

      optarg("-format") g_format_out = thearg(); endarg()
      optarg("-record-format") g_format_in = thearg(); endarg()

      uniarg("--stingier") g_stingy++; endarg()
      optarg("-peek")g_ntest  = atoi(thearg()); endarg()
      uniarg("--fastqx-in") g_format_in  = fastqxformatin; endarg()
      uniarg("--fastax-in") g_format_in  = fastaxformatin; endarg()
      uniarg("--fasta-in")  g_format_in  = fastaformatin;  endarg()
      uniarg("--fasta-out") g_format_out = fastaformatout; endarg()

      optarg("-hsx") int a = atoi(thearg()); if (a >= 2  && a <= 32) loghs = a; endarg()
      optarg("-dsx") int a = atoi(thearg()); if (a >= 2  && a <= 36) logds = a; endarg()
      optarg("-hsd") loghs_delta = atoi(thearg()); endarg()
      optarg("-dsd") logds_delta = atoi(thearg()); endarg()
      optarg("-hs")
         int a = atoi(thearg()); if (a >= HASH_SIZE_LOG_MIN && a <= HASH_SIZE_LOG_MAX) loghs = a;
      endarg()
      optarg("-ds")
         int a = atoi(thearg()); if (a >= DATA_SIZE_LOG_MIN && a <= DATA_SIZE_LOG_MAX) logds = a;
      endarg()

      optarg("-l")   lower = atoi(thearg()); endarg()
      optarg("-zip-factor")
         gzip_factor = atof(thearg());
         if (gzip_factor < 1.0) gzip_factor = 1.0;
         if (gzip_factor > 10.0)gzip_factor = 10.0;
      endarg()
      optarg("-u")   upper = atoi(thearg()); endarg()
      optarg("-tri")
         trintco = atoi(thearg());
         if (trintco < 0 || trintco > 100) {
            argh("tally", "-trint option should take argument in [0..100]");
            trintco = 100;
         }
      endarg()

      optarg("-give") limit2 = atoi(thearg()); endarg()

      uniarg("--nozip") zippit = 0; endarg()
      uniarg("--with-quality") with_quality = 1; endarg()
      uniarg("--peek")  tally_modes |= TALLY_PEEK; endarg()
      uniarg("--no-auto")  tally_modes |= TALLY_CLEVER; tally_modes ^= TALLY_CLEVER; endarg()
      uniarg("--cx") use_compression = 0; endarg()
      uniarg("--noput") do_output = 0; endarg()
      uniarg("--plenty") g_stingy = 0; endarg()
      uniarg("--unsorted") do_sort = 0; endarg()
      uniarg("--pair-by-offset") tally_modes |= TALLY_FORCE_PAIRS; endarg()
      uniarg("--no-tally") tally_modes |= TALLY_NO_TALLY; endarg()
      uniarg("--dirty-mem") g_gauge_memory = 1; endarg()
      optarg("-v")
         if (strstr(thearg(), "cmp")) tally_modes |= TALLY_VERBOSE_CMP;
         if (strstr(thearg(), "work")) g_verbose = 1;
      endarg()
      uniarg("--R") called_from_R = 1;    endarg()

      uniarg("-h")
puts("-i <fname>        input stream (gzipped file allowed) (default STDIN)");
puts("-o <fname>        (gzipped!) output stream (default out.tally.gz)");
fprintf(stdout, "--fasta-in        expect FASTA format (same as -record-format '%s')\n", fastaformatin);
fprintf(stdout, "--fasta-out       write FASTA format (same as -format '%s')\n", fastaformatout);
puts("");

puts("--with-quality    pass quality scores along, collate by taking per-base max");
puts("--no-auto         do not peek in input file and set memory parameters automatically");
puts("--peek            peek in input file and output estimated memory parameters");
puts("-zip-factor <num> assume compression factor <num> (use 1.0 for uncompressed files");
puts("-l <int>          require read length >= <int>");
puts("-u <int>          require read length <= <int>");
puts("-tri <int>        required tri-nucleotide score <= <int>");
puts("-si <int>         strip <int> bases from start of read before uniquifying");
puts("-dsi <int>        as -si but after uniquifying (degenerate sequence insert)");
puts("-sumstat  <fname> output file with counts of discarded categories");
puts("");

puts("PAIRED END functionality");
puts("-j <fname>        second paired end input stream");
puts("               -> (requires -record-format with %J or --fastqx-in or --fastax-in)");
puts("-p <fname>        (gzipped!) second output stream for second paired end (cf -j)");
fprintf(stdout, "--fastax-in       expect reaper --fastax-out format (same as -record-format '%s')\n", fastaxformatin);
fprintf(stdout, "--fastqx-in       expect reaper --fastqx-out format (same as -record-format '%s')\n", fastqxformatin);
fprintf(stdout, "               -> these two options are for re-pairing individually processed paired-end files\n");
puts("");

puts("-hsd [1,-1,2,-2]  increase or decrease hash size relative to default");
puts("-dsd [1,-1,2,-2]  increase or decrease data size relative to default");
puts("-hs k             k in 14..32 specifies hash size 2 ** k");
puts("-ds k             k in 14..31 specifies storage size 2 ** k");
puts("--unsorted        do not sort output sequences");
puts("--cx              do not compress sequence (unit testing)");
puts("   NOTE with --cx output will not be sorted in the same way");
puts("   because hash values change. To compare, sort outputs");
puts("--noput           do not output uniquified sequences");
puts("-v <work|cmp>     turn on verbosity settings");
puts("      cmp         with cmp paired end identifier mismatches will be reported");
puts("-record-format    specify input format");
puts("     The same syntax as documented under reaper --record-format,");
puts("     Additionally %J is accepted and assumes a numerical ID that");
puts("     will be strictly increasing.");
puts("     If -j is used this ID is required and will be used to match reads.");
puts("     This can be used in conjunction with the reaper %J format directive.");
puts("-format <format>  output format specification, syntax below");
puts("     %R  read");
puts("     %L  length");
puts("     %C  number of occurrences");
puts("     %T  trinucleotide score");
puts("     %I  read identifier - numerical identifier constructed on output");
puts("         CAVEAT read identifier could differ between runs depending on options");
puts("         CAVEAT read identifier is not tied to the read sequence");
puts("     %t  tab");
puts("     %s  tab");
puts("     %n  newline");
puts("     %%  percentage character");
puts("  Anything else is copied verbatim");
puts("  Paired-end reads (activated with -j option) can be output with -o and");
puts("     %R  reads concatenated from -i and -j input files, matched by %J");
puts("     %Q  quality concatenated from -i and -j input files");
puts("     %A  matched read from -i input file");
puts("     %B  matched read from -j input file");
puts("     %U  quality from -i input file");
puts("     %V  quality from -j input file");
puts("     %K  length of first paired end");
puts("  If -p option is used, %R becomes context aware and refers to %A (-o) or %B (-p)");
puts("  If -p option is used, %Q becomes context aware and refers to %A (-o) or %B (-p)");
puts("--no-tally        reads are output as they are processed (use --with-quality to retain quality)");
puts("  This can be useful for matching up paired-end files with missing reads.");
puts("  It works only if record offset information was preserved and is read back in using %J");
puts("  The supported output directives (besides %n %s %t and %%) are these:");
puts("     %R  read");
puts("     %Q  quality");
puts("     %I  identifier");
puts("     %J  output offset");
puts("  Filter options such as -tri and -si are NOT active");
puts("--pair-by-offset  assume the -i and -j input files match record-by-record");
puts("     With this option the %J directive is not needed");
exit(0);
      endarg()

      uniarg("--version")
         fprintf(stdout, "Tally version: %s\n", tally_tag);
         exit(0);
      endarg()

      failarg()

      arg_done()

   ;  if (g_format_out == fastqformatoutnoq && with_quality)
      g_format_out = fastqformatout

   ;  if (tally_modes & TALLY_NO_TALLY)
      {  if (g_format_out == fastqformatoutnoq)
         g_format_out = fastqformatoutntnoq
      ;  else if (g_format_out == fastqformatout)
         g_format_out = fastqformatoutnt
      ;  else if (g_format_out == fastaformatout)
         g_format_out = fastaformatoutnt
   ;  }

      if (called_from_R)
      argh("tally", "R is calling")

   ;  if (g_have_paired)
      usize = sizeof(struct unit2)
   ;  if (with_quality)
      usize = sizeof(struct unit3)

   ;  if
      (  g_fnin2
      && !(tally_modes & TALLY_FORCE_PAIRS)
      && g_format_in != fastqxformatin
      && (!g_format_in || !strstr(g_format_in, "%J"))
      )
      {  argh("format error", "paired end reads require one of these")
      ;  fprintf(stderr, "  1)  --fastqx-in applied to reaper output produced with reaper --fastqx-out\n")
      ;  fprintf(stderr, "  2)  A -record-format argument using the %%J tag (to retrieve record offset)\n")
      ;  fprintf(stderr, "  3)  The --pair-by-offset option (files have line-by-line matching records)\n")
      ;  X_ERR_JUMP(DONE, "Please use on of the above when using the -j option")
   ;  }
      if (!(input = myfopen(g_fnin, "r", 1)))
      X_ERR_JUMP(DONE, "bailing out (cannot read -i file)")
   ;  if (g_fnin2 && !(input2 = myfopen(g_fnin2, "r", 1)))
      X_ERR_JUMP(DONE, "bailing out (cannot read -j file)")


   ;  if (tally_modes & (TALLY_CLEVER | TALLY_PEEK))
      {  struct stat st = { 0 }
      ;  int stok = 0 == stat(g_fnin, &st)
#if WE_USE_ZLIB
      ;  if (!strcmp(g_fnin, "-") || !stok || gzseek(input, 0L, SEEK_CUR))
#else
      ;  if (!strcmp(g_fnin, "-") || !stok || fseek(input, 0L, SEEK_CUR))
#endif
         argh("tally", "no peek, never clever (memory parameters could not be estimated)")
      ;  else
         {  inputfilesize = st.st_size
         ;  X_CASCADE(
               gauge_mem_pams
               (  input
               ,  g_format_in
               ,  usize
               ,  inputfilesize
               ,  gzip_factor
               ,  with_quality
               ,  logds_delta
               ,  loghs_delta
               ,  &logds
               ,  &loghs
               )  )
         ;  if (tally_modes & TALLY_PEEK)
            argh("tally", "recommended settings are:    -hs %d -ds %d", loghs, logds)

#if WE_USE_ZLIB
         ;  if (!(tally_modes & TALLY_PEEK) && 0 > gzseek(input, 0L, SEEK_SET))
#else
         ;  if (!(tally_modes & TALLY_PEEK) && 0 > fseek(input, 0L, SEEK_SET))
#endif
            X_ERR_JUMP(DONE, "rewind to start failed")
      ;  }
         if (tally_modes & TALLY_PEEK)
         exit(0)
   ;  }

      if (!logds)
      logds = DATA_SIZE_LOG_DEFAULT
   ;  if (!loghs)
      loghs = HASH_SIZE_LOG_DEFAULT

   ;  storage_size = 1L << logds
   ;  hsize        = 1L << loghs

   ;  argh
      (  "tally"
      ,  "data log size %d (%.2fG) hash log size %d (%.2fG) unit size %d"
      ,  (int) logds
      ,  storage_size * 1.0 / (1 << 30)
      ,  (int) loghs
      ,  hsize * usize * 1.0 / (1 << 30)
      ,  (int) usize
      )
   ;

      {  int len = strlen(g_fnout)
      ;  if (0 && called_from_R && zippit && (len < 4 || strcmp(g_fnout + len -3, ".gz")))
            my_fnout = stringle("%s.gz", g_fnout)
         ,  changed_name = 1
      ;  else
         my_fnout = stringle("%s", g_fnout)
      ;  if (do_output && !(fpo = myfopen(my_fnout, gzopenmode, zippit)))
         X_ERR_JUMP(DONE, "cannot open output file")
      ;  if (do_output && g_fnout2 && !(fpo2 = myfopen(g_fnout2, gzopenmode, zippit)))
         X_ERR_JUMP(DONE, "cannot open second output file")
   ;  }

      if (g_fnouth && !(fpoh = myfopen(g_fnouth, "w", 0)))
      X_ERR_JUMP(DONE, "bailing out")

   ;  if (assoc_init(lsohash, usize, hsize))
      goto DONE
   ;  cp = storage = malloc(storage_size)
   ;  g_memusage += storage_size

   ;  status = STATUS_DATA_PROCESSING

                                                      /* funcify? */
   ;  {  int left_n = 0
      ;  int lus   = lower || upper || sinsert
      ;  unsigned active = g_have_paired ? 3 : 1          /* bit field */
      ;  unsigned force  = tally_modes & TALLY_FORCE_PAIRS
      ;  unsigned long sumreadsize = 0, sumrecsize = 0

      ;  while (1)
         {  unsigned n_compress = 0, hstatus, n_receive = 0
         ;  unsigned readme = 1        /* bits: 1 == read from input, 2 == read from input2 */

         ;  if (g_have_paired)
            {  readme = (r1->ID == r2->ID || force ? 3 : r1->ID < r2->ID ? 1 : 2) & active
            ;  switch(readme)
               {  case 1: r1->discarded_unmatched += r1->count; break
               ;  case 2: r2->discarded_unmatched += r2->count; break
            ;  }
            }

            if (readme & 1)
            {  unsigned tallystat = read_record3(input, &fb1, g_format_in, r1)
            ;  if (tallystat == TALLY_DONE)
               {  active ^= 1
               ;  r1->ID = UINT_MAX       /* so that the other stream will be instructed to read */
            ;  }
               else if (tallystat == TALLY_NOMEM || tallystat == TALLY_ERROR)
               goto DONE
            ;  else
                  sumreadsize += r1->seq_n
               ,  sumrecsize += r1->bytesize

            ;  if (with_quality && r1->q_n != r1->seq_n)
               {  if (r1->q_n < r1->seq_n)
                  memset(r1->q + r1->q_n, '~', r1->seq_n - r1->q_n)
               ;  r1->q_n = r1->seq_n
               ;  r1->q[r1->q_n] = '\0'
               ;  r1->error_quality += r1->count
            ;  }

               if (r1->ID && r1->ID <= id_prev)
               X_ERR_JUMP(DONE, "i-file %s record IDs not increasing (%u follows %u)\n", g_fnin, r1->ID, id_prev)
            ;  id_prev = r1->ID
         ;  }
            if (readme & 2)
            {  unsigned tallystat = read_record3(input2, &fb2, g_format_in, r2)
            ;  if (tallystat == TALLY_DONE)
               {  active ^= 2
               ;  r2->ID = UINT_MAX       /* so that the other stream will be instructed to read */
            ;  }
               else if (tallystat == TALLY_NOMEM || tallystat == TALLY_ERROR)
               goto DONE

            ;  if (with_quality && r2->q_n != r2->seq_n)
               {  if (r2->q_n < r2->seq_n)
                  memset(r2->q + r2->q_n, '~', r2->seq_n - r2->q_n)
               ;  r2->q_n = r2->seq_n
               ;  r2->q[r2->q_n] = '\0'
               ;  r2->error_quality += r2->count
            ;  }

               if (r2->ID && r2->ID <= id_prev2)
               X_ERR_JUMP(DONE, "j-file %s record IDs not increasing (%u follows %u)\n", g_fnin2, r2->ID, id_prev2)
            ;  id_prev2 = r2->ID
         ;  }

            if (!active)
            break

         ;  if (g_have_paired)
            {  if (r1->ID == r2->ID)
               {  r1->n_paired++
               ;  if (r1->id_n || r2->id_n)
                  {  if (strcmp(r1->id, r2->id))
                     {  n_pair_error++
                     ;  if (tally_modes & TALLY_VERBOSE_CMP)
                        fprintf(stderr, "cmp %s %s\n", r1->id, r2->id)
                  ;  }
                     n_pair_id++
               ;  }
               }
               else
               continue
         ;  }

            {  struct record *rp

            ;  for (rp = whoopsie; rp < whoopsie+2; rp++)
               {  int n_dna = 0, n_other = 0, i
               ;  rp->skip = 0
               ;  for (i=0;i<rp->seq_n;i++)
                  {  if (themap[(unsigned char) rp->seq[i]] < 5)
                     n_dna++
                  ;  else
                     n_other++
               ;  }
                  if (n_dna < n_other)
                  {  rp->skip++
                  ;  rp->discarded_alien += rp->count
               ;  }

                  if (!rp->skip && lus)
                  {  if
                     (  sinsert > rp->seq_n
                     || (lower && rp->seq_n < lower + early * sinsert)
                     || (upper && rp->seq_n > upper + early * sinsert)
                     )
                     {  rp->skip++
                     ;  rp->discarded_length += rp->count
                  ;  }
                  }

                                          /* below depends on bound check above */
                  if (!rp->skip && trintco)
                  {  int trint = trintscore(rp->seq + sinsert, rp->seq_n - sinsert)
                  ;  if (trint >= trintco)
                     {  rp->discarded_trint += rp->count
                     ;  rp->skip++
                  ;  }
                  }

                  if (!g_have_paired)
                  break
            ;  }
               if (g_have_paired)
               {  if (r1->skip && !r2->skip)
                  r2->discarded_byother += r2->count
               ;  else if (!r1->skip && r2->skip)
                  r1->discarded_byother += r1->count
            ;  }

               if (r1->skip + r2->skip)
               {  continue
            ;  }
            }

                                                   /* fixme: this does not do sequence insert */
            if (tally_modes & TALLY_NO_TALLY)
            {  if (do_output)
               {  writerecord(r1, g_format_out, fpo, zippit)
               ;  g_nt_left_out += r1->seq_n * r1->count
            ;  }

               n_pass_unique++
            ;  n_pass += r1->count
            ;  if (do_output && g_have_paired)
               {  writerecord(r2, g_format_out, fpo2, zippit)
               ;  g_nt_right_out += r2->seq_n * r2->count
            ;  }
               continue
         ;  }

                                                /* below assignment needs to be made before move_sequences;
                                                 * bit questionable design (stick left_n in struct record
                                                 * and let move_sequences deal with it would be neater).
                                                */
            left_n = r1->seq_n - early * sinsert

         ;  g_nt_left_in  += r1->seq_n * r1->count
         ;  if (g_have_paired)
            g_nt_right_in += r2->seq_n * r2->count

                                                /* this updates r1->seq_n and r1->q_n
                                                 * caters for both paired end and single file.
                                                */
         ;  if (move_sequences(r1, r2, early * sinsert))
            goto DONE

         ;  n_receive = r1->seq_n

                                                /* compression of NANANA etc doubles size; factor 2
                                                 * quality scores; factor 1 (char per base) plus zero terminator.
                                                 * Not sure whether this factors in NANANAN.
                                                */
         ;  if (storage_used + 3*n_receive + 3 >=storage_size)
            {  if (!(storage = malloc(storage_size)))
               X_ERR_JUMP(DONE, "no memory available")
            ;  cp = storage
            ;  storage_used = 0

            ;  t1 = clockit(t1, &time_task, &time_iput)
            ;  argh
               (  "tally"
               ,  "new chunk (size %lu, time %.1f, hash fill %.3f, %lu sequences)"
               ,  storage_size
               ,  time_task
               ,  current_hash->n_used * 1.0 / current_hash->capacity
               ,  old_hash_num + current_hash->n_used
               )
#if TRACING_ON
; {   static unsigned long g_cfl_prev = 0
  ;   int m
  ;   argh("trace", "conflicts %lu", (g_cfl - g_cfl_prev))
  ;   g_cfl_prev = g_cfl
  ;   for (m=0; m<256; m++)
      fprintf(stderr, " %.1f", log((1+g_hist[m])/log(10.0)))
  ;   puts("")
; }
#endif
            ;  n_chunk++
            ;  g_memusage += storage_size
         ;  }

            {  void* uptr
            ;  unsigned seq_delta = 0

            ;  if (use_compression)
                  n_compress = compress_sequence(r1->seq, r1->seq_n, cp)
               ,  len_compressed += n_compress
               ,  hstatus = stack_insert(lsohash, cp, n_compress, r1, left_n, &uptr)
               ,  seq_delta = n_compress
            ;  else
                  memcpy(cp, r1->seq, r1->seq_n)
               ,  hstatus = stack_insert(lsohash, cp, r1->seq_n, r1, left_n, &uptr)
               ,  seq_delta = r1->seq_n
            ;  len_uncompressed += r1->seq_n
            ;  if (hstatus == 'i')
               {  storage_used += seq_delta
               ;  cp += seq_delta
            ;  }

               if (!hstatus)
               X_ERR_JUMP(DONE, "hash failure -- increase -hs argument")

            ;  if (with_quality)
               {  if (hstatus == 'i')
                  {  memcpy(cp, r1->q, r1->q_n)
                  ;  cp[r1->q_n] = '\0'
                  ;  ((struct unit3*) uptr)->qp = cp
                  ;  cp += r1->q_n + 1
                  ;  storage_used += r1->q_n + 1
               ;  }
                  else                                /* if (thecount < 256) */
                  {  int i
                  ;  struct unit3* u3 = uptr
                  ;  char* q = u3->qp
                  ;  for (i=0;i<r1->q_n;i++)
                     {  unsigned qi = (unsigned char) q[i]
                     ;  unsigned n  = (unsigned char) r1->q[i]
                     ;  if (n > qi)
                        q[i] = n
                  ;  }
                  }
               }
            }
            if (limit && current_hash->n_used + old_hash_num >= limit)
            break
      ;  }

         if (g_gauge_memory)
         {  FILE* fmem = fopen("/proc/self/smaps", "r")
         ;  int c
         ;  if (!fmem)
            argh("tally", "no luck in getting private/diry estimate")
         ;  else
            while(EOF != (c = fgetc(fmem)))
            fputc(c, stdout)
      ;  }

         t1 = clockit(t1, &time_task, &time_iput)
      ;  argh("input", "done (time %.1f)", time_task)

      ;  if (input)  myfzclose(input, 1)
      ;  if (input2) myfzclose(input2, 1)

      ;  {  long unsigned nunique = 0
         ;  int i

         ;  for (i=0;lsohash[i].n_used > 0 && i < HCL_NHASH;i++)
            nunique += lsohash[i].n_used

         ;  n_hash_used = i

         ;  if (g_verbose)
            {  if (r1->roffset)
               argh("observed", "average read length: %.0f", sumreadsize * 1.0 / r1->roffset)
            ;  argh("observed", "read count: %.1fM", r1->roffset * 1.0 / (1 << 20))
            ;  if (sumrecsize)
               argh("observed", "read/record ratio: %.2f", sumreadsize * 1.0 / sumrecsize)
            ;  if (!g_have_paired)
               argh("observed", "nonredundant read count: %.1fM", nunique * 1.0 / (1 << 20))
            ;  if (inputfilesize)
               argh("observed", "gzip compression: %.1f", sumrecsize * 1.0 / inputfilesize)
         ;  }
         }
      }


      if (tally_modes & TALLY_NO_TALLY)
      {  status = 0
      ;  if (g_have_paired && n_pair_error)
         argh("_____", "WARNING %u records had non-matching identifiers", (unsigned) n_pair_error)
      ;  if (g_have_paired)
         argh("tally", "stringID error/compared/paired %d/%d/%d", (int) n_pair_error, (int) n_pair_id, (int) r1->n_paired)
      ;  argh("tally", "paired %d reads (out of %d for -i and %d for -j argument)", (int) r1->n_paired, (int) r1->roffset, (int) r2->roffset)
      ;  goto DONE
   ;  }

      {  if (g_have_paired)
         argh("tally", "skipped %u/%u (left/right)", r1->discarded_unmatched, r2->discarded_unmatched)
      ;  if (g_have_paired && n_pair_error)
         argh("_____", "WARNING %u records had non-matching identifiers", (unsigned) n_pair_error)
      ;  if (g_have_paired)
         argh("tally", "stringID error/compared/paired %d/%d/%d", (int) n_pair_error, (int) n_pair_id, (int) r1->n_paired)
      ;  if (g_tally_fieldtosmall)
         argh("_____", "WARNING %lu cases where field was truncated", (long unsigned) g_tally_fieldtosmall)

      ;  if (use_compression)
         argh
         (  "memory"
         ,  "DNA compression factor %.1f"
         ,  (double) (len_compressed ? len_uncompressed * 1.0 / len_compressed : 0.0)
         )
   ;  }

                                                      /* destroy hash; create array at start position */
                                                      /* funcify */
      {  unsigned i = 0, j
      ;  t1 = clock()
      ;  argh("array from hash", "start")
      ;  for (i=0;i < n_hash_used; i++)
         {  unsigned dest = 0
         ;  unsigned totalcount = 0
         ;  for (j=0;j<lsohash[i].capacity;j++)
            {  struct unit1* n = getunit1(lsohash[i].ls, j, usize)
            ;  totalcount += n->count
            ;  if (n->nt)
               {  if (n->count > 1)
                  lsohash[i].n_gtone++
               ;  if (j > dest)
                  memcpy(getcharp(lsohash[i].ls, dest, usize), getcharp(lsohash[i].ls, j, usize), usize)
               ;  dest++
            ;  }
            }
         ;  if (dest != lsohash[i].n_used)
               argh("assertion failed", "hash corruption detected (%u vs %u) -- file bug please", ju dest, ju lsohash[i].n_used)
            ,  argh("assertion failed", "continuing, fingers crossed")
         ;  if (dest < lsohash[i].capacity)
            memset(getcharp(lsohash[i].ls, dest, usize), 0, usize)
         ;  if (lsohash[i].n_used)
            argh
            (  "hash-stats"
            ,  "H%02d lucky %.3f fill %.3f unique, total, repeated %.1fM, %.1f%%, %.1f%%"
            ,  (int) i
            ,  lsohash[i].n_noconflict * 1.0/ lsohash[i].n_used
            ,  lsohash[i].n_used * 1.0 / lsohash[i].capacity
            ,  lsohash[i].n_used * 1.0 / (1 << 20)
            ,  totalcount * 100.0 / lsohash[i].n_used
            ,  lsohash[i].n_gtone * 100.0 / lsohash[i].n_used
            )
      ;  }
         t1 = clockit(t1, &time_task, &time_process)
      ;  argh
         (  "data-stats"
         ,  "%d chunks, used %.2f%%"
         ,  (int) n_chunk
         ,  100.0 * (storage_used + (n_chunk - 1) * storage_size) * 1.0 / (n_chunk * 1.0 * storage_size)
         )
      ;  argh("array from hash", "done (time %.1f)", time_task)
   ;  }

                                                      /* funcify */
      if (do_sort == 'C')
      {  unsigned i = 0

      ;  argh("singletons", "start - moving unique reads to the end")
      ;  stack_separate_singletons(lsohash, n_hash_used)

      ;  t1 = clockit(t1, &time_task, &time_process)
      ;  argh("singletons", "done (time %.1f)", time_task)
      ;  argh("sorting", "start")

      ;  for (i=0; i < n_hash_used; i++)
         {  unsigned ngtone = lsohash[i].n_gtone
         ;  unsigned nused  = lsohash[i].n_used
         ;  qsort(lsohash[i].ls, lsohash[i].n_gtone, usize, unit_cmp_count)
         ;  if (ngtone > 0 && getunit1(lsohash[i].ls, ngtone-1, usize)->count < 2)
            argh("_____", "a blemish, I am sad to say")
         ;  if (nused > ngtone && getunit1(lsohash[i].ls, ngtone, usize)->count > 1)
            argh("_____", "a jarring imperfection, sadly")
      ;  }

         t1 = clockit(t1, &time_task, &time_process)
      ;  argh("sorting", "done (time %.1f)", time_task)

      ;  if (n_hash_used > 1)
         {  argh("merging", "start - merging hashes for correct sorting")
         ;  stack_sortall_bycount(lsohash, n_hash_used)

         ;  t1 = clockit(t1, &time_task, &time_process)
         ;  argh("merging", "done (time %.1f)", time_task)
      ;  }
      }


      if (do_output)
      {  unsigned i, j, hh = 0         /* hh: hairy hack */
      ;  char uncompress[MAXFIELDSIZE]
      ;  int late = early ? 0 : 1
      ;  int delta = late * sinsert
      ;  t1 = clock()
      ;  argh("output", "start")
      ;  unsigned output_id = 1

      ;  for (i=0; i < n_hash_used ;i++)
         {  for (hh=0; hh < 2; hh++)
            {  unsigned thestart =  0
            ;  unsigned theend   =  lsohash[i].n_gtone

            ;  if (hh)                                      /* in second batch do all unique reads */
                  thestart = lsohash[i].n_gtone
               ,  theend   = lsohash[i].n_used

            ;  for (j=thestart;j<theend;j++)
               {  void* v  =  getcharp(lsohash[i].ls, j, usize)
               ;  struct unit1* u1 = v
               ;  struct unit2* u2 = v
               ;  struct unit3* u3 = v
               ;  int n_uncompress = 0, trint1 = 0, trint2 = 0
               ;  unsigned len1    =   g_have_paired ? u2->left_n : u1->length

               ;  if (use_compression)
                  n_uncompress = uncompress_sequence(u1->nt, u1->length, uncompress)

               ;  if (n_uncompress < 0)
                  continue

               ;  trint1   =     use_compression
                           ?  trintscore(uncompress + delta, n_uncompress - delta)
                           :  trintscore(u1->nt + delta, len1 - delta)
               ;  trint2   =     0

               ;  if (g_have_paired)
                  trint2   =     use_compression
                           ?  trintscore(uncompress + len1 + delta + SEP_PAIR_WIDTH, n_uncompress - len1 - delta - SEP_PAIR_WIDTH)
                           :  trintscore(u1->nt + len1 + delta + SEP_PAIR_WIDTH, u1->length - len1 - delta - SEP_PAIR_WIDTH)

               ;  n_pass_unique++
               ;  n_pass += u1->count

               ;  {  const char* s  =  uncompress
                  ;  unsigned ns    =  n_uncompress
                  ;  const char* q  =  NULL
                  ;  if (with_quality)
                     q = u3->qp

                  ;  if (!use_compression)
                        s  = u1->nt
                     ,  ns = u1->length
                  ;  if (fpo)
                     writeline(s, ns, delta, u1->count, len1, trint1, output_id, g_format_out, fpo, zippit, fpo2 ? 1 : 0, q)
                  ;  if (fpo2)
                     writeline(s, ns, delta, u1->count, len1, trint2, output_id, g_format_out, fpo2, zippit, 2, q)
               ;  }
                  output_id++
               ;  if (limit2 && output_id >= limit2)        /* bit ugly this */
                  break
            ;  }
               if (limit2 && output_id >= limit2)           /* bit uglier even */
               break
         ;  }
            if (limit2 && output_id >= limit2)              /* I really rest my case */
            break
      ;  }

         t1 = clockit(t1, &time_task, &time_oput)
      ;  argh("output", "done (time %.1f)", time_task)

      ;  {  double time_total = time_process + time_iput + time_oput
         ;  if (time_total)
            argh
            (  "time"
            ,  "time spent in input: %2d%% compute: %2d%% output: %2d%%"
            ,  (int) (100.0 * time_iput / time_total)
            ,  (int) (100.0 * time_process / time_total)
            ,  (int) (100.0 * time_oput / time_total)
            )
      ;  }
      }

      status = 0
   ;  DONE
      :
                        /* this must follow DONE above, to accommodate NO_TALLY mode */
      if (do_output)
      {  if (fpo)  myfzclose(fpo, zippit)
      ;  if (fpo2) myfzclose(fpo2, zippit)
   ;  }

      if (status != STATUS_CLINE_PROCESSING && do_output)
      {  FILE* thef = fpoh ? fpoh : stderr
      ;  const char* sep = g_have_paired ? "_left" : ""
      ;  fprintf(thef, "discarded%s_unmatched=%u\n", sep, r1->discarded_unmatched)
      ;  fprintf(thef, "discarded%s_alien=%u\n", sep, r1->discarded_alien)
      ;  fprintf(thef, "discarded%s_length=%u\n", sep, r1->discarded_length)
      ;  fprintf(thef, "discarded%s_trint=%u\n", sep, r1->discarded_trint)

      ;  if (g_have_paired)
         fprintf(thef, "discarded_left_byright=%u\n", r1->discarded_byother)
      ;  if (with_quality)
         fprintf(thef, "error%s_quality=%u\n", sep, r1->error_quality)

      ;  if (!g_have_paired)
         {  const char* warn = ""
         ;  if (g_nt_left_in - n_pass * sinsert != g_nt_left_out)
            warn = "______"
         ;  fprintf(thef, "%snt_in=%lu\n", warn, jlu g_nt_left_in)
         ;  fprintf(thef, "%snt_out=%lu\n", warn, jlu g_nt_left_out)
      ;  }
         else
         {  const char* warn_left = "", *warn_right = ""

         ;  if (g_nt_left_in - n_pass * sinsert + g_nt_left_bogon_diff != g_nt_left_out)
            warn_left = "______"

         ;  if (g_nt_right_in - n_pass * sinsert + - g_nt_left_bogon_diff != g_nt_right_out)
            warn_right = "______"

         ;  fprintf(thef, "discarded_right_unmatched=%u\n", r2->discarded_unmatched)
         ;  fprintf(thef, "discarded_right_alien=%u\n", r2->discarded_alien)
         ;  fprintf(thef, "discarded_right_length=%u\n", r2->discarded_length)
         ;  fprintf(thef, "discarded_right_trint=%u\n", r2->discarded_trint)
         ;  fprintf(thef, "discarded_right_byleft=%u\n", r2->discarded_byother)

         ;  fprintf(thef, "%snt_left_in=%lu\n", warn_left, jlu g_nt_left_in)
         ;  fprintf(thef, "%snt_left_out=%lu\n", warn_left, jlu g_nt_left_out)
         ;  fprintf(thef, "%snt_right_in=%lu\n", warn_right, jlu g_nt_right_in)
         ;  fprintf(thef, "%snt_right_out=%lu\n", warn_right, jlu g_nt_right_out)

         ;  if (n_bogon)
            {  long diff = g_nt_left_bogon_diff
            ;  if (diff < 0)
               diff = -diff
            ;  fprintf(thef, "%s%snt_hash_ambiguity_difference=%ld\n", warn_left, warn_right, (long) diff)
            ,  fprintf(thef, "%s%shash_ambiguity_count=%lu\n", warn_left, warn_right, jlu n_bogon)
         ;  }
            if (with_quality)
            fprintf(thef, "error_right_quality=%u\n", r2->error_quality)
      ;  }
         fprintf(thef, "passed_unique=%u\n", n_pass_unique)
      ;  fprintf(thef, "passed_total=%u\n", n_pass)
      ;  if (g_have_paired)
            fprintf(thef, "num_records_left=%u\n", r1->roffset)
         ,  fprintf(thef, "num_records_right=%u\n", r2->roffset)
      ;  else
         fprintf(thef, "num_records=%u\n", r1->roffset)

      ;  {  unsigned long lefttally
         =     r1->discarded_unmatched
            +  r1->discarded_alien
            +  r1->discarded_length
            +  r1->discarded_trint
            +  r1->discarded_byother
            +  n_pass

         ;  unsigned long righttally
         =     r2->discarded_unmatched
            +  r2->discarded_alien
            +  r2->discarded_length
            +  r2->discarded_trint
            +  r2->discarded_byother
            +  n_pass

         ;  if (lefttally != r1->roffset_counted)
            argh("_____", "count discrepancy in %sinput file (tally %lu counted offset %u)", g_have_paired ? "left " : "", lefttally, r1->roffset_counted)

         ;  if (g_have_paired && righttally !=  r2->roffset_counted)
            argh("_____", "count discrepancy in right input file (tally %lu counted offset %u)", righttally, r2->roffset_counted)
      ;  }
         argh("memusage", "%lu bytes", g_memusage)
   ;  }

      if (changed_name)
      argh("tally", "\n\nBeware! output name changed to have .gz suffix, now %s", my_fnout)

   ;  {  unsigned i
      ;  for (i=0; i < n_hash_used ;i++)
         free(lsohash[i].ls)
   ;  }

      if (fpoh) myfzclose(fpoh, 0)
   ;  if (my_fnout) free(my_fnout)

   ;  return status
;  }


#ifndef BUILD_R_BINDINGS

int main
(  int argc
,  const char* argv[]
)
   {  return tally_main(argc, argv)
;  }

#endif



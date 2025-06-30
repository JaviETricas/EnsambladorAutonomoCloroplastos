
/*
 * (C) Copyright 2011, 2012, 2013, 2014, 2015 European Molecular Biology Laboratory.
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

/* TODO:
 * Sparse hash implementation. ls member of struct df
 * n_max_sources unlimit at some point.
 *
 *    trace_sequences -> unroll_sequence1 -> test_sequence
 * basic logic is sane, but embellished, frayed, and cluttered.
 */


#ifndef DEBUG_ON
#define DEBUG_ON 0
#endif


#if 0
#define KMER_DEBUG 43363288
#endif



#include "version.h"
#include "trint.h"
#include "slib.h"
#include "sw.h"

#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>


#if 0
#define BASEMAP(b)   (b == 'A' ? 0 : b == 'C' ? 1 : b == 'G' ? 2 : b == 'T' ? 3 : 4)
#else
#define BASEMAP(b)   themap[b]
#endif
                              /* A C G T are set in main() */
unsigned themap[256] =
   {0,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5
   ,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5
   ,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5
   ,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5
   ,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5
   ,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5
   ,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5
   ,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5
   };

char dna[4] = { 'A', 'C', 'G', 'T' };

static long unsigned g_n_reads = 0;
static int g_minion_id = 1;
static int g_test = 0;
ZFILE  g_pairfp = NULL;

typedef long unsigned kun;         /* k-mer unsigned */
typedef unsigned ctr;         /* counter */
#define lu (long unsigned)
#define ju (unsigned)         /* just unsigned */
#define ji (int)              /* just int */


enum
{  MODE_ADAPTER = 1
,  MODE_SEQUENCES
,  MODE_GAUGE_ADAPTER
}  ;


unsigned hash_sequence
(  char* s
,  unsigned length
)
   {  unsigned i, h=0
   ;  for (i=0;i<length;i++)
      h = 0 + ((h << 6) ^ (h << 21) ^ (h >> 13) ^ (s[i] * 71523))
   ;  return h >> 6
;  }


struct adapter_candidate
{  char* seq
;  unsigned seq_n

;  double score_fanout
;  double density_start

;  int score_fanout_rank
;  int density_start_rank

;  double density_prefix
;  double prefix_fanout
;  double id_len
;  int id
;  int isaprefix
;
}  ;



static struct adapter_candidate* ac_g = NULL;
static unsigned ac_N = 0;
static unsigned ac_i = 0;

static int mode_g = 0;

int sort_ac_byscore
(  const void* xv
,  const void* yv
)
   {  const struct adapter_candidate * x = xv
   ;  const struct adapter_candidate * y = yv
   ;  if (x->isaprefix + y->isaprefix == 1)
      return x->isaprefix ? 1 : -1
   ;  return x->score_fanout > y->score_fanout ? -1 : x->score_fanout < y->score_fanout ? 1 : 0
;  }


int sort_ac_bydensity
(  const void* xv
,  const void* yv
)
   {  const struct adapter_candidate * x = xv
   ;  const struct adapter_candidate * y = yv
   ;  if (x->isaprefix + y->isaprefix == 1)
      return x->isaprefix ? 1 : -1
   ;  return x->density_start > y->density_start ? -1 : x->density_start < y->density_start ? 1 : 0
;  }

int sort_ac_byid
(  const void* xv
,  const void* yv
)
   {  const struct adapter_candidate * x = xv
   ;  const struct adapter_candidate * y = yv
   ;  if (x->isaprefix + y->isaprefix == 1)
      return x->isaprefix ? 1 : -1
   ;  return x->id_len > y->id_len ? -1 : x->id_len < y->id_len ? 1 : 0
;  }


struct kmer
{  kun kmer
;  ctr ttr[4]                /* counts of the four possible kmers 'to the right' */
;  ctr nout                  /* sum of ttr. */
;  ctr nin                   /* sum of incoming edges */
;  int heap_offset
#if 0
;  float read_ofs            /* average starting position in read */
;  float read_var
#endif
;
}  ;



#define X_LEFT_A    ( 1 <<  0 )
#define X_LEFT_C    ( 1 <<  1 )
#define X_LEFT_G    ( 1 <<  2 )
#define X_LEFT_T    ( 1 <<  3 )

#define X_RIGHT_A   ( 1 <<  4 )
#define X_RIGHT_C   ( 1 <<  5 )
#define X_RIGHT_G   ( 1 <<  6 )
#define X_RIGHT_T   ( 1 <<  7 )

#define X_LEFTRIGHT ((1 <<  8 ) -1)
#define X_LEFT_STAR ((1 <<  4 ) -1)
#define X_RIGHT_STAR (X_LEFT_STAR << 4)

#define X_SEEN      ( 1 <<  8 )
#define X_SOURCE    ( 1 <<  9 )
#define X_BUSY      ( 1 << 10 )

#define B_LEFT_DROP             11
#define X_LEFT_DROP       (1 << B_LEFT_DROP)
#define B_LEFT_FREQ             12
#define X_LEFT_FREQ       (1 << B_LEFT_FREQ)
#define B_LEFT_DEPTH            13
#define X_LEFT_DEPTH      (1 << B_LEFT_DEPTH)

#define B_RIGHT_DROP             14
#define X_RIGHT_DROP       (1 << B_RIGHT_DROP)
#define B_RIGHT_FREQ             15
#define X_RIGHT_FREQ       (1 << B_RIGHT_FREQ)
#define B_RIGHT_DEPTH            16
#define X_RIGHT_DEPTH      (1 << B_RIGHT_DEPTH)

#define X_LEFT_REJECT     (X_LEFT_DROP | X_LEFT_FREQ | X_LEFT_DEPTH)
#define X_RIGHT_REJECT   (X_RIGHT_DROP | X_RIGHT_FREQ | X_RIGHT_DEPTH)

#define X_CYCLE     ( 1 << 17 )


void fillbits
(  char* buf
,  unsigned bits
)
   {  int i, delta = 0
   ;  strncpy(buf, "ACGT,ACGT,@sbt", 16)
   ;  for (i=0;buf[i];i++)
      {  if (i/4 > delta)
         delta++
      ;  if (!(bits & (1 << i)))
         buf[i+delta] = ' '
   ;  }
   }


struct kannot                    /* kmer annotation */
{  kun kmer
;  double in_out_ratio
;  unsigned ccid                 /* connected component ID */
;  unsigned flags
;
}  ;


struct cc                        /* connected component */
{  unsigned  id
;  ctr       N
;  ctr       N_sink
;  ctr       N_source
;  ctr       N_bpoint            /* branching point */
;  unsigned  skip
;  ctr       N_io_fmin
;  ctr       N_io_xmax
;  ctr       N_io_tmin
;
}  ;


#define OPT_ACCEPT_VERBOSE    1 << 0

struct df                        /* data frame */
{  struct kmer*   ls
;  unsigned       n_ls
;  struct kannot* ws             /* this represents the filtered set of kmers, the working set */
;  unsigned       n_ws           /* n working space */
;  kun*           ccstack        /* size n_ws, working space */
;  ctr            n_sources      /* sources are put at end of ccstack */
;  unsigned       k
;  kun            kmask
;  unsigned       options
;  FILE*          fpaccept
;  FILE*          fpextra        /* fasta output for adapter */
;  FILE*          fpdebug
;  ZFILE          fpdata
                              /* fixme: bundle magic parameters */
;  double         io_xmax        /* maximum drop-off */
;  double         io_fmin        /* frequency of base */
;  double         io_tmin        /* transition depth (absolute number of base) */
;  unsigned       io_imin        /* absolute number of total kmer output */
;  unsigned       io_omin        /* absolute number of total kmer input */
;  unsigned       n_variant_source_max
;  unsigned       n_variant_component_max
;  unsigned       n_max_sources
;  unsigned       tailtri
;  unsigned       kmertri
;  struct cc*     ccannot        /* size n_ws, working space */
;  unsigned       n_show
;
}  ;


struct kmer* ls_g = NULL;



int unrepeatiness
(  unsigned k
,  kun kmer
,  kun kmask
)
   {  int j
   ;  int n_diff_min = k
   ;  for (j=1;j<=3;j++)
      {  int i, n_diff = 0
      ;  kun shiftmer
         =  (  (  (kmer << (2 * j))
               |  (kmer >> (2 * (k-j)))
               )
            &  kmask
            )
      ;  kun xor = shiftmer ^ kmer
      ;  for (i=0;i<k;i++)
         if ((xor >> (2 * i)) & 3)
         n_diff++
      ;  if (n_diff_min > n_diff)
         n_diff_min = n_diff
   ;  }
      return n_diff_min
;  }



void fillbuf
(  char* buf
,  kun   kmer
,  unsigned k
)
   {  buf[k] = '\0'
   ;  while (k--)
      {  buf[k] = dna[kmer & 3]
      ;  kmer = kmer >> 2
   ;  }
   }


void get_input_transition_depth
(  struct df* df
,  unsigned k
,  kun kmer
,  unsigned* output
,  unsigned depth
,  unsigned d           /* should move up to depth */
)
   {  unsigned i
   ;  for (i=0; i< 4; i++)
      {  kun shiftmer = (i << (2 * k - 2)) | (kmer >> 2)
      ;  if (d >= depth)
         {  unsigned entry = (shiftmer >> (2 * (k - depth))) & ((1 << ( 2 * depth)) -1)
         ;  output[entry] = df->ls[shiftmer].nout
      ;  }
         else
         get_input_transition_depth(df, k, shiftmer, output, depth, d+1)
   ;  }
   }



void get_input_transition
(  struct df* df
,  unsigned k
,  kun kmer
,  unsigned output[4]
)
   {  unsigned i
   ;  for (i=0; i< 4; i++)
      {  kun shiftmer = (i << (2 * k - 2)) | (kmer >> 2)
      ;  output[i] = df->ls[shiftmer].ttr[kmer & 3]
   ;  }
;  }


unsigned sift
(  struct df* df
,  ctr* siftinfo
)
   {  int i
   ;  ctr npass = 0
   ;  struct kmer *ls = df->ls
   ;  unsigned dictsize = df->n_ls

   ;  for (i=0; i<dictsize; i++)
      {  int pass = 0
      ;  if (ls[i].nout > 0 || ls[i].nin > 0)
         {  if (ls[i].nout < df->io_omin && ls[i].nin < df->io_imin)
            {  siftinfo[3]++
            ;  siftinfo[2] += ls[i].nout
         ;  }
            else
            {  pass = 1
            ;  siftinfo[5]++
            ;  siftinfo[4] += ls[i].nout
         ;  }
         }

         if (pass)
         npass++
      ;  else
            ls[i].nout = 0
         ,  ls[i].nin  = 0
   ;  }
      return npass
;  }


                           /* noteme added extra code to catch kmers with no out */
int cmp_kmer_out
(  const void* a
,  const void* b
)
   {  kun u1 = ((const struct kannot*) a)[0].kmer
   ;  kun u2 = ((const struct kannot*) b)[0].kmer
   ;  return ls_g[u1].nout < ls_g[u2].nout ? 1 : ls_g[u1].nout > ls_g[u2].nout ? -1 : 0
;  }

               /* use heap insert ...
                * top of heap is the smallest of the current set of values
                * we start with all zeroes 
               */
void get_top_by_out
(  struct df* df
,  unsigned topsize
)
   {  unsigned i
   ;  struct kmer *ls = df->ls
   ;  unsigned dictsize = df->n_ls
   ;  unsigned mytopsize = topsize | 1       /* make sure it's odd */
   ;  struct kannot* heap = calloc(mytopsize, sizeof heap[0])
   ;  struct kmer ls0 = ls[0]

   ;  df->ws   = heap

   ;  ls[0].nout  = 0             /* heap initialised to all-zero; points to AAA..AAA, temporarily set to zero */
   ;  ls[0].nin   = 0             /* heap initialised to all-zero; points to AAA..AAA, temporarily set to zero */

   ;  for (i=0;i<dictsize;i++)
      {  if (!ls[i].nout && ls[i].nin)                    /* dangersign invariant broken */
         ls[i].nout = 1
      ;  if (ls[i].nout > ls[heap[0].kmer].nout)
         {  unsigned p = 0                                /* parent */
         ;  unsigned c = 1
         ;  heap[0].kmer = i
         ;  while (c+1 < mytopsize)
            {  if (ls[heap[c+1].kmer].nout < ls[heap[c].kmer].nout)   /* take smallest of children */
               c++
            ;  if (ls[heap[p].kmer].nout > ls[heap[c].kmer].nout)
               {  heap[p] = heap[c]                       /* move smallest child up */
               ;  heap[c].kmer = i
               ;  p = c                                   /* next parent will be place of smallest child */
               ;  c = 2*p + 1
            ;  }
               else
               break
         ;  }
         }
      }
      ls_g = df->ls
   ;  ls[0] = ls0
   ;  qsort(df->ws, mytopsize, sizeof df->ws[0], cmp_kmer_out)
   ;  df->n_ws = topsize
   ;  argh("minion", "set workspace to size %d", (int) topsize)
;  }


ctr g_fieldtosmall = 0;

#define MAXFIELDSIZE 511


struct record
{  char seq       [MAXFIELDSIZE+1]
;  char q         [MAXFIELDSIZE+1]
;  char discard   [MAXFIELDSIZE+1]
;  char id        [MAXFIELDSIZE+1]

;  unsigned seq_n
;  unsigned q_n
;  unsigned discard_n
;  unsigned id_n

;  unsigned count
;  unsigned ID
;
}  ;


#define MINION_ERROR 3
#define MINION_NOMEM 2
#define MINION_DONE  1
#define MINION_OK    0


#define izblank(c) ((unsigned char) (c) == ' ' || (unsigned char) (c) == '\t')


int cpytofield
(  char* dest
,  char* src
,  int   n
)
   {  if (n > MAXFIELDSIZE)
         n = MAXFIELDSIZE
      ,  g_fieldtosmall++
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
   {  const char* fmtp = format
   ;  unsigned rlstat = 0
#define LINE_LIMIT 8191
   ;  char buf[LINE_LIMIT+1]
   ;  unsigned n_received = 0
   ;  unsigned n_truncated = 0
   ;  char* bufp = NULL, *curp = NULL, *bufz = NULL
   ;  int esc = 0
   ;  unsigned n_lines = 0

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

   ;  while (fmtp[0])
      {  if (!bufp)
         {  rlstat |= kraken_hookline(ip, fb, buf, LINE_LIMIT, &n_received, &n_truncated)
         ;  if (n_truncated || rlstat & (RL_DONE | RL_ERROR))
            break
         ;  bufp = buf
         ;  bufz = buf + n_received
         ;  esc  = 0
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

     case 'H':
               while (bufp < bufz && (unsigned char) bufp[0] != (unsigned char) fmtp[1])
               bufp++
            ;  if (bufp == bufz)
               goto DONE
            ;  rec->discard_n = cpytofield(rec->discard, curp, bufp -curp)
            ;  fmtp++
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
            case 'Q':
                  while (bufp < bufz && !izblank((unsigned char) bufp[0]))
                  bufp++
               ;  rec->q_n = cpytofield(rec->q, curp, bufp - curp)
               ;  break
               ;
            case 'I':
                  while (bufp < bufz && !izblank((unsigned char) bufp[0]))
                  bufp++
               ;  rec->id_n = cpytofield(rec->id, curp, bufp - curp)
               ;  break
               ;
            case 'R':
                  while (bufp < bufz && isalpha((unsigned char) bufp[0]))
                     bufp[0] = toupper((unsigned char) bufp[0])
                  ,  bufp++
               ;  rec->seq_n = cpytofield(rec->seq, curp, bufp - curp)
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

      if (n_truncated)
         argh("minion", "line too long")
      ,  rlstat |= RL_ERROR

   ;  if (n_lines && (fmtp[0] || (bufp && bufp - buf != n_received)))
      {  argh("minion", "parse error at line %lu (remaining format string [%s], buffer [%s])", jlu n_lines, fmtp, bufp)
      ;  return MINION_ERROR
   ;  }

      if (rlstat & RL_DONE)
      return MINION_DONE

   ;  if (rlstat & RL_ERROR)
      return MINION_ERROR

   ;  if (rlstat & RL_NOMEM)           /* fixme not checked yet */
      return MINION_NOMEM

   ;  return MINION_OK
;  }


unsigned kmer_from_buf(const char* buf, unsigned K, unsigned* m)
   {  unsigned hid = 0
   ;  unsigned mask = 0
   ;  int i
   ;  for (i=0;i<K;i++)
      {  unsigned b = themap[(unsigned char) buf[i]]
      ;  mask = mask << 1
      ;  if (b >= 4)
         mask |= 1
      ;  else
         hid |= b << 2*(K-1-i)
   ;  }
      *m = mask
   ;  return hid
;  }


unsigned parse_sequence
(  struct df* df
,  const char* s
,  int n
,  unsigned duplicity
,  unsigned seq_id
,  ctr count[2]
)
   {  int i = 0
   ;  struct kmer* ls = df->ls
   ;  unsigned k      = df->k
   ;  kun kmask  = df->kmask
   ;  kun kmer  = 0
   ;  unsigned nbits = 0
   ;  unsigned nmask = (1 << (k+1)) - 1         /* note: mask for k+1-mer */
   ;  unsigned index[64] = { 0 }
   ;  unsigned n_tri = 0, n_skipped = 0

   ;  for (i=0;i<n;i++)
      {  unsigned b = BASEMAP((unsigned char) s[i])
      ;  unsigned tri_in  = ((kmer << 2) | b) & 63
      ;  unsigned tri_out = (kmer >> (2*k-6)) & 63

      ;  nbits = nbits << 1
      ;  if (b > 3)
         nbits |= 1

      ;  if (i >= 2 && !(nbits & 7))
         n_tri += (index[tri_in]++ == 0)

      ;  if (i >= k && !((nbits >> (k-2)) & 7))
         n_tri -= (index[tri_out]-- == 1)

;if(0)fprintf(stderr, "ntri %d (i/o %d %d N %d %d)\n", ji n_tri, ji tri_in, ji tri_out, ji (nbits & 7), ji ((nbits >> (k-2)) & 7))
      ;  if (i >= k)
         {  if (!(nmask & nbits) && n_tri >= df->kmertri)
            {  unsigned h = kmer & kmask
            ;  ls[h].ttr[b] += duplicity
            ;  count[0] += duplicity
            ;  count[1] += ls[h].nout == 0
            ;  ls[h].nout += duplicity
            ;  ls[((kmer << 2) | b) & kmask].nin += duplicity
#if 0
            ;  float this_ofs = (i-k+1) * 1.0 * duplicity
            ;  ls[h].read_ofs += this_ofs
            ;  ls[h].read_var += pow(this_ofs - ls[h].read_ofs * 1.0 / ls[h].nout, 2.0)
#endif
         ;  }
            else if (n_tri < df->kmertri)
            n_skipped++
      ;  }
         kmer = (kmer << 2) | b           /* screws kmer up if N etc, but will be shifted out */
   ;  }
      return n_skipped
;  }


unsigned nroffbits
(  unsigned v
)
   {  const unsigned char lut[16] = {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4}
   ;  return lut[v & 0x0f] + lut[v >> 4]
;  }


                                       /* fixme: defensive test, make sure
                                        * that "X_RIGHT_? is set" ==> heap_offset >= 0
                                       */
unsigned consider_shift
(  kun kmer
,  struct df* df
,  unsigned char nb[4]
)
   {  int ho            =  df->ls[kmer].heap_offset
   ;  ctr* ttr          =  df->ls[kmer].ttr
   ;  struct kannot* ka =  ho < 0 ? NULL : df->ws+ho
   ;  unsigned n = 0
   ;  unsigned char nb2[4]

   ;  if (!ka)
      return 0

   ;  if (ka->flags & X_RIGHT_A)
      nb[n++] = 0
   ;  if (ka->flags & X_RIGHT_C)
      nb[n++] = 1
   ;  if (ka->flags & X_RIGHT_G)
      nb[n++] = 2
   ;  if (ka->flags & X_RIGHT_T)
      nb[n++] = 3

   ;  if (n <= 1) return n
   ;  memcpy(nb2, nb, 4)

   ;  if (ttr[nb[1]] > ttr[nb[0]])
         nb[1] = nb[0]
      ,  nb[0] = nb2[1]

   ;  if (n == 2) return 2
                                    /* now first two elements sorted */
   ;  if (ttr[nb[2]] > ttr[nb[0]])
      {  memmove(nb+1, nb, 2)
      ;  nb[0] = nb2[2]
   ;  }
      else if (ttr[nb[2]] > ttr[nb[1]])
         nb[2] = nb[1]
      ,  nb[1] = nb2[2]

   ;  if (n == 3) return 3
                                    /* now first three elements sorted */
   ;  {  int i
      ;  for (i=3; i>0; i--)
         if (ttr[nb[i-1]] > ttr[nb[3]])
         break
      ;  if (i != 3)
            memmove(nb+i+1, nb+i, 3-i)
         ,  nb[i] = nb2[3]
   ;  }
      return 4
;  }




   /* return value: the number of matching nucleotides.
   */
int match_kmers
(  kun left
,  kun right
,  unsigned k
)
   {  int i = 0
   ;  kun kmask = (1 << (2*k)) - 1
   ;  for (i=0;i<k;i++)
      if (  (left & (kmask >> 2*i)) == ((right >> 2*i)  & kmask))
      break
   ;  return k-i
;  }


int output_assemble
(  const char* msg
,  char* buf
,  int buf_n
,  int i_variant
,  int n_cycle
,  ctr max_seen
,  kun kmer_start
,  kun kmer_end
,  struct df* df
,  struct kannot* annot
)
   {  struct cc* ccannot = df->ccannot + df->ws[df->ls[kmer_start].heap_offset].ccid
   ;  fprintf(df->fpaccept, ">minion%d depth=%lu group=%u", g_minion_id++, lu max_seen, ju df->ws[df->ls[kmer_start].heap_offset].ccid)
   ;  if (slib_verbose_level > 0)
      fprintf
      (  df->fpaccept
      ,  " [source,branch,sink %u %u %u] [rejected min-f min-t max-x %u %u %u] [maxcount=%lu variant=%d seqlen=%d] [start,max,end %lu/%lu %lu %lu/%lu] [type=%s] [cycle=%d]"

      ,  ju ccannot->N_source
      ,  ju ccannot->N_bpoint
      ,  ju ccannot->N_sink

      ,  ju ccannot->N_io_fmin
      ,  ju ccannot->N_io_tmin
      ,  ju ccannot->N_io_xmax

      ,  lu max_seen
      ,  ji i_variant
      ,  ji buf_n

      ,  lu df->ls[kmer_start].nin
      ,  lu df->ls[kmer_start].nout
      ,  lu max_seen
      ,  lu df->ls[kmer_end].nin
      ,  lu df->ls[kmer_end].nout

      ,  msg
      ,  ji n_cycle
      )
   ;  fprintf(df->fpaccept, "\n%s\n", buf)
   ;  return 1
;  }


int test_sequence
(  const char* msg
,  char* buf
,  int buf_n
,  int i_variant
,  int n_cycle
,  ctr max_seen
,  kun kmer_start
,  kun kmer_end
,  struct df* df
,  struct kannot* annot
)
   {  unsigned k = df->k
   ;  int vb = 0
   ;

      if (0 && vb)
      {  int p = k, i
      ;  kun mymer
      ;  char az[] = "0|1)2>3=4+5*6#7"

      ;  fprintf(df->fpaccept, "acc\n")
      ;  fprintf(df->fpaccept, "accepti ---------- ")    /* kmer input */
      ;  for (mymer = kmer_start, p=k;p<=buf_n;p++)
         {  kun mymer_next    =  ((mymer << 2) | BASEMAP((unsigned char) buf[p])) & df->kmask
         ;  unsigned logit    =  0.5 + log(df->ls[mymer].nin * 1.0) / log(sqrt(10.0))
         ;  if (logit > 14)
            logit = 14
         ;  fputc(az[logit], df->fpaccept)
         ;  mymer = mymer_next
      ;  }
         fputc('\n', df->fpaccept)

      ;  mymer = kmer_start
      ;  fprintf(df->fpaccept, "accepto ---------- ")    /* kmer output */
      ;  for (mymer = kmer_start, p=k;p<=buf_n;p++)
         {  kun mymer_next    =  ((mymer << 2) | BASEMAP((unsigned char) buf[p])) & df->kmask
         ;  unsigned logit    =  0.5 + log(df->ls[mymer].nout * 1.0) / log(sqrt(10.0))
         ;  if (logit > 14)
            logit = 14
         ;  fputc(az[logit], df->fpaccept)
         ;  mymer = mymer_next
      ;  }
         fputc('\n', df->fpaccept)

      ;  mymer = kmer_start
      ;  fprintf(df->fpaccept, "acceptx ---------- ")    /* drop-off and skew */
      ;  for (i=0;i<k-1;i++) fputc(' ', df->fpaccept)
      ;  for (mymer = kmer_start, p=k;p<buf_n;p++)
         {  kun mymer_next    =  ((mymer << 2) | BASEMAP((unsigned char) buf[p])) & df->kmask
         ;  unsigned char c   =  '0' + (int) (0.5 + (df->ls[mymer].ttr[mymer_next & 3] * 10.0 / df->ls[mymer].nout))
         ;  fputc(c, df->fpaccept)
#if 0
         ;  int ho            =  df->ls[mymer].heap_offset
         ;  unsigned fl       =  ho < 0 ? 1 : df->ws[ho].flags
         ;  unsigned char c
            =     fl == 1
               ? '!'
               :  fl & X_RIGHT_DEPTH
               ?  'd'
               :  fl == (X_RIGHT_DROP | X_RIGHT_FREQ)
               ?  'e'
               :  fl == X_RIGHT_FREQ
               ?  'f'
               :  fl == X_RIGHT_DROP
               ?  'x'
               :  '-'
         ;  fputc(c, df->fpaccept)
#endif
         ;  mymer = mymer_next
      ;  }
         fputs("!\n", df->fpaccept)

      ;  mymer = kmer_start
      ;  fprintf(df->fpaccept, "acceptr ---------- ")    /* right looking */
      ;  for (i=0;i<k-1;i++) fputc(' ', df->fpaccept)
      ;  for (mymer = kmer_start, p=k;p<=buf_n;p++)
         {  kun mymer_next    =  ((mymer << 2) | BASEMAP((unsigned char) buf[p])) & df->kmask
         ;  int ho            =  df->ls[mymer].heap_offset
         ;  unsigned n_consider = ho < 0 ? 0 : nroffbits(df->ws[ho].flags & X_RIGHT_STAR)
;if(0)fputc(n_consider == 0 ? '!' : n_consider == 1 ? '-' : '0' + n_consider, df->fpaccept)
         ;  mymer = mymer_next
      ;  }
         fputc('\n', df->fpaccept)

      ;  mymer = kmer_start
      ;  fprintf(df->fpaccept, "acceptx ---------- !")    /* drop-off and skew */
      ;  for (mymer = kmer_start, p=k;p<buf_n;p++)
         {  kun mymer_next    =  ((mymer << 2) | BASEMAP((unsigned char) buf[p])) & df->kmask
         ;  unsigned char c   =  '0' + (int) (0.5 + (df->ls[mymer].ttr[mymer_next & 3] * 10.0 / df->ls[mymer_next].nin))
         ;  fputc(c, df->fpaccept)
#if 0
         ;  int ho            =  df->ls[mymer].heap_offset
         ;  unsigned fl       =  ho < 0 ? 1 : df->ws[ho].flags
         ;  unsigned char c
            =     fl == 1
               ? '!'
               :  fl & X_LEFT_DEPTH
               ?  'd'
               :  fl == (X_LEFT_DROP | X_LEFT_FREQ)
               ?  'e'
               :  fl == X_LEFT_FREQ
               ?  'f'
               :  fl == X_LEFT_DROP
               ?  'x'
               :  '-'
         ;  fputc(c, df->fpaccept)
#endif
         ;  mymer = mymer_next
      ;  }
         fputc('\n', df->fpaccept)

      ;  mymer = kmer_start
      ;  fprintf(df->fpaccept, "acceptl ---------- ")    /* left looking */
      ;  for (mymer = kmer_start, p=k;p<=buf_n;p++)
         {  kun mymer_next    =  ((mymer << 2) | BASEMAP((unsigned char) buf[p])) & df->kmask
         ;  int ho            =  df->ls[mymer].heap_offset
         ;  unsigned n_consider = ho < 0 ? 0 : nroffbits(df->ws[ho].flags & X_LEFT_STAR)
         ;  fputc(n_consider == 0 ? '!' : n_consider == 1 ? '-' : '0' + n_consider, df->fpaccept)
         ;  mymer = mymer_next
      ;  }
         fputc('\n', df->fpaccept)

   ;  }

      if (mode_g == MODE_ADAPTER)
      {  int p = k, i
      ;  kun mymer
      ;  unsigned lpc[4] = { 0 }    /* left prefix counts */
      ;  unsigned transitions[64] = { 0 }
      ;  unsigned isa_prefix = 0
      ;  double prefix_sumsq = 0.0, prefix_count  = 0.0

      ;  struct adapter_candidate* ac = ac_g + ac_i
      ;  if (ac_i >= ac_N)
         return 1
      ;  ac_i++

      ;  get_input_transition(df, k, kmer_start, lpc)
      ;  get_input_transition_depth(df, k, kmer_start, transitions, 3, 1)

      ;  for (i=0; i< 64; i++)
         {  prefix_sumsq += pow(transitions[i], 2.0)
         ;  prefix_count += transitions[i]
      ;  }

         {  double thespread = prefix_count ? prefix_sumsq / pow(prefix_count, 2.0) : 1.0
         ;  double density_prefix = prefix_count * 100.0 / g_n_reads
         ;  double density_start = df->ls[kmer_start].nout * 100.0 / g_n_reads
         ;  double score_nominator = density_start * thespread

         ;  ac->score_fanout = score_nominator ? density_prefix / score_nominator : 0.0

               /* if density_prefix / density_start < 1, then correct */
         ;  if (density_prefix && density_start > density_prefix)
            ac->score_fanout *= density_start * 1.0 / density_prefix

         ;  ac->prefix_fanout = thespread
         ;  ac->density_prefix = density_prefix
         ;  ac->density_start = density_start
         ;  ac->density_start_rank = ac_i
         ;  ac->seq = malloc(buf_n+1)
         ;  if (ac->seq)
               memcpy(ac->seq, buf, buf_n)
            ,  ac->seq_n = buf_n
            ,  ac->seq[ac->seq_n] = '\0'
      ;  }

         for (mymer = kmer_start, p=k;p<=buf_n;p++)
         {  kun mymer_next    =  ((mymer << 2) | BASEMAP((unsigned char) buf[p])) & df->kmask
         ;  unsigned count_out =  df->ls[mymer].nout
         ;  double density     = count_out * 100.0 / g_n_reads
         ;  if (0)
            {  fprintf
               (  df->fpaccept
               ,  "%-10.3f  %.*s %12u"
               ,  density
               ,  (int) k, buf+p-k
               ,  (unsigned) mymer
               )
            ;  if (p < buf_n)
               fprintf(df->fpaccept, " (%c: %.3f)\n", (int) buf[p], df->ls[mymer].ttr[BASEMAP((unsigned char) buf[p])] * 100.0 / count_out)
         ;  }
            mymer = mymer_next               /* note: when p == buf_n, a 0 was shifted in */
      ;  }
         for (i=0;i<4;i++)
         {  kun mymer_test = (mymer | i) & df->kmask
         ;  unsigned count_in    =  df->ls[mymer_test].nin
         ;  if(0)fprintf(df->fpaccept, " (%c: %.3f)", dna[i], count_in * 100.0 / g_n_reads)
         ;  isa_prefix |=  (df->ws[df->ls[mymer_test].heap_offset].flags & X_SOURCE)
      ;  }
         ac->isaprefix = isa_prefix > 0
   ;  }

      if (mode_g == MODE_SEQUENCES)            /* sequence */
      {  struct cc* ccannot = df->ccannot + df->ws[df->ls[kmer_start].heap_offset].ccid
      ;  fprintf(df->fpaccept, ">minion%d depth=%lu group=%u", g_minion_id++, lu max_seen, ju df->ws[df->ls[kmer_start].heap_offset].ccid)
      ;  if (slib_verbose_level > 0)
         fprintf
         (  df->fpaccept
         ,  " [source,branch,sink %u %u %u] [rejected min-f min-t max-x %u %u %u] [maxcount=%lu variant=%d seqlen=%d] [start,max,end %lu/%lu %lu %lu/%lu] [type=%s] [cycle=%d]"

         ,  ju ccannot->N_source
         ,  ju ccannot->N_bpoint
         ,  ju ccannot->N_sink

         ,  ju ccannot->N_io_fmin
         ,  ju ccannot->N_io_tmin
         ,  ju ccannot->N_io_xmax

         ,  lu max_seen
         ,  ji i_variant
         ,  ji buf_n

         ,  lu df->ls[kmer_start].nin
         ,  lu df->ls[kmer_start].nout
         ,  lu max_seen
         ,  lu df->ls[kmer_end].nin
         ,  lu df->ls[kmer_end].nout

         ,  msg
         ,  ji n_cycle
         )
      ;  fprintf(df->fpaccept, "\n%s\n", buf)
   ;  }
      return 1
;  }

                                    /* fixme; bit expensive perhaps to recurse for traversal.
                                     * It's very convenient though.
                                    */
ctr unroll_sequence1
(  char* buf
,  int n_buf
,  int i_buf            /* write next base here */
,  int* n_variants
,  int* n_cycle
,  ctr* max_seen
,  kun root
,  kun kmer
,  struct df* df
)
   {  struct kmer* ls = df->ls
   ;  struct kannot* ws = df->ws
   ;  kun kmask = df->kmask
   ;  unsigned char nb[4] = { 0 }
   ;  int ret = 0
   ;  int i = 0

   ;  unsigned n_consider = consider_shift(kmer, df, nb)

#ifdef KMER_DEBUG
;if(kmer == KMER_DEBUG)fprintf(stderr, "%d %d %d %d %d %d flags %d\n", KMER_DEBUG, n_consider, ji nb[0], ji nb[1], ji nb[2], ji nb[3], ji df->ws[df->ls[kmer].heap_offset].flags)
#endif

   ;  if (n_variants[0] >= df->n_variant_source_max)                  /* fixme magic constant */
      return 0

   ;  buf[i_buf] = '\0'

                                                /* hierverder: contorted logic. */
   ;  if (i_buf == n_buf || !n_consider)
      {  const char* msg = !n_consider ? "curtains" : "cutoff"
      ;  if (test_sequence(msg, buf, i_buf, n_variants[0], n_cycle[0], max_seen[0], root, kmer, df, NULL))
         {  n_variants[0]++
         ;  ret++ 
      ;  }
         if (i_buf == n_buf)
         return 0
   ;  }

                              /* fixme: if accepted with prefix, the prefix can be something
                               * previously discarded -- ideally check on accepted only
                              */
      for (i=0; i< n_consider; i++)
      {  unsigned j = nb[i]
      ;  kun shiftmer = (kmer << 2 | j) & kmask 
      ;  ctr max_seen2 = max_seen[0]
      ;  struct kannot* shift_annot = ws+ls[shiftmer].heap_offset

      ;  if (shift_annot->flags & X_BUSY)           /* already seen in this consensus read */
         {  n_cycle[0]++
;if (df->fpdebug)fprintf(df->fpaccept, "cycle\n")
         ;  continue
      ;  }
         shift_annot->flags  |= X_BUSY

      ;  if (df->ls[shiftmer].nout > max_seen2)
         max_seen2 = df->ls[shiftmer].nout

      ;  buf[i_buf] = dna[j]
      ;  buf[i_buf+1] = '\0'
                                                   /* start of consensus sequence already accepted */
                                                   /* fixme: expand_component should prevent this?
                                                    * note that io tilt introduces *extra* breaks.
                                                   */
      ;  if (shift_annot->flags & X_SOURCE)
         {  if (test_sequence("prefix", buf, i_buf+1, n_variants[0], n_cycle[0], max_seen2, root, shiftmer, df, shift_annot))
            {  n_variants[0]++
            ;  ret++
         ;  }
         }
         else
         {  shift_annot->flags |= X_SEEN
         ;  ret += unroll_sequence1(buf, n_buf, i_buf+1, n_variants, n_cycle, &max_seen2, root, shiftmer, df)
      ;  }
         shift_annot->flags ^= X_BUSY
   ;  }
      return ret
;  }



struct branchpoint {
   kun   kmer;                      /* current kmer before shifting in next */
   int   n_branch;                  /* The number of options */
   int   i_branch;                  /* The current option */
   int   depth;                     /* Highest k-mer count in this sequence */
   unsigned char branches[4];       /* the options */
}  ;


void dump_bp
(  FILE* fp
,  struct branchpoint* bp
,  int k
,  int n
,  int printall
)
   {  int i, j
   ;  char buf[32]
   ;  for (i=0; i<n; i++)
      {  if (printall || bp[i].n_branch > 1)
         {  fillbuf(buf, bp[i].kmer, k)
         ;  fprintf(fp, "  %d", i)
         ;  for (j=bp[i].i_branch;j<bp[i].n_branch;j++)
            fprintf(fp, " %c", (unsigned int) dna[bp[i].branches[j]])
      ;  }
      }
      fputc('\n', fp)
;  }


/*
 *    nozashum: busy not set.
 *    zashum: busy set, corresponds with single k-mers.
*/
ctr expand_source
(  char* buf
,  int n_buf
,  struct branchpoint* branchpoints
,  int i_buf            /* write next base here */
,  kun start
,  struct df* df
)
   {  struct kmer* ls = df->ls
   ;  struct kannot* ws = df->ws, *current_annot
   ;  kun kmask = df->kmask
   ;  struct branchpoint *bp = branchpoints + 1             /* index 0 is sentinel value */
   ;  ctr n_cycle  = 0
   ;  ctr n_variants = 0
   ;  kun current = start
   ;  char printbuf[32]
   ;  struct branchpoint sentinel = { 0 }
   ;  int n_real_branches = 0
   ;  ctr n_ttt = 0

   ;  fillbuf(buf, start, df->k)
   ;  i_buf = df->k
   ;  buf[i_buf] = '\0'

   ;  bp->i_branch = 0
   ;  bp->n_branch = 1
   ;  bp->kmer = start

   ;  branchpoints[0] = sentinel

   ;  current_annot = ws+ls[start].heap_offset
   ;  current_annot->flags |= X_BUSY

#define TYPE_DEFAULT "leaf"
#define TYPE_CYCLE   "cycle"

   ;  while (1)
      {  unsigned char base = 0
      ;  unsigned n_consider = consider_shift(current, df, bp->branches)
      ;  const char* type = TYPE_DEFAULT

;fillbuf(printbuf, current, df->k)

;if(0)fprintf(stderr, "height %d real %d\n", (int) (bp - branchpoints), n_real_branches)
;if(0)fprintf(stderr, "%s consider %d\n", printbuf, (int) n_consider)
      ;  if (!n_consider)           /* go down the stack to find a fresh branch point. */
         {  A_STEP_BACK
         :  n_variants++

;  if (n_real_branches > 0)
   dump_bp(df->fpaccept, branchpoints, df->k, bp - branchpoints, 0)

         ;  if (strcmp(type, TYPE_DEFAULT))
            fprintf(df->fpaccept, "type %s at offset %d\n", type, ji i_buf)

         ;  if (mode_g != MODE_ADAPTER)
            output_assemble(type, buf, i_buf, n_variants, n_cycle, bp->depth, start, current, df, ws+ls[current].heap_offset)

                                    /* TODO: funcify below clause */
         ;  if (mode_g == MODE_ADAPTER)
            {  unsigned lpc[4] = { 0 }    /* left prefix counts */
            ;  unsigned transitions[64] = { 0 }
            ;  double prefix_sumsq = 0.0, prefix_count  = 0.0
            ;  unsigned i

            ;  struct adapter_candidate* ac = ac_g + ac_i
            ;  if (ac_i >= ac_N)
               return 1
            ;  ac_i++

            ;  get_input_transition(df, df->k, start, lpc)
            ;  get_input_transition_depth(df, df->k, start, transitions, 3, 1)

            ;  for (i=0; i< 64; i++)
               {  prefix_sumsq += pow(transitions[i], 2.0)
               ;  prefix_count += transitions[i]
            ;  }

               {  double thespread = prefix_count ? prefix_sumsq / pow(prefix_count, 2.0) : 1.0
               ;  double density_prefix = prefix_count * 100.0 / g_n_reads
               ;  double density_start = df->ls[start].nout * 100.0 / g_n_reads
               ;  double score_nominator = density_start * thespread

               ;  ac->score_fanout = score_nominator ? density_prefix / score_nominator : 0.0

                     /* if density_prefix / density_start < 1, then correct */
               ;  if (density_prefix && density_start > density_prefix)
                  ac->score_fanout *= density_start * 1.0 / density_prefix

               ;  ac->prefix_fanout = thespread
               ;  ac->density_prefix = density_prefix
               ;  ac->density_start = density_start
               ;  ac->density_start_rank = ac_i
               ;  ac->seq = malloc(i_buf+1)
               ;  if (ac->seq)
                     memcpy(ac->seq, buf, i_buf)
                  ,  ac->seq_n = i_buf
                  ,  ac->seq[ac->seq_n] = '\0'
            ;  }
            }

         ;  while (bp > branchpoints)
            {  n_ttt++
            ;  bp->i_branch++

            ;  if (bp->i_branch < bp->n_branch)       /* Note: bp still set to busy, as needed */
               break

            ;  if (ws[ls[bp->kmer].heap_offset].flags & X_CYCLE)     /* Do not flip X_BUSY, a lower branchpoint still points to this kmer */
               ws[ls[bp->kmer].heap_offset].flags ^= X_CYCLE
            ;  else if (!(ws[ls[bp->kmer].heap_offset].flags & X_BUSY))
               {  fillbuf(buf, bp->kmer, df->k)
               ;  dump_bp(df->fpaccept, branchpoints, df->k, bp-branchpoints, 1)
               ;  fprintf(stderr, "Error level %d %u %s %d type %s\n", (int) (bp-branchpoints), (unsigned) bp->kmer, buf, (int) bp->n_branch, type)
               ;  exit(1)
            ;  }
               else
               ws[ls[bp->kmer].heap_offset].flags ^= X_BUSY

            ;  if (bp->n_branch > 1)
               n_real_branches--

            ;  bp->i_branch = 0
            ;  bp--
            ;  i_buf--
         ;  }
            if (bp == branchpoints)
            break
      ;  }
         else
         {  bp->i_branch = 0
         ;  bp->n_branch = n_consider
         ;  if (n_consider > 1)
            n_real_branches++
      ;  }

         base     =  bp->branches[bp->i_branch]
      ;  current  = (bp->kmer << 2 | base) & kmask

      ;  bp++
      ;  bp[0]    =  sentinel       /* note: no branches yet */
      ;  bp->kmer =  current
      ;  bp->depth=  bp[-1].depth
      ;  current_annot = ws+ls[current].heap_offset
      ;  if (df->ls[current].nout > bp->depth)
         bp->depth = df->ls[current].nout

                                                           /* already seen in this consensus read
                                                            * Note that X_BUSY is attached to the kmer, not to
                                                            * the branchpoint.
                                                           */
      ;  if (current_annot->flags & X_BUSY)
         {  n_cycle++
         ;  current_annot->flags |= X_CYCLE
;fillbuf(printbuf, current, df->k)
;fprintf(df->fpaccept, "cycle detected at %s (level %d)\n", printbuf, (int) (bp - branchpoints))
         ;  type = TYPE_CYCLE
         ;  goto A_STEP_BACK
      ;  }

         current_annot->flags  |= X_BUSY     /* hierverder */

      ;  buf[i_buf++] = dna[base]
      ;  buf[i_buf] = '\0'

      ;  if (current_annot->flags & X_SOURCE)    /* expand_component should prevent this -- check */
         {  type = "fixme hit source"
         ;  goto A_STEP_BACK
      ;  }
         if (i_buf >= n_buf)
         {  type = "length limit"
         ;  goto A_STEP_BACK
      ;  }
      }

#if 0
      fillbuf(buf, start, df->k)
   ;  if (df->ws[df->ls[start].heap_offset].flags & X_BUSY)
      fprintf(df->fpaccept, "%d zashum %s (%d %d)\n", (int) g_minion_id, buf, (int) (bp - branchpoints), n_real_branches)
   ;  else
      fprintf(df->fpaccept, "\n%d nozashum %s (%d %d)\n\n", (int) g_minion_id, buf, (int) (bp - branchpoints), n_real_branches)
#endif
   ;  return n_variants
;  }



void annotate_top
(  struct df* df
)
   {  ZFILE fp = df->fpdata
   ;  struct kmer* ls = df->ls
   ;  unsigned dictsize = df->n_ls
   ;  struct kannot* ws = df->ws
   ;  unsigned topsize  = df->n_ws
   ;  unsigned k        = df->k

   ;  int i
   ;  char buf[128]

   ;  for (i=0; i< dictsize;i++)
      ls[i].heap_offset = -1

   ;  if (fp)
      gzprintf(fp, "kmer\tnmer\t#in\t#out\t#leftA\tA\t#leftC\tC\t#leftG\tG\t#leftT\tT\t#rightA\tA\t#rightC\tC\t#rightG\tG\t#rightT\tT\tior\tpskew\tsskew\n")

   ;  for (i=0;i<topsize;i++)
      {  kun kmer = ws[i].kmer
      ;  struct kannot* kan = ws+i
      ;  int j

      ;  ls[kmer].heap_offset = i

      ;  fillbuf(buf, kmer, k)
      ;  if (fp)
         gzprintf(fp, "%s\t%d\t%d\t%d", buf, ji kmer, ji ls[kmer].nin, ji ls[kmer].nout)
      ;  for (j=0;j<4;j++)
         {  kun prefix = (j << (2 *k -2)) | (kmer >> 2) 
         ;  if (fp)
            gzprintf(fp, "\t%d\t%c", ls[prefix].ttr[kmer & 3], dna[j])
      ;  }
         for (j=0;j<4;j++)
         {  if (fp)
            gzprintf(fp, "\t%d\t%c", ls[kmer].ttr[j], dna[j])
      ;  }
         kan->in_out_ratio = (1+ls[kmer].nin) * 1.0 / (1+ls[kmer].nout)
      ;  if (fp)
         gzprintf(fp, "\t%.1f\n", kan->in_out_ratio)
   ;  }
   }


                              /* hierverder.
                                 The way forward is to beat this into shape.
                                 Good logic for marking as 'done'
                                 while making sure edges are marked on both sides (left and right kmer)
                              */
void expand_component
(  struct df *df
,  kun kmer
,  unsigned id
,  struct cc* ccannot
)
   {  int i_stack = 1, n_seen = 1
   ;  kun kmask = df->kmask
   ;  unsigned k     = df->k
   ;  char buf[64]
   ;  char buf2[64]
   ;  kun* stack = df->ccstack    /* size df->n_ws */
   ;  int debug = 0

   ;  stack[0] = kmer

   ;  fillbuf(buf, kmer, k)

#if DEBUG_ON
   ;  if (df->fpdebug)
      fprintf(df->fpdebug, "\nstart %s in component %d\n", buf, ji id)
#endif

   ;  while (i_stack > 0)
      {  kun thismer  =  stack[i_stack-1]
      ;  unsigned n_left   =  0
      ;  unsigned n_right  =  0
      ;  int ho            =  df->ls[thismer].heap_offset
      ;  struct kannot* an =  ho < 0 ? NULL : df->ws+ho
      ;  struct kmer*   km =  df->ls+thismer        /* kmer middle */
      ;  int j

      ;  fillbuf(buf, thismer, k)
      ;  i_stack--

      ;  if (!an || an->ccid > 0)                   /* !an should not happen */
         continue

      ;  an->ccid  = id


;if(0)an->flags = 0

;if(0)fprintf(stderr, "stack %d top %s id %d\n", ji i_stack, buf, an->ccid)
;if (0 && debug++ > 40)
exit(1)
      ;  for (j=0;j<4;j++)
         {  kun ttr           =  ((thismer << 2) | j) & kmask                    /* to the right */
         ;  kun ttl           =  ((thismer >> 2) | (j << (2*k-2))) & kmask       /* to the left  */
         ;  struct kmer* kr   =  df->ls+ttr
         ;  struct kmer* kl   =  df->ls+ttl
         ;  int hr            =  kr->heap_offset
         ;  int hl            =  kl->heap_offset
         ;  struct kannot* ar =  hr >= 0 ? df->ws+hr : NULL
         ;  struct kannot* al =  hl >= 0 ? df->ws+hl : NULL

         ;  if (ar)
            {  const char* msg = "rskip"
            ;  unsigned fmin_ok
               =     km->ttr[j] >= km->nout * df->io_fmin
                  && km->ttr[j] >= kr->nin * df->io_fmin
            ;  unsigned xmax_ok
               =     kr->nin * df->io_xmax >= km->nout
                  && kr->nin <= km->nout * df->io_xmax
            ;  unsigned tmin_ok = km->ttr[j] >= df->io_tmin

            ;  an->flags |= ((1-fmin_ok) << B_RIGHT_FREQ) | ((1-xmax_ok) << B_RIGHT_DROP) | ((1-tmin_ok) << B_RIGHT_DEPTH)

            ;  ccannot->N_io_fmin += fmin_ok == 0
            ;  ccannot->N_io_xmax += xmax_ok == 0
            ;  ccannot->N_io_tmin += tmin_ok == 0

            ;  if (fmin_ok * xmax_ok * tmin_ok)
               {  n_seen++
               ;  n_right++
               ;  msg = "pickr"
               ;  an->flags |= 1 << (4+j)                /* danger sign; X_RIGHT_A ... X_RIGHT_T */
               ;  if (!ar->ccid)
                  stack[i_stack++] = ttr
               ;  else if (ar->ccid != id)
                  argh("_^^_", "no right symmetry :-(")
               ;  if (g_pairfp)
                     fillbuf(buf2, ttr, k)
                  ,  gzprintf(g_pairfp, "%s\t%s\n", buf, buf2)
            ;  }


#if DEBUG_ON
               fillbuf(buf2, ttr, k)
            ;  if (df->fpdebug)
               fprintf
               (  df->fpdebug
               ,  "%s from %s to %s test %d %d %d %d\n"
               ,  msg, buf, buf2
               ,  km->ttr[j] >= km->nout * df->io_fmin
               ,  km->ttr[j] >= kr->nin * df->io_fmin
               ,  kr->nin * df->io_xmax >= km->nout
               ,  kr->nin <= km->nout * df->io_xmax
               )
#endif
         ;  }

            if (al)
            {  const char* msg = "lskip"
            ;  unsigned fmin_ok
               =     kl->ttr[thismer & 3] >= kl->nout * df->io_fmin 
                  && kl->ttr[thismer & 3] >= km->nin  * df->io_fmin 
            ;  unsigned xmax_ok
               =     kl->nout * df->io_xmax >= km->nin
                  && kl->nout <= km->nin * df->io_xmax
            ;  unsigned tmin_ok = kl->ttr[thismer & 3] >= df->io_tmin

            ;  an->flags |= ((1-fmin_ok) << B_LEFT_FREQ) | ((1-xmax_ok) << B_LEFT_DROP) | ((1-tmin_ok) << B_LEFT_DEPTH)

            ;  ccannot->N_io_fmin += fmin_ok == 0
            ;  ccannot->N_io_xmax += xmax_ok == 0
            ;  ccannot->N_io_tmin += tmin_ok == 0

            ;  if (fmin_ok * xmax_ok * tmin_ok)
               {  n_seen++
               ;  n_left++
               ;  msg = "pickl"
               ;  an->flags |= 1 << j                    /* danger sign; X_LEFT_A ... X_LEFT_T */
               ;  if (!al->ccid)
                  stack[i_stack++] = ttl
               ;  else if (al->ccid != id)
                  argh("_^^_", "no left symmetry :-(")
               ;  if (g_pairfp)
                     fillbuf(buf2, ttl, k)
                  ,  gzprintf(g_pairfp, "%s\t%s\t%d\n", buf, buf2, (int) kl->ttr[thismer & 3])
            ;  }
#if DEBUG_ON
               fillbuf(buf2, ttl, k)
            ;  if (df->fpdebug)
               fprintf
               (  df->fpdebug
               ,  "%s from %s to %s test %d %d %d %d\n"
               ,  msg, buf, buf2
               ,  kl->ttr[thismer & 3] >= kl->nout * df->io_fmin 
               ,  kl->ttr[thismer & 3] >= km->nin  * df->io_fmin 
               ,  kl->nout * df->io_xmax >= km->nin
               ,  kl->nout <= km->nin * df->io_xmax
               )
#endif
         ;  }
         }

      ;  if (!n_right)
         {  ccannot->N_sink++
         ;  if (g_pairfp)
               gzprintf(g_pairfp, "//NODECLASS\t%s\tsink\tdebruijn\n", buf)
            ,  gzprintf(g_pairfp, "//NODESIZE\t%s\t16\n", buf)
      ;  }

                                    /* fixme: case no left due to cycle */
      ;  if (!n_left)
         {  stack[df->n_ws - ++df->n_sources] = thismer
         ;  ccannot->N_source++
         ;  if (g_pairfp)
               gzprintf(g_pairfp, "//NODECLASS\t%s\tsource\tdebruijn\n", buf)
            ,  gzprintf(g_pairfp, "//NODESIZE\t%s\t16\n", buf)
      ;  }

         ccannot->N_bpoint += (n_left > 1) || (n_right > 1)

#if DEBUG_ON
      ;  if (df->fpdebug)
         {  fillbuf(buf, thismer, k)
         ;  if (!n_right)
            fprintf(df->fpdebug, "  sink %s %d\n", buf, n_left)
         ;  if (!n_left)
            fprintf(df->fpdebug, "  source %s %d\n", buf, n_right)
         ;  if (n_left && n_right)
            fprintf(df->fpdebug, "  node %s %d %d\n", buf, n_left, n_right)
      ;  }
#endif

   ;  }
      if (!ccannot->N_source)
      {  fillbuf(buf, kmer, k)
      ;  argh("_()_", "lead %s yields no source", buf)
   ;  }
      ccannot->N = n_seen
;  }


#define N_TRACK_CCANNOT 20000


unsigned compute_cc
(  struct df* df
)
   {  unsigned i = 0, i_cc = 1               /* why at 1? just for reporting? */
   ;  unsigned topsize = df->n_ws
   ;  df->ccstack = calloc(topsize, sizeof df->ccstack[0])     /* this is a list of kmers */
   ;  char buf[64]

   ;  if (df->n_max_sources > N_TRACK_CCANNOT)
         argh("minion", "number of leads reduced to %d", (int) N_TRACK_CCANNOT)
      ,  df->n_max_sources = N_TRACK_CCANNOT

   ;  while (i_cc <= df->n_max_sources)
      {  struct cc *ccannot = df->ccannot+i_cc
      ;  while (i < topsize && df->ws[i].ccid)         /* seen as sequence or in cc */
         i++
      ;  if (i >= topsize)
         break

      ;  ccannot->id = i_cc                            /* vger-alloc, N_TRACK_CCANNOT */
      ;  expand_component(df, df->ws[i].kmer, i_cc, ccannot)

      ;  if (i_cc < 100)
            fillbuf(buf, df->ws[i].kmer, df->k)
         ,  vargh
            (  1
            ,  "lead"
            ,  "%s density %.2f%% component %d size %d source %d sink %d"
            ,  buf
            ,  ju df->ls[df->ws[i].kmer].nout * 100.0 / g_n_reads
            ,  ji i_cc
            ,  ji ccannot->N
            ,  ji ccannot->N_source
            ,  ji ccannot->N_sink
            )
      ;  else if (i_cc == 100)
         argh("minion", "further leads suppressed")

      ;  i_cc++
   ;  }
      return i_cc
;  }


                        /* hierverder; perhaps make sure that we are working in a single component?
                         * expand_component defines components, sources, sinks.
                         * missing sinks may have to do with io_ratio .. or similar criterion.
                         * trace_sequences + unroll_sequence1 work in parallel, ish.
                        */

#define MAX_SEQUENCE_LENGTH (1 << 20)
void trace_sequences
(  struct df* df
)
   {  unsigned i  = 0
   ;  ctr n_sequences = 0
   ;  ctr n_sources   = 0
   ;  ctr n_components = 0
   ;  ctr n_component_variants = 0
   ;  char* buf = myalloc(MAX_SEQUENCE_LENGTH)
   ;  struct branchpoint* bp = myalloc(MAX_SEQUENCE_LENGTH * sizeof bp[0])
   ;  kun ccid_prev = 0

   ;  vargh(1, "minion", "rolling in the deep (%u sources)", ju df->n_sources)

   ;  if (df->options & OPT_ACCEPT_VERBOSE)
      argh("minion", "log-scale: 0|1)2>3=4+5*6#7")

   ;  for (i=0; i<df->n_sources; i++)
      {  int ho = -1
      ;  int n_cycle = 0
      ;  kun start = df->ccstack[df->n_ws-i-1]
      ;  ctr max_seen = df->ls[start].nout
      ;  struct cc* can
      ;  int n_variants_deprecated = 0

      ;  ho = df->ls[start].heap_offset
      ;  fillbuf(buf, start, df->k)

      ;  vargh(2, "minion", "source %s", buf)

      ;  if (ho < 0)
         {  argh("minion", "skip %s", buf)
         ;  continue
      ;  }

         if (df->ws[ho].flags & X_BUSY)
         fprintf(stderr, "too busy: %s\n", buf)

      ;  df->ws[ho].flags |= X_SOURCE
      ;  if (g_test)
         df->ws[ho].flags |= X_BUSY
      ;  can = df->ccannot + df->ws[ho].ccid

      ;  if (can->id != ccid_prev)
            n_component_variants = 0
         ,  n_components++
         ,  ccid_prev = can->id

      ;  if (n_component_variants < df->n_variant_component_max)
         {  int n
         ;  if (df->options & OPT_ACCEPT_VERBOSE)
            fprintf
            (  df->fpaccept
            ,  "search %s (io %d/%d) in component %u\n"
            ,  buf
            ,  ji df->ls[start].nin
            ,  ji df->ls[start].nout
            ,  ju df->ws[ho].ccid
            )
         ;  if (g_test)
            n = unroll_sequence1(buf, MAX_SEQUENCE_LENGTH-1, df->k, &n_variants_deprecated, &n_cycle, &max_seen, start, start, df)
         ;  else
            n = expand_source(buf, MAX_SEQUENCE_LENGTH-1, bp, df->k, start, df)
         ;  n_component_variants += n
         ;  n_sequences += n
      ;  }

                                             /* fixme: n_variants, n_hit, n_variant_component_max, n_variant_source_max */
         if (df->ws[ho].flags & X_BUSY)
            fprintf(stderr, "busy: %s\n", buf)
         ,  df->ws[ho].flags ^= X_BUSY

      ;  if (df->n_max_sources && ++n_sources >= df->n_max_sources)
         break
   ;  }
      vargh
      (  1
      ,  "minion"
      ,  "luctor et emergo (%u, %u, %u sequences, sources, components)"
      ,  ju n_sequences
      ,  ju n_sources
      ,  ju n_components
      )
;  }


void output_settings
(  struct df* df
)
   {  if (!slib_verbose_level)
      return
   ;  argh("minion", "settings used:")
   ;  fprintf(stderr, " -max-x       %8g  maximum drop-off between k-mers\n", df->io_xmax)
   ;  fprintf(stderr, " -min-i       %8d  minimal required input of a node\n", df->io_imin)
   ;  fprintf(stderr, " -min-o       %8d  minimal required output of a node\n", df->io_omin)
   ;  fprintf(stderr, " -min-t       %8g  minimal required depth for any transition\n", df->io_tmin)
   ;  fprintf(stderr, " -min-f       %8g  minimal required base frequency for variant junctions\n", df->io_fmin)
   ;  fprintf(stderr, " -r           %8d  minimal required trimer count\n", df->kmertri)
   ;  fprintf(stderr, " -N           %8d  maximal number of components to analyse\n", df->n_max_sources)
   ;  fprintf(stderr, " -dust-suffix %8d  dust-type complexity criterion for stripping read suffixes\n", df->tailtri)
   ;  fprintf(stderr, " -k           %8d  word length\n", df->k)
;  }



void list_adapter_sequences
(  struct df* df
,  const char* theaptr
)
#define MATRIX_SIZE 8192
   {  SWNUM data[MATRIX_SIZE] = { 0 }
   ;  struct sw_alninfo ai = { 0 }
   ;  struct sw_param swp = { 0 }
   ;  int indel_allowed = 0
   ;  int i

   ;  swp.cost_gapleft  =  1
   ;  swp.cost_gapright =  1
   ;  swp.cost_subst    =  1
   ;  swp.gain_match    =  2
   ;  swp.left_skip     =  0
   ;  swp.right_limit   =  0
   ;  swp.flags         =  0

;fprintf(stderr, "have %d\n", (int) ac_i)
   ;  qsort(ac_g, ac_i, sizeof ac_g[0], sort_ac_bydensity)
   ;  for (i=0; i< ac_i;i++)
      ac_g[i].density_start_rank = i+1

   ;  qsort(ac_g, ac_i, sizeof ac_g[0], sort_ac_byscore)
   ;  for (i=0; i< ac_i;i++)
      ac_g[i].score_fanout_rank = i+1

   ;  if (theaptr)
      {  for (i=0; i<ac_i; i++)
         {  if (sw_fill(&ai, data, MATRIX_SIZE, ac_g[i].seq, theaptr, indel_allowed, &swp))
               ac_g[i].id_len = 0.0
            ,  ac_g[i].id = 0.0
         ;  sw_trace_trim(&ai, &swp, ai.max_ij, SW_TRIMLEFT |  SW_TRIMRIGHT)
         ;  {  double id_right = ai.n_match * 1.0 / (ai.nj-1+ai.n_insr)
            ;  double id_left  = ai.n_match * 1.0 / (ai.ni-1+ai.n_insl)
            ;  double id = id_right > id_left ? id_right : id_left
            ;  ac_g[i].id_len = id * ai.n_match
            ;  ac_g[i].id = 100.5 * id
         ;  }
         }
         qsort(ac_g, ac_i, sizeof ac_g[0], sort_ac_byid)
      ;  for (i=0; i<ac_i; i++)
         {  if (sw_fill(&ai, data, MATRIX_SIZE, ac_g[i].seq, theaptr, indel_allowed, &swp))
            {  fprintf(df->fpaccept, "Skipping sequence %s\n", ac_g[i].seq)
            ;  continue
         ;  }
            sw_trace_trim(&ai, &swp, ai.max_ij, SW_TRIMLEFT |  SW_TRIMRIGHT)

         ;  fputs("\n\n", df->fpaccept)
         ;  fputs("(predicted sequence)\n", df->fpaccept)
         ;  sw_pp_sparse(&ai, df->fpaccept, 0)
         ;  fputs("(query sequence)\n", df->fpaccept)
         ;  fprintf(df->fpaccept, "---\nmatch-score=%d\nmatch-count=%d\n", (int) ac_g[i].id, (int) ai.n_match)

         ;  fprintf
            (  df->fpaccept
            ,  "sequence-density=%.2f\n"
                  "sequence-density-rank=%d\n"
                     "fanout-score=%.2f\n"
                        "fanout-score-rank=%d\n"
                              "prefix-density=%.2f\n"
                                    "prefix-fanout=%.0f\n"
                                          "sequence=%s\n"
            ,  ac_g[i].density_start
            ,  (int) ac_g[i].density_start_rank
            ,  ac_g[i].score_fanout
            ,  (int) ac_g[i].score_fanout_rank
            ,  ac_g[i].density_prefix
            ,  ac_g[i].prefix_fanout ? 0.5 + 1.0 / ac_g[i].prefix_fanout : 1.0
            ,  ac_g[i].seq
            )
         ;  if (!df->n_show)
            fputs("Only the best match was shown (use -show <num> to see more)\n", df->fpaccept)
         ;  if (i+1 >= df->n_show)
            break
      ;  }
      }
      else
      {  int toggle = 0

      ;  for (toggle =0; toggle < 2; toggle++)
         {  qsort(ac_g, ac_i, sizeof ac_g[0], !toggle ? sort_ac_bydensity : sort_ac_byscore)
         ;  for (i=0; i<ac_i; i++)
            {  if (!toggle && ac_g[i].density_start_rank > df->n_show)
               break
            ;  if (toggle && ac_g[i].density_start_rank <= df->n_show)
               continue
            ;  if (toggle && ac_g[i].score_fanout_rank > df->n_show)
               break

            ;  fputs("\n\n", df->fpaccept)
            ;  if (df->fpextra)
               fprintf
               (  df->fpextra
               ,  ">minion-%d-%d density=%.0f dr=%d fr=%d\n%s\n"
               ,  ac_g[i].density_start_rank
               ,  ac_g[i].score_fanout_rank
               ,  ac_g[i].density_start
               ,  ac_g[i].density_start_rank
               ,  ac_g[i].score_fanout_rank
               ,  ac_g[i].seq
               )
            ;  fprintf
               (  df->fpaccept
               ,  "criterion=%s\n"
                     "sequence-density=%.2f\n"
                        "sequence-density-rank=%d\n"
                           "fanout-score=%.2f\n"
                              "fanout-score-rank=%d\n"
                                 "prefix-density=%.2f\n"
                                    "prefix-fanout=%.1f\n"
                                       "sequence=%s\n"
               ,  !toggle ? "sequence-density" : "fanout-score"
                  ,  ac_g[i].density_start
                     ,  (int) ac_g[i].density_start_rank
                        ,  ac_g[i].score_fanout
                           ,  (int) ac_g[i].score_fanout_rank
                              ,  ac_g[i].density_prefix
                                 ,  ac_g[i].prefix_fanout ? 1.0 / ac_g[i].prefix_fanout : 1.0
                                    ,  ac_g[i].seq
               )
            ;  if (i+1 >= df->n_show)
               break
         ;  }
         }
      }
   }


int minion_main
(  int argc
,  const char* argv[]
)
   {  ZFILE input = NULL
   ;  const char* g_fnin = "-"
   ;  const char* fndata = NULL, *fnout = "-", *fndebug = NULL, *fnextra = NULL
   ;  const char* theaptr = NULL
   ;  const char* g_format_in = "@%#%R%n%#%Q%n"
   ;  unsigned n_truncated = 0
   ;  unsigned nread = 0
   ;  struct df df = { 0 }
   ;  int status = 31
   ;  const char* gzopenmode = "w1"
   ;  unsigned id_prev = 0
   ;  int zippit = 1
   ;  int limit = 0

   ;  struct record rec    =  { { 0 }, { 0 }, { 0 }, { 0 }, 0 }

   ;  struct kmer* ls = NULL

   ;  struct file_buffer fb1 = { { 0 }, 0 }

   ;  kraken_readline(NULL, NULL, 0, NULL, &n_truncated)            /* resets static buffers */
   ;  kraken_hookline(NULL, &fb1, NULL, 0, NULL, &n_truncated)      /* resets buffers */

   ;  themap['A'] = 0
   ;  themap['C'] = 1
   ;  themap['G'] = 2
   ;  themap['T'] = 3
   ;  themap['N'] = 4
   ;  themap['X'] = 4

   ;  if
      (  argc < 2
      || (  strcmp(argv[1], "search-adapter")
         && strcmp(argv[1], "assemble")
         && strcmp(argv[1], "nofilter")
         && strcmp(argv[1], "lax-test")
         && strcmp(argv[1], "tiny-test")
         && strcmp(argv[1], "gauge-adapter")
         && strcmp(argv[1], "help")
         && !strstr(argv[1], "version")
         )
      )
      {  argh("minion", "Usage: minion <search-adapter|help|version>")
      ;  argh("example", "minion help | search-adapter | gauge-adapter | assemble | tiny-test | lax-test")
      ;  argh("example", "minion search-adapter -i FASTQFILE")
      ;  exit(0)
   ;  }

      df.tailtri  = 20
   ;  df.kmertri  = 4
   ;  df.io_xmax  = 5.0
   ;  df.n_variant_source_max = 10
   ;  df.n_variant_component_max = 100
   ;  df.n_max_sources = 0

   ;  if (strstr(argv[1], "version"))
      {  fprintf(stdout, "Minion version: %s\n", minion_tag)
      ;  exit(0)
   ;  }
      if (!strcmp(argv[1], "search-adapter"))
      {  df.k        =  12
      ;  df.io_fmin  =  0.6
      ;  df.io_xmax  =  2.0
      ;  df.io_tmin  =  100
      ;  df.io_imin  =  100
      ;  df.io_omin  =  100
      ;  df.n_max_sources =  50
      ;  df.n_show   =  1
      ;  limit       =  2000000
      ;  mode_g      =  MODE_ADAPTER
   ;  }
      else if (!strcmp(argv[1], "lax-test"))
      {  df.k       =   12
      ;  df.io_fmin =   0.05
      ;  df.io_xmax =   4.0
      ;  df.io_tmin =   5
      ;  df.io_imin =   5
      ;  df.io_omin =   5
      ;  df.n_max_sources =  1000
      ;  mode_g = MODE_SEQUENCES
   ;  }
      else if (!strcmp(argv[1], "nofilter"))
      {  df.k       =   12
      ;  df.io_fmin =   0.01
      ;  df.io_xmax =   100.0
      ;  df.io_tmin =   1
      ;  df.io_imin =   1
      ;  df.io_omin =   1
      ;  df.n_max_sources =  1000
      ;  mode_g = MODE_SEQUENCES
   ;  }
      else if (!strcmp(argv[1], "assemble"))
      {  df.k       =   14
      ;  df.io_fmin =   0.1
      ;  df.io_xmax =   2.0
      ;  df.io_tmin =   10
      ;  df.io_imin =   20
      ;  df.io_omin =   20
      ;  df.n_max_sources =  1000
      ;  mode_g = MODE_SEQUENCES
   ;  }
      else if (!strcmp(argv[1], "tiny-test"))
      {  df.k       =   7
      ;  df.io_fmin =   0.2
      ;  df.io_tmin  =  0
      ;  df.io_imin =   0
      ;  df.io_omin =   0
      ;  df.kmertri =   0
      ;  df.n_max_sources =  200
      ;  mode_g = MODE_SEQUENCES
   ;  }
      else if (!strcmp(argv[1], "gauge-adapter"))
      {  df.kmertri  =  0
      ;  df.io_tmin  =  3
      ;  df.io_imin  =  3
      ;  df.io_omin  =  3
      ;  mode_g = MODE_GAUGE_ADAPTER
      ;  df.k        =  7
   ;  }
      else if (!strcmp(argv[1], "help"))
      {
puts("-o <fname>        output stream (default STDOUT)");
puts("-i <fname>        input stream (gzipped file allowed) (default STDIN)");
puts("-d <fname>        (gzipped!) data stream (de Bruijn graph representation)");
puts("-e <fname>        debug output stream");
puts("-k <int>          kmer size (k)");
puts("-r <int>          minimal required k-mer tri-nucleotide content");
puts("-dust-suffix <int>  dust-type criterion for trimming read tails");
puts("-nvar <int>       maximum number of variants to track per graph component");
puts("-leads <int>      use no more than <int> leads (deBr graph components)");
puts("--fasta           expect fasta input");
puts("-record-format <fmt>  specify input format");
puts("-v                increase verbosity output levels");
puts("------ The options below have mode-dependent defaults ------");
puts("-nos <int>        maximum number of sources to consider");
puts("-min-i <int>      minimum required input of a kmer");
puts("-min-o <int>      minimum required output of a kmer");
puts("-max-x <num>      maximum allowed drop-off between nodes");
puts("-min-f <num>      minimum required base frequency between variant junctions");
puts("------ In search-adapter mode ------");
puts("-show <int>       show <int> top candidates for two criteria");
puts("-adapter <oligo>  cross-reference this adapter (consider -write-fasta and using swan)");
puts("-write-fasta <fname>    write output in FASTA format");
exit(0);
       }
       else
       { argh("minion", "no such mode: %s", argv[1])
       ; exit(1)
    ;  }

         arg_switch2()

      optarg("-record-format") g_format_in = thearg(); endarg()
      uniarg("--fasta") g_format_in = ">%#%R%n"; endarg()
      optarg("-o")   fnout   = thearg(); endarg()
      optarg("-do")  limit   = atoi(thearg());   endarg()
      optarg("-i")   g_fnin  = thearg(); endarg()
      optarg("-d")   fndata  = thearg(); endarg()
      optarg("-edges")   g_pairfp  = myfopen(thearg(), gzopenmode, 1); endarg()
      optarg("-e")   fndebug = thearg(); endarg()
      optarg("-adapter") theaptr = thearg(); endarg()
      optarg("-write-fasta") fnextra = thearg(); endarg()
      optarg("-k")   df.k    = atoi(thearg()); endarg()
      optarg("-r")   df.kmertri        = atoi(thearg()); endarg()
      optarg("-N")   df.n_max_sources       = atoi(thearg()); endarg()
      optarg("-min-i")  df.io_imin    = atoi(thearg()); endarg()
      optarg("-min-o")  df.io_omin    = atoi(thearg()); endarg()
      optarg("-min-t")  df.io_tmin    = atoi(thearg()); endarg()
      optarg("-min-io") df.io_omin = df.io_imin = atoi(thearg()); endarg()
      optarg("-max-x")  df.io_xmax = atof(thearg()); endarg()
      optarg("-min-f")  df.io_fmin    = atof(thearg()); endarg()
      optarg("-nvar")   df.n_variant_source_max = atoi(thearg()); endarg()
      optarg("-dust-suffix") df.tailtri  = atoi(thearg()); endarg()
      optarg("-show") df.n_show = atoi(thearg()); endarg()
      optarg("-leads") df.n_max_sources = atoi(thearg()); endarg()
      uniarg("-z") slib_verbose_level++; output_settings(&df); exit(0); endarg()
      uniarg("--test") g_test = 1; endarg()

      uniarg("-v")  slib_verbose_level++;  endarg()

      uniarg("--nozip") zippit = 0; endarg()
      uniarg("--version")
fprintf(stdout, "minion version: %s\n", minion_tag);
status = 0;
goto DONE;
      endarg()

      failarg()
      arg_done()
      ;

      if (mode_g == MODE_ADAPTER)
      {  ac_N =  df.n_max_sources
      ;  ac_g =  calloc(ac_N, sizeof ac_g[0])
   ;  }

      if (mode_g == MODE_GAUGE_ADAPTER && (!theaptr || strlen(theaptr) < df.k))
      X_ERR_JUMP(DONE, "mode gauge-adapter requires -adapter option")

   ;  if (!(input = myfopen(g_fnin, "r", 1)))
      X_ERR_JUMP(DONE, "bailing out (cannot read -i file)")

   ;  if (fnextra && !(df.fpextra = myfopen(fnextra, gzopenmode, 0)))
      X_ERR_JUMP(DONE, "cannot open FASTA output file")

   ;  if (fndata && !(df.fpdata = myfopen(fndata, gzopenmode, 1)))
      X_ERR_JUMP(DONE, "cannot open data output file")

   ;  if (!(df.fpaccept = myfopen(fnout, "w", 0)))
      X_ERR_JUMP(DONE, "cannot open main output file")

   ;  if (fndebug && !(df.fpdebug = myfopen(fndebug, "w", 0)))
      X_ERR_JUMP(DONE, "cannot open debug output file")

   ;  output_settings(&df)

#if 0
   ;  if (mode_g == MODE_SEQUENCES && df.n_show)
      df.n_max_sources = df.n_show
#endif

   ;  {  int i
      ;  df.n_ls = 1 << (2 * df.k)
      ;  if (!(ls = calloc(df.n_ls, sizeof ls[0])))
         X_ERR_JUMP(DONE, "cannot alloc kmer array")
      ;  vargh
         (  1
         ,  "minion"
         ,  "malloc'ed %.1fG, unit size %d"
         ,  (double) ((df.n_ls * 1.0 * sizeof ls[0]) / (1.0 * (1 << 30)))
         ,  ji sizeof ls[0]
         )
      ;  for (i=0;i<df.n_ls;i++)
         ls[i].kmer = i
      ;  df.kmask = (1 << 2*df.k) - 1
      ;  df.ls = ls
      ;  vargh(1, "minion", "now possess %u fine fields, %lu", ju df.n_ls, lu df.kmask)
   ;  }

      {  unsigned n_bases_read = 0, n_bases_removed = 0
      ;  ctr countinfo[2] = { 0 }
      ;  unsigned n_broken_reads = 0, n_kmers_skipped = 0
      ;  argh("minion", "reading reads")
      ;  while (1)
         {  unsigned n_skipped = 0

         ;  if (g_format_in)
            {  unsigned minionstat = read_record3(input, &fb1, g_format_in, &rec)

            ;  if (rec.ID && rec.ID <= id_prev)
               X_ERR_JUMP(DONE, "i-file %s record IDs not increasing (%u follows %u)\n", g_fnin, ju rec.ID, ju id_prev)
            ;  id_prev = rec.ID
            ;  if (minionstat == MINION_DONE)
               break
            ;  else if (minionstat == MINION_NOMEM || minionstat == MINION_ERROR)
               X_ERR_JUMP(DONE, "read error")
            ;  nread++
            ;  n_bases_read += rec.seq_n

            ;  if (++g_n_reads % 1000000 == 0)
               fprintf(stderr, "%3d\n", (int) g_n_reads/1000000)
            ;  else if (g_n_reads % 20000 == 0)
               fputc('.', stderr)

            ;  if (df.tailtri)
               {  int max = 0, maxid = 0
               ;  do
                  {  dustscore_tail(rec.seq, rec.seq_n, NULL, &max, &maxid)
                  ;  if (max >= df.tailtri)
                        rec.seq[maxid] = '\0'
                     ,  n_bases_removed += (rec.seq_n - maxid)
                     ,  rec.seq_n = maxid
               ;  }
                  while (max >= df.tailtri)
;if(0)fprintf(stdout, "read seq [%s] %d:%d\n", rec.seq, ji maxid, ji max)
            ;  }
               n_skipped = parse_sequence(&df, rec.seq, rec.seq_n, rec.count, nread, countinfo)
            ;  n_broken_reads  += (n_skipped > 0)
            ;  n_kmers_skipped += n_skipped
            ;  if (g_n_reads == limit)
               break
         ;  }
         }
         if (g_n_reads % 1000000) fputc('\n', stderr)
      ;  vargh
         (  1
         ,  "minion"
         ,  "have %u reads, removed %u bases by dusting (%.1f%% of total)"
         ,  ju g_n_reads
         ,  ju n_bases_removed
         ,  n_bases_removed * 100.0 / n_bases_read
         )
      ;  vargh
         (  1
         ,  "minion"
         ,  "skipped %u kmers in %u reads (%.1f average)"
         ,  ju n_kmers_skipped
         ,  ju n_broken_reads
         ,  n_broken_reads ? n_kmers_skipped * 1.0 / n_broken_reads : 0.0
         )
      ;  vargh
         (  1
         ,  "minion"
         ,  "have %.1fM total and %.1fM unique kmers (average %.1f)"
         ,  countinfo[0] / 1000000.0
         ,  countinfo[1] / 1000000.0
         ,  countinfo[0] * 1.0 / countinfo[1]
         )
   ;  }

      if (mode_g == MODE_GAUGE_ADAPTER)
      {  int i = 0
      ;  int l = strlen(theaptr)
      ;  unsigned mask = 0
      ;  unsigned kmer_prev = 0

;fprintf(stderr, "%i\n", l)
      ;  for (i=0; i<=l-df.k; i++)
         {  unsigned kmer = kmer_from_buf(theaptr+i, df.k, &mask)
         ;  unsigned nextbase = 0, lastbase = themap[(unsigned char) theaptr[i+df.k-1]]
         ;  unsigned a_kmer =   kmer >> 2
         ;  unsigned c_kmer =  (kmer >> 2 ) | (1 << (2 * df.k -2))
         ;  unsigned g_kmer =  (kmer >> 2 ) | (2 << (2 * df.k -2))
         ;  unsigned t_kmer =  (kmer >> 2 ) | (3 << (2 * df.k -2))

         ;  if (i < l-df.k)
            nextbase = themap[(unsigned char) theaptr[i+df.k]]

;if (nextbase > 3)
 die(1, "minion")
         ;  fprintf
            (  stderr
            ,  "%.*s to %.*u [ %7d %7d %7d %7d ]  [ %7d %7d %7d %7d ]  %u %u\n"
            ,  df.k
            ,  theaptr+i
            ,  df.k
            ,  kmer
            ,  (int) (0.5 + (1000000.0 * df.ls[a_kmer].ttr[lastbase]) / df.ls[kmer_prev].ttr[lastbase])
            ,  (int) (0.5 + (1000000.0 * df.ls[c_kmer].ttr[lastbase]) / df.ls[kmer_prev].ttr[lastbase])
            ,  (int) (0.5 + (1000000.0 * df.ls[g_kmer].ttr[lastbase]) / df.ls[kmer_prev].ttr[lastbase])
            ,  (int) (0.5 + (1000000.0 * df.ls[t_kmer].ttr[lastbase]) / df.ls[kmer_prev].ttr[lastbase])

            ,  (int) (0.5 + (1000000.0 * df.ls[kmer].ttr[0]) / df.ls[kmer].ttr[nextbase])
            ,  (int) (0.5 + (1000000.0 * df.ls[kmer].ttr[1]) / df.ls[kmer].ttr[nextbase])
            ,  (int) (0.5 + (1000000.0 * df.ls[kmer].ttr[2]) / df.ls[kmer].ttr[nextbase])
            ,  (int) (0.5 + (1000000.0 * df.ls[kmer].ttr[3]) / df.ls[kmer].ttr[nextbase])
            ,  df.ls[kmer].nin
            ,  df.ls[kmer].nout
            )
         ;  kmer_prev = kmer
      ;  }
         exit(0)
   ;  }

      {  ctr siftinfo[6] = { 0 }
      ;  if (df.kmertri || df.io_omin || df.io_imin)
         {  sift(&df, siftinfo)
         ;  vargh
            (  1
            ,  "minion"
            ,  "removed %.1fM total and %.1fM unique kmers below threshold (average %.1f)"
            ,  siftinfo[2] / 1000000.0
            ,  siftinfo[3] / 1000000.0
            ,  siftinfo[2] * 1.0 / siftinfo[3]
            )
         ;  vargh
            (  1
            ,  "minion"
            ,  "retained %.1fM total and %.1fM unique kmers (average %.1f)"
            ,  siftinfo[4] / 1000000.0
            ,  siftinfo[5] / 1000000.0
            ,  siftinfo[4] * 1.0 / siftinfo[5]
            )
      ;  }

         vargh(1, "minion", "sorting %lu kmers by frequency", (long unsigned) siftinfo[5])
      ;  get_top_by_out(&df, siftinfo[5])                /* fixme; heap no longer necessary. */
                                                         /* gaah: ^ what is this stupid comment */

      ;  annotate_top(&df)

                                                         /* vger-alloc */
      ;  df.ccannot = calloc(N_TRACK_CCANNOT, sizeof df.ccannot[0])

      ;  {        vargh(0, "minion", "connected component analysis")
         ;  compute_cc(&df)
         ;        vargh(0, "minion", "building consensus sequences")
         ;  trace_sequences(&df)

         ;  if (mode_g == MODE_ADAPTER)
            list_adapter_sequences(&df, theaptr)
      ;  }
      }

      status = 0
   ;  DONE
      :

      if (input)
      myfzclose(input, 1)

   ;  myfzclose(g_pairfp, 1)
   ;  myfzclose(df.fpdata, zippit)
   ;  myfzclose(df.fpaccept, 0)
   ;  myfzclose(df.fpdebug, 0)
   ;  myfzclose(df.fpextra, 0)

   ;  return status
;  }


int main
(  int argc
,  const char* argv[]
)
   {  return minion_main(argc, argv)
;  }






/*
 * (C) Copyright 2014, 2015 European Molecular Biology Laboratory.
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


/* TODO: memory read_fasta_file
 * kmer indexing: at the moment repeated kmers in query sequence are counted
 * as independent matches. Ideally we'd keep track of that.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include <errno.h>

#include "sw.h"
#include "slib.h"
#include "dna.h"
#include "version.h"



#define PRINT_KEY_VALUE 1
#define PRINT_EXCISE    2

#define GREP_MATCH      1
#define GREP_NOMATCH    2

enum {
   GREPV_FORMAT_IDLIST = 1,
   GREPV_FORMAT_FASTA
}  ;


static const char* g_format_grep = ">%I%n%Q%n";
static int g_grepv_format = GREPV_FORMAT_FASTA;
static int g_grep = 0;
static int g_cell = 0;
static int g_mode_print = 0;
static int g_qlen = 0;
static int g_rlen = 0;
static int g_indel_allowed = 1;
static int g_dump_matrix = 0;
static int g_compare_name = 0;
static unsigned g_n_seeds = 1;
static unsigned g_w_seeds = 0;

static unsigned g_n_alignments = 0;
static unsigned g_n_seed_skipped = 0;


struct req
{  char* s
;  int slen
;  char* annot
;  unsigned busy
;  unsigned last_kmer_independent_pos         /* for counting index hits */
;  unsigned last_kmer_pos
;  unsigned n_kmer_independent_hits    /* for counting index hits */
;  unsigned first_kmer_position
;
}  ;


#ifndef SWAN_K_TYPE
#     define SWAN_K_TYPE unsigned long
#endif

#define ktype SWAN_K_TYPE

struct kmer_x_req
{  ktype kmer
;  unsigned  req_index
;
}  ;


struct kmer_index
{  unsigned* kxr_offset
;  unsigned* kxr_length
;  struct kmer_x_req* kxr
;  unsigned k              /* defines size of the above */
;
}  ;


struct kmer_gen
{  unsigned cur_kmer
;  int last_N
;  int i
;  const char* s
;  unsigned slen
;  unsigned k
;
}  ;


void kmer_gen_init
(  struct kmer_gen* kg
,  const char* s
,  unsigned slen
,  unsigned k
)
   {  kg->cur_kmer = 0
   ;  kg->last_N   = -1
   ;  kg->i        = 0
   ;  kg->s        = s
   ;  kg->slen     = slen
   ;  kg->k        = k
;  }


unsigned* kmer_gen_next
(  struct kmer_gen* kg
)
   {  int i = kg->i
   ;  unsigned base = BASEMAP((unsigned char) kg->s[i])        /* will hit/need \0 */

   ;  if (i >= kg->slen)
      return NULL
   ;  kg->i++              /* note, i contains previous value */

   ;  kg->cur_kmer = (kg->cur_kmer << 2)
   ;  if (base == 4)
      kg->last_N = i
   ;  else if (base < 4)
      kg->cur_kmer |= base
   ;  kg->cur_kmer &= ((1 << (2*kg->k)) - 1)

   ;  if (kg->last_N + kg->k <= i)
      return &kg->cur_kmer
   ;  return kmer_gen_next(kg)
;  }


void sw_printaln5
(  struct sw_alninfo* ai
,  int accept
,  int score
,  void* fp
,  int zippit
,  int recno
)
   {  unsigned o = ai->aln_ofs
   ;  char space[512]
   ;  char buf[8192]
   ;  int n
   ;  int leftflush = ai->lft_start > ai->rgt_start
   ;  int shift = (int) ai->lft_start - (int) ai->rgt_start

   ;  memset(space, ' ', 512)
   ;  space[511] = '\0'

   ;  if (shift < 0)
      shift = -shift

   ;  n
   =  snprintf
      (  buf
      ,  8192
      ,  "%.*s%.*s%s%s : [%d,%d]\n%.*s%.*s%s : score %d\n%.*s%.*s%s%s : [%d,%d] %d/%d recno=%d\n.\n"
       /* 1 1 1 1 1 1     2  2    3 3 3 3 3          3   4 4 4 4 4 4    5  5   6  6      7 */
/* 1 */
      ,  leftflush ? 0 : shift         /* this much */
      ,  space
      ,  (int) (ai->lft_start-1)       /* this much */
      ,  ai->left
      ,  ai->aln_lft+o
      ,  ai->left + (ai->lft_end)
/* 2 */
      ,  (int) ai->lft_start
      ,  (int) ai->lft_end
/* 3 */
      ,  shift
      ,  space
      ,  (int) ((leftflush ? ai->rgt_start : ai->lft_start) - 1)
      ,  space
      ,  ai->aln_aln+o
      ,  (int) ai->data[ai->max_ij]
/* 4 */
      ,  leftflush ? shift : 0
      ,  space
      ,  (int) (ai->rgt_start -1)
      ,  ai->right
      ,  ai->aln_rgt+o
      ,  ai->right + (ai->rgt_end)
/* 5 */
      ,  (int) ai->rgt_start
      ,  (int) ai->rgt_end
/* 6 */
      ,  (int) accept
      ,  (int) score
/* 7 */
      ,  (int) recno
      )
#if WE_USE_ZLIB
   ;  if (n > 0)
      {
         if (zippit) gzwrite(fp, buf, n)
      ;  else        fputs(buf, fp)
#else
         fputs(buf, fp)
#endif
   ;  }
   }



   /* A stopgap solution.
   */
SWNUM* wrap_sw_fill
(  struct sw_alninfo* ai
,  const char *left           /* fixme, length-encode */
,  const char *right          /* fixme, length-encode */
,  int indel_allowed
,  struct sw_param* swp
)
   {  SWNUM* allocdata  = NULL
   ;  unsigned leftlen  = strlen(left)
   ;  unsigned rightlen = strlen(right)
   ;  unsigned long thesize  = (leftlen+1) * (rightlen+1)
   ;  int ret = 0
   
   ;  allocdata = myalloc(thesize * sizeof allocdata[0])

   ;  ret
      =  swp->flags & (SW_NW_FILL | SW_NW_CODE)
      ?  sw_fill_nw(ai, allocdata, thesize, left, right, indel_allowed, swp)
      :  sw_fill(ai, allocdata, thesize, left, right, indel_allowed, swp)

   ;  if (ret)
      {  free(allocdata)
      ;  die(1, "memory allocation failure, most likely")
   ;  }
      return allocdata
;  }


int pp_cb                       /* pretty print callback */
(  const struct sw_alninfo* ai
,  char* buf
,  unsigned bufsize
)
   {  return snprintf(buf, bufsize, " ref-offset=%d query-offset=%d aln={%s}", ai->lft_start, ai->rgt_start, ai->aln_aln+ai->aln_ofs)
;  }


int do_align
(  FILE* fpo
,  const char* seq_ref
,  const char* annot_ref
,  const char* seq_q
,  const char* annot_right
,  int recno_ref
,  int recno_q
,  struct sw_param* swp
,  const char* format
)
   {  SWNUM* thedata = NULL
   ;  struct sw_alninfo ai = { 0 }
   ;  int cell = 0

   ;  thedata = wrap_sw_fill(&ai, seq_ref, seq_q, g_indel_allowed, swp)

   ;  if (g_dump_matrix)
      sw_dump(&ai)

   ;  cell = ai.max_ij

   ;  if (g_cell && !(g_cell < ai.nj+1 || g_cell >= ai.ni * ai.nj))
      cell = g_cell

;if(1)fprintf(stderr, "DO ALIGN %d (%d) (%d)\n", (int) cell, (int) ai.nj+1, (int) (ai.ni * ai.nj -1))
   ;  if (swp->flags & (SW_NW_TRACE | SW_NW_CODE))
      sw_trace_nw(&ai, swp, cell)
   ;  else sw_trace(&ai, swp, cell)

   ;  if (format)
      sw_format(annot_ref, annot_right, cell, &ai, fpo, 0, recno_ref, NULL, format)
   ;  else if (g_mode_print & PRINT_KEY_VALUE)
      sw_lp(annot_ref, annot_right, cell, &ai, fpo, 0, recno_ref, recno_q,NULL)
   ;  else if (g_mode_print & PRINT_EXCISE)
      sw_pp_excise(annot_ref, annot_right, cell, &ai, fpo, 0, recno_ref, NULL)
   ;  else
      sw_pp2(annot_ref, annot_right, cell, &ai, fpo, 0, recno_ref, recno_q, NULL)

   ;  if (thedata)
      free(thedata)

   ;  return 0
;  }


   /* controlling k-mers.
    * We may find the same k-mer in the query sequence. Hmmmm.
    * Perhaps best have a repetitive sequence filter?
   */

int search_best_match
(  FILE* fpo
,  struct kmer_index* kidx
,  struct sw_param *swp
,  struct req* theref            /* many */
,  struct req* thequery          /* one  */
,  int* overhang
,  unsigned* todo
,  unsigned n_todo
,  unsigned id_threshold
,  int* n_foundp
,  int  recno_q
)
   {
#define MATRIX_SIZE 8192
   ;  struct sw_alninfo ai = { 0 }
   ;  struct kmer_gen kg
   ;  int best_i = -1, i, n_found = 0
   ;  int best_identity = 0
   ;  int best_overhang = 1 << 20
   ;  unsigned* kmerp
   ;  int with_index = !n_todo
   ;  int skip_by_grep = 0

   ;  if (with_index)
      {  kmer_gen_init(&kg, thequery->s, strlen(thequery->s), kidx->k)
      ;  unsigned offset = 0
      ;  while ((kmerp = kmer_gen_next(&kg)))
         {  int i
         ;  unsigned thekmer = kmerp[0]
;if (0 && kg.i == kidx->k)
fprintf(stderr, "(First k-mer %d %.*s seen %d times in index)\n", (int) thekmer, (int) kidx->k, thequery->s, (int) kidx->kxr_length[thekmer])
         ;  for (i=0; i<kidx->kxr_length[thekmer]; i++)
            {  int reqid = kidx->kxr[kidx->kxr_offset[thekmer]+i].req_index
            ;  if (!theref[reqid].busy)
               {  todo[n_todo++] = reqid
               ;  theref[reqid].busy = 1
               ;  theref[reqid].last_kmer_independent_pos = offset
               ;  theref[reqid].last_kmer_pos = offset
               ;  theref[reqid].first_kmer_position = offset
               ;  theref[reqid].n_kmer_independent_hits = 1
            ;  }
               else
               {  if (offset - theref[reqid].last_kmer_independent_pos >= kidx->k)
                  {  theref[reqid].n_kmer_independent_hits++
                  ;  theref[reqid].last_kmer_independent_pos = offset
;if(0)fprintf(stderr, "--- hit %d\n", (int) offset)
               ;  }
                  theref[reqid].last_kmer_pos = offset
            ;  }
            }
            offset++
      ;  }
;if(0)fprintf(stderr, "(Query hits %d index targets)\n", (int) n_todo);
      }

if (0 && with_index)
fprintf(stderr, "-- %s %d\n", thequery->annot, (int) n_todo);
      for (i=0; i<n_todo; i++)
      {  struct req* rs = theref+todo[i]
      ;  SWNUM* thedata = NULL
      ;  int skip_by_index
         =  (  with_index
            && rs->n_kmer_independent_hits < g_n_seeds
            && (!g_w_seeds || rs->last_kmer_pos - rs->first_kmer_position + kidx->k < g_w_seeds)
            )

;if(0)fprintf(stderr, "recs %d %d %d %d\n", rs->first_kmer_position, rs->last_kmer_independent_pos, rs->last_kmer_pos, rs->n_kmer_independent_hits)
                           /* this is a weak check ... need to worry about busy not being turned off */
      ;  if (with_index && !rs->busy)
         die(1, "impossibly not busy")


      ;  int skip =
         (  rs == thequery
         || (g_compare_name && !strcmp(thequery->annot, rs->annot))
         || skip_by_index
         || skip_by_grep
         )

      ;  g_n_seed_skipped += skip_by_index && rs->n_kmer_independent_hits > 0

      ;  rs->last_kmer_independent_pos = 0
      ;  rs->n_kmer_independent_hits = 0
      ;  rs->busy = 0

      ;  if (skip)
         continue

      ;  g_n_alignments++
      ;  thedata = wrap_sw_fill(&ai, rs->s, thequery->s, g_indel_allowed, swp)
;if (g_dump_matrix)
   sw_dump(&ai)

      ;  if (swp->flags & (SW_NW_TRACE | SW_NW_CODE))
         sw_trace_nw(&ai, swp, ai.ni * ai.nj -1)
      ;  else
         sw_trace(&ai, swp, ai.max_ij)

;if(0)fprintf(stderr, "found id %d\n", (int) ai.aln_identity)
      ;  if (id_threshold)
         {  if (id_threshold && ai.aln_identity >= id_threshold)
            {  n_found++
            ;  best_i = todo[i]
                                 /* we list all matches exceeding identity if not grepping */
            ;  if (!g_grep)
               do_align(fpo, rs->s, rs->annot, thequery->s, thequery->annot, todo[i]+1, recno_q, swp, NULL)
            ;  else if (g_grep & GREP_MATCH)
               skip_by_grep = 1
         ;  }
         }
         else if
         (  ai.aln_identity > best_identity
         || (ai.aln_identity == best_identity && ai.aln_overhang < best_overhang)
         )
         {  best_i = todo[i]
         ;  best_identity = ai.aln_identity
         ;  best_overhang = ai.aln_overhang
         ;  n_found = 1
      ;  }
         if (thedata)
         free(thedata)
   ;  }

      overhang[0] = best_overhang
   ;  n_foundp[0] = n_found

   ;  return best_i
;  }


const char* to_ACGT
(  const char* buf
,  int buf_n
,  int* n_nomatch
)
   {  char* dst = myalloc(buf_n+1)
   ;  int i, n_miss = 0

   ;  for (i=0; i<buf_n; i++)
      {  int c = ((unsigned char*) buf)[i]
      ;  if (isalpha(c))
         c = toupper(c)
      ;  if (c == 'U')
         c = 'T'
      ;  if (BASEMAP(c) > 4)
            n_miss++
         ,  c = 'N'
      ;  dst[i] = c
   ;  }
      dst[buf_n] = '\0'
   ;  n_nomatch[0] = n_miss
   ;  return dst
;  }


const char* warn_to_ACGT
(  const char* buf
,  int buf_n
)
   {  int n_nomatch = 0
   ;  const char* dst = to_ACGT(buf, buf_n, &n_nomatch)
   ;  if (n_nomatch)
      argh("swan", "strange characters in input oligo [%s], proceeding with [%s]", buf, dst)
   ;  return dst
;  }


                  /* requires 0 terminated list */
int find_dna
(  char* buf
,  int buf_n
,  int* start
)
   {  int best_start = 0, best_len = 0, s = 0, len = 0
   ;  char* a = buf, *z = buf + buf_n
   ;  while (a <= z)
      {  int c = ((unsigned char*) a)[0]
      ;  if (c == 'U')
         c = a[0] = 'T'
      ;  if (BASEMAP(c) < 4)
         {  if (!len)
            s = a - buf
         ;  len++
      ;  }
         else if (len > best_len)
         {  best_start = s
         ;  best_len = len
         ;  len = 0
      ;  }
         a++
   ;  }
      start[0] = best_start
   ;  return best_len
;  }


void printseq
(  struct req* rs
)
   {  fprintf(stderr, ">%s\n%s\n", rs->annot, rs->s)
;  }


unsigned req_total_kmer_length
(  struct req* rs
,  unsigned k
)
   {  unsigned tl = 0
   ;  while (rs->s)
      {  if (rs->slen >= k)
         tl += rs->slen + 1 -k
      ;  rs++
   ;  }
      return tl
;  }


int kxr_cmp_req
(  const void*  x
,  const void* y
)
   {  const struct kmer_x_req* xx = x
   ;  const struct kmer_x_req* yy = y
   ;  if (xx->req_index < yy->req_index)
      return -1
   ;  if (xx->req_index > yy->req_index)
      return 1
   ;  return 0
;  }


int kxr_cmp_kmr
(  const void*  x
,  const void* y
)
   {  const struct kmer_x_req* xx = x
   ;  const struct kmer_x_req* yy = y
   ;  if (xx->kmer < yy->kmer)
      return -1
   ;  if (xx->kmer > yy->kmer)
      return 1
   ;  return 0
;  }

#define KMASK(thek) ((((ktype) 1ULL) << 2*((ktype) thek)) -1)


unsigned destructive_quash       /* used with concordance; destruction is ok */
(  struct kmer_x_req* kidx
,  unsigned n
)
   {  unsigned i, target = 0
   ;  for (i=1; i<n; i++)
      {  if (kidx[i].req_index != kidx[target].req_index)
         {  target++ 
         ;  if (target != i)
            kidx[target] = kidx[i]
      ;  }
      }
      return target+1
;  }


unsigned long make_index
(  struct req* rs
,  struct kmer_index* kidx
,  unsigned k
,  int build_all
)
   {  unsigned total_length = req_total_kmer_length(rs, k)

   ;  struct req* r = rs
   ;  struct kmer_x_req* kxr
   ;  unsigned long xi = 0

   ;  kidx->k = k
   ;  if (!k)
      return 0

   ;  if (build_all)
      {  kidx->kxr_offset = myalloc((KMASK(k)+1) * sizeof kidx->kxr_offset[0])
      ;  kidx->kxr_length = myalloc((KMASK(k)+1) * sizeof kidx->kxr_length[0])
   ;  }
      kxr = kidx->kxr  = myalloc(total_length * sizeof kxr[0])

   ;  argh("swan", "building %d-mer index (use -index 0 to run without index) (kmer/llu/int18=%d/%d/%d)", (int) k, sizeof(ktype), sizeof(long long unsigned), sizeof( __int128))

   ;  while (r->s)
      {  unsigned i
      ;  ktype kmer = 0
      ;  int last_N = -1
      ;  for ( i=0; i < r->slen; i++ )
         {  unsigned base = BASEMAP((unsigned char) r->s[i])
         ;  kmer = (kmer << (ktype)2)
         ;  if (base == 4)
            last_N = i
         ;  else if (base < 4)
            kmer |= (ktype) base
         ;  if (last_N + k <= i)
            {  kxr[xi].kmer = kmer & KMASK(k)
            ;  kxr[xi].req_index = r - rs
;if(0)fprintf(stderr, "%llu %u\n", (long long unsigned) kxr[xi].kmer, kxr[xi].req_index)
            ;  xi++
         ;  }
         }
         r++
   ;  }

      if (build_all)
      {  int i
      ;  for (i=0; i<=KMASK(k); i++)         /* all counts to zero */
         kidx->kxr_length[i] = 0

      ;  for (i=0; i<xi; i++)                /* count what's there .. */
         kidx->kxr_length[kxr[i].kmer]++

      ;  kidx->kxr_offset[0] = 0             /* set offsets */
      ;  for (i=1; i<=KMASK(k); i++)
         kidx->kxr_offset[i] = kidx->kxr_offset[i-1] + kidx->kxr_length[i-1]

   ;  }
      qsort(kxr, xi, sizeof kxr[0], kxr_cmp_kmr)
   ;  if (build_all)
      {  int i
      ;  for (i=0; i<=KMASK(k); i++)         /* sort by sequence identifier/record-id */
         qsort(kxr+kidx->kxr_offset[i], kidx->kxr_length[i], sizeof kxr[0], kxr_cmp_req)
   ;  }

      argh("swan", "done building (%lu k-mers total)", xi)
   ;  return xi
;  }



   /* sequences can be arbitrarily long
    * Entire file is read into one chunk of memory,
    * then parsed into its constituent sequences.
   */ 

struct req* read_fasta_file
(  ZFILE fp
,  unsigned*  n_ref
,  const char* substr
,  int minlen
)
   {  unsigned long nbytes
   ;  char* fdata = (char*) read_a_file(fp, &nbytes)
   ;  char* a = fdata
   ;  char* start = NULL
   ;  unsigned nseq, i_rec = 0
   ;  struct req* rs = NULL
   ;  struct req rs0 = { 0 }
   ;  int n_skipped = 0

   ;  if (!fdata)
      exit(1)

   ;  nseq = a[0] == '>'

   ;  if (a[0] != '>')
      {  argh("fasta read", "file does not start with >, worryingly")
      ;  a = strstr(a, "\n>")
      ;  if (!a)
         a = fdata+nbytes
      ;  else
            nseq = 1
         ,  a++
   ;  }
      start = a

   ;  while ((a = strstr(a, "\n>")))
         nseq++
      ,  a++

   ;  fprintf(stderr, "read %d sequences\n", (int) nseq)
   ;  if (!(rs = malloc((nseq+1) * sizeof rs[0])))       /* add sentinel value */
      exit(1)

   ;  a = start
   ;  while (a && a[0] == '>')
      {  char* b = strstr(a, "\n"), *c = b, *x, *dest

      ;  if (!b)                       /* should not happen, really */
            argh("fasta read", "error")
         ,  exit(1)

      ;  rs[i_rec] = rs0
      ;  rs[i_rec].annot = a+1
      ;  c = strstr(b, "\n>")          /* these two are very much order-dependent */
      ;  rs[i_rec].annot[b-a-1] = '\0'

      ;  rs[i_rec].s = b+1
      ;  rs[i_rec].slen = c ? c - b - 1 : (nbytes - 1) - ((b+1) - fdata)
      ;  rs[i_rec].s[rs[i_rec].slen] = '\0'
      ;  dest = rs[i_rec].s

      ;  for (x = rs[i_rec].s; x < rs[i_rec].s + rs[i_rec].slen; x++)
         {  if (isalpha(x[0]))
            x[0] = toupper(x[0])
         ;  if (x[0] == 'U')
            x[0] = 'T'

         ;  if (BASEMAP((unsigned char) x[0]) < 4)
            /* NOTHING */
         ;  else if (isspace(x[0]))
            x[0] = '\0'
         ;  else
            x[0] = 'N'

         ;  if (x[0])
            dest++[0] = x[0]
      ;  }

         dest[0] = '\0'
      ;  rs[i_rec].slen = dest - rs[i_rec].s

      ;  a = c ? c+1 : NULL
      ;  if (substr && !strstr(rs[i_rec].annot, substr))
         n_skipped++
      ;  else if (minlen && rs[i_rec].slen < minlen)
         n_skipped++
      ;  else
         i_rec++
   ;  }
      if (i_rec + n_skipped == nseq)
      rs[i_rec] = rs0
   ;  else
      argh("fasta read", "record count discrepancy (%d/%d) - demons", (int) i_rec, (int) nseq)

   ;  n_ref[0] = i_rec
   ;  return rs
;  }


static char getbase[4] = { 'A', 'C', 'G', 'T' };
#define LSIZE(k)     (((ktype) 1) << (2*((ktype) k)))    /* language size, 4096, 16384          */
#define LOMEGA(k)    (LSIZE((ktype) k)-1)      /* last (all-one) word, 4095, 16383,   */

         /* Get the word that corresponds to a hash
         */
char* get_sylmer
(  int k
,  ktype sylid
,  char buf[101]
)
   {  int i

   ;  if (sylid > LOMEGA(k))
      {  for (i=0;i<k;i++)
         buf[i] = 'X'
   ;  }
      else
      for (i=0;i<k;i++)
      buf[i] = getbase[(sylid >> ((ktype)2*(k-i-1))) & (ktype)3]

   ;  buf[k] = '\0'
   ;  return buf
;  }



int main
(  int argc
,  char* argv[]
)
   {  FILE* fpo = stdout
   ;  FILE* fpo_grepv = stdout
   ;  ZFILE fpref = NULL, fpquery = NULL
   ;  const char* g_fnout = "-", *g_fnout_grepv = NULL
   ;  const char* g_string_query = NULL
   ;  const char* g_string_ref = NULL

   ;  struct sw_param swp = { 0 }
   ;  unsigned n_truncated = 0
   ;  unsigned index_kmer = 0
   ;  unsigned theidentity = 0
   ;  unsigned limit = 0
   ;  unsigned B_concordance = 0

   ;  const char* singleref = NULL, *singlequery = NULL
   ;  const char* fnref = NULL, *fnquery = NULL

   ;  swp.cost_gapleft      =  3
   ;  swp.cost_gapright     =  3
   ;  swp.cost_subst        =  1
   ;  swp.gain_match        =  4
   ;  swp.left_skip         =  0
   ;  swp.right_limit       =  0
   ;  swp.flags             =  SW_TRIMLEFT | SW_TRIMRIGHT

   ;  kraken_readline(NULL, NULL, 0, NULL, &n_truncated)            /* resets static buffers */
   ;  themap_init()

/* enter macromagical option world */
      
   ;  arg_switch()

      uniarg("--version")
fprintf(stdout, "Swan version: %s\n", swan_tag);
exit(0);
      endarg()

      uniarg("-h")
puts("-o                output file name (STDOUT)");
puts("-r FASTA-file     fasta file for reference");
puts("-q FASTA-file     fasta file for query");
puts("-rs DNA-string    reference string to align (displayed on top)");
puts("-qs DNA-string    query string to align (displayed below)");
puts("-q-len <int>      only consider sequences at least this long");
puts("-r-len <int>      only consider sequences at least this long");
puts("-q-string <string> e.g. hsa, mmu; only matching identifiers are considered");
puts("-r-string <string> e.g. hsa, mmu; only matching identifiers are considered");
puts("-identity <int>   require matches with at least <int> identity (0-100)");
puts("-index <int>      k-mer size to build index on (suggest 8 to 12; filters on k-mer match!)");
puts("-concordance <int> as -index, build index and output as concordance");
puts("-n-seeds <int>    require <int> independent k-mer hits for a match to be considered (overlap not allowed)");
puts("-w-seeds <int>    require two seeds spanning at least <int> bases (overlap allowed)");
puts("--grep            output sequences that match the reference (requires -identity)");
puts("--grepv           output sequences (see -grepv-o) that do not match the reference (requires -identity)");
puts("-grepv-o <fname>  output file for non-matching sequences");
puts("-swp M/S/G        match/substitution/gap : gain/cost/cost");
puts("-lsrl L/R         reference/left-skip / query/right-limit (adapter specific)");
puts("--noindel         do not consider indels while aligning");
puts("--skip-same-name  do not compare sequences with identical names");
puts("--matrix          dump alignment matrix");
puts("--key-value       output easily parseable line-based key-value output");
puts("--excise          excise the aligned part when printing");
puts("--nw              Needleman-Wunsch fill and trace (EXPERIMENTAL)");
puts("--nw-fill         Needleman-Wunsch fill (EXPERIMENTAL)");
puts("--nw-trace        Needleman-Wunsch trace (EXPERIMENTAL)");
puts("--no-trimleft     Do not trim alignment");
puts("--no-trimright    Do not trim alignment");
/* puts("--nw-code         Use NW code without triggering NW paths (PRIOR TO MERGING NW CODE)"); */
puts("--debug           Output diagnostic information");
puts("-do <int>         process the top <int> entries from the reference file");
puts("-cell <int>       align from cell <int>");
printf("(alignment unit: %s)\n", STRINGIFY(SWNUM));
exit(0);
      endarg()

      uniarg("--verbose") swp.flags |= SW_VERBOSE; endarg()
      uniarg("--no-trimleft") swp.flags |= SW_TRIMLEFT; swp.flags ^= SW_TRIMLEFT; endarg()
      uniarg("--notrim") swp.flags |= SW_TRIM; swp.flags ^= SW_TRIM; endarg()
      uniarg("--no-trimright") swp.flags |= SW_TRIMRIGHT; swp.flags ^= SW_TRIMRIGHT; endarg()
      uniarg("--skip-same-name")  g_compare_name = 1; endarg()
      uniarg("--nw")  swp.flags |= ( SW_NW_FILL | SW_NW_TRACE ); endarg()
      uniarg("--nw-fill") swp.flags |= SW_NW_FILL; endarg()
      uniarg("--nw-trace") swp.flags |= SW_NW_TRACE; endarg()
      uniarg("--nw-code") swp.flags |= SW_NW_CODE; endarg()
      optarg("-o")  g_fnout = thearg(); endarg()
      optarg("-r-string") g_string_ref = thearg(); endarg()
      optarg("-q-string") g_string_query = thearg(); endarg()
      optarg("-q-len") g_qlen = atoi(thearg()); endarg()
      optarg("-r-len") g_rlen = atoi(thearg()); endarg()
      optarg("-cell") g_cell = atoi(thearg()); endarg()
      optarg("-swp")
         if
         (  3 != sscanf(thearg(), "%u/%u/%u", &swp.gain_match, &swp.cost_subst, &swp.cost_gapleft)
         )
         die(1, "-swp requires %%u/%%u/%%u format");
         swp.cost_gapright = swp.cost_gapleft;
      endarg()
      optarg("-lsrl")
         if
         (  2 != sscanf(thearg(), "%u/%u", &swp.left_skip, &swp.right_limit)
         )
         die(1, "-lsrl requires %%u/%%u format");
      endarg()
      uniarg("--noindel") g_indel_allowed = 0; endarg()
      uniarg("--matrix") g_dump_matrix = 1; endarg()
      uniarg("--key-value") g_mode_print |= PRINT_KEY_VALUE; endarg()
      uniarg("--excise") g_mode_print |= PRINT_EXCISE; endarg()
      optarg("-do")  limit  = atoi(thearg());  endarg()
      optarg("-qs")  singlequery  = warn_to_ACGT(thearg(), strlen(thearg()));  endarg()
      optarg("-rs")  singleref    = warn_to_ACGT(thearg(), strlen(thearg()));  endarg()
      optarg("-identity")       theidentity  = atoi(thearg()); endarg()
      uniarg("--grep")  g_grep |= 1; endarg()
      optarg("-grep-format")  g_format_grep = thearg(); endarg()
      uniarg("--grepv")  g_grep |= 2; endarg()
      uniarg("--grepv-identifiers") g_grepv_format = GREPV_FORMAT_IDLIST; endarg()
      optarg("-grepv-o")      g_fnout_grepv = thearg(); endarg()
      optarg("-index") index_kmer = atoi(thearg()); if (index_kmer > 15) die(1, "index cannot exceed 15"); endarg()
      optarg("-concordance") index_kmer = atoi(thearg()); B_concordance = 1; if (index_kmer > 63) die(1, "concordance index cannot exceed 63"); endarg()
      optarg("-concordancex") index_kmer = atoi(thearg()); B_concordance = 3; if (index_kmer > 63) die(1, "concordance index cannot exceed 63"); endarg()
      optarg("-n-seeds") g_n_seeds = atoi(thearg()); endarg()
      optarg("-w-seeds") g_w_seeds = atoi(thearg()); endarg()
      optarg("-q")  fnquery  = thearg();  endarg()
      optarg("-r")  fnref   = thearg();  endarg()
      failarg()
      arg_done()

/* exit macromagicalitaciousness */

   ;  if (B_concordance && index_kmer >= 4 * sizeof(ktype))
      die(1, "This swan was built to index k-mers up to maximum size %d", 4 * sizeof(ktype) -1)

   ;  if (swp.flags & ( SW_NW_FILL | SW_NW_TRACE ))
      {  swp.flags |= SW_TRIM
      ;  swp.flags ^= SW_TRIM
   ;  }

      if (g_grep && !theidentity)
      die(1, "grep functionality requires setting an identity threshold")

   ;  if (singleref && fnquery)
      die(1, "Single reference and file query not supported. Swap reference and query or use -r <(echo -e \">dummy\\nACGTCGT\")")

   ;  if ((g_grep & 2) && !g_fnout_grepv)
      die(1, "--grepv requires specifying a file name with -grepv-o")

   ;  if (!fnquery && singlequery && strlen(singlequery) < index_kmer)
      {  argh("swan", "query length is very small; index disabled accordingly")
      ;  index_kmer = 0
   ;  }

      if (!(fpo = myfopen(g_fnout, "w", 0)))
      exit(1)
   ;  if (fnref && !(fpref = myfopen(fnref, "r", 1)))
      exit(1)
   ;  if (fnquery && !(fpquery = myfopen(fnquery, "r", 1)))
      exit(1)
   ;  if (g_fnout_grepv && !(fpo_grepv = myfopen(g_fnout_grepv, "w", 0)))
      exit(1)

   ;  fprintf(stderr, "Identity set to %d, index to %d\n", (int) theidentity, (int) index_kmer)

   ;  if (fpref)
      {  unsigned nref  = 0, nquery = 0
      ;  int i, overhang = 0, n_found = 0, best_i = 0
      ;  struct req* ls_ref = read_fasta_file(fpref, &nref, g_string_ref, g_rlen)
      ;  struct req* ls_query
      ;  struct req* q
      ;  struct req  query_single[2] = { { 0 }, { 0 } }

      ;  unsigned* todo = myalloc(nref * sizeof todo[0])
      ;  unsigned n_todo = 0
      ;  unsigned long n_index_items = 0;

      ;  struct kmer_index ref_index = { 0 }

      ;  if (!index_kmer)
         {  n_todo = nref
         ;  for (i=0; i<n_todo; i++)
            todo[i] = i
      ;  }

                                          /* some scary stuff going on here.
                                           * if B_concordance, then some big data structures are not built.
                                           * below works always, but no index is made when k == 0.
                                          */
         n_index_items = make_index(ls_ref, &ref_index, index_kmer, !B_concordance)

      ;  if (B_concordance)
         {  unsigned long xi, offset = 0, stretch = 1, j                              /* note stretch = 1 ... */
         ;  int multi_only = (B_concordance & 2)
				 ;	char buf[101]
         ;  struct kmer_x_req* kxr = ref_index.kxr
         ;  for (xi=1; xi<=n_index_items; xi++)                                       /* .. combines with xi=1 start */
            {  int end = (xi == n_index_items)                                        /* our index runs over the array size .. */
            ;  int change = !end && kxr[xi].kmer != kxr[xi-1].kmer                    /* .. so we need && shortcut here; chicanery */

            ;  if (change || end)
               {  int nq
               ;  qsort(kxr+offset, stretch, sizeof kxr[0], kxr_cmp_req)
               ;  nq = destructive_quash(kxr+offset, stretch)
               ;  if (!multi_only || nq > 1)
                  {  fprintf(fpo, "%s\t", get_sylmer(index_kmer, kxr[offset].kmer, buf));
                  ;  for (j=0; j<nq; j++)
                     fprintf(fpo, "%s%s", (j ? "," : ""), ls_ref[kxr[offset+j].req_index].annot)
                  ;  putc('\n', fpo)
               ;  }
                  offset = xi
               ;  stretch = 1
            ;  }
               else
               stretch++
         ;  }
            fclose(fpo)
         ;  exit(0)
      ;  }

         if (singlequery)
         {  query_single[0].s = stringle("%s", singlequery)
         ;  query_single[0].annot = stringle("(query)")
         ;  ls_query = query_single
         ;  nquery = 1
      ;  }
         else if (fpquery)
         ls_query = read_fasta_file(fpquery, &nquery, g_string_query, g_qlen)
      ;  else
         ls_query = ls_ref

                  /* hierverder. In the middle of shit */
      ;  q = ls_query
      ;  while (q->s)
                              /* within search_best_match we may also call do_align, namely this case:
                               * if (id_threshold && ai.aln_identity >= id_threshold)
                              */
         {  int recno_q = q - ls_query + 1
         ;  best_i = search_best_match(fpo, &ref_index, &swp, ls_ref, q, &overhang, todo, n_todo, theidentity, &n_found, recno_q)
;if(0)fprintf(stderr, "%d best %d found %s\n", (int) best_i, (int) n_found, q->s)
;if (best_i >= 0 && !n_found) die(1, "impossible (%d)", best_i)
         ;  if (g_grep)
            {  if (n_found && (g_grep & GREP_MATCH))
               do_align(fpo, ls_ref[best_i].s, ls_ref[best_i].annot, q->s, q->annot, best_i+1, recno_q, &swp, g_format_grep)
            ;  else if (!n_found && (g_grep & GREP_NOMATCH))
               {  if (g_grepv_format == GREPV_FORMAT_FASTA)
                  fprintf(fpo_grepv, ">%s\n%s\n", q->annot, q->s)
               ;  else
                  fprintf(fpo_grepv, "%s\n", q->annot)
            ;  }
            }
            else if (best_i >= 0 && !theidentity)
            {  do_align(fpo, ls_ref[best_i].s, ls_ref[best_i].annot, q->s, q->annot, best_i+1, recno_q, &swp, NULL)
         ;  }
            q++
         ;  if (limit && q-ls_query >= limit)
            break
      ;  }
         free(todo)
      ;  if (index_kmer)
         argh("swan", "Aligned %u pairs skipped %u deficient seed matches", g_n_alignments, g_n_seed_skipped)
   ;  }                             /* hierverder: no fpquery + single reference apparently */
      else if (singleref)
      do_align(fpo, singleref, "reference", singlequery ? singlequery : singleref, "query", 1, 1, &swp, NULL)
   ;  else
      die(1, "no reference!")

   ;  fclose(fpo)
   ;  if (fpo_grepv)
      fclose(fpo_grepv)

   ;  return 0
;  }








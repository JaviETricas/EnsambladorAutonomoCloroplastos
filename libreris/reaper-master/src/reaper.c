

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

/*

This is reaper, a flexible and fast program for pre-processing short-read sequencing data.

-  Match read with 3p-adaptor
   -  uses Smith-Waterman local alignment over full read and adaptor
      1. The best local alignment is considered and tested against user-set thresholds
      2. Optionally weaker matches between read-suffix and adaptor-prefix can be considered as well
   -  read/adaptor match is user controlled (alignment identity threshold)
   -  handling of poly-A (optional)
   -  handling of N in sequence (optional)
   -  handling of B quality scores (optional)

-  Check for ligation contamination
   -  Smith-Waterman local alignment

-  Handles barcodes
   -  Output split into per-barcode clean files and annotated discard files
   -  Recognises different layouts of barcodes and sequence inserts

-  Sane defaults and adaptable.
   -  Outputs per-barcode summary statistics for comprehensive QC plots
   -  Output fields can be chosen and formatted
   -  Per-read characteristics available (e.g. alignment, trinucleotide score)

-  One read-at-time, written in C
   -  Small memory footprint (order of megabytes)
   -  1-4 million reads per minute (depending on read length, contamination check)
   -  Twinned with 'tally' for fast uniquifying of reads


NOTES

*  Implementation caveat I.
   param->buck[buck_n] exists, contrary to normal C convention.
   qc_tally_in and qc_tally_out write in this bucket - it is used for
   grand-total reporting.  See /specialbucket/.
   This has the advantage that normal bucket processing code looks normal, and
   will not screw up the special bucket - only memory management looks peculiar.

*  Implementation caveat II.
   The code mixes 0-based ofsets and 1-based ofsets. The first derive
   from normal C-strings, the second from alignment matrices -- this
   is because they are zero-bordered; the first base in e.g. the adaptor
   maps to the second column of the alignment matrix, indexed 1.
   See /stringoffsets/.

*  Implementation Note.
   do_adaptor_3p and do_barcode_3p are quite different. The former allows the
   absence of a match, the latter requires exactly one match. The former allows
   more laxness when having an end-to-start match, and additionally allows a
   3p-adaptor prefix match if it is folllowed by low complexity code.


*/


#define REAPER_FAIL 29
#define REAPER_OK   0


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <signal.h>
#include <limits.h>
#include <time.h>
#include <zlib.h>
#include <ctype.h>


#include "reaper.h"
#include "version.h"
#include "sw.h"
#include "slib.h"
#include "trint.h"
#include "table.h"
#include "dna.h"


#ifndef MAXFIELDSIZE
#  define MAXFIELDSIZE  16384
#endif


int g_debug = 0;
int g_guard = 0;
int g_restrict = 18;

unsigned long g_fieldtosmall = 0;

                     /* NOTE these are coupled to the below array */
enum program_mode
{  CASE_NOTSPECIFIED = 0
                     /* Caveat: we use the test g_mode > CASE_NOBC in the program */
,  CASE_NOBC
,  CASE_3P_BC
,  CASE_5P_BC
}  ;

                     /* NOTE these are coupled to the above enum */
const char* G_MODE[6]
=  {  "NA"

   ,  "no-bc"
   ,  "3p-bc"
   ,  "5p-bc"
   };

enum program_mode g_mode = CASE_NOTSPECIFIED;


const char* x_3pgb[2] = { "3p-global -", "3p-global +" };
const char* x_tabu[2] = { "tabu -", "tabu +" };
const char* x_3ppx[2] = { "3p-prefix -", "3p-prefix +" };
const char* x_3pbc[2] = { "3p-barcode -", "3p-barcode +" };
const char* x_5pbc[2] = { "5p-barcode -", "5p-barcode +" };
const char* x_5psi[2] = { "3p-sinsert -", "3p-sinsert +" };



int dna[5] = { 'A', 'C', 'G', 'T', 'N' };


#define MATRIX_SIZE 8192


#define MAXQUALITY 256

#define dim unsigned

struct qc_composition
{  dim A[MAXFIELDSIZE]
;  dim C[MAXFIELDSIZE]
;  dim G[MAXFIELDSIZE]
;  dim T[MAXFIELDSIZE]
;  dim N[MAXFIELDSIZE]
;
}  ;


struct atabu
{  char*       tabuseq
;  unsigned    tabulen
;  struct match_requirement mr_tabu
;  struct atabu* next
;
}  ;


         /* Initially only used for barcode funneling, but can also be used for
          * other purposes.  This bucket design with multiple purposes does not
          * yet have a transparent design: the last bucket in struct param is
          * the only non-barcode bucket, and is special.
          * -  barcodes_io_open only opens streams for the barcode buckets.
         */
struct bucket
{  const char* barcode3p
;  const char* sinsert3p
;  const char* adaptor3p

;  const char* barcode5p
;  const char* sinsert5p

;  const char* cat3p    /* new */
;  int sinsert3p_n      /* new */
;  int sinsert5p_n      /* new */

;  const char* filetag  /* use instead of barcode in file name */

;  struct atabu tabu2   /* fixme, linked list with first link not freeable, designtabu
                         * The '2' in tabu2 has no meaning -- artefact from rewrite,
                         * kept around because it makes searching easier.
                         * Too lazy to think of something better.
                        */

;  int barcode_n
;  int barcode_start
;  int barcode_end

;  int cache_lint_offset
;  int cache_score
;  int cache_accept

;  char* fname

;  struct qc_composition qc_cpsin      /* input read composition */
;  struct qc_composition qc_cpsout     /* clean read composition */

;  unsigned global_in_max_size
;  unsigned global_out_max_size

;  dim qc_lenout[MAXFIELDSIZE+1]       /* clean read length */
;  dim * qc_bcq                        /* base call quality, one block is one read offet
                                        *    read-offset 0:    q0, q1, q2 ..... q255
                                        *    read-offset 1:    q0, q1, q2 ..... q255
                                        *       .....
                                        *    read-offset N-1:  q0, q1, q2 ..... q255
                                       */
;  ZFILE fp
;
}  ;


struct match_account
{  unsigned match_type        /* n(one) g(lobal) p(refix), t(ail) */
;  double suffix_complexity   /* for prefix match */
;  unsigned n_perfect         /* head to tail length */
;  int nnn_offset
;  int bbb_offset
;  int aaa_offset
;  int qqq_offset
;  int dust_offset
;
}  ;


struct param
{  unsigned cleanlength
#define VB_ALIGN   1 << 0
#define VB_LINT    1 << 1
#define VB_CLEAN   1 << 2
;  unsigned verbose
#define MODE_QC_EARLY      1 << 0
#define MODE_QC_LATE       1 << 1
#define MODE_QC_FULL       1 << 3
#define MODE_FULLLENGTH    1 << 4
#define MODE_KEEPALL       1 << 5
#define MODE_DUSTPOLYA     1 << 6
#define MODE_STREAMCLEAN   1 << 7
;  unsigned modes
;  unsigned sample_lint
;  unsigned sample_clean

                           /* we kludge/fudge this a bit. buck[buck_n] will also be utilised
                            * contrary to common conventions. It's an 'overflow' bucket
                            * for all reads that could not be classified in the barcode cases,
                            * or just all reads in the case where no barcode is present.
                           */
;  struct bucket *buck
;  int  buck_n             /* will be -1 initalised, so that loops 0 <= buck_n work */

;  unsigned p_3p_tth       /* prefix match, minimal tail to head perfect match */
;  const char* fnmeta
;  ZFILE fpclean
;  ZFILE fplint
;  int zippit
;  int polya               /* strip poly-As of at least this size */
;  unsigned taildust       /* strip tail if dust-like criterion exceeds this */
;  unsigned taildustlate   /* after adapter matching, strip tail if dust-like criterion exceeds this */
;  unsigned taildust_bases
;  unsigned taildustlate_bases
;  int trintco             /* if != 0, remove reads with tri-nt score >= trintco */
;  int length_co
;  int q_winlen
;  int q_winofs
;  int q_readofs           /* skip first readofs bases */
;  int q_cutoff
;  const char* adaptor_3p
;  const char* sinsert_5p
;  const char* adaptor_test   /* experimental */
;  struct atabu the_tabu2
;  const char* basename
;  struct sw_param swp_3p
;  struct sw_param swp_5p
;  struct sw_param swp_tabu
;  struct sw_param swp_raptor                   /* reaper - adaptor testing   */

;  struct match_requirement mr_3p_global        /* used without clause        */
;  struct match_requirement mr_3p_prefix        /* only used if mr_minlen > 0 */
;  struct match_requirement mr_3p_barcode       /* used without clause        */
;  struct match_requirement mr_5p_require       /* used without clause        */
;  struct match_requirement mr_tabu             /* used without clause        */
;  struct match_requirement mr_5p_barcode

;  double suffix_complexity_max                 /* in conjunction with prefix match */
;  double loco_check

;  unsigned nnn_count
;  unsigned nnn_outof

;  int anchor_offset
;  int length_limit
;  const char* anchor_string

;  unsigned kmer                                /* to analyse distance to adapter */

;  const char* format_clean
;  const char* format_lint
;  struct table metadata
;
}  ;


void param_init
(  struct param * pam
)                                            /*  1. minlen
                                                  | 2. maxedit
                                                  |  | 3. maxgap
                                                  |  |  | 4. offset
                                                  |  |  |  |        */
   {  struct match_requirement mr_3p_global = { 14, 2, 1, 0 }
   ;  struct match_requirement mr_3p_prefix = {  8, 2, 0, 2 }
   ;  struct match_requirement mr_3p_barcode= {  0, 6, 1, 0 }
   ;  struct match_requirement mr_5p_require= {  0, 2, 1, 10}
   ;  struct match_requirement mr_tabu      = { 16, 2, 1, 0 }
   ;  struct match_requirement mr_5p_barcode= {  0, 0, 0, 2 }  /* bit field 1 = zip 2 = normal alignment */

   ;  struct atabu theboo           =  { "", 0, { 0 }, NULL }
   
   ;  pam->cleanlength              =  0
   ;  pam->verbose                  =  0
   ;  pam->modes                    =  MODE_QC_FULL
   ;  pam->sample_clean             =  1
   ;  pam->sample_lint              =  1
   ;  pam->polya                    =  0
   ;  pam->taildust                 =  0
   ;  pam->taildustlate             =  0
   ;  pam->taildust_bases           =  0
   ;  pam->taildustlate_bases       =  0

   ;  pam->p_3p_tth                 =  0   /* strip perfect (read) tail to head (3p adaptor) match of any length >= 1 */

   ;  pam->mr_3p_global             =  mr_3p_global
   ;  pam->mr_3p_prefix             =  mr_3p_prefix
   ;  pam->mr_3p_barcode            =  mr_3p_barcode
   ;  pam->mr_5p_require            =  mr_5p_require
   ;  pam->mr_tabu                  =  mr_tabu
   ;  pam->mr_5p_barcode            =  mr_5p_barcode

   ;  pam->buck                     =  NULL
   ;  pam->buck_n                   =  -1

   ;  pam->fnmeta                   =  NULL
   ;  pam->adaptor_3p               =  NULL
   ;  pam->adaptor_test             =  NULL
   ;  pam->sinsert_5p               =  NULL
   ;  pam->basename                 =  "out"
   ;  pam->trintco                  =  0
   ;  pam->length_co                =  0
   ;  pam->q_winlen                 =  0
   ;  pam->q_winofs                 =  0
   ;  pam->q_readofs                =  10 
   ;  pam->q_cutoff                 =  0
   ;  pam->zippit                   =
#if   WE_USE_ZLIB
   1
#else
   0
#endif

   ;  pam->swp_3p.cost_gapleft      =  3
   ;  pam->swp_3p.cost_gapright     =  3
   ;  pam->swp_3p.cost_subst        =  1
   ;  pam->swp_3p.gain_match        =  4
   ;  pam->swp_3p.left_skip         =  0
   ;  pam->swp_3p.right_limit       =  0
   ;  pam->swp_3p.flags             =  0

   ;  pam->swp_5p                   =  pam->swp_3p
   ;  pam->swp_tabu                 =  pam->swp_3p

   ;  pam->swp_raptor               =  pam->swp_3p
   ;  pam->swp_raptor.cost_gapleft  =  1
   ;  pam->swp_raptor.cost_gapright =  1
   ;  pam->swp_raptor.cost_subst    =  1
   ;  pam->swp_raptor.gain_match    =  2

   ;  pam->the_tabu2                =  theboo

   ;  pam->nnn_count                =  0
   ;  pam->nnn_outof                =  0

   ;  pam->anchor_offset            =  -1
   ;  pam->anchor_string            =  NULL

   ;  pam->kmer                     =  10

   ;  pam->suffix_complexity_max    =  0.25
   ;  pam->loco_check               =  0.0

   ;  pam->format_clean             =  "@%I%n%C%n+%n%Q%n"
   ;  pam->format_lint              =  "@%I%tmsg=%M%n%R%n+%n%Y%n"

   ;  pam->length_limit             =  0

   ;  pam->fpclean                  =  NULL
   ;  pam->fplint                   =  NULL

   ;  memset(&pam->metadata, 0, sizeof pam->metadata)
;  }



unsigned get_mdta       /* minimum distance to adapter, per k-mer varied along read */
(  char* buf_distance
,  char* buf_offset
,  const char* seq
,  unsigned seq_n
,  struct param* pam
)
   {  const char* theaptr = pam->adaptor_test
   ;  unsigned k = pam->kmer
   ;  struct sw_param* swp = &pam->swp_raptor
   ;  static unsigned char* mdta = NULL
   ;  unsigned mask = 0
   ;  unsigned kmer = 0
   ;  char seq2[MAXFIELDSIZE+1]
   ;  int i = 0
   ;  if (!mdta)
      mdta = calloc(1 << 2*k, 1)

   ;  buf_distance[0] = '\0'
   ;  buf_offset[0] = '\0'

   ;  strncpy(seq2, seq, MAXFIELDSIZE)

   ;  if (!theaptr)
      theaptr = pam->adaptor_3p

   ;  if (!theaptr)           /* give up */
      seq_n = 0

   ;  for (i=0; i < seq_n; i++)
      {  int dist = 7, offset = 31

      ;  unsigned b = themap[(unsigned char) seq2[i]]
      ;  mask = (mask << 1) & ((1 << k) -1)
      ;  kmer = (kmer << 2) & ((1 << 2*k) -1)
      ;  if (b > 3)
         mask |= 1
      ;  else
         kmer |= b
;if(0)fprintf(stderr, "[%s]\n", seq2+i)

      ;  if (i+1 < k)
         {  buf_distance[i] = '.'
         ;  buf_distance[seq_n+i+1] = '.'
         ;  continue
      ;  }
                                         /* store distance in 3 bits (0..7), offset in 5 bits (1..31) */
         if (mask || !mdta[kmer])
         {  SWNUM data[MATRIX_SIZE] = { 0 }
         ;  struct sw_alninfo ai = { 0 }
         ;  char bu = seq2[i+1]
         ;  seq2[i+1] = '\0'             /* sw_fill wants 0-terminated. blah. */
         ;  if
            ( !sw_fill
               (  &ai
               ,  data
               ,  MATRIX_SIZE
               ,  seq2 + i+1-k
               ,  theaptr
               ,  swp->cost_gapleft + swp->cost_gapright
               ,  swp
            )  )
            {  sw_trace_trim(&ai, swp, ai.max_ij, 0)         /* do_barcode_3p */
            ;  dist = k - ai.n_match + ai.n_insr + ai.n_insl
            ;  offset = ai.rgt_start
         ;  }
            if (dist > 7)     dist = 7
         ;  if (offset > 31)  offset = 31
         ;  if (!mask)
            mdta[kmer] = dist | (offset << 3)
         ;  seq2[i+1] = bu
      ;  }
         else
            dist = mdta[kmer] & 7
         ,  offset = mdta[kmer] >> 3

      ;  buf_distance[i] = '0' + (dist > 9 ? 9 : dist)
      ;  buf_distance[seq_n+i+1] = offset["0123456789ABCDEFGHIJKLMNOPQRSTUV"]
   ;  }


   ;  if (seq_n)
      {  buf_distance[i] = '\n'
      ;  buf_distance[2 * seq_n+1] = '\0'
      ;  buf_distance[0] = 'D'
      ;  buf_distance[seq_n+1] = 'O'
   ;  }
      else
      {  memcpy(buf_distance, "D\nO\0", 4)
      ;  return 3
   ;  }
      return 2 * seq_n + 1
;  }



struct read_stats_global
{  unsigned N
;  unsigned N_err       /* globally set */
;  unsigned N_truncated    /* */
;  unsigned N_alien     /* not ACGTUNX */
;  unsigned N_cflr      /* conflict resolved, per case */
;  unsigned N_clean     /* per case */

;  unsigned N_bb        /* globally set */
;  unsigned N_qq        /* globally set */
;  unsigned N_anchor    /* globally set */
;  unsigned N_loco      /* globally set */
;  unsigned N_nomsg     /* globally set */
;  unsigned N_cfl       /* per case */
;  unsigned N_nn        /* per case */
;  unsigned N_aa        /* per case */
;  unsigned N_tri       /* per case */
;  unsigned N_nomatch   /* per case */
;  unsigned N_lenlo     /* per case */
;  unsigned N_tabu      /* per case */
;  unsigned N_5p_nosinsert /* per case */          /* This group should add up to */

;  unsigned N_lint      /* per case */          /* to this */

;  unsigned long B      /* */
;  unsigned long N_cells  /*    total size of alignment matrices */
;  unsigned B_BB        /* bases lost due to BBB, globally set */
;  unsigned B_tri       /* bases lost due to tri-tail trimming, globally set */
;  unsigned B_NN        /* bases lost due to NNN, globally set */
;  unsigned B_AA        /* bases lost due to AAA, (after adaptor stripping) */
;  unsigned B_QQ        /* bases lost due to QQQ, (after adaptor stripping) */
;  unsigned B_lco       /* bases lost due to user specified length cutoff */

/* hierverder; implement below */
;  unsigned B_3pA       /* bases lost due to 3pA stripping */
;  unsigned B_5pMM      /* mismatches in required 5p adaptor */
;
}  ;


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


struct record
{  char id        [MAXFIELDSIZE+1]
;  char seq       [MAXFIELDSIZE+1]
;  char q         [MAXFIELDSIZE+1]
;  char discard   [MAXFIELDSIZE+1]
;  char annot     [MAXFIELDSIZE+1]

;  unsigned id_n
;  unsigned seq_n
;  unsigned q_n
;  unsigned discard_n
;  unsigned annot_n

;  int clean_n

;  unsigned count
;  const char* out_message
;  struct read_stats_global rs            /* NOTE maintains state between different read_record2 calls */
;
}  ;



void validate_sequence
(  struct record* rec
)
   {  int i
   ;  for (i=0;i<rec->seq_n;i++)
      {  rec->seq[i] = toupper((unsigned char) rec->seq[i])
      ;  if (BASEMAP((unsigned char) rec->seq[i]) > 4)
            rec->rs.N_alien++
         ,  rec->seq[i] = 'N'
   ;  }
   }


#define RR_ERROR 3
#define RR_NOMEM 2
#define RR_DONE  1
#define RR_OK    0


         /* At the moment this keeps reading in case of RL_NOMEM,
          * but readline does not flush overly long lines yet,
          * so there is not yet a mechanism to just ignore
          * a record with lines that are too long.
          * todo fixme
         */
int read_record2
(  ZFILE  ip
,  const char* format
,  struct record* rec
,  int length_limit
)
   {  unsigned char* a = (unsigned char*) format
   ;  unsigned rlstat = 0
   ;  int q_received =  0
   ;  unsigned truncated  =  0
   ;  int maxfieldsize = MAXFIELDSIZE

   ;  if (length_limit && length_limit < MAXFIELDSIZE)
      maxfieldsize = length_limit

   ;  rec->id_n      =  0
   ;  rec->seq_n     =  0
   ;  rec->annot_n   =  0
   ;  rec->discard_n =  0
   ;  rec->clean_n   =  0
   ;  rec->q_n       =  0
   ;  rec->count     =  1
   ;  rec->id[0]     =  '\0'
   ;  rec->seq[0]    =  '\0'
   ;  rec->discard[0]=  '\0'
   ;  rec->q[0]      =  '\0'
   ;  rec->annot[0]  =  '\0'

   ;  while (a[0])
      {  char* receive = NULL
      ;  unsigned *np_receive = 0
      ;  switch(a[0])
         {  case 'I':   receive = rec->id;      np_receive = &rec->id_n;        break
         ;  case 'R':   receive = rec->seq;     np_receive = &rec->seq_n;       break
         ;  case 'D':   receive = rec->discard; np_receive = &rec->discard_n;   break
         ;  case 'A':   receive = rec->annot;   np_receive = &rec->annot_n;     break
         ;  case 'Q':   receive = rec->q;       np_receive = &rec->q_n; q_received = 0; break
      ;  }

         rlstat |= kraken_readline(ip, receive, maxfieldsize, np_receive, &truncated)
      ;  if (a[0] == 'R')
         rec->rs.N_truncated += truncated

      ;  if (rlstat & (RL_DONE | RL_ERROR))
         break
      ;  a++
   ;  }

      if (rlstat & (RL_DONE | RL_ERROR))
      {  if (a != (unsigned char*) format && a[0])
         argh("reaper", "last record was incomplete")
      ;  return RR_DONE
   ;  }

      if (rlstat & RL_NOMEM)           /* fixme not checked yet */
      return RR_NOMEM

   ;  rec->rs.N++

   ;  validate_sequence(rec)

   ;  if (q_received && rec->q_n != rec->seq_n)
      {  rec->rs.N_err++
      ;  rec->q_n = 0
      ;  rec->q[0] = '\0'
   ;  }

      return RR_OK
;  }


                                       /*
                                          ;  case 'I':   receive = rec->id;      np_receive = &rec->id_n;        break
                                          ;  case 'R':   receive = rec->seq;     np_receive = &rec->seq_n;       break
                                          ;  case 'D':   receive = rec->discard; np_receive = &rec->discard_n;   break
                                          ;  case 'Q':   receive = rec->q;       np_receive = &rec->q_n; q_received = 0; break
                                       */

#define izblank(c) ((unsigned char) (c) == ' ' || (unsigned char) (c) == '\t')

int read_recordx
(  ZFILE  ip
,  const char* format
,  struct record* rec
)
   {  const char* fmtp = format, *fmtz = format + strlen(format)
   ;  unsigned rlstat = 0
#define THE_BUFSIZE 8191
   ;  char buf[THE_BUFSIZE+1]
   ;  unsigned n_received = 0, n_lines = 0
   ;  unsigned truncated = 0
   ;  char* bufp = NULL, *curp = NULL, *bufz = NULL
   ;  int esc = 0

   ;  rec->id_n      =  0
   ;  rec->seq_n     =  0
   ;  rec->discard_n =  0
   ;  rec->q_n       =  0
   ;  rec->clean_n   =  0
   ;  rec->count     =  1
   ;  rec->id[0]     =  '\0'
   ;  rec->seq[0]    =  '\0'
   ;  rec->discard[0]=  '\0'
   ;  rec->q[0]      =  '\0'

   ;  while (fmtp < fmtz)
      {  if (!bufp)
         {  rlstat |= kraken_readline(ip, buf, THE_BUFSIZE, &n_received, &truncated)
         ;  if (truncated || rlstat & (RL_DONE | RL_ERROR))
            break
         ;  bufp = buf
         ;  bufz = buf + n_received
         ;  n_lines++
         ;  esc  = 0
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
         case 'C': case 'X':
                  {  int n_scanned = 0
                  ;  if (sscanf(bufp, "%u%n", &rec->count, &n_scanned) < 1)
                     goto DONE
                  ;  bufp += n_scanned
               ;  }
                  break
               ;
            case 'G':
                  while (bufp < bufz && !izblank((unsigned char) bufp[0]))
                  bufp++
               ;  rec->discard_n = cpytofield(rec->discard, curp, bufp -curp)
               ;  break
               ;
            case 'F':
                  while (bufp < bufz && '\t' != (unsigned char) bufp[0])
                  bufp++
               ;  rec->discard_n = cpytofield(rec->discard, curp, bufp -curp)
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
            case 'a':
                  rec->annot_n = cpytofield(rec->annot, rec->discard, rec->discard_n)
               ;  break
               ;
            case 'A':
                  while (izblank((unsigned char) curp[0]))
                  curp++
               ;  bufp = curp
               ;  while (bufp < bufz)
                  bufp++
               ;  rec->annot_n = cpytofield(rec->annot, curp, bufp -curp)
               ;  break
               ;
            case 'I':
                  while (bufp < bufz && !izblank((unsigned char) bufp[0]))
                  bufp++
               ;  rec->id_n = cpytofield(rec->id, curp, bufp -curp)
               ;  break
               ;
            case 'Q':
                  while (bufp < bufz && !izblank((unsigned char) bufp[0]))
                  bufp++
               ;  rec->q_n = cpytofield(rec->q, curp, bufp -curp)
               ;  break
               ;
            case 'R':
                  while (bufp < bufz && isalpha((unsigned char) bufp[0]))
                  bufp++
               ;  rec->seq_n = cpytofield(rec->seq, curp, bufp -curp)
               ;  validate_sequence(rec)
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
      {  argh("reaper", "parse error (remaining format string %s, buffer [%s])", fmtp, bufp)
      ;  return RR_ERROR
   ;  }

      if (rlstat & RL_DONE)
      return RR_DONE

   ;  if (rlstat & RL_ERROR)
      return RR_ERROR

   ;  if (rlstat & RL_NOMEM)           /* fixme not checked yet */
      return RR_NOMEM

   ;  rec->rs.N++

   ;  return RR_OK
;  }




                                                         /*
/stringoffsets/
=========================================================
alnspace

   alignment space:
               0 1 2 3 4 5 6 7 8 9 10        coordinates
                 G G G A T C C T A G         read
                             C T A G G C C   adaptor
left_start  : 7
lint_offset : left_start - 1 == 6 == clean_n

=========================================================

   string space:
               0 1 2 3 4 5 6 7 8 9 10        coordinates
               G G G A N C N T N G           read
                       ^                     N offset
 
N offset is 4, we return 5, so that
lint_offset == NNN_offset - 1 == 4 == clean_n

========================================================= */



int get_nnn_offset
(  const char* s
,  unsigned n
,  unsigned count
,  unsigned outof
)
   {  unsigned hist[SW_ALN_MAXSIZE] = { 0 }

   ;  if (!outof || count > outof)
      return -1

   ;  int i, nn = 0
   ;  for (i=0;i<n;i++)
      {  unsigned char a = s[i]
      ;  unsigned char e = hist[i % outof]
      ;  hist[i % outof] = a
      ;  nn += (a == 'N') - (e == 'N')
      ;  if (nn >= count)
         {  while (nn > 0)     /* now check which ones are not N */
            {  if (hist[i-- % outof] == 'N')
               nn--
         ;  }
            return i+2      /* noteme, alignment space coordinate
                             * i+1 is the string-space NNNN offset,
                             * but we need to return offset in aligment space (see above)
                            */
      ;  }
      }
      return -1
;  }


                                                /* same return semantics as get_nnn_offset */
void get_dust_info
(  const char* s
,  unsigned n
,  int* maxp                                    /* dust score for tail */
,  int* maxidp                                  /* offset where this is attained */
,  int* dominant_base                           /* dominant base */
,  int* dominant_base_count                     /* and its count */
)
   {  int max = 0, maxid =0, i
   ;  unsigned count[6] = { 0 }

   ;  dustscore_tail(s, n, NULL, &max, &maxid)

   ;  maxp[0]     =  max
   ;  maxidp[0]   =  maxid

   ;  dominant_base[0]       = 4
   ;  dominant_base_count[0] = 0

   ;  for (i=maxid; i<n; i++)
      count[themap[(unsigned char) s[i]]]++

   ;  for (i=0;i<4;i++)
      if (count[i] > dominant_base_count[0])
         dominant_base[0] = i
      ,  dominant_base_count[0] = count[i]
;  }

                                                /* same return semantics as get_nnn_offset */
int get_dust_offset
(  const char* s
,  unsigned n
,  unsigned threshold
,  unsigned bases
)
   {  int max = 0, maxid =0, found = 0
   ;  dustscore_tail(s, n, NULL, &max, &maxid)
   ;  found = threshold && max >= threshold && maxid+6 <= n

   ;  if (found && bases)
      {  int i
      ;  unsigned count[6] = { 0 }
      ;  for (i=maxid; i<n; i++)                    /* fixme: i=0 should that be i=maxid ? */
         count[themap[(unsigned char) s[i]]]++
      ;  for (i=0;i<4;i++)
         {  if (count[i] * 2 > (n - maxid) && (bases & (1 << i)))
            break
      ;  }
         found = i < 4
   ;  }
      return found ? (maxid+1) : -1
;  }


                                                /* same return semantics as get_nnn_offset */
int get_polya_offset
(  const char* s
,  unsigned n
,  int minlen
)
   {  int i = n                                 /* Note: should work for q_n == 0 */
   ;  int polya = 0

   ;  while (--i >= 0 && s[i] == 'A')
      polya++
   ;  return minlen && polya >= minlen ? n + 1 - polya : -1
;  }

int cmp_int
(  const void* a
,  const void* b
)
   {  return ((int*)a)[0] - ((int*)b)[0] 
;  }


int get_qqq_offset
(  struct record* rec
,  int cut                                /* cut when median drops below this */
,  unsigned winlen
,  unsigned winofs
,  unsigned readofs
)
   {  unsigned char hist[256] = { 0 }
   ;  int window[255]                     /* signed so that we can use a-b for comparison */
   ;  unsigned winmid
   ;  int i = 0, n_lt = 0, n_gt = 0
   ;  int median = 0
   ;  int lastqual = 0

   ;  if (!winlen || readofs + winlen > rec->q_n)
      return -1

   ;  lastqual = rec->q[rec->q_n-1]

   ;  winmid = (winlen+1) / 2             /* e.g. 5 -> 3; n_smaller + n_equal must be >= winmid 
                                           * and also n_larger + n_equal must be >= winmid
                                          */
   ;  for (i=readofs;i<winlen+readofs; i++)
      window[i % winlen] = rec->q[i]

   ;  qsort(window, winlen, sizeof window[0], cmp_int)
   ;  median = window[winmid-1]

   ;  for (i=readofs;i<winlen+readofs; i++)
      {  unsigned i_out = i % winlen
      ;  window[i_out] = rec->q[i]
      ;  hist[(unsigned char) window[i_out]]++
      ;  n_lt += window[i_out] < median
      ;  n_gt += window[i_out] > median
   ;  }

      i = readofs + winlen - 1                                 /* necessary for median >= cut logic */

   ;  while (median >= cut && ++i < rec->q_n + winlen)
      {  unsigned i_out = i % winlen
      ;  unsigned rgt = i < rec->q_n ? rec->q[i] : lastqual    /* shift this value in */
      ;  unsigned lft = window[i_out]                          /* shift this value out */
      ;  window[i_out] = rgt

      ;  hist[rgt]++
      ;  hist[lft]--

      ;  if (rgt < median) n_lt++
      ;  if (rgt > median) n_gt++
      ;  if (lft < median) n_lt--
      ;  if (lft > median) n_gt--

      ;  while (!hist[median] || n_lt + hist[median] < winmid)
         {  n_lt += hist[median]
         ;  n_gt -= hist[median+1]
         ;  median++
      ;  }
         while (!hist[median] || n_gt + hist[median] < winmid)
         {  n_gt += hist[median]
         ;  n_lt -= hist[median-1]
         ;  median--
      ;  }
      }
      {  int j = (i+1+winofs)-winlen
      ;  if (j > rec->q_n)
         j = rec->q_n
      ;  if (i < (rec->q_n + winlen))
         return j
   ;  }
      return -1
;  }



int get_bbb_offset
(  struct record* rec
)
   {  int i = rec->q_n                          /* Note: should work for q_n == 0 */

   ;  while (--i >= 0)
      {  if (rec->q[i] != 'B')
         break
   ;  }
      return i+1 == rec->q_n ? -1 : i+2         /* same return semantics as get_nnn_offset */
;  }


double string_complexity
(  const char* a
,  unsigned n
)
   {  int i, n_pairs = 0
   ;  if (!n)
      return 1.0
   ;  for (i=0;i+1<n;i++)
      n_pairs += a[i] == a[i+1]
   ;  return (1.0 - n_pairs * 1.0 / n)
;  }


#if 0
   /* fixme: n_stretch in callers often computed using SW_ALN_MAXSIZE, but trimming changes this */
   /* returns number of matches if parameters satisfied, 0 otherwise */
int alignment_edit_ok
(  const unsigned char* aln
,  unsigned aln_length
,  unsigned n_stretch
,  unsigned n_edit_max
,  unsigned n_gap_max
)
   {  int i, n_edit = 0, n_gap = 0, n_subst = 0, n_match = 0, ok = 0
   ;  if (!n_stretch)
      return 0
   ;  unsigned hist[SW_ALN_MAXSIZE] = { 0 }
;if(0)fprintf(stderr, "length %d n_stretch %d e %d\n", (int) aln_length ,  (int) n_stretch ,  (int) n_edit_max)
   ;  for (i=0;i<aln_length;i++)
      {  unsigned char a = aln[i]
      ;  unsigned char e = hist[i % n_stretch]
      ;  hist[i % n_stretch] = a
      ;  n_subst += (a == 'x') - (e == 'x')
      ;  n_gap   += (a == '-') - (e == '-')
      ,  n_edit   = n_subst + n_gap
      ;  n_match += a == '|'
      ;  if (  n_edit <= n_edit_max
            && i+1 + n_edit_max >= n_stretch + n_edit
            && n_gap <= n_gap_max
            )
         ok = 1
   ;  }
      return ok ? n_match : 0
;  }
#endif


                        /* the prefix match is complicated.
                         * We must be able to find a good match at the end
                         * that's not beaten by a longer shoddy match at start
                         * for which the suffix complexity criterion fails. This
                         * probably, ideally, requires using the last-row max_ej cell.
                         * The max_ej cell will not, however, find prefix matches
                         * followed by low complexity reads.
                         */

int sw_match_end
(  struct sw_alninfo* ai
,  struct match_account* ma
,  struct sw_param* swp
,  unsigned prefix_perfect_minlen
,  struct match_requirement* mr_prefix
)
   {  unsigned n_gap  = ai->n_insl + ai->n_insr
   ;  unsigned n_edit = n_gap + ai->n_subst
   ;  unsigned prefix_length = ai->rgt_end + 1 - ai->rgt_start

   ;  ma->suffix_complexity = 0.0

                        /* This code uses the globally highest scoring alignment,
                         * and applies the prefix match requirement to the entire
                         * alignment, not just a subpart.
                         * This is so that a case such as
GATAGATCGAAACCATCCTCTGCTAGTTTTC-GTATG : [16,36]
               ||x|x|||-x||x|||-|xx|| : score 56
               TCGTATGC-CGTCTTCTGCTTGT : [1,21] 0/0 recno=1
                         * is not considered a match.
                        */
   ;  unsigned prefix_match
      =  (  (!mr_prefix->mr_offset || (ai->rgt_start <= mr_prefix->mr_offset))
         && prefix_length >= mr_prefix->mr_minlen
         && n_edit <= mr_prefix->mr_maxedit
         && n_gap  <= mr_prefix->mr_maxgap
         )

                        /* now we just inspect the end. */
   ;  if (!prefix_match && prefix_perfect_minlen)
      {  int ofs = (ai->ni-1) * ai->nj  /* this sits on the left column of zeros, last entry */
      ;  unsigned j, bestj = 0
      ;  for (j=prefix_perfect_minlen;j<ai->nj;j++)
         {  if (ai->data[ofs+j]  == j * swp->gain_match)
            bestj = j
      ;  }
         ma->n_perfect = bestj
      ;  return bestj ? 1 : 0
   ;  }

                        /* NOTE. The prefix match may fall just short of the read,
                         * for example miss the last one or two bases. The suffix
                         * complexity will still be small, as a result of
                         * prefix_length (the matching part of the adaptor)
                         * being part of the formula
                        */
      if (prefix_match)
      {  if (ai->lft_end + 1 < ai->ni)
         {  double left_trail = ai->ni - ai->lft_end - 1
         ;  double c = string_complexity(ai->left+ai->lft_end+1, ai->ni - ai->lft_end-1)
         ;  ma->suffix_complexity = c * (left_trail) /  (1.0 * prefix_length + left_trail)
      ;  }
         return 1
   ;  }
      return 0
;  }


         /* barcode code */
int sw_match_bc3p
(  struct sw_alninfo* ai
,  struct match_requirement* mr_global
,  struct match_requirement* mr_prefix
,  struct match_requirement* mr_barcode
,  int  sinsert_n       /* number of bases before the barcode, not yet used. */
,  int  barcode_n       /* length of barcode, not yet used. */
,  int* n_match
)
   {  int n_length = mr_global->mr_minlen
   ;  int n_edit_max = mr_global->mr_maxedit
   ;  int n_gap_max  = mr_global->mr_maxgap
   ;  n_match[0] = 0
   ;
                        /* /barcodehardcoded/
                         * note: 3p-sinsert+barcode start hardcoded to be in pos 1 or 2
                         * (beware of alignment space, which uses offset 1)
                         * Now we have to be careful wrt to the barcode match requirements;
                         * How do we locate the right spot in the alignment?
                        */
      if (ai->rgt_start != 1 && ai->rgt_start != 2)
      {  if(0)fprintf(stderr, "bail-out 2 (%d)\n", (int) ai->rgt_start)
      ;  return 0
   ;  }
   
                              /* reduce global range alignment length parameter if match is at end of range */
      if (ai->ni - ai->lft_start < n_length && mr_prefix->mr_minlen)
      {  n_length = ai->ni - ai->lft_start
      ;  if (n_length < mr_prefix->mr_minlen)
         return 0
      ;  n_edit_max = mr_prefix->mr_maxedit
      ;  n_gap_max = mr_prefix->mr_maxgap
   ;  }

      return
      (  (  !mr_barcode->mr_minlen     /* check barcode range alignment if specied */
         || sw_alignment_edit_ok
            (  ai->aln_aln + ai->aln_ofs
                                       /* fixme, better control over barcode */
            ,  mr_barcode->mr_minlen   /* fixed length search: can this overrun the alignment? document! */
            ,  mr_barcode->mr_minlen
            ,  mr_barcode->mr_maxedit
            ,  mr_barcode->mr_maxgap
            )
         )
      &&
                              /* second argument, alignment length: this is not
                               * strictly correct, but the remainder of the
                               * alignment is zero-filled.  Note: our search
                               * space (right_limit in sw_fill) has been
                               * restricted in do_barcode_3p.
                              */

#define WOTCHA 0
                              /* check global range alignment */
         (  n_match[0]
         =  sw_alignment_edit_ok
            (  ai->aln_aln + ai->aln_ofs
            ,  ai->aln_end - ai->aln_ofs
            ,  n_length
            ,  n_edit_max
            ,  n_gap_max
         )  )
      )
;  }


/* aln_ofs == SW_ALN_MAXSIZE + 1 */


int sw_match_any
(  struct sw_alninfo* ai
,  struct match_requirement* mr
)
   {
      return
         (  !mr->mr_offset
         || (mr->mr_offset > 0 && ai->rgt_start <= mr->mr_offset)
         || (mr->mr_offset < 0 && (ai->nj - ai->rgt_end) <= -mr->mr_offset)
         )
      && sw_alignment_edit_ok
         (  ai->aln_aln + ai->aln_ofs
         ,  ai->aln_end - ai->aln_ofs
         ,  mr->mr_minlen
         ,  mr->mr_maxedit
         ,  mr->mr_maxgap
         )
;  }


int sw_printaccount
(  const struct match_account* ma
,  char* buf
,  unsigned bufsize
)
   {  return
      snprintf
      (  buf
      ,  bufsize
      ,  "mt=%c,sc=%.2f,ht=%d,nn=%d,bb=%d,aa=%d,qq=%d"
      ,  (int) ma->match_type
      ,        ma->suffix_complexity
      ,  (int) ma->n_perfect
      ,  (int) ma->nnn_offset
      ,  (int) ma->bbb_offset
      ,  (int) ma->aaa_offset
      ,  (int) ma->qqq_offset
      )
;  }


int sw_printaln3
(  const struct sw_alninfo* ai
,  char* buf
,  unsigned bufsize
)
   {  unsigned o = ai->aln_ofs
   ;  return
      snprintf
      (  buf
      ,  bufsize
      ,  "ls=%d,le=%d,rs=%d,re=%d,score=%d,left=%s,aln=%s,right=%s"
      ,  (int) ai->lft_start
      ,  (int) ai->lft_end
      ,  (int) ai->rgt_start
      ,  (int) ai->rgt_end
      ,  (int) ai->data[ai->max_ij]
      ,  ai->aln_lft+o
      ,  ai->aln_aln+o
      ,  ai->aln_rgt+o
      )
;  }


void dump_read
(  int cleaned
,  void* fpo
,  const struct record* rec
,  const struct sw_alninfo* ai
,  const struct match_account* ma
,  struct param* pam
)
#define THEBUFSIZE 1023
   {  char buf[THEBUFSIZE+1]
   ;  char buf2[THEBUFSIZE+1]
#if WE_USE_ZLIB
   ;  int zippit = pam->zippit
#endif
   ;  const char* format = cleaned ? pam->format_clean : pam->format_lint
   ;  const char* f = format
   ;  const char* z = format + strlen(format)
   ;  char dust [MAXFIELDSIZE+1]
   ;  int escape = 0, tnt
   ;  int dustidx = 0, dustscore = 0

   ;  while (f < z)
      {  unsigned c = (unsigned char) f[0]
      ;  const char* p = buf
      ;  int n = 0
      ;  if (!escape)
         {  if (c == '%')
            escape = 1
         ;  else
               buf[0] = c
            ,  n = 1
      ;  }
         else
         {  int cleanlen = cleaned ? rec->clean_n : rec->seq_n
         ;  switch(c)
            {  case 'R': p = rec->seq; n = rec->seq_n;
break;
            ;  case 'C': case 'E': p = rec->seq; n = rec->clean_n;
                                 if (c == 'E' && n == 0) { p = "N"; n = 1; }
break;
            ;  case 'V': revcompl(rec->seq, rec->clean_n, buf); p = buf; n = rec->clean_n;
break;
            ;  case 'Z': { int themax = rec->seq_n;
                           if (pam->length_co > themax) themax = pam->length_co;
                           n = snprintf(buf, THEBUFSIZE+1, "%.*s", (int) rec->clean_n, rec->seq);
                           while (n < themax && n < THEBUFSIZE)
                              buf[n++] = 'N';
                           buf[n] = '\0';
break;
                         }
            ;  case 'I': p = rec->id; n  = rec->id_n; break;
            ;  case 'L': n = snprintf(buf, THEBUFSIZE+1, "%d", (int) rec->clean_n); break;
            ;  case 'X': n = snprintf(buf, THEBUFSIZE+1, "%u", (unsigned) rec->count); break;
            ;  case 'Q': p = rec->q; n = cleanlen; break;
            ;  case 'Y': p = rec->q; n = rec->q_n; break;
            ;  case 'U': {  int score, score_offset, base, base_count
                         ;  int seqlen = cleanlen
                         ;  get_dust_info(rec->seq, seqlen, &score, &score_offset, &base, &base_count)
                         ;  n = snprintf(buf, THEBUFSIZE+1, "score=%d,offset=%d,len=%d,base=%c,count=%d",
                                    score, score_offset+1, seqlen - score_offset, dna[base], base_count)
                       ; }
;  break
            ;  case 'T':   tnt = trintscore(rec->seq, cleanlen)
                        ;  n = snprintf(buf, THEBUFSIZE+1, "%d", tnt)
;  break
            ;  case 'D':   dustscore_tail(rec->seq, rec->clean_n, dust, &dustscore, &dustidx)
;  p = dust; n = rec->clean_n; break
            ;  case '_':  if (!dustscore)
                            dustscore_tail(rec->seq, rec->seq_n, NULL, &dustscore, &dustidx)
                        ;   n = snprintf(buf, THEBUFSIZE+1, "%d:%d", dustidx+1, dustscore);
break
            ;  case 'M': p = rec->out_message; n = strlen(p); break
            ;  case 'i': case 'J': n = snprintf(buf, THEBUFSIZE+1, "%d", (int) rec->rs.N); break;
            ;  case 'f': n = snprintf(buf, THEBUFSIZE+1, "%d", (int) (4 * (rec->rs.N-1) + 2)); break;
            ;  case '3': if (ai) n = sw_printaln3(ai, buf, THEBUFSIZE+1);  break
            ;  case 'A': p = rec->annot; n = rec->annot_n; break
            ;  case '?': n = get_mdta(buf, buf2, rec->seq, rec->seq_n, pam); break
            ;  case '=': if (ma) n = sw_printaccount(ma, buf, THEBUFSIZE+1);  break
            ;  case 'n':
               case 't':
               case 's':
               case '%':   buf[0] = c == 'n' ? '\n' : c == 't' ? '\t' : c == 's' ? ' ' : '%'
                        ;  n = 1
;  break
            ;  case 'q': { int a = (unsigned char) f[1];
                           int themax = cleanlen;
                           if (pam->length_co > themax) themax = pam->length_co;
                           f++;
                           n = snprintf(buf, THEBUFSIZE+1, "%.*s", (int) rec->clean_n, rec->q);
                                       while (n < themax && n < THEBUFSIZE)
                                          buf[n++] = a;
                                       buf[n] = '\0';
break;
                         }
         ;  }
            escape = 0
      ;  }
         if (n > 0)
         {
#if WE_USE_ZLIB
            if (zippit) gzwrite(fpo, p, n)
         ;  else        fwrite(p, n, 1, fpo)
#else
         fwrite(p, n, 1, fpo)
#endif
      ;  }
         f++
   ;  }
   }
#undef THEBUFSIZE



                           /* low complexity */
int loco
(  struct record* rec
,  struct param* pam
,  double loco_check
)
   {  if (string_complexity(rec->seq, rec->seq_n) < loco_check)
      {  rec->out_message = "loco"
      ;  dump_read(0, pam->fplint, rec, NULL, NULL, pam)
      ;  return 1
   ;  }
      return 0
;  }


/* hierverder segfault if buckid == buck_n or if maxfieldsize == MAXFIELDSIZE.
   should be plenty, no?
*/


void  qc_tally_in
(  struct record* rec
,  struct param* pam
,  int buckid
,  int N
)
   {  int j
   ;  for (j=0;j<N;j++)
      {  switch(rec->seq[j])
         {  case 'A': pam->buck[buckid].qc_cpsin.A[j]++; break
         ;  case 'C': pam->buck[buckid].qc_cpsin.C[j]++; break
         ;  case 'G': pam->buck[buckid].qc_cpsin.G[j]++; break
         ;  case 'T': pam->buck[buckid].qc_cpsin.T[j]++; break
         ;  case 'N': pam->buck[buckid].qc_cpsin.N[j]++; break
      ;  }
         if (rec->q_n)
         {  int kk = MAXQUALITY * j + ((unsigned char) rec->q[j])
;if(0)fprintf(stderr, "testing2 recid=%d buckid=%d qid=%d nt=%d val=%d rootp=%p kk=%d seq=%s\n"
,  (int) rec->rs.N
,  buckid
,  (int) (7 & (unsigned char) rec->q[j])
,  (int) j
,  0
,  (void*) (pam->buck[buckid].qc_bcq)
,  kk
,  rec->seq
)
         ;  pam->buck[buckid].qc_bcq[kk] += 1     /* one block is one read-offset */
                                                  /*  [ nt0 size MAXQUALITY ] [ nt1 size MAXQUALITY] */
;        }
      }
      if (N > pam->buck[buckid].global_in_max_size)
      pam->buck[buckid].global_in_max_size = N
;  }



void  qc_tally_out
(  struct record* rec
,  struct param* pam
,  int buckid
,  int N
)
   {  int j
   ;  pam->buck[buckid].qc_lenout[N]++
   ;  for (j=0;j<N;j++)
      {  switch(rec->seq[j])
         {  case 'A': pam->buck[buckid].qc_cpsout.A[j]++; break
         ;  case 'C': pam->buck[buckid].qc_cpsout.C[j]++; break
         ;  case 'G': pam->buck[buckid].qc_cpsout.G[j]++; break
         ;  case 'T': pam->buck[buckid].qc_cpsout.T[j]++; break
         ;  case 'N': pam->buck[buckid].qc_cpsout.N[j]++; break
      ;  }
      }
      if (N > pam->buck[buckid].global_out_max_size)
      pam->buck[buckid].global_out_max_size = N
;  }


void qc_brief
(  struct param *pam
,  struct record *rec
,  enum program_mode mode
)
   {  argh
      (  "reaper"
      ,  "check %d errors, %d reads truncated, %d clean, %d lint, %d total"
      ,  (int) rec->rs.N_err
      ,  (int) rec->rs.N_truncated
      ,  (int) rec->rs.N_clean
      ,  (int) rec->rs.N_lint
      ,  (int) rec->rs.N
      )
   ;  if (rec->rs.N_alien)
      argh
      (  "reaper"
      ,  "found %d characters outside [ACGTUNX], average %.1f per read"
      ,  rec->rs.N_alien
      ,  rec->rs.N_alien * 1.0 / rec->rs.N
      )
   ;  if (pam->buck_n > 1)
      {  int i, j
      ;  for (i=0;i<pam->buck_n;i++)
         {  dim sum = 0
         ;  const char* bc = pam->buck[i].barcode5p
         ;  if (!bc)
            bc = pam->buck[i].barcode3p
         ;  for (j=0; j<MAXFIELDSIZE; j++)
            sum += pam->buck[i].qc_lenout[j]
         ;  argh("reaper", "barcode [%s] has %lu clean reads", bc, (long unsigned) sum)
      ;  }
      }
   }


   /* this is the quality control report function. mat-happy! */

void  qc_report
(  struct param* pam
,  struct record* rec
,  enum program_mode mode
)
   {  int i, j, buckid
   ;  for (buckid=0;buckid<=pam->buck_n;buckid++)
      {  struct bucket* buck = pam->buck+buckid
      ;  const char* stage2 = "clean"
      ;  const char* thetag = buck->filetag
      
      ;  if (!thetag)
         thetag = buck->barcode3p
      ;  if (!thetag)                                /* hack, but OK  */
         thetag = buck->barcode5p

      ;  if (!thetag)                                /* hack hack hack */
         thetag = "lane"

/* todo: qc_cpsin and qc_cpsout need to keep track of maximum read lengths */

                  /* qc hack again */
      ;  if (buckid == pam->buck_n && (mode == CASE_3P_BC || mode == CASE_5P_BC))
            thetag = "total"
         ,  stage2 = "miss"

      ;  char* fn_pre_nt  = stringle("%s.%s.report.input.nt", pam->basename, thetag)
      ;  char* fn_pre_q   = stringle("%s.%s.report.input.q", pam->basename, thetag)
      ;  char* fn_post_nt = stringle("%s.%s.report.%s.nt", pam->basename, thetag, stage2)
      ;  char* fn_post_q  = stringle("%s.%s.report.%s.q", pam->basename, thetag, stage2)
      ;  char* fn_post_len= stringle("%s.%s.report.%s.len", pam->basename, thetag, stage2)

      ;  FILE* fp_pre_nt, *fp_pre_q
            ,* fp_post_nt, *fp_post_len

      ;  dim *cpsout[5]
         =  {  buck->qc_cpsout.A
            ,  buck->qc_cpsout.C
            ,  buck->qc_cpsout.G
            ,  buck->qc_cpsout.T
            ,  buck->qc_cpsout.N
            }
      ,  *cpsin[5]
         =  {  buck->qc_cpsin.A
            ,  buck->qc_cpsin.C
            ,  buck->qc_cpsin.G
            ,  buck->qc_cpsin.T
            ,  buck->qc_cpsin.N
            }

      ;  if ((fp_pre_nt = myfopen(fn_pre_nt, "w", 0)))
         {  fprintf(fp_pre_nt, "pos\tA\tC\tG\tT\tN\n")
         ;  for (j=0;j<buck->global_in_max_size;j++)
            fprintf
            (  fp_pre_nt
            ,  "%d\t%d\t%d\t%d\t%d\t%d\n"
            ,  (int) (j+1)
            ,  (int) cpsin[0][j]
            ,  (int) cpsin[1][j]
            ,  (int) cpsin[2][j]
            ,  (int) cpsin[3][j]
            ,  (int) cpsin[4][j]
            )
         ;  myfzclose(fp_pre_nt, 0)
      ;  }

         if ((fp_post_nt = myfopen(fn_post_nt, "w", 0)))
         {  fprintf(fp_post_nt, "pos\tA\tC\tG\tT\tN\n")
         ;  for (j=0;j<buck->global_out_max_size;j++)
            fprintf
            (  fp_post_nt
            ,  "%d\t%d\t%d\t%d\t%d\t%d\n"
            ,  (int) (j+1)
            ,  (int) cpsout[0][j]
            ,  (int) cpsout[1][j]
            ,  (int) cpsout[2][j]
            ,  (int) cpsout[3][j]
            ,  (int) cpsout[4][j]
            )
         ;  myfzclose(fp_post_nt, 0)
      ;  }

         if ((fp_post_len = myfopen(fn_post_len, "w", 0)))
         {  fprintf(fp_post_len, "length\tcount\n")
         ;  for (i=0;i<=buck->global_out_max_size;i++)
            fprintf
            (  fp_post_len
            ,  "%d\t%d\n"
            ,  (int) i
            ,  (int) pam->buck[buckid].qc_lenout[i]
            )
         ;  myfzclose(fp_post_len, 0)
      ;  }

         if ((fp_pre_q = myfopen(fn_pre_q, "w", 0)))
         {  fprintf(fp_pre_q, "pos\tq0\tq10\tq25\tq50\tq75\tq90\tq100\n")
         ;  for (i=0;i<pam->buck[buckid].global_in_max_size;i++)
            {  dim total = 0, total2 = 0, j           /* fixme overflow */
            ;  dim q0=0, q10=0, q25=0, q50=0, q75=0, q90=0, q100=0

            ;  for (j=1;j<MAXQUALITY;j++)             /* offset 1 because of q0 check */
               total += buck->qc_bcq[i * MAXQUALITY + j]

            ;  for (j=1;j<MAXQUALITY;j++)
               {  int kk = i * MAXQUALITY + j
               ;  total2 += buck->qc_bcq[kk]
               ;  if (!q0  && total * 0.0  <  total2)
                  q0 =  j
               ;  if (!q10 && total * 0.10 <= total2)
                  q10 = j
               ;  if (!q25 && total * 0.25 <= total2)
                  q25 = j
               ;  if (!q50 && total * 0.50 <= total2)
                  q50 = j
               ;  if (!q75 && total * 0.75 <= total2)
                  q75 = j
               ;  if (!q90 && total * 0.90 <= total2)
                  q90 = j
               ;  if (!q100 && total * 0.9999 <= total2)
                  q100 = j
            ;  }
               fprintf
               (  fp_pre_q, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
               ,  (int) (i+1)
               ,  (int) q0
               ,  (int) q10
               ,  (int) q25
               ,  (int) q50
               ,  (int) q75
               ,  (int) q90
               ,  (int) q100
               )
         ;  }
            myfzclose(fp_pre_q, 0)
      ;  }
         free(fn_pre_nt)
      ;  free(fn_pre_q)
      ;  free(fn_post_nt)
      ;  free(fn_post_q)
      ;  free(fn_post_len)
   ;  }
      {  char* fn_post_tally= stringle("%s.sumstat", pam->basename)
      ;  FILE* fp_post_tally
      ;  if ((fp_post_tally = myfopen(fn_post_tally, "w", 0)))
         {  fprintf(fp_post_tally, "discarded_no_barcode_match=%u\n", rec->rs.N_nomatch)
         ;  fprintf(fp_post_tally, "discarded_tabu_contamination=%u\n", rec->rs.N_tabu)
         ;  fprintf(fp_post_tally, "discarded_no_5p_sinsert=%u\n", rec->rs.N_5p_nosinsert)
         ;  fprintf(fp_post_tally, "discarded_no_anchor=%u\n", rec->rs.N_anchor)
         ;  fprintf(fp_post_tally, "discarded_NNN=%u\n", rec->rs.N_nn)
         ;  fprintf(fp_post_tally, "discarded_BBB=%u\n", rec->rs.N_bb)
         ;  fprintf(fp_post_tally, "discarded_AAA=%u\n", rec->rs.N_aa)
         ;  fprintf(fp_post_tally, "discarded_QQQ=%u\n", rec->rs.N_qq)
         ;  fprintf(fp_post_tally, "discarded_tri=%u\n", rec->rs.N_tri)
         ;  fprintf(fp_post_tally, "discarded_length_cutoff=%u\n", rec->rs.N_lenlo)
         ;  fprintf(fp_post_tally, "discarded_barcode_conflict=%u\n", rec->rs.N_cfl)
         ;  fprintf(fp_post_tally, "discarded_low_complexity=%u\n", rec->rs.N_loco)
         ;  fprintf(fp_post_tally, "no_message=%u\n", rec->rs.N_nomsg)
         ;  fprintf(fp_post_tally, "total_discarded=%u\n", rec->rs.N_lint)
         ;  fprintf(fp_post_tally, "total_error=%u\n", rec->rs.N_err)
         ;  fprintf(fp_post_tally, "total_accepted=%u\n", rec->rs.N_clean)
         ;  fprintf(fp_post_tally, "total_input=%u\n", rec->rs.N)
         ;  fprintf(fp_post_tally, "conflicts_resolved=%u\n", rec->rs.N_cflr)
         ;  fprintf(fp_post_tally, "bases_removed_by_quality=%u\n", rec->rs.B_QQ)
         ;  fprintf(fp_post_tally, "bases_removed_by_qualityb=%u\n", rec->rs.B_BB)
         ;  fprintf(fp_post_tally, "bases_removed_by_tritail=%u\n", rec->rs.B_tri)
         ;  fprintf(fp_post_tally, "bases_removed_by_nnn=%u\n", rec->rs.B_NN)
         ;  fprintf(fp_post_tally, "bases_removed_by_polya=%u\n", rec->rs.B_AA)
         ;  myfzclose(fp_post_tally, 0)
         ;  free(fn_post_tally)
      ;  }
      }
   }


int buckets_prepare_qc
(  struct bucket* buck
,  int buck_n
,  enum program_mode mode
)
   {  int i
   ;  for (i=0;i<=buck_n;i++)
      {  buck[i].cache_lint_offset = 0
      ;  buck[i].cache_score = 0
      ;  buck[i].cache_accept = 0
      ;  memset(buck[i].qc_cpsin.A, 0, sizeof buck[i].qc_cpsin.A)
      ;  memset(buck[i].qc_cpsin.C, 0, sizeof buck[i].qc_cpsin.C)
      ;  memset(buck[i].qc_cpsin.G, 0, sizeof buck[i].qc_cpsin.G)
      ;  memset(buck[i].qc_cpsin.T, 0, sizeof buck[i].qc_cpsin.T)
      ;  memset(buck[i].qc_cpsin.N, 0, sizeof buck[i].qc_cpsin.N)
      ;  memset(buck[i].qc_cpsout.A, 0, sizeof buck[i].qc_cpsout.A)
      ;  memset(buck[i].qc_cpsout.C, 0, sizeof buck[i].qc_cpsout.C)
      ;  memset(buck[i].qc_cpsout.G, 0, sizeof buck[i].qc_cpsout.G)
      ;  memset(buck[i].qc_cpsout.T, 0, sizeof buck[i].qc_cpsout.T)
      ;  memset(buck[i].qc_cpsout.N, 0, sizeof buck[i].qc_cpsout.N)
      ;  memset(buck[i].qc_lenout, 0, sizeof buck[i].qc_lenout)

                           /* fixme memrelease */
      ;  // if (!(buck[i].qc_bcq = calloc(1 + (MAXQUALITY * MAXFIELDSIZE), sizeof(dim))))
      ;  if (!(buck[i].qc_bcq = calloc((MAXQUALITY * MAXFIELDSIZE) , sizeof(dim))))
         X_FINISHED("cannot allocate memory for quality tracking")
      ;  buck[i].global_in_max_size = 0
      ;  buck[i].global_out_max_size = 0
      ;  buck[i].barcode_start = 0
      ;  buck[i].barcode_end = 0
   ;  }
                  /* fixme very ugly; init code only for last "special" bucket */
                  /* rip and refactor at some point */
      buck[buck_n].barcode3p = NULL
   ;  buck[buck_n].adaptor3p = "the-3-prime-adaptor"
   ;  buck[buck_n].sinsert3p = NULL

   ;  buck[buck_n].barcode5p = NULL
   ;  buck[buck_n].sinsert5p = NULL
   ;  buck[buck_n].cat3p = NULL
   ;  buck[buck_n].filetag = NULL
   ;  return 0
;  }



void data_release
(  struct param* pam
)
   {  int i
   ;  if (pam->buck)
      {  for (i=0;i<=pam->buck_n;i++)
         {  struct bucket* buck = pam->buck+i
         ;  struct atabu* boo = &buck->tabu2, *boop = NULL     /* boop'revious */
         ;  int n = 0
         ;  while (boo)                   /* hierverder: corrupt if meta file cannot be opened! */
            {  if (boo->tabuseq)
               free(boo->tabuseq)
            ;  boop = boo->next
            ;  if (n++)
               free(boo)
            ;  boo = boop
         ;  }
            if (buck->qc_bcq) free(buck->qc_bcq)
         ;  if (buck->fname) free(buck->fname)
      ;  }
         free(pam->buck)
   ;  }

      pam->buck = NULL
   ;  pam->buck_n = -1
   ;  table_release(&pam->metadata)
;  }


int barcodes_io_open
(  struct bucket* buck
,  const char* basename
,  int buck_n
,  int zippit
)
   {  int i
   ;  for (i=0;i<buck_n;i++)
      {  const char* thetag = buck[i].filetag
      ;  if (!thetag)
         thetag = buck[i].barcode3p
      ;  if (!thetag)
         thetag = buck[i].barcode5p
      ;  buck[i].fname = stringle("%s.%s.%s", basename, thetag, zippit ? "clean.gz" : "clean")
      ;  if (!(buck[i].fp = myfopen(buck[i].fname, ZWMODE, zippit)))
         {  argh("reaper", "cannot open file %s", buck[i].fname)
         ;  return 1
      ;  }
      }
      return 0
;  }


int barcodes_io_close
(  struct bucket* buck
,  int buck_n
,  int zippit
)
   {  int i
   ;  for (i=0;i<buck_n;i++)
      myfzclose(buck[i].fp, zippit)
   ;  return 0
;  }

               /* fixme, could have seqp and qp pointers */
void rec_remove_prefix
(  struct record* rec
,  int rm
)
   {  if (rm > rec->seq_n)
      rm = rec->seq_n
   ;  memmove(rec->seq, rec->seq+rm, rec->seq_n - rm)
   ;  rec->seq_n -= rm
   ;  rec->seq[rec->seq_n] = '\0'
   ;  if (rec->q_n)
         memmove(rec->q, rec->q+rm, rec->q_n - rm)
      ,  rec->q_n   -= rm
      ,  rec->q[rec->q_n] = '\0'
;  }


int column_require
(  int i
,  const char* str
)
   {  if (i < 0)
      X_FINISHED("need column %s", str)
   ;  return 0
;  }


      /* fixme: seem to have multiple places with barcode instantiation.
       * deduplicate.
      */
int adaptors_set_case_3p
(  struct bucket* buck
,  int* buck_n
,  struct param* pam
)
   {  buck[0].barcode5p = NULL

   ;  buck[0].barcode3p = NULL
   ;  buck[0].sinsert3p = NULL
   ;  buck[0].cat3p     = NULL
   ;  buck[0].sinsert5p = pam->sinsert_5p
   ;  buck[0].adaptor3p = pam->adaptor_3p
   ;  buck[0].tabu2     = pam->the_tabu2

   ;  if (g_mode == CASE_NOBC && !pam->adaptor_3p)
      X_FINISHED("-geom 3p need a 3p adaptor")

   ;  if
      (  g_restrict > 0
      && strlen(pam->adaptor_3p) > g_restrict
/* fixme: old code had
      && strlen(pam->adaptor_3p) < g_restrict
   and led to different results sometimes.
*/
      )
      ((char*) buck[0].adaptor3p)[g_restrict] = '\0'

;fprintf(stderr, "adapter %s\n", buck[0].adaptor3p)
   ;  buck_n[0] = 1
   ;  return 0
;  }


                                                /* ugly, designtabu */
void parse_tabu
(  struct atabu* boo
,  const char*  s
,  int extent
)
   {  struct match_requirement mr_tabu      = { 16, 2, 1, 0 }
   ;  const char* tok1 = NULL
   ;  do
      {  char* s2 = NULL
      ;  tok1 = strchr(s, ',')
      ;  s2 = stringle("%.*s", (int) (tok1 ? (tok1 - s) : strlen(s)), s)
      ;  boo->tabuseq = s2
      ;  boo->tabulen = strlen(s2)
      ;  boo->mr_tabu = mr_tabu
      ;  boo->next = NULL

      ;  if (extent && boo->tabulen > extent)
            boo->tabulen = extent
         ,  boo->tabuseq[extent] = '\0'

      ;  if (tok1)
         {  boo->next = myalloc(sizeof boo[0])
         ;  boo = boo->next
         ;  s = tok1 + 1
      ;  }
      }
      while(tok1)
   ;  return
;  }


int buckets_alloc
(  struct param* pam
,  int n
)
   {  if (n > 64)
      argh("reaper", "that is a lot of barcodes, that is. %d of'em", n)

   ;  if (!(pam->buck = calloc((n+1), sizeof pam->buck[0])))
      X_FINISHED("unable to allocate data for %d barcodes", n)
   ;  pam->buck_n = n
   ;  return 0
;  }


int config_read
(  const char* fn
,  struct param* pam
,  enum program_mode mode
,  int extent
)
   {  int id_barcode    ,  doc_barcode =  0
   ;  int id_3p_ad      ,  doc_3p_ad   =  1
   ;  int id_3p_si      ,  doc_3p_si   =  2

   ;  int id_5p_si      ,  doc_5p_si   =  3
   ;  int id_tabu       ,  doc_tabu    =  4
   ;  int id_filetag    ,  doc_filetag =  5

   ;  int n_columns_found
   ;  int n_barcodes = 0
   ;  int barcodes_length = 0
   ;  struct table* thetable = &pam->metadata
   ;  ZFILE ip
   ;  struct bucket* buck = NULL

   ;  if (!(ip = myfopen(fn, ZRMODE, 1)))
      return 1

   ;  const char* dochandle[6] = { "barcode", "3p-ad", "3p-si", "5p-si", "tabu", "filetag" }

   ;  X_CASCADE_CLEAN(table_read(ip, thetable, 's', TABLE_WITH_COLNAMES), myfzclose(ip, 1))
   ;  X_CASCADE(buckets_alloc(pam, thetable->n_rows))
   ;  buck = pam->buck

   ;  id_barcode  =  table_column_index(thetable, dochandle[doc_barcode])
   ;  id_3p_ad    =  table_column_index(thetable, dochandle[doc_3p_ad])
   ;  id_3p_si    =  table_column_index(thetable, dochandle[doc_3p_si])
   ;  id_5p_si    =  table_column_index(thetable, dochandle[doc_5p_si])
   ;  id_tabu     =  table_column_index(thetable, dochandle[doc_tabu])
   ;  id_filetag  =  table_column_index(thetable, dochandle[doc_filetag])

   ;  n_columns_found
      =     (  id_3p_ad >= 0     )
         +  (  id_barcode >= 0   )
         +  (  id_3p_si >= 0     )
         +  (  id_5p_si >= 0     )
         +  (  id_tabu >= 0      )
         +  (  id_filetag >= 0   )

   ;  switch(mode)
      {  case CASE_NOBC
         :  {  X_CASCADE(column_require(id_3p_ad, dochandle[doc_3p_ad]))
            ;  X_CASCADE(column_require(id_tabu, dochandle[doc_tabu]))
            ;  if (n_columns_found > 2)
               argh("reaper", "only using columns %s and %s", dochandle[doc_3p_ad], dochandle[doc_tabu])
            ;  break
         ;  }

         case CASE_3P_BC
         :  {  X_CASCADE(column_require(id_3p_si, dochandle[doc_3p_si]))
            ;  X_CASCADE(column_require(id_barcode, dochandle[doc_barcode]))
            ;  X_CASCADE(column_require(id_tabu, dochandle[doc_tabu]))
            ;  X_CASCADE(column_require(id_3p_ad, dochandle[doc_3p_ad]))
            ;  if (n_columns_found > 4)
               argh("reaper", "only using columns %s, %s and %s", dochandle[doc_barcode], dochandle[doc_3p_si], dochandle[doc_tabu])
            ;  break
         ;  }

         case CASE_5P_BC
         :  {  X_CASCADE(column_require(id_barcode, dochandle[doc_barcode]))
            ;  X_CASCADE(column_require(id_5p_si, dochandle[doc_5p_si]))
            ;  X_CASCADE(column_require(id_3p_ad, dochandle[doc_3p_ad]))
            ;  X_CASCADE(column_require(id_tabu, dochandle[doc_tabu]))
            ;  if (n_columns_found > 4)
               argh
               (  "reaper"
               ,  "only using columns %s, %s, %s and %s"
               ,  dochandle[doc_barcode]
               ,  dochandle[doc_5p_si]
               ,  dochandle[id_3p_ad], dochandle[doc_tabu]
               )
            ;  break
         ;  }

         default
         :  ;
      }

   ;  for (n_barcodes=0; n_barcodes < thetable->n_rows; n_barcodes++)
      {  int n = n_barcodes
      ;  buck[n].barcode5p =  NULL
      ;  buck[n].sinsert5p =  NULL
      ;  buck[n].sinsert5p_n = 0

      ;  buck[n].barcode3p =  NULL
      ;  buck[n].sinsert3p =  NULL
      ;  buck[n].adaptor3p =  NULL
      ;  buck[n].sinsert3p_n = 0
      ;  buck[n].cat3p     =  NULL
      ;  buck[n].filetag   =  NULL

      ;  buck[n].barcode_n = 0

#define hte(s) (strcmp(s, "-") ? s : "")           /* hyphen-to-empty */

      ;  if (id_filetag >= 0)
         buck[n].filetag = hte(table_elem_string(thetable, n, id_filetag))

      ;  if (mode == CASE_NOBC)
         {  buck[n].adaptor3p =  hte(table_elem_string(thetable, n, id_3p_ad))
         ;  parse_tabu(&buck[n].tabu2, hte(table_elem_string(thetable, n, id_tabu)), extent)
         ;  if (id_5p_si >= 0)
               buck[n].sinsert5p = hte(table_elem_string(thetable, n, id_5p_si))
            ,  buck[n].sinsert5p_n = strlen(buck[n].sinsert5p)                   /* fixme do length assignments in unified place */
      ;  }
         else if (mode == CASE_3P_BC)
         {  buck[n].barcode3p = table_elem_string(thetable, n, id_barcode)
         ;  buck[n].adaptor3p = hte(table_elem_string(thetable, n, id_3p_ad))
         ;  buck[n].sinsert3p = hte(table_elem_string(thetable, n, id_3p_si))
         ;  parse_tabu(&buck[n].tabu2, hte(table_elem_string(thetable, n, id_tabu)), extent)

         ;  buck[n].sinsert3p_n = strlen(buck[n].sinsert3p)
         ;  buck[n].barcode_n = strlen(buck[n].barcode3p)
         ;  buck[n].cat3p     =  stringle("%s%s%s", buck[n].sinsert3p, buck[n].barcode3p, buck[n].adaptor3p)

         ;  if (buck[n].sinsert3p_n > 1)
            X_FINISHED("3p sequence insert of size > 1 not yet supported; ask Stijn")
         ;  if (n == 0)
            barcodes_length = strlen(buck[n].barcode3p)
         ;  else if (barcodes_length != strlen(buck[n].barcode3p))
            argh("reaper", "barcodes with different lengths")
            /* X_FINISHED("barcodes with different lengths") */
      ;  }
         else if (mode == CASE_5P_BC)
         {  buck[n].barcode5p = table_elem_string(thetable, n, id_barcode)
         ;  buck[n].sinsert5p = hte(table_elem_string(thetable, n, id_5p_si))
         ;  buck[n].adaptor3p = hte(table_elem_string(thetable, n, id_3p_ad))
         ;  parse_tabu(&buck[n].tabu2, hte(table_elem_string(thetable, n, id_tabu)), extent)

         ;  buck[n].sinsert5p_n = strlen(buck[n].sinsert5p)
         ;  buck[n].barcode_n = strlen(buck[n].barcode5p)
         ;  if (n == 0)
            barcodes_length = strlen(buck[n].barcode5p)
         ;  else if (barcodes_length != strlen(buck[n].barcode5p))
            argh("reaper", "barcodes with different lengths")
            /* X_FINISHED("barcodes with different lengths") */
      ;  }

         if (extent > 0)
         {  if (buck[n].adaptor3p && strlen(buck[n].adaptor3p) > extent)
            ((char*) buck[n].adaptor3p)[extent] = '\0'
      ;  }

         argh
         (  "reaper"
         ,  "scanned line %d sinsert3p [%s] barcode3p [%s] adaptor3p [%s] cat3p [%s]"
            " barcode5p [%s] sinsert5p [%s] tabu [%s]"
         ,  (int) (n+1)
         ,  buck[n].sinsert3p ? buck[n].sinsert3p : "NA"
         ,  buck[n].barcode3p ? buck[n].barcode3p : "NA"
         ,  buck[n].adaptor3p ? buck[n].adaptor3p : "NA"
         ,  buck[n].cat3p     ? buck[n].cat3p     : "NA"
         ,  buck[n].barcode5p ? buck[n].barcode5p : "NA"
         ,  buck[n].sinsert5p ? buck[n].sinsert5p : "NA"
         ,  buck[n].tabu2.tabuseq ? buck[n].tabu2.tabuseq : ""
         )
      ;  {  struct atabu* boo = buck[n].tabu2.next
         ;  while (boo)
            {  fprintf(stderr, " [%s]", boo->tabuseq)
            ;  boo = boo->next
         ;  }
         }
      }
   ;  return 0
;  }


            /* check: below comment up to date? */
            /* The second condition is because we only fill part of the alignment matrix.
             * in the equality case the matrix will be emtpy, and sw_trace_trim will crash.
             * The easiest way is to check here (BUT, fixme nevertheless; sw_trace_trim should
             * do better; it's probably crashing because max_ij is never increased and starts at 0).
             * Note in barcode mode we need at least the barcode (plus a little bit more,
             * likely), so an initial read of length cleanlength (before aligning)
             * is never going to make it.
            */

int check_for_tabu
(  struct record* rec
,  unsigned len
,  struct param* pam
,  struct atabu* boo
)
   {  struct sw_alninfo ai = { 0 }
   ;  char* s = rec->seq
   ;  char* q = rec->q
   ;  int indel_allowed = pam->mr_tabu.mr_maxgap > 0
   ;  char c_s = s[len]                      /* hack/kludge; note that len is user specified */
   ;  char c_q = q[len]                      /* hack/kludge NOTE harmless in case q_n == 0 */
   ;  int found_alignment = 0
   ;  SWNUM data[MATRIX_SIZE] = { 0 }

   ;  if
      (  !boo
      || !boo->tabulen
      || !pam->mr_tabu.mr_minlen
      )
      return 0

   ;  s[len] = '\0'
   ;  q[len] = '\0'

   ;  do
      {  if (sw_fill(&ai, data, MATRIX_SIZE, s, boo->tabuseq, indel_allowed, &pam->swp_tabu))
         X_IGNORE("sw error in check_for_tabu")
      ;  rec->rs.N_cells += ai.n_cells_used
      ;  sw_trace_trim(&ai, &pam->swp_tabu, ai.max_ij, SW_TRIM)       /* check_for_tabu */

      ;  if (sw_match_any(&ai, &pam->mr_tabu))
         {  found_alignment = 1
         ;  break
      ;  }
         boo = boo->next
   ;  }
      while (boo && boo->tabulen)

         /* in check_for_tabu, fplint; this one printed before read */
   ;  if
      (  pam->verbose & VB_ALIGN
      && rec->rs.N % pam->sample_lint == 0
      )
      sw_pp2("read", x_tabu[found_alignment], -1, &ai, pam->fplint, pam->zippit, rec->rs.N, 0, NULL)
         /* in check_for_tabu */

   ;  s[len] = c_s                           /* hack/kludge */
   ;  q[len] = c_q                           /* hack/kludge */
   ;  return found_alignment
;  }


         /* note: mr_minlen not used! */
         /* note: sinsert is left in the alignment */
int sw_match_5p_required
(  struct sw_alninfo* ai
,  int adaptor_length
,  struct match_requirement* mr
)
   {  int n_edit =   ai->n_insl + ai->n_insr + ai->n_subst
   ;  int n_border_left = THEMAX(ai->lft_start, ai->rgt_start) - 1
   ;  int n_border_right = adaptor_length - ai->lft_end
   ;  int n_border = n_border_left + n_border_right

;if(0)fprintf(stderr, "border-left: %d border-right: %d\n", n_border_left, n_border_right)

   ;  return
         (  ai->n_match + n_edit + n_border >= adaptor_length     /* cover entire adapter */
         && n_edit + n_border <= mr->mr_maxedit                   /* edit requriement */
         && ai->n_insl + ai->n_insr <= mr->mr_maxgap              /* gap requirement  */
         )
      ?  1
      :  0
;  }


                           /* only called by ntnt */
int check_5p_sinsert
(  struct record* rec
,  struct param* pam
,  ZFILE fpclean          /* clean file is barcode specific */
,  const char** msgp
)
   {  int keepall = pam->modes & MODE_KEEPALL
   ;  if
      (  !pam->sinsert_5p
      || !strlen(pam->sinsert_5p)
      )
      return 0

;if(0)fprintf(stderr, "insert %s\n", pam->sinsert_5p)
   ;  pam->swp_5p.right_limit =  strlen(pam->sinsert_5p) + pam->mr_5p_require.mr_offset

   ;  {  struct sw_alninfo ai = { 0 }
      ;  SWNUM data[MATRIX_SIZE] = { 0 }
      ;  int adaptor_length = strlen(pam->sinsert_5p)
      ;  int indel_allowed = pam->mr_5p_require.mr_maxgap > 0
      ;  int found_alignment = 0

      ;  if (sw_fill(&ai, data, MATRIX_SIZE, pam->sinsert_5p, rec->seq, indel_allowed, &pam->swp_5p))
         X_FINISHED_VAL(-1, "sw error in check_5p_sinsert")
      ;  rec->rs.N_cells += ai.ni * ai.nj
      ;  sw_trace_trim
         (  &ai
         ,  &pam->swp_5p
         ,  ai.max_ij
         ,  indel_allowed ? SW_TRIM : 0
         )         /* check_5p_sinsert */

      ;  found_alignment = sw_match_5p_required(&ai, adaptor_length, &pam->mr_5p_require)

;if(0)
 fprintf(stderr, "%s %s %d\n", pam->sinsert_5p, rec->seq, found_alignment)
,sw_pp2("read", x_5psi[found_alignment], -1, &ai, stdout, 0, 0, rec->rs.N, NULL)

            /* in check_5p_sinsert, fpclean, fplint; caller ntnt returns immediately */
      ;  if
         (  pam->verbose & VB_ALIGN
         && rec->rs.N % pam->sample_lint == 0
         )
         sw_pp2("read", x_5psi[found_alignment], -1, &ai, pam->fplint, pam->zippit, rec->rs.N, 0, NULL)
            /* in check_5p_sinsert */

      ;  if (found_alignment)
         return ai.rgt_end + (ai.ni-1-ai.lft_end)        /* ai.ni-1 is the length of left sequence */
      ;  else
         {  if (keepall)
            {  msgp[0] = "5p-sinsert"
            ;  return rec->seq_n
         ;  }
            rec->rs.N_5p_nosinsert++
         ;  rec->rs.N_lint++
         ;  rec->out_message = "5p-absent"
         ;  dump_read(0, pam->fplint, rec, NULL, NULL, pam)
         ;  return -1
      ;  }
      }
      return 0
;  }


/* fixme; what if the input read has length smaller than cleanlength and no
 * adaptor or tabu is found? detected?
*/

int ntnt                  /* nucleotide nip and tuck */
(  struct record* rec
,  struct param* pam       /* lint file is taken from pam, non-specific */
,  ZFILE fpclean          /* clean file is barcode specific */
)
   {  int clean = 0, lint_offset = 1, p_5p_offset = 0
   ;  int indel_allowed = pam->mr_3p_global.mr_maxgap > 0
   ;  int keepall = pam->modes & MODE_KEEPALL
   ;  const char* msg = "clean"
   ;  int found_alignment = 0

   ;  struct sw_alninfo ai = { 0 }
   ;  struct match_account ma = { 'n', 0, 0, -1, -1 }
   ;  struct sw_param myswp = pam->swp_3p

   ;  SWNUM data[MATRIX_SIZE] = { 0 }

                           /* returns the offset of the first good bit */
                           /* updates rec->rs                          */
   ;  p_5p_offset = check_5p_sinsert(rec, pam, fpclean, &msg)

   ;  if (p_5p_offset < 0)
      return 1
   ;  else if (p_5p_offset > 0)
      rec_remove_prefix(rec, p_5p_offset)

   ;  if (g_guard > 0)
      myswp.left_skip = g_guard

                           /* hierverder: if empty read, everything below is still acted
                            * upon; do we hit any of those conditions?
                           */

                           /* CASE_NOBC may or may not set pam->adaptor_3p and pam->tabu2
                            * fixme check.
                           */
   ;  if (pam->adaptor_3p)
      {  if
         (  sw_fill
            (  &ai
            ,  data
            ,  MATRIX_SIZE
            ,  rec->seq
            ,  pam->adaptor_3p
            ,  indel_allowed
            ,  &myswp
         )  )
         X_FINISHED("sw error in nip_and_tuck")
      ;  rec->rs.N_cells += ai.ni * ai.nj
      ;  sw_trace_trim(&ai, &myswp, ai.max_ij, SW_TRIMLEFT)      /* decideme: why not SW_TRIMRIGHT? */
   ;  }

      ma.nnn_offset = get_nnn_offset(ai.left, ai.ni-1, pam->nnn_count, pam->nnn_outof)
   ;  ma.qqq_offset = get_qqq_offset(rec, pam->q_cutoff, pam->q_winlen, pam->q_winofs, pam->q_readofs)
   ;  ma.bbb_offset = pam->modes & MODE_QC_LATE ? get_bbb_offset(rec) : -1
   ;  ma.aaa_offset = get_polya_offset(ai.left, ai.ni-1, pam->polya)
   ;  ma.dust_offset= get_dust_offset(ai.left, ai.ni-1, pam->taildust, pam->taildust_bases)

         /* NOTE everything in alignment space is offset one relative to string
          * space (C's normal 0-based indexing).  In alignment space the
          * alignment starts at ai.lft_start, so the good bit ends at
          * ai.lft_start-1. In string space that's at ai.lft_start-2. The
          * length of the good bit is thus ai.lft_start - 1: that's what we
          * need to compare with pam->cleanlength.

          * NNN space is just string-based, so lint_offset is not adjusted when
          * set to ma.nnn_offset. The same holds for BBB space.

          * Finally, if 0 == pam->cleanlength, then we should always accept the
          * read, even though the length cutoff comparison might *still* fail
          * due to the ai.rgt_start offsetting -- we may have a clean
          * readlength of -2 if for instance 3 == rgt_start.
         */

   ;  if (pam->adaptor_3p)
      found_alignment = sw_match_any(&ai, &pam->mr_3p_global)

   ;  if
      (  pam->adaptor_3p[0]               /* fixme/docme */
      && pam->verbose & VB_ALIGN
      && rec->rs.N % pam->sample_lint == 0
      )
      sw_pp2("read", x_3pgb[found_alignment], -1, &ai, pam->fplint, pam->zippit, rec->rs.N, 0, NULL)

   ;  if (found_alignment)
      {  ma.match_type = 'g'
      ;  clean = (ai.lft_start + 1) > pam->cleanlength + ai.rgt_start
      ;  msg = "len-global"
      ;  if (clean)
         lint_offset = (ai.lft_start + 1) - ai.rgt_start
      ;  else if (!pam->cleanlength && ai.lft_start < ai.rgt_start)
            clean = 1
         ,  lint_offset = 1
      ;  else
         rec->rs.N_lenlo++
      ;
      }

                              /* below, first use the overal best alignment.
                               * This may find a prefix match followed by low complexity sequence.
                               * Then look at the tail of the read only.
                              */
      else if (pam->adaptor_3p && pam->mr_3p_prefix.mr_minlen)
      {  unsigned z1 = sw_match_end(&ai, &ma, &myswp, pam->p_3p_tth, &pam->mr_3p_prefix)
;if(0)fprintf(stderr, "z1 %d\n", (int) z1)
      ;  found_alignment = z1 && ma.suffix_complexity < pam->suffix_complexity_max
      ;  if (!found_alignment)
         {  sw_trace_trim(&ai, &myswp, ai.max_ej, SW_TRIMLEFT)     /* DON'T trim right; it needs to stay anchored */
         ;  found_alignment = sw_match_end(&ai, &ma, &myswp, pam->p_3p_tth, &pam->mr_3p_prefix)
         ;  if
            (  pam->adaptor_3p[0]               /* fixme/docme */
            && pam->verbose & VB_ALIGN
            && rec->rs.N % pam->sample_lint == 0
            )
            sw_pp2("read", x_3ppx[found_alignment], -1, &ai, pam->fplint, pam->zippit, rec->rs.N, 0, NULL)
;if(0)fprintf(stderr, "found_alignment %d\n", (int) found_alignment)
      ;  }

         if (found_alignment)
         {  if (ma.n_perfect)
               clean = ai.ni > pam->cleanlength + ma.n_perfect
            ,  ma.match_type = 't'
         ;  else
               clean = (ai.lft_start + 1) > pam->cleanlength + ai.rgt_start
            ,  ma.match_type = 'p'

         ;  msg = "len-suffix"
         ;  if (clean)
            {  lint_offset
               =  ma.n_perfect
               ?  ai.ni - ma.n_perfect
               : (ai.lft_start + 1) - ai.rgt_start
         ;  }
            else if (!pam->cleanlength && ai.lft_start < ai.rgt_start)
            {  clean = 1
            ;  lint_offset = 1
         ;  }
            else
            rec->rs.N_lenlo++
      ;  }
         else if (ai.ni > pam->cleanlength)
            clean = 1
         ,  lint_offset = ai.ni
   ;  }

      else
         clean = 1
      ,  lint_offset = rec->seq_n + 1

                        /* This may happen due to the correction for right start
                         * (the alignment can start at base 2 in the adapter,
                         * base 0 in the read).
                         * Set to 1 by way of defensive programming.
                        */
   ;  if (!lint_offset)
      lint_offset = 1

                        /* fixme: funcify */
   ;  if (clean && ma.nnn_offset > 0)
      {  rec->rs.B_NN += rec->seq_n+1 - ma.nnn_offset
      ;  if (ma.nnn_offset <= pam->cleanlength)
         {  clean = 0
         ;  rec->rs.N_nn++
         ;  msg = "NNN"
      ;  }
         else if (ma.nnn_offset < lint_offset)
         lint_offset = ma.nnn_offset
   ;  }

                        /* fixme: funcify */
      if (clean && ma.qqq_offset > 0)
      {  rec->rs.B_QQ += rec->seq_n+1 - ma.qqq_offset
      ;  if (ma.qqq_offset <= pam->cleanlength)
         {  clean = 0
         ;  rec->rs.N_qq++
         ;  msg = "lowmq"
      ;  }
         else if (ma.qqq_offset < lint_offset)
         lint_offset = ma.qqq_offset
   ;  }

                                       /* ach, require 3 bases, why not */
      if (pam->taildustlate && clean && lint_offset > 3)
      {  int dust_offset = get_dust_offset(ai.left, lint_offset-1, pam->taildustlate, pam->taildustlate_bases)
      ;  if (dust_offset >= 0)
         {  if (ma.dust_offset < 0 || dust_offset < ma.dust_offset)
            ma.dust_offset = dust_offset
      ;  }
      }

                        /* fixme: funcify */
      if (clean && ma.dust_offset > 0)
      {  rec->rs.B_tri += rec->seq_n+1 - ma.dust_offset
      ;  if (ma.dust_offset <= pam->cleanlength)
         {  clean = 0
         ;  rec->rs.N_tri++
         ;  msg = "tritail"
      ;  }
         else if (ma.dust_offset < lint_offset)
         lint_offset = ma.dust_offset
   ;  }


                        /* fixme: funcify */
      if (clean && ma.bbb_offset > 0)
      {  rec->rs.B_BB += rec->seq_n+1 - ma.bbb_offset
      ;  if (ma.bbb_offset <= pam->cleanlength)
         {  clean = 0
         ;  rec->rs.N_bb++
         ;  msg = "lowq"
      ;  }
         else if (ma.bbb_offset < lint_offset)
         lint_offset = ma.bbb_offset
   ;  }


                        /* fixme: funcify
                         * hierverder: is B_AA incremenent OK? (doublecounting)
                        */
      if (clean && ma.aaa_offset > 0)
      {  rec->rs.B_AA += rec->seq_n+1 - ma.aaa_offset
      ;  if (ma.aaa_offset <= pam->cleanlength)
         {  clean = 0
         ;  rec->rs.N_aa++
         ;  msg = "polya"
      ;  }
         else if (ma.aaa_offset < lint_offset)
         lint_offset = ma.aaa_offset
;if(0)fprintf(stderr, "polya %d lint %d clean %d seq %s %d\n", (int) ma.aaa_offset, (int) lint_offset, (int) clean, rec->seq, (int) rec->rs.N)
   ;  }


                        /* fixme: funcify */
      if (clean && pam->length_co > 0)
      {  rec->rs.B_lco += rec->seq_n - pam->length_co
      ;  if (pam->length_co+1 < lint_offset)
         lint_offset = pam->length_co+1
   ;  }


      if (clean && pam->trintco && pam->trintco < trintscore(rec->seq, lint_offset-1))
      {  clean = 0
      ;  rec->rs.N_tri++
      ;  msg = "trint"
   ;  }

                       /* caller sets pam->the_tabu2 */
      if
      (  clean
      && pam->the_tabu2.tabulen
      && check_for_tabu(rec, lint_offset-1, pam, &pam->the_tabu2)
      )
      {  if (keepall)
         lint_offset = 1
      ;  else
         {  clean = 0
         ;  rec->rs.N_tabu++
      ;  }
         msg = "tabu"
   ;  }

      clean |= keepall

   ;  if (clean)
      {  rec->rs.N_clean++
      ;  rec->clean_n = lint_offset - 1
      ;  rec->out_message = msg
      ;  dump_read(1, fpclean, rec, &ai, &ma, pam)
   ;  }
      else
      {  rec->rs.N_lint++
      ;  if (!msg)
            msg = "no message (weird)"
         ,  rec->rs.N_nomsg++
      ;  rec->out_message = msg
      ;  dump_read(0, pam->fplint, rec, &ai, &ma, pam)
      ;  rec->clean_n = -1
   ;  }

      return 0
;  }


                              /* fixfixme: adaptor settings */
int do_adaptor_3p
(  struct record* rec
,  struct param* pam
)
   {  pam->adaptor_3p = pam->buck[0].adaptor3p
   ;  pam->the_tabu2  = pam->buck[0].tabu2
   ;  pam->sinsert_5p = pam->buck[0].sinsert5p
#if 0
         =    (pam->buck[bestid].sinsert5p_n)
            ?  pam->buck[bestid].sinsert5p
            :  NULL
#endif
   ;  X_CASCADE(ntnt(rec, pam, pam->fpclean))

   ;  if (rec->clean_n >= 0)
      qc_tally_out(rec, pam, pam->buck_n, rec->clean_n)     /* /specialbucket/ */
   ;  return 0
;  }


int do_barcode_3p
(  struct record *rec
,  struct param *pam
)
   {  int NNN_offset = get_nnn_offset(rec->seq, rec->seq_n, pam->nnn_count, pam->nnn_outof)
   ;  int BBB_offset = pam->modes & MODE_QC_LATE ? get_bbb_offset(rec) : -1

   ;  int lint_offset = 0
   ;  int indel_allowed = pam->mr_3p_global.mr_maxgap > 0
   ;  SWNUM data[MATRIX_SIZE] = { 0 }
   ;  struct sw_alninfo ai = { 0 }
   ;  int buckid, id_clean = 0, n_clean_best = 0
   ;  const char* msg = "nomatch"
   ;  int score_best = 0, n_pass = 0

; const char* cfl[10] = { "cfl0", "cfl1",  "cfl2",  "cfl3",  "cfl4",  "cfl5",  "cfl6",  "cfl7",  "cfl8", "cflZ" }

   ;  for (buckid = 0; buckid < pam->buck_n; buckid++)
      {  int score = 0
                  /* NOTE:
                   * We used to do restriction of search space
                   * (e.g. rec->seq + pam->cleanlength, limit
                   * search to start of adaptor). That restriction was dropped,
                   * but the reasons for dropping were not noted unfortunately.
                  */
      ;  if
         (  sw_fill
            (  &ai
            ,  data
            ,  MATRIX_SIZE
            ,  rec->seq
            ,  pam->buck[buckid].cat3p
            ,  indel_allowed
            ,  &pam->swp_3p
         )  )
         X_FINISHED("sw error in do_barcode_3p")

      ;  rec->rs.N_cells += ai.ni * ai.nj
      ;  sw_trace_trim(&ai, &pam->swp_3p, ai.max_ij, 0)         /* do_barcode_3p */

      ;  lint_offset = 1
      ;  if ((ai.lft_start + 1) > ai.rgt_start)
         lint_offset = (ai.lft_start + 1) - ai.rgt_start

                                 /* donow: ship barcode offset and length */
      ;  if
         (  sw_match_bc3p
            (  &ai
            ,  &pam->mr_3p_global
            ,  &pam->mr_3p_prefix
            ,  &pam->mr_3p_barcode
            ,  pam->buck[buckid].sinsert3p_n
            ,  pam->buck[buckid].barcode_n
            ,  &score            /* ignored for now */
         )  )
         {  score = ai.n_match - ai.n_insl - ai.n_insr - ai.n_subst
         ;  n_pass++
         ;  if (score > score_best)
            {  id_clean = buckid
            ;  n_clean_best = 1
            ;  score_best = score
         ;  }
            else if (score == score_best)
            n_clean_best++
         ;  pam->buck[buckid].cache_score = score
         ;  pam->buck[buckid].cache_accept = 1
      ;  }
         else
            lint_offset = ai.ni
         ,  pam->buck[buckid].cache_score = score
         ,  pam->buck[buckid].cache_accept = 0

      ;  pam->buck[buckid].cache_lint_offset = lint_offset
   ;  }

      if (n_clean_best == 0)
      rec->rs.N_nomatch++

   ;  if (n_clean_best == 1 && NNN_offset > 0)
      {  rec->rs.B_NN += rec->seq_n+1 - NNN_offset
      ;  if (NNN_offset <= pam->cleanlength)
         {  n_clean_best = 0
         ;  rec->rs.N_nn++
         ;  msg = "NNN"
      ;  }
         else if (NNN_offset < lint_offset)
         lint_offset = NNN_offset
   ;  }

      if (n_clean_best == 1 && BBB_offset > 0)
      {  rec->rs.B_BB += rec->seq_n+1 - BBB_offset
      ;  if (BBB_offset <= pam->cleanlength)
         {  n_clean_best = 0
         ;  rec->rs.N_bb++
         ;  msg = "lowq"
      ;  }
         else if (BBB_offset < lint_offset)
         lint_offset = BBB_offset
   ;  }

               /* The check for pam->cleanlength is so that we do
                * not drop reads that have a virtual 'negative'
                * offset of the first base of aligned adaptor.
               */
      if
      (  n_clean_best == 1
      && pam->cleanlength
      && ai.lft_start + 1 <= pam->cleanlength + ai.rgt_start
      )
      {  n_clean_best = 0
      ;  msg = "len"
      ;  rec->rs.N_lenlo++
   ;  }

               /* Soft assertion: see comment below */
      if (n_clean_best == 1 && !lint_offset)
      {  argh("reaper", "nasty zero lint offset, please fix code")
      ;  lint_offset = 1
   ;  }

      if
      (  n_clean_best == 1
      && check_for_tabu(rec, lint_offset-1, pam, &pam->buck[id_clean].tabu2)

               /* warning: in the above check (and further below as well)
                * we need that lint_offset > 0.
                * this should be garantueed by having a match and
                * by restricting the search space, but
                * with an all-zero alignment matrix it could
                * still be icky. So there's an assertion above.
               */
      )
      {  n_clean_best = 0
      ;  rec->rs.N_tabu++
      ;  msg = "tabu"
   ;  }

      if (n_clean_best == 1)
      {  rec->clean_n = pam->buck[id_clean].cache_lint_offset - 1
      ;  rec->out_message = NULL
      ;  dump_read(1, pam->buck[id_clean].fp, rec, NULL, NULL, pam)
      ;  rec->rs.N_clean++
#if 0
      ;  if (n_pass > 1)
         rec->rs.N_cflr++
#endif
      ;  qc_tally_in(rec, pam, id_clean, rec->seq_n)
      ;  qc_tally_out(rec, pam, id_clean, rec->clean_n)
   ;  }
      else
      {  rec->out_message = n_clean_best ? cfl[n_clean_best > 9 ? 9 : n_clean_best] : msg
      ;  dump_read(0, pam->fplint, rec, NULL, NULL, pam)
      ;  qc_tally_out(rec, pam, pam->buck_n, rec->seq_n)    /* specialbucket */
      ;  rec->rs.N_lint++
      ;  if (n_clean_best > 1)
         rec->rs.N_cfl++
   ;  }

            /* in do_barcode_3p, buck[].fp, fplint  */
      if
      (  pam->verbose & VB_ALIGN
      && rec->rs.N % pam->sample_lint == 0
      )
      {  for (buckid = 0; buckid < pam->buck_n; buckid++)
         {  int score = 0
         ;  if
            (  sw_fill
               (  &ai
               ,  data
               ,  MATRIX_SIZE
               ,  rec->seq
               ,  pam->buck[buckid].cat3p
               ,  indel_allowed
               ,  &pam->swp_3p
            )  )
            X_FINISHED("sw error in do_barcode_3p")
         ;  rec->rs.N_cells += ai.ni * ai.nj
         ;  sw_trace_trim(&ai, &pam->swp_3p, ai.max_ij, 0)    /* do_barcode_3p */
         ;  sw_match_bc3p
            (  &ai
            ,  &pam->mr_3p_global
            ,  &pam->mr_3p_prefix
            ,  &pam->mr_3p_barcode
            ,  pam->buck[buckid].sinsert3p_n
            ,  pam->buck[buckid].barcode_n
            ,  &score         /* ignored */
            )
         ;  sw_pp2
            (  "read"
            ,  x_3pbc[pam->buck[buckid].cache_accept]
            ,  pam->buck[buckid].cache_score
            ,  &ai
            ,  pam->fplint
            ,  pam->zippit
            ,  rec->rs.N
            ,  0
            ,  NULL
            )
            /* in do_barcode_3p */
      ;  }
      }
      return 0
;  }


int n_zip_mismatch
(  const char* s1
,  const char* s2
)
   {  int n = 0
   ;  const char* a, *b
   ;  for (a=s1,b=s2; *a && *b; a++, b++)
      n += *a != *b
   ;  return n
;  }


int do_barcode_5p
(  struct record *rec
,  struct param *pam
)
   {  int buckid = 0
   ;  int bestid = -1
   ;  struct match_requirement mr = pam->mr_5p_barcode
   ;  int indel_allowed = mr.mr_maxgap > 0
   ;  int bestend = -1
   ;  int conflict = 0

   ;  SWNUM data[MATRIX_SIZE] = { 0 }
   ;  struct sw_alninfo ai = { 0 }

   ;  {  for (buckid = 0; buckid < pam->buck_n; buckid++)
         {  if (!strncmp(pam->buck[buckid].barcode5p, rec->seq, pam->buck[buckid].barcode_n))
            {  bestend = pam->buck[buckid].barcode_n
            ;  break
         ;  }
         }
         if (bestend > 0)
         bestid = buckid
   ;  }

               /* mr.mr_offset & 1 indicates a 'zip' alignment,
                * toe to toe and head to head aligned with only substitutions allowed.
               */
      if (bestid < 0 && mr.mr_maxedit > 0 && (mr.mr_offset & 1))
      {  for (buckid = 0; buckid < pam->buck_n; buckid++)
         {  int z = INT_MAX
         ;  if (pam->buck[buckid].barcode_n <= rec->seq_n)
            z = n_zip_mismatch(pam->buck[buckid].barcode5p, rec->seq)
#if WE_USE_ZLIB
;if(0)gzprintf(pam->fplint, "code %s seq %s zip %d\n", pam->buck[buckid].barcode5p, rec->seq, z)
#endif
         ;  if (z <= mr.mr_maxedit)
            {  if (bestid >= 0)
               {  conflict = 1
#if WE_USE_ZLIB
               ;  gzprintf
                  (  pam->fplint
                  ,  "> 5p-conflict-zip %s %s [%s/%d] recno %d\n"
                  ,  pam->buck[bestid].barcode5p
                  ,  pam->buck[buckid].barcode5p
                  ,  rec->seq
                  ,  (int) rec->seq_n
                  ,  (int) rec->rs.N
                  )
#endif
               ;  bestid = -1
               ;  break
            ;  }
               else
               bestid = buckid
         ;  }
         }
         if (!conflict && bestid >= 0)
         bestend = pam->buck[bestid].barcode_n
   ;  }


      if (!conflict && bestid < 0 && (mr.mr_offset & 2) && mr.mr_maxedit)
      {  int n_edit = 0
      ;  int n_gap = 0
      ;  for (buckid = 0; buckid < pam->buck_n; buckid++)
         {  int found_alignment = 0
         ;  pam->swp_5p.right_limit  = pam->buck[buckid].barcode_n + mr.mr_maxgap
         ;  if (rec->seq_n < pam->buck[buckid].barcode_n)
            continue
         ;  if
            (  sw_fill
               (  &ai
               ,  data
               ,  MATRIX_SIZE
               ,  pam->buck[buckid].barcode5p
               ,  rec->seq
               ,  indel_allowed
               ,  &pam->swp_5p                     /* governs alignment algorithm, settings are very stable */
            )  )
            X_FINISHED("sw error in do_barcode_5p")

         ;  rec->rs.N_cells += ai.ni * ai.nj

         ;  sw_trace_trim(&ai, &pam->swp_5p, ai.max_ij, 0)    /* do_barcode_5p */
         ;  n_gap = ai.n_insl + ai.n_insr  
         ;  n_edit = ai.n_subst + n_gap

         ;  found_alignment
            =  (  (ai.lft_end - ai.lft_start) + 1 + mr.mr_maxedit >= pam->buck[buckid].barcode_n + n_edit   /* + (ai.rgt_start - 1) */
               && mr.mr_maxedit >= n_edit
               && mr.mr_maxgap >= n_gap
               )

;if(0)fprintf(stderr, "%d\n", (int) ((ai.lft_end - ai.lft_start) + 1))
;if(0)fprintf(stderr, "found %d bc %s (%d %d)\n", (int) found_alignment, pam->buck[buckid].barcode5p, (int) n_edit, (int) (ai.rgt_start -1))
               /* in do_barcode_5p, buck[].fp, fplint */
         ;  if
            (  pam->verbose & VB_ALIGN
            && rec->rs.N % pam->sample_lint == 0
            )
            sw_pp2("read", x_5pbc[found_alignment], -1, &ai, pam->fplint, pam->zippit, rec->rs.N, 0, NULL)
               /* in do_barcode_5p */

         ;  if (found_alignment)
            {  if (bestid >= 0)
               {
#if WE_USE_ZLIB
;  gzprintf(pam->fplint, "> 5p-conflict-align %s %s %s\n", pam->buck[bestid].barcode5p, pam->buck[buckid].barcode5p, rec->seq)
#endif
               ;  bestid = -1
               ;  conflict = 1
               ;  break
            ;  }
               else
                  bestid = buckid
               ,  bestend = ai.rgt_end + (pam->buck[bestid].barcode_n - ai.lft_end)
         ;  }
         }
      }

      if (bestid < 0)
      {  rec->rs.N_lint++
      ;  rec->out_message = conflict ? "conflict" : "nobc"
      ;  qc_tally_out(rec, pam, pam->buck_n, rec->seq_n)    /* specialbucket */
      ;  dump_read(0, pam->fplint, rec, NULL, NULL, pam)
      ;  if (conflict)
         rec->rs.N_cfl++
      ;  else
         rec->rs.N_nomatch++
   ;  }
      else
      {  qc_tally_in(rec, pam, bestid, rec->seq_n)          /* tally in per-barcode */
#if WE_USE_ZLIB
;if(0)gzprintf(pam->fplint, "# %s %s\n", pam->buck[bestid].barcode5p, rec->seq)
#endif
      ;  rec_remove_prefix(rec, bestend)                    /* inclusive end in 1-based offset == length [12345] */

      ;  pam->sinsert_5p
         =    (pam->buck[bestid].sinsert5p_n)
            ?  pam->buck[bestid].sinsert5p
            :  NULL
      ;  pam->adaptor_3p = pam->buck[bestid].adaptor3p
      ;  pam->the_tabu2  = pam->buck[0].tabu2

      ;  X_CASCADE(ntnt(rec, pam, pam->buck[bestid].fp))

      ;  if (rec->clean_n >= 0)
         qc_tally_out(rec, pam, bestid, rec->clean_n)
   ;  }
      return 0
;  }



int parse_mr
(  const char* form
,  struct match_requirement* mr
,  int maxnums
)
   {  int n_consumed, ok = 1
   ;  mr->mr_maxgap = 0

   ;  do
      {  if
         (  maxnums >= 4
         && 4 == sscanf(form, "%u/%u/%u/%d%n",
            &mr->mr_minlen, &mr->mr_maxedit, &mr->mr_maxgap, &mr->mr_offset, &n_consumed)
         )
         break
      ;  if
         (  maxnums >= 3
         && 3 == sscanf(form, "%u/%u/%u%n",
            &mr->mr_minlen, &mr->mr_maxedit, &mr->mr_maxgap, &n_consumed)
         )
         break
      ;  if
         (  maxnums >= 2
         && 2 == sscanf(form, "%u/%u%n",
            &mr->mr_minlen, &mr->mr_maxedit, &n_consumed)
         )
         break
      ;  ok = 0
   ;  }
      while (0)

   ;  if (ok && mr->mr_maxedit < mr->mr_maxgap)
      {  argh
         (  "error"
         ,  "maximum gap count %d should not exceed maximum edit distance %d in %s"
         ,  ji mr->mr_maxgap
         ,  ji mr->mr_maxedit
         ,  form
         )
      ;  mr->mr_maxedit = mr->mr_maxgap
      ;  argh("fixit", "setting maximum edit distance to %d", ji mr->mr_maxedit)
   ;  }
      return ok && n_consumed == strlen(form)
;  }


int parse_anchor
(  const char* form
,  int *d
,  const char **anchor
)
   {  int n_consumed = 0
   ;  if (1 != sscanf(form, "%d/%n", d, &n_consumed))
      return 0
   ;  if (n_consumed == 0 || form[n_consumed-1] != '/')
      return 0
   ;  d[0] -= 1
   ;  anchor[0] = form + n_consumed
   ;  return 1
;  }


int parse_pair
(  const char* form
,  unsigned *u1
,  unsigned *u2
)
   {  int n_consumed
   ;  if (2 != sscanf(form, "%u/%u%n", u1, u2, &n_consumed))
      return 0
   ;  if (n_consumed != strlen(form))
      return 0
   ;  return 1
;  }


int parse_dust_spec
(  const char* form
,  unsigned *cutoff
,  unsigned *bases
)
   {  int n_converted = 0
   ;  if
      (   2 != sscanf(form, "%u/%n", cutoff, &n_converted) 
      &&  1 != sscanf(form, "%u", cutoff)
      )
      return 0
   ;  if (n_converted > 0)
      {  const char* o = form + n_converted
      ;  if (strchr(o, 'A')) bases[0] |= 1 << 0
      ;  if (strchr(o, 'C')) bases[0] |= 1 << 1
      ;  if (strchr(o, 'G')) bases[0] |= 1 << 2
      ;  if (strchr(o, 'T')) bases[0] |= 1 << 3
   ;  }
      return 1
;  }


static volatile sig_atomic_t g_abort_loop = 0;


void reaper_sig_catch
(  int sig
)
   {  if (sig == SIGINT)
      g_abort_loop = 1
;  }


int reaper_main
(  int argc
,  const char* argv[]
)
   {  ZFILE input = NULL

   ;  const char* g_fnin = "-"
   ;  char* fnclean = NULL, *fnlint = NULL
   ;  unsigned g_limit = 0

   ;  const char* g_format_in   =   "@%I%A%n%R%n+%#%Q%n"
   ;  const char* g_format_in2  =   NULL
   ;  int status = 19
   ;  int called_from_R = 0

                                                      /* below initialises read_stats_global (rec.rs) */
   ;  struct record rec = { { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, 0 }
   ;  struct param  pam = { 0 }

#ifdef WINDOWS_BUILD
   ;  set_argh("reaperlog.txt", 1)
#endif

   ;  param_init(&pam)
   ;  g_abort_loop = 0
   ;  signal(SIGINT, reaper_sig_catch)

   ;  themap_init()
   ;

/* enter macromagicalitaciousness */

      arg_switch()

      optarg("-i")            g_fnin = thearg();               endarg()
      optarg("-fastq")        g_fnin = thearg();               endarg()
      optarg("-clean-length") pam.cleanlength = atoi(thearg());  endarg()
      optarg("-polya")        pam.polya = atoi(thearg());      endarg()
      optarg("-dust-suffix")
         if (!parse_dust_spec(thearg(), &pam.taildust, &pam.taildust_bases))
         X_ERR_JUMP(DONE, "-dust-suffix argument needs <cut-off> or <cutoff>/<ACGT> pattern");
      endarg()
      optarg("-dust-suffix-late")
         if (!parse_dust_spec(thearg(), &pam.taildustlate, &pam.taildustlate_bases))
         X_ERR_JUMP(DONE, "-dust-suffix-late argument needs <cut-off> or <cutoff>/<ACGT> pattern");
      endarg()
      optarg("-qqq-check")
         if
         (  4 != sscanf(thearg(), "%u/%u/%u/%u", &pam.q_cutoff, &pam.q_winlen, &pam.q_winofs, &pam.q_readofs)
         && 3 != sscanf(thearg(), "%u/%u/%u", &pam.q_cutoff, &pam.q_winlen, &pam.q_winofs)
         && 2 != sscanf(thearg(), "%u/%u", &pam.q_cutoff, &pam.q_winlen)
         )
         X_ERR_JUMP(DONE, "-medq argument needs <cutoff>/<winlen>[/<winoffset>] pattern");
         if (pam.q_winofs > pam.q_winlen) pam.q_winofs = pam.q_winlen;
         if (pam.q_winlen > 255)          pam.q_winlen = 255;
         if (pam.q_winlen % 2 == 0)       pam.q_winlen++;
         if (!pam.q_winofs)               pam.q_winofs = (pam.q_winlen + 1) / 2;
      endarg()
      optarg("-anchor")
         if (!parse_anchor(thearg(), &pam.anchor_offset, &pam.anchor_string))
         X_ERR_JUMP(DONE, "-anchor argument does not follow <OFFSET>/<BASES> pattern");
         fprintf(stderr, "filtering on substring %s at position %d\n", pam.anchor_string, pam.anchor_offset+1);
      endarg()
      optarg("-tri")
         pam.trintco = atoi(thearg());
         if (pam.trintco < 0 || pam.trintco > 100) {
            argh("reaper", "-trint option should take argument in [0..100]");
            pam.trintco = 0;
         }
      endarg()
      optarg("-mr-tabu")
         if (!parse_mr(thearg(), &pam.mr_tabu, 4))
         X_ERR_JUMP(DONE, "-mr-tabu argument needs x/y or x/y/z or x/y/z/u pattern");
      endarg()
      optarg("-5p-barcode")
         if (!parse_mr(thearg(), &pam.mr_5p_barcode, 4))
         X_ERR_JUMP(DONE, "-5p-barcode argument needs x/y or x/y/z or x/y/z/u pattern");
      endarg()
      optarg("-5p-sinsert")
         if (!parse_mr(thearg(), &pam.mr_5p_require, 4))
         X_ERR_JUMP(DONE, "-5p-sinsert argument needs x/y or x/y/z or x/y/z/u pattern");
      endarg()
      optarg("-3p-global")
         if (!parse_mr(thearg(), &pam.mr_3p_global, 4))
         X_ERR_JUMP(DONE, "-3p-global argument needs x/y or x/y/z or x/y/z/u pattern");
      endarg()
      optarg("-3p-prefix")
         if (!parse_mr(thearg(), &pam.mr_3p_prefix, 4))
         X_ERR_JUMP(DONE, "-3p-prefix argument needs x/y or x/y/z or x/y/z/u pattern");
      endarg()
      optarg("-3p-barcode")
         if (!parse_mr(thearg(), &pam.mr_3p_barcode, 3))
         X_ERR_JUMP(DONE, "-3p-barcode argument needs x/y or x/y/z pattern");
      endarg()
      optarg("-3p-head-to-tail")   pam.p_3p_tth = atoi(thearg());    endarg()    /* suffix perfect match */
      optarg("-swp")
         if (3 != sscanf(thearg(), "%u/%u/%u",
            &pam.swp_3p.gain_match, &pam.swp_3p.cost_subst, &pam.swp_3p.cost_gapleft)
            )
         X_ERR_JUMP(DONE, "-swp argument needs MATCH/SUBS/GAP pattern (gain/cost/cost)");
         pam.swp_3p.cost_gapright = pam.swp_3p.cost_gapleft;
         pam.swp_5p    =  pam.swp_3p;
         pam.swp_tabu  =  pam.swp_3p;
         pam.swp_raptor=  pam.swp_3p;
      endarg()

      optarg("-loco")
         pam.loco_check = atof(thearg());
      endarg()

      optarg("-geom")
         if (!strcmp(thearg(), "3p-bc"))
         g_mode = CASE_3P_BC
      ;  else if (!strcmp(thearg(), "5p-bc"))
         g_mode = CASE_5P_BC
      ;  else if (!strcmp(thearg(), "no-bc"))
         g_mode = CASE_NOBC
      ;  else
         X_ERR_JUMP(DONE, "unsupported geometry [%s]", thearg())
      ;
      endarg()

      optarg("-meta")    pam.fnmeta = thearg();           endarg()
      optarg("-record-format2")    g_format_in2 = thearg();       endarg()
      optarg("-record-format")     g_format_in = thearg();      endarg()
      uniarg("--fasta-in")         g_format_in = ">%I%A%n%R%n";   endarg()
      uniarg("--fastqx-out")
         pam.format_clean  = "@%I%trecno=%J%n%C%n+%n%Q%n";
         pam.format_lint   = "@%I%trecno=%J%tmsg=%M%n%R%n+%n%Y%n";
      endarg()
      uniarg("--fasta-out")
         pam.format_clean  = ">%I%n%C%n";
         pam.format_lint   = ">%I%tmsg=%M%n%R%n";
      endarg()
      uniarg("--fastax-out")
         pam.format_clean  = ">%I%trecno=%J%n%C%n";
         pam.format_lint   = ">%I%trecno=%J%tmsg=%M%n%R%n";
      endarg()
      optarg("-sc-max")       pam.suffix_complexity_max = atof(thearg());  endarg()
      optarg("-nnn-check")
         if (!parse_pair(thearg(), &pam.nnn_count, &pam.nnn_outof))
         X_ERR_JUMP(DONE, "-nnn-check argument does not follow x/y pattern");
      endarg()
      optarg("-sample")
         if (!parse_pair(thearg(), &pam.sample_clean, &pam.sample_lint))
         X_ERR_JUMP(DONE, "-sample argument does not follow x/y pattern");
      endarg()
      optarg("-debug")
         if (strchr(thearg(), 'a')) pam.verbose |= VB_ALIGN;
         if (strchr(thearg(), 'c')) pam.verbose |= VB_CLEAN;
         if (strchr(thearg(), 'l')) pam.verbose |= VB_LINT;
      endarg()
      optarg("-do")           g_limit = atoi(thearg());        endarg()
      uniarg("--nozip")       pam.zippit = 0;                  endarg()
      uniarg("--R")           called_from_R = 1;               endarg()
      uniarg("--noqc")        BIT_OFF(pam.modes, MODE_QC_FULL);endarg()
      uniarg("--bcq-early")   BIT_ON(pam.modes, MODE_QC_EARLY);endarg()
      uniarg("--bcq-late")    BIT_ON(pam.modes, MODE_QC_LATE); endarg()
      uniarg("--full-length") BIT_ON(pam.modes, MODE_FULLLENGTH ); endarg()
      uniarg("--keep-all")    BIT_ON(pam.modes, MODE_KEEPALL ); endarg()
      optarg("-3pa")          pam.adaptor_3p = hte(thearg());       endarg()
      optarg("-raptor")       pam.adaptor_test = thearg();       endarg()
      optarg("-kmer")         pam.kmer = atoi(thearg());       endarg()
      optarg("-5psi")         pam.sinsert_5p = hte(thearg());       endarg()
      optarg("-tabu") parse_tabu(&pam.the_tabu2, hte(thearg()), g_restrict);    endarg()
      optarg("-basename")     pam.basename = thearg();         endarg()
      optarg("-format-clean") pam.format_clean = thearg();     endarg()
      optarg("-length-limit") pam.length_limit = atoi(thearg());     endarg()
      optarg("-format-lint")  pam.format_lint  = thearg();      endarg()
      optarg("-trim-length")  pam.length_co    = atoi(thearg());     endarg()

   uniarg("--stream-clean")   BIT_ON(pam.modes, MODE_STREAMCLEAN); endarg()
   uniarg("--record-format")
fprintf(stdout, "-record-format <format string> (extended record description syntax)\n");
puts("     %R  expect read (longest sequence found over [a-zA-Z]*) - empty read allowed");
puts("     %I  expect identifier (longest sequence of non-blank)");
puts("     %A  expect overflow annotation field (initial blanks skipped, rest of line)");
puts("     %a  use previous %F or %G or %H field as annotation");
puts("     %B  expect annotation field (longest sequence of non-blank)");
puts("     %Q  expect quality (longest sequence of non-blank)");
puts("     %X  expect count (a nonnegative integer number)");
puts("     %F  expect and discard field (longest sequence of non-tab)");
puts("     %G  expect and discard field (longest sequence of non-blank)");
puts("     %Hx (x is a placeholder) expect and discard up to x");
puts("     %#  discard everything until end of line");
puts("     %b  expect run of blanks (space or tab)");
puts("     %n  expect end of line match");
puts("     %.  expect and discard any character");
puts("     %s  expect a space");
puts("     %t  expect a tab");
puts("     %%  expect a percent sign");
puts("Anything else requires a literal match, see the FASTA example just below");
puts("The string '>%I%#%R%n' will parse FASTA-ish input (but use --fasta-in instead)");
puts("The string '%F%t%R%#' retrieves the second field from each line");
status = 0;
goto DONE;
   endarg()

   optarg("-restrict")
      g_restrict = atoi(thearg());
   endarg()

   optarg("-guard")
      g_guard = atoi(thearg());
   endarg()

   uniarg("-z")
fprintf(stdout, "Struct param:  %d\n", (int) sizeof(struct param));
fprintf(stdout, "Struct bucket: %d\n", (int) sizeof(struct bucket));
status = 0;
goto DONE;
   endarg()

   uniarg("--version")
fprintf(stdout, "Reaper version: %s\n", reaper_tag);
status = 0;
goto DONE;
   endarg()

   uniarg("-h")
puts("\nRequired options");
puts("-geom <mode>   mode in {no-bc, 3p-bc, 5p-bc}");
puts("-meta <fname>  file with geometry-dependent format. Required columns:");
puts("   Geometry    Columns:");
puts("      no-bc          3p-ad     -       -      -    tabu");
puts("      3p-bc          3p-ad  barcode  3p-si    -    tabu");
puts("      5p-bc          3p-ad  barcode    -    5p-si  tabu");
puts("bc=barcode, ad=adaptor, si=sequence insert");
puts("Columns 3p-si, 5p-si, 3p-ad and tabu may all be empty");
puts("Alternatively, to express absence, a single hyphen may be used");

puts("\nImportant options");
puts("-i <fname>        input stream (gzipped file allowed) (default STDIN)");
fprintf(stdout, "-clean-length <int> minimum allowed clean length (default %d)\n", (int) pam.cleanlength);
puts("-guard <int>      protect first <int> bases in read from adapter and tabu matching");
fprintf(stdout, "-restrict <int>   only use the first <int> bases of adapter or tabu sequence (default %d)\n", (int) g_restrict);
puts("                  This is to avoid false positive matches");
puts("-tri <threshold>  filter out reads with tri-nt score > threshold");
puts("                  a reasonable <threshold> is 35");
fprintf(stdout, "-qqq-check  <val>/<winlen>  cut sequence when median quality drops below <val>\n");
fprintf(stdout, "-qqq-check  <val>/<winlen>/<winofs> as above, cut at <winofs> (default %d)\n", (int) pam.q_winofs);
fprintf(stdout, "-qqq-check  <val>/<winlen>/<winofs>/<readofs> as above, start at <readofs>\n");
fprintf(stdout, "-dust-suffix <threshold> dust theshold for read suffix (default %d, suggested 20)\n", (int) pam.taildust);
fprintf(stdout, "-nnn-check <count>/<outof> (default %d/%d)\n", (int) pam.nnn_count, (int) pam.nnn_outof);
puts("      disregard read onwards from seeing <count> N's in <outof> bases");

puts("\nAlignment options");
puts("Options to specify when part of an alignment triggers a match:");
#define local_EXPAND_MR(mr) (int) mr.mr_minlen, (int) mr.mr_maxedit, (int) mr.mr_maxgap, (int) mr.mr_offset
fprintf(stdout, "-3p-global  l/e[/g[/o]]  (default %d/%d/%d/%d)\n", local_EXPAND_MR(pam.mr_3p_global));
fprintf(stdout, "-3p-prefix  l/e[/g[/o]]  (default %d/%d/%d/%d)\n", local_EXPAND_MR(pam.mr_3p_prefix));
fprintf(stdout, "-3p-barcode l/e[/g[/o]]  (default %d/%d/%d/%d)\n", local_EXPAND_MR(pam.mr_3p_barcode));
fprintf(stdout, "-5p-barcode l/e[/g[/o]]  (default %d/%d/%d/%d)\n", local_EXPAND_MR(pam.mr_5p_barcode));
fprintf(stdout, "-5p-sinsert l/e[/g[/o]]  (default %d/%d/%d/%d)\n", local_EXPAND_MR(pam.mr_5p_require));
fprintf(stdout, "-mr-tabu    l/e[/g[/o]]  (default %d/%d/%d/%d)\n", local_EXPAND_MR(pam.mr_tabu));
fprintf(stdout, "-3p-head-to-tail l minimal trailing perfect match length (default %d)\n", (int) pam.p_3p_tth);
#undef local_EXPAND_MR
puts("   syntax used in the above:");
puts("      l  <int> minimum length required to count sub-alignment as match");
puts("      e  <int> maximum allowed edit distance");
puts("      g  <int> [optional, not active when set to 0] maximum allowed number of gaps");
puts("      o  <int> [optional, not active when set to 0] offset:");
puts("            o= 5  requires alignment to start in the first five bases of adaptor");
puts("            o=-5  requires alignment to end in the last five bases of adaptor");
fprintf(stdout, "-swp M/S/G match/substitution/gap gain/cost/cost (default %u/%u/%u)\n",
   pam.swp_3p.gain_match, pam.swp_3p.cost_subst, pam.swp_3p.cost_gapleft);
puts("\nInput/output options");
puts("--fasta-in      read FASTA input");
fprintf(stdout, "-record-format <format string> (record description, default %s)\n", g_format_in);
puts("   [ -record-format syntax is output when supplying --record-format ]");
puts("-record-format2 <format string> (simple line formats, one field per line):");
puts("      R  read");
puts("      I  read identifier");
puts("      Q  quality scores");
puts("      D  discard field");
puts("");
puts("-basename <pfx>   pfx.lint.gz, pfx.clean.gz pfx.report etc will be constructed");
puts("-format-clean <format string> (output for clean reads)");
puts("-format-lint <format string> (output for filtered reads)");
puts("   -format-clean/lint specification syntax:");
puts("      %R  read");
puts("      %C  clean read");
puts("      %Z  clean read padded with Ns if necessary");
puts("      %V  reverse complement of clean read");
puts("      %I  read identifier");
puts("      %Q  clean or input read quality (for clean / lint file respectively)");
puts("      %X  read count (only applicable if -record-format is used)");
puts("      %Y  input read quality");
puts("      %q<c>  clean input read quality padded with character <c>");
puts("      %A  annotation field");
puts("      %L  clean read length");
puts("      %M  message describing cause for filtering (lint file)");
puts("      %T  trinucleotide complexity score (clean/lint file)");
puts("      %U  dUst sUffix complexity information");
puts("      %3  best read/3p-adaptor alignment");
puts("      %=  alignment characteristics");
puts("            mt=matchtype");
puts("            sc=suffix-complexity");
puts("            ht=head-tail-match");
puts("            nn=N-match-offset");
puts("            bb=B-match-offset");
puts("            aa=Polya-offset");
puts("            qq=Quality-trim-offset");
puts("      %n  newline");
puts("      %J  record offset, unique for each read. Use to match paired-end reads");
puts("      %f  fastq line number based on standard fastq format");
puts("      %t  tab");
puts("      %%  percent sign");
puts("   Anything else is copied verbatim");
puts("-debug [acl]+     a=alignments c=clean l=lint");
puts("-sample c/l       if debug, sample every c/l clean/lint read");
puts("--nozip           do not output gzipped files");
puts("--noqc            do not output quality report files");

puts("\nMiscellaneous options");
puts("--bcq-early       perform early 'B' quality filtering (when reading sequences)");
puts("--bcq-late        perform late 'B' quality filtering (before outputting sequences)");
puts("--full-length     only allow reads not shortened in any filter step");
puts("--keep-all        delete rather than discard reads (e.g. tabu match, missing 5p-sinsert)");
puts("-trim-length <int>     cut reads back to length <int>");
fprintf(stdout, "-polya <int>      remove trailing A's if length exceeds <int>\n");
fprintf(stdout, "-sc-max <f>       threshold for complexity of suffix following prefix match (default %.2f)\n", pam.suffix_complexity_max);
puts("\nOptions for use when running reaper with -geom no-bc");
puts("-3pa <three prime adaptor>");
puts("-5psi <five prime sequence insert>");
puts("-tabu <tabu sequence>");
puts("");
status = 0;
X_ERR_JUMP(DONE, "DON'T PANIC");
   endarg()

      failarg()
      arg_done()

/* curtains macromagicalitaciousness */


   ;  {  unsigned dummy
      ;  kraken_readline(NULL, NULL, 0, NULL, &dummy)      /* resets static buffers */
   ;  }

      if (called_from_R)
      argh("reaper", "R is calling")

   ;  if (g_mode == CASE_NOTSPECIFIED)
      X_ERR_JUMP(DONE, "-geom option required (see -h for available options)")
   ;  else if (g_mode > CASE_NOBC && !pam.fnmeta)
      X_ERR_JUMP(DONE, "-meta option required with geometry %s", G_MODE[g_mode])

                           /* fixme: the logic of barcodes_io_open
                            * versus fpclean is somewhat entangled
                            * and dispersed throughout.
                           */
   ;  if (pam.fnmeta)
      {  if (config_read(pam.fnmeta, &pam, g_mode, g_restrict))
         X_ERR_JUMP(DONE, "failed to parse barcode spec file")

      ;  if
         (  g_mode != CASE_NOBC
         && barcodes_io_open(pam.buck, pam.basename, pam.buck_n, pam.zippit)
         )
         X_ERR_JUMP(DONE, "error opening barcode bucket files")

      ;  X_TRY_JUMP(DONE, buckets_prepare_qc(pam.buck, pam.buck_n, g_mode))

      ;  if (g_mode == CASE_3P_BC)
         {  pam.swp_3p.left_skip  = pam.cleanlength
         ;  pam.swp_3p.right_limit = pam.mr_3p_global.mr_minlen       /* fixme document */
         ;  pam.swp_3p.flags |= SW_ENDTOSTART
         ;  pam.p_3p_tth = 0
      ;  }
                                          /* retireme/fixme: the below is spurious given call above. */
#if 0
                                                            /* fixme: this is a bit of a kludge: separate tally/qc stuff  */
                                          ;  if (g_mode == CASE_NOBC)
                                             X_TRY_JUMP(DONE, buckets_prepare_qc(pam.buck, pam.buck_n, g_mode))
#endif
   ;  }
      else if (g_mode == CASE_NOBC)
      {  X_TRY_JUMP(DONE, buckets_alloc(&pam, 1))
      ;  X_TRY_JUMP(DONE, adaptors_set_case_3p(pam.buck, &pam.buck_n, &pam))
      ;  X_TRY_JUMP(DONE, buckets_prepare_qc(pam.buck, pam.buck_n, g_mode))
                                          /* note: in this branch we do not open io for (the one) bucket */
   ;  }
      else
      X_ERR_JUMP(DONE, "-meta option is required for all modes except '3p-ad'")


   ;  fnclean
      =  stringle
         (  "%s.%s.%s"
         ,  pam.basename
         ,  (g_mode == CASE_NOBC) ? "lane" : "ocean"    /* ocean case should not materialise */
         ,  pam.zippit ? "clean.gz" : "clean"
         )
   ;  fnlint = stringle("%s.%s", pam.basename, pam.zippit ? "lint.gz" :  "lint")
   ;

      if (!(input = myfopen(g_fnin, "r", 1)))
      X_ERR_JUMP(DONE, "bailing out")

   ;  if (g_mode == CASE_NOBC)
      {  if (pam.modes & MODE_STREAMCLEAN)
         fnclean = stringle("-")
      ;  if (!(pam.fpclean = myfopen(fnclean, ZWMODE, pam.zippit)))
         X_ERR_JUMP(DONE, "bailing out")
   ;  }

      if (!(pam.fplint = myfopen(fnlint, ZWMODE, pam.zippit)))
      X_ERR_JUMP(DONE, "bailing out")

   ;  {  clock_t t1 = clock(), t2
      ;  struct read_stats_global rsprev = { 0 }     /* read stats per million */
      ;  int anchor_len = 0
      ;  if (pam.anchor_string)
         anchor_len = strlen(pam.anchor_string)

      ;  fprintf(stderr, "---\n")
      ;  fprintf(stderr, "mRpm   million reads per minute\n")
      ;  fprintf(stderr, "mNpm   million nucleotides per minute\n")
      ;  fprintf(stderr, "mCps   million alignment cells per second\n")
      ;  fprintf(stderr, "lint   total removed reads (per 10K), sum of columns to the left\n")
      ;  fprintf
         (  stderr
         ,  "%-37s %7s %3s %4s %4s %4s %6s%3s %4s %4s %4s %4s %4s %4s %4s %4s %4s%s\n"
         ,  "25K reads per dot, 1M reads per line"
         ,  "seconds",  "mr", "mRpm", "mNpm", "mCps"
         ,  "{error",  "qc",  "low", "len", "NNN", "tabu", "nobc", "cflr", "cfl", "lint", "OK", "} per 10K"
         )
      ;  do
         {  unsigned rrstat =    g_format_in2
                              ?  read_record2(input, g_format_in2, &rec, pam.length_limit)
                              :  read_recordx(input, g_format_in, &rec)

         ;  if (rec.seq_n > MATRIX_SIZE / 24)
               rec.seq_n = MATRIX_SIZE / 24
            ,  rec.seq[rec.seq_n] = '\0'

         ;  if (rrstat == RR_DONE)
            break
         ;  else if (rrstat == RR_NOMEM || rrstat == RR_ERROR)
            goto DONE

         ;  if (pam.modes & MODE_FULLLENGTH)
            pam.cleanlength = rec.seq_n

         ;  rec.rs.B += rec.seq_n

         ;  do
            {  qc_tally_in(&rec, &pam, pam.buck_n, rec.seq_n)        /* specialbucket */
            ;  if (pam.anchor_offset >= 0 && pam.anchor_string)
               {  if
                  (  pam.anchor_offset <= rec.seq_n
                  && strncmp(rec.seq+pam.anchor_offset, pam.anchor_string, anchor_len)
                  )
                  {  rec.rs.N_anchor++
                  ;  rec.rs.N_lint++
                  ;  break
               ;  }
               }
               if (pam.modes & MODE_QC_EARLY)
               {  int bbb_ofs
               ;  if ((bbb_ofs = get_bbb_offset(&rec)) > 0)
                  {  rec.rs.B_BB += rec.seq_n + 1 - bbb_ofs
                  ;  if (bbb_ofs <= pam.cleanlength)
                     {  rec.rs.N_bb++; rec.rs.N_lint++; break ; }
                     else
                     {  rec.seq_n = bbb_ofs-1
                     ;  rec.q_n = bbb_ofs-1
                     ;  rec.seq[bbb_ofs-1] = '\0'
                     ;  rec.q[bbb_ofs-1] = '\0'
                  ;  }
                  }
               }

               if (pam.loco_check && loco(&rec, &pam, pam.loco_check))
                  {  rec.rs.N_loco++; rec.rs.N_lint++; break ;  }

               switch(g_mode)
               {  case CASE_NOBC
                : do_adaptor_3p(&rec, &pam); break
                ;
                  case CASE_3P_BC
                : do_barcode_3p(&rec, &pam); break
                ;
                  case CASE_5P_BC
                : do_barcode_5p(&rec, &pam); break
                ;
                  case CASE_NOTSPECIFIED
                : X_ERR_JUMP(DONE, "impossible")
            ;  }
               /* pam->verbose & VB_ALIGN */
         ;  }
            while (0)

         ;  if (25000 * (rec.rs.N / 25000) == rec.rs.N)
            {  fputc('.', stderr)
            ;  if (rec.rs.N % 1000000 == 0)
               {  double nsecs
               ;  t2 = clock()
               ;  unsigned diffn  = (rec.rs.N - rsprev.N)
               ;  unsigned diffnt = (rec.rs.B - rsprev.B)
               ;  unsigned diffsw = (rec.rs.N_cells - rsprev.N_cells)
               ;  nsecs = (t2 - t1) * 1.0 / CLOCKS_PER_SEC
               ;  fprintf
                  (  stderr
                  ,  " %4.0f %3d %4.1f %4.0f %4.0f"
                  ,  nsecs
                  ,  (int) (rec.rs.N / 1000000)
                  ,  (double) (60.0 / nsecs)
                  ,  (double) (diffnt * 1.0 / (nsecs * 1000000.0 / 60.0))
                  ,  (double) (diffsw * 1.0 / (nsecs * 1000000.0))
                  )
               ;  fprintf
                  (  stderr
                  ,  " %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d\n"
                  ,  (int)  (10000.5 * (rec.rs.N_err - rsprev.N_err) * 1.0 / diffn)
                  ,  (int)  (10000.5 * (rec.rs.N_bb  - rsprev.N_bb ) * 1.0 / diffn)
                  ,  (int)  (10000.5 * (rec.rs.N_loco - rsprev.N_loco) * 1.0 / diffn)
                  ,  (int)  (10000.5 * (rec.rs.N_lenlo - rsprev.N_lenlo) * 1.0 / diffn)
                  ,  (int)  (10000.5 * (rec.rs.N_nn - rsprev.N_nn) * 1.0 / diffn)
                  ,  (int)  (10000.5 * (rec.rs.N_tabu - rsprev.N_tabu) * 1.0 / diffn)
                  ,  (int)  (10000.5 * (rec.rs.N_nomatch - rsprev.N_nomatch) * 1.0 / diffn)
                  ,  (int)  (10000.5 * (rec.rs.N_cflr - rsprev.N_cflr) * 1.0 / diffn)
                  ,  (int)  (10000.5 * (rec.rs.N_cfl - rsprev.N_cfl) * 1.0 / diffn)
                  ,  (int)  (10000.5 * (rec.rs.N_lint - rsprev.N_lint) * 1.0 / diffn)
                  ,  (int)  (10000.5 * (rec.rs.N_clean - rsprev.N_clean) * 1.0 / diffn)
                  )
               ;  t1 = t2
               ;  rsprev = rec.rs
            ;  }
            }
         }
         while (!g_abort_loop && (!g_limit || rec.rs.N < g_limit))
      ;  fputc('\n', stderr)

      ;  if
         (     rec.rs.N_bb    +  rec.rs.N_loco
            +  rec.rs.N_lenlo +  rec.rs.N_nn
            +  rec.rs.N_tabu  +  rec.rs.N_cfl 
            +  rec.rs.N_5p_nosinsert  + rec.rs.N_nomatch
            +  rec.rs.N_aa
            +  rec.rs.N_tri   +  rec.rs.N_qq
            +  rec.rs.N_anchor
         != rec.rs.N_lint
         )
         argh("reaper", "(not overly important) discarded category counts do not add up to lint total")

      ;  if (rec.rs.N_clean + rec.rs.N_lint + rec.rs.N_err != rec.rs.N)
         argh("reaper", "(not overly important) dispatch category counts do not add up to input total")

      ;  if (rec.rs.N_nomsg)
         argh("reaper", "(not overly important) there were %u discarded reads with no message", rec.rs.N_nomsg)

      ;  qc_brief(&pam, &rec, g_mode)

      ;  if (g_fieldtosmall)
         argh("reaper", "___ WARNING ___ %lu cases where field was truncated", (long unsigned) g_fieldtosmall)

      ;  if (pam.modes & MODE_QC_FULL)
         qc_report(&pam, &rec, g_mode)
   ;  }

      status = 0
   ;  DONE:

   ;  myfzclose(input, 1)

   ;  myfzclose(pam.fpclean, pam.zippit)
   ;  myfzclose(pam.fplint, pam.zippit)

   ;  barcodes_io_close(pam.buck, pam.buck_n, pam.zippit)

   ;  if (0) data_release(&pam)

   ;  if (fnclean) free(fnclean)
   ;  if (fnlint) free(fnlint)

   ;  close_argh()

   ;  return status
;  }



#ifndef BUILD_R_BINDINGS

int main
(  int argc
,  const char* argv[]
)
   {  return reaper_main(argc, argv)
;  }

#endif



#if 0
#ifdef RREAPERR


SEXP R_dispatchee
(  SEXP list
,  int (*themain)(int , const char* [])
,  const char* callee
)
   {  /* SEXP elmt = R_NilValue */
   ;  SEXP result
   ;  SEXP names = getAttrib(list, R_NamesSymbol)
   ;  int argc = 1
   ;  const char* argv[100] = { 0 }

   ;  argv[0] = "dummy-internal"

   ;  int *res = NULL
   ;  int i = 0

   ;  PROTECT(result = allocVector(INTSXP, 1))
   ;  res = INTEGER(result)

   ;  if (length(list) <= 49)                   /* space for argv[0] and final NULL */
      for (i = 0; i < length(list); i++)
      {  const char* key = CHAR(STRING_ELT(names, i))
      ;  const char* val = CHAR(STRING_ELT(list, i))
      ;  if (!strlen(key))
         {  if (!strncmp(val, "--", 2))
            argv[argc++] = val
         ;  else
            {  argh(callee, "argument %s is not of type \"--mode\"", val)
            ;  break
         ;  }
         }
         else
         {  if (!strncmp(key, "--", 2))
            {  if (strlen(val) && strcmp(val, "on") && strcmp(val, "ON"))
               {  argh(callee, "option %s=%s is not recognised as \"--mode=on\"", key, val)
               ;  break
            ;  }
               argv[argc++] = key
         ;  }
            else
            {  argv[argc++] = key
            ;  argv[argc++] = val
         ;  }
         }
      }

   ;  argv[argc] = NULL

   ;  if (i != length(list))
      *res = 29
   ;  else
      {  int j
      ;  fprintf(stderr, "Passing to %s:", callee)
      ;  for (j=0;j<argc;j++)
         fprintf(stderr, " %s", argv[j])
      ;  fputc('\n', stderr)
      ;  *res=themain(argc, argv)
   ;  }
      fputc('\n', stderr)

   ;  UNPROTECT(1)
   ;  return result
;  }



SEXP reaperC
(  SEXP list
)
   {  return R_dispatchee(list, reaper_main, "reaper")
;  }



#endif
#endif








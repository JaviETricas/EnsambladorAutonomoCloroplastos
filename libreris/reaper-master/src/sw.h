
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


#ifndef  include_sw
#define  include_sw

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/* CAVEAT:
   Caller supports missing q (and the q we get might be garbage, although
   certain to be in array positions that exist but could have 0 bytes etc),
   but we don't check this. It's callers responsibility.
*/


/* the alignment matrix is layed out as follows,

                x    0  1  2  .. .. ..  nj-1 ("right")
                  +----^--^--^--------^---
               0  |  0  0  0  0  ..  0  0
               1  |  0  .  .. .. .. .. ..
               2  |  0  :  .           ..
               :  |  0  :     .        ..
               :  |  0  :        .     ..
               :  |  0  :           .  ..
   ("left")  ni-1 |  0  .  .. .. .. .. ..

      nj == strlen(right) + 1
      ni == strlen(left) + 1

   In memory this is a concatenation of rows.  The first row and column are set
   to all-zero.  sometimes we use elem(ai, i, j) (in sw.c), but other times we
   use ai->data[ij] and compute dependent i and j as

      i = ij / ai->nj
      j = ij - (i * ai->nj)

   Fixme:
      document left_skip and right_limit
*/



#define STRINGIFY0(v) #v
#define STRINGIFY(v) STRINGIFY0(v)

                  /* defining SWNUM as double has been the default, so well tested.
                   * 'unsigned' is almost two times faster.
                  */
#ifndef SWNUM
#     define SWNUM double
#endif




struct sw_alninfo
{  unsigned lft_start   /* NOTE this is in alignment space: first char corresponds to position 1 */
;  unsigned lft_end     /* NOTE inclusive */
#define SW_LFT_OHL(ai)  ((ai)->lft_start - 1)
#define SW_LFT_OHR(ai)  ((ai)->ni - 1 - (ai)->lft_end)
#define SW_LFT_OHLR(ai) (SW_LFT_OHL(ai) + SW_LFT_OHR(ai))
#define SW_RGT_OHL(ai)  ((ai)->rgt_start - 1)
#define SW_RGT_OHR(ai)  ((ai)->nj - 1 - (ai)->rgt_end)
#define SW_RGT_OHLR(ai) (SW_RGT_OHL(ai) + SW_RGT_OHR(ai))
#define SW_OVERHANG(ai) (SW_LFT_OHLR(ai) + SW_RGT_OHLR(ai))
;  unsigned rgt_start   /* NOTE this is in alignment space: first char corresponds to position 1 */
;  unsigned rgt_end     /* NOTE inclusive */
;  unsigned n_match
;  unsigned n_subst
;  unsigned n_insl
;  unsigned n_insr
;  SWNUM* data
;  const char* left
;  const char* right
;  const char* q
;  unsigned ni
;  unsigned nj
;  unsigned max_ij       /* i in left, j in right with max value in matrix */
;  unsigned max_ej       /* i is last in left, j somewhere in right with max value in row e(nd) */
#define SW_ALN_MAXSIZE 255
;  unsigned char aln_lft[SW_ALN_MAXSIZE+1]
;  unsigned char aln_aln[SW_ALN_MAXSIZE+1]
;  unsigned char aln_rgt[SW_ALN_MAXSIZE+1]
;  unsigned aln_ofs      /* start of alignment, inclusive */
;  unsigned aln_end      /* end of alignment, exclusive */
;  unsigned n_cells_used

;  unsigned aln_identity               /* per base, in aligned stretch, best of ref and query */
;  int aln_overhang
;  int tag_a             /* encode user information */
;  int tag_b             /* encode user information */
;  int tag_c             /* encode user information */
;  int tag_d             /* encode user information */
;
}  ;


#define SW_VERBOSE        ( 1 <<  0 )
#define SW_ENDTOSTART     ( 1 <<  1 )
#define SW_INDELALLOWED   ( 1 <<  2 )
#define SW_TRIMLEFT       ( 1 <<  3 )
#define SW_TRIMRIGHT      ( 1 <<  4 )
#define SW_TRIM           (SW_TRIMLEFT | SW_TRIMRIGHT)
#define SW_NW_FILL        ( 1 <<  5 )
#define SW_NW_TRACE       ( 1 <<  6 )
#define SW_NW_CODE        ( 1 <<  7 )


struct sw_param
{  int cost_gapleft
;  int cost_gapright
;  int gain_match
;  int cost_subst
;  unsigned left_skip
;  unsigned right_limit
;  unsigned flags
;
}  ;


unsigned sw_fit
(  struct sw_alninfo* ai
,  unsigned ij
)  ;


int sw_fill
(  struct sw_alninfo* ai
,  SWNUM* data
,  unsigned data_size
,  const char *left           /* fixme, length-encode */
,  const char *right          /* fixme, length-encode */
,  int indel_allowed
,  struct sw_param* swp
)  ;


int sw_fill_nw
(  struct sw_alninfo* ai
,  SWNUM* data
,  unsigned data_size
,  const char *left           /* fixme, length-encode */
,  const char *right          /* fixme, length-encode */
,  int indel_allowed
,  struct sw_param* swp
)  ;


void sw_trace
(  struct sw_alninfo* ai
,  struct sw_param* swp
,  unsigned ij
)  ;


void sw_trace_trim
(  struct sw_alninfo* ai
,  struct sw_param* swp
,  unsigned ij
,  unsigned trim_flags
)  ;


void sw_trace_nw
(  struct sw_alninfo* ai
,  struct sw_param* swp
,  unsigned ij
)  ;

void sw_printaln
(  struct sw_alninfo* ai
,  FILE* out
)  ;

#if 0
void sw_pp
(  struct sw_alninfo* ai
,  int accept
,  int score
,  FILE* fp
,  int zippit
,  int recno
)  ;
#endif

void sw_pp2
(  const char* annot_left
,  const char* annot_right
,  int cell       /* use -1 to get default cell for maximum score */
,  struct sw_alninfo* ai
,  void* fp
,  int zippit
,  int recno_ref
,  int recno_q
,  int(*cb)(const struct sw_alninfo* ai, char* buf, unsigned bufsize)
)  ;

   /* only print excised part */
void sw_pp_excise
(  const char* annot_left
,  const char* annot_right
,  int cell
,  struct sw_alninfo* ai
,  void* fp
,  int zippit
,  int recno
,  int(*cb)(const struct sw_alninfo* ai, char* buf, unsigned bufsize)
)  ;


void sw_lp                       /* line print */
(  const char* annot_left
,  const char* annot_right
,  int cell
,  struct sw_alninfo* ai
,  void* fp
,  int zippit
,  int recno_ref
,  int recno_q
,  int(*cb)(const struct sw_alninfo* ai, char* buf, unsigned bufsize)
)  ;


void sw_format                   /* format */
(  const char* annot_left
,  const char* annot_right
,  int cell
,  struct sw_alninfo* ai
,  void* fp
,  int zippit
,  int recno
,  int(*cb)(const struct sw_alninfo* ai, char* buf, unsigned bufsize)
,  const char* format
)  ;


void sw_pp_sparse
(  struct sw_alninfo* ai
,  void* fp
,  int zippit
)  ;


void sw_dump
(  struct sw_alninfo* ai
)  ;


               /* not yet used in sw.c, but this is probably the best place for future developments */
struct match_requirement
{  unsigned mr_minlen
;  unsigned mr_maxedit
;  unsigned mr_maxgap
;  int mr_offset
;
}  ;


int sw_alignment_edit_ok
(  const unsigned char* aln
,  unsigned aln_length
,  unsigned n_stretch
,  unsigned n_edit_max
,  unsigned n_gap_max
)  ;


#endif


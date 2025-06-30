
/*
 * STALLED DEVELOPMENT.
 * matecare4 was a stage in the separation of tally into two different
 * programs, namely (1) core tally functionality (tallying) (2) re-matching
 * read pairs (this program).
*/

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


#define SEP_PAIR_WIDTH 1

#ifndef TRACING_ON
#define TRACING_ON 0
#endif

#include "version.h"
#include "slib.h"
#include "trint.h"
#include "dna.h"
#include "sw.h"

#include <math.h>
#include <stdio.h>
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


unsigned g_verbose = 0;
unsigned g_ntest = 100000;

unsigned n_bogon = 0;
long g_nt_left_bogon_diff = 0;    /* number of nt lost or gained due to bogons */
size_t g_nt_left_in     = 0;
size_t g_nt_right_in    = 0;
size_t g_nt_left_out    = 0;
size_t g_nt_right_out   = 0;


struct sw_param swp_circ;
struct match_requirement mr_circ = { 16, 2, 1, 0 };

const char* me = "matecare4";



char* parse_circ_spec
(  const char* arg
)
   {  int n_converted = 0
   ;  unsigned min_length = 0, max_edit = 0, max_gap = 0
   ;  char* seq = stringle("%s", arg)
   ;  if
      (   4 != sscanf(arg, "%[ACGT]/%u/%u/%u%n", seq, &max_edit, &max_gap, &min_length, &n_converted)
      &&  3 != sscanf(arg, "%[ACGT]/%u/%u%n", seq, &max_edit, &max_gap, &n_converted)
      &&  2 != sscanf(arg, "%[ACGT]/%u%n", seq, &max_edit, &n_converted)
      &&  1 != sscanf(arg, "%[ACGT]%n", seq, &n_converted)
      )
      return NULL

   ;  swp_circ.cost_gapleft      =  3
   ;  swp_circ.cost_gapright     =  3
   ;  swp_circ.cost_subst        =  1
   ;  swp_circ.gain_match        =  4
   ;  swp_circ.left_skip         =  0
   ;  swp_circ.right_limit       =  0
   ;  swp_circ.flags             =  0

   ;  mr_circ.mr_minlen          =  min_length ? min_length : strlen(seq)
   ;  mr_circ.mr_maxedit         =  max_edit
   ;  mr_circ.mr_maxgap          =  max_gap
   ;  mr_circ.mr_offset          =  0

   ;  return seq
;  }




#define MATECARE_ERROR KRAK_READ_ERROR
#define MATECARE_NOMEM KRAK_READ_NOMEM
#define MATECARE_DONE  KRAK_READ_DONE
#define MATECARE_OK    KRAK_READ_OK


void writerecord
(  struct record_tally* rec
,  const char* format
,  void* fpo
,  int zippit
)
   {  char output[KRAK_MAXFIELDSIZE+1]
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
                           n = snprintf(output, KRAK_MAXFIELDSIZE, "%u", trint); }
                           break
            ;  case 'q':
                         { int a = (unsigned char) f[1], i; if (!a) a = '~';
                           for (i=0;i<rec->seq_n && i < KRAK_MAXFIELDSIZE; i++)
                              output[i] = a;
                           output[i] = '\0';
                           n = i;
                           f++;
                           break;
                         }
            ;  case 'I': p = rec->id ;   n = rec->id_n ;  break
            ;  case 'J': n = snprintf(output, KRAK_MAXFIELDSIZE, "%u", rec->ID) ;  break
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


#define MATRIX_SIZE (1 << 14)


int sw_match_any
(  struct sw_alninfo* ai
,  struct match_requirement* mr
)
   {  return
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


int decirc
(  struct record_tally* rec
,  const char* circ
,  const char* circ_rc
)
   {  struct sw_alninfo ai[4] = { { 0 }, { 0 }, { 0 }, { 0 } }
   ;  int indel_allowed = mr_circ.mr_maxgap > 0
   ;  int hit[4] = { 0 }
   ;  int idx_0f = 0, idx_0r = 1, idx_1f = 2, idx_1r = 3
   ;  int i
   ;  int n_hits = 0
   ;  int dbg = 0

   ;  for (i=0;i<2;i++)
      {  struct record_tally* rp = rec+i
      ;  char* s = rp->seq
      ;  int idx_f = 2*i + 0
      ;  int idx_r = 2*i + 1
      ;  SWNUM data[MATRIX_SIZE] = { 0 }

      ;  if (sw_fill2(ai+idx_f, data, MATRIX_SIZE, s, circ, indel_allowed, &swp_circ))
         X_ERR_JUMP(DONE, "sw error in circularisation tag search I")
      ;  sw_trace3(ai+idx_f, &swp_circ, ai[idx_f].max_ij, SW_TRIM)

      ;  if (sw_fill2(ai+idx_r, data, MATRIX_SIZE, s, circ_rc, indel_allowed, &swp_circ))
         X_ERR_JUMP(DONE, "sw error in circularisation tag search II")
      ;  sw_trace3(ai+idx_r, &swp_circ, ai[idx_r].max_ij, SW_TRIM)
   ;  }

      hit[idx_0f] = sw_match_any(ai+idx_0f, &mr_circ)
   ;  hit[idx_0r] = sw_match_any(ai+idx_0r, &mr_circ)
   ;  hit[idx_1f] = sw_match_any(ai+idx_1f, &mr_circ)
   ;  hit[idx_1r] = sw_match_any(ai+idx_1r, &mr_circ)
   ;

      for (i=0; i<4; i++)
      {  struct sw_alninfo* ap = ai+i
      ;  struct record_tally* rp = rec + i/2
      ;  int found_alignment = hit[i]
      ;  int j
      ;  n_hits += hit[i] != 0
      ;  if (found_alignment)
         {  ap->tag_a = (ap->lft_start-1) - SW_RGT_OHL(ap)        /* tag_a and tag_b are in C 0-offset string space with tag_b inclusive */
         ;  ap->tag_b = (ap->lft_end-1) + SW_RGT_OHR(ap)
         ;  if (ap->tag_a < 0)
            ap->tag_a = 0
         ;  if (ap->tag_b >= ap->ni )
            ap->tag_b = ap->ni-1

         ;  for (j=ap->tag_a; j<=ap->tag_b; j++)
            rp->seq[j] = i % 2 ? 'R' : 'F'

;if(dbg)fprintf(stderr, "ummmm %d %d\n", (int) ap->lft_start, (int) ap->lft_end)
;if(dbg)fprintf(stderr, "%d (%d) %d %d %d %d %.*s\n"
         ,  found_alignment
         ,  i
         ,  ap->tag_a
         ,  ap->tag_b
         ,  ap->rgt_start
         ,  ap->rgt_end
         ,  ap->tag_b - ap->tag_a +1
         ,  rp->seq+ap->tag_a
   )
         ;  if (ap->tag_b - ap->tag_a + 1 < 5)
            {  fprintf(stderr, "ummmm %d %d ni=%d ohl=%d ohr=%d\n", (int) ap->lft_start, (int) ap->lft_end, (int) ap->ni, (int) SW_RGT_OHL(ap), (int) SW_RGT_OHR(ap))
            ;  exit(23)
         ;  }
      ;  }
      }

if (dbg && n_hits)
fprintf(stderr, "\n");
      DONE
   :  return 0
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


#define MATECARE_VERBOSE_CMP  1 << 0
#define MATECARE_FORCE_PAIRS  1 << 2


#define STATUS_CLINE_PROCESSING 22
#define STATUS_DATA_PROCESSING  99


int matecare_main
(  int argc
,  const char* argv[]
)
   {  int status = STATUS_CLINE_PROCESSING
   ;  int called_from_R = 0
   ;  int changed_name = 0
   ;  int lower = 0, upper = 0, trintco = 0
   ;  unsigned n_pass = 0

   ;  const char* gzopenmode = ZWMODE

   ;  unsigned int matecare_modes = 0

   ;  ZFILE fpo = NULL, input = NULL, input2 = NULL, fpo2 = NULL
   ;  FILE* fpsumstat = NULL
   ;  const char* g_fnin = "-", *g_fnin2 = NULL
   ;  const char* g_fnout = "out.mate1.gz", *g_fnout2 = "out.mate2.gz"
   ;  const char* g_fnsumstat = NULL

   ;  char* g_circ_seq = NULL, *g_circ_rc = NULL

   ;  const char* fastqformatin     =  "@%I%#%R%n+%#%Q%n"
   ;  const char* fastqxformatin    =  "@%I%brecno=%J%#%R%n+%#%Q%n"

   ;  const char* fastaformatin     =  ">%I%#%R%n"
   ;  const char* fastaxformatin     =  ">%I%brecno=%J%#%R%n"

   ;  const char* fastqformatout    =  "@" "%I" "%n%R%n+%n%Q%n"
   ;  const char* fastaformatout    =  ">" "%I" "%n%R%n"

   ;  const char* g_format_in       =  fastqformatin
   ;  const char* g_format_out      =  fastqformatout

   ;  char* my_fnout   = NULL
   ;  unsigned limit   = 0
   ;  unsigned limit2  = 0
   ;  int zippit = 1
   ;  int do_output = 1
   ;  unsigned truncated = 0
   ;  unsigned n_pair_error = 0, n_pair_id = 0

   ;  struct record_tally whoopsie[2]  =  { { { 0 }, { 0 }, { 0 }, { 0 }, 0 }
                                          , { { 0 }, { 0 }, { 0 }, { 0 }, 0 }
                                          }
   ;  struct record_tally* r1 = whoopsie+0
   ;  struct record_tally* r2 = whoopsie+1

   ;  struct file_buffer fb1 = { { 0 }, 0 }
   ;  struct file_buffer fb2 = { { 0 }, 0 }

   ;  unsigned id_prev = 0
   ;  unsigned id_prev2 = 0

#ifdef WINDOWS_BUILD
   ;  set_argh("reaperlog.txt", 1)
#endif

   ;  kraken_readline(NULL, NULL, 0, NULL, &truncated)            /* resets static buffers */
   ;  kraken_hookline(NULL, &fb1, NULL, 0, NULL, &truncated)      /* resets buffers */
   ;  kraken_hookline(NULL, &fb2, NULL, 0, NULL, &truncated)      /* resets buffers */

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
      optarg("-j")   g_fnin2  = thearg(); endarg()

      optarg("-circ")   g_circ_seq = parse_circ_spec(thearg());
         if (!g_circ_seq)
         X_ERR_JUMP(DONE, "error parsing -circ argument as SEQ[/edit[/gap]]");
         if (!(g_circ_rc = malloc(strlen(g_circ_seq)+1)))
         X_ERR_JUMP(DONE, "error mallocing circ reversec complement");
         revcompl(g_circ_seq, strlen(g_circ_seq), g_circ_rc);
         fprintf(stderr, "seq [%s] [%s] %d %d\n", g_circ_seq, g_circ_rc, mr_circ.mr_maxedit, mr_circ.mr_maxgap);
      endarg()

      optarg("-do")  limit    = atoi(thearg()); endarg()
      optarg("-gzwrite")  gzopenmode = thearg(); endarg()
      optarg("-sumstat") g_fnsumstat = thearg(); endarg()

      optarg("-format") g_format_out = thearg(); endarg()
      optarg("-record-format") g_format_in = thearg(); endarg()

      uniarg("--fastqx-in") g_format_in  = fastqxformatin; endarg()
      uniarg("--fastax-in") g_format_in  = fastaxformatin; endarg()
      uniarg("--fasta-in")  g_format_in  = fastaformatin;  endarg()
      uniarg("--fasta-out") g_format_out = fastaformatout; endarg()

      optarg("-l")   lower = atoi(thearg()); endarg()
      optarg("-u")   upper = atoi(thearg()); endarg()
      optarg("-tri")
         trintco = atoi(thearg());
         if (trintco < 0 || trintco > 100) {
            argh(me, "-trint option should take argument in [0..100]");
            trintco = 100;
         }
      endarg()

      optarg("-give") limit2 = atoi(thearg()); endarg()

      uniarg("--nozip") zippit = 0; endarg()
      uniarg("--noput") do_output = 0; endarg()

      uniarg("--pair-by-offset") matecare_modes |= MATECARE_FORCE_PAIRS; endarg()
      optarg("-v")
         if (strstr(thearg(), "cmp")) matecare_modes |= MATECARE_VERBOSE_CMP;
         if (strstr(thearg(), "work")) g_verbose = 1;
      endarg()
      uniarg("--R") called_from_R = 1;    endarg()

      uniarg("-h")
puts("-i <fname>        first mate input stream (gzipped file allowed)");
puts("-j <fname>        second mate input stream (gzipped file allowed)");
puts("-o <fname>        first mate (gzipped!) output stream (default out.mate1.gz)");
puts("-p <fname>        second mate (gzipped!) output stream (default out.mate2.gz)");
fprintf(stdout, "--fasta-in        expect FASTA format (same as -record-format '%s')\n", fastaformatin);
fprintf(stdout, "--fasta-out       write FASTA format (same as -format '%s')\n", fastaformatout);
fprintf(stdout, "--fastax-in       expect FASTA format augmented with record offset (same as -record-format '%s')\n", fastaxformatin);
fprintf(stdout, "--fastqx-in       expect FASTQ format augmented with record offset (same as -format '%s')\n", fastqxformatin);
puts("");

puts("-l <int>          require read length >= <int>");
puts("-u <int>          require read length <= <int>");
puts("-sumstat  <fname> output file with counts of discarded categories");
puts("");

puts("PAIRING OF MATES");
puts("   There are two modes. In the first the input is assumed to be paired");
puts("   by offset; both files contain the same number of records and there is");
puts("   a record-by-record correspondence. In this case, use --pair-by-offset.");
puts("");

puts("   In the second mode, the original pairing is no longer supposed to be present.");
puts("   This can happen when the sequence files were independently processed");
puts("   and reads were filtered out.");
puts("   This is the case for example when reaper has been used to discard reads");
puts("   using criteria such as read quality and sequence complexity.");
puts("   In this case it is necessary that the original record offset information");
puts("   is preserved in the filtered read files that are the input to matecare4.");
puts("   Reaper will include a tag recno=<NUMBER> in its output when supplied");
puts("   with either of the options --fastax-out or --fastqx-out. matecare4");
puts("   will parse these formats if supplied with one of the corresponding options");
puts("   --fastax-in or --fastqx-in (as can be seen in the description of these options above).");
puts("   Alternatively, custom formats containing such a record offset number");
puts("   can be specified with the -record-format option and by including the %J");
puts("   identifier in the record format (see below).");
puts("");

puts("-record-format    specify input format");
puts("     The same syntax as documented under reaper --record-format,");
puts("     Additionally %J is accepted and assumes a numerical ID that");
puts("     will be strictly increasing. This ID will be used to match reads.");
puts("--pair-by-offset  assume the -i and -j input files match record-by-record");
puts("     With this option the %J directive is not needed");
puts("-v <work|cmp>     turn on verbosity settings");
puts("      cmp         with cmp paired end identifier mismatches will be reported");
puts("-format <format>  output format specification, syntax below");
puts("     %R  read");
puts("     %Q  quality");
puts("     %L  length");
puts("     %C  number of occurrences");
puts("     %T  trinucleotide score");
puts("     %q<CHAR> output mono-quality score <CHAR>");
puts("     %I  copy input read identifier");
puts("     %J  copy input ordinal identifier (as output by reaper or other program)");
puts("     %t  tab");
puts("     %s  tab");
puts("     %n  newline");
puts("     %%  percentage character");
exit(0);
      endarg()

      uniarg("--version")
         fprintf(stdout, "matecare version: %s\n", matecare_tag);
         exit(0);
      endarg()

      failarg()

      arg_done()

   ;  if (called_from_R)
      argh(me, "R is calling")

   ;  if
      (  g_fnin2
      && !(matecare_modes & MATECARE_FORCE_PAIRS)
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

      if (g_fnsumstat && !(fpsumstat = myfopen(g_fnsumstat, "w", 0)))
      X_ERR_JUMP(DONE, "bailing out")


   ;  status = STATUS_DATA_PROCESSING
   ;
                                                      /* funcify? */
      {  int lus   = lower || upper
      ;  unsigned active = 3                          /* bit field */
      ;  unsigned force  = matecare_modes & MATECARE_FORCE_PAIRS

      ;  while (1)
         {  unsigned readme = 1        /* 1 == skip input1, 2 == skip input2, 3 == read both */

         ;  readme = (r1->ID == r2->ID || force ? 3 : r1->ID < r2->ID ? 1 : 2) & active

         ;  if (limit && (r1->roffset > limit || r2->roffset >= limit))
            break

         ;  switch(readme)
            {  case 1: r1->discarded_unmatched += r1->count; break
            ;  case 2: r2->discarded_unmatched += r2->count; break
            /* case 3: read from both streams, meaning that r1 and r2 matched. */
         ;  }

            if (readme & 1)
            {  unsigned thestat = read_record_tally(input, &fb1, g_format_in, r1)
            ;  if (thestat == KRAK_READ_DONE)
               {  active ^= 1
               ;  r1->ID = UINT_MAX       /* so that the other stream will be instructed to read */
            ;  }
               else if (thestat == KRAK_READ_NOMEM || thestat == KRAK_READ_ERROR)
               goto DONE

            ;  if (r1->q_n != r1->seq_n)
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
            {  unsigned thestat = read_record_tally(input2, &fb2, g_format_in, r2)
            ;  if (thestat == KRAK_READ_DONE)
               {  active ^= 2
               ;  r2->ID = UINT_MAX       /* so that the other stream will be instructed to read */
            ;  }
               else if (thestat == KRAK_READ_NOMEM || thestat == KRAK_READ_ERROR)
               goto DONE

            ;  if (r2->q_n != r2->seq_n)
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

         ;  if (r1->ID == r2->ID)
            {  r1->n_paired++
            ;  if (r1->id_n || r2->id_n)
               {  if (strcmp(r1->id, r2->id))
                  {  n_pair_error++
                  ;  if (matecare_modes & MATECARE_VERBOSE_CMP)
                     fprintf(stderr, "cmp %s %s\n", r1->id, r2->id)
               ;  }
                  n_pair_id++
            ;  }
            }
            else
            continue

         ;  {  struct record_tally *rp

            ;  if (g_circ_seq)
               decirc(whoopsie, g_circ_seq, g_circ_rc)

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
                     (  (lower && rp->seq_n < lower)
                     || (upper && rp->seq_n > upper)
                     )
                     {  rp->skip++
                     ;  rp->discarded_length += rp->count
                  ;  }
                  }

                                          /* below depends on bound check above */
                  if (!rp->skip && trintco)
                  {  int trint = trintscore(rp->seq, rp->seq_n)
                  ;  if (trint >= trintco)
                     {  rp->discarded_trint += rp->count
                     ;  rp->skip++
                  ;  }
                  }
               }

               if (r1->skip && !r2->skip)
               r2->discarded_byother += r2->count
            ;  else if (!r1->skip && r2->skip)
               r1->discarded_byother += r1->count

            ;  if (r1->skip + r2->skip)
               continue
         ;  }

            if (do_output)
            {  writerecord(r1, g_format_out, fpo, zippit)
            ;  g_nt_left_out += r1->seq_n * r1->count
            ;  writerecord(r2, g_format_out, fpo2, zippit)
            ;  g_nt_right_out += r2->seq_n * r2->count
         ;  }

            n_pass += r1->count

         ;  g_nt_left_in  += r1->seq_n * r1->count
         ;  g_nt_right_in += r2->seq_n * r2->count

         ;  if (limit2 && n_pass >= limit2)
            break
      ;  }

         if (input)  myfzclose(input, 1)
      ;  if (input2) myfzclose(input2, 1)
   ;  }


      {  status = 0
      ;  if (n_pair_error)
         argh("_____", "WARNING %u records had non-matching identifiers", (unsigned) n_pair_error)
      ;  argh(me, "stringID error/compared/paired %d/%d/%d", (int) n_pair_error, (int) n_pair_id, (int) r1->n_paired)
      ;  argh(me, "paired %d reads (out of %d for -i and %d for -j argument)", (int) r1->n_paired, (int) r1->roffset, (int) r2->roffset)
      ;  goto DONE
   ;  }

      status = 0
   ;  DONE
      :

      if (do_output)
      {  if (fpo)  myfzclose(fpo, zippit)
      ;  if (fpo2) myfzclose(fpo2, zippit)
   ;  }

      if (status != STATUS_CLINE_PROCESSING && do_output)
      {  FILE* thef = fpsumstat ? fpsumstat : stderr
      ;  const char* warn_left = "", *warn_right = ""
      ;  fprintf(thef, "discarded_left_unmatched=%u\n", r1->discarded_unmatched)
      ;  fprintf(thef, "discarded_left_alien=%u\n", r1->discarded_alien)
      ;  fprintf(thef, "discarded_left_length=%u\n", r1->discarded_length)
      ;  fprintf(thef, "discarded_left_trint=%u\n", r1->discarded_trint)

      ;  fprintf(thef, "discarded_left_byright=%u\n", r1->discarded_byother)
      ;  fprintf(thef, "error_left_quality=%u\n", r1->error_quality)

      ;  if (g_nt_left_in + g_nt_left_bogon_diff != g_nt_left_out)
         warn_left = "______"

      ;  if (g_nt_right_in + - g_nt_left_bogon_diff != g_nt_right_out)
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

         fprintf(thef, "error_right_quality=%u\n", r2->error_quality)
      ;  fprintf(thef, "passed_total=%u\n", n_pass)
      ;  fprintf(thef, "num_records_left=%u\n", r1->roffset)
      ;  fprintf(thef, "num_records_right=%u\n", r2->roffset)

      ;  {  unsigned long leftsum
         =     r1->discarded_unmatched
            +  r1->discarded_alien
            +  r1->discarded_length
            +  r1->discarded_trint
            +  r1->discarded_byother
            +  n_pass

         ;  unsigned long rightsum
         =     r2->discarded_unmatched
            +  r2->discarded_alien
            +  r2->discarded_length
            +  r2->discarded_trint
            +  r2->discarded_byother
            +  n_pass

         ;  if (leftsum != r1->roffset_counted)
            argh("_____", "count discrepancy in left input file (sumstat sum %lu counted offset %u)", leftsum, r1->roffset_counted)

         ;  if (rightsum !=  r2->roffset_counted)
            argh("_____", "count discrepancy in right input file (sumstat sum %lu counted offset %u)", rightsum, r2->roffset_counted)
      ;  }
      }

      if (changed_name)
      argh(me, "\n\nBeware! output name changed to have .gz suffix, now %s", my_fnout)

   ;  if (fpsumstat) myfzclose(fpsumstat, 0)
   ;  if (my_fnout) free(my_fnout)

   ;  return status
;  }


#ifndef BUILD_R_BINDINGS

int main
(  int argc
,  const char* argv[]
)
   {  return matecare_main(argc, argv)
;  }

#endif



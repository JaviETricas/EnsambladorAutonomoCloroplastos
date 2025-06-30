

/* Author: Stijn van Dongen

 * features:
    - close to doing the very minimal amount of work required
    - no limits on matrix size; limited by available memory
    - transpose tabular data without using in-memory transpose
    - read data as string, write transpose by following array of pointers
    - record and field separator can be changed.
    - can read and write gzipped files.
    - can handle headers lacking dummy lead field.
    - in pad mode (--pad-ragged) it will padd the transposed table with empty fields
      to make it rectangular, should the input have rows with different amount of cells.
    - In normal operation, a ragged array will be considered an error and lead to termination.
    - check mode

 * caveats/TODO:
    - A file such as this:
\thdr2\thdr3\n
\trowname\tcell1\tcell2\tcell3\n
      seems to have two columns + leading tab in the header.  transpose will
      see two tabs, consider this three columns (the first of which is empty), view
      this as a three column header which can be
      padded by adding an empty column, resulting in two empty columns.

 * TODO
    - enable format-preservation without adding leading tab if header line is short.
    - allow mapping of empty/NaN/NA/-/ fields.
    - encapsulate input reading
    ? Open output file only if tests are passed.
    ? write row, col labels
    / test with empty lines (rec_offset stress-test)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>

#include "slib.h"


unsigned char g_tab = '\t';
unsigned char g_nl  = '\n';
int removed_last_newline = 0;



#define DEBUG_MATRIX     (1 << 0)
#define DEBUG_ROW        (1 << 1)
#define DEBUG_CELL       (1 << 2)

unsigned g_debug = 0;
int g_check = 0, g_pad = 0, g_normalise = 1, g_verbose = 1;



struct offset
{  unsigned long this            /* when writing a cell, this is where we start to write */
;  unsigned long next            /* and this is where it ends. Next loop this is set to next */
;  unsigned long n_sep
;
}  ;


      /* Returns 0 if we need to add a fake tab, 1 otherwise
      */
int field_check_ok
(  unsigned char* data
,  unsigned long N
,  struct offset* offsets
,  unsigned long n_offsets
)
   {  int delta = 1
   ;  unsigned long i = 0
   ;  unsigned long n_nok = 0

   ;  for (i=0;i<n_offsets-1;i++)
      {  unsigned char* o = data+offsets[i].this, *p = o, *z = o + offsets[i].next
      ;  while (p < z && (p = memchr(p, g_tab, offsets[i].next - offsets[i].this - (p-o))))
            offsets[i].n_sep++
         ,  p++
;if (g_debug & DEBUG_ROW)
 fprintf(stderr, "line %lu this %lu next %lu nsep %lu\n", i, offsets[i].this, offsets[i].next, offsets[i].n_sep)
;     }

                           /* In check/strict mode we compare with the *second* record.
                            * The first record may have things going on with the first column
                            * -- it may be meaningless and hence missing. Additionally,
                            * the separator may be mising as well.
                           */
      if (n_offsets > 2)
      {  unsigned long reference = offsets[1].n_sep, n_header = offsets[0].n_sep-1

      ;  if (g_debug & DEBUG_MATRIX)
         fprintf(stderr, "reference count %lu, first line %lu\n", reference, n_header)

      ;  if (n_header != reference && (n_header+1) != reference)
         {  if (g_verbose)
            arrr("Header line and second line do not match up (%lu versus %lu)", n_header, reference)
         ;  if (g_check || !g_pad)
            exit(1)
      ;  }

         for (i=2;i<n_offsets-1;i++)
         if (offsets[i].n_sep != reference)
         {  if (g_verbose)
            fprintf(stderr, " %lu/%lu", (i+1), (offsets[i].n_sep+1))
         ;  n_nok++
      ;  }

         if (n_nok)
         {  if (g_verbose)
            fprintf(stderr, " (line no/field count)\n")
         ;  if (g_verbose)
            arrr
            ( "Count discrepancies found (%lu), used second line as reference with %lu fields, %lu lines total"
            ,  n_nok
            ,  reference+1
            , (n_offsets-1)
            )
         ;  if(g_check || !g_pad)
            exit(1)
      ;  }

         if (g_verbose && !n_nok)
         arrr("Have %lu x %lu matrix", (n_offsets-1), (offsets[1].n_sep+1))
      ;  else if (g_verbose)
         arrr("Have ragged matrix with %lu rows", (n_offsets-1))

      ;  if (g_check)
         exit(0)
   ;  }
                           /* NOTE for the first line an artificial separator was added (add_fake_tab)
                            * Hence we do not have offsets[0].n_sep+1 as above
                           */
      else if (n_offsets == 2 && g_verbose)
      arrr("Have %lu x %lu matrix", (n_offsets-1), (offsets[0].n_sep))
   ;  else
      arrr("Have 0 x 0 matrix")     /* This is unreachable code; an empty file will be considered
                                     * to have a single cell with an empty string in it.
                                    */

   ;  if (g_check)
      exit(0)
   ;

                     /* (add_fake_tab) Below comparison looks crazy. That's because we also counted the
                      * fake tab that exists as a precaution. It was easiest to do that the
                      * way it is done, and accept this bizarre comparison -- bizarre and crazy
                      * because if the identity is true we have in fact an off-by-one situation,
                      * indicating a header line with one fewer element than the second line.
                     */
      if (i >= 2 && offsets[0].n_sep == offsets[1].n_sep && n_nok == 0 && g_normalise)
      {  if (g_verbose)
         arrr
         (  "first and second line have %lu and %lu separators; I'll add a leading separator"
         ,  offsets[0].n_sep-1
         ,  offsets[1].n_sep
         )
      ;  delta = 0
   ;  }

      return delta
;  }


void do_transpose
(  void* fpo
,  unsigned char* data
,  unsigned long N
,  struct offset* offsets
,  unsigned long n_offsets
,  int zippit
)
   {  int n_active, n_rows = 0
   ;  if (field_check_ok(data, N, offsets, n_offsets))
         offsets[0].this = 1       /* No need to parse the fake tab */
      ,  data[1] = g_nl            /* Not meaningful; still need to test this however. */
   ;  do
      {  unsigned long i
      ;  n_active = 0           /* this keeps track whether an unfinished input column still exists */
      ;  for (i=0; i<n_offsets-1;i++)
         {  unsigned long o = offsets[i].this
         ;  unsigned remaining = offsets[i].next - offsets[i].this
         ;  unsigned char* p = remaining ? memchr(data+o, g_tab, remaining) : NULL
         ;  if (remaining && !p)
            p = data+offsets[i].next
         ;  if (i) {
#if WE_USE_ZLIB
               if (zippit)
               gzputc(fpo, g_tab)
            ;  else
#endif
               fputc(g_tab, fpo)
         ;  }
            if (p)
            {  p[0] = '\0'
#if WE_USE_ZLIB
            ;  if (zippit)
               gzwrite(fpo, data+o+1, p-(data+o)-1)
            ;  else
#endif
               fwrite(data+o+1, 1, p-(data+o)-1,  fpo)
            ;  offsets[i].this = p-data
         ;  }
            if (offsets[i].this < offsets[i].next)
            n_active++
      ;  }
         if (n_active || g_nl == '\n') {
#if WE_USE_ZLIB
            if (zippit)
            gzputc(fpo, g_nl)
         ;  else
#endif
            putc(g_nl, fpo)
         ;  n_rows++
      ;  }
      }
      while (n_active > 0)

   ;  if (removed_last_newline) {
#if WE_USE_ZLIB
         if (zippit)
         gzputc(fpo, '\n')
      ;  else
#endif
         putc('\n', fpo)
   ;  }
   }


int main
(  int argc
,  char* argv[]
)
   {  unsigned long n_data_alloc = 1<<20
   ;  unsigned char* data = myalloc(n_data_alloc)  /* must be bigger than BUFSIZE + 1 */
   ;  void* fpo = NULL
   ;  ZFILE input
#define BUFSIZE 4096
   ;  unsigned long N = 0, n_offsets = 0, n_offsets_alloc = 1<<10
   ;  long n = 0
   ;  struct offset *offsets = myalloc(n_offsets_alloc * sizeof offsets[0])
   ;  int zippit = 1
   ;  const char* g_fnin = "-"
   ;  const char* g_fnout = "-"
   ;

/* enter macromagical option world */

      arg_switch()

      optarg("-i")  g_fnin = thearg();  endarg()
      optarg("-o")  g_fnout = thearg();  endarg()
      optarg("-t")  g_tab  = thearg()[0];  endarg()
      optarg("-n")  g_nl   = thearg()[0];  endarg()
      uniarg("--debug")  g_debug |= DEBUG_MATRIX; endarg()
      uniarg("-z")  g_debug = 1 | (g_debug << 1); endarg()
#if WE_USE_ZLIB
      uniarg("--nozip")  zippit = 0;  endarg()
#endif
      uniarg("--check")          g_check = 1; endarg()
      uniarg("--pad-ragged")     g_pad = 1; g_normalise = 0; endarg()
      uniarg("--quiet")          g_verbose = 0; endarg()
   uniarg("-h")
puts("-i <fname>        input stream (gzipped file allowed)");
puts("-o <fname>        output stream (gzipped file allowed)");
puts("-t <CHAR>         field delimiter (default TAB)");
puts("-n <CHAR>         record delimiter (default NEWLINE)");
puts("--pad-ragged      allow ragged table input, pad output");
puts("--debug           output some more info (use -z -z for record offset information)");
puts("--quiet           no messages");
puts("--check           check table consistency (rectangulosity) and exit");
#if WE_USE_ZLIB
puts("--nozip           do not gzip the output");
#endif
exit(0);
   endarg()

      failarg()
      arg_done()

/* exit macromagicalitaciousness */

   ;  if (!(input = myfopen(g_fnin, "r", 1)))
      exit(1)
   ;  if (!g_check && !(fpo = myfopen(g_fnout, "w", zippit)))
      exit(1)

   ;  offsets[0].this = 0
   ;  offsets[0].next = -1
   ;  offsets[0].n_sep = 0
   ;  n_offsets = 1

   ;  data[N++] = g_nl   /* should have documented this .. we don't seem to memchr it */
   ;  data[N++] = g_tab  /* in case we need to add one when header line is one short: add_fake_tab */
                         /* This *is* memchr'ed by the field_check_ok code; hence we get a
                          * discrepancy for the first line.
                         */

#  if WE_USE_ZLIB
   ;  while ((n = gzread(input, data+N, BUFSIZE)) > 0)
#  else
   ;  while ((n = fread(data+N, 1, BUFSIZE, input)) > 0)
#endif
      {  if (N+n+2+BUFSIZE > n_data_alloc)     /* extra space to write '\0', optionally '\n\0' */
         {  n_data_alloc *= 1.414
         ;  if (n_data_alloc < N+n+2+BUFSIZE)
            n_data_alloc = N+n+2+10*BUFSIZE
         ;  data = myrealloc(data, n_data_alloc)
      ;  }

         data[N+n] = '\0'        /* Should not be necessary actually */

      ;  {  unsigned char* o = data+N, *p = o
         ;  while (p-o < n && (p = memchr(p, g_nl, n-(p-o))))
            {  if (n_offsets+1 > n_offsets_alloc)       /* extra space in case no trailing g_nl */
               {  n_offsets_alloc *= 1.414
               ;  offsets = myrealloc(offsets, n_offsets_alloc * sizeof offsets[0])
            ;  }
               offsets[n_offsets].this = p-data
            ;  offsets[n_offsets-1].next = p-data
            ;  offsets[n_offsets-1].n_sep = 0
            ;  n_offsets++
            ;  p++
         ;  }
         }
         N += n
   ;  }

      myfzclose(input, 1)
   ;  if (g_verbose)
      arrr("Read %lu bytes", N-2)

   ;  if (g_nl != '\n' && data[N-1] == '\n')
         data[--N] = '\0'
      ,  removed_last_newline = 1

   ;  offsets[n_offsets-1].next  = N
   ;  offsets[n_offsets-1].n_sep = 0

   ;  if (N > 0 && data[N-1] != g_nl)
      {  data[N]  = g_nl
      ;  data[N+1] = '\0'
      ;  offsets[n_offsets].next = N
      ;  offsets[n_offsets].this   = N
      ;  offsets[n_offsets].n_sep  = 0
      ;  n_offsets++
      ;  N++
   ;  }

      offsets[n_offsets-1].next = -1

;if(0)fprintf(stdout, "[%s], %d offsets\n", data, (int) n_offsets)

   ;  do_transpose(fpo, data, N, offsets, n_offsets, zippit)
   ;  myfzclose(fpo, zippit)
   ;  free(data)
   ;  free(offsets)
   ;  return 0
;  }








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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>
#include <sys/stat.h>
#include <ctype.h>

#include "slib.h"


int slib_argv_i = 0, slib_noargallowed = 0;
const char* slib_argh_lead = "";
int g_krak_fieldtosmall = 0;

unsigned slib_verbose_level = 0;

FILE* fpargh = NULL;


int set_argh(const char* fname, int die)
   {  FILE* fp = myfopen(fname, "a", 0)
   ;  if (!fp && die)
      exit(11)
   ;  if (fp)
      fpargh = fp
   ;  else
      argh("argh", "cannot argh to file %s", fname)
   ;  return fp ? 0 : 1
;  }


int close_argh()
   {  if (fpargh && fpargh != stderr)
      fclose(fpargh)
   ;  return 0
;  }


void vargh
(  unsigned level
,  const char* tag
,  const char* fmt
,  ...
)
   {  va_list  args
   ;  FILE* fp = fpargh
   ;  if (slib_verbose_level < level)
      return
   ;  if (!fp)
      fp = stderr
   ;  va_start(args, fmt)
   ;  fprintf(fp, "%s[%s] ", slib_argh_lead, tag)
   ;  vfprintf(fp, fmt, args)
   ;  fputc('\n', fp)
   ;  va_end(args)
;  }


int argh
(  const char* tag
,  const char* fmt
,  ...
)
   {  va_list  args
   ;  FILE* fp = fpargh
   ;  if (!fp)
      fp = stderr
   ;  va_start(args, fmt)
   ;  fprintf(fp, "%s[%s] ", slib_argh_lead, tag)
   ;  vfprintf(fp, fmt, args)
   ;  fputc('\n', fp)
   ;  va_end(args)
   ;  return 29         /* Not sure what I was thinking */
;  }


int arrr
(  const char* fmt
,  ...
)
   {  va_list  args
   ;  FILE* fp = fpargh
   ;  if (!fp)
      fp = stderr
   ;  va_start(args, fmt)
   ;  vfprintf(fp, fmt, args)
   ;  fputc('\n', fp)
   ;  va_end(args)
   ;  return 29         /* Not sure what I was thinking */
;  }


void die
(  int e
,  const char* fmt
,  ...
)
   {  va_list  args
   ;  FILE* fp = fpargh
   ;  if (!fp)
      fp = stderr
   ;  va_start(args, fmt)
   ;  vfprintf(fp, fmt, args)
   ;  fputc('\n', fp)
   ;  va_end(args)
   ;  exit(e)
;  }


void* myrealloc(void* mem, size_t n)
   {  mem = realloc(mem, n)
   ;  if (!mem && n)
      {  arrr("cannae realloc %lu to bytes", (unsigned long) n)
      ;  exit(ENOMEM)
   ;  }
      return mem
;  }


char* stringle
(  const char* fmt
,  ...
)
   {  va_list  args
#define MAX_FILENAME_LENGTH 4095
   ;  char buf[MAX_FILENAME_LENGTH+1]
   ;  char* fn
   ;  int written

   ;  va_start(args, fmt)
   ;  written = vsnprintf(buf, MAX_FILENAME_LENGTH, fmt, args)
   ;  va_end(args)

   ;  if (written >= MAX_FILENAME_LENGTH || written < 0 || !(fn = malloc(written+1)))
         arrr("cannot construct full character string")
      ,  exit(1)
   ;  memcpy(fn, buf, written)
   ;  fn[written] = '\0'
   ;  return fn
;  }


void myfzclose(void* fpp, int zippit)
   {  if(fpp)
      {  int ret
#if WE_USE_ZLIB
      ;  gzflush(fpp, Z_FINISH)
      ;  ret = zippit ? gzclose(fpp) : fclose(fpp)
#else
      ;  ret = fclose(fpp)
#endif
   ;  if (ret)
      fprintf(stderr, "[slib] return code %d\n", (int) ret)
   ;  }
   }


#if ON_WINDOWS

   /* Not sure whether the do_read insertion below is actually
    * really needed under windows.
   */

void* myfopen (const char* fname, const char* mode, int z)
   {  FILE* fp = NULL
   ;  int do_read = strchr(mode, 'r') != NULL      /* this was at some point inserted to do windows compiles */
   ;  if (!strcmp(fname, "-"))
      {  if (strchr(mode, 'w'))
         fp = stdout
      ;  else if (strchr(mode, 'r'))
         fp = stdin
      ;  else
         arrr("fire: unrecognized open mode <%s>", mode)     /* fp now NULL */
   ;  }
      else if (!strcmp(fname, "stderr") && !strchr(mode, 'w'))
      fp = stderr
   ;  else
      {  fp = fopen(fname, do_read ? "r" : "wb")
      ;  if (!fp)
         arrr("cannae open file <%s> in mode <%s>", fname, mode)
   ;  }

#if WE_USE_ZLIB
      if (z && fp)
      {  int fn = fileno(fp)
      ;  if (fn < 0)
         arrr("bad fileno (file %s)", fname)
      ;  return fn < 0 ? NULL : gzdopen(fn, mode)
   ;  }
#endif
      return fp
;  }

#else

void* myfopen (const char* fname, const char* mode, int z)
   {  FILE* fp = NULL
   ;  /*   int do_read = strchr(mode, 'r') != NULL   */
   ;  if (!strcmp(fname, "-"))
      {  if (strchr(mode, 'w'))
         fp = stdout
      ;  else if (strchr(mode, 'r'))
         fp = stdin
      ;  else
         arrr("fire: unrecognized open mode <%s>", mode)     /* fp now NULL */
   ;  }
      else if (!strcmp(fname, "stderr") && !strchr(mode, 'w'))
      fp = stderr
   ;  else
      {  fp = fopen(fname, mode)
      ;  if (!fp)
         arrr("cannae open file <%s> in mode <%s>", fname, mode)
   ;  }

#if WE_USE_ZLIB
      if (z && fp)
      {  int fn = fileno(fp)
      ;  if (fn < 0)
         arrr("bad fileno (file %s)", fname)
      ;  return fn < 0 ? NULL : gzdopen(fn, mode)
   ;  }
#endif
      return fp
;  }

#endif


void* myalloc(size_t n)
   {  return myrealloc(NULL, n)
;  }

#define RL_OK     0
#define RL_DONE   1
#define RL_ERROR  2
#define RL_NOMEM  4


int kraken_hookline
(  ZFILE ip
,  struct file_buffer* fb
,  char* dest
,  unsigned  dest_capacity
,  unsigned* dest_written
,  unsigned* truncated
)
   {  unsigned n_copied             =  0
   ;  int have_newline              =  0
   ;  int have_cr                   =  0
   ;  int n_add                     =  0
   ;  int n_added                   =  0

   ;  if (!dest)                    /* special initialiser invocation */
      {  fb->rl_ct = 0
      ;  fb->rl_offset = 0
      ;  fb->rl_available = 0
      ;  fb->rl_flush = 0
      ;  return RL_OK
   ;  }

      fb->rl_ct++
   ;  truncated[0] = 0

                           /* it's a slight blemish that below reads will be attempted
                            * while feof(ip) or gzeof(ip) is already true.
                            * Solution would be to keep track of eof() if reads return 0 bytes,
                            * or to restructure the code.
                           */
   ;  do
      {  const unsigned char* p
      ;  if (!fb->rl_available)
         {  fb->rl_offset = 0
         ;  do
            {
#if WE_USE_ZLIB
               fb->rl_available = gzread(ip, fb->rl_buf, RL_BUFSIZE)
#else
               fb->rl_available = fread(fb->rl_buf, 1, RL_BUFSIZE, ip)
#endif
            ;  if (fb->rl_flush && (p = memchr(fb->rl_buf, '\n', fb->rl_available)))
               {  fb->rl_offset = 1 + p - fb->rl_buf
               ;  fb->rl_available -= 1 + p - fb->rl_buf
               ;  fb->rl_flush = 0
            ;  }
            }
            while (fb->rl_flush && fb->rl_available)
      ;  }

         if (!fb->rl_available)
         return RL_DONE

      ;  if ((p = memchr(fb->rl_buf + fb->rl_offset, '\n', fb->rl_available)))
         have_newline = 1
      ;  else
         p = fb->rl_buf + fb->rl_offset + fb->rl_available

      ;  if (p > fb->rl_buf && p[-1] == '\r')
         have_cr = 1

      ;  n_added = n_add = (p - fb->rl_buf) - fb->rl_offset

      ;  if (n_copied + n_added > dest_capacity)
         {  truncated[0] = 1
         ;  n_added = dest_capacity - n_copied        /* adjust number to copy */
      ;  }

         if (n_added > fb->rl_available)
         {  arrr("panic n_added=%d rl_available=%d", (int) n_added, (int) fb->rl_available)
         ;  return RL_ERROR
      ;  }

         memcpy(dest+n_copied, fb->rl_buf + fb->rl_offset, n_added)
#if 0
;dest[n_copied+n_added] = '\0'
;fprintf(stderr, "copy [%d] offset [%d] available [%d] line [%d] newline [%d / %d]\n"
, (int) n_added, (int) fb->rl_offset, (int) fb->rl_available, (int) fb->rl_ct, have_newline, (int) (p - fb->rl_buf))
#endif
      ;  n_copied     += n_added
      ;  fb->rl_offset    += n_added + have_newline
      ;  fb->rl_available -= n_added + have_newline
   ;  }
      while (n_added == n_add && !have_newline)

   ;  if (n_added < n_add)
      {  if (have_newline)                /* already found next offset */
            fb->rl_offset += n_add - n_added
         ,  fb->rl_available -= n_add - n_added
      ;  else                             /* this is the annoying case */
            fb->rl_flush = 1
         ,  fb->rl_offset = 0
         ,  fb->rl_available = 0
   ;  }

      if (have_cr)
      n_copied--

   ;  dest[n_copied] = '\0'               /* fixme document why this is safe */
   ;  dest_written[0] = n_copied
   ;  return RL_OK
;  }



   /* -  only allows lines of length up to dest_capacity
    * -  to reset static counters, invoke with kraken_readline(NULL, NULL, 0, 0)
   */

int kraken_readline
(  ZFILE ip
,  char* dest
,  unsigned  dest_capacity
,  unsigned* dest_written
,  unsigned* truncated
)
   {  static unsigned char rl_buf[RL_BUFSIZE]
   ;  static unsigned rl_offset     =  0
   ;  static int rl_available       =  0
   ;  static int rl_ct              =  0
   ;  static int rl_flush           =  0
   ;  static int rl_bufsize         =  RL_BUFSIZE

   ;  unsigned n_copied             =  0
   ;  int have_newline              =  0
   ;  int have_cr                   =  0
   ;  int n_add                     =  0
   ;  int n_added                   =  0

   ;  if (!ip && !dest)    /* special initialiser invocation */
      {  rl_ct = 0
      ;  rl_offset = 0
      ;  rl_available = 0
      ;  rl_flush = 0
      ;  return RL_OK
   ;  }

      rl_ct++
   ;  truncated[0] = 0

                           /* it's a slight blemish that below reads will be attempted
                            * while feof(ip) or gzeof(ip) is already true.
                            * Solution would be to keep track of eof() if reads return 0 bytes,
                            * or to restructure the code.
                           */
   ;  do
      {  const unsigned char* p
      ;  if (!rl_available)
         {  rl_offset = 0
         ;  do
            {
#if WE_USE_ZLIB
               rl_available = gzread(ip, rl_buf, rl_bufsize)
#else
               rl_available = fread(rl_buf, 1, rl_bufsize, ip)
#endif
            ;  if (rl_flush && (p = memchr(rl_buf, '\n', rl_available)))
               {  rl_offset = 1 + p - rl_buf
               ;  rl_available -= 1 + p - rl_buf
               ;  rl_flush = 0
            ;  }
            }
            while (rl_flush && rl_available)
      ;  }

         if (!rl_available)
         return RL_DONE

      ;  if ((p = memchr(rl_buf + rl_offset, '\n', rl_available)))
         have_newline = 1
      ;  else
         p = rl_buf + rl_offset + rl_available

      ;  if (p[-1] == '\r')
         have_cr = 1

      ;  n_added = n_add = (p - rl_buf) - rl_offset

      ;  if (n_copied + n_added > dest_capacity)
         {  truncated[0] = 1
         ;  n_added = dest_capacity - n_copied        /* adjust number to copy */
      ;  }

         if (n_added > rl_available)
         {  arrr("panic n_added=%d rl_available=%d", (int) n_added, (int) rl_available)
         ;  return RL_ERROR
      ;  }

         memcpy(dest+n_copied, rl_buf + rl_offset, n_added)
#if 0
;dest[n_copied+n_added] = '\0'
;fprintf(stderr, "copy [%d] offset [%d] available [%d] line [%d] newline [%d / %d]\n"
, (int) n_added, (int) rl_offset, (int) rl_available, (int) rl_ct, have_newline, (int) (p - rl_buf))
#endif
      ;  n_copied     += n_added
      ;  rl_offset    += n_added + have_newline
      ;  rl_available -= n_added + have_newline
   ;  }
      while (n_added == n_add && !have_newline)

   ;  if (n_added < n_add)
      {  if (have_newline)                /* already found next offset */
            rl_offset += n_add - n_added
         ,  rl_available -= n_add - n_added
      ;  else                             /* this is the annoying case */
            rl_flush = 1
         ,  rl_offset = 0
         ,  rl_available = 0
   ;  }

      if (have_cr)
      n_copied--

   ;  dest[n_copied] = '\0'               /* fixme document why this is safe */
;if(0)fprintf(stderr, "Yup [%s] [%.*s] offset %d cr %d hn %d\n", dest, 10, rl_buf + rl_offset, (int) rl_offset, have_cr, have_newline)
   ;  dest_written[0] = n_copied
   ;  return RL_OK
;  }



unsigned char* read_a_file
(  ZFILE ip
,  unsigned long* sizep
)
   {
#define BUFSIZE 4096
   ;  unsigned long n_data_alloc = 4096
   ;  unsigned long N = 0
   ;  long n = 0
   ;  unsigned char* d

#  if WE_USE_ZLIB
   ;  d = malloc(n_data_alloc)
   ;  while ((n = gzread(ip, d+N, BUFSIZE)) > 0)
#  else
   ;  {  struct stat mystat
      ;  if (fstat(fileno(ip), &mystat))
         fprintf(stderr,  "[read_a_file] cannae stat file\n")
      ;  else if (mystat.st_size > n_data_alloc)
         n_data_alloc = mystat.st_size + 1
   ;  }
      d = malloc(n_data_alloc)
   ;  while ((n = fread(d+N, 1, BUFSIZE, ip)) > 0)
#endif
      {  if (N+n+2+4096 > n_data_alloc)     /* extra space to write \n \0 at EOF */
         {  n_data_alloc *= 1.414
         ;  if (n_data_alloc < N+n+2+4096)
            n_data_alloc = N+n+2+10*4096
         ;  d = realloc(d, n_data_alloc)
      ;  }
         N += n
   ;  }

      if (!N || d[N-1] != '\n')
      d[N++] = '\n'

   ;  sizep[0] = N
   ;  d[N] = '\0'
   ;  return d
#undef BUFSIZE
;  }



void revcompl
(  const char* seq
,  int   n
,  char* buf
)
   {  int i
   ;  for (i=0; i<=n/2; i++)
      {  int base_lft = (unsigned char) seq[i]
      ;  int base_rgt = (unsigned char) seq[n-i-1]
      ;  buf[n-i-1] = 'A' + (base_lft > 'A' && base_lft < 'T' ? 'I' - base_lft : 'T' - base_lft)
      ;  buf[i]     = 'A' + (base_rgt > 'A' && base_rgt < 'T' ? 'I' - base_rgt : 'T' - base_rgt)
   ;  }
      buf[n] = '\0'
;  }


#if 0
unsigned char* read_file
(  ZFILE ip
,  unsigned* nbytes
)
   {  long n = 0
#define BUFSIZE 4096
   ;  unsigned long n_data_alloc = 1 << 20            /* 1M */
   ;  unsigned long N = 0
   ;  unsigned char* data = myalloc(n_data_alloc)
#  if WE_USE_ZLIB
   ;  while ((n = gzread(ip, data+N, BUFSIZE)) > 0)
#  else
   ;  while ((n = fread(data+N, 1, BUFSIZE, ip)) > 0)
#endif
      {  if (N+n+2+BUFSIZE > n_data_alloc)     /* extra space to write '\0', optionally '\n\0 */
         {  n_data_alloc *= 1.414
         ;  if (n_data_alloc < N+n+2+BUFSIZE)
            n_data_alloc = N+n+2+10*BUFSIZE
         ;  data = myrealloc(data, n_data_alloc)
      ;  }

         data[N+n] = '\0'        /* Should not be necessary actually */
      ;  N += n
   ;  }
      nbytes[0] = N
   ;  return data
#undef BUFSIZE
;  }
#endif


#define izblank(c) ((unsigned char) (c) == ' ' || (unsigned char) (c) == '\t')


int krak_cpytofield
(  char* dest
,  char* src
,  int   n
)
   {  if (n > KRAK_MAXFIELDSIZE)
         n =  KRAK_MAXFIELDSIZE
      ,  g_krak_fieldtosmall++
   ;  memcpy(dest, src, n)
   ;  dest[n] = '\0'
   ;  return n
;  }


int read_record_tally
(  ZFILE  ip
,  struct file_buffer* fb
,  const char* format
,  struct record_tally* rec
)
   {  const char* fmtp = format, *fmtz = format + strlen(format)
   ;  unsigned rlstat = 0
#ifndef READ_RECORD_LINE_LIMIT
#define READ_RECORD_LINE_LIMIT ((1<<15)-1)
#endif
   ;  char buf[READ_RECORD_LINE_LIMIT+1]
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
         {  rlstat |= kraken_hookline(ip, fb, buf, READ_RECORD_LINE_LIMIT, &n_received, &truncated)
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
               ;  rec->q_n = krak_cpytofield(rec->q, curp, bufp - curp)
               ;  break
               ;
            case 'I':
                  while (bufp < bufz && !izblank((unsigned char) bufp[0]))
                  bufp++
               ;  rec->id_n = krak_cpytofield(rec->id, curp, bufp - curp)
               ;  break
               ;
            case 'R':
                  while (bufp < bufz && !izblank((unsigned char) bufp[0]))
                     bufp[0] = toupper((unsigned char) bufp[0])
                  ,  bufp++
               ;  rec->seq_n = krak_cpytofield(rec->seq, curp, bufp - curp)
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
      ;  return KRAK_READ_ERROR
   ;  }

      if (rlstat & RL_DONE)
      return KRAK_READ_DONE

   ;  if (rlstat & RL_ERROR)
      return KRAK_READ_ERROR

   ;  if (rlstat & RL_NOMEM)           /* fixme not checked yet */
      return KRAK_READ_NOMEM

   ;  rec->roffset++
   ;  rec->roffset_counted += rec->count
   ;  return KRAK_READ_OK
;  }




int record_tally_shift_sequence
(  struct record_tally* rec
,  int delta
)
   {  int carry_quality = rec->q_n == rec->seq_n

   ;  if (delta)
      {  memmove(rec->seq, rec->seq+delta, rec->seq_n - delta)
      ;  rec->seq_n -= delta
      ;  if (carry_quality)
         {  memmove(rec->q, rec->q+delta, rec->q_n - delta)
         ;  rec->q_n -= delta
      ;  }
      }

      rec->q[rec->q_n] = '\0'
   ;  rec->seq[rec->seq_n] = '\0'
   ;  return 0
;  }




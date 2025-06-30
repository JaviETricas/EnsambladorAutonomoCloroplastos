
/*
 * (C) Copyright 2011, 2012, 2013, 2014, 2015 European Molecular Biology Laboratory.
 * Author: Stijn van Dongen <micans@gmail.com>
 * Contact: <micans@gmail.com>
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


/* This little library is not pretty, but it is not sprawling either.
 * It was written to be maximally minimalistic.
*/


#ifndef  include_slib
#define  include_slib

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>


#define KRAK_READ_ERROR 3
#define KRAK_READ_NOMEM 2
#define KRAK_READ_DONE  1
#define KRAK_READ_OK    0


#define WE_USE_ZLIB 1
#ifndef WE_USE_ZLIB
#endif

#ifndef ON_WINDOWS
#define ON_WINDOWS 0
#endif


#if ON_WINDOWS
#define COMPRESSION_LEVEL ""
#else
#define COMPRESSION_LEVEL "1"
#endif


#if WE_USE_ZLIB   /* the ZFILE abstraction carries some idiosyncratic intricacies. It's ugly, but there are reasons */
#include <zlib.h> /* Sadly this margin etc etc */
typedef gzFile ZFILE;
#define  ZWMODE "wb" COMPRESSION_LEVEL
#define  ZRMODE "rb"
#else
typedef FILE*  ZFILE;
#define  ZWMODE "w" COMPRESSION_LEVEL
#define  ZRMODE "r"
#endif

#define BIT_OFF(var, bit) (var) |= (bit), (var) ^= (bit)
#define BIT_ON(var, bit) (var) |= (bit)


extern int slib_argv_i, slib_noargallowed;

extern unsigned slib_verbose_level;

extern int g_krak_fieldtosmall;


#define thearg()  (  slib_noargallowed      \
                  ?  arrr              \
                     (  "programmer brain damage, do not use option %s"    \
                     ,  argv[slib_argv_i]   \
                     )                 \
                     ,  "argument-not-available" : argv[slib_argv_i]     \
                  )
#define arg_switch()          slib_argv_i = 0; slib_noargallowed = 0; while (++slib_argv_i < argc) { if (0) {}
#define arg_switch2()         slib_argv_i = 1; slib_noargallowed = 0; while (++slib_argv_i < argc) { if (0) {}
#define incopt(f)             if (++slib_argv_i >= argc) return arrr("option %s requires argument", f)
#define optarg(f)             else if (!strcmp(thearg(), f)) { incopt(f);
#define uniarg(f)             else if (!strcmp(thearg(), f)) { slib_noargallowed = 1;
#define uniargs(f,g)          else if (!strcmp(thearg(), f) || !strcmp(thearg(), g)) { slib_noargallowed = 1;
#define failarg()             else { return arrr("unrecognised option <%s>", thearg()); }
#define endarg()              slib_noargallowed = 0; }
#define arg_done()            }

#define X_IGNORE(...)    do {  fprintf(stderr, __VA_ARGS__); fputc('\n', stderr); return 0; } while(0)
#define X_FINISHED(...)  do {  fprintf(stderr, __VA_ARGS__); fputc('\n', stderr); return 1; } while(0)
#define X_FINISHED_VAL(val, ...)    do {  fprintf(stderr, __VA_ARGS__); fputc('\n', stderr); return val; } while(0)
#define X_CASCADE(expression) do { int ret = (expression); if (ret) return ret; } while (0)
#define X_CASCADE_CLEAN(expression_test, expression_clean) do { int ret = (expression_test); if (ret) { expression_clean; return ret; } } while (0)

#define X_ERR_JUMP(label, ...) do {  fprintf(stderr, __VA_ARGS__); fputc('\n', stderr); goto label; } while (0)
#define X_TRY_JUMP(label, expression) do {  int ret = (expression); if (ret) { goto label; } } while (0)

#define ji  (int)
#define ju  (unsigned)
#define jl  (long)
#define jlu (unsigned long)

#define THEMAX(a,b) ((a) > (b) ? (a) : (b))


   /* prints to STDERR by default, settable using fpargh */
int argh
(  const char* tag
,  const char* fmt
,  ...
)  ;

void vargh
(  unsigned verbosity_level
,  const char* tag
,  const char* fmt
,  ...
)  ;
extern FILE* fpargh;
int set_argh(const char* fname, int die);
int close_argh(void);


int arrr
(  const char* fmt
,  ...
)  ;

void die
(  int e
,  const char* fmt
,  ...
)  ;

void* myalloc
(  size_t n
)  ;

void* myrealloc
(  void* mem
,  size_t n
)  ;

char* stringle
(  const char* fmt
,  ...
)  ;

void myfzclose
(  void* fpp
,  int zippit
)  ;

   /* return zippit ? ZFILE* : *FILE */
void* myfopen
(  const char* fname
,  const char* mode
,  int zippit
)  ;

unsigned char* read_a_file
(  ZFILE fp
,  unsigned long* sizep
)  ;


#define RL_OK     0
#define RL_DONE   1
#define RL_ERROR  2
#define RL_NOMEM  4


#define READSIZE    (1 << 15)

#define RL_BUFSIZE    (1 << 15)           /* 8192 */


struct file_buffer
{  unsigned char rl_buf[RL_BUFSIZE]
;  unsigned rl_offset
;  int rl_available
;  int rl_ct
;  int rl_flush
;
}  ;


int kraken_hookline
(  ZFILE ip
,  struct file_buffer* fb
,  char* dest
,  unsigned  dest_capacity
,  unsigned* dest_written
,  unsigned* truncated
)  ;


   /* -  uses static variables. Do not use two instances simultaneously (use hookline for that).
      -  only allows lines of length up to dest_capacity
      -  is not yet able to truncate or flush lines that are longer
   */
int kraken_readline
(  ZFILE ip
,  char* dest
,  unsigned  dest_capacity
,  unsigned* dest_written
,  unsigned* truncated
)  ;


void revcompl
(  const char* seq
,  int   n
,  char* buf
)  ;



#ifndef KRAK_MAXFIELDSIZE
#define KRAK_MAXFIELDSIZE ((1 << 15)-1)
#endif


struct record_tally
{  char seq       [KRAK_MAXFIELDSIZE+1]
;  char q         [KRAK_MAXFIELDSIZE+1]
;  char discard   [KRAK_MAXFIELDSIZE+1]
;  char id        [KRAK_MAXFIELDSIZE+1]

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


int read_record_tally
(  ZFILE  ip
,  struct file_buffer* fb
,  const char* format
,  struct record_tally* rec
)  ;


int record_tally_shift_sequence
(  struct record_tally* rec
,  int delta
)  ;


#endif


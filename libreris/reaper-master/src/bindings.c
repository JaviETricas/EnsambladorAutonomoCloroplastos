
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


#include <R.h>
#include <Rinternals.h>

#include "reaper.h"
#include "tally.h"
#include "slib.h"


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

      argv[argc] = NULL

   ;  {  int j
      ;  fprintf(stderr, "Passing to %s:", callee)
      ;  for (j=0;j<argc;j++)
         fprintf(stderr, " %s", argv[j])
      ;  fputc('\n', stderr)
   ;  }
      fputc('\n', stderr)

   ;  if (i != length(list))
      {  int l = length(list)
      ;  Rprintf("length discrepancy %d %d\n", i, l)
   ;  }

      *res=themain(argc, argv)

   ;  UNPROTECT(1)
   ;  return result
;  }



SEXP reaperC
(  SEXP list
)
   {  return R_dispatchee(list, reaper_main, "reaper")
;  }



SEXP tallyC
(  SEXP list
)
   {  return R_dispatchee(list, tally_main, "tally")
;  }




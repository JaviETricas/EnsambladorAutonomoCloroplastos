
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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include <errno.h>

#include "slib.h"
#include "version.h"
#include "table.h"




int main
(  int argc
,  char* argv[]
)
   {  FILE* fpout = stdout
   ;  char* fnout = "-"
   ;  char* fnin  = "-"
   ;  ZFILE input = NULL
   ;  unsigned truncated = 0
   ;  unsigned rlstat = 0
   ;  unsigned nc = 0
   ;  struct table thetable = { NULL }

#define READLINE_MAX (1<<15)
   ;  char buf[READLINE_MAX+1]
   ;  kraken_readline(NULL, NULL, 0, NULL, NULL)            /* resets static buffers */

/* enter macromagical option world */

   ;  arg_switch()

      uniarg("--version")
fprintf(stdout, "Sword version: %s\n", sword_tag);
exit(0);
      endarg()

      uniarg("-h")
puts("");
exit(0);
      endarg()

      failarg()
      arg_done()

/* exit macromagicalitaciousness */

   ;  if (fnout && !(fpout = myfopen(fnout, "r", 1)))
      exit(1)

   ;  if (fnin && !(input = myfopen(fnin, "r", 1)))
      exit(1)

   ;  table_read(input, &thetable, 's', 0)
   ;  myfzclose(input, 1)

   ;  fclose(fpout)

   ;  return 0
;  }







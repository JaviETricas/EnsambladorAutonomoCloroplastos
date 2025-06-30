
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

/* TODO
 *    support zlib. Note: would be nice to know size of uncompressed data.
*/

#include "table.h"
#include "slib.h"


int table_check_data
(  unsigned char* d
,  unsigned long  dsize
,  unsigned long* n_col1
,  unsigned long* n_col2
,  unsigned long* n_row
)
   {  unsigned long c1 = 0, c2 = 0, r = 0, nsep = 0, i
   ;  for (i=0; i <= dsize; i++)
      {  if (d[i] == '\t')
         nsep++
      ;  else if (d[i] == '\n' || i == dsize)
         {  r++
         ;  if (r == 1)
            c1 = nsep+1
         ;  else if (r == 2)
            c2 = nsep+1
         ;  else if (c2 != nsep+1)
            {  fprintf(stderr, "column count mismatch at line %d\n", (int) r)
            ;  break
         ;  }
            if (r == 2 && c1 != c2 && c1 + 0 != c2)      /* fixme: allow c1 + 1 == c2 */
            {  fprintf(stderr, "column count mismatch at line 2\n")
            ;  break
         ;  }
            nsep = 0
      ;  }
      }
      n_col1[0] = c1
   ;  n_col2[0] = c2
   ;  n_row[0]  = r
   ;  return i == dsize + 1
;  }


void table_print_tp
(  struct table* t
)
   {  unsigned r = 0, c = 0
   ;  for (c=0;c<t->n_cols;c++)
      {  for (r=0;r<t->n_rows;r++)
         {  char* field = NULL
         ;  memcpy(&field, table_elem(t,r,c), sizeof field)
         ;  fprintf(stdout, "%s%s", r ? "\t" : "", ((char**) table_elem(t,r,c))[0])
      ;  }
         fputc('\n', stdout)
   ;  }
;  }


void table_print
(  struct table* t
)
   {  unsigned r = 0, c = 0
   ;  int rowlead = t->flags & TABLE_WITH_ROWNAMES
   ;  int collead = t->flags & TABLE_WITH_COLNAMES

   ;  if (collead)
      {  fprintf(stdout, "%s", t->dummy)
      ;  for (c=0;c<t->n_cols;c++)
         fprintf(stdout, "\t%s", t->colnames[c])
      ;  fputc('\n', stdout)
   ;  }

      for (r=0;r<t->n_rows;r++)
      {  if (rowlead)
         fputs(t->rownames[r], stdout)
      ;  for (c=0;c<t->n_cols;c++)
         {  char* field = NULL
         ;  memcpy(&field, table_elem(t,r,c), sizeof field)
         ;  fprintf(stdout, "%s%s", c || collead ? "\t" : "", ((char**) table_elem(t,r,c))[0])
      ;  }
         fputc('\n', stdout)
   ;  }
;  }


int table_fill
(  struct table* t
)
   {  unsigned long rowid = 0, colid = 0, sof = 0, i
   ;  unsigned expect_rnames = t->flags & TABLE_WITH_ROWNAMES ? 1 : 0
   ;  unsigned expect_cnames = t->flags & TABLE_WITH_COLNAMES ? 1 : 0
   ;  unsigned char* d = t->rawdata
   ;  unsigned long dsize = t->rawdata_n

   ;  for (i=0; i<= dsize; i++)
      {  const unsigned char* field = d+sof
      ;  if (i != dsize && d[i] != '\t' && d[i] != '\n')
         continue

      ;  if (memchr(d+sof, 0, i-sof))
         X_FINISHED("file contains nul bytes")

      ;  if (!rowid && expect_cnames)                 /* mq future header column count mismatch support criticalaciousness */
         {  if (expect_rnames && !colid)
            t->dummy = (char*) d+sof
         ;  else
            memcpy(t->colnames + colid - expect_rnames, &field, sizeof field)
      ;  }
         else if (!colid && expect_rnames)
         memcpy(t->rownames + rowid - expect_cnames, &field, sizeof field)

      ;  else if (t->flags & TABLE_BY_ROW_MAJOR)
         memcpy(table_elem_byrow(t, rowid - expect_cnames, colid - expect_rnames), &field, sizeof field)
      ;  else
         memcpy(table_elem_bycol(t, rowid - expect_cnames, colid - expect_rnames), &field, sizeof field)

      ;  colid++
      ;  if (d[i] == '\n')
            rowid++
         ,  colid = 0
      ;  d[i] = '\0'
      ;  sof = i+1
   ;  }
      if (colid != t->n_cols + expect_rnames || (rowid+1) != t->n_rows + expect_cnames)
      X_FINISHED("[table fill] assertion on dimensions failed")
   ;  return 0
;  }



   /* d:          t.......t........t......t......n......t......t
      t->data     [   ][   ][   ]
   */


int table_read
(  ZFILE fp
,  struct table* t
,  unsigned type
,  unsigned flags
)
   {  unsigned long n_col1 = 0, n_col2 = 0, n_row = 0
   ;  unsigned expect_rnames = flags & TABLE_WITH_ROWNAMES ? 1 : 0
   ;  unsigned expect_cnames = flags & TABLE_WITH_COLNAMES ? 1 : 0

   ;  long type_long
   ;  double type_double
   ;  unsigned char* type_string 
   ;  int status = 1

   ;  t->rawdata = read_a_file(fp, &t->rawdata_n)

   ;  if (!t->rawdata_n)
      X_ERR_JUMP(FAIL, "empty table?")

   ;  if (t->rawdata[t->rawdata_n-1] == '\n')
         t->rawdata_n--
      ,  t->rawdata[t->rawdata_n] = '\0'

   ;  if (!table_check_data(t->rawdata, t->rawdata_n, &n_col1, &n_col2, &n_row))
      X_ERR_JUMP(FAIL, "table error")

   ;  t->flags = flags
   ;  t->type  = type

   ;  if (type == 'l')
      t->elem_size = sizeof type_long
   ;  else if (type == 'd')
      t->elem_size = sizeof type_double
   ;  else if (type == 's')
      t->elem_size = sizeof type_string

   ;  if (n_row < 1 + expect_rnames)
      X_ERR_JUMP(FAIL, "need at least one row")
   ;  if (n_col2 < 1 + expect_cnames)
      X_ERR_JUMP(FAIL, "need at least one column")

   ;  t->n_rows = n_row - expect_cnames
   ;  t->n_cols = n_col2 - expect_rnames

;if(0)fprintf(stderr, "table with %d rows, %d columns elem size %d\n", (int) t->n_rows, (int) t->n_cols, (int) t->elem_size)

   ;  t->data = malloc(t->elem_size * t->n_rows * t->n_cols)
   ;  if (expect_rnames)
      t->rownames = malloc((n_row - expect_cnames) * sizeof t->rownames[0])
   ;  if (expect_cnames)
      t->colnames = malloc((n_col2 - expect_rnames) * sizeof t->colnames[0])

   ;  if (type != 's')
      X_ERR_JUMP(FAIL, "only type string supported yet")

   ;  X_TRY_JUMP(FAIL, table_fill(t))

   ;  status = 0
   ;
      FAIL
   :
;if(0)fprintf(stderr, "table with %d rows, %d columns elem size %d\n", (int) t->n_rows, (int) t->n_cols, (int) t->elem_size)
   ;  return status
;  }


int table_column_index
(  struct table* t
,  const char* id
)
   {  int hit = -1, i
   ;  for (i=0;i<t->n_cols;i++)
      {  if (!strcmp(t->colnames[i], id))
         {  if (hit >= 0)
            arrr("column %s is present more than once", id)
         ;  else
            hit = i
      ;  }
      }
   ;  return hit
;  }


void table_release
(  struct table* t
)
   {  if (t->data) free(t->data)
   ;  if (t->colnames) free(t->colnames)
   ;  if (t->rownames) free(t->rownames)
   ;  if (t->rawdata) free(t->rawdata)
   ;  t->data = NULL
   ;  t->rawdata = NULL
   ;  t->colnames = NULL
   ;  t->rownames = NULL
;  }




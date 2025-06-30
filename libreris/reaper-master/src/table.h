
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


#include "slib.h"


struct table
{  void* data
;  unsigned char* rawdata
;  unsigned long rawdata_n
;  unsigned n_rows
;  unsigned n_cols
;  char** colnames
;  char** rownames
;  char*  dummy         /* the colname for the rownames, if present    */
;  unsigned type        /* 'd'ouble, 'l'ong, long 'u'nsigned, 's'tring */
;  unsigned elem_size

#define TABLE_BY_ROW_MAJOR    1 << 0
#define TABLE_WITH_ROWNAMES   1 << 1
#define TABLE_WITH_COLNAMES   1 << 2

#define table_elem_byrow(t,i,j) (((char*) (t)->data) + (t)->elem_size * ((i) * (t)->n_cols   +   (j)))
#define table_elem_bycol(t,i,j) (((char*) (t)->data) + (t)->elem_size * ((i)   +   (t)->n_rows * (j)))
#define table_elem(t,i,j) (((t)->flags & TABLE_BY_ROW_MAJOR) ? table_elem_byrow(t,i,j) : table_elem_bycol(t,i,j))
#define table_elem_string(t,i,j) (((char**) table_elem(t,i,j))[0])

;  unsigned flags
;
}  ;



int table_check_data
(  unsigned char* d
,  unsigned long  dsize
,  unsigned long* n_col1
,  unsigned long* n_col2
,  unsigned long* n_row
)  ;


void table_print_tp
(  struct table* t
)  ;


void table_print
(  struct table* t
)  ;


int table_fill
(  struct table* t
)  ;


int table_read
(  ZFILE fp
,  struct table* t
,  unsigned type
,  unsigned flags
)  ;


int table_column_index
(  struct table* t
,  const char* id
)  ;


void table_release
(  struct table* t
)  ;



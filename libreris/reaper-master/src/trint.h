
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

#ifndef  include_trint
#define  include_trint


/* Herve Pages-type criterion (weighted frequencies) */

int trintscore
(  const char* seq
,  int seqlen
)  ;

int trintscore_tail
(  const char* seq
,  int seqlen
,  char* scratch
,  int* scorep
,  int* maxidp
)  ;


/* Dust-type criterion without the maximum-across-nested interval requirement */

int dustscore_tail
(  const char* seq
,  int seqlen
,  char* scratch
,  int* scorep
,  int* maxidp
)  ;

int dustscore
(  const char* seq
,  int seqlen
)  ;

#endif



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

#include <math.h>
#include <stdio.h>

#include "trint.h"


/* Following a suggestion by Hervé Pagès
   at https://stat.ethz.ch/pipermail/bioc-sig-sequencing/2009-February/000170.html

   tnf <- trinucleotideFrequency2(dict0)
   tnf2 <- tnf - 1L
   tnf2[tnf2 < 0] <- 0L
   scores <- rowSums(tnf2 * tnf2)

   We ignore Ns -- affected trinucleotides are ignored, but new
   ones come into existence. The result is the same as if all Ns
   were deleted from the sequence.
*/


const char* trimer_store = "AAAAACAAGAATACAACCACGACTAGAAGCAGGAGTATAATCATGATTCAACACCAGCATCCACCCCCGCCTCGACGCCGGCGTCTACTCCTGCTTGAAGACGAGGATGCAGCCGCGGCTGGAGGCGGGGGTGTAGTCGTGGTTTAATACTAGTATTCATCCTCGTCTTGATGCTGGTGTTTATTCTTGTTT";


int trintscore
(  const char* seq
,  int seqlen
)
#define BASEMAP(b)   (b == 'A' ? 0 : b == 'C' ? 1 : b == 'G' ? 2 : b == 'T' ? 3 : 4)
   {  unsigned index[64] = { 0 }
   ;  int i
   ;  double score = 0.0
   ;  unsigned tri = 0, n_used = 0
   ;  if (seqlen < 4)                        /* not strictly necessary, but hey */
      return 0

   ;  for (i=0;i<seqlen;i++)
      {  unsigned b = BASEMAP(seq[i])
      ;  if (b > 3)
         continue
      ;  n_used++
      ;  tri = ((tri << 2) | b) & 0x3f       /* x3f == d63 == b111111 */
      ;  if (n_used > 2)
         index[tri]++
   ;  }
      if (n_used > 3)                        /* that means at least two different trint */
      for (i=0;i<64;i++)
      if (index[i])
      score += ((index[i]-1) * 1.0 / (n_used-3)) * ((index[i]-1) * 1.0 / (n_used-3))
   ;  return (int) (0.5 + sqrt(score) * 100.0)
;  }


int dustscore_tail
(  const char* seq
,  int seqlen
,  char* scratch
,  int* maxp
,  int* maxidp
)
   {  int maxidx = seqlen
   ;  double thescore = 0.0, maxscore = 0.0
   ;  int i

   ;  maxp[0] = 0
   ;  if (scratch)
      scratch[seqlen] = '\0'
   ;

      {  unsigned index[64] = { 0 }
      ;  unsigned tri = 0, n_used = 0, n_n = 0
      ;  int sumhsq = 0                            /* sum of half squares */
      ;  for (i=1;i<=seqlen;i++)
         {  unsigned b = BASEMAP(seq[seqlen - i])
         ;  if (b > 3)                             /* 0-3 is ACGT */
            n_n++
         ;  else
            {  n_used++
            ;  tri = ((tri << 2) | b) & 0x3f       /* x3f == d63 == b111111 */
            ;  if (n_used > 2)
               {  index[tri]++
               ;  sumhsq += index[tri] - 1
            ;  }
            }

            {  thescore
               =     n_used > 3
                  ?  (10.0 * (sumhsq + (n_n * (n_n + 1))) / (1.0 * (n_used-3)))
                  :  0.0

            ;  if (thescore >= maxscore)
                  maxscore = thescore
               ,  maxidx = seqlen-i
            ;  if (scratch)
               scratch[seqlen-i] = thescore > 93 ? 126 : thescore + 33
         ;  }
         }
      }

      maxp[0] = (int) (0.5 + maxscore)
   ;  maxidp[0] = maxidx
   ;  return thescore
;  }


int dustscore
(  const char* seq
,  int seqlen
)
   {  int max = 0, maxid = 0
   ;  return dustscore_tail(seq, seqlen, NULL, &max, &maxid)
;  }


int trintscore_tail
(  const char* seq
,  int seqlen
,  char* scratch
,  int* maxp
,  int* maxidp
)
   {  int thescore = 0, maxscore = 0, maxidx = 0
   ;  int i

   ;  maxp[0] = 0
   ;  if (scratch)
      scratch[seqlen] = '\0'
   ;

      {  unsigned index[64] = { 0 }
      ;  unsigned tri = 0, n_used = 0, n_n = 0
      ;  int cplty = 0                             /* complexity */
      ;  for (i=1;i<=seqlen;i++)
         {  unsigned b = BASEMAP(seq[seqlen - i])
         ;  if (b > 3)                             /* 0-3 is ACGT */
            n_n++
         ;  else
            {  n_used++
            ;  tri = ((tri << 2) | b) & 0x3f       /* x3f == d63 == b111111 */
            ;  if (n_used > 2)
                                                   /* this leads to the Hervé Pagès criterion:
                                                    * #tri  desired score     increment
                                                    *    1           0        -
                                                    *    2           1        1
                                                    *    3           4        3
                                                    *    4           9        5
                                                    * et cetera
                                                   */
               {  if (index[tri] > 0)
                  cplty += index[tri] + index[tri] - 1
               ;  index[tri]++
            ;  }
            }

            {  thescore
               =     n_used > 3
                  ?  (int) (0.5 + 100.0 * sqrt((cplty + n_n * (n_n + 1)) * 1.0 / ((n_used-3) * (n_used-3))))
                  :  0

            ;  if (thescore >= maxscore)
                  maxscore = thescore
               ,  maxidx = seqlen-i

            ;  if (scratch)
               scratch[seqlen-i] = thescore > 93 ? 126 : thescore + 33
         ;  }
         }
      }

      maxp[0] = maxscore
   ;  maxidp[0] = maxidx
   ;  return thescore
;  }





#if 0
   ;  {  unsigned index[64] = { 0 }
      ;  unsigned tri = 0, n_used = 0, n_n = 0, sumhsq = 0     /* sum of half squares */
      ;  for (i=0;i<seqlen;i++)
         {  unsigned b = BASEMAP(seq[i])
         ;  if (b > 3)
            n_n++
         ;  else
            {  n_used++
            ;  tri = ((tri << 2) | b) & 0x3f       /* x3f == d63 == b111111 */
            ;  index[tri]++
            ;  sumhsq += index[tri]-1
         ;  }
            {  unsigned s = (int) (10.0 * (sumhsq + (n_n * (n_n + 1))) / (1.0 * i))
            ;  if (s > 255)
               s = 255
            ;  if (scratch)
               scratch[i] = s
         ;  }
         }
      }
#endif

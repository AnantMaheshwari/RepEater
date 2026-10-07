
/*  File: tanbed.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description: make a bed file from a FasTAN .1aln file
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 16 10:12 2025 (rd109)
 * Created: Thu Oct 16 02:48:43 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "alntools.h"

int main (int argc, char *argv[])
{
  storeCommandLine (argc, argv) ;
  --argc ; ++argv ;
  if (argc != 1) die ("Usage: tanbed <.1aln file>") ;

  OneSchema *schema = oneSchemaCreateFromText (schemaText) ;
  OneFile *of = oneFileOpenRead (*argv, schema, "aln", 1) ;
  if (!of) die ("failed to open .1aln file %s", *argv) ;
  Gdb *gdb = readGdb (of, 1, stderr) ;
  if (of->lineType != 'A') die ("unexpected line type %c", of->lineType) ;

  I64 nAlign = 0 ;
  oneStats (of, 'A', &nAlign, 0, 0) ;
  Array ab = arrayCreate (nAlign, TanLine) ;
  I64 totAlign = 0 ;
  while (of->lineType == 'A')
    { TanLine *b = arrayp(ab, arrayMax(ab), TanLine) ;
      int ctg = oneInt(of,0) ;
      if (ctg != oneInt(of,3))
	die ("target mismatch line %lld - not a TAN file?", (long long)of->line) ;
      b->seq     = ctg2seq(gdb,ctg) ;
      b->start   = ctg2pos(gdb,ctg,oneInt(of,4)) ;
      b->end     = ctg2pos(gdb,ctg,oneInt(of,2)) ;
      double len = oneInt(of,2) - oneInt(of,4) ;
      totAlign  += oneInt(of,2) - oneInt(of,4) ;
      while (oneReadLine(of) && of->lineType != 'A')
	if (of->lineType == 'D') b->score = (int)(1000*(1.0 - oneInt(of,0)/len)) ;
	else if (of->lineType == 'U') b->unit = oneInt(of,0) ;
    }
  oneFileClose (of) ;
  arraySort (ab, tanSort) ;

  int i ;
  for (i = 0 ; i < arrayMax(ab) ; ++i)
    { TanLine *b = arrp(ab, i, TanLine) ;
      printf ("%s\t", dictName(gdb->seqDict, b->seq)) ;
      printf ("%lld\t%lld\t%d\t%d\n", (long long)b->start, (long long)b->end, b->unit, b->score) ;
    }
  
  I64 totSeq = 0 ; for (i = 0 ; i < gdb->nSeq ; ++i) totSeq += gdb->seqLen[i] ;
  fprintf (stderr, "processed %lld alignments total length %lld from %s length %lld (%.1f %%)\n",
	   nAlign, totAlign, *argv, totSeq, totAlign/(0.01*totSeq)) ;
}

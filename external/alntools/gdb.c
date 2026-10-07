/*  File: gdb.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 16 10:04 2025 (rd109)
 * Created: Sun Oct 19 21:45:26 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "alntools.h"

void reportGdb (Gdb *gdb, FILE *f)
{ fprintf (f, "%d seqs %d contigs (%d gaps)", gdb->nSeq, gdb->nCtg, gdb->nGap) ;
  fprintf (f, ", totSeq %lld totCtg %lld (%.3f%%)", (long long) gdb->totSeq,
	   (long long) gdb->totCtg, gdb->totCtg/(0.01*gdb->totSeq)) ;
  if (gdb->totMask)
    fprintf (f, ", %d masks totMask %lld (%.1f%%)",
	     (int) gdb->maxMask/2, (long long) gdb->totMask, gdb->totMask/(0.01*gdb->totSeq)) ;
  fputc ('\n', f) ;
}

Gdb *readGdb (OneFile *of, int k, FILE *report) // don't report if !report
{
  int i ;
  Gdb *gdb = new0 (1, Gdb) ;

  for (i = 0 ; i < oneReferenceCount(of) ; ++i)
    if (of->reference[i].count == k) gdb->seqFileName = of->reference[i].filename ;
    else if (of->reference[i].count == 3) gdb->seqPathName = of->reference[i].filename ;
  
  oneStats (of, 'S', &gdb->maxSeq, 0, 0) ;
  oneStats (of, 'C', &gdb->maxCtg, 0, 0) ;
  oneStats (of, 'M', 0, 0, &gdb->maxMask) ;
  gdb->seqDict = dictCreate (gdb->maxSeq) ;
  gdb->seqLen = new0 (gdb->maxSeq, I64) ;
  gdb->ctgLen = new0 (gdb->maxCtg, I64) ;
  gdb->ctgSeq = new0 (gdb->maxCtg, int) ;
  gdb->ctgPos = new0 (gdb->maxCtg, I64) ;
  if (gdb->maxMask)
    { gdb->ctgMaskCount = new0 (gdb->maxCtg, int) ;
      gdb->ctgMaskStart = new0 (gdb->maxCtg, int) ;
      gdb->mask = new0 (gdb->maxMask, I64) ;
    }

  if (!strcmp (of->fileType, "gdb"))	// stand alone
    oneGoto (of, 'S', 0) ; 		// go to the start of the data
  else                              	// a skeleton inside another file type
    { if (!oneGoto (of, 'g', k)) 	// need to go to the start of the k'th embedded GDB
	die ("failed to go to GDB %d in %s", k, oneFileName (of)) ;
      oneReadLine (of) ;         	// swallow the g line
    }

  I64  end = 0, inMask = 0 ;
  bool isDone = false ;
  while (!isDone && oneReadLine (of))
    switch (of->lineType)
      {
      case 'f':
	gdb->fA = oneReal(of,0) ; gdb->fC = oneReal(of,1) ;
	gdb->fG = oneReal(of,2) ; gdb->fT = oneReal(of,3) ;
	break ;
      case 'u':
	gdb->isUpper = true ;
	break ;
      case 'S':
	dictAdd (gdb->seqDict, oneString(of), 0) ;
	if (gdb->nSeq > 0) gdb->seqLen[gdb->nSeq-1] = end ;
	end = 0 ; // reset counter
	++gdb->nSeq ;
	break ;
      case 'G':
	if (!gdb->nSeq) die ("G line before S line in GDB") ;
	end += oneInt(of,0) ;
	gdb->totSeq += oneInt(of,0) ;
	++gdb->nGap ;
	break ;
      case 'C':
	if (!gdb->nSeq) die ("C line before S line in GDB") ;
	gdb->ctgSeq[gdb->nCtg] = gdb->nSeq-1 ;
	gdb->ctgLen[gdb->nCtg] = oneInt(of,0) ;
	gdb->ctgPos[gdb->nCtg] = end ;
	end += oneInt(of,0) ;
	gdb->totSeq += oneInt(of,0) ;
	gdb->totCtg += oneInt(of,0) ;
	++gdb->nCtg ;
	break ;
      case 'M':
	if (!gdb->nCtg) die ("M line before C line in GDB") ;
	if (oneLen(of) % 2) die ("size of Mask list must be even") ;
	if (gdb->ctgMaskCount[gdb->nCtg-1]) die ("> 1 M lines per C line in GDB") ;
	gdb->ctgMaskCount[gdb->nCtg-1] = oneLen(of) ;
	gdb->ctgMaskStart[gdb->nCtg-1] = inMask ;
	memcpy (&gdb->mask[inMask], oneIntList(of), oneLen(of)*sizeof(I64)) ;
	inMask += oneLen(of) ;
	I64 *m = oneIntList(of) ;
	for (i = 0 ; i < oneLen(of) ; i += 2) gdb->totMask += m[i+1] - m[i] ;
	break ;
      default: // anything else (including another 'g') is the end of this GDB
	isDone = true ;
	break ;
      }
  // must close the final sequence
  if (gdb->nSeq > 0) gdb->seqLen[gdb->nSeq-1] = end ;

  if (report)
    { fprintf (report, "read GDB from %s : ", oneFileName(of)) ;
      reportGdb (gdb, report) ;
    }
  
  return gdb ;
}

void writeGdb (OneFile *of, Gdb *gdb, int k, FILE *report)
{
  if (gdb->seqFileName) oneAddReference (of, gdb->seqFileName, k) ;
  if (gdb->seqPathName) oneAddReference (of, gdb->seqPathName, 3) ;

  if (strcmp (of->fileType, "gdb")) oneWriteLine (of, 'g', 0, 0) ;
  if (gdb->fA + gdb->fC + gdb->fG + gdb->fT)
    { oneReal(of,0) = gdb->fA ; oneReal(of,1) = gdb->fC ;
      oneReal(of,2) = gdb->fG ; oneReal(of,3) = gdb->fT ;
      oneWriteLine (of, 'f', 0, 0) ;
    }
  if (gdb->isUpper) oneWriteLine (of, 'u', 0, 0) ;
  int i, j = 0 ;
  I64 end ;
  for (i = 0 ; i < gdb->nSeq ; ++i)
    { char *name = dictName(gdb->seqDict,i) ;
      oneWriteLine (of, 'S', strlen(name), name) ;
      end = 0 ;
      while (j < gdb->nCtg && gdb->ctgSeq[j] == i)
	{ if (gdb->ctgPos[j] > end)
	    { oneInt(of,0) = gdb->ctgPos[j] - end ;
	      oneWriteLine (of, 'G', 0, 0) ;
	      end = gdb->ctgPos[j] ;
	    }
	  oneInt(of,0) = gdb->ctgLen[j] ;
	  oneWriteLine (of, 'C', 0, 0) ;
	  end += gdb->ctgLen[j] ;
	  if (gdb->totMask && gdb->ctgMaskCount[j])
	    oneWriteLine (of, 'M', gdb->ctgMaskCount[j], &gdb->mask[gdb->ctgMaskStart[j]]) ;
	  ++j ;
	}
      if (end < gdb->seqLen[i])
	{ oneInt(of,0) = gdb->seqLen[i] - end ;
	  end += gdb->seqLen[i] - end ; // not needed, but nice to have it here
	  oneWriteLine (of, 'G', 0, 0) ;
	}
    }

  if (report)
    { fprintf (report, "wrote GDB to %s : ", oneFileName(of)) ;
      reportGdb (gdb, report) ;
    }
}

void gdbDestroy (Gdb *gdb)
{
  dictDestroy (gdb->seqDict) ;
  newFree (gdb->seqLen, gdb->maxSeq, I64) ;
  newFree (gdb->ctgLen, gdb->maxCtg, I64) ;
  newFree (gdb->ctgSeq, gdb->maxCtg, int) ;
  newFree (gdb->ctgPos, gdb->maxCtg, I64) ;
  if (gdb->maxMask)
    { newFree (gdb->ctgMaskCount, gdb->maxCtg, int) ;
      newFree (gdb->ctgMaskStart, gdb->maxCtg, int) ;
      newFree (gdb->mask, gdb->maxMask, I64) ;
    }
  newFree (gdb, 1, Gdb) ;
}

/********************** main function for gdbmask ******************************/

#ifdef GDB_MASK

int main (int argc, char *argv[])
{
  storeCommandLine (argc, argv) ;
  --argc ; ++argv ;

  char *outFileName = 0 ;
  if (argc >= 2 && !strcmp(*argv, "-o"))
    { outFileName = argv[1] ; argc -= 2 ; argv += 2 ; }

  if (argc != 2)
    { fprintf (stderr, "Usage: gdbmask [-o <outfile>] <XX.1gdb> <bed file>\n"
	       "  default is to overwrite the input .1gdb - stash beforehand or use -o to keep\n") ;
      exit (1) ;
    }
     
  if (!outFileName) outFileName = argv[0] ; // seems to be OK

  OneSchema *schema = oneSchemaCreateFromText (schemaText) ;
  OneFile *ofIn = oneFileOpenRead (argv[0],schema, "gdb", 1) ;
  if (!ofIn) die ("failed to open %s as a gdb ONEcode file", argv[0]) ;
  Gdb *gdb = readGdb (ofIn, 1, stdout) ;

  if (gdb->maxMask) // throw away the current mask
    { newFree (gdb->mask, gdb->maxMask, I64) ;
      gdb->maxMask = 0 ;
      gdb->totMask = 0 ;
      memset (gdb->ctgMaskCount, 0, gdb->maxCtg*sizeof(int)) ;
      memset (gdb->ctgMaskStart, 0, gdb->maxCtg*sizeof(int)) ;
    }
  else
    { gdb->ctgMaskCount = new0 (gdb->maxCtg, int) ;
      gdb->ctgMaskStart = new0 (gdb->maxCtg, int) ;
    }

  // next read the bed file
  FILE *f = fopen (argv[1], "r") ;
  if (!f) die ("failed to open %s for reading", argv[1]) ;
  Array    ab = arrayCreate (4096, TanLine) ;
  TanLine *bl = arrayp (ab, arrayMax(ab), TanLine) ;
  char     nameBuf[4096] ;
  while (fscanf (f, "%s\t%lld\t%lld\t%d\t%d\n",
		 nameBuf, &bl->start, &bl->end, &bl->unit, &bl->score) == 5)
    { if (!dictFind (gdb->seqDict, nameBuf, &bl->seq))
	die ("failed to find %s in sequence names for %s", nameBuf, argv[0]) ;
      if (bl->start < 0 || bl->start > bl->end)
	die ("illegal start %lld end %lld line %lld",
	     (long long)bl->start, (long long)bl->end, (long long)arrayMax(ab)) ;
      if (bl->end > gdb->seqLen[bl->seq])
	die ("end %lld off the end %lld of the sequence",
	     (long long)bl->end, (long long) gdb->seqLen[bl->seq]) ;
      bl = arrayp (ab, arrayMax(ab), TanLine) ;
    }
  --arrayMax(ab) ;
  arraySort (ab, tanSort) ; // need to ensure it is sorted
  fprintf (stdout, "read %lld bed lines from %s\n", (long long)arrayMax(ab), argv[1]) ;
  fclose (f) ;

  // now make the new mask, but first need to sort
  gdb->maxMask = 2*arrayMax(ab) ;
  gdb->mask = new (gdb->maxMask, I64) ;
  I64 *m = gdb->mask ;
  TanLine *bi = arrp(ab, 0, TanLine) ;
  TanLine *bend = arrp(ab, arrayMax(ab), TanLine) ; // will not deref bend - just an end marker
  while (bi < bend)
    { TanLine *bj = bi ;
      I64 *m0 = m ;
      while (bj < bend && bj->seq == bi->seq)
	{ *m++ = bj->start ;
	  *m++ = bj->end ;
	  gdb->totMask += (bj->end - bj->start) ;
	  bj++ ;
	}
      gdb->ctgMaskStart[bi->seq] = m0 - gdb->mask ;
      gdb->ctgMaskCount[bi->seq] = m - m0 ;
      bi = bj ;
    }
  arrayDestroy (ab) ;
  
  OneFile *ofOut = oneFileOpenWriteNew (outFileName, schema, "gdb", true, 1) ;
  if (!ofOut) die ("failed to open %s as output", outFileName) ;
  oneInheritProvenance (ofOut, ofIn) ;
  oneAddProvenance (ofOut, "gdbmask", VERSION, getCommandLine()) ;
  writeGdb (ofOut, gdb, 1, stdout) ;
  oneFileClose (ofOut) ;
  if (!strcmp (outFileName, argv[0])) printf ("!! rerun GIXmake %s to apply new mask\n", argv[0]) ;

  gdbDestroy (gdb) ;
  oneFileClose (ofIn) ;

  destroyCommandLine () ;
  return 0 ; 
}

#endif

// I wrote the next routine, but never used it, so have not debugged it yet

/*
OneFile *gdbFile (OneFile *ofAln, int number)
{
  char *source = 0 ;
  char *cpath = 0 ;
  int i ;
  for (i = 0 ; i < oneReferenceCount(of) ; ++i)
    if (ofAln->reference[i].count == number) source = ofAln->reference[i].filename ;
    else if (ofAln->reference[i].count == 3) cpath = ofAln->reference[i].filename ;
  if (!source) die ("failed to find information on source %d in %s", number, of->fileName) ;
  OneFile *ofGdb = oneFileOpenRead (source, schema, "gdb", 1) ;
  if (ofGdb) return ofGdb ;
  
  char buf1[2048] ;
  if (snprintf (buf, 2040, "FAtoGDB %s", source) > 2040)
    die ("name %s too long for FAtoGDB call (>2040 chars)", source) ;
  if (system (buf)) // failed
    { if (!cpath) die ("failed to find source directory in %s", of->fileName) ;
      if (snprintf (buf, 2040, "FAtoGDB %s/%s", cpath, source) > 2040)
	die ("name %s/%s too long for FAtoGDB call (>2040 chars)", cpath, source) ;
      if (system (buf)) // failed again
	die ("The .1aln file %s referred to a sequence file %s but we need\n"
	     "to mask a .1gdb file. We tried to run FAtoGDB to make the .1gdb, but that failed\n"
	     "both on %s and %s/%s. Please make your own .1gdb file and rerun FasTAN on that.\n",
	     source, source, cpath, source) ;
    }

  // now find the GDB that we made
  source = buf + 8 ;                                   // after "FAtoGDB "
  char *s = source ; while (*s) s++ ;                  // go to end
  if (s - source > 3 && !strcmp (s-3, ".gz")) s -= 3 ; // strip off .gz
  while (s > source && s != '.') --s ;                 // go back to previous .
  sprintf (s, ".1gdb") ;                               // add .1gdb on the end
  *s = '.' ;
  ofGdb = oneFileOpenRead (buf, schema, "gdb", 1) ;
  if (ofGdb) return ofGdb ;
  die ("The .1aln file %s referred to a sequence file %s but we need\n"
       "to mask a .1gdb file. We successfully ran FAtoGDB to make the .1gdb, but can't find it.\n"
       "Please rerun FasTAN on the .1gdb file and rerun gdbmask on that.\n") ;
  return 0 ; // for compiler happiness
}
*/

/******************* end of file *******************/

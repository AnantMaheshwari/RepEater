/*  File: taco.c
 *  Authors: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description: TAndem COmpressor
 *   Subcommands: compress, collapse, reconstruct, chain, merge, lift, info
 * HISTORY:
 * Last edited: Nov 16 12:57 2025 (rd109)
 * Created: Sat Nov 15 23:37:43 2025 (rd109)
 * Update: Mon Apr 6 20:20 2025 (am3320) - composable .1taco + .1seq system
 * Update: Aug 23 2026 (am3320) - compress: tolerate trailing N absent from the GDB
 *-------------------------------------------------------------------
 */

#include "alntools.h"
#include "seqio.h"
#include "salsa.h"

/****************** .1seq schema with scaffold structure ******************/
/* Matches FasTAN/FastLTR's expected layout so the FA reader and .1seq reader
   produce the same contig structure: scaffolds split at runs of N. */

static char *scaffoldSeqSchemaText =
  "1 3 def 2 1\n"
  ".\n"
  "P 3 seq                     SEQUENCE\n"
  "O s 2 3 INT 6 STRING        length and id for a scaffold (group of contigs)\n"
  "G S                         scaffolds (s) group sequence objects (S)\n"
  "D n 2 4 CHAR 3 INT          non-acgt run between contigs within a scaffold\n"
  "O S 1 3 DNA                 sequence (one contig)\n" ;

/* Write one scaffold's compressed sequence as s + alternating S/n lines,
   splitting at runs of 'N' or 'n' so the .1seq matches what a FASTA reader
   would build via N-gap splitting.  Caller is responsible for any -f FASTA. */
static void writeScaffoldSeq (OneFile *ofSeq, char *seq, I64 len, char *id)
{
  oneInt(ofSeq, 0) = len ;
  oneWriteLine (ofSeq, 's', strlen(id), id) ;

  I64 i = 0 ;
  while (i < len)
    { I64 contigStart = i ;
      while (i < len && seq[i] != 'N' && seq[i] != 'n') ++i ;
      I64 contigLen = i - contigStart ;
      if (contigLen > 0)
        oneWriteLine (ofSeq, 'S', contigLen, seq + contigStart) ;

      I64 gapStart = i ;
      char gapChar = (i < len) ? seq[i] : 'N' ;  // preserve case from input
      while (i < len && (seq[i] == 'N' || seq[i] == 'n')) ++i ;
      I64 gapLen = i - gapStart ;
      if (gapLen > 0)
        { oneChar(ofSeq, 0) = gapChar ;
          oneInt(ofSeq, 1) = gapLen ;
          oneWriteLine (ofSeq, 'n', 0, 0) ;
        }
    }
}

/****************** compress subcommand ******************/

static int tacoCompress (int argc, char *argv[])
{
  char *outPrefix = 0 ;
  bool emitFasta = false ;
  while (argc >= 2 && **argv == '-')
    { if (!strcmp(*argv, "-o")) { outPrefix = argv[1] ; argc -= 2 ; argv += 2 ; }
      else if (!strcmp(*argv, "-f")) { emitFasta = true ; --argc ; ++argv ; }
      else break ;
    }

  if (argc != 2)
    { fprintf (stderr, "Usage: taco compress [-o <prefix>] [-f] <input.1aln> <seqFile>\n"
	       "  Generalized repeat (e.g. detected by FasTAN/FastLTR) compression and composition system.\n"
	       "  Outputs: prefix.1taco (coordinate map) + prefix.1seq (sequences)\n"
	       "  -f: also emit FASTA (prefix-taco.fa.gz)\n ") ;
      return 1 ;
    }

  /* open the input .1aln file */
  OneSchema *schema = oneSchemaCreateFromText (schemaText) ;
  OneFile *ofIn = oneFileOpenRead (argv[0], schema, "aln", 1) ;
  if (!ofIn) die ("failed to open %s as a .1aln file", argv[0]) ;
  Gdb *gdb = readGdb (ofIn, 1, stderr) ;
  if (ofIn->lineType != 'A') die ("unexpected line type %c", ofIn->lineType) ;

  /* open the input sequence */
  SeqIO *inIO = seqIOopenRead (argv[1], dna2textConv, false) ;
  if (!inIO) die ("failed to open %s to read as a sequence file", argv[1]) ;

  /* determine output prefix */
  char prefixBuf[1024] ;
  if (!outPrefix)
    { strncpy (prefixBuf, argv[1], sizeof(prefixBuf) - 20) ;
      char *s = prefixBuf + strlen(prefixBuf) ;
      if (s > prefixBuf+3 && !strcmp(s-3, ".gz")) { s -= 3 ; *s = 0 ; }
      while (s > prefixBuf && *s != '.') --s ;
      if (s == prefixBuf) s += strlen(prefixBuf) ; else *s = 0 ;
      outPrefix = prefixBuf ;
    }

  /* open .1taco output */
  char tacoFileName[1024] ;
  snprintf (tacoFileName, sizeof(tacoFileName), "%s.1taco", outPrefix) ;
  OneSchema *mapSchema = oneSchemaCreateFromText (salsaSchemaText) ;
  OneFile *ofTaco = oneFileOpenWriteNew (tacoFileName, mapSchema, "taco", true, 1) ;
  oneSchemaDestroy (mapSchema) ;
  if (!ofTaco) die ("failed to open %s to write", tacoFileName) ;
  oneAddProvenance (ofTaco, "taco", VERSION, getCommandLine()) ;

  /* open .1seq output (scaffold-aware: s + S/n alternating) */
  char seqFileName[1024] ;
  snprintf (seqFileName, sizeof(seqFileName), "%s.1seq", outPrefix) ;
  OneSchema *seqSchema = oneSchemaCreateFromText (scaffoldSeqSchemaText) ;
  OneFile *ofSeq = oneFileOpenWriteNew (seqFileName, seqSchema, "seq", true, 1) ;
  oneSchemaDestroy (seqSchema) ;
  if (!ofSeq) die ("failed to open %s to write", seqFileName) ;
  oneAddProvenance (ofSeq, "taco", VERSION, getCommandLine()) ;

  /* optionally open FASTA output */
  SeqIO *fastaIO = 0 ;
  if (emitFasta)
    { char fastaName[1024] ;
      snprintf (fastaName, sizeof(fastaName), "%s-taco.fa.gz", outPrefix) ;
      fastaIO = seqIOopenWrite (fastaName, 0, dna2textConv, 0) ;
      if (!fastaIO) die ("failed to open %s to write", fastaName) ;
    }

  /* read the .1aln file */
  I64 nAlign = 0 ;
  oneStats (ofIn, 'A', &nAlign, 0, 0) ;
  Array at = arrayCreate (nAlign, TanLine) ;
  while (ofIn->lineType == 'A')
    { TanLine *t = arrayp(at, arrayMax(at), TanLine) ;
      int ctg = oneInt(ofIn,0) ;
      if (ctg != oneInt(ofIn,3))
	die ("target mismatch line %lld - not a TAN file?", (long long)ofIn->line) ;
      t->seq     = ctg2seq(gdb,ctg) ;
      t->start   = ctg2pos(gdb,ctg,oneInt(ofIn,4)) ;
      t->end     = ctg2pos(gdb,ctg,oneInt(ofIn,2)) ;
      double len = oneInt(ofIn,2) - oneInt(ofIn,4) ;
      while (oneReadLine(ofIn) && ofIn->lineType != 'A')
	if (ofIn->lineType == 'D') t->score = (int)(1000*(1.0 - oneInt(ofIn,0)/len)) ;
	else if (ofIn->lineType == 'U') t->unit = oneInt(ofIn,0) ;
    }
  oneFileClose (ofIn) ;
  arraySort (at, tanSort) ;

  /* build seqStart lookup */
  int i ;
  I32 *seqStart = new0 (gdb->nSeq, I32) ;
  I32 lastSeq = -1 ;
  for (i = 0 ; i < arrayMax(at) ; ++i)
    { I32 seq = arrp(at, i, TanLine)->seq ;
      if (seq != lastSeq) { seqStart[seq] = i ; lastSeq = seq ; }
    }

  /* allocate buffers */
  I64 maxSeqLen = 0 ;
  for (i = 0 ; i < gdb->nSeq ; ++i)
    if (gdb->seqLen[i] > maxSeqLen) maxSeqLen = gdb->seqLen[i] ;
  char *tacoSeq = new (maxSeqLen, char) ;
  Array seqLift = arrayCreate (64, SalsaEntry) ;

  /* process each sequence */
  while (seqIOread (inIO))
    { U32 seq ;
      char *id = sqioId(inIO) ;
      if (inIO->descLen)
	{ id = malloc (strlen(id) + strlen(sqioDesc(inIO)) + 2) ;
	  sprintf (id, "%s %s", sqioId(inIO), sqioDesc(inIO)) ;
	}
      if (!dictFind (gdb->seqDict, id, &seq))
	die ("failed to match name %s in %s to %s", sqioId(inIO), argv[1], argv[0]) ;
      if (id != sqioId(inIO)) free (id) ;
      /* FasTAN's GDB builder closes a scaffold at its last contig, so a run of N at
	 the end of a scaffold is absent from gdb->seqLen (a leading run is kept, as a
	 zero-length contig plus a gap).  Tolerate a source longer than the GDB when
	 the surplus is all N: compress just the first gdb->seqLen[seq] bases, which is
	 what the GDB and all downstream coordinates are relative to. */
      if (inIO->seqLen > (U64) gdb->seqLen[seq])
	{ char *s = sqioSeq(inIO) ;
	  I64   j ;
	  for (j = gdb->seqLen[seq] ; j < (I64) inIO->seqLen ; ++j)
	    if (s[j] != 'N' && s[j] != 'n') break ;
	  if (j < (I64) inIO->seqLen)
	    die ("length mismatch for seq %s: %lld in %s, %lld in %s - surplus is not all N",
		 sqioId(inIO), (long long) gdb->seqLen[seq], argv[0],
		 (long long) inIO->seqLen, argv[1]) ;
	  warn ("seq %s: ignoring %lld trailing N in %s that are absent from the GDB in %s",
		sqioId(inIO), (long long) (inIO->seqLen - gdb->seqLen[seq]), argv[1], argv[0]) ;
	}
      else if (inIO->seqLen != (U64) gdb->seqLen[seq])
	die ("length mismatch for seq %s: %lld in %s, %lld in %s",
	     sqioId(inIO), (long long) gdb->seqLen[seq], argv[0],
	     (long long) inIO->seqLen, argv[1]) ;

      U64   nNew = 0, nOld = 0 ;
      char *oldSeq = sqioSeq(inIO) ;
      arrayMax(seqLift) = 0 ;

      for (i = seqStart[seq] ; i < arrayMax(at) && arrp(at,i,TanLine)->seq == seq ; ++i)
	{ TanLine *t = arrp(at,i,TanLine) ;
	  if (t->start < nOld)
	    { t->start = nOld ;
	      if (t->start + t->unit > t->end)
		warn ("buried repeat %s size %d ending %lld when last ended at %lld",
		      sqioId(inIO), t->unit, (long long) t->end, (long long) nOld) ;
	      continue ;
	    }
	  SalsaEntry *e = arrayp (seqLift, arrayMax(seqLift), SalsaEntry) ;
	  e->orig_start = t->start ;
	  e->orig_end   = t->end ;
	  e->taco_pos   = nNew + (t->start - nOld) ;
	  e->unit_len   = t->unit ;
	  /* no need to copy removed DNA here; oldSeq stays valid until end of seqIOread */
	  memcpy (tacoSeq+nNew, oldSeq+nOld, t->start - nOld + t->unit) ;
	  nNew += t->start - nOld + t->unit ;
	  nOld = t->end ;
	}
      memcpy (tacoSeq+nNew, oldSeq+nOld, gdb->seqLen[seq] - nOld) ;
      nNew += gdb->seqLen[seq] - nOld ;

      /* write .1taco: c line + I line + L lines */
      oneInt(ofTaco,0) = gdb->seqLen[seq] ;
      oneInt(ofTaco,1) = nNew ;
      oneWriteLine (ofTaco, 'c', 0, 0) ;
      oneWriteLine (ofTaco, 'I', strlen(sqioId(inIO)), sqioId(inIO)) ;

      /* write .1seq: scaffold (s) + contigs (S) + gaps (n) */
      writeScaffoldSeq (ofSeq, tacoSeq, nNew, sqioId(inIO)) ;

      /* write L events + removed DNA S lines to .1taco */
      int j ;
      for (j = 0 ; j < arrayMax(seqLift) ; ++j)
	{ SalsaEntry *e = arrp(seqLift, j, SalsaEntry) ;
	  oneInt(ofTaco,0) = e->orig_start ;
	  oneInt(ofTaco,1) = e->orig_end ;
	  oneInt(ofTaco,2) = e->taco_pos ;
	  oneInt(ofTaco,3) = e->unit_len ;
	  oneWriteLine (ofTaco, 'L', 0, 0) ;
	  I64 removedLen = e->orig_end - e->orig_start - e->unit_len ;
	  if (removedLen > 0)
	    oneWriteLine (ofTaco, 'S', removedLen, oldSeq + e->orig_start + e->unit_len) ;
	}

      /* optionally write FASTA */
      if (fastaIO)
	seqIOwrite (fastaIO, sqioId(inIO), inIO->descLen ? sqioDesc(inIO) : 0, nNew, tacoSeq, 0) ;
    }

  /* cleanup */
  arrayDestroy (seqLift) ;
  newFree (tacoSeq, maxSeqLen, char) ;
  newFree (seqStart, gdb->nSeq, I32) ;
  arrayDestroy (at) ;
  gdbDestroy (gdb) ;
  seqIOclose (inIO) ;
  oneFileClose (ofTaco) ;
  oneFileClose (ofSeq) ;
  if (fastaIO) seqIOclose (fastaIO) ;

  return 0 ;
}

/****************** collapse subcommand ******************/

typedef struct {
  I32 seq ;
  I64 start, end ;
  I32 unit ;
} CollapseEvent ;

static int eventSort (const void *a, const void *b)
{
  CollapseEvent *ea = (CollapseEvent*)a, *eb = (CollapseEvent*)b ;
  if (ea->seq != eb->seq) return ea->seq - eb->seq ;
  if (ea->start != eb->start) return (ea->start < eb->start) ? -1 : 1 ;
  return (ea->end < eb->end) ? -1 : (ea->end > eb->end) ? 1 : 0 ;
}

static int tacoCollapse (int argc, char *argv[])
{
  char *outPrefix = 0 ;
  bool emitFasta = false ;
  while (argc >= 2 && **argv == '-')
    { if (!strcmp(*argv, "-o")) { outPrefix = argv[1] ; argc -= 2 ; argv += 2 ; }
      else if (!strcmp(*argv, "-f")) { emitFasta = true ; --argc ; ++argv ; }
      else break ;
    }

  if (argc != 2)
    { fprintf (stderr, "Usage: taco collapse [-o <prefix>] [-f] <spec.tsv> <input.fa>\n"
	       "  spec.tsv: tab-separated lines: seq_index orig_start orig_end unit_len\n"
	       "    unit_len > 0: tandem compression (keep first unit_len bases)\n"
	       "    unit_len = 0: full deletion (remove entire interval)\n"
	       "  Outputs: prefix.1taco + prefix.1seq\n"
	       "  -f: also emit FASTA (prefix_collapsed.fa.gz)\n") ;
      return 1 ;
    }

  /* read the TSV collapse spec */
  FILE *specFile = fopen (argv[0], "r") ;
  if (!specFile) die ("failed to open spec file %s", argv[0]) ;
  Array events = arrayCreate (1024, CollapseEvent) ;
  char line[1024] ;
  int nSeq = 0 ;
  while (fgets (line, sizeof(line), specFile))
    { if (line[0] == '#' || line[0] == '\n') continue ;
      CollapseEvent *e = arrayp (events, arrayMax(events), CollapseEvent) ;
      if (sscanf (line, "%d %lld %lld %d", &e->seq, (long long*)&e->start,
		  (long long*)&e->end, &e->unit) != 4)
	die ("bad line in spec file: %s", line) ;
      if (e->seq >= nSeq) nSeq = e->seq + 1 ;
    }
  fclose (specFile) ;
  arraySort (events, eventSort) ;

  /* open input sequence file */
  SeqIO *inIO = seqIOopenRead (argv[1], dna2textConv, false) ;
  if (!inIO) die ("failed to open %s to read as a sequence file", argv[1]) ;

  /* determine output prefix */
  char prefixBuf[1024] ;
  if (!outPrefix)
    { strncpy (prefixBuf, argv[1], sizeof(prefixBuf) - 20) ;
      char *s = prefixBuf + strlen(prefixBuf) ;
      if (s > prefixBuf+3 && !strcmp(s-3, ".gz")) { s -= 3 ; *s = 0 ; }
      while (s > prefixBuf && *s != '.') --s ;
      if (s == prefixBuf) s += strlen(prefixBuf) ; else *s = 0 ;
      outPrefix = prefixBuf ;
    }

  /* open .1taco output */
  char tacoFileName[1024] ;
  snprintf (tacoFileName, sizeof(tacoFileName), "%s.1taco", outPrefix) ;
  OneSchema *mapSchema = oneSchemaCreateFromText (salsaSchemaText) ;
  OneFile *ofTaco = oneFileOpenWriteNew (tacoFileName, mapSchema, "taco", true, 1) ;
  oneSchemaDestroy (mapSchema) ;
  if (!ofTaco) die ("failed to open %s to write", tacoFileName) ;
  oneAddProvenance (ofTaco, "taco", VERSION, getCommandLine()) ;

  /* open .1seq output (scaffold-aware: s + S/n alternating) */
  char seqFileName[1024] ;
  snprintf (seqFileName, sizeof(seqFileName), "%s.1seq", outPrefix) ;
  OneSchema *seqSchema = oneSchemaCreateFromText (scaffoldSeqSchemaText) ;
  OneFile *ofSeq = oneFileOpenWriteNew (seqFileName, seqSchema, "seq", true, 1) ;
  oneSchemaDestroy (seqSchema) ;
  if (!ofSeq) die ("failed to open %s to write", seqFileName) ;
  oneAddProvenance (ofSeq, "taco", VERSION, getCommandLine()) ;

  /* optionally open FASTA output */
  SeqIO *fastaIO = 0 ;
  if (emitFasta)
    { char fastaName[1024] ;
      snprintf (fastaName, sizeof(fastaName), "%s_collapsed.fa.gz", outPrefix) ;
      fastaIO = seqIOopenWrite (fastaName, 0, dna2textConv, 0) ;
      if (!fastaIO) die ("failed to open %s to write", fastaName) ;
    }

  /* build seqStart lookup */
  int i ;
  I32 *seqStartArr = new0 (nSeq, I32) ;
  I32 lastSeq = -1 ;
  for (i = 0 ; i < arrayMax(events) ; ++i)
    { I32 seq = arrp(events, i, CollapseEvent)->seq ;
      if (seq != lastSeq) { seqStartArr[seq] = i ; lastSeq = seq ; }
    }

  Array seqLift = arrayCreate (64, SalsaEntry) ;

  /* process sequences */
  int seqIdx = 0 ;
  while (seqIOread (inIO))
    { I64 seqLen = inIO->seqLen ;
      char *outBuf = new (seqLen, char) ;
      char *inSeq = sqioSeq(inIO) ;
      U64 nNew = 0, nOld = 0 ;
      arrayMax(seqLift) = 0 ;

      if (seqIdx < nSeq)
	for (i = seqStartArr[seqIdx] ;
	     i < arrayMax(events) && arrp(events,i,CollapseEvent)->seq == seqIdx ; ++i)
	  { CollapseEvent *ce = arrp(events, i, CollapseEvent) ;
		    if (ce->start < 0 || ce->end > seqLen || ce->start >= ce->end)
		      { warn ("out-of-range event seq %d [%lld,%lld) vs scaffold len %lld; skipping",
			      seqIdx, (long long)ce->start, (long long)ce->end, (long long)seqLen) ;
			continue ;
		      }
	    if (ce->start < (I64)nOld)
	      { if (ce->start + ce->unit > ce->end)
		  warn ("buried event seq %d at %lld when cursor at %lld",
			seqIdx, (long long)ce->start, (long long)nOld) ;
		continue ;
	      }
	    SalsaEntry *se = arrayp (seqLift, arrayMax(seqLift), SalsaEntry) ;
	    se->orig_start = ce->start ;
	    se->orig_end   = ce->end ;
	    se->taco_pos   = nNew + (ce->start - nOld) ;
	    se->unit_len   = ce->unit ;
	    /* no need to copy removed DNA; inSeq stays valid until end of seqIOread */
	    memcpy (outBuf+nNew, inSeq+nOld, ce->start - nOld) ;
	    nNew += ce->start - nOld ;
	    if (ce->unit > 0)
	      { memcpy (outBuf+nNew, inSeq+ce->start, ce->unit) ;
		nNew += ce->unit ;
	      }
	    nOld = ce->end ;
	  }

      memcpy (outBuf+nNew, inSeq+nOld, seqLen - nOld) ;
      nNew += seqLen - nOld ;

      /* write .1taco */
      oneInt(ofTaco,0) = seqLen ;
      oneInt(ofTaco,1) = nNew ;
      oneWriteLine (ofTaco, 'c', 0, 0) ;
      oneWriteLine (ofTaco, 'I', strlen(sqioId(inIO)), sqioId(inIO)) ;

      /* write .1seq: scaffold (s) + contigs (S) + gaps (n) */
      writeScaffoldSeq (ofSeq, outBuf, nNew, sqioId(inIO)) ;

      /* write L events + removed DNA S lines to .1taco */
      int j ;
      for (j = 0 ; j < arrayMax(seqLift) ; ++j)
	{ SalsaEntry *se = arrp(seqLift, j, SalsaEntry) ;
	  oneInt(ofTaco,0) = se->orig_start ;
	  oneInt(ofTaco,1) = se->orig_end ;
	  oneInt(ofTaco,2) = se->taco_pos ;
	  oneInt(ofTaco,3) = se->unit_len ;
	  oneWriteLine (ofTaco, 'L', 0, 0) ;
	  I64 removedLen = se->orig_end - se->orig_start - se->unit_len ;
	  if (removedLen > 0)
	    oneWriteLine (ofTaco, 'S', removedLen, inSeq + se->orig_start + se->unit_len) ;
	}

      if (fastaIO)
	seqIOwrite (fastaIO, sqioId(inIO), inIO->descLen ? sqioDesc(inIO) : 0, nNew, outBuf, 0) ;

      newFree (outBuf, seqLen, char) ;
      ++seqIdx ;
    }

  arrayDestroy (seqLift) ;
  if (nSeq) newFree (seqStartArr, nSeq, I32) ;
  arrayDestroy (events) ;
  seqIOclose (inIO) ;
  oneFileClose (ofTaco) ;
  oneFileClose (ofSeq) ;
  if (fastaIO) seqIOclose (fastaIO) ;

  return 0 ;
}

/****************** reconstruct subcommand ******************/

static int tacoReconstruct (int argc, char *argv[])
{
  char *outFileName = 0 ;
  while (argc >= 2 && **argv == '-')
    { if (!strcmp(*argv, "-o")) { outFileName = argv[1] ; argc -= 2 ; argv += 2 ; }
      else break ;
    }

  if (argc != 2)
    { fprintf (stderr, "Usage: taco reconstruct [-o <output.fa>] <input.1taco> <input.1seq>\n"
	       "  Reconstruct original FASTA from .1taco + .1seq pair.\n") ;
      return 1 ;
    }

  /* open .1taco */
  OneSchema *mapSchema = oneSchemaCreateFromText (salsaSchemaText) ;
  OneFile *ofTaco = oneFileOpenRead (argv[0], mapSchema, "taco", 1) ;
  oneSchemaDestroy (mapSchema) ;
  if (!ofTaco) die ("failed to open %s", argv[0]) ;
  I64 nSeq = 0 ;
  oneStats (ofTaco, 'c', &nSeq, 0, 0) ;

  /* open .1seq */
  OneSchema *seqSchema = oneSchemaCreateFromText (seqioSchemaText) ;
  OneFile *ofSeq = oneFileOpenRead (argv[1], seqSchema, "seq", 1) ;
  oneSchemaDestroy (seqSchema) ;
  if (!ofSeq) die ("failed to open %s", argv[1]) ;

  /* read coordinate map (with removed DNA) + compressed seqs */
  Salsa *salsas = salsaReadWithSeq (ofTaco, ofSeq, (int)nSeq) ;
  oneFileClose (ofTaco) ;
  oneFileClose (ofSeq) ;

  /* open output */
  if (!outFileName) outFileName = "-" ; /* stdout */
  SeqIO *outIO = seqIOopenWrite (outFileName, 0, dna2textConv, 0) ;
  if (!outIO) die ("failed to open %s to write", outFileName) ;

  int i ;
  for (i = 0 ; i < (int)nSeq ; ++i)
    { Salsa *s = &salsas[i] ;
      char *orig = salsaReconstruct (s) ;
      seqIOwrite (outIO, s->seqName ? s->seqName : "unknown", 0, s->origLen, orig, 0) ;
      free (orig) ;
    }

  seqIOclose (outIO) ;
  salsaDestroy (salsas, (int)nSeq) ;
  return 0 ;
}

/****************** chain subcommand ******************/

static int tacoChain (int argc, char *argv[])
{
  char *outPrefix = 0 ;
  while (argc >= 2 && **argv == '-')
    { if (!strcmp(*argv, "-o")) { outPrefix = argv[1] ; argc -= 2 ; argv += 2 ; }
      else break ;
    }

  if (argc < 4 || argc % 2 != 0)
    { fprintf (stderr, "Usage: taco chain [-o <prefix>] <a.1taco> <a.1seq> <b.1taco> <b.1seq> [...]\n"
	       "  Flatten multiple .1taco+.1seq layers into a single pair.\n"
	       "  Files given in application order (a applied first, b next, etc.).\n") ;
      return 1 ;
    }

  int nLayers = argc / 2 ;

  /* read first layer */
  OneSchema *mapSchema = oneSchemaCreateFromText (salsaSchemaText) ;
  OneSchema *seqSchema = oneSchemaCreateFromText (seqioSchemaText) ;

  OneFile *ofTaco = oneFileOpenRead (argv[0], mapSchema, "taco", 1) ;
  if (!ofTaco) die ("failed to open %s", argv[0]) ;
  I64 nSeq = 0 ;
  oneStats (ofTaco, 'c', &nSeq, 0, 0) ;
  OneFile *ofSeq = oneFileOpenRead (argv[1], seqSchema, "seq", 1) ;
  if (!ofSeq) die ("failed to open %s", argv[1]) ;
  Salsa *current = salsaReadWithSeq (ofTaco, ofSeq, (int)nSeq) ;
  int currentN = (int)nSeq ;
  oneFileClose (ofTaco) ;
  oneFileClose (ofSeq) ;

  /* compose each subsequent layer */
  int k ;
  for (k = 1 ; k < nLayers ; ++k)
    { oneSchemaDestroy (mapSchema) ; mapSchema = oneSchemaCreateFromText (salsaSchemaText) ;
      oneSchemaDestroy (seqSchema) ; seqSchema = oneSchemaCreateFromText (seqioSchemaText) ;
      ofTaco = oneFileOpenRead (argv[k*2], mapSchema, "taco", 1) ;
      if (!ofTaco) die ("failed to open %s", argv[k*2]) ;
      I64 nSeq2 = 0 ;
      oneStats (ofTaco, 'c', &nSeq2, 0, 0) ;
      if (!nSeq2) oneStats (ofTaco, 's', &nSeq2, 0, 0) ;
      if (nSeq2 != currentN) die ("sequence count mismatch: %d vs %d at layer %d", currentN, (int)nSeq2, k) ;
      ofSeq = oneFileOpenRead (argv[k*2+1], seqSchema, "seq", 1) ;
      if (!ofSeq) die ("failed to open %s", argv[k*2+1]) ;
      Salsa *outer = salsaReadWithSeq (ofTaco, ofSeq, (int)nSeq2) ;
      oneFileClose (ofTaco) ;
      oneFileClose (ofSeq) ;

      /* compose for each sequence */
      Salsa *composed = new0 (currentN, Salsa) ;
      int i ;
      for (i = 0 ; i < currentN ; ++i)
	composed[i] = salsaCompose (&current[i], &outer[i]) ;
      salsaDestroy (current, currentN) ;
      salsaDestroy (outer, (int)nSeq2) ;
      current = composed ;
    }
  oneSchemaDestroy (mapSchema) ;
  oneSchemaDestroy (seqSchema) ;

  /* write output */
  if (!outPrefix) outPrefix = "chained" ;
  char tacoFileName[1024], seqFileName2[1024] ;
  snprintf (tacoFileName, sizeof(tacoFileName), "%s.1taco", outPrefix) ;
  snprintf (seqFileName2, sizeof(seqFileName2), "%s.1seq", outPrefix) ;

  mapSchema = oneSchemaCreateFromText (salsaSchemaText) ;
  ofTaco = oneFileOpenWriteNew (tacoFileName, mapSchema, "taco", true, 1) ;
  oneSchemaDestroy (mapSchema) ;
  if (!ofTaco) die ("failed to open %s to write", tacoFileName) ;
  oneAddProvenance (ofTaco, "taco", VERSION, getCommandLine()) ;

  seqSchema = oneSchemaCreateFromText (seqioSchemaText) ;
  ofSeq = oneFileOpenWriteNew (seqFileName2, seqSchema, "seq", true, 1) ;
  oneSchemaDestroy (seqSchema) ;
  if (!ofSeq) die ("failed to open %s to write", seqFileName2) ;
  oneAddProvenance (ofSeq, "taco", VERSION, getCommandLine()) ;

  salsaWriteWithSeq (ofTaco, ofSeq, current, currentN) ;
  oneFileClose (ofTaco) ;
  oneFileClose (ofSeq) ;
  salsaDestroy (current, currentN) ;

  return 0 ;
}

/****************** merge subcommand ******************/

static int tacoMergeCmd (int argc, char *argv[])
{
  char *outPrefix = 0 ;
  while (argc >= 2 && **argv == '-')
    { if (!strcmp(*argv, "-o")) { outPrefix = argv[1] ; argc -= 2 ; argv += 2 ; }
      else break ;
    }

  if (argc != 4)
    { fprintf (stderr, "Usage: taco merge [-o <prefix>] <a.1taco> <a.1seq> <b.1taco> <b.1seq>\n"
	       "  Merge two non-overlapping .1taco+.1seq pairs on the same original.\n") ;
      return 1 ;
    }

  OneSchema *mapSchema = oneSchemaCreateFromText (salsaSchemaText) ;
  OneSchema *seqSchema = oneSchemaCreateFromText (seqioSchemaText) ;

  /* read a */
  OneFile *ofTaco = oneFileOpenRead (argv[0], mapSchema, "taco", 1) ;
  if (!ofTaco) die ("failed to open %s", argv[0]) ;
  I64 nSeqA = 0 ;
  oneStats (ofTaco, 'c', &nSeqA, 0, 0) ;
  if (!nSeqA) oneStats (ofTaco, 's', &nSeqA, 0, 0) ;
  OneFile *ofSeq = oneFileOpenRead (argv[1], seqSchema, "seq", 1) ;
  if (!ofSeq) die ("failed to open %s", argv[1]) ;
  Salsa *sA = salsaReadWithSeq (ofTaco, ofSeq, (int)nSeqA) ;
  oneFileClose (ofTaco) ;
  oneFileClose (ofSeq) ;

  /* read b */
  oneSchemaDestroy (mapSchema) ; mapSchema = oneSchemaCreateFromText (salsaSchemaText) ;
  oneSchemaDestroy (seqSchema) ; seqSchema = oneSchemaCreateFromText (seqioSchemaText) ;
  ofTaco = oneFileOpenRead (argv[2], mapSchema, "taco", 1) ;
  if (!ofTaco) die ("failed to open %s", argv[2]) ;
  I64 nSeqB = 0 ;
  oneStats (ofTaco, 'c', &nSeqB, 0, 0) ;
  if (!nSeqB) oneStats (ofTaco, 's', &nSeqB, 0, 0) ;
  if (nSeqA != nSeqB) die ("sequence count mismatch: %d vs %d", (int)nSeqA, (int)nSeqB) ;
  ofSeq = oneFileOpenRead (argv[3], seqSchema, "seq", 1) ;
  if (!ofSeq) die ("failed to open %s", argv[3]) ;
  Salsa *sB = salsaReadWithSeq (ofTaco, ofSeq, (int)nSeqB) ;
  oneFileClose (ofTaco) ;
  oneFileClose (ofSeq) ;
  oneSchemaDestroy (mapSchema) ;
  oneSchemaDestroy (seqSchema) ;

  /* merge */
  Salsa *merged = new0 ((int)nSeqA, Salsa) ;
  int i ;
  for (i = 0 ; i < (int)nSeqA ; ++i)
    merged[i] = salsaMerge (&sA[i], &sB[i]) ;
  salsaDestroy (sA, (int)nSeqA) ;
  salsaDestroy (sB, (int)nSeqB) ;

  /* write output */
  if (!outPrefix) outPrefix = "merged" ;
  char tacoFileName[1024], seqFileName2[1024] ;
  snprintf (tacoFileName, sizeof(tacoFileName), "%s.1taco", outPrefix) ;
  snprintf (seqFileName2, sizeof(seqFileName2), "%s.1seq", outPrefix) ;

  mapSchema = oneSchemaCreateFromText (salsaSchemaText) ;
  ofTaco = oneFileOpenWriteNew (tacoFileName, mapSchema, "taco", true, 1) ;
  oneSchemaDestroy (mapSchema) ;
  oneAddProvenance (ofTaco, "taco", VERSION, getCommandLine()) ;

  seqSchema = oneSchemaCreateFromText (seqioSchemaText) ;
  ofSeq = oneFileOpenWriteNew (seqFileName2, seqSchema, "seq", true, 1) ;
  oneSchemaDestroy (seqSchema) ;
  oneAddProvenance (ofSeq, "taco", VERSION, getCommandLine()) ;

  salsaWriteWithSeq (ofTaco, ofSeq, merged, (int)nSeqA) ;
  oneFileClose (ofTaco) ;
  oneFileClose (ofSeq) ;
  salsaDestroy (merged, (int)nSeqA) ;

  return 0 ;
}

/****************** lift subcommand ******************/

static int tacoLift (int argc, char *argv[])
{
  bool reverse = false ;
  while (argc >= 1 && **argv == '-')
    { if (!strcmp(*argv, "-r")) { reverse = true ; --argc ; ++argv ; }
      else break ;
    }

  if (argc < 2)
    { fprintf (stderr, "Usage: taco lift [-r] <annotations.bed> <map.1taco> [map2.1taco ...]\n"
	       "  default: final -> original (lift collapsed coords back to original)\n"
	       "  -r: original -> final (project original coords to collapsed space)\n"
	       "  annotations.bed: tab-separated: seq_index start end [extra columns...]\n"
	       "  .1taco files in round order (round 0 first)\n") ;
      return 1 ;
    }

  char *bedFile = argv[0] ;
  int nMaps = argc - 1 ;
  char **mapFiles = argv + 1 ;

  /* load the salsa stack */
  SalsaStack *stack = salsaStackLoad (mapFiles, nMaps) ;

  /* process annotations */
  FILE *fp = strcmp(bedFile, "-") ? fopen (bedFile, "r") : stdin ;
  if (!fp) die ("failed to open %s", bedFile) ;

  char line[4096] ;
  while (fgets (line, sizeof(line), fp))
    { if (line[0] == '#' || line[0] == '\n')
	{ fputs (line, stdout) ;
	  continue ;
	}
      int seq ;
      long long start, end ;
      int nParsed ;
      char *rest = 0 ;

      if (sscanf (line, "%d\t%lld\t%lld%n", &seq, &start, &end, &nParsed) < 3)
	{ fprintf (stderr, "warning: skipping unparseable line: %s", line) ;
	  continue ;
	}
      rest = line + nParsed ;

      I64 oa, ob ;
      if (reverse)
	salsaStackOrigToFinalInterval (stack, seq, (I64)start, (I64)end, &oa, &ob) ;
      else
	salsaStackFinalToOrigInterval (stack, seq, (I64)start, (I64)end, &oa, &ob) ;

      printf ("%d\t%lld\t%lld", seq, (long long)oa, (long long)ob) ;
      if (rest && *rest) fputs (rest, stdout) ;
      else putchar ('\n') ;
    }

  if (fp != stdin) fclose (fp) ;
  salsaStackDestroy (stack) ;

  return 0 ;
}

/****************** info subcommand ******************/

static int tacoInfo (int argc, char *argv[])
{
  if (argc != 1)
    { fprintf (stderr, "Usage: taco info <input.1taco>\n"
	       "  Print summary statistics of a .1taco file.\n") ;
      return 1 ;
    }

  OneSchema *schema = oneSchemaCreateFromText (salsaSchemaText) ;
  OneFile *of = oneFileOpenRead (argv[0], schema, "taco", 1) ;
  oneSchemaDestroy (schema) ;
  if (!of) die ("failed to open %s", argv[0]) ;

  I64 nSeq = 0, nEvents = 0 ;
  oneStats (of, 'c', &nSeq, 0, 0) ;
  oneStats (of, 'L', &nEvents, 0, 0) ;

  printf ("File: %s\n", argv[0]) ;
  printf ("Sequences: %lld\n", (long long)nSeq) ;
  printf ("Lift events: %lld\n", (long long)nEvents) ;

  /* read and print per-sequence details */
  Salsa *salsas = salsaRead (of, (int)nSeq) ;
  oneFileClose (of) ;

  I64 totOrig = 0, totTaco = 0, totRemoved = 0 ;
  int i, j ;
  for (i = 0 ; i < (int)nSeq ; ++i)
    { Salsa *s = &salsas[i] ;
      printf ("  seq %d", i) ;
      if (s->seqName) printf (" (%s)", s->seqName) ;
      printf (": origLen=%lld tacoLen=%lld events=%d compression=%.1f%%\n",
	      (long long)s->origLen, (long long)s->tacoLen, s->n,
	      s->origLen > 0 ? 100.0 * (1.0 - (double)s->tacoLen / s->origLen) : 0.0) ;
      totOrig += s->origLen ;
      totTaco += s->tacoLen ;
      for (j = 0 ; j < s->n ; ++j)
	totRemoved += s->entry[j].orig_end - s->entry[j].orig_start - s->entry[j].unit_len ;
    }
  printf ("Total: origLen=%lld tacoLen=%lld removed=%lld compression=%.1f%%\n",
	  (long long)totOrig, (long long)totTaco, (long long)totRemoved,
	  totOrig > 0 ? 100.0 * (1.0 - (double)totTaco / totOrig) : 0.0) ;

  salsaDestroy (salsas, (int)nSeq) ;
  return 0 ;
}

/****************** windows subcommand ******************/

static int tacoWindows (int argc, char *argv[])
{
  I64 flank = 0 ;
  while (argc >= 2 && **argv == '-')
    { if (!strcmp(*argv, "-f")) { flank = strtoll (argv[1], 0, 10) ; argc -= 2 ; argv += 2 ; }
      else break ;
    }

  if (argc != 1)
    { fprintf (stderr, "Usage: taco windows [-f <flank>] <input.1taco>\n"
	       "  Emit tab-separated scan windows (name<TAB>beg<TAB>end) around each\n"
	       "  compression/collapse scar, in the compressed (.1seq) coordinate space,\n"
	       "  padded by <flank> on each side and merged per sequence.  Feeds -r of\n"
	       "  FasTAN/FastLTR so the next round only re-scans changed neighborhoods.\n") ;
      return 1 ;
    }

  OneSchema *schema = oneSchemaCreateFromText (salsaSchemaText) ;
  OneFile *of = oneFileOpenRead (argv[0], schema, "taco", 1) ;
  oneSchemaDestroy (schema) ;
  if (!of) die ("failed to open %s", argv[0]) ;

  I64 nSeq = 0 ;
  oneStats (of, 'c', &nSeq, 0, 0) ;
  Salsa *salsas = salsaRead (of, (int)nSeq) ;
  oneFileClose (of) ;

  //  Emit each scar's [taco_pos, taco_pos+unit_len) padded by flank, clamped and
  //  merged per sequence (entry[] is sorted by taco_pos, so windows are ascending).
  int i, j ;
  for (i = 0 ; i < (int)nSeq ; ++i)
    { Salsa *s = &salsas[i] ;
      if (s->seqName == NULL) continue ;
      I64 curBeg = -1, curEnd = -1 ;
      for (j = 0 ; j < s->n ; ++j)
	{ I64 b = s->entry[j].taco_pos - flank ;
	  I64 e = s->entry[j].taco_pos + s->entry[j].unit_len + flank ;
	  if (b < 0) b = 0 ;
	  if (e > s->tacoLen) e = s->tacoLen ;
	  if (b >= e) continue ;
	  if (curBeg < 0) { curBeg = b ; curEnd = e ; }
	  else if (b <= curEnd) { if (e > curEnd) curEnd = e ; }
	  else { printf ("%s\t%lld\t%lld\n", s->seqName, (long long)curBeg, (long long)curEnd) ;
	         curBeg = b ; curEnd = e ; }
	}
      if (curBeg >= 0)
	printf ("%s\t%lld\t%lld\n", s->seqName, (long long)curBeg, (long long)curEnd) ;
    }

  salsaDestroy (salsas, (int)nSeq) ;
  return 0 ;
}

/****************** main dispatcher ******************/

int main (int argc, char *argv[])
{
  storeCommandLine (argc, argv) ;
  --argc ; ++argv ;

  if (argc < 1)
    { fprintf (stderr,
	       "Usage: taco <subcommand> [options] [args...]\n"
	       "\n"
	       "Subcommands:\n"
	       "  compress     .1aln + FASTA -> .1taco + .1seq (tandem compression)\n"
	       "  collapse     TSV spec + FASTA -> .1taco + .1seq (generic collapse)\n"
	       "  reconstruct  .1taco + .1seq -> FASTA (inverse / undo)\n"
	       "  chain        flatten multiple .1taco + .1seq layers into one\n"
	       "  merge        merge non-overlapping .1taco + .1seq on same original\n"
	       "  lift         coordinate liftover through .1taco file(s)\n"
	       "  info         print summary of .1taco file\n"
		       "  windows      .1taco -> scan windows around scars (feeds FasTAN -r)\n") ;
      return 1 ;
    }

  if      (!strcmp(argv[0], "compress"))    return tacoCompress    (argc-1, argv+1) ;
  else if (!strcmp(argv[0], "collapse"))    return tacoCollapse    (argc-1, argv+1) ;
  else if (!strcmp(argv[0], "reconstruct")) return tacoReconstruct (argc-1, argv+1) ;
  else if (!strcmp(argv[0], "chain"))       return tacoChain       (argc-1, argv+1) ;
  else if (!strcmp(argv[0], "merge"))       return tacoMergeCmd    (argc-1, argv+1) ;
  else if (!strcmp(argv[0], "lift"))        return tacoLift        (argc-1, argv+1) ;
  else if (!strcmp(argv[0], "info"))        return tacoInfo        (argc-1, argv+1) ;
  else if (!strcmp(argv[0], "windows"))     return tacoWindows     (argc-1, argv+1) ;
  else
    { fprintf (stderr, "Unknown subcommand: %s\n", argv[0]) ;
      return 1 ;
    }
}

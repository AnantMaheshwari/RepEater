/*  File: salsa.c
 *  Author: Anant Maheshwari (am3320@cam.ac.uk)
 *-------------------------------------------------------------------
 * Description: Sorted Array for Liftover of Sequence Addresses
 *   Coordinate map between original and taco-compressed sequences.
 *
 * HISTORY:
 * Created: 2026 (am3320)
 *-------------------------------------------------------------------
 */

#include "salsa.h"
#include <string.h>

/****************** internal binary search helpers ******************/

/* Return index of last entry with entry[i].orig_start <= pos, or -1 if none. */
static int bsearchOrigStart (SalsaEntry *entry, int n, I64 pos)
{
  int lo = 0, hi = n - 1, result = -1 ;
  while (lo <= hi)
    { int mid = (lo + hi) / 2 ;
      if (entry[mid].orig_start <= pos) { result = mid ; lo = mid + 1 ; }
      else hi = mid - 1 ;
    }
  return result ;
}

/* Return index of last entry with entry[i].taco_pos <= pos, or -1 if none. */
static int bsearchTacoPos (SalsaEntry *entry, int n, I64 pos)
{
  int lo = 0, hi = n - 1, result = -1 ;
  while (lo <= hi)
    { int mid = (lo + hi) / 2 ;
      if (entry[mid].taco_pos <= pos) { result = mid ; lo = mid + 1 ; }
      else hi = mid - 1 ;
    }
  return result ;
}

/****************** build ******************/

void salsaFinalize (Salsa *salsa, Array entryArr, I64 origLen, I64 tacoLen)
{
  int i ;
  salsa->n       = arrayMax (entryArr) ;
  salsa->origLen = origLen ;
  salsa->tacoLen = tacoLen ;
  salsa->tacoSeq = 0 ;
  salsa->seqName = 0 ;
  salsa->entry = salsa->n ? new (salsa->n, SalsaEntry) : 0 ;
  if (salsa->n)
    memcpy (salsa->entry, arrp(entryArr, 0, SalsaEntry), salsa->n * sizeof(SalsaEntry)) ;
  salsa->cumShift = new (salsa->n + 1, I64) ;
  salsa->cumShift[0] = 0 ;
  for (i = 0 ; i < salsa->n ; ++i)
    salsa->cumShift[i+1] = salsa->cumShift[i]
      + (salsa->entry[i].orig_end - salsa->entry[i].orig_start - salsa->entry[i].unit_len) ;
}

void salsaDestroy (Salsa *salsas, int nSeq)
{
  int i, j ;
  for (i = 0 ; i < nSeq ; ++i)
    { if (salsas[i].entry)
	{ for (j = 0 ; j < salsas[i].n ; ++j)
	    if (salsas[i].entry[j].removed) free (salsas[i].entry[j].removed) ;
	  newFree (salsas[i].entry, salsas[i].n, SalsaEntry) ;
	}
      if (salsas[i].cumShift) newFree (salsas[i].cumShift, salsas[i].n + 1, I64) ;
      if (salsas[i].tacoSeq)  free (salsas[i].tacoSeq) ;
      if (salsas[i].seqName)  free (salsas[i].seqName) ;
    }
  newFree (salsas, nSeq, Salsa) ;
}

/****************** I/O ******************/

void salsaWrite (OneFile *of, Salsa *salsas, int nSeq)
{
  int i, j ;
  for (i = 0 ; i < nSeq ; ++i)
    { Salsa *s = &salsas[i] ;
      oneInt(of,0) = s->origLen ;
      oneInt(of,1) = s->tacoLen ;
      oneWriteLine (of, 'c', 0, 0) ;
      if (s->seqName)
	oneWriteLine (of, 'I', strlen(s->seqName), s->seqName) ;
      for (j = 0 ; j < s->n ; ++j)
        { SalsaEntry *e = &s->entry[j] ;
          oneInt(of,0) = e->orig_start ;
          oneInt(of,1) = e->orig_end ;
          oneInt(of,2) = e->taco_pos ;
          oneInt(of,3) = e->unit_len ;
          oneWriteLine (of, 'L', 0, 0) ;
          I64 removedLen = e->orig_end - e->orig_start - e->unit_len ;
          if (e->removed && removedLen > 0)
            oneWriteLine (of, 'S', removedLen, e->removed) ;
        }
    }
}

Salsa *salsaRead (OneFile *of, int nSeq)
{
  Salsa *salsas = new0 (nSeq, Salsa) ;
  Array entryArr = arrayCreate (64, SalsaEntry) ;
  int seq = -1, i, j ;

  while (oneReadLine (of))
    switch (of->lineType)
      {
      case 'c':
        if (seq >= 0)
          { salsas[seq].n = arrayMax (entryArr) ;
            salsas[seq].entry = salsas[seq].n ? new (salsas[seq].n, SalsaEntry) : 0 ;
            if (salsas[seq].n)
              memcpy (salsas[seq].entry, arrp(entryArr, 0, SalsaEntry),
                      salsas[seq].n * sizeof(SalsaEntry)) ;
          }
        if (++seq >= nSeq) die ("too many sequence lines in .1taco file") ;
        salsas[seq].origLen = oneInt(of,0) ;
        salsas[seq].tacoLen = oneInt(of,1) ;
        salsas[seq].tacoSeq = 0 ;
        salsas[seq].seqName = 0 ;
        arrayMax(entryArr) = 0 ;
        break ;
      case 'I':
	if (seq >= 0)
	  { char *s = oneString(of) ;
	    salsas[seq].seqName = new (strlen(s) + 1, char) ;
	    strcpy (salsas[seq].seqName, s) ;
	  }
	break ;
      case 'L':
        { SalsaEntry *e = arrayp (entryArr, arrayMax(entryArr), SalsaEntry) ;
          e->orig_start = oneInt(of,0) ;
          e->orig_end   = oneInt(of,1) ;
          e->taco_pos   = oneInt(of,2) ;
          e->unit_len   = oneInt(of,3) ;
          e->removed    = 0 ;
        }
        break ;
      case 'S':
        /* removed DNA for the immediately preceding L event */
        if (arrayMax(entryArr) > 0)
          { SalsaEntry *e = arrp (entryArr, arrayMax(entryArr) - 1, SalsaEntry) ;
            I64 len = oneLen(of) ;
            e->removed = new (len + 1, char) ;
            memcpy (e->removed, oneDNAchar(of), len) ;
            e->removed[len] = 0 ;
          }
        break ;
      }

  /* finalise last sequence */
  if (seq >= 0)
    { salsas[seq].n = arrayMax (entryArr) ;
      salsas[seq].entry = salsas[seq].n ? new (salsas[seq].n, SalsaEntry) : 0 ;
      if (salsas[seq].n)
        memcpy (salsas[seq].entry, arrp(entryArr, 0, SalsaEntry),
                salsas[seq].n * sizeof(SalsaEntry)) ;
    }

  /* compute cumShift for all sequences */
  for (i = 0 ; i <= seq ; ++i)
    { Salsa *s = &salsas[i] ;
      s->cumShift = new (s->n + 1, I64) ;
      s->cumShift[0] = 0 ;
      for (j = 0 ; j < s->n ; ++j)
        s->cumShift[j+1] = s->cumShift[j]
          + (s->entry[j].orig_end - s->entry[j].orig_start - s->entry[j].unit_len) ;
    }

  arrayDestroy (entryArr) ;
  return salsas ;
}

/****************** I/O with companion .1seq ******************/

void salsaWriteWithSeq (OneFile *ofTaco, OneFile *ofSeq, Salsa *salsas, int nSeq)
{
  /* coordinate map + removed DNA now all live in .1taco via salsaWrite */
  salsaWrite (ofTaco, salsas, nSeq) ;
  /* .1seq holds only the compressed sequence: one S line per sequence */
  int i ;
  for (i = 0 ; i < nSeq ; ++i)
    oneWriteLine (ofSeq, 'S', salsas[i].tacoLen, salsas[i].tacoSeq) ;
}

Salsa *salsaReadWithSeq (OneFile *ofTaco, OneFile *ofSeq, int nSeq)
{
  /* read coordinate map + removed DNA from .1taco */
  Salsa *salsas = salsaRead (ofTaco, nSeq) ;

  /* read one compressed sequence per entry from .1seq.
     Two formats are accepted:
       (a) scaffold form: s + (S|n)* per scaffold; gaps reconstructed as runs of 'N'
       (b) legacy single-S form: optional I + one S per sequence */
  int seq = 0 ;
  if (!oneReadLine (ofSeq))
    { if (nSeq == 0) return salsas ;
      die ("salsaReadWithSeq: empty .1seq, expected %d sequences", nSeq) ;
    }

  while (seq < nSeq)
    { if (ofSeq->lineType == 's')
        { I64 scafLen = oneInt(ofSeq, 0) ;
          salsas[seq].tacoSeq = new (scafLen + 1, char) ;
          I64 pos = 0 ;
          bool atNextScaf = false ;
          while (oneReadLine (ofSeq))
            { if (ofSeq->lineType == 'S')
                { I64 sl = oneLen(ofSeq) ;
                  if (pos + sl > scafLen) die ("scaffold %d overflow", seq) ;
                  memcpy (salsas[seq].tacoSeq + pos, oneDNAchar(ofSeq), sl) ;
                  pos += sl ;
                }
              else if (ofSeq->lineType == 'n')
                { I64 nl = oneInt(ofSeq, 1) ;
                  if (pos + nl > scafLen) die ("scaffold %d gap overflow", seq) ;
                  memset (salsas[seq].tacoSeq + pos, 'N', nl) ;
                  pos += nl ;
                }
              else if (ofSeq->lineType == 's')
                { atNextScaf = true ; break ; }
            }
          salsas[seq].tacoSeq[pos] = 0 ;
          if (pos != scafLen)
            die ("salsaReadWithSeq: scaffold %d length mismatch: declared %lld, got %lld",
                 seq, (long long)scafLen, (long long)pos) ;
          ++seq ;
          if (!atNextScaf && seq < nSeq)
            die ("salsaReadWithSeq: expected %d sequences in .1seq, got %d", nSeq, seq) ;
        }
      else if (ofSeq->lineType == 'S')
        { I64 len = oneLen(ofSeq) ;
          salsas[seq].tacoSeq = new (len + 1, char) ;
          memcpy (salsas[seq].tacoSeq, oneDNAchar(ofSeq), len) ;
          salsas[seq].tacoSeq[len] = 0 ;
          ++seq ;
          if (seq < nSeq && !oneReadLine (ofSeq))
            die ("salsaReadWithSeq: expected %d sequences in .1seq, got %d", nSeq, seq) ;
        }
      else
        { /* skip unrelated lines (e.g. legacy I lines) */
          if (!oneReadLine (ofSeq))
            die ("salsaReadWithSeq: expected %d sequences in .1seq, got %d", nSeq, seq) ;
        }
    }
  return salsas ;
}

/****************** reconstruction ******************/

char *salsaReconstruct (Salsa *salsa)
{
  char *result = new (salsa->origLen + 1, char) ;
  I64 tacoPos = 0, origPos = 0 ;
  int i ;

  for (i = 0 ; i < salsa->n ; ++i)
    { SalsaEntry *e = &salsa->entry[i] ;
      /* copy 1:1 gap region from tacoSeq */
      I64 gap = e->taco_pos - tacoPos ;
      if (gap > 0) { memcpy (result + origPos, salsa->tacoSeq + tacoPos, gap) ; origPos += gap ; tacoPos += gap ; }
      /* copy the kept unit from tacoSeq */
      if (e->unit_len > 0)
	{ memcpy (result + origPos, salsa->tacoSeq + tacoPos, e->unit_len) ;
	  origPos += e->unit_len ; tacoPos += e->unit_len ;
	}
      /* insert the removed DNA */
      I64 removedLen = e->orig_end - e->orig_start - e->unit_len ;
      if (removedLen > 0 && e->removed)
	{ memcpy (result + origPos, e->removed, removedLen) ; origPos += removedLen ; }
    }
  /* copy trailing 1:1 region */
  I64 tail = salsa->tacoLen - tacoPos ;
  if (tail > 0) memcpy (result + origPos, salsa->tacoSeq + tacoPos, tail) ;
  origPos += tail ;
  result[origPos] = 0 ;
  return result ;
}

/****************** merge and compose ******************/

/* helper: recompute cumShift and taco_pos for a sorted event array */
static void recomputeShifts (Salsa *salsa)
{
  int i ;
  salsa->cumShift = new (salsa->n + 1, I64) ;
  salsa->cumShift[0] = 0 ;
  I64 shift = 0 ;
  for (i = 0 ; i < salsa->n ; ++i)
    { salsa->entry[i].taco_pos = salsa->entry[i].orig_start - shift ;
      I64 removed = salsa->entry[i].orig_end - salsa->entry[i].orig_start - salsa->entry[i].unit_len ;
      shift += removed ;
      salsa->cumShift[i+1] = shift ;
    }
  salsa->tacoLen = salsa->origLen - shift ;
}

/* helper: build tacoSeq from original sequence and event list */
static char *buildTacoSeq (char *origSeq, I64 origLen, SalsaEntry *entry, int n, I64 tacoLen)
{
  char *taco = new (tacoLen + 1, char) ;
  I64 nOld = 0, nNew = 0 ;
  int i ;
  for (i = 0 ; i < n ; ++i)
    { SalsaEntry *e = &entry[i] ;
      I64 gap = e->orig_start - nOld ;
      if (gap > 0) { memcpy (taco + nNew, origSeq + nOld, gap) ; nNew += gap ; }
      if (e->unit_len > 0) { memcpy (taco + nNew, origSeq + e->orig_start, e->unit_len) ; nNew += e->unit_len ; }
      nOld = e->orig_end ;
    }
  I64 tail = origLen - nOld ;
  if (tail > 0) { memcpy (taco + nNew, origSeq + nOld, tail) ; nNew += tail ; }
  taco[nNew] = 0 ;
  return taco ;
}

static int entryCompare (const void *a, const void *b)
{ SalsaEntry *ea = (SalsaEntry*)a, *eb = (SalsaEntry*)b ;
  if (ea->orig_start != eb->orig_start) return (ea->orig_start < eb->orig_start) ? -1 : 1 ;
  return (ea->orig_end < eb->orig_end) ? -1 : (ea->orig_end > eb->orig_end) ? 1 : 0 ;
}

Salsa salsaMerge (Salsa *a, Salsa *b)
{
  Salsa result ;
  memset (&result, 0, sizeof(Salsa)) ;
  result.origLen = a->origLen ;
  result.n = a->n + b->n ;
  result.entry = result.n ? new (result.n, SalsaEntry) : 0 ;

  int i ;
  if (a->n) memcpy (result.entry, a->entry, a->n * sizeof(SalsaEntry)) ;
  if (b->n) memcpy (result.entry + a->n, b->entry, b->n * sizeof(SalsaEntry)) ;

  for (i = 0 ; i < a->n ; ++i)
    if (a->entry[i].removed)
      { I64 len = a->entry[i].orig_end - a->entry[i].orig_start - a->entry[i].unit_len ;
	result.entry[i].removed = new (len + 1, char) ;
	memcpy (result.entry[i].removed, a->entry[i].removed, len) ;
	result.entry[i].removed[len] = 0 ;
      }
  for (i = 0 ; i < b->n ; ++i)
    if (b->entry[i].removed)
      { I64 len = b->entry[i].orig_end - b->entry[i].orig_start - b->entry[i].unit_len ;
	result.entry[a->n + i].removed = new (len + 1, char) ;
	memcpy (result.entry[a->n + i].removed, b->entry[i].removed, len) ;
	result.entry[a->n + i].removed[len] = 0 ;
      }

  /* sort by orig_start */
  qsort (result.entry, result.n, sizeof(SalsaEntry), entryCompare) ;

  /* verify no overlaps */
  for (i = 1 ; i < result.n ; ++i)
    if (result.entry[i].orig_start < result.entry[i-1].orig_end)
      die ("salsaMerge: overlapping events at orig_start %lld and orig_end %lld",
	   (long long)result.entry[i].orig_start, (long long)result.entry[i-1].orig_end) ;

  /* recompute taco_pos and cumShift */
  recomputeShifts (&result) ;

  /* build tacoSeq if both inputs have it */
  if (a->tacoSeq && b->tacoSeq)
    { char *orig = salsaReconstruct (a) ; /* reconstruct from a (same original as b) */
      result.tacoSeq = buildTacoSeq (orig, result.origLen, result.entry, result.n, result.tacoLen) ;
      free (orig) ;
    }

  /* copy seqName from a */
  if (a->seqName)
    { result.seqName = new (strlen(a->seqName) + 1, char) ;
      strcpy (result.seqName, a->seqName) ;
    }

  return result ;
}

Salsa salsaCompose (Salsa *inner, Salsa *outer)
{
  /* inner: orig -> A-space.  outer: A-space -> B-space.  result: orig -> B-space. */
  Salsa result ;
  memset (&result, 0, sizeof(Salsa)) ;
  result.origLen = inner->origLen ;

  /* Step 1: reconstruct the original sequence from inner */
  char *origSeq = salsaReconstruct (inner) ;

  /* Step 2: lift each outer event back to original-space and collect all events */
  Array events = arrayCreate (inner->n + outer->n, SalsaEntry) ;

  
  /* An inner event is absorbed if an outer event spans across it in A-space. */
  bool *innerAbsorbed = new0 (inner->n, bool) ;
  int lo = 0 ; /* first inner event that might still be absorbed */

  int j ;
  for (j = 0 ; j < outer->n ; ++j)
    { SalsaEntry *o = &outer->entry[j] ;
      /* lift outer event interval [orig_start, orig_end) from A-space to original-space */
      I64 origA, origB ;
      salsaTacoToOrigInterval (inner, o->orig_start, o->orig_end, &origA, &origB) ;

      /* advance lo past inner events that end before this window */
      while (lo < inner->n && inner->entry[lo].orig_end <= origA) ++lo ;

      /* mark inner events absorbed by this outer event */
      int k ;
      for (k = lo ; k < inner->n && inner->entry[k].orig_start < origB ; ++k)
        if (!innerAbsorbed[k] && inner->entry[k].orig_end <= origB)
          innerAbsorbed[k] = true ;

      /* build the removed DNA for this outer event in original-space */
      I64 removedLen = (origB - origA) - o->unit_len ;
      char *removed = 0 ;
      if (removedLen > 0)
        { removed = new (removedLen + 1, char) ;
          memcpy (removed, origSeq + origA + o->unit_len, removedLen) ;
          removed[removedLen] = 0 ;
        }

      SalsaEntry *e = arrayp (events, arrayMax(events), SalsaEntry) ;
      e->orig_start = origA ;
      e->orig_end   = origB ;
      e->unit_len   = o->unit_len ;
      e->taco_pos   = 0 ; /* will be recomputed */
      e->removed    = removed ;
    }

  /* Add non-absorbed inner events */
  int i ;
  for (i = 0 ; i < inner->n ; ++i)
    if (!innerAbsorbed[i])
      { SalsaEntry *src = &inner->entry[i] ;
        SalsaEntry *e = arrayp (events, arrayMax(events), SalsaEntry) ;
        e->orig_start = src->orig_start ;
        e->orig_end   = src->orig_end ;
        e->unit_len   = src->unit_len ;
        e->taco_pos   = 0 ;
        I64 removedLen = src->orig_end - src->orig_start - src->unit_len ;
        if (removedLen > 0 && src->removed)
          { e->removed = new (removedLen + 1, char) ;
            memcpy (e->removed, src->removed, removedLen) ;
            e->removed[removedLen] = 0 ;
          }
        else
          e->removed = 0 ;
      }
  newFree (innerAbsorbed, inner->n, bool) ;

  /* sort by orig_start */
  result.n = arrayMax (events) ;
  result.entry = result.n ? new (result.n, SalsaEntry) : 0 ;
  if (result.n)
    { memcpy (result.entry, arrp(events, 0, SalsaEntry), result.n * sizeof(SalsaEntry)) ;
      qsort (result.entry, result.n, sizeof(SalsaEntry), entryCompare) ;
    }
  arrayDestroy (events) ;

  /* recompute taco_pos and cumShift */
  recomputeShifts (&result) ;

  /* build the flat tacoSeq */
  result.tacoSeq = buildTacoSeq (origSeq, result.origLen, result.entry, result.n, result.tacoLen) ;
  free (origSeq) ;

  /* copy seqName from inner */
  if (inner->seqName)
    { result.seqName = new (strlen(inner->seqName) + 1, char) ;
      strcpy (result.seqName, inner->seqName) ;
    }

  return result ;
}

/****************** coordinate transforms ******************/

/*
 * Invariant exploited throughout:
 *   entry[i].taco_pos == entry[i].orig_start - cumShift[i]
 */

I64 salsaOrigToTacoPos (Salsa *salsa, I64 pos)
{
  int i = bsearchOrigStart (salsa->entry, salsa->n, pos) ;
  if (i < 0)
    return pos ;                          /* before all events: identity */
  if (pos < salsa->entry[i].orig_end)    /* inside event i */
    { if (salsa->entry[i].unit_len == 0)   /* full deletion: map to boundary */
	return salsa->entry[i].taco_pos ;
      return salsa->entry[i].taco_pos
	+ (pos - salsa->entry[i].orig_start) % salsa->entry[i].unit_len ;
    }
  return pos - salsa->cumShift[i+1] ;    /* in 1:1 gap after event i */
}

I64 salsaTacoToOrigPos (Salsa *salsa, I64 pos)
{
  int i = bsearchTacoPos (salsa->entry, salsa->n, pos) ;
  if (i < 0)
    return pos ;
  if (pos < salsa->entry[i].taco_pos + salsa->entry[i].unit_len) /* inside event i */
    return salsa->entry[i].orig_start + (pos - salsa->entry[i].taco_pos) ;
  return pos + salsa->cumShift[i+1] ;
}

/* Interval [a,b): taco -> orig  (expanding) */
void salsaTacoToOrigInterval (Salsa *salsa, I64 a, I64 b, I64 *oa, I64 *ob)
{
  int i, j ;

  /* left endpoint */
  i = bsearchTacoPos (salsa->entry, salsa->n, a) ;
  if (i < 0)
    *oa = a ;
  else if (a < salsa->entry[i].taco_pos + salsa->entry[i].unit_len)
    *oa = salsa->entry[i].orig_start ;
  else
    *oa = a + salsa->cumShift[i+1] ;

  /* right endpoint: operate on b-1 (last position inside [a,b)) */
  j = (b > 0) ? bsearchTacoPos (salsa->entry, salsa->n, b - 1) : -1 ;
  if (j < 0)
    *ob = b ;
  else if (b - 1 < salsa->entry[j].taco_pos + salsa->entry[j].unit_len)
    *ob = salsa->entry[j].orig_end ;
  else
    *ob = b + salsa->cumShift[j+1] ;
}

/* Interval [a,b): orig -> taco  (contracting) */
void salsaOrigToTacoInterval (Salsa *salsa, I64 a, I64 b, I64 *ta, I64 *tb)
{
  int i, j ;

  /* left endpoint */
  i = bsearchOrigStart (salsa->entry, salsa->n, a) ;
  if (i < 0)
    *ta = a ;
  else if (a < salsa->entry[i].orig_end)
    *ta = salsa->entry[i].taco_pos ;
  else
    *ta = a - salsa->cumShift[i+1] ;

  /* right endpoint: operate on b-1 */
  j = (b > 0) ? bsearchOrigStart (salsa->entry, salsa->n, b - 1) : -1 ;
  if (j < 0)
    *tb = b ;
  else if (b - 1 < salsa->entry[j].orig_end)
    *tb = salsa->entry[j].taco_pos + salsa->entry[j].unit_len ;
  else
    *tb = b - salsa->cumShift[j+1] ;
}

/****************** SalsaStack: multi-layer composition ******************/

SalsaStack *salsaStackCreate (int nLayer)
{
  SalsaStack *stack = new0 (1, SalsaStack) ;
  stack->nLayer = nLayer ;
  stack->nSeq   = new0 (nLayer, int) ;
  stack->layer  = new0 (nLayer, Salsa*) ;
  return stack ;
}

void salsaStackSetLayer (SalsaStack *stack, int k, Salsa *salsas, int nSeq)
{
  stack->layer[k] = salsas ;
  stack->nSeq[k]  = nSeq ;
}

SalsaStack *salsaStackLoad (char **filenames, int nFiles)
{
  SalsaStack *stack = salsaStackCreate (nFiles) ;
  int k ;
  for (k = 0 ; k < nFiles ; ++k)
    { OneSchema *schema = oneSchemaCreateFromText (salsaSchemaText) ;
      OneFile *of = oneFileOpenRead (filenames[k], schema, "taco", 1) ;
      oneSchemaDestroy (schema) ;
      if (!of) die ("salsaStackLoad: failed to open %s", filenames[k]) ;
      I64 nSeq = 0 ;
      oneStats (of, 'c', &nSeq, 0, 0) ;
      Salsa *salsas = salsaRead (of, (int)nSeq) ;
      oneFileClose (of) ;
      salsaStackSetLayer (stack, k, salsas, (int)nSeq) ;
    }
  return stack ;
}

void salsaStackDestroy (SalsaStack *stack)
{
  int k ;
  for (k = 0 ; k < stack->nLayer ; ++k)
    salsaDestroy (stack->layer[k], stack->nSeq[k]) ;
  newFree (stack->nSeq, stack->nLayer, int) ;
  newFree (stack->layer, stack->nLayer, Salsa*) ;
  newFree (stack, 1, SalsaStack) ;
}

I64 salsaStackFinalToOrigPos (SalsaStack *stack, int seq, I64 pos)
{
  int k ;
  for (k = stack->nLayer - 1 ; k >= 0 ; --k)
    pos = salsaTacoToOrigPos (&stack->layer[k][seq], pos) ;
  return pos ;
}

I64 salsaStackOrigToFinalPos (SalsaStack *stack, int seq, I64 pos)
{
  int k ;
  for (k = 0 ; k < stack->nLayer ; ++k)
    pos = salsaOrigToTacoPos (&stack->layer[k][seq], pos) ;
  return pos ;
}

void salsaStackFinalToOrigInterval (SalsaStack *stack, int seq, I64 a, I64 b, I64 *oa, I64 *ob)
{
  int k ;
  for (k = stack->nLayer - 1 ; k >= 0 ; --k)
    salsaTacoToOrigInterval (&stack->layer[k][seq], a, b, &a, &b) ;
  *oa = a ; *ob = b ;
}

void salsaStackOrigToFinalInterval (SalsaStack *stack, int seq, I64 a, I64 b, I64 *fa, I64 *fb)
{
  int k ;
  for (k = 0 ; k < stack->nLayer ; ++k)
    salsaOrigToTacoInterval (&stack->layer[k][seq], a, b, &a, &b) ;
  *fa = a ; *fb = b ;
}

/************************ end of file *************************/

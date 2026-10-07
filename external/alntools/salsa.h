/*  File: salsa.h
 *  Author: Anant Maheshwari (am3320@cam.ac.uk)
 *-------------------------------------------------------------------
 * Description: Sorted Array for Liftover of Sequence Addresses
 *   Coordinate map between original and taco-compressed sequences.
 *   Stores one Salsa per sequence.  Each Salsa holds a sorted array of
 *   SalsaEntry records, one per tandem-repeat compression event, plus a
 *   cumulative-shift array for O(log n) coordinate transforms.
 *
 *   Coordinate conventions: all positions are 0-based; intervals are half-open [a,b).
 *   Output files use the .1taco OneCode file type.
 *
 * HISTORY:
 * Created: 2026 (am3320)
 *-------------------------------------------------------------------
 */

#ifndef SALSA_DEFINED
#define SALSA_DEFINED

#include "ONElib.h"
#include "utils.h"
#include "array.h"

/****************** schema for the .1taco file type ******************/

static char salsaSchemaText[] =
  "1 3 def 2 1\n"
  ".\n"
  "P 4 taco                    TACO MAP\n"
  "O c 2 3 INT 3 INT           sequence: origLen tacoLen\n"
  "D I 1 6 STRING              sequence identifier\n"
  "G L                         sequence (c) groups lift entries (L)\n"
  "O L 4 3 INT 3 INT 3 INT 3 INT  lift event: orig_start orig_end taco_pos unit_len\n"
  "D S 1 3 DNA                 removed DNA for this event (len = orig_end-orig_start-unit_len)\n"
;

/****************** data structures ******************/

typedef struct {
  I64 orig_start ;  // start of tandem repeat in original coords (0-based, inclusive)
  I64 orig_end ;    // end of tandem repeat in original coords (exclusive)
  I64 taco_pos ;    // start of the compressed unit in taco coords
  I32 unit_len ;    // length of one repeat unit (the single copy kept)
  char *removed ;   // removed DNA (NULL if not loaded); len = orig_end - orig_start - unit_len
} SalsaEntry ;

typedef struct {
  int        n ;        // number of compression events
  SalsaEntry *entry ;   // array of n events, sorted by orig_start (== sorted by taco_pos)
  I64        *cumShift ; // cumShift[i] = total bases removed by events 0..i-1;  size n+1
  I64         origLen ; // original sequence length
  I64         tacoLen ; // compressed sequence length
  char       *tacoSeq ; // compressed sequence (NULL if not loaded); len = tacoLen
  char       *seqName ; // sequence identifier (NULL if not loaded)
} Salsa ;

/****************** function declarations ******************/

/* Build --------------------------------------------------------------------- */

// Finalise one sequence's Salsa after its compression loop.
// entryArr is a temporary Array of SalsaEntry built during the loop (see taco.c).
void salsaFinalize (Salsa *salsa, Array entryArr, I64 origLen, I64 tacoLen) ;

void salsaDestroy (Salsa *salsas, int nSeq) ;

/* I/O ----------------------------------------------------------------------- */

// Write all nSeq salsas to an open .1taco OneFile (opened for "taco").
void salsaWrite (OneFile *of, Salsa *salsas, int nSeq) ;

// Read from an open .1taco OneFile.  nSeq must match the number of 's' lines.
Salsa *salsaRead (OneFile *of, int nSeq) ;

/* I/O with companion .1seq DNA file ----------------------------------------- */

// Write nSeq salsas to .1taco + companion .1seq (compressed seqs + removed DNA).
// Requires salsa[i].tacoSeq and entry[j].removed to be set.
void salsaWriteWithSeq (OneFile *ofTaco, OneFile *ofSeq, Salsa *salsas, int nSeq) ;

// Read paired .1taco + .1seq files.  Populates tacoSeq and entry[j].removed.
Salsa *salsaReadWithSeq (OneFile *ofTaco, OneFile *ofSeq, int nSeq) ;

/* Reconstruction ------------------------------------------------------------ */

// Reconstruct original sequence from a Salsa with tacoSeq and removed DNA.
// Returns newly allocated char[origLen].
char *salsaReconstruct (Salsa *salsa) ;

/* Composition --------------------------------------------------------------- */

// Flatten two layers: inner (orig->A) and outer (A->B) into one (orig->B).
// Full general case: handles overlapping events across layers.
// Both must have tacoSeq and removed DNA loaded.
Salsa salsaCompose (Salsa *inner, Salsa *outer) ;

// Merge non-overlapping events on the same original sequence.
// Dies if events overlap.  Both must have tacoSeq and removed DNA loaded.
Salsa salsaMerge (Salsa *a, Salsa *b) ;

/* Coordinate transforms ----------------------------------------------------- */

// Single-position: orig -> taco.
//   Positions inside a repeat map via modulo into the kept copy.
I64 salsaOrigToTacoPos (Salsa *salsa, I64 pos) ;

// Single-position: taco -> orig.
//   Positions inside a compressed unit map to the corresponding offset in the first copy.
I64 salsaTacoToOrigPos (Salsa *salsa, I64 pos) ;

// Half-open interval [a,b): taco -> orig  (expanding).
//   If an endpoint falls inside a compressed unit it is pushed out to orig_start / orig_end
//   so the result spans all original sequence covered by the taco interval.
void salsaTacoToOrigInterval (Salsa *salsa, I64 a, I64 b, I64 *oa, I64 *ob) ;

// Half-open interval [a,b): orig -> taco  (contracting).
//   If an endpoint falls inside a repeat it pulls in to taco_pos / taco_pos+unit_len.
void salsaOrigToTacoInterval (Salsa *salsa, I64 a, I64 b, I64 *ta, I64 *tb) ;

/* SalsaStack: multi-layer composition for iterative collapse rounds ---------- */

typedef struct {
  int     nLayer ;     // number of .1taco layers
  int    *nSeq ;       // nSeq[k] = sequences in layer k
  Salsa **layer ;      // layer[k] = Salsa array for round k
} SalsaStack ;

SalsaStack *salsaStackCreate  (int nLayer) ;
void        salsaStackSetLayer (SalsaStack *stack, int k, Salsa *salsas, int nSeq) ;
SalsaStack *salsaStackLoad    (char **filenames, int nFiles) ;
void        salsaStackDestroy (SalsaStack *stack) ;

I64  salsaStackFinalToOrigPos      (SalsaStack *stack, int seq, I64 pos) ;
I64  salsaStackOrigToFinalPos      (SalsaStack *stack, int seq, I64 pos) ;
void salsaStackFinalToOrigInterval (SalsaStack *stack, int seq, I64 a, I64 b, I64 *oa, I64 *ob) ;
void salsaStackOrigToFinalInterval (SalsaStack *stack, int seq, I64 a, I64 b, I64 *fa, I64 *fb) ;

#endif /* SALSA_DEFINED */

/************************ end of file *************************/

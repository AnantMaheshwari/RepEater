/*  File: alntools.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 16 10:18 2025 (rd109)
 * Created: Thu Oct 16 02:49:56 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "dict.h"
#include "utils.h"
#include "ONElib.h"

#define VERSION "0.1"

static char schemaText[] =
  "1 3 def 2 1                 schema for aln and FastGA\n"
  ".\n"
  ".\n"
  "P 3 gdb                             GDB\n"
  "D f 4 4 REAL 4 REAL 4 REAL 4 REAL   global: base frequency vector\n"
  "D u 0                               global: upper case when displayed\n"
  "O S 1 6 STRING                      id for a scaffold\n"
  "D G 1 3 INT                         gap of given length\n"
  "D C 1 3 INT                         contig of given length\n"
  "D M 1 8 INT_LIST                    mask pair list for a contig\n"
  ".\n"
  "P 3 seq                     SEQUENCE\n"
  "O s 2 3 INT 6 STRING        length and id for group of sequences = a scaffold\n"
  "G S                         scaffolds (s) group sequence objects (S)\n"
  "D n 2 4 CHAR 3 INT          non-acgt chars outside (between) sequences within scaffold\n"
  "O S 1 3 DNA                 sequence\n"
  "D I 1 6 STRING              identifier of sequence\n"
  ".\n"
  "P 3 aln                     ALIGNMENTS\n"
  "D t 1 3 INT                 trace point spacing in a - global\n"
  ".                           GDB skeleton (may not be presend)\n"
  "O g 0                       groups scaffolds into a GDB skeleton\n"
  "G S                         collection of scaffolds constituting a GDB\n"
  "O S 1 6 STRING              id for a scaffold\n"
  "D G 1 3 INT                 gap of given length\n"
  "D C 1 3 INT                 contig of given length\n"
  "D M 1 8 INT_LIST            mask pair list for a contig\n"
  ".\n"
  "O a 0                       groups A's into a colinear chain\n"
  "G A                         chains (a) group alignment objects (A)\n"
  "D p 2 3 INT 3 INT           spacing in a,b between end of previous alignment and start of next\n"
  ".                           alignment: a_read[beg..end] b_read[beg..end], 0-indexing\n"
  "O A 6 3 INT 3 INT 3 INT 3 INT 3 INT 3 INT\n"
  "D L 2 3 INT 3 INT           lengths of sequences a and b\n"
  "D R 0                       flag: reverse-complement sequence b\n"
  "D D 1 3 INT                 differences: number of diffs = substitions + indels\n"
  "D T 1 8 INT_LIST            trace points in b\n"
  "D X 1 8 INT_LIST            number of differences in alignment per trace interval\n"
  "D Q 1 3 INT                 quality: alignment confidence in phred units (currently unused)\n"
  "D E 1 3 INT                 match: number of equal bases (currently unused)\n"
  "D Z 1 6 STRING              cigar string: encodes precise alignment (currently unused)\n"
  "D U 1 3 INT                 putative unit size of a TR alignment (FASTAN)\n"
  "P 4 taco                    TACO MAP\n"
  "O c 2 3 INT 3 INT           sequence: origLen tacoLen\n"
  "D I 1 6 STRING              sequence identifier\n"
  "G L                         sequence (c) groups lift entries (L)\n"
  "O L 4 3 INT 3 INT 3 INT 3 INT  lift event: orig_start orig_end taco_pos unit_len\n"
  ;

typedef struct {
  char  *seqFileName, *seqPathName ;  // from of->reference - not owned by the Gdb
  double fA, fC, fG, fT ;	// 'f' frequences of A,C,G,T
  bool   isUpper ;	 	// 'u' upper case for non-masked sequence
  int    nSeq, nCtg, nGap ;
  I64    maxSeq, maxCtg, maxMask ; // used to allocate arrays below
  I64    totSeq, totCtg, totMask ; // total length including gaps, without gaps, of mask
  DICT  *seqDict ;       	// names of sequences
  I64   *seqLen ;	 	// lengths of sequences
  I64   *ctgLen ;	 	// contig lengths
  int   *ctgSeq ;	 	// parent sequence for each contig
  I64   *ctgPos ;	 	// offset in parent of each contig
  int   *ctgMaskCount ;  	// number of masks in each contig
  int   *ctgMaskStart ;         // start of contigs's mask entries in ->mask
  I64   *mask ;		 	// mask positions (two per masked region)
} Gdb ;

static inline int ctg2seq (Gdb *gdb, int ctg) { return gdb->ctgSeq[ctg] ; }
static inline I64 ctg2pos (Gdb *gdb, int ctg, I64 x) { return gdb->ctgPos[ctg] + x ; }

/****************** BED file structures ***********************/

typedef struct {
  U32 seq ;          // for bed output need seq
  int ctg ;          // for masking need contig
  I64 start, end ;
  int unit ;
  int score ;
} TanLine ;

static int tanSort (const void *a, const void *b)
{
  TanLine *ta = (TanLine*)a, *tb = (TanLine*)b ;
  if (ta->seq != tb->seq) return ta->seq - tb->seq ;
  if (ta->start != tb->start) return ta->start - tb->start ;
  return ta->end - tb->end ;
}

/************ in gdb.c ************/

Gdb *readGdb (OneFile *of, int k, FILE *report) ;
// read from either a .1gdb file (set k == 1) or the k'th skeleton within of
// writes a summary line if report != NULL

void writeGdb (OneFile *of, Gdb *gdb, int k, FILE *report) ;

void gdbDestroy (Gdb *gdb) ;

OneFile *gdbFile (OneFile *ofAln, int number) ;

/************************ end of file *************************/

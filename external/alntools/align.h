/*  File: align.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2024
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Aug 10 00:26 2024 (rd109)
 * Created: Fri Aug  9 14:21:12 2024 (rd109)
 *-------------------------------------------------------------------
 */

typedef unsigned char      uint8;
typedef unsigned short     uint16;
typedef unsigned int       uint32;
typedef unsigned long long uint64;
typedef signed char        int8;
typedef signed short       int16;
typedef signed int         int32;
typedef signed long long   int64;
typedef float              float32;
typedef double             float64;

typedef struct
  { void     *trace;
    int       tlen;
    int       diffs;
    int       abpos, bbpos; // inclusive, 0-indexed
    int       aepos, bepos; // exclusive, 0-indexed
  } Path;

#define COMP_FLAG  0x1	// b sequence should be reverse-complemented
                        // Gene had other flags, but none of these are used in .1aln files
#define COMP(x)   ((x) & COMP_FLAG)

typedef struct {
  Path    path;         // Path: begin- and end-point of alignment + diffs
  uint32  flags;        // Pipeline status and complementation flags
  int     aread;        // Id # of A sequence - 0 indexed
  int     bread;        // Id # of B sequence - 0 indexed
} Overlap;

/******** end of file *********/

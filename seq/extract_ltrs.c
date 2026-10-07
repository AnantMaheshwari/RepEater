/* extract_ltrs.c
 *
 * Read a FastLTR .1aln alignment file and a genome FASTA, write an
 * LTR-sequence FASTA with the structural header schema consumed by
 * engine/feature_cluster.py:
 *
 *   >LTR_<n> flank_l=<a-b> repeat1=<a-b> repeat2=<a-b> flank_r=<a-b>
 *
 * Coordinates are scaffold-absolute. Alignments are sorted by
 * (scaffold_idx, r1s, r2s) for deterministic output regardless
 * of .1aln record order.
 *
 * Links against the SEQUENCE_UTILITIES copies of ONElib / seqio / utils /
 * array / dict. See algorithms/build_extract_ltrs.sh.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "array.h"
#include "dict.h"
#include "seqio.h"
#include "ONElib.h"

#define FLANK 100

typedef struct {
  I64 scaffold_idx;
  I64 offset;
} Contig;

typedef struct {
  I64 scaffold_idx;
  I64 r1s, r1e;   /* repeat1 (b copy) */
  I64 r2s, r2e;   /* repeat2 (a copy) */
} Align;

typedef struct {
  char *seq;
  I64   len;
} ScafSeq;

static int cmp_align(const void *a, const void *b) {
  const Align *x = (const Align *)a;
  const Align *y = (const Align *)b;
  if (x->scaffold_idx != y->scaffold_idx)
    return (x->scaffold_idx < y->scaffold_idx) ? -1 : 1;
  if (x->r1s != y->r1s) return (x->r1s < y->r1s) ? -1 : 1;
  if (x->r2s != y->r2s) return (x->r2s < y->r2s) ? -1 : 1;
  return 0;
}

static void usage(void) {
  fprintf(stderr,
          "usage: extract_ltrs <genome.fa[.gz]> <in.1aln> <out.fasta>\n");
  exit(1);
}

int main(int argc, char *argv[]) {
  if (argc != 4) usage();
  const char *genome_path = argv[1];
  const char *aln_path    = argv[2];
  const char *out_path    = argv[3];

  /* ---- parse .1aln --------------------------------------------------- */
  OneFile *of = oneFileOpenRead(aln_path, 0, 0, 1);
  if (!of) die("failed to open %s", aln_path);

  DICT  *scaf_dict = dictCreate(4096);
  Array  contigs   = arrayCreate(8192, Contig);
  Array  aligns    = arrayCreate(8192, Align);

  I64 cur_scaf    = -1;
  I64 scaf_offset = 0;
  Align cur = {0};
  bool have_align = false;

  while (oneReadLine(of)) {
    char t = of->lineType;
    if (t == 'S') {
      char *name = oneString(of);
      U64 k;
      dictAdd(scaf_dict, name, &k);
      cur_scaf    = (I64)k;
      scaf_offset = 0;
    } else if (t == 'C') {
      I64 clen = oneInt(of, 0);
      Contig *c = arrayp(contigs, arrayMax(contigs), Contig);
      c->scaffold_idx = cur_scaf;
      c->offset       = scaf_offset;
      scaf_offset    += clen;
    } else if (t == 'G') {
      scaf_offset += oneInt(of, 0);
    } else if (t == 'A') {
      I64 a_contig = oneInt(of, 0);
      I64 a_start  = oneInt(of, 1);
      I64 a_end    = oneInt(of, 2);
      /* b_contig = oneInt(of, 3);  assumed == a_contig for LTRs */
      I64 b_start  = oneInt(of, 4);
      I64 b_end    = oneInt(of, 5);
      if (a_contig < 0 || a_contig >= (I64)arrayMax(contigs))
        die("A line references out-of-range contig %lld", (long long)a_contig);
      Contig *c = arrp(contigs, a_contig, Contig);
      cur.scaffold_idx = c->scaffold_idx;
      cur.r1s = c->offset + b_start;
      cur.r1e = c->offset + b_end;
      cur.r2s = c->offset + a_start;
      cur.r2e = c->offset + a_end;
      have_align = true;
    } else if (t == 'U') {
      if (have_align &&
          (cur.r1e - cur.r1s) != 0 &&
          (cur.r2e - cur.r2s) != 0) {
        array(aligns, arrayMax(aligns), Align) = cur;
      }
      have_align = false;
      memset(&cur, 0, sizeof(cur));
    }
    /* D and any other types ignored */
  }
  oneFileClose(of);

  fprintf(stderr,
          "Parsed %llu scaffolds, %llu contigs, %llu alignments from %s\n",
          (unsigned long long)dictMax(scaf_dict),
          (unsigned long long)arrayMax(contigs),
          (unsigned long long)arrayMax(aligns),
          aln_path);

  /* ---- stream genome, capture scaffolds we care about --------------- */
  U64 n_scaf = dictMax(scaf_dict);
  ScafSeq *scaf_seqs = new0(n_scaf, ScafSeq);

  SeqIO *si = seqIOopenRead((char *)genome_path, dna2textConv, 0);
  if (!si) die("failed to open %s", genome_path);
  while (seqIOread(si)) {
    U64 k;
    bool found = false;
    if (si->descLen > 0) {
      char *fullname = (char *)malloc(si->idLen + 1 + si->descLen + 1);
      if (!fullname) die("out of memory building scaffold name");
      memcpy(fullname, sqioId(si), si->idLen);
      fullname[si->idLen] = ' ';
      memcpy(fullname + si->idLen + 1, sqioDesc(si), si->descLen);
      fullname[si->idLen + 1 + si->descLen] = '\0';
      found = dictFind(scaf_dict, fullname, &k);
      free(fullname);
    }
    if (!found && !dictFind(scaf_dict, sqioId(si), &k)) continue;
    scaf_seqs[k].len = si->seqLen;
    scaf_seqs[k].seq = (char *)malloc(si->seqLen);
    if (!scaf_seqs[k].seq) die("out of memory copying scaffold %s", sqioId(si));
    memcpy(scaf_seqs[k].seq, sqioSeq(si), si->seqLen);
  }
  seqIOclose(si);

  /* ---- sort alignments for deterministic output ---------------------- */
  qsort(arrp(aligns, 0, Align), arrayMax(aligns), sizeof(Align),
        cmp_align);

  /* ---- emit fasta --------------------------------------------------- */
  FILE *out = fopen(out_path, "w");
  if (!out) die("failed to open %s for writing", out_path);

  U64 written = 0, skipped = 0;
  for (U64 i = 0; i < arrayMax(aligns); ++i) {
    Align *a = arrp(aligns, i, Align);
    ScafSeq *s = &scaf_seqs[a->scaffold_idx];
    if (!s->seq) {
      fprintf(stderr, "warn: scaffold %s not found in genome\n",
              dictName(scaf_dict, a->scaffold_idx));
      ++skipped;
      continue;
    }
    I64 ltr_start = a->r1s;
    I64 ltr_end   = a->r2e;
    if (ltr_end > s->len) {
      fprintf(stderr,
              "warn: coord %lld-%lld exceeds scaffold %s len %lld\n",
              (long long)ltr_start, (long long)ltr_end,
              dictName(scaf_dict, a->scaffold_idx), (long long)s->len);
      ++skipped;
      continue;
    }
    I64 es = ltr_start > FLANK ? ltr_start - FLANK : 0;
    I64 ee = ltr_end + FLANK < s->len ? ltr_end + FLANK : s->len;

    ++written;
    fprintf(out,
            ">LTR_%llu scaffold=%s flank_l=%lld-%lld repeat1=%lld-%lld "
            "repeat2=%lld-%lld flank_r=%lld-%lld\n",
            (unsigned long long)written,
            dictName(scaf_dict, a->scaffold_idx),
            (long long)es,       (long long)a->r1s,
            (long long)a->r1s,   (long long)a->r1e,
            (long long)a->r2s,   (long long)a->r2e,
            (long long)a->r2e,   (long long)ee);
    fwrite(s->seq + es, 1, (size_t)(ee - es), out);
    fputc('\n', out);
  }
  fclose(out);

  fprintf(stderr, "Wrote %llu LTR sequences to %s (%llu skipped)\n",
          (unsigned long long)written, out_path,
          (unsigned long long)skipped);

  /* ---- cleanup (cosmetic; process is exiting) ----------------------- */
  for (U64 i = 0; i < n_scaf; ++i) if (scaf_seqs[i].seq) free(scaf_seqs[i].seq);
  free(scaf_seqs);
  arrayDestroy(aligns);
  arrayDestroy(contigs);
  dictDestroy(scaf_dict);
  return 0;
}

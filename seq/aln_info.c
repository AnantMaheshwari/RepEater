/* aln_info.c
 *
 * Walk a .1aln file using the ONElib API and print summary statistics 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "array.h"
#include "dict.h"
#include "ONElib.h"

typedef struct {
    I64 scaf_idx;     /* index into scaf_dict */
    I64 len;          /* contig length (from C line) */
} Contig;

typedef struct {
    I64  total_len;   /* sum of contig lengths for this scaffold */
    char has_aln;     /* 1 if any A line references a contig of this scaffold */
} ScafInfo;

static DICT *load_names(const char *path) {
    FILE *fp = fopen(path, "r");
    if (!fp) die("failed to open %s", path);
    DICT *d = dictCreate(256);
    char buf[8192];
    while (fgets(buf, sizeof(buf), fp)) {
        size_t n = strlen(buf);
        while (n > 0 && (buf[n-1] == '\n' || buf[n-1] == '\r' ||
                         buf[n-1] == ' '  || buf[n-1] == '\t'))
            buf[--n] = '\0';
        if (n > 0) {
            U64 k;
            dictAdd(d, buf, &k);
        }
    }
    fclose(fp);
    return d;
}

int main(int argc, char *argv[]) {
    if (argc < 2 || argc > 4) {
        fprintf(stderr, "usage: aln_info <file.1aln> [--names <names.txt>]\n");
        return 1;
    }
    const char *aln_path = argv[1];
    const char *names_path = NULL;
    if (argc == 4) {
        if (strcmp(argv[2], "--names") != 0) {
            fprintf(stderr, "aln_info: unknown option %s\n", argv[2]);
            return 1;
        }
        names_path = argv[3];
    }

    DICT *restrict_names = names_path ? load_names(names_path) : NULL;

    OneFile *of = oneFileOpenRead(aln_path, 0, 0, 1);
    if (!of) die("failed to open %s", aln_path);

    DICT  *scaf_dict = dictCreate(4096);
    Array  contigs   = arrayCreate(8192, Contig);
    Array  scafs     = arrayCreate(4096, ScafInfo);
    I64    cur_scaf  = -1;
    I64    aln_count = 0;

    while (oneReadLine(of)) {
        char t = of->lineType;
        if (t == 'S') {
            char *name = oneString(of);
            U64 k;
            dictAdd(scaf_dict, name, &k);
            cur_scaf = (I64)k;
            /* ensure scafs has a slot for this index, zero-initialised */
            if ((I64)arrayMax(scafs) <= cur_scaf) {
                ScafInfo *s = arrayp(scafs, cur_scaf, ScafInfo);
                s->total_len = 0;
                s->has_aln   = 0;
            }
        } else if (t == 'C') {
            I64 clen = oneInt(of, 0);
            Contig *c = arrayp(contigs, arrayMax(contigs), Contig);
            c->scaf_idx = cur_scaf;
            c->len      = clen;
            if (cur_scaf >= 0) {
                ScafInfo *s = arrp(scafs, cur_scaf, ScafInfo);
                s->total_len += clen;
            }
        } else if (t == 'A') {
            aln_count++;
            I64 a_contig = oneInt(of, 0);
            if (a_contig >= 0 && a_contig < (I64)arrayMax(contigs)) {
                Contig *c = arrp(contigs, a_contig, Contig);
                if (c->scaf_idx >= 0) {
                    ScafInfo *s = arrp(scafs, c->scaf_idx, ScafInfo);
                    s->has_aln = 1;
                }
            }
        }
    }
    oneFileClose(of);

    /* Pass 2: aggregate totals + active-scaffold bp + (optional) restricted bp */
    I64 total_bp = 0, active_bp = 0, restricted_bp = 0;
    U64 n_scafs = arrayMax(scafs);
    for (U64 i = 0; i < n_scafs; i++) {
        ScafInfo *s = arrp(scafs, i, ScafInfo);
        total_bp += s->total_len;
        if (s->has_aln) active_bp += s->total_len;
        if (restrict_names) {
            char *name = (char *)dictName(scaf_dict, i);
            U64 k;
            if (name && dictFind(restrict_names, name, &k))
                restricted_bp += s->total_len;
        }
    }

    printf("alignments=%lld\n",         (long long)aln_count);
    printf("total_scaffold_bp=%lld\n",  (long long)total_bp);
    printf("active_scaffold_bp=%lld\n", (long long)active_bp);
    if (restrict_names)
        printf("restricted_scaffold_bp=%lld\n", (long long)restricted_bp);
    printf("active_scaffolds:\n");
    for (U64 i = 0; i < n_scafs; i++) {
        ScafInfo *s = arrp(scafs, i, ScafInfo);
        if (s->has_aln) {
            char *name = (char *)dictName(scaf_dict, i);
            if (name) printf("%s\n", name);
        }
    }

    arrayDestroy(scafs);
    arrayDestroy(contigs);
    dictDestroy(scaf_dict);
    if (restrict_names) dictDestroy(restrict_names);
    return 0;
}

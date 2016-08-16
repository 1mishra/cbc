
#include "edit.h"
#include <iostream>
#include <vector>
#include "sam_block.h"

#define UINT32_MAX ((uint32_t) - 1)
#define INT32_MIN (-0x7fffffff - 1)
#define min(a,b) \
  ({ __typeof__ (a) _a = (a); \
   __typeof__ (b) _b = (b); \
   _a < _b ? _a : _b; })

#define max(a,b) \
  ({ __typeof__ (a) _a = (a); \
   __typeof__ (b) _b = (b); \
   _a > _b ? _a : _b; })

#define DEBUG false
#define VERIFY false

using namespace std;

static uint32_t const EDITS = 3;


struct entry {
  int32_t value;
  enum {DIAGONAL, LEFT, UP, NONE} backpointer;
};

static uint32_t min_index(uint32_t *array, uint32_t size) {
  uint32_t min = UINT32_MAX;
  uint32_t index = 0;
  for (uint32_t i = 0; i < size; i++) {
    if (array[i] < min) {
      min = array[i];
      index = i;
    }
  }
  return index;
}

static uint32_t max_index(int32_t *array, uint32_t size) {
  int32_t max = INT32_MIN;
  uint32_t index = 0;
  for (uint32_t i = 0; i < size; i++) {
    if (array[i] > max) {
      max = array[i];
      index = i;
    }
  }
  return index;
}

void init_sequence(struct sequence *seq, uint32_t *Dels, struct ins *Insers, struct snp *SNPs) {
  seq->n_dels = 0, seq->n_ins = 0, seq->n_snps = 0;
  seq->Dels = Dels;
  seq->Insers = Insers;
  seq->SNPs = SNPs;
}

static uint32_t edit_dist_helper(char *str1, char *str2, uint32_t s1, uint32_t s2, vector<vector<uint32_t> > matrix) {
  for (uint32_t i = 0; i <= s1; i++) {
    matrix[i][0] = i;
  }
  for (uint32_t j = 0; j <= s2; j++) {
    matrix[0][j] = j;
  }
  matrix[0][0] = 0;
  for (uint32_t i = 1; i <= s1; i++) {
    for (uint32_t j = 1; j <= s2; j++) {
      if (str1[i-1] == str2[j-1]) {
        matrix[i][j] = matrix[i-1][j-1];
      } else {
        matrix[i][j] = min(matrix[i-1][j-1] + 1, min(matrix[i-1][j] + 1, matrix[i][j-1] + 1));
      }
    }
  }
  return matrix[s1][s2];
}

/************************
 * Compute edit distance
 **********************/

uint32_t edit_dist(char *str1, char *str2, uint32_t s1, uint32_t s2) {

  vector<vector<uint32_t>> matrix(s1+1, vector<uint32_t>(s2 + 1));
  uint32_t dist = edit_dist_helper(str1, str2, s1, s2, matrix);
  return dist;
}

static void fill_target(char *ref, char *target, int prev_pos, int cur_pos, uint32_t *ref_pos, uint32_t *Dels, uint32_t *dels_pos, uint32_t numDels) {

  uint32_t ref_start = *ref_pos;
  if (prev_pos == cur_pos) {
    return;
  }
  for (int i = prev_pos; i < cur_pos; i++) {
    target[i] = ref[*ref_pos];
    (*ref_pos)++;
    while (*dels_pos < numDels && *ref_pos >= Dels[*dels_pos]) {
      if (DEBUG) printf("DELETE %d\n", Dels[*dels_pos]);
      (*ref_pos)++; 
      (*dels_pos)++;
    }

  }
  if (VERIFY) assert(*dels_pos <= numDels);
  if (DEBUG) printf("MATCH [%d, %d), ref [%d, %d)\n", prev_pos, cur_pos, ref_start, *ref_pos);
}

void reconstruct_read_from_ops(struct sequence *seq, char *ref, char *target, uint32_t len) {
  uint32_t start_copy = 0, ref_pos = 0;

  uint32_t ins_pos = 0, snps_pos = 0, dels_pos = 0;
  uint32_t numIns = seq->n_ins, numSnps = seq->n_snps, numDels = seq->n_dels;
  //printf("snps %d, dels %d, ins %d\n", numSnps, numDels, numIns);
  struct ins *Insers = seq->Insers;
  struct snp *SNPs = seq->SNPs;
  uint32_t *Dels = seq->Dels;

  uint32_t buf[2];

  while (numDels > 0 && dels_pos < numDels && ref_pos >= Dels[dels_pos]) {
    if (DEBUG) printf("DELETE %d\n", Dels[dels_pos]);
    ref_pos++;
    dels_pos++;
  }
  for (uint32_t i = 0; i < numIns + numSnps; i++) {
    buf[0] = (ins_pos < numIns)  ? Insers[ins_pos].pos  : UINT32_MAX;
    buf[1] = (snps_pos < numSnps) ? SNPs[snps_pos].pos : UINT32_MAX;
    uint32_t index = min_index(buf, 2);
    if (index == 0) {
      fill_target(ref, target, start_copy, Insers[ins_pos].pos, &ref_pos, Dels, &dels_pos, numDels);
      if (DEBUG) printf("Insert %c at %d\n", basepair2char(Insers[ins_pos].targetChar), Insers[ins_pos].pos);
      target[Insers[ins_pos].pos] = basepair2char(Insers[ins_pos].targetChar);
      start_copy = Insers[ins_pos].pos + 1;
      ins_pos++;
    } else {
      fill_target(ref, target, start_copy, SNPs[snps_pos].pos, &ref_pos, Dels, &dels_pos, numDels);
      if (DEBUG) printf("Replace %c with %c at %d\n", basepair2char(SNPs[snps_pos].refChar), basepair2char(SNPs[snps_pos].targetChar), SNPs[snps_pos].pos);
      target[SNPs[snps_pos].pos] = basepair2char(SNPs[snps_pos].targetChar);
      start_copy = SNPs[snps_pos].pos + 1;
      snps_pos++;
      ref_pos++;
    }
  }
  fill_target(ref, target, start_copy, len, &ref_pos, Dels, &dels_pos, numDels);
  if (VERIFY) assert(snps_pos == numSnps);
  if (VERIFY) assert(dels_pos == ins_pos);
}

// sequence transforms str2 into str1
// callee allocates a large enough array (>= sizeof(operation) * max(s1, 22)) 
// returns number of operations used 
uint32_t edit_sequence(char *str1, char *str2, uint32_t s1, uint32_t s2, struct sequence *seq) {

  vector<vector<uint32_t> > matrix(s1+1, vector<uint32_t>(s2 + 1));

  edit_dist_helper(str1, str2, s1, s2, matrix);

  uint32_t n_dels_tmp = 0, n_ins_tmp = 0, n_snps_tmp = 0; 
  uint32_t Dels_tmp[max(s1,s2)];
  struct ins Insers_tmp[max(s1, s2)];
  struct snp SNPs_tmp[max(s1, s2)];

  uint32_t i = s1;
  uint32_t j = s2;
  uint32_t dist = matrix[i][j];

  while (i != 0 || j != 0) {
    uint32_t edits[EDITS];      
    edits[0] = (i > 0 && j > 0) ? matrix[i-1][j-1] : UINT32_MAX;    // REPLACE
    edits[1] = (j > 0) ? matrix[i][j-1] : UINT32_MAX;               // DELETE
    edits[2] = (i > 0) ? matrix[i-1][j] : UINT32_MAX;               // INSERT
    uint32_t index = min_index(edits, EDITS);
    switch (index) {
      case 0: {
                if (str1[i-1] != str2[j-1]) {
                  SNPs_tmp[n_snps_tmp].targetChar = char2basepair(str1[i-1]);
                  SNPs_tmp[n_snps_tmp].refChar = char2basepair(str2[j-1]);
                  //printf("Replace %c with %c, Reference: %d, targetPos: %d\n", str2[j-1], str1[i-1], (j-1), (i-1));
                  SNPs_tmp[n_snps_tmp].pos = i - 1;
                  n_snps_tmp++;
                }
                i--;
                j--;
                break;
              } 
      case 1: {
                Dels_tmp[n_dels_tmp] = j - 1;
                n_dels_tmp++;
                j--;
                break;
              }
      case 2: {
                Insers_tmp[n_ins_tmp].targetChar = char2basepair(str1[i-1]);
                Insers_tmp[n_ins_tmp].pos = i - 1;
                n_ins_tmp++;
                if (VERIFY) assert(isalpha(str1[i-1]));
                i--;
                break;
              }
    }
  }

  seq->n_dels = n_dels_tmp;
  seq->n_ins = n_ins_tmp;
  seq->n_snps = n_snps_tmp;
  for (uint32_t i = 0; i < n_dels_tmp; i++) {
    seq->Dels[i] = Dels_tmp[n_dels_tmp - i - 1];  
  }
  for (uint32_t i = 0; i < n_ins_tmp; i++) {
    seq->Insers[i] = Insers_tmp[n_ins_tmp - i - 1];  
  }
  for (uint32_t i = 0; i < n_snps_tmp; i++) {
    seq->SNPs[i] = SNPs_tmp[n_snps_tmp - i - 1];  
  }
  if (DEBUG) {
    printf("snps %d, dels %d, ins %d\n", n_snps_tmp, n_dels_tmp, n_ins_tmp);
  }


  if (VERIFY || DEBUG) {
    char buf[1024];
    reconstruct_read_from_ops(seq, str2, buf, s1);
    //assert(edit_dist(str1, buf, s1, s1) == 0);
    if (DEBUG) printf("Ref: %.100s\nTar: %.100s\nAtt: %.100s\n", str2, str1, buf);
    //if (DEBUG) printf("Distance between the reconstructed and ref: %d\n", edit_dist(str1, buf, s1, s1));
    //if (DEBUG) printf("Computed dist: %d\n", (int) dist);
  }

  return dist;
}

/*
static uint32_t glocal_align(char *str1, char *str2, uint32_t s1, uint32_t s2, struct entry matrix[s1+1][s2+1], uint32_t *i_max, uint32_t *j_max) {

  for (uint32_t i = 0; i <= s1; i++) {
    for (uint32_t j = 0; j <= s2; j++) {
      struct entry ent;
      ent.backpointer = NONE;
      ent.value = 0;
      if (i == 0) {
        matrix[i][j] = ent;
      } else if (j == 0) {
        matrix[i][j] = ent;
      } else {
        int32_t edits[4];
        bool match = str1[i-1] == str2[j-1];
        edits[0] = (match) ? matrix[i-1][j-1].value + 1 : matrix[i-1][j-1].value - 1;
        edits[1] = matrix[i-1][j].value - 2;
        edits[2] = matrix[i][j-1].value - 2;
        uint32_t index = max_index(edits, 3);
        ent.value = edits[index];
        if (index == 0) ent.backpointer = DIAGONAL;
        else if (index == 1) ent.backpointer = UP;
        else if (index == 2) ent.backpointer = LEFT;
        matrix[i][j] = ent;
      }
      if (i == s1 && j >= s1 && matrix[i][j].value > matrix[*i_max][*j_max].value) {
        *i_max = i;
        *j_max = j;
      }
    }
  }
  return matrix[*i_max][*j_max].value;
}

uint32_t glocal_sequence(char *str1, char *str2, uint32_t s1, uint32_t s2, struct sequence *seq) {

  struct entry matrix[s1+1][s2+1];

  uint32_t i = s1;
  uint32_t j = s1;

  uint32_t dist = glocal_align(str1, str2, s1, s2, matrix, &i, &j);
  uint32_t i_copy = i;
  uint32_t j_copy = j;
  if (DEBUG) printf("i: %d, j: %d\n", i, j);


  uint32_t n_dels_tmp = 0, n_ins_tmp = 0, n_snps_tmp = 0; 
  uint32_t Dels_tmp[max(s1,s2)];
  struct ins Insers_tmp[max(s1, s2)];
  struct snp SNPs_tmp[max(s1, s2)];

  while (i != 0 || j != 0) {
    printf("%d, %d\n", i, j);
    switch (matrix[i][j].backpointer) {
      case DIAGONAL: {
                if (str1[i-1] != str2[j-1]) {
                  SNPs_tmp[n_snps_tmp].targetChar = char2basepair(str1[i-1]);
                  SNPs_tmp[n_snps_tmp].refChar = char2basepair(str2[j-1]);
                  //printf("Replace %c with %c, Reference: %d, targetPos: %d\n", str2[j-1], str1[i-1], (j-1), (i-1));
                  SNPs_tmp[n_snps_tmp].pos = i - 1;
                  n_snps_tmp++;
                }
                i--;
                j--;
                break;
              } 
      case LEFT: {
                Dels_tmp[n_dels_tmp] = j - 1;
                n_dels_tmp++;
                j--;
                break;
              }
      case UP: {
                Insers_tmp[n_ins_tmp].targetChar = char2basepair(str1[i-1]);
                Insers_tmp[n_ins_tmp].pos = i - 1;
                n_ins_tmp++;
                if (VERIFY) assert(isalpha(str1[i-1]));
                i--;
                break;
              }
      default: 
               printf("i: %d, j:%d\n", i, j);
    }
  }
  assert(j == 0);

  for (uint32_t k = 0; k < n_dels_tmp; k++) {
    seq->Dels[k] = Dels_tmp[n_dels_tmp - k - 1];  
  }
  for (uint32_t k = 0; k < n_ins_tmp; k++) {
    seq->Insers[k] = Insers_tmp[n_ins_tmp - k - 1];  
  }
  for (uint32_t k = 0; k < n_snps_tmp; k++) {
    seq->SNPs[k] = SNPs_tmp[n_snps_tmp - k - 1];  
  }

  seq->n_dels = n_dels_tmp;
  seq->n_ins = n_ins_tmp;
  seq->n_snps = n_snps_tmp;

  if (DEBUG) {
    printf("snps %d, dels %d, ins %d\n", n_snps_tmp, n_dels_tmp, n_ins_tmp);
  }


  if (VERIFY || DEBUG) {
    char buf[1024];
    reconstruct_read_from_ops(seq, str2, buf, s1);
    if (DEBUG) printf("Ref: %.15s\nTar: %.5s\nAtt: %.5s\n", str2, str1, buf);
    assert(edit_dist(str1, buf, s1, s1) == 0);
  }

  return dist;
}*/


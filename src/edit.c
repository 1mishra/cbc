
#include "edit.h"
#include <assert.h>
#include <ctype.h>

#define min(a,b) \
 ({ __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b; })

#define max(a,b) \
 ({ __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b; })

#define DEBUG true

static uint32_t const EDITS = 3;

uint32_t min_index(uint32_t *array, uint32_t size) {
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

static uint32_t edit_dist_helper(char *str1, char *str2, uint32_t s1, uint32_t s2, uint32_t matrix[s1+1][s2+1]) {
  for (uint32_t i = 0; i <= s1; i++) {
    for (uint32_t j = 0; j <= s2; j++) {
      if (i == 0) {
        matrix[i][j] = j;
      } else if (j == 0) {
        matrix[i][j] = i;
      } else if (str1[i-1] == str2[j-1]) {
        matrix[i][j] = matrix[i-1][j-1];
      } else {
        matrix[i][j] = min(matrix[i-1][j-1], min(matrix[i-1][j], matrix[i][j-1])) + 1;
      }
    }
  }
  return matrix[s1][s2];
}

/************************
 * Compute edit distance
 **********************/

uint32_t edit_dist(char *str1, char *str2, uint32_t s1, uint32_t s2) {

  uint32_t matrix[s1+1][s2+1];
  uint32_t dist = edit_dist_helper(str1, str2, s1, s2, matrix);
  return dist;
}

void reconstruct_read_from_ops(struct operation *ops, uint32_t ops_len, char *ref, char *target) {
  uint32_t ref_ctr = 0, target_ctr = 0;
  for (uint32_t i = 0; i < ops_len; i++) {
    switch (ops[i].edit_op) {
      case MATCH: {
        for (int k = 0; k < ops[i].value; k++) {
          target[target_ctr] = ref[ref_ctr];
          target_ctr++;
          ref_ctr++;
        }
        break;
                  }
      case REPLACE: {
        target[target_ctr] = (char) ops[i].value;
        target_ctr++;
        ref_ctr++;
        break;
                    }
      case DELETE: {
        ref_ctr += ops[i].value;
        break;
                   }
      case INSERT:
        target[target_ctr] = (char) ops[i].value;
        target_ctr++;
        break;

    }
  }
}

static uint32_t compact_seq(struct operation *tmp_seq, uint32_t count, struct operation *seq) {
  uint32_t run = 0;

  uint32_t length = 0;
  uint32_t total_structs = 0;

  uint32_t num_dels = 0;
  uint32_t num_ins = 0;

  edit edit_type = DELETE;
  for (int32_t i = count - 1; i >= 0; i--) {
    if (run != 0 && edit_type != tmp_seq[i].edit_op) {
      seq->edit_op = edit_type;
      seq->value = run;
      
      if (edit_type == MATCH) {
        if (DEBUG) printf("MATCH %d\n", (int) seq->value);
        length += seq->value;
      } else {
        num_dels += seq->value;
        if (DEBUG) printf("DELETE %d\n", (int) seq->value);
      }
      run = 0;
      seq++;
      total_structs++;
    }
    switch (tmp_seq[i].edit_op) {
      case MATCH:
        run++;
        edit_type = MATCH;
        break;
      case REPLACE: 
        seq->edit_op = REPLACE;
        seq->value = tmp_seq[i].value;
        if (DEBUG) {
          printf("REPLACE %c\n", (char) tmp_seq[i].value);
        }
        total_structs++;
        length++;
        seq++;
        break;
      case INSERT:
        seq->edit_op = INSERT;
        seq->value = tmp_seq[i].value;
        if (DEBUG) {
          printf("INSERT %c\n", (char) tmp_seq[i].value);
        }
        num_ins++;

        total_structs++;
        length++;
        seq++;
        break;
      case DELETE:
        run++;
        edit_type = DELETE;
        break;
    }
  }
  if (run > 0) {
    seq->edit_op = edit_type;
    seq->value = run;
    if (edit_type == MATCH) {
      length += seq->value;
      if (DEBUG) printf("MATCH %d\n", (int) seq->value);
    } else {
      num_dels += seq->value;
      if (DEBUG) printf("DELETE %d\n", (int) seq->value);
    }
    total_structs++;
  }
  return total_structs;
}

// sequence transforms str2 into str1
// callee allocates a large enough array (>= sizeof(operation) * max(s1, 22)) 
// returns number of operations used 
uint32_t edit_sequence(char *str1, char *str2, uint32_t s1, uint32_t s2, struct operation *seq) {
  uint32_t matrix[s1+1][s2+1];

  uint32_t dist = edit_dist_helper(str1, str2, s1, s2, matrix);

  // change this later to be more robust
  struct operation tmp_seq[1024];

  uint32_t i = s1;
  uint32_t j = s2;
  size_t count = 0;


  while (i != 0 || j != 0) {
    uint32_t edits[EDITS];      
    edits[0] = (i > 0 && j > 0) ? matrix[i-1][j-1] : UINT32_MAX;    // REPLACE
    edits[1] = (j > 0) ? matrix[i][j-1] : UINT32_MAX;               // DELETE
    edits[2] = (i > 0) ? matrix[i-1][j] : UINT32_MAX;               // INSERT
    uint32_t index = min_index(edits, EDITS);
    switch (index) {
      case 0: {
                if (str1[i-1] != str2[j-1]) {
                  tmp_seq[count].value = str1[i-1];
                  tmp_seq[count].edit_op = REPLACE;
                } else { 
                  tmp_seq[count].edit_op = MATCH;
                }
                i--;
                j--;
                break;
              } 
      case 1: {
                tmp_seq[count].edit_op = DELETE;
                j--;
                break;
              }
      case 2: {
                tmp_seq[count].edit_op = INSERT;
                tmp_seq[count].value = str1[i-1];
                assert(isalpha(str1[i-1]));
                i--;
                break;
              }
    }
    count++;
  }
  if (count > 1024) {
    // fail fast until refactoring 
    assert(count <= 1024);
  }

  uint32_t ops_len = compact_seq(tmp_seq, count, seq);
  char buf[1024];
  reconstruct_read_from_ops(seq, ops_len, str2, buf);
  printf("Ori: %.100s\nRef: %.100s\nAtt: %.100s\n", str2, str1, buf);
  printf("Computed dist: %d\n", (int) dist);
  exit(1);
  return dist;
}

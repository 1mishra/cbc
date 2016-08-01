
#include "edit.h"

#define min(a,b) \
 ({ __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b; })

#define max(a,b) \
 ({ __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b; })

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

uint32_t edit_dist_helper(char *str1, char *str2, uint32_t s1, uint32_t s2, uint32_t **matrix) {
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

  uint32_t **matrix = (uint32_t **) malloc((s1 + 1) * sizeof(uint32_t*));
  for (uint32_t i = 0; i <= s1; i++) {
    matrix[i] = (uint32_t *) malloc((s2+1) * sizeof(uint32_t));
  }
  uint32_t dist = edit_dist_helper(str1, str2, s1, s2, matrix);

  for (uint32_t i = 0; i <= s1; i++) {
    free(matrix[i]);
  }
  free(matrix);
  return dist;
}

static void compact_seq(struct operation *tmp_seq, uint32_t count, struct operation *seq) {
  uint32_t run = 0;
  edit edit_type = INSERT;
  for (int32_t i = count - 1; i >= 0; i--) {
    if (run != 0 && edit_type != tmp_seq[i].edit_op) {
      seq->edit_op = edit_type;
      seq->value = run;
      run = 0;
      seq++;
    }
    switch (tmp_seq[i].edit_op) {
      case MATCH:
        run++;
        edit_type = MATCH;
        break;
      case REPLACE: case INSERT:
        seq->edit_op = tmp_seq[i].edit_op;
        seq->value = tmp_seq[i].value;
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
  }
}

// sequence transforms str2 into str1
// callee allocates a large enough array (>= sizeof(operation) * max(s1, 22)) 
// returns number of operations used 
uint32_t edit_sequence(char *str1, char *str2, uint32_t s1, uint32_t s2, struct operation *seq) {
  uint32_t **matrix = (uint32_t **) malloc((s1 + 1) * sizeof(uint32_t*));
  for (uint32_t i = 0; i <= s1; i++) {
    matrix[i] = (uint32_t *) malloc((s2+1) * sizeof(uint32_t));
  }
  uint32_t dist = edit_dist_helper(str1, str2, s1, s2, matrix);

  struct operation *tmp_seq = (struct operation *) malloc(max(s1,s2) * sizeof(struct operation));

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
                i--;
                break;
              }
    }
    count++;
  }
  for (uint32_t i = 0; i <= s1; i++) {
    free(matrix[i]);
  }
  free(matrix);

  compact_seq(tmp_seq, count, seq);
  free(tmp_seq);
  return dist;
}

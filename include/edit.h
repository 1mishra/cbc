#ifndef _EDIT_H_
#define _EDIT_H_

#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>

typedef enum {INSERT, DELETE, REPLACE, MATCH} edit;

struct operation {
  edit edit_op;
  // cast this to a char if REPLACE/INSERT
  // otherwise contains the number of deletions/matches in a row
  int value;
};

static uint32_t edit_dist_helper(char *str1, char *str2, uint32_t s1, uint32_t s2, uint32_t **matrix);
uint32_t edit_dist(char *str1, char *str2, uint32_t s1, uint32_t s2);

void reconstruct_read_from_ops(struct operation *ops, uint32_t ops_len, char *ref, char *target);
static uint32_t compact_seq(struct operation *tmp_seq, uint32_t count, struct operation *seq);
uint32_t edit_sequence(char *str1, char *str2, uint32_t s1, uint32_t s2, struct operation *seq);
#endif 

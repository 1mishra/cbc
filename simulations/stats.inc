//
//  stats.c
//  samSecComp
//
//  Created by Mikel Hernaez on 6/2/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//
//
//  Created by Mikel Hernaez on 6/2/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License version 3,
// as published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Affero General Public License for more details.
// 
// You should have received a copy of the GNU Affero General Public
// License along with this program. If not, see
// <http://www.gnu.org/licenses/>.

#include <stdio.h>
#include <stdlib.h>

#define MAX_CARDINALITY 1000000
#define MAX_ALPHA 50000000

typedef struct stream_stats_t{
    
    int32_t *alphabet;
    uint32_t *counts;
    int32_t *alphaMap;
    uint8_t *alphaExist;
    uint32_t alphabetCard;
    uint32_t step;
    uint32_t n;
    
}*stream_stats;

//////////////////////////////////////////////////////////////////////////////////////////
//                                                                                      //
//                                                                                      //
//                                  INITIALIZATION                                      //
//                                                                                      //
//                                                                                      //
//////////////////////////////////////////////////////////////////////////////////////////

stream_stats* initialize_stream_stats_F(){
    
    stream_stats *s = (stream_stats*) calloc(1, sizeof(stream_stats));
    
    s[0] = (stream_stats) calloc(1, sizeof(struct stream_stats_t));
    
    // Allocate memory
    s[0]->counts = (uint32_t*) calloc(8, sizeof(uint32_t));
    
    // An extra for the cumcounts
    s[0]->counts += 1;
    
    s[0]->alphabetCard = 2;
    
    s[0]->counts[0] = 1;
    s[0]->counts[1] = 1;
    
    s[0]->n = 2;
    
    // STEP
    s[0]->step = 1;
    
    return s;
    
}

stream_stats* initialize_stream_stats_P(){
    
    
    
    stream_stats *s;
    
    s = (stream_stats*) calloc(1, sizeof(stream_stats));
    
    s[0] = (stream_stats) calloc(1, sizeof(struct stream_stats_t));
    
    // Allocate memory
    s[0]->alphabet = (int32_t *) calloc(MAX_CARDINALITY + 2, sizeof(int32_t));
    s[0]->counts = (uint32_t*) calloc(MAX_CARDINALITY + 3, sizeof(uint32_t));
    s[0]->alphaExist = (uint8_t*) calloc(MAX_ALPHA + 2, sizeof(uint8_t));
    s[0]->alphaMap = (int32_t*) calloc(MAX_ALPHA + 2, sizeof(int32_t));
    
    // We shift the pointers to consider -1 (new alpha).
    s[0]->alphaMap += 1;
    s[0]->alphaExist += 1;
    s[0]->alphabet += 1;
    // An extra for the cumcounts
    s[0]->counts += 2;
    
    
    s[0]->alphabetCard = 0;
    s[0]->alphabet[-1] = -1;
    s[0]->alphaMap[-1] = -1;
    
    s[0]->alphaExist[-1] = 1;
    s[0]->counts[-1] = 1;
    s[0]->n = 1;
    
    // STEP
    s[0]->step = 10;
    
    return s;
    
}

stream_stats* initialize_stream_stats_M(){
    
    stream_stats *s = (stream_stats*) calloc(256, sizeof(stream_stats));
    
    uint32_t i = 0;
    
    for (i = 0; i < 256; i++) {
        
        s[i] = (stream_stats) calloc(1, sizeof(struct stream_stats_t));
        
        // Allocate memory
        s[i]->counts = (uint32_t*) calloc(8, sizeof(uint32_t));
        
        // An extra for the cumcounts
        s[i]->counts += 1;
        
        s[i]->alphabetCard = 2;
        
        s[i]->counts[0] = 1;
        s[i]->counts[1] = 1;
        
        s[i]->n = 2;
        
        // STEP
        s[i]->step = 1;
    }
    
    
    
    
    
    return s;
    
}

stream_stats* initialize_stream_stats_S(uint32_t readLength){
    
    stream_stats *s;
    
    uint32_t i = 0;
    
    s = (stream_stats*) calloc(1, sizeof(stream_stats));
    
    s[0] = (stream_stats) calloc(1, sizeof(struct stream_stats_t));
    
    // Allocate memory
    s[0]->counts = (uint32_t*) calloc(readLength + 2, sizeof(uint32_t));
    
    // An extra for the cumcounts
    s[0]->counts += 1;
    
    s[0]->alphabetCard = readLength;
    
    s[0]->n = 0;
    for (i = 0; i < readLength; i++) {
        s[0]->counts[i] = 1;
        s[0]->n += s[0]->counts[i];
    }
    
    // STEP
    s[0]->step = 10;
    
    return s;
    
}

stream_stats* initialize_stream_stats_I(uint32_t readLength){
    
    stream_stats *s;
    
    s = (stream_stats*) calloc(1, sizeof(stream_stats));
    
    uint32_t i = 0;
    
    s[0] = (stream_stats) calloc(1, sizeof(struct stream_stats_t));
    
    // Allocate memory
    s[0]->counts = (uint32_t*) calloc(readLength + 2, sizeof(uint32_t));
    
    // An extra for the cumcounts
    s[0]->counts += 1;
    
    s[0]->alphabetCard = readLength;
    
    s[0]->n = 0;
    for (i = 0; i < readLength; i++) {
        s[0]->counts[i] = 1;
        s[0]->n += s[0]->counts[i];
    }
    
    // STEP
    s[0]->step = 16;
    
    return s;
    
}

stream_stats* initialize_stream_stats_pos(uint32_t readLength){
    
    stream_stats* s;
    
    uint32_t i = 0, j = 0;
    
    uint32_t num_models = 65536;
    
    s = (stream_stats*) calloc(num_models, sizeof(stream_stats));
    
    for (j = 0; j < num_models; j++) {
    
        s[j] = (stream_stats) calloc(1, sizeof(struct stream_stats_t));
    
        // Allocate memory
        s[j]->counts = (uint32_t*) calloc(readLength + 2, sizeof(uint32_t));
    
        // An extra for the cumcounts
        s[j]->counts += 1;
    
        s[j]->alphabetCard = readLength;
    
        s[j]->n = 0;
        for (i = 0; i < readLength; i++) {
            s[j]->counts[i] = 1;
            s[j]->n += s[0]->counts[i];
        }
    
        // STEP
        s[j]->step = 10;
    }
    
    return s;
    
}

stream_stats* initialize_stream_stats_char(){
    
    stream_stats* s;
    
    s = (stream_stats*) calloc(6, sizeof(stream_stats));
    
    uint32_t i = 0, j;
    
    for (j = 0; j < 6; j++) {
        
        s[j] = (stream_stats) calloc(1, sizeof(struct stream_stats_t));
        
        // Allocate memory
        s[j]->counts = (uint32_t*) calloc(16, sizeof(uint32_t));
        
        // An extra for the cumcounts
        s[j]->counts += 1;
        
        s[j]->alphabetCard = 5;
        
        
        s[j]->n = 0;
        for (i = 0; i < 4; i++) {
            s[j]->counts[i] = (i == j)? 0:8;
            s[j]->n += s[j]->counts[i];
        }
        //{A,C,G,T} to N
        s[j]->counts[4] = 1;
        s[j]->n ++;
        
        // STEP
        s[j]->step = 8;
    }
 
    s[0]->counts[1] += 8;
    s[0]->counts[2] += 8;
    s[0]->n += 16;
    
    s[1]->counts[0] += 8;
    s[1]->counts[3] += 8;
    s[1]->n += 16;
    
    s[2]->counts[0] += 8;
    s[2]->counts[3] += 8;
    s[2]->n += 16;
    
    s[3]->counts[1] += 8;
    s[3]->counts[2] += 8;
    s[3]->n += 16;
    
    // Allow N to N
    // s[4]->counts[4] = 1, s[4]->n++;
    
    // We need to increase the total count of the insertion case
    // s[5]->n++;
    
    return s;
    
}





//////////////////////////////////////////////////////////////////////////////////////////
//                                                                                      //
//                                                                                      //
//                                  RESCALE                                             //
//                                                                                      //
//                                                                                      //
//////////////////////////////////////////////////////////////////////////////////////////

uint32_t rescale_stats(stream_stats s){
    
    int32_t i = 0;
    
    s->n = 0;
    for (i = 0; i < (int32_t) s->alphabetCard; i++){
        s->counts[i] >>= 1, s->counts[i]++;
        s->n += s->counts[i];
    }
    if (s->counts[-1] != 0) {
        s->counts[-1] >>= 1, s->counts[-1]++;
        s->n += s->counts[-1];
    }
    
    return s->n;
}

uint32_t rescale_stats_var(stream_stats s){
    
    int32_t i = 0, tempCard = 0;
    
    s->n = 0;
    for (i = -1; i < (int32_t) s->alphabetCard; i++){
        s->counts[i] >>= 1;
        if (s->counts[i] == 0) {
            s->alphaExist[s->alphabet[i]] = 0;
            tempCard--;
        }
        else s->n += s->counts[i];
    }
    
    return s->n;
}

//////////////////////////////////////////////////////////////////////////////////////////
//                                                                                      //
//                                                                                      //
//                                  UPDATE                                              //
//                                                                                      //
//                                                                                      //
//////////////////////////////////////////////////////////////////////////////////////////

uint32_t update_stats_F(stream_stats stats, int32_t x, uint32_t m){
    
    // Update the statistics
    stats->counts[x]+= stats->step, stats->n+= stats->step;
    
    // Rescale if necessary
    if (stats->n >= (1 << (m - 3))){
        
        stats->n = 0;
        stats->counts[0] >>= 1, stats->counts[0]++, stats->n += stats->counts[0];
        stats->counts[1] >>= 1, stats->counts[1]++, stats->n += stats->counts[1];
    }
    
    return 1;
}

uint32_t update_stats_P(stream_stats stats, int32_t x, uint32_t m, uint8_t stepSize){
    
    int32_t i = 0;
    
    // Update the statistics
    if (stepSize)
        stats->counts[x]+= stats->step, stats->n+= stats->step;
    else
        stats->counts[x]+= 1, stats->n+= 1;
    
    // Rescale if necessary
    if (stats->n >= (1 << (m - 3))){
        
        stats->n = 0;
        for (i = -1; i < (int32_t) stats->alphabetCard; i++){
            stats->counts[i] >>= 1, stats->counts[i]++;
            stats->n += stats->counts[i];
        }
    }
    
    return 1;
}

uint32_t update_stats_M(stream_stats stats, int32_t x, uint32_t m){
    
    // Update the statistics
    stats->counts[x]+= stats->step, stats->n+= stats->step;
    
    // Rescale if necessary
    if (stats->n >= (1 << (m - 3))){
        
        stats->n = 0;
        stats->counts[0] >>= 1, stats->counts[0]++, stats->n += stats->counts[0];
        stats->counts[1] >>= 1, stats->counts[1]++, stats->n += stats->counts[1];
    }
    
    return 1;
}

uint32_t update_stats_S(stream_stats stats, int32_t x, uint32_t m){
    
    int32_t i = 0;
    // Update the statistics
    stats->counts[x]+= stats->step, stats->n+= stats->step;
    
    // Rescale if necessary
    if (stats->n >= (1 << (m - 3))){
        
        stats->n = 0;
        for (i = 0; i < (int32_t) stats->alphabetCard; i++){
            stats->counts[i] >>= 1, stats->counts[i]++;
            stats->n += stats->counts[i];
        }
    }
    
    return 1;
}

uint32_t update_stats_I(stream_stats stats, int32_t x, uint32_t m){
    
    int32_t i = 0;
    // Update the statistics
    stats->counts[x]+= stats->step, stats->n+= stats->step;
    
    // Rescale if necessary
    if (stats->n >= (1 << (m - 3))){
        
        stats->n = 0;
        for (i = 0; i < (int32_t) stats->alphabetCard; i++){
            stats->counts[i] >>= 1, stats->counts[i]++;
            stats->n += stats->counts[i];
        }
    }
    
    return 1;
}

uint32_t update_stats_pos(stream_stats stats, int32_t x, uint32_t m){
    
    int32_t i = 0;
    // Update the statistics
    stats->counts[x]+= stats->step, stats->n+= stats->step;
    
    // Rescale if necessary
    if (stats->n >= (1 << (m - 3))){
        
        stats->n = 0;
        for (i = 0; i < (int32_t) stats->alphabetCard; i++){
            stats->counts[i] >>= 1, stats->counts[i]++;
            stats->n += stats->counts[i];
        }
    }
    
    return 1;
}

uint32_t update_stats_char(stream_stats stats, int32_t x, uint32_t m){
    
    int32_t i = 0;
    // Update the statistics
    stats->counts[x]+= stats->step, stats->n+= stats->step;
    
    // Rescale if necessary
    if (stats->n >= (1 << (m - 3))){
        
        stats->n = 0;
        for (i = 0; i < (int32_t) stats->alphabetCard; i++){
            stats->counts[i] >>= 1, stats->counts[i]++;
            stats->n += stats->counts[i];
        }
    }
    
    return 1;
}

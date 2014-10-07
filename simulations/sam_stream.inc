//
//  sam_stream.c
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
#include <string.h>
#include "arith.inc"

#define mF 18
#define mP 22
#define mM 14
#define mS 18
#define mI 18
#define mp 18
#define mc 13



typedef struct samStream_t{
    stream_stats* stats;
    Arithmetic_code a;
    osStream os;
}*samStream;

typedef struct sam_compressor_t{
    samStream F;
    samStream P;
    samStream M;
    samStream S;
    samStream I;
    samStream p;
    samStream c;
}*sam_compressor;



samStream initialize_samStream_F(char* osPath, uint8_t decompressor_flag){
    
    samStream sStream;
    FILE *fos;
    
    uint32_t osPathLength = (uint32_t)strlen(osPath);
    
    strcat(osPath, "_F.ido");
    fos = (decompressor_flag)? fopen(osPath, "r"):fopen(osPath, "w");
    
    sStream = (samStream) calloc(1, sizeof(struct samStream_t));
    sStream->stats = initialize_stream_stats_F();
    sStream->a = initialize_arithmetic_encoder(mF);
    sStream->os = initialize_osStream(1, fos, NULL, decompressor_flag);
    sStream->a->t = (decompressor_flag)? read_uint32_from_stream(sStream->a->m, sStream->os):0;
    
    *(osPath + osPathLength) = 0;
    
    return sStream;
    
}

samStream initialize_samStream_P(char* osPath, uint8_t decompressor_flag){
    
    samStream sStream;
    FILE *fos, *fos_A;
    
    uint32_t osPathLength = (uint32_t)strlen(osPath);
    
    strcat(osPath, "_P.ido");
    fos = (decompressor_flag)? fopen(osPath, "r"):fopen(osPath, "w");
    
    *(osPath + osPathLength) = 0;
    
    strcat(osPath, "_P_A.ido");
    fos_A = (decompressor_flag)? fopen(osPath, "r"):fopen(osPath, "w");
    
    sStream = (samStream) calloc(1, sizeof(struct samStream_t));
    sStream->stats = initialize_stream_stats_P();
    sStream->a = initialize_arithmetic_encoder(mP);
    sStream->os = initialize_osStream(1, fos, fos_A, decompressor_flag);
    sStream->a->t = (decompressor_flag)? read_uint32_from_stream(sStream->a->m, sStream->os):0;
    
    *(osPath + osPathLength) = 0;
    
    return sStream;
    
}

samStream initialize_samStream_M(char* osPath, uint8_t decompressor_flag){
    
    samStream sStream;
    FILE *fos;
    
    uint32_t osPathLength = (uint32_t)strlen(osPath);
    
    strcat(osPath, "_M.ido");
    fos = (decompressor_flag)? fopen(osPath, "r"):fopen(osPath, "w");
    
    sStream = (samStream) calloc(1, sizeof(struct samStream_t));
    sStream->stats = initialize_stream_stats_M();
    sStream->a = initialize_arithmetic_encoder(mM);
    sStream->os = initialize_osStream(1, fos, NULL, decompressor_flag);
    sStream->a->t = (decompressor_flag)? read_uint32_from_stream(sStream->a->m, sStream->os):0;
    
    *(osPath + osPathLength) = 0;
    
    return sStream;
    
}

samStream initialize_samStream_S(char* osPath, uint8_t decompressor_flag, uint32_t readLength){
    
    samStream sStream;
    FILE *fos;
    
    uint32_t osPathLength = (uint32_t)strlen(osPath);
    
    strcat(osPath, "_S.ido");
    fos = (decompressor_flag)? fopen(osPath, "r"):fopen(osPath, "w");
    
    sStream = (samStream) calloc(1, sizeof(struct samStream_t));
    sStream->stats = initialize_stream_stats_S(readLength);
    sStream->a = initialize_arithmetic_encoder(mS);
    sStream->os = initialize_osStream(1, fos, NULL, decompressor_flag);
    sStream->a->t = (decompressor_flag)? read_uint32_from_stream(sStream->a->m, sStream->os):0;
    
    *(osPath + osPathLength) = 0;
    
    return sStream;
    
}

samStream initialize_samStream_I(char* osPath, uint8_t decompressor_flag, uint32_t readLength){
    
    samStream sStream;
    FILE *fos;
    
    
    uint32_t osPathLength = (uint32_t)strlen(osPath);
    
    strcat(osPath, "_I.ido");
    fos = (decompressor_flag)? fopen(osPath, "r"):fopen(osPath, "w");
    
    sStream = (samStream) calloc(1, sizeof(struct samStream_t));
    sStream->stats = initialize_stream_stats_I(readLength);
    sStream->a = initialize_arithmetic_encoder(mI);
    sStream->os = initialize_osStream(1, fos, NULL, decompressor_flag);
    sStream->a->t = (decompressor_flag)? read_uint32_from_stream(sStream->a->m, sStream->os):0;
    
    *(osPath + osPathLength) = 0;
    
    return sStream;
    
}

samStream initialize_samStream_pos(char* osPath, uint8_t decompressor_flag, uint32_t readLength){
    
    samStream sStream;
    FILE *fos;
    
    uint32_t osPathLength = (uint32_t)strlen(osPath);
    
    strcat(osPath, "_pos.ido");
    fos = (decompressor_flag)? fopen(osPath, "r"):fopen(osPath, "w");
    
    sStream = (samStream) calloc(1, sizeof(struct samStream_t));
    sStream->stats = initialize_stream_stats_pos(readLength);
    sStream->a = initialize_arithmetic_encoder(mp);
    sStream->os = initialize_osStream(1, fos, NULL, decompressor_flag);
    sStream->a->t = (decompressor_flag)? read_uint32_from_stream(sStream->a->m, sStream->os):0;
    
    *(osPath + osPathLength) = 0;
    
    return sStream;
}

samStream initialize_samStream_char(char* osPath, uint8_t decompressor_flag){
    
    samStream sStream;
    FILE *fos;
    
    uint32_t osPathLength = (uint32_t)strlen(osPath);
    
    strcat(osPath, "_char.ido");
    fos = (decompressor_flag)? fopen(osPath, "r"):fopen(osPath, "w");
    
    sStream = (samStream) calloc(1, sizeof(struct samStream_t));
    sStream->stats = initialize_stream_stats_char();
    sStream->a = initialize_arithmetic_encoder(mc);
    sStream->os = initialize_osStream(1, fos, NULL, decompressor_flag);
    sStream->a->t = (decompressor_flag)? read_uint32_from_stream(sStream->a->m, sStream->os):0;
    
    *(osPath + osPathLength) = 0;
    
    return sStream;
    
}

sam_compressor initialize_sam_compressor(char osPath[], uint8_t streamDirection, uint32_t *readLength){
    
    sam_compressor s;
    
    s = calloc(1, sizeof(struct sam_compressor_t));
    
    
    s->F = initialize_samStream_F(osPath, streamDirection);
    s->P = initialize_samStream_P(osPath, streamDirection);
    if (streamDirection == 1){
        // Read READ_LENGTH
        fscanf(s->P->os->fosA, "%d", readLength);
    }
    s->M = initialize_samStream_M(osPath, streamDirection);
    s->S = initialize_samStream_S(osPath, streamDirection, *readLength);
    s->I = initialize_samStream_I(osPath, streamDirection, *readLength);
    s->p = initialize_samStream_pos(osPath, streamDirection, *readLength);
    s->c = initialize_samStream_char(osPath, streamDirection);
    
    return s;
}

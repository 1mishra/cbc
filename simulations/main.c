//
//  main.cpp
//  samComp
//
//  Created by Mikel Hernaez on 4/13/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//
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


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>


#include <stdint.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <ctype.h>

#include "input_file.inc"
#include "sam_stream.inc"

#define isnumber isdigit
#define MAX_READ_LENGTH 1024
#define MAX_PATH_LENGTH 256
#define MAX_BP_CHR 300000000
#define MAX_HEADER 512
#define COMPRESSION 0
#define DECOMPRESSION 1

#define BITS_DELTA 7

uint32_t READ_LENGTH;
char _readLength[4];

char *reference;

uint8_t snpInRef[MAX_BP_CHR] = {0};
uint32_t cumsumP = 0;

//**************************************************************//
//                                                              //
//                  STORE REFERENCE IN MEMORY                   //
//                                                              //
//**************************************************************//
int store_reference_in_memory(FILE* refFile){
    uint32_t ch, letterCount;
    char *H;
    
    reference = (char *) malloc(MAX_BP_CHR*sizeof(char));
    
    // ******* Read and Store Reference****** //
    letterCount = 0;
    
    ch = getc(refFile);
    if (ch == '>' || ch == '@'){
        H = (char*)malloc(MAX_HEADER * sizeof(char));
        fgets(H, MAX_HEADER, refFile);
        free(H);
    }
    else{
        reference[letterCount++] = toupper(ch);
    }
    
    while (EOF != (ch=getc(refFile)))
        if (ch !='\n')
            reference[letterCount++] = toupper(ch);
    
    reference[letterCount] = '#';
    reference[letterCount+1] = '\0';
    
    reference = (char *) realloc(reference, letterCount+2);
    
    return letterCount+2;
    
}


// To store the stats of the chars both in ref and target

enum BASEPAIR {
    A,
    C,
    G,
    T,
    N,
    O
};

typedef struct ch_t{
    char refChar;
    char targetChar;
}ch_t;

typedef struct snp{
    uint32_t pos;
    enum BASEPAIR refChar;
    enum BASEPAIR targetChar;
    uint32_t ctr;
}snp;

typedef struct ins{
    uint32_t pos;
    enum BASEPAIR targetChar;
}ins;

int help_enc(){
    
    printf("Usage: <inputFile> <compressedFilesPath>\n");
    return 1;
    
}

int help_dec(){
    
    printf("Usage: <outputFile> <FileContainingPathToReferences> <compressedFile.ido>\n");
    return 1;
    
}


//////////////////////////////////////////////////////
//                                                  //
//                                                  //
//                  GENERAL FUNCTIONS               //
//                                                  //
//                                                  //
//////////////////////////////////////////////////////

int char2basepair(char c)
{
    switch(c)
    {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return 4;
    }
}

int basepair2char(enum BASEPAIR c)
{
    switch(c)
    {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return 'N';
    }
}

char bp_complement(char c){
    
    switch(c)
    {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return c;
    }
}






//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
//                                                                                                          //
//                                         COMPRESSION                                                      //
//                                                                                                          //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


uint32_t compute_delta_to_first_snp(prevPos, readLen){
    
    uint32_t deltaOut;
    uint32_t j = 0;
    
    deltaOut = readLen + 2;
    
    for (j=0;j<readLen - prevPos; j++){
        if (snpInRef[cumsumP - 1 + j + prevPos] == 1){
            deltaOut = j;
            break;
        }
    }
    
    return deltaOut;
}



uint32_t compress_F(samStream F, uint16_t flag){
    
    int isReversed = 0;
    
    isReversed = flag << 11;
    isReversed >>= 15;
    
    // Send match to the arithmetic encoder
    arithmetic_encoder_step(F->a, F->stats[0], isReversed, F->os);
    // Update Stats
    update_stats_F(F->stats[0], isReversed, F->a->m);
    
    return isReversed;
    
}

uint32_t decompress_F(samStream F){
    
    uint32_t isReversed = 0;
    
    // Decode the match
    isReversed = arithmetic_decoder_step(F->a, F->stats[0], F->os);
    // Update Stats
    update_stats_F(F->stats[0], isReversed, F->a->m);
    
    return isReversed;
}

uint32_t compress_P(samStream P, uint32_t pos){
    
    static uint32_t prevPos = 0;
    enum {SMALL_STEP = 0, BIG_STEP = 1};
    int i = 0;
    int32_t x = 0;
    
    // Check if we are changing chromosomes.
    if (pos < prevPos) {
        
        cumsumP = 0;
        for (i=0; i<MAX_BP_CHR; i++){
            snpInRef[i] = 0;
        }
        
        
        // Send a -1 to indicate that you should read from the A,
        arithmetic_encoder_step(P->a, P->stats[0], -1, P->os);
        // Update the statistics of the -1
        update_stats_P(P->stats[0], -1, P->a->m, SMALL_STEP);
        
        // Send a -1 to A,
        fprintf(P->os->fosA, "%d\n", -1);
        
        // Compress the position
        x = pos;
        if (P->stats[0]->alphaExist[x]){
            // Send x to the arithmetic encoder
            arithmetic_encoder_step(P->a, P->stats[0], P->stats[0]->alphaMap[x], P->os);
            // Update the statistics
            update_stats_P(P->stats[0], P->stats[0]->alphaMap[x], P->a->m, BIG_STEP);
        }
        else{
            
            // Send a -1 to the arithmetic encoder
            arithmetic_encoder_step(P->a, P->stats[0], -1, P->os);
            // Update the statistics for the -1
            update_stats_P(P->stats[0], -1, P->a->m, SMALL_STEP);
            
            // Send the new letter to the alphabet stream
            fprintf(P->os->fosA, "%d\n", x);
            
            // Update the statistics of the alphabet for x
            P->stats[0]->alphabet[P->stats[0]->alphabetCard] = x;
            P->stats[0]->alphaMap[x] = P->stats[0]->alphabetCard;
            P->stats[0]->alphaExist[x] = 1;
            update_stats_P(P->stats[0], P->stats[0]->alphabetCard++, P->a->m, SMALL_STEP);
        }
        
        
        prevPos = pos;
        
        return x;
    }
    
    
    // Compress the position diference
    x = pos - prevPos;
    
    if (P->stats[0]->alphaExist[x]){
        // Send x to the arithmetic encoder
        arithmetic_encoder_step(P->a, P->stats[0], P->stats[0]->alphaMap[x], P->os);
        // Update the statistics
        update_stats_P(P->stats[0], P->stats[0]->alphaMap[x], P->a->m, BIG_STEP);
    }
    else{
        
        // Send a -1 to the arithmetic encoder
        arithmetic_encoder_step(P->a, P->stats[0], -1, P->os);
        // Update the statistics for the -1
        update_stats_P(P->stats[0], -1, P->a->m, SMALL_STEP);
        
        // Send the new letter to the alphabet stream
        fprintf(P->os->fosA, "%d\n", x);
        
        // Update the statistics of the alphabet for x
        P->stats[0]->alphabet[P->stats[0]->alphabetCard] = x;
        P->stats[0]->alphaMap[x] = P->stats[0]->alphabetCard;
        P->stats[0]->alphaExist[x] = 1;
        update_stats_P(P->stats[0], P->stats[0]->alphabetCard++, P->a->m, SMALL_STEP);
    }
    
    prevPos = pos;
    
    return x;
}

int32_t decompress_P(samStream P){
    
    static uint32_t prevPos = 0;
    
    int32_t pos, alphaMapX = 0, x = 0;
    
    enum {SMALL_STEP = 0, BIG_STEP = 1};
    
    alphaMapX = arithmetic_decoder_step(P->a, P->stats[0], P->os);
    
    x = P->stats[0]->alphabet[alphaMapX];
    
    if (x == -1) {
        
        // Update the statistics
        update_stats_P(P->stats[0], alphaMapX, P->a->m, SMALL_STEP);
        
        // Read from the osA to get the unknown alphabet letter alpha
        fscanf(P->os->fosA, "%d", &x);
        
        // Check if we are changing chromosomes.
        if (x == -1) {
            
            prevPos = 0;
            return -1;
        }
        
        // Update the statistics of the alphabet for x
        P->stats[0]->alphabet[P->stats[0]->alphabetCard] = x;
        P->stats[0]->alphaMap[x] = P->stats[0]->alphabetCard;
        P->stats[0]->alphaExist[x] = 1;
        update_stats_P(P->stats[0], P->stats[0]->alphabetCard++, P->a->m, SMALL_STEP);
    }
    
    else
    {
        // Update the statistics
        update_stats_P(P->stats[0], alphaMapX, P->a->m, BIG_STEP);
        
    }
    
    pos = prevPos + x;
    
    prevPos = pos;
    
    return pos;
}

uint32_t compress_M(samStream M, uint8_t match, uint32_t P){
    
    uint32_t ctx = 0;
    static uint8_t  prevM = 0;
    
    // Compute Context
    P = (P != 0)? 0:1;
    //prevP = (prevP > READ_LENGTH)? READ_LENGTH:prevP;
    //prevP = (prevP > READ_LENGTH/4)? READ_LENGTH:prevP;
    
    ctx = (P << 1) | prevM;
    
    //ctx = 0;
    // Send match to the arithmetic encoder
    arithmetic_encoder_step(M->a, M->stats[ctx], match, M->os);
    
    // Update Stats
    update_stats_M(M->stats[ctx], match, M->a->m);
    
    prevM = match;
    
    return 1;
}

uint32_t decompress_M(samStream M, uint32_t P){
    
    uint32_t match, ctx = 0;
    static uint8_t prevM = 0;
    
    P = (P != 0)? 0:1;
    
    ctx = (P << 1) | prevM;
    
    // Decode the match
    match = arithmetic_decoder_step(M->a, M->stats[ctx], M->os);
    
    // Update Stats
    update_stats_M(M->stats[ctx], match, M->a->m);
    
    prevM = match;
    
    return match;
}

uint32_t compress_S(samStream S, uint8_t numSnps){
    
    // Send match to the arithmetic encoder
    arithmetic_encoder_step(S->a, S->stats[0], numSnps, S->os);
    // Update Stats
    update_stats_S(S->stats[0], numSnps, S->a->m);
    
    return 1;
    
}

uint32_t decompress_S(samStream S){
    
    int32_t numSnps;
    // Decode the symbol
    numSnps = arithmetic_decoder_step(S->a, S->stats[0], S->os);
    // Update Stats
    update_stats_S(S->stats[0], numSnps, S->a->m);
    
    return numSnps;
}

uint32_t compress_I(samStream I, uint8_t numIndels){
    
    // Send numIndels to the arithmetic encoder
    arithmetic_encoder_step(I->a, I->stats[0], numIndels, I->os);
    // Update Stats
    update_stats_I(I->stats[0], numIndels, I->a->m);
    return 1;
    
}

uint32_t decompress_I(samStream I){
    
    int32_t numIndels;
    // Decode the symbol
    numIndels = arithmetic_decoder_step(I->a, I->stats[0], I->os);
    // Update Stats
    update_stats_S(I->stats[0], numIndels, I->a->m);
    
    return numIndels;
}

uint32_t compress_pos(samStream p, uint32_t pos, uint32_t prevPos, uint32_t flag){
    
    uint32_t ctx = 0;
    
    //flag = 0;
    ctx = prevPos << 1 | flag;
    
    // Send pos to the arithmetic encoder
    arithmetic_encoder_step(p->a, p->stats[ctx], pos, p->os);
    // Update Stats
    update_stats_pos(p->stats[ctx], pos, p->a->m);
    
    return 1;
    
}

uint32_t decompress_pos(samStream p, uint32_t prevPos, uint32_t flag){
    
    uint32_t pos, ctx = 0;
    
    ctx = prevPos << 1| flag;
    
    // Decode the pos
    pos = arithmetic_decoder_step(p->a, p->stats[ctx], p->os);
    
    // Update Stats
    update_stats_pos(p->stats[ctx], pos, p->a->m);
    
    return pos;
    
    
}

uint32_t compress_chars(samStream s, enum BASEPAIR ref, enum BASEPAIR target){
    
    // Send target to the arithmetic encoder (stats based on ref)
    arithmetic_encoder_step(s->a, s->stats[ref], target, s->os);
    // Update Stats
    update_stats_char(s->stats[ref], target, s->a->m);
    
    return 1;
    
}

uint8_t decompress_char(samStream s, enum BASEPAIR ref){
    
    uint32_t target;
    
    target = arithmetic_decoder_step(s->a, s->stats[ref], s->os);
    
    // Update Stats
    update_stats_char(s->stats[ref], target, s->a->m);
    
    return basepair2char( (enum BASEPAIR)target);
    
}




int add_snps_to_array(char* edits, snp* SNPs, unsigned int *numSnps, unsigned int insertionPos, char *read){
    
    static unsigned int prevEditPtr = 0, cumPos = 0;
    
    int pos = 0, tempPos = 0, ctr;
    char ch = 0;
    
    uint8_t flag = 0;
    
    edits += prevEditPtr;
    
    while (*edits != 0 ) {
        
        pos = atoi(edits);
        tempPos = pos;
        
        // if there are deletions after pos, we need to add those positions that come after the deletions
        ctr = 0;
        if ( isnumber(*(edits+1)) ) ctr = 2;
        else ctr = 1;
        ch = *(edits+ctr);
        ctr++;
        while (ch == '^'){
            while (isnumber(*(edits+ctr)) == 0) {
                ctr++;
            }
            tempPos += atoi(edits + ctr);
            if ( isnumber(*(edits+ctr+1)) ) ctr += 2;
            else ctr += 1;
            ch = *(edits+ctr);
            ctr++;
            if (ch == '\0'){
                flag = 1;
                break;
            }
        }
        
        if (flag == 1){
            flag = 0;
            break;
        }
        
        
        if (cumPos + tempPos >= insertionPos){
            cumPos++;
            return cumPos;
        }
        
        if ( isnumber(*(edits+1)) ) edits += 2, prevEditPtr += 2;
        else edits += 1, prevEditPtr += 1;
        
        ch = *edits++, prevEditPtr++;
        
        while (ch == '^'){
            while (isnumber(*edits) == 0)
                edits++, prevEditPtr++;
            pos += atoi(edits);
            if ( isnumber(*(edits+1)) ) edits += 2, prevEditPtr += 2;
            else edits += 1, prevEditPtr += 1;
            ch = *edits++, prevEditPtr++;
        }
        
        if (ch == '\0')
            break;
        cumPos += pos;
        
        SNPs[*numSnps].pos = pos;
        SNPs[*numSnps].refChar = char2basepair(ch);
        SNPs[*numSnps].targetChar = char2basepair(read[cumPos]);
        (*numSnps)++;
        cumPos++;
        if (*edits == 0)
            break;
    }
    
    prevEditPtr = 0;
    cumPos = 0;
    return 0;
}

uint32_t compress_edits(sam_compressor s, char *edits, char *cigar, char *read, uint32_t P, uint8_t flag){
    
    unsigned int numIns = 0, numDels = 0, numSnps = 0, lastSnp = 1;
    int i = 0, M = 0, I = 0, D = 0, pos = 0, ctr = 0, prevPosI = 0, prevPosD = 0, ctrS = 0, S = 0;
    uint32_t delta = 0;
    
    // pos in the reference
    cumsumP = cumsumP + P;
    
    uint32_t Dels[MAX_READ_LENGTH];
    ins Insers[MAX_READ_LENGTH];
    snp SNPs[MAX_READ_LENGTH];
    
    uint8_t firstCase = 1;
    
    uint32_t prev_pos = 0;
    
    if(strcmp(edits, _readLength) == 0){
        // The read matches perfectly.
        compress_M(s->M, 1, P);
        return 0;
    }
    
    compress_M(s->M, 0, P);
    
    // The read does not match perfectly
    // Compute the edits
    while (*cigar != 0){
        if ( isdigit( *(cigar + i) ) == 0 )
            switch ( *(cigar + i) ){
                    // compute the position of the edit
                case 'M':
                    M += atoi(cigar);
                    
                    firstCase = 0;
                    cigar = cigar + i + 1;
                    i = -1;
                    break;
                    
                    // Store a insertion and all the previous snps
                case 'I':
                    I = atoi(cigar);
                    for (ctr = 0; ctr < I ; ctr++) {
                        pos = M;
                        if (lastSnp != 0)
                            lastSnp = add_snps_to_array(edits, SNPs, &numSnps, pos + numIns, read);
                        
                        Insers[numIns].pos = pos - prevPosI;
                        Insers[numIns].targetChar = char2basepair(read[pos+numIns]);
                        prevPosI = pos;
                        numIns++;
                    }
                    
                    firstCase = 0;
                    cigar = cigar + i + 1;
                    i = -1;
                    break;
                    
                    // Store the deletion
                case 'D':
                    D = atoi(cigar);
                    for (ctr = 0; ctr < D ; ctr++) {
                        pos = M;
                        Dels[numDels] = pos - prevPosD;
                        prevPosD = pos;
                        numDels++;
                    }
                    
                    firstCase = 0;
                    cigar = cigar + i + 1;
                    i = -1;
                    break;
                    
                case '*':
                    return 1;
                case 'S':
                    if (firstCase == 1) {
                        // S is the first thing we see
                        S = atoi(cigar);
                        for (ctrS = 0; ctrS < S; ctrS++){
                            if (lastSnp != 0)
                                lastSnp = add_snps_to_array(edits, SNPs, &numSnps, numIns, read);
                            Insers[numIns].pos = 0;
                            Insers[numIns].targetChar = char2basepair(read[ctrS]);
                            numIns++;
                        }
                    }else{
                        // We are at the end
                        S = atoi(cigar);
                        for (ctr = 0; ctr < S ; ctr++) {
                            pos = M;
                            Insers[numIns].pos = pos - prevPosI;
                            Insers[numIns].targetChar = char2basepair(read[pos+numIns]);
                            prevPosI = pos;
                            numIns++;
                        }
                    }
                    
                    firstCase = 0;
                    cigar = cigar + i + 1;
                    i = -1;
                    break;
                default:
                    break;
                    printf("Something besides MIDS appeared in the cigar\n");
                    printf("%c\n", *(cigar + i));
            }
        i++;
    }
    
    if (lastSnp != 0)
        add_snps_to_array(edits, SNPs, &numSnps, READ_LENGTH + 1, read);
    
    
    // Compress the edits
    if ((numDels | numIns) == 0) {
        compress_S(s->S, numSnps);
    }
    else{
        compress_S(s->S, 0);
        compress_I(s->I, numSnps);
        compress_I(s->I, numDels);
        compress_I(s->I, numIns);
    }
    
    // Store the positions and Chars in the corresponding vector
    prev_pos = 0;
    for (i = 0; i < numDels; i++){
        compress_pos(s->p, Dels[i], prev_pos, flag);
        prev_pos += Dels[i];
    }
    prev_pos = 0;
    for (i = 0; i < numSnps; i++){
        
        // compute delta to next snp
        delta = compute_delta_to_first_snp(prev_pos, READ_LENGTH);
        /*delta = READ_LENGTH + 2;
         for (j=0;j<READ_LENGTH - prev_pos; j++){
         if (snpInRef[cumsumP - 1 + j] == 1){
         delta = j;
         break;
         }
         }*/
        
        delta = (delta << BITS_DELTA);
        compress_pos(s->p, SNPs[i].pos, delta + prev_pos, flag);
        prev_pos += SNPs[i].pos + 1;
        snpInRef[cumsumP + prev_pos - 1 - 1] = 1;
        
        compress_chars(s->c, SNPs[i].refChar, SNPs[i].targetChar);
        
    }
    prev_pos = 0;
    for (i = 0; i < numIns; i++){
        compress_pos(s->p, Insers[i].pos, prev_pos, flag);
        prev_pos += Insers[i].pos;
        
        compress_chars(s->c, O, Insers[i].targetChar);
    }
    
    return cumsumP;
    
}

uint32_t reconstruct_read(sam_compressor s, uint32_t pos, uint8_t invFlag, uint8_t **read, FILE *fastqFile){
    
    unsigned int numIns = 0, numDels = 0, numSnps = 0, delPos = 0, ctrPos = 0, snpPos = 0, insPos = 0;
    uint32_t currentPos = 0, prevIns = 0, prev_pos = 0, delta = 0, deltaPos = 0;
    
    static uint32_t prevPos = 0;
    
    unsigned int ctrDels = 0;
    int i = 0, j = 0;
    
    uint8_t match;
    
    enum BASEPAIR refbp;
    
    if (pos < prevPos){
        deltaPos = pos;
    }else{
        deltaPos = pos - prevPos;
    }
    prevPos = pos;
    
    // The read matches perfectly.
    match = decompress_M(s->M, deltaPos);
    
    // cumsumP is equal to pos
    cumsumP = pos;
    
    // If there is a match, reconstruct the read
    if (match){
        
        switch (invFlag) {
            case 0:
                for (ctrPos=0; ctrPos<READ_LENGTH; ctrPos++)
                    fputc(reference[pos + ctrPos - 1],fastqFile);
                
                fprintf(fastqFile, "\n");
                break;
                
            case 1:
                for (ctrPos=0; ctrPos<READ_LENGTH; ctrPos++)
                    fputc(bp_complement( reference[pos + READ_LENGTH - 1 - ctrPos - 1] ),fastqFile);
                
                fprintf(fastqFile, "\n");
                break;
            default:
                printf("ERROR: invFlag must be 0 or 1\n");
                
        }
        
        
        return 1;
    }
    
    // There is no match, retreive the edits
    else{
        numSnps = decompress_S(s->S);
        
        
        if (numSnps == 0){
            numSnps = decompress_I(s->I);
            numDels = decompress_I(s->I);
            numIns = decompress_I(s->I);
        }
    }
    
    // Reconstruct the read
    
    // Deletions
    prev_pos = 0;
    for (ctrDels = 0; ctrDels < numDels; ctrDels++){
        
        delPos = decompress_pos(s->p, prev_pos, invFlag);
        prev_pos += delPos;
        
        // Do not take the deleted characters from the reference
        for (ctrPos = 0; ctrPos<delPos; ctrPos++){
            (*read)[currentPos] = reference[pos + currentPos - 1 + ctrDels];
            currentPos++;
        }
    }
    
    
    // Fill up the rest of the read up to numIns
    for (ctrPos = currentPos; ctrPos<READ_LENGTH - numIns; ctrPos++){
        (*read)[currentPos] = reference[pos + currentPos - 1 + numDels];
        currentPos++;
    }
    
    // SNPS
    currentPos = 0;
    prev_pos = 0;
    for (i = 0; i < numSnps; i++){
        
        // compute delta to next snp
        delta = compute_delta_to_first_snp(prev_pos, READ_LENGTH);
        delta = (delta << BITS_DELTA);
        
        snpPos = decompress_pos(s->p, delta + prev_pos, invFlag);
        prev_pos += snpPos + 1;
        snpInRef[cumsumP + prev_pos - 1 - 1] = 1;
        
        refbp = char2basepair( (*read)[currentPos + snpPos] );
        (*read)[currentPos + snpPos] = decompress_char(s->c, refbp);
        currentPos = currentPos + snpPos + 1;
        
    }
    currentPos = 0;
    
    // Insertions and write the read
    switch (invFlag) {
            // There is NO inversion
        case 0:
            // Write the insertions
            prev_pos = 0;
            for (i = 0; i < numIns; i++){
                
                insPos = decompress_pos(s->p, prev_pos, invFlag);
                prev_pos += insPos;
                
                // write the read up to the insertion
                for (ctrPos=0; ctrPos<insPos; ctrPos++)
                    fputc((*read)[currentPos],fastqFile), currentPos++;
                
                // write the insertion
                fputc(decompress_char(s->c, O), fastqFile);
            }
            
            // write the rest of the read
            for (ctrPos=currentPos; ctrPos<READ_LENGTH - numIns; ctrPos++)
                fputc((*read)[currentPos],fastqFile), currentPos++;
            
            fprintf(fastqFile, "\n");
            return 0;
            
            // There is an inversion
        case 1:
            prevIns = 0;
            prev_pos = 0;
            for (i = 0; i < numIns; i++){
                
                insPos = decompress_pos(s->p, prev_pos, invFlag);
                prev_pos += insPos;
                insPos += prevIns;
                // move the read one position to the left (make room for the insertion)
                for (ctrPos = READ_LENGTH - 1; ctrPos > insPos; ctrPos--)
                    (*read)[ctrPos] = (*read)[ctrPos - 1];
                // Add the insertion
                (*read)[insPos] = decompress_char(s->c, O);
                prevIns = insPos + 1;
            }
            
            // write the read inverted
            for (ctrPos=0; ctrPos<READ_LENGTH; ctrPos++)
                fputc(bp_complement((*read)[READ_LENGTH - 1 - ctrPos]),fastqFile);
            fprintf(fastqFile, "\n");
            
            return 1;
            
            
        default:
            printf("ERROR: invFLAG different from 0 and 1\n");
            return 2;
    }
}

int end_compression(sam_compressor s){
    
    encoder_last_step(s->F->a, s->F->os);
    
    encoder_last_step(s->P->a, s->P->os);
    fclose(s->P->os->fosA);
    
    encoder_last_step(s->M->a, s->M->os);
    
    encoder_last_step(s->S->a, s->S->os);
    
    encoder_last_step(s->I->a, s->I->os);
    
    encoder_last_step(s->p->a, s->p->os);
    
    encoder_last_step(s->c->a, s->c->os);
    
    
    return 1;
}

uint32_t decompress_chromosome(sam_compressor sComp, FILE *fout){
    
    unsigned int n = 0, prevPos = 0;
    
    uint32_t invFlag;
    int pos;
    uint8_t *read;
    
    read = (uint8_t*)calloc(READ_LENGTH, sizeof(uint8_t));
    
    while (1)
    {
        
        pos = decompress_P(sComp->P);
        
        if (pos == -1){
            free(read);
            return n;
        }
        
        if (pos == 0){
            free(read);
            return 0;
        }
        
        
        
        invFlag = decompress_F(sComp->F);
        
        reconstruct_read(sComp, pos, invFlag, &read, fout);
        
        prevPos = pos;
        
        n++;
    }
    
    
    printf("number of reads: %u\n", n);
    
    return n;
    
}

uint32_t start_sam_compression(FILE *fin, char* osPath){
    
    unsigned int n = 0, tempP = 0, tempF = 0;
    sam_compressor sComp;
    sam_line samLine;
    
    samLine = (sam_line) calloc(1, sizeof(struct sam_line_t));
    samLine->identifier = (char*) calloc(1, 4*MAX_READ_LENGTH);
    samLine->refname = (char*) calloc(1, 4*MAX_READ_LENGTH);
    samLine->cigar = (char*) calloc(1, 4*MAX_READ_LENGTH);
    samLine->edits = (char*) calloc(1, 4*MAX_READ_LENGTH);
    samLine->read = (char*) calloc(1, 4*MAX_READ_LENGTH);
    
    // Loop over the lines of the sam file
    n = 0;
    while ( read_line_from_sam(samLine, fin) )
    {
        // This is the first time we read from the SAM file
        if (n == 0){
            // Get the read_length
            while (samLine->read[n++] != 0)
                READ_LENGTH++;
            sprintf(_readLength, "%d", READ_LENGTH);
            n = 0;
            
            // Initialize the compressor
            sComp = initialize_sam_compressor(osPath, COMPRESSION, &READ_LENGTH);
            
            // Write READ_LENGTH to P->os->fosA
            fprintf(sComp->P->os->fosA, "%d\n", READ_LENGTH);
            
        }
        
        if (n%1000000 == 0) {
            printf("%d Million reads compressed.\n", n/1000000);
        }
        
        // Compress sam line
        tempF = compress_F(sComp->F, samLine->invFlag);
        tempP = compress_P(sComp->P, samLine->pos);
        compress_edits(sComp, samLine->edits, samLine->cigar, samLine->read, tempP, tempF);
        
        n++;
    }
    compress_P(sComp->P, 0);
    end_compression(sComp);
    
    printf("Number of reads compressed: %u\n", n);
    
    // free(samLine->cigar), free(samLine.edits), free(samLine.read_), free(samLine.identifier), free(samLine.refname);
    fclose(fin);
    
    return n;
    
}

uint32_t start_sam_decompression(FILE* fout, char* osPath, FILE* fref){
    
    char refPath[MAX_PATH_LENGTH];
    sam_compressor samDecomp;
    FILE* fchr;
    
    uint32_t n = 0, chrCtr = 0, cumctr = 0, i = 0;
    
    samDecomp = initialize_sam_compressor(osPath, DECOMPRESSION, &READ_LENGTH);
    
    // fscanf(samDecomp->P->os->fosA, "%d", &READ_LENGTH);
    
    while (1) {
        
        // Read path of Reference from refFile
        fscanf(fref, "%s", refPath);
        
        // Open the Ref file
        fchr = fopen ( refPath , "r" );
        if (fchr==NULL) {fputs ("Chromosome (ref) File error\n",stderr); exit (1);}
        
        // Store Ref sequence in memory
        store_reference_in_memory(fchr);
        fclose (fchr);
        
        printf("reference stored\n");
        
        // Clean snpInRef vector and reset cumsumP
        cumsumP = 0;
        for (i=0; i<MAX_BP_CHR; i++){
            snpInRef[i] = 0;
        }
        
        n = decompress_chromosome(samDecomp, fout);
        cumctr += n;
        
        if (n == 0)
            break;
        
        printf("Chromosome %d decompressed (%d reads).\n", ++chrCtr, n);
        
    }
    
    return cumctr;
    
    
}

///////////////////////////////
int help(){
  
  printf("To compress:\n");
  printf("./main.run samFile outputPrefix\n");
  printf("To decompress:\n");
  printf("./main.run readsFileForOutput reference.txt outputPrefix\n");
  printf("The reference.txt is a file that contains the path to chrX in the Xth line\n");
  printf("see README for more detailed information on how to run the program\n");

  return 1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
//                                                                                                          //
//                                         MAIN                                                             //
//                                                                                                          //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char * argv[])
{
    FILE *fin, *fout, *fref;
    
    char osPath[MAX_PATH_LENGTH];
    
    
    
    
    uint32_t numReads = 0;
    
    clock_t begin;
    clock_t end;
    
    // Check that the number of input parameters is correct before open the files!
    if (argc == 3){
        
        fin = fopen( argv[1] , "r");
        
        strcpy(osPath, argv[2]);
        
        // Start the compression
        begin = clock();
        numReads = start_sam_compression(fin, osPath);
        end = clock();
        
        printf("Time compressing the SAM file: %.02f seconds\n\n", (double)(end - begin)/ CLOCKS_PER_SEC);
        
    }
    
    else if (argc == 4)
    {
        
        fout = fopen(argv[1], "w");
        fref = fopen(argv[2], "r");
        
        strcpy(osPath, argv[3]);
        
        begin = clock();
        start_sam_decompression(fout, osPath, fref);
        end = clock();
        
        printf("Time decompressing the SAM file: %.02f seconds\n\n", (double)(end - begin)/ CLOCKS_PER_SEC);
    }else{
      help();
      return 1;
    }
    
    
    
    return 0;
}

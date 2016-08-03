//
//  read_decompression.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/13/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include "read_compression.h"
#include <ctype.h>
#include <stdint.h>

#define DEBUG true
//**************************************************************//
//                                                              //
//                  STORE REFERENCE IN MEMORY                   //
//                                                              //
//**************************************************************//
int store_reference_in_memory(FILE* refFile){
    uint32_t letterCount, endoffile = 1;
    char header[1024];
    
    reference = (char *) malloc(MAX_BP_CHR*sizeof(char));
    
    // ******* Read and Store Reference****** //
    letterCount = 0;
    
    // Remove the header
    fgets(header, sizeof(header), refFile);
    
    while (fgets(&reference[letterCount], 1024, refFile))
    {
        
        if(reference[letterCount] == '>' || reference[letterCount] == '@'){
            endoffile = 0;
            break;
        }
       
        reference[letterCount] = toupper(reference[letterCount]);
        while (reference[letterCount++] != '\n' ) { 
          reference[letterCount] = toupper(reference[letterCount]);
        }
        letterCount--;
        
    }
        
    reference[letterCount] = '\0';
    
    reference = (char *) realloc(reference, letterCount);
    
    if (endoffile)
        return END_GENOME_FLAG;
    
    return letterCount;
    
}


/************************
 * Decompress the read
 **********************/
uint32_t decompress_read(Arithmetic_stream as, sam_block sb, uint8_t chr_change, struct sam_line_t *sline){
    
    int invFlag, tempP, k;
    uint32_t readLen;
    uint16_t maskedReadVal;
    
    read_models models = sb->reads->models;
    
    
    // Decompress read length (4 uint8)
    readLen = 0;
    for (k=0;k<4;k++) {
        maskedReadVal = decompress_uint8t(as, models->rlength[k]);
        readLen = readLen | maskedReadVal<<(k*8);
    }
    
    
    
    // Decompress the read
    tempP = decompress_pos(as, models->pos, models->pos_alpha, chr_change, &sline->pos);
    
    invFlag = decompress_flag(as, models->flag, &sline->flag);
    
    reconstruct_read(as, models, tempP, invFlag, sline->read, readLen, sline->cigar);
    
    return invFlag;
}


/************************
 * Decompress the cigar
 **********************/
uint32_t decompress_cigar(Arithmetic_stream as, sam_block sb, struct sam_line_t *sline)
{
    uint8_t cigarFlags;
    int cigarCtr,cigarLen;
    
    read_models models = sb->reads->models;
    
    
    // Decompress cigarFlags
    cigarFlags = decompress_uint8t(as,models->cigarFlags[0]);
    
    if(cigarFlags==0) {
        //The cigar is not the recCigar.
        cigarLen = decompress_uint8t(as,models->cigar[0]);
        for(cigarCtr = 0; cigarCtr<cigarLen; cigarCtr++)
            sline->cigar[cigarCtr] = decompress_uint8t(as,models->cigar[0]);
        
        sline->cigar[cigarCtr] = '\0';
    }
    
    return 1;
}



/***********************
 * Decompress the Flag
 **********************/
uint32_t decompress_flag(Arithmetic_stream a, stream_model *F, uint32_t *flag){
    
    
    // In this case we are just compressing the binary information of whether the
    // read is in reverse or not. we use F[0] as there is no context for the flag.
    uint16_t x;
    // Read the value from the Arithmetic Stream
    x = read_value_from_as(a, F[0]);
    
    // Update model
    update_model(F[0], x);
    
    *flag = x;
    
    x = x & 16;
    x >>= 4;
    
    return x;
    
}

/***********************************
 * Decompress the Alphabet of Position
 ***********************************/
uint32_t decompress_pos_alpha(Arithmetic_stream as, stream_model *PA){
    
    uint32_t Byte = 0, x = 0;
    
    // we encode byte per byte i.e. x = [B0 B1 B2 B3]
    
    // Read B0 from the Arithmetic Stream using the alphabet model
    Byte = read_value_from_as(as, PA[0]);
    // Update model
    update_model(PA[0], Byte);
    // Reconstruct B0
    x |= (Byte << 24);
    
    // Read B1 from the Arithmetic Stream using the alphabet model
    Byte = read_value_from_as(as, PA[1]);
    // Update model
    update_model(PA[1], Byte);
    // Reconstruct B1
    x |= (Byte << 16);
    
    // Send B2 to the Arithmetic Stream using the alphabet model
    Byte = read_value_from_as(as, PA[2]);
    // Update model
    update_model(PA[2], Byte);
    // Reconstruct B2
    x |= (Byte << 8);
    
    // Send B3 to the Arithmetic Stream using the alphabet model
    Byte = read_value_from_as(as, PA[3]);
    // Update model
    update_model(PA[3], Byte);
    // Reconstruct B3
    x |= (Byte);
    
    return x;
    
    
}

/**************************
 * Decompress the Position
 *************************/
uint32_t decompress_pos(Arithmetic_stream as, stream_model *P, stream_model *PA, uint8_t chr_change, uint32_t *p){
    
    static uint32_t prevPos = 0;
    
    int32_t pos, alphaMapX = 0, x = 0;
    
    enum {SMALL_STEP = 0, BIG_STEP = 1};
    
    // Check if we are changing chromosomes.
    if (chr_change)
        prevPos = 0;
    
    // Read from the AS and get the position
    alphaMapX = read_value_from_as(as, P[0]);
    
    x = P[0]->alphabet[alphaMapX];
    
    // Update the statistics
    update_model(P[0], alphaMapX);
    
    // A new value of pos
    if (x == -1) {
        
        // Read from the AS to get the unknown alphabet letter alpha
        x = decompress_pos_alpha(as, PA);
        
        // Update the statistics of the alphabet for x
        P[0]->alphaExist[x] = 1;
        P[0]->alphaMap[x] = P[0]->alphabetCard; // We reserve the bin 0 for the new symbol flag
        P[0]->alphabet[P[0]->alphabetCard] = x;
        
        update_model(P[0], P[0]->alphabetCard++);
    }
    
    // Decompress the position diference (+ 1 to reserve 0 for new symbols)
    pos = prevPos + x - 1;
    
    *p = pos;
    
    prevPos = pos;
    
    return pos;
}

/****************************
 * Decompress the match
 *****************************/
uint32_t decompress_match(Arithmetic_stream a, stream_model *M, uint32_t P){
    
    uint32_t ctx = 0;
    static uint8_t  prevM = 0;
    
    uint8_t match = 0;
    
    
    // Compute Context
    P = (P != 1)? 0:1;
    //prevP = (prevP > READ_LENGTH)? READ_LENGTH:prevP;
    //prevP = (prevP > READ_LENGTH/4)? READ_LENGTH:prevP;
    
    ctx = (P << 1) | prevM;
    
    //ctx = 0;
    
    // Read the value from the Arithmetic Stream
    match = read_value_from_as(a, M[ctx]);
    
    // Update model
    update_model(M[ctx], match);
    
    prevM = match;
    
    return match;
}

/*************************
 * Decompress the snps
 *************************/
uint32_t decompress_snps(Arithmetic_stream a, stream_model *S){
    
    uint8_t numSnps = 0;
    // No context is used for the numSnps for the moment.
    
    // Send the value to the Arithmetic Stream
    numSnps = read_value_from_as(a, S[0]);
    
    // Update model
    update_model(S[0], numSnps);
    
    return numSnps;
    
}


/********************************
 * Decompress the indels
 *******************************/
uint32_t decompress_indels(Arithmetic_stream a, stream_model *I){
    
    uint8_t numIndels = 0;
    // No context is used for the numIndels for the moment.
    
    // Read the value from the Arithmetic Stream
    numIndels = read_value_from_as(a, I[0]);
    
    // Update model
    update_model(I[0], numIndels);
    
    return numIndels;
    
}

/*******************************
 * Decompress the variations
 ********************************/
uint32_t decompress_var(Arithmetic_stream a, stream_model *v,  uint32_t prevPos, uint32_t flag){
    
    uint32_t ctx = 0;
    uint32_t pos = 0;
    
    //flag = 0;
    ctx = prevPos << 1 | flag;
    
    // Read the value from the Arithmetic Stream
    pos = read_value_from_as(a, v[ctx]);
    
    // Update model
    update_model(v[ctx], pos);
    
    return pos;
    
}

/*****************************************
 * Decompress the chars
 ******************************************/
uint8_t decompress_chars(Arithmetic_stream a, stream_model *c, enum BASEPAIR ref){
    
    uint32_t target = 0;
    
    // Read the value from the Arithmetic Stream
    target = read_value_from_as(a, c[ref]);
    
    // Update model
    update_model(c[ref], target);
    
    return basepair2char((enum BASEPAIR)target);
    
}

uint32_t argmin(uint32_t *arr, uint32_t len) {
  uint32_t min = UINT32_MAX;
  uint32_t index = 0;
  for (uint32_t i = 0; i < len; i++) {
    if (min > arr[i]) {
      min = arr[i]; 
      index = i;
    }
  }
  return index;
}

/*****************************************
 * Reconstruct the read
 ******************************************/
uint32_t reconstruct_read(Arithmetic_stream as, read_models models, uint32_t pos, uint8_t invFlag, char *read, uint32_t readLen, char* recCigar){
    
    unsigned int numIns = 0, numDels = 0, numSnps = 0, delPos = 0, ctrPos = 0, snpPos = 0, insPos = 0;
    uint32_t currentPos = 0, prevIns = 0, prev_pos = 0, delta = 0, deltaPos = 0;

    
    uint32_t Dels[MAX_READ_LENGTH];
    ins Insers[MAX_READ_LENGTH];
    snp SNPs[MAX_READ_LENGTH];
    
    static uint32_t prevPos = 0;
    
    unsigned int ctrDels = 0, readCtr = 0;
    int i = 0;
    uint8_t tmpChar;
    unsigned int returnVal;
    
    uint8_t match;
    
    enum BASEPAIR refbp;
    
    if (pos < prevPos){
        deltaPos = pos;
    }else{
        deltaPos = pos - prevPos + 1;// deltaPos is 1-based.
    }
    prevPos = pos;
    
    // The read matches perfectly.
    match = decompress_match(as, models->match, deltaPos);
    
    // cumsumP is equal to pos
    cumsumP = pos;
    
    // If there is a match, reconstruct the read
    if (match) {
      for (ctrPos=0; ctrPos<models->read_length; ctrPos++)
        read[readCtr++] = reference[pos + ctrPos - 1];
        return 1;
    }
    // There is no match, retreive the edits
    else{
        numSnps = decompress_snps(as, models->snps);
        
        
        if (numSnps == 0){
            numSnps = decompress_indels(as, models->indels);
            numDels = decompress_indels(as, models->indels);
            numIns = decompress_indels(as, models->indels);
        } else {
            numDels = 0;
            numIns = 0;
        }
    }
    
    //printf("snps %d, dels %d, ins %d\n", numSnps, numDels, numIns);

    // Reconstruct the read
    
    struct sequence seq;
    init_sequence(&seq, Dels, Insers, SNPs);
    seq.n_ins = numIns;
    seq.n_snps = numSnps;
    seq.n_dels = numDels;

    char *tempRead = (char*)alloca(models->read_length*sizeof(char) + 2);
    // Deletions
    prev_pos = 0;
    for (ctrDels = 0; ctrDels < numDels; ctrDels++){
        delPos = decompress_var(as, models->var, prev_pos, invFlag);
        //printf("Delete ref at %d, prev %d\n", delPos, prev_pos);
        Dels[ctrDels] = delPos + prev_pos;
        prev_pos += delPos;
    }

    currentPos = 0;
    
    prev_pos = 0;
    for (i = 0; i < numIns; i++){
        insPos = decompress_var(as, models->var, prev_pos, invFlag);
        Insers[i].pos = prev_pos + insPos;
        Insers[i].targetChar = char2basepair(decompress_chars(as, models->chars, O));
        //printf("Insert %c at offset %d, prev_pos %d\n", basepair2char(Insers[i].targetChar), insPos, prev_pos);
        prev_pos += insPos;
    }

    // SNPS
    prev_pos = 0;
    uint32_t ref_pos = 0;
    for (i = 0; i < numSnps; i++){
        
        assert(prev_pos <  models->read_length);
        snpPos = decompress_var(as, models->var, prev_pos, invFlag);

        // Get adjusted ref_pos
        // TODO FIX THIS. it's buggy
        ref_pos = snpPos + prev_pos; 

        int add = 0;
        for (int i = 0; i < numDels; i++) {
          if (Dels[i] <= ref_pos) add++;
        }
        int sub = 0;
        for (int i = numIns - 1; i >= 0; i--) {
          if (Insers[i].pos <= ref_pos) sub++;
        }
        ref_pos = ref_pos + add - sub;

        refbp = char2basepair( reference[pos + ref_pos - 1] );
        SNPs[i].refChar = refbp;
        SNPs[i].targetChar = char2basepair(decompress_chars(as, models->chars, refbp));
        SNPs[i].pos = prev_pos + snpPos;

        //printf("Replace %c with %c, Reference: %d, offset: %d, prev_pos = %d\n", basepair2char(refbp), basepair2char(SNPs[i].targetChar), ref_pos, snpPos, prev_pos);
        prev_pos += snpPos;

    }
    reconstruct_read_from_ops(&seq, &(reference[pos - 1]), read, models->read_length);

    if (invFlag != 1 && invFlag != 0) {
      printf("ERROR: invFLAG different from 0 and 1\n");
      return 2;
    } 
    /*
    if (invFlag == 1) {
      char tmp;
      for (int i = 0; i < models->read_length / 2; i++) {
        int k = models->read_length - i - 1;
        tmp = read[k];
        read[k] = bp_complement(read[i]);
        read[i] = bp_complement(tmp);
      }
    }*/
    return returnVal;
}

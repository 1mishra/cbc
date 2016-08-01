//
//  reads_compression.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/5/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//


#include "read_compression.h"

/************************
 * Compress the read
 **********************/
uint32_t compress_read(Arithmetic_stream as, read_models models, read_line samLine, uint8_t chr_change){
    
    int tempF, PosDiff, chrPos, k;
    uint32_t mask;
    uint16_t maskedReadVal;
    // For now lets skip the unmapping ones
    if (samLine->invFlag & 4) {
//        assert(strcmp("*", samLine->cigar) == 0);
        return 1;
    }
    // compress read length (assume int)
    for (k=0;k<4;k++) {
        mask = 0xFF<<(k*8);
        maskedReadVal = (uint8_t)(models->read_length & mask)>>(k*8);
        compress_uint8t(as, models->rlength[k], maskedReadVal);
    }
    
    // Compress sam line
    PosDiff = compress_pos(as, models->pos, models->pos_alpha, samLine->pos, chr_change);
    tempF = compress_flag(as, models->flag, samLine->invFlag);
    //tempF = compress_flag(as, models->flag, 0);
    chrPos = compress_edits(as, models, samLine->edits, samLine->cigar, samLine->read, PosDiff, tempF, &(samLine->cigarFlags));
    
    assert(samLine->pos  == chrPos);

    return 1;
}


/***********************
 * Compress the Flag
 **********************/
uint32_t compress_flag(Arithmetic_stream a, stream_model *F, uint16_t flag){
    
    
    // In this case we need to compress the whole flag, althugh the binary information of whether the
    // read is in reverse or not is the most important. Thus, we return the binary info.
    //we use F[0] as there is no context for the flag.
    
    uint16_t x = 0;
    
    x = flag << 11;
    x >>= 15;
    
    // Send the value to the Arithmetic Stream
    send_value_to_as(a, F[0], flag);
    
    // Update model
    update_model(F[0], flag);
    
    return x;
    
}

/***********************************
 * Compress the Alphabet of Position
 ***********************************/
uint32_t compress_pos_alpha(Arithmetic_stream as, stream_model *PA, uint32_t x){
    
    uint32_t Byte = 0;
    
    // we encode byte per byte i.e. x = [B0 B1 B2 B3]
    
    // Send B0 to the Arithmetic Stream using the alphabet model
    Byte = x >> 24;
    send_value_to_as(as, PA[0], Byte);
    // Update model
    update_model(PA[0], Byte);
    
    // Send B1 to the Arithmetic Stream using the alphabet model
    Byte = (x & 0x00ff0000) >> 16;
    send_value_to_as(as, PA[1], Byte);
    // Update model
    update_model(PA[1], Byte);
    
    // Send B2 to the Arithmetic Stream using the alphabet model
    Byte = (x & 0x0000ff00) >> 8;
    send_value_to_as(as, PA[2], Byte);
    // Update model
    update_model(PA[2], Byte);
    
    // Send B3 to the Arithmetic Stream using the alphabet model
    Byte = (x & 0x000000ff);
    send_value_to_as(as, PA[3], Byte);
    // Update model
    update_model(PA[3], Byte);
    
    return 1;
    
    
}

/***********************
 * Compress the Position
 **********************/
uint32_t compress_pos(Arithmetic_stream as, stream_model *P, stream_model *PA, uint32_t pos, uint8_t chr_change){
    
    static uint32_t prevPos = 0;
    enum {SMALL_STEP = 0, BIG_STEP = 1};
    int32_t x = 0;
    
    // TODO diferent update models for updating -1 and already seen symbols
    // i.e., SMALL_STEP and BIG_STEP
    
    // Check if we are changing chromosomes.
    if (chr_change)
        prevPos = 0;
    
    
    // Compress the position diference (+ 1 to reserve 0 for new symbols)
    x = pos - prevPos + 1;
    
    if (P[0]->alphaExist[x]){
        // Send x to the Arithmetic Stream
        send_value_to_as(as, P[0], P[0]->alphaMap[x]);
        // Update model
        update_model(P[0], P[0]->alphaMap[x]);
    }
    else{
        
        // Send 0 to the Arithmetic Stream
        send_value_to_as(as, P[0], 0);
        
        // Update model
        update_model(P[0], 0);
        
        // Send the new letter to the Arithmetic Stream using the alphabet model
        compress_pos_alpha(as, PA, x);
        
        // Update the statistics of the alphabet for x
        P[0]->alphaExist[x] = 1;
        P[0]->alphaMap[x] = P[0]->alphabetCard; // We reserve the bin 0 for the new symbol flag
        P[0]->alphabet[P[0]->alphabetCard] = x;
        
        // Update model
        update_model(P[0], P[0]->alphabetCard++);
    }
    
    prevPos = pos;
    
    return x;
}

/****************************
 * Compress the match
 *****************************/
uint32_t compress_match(Arithmetic_stream a, stream_model *M, uint8_t match, uint32_t P){
    
    uint32_t ctx = 0;
    static uint8_t  prevM = 0;
    
    
    // Compute Context
    P = (P != 1)? 0:1;
    //prevP = (prevP > READ_LENGTH)? READ_LENGTH:prevP;
    //prevP = (prevP > READ_LENGTH/4)? READ_LENGTH:prevP;
    
    ctx = (P << 1) | prevM;
    
    //ctx = 0;
    
    // Send the value to the Arithmetic Stream
    send_value_to_as(a, M[ctx], match);
    
    // Update model
    update_model(M[ctx], match);
    
    prevM = match;
    
    return 1;
}

/*************************
 * Compress the snps
 *************************/
uint32_t compress_snps(Arithmetic_stream a, stream_model *S, uint8_t numSnps){
    
    
    // No context is used for the numSnps for the moment.
    
    // Send the value to the Arithmetic Stream
    send_value_to_as(a, S[0], numSnps);
    
    // Update model
    update_model(S[0], numSnps);
    
    return 1;
    
}


/********************************
 * Compress the indels
 *******************************/
uint32_t compress_indels(Arithmetic_stream a, stream_model *I, uint8_t numIndels){
    
    
    // Nos context is used for the numIndels for the moment.
    
    // Send the value to the Arithmetic Stream
    send_value_to_as(a, I[0], numIndels);
    
    // Update model
    update_model(I[0], numIndels);
    
    return 1;
    
}

/*******************************
 * Compress the variations
 ********************************/
uint32_t compress_var(Arithmetic_stream a, stream_model *v, uint32_t pos, uint32_t prevPos, uint32_t flag){
    
    uint32_t ctx = 0;
    
    //flag = 0;
    ctx = prevPos << 1 | flag;
    
    // Send the value to the Arithmetic Stream
    send_value_to_as(a, v[ctx], pos);
    
    // Update model
    update_model(v[ctx], pos);
    
    return 1;
    
}

/*****************************************
 * Compress the chars
 ******************************************/
uint32_t compress_chars(Arithmetic_stream a, stream_model *c, enum BASEPAIR ref, enum BASEPAIR target){
    
    // Send the value to the Arithmetic Stream
    send_value_to_as(a, c[ref], target);
    
    // Update model
    update_model(c[ref], target);
    
    return 1;
    
}

/*****************************************
 * Compress the edits
 ******************************************/
uint32_t compress_edits(Arithmetic_stream as, read_models rs, char *edits, char *cigar, char *read, uint32_t deltaP, uint8_t flag, uint8_t* cigarFlags){
    
    unsigned int numIns = 0, numDels = 0, numSnps = 0, lastSnp = 1;
    int i = 0, tmpi = 0, M = 0, I = 0, D = 0, pos = 0, ctr = 0, prevPosI = 0, prevPosD = 0, ctrS = 0, S = 0;
    uint32_t delta = 0;
    char recCigar[MAX_CIGAR_LENGTH];

    
    //uint32_t k, tempValue, tempSum = 0;
    
    uint32_t posRef, posRead, tmpM, tmpI, tmpD, tmpS, match;
    char *tmpcigar, *tmpEdits;
    char *origCigar = cigar;
    
    
    // pos in the reference
    cumsumP = cumsumP + deltaP - 1;// DeltaP is 1-based
    
    uint32_t Dels[MAX_READ_LENGTH];
    ins Insers[MAX_READ_LENGTH];
    snp SNPs[MAX_READ_LENGTH];
    
    uint8_t firstCase = 1;
    
    uint32_t prev_pos = 0;
    
    
    //ALERTA AQUI y en su analogo descomp.: si pasamos por aqui no hacemos nada con el cigar... (arreglado ya?)
    if(strcmp(edits, rs->_readLength) == 0){
        // The read matches perfectly.
        compress_match(as, rs->match, 1, deltaP);
        
        return cumsumP;
    }
    
    compress_match(as, rs->match, 0, deltaP);
    
    if (strcmp("0A0T0A0A0A96", edits)==0) {
        ;
    }
    
    if (lastSnp != 0)
        add_snps_to_array(edits, SNPs, &numSnps, rs->read_length + 1, read);
    
    
    // Compress the edits
    if ((numDels | numIns) == 0) {
        compress_snps(as, rs->snps, numSnps);
    }
    else{
        compress_snps(as, rs->snps, 0);
        compress_indels(as, rs->indels, numSnps);
        compress_indels(as, rs->indels, numDels);
        compress_indels(as, rs->indels, numIns);
    }
    
    // Store the positions and Chars in the corresponding vector
    prev_pos = 0;
    for (i = 0; i < numDels; i++){
        compress_var(as, rs->var, Dels[i], prev_pos, flag);
        prev_pos += Dels[i];
    }
    prev_pos = 0;
    for (i = 0; i < numSnps; i++){
        
        // compute delta to next snp
        delta = compute_delta_to_first_snp(prev_pos, rs->read_length);
        /*delta = READ_LENGTH + 2;
         for (j=0;j<READ_LENGTH - prev_pos; j++){
         if (snpInRef[cumsumP - 1 + j] == 1){
         delta = j;
         break;
         }
         }*/
        
        delta = (delta << BITS_DELTA);
        compress_var(as, rs->var, SNPs[i].pos, delta + prev_pos, flag);
        prev_pos += SNPs[i].pos + 1;
        snpInRef[cumsumP + prev_pos - 1 - 1] = 1;
        
        compress_chars(as, rs->chars, SNPs[i].refChar, SNPs[i].targetChar);
        
    }
    prev_pos = 0;
    for (i = 0; i < numIns; i++){
        compress_var(as, rs->var, Insers[i].pos, prev_pos, flag);
        prev_pos += Insers[i].pos;
        
        compress_chars(as, rs->chars, O, Insers[i].targetChar);
    }
    
    

    return cumsumP;
    
}



/******************************************
 * Function to look for snps in the cigar
 ****************************************/
int add_snps_to_array(char* edits, snp* SNPs, unsigned int *numSnps, unsigned int insertionPos, char *read){
    
    static unsigned int prevEditPtr = 0, cumPos = 0;
    
    int pos = 0, tempPos = 0, ctr;
    char ch = 0;
    
    uint8_t flag = 0;
    
    edits += prevEditPtr;
    
    while (*edits != 0 ) {
        
        pos = atoi(edits);
        
        tempPos = pos;
        
        ctr = compute_num_digits(pos);
        ch = *(edits+ctr);
        ctr++;
        
        // if there are deletions after pos, we need to add those positions that come after the deletions
        while (ch == '^'){
            while (isdigit(*(edits+ctr)) == 0) {
                ctr++;
            }
            tempPos += atoi(edits + ctr);
            
            ctr += compute_num_digits(atoi(edits + ctr));
            
            ch = *(edits+ctr);
            ctr++;
            if (ch == '\0'){
                flag = 1;
                break;
            }
        }
        
        if (flag == 1){
            break;
        }
        
        
        if (cumPos + tempPos >= insertionPos){
            cumPos++;
            return cumPos;
        }
        
        tempPos = atoi(edits);
        edits += compute_num_digits(tempPos);
        prevEditPtr += compute_num_digits(tempPos);
        
        //if ( isnumber(*(edits+1)) ) edits += 2, prevEditPtr += 2;
        //else edits += 1, prevEditPtr += 1;
        
        ch = *edits++, prevEditPtr++;
        
        while (ch == '^'){
            while (isdigit(*edits) == 0)
                edits++, prevEditPtr++;
            pos += atoi(edits);
            
            tempPos = atoi(edits);
            edits += compute_num_digits(tempPos);
            prevEditPtr += compute_num_digits(tempPos);
            
            //if ( isnumber(*(edits+1)) ) edits += 2, prevEditPtr += 2;
            //else edits += 1, prevEditPtr += 1;
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

uint32_t compute_delta_to_first_snp(uint32_t prevPos, uint32_t readLen){
    
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

uint32_t compute_num_digits(uint32_t x){
    
    //Get the number of digits (We assume readLength < 1000)
    
    if (x < 10)
        return 1;
    else if (x < 100)
        return 2;
    else if (x < 1000)
        return 3;
    else if (x < 10000)
        return 4;
    else if (x < 100000)
        return 5;
    else if (x < 1000000)
        return 6;
    else if (x < 10000000)
        return 7;
    else if (x < 100000000)
        return 8;
    else
        return 9;
    
}



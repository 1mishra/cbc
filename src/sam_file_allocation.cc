//
//  sam_file_allocation.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/18/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include "sam_block.h"
#include <sstream>
#include <string>
#include <iostream>

using namespace std;

uint32_t get_read_length(FILE *f){
    
    int ch, header_bytes = 0;
    char buffer[2048]; // 2 KB buffer
    
    // We use this opportunity to remove the headers
    // getc coge caracter, fgets coge 2048 o hasta salto linea.
    // con esto simplemente calculamos los bytes de header (empiezan con @),
    // el contenido de buffer nos da igual (se override cada vez).
    while ((ch = getc(f)) == '@') {
        fgets(buffer, 2048, f);
        header_bytes += strlen(buffer) + 1; // +1 to account for the @
    }
    
    // rewind the file pointer to be at the beginning of the first read
    fseek(f, header_bytes, SEEK_SET);
    
    // Extract the first read
    fscanf(f, "%*s %*d %*s %*d %*d %*s %*s %*d %*d %s", buffer);
    
    // rewind the file pointer to be at the beginning of the first read
    fseek(f, header_bytes, SEEK_SET);
    
    
    return (uint32_t)strlen(buffer);
    
}

/**
 *
 */
rname_block alloc_rname_block(){
    uint32_t i = 0;
    
    rname_block rtn = (rname_block) calloc(1, sizeof(struct rname_block_t));
    
    rtn->block_length = MAX_LINES_PER_BLOCK;
    
    rtn->rnames = (char**) calloc(rtn->block_length, sizeof(char*));
    
    // allocate the memory for each of the lines
    for (i = 0; i < rtn->block_length; i++) {
        
        rtn->rnames[i] = (char*) calloc(MAX_READ_LENGTH, sizeof(char));
    }
    
    // Allocate (and initialize) the models for the rnames
    rtn->models = alloc_rname_models_t();
    
    return rtn;
}

/**
 *
 */
mapq_block alloc_mapq_block(){
    
    mapq_block rtn = (mapq_block) calloc(1, sizeof(struct mapq_block_t));
    
    rtn->block_length = 1;
    
    rtn->mapq = (char*) calloc(rtn->block_length, sizeof(char));
    
    // Allocate (and initialize) the models for the rnames
    rtn->models = alloc_mapq_models_t();
    
    return rtn;
}



/**
 *
 */
rnext_block alloc_rnext_block(){
    uint32_t i = 0;
    
    rnext_block rtn = (rnext_block) calloc(1, sizeof(struct rnext_block_t));
    
    rtn->block_length = 1;
    
    rtn->rnext = (char**) calloc(rtn->block_length, sizeof(char*));
    
    // allocate the memory for each of the lines
    for (i = 0; i < rtn->block_length; i++) {
        
        rtn->rnext[i] = (char*) calloc(MAX_READ_LENGTH, sizeof(char));
    }
    
    // Allocate (and initialize) the models for the rnames
    rtn->models = alloc_rnext_models_t();
    
    return rtn;
}

/**
 *
 */
pnext_block alloc_pnext_block(){
    
    pnext_block rtn = (pnext_block) calloc(1, sizeof(struct pnext_block_t));
    
    rtn->block_length = 1;
    
    rtn->pnext = (uint32_t *) calloc(rtn->block_length, sizeof(uint32_t*));
    
    // Allocate (and initialize) the models for the rnames
    rtn->models = alloc_pnext_models_t();
    
    return rtn;
}

/**
 *
 */
tlen_block alloc_tlen_block(){
    
    tlen_block rtn = (tlen_block) calloc(1, sizeof(struct tlen_block_t));
    
    rtn->block_length = 1;
    
    rtn->tlen = (int32_t *) calloc(rtn->block_length, sizeof(int32_t*));
    
    // Allocate (and initialize) the models for the rnames
    rtn->models = alloc_tlen_models_t();
    
    return rtn;
}



/**
 *
 */
id_block alloc_id_block(){
    
    uint32_t i = 0;
    
    id_block rtn = (id_block) calloc(1, sizeof(struct id_block_t));
    
    rtn->block_length = 1;
    
    rtn->IDs = (char**) calloc(rtn->block_length, sizeof(char*));
    
    // allocate the memory for each of the lines
    for (i = 0; i < rtn->block_length; i++) {
        
        rtn->IDs[i] = (char*) calloc(MAX_READ_LENGTH, sizeof(char));
    }
    
    // Allocate (and initialize) the models for the IDs
    rtn->models = alloc_id_models_t();
    
    return rtn;
}

/**
 *
 */
aux_block alloc_aux_block() {
    
    uint32_t i = 0;
    
    aux_block rtn = (aux_block) calloc(1, sizeof(struct aux_block_t));
    
    rtn->block_length = 10;
    
    rtn->aux_str = (char**) calloc(MAX_AUX_FIELDS, sizeof(char*));
    
    // allocate the memory for each of the lines
    for (i = 0; i < MAX_AUX_FIELDS; i++) {
        rtn->aux_str[i] = (char*) calloc(MAX_AUX_LENGTH, sizeof(char));
    }
    
    rtn->most_common = (char**) calloc(MOST_COMMON_LIST_SIZE, sizeof(char*));
    // allocate the memory for each of the lines
    for (i = 0; i < MOST_COMMON_LIST_SIZE; i++) {
        rtn->most_common[i] = (char*) calloc(MAX_AUX_LENGTH, sizeof(char));
    }

    // Allocate (and initialize) the models for the aux
    rtn->models = alloc_aux_models_t();
    
    return rtn;
}

/**
 *
 */
read_block alloc_read_block_t(uint32_t read_length){
    
    read_block rf = (read_block) calloc(1, sizeof(struct read_block_t));
    
    rf->block_length = 1;
    
    rf->lines = (read_line) calloc(rf->block_length, sizeof(struct read_line_t));
        
    rf->lines->cigar = (char*) calloc(1, 2*read_length);
    rf->lines->edits = (char*) calloc(1, 2*read_length);
    rf->lines->read = (char*) calloc(1, read_length + 3);
    
    // Allocate (and initialize) the models for the reads
    rf->models = alloc_read_models_t(read_length);
    
    return rf;
}

sam_block alloc_sam_models(Arithmetic_stream as, FILE * fin, FILE *fref, struct qv_options_t *qv_opts, uint8_t mode){
    
    uint32_t i = 0;
    
    sam_block sb = (sam_block) calloc(1, sizeof(struct sam_block_t));
    
    sb->fs = fin;
    
    sb->fref = fref;
    
    // initialize the codebook_model
    uint32_t rescale = 1 << 20;
    sb->codebook_model = initialize_stream_model_codebook(rescale);
    
    // Get the Read Length
    if (mode == DECOMPRESSION || mode == DOWNLOAD) {
        // get the readLength from the ios buffer
        sb->read_length =  decompress_int(as, sb->codebook_model);
    }
    else{
        // get the read length from input file and move file pointer after headers
        sb->read_length = get_read_length(sb->fs);
        // write readLength directly to AS using the codebook model
        compress_int(as, sb->codebook_model, sb->read_length);
    }
    
    
    // Allocate the memory:
    
    //READS,
    sb->reads = alloc_read_block_t(sb->read_length);
    
    //IDs
    sb->IDs = alloc_id_block();
    
    //aux
    sb->aux = alloc_aux_block();
    
    //RNAMEs
    sb->rnames = alloc_rname_block();
    
    //MAPQ
    sb->mapq = alloc_mapq_block();
    
    //RNEXT
    sb->rnext = alloc_rnext_block();
    
    //PNEXT
    sb->pnext = alloc_pnext_block();
    
    //TLEN
    sb->tlen = alloc_tlen_block();
    
    return sb;
    
}


uint32_t load_sam_line(sam_block sb){
    
    int32_t j = 0;
    read_line rline = sb->reads->lines;
    qv_line qvline = sb->QVs->qv_lines;
    
/*
    for (i = 0; i < QV->block_length; i++) {
        
        stringstream ss;
        char c;
        while (1) {
            c = fgetc(QV->fs);
            if (c == EOF) return 1;
            if (c == '\n') break;
            ss << c;
        }
        
        
        size_t counter = 0;

        string token;
        while (getline(ss, token, '\t')) {
            counter++;
            // check flag
            if (counter == 2) {
            
            // reached quality values 
            } else if (counter == 11) {
                if (invFlag & 16) { 
                    // The read is inverted, so we should read the quality 
                    // values backwards
                    for (j = QV->columns - 1; j >= 0; j--) {
                        qvline[i].data[j] = (*ptr) - 33, ptr++;
                    }
                }
            } else { 
                // The read is not inversed, so we should read the
                // quality values normally (forward)
                for (j = 0; j < QV->columns; j++) {
                    qvline[i].data[j] = (*ptr) - 33, ptr++;
                }
            }
        }
    }*/
    
    char buffer[1024];
    char *ptr;
    char *ID_line = *sb->IDs->IDs;
    char *rname_line = *sb->rnames->rnames;
    char *rnext = *sb->rnext->rnext;
    
    char **aux_fields = sb->aux->aux_str;

    rline->read_length = sb->read_length;
    // Read compulsory fields
    if (fgets(buffer, 1024, sb->fs)) {
        // ID
        ptr = strtok(buffer, "\t");
        strcpy(ID_line, ptr);
        // FLAG
        ptr = strtok(NULL, "\t");
        rline->invFlag = atoi(ptr);
        // RNAME
        ptr = strtok(NULL, "\t");
        strcpy(rname_line, ptr);
        // POS
        ptr = strtok(NULL, "\t");
        rline->pos = atoi(ptr);
        // MAPQ
        ptr = strtok(NULL, "\t");
        *sb->mapq->mapq = atoi(ptr);
        // CIGAR
        ptr = strtok(NULL, "\t");
        strcpy(rline->cigar, ptr);
        // RNEXT
        ptr = strtok(NULL, "\t");
        strcpy(rnext, ptr);
        // PNEXT
        ptr = strtok(NULL, "\t");
        *sb->pnext->pnext = atoi(ptr);
        // TLEN
        ptr = strtok(NULL, "\t");
        *sb->tlen->tlen = atoi(ptr);
        // SEQ
        ptr = strtok(NULL, "\t");
        strcpy(rline->read, ptr);
        // QUAL
        ptr = strtok(NULL, "\t");
        
        // Read the AUX fields until end of line, and store the MD field
        int auxCnt = 0;
        while( NULL != (ptr = strtok(NULL, "\t")) ){
            //MD:_:_ is a special case
            if (*ptr == 'M' && *(ptr+1) == 'D'){
                // skip MD:Z:
                ptr += 5;
                strcpy(rline->edits, ptr);
                //break;
            } else {
                strcpy(aux_fields[auxCnt], ptr);
                auxCnt++;
                //if we have reached the max. allowed aux fields, break.
                if(auxCnt==MAX_AUX_FIELDS) break;
            }
        }
        
        //awful hack to check if the last field has a \n.
        if(auxCnt!=0) if(aux_fields[auxCnt-1][strlen(aux_fields[auxCnt-1])-1]=='\n') aux_fields[auxCnt-1][strlen(aux_fields[auxCnt-1])-1]=0;
        
        sb->aux->aux_cnt = auxCnt;
        
        return 0;
    }
    else
        return 1;
}


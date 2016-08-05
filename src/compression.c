//
//  compression.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 12/4/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include <stdio.h>
#include <pthread.h>
#include <time.h>
#include "sam_block.h"
#include "read_compression.h"

int print_line(struct sam_line_t *sline, uint8_t print_mode, FILE *fs){
    
    char foo[] = "CIGAR";
    
    int32_t i = 0;
    
    switch (print_mode) {
        case 0:
            if ((sline->flag & 16) == 16) {

                for (i = sline->readLength - 1; i >=0 ; --i)
                    fputc(bp_complement(sline->read[i]), fs);
                
                fputc('\n', fs);
            }
            else{
                fprintf(fs,"%s\n",sline->read);
            }
            break;
            
        default:
            break;
    }
    return 0;
}

int compress_line(Arithmetic_stream as, sam_block samBlock, uint8_t lossiness)  {
    
    uint8_t chr_change;
    
    // Load the data from the file
    if(load_sam_line(samBlock))
        return 0;
    
    if ( (samBlock->reads->lines->invFlag & 4) == 4) {
        return 1;
    }
        
    // Compress sam line
    
    chr_change = compress_rname(as, samBlock->rnames->models, *samBlock->rnames->rnames);
        
    if (chr_change == 1){

        // Store Ref sequence in memory
        store_reference_in_memory(samBlock->fref);
        // Reset cumsumP
        cumsumP = 0;
        memset(snpInRef, 0, MAX_BP_CHR);
    }
    
    compress_read(as, samBlock->reads->models, samBlock->reads->lines, chr_change);
    
    return 1;
}

int decompress_line(Arithmetic_stream as, sam_block samBlock, uint8_t lossiness) {
    
    int32_t chr_change = 0;
    
    uint32_t decompression_flag = 0;
    
    struct sam_line_t sline;
    
    
    //This is only for fixed length? i think so.
    sline.readLength = samBlock->read_length;
    
    //printf("Decompressing the block...\n");
    // Loop over the lines of the sam block
        
    chr_change = decompress_rname(as, samBlock->rnames->models, sline.rname);
        
    if (chr_change == -1)
        return 0;
        
    if (chr_change == 1){
            
        //printf("Chromosome %d decompressed.\n", ++chrCtr);
            
        // Store Ref sequence in memory
        store_reference_in_memory(samBlock->fref);
            
        // reset cumsumP
        cumsumP = 0;

        memset(snpInRef, 0, MAX_BP_CHR);

    }
        
    decompression_flag = decompress_read(as,samBlock, chr_change, &sline);
    
    print_line(&sline, 0, samBlock->fs);
    
    return 1;
}



void* compress(void *thread_info){
    
    uint64_t compress_file_size = 0;
    clock_t begin;
    clock_t ticks;
    
    uint32_t lineCtr = 0;
    
    printf("Compressing...\n");
    begin = clock();
    
    struct compressor_info_t info = *((struct compressor_info_t *)thread_info);
    
    struct qv_options_t opts = *(info.qv_opts);
    
    // Allocs the Arithmetic and the I/O stream
    Arithmetic_stream as = alloc_arithmetic_stream(info.mode, info.fcomp);
    
    // Allocs the different blocks and all the models for the Arithmetic
    sam_block samBlock = alloc_sam_models(as, info.fsam, info.fref, &opts, info.mode);
    
    
    if (info.lossiness == LOSSY) {
        compress_int(as, samBlock->codebook_model, LOSSY);
        initialize_qv_model(as, samBlock->QVs, COMPRESSION);
    }
    else
        compress_int(as, samBlock->codebook_model, LOSSLESS);
    
    while (compress_line(as, samBlock, info.lossiness)) {
        ++lineCtr;
        if (lineCtr % 1000000 == 0) {
          printf("[cbc] compressed %" PRIu32 " lines\n", lineCtr);
        }
    }
    // Load and compress the blocks
    //while(compress_block(as, samBlock)){
    //    reset_QV_block(samBlock->QVs, info.mode);
    
    //   n += samBlock->block_length;
    //}
    
    // Check if we are in the last block
    compress_rname(as, samBlock->rnames->models, "\n");
    
    //end the compression
    compress_file_size = encoder_last_step(as);
    
    printf("Final Size: %lld\n", compress_file_size);
    //printf("%f Million reads compressed using %f MB.\n", (double)n/1000000.0, (double)compress_file_size/1000000.0);
    
    // free(samLine->cigar), free(samLine.edits), free(samLine.read_), free(samLine.identifier), free(samLine.refname);
    
    fclose(info.fsam);
    
    ticks = clock() - begin;
    
    printf("Compression took %f\n", ((float)ticks)/CLOCKS_PER_SEC);
    
    //pthread_exit(NULL);
    return NULL;
}


void* decompress(void *thread_info){
    
    uint64_t n = 0;
    clock_t begin = clock();
    clock_t ticks;
    
    struct compressor_info_t *info = (struct compressor_info_t *)thread_info;
    
    Arithmetic_stream as = alloc_arithmetic_stream(info->mode, info->fcomp);
    
    sam_block samBlock = alloc_sam_models(as, info->fsam, info->fref, info->qv_opts, DECOMPRESSION);
    
    info->lossiness = decompress_int(as, samBlock->codebook_model);
    
    // Start the decompression
    // initialize the QV model
    if (info->lossiness == LOSSY) {
        initialize_qv_model(as, samBlock->QVs, DECOMPRESSION);
    }
    
    // Decompress the blocks
    while(decompress_line(as, samBlock, info->lossiness)){
        //reset_QV_block(samBlock->QVs, DECOMPRESSION);
        n++;
    }
    
    n += samBlock->block_length;
    
    //end the decompression
    //compress_file_size = encoder_last_step(as);
    
    ticks = clock() - begin;
    
    printf("Decompression took %f\n", ((float)ticks)/CLOCKS_PER_SEC);
    
    //printf("%f Million reads decompressed.\n", (double)n/1000000.0);
    
    // free(samLine->cigar), free(samLine.edits), free(samLine.read_), free(samLine.identifier), free(samLine.refname);
    
    fclose(info->fsam);
    fclose(info->fref);
    
    return NULL;
}

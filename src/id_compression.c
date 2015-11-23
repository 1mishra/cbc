//
//  id_compression.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 12/10/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

// Compression of the IDs -- work based on the compression of ID in Samcomp by Mahonney and Bonfiled (2012)

#include <stdio.h>
#include "sam_block.h"

uint8_t decompress_uint8t(Arithmetic_stream as, stream_model model){
    
    // Send the value to the Arithmetic Stream
    uint8_t c = read_value_from_as(as, model);
    
    // Update model
    update_model(model, c);
    
    return c;
    
}


int compress_uint8t(Arithmetic_stream as, stream_model model, uint8_t c){
    
    // Send the value to the Arithmetic Stream
    send_value_to_as(as, model, c);
    
    // Update model
    update_model(model, c);
    
    return 1;
    
}

int compress_rname(Arithmetic_stream as, rname_models models, char *rname){
    
    static char prev_name[1024] = {0};
    static int prevChar = 0;
    
    uint32_t ctr = 0;
    
    if(strcmp(rname, prev_name) == 0){
        
        compress_uint8t(as, models->same_ref[0], 0);
        return 0;
        
    }
    
    else{
        compress_uint8t(as, models->same_ref[0], 1);
        while (*rname) {
            compress_uint8t(as, models->rname[prevChar], *rname);
            prev_name[ctr++] = *rname;
            prevChar = *rname++;
        }
        compress_uint8t(as, models->rname[prevChar], 0);
        prev_name[ctr] = 0;
        return 1;
    }
    
}

int decompress_rname(Arithmetic_stream as, rname_models models, char *rname){
    
    static char prev_name[1024] = {0};
    static int prevChar = 0;
    
    uint8_t chr_change = 0, ch;
    
    uint32_t ctr = 0;
    
    chr_change = decompress_uint8t(as, models->same_ref[0]);
    
    if (chr_change) {
        
        while ( (ch = decompress_uint8t(as, models->rname[prevChar])) ) {
            
            if (ch == '\n') {
                return -1;
            }
            prev_name[ctr++] = ch;
            prevChar = ch;
            *rname = ch, rname++;
        }
        
    }
    
    return chr_change;
    
}

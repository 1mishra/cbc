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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>

typedef struct sam_line_t{
    char *identifier;
    char *refname;
    char *cigar;
    char *edits;
    char *read;
    uint16_t invFlag;
    // int invFlag;
    int pos;
}*sam_line;


uint8_t read_line_from_sam(struct sam_line_t *samLine, FILE *f){

    static uint8_t firstCall = 1;    
    uint8_t thereIsLine = 0;
    char *token;
    char tmp[1000];
    int ch;
   
    if (firstCall == 1){ 
      char line[100000] = {0};

      while ( fgets(line, 100000, f) != NULL) {
          if (line[0] == '@'){
              // Line is a header
              
          }else{
              firstCall = 0;
              // Line belongs to a read
              thereIsLine = 1;
          
              // Identifier
              token = strtok(line, "\t");
              strcpy(samLine->identifier, token);
          
              // Flag
              token = strtok(NULL, "\t");
              samLine->invFlag = atoi(token);
          
              // Refname
              token = strtok(NULL, "\t");
              strcpy(samLine->refname, token);
          
              // Pos
              token = strtok(NULL, "\t");
              samLine->pos = atoi(token);
          
              // MapQ
              token = strtok(NULL, "\t");
          
              // Cigar
              token = strtok(NULL, "\t");
              strcpy(samLine->cigar, token);
          
              // Rnext, Pnext, Tlen
              token = strtok(NULL, "\t");
              token = strtok(NULL, "\t");
              token = strtok(NULL, "\t");
          
              // Read
              token = strtok(NULL, "\t");
              strcpy(samLine->read, token);
          
              // QV sequence
              token = strtok(NULL, "\t");
          
              /* walk through AUX fields */
              samLine->edits[0] = '$';
              token = strtok(NULL, "\t");
              while( token != NULL )
              {
                  if (strncmp(token, "MD:Z:", 5) == 0) {
                      strcpy(samLine->edits, token + 5);
                  }
                  token = strtok(NULL, "\t");
              }
              break;
          }
      }
    }else{
      // firstCall = 0
      // Read compulsory fields
      if (EOF!=fscanf(f, "%s %"SCNu16" %s %d %*d %s %*s %*d %*d %s %*s", samLine->identifier, &samLine->invFlag, samLine->refname, &samLine->pos, samLine->cigar, samLine->read)){
      thereIsLine = 1;
      // Read the AUX fields until end of line, and store the MD field
      while('\n'!=(ch=fgetc(f))){
        // Do something
        if (ch == 'M'){
          ch = fgetc(f);
          if (ch == 'D'){
            // Read :Z:
            ch = fgetc(f);
            ch = fgetc(f);
            ch = fgetc(f);
            fscanf(f, "%s", samLine->edits);
          }
        }

      }

    }


    }
    
    return thereIsLine;

}


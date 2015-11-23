//
//  aux_compression.c
//  samFastqIO
//
//  Created by Carlos Navarro Astiasarán on 29/08/15.
//  Copyright (c) 2015 Carlos Navarro Astiasarán. All rights reserved.
//

#include "sam_block.h"


//Converts a-z to 0-25, A-Z to 26-51, 0-9 to 52-61. (so that 6bits is enough)
uint8_t charMap(char c) {
    if (c >=97 && c<=122) return c-97;
    else if(c>=65 && c<=90) return c-39;
    else if(c>=48 && c<=57) return c+4;
    else return 63;
}

//Inverse of charMap
char inverseCharMap(uint8_t c) {
    if(c>=0 && c<=25) return c+97;
    else if(c>=26 && c<=51) return c+39;
    else if(c>=52 && c<=61) return c-4;
    else return 0;
}

uint8_t typeMap(char t) {
    int index;
    for(index=0;index<TYPELUTLENGTH;index++)
        if(typeLUT[index]==t) return index;
    return TYPELUTLENGTH;
}


/*   ptr: Points to the first char of the tag.
     The uint16_t output has the mapped data.
     rawType just returns the type.
 */
uint16_t preprocessTagType(char *ptr, uint8_t* rawType) {
    char tagChar1,tagChar2,type;
    
    tagChar1 = *ptr;
    tagChar2 = *(ptr+1);
    type     = *(ptr+3);
    *rawType = type;
    
    //We first check whether this TAG:TYPE appears on the LUT
    char tagType[4];
    tagType[0]=tagChar1; tagType[1]=tagChar2; tagType[2]=type; tagType[3]=0;
    int index;
    for(index=0;index<TAGTYPELUTLENGTH;index++)
        if(strcmp(tagTypeLUT[index],tagType)==0) break;
    
    //index=TAGTYPELUTLENGTH; //force to use NOT FOUND (just 2 check)
    if(index==TAGTYPELUTLENGTH) {
        //not found:
        //1st bit flag = 1,
        //2nd to 7th bit = first char of tag (mapped w/charMap),
        //8th to 13th bit = second char of tag (mapped w/charMap),
        //14th to 16th bit = type (mapped w/LUT).
        return NOTFOUNDLUTFLAG | charMap(tagChar1)<<9 | charMap(tagChar2)<<3 | typeMap(type); //what if typeMap(type) is not mapped to anything? (unknown type) then typeMap(type) = TYPELUTLENGTH
        
    } else {
        //found:
        //1st bit flag = 0,
        //2nd to 10th bit = 0,
        //11th to 16th bit = index.
        return index;
    }
}

uint8_t inversePreprocessTagType(char byte1, char byte2, char *ptr)
{
    uint8_t mappedTagChar1,mappedTagChar2,mappedType;
    char tagChar1,tagChar2,type;
    
    uint8_t tagTypeIndex;
    
    uint16_t d = (byte1<<8) | byte2;
    
    //appears on LUT?
    if (d & NOTFOUNDLUTFLAG) {
        //No
        mappedTagChar1 = (d & MASKTAGCHAR1)>>9;
        mappedTagChar2 = (d & MASKTAGCHAR2)>>3;
        mappedType = (d & MASKTYPE);
        tagChar1 = inverseCharMap(mappedTagChar1);
        tagChar2 = inverseCharMap(mappedTagChar2);
        type = typeLUT[mappedType];
        *ptr=tagChar1;
        *(ptr+1)=tagChar2;
        *(ptr+2)=':';
        *(ptr+3)=type;
        *(ptr+4)=':';
        return 1;
    }
    
    //Yes
    tagTypeIndex = d & 0xFF;
    *ptr=(char)tagTypeLUT[tagTypeIndex][0];
    *(ptr+1)=(char)tagTypeLUT[tagTypeIndex][1];
    *(ptr+2)=':';
    *(ptr+3)=(char)tagTypeLUT[tagTypeIndex][2];
    *(ptr+4)=':';
    return 0;
}

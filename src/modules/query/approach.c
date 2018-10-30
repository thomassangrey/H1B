
#ifndef APPROACH_
#define APPROACH_

#include <stdio.h>
#include <stdlib.h>
#include "../utilities/C_arr_L.h"
#include "../utilities/bin_struct.h"
#include "approach.h"
#include <string.h>
#include "../utilities/bin_struct.h"


long* approach(bin_struct BS, double val,  long * labels, 
              char *relational) {
    
    long relational_val;
    long *index;
    index = (long *) malloc(2*sizeof(long)*2);
    index[0] = -1;
    index[1] = -1;


    // Handle various search cases: =, <=, >=, <, >, Nearest

    if (!strncmp(relational, "E", 1) ) {
        relational_val = 1;
    }
    else if (!strncmp(relational, "GE", 2) ) {
        relational_val = 2;
    }
    else if (!strncmp(relational, "LE", 2) ) {
        relational_val = 3;
    }
    else if (!strncmp(relational, "G", 1) ) {
        relational_val = 4;
    }
    else if (!strncmp(relational, "L", 1) ) {
        relational_val = 5;
    }
    else if (!strncmp(relational, "N", 1) ) {
        relational_val = 6;
    }

    else {
        printf("Incorrect relational string\n");
        relational_val = 0;
    }

    switch (relational_val) {

        case 1 :
            
            if ( (labels[1] != -1) && (labels[0] != -1) ) {
                if (BS.labels[labels[0]+1] == val ) {
                    index[0] = BS.positions[labels[0]] + BS.run_lengths[labels[0]] ;
                    index[1] = BS.positions[labels[1]] -1;
    
                    return index;
                }
                else {
                    index[0] = -1;
                    index[1] = -1;
                    return index;
                }
            }
            else if ((labels[0] == -1) & (labels[1] != -1)){
                if (BS.labels[0] == val) {
                    index[0] = BS.positions[0];
                    index[1] = BS.positions[0] + BS.run_lengths[0] -1;
                    return index;
                }
                else {
                    index[0] = -1;
                    index[1] = -1; 
                    return index;
                }
            }
            else if ((labels[0] != -1) & (labels[1] == -1)){
                if (BS.labels[BS.len - 1] == val) {
                    index[0] = BS.positions[labels[0] + 1];
                    index[1] = BS.positions[labels[0] + 1] + \
                                   BS.run_lengths[labels[0] + 1] - 1;
                    return index;
                }
                else {
                    index[0] = -1;
                    index[1] = -1; 
                    return index;
                }
            }
            else {
                index[0] = -1;
                index[1] = -1; 
                return index;
            }
            

            
        case 2 :
            if ( ((labels[1]) != -1) && (labels[0] != -1) ) {
                
                if (BS.labels[labels[0]+1] == val ) {
                    index[0] = BS.positions[labels[0] + 1] ;
                    index[1] = BS.positions[labels[0] + 1] + \
                               BS.run_lengths[labels[0] + 1] -1;
                    return index;
                }
                else {
                    index[0] = BS.positions[labels[1]];
                    index[1] = BS.positions[labels[1]] + BS.run_lengths[labels[1]] - 1;
                    return index;
                }
            }
            else if (labels[0] == -1){
                index[0] = BS.positions[0];
                index[1] = BS.positions[0] + BS.run_lengths[0] - 1;
                return index;
            }
            else if (labels[1] == -1) {
                if (BS.labels[BS.len - 1] >= val) {
                    index[0] = BS.positions[BS.len - 1];
                    index[1] = BS.positions[BS.len - 1] + BS.run_lengths[BS.len - 1] -1;
                    return index;
                }
                else {
                    index[0] = -1;
                    index[1] = -1; 
                    return index;
                }
            }


        case 3 :  

            if ( ((labels[1]) != -1) && (labels[0] != -1) ) {
                
                if (BS.labels[labels[0]+1] == val ) {
                    index[0] = BS.positions[labels[0] + 1] ;
                    index[1] = BS.positions[labels[0] + 1] + \
                               BS.run_lengths[labels[0] + 1] -1;
                    return index;
                }
                else {
                    index[0] = BS.positions[labels[0]];
                    index[1] = BS.positions[labels[0]] + BS.run_lengths[labels[0]] - 1;
                    return index;
                }
            }
            else if (labels[0] == -1){
                index[0] = -1;
                index[1] = -1; 
                return index;                
            }
            else if (labels[1] == -1) {
                    index[0] = BS.positions[BS.len - 1];
                    index[1] = BS.positions[BS.len - 1] + BS.run_lengths[BS.len - 1] -1;
                    return index;
                }


        case 4 :
            if ( ((labels[1]) != -1) && (labels[0] != -1) ) {
            
                index[0] = BS.positions[labels[1]] ;
                index[1] = BS.positions[labels[1]] + BS.run_lengths[labels[1]] - 1;
                return index;
            }
            else if (labels[0] == -1){
                index[0] = BS.positions[labels[1]] ;
                index[1] = BS.positions[labels[1]] + BS.run_lengths[labels[1]] - 1;
                return index;
            }
            else if (labels[1] == -1) {

                index[0] = -1;
                index[1] = -1; 
                return index;
            }

        case 5 :
            if ( ((labels[1]) != -1) && (labels[0] != -1) ) {
            
                index[0] = BS.positions[labels[0]] ;
                index[1] = BS.positions[labels[0]] + BS.run_lengths[labels[0]] - 1;
                printf("index = [%ld,%ld]",index[0], index[1]);

                return index;
            }
            
            else if (labels[0] == -1) {

                index[0] = 0;
                if ( (BS.positions[labels[1]] - 1) > 0 )
                    index[1] = BS.positions[labels[1]] - 1;
                else
                    index[1] = 0;
                printf("index = [%ld,%ld]",index[0], index[1]);
                printf("labels = [%ld,%ld]",labels[0], labels[1]);
                return index; 
            }
            
            else if (labels[1] == -1){
                index[0] = BS.positions[labels[0]] ;
                index[1] = BS.positions[labels[0]] + BS.run_lengths[labels[0]] - 1;
                printf("index = [%ld,%ld]",index[0], index[1]);
                return index;
            }

        case 6 :
            if ( ((labels[1]) != -1) && (labels[0] != -1) ) {
                
                if (BS.labels[labels[0]+1] == val ) {
                    index[0] = BS.positions[labels[0]] + BS.run_lengths[labels[0]];
                    index[1] = BS.positions[labels[1]] -1;
                    return index;
                }
                else if ( (val - BS.labels[labels[0]]) > (BS.labels[labels[1]] - val) ) {
                        index[0] = BS.positions[labels[1]];
                        index[1] = BS.positions[labels[1]] + BS.run_lengths[labels[1]] - 1;
                        return index;
                }
                else if ( (val - BS.labels[labels[0]]) < (BS.labels[labels[1]] - val) ) {
                        index[0] = BS.positions[labels[0]];
                        index[1] = BS.positions[labels[0]] + BS.run_lengths[labels[0]] - 1;
                        return index;
                }
            }
            else if (labels[0] == -1) {
                index[0] = BS.positions[0];
                index[1] = BS.positions[0] + BS.run_lengths[0] -1;
                return index;
            }

            else if (labels[1] == -1) {
                index[0] = BS.positions[BS.len - 1];
                index[1] = BS.positions[BS.len - 1] + BS.run_lengths[BS.len - 1] -1;
                return index;
            }


        default :
            index[1]  = -1; 
            index[0]   = -1;
            return index;

    }

}
#endif




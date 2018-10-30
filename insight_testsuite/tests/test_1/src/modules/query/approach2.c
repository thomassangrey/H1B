
#ifndef APPROACH2_
#define APPROACH2_

#include <stdio.h>
#include <stdlib.h>
#include"../utilities/C_arr_L.h"
#include "../utilities/bin_struct.h" 
#include "approach.h"
#include "approach.c" 



static C_arr_L
LTE_case(long *index, C_arr_L lookup_array, bin_struct BS, C_arr_L temp_array) {
    long idx, Hi, Lo;
    
    Lo = 0;
    if (index[1] == -1) {
        Hi = lookup_array.len -1;
    }
    else {
        Hi = index[1];
    }
    temp_array.C_list_L = (long *) realloc(temp_array.C_list_L, \
                                      2*(Hi - Lo + 1)*sizeof(long));
    for( idx = Lo; idx <= Hi; idx++) {
        temp_array.C_list_L[idx - Lo] = lookup_array.C_list_L[idx];
    }
    temp_array.len = idx - Lo;
    // printf("temp_array.len = %ld\n",temp_array.len);
    return temp_array;
}

static C_arr_L
GTE_case(long *index, C_arr_L lookup_array, bin_struct BS, C_arr_L temp_array) {
    long idx, Hi, Lo;
    
    Hi = lookup_array.len - 1;
    if (index[0] == -1) {
        Lo = 0;
    }
    else {
        Lo = index[0];
    }
    temp_array.C_list_L = (long *) realloc(temp_array.C_list_L, \
                                      2*(Hi - Lo + 1)*sizeof(long));
    for( idx = Lo; idx <= Hi; idx++) {
        temp_array.C_list_L[idx - Lo] = lookup_array.C_list_L[idx];
    }
    temp_array.len = idx - Lo;
    // printf("temp_array.len = %ld\n",temp_array.len);
    return temp_array;
}
static C_arr_L
LT_case(long *index, C_arr_L lookup_array, bin_struct BS, C_arr_L temp_array) {
    long idx, Hi, Lo;
    
    Lo = 0;
    if (index[1] == -1) {
        Hi = lookup_array.len -1;
    }
    else {
        Hi = index[1];
    }
    printf("Hi = %ld, Lo = %ld\n",Hi, Lo);
    temp_array.C_list_L = (long *) realloc(temp_array.C_list_L, \
                                      2*(Hi - Lo + 1)*sizeof(long));
    for( idx = Lo; idx <= Hi; idx++) {
        temp_array.C_list_L[idx - Lo] = lookup_array.C_list_L[idx];
    }
    temp_array.len = idx - Lo;
    return temp_array;
}
static C_arr_L
GT_case(long *index, C_arr_L lookup_array, bin_struct BS, C_arr_L temp_array) {
    
    long idx, Hi, Lo;     
    Hi = lookup_array.len - 1;
    if (index[0] == -1) {
        Lo = 0;
    }
    else {
        Lo = index[0];
    }

    temp_array.C_list_L = (long *) realloc(temp_array.C_list_L, \
                                      2*(Hi - Lo + 1)*sizeof(long));
    for( idx = Lo; idx <= Hi; idx++) {
        temp_array.C_list_L[idx - Lo] = lookup_array.C_list_L[idx];
    }
    temp_array.len = idx - Lo;
    return temp_array;
} 

C_arr_L approach2(long *labels, double val, C_arr_L lookup_array, \
                   bin_struct BS, long relational_val, C_arr_L temp_array) {


char * relational;
long * index;

index = malloc(2*2*sizeof(long));

switch ((int) relational_val) {    
        case 1 :
            relational = "LE";
            index = approach( BS, val, labels, relational); 
            temp_array = LTE_case(index, lookup_array, BS, temp_array);
            break;
        case 2 :
            relational = "GE";
            index = approach( BS, val, labels, relational); 
            temp_array = GTE_case(index, lookup_array, BS, temp_array);
            break;
        case 3 :
            relational = "L";
            index = approach( BS, val, labels, relational); 
            temp_array = LT_case(index, lookup_array, BS, temp_array);
            break;
        case 4 :
            relational = "G";
            index = approach( BS, val, labels, relational); 
            temp_array = GT_case(index, lookup_array, BS, temp_array);
            break;
        
        default:
            temp_array.C_list_L = NULL;
            temp_array.len = -1;
            printf("Approach2: Default switch case. This shouldn't happen...\n");
    }
free(index);
return temp_array;
}
#endif






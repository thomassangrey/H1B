
#include <stdio.h>
#include <stdlib.h>
#include "HiLo.h"
#include "approach.c"
#include <string.h>
#include "time.h"
#include "../sort/sortbymedian.c"


int main() {
    int idx, length, start, middle, end, *index, K, *lookup, median, \
        *temp_lookup;
    double *sorted_series, *random_series, *temp_series, val;
    char * approach;
    
    
    
    srand(time(0));
    index = (int *) malloc(sizeof(int)*2);
    length = 10;
    middle = 0;
    start = -1;
    end = length - 1;
    K=-1;
    median = 0;

    index[0] = 0;
    index[1] = 1;

     val = 8;

    sorted_series = (double *) malloc(sizeof(double)*length);
    random_series = (double *) malloc(sizeof(double)*length);
    temp_series = (double *) malloc(sizeof(double)*length);
    lookup = (int *) malloc(sizeof(int)*length);
    temp_lookup = (int *) malloc(sizeof(int)*length);
    
    length = 10;
        for(idx = 0; idx < length; idx++) {
            random_series[idx] = (int) (length * ((double) rand()/ (double) RAND_MAX) ); 
            sorted_series[idx] = random_series[idx]; 
            temp_series[idx] = random_series[idx]; 
            
            printf("random_series[%d] = %f\n", idx, random_series[idx]);
        }
    
    sort(start, end, K, median, length, sorted_series, temp_series, lookup, temp_lookup);
    
    printf("\n");
    for(idx = 0; idx < length; idx++) {
        printf("sorted_series[%d] = %f\n", idx, sorted_series[idx]);
    }
    approach = "E";
    index = HiLo(sorted_series, val, start, \
                    middle, end, length, index);
    printf("uncorrected index = [%d, %d]\n", index[0],index[1]);
    index = relational(sorted_series, val, length, index, approach); 

    printf("index = [%d, %d]\n", index[0],index[1]);
}

int* HiLo(double * sorted_series, double val, int start, \
            int middle, int end, int length, int * index) {
    
int temp_val;
    if (start == -1) {
        start = 0;
        end = length - 1;
        middle = (int) ((end + start)/2);
        index[0] = 0;
        index[1] = 1;
    }
    
    if (middle == start) {
         printf("pos 0\n");
         index[0] = middle;
         index[1] = middle;
         return index;
     }
    
    else {
        if (val > sorted_series[middle] ){
            start = middle;
            middle = 1+ (int) ((end + start)/2);
            printf("pos A\n");
            HiLo(sorted_series, val, start, middle, end, length, index);
            
        }
        else if (val < sorted_series[middle] ) {
            end = middle;
            middle = (int) ((end + start)/2);
            printf("pos B\n");
            HiLo(sorted_series, val, start, middle, end, length, index);
        }
        else { //2
            if ( sorted_series[start] >= val ) {
                if ( sorted_series[end] <= val ) {
                    index[1] = end;
                    temp_val = sorted_series[end];
                    while (sorted_series[end] == temp_val ) {
                        end -= 1;
                    }
                    index[0] = end+1;
                    printf("pos 1\n");
                }
                else {
                    index[0] = start;
                    temp_val = sorted_series[start];
                    while (sorted_series[start] == temp_val ) {
                        start += 1;
                    }
                    index[1] = start - 1;
                     printf("pos 2\n");
                }
            }

            else {  //1
                if ( sorted_series[end ] <= val ) {
                    index[1] = end;
                    temp_val = sorted_series[end];
                    while ( sorted_series[end - 1]  == temp_val ) {
                        end -= 1;
                    }
                    index[0] = end;
                     printf("pos 3\n");
                }
                else {
                    while ( sorted_series[end - 1]  > val ) {
                        end -= 1;
                    }
                    index[1] = end;
                    while ( sorted_series[start + 1]  < val ) {
                        start += 1;
                    }
                    index[0] = start;
                     printf("pos 4\n");
                }
            }
        }

    }
    return index;

}

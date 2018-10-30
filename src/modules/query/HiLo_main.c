#include "../utilities/bin_struct.h"
#include <string.h>
#include "HiLo2.c"
#include "HiLo2.h"
#include "time.h"
#include "../sort/sortbymedian.c"
#include "limits.h"
#include "../utilities/bin_struct.h"
#include "../utilities/BinSeries.c"
#include "approach.h"
#include "approach.c"

int main() {
    long idx, length, start, middle, end, *index, K, *lookup, median, \
        *temp_lookup, interval, *labels;
    double *sorted_series, *random_series, *temp_series, val;
    char * relational;
    bin_struct BS;
    
    
    
    srand(time(0));
    
    index = (long *) malloc(2*sizeof(long)*2);
    labels = (long *) malloc(2*sizeof(long)*2);
 
    length = 20;
    middle = 0;
    start = -1;
    end = length - 1;
    K=-1;
    median = 0;
    interval = INT_MAX;

    index[0] = -1;
    index[1] = -1;
    index[0] = -1;
    index[1] = -1;

     val = 0;

    sorted_series = (double *) malloc(2*sizeof(double)*length);
    random_series = (double *) malloc(2*sizeof(double)*length);
    temp_series = (double *) malloc(2*sizeof(double)*length);
    lookup = (long *) malloc(2*sizeof(long)*length);
    temp_lookup = (long *) malloc(2*sizeof(long)*length);
    
        for(idx = 0; idx < length; idx++) {
            // random_series[idx] = (long) ((length/10) * ((double) rand()/ (double) RAND_MAX) ); 
            random_series[idx] = (long) (0); 
            sorted_series[idx] = random_series[idx]; 
            temp_series[idx] = random_series[idx]; 
            
            printf("random_series[%ld] = %f\n", idx, random_series[idx]);
        }
    
    sort(start, end, K, median, length, sorted_series, temp_series, lookup, temp_lookup);
    
    printf("\n");
    for(idx = 0; idx < length; idx++) {
        printf("sorted_series[%ld] = %f\n", idx, sorted_series[idx]);
    }
    BS = BinSeries( sorted_series,  length);
    // for(idx = 0; idx < BS.len; idx++) {
    //     printf("BS.labels[%ld] = %f\n", idx, BS.labels[idx]);
    // }
    // for(idx = 0; idx < BS.len; idx++) {
    //     printf("BS.positions[%ld] = %ld\n", idx, BS.positions[idx]);
    // }
    for(idx = 0; idx < BS.len; idx++) {
        printf("BS.run_lengths[%ld] = %ld\n", idx, BS.run_lengths[idx]);
    }
    printf("BS.run_lengths[%d] = %ld\n", 0, BS.run_lengths[0]);
    relational = "E";
    start = -1;
    labels = HiLo(BS.labels, val, start, middle, BS.len - 1, \
    	          BS.len, labels, interval);
    printf("labels = [%ld, %ld]\n", labels[0],labels[1]);
    index = approach(BS, val, labels, relational); 

    printf("raw index = [%ld, %ld]\n", BS.positions[labels[0]],BS.positions[labels[1]]);
    printf("index = [%ld, %ld]\n", index[0],index[1]);
}
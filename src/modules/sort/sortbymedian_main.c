#include <stdio.h>
#include <stdlib.h>
#include "sortbymedian.h"
#include "sortbymedian.c"
#include "time.h"


int main() {
    long start, length, end, K, i;
    double median;
    median      =0.0;
    start       = 0;
    length      = 30000;
    end         = length -1;
    K           = -1;
 
    srand(time(0));

    double *tempseries =  malloc(2*length * sizeof (double));
    double *series =  malloc(2*length * sizeof (double));
    long *lookup =  malloc(2*length * sizeof (long));
    long *templookup =  malloc(2*length * sizeof (long));

    
    if ( (tempseries == NULL) | (series == NULL)) {
      printf("in NULL test\n");
    /* Memory could not be allocated, the program should
       handle the error here as appropriate. */
    } else {
        for (i = start; i < length; i++) {
            tempseries[i] = i;//(long) (length * ((double) rand())/( (double) RAND_MAX )) ;
            lookup[i] = i;
        }
        for (i = start; i < length; i++) {
            series[i] = tempseries[i];
            templookup[i] = i;
        }
    }

    // for (i = 0; i < length; i++) {
    //     printf("series[%ld] = %f\n\n",i,series[i]);
    // }

    start = -1;
    sort(start, end, K, median, length, series, tempseries, lookup, templookup);
    
    
    for (i = 1; i < length; i++) {
         if ((series[i] -series[i-1]) < 0 ) {
          printf("Error\n");
         }
         
        // printf("series[%ld] = %f\tlookup[%ld] = %ld\n",i,series[i], i, lookup[i]);
    }

    
    free(series);  
    free(lookup); 
    free(templookup);  
    free(tempseries); 

    
    printf("goodbye...");
    return 0;
}
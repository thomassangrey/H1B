#include <stdio.h>
#include <stdlib.h>
#include "sortbymedian.c"
#include "time.h"

int sort(int start, int end, int K, double median, int length, \
                    double *series, double *tempseries, \
                    int *lookup, int *templookup);
int main() {
    int start = 0;
    int length = 10;
    int end = length -1;
    int K = -1;
    double median = 0.0;
    // double *series;
    // double *tempseries;
    int i;

    srand(time(0));

    double *tempseries = (double *) malloc(length * sizeof (double));
    double *series = (double *) malloc(length * sizeof (double));
    int *lookup = (int *) malloc(length * sizeof (int));
    int *templookup = (int *) malloc(length * sizeof (int));

    
    if ( (tempseries == NULL) | (series == NULL)) {
    /* Memory could not be allocated, the program should
       handle the error here as appropriate. */
    } else {
        for (i = start; i <= end + 1; i++) {
            tempseries[i] = (int) (length * ((double) rand())/( (double) RAND_MAX ));
            lookup[i] = i;
        }
        for (i = start; i <= end + 1; i++) {
            series[i] = tempseries[i];
            templookup[i] = i;
        }
    }

    for (i = 0; i <= end; i++) {
        printf("series[%d] = %f\n\n",i,series[i]);
    }

    start = -1;
    sort(start, end, K, median, length, series, tempseries, lookup, templookup);
    

    for (i = 0; i <= end; i++) {
        printf("series[%d] = %f\tlookup[%d] = %d\n",i,series[i], i, lookup[i]);
    }
    free(tempseries);  /* We are done with the int objects, and
                   free the associated pointer. */
    tempseries = NULL; /* The pointed-to-data  must not be used again,
                   unless re-assigned by using malloc
                   again. */
    free(series);  /* We are done with the int objects, and
                   free the associated pointer. */
    series = NULL; /* The pointed-to-data  must not be used again,
                   unless re-assigned by using malloc
                   again. */
    free(lookup);  /* We are done with the int objects, and
                   free the associated pointer. */
    tempseries = NULL; /* The pointed-to-data  must not be used again,
                   unless re-assigned by using malloc
                   again. */
    free(templookup);  /* We are done with the int objects, and
                   free the associated pointer. */
    series = NULL; /* The pointed-to-data  must not be used again,
                   unless re-assigned by using malloc
                   again. */

    
    printf("goodbye...");
    
}
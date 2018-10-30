#include <stdio.h>
#include <stdlib.h>
#include "test.h"
int *   test(int sum, int length, int counter, int * tempseries) {
    counter += 1;
    printf("tempseries[%d] = %d (before)\n", counter,tempseries[counter]);
    tempseries[counter] = sum;
    sum += counter;
    printf("tempseries[%d] is = %d (after)\n",counter, tempseries[counter] );
    if ( counter < length ) {
                test(sum, length, counter, tempseries);        
    }

    return tempseries;
    }


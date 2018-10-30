#include <stdio.h>
#include <stdlib.h>
#include "HiLo2.h"


long* HiLo(double * sorted_series, double val, long start, \
            long middle, long end, long length, long * index, long interval) {
    


    if (start == -1) {
        start = 0;
        end = length - 1;
        middle = (long) ((end + start)/2);
        index[0] = -1;
        index[1] = -1;
    }

    if (val > sorted_series[end]) {
        index[0] = end;
        index[1] = -1;
        return index;
    }
    if (val < sorted_series[0]) {
        index[0] = -1;
        index[1] = start;
        return index;
    }
    if ( (end - start + 1) < interval ){
        interval = (end - start + 1);
        if (sorted_series[middle] > val) {
            end = middle;
            middle = (long) ((end + start)/2);
            index[1] = end;
            HiLo(sorted_series, val, start, middle, end, length, \
                               index, interval);
        }
        else if (sorted_series[middle] < val) {
            start = middle;
            middle = (long) ((end + start)/2);
            index[0] = start;
            HiLo(sorted_series, val, start, middle, end, length, \
                                index, interval);
        }
        else {  
            index[0] = middle - 1;
            index[1] = middle + 1;
        }

    }
    return index;
}

#include <stdio.h>
#include <stdlib.h>
#include "bin_struct.h"
bin_struct BinSeries( double * sorted_series,  long length) {

// A simple data compression but I only need it to simplify search methods
// Bin the series creating three arrays: i) unique values labeling the bins,
// ii) number of values in each bin, iii) original positions within sorted_series
// of the first appearance in that series of each unique binned value.
// Retruns a struct containing i), ii), and iii) and their lengths

bin_struct BS;
long idx;
long bins_idx = -1;
long runs_idx;

BS.labels      = (double *) malloc(2*length*sizeof(double));
BS.run_lengths = (long *) malloc(2*length*sizeof(long));
BS.positions   = (long *) malloc(2*length*sizeof(long));

    idx = -1;
	while (idx < length)  {
		idx++;
		bins_idx++;
		BS.labels[bins_idx] = sorted_series[idx];
		BS.positions[bins_idx] = idx;
		runs_idx = 1;
		if(sorted_series[idx+1] == sorted_series[idx]) {
            while (sorted_series[idx+1] == sorted_series[idx]){
            	idx++;
            	runs_idx++;
            }
		}
		BS.run_lengths[bins_idx] = runs_idx;
	}
BS.labels      = (double *) realloc(BS.labels, (bins_idx+1)*sizeof(double));
BS.run_lengths = (long *) realloc(BS.run_lengths, (bins_idx+1)*sizeof(long));
BS.positions   = (long *) realloc(BS.positions, (bins_idx+1)*sizeof(long));
BS.len = bins_idx;
return BS;
}
    

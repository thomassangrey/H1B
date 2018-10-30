#include<stdlib.h>
#include<stdio.h>
#include "../utilities/C_arr_L.h"
#include "intersection.c"
#include "time.h"

int main() {
	long  idx, srn_idx, num_records, B;
	C_arr_L records, temp_array;
 
    srand(time(0));
    temp_array.len =10;
    num_records = 10;
    records.len = num_records;

    records.C_list_L =  malloc(2*num_records*sizeof(long));
    temp_array.C_list_L = malloc(2*temp_array.len*sizeof(long));
    
	for(idx = 0;idx < records.len; idx++) {
		srn_idx = (long) (num_records * ((double) rand()/ (double) RAND_MAX) );   
	    records.C_list_L[idx] = idx;//srn_idx;
	}

	for(idx = 0;idx < temp_array.len; idx++) {
	    srn_idx = (long) (temp_array.len * ((double) rand()/ (double) RAND_MAX) );
	    temp_array.C_list_L[idx] = srn_idx;
	}
	
	for(idx = 0;idx < records.len; idx++) {
	    printf("\noldRecords[%ld] = %ld\n", idx, records.C_list_L[idx]);
	}
	for(idx = 0;idx < temp_array.len; idx++) {
	    printf("\ntemp[%ld] = %ld\n", idx, temp_array.C_list_L[idx]);
	}

	records = intersection(temp_array, records);

	for(idx = 0;idx < records.len; idx++) {
	    printf("\nnewRecords[%ld] = %ld\n", idx, records.C_list_L[idx]);
	}


	for (idx = 0; idx < records.len; idx++) {
        idx++;
        if(idx >0) {
            B = (records.C_list_L[idx] - records.C_list_L[idx-1]);
            if (B<0) {
        	    printf("ERROR");
            }
        }
    }
	free(records.C_list_L);
	free(temp_array.C_list_L);
	
	return 0;
}
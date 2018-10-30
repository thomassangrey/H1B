
#include<stdlib.h>
#include<stdio.h>
#include "../utilities/C_arr_L.h"
#include"../sort/sortbymedian.c"
#include"../sort/sortbymedian.h"
#include"limits.h"


C_arr_L intersection(C_arr_L temp_array, C_arr_L records) {
	long     tmp_idx, idx, jdx,  kdx, NR_maxsize, num_records;
    long     cnt, tmp, start, end, length, K, median;
    long    *lookup, *templookup;
    double  *tempseries, *dbl_temp, *dbl_records;
	C_arr_L  new_records;
    
    num_records = records.len;
    if (temp_array.len < num_records) {
        NR_maxsize = temp_array.len;
    }
    else {
        NR_maxsize = num_records;
    }
    
    dbl_records  = malloc(2*records.len*sizeof(double));
    dbl_temp     = malloc(2*temp_array.len*sizeof(double));

    for (kdx = 0; kdx < records.len; kdx ++) {
        dbl_records[kdx] = records.C_list_L[kdx];
    }
    for (kdx = 0; kdx < temp_array.len; kdx ++) {
        dbl_temp[kdx] = temp_array.C_list_L[kdx];
    }

    
    new_records.C_list_L = (long *) malloc(2*NR_maxsize*sizeof(long));
	
    tempseries = malloc(2*records.len*sizeof(double));
    lookup     = malloc(2*records.len*sizeof(long));
    templookup = malloc(2*records.len*sizeof(long));
    
    start = -1;
    end = records.len - 1;
    K = -1;
    median = -1;
    length = records.len; 
   
    sort(start, end, K, median, records.len, dbl_records, \
                 tempseries, lookup, templookup );
    free(tempseries);
    free(lookup);
    free(templookup);

    tempseries = malloc(2*temp_array.len*sizeof(double));
    lookup     = malloc(2*temp_array.len*sizeof(long));
    templookup = malloc(2*temp_array.len*sizeof(long));

    start = -1;
    end = records.len - 1;
    K = -1;
    median = -1;
    length = records.len; 

    sort(-1, temp_array.len - 1, -1, 0, temp_array.len, dbl_temp, \
                 tempseries, lookup, templookup );
    
    free(tempseries);
    free(lookup);
    free(templookup);

 
    for (kdx = 0; kdx < records.len; kdx ++) {
        records.C_list_L[kdx] = (long) dbl_records[kdx];
    }
    for (kdx = 0; kdx < temp_array.len; kdx ++) {
        temp_array.C_list_L[kdx] = (long) dbl_temp[kdx];
    }



    idx = 0;
    jdx = 0;
    cnt = -1;
    

    while (dbl_records[idx] == temp_array.C_list_L[jdx]) {
        cnt++;
        new_records.C_list_L[cnt] = dbl_records[idx];
        jdx++;
        idx++;
    }
    tmp = dbl_records[idx];
    tmp_idx = idx;
    dbl_records[idx] = -INT_MAX;

    while ((jdx < temp_array.len) & (idx < records.len) ) {
        while ( (temp_array.C_list_L[jdx] < dbl_records[idx]) & \
                              (jdx < temp_array.len) ){
            jdx++;
        }
        if ((jdx < temp_array.len) & (idx < records.len) ) {
            if (dbl_records[idx] == temp_array.C_list_L[jdx]) {
                    cnt++;
                    new_records.C_list_L[cnt] = dbl_records[idx];
                    jdx++;
                    // printf("new_records.C_list_L[%ld] = %ld\n",cnt,new_records.C_list_L[cnt]);
            }
        }
        while ( (idx < records.len) & \
            (dbl_records[idx] < temp_array.C_list_L[jdx]) ) {
                idx++;
        }
        if ((jdx < temp_array.len) & (idx < records.len) ) {
            if (dbl_records[idx] == temp_array.C_list_L[jdx]) {
                    cnt++;
                    new_records.C_list_L[cnt] = dbl_records[idx];
                    idx++;
                    // printf("new_records.C_list_L[%ld] = %ld\n",cnt,new_records.C_list_L[cnt]);
            }
        }
    }            
        
        
    dbl_records[tmp_idx] = tmp;


    records.len = cnt + 1;
    
    records.C_list_L = realloc(records.C_list_L, 2*records.len*sizeof(long));
    for(idx = 0; idx < records.len; idx++) {
        records.C_list_L[idx] = new_records.C_list_L[idx];
    }

    free(new_records.C_list_L);
    free(dbl_records);
    free(dbl_temp);
    
    return records;

}
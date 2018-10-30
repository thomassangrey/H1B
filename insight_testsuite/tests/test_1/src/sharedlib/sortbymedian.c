// #include <stdio.h>
// #include <stdlib.h>
// double * sort(int start, int end, int K, double median, int length, double *series, double *tempseries);

// int main() {
//     int start = 0;
//     int length = 10;
//     int end = length -1;
//     int K = 0;
//     double median = 0.0;
//     // double *series;
//     // double *tempseries;
//     int i;
//     double *result;
//     double res;

//     // for (i = start; i < end + 1; i++) {
//     //     // tempseries[i] = 1.0+56*i -9*i*i;
//     //     res = 1+56*i -9.0*i*i;
//     // }


//     double *tempseries = (double *) malloc(length * sizeof (double));
//     double *series = (double *) malloc(length * sizeof (double));

    
//     if ( (tempseries == NULL) | (series == NULL)) {
//     /* Memory could not be allocated, the program should
//        handle the error here as appropriate. */
//     } else {
//     for (i = start; i < end + 1; i++) {
//         tempseries[i] = 10-i;
//         res = 10-i;
//     for (i = start; i < end + 1; i++) {
//         series[i] = 10-i;
//         res = 10-i;
//     }
//     start =-1;
//     series = sort(start, end, K, median, length, series, tempseries);
    
//     }
//     for (i = start+1; i < end + 1; i++) {
//         printf("series[%d] = %f\n",i,series[i]);
//     }
//     free(tempseries);  /* We are done with the int objects, and
//                    free the associated pointer. */
//     tempseries = NULL; /* The pointed-to-data  must not be used again,
//                    unless re-assigned by using malloc
//                    again. */
//     free(series);  /* We are done with the int objects, and
//                    free the associated pointer. */
//     series = NULL; /* The pointed-to-data  must not be used again,
//                    unless re-assigned by using malloc
//                    again. */
//     }

    
//     printf("hello...");
    
// }

double * sort(int start, int end, int K, double median, int length, \
	                double *series, double *tempseries, \
                    double *lookup, double *templookup) {
    int i;
    int sum;
    int Ldx;
    int Hdx;
    int RB;
    int L_start;
    int L_end;
    int R_start;
    int R_end;
    int bisected;
    int nonsingleton;
    
    if(start== -1) {
        start = 0;
    	end = length -1;
        for (i = start; i < end + 1; i++) {
            tempseries[i] = series[i];
            lookup[i] = i;
            templookup[i] = lookup[i];
        }
             

    K += 1;
     

    }
    sum = 0;
    for (i = start; i < end + 1; i++) {
        sum += tempseries[i];
        // printf("tempseries[i] = %f \n", tempseries[i]);
    
    }
    median = sum/(end+1 -start);

    //printf("the start is %d, end is %d, median is %f \n", start, end, median);
    Ldx = start - 1;
    Hdx = end + 1;
    for (i = start; i < end + 1; i++) {
    	if ( tempseries[i] <= median ){
            Ldx += 1;
            series[Ldx] = tempseries[i];
            lookup[Ldx] = templookup[Ldx];
    	}
    	else if (tempseries[i] > median) {
    		Hdx -= 1;
    		series[Hdx] = tempseries[i];
            lookup[Hdx] = templookup[Hdx];
    	}
    }


    if (Ldx < start) {
        RB = start;
        if (K == 0){
            return lookup;
        }
        for (i=0; i<=length-1; i++) {
            printf("s%i %f: ",i,series[i]);
        }
        printf("\nSeries %d: \n",K);
        return (double *) series;
    }
    else {
    	    RB = Ldx;
    }
        
    /////LEFT FORK
    L_start = start;
    L_end = RB;
    

    for (i = L_start; i <= L_end; i++) {
        tempseries[i] = series[i];
        templookup[i] = lookup[i];
    }
    
    nonsingleton = ( (L_end - L_start) >= 1 );
    bisected =  ( (end - L_end) > 0 );
    if (bisected && nonsingleton) {
    	length = L_end - L_start +1;
    	sum = 0;
    	for (i = start; i <= end + 1; i++) {
            sum += tempseries[i];
        }
        median = sum/length;
        lookup = sort(L_start, L_end, K, median, length, \
                 series, tempseries, lookup. templookup);    
    }
    

    ///// RIGHT FORK
    // end unchanged
    R_start = RB + 1;
    R_end = end;

    for (i = R_start; i <= R_end; i++) {
        tempseries[i] = series[i];
        templookup[i] = lookup[i];
    }

    nonsingleton = ( (R_end - R_start) >= 1 );
    bisected =  ( (end - R_start) >= 0 );
    if (bisected && nonsingleton) {
    	sum = 0;
    	length = R_end - R_start +1;
    	for (i = R_start; i <= R_end + 1; i++) {
            sum += tempseries[i];
        }
        median = sum/length;
        lookup = sort(R_start, R_end, K, median, length, \
                 series, tempseries, lookup, tmeplookup);
    }
    for (i=0; i<=length-1; i++) {
        printf("s%i %f: ",i,series[i]);
    }
    printf("\nseries %d: \n",K);
        return (double *) series;
    return lookup;
}
    




int sort(int start, int end, int K, double median, int length, \
	                double *series, double *tempseries, \
                    int *lookup, int *templookup) {
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

    Ldx = start - 1;
    Hdx = end + 1;
    for (i = start; i < end + 1; i++) {
    	if ( tempseries[i] <= median ){
            Ldx += 1;
            series[Ldx] = tempseries[i];
            lookup[Ldx] = templookup[i];
    	}
    	else if (tempseries[i] > median) {
    		Hdx -= 1;
    		series[Hdx] = tempseries[i];
            lookup[Hdx] = templookup[i];
    	}
    }


    if (Ldx < start) {
        RB = start;
        if (K == 0){
            return 1;
        }
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
        sort(L_start, L_end, K, median, length, \
                 series, tempseries, lookup, templookup);    
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
        sort(R_start, R_end, K, median, length, \
                 series, tempseries, lookup, templookup);
    }


    return 1;
}
    


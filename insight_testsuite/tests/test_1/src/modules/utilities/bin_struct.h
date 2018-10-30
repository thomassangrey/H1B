#ifndef BIN_STRUCT
#define BIN_STRUCT

typedef struct {
        long        len;
        double     *labels;
        long       *run_lengths;
        long       *positions;
    } bin_struct;
#endif
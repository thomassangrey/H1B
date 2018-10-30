
####  Getting started with python ctypes

# sum.c is:

# int our_function(int num_numbers, int *numbers) {
#     int i;
#     int sum = 0;
#     for (i = 0; i < num_numbers; i++) {
#         sum += numbers[i];
#     }
#     return sum;
# }
# cc -fPIC -shared -o libsum.so sum.c
# not sure if the following line is necessary:
# export LD_LIBRARY_PATH=.


import ctypes
import random

S = []
for i in range(10):
    S.append(random.randint(1,10))
lenS = len(S)
_sortbymedian = ctypes.CDLL('sharedlib/libsortbymedian.so')
_sortbymedian.argtype = (ctypes.c_int,ctypes.c_int, ctypes.c_int, \
ctypes.c_double,ctypes.c_int, ctypes.POINTER(ctypes.c_double),\
ctypes.POINTER(ctypes.c_double))
# _sortbymedian._FuncPtr().argtype = (ctypes.c_int, \
#            	ctypes.c_int, ctypes.c_int, ctypes.c_double, \
#            	ctypes.c_int, ctypes.POINTER(ctypes.c_double),\
# 	        ctypes.POINTER(ctypes.c_double))
# _sortbymedian.restype = ctypes.c_double * lenS
_sortbymedian.restype = ctypes.POINTER(ctypes.c_double)
def sort(start, end, K, median, length, series, tempseries):
    array_type = ctypes.c_double * lenS
    series = _sortbymedian.sort(ctypes.c_int(start), \
    ctypes.c_int(end), ctypes.c_int(K), ctypes.c_double(median),\
    ctypes.c_int(length), array_type( * series), array_type( * tempseries))
    return series



series = sort(-1,-1,0,0,10,S,S)

#ifndef PYLISTL_GETLIST
#define PYLISTL_GETLIST
#include "Python.h"
#include "C_arr_L.h"



static C_arr_L 
PyListL_GetList(PyObject *self, PyObject *args) {
    long     n, idx, isLong;
    long     C_list_itm;
    PyObject   *Py_list = NULL, *Py_list_itm = NULL;
    C_arr_L      C_array_L = {-1, NULL};
    
    
    if ( (PyArg_ParseTuple(args, "O", &Py_list)) && PyList_Check(Py_list) ) {
        
        n = PyObject_Length(Py_list);
        C_array_L.len = n;
    
    C_array_L.C_list_L = (long *) malloc(2* n*sizeof(long) );

        for(idx = 0; idx < n; idx++) {
            
            Py_list_itm = PyList_GetItem(Py_list, idx);
            Py_XINCREF(Py_list_itm);
            isLong   = PyLong_Check(Py_list_itm);
            
            if ( isLong > 0 ) {
                C_list_itm = PyLong_AsLong(Py_list_itm);
                Py_DECREF(Py_list_itm);
                C_array_L.C_list_L[idx] = C_list_itm;
            }
            else {
                Py_DECREF(Py_list_itm);
                Py_DECREF(Py_list);
                free(C_array_L.C_list_L); 
                printf("Error in PyListL_GetList type check...");
                return  C_array_L;                    
                }
            }
        } 

    else {
        printf("Error in PyListL_GetList... at List_check or ParseTuple\n");
        
        return C_array_L;
    }
    // printf("Exiting PyListD_GetList...\n");
    return C_array_L;
}

#endif


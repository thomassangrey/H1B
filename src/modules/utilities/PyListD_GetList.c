
#ifndef PYLISTD_GETLIST
#define PYLISTD_GETLIST

#include "Python.h"
#include "C_arr_D.h"



static C_arr_D 
PyListD_GetList(PyObject *self, PyObject *args) {
    long       n, idx, isFloat, isLong, conversion = 0;
    double     C_list_itm;
    PyObject   *Py_list = NULL, *Py_list_itm = NULL;
    C_arr_D     C_array_D = {-1, NULL};
    
    if ( (PyArg_ParseTuple(args, "O", &Py_list)) && PyList_Check(Py_list) ) {

        n = PyObject_Length(Py_list);
        C_array_D.len = n;
    
        C_array_D.C_list_D = (double *) malloc(n*sizeof(double) );

        for(idx = 0; idx < n; idx++) {
            
            Py_list_itm = PyList_GetItem(Py_list, idx);
            // printf("1: PyList_itm refcnt = %lu",Py_list_itm->ob_refcnt);
            Py_XINCREF(Py_list_itm);
            isFloat = PyFloat_Check(Py_list_itm);
            
            
            if ( isFloat > 0 ) {
                C_list_itm = PyFloat_AsDouble(Py_list_itm);
                Py_DECREF(Py_list_itm);
                C_array_D.C_list_D[idx] = C_list_itm;
                // printf("2: PyList_itm refcnt = %lu",Py_list_itm->ob_refcnt);
            }
            else {
                PyErr_Clear();
                isLong = PyLong_Check(Py_list_itm);
                if (isLong) {
                    C_list_itm = PyLong_AsDouble(Py_list_itm); 
                    Py_DECREF(Py_list_itm);
                    C_array_D.C_list_D[idx] = C_list_itm;
                    conversion = 1;
                    // printf("3: PyList_itm refcnt = %lu",Py_list_itm->ob_refcnt);
                }
                else {
                    Py_DECREF(Py_list_itm);
                    Py_DECREF(Py_list);
                    free(C_array_D.C_list_D); 
                    printf("Error in PyListD_GetList... Unknown type\n");
                    return  C_array_D;         
                }        
            }
        }
    }

    else {
        printf("Error in PyListD_GetList... at List_check or ParseTuple\n");
        return C_array_D;
    }
    if(conversion) {
        // printf("Forced conversion of PyLong to Double...\n");   
    }
    // printf("Exiting PyListD_GetList...\n");
    return C_array_D;
}

#endif

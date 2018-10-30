#include <stdio.h>
#include <stdlib.h>
#include "Python.h"
#include "C_arr.h"
#include "PyListD_GetList.c"
#include "sortbymedian.h"
#include "sortbymedian.c"
#include "time.h"
 
static char module_docstring[] = "Sort Map heap sorts and returns a lookup";

static char sort_docstring[] = \
    "Sorts a Python list of doubles and returns a sort order lookup map. O(n*ln(n)) ";

static char randomize_docstring[] = \
    "randomizes a list in place. Returns List for success, NULL for fail.";

static PyObject * 
sort_map_sort(PyObject *self, PyObject *args) {
    int           n, idx, start, end, K,length, *lookup, *templookup;
    double       *tempseries, median;
    PyObject     *Py_list_out, *Py_list;
    C_arr         C_array;
     
    C_array = PyListD_GetList(self, args);
    n = C_array.len;   
    // printf("NEW SORTMAP\n");
     
    tempseries = (double *)malloc( n*sizeof(double) );
    templookup = (int *)malloc( n*sizeof(int) );
    lookup =     (int *)malloc( n*sizeof(int) );
    start = -1;
    end = n - 1; 
    length = n; 
    median = 0.0;
    K= -1;

    for (idx = 0; idx < n; idx++) {
        tempseries[idx] = C_array.C_list[idx];
        lookup[idx] = idx;
        templookup[idx] = lookup[idx];
    }
    
    sort(start, end, K, median, length, C_array.C_list, tempseries, lookup, templookup);
    
    // printf("n = %i elements\n",n);
    Py_list_out = PyList_New(n);
    
    for (idx = 0; idx < n; idx++) {
        if ( PyList_SetItem( Py_list_out, idx, \
                             PyLong_FromLong(lookup[idx]) ) < 0 ) {
            printf("Error in PyList_out_SetItem");
            Py_DECREF(Py_list_out);
            free(C_array.C_list);
            free(tempseries);
            free(lookup);
            free(templookup);
            return NULL;
        }
    }
    if ( PyArg_ParseTuple(args, "O", &Py_list) ) {
        
        for (idx = 0; idx < n; idx++) {
            if ( PyList_SetItem( Py_list, idx, \
                                 PyFloat_FromDouble(C_array.C_list[idx]) ) < 0 ) {
                printf("Error in PyList_SetItem");
                Py_DECREF(Py_list);
                free(C_array.C_list);
                free(tempseries);
                free(lookup);
                free(templookup);
                return NULL;
            }
        }
    }
    
    // printf("Exiting sort...\n"); 
    free(C_array.C_list);
    free(tempseries);  
    free(lookup);
    free(templookup);
    return Py_list_out;
}

static PyObject *
sort_map_randomize(PyObject *self, PyObject *args) {
    int           n, idx,  srn_idx;
    PyObject     *Py_list, *result;
    C_arr         C_array;
    srand(time(0));
    C_array = PyListD_GetList(self, args);
    n = C_array.len;

    if ( PyArg_ParseTuple(args, "O", &Py_list) ) {
        Py_XINCREF(Py_list);
        for (idx = 0; idx < n; idx++) {
            srn_idx = (int) (n * ((double) rand()/ (double) RAND_MAX) );
            if ( PyList_SetItem( Py_list, idx, \
                                 PyFloat_FromDouble(C_array.C_list[srn_idx]) ) < 0 ) {
                printf("Error in PyList_SetItem");
                Py_DECREF(Py_list); 
                free(C_array.C_list);
                return NULL;
            }
        }
    }
    else {
        printf("Error in PyList_GetList");
        free(C_array.C_list);
        Py_DECREF(Py_list); 
        return NULL;
    }
    // printf("sort_map.randomize exiting...\n");
    free(C_array.C_list);
    // return Py_list; 
    Py_DECREF(Py_list);
    Py_INCREF(Py_None);
    result = Py_None;
    return result;
}

static PyMethodDef module_methods[] = {
    {"sort", (PyCFunction)sort_map_sort, METH_VARARGS, sort_docstring},
    {"randomize", (PyCFunction)sort_map_randomize, METH_VARARGS, randomize_docstring},
    {NULL, NULL, -1, NULL} 
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "sort_map",
        module_docstring,
        -1,
        module_methods,
        NULL,
        NULL,
        NULL,
        NULL,
};

PyMODINIT_FUNC
PyInit_sort_map(void)
{
    PyObject *module = PyModule_Create(&moduledef);
    return module;

}


// /* Boilerplate to return "None" */
//         Py_INCREF(Py_None);
//         result = Py_None;

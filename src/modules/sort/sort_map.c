#include <stdio.h>
#include <stdlib.h>
#include "Python.h"
#include "../utilities/C_arr_D.h"
#include "../utilities/C_arr_L.h"
#include "../utilities/PyListD_GetList.c"
#include "sortbymedian.h"
#include "sortbymedian.c"
#include "time.h"
#include "../utilities/bin_struct.h"
#include "../utilities/BinSeries.c"  

static char module_docstring[] = \
"Sort Map module creates a heap sort of a list of longs or doubles \
and returns a lookup table.";

static char sort_docstring[] = \
    "Sorts a Python list of doubles or longs and returns a sort \
    order lookup map. O(n*ln(n)) ";
 
static char randomize_docstring[] = \
    "randomizes a list in place. Returns List for success, NULL for fail.";

static char build_bin_docstring[] = \
    "Bins a list. returns 3 lists: 1) unique element values in ascending \
     order, 2) the number of occurences each unique element (run-length), \
     and 3) the starting position of each run";


static PyObject *  
sort_map_build_bin(PyObject *self, PyObject *args) {
// Bins the series to build a PyObject called Py_dict_out. Has four list items:
// 1) len: length: the number of unique elements, 2) labels: unique values,
// 3) the starting positions for each unique value (built from an ordered list),
// and 4) the number of unique elements.
    long          idx;
    PyObject     *Py_dict_out, *Py_list_labels, *Py_list_positions, \
                 *Py_list_run_lengths;
    bin_struct    BS;
    C_arr_D       SS_array; 

    SS_array = PyListD_GetList(self, args);

    BS = BinSeries(SS_array.C_list_D, SS_array.len);  
    Py_list_labels = PyList_New(BS.len); 
    Py_list_positions = PyList_New(BS.len);
    Py_list_run_lengths = PyList_New(BS.len);
    
    for(idx = 0; idx < BS.len; idx++) {
        PyList_SetItem(Py_list_labels, idx, PyFloat_FromDouble(BS.labels[idx]));
        PyList_SetItem(Py_list_positions, idx, PyLong_FromLong(BS.positions[idx]));
        PyList_SetItem(Py_list_run_lengths, idx, PyLong_FromLong(BS.run_lengths[idx]));
    }

    Py_dict_out = PyDict_New(); 
    PyDict_SetItemString(Py_dict_out, "len", PyLong_FromLong(BS.len));
    PyDict_SetItemString(Py_dict_out, "labels", Py_list_labels);
    PyDict_SetItemString(Py_dict_out, "positions", Py_list_positions);
    PyDict_SetItemString(Py_dict_out, "run_lengths", Py_list_run_lengths);
    free(SS_array.C_list_D);
    return Py_dict_out; 
}

static PyObject * 
sort_map_sort(PyObject *self, PyObject *args) {
    long           n, idx, start, end, K,length, *lookup, *templookup;
    double       *tempseries, median;
    PyObject     *Py_list_out, *Py_list;
    C_arr_D       C_array_D;
     
    C_array_D = PyListD_GetList(self, args); 
    n = C_array_D.len;    
    
     
    tempseries = (double *)malloc(2* n*sizeof(double) );
    templookup = (long *)malloc(2* n*sizeof(long) );
    lookup =     (long *)malloc(2* n*sizeof(long) );
    start = -1;
    end = n - 1; 
    length = n; 
    median = 0.0;
    K= -1;
 
    for (idx = 0; idx < n; idx++) {
        tempseries[idx] = C_array_D.C_list_D[idx];
        lookup[idx] = idx; 
        templookup[idx] = lookup[idx]; 
    }
     
    sort(start, end, K, median, length, C_array_D.C_list_D, tempseries, lookup, templookup);
    
    // printf("n = %ld elements\n",n);
    Py_list_out = PyList_New(n);
    
    for (idx = 0; idx < n; idx++) { 
        if ( PyList_SetItem( Py_list_out, idx, \
                             PyLong_FromLong(lookup[idx]) ) < 0 ) {
            printf("Error in PyList_out_SetItem");
            Py_DECREF(Py_list_out);
            free(C_array_D.C_list_D);
            free(tempseries);
            free(lookup);
            free(templookup); 
            return NULL;
        }
    }
    if ( PyArg_ParseTuple(args, "O", &Py_list) ) {
         
        for (idx = 0; idx < n; idx++) {
            if ( PyList_SetItem( Py_list, idx, \
                                 PyFloat_FromDouble(C_array_D.C_list_D[idx]) ) < 0 ) {
                printf("Error in PyList_SetItem");
                Py_DECREF(Py_list);
                free(C_array_D.C_list_D);
                free(tempseries);
                free(lookup);
                free(templookup);
                return NULL; 
            } 
        }  
    } 
    
    // printf("Exiting sort...\n"); 
    free(C_array_D.C_list_D);
    free(tempseries);  
    free(lookup);
    free(templookup);
    return Py_list_out;
}

static PyObject *
sort_map_randomize(PyObject *self, PyObject *args) {
    long           n, idx,  srn_idx;
    PyObject     *Py_list, *result;
    C_arr_D       C_array_D;
    srand(time(0));

    C_array_D = PyListD_GetList(self, args);
    n = C_array_D.len;

    if ( PyArg_ParseTuple(args, "O", &Py_list) ) {
        Py_XINCREF(Py_list);
        for (idx = 0; idx < n; idx++) {
            srn_idx = (long) (n * ((double) rand()/ (double) RAND_MAX) );
            if ( PyList_SetItem( Py_list, idx, \
                                 PyFloat_FromDouble(C_array_D.C_list_D[srn_idx]) ) < 0 ) {
                printf("Error in PyList_SetItem");
                Py_DECREF(Py_list); 
                free(C_array_D.C_list_D);
                return NULL;
            }
        } 
    }
    else {
        printf("Error in PyList_GetList");
        free(C_array_D.C_list_D);
        Py_DECREF(Py_list); 
        return NULL;
    }
    // printf("sort_map.randomize exiting...\n");
    free(C_array_D.C_list_D);
    // return Py_list; 
    Py_DECREF(Py_list);
    Py_INCREF(Py_None);
    result = Py_None;
    return result;
}

static PyMethodDef module_methods[] = {
    {"sort", (PyCFunction)sort_map_sort, METH_VARARGS, sort_docstring},
    {"randomize", (PyCFunction)sort_map_randomize, METH_VARARGS, randomize_docstring},
    {"build_bin", (PyCFunction)sort_map_build_bin, METH_VARARGS, build_bin_docstring},
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

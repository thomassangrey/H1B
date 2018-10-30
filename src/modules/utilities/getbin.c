
#ifndef GETBIN
#define GETBIN

#include "Python.h"
#include "bin_struct.h"
#include "PyListD_GetList.c"
#include "PyListL_GetList.c"

static bin_struct 
getbin(PyObject *self, PyObject *args) {
    PyObject     *Py_dict, *labels_list, *run_lengths_list, *positions_list;
    bin_struct    BS;
    if ( (PyArg_ParseTuple(args, "O", &Py_dict)) && PyDict_Check(Py_dict) ){
        BS.len = (long) PyLong_AsLong(PyDict_GetItemString(Py_dict, "len"));
        BS.labels      = (double *) malloc(2*BS.len*sizeof(double));
        BS.positions   = (long *) malloc(2*BS.len*sizeof(long));
        BS.run_lengths = (long *) malloc(2*BS.len*sizeof(long));

        labels_list      = Py_BuildValue("(O)", PyDict_GetItemString(Py_dict, "labels"));
        positions_list   = Py_BuildValue("(O)", PyDict_GetItemString(Py_dict, "positions"));
        run_lengths_list = Py_BuildValue("(O)", PyDict_GetItemString(Py_dict, "run_lengths"));

        BS.labels =      (PyListD_GetList(self, labels_list)).C_list_D;
        BS.positions =   ( (PyListL_GetList(self, positions_list)).C_list_L );  
        BS.run_lengths = ( (PyListL_GetList(self, run_lengths_list)).C_list_L );//
        
        return BS;
    }
    else {
        BS.len = -1;
        BS.labels      = NULL;
        BS.run_lengths = NULL;
        BS.positions   = NULL;
        printf("Error in getbin.c PyArg_ParseTuple\n");
        return BS;
    }
}

#endif

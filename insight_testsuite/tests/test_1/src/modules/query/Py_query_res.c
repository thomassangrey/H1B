#include <stdio.h>
#include <stdlib.h>
#include "Python.h"
#include "../utilities/C_arr_D.h"
#include "../utilities/PyListD_GetList.c"
#include "../utilities/C_arr_L.h"
#include "../utilities/PyListL_GetList.c"
#include "HiLo2.h"   
#include "HiLo2.c"
#include "intersection.h"
#include "intersection.c"
#include <string.h> 
#include "limits.h"
#include "../utilities/bin_struct.h"  
#include "../utilities/getbin.c"
#include "approach2.h"
#include "approach2.c" 
#include "approach.c" 
#include "approach.h"


static char module_docstring[] = "creates a subspace matching keyed values. ";



static PyObject * 
query_find(PyObject *self, PyObject *args) {
    long          idx, qdx, pdx, tdx, vdx, temp_idx, new, \
                  start, middle, end, length, interval, \
                  n_S, n_V, n_SS, n_L, inp_grp_num;
    long         *index, *labels, relational_val;
    PyObject     *Py_list_out, *Py_list_V, *Py_list_S, *Py_list_SS, \
                 *Py_list_R, *Py_list_L, *Py_dict_BS, *Py_relational;
    C_arr_D       series_array, sorted_series_array, val_array;  
    C_arr_L       lookup_array, records, temp_array;
    Py_ssize_t    tuple_size;
    bin_struct    BS; 
    char *        relational;

    tuple_size = PyTuple_Size(args); 
    inp_grp_num = (long) tuple_size/6;
    index = (long *) malloc( 2*sizeof(long) ); 
    labels = (long *) malloc( 2*sizeof(long) ); 

    Py_relational = Py_BuildValue("(O)", PyTuple_GetItem(args, 0));
    PyArg_ParseTuple(Py_relational, "z", &relational);
 
    if (!strncmp(relational, "LTE", 3) ) {
        relational_val = 1;
    }
    else if (!strncmp(relational, "GTE", 3) ) {
        relational_val = 2;
    }
    else if (!strncmp(relational, "LT", 2) ) {
        relational_val = 3;
    }
    else if (!strncmp(relational, "GT", 2) ) {
        relational_val = 4;
    }
    else {
        relational_val = 0;
    }  
    
    if (tuple_size == 7) {
        Py_list_R  = Py_BuildValue("(O)", PyTuple_GetItem(args, 6));
    }  

 // qdx indexes input group number for groups of searches. Not currently
//  implemented elsewhere in Python program

    for (qdx = 0; qdx < inp_grp_num; qdx++) {
        Py_list_V = Py_BuildValue("(O)", PyTuple_GetItem(args, qdx*6 +1));
        Py_list_S = Py_BuildValue("(O)", PyTuple_GetItem(args, qdx*6 + 2));
        Py_list_SS = Py_BuildValue("(O)", PyTuple_GetItem(args, qdx*6 + 3));
        Py_list_L  = Py_BuildValue("(O)", PyTuple_GetItem(args, qdx*6 + 4));
        Py_dict_BS = Py_BuildValue("(O)", PyTuple_GetItem(args, qdx*6 + 5));
      
        val_array = PyListD_GetList(self, Py_list_V);
        series_array = PyListD_GetList(self, Py_list_S);
        sorted_series_array = PyListD_GetList(self, Py_list_SS);
        lookup_array = PyListL_GetList(self, Py_list_L);
        BS = getbin(self, Py_dict_BS);

        n_V = val_array.len; 
        n_S = series_array.len;
        n_SS = sorted_series_array.len;
        n_L = lookup_array.len;
        
        temp_array.C_list_L= (long *) malloc(2*n_S*sizeof(long));
        start = -1;
        end  = n_S - 1;
        length = n_S; 
        middle =0 ;
        interval = INT_MAX;
        temp_array.len = 0; 
        temp_idx = -1;  

// ?? Put filter in place for handling len(val_array) > 1 AND 
// relatioanl_val == 0 ??
        for (vdx = 0; vdx < n_V; vdx++){
            index[0]  = -1; 
            index[1]  = -1;
            labels[0] = -1; 
            labels[1] = -1;
            
            labels = HiLo(BS.labels, val_array.C_list_D[vdx], start, \
                              middle, end, BS.len, labels, interval);
        
            if (relational_val == 0) {
                index = approach( BS, val_array.C_list_D[vdx], \
                                  labels, relational); 
                
                if ( (index[0] > -1) && (index[1] > -1) ) {
                    
                    for(pdx = index[0]; pdx <=index[1]; pdx++){
                        new = 1;
                        for(tdx = 0; tdx < temp_array.len; tdx++){ //CHANGED < to <=
     
                            if ( temp_array.C_list_L[tdx] == (long) lookup_array.C_list_L[pdx]) {
                                new = 0;
                                break;  
                            } 
                        }
                        if (new) {
                            
                            temp_idx = temp_idx + 1;
                            temp_array.C_list_L[temp_idx] = (long) lookup_array.C_list_L[pdx]; 
                            } 
                    }
                    temp_array.len = temp_idx + 1;
                    
                }
                else {
                     printf("index[0] or index[1] is -1 at vdx = %ld\n", vdx); 
                     printf("index = [%ld, %ld]", index[0], index[1]);
                     printf("qdx = %ld\n",qdx);
                }
            }
            else {
                temp_array = approach2(labels, val_array.C_list_D[vdx], \
                                         lookup_array, BS, \
                                         relational_val, temp_array);

            } 
        }
 
      // Intersect temp_array.C_list_Lwith records
        if (qdx == 0) {
            if (tuple_size == 7) {
                records= PyListL_GetList(self, Py_list_R);
                records = intersection(temp_array, records);
                }
            else {
                records.C_list_L = (long *) malloc(2*temp_array.len*sizeof(long));
                
                for(idx = 0; idx < temp_array.len; idx++) {
                    records.C_list_L[idx] =  temp_array.C_list_L[idx];
               }
                records.len = temp_array.len;  
            }
            
        }
        else {
            
            records = intersection(temp_array, records);
        }

        free(temp_array.C_list_L); 
        free(val_array.C_list_D);
        free(series_array.C_list_D);
        free(sorted_series_array.C_list_D);
        free(lookup_array.C_list_L);
        free(BS.labels);
        free(BS.positions);
        free(BS.run_lengths);
        free(labels);
        

    
        Py_DECREF(Py_list_V);
        Py_DECREF(Py_list_S);
        Py_DECREF(Py_list_SS);
        Py_DECREF(Py_list_L);
        Py_DECREF(Py_dict_BS);
        Py_DECREF(Py_relational);
        if (tuple_size == 7) {
            Py_DECREF(Py_list_R);
        }    
    }

    // compare record_tmp to current records and form the intersection


    Py_list_out = PyList_New(records.len);
                    
    if (records.len > 0 ) {
        for (idx = 0; idx < records.len; idx++) {

            if ( PyList_SetItem( Py_list_out, idx, \
                                 PyLong_FromLong(records.C_list_L[idx])) < 0 ) {
                printf("Error in PyList_out_SetItem");
                Py_DECREF(Py_list_out);
                free(records.C_list_L);
                free(index);
                printf("Exiting query...\n");   
                return NULL;
            }
        } 
    }
    else {
        Py_INCREF(Py_None);
        Py_list_out = Py_None;
    }
    free(records.C_list_L);
    free(index);
    printf("Exiting query...\n");
    // printf("Py_list_S refcnt = %lu\n",Py_list_S->ob_refcnt);
    // printf("Py_list_SS refcnt = %lu\n",Py_list_SS->ob_refcnt);
    // printf("Py_list_L refcnt = %lu\n",Py_list_L->ob_refcnt);
    // printf("Py_dict_BS refcnt = %lu\n",Py_dict_BS->ob_refcnt);
    // printf("Py_relational refcnt = %lu\n",Py_relational->ob_refcnt);
    // printf("Py_list_out refcnt = %lu\n",Py_list_out->ob_refcnt);
    
    return Py_list_out;
}




static PyMethodDef module_methods[] = {
    {"find", (PyCFunction)query_find, METH_VARARGS, ""},
    // {"getbin", (PyCFunction)query_getbin, METH_VARARGS, "getbin"},
    {NULL, NULL, -1, NULL} 
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "query",
        module_docstring,
        -1,
        module_methods,
        NULL,
        NULL,
        NULL,
        NULL,
};

PyMODINIT_FUNC
PyInit_query(void)
{
    PyObject *module = PyModule_Create(&moduledef);
    return module;

}



#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL ASE_ARRAY_API
#include <numpy/arrayobject.h>

#define PY3 (PY_MAJOR_VERSION >= 3)

// Holonomic constraints
PyObject* adjust_positions(PyObject *self, PyObject *args);
PyObject* adjust_momenta(PyObject *self, PyObject *args);
PyObject* adjust_positions_general(PyObject *self, PyObject *args);
PyObject* adjust_momenta_general(PyObject *self, PyObject *args);
// TIP3P forces
PyObject* calculate_forces_H2O(PyObject *self, PyObject *args);


static PyMethodDef functions[] = {
    {"adjust_positions", adjust_positions, METH_VARARGS, 0},
    {"adjust_momenta", adjust_momenta, METH_VARARGS, 0},
    {"adjust_positions_general", adjust_positions_general, METH_VARARGS, 0},
    {"adjust_momenta_general", adjust_momenta_general, METH_VARARGS, 0},
    {"calculate_forces_H2O", calculate_forces_H2O, METH_VARARGS, 0},
    {0, 0, 0, 0}
};

#if PY3
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_ase",
    "C-extension for ASE",
    -1,
    functions,
    NULL,
    NULL,
    NULL,
    NULL
};
#endif

static PyObject* moduleinit(void)
{
#if PY3
    PyObject* m = PyModule_Create(&moduledef);
#else
    PyObject* m = Py_InitModule3("_ase", functions,
                                 "C-extension for ASE\n\n...\n");
#endif

    if (m == NULL)
        return NULL;

#ifndef PARALLEL
    // ase-python needs to import arrays at the right time, so this is
    // done in gpaw_main().  In serial, we just do it here:
    import_array1(0);
#endif
    return m;
}

#ifndef ASE_INTERPRETER

#if PY3
PyMODINIT_FUNC PyInit__ase(void)
{
    return moduleinit();
}
#else
PyMODINIT_FUNC init_ase(void)
{
    moduleinit();
}
#endif

#else // ifndef GPAW_INTERPRETER

#if PY3
#define moduleinit0 moduleinit
#else
void moduleinit0(void) { moduleinit(); }
#endif


int
ase_main()
{
    int status = -1;

    PyObject *ase_mod = NULL, *pymain = NULL;

    ase_mod = PyImport_ImportModule("ase");
    if(ase_mod == NULL) {
        status = 3;  // Basic import failure
    } else {
        pymain = PyObject_GetAttrString(ase_mod, "main");
    }

    if(pymain == NULL) {
        status = 4;  // gpaw.main does not exist for some reason
        //PyErr_Print();
    } else {
        // Returns Py_None or NULL (error after calling user script)
        // We already imported the Python parts of numpy.  If we want, we can
        // later attempt to broadcast the numpy C API imports, too.
        // However I don't know how many files they are, and we need to
        // figure out how to broadcast extension modules (shared objects).
        import_array1(0);
        PyObject *pyreturn = PyObject_CallFunction(pymain, "");
        status = (pyreturn == NULL);
        Py_XDECREF(pyreturn);
    }

    Py_XDECREF(pymain);
    Py_XDECREF(ase_mod);
    return status;
}


int
main(int argc, char **argv)
{
#if PY3
#define PyChar wchar_t
    wchar_t* wargv[argc];
    wchar_t* wargv2[argc];
    for (int i = 0; i < argc; i++) {
        int n = 1 + mbstowcs(NULL, argv[i], 0);
        wargv[i] = (wchar_t*)malloc(n * sizeof(wchar_t));
        wargv2[i] = wargv[i];
        mbstowcs(wargv[i], argv[i], n);
    }
#else
#define PyChar char
    char** wargv = argv;
#endif

    Py_SetProgramName(wargv[0]);
    PyImport_AppendInittab("_ase", &moduleinit0);
    Py_Initialize();

    PySys_SetArgvEx(argc, wargv, 0);
    int status = ase_main();

    if(status != 0) {
        PyErr_Print();
    }

    Py_Finalize();

#if PY3
    for (int i = 0; i < argc; i++)
        free(wargv2[i]);
#endif

    return status;
}
#endif // ASE_INTERPRETER

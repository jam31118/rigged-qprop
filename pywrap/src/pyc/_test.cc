#include <Python.h>
#include "imag-prop.h"

// Declare each module method
static PyObject *_test_imag_prop(PyObject *self, PyObject *args);

// Define an array of module methods
static PyMethodDef module_methods[] = {
	{"imag_prop", _test_imag_prop, METH_VARARGS, "An imaginary propagation"},
	{NULL, NULL, 0, NULL}
};

// Define this module
static struct PyModuleDef _test = {
	PyModuleDef_HEAD_INIT,
	"_test",
	"DOCUMENTATION: A set of arithmetics",
	-1,
	module_methods
};

// Create module object
PyMODINIT_FUNC PyInit__test(void) {
	return PyModule_Create(&_test);
}

// Define module methods
static PyObject *_test_imag_prop(PyObject *self, PyObject *args) {
	if (imag_prop() != 0) {
		PyErr_SetString(PyExc_Exception, "Failed to execute imaginary propagation.");
		return NULL;
	}

	Py_INCREF(Py_None);
	return Py_None;
}


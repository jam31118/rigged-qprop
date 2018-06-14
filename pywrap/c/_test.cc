#include <Python.h>
#include "imag-prop.h"
#include "real-prop.h"


// Declare each module method
static PyObject *_test_imag_prop(PyObject *self, PyObject *args);
static PyObject *_test_real_prop(PyObject *self, PyObject *args);

// Define an array of module methods
static PyMethodDef module_methods[] = {
	{"imag_prop", _test_imag_prop, METH_VARARGS, "An imaginary propagation"},
	{"real_prop", _test_real_prop, METH_VARARGS, "An real propagation"},
	{NULL, NULL, 0, NULL}
};

// Define this module
static struct PyModuleDef _test = {
	PyModuleDef_HEAD_INIT,
	"_test",
	"Python Wrapper of QPROP",
	-1,
	module_methods
};

// Create module object
PyMODINIT_FUNC PyInit__test(void) {
	return PyModule_Create(&_test);
}

// Define module methods
static PyObject *_test_imag_prop(PyObject *self, PyObject *args) {
	int return_code = imag_prop(1, NULL);
	if (return_code != 0) {
		PyErr_SetString(PyExc_Exception, "Failed to execute imaginary propagation.");
		return NULL;
	}

	PyObject *result_tuple = Py_BuildValue("i", return_code);
	return result_tuple;
}

static PyObject *_test_real_prop(PyObject *self, PyObject *args) {
	int return_code = real_prop(1, NULL);
	if (return_code != 0) {
		PyErr_SetString(PyExc_Exception, "Failed to execute real propagation.");
		return NULL;
	}

	PyObject *result_tuple = Py_BuildValue("i", return_code);
	return result_tuple;
}


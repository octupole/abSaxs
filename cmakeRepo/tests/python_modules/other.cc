#include "../../../cmakeRepo/tests/python_modules/other.h"

PyObject * meaning_of_life(PyObject *_module) {
    return PyLong_FromLong(42);
}

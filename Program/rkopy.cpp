#define PYBIND11_DETAILED_ERROR_MESSAGES 1
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include "rko.h"

PYBIND11_MODULE(rkopy, m) {

    m.def("solve", &solve);

}
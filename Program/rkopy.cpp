#include "rko.h"
// #include "teste.h"

#define PYBIND11_DETAILED_ERROR_MESSAGES 1
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

// namespace py = pybind11;


// PYBIND11_MAKE_OPAQUE(std::vector<double>)

PYBIND11_MODULE(rkopy, m) {

    m.def("solve", &solve);
    // m.def("teste", &teste);

}
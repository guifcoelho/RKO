from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension

ext_modules = [
    Pybind11Extension(
        "rkopy",
        ["Program/rkopy.cpp"],
        extra_compile_args=["/openmp"]
    ),
]

setup(name="rkopy", ext_modules=ext_modules)

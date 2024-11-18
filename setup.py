from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension

ext_modules = [
    Pybind11Extension(
        "rkopy",
        ["Program/rkopy.cpp"]
    ),
]

setup(name="rkopy", ext_modules=ext_modules)

# distutils: language = c++
# cython: c_string_type=unicode, c_string_encoding=utf8

cimport RealValuedGOMEA 
from RealValuedGOMEA import pyRealValuedGOMEA
cimport DiscreteGOMEA
from DiscreteGOMEA import pyDiscreteGOMEA

include "Fitness.pyx"

def RealValuedGOMEA(**kwargs):
    return pyRealValuedGOMEA(**kwargs)

def DiscreteGOMEA(**kwargs):
    return pyDiscreteGOMEA(**kwargs)


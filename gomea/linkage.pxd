import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "gomea/src/common/linkage_model.hpp" namespace "gomea":
    cdef cppclass linkage_config_t:
        linkage_config_t() except +
        linkage_config_t(size_t) except +
        linkage_config_t(int,bool,int,int) except +
        linkage_config_t(int,bool,bool) except +
        linkage_config_t(vector[vector[int]]) except +
        linkage_config_t(string) except +

cdef class LinkageModel:
    cdef linkage_config_t *c_inst

cdef class cLinkageTree(LinkageModel):
    pass

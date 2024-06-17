import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "gomea/src/common/linkage_config.hpp" namespace "gomea":
    cdef cppclass linkage_config_t:
        linkage_config_t() except +
        linkage_config_t(int,bool) except +
        linkage_config_t(int,bool,int,bool) except +
        linkage_config_t(int,bool,bool) except +
        linkage_config_t(vector[vector[int]]) except +
        linkage_config_t(string) except +

cdef class LinkageModel:
    cdef linkage_config_t *c_inst

cdef class Univariate(LinkageModel):
    pass

cdef class BlockMarginalProduct(LinkageModel):
    pass

cdef class Full(LinkageModel):
    pass

cdef class LinkageTree(LinkageModel):
    pass

cdef class StaticLinkageTree(LinkageModel):
    pass

cdef class Conditional(LinkageModel):
    pass

cdef class Custom(LinkageModel):
    pass
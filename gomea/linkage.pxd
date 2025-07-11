import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "gomea/src/common/linkage_config.hpp" namespace "gomea":
    cdef cppclass linkage_config_t:
        @staticmethod
        linkage_config_t *constructor_UNI()
        @staticmethod
        linkage_config_t *constructor_MPM(int)
        @staticmethod
        linkage_config_t *constructor_LT(string,bool,int,bool)
        @staticmethod
        linkage_config_t *constructor_COND(int,bool,bool)
        @staticmethod
        linkage_config_t *constructor_CUSTOM(vector[vector[int]])
        @staticmethod
        linkage_config_t *constructor_FILE(string)

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

#cdef class Conditional(LinkageModel):
#    pass

cdef class Custom(LinkageModel):
    pass
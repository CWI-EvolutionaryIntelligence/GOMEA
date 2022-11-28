import numpy as np
cimport numpy as np

cdef extern from "gomea/src/common/linkage_model.hpp" namespace "gomea":
    cdef cppclass linkage_config_t:
        linkage_config_t() except +
        linkage_config_t(size_t) except +

cdef class LinkageModel:
    cdef linkage_config_t *c_inst

cdef class Univariate(LinkageModel):
    pass

cdef class MarginalProductModel(LinkageModel):
    pass

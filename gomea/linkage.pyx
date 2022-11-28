cdef class Univariate(LinkageModel):
    def __cinit__(self):
        self.c_inst = new linkage_config_t()

cdef class MarginalProductModel(LinkageModel):
    def __cinit__(self, block_size : size_t ):
        self.c_inst = new linkage_config_t(block_size)


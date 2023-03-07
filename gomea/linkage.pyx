from libcpp.vector cimport vector
from libcpp cimport bool

cdef class LinkageModel:
    def __dealloc__(self):
        del self.c_inst 

cdef class Univariate(LinkageModel):
    def __cinit__(self):
        self.c_inst = new linkage_config_t()

cdef class MarginalProductModel(LinkageModel):
    def __cinit__(self,
        block_size : size_t
    ):
        self.c_inst = new linkage_config_t(block_size)

cdef class Full(LinkageModel):
    def __cinit__(self):
        self.c_inst = new linkage_config_t(0)

cdef class StaticLinkageTree(LinkageModel):
    def __cinit__(self,
        similarity_measure : int = 2,
        filtered : bool = True,
        maximum_set_size : int = -1
    ):
        cdef bool is_static = True 
        self.c_inst = new linkage_config_t(similarity_measure, filtered, maximum_set_size, is_static)

cdef class LinkageTree(LinkageModel):
    def __cinit__(self,
        similarity_measure : int = 0,
        filtered : bool = False,
        maximum_set_size : int = -1
    ):
        cdef bool is_static = False
        self.c_inst = new linkage_config_t(similarity_measure, filtered, maximum_set_size, is_static)

cdef class Conditional(LinkageModel):
    def __cinit__(self,
        max_clique_size : int = 1,
        cliques_as_fos_elements : bool = True,
        include_full_fos_element : bool = True
    ):
        self.c_inst = new linkage_config_t(max_clique_size,cliques_as_fos_elements,include_full_fos_element)

cdef class Custom(LinkageModel):
    def __cinit__(self,
        vector[vector[int]] FOS
    ):
        self.c_inst = new linkage_config_t(FOS)
        
cdef class FromFile(LinkageModel):
    def __cinit__(self,
        filename: string
    ):
        #f = str.encode(filename)
        cdef string f = filename.decode()
        self.c_inst = new linkage_config_t(f)

def UCondGG():
    return Conditional(1,False,True)

def UCondFG():
    return Conditional(1,True,False)

def UCondHG():
    return Conditional(1,True,True)

def MCondHG(*args,**kwargs):
    if len(args) == 1:
        return Conditional(args[0],True,True)
    elif "max_clique_size" in kwargs:
        return Conditional(kwargs["max_clique_size"],True,True)
    else:
        raise RuntimeError("Invalid arguments.")


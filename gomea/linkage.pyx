from libcpp.vector cimport vector
from libcpp cimport bool

cdef class LinkageModel:
    def __dealloc__(self):
        del self.c_inst 

cdef class Univariate(LinkageModel):
    def __cinit__(self):
        self.c_inst = linkage_config_t.constructor_UNI()

cdef class BlockMarginalProduct(LinkageModel):
    def __cinit__(self,
        block_size : int = -1
    ):
        self.c_inst = linkage_config_t.constructor_MPM(block_size)

cdef class Full(LinkageModel):
    def __cinit__(self):
        self.c_inst = linkage_config_t.constructor_MPM(0)

cdef class StaticLinkageTree(LinkageModel):
    def __cinit__(self,
        similarity_measure = "VIG",
        filtered : bool = True,
        maximum_set_size : int = -1
    ):
        self.c_inst = linkage_config_t.constructor_LT(str.encode(similarity_measure), filtered, maximum_set_size, True)

cdef class LinkageTree(LinkageModel):
    def __cinit__(self,
        similarity_measure = "MI",
        filtered : bool = False,
        maximum_set_size : int = -1
    ):
        self.c_inst = linkage_config_t.constructor_LT(str.encode(similarity_measure), filtered, maximum_set_size, False)

#cdef class Conditional(LinkageModel):
#    def __cinit__(self,
#        max_clique_size : int = 1,
#        cliques_as_fos_elements : bool = True,
#        include_full_fos_element : bool = True
#    ):
#        if not cliques_as_fos_elements and not include_full_fos_element:
#            raise AssertionError("At least one of input parameters 'cliques_as_fos_elements' or 'include_full_fos_element' must be True.")
#        self.c_inst = new linkage_config_t(max_clique_size,cliques_as_fos_elements,include_full_fos_element)

cdef class Custom(LinkageModel):
    def __cinit__(self,
        *args,
        **kwargs
    ):
        cdef string file = str.encode("")
        cdef vector[vector[int]] fos = []
        if 'file' in kwargs:
            file = str.encode(kwargs['file'])
        if 'fos' in kwargs:
            fos = kwargs['fos']
        if(len(file) is 0 and len(fos) is 0):
            raise AssertionError("Constructor requires exactly 1 argument.")
        if(len(file) > 0 and len(fos) > 0):
            raise AssertionError("Constructor requires exactly 1 argument.")
        if( len(file) > 0 ):
            self.c_inst = linkage_config_t.constructor_FILE(file)
        elif( len(fos) > 0 ):
            self.c_inst = linkage_config_t.constructor_CUSTOM(fos)
        else:
            raise AssertionError("Constructor requires exactly 1 argument.")

#def UCondGG():
#    return Conditional(1,False,True)

#def UCondFG():
#    return Conditional(1,True,False)

#def UCondHG():
#    return Conditional(1,True,True)

#def MCondHG(*args,**kwargs):
#    if len(args) == 1:
#        return Conditional(args[0],True,True)
#    elif "max_clique_size" in kwargs:
#        return Conditional(kwargs["max_clique_size"],True,True)
#    else:
#        raise RuntimeError("Invalid arguments.")


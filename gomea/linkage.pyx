from libcpp.vector cimport vector

cdef class cUnivariate(LinkageModel):
    def __cinit__(self):
        self.c_inst = new linkage_config_t()

cdef class cMarginalProductModel(LinkageModel):
    def __cinit__(self,
        block_size : size_t
    ):
        self.c_inst = new linkage_config_t(block_size)

cdef class cFull(LinkageModel):
    def __cinit__(self):
        self.c_inst = new linkage_config_t(0)

cdef class cLinkageTree(LinkageModel):
    def __cinit__(self,
        similarity_measure : int = 0,
        filtered : bool = False,
        maximum_set_size : int = -1
    ):
        print(similarity_measure,filtered,maximum_set_size)
        self.c_inst = new linkage_config_t(similarity_measure, filtered, maximum_set_size, 0)

cdef class cConditional(LinkageModel):
    def __cinit__(self,
        max_clique_size : int = 1,
        cliques_as_fos_elements : bool = True,
        include_full_fos_element : bool = True
    ):
        self.c_inst = new linkage_config_t(max_clique_size,cliques_as_fos_elements,include_full_fos_element)

cdef class cCustom(LinkageModel):
    def __cinit__(self,
        vector[vector[int]] FOS
    ):
        self.c_inst = new linkage_config_t(FOS)
        
cdef class cFromFile(LinkageModel):
    def __cinit__(self,
        filename: string
    ):
        self.c_inst = new linkage_config_t(filename)

def Univariate():
    return cUnivariate()

def MarginalProductModel(*args, **kwargs):
    if len(args) == 1:
        return cMarginalProductModel(*args)
    elif "block_size" in kwargs:
        return cMarginalProductModel(**kwargs) 
    else:
        raise RuntimeError("Invalid arguments.")

def Full():
    return cFull()

def LinkageTree(*args,**kwargs):
    if len(args) == 3:
        return cLinkageTree(*args)
    elif len(args) == 0:
        return cLinkageTree(**kwargs)
    else:
        raise RuntimeError("Invalid arguments.")

def UCondGG():
    return cConditional(1,False,True)

def UCondFG():
    return cConditional(1,True,False)

def UCondHG():
    return cConditional(1,True,True)

def MCondHG(*args,**kwargs):
    if len(args) == 1:
        return cConditional(args[0],True,True)
    elif "max_clique_size" in kwargs:
        return cConditional(kwargs["max_clique_size"],True,True)
    else:
        raise RuntimeError("Invalid arguments.")

def Conditional(*args,**kwargs):
    if len(args) == 3:
        return cConditional(args)
    elif len(args) == 0:
        return cConditional(**kwargs)
    else:
        raise RuntimeError("Invalid arguments.")

def Custom(*args,**kwargs):
    if len(args) == 1:
        return cConditional(args)
    elif "FOS" in kwargs:
        return cCustom(**kwargs)
    else:
        raise RuntimeError("Invalid arguments.")

def FromFile(*args,**kwargs):
    if len(args) == 1:
        return cFromFile(str.encode(args[0]))
    elif "filename" in kwargs:
        kwargs["filename"] = str.encode(kwargs["filename"])
        return cFromFile(**kwargs)
    else:
        raise RuntimeError("Invalid arguments.")


# distutils: language = c++

cdef int fitness_test() except -1:
    print("pyfitness - TEST")
    return 1

cdef public int fitness_embedded() except -1:
    print("pyFitness - EMBED TEST")
    return 1


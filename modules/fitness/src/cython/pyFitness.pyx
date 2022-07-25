# distutils: language = c++

cdef public int fitness_embedded() except -1:
    print("EMBED TEST")
    return 1

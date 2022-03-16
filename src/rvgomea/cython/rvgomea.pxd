cdef extern from "rvgomea.cpp":
    pass

# Declare the class with cdef
cdef extern from "rvgomea.h" namespace "shapes":
    cdef cppclass rvgomea:
        rvgomea() except +
        rvgomea(int, int, int, int) except +
        int x0, y0, x1, y1
        int getArea()
        void getSize(int* width, int* height)
        void move(int, int)

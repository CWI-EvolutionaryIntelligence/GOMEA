from libcpp.map cimport map
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "gomea/src/common/output_statistics.hpp" namespace "gomea":
    cdef cppclass output_statistics_t:
        output_statistics_t() except +
        vector[string] getGenerationalMetricNames() except +
        vector[T] getGenerationalMetricValues[T](string) except +
        vector[string] getFinalMetricNames() except +
        T getFinalMetricValue[T](string) except +

cdef class OutputStatisticsWrapper:
    cdef output_statistics_t *c_ptr
    
    @staticmethod
    cdef OutputStatisticsWrapper from_ptr(output_statistics_t *_ptr)

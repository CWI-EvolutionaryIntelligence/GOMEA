from gomea.output cimport output_statistics_t
from ctypes import * 

cdef class OutputStatisticsWrapper:
    def __cinit__(self):
        pass

    def __init__(self):
        raise TypeError("This class cannot be instantiated directly.")

    @staticmethod
    cdef OutputStatisticsWrapper from_ptr(output_statistics_t *_ptr):
        cdef OutputStatisticsWrapper wrapper = OutputStatisticsWrapper.__new__(OutputStatisticsWrapper)
        wrapper.c_ptr = _ptr
        return wrapper

class OutputStatistics:
    def __init__(self, stats : OutputStatisticsWrapper):
        self.generational_metrics = stats.c_ptr[0].getGenerationalMetricNames()
        self.generational_metrics_dict = {}
        self.final_metrics = stats.c_ptr[0].getFinalMetricNames()
        self.final_metrics_dict = {}

        def getGenerationalMetricValues( metric_name ):
            try:
                return stats.c_ptr[0].getGenerationalMetricValues[int](metric_name)
            except Exception as err:
                pass
            try:
                return stats.c_ptr[0].getGenerationalMetricValues[double](metric_name)
            except Exception as err:
                pass
            try:
                b_result = stats.c_ptr[0].getGenerationalMetricValues[string](metric_name)
                return [byte_string.decode('utf-8') for byte_string in b_result]
            except Exception as err:
                pass

        def getFinalMetricValue( metric_name ):
            try:
                return stats.c_ptr[0].getFinalMetricValue[int](metric_name)
            except Exception as err:
                pass
            try:
                return stats.c_ptr[0].getFinalMetricValue[double](metric_name)
            except Exception as err:
                pass
            try:
                b_result = stats.c_ptr[0].getFinalMetricValue[string](metric_name)
                return b_result.decode('utf-8')
            except Exception as err:
                pass

        for metric_name in self.generational_metrics:
            s_metric_name = metric_name.decode()
            self.generational_metrics_dict[s_metric_name] = getGenerationalMetricValues(metric_name)
        self.generational_metrics = [s.decode() for s in self.generational_metrics]
        for metric_name in self.final_metrics:
            s_metric_name = metric_name.decode()
            self.final_metrics_dict[s_metric_name] = getFinalMetricValue(metric_name)
        self.final_metrics = [s.decode() for s in self.final_metrics]

    def __getitem__(self, metric):
        assert( metric in self.generational_metrics, str(metric)+" is not a valid metric." )
        return self.generational_metrics_dict[metric]

    def getFinalStatistics(self):
        return self.final_metrics_dict.copy()

    def getGenerationalStatistics(self):
        return self.generational_metrics_dict.copy()
    
    def printFinalStatistics(self):
        print("Final statistics:")
        for k,v in self.final_metrics_dict.items():
            print(k,":",v)

    def printGenerationalStatistics(self):
        print("Generational statistics:")
        for k,v in self.generational_metrics_dict.items():
            print(k,":",v)
    
    def printAllStatistics(self):
        self.printGenerationalStatistics()
        self.printFinalStatistics()

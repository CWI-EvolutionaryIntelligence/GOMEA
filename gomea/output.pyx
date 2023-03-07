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
        self.metrics = stats.c_ptr[0].getAllMetricNames()
        self.metrics_map = {}

        def getMetricValues( metric_name ):
            try:
                return stats.c_ptr[0].getMetricValues[int](metric_name)
            except Exception as err:
                pass
            try:
                return stats.c_ptr[0].getMetricValues[double](metric_name)
            except Exception as err:
                pass

        for metric_name in self.metrics:
            s_metric_name = metric_name.decode()
            self.metrics_map[s_metric_name] = getMetricValues(metric_name)
        self.metrics = [s.decode() for s in self.metrics]
        #self.elitist_genotype = stats.c_ptr[0].getElitistGenotype()
        #self.elitist_fitness = stats.c_ptr[0].getElitistFitness()

    def __getitem__(self, metric):
        assert( metric in self.metrics, str(metric)+" is not a valid metric." )
        return self.metrics_map[metric][-1]

    def getFinalStatistics(self):
        final_stats = {}
        for metric in self.metrics:
            final_stats[metric] = self.metrics_map[metric][-1]
        return final_stats
    
    def printFinalStatistics(self):
        final_stats = self.getFinalStatistics()
        print("Final statistics:")
        for k,v in final_stats.items():
            print(k,":",v)
    
    def printAllStatistics(self):
        print("All statistics:")
        for k,v in self.metrics_map.items():
            print(k,":",v)

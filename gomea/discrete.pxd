from libcpp.string cimport string
from libcpp cimport bool
from gomea.output cimport output_statistics_t
from gomea.fitness cimport fitness_t
from gomea.linkage cimport linkage_config_t 

cdef extern from "gomea/src/discrete/Config.hpp" namespace "gomea::discrete":
    cdef cppclass Config:
        Config() except +

        int maximumNumberOfGOMEAs, IMSsubgenerationFactor, basePopulationSize, maxArchiveSize, maximumNumberOfEvaluations, maximumNumberOfGenerations, AnalyzeFOS
        long long randomSeed
        bool fix_seed, generational_statistics, generational_solution, verbose
        output_frequency_t output_frequency
        double maximumNumberOfSeconds
        string folder, problemInstancePath
        fitness_t[char] *fitness
        linkage_config_t *linkage_config

cdef extern from "gomea/src/discrete/gomeaIMS.hpp" namespace "gomea::discrete":
    cdef cppclass gomeaIMS:
        gomeaIMS() except +
        gomeaIMS(Config*) except +
        
        void run() except +
        void runGeneration() except +
        bool checkTermination() except +
        double getProgressUntilTermination() except +
        output_statistics_t output

cdef extern from "src/common/gomea_defs.hpp" namespace "gomea":
    ctypedef enum output_frequency_t:
        GEN, IMS_GEN, NEW_ELITE

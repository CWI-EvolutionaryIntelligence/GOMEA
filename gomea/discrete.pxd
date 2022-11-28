from libcpp.string cimport string
from libcpp cimport bool
from gomea.fitness cimport fitness_t
from gomea.linkage cimport linkage_config_t 

cdef extern from "gomea/src/discrete/Config.hpp" namespace "gomea::discrete":
    cdef cppclass Config:
        Config() except +

        int maximumNumberOfGOMEAs, IMSsubgenerationFactor, basePopulationSize, maxArchiveSize, maximumNumberOfEvaluations, maximumNumberOfGenerations
        long long randomSeed
        bool fix_seed
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

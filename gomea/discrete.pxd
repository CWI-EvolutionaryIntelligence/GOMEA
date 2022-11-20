from libcpp.string cimport string
from libcpp cimport bool
from gomea.fitness cimport fitness_t

cdef extern from "gomea/src/discrete/Config.hpp" namespace "gomea::discrete":
    cdef cppclass Config:
        Config() except +

        int maximumNumberOfGOMEAs, IMSsubgenerationFactor, basePopulationSize, maxArchiveSize, maximumNumberOfEvaluations, maximumNumberOfGenerations
        double maximumNumberOfSeconds
        string folder, problemInstancePath
        fitness_t[char] *fitness

cdef extern from "gomea/src/discrete/gomeaIMS.hpp" namespace "gomea::discrete":
    cdef cppclass gomeaIMS:
        gomeaIMS() except +
        gomeaIMS(Config*) except +
        
        void run() except +
        void runGeneration() except +
        bool checkTermination() except +
        double getProgressUntilTermination() except +

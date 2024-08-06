# distutils: language = c++
# cython: c_string_type=unicode, c_string_encoding=utf8

from gomea.output cimport OutputStatisticsWrapper
from gomea.output import OutputStatistics
from gomea.discrete cimport gomeaIMS, Config
from gomea.linkage cimport LinkageModel, StaticLinkageTree
from libcpp.string cimport string
from libcpp cimport bool
from cpython.exc cimport PyErr_CheckSignals
from tqdm import tqdm
import inspect

from gomea.fitness cimport FitnessFunction, GBOFitnessFunction, BBOFitnessFunction
include "EmbeddedFitness.pxi"

# Create a Cython extension type which holds a C++ instance
# as an attribute and create a bunch of forwarding methods
# Python extension type.
cdef class DiscreteGOMEA:
    cdef gomeaIMS c_inst  # Hold a C++ instance which we're wrapping
    cdef Config c_config

    def __cinit__(self,
        # Optimization problem settings (required)
        fitness: FitnessFunction, 
        # GOMEA settings
        linkage_model : LinkageModel = StaticLinkageTree(),
        max_number_of_populations : int = 25,
        IMS_subgeneration_factor: int = 4,
        base_population_size: int = 2,
        max_number_of_evaluations : int = -1,
        max_number_of_generations : int = -1,
        max_number_of_seconds : float = -1.0,
        random_seed : int = -1,
    ):

        # Initialize attributes 
        #_, _, _, values = inspect.getargvalues(inspect.currentframe())
        #for arg, val in values.items():
            #setattr(self, arg, val)

        self.c_config = Config()
        assert( (<FitnessFunction?>fitness).c_inst_discrete != NULL, "FitnessFunction is not discrete." )
        self.c_config.fitness = (<FitnessFunction?>fitness).c_inst_discrete
        self.c_config.linkage_config = linkage_model.c_inst
        self.c_config.folder = string(b"output_discrete_gomea")
        self.c_config.maximumNumberOfGOMEAs = max_number_of_populations
        self.c_config.IMSsubgenerationFactor = IMS_subgeneration_factor
        self.c_config.basePopulationSize = base_population_size
        self.c_config.problemInstancePath = string(b"") 
        self.c_config.maxArchiveSize = 1000 
        self.c_config.maximumNumberOfEvaluations = max_number_of_evaluations
        self.c_config.maximumNumberOfGenerations = max_number_of_generations
        self.c_config.maximumNumberOfSeconds = max_number_of_seconds
        self.c_config.AnalyzeFOS = 0
        #if analyze_fos:
        #    self.c_config.AnalyzeFOS = 1
        #self.c_config.verbose = False 
        self.c_config.fix_seed = False
        if random_seed != -1:
            self.c_config.randomSeed = random_seed
            self.c_config.fix_seed = True

        # Initialize C++ instance
        self.c_inst = gomeaIMS(&self.c_config)

    def get_progress(self):
        return self.c_inst.getProgressUntilTermination()

    def check_termination(self):
        cdef bool t = self.c_inst.checkTermination()
        return t

    def run_generation(self):
        self.c_inst.runGeneration()

    def init_progress_bar(self):
        progress_bar = tqdm(desc="DiscreteGOMEA",unit='%',total=100)
        return progress_bar

    def update_progress_bar(self, progress_bar):
        progress_bar.n = self.get_progress()
        progress_bar.refresh()
    
    def run(self):
        self.c_inst.run()
        return OutputStatistics(OutputStatisticsWrapper.from_ptr(&self.c_inst.output))

    def run_with_progress(self):
        with self.init_progress_bar() as progress_bar:
            while( not self.check_termination() ):
                self.run_generation()
                self.update_progress_bar(progress_bar)
                #print("PyErr:",PyErr_CheckSignals() )
                if( PyErr_CheckSignals() == -1 ):
                    break
            self.update_progress_bar(progress_bar)
        #self.c_inst.run()

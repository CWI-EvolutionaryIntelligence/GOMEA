# distutils: language = c++
# cython: c_string_type=unicode, c_string_encoding=utf8

from gomea.output cimport OutputStatisticsWrapper
from gomea.output import OutputStatistics
from gomea.discrete cimport gomeaIMS, Config
from gomea.linkage cimport LinkageModel, LinkageTree
from libcpp.string cimport string
from libcpp cimport bool
from cpython.exc cimport PyErr_CheckSignals
from tqdm import tqdm
import inspect

from gomea.fitness cimport FitnessFunction, GBOFitnessFunction
include "gomea/EmbeddedFitness.pxi"

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
        linkage_model : LinkageModel = LinkageTree(),
        folder: string=string(b"output_discrete_gomea"),
        maximum_number_of_GOMEAs: int=25,
        IMS_subgeneration_factor: int=4,
        base_population_size: int=2,
        problem_instance_path: string=string(b""),
        max_archive_size: int=1000,
        max_evals : int=-1,
        max_gens : int=-1,
        max_time : double=-1.0,
        random_seed : int=-1,
        # Other
        verbose : bool=False,
        analyze_fos : bool=False
    ):

        # Initialize attributes 
        #_, _, _, values = inspect.getargvalues(inspect.currentframe())
        #for arg, val in values.items():
            #setattr(self, arg, val)

        self.c_config = Config()
        self.c_config.fitness = (<FitnessFunction?>fitness).c_inst_discrete
        self.c_config.linkage_config = linkage_model.c_inst
        self.c_config.folder = folder
        self.c_config.maximumNumberOfGOMEAs = maximum_number_of_GOMEAs
        self.c_config.IMSsubgenerationFactor = IMS_subgeneration_factor
        self.c_config.basePopulationSize = base_population_size
        self.c_config.problemInstancePath = problem_instance_path
        self.c_config.maxArchiveSize = max_archive_size
        self.c_config.maximumNumberOfEvaluations = max_evals
        self.c_config.maximumNumberOfGenerations = max_gens
        self.c_config.maximumNumberOfSeconds = max_time
        self.c_config.AnalyzeFOS = 0
        if analyze_fos:
            self.c_config.AnalyzeFOS = 1
        #self.c_config.verbose = verbose
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

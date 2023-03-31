# distutils: language = c++
# cython: c_string_type=unicode, c_string_encoding=utf8

from gomea.real_valued cimport rvg_t, Config
from gomea.output cimport OutputStatisticsWrapper
from gomea.output import OutputStatistics
from libcpp.string cimport string
from libcpp cimport bool
import inspect

from gomea.fitness cimport FitnessFunction, GBOFitnessFunction
from gomea.linkage cimport LinkageModel, StaticLinkageTree
include "gomea/EmbeddedFitness.pxi"

# Create a Cython extension type which holds a C++ instance
# as an attribute and create a bunch of forwarding methods
# Python extension type.
cdef class RealValuedGOMEA:
    cdef rvg_t c_inst  # Hold a C++ instance which we're wrapping
    cdef Config c_config

    def __cinit__(self,
        # Optimization problem settings (required)
        fitness: FitnessFunction, 
        # GOMEA parameters
        linkage_model : LinkageModel = StaticLinkageTree(),
        lower_init_range: double=0.0,
        upper_init_range: double=1.0,
        random_seed: int=-1,
        # IMS settings (optional)
        max_number_of_populations: int=25,
        max_number_of_generations: int=-1,
        IMS_subgeneration_factor: int=8,
        base_population_size: int=10,
        # Termination settings (optional)
        max_number_of_evaluations : double=-1,
        max_number_of_seconds : double=-1,
    ):
        
        # Initialize attributes 
        self.c_config = Config()
        self.c_config.problem_index = 0 
        self.c_config.fitness = (<FitnessFunction?>fitness).c_inst_realvalued
        self.c_config.linkage_config = linkage_model.c_inst
        self.c_config.use_vtr = True 
        self.c_config.lower_user_range = lower_init_range
        self.c_config.upper_user_range = upper_init_range
        self.c_config.tau = 0.35
        self.c_config.distribution_multiplier_decrease = 0.9
        self.c_config.st_dev_ratio_threshold = 1.0
        self.c_config.fitness_variance_tolerance = 0.0
        self.c_config.maximum_no_improvement_stretch = 100
        self.c_config.selection_during_gom = True
        self.c_config.update_elitist_during_gom = True 
        self.c_config.maximum_number_of_populations = max_number_of_populations
        self.c_config.number_of_subgenerations_per_population_factor = IMS_subgeneration_factor
        self.c_config.base_population_size = base_population_size
        self.c_config.maximum_number_of_generations = max_number_of_generations
        self.c_config.maximum_number_of_evaluations = max_number_of_evaluations
        self.c_config.maximum_number_of_seconds = max_number_of_seconds
        self.c_config.black_box_evaluations = False
        self.c_config.write_generational_statistics = False
        self.c_config.write_generational_solutions = False
        self.c_config.verbose = False
        self.c_config.print_verbose_overview = False 
        self.c_config.fix_seed = 0
        if random_seed != -1:
            self.c_config.random_seed = random_seed
            self.c_config.fix_seed = 1

        # Initialize C++ instance
        self.c_inst = rvg_t(&self.c_config)
    
    def print_usage(self):
        self.c_inst.printUsage()

    def run(self):
        self.c_inst.run()
        return OutputStatistics(OutputStatisticsWrapper.from_ptr(&self.c_inst.output))

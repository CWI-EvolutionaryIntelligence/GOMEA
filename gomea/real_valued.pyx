# distutils: language = c++
# cython: c_string_type=unicode, c_string_encoding=utf8

from gomea.real_valued cimport rvg_t, Config
from gomea.output cimport OutputStatisticsWrapper
from gomea.output import OutputStatistics
from libcpp.string cimport string
from libcpp cimport bool
import inspect

from gomea.fitness cimport FitnessFunction, GBOFitnessFunction, BBOFitnessFunction
from gomea.linkage cimport LinkageModel, StaticLinkageTree
include "EmbeddedFitness.pxi"

# Create a Cython extension type which holds a C++ instance
# as an attribute and create a bunch of forwarding methods
# Python extension type.
cdef class RealValuedGOMEA:
    cdef rvg_t c_inst  # Hold a C++ instance which we're wrapping
    cdef Config c_config
    cdef FitnessFunction fitness

    def __cinit__(self,
        # Optimization problem settings (required)
        fitness: FitnessFunction, 
        # GOMEA settings (optional)
        linkage_model : LinkageModel = StaticLinkageTree(),
        lower_init_range: float = 0.0,
        upper_init_range: float = 1.0,
        random_seed: int=-1,
        # IMS settings (optional)
        max_number_of_populations: int = 25,
        IMS_subgeneration_factor: int = 8,
        base_population_size: int = 10,
        # Termination settings (optional)
        max_number_of_generations: int = -1,
        max_number_of_evaluations : float = -1,
        max_number_of_seconds : float = -1,
        # Output settings (optional)
        generational_statistics = True,
        generational_solution = False,
        #output_frequency = 'GEN',
        verbose = False,
    ):
        
        # Initialize attributes 
        self.c_config = Config()
        self.c_config.problem_index = 0
        assert( (<FitnessFunction?>fitness).c_inst_realvalued != NULL, "FitnessFunction is not real-valued." )
        self.c_config.fitness = (<FitnessFunction?>fitness).c_inst_realvalued
        # as we only store a pointer to a field of the FitnessFunction (which actually manages the memory around the pointer)
        # there is the risk of the FitnessFunction being removed & its contents deallocated. This leaves the above pointer dangling.
        # This is particularly likely if the FitnessFunction was constructed as an argument to this class -
        # in which case it & corresponding memory will be deallocated after this initialization procedure has completed.
        # Fix this by keeping the reference to the FitnessFunction around for as long as this instance exists.
        self.fitness = fitness

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
        self.c_config.generational_statistics = generational_statistics 
        self.c_config.generational_solution = generational_solution
        self.c_config.verbose = verbose
        self.c_config.print_verbose_overview = False 
        output_frequency = 'GEN'
        if output_frequency == 'GEN':
            self.c_config.output_frequency = GEN
        elif output_frequency == 'IMS_GEN':
            self.c_config.output_frequency = IMS_GEN
        elif output_frequency == 'NEW_ELITE':
            self.c_config.output_frequency = NEW_ELITE
        else:
            raise ValueError("output_frequency must be one of 'GEN', 'IMS_GEN', 'NEW_ELITE'")
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

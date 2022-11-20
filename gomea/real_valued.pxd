from libcpp.string cimport string
from libcpp cimport bool
from gomea.fitness cimport fitness_t

cdef extern from "src/real_valued/Config.hpp" namespace "gomea::realvalued":
    cdef cppclass Config:
        Config() except +

        int problem_index, maximum_no_improvement_stretch, FOS_element_size, maximum_number_of_populations, number_of_subgenerations_per_population_factor, base_population_size 
        double lower_user_range, upper_user_range, tau, distribution_multiplier_decrease, st_dev_ratio_threshold, fitness_variance_tolerance, maximum_number_of_evaluations, maximum_number_of_seconds, vtr
        short use_conditional_sampling, fix_seed, use_vtr, write_generational_statistics, write_generational_solutions, black_box_evaluations, selection_during_gom, update_elitist_during_gom, verbose, print_verbose_overview
        fitness_t[double] *fitness

cdef extern from "src/real_valued/rv-gomea.hpp" namespace "gomea::realvalued":
    cdef cppclass rvg_t:
        rvg_t() except +
        rvg_t(Config*) except +
        
        void printUsage() except +

        void run() except +
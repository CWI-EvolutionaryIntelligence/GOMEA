from libcpp.string cimport string
from libcpp cimport bool
from libcpp.memory cimport shared_ptr
from gomea.output cimport output_statistics_t
from gomea.fitness cimport fitness_t
from gomea.linkage cimport linkage_config_t

cdef extern from "src/real_valued/Config.hpp" namespace "gomea::realvalued":
    cdef cppclass Config:
        Config() except +

        int problem_index, maximum_no_improvement_stretch, maximum_number_of_populations, number_of_subgenerations_per_population_factor, base_population_size, maximum_number_of_generations 
        double lower_user_range, upper_user_range, tau, distribution_multiplier_decrease, st_dev_ratio_threshold, fitness_variance_tolerance, maximum_number_of_evaluations, maximum_number_of_seconds 
        bool fix_seed, use_vtr, write_generational_statistics, write_generational_solutions, black_box_evaluations, selection_during_gom, update_elitist_during_gom, verbose, print_verbose_overview
        long long random_seed
        linkage_config_t *linkage_config
        shared_ptr[fitness_t[double]] fitness

cdef extern from "src/real_valued/rv-gomea.hpp" namespace "gomea::realvalued":
    cdef cppclass rvg_t:
        rvg_t() except +
        rvg_t(Config*) except +
        
        void run() except +
        
        output_statistics_t output

#pragma once

#include <vector>
#include <set>
#include <string>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <deque>
#include <random>
#include <memory>

#include "gomea/src/common/solution.hpp"
#include "gomea/src/utils/time.hpp"
#include "gomea/src/utils/tools.hpp"

namespace gomea{

namespace linkage
{
typedef enum{
	UNIVARIATE = 0,
    FULL = 1,
    MPM = 2, // marginal product model
	LINKAGE_TREE = 3,
    CONDITIONAL = 4,
    FROM_FILE = 5,
    CUSTOM_LM = 6
} linkage_model_type;

typedef enum{
    MI = 0, // mutual information
    NMI = 1, // normalized mutual information
    VIG = 2, // variable interaction graph
    TIGHT = 3,  // tightly linked -- based on distance in genotype
    RANDOM = 4, // randomly uniformly distributed
} similarity_measure_type;
}

class linkage_config_t{
    public:
        linkage_config_t();
        linkage_config_t( int block_size_ );
        linkage_config_t( std::string similarityMeasure, bool filtered, int maximumSetSize, bool is_static );
        linkage_config_t( int max_clique_size_, bool include_cliques_as_fos_elements_, bool include_full_fos_element_);
        linkage_config_t( const vec_t<vec_t<int>> &FOS );
        linkage_config_t( std::string filename );

        // Required for Cython to distinguish between different constructors
        static inline linkage_config_t* constructor_UNI() {
            return new linkage_config_t();
        }
        static inline linkage_config_t* constructor_MPM( int block_size_ ) {
            return new linkage_config_t( block_size_ );
        }
        static inline linkage_config_t* constructor_LT( std::string similarityMeasure, bool filtered, int maximumSetSize, bool is_static ) {
            return new linkage_config_t( similarityMeasure, filtered, maximumSetSize, is_static );
        }
        static inline linkage_config_t* constructor_COND(int max_clique_size_, bool include_cliques_as_fos_elements_, bool include_full_fos_element_) {
            return new linkage_config_t(max_clique_size_, include_cliques_as_fos_elements_, include_full_fos_element_);
        }
        static inline linkage_config_t* constructor_CUSTOM(const vec_t<vec_t<int>> &FOS) {
            return new linkage_config_t(FOS);
        }
        static inline linkage_config_t* constructor_FILE( std::string filename_ ) {
            return new linkage_config_t( filename_ );
        }
    
        linkage::linkage_model_type type;
        linkage::similarity_measure_type lt_similarity_measure = linkage::similarity_measure_type::MI;

        bool is_mpm = false;
        int mpm_block_size = -1;
        
        bool lt_filtered = false;
        int lt_maximum_set_size = -1;
        bool lt_is_static = false;

		int cond_max_clique_size = 1;
        bool cond_include_cliques_as_fos_elements = true,
             cond_include_full_fos_element = true;

        std::string filename = "";
        
        vec_t<vec_t<int>> FOS;

        static linkage::similarity_measure_type parseSimilarityMeasure(std::string _similarityMeasure)
        {
            std::string similarityMeasure(_similarityMeasure);
            similarityMeasure = utils::toUpper(similarityMeasure);
            if( similarityMeasure == "MI" )
                return linkage::similarity_measure_type::MI;
            else if( similarityMeasure == "NMI" )
                return linkage::similarity_measure_type::NMI;
            else if( similarityMeasure == "VIG" || similarityMeasure == "WVIG")
                return linkage::similarity_measure_type::VIG;
            else if( similarityMeasure == "TIGHT" )
                return linkage::similarity_measure_type::TIGHT;
            else if( similarityMeasure == "RANDOM" )
                return linkage::similarity_measure_type::RANDOM;
            else
                throw std::runtime_error("Unknown similarity measure: " + similarityMeasure);
        }
};

}
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
    MPM = 2,
	LINKAGE_TREE = 3,
    CONDITIONAL = 4,
    FROM_FILE = 5,
    CUSTOM_LM = 6
} linkage_model_type;
}

class linkage_config_t{
    public:
        linkage_config_t();
        linkage_config_t( bool is_mpm_, int block_size_ ); // needs bool parameter to distinguish from FOS matrix constructor (in Cython)
        linkage_config_t( int similarityMeasure, bool filtered, int maximumSetSize, bool is_static );
        linkage_config_t( int max_clique_size_, bool include_cliques_as_fos_elements_, bool include_full_fos_element_);
        linkage_config_t( const vec_t<vec_t<int>> &FOS );
        linkage_config_t( std::string filename );
    
        linkage::linkage_model_type type;

        bool is_mpm = false;
        int mpm_block_size = -1;
        
        int lt_similarity_measure = 0;
        bool lt_filtered = false;
        int lt_maximum_set_size = -1;
        bool lt_is_static = false;

		int cond_max_clique_size = 1;
        bool cond_include_cliques_as_fos_elements = true,
             cond_include_full_fos_element = true;

        std::string filename = "";
        
        vec_t<vec_t<int>> FOS;
};

}
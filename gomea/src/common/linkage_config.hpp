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
	UNIVARIATE,
    FULL,
    MPM,
	LINKAGE_TREE,
    CONDITIONAL,
    FROM_FILE,
    CUSTOM_LM
} linkage_model_type;
}

class linkage_config_t{
    public:
        linkage_config_t();
        linkage_config_t( size_t block_size ); 
        linkage_config_t( int similarityMeasure, bool filtered, int maximumSetSize, bool is_static );
        linkage_config_t( int max_clique_size_, bool include_cliques_as_fos_elements_, bool include_full_fos_element_);
        linkage_config_t( const vec_t<vec_t<int>> &FOS );
        linkage_config_t( std::string filename );
    
        linkage::linkage_model_type type;

        size_t mpm_block_size = -1;
        
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
#include "gomea/src/common/linkage_config.hpp"

namespace gomea{
        
linkage_config_t::linkage_config_t()
{
	type = linkage::UNIVARIATE;
}

linkage_config_t::linkage_config_t( int block_size_ ) : mpm_block_size(block_size_)
{
	type = linkage::MPM;
	is_mpm = true;
}

linkage_config_t::linkage_config_t(std::string similarityMeasure_, bool filtered_, int maximumSetSize_, bool is_static ) 
	: lt_filtered(filtered_), lt_maximum_set_size(maximumSetSize_), lt_is_static(is_static)
{
	type = linkage::LINKAGE_TREE;
	lt_similarity_measure = parseSimilarityMeasure(similarityMeasure_);
}

linkage_config_t::linkage_config_t( const vec_t<vec_t<int>> &FOS_ ) : FOS(FOS_)
{
	type = linkage::CUSTOM_LM;
}

linkage_config_t::linkage_config_t( std::string filename_ ) : filename(filename_)
{
	type = linkage::FROM_FILE;
}

linkage_config_t::linkage_config_t( int max_clique_size_, bool include_cliques_as_fos_elements_, bool include_full_fos_element_)
	: cond_max_clique_size(max_clique_size_), cond_include_cliques_as_fos_elements(include_cliques_as_fos_elements_), cond_include_full_fos_element(include_full_fos_element_)
{
	type = linkage::CONDITIONAL;
}

} // namespace gomea
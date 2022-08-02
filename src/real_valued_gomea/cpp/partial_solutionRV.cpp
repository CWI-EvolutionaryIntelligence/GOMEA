#include "partial_solutionRV.hpp"

namespace gomea{
namespace realvalued{

template<class T>
partial_solution_t<T>::partial_solution_t( vec_t<T> &touched_variables, vec_t<T> &sample_zs, vec_t<int> &touched_indices ) : partial_solution_t(touched_variables,touched_indices)
{
	this->sample_zs = sample_zs;
}

template<class T>
void partial_solution_t<T>::setSampleMean( vec_t<T> &means )
{
	this->sample_means = means;
}

template class partial_solution_t<float>;
template class partial_solution_t<double>;

}}

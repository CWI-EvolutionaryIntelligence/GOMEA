#pragma once

#include <stdio.h>
#include <math.h>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <cassert>
#include <stdexcept>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iostream> 
#include <random>
#include <variant>
#include <memory>
#include <iomanip>

#define MY_PI 3.141592653589793238462643383279
namespace gomea{

template<class T>
using vec_t = std::vector<T>;
template<class T>
using vec_pt = std::shared_ptr<vec_t<T>>;
typedef std::variant<char, int, float, double> genotype_t;
typedef std::map<int,std::set<int>> graph_t;

}

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
#include <iomanip>

#define MY_PI 3.141592653589793238462643383279
namespace gomea{

template<class T>
using vec_t = std::vector<T>;
typedef std::variant<char, int, float, double> genotype_t;
typedef std::map<int,std::set<int>> graph_t;
typedef enum output_frequency_t { GEN, IMS_GEN, NEW_ELITE } output_frequency_t;

}

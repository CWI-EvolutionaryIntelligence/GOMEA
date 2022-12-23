#pragma once

#include <stdio.h>
#include <math.h>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <cassert>
#include <stdexcept>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iostream> 
#include <random>
#include <variant>

namespace gomea{

template<class T>
using vec_t = std::vector<T>;
typedef std::variant<char, int, float, double> genotype_t;

}

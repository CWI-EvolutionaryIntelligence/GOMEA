#pragma once

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <vector>
#include <algorithm>
#include <set>
#include <unordered_map>
#include <map>
#include <memory>
#include <cassert>
#include <limits>

#include "gomea/src/utils/tools.hpp"
#include "gomea/src/utils/linalg.hpp"

namespace gomea{
namespace realvalued{

int *getRanks(double *array, int array_size );
int *getRanksFromSorted(int *sorted, int array_size );

vecE random1DNormalUnitVector( int length );

}}

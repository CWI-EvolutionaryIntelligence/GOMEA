#pragma once

#include <string> 
#include <iostream>

#include "gomea/src/common/output_statistics.hpp"
#include "gomea/src/discrete/Config.hpp"

namespace gomea{
namespace discrete{

class GOMEA
{
public:
    virtual void run() = 0;
    virtual ~GOMEA(){};

    double readVTR(Config *config);

	output_statistics_t output;
    Config *config;
};

}}
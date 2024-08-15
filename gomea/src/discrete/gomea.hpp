#pragma once

#include <string> 
#include <iostream>

#include "gomea/src/discrete/config.hpp"

namespace gomea{
namespace discrete{

class GOMEA
{
public:
    virtual void run() = 0;
    virtual ~GOMEA(){};

    double readVTR(config_t *config);
};

}}
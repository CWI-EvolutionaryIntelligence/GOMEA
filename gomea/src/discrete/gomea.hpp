#pragma once

#include <string> 
#include <iostream>

#include "gomea/src/discrete/Config.hpp"

namespace gomea{
namespace discrete{

class GOMEA
{
public:
    virtual void run() = 0;
    virtual ~GOMEA(){};

    double readVTR(Config *config);
};

}}
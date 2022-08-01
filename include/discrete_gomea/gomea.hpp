#pragma once

#include <string> 
#include <iostream>
using namespace std;

#include <Config.hpp>

class GOMEA
{
public:
    virtual void run() = 0;
    virtual ~GOMEA(){};

    double readVTR(Config *config);
};
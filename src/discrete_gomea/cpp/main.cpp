#include "discrete_gomea/Config.hpp"
#include "discrete_gomea/gomea.hpp"
#include "discrete_gomea/gomeaIMS.hpp"


int main(int argc, char **argv)
{
    Config *config = new Config();
    config->parseCommandLine(argc, argv);
    config->checkOptions();
    config->printOverview();

    config->rng.seed(config->randomSeed);

    GOMEA *gomeaInstance = new gomeaIMS(config);

    try
    {
        gomeaInstance->run();
    }
    catch (customException &ex)
    {}

    delete gomeaInstance;
    delete config;
    
    return 0;
}
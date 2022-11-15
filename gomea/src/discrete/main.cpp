#include "gomea/src/discrete/Config.hpp"
#include "gomea/src/discrete/gomea.hpp"
#include "gomea/src/discrete/gomeaIMS.hpp"

namespace gomea{
namespace discrete{

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

}}
#include "gomea/src/discrete/config.hpp"
#include "gomea/src/discrete/gomea.hpp"
#include "gomea/src/discrete/gomeaIMS.hpp"

using gomea::discrete::config_t;
using gomea::discrete::GOMEA;
using gomea::discrete::gomeaIMS;
int main(int argc, char **argv)
{
    config_t *config = new config_t();
    config->parseCommandLine(argc, argv);
    config->checkOptions();
    config->printOverview();

    GOMEA *gomeaInstance = new gomeaIMS(config);

    try
    {
        gomeaInstance->run();
    }
    catch (gomea::utils::terminationException &ex)
    {}

    delete gomeaInstance;
    delete config;
    
    return 0;
}

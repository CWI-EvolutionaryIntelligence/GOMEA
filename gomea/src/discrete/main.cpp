#include "gomea/src/discrete/Config.hpp"
#include "gomea/src/discrete/gomea.hpp"
#include "gomea/src/discrete/gomeaIMS.hpp"

using namespace gomea;
int main(int argc, char **argv)
{
    discrete::Config *config = new discrete::Config();
    config->parseCommandLine(argc, argv);
    config->checkOptions();
    config->printOverview();

    discrete::GOMEA *gomeaInstance = new discrete::gomeaIMS(config);

    try
    {
        gomeaInstance->run();
    }
    catch (utils::customException &ex)
    {}

    delete gomeaInstance;
    delete config;
    
    return 0;
}

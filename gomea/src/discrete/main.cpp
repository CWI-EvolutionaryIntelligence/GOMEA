#include "gomea/src/discrete/Config.hpp"
#include "gomea/src/discrete/gomea.hpp"
#include "gomea/src/discrete/gomeaIMS.hpp"

using gomea::discrete::Config;
using gomea::discrete::GOMEA;
using gomea::discrete::gomeaIMS;
int main(int argc, char **argv)
{
    Config *config = new Config();
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

    gomeaInstance->output.writeToFile("out.dat");

    delete gomeaInstance;
    delete config;
    
    return 0;
}

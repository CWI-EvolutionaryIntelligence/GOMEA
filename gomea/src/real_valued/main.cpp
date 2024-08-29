#include "gomea/src/real_valued/Config.hpp"
#include "gomea/src/real_valued/rv-gomea.hpp"

using gomea::realvalued::Config;
using gomea::realvalued::rvg_t;
int main(int argc, char **argv)
{
    //Config *config = new Config();
    //config->parseCommandLine(argc, argv);
    //config->checkOptions();
    //config->printOverview();

    rvg_t *gomeaInstance = new rvg_t(argc, argv);

    try
    {
        gomeaInstance->run();
    }
    catch (gomea::utils::customException &ex)
    {}

    delete gomeaInstance;
    //delete config;
    
    return 0;
}

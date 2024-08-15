#include "gomea/src/discrete/gomea.hpp"

namespace gomea{
namespace discrete{

double GOMEA::readVTR(config_t *config)
{
    std::string filename = config->folder + "/vtr.txt";
    std::ifstream inFile;
    inFile.open(filename);

    std::string vtr_str;
    inFile >> vtr_str;
    double vtr = std::stod(vtr_str);

    inFile.close();

    return vtr;
}

}}
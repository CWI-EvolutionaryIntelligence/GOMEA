#include "gomea/src/discrete/gomea.hpp"

namespace gomea{
namespace discrete{

double GOMEA::readVTR(Config *config)
{
    string filename = config->folder + "/vtr.txt";
    ifstream inFile;
    inFile.open(filename);

    string vtr_str;
    inFile >> vtr_str;
    double vtr = stod(vtr_str);

    inFile.close();

    return vtr;
}

}}
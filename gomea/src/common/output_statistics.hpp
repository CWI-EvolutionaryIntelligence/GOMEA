#pragma once

#include "gomea/src/common/gomea_defs.hpp"
#include "gomea/src/common/solution.hpp"

namespace gomea{

typedef std::variant<int, double> metric_t;

class output_statistics_t
{
    public:
        output_statistics_t(){};
        
        void addMetricValue( std::string metric_name, int key, int value );
        void addMetricValue( std::string metric_name, int key, double value );
        template<class T>
        void addMetricValueGeneric( std::string metric_name, int key, T value );
        
        vec_t<std::string> getAllMetricNames();
        template<class T>
        vec_t<T> getMetricValues( std::string metric_name );
        template<class T>
        vec_t<T> getMetricValuesForKeys( std::string metric_name, vec_t<int> input_keys );
        template<class T>
        T getMetricValue( std::string metric_name, int key );

        void writeToFile( std::string filename = "statistics.dat" );
        
        std::map<std::string, std::map<int,metric_t>> metrics_map;
        std::set<int> all_keys;
        solution_t<genotype_t> *elitist_solution;
};

}
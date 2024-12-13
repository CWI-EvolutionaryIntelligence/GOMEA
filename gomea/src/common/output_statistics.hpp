#pragma once

#include "gomea/src/common/gomea_defs.hpp"
#include "gomea/src/common/solution.hpp"

namespace gomea{

typedef std::variant<int, double, std::string> metric_t;

class output_statistics_t
{
    public:
        output_statistics_t(){};
        
        void addGenerationalMetricValue( std::string metric_name, int key, int value );
        void addGenerationalMetricValue( std::string metric_name, int key, double value );
        void addGenerationalMetricValue( std::string metric_name, int key, std::string value );
        template<class T>
        void addGenerationalMetricValueGeneric( std::string metric_name, int key, T value );
        
        vec_t<std::string> getGenerationalMetricNames();
        template<class T>
        vec_t<T> getGenerationalMetricValues( std::string metric_name );
        template<class T>
        vec_t<T> getGenerationalMetricValuesForKeys( std::string metric_name, vec_t<int> input_keys );
        template<class T>
        T getGenerationalMetricValue( std::string metric_name, int key );

        void addFinalMetricValue( std::string metric_name, int value );
        void addFinalMetricValue( std::string metric_name, double value );
        void addFinalMetricValue( std::string metric_name, std::string value );
        template<class T>
        void addFinalMetricValueGeneric( std::string metric_name, T value );
        
        vec_t<std::string> getFinalMetricNames();
        template<class T>
        T getFinalMetricValue( std::string metric_name );

        void writeToFile( std::string filename = "statistics.dat" );
        
        std::unordered_map<std::string, std::unordered_map<int,metric_t>> generational_metrics_map;
        std::unordered_map<std::string, metric_t> final_metrics_map;
        std::set<int> generational_keys;
        solution_t<genotype_t> *elitist_solution;
        int number_of_writes = 0;
};

}
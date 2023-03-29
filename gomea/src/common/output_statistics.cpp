#include "gomea/src/common/output_statistics.hpp"

namespace gomea{

void output_statistics_t::addMetricValue( std::string metric_name, int key, int value ){addMetricValueGeneric<int>(metric_name,key,value);}
void output_statistics_t::addMetricValue( std::string metric_name, int key, double value ){addMetricValueGeneric<double>(metric_name,key,value);}
template<class T>
void output_statistics_t::addMetricValueGeneric( std::string metric_name, int key, T value )
{
    metric_t metric = value;
    if( metrics_map.find(metric_name) == metrics_map.end() )
    {
        metrics_map[metric_name] = std::unordered_map<int,metric_t>();
    }
    metrics_map[metric_name][key] = metric;
    all_keys.insert(key);
}

template<class T>
T output_statistics_t::getMetricValue( std::string metric_name, int key )
{
    assert( metrics_map.find(metric_name) != metrics_map.end() );
    assert(metrics_map[metric_name].find(key) != metrics_map[metric_name].end());
    return( std::get<T>(metrics_map[metric_name][key]) );   
}
template int output_statistics_t::getMetricValue( std::string metric_name, int key );
template double output_statistics_t::getMetricValue( std::string metric_name, int key );
template<>
metric_t output_statistics_t::getMetricValue( std::string metric_name, int key )
{
    assert( metrics_map.find(metric_name) != metrics_map.end() );
    assert( metrics_map[metric_name].find(key) != metrics_map[metric_name].end() );
    return( metrics_map[metric_name][key] );   
}

template<class T>
vec_t<T> output_statistics_t::getMetricValues( std::string metric_name ) 
{
    assert( metrics_map.find(metric_name) != metrics_map.end() );
    vec_t<T> metric_vec;
    for (int k : all_keys)
    {
        assert(metrics_map[metric_name].find(k) != metrics_map[metric_name].end());
        metric_vec.push_back(getMetricValue<T>(metric_name, k));
    }
    return( metric_vec );
}
template vec_t<int> output_statistics_t::getMetricValues( std::string metric_name ); 
template vec_t<double> output_statistics_t::getMetricValues( std::string metric_name ); 

template<class T>
vec_t<T> output_statistics_t::getMetricValuesForKeys( std::string metric_name, vec_t<int> input_keys )
{
    assert( metrics_map.find(metric_name) != metrics_map.end() );
    vec_t<T> metric_vec;
    std::sort(input_keys.begin(), input_keys.end());
    for (int k : input_keys)
    {
        assert(metrics_map[metric_name].find(k) != metrics_map[metric_name].end());
        metric_vec.push_back(getMetricValue<T>(metric_name, k));
    }
    return( metric_vec );
}
template vec_t<int> output_statistics_t::getMetricValuesForKeys( std::string metric_name, vec_t<int> input_keys );
template vec_t<double> output_statistics_t::getMetricValuesForKeys( std::string metric_name, vec_t<int> input_keys );

vec_t<std::string> output_statistics_t::getAllMetricNames()
{
    vec_t<std::string> metrics;
    for( auto const &item : metrics_map )
    {
        metrics.push_back(item.first);
    }
    return metrics;
}

void output_statistics_t::writeToFile(std::string filename) 
{
    char    string[1000];
    FILE *file = fopen( filename.c_str(), "w" );

    if( file != NULL )
    {
        sprintf(string, "# ");
        fputs(string, file);
        for (auto const &v : metrics_map)
        {
            sprintf( string, "%20s ", v.first.c_str());
            fputs(string, file);
        }
        sprintf(string, "\n");
        fputs(string, file);

        for( int k : all_keys )
        {
            for (auto const &item : metrics_map)
            {
                auto vals = item.second;
                std::string strval = std::visit([](auto&& arg){return std::to_string(arg);}, vals[k]);
                sprintf(string, "%20s ", strval.c_str());
                fputs(string, file);
            }
            sprintf(string, "\n");
            fputs(string, file);
        }

        fclose(file);
    }
}

}
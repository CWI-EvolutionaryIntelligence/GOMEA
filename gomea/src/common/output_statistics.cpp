#include "gomea/src/common/output_statistics.hpp"

namespace gomea{

void output_statistics_t::addGenerationalMetricValue( std::string metric_name, int key, int value) { addGenerationalMetricValueGeneric<int>(metric_name, key, value); }
void output_statistics_t::addGenerationalMetricValue( std::string metric_name, int key, double value ) { addGenerationalMetricValueGeneric<double>(metric_name,key,value); }
void output_statistics_t::addGenerationalMetricValue( std::string metric_name, int key, std::string value ) { addGenerationalMetricValueGeneric<std::string>(metric_name,key,value); }
template<class T>
void output_statistics_t::addGenerationalMetricValueGeneric( std::string metric_name, int key, T value )
{
    metric_t metric = value;
    if( generational_metrics_map.find(metric_name) == generational_metrics_map.end() )
    {
        generational_metrics_map[metric_name] = std::unordered_map<int,metric_t>();
    }
    generational_metrics_map[metric_name][key] = metric;
    generational_keys.insert(key);
}

template<class T>
T output_statistics_t::getGenerationalMetricValue( std::string metric_name, int key )
{
    assert( generational_metrics_map.find(metric_name) != generational_metrics_map.end() );
    assert(generational_metrics_map[metric_name].find(key) != generational_metrics_map[metric_name].end());
    return( std::get<T>(generational_metrics_map[metric_name][key]) );   
}
template int output_statistics_t::getGenerationalMetricValue( std::string metric_name, int key );
template double output_statistics_t::getGenerationalMetricValue( std::string metric_name, int key );
template std::string output_statistics_t::getGenerationalMetricValue( std::string metric_name, int key );
template<>
metric_t output_statistics_t::getGenerationalMetricValue( std::string metric_name, int key )
{
    assert( generational_metrics_map.find(metric_name) != generational_metrics_map.end() );
    assert( generational_metrics_map[metric_name].find(key) != generational_metrics_map[metric_name].end() );
    return( generational_metrics_map[metric_name][key] );   
}

template<class T>
vec_t<T> output_statistics_t::getGenerationalMetricValues( std::string metric_name ) 
{
    assert( generational_metrics_map.find(metric_name) != generational_metrics_map.end() );
    vec_t<T> metric_vec;
    for ( int k : generational_keys )
    {
        assert(generational_metrics_map[metric_name].find(k) != generational_metrics_map[metric_name].end());
        metric_vec.push_back(getGenerationalMetricValue<T>(metric_name, k));
    }
    return( metric_vec );
}
template vec_t<int> output_statistics_t::getGenerationalMetricValues( std::string metric_name ); 
template vec_t<double> output_statistics_t::getGenerationalMetricValues( std::string metric_name );
template vec_t<std::string> output_statistics_t::getGenerationalMetricValues( std::string metric_name );

template<class T>
vec_t<T> output_statistics_t::getGenerationalMetricValuesForKeys( std::string metric_name, vec_t<int> input_keys )
{
    assert( generational_metrics_map.find(metric_name) != generational_metrics_map.end() );
    vec_t<T> metric_vec;
    std::sort(input_keys.begin(), input_keys.end());
    for (int k : input_keys)
    {
        assert(generational_metrics_map[metric_name].find(k) != generational_metrics_map[metric_name].end());
        metric_vec.push_back(getGenerationalMetricValue<T>(metric_name, k));
    }
    return( metric_vec );
}
template vec_t<int> output_statistics_t::getGenerationalMetricValuesForKeys( std::string metric_name, vec_t<int> input_keys );
template vec_t<double> output_statistics_t::getGenerationalMetricValuesForKeys( std::string metric_name, vec_t<int> input_keys );
template vec_t<std::string> output_statistics_t::getGenerationalMetricValuesForKeys( std::string metric_name, vec_t<int> input_keys );

vec_t<std::string> output_statistics_t::getGenerationalMetricNames()
{
    vec_t<std::string> metrics;
    for( auto const &item : generational_metrics_map )
    {
        metrics.push_back(item.first);
    }
    return metrics;
}

void output_statistics_t::addFinalMetricValue( std::string metric_name, int value) { addFinalMetricValueGeneric<int>(metric_name, value); }
void output_statistics_t::addFinalMetricValue( std::string metric_name, double value ) { addFinalMetricValueGeneric<double>(metric_name, value); }
void output_statistics_t::addFinalMetricValue( std::string metric_name, std::string value ) { addFinalMetricValueGeneric<std::string>(metric_name, value); }
template<class T>
void output_statistics_t::addFinalMetricValueGeneric( std::string metric_name, T value )
{
    metric_t metric = value;
    final_metrics_map[metric_name] = metric;
}

template<class T>
T output_statistics_t::getFinalMetricValue( std::string metric_name )
{
    assert( final_metrics_map.find(metric_name) != final_metrics_map.end() );
    return( std::get<T>(final_metrics_map[metric_name]) );   
}
template int output_statistics_t::getFinalMetricValue( std::string metric_name );
template double output_statistics_t::getFinalMetricValue( std::string metric_name );
template std::string output_statistics_t::getFinalMetricValue( std::string metric_name );

vec_t<std::string> output_statistics_t::getFinalMetricNames()
{
    vec_t<std::string> metrics;
    for( auto const &item : final_metrics_map )
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
        for (auto const &v : generational_metrics_map)
        {
            sprintf( string, "%20s ", v.first.c_str());
            fputs(string, file);
        }
        sprintf(string, "\n");
        fputs(string, file);

        for( int k : generational_keys )
        {
            for (auto const &item : generational_metrics_map)
            {
                auto vals = item.second;
                std::string strval = std::visit([](auto&& arg){
                    if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, std::string>) {
                        return arg;
                    } else {
                        return std::to_string(arg);
                    }}, vals[k]);
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
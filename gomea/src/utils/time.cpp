#include "gomea/src/utils/time.hpp"

namespace gomea{
namespace utils{

time_t start_time;
std::unordered_map<std::string,double> timers;

void initStartTime()
{
	start_time = getTimestamp();
}

void printTimestamp( time_t timestamp )
{
	//std::chrono::time_point<std::chrono::system_clock, std::chrono::microseconds> time_point;
	//std::time_t tc = std::chrono::high_resolution_clock::to_time_t(timestamp);
	//std::cout << "time: " << timestamp << std::endl;
	//std::cout << "time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(timestamp.time_since_epoch()).count();
	std::cout << "time: " << timestamp.time_since_epoch().count() << std::endl;
}

time_t getTimestamp()
{
	//std::chrono::time_point<std::chrono::system_clock,std::chrono::microseconds> time_point;
	return std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::system_clock::now());
}

void addToTimer( std::string label, time_t startTimestamp )
{
	if( timers.count(label) == 0 )
		timers[label] = 0.0;
	timers[label] += getElapsedTimeSeconds(startTimestamp);
}

void clearTimers()
{
	timers.clear();
}

double getTimer( std::string label )
{
	if( timers.count(label) == 0 )
		return 0.0;
	return timers[label];
}

long long getElapsedTimeMicroseconds(time_t startTimestamp)
{
    auto timestamp_now = getTimestamp();
	auto diff = std::chrono::duration_cast<std::chrono::microseconds>(timestamp_now - startTimestamp);
	return( diff.count() );
}

double getElapsedTimeMilliseconds(time_t startTimestamp)
{
	return( getElapsedTimeMicroseconds(startTimestamp) / 1e3 );
}

double getElapsedTimeSeconds(time_t startTimestamp)
{
	return( getElapsedTimeMicroseconds(startTimestamp) / 1e6 );
}
        
double getElapsedTimeSinceStartSeconds()
{
	return getElapsedTimeSeconds(start_time);
}

double getElapsedTimeSinceStartMilliseconds()
{
	return getElapsedTimeMicroseconds(start_time)/1000.0;
}

}}


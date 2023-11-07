#pragma once

#include <iostream>
#include <chrono>
#include <cstddef>
#include <string>
#include <unordered_map>

namespace gomea{
    typedef std::chrono::time_point<std::chrono::system_clock,std::chrono::microseconds> time_t;
    //typedef std::chrono::high_resolution_clock::time_point time_t;
    namespace utils{
        extern time_t start_time;

        void initStartTime();
        void addToTimer( std::string label, time_t startTimestamp );
        void clearTimers();
        double getTimer(std::string label);
        time_t getTimestamp();
        void printTimestamp( time_t timestamp );
        long long getElapsedTimeMicroseconds(time_t startTimestamp);
        double getElapsedTimeMilliseconds(time_t startTimestamp);
        double getElapsedTimeSeconds(time_t startTimestamp);
        double getElapsedTimeSinceStartMilliseconds();
        double getElapsedTimeSinceStartSeconds();
        double getElapsedTimeSinceStartSeconds(time_t startTimestamp);
    }
}

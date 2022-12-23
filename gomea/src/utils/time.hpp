#pragma once

#include <chrono>
#include <cstddef>

namespace gomea{
    typedef std::chrono::high_resolution_clock::time_point time_t;
    namespace utils{
        extern time_t start_time;

        void initStartTime();
        time_t getTimestamp();
        long long getElapsedTimeMilliseconds(time_t startTimestamp);
        double getElapsedTimeSeconds(time_t startTimestamp);
        long long getElapsedTimeSinceStartMilliseconds();
        double getElapsedTimeSinceStartSeconds();
    }
}

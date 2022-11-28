#pragma once

#include <chrono>
#include <cstddef>

namespace gomea{
    typedef std::chrono::high_resolution_clock::time_point time_t;
    namespace utils{
        time_t getTimestamp();
        long long getElapsedTimeMilliseconds(time_t startTimestamp);
        double getElapsedTimeSeconds(time_t startTimestamp);
    }
}

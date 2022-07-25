#include "time.hpp"

long long getTimestamp()
{
    struct timeval tv;
    long long result;

    gettimeofday (&tv, NULL);
    result = (tv.tv_sec * 1000) + (tv.tv_usec / 1000);

    return  result;
}

long long getTime(long long startTimestamp)
{
    long long timestamp_now, difference;

    timestamp_now = getTimestamp();

    difference = timestamp_now-startTimestamp;

    return difference;
}

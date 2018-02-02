#ifndef TIME_KEEPER_H
#define TIME_KEEPER_H

#include <time.h>

class TimeKeeper
{
public:
    TimeKeeper(){initial_time = clock();}
    void StartMeasurement(){initial_time = clock();}
    float GiveMeasurement(){return float(clock()-initial_time)/CLOCKS_PER_SEC;}
private:
    clock_t initial_time;
};
#endif

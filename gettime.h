#ifndef _BENCH_GETTIME_INCLUDED
#define _BENCH_GETTIME_INCLUDED

#include <stdlib.h>
#include <sys/time.h>
#include <iomanip>
#include <iostream>

struct timer {
    double totalTime;
    double lastTime;
    double totalWeight;
    bool on;
    struct timezone tzp;
    timer() {
        struct timezone tz = {0, 0};
        totalTime=0.0; 
        totalWeight=0.0;
        on=0; tzp = tz;}
        double getTime() {
        timeval now;
        gettimeofday(&now, &tzp);
        return ((double) now.tv_sec) + ((double) now.tv_usec)/1000000.;
    }
    void start () {
        on = 1;
        lastTime = getTime();
    } 
    double stop () {
        on = 0;
        double d = (getTime()-lastTime);
        totalTime += d;
        return d;
    } 
    double stop (double weight) {
        on = 0;
        totalWeight += weight;
        double d = (getTime()-lastTime);
        totalTime += weight*d;
        return d;
    } 

    double total() {
        if (on) 
          return totalTime + getTime() - lastTime;
        else 
          return totalTime;
    }

    double next() {
        if (!on) 
          return 0.0;
        double t = getTime();
        double td = t - lastTime;
        totalTime += td;
        lastTime = t;
        return td;
    }

    void reportT(double time) {
        std::cout << "Algorithm Running Time: " << std::setprecision(3) << time <<  std::endl;;
    }
};

static timer timer_var;
#define startTime() timer_var.start();
#define nextTime(_string) timer_var.reportNext(_string);
#define nextTimeN() timer_var.reportT(timer_var.next());

#endif 


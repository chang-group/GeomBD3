#ifndef _Timer_h_
#define _Timer_h_

class Timer {
  public:
    timeval startTime;
    double duration;

    void start(){
      gettimeofday(&startTime, NULL);
    }

    void stop(){
      timeval endTime;
      long seconds, useconds;

      gettimeofday(&endTime, NULL);

      seconds  = endTime.tv_sec  - startTime.tv_sec;
      useconds = endTime.tv_usec - startTime.tv_usec;

      duration = seconds + useconds/1000000.0;
    }

    void print(ostream *stream){
      double seconds = duration;

      // calculate human time
      int minutes = seconds / 60;
      seconds -= (minutes * 60.);
      int hours = minutes / 60;
      minutes -= (hours * 60);

      // output timing information
      *stream << "time: ";
      if(hours > 0.) {
        *stream << hours << " hours ";
      }
      if(minutes > 0.) {
        *stream << minutes << " minutes ";
      }
      *stream << seconds << " seconds elapsed." << endl;
    }

    void log_current(ostream *stream) {
      stop(); //doesn't actually stop anything, just grab current time

      double seconds = duration;

      // calculate human time
      int minutes = seconds / 60;
      seconds -= (minutes * 60.);
      int hours = minutes / 60;
      minutes -= (hours * 60);

      // output timing information
      if(hours > 0.) {
        *stream << hours << ":";
      }
      if(minutes > 0.) {
        *stream << minutes << ":";
      }
      *stream << seconds;
    }
};

#endif

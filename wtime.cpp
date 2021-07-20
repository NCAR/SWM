#include <sys/time.h>

double wtime()
{
   /* Use a generic timer */
   static int sec = -1;
   struct timeval tv;
   struct timezone tz;
   gettimeofday(&tv, &tz);
   if (sec < 0) sec = tv.tv_sec;
   return (tv.tv_sec - sec) + 1.0e-6*tv.tv_usec;
}

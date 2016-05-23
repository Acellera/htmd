/* Modifications by Acellera Ltd */

#ifdef TIMEBOMB
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

//Default is valid for 60 days
#ifndef LAPSE
  #define LAPSE (86400*60)
//39312000
#endif
//Warning of expiration 15 days before
#define WARNTIME 1296000 

void timebomb() {
  unsigned int cdate = CDATE;
  unsigned int lapse = LAPSE;
  time_t curr_time = time(NULL);
//  printf("Unix compilation time %d, time %d\n", cdate, curr_time);  
  if (labs(curr_time - cdate) > lapse - WARNTIME) {//abs in case they set it to very old time           
    fprintf ( stderr, "# WARNING. This build is about to expire. Please write to info@acellera.com for more information.\n" );
    fprintf ( stdout, "# WARNING. This build is about to expire. Please write to info@acellera.com for more information.\n" );
  }
  if (labs(curr_time - cdate) > lapse) {//abs in case they set it to very old time
    fprintf ( stderr, "# WARNING. This build has expired. Please write to info@acellera.com for more information.\n" );
    fprintf ( stdout, "# WARNING. This build has expired. Please write to info@acellera.com for more information.\n" );
    exit(1);
  }
}
#endif


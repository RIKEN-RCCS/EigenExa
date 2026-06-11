#include <stdio.h>
#include <stdlib.h>

#if __IBM_REGISTER_VARS
#	define	GET_DELTA	get_delta
#else
#	define	GET_DELTA	get_delta_
#endif

extern "C"{
void GET_DELTA(void* a1, void* a2, int* i);
}

void
GET_DELTA (void *a1, void *a2, int *i)
{
  long b1 = (long) a1;
  long b2 = (long) a2;
  long delta;

  delta = b1 - b2;
  *i = (int) delta;
}

#if __IBM_REGISTER_VARS
void
flush (int *unit)
{
  fflush (stdout);
  flush_ (unit);
  sleep (0);
}
#endif


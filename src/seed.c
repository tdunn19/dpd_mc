#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

main()
{
	srand(time(NULL) + getpid());
	int aNumber = rand();
	printf("%d\t\tiseed\n", aNumber);
}


/* file: timestamp.c */
/* Print current time. */

#include <stdio.h>
#include <time.h>
#include "timestamp.h"

void timestamp()
{
	char time_buffer[TIME_STAMP_SIZE];
	const struct tm *tm;
	time_t now;

	now = time(NULL);
	tm = localtime(&now);

	strftime(time_buffer, TIME_STAMP_SIZE, "%d %B %Y %I:%M:%S %p", tm);

	fprintf(stdout, "%s\n", time_buffer);

	return;
}

// File: errorlib.c

/*
 * Subroutines to handle errors and warnings.
 * This provides a straightforward error flag implementation.
 *
 * Subroutines:
 *  error_flag
 *  error_message
 */

#include <stdio.h>
#include <stdlib.h>
#include "errorlib.h"

/*
 * error_flag: Alert user to function flagging an error.
 */
void error_flag(int rank, int error_value, char *function_name)
{
        fprintf(stderr, "*** Error: rank = %d, ", rank);
        fprintf(stderr, "%s: %d\n", function_name, error_value);
        fflush(stderr);
}

/*
 * error_message: Alert user to function flagging an error and print
 * error message.
 */
void error_message(int rank, char *emessage, char *function_name)
{
        fprintf(stderr, "*** Error: rank = %d, ", rank);
        fprintf(stderr, "%s: %s\n", function_name, emessage);
        fflush(stderr);
}

// File: errorlib.c

/*
 * Subroutines to handle errors and warnings.
 * This provides a straightforward error flag implementation.
 *
 * Subroutines:
 *  error_flag
 */

#include <stdio.h>
#include <stdlib.h>
#include "errorlib.h"

/*
 * error_flag: Alert user to function flagging an error.
 */
void error_flag(int error_value, char *function_name)
{
    fprintf(stderr, "*** Error: ");
    fprintf(stderr, "%s: %d\n", function_name, error_value);
    return;
}

/*
 * error_message: Alert user to function flagging and error, printing
 * a warning message. 
 */
void error_message(char *function_name, char *message)
{
        fprintf(stderr, "*** Error: ");
        fprintf(stderr, "%s: %s\n", function_name, message);
        return;
}

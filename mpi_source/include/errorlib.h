// File: errorlib.h

#ifndef errorlib_h
#define errorlib_h

/*
 * error_flag: Alert user to function flagging an error.
 */
void error_flag(
        int rank,           /* Process rank */
        int error_value,    /* Value of error flag */
        char *function_name /* Name of function flagging error */
    );

/*
 * error_message: Alert user to function flaggin an error and print
 * error message.
 */
void error_message(
        int rank,           /* Process rank */
        char *emessage,     /* Error message to print */
        char *function_name /* Function name */
        );

#endif

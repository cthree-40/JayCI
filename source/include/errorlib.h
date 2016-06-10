// File: errorlib.h

#ifndef errorlib_h
#define errorlib_h

/*
 * error_flag: Alert user to function flagging an error.
 */
void error_flag(
    int error_value, /* Value of error flag */
    char *function_name /* Name of function flagging error */
    );

#endif

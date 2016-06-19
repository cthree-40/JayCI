// File: cmdline.h

/*
 * Process command line arguments.
 */

#ifndef cmdline_h
#define cmdline_h

/*
 * processcmdlineargs: processes command line arguments. Output of function
 * is error flag. Returns amoutn of memory to be allocated in bytes.
 */
int processcmdlineargs(
    int argc, /* Number of arguments */
    char *argv[], /* Value of arguments */
    long long int *memory /* Memory to be allocated. */
    );

#endif

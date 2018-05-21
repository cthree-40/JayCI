// File: binary_strings.h
/*
 * Utilities for binary string representation of electron orbital strings.
 */
#ifndef binary_strings_h
#define binary_strings_h

/*
 * struct eostring: electron occupation string
 */
struct eostring {
        int doccx;
        int actvx;
        int index;
        struct occstr string;
};

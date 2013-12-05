#!/usr/bin/perl

#####################################################################
# This script generates the following input files:
#     det.list
#     pcross.ref
#     qcross.ref
#     pstep.list
#     qstep.list
#     pstring.list
#     qstring.list
#####################################################################

# Only to be used for FCI calculations

#####################################################################

use strict;
use warnings;

print "Beginning...";

my $dtrmnts = 4008004;
my $astrngs = 2002;
my $bstrngs = 2002;

# Open the files pcross.ref and qcross.ref
open(<PXREF>, ">pcross.ref");
my $i=0;
my $j=0;
for ( $i = 1; $i <= $astrings; $i++ ){ 
  for ( $j = 1; $j <= $bstrings; $j++ ){
    print PXREF " 

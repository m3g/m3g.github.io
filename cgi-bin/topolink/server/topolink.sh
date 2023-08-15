#!/bin/bash

inputfile=$1
pdbfile=$2
output=$3
#email=$4

dir=$(dirname $output)

./runback.sh $inputfile $pdbfile $output >& /dev/null &
#./runback.sh $inputfile $pdbfile $output $email >& /dev/null &

exit


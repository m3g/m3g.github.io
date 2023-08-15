#!/bin/bash
echo "Executing ... please wait ... "
gfortran -o ./agriprec2 ../../agriprec/agriprec2.for
./agriprec2 < ./files/data.dat > ./files/result.dat

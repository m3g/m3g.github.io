#!/bin/bash

inputfile=$1
pdbfile=$2
output=$3
#email=$4

dir=$(dirname $output)

echo $inputfile > $dir/runlog.txt
echo $pdbfile >> $dir/runlog.txt
echo $output >> $dir/runlog.txt
#echo $email >> $dir/runlog.txt

tmpdir=`basename $dir`

printlinks=`grep printlinks $inputfile`

pdbfile_name=$(basename $pdbfile)

./wait_for.tcl www-data topolink-web 2 2

./topolink-web $inputfile $pdbfile >& $output 

./write_index_final.py $dir $(basename $pdbfile) $(basename $inputfile)

if [[ "$printlinks" == *"yes"* ]] ; then
  zip -j $dir/links.zip $dir/links/* >& $dir/zip.log
fi

#echo " 
#
#Please find the output your TopoLink run at: 
#
#http://leandro.iqm.unicamp.br/topolink/serverfiles/$tmpdir/index.shtml
#
#This link will be available for about 3 days.
#
#If you have any question, please contact Leandro Martinez, at leandro@iqm.unicamp.br
#
#If TopoLink was useful to you, please remember to reference it by citing:
#
#A. Ferrari, F. C. Gozzo, L. Martinez, 
#Evaluation of structural models using chemical crosslinking distance constraints.
#2018.
#
#" \
#| mutt -F ./muttrc -s "TopoLink result for: $pdbfile_name" -- $email >& $dir/mail.log



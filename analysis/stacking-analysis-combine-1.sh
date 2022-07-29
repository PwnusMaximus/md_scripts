#!/bin/bash

#This script is intended to be run in the "upper" stacking
#folder. It will output a final "FINAL-FINAL.dat" file in the 
#Directory where this is run from 
#
#when done, take the FINAL-FINAL.dat file, name it for the system and move
#use 'stacking-final-sort-merge.sh' to combine the files. 


for i in stacking-* 
	do
	cd $i
	echo "Moved to "$i"" 
	echo "working on "$i""
	sed -i '/Undefined/d' ALL-stacking*.dat 
	sed -i 's/stacked//g' ALL-stacking*.dat
	sed -i 's/==>//g' ALL-stacking*.dat
	sed -i 's/<==//g' ALL-stacking*.dat
	sed -i 's/stacking-count-//g' ALL-stacking*.dat
	sed -i 's/.dat//g' ALL-stacking*.dat
	grep --no-filename '[0-9][A-Z]' ALL-stacking*.dat > NAMES.dat
#	cat NAMES.dat 
	sed -i "s/^/${i}/" NAMES.dat
	sed -i 's/stacking-//g' NAMES.dat 
	grep --no-filename '[0-9][0-9][0-9][0-9]' ALL-stacking*.dat > DATA.dat
#	cat DATA.dat
	paste -d , NAMES.dat DATA.dat > FINAL.dat 
#	cat FINAL.dat
	cd ../
	echo "Done with "$i"" 
done 

#search collect and combine all files just edited
find . -type f -name 'FINAL.dat' -exec cat {} >> FINAL-FINAL.dat \;
sed -i 's/ /-/' FINAL-FINAL.dat #make the first column one single string
#cat  FINAL-FINAL.dat 

echo "done" 

#!/bin/bash

# This script looks for all aromatic stacking interactions with a reference residue. Interaction present when the aromatic center of masses are < 5 Ang and have a interplanar angle < 30. These calues can be changed with the appropriate variables in the header.
# usage: replace the lines below (7-9) specifying the name of the prmtop, trajectory and the residue you want to search interactions for. Then run the script ('bash run_find_stacking.sh')
#works on Residues A,DA,G,DG,C,DG,T,DT,U,PHE,TYR,Phe,Tyr,HIS,His,HIP,Hip,HID,Hid,TRP and Trp. Any other modified residues have to be added to the script (template lines 60-65) or done by hand. 
# EDIT HERE and the sequence section below
PRMTOP="../../stripped.prmtop"
TRAJFILE="trajin ../../ALL-CONCAT.nc 1 last 40" 
REF_NUM="2"
REF_RES=":2@N1,C2,N3,C4,C5,C6"; # Change this for the atoms you are wanting to use as reference
AROMATIC_RES=":I2G,I2A,I2C,I2U,I1G,I1A,I1C,I1U,DT5,DT3,G5,G3,A5,A3,U5,U3,C5,C3,A,DA,G,DG,DC,C,T,DT,U,PHE,TYR,Phe,Tyr,HIS,His,HIP,Hip,HID,Hid,TRP,Trp" #Theses are the residues that it is built to search.
cutoff_distance1="5" # cutoff of when any atoms are interacting
MINDIST="0" #minimum distance cutoff for stacking
MAXDIST="5" #maximum distance cutoff for stacking
MINANG1="0" #stacking definition 1, minimum angle
MAXANG1="30" #stacking definition 1, maximum angle
MINANG2="150" #'flipped' stacking definition 2, minimum angle
MAXANG2="180" #'flipped' stacking definition 2, maximum angle


# create file to get frames file
echo "getting residues in range..."
echo "parm $PRMTOP
$TRAJFILE
mask ($REF_RES<@$cutoff_distance1)&($AROMATIC_RES)&!($(echo "$REF_RES" | grep -o ":[0-9]\+")) maskout wedge_frames.dat
run
quit" > get_distance_frames.in


# use mask to get frames based on the distance cutoffs
cpptraj -i get_distance_frames.in  > /dev/null # output is not printed to terminal

#get resn and resid
awk '{ print $4$5 }' wedge_frames.dat | sed '1d' | sort | uniq > stacking_residues.log


#starts file to calculate the vectors and distances
echo "determining if stacked..."
echo "parm $PRMTOP
$TRAJFILE
vector VECT-ref-ring corrplane $REF_RES" > stacking.in

#for each residue that is within the cutoff distance; get vector and distance lines
for res in $(cat stacking_residues.log) #get each residue within the distance cutoff
do
	if [[ $( echo "$res"  | grep -o "^[0-9].*" | sed '0,/[0-9]*/{s/^[0-9]*//g}' ) == "I2U" ]] || [[ $( echo "$res"  | grep -o "^[0-9].*" | sed '0,/[0-9]*/{s/^[0-9]*//g}' ) == "I2C" ]] || [[ $( echo "$res"  | grep -o "^[0-9].*" | sed '0,/[0-9]*/{s/^[0-9]*//g}' ) == "I1U" ]] || [[ $( echo "$res"  | grep -o "^[0-9].*" | sed '0,/[0-9]*/{s/^[0-9]*//g}' ) == "I1C" ]] || [[ $( echo "$res"  | grep -o "^[0-9].*" | sed '0,/[0-9]*/{s/^[0-9]*//g}' ) == "DT3" ]] || [[ $( echo "$res"  | grep -o "^[0-9].*" | sed '0,/[0-9]*/{s/^[0-9]*//g}' ) == "DT5" ]]  || [[ $( echo "$res"  | grep -o "^[0-9].*" | sed '0,/[0-9]*/{s/^[0-9]*//g}' ) == "C3" ]] || [[ $( echo "$res"  | grep -o "^[0-9].*" | sed '0,/[0-9]*/{s/^[0-9]*//g}' ) == "C5" ]] || [[ $( echo "$res"  | grep -o "^[0-9].*" | sed '0,/[0-9]*/{s/^[0-9]*//g}' ) == "U3" ]] || [[ $( echo "$res"  | grep -o "^[0-9].*" | sed '0,/[0-9]*/{s/^[0-9]*//g}' ) == "U5" ]] || [[ $( echo "$res"  | grep -o "^[0-9].*" "[A-Z]\+" ) == "DC" ]] || [[ $( echo "$res"  | grep -o "^[0-9].*" "[A-Z]\+" ) == "C" ]] || [[ $( echo "$res"  | grep -o "^[0-9].*" "[A-Z]\+" ) == "DT" ]] || [[ $( echo "$res"  | grep -o "^[0-9].*" "[A-Z]\+" ) == "T" ]] || [[ $( echo "$res"  | grep -o "^[0-9].*" "[A-Z]\+" ) == "U" ]]; #checks for Pyrimidines
	then
		echo "vector VECT-$res-ring corrplane :$(echo "$res" | grep -o "^\s*[0-9]+")@N1,C2,N3,C4,C5,C6" >> stacking.in
		echo "distance DIST-$res :$(echo "$res" | grep -o "[0-9]\+")@N1,C2,N3,C4,C5,C6 $REF_RES out stacking-$res.out" >> stacking.in
				
	elif [[ $( echo "$res"  | grep -o "^[0-9].*" | sed '0,/[0-9]*/{s/^[0-9]*//g}' ) == "I2A" ]] || [[ $( echo "$res"  | grep -o "^[0-9].*" | sed '0,/[0-9]*/{s/^[0-9]*//g}' ) == "I2G" ]] || [[ $( echo "$res"  | grep -o "^[0-9].*" | sed '0,/[0-9]*/{s/^[0-9]*//g}' ) == "I1A" ]] || [[ $( echo "$res"  | grep -o "^[0-9].*" | sed '0,/[0-9]*/{s/^[0-9]*//g}' ) == "I1G" ]] || [[ $( echo "$res"  | grep -o "^[0-9].*" | sed '0,/[0-9]*/{s/^[0-9]*//g}' ) == "A3" ]] || [[ $( echo "$res"  | grep -o "^[0-9].*" | sed '0,/[0-9]*/{s/^[0-9]*//g}' ) == "A5" ]] || [[ $( echo "$res"  | grep -o "^[0-9].*" | sed '0,/[0-9]*/{s/^[0-9]*//g}' ) == "G3" ]] || [[ $( echo "$res"  | grep -o "^[0-9].*" | sed '0,/[0-9]*/{s/^[0-9]*//g}' ) == "G5" ]] || [[ $( echo "$res"  | grep -o "^[0-9].*" "[A-Z]\+" ) == "DA" ]] || [[ $( echo "$res"  | grep -o "^[0-9].*" "[A-Z]\+" ) == "A" ]] || [[ $( echo "$res"  | grep -o "^[0-9].*" "[A-Z]\+" ) == "DG" ]] || [[ $( echo "$res"  | grep -o "^[0-9].*" "[A-Z]\+" ) == "G" ]]; # checks for purines
	then
		echo "vector VECT-$res-ring corrplane :$(echo "$res" | grep -o "^\s*[0-9]+")@N1,C2,N3,C4,C5,C6,N7,C8,N9" >> stacking.in
		echo "distance DIST-$res :$(echo "$res" | grep -o "[0-9]\+")@N1,C2,N3,C4,C5,C6,N7,C8,N9 $REF_RES out stacking-$res.out" >> stacking.in
	elif [[  $( echo "$res"  | grep -o "[A-Z]\+" ) == "PHE" ]] || [[ $( echo "$res"  | grep -o "[A-Z]\+" ) == "Phe" ]] || [[ $( echo "$res"  | grep -o "[A-Z]\+" ) == "Tyr" ]] || [[ $( echo "$res"  | grep -o "[A-Z]\+" ) == "TYR" ]]; # checks for Phe and Tyr
	then
		echo "vector VECT-$res-ring corrplane :$(echo "$res" | grep -o "^[0-9]\+")@CG,CD1,CD2,CE1,CE2,CZ" >> stacking.in
		echo "distance DIST-$res :$(echo "$res" | grep -o "[0-9]\+")@CG,CD1,CD2,CE1,CE2,CZ $REF_RES out stacking-$res.out" >> stacking.in
	elif [[  $( echo "$res"  | grep -o "[A-Z]\+" ) == "HIS" ]] || [[ $( echo "$res"  | grep -o "[A-Z]\+" ) == "His" ]] || [[  $( echo "$res"  | grep -o "[A-Z]\+" ) == "HIP" ]] || [[ $( echo "$res"  | grep -o "[A-Z]\+" ) == "Hip" ]] || [[  $( echo "$res"  | grep -o "[A-Z]\+" ) == "HID" ]] || [[ $( echo "$res"  | grep -o "[A-Z]\+" ) == "Hid" ]]
	then
		echo "vector VECT-$res-ring corrplane :$(echo "$res" | grep -o "[0-9]\+")@CG,ND1,CD2,CE1,NE2" >> stacking.in
		echo "distance DIST-$res :$(echo "$res" | grep -o "[0-9]\+")@CG,ND1,CD2,CE1,NE2 $REF_RES out stacking-$res.out" >> stacking.in
	elif [[  $( echo "$res"  | grep -o "[A-Z]\+" ) == "TRP" ]] || [[ $( echo "$res"  | grep -o "[A-Z]\+" ) == "Trp" ]]
	then
		echo "vector VECT-$res-ring corrplane :$(echo "$res" | grep -o "[0-9]\+")@CG,CD1,CD2,NE1,CE2,CE3,CZ2,CZ3,CH2" >> stacking.in
		echo "distance DIST-$res :$(echo "$res" | grep -o "[0-9]\+")@CG,CD1,CD2,NE1,CE2,CE3,CZ2,CZ3,CH2 $REF_RES out stacking-$res.out" >> stacking.in
	elif [[  $( echo "$res"  | grep -o "[A-Z]\+" ) == "ARG" ]] || [[ $( echo "$res"  | grep -o "[A-Z]\+" ) == "Arg" ]]
	then
		echo "vector VECT-$res-ring corrplane :$(echo "$res" | grep -o "[0-9]\+")@NE,CZ,NH1,NH2" >> stacking.in
		echo "distance DIST-$res :$(echo "$res" | grep -o "[0-9]\+")@NE,CZ,NH1,NH2 $REF_RES out stacking-$res.out" >> stacking.in
#	elif [[  $( echo "$res"  | grep -o "[A-Z]\+" ) == "THREE_LETTER_CODE" ]] || [[ $( echo "$res"  | grep -o "[A-Z]\+" ) == "ALTERNATIVE_THREE_LETTER_CODE" ]]
#	then
#		echo "vector VECT$res-ring corrplane :$(echo "$res" | grep -o "[0-9]\+")@CG,CD1,CD2,NE1,CE2,CE3,CZ2,CZ3,CH2" >> stacking.in #Change atom names to be the aromatic portion of the resid 
#		echo "distance DIST-$res :$(echo "$res" | grep -o "[0-9]\+")@CG,CD1,CD2,NE1,CE2,CE3,CZ2,CZ3,CH2 $REF_RES out stacking-$res.out" > #Change atom names here to the same ones 
	fi
	echo "vectormath vec1 VECT-ref-ring vec2 VECT-$res-ring out stacking-$res.out dotangle name $res-DOTANGLE" >> stacking.in
#Determines if the distance and angle qualify as 'stacked' for any given interaction and determines occupancy
	echo "calcstate state stacked,DIST-$res,$MINDIST,$MAXDIST,$res-DOTANGLE,$MINANG1,$MAXANG1 state stacked,DIST-$res,$MINDIST,$MAXDIST,$res-DOTANGLE,$MINANG2,$MAXANG2 out stacking-$res.dat countout stacking-count-$res.dat" >> stacking.in
done

#finish the vector and distance command
echo "run
quit" >> stacking.in

#execute created cpptraj input file
cpptraj -i stacking.in > /dev/null

#rearrange and output the final data. 
echo "final output file..."
gawk -i inplace '{print $3, $4}' stacking-count-*.dat
tail -n +1 stacking-count-*.dat > ALL-stacking-$REF_NUM.dat
sed -i '/^$/d' ALL-stacking-$REF_NUM.dat
sed -i '/State_/d' ALL-stacking-$REF_NUM.dat

echo "$REF_RES ...Done"


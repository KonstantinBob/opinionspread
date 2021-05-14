#!/usr/bin/env bash 

# constant
usage="Resultscript exeName logName flagMake flagMode \n 	make flags: -s speed -d debug \n 	mode flags:  -g dot pictures "

# clear screen
clear

# check arguments
if [ "$#" -ne 4 ];
then
	echo -e $usage
	exit 1;
fi

# save arguments
exeName=$1
logName=$2
makeFlag=$3
modeFlag=$4

# do profile
if [ $makeFlag == "-d" ];
then
	gprof $exeName gmon.out > analysis_$logName.txt
fi

# make the graph plots
if [ $modeFlag == "-g" ];
then
	#format selection
	#format="pdf"
	format="png"
	
	for file in *.dot
	do
		neato -T$format $file -o $(basename "$file" | cut -d. -f1).$format 
	done
fi

# make the  plots
if [ $modeFlag == "-p" ];
then
	for file in f_*
	do
		Histoplot.py $file 
	done
fi

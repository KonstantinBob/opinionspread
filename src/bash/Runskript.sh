#!/usr/bin/env bash 

# constant
usage="Runscript executableName outputFilename flagMake flagRun pVal qVal nodes numberNet\n 	make flags: -s speed -d debug \n 	run flags:  -o output -h nohup"

# clear screen
clear

#check arguments
if [ "$#" -ne 8 ];
then
	echo -e $usage
	exit 1;
fi

#save arguments
exeName=$1
logName=$2
makeFlag=$3
runFlag=$4
pVal=$5
qVal=$6
nodes=$7
numberNet=$8
#target directory for data
target=$(pwd)

# make shure that you are in the right directory
cd ~/Thesis/Code

# make the executable
if [ $makeFlag == "-s" ];
then
	make -B speed;
elif [ $makeFlag == "-d" ];
then
	make debug;
else
	echo "Unknown argument" $3
	exit 1;
fi

if [ $? -eq 0 ];
then
	# move it to data directory
	mv Echse $target/$exeName
			
	# change directory to there
	cd $target
	
	#clear again
	clear
	
	# run production
	if [ $runFlag == "-o" ];
	then
		./$exeName ${pVal} ${qVal} ${nodes} ${numberNet} ;
	elif [ $runFlag == "-h" ];
	then
		stdbuf -oL nohup ./$exeName ${pVal} ${qVal} ${nodes} ${numberNet} </dev/null > $logName.log 2>&1 &
	else
		echo -e $usage
		exit 1;
	fi
		
fi

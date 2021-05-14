#!/usr/bin/env bash 

usage="BatchSkript executableName outputFilename flagMake flagRun netType \n 	make flags: -s speed -d debug \n 	run flags:  -o output -h nohup"

#check arguments
if [ "$#" -ne 5 ];
then
	echo -e $usage
	exit 1;
fi

#save arguments
exeName=$1
logName=$2
makeFlag=$3
runFlag=$4		
netType=$5

# dictionary for network names
declare -A networkNames
networkNames=( ["0"]="FF" ["1"]="ER" ["2"]="SK" ["3"]="WS" ["4"]="BA")
	
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
	# clear make output
	clear

	# move it to data directory
	mv Echse $target/$exeName

	# change directory to there
	cd $target

	# execution					
	nice -n 0 parallel  -j 8 "stdbuf -oL nohup ./$exeName 0.95 0.95 {1} ${netType} </dev/null > ${networkNames[${netType}]}_{1}_0.95_0.95_.log 2>&1 " :::  125 160 250 325 500 750 1000 1500 2000 3000 4000 6000 8000
fi

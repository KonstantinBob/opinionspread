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

	#set h 
	h=0.05
	
	# 1 -> beta
	# 2-> size
	
	# h->p
	# beta -> q
					
	nice -n 18 parallel  -j 8 "stdbuf -oL nohup ./$exeName $h {1} {2} $netType </dev/null >${logName}_{2}_${h}_{1}_.log 2>&1 " ::: 0.001 0.003 0.005 0.008 0.01 0.03 0.05 0.08 0.1 0.13 0.18 ::: 125 250 500 1000 2000 4000 8000 16000
	

fi

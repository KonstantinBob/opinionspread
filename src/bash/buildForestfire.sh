#!/usr/bin/env bash 

# make shure that you are in the right directory
cd ~/Thesis/Snap/Snap-3.0

# build the new executeable
make all

if [ $? -eq 0 ];
then
	# clear screen
	clear
	
	# move it to data directory
	cp ~/Thesis/Snap/Snap-3.0/examples/forestfire/forestfire ~/Thesis/Data/ForestFireEffDiaCheck/CheckEffDiaDistr
	
fi
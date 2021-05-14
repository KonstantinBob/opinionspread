#!/usr/bin/env bash 

# check the path variable
for path in $(echo $PATH | cut -d:)
do
	echo $path
done

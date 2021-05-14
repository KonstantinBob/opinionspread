#!/usr/bin/env bash 

# loop over packages
for package in $(cat ~/Thesis/TeX/usedPackages.txt)
do
	apt-cache search $package | cut -f1 -d' ' | egrep "texlive*" >> debPackages
done

# purify the results
cat debPackages | sort | uniq  > debPackagesSorted 

# check if an update is needed 

for package in $(cat ~/Thesis/TeX/debPackagesSorted)
do
	sudo apt-get install $package
done

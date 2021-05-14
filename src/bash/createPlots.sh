#!/usr/bin/env bash 

# constant
usage="createPlots  outputName datafiles"

# check arguments
if [ "$#" -lt 2  ];
then
	echo -e $usage
	exit 1;
fi

# parse input
output=$1

# dictionary for network names
declare -A networkNames
networkNames=( ["ER"]="Erdos" ["WS"]="Watts")

# loop over datafiles
for file  in "${@:2}"
do

	# parse the net, size, p and q
	net=$( echo ${file}  | cut -d_ -f1)
	size=$( echo ${file}  | cut -d_ -f2)
	p=$( echo ${file}  | cut -d_ -f3)
	q=$( echo ${file}  | cut -d_ -f4)

	# output tex code
	echo "\begin{figure}" >>${output}
	echo "	\centering"  >>${output}
	echo "	\includegraphics[width=\textwidth]{${file}}" >>${output}
	echo "	\caption[P]{${networkNames[${net}]} with \$\\no=${size}\$, \$\\p=${p}\$ and \$\\q=${q}\$}" >>${output}
	echo "\end{figure}" >>${output}
done

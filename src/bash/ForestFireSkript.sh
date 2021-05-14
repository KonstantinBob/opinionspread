#!/usr/bin/env bash 

for i in `seq 1 250`;
        do
           ./forestfire -n:1698 -f:0.37 -b:0.32  | tail -n 1 >> effDiaDistrSnap.dat	    	  
        done
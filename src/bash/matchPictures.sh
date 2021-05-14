#!/usr/bin/env bash 

for file in custom*.png;
        do
            # get the number of the file
            filenumber=$(echo $file | sed -e s/[^0-9]//g)
            convert $file graphFF$filenumber.png +append result$filenumber.png
        done
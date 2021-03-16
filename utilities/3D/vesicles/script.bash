#!/bin/bash

files=("ml_getlinept_mex.mexa64" "ml_hitmiss.mexa64" "ml_moments_1.mexa64" "ml_multout.mexa64" 
"ml_texture.mexa64" "ml_Znl.mexa64")

for file in ${files[*]}
do
    find /home/icaoberg.projects/slic/source -name "$file" -exec cp {} . \;
done


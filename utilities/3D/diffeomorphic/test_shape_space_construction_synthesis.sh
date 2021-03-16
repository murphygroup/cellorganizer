#!/bin/bash
command="test_shape_space_construction_synthesis"
#PBS -N $command -l nodes=1:ppn=1
##PBS -N $command -l nodes=1:ppn=12
##PBS -l walltime=10:00:00
#PBS -q model1
##PBS -q pool1
##PBS -l mem=4g
cd $PBS_O_WORKDIR
# unset DISPLAY
/usr/local/bin/matlab -nosplash -nodisplay -nodesktop -r "try; ${command}; catch err; getReport(err, 'extended'), end; exit"


#!/bin/bash

source ./bin/activate

make clean
make html
rsync -ruv source/images/ build/html/_images/ 
deactivate

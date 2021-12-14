#!/bin/bash

python3 -m venv .

source ./bin/activate
pip3 install numpy scipy matplotlib
pip3 install -I sphinx
pip3 install xlrd
pip3 install pandas==0.22.0
pip3 install tabulate
pip3 install sphinx_rtd_theme
pip3 install sphinxcontrib.gist
pip3 install sphinxcontrib-programoutput
pip3 install -e git+git://github.com/icaoberg/sphinxcontrib-pyexec.git@master#egg=sphinxcontrib-pyexec
deactivate

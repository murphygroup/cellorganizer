About CellOrganizer for Matlab
******************************
CellOrganizer for Matlab is the most flexible and powerful of the CellOrganizer suite, since it interfaces directly with Matlab which facilities customization and pre- and post-processing.

Installing CellOrganizer locally
********************************

Requirements
------------
* Matlab 2018a or newer
	* Bioinformatics Toolbox
	* Computer Vision System Toolbox
	* Control System Toolbox
	* Curve Fitting Toolbox
	* Image Processing Toolbox
	* Mapping Toolbox
	* Optimization Toolbox
	* Robust Control Toolbox
	* Signal Processing Toolbox
	* Simulink
	* Simulink Design Optimization
	* Statistics and Machine Learning Toolbox
	* System Identification Toolbox
	* Wavelet Toolbox

Downloading CellOrganizer
-------------------------
To download the latest CellOrganizer for Matlab distribution go to the `download page <http://www.cellorganizer.org/cellorganizer-2-9-0/>`_. Afterwards, extract the contents of the release into a local directory of your preference. 

For example,

.. code-block:: bash

	cd ~/
	wget -nc http://cellorganizer.org/downloads/v2.9/cellorganizer_v2.9.0_and_image_collection.tgz
	tar -xvf cellorganizer_v2.9.0_and_image_collection.tgz
	rm -fv cellorganizer_v2.9.0_and_image_collection.tgz

The commands above will download and extract to disk the contents of CellOrganizer v2.9.0.

Starting CellOrganizer
----------------------
To start using CellOrganizer, `start a Matlab session <https://www.mathworks.com/help/matlab/matlab_env/start-matlab-on-linux-platforms.html>`_ and `change the directory <https://www.mathworks.com/help/matlab/ref/cd.html>`_ to the location of CellOrganizer folder and run `setup.m`.

In the Matlab, type

.. code-block:: matlab

	>> cd( ‘/path/to/folder/cellorganizer’ );
	>> setup();

If you were successful you should see a message like

.. code-block:: bash

	>> setup
	Checking for new stable version... Version is up to date.

You are now ready to use CellOrganizer for Matlab.

Installing CellOrganizer on a HPC
*********************************

Requirements
------------

.. IMPORTANT::
   Remember to contact your HPC cluster managers to know if you have access to a Matlab license that would allow you to run `Matlab <https://www.mathworks.com/products/matlab.html>`_  in your cluster.

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

.. IMPORTANT::
   Sometimes, compute nodes do not access the web. Make sure to download the tarball from the front end or contact your HPC managers for further instructions.

CellOrganizer for Matlab is the most flexible and powerful of the CellOrganizer suite, since it interfaces directly with Matlab which facilities customization and pre- and post-processing.

To download the latest CellOrganizer for Matlab distribution go to the `download page <http://www.cellorganizer.org/cellorganizer-2-9-0/>`_ and download the latest release. Afterwards, extract the contents of the release into a local directory of your preference. 

For example,

.. code-block:: bash

	cd ~/
	wget -nc http://cellorganizer.org/downloads/v2.9/cellorganizer_v2.9.0_and_image_collection.tgz
	tar -xvf cellorganizer_v2.9.0_and_image_collection.tgz
	rm -fv cellorganizer_v2.9.0_and_image_collection.tgz

The commands above will download and extract to disk the contents of CellOrganizer v2.9.0.

Starting CellOrganizer
----------------------

The next instructions assume the HPC cluster you have access to uses `SLURM <https://slurm.schedmd.com/>`_ as its default scheduler. 

.. IMPORTANT::
	Using CellOrganizer for Matlab on all possible popular schedulers is beyond the scope of this document, however feel free to contact us through the mailing list and we will do our best to help you.

Using CellOrganizer for Matlab interactively
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use `srun <https://slurm.schedmd.com/srun.html>`_ or `salloc <https://slurm.schedmd.com/salloc.html>`_ to allocate resources to start Matlab. 

For example 

.. code-block:: bash

	srun -p pool --mem=8Gb --pty /bin/bash

To start using CellOrganizer, `start a Matlab session <https://www.mathworks.com/help/matlab/matlab_env/start-matlab-on-linux-platforms.html>`_ and `change the directory <https://www.mathworks.com/help/matlab/ref/cd.html>`_ to the location of CellOrganizer folder and run `setup.m`. 

In the Matlab, type

.. code-block:: matlab

	cd( ‘/path/to/folder/cellorganizer’ );
	setup();

If you were successful you should see a message like

.. code-block:: bash

	>> setup
	Checking for new stable version... Version is up to date.

You are now ready to use CellOrganizer for Matlab within that allocation or interactive session.

Submitting a job for CellOrganizer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some tasks in CellOrganizer, with special attention on training, require considerable resources. At times, it will be more efficient to submit a job to a scheduler rather than doing it interactively.

For example, to run `demo3D00`, you could create a file called `script.sh` whose contents are

.. code-block:: matlab

	#!/bin/bash
	#
	#$ -j y
	#$ -S /bin/bash
	#$ -cwd

	## the next line selects the partition/queue
	#SBATCH -p pool

	## the next line selects the number of cores
	#SBATCH -n 4

	## the next line selects the memory size
	#SBATCH --mem=8G

	## the next line selects the walltime
	#SBATCH -t 00:30:00

	cd /path/to/cellorganizer/folder
	matlab -nodesktop -nosplash -r "setup(); demo3D00(), exit;"

Then use the command `sbatch <https://slurm.schedmd.com/sbatch.html>`_ to submit it to the scheduler

.. code-block:: bash

	sbatch script.sh

to add the job to the scheduler.

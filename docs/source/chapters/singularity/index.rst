About CellOrganizer for Singularity
***********************************
Singularity enables users to have full control of their environment. Singularity containers can be used to package entire scientific workflows, software and libraries, and even data.

CellOrganizer for Singularity is an image with compiled binaries from CellOrganizer functions

- **img2slml**, the top-level function to train generative models of cells, and
- **slml2img**, the top-level function to generate simulated instances from a trained generative model.
- **slml2info**, the top-level function to generate a report from information extracted from a single generative model.
- **slml2report**, the top-level function to generate a report from comparing generative models.
- **slml2slml**, the top-level function to combine models into a single model file.

Installing CellOrganizer for Singularity
****************************************

About Singularity
=================
Singularity performs operating-system-level virtualization. To learn about Singularity and how to use it, click `here <https://www.sylabs.io/guides/2.6/user-guide/index.html>`_

Setup
-----
The following instructions describe

* How to install Singularity, the virtualization engine that will run the container
* How to download the latest Cellorganizer image
* How to start a shell container from the Singularity image
* How to run some of the demos included in the container

Source Code
===========
The source code to build the Singularity image can be found at `here <https://github.com/icaoberg/singularity-cellorganizer/>`_.

For convenience, the Singularity image can be found at `CellOrganizer <http://www.cellorganizer.org/singularity/>`_.

Installing Singularity
======================
Before downloading the image, you need to install Singularity. Installing Singularity is beyond the scope of this document.

To follow the official instructions to install Singularity, click `here <https://sylabs.io/docs//>`_.

**Download the repository and build the image using Singularity**

Open terminal and enter the commands::

	git clone https://github.com/icaoberg/singularity-cellorganizer
	cd singularity-cellorganizer
	bash ./script.sh .

Download the most recent image using Singularity command line (Recommended)
---------------------------------------------------------------------------

Open terminal and enter the command::

	singularity pull shub://murphygroup/singularity-cellorganizer

Pulling the image from Singularity Hub should produce output similar to

.. raw:: html

	<script id="asciicast-JDe3ygUxbchBTRLjrWVnBlQy3" src="https://asciinema.org/a/JDe3ygUxbchBTRLjrWVnBlQy3.js" async></script>

Once the download is complete, you can confirm the image was downloaded by listing the singularity images of the repository::

    ls murphygroup-singularity-cellorganizer-master-latest.simg 

Demos
=====

There are several demos included within the CellOrganizer software bundle. These demos are intended to illustrate CellOrganizer's functionality, and should be used to familiarize the user with the input/output format of various top-level functions such as **img2slml** and **slml2img**. Certains demos have been deprecated and will be removed in future versions of CellOrganizer.

+----------+------------+-------------+-----------+-------------+
|Demo      | Training   | Synthesis   | Other     | Deprecated  |
+==========+============+=============+===========+=============+
| demo2D00 |            | True        |           |             | 
+----------+------------+-------------+-----------+-------------+
| demo2D01 | True       |             |           |             |
+----------+------------+-------------+-----------+-------------+
| demo2D02 |            | True        |           |             |
+----------+------------+-------------+-----------+-------------+
| demo2D03 | True       |             |           | v2.8.1      |
+----------+------------+-------------+-----------+-------------+
| demo2D04 | True       |             |           |             |
+----------+------------+-------------+-----------+-------------+
| demo2D05 | True       |             |           |             |
+----------+------------+-------------+-----------+-------------+
| demo2D06 |            |  True       |           |             |
+----------+------------+-------------+-----------+-------------+
| demo2D07 |            |  True       |           |             |
+----------+------------+-------------+-----------+-------------+
| demo2D08 |  True      |             |           |             | 
+----------+------------+-------------+-----------+-------------+
| demo2D09 |  True      |             |           |             |
+----------+------------+-------------+-----------+-------------+
| demo3D00 |            | True        |           |             |
+----------+------------+-------------+-----------+-------------+
| demo3D01 |            | True        |           |             |
+----------+------------+-------------+-----------+-------------+
| demo3D04 |            | True        |           |             |
+----------+------------+-------------+-----------+-------------+
| demo3D05 |            | True        |           |             |
+----------+------------+-------------+-----------+-------------+
| demo3D06 |            | True        |           | v2.8.1      |
+----------+------------+-------------+-----------+-------------+
| demo3D07 |            | True        |           |             |
+----------+------------+-------------+-----------+-------------+
| demo3D08 |            | True        |           |             |
+----------+------------+-------------+-----------+-------------+
| demo3D09 |            | True        |           |             |
+----------+------------+-------------+-----------+-------------+
| demo3D10 |            | True        |           |             |
+----------+------------+-------------+-----------+-------------+
| demo3D11 | True       |             |           |             |
+----------+------------+-------------+-----------+-------------+
| demo3D12 | True       |             |           |             |
+----------+------------+-------------+-----------+-------------+
| demo3D15 |            | True        |           |             |
+----------+------------+-------------+-----------+-------------+
| demo3D17 |            | True        |           |             |
+----------+------------+-------------+-----------+-------------+
| demo3D19 | True       |             |   Report  |             | 
+----------+------------+-------------+-----------+-------------+
| demo3D20 | True       |             |   Plot    |             | 
+----------+------------+-------------+-----------+-------------+
| demo3D25 |            | True        |           |             |
+----------+------------+-------------+-----------+-------------+
| demo3D29 | True       |             |   Info    |             | 
+----------+------------+-------------+-----------+-------------+
| demo3D34 |            | True        |           |             |   
+----------+------------+-------------+-----------+-------------+
| demo3D35 | True       |             |  Info     |             | 
+----------+------------+-------------+-----------+-------------+
| demo3D42 | True       |             |           |             |
+----------+------------+-------------+-----------+-------------+
| demo3D44 |            | True        |           |             |
+----------+------------+-------------+-----------+-------------+
| demo3D47 |            |             |  Model    |             |    
+----------+------------+-------------+-----------+-------------+
| demo3D48 | True       |             |           |             |
+----------+------------+-------------+-----------+-------------+
| demo3D50 | True       |             |           |             |
+----------+------------+-------------+-----------+-------------+
| demo3D51 | True       |             |  Plot     |             | 
+----------+------------+-------------+-----------+-------------+
| demo3D52 | True       |             |           |             |
+----------+------------+-------------+-----------+-------------+
| demo3D53 |            | True        |           |             |
+----------+------------+-------------+-----------+-------------+
| demo3D55 |            | True        |  Plot     |             |    
+----------+------------+-------------+-----------+-------------+ 

Running CellOrganizer for Singularity
*************************************

List all applications in the cellorganizer-singularity image
============================================================
To list the CellOrganizer functions included in the image, open Terminal and enter the command::

	singularity apps  singularity-cellorganizer/murphygroup-singularity-cellorganizer-master-latest.simg

This will display these functions

* img2slml
* slml2img
* slml2info
* slml2report
* slml2slml

Run a demo that invokes img2slml
================================
An example of a demo that trains a generative model from a series of `.tif` image files is `demo2D01`. To run this demo, change your current directory to `~/singularity-cellorganizer/demos/2D/demo2D01` by entering::

	cd demos/2D/demo2D01

You should find the shell script **demo2D01.sh**. To run the demo enter the command::

	singularity run -a img2slml ../../../murphygroup-singularity-cellorganizer-master-latest.simg demo2D01.sh

The '-a' flag allows us to specify the function binary that we will use in the script. This demo will save a folder `param` containing .mat files as well as a `.mat` file `lamp2.mat` to the same directory (`~/singularity-cellorganizer/demos/2D/demo2D01`). These `.mat` files contain information characterizing the trained generative model.

Running the demo in the container should produce results similar to

Run a demo that invokes slml2img
================================
An example of a demo that produces simulated images from a trained generative model is `demo2D02`. To run this demo, change your current directory to `/home/singularity-cellorganizer/demos/2D/demo2D02` by entering from your home directory::

	cd demos/2D/demo2D02

You should find the shell script `demo2D02.sh`. To run the demo, enter the command::

	singularity run -a slml2img ../../../murphygroup-singularity-cellorganizer-master-latest.simg demo2D02.sh

This demo will save a folder `img` containing these simulated images to the same directory.

Run custom script that invokes img2slml
=======================================
An example running custom function parameters for img2slml stored within a .txt file. Within this directory (i.e. `/path/to/input.txt`), you can run the command::

	singularity run -a img2slml ~/singularity-cellorganizer/murphygroup-singularity-cellorganizer-master-latest.simg img2slml input.txt

Running CellOrganizer for Docker through Jupyter-Notebook
*********************************************************
Once you have the Jupyter running, now you can start playing around with the demos. We have a set of sample demos that
can help you get started.

Your Jupyter folder structure looks like this.

.. figure:: ../../source/chapters/docker/notebook.png
    :width: 500px


Downloading the sample notebooks, images and models
----------------------------------------------------
The run script had created a folder named mmbios2021 on your desktop. Use this shared folder to transfer files (if needed) between
your system and the docker container.

* Use the Download_files notebook to get the latest images, models and notebooks.

* These are fetched from the web and kept under mmbios2021 folder on your desktop.

.. figure:: ../../source/chapters/docker/download_files.png
    :width: 500px

The demos use various image sets and models located inside mmbios2021/images and mmbios2021/models folder.

The demo notebooks are located under mmbios2021/notebooks and are further divided as:

* demo_notebooks
    2D
    A collection of 2D jupyter notebooks demos

    3D
    A collection of 3D jupyter notebooks demos

* workshop_demos
    A collection of jupyter notebooks used during mmbios2021 workshop

Run a demo that invokes img2slml
--------------------------------
An example of a demo that trains a generative model from a series of `.tif` image files is `demo2D01`. To run a demo, simply click the run button at the top of the notebook.

.. figure:: ../../source/chapters/docker/Run_Button.png
    :width: 500px

This demo will save a folder `param` containing .mat files as well as a `.mat` file `lamp2.mat` to the same directory (`/home/cellorganizer/demos/2D/demo2D01`). These `.mat` files contain information characterizing the trained generative model.

Run a demo that invokes slml2img
--------------------------------
An example of a demo that produces simulated images from a trained generative model is `demo2D02`.

This demo will save a folder `img` containing these simulated images to the same directory.

Export generated data out of the container
------------------------------------------
To export generated data out of the container, click the files in the directory that will be exported and click download.

.. figure:: ../../source/chapters/docker/Download.png
    :width: 500px


Important Docker Configurations
=================================
Some demos are compute heavy and might need more compute resource than which is allowed by default via docker.
To increase RAM memory allocation and number of CPUs. Do the following steps :

* alter these two flags --memory="8g" --cpus=2 to change the memory and cpu allocations.
Mac
===
 Docker > Preference > Resources and change the default allocations.

.. figure:: ../../source/chapters/docker/preferences_docker.png
    :width: 500px

Ubuntu
======
Check the docker documentation to see, how you can use commands to set various parameters.
https://docs.docker.com/config/containers/resource_constraints/


Windows
=======
Check the docker documentation to see, how you can use commands to set various parameters.
https://docs.docker.com/docker-for-windows/

Building Docker image
=====================
* Clone the repo https://github.com/murphygroup/docker-cellorganizer-jupyter-notebook.

* Use the build script to create the docker image ( you will need the matlab binaries of CellOrganizer or it will download the latest from web)

* Once the image builds successfully ( docker images command can be used to check if image is created), you can go ahead to make the container.

* Use the run script to run the container with this docker image.

* The container starts with the url of accessing Jupyter notebook.

Visualizing Results
=====================

* To visualize a result, one can create a notebook cell (the notebook must reside in same location as index.html)

| from IPython.core.display import HTML
| HTML(filename="./index.html")

# Same folder location, ensures the embedded images don't break while rendering.

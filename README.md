CellOrganizer
-------------
I. CellOrganizer takes as input 2D-4D images of cells and creates generative models of many aspects of cell size, shape and internal organization.  It also contains utilities for evaluating, comparing and otherwise using this models.

In addition to the Matlab version, CellOrganizer is also available as compiled executables in a Docker image with Jupyter Notebook support.

The full documentation for CellOrganizer is available at 

https://cellorganizer.readthedocs.io/en/latest/

An extensive tutorial on the Docker/Jupyter version of CellOrganizer can be found at

https://docs.google.com/document/d/10LFLwMmpuQrN8880Sf0ndh7VYeWZDJ7E/

The release notes for v2.10 can be found at

https://docs.google.com/document/d/1D7AAdWc-_xvBu3YCtOfURevpjlCZrZrRAROhJhipCHQ/edit?usp=sharing

Additional information on CellOrganizer is available at

http://www.cellOrganizer.org


II. Using CellOrganizer from the Matlab command line

Enter "startup" at the Matlab command prompt to add the relevant paths.  The main CellOrganizer commands are

> img2slml

which trains a model from a set of images

and

> slml2img

which synthesizes one or more images from a trained model.
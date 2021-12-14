Input Type
-------------
In order to train a model, CellOrganizer requires you to have an image dataset. CellOrganizer provides sample image datasets for you to explore the software with.
These sample image datasets are further described `here <http://murphylab.web.cmu.edu/data/>`_. Once comfortable with the provided datasets, you are encouraged to upload your own image dataset for further
computation given that all images in your dataset will meet certain parameters.


TIFF vs OME.TIFF
-----------------
A TIFF image is a format that retains more pixels and provides a better resolution than a regular JPEG image. However, an OME.TIFF image retains the pixels of the image plus the metadata contained with it in an XML piece within the image.


What is an OME.TIFF format?
----------------------------
OME.TIFF is a bio-format that retains the pixels in multipage TIFF format, and simultaneously contains XML metadata because of the TIFF format it allows compatibility with many more applications.
This image type combines both the TIFF format with the OME-XML format.


Image Specifications
---------------------

CellOrganizer training demos on MATLAB accept the following image file formats:
 * TIFF
 * OME-TIFF.
 * The minimum number of images to run a model is 20.
 * Images need to be processed as OME.TIFF format before they can be uploaded
 * Dataset if it is a collection should be a mix of HeLa Cells
 * Images must be static 2D or 3D dimensionality
 * Type (file): .mat, .OME.TIFF, .pdf, tabular, .tiff, .txt, vcml, xlsx, xml, zip
 * Import image from a URL - if images are in a repository just add the url for it.

.. IMPORTANT::
    if you need to change the format of your images to OME.TIFF check `here <https://www-legacy.openmicroscopy.org/site/products/ome-tiff>`_ for more on what you will need to do. In most demos, CellOrganizer anticipates that you have a unique image dataset for each of the following three channels: cellular shapes, nucleoli, and proteins.
    You must set each of these variables individually and can choose to remove or add variables to comply with your own dataset. To set a variable, you are expected to either provide a list of the image filenames or a pattern of the image filenames.
    **All images in the dataset are to be binarized and contrast – stretched prior to the main processing step.**

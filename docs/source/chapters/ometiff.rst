OME-TIFF
--------

The software supports saving synthetic images as valid `OME.TIFF <https://docs.openmicroscopy.org/ome-model/5.6.3/#ome-tiff>`_. Currently, the software can only save images with one protein pattern as OME.TIFFs. To do this, use the output flag

.. code:: matlab

	options.output.OMETIFF = true;

For an example, investigate and run demo3D34.

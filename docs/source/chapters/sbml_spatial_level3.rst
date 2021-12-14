SBML Spatial Level 3
--------------------

The software supports saving synthetic images as valid `SBML Level 3 Spatial <http://sbml.org/Documents/Specifications/SBML_Level_3/Packages/spatial>`_ instances. Currently, the software can only save images with one vesicle pattern as OME.TIFFs. To do this, use the output flag

.. code:: matlab

	options.output.SBMLSpatial = true;

For an example, investigate and run demo3D34.

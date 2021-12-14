slml2img
********
This function synthesizes an image from a list of SLML models.
Instances may be saved in the following forms:

#. tiff stacks: a 3D tiff image stack for each pattern generated using the input models
#. OME.TIFF: saves synthetic images as valid `OME.TIFF <https://docs.openmicroscopy.org/ome-model/5.6.3/#ome-tiff>`_
#. indexed images: a single 3D tiff image stack where each pattern is represented by a number 1-n
#. object mesh: a .obj mesh file for each pattern generated using the input models (blenderfile option)
#. SBML-Spatial file: a Systems Biology Markup Language (SBML) instance XML file utilizing the Spatial extension in level 3 version 1
#. Virtual Cell Markup-Language (VCML): a Virtual Cell Markup Language and the native Virtual Cell format.

**Example:**

*slml2img(models, options)*

=======================  ========================================
List Of Input Arguments  Descriptions
=======================  ========================================
models                   A cell array of filenames
options                  A structure holding the function options
=======================  ========================================



General Options
==============

*options.synthesis* (mandatory)
    * Synthesis parameter that allows to synthesize 'nucleus', 'cell', 'framework' or 'all'.

*options.debug* (optional) **[false]**
    * If set to true, then the function will (1) keep temporary results folder, (2) will print information useful for debugging.

*options.targetDirectory* (optional) **[current]**
    * Directory where the images are going to be saved.

*options.prefix* (optional) **['demo']**
    * Filename prefix for the synthesized images.

*options.numberOfSynthesizedImages* (optional) **[1]**
    * Number of synthesized images.

*options.image_size* (optional) **[1024 1024]**
    * The image size is [1024 1024] for both 2D and 3D in x and y.

*options.compression* (optional) **['lzw']**
    * Compression of tiff, i.e. 'none', 'lzw' and 'packbits'

*options.microscope* (optional) **['none']**
    * Microscope model from which we select a point spread function.

*options.sampling.method* (optional) **['trimmed']**
    * Can be 'disc', 'sampled' or 'trimmed'.

*options.resolution.cell* (optional)
    * The resolution of the cell and nucleus that are being passed in

*options.instance.cell* (optional) **[empty]**
    * A binary cell image to be filled with objects.

*options.instance.nucleus* (optional) **[empty]**
    * A binary nuclear image to be filled with objects.


Protein Submodel Options
^^^^^^^^^^^^^^^^^^^^^^^^
*options.overlapsubsize* (optional) **[0.3]**
    * Defines the downsampling fraction to perform during object overlap avoidance.

*options.overlapthresh* (optional) **[0]**
    * Defines the amount of overlap that is allowed between objects.

*options.protein.cytonuclearflag* (optional) **[cytonuclearflag included in the model]**
    * Defines the allowable region for protein placement.

*options.oobbuffer* (optional) **[0]**
    * The thickness in microns of an additional buffer zone inside the boundary of a cell in which an object cannot be placed.

*options.rendAtStd* (optional) **[2]**
    * Defines the number of standard deviations to render Gaussian objects at.

*options.resolution.objects* (optional)
    * The resolution of the object model being synthesized


Model Specific Options
======================

2D PCA
^^^^^^^^
learn more `here <https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty983/5232995>`_

*options.model.pca.pca_synthesis_method* (mandatory) **['reconstruction' or 'random sampling']**
    * The method in which the generated image is synthesized.

*options.model.pca.imageSize* (mandatory) **[1024, 1024]**
    * image size of the resulting synthesized image


3D SPHARM-RPDM
^^^^^^^^^^^^^^^
learn more `here <https://link.springer.com/protocol/10.1007%2F978-1-4939-9102-0_11>`_

*options.model.spharm_rpdm.synthesis_method* (mandatory) **['reconstruction' or 'random sampling']**


3D T-Cell Distribution
^^^^^^^^^^^^^^^^^^^^^^
learn more `here <https://link.springer.com/protocol/10.1007/978-1-4939-6881-7_25>`_

*options.model.tcell.results_location* (mandatory)
    * File path for where the results should be saved.

*options.model.tcell.named_option_set* (mandatory)
    * The running choice for CellOrganizer and one sensor of two-point annotation

*options.model.tcell.sensor* (mandatory)
    * Set up protein name

*options.model.tcell.model_type_to_include* (mandatory)
    * Set up for model to include

*options.model.tcell.use_two_point_synapses* (optional)
    * Set up the mode of synapse to use, as a default, we use one-point, if needed you can use two-point by set up the option as true

*options.model.tcell.timepoints_to_include* (optional)
    * If creation of models for only a subset of the time points is desired, edit to specify which time points to include

Output Options
==============
OMETIFF
^^^^^^^
*options.output.ometiff* (optional) **[false]**
    * Boolean flag specifying whether to write out an (.ome.tif) OME TIFF.
SBML
^^^^
*options.output.SBML* (optional) **[false]**
    * Boolean flag specifying whether to write out (.xml) files with SBML-Spatial 2 representations of geometries. Default is false.
    
*options.output.SBML.downsampling* (optional) **[1]**
    * Downsampling fraction for the creation of SBML Spatial files when output.SBML or output.SBMLSpatial are true (1 means no downsampling, 1/5 means 1/5 the size).
    
*options.output.SBML.spatial* (optional) **[false]**
    * Boolean flag specifying whether to write out (.xml) file with SBML-Spatial 3 representations of geometries. Default is false.
    
*options.output.SBML.spatialimage* (optional) **[false]**
    * Boolean flag specifying whether SBML-Spatial 3 output represents geometries with image volumes instead of meshes. Meshes are not supported by Virtual Cell. Default is false.
    
*options.output.SBML.spatialusecompression* (optional) **[true]**
    * Boolean flag specifying whether to write SBML Spatial output using compression. Default is true.
    
*options.output.SBML.spatialuseanalytic_meshes* (optional) **[false]**
    * Boolean flag specifying whether to use analytic meshes instead of isosurfaces of rasterized shapes. Default is false.
    
*options.output.SBML.spatialvcellcompatible* (optional) **[false]**
    * Boolean flag specifying whether to write SBML Spatial output compatible with Virtual Cell but not the Level 3 Version 1 Release 0.90 draft specifications. Default is false.

*options.output.SBML.translations* (optional) **[{}]**
    * N x 2 cell array of strings (first column) to be replaced by other strings (second column) in CellOrganizer-generated SBML.

VCML
^^^^
*options.output.VCML.writeVCML* (optional) **[false]**
    * Boolean flag specifying whether to write out VCML files for use with Virtual Cell.
    
*options.output.VCML.input_filename* (optional) **[false]**
    * String specifying Virtual Cell VCML file with biochemistry which will be combined with generated geometry in output file.

*options.output.VCML.downsampling* (optional) **[1]**
    * Downsampling fraction for the creation of object files (1 means no downsampling, 1/5 means 1/5 the size).

*options.output.VCML.addTranslocationIntermediates* (optional) **[true]**
    * Boolean flag specifying whether to create intermediate species and reactions for reactions involving non-adjacent translocations, which are valid in cBNGL but not Virtual Cell.

*options.output.VCML.numSimulations* (optional)  **[1]**
    * Number of simulations in VCML file.

*options.output.VCML.translations* (optional) **[{0,2}]**
    * N x 2 cell array of strings (first column) to be replaced by other strings (second column).

*options.output.VCML.defaultDiffusionCoefficient* (optional) **[1.0958e-11]**
    * Double specifying diffusion coefficient in meters squared per second.

*options.output.VCML.NET.filename* (optional) **[' ']**
    * String specifying BioNetGen network file to include in VCML files for use with Virtual Cell.

*options.output.VCML.NET.units.concentration* (optional) **['uM']**
    * String specifying concentration units in NET file.

*options.output.VCML.NET.units.length* (optional) **['um']**
    * String specifying length units in NET file.

*options.output.VCML.NET.units.time* (optional) **['s']**
    * String specifying time units in NET file.

*options.output.VCML.NET.effectiveWidth* (optional) **[3.8775e-9]**
    * Double specifying surface thickness in meters.

*options.output.VCML.NET.useImageAdjacency* (optional) **[true]**
    * Boolean specifying whether to derive compartment adjacency from the synthetic image. Can break Virtual Cell compatibility due to inclusion of BioNetGen representation of translocation between non-adjacent compartments.

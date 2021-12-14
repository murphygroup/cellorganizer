.. api:

Main Functions
==============


.. note::

   This section only includes the headers of the most important functions
   in CellOrganizer.

img2slml
********

Method header::

  function answer = img2slml( varargin )
  % IMG2SLML Trains a generative model of subcellular location from a
  % collection of images and saves the model to disk.
  %
  % A CellOrganizer model consists of four components,
  %
  % 1) a (optional) documentation component
  % 2) a nuclear membrane model,
  % 3) a cell membrane model and,
  % 4) a protein pattern model.
  %
  % ┌────────────────────────┐
  % │List Of Input Parameters│
  % └────────────────────────┘
  % Inputs                      Descriptions
  % ------                      ------------
  % dimensionality              2D/3D
  % dnaImagesDirectoryPath      DNA images collection directory, list of files or pattern
  % cellImagesDirectoryPath     Cell images collection directory, list of files or pattern
  % proteinImagesDirectoryPath  Protein images collection directory, list of files or pattern
  % options                     List of options
  %
  % The input argument options holds the valid parameters for all of these components.
  %
  % ┌───────────────┐
  % │List Of Options│
  % └───────────────┘
  % Mandatory options         Descriptions
  % -----------------         ------------
  % model.resolutions         Any double 1x2/1x3 double vector.
  %                           (microns/voxel).
  %
  % Generic model options     Descriptions
  % ---------------------     ------------
  % masks                     (optional) Masks collection directory.
  %
  % train.flag                (optional) Selects what model is going to be trained ('nuclear',
  %                           'framework', or 'all'). Default is 'all'.
  %
  % model.name                (optional) Holds the name of the model. Default is empty.
  % model.id                  (optional) Holds the id of the model. Default is a randomly generated string.
  % model.filename            Holds the output filename.
  % downsampling              Downsampling vector to be used during preprocessing.
  %
  % Debugging options
  % -----------------
  % debug                     If set to true, then the function will (1) keep temporary results folder, (2) will
  %                           print information useful for debugging. Default is false.
  % display                   If set to true, then plots useful for debugging with be open. This functionality is
  %                           meant for debugging only, setting this to true will considerably slow down
  %                           computation. Default is false;
  % save_segmentations        Will save the segmentations to the model file. Setting this option to true will create
  %                           a considerably large file.
  %
  % Nuclear shape model options  Descriptions
  % ---------------------------  ------------
  % nucleus.class                (mandatory) Holds the nuclear membrane model class.
  % nucleus.type                 (mandatory) Holds the nuclear membrane model type.
  % nucleus.name                 (optional) Holds the name of the nuclear model. Default is empty.
  % nucleus.id                   (optional) Holds the id of the nuclear model. Default is a randomly generated string.
  %
  % Cell shape model options  Descriptions
  % ------------------------  ------------
  % cell.class                (mandatory) Holds the cell membrane model class.
  % cell.type                 (mandatory) Holds the cell membrane model type.
  % cell.name                 (optional) Holds the name of the cell model. Default is empty.
  % cell.id                   (optional) Holds the id the cell model. Default is empty.
  %
  % Protein shape model options  Descriptions
  % ---------------------------  ------------
  % protein.class                (mandatory) Holds the protein model class.
  % protein.type                 (mandatory) Holds the protein model type.
  % protein.name                 (optional) Holds the name of the protein model. The default is empty.
  % protein.id                   (optional) Holds the id of the protein model. Default is a randomly generated string.
  % protein.cytonuclearflag      (optional) Determines whether the protein pattern will be generated in
  %                              the cytosolic space ('cyto'), nuclear space ('nuc') or everywhere ('all').
  %                              Default is cyto.
  %
  % ┌────────────────────────────────────┐
  % │List Of Options per Model class/type│
  % └────────────────────────────────────┘
  % 2D PCA model options
  % --------------------
  % model.pca.latent_dim      (optional) This specifies how many latent dimensions should be used for modeling
  %                           the shape space. Valid values arepositive integers. The default is 15.
  %
  % 2D diffeomorphic model options
  % ------------------------------
  % model.diffeomorphic.distance_computing_method     (optional) ‘faster'
  % model.diffeomorphic.com_align                     (optional) 'nuc'
  %
  % T cell distribution model options
  % ---------------------------------
  % model.tcell.synapse_location            (mandatory) File path to annotation of the synapse positions of the T cells as input.
  % model.tcell.results_location            (mandatory) File path for where the results should be saved.
  % model.tcell.named_option_set            (mandatory) The running choice for CellOrganizer and one sensor of two-point annotation.
  % model.tcell.use_two_point_synapses      (optional) Set up the mode of synapse to use, as a default, we use one-point,
  %                                         if needed you can use two-point by set up the option as true.
  % model.tcell.sensor                      Set up protein name.
  % model.tcell.timepoints_to_include       (optional) If creation of models for only a subset of the time points is desired,
  %                                         edit to specify which time points to include.
  % model.tcell.model_type_to_include       (mandatory) Set up for model to include.
  % model.tcell.infer_synapses              (mandatory) true or false.
  % model.tcell.adjust_one_point_alignment  (optional) Set up alignment adjustment true or false.
  % model.tcell.ometiff                     (optional) If true, then it assumes images are OME.TIFFs with annotations. Default is false.
  %
  % 3D SPHARM-RPDM model options
  % ----------------------------
  % model.spharm_rpdm.alignment_method 	(optional) method by which cells willbe aligned when producing shape descriptors
  %                                       The possible values are 'major_axis' (defaut) or 'foe'.
  % model.spharm_rpdm.rotation_plane 	(optional)  Dimensions of image that will used for alignment.
  %                                       The possible values are 'xy' (defaut), 'xz', 'yz' or ‘xyz'. For example,
  %                                       if ‘xy‘ is specified, each cellwill be rotated in the 	xy plane (around the z axis).
  % model.spharm_rpdm.postprocess         (optional) This specifies whether alignment and size normalization
  %                                       should be done after parameterization.  The values are ‘true’ (default) and ‘false’.
  % model.spharm_rpdm.maxDeg              (optional) This specifies the degree up to which spherical harmonics
  %                                       should be calculated.  Valid values are positive integers.  The default is 31.
  % model.spharm_rpdm.components          (mandatory) This specifies which components should be included in the
  %                                       shape model.  The valid values are {'cell'}, {'nuc'}, or {'cell', 'nuc'}.
  % model.spharm_rpdm.latent_dim          (optional) This specifies how many latent dimensions should be used for
  %                                       modeling the shape space.  Valid values are positive integers.  The default is 15.
  %
  % ┌─────────────┐
  % │Documentation│
  % └─────────────┘
  % This is an optional structure with multiple elements that holds documentation about this model.
  %
  % documentation.<name>      Holds the value of variable <name>. This is meant to be meta information. Default is empty.
  %
  % Helper Options
  % -------------
  % verbose                   (optional) Displays messages to screen. Default is true.
  % debug                     (optional) Reports errors and warnings. Default is false.

slml2info
*********

Method header::

  function answer = slml2info( varargin )
  % SLML2INFO Generate a report from information extracted from a genearative model file
  %
  % List Of Input Arguments  Descriptions
  % -----------------------  ------------
  % filenames                List of files
  % options                  Options structure
  %
  % Example
  % > filenames = {'/path/to/model/file/model.mat'};
  % > answer = slml2info( filenames );

slml2img
********

Method header::

  % SLML2IMG Synthesizes an image from a list of SLML models.
  %
  % Instances may be saved in the following forms:
  % a) tiff stacks: a 3D tiff image stack for each pattern generated using the input models
  % b) indexed images: a single 3D tiff image stack where each pattern is represented by a number 1-n
  % c) object mesh: a .obj mesh file for each pattern generated using the input models (blenderfile option)
  % d) SBML-Spatial file: a Systems Biology Markup Language (SBML) instance XML file utilizing the Spatial extension in level 3 version 1
  %
  %
  % List Of Input Arguments  Descriptions
  % -----------------------  ------------
  % models                   A cell array of filenames
  % options                  A structure holding the function options
  %
  % The shape of options is described
  %
  % List Of Options           Descriptions
  % ---------------           ------------
  % targetDirectory           (optional) Directory where the images are going to be saved. Default is current directory.
  % prefix                    (optional) Filename prefix for the synthesized images. Default is 'demo'
  % numberOfSynthesizedImages (optional) Number of synthesized images. Default is 1.
  % compression               (optional) Compression of tiff, i.e. 'none', 'lzw' and 'packbits'
  % microscope                (optional) Microscope model from which we select a point spread function. Default is 'none'
  % synthesis                 (optional) Synthesis parameter that allows to
  %                                      synthesize 'nucleus', 'framework' or 'all'. Default is 'all'
  % protein.cytonuclearflag   (optional) Defines the allowable region for protein placement.
  %                                      The default is the cytonuclearflag included in the model.
  % sampling.method           (optional) Can be 'disc', 'sampled' or 'trimmed'. Default is trimmed
  % savePDF                   (optional) Saves the probability density function for a given pattern during 2D synthesis. Default is false.
  % spherical_cell            (optional) Boolean flag that indicates whether a cell is spherical. Default is false.
  % overlapsubsize            (optional) Defines the downsampling fraction to perform during object overlap avoidance. Default is 0.3.
  % overlapthresh             (optional) Defines the amount of overlap that is allowed between objects. Default is 1.
  % rendAtStd                 (optional) Defines the number of standard deviations to render Gaussian objects at. Default is 2.
  % sampling.method.density   (optional) An integer. Default is empty.
  % protein.cytonuclearflag   (optional) Can 'cyto', 'nucleus' or 'all'. Default is all.
  % resolution.cell           (optional) The resolution of the cell and nucleus that are being passed in
  % resolution.objects        (optional) The resolution of the object model being synthesized
  % instance.cell             (optional) A binary cell image to be filled with objects. Default is empty.
  % instance.nucleus          (optional) A binary nuclear image to be filled with objects. Default is empty.
  % image_size                (optional) The image size. Default is [1024 1024] for both 2D and 3D in x and y
  % synthesis.diffeomorphic.maximum_iterations (optional) Integer defining the maximum number of iterations during diffeo inference. Default is 100.
  %
  % Random walk options
  % -------------------
  % randomwalk                (optional) Boolean flag of whether to perform a shape space walk. Default is False.
  % framefolder               (optional) The folder in which to look for completed frames and save finished frames from the diffeomorphic synthesis.
  %                                      The default is './frames/'.
  % walksteps                 (optional) The integer number of steps to walk during a shape space walk. Default is 1.
  % walk_type                 (optional) Type of random walk to perform. Default is 'willmore'.
  %
  % Helper options
  % --------------
  %
  % debug                     (optional) Keeps temporary results and catches
  %                           errors with full reports. Default is false;
  % display                   (optional) Will make pretty plots. Turning this
  %                           flag on will slow down synthesis. Default is
  %                           false.
  % verbose                   (optional) Print the intermediate steps to screen. Default is false.
  %
  % Outputs
  % -------
  % output.tifimages           (optional) Boolean flag specifying whether to write out tif images. Default is true.
  % output.indexedimage        (optional) Boolean flag specifying whether to write out indexed image. Default is false.
  % output.blenderfile         (optional) Boolean flag specifying whether to write out (.obj) files for use in blender. Default is false;
  % output.blender.downsample  (optional) ownsampling fraction for the creation of object files (1 means no downsampling, 1/5 means 1/5 the size).
  % output.SBML                (optional) boolean flag specifying whether to write out (.xml) files with SBML-Spatial representations of geometries. Default is false;

slml2report
***********

Method header::

  % SLML2REPORT Generate a report comparing two SLML generative models
  %
  % List Of Input Arguments  Descriptions
  % -----------------------  ------------
  % model1                   A generative model filename 
  % model2                   A generative model filename
  %
  % Example
  % > filename1 = '/path/to/model/model1.mat';
  % > filename2 = '/path/to/model/model2.mat';
  % answer = slml2report( filename1, filename2 );

slml2slml
*********

Method header::

  function model = slml2slml( varargin )
  % SLML2SLML Combines multiple SLML files into a single model file.
  %
  % List Of Input Arguments     Descriptions
  % -----------------------     ------------
  % files                       list of paths of models need be combined
  % options                     Options structure
  %
  % The input argument options holds the valid parameters for these components.
  % The shape of options is described below
  %
  % List Of Parameters        Descriptions
  % ------------------        ------------
  % output_filename           (optional)the file name of output model,
  %                           default is "model.mat"
  % seletion
  % name
  %
  % Documentation (optional)
  % ------------------------
  % This is an optional structure with multiple elements that holds documentation about this model.
  % If the decumentation isn't input, function will inherit documentation from first
  % model in list if model is present


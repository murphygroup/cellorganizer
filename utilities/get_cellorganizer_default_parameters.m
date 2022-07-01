function param = get_cellorganizer_default_parameters( option, param )
%GET_CELLORGANIZER_DEFAULT_PARAMETERS This helper function helps set up the
%default parameters for CellOrganizer

% Ivan E. Cao-Berg
%
% Copyright (C) 2015-2020 Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

% icaoberg 2018/10/26 Added checks for param.train.flag = 'tcell'
% xruan 05/09/2019 Added synthesis resolution for spharm model

%these are the default baby parameters for cellorganizer
param = ml_initparam( param, struct( ...
    'debug', false, ...
    'display', false , ...
    'verbose', true ));

param = ml_initparam( param, struct( ...
    'slml', struct( 'level', '2' )));

param = ml_initparam( param, struct( ...
    'slml', struct( 'version', '8.0' )));

param = ml_initparam( param, struct( ...
    'randomwalk', false ));

if isfield( param, 'model' ) && ...
        ( ~isfield( param.model, 'id' ) || ...
        ( isfield( param.model, 'id' ) && isempty( param.model.id ) ) )
    param.model.id = uuidgen();
end

if strcmpi( option, 'training' )
    param = training( param );
end

if strcmpi( option, 'synthesis' )
    param = synthesis( param );
end
end%get_cellorganizer_default_parameters

function param = synthesis( param )
    param = ml_initparam( param, struct( ...
        'temporary_results', [ pwd filesep 'temp' ] ));

    %these are optional parameters that are generally hidden from users
    param = ml_initparam( param, struct( ...
        'overwrite_synthetic_instances', true, ...
        'clean_synthetic_instances', true, ...
        'numberOfSynthesizedImages', 1, ...
        'targetDirectory', pwd, ...
        'compression', 'none', ...
        'prefix', 'demo', ...
        'synthesis', 'all' ));

    param = ml_initparam( param, struct( ...
        'overlapsubsize', 0.3, ...
        'overlapthresh', 1, ...
        'rendAtStd', 2 ));

    param = ml_initparam( param, struct( ...
        'oobthresh', 0, ...
        'oobbuffer', 0 ));

    param = ml_initparam( param, struct( ...
        'framework_cropping', true ...
        )); % Debug

    if param.numberOfSynthesizedImages <= 0
        warning( ['Number of synthesized images is ' ...
            param.numberOfSynthesizedImages '. Setting default value to 1.'] );
        param.numberOfSynthesizedImages = 1;
    end

    if isempty(param.prefix)
        warning('Prefix cannot be empty. Setting default value to ''demo''');
        param.prefix = 'demo';
    end

    if ~isfield( param, 'output' )
        param.output = struct();
    end
    
    % reset SBML synthesis options into Output
    if isfield(param,'SBML')
        if isfield(param.SBML,'SBML')
            param.output.SBML = param.SBML.SBML;
        end

        if isfield(param.SBML,'downsampling')
            param.output.SBMLDownsampling = param.SBML.downsampling;
        end

        if isfield(param.SBML,'spatial_image')
            param.output.SBMLSpatialImage = param.SBML.spatial_image;
        end

        if isfield(param.SBML,'spatial_use_compression')
            param.output.SBMLSpatialUseCompression = param.SBML.spatial_use_compression;
        end

        if isfield(param.SBML,'spatial_use_analytic_meshes')
            param.output.SBMLSpatialUseAnalyticMeshes = param.SBML.spatial_use_analytic_meshes;
        end

        if isfield(param.SBML,'spatial_vcell_compatible')
            param.output.SBMLSpatialVCellCompatible = param.SBML.spatial_vcell_compatible;
        end

        if isfield(param.SBML,'include_ec')
            param.output.SBMLIncludeEC = param.SBML.include_ec;
        end

        if isfield(param.SBML,'ec_scale')
            param.output.SBMLECScale = param.SBML.ec_scale;
        end

        if isfield(param.SBML,'flip_x_to_align')
            param.output.SBMLFlipXToAlign = param.SBML.flip_x_to_align;
        end

        if isfield(param.SBML,'translations')
            param.output.SBMLTranslations = param.SBML.translations;
        end

        if isfield(param.SBML,'end_time')
            param.output.SBMLEndTime = param.SBML.end_time;
        end

        if isfield(param.SBML,'default_time_step')
            param.output.SBMLDefaultTimeStep = param.SBML.default_time_step;
        end

        if isfield(param.SBML,'max_time_step')
            param.output.SBMLMaxTimeStep = param.SBML.max_time_step;
        end

        if isfield(param.SBML,'output_time_step')
            param.output.SBMLOutputTimeStep = param.SBML.output_time_step;
        end
    end
    
    if isfield(param, 'output') && isfield(param.output, 'VCML') && isa(param.output.VCML,'logical')
        temp_VCML = param.output.VCML;
        param.output = rmfield(param.output,'VCML');
        param.output.VCML.writeVCML = temp_VCML;
    end

    if isfield(param, 'VCML')
        if ~isfield( param.output, 'VCML' )
            param.output.VCML = struct();
        end
        param.output.VCML = ml_initparam( param.output.VCML, param.VCML);
        param = rmfield(param, 'VCML');
    end
    if isfield(param.output, 'writeVCML')
        param.output.VCML.writeVCML = param.output.writeVCML;
    end

    if ~isfield( param.output, 'MCellMDL' )
        param.output.MCellMDL = struct();
    end
    if isfield(param, 'MCellMDL')
        param.output.MCellMDL = ml_initparam( param.output.MCellMDL, param.MCellMDL);
        param = rmfield(param, 'MCellMDL');
    end
    if isfield(param.output, 'writeMCellMDL')
        param.output.MCellMDL.writeMCellMDL = param.output.writeMCellMDL;
    end
    
    if isfield(param, 'NET')
        if isfield(param.NET, 'filename')
            param.output.NET.filename = param.NET.filename;
        end
        if isfield(param.NET, 'translations')
            param.output.NET.translations = param.NET.translations;
        end
        if ~isfield( param.output, 'NET' )
            param.output.NET = struct();
        end
        if isfield(param.NET, 'units.concentration')
            param.output.NET.units.concentration = param.NET.units.concentration;
        end
        if isfield(param.NET, 'units.time')
            param.output.NET.units.time = param.NET.units.time;
        end
        if isfield(param.NET, 'effectiveWidth')
            param.output.NET.effectiveWidth = param.NET.effectiveWidth;
        end
        param.output.NET = ml_initparam( param.output.NET, param.NET);
        param = rmfield(param, 'NET');
    end
    
    for temp = {'framework_min_clearance', 'framework_clearance_n_max_filter_rounds', 'intersecting_mesh_object_policy'}'
        temp = temp{1};
        if isfield(param, temp)
            param.output.(temp) = param.(temp);
            param = rmfield(param, temp);
        end
    end

    %%%% Default Output Options %%%%%
    % output.tifimages                                    (optional) Boolean flag specifying whether to write out tif images. Default is true.
    % output.indexedimage                                 (optional) Boolean flag specifying whether to write out indexed image. Default is false.
    % output.blenderfile                                  (optional) Boolean flag specifying whether to write out (.obj) files for use in blender. Default is false.
    % output.meshes                                       (optional) Boolean flag specifying whether to write out (.obj) files of analytic meshes (if available, does not work with every model type). Default is false.
    % output.shape_space_coords                           (optional) Boolean flag specifying whether to write out (.mat, .txt) files containing shape space coordinates. Currently only for SPHARM geometry. Default is false.
    % output.framework_min_clearance                      (optional) double specifying the minimum distance in μm to impose between nucleus and cell after synthesis. -inf to disable. Currently only used for framework meshes assuming corresponding vertices by instance2MCellMDL. Default is -inf.
    % output.framework_clearance_n_max_filter_rounds      (optional) integer specifying the number of rounds of maximum filter to apply to the projections of cell vertices onto nucleus normals among the immediate neighbors of each vertex. 0 to disable. Currently only used for framework meshes assuming corresponding vertices by instance2MCellMDL. Default is 1.
    % output.intersecting_mesh_object_policy              (optional) string specifying policy for checking framework and objects for intersection and whether to remove objects or reject the synthesized cell entirely. Currently untested for values other than 'ignore'. Currently only used for framework meshes assuming corresponding vertices by instance2MCellMDL. Default is 'ignore'.
    % output.blender.downsample                           (optional) downsampling fraction for the creation of object files (1 means no downsampling, 1/5 means 1/5 the size).
    % output.SBML                                         (optional) boolean flag specifying whether to write out (.xml) files with SBML-Spatial 2 representations of geometries. Default is false.
    % output.SBMLDownsampling                             (optional) downsampling fraction for the creation of SBML Spatial files when output.SBML or output.SBMLSpatial are true (1 means no downsampling, 1/5 means 1/5 the size). Default is 1.
    % output.SBMLSpatial                                  (optional) boolean flag specifying whether to write out (.xml) file with SBML-Spatial 3 representations of geometries. Default is false.
    % output.SBMLSpatialImage                             (optional) boolean flag specifying whether SBML-Spatial 3 output represents geometries with image volumes instead of meshes. Meshes are not supported by Virtual Cell. Default is false.
    % output.SBMLSpatialUseCompression                    (optional) boolean flag specifying whether to write SBML Spatial output using compression. Default is true.
    % output.SBMLSpatialUseAnalyticMeshes                 (optional) boolean flag specifying whether to use analytic meshes instead of isosurfaces of rasterized shapes. Default is false.
    % output.SBMLSpatialVCellCompatible                   (optional) boolean flag specifying whether to write SBML Spatial output compatible with Virtual Cell but not the Level 3 Version 1 Release 0.90 draft specifications. Default is false.
    % output.SBMLSpatialImageDownsampling                 (optional) downsampling fraction for the creation of SBML Spatial files when output.SBMLSpatialImage is true (1 means no downsampling, 1/5 means 1/5 the size). Default is 1.
    % output.SBMLTranslations                             (optional) N x 2 cell array of strings (first column) to be replaced by other strings (second column) in CellOrganizer-generated SBML. Default is {}.
    % output.SBMLIncludeEC                                (optional) boolean flag specifying whether to include an extracellular region in SBML Spatial output. Default is false.
    % output.SBMLECScale                                  (optional) scaling for extracellular region in SBML Spatial output relative to bounding box for cell shape. Default is 1.
    % output.SBMLDefaultDiffusionCoefficient              (optional) double specifying diffusion coefficient in meters squared per second. Default is 1.0958e-11.
    % output.VCML.writeVCML                               (optional) boolean flag specifying whether to write out VCML files for use with Virtual Cell. Default is false.
    % output.VCML.input_filename                          (optional) string specifying Virtual Cell VCML file with biochemistry which will be combined with generated geometry in output file. Default is empty string.
    % output.VCML.downsampling                            (optional) downsampling fraction for the creation of object files (1 means no downsampling, 1/5 means 1/5 the size). Default is 1.
    % output.VCML.add_translocation_intermediates         (optional) boolean flag specifying whether to create intermediate species and reactions for reactions involving non-adjacent translocations, which are valid in cBNGL but not Virtual Cell. Default is true.
    % output.VCML.translations                            (optional) N x 2 cell array of strings (first column) to be replaced by other strings (second column).
    % output.VCML.default_diffusion_coefficient           (optional) double specifying diffusion coefficient in meters squared per second. Default is 1.0958e-11.
    % output.VCML.NET.filename                            (optional) string specifying BioNetGen network file to include in VCML files for use with Virtual Cell. Default is empty string.
    % output.VCML.NET.units.concentration                 (optional) string specifying concentration units in NET file. Default is 'uM'.
    % output.VCML.NET.units.length                        (optional) string specifying length units in NET file. Default is 'um'.
    % output.VCML.NET.units.time                          (optional) string specifying time units in NET file. Default is 's'.
    % output.VCML.NET.effectiveWidth                      (optional) double specifying surface thickness in meters. Default is 3.8775e-9.
    % output.VCML.NET.useImageAdjacency                   (optional) boolean specifying whether to derive compartment adjacency from the synthetic image. Can break Virtual Cell compatibility due to inclusion of BioNetGen representation of translocation between non-adjacent compartments. Default is true.
    % output.VCML.num_simulations                         (optional) number of simulations in VCML file.
    % output.VCML.delete_input_simulations                (optional) boolean specifying whether to delete simulations in VCML file or to modify them using `output.VCML.end_time`, etc. Default is false (behavior changed since version 2.9.0).
    % output.MCellMDL.writeMCellMDL                       (optional) boolean flag specifying whether to write out MCellMDL files for use with MCell. Default is false.
    % output.MCellMDL.downsampling                        (optional) downsampling fraction for the creation of object files (1 means no downsampling, 1/5 means 1/5 the size). Default is 1.
    % output.MCellMDL.addTranslocationIntermediates       (optional) boolean flag specifying whether to create intermediate species and reactions for reactions involving non-adjacent translocations, which are valid in cBNGL but not MCell. Default is true.
    % output.MCellMDL.numSimulations                      (optional) number of simulations in MCellMDL file.
    % output.MCellMDL.translations                        (optional) N x 2 cell array of strings (first column) to be replaced by other strings (second column).
    % output.MCellMDL.defaultDiffusionCoefficient         (optional) double specifying diffusion coefficient in meters squared per second. Default is 1.0958e-11.
    % output.MCellMDL.input_filename_pattern              (optional) string specifying pattern matching a set of MCell MDL files to be combined with generated MDL files. This should be empty or `[path][prefix].*.mdl`. If not empty, CellOrganizer will only generate the geometry file and will copy the other files matching the pattern to the output directory, and it is the user's responsibility to ensure compatibility between the input and CellOrganizer's output. Default is `''`.
    % output.NET.filename                                 (optional) string specifying BioNetGen network file to include in VCML or MCell MDL files for use with Virtual Cell or MCell MDL files for MCell. Default is `''`. Only one of `output.MCellMDL.input_filename_pattern` and `output.NET.filename` can be non-empty.
    % output.NET.units.concentration                      (optional) string specifying concentration units in NET file. Default is 'uM'.
    % output.NET.units.length                             (optional) string specifying length units in NET file. Default is 'um'.
    % output.NET.units.time                               (optional) string specifying time units in NET file. Default is 's'.
    % output.NET.translations                             (optional) N x 2 cell array of strings (first column) to be replaced by other strings (second column).
    % output.NET.effective_width                          (optional) double specifying surface thickness in meters. Default is 3.8775e-9.
    % % output.NET.use_image_adjacency                    (optional) boolean specifying whether to derive compartment adjacency from the synthetic image. Can break Virtual Cell compatibility due to inclusion of BioNetGen representation of translocation between non-adjacent compartments. Default is true.
    % output.NET.downsampling                             (optional) downsampling fraction for the creation of object files (1 means no downsampling, 1/5 means 1/5 the size). Default is 1.
    % output.NET.add_translocation_intermediates          (optional) boolean flag specifying whether to create intermediate species and reactions for reactions involving non-adjacent translocations, which are valid in cBNGL but not Virtual Cell. Default is true.
    % output.NET.default_diffusion_coefficient            (optional) double specifying diffusion coefficient in meters squared per second. Default is 1.0958e-11.
    % output.OMETIFF                                      (optional) boolean flag specifying whether to write out an (.ome.tif) OME TIFF. Default is false.
    
    % From Luby-Phelps 2000, "Cytoarchitecture and physical properties of cytoplasm: volume, viscosity, diffusion, intracellular surface area," Table I, column "Cytoplasmic D (μm2/sec)," excluding entries containing only upper bounds and taking the average of the bounds of entries containing ranges:
    default_diffusion_coefficient = (0.9 + 6.9 + 43 + 27 + 5.9 + 1.7 + 6.8 + (7+11)/2 + 13.5 + (6+11)/2 + 6.7 + 1.6) / 12; % um2.s-1
    
    output_defaults = struct();
    output_defaults.tifimages = true;
    output_defaults.indexedimage = false;
    output_defaults.blenderfile = false;
    output_defaults.meshes = false;
    output_defaults.shape_space_coords = false;
    output_defaults.framework_min_clearance = -inf;
    output_defaults.framework_clearance_n_max_filter_rounds = 1;
    output_defaults.intersecting_mesh_object_policy = 'ignore';
    output_defaults.SBML = false;
    output_defaults.SBMLDownsampling = 1;
    output_defaults.SBMLSpatial = false;
    output_defaults.SBMLSpatialImage = false;
    output_defaults.SBMLSpatialUseCompression = true;
    output_defaults.SBMLSpatialUseAnalyticMeshes = false;
    output_defaults.SBMLSpatialVCellCompatible = false;
    output_defaults.SBMLSpatialImageDownsampling = 1;
    output_defaults.SBMLTranslations = {};
    output_defaults.SBMLIncludeEC = true;
    output_defaults.SBMLECScale = 1;
    output_defaults.SBMLDefaultDiffusionCoefficient = default_diffusion_coefficient;
    output_defaults.SBMLEndTime = 20;
    output_defaults.SBMLDefaultTimeStep = 0.05;
    output_defaults.SBMLMaxTimeStep = 0.1;
    output_defaults.SBMLOutputTimeStep = 0.5;
    output_defaults.SBMLFlipXToAlign = true; % Debug
    output_defaults.VCML = struct();
    output_defaults.VCML.writeVCML = false;
    output_defaults.VCML.input_filename = '';
    output_defaults.VCML.downsampling = 1;
    output_defaults.VCML.add_translocation_intermediates = true;
    output_defaults.VCML.objects_always_present = true;
    output_defaults.VCML.translations = cell(0, 2);
    output_defaults.VCML.default_diffusion_coefficient = unit_convert('um2.s-1', 'm2.s-1', default_diffusion_coefficient);
    output_defaults.VCML.num_simulations = 1;
    output_defaults.VCML.delete_input_simulations = false;
    output_defaults.VCML.end_time = 20;
    output_defaults.VCML.default_time_step = 0.05;
    output_defaults.VCML.min_time_step = 0.0;
    output_defaults.VCML.max_time_step = 0.1;
    output_defaults.VCML.output_time_step = 0.5;
    output_defaults.VCML.absolute_tolerance = 1e-9;
    output_defaults.VCML.relative_tolerance = 1e-7;
    output_defaults.MCellMDL = output_defaults.VCML;
    output_defaults.MCellMDL = rmfield(output_defaults.MCellMDL, {'writeVCML', 'min_time_step'});
    output_defaults.MCellMDL.writeMCellMDL = false;
    output_defaults.MCellMDL.interaction_radius = nan;
    output_defaults.MCellMDL.max_real_time = '1:0:0:0';
    output_defaults.MCellMDL.use_reaction_rate_hack = false;
    output_defaults.MCellMDL.input_filename_pattern = '';
    output_defaults.NET.filename = '';
    output_defaults.NET.units = struct();
    output_defaults.NET.units.concentration = 'uM';
    output_defaults.NET.units.length = 'um';
    output_defaults.NET.units.time = 's';
    output_defaults.NET.translations = cell(0, 2);
    % From Mitra et al. 2004, "Modulation of the bilayer thickness of exocytic pathway membranes by membrane proteins rather than cholesterol," using values from abstract
    output_defaults.NET.effectiveWidth = unit_convert('angstrom', 'm', (37.5+39.5+35.6+42.5) / 4);
    output_defaults.NET.effective_width = unit_convert('angstrom', 'm', (37.5+39.5+35.6+42.5) / 4);
    % output_defaults.NET.use_image_adjacency = true;
    output_defaults.NET.downsampling = 1;
    output_defaults.NET.add_translocation_intermediates = true;
    % From Luby-Phelps 2000, "Cytoarchitecture and physical properties of cytoplasm: volume, viscosity, diffusion, intracellular surface area," Table I, column "Cytoplasmic D (�m2/sec)," excluding entries containing only upper bounds and taking the average of the bounds of entries containing ranges:
    output_defaults.NET.default_diffusion_coefficient = unit_convert('um2.s-1', 'm2.s-1', (0.9 + 6.9 + 43 + 27 + 5.9 + 1.7 + 6.8 + (7+11)/2 + 13.5 + (6+11)/2 + 6.7 + 1.6) / 12);
    output_defaults.OMETIFF = false;

    if ~isfield( param, 'output' )
        param.output = output_defaults;
    else
        param.output = ml_initparam( param.output, output_defaults );
    end

    %T cell model options
    if isfield( param, 'model' ) && isfield( param.model, 'tcell' )
        if isfield( param.model.tcell, 'use_two_point_synapses' )
            param.use_two_point_synapses = param.model.tcell.use_two_point_synapses;
        end
        if isfield( param.model.tcell, 'timepoints_to_include' )
            param.timepoints_to_include = param.model.tcell.timepoints_to_include;
        end
        if isfield( param.model.tcell, 'conditions_to_include' )
            param.conditions_to_include = param.model.tcell.conditions_to_include;
        end
        if isfield( param.model.tcell, 'model_type_to_include' )
            param.model_type_to_include = param.model.tcell.model_type_to_include;
        end
        param.model = rmfield( param.model, 'tcell' );
    end

    %PCA model options
    if isfield( param, 'model' ) && isfield( param.model, 'pca' )
        if isfield( param.model.pca, 'pca_synthesis_method' )
            param.pca_synthesis_method = param.model.pca.pca_synthesis_method;
        end

        % param.model = rmfield( param.model, 'pca' );
    end

    %SPHARM model options
    if isfield( param, 'model' ) && isfield(param.model, 'spharm_rpdm')
        if isfield(param.model.spharm_rpdm, 'synthesis_method')
            param.spharm_rpdm.synthesis_method = param.model.spharm_rpdm.synthesis_method;
        end
        % xruan 05/09/2019
        if isfield(param.model.spharm_rpdm, 'synthesis_resolution')
            param.spharm_rpdm.synthesis_resolution = param.model.spharm_rpdm.synthesis_resolution;
        end
        if isfield(param.model.spharm_rpdm, 'imageSize')
            param.imageSize = param.model.spharm_rpdm.imageSize;
        end
        param.model = rmfield( param.model, 'spharm_rpdm' );
    end
end%synthesis

function param = training( param )
    param = ml_initparam( param, struct( ...
        'temporary_results', [ pwd filesep 'temp' ] ));

    %Python3
    if ~isfield( param, 'python_path' )
        python_path = '/usr/local/bin/python3';
        if exist( python_path )
            param = ml_initparam( param, struct( ...
                'python_path', python_path ));
        end
    end

    %these are optional parameters that are generally hidden from users
    param = ml_initparam( param, struct( ...
        'targetDirectory', pwd, ...
        'save_helper_figures', false, ...
        'save_helper_movies', false, ...
        'interactive', false, ...
            'time_to_pause', 0.001, ...
            'skip_preprocessing', false));

    param = ml_initparam(param, ...
        struct('saveIntermediateResults', true ) );

    param.model = ml_initparam(param.model, ...
        struct('filename', 'model.xml' ) );

        param = ml_initparam(param, ...
            struct('cytonuclearflag', 'cyto' ) );

param = ml_initparam(param, ...
    struct('save_segmentations', false ) );

    %icaoberg 2018/10/26 setup default values for param.train.flag
    if ~isfield( param, 'train' ) || ~isfield( param.train, 'flag' )
        disp( ['Model training flag not set. Defaulting to ''tcell''.' ]);
        param.train.flag = 'nuclear';
    else
        if ~( strcmpi( param.train.flag, 'nuclear' ) || ...
                strcmpi( param.train.flag, 'framework' ) || ...
                strcmpi( param.train.flag, 'cell' ) || ...
                strcmpi( param.train.flag, 'protein' ) || ...
                strcmpi( param.train.flag, 'tcell' ) || ...
                strcmpi( param.train.flag, 'all' ) )
            disp( ['Unrecognized training: ' param.train.flag '. Defaulting to ''framework''.']);
            param.train.flag = 'all';
        end
    end

    if ismember( param.train.flag, {'protein', 'all'} )
        if ~(isfield( param, 'protein' ) && isfield( param.protein, 'id' ))
            param.protein.id = uuidgen();
        end

        if ~(isfield( param, 'protein' ) && isfield( param.protein, 'name' ))
            param.protein.name = 'unset';
        end
    end

    if ismember( param.train.flag, {'cell', 'framework', 'all'} )
        if ~(isfield( param, 'cell' ) && isfield( param.cell, 'id' ))
            param.cell.id = uuidgen();
        end

        if ~(isfield( param, 'cell' ) && isfield( param.cell, 'name' ))
            param.cell.name = 'unset';
        end
    end

    if ismember( param.train.flag, {'nucleus', 'framework', 'all'} )
        if ~(isfield( param, 'nucleus' ) && isfield( param.nucleus, 'id' ))
            param.nucleus.id = uuidgen();
        end

        if ~(isfield( param, 'nucleus' ) && isfield( param.nucleus, 'name' ))
            param.nucleus.name = 'unset';
        end
    end

    %T-Cell model options
    if isfield( param.model, 'tcell' )
        if isfield( param.model.tcell, 'synapse_location' )
            param.synapse_location = param.model.tcell.synapse_location;
        end
        if isfield( param.model.tcell, 'named_options_set' )
            param.named_options_set = param.model.tcell.named_options_set;
        end
        if isfield( param.model.tcell, 'use_two_point_synapses' )
            param.use_two_point_synapses = param.model.tcell.use_two_point_synapses;
        end
        if isfield( param.model.tcell, 'sensor' )
            param.sensor = param.model.tcell.sensor;
        end
        if isfield( param.model.tcell, 'timepoints_to_include' )
            param.timepoints_to_include = param.model.tcell.timepoints_to_include;
        end
        if isfield( param.model.tcell, 'model_type_to_include' )
            param.model_type_to_include = param.model.tcell.model_type_to_include;
        end
        if isfield( param.model.tcell, 'ometiff' )
            param.tcell.ometiff = param.model.tcell.ometiff;
        end
        if isfield( param.model.tcell, 'adjust_one_point_alignment' )
            param.adjust_one_point_alignment = param.model.tcell.adjust_one_point_alignment;
        end
        if isfield( param.model.tcell, 'infer_synapses' )
            param.infer_synapses = param.model.tcell.infer_synapses;
        end
        param.model = rmfield( param.model, 'tcell' );
    end

    if (isfield( param, 'protein' ) && ...
            isfield( param.protein, 'class' ) && ...
            strcmpi( param.protein.class, 'standardized_voxels' )) && ...
            (isfield( param, 'protein' ) && ...
            isfield( param.protein, 'type' ) && ...
            strcmpi( param.protein.type, 'standardized_map_half-ellipsoid' ))
        if isfield(param, 'tcell') && isfield(param.tcell, 'ometiff') && ...
                (param.tcell.ometiff == true)
            param.synapse_location = {'./temp_for_ometiff/annotations/*.csv'};
        end
    end

    %PCA model options
    if isfield( param.model, 'pca' )
        if isfield( param.model.pca, 'latent_dim' )
            param.latent_dim = param.model.pca.latent_dim;
        end

        if isfield( param.model.pca, 'imageSize' )
            param.imageSize = param.model.pca.imageSize;
        end
    end

    if ismember( param.train.flag, {'protein', 'all'} )
        if ~(isfield( param, 'protein' ) && isfield( param.protein, 'id' ))
            param.protein.id = uuidgen();
        end

        if ~(isfield( param, 'protein' ) && isfield( param.protein, 'name' ))
            param.protein.name = 'unset';
        end
    end

    if ismember( param.train.flag, {'cell', 'framework', 'all'} )
        if ~(isfield( param, 'cell' ) && isfield( param.cell, 'id' ))
            param.cell.id = uuidgen();
        end

        if ~(isfield( param, 'cell' ) && isfield( param.cell, 'name' ))
            param.cell.name = 'unset';
        end
    end

    if ismember( param.train.flag, {'nucleus', 'framework', 'all'} )
        if ~(isfield( param, 'nucleus' ) && isfield( param.nucleus, 'id' ))
            param.nucleus.id = uuidgen();
        end

        if ~(isfield( param, 'nucleus' ) && isfield( param.nucleus, 'name' ))
            param.nucleus.name = 'unset';
        end
    end

    %T-Cell model options
    if isfield( param.model, 'tcell' )
        if isfield( param.model.tcell, 'synapse_location' )
            param.synapse_location = param.model.tcell.synapse_location;
        end
        if isfield( param.model.tcell, 'named_options_set' )
            param.named_options_set = param.model.tcell.named_options_set;
        end
        if isfield( param.model.tcell, 'use_two_point_synapses' )
            param.use_two_point_synapses = param.model.tcell.use_two_point_synapses;
        end
        if isfield( param.model.tcell, 'sensor' )
            param.sensor = param.model.tcell.sensor;
        end
        if isfield( param.model.tcell, 'timepoints_to_include' )
            param.timepoints_to_include = param.model.tcell.timepoints_to_include;
        end
        if isfield( param.model.tcell, 'model_type_to_include' )
            param.model_type_to_include = param.model.tcell.model_type_to_include;
        end
        if isfield( param.model.tcell, 'ometiff' )
            param.tcell.ometiff = param.model.tcell.ometiff;
        end
        if isfield( param.model.tcell, 'adjust_one_point_alignment' )
            param.adjust_one_point_alignment = param.model.tcell.adjust_one_point_alignment;
        end
        if isfield( param.model.tcell, 'infer_synapses' )
            param.infer_synapses = param.model.tcell.infer_synapses;
        end
        param.model = rmfield( param.model, 'tcell' );
    end

    if (isfield( param, 'protein' ) && ...
            isfield( param.protein, 'class' ) && ...
            strcmpi( param.protein.class, 'standardized_voxels' )) && ...
            (isfield( param, 'protein' ) && ...
            isfield( param.protein, 'type' ) && ...
            strcmpi( param.protein.type, 'standardized_map_half-ellipsoid' ))
        if isfield(param, 'tcell') && isfield(param.tcell, 'ometiff') && ...
                (param.tcell.ometiff == true)
            param.synapse_location = {'./temp_for_ometiff/annotations/*.csv'};
        end
    end

    %PCA model options
    if isfield( param.model, 'pca' )
        if isfield( param.model.pca, 'latent_dim' )
            param.latent_dim = param.model.pca.latent_dim;
        end

        if isfield( param.model.pca, 'imageSize' )
            param.imageSize = param.model.pca.imageSize;
        end

        param.model = rmfield( param.model, 'pca' );
    end

    %SPHARM model options
    if isfield(param.model,'spharm_rpdm')
        if isfield(param.model.spharm_rpdm, 'alignment_method')
            param.spharm_rpdm.alignment_method = param.model.spharm_rpdm.alignment_method;
        end
        if isfield(param.model.spharm_rpdm, 'rotation_plane')
            param.spharm_rpdm.rotation_plane = param.model.spharm_rpdm.rotation_plane;
        end
        if isfield(param.model.spharm_rpdm, 'postprocess')
            param.spharm_rpdm.postprocess = param.model.spharm_rpdm.postprocess;
        end
        if isfield(param.model.spharm_rpdm, 'maxDeg')
            param.spharm_rpdm.maxDeg = param.model.spharm_rpdm.maxDeg;
        end
        if isfield(param.model.spharm_rpdm, 'components')
            param.spharm_rpdm.components = param.model.spharm_rpdm.components;
        end
        if isfield(param.model.spharm_rpdm, 'latent_dim')
            param.latent_dim = param.model.spharm_rpdm.latent_dim;
        end
        if isfield(param.model.spharm_rpdm, 'segminnucfraction')
            param.segminnucfraction = param.model.spharm_rpdm.segminnucfraction;
        end
    end
if ~isfield( param, 'labels' )
       param.labels = {};
end
end%training

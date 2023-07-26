function answer = image2SPHARMparameterization(varargin)

% July 15, 2023 R.F. Murphy add options.spharmrotate to do rotation by various methods

if isdeployed
    
    filename = is_deployed(varargin{1});
    load(filename);
    load(image_path);
       
else

    cur_image = varargin{1};
    options = varargin{2};

end

% set default options
if ~isfield(options, 'NMfirsttry_maxiter')
    options.NMfirsttry_maxiter = 300;
end
if ~isfield(options, 'NMretry_maxiter')
    options.NMretry_maxiter = 100;
end
if ~isfield(options, 'NMretry_maxiterbig')
    options.NMretry_maxiterbig = 300;
end
if ~isfield(options, 'NMcost_tol')
    options.NMcost_tol = 1e-7;
end
if ~isfield(options, 'NMlarge_tol') 
    options.NMlarge_tol = 1e-7;
end
if ~isfield(options, 'maxDeg')
    options.maxDeg = 31;
end
if ~isfield(options, 'hd_thresh')
    options.hd_thresh = 10;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cur_image = ml_readimage(image_path);
% disp(size(cur_image));
% disp(cur_image);
% cur_image = loadImage(cur_image, options.downsampling);
param_output = spharm_rpdm_image_parameterization(cur_image, options);

% rotationally align if desired
if options.spharmrotate
    param_post = spherical_parameterization_postprocess(param_output.vertices, param_output.faces, ...
        param_output.sph_verts, param_output.fvec, param_output.final_hd, param_output.jaccard_index, options);
    param_output.fvec = param_post.fvec;
    param_output.vertices = param_post.vertices;
    param_output.faces = param_post.faces;
    param_output.sph_verts = param_post.sph_verts;
%    if isfield('param_post','R')
    try
        param_output.R = param_post.R;
        param_output.rotation_center = param_post.rotation_center;
    catch
    end
end

%save output if deployed
% if isdeployed
disp('saving parameterizations...');
% output_dir = join([options.output_dir, '/param_output.mat']);
save(options.output_filepath, 'param_output');

if exist( [options.output_filepath(1:end-3) 'mat'], 'file' )
    answer = true;
else
    answer = false;
end


if isdeployed
    close all
end

end
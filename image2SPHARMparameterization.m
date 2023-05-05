function answer = image2SPHARMparameterization(varargin)

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
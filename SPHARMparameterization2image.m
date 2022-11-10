function answer = SPHARMparameterization2image(varargin)

if isdeployed
    
    filename = is_deployed(varargin{1});
    load(filename);
    load(model_path);
    

else
    
    param_output = varargin{1};
    options = varargin{2};
    
end

%set default options
if ~isfield(options, 'cropping')
    options.cropping = 'tight';
end
if ~isfield(options, 'oversampling_scale')
    options.oversampling_scale = 1;
end
if ~isfield(options, 'debug')
    options.debug = false;
end
%%%%%%%%%%%%%%%%%%%%%

deg = param_output.deg;
fvec = param_output.fvec;
img = spharm2image(deg, fvec, options);

disp('saving image');
save(options.output_filepath, 'img');

if exist( [options.output_filepath(1:end-3) 'mat'], 'file' )
    answer = true;
else
    answer = false;
end


end
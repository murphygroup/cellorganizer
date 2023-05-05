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
if ~isfield(options, 'figtitle')
    options.figtitle = [];
end
if ~isfield(options, 'plot')
    options.plot = 0;
end
if ~isfield(options, 'dpi')
    options.dpi = 150;
end
if ~isfield(options, 'filename')
    options.filename = [];
end
if ~isfield(options,'meshtype.type')
    options.meshtype.type = 'triangular';
end
if ~isfield(options,'meshtype.nPhi')
    options.meshtype.nPhi = 64;
end
if ~isfield(options,'meshtype.nTheta')
    options.meshtype.nTheta = 32;
end
if ~isfield(options,'nVertices')
    options.nVertices = 4002;
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

if isdeployed
    close all
end

end
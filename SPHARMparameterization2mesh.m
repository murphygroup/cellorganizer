function answer = spharmparameterization2mesh(varargin)

% Oct. 11, 2022 R.F.Murphy fix argument handling for non-deployed
if isdeployed

    filename_deployed = is_deployed(varargin{1});
    load(filename_deployed);
    
    % this reads in param_output
    load(model_path);
    
else
    
    param_output = varargin{1};
    options = varargin{2};
    
end
    
deg = param_output.deg;
fvec = param_output.fvec;
meshtype = options.meshtype;
plot = options.plot;
figtitle = options.figtitle;
filename = options.filename;
dpi = options.dpi;
[Zvert, fs] = spharm2meshfigure(deg,fvec,meshtype,plot,figtitle,filename,dpi);
mesh_out = [];
mesh_out.Zvert = Zvert;
mesh_out.fs = fs;

disp('saving mesh...');
% output_dir = join([options.output_dir, '/mesh_output.mat']);
save(options.output_filepath, 'mesh_out');

if exist( [options.output_filepath(1:end-3) 'mat'], 'file' )
    answer = true;
else
    answer = false;
end

end
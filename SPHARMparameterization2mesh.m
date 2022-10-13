function answer = spharmparameterization2mesh(varargin)

% Oct. 11, 2022 R.F.Murphy fix argument handling for non-deployed
if isdeployed
    disp('Running deployed version of reconstruct_spharm_descriptor_to_mesh...');

    %getting info read into matlab
    %when method is deployed
    
    text_file = varargin{1};

    [filepath, name, ext] = fileparts(text_file);

    if ~exist(text_file, 'file')
        warning('Input file does not exist. Exiting method.');
        return
    end

    disp(['Attempting to read input file ' text_file]);
    fid = fopen(text_file, 'r' );

    disp('Evaluating lines from input file');
    while ~feof(fid)
        line = fgets(fid);
        disp(line);
        try
            eval(line);
        catch err
            disp('Unable to parse line');
            getReport(err)
            return
        end
    end
    fclose(fid);
    if ~exist('options', 'var')
        options = {};
    end
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
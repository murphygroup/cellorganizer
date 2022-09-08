function answer = spharmparameterization2mesh(varargin)

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

    model = load(model_path);
    descriptor = model.cellShapeModel.all_spharm_descriptors;
    components = model.cellShapeModel.components;

    cellind = options.cellind;

    if numel(model.components) > 1
        descriptors = reshape(descriptors, [], 2, size(descriptors, 2), size(descriptors, 3));
        descriptors = permute(descriptors, [1, 3, 2, 4]);
        descriptors = descriptors(: , :, :, cellind);
    end

    descriptors = descriptors(: , :, cellind);
    
else
    
    descriptor = varargin{1};
    components = varargin{2};
    
end
    
mesh_out = reconstruct_spharm_descriptor_to_mesh(descriptor, components);

disp('saving mesh...');
% output_dir = join([options.output_dir, '/mesh_output.mat']);
save(options.output_filepath, 'mesh_out');

if exist( [options.output_filepath(1:end-3) 'mat'], 'file' )
    answer = true;
else
    answer = false;
end

end
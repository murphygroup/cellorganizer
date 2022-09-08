function answer = spharmparameterization2image(varargin)

if isdeployed
    disp('Running deployed version of spharm_rpdm_sample_or_reconstruct_images...');

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
    

else
    
    model = varargin{1};
    options = varargin{2};
    
end

image_mat = spharm_rpdm_sample_or_reconstruct_images(model, options);

disp('saving image(s)');
save(options.output_filepath, 'image_mat');

if exist( [options.output_filepath(1:end-3) 'mat'], 'file' )
    answer = true;
else
    answer = false;
end


end
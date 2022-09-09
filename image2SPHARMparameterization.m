function answer = image2SPHARMparameterization(varargin)

if isdeployed
    disp('Running deployed version of image2SPHARMparameterization');
    %getting info read into matlab
    %when method is deployed
    
    %init path to mat file
    image_path = '';

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

    load(image_path);
    

    %read in image
%         img = ml_readimage(cur_image);
%         cur_image = loadImage(img, options.downsampling);
       
else

    cur_image = varargin{1};
    options = varargin{2};

end

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

end
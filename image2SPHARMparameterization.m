function answer = image2SPHARMparameterization(varargin)

if isdeployed
    
    is_deployed(varargin)
    load(image_path);
       
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
function answer = SPHARMparameterization2image(varargin)

if isdeployed
    
    is_deployed(varargin)

    load(model_path);
    

else
    
    param_output = varargin{1};
    options = varargin{2};
    
end

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
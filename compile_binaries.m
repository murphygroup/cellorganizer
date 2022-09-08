% Setup all the necessary paths.
setup

% If cellorganizer-binaries directory already exists, then remove it.
if exist('cellorganizer-binaries', 'dir')
    rmdir cellorganizer-binaries s
end

% Create the cellorganizer-binaries directory
mkdir cellorganizer-binaries

% cd into cellorganizer-binaries directory
cd cellorganizer-binaries

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% img2slml %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try 
    mcc -a ../utilities/bfmatlab/bioformats_package.jar -m img2slml.m -a calc_sine_2_jacobian_sphere.m -a calc_sine_1_jacobian_sphere.m -a calc_sine_3_jacobian_sphere.m -a calc_sine_4_jacobian_sphere.m -o img2slml
catch exception
    warning('Failed to compile img2slml.m');
    msgText = getReport(exception);
    disp(msgText);
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% slml2img %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    mcc -a ../utilities/bfmatlab/bioformats_package.jar -m slml2img.m -o slml2img
catch exception
    warning('Failed to compile slml2img.m');
    msgText = getReport(exception);
    disp(msgText);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% img2shapespace %%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    mcc -a ../utilities/bfmatlab/bioformats_package.jar -m img2shapespace.m -o img2shapespace -a pad_imfunc_to_window_size.m -a diffeo_img_function.m
catch exception
    warning('Failed to compile img2shapespace.m');
    msgText = getReport(exception);
    disp(msgText);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% slml2info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    mcc -m slml2info.m -o slml2info -a pad_imfunc_to_window_size.m -a diffeo_img_function.m
catch exception
    warning('Failed to compile slml2info.m');
    msgText = getReport(exception);
    disp(msgText);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% slml2report %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    mcc -m slml2report.m -o slml2report
catch exception
    warning('Failed to compile slml2report.m');
    msgText = getReport(exception);
    disp(msgText);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% slml2slml %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try 
    mcc -m slml2slml.m -o slml2slml
catch exception
    warning('Failed to compile slml2slml.m');
    msgText = getReport(exception);
    disp(msgText);
end

%%%%%%%%% image2SPHARMparameterization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    mcc -a ../utilities/bfmatlab/bioformats_package.jar -m image2spharmparameterization.m -a calc_sine_2_jacobian_sphere.m -a calc_sine_1_jacobian_sphere.m -a calc_sine_3_jacobian_sphere.m -a calc_sine_4_jacobian_sphere.m -o image2spharmparameterization
catch exception
    warning('Failed to compile image2SPHARMparameterization.m');
    msgText = getReport(exception);
    disp(msgText);
end

%%%%%%%%% SPHARMparameterization2image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    mcc -a ../utilities/bfmatlab/bioformats_package.jar -m spharmparameterization2image.m -o spharmparameterization2image
catch exception
    warning('Failed to compile spharmparameterization2image.m');
    msgText = getReport(exception);
    disp(msgText);
end

%%%%%%%%% SPHARMparameterization2mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    mcc -m spharmparameterization2mesh.m -o spharmparameterization2mesh
catch exception
    warning('Failed to compile spharmparameterization2mesh.m');
    msgText = getReport(exception);
    disp(msgText);
end


% Remove all unnecessary files.
delete mccExcludedFiles.log
delete readme.txt
delete requiredMCRProducts.txt
delete run_img2shapespace.sh
delete run_img2slml.sh
delete run_slml2img.sh
delete run_slml2slml.sh
delete run_slml2info.sh
delete run_slml2report.sh
delete run_image2spharmparameterization.sh
delete run_spharmparameterization2image.sh
delete run_spharmparameterization2mesh.sh

% Return to previous working directory
cd ..

% Tar and compress cellorganizer-binaries directory
% Remove cellorganizer-binaries directory
try
    tar('cellorganizer-binaries.tgz', 'cellorganizer-binaries');
    rmdir cellorganizer-binaries s
catch
    warning('Failed to tar cellorganizer-binaries directory');
end

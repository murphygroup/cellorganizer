%set options
options = [];
options.NMfirsttry_maxiter = 300;
options.NMretry_maxiter = 100;
options.NMretry_maxiterbig = 300;

% decreasing these numbers decreases compute time but potential reduces model quality
options.NMcost_tol = 1e-7;
options.NMlargr_tol = 1e-7;

% degree of spherical harmonic descriptor
options.maxDeg = 31;
% if the error in the parameterization for a given cell is higher than this, discard that cell
% (note that there is a separate option hd_threshold that controls which cells are in reports)
options.hd_thresh = 10

% set where the output parameterization should be saved
options.output_filepath = '/tmp/testSPHARMparameterization.mat';

%read in selected image
global CELLORGANIZER_ROOT_PATH
path = [CELLORGANIZER_ROOT_PATH '/images/HeLa/3D/processed/LAM_cell1_ch1_t1.tif'];
iinfo=imfinfo(path);
for islice=1:length(iinfo)
    img(:,:,islice) = imread(path,'Index',islice);
end
size(img)

%downsample
xyreduced = imresize(img,1/8);
for zslice = 2:2:size(xyreduced,3)
    oz = zslice/2;
    data(:,:,oz) = (xyreduced(:,:,zslice-1)+xyreduced(:,:,zslice))/2;
end
%make sure data is shaped correctly
size(data)

image2SPHARMparameterization(data, options); %return back descriptors as dict
% read in the outputs
load(options.output_filepath);
param_output

% generate reconstructed mesh from parameterization
options_mesh = [];
options_mesh.meshtype.type = 'triangular';
options_mesh.meshtype.nVertices = 4002;
options_mesh.figtitle = [];
options_mesh.filename = [];
options_mesh.plot = 0; %don't show mesh figure
options_mesh.dpi = 150;
options_mesh.output_filepath = '/tmp/SPHARMparameterization2mesh.mat';
SPHARMparameterization2mesh(param_output, options_mesh);
load(options_mesh.output_filepath);
mesh_out

% generate reconstructed image from paramterization (regenerates mesh)
options_image = options_mesh;
options_image.output_filepath = '/tmp/SPHARMparameterization2image.mat';
options_image.cropping = 'tight';
options_image.oversampling_scale = 1;
options_image.debug = true;

SPHARMparameterization2image(param_output, options_image);
abc=load(options_image.output_filepath);
recon_img = abc.img;
size(recon_img)
figure
imshow(reshape_2d(recon_img, -1), []);

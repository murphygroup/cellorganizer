function answer = test2A ( options )


model_name = "Module2A";
module_name = model_name;
% topdir = "/home/murphylab/cellorganizer/local/results"
% outputdir = topdir + "/" + module_name
% if not os.path.exists(topdir):
%     os.makedirs(topdir)
% os.chdir(topdir)
% os.system("ls")
% if not os.path.exists(outputdir):
%     os.makedirs(outputdir)
% os.chdir(outputdir)
% os.system("ls")
% if not os.path.exists(model_name):
%     os.makedirs(model_name)
% os.chdir(model_name)
% os.system("ls")
% 
% options = {}

% #set general options
options.model.name = module_name;
options.model.filename = module_name + '.mat';
options.output_filename = module_name; 
options.verbose = 0;
options.debug = 0;
options.display = 0;
options.model.id = 123;
options.downsampling = [5, 5, 1];
options.model.resolution = [0.049, 0.049, 0.2000];
options.if_skip_cell_nuclear_model = 0;

% #set options that control what kind of model is created
options.train.flag = 'framework';
options.cell.class = 'cell_membrane';
options.cell.type = 'ratio';
% # options.nucleus.class = 'nuclear_membrane';
% # options.nucleus.type = 'cylindrical_surface';
% #options.model.spharm_rpdm.components = {'cell'};
% # if we are going to include nuclear/DNA images;
options.nucleus.class = 'nuclear_membrane';
options.nucleus.type = 'cylindrical_surface';

options.masks = {};

directory = 'images/HeLa/3D/processed';
dna = {}; cellm = {}; ptions.labels = {};
for i = 1:10
    dna{i} = [directory filesep 'LAM_cell1' num2str(i-1) '_ch0_t1.tif'];
    cellm{i} = [directory filesep 'LAM_cell1' num2str(i-1) '_ch1_t1.tif'];
%     protein{i} = [directory filesep 'LAM_cell1' num2str(i) '_ch2_t1.tif'];
%     options.labels{length(options.labels)+1} = 'Nucleoli';
    options.masks{i} = [directory filesep 'LAM_cell1' num2str(i-1) '_mask_t1.tif'];
end
% 
% directory = '/home/murphylab/cellorganizer/local/images/HeLa/3D/processed/'
% dnaImagesDirectoryPath = []
% cellImagesDirectoryPath = []
% 
% file_pattern = 'LAM_cell1?_mask_t1.tif'
% for name in glob.glob(directory + file_pattern):
%     options["masks"].append(name)
% 
% file_pattern = 'LAM_cell1?_ch1_t1.tif'
% for name in glob.glob(directory + file_pattern):
%     cellImagesDirectoryPath.append(name)
%     
% file_pattern = 'LAM_cell1?_ch0_t1.tif'
% for name in glob.glob(directory + file_pattern):
%     dnaImagesDirectoryPath.append(name)
% 
% cellImagesDirectoryPath.sort()
% dnaImagesDirectoryPath.sort()
% options["masks"].sort()

answer = img2slml('3D', dna, cellm, [], options);


% #now set options for img2info
% # options = {}
         
answer = slml2info([module_name + ".mat"],options);

% HTML(filename="index.html")



end
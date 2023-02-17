options = [];
options.shape_evolution = 'null';
options.labels = 'unique';
options.subsize = 400; %controls the size of the shape (larger number is smaller)
options.includenuclear = 0;
options.includecell = 0;
options.includeprot = 1;
options.paired = 0;
model_folder = "./"
model_files_1 = 'golgi_GT_model_s.mat';
model_files_2 = 'nucleoli_GT_model_s.mat';
slml2report(model_files_1 , model_files_2, options)


options = [];
options.shape_evolution = 'null';
options.labels = 'unique';
options.subsize = 400; %controls the size of the shape (larger number is smaller)
%options.includenuclear = 0;
options.paired = 0;
model_folder = "./"
model_files_1 = 'Module2A1.mat';
model_files_2 = 'Module2A2.mat';
slml2report(model_files_1 , model_files_2, options)

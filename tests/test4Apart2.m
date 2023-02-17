options = [];
options.shape_evolution = 'none';
options.labels = 'unique';
options.subsize = 400; %controls the size of the shape (larger number is smaller)
options.viewangle = [0,90]; % down z axis
%options.viewangle = [90,0]; % side view
options.hd_threshold = 10 % filter out objects with Hausdorff distance greater than this
        
answer = slml2info({'demo3D61.mat'}, options);

clear all
close all

addpath slicfiles/
cellnums =  [1 2 4 6 10 12:42 46:49 51:52]
K = 1;

%subfolder = 1;
%subfolder = 2;
subfolder = 3;  %%

featType = 'all';
imCategory = 'general';
batchno = 1;


w = [1,1,1,1,1,1,1,0,0,0,0,0,0
     0,0,0,0,1,1,1,1,1,0,0,0,0
     0,0,0,0,1,1,1,1,1,0,1,1,1
     0,0,0,0,1,1,1,1,1,1,1,1,1
     1,1,1,1,1,1,1,1,1,1,1,1,1];

for i = 1:size(w,1)
finTable = [];
K = 1;
for cellnum = cellnums
       disp([num2str(w(i,:)),' cell ',num2str([cellnum])])
	finTable(K,:) = [recover_parameters_from_real_image_ownlib(subfolder,batchno,cellnum,w(i,:)),cellnum];
	K = K + 1;
end

if ~exist(['outputs_' num2str(subfolder) '/results/'],'dir')
   mkdir(['outputs_' num2str(subfolder) '/results/']);
end
%save(['outputs_' num2str(subfolder) '/results/result_new.mat']);
save(['outputs_' num2str(subfolder) '/results/result_new_' regexprep(num2str(w(i,:)),'  ','_') '.mat']);

end


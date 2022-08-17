function perobjmdl(savepath,resolution)
%This function saves .mdl files for simulating object meshes in MCell.
%currently it saves a single .mdl file for each endosome with only size
%data. 
%
%Inputs: 
%savepath = string pointing to the desired folder to save the .mdl files
%
%Outputs:
%set of .mdl files 
%
%Author: Devin Sullivan July 23, 2013
mkdir(savepath)

if isempty(resolution)
    resolution = [1,1,1];
end

%first load temp results 
load([pwd filesep 'temp' filesep 'primitives']);
load([pwd filesep 'temp' filesep 'OriginalObjects.mat']);
%         'GaussObjects','pos');
% objsizevec = objsizevec.*repmat(resolution,[size(objsizevec,1),1]);


%Create an image for each object
for i = 1:size(objsizevec,1)
    %initialize image to zeros
    img = zeros(ceil(objsizevec(i,1)),ceil(objsizevec(i,2)),ceil(objsizevec(i,3)));
    %put objects into image
    protimg = ml_imaddobj2(blankimg,GaussObjects,...
        struct('method','replace','pos',newpos,'objectmethod',param.sampling.method));
    
    

end



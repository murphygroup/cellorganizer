function FV = makeIso2mesh(volimage,param)

if nargin<2
    param.display = 1;
end

if ischar(volimage)
    try
        volimage = ml_readimage(volimage);
    catch
        %It is likely a "preSBML.mat" file
        load(volimage);
        volimage = frameworkSBML.list(2).img;
    end
end

% volimage is a volumetric image such as an X-ray or MRI image
% A,b are registration matrix and vector, respectively
%% perform mesh generation

%This was the default values used in the example file from iso2mesh
% [node,elem,face]=vol2mesh(volimage>0.05,1:size(volimage,1),1:size(volimage,2),...
%     1:size(volimage,3),2,2,1);
%Modified numbers 
[node,elem,face]=vol2mesh(volimage>0.05,1:size(volimage,1),1:size(volimage,2),...
    1:size(volimage,3),10,10,1);

%% alternatively, one can use the following cmd as a less robust approach
% [node,elem,face]=vol2mesh(volimage>0.05,1:size(volimage,1),1:size(volimage,2),...
%                           1:size(volimage,3),0.2,2,1,'simplify');
FV.vertices = node;
FV.vertices(:,1:2) = node(:,2:-1:1);
FV.faces = face;
%% visualize the resulting mesh
if param.display
    plotmesh(node,face);
    axis equal;
end
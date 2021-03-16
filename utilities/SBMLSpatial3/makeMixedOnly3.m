function makeMixedOnly3

%This function makes the SBML-spatial standard 'MixedOnly.xml' file. 
%
%The function uses pre-made matlab structures from CellOrganizer. To
%completely re-create the instance please run demo3DSBML_Standards
%
%Created by: Devin P Sullivan 9/5/14


if nargin<1
    preSBML = 'preSBML.mat';
end

load(preSBML)

% simpleCube = zeros(size(img));
simpleCube = zeros(10,10,10);
%simpleCube(100:size(img,1)-100,100:size(img,2)-100,:) = 1;
border = 2;
simpleCube(border:size(simpleCube,1)-border,border:size(simpleCube,2)-border,:) = 1;

simpleCubeSBML = createSBMLFrameworkstruct3({simpleCube},param);



outdir = [pwd, filesep, 'MixedOnly.xml'];

instance2SBML_mod(primitives,simpleCubeSBML,outdir,[],model.proteinModel.resolution)
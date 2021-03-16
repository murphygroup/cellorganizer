function makeCSGOnly(preSBML)
%This function makes the SBML-spatial standard 'CSGOnly.xml' file. 
%
%The function uses pre-made matlab structures from CellOrganizer. To
%completely re-create the instance please run demo3DSBML_Standards
%
%Created by: Devin P Sullivan 9/5/14

if nargin<1
    preSBML = 'preSBML.mat';
end

load(preSBML)

outdir = [pwd, filesep, 'CSGOnly'];

instance2SBML_mod(primitives,[],[outdir,'.xml'],[],model.proteinModel.resolution)
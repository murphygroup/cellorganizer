function loadAllRegionFeaturesSimple(readdir,rootmeta,rootdir,outfilename)

% readdir is the directory storing features.
% rootmeta is the directory storing metadata for HPA images.
% rootdir is the root directory storing all images.

% Copyright (C) 2010  Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

% 10 Jan 10 - Jieyue Li
% 2011-11-12 tebuck: copied from loadAllRegionFeatures.m,
% modifying to use no metadata.


if ~exist('readdir','var')
    readdir = './data/features/region/';
end

if ~exist('rootdir','var')
    rootdir = '/home/jieyuel/images/HPA/images/IFconfocal/';
end

tmpidx = strfind(imagelist,rootdir); %%
tmpidx = cellfun(@isempty,tmpidx,'UniformOutput',false);
tmpidx = cell2mat(tmpidx);
tmpidx = ~tmpidx;
imagelist = imagelist(tmpidx);
antibodyids = antibodyids(tmpidx);
classlabels = classlabels(tmpidx);
cellabels = cellabels(tmpidx);
specificity = specificity(tmpidx);
staining = staining(tmpidx);

staining = cellfun(@(x) x(x(2)),staining,'UniformOutput',false); %%
staining = cell2mat(staining); %%

ulist = listunixformat( imagelist);


featsets = {...
    'mutualInfo', ...
    'mutualInfox2', ...
    'nonObjFluor', ...
    'nonObjFluorx2', ...
    'obj', ...
    'objx2', ...
    'overlap', ...
    'overlapx2', ...
    'texture', ...
    'texturex2', ...
    'texturex4', ...
    'tas', ...
    'tasx2', ...
    'objRegion',...
    'objRegionx2',...
  };

numfeats = [30, 30, 4, 4, 18, 18, 39, 39, 156, 156, 156, 18, 18, 14, 14];


datasettype = 'Region';


encount = cumsum(numfeats);
stcount = [1 encount(1:end-1)+1];

for i=1:length(featsets)
    tmplist = findreplacestring( ulist, '/','_');
    rootdirtmp = findreplacestring( rootdir, '/','_');
    tmplist = findreplacestring( tmplist, rootdirtmp, readdir);
    flists{i} = findreplacestring( tmplist, '.tif', ['_' featsets{i} '.mat']);
    featlists{i} = listmatlabformat( flists{i})';
end
%keyboard

allfeatures = zeros(length(imagelist)*20,encount(end));
idx = zeros(length(imagelist)*20,1);
allclasslabels = [];
allclasslabels{length(imagelist)*20} = [];
allantibodyids = idx;
allcellids = idx;
allspecificity = idx;
allstaining = idx;
allimagelist = allclasslabels;
counter = 1;
fiddle=fopen('regionload.txt','w');
fwrite(fiddle, num2str(length(imagelist)));
fwrite(fiddle,char(10));

disp( length(imagelist));
for i=1:length(imagelist)
    allnames = [];
    breaker = 0;
    features = zeros(1,encount(end));
    for j=1:length(featsets)
        if ~exist(featlists{j}{i},'file')
            breaker = 1;
            continue;
        end

        errch=0;
        eval('load(featlists{j}{i});','disp(i);disp(featlists{j}{i});delete(featlists{j}{i});errch=1;');
        if errch
            continue;
        end

        load(featlists{j}{i});
        allnames = [allnames names];
        if size(feats,1)~=size(features,1)
            features = repmat(features, [size(feats,1) 1]);
        end
        features(:,stcount(j):encount(j)) = feats;
    end
    
    if ~breaker
        idx(counter:counter+size(features,1)-1) = 1;
        allfeatures(counter:counter+size(features,1)-1,:) = features;
        allclasslabels(counter:counter+size(features,1)-1) = classlabels(i);
        allantibodyids(counter:counter+size(features,1)-1) = antibodyids(i);
        allcellids(counter:counter+size(features,1)-1) = cellabels(i);
        allimagelist(counter:counter+size(features,1)-1) = imagelist(i);
        allspecificity(counter:counter+size(features,1)-1) = specificity(i);
        allstaining(counter:counter+size(features,1)-1) = staining(i);
        counter = counter+size(features,1);
    end
fwrite(fiddle, [ num2str(i) '  ']);
end
char(10);

fwrite(fiddle,char(10));
fclose(fiddle);

names = allnames;

ind = find(idx);
allfeatures = allfeatures(ind,:);
imagelist = allimagelist(ind);
classlabels = allclasslabels(ind);
antibodyids = allantibodyids(ind);
cellabels = allcellids(ind);
specificity = allspecificity(ind);
staining = allstaining(ind);


[classes Icl classlabels] = unique(classlabels);
[proteins Ipr protlabels] = unique(antibodyids);

if ~exist('outfilename','var')
   outfilename = './data/features/regionfeatures_all.mat';
end
save( outfilename,...
    'allfeatures','classlabels','antibodyids','imagelist','specificity',...
    'classes','staining','cellabels','names');

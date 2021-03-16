function [features_filename, feature_computation_time] = loadAllRegionFeatures(readdir,rootmeta,rootdir,outfilename, featsets, old_feature_set_behavior)
%function [features_filename] = loadAllRegionFeatures(readdir,rootmeta,rootdir,outfilename, featsets, old_feature_set_behavior)
%%function [features_filename] = loadAllRegionFeatures(readdir,rootmeta,rootdir,outfilename, channel_as_protein, featsets, old_feature_set_behavior)

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
% 2011-11-15 tebuck: adding slfnames to files.
% 2011-11-16 tebuck: channel_as_protein sets the protein image to
% that channel's image and blanks the channel's image.
% 2011-11-18 tebuck: made featsets an arg. INCOMPLETE, NEEDS
% ACCOMODATIONS LIKE WITH numfeats. Depends on DataHash.
% 2011-12-30 tebuck: channel_as_protein removed here as well.



if ~exist('readdir','var')
    readdir = './data/features/region/';
end

if ~exist('rootmeta','var')
    rootmeta = './images/';
end

if ~exist('rootdir','var')
    rootdir = '/home/jieyuel/images/HPA/images/IFconfocal/';
end

if nargin > 1
  feature_computation_time = 0.; 
end


% $$$ if ~exist('channel_as_protein','var')
% $$$     channel_as_protein = 'protein';
% $$$ end


load([rootmeta,'hpalistsall.mat']);

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



original_featsets = {...
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
original_numfeats = [30, 30, 4, 4, 18, 18, 39, 39, 156, 156, 156, 18, 18, 14, 14];
% tebuck: custom features:
all_featsets = original_featsets; 
all_numfeats = original_numfeats; 
all_featsets = [all_featsets, {'nucStats'}]; 
all_numfeats = [all_numfeats, 20]; 
all_featsets = [all_featsets, {'protTotalIntensity'}]; 
all_numfeats = [all_numfeats, 1]; 
all_featsets = [all_featsets, {'protMeanIntensity'}]; 
all_numfeats = [all_numfeats, 1]; 

% $$$ if ~exist('featsets','var')
% $$$   %featsets = original_featsets;
% $$$   old_feature_set_behavior = true;
% $$$ end


% Build a dictionary of number of features per feature set:
all_numfeats_cell = num2cell(all_numfeats);
feature_length_dictionary = struct();
for field_index = 1:length(all_featsets)
  feature_length_dictionary.(all_featsets{field_index}) = all_numfeats(field_index);
end
%feature_length_dictionary


feature_set_hash = ''; 
% $$$ if exist('old_feature_set_behavior', 'var')
if ~exist('featsets','var') || isempty(featsets)
  numfeats = original_numfeats;
  featsets = original_featsets;
else
  %numfeats = getfield(feature_length_dictionary,
  %original_featsets, {':'})
  featsets = sort(featsets);
  numfeats = zeros(1, length(featsets));
  % Use a hash of the feature sets to avoid overwriting files with
  % high probability:
  %feature_set_hash = DataHash(featsets, struct('Method',
  %'SHA-512'));
  % 2^(-4 * 40) ~ 6.84e-49 prior probability of collision:
  feature_set_hash = DataHash(featsets, struct('Method', 'SHA-1'));
  for field_index = 1:length(featsets)
    if ~any(strcmp(featsets{field_index}, all_featsets))
      error('%s not a valid feature set', featsets{field_index})
    end
    numfeats(field_index) = feature_length_dictionary.(featsets{field_index});
  end
end
%numfeats

datasettype = 'Region';


encount = cumsum(numfeats);
stcount = [1 encount(1:end-1)+1];

for i=1:length(featsets)
    tmplist = findreplacestring( ulist, '/','_');
    rootdirtmp = findreplacestring( rootdir, '/','_');
    tmplist = findreplacestring( tmplist, rootdirtmp, readdir);
% $$$     flists{i} = findreplacestring( tmplist, '.tif', ['_' featsets{i} '.mat']);
    feature_set_suffix = ['_', featsets{i}];
% $$$     if ~strcmpi(channel_as_protein, 'protein')
% $$$       feature_set_suffix = [feature_set_suffix, '_', channel_as_protein, '-focus'];
% $$$     end
    flists{i} = findreplacestring( tmplist, '.tif', [feature_set_suffix, '.mat']);
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
    allslfnames = [];
    breaker = 0;
    features = zeros(1,encount(end));
    for j=1:length(featsets)
        if ~exist(featlists{j}{i},'file')
            check_for_file = featlists{j}{i}
            breaker = 1
            continue;
        end

        errch=0;
        eval('load(featlists{j}{i});','disp(i);disp(featlists{j}{i});delete(featlists{j}{i});errch=1;');
        if errch
            continue;
        end

        load(featlists{j}{i});
        allnames = [allnames names];
        allslfnames = [allslfnames slfnames];
        %fprintf('featsets{j} = %s, size(feats) = [%d, %d], all_numfeats[j] = %d\n', featsets{j}, size(feats), any(all_numfeats == size(feats,2)))
        %fprintf('featsets{j} = %s, size(feats) = [%d, %d], all_numfeats(j) = %d\n', featsets{j}, size(feats), all_numfeats(j))
        %fprintf('featsets{j} = %s, size(feats) = [%d, %d], numfeats should be %d\n', featsets{j}, size(feats), feature_length_dictionary.(featsets{j}))
        if size(feats, 2) ~= feature_length_dictionary.(featsets{j})
          %fprintf('featsets{j} = %s, size(feats) = [%d, %d], numfeats should be %d\n', featsets{j}, size(feats), feature_length_dictionary.(featsets{j}))
          warning('featsets{j} = %s, size(feats) = [%d, %d], numfeats should be %d\n', featsets{j}, size(feats), feature_length_dictionary.(featsets{j}))
        end
        
        if size(feats,1)~=size(features,1)
            features = repmat(features, [size(feats,1) 1]);
        end
        features(:,stcount(j):encount(j)) = feats;
        if exist('feature_computation_time', 'var') && exist('computation_time', 'var')
          feature_computation_time = feature_computation_time + computation_time; 
        end
        
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
slfnames = allslfnames;

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
   outfilename = './data/features/regionfeatures_all';
   outfilename = [outfilename, '.mat'];
end
if ~isempty(feature_set_hash)
  outfilename = [outfilename(1:end - 4), '_', feature_set_hash, outfilename(end - 3:end)];
end
save( outfilename,...
    'allfeatures','classlabels','antibodyids','imagelist','specificity',...
    'classes','staining','cellabels','names','slfnames');

if nargout > 0
  features_filename = outfilename;
end

% function loadAllFieldFeatures
clear
addpath lib/system

load data/metadata/hpalistsall.mat

ulist = listunixformat( imagelist);


featsets = {...
    'mutualInfo16', ...
    'nonObjFluor', ...
    'obj', ...
    'overlap', ...
    'nonObjFluorMulti', ...
    'objMulti', ...
    'overlapMulti', ...
    'tas', ...
    'texture64', ...
    'texture64x2', ...
    'texture64x4', ...
  };

numfeats = [6,  2,  9, 24,  2, 9, 24, ...
    18, 26, 26, 26];



datasettype = 'Field';
for i=1:length(featsets)
    featsets{i} = [featsets{i} datasettype];
end



encount = cumsum(numfeats);
stcount = [1 encount(1:end-1)+1];

for i=1:length(featsets)
    tmplist = findreplacestring( ulist, '/','_');
    tmplist = findreplacestring( tmplist, '_images_HPA_images_IFconfocal_', ['data/features/' lower(datasettype) '/']);
    flists{i} = findreplacestring( tmplist, '.tif', ['_' featsets{i} '.mat']);
    featlists{i} = listmatlabformat( flists{i})';
end


allfeatures = zeros(length(imagelist),encount(end));
allslfnames = [];
allnames = [];
idx = zeros(length(imagelist),1);
fiddle=fopen([lower(datasettype) 'load2.txt'],'w');
fwrite(fiddle, num2str(length(imagelist)));
fwrite(fiddle,char(10));
for i=1:length(imagelist)
    breaker = 0;
    features = zeros(1,encount(end));
    allslfnames = [];
    allnames = [];
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
% featlists{j}{i}
        features(stcount(j):encount(j)) = feats;
        allslfnames = [allslfnames slfnames];
        allnames = [allnames names];
    end
    
    if ~breaker
        idx(i) = 1;
        allfeatures(i,:) = features;
    else
        disp(i);
    end
fwrite(fiddle, [ num2str(i) '  ']);
end
fwrite(fiddle,char(10));
fclose(fiddle);
ind = find(idx);
allfeatures = allfeatures(ind,:);
imagelist = imagelist(ind);
classlabels = classlabels(ind);
antibodyids = antibodyids(ind);
specificity = specificity(ind);
staining = staining(ind);
cellabels = cellabels(ind);

[classes Icl classlabels] = unique(classlabels);
[proteins Ipr protlabels] = unique(antibodyids);

names = allnames;
slfnames = allslfnames;
save( ['data/tmp/' lower(datasettype) 'features_all2.mat'],...
    'allfeatures','classlabels','antibodyids','imagelist','specificity',...
    'classes','staining','cellabels','names','slfnames');
    
% filtertrainfeats_script

% save( 'data/tmp/fieldfeatures_train.mat','allfeatures','classlabels','classes','imagelist','antibodyids');
% save( 'data/tmp/fieldfeatures_train.mat','allfeatures','classlabels','classes','imagelabels','antibodyids');


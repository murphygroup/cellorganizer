clear

addpath lib/classification
addpath lib/classification/libsvm-mat-2.84-1
addpath lib/results
addpath lib/system

load data/metadata/hpalistsall.mat
ulist = listunixformat( imagelist);

featsets = { 'spatialRegion'};
numfeats = 120;


encount = cumsum(numfeats);
stcount = [1 encount(1:end-1)+1];

for i=1:length(featsets)
    tmplist = findreplacestring( ulist, '/','_');
    tmplist = findreplacestring( tmplist, '_images_HPA_images_IFconfocal_', 'data/features/region/');
    flists{i} = findreplacestring( tmplist, '.tif', ['_' featsets{i} '.mat']);
    featlists{i} = listmatlabformat( flists{i})';
end

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
fprintf( [num2str(i) '  ']);
    breaker = 0;
    features = zeros(1,encount(end));
    for j=1:length(featsets)
        if ~exist(featlists{j}{i},'file')
            breaker = 1;
            continue;
        end
        load(featlists{j}{i});
        if size(allfeats,1)~=size(features,1)
            features = repmat(features, [size(allfeats,1) 1]);
        end
        features(:,stcount(j):encount(j)) = allfeats;
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

outfilename = 'data/tmp/histfeatures_all.mat';
save( outfilename,...
    'allfeatures','classlabels','antibodyids','imagelist','specificity',...
    'classes','staining','cellabels');

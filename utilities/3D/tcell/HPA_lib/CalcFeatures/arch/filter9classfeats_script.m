
% choose only 1 location pattern images
idx = find(specificity~=1);
allfeatures(idx,:) = [];
classlabels(idx) = [];
antibodyids(idx) = [];
imagelist(idx) = [];
specificity(idx) = [];
staining(idx) = [];
cellabels(idx) = [];

% consider all staining
idx = find(staining==4);
allfeatures(idx,:) = [];
classlabels(idx) = [];
antibodyids(idx) = [];
imagelist(idx) = [];
specificity(idx) = [];
staining(idx) = [];
cellabels(idx) = [];


if 1==1
% remove cytoplasm class
idx = find(strcmp( 'loc_cytoplasm',classes));
idx = find(classlabels==idx);
allfeatures(idx,:) = [];
classlabels(idx) = [];
antibodyids(idx) = [];
imagelist(idx) = [];
specificity(idx) = [];
staining(idx) = [];
cellabels(idx) = [];
end

% merge nuclear classes
idx1 = find(strcmp( 'loc_nucleus_without_nucleoli',classes));
idx2 = find(strcmp( 'loc_nucleus',classes));
classlabels(classlabels==idx1) = idx2;

[U I classlabels] = unique(classlabels);
classes = classes(U);

% remove classes with less than N representative proteins
[proteins Ipr protlabels] = unique(antibodyids);
[c b] = hist(classlabels(Ipr),1:1:max(classlabels));
N = 10;
idx = find(c<N);
for i=1:length(idx)
    ind = find(classlabels==idx(i));
    allfeatures(ind,:) = [];
    classlabels(ind) = [];
    antibodyids(ind) = [];
    imagelist(ind) = [];
    specificity(ind) = [];
    staining(ind) = [];
    cellabels(ind) = [];
end
[proteins Ipr protlabels] = unique(antibodyids);
[U I classlabels] = unique(classlabels);
classes = classes(U);

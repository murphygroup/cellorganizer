function [allfeatures, labels, filters] = filterfeats( allfeatures, labels, filters)
% function [allfeatures, labels, filters] = filterfeats( allfeatures, labels, filters)
% 
% labels and filters are structs

classlabels = labels.classlabels;
antibodyids = labels.antibodyids;
imagelist = labels.imagelist;
classes = labels.classes;

specificity = filters.specificity;
staining = filters.staining;
cellabels = filters.cellabels;
filtertype = filters.filtertype;

% FILTER ON NUMBER PATTERNS
switch filtertype
  case 'multiacross'
    % remove only 1 location pattern images
    % idx = find(specificity==1);
    idx = find(specificity>2);
  case 'multia431'
    % remove 3+ location pattern images
    idx = find(specificity>2);
  case 'multiu2os'
    % remove 3+ location pattern images
    idx = find(specificity>2);
  case 'multiu251mg'
    % remove 3+ location pattern images
    idx = find(specificity>2);
  otherwise
    % choose all non-1 location pattern images
    idx = find(specificity~=1);
end
allfeatures(idx,:) = [];
classlabels(idx) = [];
antibodyids(idx) = [];
imagelist(idx) = [];
specificity(idx) = [];
staining(idx) = [];
cellabels(idx) = [];


% FILTER ON STAINING
% consider all staining
idx = find(staining==4);
allfeatures(idx,:) = [];
classlabels(idx) = [];
antibodyids(idx) = [];
imagelist(idx) = [];
specificity(idx) = [];
staining(idx) = [];
cellabels(idx) = [];


% FILTER ON CELL TYPE
idx = [];
switch lower(filtertype)
  case 'multia431'
    idx = find(cellabels~=1);
  case 'multiu2os'
    idx = find(cellabels~=2);
  case 'multiu251mg'
    idx = find(cellabels~=3);
  case 'a431'
    idx = find(cellabels~=1);
  case 'u2os'
    idx = find(cellabels~=2);
  case 'u251mg'
    idx = find(cellabels~=3);
  otherwise
    % remove cytoplasm class
    idx = [];
end
allfeatures(idx,:) = [];
classlabels(idx) = [];
antibodyids(idx) = [];
imagelist(idx) = [];
specificity(idx) = [];
staining(idx) = [];
cellabels(idx) = [];


% FILTER ON CYTOPLASM CLASS
idx = find(strcmp( 'loc_cytoplasm',classes));
switch lower(filtertype)
  case 'cyto'
    % remove all non cytoplasm class
    idx = find(classlabels~=idx);
  otherwise
    % keep cytoplasm class
    idx = [];
    % idx = find(classlabels==idx);
end
allfeatures(idx,:) = [];
classlabels(idx) = [];
antibodyids(idx) = [];
imagelist(idx) = [];
specificity(idx) = [];
staining(idx) = [];
cellabels(idx) = [];


if 1==1
% MERGE CLASSES
% merge nuclear classes
idx1 = find(strcmp( 'loc_nucleus_without_nucleoli',classes));
idx2 = find(strcmp( 'loc_nucleus',classes));
classlabels(classlabels==idx1) = idx2;
end

[U I classlabels] = unique(classlabels);
classes = classes(U);


% REMOVE SMALL CLASSES
% remove classes with less than N representative proteins
[proteins Ipr protlabels] = unique(antibodyids+100000*cellabels);
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


% filter rows with isnan features
[idx cidx] = find(isnan(allfeatures));
allfeatures(idx,:) = [];
classlabels(idx) = [];
antibodyids(idx) = [];
imagelist(idx) = [];
specificity(idx) = [];
staining(idx) = [];
cellabels(idx) = [];


labels.classlabels = classlabels;
labels.antibodyids = antibodyids;
labels.imagelist = imagelist;
labels.classes = classes;

filters.specificity = specificity;
filters.staining = staining;
filters.cellabels = cellabels;
filters.filtertype = filtertype;


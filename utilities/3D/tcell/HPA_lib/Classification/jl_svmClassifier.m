function [predict_label accuracy weights] = jl_svmClassifier(allfeatures,modelfile,classlabels)

%classify new cells or images with SLF features using trained svm model.

if ~exist('modelfile','var')
   modelfile = './HPA_lib/trainedModel/IF_atlas_5_regionfeatures_single_A431R2_5folds_svmclassifier.mat';
end
load(modelfile,'model','idx_feat','mu','st');


%allfeatures = single(allfeatures);

[rows,cols] = find(isnan(allfeatures));
allfeatures(rows,:) = [];

%standardizing
allfeatures = (allfeatures - repmat( mu, [size(allfeatures,1) 1]))./ repmat( st, [size(allfeatures,1) 1]);

% normalizing
allfeatures = normalizeFeatures( allfeatures);

%feature selected
allfeatures = allfeatures(:,idx_feat);


%classifying
if ~exist('classlabels','var')
   classlabels = zeros(size(allfeatures,1),1);
else
   classlabels = classlabels(:);
   classlabels(rows) = [];
end

[predict_label accuracy weights] = svmpredict(classlabels, double(allfeatures), model, '-b 1'); 


function allweights = svmCrossValidate( tmpname, allfeatures, classlabels, antibodyids, idx_sda, c, g, Ninit)
% function [protlabels,predlabels,proteins,newallweights] = ...
%     svmCrossValidate( allfeatures, classlabels, antibodyids, idx_sda)

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


if isempty(idx_sda)
   idx_sda = [1:size(allfeatures,2)];
end
length(idx_sda)


if ~isempty(tmpname) && exist(['./data/classification/',tmpname,'_',num2str(Ninit),'folds.mat'],'file')
   ['./data/classification/',tmpname,'_',num2str(Ninit),'folds.mat']
   load(['./data/classification/',tmpname,'_',num2str(Ninit),'folds.mat'],'trainidx','testidx','protlabels','N');
else
   [trainidx,testidx,protlabels,N] = splitbyproteins(classlabels,antibodyids, Ninit);
end
N


U = unique(classlabels);


if ~exist('Nfolds','var')
    Nfolds = 1:N;
end


predlabels = zeros(length(classlabels),1);
allweights = zeros(length(classlabels),length(U));
committee = [];
for i=1:N
    trainind = find(ismember(protlabels,trainidx{i}));
    testind = find(ismember(protlabels,testidx{i}));

%     trainind = trainidx{i};
%     testind = testidx{i};

    traindata = allfeatures( trainind,:);
    testdata = allfeatures( testind,:);
    
    trainlabels = classlabels(trainind);
    testlabels = classlabels(testind);
    
    % zscore standardization
    [traindata, testdata, mu, st] = standardizeFeatures( traindata, testdata);

    % normalizing
    traindata = normalizeFeatures( traindata);
    testdata = normalizeFeatures( testdata);
    
    ut = clock;
    ut = [num2str(ut(end-2)) num2str(ut(end-1)) num2str(ut(end))];
    ut(ut=='.') = '_';

    if exist('idx_sda','var')
        idx_feat = idx_sda{i};
    else
        idx_feat = featureSelection( traindata, trainlabels, 'sda', ut);
    end

    idx = unique(idx_feat);

    traindata = traindata(:,idx);
    testdata = testdata(:,idx);
    
    wopt = classweights( trainlabels);
%     options = ['-s 0 -m 1000 -t 2 -c 2048 -g 4 -b 1' wopt];
%     options = ['-s 0 -t 2 -c 256 -g 16 -b 1' wopt];
%    options = ['-s 0 -t 2 -c 256 -g ',num2str(g),' -b 1' wopt];
    options = ['-s 0 -t 2 -c ',num2str(c(i)),' -g ',num2str(g(i)),' -b 1' wopt];

    model = svmtrain( trainlabels, double(traindata), options);
    [predict_label accuracy weights] = svmpredict( testlabels, double(testdata), model, '-b 1');

    committee(i).model = model;
    committee(i).idx_feat = idx;
    committee(i).options = options;
    
    [a ind] = max(weights,[],2);
    predlabels(testind) = model.Label(ind);
    allweights(testind,model.Label) = weights;
end

%[cc uu dd ww] = conmatrix( classlabels, predlabels);


% proteins = unique(antibodyids);
% predlabels = zeros(size(proteins));
% protlabels = zeros(size(proteins));
% newallwegihts = zeros(length(proteins),size(allweights,2));
% for i=1:length(proteins)
%     idx = find(proteins(i)==antibodyids);
%     weights = sum(allweights(idx,:),1);
%     [a label] = max(weights,[],2);
    
%     predlabels(i) = label;
    
%     protlabels(i) = mode(classlabels(idx));
%     newallweights(i,:) = weights/length(idx);
% end

% [cc uu dd ww] = conmatrix( protlabels,predlabels);


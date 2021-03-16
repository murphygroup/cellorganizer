function jl_CV_split(dataname,Ninit,Nclass,Nfolds)

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


if ~exist('dataname','var')
   dataname = 'regionfeatures_single_A431.mat';
end
dataname

if ~exist('Ninit','var')
   Ninit = 5;
end
Ninit

if ~exist('Nclass','var')
   Nclass = Ninit;
elseif Nclass<Ninit
   warning('Nclass cannot be smaller than Ninit and has been increased to equal to Ninit!');
   Nclass = Ninit;
end
Nclass

load(dataname)
classlabels = classlabels(:);
imagelist = imagelist(:);
classes = classes(:);

[rows,cols] = find(isnan(allfeatures));
allfeatures(rows,:) = [];
classlabels(rows,:) = [];
antibodyids(rows,:) = [];
imagelist(rows,:) = [];
specificity(rows,:) = [];
staining(rows,:) = [];

%iidx = find(isinf(allfeatures));  %%
%allfeatures(iidx) = sign(allfeatures(iidx)).*realmax('single');
[iidx,iidy] = find(isinf(allfeatures));
allfeatures(:,iidy) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%delete classes with small number of antibodies.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%class_set=unique(classlabels);
%for i=1:length(class_set)
%c(i)=length(unique(antibodyids(find(classlabels==class_set(i)))));
%end

class_del = [];
class_set = unique(classlabels);
classes_idx = 1:length(classes);
if length(unique(classlabels))<length(classes)
   warning('Some classes are not in the classlabels. Delete them... Check data!');
   tmpind = ~ismember(classes_idx,class_set); %%
   class_del = [class_del,classes_idx(tmpind)];
elseif length(unique(classlabels))>length(classes)
   warning('Some classlabels are not in the classes. Delete them... Problematic and Data-Check recommended!');
   tmpind = ~ismember(class_set,classes_idx); %%
   class_del = [class_del,class_set(tmpind)'];
end

for i=1:length(class_set)
    numc(i)=length(unique(antibodyids(find(classlabels==class_set(i)))));
    if numc(i)<Nclass
       class_del = [class_del,class_set(i)];
    end
end
tmpidx = ismember(classlabels,class_del);
allfeatures(tmpidx,:) = [];
antibodyids(tmpidx) = [];
classlabels(tmpidx) = [];
classes(class_del) = [];
class_del = sort(class_del);
for i=1:length(class_del)
    tmpdex = find(classlabels>class_del(i));
    classlabels(tmpdex) = classlabels(tmpdex)-1;
    class_del = class_del-1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%addpath('./libsvm-mat-2.88-1/');
%addpath('./lib_Justin/lib/classification');
%addpath('./rtag/randtag_scripts/SLIC/SDA/matlab');

%addpath('./lib_Justin/lib/classification');
%addpath('./lib_Justin/lib/classify');
%addpath('./lib_Justin/lib/results');


boolvar = 0;
if exist('Nfolds','var')
    boolvar = 1;
end

abids2 = antibodyids;

%Ninit = 5;
[trainidx,testidx, protlabels,N] = splitbyproteins( classlabels, abids2, Ninit);
N


if ~exist('Nfolds','var')
    Nfolds = 1:N;
end
Nfolds


traintrainidx = cell(N,1);
traintestidx = cell(N,1);
trainprotlabels = cell(N,1);
trainN = zeros(N,1);
for i=Nfolds
    disp(['evaluating fold: ' num2str(i)]); 
    trainind = find(ismember(protlabels,trainidx{i}));
    testind = find(ismember(protlabels,testidx{i}));
    
    trainlabels = classlabels(trainind);
    testlabels = classlabels(testind);

    trainprots = antibodyids(trainind);

    trainabids2 = trainprots;

    [tmp_traintrainidx, tmp_traintestidx, tmp_trainprotlabels, tmp_trainN] = splitbyproteins( trainlabels, trainabids2, Ninit-1); %%
    traintrainidx{i} = tmp_traintrainidx;
    traintestidx{i} = tmp_traintestidx;
    trainprotlabels{i} = tmp_trainprotlabels;
    trainN(i)  = tmp_trainN;
end



tmpidx1 = strfind(dataname,'.');
tmpidx2 = strfind(dataname,'/');
tmpname=dataname;
if isempty(tmpidx2)
   tmpidx2 = 0;
   tmpname = tmpname(1:tmpidx1(end)-1);
else
   tmpname(tmpidx2(end))='_';
   tmpname = tmpname(tmpidx2(end-1)+1:tmpidx1(end)-1);
end
save(['./data/classification/',tmpname,'_',num2str(Ninit),'folds.mat'],'trainidx','testidx','protlabels','N','classes','allfeatures',...
     'classlabels','antibodyids','dataname','Ninit','Nclass','Nfolds','class_del','traintrainidx','traintestidx','trainprotlabels','trainN','boolvar','numc');



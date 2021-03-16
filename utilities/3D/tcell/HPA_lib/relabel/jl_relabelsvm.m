function jl_relabelsvm(dataname,Ninit2,Nclass2,topN,numLevels)

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

tmpidx1 = strfind(dataname,'.');
tmpidx2 = strfind(dataname,'/');
tmpname=dataname;
if isempty(tmpidx2)
   tmpidx2 = 0;
else
   tmpname(tmpidx2(end))='_';
   tmpname = tmpname(tmpidx2(end-1)+1:tmpidx1(end)-1);
end

load(['./data/classification/svm',tmpname])


if ~exist('Ninit2','var')
   Ninit2 = 5;
end
if ~exist('Nclass2','var')
   Nclass2 = 5;
end
if ~exist('Ninit','var')||~exist('Nclass','var')
   Ninit = Ninit2;
   Nclass = Nclass2;
end
if Nclass<Ninit
   warning('Nclass cannot be smaller than Ninit and has been increased to equal to Ninit!');
   Nclass = Ninit;
end
Ninit
Nclass

if ~exist('allfeatures','var')||~exist('classes','var')||~exist('classlabels','var')||~exist('antibodyids','var')
   if exist(['./data/classification/',tmpname,'_',num2str(Ninit),'folds.mat'],'file')
      ['./data/classification/',tmpname,'_',num2str(Ninit),'folds.mat']
      load(['./data/classification/',tmpname,'_',num2str(Ninit),'folds.mat'],'allfeatures','classes','classlabels','antibodyids')
   else
      load(dataname,'allfeatures','classes','classlabels','antibodyids');
      [rows,cols] = find(isnan(allfeatures));
      allfeatures(rows,:) = [];
      classlabels(rows,:) = [];
      antibodyids(rows,:) = [];

      %allfeatures = single(allfeatures); %%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %delete classes with small number of antibodies.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
   end
end

   
idx_sda = cell(length(committee),1);
c = [];
g = [];
for i=1:length(committee)
    idx_sda{i} = committee(i).idx_feat;

    tmp_ind1 = strfind(committee(i).options,'-g');
    tmp_ind2 = strfind(committee(i).options,' ');
    tmp_ind3 = find(tmp_ind2-(tmp_ind1+3)>=0);
    tmp_ind4 = tmp_ind2(tmp_ind3(1));
    g = [g str2num(committee(i).options(tmp_ind1+3:tmp_ind4-1))];

    tmp_ind5 = strfind(committee(i).options,'-c');
    tmp_ind6 = find(tmp_ind2-(tmp_ind5+3)>=0);
    tmp_ind7 = tmp_ind2(tmp_ind6(1));
    c = [c str2num(committee(i).options(tmp_ind5+3:tmp_ind7-1))];
end
%idx_sda = unique(idx_sda);
%length_idx_sda = length(idx_sda)
%g = mode(g);
g
%c = mode(c);
c


%allfeatures = allfeatures(:,idx_sda); %%allfeatures has already been selected by SDA.
allfeatures = single(allfeatures); %%

%keyboard

all_cc = [];
all_uu = [];
all_dd = [];
all_ww = [];
[truelabels, predlabels, protlabels, allweights] = voteAcrossProtein( allweights, classlabels, antibodyids);
[cc uu dd ww] = conmatrix( truelabels, predlabels)
all_cc = [all_cc;{cc}];
all_uu = [all_uu;{uu}];
all_dd = [all_dd;{dd}];
all_ww = [all_ww;{ww}];

val = max(allweights,[],2);

index = find(predlabels~=truelabels);

[s index2] = sort(val(index));

rankedlist = protlabels(index(index2));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pseudorecurse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('topN','var')
   topN = 5;
end
topN
if ~exist('numLevels','var')
   numLevels = 5;
end
numLevels

%[newrankedlist] = pseudorecurse(allfeatures,idx_sda,g,Ninit,rankedlist,classlabels,antibodyids,protlabels,predlabels,topN,numLevels); %%

numclassif = pascal(numLevels + 1,2);
numclassif = abs(numclassif(:,1));

indexer = 1;
tmp_rankedlist = [];
tmp_rankedlist{1} = rankedlist;

%newclasslabels = classlabels; %%****************1 uses only this in save!!!!!!!!!!!!!!!!!!
for i=2:length(numclassif)
    newclasslabels = classlabels; %%*************2 uses only this in save!!!!!!!!!!!!!!!!!!Default, original and (maybe) useful.
    for j=indexer:topN-1+indexer
        % tmp = rankedlist;
        tmp = tmp_rankedlist{i-1};
        idx = find(tmp(j)==antibodyids);
        idx2 = find(tmp(j)==protlabels);
        newclasslabels(idx) = predlabels(idx2);
    end
    % indexer = indexer+topN;


    %[newrankedlist] = rerank(allfeatures,newclasslabels,antibodyids,idx_sda,g,Ninit); %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rerank
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %newallweights = svmCrossValidate( allfeatures, newclasslabels, antibodyids, idx_sda, g, Ninit); %%
    %newallweights = svmCrossValidate( allfeatures, newclasslabels, antibodyids, [1:size(allfeatures,2)], g, Ninit); %%
    %newallweights = svmCrossValidate( allfeatures, newclasslabels, antibodyids, [], c, g, Ninit); %%allfeatures has already been selected by SDA.
    newallweights = svmCrossValidate( tmpname, allfeatures, newclasslabels, antibodyids, idx_sda, c, g, Ninit); %%
    % correctHoles
    [newtruelabels, newpredlabels, newprotlabels, newallweights] = voteAcrossProtein( newallweights, newclasslabels, antibodyids);
    [cc uu dd ww] = conmatrix( newtruelabels, newpredlabels)
    all_cc = [all_cc;{cc}];
    all_uu = [all_uu;{uu}];
    all_dd = [all_dd;{dd}];
    all_ww = [all_ww;{ww}];

    newval = max(newallweights,[],2);

    pwrong = 1-newval;

    index = find(newpredlabels~=newtruelabels);

    [s index2] = sort(-pwrong(index));

    newrankedlist = newprotlabels(index(index2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    tmp_rankedlist{i} = newrankedlist;
    disp('.');
end
clear allfeatures

all_rankedlist = [];
for i=1:length(numclassif)
    %all_rankedlist = [all_rankedlist; repmat( tmp_rankedlist{i}, [numclassif(i) 1])]; %%
    all_rankedlist = [all_rankedlist; tmp_rankedlist{i}];
end

[c b] = hist(all_rankedlist,1:1:max(all_rankedlist));
counts = c(c>0);
prots = b(c>0);

newrankedlist = prots(counts==max(counts))';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



newlabels = [];
for i=1:length(newrankedlist)
    newlabels{i,1} = num2str(newrankedlist(i));
    idx = find(newrankedlist(i)==protlabels);
    tmp = classes{truelabels(idx)};
    newlabels{i,2} = tmp;
    tmp = classes{predlabels(idx)};
    newlabels{i,3} = tmp;
end

oldlabels = [];
for i=1:length(rankedlist)
    oldlabels{i,1} = num2str(rankedlist(i));
    idx = find(rankedlist(i)==protlabels);
    tmp = classes{truelabels(idx)};
    oldlabels{i,2} = tmp;
    tmp = classes{predlabels(idx)};
    oldlabels{i,3} = tmp;
end



%{
save(['./data/classification/relabelsvm2_',num2str(numLevels),dataname],'Nclass','Nclass2','Ninit','Ninit2','all_cc','all_dd','all_uu','all_ww','all_rankedlist',...
     'antibodyids','classlabels','class_del','classes','committee','counts','dataname','g','b','c','newclasslabels','newlabels','newrankedlist','numLevels','numc',...
     'numclassif','oldlabels','predlabels','protlabels','rankedlist','tmp_rankedlist','topN','truelabels');
%}
save(['./data/classification/relabelsvm_',num2str(numLevels),tmpname],'Nclass','Nclass2','Ninit','Ninit2','all_cc','all_dd','all_uu','all_ww','all_rankedlist',...
     'antibodyids','classlabels','classes','committee','counts','dataname','g','b','c','newclasslabels','newlabels','newrankedlist','numLevels',...
     'numclassif','oldlabels','predlabels','protlabels','rankedlist','tmp_rankedlist','topN','truelabels');


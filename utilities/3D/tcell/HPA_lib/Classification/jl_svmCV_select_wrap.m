function jl_svmCV_select_wrap(dataname,CV,flag,pattern,rootdir)

%'Full' two-level cross validation; 
%in outer level, each test fold could use different pairs of parameters which is not the case in jl_svmCV_select (it will select 
%only one pair of parameters for all 5 'outer test' folds).

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
   tmpname = tmpname(1:tmpidx1(end)-1);
else
   tmpname(tmpidx2(end))='_';
   tmpname = tmpname(tmpidx2(end-1)+1:tmpidx1(end)-1);
end
tmpname

tmpidx1 = strfind(dataname,'.');
tmpidx3 = strfind(dataname,'Dist');
tmpname2 = dataname;
boottag = tmpname2(tmpidx3+4:tmpidx1(end)-1)
tmpname2 = [tmpname2(1:tmpidx3-2),'.mat']


if ~exist('CV','var')
   CV = 5;
end
CV

if ~exist('flag','var')
   flag = 0;
end
flag

if ~exist('rootdir','var')
   rootdir = './data/classification/';
end
rootdir

if ~exist('pattern','var')
   pattern = 'svm';
end
pattern



if flag == 0

currentdir = pwd;
cd(rootdir)

funcs = cell(CV,1);
jobs = cell(CV,1);

for k = 1:CV

%svmcg_files = dir([pattern,'*_',tmpname,'.mat']); %%
%svmcg_files = dir([pattern,'*',tmpname,'.mat']); %%
svmcg_files = dir([pattern,'*','-',num2str(k),tmpname,'.mat']); %%
length(svmcg_files);

%keyboard
accuracies_cg = [];
for i=1:length(svmcg_files)
%    load(svmcg_files(i).name,'uu3');
%    uu3
%    accuracies_cg = [accuracies_cg;uu3];

    %{
    load(svmcg_files(i).name,'acc_tuning_ave');
    acc_tuning_ave
    accuracies_cg = [accuracies_cg;acc_tuning_ave];
    %}
    load(svmcg_files(i).name,'acc_tuning');
    acc_tuning;
    accuracies_cg = [accuracies_cg;acc_tuning(k)];  %%
end
[svmcg_accuracy,svmcg_idx] = max(accuracies_cg);
svmcg_accuracy;
svmcg_idx;
svmcg_files(svmcg_idx);
load(svmcg_files(svmcg_idx).name,'c_select','g_select');
c_select = c_select(k)
if exist('g_select','var')
   g_select = g_select(k)
end
length_svmcg_files = length(svmcg_files)

disp(['jl_svmCV(',num2str(c_select),',',num2str(g_select),',','[]',',','0,','''',dataname,'''',',',num2str(k),')']);
disp(['qarray',' ','-m',' ','"','jl_svmCV(',num2str(c_select),',',num2str(g_select),',','[]',',','0,','''',dataname,'''',',',num2str(k),');','"']);

funcs{k} = ['jl_svmCV(',num2str(c_select),',',num2str(g_select),',','[]',',','0,','''',dataname,'''',',',num2str(k),')'];
jobs{k} = ['qarray',' ','-m',' ','"','jl_svmCV(',num2str(c_select),',',num2str(g_select),',','[]',',','0,','''',dataname,'''',',',num2str(k),');','"'];


%%delete unuseful files
tmpfiles = svmcg_files;
tmpfiles(svmcg_idx) = [];
%{
for i=1:length(tmpfiles)
    delete(tmpfiles(i).name);
end
%}




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the final accuracy for report
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%run job
%qarray -m "jl_svmCV(2^10,2^0,29,0,'./data/IF_atlas_5/regionfeatures_single_A431.mat')"
%Or, 
%qarray -m "jl_svmCV(2^10,2^0,[],0,'./data/IF_atlas_5/regionfeatures_single_A431.mat')"




%{
%keyboard

disp(['jl_svmDistCV(',num2str(mode(c_select)),',','''',boottag,'-','''',',','0,','''',tmpname2,'''',')']);
%keyboard
jl_svmDistCV(mode(c_select),[boottag,'-'],0,tmpname2);
%}

%keyboard

end

save([pattern tmpname '.mat'], 'funcs','jobs');



else

addpath(genpath('/home/jieyuel/HPA_Project/lib_Justin/lib/classify'));
addpath('/home/jieyuel/HPA_Project/lib_Justin/lib/results');

addpath('/home/jieyuel/HPA_Project/lib_Justin/HPAconfocal_lib/classification')

currentdir = pwd;
cd(rootdir)

load([tmpname,'_',num2str(CV),'folds.mat'],'classlabels','antibodyids','classes','Ninit','Nclass','Nfolds');
%keyboard
predlabelsall = zeros(length(classlabels),1);
allweightsall = zeros(length(classlabels),length(unique(classlabels)));
committeeall = struct([]);
%svmcg_files = dir([pattern,'-','*',tmpname,'.mat']); %%
for k = 1:CV
    CV_file = [pattern,'-',num2str(k),tmpname,'.mat']; %%
    load(CV_file,'predlabels','allweights','SDA','committee')
    predlabelsall = predlabelsall+predlabels;
    allweightsall = allweightsall+allweights;
    committeeall = [committeeall;committee(k)];
end
predlabels = predlabelsall;
allweights = allweightsall;
clear committee
committee = committeeall;



[cc1 uu1 dd1 ww1] = conmatrix( classlabels, predlabels)


[protlabels,protpredlabels,proteins,newallweights] = ...
    modeAcrossProtein( allweights, classlabels, antibodyids);

[cc2 uu2 dd2 ww2] = conmatrix( protlabels, protpredlabels);
cc2
dd2


[protlabels,protpredlabels,proteins,newallweights] = ...
    voteAcrossProtein( allweights, classlabels, antibodyids);

[cc3 uu3 dd3 ww3] = conmatrix( protlabels, protpredlabels);
cc3
dd3


[tmp_protlabels,tmp_protpredlabels,proteins] = ...
    jl_voteAcrossProtein( classlabels, antibodyids, predlabels);

[tmp_cc3 tmp_uu3 tmp_dd3 tmp_ww3] = conmatrix( tmp_protlabels, tmp_protpredlabels)




save([pattern tmpname '.mat'], 'allweights','antibodyids','classlabels','committee',...
      'cc1','uu1','dd1','ww1','cc2','uu2','dd2','ww2','cc3','uu3','dd3','ww3','classes','Ninit','Nclass','Nfolds','predlabels',...
      'protlabels','protpredlabels','SDA','tmp_cc3','tmp_uu3','tmp_dd3','tmp_ww3','tmp_protlabels','tmp_protpredlabels');

end

cd(currentdir)





function [protlabels,protpredlabels,proteins] = jl_voteAcrossProtein(classlabels,antibodyids,predlabels)

proteins = unique(antibodyids);
protpredlabels = zeros(size(proteins));
protlabels = zeros(size(proteins));
for i=1:length(proteins)
    idx = find(proteins(i)==antibodyids);

    protpredlabels(i) = mode(predlabels(idx));
    protlabels(i) = mode(classlabels(idx));
end


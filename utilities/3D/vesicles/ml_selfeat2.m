function [idx,scores] = ml_selfeat2(feats1,feats2,param)
%ML_SELFEAT2 Select similar features between two data sets
%   IDX = ML_SELFEAT2(FEATS1,FEATS2) returns the indices of the features that
%   are supposed to similar between the [MCF] FEATS1 and [MCF] FEATS2. The
%   indices are sorted with the descending order of similarities. FEATS1 and
%   FEATS2 should contain the same number of [feature matrix]s.
%   
%   IDX = ML_SELFEAT2(FEATS1,FEATS2,PARAM) specifies how to select the features
%   by the structure PARAM, which has the following fields:
%       'method' - methods of selecting and sorting the features
%           'avgd' : average z-score distance
%           'avgr' : average rank
%           'avgt' : average t-test scores
%           'sda1' : sda from FEATS1
%           'sda2' : sda from FEATS2
%           'sdap' : sda from pooled data (default)
%           'sdas' : sda from pooled data of common classes
%           'tsda' : t-test first then sda
%           'mahd' : mahalanobis distance
%           'ldac' : lda classification
%      'classes1' - class names for FEATS1. This is a [string array].
%      'classes2' - class names for FEATS2. This is a [string array].
%      If either 'classes1' or 'classes2' does not exist, FEATS1 and FEATS2
%      must have the same number of [feature matrix]s so that the elements at
%      the same position correspond to the same class.
%
%   [IDX,SCORES] = ML_SELFEAT2(...) returns the sorted scores as will.
%
%   See also

%   10-Jan-2007 Initial write T. Zhao
%   Copyright (c) 2007 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.
%   
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu


if nargin < 2
    error('2 or 3 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct('method','sdap'));

if ~isfield(param,'classes1') |  ~isfield(param,'classes2')
    if length(feats1)~=length(feats2)
        error(['The first two inputs must have the same number of ' ...
               'elements while the classes are not specified']);
    else
        classidx = 1:length(feats1);      
    end
else
    for i=1:length(param.classes2)
        classidx(i) = strmatch(param.classes2{i},param.classes1);
    end
end

if ismember(param.method,{'avgd','avgr'})
    sortresult.indices = [];
    sortresult.dists = [];
    sortresult.rank = [];
    for i=1:length(classidx)
        [idx,dists,r] = ml_sortfeatdist(feats2{i},feats1{classidx(i)});
        sortresult.indices(i,:) = idx;
        sortresult.dists(i,idx) = dists;
        sortresult.rank(i,:) = r;      
    end
end

switch param.method
    case 'rand'
        idx = randperm(size(feats1{1},2));
        scores = 1:length(idx);
    case 'ldac'
        nfeat = size(feats1{1},2);
        avgaccs = [];
        [combfeats1,combclass1] = ml_mcf2combfeats(feats1);
        [combfeats2,combclass2] = ml_mcf2combfeats(feats2);
        errs = [];
        for i=1:nfeat
            sample = combfeats2(:,i);
            training = combfeats1(:,i);
            [testidx,err] = classify(sample,training,combclass1);
            avgaccs(i) = sum(combclass2==testidx)/length(combclass2);
        end
%         [scores,idx] = sort(errs);
        [scores,idx] = sort(avgaccs,2,'descend');
    case 'avgd'
        meandists = mean(sortresult.dists,1);
        [scores,idx] = sort(meandists);
    case 'avgr'
        meanrank = mean(sortresult.rank,1);
        [scores,idx] = sort(meanrank);
    case 'avgt'
        param = ml_initparam(param,struct('thresh',[]));
        nfeat = size(feats1{1},2);
        dists = [];
        ps = [];
        for i=1:length(classidx)
            for j=1:nfeat
                f1 = feats1{classidx(i)}(:,j);
                f2 = feats2{i}(:,j);
                %dists(i,j) = abs(mean(f1)-mean(f2))/sqrt(var(f1)+var(f2));
                [ps(i,j),ts] = ml_ttest2(f1,f2);
                ts(i,j) = abs(ts);
            end
        end
        meanps = mean(ps,1);
        [scores,idx] = sort(-meanps);
        scores = -scores;
        if ~isempty(param.thresh)
            idx(-sortedps<param.thresh) = [];
            scores(-sortedps<param.thresh) = [];
        end
        
    case 'sda1'
        idx = ml_stepdisc(feats1);
        scores = [];
    case 'sda2'
        idx = ml_stepdisc(feats2);
        scores = [];
    case 'sdap'
        combfeats = feats1;
        for i=1:length(classidx)
            combfeats{classidx(i)} = [feats1{classidx(i)};feats2{i}];
        end
        idx = ml_stepdisc(combfeats);
        scores = [];
    case 'sdas'
        combfeats = {};
        for i=1:length(classidx)
            combfeats{i} = [feats1{classidx(i)};feats2{i}];
        end
        idx = ml_stepdisc(combfeats);
        scores = [];
    case 'tsda'
        param2 = param;
        param.method = 'avgt';
        [idx2,scores] = ml_selfeat2(feats1,feats2,param);
        nfeat = size(feats1{1},2);
        idx2 = idx2(1:round(nfeat/2));
        feats1 = ml_selfeat_mcf(feats1,idx2);
        idx3 =ml_stepdisc(feats1);
        idx = idx2(idx3);
        scores = scores(idx3);
    case 'mahd'
        nfeat = size(feats1{1},2);
        dists = [];
        nsamples = [];
        for i=1:length(classidx)
            nsamples(i) = size(feats1{classidx(i)},1);
        end
        
        for i=1:nfeat
            ufeats1 = [];
            ufeats2 = [];
            for j=1:length(classidx)
                ufeats1(:,j) = feats1{classidx(j)}(1:min(nsamples),i);
                ufeats2(:,j) = feats2{j}(:,i);
            end
            dists(i) = ml_twomaha(ufeats1,ufeats2,0,1);
        end
        [scores,idx] = sort(dists);
    otherwise
        error(['Unrecognized method: ' param.method]);
end

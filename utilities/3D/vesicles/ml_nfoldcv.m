function [avgcm,avgacc,ncm,fveclen,otherinfo] = ...
    ml_nfoldcv(features,n,regmeth,t)
%ML_NFOLDCV N-fold cross validation of classification.
%   AVGCM = ML_NFOLDCV(FEATURES,N,REGMETH,T) returns the averaged
%   confusion matrix of N fold cross validation of classifying FEATURES,
%   which is a cell array of feature matrices. Each matrix in the cell
%   array FEATURES contains features from a unique class. REGMETH is
%   the classification method and T is the structure of classificaion 
%   parameters. If T is empty, default parameters will be used. 
%   See ML_REGRESS for more details.
%
%   [AVGCM,AVGACC,NCM,FVECLEN,OTHERINFO] = ML_NFOLDCV(...) also returns 
%   overall accuracy AVGACC and confusion matrix with numbers, NCM. FVECLEN
%   is the number of features wich give the best accuracy. OTHERINFO
%   contains information about which samples are not correctly classified.

%   22-Jun-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

% Copyright (C) 2007  Murphy Lab
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

if nargin < 4
    error('Exactly 4 arguments are required')
end

nclass=length(features);

rand('state',0);
for i=1:nclass
    norgset(i)=size(features{i},1);
    perm{i}=randperm(norgset(i));
    permfold(i,:)=[0,cumsum(ml_redistrnum(norgset(i),n))];
end

test_sda_len = size(features{1},2);
stop_len = 15;

if isempty(t)
    t = struct([]);
end
if ~isfield(t,'sda')
    t = ml_initparam(t,struct('sda',[]));
end

if isempty(t.sda)
    avgaccs = zeros(1,1);
else
    %stop_len = 15;  Y.H.: Mar 07
    avgaccs = zeros(1,test_sda_len);
end

for len=1:length(avgaccs)
    for k=1:n
        k
        for i=1:nclass
            pf=perm{i};
            testsel{k,i}=pf(permfold(i,k)+1:permfold(i,k+1));            
            testset{i}=features{i}(testsel{k,i},1:test_sda_len-len+1);
            pf(permfold(i,k)+1:permfold(i,k+1))=[];
            trainsel{k,i} = pf;
            trainset{i}=features{i}(pf,1:test_sda_len-len+1);
        end
        [trainfeats,trainclass,traincellidx]=ml_mcf2combfeats(trainset);
        [testfeats,testclass,testcellidx]=ml_mcf2combfeats(testset);
        if regmeth==1
            for u=1:max(trainclass)
                if sum(trainclass==u)<=size(trainfeats,2)
                    error(['Cross validation error - ' ...
                        'The number of images must exceed ' ...
                        'the number features for LDA.']);
                end
            end
        end
        regmodel=ml_trainclassif(trainfeats,trainclass,t,regmeth);
        tg=ml_evalreg(testfeats,regmodel);
        [ncvcm{k},pcvcm{k},cvacc(k),errinfo{k}] = ...
            ml_summaryclassif(testclass,tg,nclass);
        if ~isempty(errinfo{k})
            errinfo{k}(:,1) = testcellidx(errinfo{k}(:,1));
        end
    end    
    [avgcms{len},avgaccs(len),ncms{len}]=ml_calcavgcm(ncvcm);

    if len >= stop_len
        break;
    end
end
[output, bestidx] = max(avgaccs);
fveclen = test_sda_len-bestidx+1;
avgcm = avgcms{bestidx};
avgacc = avgaccs(bestidx);
ncm = ncms{bestidx};
otherinfo = struct('errinfo', {errinfo},'trainsel',{trainsel}, ...
    'testsel',{testsel});


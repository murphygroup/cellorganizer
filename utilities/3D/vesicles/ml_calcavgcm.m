function [avgcm,avgacc,ncm,kappa,ua,pa]=ml_calcavgcm(cvcm,merge)
%ML_CALCAVGCM Summarize confusion matrices from cross validation.
%   AVGCM = ML_CALCAVGCM(CVCM) returns the average confusion matrix from
%   the cell array of confusion matrix, CVCM.
%   
%   AVGCM = ML_CALCAVGCM(CVCM,MERGE) will merge classes with the indices
%   in the vector MERGE. For example, if MERGE is [1 3 4], classes 1, 3
%   and 4 will be merged.
%
%   [AVGCM,AVGACC,NCM,KAPPA,UA,PA] = ML_CALCAVGCM(...) also returns
%   overall accuracy AVGACC, confusion matrix with numbers NCM, kappa
%   values KAPPA, user's accuracy UA, producers' accuracy PA.
%   The Kappa is calculated as 
%       (Observed agreement - Chance agreement)/(1-Chance agreement).
%   It indicates how good the agreement is (from Landis and Koch, 1977):
%       0.00-poor-0.01-slight-0.20-fair-0.40-moderate-0.60-substantial-
%       0.80-almost perfect-1.00

%   ??-???-???? Initial write T. Zhao
%   17-MAY-2004 Modified T. Zhao
%       - change parameter meaning
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

if nargin < 1
    error('1 or 2 arguments are required')
end

if ~exist('merge','var')
    merge=[];
end

%initialization
for i=1:length(cvcm)
    if ~isempty(cvcm{i})
        ncm=zeros(size(cvcm{i}));
        break;
    end
end

%add confusion matrices of all folds
for i=1:length(cvcm)
    if ~isempty(cvcm{i})
        ncm=ncm+cvcm{i};
    end
end

%merge classes
if ~isempty(merge)
    if max(merge)>length(ncm)
        error('Merging indices exceed the number of confusion matrix');
    end
        
    ncm2=ncm;
    ncm2(merge(1),:)=sum(ncm(merge,:),1);
    ncm2(merge(2:end),:)=[];
    
    ncm2(:,merge(1))=sum(ncm2(:,merge),2);
    ncm2(:,merge(2:end))=[];
    
    ncm=ncm2;
end

%calculate accurate rate
avgcm=ncm./repmat(sum(ncm,2),1,size(ncm,2))*100;
avgacc=sum(diag(ncm))/sum(ncm(:));

srow=sum(ncm,2);
scol=sum(ncm,1)';

%users accuracy and producers accuracy
if scol==0
    ua = 0;
else
    ua=(diag(ncm)./scol)';
end
if srow==0
    pa = 0;
else
    pa=(diag(ncm)./srow)';
end

%kappa
if length(ncm)==1
    kappa = 1;
else
    n=sum(ncm(:));
    po = n*sum(diag(ncm));
    pc = sum(srow.*scol);
    kappa = (po-pc)/(n^2-pc);
end

% kappa = (n*sum(diag(ncm))-sum(srow.*scol));
% if kappa~=0
%     kappa=(n*sum(diag(ncm))-sum(srow.*scol))/(n^2-sum(srow.*scol));
% end


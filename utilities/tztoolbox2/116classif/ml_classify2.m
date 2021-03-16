function result = ml_classify2(feats1,feats2,param)
%ML_CLASSIFY2 Train classifiers on one data set and classify the other.
%   RESULT = ML_CLASSIFY2(FEATS1,FEATS2) trains a classifier on FEATS1 and test
%   it on FEATS2. FEATS1 and FEATS2 are both [MCF]s. It returns a structure
%   containing the classification results.
%   
%   RESULT = ML_CLASSIFY2(FEATS1,FEATS2,PARAM) lets users customize the 
%   parameters by the structure PARAM, which could have the following fields:
%       'classifierid' - ID of the classifier. See ML_REGRESS for more details.
%       't' - additinoal parameters of the classifier. See ML_REGRESS for more 
%           details.
%       'cv1' - do cross validation for FEATS1 if it is 1.
%       'cv2' - do cross validation for FEATS2 if it is 1.
%
%   See also

%   11-Jan-2007 Initial write T. Zhao
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

param = ml_initparam(param,struct('classifierid',2,'cv1',1,'cv2',1,'t',[]));

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

[combfeats1,labels] = ml_mcf2combfeats(feats1);

classifier = ml_trainclassif(combfeats1, ...
         labels,param.t,param.classifierid);

if param.cv1==1
    [result.avgcm1,result.avgacc1,result.ncm1] = ml_nfoldcv(feats1, ...
           10,param.classifierid,param.t);
end

if param.cv2==1
    [result.avgcm2,result.avgacc2,result.ncm2] = ml_nfoldcv(feats2, ...
              10,param.classifierid,param.t);
end

result.tcm = [];
result.tacc = [];
for i=1:length(classidx)
    ty = classidx(i);
    ey = ml_evalreg(feats2{i},classifier);

    %average accuracy
    acc = sum(ey==ty)/length(ey)
    
    %count numbers in the training classes
    s = ml_countnum(ey);
    
    if min(ey)>1
        s = [zeros(1,min(ey)-1) s];
    end
    if max(ey)<10
        s = [s zeros(1,length(feats1)-max(ey))];
    end
    result.tcm = [result.tcm;s];
    result.tacc(i) = acc;
end

result.classidx = classidx;
result.classifier = classifier;

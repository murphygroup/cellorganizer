function [ncm,pcm,avgacc,errinfo]=ml_summaryclassif(tlabel,elabel,nclass)
%ML_SUMMARYCLASSIF Summarize classification results.
%   NCM = ML_SUMMARYCLASSIF(TLABEL,ELABEL) returns a confusion matrix
%   with numbers by checking true class labels TLABEL and predicted class
%   labels ELABEL. TLABEL and ELABEL must be column vectors of positive
%   integers with the same length.
%
%   NCM = ML_SUMMARYCLASSIF(TLABEL,ELABEL,NCLASS) specifies that there are
%   NCLASS classes. So the size of NCM will be NCLASS X NCLASS.
%
%   [NCM,PCM,AVGACC] = ML_SUMMARYCLASSIF(TLABEL,ELABEL) also returns
%   confusion matrix with percentage PCM and overall accuracy AVGACC.
%
%   See also ML_SUMMARYTESTCLASSIF

%   17-MAY-2004 Modified T. Zhao
%       - return confusion matrix of numbers 
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

if nargin < 2
    error('Exactly 2 arguments are required')
end

if size(tlabel,2) > 1
    error('The 1st argument is not a column vector.');
end

if size(elabel,2) > 1
    error('The 2nd argument is not a column vector.');
end

if length(tlabel) ~= length(elabel)
    error('The two arguments do not have the same length.');
end

%find distinguished classes
if ~exist('nclass','var')
    nclass = max([tlabel; elabel]);
end

% caclass=ml_findclass(tlabel); 
% nclass=length(caclass);

nsample=length(elabel);

ncm=zeros(nclass);
errinfo = [];

for i=1:nsample
    ncm(tlabel(i),elabel(i))=ncm(tlabel(i),elabel(i))+1;
    if tlabel(i)~=elabel(i)
        errinfo = [errinfo;i tlabel(i) elabel(i)];
    end
end

pcm=ml_normrow(ncm)*100;
avgacc=sum(tlabel==elabel)/nsample;
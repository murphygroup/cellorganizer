function out = tz_ncclassif(x,y,param)
%TZ_NCCLASSIF Train a nearest center classifier.
%   OUT = TZ_NCCLASSIF(X,Y) returns a nearest center classifier trained
%   from the [feature matrix]  X and [label vector] Y. OUT is a structure
%   containing information of the classiifer. It has the following fields:
%       'centers' - it is a matrix. The ith row is the center of a group
%           with label i.
%   
%   OUT = TZ_NCCLASSIF(X,Y,PARAM) also allows setting paramters. PARAM is a
%   structure and it has the following fields:
%       NOT AVAILABLE YET
%   
%   See also

%   03-Mar-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('2 or 3 arguments are required')
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct([]));

%Group data into a cell array
x2 = ml_combfeats2mcf(x,y);

%Calculate centers
for i=1:length(x2)
    out.centers(i,:) = mean(x2{i},1);
end
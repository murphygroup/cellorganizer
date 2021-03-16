function y = tz_data2pca(x,param)
%TZ_DATA2PCA Convert data into principle component representation.
%   Y = TZ_DATA2PCA(X) returns a [feature matirx]  which is the principle
%   component representation of the [feature matrix] X. 
%   
%   Y = TZ_DATA2PCA(X,PARAM) let user specify parameters such as number of
%   components or percentage:
%       'reduce' - how to reduce the data
%           'ncmp' : number of components
%           'perc' : percentage
%           'none' : no reductant
%       'reduceparam' - number of components if PARAM.reduce is 'ncmp' and
%           percentage if PARAM.reduce is 'perc'
%   
%   See also

%   29-May-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('1 or 2 arguments are required')
end

if size(x,1)==1
    warning('There is only one data point. PCA does not make any change.');
    y = x;
    return
end

[pcavec,coeff] = princomp(x);
switch param.reduce
    case 'ncmp'
        y = coeff(:,1:param.reduceparam);
    case 'perc'
        error('Sorry, option ''perc'' is not available currently!');
    case 'none'
        y = coeff;
    otherwise
        error('Unrecognized reduce method');
end


function tz_plotgroup(x,group,dims)
%TZ_PLOTGROUP Plot grouped data.
%   TZ_PLOTGROUP(X,GROUP) plot data matrix X with different colors and
%   markes according to their group membership in the integer vector
%   GROUP. If the number of columns is greater than 2, the first two
%   principal components will be plotted.
%   
%   TZ_PLOTGROUP(X,GROUP,DIMS) plots the data with specified dimensions
%   in DIMS. PCA is also used if the number of dimensions is greater
%   than 2.

%   03-Jun-2005  Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University
  
if nargin < 2
    error('2 or 3 arguments are required')
end

if exist('dims','var')
    x=x(:,dims);
end

if size(x,2)==1
    x(:,2) = zeros(size(x,1),1);
end

if size(x,2)>2
    [coeff,score]=princomp(x);
    x=score(:,1:2);
end

gcolors='ykrbgmc';
gmarks='o.x+*sdhp';
ncolors=length(gcolors);
nmarks=length(gmarks);

if iscell(group)
    ngroup=length(group);
else
    ngroup=max(group);
end

for i=1:ngroup
%     [marksel,colorsel]=ind2sub([nmarks,ncolors],i);
    marksel=mod(i,nmarks)+1;
    colorsel=mod(i,ncolors)+1;
    if iscell(group)
        groupsel=group{i};
    else
        groupsel=find(group==i);
    end
    plot(x(groupsel,1),x(groupsel,2),[gcolors(colorsel) gmarks(marksel)]);
    hold on
end
hold off

function osm = tz_trainosm(combclass,combcellidx,combdist,clstlabel)
%TZ_TRAINOSM Train object sampling model (OSM)
%   OSM = TZ_TRAINOSM(COMBCLASS,COMBCELLIDX,COMBDIST,CLSTLABEL) returns
%   an object sampling model for multiple classes. OSM is an structure
%   array containing the following fields:
%       'occurclst' - a vector of occured clusters in each class
%       'nmean' - a vector of mean of object number
%       'nvar' -  the covariance of object number
%       'distpara' - parameters for distances
%
%   See also TZ_GENOSMPARAM

%   22-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 4
    error('Exactly 4 arguments are required')
end

cellobjcom = tz_cellobjcom(combclass,combcellidx, ...
    ml_label2post(clstlabel),max(clstlabel),'#$objnum', ...
    [],[],[],{},{});

for i=1:length(cellobjcom)
    objnums = cellobjcom{i};
    occuridx{i} = find(sum(objnums,1)>0);
    occurnums = objnums(:,occuridx{i});
    nmean{i}=mean(occurnums,1);
    nvar{i}=cov(occurnums);
end

distpara = reshape(tz_clstfeat(combdist,clstlabel, ...
    {'mean','std'}),2,[])';

osm = struct('occurclst',occuridx,'nmean',nmean,'nvar',nvar, ...
    'distpara',distpara);

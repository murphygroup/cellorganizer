function [feats,names,slfnames] = tz_regionfeatmat(img,dnaimg,masks)
%TZ_REGIONFEATMAT Calculate feature matrix of regions in an image.
%   FEATS = TZ_REGIONFEATMAT(IMG,DNAIMG,MASKS) returns a feature matrix of 
%   the regions in IMG. DNAIMG is an image from dna channel and used in 
%   feature calculation. See TZ_MASKIMGS for more information about the two
%   parameters IMG and MASKS.
%   
%   [FEATS,NAMES,SLFNAMES] = TZ_REGIONFEATMAT(...) also returns feature
%   names and SLF names.
%
%   See also

%   11-Nov-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end

regions = tz_maskimgs(img,masks);

if ~isempty(dnaimg)
    dnas = tz_maskimgs(dnaimg,masks);
end

feats = [];
names = [];
slfnames = [];

for k=1:length(regions)
    if ~isempty(dnas)
        dnaRegion = dnas{k};
    else
        dnaRegion = [];
    end
    
    [names, feats(k,:), slfnames] = ml_featset( double(regions{k}), [], ...
        double(dnaRegion), 'all180',0,0,'yesbgsub','rc');  
end

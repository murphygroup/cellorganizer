function featname = tz_yeastimg2featname(imgpath)
%TZ_YEASTIMG2FEATNAME Image file path to feature file name
%   FEATNAME = TZ_YEASTIMG2FEATNAME(IMGPATH) converts the path of a [yeast
%   image] into a unique file name with the information of chromosome and
%   image file name.
%   
%   See also

%   12-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

imagePathTokens = tz_strtok(imgpath,'/');
chromosomeName = imagePathTokens{end-1};
filename = strtok(imagePathTokens{end},'.'); %without extension

featname = [chromosomeName '_' filename];


function ey = tz_testncclassif(x,classifier)
%TZ_TESTNCCLASSIF Test a trained nearest center classifier.
%   EY = TZ_TESTNCCLASSIF(X,CLASSIFIER) returns the testing [label vector] 
%   of the [feature matrix] x from the nearest center classifier (NCC)
%   CLASSIFIER.
%   
%   See also

%   03-Mar-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end

ds = tz_pdist2(x,classifier.centers);

[mindist,ey] = min(ds,[],2);
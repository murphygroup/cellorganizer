function classobjcom=tz_classobjcom(post,combclass)
%TZ_CLASSOBJCOM Sort the combined clusters into cell arrays.
%   CLASSOBJCOM = TZ_CLASSOBJCOM(POST,COMBCLASS) returns a cell array of
%   matrices. Actually it seems that it is the same as ML_COMBFEATS2MCF.

%   03-MAR-2004 Initial write T. Zhao
%   03-NOV-2004 Modified T. Zhao
%       - add comments
%       - change parameters (classes,post) --> (post,combclass)
%   Copyright (c) Murphy Lab, Carnegie Mellon University

% warning('the function has been updated')

nclass=max(combclass);

for i=1:nclass
    classobjcom{i}=post(combclass==i,:);
end

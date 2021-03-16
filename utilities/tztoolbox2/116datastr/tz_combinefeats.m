function feats = tz_combinefeats(feat1,feat2)
%TZ_COMBINEFEATS Combine two feature matrices.
%   FEATS = TZ_COMBINEFEATS(FEAT1,FEAT2) combines two feature matrices
%   FEAT1 and FEAT2 along column and return the combined feature matrix
%   FEATS. This is supposed to combine different feature set on the same
%   data set.
%   
%   See also TZ_COMBFEATS2_MCF

%   21-MAY-2003 Initial write  TINGZ
%   14-DEC-2003 Modified T. Zhao
%   19-MAY-2004 Modified T. Zhao
%       - warning to error
%   31-OCT-2004 Modified T. Zhao
%       - add description
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

if size(feat1,1) ~= size(feat2,1)
    error('Sample size unmatched. Combining failed')
end

feats=[feat1 feat2];

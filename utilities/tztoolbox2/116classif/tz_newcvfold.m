function [testsel,trainsel] = tz_newcvfold(nsample,nfold,randstate)

%TZ_NEWCVFOLD Creates permutations for crss validation.
%   TESTSEL = TZ_NEWCVFOLD(NSAMPLE,NFOLD) returns the indices of testing
%   set in NFOLD fold cross validation. The indices are generated
%   randomly. TESTSEL is a cell array, which contains NFOLD vectors,
%   corresponding to each fold.
%   
%   TESTSEL = TZ_NEWCVFOLD(NSAMPLE,NFOLD,RANDSTATE) lets user specify
%   random state by RANDSTATE, which is used for rand('state',RANDSATE)
%   
%
%   [TESTSEL,TRAINSEL] = TZ_NEWCVFOLD(...) also returns indices of 
%   training set in TRAINSEL, which is also a cell array.
%   
%   See also TZ_CVPERMUTE

%   15-May-2005 Initial write Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University


if nargin < 2
    error('2 or 3 arguments are required')
end

if ~exist('randstate','var')
    sampleperm=randperm(nsample);
else
    if randstate>=0
        rand('state',randstate);
        sampleperm=randperm(nsample);
    else
        sampleperm=1:nsample;
    end
end

folds=ml_redistrnum(nsample,nfold);
stepfolds=cumsum([0,folds]);
fullsel=1:nsample;
for i=1:length(stepfolds)-1
    testsel{i}=sort(sampleperm(stepfolds(i)+1:stepfolds(i+1)));
    trainsel{i}=fullsel;
    trainsel{i}(testsel{i})=[];
end
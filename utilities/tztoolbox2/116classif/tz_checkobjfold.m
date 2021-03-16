function checkcode = tz_checkobjfold(dir1,dir2,nfold)
%TZ_CHECKOBJFOLD Compare two permutations for objects.
%   CHECKCODE = TZ_CHECKOBJFOLD(DIR1,DIR2,NFOLD) compares object cv
%   permutations under DIR1 and DIR2 and returns an integer value
%   CHECKCODE:
%       0 - exactly the same
%       1 - different permutation, same data size
%       2 - different data size
%   NFOLD is the number of folds.
%
%   See also TZ_CHECKCVFOLD

%   02-Apr-2005  Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('At least 3 arguments are required')
end

checkcode=0;
for i=1:nfold
    perm1=load([dir1 '/' 'sel' num2str(i) 'fold.mat']);
    perm2=load([dir2 '/' 'sel' num2str(i) 'fold.mat']);
    if (length(perm1.trainsel) ~= length(perm2.trainsel)) |...
            (length(perm1.testsel) ~= length(perm2.testsel))
        checkcode=2;
        break;
    end
    
    if (any(perm1.trainsel ~= perm2.trainsel) |...
            (any(perm1.testsel ~= perm2.testsel)))
        checkcode=1;
        break;
    end
end

switch checkcode
case 0
    disp('The two permutations are exactly the same');
case 1
    disp('The two permutations are different, but have the same data size');
case 2
    disp('The data sizes are different');
end
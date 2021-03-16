function checkcode = tz_checkcvfold(cv1,cv2)
%TZ_CHECKCVFOLD Compare two permutations for cross validation.
%   CHECKCODE = TZ_CHECKCVFOLD(CV1,CV2) returns an integer indicating
%   if the two permutations CV1 and CV2 are the same. If CHECKCODE is:
%       0 - exactly the same
%       1 - same size, different numbers of fold
%       2 - same size, different permutation
%       3 - different size
%   CV1 and CV2 are cell arrays of vectors of indices.    

%   25-May-2005  Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('At least 2 arguments are required')
end

n1=0;
n2=0;
for i=1:length(cv1)
    n1=n1+length(cv1{i});
end

for i=1:length(cv2)
    n2=n2+length(cv2{i});
end

checkcode=0;

if n1~=n2
    checkcode=3;
else
    if length(cv1)~=length(cv2)
        checkcode=1;
    else  
        for i=1:length(cv1)
            if length(cv1{i})~=length(cv2{i})
                checkcode=2;
                break;
            end
            if any(sort(cv1{i})~=sort(cv2{i}))
                checkcode=2;
                break;
            end
        end
    end
end

%       0 - exactly the same
%       1 - same size, different numbers of fold
%       2 - same size, different permutation
%       3 - different size

switch checkcode
case 0
    disp('The two permutations are exactly the same');
case 1
    disp('The two permutations have different numbers of folds, but have the same data size');
case 2
    disp('The two permutations are different, but have the same data size');
case 3
    disp('The data sizes are different');
end

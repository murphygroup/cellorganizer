function [celltrainsel,celltestsel] = tz_obj2cellperm(combcellidx, ...
    combclass,nfold,option,loadpath)
%TZ_OBJ2CELLPERM Convert object permutation to cell permutation.
%   CELLTRAINSEL = 
%       TZ_OBJ2CELLPERM(COMBCELLIDX,COMBCLASS,NFOLD,OPTION,LOADPATH)
%   returns indices of training set converted from permutation files
%   under LOADPATH. CELLTRAINSEL{J}{I} is a vector of indices in class J
%   and fold I. See TZ_OBJCVPERMUTE for details of  COMBCELLIDX, 
%   COMBCLASS, NFOLD and OPTION. 
%    
%   [CELLTRAINSEL,CELLTESTSEL] = TZ_OBJ2CELLPERM(...) also returns
%   indices of testing set.
%
%   See also TZ_CELL2OBJPERM

%   31-May-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

 if nargin < 5
    error('At least 5 arguments are required')
end 

nclass=max(combclass);

for i=1:nfold
    load([loadpath '/sel' num2str(i) 'fold.mat']);
    testobjclassidx=combclass(testsel==1);
    testobjcellidx=combcellidx(testsel==1);
    trainobjclassidx=combclass(trainsel==1);
    trainobjcellidx=combcellidx(trainsel==1);
    for j=1:nclass
        objcellidx=testobjcellidx(testobjclassidx==j);
        celltestsel{j}{i}=tz_getuniquenum(objcellidx);
        objcellidx=trainobjcellidx(trainobjclassidx==j);
        celltrainsel{j}{i}=tz_getuniquenum(objcellidx);
    end
end

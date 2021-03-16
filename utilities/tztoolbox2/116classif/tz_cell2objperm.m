function tz_cell2objperm(celltestsel,combcellidx,combclass,option,savepath)
%TZ_OBJCVPERMUTE Convert cell permutation to object permutation
%   TZ_OBJCVPERMUTE(CELLTESTSEL,COMBCELLIDX,COMBCLASS,OPTION,SAVEPATH)
%   converts testing indices of cell cross validation into object cross
%   validation. CELLTESTSEL is a cell array of cell array. 
%   CELLTESTSEL{J}{I} returns a vector of indices in class J and fold I.
%   COMBCELLIDX and COMBCLASS are column vectors of cell indices and 
%   class labels repectively. They must have the same number of rows.
%   OPTINS is the same as OPTION in TZ_OBJCVPERMUTE. The results will
%   be saved under SAVEPATH. This is the same as those from TZ_OBJCVPERMUTE.
%
%   See also TZ_OBJ2CELLPERM

%   30-May-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 5
    error('At least 5 arguments are required')
end

global glMtpath;

nclass=max(combclass);
nobj=length(combclass);

if isempty(savepath)
    savepath='results';
end

load permid.mat
createtime=clock;

for i=1:length(celltestsel{1})
    foldname=['sel' num2str(i) 'fold' option '.mat']
    
    testsel=zeros(nobj,1);
    
    for j=1:nclass
        for k=celltestsel{j}{i}
            testsel=testsel+((combclass==j)&(combcellidx==k));
        end
    end
    
    trainsel=ones(nobj,1)-testsel;
    
    savefile=[savepath '/' foldname];
    save(savefile,'trainsel','testsel','permid','createtime');
end

permid=permid+1;
save([glMtpath '/permid.mat'],'permid');

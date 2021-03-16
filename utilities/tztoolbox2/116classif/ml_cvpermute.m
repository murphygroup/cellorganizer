function ml_cvpermute(combobj,combcellidx,combclass,nfold,savepath, ...
    randstate)
%ML_CVPERMUTE Permuates for cross validation of object level classification
%   ML_CVPERMUTE(COMBOBJ,COMBCELLIDX,COMBCLASS,NFOLD,SAVEPATH) generates a
%   mat data file under SAVEPATH. The file contains indices of training 
%   sets and testing sets for NFOLD-fold cross validation. It also contains
%   permutation ID and created time. COMBOBJ is the [feature matrix] of all
%   objects. COMBCELLIDX is the column vector of cell indices and COMBCLASS
%   is the column vector of class labels.
%   
%   ML_CVPERMUTE(COMBOBJ,COMBCELLIDX,COMBCLASS,NFOLD,SAVEPATH,RANDSTATE)
%   uses RANDSTATE to set up random state for permutation. See RAND for
%   details: RAND('state',RANDSTATE).
%   
%   See also

%   09-Feb-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU


oldseed=rand('state');

if ~exist('randstate','var')
    randstate=[];
end

if ~isempty(randstate)
    rand('state',randstate);
end

k=1;
nclass=max(combclass);
nobj=size(combobj,1);

for i=1:nclass
    norgset(i)=max(combcellidx(combclass==i));
    perm{i}=randperm(norgset(i));
end

for i=1:nclass
    permfold(i,:)=[0,cumsum(ml_redistrnum(norgset(i),nfold))];
end

load permid.mat
createtime=clock;

for i=1:nfold
    foldname=['sel' num2str(i) 'fold' '.mat']
    testobjset=[];
    trainobjset=[];
    testsel=zeros(nobj,1);

    if isempty(savepath)
        savepath='results';
    end
    
    %generating training set and testing set
    
    %generate testing set
    for j=1:nclass
        for k=perm{j}(permfold(j,i)+1:permfold(j,i+1))
            %select class j and cell k
            testsel=testsel+((combclass==j)&(combcellidx==k));
        end
    end
    
    trainsel=ones(nobj,1)-testsel;
    
    savefile=[savepath '/' foldname];
    save(savefile,'trainsel','testsel','permid','createtime');
end

permid=permid+1;

save('/home/tpeng/matlab/permid.mat','permid');

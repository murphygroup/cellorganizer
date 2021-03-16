function tz_objcvpermute(combobj,combcellidx,combclass,nfold,option, ...
    savepath,randstate)
%TZ_OBJCVPERMUTE Generate cross validation permutation for objects.
%   TZ_OBJCVPERMUTE(COMBOBJ,COMBCELLIDX,COMBCLASS,NFOLD,OPTION,SAVEPATH)
%   generates permutations for object level analysis. COMBOBJ is the
%   feature matrix of objects. COMBCELLIDX is a column vector of indices
%   indicating which cell an object is located. COMBCLASS is a column
%   vector of class labels of objects. COMOBJ, COMBCELLIDX AND COMBCLASS
%   must have the same number of rows. NFOLD is the number of folds.
%   OPTION is for some special processiong. Currently there is only
%   one option called 'merge34', which does merging class 3 and 4.
%   COMBOBJ and OPTION might be removed in the future version. The
%   permutation results will be saved into files under SAVEPATH
%   with permutaion id.
%   
%   TZ_OBJCVPERMUTE(COMBOBJ,COMBCELLIDX,COMBCLASS,NFOLD,OPTION,SAVEPATH,
%   RANDSTATE) also lets users specifys the random state so that they
%   can repeat their results.

%   28-Aug-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 6
    error('6 or 7 arguments are required')
end

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

%%%%%%%%%%%%%%%merge class 3 4%%%%%%%%%%%%%%%%%
if strcmp(option,'merge34')
    nclass3=max(combcellidx(combclass==3));
    combcellidx(combclass==4)=combcellidx(combclass==4)+nclass3;
    combclass(combclass==4)=3;
    
    for i=5:nclass
        combclass(combclass==i)=i-1;
    end
    nclass=nclass-1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%remove class 3%%%%%%%%%%%%%%%%%
if strcmp(option,'remove3')
    allcombcellidx=combcellidx;
    allcombobj=combobj;
    allcombclass=combclass;
    newcombcellidx=combcellidx(combclass==3);
    newcombobj=combobj(combclass==3,:);
    newcombclass=combclass(combclass==3);
    
    combcellidx(combclass==3)=[];
    combobj(combclass==3,:)=[];
    combclass(combclass==3)=[];
    
    for i=4:nclass
        combclass(combclass==i)=i-1;
    end
    nclass=nclass-1;
    nobj=size(combobj,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:10
    norgset(i)=max(combcellidx(combclass==i));
    perm{i}=randperm(norgset(i));
end

for i=1:10
    permfold(i,:)=[0,cumsum(tz_redistrnum(norgset(i),nfold))];
end

load permid.mat
createtime=clock;

for i=1:nfold
    foldname=['sel' num2str(i) 'fold' option '.mat']
    testobjset=[];
    trainobjset=[];
    testsel=zeros(nobj,1);

    if isempty(savepath)
        savepath='results';
    end
    
    %generating training set and testing set
    for j=1:nclass
        for k=perm{j}(permfold(j,i)+1:permfold(j,i+1))
            testsel=testsel+((combclass==j)&(combcellidx==k));
        end
    end
    
    trainsel=ones(nobj,1)-testsel;
    
    savefile=[savepath '/' foldname];
    save(savefile,'trainsel','testsel','permid','createtime');
end

permid=permid+1;

save('/home/tingz/matlab/permid.mat','permid');

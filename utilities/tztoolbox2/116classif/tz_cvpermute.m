function tz_cvpermute(combobj,combcellidx,combclass,nfold,option,savepath,randstate)
%TZ_CVPERMUTE Obsolete.

%function tz_cvpermute(combobj,combcellidx,combclass,nfold,option,savepath,randstate)
%
%OVERVIEW:
%   create permutation for cross validation

error('Function tz_cvpermute is out of date. Please use tz_objcvpermute.');

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

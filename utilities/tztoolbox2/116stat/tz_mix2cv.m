function b=tz_mix2cv(combobj,combcellidx,combclass,nfold,clstmeth, ...
    clustk,classifmeth1,t1,option,loadpath,other)

%function b=tz_mix2cv(combobj,combcellidx,combclass,nfold,clstmeth,clustk,classifmeth1,t1,option,loadpath,other)
%
%OVERVIEW:
%   cross validation decompose 2-mixture 
%PARAMETERS:
%   combobj - combined objects
%   combcellidx - combined cell idx
%   combclass - combined class label
%   nfold - n fold
%   clstmeth - clustering method
%   clustk - cluster number
%   classifmeth1 - first level classifier
%   t1 - parameters for first level classifer
%   option - option of reducing features
%   loadpath - directory of loading clustering results
%   other - another option
%RETURN:
%   b - decomposing results, cell array, each element is a matrix,
%DESCRIPTION:
%   
%HISTORY:
%   ??-???-2004 Initial write TINGZ
%   05-NOV-2004 Modified TINGZ
%       - add comments

clustername=[clstmeth option other];


nclass=max(combclass);
nobj=size(combobj,1);

clustername=[clstmeth option other];

%%%%%%%%%%%%%%%remove class 3%%%%%%%%%%%%%%%%%
if strcmp(option,'remove3')
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

%extract fluorescence fraction
%flindex=combobj(:,4);

% select features
combobj=combobj(:,[1 3 5:end]);

for kk=clustk
    for fold=1:nfold
        foldname=['sel' num2str(fold) 'fold' option '.mat']

        testobjset=[];
        trainobjset=[];
        testsel=zeros(nobj,1);


        clusterfilename=[num2str(kk) clustername num2str(fold) 'fold.mat'];
        %clustering
        if isempty(loadpath)
            error('Empty loadpath');
        else
            load([loadpath '/' num2str(kk) clustername num2str(fold) 'fold.mat'])
            [trainnorm,testnorm]=mb_featurenorm(combobj(find(trainsel),:),combobj(find(testsel),:));
        end
        
        if (strcmp(clstmeth,'kmeans') | strcmp(clstmeth,'kmeansmahal')) & kk~=1
%             trainset=tz_cellobjcom(combclass(find(trainsel)),combcellidx(find(trainsel)),trainpost);
            [m,objidcs]=max(trainpost');
            objidcs=objidcs';
        end        
        
        %reassign labels
        if ~strcmp(classifmeth1,'classcluster') & kk==1
            objidcs=ones(sum(trainsel),1);
            allpost=ones(nobj,1);
        else
            if (strcmp(classifmeth1,'lknn') & t1==1) | (strcmp(clstmeth,'kmeans') & strcmp(classifmeth1,'dist')) | strcmp(other,'leavetrain')
                allpost=[objidcs;tz_classify(testnorm,trainnorm,objidcs,classifmeth1,0,t1)];
            else
                allpost=tz_classify([trainnorm;testnorm],trainnorm,objidcs,classifmeth1,0,t1);
            end
        end
        
       
        ncluster=kk;
        
        if strcmp(clstmeth,'classcluster')
            ncluster=max(objidcs);
        end
        
        ncluster
        
        %Calculate cell features
        trainset=tz_cellobjcom(combclass(find(trainsel)),combcellidx(find(trainsel)),allpost(1:length(objidcs)),ncluster,'objnum');
        testset=tz_cellobjcom(combclass(find(testsel)),combcellidx(find(testsel)),allpost((length(objidcs)+1):end),ncluster,'objnum');
        
        [objsq,objnum]=tz_trainobjcom(trainset,[]);
        
        
        k=1;
        for c1=1:nclass
            for c2=(c1+1):nclass
                mixp{k}=[c1 c2;sum(objnum{c1}(:,2)),sum(objnum{c2}(:,2))];
                k=k+1;
            end
        end
        b{fold}=[];
        for k=1:length(mixp)
            
           
            
            m=1;
            cperm=mixp{k}(1,:);
            
            for i=1:size(testset{cperm(1)},1)
                for j=1:size(testset{cperm(2)},1)
                    cperm=mixp{k}(1,:);
                    
                    logl=tz_maxloglmnom2(objsq,objnum,testset{cperm(1)}(i,:)+testset{cperm(2)}(j,:),mixp);
                    [a,b{fold}(k,m)]=max(logl);
                    m=m+1;
                end
                i
            end
        end
       
    end
end
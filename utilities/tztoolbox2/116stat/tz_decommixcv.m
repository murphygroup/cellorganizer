function [ress,mress]=tz_decommixcv(combobj,weights,combcellidx, ...
    combclass,nfold,clstmeth,clustk,clstsda,classifmeth1,t1,workdir, ...
    other,ntrial,classnum,cellnums,evalmeths,fitmeth,minalpha,estnum)
%TZ_DECOMMIXCV Evaluate mixture models by cross validation
%   RESS = TZ_DECOMMIXCV(COMBOBJ,WEIGHTS,COMBCELLIDX,COMBCLASS,NFOLD,
%   CLSTMETH,CLUSTK,CLSTSDA,CLASSIFMETH1,T1,WORKDIR,OTHER,NTRIAL,CLASSNUM,
%   CELLNUMS,EVALMETHS,FITMETH,MINALPHA,ESTNUM)
%   
%   [RESS,MRESS] = TZ_DECOMMIXCV(...)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function res=tz_evalmix(combobj,combcellidx,combclass,nfold,clstmeth,clustk,classifmeth1,t1,option,loadpath,other,ntrial,evalmeth,fitmeth)
%
% function [avgcm,avgacc,aic]=tz_objclassif(combobj,combcellidx,combclass,combcoords,combobjects,
%fluofrac,cofdist,...
%     nfold,clstmeth,clustk,clstsda,classifmeth1,t1,classifmeth2,preprocess,t2,featset,workdir,redo,...
%     other)

%OVERVIEW:
%   Evaluate the goodness-of-fit of mixture models by cross validation
%PARAMETERS:
%   same as tz_objclassif except
%   weights - 
%   evalmeth - evaluation method
%   fitmeth - mixture fitting method
%RETURN:
%   res - overall error
%DESCRIPTION:
%
%HISTORY:
%   28-MAY-2004 Initial write Tingz
%   31-MAY-2004 Modified TINGZ
%       - add evalmeth
%   14-JUN-2004 Modified TINGZ
%       - add fitmeth
%   04-NOV-2004 Modified TINGZ
%       - revise comments
%       - add 'per' evalmeth

if ~exist('fitmeth','var')
    fitmeth='Hs';
end

% clustername=[clstmeth num2str(clstsda) option other];


nclass=max(combclass);
nobj=size(combobj,1);

%for saving clustering results
clustername=[clstmeth num2str(clstsda) other{1}];

% featset='#$objnum';
if isempty(weights)
  featset='#!objnum';
else
  featset='#!featsum';
end

%directory for saving object type learning results
objtypedirname=[classifmeth1 tz_cell2str(t1,[]) other{2}]

%directory for saving cell level feature results
cellfeatdirname=['featset_' featset other{3}]
    
%directory for the 2nd classification results
classifdirname=['2nd' classifmeth2 num2str(preprocess) tz_cell2str(t2,[]) other{4}]

for kk=clustk
    inc=1;
    %directory for saving clustering results
    clusterdirname=[num2str(kk) clustername];
    clusterpath=[workdir '/' clusterdirname];
    objtypepath=[clusterpath '/' objtypedirname];
    cellfeatpath=[objtypepath '/' cellfeatdirname];
    
    for fold=1:nfold
        res(fold)=0;
%         foldname=['sel' num2str(fold) 'fold' option '.mat']
        foldname=['sel' num2str(i) 'fold' '.mat']
        
        testobjset=[];
        trainobjset=[];
        testsel=zeros(nobj,1);
        
        %number of clusters + [cluster_method sda_option option other] + fold + 'fold.mat'
        clusterfilename=[num2str(kk) clustername num2str(i) 'fold.mat'];
        
        %learned object type file name
        objtypefilename=[classifmeth1 tz_cell2str(t1,[]) clusterfilename];
        
        %cell level feature file name
        cellfeatfilename=[featset objtypefilename];
                
        load([workdir '/' foldname]);
        
%         clusterfilename=[num2str(kk) clustername num2str(fold) 'fold.mat'];
%         objtypefilename=[classifmeth1 num2str(t1) clusterfilename];
        objnumfilename=['#!objnum' objtypefilename];
        
        if isempty(weights)
            cellfeatfilename=[featset objtypefilename clusterfilename];
        else
            cellfeatfilename=['w' num2str(round(mean(weights))) featset objtypefilename clusterfilename];
        end
        
        %clustering
        if isempty(loadpath)
            error('Empty loadpath');
        else
            load([loadpath '/' clusterfilename])
            if exist('clstfeatsel','var')
                if ~isempty(clstfeatsel)
                    combobj=combobj(:,clstfeatsel);
                end
            end
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
            if exist([loadpath '/' objtypefilename],'file')
                load([loadpath '/' objtypefilename]);
            else
                if (strcmp(classifmeth1,'lknn') & t1==1) | (strcmp(clstmeth,'kmeans') & strcmp(classifmeth1,'dist')) | strcmp(other,'leavetrain')
                    allpost=[objidcs;tz_classify(testnorm,trainnorm,objidcs,classifmeth1,0,t1)];
                else
                    allpost=tz_classify([trainnorm;testnorm],trainnorm,objidcs,classifmeth1,0,t1);
                end
            end
        end
        
        
        ncluster=kk;
        
        if strcmp(clstmeth,'classcluster')
            ncluster=max(objidcs);
        end
        
        ncluster
        

        if estnum==1
            load([loadpath '/' objnumfilename]);
            trainset2=trainset;
            testset2=testset;
            
        end
        
        if exist([loadpath '/' cellfeatfilename],'file')
            load([loadpath '/' cellfeatfilename]);
        else
            %Calculate cell features
            if ~isempty(weights)
                trainset=tz_cellobjcom(combclass(find(trainsel)),combcellidx(find(trainsel)),allpost(1:length(objidcs)),ncluster,'#$featsum',weights(find(trainsel),:));
                testset=tz_cellobjcom(combclass(find(testsel)),combcellidx(find(testsel)),allpost((length(objidcs)+1):end),ncluster,'#$featsum',weights(find(testsel),:));
            else
                trainset=tz_cellobjcom(combclass(find(trainsel)),combcellidx(find(trainsel)),allpost(1:length(objidcs)),ncluster,'#$objnum',weights(find(trainsel),:));
                testset=tz_cellobjcom(combclass(find(testsel)),combcellidx(find(testsel)),allpost((length(objidcs)+1):end),ncluster,'#$objnum',weights(find(testsel),:));
                
            end
            save([loadpath '/' cellfeatfilename],'trainset','testset');
        end
        
        
        [objsq,objnum]=tz_trainobjcom(trainset,[]);
        p=tz_normobjcom(objsq);
        pw=sum(objsq,2);
        
        if estnum==1
            [objsq2,objnum2]=tz_trainobjcom(trainset2,[]);
            p2=tz_normobjcom(objsq2);
            pw2=sum(objsq2,2);
        end
        tmptrainset=trainset;
        tmptestset=testset;
        
        for cellnum=cellnums
            
            for ne=1:length(evalmeths)
                evalmeth=evalmeths{ne};
                trainset=tmptrainset;
                testset=tmptestset;
                switch evalmeth
                case 'mp2'
                    trainset{3}=[trainset{3};trainset{4}];
                    trainset(4)=[];
                    testset{3}=[testset{3};testset{4}];
                    testset(4)=[];
                    [objsq,objnum]=tz_trainobjcom(trainset,[]);
                    p=tz_normobjcom(objsq);
                    pw=sum(objsq,2);
                case 'mpe'
                    trainset{3}=[trainset{3};trainset{4}];
                    trainset{5}=[trainset{5};trainset{9}];
                    trainset(4)=[];
                    trainset(8)=[];
                    testset{3}=[testset{3};testset{4}];
                    testset{5}=[testset{5};testset{9}];
                    testset(4)=[];
                    testset(8)=[];
                    [objsq,objnum]=tz_trainobjcom(trainset,[]);
                    p=tz_normobjcom(objsq);
                    pw=sum(objsq,2);
                end
                nclass=length(trainset);    
                res=0;
                for i=1:ntrial
                    scell=[];
                    scell2=[];
                    %             sc=binornd(1,0.5,1,nclass);
                    %             while(sum(sc)==0)
                    %                 sc=binornd(1,0.5,1,nclass);
                    %             end
                    if classnum<=0
                        snclass=randperm(nclass);
                    else
                        snclass=classnum;
                    end
                    
                    sclass=randperm(nclass);
                    
                    
                    sc=zeros(1,nclass);
                    sc(sclass(1:snclass(1)))=1;
                    
                    for j=1:nclass
                        
                        if sc(j)==1
                            perm=randperm(size(testset{j},1));
                            tmpcell=testset{j}(perm(1),:);
                            if estnum==1
                                tmpcell2=testset2{j}(perm(1),:);
                            end
                            for k=2:cellnum
                                tmpcell=tmpcell+testset{j}(perm(k),:);
                                if estnum
                                    tmpcell2=tmpcell2+testset2{j}(perm(k),:);
                                end
                            end
                            scell=[scell;tmpcell];
                            if estnum
                                scell2=[scell2;tmpcell2];
                            end
                        end
                        
                    end
                    y=sum(scell,1);
                    if estnum
                        y2=sum(scell2,1);
                    end
                    %             [alpha,selmodel,aic]=tz_aicfitmnomk(y,p);
                    %             ealpha=zeros(1,nclass);
                    %             ealpha(selmodel)=alpha;
                    talpha=zeros(1,nclass);
                    
                    %                     talpha1=sum(scell,2)'/sum(y);
                    
                    talpha1=talpha;
                    talpha1(find(sc))=cellnum;
                    talpha1=talpha1/sum(talpha1);
                    
                    talpha2=talpha;
                    talpha2(find(sc))=sum(scell,2)';
                    talpha2=talpha2/sum(talpha2);
                    
                    %                     if isempty(weights)
                    %                         talpha(find(sc))=sum(scell,2)'/sum(y);
                    %                     else
                    %                         talpha(find(sc))=cellnum;
                    %                         talpha=talpha/sum(talpha);
                    %                     end
                    
                    switch fitmeth
                    case 'Hx'
                        ealpha=tz_fitmnomk(y2,p2,minalpha);
                        ealpha2=tz_fitmnomk(y,p(ealpha>0,:),0);
                        ealpha(ealpha>0)=ealpha2;
                        talpha=talpha2;
                    case 'Hf'
                        succ=0;
                        while succ==0
                            [ealpha,loglk,succ]=tz_fitmnomk_fix(y,p,snclass(1),[],minalpha);
                        end
                        talpha=talpha2;
                    case 'Hs'
                        succ=0;
                        while succ==0
                            [ealpha,loglk,succ]=tz_fitmnomk(y,p,minalpha);
                        end
                        talpha=talpha2;
                    case 'Hm'  
                        succ=0;
                        while succ==0
                            [ealpha,loglk,succ]=tz_fitmnomk(y*10/max(y),p,minalpha);
                        end
                        talpha=talpha2;
                        
                    case 'Hn'
                        sobj=sum(objsq,2);
                        [ealpha,loglk,succ]=tz_fitmnomk(y,p,minalpha);
                        ealpha=ealpha./sobj';
                        ealpha=ealpha/sum(ealpha);
                        talpha=talpha1;
%                         talpha(find(sc))=sum(scell,2)'/sum(y);
                    case 'Ht'
                        succ=0;
                        while succ==0
                            [ealpha,loglk,succ]=tz_fitmnomk(y,p,minalpha);
                        end
                        talpha=talpha2;
                        
                    case 'Rd'
                        [fitx,logll]=tz_redistrmnom(y,p);
                        ealpha=sum(fitx,2);
                        ealpha=ealpha/sum(ealpha);
                        ealpha=ealpha';
                        talpha=talpha2;
%                         talpha(find(sc))=sum(scell,2)'/sum(y);
                    case 'Ln'
                        ealpha=tz_solvemix(y,objsq);
                        talpha=talpha1;
                    case 'Lo'
                        ealpha=tz_solvemix(y,p);
                        talpha=talpha2;
                        %                         talpha(find(sc))=sum(scell,2)'/sum(y);
                    case 'Lf'
                        
                        ealpha=tz_solvemix_fix(y,p,snclass(1),[]);
                        
                        talpha=talpha2;
                    end
                    
                    
                    ealpha
                    talpha
                    switch evalmeth
                    case 'sig'
                        err=tz_evalmixsig(ealpha,talpha,1000);
                    case 'sig2'
                        err=tz_evalmixsig(ealpha,talpha,1000,'nonp');
                    case 'abs'
                        err=sum(abs(ealpha-talpha));
                    case 'sq2'
                        err=sum((ealpha-talpha).^2);
                    case 'per'
                        err=1-sum(min([ealpha;talpha]));
                    case 'mpe'
%                         ealpha(3)=ealpha(3)+ealpha(4);
%                         talpha(3)=talpha(3)+talpha(4);
%                         ealpha(4)=0;
%                         talpha(4)=0;
%                         ealpha(5)=ealpha(5)+ealpha(9);
%                         talpha(5)=talpha(5)+talpha(9);
%                         ealpha(9)=0;
%                         talpha(9)=0;
                        err=1-sum(min([ealpha;talpha]));
                    case 'mp2'
%                         ealpha(3)=ealpha(3)+ealpha(4);
%                         talpha(3)=talpha(3)+talpha(4);
%                         ealpha(4)=0;
%                         talpha(4)=0;
%                         
                        err=1-sum(min([ealpha;talpha]));
                    end
                    
                    res=res+err;
                end
                
                res=res/ntrial;
                ress{inc}=struct('res',res,'cellnum',cellnum,'evalmeth',evalmeth,'fold',fold);
                inc=inc+1;
            end
            
        end
    end
end
inc=1;
for i=1:nfold
    for j=1:length(cellnums)
        for k=1:length(evalmeths)
            mress(i,j,k)=ress{inc}.res;
            inc=inc+1;
        end
    end
end
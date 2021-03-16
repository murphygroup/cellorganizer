function [avgcm,avgacc,aic]=tz_testcv(combobj,combcellidx,combclass,combcoords,combobjects,nfold,...
    clstmeth,clustk,clstsda,classifmeth1,t1,classifmeth2,preprocess,t2,featset,option,loadpath,...
    permpath,savepath,other,comments,recalc)
%TZ_TESTCV Obsolete. Please use TZ_OBJCLASSIF intead.
% function [avgcm,avgacc,aic]=tz_testcv(combobj,combcellidx,combclass,combcoords,combobjects,nfold,...
%     clstmeth,clustk,clstsda,classifmeth1,t1,classifmeth2,preprocess,t2,featset,option,loadpath,...
%     permpath,savepath,other,comments,recalc)
%
%OVERVIEW:
%   stratified cross validation for object level classification
%PARAMETERS:
%   combobj - object-level features nobjxp
%   combcelledx - cell indices for objects nobjx1
%   combclass - class indices for objects nobjx1
%   nfold - nfold-fold cross validation
%   clstmeth - clustering method
%   clustk - cluster numbers
%   clstsda - feature selection before clustering
%       0 - no feature selection
%       <0 - all sda selected feautres
%       >0 - best clstsda features
%   classifmeth1 - first level classification
%   t1 - additional parameters for first level classification
%   classifmeth2 - second level classification
%   t2 - additional parameters for second level classification
%   featset - cell-level features
%   option - options for preprocessing
%       'merge34'
%       'remove3'
%   loadpath - loading path for clustering results
%   permpath - loading path for cv partitions
%   savepath - path for saving results
%   other - some weired options
%   comments - comments for saved results
%   recalc - recalculate object types and cell level features or not
%DESCRIPTION
%   
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   08-APR-2004 Modified TINGZ
%   18-APR-2004 Modified TINGZ
%       - Add classification methods option
%   19-APR-2004 Modified TINGZ
%       - Add 1st and 2nd classification methods option
%   16-MAY-2004 Modified TINGZ
%       - add savepath and save comments
%   18-MAY-2004 Modified TINGZ
%       - Support nonexist files
%   23-MAY-2004 Modified TINGZ
%       - record while saving files
%   14-NOV-2004 Modified TINGZ
%       - remove feature selection
%   29-NOV-2004 Modified TINGZ
%       - give options for feature selection
%   29-NOV-2004 Modified TINGZ
%       - dist classifier

oldseed=rand('state');

%initialize random state to make the results repeatable
rand('state',0);

%for saving classification results
%if isempty(t3)
setname=['cluster' clstmeth num2str(clstsda) '1st' classifmeth1 num2str(t1) '2nd' classifmeth2 num2str(preprocess) num2str(t2) featset option other]; 
    %else
%     setname=['cluster' clstmeth num2str(clstsda) '1st' classifmeth1 num2str(t1) '2nd' classifmeth2 num2str(preprocess) num2str(t2) 't2' num2str(t3) featset option other]; 
% end

%for saving clustering results
clustername=[clstmeth num2str(clstsda) option other];

%nfold=10;
k=1;
nclass=max(combclass);
nobj=size(combobj,1);

%Do not normalize the feature if it's fluorescence fraction
% if strcmp(featset,'fluofrac')% | strcmp(featset,'objnum')
%     preprocess=1;
% else
%     preprocess=1;
% end

if isempty(preprocess)
    preprocess=1;
end

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
    combcellidx(combclass==3)=[];
    combobj(combclass==3,:)=[];
    if ~isempty(combcoords)
        combcoords(combclass==3,:)=[];
    end
    
    combclass(combclass==3)=[];
    

    for i=4:nclass
        combclass(combclass==i)=i-1;
    end
    nclass=nclass-1;
    nobj=size(combobj,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% select features
if isempty(other)
    %extract fluorescence fraction
    flindex=combobj(:,4);
    cdindex=combobj(:,2);
    combobj=combobj(:,[1 3 5:end]);
else
    if other(1)~='%'
        %extract fluorescence fraction
        flindex=combobj(:,4);
        cdindex=combobj(:,2);
        combobj=combobj(:,[1 3 5:end]);
    else
        flindex=ones(size(combobj,1),1);
        cdindex=ones(size(combobj,1),1);
    end
end

%%%%%%%%%%%%%%%%%%n fold permutation%%%%%%%%%%%%%%
% if isempty(loadpath)
%     if isempty(permpath)
        for i=1:nclass
            norgset(i)=max(combcellidx(combclass==i));
            perm{i}=randperm(norgset(i));
        end
        
        for i=1:nclass
            permfold(i,:)=[0,cumsum(tz_redistrnum(norgset(i),nfold))];
        end
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for kk=clustk
   
    for i=1:nfold
        foldname=['sel' num2str(i) 'fold' option '.mat']

        testobjset=[];
        trainobjset=[];
        testsel=zeros(nobj,1);

%         if strcmp(clstmeth,'kmeansmahal')
%             % Seed the random number generator
%             options(5) = 1; % let seeds be picked randomly by rm_mahalkmeans
%             seeds = ones(kk,size(trainnorm,2));
%         end
        if isempty(savepath)
            savepath='results';
        end
        
        %number of clusters + [cluster_method sda_option option other] + fold + 'fold.mat'
        clusterfilename=[num2str(kk) clustername num2str(i) 'fold.mat'];
        
        %learned object type file name
        objtypefilename=[classifmeth1 num2str(t1) clusterfilename];
        
        %cell level feature file name
        cellfeatfilename=[featset objtypefilename clusterfilename];
        
        newclustering=1;
        
        if ~isempty(loadpath)
            loadfile=[loadpath '/' clusterfilename]
            newclustering=0;
            if(isempty(dir(loadfile)))
                warning('specified clustering file not found. New clustering will be done');
                newclustering=1;
            end
        end
             
        %clustering
        if newclustering==1    
            if isempty(permpath)
                %generating training set and testing set
                for j=1:nclass
                    for k=perm{j}(permfold(j,i)+1:permfold(j,i+1))
                        testsel=testsel+((combclass==j)&(combcellidx==k));
                    end
                end
                
                trainsel=ones(nobj,1)-testsel;
                
                savefile=[savepath '/' foldname];
                save(savefile,'trainsel','testsel','comments');
                tz_updaterecord([savefile ' saved'],'results');
                
            else
                load([permpath '/' foldname]);
            end
            
            %feature selection
            if ~isempty(clstsda)
                if clstsda~=0
                    clstfeatsel=tz_sda(combobj(find(trainsel),:),combclass(find(trainsel),:));
                    if clstsda<0
                        combobj=combobj(:,clstfeatsel); 
                    else
                        combobj=combobj(:,clstfeatsel(1:clstsda));
                    end
                else
                    clstfeatsel=[];
                end
            else
                clstfeatsel=[];
            end
            
            [trainnorm,testnorm]=ml_featurenorm(combobj(find(trainsel),:),combobj(find(testsel),:));
            if strcmp(clstmeth,'kmeans') | strcmp(clstmeth,'kmeansmahal') |strcmp(clstmeth,'tzkmeans')
                if strcmp(other,'35randseeds')
                    rand('state',[1:35]'+2);
                    randperm(size(trainnorm,1));
                else
                    rand('state',0);
                end
                
                randidx = randperm(size(trainnorm,1));
                seeds = trainnorm(randidx(1:kk),:);
                
                % k-means clustering options
                options = zeros(1,14);
            end
            
            
            
            disp('clustering...')
            
            %start clustering
            switch clstmeth
            case 'kmeans'
                if kk~=1
                    [centers,options,trainpost,errlog]=netlab_kmeans(seeds,trainnorm,options);
                    savefile=[savepath '/' clusterfilename];
                    save(savefile,'centers','clstfeatsel','trainpost','trainsel','testsel','comments');
                    tz_updaterecord([savefile ' saved'],'results');
                else
                    trainpost=ones(sum(trainsel),1);
                end
            case 'tzkmeans' %remove small clsuters
                if kk~=1
                    [centers,options,trainpost,errlog]=tz_kmeans(seeds,trainnorm,options,-1);
                    savefile=[savepath '/' clusterfilename];
                    save(savefile,'centers','clstfeatsel','trainpost','trainsel','testsel','comments');
                    tz_updaterecord([savefile ' saved'],'results');
                else
                    trainpost=ones(sum(trainsel),1);
                end
            case 'kmeansmahal'
                if kk~=1
                    options(5)=1;
                    [centers,options,trainpost,errlog]=rm_mahalkmeans(seeds,trainnorm,options);
                    savefile=[savepath '/' clusterfilename];
                    save(savefile,'centers','clstfeatsel','trainpost','trainsel','testsel','comments');
                    tz_updaterecord([savefile ' saved'],'results');
                else
                    trainpost=ones(sum(trainsel),1);
                end
            case 'classcluster' %classwise clustering
                if kk==1
                    objidcs=combclass(find(trainsel));
                else
                    rand('state',0);
                    objidcs=tz_classcluster(trainnorm,combclass(find(trainsel)),kk,0);   
                end
                trainpost=tz_label2post(objidcs);
                savefile=[savepath '/' clusterfilename];
                save(savefile,'objidcs','clstfeatsel','trainpost','trainsel','testsel','comments');                
                tz_updaterecord([savefile ' saved'],'results');
            end
            disp('done!')
        else
            %load results if there exist clustering files
            load(loadfile)
            if ~isempty(clstsda)
                if clstsda~=0
                    if ~isempty(clstfeatsel)
                        if clstsda<0
                            combobj=combobj(:,clstfeatsel); 
                        else
                            combobj=combobj(:,clstfeatsel(1:clstsda));
                        end
                        
                    end
                end
            else
                clstfeatsel=[];
            end
            [trainnorm,testnorm]=ml_featurenorm(combobj(find(trainsel),:),combobj(find(testsel),:));
        end
        
        if (strcmp(clstmeth,'kmeans') | strcmp(clstmeth,'kmeansmahal') | strcmp(clstmeth,'tzkmeans')) & kk~=1
%             trainset=tz_cellobjcom(combclass(find(trainsel)),combcellidx(find(trainsel)),trainpost);
            [m,objidcs]=max(trainpost');
            objidcs=objidcs';
        end        
        
        %reassign labels
        if ~strcmp(clstmeth,'classcluster') & kk==1
            objidcs=ones(sum(trainsel),1);
            allpost=ones(nobj,1);
        else

            typelearn=1;
            if ~isempty(loadpath)  & recalc==0
                if exist([loadpath '/' objtypefilename],'file')
                    load([loadpath '/' objtypefilename]);
                    typelearn=0;
                end
            end
            
            if typelearn==1
                if (strcmp(classifmeth1,'lknn') & t1==1) | strcmp(other,'leavetrain')
                    allpost=[objidcs;tz_classify(testnorm,trainnorm,objidcs,classifmeth1,0,t1)];
                else
                    allpost=tz_classify([trainnorm;testnorm],trainnorm,objidcs,classifmeth1,0,t1);
                end
                save([savepath '/' objtypefilename],'allpost');
            end
            %             end
        end
        
        [aic(i), clusterinfo] = rm_akaike_euclid( trainnorm, trainpost);
        
        ncluster=kk;
        
        if strcmp(clstmeth,'classcluster')
            ncluster=max(objidcs);
        end
        
        ncluster
        if strcmp(classifmeth2,'nnobj')
            [training,testing]=ml_featurenorm(combobj(find(trainsel),:),combobj(find(testsel),:));
            traincell=[training,combcellidx(find(trainsel))];
            if t2==1
                testing=[testing,flindex(find(testsel))];
            end
            rlabels = combclass(find(trainsel));
            [testcell,tlabels]=tz_combfeat2cell(testing,combclass(find(testsel)),combcellidx(find(testsel)));
        else
            %Calculate cell features
            if strcmp(featset,'newset')
                trainset=tz_cellobjfeats_clst(combclass(find(trainsel)),combcellidx(find(trainsel)),combobj(find(trainsel),:),allpost(1:length(objidcs)),ncluster);
                testset=tz_cellobjfeats_clst(combclass(find(testsel)),combcellidx(find(testsel)),combobj(find(testsel),:),allpost((length(objidcs)+1):end),ncluster);
            else
                if ~isempty(combcoords)
                    trcombcoords=combcoords(find(trainsel),:);
                    tscombcoords=combcoords(find(testsel),:);
                else
                    trcombcoords=[];
                    tscombcoords=[];
                end
                
                if ~isempty(combobjects)
                    trcombobjects=combobjects(find(trainsel));
                    tscombobjects=combobjects(find(testsel));
                else
                    trcombobjects={};
                    tscombobjects={};
                end
                featcal=1;
                if ~isempty(loadpath) & recalc==0
                    if exist([loadpath '/' cellfeatfilename],'file')
                        load([loadpath '/' cellfeatfilename]);
                        featcal=0;
                    end
                end
                
                if featcal==1
                    trainset=tz_cellobjcom(combclass(find(trainsel)),combcellidx(find(trainsel)),allpost(1:length(objidcs)),ncluster,featset,combobj(find(trainsel),:),flindex(find(trainsel)),cdindex(find(trainsel)),trcombcoords,tscombobjects);
                    testset=tz_cellobjcom(combclass(find(testsel)),combcellidx(find(testsel)),allpost((length(objidcs)+1):end),ncluster,featset,combobj(find(testsel),:),flindex(find(testsel)),cdindex(find(testsel)),trcombcoords,tscombobjects);
                    save([savepath '/' cellfeatfilename],'trainset','testset');
                end
            end
            %Real labels for training set and testing set
            traincell=[];
            rlabels=[];
            testcell=[];
            tlabels=[];
            for m=1:nclass
                tt=testset{m};
                testcell=[testcell;tt];
                tlabels=[tlabels;ones(size(tt,1),1)*m];
                tt=trainset{m};
                traincell=[traincell;tt];
                rlabels=[rlabels;ones(size(tt,1),1)*m];          
            end
            
        end        
        
        %     objmodel=tz_normobjcom(objsq);
        %     %logl=tz_testlogl(objsq,testobj);
        %     mixobj=zeros(1,size(objmodel,2))
        %     for m=1:length(testset)
        %         mixobj=mixobj+sum(testset{m}(:,1:(end-1)),1);
        %     end
        %     [x,selmodel,aic]=tz_aicfitmnom(mixobj,objmodel)
        %        
        %     sum(x,2)
        %[a,b]=max(logl')
        
        %     cvcm1={};
        %     cvcm2={};
        
        %Classification
        %[cvcm1{i},avgacc1(i)]=tz_loglclassif(objsq,[],testobj,labels);
        %[cvcm2{i},avgacc2(i)]=tz_loglclassif(objsq,objnum,testobj,labels);
        
           
        if strcmp(classifmeth2,'comp')
            rejs=zeros(nclass,nclass);
            for ci=1:nclass
                ncell=sum(tlabels==ci);
                ccells=testcell(tlabels==ci,:);
                inc=1;
                %intra-class comparison
                for cj=1:ncell
                    if inc>10
                        break;
                    end
                    for ck=cj:ncell-cj
                        s1=ccells(ck,:);
                        s2=ccells(ck+cj,:);
                        pvalue=tz_chi2test2(s1,s2);
                        rejs(ci,ci)=rejs(ci,ci)+(pvalue<0.05);
                        inc=inc+1;
                        if inc>10
                            break;
                        end
                    end
                end
                %inter-class comparison
                for cj=(ci+1):nclass
                    ncell2=sum(tlabels==cj);
                    ccells2=testcell(tlabels==cj,:);
                    npair=ncell*ncell2;
                    mpair=reshape(1:npair,ncell,ncell2);
                    randpair=randperm(npair);
                    
                    for ck=1:10
                        [c1,c2]=find(mpair==randpair(ck))
                        
%                         if c2==10
%                             ncell2
%                         end
                        s1=ccells(c1,:);
                        s2=ccells2(c2,:);
                        pvalue=tz_chi2test2(s1,s2);
                        rejs(ci,cj)=rejs(ci,cj)+(pvalue<0.05);
                    end
                end
            end
        else
            tg=tz_classify(testcell,traincell,rlabels,classifmeth2,preprocess,t2,ncluster);
            if preprocess>=2
                load('/home/tingz/tmp/featsel.mat');
                sda{i}=featsel;
            else
                sda{i}=NaN;
            end
            [ncvcm{i},pcvcm{i},cvacc(i)]=tz_summaryclassif(tlabels,tg);
        end
    end
    
    if strcmp(classifmeth2,'comp')
        avgcm=rejs;
        avgacc=0;
        savefile=[savepath '/' 'results' num2str(kk) setname num2str(i) 'fold.mat'];
        save(savefile,'rejs','comments');
        tz_updaterecord([savefile ' saved'],'results');
    else
        [avgcm,avgacc]=tz_calcavgcm(ncvcm);
        
        savefile=[savepath '/' 'results' num2str(kk) setname num2str(i) 'fold.mat'];
        save(savefile,'ncvcm','pcvcm','cvacc','aic','sda','comments');
        tz_updaterecord([savefile ' saved'],'results');
    end
   
    display(setname);
    
%     rand('state',0);
end

rand('state',oldseed);
function [avgcm,avgacc,aic,ncm,kappa,ua,pa]=...
    tz_objclassif(combobj,combcellidx,combclass,combcoords,combobjects,...
    fluofrac,cofdist,nfold,pfold,clstmeth,clustk,clstsda,...
    classifmeth1,t1,classifmeth2,preprocess,t2,featset,workdir,redo,other)
%TZ_OBJCLASSIF Object-based classification.
%   AVGCM = TZ_OBJCLASSIF(COMBOBJ,COMBCELLIDX,COMBCLASS,COMBCOORDS,
%       COMBOBJECTS,FLUOFRAC,COFDIST,NFOLD,PFOLD,CLSTMETH,CLUSTK,CLSTSDA,
%       CLASSIFMETH1,T1,CLASSIFMETH2,PREPROCESS,T2,FEATSET,WORKDIR,REDO,
%       OTHER) returns confusion matrices of object level classification
%       by NFOLD-fold stratified cross validation. COMBOBJ is an nxp matrix 
%       of the object level features if there are n objects and p features.
%       COMBCELLIDX is a nx1 vector specifying which cell each object is
%       from and COMBCLASS is also a nx1 vector specifying the cell class
%       membership of each object. Both combcellidx and combclass must be
%       integer vectors and the integers are contiguous and starting from 1.
%       Only folds in the vector PFLOD will be executed if PFOLD is not
%       empty. This is especially for parallel processing. CLSTMETH
%       specifies the clustering method for learning object types:
%           'kmeans' - batch-kmeans clustering      
%           'tz_kmeans' - batch-kmeans, but small clusters will be removed
%           'kmeansmahal' - kmeans clustering with mahalanobosis dist.
%           'classcluster' - classwise clustering
%       CLUSTK is a vector of clustering numbers that will be tried. 
%       Notice that for classwise clustering the numbers are actually
%       CLUSTK*(number of class). If CLUSTK is -1, then the clustering
%       will search the optimal cluster number automatically. But before
%       clustering, some features will be selected based on SDA if 
%       CLUSTSDA is not empty or 0. If CLSTSDA<0, all significant features
%       determined by SDA will be used. If CLSTSDA>0, only the first
%       CLSTK SDA-sorted features will be used. After clustering, all
%       objects will be assigned a type by the first level classifier,
%       which is specified by CLASSIFMETH1. T1 is a cell array of the 
%       parameters for classification. See TZ_CLASSIFY for more details.
%       CLASSIFMETH2 specifies the second level classification, and
%       PREPROCESS and T2 are parameters for this classification. Also see
%       TZ_CLASSIFY for more details. FFEATSET is a string specifying 
%       the set of cell-level features. See TZ_CELLOBJCOM for more
%       details. WORKDIR is the directory for saving or loading results.
%       REDO is a number to specify how to redo calculation, the order is
%       [arguments,permute,clustering,1st classif,2nd classif]. OTHER t
%       akes some strange options. Currently it has 
%       [clustering,1st classifier,2nd classifier,feat].
%       The intermediate results will be saved in the disk with the 
%       structure:
%       workdir/ (specified)
%           |
%           cluster/
%                  |
%                  1st classif/
%                             |
%                             2nd classif/
%                                        |
%   It is a new version of TZ_TESTCV.
%
%   Notice:
%       If CLASSIFMETH2 is started with 'mix', the whole function will 
%       become mixture decomposition. T2 will be a cell array contating
%       {'ntrial','classnum','cellnums','evalmeths','fitmeth','minalpha'.
%       It is related to the structured parameter for the function
%       TZ_DECOMMIXTR. See TZ_DECOMMIXTR for more details.
%   
%   [AVGCM,AVGACC,AIC,NCM,KAPPA,UA,PA] = TZ_OBJCLASSIF(...) also returns
%   overall accuracy, AIC values, confusion matrices with number, Kappa,
%   user's accuracy and producer's accuracy.

%   ??-???-???? Initial write T.Zhao
%   08-APR-2004 Modified T. ZHAO
%   18-APR-2004 Modified T. ZHAO
%       - Add classification methods option
%   19-APR-2004 Modified T. ZHAO
%       - Add 1st and 2nd classification methods option
%   16-MAY-2004 Modified T. ZHAO
%       - add savepath and save comments
%   18-MAY-2004 Modified T. ZHAO
%       - Support nonexist files
%   23-MAY-2004 Modified T. ZHAO
%       - record while saving files
%   14-NOV-2004 Modified T. ZHAO
%       - remove feature selection
%   29-NOV-2004 Modified T. ZHAO
%       - give options for feature selection
%   29-NOV-2004 Modified T. ZHAO
%       - dist classifier
%   18-FEB-2004 Modified T. ZHAO
%       - updated from tz_testcv    
%   03-MAR-2004 Modified T. ZHAO 
%       - add 'Gl' for mixture 
%   Copyright (c) Murphy Lab, Carnegie Mellon University


if ~exist(workdir,'dir')
    mkdir(workdir);
end

if isempty(other)
    other = {[],[],[],[]};
else
    other={other{:},[],[],[],[]};
end

datasize=size(combobj);
nclass=max(combclass);
% k=1;
nobj=size(combobj,1);

%number of cells in each class
ncells = zeros(1,nclass);
for i=1:nclass
    ncells(i) = max(combcellidx(combclass==i));
end

if redo>0
    %check whether the data are matched
    if ~exist([workdir '/' 'arg.mat'],'file')
        save([workdir '/arg.mat'],'datasize','nclass','ncells','nfold');
        ml_cvpermute(combobj,combcellidx,combclass,nfold,workdir);
    else
        paras=load([workdir '/' 'arg.mat']);
        if ~all([datasize,nclass,ncells,nfold] == ...
                [paras.datasize,paras.nclass,paras.ncells,paras.nfold])
            disp([datasize,nclass,ncells,nfold]- ...
                [paras.datasize,paras.nclass,paras.ncells,paras.nfold]);
            error('data not matched');
        end
    end
else %new permutation
    %save input arguments
    save([workdir '/arg.mat'],'datasize','nclass','ncells','nfold');
    ml_cvpermute(combobj,combcellidx,combclass,nfold,workdir);
end

%for saving classification results
if isstruct(t2)
    t2 = orderfields(t2);
    t2cell = struct2cell(t2);
else
    t2cell = t2;
end

setname=['cluster' clstmeth num2str(clstsda) '1st' classifmeth1 ...
    tz_cell2str(t1,[])...
    '2nd' classifmeth2 num2str(preprocess) tz_cell2str(t2cell,[]) ...
    featset tz_cell2str(other,[])];

%for saving clustering results
clustername=[clstmeth num2str(clstsda) other{1}];


if isempty(fluofrac)
    fluofrac=ones(size(combobj,1),1);
end
if isempty(cofdist)
    cofdist=ones(size(combobj,1),1);
end

%directory for saving object type learning results
objtypedirname=[classifmeth1 tz_cell2str(t1,[]) other{2}]

%directory for saving cell level feature results
cellfeatdirname=['featset_' featset other{3}]
    
%directory for the 2nd classification results
classifdirname=['2nd' classifmeth2 num2str(preprocess) ...
    tz_cell2str(t2cell,[]) other{4}]

if isempty(pfold)
    pfold=1:nfold;
end
ki=1;
for kk=clustk
    %directory for saving clustering results
    clusterdirname=[num2str(kk) clustername];
    clusterpath=[workdir '/' clusterdirname];
    objtypepath=[clusterpath '/' objtypedirname];
    cellfeatpath=[objtypepath '/' cellfeatdirname];
    
    for i=pfold
        foldname=['sel' num2str(i) 'fold' '.mat']

        %intialize
%         testobjset=[];
%         trainobjset=[];
        testsel=zeros(nobj,1); 
        
        %number of clusters + [cluster_method sda_option option other] + fold + 'fold.mat'
        clusterfilename=[num2str(kk) clustername num2str(i) 'fold.mat'];
        
        %learned object type file name
        objtypefilename=[classifmeth1 tz_cell2str(t1,[]) clusterfilename];
        
        %cell level feature file name
        cellfeatfilename=[featset objtypefilename];
                
        load([workdir '/' foldname]);
        
        if ~exist(clusterpath,'dir')
            mkdir(workdir,clusterdirname);
        end
        
        if redo<=2 %do clustering
            newclustering=1;
        else
            newclustering=0;
             
            if ~exist([clusterpath '/' clusterfilename],'file')
                warning('specified clustering file not found. New clustering will be done');
                newclustering=1;
            end  
            
        end
        
        clstfeatsel = [];
        
        if newclustering==1
            if ~isempty(clstsda)
                clstfeatsel = ml_sda(trainobj,trainclass);
            end
        else
            load([clusterpath '/' clusterfilename])            
        end
        
        if ~isempty(clstsda)
            if clstsda~=0
                clstfeatsel = ml_sda(trainobj,trainclass);
                combobj=combobj(:,clstfeatsel);
            else
                combobj=combobj(:,clstfeatsel(1:clstsda));
            end
        end
        
        trainidx = find(trainsel);
        testidx = find(testsel);  
        
        if strcmp(other{1},'eval')            
            [trainobj,trainclass,evalobj,evalclass] = ml_splitset( ...
                [combobj(trainidx,:) combcellidx(trainidx,:) trainidx], ...
                combclass(trainidx,:));
            trainidx = trainobj(:,end);
            evalidx = evalobj(:,end);
            traincellidx = trainobj(:,end-1);
            evalcellidx = evalobj(:,end-1);            
            trainobj = trainobj(:,1:end-2);
            evalobj = evalobj(:,1:end-2);
        else
            trainobj = combobj(trainidx,:);
            evalobj = trainobj;
            trainclass = combclass(trainidx,:);
            evalclass = trainclass;
            evalidx = trainidx;
            traincellidx = combcellidx(trainidx,:);
            evalcellidx = traincellidx;
            evalidx = trainidx;
        end  
        testobj = combobj(testidx,:);
        testcellidx = combcellidx(testidx,:);
        testclass = combclass(testidx,:);

        [trainnorm,testnorm]=ml_featurenorm(trainobj,testobj);
        if strcmp(other{1},'eval')
            [trainnorm,evalnorm] = ml_featurenorm( ...
                trainobj,evalobj);
        else
            evalnorm = trainnorm;
        end
                
        %clustering
        if newclustering==1                              
            if strcmp(clstmeth,'kmeans') | ...
                    strcmp(clstmeth,'kmeansmahal') | ...
                    strcmp(clstmeth,'tzkmeans')               
                rand('state',0);
                
                randidx = randperm(size(trainnorm,1));
                
                if strcmp(other{1},'35randseeds')
%                     rand('state',[1:35]'+2);
%                     randperm(size(trainnorm,1));
                    randidx = 1:kk;
                else
                    if ~isempty(strmatch('classeeds',other{1}))
                        seedk = str2num(other{1}(10:end));
                        seedclassidx = find(trainclass== ...
                                            seedk);
                        randidx = randsample(seedclassidx, ...
                                             length(seedclassidx));
                        
                    end
                    

                end
                
                seeds = trainnorm(randidx(1:kk),:);
                
                
                % k-means clustering options
                options = zeros(1,14);
            end
           
            disp('clustering...')
            savefile=[clusterpath '/' clusterfilename];
            %start clustering
            switch clstmeth
                case 'kmeans'
                    if kk~=1
                        [centers,options,trainpost] = ...
                            netlab_kmeans(seeds,trainnorm,options);
                        %savefile=[clusterpath '/' clusterfilename];
                        save(savefile,'centers','clstfeatsel','trainpost');
                        %tz_updaterecord([savefile ' saved'],'results');
                    else
                        trainpost=ones(length(trainclass),1);
                    end
                case 'xmeans'
                    if kk~=1
                        objidcs = tz_xmeans(struct('max_ctrs',kk));
                        trainpost=ml_label2post(objidcs);
                        save(savefile,'centers','clstfeatsel','trainpost');
                    else
                        trainpost=ones(length(trainclass),1);
                    end
                case 'tzkmeans' %remove small clsuters
                    if kk~=1
                        [centers,options,trainpost] = ...
                            tz_kmeans(seeds,trainnorm,options,-1);
                        %savefile=[clusterpath '/' clusterfilename];
                        save(savefile,'centers','clstfeatsel','trainpost');
                        %tz_updaterecord([savefile ' saved'],'results');
                    else
                        trainpost=ones(length(trainclass),1);
                    end
                case 'kmeansmahal'
                    if kk~=1
                        options(5)=1;
                        [centers,options,trainpost] = ...
                            rm_mahalkmeans(seeds,trainnorm,options);
                        %savefile=[clusterpath '/' clusterfilename];
                        save(savefile,'centers','clstfeatsel','trainpost');
%                     tz_updaterecord([savefile ' saved'],'results');
                    else
                        trainpost=ones(length(trainclass),1);
                    end
                case 'classcluster' %classwise clustering
                    if kk==1
                        objidcs=trainclass;
                    else
                        rand('state',0);
                        if kk<-2 %set number of cluters to search
                            ccParam = [1:(-kk) 10];
                        else
                            ccParam = [];
                        end
                        
                        objidcs = ...
                            tz_classcluster(trainnorm, ...
                            trainclass,kk,0,ccParam);   
                    end
                    trainpost=ml_label2post(objidcs);
                    %savefile=[clusterpath '/' clusterfilename];
                    save(savefile,'objidcs','clstfeatsel','trainpost');                
%                 tz_updaterecord([savefile ' saved'],'results');
            end
            disp('done!')
        end
        
        if (strcmp(clstmeth,'kmeans') | strcmp(clstmeth,'xmeans') | ...
                strcmp(clstmeth,'kmeansmahal') | ...
                strcmp(clstmeth,'tzkmeans')) & kk~=1
%             trainset=tz_cellobjcom(combclass(find(trainsel)),combcellidx(find(trainsel)),trainpost);
            [m,objidcs]=max(trainpost');
            objidcs=objidcs';
        end        
        
        %reassign labels
        if ~strcmp(clstmeth,'classcluster') & kk==1
            if ~exist(objtypepath,'dir')
                mkdir(clusterpath,objtypedirname);
            end
            objidcs=ones(length(trainclass),1);
            allpost=ones(nobj,1);
        else
            typelearn=0;
            if ~exist(objtypepath,'dir')
                mkdir(clusterpath,objtypedirname);
            end
            
            if redo<=3 %do first level classification
                typelearn=1;
            else
                if ~exist([objtypepath '/' objtypefilename],'file')
                    typelearn=1;
                    warning('specified object type file not found. New object type assigning will be done');
                else
                    load([objtypepath '/' objtypefilename]);
                    typelearn=0;
                end
                
            end
            if typelearn==1
                if ((strcmp(classifmeth1,'lknn') & t1==1) | ...
                        strcmp(other{2},'leavetrain')) & ...
                        strcmp(other{1},'eval')
                    allpost = [objidcs;tz_classify(testnorm, ...
                            trainnorm,objidcs,classifmeth1,0,t1)];
                else
                    allpost = tz_classify([evalnorm;testnorm], ...
                        trainnorm,objidcs,classifmeth1,0,t1{:});
                end
                save([objtypepath '/' objtypefilename],'allpost');
            end
            %             end
        end
        
        evalpost = allpost(1:length(evalidx));
        testpost = allpost(length(evalidx)+1:end);
        
        %tz- 12-Oct-2006
        % [aic(i), clusterinfo] = rm_akaike_euclid( trainnorm, trainpost);
        %tz--
        
        %tz+ 12-Oct-2006
        aic(i) = tz_aicbic(trainnorm,ml_post2label(trainpost),'euc');
        %tz++
        
        ncluster=kk;
        
        if strcmp(clstmeth,'classcluster')
            ncluster=max(objidcs);
        end
        
        ncluster
        if strcmp(classifmeth2,'nnobj')
            %[training,testing]=ml_featurenorm( ...
            %    combobj(find(trainsel),:),combobj(find(testsel),:));
            training = evalnorm;
            testing  = testnorm;
            traincell=[training,traincellidx];
            if t2==1
                testing=[testing,flindex(testcellidx)];
            end
            rlabels = evalclass;
            [testcell,tlabels]=tz_combfeat2cell(testing, ...
                testclass,testcellidx);
        else
            %Calculate cell features
            if strcmp(featset,'newset')
                trainset=tz_cellobjfeats_clst( ...
                    evalclass, ...
                    evalcellidx, ...
                    evalobj, ...
                    evalpost,ncluster);
                testset=tz_cellobjfeats_clst(testclass, ...
                    testcellidx, ...
                    testobj, ...
                    testpost,ncluster);
            else
                if ~exist(cellfeatpath,'dir')
                    mkdir(objtypepath,cellfeatdirname);
                end

                if redo<=4 %do cell-level feature calculation
                    featcalc=1;
                else
                    if ~exist([cellfeatpath '/' cellfeatfilename],'file')
                        featcalc=1;
                        warning(['Specified feature file not found. '...
                            'Features will be calculated.']);
                    else
                        load([cellfeatpath '/' cellfeatfilename]);
                        featcalc=0;
                    end
                    
                end
                
                if featcalc==1
                    if ~isempty(combcoords)
                        trcombcoords=combcoords(evalidx,:);
                        tscombcoords=combcoords(find(testsel),:);
                    else
                        trcombcoords=[];
                        tscombcoords=[];
                    end
                    
                    if ~isempty(combobjects)
                        trcombobjects=combobjects(evalidx);
                        tscombobjects=combobjects(testidx);
                    else
                        trcombobjects={};
                        tscombobjects={};
                    end
               
                    trainset=tz_cellobjcom(evalclass, ...
                        evalcellidx,...
                        evalpost,...
                        ncluster,featset,evalobj, ...
                        fluofrac(evalidx), ...
                        cofdist(evalidx), ...
                        trcombcoords,trcombobjects);
                    testset=tz_cellobjcom(testclass, ...
                        testcellidx, ...
                        testpost,...
                        ncluster,featset,testobj, ...
                        fluofrac(testidx),cofdist(testidx), ...
                        tscombcoords,tscombobjects);
                    save([cellfeatpath '/' cellfeatfilename], ...
                        'trainset','testset');
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
        
        classif=0;
        classifpath=[cellfeatpath '/' classifdirname];
        
        savefile=['results' num2str(kk) setname num2str(i) 'fold.mat'];
        
        if ~exist(classifpath,'dir')
            mkdir(cellfeatpath,classifdirname);
        end
        
        if redo<=5 %do second level classificaiton
            classif=1;
        else
            if ~exist([classifpath '/' savefile],'file')
                classif=1;
                warning(['specified classification file not found. ' ...
                    'It is doing classification now ...']);
            end
        end
        
        if classif==1
            if ~strcmp(classifmeth2(1:3),'mix')
                tg=tz_classify(testcell,traincell,rlabels, ...
                    classifmeth2,preprocess,t2{:});
                [ncvcm,pcvcm,cvacc]=ml_summaryclassif(tlabels,tg);
                save([classifpath '/' savefile], ...
                    'ncvcm','pcvcm','cvacc','aic');
            else
                if isfield(t2,'fitmeth')
                    fitmeth = t2.fitmeth;
                else
                    fitmeth = t2{5};
                end
                
                if strmatch(fitmeth,{'Gl','Gw'})
                    trobj=[evalobj, ...
                        evalclass,...
                        evalcellidx, ...
                        evalpost];
                    tsobj=[testobj, ...
                        testclass, ...
                        testcellidx, ...
                        testpost];
                    for si=1:nclass
                        trainset{si}=trobj(trobj(:,end-2)==si,:);
                        testset{si}=tsobj(tsobj(:,end-2)==si,:);
                    end
                end

                
                [ress,mress]=tz_decommixtr(trainset,testset, ...
                    struct('ntrial',t2{1},'classnum',t2{2}, ...
                    'cellnums',t2{3},'evalmeths',t2{4},'fitmeth',t2{5}, ...
                    'minalpha',t2{6}));
                save([classifpath '/' savefile],'ress','mress');
            end
        else
            load([classifpath '/' savefile]) 
        end
        if ~strcmp(classifmeth2(1:3),'mix')
            ncvcms{i}=ncvcm;
        else
            avgcm(:,:,i)=1-mress;
        end
            
    end
    
    if ~strcmp(classifmeth2(1:3),'mix')
        [avgcm{ki},avgacc(ki),ncm{ki},kappa(ki),tua,tpa] = ...
            ml_calcavgcm(ncvcms);
        ua(ki,:)=tua;
        pa(ki,:)=tpa;
    else
        avgacc(ki)=mean(avgcm,3);
    end
    ki=ki+1;
    display(setname);
end

function tz_createnewfold(combobj,combcellidx,combclass,nfold,option,savepath,randstate)

%function tz_createnewfold(combobj,combcellidx,combclass,nfold,clstmeth,clustk,classifmeth1,t1,classifmeth2,t2,featset,option,loadpath,permpath,savepath,other,comments)
%
%OVERVIEW:
%   create permutation for cross validation for object level classification

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

%extract fluorescence fraction
flindex=combobj(:,4);

% select features
combobj=combobj(:,[1 3 5:end]);

%%%%%%%%%%%%%%%%%%n fold permutation%%%%%%%%%%%%%%
% if isempty(loadpath)
%     if isempty(permpath)
        for i=1:10
            norgset(i)=max(combcellidx(combclass==i));
            perm{i}=randperm(norgset(i));
        end
        
        for i=1:10
            permfold(i,:)=[0,cumsum(tz_redistrnum(norgset(i),nfold))];
        end
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for kk=clustk
    for i=1:nfold
        foldname=['sel' num2str(i) 'fold' option '.mat']
        newfoldname=['sel' num2str(i) 'fold' '.mat']

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
        
        clusterfilename=[num2str(kk) clustername num2str(i) 'fold.mat'];
        %clustering
        newclustering=1;
        
        if ~isempty(loadpath)
            loadfile=[loadpath '/' num2str(kk) clustername num2str(i) 'fold.mat']
            newclustering=0;
            if(isempty(dir(loadfile)))
                warning('specified clustering file not found. New clustering will be done');
                newclustering=1;
            end
        end
                      
        if newclustering==1    
            %if isempty(permpath)
                %generating training set and testing set
                for j=1:nclass
                    for k=perm{j}(permfold(j,i)+1:permfold(j,i+1))
                        testsel=testsel+((combclass==j)&(combcellidx==k));
                    end
                end
                
                trainsel=ones(nobj,1)-testsel;
                newtrainsel=trainsel;
                newtestsel=testsel;
                load([permpath '/' foldname]);
                
                s=1;
                for j=[1:2 4:10]
                    newtrainsel(allcombclass==j)=trainsel(combclass==s);
                    s=s+1;
                end
                trainsel=newtrainsel;
                testsel=newtestsel;
                
                save([savepath '/' newfoldname],'trainsel','testsel','comments');
                %else
                
                %end
            
            [trainnorm,testnorm]=mb_featurenorm(combobj(find(trainsel),:),combobj(find(testsel),:));
            if strcmp(clstmeth,'kmeans') | strcmp(clstmeth,'kmeansmahal')
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
            switch clstmeth
            case 'kmeans'
                if kk~=1
                    [centers,options,trainpost,errlog]=netlab_kmeans(seeds,trainnorm,options);
                    save([savepath '/' clusterfilename],'centers','trainpost','trainsel','testsel','comments');
                else
                    trainpost=ones(sum(trainsel),1);
                end
            case 'kmeansmahal'
                if kk~=1
                    options(5)=1;
                    [centers,options,trainpost,errlog]=rm_mahalkmeans(seeds,trainnorm,options);
                    save([savepath '/' clusterfilename],'centers','trainpost','trainsel','testsel','comments');
                else
                    trainpost=ones(sum(trainsel),1);
                end
            case 'classcluster'
                if kk==1
                    objidcs=combclass(find(trainsel));
                else
                    rand('state',0);
                    objidcs=tz_classcluster(trainnorm,combclass(find(trainsel)&combclass==3),kk,0);   
                end
                trainpost=tz_label2post(objidcs);
                save([savepath '/' clusterfilename],'objidcs','trainpost','trainsel','testsel','comments');                
            end
            disp('done!')
        else
            load(loadfile)
            [trainnorm,testnorm]=mb_featurenorm(combobj(find(trainsel),:),combobj(find(testsel),:));
        end
        
        if (strcmp(clstmeth,'kmeans') | strcmp(clstmeth,'kmeansmahal')) & kk~=1
%             trainset=tz_cellobjcom(combclass(find(trainsel)),combcellidx(find(trainsel)),trainpost);
            [m,objidcs]=max(trainpost');
            objidcs=objidcs';
        end        
        
        %reassign labels
        if ~strcmp(clstmeth,'classcluster') & kk==1
            objidcs=ones(sum(trainsel),1);
            allpost=ones(nobj,1);
        else
            if (strcmp(classifmeth1,'lknn') & t1==1) | (strcmp(clstmeth,'kmeans') & strcmp(classifmeth1,'dist')) | strcmp(other,'leavetrain')
                allpost=[objidcs;tz_classify(testnorm,trainnorm,objidcs,classifmeth1,0,t1)];
            else
                allpost=tz_classify([trainnorm;testnorm],trainnorm,objidcs,classifmeth1,0,t1);
            end
        end
        
        [aic(i), clusterinfo] = rm_akaike_euclid( trainnorm, trainpost);
        
        ncluster=kk;
        
        if strcmp(clstmeth,'classcluster')
            ncluster=max(objidcs);
        end
        
        ncluster
        
        %Calculate cell features
        if strcmp(featset,'newset')
            trainset=tz_cellobjfeats_clst(combclass(find(trainsel)),combcellidx(find(trainsel)),combobj(find(trainsel),:),allpost(1:length(objidcs)),ncluster);
            testset=tz_cellobjfeats_clst(combclass(find(testsel)),combcellidx(find(testsel)),combobj(find(testsel),:),allpost((length(objidcs)+1):end),ncluster);
        else
            trainset=tz_cellobjcom(combclass(find(trainsel)),combcellidx(find(trainsel)),allpost(1:length(objidcs)),ncluster,featset,combobj(find(trainsel),:),flindex(find(trainsel)));
            testset=tz_cellobjcom(combclass(find(testsel)),combcellidx(find(testsel)),allpost((length(objidcs)+1):end),ncluster,featset,combobj(find(testsel),:),flindex(find(testsel)));
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
        tg=tz_classify(testcell,traincell,rlabels,classifmeth2,preprocess,t2);
        [ncvcm{i},pcvcm{i},cvacc(i)]=tz_summaryclassif(tlabels,tg);
    end
    
    [avgcm,avgacc]=tz_calcavgcm(ncvcm);
    save([savepath '/' 'results' num2str(kk) setname num2str(i) 'fold.mat'],'ncvcm','pcvcm','cvacc','aic','comments');
   
    display(setname);
    
%     rand('state',0);
end

rand('state',oldseed);
function gnf_objtrain(combobj,combcellidx,combclass,clstmeth,clustk,classifmeth,t,featset,workdir)

if ~exist(workdir,'dir')
    mkdir(workdir);
end

datasize=size(combobj);
nclass=max(combclass);
nobj=size(combobj,1);
% number of cells in each class
ncells = zeros(1,nclass);
for i=1:nclass
    ncells(i) = max(combcellidx(combclass==i));
end

% ------------------------- %
% SDA code %
clstfeatsel = [];
% ------------------------- %

% directory for saving clustering results
clustername = clstmeth;%[clstmeth num2str(clstsda) other{1}];

% directory for saving object type learning results
objtypedirname = classifmeth;%[classifmeth1 tz_cell2str(t1,[]) other{2}]

[trainnorm,testnorm] = ml_featurenorm(combobj,combobj);

ki = 1;
for kk=clustk
    
    %directory for saving clustering results
    clusterdirname=[num2str(kk) clustername];
    clusterpath=[workdir '/' clusterdirname];
    objtypepath=[clusterpath '/' objtypedirname];
    if ~exist(clusterpath,'dir')
            mkdir(workdir,clusterdirname);
    end
    clusterfilename=[num2str(kk) clustername '.mat'];
    
    if strcmp(clstmeth,'kmeans') | ...
            strcmp(clstmeth,'kmeansmahal') | ...
            strcmp(clstmeth,'tzkmeans')
    
        rand('state',0);
        randidx = randperm(size(trainnorm,1));
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
                %tz_updaterecord([savefile ' saved'],'results');
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
            %tz_updaterecord([savefile ' saved'],'results');
    end
    disp('done!')
end
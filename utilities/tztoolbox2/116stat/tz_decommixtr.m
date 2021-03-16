function [ress,mress]=tz_decommixtr(trainset,testset,param)
%TZ_DECOMMIXTR Mixture decomposition.
%   RESS = TZ_DECOMMIXTR(TRAINSET,TESTSET,PARAM) returns the average error 
%   of mixture decomposition. TRAINSET is a cell array of training feature
%   matrices. Each matrix corresponds to one component. TESTSET is also a
%   cell array of feature matrices. During testing, they will be randomly
%   mixed. PARAM is a structure of parameters. It has the following fields:
%       'ntrial' - number of trials. default 10.
%       'classnum' - number of classes to mix. If it is 0, the number will
%           generated randomly. default -1.
%       'cellnums' - number of cells from each class to mix. It is a vector,
%           and all numbers in the vector will be tried. default 1.
%       'evalmeths' - significance evaluation methods. It is an integer vector
%           corresponding one or more of the following strings:
%               'sig' : signifcance level calcuated by TZ_EVALMIXSIG 
%                       with the method 'norm'
%               'sig2': signifcance level calcuated by TZ_EVALMIXSIG 
%                       with the method 'nonp'
%               'abs' : total absolute error
%               'sq2' : total squre error
%               'per' : error rate
%               'mpe' : error rate after mering class 3 and 4, 5 and 9. (for
%                       merging giantin and gpp, lysosome and endosome)
%               'mp2' : error rate after merging class 3 and 4. (for merging
%                       giantin and gpp)
%       'fitmeth' - decomposition method. 
%               'Hf' : multinomial on fixed number of components
%               'Hs' : multinomial decomposition
%               'Hm' : multinomial decomposition with adjusted data
%               'Hn' : multinomial decomposition with object size weights
%               'Ht' : multinomial decomposition  with success
%               'Ln' : linear regression
%               'Lo' : linear regression with adjusted data
%               'Lf' : linear regression with fixed number of components       
%               'Gl' : maximum likelihood estimation for general distribution
%               'Gw' : weighted maximum likelihood estimation for general 
%                      distribution
%               'nn' : neural network estimation
%
%   RESS = TZ_DECOMMIXTR(TRAINSET,TESTSET) uses the default parameters.
%
%   [RESS,MRESS] = TZ_DECOMMIXTR(...)
%   
%   See also TZ_DECOMMIX_LK

%   19-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param, ...
    struct('ntrial',10,'classnum',0,'cellnums',1, ...
    'evalmeths',5,'fitmeth','Hs','minalpha',[]));

inc=1;

%backup original data
tmptrainset=trainset;
tmptestset=testset;

%significance evaluation methods
levalmeths={'sig','sig2','abs','sq2','per','mpe','mp2'};
evalmeths = param.evalmeths;
evalmeths=levalmeths(evalmeths);

fitmeth = param.fitmeth;
ntrial = param.ntrial;
minalpha = param.minalpha;
cellnums = param.cellnums;
classnum = param.classnum;

for cellnum=cellnums 
    for ne=1:length(evalmeths)
        evalmeth=evalmeths{ne};
        trainset=tmptrainset;
        testset=tmptestset;
        switch evalmeth
        case 'mp2'  %merge class 3 4
            trainset{3}=[trainset{3};trainset{4}];
            trainset(4)=[];
            testset{3}=[testset{3};testset{4}];
            testset(4)=[];          
        case 'mpe'  %merge class 3,4 and 5 9
            trainset{3}=[trainset{3};trainset{4}];
            trainset{5}=[trainset{5};trainset{9}];
            trainset(4)=[];
            trainset(8)=[];
            testset{3}=[testset{3};testset{4}];
            testset{5}=[testset{5};testset{9}];
            testset(4)=[];
            testset(8)=[];
        end
        
        if ~strcmp(fitmeth,'Gl') & ~strcmp(fitmeth,'nn') ...
                & ~strcmp(fitmeth,'Gw')
            [objsq,objnum]=tz_trainobjcom(trainset,[]);
            p=tz_normobjcom(objsq);
            pw=sum(objsq,2);
        end
        
        if strcmp(fitmeth,'nn')
            training=[];
            group=[];
            for i=1:trainset
                training=[training;trainset{i}];
                group=[group;i+zeros(size(trainset{i},1),1)];
            end
            regmodel = tz_bpnnunmix_clst(training,group);
        end
        
        nclass=length(trainset);    
        res=0;
        for i=1:ntrial
            scell=[];
            scell2=[];

            %decide number of basic components
            if classnum<=0
                snclass=randperm(nclass);
            else
                snclass=classnum;
            end
            
            %decide classes of basic components
            sclass=randperm(nclass);       
            sc=zeros(1,nclass);
            sc(sclass(1:snclass(1)))=1;
            
            %randomly select cells to form testing set
            for j=1:nclass
                if sc(j)==1
                    if strcmp(fitmeth,'Gl') | strcmp(fitmeth,'Gw')
                        % for general likelihood
                        testcell=ml_findclass(testset{j}(:,1:end-1));         
                        perm=randperm(length(testcell));
                        tmpcell=testcell{perm(1)};
                    else
                        perm=randperm(size(testset{j},1));
                        tmpcell=testset{j}(perm(1),:);
                    end

                    
                    for k=2:cellnum
                        if strcmp(fitmeth,'Gl') | strcmp(fitmeth,'Gw')
                            tmpcell=[tmpcell;testcell{perm(k)}];
                        else
                            tmpcell=tmpcell+testset{j}(perm(k),:);
                        end
                    end
                    scell=[scell;tmpcell];
                end
            end
            
            if ~strcmp(fitmeth,'Gl') | strcmp(fitmeth,'Gw')
                y=sum(scell,1);
            end
            
            talpha=zeros(1,nclass);
                        
            
            
            if strcmp(fitmeth,'Gl') | strcmp(fitmeth,'Gw')
                for i=1:nclass
                    talpha(i)=sum(scell(:,end-2)==i);
                end
                talpha=talpha/sum(talpha);
            else
                talpha1=talpha;
                talpha1(find(sc))=cellnum;
                talpha1=talpha1/sum(talpha1);
                
                talpha2=talpha;
                talpha2(find(sc))=sum(scell,2)';
                talpha2=talpha2/sum(talpha2);
            end
            
            switch fitmeth
%             case 'Hx'
%                 ealpha=tz_fitmnomk(y2,p2,minalpha);
%                 ealpha2=tz_fitmnomk(y,p(ealpha>0,:),0);
%                 ealpha(ealpha>0)=ealpha2;
%                 talpha=talpha2;
            case 'Hf'
                succ=0;
                while succ==0
                    [ealpha,loglk,succ]=tz_fitmnomk_fix(y,p,snclass(1),[],minalpha);
                end
                talpha=talpha2;
            case 'Hs'
                succ=0;
                %%%%%%%%%
%                 y=y(:,[1:15, 17:29,31:end]);
                %%%%%%%%%%%
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
                
%             case 'Rd'
%                 [fitx,logll]=ml_redistrmnom(y,p);
%                 ealpha=sum(fitx,2);
%                 ealpha=ealpha/sum(ealpha);
%                 ealpha=ealpha';
%                 talpha=talpha2;
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
            case 'Gl'
                for i=1:nclass
                    %training{i}=trainset{i}(:,[1:end-3 end]);
                    training{i}=trainset{i}(:,[1:end-3]);
                    weights{i} = trainset{i}(:,1);
%                     objtypes{i}=trainset{i}(:,end);
                end
%                 ealpha=tz_decommix_lk(scell(:,1:end-3),training, ...
%                     'mvno',[],scell(:,1),{weights});
                ealpha = ml_decommixlk(scell(:,1:end-3),training);
                %scell(:,1),{});
            case 'Gw'
                for i=1:nclass
                    training{i}=trainset{i}(:,[1:end-3]);
                    weights{i} = trainset{i}(:,1);
%                     objtypes{i}=trainset{i}(:,end);
                end
%                 ealpha=tz_decommix_lk(scell(:,1:end-3),training, ...
%                     'mvnw',[],scell(:,1),{weights});   
%                 ealpha=tz_decommix_lk(scell(:,1:end-3),training, ...
%                     'mvnw',[],scell(:,1),[]);  
                ealpha = ml_decommixlk(scell(:,1:end-3),training, ...
                    struct('weights',scell(:,1)));
            case 'nn'
                ealpha=tz_evalbpnnreg(y,regmodel);
                ealpha=tz_normobjcom(ealpha,0);
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
                err=1-sum(min([ealpha;talpha]));
            case 'mp2'
                err=1-sum(min([ealpha;talpha]));
            end
            
            res=res+err;
        end
        
        res=res/ntrial;
        ress{inc}=struct('res',res,'cellnum',cellnum,'evalmeth',evalmeth);
        inc=inc+1;
    end
    
end

inc=1;

for j=1:length(cellnums)
    for k=1:length(evalmeths)
        mress(j,k)=ress{inc}.res;
        inc=inc+1;
    end
end

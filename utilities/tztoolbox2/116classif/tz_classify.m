function class = tz_classify(sample,training,group,classmeth, ...
    preprocess,varargin)
%TZ_CLASSIFY Supervised learning. It will be replaced in the future.
%   CLASS = TZ_CLASSIFY(SAMPLE,TRAINING,GROUP,CLASSMETH,PREPROCESS) returns the
%   predicted classes of the testing [feature matrix] SAMPLE from the
%   classifier trained on the training [feature matrix] TRAINING and its [label
%   vector] GROUP. CLASSMETH is set to specify what classifier will be trained:
%       'lda' - Linear discriminant analysis (LDA). See CLASSIFY.
%       'svm' - Support vector machine (SVM). See TZ_SVMCLASSIFY.
%       'mnlogl' - Clasissification based on multinomial models.
%           See TZ_MNLOGLCLASSIFY.
%       'logl' - maximum log likehood. See TZ_LOGLCLASSIFY.
%       'dist' - Nearest center classifier (NCC). See TZ_DISTCLASSIFY.
%       'knn' - K-nearest neighbor (KNN) classifier. See TZ_KNNCLASSIFY.
%       'lknn' - KNN for large samples. See TZ_LKNNCLASSIFY.
%       'bpnn' - Back propagation neural network (BPNN). See TZ_BPNNCLASSIFY.
%       'nnobj' - Classification based on object neighbors. 
%           See TZ_NNOBJCLASSIFY.
%       'kde' - Kernel density estimation. See TZ_KDECLASSIFY.
%   PREPROCESS is a number to specify how to preprocess the data:
%       0 - no preprocessing
%       1 - z-score
%       2 - sda and z-score
%       3 - sda for object level classification
%       4 - sda
%
%   CLASS = TZ_CLASSIFY(SAMPLE,TRAINING,GROUP,CLASSMETH,PREPROCESS,T) 
%   customizes parameters for the selected classifier. See the correspoding
%   classification function for more details.
%   
%   CLASS = TZ_CLASSIFY(SAMPLE,TRAINING,GROUP,CLASSMETH,PREPROCESS,T,T2) is
%   only useful for 'nnobj'. T2 is the number of clusters.
%   
%   See also ML_TRAINCLASSIF

%   ??-???-???? Initial write TINGZ
%   17-MAY-2004 Add comments TINGZ
%   09-FEB-2005 Add comments TINGZ
%   09-FEB-2005 Modified TINGZ
%       - change error message

if nargin<7
    varargin{2}=1;
end

if nargin<6
    varargin{1}=1;
end

t=varargin{1};
t2=varargin{2};

if ~strcmp(classmeth,'mnlogl') & ~strcmp(classmeth,'nnobj')
    switch preprocess
    
        case 1 %normalization
            [training,sample]=ml_featurenorm(training,sample);
        case 2 %sda and normalization
            ntotal=size(training,2);
            fullidx=1:ntotal;
            
            constidx=tz_constidx(training);
            if ~isempty(constidx)
                training(:,constidx)=[];
                sample(:,constidx)=[];
                fullidx(constidx)=[];
                disp(['constant features removed:' num2str(constidx)]);
            end
            
            featsel=tz_sda(training,group);
            
            training=training(:,featsel);
            sample=sample(:,featsel);
            featsel=fullidx(featsel);
            save('/home/tingz/tmp/featsel.mat','featsel');
            [training,sample]=ml_featurenorm(training,sample);
        case 4 %sda only
            ntotal=size(training,2);
            fullidx=1:ntotal;
            
            constidx=ml_constidx(training);
            if ~isempty(constidx)
                training(:,constidx)=[];
                sample(:,constidx)=[];
                fullidx(constidx)=[];
                disp(['constant features removed:' num2str(constidx)]);
            end
            
            featsel=tz_sda(training,group);
            
            training=training(:,featsel);
            sample=sample(:,featsel);
            featsel=fullidx(featsel);
            save('/home/tingz/tmp/featsel.mat','featsel');
        case 3 %sda for repeated features, especially for object level classification
            ntotal=size(training,2);
            fullidx=1:ntotal;
            subidx=1:ntotal;
            
            featmeannum=ntotal/t2-2; %number of mean features. 
                                     %t2 is the number of clusters
          
            subtrainset=training(:,1:t2); %partial features
            constidx=ml_constidx(subtrainset); %constant features
            
            if ~isempty(constidx) %remove const features
                subtrainset(:,constidx)=[];
                subidx(constidx)=[];          
                disp(['constant features removed:' num2str(constidx)]);
            end
            
            featsel=ml_sda(subtrainset,group); %feature selection
            subidx=[subidx(featsel) fullidx((t2+1):end)];
            
            rmidx=fullidx;
            rmidx(subidx)=[]; %indices of features to be removed
            
            rmidx2=t2+rmidx; %all indices of features to be removed
            for i=1:length(rmidx)
                rmidx2=[rmidx2 t2*2+(rmidx(i)-1)*featmeannum+(1:featmeannum)];
            end
            
            rmidx=[rmidx,rmidx2];
            subidx=fullidx;
            subidx(rmidx)=[];
        
            subtrainset=training(:,subidx);
            
            constidx=ml_constidx(subtrainset);
            
            if ~isempty(constidx)
                subtrainset(:,constidx)=[];
                subidx(:,constidx)=[];          
                disp(['constant features removed:' num2str(constidx)]);
            end
            
            featsel=ml_sda(subtrainset,group);
            featsel=subidx(featsel);
            
            save('/home/tingz/tmp/featsel.mat','featsel');
            training=training(:,featsel);
            sample=sample(:,featsel);
            [training,sample]=ml_featurenorm(training,sample);
    end
end

switch classmeth
    case 'lda'
        if ~ischar(t)
            class = classify(sample,training,group);
        else
            class = classify2(sample,training,group,t);
        end
    case 'svm'
        class = tz_svmclassify(sample,training,group,t);
    case 'mnlogl'
        class = tz_mnloglclassify(sample,training,group,t);
    case 'logl'
        class = tz_loglclassify(sample,training,group,t);
    case 'dist'
        class = tz_distclassify(sample,training,group);
    case 'knn'
        class = tz_knnclassify(sample,training,group,t);
    case 'lknn'
        class = tz_lknnclassify(sample,training,group,t);
    case 'bpnn'
        class = tz_bpnnclassify(sample,training,group,t);
    case 'nnobj'
        class = tz_nnobjclassify(sample,training,group,t);
    case 'kde'
        tmp=version;
        if str2num(tmp(1:3))<6.5
            error('the option kde can not work for version lower than 6.5')
        end    
        class = tz_kdeclassify(sample,training,group,t);    
    case 'libsvm'
        if t==0 %search parameters automatically
            classifier = ml_trainclassif(training,group,struct('args', ...
                                                              []),4);
        else
            classifier = ml_trainclassif(training,group,struct([]),4);
        end
        class = ml_evalreg(sample,classifier);
    otherwise
        disp('The classification method should be one of these:')
        disp('lda svm mnlogl dist knn lknn bpnn nnobj kde');
        error('invalid classification method');
end

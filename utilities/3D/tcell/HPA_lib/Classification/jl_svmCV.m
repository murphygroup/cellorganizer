function [uwacc,wacc,macc] = jl_svmCV(c,g,filetag,tuning,dataname,predicting,Nfolds,Ninit,Nclass,SDA)

%addpath(genpath('./lib_Justin/lib/classify'));
%addpath('./lib_Justin/lib/results');

%addpath('./lib_Justin/HPAconfocal_lib/classification')
%addpath('./lib_Justin/HPAconfocal_lib/classification/libsvm-mat-2.84-1')

%addpath('./lib_Justin/lib/classification/libsvm-mat-2.84-1')


if nargin==1 %%
   c
   dataname = c;
   clear c
end

if ~exist('tuning','var')
   tuning = 0;
end
tuning

if ~exist('predicting','var')
   predicting = 0;
end
predicting

if ~exist('dataname','var')
   %dataname = 'regionfeatures_single_A431.mat';
   dataname = 'regionfeatures_all_A431.mat'; %%default
end
dataname

if ~exist('Ninit','var')
   Ninit = 5;
end
Ninit

if ~exist('Nclass','var')
   Nclass = Ninit;
elseif Nclass<Ninit
   warning('Nclass cannot be smaller than Ninit and has been increased to equal to Ninit!');
   Nclass = Ninit;
end
Nclass

if ~exist('SDA','var')
   SDA = 1;
end
SDA

clear tmp_Nfolds
if exist('Nfolds','var')
   tmp_Nfolds = Nfolds;
end


tmpidx1 = strfind(dataname,'.');
tmpidx2 = strfind(dataname,'/');
tmpname=dataname;
if isempty(tmpidx2)
   tmpidx2 = 0;
   if ~isempty(tmpidx1)
      tmpname = tmpname(1:tmpidx1(end)-1);
   end
else
   tmpname(tmpidx2(end))='_';
   tmpname = tmpname(tmpidx2(end-1)+1:tmpidx1(end)-1);
end
tmpname

if exist(['./data/classification/',tmpname,'_',num2str(Ninit),'folds.mat'])
   ['./data/classification/',tmpname,'_',num2str(Ninit),'folds.mat']
   load(['./data/classification/',tmpname,'_',num2str(Ninit),'folds.mat']);
else

tmp_g = g;
load(dataname)

if exist('f','var')
%   f = f(1:10:end,:); %
%   classlabels = classlabels(1:10:end,:); %
%   antibodyids = antibodyids(1:10:end,:); %
   allfeatures = f;
   clear f g
   g = tmp_g;
   SDA = 0;
   allfeatures(isnan(allfeatures)) = 0;
end

classlabels = classlabels(:);
imagelist = imagelist(:);
classes = classes(:);

[rows,cols] = find(isnan(allfeatures));
allfeatures(rows,:) = [];
classlabels(rows,:) = [];
antibodyids(rows,:) = [];
imagelist(rows,:) = [];
specificity(rows,:) = [];
staining(rows,:) = [];

tmpmax = max(allfeatures,[],1);
tmpmin = min(allfeatures,[],1);
tmpidx = find((tmpmax-tmpmin)<=1e-2);
allfeatures(:,tmpidx) = [];

%iidx = find(isinf(allfeatures));  %%
%allfeatures(iidx) = sign(allfeatures(iidx)).*realmax('single');
[iidx,iidy] = find(isinf(allfeatures));
allfeatures(:,iidy) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%delete classes with small number of antibodies.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%class_set=unique(classlabels);
%for i=1:length(class_set)
%c(i)=length(unique(antibodyids(find(classlabels==class_set(i)))));
%end

class_del = [];
class_set = unique(classlabels);
classes_idx = 1:length(classes);
if length(unique(classlabels))<length(classes)
   warning('Some classes are not in the classlabels. Delete them... Check data!');
   tmpind = ~ismember(classes_idx,class_set); %%
   class_del = [class_del,classes_idx(tmpind)];
elseif length(unique(classlabels))>length(classes)
   warning('Some classlabels are not in the classes. Delete them... Problematic and Data-Check recommended!');
   tmpind = ~ismember(class_set,classes_idx); %%
   class_del = [class_del,class_set(tmpind)'];
end

for i=1:length(class_set)
    numc(i)=length(unique(antibodyids(find(classlabels==class_set(i)))));
    if numc(i)<Nclass
       class_del = [class_del,class_set(i)];
    end
end
tmpidx = ismember(classlabels,class_del);
allfeatures(tmpidx,:) = [];
antibodyids(tmpidx) = [];
classlabels(tmpidx) = [];
classes(class_del) = [];
class_del = sort(class_del);
for i=1:length(class_del)
    tmpdex = find(classlabels>class_del(i));
    classlabels(tmpdex) = classlabels(tmpdex)-1;
    class_del = class_del-1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%keyboard


boolvar = 0;
if exist('Nfolds','var')
    boolvar = 1;
end

abids2 = antibodyids;

%Ninit = 5;
[trainidx,testidx, protlabels,N] = splitbyproteins( classlabels, abids2, Ninit);

end %%
N


if ~exist('Nfolds','var')
    Nfolds = 1:N;
end
Nfolds

if exist('tmp_Nfolds','var')
   disp(['Only run the fold ',num2str(tmp_Nfolds),'.']);
   Nfolds = tmp_Nfolds;
end
Nfolds


U = unique(classlabels);

if ~exist('c','var')
   c = 2.^(2:2:8);
end
if ~exist('g','var')
   g = 2.^(-2:2:4);
end
if length(c)==1 && length(g)==1
   if ~exist('filetag')
      warning('We need a filetag to distinguish result from each different set of c and g or no parameter tuning!!!'); %%
   end
else
   filetag = []; %%
end 

predlabels = zeros(length(classlabels),1);
allweights = zeros(length(classlabels),length(U));
committee = [];
%allfeatures = double(allfeatures);
allfeatures = single(allfeatures);
%keyboard
c_select = zeros(length(Nfolds),1);
g_select = zeros(length(Nfolds),1);
acc_tuning = zeros(length(Nfolds),1);
wopt = [];
for i=Nfolds
    disp(['evaluating fold: ' num2str(i)]); 
    if exist('protlabels','var')
       trainind = find(ismember(protlabels,trainidx{i}));
       testind = find(ismember(protlabels,testidx{i}));
    elseif exist('protlabelsS','var') && exist('protlabelsM','var')
       trainind = ind_single(ismember(protlabelsS,trainidxS{i}));
       testind = [ind_single(ismember(protlabelsS,testidxS{i}));ind_multi(ismember(protlabelsM,testidxM{i}))];
    else
       error('No protlabels from folds!');
    end


    traindata = allfeatures( trainind,:);
    testdata = allfeatures( testind,:);
    
    trainlabels = classlabels(trainind);
    testlabels = classlabels(testind);

    [traindata, testdata, mu, st] = standardizeFeatures( traindata, testdata);

    % normalizing
    traindata = normalizeFeatures( traindata);
    testdata = normalizeFeatures( testdata);

    trainprots = antibodyids(trainind);
    
    trainabids2 = trainprots;

    if SDA == 1
       idx = featureSelection( traindata, trainlabels, 'sda', [tmpname,num2str(filetag),num2str(i)]);
    else
       idx = 1:size(allfeatures,2);
    end
    idx = unique(idx);
    length(idx)
    traindata = traindata(:,idx);
    testdata = testdata(:,idx);

    %keyboard
    
    %if length(c)>1 || length(g)>1 %%
    if tuning==1
    
    if ~exist(['./data/classification/',tmpname,'_',num2str(Ninit),'folds.mat']) %%
       [traintrainidx, traintestidx, trainprotlabels, trainN] = splitbyproteins( trainlabels, trainabids2, Ninit-1); %%
    end

    %all_acc = zeros(length(traintrainidx),length(c),length(g));
    all_acc = zeros(length(traintrainidx{i}),length(c),length(g));  %%
    for zed = 1:length(traintrainidx{i})
        traintrainind = find(ismember(trainprotlabels{i},traintrainidx{i}{zed}));
        traintestind = find(ismember(trainprotlabels{i},traintestidx{i}{zed}));

        traintraindata = traindata( traintrainind,:);
        traintestdata = traindata( traintestind,:);
    
        traintrainlabels = trainlabels(traintrainind);
        traintestlabels = trainlabels(traintestind);

        ma = 0;
        wopt = classweights( traintrainlabels);  %uw
        for j=1:length(c)
            for k=1:length(g)
                disp(['Tuning for parameters: C: ',num2str(c(j)),' g: ',num2str(k),'.']);
                options = ['-s 0 -t 2 -c ' num2str(c(j)) ' -g ' num2str(g(k)) ' -b 0 ' wopt];
                
                %keyboard
                model = svmtrain( traintrainlabels, double(traintraindata), options); %traintraindata(:,idx)
                [predict_label accuracy weights] = svmpredict( traintestlabels, double(traintestdata), model, '-b 0'); %traintestdata(:,idx)
                %[ccc uuu ddd www] = conmatrix(traintestlabels, predict_label);
                all_acc(zed,j,k) = sum(predict_label(:)==traintestlabels(:));
            end
        end
        clear traintraindata traintestdata
    end
    %s(:,:) = mean(all_acc,1);
    s(:,:) = sum(all_acc,1)/length(trainlabels);
    [row col] = find(s==max(s(:)));
    acc_tuning(i) = max(s(:))
    
    if length(row)>1
        row = row(floor(length(row)/2));
        col = col(floor(length(col)/2));
    end
    c(row)
    c_select(i)=c(row);
    g(col)
    g_select(i)=g(col);

    end

    if predicting==1

    if (length(c)>1 || length(g)>1) && (tuning==0)
       error('When tuning==0 it means that we are calculating the final accuracy; c and g should both be single chosen values.');
    end
    row = 1;
    col = 1;
    %end
    c(row)
    c_select(i)=c(row);
    g(col)
    g_select(i)=g(col);

    
    wopt = classweights( trainlabels);  %uw
    options = ['-s 0 -t 2 -c ' num2str(c(row)) ' -g ' num2str(g(col)) ' -b 1 ' wopt];
    %options = ['-s 0 -t 0 -c ' num2str(c(row)) ' -b 1 ' wopt];

    model = svmtrain( trainlabels, double(traindata), options); %traindata(:,idx)
    [predict_label accuracy weights] = svmpredict( testlabels, double(testdata), model, '-b 1'); %testdata(:,idx)

    [o1 o2 o3 o4] = conmatrix( testlabels, predict_label)

    committee(i).model = model;
    committee(i).idx_feat = idx;
    committee(i).options = options;
    committee(i).mu = mu;
    committee(i).st = st;
    
    [a ind] = max(weights,[],2);
    predlabels(testind) = model.Label(ind);
    allweights(testind,model.Label) = weights;
    end
    clear traindata testdata model
end
acc_tuning_ave = mean(acc_tuning);

[cc1 uu1 dd1 ww1] = conmatrix( classlabels, predlabels)

Nstr = [];
for i=1:length(Nfolds)
    Nstr = [Nstr num2str(Nfolds(i)) '_'];
end
%Nstr(end) = [];

[protlabels,protpredlabels,proteins,newallweights] = ...
    modeAcrossProtein( allweights, classlabels, antibodyids);

[cc2 uu2 dd2 ww2] = conmatrix( protlabels, protpredlabels);
cc2
dd2



[protlabels,protpredlabels,proteins,newallweights] = ...
    voteAcrossProtein( allweights, classlabels, antibodyids);

[cc3 uu3 dd3 ww3] = conmatrix( protlabels, protpredlabels);
cc3
dd3



[tmp_protlabels,tmp_protpredlabels,proteins] = ...
    jl_voteAcrossProtein( classlabels, antibodyids, predlabels);

[tmp_cc3 tmp_uu3 tmp_dd3 tmp_ww3] = conmatrix( tmp_protlabels, tmp_protpredlabels)



if ~isstr(filetag)
   filetag = num2str(filetag);
end
if exist('tmp_Nfolds','var')
   filetag = [filetag,'-',num2str(tmp_Nfolds)];
end

%keyboard
if nargout == 0
   if boolvar
      save(['./data/classification/' 'svm' Nstr filetag tmpname], 'allweights','committee','filetag','Nstr','antibodyids','acc_tuning','classlabels',...
      'cc1','uu1','dd1','ww1','cc2','uu2','dd2','ww2','cc3','uu3','dd3','ww3','classes','Ninit','Nclass','Nfolds','predlabels','c','g','c_select','g_select','class_del','numc',...
      'protlabels','protpredlabels','acc_tuning_ave','SDA','tmp_cc3','tmp_uu3','tmp_dd3','tmp_ww3','tmp_protlabels','tmp_protpredlabels');
   else
      save(['./data/classification/' 'svm' filetag tmpname], 'allweights','committee','filetag','Nstr','antibodyids','acc_tuning','classlabels',...
      'cc1','uu1','dd1','ww1','cc2','uu2','dd2','ww2','cc3','uu3','dd3','ww3','classes','Ninit','Nclass','Nfolds','predlabels','c','g','c_select','g_select','class_del','numc',...
      'protlabels','protpredlabels','acc_tuning_ave','SDA','tmp_cc3','tmp_uu3','tmp_dd3','tmp_ww3','tmp_protlabels','tmp_protpredlabels');  %svmuw
   end
else
   uwacc = uu1
   wacc = uu3
   macc = uu2
end





function [protlabels,protpredlabels,proteins] = jl_voteAcrossProtein(classlabels,antibodyids,predlabels)

proteins = unique(antibodyids);
protpredlabels = zeros(size(proteins));
protlabels = zeros(size(proteins));
for i=1:length(proteins)
    idx = find(proteins(i)==antibodyids);

    protpredlabels(i) = mode(predlabels(idx));
    protlabels(i) = mode(classlabels(idx));
end


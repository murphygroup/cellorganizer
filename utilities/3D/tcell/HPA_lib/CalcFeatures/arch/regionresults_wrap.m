clear

addpath lib/classification
addpath lib/classification/libsvm-mat-2.84-1
addpath lib/results

load regionlrclass.mat

[a predlabels] = max(allweights,[],2);

[C U D W] = conmatrix( classlabels, predlabels);

Ntot = sum(C(:));
Ncorr = U*Ntot;

D = round(D*1000)/10
U = round(U*1000)/10
W = round(W*1000)/10

% classes(idx_rmclass) = [];
% classes'

fid = fopen('results_region.csv','w');

fwrite(fid,['Single region confusion matrix,' char(10)]);
fwrite(fid,',');
for i=1:length(classes)
    fwrite(fid,classes{i});
    fwrite(fid,',');
end
fwrite(fid,['#samples,' char(10)]);

for i=1:length(classes)
    fwrite(fid,classes{i});
    fwrite(fid,',');
    for j=1:length(classes)
        fwrite(fid,num2str(D(i,j)));
        fwrite(fid,',');
    end
    fwrite(fid, [num2str(sum(C(i,:))) ',' char(10)]);
end
fwrite(fid,['Accuracy,' num2str(U) ',' char(10)]);
fwrite(fid,['Class weighted accuracy,' num2str(W) ',' char(10)]);



[val predlabels] = max(allweights, [], 2);
[s ind] = sort(val);

re = zeros(length(s),1);
for i=length(s):-1:1
    prl = predlabels(ind(i:end));
    cll = classlabels(ind(i:end))';
    pr(i) = sum(prl==cll)/length(prl);
    re(i) = pr(i)*length(prl)/Ntot;
end

% plot(re,pr,'.');
% xlabel('Recall');
% ylabel('Precision');

rthr = 0.6;
idx = find(re>rthr);
i = idx(end);
pthr = s(i);

prl = predlabels(ind(i:end));
cll = classlabels(ind(i:end));

fwrite(fid,[char(10) char(10) char(10)]);

[Cpr Upr Dpr Wpr] = conmatrix( cll, prl);

Dpr = round(Dpr*1000)/10
Upr = round(Upr*1000)/10
Wpr = round(Wpr*1000)/10

fwrite(fid,['Single region with recall = ' num2str(rthr*100) '% / prob. thr = ' num2str(pthr) ',' char(10)]);
fwrite(fid,',');
for i=1:length(classes)
    fwrite(fid,classes{i});
    fwrite(fid,',');
end
fwrite(fid,['#samples,' char(10)]);

for i=1:length(classes)
    fwrite(fid,classes{i});
    fwrite(fid,',');
    for j=1:length(classes)
        fwrite(fid,num2str(Dpr(i,j)));
        fwrite(fid,',');
    end
    fwrite(fid, [num2str(sum(Cpr(i,:))) ',' char(10)]);
end
fwrite(fid,['Accuracy,' num2str(Upr) ',' char(10)]);
fwrite(fid,['Class weighted accuracy,' num2str(Wpr) ',' char(10)]);





% classif by images
[Ui ii imagelabels] = unique(imagelist);
Ui = unique(imagelabels);
zed = zeros(length(imagelist),1);
ppredlabels = zeros(length(imagelist),1);
pclasslabels = zeros(length(imagelist),1);
allvals = zeros(length(Ui), size(allweights,2));
for i=1:length(Ui)
    idx = find(imagelabels==Ui(i));
    val = prod(allweights(idx,:),1);
    allvals(i,:) = val;
    [a nlabel] = max(val);
    ppredlabels(i) = nlabel;
    pclasslabels(i) = classlabels(idx(1));
end

[C U D W] = conmatrix( pclasslabels, ppredlabels);
Ntot = sum(C(:));
Ncorr = U*Ntot;

D = round(D*1000)/10
U = round(U*1000)/10
W = round(W*1000)/10




fwrite(fid,[char(10) char(10) char(10)]);
fwrite(fid,['Field level confusion matrix,' char(10)]);
fwrite(fid,',');
for i=1:length(classes)
    fwrite(fid,classes{i});
    fwrite(fid,',');
end
fwrite(fid,['#samples,' char(10)]);

for i=1:length(classes)
    fwrite(fid,classes{i});
    fwrite(fid,',');
    for j=1:length(classes)
        fwrite(fid,num2str(D(i,j)));
        fwrite(fid,',');
    end
    fwrite(fid, [num2str(sum(C(i,:))) ',' char(10)]);
end
fwrite(fid,['Accuracy,' num2str(U) ',' char(10)]);
fwrite(fid,['Class weighted accuracy,' num2str(W) ',' char(10)]);



[val predlabels] = max(allvals, [], 2);
[s ind] = sort(val);

re = zeros(length(s),1);
pr = zeros(length(s),1);
for i=length(s):-1:1
    prl = ppredlabels(ind(i:end));
    cll = pclasslabels(ind(i:end));
    pr(i) = sum(prl==cll)/length(prl);
    re(i) = pr(i)*length(prl)/Ntot;
end

xmin = min(re)*100;
xmax = max(re)*100;
xmax = xmax + 0.1*(xmax-xmin);
ymin = min(pr)*100;
ymax = max(pr)*100;
ymax = ymax + 0.1*(ymax-ymin);

break
figure1 = figure('Color',[1 1 1],'Position', [100 100 650 650]);
axes('Parent',figure1,'FontSize',28,'LineWidth',3, 'XTick',[0 0.2 0.4 0.6 0.8 1]*100, 'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]*100, 'Position',[.20 .20 .75 .75]);
xpos = [(xmax-xmin)/2+xmin ymin-(ymax-ymin)*.15];
ypos = [-10 (ymax-ymin)/2+ymin];

plot(re*100,pr*100,'k^','LineWidth',2,'Color',[0 0 0]);
xlabel('Recall (%)','FontSize',32,'Position',xpos);% ,'FontName','Kartika');
ylabel('Precision (%)','FontSize',32,'Position',ypos);% ,'FontName','Kartika');
axis([xmin xmax ymin ymax]);
print -f1 -depsc PR.eps


rthr = 0.6;
idx = find(re>rthr);
i = idx(end);
pthr = s(i);

prl = ppredlabels(ind(i:end));
cll = pclasslabels(ind(i:end));

fwrite(fid,[char(10) char(10) char(10)]);

[Cpr Upr Dpr Wpr] = conmatrix( cll, prl);

Dpr = round(Dpr*1000)/10
Upr = round(Upr*1000)/10
Wpr = round(Wpr*1000)/10

fwrite(fid,['Field level with recall = ' num2str(rthr*100) '% / prob. thr = ' num2str(pthr) ',' char(10)]);
fwrite(fid,',');
for i=1:length(classes)
    fwrite(fid,classes{i});
    fwrite(fid,',');
end
fwrite(fid,['#samples,' char(10)]);

for i=1:length(classes)
    fwrite(fid,classes{i});
    fwrite(fid,',');
    for j=1:length(classes)
        fwrite(fid,num2str(Dpr(i,j)));
        fwrite(fid,',');
    end
    fwrite(fid, [num2str(sum(Cpr(i,:))) ',' char(10)]);
end
fwrite(fid,['Accuracy,' num2str(Upr) ',' char(10)]);
fwrite(fid,['Class weighted accuracy,' num2str(Wpr) ',' char(10)]);
fclose(fid)



break












[val predlabels] = max(allvals, [], 2);
[s ind] = sort(val);

re = zeros(length(s),1);
for i=length(s):-1:1
    prl = predlabels(ind(i:end));
    cll = classlabels(ind(i:end))';
    pr(i) = sum(prl==cll)/length(prl);
    re(i) = pr(i)*length(prl)/Ntot;
end

plot(re,pr,'.');
xlabel('Recall');
ylabel('Precision');




[Cpr Upr Dpr Wpr] = conmatrix( cllp, prlp);

Dpr = round(Dpr*1000)/10
Upr = round(Upr*1000)/10
Wpr = round(Wpr*1000)/10




idx = find(re>0.6);
i = idx(end);

prl = predlabels(ind(i:end));
cll = classlabels(ind(i:end));

if 1==2
idx = find(val>0.6);
prl = predlabels(idx);
cll = classlabels(idx);
end

fwrite(fid,[char(10) char(10) char(10)]);


[C U D W] = conmatrix( cll, prl);
C(idx_rmclass,:) = [];
C(:,idx_rmclass) = [];
D(idx_rmclass,:) = [];
D(:,idx_rmclass) = [];

W = D.*eye(size(D));
W = sum(W(:)) / size(D,1);

C
U = round(U*1000)/10
D = round(D*1000)/10
W = round(W*1000)/10

fwrite(fid,['Unweighted results with recall set at 60%,' char(10)]);
fwrite(fid,',');
for i=1:length(classes)
    fwrite(fid,classes{i});
    fwrite(fid,',');
end
fwrite(fid,char(10));

for i=1:length(classes)
    fwrite(fid,classes{i});
    fwrite(fid,',');
    for j=1:length(classes)
        fwrite(fid,num2str(C(i,j)));
        fwrite(fid,',');
    end
    fwrite(fid,char(10));
end
fwrite(fid,['Accuracy,' num2str(U) ',' char(10) char(10) char(10)]);


fwrite(fid,['Class weighted results with recall set at 60%,' char(10)]);
fwrite(fid,',');
for i=1:length(classes)
    fwrite(fid,classes{i});
    fwrite(fid,',');
end
fwrite(fid,char(10));

for i=1:length(classes)
    fwrite(fid,classes{i});
    fwrite(fid,',');
    for j=1:length(classes)
        fwrite(fid,num2str(D(i,j)));
        fwrite(fid,',');
    end
    fwrite(fid,char(10));
end
fwrite(fid,['Accuracy,' num2str(W) ',' char(10)]);



figure

if 1==1

[C U D W] = conmatrix( pclasslabels, ppredlabels);


[val predlabels] = max(allvals, [], 2);
[s ind] = sort(val);


Ntot = sum(C(:));
Ncorr = U*Ntot;
re = zeros(length(s),1);
pr = zeros(length(s),1);
for i=length(s):-1:1
    prl = ppredlabels(ind(i:end));
    cll = pclasslabels(ind(i:end));
    pr(i) = sum(prl==cll)/length(prl);
    re(i) = pr(i)*length(prl)/Ntot;
end

plot(re,pr,'.');
xlabel('Recall');
ylabel('Precision');

idx = find(re>0.6);
i = idx(end);

prl = ppredlabels(ind(i:end));
cll = pclasslabels(ind(i:end));

[C U D W] = conmatrix( cll,prl);

C(end+1,:) = 0;
C(:,end+1) = 0;
D(end+1,:) = 0;
D(:,end+1) = 0;

C(idx_rmclass,:) = [];
C(:,idx_rmclass) = [];
D(idx_rmclass,:) = [];
D(:,idx_rmclass) = [];

W = D.*eye(size(D));
W = sum(W(:)) / size(D,1);

C
U = round(U*1000)/10
D = round(D*1000)/10
W = round(W*1000)/10

fwrite(fid,['Unweighted results with recall set at 60%,' char(10)]);
fwrite(fid,',');
for i=1:length(classes)
    fwrite(fid,classes{i});
    fwrite(fid,',');
end
fwrite(fid,char(10));

for i=1:length(classes)
    fwrite(fid,classes{i});
    fwrite(fid,',');
    for j=1:length(classes)
        fwrite(fid,num2str(C(i,j)));
        fwrite(fid,',');
    end
    fwrite(fid,char(10));
end
fwrite(fid,['Accuracy,' num2str(U) ',' char(10) char(10) char(10)]);


fwrite(fid,['Class weighted results with recall set at 60%,' char(10)]);
fwrite(fid,',');
for i=1:length(classes)
    fwrite(fid,classes{i});
    fwrite(fid,',');
end
fwrite(fid,char(10));

for i=1:length(classes)
    fwrite(fid,classes{i});
    fwrite(fid,',');
    for j=1:length(classes)
        fwrite(fid,num2str(D(i,j)));
        fwrite(fid,',');
    end
    fwrite(fid,char(10));
end
fwrite(fid,['Accuracy,' num2str(W) ',' char(10)]);





end



fclose(fid);








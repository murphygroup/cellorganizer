clear

addpath lib/classification
addpath lib/classification/libsvm-mat-2.84-1
addpath lib/results

load fieldclass.mat

[C U D W] = conmatrix( classlabels, predlabels);
% C(idx_rmclass,:) = [];
% C(:,idx_rmclass) = [];
% D(idx_rmclass,:) = [];
% D(:,idx_rmclass) = [];
% W = D.*eye(size(D));
% W = sum(W(:)) / size(D,1);
% classes(idx_rmclass) = [];

Ntot = sum(C(:));
Ncorr = U*Ntot;

classes
D = round(D*1000)/10
U = round(U*1000)/10
W = round(W*1000)/10

fid = fopen('results_field.csv','w');

fwrite(fid,['Confusion matrix for all classified data,' char(10)]);
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
    cll = classlabels(ind(i:end));
    pr(i) = sum(prl==cll)/length(prl);
    re(i) = pr(i)*length(prl)/Ntot;
end

break
xmin = min(re)*100;
xmax = max(re)*100;
xmax = xmax + 0.1*(xmax-xmin);
ymin = min(pr)*100;
ymax = max(pr)*100;
ymax = ymax + 0.1*(ymax-ymin);

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

prl = predlabels(ind(i:end));
cll = classlabels(ind(i:end));




fwrite(fid,[char(10) char(10) char(10)]);

[Cpr Upr Dpr Wpr] = conmatrix( cll, prl);
% Cpr(idx_rmclass,:) = [];
% Cpr(:,idx_rmclass) = [];
% Dpr(idx_rmclass,:) = [];
% Dpr(:,idx_rmclass) = [];
% Wpr = Dpr.*eye(size(Dpr));
% Wpr = sum(Wpr(:)) / size(Dpr,1);

Dpr = round(Dpr*1000)/10
Upr = round(Upr*1000)/10
Wpr = round(Wpr*1000)/10

fwrite(fid,['Confusion matrix with recall = ' num2str(rthr*100) '% / prob. thr = ' num2str(pthr) ',' char(10)]);
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
fclose(fid);

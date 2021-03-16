function [slffeats, slfnames, histfeats, histnames] = spatialStatistics( img1, img2, DSF)

s = size(img1);

mask = zeros(size(img1));
mask(1:DSF:s(1),1:DSF:s(2)) = 1;
[r1,c1] = find(img1 & mask);
[r2,c2] = find(img2 & mask);


coord1 = [r1,c1];
coord2 = [r2,c2];

odist = pdist( [coord1; coord2]);
odist = squareform(odist);
% odist = odist(1:size(coord1,1),size(coord1,1)+1:end);



% odist = triu(odist+1, 1) - 1;
% odist(odist==-1) = [];





NBINS = 39;
bsize = 1/NBINS;
rmout = 0.05;


if isempty(odist)
    slffeats = zeros(1,4);
    histfeats = zeros(1,80);
    slfnames = [];
    histnames = [];
    return
end

[s ind] = sort(odist);
qt = floor(length(odist)*rmout/2);
ind = ind(1:end-qt+1);
odist = odist(ind);


odist = odist - min(odist);
odist = odist / max(odist);

[c b] = hist(odist,0:bsize:1);

histfeats1 = c / sum(c);
histfeats2 = c / sum(bsize*c);%*(1-rmout);
histfeats3 = c / sqrt(sum(c.^2));

mu = mean(odist);
st = std(odist);

al = mu*(mu*(1-mu)/st^2-1);
be = (1-mu)*(mu*(1-mu)/st^2-1);

slffeats = [mu st al be];
slfnames = {'Mean','St Dev',' Alpha','Beta'};


histfeats = [histfeats1 histfeats2 histfeats3];
histnames = [repmat({'histfeatsNormCounts'},[1 NBINS+1]) ...
             repmat({'histfeatsNormArea'},[1 NBINS+1]) ...
             repmat({'histfeatsNormMagArea'},[1 NBINS+1])];
for i=1:NBINS+1
    idx = i:NBINS+1:length(histfeats);
    for j=1:length(idx)
        histnames{idx(j)} = [histnames{idx(j)} '_' num2str(b(i))];
    end
end

return

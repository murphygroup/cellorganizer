clear

d = dir('data/maskvis/*jpg');
rand('seed',1);

rp = randperm(length(d));
ind = rp(1:10);

allcounts = zeros(length(ind),6);
for i=1:length(ind)
    str = ['data/maskvis/' d(ind(i)).name];
    img = imread(str);

    img = img(:,:,[3 2 1]);
    
    manualstr = ['data/manualmasks/' d(ind(i)).name(1:end-3) 'png'];
    
    if exist(manualstr,'file')
        manual = imread(manualstr);
    else
        manual = manualcrop(img);
        imwrite(manual,manualstr);
    end

    str2 = ['data/masks/' d(ind(i)).name(1:end-3) 'png'];
    water = imread(str2);
    
    [cPrecision,cRecall,pPrecision,pRecall,counts] = evalseg( manual, water);
    allcounts(i,:) = counts;
end

counts = sum(allcounts,1);
cTP = counts(1);
cFP = counts(2);
cFN = counts(3);
pTP = counts(4);
pFP = counts(5);
pFN = counts(6);

cRecall = cTP/(cTP + cFN);
cPrecision = cTP/(cTP+cFP);

pRecall = pTP/(pTP + pFN);
pPrecision = pTP/(pTP+pFP);

[cPrecision cRecall pPrecision pRecall]

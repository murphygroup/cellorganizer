function objects=tz_wavfindobj(img,cropimg,dnaproc,dnat,nobj,wavname, ...
    level,featset,alpha,option,t)
%TZ_WAVFINDOBJ Wavelet-based object detection.
%   OBJECTS = TZ_WAVFINDOBJ(IMG,CROPIMG,DNAPROC,DNAT,NOBJ,WAVNAME,
%       LEVEL,FEATSET,ALPHA,OPTION,T) returns the cell array of objects
%   detected from IMG based on wavelet-based features of each pixel.
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function tz_wavfindobj(img,nobj,wavname)
%
%OVERVIEW:
%   wavelet-based object detection
%PARAMETERS:
%   img - input image
%   cropimg - mask image
%   dnaproc - processed dna image
%   nobj - object number
%   wavname - wavelet basis
%   featset - pixel feat set
%   alpha - parameter for feature calculation
%   option - preprocessing
%   t - parameters for smoothing
%RETURN:
%   objects - object cell array
%DESCRIPTION:
%
%HISTORY:
%   11-NOV-2004 Initial write TINGZ
%   04-DEC-2004 Modified TINGZ
%       - add dnaproc

imgsize=size(img);

switch option
case 1
    %%% subtract background
    testimg = mb_imgbgsub(img, 'common');
    
    %remove background
    pfeats=tz_wavpfeats(testimg,dnaproc,dnat,1,wavname,featset,alpha);
    np=size(pfeats,1);
    
    rand('state',0);
    seeds=tz_randsel(np,2);
    selfeats=pfeats;
    centers=selfeats(seeds,:);
    options=zeros(1,14);
    options(14)=50;
    [centers,options,post,errlog] = netlab_kmeans(centers,selfeats,options);
    label=tz_post2label(post);
    
    if sum(label==1)<sum(label==2)
        lablel(label==2)=0;
    else
        label=label-1;
    end
    label(cropimg(:)==0)=0;
case 2
    testimg=ml_preprocess(img, cropimg, 'ml', 'yes');
    pfeats=tz_wavpfeats(img,dnaproc,dnat,level,wavname,featset,alpha);
    np=size(pfeats,1);
    label=testimg(:)>0;
case 3
    ker=eval(['fspecial(' tz_cell2str(t) ')']);
    img2=filter2(ker,img,'same');
    testimg=ml_preprocess(img, cropimg, 'ml', 'yes');
    pfeats=tz_wavpfeats(img2,dnaproc,dnat,level,wavname,featset,alpha);
    np=size(pfeats,1);
    label=testimg(:)>0;
end

objpfeats=pfeats(label==1,:);

rand('state',0)
selfeats=objpfeats;

if nobj>=1
    np=size(objpfeats,1);
    seeds=tz_randsel(np,nobj);
    centers=selfeats(seeds,:);
    options=zeros(1,14);
    options(14)=50;
    
    [centers,options,post,errlog] = netlab_kmeans(centers,selfeats,options);
else
    
    [aics,centers,posts]=tz_kmeansaic(selfeats,1:10,2,0,0);
    [minaic,i,j]=tz_min(aics);
    nobj=i;
    post=posts{i,j};
end

labels=tz_post2label(post);

img2=zeros(imgsize);
img2(label==1)=labels;
%figure
% imshow(img2,[]);
image(img2*50/max(img2(:)))

for i=1:nobj
    intensities{i}=testimg(img2==i);
    if ~isempty(intensities{i})
        meani(i)=mean(intensities{i});
    else
        warning('empty cluster');
        meani(i)=0;
    end
end

ranklabels=tz_rank(meani);
for i=1:nobj
    [r,c]=find(img2==i);    
    objects{ranklabels(i)}=[r,c,intensities{i}];
end

objects=rmempty(objects);

drawnow

function c2=rmempty(c)

k=1;
for i=1:length(c)
    if ~isempty(c{i})
        c2{k}=c{i};
        k=k+1;
    end
end

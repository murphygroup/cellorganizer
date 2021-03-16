function pfeats=tz_wavpfeats(img,dnaproc,t,level,wname,featset,alpha)
%TZ_WAVPFEATS Wavelet-based features for each pixel in an image.
%   PFEATS = TZ_WAVPFEATS(IMG,DNAPROC,T,LEVEL,WNAME,FEATSET,ALPHA) returns
%   the wavelet-based features of all pixels in the [image] IMG. DNAPROC is
%   the processed image of DNA channel. T is the option for calculating DNA
%   related features. LEVEL is the wavelet decomposition level. WNAME
%   specifies the name of wavelet bases. FEATSET is a [string array] of
%   feature set names. ALPHA is the nonlinear transformation parameter. The
%   returned value PFEATS is a [feature matrix].

%   07-NOV-2004 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 7
    error('Exactly 7 arguments are required')
end

if length(wname)==5
    switch wname
    case 'gabor'
        [out,avgout,freqs,thetas]=tz_gaborfiltbank(img,level);
        switch featset{1}
        case 'all'
            
            for i=1:length(out)
                img2=out{i};
                pfeats(:,i)=img2(:);
            end
        case 'avg'
            for i=1:length(avgout)
                img2=avgout{i};
                pfeats(:,i)=img2(:);
            end
        end
    case 'fhaar'
        filters={[1 2 1;2 4 2;1 2 1],[-1 -2 -1;2 4 2;-1 -2 -1;], ...
            [-1 2 -1;-2 4 -2;-1 2 -1],[1 -2 1;-2 4 -2;1 -2 1]};
        switch featset{1}
        case 'all'
            featsel=1:length(filters);
        case 'avg'
            featsel=1:length(filters);
        case 'mean'
            featsel=1:length(filters);
        case 'low'
            featsel=1;
        case 'high'
            featsel=2:length(filters);
        end
        
        for i=featsel
            img2=conv2(img,filters{i},'same');
            pfeats(:,i)=img2(:);
        end
        
        if strcmp(featset{1},'avg')
            
            pfeats(:,2)=mean(pfeats(:,2:end),2);

            pfeats(:,3:end)=[];
        end
        if strcmp(featset{1},'sum')
            
            pfeats(:,2)=sum(pfeats(:,2:end),2);

            pfeats(:,3:end)=[];
        end
    case 'inten'
        pfeats=img(:);
    end
   
else
    if level==0
        pfeats=img(:);
        return;
    end
    
    %wavelet decomposition
    [c,s]=wavedec2(img,level,wname);
    
    %approximation
    a=c(1:prod(s(1,:)));
    
    %details
    for i=1:level
        [h{i},v{i},d{i}]=detcoef2('all',c,s,i);
    end
    
    pfeats=[];
    for i=1:length(featset)
        switch(featset{i})
        case 'coef'
            pfeats=[pfeats,img(:)];
            for i=1:level
                sh=imresize(h{i},size(img),'nearest');
                sv=imresize(v{i},size(img),'nearest');
                sd=imresize(d{i},size(img),'nearest');
                pfeats=[pfeats,sh(:),sv(:),sd(:)];
            end
            sa=imresize(reshape(a,s(1,:)),size(img),'nearest');
            pfeats=[pfeats,sa(:)];
        case 'app'
            x=idwt2(reshape(a,s(1,:)),[],[],[],'haar');
            pfeats=[pfeats,x(:)];
        case 'pos'
            pfeats=[pfeats,repmat((1:size(img,1))',size(img,2),1)/50, ...
                repmat((1:size(img,2))',size(img,1),1)/50];
        end
    end
end

if alpha>0
    pfeats=tz_filtimg2featimg(pfeats,alpha);
end

if strcmp(featset{end},'dna')
    bwdna=dnaproc>0;
    filldna=imfill(bwdna,'holes');
    switch(t(1))
    case 1
        dnapfeats=ones(size(dnaproc));
        dnapfeats(filldna~=0)=-dnapfeats(filldna~=0);
        dnapfeats=dnapfeats(:);
    case 2
        dnaedge=bwperim(filldna);
        dnadist=bwdist(dnaedge);
        dnapfeats=dnadist;
        dnapfeats(filldna~=0)=-dnapfeats(filldna~=0);
        dnapfeats=dnapfeats(:);
    case 3
        dnaedge=bwperim(filldna);
        dnadist=bwdist(dnaedge);
        dnapfeats=dnadist;
        dnapfeats(filldna~=0)=-dnapfeats(filldna~=0);
        dnapfeats=dnapfeats(:);
        dnapfeats=tz_filtimg2featimg(dnapfeats,t(3));
        
    end
    
    dnapfeats=dnapfeats*max(pfeats(:))/t(2);
    pfeats=[pfeats,dnapfeats];    
end


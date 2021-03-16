function alpha=tz_decommix_lk(sample,training,method,minalpha,weights,t)
%TZ_DECOMMIX_LK Mixture decomposition by general likelihood.
%   ALPHA = TZ_DECOMMIX_LK(SAMPLE,TRAINING,METHOD,MINALPHA,WEIGHTS,T)
%   returns a vector of estimated coefficients of the fundamental
%   components in the mixture pattern SAMPLE, which is a data matrix.
%   The parameters of fundamental patterns are learned from data TRAINING.
%   Usually TRAING is a cell array of data matrices. Each matrix is for
%   one component. But TRAINING itself will also be the the parameters if 
%   METHOD is 'ready'. For more options about METHOD and T, see TZ_TRAINLK.
%   MINALPHA is the minimal possible coefficient that is greater than 0.
%   WEIGHTS specifies weights for each sample. Higher weight means higher
%   influence. If WEIGHTS is empty ([]), it means equal weights.
%   
%   See also TZ_TRAINLK.

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University


if nargin < 6
    error('Exactly 6 arguments are required')
end

if isempty(weights)
    weights = ones(size(sample,1),1);
else
    weights=weights/sum(weights);
end

%number of base components
m=length(training);
tm=m;

if strcmp(method,'ready')
    params=training; 
else
    for i=1:m
        if ~isempty(t)
            trainargs = {t{:},i};
        else
            trainargs = t;
        end
        params{i}=tz_trainlk(training{i},method,trainargs);
    end
end

for i=1:m
    lk(:,i)=tz_evalk(sample,params{i});
%     if ~isempty(weights)
%         lk(:,i)=lk(:,i).^weights;
%     end
    %     lk(:,i)=lk(:,i)/sum(lk(:,i));
end

invalidLkIdx = sum(abs(lk),2)==0;
lk(invalidLkIdx,:)=[];
if ~isempty(weights)
    weights(invalidLkIdx)=[];
end

%minimal weight. any weight less than minahpha will be set to zero
if isempty(minalpha)
    minalpha=.1/(size(sample,1));
end


%intialize
maxiter=100;
mine=1e-5;
sel=1:m;
succ=1;
i=1;

%Initialize alpha
alpha=ones(1,m)/m;

loglk=[];

halfweights = weights.^(1/2);
halfweightsmat = repmat(halfweights,1,size(lk,2)-1);

while i<maxiter
    
    %%%%%%%%%%%%%loglk%%%%%%%%%%%%
%     length(weights)
%     size(lk)
    loglk=[loglk,sum(weights.*log(lk*alpha'))];
    plot(loglk);
    drawnow
    %%%%%%%%%%%%%%%%%%%%
    
    
    dvm=lk*alpha';
    dm=[];
    phm=[];
    dm = tz_addcol(lk(:,1:end-1),-lk(:,end));
%     for j=1:size(lk,1)
%         dm=[dm;lk(j,1:end-1)-lk(j,end)];
%         phm(:,:,j)=weights(j)*dm(j,:)'*dm(j,:)/dvm(j)^2;
%     end
    
    grd=[];
    for k=1:m-1
        grd(k)=sum(weights.*dm(:,k)./dvm);
    end
    
    dvmmat = repmat(dvm,1,size(lk,2)-1);
    wdm = halfweightsmat(:,1:size(dm,2)).*dm./dvmmat;
    hm = -wdm'*wdm;
    
    if det(hm)==0 | size(sample,1)<m-1
        hm=-eye(size(hm,1));
    end

%     hm=-sum(phm,3);

    dalpha=grd*inv(hm)';
    
    tmpalpha=alpha(1:end-1)-dalpha;
    while any(tmpalpha<0) | sum(tmpalpha)>1
        dalpha=dalpha/2;
        tmpalpha=alpha(1:end-1)-dalpha;
    end
    
    alpha=[tmpalpha,1-sum(tmpalpha)];
    tmpalpha=alpha;
    
    i=i+1;
    if any(isnan(alpha))
        warning('invalid alpha');
        alpha=ones(1,tm)/tm;
        succ=0;
        break;
    end
    
    if i==maxiter | all(abs(dalpha)<mine)
        %         tmpalpha
        %         minalpha
        if(any(tmpalpha<minalpha))
            [minea,minind]=min(tmpalpha);
            alpha(minind)=[];
            alpha=alpha/sum(alpha);
            sel(minind)=[];
            lk(:,minind)=[];
            
            m=size(lk,2);
            i=1;
            if length(alpha)==1
                alpha=1;
                break
            end
        else
            if i>=maxiter
                disp('max iteration reached');
            end
            break;
        end
        
    end
end

if succ==1
    m=length(training);
    ealpha=alpha;
    alpha=zeros(1,m);
    alpha(sel)=ealpha;
else
    warning('fitting failed');
end

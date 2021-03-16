function alpha = ml_decommixlk(sample,training,param)
%ML_DECOMMIXLK
%   ALPHA = ML_DECOMMIXLK(SAMPLE,TRAINING)
%   
%   ALPHA = ML_DECOMMIXLK(SAMPLE,TRAINING,PARAM)
%   
%   See also

%   15-Dec-2006 Initial write T. Zhao
%   Copyright (c) 2006 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.
%   
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu


if nargin < 2
    error('2 or 3 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct('weights',[], ...
    'f',struct('name','mvn'),'minalpha',[]));

if isempty(param.weights)
    weights = ones(size(sample,1),1);
else
    weights=param.weights/sum(param.weights);
end

%number of base components
m=length(training);
tm=m;

for i=1:m
    f{i} = ml_estpdf(training{i},param.f);
end

for i=1:m
    lk(:,i) = ml_pdf(sample,f{i});
end
                  
invalidLkIdx = find(sum(abs(lk),2)==0);
lk(invalidLkIdx,:)=[];
if ~isempty(weights)
    weights(invalidLkIdx)=[];
end

%minimal weight. any weight less than minahpha will be set to zero
minalpha = param.minalpha;
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
    length(weights)
    size(lk)
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

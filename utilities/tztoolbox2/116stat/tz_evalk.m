function lk=tz_evalk(Y,params,islog)
%TZ_EVALK Calculate likelihood or log likelihood.
%   LK = TZ_EVALK(Y,PARAMS) returns the likelihood of each data point in
%   the [feature matrix] Y based on the distribution spcefified by PARAMS,
%   which has the following fields:
%   'method' - what kind of distribution
%       'mnom' : multinomial distribution
%           'p' - multinomial probability vector      
%       'mnomn' : multinomial with numbers
%           'p' - same as 'mnom'
%           'objnum' - histogram of object types. See TZ_DISCRDENS for more details.
%       'mvn' or 'mvnw' : multivariate normal distribution
%           'mu' : mean
%           'sigma' : covariance matrix
%       'mvnm' or 'mvno' : Gaussian mixture
%           'mus' : a matrix. each column is a mean.
%           'covs' : a three dimensional matrix. covs(:,:,i) is the covariance matrix of the ith component.
%
%   LK = TZ_EVALK(Y,PARAMS,1) returns the log likelihood. TZ_EVALK(Y,PARAMS,0) is the same as TZ_EVALK(Y,PARAMS).
%   
%
%   See also TZ_TRAINLK

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

% function logl=tz_evalk(params,Y)

if ~exist('islog','var')
    islog=0;
end

if length(params.method)>=4
    if strcmp(params.method(1:4),'mnom')
        for i=1:size(Y,1)
            %remove zeros from parameters
            s=Y(i,params.p~=0);
            tmpp=params.p(params.p~=0);
            
            %Special case: 0 probability
            if any(params.p==0)
                if any(Y(i,params.p==0)~=0)
                    if islog==0
                        lk(i)=0;
                    else
                        lk(i)=-Inf;
                    end
                    continue
                end
            end
            
            %likelihood
            if islog==0
                lk(i)=prod(tmpp.^s);
            else
                lk(i)=sum(s.*log(tmpp));
            end
        end 
    end
end

switch params.method
case 'mnom'
    lk=lk';
case 'mnomn'    %mutinomial with numbers
    if islog==0
        lk=lk.*tz_discrdens(params.objnum,50,sum(Y,2));
    else
        lk=lk+log(tz_discrdens(params.objnum,50,sum(Y,2)));
    end
    lk=lk';
case {'mvn','mvnw'}      %normal distribution
    
    tY=Y;
    if ~isempty(params.constidx)
        tY(:,params.constidx)=[];
    end

    lk=mvnpdf(tY,params.mu,params.sigma);
    if ~isempty(params.constidx)
        for i=1:size(Y,1)     
            if any(params.constx~=Y(i,params.constidx))
                lk(i)=0;
            end
        end
    end
    
    if islog~=0
        lk=log(lk);
    end
    
case {'mvnm','mvno'} %gaussian mixture
    constidx=params.constidx;
    tY=Y;
    if ~isempty(constidx)
        tY(:,constidx)=[];
    end
    
    lk=[];
    for i=1:params.m
        lk=[lk,mvnpdf(tY,params.mus(:,i)',params.covs(:,:,i))];
    end
    
    lk=lk*params.pp';

    if ~isempty(constidx)
        for i=1:size(Y,1)
            if any(params.constx~=Y(i,constidx))
                logl(i)=0;
            end
        end
    end
    
    if islog~=0
        lk=log(lk);
    end
end  

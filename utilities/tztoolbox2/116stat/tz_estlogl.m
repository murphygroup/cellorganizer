function logl=tz_estlogl(X,Y,method)
%TZ_ESTLOGL Log-likelihood estimation.
%   LOGL = TZ_ESTLOGL(X,Y,METHOD)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%TZ_ESTLOGL: calculate the log-likelihood of Y under the density estimated from X
%

if length(method)>=4
    if strcmp(method(1:4),'mnom')
        totalnum=sum(X,1);
        p=totalnum/sum(totalnum);
        for i=1:size(Y,1)
            s=Y(i,p~=0);
            tmpp=p(p~=0);
            if any(p==0)
                if any(Y(i,p==0)~=0)
                    logl(i)=-Inf;
                    continue
                end
            end
            logl(i)=sum(s.*log(tmpp));
        end 
    end
end

switch(method)
    
case 'mnomn'    %mutinomial with numbers
    objnum=[sum(X,2),ones(size(X,1),1)];
    logl=logl+log(tz_discrdens(objnum,50,sum(Y,2)));
case 'mvn'      %normal distribution
    constidx=tz_constidx(X);
    tX=X;
    tY=Y;
    if ~isempty(constidx)
        tX(:,constidx)=[];
        tY(:,constidx)=[];
    end
    mu=mean(tX,1);
    sigma=cov(tX);
    logl=log(mvnpdf(tY,mu,sigma))';
    if ~isempty(constidx)
        for i=1:size(Y,1)
            
            if any(X(1,constidx)~=Y(i,constidx))
                logl(i)=-Inf;
            end
            
        end
    end
case 'mvnm' %gaussian mixture
    constidx=tz_constidx(X);
    tX=X;
    tY=Y;
    if ~isempty(constidx)
        tX(:,constidx)=[];
        tY(:,constidx)=[];
    end
    
    [bestk,bestpp,bestmu,bestcov,dl,countf] = mixtures4(tX',1,4,0.01,1e-6,0)
    lk=[];
    for i=1:bestk
        lk=[lk,mvnpdf(tY,bestmu(:,i)',bestcov(:,:,1))];
    end
    lk=lk*bestpp';
    logl=log(lk)';
    if ~isempty(constidx)
        for i=1:size(Y,1)
            
            if any(X(1,constidx)~=Y(i,constidx))
                logl(i)=-Inf;
            end
            
        end
    end
end
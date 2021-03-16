function [B,res]=tz_transformfeats(X,Y,option,t)

%TZ_TRANSFORM: do multivariate regression (or prediction)

projx=sum(abs(X),1);
projy=sum(abs(X),2);
if any(projx==0)
    warning('zero columns are removed')
    X(:,projx==0)=[];
end

if any(projy==0)
    warning('zero rows are removed')
    X(projy==0,:)=[];
    Y(projy==0,:)=[];
end

switch option
case 'ln'
    X=[ones(size(X,1),1),X];
    B=inv(X'*X)*X'*Y;
case 'sv'
    X=[ones(size(X,1),1),X];
    [V,S,U]=svd(cov([X,Y]));
    V=U(1:size(X,2),:);
    W=U(size(X,2)+1:end,:);
    B=V*W';
case 'pl'
    X=[ones(size(X,1),1),X];
    nF=min(size(X))/2;
    [T,P,W,U,Q,B,ssq,Ro,Rv,Lo,Lv] = mypls(X,Y,round(nF));
    B=W*Q';
case 'pr'
    if ~ischar(t)
        if t==1
            t=size(X,2);
        else
            [pc,score,latent] = princomp(X);
            latent=cumsum(latent)/sum(latent);
            t=find(latent>=t);
            t=t(1);
        end
    else
        t=str2num(t);
    end
    
    result=cpcr(X,Y,'plots',0,'k',t);
        
    B=[result.int;result.slope];
%     X=[ones(size(X,1),1),X];
end

Yh=tz_evaltransf(X,B);
res=sum(sum((Yh-Y).^2));
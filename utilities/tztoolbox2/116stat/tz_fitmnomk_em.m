function [alpha,loglk,succ]=tz_fitmnomk_em(y,p)

%function [alpha,loglk,succ]=tz_fitmnomk(y,p)
%OVERVIEW:
%   Solve multinomial mixture of m components, p: mxk matrix
%PARAMETERS:
%   y - data
%   p - parameters of components
%RETURN:
%   alpha - weigths
%   loglk - log likelihood
%   succ - failed or not
%DESCRIPTION:
%   model: g(x)=alpha1*f1(x)+alpha2*f2(x)+...+alphan*fn(x); alpha1+...+alphan=1
%           fi ~ multinomial
%   Using EM algorithm to find maximum point
%
%HISTORY:
%   15-MAY-2004 Initial write TINGZ
%   03-JUN-2004 Modified TINGZ
%       - debug singular hessian matrix

succ=1;
n=sum(y);
if(n==0)
    alpha=nan;
    loglk=-Inf;
    succ=0;
    warning('invalid data');
    return;
end

%keep original parameters and data
orgy=y;
orgp=p;
loglk=0;

%minimal weight. any weight less than minahpha will be set to zero
minalpha=0.1/(sum(y));

tp=sum(orgp,1);

%The logl likelihood is -Inf if there is any sample drawn from the probalbity 0
if any(tp==0)
    if any(orgy(tp==0)~=0)
        loglk=-Inf;
    end
end

%remove 0s
% p(:,tp==0 | y==0)=[];
% y(tp==0 | y==0)=[];

[m,k]=size(orgp);

%Initiate alpha
alpha=ones(1,m)/m;
ialpha=alpha;
sel=1:m;

if(any(sum(p,2)==0))
    alpha(sum(p,2)==0)=[];
    sel(sum(p,2)==0)=[];
    p(sum(p,2)==0,:)=[];
    m=size(p,1);
end

if m<=1
    if(isempty(p))
        succ=0;
        alpha=ialpha;
        loglk=-Inf;
        warning('no data or parameter left');
        return;
    end
    
    alpha=1;
    
    if loglk==0
        loglk=sum(y.*log(p));
    end
else    
    maxiter=500;
    mine=1e-5;
    
    %%%%%%%%%%%%%Plot loglk surface - for debugging%%%%%%%%
    % [a1,a2]=meshgrid(0:0.01:1,0:0.01:1);
    % 
    % for j=1:101
    %     for i=1:101
    %         allalpha=[a1(i,j),a2(i,j),1-a1(i,j)-a2(i,j)];
    %         if any(allalpha<0) | sum(allalpha)>1
    %             sloglk(i,j)=NaN;
    %         else
    %             factor1=allalpha*p;
    %             sloglk(i,j)=sum(y.*log(factor1));
    %         end
    %     end
    % end
    % surf(a1,a2,sloglk);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    i=1;
    while i<maxiter
        if(isempty(p))
            succ=0;
            warning('no data or parameter left');
            break;
        end
        
        if(size(p,2))==1
            [maxv,maxi]=max(alpha);
            alpha=alpha*0;
            alpha(maxi)=1;
            break;
        end
        
        
        tmpalpha=alpha;
        nalpha=length(alpha);
        mp=(n*tmpalpha*p);
        for u=1:nalpha
            alpha(u)=sum(y.*p(u,:)*tmpalpha(u)./mp);
        end
        
        dalpha=alpha-tmpalpha;
%         min(abs(dalpha))
        i=i+1;
        
     
        
        if i==maxiter | all(abs(dalpha)<mine)
            if(any(tmpalpha<minalpha))
                alpha(tmpalpha<minalpha)=[];
                alpha=alpha/sum(alpha);
                sel(tmpalpha<minalpha)=[];
                p(tmpalpha<minalpha,:)=[];
                
                ntp=sum(p,1);
                if any(ntp==0)
                    p(:,ntp==0)=[];
                    y(ntp==0)=[];
                end
                
                alpha(sum(p,2)==0)=[];
                sel(sum(p,2)==0)=[];
                p(sum(p,2)==0,:)=[];
                m=size(p,1);
                i=1;
            else
                if i>=maxiter
                    disp('max iteration reached');
                end
                break;
            end
            
        end
    end
end

if succ==1
    m=size(orgp,1);
    
    ealpha=alpha;
    alpha=zeros(1,m);
    alpha(sel)=ealpha;
    
    %if m~=1
    %    alpha(m)=1-sum(alpha);
    %end
    
    if loglk==0
        loglk=sum(y.*log(ealpha*p));
    end
else
    warning('fitting failed');
end
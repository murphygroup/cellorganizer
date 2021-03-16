function [alpha,loglk,succ]=tz_fitmnomk(y,p,minalpha)

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
%   Using Newton's method and gradient ascent to find maximum point
%
%HISTORY:
%   15-MAY-2004 Initial write TINGZ
%   03-JUN-2004 Modified TINGZ
%       - debug singular hessian matrix

succ=1;

if(sum(y)==0)
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
if ~exist('minalpha','var')
    minalpha=.1/(sum(y));
else
    if isempty(minalpha)
        minalpha=.1/(sum(y));
    end
end

tp=sum(orgp,1);

%The logl likelihood is -Inf if there is any sample drawn from the probalbity 0
if any(tp==0)
    if any(orgy(tp==0)~=0)
        loglk=-Inf;
    end
end

%remove 0s
p(:,tp==0 | y==0)=[];
y(tp==0 | y==0)=[];

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
end

if m==2
    alpha=tz_fitmnom2(y,p(1,:),p(2,:));
    alpha=[alpha,1-alpha];
    
    if(loglk==0)
        loglk=sum(y.*log(alpha*p));
    end
end



%remove 0s
% zeroc=find(sum(p)==0);
% y(zeroc)=[];
% p(:,zeroc)=[];

if m>2

    
    maxiter=100;
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
        %allalpha=[alpha,1-sum(alpha)];
        factor1=alpha*p;
        
        if any(factor1==0)
            warning('factor 1 is 0');
        end
        
        %loglk(i)=sum(y.*log(factor1));
        d1=[];
        d2=[];
        
        %     yy=y;
        %     factor1(yy==0)=[];
        %     pp(yy==0)
        %     yy(yy==0)=[];
        
        %Calculate gradients
        for j=1:(m-1)
            d1(j)=sum(y.*(p(j,:)-p(m,:))./factor1);  
        end
        
        if sum(y~=0)>=m-1
            %Calculate Hessian matrix
            d2=zeros(m-1,m-1);
            for r=1:m-1
                for s=1:r
                    d2(r,s)=-sum(y.*(p(r,:)-p(m,:)).*(p(s,:)-p(m,:))./factor1.^2);
                end
            end
            d2=d2+(d2-diag(diag(d2)))';
            
            iH=inv(d2);
            if any(isinf(iH))
                iH=-eye(m-1);
                %break;
            end
        else
            iH=-eye(m-1);
        end
        

        
        dalpha=d1*iH;
        
        if any(isinf(dalpha)) | any(isnan(dalpha))
            succ=0;
            break;
        end
        
%         if all(abs(dalpha)<mine)
%             break
%         end 
        
        tmpalpha=alpha(1:end-1)-dalpha;
        while any(tmpalpha<0) | sum(tmpalpha)>1
            dalpha=dalpha/2;
            tmpalpha=alpha(1:end-1)-dalpha;
        end
        
        alpha=[tmpalpha,1-sum(tmpalpha)];
        tmpalpha=alpha;
        
        if any(isnan(alpha))
            warning('invalid alpha');
        end
        
        %alpha(tmpalpha<minalpha)=0;
        
        %  p=tz_normobjcom(p);
        
        i=i+1;
        if i==maxiter | all(abs(dalpha)<mine)
            if(any(tmpalpha<minalpha))
                [minea,minind]=min(tmpalpha);
                alpha(minind)=[];
                alpha=alpha/sum(alpha);
                sel(minind)=[];
                p(minind,:)=[];
                
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
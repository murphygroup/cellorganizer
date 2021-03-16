function [x,logll]=tz_redistrmnom(y,p,x0)

%[x,logll]=tz_redistrmnom(y,p)
%
%OVERVIEW:
%   Estimate multinoial mixture model
%PARAMETERS:
%   y: data
%   p: kxn matrix for k mixture models with n nomial
%   x0 - intial guess
%RETURN:
%   x - components. kxn matrix
%DESCRIPTION
%
%HISTORY:
%   24-MAR-2004 Initial write TINGZ
%   31-MAR-2004 Modified TINGZ
%   05-NOV-2004 Modified TINGZ
%       - add comments

n=size(p,2);
nmix=size(p,1);

if ~exist('x0','var')
    x0=zeros(nmix,n);
    for i=1:n 
        x0(:,i)=tz_redistrnum(y(i),nmix)';
    end
end

x=x0;

%condition for iteration
incflag=1;

%precalculation
for i=1:nmix
    for j=1:nmix
        ratiop{i,j}=p(i,:)./p(j,:);   
        ratiop{i,j}(p(i,:)==0)=0;
    end
end

%index for logll recording
% k=1;
% logll(k)=tz_mnomlogl(x,p)

while(incflag==1)
    incflag=0;
    
    %loop for checking all slots
    for i=1:nmix-1
        for j=(i+1):nmix
            %up direction
            x1=x(i,:);
            x2=x(j,:);
            n1=sum(x1);
            n2=sum(x2);
            if n1==0
                n1=1;
            end
            x1(x1==0)=1;
            incp=((n2+1)/n1)*ratiop{j,i}.*x1./(x2+1);
            incp(x(i,:)==0)=0;
            [maxincp1,pos1]=max(incp);
            if(maxincp1>1 & x(i,pos1)>0 &n1>0)
                x(i,pos1)=x(i,pos1)-1;
                x(j,pos1)=x(j,pos1)+1;
                incflag=1;
            end

            %down direction
            x1=x(j,:);
            x2=x(i,:);
            n1=sum(x1);
            n2=sum(x2);
            if n1==0
                n1=1;
            end
            x1(x1==0)=1;
            incp=((n2+1)/n1)*ratiop{i,j}.*x1./(x2+1);
            incp(x(j,:)==0)=0;
            [maxincp1,pos1]=max(incp);
            if(maxincp1>1 & x(j,pos1)>0 &n1>0)
                x(j,pos1)=x(j,pos1)-1;
                x(i,pos1)=x(i,pos1)+1;
                incflag=1;
            end
        end
    end
    
    %update log likelihood
%     k=k+1;
%     logll(k)=tz_mnomlogl(x,p);
%     
    if(incflag==0)
        break;
    end
end



logll=tz_mnomlogl(x,p);

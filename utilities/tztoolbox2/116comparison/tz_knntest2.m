function [pvalue,ts]=tz_knntest2(X1,X2,k,s)
%TZ_KNNTEST2 Obsolete.
%
%See also ML_KNNTEST2

%function [pvalue,ts]=tz_knntest2(X1,X2,k,s)
%
%OVERVIEW:
%   Nearest neighbor two smaple test 
%PARAMETERS:
%   X1 - the first group of samples n1xp
%   X2 - the second group of samples n2xp
%   k - k nearest neighbor
%   s - distance fuction for finding neighbors
%RETURN:
%   pvalue - p-value
%   ts - test statistic
%DESCRPTION:
%   Original paper - A multivariate two-sample test based on the number of nearest neighbor type coincidences
%                    Norbert Henze, 1988
%
%HISTORY:
%   03-APR-2004 Initial Write TINGZ
%   17-MAY-2004 Add comments TINGZ
%   17-JUL-2004 Modified TINGZ
%       - change ts

error(tz_genmsg('of','tz_knntest2','ml_knntest2');

if ~exist('s','var')
    s='eu';
end

if k==-1
    k=min([size(X1,1),size(X2,1)])-1;
end

X=zscore([X1;X2]);

n1=size(X1,1);
n2=size(X2,1);
n=n1+n2;

pwdist=squareform(tz_pdist(X,s,n1));
mindist=min(pwdist(pwdist~=0));
[index1,index2]=find(pwdist==0);
for i=1:length(index1)
    if index1(i)~=index2(i)
        pwdist(index1(i),index2(i))=unifrnd(0,mindist);
    end
end

knnmatrix=tz_rank(pwdist);
knnmatrix=knnmatrix-1;
B=[[ones(n1,n1),zeros(n1,n2)];[zeros(n2,n1),ones(n2,n2)]];
t=zeros(size(B));
for r=1:k
    a(:,:,r)=(knnmatrix==r);
    t=t+a(:,:,r).*B;
end

ts=sum(t(:));

ap=sum(a,3);



d=sum(ap,1);
c=sum((d-k).^2)/n/k;
v=sum(sum(ap.*ap'))/n/k;
m=(n1*(n1-1)+n2*(n2-1))/(n-1);
q=4*(n1-1)*(n2-1)/((n-2)*(n-3));
EL=k*m;
VarL=k*n1*n2*(q*(1+v-2*k/(n-1))+(1-q)*c)/(n-1);

pvalue=1-normcdf((ts-EL)/sqrt(VarL));

ts=ts/n/k;
function [pvalue,ts,res,stdres]=tz_testindep(data,option)

%function [pvalue,ts]=tz_testindep(data,option)
%
%OVERVIEW:
%   test indpendence of 2D discrete table
%PARAMETERS:
%   data - IxJ table
%   option - method of testing
%RETURN:
%   pvalue - p-value
%   ts - test statistic
%   res - residuals
%   stdres - standarized residuals
%DESCRIPTION:
%
%HISTORY:
%   04-NOV-2004 Initial write TINGZ

dataj=sum(data,1);
datai=sum(data,2);
ndata=sum(data(:));

pi=datai/ndata;
pj=dataj/ndata;
pe=pi*pj;
ne=pe*ndata;

switch(option)
case 'c' %chi square
    ts=sum((data(:)-ne(:)).^2./ne(:));
case 'g' %G square
    ts=sum(2*data(:).*log(data(:)./ne(:)));
case 'o' %order
    pdata=data/sum(data(:));
    for i=1:size(data,1)
        for j=1:size(data,2)
            subdata=pdata;
            subdata(1:i,:)=[];
            subdata(:,j:end)=[];
            pid(i,j)=sum(subdata(:));
%             subdata=pdata;
%             subdata(i:end,:)=[];
%             subdata(:,1:j)=[];
%             pid(i,j)=pid(i,j)+sum(subdata(:));
            
            subdata=pdata;
            subdata(i:end,:)=[];
            subdata(:,j:end)=[];
            pic(i,j)=sum(subdata(:));
%             subdata=pdata;
%             subdata(1:i,:)=[];
%             subdata(:,1:j)=[];
%             pic(i,j)=pic(i,j)+sum(subdata(:));
        end
    end
    pc=2*sum(sum(pdata.*pic));
    pd=2*sum(sum(pdata.*pid));
    ts=(pc-pd)/(pc+pd);
    se=sqrt(16*sum(sum(pdata.*(pd*pic-pc*pid).^2))/(pc+pd)^4);
    pvalue=(1-normcdf(abs(ts),0,se))*2;
end

res=(data-ne)./sqrt(ne);
stdres=res./sqrt((1-pi)*(1-pj));
if option~='o'
    df=(size(data,1)-1)*(size(data,2)-1);
    pvalue=1-chi2cdf(ts,df);
end
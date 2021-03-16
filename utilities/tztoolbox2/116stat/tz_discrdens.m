function density=tz_discrdens(data,winsize,n)
%TZ_DISCRDENS Discrete density estimation by kernel method.
%   DENSITY = TZ_DISCRDENS(DATA,WINSIZE,N) returns the probability densisties of the integer vector N. DATA is a 2-column matrix for density estimation. The first column contains values and the second column contains numbers. WNDSIZE is the window size of the binomial kernel.
%   
%   See also

%   15-MAR-2004 Initial write T. Zhao
%   28-MAR-2004 Modified TINGZ
%       - add comments   
%   Copyright (c) Murphy Lab, Carnegie Mellon University

kernel=binopdf(0:winsize,winsize,0.5);
halfwin=round(winsize/2);
data(:,1)=data(:,1)-1;
data(:,2)=data(:,2)/sum(data(:,2));

for i=1:length(n)
    for j=1:size(data,1)
        if data(j,1)+1<=halfwin
            weight(j)=data(j,2)*binopdf(n(i)-1,winsize,(data(j,1)+0.5)/winsize);    
        else
            weight(j)=data(j,2)*binopdf(n(i)-data(j,1)-1+halfwin,winsize,0.5);
        end
    end
    density(i)=sum(weight);
        
end

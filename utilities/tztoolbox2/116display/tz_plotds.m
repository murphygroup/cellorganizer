function tz_plotds(x,sn)
%TZ_PLOTDS Plot digital signal.
%   TZ_PLOTDS(X) plots the column vector X as digital signal.
%   
%   TZ_PLOTDS(X,SN) specifies the starting position of the digital signal
%   so that the first element of X will start at SN.

%   16-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('1 or 2 arguments are required')
end

if ~exist('sn','var')
    sn=0;
end

n=length(x);
ax=[sn+(0:n-1);sn+(0:n-1)];
ay=[zeros(1,n);x];

plot(ax(2,:),ay(2,:),'ro');

line(ax,ay,'color','b','linewidth',2);
line([sn-1 n+sn],[0 0],'color','k');

axis([sn-1 n+sn 0 1])
axis 'auto y'
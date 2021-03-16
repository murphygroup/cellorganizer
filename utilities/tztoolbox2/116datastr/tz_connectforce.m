function f=tz_connectforce(x1,x2,option)
%TZ_CONNECTFORCE Calculate attractive force between two vectors.
%   F = TZ_CONNECTFORCE(X1,X2,OPTION) returns the attractive force btween
%   two vectors X1 and X2, which should have the same size. OPTIONS
%   specifies the method of calculating the force:
%       'grv' - gravity
%       'crr' - correlation 

%   15-Sep-2005 Initial write T. Zhao
%   ??-???-2004 Initial write T. Zhao
%   04-NOV-2004 Modified T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('Exactly 3 arguments are required')
end

switch option
case 'grv'
    ms=x1.*x2;
    x1(ms==0)=[];
    x2(ms==0)=[];
    ms(ms==0)=[];
    if isempty(ms)
        f=0;
    else
        f=sum(ms./(x2-x1+1).^2);
    end
case 'crr'
    f=corrcoef(x1,x2);
    f=f(2);
end
        

    
